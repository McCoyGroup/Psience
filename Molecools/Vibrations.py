from McUtils.Coordinerds.CoordinateSystems import CoordinateSystem, CartesianCoordinates3D, CoordinateSet
from McUtils.Data import AtomData, UnitsData
from .Molecule import Molecule

import numpy as np, scipy.linalg as slag

__all__ = [
    "MolecularVibrations",
    "NormalModeCoordinates",
]

class MolecularVibrations:

    def __init__(self, molecule, basis, freqs = None, init = None):
        """Sets up a vibration for a Molecule object over the CoordinateSystem basis

        :param molecule:
        :type molecule: Molecule
        :param init:
        :type init: None | CoordinateSet
        :param basis:
        :type basis: CoordinateSystem
        """
        self._mol = molecule
        self._coords = init if init is not None else self._mol._coords
        self._basis = basis
        if freqs is None and hasattr(basis, "freqs"):
            freqs = basis.freqs
        self.freqs = freqs

    @property
    def coords(self):
        return self._coords
    @property
    def basis(self):
        return self._basis

    def __len__(self):
        return self._basis.matrix.shape[0]

    def displace(self, displacements = None, amt = .1, n = 1, which = 0):

        if displacements is None:
            displacements = self._basis.displacement(np.arange(1, n+1)*amt)
        displacements = np.asarray(displacements)
        if displacements.ndim == 1:
            disp_coords = CoordinateSet(np.zeros((n,) + self._basis.coordinate_shape), self._basis)
            disp_coords[:, which] = displacements
        else:
            disp_coords = CoordinateSet(displacements, self._basis)
        coords = np.broadcast_to(self._coords, (n, ) + self._coords.shape)
        disp_coords = disp_coords.convert(self._coords.system)
        displaced = coords + disp_coords
        return displaced

    def visualize(self, step_size = .1, steps = (5, 5), which = 0, anim_opts = None, mode = 'fast', **plot_args):
        from McUtils.Plots import Animator

        if isinstance(steps, (int, np.integer)):
            steps = [steps, steps]

        left  = np.flip(self.displace(amt=-step_size, n=steps[0], which=which), axis=0)
        right = self.displace(amt=step_size, n=steps[1], which=which)
        all_geoms = np.concatenate((left, np.broadcast_to(self._coords, (1, ) + self._coords.shape), right))

        figure, atoms, bonds = self._mol.plot(*all_geoms, objects = True, mode = mode, **plot_args)

        def animate(*args, frame = 0, _atoms = atoms, _bonds = bonds, _figure = figure):

            my_stuff = []
            nframes = len(_atoms)
            forward = ( frame // nframes ) % 2 == 0
            frame = frame % nframes
            if not forward:
                frame = -frame - 1

            if _atoms[frame] is not None:
                for a in _atoms[frame]:
                    coll = a.plot(figure)
                    my_stuff.append(coll)
            if _bonds[frame] is not None:
                for b in _bonds[frame]:
                    for bb in b:
                        p = bb.plot(figure)
                        try:
                            my_stuff.extend(p)
                        except ValueError:
                            my_stuff.append(p)

            return my_stuff

        if anim_opts is None:
            anim_opts = {}
        return Animator(figure, None, plot_method = animate, **anim_opts)

class NormalModeCoordinates(CoordinateSystem):
    """A prettied up version of a Coordinerds CoordinateSystem object
    Has function for generating, though, too
    """
    name="MolecularNormalModes"
    def __init__(self, molecule, coeffs, name=None, freqs=None):
        super().__init__(
            matrix=coeffs,
            name=self.name if name is None else name,
            basis=molecule.sys
        )
        self.freqs = freqs

    @classmethod
    def from_force_constants(cls,
                             molecule,
                             fcs,
                             atoms = None,
                             masses = None,
                             mass_units = "AtomicMassUnits",
                             inverse_mass_matrix = False,
                             energy_units = "Wavenumbers",
                             remove_transrot = True,
                             normalize = True,
                             **opts
                             ):
        """Generates normal modes from the specified force constants

        :param fcs:
        :type fcs:
        :param atoms:
        :type atoms:
        :param masses:
        :type masses:
        :param target_units:
        :type target_units:
        :param opts:
        :type opts:
        :return:
        :rtype:
        """

        if atoms is not None and masses is None:
            masses = np.array([AtomData[a, "Mass"] if isinstance(a, str) else a for a in atoms])

        if mass_units != "AtomicUnitOfMass" and mass_units != "ElectronMass":
            mass_conv = UnitsData.convert(mass_units, "AtomicUnitOfMass"),
        else:
            mass_conv = 1

        if masses is not None:
            masses = np.asarray(masses)
            masses = masses*mass_conv
            if masses.ndim == 1:
                masses = np.broadcast_to(masses, (len(masses), 3)).T.flatten()
                masses = np.diag(masses)
                inverse_mass_matrix = True
        else:
            masses = np.eye(len(fcs))

        freqs, modes = slag.eigh(fcs, masses, type=(1 if inverse_mass_matrix else 3))
        if normalize:
            normalization = np.broadcast_to(1/np.linalg.norm(modes, axis=0), modes.shape)
            modes = modes * normalization

        freqs = np.sign(freqs) * np.sqrt(np.abs(freqs))
        sorting = np.argsort(freqs)
        if remove_transrot:
            sorting = sorting[6:]

        freqs = freqs[sorting]
        if energy_units is not None and energy_units != "Hartrees":
            freqs = freqs * UnitsData.convert("Hartrees", energy_units)
        modes = modes[:, sorting]

        return cls(molecule, modes, freqs = freqs, **opts)