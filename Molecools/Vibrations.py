from ..Coordinerds.CoordinateSystems import CoordinateSystem, CartesianCoordinates3D, CoordinateSet
from McUtils.Data import AtomData, UnitsData
from .Molecule import Molecule

import numpy as np

__all__ = [
    "NormalModeCoordiantes",
    "GMatrix"
]

class VibrationalModes:

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
        self.freqs = freqs if freqs is not None else basis.freqs

    def __len__(self):
        return self._basis.matrix.shape[0]

    def displace(self, amt = .1, n = 1, which = 0):

        displacements = CoordinateSet(np.zeros((n, self._basis.dimension)), self._basis)
        displacements[:, which] = self._basis.displacement(np.arange(1, n+1)*amt)
        coords = np.broadcast_to(self._coords, (n, ) + self._coords.shape)
        displacements = displacements.convert(self._coords.system)
        displaced = displacements + coords
        return displaced

    def visualize(self, step_size = .1, steps = (5, 5), which = 0, anim_opts = None, **plot_args):
        from McUtils.Plots import Animator

        if isinstance(steps, (int, np.integer)):
            steps = [steps, steps]

        left  = np.flip(self.displace( -step_size, steps[0], which), axis=0)
        right = self.displace(  step_size, steps[1], which)
        all_geoms = np.concatenate((left, np.broadcast_to(self._coords, (1, ) + self._coords.shape), right))

        figure, atoms, bonds = self._mol.plot(*all_geoms, objects = True, **plot_args)

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
                        my_stuff.extend(bb.plot(figure))

            return my_stuff

        if anim_opts is None:
            anim_opts = {}
        return Animator(figure, None, plot_method = animate, **anim_opts)

class NormalModeCoordiantes(CoordinateSystem):
    """A prettied up version of a Coordinerds CoordinateSystem object
    Has function for generating, though, too
    """
    def __init__(self, coeffs, basis = CartesianCoordinates3D, name = None, freqs = None):
        super().__init__(
            matrix = coeffs,
            name = name,
            basis = basis
        )
        self.freqs = freqs

    @classmethod
    def from_force_constants(cls, fcs, atoms = None, target_units = "Wavenumbers", **opts):
        if atoms is None:
            # intended for the internal coordinate case where it's hard to define masses in general
            weighted_fcs = fcs
        else:
            masses = np.array([
                AtomData[a, "Mass"] * UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
                    if isinstance(a, str) else
                a
                    for a in atoms
            ])
            masses = np.broadcast_to(masses, (len(masses), 3)).T.flatten()
            g = 1/np.sqrt(masses)
            g = np.outer(g, g)
            weighted_fcs = g * fcs

        freqs, modes = np.linalg.eigh(weighted_fcs)
        freqs = np.sign(freqs) * np.sqrt(np.abs(freqs))
        sorting = np.argsort(freqs)
        freqs = freqs[sorting]
        if target_units is not None and target_units == "Wavenumbers":
            freqs = freqs * UnitsData.convert("Hartrees", "Wavenumbers")
        modes = modes.T[sorting]

        return cls(modes, freqs = freqs, **opts)

#TODO: make it possible to just extract certain G-matrix elements without computing the whole thing
class GMatrix:
    """Represents Wilson's G Matrix between two coordinate systems"""
    def __init__(self, system1, system2, masses, **fd_opts):
        self.sys1 = system1
        self.sys2 = system2
        self.masses = masses,
        self.opts = fd_opts
        self._jacobian = None

    @property
    def array(self):
        """Returns the numpy array form of the G-matrix

        :return:
        :rtype: np.ndarray
        """
        return self.asarray()
    @property
    def jacobian(self):
        if self._jacobian is None:
            self._jacobian = self.sys1.jacobian(self.sys2, **self.opts)
        return self._jacobian
    def asarray(self, **opts):
        if len(opts) == 0:
            jacobian = self.jacobian
        else:
            opts = dict(self.opts, **opts)
            jacobian = self.sys1.jacobian(self.sys2, **opts)
        jj = np.matmul(jacobian, jacobian.T)
        if self.masses is not None:
            jj = self.masses * jj # mass weight by the rows
        return jj


