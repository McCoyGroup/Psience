"""
Provides classes that manage molecular vibrations
"""

import numpy as np, scipy.linalg as slag
from McUtils.Coordinerds import CoordinateSystem, CoordinateSet
from McUtils.Data import AtomData, UnitsData
import McUtils.Numputils as nput

from .MoleculeInterface import AbstractMolecule
from .Transformations import MolecularTransformation

__all__ = [
    "MolecularVibrations",
    "MolecularNormalModes",
]

__reload_hook__ = [".MoleculeInterface", ".Transformations", '..Modes']

class MolecularVibrations:

    def __init__(self,
                 molecule,
                 basis,
                 freqs=None,
                 init=None
                 ):
        """
        Sets up a vibration for a Molecule object over the CoordinateSystem basis

        :param molecule:
        :type molecule: AbstractMolecule
        :param init:
        :type init: None | CoordinateSet
        :param basis:
        :type basis: MolecularNormalModes
        """
        self._mol = molecule
        self._coords = init
        self._basis = basis
        self._freqs = freqs
        self._widg = None

    def change_mol(self, mol):
        """
        **LLM Docstring**

        Rebind this set of vibrations to a different molecule, re-expressing the basis for the new molecule while keeping the same frequencies and reference coordinates.

        :param mol: the new molecule to associate with the vibrations
        :type mol: AbstractMolecule
        :return: a new `MolecularVibrations` built from `mol`, with `basis` changed via `basis.change_mol(mol)` and the same `freqs`/`init` as `self`
        :rtype: MolecularVibrations
        """
        return type(self)(
            mol,
            self._basis.change_mol(mol),
            freqs=self._freqs,
            init=self._coords
        )

    @property
    def basis(self):
        """
        **LLM Docstring**

        Property getter/setter for the underlying vibrational basis (a `MolecularNormalModes` object) used to describe the vibrations.

        :param basis: (setter only) the new basis to store
        :type basis: MolecularNormalModes
        :return: (getter) the stored basis
        :rtype: MolecularNormalModes
        """
        return self._basis
    @basis.setter
    def basis(self, basis):
        """
        **LLM Docstring**

        Property getter/setter for the underlying vibrational basis (a `MolecularNormalModes` object) used to describe the vibrations.

        :param basis: (setter only) the new basis to store
        :type basis: MolecularNormalModes
        :return: (getter) the stored basis
        :rtype: MolecularNormalModes
        """
        self._basis = basis
    @property
    def molecule(self):
        """
        **LLM Docstring**

        Property getter/setter for the molecule associated with these vibrations. Setting it also propagates the new molecule to `self.basis.molecule`.

        :param mol: (setter only) the new molecule to associate
        :type mol: AbstractMolecule
        :return: (getter) the stored molecule
        :rtype: AbstractMolecule
        """
        return self._mol
    @molecule.setter
    def molecule(self, mol):
        """
        **LLM Docstring**

        Property getter/setter for the molecule associated with these vibrations. Setting it also propagates the new molecule to `self.basis.molecule`.

        :param mol: (setter only) the new molecule to associate
        :type mol: AbstractMolecule
        :return: (getter) the stored molecule
        :rtype: AbstractMolecule
        """
        self._mol = mol
        self._basis.molecule = mol

    @property
    def freqs(self):
        """
        **LLM Docstring**

        Frequencies associated with the vibrations. Returns the explicitly stored frequencies if present; otherwise falls back to `self.basis.freqs` when the basis defines that attribute.

        :return: the vibrational frequencies, or `None` if neither `self` nor `self.basis` has them
        :rtype: np.ndarray | None
        """
        freqs = self._freqs
        if freqs is None and hasattr(self.basis, "freqs"):
            freqs = self.basis.freqs
        return freqs

    @property
    def coords(self):
        """
        :return:
        :rtype: CoordinateSet
        """
        if self._coords is None:
            if self._basis.in_internals:
                return self._mol.internal_coordinates
            else:
                return self._mol.coords
        else:
            return self._coords

    def __len__(self):
        """
        **LLM Docstring**

        Number of vibrational modes, taken from the number of columns of the basis's mode matrix.

        :return: the number of vibrational modes
        :rtype: int
        """
        return self._basis.matrix.shape[1]

    def displace(self, displacements=None, amt=.1, n=1, which=0):
        """
        Displaces along the vibrational mode specified by `which`
        :param displacements:
        :type displacements:
        :param amt:
        :type amt:
        :param n:
        :type n:
        :param which:
        :type which:
        :return:
        :rtype:
        """
        from .CoordinateSystems import MolecularCartesianCoordinateSystem

        if displacements is None:
            displacements = self._basis.displacement(np.arange(1, n+1)*amt)
        displacements = np.asarray(displacements)
        if displacements.shape == ():
            displacements = np.array([displacements])

        coords = self.coords
        if displacements.ndim == 1:
            coord_shape = self._basis.coordinate_shape
            if coord_shape is None:
                coord_shape = (len(self.freqs),)

            disp_coords = CoordinateSet(np.zeros((n,) + coord_shape), self._basis)
            disp_coords[:, which] = displacements
        else:
            disp_coords = CoordinateSet(displacements, self._basis)
        displaced = disp_coords.convert(self.coords.system)

        return displaced

    def visualize(self, step_size=5, steps=(2, 2), which=0, anim_opts=None, mode='fast', **plot_args):
        """
        :param step_size:
        :type step_size:
        :param steps:
        :type steps:
        :param which:
        :type which:
        :param anim_opts:
        :type anim_opts:
        :param mode:
        :type mode:
        :param plot_args:
        :type plot_args:
        :return:
        :rtype:
        """
        from McUtils.Plots import Animator

        coords = self.coords
        if coords is None:
            raise ValueError("can't visualize vibrations if no starting structure is provided")

        if isinstance(steps, (int, np.integer)):
            steps = [steps, steps]

        left  = np.flip(self.displace(amt=-step_size, n=steps[0], which=which), axis=0)
        right = self.displace(amt=step_size, n=steps[1], which=which)

        if mode=='jupyter':

            all_geoms = np.concatenate((left, np.broadcast_to(coords, (1,) + coords.shape), right,
                                        np.flip(right, axis=0)[1:], np.broadcast_to(coords, (1,) + coords.shape), np.flip(left, axis=0)[:-1]
                                        ))

            from McUtils.Jupyter import MoleculeGraphics
            try:
                ats = self._mol.atoms
            except:
                raise Exception(self._mol)
            return MoleculeGraphics(self._mol.atoms,
                                    all_geoms,
                                    bonds=self._mol.bonds
                                    )
        else:

            all_geoms = np.concatenate((left, np.broadcast_to(coords, (1,) + coords.shape), right))
            figure, atoms, bonds = self._mol.plot(*all_geoms, objects=True, mode=mode, **plot_args)

            def animate(*args, frame=0, _atoms=atoms, _bonds=bonds, _figure=figure):
                """
                Animates a vibration. Currently uses Matplotlib and so can be _verrry_ slow
                """

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
            return Animator(figure, None, plot_method=animate, **anim_opts)

    def _ipython_display_(self):
        """
        **LLM Docstring**

        IPython/Jupyter display hook: builds the interactive widget via `to_widget` and displays it.

        :return: the result of displaying the widget
        :rtype: object
        """
        return self.to_widget().display()

    def to_widget(self):
        """
        **LLM Docstring**

        Build (and cache) an interactive Jupyter widget for browsing through the vibrational modes: a menu to select which mode (`which`), paired with a live display that calls `self.visualize` for the selected mode.

        :return: the constructed `JHTML.Div` widget, or `None` if a widget was already built and cached in `self._widg`
        :rtype: JHTML.Div | None
        """
        from McUtils.Jupyter import JHTML, MenuSelect, ButtonGroup, FunctionDisplay, VariableNamespace

        def get_which(which):
            """
            **LLM Docstring**

            Parse the mode-selector menu's string value into a zero-based mode index, defaulting to `0` if the value can't be parsed as an integer.

            :param which: the raw string value from the menu widget (expected to be a 1-based mode number as text)
            :type which: str
            :return: the zero-based mode index
            :rtype: int
            """
            try:
                which = int(which.strip())
            except:
                which = 0
            else:
                which = which - 1
            return which

        if self._widg is None:
            with VariableNamespace():
                return JHTML.Div(
                    MenuSelect('which', [str(x+1) for x in range(len(self))], value='1', menu_type=ButtonGroup),
                    FunctionDisplay(
                        lambda which='1',**kw:self.visualize(which=get_which(which), mode='jupyter'),
                        ['which']
                    )
                )

    def __getitem__(self, item):
        """
        Takes a slice of the modes

        :param item:
        :type item:
        :return:
        :rtype:
        """

        if nput.is_numeric(item):
            item = (item,)
        elif not nput.is_numeric(item[0]):
            item = tuple(item[0])

        m = self._basis[item]
        if self.freqs is not None:
            f = self.freqs[item,]
        else:
            f = None
        if isinstance(m, MolecularNormalModes):
            return type(self)(
                self._mol,
                m,
                freqs=f,
                init=self._coords
            )
        else:
            raise ValueError("{}: can't take element {} from {}".format(
                type(self).__name__,
                item,
                self
            ))

    def embed(self, frame):
        """

        :param frame:
        :type frame: MolecularTransformation
        :return:
        :rtype:
        """

        basis = self.basis.embed(frame)
        if self._coords is not None:
            raise ValueError("currently not reembedding coords not coming from a Molecule object")

        return type(self)(self.molecule, basis)

    def rescale(self, scaling):
        """
        Multiplies each mode by some scaling factor
        :param phases:
        :type phases:
        :return:
        :rtype:
        """

        return type(self)(self._mol,
                          self.basis.rescale(scaling),
                          freqs=self._freqs,
                          init=self._coords
                          )
    def rotate(self, scaling):
        """
        Multiplies each mode by some scaling factor
        :param phases:
        :type phases:
        :return:
        :rtype:
        """

        print(scaling)

        return type(self)(self._mol,
                          self.basis.rotate(scaling),
                          freqs=self._freqs,
                          init=self._coords
                          )

    def __repr__(self):
        """
        **LLM Docstring**

        Debug string representation showing the class name, basis, and molecule.

        :return: string of the form `ClassName(basis, molecule)`
        :rtype: str
        """
        return "{}({}, {})".format(type(self).__name__, self.basis, self.molecule)
class MolecularNormalModes(CoordinateSystem):
    """
    A Coordinerds CoordinateSystem object that manages all of the data needed to
     work with normal mode coordinates + some convenience functions for generating and whatnot
    """
    name="MolecularNormalModes"
    def __init__(self,
                 molecule, coeffs,
                 name=None, freqs=None,
                 internal=False, origin=None, basis=None, inverse=None
                 ):
        """
        :param molecule:
        :type molecule: AbstractMolecule
        :param coeffs:
        :type coeffs:
        :param name:
        :type name:
        :param freqs:
        :type freqs:
        :param internal:
        :type internal:
        :param origin:
        :type origin:
        :param basis:
        :type basis:
        :param inverse:
        :type inverse:
        """
        coeffs = np.asanyarray(coeffs)
        if freqs is None:
            freqs = np.diag(coeffs.T@coeffs)
        freqs = np.asanyarray(freqs)
        if inverse is None:
            if freqs is not None:
                inverse = coeffs.T/freqs[:, np.newaxis]
            else:
                inverse = None
        self._molecule = molecule
        self.in_internals = internal
        # if origin is None:
        #     origin = molecule.coords
        if basis is None:
            basis = molecule.coords.system #type: MolecularCartesianCoordinateSystem | MolecularZMatrixCoordinateSystem

        # norms = np.linalg.norm(
        #                 coeffs,
        #                 axis=0
        #             )
        # print(norms)
        # if norms[1] < .008:
        #     raise Exception('wat')
        super().__init__(
            matrix=coeffs,
            inverse=inverse,
            name=self.name if name is None else name,
            basis=basis,
            dimension=(len(freqs),),
            origin=origin
        )
        self.freqs = freqs
    @property
    def molecule(self):
        """
        **LLM Docstring**

        Property getter/setter for the molecule these normal modes belong to. Setting it also propagates the new molecule to `self.basis.molecule` before updating `self._molecule`.

        :param mol: (setter only) the new molecule to associate
        :type mol: AbstractMolecule
        :return: (getter) the stored molecule
        :rtype: AbstractMolecule
        """
        return self._molecule
    @molecule.setter
    def molecule(self, mol):
        """
        **LLM Docstring**

        Property getter/setter for the molecule these normal modes belong to. Setting it also propagates the new molecule to `self.basis.molecule` before updating `self._molecule`.

        :param mol: (setter only) the new molecule to associate
        :type mol: AbstractMolecule
        :return: (getter) the stored molecule
        :rtype: AbstractMolecule
        """
        # thought I wanted this to reset the origin and things..
        self.basis.molecule = mol
        self._molecule = mol

    def change_mol(self, mol):
        """
        **LLM Docstring**

        Rebind these normal modes to a different molecule, keeping the mode matrix, name, frequencies, internal-coordinate flag, origin, basis, and inverse the same.

        :param mol: the new molecule to associate with the normal modes
        :type mol: AbstractMolecule
        :return: a new `MolecularNormalModes` for `mol` built from this object's `matrix`, `name`, `freqs`, `in_internals`, `_origin`, `basis`, and `inverse`
        :rtype: MolecularNormalModes
        """
        return type(self)(
            mol,
            self.matrix,
            name=self.name,
            freqs=self.freqs,
            internal=self.in_internals,
            origin=self._origin,
            basis=self.basis,
            inverse=self.inverse
        )

    @property
    def coords_by_modes(self):
        """
        **LLM Docstring**

        Matrix mapping mode displacements back to coordinate displacements (i.e. the inverse transformation).

        :return: the stored inverse matrix, `self.inverse`
        :rtype: np.ndarray
        """
        return self.inverse
    @property
    def modes_by_coords(self):
        """
        **LLM Docstring**

        Matrix mapping coordinate displacements onto normal-mode displacements.

        :return: the stored mode matrix, `self.matrix`
        :rtype: np.ndarray
        """
        return self.matrix

    # also need a Cartesian equivalent of this
    def to_internals(self, intcrds=None, dYdR=None, dRdY=None):
        """
        **LLM Docstring**

        Intended to convert Cartesian-basis normal modes into an internal-coordinate representation using the supplied Jacobians. The method immediately raises `NotImplementedError` directing callers to use the newer `NormalModes` object instead, so the remaining body (computing `dQdR`/`dRdQ` and constructing a new internal-coordinate `MolecularNormalModes`) is dead code that never executes.

        :param intcrds: internal-coordinate system to convert into; if `None`, taken from `self.molecule.internal_coordinates`
        :type intcrds: CoordinateSet | None
        :param dYdR: Jacobian of mass-weighted Cartesians with respect to internal coordinates; computed from `intcrds`/`molecule` if not given (unreachable)
        :type dYdR: np.ndarray | None
        :param dRdY: Jacobian of internal coordinates with respect to mass-weighted Cartesians; computed if not given (unreachable)
        :type dRdY: np.ndarray | None
        :return: never returns; always raises
        :rtype: None
        :raises NotImplementedError: always, instructing the caller to use the new `NormalModes` object
        """
        raise NotImplementedError("use new `NormalModes` object")
        if self.in_internals:
            return self
        else:
            if intcrds is None:
                intcrds = self.molecule.internal_coordinates
                if intcrds is None:
                    raise ValueError("{}.{}: can't convert to internals when molecule {} has no internal coordinate specification".format(
                        type(self).__name__,
                        'to_internals',
                        self.molecule
                    ))
            if dRdY is None or dYdR is None:
                internals = intcrds.system
                ccoords = self.molecule.coords
                carts = ccoords.system
                ncrds = self.matrix.shape[0]
                dXdR = intcrds.jacobian(carts, 1).reshape(ncrds, ncrds)
                dRdX = ccoords.jacobian(internals, 1).reshape(ncrds, ncrds)
                masses = self.molecule.masses
                mass_conv = np.sqrt(np.broadcast_to(masses[:, np.newaxis], (3, len(masses))).flatten())
                dYdR = dXdR * mass_conv[np.newaxis]
                dRdY = dRdX / mass_conv[:, np.newaxis]

            # get the new normal modes
            dQdR = dYdR@self.matrix
            dRdQ = self.inverse@dRdY

            return type(self)(self.molecule, dQdR,
                              basis = intcrds.system, origin=intcrds, inverse=dRdQ,
                              internal=True, freqs=self.freqs
                              )
    @property
    def origin(self):
        """
        **LLM Docstring**

        Reference geometry the normal modes are expanded about. If no explicit origin was stored, falls back to the molecule's internal coordinates (if `in_internals`) or Cartesian coordinates otherwise.

        :return: the origin coordinates
        :rtype: CoordinateSet
        """
        if self._origin is None:
            if self.in_internals:
                return self.molecule.internal_coordinates
            else:
                return self.molecule.coords
        else:
            return self._origin

    def embed(self, frame):
        """

        :param frame:
        :type frame: MolecularTransformation
        :return:
        :rtype:
        """

        # raise NotImplementedError("Haven't finished doing this... :)")
        # import copy

        if self.in_internals:
            raise ValueError("Internal coordinate normals modes can't be re-embedded")

        tmat = frame.transformation_function.transform

        if not np.allclose((tmat @ tmat.T), np.eye(3)):
            raise ValueError("embedding matrix {} isn't a proper rotation".format(tmat))
        if  np.round(np.linalg.det(tmat), 7) != 1:
            raise ValueError("embedding matrix {} isn't a proper rotation; determinant is {}".format(tmat, np.linalg.det(tmat)))

        mat = self.matrix.reshape((self.matrix.shape[0]//3, 3, -1))
        mat_2 = np.moveaxis(
            np.tensordot(tmat, mat, axes=[1, 1]),
            0,
            1
        )
        mat = mat_2.reshape(self.matrix.shape)

        if self._inv is not None:
            inv = self.inverse.reshape((-1, self.matrix.shape[0] // 3, 3))
            inv = np.tensordot(inv, np.linalg.inv(tmat), axes=[2, 0])
            inv = inv.reshape(self.inverse.shape)
        else:
            inv = None

        # norms = np.linalg.norm(
        #         self.matrix,
        #         axis=0
        #     )
        # norms_2 = np.linalg.norm(
        #         mat,
        #         axis=0
        #     )
        # raise Exception(norms, norms_2)

        if self._origin is not None:
            # raise NotImplementedError("transforming explicit origin needs to be supported")
            raise Exception(self.origin)
            orig = frame(self.origin)
            # orig = np.tensordot(tmat, orig, axes=[1, 1])
        else:
            orig = None

        return type(self)(
            self.molecule,
            mat,
            inverse=inv,
            origin=orig,
            freqs=self.freqs
        )

    def insert(self, val, where):
        """
        Inserts values into the appropriate positions in the mode matrix

        :param val:
        :type val:
        :param where:
        :type where:
        :return:
        :rtype:
        """
        mat = self.matrix.reshape((self.matrix.shape[0] // 3, 3, self.matrix.shape[1]))
        mat = np.insert(mat, where, val, axis=0)
        mat = np.reshape(mat, (-1, mat.shape[2]))
        # import McUtils.Plots as plt
        # plt.ArrayPlot(self.matrix)
        # plt.ArrayPlot(mat).show()
        inv = self.inverse.reshape((self.inverse.shape[0], self.inverse.shape[1] // 3, 3))
        inv = np.insert(inv, where, val, axis=1)
        inv = np.reshape(inv, (inv.shape[0], -1))
        # import McUtils.Plots as plt
        # plt.ArrayPlot(self.inverse)
        # plt.ArrayPlot(inv).show()

        return type(self)(
            self.molecule, mat,
            name=self.name, freqs=self.freqs,
            internal=self.in_internals, origin=self._origin,
            basis=self.basis, inverse=inv
        )


    # def rescale(self, scaling_factors, in_place=False):
    #     """
    #     Rescales each mode in the expansion matrix
    #     by the passed `scaling_factors`
    #
    #     :param scaling_factors:
    #     :type scaling_factors:
    #     :return:
    #     :rtype:
    #     """
    #
    #     mat = self.matrix * scaling_factors[np.newaxis, :]
    #     inv = self.inverse / scaling_factors[:, np.newaxis]
    #
    #     return type(self)(
    #         self.molecule, mat,
    #         name=self.name, freqs=self.freqs,
    #         internal=self.in_internals, origin=self._origin,
    #         basis=self.basis, inverse=inv
    #     )

    def to_new_modes(self):
        """
        Converts to the new generalized normal modes

        :return:
        """
        from ..Modes import NormalModes
        if self.in_internals:
            basis = self.molecule.internal_coordinates.system
        else:
            basis = self.molecule.coords.system
        return NormalModes(
            basis,
            self.inverse.T,
            inverse=self.matrix.T,
            freqs=self.freqs,
            origin=self.origin,
            masses=self.molecule.atomic_masses,
            mass_weighted=False
        )
    @classmethod
    def from_new_modes(cls, mol, modes):
        """
        Converts to the new generalized normal modes

        :return:
        """
        modes = modes.remove_mass_weighting()
        return cls(
            mol,
            modes.inverse.T,
            inverse=modes.matrix.T,
            freqs=modes.freqs,
            origin=modes.origin
        )

    @classmethod
    def from_force_constants(cls,
                             molecule,
                             fcs,
                             *,
                             atoms=None,
                             masses=None,
                             mass_units="AtomicMassUnits",
                             inverse_mass_matrix=False,
                             remove_transrot=True,
                             dimensionless=False,
                             mass_weighted=False,
                             normalize=False,
                             **opts
                             ):
        """
        Generates normal modes from the specified force constants

        :param molecule:
        :type molecule: AbstractMolecule
        :param fcs: force constants array
        :type fcs: np.ndarray
        :param atoms: atom list
        :type atoms: Iterable[str]
        :param masses: mass list
        :type masses: Iterable[float]
        :param mass_units: units for the masses...not clear if this is useful or a distraction
        :type mass_units: str
        :param inverse_mass_matrix: whether or not we have G or G^-1 (default: `False`)
        :type inverse_mass_matrix: bool
        :param remove_transrot: whether or not to remove the translations and rotations (default: `True`)
        :type remove_transrot: bool
        :param normalize: whether or not to normalize the modes (default: `True`)
        :type normalize: bool
        :param opts:
        :type opts:
        :return:
        :rtype: MolecularNormalModes
        """

        from ..Modes import NormalModes

        if atoms is not None and masses is None:
            masses = atoms

        freqs, modes, inv = NormalModes.get_normal_modes(
            fcs,
            masses,
            # mass_units="AtomicMassUnits",
            remove_transrot=remove_transrot,
            dimensionless=dimensionless,
            mass_weighted=mass_weighted
        )

        return cls(molecule, inv.T, inverse=modes.T, freqs=freqs, **opts)


    def __getitem__(self, item):
        """
        Takes a slice of the modes
        :param item:
        :type item:
        :return:
        :rtype:
        """

        if nput.is_numeric(item):
            item = (item,)
        elif not nput.is_numeric(item[0]):
            item = tuple(item[0])

        sub_modes = self.matrix[:, item]
        inv = self._inv
        if inv is not None:
            # i0 = inv
            inv = inv[item, :]
            # raise Exception([item, sub_modes.shape, inv.shape, i0.shape])
        freq = self.freqs[item,]
        return type(self)(
            self.molecule,
            sub_modes,
            name=self.name,
            freqs=freq,
            internal=self.in_internals,
            origin=self._origin,
            basis=self.basis,
            inverse=inv
        )