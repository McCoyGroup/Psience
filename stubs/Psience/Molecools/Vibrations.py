"""
Provides classes that manage molecular vibrations
"""
import numpy as np, scipy.linalg as slag
from McUtils.Coordinerds import CoordinateSystem, CoordinateSet
from McUtils.Data import AtomData, UnitsData
import McUtils.Numputils as nput
from .MoleculeInterface import AbstractMolecule
from .Transformations import MolecularTransformation
__all__ = ['MolecularVibrations', 'MolecularNormalModes']
__reload_hook__ = ['.MoleculeInterface', '.Transformations', '..Modes']

class MolecularVibrations:

    def __init__(self, molecule, basis, freqs=None, init=None):
        """
        Sets up a vibration for a Molecule object over the CoordinateSystem basis

        :param molecule:
        :type molecule: AbstractMolecule
        :param init:
        :type init: None | CoordinateSet
        :param basis:
        :type basis: MolecularNormalModes
        """
        ...

    def change_mol(self, mol):
        """
        **LLM Docstring**

        Rebind this set of vibrations to a different molecule, re-expressing the basis for the new molecule while keeping the same frequencies and reference coordinates.

        :param mol: the new molecule to associate with the vibrations
        :type mol: AbstractMolecule
        :return: a new `MolecularVibrations` built from `mol`, with `basis` changed via `basis.change_mol(mol)` and the same `freqs`/`init` as `self`
        :rtype: MolecularVibrations
        """
        ...

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
        ...

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
        ...

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
        ...

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
        ...

    @property
    def freqs(self):
        """
        **LLM Docstring**

        Frequencies associated with the vibrations. Returns the explicitly stored frequencies if present; otherwise falls back to `self.basis.freqs` when the basis defines that attribute.

        :return: the vibrational frequencies, or `None` if neither `self` nor `self.basis` has them
        :rtype: np.ndarray | None
        """
        ...

    @property
    def coords(self):
        """
        :return:
        :rtype: CoordinateSet
        """
        ...

    def __len__(self):
        """
        **LLM Docstring**

        Number of vibrational modes, taken from the number of columns of the basis's mode matrix.

        :return: the number of vibrational modes
        :rtype: int
        """
        ...

    def displace(self, displacements=None, amt=0.1, n=1, which=0):
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
        ...

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
        ...

    def _ipython_display_(self):
        """
        **LLM Docstring**

        IPython/Jupyter display hook: builds the interactive widget via `to_widget` and displays it.

        :return: the result of displaying the widget
        :rtype: object
        """
        ...

    def to_widget(self):
        """
        **LLM Docstring**

        Build (and cache) an interactive Jupyter widget for browsing through the vibrational modes: a menu to select which mode (`which`), paired with a live display that calls `self.visualize` for the selected mode.

        :return: the constructed `JHTML.Div` widget, or `None` if a widget was already built and cached in `self._widg`
        :rtype: JHTML.Div | None
        """
        ...

    def __getitem__(self, item):
        """
        Takes a slice of the modes

        :param item:
        :type item:
        :return:
        :rtype:
        """
        ...

    def embed(self, frame):
        """

        :param frame:
        :type frame: MolecularTransformation
        :return:
        :rtype:
        """
        ...

    def rescale(self, scaling):
        """
        Multiplies each mode by some scaling factor
        :param phases:
        :type phases:
        :return:
        :rtype:
        """
        ...

    def rotate(self, scaling):
        """
        Multiplies each mode by some scaling factor
        :param phases:
        :type phases:
        :return:
        :rtype:
        """
        ...

    def __repr__(self):
        """
        **LLM Docstring**

        Debug string representation showing the class name, basis, and molecule.

        :return: string of the form `ClassName(basis, molecule)`
        :rtype: str
        """
        ...

class MolecularNormalModes(CoordinateSystem):
    """
    A Coordinerds CoordinateSystem object that manages all of the data needed to
     work with normal mode coordinates + some convenience functions for generating and whatnot
    """
    name = 'MolecularNormalModes'

    def __init__(self, molecule, coeffs, name=None, freqs=None, internal=False, origin=None, basis=None, inverse=None):
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
        ...

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
        ...

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
        ...

    def change_mol(self, mol):
        """
        **LLM Docstring**

        Rebind these normal modes to a different molecule, keeping the mode matrix, name, frequencies, internal-coordinate flag, origin, basis, and inverse the same.

        :param mol: the new molecule to associate with the normal modes
        :type mol: AbstractMolecule
        :return: a new `MolecularNormalModes` for `mol` built from this object's `matrix`, `name`, `freqs`, `in_internals`, `_origin`, `basis`, and `inverse`
        :rtype: MolecularNormalModes
        """
        ...

    @property
    def coords_by_modes(self):
        """
        **LLM Docstring**

        Matrix mapping mode displacements back to coordinate displacements (i.e. the inverse transformation).

        :return: the stored inverse matrix, `self.inverse`
        :rtype: np.ndarray
        """
        ...

    @property
    def modes_by_coords(self):
        """
        **LLM Docstring**

        Matrix mapping coordinate displacements onto normal-mode displacements.

        :return: the stored mode matrix, `self.matrix`
        :rtype: np.ndarray
        """
        ...

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
        ...

    @property
    def origin(self):
        """
        **LLM Docstring**

        Reference geometry the normal modes are expanded about. If no explicit origin was stored, falls back to the molecule's internal coordinates (if `in_internals`) or Cartesian coordinates otherwise.

        :return: the origin coordinates
        :rtype: CoordinateSet
        """
        ...

    def embed(self, frame):
        """

        :param frame:
        :type frame: MolecularTransformation
        :return:
        :rtype:
        """
        ...

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
        ...

    def to_new_modes(self):
        """
        Converts to the new generalized normal modes

        :return:
        """
        ...

    @classmethod
    def from_new_modes(cls, mol, modes):
        """
        Converts to the new generalized normal modes

        :return:
        """
        ...

    @classmethod
    def from_force_constants(cls, molecule, fcs, *, atoms=None, masses=None, mass_units='AtomicMassUnits', inverse_mass_matrix=False, remove_transrot=True, dimensionless=False, mass_weighted=False, normalize=False, **opts):
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
        ...

    def __getitem__(self, item):
        """
        Takes a slice of the modes
        :param item:
        :type item:
        :return:
        :rtype:
        """
        ...