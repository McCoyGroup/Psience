"""
A collection of methods used in computing molecule properties
"""
import itertools
import numpy as np, scipy.sparse as sp, itertools as ip, os, abc
from collections import namedtuple
import McUtils.Numputils as nput
from McUtils.Coordinerds import CoordinateSet
from McUtils.ExternalPrograms import OpenBabelInterface
from McUtils.GaussianInterface import GaussianFChkReader, GaussianFChkReaderException
from McUtils.Data import AtomData, UnitsData, BondData
from McUtils.Zachary import FiniteDifferenceDerivative, TensorDerivativeConverter
from .MoleculeInterface import AbstractMolecule
from .Transformations import MolecularTransformation
from .Vibrations import MolecularVibrations, MolecularNormalModes
__all__ = ['StructuralProperties', 'BondingProperties', 'MolecularProperties', 'MolecularPropertyError', 'OpenBabelMolManager', 'DipoleSurfaceManager', 'PotentialSurfaceManager', 'NormalModesManager']
__reload_hook__ = ['.MoleculeInterface', '.Vibrations', '.Transformations']

class MolecularPropertyError(Exception):
    """
    General error class for MolecularProperties
    """

class StructuralProperties:
    """
    The set of molecular properties that depend on its coordinates/configuration.
    Slowly trying to move code out of this and into numputils/Hamiltonian/Evaluator
    """

    @classmethod
    def get_prop_mass_weighted_coords(cls, coords, masses):
        """Gets the mass-weighted coordinates for the system

        :param coords:
        :type coords: CoordinateSet
        :param masses:
        :type masses:
        :return:
        :rtype:
        """
        ...

    @classmethod
    def get_prop_center_of_mass(cls, coords, masses):
        """Gets the center of mass for the coordinates

        :param coords:
        :type coords: CoordinateSet
        :param masses:
        :type masses:
        :return:
        :rtype:
        """
        ...

    @classmethod
    def get_prop_inertia_tensors(cls, coords, masses):
        """
        Computes the moment of intertia tensors for the walkers with coordinates coords (assumes all have the same masses)

        :param coords:
        :type coords: CoordinateSet
        :param masses:
        :type masses: np.ndarray
        :return:
        :rtype:
        """
        ...

    @classmethod
    def get_prop_inertial_frame_derivatives(cls, crds, mass, sel=None):
        """
        **LLM Docstring**

        Compute the first and second derivatives of the (mass-weighted) inertia tensor with respect to mass-weighted Cartesian displacements, using closed-form tensor expressions rather than finite differences. Handles batched inputs (multiple geometries at once) and can restrict to a subset of atoms via `sel`.

        :param crds: Cartesian coordinates of the atoms, with an optional leading batch dimension
        :type crds: np.ndarray
        :param mass: atomic masses (only positive entries are treated as real atoms), with an optional leading batch dimension
        :type mass: np.ndarray
        :param sel: optional subset of atom indices to restrict the calculation to; intersected with the real (mass > 0) atoms
        :type sel: np.ndarray | None
        :return: `[I0Y, I0YY]`, the first derivative tensor (shape `(..., 3*nAt, 3, 3)`) and second derivative tensor (shape `(..., 3*nAt, 3*nAt, 3, 3)`) of the inertia tensor with respect to mass-weighted Cartesian displacements
        :rtype: list[np.ndarray]
        """
        ...

    @classmethod
    def get_prop_moments_of_inertia(cls, coords, masses):
        """
        Computes the moment of inertia tensor for the walkers with coordinates coords (assumes all have the same masses)

        :param coords:
        :type coords: CoordinateSet
        :param masses:
        :type masses: np.ndarray
        :return:
        :rtype:
        """
        ...

    @classmethod
    def get_prop_principle_axis_rotation(cls, coords, masses, sel=None, inverse=False, shift_sel=True):
        """
        Generates the principle axis transformation for a set of coordinates and positions

        :param coords:
        :type coords: CoordinateSet
        :param masses:
        :type masses: np.ndarray
        :return:
        :rtype: MolecularTransformation | List[MolecularTransformation]
        """
        ...

    @classmethod
    def get_principle_axis_embedded_coords(cls, coords, masses, sel=None):
        """
        Returns coordinate embedded in the principle axis frame

        :param coords:
        :type coords:
        :param masses:
        :type masses:
        :return:
        :rtype:
        """
        ...

    @classmethod
    def get_prop_principle_axis_data(cls, coords, masses):
        """
        Generates the principle axis transformation for a set of coordinates and positions

        :param coords:
        :type coords: CoordinateSet
        :param masses:
        :type masses: np.ndarray
        :return:
        :rtype: MolecularTransformation | List[MolecularTransformation]
        """
        ...

    @classmethod
    def get_prop_translation_rotation_eigenvectors(cls, coords, masses):
        """
        Returns the eigenvectors corresponding to translations and rotations
        in the system

        :param coords:
        :type coords:
        :param masses:
        :type masses:
        :return:
        :rtype:
        """
        ...
    planar_ref_tolerance = 1e-06

    @classmethod
    def get_eckart_rotations(cls, masses, ref, coords, sel=None, in_paf=False, planar_ref_tolerance=None, proper_rotation=False):
        """
        Generates the Eckart rotation that will align ref and coords, assuming initially that `ref` and `coords` are
        in the principle axis frame

        :param masses:
        :type masses:
        :param ref:
        :type ref:
        :param coords:
        :type coords: np.ndarray
        :return:
        :rtype:
        """
        ...
    EmbeddingData = namedtuple('PrincipleAxisData', ['coords', 'com', 'axes'])
    EckartData = namedtuple('EckartData', ['rotations', 'reference_data', 'coord_data'])

    @classmethod
    def get_eckart_embedding_data(cls, masses, ref, coords, sel=None, in_paf=False, planar_ref_tolerance=None, proper_rotation=False):
        """
        Embeds a set of coordinates in the reference frame

        :param masses:
        :type masses: np.ndarray
        :param ref:
        :type ref: CoordinateSet
        :param coords:
        :type coords: CoordinateSet
        :return:
        :rtype:
        """
        ...

    @classmethod
    def get_prop_eckart_transformation(cls, masses, ref, coords, sel=None, inverse=False, reset_com=False, in_paf=False, planar_ref_tolerance=None, proper_rotation=False):
        """
        Computes Eckart transformations for a set of coordinates

        :param masses:
        :type masses: np.ndarray
        :param ref:
        :type ref: CoordinateSet
        :param coords:
        :type coords: CoordinateSet
        :return:
        :rtype: MolecularTransformation | List[MolecularTransformation]
        """
        ...

    @classmethod
    def get_eckart_embedded_coords(cls, masses, ref, coords, reset_com=True, in_paf=False, sel=None, planar_ref_tolerance=None, proper_rotation=False):
        """
        Embeds a set of coordinates in the reference frame

        :param masses:
        :type masses: np.ndarray
        :param ref:
        :type ref: CoordinateSet
        :param coords:
        :type coords: CoordinateSet
        :return:
        :rtype:
        """
        ...

    @classmethod
    def get_prop_g_matrix(cls, masses, coords, internal_coords):
        """
        Gets the molecular g-matrix
        :param masses:
        :type masses: np.ndarray
        :param coords:
        :type coords: CoordinateSet
        :param internal_coords:
        :type internal_coords: CoordinateSet
        :return:
        :rtype:
        """
        ...

    @classmethod
    def get_prop_coriolis_constants(cls, carts, modes, masses):
        """
        **LLM Docstring**

        Compute the Coriolis zeta constants for a set of normal modes, by projecting the mode displacement vectors into the molecule's inertial (principal-axis) frame and forming the antisymmetric combination that gives the coupling between each pair of modes about each principal axis.

        :param carts: Cartesian coordinates of the atoms, with an optional leading batch dimension
        :type carts: np.ndarray
        :param modes: the (mass-weighted) normal-mode displacement vectors
        :type modes: np.ndarray
        :param masses: the atomic masses
        :type masses: np.ndarray
        :return: a `(zeta, (mom_i, eigs))` pair, where `zeta` has shape `(..., 3, nmodes, nmodes)` giving the zeta constant about each of the 3 principal axes for every mode pair, `mom_i` is the moments of inertia, and `eigs` is the principal-axis frame
        :rtype: tuple[np.ndarray, tuple[np.ndarray, np.ndarray]]
        """
        ...

class BondingProperties:
    """
    The set of properties that depend only on bonding
    and that kind of format
    """

    @classmethod
    def get_prop_adjacency_matrix(cls, atoms, bonds):
        """
        Returns the adjacency matrix for the molecule

        :param bonds:
        :type bonds: Iterable[int]
        :return:
        :rtype:
        """
        ...

    @classmethod
    def get_prop_connectivity(cls, atoms, bonds):
        """
        Returns the adjacency matrix for the molecule

        :param bonds:
        :type bonds: Iterable[int]
        :return:
        :rtype:
        """
        ...

    @classmethod
    def get_prop_fragments(cls, atoms, bonds):
        """
        Returns the fragments for the molecule

        :param bonds:
        :type bonds: Iterable[int]
        :return:
        :rtype:
        """
        ...

    @classmethod
    def get_prop_zmat_ordering(cls, atoms, bonds):
        """
        Gets a guessed Z-matrix ordering for the molecule with connectivity defined by bonds based on the following:
            1. Fragments are separated out
            2. The atom with the highest degree of connectivity in each fragment is chosen as the fragment "label"
            3. Fragments are ordered by connectivity of the label from high to low
            4. Fragment labels reference each other with:
                a) the second label on the x-axis
                b) the 3rd in the xy-plane
                c) all others relative to the first three
            5. All other atoms are sorted first by fragment label, then by connection to the fragment label, and then by connectivity
            6. Atoms reference each other based on the following:
                a) if the atom has one bond:
                    i)   the atom it is bound to
                    ii)  the lowest-connectivity atom that one is bound to
                    iii) the second-lowest-connectivity atom OR the next fragment label
                b) if the atom has two bonds:
                    i)   the highest-connectivity atom it is bound to
                    ii)  the lowest-connectivity atom it is bound to
                    iii) the lowest-connectivity atom (i) is bound to
                c) if the atom has three bonds:
                    i)   the highest-connectivity atom it is bound to
                    ii)  the lowest-connectivity atom it is bound to
                    iii) the second-highest connectivity atom it is bound to
              if any of these atoms do not exist, the fragment labels will be used in their place

        :param bonds:
        :type bonds:
        :return:
        :rtype:
        """
        ...

    @classmethod
    def get_prop_guessed_bonds(cls, mol, tol=1.05, guess_type=True, covalent_radius_scaling=1.1):
        """
        Guesses the bonds for the molecule by finding the ones that are less than some percentage of a single bond for that
        pair of elements

        :return:
        :rtype:
        """
        ...

class MolecularProperties:
    """
    An object whose sole purpose in life is to get molecular properties
    A property should be implemented in two parts:
        1) a classmethod called get_prop_<prop name> that takes raw inputs and uses them to compute a property
        2) a classmethod called <prop name> that extracts the property from a passed molecule

    All properties should be appropriately vectorized to work on a single configuration or a set of configurations
    """

    @classmethod
    def mass_weighted_coords(cls, mol):
        """
        Computes the moments of inertia

        :param mol:
        :type mol: AbstractMolecule
        :return:
        :rtype:
        """
        ...

    @classmethod
    def g_matrix(cls, mol):
        """
        :param mol:
        :type mol: AbstractMolecule
        :return:
        :rtype:
        """
        ...

    @classmethod
    def center_of_mass(cls, mol):
        """
        Computes the moments of inertia

        :param mol:
        :type mol: AbstractMolecule
        :return:
        :rtype:
        """
        ...

    @classmethod
    def inertia_tensor(cls, mol):
        """
        Computes the inertia tensors for the stored geometries

        :param mol:
        :type mol: AbstractMolecule
        :return:
        :rtype:
        """
        ...

    @classmethod
    def moments_of_inertia(cls, mol):
        """
        Computes the moments of inertia

        :param mol:
        :type mol: AbstractMolecule
        :return:
        :rtype:
        """
        ...

    @classmethod
    def principle_axis_data(cls, mol, sel=None):
        """
        Generates the center of masses and inertial axes

        :param mol:
        :type mol: AbstractMolecule
        :return:
        :rtype:
        """
        ...

    @classmethod
    def principle_axis_transformation(cls, mol, sel=None, inverse=False):
        """
        Generates the principle axis transformation for a Molecule

        :param mol:
        :type mol: AbstractMolecule
        :return:
        :rtype:
        """
        ...

    @classmethod
    def eckart_embedding_data(cls, mol, coords, sel=None, in_paf=False, planar_ref_tolerance=None, proper_rotation=False):
        """

        :param mol:
        :type mol: AbstractMolecule
        :param coords:
        :type coords:
        :param sel:
        :type sel:
        :return:
        :rtype:
        """
        ...

    @classmethod
    def eckart_transformation(cls, mol, ref_mol, sel=None, inverse=False, planar_ref_tolerance=None, reset_com=True, proper_rotation=False):
        """

        :param ref_mol: reference geometry
        :type ref_mol: AbstractMolecule
        :param mol: molecules to get Eckart embeddings for
        :type mol: AbstractMolecule
        :param sel: coordinate selection to use when doing the Eckart stuff
        :type mol:
        :return:
        :rtype:
        """
        ...

    @classmethod
    def eckart_embedded_coords(cls, mol, coords, sel=None, in_paf=False, reset_com=True, planar_ref_tolerance=None, proper_rotation=False):
        """

        :param mol:
        :type mol: AbstractMolecule
        :param coords:
        :type coords:
        :param sel:
        :type sel:
        :return:
        :rtype:
        """
        ...

    @classmethod
    def coriolis_constants(cls, mol):
        """
        **LLM Docstring**

        Compute the Coriolis zeta constants for a molecule's own normal modes, using its Cartesian coordinates, normal-mode matrix, and atomic masses.

        :param mol: the molecule to compute Coriolis constants for
        :type mol: AbstractMolecule
        :return: the zeta constants, as the first element of `StructuralProperties.get_prop_coriolis_constants`
        :rtype: np.ndarray
        """
        ...

    @classmethod
    def translation_rotation_eigenvectors(cls, mol, sel=None):
        """

        :param mol: molecules to get eigenvectors for
        :type mol: AbstractMolecule
        :param sel: coordinate selection to use when doing the rotation/translation calculations
        :type mol:
        :return:
        :rtype:
        """
        ...

    @classmethod
    def fragments(cls, mol):
        """

        :param mol:
        :type mol: AbstractMolecule
        :return:
        :rtype:
        """
        ...

    @classmethod
    def fragment_indices(cls, mol):
        """

        :param mol:
        :type mol: AbstractMolecule
        :return:
        :rtype:
        """
        ...

    @classmethod
    def edge_graph(cls, mol):
        """

        :param mol:
        :type mol: AbstractMolecule
        :return:
        :rtype:
        """
        ...

    @classmethod
    def guessed_bonds(cls, mol, tol=1.05, guess_type=True):
        """
        Guesses the bonds for the molecule by finding the ones that are less than some percentage of a single bond for that
        pair of elements

        :return:
        :rtype:
        """
        ...

    @classmethod
    def get_prop_chemical_formula(cls, atoms):
        """

        :param atoms:
        :type atoms: Tuple[str]
        :return:
        :rtype:
        """
        ...

    @classmethod
    def chemical_formula(cls, mol):
        """

        :param mol:
        :type mol: AbstractMolecule
        :return:
        :rtype:
        """
        ...

class StringFormatHandler(metaclass=abc.ABCMeta):
    """
    A base class to handle converting to/from string formats
    mostly just here to implement SDF so OpenBabel can read off that
    """

    def __init__(self, mol):
        """
        **LLM Docstring**

        Store the molecule this format handler will convert to/from a string representation.

        :param mol: the molecule to convert
        :type mol: AbstractMolecule
        :return: None
        :rtype: None
        """
        ...

    @abc.abstractmethod
    def get_string(self, atoms, coords, bonds, meta):
        """
        Converts the molecular info to string format

        :param atoms:
        :type atoms: tuple[str]
        :param coords:
        :type coords: CoordinateSet
        :param bonds:
        :type bonds: tuple[tuple]
        :param meta:
        :type meta: dict
        :return:
        :rtype: str
        """
        ...

    def convert(self):
        """
        **LLM Docstring**

        Convert the stored molecule to its string representation by extracting its atoms, coordinates, bonds, and metadata and passing them to the (subclass-defined) `get_string`.

        :return: the string representation of the molecule
        :rtype: str
        """
        ...

    @classmethod
    @abc.abstractmethod
    def parse_string(cls, str):
        """

        :param str:
        :type str:
        :return: (atoms, coords, bonds, meta)
        :rtype: (list, np.ndarray, list, dict)
        """
        ...

    @classmethod
    def parse(cls, string):
        """
        **LLM Docstring**

        Parse a string representation back into a `Molecule` object, using the (subclass-defined) `parse_string` to extract the atoms, coordinates, bonds, and metadata.

        :param string: the string to parse
        :type string: str
        :return: the reconstructed molecule
        :rtype: Molecule
        """
        ...

class SDFFormatHandler(StringFormatHandler):
    misc_useless_structural_data_header = ' 0     0  0  0  0  0  0999 V2000'
    program = 'Psience.Molecools'

    def __init__(self, mol):
        """
        :param mol:
        :type mol: AbstractMolecule
        :param program:
        :type program:
        :param comment:
        :type comment:
        """
        ...

    def convert_header(self, comment=None):
        """
        **LLM Docstring**

        Build the 3-line SDF header block: the molecule's name, the generating program name, and a comment line.

        :param comment: extra text appended to the comment line, after `self.comment`
        :type comment: str | None
        :return: the newline-joined header block
        :rtype: str
        """
        ...

    def convert_counts_line(self, atoms, bonds):
        """
        **LLM Docstring**

        Build the SDF "counts" line giving the number of atoms and bonds, plus the fixed structural-data suffix required by the SDF format.

        :param atoms: the atoms of the molecule
        :type atoms: tuple[str]
        :param bonds: the bonds of the molecule
        :type bonds: tuple[tuple]
        :return: the formatted counts line
        :rtype: str
        """
        ...

    def convert_coordinate_block(self, atoms, coords):
        """
        **LLM Docstring**

        Build the SDF atom block: one fixed-width line per atom giving its Cartesian coordinates and element symbol (plus the required trailing zero fields).

        :param atoms: the atoms of the molecule
        :type atoms: tuple[str]
        :param coords: the Cartesian coordinates, one row per atom
        :type coords: CoordinateSet
        :return: the newline-joined atom block
        :rtype: str
        """
        ...

    def convert_bond_block(self, bonds):
        """
        **LLM Docstring**

        Build the SDF bond block: one fixed-width line per bond giving the 1-indexed atom pair and bond order (defaulting to a single bond if none is specified).

        :param bonds: the bonds, each as an `(atom1, atom2[, order])` tuple of 0-indexed atom indices
        :type bonds: tuple[tuple]
        :return: the newline-joined bond block
        :rtype: str
        """
        ...

    def convert_metadata(self, meta):
        """
        **LLM Docstring**

        Build the SDF metadata block, formatting each key/value pair as an SDF data-item tag (`>  <key>` followed by its value).

        :param meta: the metadata to serialize
        :type meta: dict
        :return: the double-newline-joined metadata block
        :rtype: str
        """
        ...

    def get_single_structure_string(self, atoms, coords, bonds, meta):
        """
        **LLM Docstring**

        Assemble the full SDF record for a single structure by combining the header, counts line, atom block, bond block, the `M  END` terminator, metadata block, and the `$$$$` record separator.

        :param atoms: the atoms of the molecule
        :type atoms: tuple[str]
        :param coords: the Cartesian coordinates
        :type coords: CoordinateSet
        :param bonds: the bonds of the molecule
        :type bonds: tuple[tuple]
        :param meta: the molecule's metadata
        :type meta: dict
        :return: the complete SDF record for this structure
        :rtype: str
        """
        ...

    def get_string(self, atoms, coords, bonds, meta):
        """
        Converts the molecular info to string format

        :param atoms:
        :type atoms: tuple[str]
        :param coords:
        :type coords: CoordinateSet
        :param bonds:
        :type bonds: tuple[tuple]
        :param meta:
        :type meta: dict
        :return:
        :rtype: str
        """
        ...

    @classmethod
    def parse_string(cls, str):
        """

        :param str:
        :type str:
        :return: (atoms, coords, bonds, meta)
        :rtype: (list, np.ndarray, list, dict)
        """
        ...

class PropertyManager(metaclass=abc.ABCMeta):
    """
    A utility base class so to make it easier to have a unified way to
    handled derived properties
    """
    name = None

    def __init__(self, mol):
        """
        **LLM Docstring**

        Store the molecule this property manager is attached to.

        :param mol: the molecule whose derived property this manager handles
        :type mol: Molecule
        :return: None
        :rtype: None
        """
        ...

    def set_molecule(self, mol):
        """
        **LLM Docstring**

        Rebind this manager to a different molecule.

        :param mol: the new molecule to associate
        :type mol: AbstractMolecule
        :return: None
        :rtype: None
        """
        ...

    @classmethod
    @abc.abstractmethod
    def from_data(cls, mol, data):
        """
        **LLM Docstring**

        Abstract constructor for building a property manager directly from raw data. Not implemented on the base class; subclasses must override it.

        :param mol: the molecule the manager will be attached to
        :type mol: AbstractMolecule
        :param data: the raw data to build the manager's property from
        :type data: object
        :return: never returns on the base class
        :rtype: PropertyManager
        :raises NotImplementedError: always, on the base class
        """
        ...

    @abc.abstractmethod
    def load(self):
        """
        Loads in the values

        :return:
        :rtype:
        """
        ...

    @abc.abstractmethod
    def update(self, val):
        """
        Updates the held values

        :param val:
        :type val:
        :return:
        :rtype:
        """
        ...

    @abc.abstractmethod
    def apply_transformation(self, transf):
        """
        Applies a transformation to the held values

        :param transf:
        :type transf: MolecularTransformation
        :return:
        :rtype: PropertyManager
        """
        ...

    @abc.abstractmethod
    def insert_atoms(self, atoms, coords, where):
        """
        Handles the insertion of new atoms into the structure

        :param atoms:
        :type atoms: tuple[str]
        :param coords:
        :type coords: CoordinateSet
        :param where:
        :type where: tuple[int]
        :return:
        :rtype:
        """
        ...

    @staticmethod
    def _insert_derivative_zeros(derivs, where, extra_shape=0):
        """
        Inserts zeros into derivative terms to account for
        the insertion of dummy atoms into a system

        :param derivs:
        :type derivs:
        :param where:
        :type where:
        :param extra_shape:
        :type extra_shape:
        :return:
        :rtype:
        """
        ...

    @abc.abstractmethod
    def delete_atoms(self, where):
        """
        Handles the deletion from the structure

        :param atoms:
        :type atoms: tuple[str]
        :param coords:
        :type coords: CoordinateSet
        :param where:
        :type where: tuple[int]
        :return:
        :rtype:
        """
        ...

    @staticmethod
    def _transform_derivatives(derivs, transf):
        """
        Handles the transformation of derivative tensors
        across coordinate systems

        :param derivs: derivative tensors starting at order 0
        :type derivs: tuple[np.ndarray]
        :param transf:
        :type transf: MolecularTransformation
        :return:
        :rtype: tuple[np.ndarray]
        """
        ...

    def __repr__(self):
        """
        **LLM Docstring**

        Debug string representation using the manager's `name` (or its class name if `name` is unset) and the associated molecule.

        :return: string of the form `Name(molecule)`
        :rtype: str
        """
        ...

    def copy(self):
        """
        **LLM Docstring**

        Make a shallow copy of this manager.

        :return: a shallow copy of `self`
        :rtype: PropertyManager
        """
        ...

class OpenBabelMolManager(PropertyManager):
    name = 'OBMol'

    def __init__(self, mol, obmol=None):
        """
        **LLM Docstring**

        Wrap an OpenBabel molecule object for a given molecule, optionally supplying an already-constructed OpenBabel molecule.

        :param mol: the molecule this manager wraps an OpenBabel representation for
        :type mol: AbstractMolecule
        :param obmol: an existing OpenBabel `OBMol` object to use, if already available
        :type obmol: object | None
        :return: None
        :rtype: None
        """
        ...

    def load(self):
        """
        **LLM Docstring**

        Lazily build (and cache) the underlying OpenBabel `OBMol` object by converting the molecule to an SDF string via `SDFFormatHandler` and reading it in with `pybel`.

        :return: the OpenBabel `OBMol` object
        :rtype: object
        """
        ...

    def update(self, val):
        """
        **LLM Docstring**

        Replace the cached OpenBabel molecule object.

        :param val: the new OpenBabel molecule object
        :type val: object
        :return: None
        :rtype: None
        """
        ...

    def apply_transformation(self, transf):
        """
        Applies a transformation to the held values

        :param transf:
        :type transf: MolecularTransformation
        :return:
        :rtype: PropertyManager
        """
        ...

    def insert_atoms(self, atoms, coords, where):
        """
        Handles the insertion of new atoms into the structure

        :param atoms:
        :type atoms: tuple[str]
        :param coords:
        :type coords: CoordinateSet
        :param where:
        :type where: tuple[int]
        :return:
        :rtype:
        """
        ...

    def delete_atoms(self, where):
        """
        Handles the deletion from the structure

        :param atoms:
        :type atoms: tuple[str]
        :param coords:
        :type coords: CoordinateSet
        :param where:
        :type where: tuple[int]
        :return:
        :rtype:
        """
        ...

class DipoleSurfaceManager(PropertyManager):
    name = 'DipoleSurface'

    def __init__(self, mol, surface=None, derivatives=None, polarizability_derivatives=None):
        """
        **LLM Docstring**

        Set up the dipole-surface manager, either copying the surface/derivatives/analytic-flag from an existing surface-like object, or storing the supplied `surface`/`derivatives`/`polarizability_derivatives` directly (unpacking `derivatives` into numerical and analytic parts if given as a dict).

        :param mol: the molecule this manager is attached to
        :type mol: AbstractMolecule
        :param surface: an existing dipole surface, or another object exposing `_surf`/`_derivs`/`_analytic_derivatives` to copy from
        :type surface: object | None
        :param derivatives: raw dipole derivative tensors, or a dict with `'numerical'`/`'analytic'` keys
        :type derivatives: list[np.ndarray] | dict | None
        :param polarizability_derivatives: raw polarizability derivative tensors
        :type polarizability_derivatives: list[np.ndarray] | None
        :return: None
        :rtype: None
        """
        ...

    @classmethod
    def from_data(cls, mol, data):
        """
        **LLM Docstring**

        Abstract constructor for building a dipole-surface manager directly from raw data. Not implemented.

        :param mol: the molecule the manager will be attached to
        :type mol: AbstractMolecule
        :param data: the raw data to build from
        :type data: object
        :return: never returns
        :rtype: DipoleSurfaceManager
        :raises NotImplementedError: always
        """
        ...

    @property
    def surface(self):
        """
        **LLM Docstring**

        The dipole surface object, lazily loaded via `load_dipole_surface` if not already set.

        :return: the dipole surface
        :rtype: object
        """
        ...

    @property
    def numerical_derivatives(self):
        """
        **LLM Docstring**

        The numerically-computed dipole derivatives, lazily loaded (together with the analytic derivatives, if provided by the same source) via `load_dipole_derivatives` if neither is already cached.

        :return: the numerical dipole derivatives, or `None` if unavailable
        :rtype: tuple | None
        """
        ...

    def get_derivatives(self, quiet=False):
        """
        **LLM Docstring**

        Fetch the analytic dipole derivatives, lazily loading them (together with the numerical derivatives, if bundled) via `load_dipole_derivatives` if not already cached.

        :param quiet: if `True`, suppresses errors from the underlying loader when the derivatives can't be found
        :type quiet: bool
        :return: the analytic dipole derivatives
        :rtype: tuple
        """
        ...

    @property
    def derivatives(self):
        """
        **LLM Docstring**

        Property getter/setter for the analytic dipole derivatives. The getter delegates to `get_derivatives()`; the setter accepts either a raw derivative tuple or a dict with `'numerical'`/`'analytic'` keys.

        :param derivatives: (setter only) the new derivatives to store
        :type derivatives: tuple | dict
        :return: (getter) the analytic dipole derivatives
        :rtype: tuple
        """
        ...

    @derivatives.setter
    def derivatives(self, derivatives):
        """
        **LLM Docstring**

        Property getter/setter for the analytic dipole derivatives. The getter delegates to `get_derivatives()`; the setter accepts either a raw derivative tuple or a dict with `'numerical'`/`'analytic'` keys.

        :param derivatives: (setter only) the new derivatives to store
        :type derivatives: tuple | dict
        :return: (getter) the analytic dipole derivatives
        :rtype: tuple
        """
        ...

    def get_pol_derivatives(self, quiet=False):
        """
        **LLM Docstring**

        Fetch the polarizability derivatives, lazily loading them via `load_polarizability_derivatives` if not already cached.

        :param quiet: if `True`, suppresses errors from the underlying loader when the derivatives can't be found
        :type quiet: bool
        :return: the polarizability derivatives
        :rtype: list
        """
        ...

    @property
    def polarizability_derivatives(self):
        """
        **LLM Docstring**

        Property getter/setter for the polarizability derivatives. The getter delegates to `get_pol_derivatives()`.

        :param derivatives: (setter only) the new polarizability derivatives to store
        :type derivatives: list
        :return: (getter) the polarizability derivatives
        :rtype: list
        """
        ...

    @polarizability_derivatives.setter
    def polarizability_derivatives(self, derivatives):
        """
        **LLM Docstring**

        Property getter/setter for the polarizability derivatives. The getter delegates to `get_pol_derivatives()`.

        :param derivatives: (setter only) the new polarizability derivatives to store
        :type derivatives: list
        :return: (getter) the polarizability derivatives
        :rtype: list
        """
        ...

    def load(self):
        """
        **LLM Docstring**

        Return whichever representation of the dipole is available: the surface object if one is set, otherwise the derivative expansion.

        :return: the dipole surface or its derivatives
        :rtype: object
        """
        ...

    def update(self, val):
        """
        Updates the held values

        :param val:
        :type val:
        :return:
        :rtype:
        """
        ...

    def load_dipole_surface(self):
        """
        **LLM Docstring**

        Not currently implemented: general (non-derivative-expansion) dipole surfaces aren't supported yet.

        :return: never returns
        :rtype: object
        :raises NotImplementedError: always, since general dipole surfaces aren't supported
        """
        ...

    def _load_gaussian_fchk_dipoles(self, file, masses=None, freqs=None):
        """
        **LLM Docstring**

        Parse dipole moment, gradient, and (if present) second/third derivatives plus their numerical counterparts out of a Gaussian FChk file, applying the necessary atomic-mass-unit unit conversions to each derivative order.

        :param file: path to the Gaussian FChk file
        :type file: str
        :param masses: unused in the current implementation
        :type masses: np.ndarray | None
        :param freqs: unused in the current implementation
        :type freqs: np.ndarray | None
        :return: a dict with `'analytic'` and `'numerical'` keys, each a tuple of the available derivative orders (moment, gradient, and higher orders where present)
        :rtype: dict
        """
        ...

    def load_dipole_derivatives(self, file=None, quiet=False):
        """
        Loads dipole derivatives from a file (or from `source_file` if set)

        :param file:
        :type file:
        :return:
        :rtype:
        """
        ...

    def load_polarizability_derivatives(self, file=None, quiet=False):
        """
        Loads dipole derivatives from a file (or from `source_file` if set)

        :param file:
        :type file:
        :return:
        :rtype:
        """
        ...

    def apply_transformation(self, transf):
        """
        **LLM Docstring**

        Apply a coordinate transformation to a copy of this dipole-surface manager, transforming the surface object (if set) and/or the derivative tensors (if set); non-linear transformations of the derivatives are not yet supported.

        :param transf: the transformation to apply, either an object exposing `transformation_function` or a `3x3` transformation matrix
        :type transf: object | np.ndarray
        :return: a new `DipoleSurfaceManager` with the transformation applied
        :rtype: DipoleSurfaceManager
        :raises NotImplementedError: if the derivatives are set and the transformation is not linear
        """
        ...

    def insert_atoms(self, atoms, coords, where):
        """
        Handles the insertion of new atoms into the structure

        :param atoms:
        :type atoms: tuple[str]
        :param coords:
        :type coords: CoordinateSet
        :param where:
        :type where: tuple[int]
        :return:
        :rtype:
        """
        ...

    def delete_atoms(self, where):
        """
        Handles the deletion from the structure

        :param atoms:
        :type atoms: tuple[str]
        :param coords:
        :type coords: CoordinateSet
        :param where:
        :type where: tuple[int]
        :return:
        :rtype:
        """
        ...

class PotentialSurfaceManager(PropertyManager):
    name = 'PotentialSurface'

    def __init__(self, mol, surface=None, derivatives=None):
        """
        **LLM Docstring**

        Set up the potential-surface manager, either copying the surface/derivatives/surface-coordinates from an existing surface-like object, or storing the supplied `surface`/`derivatives` directly.

        :param mol: the molecule this manager is attached to
        :type mol: AbstractMolecule
        :param surface: an existing potential surface, or another object exposing `_surf`/`_derivs`/`_surface_coords` to copy from
        :type surface: object | None
        :param derivatives: raw potential derivative tensors
        :type derivatives: list[np.ndarray] | None
        :return: None
        :rtype: None
        """
        ...

    @classmethod
    def from_data(cls, mol, data):
        """
        **LLM Docstring**

        Abstract constructor for building a potential-surface manager directly from raw data. Not implemented.

        :param mol: the molecule the manager will be attached to
        :type mol: AbstractMolecule
        :param data: the raw data to build from
        :type data: object
        :return: never returns
        :rtype: PotentialSurfaceManager
        :raises NotImplementedError: always
        """
        ...

    @property
    def surface(self):
        """
        **LLM Docstring**

        The potential surface object, lazily loaded via `load_potential_surface(self.surface_coords)` if not already set.

        :return: the potential surface
        :rtype: object
        """
        ...

    @property
    def surface_coords(self):
        """
        **LLM Docstring**

        Property getter/setter for the coordinate specification used when loading the potential surface (e.g. bond/angle/dihedral tuples).

        :param coords: (setter only) the new coordinate specification
        :type coords: object
        :return: (getter) the stored coordinate specification
        :rtype: object
        """
        ...

    @surface_coords.setter
    def surface_coords(self, coords):
        """
        **LLM Docstring**

        Property getter/setter for the coordinate specification used when loading the potential surface (e.g. bond/angle/dihedral tuples).

        :param coords: (setter only) the new coordinate specification
        :type coords: object
        :return: (getter) the stored coordinate specification
        :rtype: object
        """
        ...

    def get_derivs(self, quiet=False):
        """
        **LLM Docstring**

        Fetch the potential-energy derivative tensors, lazily loading them via `load_potential_derivatives` if not already cached.

        :param quiet: if `True`, suppresses errors from the underlying loader when the derivatives can't be found
        :type quiet: bool
        :return: the potential derivative tensors
        :rtype: tuple
        """
        ...

    @property
    def derivatives(self):
        """
        **LLM Docstring**

        Property getter/setter for the potential-energy derivative tensors. The getter delegates to `get_derivs()`.

        :param v: (setter only) the new derivative tensors to store
        :type v: tuple
        :return: (getter) the potential derivative tensors
        :rtype: tuple
        """
        ...

    @derivatives.setter
    def derivatives(self, v):
        """
        **LLM Docstring**

        Property getter/setter for the potential-energy derivative tensors. The getter delegates to `get_derivs()`.

        :param v: (setter only) the new derivative tensors to store
        :type v: tuple
        :return: (getter) the potential derivative tensors
        :rtype: tuple
        """
        ...

    @property
    def force_constants(self):
        """
        **LLM Docstring**

        The force-constant (second-derivative/Hessian) tensor from the potential derivative expansion.

        :return: `self.derivatives[1]`
        :rtype: np.ndarray
        """
        ...

    def load_potential_derivatives(self, file=None, mode=None, quiet=False):
        """
        Loads potential derivatives from a file (or from `source_file` if set)

        :param file:
        :type file:
        :return:
        :rtype:
        """
        ...

    def load_potential_surface(self, coordinates):
        """
        **LLM Docstring**

        Build a `PotentialSurface` object by loading a Gaussian scan log file (`self.mol.source_file`) and, if `coordinates` gives an iterable of atom-index tuples (bond/angle/dihedral specs), constructing the coordinate-transform function used to interpret the scan coordinates; otherwise `coordinates` itself is treated as that transform function.

        :param coordinates: either a callable that maps Cartesian coordinates to scan coordinates, or an iterable of 2/3/4-atom index tuples specifying bond lengths/angles/dihedrals to use as scan coordinates
        :type coordinates: callable | Iterable[tuple]
        :return: the loaded potential surface
        :rtype: PotentialSurface
        :raises ValueError: if an entry in `coordinates` isn't a 2, 3, or 4-atom tuple
        """
        ...

    def load(self, coordinates=None):
        """
        **LLM Docstring**

        Load the potential surface using `coordinates` if one is given and no surface is cached yet, then return whichever representation (surface or derivatives) is available.

        :param coordinates: coordinate specification forwarded to `load_potential_surface` if the surface still needs to be loaded
        :type coordinates: object | None
        :return: the potential surface or its derivatives
        :rtype: object
        """
        ...

    def update(self, val):
        """
        Updates the held values

        :param val:
        :type val:
        :return:
        :rtype:
        """
        ...

    def apply_transformation(self, transf):
        """
        **LLM Docstring**

        Apply a coordinate transformation to a copy of this potential-surface manager, transforming the surface object (if set) and/or the derivative tensors (if set) via either a direct linear transform or a general `TensorDerivativeConverter` re-expansion.

        :param transf: the transformation to apply, either an object exposing `transformation_function`, a `3x3` transformation matrix, or a general derivative-tensor sequence
        :type transf: object | np.ndarray | list
        :return: a new `PotentialSurfaceManager` with the transformation applied
        :rtype: PotentialSurfaceManager
        """
        ...

    def insert_atoms(self, atoms, coords, where):
        """
        Handles the insertion of new atoms into the structure

        :param atoms:
        :type atoms: tuple[str]
        :param coords:
        :type coords: CoordinateSet
        :param where:
        :type where: tuple[int]
        :return:
        :rtype:
        """
        ...

    def delete_atoms(self, where):
        """
        Handles the deletion from the structure

        :param atoms:
        :type atoms: tuple[str]
        :param coords:
        :type coords: CoordinateSet
        :param where:
        :type where: tuple[int]
        :return:
        :rtype:
        """
        ...

class NormalModesManager(PropertyManager):

    def __init__(self, mol, normal_modes=None):
        """
        **LLM Docstring**

        Set up the normal-modes manager, converting `normal_modes` into a `MolecularVibrations` object from a dict spec, rebinding it to `mol` if it already carries mode data, or leaving it as `None`.

        :param mol: the molecule this manager is attached to
        :type mol: AbstractMolecule
        :param normal_modes: the normal modes, given as a `MolecularVibrations`-like object, a dict with `'matrix'`/`'freqs'`/(optionally) `'inverse'`/`'origin'` keys, or `None`
        :type normal_modes: object | dict | None
        :return: None
        :rtype: None
        """
        ...

    @classmethod
    def from_data(cls, mol, data):
        """
        **LLM Docstring**

        Build a `NormalModesManager` from raw normal-mode data of various shapes: an existing `MolecularVibrations` (or `None`), a dict of matrix/frequency data, or generic mode data understood by `NormalModes.prep_modes`.

        :param mol: the molecule the manager will be attached to
        :type mol: AbstractMolecule
        :param data: the raw normal-mode data
        :type data: MolecularVibrations | dict | object | None
        :return: the constructed manager
        :rtype: NormalModesManager
        """
        ...

    def set_molecule(self, mol):
        """
        **LLM Docstring**

        Rebind this manager (and its stored modes, if any) to a different molecule.

        :param mol: the new molecule to associate
        :type mol: AbstractMolecule
        :return: None
        :rtype: None
        """
        ...

    def get_modes(self, quiet=False, allow_compute=True):
        """
        **LLM Docstring**

        Fetch the normal modes, computing them via `get_normal_modes` if not already cached.

        :param quiet: if `True`, suppresses errors when the modes can't be found/computed
        :type quiet: bool
        :param allow_compute: whether force constants may be computed from an energy evaluator if not otherwise available
        :type allow_compute: bool
        :return: the normal modes
        :rtype: MolecularVibrations
        """
        ...

    @property
    def modes(self):
        """

        :return:
        :rtype: MolecularVibrations
        """
        ...

    @modes.setter
    def modes(self, modes):
        """
        **LLM Docstring**

        Setter for the `modes` property: coerces `modes` into a `MolecularVibrations` object via `construct_normal_modes` if it isn't already one, then stores it.

        :param modes: the new normal modes, in any form accepted by `construct_normal_modes`
        :type modes: MolecularVibrations | dict | np.ndarray | object
        :return: None
        :rtype: None
        :raises TypeError: if the coerced value still isn't a `MolecularVibrations` instance
        """
        ...

    def construct_normal_modes(self, modes):
        """
        **LLM Docstring**

        Coerce a variety of mode representations (a dict of matrix/kwargs, a raw `np.ndarray`, an object exposing `modes_by_coords`/`coords_by_modes`, or an existing `MolecularNormalModes`) into a `MolecularVibrations` object bound to this manager's molecule.

        :param modes: the mode data to coerce
        :type modes: dict | np.ndarray | MolecularNormalModes | object
        :return: the constructed (or passed-through) `MolecularVibrations` object
        :rtype: MolecularVibrations
        :raises ValueError: if `modes` doesn't match any recognized representation
        """
        ...

    def load(self):
        """
        **LLM Docstring**

        Return the currently cached normal modes without attempting to compute them.

        :return: the cached modes, or `None` if none have been set
        :rtype: MolecularVibrations | None
        """
        ...

    def update(self, modes):
        """

        :return:
        :rtype:
        """
        ...

    @classmethod
    def _parse_orca_file_modes(cls, mol, file, quiet=False):
        """
        **LLM Docstring**

        Parse vibrational frequencies and normal modes from an ORCA log file, dropping the 6 translational/rotational modes, mass-weighting/renormalizing the mode vectors using the molecule's G-matrix, and converting frequencies to Hartrees.

        :param mol: the molecule the modes belong to (used for its G-matrix)
        :type mol: AbstractMolecule
        :param file: path to the ORCA log file
        :type file: str
        :param quiet: accepted for interface consistency with sibling parsers but not used in this method's body
        :type quiet: bool
        :return: `(freqs, modes, inv)` -- the vibrational frequencies, mode matrix, and its (pseudo)inverse
        :rtype: tuple[np.ndarray, np.ndarray, np.ndarray]
        """
        ...

    @classmethod
    def _parse_log_file_modes(cls, file, quiet=False):
        """
        **LLM Docstring**

        Parse vibrational frequencies and normal modes from a Gaussian `.log` file, using the parsed reduced masses and atomic masses to renormalize and mass-weight the mode vectors.

        :param file: path to the Gaussian log file
        :type file: str
        :param quiet: if `True`, returns `None` instead of raising when coordinates or normal modes are missing from the file
        :type quiet: bool
        :return: `(freqs, modes, inv)` -- the vibrational frequencies (in Hartrees), mode matrix, and its inverse; or `None` if `quiet` and the required data is missing
        :rtype: tuple[np.ndarray, np.ndarray, np.ndarray] | None
        :raises ValueError: if `quiet` is `False` and coordinates or normal modes are missing from the file
        """
        ...

    @classmethod
    def _parse_fchk_file_modes(cls, mol, file, rephase=True, quiet=False):
        """
        **LLM Docstring**

        Parse vibrational frequencies and normal modes from a Gaussian FChk file, mass-weighting/renormalizing them using the reduced masses and atomic weights, and (optionally) rephasing them to match the dipole-derivative-based convention.

        :param mol: the molecule the modes belong to
        :type mol: AbstractMolecule
        :param file: path to the Gaussian FChk file
        :type file: str
        :param rephase: whether to rephase the modes using `get_fchk_normal_mode_rephasing`
        :type rephase: bool
        :param quiet: accepted for interface consistency with sibling parsers but not used in this method's body
        :type quiet: bool
        :return: the parsed (and possibly rephased) normal modes
        :rtype: MolecularNormalModes
        """
        ...

    @classmethod
    def _parse_molpro_file_modes(cls, mol, file, quiet=False):
        """
        **LLM Docstring**

        Parse vibrational frequencies and normal modes from a MOLPRO log file, discarding near-zero (translational/rotational) frequencies, correcting any imaginary-frequency sign discontinuities, sorting by frequency, and mass-weighting the mode vectors using the molecule's G-matrix.

        :param mol: the molecule the modes belong to (used for its G-matrix)
        :type mol: AbstractMolecule
        :param file: path to the MOLPRO log file
        :type file: str
        :param quiet: if `True`, returns `None` instead of raising when the normal-modes block is missing
        :type quiet: bool
        :return: `(freqs, modes, inv)` -- the vibrational frequencies (in Hartrees), mode matrix, and its inverse
        :rtype: tuple[np.ndarray, np.ndarray, np.ndarray] | None
        :raises ValueError: if `quiet` is `False` and the normal-modes block is missing from the file
        """
        ...
    recalc_normal_mode_tolerance = 1e-08

    def load_normal_modes(self, file=None, mode=None, rephase=True, recalculate=False, quiet=False):
        """
        Loads potential derivatives from a file (or from `source_file` if set)

        :param file:
        :type file:
        :param rephase: whether to rephase FChk normal modes or not
        :type rephase: bool
        :return:
        :rtype:
        """
        ...

    def get_normal_modes(self, quiet=False, compute_force_constants=True, **kwargs):
        """
        Loads normal modes from file or calculates
        from force constants

        :param kwargs:
        :type kwargs:
        :return:
        :rtype:
        """
        ...

    def get_force_constants(self, compute_force_constants=True, quiet=False):
        """
        **LLM Docstring**

        Fetch the potential-energy Hessian (force constants), computing it from the molecule's energy evaluator via a numerical order-2 calculation if not already stored on the molecule.

        :param compute_force_constants: whether to allow computing the force constants via `self.mol.calculate_energy` if not already available
        :type compute_force_constants: bool
        :param quiet: if `True`, returns `None` instead of raising when force constants aren't available and can't be computed
        :type quiet: bool
        :return: the force-constant (Hessian) matrix
        :rtype: np.ndarray | None
        :raises ValueError: if `quiet` is `False`, no derivatives are stored, and no energy evaluator is available (or `compute_force_constants` is `False`)
        """
        ...

    @classmethod
    def get_dipole_derivative_based_rephasing(cls, modes, analytic_dipoles, numerical_dipoles, strict=True, allow_swaps=False):
        """
        **LLM Docstring**

        Derive the sign (and, in limited cases, swap) corrections needed to align a set of normal modes' numerically-computed dipole derivatives with their analytically-computed counterparts, by comparing the inner-product matrix between the two (each normalized) and flagging modes whose overlap magnitude falls below a threshold as needing special handling.

        :param modes: the normal modes whose phase convention is being checked/corrected
        :type modes: MolecularNormalModes
        :param analytic_dipoles: the analytic dipole derivative expansion (only the first-derivative term is used)
        :type analytic_dipoles: tuple | None
        :param numerical_dipoles: the numerical dipole derivative expansion (only the first-derivative term is used)
        :type numerical_dipoles: tuple | None
        :param strict: if `True`, raise on shape mismatches between the numerical and analytic derivatives rather than returning `None`
        :type strict: bool
        :param allow_swaps: if `True`, attempt to resolve ambiguous (nearly-degenerate) modes by finding candidate index swaps instead of only sign corrections (currently raises with the candidate mapping rather than resolving it automatically)
        :type allow_swaps: bool
        :return: the rephasing matrix to apply to the modes, or `None` if the necessary dipole data isn't available
        :rtype: np.ndarray | None
        :raises ValueError: if `strict` is `True` and the numerical/analytic derivative shapes don't align, or if more than one mode needs rephasing and `allow_swaps` is `False`
        """
        ...

    @classmethod
    def get_fchk_normal_mode_rephasing(cls, mol, modes, use_dipoles=False):
        """
        Returns the necessary rephasing to make the numerical dipole derivatives
        agree with the analytic dipole derivatives as pulled from a Gaussian FChk file
        :return:
        :rtype:
        """
        ...

    @classmethod
    def _handle_rephasing_update(cls, together, opposite, rephasing_targets, non_rephased, strict=True):
        """
        **LLM Docstring**

        Update the running record of which modes must be rephased together versus oppositely, merging the new `together`/`opposite` constraint sets into `rephasing_targets` for every mode involved (skipping modes already marked as non-rephased).

        :param together: mode indices that must share the same sign correction as each other
        :type together: list[int]
        :param opposite: mode indices that must have the opposite sign correction to `together`
        :type opposite: list[int]
        :param rephasing_targets: mapping from mode index to a `(same_sign_set, opposite_sign_set)` pair, updated in place
        :type rephasing_targets: dict
        :param non_rephased: set of mode indices already excluded from rephasing
        :type non_rephased: set
        :param strict: accepted for interface consistency but not used in this method's current body (the contradiction-checking logic that used it is commented out)
        :type strict: bool
        :return: None (updates `rephasing_targets` in place)
        :rtype: None
        """
        ...

    @classmethod
    def get_partial_cubic_based_rephasing(cls, modes, partial_cubics: np.ndarray, degenerate_freq_threshold=1e-06, equivalent_threshold=0.01, nonzero_threshold=4.5e-06, ignore_single_zeros=False, frequency_scale=True):
        """
        **LLM Docstring**

        Derive a rephasing (sign/permutation) matrix for a set of nearly-degenerate normal modes using consistency of the partial cubic force constants: transforms `partial_cubics` into the mode basis (optionally frequency-scaled), groups modes into degenerate blocks by frequency, resolves pairwise sign/swap ambiguities within each block by comparing triples of cubic force constants that should agree up to the rephasing, and returns the resulting rephasing transformation.

        :param modes: the normal modes to rephase
        :type modes: MolecularNormalModes
        :param partial_cubics: the partial cubic force-constant tensor used to disambiguate the phase/ordering of degenerate modes
        :type partial_cubics: np.ndarray
        :param degenerate_freq_threshold: maximum frequency difference for two modes to be treated as degenerate
        :type degenerate_freq_threshold: float
        :param equivalent_threshold: relative-difference threshold for treating two cubic force-constant magnitudes as equal
        :type equivalent_threshold: float
        :param nonzero_threshold: absolute magnitude below which a cubic force constant is treated as zero
        :type nonzero_threshold: float
        :param ignore_single_zeros: accepted parameter but not referenced in the current method body
        :type ignore_single_zeros: bool
        :param frequency_scale: whether to frequency-scale the cubic tensor before comparing terms
        :type frequency_scale: bool
        :return: the rephasing (sign/permutation) matrix to apply to the modes
        :rtype: np.ndarray
        :raises ValueError: if a degenerate block can't be resolved with a single pairwise swap, or if inconsistent (asymmetric) rephasing requirements are detected among a triple of modes
        """
        ...

    def apply_transformation(self, transf):
        """
        **LLM Docstring**

        Apply a coordinate transformation to a copy of this normal-modes manager by embedding its modes through the transformation.

        :param transf: the transformation to apply to the modes
        :type transf: object
        :return: a new `NormalModesManager` with the transformed modes
        :rtype: NormalModesManager
        """
        ...

    def insert_atoms(self, atoms, coords, where):
        """
        Handles the insertion of new atoms into the structure

        :param atoms:
        :type atoms: tuple[str]
        :param coords:
        :type coords: CoordinateSet
        :param where:
        :type where: tuple[int]
        :return:
        :rtype:
        """
        ...

    def delete_atoms(self, where):
        """
        Handles the deletion from the structure

        :param atoms:
        :type atoms: tuple[str]
        :param coords:
        :type coords: CoordinateSet
        :param where:
        :type where: tuple[int]
        :return:
        :rtype:
        """
        ...