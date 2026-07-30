"""
Provides a simple Molecule class that we can adapt as we need
Uses AtomData to get properties and whatnot
"""
from __future__ import annotations
import io
import os, numpy as np
import tempfile
import typing
import collections
from McUtils.Data import AtomData, UnitsData, BondData
from McUtils.Coordinerds import CoordinateSet, ZMatrixCoordinates, CartesianCoordinates3D, PrimitiveCoordinatePicker, RedundantCoordinateGenerator
import McUtils.Numputils as nput
import McUtils.Iterators as itut
import McUtils.Devutils as dev
import McUtils.Coordinerds as coordops
import McUtils.Zachary as zach
from McUtils.Zachary import Mesh
from McUtils.Graphs import EdgeGraph
import McUtils.Plots as plt
from McUtils.ExternalPrograms import RDMolecule, ExternalProgramJob
from McUtils.Scaffolding import Logger
import McUtils.Symmetry as symm
from .MoleculeInterface import *
from .CoordinateSystems import MolecularEmbedding, ModeEmbedding
from .Evaluator import MolecularEvaluator, EnergyEvaluator, DipoleEvaluator, ChargeEvaluator, ReducedDimensionalPotentialHandler, DipolePolarizabilityEvaluator
from .Hamiltonian import MolecularHamiltonian
from .Properties import *
from .Serializers import MoleculePropertyCache
__all__ = ['Molecule', 'MolecoolException']
__reload_hook__ = ['..Modes', '.MoleculeInterface', '.CoordinateSystems', '.Hamiltonian', '.Evaluator', '.Properties']
from .Transformations import MolecularTransformation

class Molecule(AbstractMolecule):
    """
    General purpose 'Molecule' class where the 'Molecule' need not be a molecule at all
    """

    def __init__(self, atoms, coords, bonds=None, masses=None, name=None, internals=None, rdmol=None, dipole_surface=None, dipole_derivatives=None, potential_surface=None, potential_derivatives=None, normal_modes=None, source_file=None, guess_bonds=True, charge=None, formal_charges=None, spin=None, display_mode=None, display_settings=None, energy=None, energy_evaluator=None, dipole_evaluator=None, charge_evaluator=None, polarizability_evaluator=None, polarizability_derivatives=None, checkpoint_file=None, **metadata):
        """
        :param atoms: atoms specified by name, either full name or short
        :type atoms: Iterable[str]
        :param coords: coordinates for the molecule, assumed to be in Bohr by default
        :type coords: np.ndarray | Iterable[Iterable[float]]
        :param bonds: bond specification for the molecule
        :type bonds: Iterable[Iterable[int]] | None
        :param obmol: OpenBabel molecule for doing conversions
        :type obmol:
        :param charge: Net charge on the molecule
        :type charge: int | None
        :param name: Name for the molecule
        :type name: str | None
        :param name: The internal coordinate specification for the molecule
        :type name: np.ndarray[int] | None
        :param dipole_surface: The dipole surface for the system
        :type dipole_surface: DipoleSurface | None
        :param dipole_derivatives: Derivatives of the dipole surface
        :type dipole_derivatives: Iterable[np.ndarray] | None
        :param potential_surface: The potential surface for the system
        :type potential_surface: PotentialSurface | None
        :param potential_derivatives: Derivatives of the potential surface
        :type potential_derivatives: Iterable[np.ndarray] | None
        :param guess_bonds: Whether or not to guess the bonding arrangement when that would be used
        :type guess_bonds: bool
        :param source_file: The data file the molecule was loaded from
        :type source_file: str
        :param kw: Other bound parameters that might be useful
        :type kw:
        """
        ...

    def modify(self, atoms=dev.default, coords=dev.default, *, internals=dev.default, masses=dev.default, bonds=dev.default, guess_bonds=dev.default, energy=dev.default, energy_evaluator=dev.default, dipole_evaluator=dev.default, charge_evaluator=dev.default, polarizability_evaluator=dev.default, charge=dev.default, spin=dev.default, rdmol=dev.default, display_mode=dev.default, display_settings=dev.default, normal_modes=dev.default, dipole_surface=dev.default, potential_surface=dev.default, dipole_derivatives=dev.default, potential_derivatives=dev.default, polarizability_derivatives=dev.default, meta=dev.default, source_file=dev.default):
        """
        **LLM Docstring**

        Build a new `Molecule` that is a copy of this one with the given fields overridden, treating any argument left at its `dev.default` sentinel as "keep the current value" (with some fields, like masses/energy/surfaces, only carried over automatically when the arguments they depend on -- e.g. `atoms`, `coords`, `energy_evaluator` -- are also left unspecified).

        :param atoms: replacement atoms, or `dev.default` to keep the current ones
        :type atoms: Iterable[str] | object
        :param coords: replacement coordinates, or `dev.default` to keep the current ones
        :type coords: np.ndarray | object
        :param internals: replacement internal-coordinate specification
        :type internals: object
        :param masses: replacement masses
        :type masses: np.ndarray | object
        :param bonds: replacement bonds
        :type bonds: object
        :param guess_bonds: replacement bond-guessing flag
        :type guess_bonds: bool | object
        :param energy: replacement cached energy value
        :type energy: float | object
        :param energy_evaluator: replacement energy evaluator
        :type energy_evaluator: object
        :param dipole_evaluator: replacement dipole evaluator
        :type dipole_evaluator: object
        :param charge_evaluator: replacement charge evaluator
        :type charge_evaluator: object
        :param polarizability_evaluator: replacement polarizability evaluator
        :type polarizability_evaluator: object
        :param charge: replacement net charge
        :type charge: int | object
        :param spin: replacement spin
        :type spin: object
        :param rdmol: replacement RDKit molecule
        :type rdmol: object
        :param display_mode: replacement display mode
        :type display_mode: str | object
        :param display_settings: replacement display settings
        :type display_settings: dict | object
        :param normal_modes: replacement normal modes
        :type normal_modes: object
        :param dipole_surface: replacement dipole surface
        :type dipole_surface: object
        :param potential_surface: replacement potential surface
        :type potential_surface: object
        :param dipole_derivatives: replacement dipole derivatives
        :type dipole_derivatives: object
        :param potential_derivatives: replacement potential derivatives
        :type potential_derivatives: object
        :param polarizability_derivatives: replacement polarizability derivatives
        :type polarizability_derivatives: object
        :param meta: replacement/merged metadata
        :type meta: dict | object
        :param source_file: replacement source file path
        :type source_file: str | object
        :return: the new, modified `Molecule`
        :rtype: Molecule
        """
        ...

    def __del__(self):
        """
        **LLM Docstring**

        Clean up the molecule's coordinate embedding (deregistering any converters it registered) when the object is garbage collected.

        :return: None
        :rtype: None
        """
        ...

    def to_state(self, serializer=None):
        """
        **LLM Docstring**

        Serialize this molecule's essential data (atoms, masses, coordinates, bonds, internal-coordinate spec, evaluators, potential/dipole derivatives, charge, spin) into a plain dict, stripping out non-serializable embedding-specific converter options first.

        :param serializer: accepted for interface consistency but not used in this method's body
        :type serializer: object | None
        :return: the serialized state dict
        :rtype: dict
        """
        ...

    @classmethod
    def from_state(cls, data, serializer=None):
        """
        **LLM Docstring**

        Reconstruct a `Molecule` from a previously serialized state dict, by passing its entries directly as constructor keyword arguments.

        :param data: the serialized state, as produced by `to_state`
        :type data: dict
        :param serializer: accepted for interface consistency but not used in this method's body
        :type serializer: object | None
        :return: the reconstructed molecule
        :rtype: Molecule
        """
        ...

    def cached_eval(self, key, generator, *, condition=None, args=(), kwargs=None):
        """
        **LLM Docstring**

        Evaluate (and cache, via this molecule's on-disk/in-memory checkpoint) a value under `key`, delegating to `self.checkpoint.cached_eval`.

        :param key: the cache key to look up or populate
        :type key: str
        :param generator: callable used to compute the value when it is not already cached
        :type generator: callable
        :param condition: optional predicate controlling whether the cached value should be recomputed
        :type condition: callable | None
        :param args: positional arguments passed to `generator` if it is called
        :type args: tuple
        :param kwargs: keyword arguments passed to `generator` if it is called
        :type kwargs: dict | None
        :return: the cached or newly computed value
        :rtype: object
        """
        ...

    @classmethod
    def _generate_auto_spec(cls, atoms, bonds, base_coords=None, **opts):
        """
        **LLM Docstring**

        Automatically pick a set of primitive internal coordinates (bond stretches, angles, dihedrals) from the bonding graph, via `PrimitiveCoordinatePicker`.

        :param atoms: the atom labels
        :type atoms: Iterable[str]
        :param bonds: the bonds to build coordinates from
        :type bonds: Iterable[Iterable[int]]
        :param base_coords: coordinates to seed/prioritize the picker with
        :type base_coords: Iterable | None
        :param opts: extra options forwarded to `PrimitiveCoordinatePicker`
        :type opts: dict
        :return: the picked coordinate specs
        :rtype: list
        """
        ...

    @classmethod
    def _generate_stretch_spec(cls, atoms, bonds, **opts):
        """
        **LLM Docstring**

        Build a "natural"-coordinate specification consisting only of the bond-stretch coordinates implied by `bonds`.

        :param atoms: the atom labels (unused directly, kept for interface consistency with `_generate_auto_spec`)
        :type atoms: Iterable[str]
        :param bonds: the bonds to build stretch coordinates from
        :type bonds: Iterable[Iterable[int]]
        :param opts: extra options, unused
        :type opts: dict
        :return: the list of bond-stretch coordinate specs
        :rtype: list
        """
        ...

    @classmethod
    def _auto_auto_spec(cls, spec_generator, atoms, coords, bonds, redundant=False, base_coordinates=None, masses=None, untransformed_coordinates=None, prune_coordinates=True, pruning_options=None, formal_charges=None, **opts):
        """
        **LLM Docstring**

        Shared driver behind `_auto_spec`/`_stretch_spec`: guesses bonds if not given, optionally sets up a redundant-coordinate specification (folding in any `untransformed_coordinates`), generates the primitive coordinate specs via `spec_generator`, and (optionally) prunes them down to a well-conditioned, non-redundant subset via `RedundantCoordinateGenerator.prune_coordinate_specs`.

        :param spec_generator: the coordinate-generating function to use (`_generate_auto_spec` or `_generate_stretch_spec`)
        :type spec_generator: callable
        :param atoms: the atom labels
        :type atoms: Iterable[str]
        :param coords: the Cartesian coordinates
        :type coords: np.ndarray
        :param bonds: the bonds to use; guessed via RDKit if `None`
        :type bonds: Iterable[Iterable[int]] | None
        :param redundant: whether to build a redundant coordinate specification
        :type redundant: bool
        :param base_coordinates: seed coordinates to prioritize/include
        :type base_coordinates: Iterable | None
        :param masses: atomic masses, used for pruning; computed from `atoms` if not given
        :type masses: np.ndarray | None
        :param untransformed_coordinates: coordinates that should remain untransformed under the redundant transformation
        :type untransformed_coordinates: Iterable | None
        :param prune_coordinates: whether to prune the generated coordinate specs down to a well-conditioned subset
        :type prune_coordinates: bool
        :param pruning_options: extra options forwarded to `RedundantCoordinateGenerator.prune_coordinate_specs`
        :type pruning_options: dict | None
        :param formal_charges: formal charges used when guessing bonds
        :type formal_charges: Iterable | None
        :param opts: extra options forwarded to `spec_generator`
        :type opts: dict
        :return: the resulting internal-coordinate specification dict (with `'specs'`, and `'redundant'`/`'untransformed_coordinates'` if applicable)
        :rtype: dict
        """
        ...

    @classmethod
    def _auto_spec(cls, atoms, coords, bonds, **opts):
        """
        **LLM Docstring**

        Build an automatically-chosen internal-coordinate specification from the bonding graph, via `_auto_auto_spec` with `_generate_auto_spec`.

        :param atoms: the atom labels
        :type atoms: Iterable[str]
        :param coords: the Cartesian coordinates
        :type coords: np.ndarray
        :param bonds: the bonds to use
        :type bonds: Iterable[Iterable[int]] | None
        :param opts: extra options forwarded to `_auto_auto_spec`
        :type opts: dict
        :return: the resulting internal-coordinate specification dict
        :rtype: dict
        """
        ...

    @classmethod
    def _stretch_spec(cls, atoms, coords, bonds, **opts):
        """
        **LLM Docstring**

        Build a "natural"/stretch-only internal-coordinate specification from the bonding graph, via `_auto_auto_spec` with `_generate_stretch_spec`.

        :param atoms: the atom labels
        :type atoms: Iterable[str]
        :param coords: the Cartesian coordinates
        :type coords: np.ndarray
        :param bonds: the bonds to use
        :type bonds: Iterable[Iterable[int]] | None
        :param opts: extra options forwarded to `_auto_auto_spec`
        :type opts: dict
        :return: the resulting internal-coordinate specification dict
        :rtype: dict
        """
        ...

    def canonicalize_internals(self, spec, atoms, coords, bonds, relocalize=True, masses=None):
        """
        **LLM Docstring**

        Normalize the many accepted forms of an internal-coordinate specification (the strings `'auto'`/`'zmatrix'`, a dict with `'primitives'`/`'specs'`/`'zmatrix'` keys where `'specs'` may itself be `'auto'`/`'natural'`, a bare Z-matrix-like array, or a bare list of primitive specs) down into the canonical dict form expected by `MolecularEmbedding`, recursively re-dispatching as needed.

        :param spec: the internal-coordinate specification to canonicalize
        :type spec: str | dict | Iterable | None
        :param atoms: the atom labels
        :type atoms: Iterable[str]
        :param coords: the Cartesian coordinates
        :type coords: np.ndarray
        :param bonds: the bonds to use when auto-generating coordinates
        :type bonds: Iterable[Iterable[int]] | None
        :param relocalize: whether redundant coordinates should be relocalized by default
        :type relocalize: bool
        :param masses: atomic masses, forwarded to the auto-generation routines
        :type masses: np.ndarray | None
        :return: the canonicalized specification
        :rtype: dict | None
        :raises ValueError: if `spec` is a string that isn't recognized (`'auto'`/`'zmatrix'`/`'natural'`)
        """
        ...

    def prep_internal_spec(self, spec, relocalize=True, masses=None):
        """
        **LLM Docstring**

        Canonicalize an internal-coordinate specification against this molecule's own atoms, coordinates, bonds, and masses, via `canonicalize_internals`.

        :param spec: the internal-coordinate specification to canonicalize
        :type spec: object
        :param relocalize: whether redundant coordinates should be relocalized by default
        :type relocalize: bool
        :param masses: atomic masses to use instead of this molecule's own
        :type masses: np.ndarray | None
        :return: the canonicalized specification
        :rtype: dict | None
        """
        ...

    @property
    def embedding(self):
        """
        **LLM Docstring**

        Property getter/setter for the molecule's `MolecularEmbedding`. Setting it also resets the cached evaluator and Hamiltonian, since both depend on the embedding.

        :param e: (setter only) the new embedding
        :type e: MolecularEmbedding
        :return: (getter) the current embedding
        :rtype: MolecularEmbedding
        """
        ...

    @embedding.setter
    def embedding(self, e):
        """
        **LLM Docstring**

        Property getter/setter for the molecule's `MolecularEmbedding`. Setting it also resets the cached evaluator and Hamiltonian, since both depend on the embedding.

        :param e: (setter only) the new embedding
        :type e: MolecularEmbedding
        :return: (getter) the current embedding
        :rtype: MolecularEmbedding
        """
        ...

    def get_evaluator(self, embedding=None, normal_modes=dev.default):
        """
        **LLM Docstring**

        Build a `MolecularEvaluator` for this molecule's embedding (or an alternate one), using either the given `normal_modes` or this molecule's own.

        :param embedding: an alternate embedding to build the evaluator for; defaults to `self.embedding`
        :type embedding: MolecularEmbedding | None
        :param normal_modes: normal modes to use instead of `self._normal_modes`
        :type normal_modes: object
        :return: the constructed evaluator
        :rtype: MolecularEvaluator
        """
        ...

    @property
    def evaluator(self):
        """
        **LLM Docstring**

        Property getter/setter for the `MolecularEvaluator` used for energy/derivative calculations. The getter lazily builds one via `get_evaluator` the first time it's needed.

        :param e: (setter only) the new evaluator
        :type e: MolecularEvaluator
        :return: (getter) the cached (or newly built) evaluator
        :rtype: MolecularEvaluator
        """
        ...

    @evaluator.setter
    def evaluator(self, e):
        """
        **LLM Docstring**

        Property getter/setter for the `MolecularEvaluator` used for energy/derivative calculations. The getter lazily builds one via `get_evaluator` the first time it's needed.

        :param e: (setter only) the new evaluator
        :type e: MolecularEvaluator
        :return: (getter) the cached (or newly built) evaluator
        :rtype: MolecularEvaluator
        """
        ...

    @property
    def coords(self):
        """
        **LLM Docstring**

        Property getter/setter for the Cartesian coordinates, delegating to `self.embedding.coords`.

        :param coords: (setter only) the new Cartesian coordinates
        :type coords: np.ndarray
        :return: (getter) the Cartesian coordinates
        :rtype: CoordinateSet
        """
        ...

    @coords.setter
    def coords(self, coords):
        """
        **LLM Docstring**

        Property getter/setter for the Cartesian coordinates, delegating to `self.embedding.coords`.

        :param coords: (setter only) the new Cartesian coordinates
        :type coords: np.ndarray
        :return: (getter) the Cartesian coordinates
        :rtype: CoordinateSet
        """
        ...

    @property
    def masses(self):
        """
        **LLM Docstring**

        Property getter/setter for the atomic masses. The setter also updates the embedding's masses (via `atomic_masses`, i.e. in atomic units).

        :param masses: (setter only) the new masses
        :type masses: np.ndarray
        :return: (getter) the atomic masses
        :rtype: np.ndarray
        """
        ...

    @masses.setter
    def masses(self, masses):
        """
        **LLM Docstring**

        Property getter/setter for the atomic masses. The setter also updates the embedding's masses (via `atomic_masses`, i.e. in atomic units).

        :param masses: (setter only) the new masses
        :type masses: np.ndarray
        :return: (getter) the atomic masses
        :rtype: np.ndarray
        """
        ...

    @property
    def internals(self):
        """
        **LLM Docstring**

        Getter for the raw (canonicalized) internal-coordinate specification, delegating to `self.embedding.internals`. (A setter with the same name separately rebuilds the embedding from a new specification.)

        :return: the internal-coordinate specification, or `None` if none is set
        :rtype: dict | None
        """
        ...

    @property
    def charge(self):
        """
        **LLM Docstring**

        Property getter/setter for the molecule's net charge, stored in its metadata (getter defaults to `0` if unset).

        :param c: (setter only) the new net charge
        :type c: int
        :return: (getter) the net charge
        :rtype: int
        """
        ...

    @charge.setter
    def charge(self, c):
        """
        **LLM Docstring**

        Property getter/setter for the molecule's net charge, stored in its metadata (getter defaults to `0` if unset).

        :param c: (setter only) the new net charge
        :type c: int
        :return: (getter) the net charge
        :rtype: int
        """
        ...

    @property
    def spin(self):
        """
        **LLM Docstring**

        Property getter/setter for the molecule's spin, stored in its metadata.

        :param c: (setter only) the new spin value
        :type c: object
        :return: (getter) the spin, or `None` if unset
        :rtype: object | None
        """
        ...

    @spin.setter
    def spin(self, c):
        """
        **LLM Docstring**

        Property getter/setter for the molecule's spin, stored in its metadata.

        :param c: (setter only) the new spin value
        :type c: object
        :return: (getter) the spin, or `None` if unset
        :rtype: object | None
        """
        ...

    @property
    def charges(self):
        """
        **LLM Docstring**

        Property getter/setter for the molecule's per-atom partial charges, stored in its metadata.

        :param c: (setter only) the new per-atom charges
        :type c: np.ndarray
        :return: (getter) the per-atom charges, or `None` if unset
        :rtype: np.ndarray | None
        """
        ...

    @charges.setter
    def charges(self, c):
        """
        **LLM Docstring**

        Property getter/setter for the molecule's per-atom partial charges, stored in its metadata.

        :param c: (setter only) the new per-atom charges
        :type c: np.ndarray
        :return: (getter) the per-atom charges, or `None` if unset
        :rtype: np.ndarray | None
        """
        ...

    @property
    def formal_charges(self):
        """
        **LLM Docstring**

        Property getter/setter for the molecule's per-atom formal charges, stored in its metadata.

        :param c: (setter only) the new per-atom formal charges
        :type c: np.ndarray
        :return: (getter) the per-atom formal charges, or `None` if unset
        :rtype: np.ndarray | None
        """
        ...

    @formal_charges.setter
    def formal_charges(self, c):
        """
        **LLM Docstring**

        Property getter/setter for the molecule's per-atom formal charges, stored in its metadata.

        :param c: (setter only) the new per-atom formal charges
        :type c: np.ndarray
        :return: (getter) the per-atom formal charges, or `None` if unset
        :rtype: np.ndarray | None
        """
        ...

    def get_charge_evaluator(self, evaluator=None, **opts):
        """
        **LLM Docstring**

        Resolve (and, if needed, instantiate from this molecule) a charge-evaluator object, defaulting to `self.charge_evaluator` if none is given explicitly.

        :param evaluator: an explicit evaluator (or evaluator-type specification) to resolve; defaults to `self.charge_evaluator`
        :type evaluator: object | None
        :param opts: extra options forwarded to the evaluator's `from_mol` constructor, if applicable
        :type opts: dict
        :return: the resolved charge-evaluator instance
        :rtype: object
        :raises ValueError: if the evaluator type can't be resolved
        """
        ...

    def calculate_charges(self, evaluator=None, order=None, **opts):
        """
        **LLM Docstring**

        Compute partial-charge values (and, optionally, their derivatives) using the resolved charge evaluator at this molecule's current geometry.

        :param evaluator: an explicit evaluator (or evaluator-type specification) to use; defaults to `self.charge_evaluator`
        :type evaluator: object | None
        :param order: the highest derivative order to compute; if `None`, only the charges themselves are returned
        :type order: int | None
        :param opts: extra options forwarded to `get_charge_evaluator`
        :type opts: dict
        :return: the charges (if `order` is `None`) or the full charge/derivative expansion
        :rtype: np.ndarray | list[np.ndarray]
        """
        ...

    @internals.setter
    def internals(self, spec):
        """
        **LLM Docstring**

        Getter for the raw (canonicalized) internal-coordinate specification, delegating to `self.embedding.internals`. (A setter with the same name separately rebuilds the embedding from a new specification.)

        :return: the internal-coordinate specification, or `None` if none is set
        :rtype: dict | None
        """
        ...

    @property
    def internal_coordinates(self):
        """
        **LLM Docstring**

        Getter for the internal coordinates at the molecule's current geometry, delegating to `self.embedding.internal_coordinates`.

        :return: the internal coordinates, or `None` if none are defined
        :rtype: CoordinateSet | None
        """
        ...

    @property
    def redundant_internal_transformation(self):
        """
        **LLM Docstring**

        Getter for the redundant-to-non-redundant internal-coordinate transformation, delegating to `self.embedding.redundant_internal_transformation`.

        :return: the redundant transformation, or `None` if not applicable
        :rtype: np.ndarray | None
        """
        ...

    @classmethod
    def _check_label(cls, label, allowed_coordinate_types=None, excluded_coordinate_types=None, allowed_ring_types=None, excluded_ring_types=None, allowed_group_types=None, excluded_group_types=None):
        """
        **LLM Docstring**

        Test whether a coordinate label passes a set of allow/exclude filters on its atom types, ring membership, and functional-group membership.

        :param label: the coordinate label to test (exposing `.atoms`, `.ring`, `.group` attributes)
        :type label: object
        :param allowed_coordinate_types: if given, `label.atoms` must be among these to pass
        :type allowed_coordinate_types: Iterable | None
        :param excluded_coordinate_types: if given, `label.atoms` must not be among these to pass
        :type excluded_coordinate_types: Iterable | None
        :param allowed_ring_types: if given, `label.ring` must be among these to pass
        :type allowed_ring_types: Iterable | None
        :param excluded_ring_types: if given, `label.ring` must not be among these to pass
        :type excluded_ring_types: Iterable | None
        :param allowed_group_types: if given, `label.group` must be among these to pass
        :type allowed_group_types: Iterable | None
        :param excluded_group_types: if given, `label.group` must not be among these to pass
        :type excluded_group_types: Iterable | None
        :return: whether the label passes every specified filter
        :rtype: bool
        """
        ...

    @classmethod
    def get_coordinate_filer(cls, allowed_coordinate_types=None, excluded_coordinate_types=None, allowed_ring_types=None, excluded_ring_types=None, allowed_group_types=None, excluded_group_types=None):
        """
        **LLM Docstring**

        Build a filter function (closing over the given allow/exclude criteria) that, given a dict of coordinate-to-label mappings, returns only the entries whose label passes `_check_label`.

        :param allowed_coordinate_types: forwarded to `_check_label`
        :type allowed_coordinate_types: Iterable | None
        :param excluded_coordinate_types: forwarded to `_check_label`
        :type excluded_coordinate_types: Iterable | None
        :param allowed_ring_types: forwarded to `_check_label`
        :type allowed_ring_types: Iterable | None
        :param excluded_ring_types: forwarded to `_check_label`
        :type excluded_ring_types: Iterable | None
        :param allowed_group_types: forwarded to `_check_label`
        :type allowed_group_types: Iterable | None
        :param excluded_group_types: forwarded to `_check_label`
        :type excluded_group_types: Iterable | None
        :return: the constructed coordinate-filtering function
        :rtype: callable
        """
        ...
    default_coordinate_pruning = 'graph'

    def get_bond_graph_internals(self, include_stretches=True, include_bends=True, include_dihedrals=True, include_fragments=True, pruning=None, fragment=None, base_internals=None, use_distance_matrix=True, concatenate=True):
        """
        **LLM Docstring**

        Build a set of internal coordinates (bond stretches, bends, dihedrals, and/or inter-fragment coordinates) directly from the bonding graph, optionally restricted to a single fragment (recursively, with the result permuted back into the full atom indexing) and/or pruned down to a well-conditioned subset.

        :param include_stretches: whether to include bond-stretch coordinates
        :type include_stretches: bool
        :param include_bends: whether to include bond-angle coordinates
        :type include_bends: bool
        :param include_dihedrals: whether to include dihedral-angle coordinates
        :type include_dihedrals: bool
        :param include_fragments: whether to include coordinates connecting separate molecular fragments
        :type include_fragments: bool
        :param pruning: whether/how to prune the resulting coordinates (`True` for the default method, or an explicit method spec), forwarded to `prune_internals`
        :type pruning: bool | str | dict | None
        :param fragment: restrict to a single fragment, given as a fragment index or an explicit list of atom indices
        :type fragment: int | Iterable[int] | None
        :param base_internals: accepted and forwarded when recursing on a fragment, but not otherwise used directly in this method's own body
        :type base_internals: object | None
        :param use_distance_matrix: whether to precompute a distance matrix for the fragment-coordinate generation
        :type use_distance_matrix: bool
        :param concatenate: whether to concatenate the different coordinate categories (stretches/bends/dihedrals/fragments) into a single list, or return them as separate groups
        :type concatenate: bool
        :return: the generated internal coordinates, as a single concatenated list or a list of category groups depending on `concatenate`
        :rtype: list
        :raises ValueError: if `pruning` is requested while `concatenate` is `False`
        """
        ...

    def prune_internals(self, coords, method='b_matrix', check_rigidity=True):
        """
        **LLM Docstring**

        Reduce a set of internal coordinates down to a non-redundant, well-conditioned subset, defaulting to a B-matrix-rank-based method (building the necessary translation/rotation-projected B-matrix generator and a sensible `max_coords` cap) if no custom method is supplied.

        :param coords: the internal-coordinate specs to prune
        :type coords: list
        :param method: the pruning method: a method-name string, or a dict of method options (with a `'method'` key defaulting to `'b_matrix'`)
        :type method: str | dict
        :param check_rigidity: whether to check that the pruned coordinate set spans a rigid (non-redundant) representation
        :type check_rigidity: bool
        :return: the pruned coordinate specs
        :rtype: list
        """
        ...

    def get_labeled_internals(self, coordinate_filter=None, allowed_coordinate_types=None, excluded_coordinate_types=None, allowed_ring_types=None, excluded_ring_types=None, allowed_group_types=None, excluded_group_types=None, include_stretches=True, include_bends=True, include_dihedrals=True, include_fragments=True, coordinate_sorting=None, pruning=False):
        """
        **LLM Docstring**

        Build the internal coordinates from the bonding graph (via `get_bond_graph_internals`) and label each one by its atom types/ring/functional-group membership (via `edge_graph.get_label_types` and `coordops.get_coordinate_label`), then filter and sort them.

        :param coordinate_filter: an explicit filter function to apply instead of building one from the allow/exclude arguments
        :type coordinate_filter: callable | None
        :param allowed_coordinate_types: forwarded to `get_coordinate_filer` if `coordinate_filter` is not given
        :type allowed_coordinate_types: Iterable | None
        :param excluded_coordinate_types: forwarded to `get_coordinate_filer`
        :type excluded_coordinate_types: Iterable | None
        :param allowed_ring_types: forwarded to `get_coordinate_filer`
        :type allowed_ring_types: Iterable | None
        :param excluded_ring_types: forwarded to `get_coordinate_filer`
        :type excluded_ring_types: Iterable | None
        :param allowed_group_types: forwarded to `get_coordinate_filer`
        :type allowed_group_types: Iterable | None
        :param excluded_group_types: forwarded to `get_coordinate_filer`
        :type excluded_group_types: Iterable | None
        :param include_stretches: whether to include bond-stretch coordinates
        :type include_stretches: bool
        :param include_bends: whether to include bond-angle coordinates
        :type include_bends: bool
        :param include_dihedrals: whether to include dihedral-angle coordinates
        :type include_dihedrals: bool
        :param include_fragments: whether to include inter-fragment coordinates
        :type include_fragments: bool
        :param coordinate_sorting: a custom sorting function to apply to the labeled coordinates instead of the default `coordops.sort_internal_coordinates`; pass a falsy value to skip sorting
        :type coordinate_sorting: callable | bool | None
        :param pruning: whether/how to prune the coordinates, forwarded to `get_bond_graph_internals`
        :type pruning: bool | str | dict
        :return: a mapping from coordinate spec to its label, filtered and sorted
        :rtype: dict
        """
        ...

    def get_mode_labels(self, internals=None, modes=None, use_redundants=True, expansions=None, return_modes=False, **internals_opts):
        """
        **LLM Docstring**

        Assign human-readable labels (e.g. "C-H stretch") to a set of normal modes by projecting them onto labeled internal coordinates, handling both redundant and non-redundant internal-coordinate expansions and both Cartesian- and internal-coordinate-basis modes.

        :param internals: the labeled internal coordinates to project onto; computed via `get_labeled_internals` if not given
        :type internals: dict | None
        :param modes: the normal modes to label; computed via `get_normal_modes` if not given
        :type modes: object | None
        :param use_redundants: whether to build a redundant-coordinate expansion (with relocalization) for the projection, rather than using the internal coordinates directly
        :type use_redundants: bool
        :param expansions: precomputed `(expansions, inverse_expansion)` internal-coordinate Jacobian data to reuse instead of recomputing it
        :type expansions: tuple | None
        :param return_modes: whether to also return the internal-coordinate-basis mode matrix alongside the labels
        :type return_modes: bool
        :param internals_opts: extra options forwarded to `get_labeled_internals` if `internals` is not given
        :type internals_opts: dict
        :return: the mode labels, or `(internal_modes, labels)` if `return_modes` is set
        :rtype: list | tuple
        """
        ...

    @property
    def mode_embedding(self):
        """
        **LLM Docstring**

        The (cached) `ModeEmbedding` combining this molecule's coordinate embedding with its normal modes, built lazily the first time it's needed.

        :return: the mode embedding
        :rtype: ModeEmbedding
        """
        ...

    def get_internals(self, coords=None, *, strip_embedding=True):
        """
        **LLM Docstring**

        Fetch internal coordinates, either the molecule's own cached ones or those for an alternate set of Cartesian `coords`, via `self.embedding.get_internals`.

        :param coords: alternate Cartesian coordinates to convert instead of using the cached internal coordinates
        :type coords: np.ndarray | None
        :param strip_embedding: whether to strip the fixed embedding coordinates from the result
        :type strip_embedding: bool
        :return: the internal coordinates, or `None` if none are defined
        :rtype: CoordinateSet | None
        """
        ...

    def get_cartesians_by_internals(self, order=None, coords=None, *, strip_embedding=False, **kw):
        """
        **LLM Docstring**

        Fetch the Cartesians-by-internals Jacobian expansion, via `self.embedding.get_cartesians_by_internals`.

        :param order: the highest derivative order to compute
        :type order: int | None
        :param coords: alternate coordinates to compute the Jacobian at
        :type coords: np.ndarray | None
        :param strip_embedding: whether to strip the fixed embedding coordinates from the result
        :type strip_embedding: bool
        :param kw: extra options forwarded to the embedding
        :type kw: dict
        :return: the Cartesians-by-internals Jacobian tensors
        :rtype: list[np.ndarray]
        """
        ...

    def get_internals_by_cartesians(self, order=None, *, coords=None, strip_embedding=False, **kw):
        """
        **LLM Docstring**

        Fetch the internals-by-Cartesians Jacobian expansion, via `self.embedding.get_internals_by_cartesians`.

        :param order: the highest derivative order to compute
        :type order: int | None
        :param coords: alternate coordinates to compute the Jacobian at
        :type coords: np.ndarray | None
        :param strip_embedding: whether to strip the fixed embedding coordinates from the result
        :type strip_embedding: bool
        :param kw: extra options forwarded to the embedding
        :type kw: dict
        :return: the internals-by-Cartesians Jacobian tensors
        :rtype: list[np.ndarray]
        """
        ...

    def get_cartesians_by_modes(self, order=None, **kw):
        """
        **LLM Docstring**

        Fetch the Cartesians-by-normal-modes Jacobian expansion, via `self.mode_embedding.get_cartesians_by_internals`.

        :param order: the highest derivative order to compute
        :type order: int | None
        :param kw: extra options forwarded to the mode embedding
        :type kw: dict
        :return: the Cartesians-by-modes Jacobian tensors
        :rtype: list[np.ndarray]
        """
        ...

    def get_modes_by_cartesians(self, order=None, **kw):
        """
        **LLM Docstring**

        Fetch the normal-modes-by-Cartesians Jacobian expansion, via `self.mode_embedding.get_internals_by_cartesians`.

        :param order: the highest derivative order to compute
        :type order: int | None
        :param kw: extra options forwarded to the mode embedding
        :type kw: dict
        :return: the modes-by-Cartesians Jacobian tensors
        :rtype: list[np.ndarray]
        """
        ...

    @property
    def dipole_surface(self):
        """
        :return:
        :rtype: DipoleSurfaceManager
        """
        ...

    @dipole_surface.setter
    def dipole_surface(self, val):
        """
        **LLM Docstring**

        Setter for the `DipoleSurfaceManager`, requiring the new value to already be one. (The getter simply returns `self._dips`.)

        :param val: the new dipole-surface manager
        :type val: DipoleSurfaceManager
        :return: None
        :rtype: None
        :raises TypeError: if `val` isn't a `DipoleSurfaceManager`
        """
        ...

    @property
    def dipole_derivatives(self):
        """
        **LLM Docstring**

        Property getter/setter for the dipole derivative tensors. The getter delegates to `self.dipole_surface.get_derivatives(quiet=True)`; the setter assigns to `self.dipole_surface.derivatives`.

        :param derivs: (setter only) the new dipole derivative tensors
        :type derivs: list[np.ndarray]
        :return: (getter) the dipole derivative tensors, or `None` if unavailable
        :rtype: list[np.ndarray] | None
        """
        ...

    @dipole_derivatives.setter
    def dipole_derivatives(self, derivs):
        """
        **LLM Docstring**

        Property getter/setter for the dipole derivative tensors. The getter delegates to `self.dipole_surface.get_derivatives(quiet=True)`; the setter assigns to `self.dipole_surface.derivatives`.

        :param derivs: (setter only) the new dipole derivative tensors
        :type derivs: list[np.ndarray]
        :return: (getter) the dipole derivative tensors, or `None` if unavailable
        :rtype: list[np.ndarray] | None
        """
        ...

    def get_cartesian_dipole_derivatives(self, order=None, evaluator=None, include_constant_term=False):
        """
        **LLM Docstring**

        Fetch the dipole derivatives in Cartesian coordinates, computing them via `calculate_dipole` (and caching the result on the molecule, if the same evaluator is configured as the default) if not already available to the requested order.

        :param order: the highest derivative order needed; if `None`, whatever is available is returned
        :type order: int | None
        :param evaluator: an explicit dipole evaluator to use instead of `self.dipole_evaluator`
        :type evaluator: object | None
        :param include_constant_term: whether to include the zeroth-order (reference dipole) term in the result
        :type include_constant_term: bool
        :return: the dipole derivative tensors (from first order, or zeroth if `include_constant_term`), or `None` if unavailable
        :rtype: list[np.ndarray] | None
        """
        ...

    def get_internal_dipole_derivatives(self, order=None, reembed=True, strip_embedding=True):
        """
        **LLM Docstring**

        Fetch the dipole derivatives re-expressed in internal coordinates, by re-expanding the Cartesian dipole derivatives through the Cartesians-by-internals Jacobian.

        :param order: the highest derivative order needed
        :type order: int | None
        :param reembed: whether to use the Eckart-reembedded Cartesians-by-internals Jacobian
        :type reembed: bool
        :param strip_embedding: whether to strip the fixed embedding coordinates from the Jacobian
        :type strip_embedding: bool
        :return: the internal-coordinate dipole derivative tensors
        :rtype: list[np.ndarray]
        """
        ...

    def get_cartesian_polarizability_derivatives(self, order=None, evaluator=None, include_constant_term=False):
        """
        **LLM Docstring**

        Fetch the dipole-polarizability derivatives in Cartesian coordinates, computing them via `calculate_dipole_polarizability` (and caching the result, if the same evaluator is configured as the default) if not already available to the requested order, or re-expanding them through the normal modes if they were stored in a mode basis smaller than the full Cartesian space.

        :param order: the highest derivative order needed; if `None`, whatever is available is returned
        :type order: int | None
        :param evaluator: an explicit polarizability evaluator to use instead of `self.polarizability_evaluator`
        :type evaluator: object | None
        :param include_constant_term: accepted for interface consistency with `get_cartesian_dipole_derivatives` but not used in this method's body
        :type include_constant_term: bool
        :return: the polarizability derivative tensors
        :rtype: list[np.ndarray]
        """
        ...

    def get_hamiltonian(self, embedding=None, potential_derivatives=None, modes=None, dipole_derivatives=None, **etc):
        """
        **LLM Docstring**

        Build a `MolecularHamiltonian` for this molecule, defaulting to its own embedding, potential derivatives, normal modes, and dipole derivatives wherever not explicitly overridden.

        :param embedding: an alternate coordinate embedding to use
        :type embedding: MolecularEmbedding | None
        :param potential_derivatives: alternate potential-energy derivative tensors to use
        :type potential_derivatives: list[np.ndarray] | None
        :param modes: alternate normal modes to use
        :type modes: object | None
        :param dipole_derivatives: alternate dipole derivative tensors to use
        :type dipole_derivatives: list[np.ndarray] | None
        :param etc: extra options forwarded to the `MolecularHamiltonian` constructor
        :type etc: dict
        :return: the constructed Hamiltonian
        :rtype: MolecularHamiltonian
        """
        ...

    @property
    def hamiltonian(self):
        """
        **LLM Docstring**

        Property getter/setter for the `MolecularHamiltonian`. The getter lazily builds one via `get_hamiltonian` the first time it's needed.

        :param e: (setter only) the new Hamiltonian
        :type e: MolecularHamiltonian
        :return: (getter) the cached (or newly built) Hamiltonian
        :rtype: MolecularHamiltonian
        """
        ...

    @hamiltonian.setter
    def hamiltonian(self, e):
        """
        **LLM Docstring**

        Property getter/setter for the `MolecularHamiltonian`. The getter lazily builds one via `get_hamiltonian` the first time it's needed.

        :param e: (setter only) the new Hamiltonian
        :type e: MolecularHamiltonian
        :return: (getter) the cached (or newly built) Hamiltonian
        :rtype: MolecularHamiltonian
        """
        ...

    @property
    def potential_surface(self):
        """
        :return:
        :rtype: PotentialSurfaceManager
        """
        ...

    @potential_surface.setter
    def potential_surface(self, val):
        """
        **LLM Docstring**

        Setter for the `PotentialSurfaceManager`, requiring the new value to already be one. (The getter simply returns `self._pes`.)

        :param val: the new potential-surface manager
        :type val: PotentialSurfaceManager
        :return: None
        :rtype: None
        :raises TypeError: if `val` isn't a `PotentialSurfaceManager`
        """
        ...

    @property
    def potential_derivatives(self):
        """
        **LLM Docstring**

        Property getter/setter for the potential-energy derivative tensors. The getter fetches them from `self.potential_surface`, normalizing missing/placeholder entries to `0` and trimming any trailing zero-padded (unset) higher-order terms; the setter assigns to `self.potential_surface.derivatives`.

        :param derivs: (setter only) the new potential derivative tensors
        :type derivs: list[np.ndarray]
        :return: (getter) the potential derivative tensors with trailing unset orders trimmed, or `None` if unavailable
        :rtype: list[np.ndarray] | None
        """
        ...

    @potential_derivatives.setter
    def potential_derivatives(self, derivs):
        """
        **LLM Docstring**

        Property getter/setter for the potential-energy derivative tensors. The getter fetches them from `self.potential_surface`, normalizing missing/placeholder entries to `0` and trimming any trailing zero-padded (unset) higher-order terms; the setter assigns to `self.potential_surface.derivatives`.

        :param derivs: (setter only) the new potential derivative tensors
        :type derivs: list[np.ndarray]
        :return: (getter) the potential derivative tensors with trailing unset orders trimmed, or `None` if unavailable
        :rtype: list[np.ndarray] | None
        """
        ...

    def get_cartesian_potential_derivatives(self, order=None, evaluator=None, use_cached=True):
        """
        **LLM Docstring**

        Fetch the potential-energy derivatives in Cartesian coordinates, computing them via `calculate_energy` (and caching the result on the molecule, if the same evaluator is configured as the default) if not already available to the requested order.

        :param order: the highest derivative order needed; defaults to `2` if a fresh calculation is required
        :type order: int | None
        :param evaluator: an explicit energy evaluator to use instead of `self.energy_evaluator`
        :type evaluator: object | None
        :param use_cached: accepted for interface consistency but not used in this method's body
        :type use_cached: bool
        :return: the potential derivative tensors, truncated to `order` if given, or `None` if unavailable
        :rtype: list[np.ndarray] | None
        """
        ...

    def get_internal_potential_derivatives(self, order=None, reembed=True, strip_embedding=True, zero_gradient=False):
        """
        **LLM Docstring**

        Fetch the potential-energy derivatives re-expressed in internal coordinates, by re-expanding the Cartesian potential derivatives through the Cartesians-by-internals Jacobian.

        :param order: the highest derivative order needed
        :type order: int | None
        :param reembed: whether to use the Eckart-reembedded Cartesians-by-internals Jacobian
        :type reembed: bool
        :param strip_embedding: whether to strip the fixed embedding coordinates from the Jacobian
        :type strip_embedding: bool
        :param zero_gradient: whether to zero out the first-order (gradient) term before re-expanding
        :type zero_gradient: bool
        :return: the internal-coordinate potential derivative tensors, or `None` if unavailable
        :rtype: list[np.ndarray] | None
        """
        ...

    @property
    def normal_modes(self):
        """
        :return:
        :rtype: NormalModesManager
        """
        ...

    @normal_modes.setter
    def normal_modes(self, val):
        """
        **LLM Docstring**

        Setter for the normal modes: wraps the given value in a `NormalModesManager` (via `NormalModesManager.from_data`) unless it already is one. (A separate getter, not shown here, returns the current modes.)

        :param val: the new normal modes, in any form accepted by `NormalModesManager.from_data`
        :type val: object
        :return: None
        :rtype: None
        """
        ...

    def get_normal_modes(self, masses=None, potential_derivatives=None, use_internals=None, project_transrot=True, **opts):
        """
        **LLM Docstring**

        Compute this molecule's normal modes (via `NormalModes.from_molecule`), optionally using alternate masses/potential derivatives and controlling whether internal coordinates and translation/rotation projection are used.

        :param masses: masses to use instead of this molecule's own
        :type masses: np.ndarray | None
        :param potential_derivatives: potential derivatives to use instead of this molecule's own
        :type potential_derivatives: list[np.ndarray] | None
        :param use_internals: whether to compute the modes in internal coordinates rather than Cartesians
        :type use_internals: bool | None
        :param project_transrot: whether to project out translational/rotational degrees of freedom
        :type project_transrot: bool
        :param opts: extra options forwarded to `NormalModes.from_molecule`
        :type opts: dict
        :return: the computed normal modes
        :rtype: NormalModes
        """
        ...

    def get_reaction_path_modes(self, masses=None, potential_derivatives=None, **opts):
        """
        **LLM Docstring**

        Compute reaction-path-following normal modes for this molecule, via `ReactionPathModes.from_molecule`.

        :param masses: masses to use instead of this molecule's own
        :type masses: np.ndarray | None
        :param potential_derivatives: potential derivatives to use instead of this molecule's own
        :type potential_derivatives: list[np.ndarray] | None
        :param opts: extra options forwarded to `ReactionPathModes.from_molecule`
        :type opts: dict
        :return: the computed reaction-path modes
        :rtype: ReactionPathModes
        """
        ...

    @property
    def metadata(self):
        """
        **LLM Docstring**

        Property getter/setter for the molecule's metadata dict. The setter requires the new value to already be a `dict`.

        :param val: (setter only) the new metadata dict
        :type val: dict
        :return: (getter) the metadata dict
        :rtype: dict
        :raises TypeError: if the setter is given something that isn't a `dict`
        """
        ...

    @metadata.setter
    def metadata(self, val):
        """
        **LLM Docstring**

        Property getter/setter for the molecule's metadata dict. The setter requires the new value to already be a `dict`.

        :param val: (setter only) the new metadata dict
        :type val: dict
        :return: (getter) the metadata dict
        :rtype: dict
        :raises TypeError: if the setter is given something that isn't a `dict`
        """
        ...

    def get_harmonic_spectrum(self, **opts):
        """
        **LLM Docstring**

        Build a harmonic IR spectrum for this molecule, via `HarmonicSpectrum.from_mol`.

        :param opts: extra options forwarded to `HarmonicSpectrum.from_mol`
        :type opts: dict
        :return: the constructed harmonic spectrum
        :rtype: HarmonicSpectrum
        """
        ...

    def get_harmonic_raman_spectrum(self, **opts):
        """
        **LLM Docstring**

        Build a harmonic Raman spectrum for this molecule, via `HarmonicSpectrum.raman_from_mol`.

        :param opts: extra options forwarded to `HarmonicSpectrum.raman_from_mol`
        :type opts: dict
        :return: the constructed harmonic Raman spectrum
        :rtype: HarmonicSpectrum
        """
        ...

    def __repr__(self):
        """
        **LLM Docstring**

        Debug string representation using the molecule's name (falling back to the source-file basename, then an SMILES string from any attached RDKit molecule, then the chemical formula) and its coordinate shape.

        :return: string of the form `ClassName('name', shape=coord_shape)`
        :rtype: str
        """
        ...

    @property
    def num_atoms(self):
        """
        **LLM Docstring**

        The number of atoms in the molecule.

        :return: the atom count
        :rtype: int
        """
        ...

    @property
    def atom_positions(self):
        """
        **LLM Docstring**

        A mapping from element symbol to the list of atom indices having that symbol.

        :return: the symbol-to-indices mapping
        :rtype: dict
        """
        ...

    @property
    def dummy_positions(self):
        """
        **LLM Docstring**

        The indices of any dummy (`'X'`) atoms in the molecule.

        :return: the dummy-atom indices, or an empty list if there are none
        :rtype: list[int]
        """
        ...

    @property
    def atoms(self):
        """
        **LLM Docstring**

        The element symbols of the molecule's atoms.

        :return: the atom symbols
        :rtype: tuple[str]
        """
        ...

    def _atomic_masses(self):
        """
        **LLM Docstring**

        The atomic masses converted to atomic units (electron masses) if they appear to be given in amu (heuristically, if the smallest mass is below 100).

        :return: the masses in atomic units
        :rtype: np.ndarray
        """
        ...

    @property
    def atomic_masses(self):
        """
        **LLM Docstring**

        The atomic masses in atomic units (electron masses), via `_atomic_masses`.

        :return: the masses in atomic units
        :rtype: np.ndarray
        """
        ...

    @property
    def bonds(self):
        """
        **LLM Docstring**

        Property getter/setter for the molecule's bonds. The getter lazily guesses the bonds (via `get_guessed_bonds`) if none are set and `self.guess_bonds` is enabled.

        :param b: (setter only) the new bonds
        :type b: list[tuple] | None
        :return: (getter) the bonds, or `None` if unset and bond-guessing is disabled
        :rtype: list[tuple] | None
        """
        ...

    @bonds.setter
    def bonds(self, b):
        """
        **LLM Docstring**

        Property getter/setter for the molecule's bonds. The getter lazily guesses the bonds (via `get_guessed_bonds`) if none are set and `self.guess_bonds` is enabled.

        :param b: (setter only) the new bonds
        :type b: list[tuple] | None
        :return: (getter) the bonds, or `None` if unset and bond-guessing is disabled
        :rtype: list[tuple] | None
        """
        ...

    def break_bonds(self, bonds, use_rdkit=False, **rdopts):
        """
        **LLM Docstring**

        Build a copy of this molecule with the specified bonds removed, either by delegating to the attached RDKit molecule's `break_bonds` or by filtering the bond list directly.

        :param bonds: the bonds to remove, each an atom-index pair
        :type bonds: Iterable[tuple]
        :param use_rdkit: whether to perform the break via the RDKit molecule instead of filtering `self.bonds` directly
        :type use_rdkit: bool
        :param rdopts: extra options forwarded to the RDKit `break_bonds` call
        :type rdopts: dict
        :return: the new molecule with the given bonds removed
        :rtype: Molecule
        """
        ...

    @property
    def formula(self):
        """
        **LLM Docstring**

        The molecule's chemical formula, via `self.prop('chemical_formula')`.

        :return: the chemical formula
        :rtype: str
        """
        ...

    @property
    def multiconfig(self):
        """
        **LLM Docstring**

        Whether this molecule holds multiple geometry configurations at once, delegating to `self.coords.multiconfig`.

        :return: whether multiple configurations are stored
        :rtype: bool
        """
        ...

    @property
    def name(self):
        """
        **LLM Docstring**

        The molecule's name, falling back to `"Unnamed"` if none was set.

        :return: the molecule's name
        :rtype: str
        """
        ...

    @property
    def source_file(self):
        """
        **LLM Docstring**

        Setter for the source file: if given a plain path string, infers the source `mode` from its extension and wraps both into a `{'file':..., 'mode':...}` dict; an already-structured dict is stored as-is.

        :param src: the new source file, as a path string or a `{'file', 'mode'}` dict
        :type src: str | dict | None
        :return: None
        :rtype: None
        """
        ...

    @source_file.setter
    def source_file(self, src):
        """
        **LLM Docstring**

        Setter for the source file: if given a plain path string, infers the source `mode` from its extension and wraps both into a `{'file':..., 'mode':...}` dict; an already-structured dict is stored as-is.

        :param src: the new source file, as a path string or a `{'file', 'mode'}` dict
        :type src: str | dict | None
        :return: None
        :rtype: None
        """
        ...

    @property
    def source_mode(self):
        """
        **LLM Docstring**

        The inferred/stored source-file mode (e.g. `'fchk'`, `'log'`), if a source file is set.

        :return: the source mode, or `None` if no source file is set or no mode was recorded
        :rtype: str | None
        """
        ...

    @property
    def shape(self):
        """
        **LLM Docstring**

        The shape of the molecule's Cartesian coordinates.

        :return: the coordinate array's shape
        :rtype: tuple[int]
        """
        ...

    def __len__(self):
        """
        **LLM Docstring**

        The number of geometry configurations held by this molecule: the leading dimension of its coordinates if `multiconfig`, otherwise `1`.

        :return: the number of configurations
        :rtype: int
        """
        ...

    def copy(self):
        """
        **LLM Docstring**

        Make a copy of this molecule by calling `modify()` with no overrides.

        :return: the copied molecule
        :rtype: Molecule
        """
        ...

    def take_submolecule(self, pos):
        """
        **LLM Docstring**

        Build a new molecule containing only the atoms at the given positions, remapping bonds to the new (sub)indexing and dropping any bonds that reference atoms outside the subset.

        :param pos: the atom indices to keep, in the desired order for the submolecule
        :type pos: Iterable[int]
        :return: the constructed submolecule
        :rtype: Molecule
        """
        ...

    def prop(self, name, *args, **kwargs):
        """
        **LLM Docstring**

        Compute a named derived molecular property by dispatching to the corresponding function on `MolecularProperties`.

        :param name: the property name, matching an attribute of `MolecularProperties`
        :type name: str
        :param args: positional arguments forwarded to the property function
        :type args: tuple
        :param kwargs: keyword arguments forwarded to the property function
        :type kwargs: dict
        :return: the computed property value
        :rtype: object
        :raises MolecularPropertyError: if `name` doesn't match a known property
        """
        ...
    bond_guessing_mode = 'rdkit'

    def get_guessed_bonds(self, mode=None, **opts):
        """
        **LLM Docstring**

        Guess the bonding arrangement for this molecule, either via RDKit (from the Cartesian geometry) or via `MolecularProperties.guessed_bonds`, depending on `mode`.

        :param mode: the bond-guessing strategy to use (`'rdkit'` or another mode understood by `MolecularProperties.guessed_bonds`); defaults to `self.bond_guessing_mode`
        :type mode: str | None
        :param opts: extra options forwarded to the underlying bond-guessing routine
        :type opts: dict
        :return: the guessed bonds
        :rtype: list[tuple]
        """
        ...

    @property
    def edge_graph(self) -> EdgeGraph:
        """
        **LLM Docstring**

        The (cached) `EdgeGraph` representation of the molecule's bonding structure, built lazily via `MolecularProperties.edge_graph`.

        :return: the edge graph
        :rtype: EdgeGraph
        """
        ...

    def find_path(self, atom1, atom2):
        """
        **LLM Docstring**

        Find a path between two atoms through the bonding graph.

        :param atom1: the starting atom index
        :type atom1: int
        :param atom2: the ending atom index
        :type atom2: int
        :return: the path between the two atoms
        :rtype: list[int]
        """
        ...

    def find_substructure(self, pattern):
        """
        **LLM Docstring**

        Find matches of a substructure pattern within the molecule, delegating to the attached RDKit molecule.

        :param pattern: the substructure pattern to search for (e.g. a SMARTS string)
        :type pattern: str
        :return: the matching substructures
        :rtype: object
        """
        ...

    def apply_smarts(self, pattern):
        """
        **LLM Docstring**

        Apply a SMARTS reaction/transformation pattern to the molecule via RDKit, returning each resulting product as a new `Molecule`.

        :param pattern: the SMARTS pattern to apply
        :type pattern: str
        :return: the resulting molecules
        :rtype: list[Molecule]
        """
        ...

    def neighborhood(self, loc, size=1):
        """
        **LLM Docstring**

        Find the atoms within a given graph-distance of a location in the bonding graph.

        :param loc: the atom index to center the neighborhood on
        :type loc: int
        :param size: the neighborhood radius (in bond-graph steps)
        :type size: int
        :return: the neighboring atom indices
        :rtype: tuple[int]
        """
        ...

    def remove_hydrogens(self, positions=None, max=None, *, hydrogen_types=None):
        """
        **LLM Docstring**

        Build a copy of this molecule with hydrogen-type atoms (`H`/`D`/`T` by default) removed, either all of them or only those neighboring specified positions (optionally capped at `max` per position).

        :param positions: atom(s) whose neighboring hydrogens should be removed; if `None`, every hydrogen-type atom in the molecule is removed
        :type positions: int | Iterable[int] | None
        :param max: maximum number of hydrogens to remove per position in `positions`
        :type max: int | None
        :param hydrogen_types: the set of element symbols treated as hydrogen isotopes; defaults to `{'H', 'D', 'T'}`
        :type hydrogen_types: set[str] | None
        :return: the molecule with the selected hydrogens removed
        :rtype: Molecule
        """
        ...

    def fragment_embedding(self, fragment_indices, ref=None, return_axes=False, view_inds=(1, 2), use_moments=False):
        """
        **LLM Docstring**

        Compute a local coordinate frame (origin, offset vector, and an up-vector or full axis set) anchored at a fragment of the molecule, used as the reference frame for attaching or orienting substituents; falls back to center-of-mass/principal-axis reference points (encoded as indices `-1`/`-2`/`-3`) when the fragment doesn't have enough atoms of its own to define a frame.

        :param fragment_indices: the atom index (or indices) defining the fragment to embed
        :type fragment_indices: int | Iterable[int]
        :param ref: reference atom(s) (outside the fragment) used to anchor the origin/frame; computed from the local neighborhood if not given
        :type ref: Iterable[int] | None
        :param return_axes: whether to return a full 3x3 axis frame instead of just an up-vector
        :type return_axes: bool
        :param view_inds: which two fragment-atom positions define the "view" direction used to build the axis frame
        :type view_inds: tuple[int, int]
        :param use_moments: whether to derive the up-vector from the fragment's moments of inertia rather than from its first three atom positions
        :type use_moments: bool
        :return: `(origin, offset, up_or_axes)` -- the reference origin point, the offset from origin to the fragment's first atom, and either an up-vector or (if `return_axes`) a full axis frame
        :rtype: tuple
        """
        ...

    def attach_functional_group(self, target_fragment, atoms, new_coords, bonds='recompute', ref=None, masses=None, distance='auto', angle=0, dihedral=0, embedding='auto', bond_order=None, use_absolue_posititions=False, group_site=None) -> 'typing.Self':
        """
        **LLM Docstring**

        Build a copy of this molecule with a new group of atoms (`atoms`/`new_coords`) attached at `target_fragment`, positioning and orienting the new group using the fragment's local reference frame (bond distance/angle/dihedral, or an explicit embedding), and splicing the corresponding bonds into the result; supports designating a `group_site` atom within the new group as its attachment point, in which case the method recurses after re-deriving the embedding/bonds relative to that site.

        :param target_fragment: the atom(s) of this molecule the new group attaches to/replaces
        :type target_fragment: int | Iterable[int]
        :param atoms: the element symbols of the atoms in the new group
        :type atoms: Iterable[str]
        :param new_coords: the (local) coordinates of the new group's atoms
        :type new_coords: np.ndarray
        :param bonds: bonds within the new group; `'recompute'` to guess them fresh, `None` to reuse `self.bonds` remapped, or an explicit bond list
        :type bonds: str | list | None
        :param ref: reference atom(s) used to anchor the attachment frame; computed automatically if not given
        :type ref: Iterable[int] | None
        :param masses: masses for the new group's atoms; looked up from `atoms` if not given
        :type masses: np.ndarray | None
        :param distance: the bond distance to place the new group at; `'auto'` to look it up from `BondData`, or `None`/a number
        :type distance: str | float | None
        :param angle: rotation angle (about the up-vector) to apply to the new group
        :type angle: float
        :param dihedral: rotation angle (about the offset axis) to apply to the new group
        :type dihedral: float
        :param embedding: the reference orientation for the new group; `'auto'` to derive it from moments of inertia, or an explicit `(origin, axes)`/axes specification
        :type embedding: str | tuple | np.ndarray | None
        :param bond_order: the bond order connecting the new group to the target fragment; defaults to `1` (or inferred when `group_site` is used)
        :type bond_order: float | None
        :param use_absolue_posititions: whether `new_coords` should be used as absolute coordinates rather than being repositioned relative to the fragment frame
        :type use_absolue_posititions: bool
        :param group_site: index (within `atoms`/`new_coords`) of the atom that should serve as the attachment point; if given, the method recurses with the group re-anchored at this site and that atom excluded from the final group
        :type group_site: int | None
        :return: the molecule with the new group attached
        :rtype: Molecule
        """
        ...

    def find_heavy_atom_backbone(self, root=None):
        """
        **LLM Docstring**

        Find the longest chain of atoms in the bonding graph (the heavy-atom backbone), via `edge_graph.find_longest_chain`.

        :param root: an atom index to force as one end of the chain
        :type root: int | None
        :return: the backbone atom indices, in chain order
        :rtype: list[int]
        """
        ...

    def find_backbone_segments(self, root=None, initial_backbone=None):
        """
        **LLM Docstring**

        Split the bonding graph into backbone-connected segments, via `edge_graph.segment_by_chains`.

        :param root: an atom index to anchor the segmentation at
        :type root: int | None
        :param initial_backbone: an initial backbone chain to build the segmentation around
        :type initial_backbone: Iterable[int] | None
        :return: the resulting chain segments
        :rtype: list
        """
        ...

    def get_backbone_zmatrix(self, root=None, segments=None, return_remainder=False, return_segments=False, required_coordinates=None, isolated_coordinates=None, root_coordinates=None, initial_backbone=None, validate=True):
        """
        **LLM Docstring**

        Build a Z-matrix for a (typically single-fragment) molecule by first segmenting its bonding graph into backbone chains (via `find_backbone_segments`, validating there are no duplicate atoms across segments), constructing the base Z-matrix graph from the bonds and segments, filling in any bonds missing from the initial graph, and (if requested) enforcing required/isolated/root coordinate constraints.

        :param root: an atom index to anchor the backbone segmentation at
        :type root: int | None
        :param segments: precomputed backbone segments to use instead of calling `find_backbone_segments`
        :type segments: list | None
        :param return_remainder: whether to also return the bonds that had to be added beyond the base backbone graph
        :type return_remainder: bool
        :param return_segments: whether to also return the backbone segments used
        :type return_segments: bool
        :param required_coordinates: internal coordinates that must appear in the resulting Z-matrix
        :type required_coordinates: Iterable | None
        :param isolated_coordinates: coordinates that must be built in isolation from others
        :type isolated_coordinates: Iterable | None
        :param root_coordinates: coordinates to anchor at the root of the Z-matrix
        :type root_coordinates: Iterable | None
        :param initial_backbone: an initial backbone chain to seed the segmentation with
        :type initial_backbone: Iterable[int] | None
        :param validate: whether to validate the Z-matrix construction at each step (duplicate atoms, valid additions)
        :type validate: bool
        :return: the Z-matrix, or a tuple additionally including the segments and/or remainder bonds depending on the `return_*` flags
        :rtype: np.ndarray | tuple
        :raises ValueError: if `validate` is set and duplicate atoms are found across backbone segments
        """
        ...

    def get_canonical_zmatrix(self, ordering=None, validate=True):
        """
        **LLM Docstring**

        Build a canonical Z-matrix ordering for the molecule from a canonical fragmentation of the bonding graph.

        :param ordering: the atom ordering to use as the basis for canonicalization; defaults to the natural `0..N` ordering
        :type ordering: np.ndarray | None
        :param validate: whether to validate each Z-matrix addition
        :type validate: bool
        :return: the canonical Z-matrix
        :rtype: np.ndarray
        """
        ...

    @classmethod
    def _filter_coordinates_by_fragments(cls, inds, frags, required_coordinates):
        """
        **LLM Docstring**

        Split a list of required internal-coordinate specs into those fully contained within a single fragment (reindexed to that fragment's local atom numbering) versus those that span multiple fragments and must be merged in afterward.

        :param inds: the atom-index lists defining each fragment
        :type inds: list[Iterable[int]]
        :param frags: the corresponding submolecules for each fragment (unused directly in the body but kept for interface/length consistency)
        :type frags: list[Molecule]
        :param required_coordinates: the coordinate specs to classify, each a tuple of atom indices in the full-molecule numbering
        :type required_coordinates: Iterable[tuple] | None
        :return: `(merge_coordinates, fragment_requireds)` -- the coordinates that couldn't be assigned to a single fragment, and a per-fragment list of the (locally reindexed) coordinates assigned to each fragment
        :rtype: tuple[list, list]
        """
        ...

    def get_bond_zmatrix(self, fragments=None, segments=None, root=None, required_coordinates=None, isolated_coordinates=None, root_coordinates=None, attachment_points=None, check_attachment_points=True, validate=True, for_fragment=None, fragment_ordering=None, connect_fragments=True, initial_backbone=None):
        """
        **LLM Docstring**

        Build a full Z-matrix for the molecule from its bonding graph, handling the single-fragment case via `get_backbone_zmatrix` directly and the multi-fragment case by building a per-fragment Z-matrix for each fragment (optionally reordering/rooting/filtering required coordinates per fragment) and then splicing them together into one connected Z-matrix (via `coordops.complex_zmatrix`) unless `connect_fragments` is `False`; can also be restricted to build the Z-matrix for just one fragment (`for_fragment`), in which case it recurses on the corresponding submolecule and reindexes the result back to the full atom numbering.

        :param fragments: explicit fragment atom-index groups to use instead of `self.fragment_indices`
        :type fragments: list[Iterable[int]] | None
        :param segments: precomputed backbone segments for the single-fragment case
        :type segments: list | None
        :param root: root atom(s) to anchor the Z-matrix construction at (per fragment, in the multi-fragment case)
        :type root: int | Iterable | None
        :param required_coordinates: internal coordinates that must appear in the resulting Z-matrix
        :type required_coordinates: Iterable | None
        :param isolated_coordinates: coordinates that must be built in isolation from others
        :type isolated_coordinates: Iterable | None
        :param root_coordinates: coordinates to anchor at the root of the Z-matrix
        :type root_coordinates: Iterable | None
        :param attachment_points: explicit inter-fragment attachment points to use when connecting fragments
        :type attachment_points: dict | Iterable | None
        :param check_attachment_points: whether to validate the attachment points used when connecting fragments
        :type check_attachment_points: bool
        :param validate: whether to validate each Z-matrix addition
        :type validate: bool
        :param for_fragment: restrict the Z-matrix construction to just this fragment (an index into `self.fragment_indices`, or an explicit list of atom indices), returning the result reindexed to the full molecule
        :type for_fragment: int | Iterable[int] | None
        :param fragment_ordering: explicit ordering to apply to the fragments before connecting them
        :type fragment_ordering: Iterable[int] | None
        :param connect_fragments: whether to splice the per-fragment Z-matrices into one connected Z-matrix, or return them as a list of separate Z-matrices
        :type connect_fragments: bool
        :param initial_backbone: an initial backbone chain to seed the segmentation of each fragment with
        :type initial_backbone: Iterable[int] | None
        :return: the (connected) Z-matrix, or a list of per-fragment Z-matrices if `connect_fragments` is `False`
        :rtype: np.ndarray | list[np.ndarray]
        """
        ...

    @property
    def fragment_indices(self):
        """
        **LLM Docstring**

        The (cached) grouping of atom indices into connected molecular fragments, computed lazily via `MolecularProperties.fragment_indices`.

        :return: the list of per-fragment atom-index groups
        :rtype: list[np.ndarray]
        """
        ...

    @property
    def fragments(self):
        """
        **LLM Docstring**

        The molecule split into its connected fragments (as separate `Molecule` objects), via `MolecularProperties.fragments`.

        :return: the list of fragment molecules
        :rtype: list[Molecule]
        """
        ...

    @property
    def mass_weighted_coords(self):
        """
        :return:
        :rtype: CoordinateSet
        """
        ...

    @property
    def center_of_mass(self):
        """
        :return:
        :rtype: CoordinateSet
        """
        ...

    @property
    def inertia_tensor(self):
        """
        :return:
        :rtype: (np.ndarray, np.ndarray)
        """
        ...

    @property
    def inertial_eigensystem(self):
        """
        :return:
        :rtype: (np.ndarray, np.ndarray)
        """
        ...

    @property
    def moments_of_inertia(self):
        """
        :return:
        :rtype: np.ndarray
        """
        ...

    @property
    def inertial_axes(self):
        """
        :return:
        :rtype: np.ndarray
        """
        ...

    @property
    def translation_rotation_modes(self):
        """
        :return:
        :rtype: np.ndarray
        """
        ...

    def get_translation_rotation_projector(self, mass_weighted=False):
        """
        **LLM Docstring**

        Build the projector matrix that removes translational and rotational displacement components from a Cartesian displacement.

        :param mass_weighted: whether the projector should act on mass-weighted coordinates
        :type mass_weighted: bool
        :return: the translation/rotation projector
        :rtype: np.ndarray
        """
        ...

    def get_translation_rotation_invariant_transformation(self, mass_weighted=False, strip_embedding=True):
        """
        **LLM Docstring**

        Build the transformation (and its inverse) that projects out the translational and rotational degrees of freedom from this molecule's Cartesian coordinates.

        :param mass_weighted: whether the transformation should act on mass-weighted coordinates
        :type mass_weighted: bool
        :param strip_embedding: whether to strip the fixed embedding coordinates from the result
        :type strip_embedding: bool
        :return: the translation/rotation-invariant transformation and its inverse
        :rtype: tuple
        """
        ...

    def _load_energy(self):
        """
        **LLM Docstring**

        Load the molecule's total energy from its source file, currently only supported for Gaussian FChk files.

        :return: the total energy, or `None` if there is no source file or it isn't a supported format
        :rtype: float | None
        """
        ...

    @property
    def energy(self):
        """
        **LLM Docstring**

        The molecule's total energy: loaded from its source file if possible, else computed via `calculate_energy` if an energy evaluator is available, and cached thereafter.

        :return: the total energy, or `None` if neither a source-file energy nor an evaluator is available
        :rtype: float | None
        """
        ...
    default_energy_evalutor = 'rdkit'

    def get_energy_evaluator(self, evaluator=None, **opts):
        """
        **LLM Docstring**

        Resolve (and, if needed, instantiate from this molecule) an energy-evaluator object, defaulting to `self.energy_evaluator` (then `self.default_energy_evalutor`) if none is given explicitly.

        :param evaluator: an explicit evaluator (or evaluator-type specification) to resolve
        :type evaluator: object | None
        :param opts: extra options forwarded to the evaluator's `from_mol` constructor, if applicable
        :type opts: dict
        :return: the resolved energy-evaluator instance
        :rtype: object
        :raises ValueError: if the evaluator type can't be resolved
        """
        ...

    def get_energy_function(self, evaluator=None, *, order=None, **opts):
        """
        **LLM Docstring**

        Build a standalone callable that evaluates this molecule's energy (and derivatives) at arbitrary coordinates using the resolved energy evaluator.

        :param evaluator: an explicit evaluator (or evaluator-type specification) to use
        :type evaluator: object | None
        :param order: the derivative order the returned function should evaluate to by default
        :type order: int | None
        :param opts: extra options forwarded to `get_energy_evaluator`
        :type opts: dict
        :return: a function `(coords, order=order) -> energy_or_expansion`
        :rtype: callable
        """
        ...

    def calculate_energy(self, coords=None, *, evaluator=None, order=None, **opts):
        """
        **LLM Docstring**

        Compute the energy (and, optionally, its derivatives) at the given (or this molecule's own) coordinates using the resolved energy evaluator.

        :param coords: the coordinates to evaluate at, in Bohr; defaults to `self.coords`
        :type coords: np.ndarray | None
        :param evaluator: an explicit evaluator (or evaluator-type specification) to use
        :type evaluator: object | None
        :param order: the highest derivative order to compute; if `None`, only the energy itself is returned
        :type order: int | None
        :param opts: extra options forwarded to `get_energy_evaluator`
        :type opts: dict
        :return: the energy, or the full energy/derivative expansion
        :rtype: float | list
        """
        ...

    def partial_force_field(self, coords=None, modes=None, *, evaluator=None, order=4, mesh_spacing=1, analytic_derivative_order=None, **opts):
        """
        **LLM Docstring**

        Compute a partial (mode-selected) force-field expansion of the potential, evaluated in the given normal modes via the resolved energy evaluator's `partial_force_field` method.

        :param coords: the coordinates to evaluate at; defaults to `self.coords`
        :type coords: np.ndarray | None
        :param modes: the normal modes to expand in; computed via `get_normal_modes` if not given
        :type modes: object | None
        :param evaluator: an explicit evaluator (or evaluator-type specification) to use
        :type evaluator: object | None
        :param order: the highest derivative order to compute
        :type order: int
        :param mesh_spacing: finite-difference step size to use for the underlying force-field evaluation
        :type mesh_spacing: float
        :param analytic_derivative_order: order up to which derivatives should be computed analytically rather than numerically
        :type analytic_derivative_order: int | None
        :param opts: extra options forwarded to `get_energy_evaluator`
        :type opts: dict
        :return: the partial force-field expansion
        :rtype: object
        """
        ...

    def optimize(self, evaluator=None, *, method=None, tol=None, max_iterations=None, logger=None, reembed=True, **opts):
        """
        **LLM Docstring**

        Geometry-optimize the molecule using the resolved energy evaluator, splitting `opts` into optimizer-specific and evaluator-specific options, re-embedding the optimized (and any trajectory) coordinates back onto the original geometry via an Eckart embedding, and returning a modified copy of the molecule at the optimized geometry.

        :param evaluator: an explicit evaluator (or evaluator-type specification) to use
        :type evaluator: object | None
        :param method: the optimization method/algorithm to use
        :type method: str | None
        :param tol: convergence tolerance for the optimizer
        :type tol: float | None
        :param max_iterations: maximum number of optimizer iterations
        :type max_iterations: int | None
        :param logger: logger to report optimization progress/warnings to
        :type logger: Logger | str | None
        :param reembed: whether to Eckart-reembed the optimized geometry (and trajectory) onto the original one
        :type reembed: bool
        :param opts: extra options, split between optimizer options and evaluator-construction options
        :type opts: dict
        :return: the optimized molecule, or `(optimized_molecule, trajectory)` if the optimizer returned a trajectory
        :rtype: Molecule | tuple[Molecule, list]
        """
        ...

    def relaxed_scan(self, scan_values, scan_coordinates, evaluator=None, *, method=None, tol=None, max_iterations=None, logger=None, reembed=False, **opts):
        """
        **LLM Docstring**

        Perform a relaxed (constrained) potential-energy scan over the given coordinates/values using the resolved energy evaluator, optionally Eckart-reembedding the resulting trajectory onto the original geometry.

        :param scan_values: the coordinate value(s) to scan over
        :type scan_values: object
        :param scan_coordinates: the coordinate(s) to hold fixed/scan along
        :type scan_coordinates: object
        :param evaluator: an explicit evaluator (or evaluator-type specification) to use
        :type evaluator: object | None
        :param method: the optimization method/algorithm to use for the constrained relaxations
        :type method: str | None
        :param tol: convergence tolerance for the optimizer
        :type tol: float | None
        :param max_iterations: maximum number of optimizer iterations
        :type max_iterations: int | None
        :param logger: logger to report progress/warnings to
        :type logger: Logger | str | None
        :param reembed: whether to Eckart-reembed the resulting trajectory onto the original geometry
        :type reembed: bool
        :param opts: extra options, split between optimizer/scan options and evaluator-construction options
        :type opts: dict
        :return: `(opts, traj, meta)` -- the resolved options, the resulting geometry trajectory, and scan metadata
        :rtype: tuple
        """
        ...

    def get_dipole_evaluator(self, evaluator=None, **opts):
        """
        **LLM Docstring**

        Resolve (and, if needed, instantiate from this molecule) a dipole-evaluator object, defaulting to `self.dipole_evaluator` if none is given explicitly.

        :param evaluator: an explicit evaluator (or evaluator-type specification) to resolve
        :type evaluator: object | None
        :param opts: extra options forwarded to the evaluator's `from_mol` constructor, if applicable
        :type opts: dict
        :return: the resolved dipole-evaluator instance
        :rtype: object
        :raises ValueError: if the evaluator type can't be resolved
        """
        ...

    def get_dipole_function(self, evaluator=None, *, order=None, **opts):
        """
        **LLM Docstring**

        Build a standalone callable that evaluates this molecule's dipole (and derivatives) at arbitrary coordinates using the resolved dipole evaluator.

        :param evaluator: an explicit evaluator (or evaluator-type specification) to use
        :type evaluator: object | None
        :param order: the derivative order the returned function should evaluate to by default
        :type order: int | None
        :param opts: extra options forwarded to `get_dipole_evaluator`
        :type opts: dict
        :return: a function `(coords, order=order) -> dipole_or_expansion`
        :rtype: callable
        """
        ...

    def calculate_dipole(self, evaluator=None, order=None, **opts):
        """
        **LLM Docstring**

        Compute the dipole moment (and, optionally, its derivatives) at this molecule's current coordinates using the resolved dipole evaluator.

        :param evaluator: an explicit evaluator (or evaluator-type specification) to use
        :type evaluator: object | None
        :param order: the highest derivative order to compute; if `None`, only the dipole itself is returned
        :type order: int | None
        :param opts: extra options forwarded to `get_dipole_evaluator`
        :type opts: dict
        :return: the dipole, or the full dipole/derivative expansion
        :rtype: np.ndarray | list
        """
        ...

    def get_polarizability_evaluator(self, evaluator=None, **opts):
        """
        **LLM Docstring**

        Resolve (and, if needed, instantiate from this molecule) a dipole-polarizability-evaluator object, defaulting to `self.polarizability_evaluator` if none is given explicitly.

        :param evaluator: an explicit evaluator (or evaluator-type specification) to resolve
        :type evaluator: object | None
        :param opts: extra options forwarded to the evaluator's `from_mol` constructor, if applicable
        :type opts: dict
        :return: the resolved polarizability-evaluator instance
        :rtype: object
        :raises ValueError: if the evaluator type can't be resolved
        """
        ...

    def get_dipole_polarizability_function(self, evaluator=None, *, order=None, **opts):
        """
        **LLM Docstring**

        Build a standalone callable that evaluates this molecule's dipole polarizability (and derivatives) at arbitrary coordinates using the resolved polarizability evaluator.

        :param evaluator: an explicit evaluator (or evaluator-type specification) to use
        :type evaluator: object | None
        :param order: the derivative order the returned function should evaluate to by default
        :type order: int | None
        :param opts: extra options forwarded to `get_polarizability_evaluator`
        :type opts: dict
        :return: a function `(coords, order=order) -> polarizability_or_expansion`
        :rtype: callable
        """
        ...

    def calculate_dipole_polarizability(self, evaluator=None, order=None, **opts):
        """
        **LLM Docstring**

        Compute the dipole polarizability (and, optionally, its derivatives) at this molecule's current coordinates using the resolved polarizability evaluator.

        :param evaluator: an explicit evaluator (or evaluator-type specification) to use
        :type evaluator: object | None
        :param order: the highest derivative order to compute; if `None`, only the zeroth-order term of each returned quantity is kept
        :type order: int | None
        :param opts: extra options forwarded to `get_polarizability_evaluator`
        :type opts: dict
        :return: the polarizability expansion
        :rtype: list
        """
        ...

    @property
    def polarizability_derivatives(self):
        """
        **LLM Docstring**

        Property getter/setter for the dipole-polarizability derivative tensors, delegating to `self.dipole_surface.polarizability_derivatives`.

        :param derivs: (setter only) the new polarizability derivative tensors
        :type derivs: list[np.ndarray]
        :return: (getter) the polarizability derivative tensors
        :rtype: list[np.ndarray]
        """
        ...

    @polarizability_derivatives.setter
    def polarizability_derivatives(self, derivs):
        """
        **LLM Docstring**

        Property getter/setter for the dipole-polarizability derivative tensors, delegating to `self.dipole_surface.polarizability_derivatives`.

        :param derivs: (setter only) the new polarizability derivative tensors
        :type derivs: list[np.ndarray]
        :return: (getter) the polarizability derivative tensors
        :rtype: list[np.ndarray]
        """
        ...

    def get_reduced_potential_generator(self):
        """
        **LLM Docstring**

        Build a `ReducedDimensionalPotentialHandler` for generating reduced-dimensionality potential slices of this molecule.

        :return: the constructed handler
        :rtype: ReducedDimensionalPotentialHandler
        """
        ...

    def get_1d_potentials(self, spec, evaluator=None, energy_expansion=None, potential_params=None, **opts):
        """
        **LLM Docstring**

        Generate 1D potential-energy slices along the specified coordinate(s), via a `ReducedDimensionalPotentialHandler`.

        :param spec: the coordinate specification(s) to slice along
        :type spec: object
        :param evaluator: an explicit energy evaluator to use for the potential evaluations
        :type evaluator: object | None
        :param energy_expansion: a precomputed energy expansion to use instead of evaluating fresh
        :type energy_expansion: list[np.ndarray] | None
        :param potential_params: extra parameters controlling the potential generation
        :type potential_params: dict | None
        :param opts: extra options forwarded to `ReducedDimensionalPotentialHandler.get_1d_potentials`
        :type opts: dict
        :return: the generated 1D potentials
        :rtype: object
        """
        ...

    def evaluate(self, func, use_internals=None, order=None, strip_embedding=False):
        """
        **LLM Docstring**

        Evaluate an arbitrary function of the molecule's coordinates (and its derivatives) via the molecule's `MolecularEvaluator`.

        :param func: the function to evaluate
        :type func: callable
        :param use_internals: whether to evaluate/differentiate in internal coordinates rather than Cartesians
        :type use_internals: bool | None
        :param order: the highest derivative order to compute
        :type order: int | None
        :param strip_embedding: whether to strip the fixed embedding coordinates from the result
        :type strip_embedding: bool
        :return: the evaluated function value(s)/derivatives
        :rtype: object
        """
        ...

    def evaluate_at(self, func, coords, use_internals=None, order=None, strip_embedding=False):
        """
        **LLM Docstring**

        Evaluate an arbitrary function of the molecule's coordinates (and its derivatives) at a specific set of coordinates, via the molecule's `MolecularEvaluator`.

        :param func: the function to evaluate
        :type func: callable
        :param coords: the coordinates to evaluate at
        :type coords: np.ndarray
        :param use_internals: whether to evaluate/differentiate in internal coordinates rather than Cartesians
        :type use_internals: bool | None
        :param order: the highest derivative order to compute
        :type order: int | None
        :param strip_embedding: whether to strip the fixed embedding coordinates from the result
        :type strip_embedding: bool
        :return: the evaluated function value(s)/derivatives
        :rtype: object
        """
        ...

    def get_displaced_coordinates(self, displacements, which=None, sel=None, axes=None, use_internals=False, coordinate_expansion=None, strip_embedding=False, shift=True):
        """
        **LLM Docstring**

        Build displaced copies of the molecule's coordinates along specified atoms/axes, via the molecule's `MolecularEvaluator`.

        :param displacements: the displacement values to apply
        :type displacements: np.ndarray
        :param which: which atoms/coordinates to displace
        :type which: object | None
        :param sel: a selection of atoms to restrict the displacement to
        :type sel: Iterable[int] | None
        :param axes: which axes to displace along
        :type axes: object | None
        :param use_internals: whether the displacements are given in internal coordinates rather than Cartesians
        :type use_internals: bool
        :param coordinate_expansion: a coordinate-transformation expansion to apply to the displacements before applying them
        :type coordinate_expansion: list[np.ndarray] | None
        :param strip_embedding: whether the coordinate expansion has had its embedding coordinates stripped
        :type strip_embedding: bool
        :param shift: whether the displacements are relative shifts (`True`) or absolute target values
        :type shift: bool
        :return: the displaced coordinates
        :rtype: np.ndarray
        """
        ...

    def get_scan_coordinates(self, domains, internals=False, modes=None, order=None, which=None, sel=None, axes=None, shift=True, coordinate_expansion=None, strip_embedding=False, return_displacements=False):
        """
        **LLM Docstring**

        Build a grid of displaced coordinates spanning the given coordinate domains, optionally re-expressed through a normal-mode (or other) coordinate expansion, via the molecule's `MolecularEvaluator`.

        :param domains: the coordinate ranges to scan over
        :type domains: object
        :param internals: whether the scan coordinates are internal coordinates rather than Cartesians
        :type internals: bool
        :param modes: normal modes to scan along instead of raw coordinates; `True` to use this molecule's own modes (with translation/rotation projection disabled)
        :type modes: object | bool | None
        :param order: the Jacobian order to use when building the mode-based coordinate expansion
        :type order: int | None
        :param which: which atoms/coordinates to scan
        :type which: object | None
        :param sel: a selection of atoms to restrict the scan to
        :type sel: Iterable[int] | None
        :param axes: which axes to scan along
        :type axes: object | None
        :param shift: whether the scan values are relative shifts or absolute target values
        :type shift: bool
        :param coordinate_expansion: an explicit coordinate-transformation expansion to combine with the mode-based one (if `modes` is given) or use directly
        :type coordinate_expansion: list[np.ndarray] | None
        :param strip_embedding: whether the coordinate expansion has had its embedding coordinates stripped
        :type strip_embedding: bool
        :param return_displacements: whether to also return the raw displacement values used
        :type return_displacements: bool
        :return: the scan coordinates, or additionally the displacements if `return_displacements` is set
        :rtype: np.ndarray | tuple
        """
        ...

    def get_nearest_displacement_atoms(self, points, sel=None, axes=None, weighting_function=None, return_distances=False):
        """
        **LLM Docstring**

        Find the atoms nearest to a set of query points (optionally restricted to a selection/axes, with a custom distance weighting), via the molecule's `MolecularEvaluator`.

        :param points: the query points to find nearest atoms for
        :type points: np.ndarray
        :param sel: a selection of atoms to restrict the search to
        :type sel: Iterable[int] | None
        :param axes: which coordinate axes to compute distances over
        :type axes: object | None
        :param weighting_function: a custom function to weight distances by
        :type weighting_function: callable | None
        :param return_distances: whether to also return the computed distances
        :type return_distances: bool
        :return: the nearest atom indices, or additionally the distances if `return_distances` is set
        :rtype: np.ndarray | tuple
        """
        ...

    def get_nearest_displacement_coordinates(self, points, sel=None, axes=None, weighting_function=None, modes_nearest=False, return_distances=False):
        """
        **LLM Docstring**

        Find the displaced coordinates nearest to a set of query points, optionally in normal-mode space, via the molecule's `MolecularEvaluator`.

        :param points: the query points to find nearest coordinates for
        :type points: np.ndarray
        :param sel: a selection of atoms to restrict the search to
        :type sel: Iterable[int] | None
        :param axes: which coordinate axes to compute distances over
        :type axes: object | None
        :param weighting_function: a custom function to weight distances by
        :type weighting_function: callable | None
        :param modes_nearest: whether to find the nearest point in normal-mode space rather than Cartesian/internal space
        :type modes_nearest: bool
        :param return_distances: whether to also return the computed distances
        :type return_distances: bool
        :return: the nearest coordinates, or additionally the distances if `return_distances` is set
        :rtype: np.ndarray | tuple
        """
        ...

    def get_nearest_scan_coordinates(self, domains, sel=None, axes=None):
        """
        **LLM Docstring**

        Find the scan-grid coordinates nearest to a set of query domains, via the molecule's `MolecularEvaluator`.

        :param domains: the coordinate domains to build the nearest scan grid for
        :type domains: object
        :param sel: a selection of atoms to restrict the search to
        :type sel: Iterable[int] | None
        :param axes: which coordinate axes to compute distances over
        :type axes: object | None
        :return: the nearest scan coordinates
        :rtype: np.ndarray
        """
        ...

    @classmethod
    def _get_atomic_radius(cls, atom_data, radius_type=None):
        """
        **LLM Docstring**

        Resolve the display radius for an atom-data record: the icon radius (falling back to the van der Waals radius if too small), or an explicitly named radius field.

        :param atom_data: the atom-data record to get a radius for
        :type atom_data: dict
        :param radius_type: the specific radius field to use; if `None`, uses the icon radius with a van der Waals fallback
        :type radius_type: str | None
        :return: the resolved atomic radius
        :rtype: float
        """
        ...

    def plot_molecule_function(self, function, *, axes, sel=None, embed=False, modes_nearest=False, domain=None, domain_padding=1, plot_points=500, weighting_function=None, mask_function=None, mask_value=0, plot_atoms=False, atom_colors=None, atom_radii=None, plotter=None, epilog=None, **plot_options):
        """
        **LLM Docstring**

        Plot an arbitrary scalar function of the molecule's geometry over a 1D or 2D grid of displaced coordinates, evaluating the function at the grid points nearest to each displacement, optionally masking out-of-range values and overlaying the atoms' positions.

        :param function: the function to evaluate, taking a batch of coordinates and returning scalar values
        :type function: callable
        :param axes: which atom/coordinate axis (or axes, up to 2) to build the plot grid over
        :type axes: int | Iterable[int]
        :param sel: a selection of atoms to restrict the displacement search to
        :type sel: Iterable[int] | None
        :param embed: whether to Eckart-embed the evaluation points before calling `function`
        :type embed: bool
        :param modes_nearest: whether to find the nearest displacement in normal-mode space
        :type modes_nearest: bool
        :param domain: explicit `(min, max)` bounds per axis; computed from the geometry's bounding box (plus `domain_padding`) if not given
        :type domain: np.ndarray | None
        :param domain_padding: padding to add around the automatically computed domain
        :type domain_padding: float | Iterable[float] | None
        :param plot_points: number of grid points per axis
        :type plot_points: int | Iterable[int]
        :param weighting_function: custom distance-weighting function used when finding nearest displacements
        :type weighting_function: callable | None
        :param mask_function: a function `(values, eval_points, dists) -> mask` used to blank out certain grid values
        :type mask_function: callable | None
        :param mask_value: the value to substitute where `mask_function` returns `True`
        :type mask_value: float
        :param plot_atoms: whether to overlay disks at the atoms' projected positions
        :type plot_atoms: bool
        :param atom_colors: per-atom overlay colors; defaults to each atom's icon color
        :type atom_colors: list | None
        :param atom_radii: per-atom overlay radii; defaults to each atom's resolved radius
        :type atom_radii: list | None
        :param plotter: the plotting class/function to use; defaults to `plt.Plot` (1D) or `plt.TriContourPlot` (2D)
        :type plotter: type | callable | None
        :param epilog: extra plot elements to draw on top of the function plot
        :type epilog: list | None
        :param plot_options: extra options forwarded to the plotting call
        :type plot_options: dict
        :return: the resulting plot
        :rtype: object
        :raises NotImplementedError: if the molecule holds more than one geometry configuration
        :raises ValueError: if more than 2 axes are requested
        """
        ...

    def get_model(self, potential_specs, dipole=None):
        """
        **LLM Docstring**

        Build an analytic `MolecularModel` (a symbolic potential/dipole surface expressed in terms of this molecule's Z-matrix internal coordinates) from a specification of per-coordinate potential (and optional dipole) function contributions.

        :param potential_specs: either a raw potential expression, or a dict mapping internal-coordinate index (or index tuples, for coupled terms) to a function-type specification (a dict with a `'function_type'` key or a single `{function_type: params}` entry) describing that coordinate's contribution
        :type potential_specs: dict | object
        :param dipole: an optional 3-element (x/y/z) sequence of dipole-component specifications, each following the same per-coordinate specification format as `potential_specs`
        :type dipole: Iterable | None
        :return: the constructed analytic molecular model
        :rtype: MolecularModel
        :raises ValueError: if the molecule has no internal coordinates, if its internal-coordinate spec isn't a plain Z-matrix, or if `dipole` isn't given exactly 3 components
        """
        ...

    def get_point_group(self, coords=None, masses=None, *, sel=None, verbose=False, return_identifier=False, **tols):
        """
        **LLM Docstring**

        Identify the molecular point group at the given (or this molecule's own) geometry, via `McUtils.Symmetry.PointGroupIdentifier`.

        :param coords: the coordinates to identify the point group for; defaults to `self.coords`
        :type coords: np.ndarray | None
        :param masses: the masses to use; defaults to `self.atomic_masses`
        :type masses: np.ndarray | None
        :param sel: a selection of atoms to restrict the identification to
        :type sel: Iterable[int] | None
        :param verbose: whether the identifier should print diagnostic output
        :type verbose: bool
        :param return_identifier: whether to also return the underlying `PointGroupIdentifier` object
        :type return_identifier: bool
        :param tols: extra tolerance options forwarded to `PointGroupIdentifier`
        :type tols: dict
        :return: the identified point group, or `(identifier, point_group)` if `return_identifier` is set
        :rtype: object | tuple
        """
        ...

    @property
    def point_group(self) -> symm.PointGroup:
        """
        **LLM Docstring**

        Property getter/setter for the molecule's (cached) point group. The getter computes it lazily via `get_point_group` the first time it's needed.

        :param point_group: (setter only) the new point group to cache
        :type point_group: symm.PointGroup
        :return: (getter) the cached (or newly computed) point group
        :rtype: symm.PointGroup
        """
        ...

    @point_group.setter
    def point_group(self, point_group):
        """
        **LLM Docstring**

        Property getter/setter for the molecule's (cached) point group. The getter computes it lazily via `get_point_group` the first time it's needed.

        :param point_group: (setter only) the new point group to cache
        :type point_group: symm.PointGroup
        :return: (getter) the cached (or newly computed) point group
        :rtype: symm.PointGroup
        """
        ...

    def get_point_group_embedded_coordinates(self, pg=None, sel=None, return_point_group=False, return_identifier=False, **tols):
        """
        **LLM Docstring**

        Compute this molecule's coordinates re-embedded into the standard reference frame of its molecular point group (translating to the center of mass, aligning to the principal axes, then aligning to the point group's own axis convention).

        :param pg: an explicit point group to embed into; identified automatically via `get_point_group` if not given
        :type pg: symm.PointGroup | None
        :param sel: a selection of atoms to restrict the point-group identification to
        :type sel: Iterable[int] | None
        :param return_point_group: whether to also return the (embedded) point group
        :type return_point_group: bool
        :param return_identifier: whether to also return the underlying `PointGroupIdentifier`
        :type return_identifier: bool
        :param tols: extra tolerance options forwarded to `PointGroupIdentifier`
        :type tols: dict
        :return: the point-group-embedded coordinates, or additionally the point group and/or identifier depending on the `return_*` flags
        :rtype: np.ndarray | tuple
        """
        ...

    def symmetrize(self, pg=None, return_identifier=False, tol=0.1, sel=None, return_coordinates=None, return_point_group=False, **tols):
        """
        **LLM Docstring**

        Symmetrize the molecule's geometry to conform exactly to a (identified or given) point group, adding/relabeling atoms as needed if the symmetrization changes the atom count, and matching the symmetrized geometry back onto the original atom ordering when the count is unchanged.

        :param pg: an explicit point group to symmetrize to; identified automatically via `get_point_group` if not given
        :type pg: symm.PointGroup | None
        :param return_identifier: whether to also return the underlying `PointGroupIdentifier`
        :type return_identifier: bool
        :param tol: tolerance used for point-group identification and symmetrization
        :type tol: float
        :param sel: a selection of atoms to restrict the symmetrization to
        :type sel: Iterable[int] | None
        :param return_coordinates: whether to return raw coordinates instead of a modified `Molecule` (when the atom count is unchanged)
        :type return_coordinates: bool | None
        :param return_point_group: whether to also return the point group used
        :type return_point_group: bool
        :param tols: extra tolerance options forwarded to `PointGroupIdentifier`
        :type tols: dict
        :return: a tuple starting with the symmetrized molecule/coordinates (or `None` if the atom count changed and no direct match was possible), followed by `(new_atoms, new_coords)`, and optionally the identifier and/or point group
        :rtype: tuple
        """
        ...

    def get_symmetrized_internals(self, point_group=None, *, internals=None, extra_internals=None, masses=None, return_expansions=False, atom_selection=None, as_characters=True, normalize=None, drop_empty_modes=None, perms=None, return_base_expansion=False, return_point_group=False, reduce_redundant_coordinates=None, ops=None, permutation_tol=0.01, **opts):
        """
        **LLM Docstring**

        Build symmetry-adapted linear combinations of internal coordinates for the molecule under a (identified or given) point group, via `McUtils.Symmetry.symmetrize_internals`.

        :param point_group: an explicit point group to symmetrize under; identified automatically via `get_point_group` (or on a submolecule, if `atom_selection` is given) if not given
        :type point_group: symm.PointGroup | None
        :param internals: the internal coordinates to symmetrize; extracted from `self.internals` (Z-matrix or generic spec) if not given
        :type internals: Iterable | None
        :param extra_internals: extra internal coordinates to prepend to the extracted set (deduplicated)
        :type extra_internals: Iterable | None
        :param masses: masses to use; defaults to `self.atomic_masses`
        :type masses: np.ndarray | None
        :param return_expansions: whether to also return the coordinate-derivative expansions used
        :type return_expansions: bool
        :param atom_selection: restrict the point-group identification to a submolecule
        :type atom_selection: Iterable[int] | None
        :param as_characters: whether to express the symmetrized coordinates in terms of irrep characters
        :type as_characters: bool
        :param normalize: whether to normalize the symmetrized coordinate combinations
        :type normalize: bool | None
        :param drop_empty_modes: whether to drop symmetry combinations that come out identically zero
        :type drop_empty_modes: bool | None
        :param perms: explicit atom permutations to use for the symmetry operations
        :type perms: Iterable | None
        :param return_base_expansion: whether to also return the base (pre-symmetrization) coordinate expansion
        :type return_base_expansion: bool
        :param return_point_group: whether to also return the point group used
        :type return_point_group: bool
        :param reduce_redundant_coordinates: whether to reduce a redundant coordinate set down to a non-redundant symmetrized basis
        :type reduce_redundant_coordinates: bool | None
        :param ops: explicit symmetry operations to use instead of the full point-group operation set
        :type ops: Iterable | None
        :param permutation_tol: tolerance used when matching atom permutations to symmetry operations
        :type permutation_tol: float
        :param opts: extra options forwarded to point-group identification/`symmetrize_internals`
        :type opts: dict
        :return: the symmetrized internal coordinates, plus any additionally requested return values, prefixed with the point group if `return_point_group` is set
        :rtype: object | tuple
        :raises ValueError: if no internal coordinates are available or given
        """
        ...

    def get_surface(self, radius_type='VanDerWaalsRadius', *, surface_type=None, radius_units='Angstroms', samples=100, radius_scaling=1, **etc):
        """
        **LLM Docstring**

        Build a surface representation of the molecule (by default, a union of per-atom spheres) using per-atom radii of the requested type.

        :param radius_type: the `AtomData` radius field to use for each atom
        :type radius_type: str
        :param surface_type: the surface class to construct; defaults to `zach.SphereUnionSurface`
        :type surface_type: type | None
        :param radius_units: the units the radii are given in before converting to Bohr
        :type radius_units: str
        :param samples: number of samples used to represent the surface
        :type samples: int
        :param radius_scaling: uniform scale factor applied to all radii
        :type radius_scaling: float
        :param etc: extra options forwarded to the surface constructor
        :type etc: dict
        :return: the constructed surface
        :rtype: zach.SphereUnionSurface
        """
        ...

    def get_surface_mesh(self, radius_type='VanDerWaalsRadius', *, surface_type=None, radius_units='Angstroms', samples=50, expansion=0.01, mesh_options=None, **etc):
        """
        **LLM Docstring**

        Build a triangulated mesh of the molecule's surface, via `get_surface` followed by `generate_mesh`.

        :param radius_type: the `AtomData` radius field to use for each atom
        :type radius_type: str
        :param surface_type: the surface class to construct; defaults to `zach.SphereUnionSurface`
        :type surface_type: type | None
        :param radius_units: the units the radii are given in before converting to Bohr
        :type radius_units: str
        :param samples: number of samples used to represent the surface
        :type samples: int
        :param expansion: expansion factor applied to the surface before meshing
        :type expansion: float
        :param mesh_options: extra options forwarded to `generate_mesh`
        :type mesh_options: dict | None
        :param etc: extra options forwarded to `get_surface`
        :type etc: dict
        :return: the generated surface mesh
        :rtype: object
        """
        ...

    def setup_AIMD(self, potential_function=None, timestep=0.5, seed=None, total_energy=None, total_energy_scaling=None, trajectories=1, sampled_modes=None, initial_energies=None, initial_displacements=None, initial_mode_directions=None, displaced_coords=None, track_kinetic_energy=False, track_velocities=False, modes=None, **etc):
        """
        **LLM Docstring**

        Set up an `AIMDSimulator` (ab initio / classical molecular dynamics) for this molecule, either propagating from an explicit set of displaced initial positions, or sampling initial normal-mode velocities/energies (from explicit mode directions, explicit per-trajectory energies, or randomly-sampled directions distributing a target total energy) when starting from the equilibrium geometry.

        :param potential_function: the energy/gradient function to use; built via `get_energy_function` if not given
        :type potential_function: callable | None
        :param timestep: the simulation timestep
        :type timestep: float
        :param seed: random seed for sampling initial mode directions
        :type seed: int | None
        :param total_energy: the total vibrational energy to distribute across sampled trajectories/modes
        :type total_energy: float | None
        :param total_energy_scaling: scale factor applied to the default total energy (sum of mode frequencies) if `total_energy` isn't given
        :type total_energy_scaling: float | None
        :param trajectories: number of trajectories to sample when randomly generating mode directions
        :type trajectories: int
        :param sampled_modes: which normal modes to sample energy into; defaults to all of them
        :type sampled_modes: Iterable[int] | None
        :param initial_energies: explicit per-trajectory, per-mode energies to use instead of sampling them
        :type initial_energies: np.ndarray | None
        :param initial_displacements: explicit initial Cartesian displacements to start the trajectories from, bypassing the normal-mode energy-sampling path entirely
        :type initial_displacements: np.ndarray | None
        :param initial_mode_directions: explicit per-trajectory mode-direction vectors to convert into initial energies
        :type initial_mode_directions: np.ndarray | None
        :param displaced_coords: which coordinates `initial_displacements` applies to
        :type displaced_coords: object | None
        :param track_kinetic_energy: whether the simulator should track kinetic energy over the trajectory
        :type track_kinetic_energy: bool
        :param track_velocities: whether the simulator should track velocities over the trajectory
        :type track_velocities: bool
        :param modes: normal modes to use for the energy-to-velocity conversion; computed via `get_normal_modes` if not given
        :type modes: object | None
        :param etc: extra options forwarded to the `AIMDSimulator` constructor
        :type etc: dict
        :return: the constructed simulator
        :rtype: AIMDSimulator
        :raises ValueError: if both `initial_energies` and `initial_mode_directions` are given
        """
        ...

    def setup_VPT(self, *, states=2, order=2, use_internals=None, potential_derivatives=None, energy_evaluator=None, dipole_derivatives=None, dipole_evaluator=None, runner='matrix', use_reaction_path=False, modes=None, projected_modes=None, mode_transformation=None, potential_terms=None, dipole_terms=None, **opts):
        """
        **LLM Docstring**

        Set up a vibrational perturbation theory (VPT) calculation for this molecule, resolving the potential/dipole derivative expansions (computing them via the configured evaluators if not already available), the normal modes (optionally localized against reaction-path-projected modes), and dispatching to the appropriate VPT runner (`VPTRunner` for the default matrix-based approach, or `AnalyticVPTRunner`).

        :param states: the target vibrational states specification (e.g. max quanta) to compute
        :type states: int | object
        :param order: the perturbation-theory order to compute the potential/dipole expansions to
        :type order: int
        :param use_internals: whether to run VPT in internal coordinates (constructing from `self.modify()`) rather than bare atoms/coordinates
        :type use_internals: bool | None
        :param potential_derivatives: explicit potential-energy derivative tensors to use instead of computing them
        :type potential_derivatives: list[np.ndarray] | None
        :param energy_evaluator: an explicit energy evaluator to use when computing the potential derivatives
        :type energy_evaluator: object | None
        :param dipole_derivatives: explicit dipole derivative tensors to use instead of computing them
        :type dipole_derivatives: list[np.ndarray] | None
        :param dipole_evaluator: an explicit dipole evaluator to use when computing the dipole derivatives
        :type dipole_evaluator: object | None
        :param runner: which VPT runner to use: `'matrix'` for `VPTRunner`, or an explicit runner class/object
        :type runner: str | object
        :param use_reaction_path: whether to project the normal modes against reaction-path modes
        :type use_reaction_path: bool
        :param modes: explicit normal modes to use instead of computing them
        :type modes: object | None
        :param projected_modes: explicit modes to localize the normal modes against
        :type projected_modes: object | None
        :param mode_transformation: an explicit mode-coordinate transformation to use
        :type mode_transformation: object | None
        :param potential_terms: precomputed potential expansion terms, bypassing derivative computation entirely
        :type potential_terms: object | None
        :param dipole_terms: precomputed dipole expansion terms, bypassing derivative computation entirely
        :type dipole_terms: object | None
        :param opts: extra options forwarded to the VPT runner's `construct` method
        :type opts: dict
        :return: the constructed VPT runner/calculation
        :rtype: object
        """
        ...

    def setup_job(self, job_type, *args, **kwargs):
        """
        **LLM Docstring**

        Set up an external-program computational job (e.g. a Gaussian/ORCA job) for this molecule, via `ExternalProgramJob.resolve` and `from_mol`.

        :param job_type: the job type to resolve (a name or job-type object)
        :type job_type: str | object
        :param args: positional arguments forwarded to the resolved job type's `from_mol`
        :type args: tuple
        :param kwargs: keyword arguments forwarded to the resolved job type's `from_mol`, merged over its default options
        :type kwargs: dict
        :return: the constructed job
        :rtype: object
        """
        ...

    def get_gmatrix(self, masses=None, coords=None, use_internals=None, power=None, **internals_opts):
        """
        **LLM Docstring**

        Compute the G-matrix (inverse effective-mass metric tensor) for this molecule, either the trivial diagonal inverse-mass form (in plain Cartesians) or the projected `B G0 B^T`-style form built from the internals-by-Cartesians Jacobian (when using internal coordinates), optionally raised to a fractional power.

        :param masses: masses to use instead of this molecule's own
        :type masses: np.ndarray | None
        :param coords: alternate coordinates to compute the G-matrix at
        :type coords: np.ndarray | None
        :param use_internals: whether to compute the G-matrix in internal coordinates; defaults to `True` if the molecule has internal coordinates defined
        :type use_internals: bool | None
        :param power: an optional power to raise the resulting G-matrix to
        :type power: float | None
        :param internals_opts: extra options forwarded to `get_internals_by_cartesians` when using internal coordinates
        :type internals_opts: dict
        :return: the G-matrix
        :rtype: np.ndarray
        """
        ...

    @property
    def g_matrix(self):
        """
        Returns the molecular g-matrix for the system
        :return:
        :rtype:
        """
        ...

    @property
    def coriolis_constants(self):
        """
        Returns the molecular g-matrix for the system
        :return:
        :rtype:
        """
        ...

    def principle_axis_frame(self, sel=None, inverse=False):
        """
        Gets the principle axis frame(s) for the molecule
        :param mol:
        :type mol:
        :param sel: selection of atoms to use when getting the Eckart frame
        :type sel:
        :param inverse: whether to return the inverse of the rotations or not
        :type inverse: bool
        :return:
        :rtype: MolecularTransformation | List[MolecularTransformation]
        """
        ...

    @property
    def principle_axis_data(self):
        """
        Gets the principle axis embedded coords and embedding parameters for the molecule

        :return:
        :rtype: MolecularTransformation | List[MolecularTransformation]
        """
        ...

    def permute_atoms(self, perm):
        """
        **LLM Docstring**

        Build a copy of this molecule with its atoms reordered according to a permutation, remapping bonds to the new indexing accordingly.

        :param perm: the new atom ordering, given as the sequence of old indices in their new order
        :type perm: Iterable[int]
        :return: the permuted molecule
        :rtype: Molecule
        """
        ...

    def apply_affine_transformation(self, transformation, load_properties=False, embed_properties=True):
        """
        **LLM Docstring**

        Apply an affine (rotation/translation) transformation to the molecule's coordinates, optionally propagating the transformation onto its normal modes, potential surface, and dipole surface as well.

        :param transformation: the affine transformation to apply, either a raw transformation matrix or an object exposing an `apply` method (e.g. a `MolecularTransformation`)
        :type transformation: np.ndarray | object
        :param load_properties: whether to force-load the normal modes/potential/dipole derivatives before transforming (so they get carried over even if not already computed); `None` to load only if cheaply available
        :type load_properties: bool | None
        :param embed_properties: whether to propagate the transformation onto already-loaded (or newly loaded) normal modes/potential/dipole surfaces
        :type embed_properties: bool
        :return: the transformed molecule (with `source_file` cleared)
        :rtype: Molecule
        """
        ...

    def apply_rotation(self, rotation_matrix, shift_com=None, load_properties=None, embed_properties=True):
        """
        **LLM Docstring**

        Apply a rotation (optionally shifted to rotate about the center of mass) to the molecule, via `apply_affine_transformation`.

        :param rotation_matrix: the rotation matrix (or full affine matrix) to apply
        :type rotation_matrix: np.ndarray
        :param shift_com: whether to first shift so the rotation is performed about the center of mass (only applied if `rotation_matrix` isn't already a 4x4 affine matrix)
        :type shift_com: bool | None
        :param load_properties: forwarded to `apply_affine_transformation`
        :type load_properties: bool | None
        :param embed_properties: forwarded to `apply_affine_transformation`
        :type embed_properties: bool
        :return: the rotated molecule
        :rtype: Molecule
        """
        ...

    def eckart_frame(self, mol, sel=None, inverse=False, planar_ref_tolerance=None, proper_rotation=False, reset_com=True):
        """
        Gets the Eckart frame(s) for the molecule

        :param mol:
        :type mol:
        :param sel: selection of atoms to use when getting the Eckart frame
        :type sel:
        :param inverse: whether to return the inverse of the rotations or not
        :type inverse: bool
        :return:
        :rtype: MolecularTransformation
        """
        ...

    def embed_coords(self, crds, sel=None, in_paf=False, planar_ref_tolerance=None, proper_rotation=False):
        """
        Embeds coords in the Eckart frame using `self` as a reference

        :param crds:
        :type crds:
        :return:
        :rtype:
        """
        ...

    def get_embedding_data(self, crds, sel=None, in_paf=False, planar_ref_tolerance=None, proper_rotation=False):
        """
        Gets the necessary data to embed crds in the Eckart frame using `self` as a reference

        :param crds:
        :type crds:
        :return:
        :rtype:
        """
        ...

    def get_embedded_molecule(self, ref=None, sel=None, planar_ref_tolerance=None, proper_rotation=False, embed_properties=True, load_properties=None, reset_com=True):
        """
        Returns a Molecule embedded in an Eckart frame if ref is not None, otherwise returns
        a principle-axis embedded Molecule
        :return:
        :rtype: Molecule
        """
        ...

    def get_rmsd(self, other: 'typing.Self | np.ndarray', sel=None, embed=True, embedding_sel=None, mass_weighted=False):
        """
        **LLM Docstring**

        Compute the (optionally mass-weighted) root-mean-square displacement between this molecule's geometry and another geometry/molecule, optionally Eckart-embedding both onto a common reference frame first and restricting the comparison to a subset of atoms.

        :param other: the geometry (or molecule) to compare against
        :type other: Molecule | np.ndarray
        :param sel: a selection of atoms to restrict the RMSD computation to
        :type sel: Iterable[int] | None
        :param embed: whether to Eckart-embed both geometries onto a common frame before comparing
        :type embed: bool
        :param embedding_sel: a selection of atoms to use for the embedding fit, if different from `sel`
        :type embedding_sel: Iterable[int] | None
        :param mass_weighted: whether to mass-weight the displacement before computing the norm
        :type mass_weighted: bool
        :return: the RMSD value(s)
        :rtype: float | np.ndarray
        """
        ...

    def align_molecule(self, other: 'typing.Self', reindex_bonds=True, permute_atoms=True, align_structures=True, sel=None, embed_properties=True, load_properties=False):
        """
        Aligns `other` with `self` by first finding the reindexing of the bonds of `other` that
        lead to the best graph overlap with `self`, then determining which atoms can be permuted based on their graph
        structures, then determining which permutation of equivalent atoms leads to the best agreement between the structures,
        and then finally finding the Eckart/min-RMSD transformation after this transformation has been applied

        :param other:
        :param reindex_bonds:
        :return:
        """
        ...

    @property
    def rdmol(self):
        """
        **LLM Docstring**

        The (cached) attached RDKit molecule representation, built lazily via `RDMolecule.from_mol` the first time it's needed (returning `None` if RDKit isn't available), and kept in sync with the current coordinates on subsequent access.

        :return: the RDKit molecule, or `None` if RDKit isn't available
        :rtype: object | None
        """
        ...

    def to_ase(self, **kwargs):
        """
        **LLM Docstring**

        Convert this molecule to an ASE (Atomic Simulation Environment) molecule object, via `ASEMolecule.from_mol`.

        :param kwargs: extra options forwarded to `ASEMolecule.from_mol`
        :type kwargs: dict
        :return: the constructed ASE molecule
        :rtype: object
        """
        ...

    @classmethod
    def from_ase(cls, ase_mol, **kwargs):
        """
        **LLM Docstring**

        Build a `Molecule` from an ASE molecule object.

        :param ase_mol: the ASE molecule to convert
        :type ase_mol: object
        :param kwargs: extra keyword arguments merged over the ASE molecule's metadata and forwarded to the constructor
        :type kwargs: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        ...

    @classmethod
    def from_zmat(cls, zmat, internals=None, axes=None, origin=None, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` from a Z-matrix specification (either a Z-matrix string or an explicit `(atoms, (ordering, coords))` tuple), converting the internal Z-matrix coordinates to Cartesians.

        :param zmat: the Z-matrix specification, as a string or an `(atoms, (ordering, coords))` tuple
        :type zmat: str | tuple
        :param internals: the internal-coordinate specification to store on the resulting molecule; defaults to the Z-matrix's own atom ordering
        :type internals: object | None
        :param axes: reference axes to use when converting to Cartesians
        :type axes: np.ndarray | None
        :param origin: reference origin to use when converting to Cartesians
        :type origin: np.ndarray | None
        :param opts: extra options forwarded to the constructor
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        ...

    @classmethod
    def from_openbabel(cls, mol, **opts):
        """

        :param mol:
        :type mol: pybel.mol
        :return:
        :rtype:
        """
        ...

    def get_obmol(self, **opts):
        """
        **LLM Docstring**

        Convert this molecule to an OpenBabel molecule object, via `OBMolecule.from_mol`.

        :param opts: extra options forwarded to `OBMolecule.from_mol`
        :type opts: dict
        :return: the constructed OpenBabel molecule
        :rtype: object
        """
        ...

    @classmethod
    def _from_log_file(cls, file, num=None, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` from the final Cartesian geometry parsed out of a Gaussian `.log` file.

        :param file: path to the Gaussian log file
        :type file: str
        :param num: which structure(s) to parse from the file (forwarded to the log reader)
        :type num: int | None
        :param opts: extra options forwarded to the constructor
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        ...

    @classmethod
    def _from_fchk_file(cls, file, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` from a Gaussian FChk file, using its atomic numbers and integer atomic weights to build isotope-labeled atom symbols and its real atomic weights as masses.

        :param file: path to the Gaussian FChk file
        :type file: str
        :param opts: extra options forwarded to the constructor
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        ...

    @classmethod
    def _from_orca_file(cls, file, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` from the final Cartesian geometry (in atomic units) parsed out of an ORCA log file.

        :param file: path to the ORCA log file
        :type file: str
        :param opts: extra options forwarded to the constructor
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        ...

    @classmethod
    def _from_molpro_file(cls, file, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` from the final Cartesian geometry parsed out of a MOLPRO log file.

        :param file: path to the MOLPRO log file
        :type file: str
        :param opts: extra options forwarded to the constructor
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        ...

    @classmethod
    def _from_hess_file(cls, file, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` from the atoms/coordinates/masses parsed out of an ORCA `.hess` file.

        :param file: path to the ORCA `.hess` file
        :type file: str
        :param opts: extra options forwarded to the constructor
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        ...

    @classmethod
    def from_rdmol(cls, rdmol, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` from an RDKit molecule (or a raw RDKit `Mol`/owning-mol object, which is first wrapped in an `RDMolecule`), carrying over its atoms, coordinates, bonds, charges, and metadata.

        :param rdmol: the RDKit molecule (or wrappable RDKit object) to convert
        :type rdmol: object
        :param opts: extra options merged over the RDKit molecule's metadata and forwarded to the constructor
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        ...

    @classmethod
    def _from_smiles(cls, smi, add_implicit_hydrogens=True, num_confs=None, conf_id=None, take_min=None, optimize=False, sanitize=False, parse_name=True, allow_cxsmiles=True, strict_cxsmiles=True, remove_hydrogens=False, replacements=None, parser_options=None, confgen_opts=None, coords=None, **opts):
        """
        **LLM Docstring**

        Build one or more `Molecule` objects from a SMILES string, generating 3D conformer(s) via RDKit and converting each resulting RDKit molecule via `from_rdmol`.

        :param smi: the SMILES string to parse
        :type smi: str
        :param add_implicit_hydrogens: whether to add implicit hydrogens
        :type add_implicit_hydrogens: bool
        :param num_confs: number of conformers to generate
        :type num_confs: int | None
        :param conf_id: which conformer id to select
        :type conf_id: int | None
        :param take_min: whether to keep only the lowest-energy conformer
        :type take_min: bool | None
        :param optimize: whether to geometry-optimize the generated conformer(s)
        :type optimize: bool
        :param sanitize: whether to sanitize the parsed RDKit molecule
        :type sanitize: bool
        :param parse_name: whether to parse a name out of the SMILES string
        :type parse_name: bool
        :param allow_cxsmiles: whether to allow CXSMILES extensions
        :type allow_cxsmiles: bool
        :param strict_cxsmiles: whether to strictly enforce CXSMILES syntax
        :type strict_cxsmiles: bool
        :param remove_hydrogens: whether to remove explicit hydrogens after parsing
        :type remove_hydrogens: bool
        :param replacements: SMILES replacement/abbreviation dict
        :type replacements: dict | None
        :param parser_options: extra low-level parser options forwarded to `RDMolecule.from_smiles`
        :type parser_options: dict | None
        :param confgen_opts: extra conformer-generation options
        :type confgen_opts: dict | None
        :param coords: explicit coordinates to use instead of generating conformers
        :type coords: np.ndarray | None
        :param opts: extra options forwarded to `from_rdmol`
        :type opts: dict
        :return: a single molecule, or a list of molecules if multiple conformers were generated
        :rtype: Molecule | list[Molecule]
        """
        ...

    @classmethod
    def _from_name(cls, name, api_key=None, add_implicit_hydrogens=True, method='pubchem', **opts):
        """
        **LLM Docstring**

        Build one or more `Molecule` objects by looking up a compound name via a chemical database API (PubChem by default, or ChemSpider) to get its SMILES string(s), then parsing those via `_from_smiles`.

        :param name: the compound name to look up
        :type name: str
        :param api_key: API key for the chosen lookup service (used for ChemSpider)
        :type api_key: str | None
        :param add_implicit_hydrogens: whether to add implicit hydrogens when parsing the resulting SMILES
        :type add_implicit_hydrogens: bool
        :param method: which lookup service to use (`'pubchem'` or `'chemspider'`)
        :type method: str
        :param opts: extra options forwarded to `_from_smiles`
        :type opts: dict
        :return: a single molecule, or a list of molecules if multiple compounds matched
        :rtype: Molecule | list[Molecule]
        :raises ValueError: if the name doesn't resolve to any known compound
        """
        ...

    @classmethod
    def _from_sdf(cls, sdf, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` from an SDF string/file, via RDKit.

        :param sdf: the SDF content to parse
        :type sdf: str
        :param opts: extra options forwarded to `from_rdmol`
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        ...

    @classmethod
    def _from_molblock(cls, sdf, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` from a MOL block string, via RDKit.

        :param sdf: the MOL block content to parse
        :type sdf: str
        :param opts: extra options forwarded to `from_rdmol`
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        ...

    @classmethod
    def _from_mol2(cls, mol, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` from a MOL2-format string/file, via RDKit.

        :param mol: the MOL2 content to parse
        :type mol: str
        :param opts: extra options forwarded to `from_rdmol`
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        ...

    @classmethod
    def _from_pdb(cls, pdb, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` from a PDB-format string/file, via RDKit.

        :param pdb: the PDB content to parse
        :type pdb: str
        :param opts: extra options forwarded to `from_rdmol`
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        ...

    @classmethod
    def _from_mrv(cls, mrv, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` from a ChemDraw/Marvin (MRV) format string/file, via RDKit.

        :param mrv: the MRV content to parse
        :type mrv: str
        :param opts: extra options forwarded to `from_rdmol`
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        ...

    @classmethod
    def _from_fasta(cls, fasta, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` from a FASTA-format sequence string, via RDKit.

        :param fasta: the FASTA content to parse
        :type fasta: str
        :param opts: extra options forwarded to `from_rdmol`
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        ...

    @classmethod
    def _from_helm(cls, helm, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` from a HELM-format string, via RDKit.

        :param helm: the HELM content to parse
        :type helm: str
        :param opts: extra options forwarded to `from_rdmol`
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        ...

    @classmethod
    def _from_cdxml(cls, cdxml, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` from a ChemDraw XML (CDXML) string/file, via RDKit.

        :param cdxml: the CDXML content to parse
        :type cdxml: str
        :param opts: extra options forwarded to `from_rdmol`
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        ...

    @classmethod
    def _from_xyz(cls, xyz, units=None, max_blocks=None, **opts):
        """
        **LLM Docstring**

        Build one or more `Molecule` objects by parsing one or more geometry blocks from XYZ-format text (or a stream/file-like object), attaching any per-block comment line as molecule metadata and converting coordinates to Bohr.

        :param xyz: the XYZ content (string, stream, or file-like object) to parse
        :type xyz: str | object
        :param units: the units the coordinates are given in; converted to Bohr (`"BohrRadius"`). If `None`, coordinates are assumed to already be in Bohr
        :type units: str | None
        :param max_blocks: maximum number of XYZ blocks (geometries) to parse; `None` parses exactly one block and returns a single `Molecule`, a negative value parses all blocks, and a positive value caps the count (returning a list)
        :type max_blocks: int | None
        :param opts: extra options forwarded to the constructor
        :type opts: dict
        :return: a single molecule (if `max_blocks` is `None`) or a list of molecules (one per parsed block)
        :rtype: Molecule | list[Molecule]
        """
        ...

    @classmethod
    def _from_xyz_file(cls, xyz_file, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` (or list of molecules) from an XYZ file, via `_from_xyz`.

        :param xyz_file: path to (or stream of) the XYZ file
        :type xyz_file: str | object
        :param opts: extra options forwarded to `_from_xyz`
        :type opts: dict
        :return: a single molecule or a list of molecules, depending on `opts['max_blocks']`
        :rtype: Molecule | list[Molecule]
        """
        ...

    @classmethod
    def _from_gspec(cls, gspec: str, charge=None, spin=None, report=None, units='Angstroms', format_options=None, **etc):
        """
        **LLM Docstring**

        Build a `Molecule` from a Gaussian "molecule specification" string (as embedded in a Gaussian job report), substituting any report-variable placeholders, normalizing formatting quirks, and dispatching to `from_string` (as a `'gspec'`-format string) or, if the string looks like a full Gaussian archive entry, to `_from_gspec_file`.

        :param gspec: the Gaussian molecule-spec string to parse
        :type gspec: str
        :param charge: the net charge; defaults to the value parsed from the spec's header if not given
        :type charge: int | None
        :param spin: the spin multiplicity; defaults to the value parsed from the spec's header if not given
        :type spin: int | None
        :param report: a parsed Gaussian job-report dict whose entries are substituted into the spec text
        :type report: dict | None
        :param units: the units the coordinates are given in
        :type units: str
        :param format_options: per-format parsing options merged with the default Z-matrix axis convention
        :type format_options: dict | None
        :param etc: extra options forwarded to `from_string`/`_from_gspec_file`
        :type etc: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        ...

    @classmethod
    def _from_gspec_file(cls, logfile, charge=None, spin=None, potential_derivatives=None, dipole_derivatives=None, target_job='Freq', **etc):
        """
        **LLM Docstring**

        Build a `Molecule` from a full Gaussian log/archive file, extracting the job report matching `target_job` (or the last report if none matches), and, for frequency jobs, the potential-energy and dipole derivative tensors (with the appropriate atomic-mass-unit conversions applied) to attach to the constructed molecule.

        :param logfile: path to (or stream of) the Gaussian log file
        :type logfile: str | object
        :param charge: the net charge to use instead of the report's own value
        :type charge: int | None
        :param spin: the spin multiplicity to use instead of the report's own value
        :type spin: int | None
        :param potential_derivatives: explicit potential derivative tensors to use instead of parsing them from the report
        :type potential_derivatives: list[np.ndarray] | None
        :param dipole_derivatives: explicit dipole derivative tensors to use instead of parsing them from the report
        :type dipole_derivatives: list[np.ndarray] | None
        :param target_job: which Gaussian job type's report to use (e.g. `'Freq'`)
        :type target_job: str
        :param etc: extra options forwarded to `from_string`
        :type etc: dict
        :return: the constructed molecule
        :rtype: Molecule
        :raises ValueError: if no job report is found in the file
        """
        ...

    @classmethod
    def from_name(cls, name, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` by looking up a compound name, via `from_string` with format `'name'`.

        :param name: the compound name to look up
        :type name: str
        :param opts: extra options forwarded to `from_string`
        :type opts: dict
        :return: the constructed molecule (or list of molecules, if multiple matches are found)
        :rtype: Molecule | list[Molecule]
        """
        ...
    _atom_strs = None

    @classmethod
    def get_atom_strings(cls):
        """
        **LLM Docstring**

        The (cached) set of up-to-2-character atomic element symbols known to `AtomData`, used for heuristically detecting SMILES/Z-matrix strings.

        :return: the set of element symbol prefixes
        :rtype: set[str]
        """
        ...

    @classmethod
    def _check_smi(cls, string, atom_types, other_syms=None):
        """
        **LLM Docstring**

        Heuristically check whether a string could be a SMILES string by stripping out all known atom-type substrings and SMILES punctuation characters and checking that nothing is left over.

        :param string: the string to check
        :type string: str
        :param atom_types: the atom-symbol substrings to strip out
        :type atom_types: Iterable[str]
        :param other_syms: the punctuation/bond-symbol characters to strip out; defaults to `cls._smi_punct`
        :type other_syms: Iterable[str] | None
        :return: whether the string is fully consumed by atom types and punctuation (i.e. looks like a SMILES string)
        :rtype: bool
        """
        ...

    @classmethod
    def _infer_str_format(cls, string: str, allow_names=False, **opts):
        """
        **LLM Docstring**

        Heuristically guess which structural-file format a raw string represents (SMILES, a compound name, MOL/SDF, Z-matrix, Gaussian molecule-spec, or XYZ), by inspecting its line count/structure and content.

        :param string: the string to classify
        :type string: str
        :param allow_names: whether a single, non-SMILES-looking word/line should be inferred as a compound name
        :type allow_names: bool
        :param opts: accepted for interface consistency but not used in this method's body
        :type opts: dict
        :return: the inferred format key (e.g. `'smi'`, `'name'`, `'mol'`, `'zmat'`, `'gspec'`, `'xyz'`)
        :rtype: str
        :raises ValueError: if no format can be inferred from the string
        """
        ...

    @classmethod
    def get_string_format_dispatchers(cls):
        """
        **LLM Docstring**

        The mapping from string-format key to the constructor method that parses that format, used by `from_string`.

        :return: the format-to-constructor mapping
        :rtype: dict
        """
        ...

    @classmethod
    def from_string(cls, string, fmt=None, allow_names=False, format_options=None, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` from a raw string in any supported structural format, inferring the format automatically if not given, dispatching to the matching in-memory parser if one exists, or otherwise writing the string to a temporary file and dispatching through `from_file`.

        :param string: the structural string to parse
        :type string: str
        :param fmt: the format to parse as; inferred via `_infer_str_format` if not given
        :type fmt: str | None
        :param allow_names: whether to allow inferring the format as a compound name
        :type allow_names: bool
        :param format_options: per-format parsing options, keyed by format, merged under `opts`
        :type format_options: dict | None
        :param opts: extra options forwarded to the format-specific parser
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        ...

    @classmethod
    def _from_ob_import(cls, file, fmt=None, **opts):
        """
        **LLM Docstring**

        Build a `Molecule` from a file via OpenBabel's format-import machinery.

        :param file: path to the file to import
        :type file: str
        :param fmt: the OpenBabel format code to use; auto-detected if not given
        :type fmt: str | None
        :param opts: extra options forwarded to `from_openbabel`
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        """
        ...

    @classmethod
    def get_file_format_dispatchers(cls):
        """
        **LLM Docstring**

        The mapping from file-format key (typically a file extension) to the constructor method that parses that format, used by `from_file`.

        :return: the format-to-constructor mapping
        :rtype: dict
        """
        ...

    @classmethod
    def from_file(cls, file, mode=None, format_options=None, use_ob_fallback=False, **opts) -> Molecule:
        """In general we'll delegate to pybel except for like Fchk and Log files

        :param file:
        :type file:
        :return:
        :rtype:
        """
        ...

    @classmethod
    def _to_smiles(cls, mol, **opts):
        """
        **LLM Docstring**

        Export a molecule to a SMILES string, via its attached RDKit molecule.

        :param mol: the molecule to export
        :type mol: Molecule
        :param opts: extra options forwarded to `rdmol.to_smiles`
        :type opts: dict
        :return: the SMILES string
        :rtype: str
        :raises ValueError: if the molecule has no attached RDKit molecule
        """
        ...

    @classmethod
    def _to_molblock(cls, mol, **opts):
        """
        **LLM Docstring**

        Export a molecule to a MOL block string, via its attached RDKit molecule.

        :param mol: the molecule to export
        :type mol: Molecule
        :param opts: extra options forwarded to `rdmol.to_molblock`
        :type opts: dict
        :return: the MOL block string
        :rtype: str
        :raises ValueError: if the molecule has no attached RDKit molecule
        """
        ...

    @classmethod
    def _to_sdf(cls, mol, **opts):
        """
        **LLM Docstring**

        Export a molecule to an SDF string, via its attached RDKit molecule.

        :param mol: the molecule to export
        :type mol: Molecule
        :param opts: extra options forwarded to `rdmol.to_sdf`
        :type opts: dict
        :return: the SDF string
        :rtype: str
        :raises ValueError: if the molecule has no attached RDKit molecule
        """
        ...

    @classmethod
    def _to_pdb(cls, mol, **opts):
        """
        **LLM Docstring**

        Export a molecule to a PDB string, via its attached RDKit molecule.

        :param mol: the molecule to export
        :type mol: Molecule
        :param opts: extra options forwarded to `rdmol.to_pdb`
        :type opts: dict
        :return: the PDB string
        :rtype: str
        :raises ValueError: if the molecule has no attached RDKit molecule
        """
        ...

    @classmethod
    def _to_cml(cls, mol, **opts):
        """
        **LLM Docstring**

        Export a molecule to a CML string, via its attached RDKit molecule.

        :param mol: the molecule to export
        :type mol: Molecule
        :param opts: extra options forwarded to `rdmol.to_cml`
        :type opts: dict
        :return: the CML string
        :rtype: str
        :raises ValueError: if the molecule has no attached RDKit molecule
        """
        ...

    @classmethod
    def _to_mrv(cls, mol, **opts):
        """
        **LLM Docstring**

        Export a molecule to a ChemDraw/Marvin (MRV) string, via its attached RDKit molecule.

        :param mol: the molecule to export
        :type mol: Molecule
        :param opts: extra options forwarded to `rdmol.to_mrv`
        :type opts: dict
        :return: the MRV string
        :rtype: str
        :raises ValueError: if the molecule has no attached RDKit molecule
        """
        ...

    @classmethod
    def _to_xyz_string(cls, mol, comment=None, units=None, num_prec=8):
        """
        **LLM Docstring**

        Format a molecule's geometry as an XYZ-format string, using a comment line (defaulting to the molecule's `repr`) and computing a column width wide enough for the largest coordinate value.

        :param mol: the molecule to export
        :type mol: Molecule
        :param comment: the comment line to use; defaults to `repr(mol)`
        :type comment: str | None
        :param units: the units to convert the coordinates into before formatting; if `None`, coordinates are left in Bohr
        :type units: str | None
        :param num_prec: number of digits after the decimal point for each coordinate
        :type num_prec: int
        :return: the formatted XYZ string
        :rtype: str
        """
        ...

    @classmethod
    def _to_zmat_string(cls, mol, units='Angstroms', float_fmt='{:8.4f}', variables=None, variable_modifications=None, **etc):
        """
        **LLM Docstring**

        Format a molecule's internal (Z-matrix) coordinates as a Z-matrix string, via `McUtils.Coordinerds.format_zmatrix_string`.

        :param mol: the molecule to export
        :type mol: Molecule
        :param units: the units to format distances in
        :type units: str
        :param float_fmt: the format string used for numeric values
        :type float_fmt: str
        :param variables: named variables to substitute for coordinate values
        :type variables: dict | None
        :param variable_modifications: modifications to apply to specific variable values
        :type variable_modifications: dict | None
        :param etc: extra options forwarded to `format_zmatrix_string`
        :type etc: dict
        :return: the formatted Z-matrix string
        :rtype: str
        :raises ValueError: if the molecule has no internal coordinates, or if they aren't a Z-matrix system
        """
        ...

    @classmethod
    def get_string_export_dispatchers(cls):
        """
        **LLM Docstring**

        The mapping from string-export format key to the exporter method that produces that format, used by `to_string`.

        :return: the format-to-exporter mapping
        :rtype: dict
        """
        ...

    def to_string(self, fmt, **opts):
        """
        **LLM Docstring**

        Export this molecule to a string in the given format, dispatching to an in-memory string exporter if one exists, otherwise round-tripping through a temporary file via `to_file`, or falling back to OpenBabel's string-export machinery.

        :param fmt: the export format key
        :type fmt: str
        :param opts: extra options forwarded to the format-specific exporter
        :type opts: dict
        :return: the exported string
        :rtype: str
        """
        ...

    @classmethod
    def get_file_export_dispatchers(cls):
        """
        **LLM Docstring**

        The mapping from file-export format key to the exporter method that writes that format to disk, used by `to_file`. Currently empty (no dedicated file exporters beyond the string-based ones and OpenBabel's fallback).

        :return: the (currently empty) format-to-exporter mapping
        :rtype: dict
        """
        ...

    def to_file(self, file, mode=None, use_ob_fallback=False, **opts):
        """
        :param file:
        :type file:
        :return:
        :rtype:
        """
        ...

    @classmethod
    def _infer_spec_format(cls, spec, **opts):
        """
        **LLM Docstring**

        Heuristically classify a raw molecule-construction `spec` (an RDKit-like object, a file path, a raw string, a dict of constructor kwargs, or an `(atoms, coords[, opts])` tuple) so `construct` can dispatch it to the right constructor.

        :param spec: the specification to classify
        :type spec: object
        :param opts: accepted for interface consistency but not used in this method's body
        :type opts: dict
        :return: `(fmt, subopts)` where `fmt` is one of `'rdmol'`, `'file'`, `'str'`, `'dict'`, or an `(atoms, coords)` pair, and `subopts` is any extra options bundled with the spec
        :rtype: tuple
        :raises ValueError: if `spec` doesn't match any recognized shape
        """
        ...

    @classmethod
    def construct(cls, spec, fmt=None, **opts):
        """
        **LLM Docstring**

        Universal `Molecule` constructor: builds a molecule from essentially any reasonable input (an existing `Molecule` to copy/modify, an RDKit/ASE object, a file path, a raw structural string, a dict of constructor kwargs, or an `(atoms, coords)`/Z-matrix pair), inferring the format automatically if not given.

        :param spec: the specification to build the molecule from
        :type spec: object
        :param fmt: an explicit format to use instead of inferring one; can also be an `(atoms, coords)` pair for direct construction
        :type fmt: str | tuple | None
        :param opts: extra options forwarded to the resolved constructor
        :type opts: dict
        :return: the constructed molecule
        :rtype: Molecule
        :raises NotImplementedError: if `spec`/`fmt` don't match any supported construction path
        """
        ...

    def plot(self, *geometries, figure=None, return_objects=False, bonds=None, bond_radius=None, atom_radius_scaling=None, atom_style=None, atom_radii=None, atom_text=None, display_atom_numbers=False, radius_type=None, bond_style=None, reconcile_bonds=True, capped_bonds=None, reflectiveness=None, vector_style=None, highlight_atoms=None, highlight_bonds=None, highlight_rings=None, highlight_styles=None, comparison_styles=None, animation_frame_styles=None, mode_vectors=None, mode_vector_origins=None, mode_vector_origin_mode='set', mode_vector_display_cutoff=0.01, principle_axes=None, principle_axes_origin=None, principle_axes_origin_mode='set', principle_axes_style=None, dipole=None, dipole_origin=None, dipole_origin_mode='set', render_multiple_bonds=None, render_fractional_bonds=None, fractional_bond_offset=None, bond_center_radius_offset=None, draw_coords=None, draw_coords_style=None, up_vector=None, multiple_bond_spacing=None, mode=None, backend=None, include_save_buttons=None, objects=False, graphics_class=None, cylinder_class=None, cylinder_options=None, sphere_class=None, sphere_options=None, arrow_class=None, arrow_options=None, line_class=None, line_options=None, disk_class=None, disk_options=None, animate=None, recording_options=None, animation_options=None, jsmol_load_script=None, include_jsmol_script_interface=None, dynamic_loading=None, units='Angstroms', label_style=None, theme='default', theme_function=None, plot_range_padding='auto', annotation_function=None, **plot_ops):
        """
        Dispatches to the appropriate `MoleculePlotter` for the resolved backend/mode.

        The full keyword surface now lives on `MoleculePlotter.plot_molecule`; calls,
        defaults, and return conventions are unchanged.
        """
        ...

    def get_animation_geoms(self, which, extent=0.35, steps=8, strip_embedding=True, units=None, coordinate_expansion=None):
        """
        **LLM Docstring**

        Build a back-and-forth looping sequence of displaced geometries for animating a single coordinate (or an arbitrary coordinate-expansion direction), via `get_scan_coordinates`.

        :param which: the coordinate index to animate (if `coordinate_expansion` is not given directly), or an explicit displacement-direction array/list of arrays (in which case this becomes the coordinate index within that expansion, defaulting to `0`)
        :type which: int | np.ndarray | list[np.ndarray]
        :param extent: the maximum displacement magnitude in each direction
        :type extent: float
        :param steps: number of steps from the equilibrium geometry out to `extent`
        :type steps: int
        :param strip_embedding: whether to strip the fixed embedding coordinates from the default Cartesians-by-internals expansion
        :type strip_embedding: bool
        :param units: units to convert the resulting geometries into; left in Bohr if `None`
        :type units: str | None
        :param coordinate_expansion: an explicit coordinate-transformation expansion to displace along, instead of the default internal-coordinate Jacobian
        :type coordinate_expansion: list[np.ndarray] | np.ndarray | None
        :return: the looping sequence of displaced geometries (out to `extent` and back)
        :rtype: np.ndarray
        """
        ...

    def animate_coordinate(self, which, extent=0.5, steps=3, return_objects=False, strip_embedding=True, units='Angstroms', backend=None, mode=None, jsmol_load_script=None, coordinate_expansion=None, **plot_opts):
        """
        **LLM Docstring**

        Build an animation of a displaced coordinate and return it as a plottable/displayable object, either as a JSMol vibration animation or as a sequence of rendered frames via `plot`.

        :param which: the coordinate index (or explicit direction) to animate, forwarded to `get_animation_geoms`/`format_animation_file`
        :type which: int | np.ndarray | list[np.ndarray]
        :param extent: the maximum displacement magnitude
        :type extent: float
        :param steps: number of steps out to `extent`
        :type steps: int
        :param return_objects: whether to return the constructed graphics objects alongside the figure (non-JSMol path only)
        :type return_objects: bool
        :param strip_embedding: whether to strip the fixed embedding coordinates from the default coordinate expansion
        :type strip_embedding: bool
        :param units: units to display the geometries in
        :type units: str
        :param backend: the rendering backend to use; defaults to `self.display_mode`
        :type backend: str | None
        :param mode: the display mode to use; defaults to `backend`
        :type mode: str | None
        :param jsmol_load_script: extra JSMol load-script text, used only in the JSMol path
        :type jsmol_load_script: str | list[str] | None
        :param coordinate_expansion: an explicit coordinate-transformation expansion to animate along
        :type coordinate_expansion: list[np.ndarray] | np.ndarray | None
        :param plot_opts: extra options forwarded to `plot`
        :type plot_opts: dict
        :return: the resulting animation figure/widget
        :rtype: object
        """
        ...

    def animate_mode(self, which, extent=0.5, steps=3, modes=None, coordinate_expansion=None, order=None, normalize=True, mass_weight=False, mass_scale=True, frequency_scale=False, **opts):
        """
        **LLM Docstring**

        Build an animation of a normal mode's displacement, converting the mode into a coordinate-expansion direction (with optional normalization, mass-weighting/scaling, and frequency scaling of the displacement extent) and delegating to `animate_coordinate`.

        :param which: the mode index to animate
        :type which: int
        :param extent: the base displacement extent, before any mass/frequency scaling
        :type extent: float
        :param steps: number of steps out to `extent`
        :type steps: int
        :param modes: the normal modes to animate; defaults to this molecule's own (converted to a fresh mode basis)
        :type modes: object | None
        :param coordinate_expansion: an additional coordinate-transformation expansion to combine with the mode-displacement direction
        :type coordinate_expansion: list[np.ndarray] | None
        :param order: the Jacobian order to use when building the mode-based coordinate expansion (only used if given)
        :type order: int | None
        :param normalize: whether to normalize the mode-displacement direction before use
        :type normalize: bool
        :param mass_weight: whether to keep the modes mass-weighted rather than removing the mass-weighting
        :type mass_weight: bool
        :param mass_scale: whether to scale the displacement extent by the mode's effective mass (only applied when not mass-weighted)
        :type mass_scale: bool
        :param frequency_scale: whether to scale the displacement extent relative to the mode's frequency
        :type frequency_scale: bool
        :param opts: extra options forwarded to `animate_coordinate`
        :type opts: dict
        :return: the resulting animation figure/widget
        :rtype: object
        """
        ...

    def _format_xyz(self, which, nat, atoms, geom, float_format='10.3f'):
        """
        **LLM Docstring**

        Format a single geometry frame as an XYZ-format block (atom count, comment line naming the frame, then one line per atom).

        :param which: the frame index/label used in the comment line
        :type which: int
        :param nat: the number of atoms
        :type nat: int
        :param atoms: the atom labels
        :type atoms: Iterable[str]
        :param geom: the Cartesian coordinates for this frame
        :type geom: np.ndarray
        :param float_format: the format spec used for each coordinate value
        :type float_format: str
        :return: the formatted XYZ block
        :rtype: str
        """
        ...

    def format_structs(self, geoms, format='xyz'):
        """
        **LLM Docstring**

        Format a batch of geometries as a concatenated multi-frame string in the given format.

        :param geoms: the geometries to format, reshaped to `(-1, natoms, 3)`
        :type geoms: np.ndarray
        :param format: the output format; currently only `'xyz'` is supported
        :type format: str
        :return: the concatenated multi-frame string
        :rtype: str
        :raises NotImplementedError: if `format` isn't `'xyz'`
        """
        ...

    def _format_jmol_displacements_files(self, coords, expansion, units='Angstroms'):
        """
        **LLM Docstring**

        Format a set of per-atom displacement vectors as JSMol-style XYZ-with-vibration blocks (base coordinates plus a zero placeholder column and the displacement vector per atom), one block per displacement direction supplied.

        :param coords: the base Cartesian coordinates
        :type coords: np.ndarray
        :param expansion: a length-1 list containing the displacement-direction tensor, reshaped to `(-1, natoms, 3)`
        :type expansion: list[np.ndarray]
        :param units: units to convert the coordinates/displacements into
        :type units: str
        :return: the list of formatted displacement blocks, one per direction
        :rtype: list[str]
        """
        ...

    def format_animation_file(self, which, format='xyz', extent=0.35, steps=8, strip_embedding=True, units='Angstroms', coordinate_expansion=None):
        """
        **LLM Docstring**

        Build a formatted animation string/block for a displaced coordinate, either as a single JSMol vibration block (`'jmol'` format) or as a looping multi-frame structure string via `get_animation_geoms`/`format_structs`.

        :param which: the coordinate index (or explicit direction) to animate
        :type which: int | np.ndarray | list[np.ndarray]
        :param format: the output format (`'jmol'` for a JSMol vibration block, otherwise forwarded to `format_structs`)
        :type format: str
        :param extent: the maximum displacement magnitude
        :type extent: float
        :param steps: number of steps out to `extent` (non-`'jmol'` formats only)
        :type steps: int
        :param strip_embedding: whether to strip the fixed embedding coordinates from the default coordinate expansion
        :type strip_embedding: bool
        :param units: units to format the geometries in
        :type units: str
        :param coordinate_expansion: an explicit coordinate-transformation expansion to animate along
        :type coordinate_expansion: list[np.ndarray] | None
        :return: the formatted animation string
        :rtype: str
        """
        ...

    def to_widget(self):
        """
        **LLM Docstring**

        Build a displayable widget for this molecule, either a JSMol applet or an X3D scene, depending on `self.display_mode`.

        :return: the constructed widget/scene object
        :rtype: object
        """
        ...
    default_display_mode = 'jsmol'

    def _ipython_display_(self):
        """
        **LLM Docstring**

        IPython/Jupyter display hook: builds the display widget via `to_widget` and displays it.

        :return: the result of displaying the widget
        :rtype: object
        """
        ...

    def _get_ob_attr(self, item):
        """
        **LLM Docstring**

        Look up an attribute on the attached `pybel` (OpenBabel) molecule, if one is set.

        :param item: the attribute name to look up
        :type item: str
        :return: the attribute's value
        :rtype: object
        :raises AttributeError: if no pybel molecule is attached
        """
        ...

class MolecoolException(Exception):
    ...