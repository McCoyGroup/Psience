from __future__ import annotations
import abc
import enum
import functools
import tempfile
import itertools
import math
import sys, os
import uuid
import time
import contextlib
import numpy as np
import warnings
import subprocess
from McUtils.Data import AtomData, UnitsData
from McUtils.Zachary import TensorDerivativeConverter, FiniteDifferenceDerivative, CoordinateFunction
import McUtils.Numputils as nput
import McUtils.Devutils as dev
import McUtils.Coordinerds as coordops
import McUtils.Iterators as itut
from McUtils.Scaffolding import Logger, read_flat_tree
from McUtils.ExternalPrograms import SingularityLauncher
from ..Data import PotentialSurface, DipoleSurface
from .CoordinateSystems import MolecularEmbedding
from .Properties import NormalModesManager
__all__ = ['MolecularEvaluator', 'EnergyEvaluator', 'RDKitEnergyEvaluator', 'AIMNet2EnergyEvaluator', 'XTBEnergyEvaluator', 'PySCFEnergyEvaluator', 'PotentialFunctionEnergyEvaluator', 'DipoleEvaluator', 'DipoleFunctionDipoleEvaluator', 'ReducedDimensionalPotentialHandler']

class MolecularEvaluator:

    def __init__(self, embedding: MolecularEmbedding, normal_modes: NormalModesManager):
        """
        **LLM Docstring**

        Store the coordinate embedding and normal-modes manager this evaluator will use to evaluate functions and build displacements.

        :param embedding: the molecular coordinate embedding to evaluate against
        :type embedding: MolecularEmbedding
        :param normal_modes: the normal-modes manager used for mode-based displacement/nearest-point queries
        :type normal_modes: NormalModesManager
        :return: None
        :rtype: None
        """
        ...

    @property
    def coords(self):
        """
        **LLM Docstring**

        The Cartesian coordinates, delegating to `self.embedding.coords`.

        :return: the Cartesian coordinates
        :rtype: CoordinateSet
        """
        ...

    @property
    def masses(self):
        """
        **LLM Docstring**

        The atomic masses, delegating to `self.embedding.masses`.

        :return: the atomic masses
        :rtype: np.ndarray
        """
        ...

    def evaluate(self, func, use_internals=None, order=None, strip_embedding=False):
        """
        **LLM Docstring**

        Evaluate a coordinate function (and, optionally, its derivatives) either directly in Cartesian coordinates or in internal coordinates with the result re-expanded back through the internals-by-Cartesians Jacobian, optionally stripping/restoring the fixed embedding coordinates around the internal-coordinate evaluation.

        :param func: the function to evaluate, expected to accept `order=` for derivative evaluation
        :type func: callable
        :param use_internals: whether to evaluate in internal coordinates; defaults to whether the embedding has internal coordinates defined
        :type use_internals: bool | None
        :param order: the highest derivative order to compute; if `None`, only the function value is returned
        :type order: int | None
        :param strip_embedding: whether to strip the fixed embedding coordinates before evaluating in internal coordinates (and restore zeroed placeholders for them in the resulting derivative tensors)
        :type strip_embedding: bool
        :return: the function value, or `[value, deriv1, deriv2, ...]` if `order` is given
        :rtype: object | list
        """
        ...

    def evaluate_at(self, func, coords, use_internals=None, order=None, strip_embedding=False):
        """
        **LLM Docstring**

        Evaluate a coordinate function at an alternate set of coordinates, by building a fresh `MolecularEvaluator` for those coordinates (keeping the same masses/internal-coordinate spec/normal modes) and calling `evaluate` on it.

        :param func: the function to evaluate
        :type func: callable
        :param coords: the coordinates to evaluate at
        :type coords: np.ndarray
        :param use_internals: whether to evaluate in internal coordinates
        :type use_internals: bool | None
        :param order: the highest derivative order to compute
        :type order: int | None
        :param strip_embedding: whether to strip the fixed embedding coordinates around the internal-coordinate evaluation
        :type strip_embedding: bool
        :return: the function value, or the value/derivative expansion if `order` is given
        :rtype: object | list
        """
        ...

    def get_displaced_coordinates(self, displacements, which=None, sel=None, axes=None, use_internals: "bool|Literal['reembed', 'convert']"=False, coordinate_expansion=None, expansion_active_positions=None, strip_embedding=False, shift=True, coords=None):
        """
        **LLM Docstring**

        Build displaced coordinates from a reference geometry by applying `displacements` -- either directly to selected flat coordinate indices/atom-axis pairs, to the full coordinate array, or (if `coordinate_expansion` is given) as a Taylor-series expansion along an arbitrary coordinate-transformation direction (e.g. a normal mode) -- optionally stripping the fixed embedding coordinates and handling the resulting internal-coordinate re-embedding/conversion convention.

        :param displacements: the displacement value(s) to apply
        :type displacements: np.ndarray
        :param which: explicit flat coordinate indices (or `(atom, axis)` pairs) to displace; mutually exclusive with `sel`/`axes`
        :type which: Iterable | None
        :param sel: a selection of atoms to restrict the displacement to (used together with `axes`)
        :type sel: Iterable[int] | None
        :param axes: which Cartesian axes to displace along (used together with `sel`)
        :type axes: Iterable[int] | None
        :param use_internals: whether the displacement (and/or resulting coordinates) should be interpreted/handled in internal coordinates; `'reembed'`/`'convert'` select how the internal-coordinate result is converted back
        :type use_internals: bool | str
        :param coordinate_expansion: a Taylor-series coordinate-transformation expansion (e.g. a normal-mode Jacobian) to displace along instead of raw coordinates
        :type coordinate_expansion: list[np.ndarray] | None
        :param expansion_active_positions: restrict `coordinate_expansion`'s output to these positions before applying
        :type expansion_active_positions: Iterable[int] | None
        :param strip_embedding: whether the coordinate expansion/target coordinates have had the fixed embedding coordinates stripped
        :type strip_embedding: bool
        :param shift: whether `displacements` are relative shifts (`True`) added to the reference, or absolute target values (`False`)
        :type shift: bool
        :param coords: alternate reference coordinates to displace from, instead of `self.embedding.coords`/`internal_coordinates`
        :type coords: np.ndarray | None
        :return: the displaced coordinates
        :rtype: CoordinateSet
        :raises ValueError: if `use_internals` is requested but the embedding has no internal-coordinate spec, or if `displacements` doesn't match the reference coordinates' shape
        """
        ...

    def get_scan_coordinates(self, domains, internals=False, which=None, sel=None, axes=None, coordinate_expansion=None, strip_embedding=False, shift=True, return_displacements=False):
        """
        **LLM Docstring**

        Build a grid of displaced coordinates spanning the given `(start, stop, num)` domains for each scanned coordinate, via `get_displaced_coordinates` applied to a meshgrid of displacement values.

        :param domains: the `(start, stop, num)` ranges to scan over, one per scanned coordinate
        :type domains: Iterable[tuple]
        :param internals: whether the scan coordinates/displacements are internal coordinates
        :type internals: bool
        :param which: explicit coordinate indices to scan; forwarded to `get_displaced_coordinates`
        :type which: Iterable | None
        :param sel: atom selection to scan over; forwarded to `get_displaced_coordinates`
        :type sel: Iterable[int] | None
        :param axes: axes to scan over; forwarded to `get_displaced_coordinates`
        :type axes: Iterable[int] | None
        :param coordinate_expansion: a coordinate-transformation expansion to scan along instead of raw coordinates
        :type coordinate_expansion: list[np.ndarray] | None
        :param strip_embedding: whether to strip the fixed embedding coordinates
        :type strip_embedding: bool
        :param shift: whether the scan values are relative shifts or absolute targets
        :type shift: bool
        :param return_displacements: whether to also return the raw displacement mesh used
        :type return_displacements: bool
        :return: the scan coordinates, or `(displacement_mesh, scan_coordinates)` if `return_displacements` is set
        :rtype: np.ndarray | tuple
        """
        ...

    def get_nearest_displacement_atoms(self, points, sel=None, axes=None, weighting_function=None, return_distances=False):
        """
        **LLM Docstring**

        For each of a set of query points, find the atom (restricted to `sel`, along `axes`) whose (optionally mass-weighted) distance to the point is smallest.

        :param points: the query points to find nearest atoms for
        :type points: np.ndarray
        :param sel: a selection of atoms to restrict the search to; defaults to all atoms
        :type sel: Iterable[int] | None
        :param axes: which Cartesian axes to compute distances over; defaults to all three
        :type axes: Iterable[int] | None
        :param weighting_function: a function of the atomic masses used to weight distances (e.g. to prefer displacing lighter atoms); defaults to `np.sqrt`
        :type weighting_function: callable | None
        :param return_distances: whether to also return the (weighted) distances to the nearest atoms
        :type return_distances: bool
        :return: the nearest atom index for each point, or `(atom_idx, distances)` if `return_distances` is set
        :rtype: np.ndarray | tuple
        """
        ...

    def get_nearest_displacement_coordinates(self, points, sel=None, axes=None, weighting_function=None, modes_nearest=False, return_distances=False):
        """
        Displaces the _nearest_ atom (in a mass-weighted sense) to the given point
        This allows for the development of functions / potentials that are in the small-displacement limit
        where everything _except_ the nearest atom to the given Cartesian position is at its equilibrium value

        :param points:
        :param sel:
        :param axes:
        :param weighting_function:
        :return:
        """
        ...

    def get_nearest_scan_coordinates(self, domains, sel=None, axes=None):
        """
        **LLM Docstring**

        Build a grid of "nearest displacement" coordinates (see `get_nearest_displacement_coordinates`) spanning the given `(start, stop, num)` domains, via a meshgrid of query points.

        :param domains: the `(start, stop, num)` ranges to scan over
        :type domains: Iterable[tuple]
        :param sel: a selection of atoms to restrict the nearest-atom search to
        :type sel: Iterable[int] | None
        :param axes: which Cartesian axes the domains correspond to
        :type axes: Iterable[int] | None
        :return: the grid of nearest-displacement coordinates
        :rtype: np.ndarray
        """
        ...

class InternalHandlingMode(enum.Enum):
    ReembedCartesians = 'reembed'
    Convert = 'convert'
    Cartesians = 'cartesians'

    @classmethod
    def resolve(cls, mode):
        """
        **LLM Docstring**

        Normalize a user-supplied internal-coordinate handling flag into an `InternalHandlingMode` enum value: `True` maps to `Convert`, `False` maps to `Cartesians`, and any other value is passed through to the enum constructor (e.g. an already-valid mode string).

        :param mode: the raw handling-mode flag to resolve
        :type mode: bool | str | InternalHandlingMode
        :return: the resolved enum member
        :rtype: InternalHandlingMode
        """
        ...

class PropertyEvaluator(metaclass=abc.ABCMeta):

    def __init__(self, embedding=None, permutation=None, use_internals=True, strip_embedding=True, flatten_internals=True, reembed_cartesians=False, supports_internals=None, **defaults):
        """
        **LLM Docstring**

        Set up the shared configuration for evaluating a molecular property: which coordinate embedding to use, whether/how to work in internal coordinates (falling back to `'reembed'` if the concrete evaluator doesn't support internals directly), whether to strip fixed embedding coordinates, whether to flatten internal-coordinate arrays, whether to Eckart-reembed Cartesians, an atom permutation to apply, and any evaluator-specific default keyword options.

        :param embedding: the coordinate embedding to evaluate the property against
        :type embedding: MolecularEmbedding | None
        :param permutation: an atom (or flattened-coordinate) permutation to apply when embedding/unembedding coordinates and derivatives
        :type permutation: np.ndarray | None
        :param use_internals: whether to evaluate in internal coordinates; if the concrete evaluator subclass doesn't support internals directly, this is coerced to `'reembed'`
        :type use_internals: bool | str
        :param strip_embedding: whether to strip the fixed embedding coordinates when converting to/from internal coordinates
        :type strip_embedding: bool
        :param flatten_internals: whether to flatten internal-coordinate arrays to a single trailing axis
        :type flatten_internals: bool
        :param reembed_cartesians: whether to Eckart-reembed Cartesian coordinates before evaluating
        :type reembed_cartesians: bool
        :param supports_internals: whether this evaluator can evaluate directly in internal coordinates; defaults to the class attribute `self.supports_internals`
        :type supports_internals: bool | None
        :param defaults: extra default keyword options passed through to `evaluate_expansion`/`evaluate_term`
        :type defaults: dict
        :return: None
        :rtype: None
        """
        ...

    def use_internal_coordinate_handlers(self):
        """
        **LLM Docstring**

        Whether this evaluator should route coordinate handling through the internal-coordinate machinery: true when `use_internals` is enabled, an embedding is set, and that embedding actually has internal coordinates defined.

        :return: whether internal-coordinate handling applies
        :rtype: bool
        """
        ...

    def embed_coords(self, coords, embed_reembedded=True):
        """
        **LLM Docstring**

        Convert raw Cartesian coordinates into whatever representation this evaluator actually evaluates in: internal coordinates (optionally stripped of fixed embedding coordinates and flattened) if `use_internals` is enabled and the embedding has internals, or Eckart-reembedded Cartesians if `reembed_cartesians` is set; then applies any configured atom/coordinate permutation.

        :param coords: the raw Cartesian coordinates to embed
        :type coords: np.ndarray
        :param embed_reembedded: whether to still perform the internal-coordinate conversion even when `use_internals` is specifically `'reembed'` (which otherwise defers embedding to elsewhere)
        :type embed_reembedded: bool
        :return: the embedded (and possibly permuted) coordinates, in whatever representation the evaluator uses
        :rtype: np.ndarray
        :raises NotImplementedError: if permuting non-flattened internal coordinates is requested
        """
        ...

    def unembed_coords(self, coords, embed_reembedded=True):
        """
        **LLM Docstring**

        Undo `embed_coords`: convert the evaluator's working coordinate representation back to plain Cartesians (via `self.embedding.get_cartesians` if internal coordinates were used), and undo any applied atom/coordinate permutation.

        :param coords: the evaluator's working-representation coordinates to unembed
        :type coords: np.ndarray
        :param embed_reembedded: whether to perform the internal-to-Cartesian conversion even when `use_internals` is specifically `'reembed'`
        :type embed_reembedded: bool
        :return: the Cartesian coordinates
        :rtype: np.ndarray
        :raises NotImplementedError: if un-permuting non-flattened internal coordinates is requested
        """
        ...

    def embed_derivs(self, coords, derivs, embed_reembedded=True):
        """
        **LLM Docstring**

        Re-express a set of derivative tensors (computed with respect to the evaluator's working coordinates) in terms of Cartesian coordinates, by re-expanding through the Cartesians-by-internals Jacobian if internal coordinates are in use.

        :param coords: the coordinates the derivatives were computed at (in the evaluator's working representation)
        :type coords: np.ndarray
        :param derivs: the derivative tensors to re-express
        :type derivs: list[np.ndarray]
        :param embed_reembedded: whether to perform the re-expansion even when `use_internals` is specifically `'reembed'`
        :type embed_reembedded: bool
        :return: the derivative tensors, re-expressed in Cartesian coordinates if applicable, or unchanged
        :rtype: list[np.ndarray]
        """
        ...

    def unembed_derivs(self, base_coords, coords, derivs, embed_reembedded=True):
        """
        **LLM Docstring**

        Undo `embed_derivs` (and any atom/coordinate permutation) on a set of Cartesian-coordinate derivative tensors, converting them back into the evaluator's working coordinate representation (internal coordinates, or Eckart-unrotated Cartesians) as appropriate.

        :param base_coords: the reference coordinates (in the evaluator's working representation) the returned derivatives should be expressed relative to
        :type base_coords: np.ndarray
        :param coords: the Cartesian coordinates the input derivatives were computed at
        :type coords: np.ndarray
        :param derivs: the Cartesian-coordinate derivative tensors to convert
        :type derivs: list[np.ndarray]
        :param embed_reembedded: whether to perform the re-expansion even when `use_internals` is specifically `'reembed'`
        :type embed_reembedded: bool
        :return: the derivative tensors, converted to the evaluator's working representation
        :rtype: list[np.ndarray]
        """
        ...

    @abc.abstractmethod
    def evaluate_term(self, coords, order, **opts):
        """
        **LLM Docstring**

        Abstract hook for evaluating the property's value (or a specific analytic derivative order) at the given coordinates. Concrete evaluator subclasses must implement this.

        :param coords: the coordinates to evaluate at, in the evaluator's working representation
        :type coords: np.ndarray
        :param order: the analytic derivative order to evaluate
        :type order: int
        :param opts: extra evaluator-specific options
        :type opts: dict
        :return: the evaluated term(s)
        :rtype: object
        """
        ...

    @classmethod
    def prep_mol_opts(cls, mol, embedding=None, charge=None, multiplicity=None, **opts):
        """
        **LLM Docstring**

        Build the keyword-argument dict used to construct an evaluator from a molecule, pulling the embedding, charge, and spin/multiplicity off the molecule wherever not explicitly overridden.

        :param mol: the molecule to build evaluator options from
        :type mol: AbstractMolecule
        :param embedding: an explicit embedding to use instead of `mol.embedding`
        :type embedding: MolecularEmbedding | None
        :param charge: an explicit charge to use instead of `mol.charge`
        :type charge: int | None
        :param multiplicity: an explicit spin multiplicity to use instead of `mol.spin`
        :type multiplicity: object | None
        :param opts: additional options passed through unchanged
        :type opts: dict
        :return: the resolved options dict
        :rtype: dict
        """
        ...

    @classmethod
    @abc.abstractmethod
    def from_mol(cls, mol, **opts):
        """
        **LLM Docstring**

        Abstract constructor for building an evaluator instance from a molecule. Concrete evaluator subclasses must implement this.

        :param mol: the molecule to build the evaluator for
        :type mol: AbstractMolecule
        :param opts: evaluator-specific construction options
        :type opts: dict
        :return: the constructed evaluator
        :rtype: PropertyEvaluator
        """
        ...

    def get_fd_opts(self, **opts):
        """
        **LLM Docstring**

        Merge the class-level `fd_defaults` with any explicitly passed finite-difference options.

        :param opts: explicit per-call finite-difference option overrides
        :type opts: dict
        :return: the merged finite-difference options
        :rtype: dict
        """
        ...

    def finite_difference_derivs(self, coords, order, batched_orders=None, displacement_generator=None, coordinate_prep=None, index_filter=None, analytic_derivative_order=None, **opts):
        """
        **LLM Docstring**

        Compute derivative tensors of the property numerically via finite differences in the evaluator's working (typically Cartesian) coordinates, wrapping `evaluate_term` as the base function fed to a `FiniteDifferenceDerivative`.

        :param coords: the reference coordinates to differentiate around
        :type coords: np.ndarray
        :param order: the highest total derivative order needed
        :type order: int
        :param batched_orders: whether `evaluate_term` returns all derivative orders in one batched call rather than being called once per order
        :type batched_orders: bool | None
        :param displacement_generator: a callable to transform the raw finite-difference displacement points before evaluation (e.g. to apply a coordinate expansion)
        :type displacement_generator: callable | None
        :param coordinate_prep: a callable to preprocess the flattened center coordinates before finite-differencing
        :type coordinate_prep: callable | None
        :param index_filter: a callable/filter restricting which derivative-tensor index combinations are computed
        :type index_filter: callable | None
        :param analytic_derivative_order: the derivative order already available analytically (via `evaluate_term`), used as the finite-difference base point rather than order 0
        :type analytic_derivative_order: int | None
        :param opts: extra options, split between `FiniteDifferenceDerivative` construction options and options forwarded to `evaluate_term`
        :type opts: dict
        :return: the finite-difference-computed derivative tensors, from `analytic_derivative_order + 1` up to `order`
        :rtype: list[np.ndarray]
        :raises NotImplementedError: if `multi_expansion_order > 0` (multi-term expansion unsplitting isn't supported yet)
        """
        ...

    def internal_finite_difference_derivs(self, coords, order, batched_orders=None, displacement_generator=None, coordinate_prep=None, index_filter=None, analytic_derivative_order=None, **opts):
        """
        **LLM Docstring**

        Compute derivative tensors of the property numerically via finite differences directly in internal coordinates, converting each finite-difference-displaced internal-coordinate point back to Cartesians (and re-embedding any analytic derivatives) before calling `evaluate_term`, when `use_internals` is set to `'reembed'`.

        :param coords: the reference coordinates (internal or Cartesian, converted to internal coordinates if needed) to differentiate around
        :type coords: np.ndarray
        :param order: the highest total derivative order needed
        :type order: int
        :param batched_orders: whether `evaluate_term` returns all derivative orders in one batched call
        :type batched_orders: bool | None
        :param displacement_generator: a callable to transform the raw finite-difference displacement points before evaluation
        :type displacement_generator: callable | None
        :param coordinate_prep: a callable to preprocess the flattened center coordinates before finite-differencing
        :type coordinate_prep: callable | None
        :param index_filter: a callable/filter restricting which derivative-tensor index combinations are computed
        :type index_filter: callable | None
        :param analytic_derivative_order: the derivative order already available analytically
        :type analytic_derivative_order: int | None
        :param opts: extra options, split between `FiniteDifferenceDerivative` construction options and options forwarded to `evaluate_term`
        :type opts: dict
        :return: the finite-difference-computed internal-coordinate derivative tensors
        :rtype: list[np.ndarray]
        """
        ...
    'Real access pattern: PropertyEvaluator.<AttrName> (7 class attributes, e.g. PropertyEvaluator.supports_internals == False). Collapsed into a dict below purely for compactness -- do not index it as a dict in real code:'
    _MEMBERS = {'supports_internals': False, 'property_units': None, 'target_property_units': None, 'distance_units': 'Angstroms', 'batched_orders': False, 'analytic_derivative_order': 1, 'multi_expansion_order': 0}

    def evaluate_expansion(self, coords, order=0, analytic_derivative_order=None, batched_orders=None, logger=None, fd_displacement_generator=None, fd_coordinate_prep=None, fd_index_filter=None, fd_handler=None, **opts):
        """
        **LLM Docstring**

        Build the property's Taylor expansion up to the requested derivative order(s), combining analytic derivatives (via `evaluate_term`, for orders up to `analytic_derivative_order`) with finite-difference derivatives (via `fd_handler`, for any higher orders requested), then applying the evaluator's energy/distance unit conversions.

        :param coords: the coordinates to evaluate the expansion at
        :type coords: np.ndarray
        :param order: the derivative order(s) to compute, either a single integer (meaning `0..order`) or an explicit sorted list of orders
        :type order: int | list[int]
        :param analytic_derivative_order: the highest order to compute analytically via `evaluate_term`; defaults to `self.analytic_derivative_order`
        :type analytic_derivative_order: int | None
        :param batched_orders: whether `evaluate_term` computes all analytic orders in one batched call
        :type batched_orders: bool | None
        :param logger: accepted for interface consistency but not used in this method's body
        :type logger: Logger | None
        :param fd_displacement_generator: forwarded to `fd_handler` as `displacement_generator`
        :type fd_displacement_generator: callable | None
        :param fd_coordinate_prep: forwarded to `fd_handler` as `coordinate_prep`
        :type fd_coordinate_prep: callable | None
        :param fd_index_filter: forwarded to `fd_handler` as `index_filter`
        :type fd_index_filter: callable | None
        :param fd_handler: the finite-difference routine to use for orders beyond `analytic_derivative_order`; defaults to `self.finite_difference_derivs`
        :type fd_handler: callable | None
        :param opts: extra options, merged with `self.defaults` and split between finite-difference construction options and options forwarded to `evaluate_term`
        :type opts: dict
        :return: the expansion terms, one per requested order (in units converted per `property_units`/`distance_units`)
        :rtype: list
        """
        ...

    def evaluate(self, coords, order=0, logger=None, fd_handler=None, analytic_derivative_order=None, unembed_derivatives=True, **opts):
        """
        **LLM Docstring**

        Top-level property-evaluation entry point: embeds the input coordinates into the evaluator's working representation, computes the Taylor expansion (temporarily disabling internal double-handling/permutation/re-embedding so `evaluate_expansion` operates on already-embedded coordinates), re-expands analytic-order derivatives back through the Cartesians-by-internals Jacobian if working in internal coordinates, and finally un-embeds the resulting derivatives back to plain Cartesian coordinates (unless `unembed_derivatives` is `False`).

        :param coords: the (Cartesian) coordinates to evaluate the property at
        :type coords: np.ndarray
        :param order: the derivative order(s) to compute
        :type order: int | list[int]
        :param logger: accepted for interface consistency but not used directly in this method's body
        :type logger: Logger | None
        :param fd_handler: an explicit finite-difference handler to use; defaults to `internal_finite_difference_derivs` if working in internal coordinates
        :type fd_handler: callable | None
        :param analytic_derivative_order: the highest order to compute analytically
        :type analytic_derivative_order: int | None
        :param unembed_derivatives: whether to convert the resulting derivatives back to plain Cartesian coordinates
        :type unembed_derivatives: bool
        :param opts: extra options forwarded to `evaluate_expansion`
        :type opts: dict
        :return: the property expansion (value plus derivatives), in Cartesian coordinates unless `unembed_derivatives` is `False`
        :rtype: list
        :raises NotImplementedError: if `multi_expansion_order > 0` and internal-coordinate re-expansion is needed
        """
        ...

    @staticmethod
    def _diag_inds(pos):
        """
        **LLM Docstring**

        Filter a list of index tuples down to just the "diagonal" ones, where every index in the tuple is the same.

        :param pos: the index tuples to filter
        :type pos: list[tuple]
        :return: the diagonal index tuples
        :rtype: list[tuple]
        """
        ...

    def partial_force_field(self, coords: np.ndarray, modes, mode_distance_units='BohrRadius', analytic_derivative_order=None, order=4, index_filter=None, fd_handler=None, mesh_spacing=0.5, **opts):
        """
        **LLM Docstring**

        Compute a partial force-field-style expansion of the property restricted to displacements along a given set of normal modes (rather than full Cartesian/internal derivatives), by generating finite-difference displacement points along the mode-coordinate directions and only requesting the "diagonal" derivative terms by default.

        :param coords: the reference coordinates to expand around
        :type coords: np.ndarray
        :param modes: the normal modes to expand the property in
        :type modes: object
        :param mode_distance_units: the units the mode displacement coordinates are given in, converted to `self.distance_units`
        :type mode_distance_units: str
        :param analytic_derivative_order: the highest order to compute analytically
        :type analytic_derivative_order: int | None
        :param order: the highest derivative order to compute
        :type order: int
        :param index_filter: which derivative-tensor index combinations to compute; defaults to `_diag_inds` (only same-mode-repeated terms)
        :type index_filter: callable | None
        :param fd_handler: an explicit finite-difference handler to use
        :type fd_handler: callable | None
        :param mesh_spacing: the finite-difference step size (normal modes typically need larger steps than Cartesian/internal coordinates)
        :type mesh_spacing: float
        :param opts: extra options forwarded to `evaluate`
        :type opts: dict
        :return: the mode-basis partial force-field expansion, with derivative orders beyond `analytic_derivative_order` rescaled back from mode units
        :rtype: list
        """
        ...
    evaluator_registry: dict

    @classmethod
    def register(cls, name, method=None):
        """
        **LLM Docstring**

        Register a named evaluator-construction method (or function) in the class's `evaluator_registry`, usable directly or as a decorator.

        :param name: the registry key to register under, or (if `method` is `None` and `name` has a `.name` attribute) the method/callable itself
        :type name: str | callable
        :param method: the evaluator-construction method/callable to register; if `None`, `register` returns a decorator instead
        :type method: callable | None
        :return: the registered method (if called directly), or a decorator function (if used as `@register(name)`)
        :rtype: callable
        """
        ...

    @classmethod
    def get_evaluators(cls):
        """
        **LLM Docstring**

        The set of built-in named evaluator constructors for this evaluator type. Returns an empty dict on the base class; concrete subclasses override this to expose their built-in evaluator options.

        :return: mapping from evaluator name to constructor
        :rtype: dict
        """
        ...

    @classmethod
    def get_evaluator_map(cls):
        """
        **LLM Docstring**

        The full set of named evaluator constructors available for this evaluator type: the built-ins from `get_evaluators`, overlaid with any user-registered ones from `evaluator_registry`.

        :return: the combined evaluator-name-to-constructor mapping
        :rtype: dict
        """
        ...

    @classmethod
    def get_evaluators_by_attributes(cls):
        """
        **LLM Docstring**

        Optional mapping used to resolve an evaluator by a set of attribute-based hints rather than by exact name, for use by the dispatch machinery. Returns `None` on the base class (no attribute-based dispatch); concrete subclasses may override this.

        :return: the attribute-based dispatch mapping, or `None`
        :rtype: dict | None
        """
        ...
    _profile_dispatch = dev.uninitialized
    default_evaluator_type = 'expansion'

    @classmethod
    def profile_generator_dispatch(cls):
        """
        **LLM Docstring**

        Lazily build (and cache) the `OptionsMethodDispatch` object used to resolve a named/attribute-described evaluator specification into a concrete evaluator constructor and options.

        :return: the cached dispatch object
        :rtype: dev.OptionsMethodDispatch
        """
        ...

    @classmethod
    @abc.abstractmethod
    def get_default_function_evaluator_type(cls):
        """
        **LLM Docstring**

        Abstract hook returning the `PropertyFunctionEvaluator`-family class to use when a bare callable is supplied as an evaluator specification. Concrete evaluator subclasses must implement this.

        :return: the default function-evaluator class
        :rtype: type
        """
        ...

    @classmethod
    def _resolve_torch_device(cls, device=None):
        """
        **LLM Docstring**

        Resolve which PyTorch device to run a model on, auto-detecting CUDA availability if not explicitly given.

        :param device: an explicit device string; if `None`, resolves to `'cuda'` if available, otherwise `'cpu'`
        :type device: str | None
        :return: the resolved device string
        :rtype: str
        """
        ...

    @classmethod
    def handle_specialization(cls, tag):
        """
        **LLM Docstring**

        Hook for resolving an evaluator "specialization" tag (the part after a `:` in an evaluator name string, e.g. `'method:tag'`) into extra options. Not supported on the base class; concrete evaluator subclasses may override this to support tags.

        :param tag: the specialization tag to handle
        :type tag: str
        :return: never returns on the base class
        :rtype: dict
        :raises NotImplementedError: always, on classes that don't override this
        """
        ...

    @classmethod
    def resolve_evaluator(cls, name):
        """
        **LLM Docstring**

        Resolve an evaluator specification -- a bare callable, an object already exposing `evaluate`, or a (possibly `'name:tag'`-qualified) name/dict/tuple understood by the evaluator-dispatch machinery -- into a concrete evaluator class/instance plus its construction options.

        :param name: the evaluator specification to resolve
        :type name: object | str | dict | tuple | list | None
        :return: `(evaluator, opts)` -- the resolved evaluator class/instance and its construction options
        :rtype: tuple
        :raises ValueError: if `name` is a non-callable, non-string/dict/tuple/list object that also doesn't expose `evaluate`
        """
        ...

    class quiet_mode:

        def __init__(self, quiet=True):
            """
            **LLM Docstring**

            Set up a context manager that optionally silences stdout for the duration of its `with` block.

            :param quiet: whether to actually silence stdout; if `False`, the context manager is a no-op
            :type quiet: bool
            :return: None
            :rtype: None
            """
            ...

        def __enter__(self):
            """
            **LLM Docstring**

            Redirect `sys.stdout` to `os.devnull` if `self.quiet` is set.

            :return: None
            :rtype: None
            """
            ...

        def __exit__(self, exc_type, exc_val, exc_tb):
            """
            **LLM Docstring**

            Restore `sys.stdout` (and close the devnull stream) if `self.quiet` is set.

            :param exc_type: the exception type, if an exception occurred in the `with` block
            :type exc_type: type | None
            :param exc_val: the exception instance, if any
            :type exc_val: BaseException | None
            :param exc_tb: the exception traceback, if any
            :type exc_tb: object | None
            :return: None
            :rtype: None
            """
            ...

class PropertyFunctionEvaluator(PropertyEvaluator):

    def __init__(self, potential_function, property_units=None, distance_units='Angstroms', batched_orders=False, analytic_derivative_order=np.inf, **defaults):
        """
        **LLM Docstring**

        Wrap a bare Python function as a `PropertyEvaluator`, storing its unit conventions and derivative-computation characteristics.

        :param potential_function: the function to wrap, expected to accept `(coords, order=..., **opts)`
        :type potential_function: callable
        :param property_units: the energy/property units the function's output is given in
        :type property_units: str | None
        :param distance_units: the distance units the function expects its input coordinates in
        :type distance_units: str
        :param batched_orders: whether the function returns all requested derivative orders in a single batched call
        :type batched_orders: bool
        :param analytic_derivative_order: the highest order the function can compute analytically; defaults to `inf` (the function is assumed to support any order directly)
        :type analytic_derivative_order: float | int
        :param defaults: extra default options forwarded to `PropertyEvaluator.__init__` (note: dispatched via `super(type(self).__bases__[0], self)`, skipping the immediate `PropertyEvaluator.__init__` override chain)
        :type defaults: dict
        :return: None
        :rtype: None
        """
        ...

    def evaluate_term(self, coords, order, **opts):
        """
        **LLM Docstring**

        Evaluate the wrapped function directly at the given coordinates/order.

        :param coords: the coordinates to evaluate at
        :type coords: np.ndarray
        :param order: the derivative order to request from the function
        :type order: int
        :param opts: extra options forwarded to the function
        :type opts: dict
        :return: the function's return value
        :rtype: object
        """
        ...
    default_property_function = None

    @classmethod
    def bind_default(cls, potential):
        """
        **LLM Docstring**

        Register a default property function to be used the next time `initialize_from_mol` is called without an explicit `property_function`.

        :param potential: the function to use as the default
        :type potential: callable
        :return: None
        :rtype: None
        """
        ...

    @classmethod
    def get_property_function(cls, prop_func, mol, **opts):
        """
        **LLM Docstring**

        Resolve the property function to wrap. On this base class, simply returns `prop_func` unchanged; subclasses may override this to do additional resolution (e.g. binding molecule-specific parameters).

        :param prop_func: the candidate property function
        :type prop_func: callable
        :param mol: the molecule the evaluator is being built for
        :type mol: AbstractMolecule
        :param opts: extra options, unused on the base class
        :type opts: dict
        :return: the resolved property function
        :rtype: callable
        """
        ...

    @classmethod
    def initializer(cls, mol, **opts):
        """
        **LLM Docstring**

        Build a bare (function-less) instance of this evaluator type; used as a lightweight constructor hook by subclasses/dispatch code that will bind a function afterward.

        :param mol: the molecule the evaluator is being built for (unused directly in this method's body)
        :type mol: AbstractMolecule
        :param opts: extra options forwarded to the constructor
        :type opts: dict
        :return: the constructed evaluator, with `potential_function=None`
        :rtype: PropertyFunctionEvaluator
        """
        ...

    @classmethod
    def initialize_from_mol(root, cls, mol, property_function=None, *, embedding=None, use_internals=None, reembed_cartesians=None, property_units=None, distance_units=None, batched_orders=None, analytic_derivative_order=None, **opts):
        """
        **LLM Docstring**

        Build a `PropertyFunctionEvaluator` for a molecule, resolving the property function to wrap (falling back to a previously `bind_default`-registered default function if none is given) and filling in construction options (embedding, unit conventions, derivative-order handling) from the molecule wherever not explicitly overridden.

        :param root: the class this method is unbound from (Python descriptor mechanics artifact; effectively unused directly)
        :type root: type
        :param cls: the evaluator class to instantiate
        :type cls: type
        :param mol: the molecule to build the evaluator for
        :type mol: AbstractMolecule
        :param property_function: the function to wrap; falls back to `cls.default_property_function` if not given
        :type property_function: callable | None
        :param embedding: an explicit embedding to use instead of `mol.embedding`
        :type embedding: MolecularEmbedding | None
        :param use_internals: whether to evaluate in internal coordinates
        :type use_internals: bool | None
        :param reembed_cartesians: whether to Eckart-reembed Cartesian coordinates
        :type reembed_cartesians: bool | None
        :param property_units: the units the function's output is given in
        :type property_units: str | None
        :param distance_units: the units the function expects its input in
        :type distance_units: str | None
        :param batched_orders: whether the function returns all derivative orders in one batched call
        :type batched_orders: bool | None
        :param analytic_derivative_order: the highest order the function computes analytically
        :type analytic_derivative_order: int | None
        :param opts: extra options forwarded to `get_property_function` and the constructor
        :type opts: dict
        :return: the constructed evaluator
        :rtype: PropertyFunctionEvaluator
        :raises ValueError: if no property function could be resolved
        """
        ...

class EnergyEvaluator(PropertyEvaluator):
    target_property_units = 'Hartrees'
    evaluator_registry = {}

    @classmethod
    def get_evaluators_by_attributes(cls):
        """
        **LLM Docstring**

        Attribute-based dispatch mapping for energy evaluators: recognizes a `potential_function` keyword as identifying a `PotentialFunctionEnergyEvaluator`.

        :return: the attribute-tuple-to-evaluator-class mapping
        :rtype: dict
        """
        ...

    @classmethod
    def get_default_function_evaluator_type(cls):
        """
        **LLM Docstring**

        The evaluator class used when a bare callable is supplied as an energy-evaluator specification.

        :return: `PotentialFunctionEnergyEvaluator`
        :rtype: type
        """
        ...

    def _modify_ase_calc(self, calc, gradient_modification_function=None, gradient_modification_mode=None, orthogonal_projection_generator=None, **opts):
        """
        **LLM Docstring**

        Wrap an ASE calculator's `calculate` method so that, if a gradient-modification function is supplied, the computed forces are intercepted, converted to a flat gradient, passed through `_modify_gradient`, and written back into the calculator's results in ASE's expected shape.

        :param calc: the ASE calculator to modify in place
        :type calc: object
        :param gradient_modification_function: the function used to modify the raw gradient/forces; if `None`, `calc` is returned unmodified
        :type gradient_modification_function: callable | None
        :param gradient_modification_mode: the modification mode forwarded to `_modify_gradient` (e.g. `'shift'`)
        :type gradient_modification_mode: str | None
        :param orthogonal_projection_generator: an optional per-geometry projector generator forwarded to `_modify_gradient`
        :type orthogonal_projection_generator: callable | None
        :param opts: extra options forwarded to `_modify_gradient`
        :type opts: dict
        :return: the (possibly modified) calculator
        :rtype: object
        """
        ...

    def to_ase(self, gradient_modification_function=None, gradient_modification_mode=None, convert_modification_distances=None, convert_modification_energies=None, orthogonal_projection_generator=None, **etc):
        """
        **LLM Docstring**

        Wrap this evaluator as an ASE `Calculator` object (via `ASECalculator`), optionally intercepting/modifying its gradient via `_modify_ase_calc`.

        :param gradient_modification_function: an optional function to modify the computed gradient/forces
        :type gradient_modification_function: callable | None
        :param gradient_modification_mode: the modification mode to use
        :type gradient_modification_mode: str | None
        :param convert_modification_distances: whether to convert distances to Bohr before calling the modification function
        :type convert_modification_distances: bool | None
        :param convert_modification_energies: whether to convert energies to Hartrees before calling the modification function
        :type convert_modification_energies: bool | None
        :param orthogonal_projection_generator: an optional per-geometry projector generator applied to the modified gradient
        :type orthogonal_projection_generator: callable | None
        :param etc: extra options forwarded to `ASECalculator`
        :type etc: dict
        :return: the constructed (and possibly gradient-modified) ASE calculator
        :rtype: object
        """
        ...

    def _modify_pysis_calc(self, calc, gradient_modification_function=None, gradient_modification_mode=None, convert_modification_distances=None, convert_modification_energies=None, **opts):
        """
        **LLM Docstring**

        Wrap a pysisyphus calculator's `get_forces` (and `get_hessian`, if present) so that, if a gradient-modification function is supplied, the computed forces/Hessian are intercepted, passed through `_modify_gradient`, and written back.

        :param calc: the pysisyphus calculator to modify in place
        :type calc: object
        :param gradient_modification_function: the function used to modify the raw forces; if `None`, `calc` is returned unmodified
        :type gradient_modification_function: callable | None
        :param gradient_modification_mode: the modification mode forwarded to `_modify_gradient`
        :type gradient_modification_mode: str | None
        :param convert_modification_distances: whether to convert distances to Bohr before calling the modification function
        :type convert_modification_distances: bool | None
        :param convert_modification_energies: whether to convert energies to Hartrees before calling the modification function
        :type convert_modification_energies: bool | None
        :param opts: extra options forwarded to `_modify_gradient`
        :type opts: dict
        :return: the (possibly modified) calculator
        :rtype: object
        """
        ...

    def to_pysis(self, gradient_modification_function=None, gradient_modification_mode=None, convert_modification_distances=None, convert_modification_energies=None, orthogonal_projection_generator=None, **etc):
        """
        **LLM Docstring**

        Wrap this evaluator as a pysisyphus calculator (via `PysisCalculator`), optionally intercepting/modifying its gradient via `_modify_pysis_calc`.

        :param gradient_modification_function: an optional function to modify the computed forces
        :type gradient_modification_function: callable | None
        :param gradient_modification_mode: the modification mode to use
        :type gradient_modification_mode: str | None
        :param convert_modification_distances: whether to convert distances to Bohr before calling the modification function
        :type convert_modification_distances: bool | None
        :param convert_modification_energies: whether to convert energies to Hartrees before calling the modification function
        :type convert_modification_energies: bool | None
        :param orthogonal_projection_generator: forwarded to `_modify_pysis_calc`
        :type orthogonal_projection_generator: callable | None
        :param etc: extra options forwarded to `PysisCalculator`
        :type etc: dict
        :return: the constructed (and possibly gradient-modified) pysisyphus calculator
        :rtype: object
        """
        ...

    def minimizer_function_by_order(self, order, allow_fd=False, modifier=None, **opts):
        """
        **LLM Docstring**

        Build a callable of the `scipy.optimize`-style `(coords, *extra) -> value` form for a specific derivative order (energy, gradient, or Hessian), using analytic evaluation when available and (optionally) finite differences otherwise, converting the result to Hartrees.

        :param order: the derivative order the function should evaluate (`0` energy, `1` gradient, `2` Hessian)
        :type order: int
        :param allow_fd: whether to fall back to `finite_difference_derivs` when the analytic derivative order isn't high enough
        :type allow_fd: bool
        :param modifier: an optional function to post-process the returned value
        :type modifier: callable | None
        :param opts: extra options forwarded to `evaluate_term`/`finite_difference_derivs`
        :type opts: dict
        :return: the constructed minimizer-compatible function, or `None` if the requested order isn't available and `allow_fd` is `False`
        :rtype: callable | None
        """
        ...

    def minimizer_func(self, **opts):
        """
        **LLM Docstring**

        Build a `scipy.optimize`-compatible energy-evaluation function, via `minimizer_function_by_order(0, ...)`.

        :param opts: extra options forwarded to `minimizer_function_by_order`
        :type opts: dict
        :return: the constructed energy function
        :rtype: callable
        """
        ...

    def minimizer_jacobian(self, **opts):
        """
        **LLM Docstring**

        Build a `scipy.optimize`-compatible gradient-evaluation function (falling back to finite differences if needed), via `minimizer_function_by_order(1, allow_fd=True, ...)`.

        :param opts: extra options forwarded to `minimizer_function_by_order`
        :type opts: dict
        :return: the constructed gradient function
        :rtype: callable
        """
        ...

    def minimizer_hessian(self, **opts):
        """
        **LLM Docstring**

        Build a `scipy.optimize`-compatible Hessian-evaluation function, via `minimizer_function_by_order(2, ...)`.

        :param opts: extra options forwarded to `minimizer_function_by_order`
        :type opts: dict
        :return: the constructed Hessian function
        :rtype: callable
        """
        ...

    def minimizer_internal_function_by_order(self, order, allow_fd=False, modifier=None, **opts):
        """
        **LLM Docstring**

        Like `minimizer_function_by_order`, but for evaluating in internal coordinates: when re-embedding, converts the working-representation coordinates back to Cartesians before calling `evaluate_term`, re-embeds any lower-order derivatives (with the appropriate Bohr-to-`distance_units` scaling), and converts the final result to Hartrees.

        :param order: the derivative order the function should evaluate
        :type order: int
        :param allow_fd: whether to fall back to `internal_finite_difference_derivs` when the analytic order isn't high enough
        :type allow_fd: bool
        :param modifier: an optional function to post-process the returned value
        :type modifier: callable | None
        :param opts: extra options forwarded to `evaluate_term`/`internal_finite_difference_derivs`
        :type opts: dict
        :return: the constructed minimizer-compatible function, or `None` if unavailable and `allow_fd` is `False`
        :rtype: callable | None
        """
        ...

    def get_internal_coordinate_indices(self, coord_spec):
        """
        **LLM Docstring**

        Resolve which internal-coordinate indices (and the total internal-coordinate dimension) correspond to a given coordinate specification, supporting both generic-internal (`'specs'`) and Z-matrix-based internal-coordinate systems.

        :param coord_spec: the coordinate(s) to look up, in whatever spec form the underlying internal-coordinate system expects
        :type coord_spec: object
        :return: `(ndim, inds)` -- the total number of internal coordinates and the resolved indices for `coord_spec`
        :rtype: tuple[int, list[int]]
        :raises ValueError: if the embedding's internal-coordinate spec isn't recognized
        """
        ...

    def get_internal_coordinate_projector(self, coord_spec, mask=None):
        """
        **LLM Docstring**

        Build a constraint function returning the orthogonal-projection matrix that removes the directions corresponding to a given set of internal coordinates (accounting for any redundant-coordinate transformation), for use in constrained optimization.

        :param coord_spec: the internal coordinate(s) to project out
        :type coord_spec: object
        :param mask: accepted for interface consistency but not used in this method's body
        :type mask: callable | None
        :return: a function `(coords) -> projector` broadcasting the projection matrix to the batch shape of `coords`
        :rtype: callable
        :raises ValueError: if the embedding's internal-coordinate spec isn't recognized
        """
        ...

    def get_internal_coordinate_constraints(self, constraints):
        """
        **LLM Docstring**

        Resolve a constraints specification for internal coordinates into a projector-generating constraint function: a dict of per-coordinate mask functions, a bare list of coordinates, a raw projection-matrix-like array, or an already-resolved constraint object passed through unchanged.

        :param constraints: the constraints specification to resolve
        :type constraints: dict | list | np.ndarray | object
        :return: the resolved constraint function/matrix
        :rtype: callable | np.ndarray | object
        """
        ...
    internal_optimizer_defaults = dict(max_displacement=0.01)

    @classmethod
    def get_optimizer_options(self):
        """
        **LLM Docstring**

        The full set of recognized keyword-argument names for `optimize`/`optimize_iterative`/`relaxed_scan`, combining the optimizer defaults, a few extra named options, and the `FiniteDifferenceDerivative` option names -- used to split a caller's `**opts` between optimizer-specific and evaluator-specific options.

        :return: the tuple of recognized option names
        :rtype: tuple[str]
        """
        ...

    def get_coordinate_projector(self, coord_spec, region_constraints=None, mask=None):
        """
        **LLM Docstring**

        Build a constraint function returning the orthogonal-projection matrix that removes the directions spanned by a general (Cartesian-space) coordinate specification, optionally restricting the projection to a specific region (via `region_constraints`) or an arbitrary mask function.

        :param coord_spec: the coordinate specification (e.g. bond/angle/dihedral definitions) whose directions should be projected out
        :type coord_spec: object
        :param region_constraints: per-coordinate `(min, max)` bounds restricting where the projection is applied
        :type region_constraints: Iterable[tuple] | None
        :param mask: an explicit mask function selecting which basis vectors to zero out, in addition to any region mask
        :type mask: callable | None
        :return: a function `(coords) -> projector` computing the (region-restricted) orthogonal projection matrix at the given coordinates
        :rtype: callable
        """
        ...

    def get_coordinate_constraints(self, constraints, region_constraints=None):
        """
        **LLM Docstring**

        Resolve a constraints specification for general (Cartesian-space) coordinates into a projector-generating constraint function: a dict of per-coordinate mask functions, a bare coordinate-list specification, a raw projection-matrix-like array, or an already-resolved constraint object passed through unchanged.

        :param constraints: the constraints specification to resolve
        :type constraints: dict | list | np.ndarray | object
        :param region_constraints: per-coordinate bounds forwarded to `get_coordinate_projector` when `constraints` is a coordinate-list spec
        :type region_constraints: Iterable[tuple] | None
        :return: the resolved constraint function/matrix
        :rtype: callable | np.ndarray | object
        """
        ...

    def _modify_gradient(self, gradient_function, modification_function, modification_mode='shift', convert_modification_distances=True, convert_modification_energies=True, orthogonal_projector=None, orthogonal_projection_generator=None, use_forces=False, coord_prep=None, grad_prep=None, grad_post=None):
        """
        **LLM Docstring**

        Build a wrapped version of a raw gradient-computing function that applies a user-supplied modification function to the gradient/forces (in `'shift'` mode, adding the modification's output to the raw gradient with appropriate unit conversions; otherwise, replacing the gradient with the modification function's output directly), then optionally projects the result through an orthogonal projector/projector-generator, converting units as needed on the way in and out.

        :param gradient_function: the underlying function that computes the raw gradient/forces
        :type gradient_function: callable
        :param modification_function: the function used to modify the gradient, called as `(coords_in_bohr, gradient)` (or with a sign flip if `use_forces`)
        :type modification_function: callable
        :param modification_mode: `'shift'` to add the modification's output to the raw gradient, or any other value to replace the gradient with it directly
        :type modification_mode: str
        :param convert_modification_distances: whether to convert distances to Bohr before calling the modification function
        :type convert_modification_distances: bool
        :param convert_modification_energies: whether to convert energies to Hartrees before calling the modification function (`'shift'` mode only)
        :type convert_modification_energies: bool
        :param orthogonal_projector: a fixed projection matrix to apply to the modified gradient
        :type orthogonal_projector: np.ndarray | None
        :param orthogonal_projection_generator: a per-geometry projector-generating function, combined with `orthogonal_projector` if both are given
        :type orthogonal_projection_generator: callable | None
        :param use_forces: whether the wrapped function operates on forces (sign-flipped gradients) rather than raw gradients
        :type use_forces: bool
        :param coord_prep: a function to extract/prepare the coordinates from the wrapped function's raw arguments
        :type coord_prep: callable | None
        :param grad_prep: a function to extract the raw gradient (and any extra state) from the wrapped function's result, or signal a short-circuit by returning `(None, state)`
        :type grad_prep: callable | None
        :param grad_post: a function to reinsert the modified gradient back into the original result structure
        :type grad_post: callable | None
        :return: the wrapped gradient-computing function
        :rtype: callable
        :raises ValueError: if `modification_mode` is `None`
        """
        ...
    scipy_no_hessian_methods = {'cg', 'bfgs'}
    scipy_no_grad_methods = {'nelder-mead'}
    use_scipy_linesearch = False

    def optimize_iterative(self, coords, coordinate_constraints=None, gradient_modification_function=None, gradient_modification_mode='shift', convert_modification_distances=True, convert_modification_energies=True, return_trajectory=False, initialization_function=None, line_search_step=None, **opts):
        """
        **LLM Docstring**

        Core geometry-optimization driver supporting several backends (a custom iterative step-minimization scheme by default, or delegating to `scipy.optimize.minimize`, ASE's optimizers, or pysisyphus, depending on `mode`), building the energy/gradient/Hessian callables (falling back to the evaluator's own `minimizer_func`/`minimizer_jacobian`/`minimizer_hessian` if not supplied), applying any coordinate constraints as an orthogonal projector, and optionally recording/returning the optimization trajectory.

        :param coords: the initial coordinates to optimize from
        :type coords: np.ndarray
        :param coordinate_constraints: coordinates to hold fixed/project out during optimization
        :type coordinate_constraints: object | None
        :param gradient_modification_function: an optional function to modify the gradient/forces during optimization (e.g. to bias the search)
        :type gradient_modification_function: callable | None
        :param gradient_modification_mode: the modification mode forwarded to `_modify_gradient`
        :type gradient_modification_mode: str
        :param convert_modification_distances: whether to convert distances to Bohr before calling the modification function
        :type convert_modification_distances: bool
        :param convert_modification_energies: whether to convert energies to Hartrees before calling the modification function
        :type convert_modification_energies: bool
        :param return_trajectory: whether to record and return the optimization trajectory
        :type return_trajectory: bool
        :param initialization_function: an optional function applied to the initial coordinates (in Bohr) before optimization begins
        :type initialization_function: callable | None
        :param line_search_step: a fixed step size to force during (non-`'scipy'`-line-search) line searches, bypassing the Wolfe-condition search
        :type line_search_step: float | None
        :param opts: the full set of optimizer options (method, tolerances, damping, region/coordinate constraints, backend-specific settings, etc.), merged with `self.optimizer_defaults`
        :type opts: dict
        :return: `(converged, coords_or_(coords,trajectory), info)` -- whether the optimization converged, the optimized (and optionally trajectory) coordinates, and backend-specific optimization metadata
        :rtype: tuple
        """
        ...

    def optimize(self, coords, use_internals=None, func=None, jacobian=None, hessian=None, logger=None, coordinate_constraints=None, orthogonal_projection_generator=None, gradient_modification_function=None, initialization_function=None, return_trajectory=False, mode=None, **opts):
        """
        **LLM Docstring**

        Geometry-optimize a structure, automatically routing through internal coordinates (building the appropriate internal-coordinate minimizer functions and constraints, then converting the result back to Cartesians) when `use_internals` is enabled and internal coordinates are available and the backend `mode` supports it, otherwise delegating directly to `optimize_iterative`.

        :param coords: the initial coordinates to optimize from
        :type coords: np.ndarray
        :param use_internals: whether to optimize in internal coordinates; defaults to `self.use_internals`
        :type use_internals: bool | None
        :param func: an explicit energy function to use instead of the default minimizer function
        :type func: callable | None
        :param jacobian: an explicit gradient function to use
        :type jacobian: callable | None
        :param hessian: an explicit Hessian function to use
        :type hessian: callable | None
        :param logger: logger to report optimization progress to
        :type logger: Logger | str | None
        :param coordinate_constraints: coordinates to hold fixed/project out during optimization
        :type coordinate_constraints: object | None
        :param orthogonal_projection_generator: an explicit per-geometry constraint projector, bypassing automatic construction from `coordinate_constraints`
        :type orthogonal_projection_generator: callable | None
        :param gradient_modification_function: an optional function to modify the gradient during optimization
        :type gradient_modification_function: callable | None
        :param initialization_function: an optional function applied to the initial coordinates before optimization
        :type initialization_function: callable | None
        :param return_trajectory: whether to record and return the optimization trajectory
        :type return_trajectory: bool
        :param mode: the optimization backend to use (`'scipy'`, `'ase'`, `'pysis'`, or `None` for the default iterative scheme); internal-coordinate optimization is skipped for `'ase'`/`'pysis'`
        :type mode: str | None
        :param opts: extra options forwarded to `optimize_iterative`
        :type opts: dict
        :return: `(converged, coords_or_(coords,trajectory), info)`
        :rtype: tuple
        """
        ...

    @classmethod
    def get_relaxed_scan_options(cls):
        """
        **LLM Docstring**

        The full set of recognized keyword-argument names for `relaxed_scan`, extending `get_optimizer_options` with the scan-specific option names.

        :return: the tuple of recognized option names
        :rtype: tuple[str]
        """
        ...

    def relaxed_scan(self, coords, displacement_specs, displacement_coords, coordinate_constraints=None, absolute_mesh=False, orthogonal_projection_generator=None, coordinate_expansion=None, expansion_active_positions=None, include_displacement_constraints=True, shift=True, adjust_displacements=True, split_scan_mesh=True, **optimization_settings):
        """
        **LLM Docstring**

        Perform a constrained ("relaxed") potential-energy scan: for each point on a mesh of displacement values along the specified coordinate(s), build a displaced starting structure (via a coordinate-index-based, coordinate-expansion-based, or user-supplied displacement-generator scheme) and re-optimize it subject to constraints holding the scanned coordinates fixed, splitting the scan mesh at the point nearest zero displacement so both directions can be scanned outward (zig-zag) from the reference geometry if `split_scan_mesh` is set.

        :param coords: the reference (starting) coordinates for the scan
        :type coords: np.ndarray
        :param displacement_specs: the `(start, stop, num)` range(s) (or explicit displacement values, if `absolute_mesh`) to scan each displaced coordinate over
        :type displacement_specs: Iterable
        :param displacement_coords: which coordinate(s) to displace along -- an internal-coordinate index/spec, a dict mapping coordinates to base displacement weights (defining a `coordinate_expansion` direction), or a callable/`(callable, callable)` pair generating custom displacement directions
        :type displacement_coords: object | dict | callable | tuple
        :param coordinate_constraints: additional coordinates to hold fixed during each relaxation, beyond the scanned ones
        :type coordinate_constraints: object | None
        :param absolute_mesh: whether `displacement_specs` gives explicit displacement values rather than `(start, stop, num)` ranges
        :type absolute_mesh: bool
        :param orthogonal_projection_generator: an explicit constraint projector generator, bypassing automatic construction
        :type orthogonal_projection_generator: callable | None
        :param coordinate_expansion: an explicit coordinate-transformation expansion direction to displace along
        :type coordinate_expansion: list[np.ndarray] | None
        :param expansion_active_positions: restrict `coordinate_expansion`'s output to these positions
        :type expansion_active_positions: Iterable[int] | None
        :param include_displacement_constraints: whether the scanned coordinates should also be added to `coordinate_constraints` (held fixed during relaxation)
        :type include_displacement_constraints: bool
        :param shift: whether displacement values are relative shifts or absolute targets
        :type shift: bool
        :param adjust_displacements: whether to convert the mesh of absolute displacement values into incremental steps between successive optimizations
        :type adjust_displacements: bool
        :param split_scan_mesh: whether to split the mesh at the point nearest zero and scan outward in both directions from the reference geometry
        :type split_scan_mesh: bool
        :param optimization_settings: extra options forwarded to `optimize` for each relaxation step
        :type optimization_settings: dict
        :return: `(did_optimize, geometries, metadata)` -- per-point convergence flags, the resulting geometries (ordered to match `displacement_specs`), and per-point optimizer metadata
        :rtype: tuple
        :raises NotImplementedError: if a Cartesian-only displacement generator is used without internal-coordinate handling, or if `coordinate_expansion` is combined with a custom displacement-generator pair
        """
        ...

    def _optimize_along_displacements(self, coords, displacement_meshes, displacement_function, zigzag=True, mesh_iterator=None, **optimization_settings):
        """
        **LLM Docstring**

        Walk a mesh of displacement values (in zig-zag or plain product order), building each successive displaced structure from the previous optimized geometry (via `displacement_function`) and re-optimizing it, accumulating the sequence of optimization results.

        :param coords: the starting coordinates for the first displacement
        :type coords: np.ndarray
        :param displacement_meshes: the per-coordinate displacement value sequences to iterate over
        :type displacement_meshes: list[np.ndarray]
        :param displacement_function: builds a displaced structure from the current coordinates and a displacement tuple
        :type displacement_function: callable
        :param zigzag: whether to traverse the mesh in a zig-zag order (rather than a plain nested-loop product) for smoother successive displacements
        :type zigzag: bool
        :param mesh_iterator: an explicit iterator over displacement tuples, bypassing the zig-zag/product construction
        :type mesh_iterator: Iterator | None
        :param optimization_settings: extra options forwarded to `optimize` at each mesh point
        :type optimization_settings: dict
        :return: the list of `(converged, coords, info)` results, one per mesh point
        :rtype: list[tuple]
        """
        ...

@EnergyEvaluator.register('rdkit')
class RDKitEnergyEvaluator(EnergyEvaluator):

    def __init__(self, rdmol, force_field='mmff', charge=None, multiplicity=None, **defaults):
        ...

    @classmethod
    def from_mol(cls, mol, **opts):
        ...
    property_units = 'Kilocalories/Mole'
    analytic_derivative_order = 1

    def evaluate_term(self, coords, order, force_field_type=None, **opts):
        ...

    def optimize(self, coords, method=None, force_field_type=None, unitary=False, tol=1e-08, max_iterations=100, **opts):
        ...

@EnergyEvaluator.register('aimnet2')
class AIMNet2EnergyEvaluator(EnergyEvaluator):
    """
    Borrows structure from AIMNet2ASE to call appropriately
    """

    def __init__(self, atoms, model='aimnet2', charge=0, device=None, model_dir=None, multiplicity=None, quiet=True, **defaults):
        ...

    @classmethod
    def handle_specialization(cls, tag):
        ...

    @classmethod
    def from_mol(cls, mol, **opts):
        ...

    def to_ase(self, gradient_modification_function=None, gradient_modification_mode='shift', convert_modification_distances=True, convert_modification_energies=True, orthogonal_projection_generator=None, **etc):
        ...

    def to_pysis(self, gradient_modification_function=None, gradient_modification_mode='shift', convert_modification_distances=False, convert_modification_energies=False, orthogonal_projection_generator=None, **etc):
        ...

    @classmethod
    def _copy_pysis(cls, base_calc, **opts):
        ...

    def _maybe_download_asset(self, file: str, url: str):
        ...

    @classmethod
    @contextlib.contextmanager
    def _overload_aimnet_modeldir(cls, root_dir):
        ...

    @classmethod
    def setup_aimnet(cls, model, device=None, model_dir=None):
        ...

    @staticmethod
    def autodiff(expr, coord, pad_dim, create_graph=False, retain_graph=False):
        ...
    property_units = 'ElectronVolts'
    batched_orders = True
    analytic_derivative_order = 2

    def prep_eval(self, coords, order, keep_graph=None, forces=None, hessian=None, **opts):
        ...

    def _apply_autodiff(self, data, root_key, deriv_key, order, property_shape, diff_var=None, key_tag='deriv', pad_coord=1, create_graph=None, retain_graph=None):
        ...

    def _replace_nans(self, e, replacement):
        ...

    def process_aimnet_derivs(self, base_shape, coords, data, root_key, order, extra_deriv_coord=None, deriv_key=None, property_shape=(), pad_coord=1, replace_nans=True, output_shape=None):
        ...

    def evaluate_term(self, coords, order, **opts):
        ...

@EnergyEvaluator.register('ase')
class ASECalcEnergyEvaluator(EnergyEvaluator):

    def __init__(self, atoms, charge=0, multiplicity=1, quiet=True, embedding=None, **defaults):
        ...

    def to_ase(self, gradient_modification_function=None, gradient_modification_mode='shift', convert_modification_distances=True, convert_modification_energies=True, orthogonal_projection_generator=None, **etc):
        ...

    @classmethod
    def from_mol(cls, mol, **opts):
        ...

    @classmethod
    def setup_calc(cls, *, calc, **settings):
        ...
    property_units = 'ElectronVolts'
    batched_orders = True
    analytic_derivative_order = 2

    def prep_ase(self, coords):
        ...

    def evaluate_term(self, coords, order, **opts):
        ...

@EnergyEvaluator.register('mace')
class MACEEnergyEvaluator(ASECalcEnergyEvaluator):
    analytic_derivative_order = 2

    @classmethod
    def handle_specialization(cls, tag):
        ...

    @classmethod
    @contextlib.contextmanager
    def _overload_mace_modeldir(cls, root_dir):
        ...

    @classmethod
    def setup_calc(cls, model='extra_large', model_type=None, device=None, model_dir=None, **settings):
        ...

@EnergyEvaluator.register('uma')
class UMAEnergyEvaluator(ASECalcEnergyEvaluator):

    @classmethod
    def handle_specialization(cls, tag):
        ...

    @classmethod
    def setup_calc(cls, model='uma-s-1p2', device=None, task_name=None, model_dir=None, login=None, hf_token=None, predict_unit_options=None, **settings):
        ...

@EnergyEvaluator.register('upet')
class UPETEnergyEvaluator(ASECalcEnergyEvaluator):

    @classmethod
    def handle_specialization(cls, tag):
        ...

    @classmethod
    def setup_calc(cls, model='pet-mad-s', device=None, version=None, login=None, hf_token=None, model_dir=None, **settings):
        ...

@EnergyEvaluator.register('orb')
class OrbModelEnergyEvaluator(ASECalcEnergyEvaluator):

    @classmethod
    def handle_specialization(cls, tag):
        ...

    @classmethod
    def setup_calc(cls, model='omol', device=None, version=None, login=None, **model_opts):
        ...

@EnergyEvaluator.register('hipnn')
class HIPNNEnergyEvaluator(EnergyEvaluator):
    property_key = 'energies'
    gradient_key = 'Grad'

    def __init__(self, atoms, model_dir, property_key=dev.default, gradient_key=dev.default, charge=0, multiplicity=None, quiet=True, **defaults):
        ...

    @classmethod
    def handle_specialization(cls, tag):
        ...

    @classmethod
    def from_mol(cls, mol, *, model_dir, **opts):
        ...

    @staticmethod
    def autodiff(expr, coord_list, pad_dim, create_graph=False, retain_graph=False):
        ...

    def process_derivatives(self, expr, coords, order):
        ...

    def setup_model(cls, model_dir, *, device=None, **opts):
        ...

    def print_model_properties(self):
        ...
    property_units = 'Kilocalories/Mole'
    batched_orders = True
    analytic_derivative_order = 1

    def process_results(self, base_shape, coords, outputs, order, property_key=None, gradient_key=None, property_shape=(), output_shape=None, **etc):
        ...

    def prep_eval(self, coords, order, **opts):
        ...

    def evaluate_term(self, coords, order, **opts):
        ...

@EnergyEvaluator.register('xtb')
class XTBEnergyEvaluator(EnergyEvaluator):
    """
    Uses XTB to calculate
    """

    def __init__(self, atoms, method='GFN2-xTB', charge=0, quiet=True, **defaults):
        ...

    @classmethod
    def handle_specialization(cls, tag):
        ...

    @classmethod
    def from_mol(cls, mol, **opts):
        ...

    def get_single_point(self, coords):
        ...
    distance_units = 'BohrRadius'
    property_units = 'Hartrees'
    batched_orders = True
    analytic_derivative_order = 1

    def evaluate_term(self, coords, order, **opts):
        ...

@EnergyEvaluator.register('pyscf')
class PySCFEnergyEvaluator(EnergyEvaluator):
    """
    """

    def __init__(self, atoms, *, level_of_theory, basis_set, caller=None, charge=None, multiplicity=None, quiet=True, molecule_options=None, level_of_theory_options=None, disp=None, **defaults):
        ...

    @classmethod
    def handle_specialization(cls, tag):
        ...
    _default_caller = 'dft'

    @property
    def caller_dispatch(self) -> dev.OptionsMethodDispatch:
        ...

    def get_callers(self):
        ...

    def resolve_caller(self, caller, lot):
        ...

    @classmethod
    def from_mol(cls, mol, **opts):
        ...

    def get_molecule(self, coords):
        ...

    def _call_dft(self, molecule, restricted=True, functional_options=None, disp=None, **opts):
        ...

    def _call_hf(self, molecule, restricted=False, **opts):
        ...

    def _call_rhf(self, molecule, **opts):
        ...

    def _call_uhf(self, molecule, **opts):
        ...

    def _call_mp2(self, molecule, reference='hf', **etc):
        ...

    def get_energy_from_calc(self, calc) -> float:
        ...

    def get_gradient_from_calc(self, calc):
        ...

    def run_calculation(self, coords, **opts):
        ...
    distance_units = 'BohrRadius'
    property_units = 'Hartrees'
    batched_orders = True
    analytic_derivative_order = 1

    def evaluate_term(self, coords, order, **opts):
        ...

class EvaluationServerEnergyEvaluator(EnergyEvaluator):

    def __init__(self, *, request_handler, **defaults):
        ...

    @abc.abstractmethod
    def setup_request(self, coords, order=None, tasks='energy', **opts) -> (list[Any], dict[Any, Any], Any):
        ...

    @abc.abstractmethod
    def process_response(self, response, state):
        ...

    @abc.abstractmethod
    def cleanup_state(self, response, state):
        ...

    def evaluate_term(self, coords, order, tasks='energy', **opts):
        ...

    def optimize(self, coords, method='server', **opts):
        ...

@EnergyEvaluator.register('mlipserver')
class MLIPServerEnergyEvaluator(EvaluationServerEnergyEvaluator):
    distance_units = 'BohrRadius'

    def __init__(self, atoms, *, container_path, energy_evaluator, session: subprocess.Popen=None, port=None, temp_dir=None, launcher_options=None, bind_sources=True, container_env=None, container_mode=None, initialization_delay_time=1, charge=None, multiplicity=None, model_dir=dev.default, pass_tokens=True, **defaults):
        ...

    @classmethod
    def from_mol(cls, mol, **opts):
        ...

    @classmethod
    def pick_port(cls):
        ...

    def get_container_env_command(self):
        ...
    dynamic_sources = ['McUtils', 'Psience', 'mlipenv']
    token_list = ['HF_TOKEN']

    def get_launcher(self):
        ...

    def resolve_address(self):
        ...

    def launch_session(self):
        ...

    def cleanup_session(self):
        ...

    def __del__(self):
        ...

    def get_io_files(self):
        ...

    def setup_request(self, coords, order=None, tasks='energy', **opts) -> (list[Any], dict[Any, Any], Any):
        ...

    def process_response(self, response, state):
        ...

    def cleanup_state(self, response, state: list[str]):
        ...

    def _run_mlip_request(self, *args, **kwargs):
        ...

@EnergyEvaluator.register('function')
class PotentialFunctionEnergyEvaluator(EnergyEvaluator):

    def __init__(self, potential_function, property_units=None, energy_units='Hartrees', distance_units='BohrRadius', **opts):
        """
        **LLM Docstring**

        Wrap a bare potential-energy function as an `EnergyEvaluator`, delegating to `PropertyFunctionEvaluator.__init__` for the actual state setup.

        :param potential_function: the function to wrap, expected to accept `(coords, order=..., **opts)`
        :type potential_function: callable
        :param property_units: the units the function's output is given in; defaults to `energy_units`
        :type property_units: str | None
        :param energy_units: the default energy units if `property_units` isn't given
        :type energy_units: str
        :param distance_units: the units the function expects its input coordinates in
        :type distance_units: str
        :param opts: extra options forwarded to `PropertyFunctionEvaluator.__init__`
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...

    def evaluate_term(self, coords, order, **opts):
        """
        **LLM Docstring**

        Evaluate the wrapped potential function, delegating to `PropertyFunctionEvaluator.evaluate_term`.

        :param coords: the coordinates to evaluate at
        :type coords: np.ndarray
        :param order: the derivative order to request from the function
        :type order: int
        :param opts: extra options forwarded to the function
        :type opts: dict
        :return: the function's return value
        :rtype: object
        """
        ...

    def use_internal_coordinate_handlers(self):
        """
        **LLM Docstring**

        Whether this evaluator should route coordinate handling through the internal-coordinate machinery, delegating to `PropertyFunctionEvaluator.use_internal_coordinate_handlers`.

        :return: whether internal-coordinate handling applies
        :rtype: bool
        """
        ...

    def evaluate(self, coords, order=0, logger=None, fd_handler=None, **opts):
        """
        **LLM Docstring**

        Evaluate the wrapped potential function's expansion, delegating to `PropertyFunctionEvaluator.evaluate` (defaulting `fd_handler` to `internal_finite_difference_derivs` when working in internal coordinates).

        :param coords: the coordinates to evaluate at
        :type coords: np.ndarray
        :param order: the highest derivative order to compute
        :type order: int
        :param logger: accepted for interface consistency but not used directly in this method's body
        :type logger: Logger | None
        :param fd_handler: an explicit finite-difference handler to use
        :type fd_handler: callable | None
        :param opts: extra options forwarded to `PropertyFunctionEvaluator.evaluate`
        :type opts: dict
        :return: the evaluated energy expansion
        :rtype: list
        """
        ...

    @classmethod
    def get_property_function(cls, prop_func, mol, **opts):
        """
        **LLM Docstring**

        Resolve the potential function to wrap; on this class, simply returns `prop_func` unchanged alongside the given options.

        :param prop_func: the candidate potential function
        :type prop_func: callable
        :param mol: the molecule the evaluator is being built for (unused directly)
        :type mol: AbstractMolecule
        :param opts: extra options, passed through unchanged
        :type opts: dict
        :return: `(prop_func, opts)`
        :rtype: tuple
        """
        ...

    @classmethod
    def from_mol(cls, mol, property_function=None, potential_function=None, **opts):
        """
        **LLM Docstring**

        Build a `PotentialFunctionEnergyEvaluator` for a molecule, resolving the potential function from either `property_function` or `potential_function` and delegating to `PropertyFunctionEvaluator.initialize_from_mol`.

        :param mol: the molecule to build the evaluator for
        :type mol: AbstractMolecule
        :param property_function: the potential function to wrap
        :type property_function: callable | None
        :param potential_function: an alias for `property_function`, used if that isn't given
        :type potential_function: callable | None
        :param opts: extra options forwarded to `PropertyFunctionEvaluator.initialize_from_mol`
        :type opts: dict
        :return: the constructed evaluator
        :rtype: PotentialFunctionEnergyEvaluator
        """
        ...
    default_property_function = None

    @classmethod
    def bind_default(cls, potential):
        """
        **LLM Docstring**

        Register a default potential function to be used the next time `from_mol` is called without an explicit potential function.

        :param potential: the function to use as the default
        :type potential: callable
        :return: None
        :rtype: None
        """
        ...

@EnergyEvaluator.register('expansion')
class PotentialExpansionEnergyEvaluator(PotentialFunctionEnergyEvaluator):

    def __init__(self, expansion, center=None, ref=None, batched_orders=None, transforms=None, transformed_derivatives=False, **opts):
        """
        **LLM Docstring**

        Build an energy evaluator from a precomputed Taylor-series potential expansion (rather than a live function), wrapping it in a `PotentialSurface` if it isn't already callable.

        :param expansion: the potential expansion, either a callable surface object or raw derivative tensors to build a `PotentialSurface` from
        :type expansion: callable | list[np.ndarray]
        :param center: the reference/expansion center for the surface, if building one from raw derivatives
        :type center: np.ndarray | None
        :param ref: the reference energy value for the surface
        :type ref: float | None
        :param batched_orders: accepted but overridden to `True` (expansions always return all requested orders in one call)
        :type batched_orders: bool | None
        :param transforms: coordinate transformation(s) to apply when building the surface from raw derivatives
        :type transforms: list | None
        :param transformed_derivatives: whether the raw derivatives are already expressed in the transformed coordinate system
        :type transformed_derivatives: bool
        :param opts: extra options forwarded to `PotentialFunctionEnergyEvaluator.__init__`
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...

    @classmethod
    def from_mol(cls, mol, property_function=None, expansion=None, **opts):
        """
        **LLM Docstring**

        Build a `PotentialExpansionEnergyEvaluator` for a molecule, resolving the expansion from either `property_function` or `expansion`.

        :param mol: the molecule to build the evaluator for
        :type mol: AbstractMolecule
        :param property_function: the potential expansion to use
        :type property_function: object | None
        :param expansion: an alias for `property_function`, used if that isn't given
        :type expansion: object | None
        :param opts: extra options forwarded to the base class's `from_mol`
        :type opts: dict
        :return: the constructed evaluator
        :rtype: PotentialExpansionEnergyEvaluator
        """
        ...

    @classmethod
    def get_property_function(cls, expansion, mol, transforms=None, **ignored):
        """
        **LLM Docstring**

        Resolve the potential expansion to use: falls back to the molecule's own `potential_derivatives` if none is given, re-expresses the derivatives in a normal-mode basis (via the molecule's normal modes) if the cubic term appears to already be in a reduced mode basis, and wraps the result in a `PotentialSurface`.

        :param expansion: the candidate expansion (raw derivative tensors, `None` to use the molecule's own, or an already-callable surface)
        :type expansion: object | None
        :param mol: the molecule the evaluator is being built for
        :type mol: AbstractMolecule
        :param transforms: an explicit coordinate transformation to build the surface with; derived from the molecule's normal modes if needed and not given
        :type transforms: list | None
        :param ignored: extra options, passed through unchanged as the returned options
        :type ignored: dict
        :return: `(expansion, ignored)` -- the resolved (possibly newly built) `PotentialSurface` and the untouched extra options
        :rtype: tuple
        """
        ...

class DipoleEvaluator(PropertyEvaluator):
    evaluator_registry = {}
    target_property_units = ('ElementaryCharge', 'BohrRadius')

    @classmethod
    def get_evaluators(cls):
        """
        **LLM Docstring**

        The built-in named dipole-evaluator constructors: an expansion-based evaluator, plus adapters for HIPNN, AIMNet2, and RDKit.

        :return: the evaluator-name-to-constructor mapping
        :rtype: dict
        """
        ...

    @classmethod
    def get_default_function_evaluator_type(cls):
        """
        **LLM Docstring**

        The evaluator class used when a bare callable is supplied as a dipole-evaluator specification.

        :return: `DipoleFunctionDipoleEvaluator`
        :rtype: type
        """
        ...

class DipoleFunctionDipoleEvaluator(DipoleEvaluator):

    def __init__(self, potential_function, property_units=('ElementaryCharge', 'BohrRadius'), distance_units='BohrRadius', **opts):
        """
        **LLM Docstring**

        Wrap a bare dipole-evaluating function as a `DipoleEvaluator`, delegating to `PropertyFunctionEvaluator.__init__` for the actual state setup.

        :param potential_function: the function to wrap, expected to accept `(coords, order=..., **opts)` and return dipole values/derivatives
        :type potential_function: callable
        :param property_units: the units the function's output is given in
        :type property_units: tuple | str
        :param distance_units: the units the function expects its input coordinates in
        :type distance_units: str
        :param opts: extra options forwarded to `PropertyFunctionEvaluator.__init__`
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...

    def embed_coords(self, coords, **kwargs):
        """
        **LLM Docstring**

        Embed coordinates into this evaluator's working representation, delegating to `PropertyFunctionEvaluator.embed_coords`.

        :param coords: the coordinates to embed
        :type coords: np.ndarray
        :param kwargs: extra options forwarded to `PropertyFunctionEvaluator.embed_coords`
        :type kwargs: dict
        :return: the embedded coordinates
        :rtype: np.ndarray
        """
        ...

    def unembed_coords(self, coords, **kwargs):
        """
        **LLM Docstring**

        Convert coordinates back to plain Cartesians, delegating to `PropertyFunctionEvaluator.unembed_coords`.

        :param coords: the coordinates to unembed
        :type coords: np.ndarray
        :param kwargs: extra options forwarded to `PropertyFunctionEvaluator.unembed_coords`
        :type kwargs: dict
        :return: the unembedded (Cartesian) coordinates
        :rtype: np.ndarray
        """
        ...

    def unembed_derivs(self, base_coords, coords, derivs, **kwargs):
        """
        **LLM Docstring**

        Convert dipole derivative tensors back to plain Cartesian coordinates, delegating to `PropertyFunctionEvaluator.unembed_derivs`.

        :param base_coords: the reference coordinates the returned derivatives should be expressed relative to
        :type base_coords: np.ndarray
        :param coords: the coordinates the input derivatives were computed at
        :type coords: np.ndarray
        :param derivs: the derivative tensors to convert
        :type derivs: list[np.ndarray]
        :param kwargs: extra options forwarded to `PropertyFunctionEvaluator.unembed_derivs`
        :type kwargs: dict
        :return: the converted derivative tensors
        :rtype: list[np.ndarray]
        """
        ...

    def evaluate_term(self, coords, order, **opts):
        """
        **LLM Docstring**

        Evaluate the wrapped dipole function, delegating to `PropertyFunctionEvaluator.evaluate_term`.

        :param coords: the coordinates to evaluate at
        :type coords: np.ndarray
        :param order: the derivative order to request from the function
        :type order: int
        :param opts: extra options forwarded to the function
        :type opts: dict
        :return: the function's return value
        :rtype: object
        """
        ...

    def evaluate(self, coords, order=0, logger=None, **opts):
        """
        **LLM Docstring**

        Evaluate the wrapped dipole function's expansion, delegating to `PropertyFunctionEvaluator.evaluate`.

        :param coords: the coordinates to evaluate at
        :type coords: np.ndarray
        :param order: the highest derivative order to compute
        :type order: int
        :param logger: accepted for interface consistency but not used directly in this method's body
        :type logger: Logger | None
        :param opts: extra options forwarded to `PropertyFunctionEvaluator.evaluate`
        :type opts: dict
        :return: the evaluated dipole expansion
        :rtype: list
        """
        ...

    @classmethod
    def get_property_function(cls, prop_func, mol, **opts):
        """
        **LLM Docstring**

        Resolve the dipole function to wrap; on this class, simply returns `prop_func` unchanged alongside the given options.

        :param prop_func: the candidate dipole function
        :type prop_func: callable
        :param mol: the molecule the evaluator is being built for (unused directly)
        :type mol: AbstractMolecule
        :param opts: extra options, passed through unchanged
        :type opts: dict
        :return: `(prop_func, opts)`
        :rtype: tuple
        """
        ...

    @classmethod
    def from_mol(cls, mol, property_function=None, dipole_function=None, **opts):
        """
        **LLM Docstring**

        Build a `DipoleFunctionDipoleEvaluator` for a molecule, resolving the dipole function from either `property_function` or `dipole_function` and delegating to `PropertyFunctionEvaluator.initialize_from_mol`.

        :param mol: the molecule to build the evaluator for
        :type mol: AbstractMolecule
        :param property_function: the dipole function to wrap
        :type property_function: callable | None
        :param dipole_function: an alias for `property_function`, used if that isn't given
        :type dipole_function: callable | None
        :param opts: extra options forwarded to `PropertyFunctionEvaluator.initialize_from_mol`
        :type opts: dict
        :return: the constructed evaluator
        :rtype: DipoleFunctionDipoleEvaluator
        """
        ...
    default_property_function = None

    @classmethod
    def bind_default(cls, potential):
        """
        **LLM Docstring**

        Register a default dipole function to be used the next time `from_mol` is called without an explicit dipole function.

        :param potential: the function to use as the default
        :type potential: callable
        :return: None
        :rtype: None
        """
        ...

class DipoleExpansionDipoleEvaluator(DipoleFunctionDipoleEvaluator):

    def __init__(self, expansion, center=None, ref=None, property_units=('ElementaryCharge', 'BohrRadius'), distance_units='BohrRadius', analytic_derivative_order=None, batched_orders=None, **opts):
        """
        **LLM Docstring**

        Build a dipole evaluator from a precomputed Taylor-series dipole expansion (rather than a live function), wrapping it in a `DipoleSurface` if it isn't already callable and defaulting the analytic derivative order to the number of terms in the expansion.

        :param expansion: the dipole expansion, either a callable surface object or raw derivative tensors to build a `DipoleSurface` from
        :type expansion: callable | list[np.ndarray]
        :param center: the reference/expansion center for the surface
        :type center: np.ndarray | None
        :param ref: the reference dipole value for the surface
        :type ref: np.ndarray | None
        :param property_units: the units the dipole values are given in
        :type property_units: tuple
        :param distance_units: the units the coordinates are given in
        :type distance_units: str
        :param analytic_derivative_order: the highest order available analytically; defaults to the number of terms in `expansion`
        :type analytic_derivative_order: int | None
        :param batched_orders: accepted but overridden to `True`
        :type batched_orders: bool | None
        :param opts: extra options forwarded to `DipoleFunctionDipoleEvaluator.__init__`
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...

    @classmethod
    def expansion_from_mol_charges(cls, mol, charges=None, origin=None):
        """
        **LLM Docstring**

        Build a first-order dipole expansion (a constant term and its Cartesian gradient) directly from a set of atomic partial charges, treating the dipole as the charge-weighted sum of displacements from a reference origin.

        :param mol: the molecule to compute the dipole expansion for
        :type mol: AbstractMolecule
        :param charges: the per-atom partial charges; computed via `mol.calculate_charges()` if not given
        :type charges: np.ndarray | None
        :param origin: the reference origin for the dipole; defaults to the molecule's center of mass
        :type origin: np.ndarray | None
        :return: `[dipole_value, gradient]` -- the zeroth- and first-order dipole expansion terms
        :rtype: list[np.ndarray]
        """
        ...

    @classmethod
    def from_mol(cls, mol, property_function=None, expansion=None, **opts):
        """
        **LLM Docstring**

        Build a `DipoleExpansionDipoleEvaluator` for a molecule, resolving the dipole expansion from either `property_function` or `expansion`.

        :param mol: the molecule to build the evaluator for
        :type mol: AbstractMolecule
        :param property_function: the dipole expansion to use
        :type property_function: object | None
        :param expansion: an alias for `property_function`, used if that isn't given
        :type expansion: object | None
        :param opts: extra options forwarded to the base class's `from_mol`
        :type opts: dict
        :return: the constructed evaluator
        :rtype: DipoleExpansionDipoleEvaluator
        """
        ...

    @classmethod
    def get_property_function(cls, expansion, mol, center=None, transforms=None, use_modes=True, **ignored):
        """
        **LLM Docstring**

        Resolve the dipole expansion to use: falls back to the molecule's own `dipole_derivatives` (or a charge-based estimate via `expansion_from_mol_charges`, if neither is available) if none is given, optionally re-expresses the expansion through the molecule's normal modes (via its Hamiltonian's `dipole_expansion`), and wraps the result in a `DipoleSurface`.

        :param expansion: the candidate dipole expansion (raw derivative tensors, `None` to use the molecule's own/charge-derived one, or an already-callable surface)
        :type expansion: object | None
        :param mol: the molecule the evaluator is being built for
        :type mol: AbstractMolecule
        :param center: the reference/expansion center to build the surface with
        :type center: np.ndarray | None
        :param transforms: an explicit coordinate transformation to use; derived from the molecule's normal-mode embedding if `use_modes` is set and none is given
        :type transforms: list | None
        :param use_modes: whether to re-express the expansion in the molecule's normal-mode basis before building the surface
        :type use_modes: bool
        :param ignored: extra options, passed through unchanged as the returned options
        :type ignored: dict
        :return: `(expansion, ignored)` -- the resolved (possibly newly built) `DipoleSurface` and the untouched extra options
        :rtype: tuple
        """
        ...

class HIPNNDipoleEvaluator(DipoleEvaluator):

    def __init__(self, atoms, model_dir, charge=0, multiplicity=None, quiet=True, **defaults):
        ...

    @classmethod
    def from_mol(cls, mol, *, model_dir, **opts):
        ...
    property_units = ('ElementaryCharge', 'Angstroms')
    distance_units = 'Angstroms'
    batched_orders = True
    analytic_derivative_order = 0

    def evaluate_term(self, coords, order, **opts):
        ...

class ChargeEvaluator(PropertyEvaluator):
    evaluator_registry = {}
    target_property_units = 'ElementaryCharge'
    default_evaluator_type = 'rdkit'

    @classmethod
    def get_evaluators(cls):
        """
        **LLM Docstring**

        The built-in named charge-evaluator constructors: adapters for RDKit and AIMNet2.

        :return: the evaluator-name-to-constructor mapping
        :rtype: dict
        """
        ...

    @classmethod
    def get_default_function_evaluator_type(cls):
        """
        **LLM Docstring**

        The evaluator class used when a bare callable is supplied as a charge-evaluator specification.

        :return: `ChargeFunctionChargeEvaluator`
        :rtype: type
        """
        ...

class ChargeFunctionChargeEvaluator(ChargeEvaluator):

    def __init__(self, potential_function, property_units='ElementaryCharge', distance_units='BohrRadius', **opts):
        """
        **LLM Docstring**

        Wrap a bare charge-evaluating function as a `ChargeEvaluator`, delegating to `PropertyFunctionEvaluator.__init__` for the actual state setup.

        :param potential_function: the function to wrap, expected to accept `(coords, order=..., **opts)` and return partial-charge values/derivatives
        :type potential_function: callable
        :param property_units: the units the function's output is given in
        :type property_units: str
        :param distance_units: the units the function expects its input coordinates in
        :type distance_units: str
        :param opts: extra options forwarded to `PropertyFunctionEvaluator.__init__`
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...

    def evaluate_term(self, coords, order, **opts):
        """
        **LLM Docstring**

        Evaluate the wrapped charge function, delegating to `PropertyFunctionEvaluator.evaluate_term`.

        :param coords: the coordinates to evaluate at
        :type coords: np.ndarray
        :param order: the derivative order to request from the function
        :type order: int
        :param opts: extra options forwarded to the function
        :type opts: dict
        :return: the function's return value
        :rtype: object
        """
        ...

    def evaluate(self, coords, order=0, logger=None, **opts):
        """
        **LLM Docstring**

        Evaluate the wrapped charge function's expansion, delegating to `PropertyFunctionEvaluator.evaluate`.

        :param coords: the coordinates to evaluate at
        :type coords: np.ndarray
        :param order: the highest derivative order to compute
        :type order: int
        :param logger: accepted for interface consistency but not used directly in this method's body
        :type logger: Logger | None
        :param opts: extra options forwarded to `PropertyFunctionEvaluator.evaluate`
        :type opts: dict
        :return: the evaluated charge expansion
        :rtype: list
        """
        ...

    @classmethod
    def from_mol(cls, mol, property_function=None, **opts):
        """
        **LLM Docstring**

        Build a `ChargeFunctionChargeEvaluator` for a molecule, delegating to `PropertyFunctionEvaluator.initialize_from_mol`.

        :param mol: the molecule to build the evaluator for
        :type mol: AbstractMolecule
        :param property_function: the charge function to wrap
        :type property_function: callable | None
        :param opts: extra options forwarded to `PropertyFunctionEvaluator.initialize_from_mol`
        :type opts: dict
        :return: the constructed evaluator
        :rtype: ChargeFunctionChargeEvaluator
        """
        ...
    default_property_function = None

    @classmethod
    def bind_default(cls, potential):
        """
        **LLM Docstring**

        Register a default charge function to be used the next time `from_mol` is called without an explicit charge function.

        :param potential: the function to use as the default
        :type potential: callable
        :return: None
        :rtype: None
        """
        ...

class RDKitChargeEvaluator(ChargeEvaluator):

    def __init__(self, rdmol, model='gasteiger', **defaults):
        ...

    @classmethod
    def from_mol(cls, mol, charge=None, multiplicity=None, **opts):
        ...
    property_units = 'ElementaryCharge'
    analytic_derivative_order = 0

    def evaluate_term(self, coords, order, model=None, **opts):
        ...

class AIMNet2ChargeEvaluator(ChargeEvaluator):

    def __init__(self, atoms, model='aimnet2', charge=0, multiplicity=None, quiet=True, **defaults):
        ...

    @classmethod
    def from_mol(cls, mol, **opts):
        ...
    property_units = 'ElementaryCharge'
    distance_units = 'Angstroms'
    batched_orders = True
    analytic_derivative_order = 6
    scaling_factor = 1

    def evaluate_term(self, coords, order, **opts):
        ...

class ChargeEvaluatorDipoleEvaluator(DipoleEvaluator):

    def __init__(self, charge_evaluator: 'ChargeEvaluator', **etc):
        ...
    batched_orders = True

    def evaluate_term(self, coords, order, **opts):
        ...

class AIMNet2DipoleEvaluator(ChargeEvaluatorDipoleEvaluator):

    @classmethod
    def from_mol(cls, mol, **opts):
        ...

class RDKitDipoleEvaluator(ChargeEvaluatorDipoleEvaluator):

    @classmethod
    def from_mol(cls, mol, **opts):
        ...

class DipolePolarizabilityEvaluator(PropertyEvaluator):
    target_property_units = 'ElementaryCharge'

    @classmethod
    def get_evaluators(cls):
        """
        **LLM Docstring**

        The built-in named dipole-polarizability-evaluator constructors: adapters for AIMNet2 and HIPNN.

        :return: the evaluator-name-to-constructor mapping
        :rtype: dict
        """
        ...

    @classmethod
    def get_default_function_evaluator_type(cls):
        """
        **LLM Docstring**

        The evaluator class used when a bare callable is supplied as a dipole-polarizability-evaluator specification.

        :return: `PolarizabilityFunctionDipolePolarizabilityEvaluator`
        :rtype: type
        """
        ...

    @abc.abstractmethod
    def evaluate_polarizability_term(self, coords, coord_order, electrostatic_order, **opts):
        """
        **LLM Docstring**

        Abstract hook for evaluating the dipole-polarizability tensor's value (or specific coordinate/electrostatic-field derivative orders) at the given coordinates. Concrete evaluator subclasses must implement this.

        :param coords: the coordinates to evaluate at
        :type coords: np.ndarray
        :param coord_order: the coordinate-derivative order to evaluate
        :type coord_order: int
        :param electrostatic_order: the electrostatic-field-derivative order to evaluate
        :type electrostatic_order: int
        :param opts: extra evaluator-specific options
        :type opts: dict
        :return: the evaluated term(s)
        :rtype: object
        """
        ...
    multi_expansion_order = 2

    def evaluate_term(self, coords, order, **opts):
        """
        **LLM Docstring**

        Evaluate the polarizability's coordinate-derivative expansion at a given coordinate order, delegating to `evaluate_polarizability_term` with a fixed electrostatic order (`self.multi_expansion_order`) unless `order` is already given as an explicit `(coord_order, electrostatic_order)` pair.

        :param coords: the coordinates to evaluate at
        :type coords: np.ndarray
        :param order: the coordinate-derivative order, or an explicit `(coord_order, electrostatic_order)` pair
        :type order: int | tuple[int, int]
        :param opts: extra options forwarded to `evaluate_polarizability_term`
        :type opts: dict
        :return: the evaluated term(s)
        :rtype: object
        """
        ...

class PolarizabilityFunctionDipolePolarizabilityEvaluator(DipolePolarizabilityEvaluator):
    ...

class FluctuatingChargeDipolePolarizabilityEvaluator(DipolePolarizabilityEvaluator):

    @classmethod
    def mu_factor(cls, total_charge, X, J_inv):
        ...

    @classmethod
    def dipole_factor_expansion(cls, coords, total_charge, X, J_inv):
        ...

    @classmethod
    def polarizability_factor(cls, coords, J_inv):
        ...

    @abc.abstractmethod
    def evaluate_electrostatic_potential_derivatives(self, coords, coord_order, electrostatic_order, **opts):
        ...
    batched_orders = True

    def evaluate_polarizability_term(self, coords, coord_order, electrostatic_order, *, charge, **opts):
        ...

    @classmethod
    def expand_electrostatic_potential(cls, total_charge, coords, potential_by_charges_by_coords):
        ...

    @classmethod
    def _distorted_harmonic_polarizabilities(cls, charges, dipole_expansion, coords, hessians):
        ...

class AIMNet2DipolePolarizabilityEvaluator(FluctuatingChargeDipolePolarizabilityEvaluator):

    def __init__(self, atoms, model='aimnet2', charge=0, multiplicity=None, quiet=True, **defaults):
        ...

    @classmethod
    def from_mol(cls, mol, **opts):
        ...

    def evaluate_polarizability_term(self, coords, coord_order, electrostatic_order, *, charge=None, **opts):
        ...

    def evaluate_electrostatic_potential_derivatives(self, coords, coord_order, electrostatic_order, **opts):
        ...

class HIPNNDipolePolarizabilityEvaluator(FluctuatingChargeDipolePolarizabilityEvaluator):

    def __init__(self, atoms, model_dir, **defaults):
        ...

    @classmethod
    def _set_extra_outputs(cls, base):
        ...

    @classmethod
    def from_mol(cls, mol, *, model_dir, embedding=None, charge=None, multiplicity=None, **opts):
        ...

    def evaluate_polarizability_term(self, coords, coord_order, electrostatic_order, *, charge=None, **opts):
        ...

    def evaluate_electrostatic_potential_derivatives(self, coords, coord_order, electrostatic_order, **opts):
        ...

class ReducedDimensionalPotentialHandler:

    def __init__(self, mol):
        """
        **LLM Docstring**

        Store the molecule this handler will build reduced-dimensionality (1D) potential slices for.

        :param mol: the molecule to build potentials for
        :type mol: AbstractMolecule
        :return: None
        :rtype: None
        """
        ...

    @classmethod
    def get_potential_params(cls, spec, re, g, local_derivs, quartic_potential_cutoff=0.05, poly_expansion_order=2, **opts):
        """
        **LLM Docstring**

        Fit a simple analytic 1D potential form (a Morse potential, or a low-order polynomial in a Duschinsky-like `w`-scaled coordinate) to a coordinate's local derivative expansion, choosing Morse when the cubic-to-quartic-derivative ratio exceeds `quartic_potential_cutoff` and a polynomial fit otherwise.

        :param spec: the coordinate specification the fit is for (passed through into the returned options for later use by a potential constructor)
        :type spec: object
        :param re: the coordinate's equilibrium/reference value
        :type re: float
        :param g: the coordinate's G-matrix element (effective inverse mass)
        :type g: float
        :param local_derivs: the local derivative expansion (value, 1st, 2nd, 3rd, ... derivatives) of the potential along this coordinate
        :type local_derivs: list[float]
        :param quartic_potential_cutoff: the cubic/quartic-derivative ratio threshold below which a polynomial (rather than Morse) fit is used
        :type quartic_potential_cutoff: float
        :param poly_expansion_order: the polynomial expansion order to keep when fitting the `'poly'` form
        :type poly_expansion_order: int
        :param opts: extra options passed through into the returned parameter dict
        :type opts: dict
        :return: the fitted potential parameters (a `'method'` key of `'morse'` or `'poly'`, plus the corresponding fit parameters and `re`)
        :rtype: dict
        """
        ...

    @classmethod
    def get_morse_potential(cls, spec, *, w, wx, g, re):
        """
        **LLM Docstring**

        Build a Morse-potential `CoordinateFunction` for a given coordinate specification and fitted parameters.

        :param spec: the coordinate specification to build the potential over
        :type spec: object
        :param w: the Morse frequency parameter
        :type w: float
        :param wx: the Morse anharmonicity parameter
        :type wx: float
        :param g: the coordinate's G-matrix element
        :type g: float
        :param re: the equilibrium/reference coordinate value
        :type re: float
        :return: the constructed Morse potential function
        :rtype: CoordinateFunction
        """
        ...

    @classmethod
    def get_poly_potential(cls, spec, *, w_coeffs=None, coeffs=None, g, re):
        """
        **LLM Docstring**

        Build a polynomial `CoordinateFunction` for a given coordinate specification, either from raw Taylor coefficients or from `w`-scaled (Duschinsky-like) coefficients that are first converted into raw coefficients.

        :param spec: the coordinate specification to build the potential over
        :type spec: object
        :param w_coeffs: `w`-scaled polynomial coefficients to convert into raw coefficients, if `coeffs` isn't given directly
        :type w_coeffs: list[float] | None
        :param coeffs: raw Taylor-series coefficients for the polynomial, used directly if given
        :type coeffs: list[float] | None
        :param g: the coordinate's G-matrix element, used to convert `w_coeffs` into raw coefficients
        :type g: float
        :param re: the equilibrium/reference coordinate value (the polynomial's expansion center)
        :type re: float
        :return: the constructed polynomial potential function
        :rtype: CoordinateFunction
        """
        ...

    @classmethod
    def get_potential_types(cls):
        """
        **LLM Docstring**

        The mapping from potential-form name to the constructor method that builds that form, used by the potential-generator dispatch.

        :return: the name-to-constructor mapping (`'morse'`, `'poly'`)
        :rtype: dict
        """
        ...

    @classmethod
    def get_potentials_by_attributes(cls):
        """
        **LLM Docstring**

        Attribute-based dispatch mapping used to infer which potential form to build from the keys present in a parameter dict (e.g. `w`/`wx` implies Morse; `coeffs`/`w_coeffs` implies polynomial).

        :return: the attribute-tuple-to-form-name mapping
        :rtype: dict
        """
        ...
    _pot_dispatch = dev.uninitialized
    default_evaluator_type = 'poly'

    @classmethod
    def potential_generator_dispatch(cls) -> 'dev.OptionsMethodDispatch':
        """
        **LLM Docstring**

        Lazily build (and cache) the `OptionsMethodDispatch` object used to resolve a parameter dict into the right potential-form constructor and its arguments.

        :return: the cached dispatch object
        :rtype: dev.OptionsMethodDispatch
        """
        ...

    def get_g(self, which, return_tf=True, order=1):
        """
        **LLM Docstring**

        Compute the G-matrix (effective inverse-mass) element(s) for a given internal-coordinate specification at the molecule's current geometry, via the mass-weighted B-matrix.

        :param which: the internal coordinate(s) to compute the G-matrix for
        :type which: object
        :param return_tf: whether to also return the raw coordinate values and coordinate-transformation Jacobians used
        :type return_tf: bool
        :param order: the Jacobian derivative order to compute
        :type order: int
        :return: the G-matrix (submatrix restricted to `which`), or `(G, (r, tf))` including the raw coordinate values/Jacobians if `return_tf` is set
        :rtype: np.ndarray | tuple
        """
        ...

    def get_anharmonic_parameters(self, which, evaluator=None, energy_expansion=None, params_handler=None, **opts):
        """
        **LLM Docstring**

        Fit 1D analytic potential parameters for each of a set of internal coordinates, by extracting each coordinate's diagonal local derivative expansion (from the full Cartesian potential expansion, re-expressed through the coordinate Jacobian) and fitting it via `params_handler` (defaulting to `get_potential_params`).

        :param which: the internal coordinate(s) to fit potentials for
        :type which: object
        :param evaluator: an explicit energy evaluator to use when computing the potential derivatives
        :type evaluator: object | None
        :param energy_expansion: a precomputed Cartesian potential-derivative expansion to use instead of computing a fresh one
        :type energy_expansion: list[np.ndarray] | None
        :param params_handler: the function used to fit each coordinate's local derivative expansion into potential parameters; defaults to `get_potential_params`
        :type params_handler: callable | None
        :param opts: extra options forwarded to `params_handler`
        :type opts: dict
        :return: the list of fitted parameter dicts, one per coordinate in `which`
        :rtype: list[dict]
        """
        ...

    def get_1d_potentials(self, which, evaluator=None, energy_expansion=None, params_handler=None, potential_params=None, coordinate_potential_handler=None, **opts):
        """
        **LLM Docstring**

        Build 1D analytic potential-energy functions for a set of internal coordinates, using any explicitly supplied per-coordinate parameters where given and automatically fitting the rest (via `get_anharmonic_parameters`), then constructing each potential via the resolved potential-form dispatch.

        :param which: the internal coordinate(s) to build potentials for, or a dict mapping coordinates to explicit parameter dicts
        :type which: object | dict
        :param evaluator: an explicit energy evaluator to use when auto-fitting missing parameters
        :type evaluator: object | None
        :param energy_expansion: a precomputed Cartesian potential-derivative expansion to reuse when auto-fitting
        :type energy_expansion: list[np.ndarray] | None
        :param params_handler: the function used to auto-fit missing coordinate parameters; forwarded to `get_anharmonic_parameters`
        :type params_handler: callable | None
        :param potential_params: explicit per-coordinate parameter dicts (list or dict form) to use instead of auto-fitting, for coordinates where available
        :type potential_params: list | dict | None
        :param coordinate_potential_handler: an explicit potential-form dispatch object to use instead of `potential_generator_dispatch()`
        :type coordinate_potential_handler: dev.OptionsMethodDispatch | None
        :param opts: extra options forwarded to `get_anharmonic_parameters` for any auto-fitted coordinates
        :type opts: dict
        :return: the list of constructed 1D potential functions, one per coordinate in `which`
        :rtype: list[CoordinateFunction]
        """
        ...