"""
Stores all of the terms used inside the VPT2 representations
"""
import numpy as np, functools as fp, itertools, time, enum, math
from McUtils.Numputils import SparseArray, levi_cevita3
import McUtils.Numputils as nput
from McUtils.Data import UnitsData
from McUtils.Scaffolding import Logger, NullLogger, Checkpointer, NullCheckpointer
from McUtils.Parallelizers import Parallelizer
from McUtils.Zachary import TensorDerivativeConverter, TensorExpansionTerms
from McUtils.Combinatorics import UniquePermutations
from ..Molecools import Molecule, MolecularVibrations, MolecularNormalModes
from .Common import PerturbationTheoryException
__all__ = ['ExpansionTerms', 'KineticTerms', 'PotentialTerms', 'CoriolisTerm', 'PotentialLikeTerm', 'DipoleTerms', 'OperatorTerms']

class DumbTensor:
    """
    A wrapper to make tensor algebra suck less
    """

    def __init__(self, tensor):
        """
        **LLM Docstring**

        Wrap a raw tensor (e.g. a `np.ndarray`) so it supports the convenience arithmetic/reshaping operations defined on `DumbTensor`.

        :param tensor: the tensor to wrap
        :type tensor: np.ndarray
        :return: None
        :rtype: None
        """
        ...

    @property
    def shape(self):
        """
        **LLM Docstring**

        Shape of the wrapped tensor.

        :return: `self.t.shape`
        :rtype: tuple
        """
        ...

    @staticmethod
    def _dot(*t, axes=None):
        """
        Flexible tensordot
        """
        ...

    def dot(self, b, *args, **kwargs):
        """
        **LLM Docstring**

        Contract this tensor with `b` using `DumbTensor._dot`, unwrapping `b` first if it is itself a `DumbTensor`.

        :param b: the other tensor (or `DumbTensor`) to contract with
        :type b: np.ndarray | DumbTensor
        :param args: extra positional arguments forwarded to `_dot`
        :type args: tuple
        :param kwargs: extra keyword arguments forwarded to `_dot`, notably `axes`
        :type kwargs: dict
        :return: a new `DumbTensor` wrapping the contraction result
        :rtype: DumbTensor
        """
        ...

    @staticmethod
    def _shift(a, *s):
        """
        **LLM Docstring**

        Apply a sequence of axis-swap transpositions to `a`. Each element of `s` is an `(i, j)` pair; for each pair, the axis at position `i` is moved to sit immediately after (or as) position `j`, shifting the intervening axes accordingly. Integers are passed through unchanged.

        :param a: the tensor to transpose, or an `int` passed through unchanged
        :type a: np.ndarray | int
        :param s: one or more `(i, j)` axis-shift pairs to apply in order
        :type s: tuple[int, int]
        :return: the transposed tensor (or the original `int`)
        :rtype: np.ndarray | int
        """
        ...

    def shift(self, *args, **kwargs):
        """
        **LLM Docstring**

        Apply `_shift` to this tensor's data and wrap the result in a new `DumbTensor`.

        :param args: `(i, j)` axis-shift pairs forwarded to `_shift`
        :type args: tuple
        :param kwargs: forwarded to `_shift`
        :type kwargs: dict
        :return: a new `DumbTensor` with the shifted axes
        :rtype: DumbTensor
        """
        ...

    def transpose(self, *perm):
        """
        **LLM Docstring**

        Transpose the wrapped tensor according to `perm` and wrap the result in a new `DumbTensor`.

        :param perm: the axis permutation to apply
        :type perm: tuple[int, ...]
        :return: a new `DumbTensor` with the transposed data
        :rtype: DumbTensor
        """
        ...

    @staticmethod
    def _contract_dim(R, targ_dim):
        """
        **LLM Docstring**

        Collapse the trailing axes of `R` (assumed to already be flattened pairwise from the end inward) down to `targ_dim` total dimensions, by repeatedly reshaping pairs of trailing axes into one.

        :param R: the array whose trailing dimensions should be merged
        :type R: np.ndarray
        :param targ_dim: the desired number of dimensions after contraction
        :type targ_dim: int
        :return: the reshaped array with `targ_dim` dimensions
        :rtype: np.ndarray
        """
        ...

    def contract_dim(self, targ_dim):
        """
        **LLM Docstring**

        Apply `_contract_dim` to this tensor's data and wrap the result in a new `DumbTensor`.

        :param targ_dim: the desired number of dimensions after contraction
        :type targ_dim: int
        :return: a new `DumbTensor` with the reduced dimensionality
        :rtype: DumbTensor
        """
        ...

    def __add__(self, other):
        """
        **LLM Docstring**

        Elementwise addition, unwrapping `other` first if it is a `DumbTensor`.

        :param other: the value to add
        :type other: np.ndarray | DumbTensor
        :return: a new `DumbTensor` wrapping `self.t + other`
        :rtype: DumbTensor
        """
        ...

    def __radd__(self, other):
        """
        **LLM Docstring**

        Reflected addition; identical to `__add__` since addition here is commutative.

        :param other: the value to add
        :type other: np.ndarray | DumbTensor
        :return: a new `DumbTensor` wrapping `self.t + other`
        :rtype: DumbTensor
        """
        ...

    def __matmul__(self, other):
        """
        **LLM Docstring**

        Operator form of `dot`; lets `@` be used to contract two `DumbTensor`s (or a `DumbTensor` and a raw array).

        :param other: the tensor to contract with
        :type other: np.ndarray | DumbTensor
        :return: a new `DumbTensor` wrapping the contraction result
        :rtype: DumbTensor
        """
        ...

    def __getitem__(self, item):
        """
        :type item: slice
        """
        ...

class MixedDerivativeHandlingModes(enum.Enum):
    Unhandled = 'unhandled'
    Numerical = 'numerical'
    Analytical = 'analytical'
    Averaged = 'averaged'
    Old = 'old'

class ImaginaryFrequencyHandlingMode(enum.Enum):
    Abs = 'abs'
    Signed = 'signed'

class JacobianKeys(enum.Enum):
    """Real access pattern: JacobianKeys.<MemberName> (this is an enum with 12 members, e.g. JacobianKeys.CartesiansByInternals == 'CartesiansByInternals'). Collapsed into a dict below purely for compactness -- do not index it as a dict in real code:"""
    _MEMBERS = {'CartesiansByInternals': 'CartesiansByInternals', 'InternalsByCartesians': 'InternalsByCartesians', 'InternalsByCartesianModes': 'InternalsByModes', 'CartesianModesByInternals': 'ModesByInternals', 'CartesiansByInternalModes': 'CartesiansByModes', 'InternalModesByCartesians': 'ModesByCartesians', 'CartesianModesByInternalModes': 'CartesianModesByInternalModes', 'InternalModesByCartesianModes': 'InternalModesByCartesianModes', 'InternalModesByInternals': 'InternalModesByInternals', 'InternalsByInternalModes': 'InternalsByInternalModes', 'CartesianModesByCartesians': 'CartesianModesByCartesians', 'CartesiansByCartesianModes': 'CartesiansByCartesianModes'}

class ExpansionTerms:
    """
    Base class for kinetic, potential, and dipole derivative terms
    """
    __props__ = ('logger', 'parallelizer', 'checkpointer', 'undimensionalize', 'numerical_jacobians', 'eckart_embed_derivatives', 'eckart_embed_planar_ref_tolerance', 'strip_dummies', 'strip_embedding', 'mixed_derivative_handling_mode', 'mixed_derivative_warning_threshold', 'mixed_derivative_handle_zeros', 'backpropagate_internals', 'direct_propagate_cartesians', 'zero_mass_term', 'internal_fd_mesh_spacing', 'internal_fd_stencil', 'cartesian_fd_mesh_spacing', 'cartesian_fd_stencil', 'cartesian_analytic_deriv_order', 'cartesian_by_internal_derivative_method', 'internal_by_cartesian_order', 'cartesian_by_internal_order', 'expansion_handling_mode', 'jacobian_warning_threshold', 'coordinate_transformations', 'coordinate_derivatives', 'imaginary_frequency_handling_mode')
    _cached_jacobians = {}

    def __init__(self, molecule, modes=None, mode_selection=None, mode_transformation=None, use_internal_modes=None, logger=None, parallelizer=None, checkpointer=None, undimensionalize=None, numerical_jacobians=True, eckart_embed_derivatives=True, eckart_embed_planar_ref_tolerance=None, strip_dummies=False, strip_embedding=True, mixed_derivative_handling_mode=None, mixed_derivative_warning_threshold=0.00025, mixed_derivative_handle_zeros=False, backpropagate_internals=False, direct_propagate_cartesians=False, zero_mass_term=10000000.0, expansion_handling_mode='old', internal_fd_mesh_spacing=0.01, internal_fd_stencil=None, cartesian_fd_mesh_spacing=0.01, cartesian_fd_stencil=None, cartesian_analytic_deriv_order=None, cartesian_by_internal_derivative_method=None, internal_by_cartesian_order=3, cartesian_by_internal_order=4, jacobian_warning_threshold=10000.0, coordinate_transformations=None, coordinate_derivatives=None, imaginary_frequency_handling_mode='abs'):
        """
        :param molecule: the molecule we're doing the expansion for
        :type molecule: Molecule
        :param modes: normal modes in Cartesian coordinates
        :type modes: MolecularVibrations
        :param mode_selection: the selection of modes to use
        :type mode_selection: None | Iterable[int]
        :param undimensionalize: whether or not we need to do some units fuckery on the modes
        :type undimensionalize: bool
        """
        ...

    @property
    def num_atoms(self):
        """
        Gets the number of atoms (excluding dummies if `strip_dummies` is `True`)

        :return:
        :rtype:
        """
        ...

    def _check_internal_modes(self, modes=None, clean=True):
        """
        **LLM Docstring**

        Determine (and cache) whether the stored mode matrix is expressed in internal coordinates rather than Cartesians, by comparing its row count against the expected Cartesian dimension (`3*natoms - 6`); if internal and `clean` is set, reshapes the mode matrices to restore any stripped embedding coordinates.

        :param modes: the modes to check; defaults to `self._modes`
        :type modes: MixtureModes | None
        :param clean: whether to reshape/pad the mode matrices via `_reshape_internal_modes` if they're found to be internal-coordinate-basis
        :type clean: bool
        :return: whether the modes are internal-coordinate-basis
        :rtype: bool
        """
        ...

    def _reshape_internal_modes(self):
        """
        **LLM Docstring**

        Pad the stored mode-by-coordinate and coordinate-by-mode matrices with zero rows/columns for the fixed embedding coordinates, if they were built with those coordinates stripped out, restoring them to the full `3*natoms`-dimensional Cartesian-adjacent shape expected elsewhere.

        :return: None
        :rtype: None
        """
        ...

    @property
    def modes(self):
        """
        **LLM Docstring**

        The stored mode object (in whatever basis -- Cartesian or internal -- it was constructed with).

        :return: the mode object
        :rtype: MixtureModes
        """
        ...

    def _tripmass(self, masses):
        """
        **LLM Docstring**

        Mass-weighting helper: optionally drops dummy-atom masses (when `self.strip_dummies` is set) or replaces them with `self.zero_mass_term`, then repeats each mass three times (once per Cartesian direction) and flattens to a length `3*n_atoms` vector.

        :param masses: per-atom masses
        :type masses: np.ndarray
        :return: the per-Cartesian-coordinate mass vector
        :rtype: np.ndarray
        """
        ...

    def get_terms(self, order=None):
        """
        Gets the terms up to the given order

        :param order:
        :type order:
        :return:
        :rtype:
        """
        ...

    def get_term(self, t):
        """
        Provides the term at order `t`

        :param t:
        :type t:
        :return:
        :rtype:
        """
        ...

    @property
    def terms(self):
        """
        **LLM Docstring**

        The (cached) full set of expansion terms, computed lazily via `get_terms()` the first time they're needed.

        :return: the expansion terms
        :rtype: list[np.ndarray]
        """
        ...

    def __getitem__(self, item):
        """
        **LLM Docstring**

        Fetch a single term at the given order, via `get_term`.

        :param item: the order to fetch
        :type item: int
        :return: the term at that order
        :rtype: np.ndarray
        """
        ...

    @staticmethod
    def _weight_derivatives(t, order=None):
        """
        **LLM Docstring**

        Apply the standard symmetric-derivative combinatorial weighting (`1/k!` on each `k`-fold-repeated-index diagonal slice) to a Taylor-expansion derivative tensor, needed because raw finite-difference/analytic derivative tensors don't yet include the Taylor-series `1/n!` factors along their repeated-index diagonals.

        :param t: the derivative tensor to weight, or an `int` passed through unchanged
        :type t: np.ndarray | int
        :param order: the derivative order (number of tensor axes) to weight; inferred from `t.shape` if not given
        :type order: int | None
        :return: the weighted tensor (or `t` unchanged, if it was an `int`)
        :rtype: np.ndarray | int
        """
        ...

    def _freq_sqrt(self, freqs, cutoff=1e-06):
        """
        **LLM Docstring**

        Compute the square root of a set of (possibly negative/imaginary) frequencies, either always taking the magnitude first (`Abs` mode) or preserving the sign of the result for frequencies whose magnitude exceeds `cutoff` (`Signed` mode), depending on `self.imaginary_frequency_handling`.

        :param freqs: the frequencies to take the square root of
        :type freqs: np.ndarray
        :param cutoff: the magnitude threshold below which a frequency is treated as exactly zero (and thus not sign-adjusted) in `Signed` mode
        :type cutoff: float
        :return: the (possibly signed) square roots
        :rtype: np.ndarray
        :raises ValueError: if `self.imaginary_frequency_handling` isn't a recognized `ImaginaryFrequencyHandlingMode`
        """
        ...

    def get_int_jacobs(self, jacs):
        """
        Gets the specified Internal->Cartesian Jacobians

        :param jacs:
        :type jacs:
        :return:
        :rtype:
        """
        ...

    def get_cart_jacobs(self, jacs):
        """
        Gets the specified Cartesian->Internal Jacobians

        :param jacs:
        :type jacs:
        :return:
        :rtype:
        """
        ...

    @property
    def inertial_frame(self):
        """
        Provides the inertial axis frame

        :return:
        :rtype:
        """
        ...

    def inertial_frame_derivatives(self):
        """
        **LLM Docstring**

        Compute the first and second derivatives of the (mass-weighted) inertia tensor with respect to mass-weighted Cartesian displacements, using closed-form tensor expressions rather than finite differences.

        :return: `[I0Y, I0YY]`, the first derivative tensor (shape `(3*nAt, 3, 3)`) and second derivative tensor (shape `(3*nAt, 3*nAt, 3, 3)`) of the inertia tensor with respect to mass-weighted Cartesian displacements
        :rtype: list[np.ndarray]
        """
        ...

    def moment_of_inertia_derivs(self, order):
        """
        **LLM Docstring**

        Compute the Taylor-series derivatives of the inverse inertia tensor with respect to the normal-mode coordinates, up to the requested order, via the recursive relation built from the first-order Cartesian inertia-tensor derivative re-expressed in mode coordinates.

        :param order: the highest derivative order to compute
        :type order: int
        :return: the reciprocal-inertia-tensor expansion terms, `order + 1` entries starting from the inertia tensor itself
        :rtype: list[np.ndarray]
        """
        ...

    def _get_embedding_coords(self):
        """
        **LLM Docstring**

        Look up the indices of the coordinates used purely for embedding (e.g. translation/rotation) from the internal coordinate system, trying `system.embedding_coords` first and falling back to `system.converter_options['embedding_coords']`.

        :return: the embedding-coordinate indices, or `None` if neither source defines them
        :rtype: np.ndarray | None
        """
        ...
    _cached_transforms = {}

    def get_coordinate_transforms(self, internal_by_cartesian_order=None, cartesian_by_internal_order=None, current_cache=None):
        """
        **LLM Docstring**

        Compute (and cache, per-molecule, both in memory and via the checkpointer) the full set of Jacobians relating Cartesian coordinates, internal coordinates, Cartesian normal modes, and internal-coordinate-basis normal modes to each other, up to the requested derivative orders: computes the internals-by-Cartesians Jacobians (mass-weighting them and warning about/zeroing any anomalously large entries), then chains them together with the mode transformation via `TensorDerivativeConverter` to populate every entry of `JacobianKeys`.

        :param internal_by_cartesian_order: derivative order (number of Cartesian derivatives) to compute for internals-by-Cartesians Jacobians; defaults to `self.internal_by_cartesian_order`
        :type internal_by_cartesian_order: int | None
        :param cartesian_by_internal_order: derivative order (number of internal derivatives) to compute for Cartesians-by-internals Jacobians; defaults to `self.cartesian_by_internal_order`
        :type cartesian_by_internal_order: int | None
        :param current_cache: an existing partial cache to extend, instead of the per-molecule cache (in-memory or checkpointed) that would otherwise be looked up
        :type current_cache: dict | None
        :return: the (possibly newly extended) cache mapping each `JacobianKeys` member to a list of Jacobian tensors by order
        :rtype: dict
        """
        ...

    @property
    def cartesian_L_matrix(self):
        """
        **LLM Docstring**

        First-order Cartesians-by-Cartesian-normal-modes transformation matrix.

        :return: the leading term of `get_cartesians_by_cartesian_modes(1)`
        :rtype: np.ndarray
        """
        ...

    def get_cartesians_by_cartesian_modes(self, order=None):
        """
        **LLM Docstring**

        Fetch the Cartesians-by-Cartesian-normal-modes Jacobians up to the requested order, computing them (via `get_coordinate_transforms`) if not already cached.

        :param order: number of derivative orders to return; if `None`, all currently cached orders are returned
        :type order: int | None
        :return: list of Jacobian tensors, one per derivative order
        :rtype: list[np.ndarray]
        :raises ValueError: if fewer cached orders are available than requested
        """
        ...

    @property
    def cartesian_L_inverse(self):
        """
        **LLM Docstring**

        First-order Cartesian-normal-modes-by-Cartesians transformation matrix.

        :return: the leading term of `get_cartesian_modes_by_cartesians(1)`
        :rtype: np.ndarray
        """
        ...

    def get_cartesian_modes_by_cartesians(self, order=None):
        """
        **LLM Docstring**

        Fetch the Cartesian-normal-modes-by-Cartesians Jacobians up to the requested order, computing them if not already cached.

        :param order: number of derivative orders to return; if `None`, all currently cached orders are returned
        :type order: int | None
        :return: list of Jacobian tensors, one per derivative order
        :rtype: list[np.ndarray]
        :raises ValueError: if fewer cached orders are available than requested
        """
        ...

    @property
    def internal_L_matrix(self):
        """
        **LLM Docstring**

        First-order internal-normal-modes-by-internals transformation matrix.

        :return: the leading term of `get_internal_modes_by_internals(1)`
        :rtype: np.ndarray
        """
        ...

    def get_internal_modes_by_internals(self, order=None, strip_embedding=True):
        """
        **LLM Docstring**

        Fetch the internal-normal-modes-by-internals Jacobians up to the requested order, optionally stripping embedding-coordinate rows from the result.

        :param order: number of derivative orders to return; if `None`, all currently cached orders are returned
        :type order: int | None
        :param strip_embedding: whether to strip embedding-coordinate rows from the result (only applied if not already stripped globally via `self.strip_embedding`)
        :type strip_embedding: bool
        :return: list of Jacobian tensors, one per derivative order
        :rtype: list[np.ndarray]
        :raises ValueError: if fewer cached orders are available than requested
        """
        ...

    @property
    def internal_L_inverse(self):
        """
        **LLM Docstring**

        First-order internals-by-internal-normal-modes transformation matrix.

        :return: the leading term of `get_internals_by_internal_modes(1)`
        :rtype: np.ndarray
        """
        ...

    def get_internals_by_internal_modes(self, order=None, strip_embedding=True):
        """
        **LLM Docstring**

        Fetch the internals-by-internal-normal-modes Jacobians up to the requested order, optionally stripping embedding-coordinate columns from the result.

        :param order: number of derivative orders to return; if `None`, all currently cached orders are returned
        :type order: int | None
        :param strip_embedding: whether to strip embedding-coordinate columns from the result (only applied if not already stripped globally via `self.strip_embedding`)
        :type strip_embedding: bool
        :return: list of Jacobian tensors, one per derivative order
        :rtype: list[np.ndarray]
        :raises ValueError: if fewer cached orders are available than requested
        """
        ...

    @property
    def cartesians_by_modes(self):
        """
        **LLM Docstring**

        All cached Cartesians-by-internal-modes Jacobians, computing the default set if not already cached.

        :return: the `JacobianKeys.CartesiansByInternalModes` entry from `get_cartesians_by_modes()`
        :rtype: list[np.ndarray]
        """
        ...

    def get_cartesians_by_modes(self, order=None):
        """
        **LLM Docstring**

        Fetch the Cartesians-by-internal-normal-modes Jacobians up to the requested order, computing them if not already cached.

        :param order: number of derivative orders to return; if `None`, all currently cached orders are returned
        :type order: int | None
        :return: list of Jacobian tensors, one per derivative order
        :rtype: list[np.ndarray]
        :raises ValueError: if fewer cached orders are available than requested
        """
        ...

    @property
    def modes_by_cartesians(self):
        """
        **LLM Docstring**

        All cached internal-normal-modes-by-Cartesians Jacobians, computing the default set if not already cached.

        :return: the `JacobianKeys.InternalModesByCartesians` entry from `get_coordinate_transforms()`
        :rtype: list[np.ndarray]
        """
        ...

    def get_modes_by_cartesians(self, order=None, strip_embedding=True):
        """
        **LLM Docstring**

        Fetch the internal-normal-modes-by-Cartesians Jacobians up to the requested order.

        :param order: number of derivative orders to return; if `None`, all currently cached orders are returned
        :type order: int | None
        :param strip_embedding: accepted for interface consistency with sibling methods but not used in this method's body
        :type strip_embedding: bool
        :return: list of Jacobian tensors, one per derivative order
        :rtype: list[np.ndarray]
        :raises ValueError: if fewer cached orders are available than requested
        """
        ...

    @property
    def cartesians_by_internals(self):
        """
        **LLM Docstring**

        All cached Cartesians-by-internals Jacobians, computing the default set if not already cached.

        :return: the `JacobianKeys.CartesiansByInternals` entry from `get_coordinate_transforms()`
        :rtype: list[np.ndarray]
        """
        ...

    def get_cartesians_by_internals(self, order=None, strip_embedding=False):
        """
        **LLM Docstring**

        Fetch the Cartesians-by-internals Jacobians up to the requested order, optionally stripping embedding coordinates from every axis but the first.

        :param order: number of derivative orders to return; if `None`, all currently cached orders are returned
        :type order: int | None
        :param strip_embedding: whether to strip embedding coordinates from the trailing axes of the result
        :type strip_embedding: bool
        :return: list of Jacobian tensors, one per derivative order
        :rtype: list[np.ndarray]
        :raises ValueError: if fewer cached orders are available than requested
        """
        ...

    @property
    def internals_by_cartesians(self):
        """
        **LLM Docstring**

        All cached internals-by-Cartesians Jacobians, computing the default set if not already cached.

        :return: the `JacobianKeys.InternalsByCartesians` entry from `get_coordinate_transforms()`
        :rtype: list[np.ndarray]
        """
        ...

    def get_internals_by_cartesians(self, order=None, strip_embedding=False):
        """
        **LLM Docstring**

        Fetch the internals-by-Cartesians Jacobians up to the requested order, optionally stripping embedding coordinates from the trailing axis.

        :param order: number of derivative orders to return; if `None`, all currently cached orders are returned
        :type order: int | None
        :param strip_embedding: whether to strip embedding coordinates from the trailing axis of the result
        :type strip_embedding: bool
        :return: list of Jacobian tensors, one per derivative order
        :rtype: list[np.ndarray]
        :raises ValueError: if fewer cached orders are available than requested
        """
        ...

    @property
    def cartesian_modes_by_internal_modes(self):
        """
        **LLM Docstring**

        All cached Cartesian-normal-modes-by-internal-normal-modes Jacobians, computing the default set if not already cached.

        :return: the `JacobianKeys.CartesianModesByInternalModes` entry from `get_coordinate_transforms()`
        :rtype: list[np.ndarray]
        """
        ...

    def get_cartesian_modes_by_internal_modes(self, order=None):
        """
        **LLM Docstring**

        Fetch the Cartesian-normal-modes-by-internal-normal-modes Jacobians up to the requested order.

        :param order: number of derivative orders to return; if `None`, all currently cached orders are returned
        :type order: int | None
        :return: list of Jacobian tensors, one per derivative order
        :rtype: list[np.ndarray]
        :raises ValueError: if fewer cached orders are available than requested
        """
        ...

    @property
    def internal_modes_by_cartesian_modes(self):
        """
        **LLM Docstring**

        All cached internal-normal-modes-by-Cartesian-normal-modes Jacobians, computing the default set if not already cached.

        :return: the `JacobianKeys.InternalModesByCartesianModes` entry from `get_coordinate_transforms()`
        :rtype: list[np.ndarray]
        """
        ...

    def get_internal_modes_by_cartesian_modes(self, order=None):
        """
        **LLM Docstring**

        Fetch the internal-normal-modes-by-Cartesian-normal-modes Jacobians up to the requested order.

        :param order: number of derivative orders to return; if `None`, all currently cached orders are returned
        :type order: int | None
        :return: list of Jacobian tensors, one per derivative order
        :rtype: list[np.ndarray]
        :raises ValueError: if fewer cached orders are available than requested
        """
        ...

class PotentialTerms(ExpansionTerms):
    """
    A helper class that can transform the derivatives of the potential from Cartesian to normal coordinates
    """
    __props__ = ExpansionTerms.__props__ + ('potential_derivatives', 'check_input_force_constants', 'hessian_tolerance', 'grad_tolerance', 'freq_tolerance')

    def __init__(self, molecule, mixed_derivs=None, modes=None, potential_derivatives=None, mode_selection=None, mode_transformation=None, full_surface_mode_selection=None, logger=None, parallelizer=None, checkpointer=None, check_input_force_constants=True, allow_higher_potential_terms=False, hessian_tolerance=0.0001, grad_tolerance=0.0001, freq_tolerance=0.002, **opts):
        """
        :param molecule: the molecule that will supply the potential derivatives
        :type molecule: Molecule
        :param mixed_derivs: whether or not the pulled derivatives are partially derivatives along the normal coords
        :type mixed_derivs: bool
        :param modes: the normal modes to use when doing calculations
        :type modes: None | MolecularVibrations
        :param mode_selection: the subset of normal modes to use
        :type mode_selection: None | Iterable[int]
        """
        ...

    @property
    def v_derivs(self):
        """
        **LLM Docstring**

        Property getter/setter for the canonicalized potential-energy derivative tensors. The getter lazily pulls the raw derivatives from `self.molecule.potential_surface.derivatives` (if none were supplied at construction) and canonicalizes them via `_canonicalize_derivs`, caching the result.

        :param v: (setter only) the new canonicalized derivative tensors to store directly
        :type v: list[np.ndarray]
        :return: (getter) the canonicalized potential derivative tensors
        :rtype: list[np.ndarray]
        """
        ...

    @v_derivs.setter
    def v_derivs(self, v):
        """
        **LLM Docstring**

        Property getter/setter for the canonicalized potential-energy derivative tensors. The getter lazily pulls the raw derivatives from `self.molecule.potential_surface.derivatives` (if none were supplied at construction) and canonicalizes them via `_canonicalize_derivs`, caching the result.

        :param v: (setter only) the new canonicalized derivative tensors to store directly
        :type v: list[np.ndarray]
        :return: (getter) the canonicalized potential derivative tensors
        :rtype: list[np.ndarray]
        """
        ...

    def _check_mode_terms(self, derivs=None):
        """
        **LLM Docstring**

        Check whether a set of derivative tensors are already fully expressed in the normal-mode basis, i.e. every tensor axis has length equal to the number of modes.

        :param derivs: the derivative tensors to check; defaults to `self.v_derivs`
        :type derivs: list[np.ndarray] | None
        :return: whether every tensor's shape matches `(n_modes,) * ndim`
        :rtype: bool
        """
        ...

    def _canonicalize_derivs(self, freqs, masses, derivs, full_mode_sel, mode_transformation):
        """
        **LLM Docstring**

        Normalize the raw potential-energy derivative tensors (gradient, Hessian, and optionally cubic/quartic force constants, possibly bundled in various tuple/object forms) into a consistent list of plain arrays, validating each tensor's shape against the expected Cartesian/internal/mode dimensions and padding/selecting mode-subset axes as needed. Returns the derivatives unchanged if they're already fully mode-basis (via `_check_mode_terms`).

        :param freqs: the mode frequencies
        :type freqs: np.ndarray
        :param masses: the atomic masses
        :type masses: np.ndarray
        :param derivs: the raw derivative data: a `(grad, fcs, fds)` triple (where `fds` may bundle third/fourth derivatives), a `(grad, fcs, thirds, fourths)` quadruple, or an arbitrary-length sequence starting with gradient/Hessian
        :type derivs: tuple | list
        :param full_mode_sel: an index selection identifying which entries of a full mode set the supplied third/fourth derivatives correspond to, used to pad them back out to the full mode dimension
        :type full_mode_sel: Iterable[int] | None
        :param mode_transformation: accepted for interface consistency but not directly used to transform the derivatives in this method's body
        :type mode_transformation: object | None
        :return: the canonicalized list of derivative tensors
        :rtype: list[np.ndarray]
        :raises PerturbationTheoryException: if any derivative tensor's shape doesn't match the expected Cartesian, internal-coordinate, or mode-basis dimensions
        """
        ...

    @classmethod
    def _symmetrize_mixed_derivatives(cls, derivs, handling_mode, mode_axes, *, logger, zero_rest=True, diagonal=True, restricted_diagonal=False, term_id=None, val_axes=0, handle_zeros=True, warning_diff=-1):
        """
        **LLM Docstring**

        Enforce/repair the expected symmetry of a mixed-derivative tensor (one computed with some axes in a different coordinate basis than others), either applying special-cased treatments for known problematic term types (`'v4_cart'`/`'v4_int'`/`'u3_cart'`/`'u3_int'` under the legacy `Old` handling mode) or, in the general case, replacing every permutation of each set of symmetry-equivalent tensor indices with a single representative value (optionally warning when the differing raw values disagree by more than `warning_diff`, and optionally treating near-zero entries specially via `handle_zeros`).

        :param derivs: the (potentially asymmetric) mixed-derivative tensor to symmetrize
        :type derivs: np.ndarray
        :param handling_mode: which symmetrization strategy to use
        :type handling_mode: MixedDerivativeHandlingModes
        :param mode_axes: accepted for interface consistency but not used directly in this method's body
        :type mode_axes: object
        :param logger: logger used to report large symmetry-violating differences
        :type logger: Logger
        :param zero_rest: accepted for interface consistency but not used directly in this method's body
        :type zero_rest: bool
        :param diagonal: accepted for interface consistency but not used directly in this method's body
        :type diagonal: bool
        :param restricted_diagonal: accepted for interface consistency but not used directly in this method's body
        :type restricted_diagonal: bool
        :param term_id: an identifying label for the term being symmetrized, used both to select the legacy special-case handling and in warning messages
        :type term_id: str | None
        :param val_axes: number of leading "value" axes (not subject to the mode-index symmetrization) in `derivs`
        :type val_axes: int
        :param handle_zeros: whether near-zero raw entries should be treated as missing/overridable by a nonzero symmetry-equivalent value rather than triggering symmetry-violation warnings
        :type handle_zeros: bool
        :param warning_diff: minimum absolute difference between two nominally-equal symmetry-related values (in Hartrees) that triggers a logged warning; disabled if not positive
        :type warning_diff: float
        :return: the symmetrized derivative tensor
        :rtype: np.ndarray
        :raises ValueError: if `handling_mode` is `Old` and no `term_id` is given
        """
        ...

    def get_terms(self, order=None, logger=None):
        """
        **LLM Docstring**

        Compute the potential-energy expansion terms in the molecule's normal-mode coordinates, zeroing out the (should-be-vanishing) gradient term, re-expanding the Cartesian (or internal-coordinate) potential derivatives through the appropriate coordinate Jacobians (handling mixed-derivative-basis terms via `_symmetrize_mixed_derivatives` where relevant), and caching the result via the checkpointer.

        :param order: the highest derivative order to compute; if `None`, uses however many terms are available in `v_derivs`
        :type order: int | None
        :param logger: logger to report progress/timing/warnings to; defaults to `self.logger`
        :type logger: Logger | None
        :return: the potential-energy expansion terms, in mode-basis coordinates, from the (zeroed) gradient upward
        :rtype: list[np.ndarray]
        """
        ...

    @classmethod
    def get_potential_optimized_coordinates(cls, V_expansion, order=2):
        """
        **LLM Docstring**

        Find, order by order, the coordinate transformation that eliminates as much of the potential-energy expansion as possible (a "potential-optimized" coordinate system), by solving at each order for the new-coordinate term that cancels the corresponding remainder of the transformed potential.

        :param V_expansion: the potential-energy expansion terms (from the quadratic/Hessian term upward, i.e. without the zeroth-order energy)
        :type V_expansion: list[np.ndarray]
        :param order: the highest order to optimize the transformation to
        :type order: int
        :return: `(forward_derivs, reverse_derivs)` -- the forward and reverse coordinate-transformation derivative tensors
        :rtype: tuple[list[np.ndarray], list[np.ndarray]]
        """
        ...

    def optimize_coordinates(self, order=2):
        """
        **LLM Docstring**

        Build the potential-optimized coordinate transformation for this potential expansion and re-express it in terms of both the internal coordinates and the Cartesians, via `get_potential_optimized_coordinates` combined with the Cartesians-by-internals/internals-by-Cartesians Jacobians.

        :param order: the highest order to optimize the transformation to
        :type order: int
        :return: `((QR, RQ), (QX, XQ))` -- the forward/reverse transformations between the optimized coordinates and both the internal coordinates (`R`) and the Cartesians (`X`)
        :rtype: tuple[tuple[list[np.ndarray], list[np.ndarray]], tuple[list[np.ndarray], list[np.ndarray]]]
        """
        ...

class KineticTerms(ExpansionTerms):
    """Represents the KE coefficients"""
    __props__ = ExpansionTerms.__props__ + ('g_derivative_threshold', 'gmatrix_tolerance', 'use_cartesian_kinetic_energy', 'check_input_gmatrixfreq_tolerance')

    def __init__(self, molecule, g_derivative_threshold=0.001, gmatrix_tolerance=1e-06, use_cartesian_kinetic_energy=False, check_input_gmatrix=True, freq_tolerance=0.002, **opts):
        """
        **LLM Docstring**

        Set up a G-matrix (kinetic-energy coefficient) expansion generator for a molecule, storing the tolerances/thresholds used to validate and warn about the computed G-matrix and its derivatives.

        :param molecule: the molecule to compute the G-matrix expansion for
        :type molecule: Molecule
        :param g_derivative_threshold: the magnitude above which a G-matrix derivative term triggers a logged warning
        :type g_derivative_threshold: float
        :param gmatrix_tolerance: the tolerance used when checking that the zeroth-order G-matrix is diagonal
        :type gmatrix_tolerance: float
        :param use_cartesian_kinetic_energy: whether to force computing the kinetic energy directly in Cartesian coordinates rather than internal coordinates
        :type use_cartesian_kinetic_energy: bool
        :param check_input_gmatrix: whether to validate the reconstructed G-matrix frequencies against the nominal mode frequencies
        :type check_input_gmatrix: bool
        :param freq_tolerance: the tolerance used for that frequency validation
        :type freq_tolerance: float
        :param opts: extra options forwarded to the base `ExpansionTerms.__init__`
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...

    def get_terms(self, order=None, logger=None, return_expressions=False):
        """
        **LLM Docstring**

        Compute the G-matrix (kinetic-energy coefficient) Taylor-expansion terms in the molecule's normal-mode coordinates, either via a simplified direct Cartesian-mode-matrix contraction (when working purely in Cartesians) or by chaining together the internals/Cartesians-by-modes Jacobians via `TensorDerivativeConverter` and differentiating iteratively, validating that the zeroth-order term is diagonal and matches the nominal mode frequencies, and warning about anomalously large higher-order derivatives.

        :param order: the highest derivative order to compute
        :type order: int | None
        :param logger: logger to report progress/timing/warnings to; defaults to `self.logger`
        :type logger: Logger | None
        :param return_expressions: whether to also return the underlying symbolic `TensorExpansionTerms` objects alongside the numeric arrays
        :type return_expressions: bool
        :return: the G-matrix expansion terms (or `(terms, expressions)` if `return_expressions` is set)
        :rtype: list[np.ndarray] | tuple
        :raises ValueError: if the zeroth-order G-matrix isn't (sufficiently) diagonal
        :raises PerturbationTheoryException: if the frequencies implied by the G-matrix don't match the nominal mode frequencies within `freq_tolerance`
        """
        ...

    @classmethod
    def _dRGQ_partition_contrib(cls, partition, R, G):
        """
        **LLM Docstring**

        Compute one term of the derivative of the reverse-transformed G-matrix with respect to the new coordinates, for a given `(r1, r2, s)` integer partition -- i.e. the contribution from taking `r1` and `r2` derivatives of the reverse transformation and `s` derivatives of the original G-matrix.

        :param partition: the `(r1, r2, s)` triple of derivative counts for this term
        :type partition: tuple[int, int, int]
        :param R: the reverse-transformation derivative tensors, indexed by order
        :type R: list[np.ndarray]
        :param G: the original G-matrix expansion terms, indexed by order
        :type G: list[np.ndarray]
        :return: the contribution of this partition to the derivative, or `0` if any required term is out of range or exactly zero
        :rtype: np.ndarray | int
        """
        ...

    @classmethod
    def _dRGQ_derivs(cls, R, G, o):
        """
        **LLM Docstring**

        Sum the contributions of every valid `(r1, r2, g)` integer partition of order `o+1` (restricted to `r1 >= rem//2` by symmetry) to get the full order-`o` derivative of the reverse-transformed G-matrix.

        :param R: the reverse-transformation derivative tensors, indexed by order
        :type R: list[np.ndarray]
        :param G: the original G-matrix expansion terms, indexed by order
        :type G: list[np.ndarray]
        :param o: the derivative order to compute
        :type o: int
        :return: the summed contribution over all valid partitions
        :rtype: np.ndarray | int
        """
        ...

    @classmethod
    def reexpress_G(self, G_expansion, forward_derivs, reverse_derivs=None, order=2):
        """
        Apply a coordinate transformation to the G-matrix

        :param forward_derivs:
        :param reverse_derivs:
        :param order:
        :return:
        """
        ...

    def reexpress(self, forward_derivs, reverse_derivs=None, order=2):
        """
        Finds a coordinate transformation the give 0 contribution to the G-matrix

        :param forward_derivs:
        :param reverse_derivs:
        :param order:
        :return:
        """
        ...

    @classmethod
    def get_kinetic_optimized_coordinates(cls, G_expansion, order=2):
        """
        **LLM Docstring**

        Intended to iteratively find the coordinate transformation that eliminates as much of the G-matrix expansion as possible order by order. As written, after computing the first-order correction `R2` and the resulting re-expressed G-matrix, the method unconditionally executes `raise Exception(new_G[1])` -- it never returns normally, so it does not currently function, and the remaining (otherwise complete-looking) loop below that line is unreachable dead code.

        :param G_expansion: the G-matrix expansion terms to optimize
        :type G_expansion: list[np.ndarray]
        :param order: the target expansion order to optimize up to
        :type order: int
        :return: never returns normally; always raises
        :rtype: None
        :raises Exception: unconditionally, carrying the (successfully computed) first-order re-expressed G-matrix derivative as its argument
        """
        ...

    def optimize_coordinates(self, order=2):
        """
        **LLM Docstring**

        Build the kinetic-energy-optimized coordinate transformation for this G-matrix expansion and re-express it in terms of both the internal coordinates and the Cartesians. Note this method calls `get_kinetic_optimized_coordinates`, which as currently implemented always raises rather than returning a transformation (see that method's docstring), so calling this method will likewise fail.

        :param order: the highest order to optimize the transformation to
        :type order: int
        :return: intended to be `((QR, RQ), (QX, XQ))`, the forward/reverse transformations between the optimized coordinates and both the internal coordinates and the Cartesians; never actually returned due to the exception raised by `get_kinetic_optimized_coordinates`
        :rtype: tuple
        :raises Exception: propagated from `get_kinetic_optimized_coordinates`
        """
        ...

class CoriolisTerm(ExpansionTerms):
    """
    Calculates the Coriolis coupling term
    """

    def get_zetas_and_momi(self):
        """
        **LLM Docstring**

        Compute the Coriolis zeta constants and the inertial-frame moments of inertia: mass-weights and frequency-descales the mode matrix, reshapes it into an atom/mode/Cartesian-axis tensor, rotates it into the inertial (principal-axis) frame, and forms the antisymmetric (Levi-Civita) combination across atoms that gives the zeta constants.

        :return: `(zeta, B_e)` -- the zeta-constant tensor (mode x mode x 3 x 3) and the inertial-frame moments of inertia
        :rtype: tuple[np.ndarray, np.ndarray]
        """
        ...

    def get_zetas(self):
        """
        **LLM Docstring**

        The Coriolis zeta-constant tensor alone, via `get_zetas_and_momi`.

        :return: the zeta-constant tensor
        :rtype: np.ndarray
        """
        ...

    def get_terms(self, order=None, J=0):
        """
        **LLM Docstring**

        Compute the Coriolis rotational-coupling operator's Taylor-expansion terms by combining the frequency-dimensioned zeta constants with the expansion of the reciprocal moment-of-inertia tensor, adding one extra coordinate axis per order.

        :param order: the highest derivative order to compute
        :type order: int | None
        :param J: the total rotational angular momentum quantum number; only `J=0` (pure vibration-rotation coupling with no external rotation) is currently supported
        :type J: int
        :return: the Coriolis expansion terms, one per order
        :rtype: list[np.ndarray]
        :raises NotImplementedError: if `J > 0`
        """
        ...

class PotentialLikeTerm(KineticTerms):
    """
    This accounts for the potential-like term.
    In Cartesian diplacement modes this is the Watson U.
    In proper internals, this is the V' term.
    """

    def get_terms(self, order=None, logger=None):
        """
        **LLM Docstring**

        Compute the potential-like (Watson `U` / internal-coordinate `V'`) correction term's Taylor-expansion terms: either directly from the trace of the reciprocal-inertia-tensor derivatives (when working purely in Cartesian modes), or, for internal coordinates, via the standard Watson pseudopotential derivation combining the G-matrix and inertia-tensor log-determinant derivatives (`d/dQ[ln(detI) - ln(detG)]`, using a matrix-cookbook identity for the log-determinant derivative).

        :param order: the highest derivative order to compute
        :type order: int | None
        :param logger: logger to report progress to; only used when computing the G-matrix terms via the base class
        :type logger: Logger | None
        :return: the potential-like-term expansion, one per order
        :rtype: list[np.ndarray]
        """
        ...

class DipoleTerms(ExpansionTerms):
    __props__ = ExpansionTerms.__props__ + ('dipole_derivatives',)

    def __init__(self, molecule, dipole_derivatives=None, mixed_derivs=None, modes=None, mode_selection=None, mode_transformation=None, full_surface_mode_selection=None, logger=None, parallelizer=None, checkpointer=None, **opts):
        """
        :param molecule: the molecule that will supply the dipole derivatives
        :type molecule: Molecule
        :param mixed_derivs: whether or not the pulled derivatives are partially derivatives along the normal coords
        :type mixed_derivs: bool
        :param modes: the normal modes to use when doing calculations
        :type modes: None | MolecularVibrations
        :param mode_selection: the subset of normal modes to use
        :type mode_selection: None | Iterable[int]
        """
        ...

    def _canonicalize_derivs(self, freqs, masses, derivs, full_mode_sel, mode_transformation):
        """
        Makes sure all of the dipole moments are clean and ready to rotate
        """
        ...

    def _check_mode_terms(self, derivs=None):
        """
        **LLM Docstring**

        Check whether a set of dipole derivative tensors are already fully expressed in the normal-mode basis, i.e. every tensor axis (aside from the trailing Cartesian-component axis where applicable) has length equal to the number of modes.

        :param derivs: the derivative tensors to check; defaults to `self.derivs[1:]`
        :type derivs: list[np.ndarray] | None
        :return: whether every tensor's shape matches `(n_modes,) * ndim`
        :rtype: bool
        """
        ...

    def get_terms(self, order=None, logger=None):
        """
        **LLM Docstring**

        Compute the dipole-moment Taylor-expansion terms (one component at a time) in the molecule's normal-mode coordinates, re-expanding the Cartesian (or internal-coordinate) dipole derivatives through the appropriate coordinate Jacobians and handling mixed-derivative-basis terms where relevant.

        :param order: the highest derivative order to compute; if `None`, uses however many terms are available in `self.derivs`
        :type order: int | None
        :param logger: logger to report progress/timing to; defaults to `self.logger`
        :type logger: Logger | None
        :return: the per-Cartesian-component dipole expansion terms
        :rtype: list
        """
        ...

class OperatorTerms(ExpansionTerms):
    """
    Literally as simple as it comes for an operator expansion.
    One dimensional, no mixed derivative stuff.
    """
    __props__ = ExpansionTerms.__props__ + ('operator_derivatives',)

    def __init__(self, molecule, operator_derivatives=None, modes=None, mode_selection=None, logger=None, parallelizer=None, checkpointer=None, **opts):
        """
        :param molecule: the molecule that will supply the dipole derivatives
        :type molecule: Molecule
        :param mixed_derivs: whether or not the pulled derivatives are partially derivatives along the normal coords
        :type mixed_derivs: bool
        :param modes: the normal modes to use when doing calculations
        :type modes: None | MolecularVibrations
        :param mode_selection: the subset of normal modes to use
        :type mode_selection: None | Iterable[int]
        """
        ...

    def _check_mode_terms(self, derivs=None):
        """
        **LLM Docstring**

        Check whether a set of operator derivative tensors are already fully expressed in the normal-mode basis, i.e. every tensor axis has length equal to the number of modes.

        :param derivs: the derivative tensors to check; defaults to `self.derivs`
        :type derivs: list[np.ndarray] | None
        :return: whether every tensor's shape matches `(n_modes,) * ndim`
        :rtype: bool
        """
        ...

    def _canonicalize_derivs(self, freqs, masses, derivs):
        """
        **LLM Docstring**

        Normalize a set of raw operator derivative tensors into a consistent list of plain (mass/frequency-undimensioned) arrays, validating each tensor's shape against the expected Cartesian or internal-coordinate dimension and applying the appropriate mass- or frequency-based unit conversion. Returns the derivatives unchanged if they're already fully mode-basis.

        :param freqs: the mode frequencies
        :type freqs: np.ndarray
        :param masses: the atomic masses
        :type masses: np.ndarray
        :param derivs: the raw operator derivative tensors, one per order
        :type derivs: list[np.ndarray]
        :return: the canonicalized (unit-converted) list of derivative tensors
        :rtype: list[np.ndarray]
        :raises PerturbationTheoryException: if any derivative tensor's shape doesn't match the expected Cartesian or internal-coordinate dimensions
        """
        ...

    def get_terms(self, order=None, logger=None):
        """
        **LLM Docstring**

        Compute the operator's Taylor-expansion terms in the molecule's normal-mode coordinates, returning the already-canonicalized `self.derivs` directly if they're already mode-basis.

        :param order: the highest derivative order to compute (only relevant when re-expansion through coordinate Jacobians is needed)
        :type order: int | None
        :param logger: logger to report progress to; defaults to `self.logger`
        :type logger: Logger | None
        :return: the operator expansion terms
        :rtype: list[np.ndarray]
        """
        ...