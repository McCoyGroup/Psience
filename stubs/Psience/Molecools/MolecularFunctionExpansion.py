import enum, numpy as np, functools as fp
from McUtils.Scaffolding import Logger, NullLogger, Checkpointer, NullCheckpointer
from McUtils.Parallelizers import Parallelizer
from McUtils.Zachary import TensorDerivativeConverter, TensorExpansionTerms

class JacobianKeys(enum.Enum):
    """Real access pattern: JacobianKeys.<MemberName> (this is an enum with 12 members, e.g. JacobianKeys.CartesiansByInternals == 'CartesiansByInternals'). Collapsed into a dict below purely for compactness -- do not index it as a dict in real code:"""
    _MEMBERS = {'CartesiansByInternals': 'CartesiansByInternals', 'InternalsByCartesians': 'InternalsByCartesians', 'InternalsByCartesianModes': 'InternalsByModes', 'CartesianModesByInternals': 'ModesByInternals', 'CartesiansByInternalModes': 'CartesiansByModes', 'InternalModesByCartesians': 'ModesByCartesians', 'CartesianModesByInternalModes': 'CartesianModesByInternalModes', 'InternalModesByCartesianModes': 'InternalModesByCartesianModes', 'InternalModesByInternals': 'InternalModesByInternals', 'InternalsByInternalModes': 'InternalsByInternalModes', 'CartesianModesByCartesians': 'CartesianModesByCartesians', 'CartesiansByCartesianModes': 'CartesiansByCartesianModes'}

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

class MolecularFunctionExpander:
    """
    Handles all necessary details to be able to
    evaluate functions in molecular internal coordinates
    and then transform back to Cartesians, and to do the
    same with respect to normal dimensionless normal mode coordinates
    """

    def __init__(self, cartesian_coordinates, internal_coordinates, transform=None, masses=None, freqs=None, parallelizer=None, logger=None, checkpointer=None, numerical_jacobians=True, eckart_embed_derivatives=True, eckart_embed_planar_ref_tolerance=None, strip_dummies=False, strip_embedding=True, internal_fd_mesh_spacing=0.001, internal_fd_stencil=None, cartesian_fd_mesh_spacing=0.01, cartesian_fd_stencil=None, cartesian_analytic_deriv_order=0, internal_by_cartesian_order=3, cartesian_by_internal_order=4, jacobian_warning_threshold=10000.0):
        """
        **LLM Docstring**

        Set up the machinery for evaluating a function in internal coordinates and converting its derivatives back to Cartesian and dimensionless normal-mode coordinates, including finite-difference settings, the mode transformation, logging, parallelization, and checkpointing.

        :param cartesian_coordinates: the Cartesian coordinates of the reference geometry
        :type cartesian_coordinates: np.ndarray
        :param internal_coordinates: the internal coordinates of the reference geometry
        :type internal_coordinates: CoordinateSet
        :param transform: normal-mode transformation to use; if given, is un-dimensionalized via `self.undimensionalize` to populate `self.modes`/`self.inverse`
        :type transform: object | None
        :param masses: atomic masses; defaults to an array of ones sized to the number of atoms
        :type masses: np.ndarray | None
        :param freqs: mode frequencies; defaults to an array of ones sized to the number of atoms
        :type freqs: np.ndarray | None
        :param parallelizer: parallelization backend to use for the Jacobian calculations
        :type parallelizer: Parallelizer | None
        :param logger: logger for warnings/diagnostics; defaults to a `NullLogger`
        :type logger: Logger | None
        :param checkpointer: checkpointing backend for caching results; defaults to a `NullCheckpointer`
        :type checkpointer: Checkpointer | None
        :param numerical_jacobians: whether all Jacobians should be computed numerically rather than analytically
        :type numerical_jacobians: bool
        :param eckart_embed_derivatives: whether derivatives should be re-embedded in the Eckart frame
        :type eckart_embed_derivatives: bool
        :param eckart_embed_planar_ref_tolerance: tolerance used when Eckart-embedding a (near-)planar reference structure
        :type eckart_embed_planar_ref_tolerance: float | None
        :param strip_dummies: whether dummy atoms should be excluded from the mass array
        :type strip_dummies: bool
        :param strip_embedding: whether embedding coordinates should be stripped out of the Jacobians by default
        :type strip_embedding: bool
        :param internal_fd_mesh_spacing: finite-difference step size used for internal-coordinate derivatives
        :type internal_fd_mesh_spacing: float
        :param internal_fd_stencil: finite-difference stencil size for internal-coordinate derivatives
        :type internal_fd_stencil: int | None
        :param cartesian_fd_mesh_spacing: finite-difference step size used for Cartesian-coordinate derivatives
        :type cartesian_fd_mesh_spacing: float
        :param cartesian_fd_stencil: finite-difference stencil size for Cartesian-coordinate derivatives
        :type cartesian_fd_stencil: int | None
        :param cartesian_analytic_deriv_order: order up to which Cartesian derivatives are computed analytically instead of by finite difference
        :type cartesian_analytic_deriv_order: int
        :param internal_by_cartesian_order: default order of internal-by-Cartesian Jacobians to compute
        :type internal_by_cartesian_order: int
        :param cartesian_by_internal_order: default order of Cartesian-by-internal Jacobians to compute
        :type cartesian_by_internal_order: int
        :param jacobian_warning_threshold: magnitude above which Jacobian entries trigger a logged warning and are zeroed out
        :type jacobian_warning_threshold: float
        :return: None
        :rtype: None
        """
        ...

    def _tripmass(self, masses):
        """
        **LLM Docstring**

        Mass-weighting helper: optionally drops dummy-atom masses (zero or negative entries, when `self.strip_dummies` is set), then repeats each mass three times (once per Cartesian direction) and flattens to a length `3*n_atoms` vector.

        :param masses: per-atom masses
        :type masses: np.ndarray
        :return: the per-Cartesian-coordinate mass vector
        :rtype: np.ndarray
        """
        ...

    def undimensionalize(self, masses, freqs, modes):
        """
        Removes units from normal modes

        :param masses:
        :type masses:
        :param modes:
        :type modes:
        :return:
        :rtype:
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

    def _get_embedding_coords(self):
        """
        **LLM Docstring**

        Look up the indices of the coordinates used purely for embedding (e.g. translation/rotation) from the internal coordinate system, trying `system.embedding_coords` first and falling back to `system.converter_options['embedding_coords']`.

        :return: the embedding-coordinate indices, or `None` if neither source defines them
        :rtype: np.ndarray | None
        """
        ...

    def get_coordinate_transforms(self, internal_by_cartesian_order=None, cartesian_by_internal_order=None):
        """
        **LLM Docstring**

        Compute (and cache in `self._cached_transforms`) the full set of Jacobians relating Cartesian coordinates, internal coordinates, Cartesian normal modes, and internal-coordinate-basis normal modes to each other, up to the requested derivative orders. Handles mass-weighting, stripping of embedding coordinates, warns about (and zeroes) anomalously large Jacobian entries, and chains the individual Jacobians together via `TensorDerivativeConverter` to populate every entry of `JacobianKeys`.

        :param internal_by_cartesian_order: derivative order (number of Cartesian derivatives) to compute for internals-by-Cartesians Jacobians; defaults to `self.internal_by_cartesian_order`
        :type internal_by_cartesian_order: int | None
        :param cartesian_by_internal_order: derivative order (number of internal derivatives) to compute for Cartesians-by-internals Jacobians; defaults to `self.cartesian_by_internal_order`
        :type cartesian_by_internal_order: int | None
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

        Fetch the internal-normal-modes-by-internals Jacobians up to the requested order, optionally stripping embedding coordinates out of the leading axis.

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

        Fetch the internals-by-internal-normal-modes Jacobians up to the requested order, optionally stripping embedding coordinates out of the trailing axis.

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

    def evaluate(self, function, deriv_order=2):
        """
        **LLM Docstring**

        Stub for evaluating `function` and its derivatives up to `deriv_order` in the transformed coordinate systems handled by this class. The body is just `...`, so this method is currently unimplemented and does nothing.

        :param function: the function to evaluate
        :type function: callable
        :param deriv_order: the highest derivative order to evaluate
        :type deriv_order: int
        :return: not implemented
        :rtype: None
        """
        ...