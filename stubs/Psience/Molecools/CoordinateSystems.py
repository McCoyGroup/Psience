"""
Defines useful extended internal coordinate frames
"""
import numpy as np
import McUtils.Devutils as dev
import McUtils.Numputils as nput
from McUtils.Parallelizers import Parallelizer
import McUtils.Coordinerds
from McUtils.Coordinerds import ZMatrixCoordinateSystem, CartesianCoordinateSystem, CoordinateSystemConverter, ZMatrixCoordinates, CartesianCoordinates3D, CoordinateSet, CompositeCoordinateSystem, GenericInternalCoordinateSystem, GenericInternalCoordinates, CartesianToGICSystemConverter, GICSystemToCartesianConverter, PrimitiveCoordinatePicker, RedundantCoordinateGenerator
from ..Modes import NormalModes
from .Properties import StructuralProperties
__all__ = ['MolecularEmbedding', 'ModeEmbedding', 'MolecularZMatrixCoordinateSystem', 'MolecularCartesianCoordinateSystem']
__reload_hook__ = ['.MoleculeInterface']

class MolecularEmbedding:

    def __init__(self, masses, coords, internals, internal_fd_opts=None, cartesian_fd_opts=None):
        """
        **LLM Docstring**

        Set up a molecule's coordinate embedding: wraps the Cartesian coordinates in a `MolecularCartesianCoordinateSystem`, canonicalizes the internal-coordinate specification (or stores it directly if already a `CoordinateSet`), and initializes the Jacobian cache and finite-difference option overrides.

        :param masses: the atomic masses
        :type masses: np.ndarray
        :param coords: the Cartesian coordinates
        :type coords: np.ndarray
        :param internals: the internal-coordinate specification (Z-matrix, generic-internal specs, a callable conversion, or an already-built `CoordinateSet`)
        :type internals: object
        :param internal_fd_opts: overrides for the internal-coordinate finite-difference defaults
        :type internal_fd_opts: dict | None
        :param cartesian_fd_opts: overrides for the Cartesian finite-difference defaults
        :type cartesian_fd_opts: dict | None
        :return: None
        :rtype: None
        """
        ...

    def get_direct_converter(self, target):
        """
        **LLM Docstring**

        Provide a converter from this embedding's coordinate system directly to plain (non-molecular) 3D Cartesian coordinates, if `target` is Cartesian-compatible.

        :param target: the coordinate system being converted to
        :type target: object
        :return: a `MolecularCartesianToRegularCartesianConverter`, or `None` if `target` isn't Cartesian-compatible
        :rtype: MolecularCartesianToRegularCartesianConverter | None
        """
        ...

    def get_inverse_converter(self, target):
        """
        **LLM Docstring**

        Provide a converter from plain (non-molecular) 3D Cartesian coordinates into this embedding's coordinate system, if `target` is Cartesian-compatible.

        :param target: the coordinate system being converted from
        :type target: object
        :return: a `RegularCartesianToMolecularCartesianConverter`, or `None` if `target` isn't Cartesian-compatible
        :rtype: RegularCartesianToMolecularCartesianConverter | None
        """
        ...

    def __del__(self):
        """
        **LLM Docstring**

        Deregister any coordinate converters this embedding registered, via `cleanup`, when the object is garbage collected.

        :return: None
        :rtype: None
        """
        ...

    def cleanup(self):
        """
        **LLM Docstring**

        Deregister every converter previously registered via `register`, if any.

        :return: None
        :rtype: None
        """
        ...

    def register(self):
        """
        **LLM Docstring**

        Register the Cartesian coordinate converters for this embedding's coordinate system with the global converter registry, if not already registered.

        :return: None
        :rtype: None
        """
        ...

    @property
    def coords(self):
        """
        **LLM Docstring**

        Property getter/setter for the Cartesian coordinates. The getter registers the coordinate converters (via `register`) before returning them. The setter accepts a raw array or an already-systemed `CoordinateSet`, invalidates the Jacobian cache and inertial-frame cache, and marks the converters as needing re-registration if the coordinate system changed.

        :param coords: (setter only) the new Cartesian coordinates
        :type coords: np.ndarray | CoordinateSet
        :return: (getter) the Cartesian coordinates
        :rtype: CoordinateSet
        """
        ...

    @coords.setter
    def coords(self, coords):
        """
        **LLM Docstring**

        Property getter/setter for the Cartesian coordinates. The getter registers the coordinate converters (via `register`) before returning them. The setter accepts a raw array or an already-systemed `CoordinateSet`, invalidates the Jacobian cache and inertial-frame cache, and marks the converters as needing re-registration if the coordinate system changed.

        :param coords: (setter only) the new Cartesian coordinates
        :type coords: np.ndarray | CoordinateSet
        :return: (getter) the Cartesian coordinates
        :rtype: CoordinateSet
        """
        ...

    @property
    def masses(self):
        """
        **LLM Docstring**

        The atomic masses, taken from the Cartesian coordinate system.

        :return: the atomic masses
        :rtype: np.ndarray
        """
        ...

    @staticmethod
    def _wrap_conv(f):
        """
        **LLM Docstring**

        Wrap a user-supplied coordinate-conversion function so it always returns a uniform `(values, opts)` pair: if the function returns a bare array, an empty options dict is paired with it; if it returns `(values, opts)`, the caller's original keyword arguments are merged underneath the returned `opts`.

        :param f: the conversion function to wrap, or `None`
        :type f: callable | None
        :return: the wrapped conversion function, or `None` if `f` is `None`
        :rtype: callable | None
        """
        ...

    @classmethod
    def canonicalize_internal_coordinate_spec(cls, spec):
        """
        **LLM Docstring**

        Normalize the many accepted forms of an internal-coordinate specification (an options dict with `'zmatrix'`/`'specs'`/`'conversion'` keys, a bare callable conversion function, a raw Z-matrix-like array, or a list of generic-internal-coordinate specs) into the single canonical dict form (`'specs'`, `'zmatrix'`, `'conversion'`, `'inverse'`, `'converter_options'`) used internally, wrapping any conversion callables via `_wrap_conv` and filling in default embedding/jacobian-prep converter options for Z-matrices.

        :param spec: the internal-coordinate specification to canonicalize
        :type spec: dict | callable | Iterable | None
        :return: the canonicalized specification dict, or `None`/the original value if `spec` is `None`
        :rtype: dict | None
        :raises ValueError: if a Z-matrix-like spec doesn't have exactly 4 columns
        """
        ...

    @property
    def internals(self):
        """
        **LLM Docstring**

        Property getter/setter for the raw (canonicalized) internal-coordinate specification. The setter re-canonicalizes the given specification and invalidates any already-computed internal coordinates.

        :param internals: (setter only) the new internal-coordinate specification, in any form accepted by `canonicalize_internal_coordinate_spec`
        :type internals: object
        :return: (getter) the canonicalized specification dict, or `None` if none is set
        :rtype: dict | None
        """
        ...

    @internals.setter
    def internals(self, internals):
        """
        **LLM Docstring**

        Property getter/setter for the raw (canonicalized) internal-coordinate specification. The setter re-canonicalizes the given specification and invalidates any already-computed internal coordinates.

        :param internals: (setter only) the new internal-coordinate specification, in any form accepted by `canonicalize_internal_coordinate_spec`
        :type internals: object
        :return: (getter) the canonicalized specification dict, or `None` if none is set
        :rtype: dict | None
        """
        ...

    @property
    def zmatrix(self):
        """
        **LLM Docstring**

        Property getter/setter for just the Z-matrix ordering array out of the internal-coordinate specification. The setter validates the Z-matrix shape, builds a fresh specification if none exists yet, and invalidates any already-computed internal coordinates.

        :param zmat: (setter only) the new Z-matrix ordering array
        :type zmat: np.ndarray | None
        :return: (getter) the stored Z-matrix array, or `None`
        :rtype: np.ndarray | None
        :raises ValueError: if `zmat` doesn't have exactly 4 columns
        """
        ...

    @zmatrix.setter
    def zmatrix(self, zmat):
        """
        **LLM Docstring**

        Property getter/setter for just the Z-matrix ordering array out of the internal-coordinate specification. The setter validates the Z-matrix shape, builds a fresh specification if none exists yet, and invalidates any already-computed internal coordinates.

        :param zmat: (setter only) the new Z-matrix ordering array
        :type zmat: np.ndarray | None
        :return: (getter) the stored Z-matrix array, or `None`
        :rtype: np.ndarray | None
        :raises ValueError: if `zmat` doesn't have exactly 4 columns
        """
        ...

    @classmethod
    def convert_to_internals(cls, coords, masses, spec):
        """
        **LLM Docstring**

        Build the internal-coordinate `CoordinateSet` described by `spec`: constructs (and registers converters for) a generic-internal, Z-matrix, or iterative-Z-matrix coordinate system as appropriate, converts `coords` into it, layers on any extra custom `conversion`/`inverse` via a `CompositeCoordinateSystem`, and returns the resulting coordinates together with the (possibly updated, e.g. with redundant-transformation info) spec.

        :param coords: the Cartesian coordinates to convert
        :type coords: CoordinateSet
        :param masses: the atomic masses
        :type masses: np.ndarray
        :param spec: the canonicalized internal-coordinate specification (as produced by `canonicalize_internal_coordinate_spec`)
        :type spec: dict
        :return: `(coords, spec)` -- the internal coordinates and the (possibly updated) specification
        :rtype: tuple[CoordinateSet, dict]
        """
        ...

    def internal_coordinates_from_spec(self, spec: dict):
        """
        **LLM Docstring**

        Build the internal coordinates for this embedding's current Cartesian coordinates and masses from a given specification, via `convert_to_internals`.

        :param spec: the canonicalized internal-coordinate specification
        :type spec: dict
        :return: `(coords, spec)`, as returned by `convert_to_internals`
        :rtype: tuple[CoordinateSet, dict]
        """
        ...

    @property
    def internal_coordinates(self):
        """
        **LLM Docstring**

        Property getter/setter for the internal coordinates. The getter lazily computes them from the stored specification (via `internal_coordinates_from_spec`) the first time they're needed. The setter requires an already-built `CoordinateSet`.

        :param ics: (setter only) the new internal coordinates
        :type ics: CoordinateSet
        :return: (getter) the internal coordinates, or `None` if no internal-coordinate specification is set
        :rtype: CoordinateSet | None
        :raises ValueError: if the setter is given something that isn't a `CoordinateSet`
        """
        ...

    @internal_coordinates.setter
    def internal_coordinates(self, ics):
        """
        **LLM Docstring**

        Property getter/setter for the internal coordinates. The getter lazily computes them from the stored specification (via `internal_coordinates_from_spec`) the first time they're needed. The setter requires an already-built `CoordinateSet`.

        :param ics: (setter only) the new internal coordinates
        :type ics: CoordinateSet
        :return: (getter) the internal coordinates, or `None` if no internal-coordinate specification is set
        :rtype: CoordinateSet | None
        :raises ValueError: if the setter is given something that isn't a `CoordinateSet`
        """
        ...

    def strip_embedding_coordinates(self, coords):
        """
        **LLM Docstring**

        Drop the fixed embedding coordinates (e.g. the 6 translation/rotation degrees of freedom implied by the Z-matrix embedding) from a coordinate array or list of derivative tensors, if the underlying internal-coordinate system defines any.

        :param coords: the coordinates (or list of derivative tensors) to strip
        :type coords: np.ndarray | list[np.ndarray]
        :return: the coordinates with embedding coordinates removed, or unchanged if there are none to strip or they're already stripped
        :rtype: np.ndarray | list[np.ndarray]
        """
        ...

    def strip_derivative_embedding(self, derivs):
        """
        **LLM Docstring**

        Drop the fixed embedding coordinates from every axis of each tensor in a list of Cartesian-derivative tensors, if the underlying internal-coordinate system defines any.

        :param derivs: the list of Cartesian-derivative tensors (order-`n` tensor at index `n-1`) to strip
        :type derivs: list[np.ndarray]
        :return: the derivative tensors with embedding coordinates removed from every relevant axis, or unchanged if there are none to strip or they're already stripped
        :rtype: list[np.ndarray]
        """
        ...

    def restore_embedding_coordinates(self, coords):
        """
        **LLM Docstring**

        Reinsert the fixed embedding coordinates (filled in from the reference internal coordinates) back into a stripped coordinate array or list, undoing `strip_embedding_coordinates`.

        :param coords: the stripped coordinates (or list) to restore
        :type coords: np.ndarray | list[np.ndarray]
        :return: the coordinates with embedding coordinates reinserted, or unchanged if there are none to restore or they're already present
        :rtype: np.ndarray | list[np.ndarray]
        """
        ...

    def restore_derivative_embedding(self, derivs):
        """
        **LLM Docstring**

        Reinsert zeroed-out placeholder entries for the fixed embedding coordinates back into every axis of each tensor in a list of stripped Cartesian-derivative tensors, undoing `strip_derivative_embedding`.

        :param derivs: the stripped list of Cartesian-derivative tensors to restore
        :type derivs: list[np.ndarray]
        :return: the derivative tensors with embedding-coordinate axes reinserted (as zeros), or unchanged if there are none to restore or they're already present
        :rtype: list[np.ndarray]
        """
        ...

    def get_internals(self, *, coords=None, strip_embedding=True):
        """
        **LLM Docstring**

        Fetch the internal coordinates, either the cached ones for this embedding's own geometry or freshly computed ones for an alternate set of Cartesian `coords`, optionally stripping the fixed embedding coordinates.

        :param coords: alternate Cartesian coordinates to convert instead of using the cached internal coordinates
        :type coords: np.ndarray | None
        :param strip_embedding: whether to strip the fixed embedding coordinates from the result
        :type strip_embedding: bool
        :return: the internal coordinates, or `None` if none are defined
        :rtype: CoordinateSet | None
        """
        ...

    def get_cartesians(self, *, coords=None, strip_embedding=True):
        """
        **LLM Docstring**

        Fetch Cartesian coordinates, either this embedding's own cached ones or the Cartesian coordinates corresponding to a given set of internal `coords`, optionally restoring any stripped embedding coordinates first.

        :param coords: internal-coordinate values to convert to Cartesians instead of returning the cached Cartesian coordinates
        :type coords: np.ndarray | None
        :param strip_embedding: whether `coords` has had its embedding coordinates stripped and needs them restored before conversion
        :type strip_embedding: bool
        :return: the Cartesian coordinates
        :rtype: CoordinateSet
        """
        ...

    @property
    def redundant_internal_transformation(self):
        """
        **LLM Docstring**

        The redundant-to-non-redundant transformation matrix associated with the current internal coordinates, if the internal-coordinate system used a redundant coordinate generator.

        :return: the redundant transformation, or `None` if not applicable
        :rtype: np.ndarray | None
        """
        ...

    @classmethod
    def _get_jacobian_storage(cls):
        """
        **LLM Docstring**

        Build a fresh, empty nested dict used to cache computed Jacobian tensors by coordinate-system pair and (for internals) by re-embedding mode.

        :return: the empty Jacobian cache structure
        :rtype: dict
        """
        ...

    def _get_internal_fd_opts(self, **opts):
        """
        **LLM Docstring**

        Merge the class-level `internal_fd_defaults`, this instance's `_internal_fd_opts` overrides, and any explicitly passed `opts`, with later sources taking precedence, to resolve the finite-difference options used for internal-coordinate Jacobians.

        :param opts: explicit per-call overrides, taking highest precedence
        :type opts: dict
        :return: the merged finite-difference options
        :rtype: dict
        """
        ...

    def _get_int_jacobs(self, jacs, coords=None, strip_embedding=False, **fd_opts):
        """
        Gets the specified dX/dRs

        :param jacs:
        :type jacs:
        :return:
        :rtype:
        """
        ...

    def _get_cart_fd_opts(self, **opts):
        """
        **LLM Docstring**

        Merge the class-level `cart_fd_defaults`, this instance's `_cart_fd_opts` overrides, and any explicitly passed `opts`, with later sources taking precedence, to resolve the finite-difference options used for Cartesian-coordinate Jacobians.

        :param opts: explicit per-call overrides, taking highest precedence
        :type opts: dict
        :return: the merged finite-difference options
        :rtype: dict
        """
        ...

    def _get_cart_jacobs(self, jacs, coords=None, strip_embedding=False, **fd_opts):
        """
        Gets the specified dR/dXs

        :param jacs:
        :type jacs:
        :return:
        :rtype:
        """
        ...

    @property
    def embedding_coords(self):
        """
        **LLM Docstring**

        The indices of the internal-coordinate system's fixed embedding coordinates (translation/rotation degrees of freedom), if any are defined.

        :return: the embedding-coordinate indices, or `None`
        :rtype: np.ndarray | None
        """
        ...

    def _get_embedding_coords(self):
        """
        **LLM Docstring**

        Look up the embedding-coordinate indices from the internal-coordinate system, trying its `embedding_coords` attribute first and falling back to its `converter_options['embedding_coords']`.

        :return: the embedding-coordinate indices, or `None` if neither source defines them
        :rtype: np.ndarray | None
        """
        ...
    cartesian_by_internals_method = 'fast'

    def get_cartesians_by_internals(self, order=None, strip_embedding=False, reembed=True, method=None, coords=None, **fd_opts):
        """
        **LLM Docstring**

        Compute the Cartesians-by-internals Jacobian expansion up to the requested `order`, either via the fast route (inverting the internals-by-Cartesians Jacobian through a translation/rotation-invariant reduction, with caching) or the classic finite-difference/analytic route, depending on `method` (auto-selected based on the internal-coordinate system type) and whether `reembed`/`strip_embedding`/explicit `coords` are requested.

        :param order: the highest derivative order to compute; if `None`, returns whatever is already cached
        :type order: int | None
        :param strip_embedding: whether to strip the fixed embedding coordinates from the result
        :type strip_embedding: bool
        :param reembed: whether to use the Eckart-reembedded (translation/rotation-invariant) formulation, for the `'fast'` method
        :type reembed: bool
        :param method: which computation strategy to use (`'fast'` or `'classic'`); auto-selected if `None`
        :type method: str | None
        :param coords: alternate Cartesian coordinates to compute the Jacobian at, instead of this embedding's own geometry
        :type coords: np.ndarray | None
        :param fd_opts: extra finite-difference options forwarded to the underlying Jacobian computation
        :type fd_opts: dict
        :return: the Cartesians-by-internals Jacobian tensors, one per order
        :rtype: list[np.ndarray]
        :raises ValueError: if fewer cached/computed orders are available than requested (classic route)
        """
        ...

    def get_internals_by_cartesians(self, order=None, strip_embedding=False, coords=None, **opts):
        """
        **LLM Docstring**

        Compute the internals-by-Cartesians Jacobian expansion up to the requested `order`, via finite difference/analytic derivatives (through `_get_cart_jacobs`), reshaping the results to `(..., ncart, ncart, ..., nint)`-style tensors and optionally stripping the fixed embedding coordinates.

        :param order: the highest derivative order to compute; if `None`, returns whatever is already cached
        :type order: int | None
        :param strip_embedding: whether to strip the fixed embedding coordinates from the result
        :type strip_embedding: bool
        :param coords: alternate Cartesian coordinates to compute the Jacobian at
        :type coords: np.ndarray | None
        :param opts: extra finite-difference options forwarded to `_get_cart_jacobs`
        :type opts: dict
        :return: the internals-by-Cartesians Jacobian tensors, one per order
        :rtype: list[np.ndarray]
        :raises ValueError: if fewer cached/computed orders are available than requested
        """
        ...

    def embed_coords(self, coords, sel=None, in_paf=False, planar_ref_tolerance=None, proper_rotation=True):
        """
        **LLM Docstring**

        Eckart-embed a set of Cartesian coordinates onto this embedding's reference geometry.

        :param coords: the coordinates to embed
        :type coords: np.ndarray
        :param sel: subset of atoms to use for the embedding fit
        :type sel: Iterable[int] | None
        :param in_paf: whether to embed into the principal-axis frame
        :type in_paf: bool
        :param planar_ref_tolerance: tolerance for detecting a (near-)planar reference structure
        :type planar_ref_tolerance: float | None
        :param proper_rotation: whether to restrict the embedding to proper (determinant +1) rotations
        :type proper_rotation: bool
        :return: the Eckart-embedded coordinates
        :rtype: np.ndarray
        """
        ...

    def unembed_derivs(self, coords, derivs, sel=None, in_paf=False, planar_ref_tolerance=None):
        """
        **LLM Docstring**

        Undo an Eckart embedding's rotation on a set of Cartesian derivative tensors, transforming them back by the combination of the embedding's axis frame and rotation.

        :param coords: the (embedded) coordinates the derivatives were computed at
        :type coords: np.ndarray
        :param derivs: the Cartesian derivative tensors to un-rotate
        :type derivs: list[np.ndarray]
        :param sel: subset of atoms used for the embedding fit
        :type sel: Iterable[int] | None
        :param in_paf: whether the embedding used the principal-axis frame
        :type in_paf: bool
        :param planar_ref_tolerance: tolerance for detecting a (near-)planar reference structure
        :type planar_ref_tolerance: float | None
        :return: the un-rotated derivative tensors
        :rtype: list[np.ndarray]
        """
        ...

    @property
    def inertia_tensor(self):
        """
        **LLM Docstring**

        The molecule's inertia tensor at its current Cartesian coordinates.

        :return: the inertia tensor
        :rtype: np.ndarray
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

        The first and second derivatives of the inertia tensor with respect to mass-weighted Cartesian displacements.

        :return: `[I0Y, I0YY]`, as returned by `StructuralProperties.get_prop_inertial_frame_derivatives`
        :rtype: list[np.ndarray]
        """
        ...

    @property
    def translation_rotation_modes(self):
        """
        **LLM Docstring**

        The (cached) translation and rotation eigenvectors of the molecule at its current geometry.

        :return: the translation/rotation eigenvectors
        :rtype: tuple
        """
        ...

    def get_translation_rotation_invariant_transformation(self, order=0, mass_weighted=True, strip_embedding=True, coords=None):
        """
        **LLM Docstring**

        Build the transformation (and its inverse) that projects out the translational and rotational degrees of freedom from a set of Cartesian coordinates.

        :param order: the derivative order of the transformation to build
        :type order: int
        :param mass_weighted: whether the transformation should act on mass-weighted coordinates
        :type mass_weighted: bool
        :param strip_embedding: whether to strip the fixed embedding coordinates from the result
        :type strip_embedding: bool
        :param coords: alternate Cartesian coordinates to build the transformation at, instead of this embedding's own geometry
        :type coords: np.ndarray | None
        :return: the translation/rotation-invariant transformation and its inverse
        :rtype: tuple
        """
        ...

class ModeEmbedding:
    """
    Provides a specialization on a `MoleculaEmbedding` to express all properties
    in terms of the attendant normal modes
    """
    modes: NormalModes

    def __init__(self, embedding: MolecularEmbedding, modes, mass_weight=None, dimensionless=None, masses=None):
        """
        **LLM Docstring**

        Set up a mode-aware specialization of a `MolecularEmbedding`: resolves whatever form `modes` was given in (a manager, a `MolecularVibrations`, or a raw normal-modes object) down to a plain normal-modes object, optionally converting it to a dimensionless or mass-weighted basis, and records whether the resulting modes are mass-weighted.

        :param embedding: the underlying molecular embedding to specialize
        :type embedding: MolecularEmbedding
        :param modes: the normal modes (or a manager/vibrations object wrapping them) to express properties in terms of
        :type modes: object
        :param mass_weight: whether to convert the modes to a mass-weighted basis
        :type mass_weight: bool | None
        :param dimensionless: whether to convert the modes to a dimensionless basis
        :type dimensionless: bool | None
        :param masses: masses to use for the mass-weighting/dimensionless conversions; defaults to the embedding's own masses
        :type masses: np.ndarray | None
        :return: None
        :rtype: None
        """
        ...

    def mw_conversion(self, strip_dummies=None):
        """
        **LLM Docstring**

        Build the diagonal mass-weighting matrix (`sqrt(mass)` per Cartesian coordinate) used to convert plain Cartesian derivatives into mass-weighted ones.

        :param strip_dummies: whether to exclude dummy (non-positive-mass) atoms from the mass vector
        :type strip_dummies: bool | None
        :return: the diagonal mass-weighting matrix
        :rtype: np.ndarray
        """
        ...

    def mw_inverse(self, strip_dummies=None):
        """
        **LLM Docstring**

        Build the diagonal inverse-mass-weighting matrix (`1/sqrt(mass)` per Cartesian coordinate) used to convert mass-weighted Cartesian derivatives back into plain ones.

        :param strip_dummies: whether to exclude dummy (non-positive-mass) atoms from the mass vector
        :type strip_dummies: bool | None
        :return: the diagonal inverse-mass-weighting matrix
        :rtype: np.ndarray
        """
        ...

    def get_mw_cartesians_by_internals(self, order=None, mass_weighted=None, coords=None, strip_embedding=True):
        """
        **LLM Docstring**

        Fetch the Cartesians-by-internals Jacobian expansion, converted to mass-weighted Cartesians if `mass_weighted` (or `self.mass_weighted` by default) is set.

        :param order: the highest derivative order to compute
        :type order: int | None
        :param mass_weighted: whether to mass-weight the result; defaults to `self.mass_weighted`
        :type mass_weighted: bool | None
        :param coords: alternate coordinates to compute the Jacobian at
        :type coords: np.ndarray | None
        :param strip_embedding: whether to strip the fixed embedding coordinates
        :type strip_embedding: bool
        :return: the (optionally mass-weighted) Cartesians-by-internals Jacobian tensors
        :rtype: list[np.ndarray]
        """
        ...

    def get_internals_by_mw_cartesians(self, order=None, mass_weighted=None, coords=None, strip_embedding=True):
        """
        **LLM Docstring**

        Fetch the internals-by-Cartesians Jacobian expansion, converted to be with respect to mass-weighted Cartesians if `mass_weighted` (or `self.mass_weighted` by default) is set.

        :param order: the highest derivative order to compute
        :type order: int | None
        :param mass_weighted: whether to mass-weight the result; defaults to `self.mass_weighted`
        :type mass_weighted: bool | None
        :param coords: alternate coordinates to compute the Jacobian at
        :type coords: np.ndarray | None
        :param strip_embedding: whether to strip the fixed embedding coordinates
        :type strip_embedding: bool
        :return: the (optionally mass-weighted) internals-by-Cartesians Jacobian tensors
        :rtype: list[np.ndarray]
        """
        ...

    def get_internals_by_cartesians(self, order=None, coords=None, strip_embedding=True):
        """
        expresses raw internals or modes (internals or Cartesian) in terms of mass-weighted Cartesians

        :param order:
        :param strip_embedding:
        :return:
        """
        ...

    def get_cartesians_by_internals(self, order=None, coords=None, strip_embedding=True):
        """
        expresses raw internals or modes (internals or Cartesian) in terms of mass-weighted Cartesians

        :param order:
        :param strip_embedding:
        :return:
        """
        ...

    def get_inertia_tensor_expansion(self, order=None, strip_embedding=True):
        """
        **LLM Docstring**

        Compute the Taylor expansion of the inertia tensor in terms of this embedding's coordinates (internal coordinates or normal modes), by re-expanding the inertial-frame derivatives through the Cartesians-by-internals Jacobian.

        :param order: the highest derivative order to compute
        :type order: int | None
        :param strip_embedding: whether to strip the fixed embedding coordinates from the underlying Jacobian
        :type strip_embedding: bool
        :return: `[I0] + [derivative terms...]`, the inertia tensor and its derivatives
        :rtype: list[np.ndarray]
        """
        ...

    def get_inertial_frame(self):
        """
        **LLM Docstring**

        The molecule's inertial (principal-axis) frame, delegated to the underlying embedding.

        :return: the inertial frame, as returned by `MolecularEmbedding.inertial_frame`
        :rtype: tuple
        """
        ...

    def get_modes_by_coords(self, mass_weighted=None, frequency_scaled=None):
        """
        **LLM Docstring**

        The modes-by-coordinates transformation matrix, optionally adjusting the mass-weighting/frequency-scaling convention of the modes first, and (if the modes are Cartesian but expressed relative to an internal-coordinate embedding) re-expressing them in terms of internal coordinates.

        :param mass_weighted: `True`/`False` to force mass-weighting on/off before extracting the matrix; `None` to leave the modes' current convention
        :type mass_weighted: bool | None
        :param frequency_scaled: `True`/`False` to adjust frequency scaling before extracting the matrix (both branches currently call `remove_frequency_scaling`); `None` to leave it unchanged
        :type frequency_scaled: bool | None
        :return: the modes-by-coordinates matrix, or `None` if no modes are set
        :rtype: np.ndarray | None
        """
        ...

    def get_coords_by_modes(self, mass_weighted=None, frequency_scaled=None):
        """
        **LLM Docstring**

        The coordinates-by-modes transformation matrix, optionally adjusting the mass-weighting/frequency-scaling convention of the modes first, and (if the modes are Cartesian but expressed relative to an internal-coordinate embedding) re-expressing them in terms of internal coordinates.

        :param mass_weighted: `True`/`False` to force mass-weighting on/off before extracting the matrix; `None` to leave the modes' current convention
        :type mass_weighted: bool | None
        :param frequency_scaled: `True`/`False` to adjust frequency scaling before extracting the matrix (both branches currently call `remove_frequency_scaling`); `None` to leave it unchanged
        :type frequency_scaled: bool | None
        :return: the coordinates-by-modes matrix, or `None` if no modes are set
        :rtype: np.ndarray | None
        """
        ...

def _get_best_axes(first_pos, axes):
    """
    Determine the best pair of inertial axes so that we don't get large-scale breakdowns from the choice of embedding

    :param first_pos:
    :type first_pos:
    :param axes:
    :type axes:
    :return:
    :rtype:
    """
    ...

class MolecularGenericInternalCoordinateSystem(GenericInternalCoordinateSystem):
    """
    Mirrors the standard ZMatrix coordinate system in _almost_ all regards, but forces an embedding
    """
    name = 'MolecularGenericInternals'

    class PassThroughRedundantGenerator:

        def __init__(self, redundant_transformation, redundant_inverse=None, masses=None, untransformed_coordinates=None, relocalize=False, **opts):
            """
            **LLM Docstring**

            Store an already-known redundant-coordinate transformation so it can be handed back unchanged (rather than computed) when `get_redundant_transformation` is called.

            :param redundant_transformation: the redundant-to-non-redundant transformation matrix to pass through
            :type redundant_transformation: np.ndarray
            :param redundant_inverse: the corresponding inverse transformation, if available
            :type redundant_inverse: np.ndarray | None
            :param masses: the atomic masses (stored as a 1-tuple due to a trailing comma in the assignment)
            :type masses: np.ndarray | None
            :param untransformed_coordinates: coordinates to leave untransformed
            :type untransformed_coordinates: object | None
            :param relocalize: whether the redundant coordinates should be relocalized
            :type relocalize: bool
            :param opts: extra options, stored but not otherwise used
            :type opts: dict
            :return: None
            :rtype: None
            """
            ...

        def get_redundant_transformation(self, base_expansions, **ignored):
            """
            **LLM Docstring**

            Return the stored redundant transformation unchanged, along with the base derivative expansions re-expanded through it.

            :param base_expansions: the derivative expansions to re-express through the redundant transformation
            :type base_expansions: list[np.ndarray]
            :param ignored: any other arguments, accepted but not used
            :type ignored: dict
            :return: `(self.tf, reexpanded_expansions)`
            :rtype: tuple[np.ndarray, list[np.ndarray]]
            """
            ...

    def __init__(self, masses, coords, /, specs, converter_options=None, redundant=False, relocalize=False, untransformed_coordinates=None, redundant_transformation=None, redundant_inverse=None, angle_ordering='ijk', **opts):
        """

        :param molecule:
        :type molecule: AbstractMolecule
        :param converter_options:
        :type converter_options:
        :param opts:
        :type opts:
        """
        ...

class MolecularZMatrixCoordinateSystem(ZMatrixCoordinateSystem):
    """
    Mirrors the standard ZMatrix coordinate system in _almost_ all regards, but forces an embedding
    """
    name = 'MolecularZMatrix'
    embedding_coords = [0, 1, 2, 4, 5, 8]

    def __init__(self, masses, coords, converter_options=None, **opts):
        """

        :param molecule:
        :type molecule: AbstractMolecule
        :param converter_options:
        :type converter_options:
        :param opts:
        :type opts:
        """
        ...

    @property
    def origins(self):
        """
        **LLM Docstring**

        The Z-matrix embedding's origin points (typically the reference center of mass), from `converter_options['origins']`.

        :return: the origin points
        :rtype: np.ndarray
        """
        ...

    @property
    def axes(self):
        """
        **LLM Docstring**

        The Z-matrix embedding's reference axes (typically two principal axes), from `converter_options['axes']`.

        :return: the reference axes
        :rtype: np.ndarray
        """
        ...

    def get_direct_converter(self, target):
        """
        **LLM Docstring**

        Provide a converter from this molecular Z-matrix system directly to the plain (non-molecular) `ZMatrix` coordinate system, if `target` is one.

        :param target: the coordinate system being converted to
        :type target: object
        :return: a `MolecularZMatrixToRegularZMatrixConverter`, or `None` if `target` isn't a `ZMatrix` system
        :rtype: MolecularZMatrixToRegularZMatrixConverter | None
        """
        ...

    def get_inverse_converter(self, target):
        """
        **LLM Docstring**

        Provide a converter from the plain `ZMatrix` coordinate system into this molecular Z-matrix system, if `target` is one.

        :param target: the coordinate system being converted from
        :type target: object
        :return: a `RegularZMatrixToMolecularZMatrixConverter`, or `None` if `target` isn't a `ZMatrix` system
        :rtype: RegularZMatrixToMolecularZMatrixConverter | None
        """
        ...

    def pre_convert_to(self, system, opts=None):
        """
        **LLM Docstring**

        Re-establish the embedding options (via `set_embedding`) before delegating to the base class's `pre_convert_to`, preserving the existing atom `ordering` if the caller supplied its own options dict.

        :param system: the coordinate system being converted to
        :type system: object
        :param opts: explicit conversion options to use instead of `self.converter_options`
        :type opts: dict | None
        :return: the resolved conversion options
        :rtype: dict
        """
        ...

    def set_embedding(self):
        """
        **LLM Docstring**

        (Re)compute and store this Z-matrix system's embedding options -- the reference origin, reference axes (chosen via `_get_best_axes` to avoid ill-conditioned choices), axis labels, masses, dummy-atom positions, and reference coordinates -- based on the current center of mass and inertial axes, fixing up the Z-matrix `ordering`'s first three rows to reference the embedding's dummy origin/axis points if an ordering is present.

        :return: None
        :rtype: None
        """
        ...

    def jacobian(self, coords, *args, reembed=None, strip_dummies=None, converter_options=None, **kwargs):
        """
        **LLM Docstring**

        Compute the Jacobian of this Z-matrix system with respect to Cartesian coordinates, handling batched/multi-frame inputs, optional dummy-atom stripping, and temporarily overriding the `reembed` converter option for the duration of the call.

        :param coords: the Cartesian coordinates to compute the Jacobian at
        :type coords: np.ndarray
        :param args: extra positional arguments forwarded to the base class's `jacobian`
        :type args: tuple
        :param reembed: whether to re-embed (Eckart-align) during the Jacobian calculation; falls back to the converter options, then defaults to `True`
        :type reembed: bool | None
        :param strip_dummies: whether to exclude dummy-atom coordinates from the Jacobian; falls back to the converter options, then defaults to `False`
        :type strip_dummies: bool | None
        :param converter_options: extra converter options merged with `self.converter_options`
        :type converter_options: dict | None
        :param kwargs: extra keyword arguments forwarded to the base class's `jacobian`
        :type kwargs: dict
        :return: the computed Jacobian tensor(s)
        :rtype: list[np.ndarray]
        """
        ...

class MolecularCartesianCoordinateSystem(CartesianCoordinateSystem):
    """
    Mirrors the standard Cartesian coordinate system in _almost_ all regards, but forces an embedding
    """
    name = 'MolecularCartesians'

    def __init__(self, masses, coords, dummy_positions=None, converter_options=None, **opts):
        """

        :param molecule:
        :type molecule: AbstractMolecule
        :param converter_options:
        :type converter_options:
        :param opts:
        :type opts:
        """
        ...

    def to_state(self, serializer=None):
        """
        **LLM Docstring**

        Serialize this coordinate system's state, adding the masses, coordinates, and dummy-atom positions on top of whatever the base class's `to_state` produces.

        :param serializer: the serializer to use, forwarded to the base class
        :type serializer: object | None
        :return: the serialized state dict
        :rtype: dict
        """
        ...

    @classmethod
    def from_state(cls, data, serializer=None):
        """
        **LLM Docstring**

        Reconstruct a `MolecularCartesianCoordinateSystem` from a previously serialized state dict.

        :param data: the serialized state, as produced by `to_state`
        :type data: dict
        :param serializer: the serializer to use, accepted for interface consistency but not used in this method's body
        :type serializer: object | None
        :return: the reconstructed coordinate system
        :rtype: MolecularCartesianCoordinateSystem
        """
        ...

    def pre_convert_to(self, system, opts=None):
        """
        **LLM Docstring**

        Ensure the masses are up to date in `converter_options`, and, if converting to a Z-matrix-family system, re-establish the embedding options (via `set_embedding`) before delegating to the base class's `pre_convert_to`.

        :param system: the coordinate system being converted to
        :type system: object
        :param opts: explicit conversion options to use instead of `self.converter_options`
        :type opts: dict | None
        :return: the resolved conversion options
        :rtype: dict
        """
        ...

    def set_embedding(self):
        """
        Sets up the embedding options...
        :return:
        :rtype:
        """
        ...

    def jacobian(self, coords, system, order=None, strip_dummies=None, converter_options=None, analytic_deriv_order=None, **kwargs):
        """
        **LLM Docstring**

        Compute the Jacobian of these Cartesian coordinates with respect to `system`, resolving the analytic-derivative order (defaulting to purely numerical for Z-matrix targets) and optionally excluding dummy-atom coordinates.

        :param coords: the coordinates to compute the Jacobian at
        :type coords: np.ndarray
        :param system: the target coordinate system
        :type system: object
        :param order: the derivative order(s) to compute
        :type order: int | list[int] | None
        :param strip_dummies: whether to exclude dummy-atom coordinates; falls back to the converter options, then defaults to `False`
        :type strip_dummies: bool | None
        :param converter_options: extra converter options merged with `self.converter_options`
        :type converter_options: dict | None
        :param analytic_deriv_order: the order up to which to compute the Jacobian analytically rather than numerically; falls back to the converter options, then defaults based on whether `system` is a Z-matrix
        :type analytic_deriv_order: int | None
        :param kwargs: extra keyword arguments forwarded to the base class's `jacobian`
        :type kwargs: dict
        :return: the computed Jacobian tensor(s)
        :rtype: list[np.ndarray]
        """
        ...

class MolecularCartesianToZMatrixConverter(CoordinateSystemConverter):
    """
    ...
    """

    def __init__(self, cart_system, zmat_system, **opts):
        """
        **LLM Docstring**

        Store the `(cart_system, zmat_system)` type pair this converter handles.

        :param cart_system: the source Cartesian coordinate system
        :type cart_system: object
        :param zmat_system: the target Z-matrix coordinate system
        :type zmat_system: object
        :param opts: extra options forwarded to the base `CoordinateSystemConverter`
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...

    @property
    def types(self):
        """
        **LLM Docstring**

        The `(source, target)` coordinate-system type pair this converter handles.

        :return: `(cart_system, zmat_system)`
        :rtype: tuple
        """
        ...

    def convert(self, coords, *, masses, dummy_positions, origins=None, axes=None, ordering=None, **kwargs):
        """
        Converts from Cartesian to ZMatrix coords, preserving the embedding
        :param coords:
        :type coords: CoordinateSet
        :param molecule:
        :type molecule:
        :param origins:
        :type origins:
        :param axes:
        :type axes:
        :param ordering:
        :type ordering:
        :param kwargs:
        :type kwargs:
        :return:
        :rtype:
        """
        ...
    base_cartesian_type = CartesianCoordinates3D
    base_internal_type = ZMatrixCoordinates

    def convert_many(self, coords, *, masses, dummy_positions, origins=None, axes=None, ordering=None, strip_embedding=True, strip_dummies=False, return_derivs=None, **kwargs):
        """
        Converts from Cartesian to ZMatrix coords, preserving the embedding

        :param coords: coordinates in Cartesians to convert
        :type coords: np.ndarray
        :param molecule:
        :type molecule: MolecularEmbedding
        :param origins: the origin for each individual structure
        :type origins: np.ndarray
        :param axes: the axes for each structure
        :type axes: np.ndarray
        :param ordering: the Z-matrix ordering spec
        :type ordering:
        :param strip_embedding: whether to strip the embedding coordinates
        :type strip_embedding:
        :param strip_dummies: whether to strip all dummy coordinates
        :type strip_dummies:
        :param kwargs:
        :type kwargs:
        :return:
        :rtype:
        """
        ...

class MolecularCartesianToRegularCartesianConverter(CoordinateSystemConverter):
    """
    ...
    """

    def __init__(self, cart_system, **opts):
        """
        **LLM Docstring**

        Store the `(CartesianCoordinates3D, cart_system)` type pair this converter handles.

        :param cart_system: the molecular Cartesian coordinate system being converted from
        :type cart_system: object
        :param opts: extra options forwarded to the base `CoordinateSystemConverter`
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...

    @property
    def types(self):
        """
        **LLM Docstring**

        The `(source, target)` coordinate-system type pair this converter handles.

        :return: `(CartesianCoordinates3D, cart_system)`
        :rtype: tuple
        """
        ...

    def convert(self, coords, **kw):
        """
        **LLM Docstring**

        Pass the coordinates through unchanged, since a molecular Cartesian system and a plain Cartesian system share the same underlying representation.

        :param coords: the coordinates to convert
        :type coords: np.ndarray
        :param kw: extra keyword arguments, returned unchanged as the options dict
        :type kw: dict
        :return: `(coords, kw)`
        :rtype: tuple
        """
        ...

    def convert_many(self, coords, **kwargs):
        """
        Converts from Cartesian to ZMatrix coords, preserving the embedding
        """
        ...

class RegularCartesianToMolecularCartesianConverter(CoordinateSystemConverter):
    """
    ...
    """

    def __init__(self, cart_system, **opts):
        """
        **LLM Docstring**

        Store the `(cart_system, CartesianCoordinates3D)` type pair this converter handles.

        :param cart_system: the molecular Cartesian coordinate system being converted to
        :type cart_system: object
        :param opts: extra options forwarded to the base `CoordinateSystemConverter`
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...

    @property
    def types(self):
        """
        **LLM Docstring**

        The `(source, target)` coordinate-system type pair this converter handles.

        :return: `(cart_system, CartesianCoordinates3D)`
        :rtype: tuple
        """
        ...

    def convert(self, coords, **kw):
        """
        **LLM Docstring**

        Pass the coordinates through unchanged, since a plain Cartesian system and a molecular Cartesian system share the same underlying representation.

        :param coords: the coordinates to convert
        :type coords: np.ndarray
        :param kw: extra keyword arguments, returned unchanged as the options dict
        :type kw: dict
        :return: `(coords, kw)`
        :rtype: tuple
        """
        ...

    def convert_many(self, coords, **kwargs):
        """
        Converts from Cartesian to ZMatrix coords, preserving the embedding
        """
        ...

class MolecularZMatrixToCartesianConverter(CoordinateSystemConverter):
    """
    ...
    """

    def __init__(self, zmat_system, cart_system, **opts):
        """
        **LLM Docstring**

        Store the `(zmat_system, cart_system)` type pair this converter handles.

        :param zmat_system: the source Z-matrix coordinate system
        :type zmat_system: object
        :param cart_system: the target Cartesian coordinate system
        :type cart_system: object
        :param opts: extra options forwarded to the base `CoordinateSystemConverter`
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...

    @property
    def types(self):
        """
        **LLM Docstring**

        The `(source, target)` coordinate-system type pair this converter handles.

        :return: `(zmat_system, cart_system)`
        :rtype: tuple
        """
        ...

    def convert(self, coords, **kw):
        """
        **LLM Docstring**

        Convert a single frame of Z-matrix coordinates to Cartesians by delegating to `convert_many` on a length-1 batch and unwrapping the result (including the per-frame `'derivs'` entries, if present).

        :param coords: the single-frame Z-matrix coordinates to convert
        :type coords: np.ndarray
        :param kw: extra keyword arguments forwarded to `convert_many`
        :type kw: dict
        :return: `(cartesian_coords, opts)` for the single frame
        :rtype: tuple
        """
        ...
    base_cartesian_type = CartesianCoordinates3D
    base_internal_type = ZMatrixCoordinates

    def convert_many(self, coords, *, masses, dummy_positions, ref_coords, origins=None, axes=None, ordering=None, reembed=False, axes_choice=None, return_derivs=None, strip_dummies=False, strip_embedding=True, planar_ref_tolerance=None, embedding_masses=None, spec=None, **kwargs):
        """
        Converts from Cartesian to ZMatrix coords, attempting to preserve the embedding
        """
        ...

class MolecularZMatrixToRegularZMatrixConverter(CoordinateSystemConverter):
    """
    ...
    """

    def __init__(self, zmat_system, **opts):
        """
        **LLM Docstring**

        Store the `(zmat_system, ZMatrixCoordinateSystem)` type pair this converter handles.

        :param zmat_system: the molecular Z-matrix coordinate system being converted from
        :type zmat_system: object
        :param opts: extra options forwarded to the base `CoordinateSystemConverter`
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...

    @property
    def types(self):
        """
        **LLM Docstring**

        The `(source, target)` coordinate-system type pair this converter handles.

        :return: `(zmat_system, ZMatrixCoordinateSystem)`
        :rtype: tuple
        """
        ...

    def convert(self, coords, **kw):
        """
        **LLM Docstring**

        Pass the coordinates through unchanged, since a molecular Z-matrix system and a plain Z-matrix system share the same underlying representation.

        :param coords: the coordinates to convert
        :type coords: np.ndarray
        :param kw: extra keyword arguments, returned unchanged as the options dict
        :type kw: dict
        :return: `(coords, kw)`
        :rtype: tuple
        """
        ...

    def convert_many(self, coords, **kwargs):
        """
        **LLM Docstring**

        Pass a batch of coordinates through unchanged, since a molecular Z-matrix system and a plain Z-matrix system share the same underlying representation.

        :param coords: the batch of coordinates to convert
        :type coords: np.ndarray
        :param kwargs: extra keyword arguments, returned unchanged as the options dict
        :type kwargs: dict
        :return: `(coords, kwargs)`
        :rtype: tuple
        """
        ...

class RegularZMatrixToMolecularZMatrixConverter(CoordinateSystemConverter):
    """
    ...
    """

    def __init__(self, zmat_system, **opts):
        """
        **LLM Docstring**

        Store the `(ZMatrixCoordinateSystem, zmat_system)` type pair this converter handles.

        :param zmat_system: the molecular Z-matrix coordinate system being converted to
        :type zmat_system: object
        :param opts: extra options forwarded to the base `CoordinateSystemConverter`
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...

    @property
    def types(self):
        """
        **LLM Docstring**

        The `(source, target)` coordinate-system type pair this converter handles.

        :return: `(ZMatrixCoordinateSystem, zmat_system)`
        :rtype: tuple
        """
        ...

    def convert(self, coords, **kw):
        """
        **LLM Docstring**

        Pass the coordinates through unchanged, since a plain Z-matrix system and a molecular Z-matrix system share the same underlying representation.

        :param coords: the coordinates to convert
        :type coords: np.ndarray
        :param kw: extra keyword arguments, returned unchanged as the options dict
        :type kw: dict
        :return: `(coords, kw)`
        :rtype: tuple
        """
        ...

    def convert_many(self, coords, **kwargs):
        """
        **LLM Docstring**

        Pass a batch of coordinates through unchanged, since a plain Z-matrix system and a molecular Z-matrix system share the same underlying representation.

        :param coords: the batch of coordinates to convert
        :type coords: np.ndarray
        :param kwargs: extra keyword arguments, returned unchanged as the options dict
        :type kwargs: dict
        :return: `(coords, kwargs)`
        :rtype: tuple
        """
        ...

class MolecularCartesianToGICConverter(CartesianToGICSystemConverter):
    """
    ...
    """

    def __init__(self, cart_system, zmat_system, **opts):
        """
        **LLM Docstring**

        Store the `(cart_system, zmat_system)` type pair this converter handles, and precompute which internal-coordinate specs are periodic (angle-like specs with more than 2 atoms) for later use in handling periodic wraparound.

        :param cart_system: the source Cartesian coordinate system
        :type cart_system: object
        :param zmat_system: the target generic-internal coordinate system
        :type zmat_system: object
        :param opts: extra options forwarded to the base converter
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...

    @property
    def types(self):
        """
        **LLM Docstring**

        The `(source, target)` coordinate-system type pair this converter handles.

        :return: `(cart_system, zmat_system)`
        :rtype: tuple
        """
        ...

    def convert_many(self, coords, *, order=0, return_derivs=None, redundant_generator: RedundantCoordinateGenerator=None, reference_internals=None, redundant_transformation=None, handle_periodicity=True, **kw):
        """
        We'll implement this by having the ordering arg wrap around in coords?
        """
        ...

class MolecularGICToCartesianConverter(GICSystemToCartesianConverter):
    """
    ...
    """

    def __init__(self, cart_system, zmat_system, **opts):
        """
        **LLM Docstring**

        Store the `(cart_system, zmat_system)` type pair this converter handles.

        :param cart_system: the target Cartesian coordinate system
        :type cart_system: object
        :param zmat_system: the source generic-internal coordinate system
        :type zmat_system: object
        :param opts: extra options forwarded to the base converter
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...

    @property
    def types(self):
        """
        **LLM Docstring**

        The `(source, target)` coordinate-system type pair this converter handles.

        :return: `(cart_system, zmat_system)`
        :rtype: tuple
        """
        ...

    def convert_many(self, coords, *, redundant_transformation=None, redundant_inverse=None, redundant_generator=None, reference_internals=None, **kw):
        """
        **LLM Docstring**

        Convert a batch of generic-internal coordinates to Cartesians, forwarding any redundant transformation (and its transpose as the corresponding inverse) as the `transformations` argument to the base converter.

        :param coords: the batch of internal coordinates to convert
        :type coords: np.ndarray
        :param redundant_transformation: the redundant-to-non-redundant transformation used to build the internal coordinates, if any
        :type redundant_transformation: np.ndarray | None
        :param redundant_inverse: accepted but not directly used (recomputed from `redundant_transformation` if needed)
        :type redundant_inverse: np.ndarray | None
        :param redundant_generator: accepted for interface consistency but not used in this method's body
        :type redundant_generator: object | None
        :param reference_internals: the reference internal coordinates the redundant transformation is defined relative to
        :type reference_internals: np.ndarray | None
        :param kw: extra keyword arguments forwarded to the base converter's `convert_many`
        :type kw: dict
        :return: `(carts, opts)`
        :rtype: tuple
        """
        ...

class RegularGICToMolecularGICConverter(CoordinateSystemConverter):
    """
    ...
    """

    def __init__(self, **opts):
        """
        **LLM Docstring**

        Store the `(GenericInternalCoordinates, MolecularGenericInternalCoordinateSystem)` type pair this converter handles.

        :param opts: extra options forwarded to the base `CoordinateSystemConverter`
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...

    @property
    def types(self):
        """
        **LLM Docstring**

        The `(source, target)` coordinate-system type pair this converter handles.

        :return: `(GenericInternalCoordinates, MolecularGenericInternalCoordinateSystem)`
        :rtype: tuple
        """
        ...

    def convert(self, coords, **kw):
        """
        **LLM Docstring**

        Pass the coordinates through unchanged, since a plain generic-internal system and a molecular generic-internal system share the same underlying representation.

        :param coords: the coordinates to convert
        :type coords: np.ndarray
        :param kw: extra keyword arguments, returned unchanged as the options dict
        :type kw: dict
        :return: `(coords, kw)`
        :rtype: tuple
        """
        ...

    def convert_many(self, coords, **kwargs):
        """
        **LLM Docstring**

        Pass a batch of coordinates through unchanged, since a plain generic-internal system and a molecular generic-internal system share the same underlying representation.

        :param coords: the batch of coordinates to convert
        :type coords: np.ndarray
        :param kwargs: extra keyword arguments, returned unchanged as the options dict
        :type kwargs: dict
        :return: `(coords, kwargs)`
        :rtype: tuple
        """
        ...

class MolecularGICConverterToRegularGIC(CoordinateSystemConverter):
    """
    ...
    """

    def __init__(self, **opts):
        """
        **LLM Docstring**

        Store the `(MolecularGenericInternalCoordinateSystem, GenericInternalCoordinates)` type pair this converter handles.

        :param opts: extra options forwarded to the base `CoordinateSystemConverter`
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...

    @property
    def types(self):
        """
        **LLM Docstring**

        The `(source, target)` coordinate-system type pair this converter handles.

        :return: `(MolecularGenericInternalCoordinateSystem, GenericInternalCoordinates)`
        :rtype: tuple
        """
        ...

    def convert(self, coords, **kw):
        """
        **LLM Docstring**

        Pass the coordinates through unchanged, since a molecular generic-internal system and a plain generic-internal system share the same underlying representation.

        :param coords: the coordinates to convert
        :type coords: np.ndarray
        :param kw: extra keyword arguments, returned unchanged as the options dict
        :type kw: dict
        :return: `(coords, kw)`
        :rtype: tuple
        """
        ...

    def convert_many(self, coords, **kwargs):
        """
        **LLM Docstring**

        Pass a batch of coordinates through unchanged, since a molecular generic-internal system and a plain generic-internal system share the same underlying representation.

        :param coords: the batch of coordinates to convert
        :type coords: np.ndarray
        :param kwargs: extra keyword arguments, returned unchanged as the options dict
        :type kwargs: dict
        :return: `(coords, kwargs)`
        :rtype: tuple
        """
        ...

class MolecularIZCoordinateSystem(MolecularZMatrixCoordinateSystem):
    name = 'MolecularIZMatrix'

class MolecularCartesianToIZConverter(MolecularCartesianToZMatrixConverter):
    """
    ...
    """
    base_internal_type = McUtils.Coordinerds.IterativeZMatrixCoordinates

class MolecularIZToCartesianConverter(MolecularZMatrixToCartesianConverter):
    """
    ...
    """
    base_internal_type = McUtils.Coordinerds.IterativeZMatrixCoordinates

    def convert_many(self, coords, *, masses, dummy_positions, ref_coords, origins=None, axes=None, ordering=None, reembed=False, axes_choice=None, return_derivs=None, strip_dummies=False, strip_embedding=True, planar_ref_tolerance=None, embedding_masses=None, fixed_atoms=None, fixed_coords=None, **kwargs):
        """
        **LLM Docstring**

        Convert a batch of iterative-Z-matrix coordinates to Cartesians, filling in sensible defaults for the extra dummy embedding masses/fixed atoms/fixed coordinates used by the iterative scheme before delegating to the base `MolecularZMatrixToCartesianConverter.convert_many`.

        :param coords: the batch of iterative-Z-matrix coordinates to convert
        :type coords: np.ndarray
        :param masses: the atomic masses
        :type masses: np.ndarray
        :param dummy_positions: indices of dummy atoms
        :type dummy_positions: list[int]
        :param ref_coords: the reference Cartesian coordinates for the embedding
        :type ref_coords: np.ndarray
        :param origins: the embedding origin points
        :type origins: np.ndarray | None
        :param axes: the embedding reference axes
        :type axes: np.ndarray | None
        :param ordering: the Z-matrix atom ordering
        :type ordering: np.ndarray | None
        :param reembed: whether to re-embed (Eckart-align) during conversion
        :type reembed: bool
        :param axes_choice: which axis pair was chosen for the embedding
        :type axes_choice: tuple | None
        :param return_derivs: which derivative orders to return
        :type return_derivs: object | None
        :param strip_dummies: whether to exclude dummy-atom coordinates
        :type strip_dummies: bool
        :param strip_embedding: whether to strip the fixed embedding coordinates
        :type strip_embedding: bool
        :param planar_ref_tolerance: tolerance for detecting a (near-)planar reference structure
        :type planar_ref_tolerance: float | None
        :param embedding_masses: masses to use for the 3 extra dummy embedding atoms, plus the real atoms; defaults to very heavy dummy masses
        :type embedding_masses: np.ndarray | None
        :param fixed_atoms: indices treated as fixed for the iterative embedding; defaults to the first 3 (dummy) atoms
        :type fixed_atoms: list[int] | None
        :param fixed_coords: coordinate indices treated as fixed; defaults to the standard embedding coordinates offset by the 3 dummy atoms
        :type fixed_coords: list[int] | None
        :param kwargs: extra keyword arguments forwarded to the base converter
        :type kwargs: dict
        :return: `(carts, opts)`, with `opts['masses']` set to the original (non-dummy-augmented) masses
        :rtype: tuple
        """
        ...

class RegularIZToMolecularIZConverter(CoordinateSystemConverter):
    """
    ...
    """

    def __init__(self, zmat_system, **opts):
        """
        **LLM Docstring**

        Store the `(IterativeZMatrixCoordinates, zmat_system)` type pair this converter handles.

        :param zmat_system: the molecular iterative-Z-matrix coordinate system being converted to
        :type zmat_system: object
        :param opts: extra options forwarded to the base `CoordinateSystemConverter`
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...

    @property
    def types(self):
        """
        **LLM Docstring**

        The `(source, target)` coordinate-system type pair this converter handles.

        :return: `(IterativeZMatrixCoordinates, zmat_system)`
        :rtype: tuple
        """
        ...

    def convert(self, coords, **kw):
        """
        **LLM Docstring**

        Pass the coordinates through unchanged, since a plain iterative-Z-matrix system and a molecular one share the same underlying representation.

        :param coords: the coordinates to convert
        :type coords: np.ndarray
        :param kw: extra keyword arguments, returned unchanged as the options dict
        :type kw: dict
        :return: `(coords, kw)`
        :rtype: tuple
        """
        ...

    def convert_many(self, coords, **kwargs):
        """
        **LLM Docstring**

        Pass a batch of coordinates through unchanged, since a plain iterative-Z-matrix system and a molecular one share the same underlying representation.

        :param coords: the batch of coordinates to convert
        :type coords: np.ndarray
        :param kwargs: extra keyword arguments, returned unchanged as the options dict
        :type kwargs: dict
        :return: `(coords, kwargs)`
        :rtype: tuple
        """
        ...

class MolecularIZToRegularIZConverter(CoordinateSystemConverter):
    """
    ...
    """

    def __init__(self, zmat_system, **opts):
        """
        **LLM Docstring**

        Store the `(zmat_system, IterativeZMatrixCoordinates)` type pair this converter handles.

        :param zmat_system: the molecular iterative-Z-matrix coordinate system being converted from
        :type zmat_system: object
        :param opts: extra options forwarded to the base `CoordinateSystemConverter`
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...

    @property
    def types(self):
        """
        **LLM Docstring**

        The `(source, target)` coordinate-system type pair this converter handles.

        :return: `(zmat_system, IterativeZMatrixCoordinates)`
        :rtype: tuple
        """
        ...

    def convert(self, coords, **kw):
        """
        **LLM Docstring**

        Pass the coordinates through unchanged, since a molecular iterative-Z-matrix system and a plain one share the same underlying representation.

        :param coords: the coordinates to convert
        :type coords: np.ndarray
        :param kw: extra keyword arguments, returned unchanged as the options dict
        :type kw: dict
        :return: `(coords, kw)`
        :rtype: tuple
        """
        ...

    def convert_many(self, coords, **kwargs):
        """
        **LLM Docstring**

        Pass a batch of coordinates through unchanged, since a molecular iterative-Z-matrix system and a plain one share the same underlying representation.

        :param coords: the batch of coordinates to convert
        :type coords: np.ndarray
        :param kwargs: extra keyword arguments, returned unchanged as the options dict
        :type kwargs: dict
        :return: `(coords, kwargs)`
        :rtype: tuple
        """
        ...