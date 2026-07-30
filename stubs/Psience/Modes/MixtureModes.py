"""
Provides support for handling modes that arise from
"""
import enum
import numpy as np, scipy.linalg as slag, collections
from McUtils.Data import AtomData, UnitsData
import numpy as np, scipy
import McUtils.Numputils as nput
from McUtils.Coordinerds import CoordinateSystem, CartesianCoordinateSystem3D, InternalCoordinateSystem
__all__ = ['MixtureModes']

class MixtureModes(CoordinateSystem):
    """
    A `McUtils.Coordinerds.CoordinateSystem` object that expresses coordinates as
    a rotation on some base set of coordinates with some associated frequencies.
    """
    name = 'MixtureModes'

    def __init__(self, basis, coeffs, freqs=None, origin=None, masses=None, inverse=None, mass_weighted=False, frequency_scaled=False, g_matrix=None, name=None):
        """
        **LLM Docstring**

        Build a `MixtureModes` object -- a `CoordinateSystem` expressing coordinates as a (possibly non-orthogonal) linear rotation of some base coordinate system, with associated frequencies, masses, and mass-weighting/frequency-scaling/G-matrix bookkeeping.

        :param basis: the base coordinate system the modes are expressed in
        :type basis: CoordinateSystem
        :param coeffs: the mode coefficient (coordinates-by-modes) matrix, or a list of such matrices/tensors representing higher-order terms of a nonlinear coordinate expansion (the first entry is used as the linear matrix)
        :type coeffs: np.ndarray | list[np.ndarray]
        :param freqs: the vibrational frequencies associated with each mode
        :type freqs: np.ndarray | None
        :param origin: the reference geometry the modes are expanded about
        :type origin: np.ndarray | None
        :param masses: the atomic masses
        :type masses: np.ndarray | None
        :param inverse: the modes-by-coordinates (inverse) matrix
        :type inverse: np.ndarray | None
        :param mass_weighted: whether the modes are mass-weighted
        :type mass_weighted: bool
        :param frequency_scaled: whether the modes are frequency-scaled (dimensionless)
        :type frequency_scaled: bool
        :param g_matrix: the G-matrix associated with the modes
        :type g_matrix: np.ndarray | None
        :param name: a name for the mode set; defaults to the class's `name` attribute
        :type name: str | None
        :return: None
        :rtype: None
        """
        ...

    def to_state(self, serializer=None):
        """
        **LLM Docstring**

        Serialize this mode set's essential data (basis, matrix, inverse, frequencies, masses, mass-weighting/frequency-scaling flags, G-matrix) into a plain dict.

        :param serializer: the serializer used to recursively serialize the `basis` object
        :type serializer: object
        :return: the serialized state dict
        :rtype: dict
        """
        ...

    @classmethod
    def from_state(cls, data, serializer=None):
        """
        **LLM Docstring**

        Reconstruct a `MixtureModes` from a previously serialized state dict.

        :param data: the serialized state, as produced by `to_state`
        :type data: dict
        :param serializer: the serializer used to recursively deserialize the `basis` object
        :type serializer: object
        :return: the reconstructed mode set
        :rtype: MixtureModes
        """
        ...

    @classmethod
    def prep_modes(cls, modes):
        """
        **LLM Docstring**

        Coerce a variety of mode representations (an already-built `MixtureModes`-like object, a dict of matrix/inverse/basis/freqs, an object exposing `.matrix`, or a raw coefficient array) into a proper `MixtureModes` instance, inferring a sensible default basis (Cartesian or generic internal) if none is given.

        :param modes: the mode data to coerce
        :type modes: object | dict | np.ndarray
        :return: the constructed (or passed-through) mode set
        :rtype: MixtureModes
        :raises ValueError: if neither a matrix nor an inverse can be determined from `modes`
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

    def modify(self, matrix=None, *, freqs=None, origin=None, masses=None, inverse=None, name=None, mass_weighted=None, frequency_scaled=None, g_matrix=None):
        """
        **LLM Docstring**

        Build a new `MixtureModes` (of the same concrete subclass) with the given fields overridden, defaulting each unspecified field to this object's own current value.

        :param matrix: replacement modes-by-coordinates matrix; defaults to `self.matrix`
        :type matrix: np.ndarray | None
        :param freqs: replacement frequencies; defaults to `self.freqs`
        :type freqs: np.ndarray | None
        :param origin: replacement reference geometry; defaults to `self.origin`
        :type origin: np.ndarray | None
        :param masses: replacement masses; defaults to `self.masses`
        :type masses: np.ndarray | None
        :param inverse: replacement coords-by-modes matrix; defaults to `self.inverse`
        :type inverse: np.ndarray | None
        :param name: replacement name; defaults to `self.name`
        :type name: str | None
        :param mass_weighted: replacement mass-weighting flag; defaults to `self.mass_weighted`
        :type mass_weighted: bool | None
        :param frequency_scaled: replacement frequency-scaling flag; defaults to `self.frequency_scaled`
        :type frequency_scaled: bool | None
        :param g_matrix: replacement G-matrix; defaults to `self.g_matrix`
        :type g_matrix: np.ndarray | None
        :return: the new, modified mode set
        :rtype: MixtureModes
        """
        ...

    def rotate(self, rot, in_place=False):
        """
        **LLM Docstring**

        Intended to rotate the mode set by a given rotation. Not implemented -- the intended semantics were judged too ambiguous.

        :param rot: the rotation to apply
        :type rot: object
        :param in_place: whether to apply the rotation in place
        :type in_place: bool
        :return: never returns
        :rtype: None
        :raises NotImplementedError: always
        """
        ...

    def transform(self, tf, inv=None, origin=None):
        """
        **LLM Docstring**

        Intended to apply a linear transformation (and its inverse) to the mode set. Not implemented -- immediately raises before reaching the (dead) implementation code below it, which is retained but unreachable.

        :param tf: the transformation matrix
        :type tf: np.ndarray
        :param inv: the inverse transformation; computed from `tf` if `tf` is square and `inv` isn't given (in the unreachable code)
        :type inv: np.ndarray | None
        :param origin: the reference geometry to transform; defaults to `self.origin` (in the unreachable code)
        :type origin: np.ndarray | None
        :return: never returns
        :rtype: None
        :raises NotImplementedError: always, before the (dead) implementation runs
        :raises ValueError: (in the unreachable code) if `inv` isn't given and `tf` isn't square
        """
        ...

    @property
    def cartesian_modes(self):
        """
        **LLM Docstring**

        Whether the modes are expressed relative to a Cartesian (as opposed to internal-coordinate) origin, inferred from the origin's dimensionality.

        :return: `True` if the origin is 2-dimensional (`(natoms, 3)`)
        :rtype: bool
        """
        ...

    def embed_coords(self, carts):
        """
        **LLM Docstring**

        Convert a batch of Cartesian displacement structures into mode coordinates, by subtracting the reference origin and projecting through the coords-by-modes (inverse) matrix.

        :param carts: the Cartesian structures to embed
        :type carts: np.ndarray
        :return: the mode coordinates
        :rtype: np.ndarray
        """
        ...

    def unembed_coords(self, mode_coords):
        """
        **LLM Docstring**

        Convert a batch of mode coordinates back into Cartesian structures, by projecting through the modes-by-coordinates matrix and adding back the reference origin.

        :param mode_coords: the mode coordinates to unembed
        :type mode_coords: np.ndarray
        :return: the Cartesian structures
        :rtype: np.ndarray
        """
        ...

    @property
    def total_transformation(self):
        """
        **LLM Docstring**

        The full (possibly multi-term, nonlinear) forward coordinate-expansion tensors this mode set was constructed with.

        :return: the list of expansion-order tensors, starting with the linear modes-by-coordinates matrix
        :rtype: list[np.ndarray]
        """
        ...

    @property
    def inverse_transformation(self):
        """
        **LLM Docstring**

        The (lazily computed and cached) inverse of `total_transformation`, i.e. the Taylor-series expansion mapping mode coordinates back to base coordinates.

        :return: the inverse expansion-order tensors
        :rtype: list[np.ndarray]
        """
        ...

    def embed_derivs(self, derivs):
        """
        **LLM Docstring**

        Re-express a set of derivative tensors (with respect to the base coordinates) in terms of mode coordinates, by re-expanding through `total_transformation`.

        :param derivs: the derivative tensors to re-express
        :type derivs: list[np.ndarray]
        :return: the re-expressed derivative tensors
        :rtype: list[np.ndarray]
        """
        ...

    def unembed_derivs(self, derivs):
        """
        **LLM Docstring**

        Re-express a set of derivative tensors (with respect to mode coordinates) back in terms of the base coordinates, by re-expanding through `inverse_transformation`.

        :param derivs: the derivative tensors to re-express
        :type derivs: list[np.ndarray]
        :return: the re-expressed derivative tensors
        :rtype: list[np.ndarray]
        """
        ...

    @property
    def is_cartesian(self):
        """
        **LLM Docstring**

        Whether the modes are expressed over a Cartesian coordinate space, either inferred from the mode-matrix row count matching `3 * len(masses)` (if masses are known) or from the base coordinate system's name.

        :return: whether the modes are Cartesian-basis
        :rtype: bool
        """
        ...

    @property
    def coords_by_modes(self):
        """
        **LLM Docstring**

        The coordinates-by-modes (inverse) transformation matrix.

        :return: `self.inverse`
        :rtype: np.ndarray
        """
        ...

    @property
    def modes_by_coords(self):
        """
        **LLM Docstring**

        The modes-by-coordinates transformation matrix.

        :return: `self.matrix`
        :rtype: np.ndarray
        """
        ...

    def _eval_G(self, masses):
        """
        **LLM Docstring**

        Compute the coordinate-space G-matrix (inverse effective-mass metric) for this mode set: for Cartesian modes, the diagonal inverse-mass matrix; for internal-coordinate modes, derived from the (already mass-weighted) inverse mode matrix.

        :param masses: the atomic masses to use (only used in the Cartesian branch)
        :type masses: np.ndarray
        :return: the G-matrix
        :rtype: np.ndarray
        :raises NotImplementedError: for non-mass-weighted internal-coordinate modes, which aren't supported
        """
        ...

    def _get_gmatrix(self, masses=None):
        """
        **LLM Docstring**

        Resolve (and cache, if computed fresh for the default masses) the G-matrix and its fractional powers, computing it via `_eval_G` if not already cached and no alternate `masses` are given.

        :param masses: alternate masses to compute a fresh (uncached) G-matrix for; if `None`, uses (and caches into) `self.g_matrix`/`self.masses`
        :type masses: np.ndarray | None
        :return: `(masses, g12, gi12)` -- the masses used, and the G-matrix's square root and inverse square root
        :rtype: tuple[np.ndarray, np.ndarray, np.ndarray]
        """
        ...

    @classmethod
    def _local_tf(cls, f, g):
        """
        **LLM Docstring**

        Compute the diagonal rescaling transformation that converts between a set of diagonal force constants `f` and G-matrix elements `g`, used to build "local mode" (Duschinsky-diagonal) representations.

        :param f: the diagonal force-constant values
        :type f: np.ndarray
        :param g: the diagonal G-matrix values
        :type g: np.ndarray
        :return: the diagonal rescaling matrix
        :rtype: np.ndarray
        """
        ...

    @classmethod
    def compute_local_transformations(cls, f, g):
        """
        **LLM Docstring**

        Compute the pair of diagonal rescaling transformations (for the force-constant matrix and its inverse-G counterpart) that define the "local mode" representation from diagonal `f`/`g` values.

        :param f: the diagonal force-constant values
        :type f: np.ndarray
        :param g: the diagonal G-matrix values
        :type g: np.ndarray
        :return: `[local_tf(f, g), local_tf(g, f)]`
        :rtype: list[np.ndarray]
        """
        ...

    @classmethod
    def compute_local_hessian(cls, f, g):
        """
        **LLM Docstring**

        Compute the rescaled ("local mode") force-constant matrix from full `f`/`g` matrices, using the diagonal rescaling from `_local_tf`.

        :param f: the force-constant (Hessian) matrix
        :type f: np.ndarray
        :param g: the G-matrix
        :type g: np.ndarray
        :return: the rescaled local Hessian
        :rtype: np.ndarray
        """
        ...

    @classmethod
    def compute_local_gmatrix(cls, f, g):
        """
        **LLM Docstring**

        Compute the rescaled ("local mode") G-matrix from full `f`/`g` matrices, using the diagonal rescaling from `_local_tf`.

        :param f: the force-constant (Hessian) matrix
        :type f: np.ndarray
        :param g: the G-matrix
        :type g: np.ndarray
        :return: the rescaled local G-matrix
        :rtype: np.ndarray
        """
        ...

    def compute_hessian(self, system='modes'):
        """
        **LLM Docstring**

        Compute the (signed-frequency-squared) Hessian matrix in either the mode basis or the underlying coordinate basis, by reconstructing it from the stored frequencies and the mode transformation.

        :param system: which basis to compute the Hessian in, `'modes'` or `'coords'`
        :type system: str
        :return: the Hessian matrix
        :rtype: np.ndarray
        :raises ValueError: if `system` is neither `'modes'` nor `'coords'`
        """
        ...

    def compute_gmatrix(self, system='modes', return_fractional=False):
        """
        **LLM Docstring**

        Compute the G-matrix in either the mode basis or the underlying coordinate basis, using whichever information is available (mass-weighting, a stored G-matrix, or Cartesian masses).

        :param system: which basis to compute the G-matrix in, `'modes'` or `'coords'`
        :type system: str
        :param return_fractional: whether to also return the G-matrix's square root and inverse square root
        :type return_fractional: bool
        :return: the G-matrix, or `(g, g12, gi12)` if `return_fractional` is set; `None` if no G-matrix information is available in the `'modes'` case
        :rtype: np.ndarray | tuple | None
        :raises ValueError: if `system` is neither `'modes'` nor `'coords'`
        """
        ...
    zero_freq_cutoff = 1e-08

    @classmethod
    def _zero_projected_modes(cls, f, g):
        """
        **LLM Docstring**

        Solve the generalized eigenvalue problem `f v = freq^2 g v`, falling back to first projecting out the null space of `g` (via its own eigendecomposition) if the direct generalized eigensolve fails due to `g` being singular.

        :param f: the force-constant matrix
        :type f: np.ndarray
        :param g: the G-matrix
        :type g: np.ndarray
        :return: `(freqs, modes)` -- the signed square-root frequencies and the corresponding eigenvectors
        :rtype: tuple[np.ndarray, np.ndarray]
        """
        ...

    def compute_freqs(self):
        """
        **LLM Docstring**

        Recompute the vibrational frequencies for this mode set from its coordinate-basis Hessian and G-matrix, via a generalized eigenvalue solve.

        :return: the signed square-root frequencies
        :rtype: np.ndarray
        """
        ...

    @property
    def local_hessian(self):
        """
        **LLM Docstring**

        The "local mode" (diagonally rescaled) force-constant matrix for this mode set.

        :return: the local Hessian
        :rtype: np.ndarray
        """
        ...

    @property
    def local_gmatrix(self):
        """
        **LLM Docstring**

        The "local mode" (diagonally rescaled) G-matrix for this mode set.

        :return: the local G-matrix
        :rtype: np.ndarray
        """
        ...

    @property
    def local_freqs(self):
        """
        **LLM Docstring**

        The diagonal entries of the local-mode Hessian, giving an approximate per-mode "local" force constant/frequency.

        :return: the local-mode diagonal values
        :rtype: np.ndarray
        """
        ...

    @property
    def local_mode_transformations(self):
        """
        **LLM Docstring**

        The pair of diagonal rescaling transformations mapping between this mode set's coordinate-basis Hessian/G-matrix and their local-mode counterparts.

        :return: `[hessian_scaling, gmatrix_scaling]`
        :rtype: list[np.ndarray]
        """
        ...

    def get_nearest_mode_transform(self, alternate_modes: np.ndarray, mass_weighted=None, atoms=None, maximum_similarity=True, unitarize=True, masses=None):
        """
        **LLM Docstring**

        Find the transformation (and its inverse) that best aligns this mode set with an externally supplied set of `alternate_modes`, either via a maximum-similarity (orthogonal Procrustes-style) alignment or via direct projection (optionally unitarized), with optional mass-(un)weighting and atom-restriction beforehand.

        :param alternate_modes: the external mode matrix to align/project against
        :type alternate_modes: np.ndarray
        :param mass_weighted: if given, coerces this mode set to (or from) mass-weighted form before comparing
        :type mass_weighted: bool | None
        :param atoms: restrict the comparison to a subset of atoms via an atom projector (Cartesian modes only)
        :type atoms: Iterable[int] | None
        :param maximum_similarity: whether to use a maximum-similarity (orthogonal Procrustes) alignment rather than direct projection
        :type maximum_similarity: bool
        :param unitarize: whether to unitarize the resulting transformation (direct-projection path only)
        :type unitarize: bool
        :param masses: masses to use for the mass-(un)weighting conversion
        :type masses: np.ndarray | None
        :return: `(tf, inv)` -- the forward and inverse transformation matrices
        :rtype: tuple[np.ndarray, np.ndarray]
        """
        ...
    ModeData = collections.namedtuple('ModeData', ['freqs', 'modes', 'inverse'])
    default_zero_freq_cutoff = 0.0001

    @classmethod
    def get_normal_modes(cls, f_matrix, mass_spec, remove_transrot=True, dimensionless=False, mass_weighted=None, zero_freq_cutoff=None, return_gmatrix=False, projector=None):
        """
        **LLM Docstring**

        Compute ordinary harmonic normal modes from a force-constant matrix and mass specification (masses, element symbols, or a full G-matrix), via a generalized eigenvalue solve on the mass-weighted Hessian, optionally removing near-zero-frequency (translation/rotation) modes, applying a custom constraint projector, mass-weighting and/or frequency-scaling (dimensionless-izing) the result, and supporting batched inputs (including "ragged" batches where a different number of modes survives per entry).

        :param f_matrix: the force-constant (Hessian) matrix/matrices
        :type f_matrix: np.ndarray
        :param mass_spec: the atomic masses (or element symbols), or an already-built G-matrix/mass matrix
        :type mass_spec: np.ndarray | list[str]
        :param remove_transrot: whether to discard near-zero-frequency modes
        :type remove_transrot: bool
        :param dimensionless: whether the resulting modes should be frequency-scaled (dimensionless)
        :type dimensionless: bool
        :param mass_weighted: whether the resulting modes should be mass-weighted; defaults to `dimensionless` when `mass_spec` is a bare mass vector
        :type mass_weighted: bool | None
        :param zero_freq_cutoff: the frequency cutoff below which modes are discarded; defaults to `cls.default_zero_freq_cutoff`
        :type zero_freq_cutoff: float | None
        :param return_gmatrix: whether to also return the (broadcast) mass/G-matrix used
        :type return_gmatrix: bool
        :param projector: an additional constraint projector to apply to the mass-weighted Hessian before diagonalizing
        :type projector: np.ndarray | None
        :return: the mode data (`ModeData(freqs, modes, inverse)`), or `(mode_data, mass_spec)` if `return_gmatrix` is set
        :rtype: object | tuple
        """
        ...
    localization_type = 'ned'
    localization_zero_freq_cutoff = 99 / 219474.56

    def get_projected_localized_mode_transformation(self, projectors, masses=None, origin=None, localization_type=None, allow_mode_mixing=False, maximum_similarity=False, unitarize=True, zero_freq_cutoff=None, orthogonal_projection=False, project_zero_gmatrix_modes=True, project_zero_gmatrix_cutoff=1e-08, atoms=None):
        """
        **LLM Docstring**

        Shared driver behind the various coordinate/atom/fragment-localized-mode constructors: diagonalizes the (optionally translation/rotation-projected) mass-weighted Hessian restricted to (or blocked by) a set of projectors, either mixing all projected subspaces together (`allow_mode_mixing`) or diagonalizing each independently and concatenating the results, then finds the transformation aligning the current modes onto the resulting localized modes via `get_nearest_mode_transform`.

        :param projectors: the projection matrix (or matrices) defining the coordinate subspace(s) to localize within
        :type projectors: list[np.ndarray]
        :param masses: masses to use instead of `self.masses`
        :type masses: np.ndarray | None
        :param origin: reference geometry to use for the translation/rotation-invariant transformation (`'direct'` localization type only)
        :type origin: np.ndarray | None
        :param localization_type: `'ned'`/other for the standard projected-Hessian approach, or `'direct'` to first factor out translation/rotation via an explicit invariant transformation (Cartesian modes only)
        :type localization_type: str | None
        :param allow_mode_mixing: whether to mix all projected subspaces together into a single diagonalization rather than treating each independently
        :type allow_mode_mixing: bool
        :param maximum_similarity: forwarded to `get_nearest_mode_transform`
        :type maximum_similarity: bool
        :param unitarize: forwarded to `get_nearest_mode_transform`
        :type unitarize: bool
        :param zero_freq_cutoff: the frequency cutoff used when diagonalizing each subspace; defaults to `self.localization_zero_freq_cutoff`
        :type zero_freq_cutoff: float | None
        :param orthogonal_projection: whether non-square entries in `projectors` should be treated as bases for an orthogonal (rather than oblique) projection
        :type orthogonal_projection: bool
        :param project_zero_gmatrix_modes: accepted for interface consistency but not used in this method's body
        :type project_zero_gmatrix_modes: bool
        :param project_zero_gmatrix_cutoff: accepted for interface consistency but not used in this method's body
        :type project_zero_gmatrix_cutoff: float
        :param atoms: accepted for interface consistency with sibling localizers but not used in this method's body
        :type atoms: object | None
        :return: `(tf, inv)`, the localizing transformation and its inverse
        :rtype: tuple[np.ndarray, np.ndarray]
        """
        ...

    def get_atom_localized_mode_transformation(self, atoms, masses=None, origin=None, localization_type='ned', allow_mode_mixing=False, maximum_similarity=False, orthogonal_projection=False, unitarize=True, zero_freq_cutoff=None):
        """
        **LLM Docstring**

        Build a localized-mode transformation that concentrates vibrational character onto (or, if `orthogonal_projection`, away from) a specified set of atoms, via `get_projected_localized_mode_transformation` with per-atom (or combined) atom-selection projectors.

        :param atoms: the atom index (or indices) to localize modes onto
        :type atoms: int | Iterable[int]
        :param masses: masses to use instead of `self.masses`
        :type masses: np.ndarray | None
        :param origin: reference geometry, forwarded to `get_projected_localized_mode_transformation`
        :type origin: np.ndarray | None
        :param localization_type: the localization scheme to use
        :type localization_type: str
        :param allow_mode_mixing: whether to build one combined projector spanning all `atoms` rather than one per atom
        :type allow_mode_mixing: bool
        :param maximum_similarity: forwarded to `get_projected_localized_mode_transformation`
        :type maximum_similarity: bool
        :param orthogonal_projection: whether to localize onto the complement of `atoms` (all other atoms) instead of `atoms` itself
        :type orthogonal_projection: bool
        :param unitarize: forwarded to `get_projected_localized_mode_transformation`
        :type unitarize: bool
        :param zero_freq_cutoff: forwarded to `get_projected_localized_mode_transformation`
        :type zero_freq_cutoff: float | None
        :return: `(tf, inv)`, the localizing transformation and its inverse
        :rtype: tuple[np.ndarray, np.ndarray]
        """
        ...

    def get_fragment_localized_mode_transformation(self, fragment, masses=None, origin=None, localization_type='ned', allow_mode_mixing=True, maximum_similarity=False, orthogonal_projection=False, unitarize=True, zero_freq_cutoff=None, **etc):
        """
        **LLM Docstring**

        Build a localized-mode transformation for a molecular fragment (a set of atoms), via `get_atom_localized_mode_transformation`.

        :param fragment: the atom indices making up the fragment
        :type fragment: Iterable[int]
        :param masses: masses to use instead of `self.masses`
        :type masses: np.ndarray | None
        :param origin: reference geometry
        :type origin: np.ndarray | None
        :param localization_type: the localization scheme to use
        :type localization_type: str
        :param allow_mode_mixing: whether to build one combined projector for the whole fragment
        :type allow_mode_mixing: bool
        :param maximum_similarity: forwarded through to the underlying localizer
        :type maximum_similarity: bool
        :param orthogonal_projection: whether to localize onto the complement of the fragment instead
        :type orthogonal_projection: bool
        :param unitarize: forwarded through to the underlying localizer
        :type unitarize: bool
        :param zero_freq_cutoff: forwarded through to the underlying localizer
        :type zero_freq_cutoff: float | None
        :param etc: extra options forwarded to `get_atom_localized_mode_transformation`
        :type etc: dict
        :return: `(tf, inv)`, the localizing transformation and its inverse
        :rtype: tuple[np.ndarray, np.ndarray]
        """
        ...

    def get_coordinate_projected_localized_mode_transformation(self, coordinate_constraints, atoms=None, masses=None, origin=None, localization_type='ned', allow_mode_mixing=False, maximum_similarity=False, orthogonal_projection=True, unitarize=True):
        """
        **LLM Docstring**

        Build a localized-mode transformation that concentrates vibrational character along a set of internal-coordinate directions (bond/angle/dihedral bases for Cartesian modes, or by coordinate index for internal-coordinate modes), optionally restricted to a subset of atoms, via `get_projected_localized_mode_transformation`.

        :param coordinate_constraints: the internal-coordinate specification(s) (for Cartesian modes) or coordinate index/indices (for internal-coordinate modes) to localize along
        :type coordinate_constraints: object | Iterable[int]
        :param atoms: restrict the resulting projector(s) to this subset of atoms (Cartesian modes only)
        :type atoms: Iterable[int] | None
        :param masses: masses to use instead of `self.masses`
        :type masses: np.ndarray | None
        :param origin: reference geometry, forwarded to `get_projected_localized_mode_transformation`
        :type origin: np.ndarray | None
        :param localization_type: the localization scheme to use
        :type localization_type: str
        :param allow_mode_mixing: whether to combine all coordinate projectors into one before localizing
        :type allow_mode_mixing: bool
        :param maximum_similarity: forwarded to `get_projected_localized_mode_transformation`
        :type maximum_similarity: bool
        :param orthogonal_projection: whether the coordinate bases define orthogonal (rather than oblique) projections
        :type orthogonal_projection: bool
        :param unitarize: forwarded to `get_projected_localized_mode_transformation`
        :type unitarize: bool
        :return: `(tf, inv)`, the localizing transformation and its inverse
        :rtype: tuple[np.ndarray, np.ndarray]
        :raises ValueError: if `coordinate_constraints` isn't recognized as coordinate indices for internal-coordinate modes
        """
        ...

    def get_internal_localized_mode_transformation(self, expansion_coordinates: 'Iterable[Iterable[int]|dict]', fixed_atoms=None, mass_weighted=False, project_transrot=True, atoms=None, maximum_similarity=False, orthogonal_projection=False, projection=False, allow_mode_mixing=False, unitarize=True, origin=None, masses=None, localization_type='ned'):
        """
        **LLM Docstring**

        Build a localized-mode transformation aligned with a set of internal-coordinate displacement directions (computed at the reference geometry), either by projecting the Hessian onto (or away from) those directions and re-diagonalizing (`projection`/`orthogonal_projection`), or by directly finding the nearest-mode alignment to the raw coordinate-derivative directions (optionally mass-weighted and/or translation/rotation-projected).

        :param expansion_coordinates: the internal coordinate(s) whose Cartesian derivative directions define the localization target
        :type expansion_coordinates: Iterable[Iterable[int]] | dict
        :param fixed_atoms: atoms to hold fixed when computing the internal-coordinate derivative tensors
        :type fixed_atoms: Iterable[int] | None
        :param mass_weighted: whether to mass-weight the coordinate-derivative directions
        :type mass_weighted: bool
        :param project_transrot: whether to project out translational/rotational components from the derivative directions (direct-alignment path only)
        :type project_transrot: bool
        :param atoms: restrict the projector(s)/alignment to a subset of atoms
        :type atoms: Iterable[int] | None
        :param maximum_similarity: forwarded to the underlying localizer
        :type maximum_similarity: bool
        :param orthogonal_projection: whether to use orthogonal (rather than oblique) projections in the projection-based path
        :type orthogonal_projection: bool
        :param projection: whether to use the projection-and-rediagonalize approach rather than direct nearest-mode alignment
        :type projection: bool
        :param allow_mode_mixing: whether to combine all coordinate directions into a single projected subspace (projection-based path)
        :type allow_mode_mixing: bool
        :param unitarize: forwarded to the underlying localizer
        :type unitarize: bool
        :param origin: reference geometry to use
        :type origin: np.ndarray | None
        :param masses: masses to use instead of `self.masses`
        :type masses: np.ndarray | None
        :param localization_type: the localization scheme to use for the projection-based path
        :type localization_type: str
        :return: `(tf, inv)`, the localizing transformation and its inverse
        :rtype: tuple[np.ndarray, np.ndarray]
        """
        ...

    def get_displacement_localized_mode_transformation(self, mode_blocks=None, atoms=None, mass_weighted=True, unitarize=True, **maximizer_opts):
        """
        **LLM Docstring**

        Build a unitary localized-mode transformation by iteratively (Jacobi-style) rotating blocks of modes to maximize a displacement-localization criterion (via `nput.jacobi_maximize`/`nput.displacement_localizing_rotation_generator`), optionally restricted to a subset of atoms and/or applied independently within specified mode blocks.

        :param mode_blocks: groups of mode indices to localize independently; a single flat list of indices localizes all modes together, and `None` localizes all modes as one block
        :type mode_blocks: Iterable[int] | list[Iterable[int]] | None
        :param atoms: restrict the localization criterion to a subset of atoms
        :type atoms: Iterable[int] | None
        :param mass_weighted: whether to localize the mass-weighted (rather than un-mass-weighted) mode matrix
        :type mass_weighted: bool
        :param unitarize: must be `True`; the Jacobi-rotation scheme only produces unitary transformations
        :type unitarize: bool
        :param maximizer_opts: extra options forwarded to `nput.jacobi_maximize`
        :type maximizer_opts: dict
        :return: `(tf_inv, tf_inv.T)`, the (unitary) localizing transformation's inverse and its transpose (used as the forward transformation)
        :rtype: tuple[np.ndarray, np.ndarray]
        :raises ValueError: if `unitarize` is `False`
        """
        ...

    def get_mass_scaled_mode_transformation(self, mass_scaling, *, atoms, localization_cutoff=0.8, num_modes=None, project_transrot=False, unitarize=True, **diag_opts):
        """
        **LLM Docstring**

        Build a localized-mode transformation by artificially scaling the masses of a set of atoms (making them effectively heavier or lighter) and re-diagonalizing, then selecting the resulting modes that are most concentrated on those atoms (by projected-displacement norm), optionally applying a localization cutoff or a fixed mode count.

        :param mass_scaling: the scale factor (or per-atom scale factors) to apply to the selected atoms' masses
        :type mass_scaling: float | np.ndarray
        :param atoms: the atom(s) whose masses should be scaled
        :type atoms: int | Iterable[int]
        :param localization_cutoff: minimum normalized displacement-on-atoms required for a mode to be selected; if `None`, the top `num_modes` (or `3*len(atoms)`) modes by that metric are taken instead
        :type localization_cutoff: float | None
        :param num_modes: the number of modes to select; defaults to `3 * len(atoms)`
        :type num_modes: int | None
        :param project_transrot: whether to project out translation/rotation using the mass-scaled masses before diagonalizing
        :type project_transrot: bool
        :param unitarize: accepted for interface consistency with sibling localizers but not used in this method's body
        :type unitarize: bool
        :param diag_opts: extra options forwarded to `NormalModes.get_normal_modes`
        :type diag_opts: dict
        :return: `(tf, tf.T)`, the selected mode transformation and its transpose, or `(None, None)` if no modes clear `localization_cutoff`
        :rtype: tuple[np.ndarray, np.ndarray] | tuple[None, None]
        """
        ...

    class LocalizationMethods(enum.Enum):
        """Real access pattern: LocalizationMethods.<MemberName> (this is an enum with 8 members, e.g. LocalizationMethods.MaximumSimilarity == 'target_modes'). Collapsed into a dict below purely for compactness -- do not index it as a dict in real code:"""
        _MEMBERS = {'MaximumSimilarity': 'target_modes', 'AtomLocalized': 'atoms', 'FragmentLocalized': 'fragment', 'DisplacmentMinimized': 'displacements', 'Internals': 'coordinates', 'CoordinateConstraints': 'constraints', 'Projected': 'projections', 'MassScaled': 'mass_scaling'}

    @property
    def localizer_dispatch(self):
        """
        **LLM Docstring**

        The mapping from `LocalizationMethods` value to the `(constructor_method, primary_argument_name)` pair used by `localize` to dispatch a localization request to the right underlying method.

        :return: the method-name-to-`(constructor, arg_name)` mapping
        :rtype: dict
        """
        ...

    def localize(self, method=None, *, atoms=None, fragment=None, target_modes=None, internals=None, mode_blocks=None, coordinate_constraints=None, projections=None, reorthogonalize=None, mass_scaling=None, unitarize=True, allow_mode_mixing=False, project_zero_gmatrix_modes=None, project_zero_gmatrix_cutoff=1e-08, **opts):
        """
        **LLM Docstring**

        Top-level entry point for building a `LocalizedModes` object: infers which localization method to use from whichever keyword argument was supplied (if `method` isn't given explicitly), dispatches to the corresponding `get_*_localized_mode_transformation` method, optionally projects out (and, if mixing is allowed, re-diagonalizes within) any modes with a near-zero effective G-matrix eigenvalue, optionally re-orthogonalizes the result, and wraps the resulting transformation in a `LocalizedModes` object.

        :param method: which localization method to use; inferred from the other arguments if not given
        :type method: str | LocalizationMethods | None
        :param atoms: atom(s) to localize onto, for the `'atoms'` method (or as a restriction for others)
        :type atoms: int | Iterable[int] | None
        :param fragment: fragment atoms to localize onto, for the `'fragment'` method
        :type fragment: Iterable[int] | None
        :param target_modes: external modes to align to, for the `'target_modes'` method
        :type target_modes: np.ndarray | None
        :param internals: internal coordinates to localize along, for the `'coordinates'` method
        :type internals: object | None
        :param mode_blocks: mode index groupings, for the `'displacements'` method
        :type mode_blocks: object | None
        :param coordinate_constraints: coordinate constraints, for the `'constraints'` method
        :type coordinate_constraints: object | None
        :param projections: explicit projectors, for the `'projections'` method
        :type projections: list[np.ndarray] | None
        :param reorthogonalize: whether to re-orthogonalize the localized modes' mass-weighted representation after localization; defaults to `not unitarize`
        :type reorthogonalize: bool | None
        :param mass_scaling: mass-scaling factor(s), for the `'mass_scaling'` method
        :type mass_scaling: float | np.ndarray | None
        :param unitarize: whether the underlying localizer should produce a unitary transformation
        :type unitarize: bool
        :param allow_mode_mixing: whether the underlying localizer is allowed to mix modes across different projected subspaces
        :type allow_mode_mixing: bool
        :param project_zero_gmatrix_modes: whether to project out (and handle) modes with a near-zero effective G-matrix eigenvalue after localizing; defaults to `allow_mode_mixing`
        :type project_zero_gmatrix_modes: bool | None
        :param project_zero_gmatrix_cutoff: the G-matrix eigenvalue magnitude below which a mode is treated as zero
        :type project_zero_gmatrix_cutoff: float
        :param opts: extra options forwarded to the dispatched localizer method
        :type opts: dict
        :return: the constructed `LocalizedModes` object, or `None` if the underlying localizer returned no transformation
        :rtype: LocalizedModes | None
        """
        ...

    def make_mass_weighted(self, masses=None):
        """
        **LLM Docstring**

        Build a mass-weighted version of this mode set by rescaling the mode/inverse matrices and origin through the G-matrix's square root/inverse-square-root; returns `self` unchanged if already mass-weighted.

        :param masses: masses to use instead of `self.masses`
        :type masses: np.ndarray | None
        :return: the mass-weighted mode set
        :rtype: MixtureModes
        """
        ...

    def remove_mass_weighting(self, masses=None):
        """
        **LLM Docstring**

        Build a non-mass-weighted version of this mode set by undoing the mass-weighting rescaling; returns `self` unchanged if not already mass-weighted.

        :param masses: masses to use instead of `self.masses`
        :type masses: np.ndarray | None
        :return: the non-mass-weighted mode set
        :rtype: MixtureModes
        """
        ...

    def _frequency_scaling(self, freqs=None):
        """
        **LLM Docstring**

        Compute the per-mode scaling factor (`sqrt(|freq|)`) used to convert between frequency-scaled and non-frequency-scaled mode bases.

        :param freqs: the frequencies to compute the scaling from; defaults to `self.freqs`
        :type freqs: np.ndarray | None
        :return: `(freqs, conv)` -- the frequencies used and their corresponding scaling factors
        :rtype: tuple[np.ndarray, np.ndarray]
        """
        ...

    def make_frequency_scaled(self, freqs=None):
        """
        **LLM Docstring**

        Build a frequency-scaled (dimensionless) version of this mode set by rescaling the mode/inverse matrices by the per-mode frequency-scaling factor; returns `self` unchanged if already frequency-scaled.

        :param freqs: the frequencies to compute the scaling from; defaults to `self.freqs`
        :type freqs: np.ndarray | None
        :return: the frequency-scaled mode set
        :rtype: MixtureModes
        """
        ...

    def remove_frequency_scaling(self, freqs=None):
        """
        **LLM Docstring**

        Build a non-frequency-scaled version of this mode set by undoing the frequency-scaling rescaling; returns `self` unchanged if not already frequency-scaled.

        :param freqs: the frequencies to compute the scaling from; defaults to `self.freqs`
        :type freqs: np.ndarray | None
        :return: the non-frequency-scaled mode set
        :rtype: MixtureModes
        """
        ...

    def make_dimensionless(self, freqs=None, masses=None):
        """
        **LLM Docstring**

        Build a fully dimensionless version of this mode set: mass-weighted and frequency-scaled.

        :param freqs: the frequencies to compute the frequency scaling from
        :type freqs: np.ndarray | None
        :param masses: the masses to compute the mass-weighting from
        :type masses: np.ndarray | None
        :return: the dimensionless mode set
        :rtype: MixtureModes
        """
        ...

    def make_dimensioned(self, freqs=None, masses=None):
        """
        **LLM Docstring**

        Build a fully dimensioned version of this mode set: not mass-weighted and not frequency-scaled.

        :param freqs: the frequencies to compute the frequency-descaling from
        :type freqs: np.ndarray | None
        :param masses: the masses to compute the mass-unweighting from
        :type masses: np.ndarray | None
        :return: the dimensioned mode set
        :rtype: MixtureModes
        """
        ...

    def apply_projection(self, proj, project_transrot=True, masses=None, origin=None):
        """
        **LLM Docstring**

        Apply an arbitrary projection matrix to the mode/inverse matrices (optionally first combining it with a translation/rotation projector), returning a new mode set restricted to (or excluding) the projected subspace.

        :param proj: the projection matrix to apply
        :type proj: np.ndarray
        :param project_transrot: whether to additionally combine `proj` with a translation/rotation projector before applying
        :type project_transrot: bool
        :param masses: masses to use for the translation/rotation projector and to store on the result
        :type masses: np.ndarray | None
        :param origin: reference geometry to use for the translation/rotation projector and to store on the result
        :type origin: np.ndarray | None
        :return: the projected mode set
        :rtype: MixtureModes
        """
        ...

    @classmethod
    def _atom_projector(cls, n, i, orthogonal_projection=False):
        """
        **LLM Docstring**

        Build a diagonal `3n x 3n` projection matrix selecting (or, if `orthogonal_projection`, excluding) the Cartesian coordinates belonging to a given set of atom indices.

        :param n: the total number of atoms
        :type n: int
        :param i: the atom index (or indices) to select
        :type i: int | Iterable[int]
        :param orthogonal_projection: whether to build the complementary projector (excluding `i` instead of selecting it)
        :type orthogonal_projection: bool
        :return: the diagonal atom-selection projection matrix
        :rtype: np.ndarray
        """
        ...

    def apply_constraints(self, coordinate_constraints, atoms=None, masses=None, origin=None, orthogonal_projection=True):
        """
        **LLM Docstring**

        Build a new mode set with the given internal coordinates constrained (projected out or onto, depending on `orthogonal_projection`), optionally restricted to a subset of atoms, via `apply_projection`.

        :param coordinate_constraints: the internal coordinate specification(s) to constrain
        :type coordinate_constraints: object
        :param atoms: restrict the constraint projector(s) to a subset of atoms
        :type atoms: Iterable[int] | None
        :param masses: masses to use instead of `self.masses`
        :type masses: np.ndarray | None
        :param origin: reference geometry to use instead of `self.origin` (un-mass-weighted)
        :type origin: np.ndarray | None
        :param orthogonal_projection: whether the constraint bases define orthogonal (rather than oblique) projections, and whether multiple constraints are combined by sequential orthogonal projection (`True`) or simple summation (`False`)
        :type orthogonal_projection: bool
        :return: the constrained mode set
        :rtype: MixtureModes
        """
        ...

    def apply_transformation(self, tf, **opts):
        """
        **LLM Docstring**

        Apply an arbitrary linear transformation to this mode set, returning the result as a `LocalizedModes` object.

        :param tf: the transformation to apply, in any form accepted by `LocalizedModes`
        :type tf: np.ndarray | tuple
        :param opts: extra options forwarded to the `LocalizedModes` constructor
        :type opts: dict
        :return: the transformed modes
        :rtype: LocalizedModes
        """
        ...

    def make_oblique(self):
        """
        **LLM Docstring**

        Build an "oblique" mode representation (see `ObliqueModeGenerator`) that makes the mode transformation as close to orthogonal as possible while still diagonalizing the mode-basis Hessian/G-matrix.

        :return: the oblique-transformed modes
        :rtype: LocalizedModes
        """
        ...