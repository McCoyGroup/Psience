__all__ = ['MolecularHamiltonian']
import numpy as np, itertools
import McUtils.Numputils as nput
import McUtils.Iterators as itut
from McUtils.Combinatorics import UniquePermutations
from McUtils.Zachary import TensorDerivativeConverter
from ..Modes import MixtureModes
from .CoordinateSystems import MolecularEmbedding, ModeEmbedding
from .Properties import PotentialSurfaceManager, DipoleSurfaceManager, NormalModesManager

class MolecularHamiltonian:
    """
    Provides a generic molecular Hamiltonian with support for expansions (for VPT),
    construction of molecular models, and other wave function generators
    """

    def __init__(self, embedding: MolecularEmbedding, potential_manager: PotentialSurfaceManager=None, potential_derivatives: list[np.ndarray]=None, modes_manager: NormalModesManager=None, modes: MixtureModes=None, dipole_manager: DipoleSurfaceManager=None, dipole_derivatives: list[np.ndarray]=None):
        """
        **LLM Docstring**

        Assemble a molecular Hamiltonian from a coordinate embedding plus normal-mode, potential, and dipole information.

        :param embedding: the molecular embedding (Cartesian/internal coordinate specification) the Hamiltonian is built over
        :type embedding: MolecularEmbedding
        :param potential_manager: an existing potential-surface manager to wrap; used in preference to `potential_derivatives` if given
        :type potential_manager: PotentialSurfaceManager | None
        :param potential_derivatives: raw potential-energy derivative tensors, used if `potential_manager` is not given
        :type potential_derivatives: list[np.ndarray] | None
        :param modes_manager: an existing normal-modes manager to use; used in preference to `modes` if given
        :type modes_manager: NormalModesManager | None
        :param modes: raw normal modes, used if `modes_manager` is not given
        :type modes: MixtureModes | None
        :param dipole_manager: an existing dipole-surface manager to wrap; used in preference to `dipole_derivatives` if given
        :type dipole_manager: DipoleSurfaceManager | None
        :param dipole_derivatives: raw dipole-surface derivative tensors, used if `dipole_manager` is not given
        :type dipole_derivatives: list[np.ndarray] | None
        :return: None
        :rtype: None
        """
        ...

    @classmethod
    def _prep_modes(cls, modes_manager, modes):
        """
        **LLM Docstring**

        Pick which normal-modes source to use: prefer `modes_manager` if supplied, otherwise fall back to the raw `modes`.

        :param modes_manager: a normal-modes manager, if available
        :type modes_manager: NormalModesManager | None
        :param modes: raw normal modes, used if no manager is given
        :type modes: MixtureModes | None
        :return: whichever of the two arguments was provided
        :rtype: NormalModesManager | MixtureModes | None
        """
        ...

    @classmethod
    def _prep_potential(cls, potential_manager, potential_derivatives):
        """
        **LLM Docstring**

        Wrap the supplied potential source in a `ScalarOperatorManager`, preferring `potential_manager` over `potential_derivatives`.

        :param potential_manager: an existing potential-surface manager
        :type potential_manager: PotentialSurfaceManager | None
        :param potential_derivatives: raw potential derivative tensors, used if no manager is given
        :type potential_derivatives: list[np.ndarray] | None
        :return: a `ScalarOperatorManager` wrapping whichever source was provided
        :rtype: ScalarOperatorManager
        """
        ...

    @classmethod
    def _prep_dipole(cls, dipole_manager, dipole_derivatives):
        """
        **LLM Docstring**

        Wrap the supplied dipole source in a `ScalarOperatorManager`, preferring `dipole_manager` over `dipole_derivatives`.

        :param dipole_manager: an existing dipole-surface manager
        :type dipole_manager: DipoleSurfaceManager | None
        :param dipole_derivatives: raw dipole derivative tensors, used if no manager is given
        :type dipole_derivatives: list[np.ndarray] | None
        :return: a `ScalarOperatorManager` wrapping whichever source was provided
        :rtype: ScalarOperatorManager
        """
        ...

    def get_VPT_expansions(self, order=2, coordinate_transformation=None):
        """
        **LLM Docstring**

        Build the set of Taylor expansions (potential, kinetic/G-matrix, and either Coriolis + Watson or pseudopotential, depending on whether the embedding uses Cartesians or internals) needed for a VPT (vibrational perturbation theory) treatment, optionally re-expressed through an extra coordinate transformation.

        :param order: either a single integer giving the base expansion order (used to derive per-term orders, with Coriolis/Watson/pseudopotential going 2 orders lower), or an explicit dict mapping term name (`'potential'`, `'kinetic'`, `'coriolis'`, `'pseudopotential'`, `'watson'`) to its order
        :type order: int | dict
        :param coordinate_transformation: an additional forward/reverse coordinate transformation to apply to every expansion, passed through `prep_transformation`
        :type coordinate_transformation: tuple | None
        :return: mapping from term name to its expansion terms
        :rtype: dict
        :raises ValueError: if `order` is a dict containing a key that isn't one of the recognized term names
        """
        ...

    @classmethod
    def prep_transformation(cls, coordinate_transformation):
        """
        **LLM Docstring**

        Normalize a coordinate transformation into a `(forward, reverse)` pair of derivative-tensor lists, computing the reverse (inverse) transformation automatically if only a forward one was given.

        :param coordinate_transformation: `None`, or a sequence of forward-transformation derivative tensors, or an already-paired `(forward, reverse)` tuple
        :type coordinate_transformation: None | list[np.ndarray] | tuple
        :return: `None` if no transformation was given, otherwise the `(forward, reverse)` derivative-tensor pair
        :rtype: tuple | None
        """
        ...

    def get_embedding(self, masses=None, modes=True, mass_weight=True, dimensionless=False):
        """
        **LLM Docstring**

        Build a `ModeEmbedding` combining this Hamiltonian's coordinate embedding with (optionally) its normal modes.

        :param masses: masses to use in the embedding; if `None`, uses the embedding's default masses
        :type masses: np.ndarray | None
        :param modes: if `True`, use `self.modes`; if `False`, use no modes (`None`); otherwise use the value passed directly as the modes object
        :type modes: bool | MixtureModes
        :param mass_weight: whether coordinates should be mass-weighted
        :type mass_weight: bool
        :param dimensionless: whether normal-mode coordinates should be made dimensionless
        :type dimensionless: bool
        :return: the constructed mode embedding
        :rtype: ModeEmbedding
        """
        ...

    def potential_expansion(self, order=None, *, dimensionless=False, embedding=None, exclude_gradient=False, **embedding_opts):
        """
        **LLM Docstring**

        Build (and optionally evaluate to a given order) a `ScalarExpansion` for the potential energy surface.

        :param order: if given, immediately evaluate the expansion terms up to this order instead of returning the `ScalarExpansion` object itself
        :type order: int | None
        :param dimensionless: whether the embedding used for the expansion should use dimensionless normal-mode coordinates
        :type dimensionless: bool
        :param embedding: an existing embedding to use; if `None`, one is built via `get_embedding`
        :type embedding: ModeEmbedding | None
        :param exclude_gradient: if `True`, zero out the first-order (gradient) term before building the expansion
        :type exclude_gradient: bool
        :param embedding_opts: extra keyword arguments forwarded to `get_embedding` when `embedding` is not given
        :type embedding_opts: dict
        :return: the `ScalarExpansion` object, or its evaluated terms if `order` was given
        :rtype: ScalarExpansion | list[np.ndarray]
        """
        ...

    def dipole_expansion(self, order=None, *, expansion=None, dimensionless=False, embedding=None, exclude_gradient=False, **embedding_opts):
        """
        **LLM Docstring**

        Build (and optionally evaluate to a given order) a `DipoleExpansion` for the dipole surface.

        :param order: if given, immediately evaluate the expansion terms up to this order instead of returning the `DipoleExpansion` object itself
        :type order: int | None
        :param expansion: raw dipole derivative tensors to use in place of `self.dipole.derivatives`
        :type expansion: list[np.ndarray] | None
        :param dimensionless: whether the embedding used for the expansion should use dimensionless normal-mode coordinates
        :type dimensionless: bool
        :param embedding: an existing embedding to use; if `None`, one is built via `get_embedding`
        :type embedding: ModeEmbedding | None
        :param exclude_gradient: accepted for interface consistency with `potential_expansion` but not used in this method's body
        :type exclude_gradient: bool
        :param embedding_opts: extra keyword arguments forwarded to `get_embedding` when `embedding` is not given
        :type embedding_opts: dict
        :return: the `DipoleExpansion` object, or its evaluated terms (with `shared=1`) if `order` was given
        :rtype: DipoleExpansion | list[np.ndarray]
        """
        ...

    def get_potential_optmizing_transformation(self, order, *, dimensionless=False, embedding=None, exclude_gradient=True, **embedding_opts):
        """
        **LLM Docstring**

        Find the coordinate transformation that optimizes (diagonalizes/simplifies) the potential expansion to the given order.

        :param order: the expansion order to optimize
        :type order: int
        :param dimensionless: whether the embedding should use dimensionless normal-mode coordinates
        :type dimensionless: bool
        :param embedding: an existing embedding to use; if `None`, one is built via `get_embedding`
        :type embedding: ModeEmbedding | None
        :param exclude_gradient: whether to zero out the gradient term before optimizing
        :type exclude_gradient: bool
        :param embedding_opts: extra keyword arguments forwarded to `get_embedding`
        :type embedding_opts: dict
        :return: the optimizing transformation, as returned by `McUtils.Numputils.optimizing_transformation`
        :rtype: object
        """
        ...

    def _get_ke_expansion(self, cls, order=None, *, coords=None, masses=None, dimensionless=False, embedding=None, **embedding_opts):
        """
        **LLM Docstring**

        Shared helper for building (and optionally evaluating) a kinetic-energy-related expansion object (`GMatrixExpansion`, `PseudopotentialExpansion`, `CoriolisRotationExpansion`, or `ReciprocalInertiaExpansion`) over a mode embedding.

        :param cls: the expansion class to instantiate with the embedding
        :type cls: type
        :param order: if given, immediately evaluate the expansion terms up to this order
        :type order: int | None
        :param coords: coordinates to pass through to `get_terms` if `order` is given
        :type coords: np.ndarray | None
        :param masses: masses to use when building the embedding
        :type masses: np.ndarray | None
        :param dimensionless: whether the embedding should use dimensionless normal-mode coordinates
        :type dimensionless: bool
        :param embedding: an existing embedding to use; if `None`, one is built via `get_embedding`
        :type embedding: ModeEmbedding | None
        :param embedding_opts: extra keyword arguments forwarded to `get_embedding` when `embedding` is not given
        :type embedding_opts: dict
        :return: the constructed expansion object, or its evaluated terms if `order` was given
        :rtype: object
        """
        ...

    def gmatrix_expansion(self, order=None, *, coords=None, dimensionless=False, embedding=None, **embedding_opts) -> 'GMatrixExpansion|np.ndarray':
        """
        **LLM Docstring**

        Build (and optionally evaluate) the G-matrix (inverse kinetic-energy tensor) expansion, via `_get_ke_expansion` with `GMatrixExpansion`.

        :param order: if given, immediately evaluate the expansion terms up to this order
        :type order: int | None
        :param coords: coordinates forwarded to `get_terms` if `order` is given
        :type coords: np.ndarray | None
        :param dimensionless: whether the embedding should use dimensionless normal-mode coordinates
        :type dimensionless: bool
        :param embedding: an existing embedding to use; if `None`, one is built via `get_embedding`
        :type embedding: ModeEmbedding | None
        :param embedding_opts: extra keyword arguments forwarded to `get_embedding`
        :type embedding_opts: dict
        :return: the `GMatrixExpansion` object, or its evaluated terms if `order` was given
        :rtype: GMatrixExpansion | np.ndarray
        """
        ...

    def get_kinetic_optmizing_transformation(self, order, *, dimensionless=False, embedding=None, **embedding_opts):
        """
        **LLM Docstring**

        Intended to iteratively construct the coordinate transformation that optimizes/decouples the kinetic-energy (G-matrix) expansion order by order. As currently written this method is unfinished/in-progress debugging code: it contains bare `print` statements, large blocks of commented-out alternatives, and unconditionally raises `Exception("???")` inside the loop body on the first iteration, so it never completes and never returns `Q, R` under normal execution.

        :param order: the target expansion order to optimize up to
        :type order: int
        :param dimensionless: whether the embedding should use dimensionless normal-mode coordinates
        :type dimensionless: bool
        :param embedding: an existing embedding to use; if `None`, one is built via `get_embedding`
        :type embedding: ModeEmbedding | None
        :param embedding_opts: extra keyword arguments forwarded to `get_embedding`
        :type embedding_opts: dict
        :return: never returns normally; always raises inside the loop
        :rtype: None
        :raises Exception: unconditionally, with message `"???"`, once `order >= 1`
        """
        ...

    def pseudopotential_expansion(self, order=None, *, dimensionless=False, embedding=None, **embedding_opts):
        """
        **LLM Docstring**

        Build (and optionally evaluate) the pseudopotential expansion, via `_get_ke_expansion` with `PseudopotentialExpansion`.

        :param order: if given, immediately evaluate the expansion terms up to this order
        :type order: int | None
        :param dimensionless: whether the embedding should use dimensionless normal-mode coordinates
        :type dimensionless: bool
        :param embedding: an existing embedding to use; if `None`, one is built via `get_embedding`
        :type embedding: ModeEmbedding | None
        :param embedding_opts: extra keyword arguments forwarded to `get_embedding`
        :type embedding_opts: dict
        :return: the `PseudopotentialExpansion` object, or its evaluated terms if `order` was given
        :rtype: PseudopotentialExpansion | list[np.ndarray]
        """
        ...

    def coriolis_expansion(self, order=None, *, dimensionless=False, embedding=None, **embedding_opts):
        """
        **LLM Docstring**

        Build (and optionally evaluate) the Coriolis-coupling expansion, via `_get_ke_expansion` with `CoriolisRotationExpansion`.

        :param order: if given, immediately evaluate the expansion terms up to this order
        :type order: int | None
        :param dimensionless: whether the embedding should use dimensionless normal-mode coordinates
        :type dimensionless: bool
        :param embedding: an existing embedding to use; if `None`, one is built via `get_embedding`
        :type embedding: ModeEmbedding | None
        :param embedding_opts: extra keyword arguments forwarded to `get_embedding`
        :type embedding_opts: dict
        :return: the `CoriolisRotationExpansion` object, or its evaluated terms if `order` was given
        :rtype: CoriolisRotationExpansion | list[np.ndarray]
        """
        ...

    def watson_expansion(self, order=None, *, dimensionless=False, embedding=None, **embedding_opts):
        """
        **LLM Docstring**

        Build (and optionally evaluate) the Watson (reciprocal-inertia) expansion, via `_get_ke_expansion` with `ReciprocalInertiaExpansion`.

        :param order: if given, immediately evaluate the expansion terms up to this order
        :type order: int | None
        :param dimensionless: whether the embedding should use dimensionless normal-mode coordinates
        :type dimensionless: bool
        :param embedding: an existing embedding to use; if `None`, one is built via `get_embedding`
        :type embedding: ModeEmbedding | None
        :param embedding_opts: extra keyword arguments forwarded to `get_embedding`
        :type embedding_opts: dict
        :return: the `ReciprocalInertiaExpansion` object, or its evaluated terms if `order` was given
        :rtype: ReciprocalInertiaExpansion | list[np.ndarray]
        """
        ...

class ScalarOperatorManager:

    def __init__(self, manager_or_derivs):
        """
        **LLM Docstring**

        Wrap a scalar-operator source (either a manager object exposing a `derivatives` attribute, or raw derivative tensors) for uniform access.

        :param manager_or_derivs: the manager object or raw derivative tensors to wrap
        :type manager_or_derivs: object | list[np.ndarray]
        :return: None
        :rtype: None
        """
        ...

    @property
    def derivatives(self):
        """
        **LLM Docstring**

        The wrapped operator's derivative tensors: delegates to `self.manager.derivatives` if `self.manager` exposes that attribute, otherwise treats `self.manager` itself as the derivative tensors.

        :return: the derivative tensors
        :rtype: list[np.ndarray]
        """
        ...

class PotentialManager(ScalarOperatorManager):

    def __init__(self, manager_or_derivs, modes: NormalModesManager):
        """
        **LLM Docstring**

        Wrap a potential source together with the normal modes needed to canonicalize its derivatives, converting the modes to a dimensionless basis if present.

        :param manager_or_derivs: the potential manager or raw derivative tensors to wrap
        :type manager_or_derivs: object | list[np.ndarray]
        :param modes: the normal-modes manager whose `.modes` (converted to a new, dimensionless basis) will be used for canonicalization
        :type modes: NormalModesManager
        :return: None
        :rtype: None
        """
        ...

    @property
    def derivatives(self):
        """
        **LLM Docstring**

        The potential's canonicalized derivative tensors: fetches the raw derivatives (from `self.manager.derivatives` if present, else `self.manager` itself) and passes them through `self.canonicalize_derivs`. Note that `PotentialManager` inherits from `ScalarOperatorManager`, which does not itself define `canonicalize_derivs` (that method lives on the unrelated `ScalarExpansion` class), so calling this property will raise an `AttributeError` unless a `canonicalize_derivs` method is supplied elsewhere (e.g. via a mixin) -- this looks like a bug/incomplete wiring rather than intended behavior.

        :return: the canonicalized derivative tensors
        :rtype: list[np.ndarray]
        :raises AttributeError: as currently written, because `self.canonicalize_derivs` is not defined on this class or its base
        """
        ...

class ScalarExpansion:
    """
    A helper class that can transform scalar derivatives from Cartesians (or Cartesian modes) to internals
    """

    def __init__(self, embedding: ModeEmbedding, derivs):
        """
        **LLM Docstring**

        Store the mode embedding and raw derivative tensors used to transform a scalar function's Taylor expansion between coordinate systems.

        :param embedding: the mode embedding describing the coordinate transformation context
        :type embedding: ModeEmbedding
        :param derivs: the raw derivative tensors of the scalar function
        :type derivs: list[np.ndarray]
        :return: None
        :rtype: None
        """
        ...

    def canonicalize_derivs(self, derivs, shared=0):
        """
        **LLM Docstring**

        Re-express derivative tensors that are given with respect to mass-weighted Cartesians into the embedding's native coordinate axis, by contracting any axis whose size matches the mass-weighted-inverse transformation's dimension with `self.embedding.mw_inverse()`. Does nothing (returns `derivs` unchanged) if the embedding has no normal modes.

        :param derivs: the derivative tensors to canonicalize
        :type derivs: list[np.ndarray]
        :param shared: number of leading axes to skip when looking for axes to transform (e.g. for expansions that share leading axes, such as dipole components)
        :type shared: int
        :return: the canonicalized derivative tensors
        :rtype: list[np.ndarray]
        """
        ...

    @classmethod
    def check_mixed_expansion(self, derivs, shared=0):
        """
        **LLM Docstring**

        Determine whether a set of derivative tensors is "mixed" -- i.e. has some tensors whose shared axis doesn't match the size of its remaining (internal-coordinate) axes, which happens when higher derivatives are supplied in a different basis than the lower ones.

        :param derivs: the derivative tensors to inspect
        :type derivs: list[np.ndarray]
        :param shared: index of the axis to compare against the trailing axis size
        :type shared: int
        :return: a `(True, mixed_axes)` pair if the expansion is mixed, where `mixed_axes` is the count of that trailing-axis size among the last tensor's shape; otherwise `(False, None)`
        :rtype: tuple[bool, int | None]
        """
        ...

    def get_scalar_transformation(self, order, transformation=None):
        """
        **LLM Docstring**

        Fetch the internals-by-Cartesians and Cartesians-by-internals Jacobian expansions from the embedding (used to transform scalar derivatives between coordinate systems), optionally re-expanding them through an extra coordinate transformation first.

        :param order: the derivative order to fetch the Jacobians to
        :type order: int
        :param transformation: an additional forward/reverse coordinate transformation to apply, passed through `MolecularHamiltonian.prep_transformation`
        :type transformation: tuple | None
        :return: the `(base_tf, rev_tf)` pair of Jacobian-derivative lists
        :rtype: tuple[list[np.ndarray], list[np.ndarray]]
        """
        ...

    def get_terms(self, order=None, *, derivs=None, transformation=None, mixed_transformation=None, mixed_derivative_handling_mode='high', modes=None, shared=0):
        """
        **LLM Docstring**

        Evaluate the scalar function's Taylor-series terms in the embedding's target coordinate system (internal coordinates or normal modes), handling canonicalization, mixed-basis derivative sets, and an optional extra coordinate transformation.

        :param order: the highest derivative order to evaluate; defaults to `len(derivs)` (or `self.derivs` if `derivs` is not given)
        :type order: int | None
        :param derivs: raw derivative tensors to use instead of `self.derivs`
        :type derivs: list[np.ndarray] | None
        :param transformation: if given, first computes the terms without it, then re-expands the result through this extra forward/reverse coordinate transformation
        :type transformation: tuple | None
        :param mixed_transformation: explicit transformation to use when handling a mixed-basis derivative set; if not given and needed, one is derived from the modes
        :type mixed_transformation: list[np.ndarray] | None
        :param mixed_derivative_handling_mode: symmetrization mode passed to `nput.symmetrize_array` when reconciling a mixed-basis expansion
        :type mixed_derivative_handling_mode: str
        :param modes: normal modes to use for building `mixed_transformation` when needed; defaults to `self.embedding.modes`
        :type modes: MixtureModes | None
        :param shared: number of leading axes treated as "shared" (not part of the coordinate re-expansion), e.g. for vector-valued expansions
        :type shared: int
        :return: the evaluated Taylor-expansion terms, one tensor per order
        :rtype: list[np.ndarray]
        """
        ...

class DipoleExpansion(ScalarExpansion):

    def __init__(self, embedding: ModeEmbedding, derivs):
        """
        **LLM Docstring**

        Set up a `ScalarExpansion` specialized for the (vector-valued) dipole surface by moving the trailing length-3 Cartesian-component axis of each derivative tensor to the front before storing them.

        :param embedding: the mode embedding describing the coordinate transformation context
        :type embedding: ModeEmbedding
        :param derivs: the raw dipole derivative tensors, each with a trailing axis of length 3 for the dipole's x/y/z components
        :type derivs: list[np.ndarray]
        :return: None
        :rtype: None
        """
        ...

    def get_terms(self, order=None, *, derivs=None, transformation=None, mixed_transformation=None, mixed_derivative_handling_mode='high', modes=None, shared=0):
        """
        **LLM Docstring**

        Evaluate the dipole's Taylor-series terms, temporarily separating the constant (reference) dipole value from the higher-order derivatives so only the derivatives are passed through the generic `ScalarExpansion.get_terms` machinery, then restoring the x/y/z component axis to its original (trailing) position in the result.

        :param order: the highest derivative order to evaluate
        :type order: int | None
        :param derivs: raw dipole derivative tensors to use instead of `self.derivs`; if given, are first reshaped to move the component axis to the front, matching the constructor's convention
        :type derivs: list[np.ndarray] | None
        :param transformation: an additional forward/reverse coordinate transformation to apply
        :type transformation: tuple | None
        :param mixed_transformation: explicit transformation for a mixed-basis derivative set
        :type mixed_transformation: list[np.ndarray] | None
        :param mixed_derivative_handling_mode: symmetrization mode used when reconciling a mixed-basis expansion
        :type mixed_derivative_handling_mode: str
        :param modes: normal modes to use for building `mixed_transformation` when needed
        :type modes: MixtureModes | None
        :param shared: number of additional leading shared axes beyond the dipole-component axis (which is always treated as shared)
        :type shared: int
        :return: the dipole expansion terms, with the x/y/z component axis restored to the trailing position
        :rtype: list[np.ndarray]
        """
        ...

class GMatrixExpansion:

    def __init__(self, embedding: ModeEmbedding):
        """
        **LLM Docstring**

        Store the mode embedding used to build the G-matrix (inverse effective mass / kinetic-energy metric tensor) expansion.

        :param embedding: the mode embedding describing the coordinate system the G-matrix is expanded in
        :type embedding: ModeEmbedding
        :return: None
        :rtype: None
        """
        ...

    def get_terms(self, order, coords=None, transformation=None):
        """
        **LLM Docstring**

        Compute the G-matrix Taylor-expansion terms in the embedding's coordinates, by contracting the internals-by-Cartesians Jacobian with itself and, for orders beyond zero, re-expanding through the Cartesians-by-internals Jacobian. If `transformation` is given, delegates to `reexpress` instead.

        :param order: the highest derivative order to compute
        :type order: int
        :param coords: coordinates forwarded to the embedding's Jacobian-fetching calls
        :type coords: np.ndarray | None
        :param transformation: if given, an additional forward/reverse coordinate transformation to reexpress the G-matrix through instead of computing it directly
        :type transformation: tuple | None
        :return: the G-matrix expansion terms, one tensor per order
        :rtype: list[np.ndarray]
        """
        ...

    @classmethod
    def _dRGQ_partition_contrib(cls, partition, R, G):
        """
        **LLM Docstring**

        Compute one term of the derivative of the reverse-transformed G-matrix with respect to the new coordinates, for a given `(r1, r2, s)` integer partition -- i.e. the contribution from taking `r1` and `r2` derivatives of the reverse transformation and `s` derivatives of the original G-matrix, including the combinatorial symmetrization over equivalent index permutations.

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

        Sum the contributions of every valid `(r1, r2, g)` integer partition of order `o+1` (restricted to `r1 >= r2` by symmetry) to get the full order-`o` derivative of the reverse-transformed G-matrix.

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
    def reexpress_G(self, G_expansion, forward_derivs, reverse_derivs, order, *, GR_expansion=None, return_GR=False):
        """
        Apply a coordinate transformation to the G-matrix

        :param forward_derivs:
        :param reverse_derivs:
        :param order:
        :return:
        """
        ...

    def reexpress(self, order, forward_derivs, reverse_derivs):
        """
        **LLM Docstring**

        Re-express this object's own G-matrix expansion (obtained via `get_terms`) through a new forward/reverse coordinate transformation, by delegating to `reexpress_G`.

        :param order: the derivative order(s) to reexpress
        :type order: int | list[int]
        :param forward_derivs: the forward-transformation derivative tensors (new coordinates as a function of old)
        :type forward_derivs: list[np.ndarray]
        :param reverse_derivs: the reverse-transformation derivative tensors (old coordinates as a function of new)
        :type reverse_derivs: list[np.ndarray]
        :return: the re-expressed G-matrix expansion terms
        :rtype: list[np.ndarray]
        """
        ...

class PseudopotentialExpansion:

    def __init__(self, g_expansion: 'ModeEmbedding|GMatrixExpansion'):
        """
        **LLM Docstring**

        Wrap (or build) a `GMatrixExpansion` used as the basis for computing the pseudopotential (Watson) correction term.

        :param g_expansion: an existing `GMatrixExpansion`, or a `ModeEmbedding` from which one will be constructed
        :type g_expansion: ModeEmbedding | GMatrixExpansion
        :return: None
        :rtype: None
        """
        ...

    @classmethod
    def lambda_expansion(cls, G_expansion, G_inv, I0_expansion, I0_inv, order):
        """
        **LLM Docstring**

        Compute the expansion of `ln(det(I0)) - ln(det(G))`-type combination (the trace terms appearing in the pseudopotential derivation) by taking traces of the tensordot-derivative products of the G-matrix and inertia-tensor expansions with their respective inverses.

        :param G_expansion: the G-matrix expansion terms
        :type G_expansion: list[np.ndarray]
        :param G_inv: the expansion terms of the G-matrix inverse
        :type G_inv: list[np.ndarray]
        :param I0_expansion: the inertia-tensor expansion terms
        :type I0_expansion: list[np.ndarray]
        :param I0_inv: the expansion terms of the inertia-tensor inverse
        :type I0_inv: list[np.ndarray]
        :param order: the derivative order to compute
        :type order: int
        :return: the term-by-term difference between the inertia-tensor and G-matrix trace expansions
        :rtype: list[np.ndarray]
        """
        ...

    @classmethod
    def get_U(cls, G_expansion, I0_expansion):
        """
        **LLM Docstring**

        Compute the pseudopotential correction term `U` from the G-matrix and inertia-tensor (`I0`) expansions, following the standard Watson pseudopotential derivation (inverting both expansions, forming the `lambda` trace expansion, and contracting it against the G-matrix expansion).

        :param G_expansion: the G-matrix expansion terms, up to order `n = len(G_expansion) - 1`
        :type G_expansion: list[np.ndarray]
        :param I0_expansion: the inertia-tensor expansion terms
        :type I0_expansion: list[np.ndarray]
        :return: the pseudopotential expansion terms
        :rtype: list[np.ndarray]
        """
        ...

    def get_terms(self, order, transformation=None):
        """
        **LLM Docstring**

        Evaluate the pseudopotential expansion terms up to `order`, either directly from the G-matrix and inertia-tensor expansions (computed two orders past `order` to account for extra derivatives needed internally), or, if `transformation` is given, via `reexpress`.

        :param order: the target expansion order
        :type order: int
        :param transformation: an additional forward/reverse coordinate transformation to reexpress the result through
        :type transformation: tuple | None
        :return: the pseudopotential expansion terms
        :rtype: list[np.ndarray]
        """
        ...

    def reexpress(self, order, forward_derivs, reverse_derivs):
        """
        Apply a coordinate transformation to the G-matrix

        :param forward_derivs:
        :param reverse_derivs:
        :param order:
        :return:
        """
        ...

class ZetaExpansion:

    def __init__(self, embedding: ModeEmbedding):
        """
        **LLM Docstring**

        Store the mode embedding used to compute the Coriolis zeta-constant expansion.

        :param embedding: the mode embedding describing the coordinate system
        :type embedding: ModeEmbedding
        :return: None
        :rtype: None
        """
        ...

    def get_terms(self, order, transformation=None):
        """
        **LLM Docstring**

        Compute the Coriolis zeta-constant expansion terms: builds the internals-by-Cartesians Jacobian (optionally reexpressed through an extra transformation), projects it into the inertial (principal-axis) frame, forms its self tensor product, and antisymmetrizes over pairs of Cartesian axes to extract the `x`, `y`, `z` zeta-constant components.

        :param order: the derivative order to compute
        :type order: int
        :param transformation: an additional forward/reverse coordinate transformation to apply to the internals-by-Cartesians Jacobian first
        :type transformation: tuple | None
        :return: the zeta expansion terms, three (x/y/z) entries per order
        :rtype: list[np.ndarray]
        """
        ...

class ReciprocalInertiaExpansion:
    """
    Pulled from Watson
    """

    def __init__(self, embedding: ModeEmbedding):
        """
        **LLM Docstring**

        Store the mode embedding used to compute the reciprocal (inverse) moment-of-inertia expansion appearing in the Watson term.

        :param embedding: the mode embedding describing the coordinate system
        :type embedding: ModeEmbedding
        :return: None
        :rtype: None
        """
        ...

    def get_terms(self, order, transformation=None):
        """
        **LLM Docstring**

        Compute the expansion of the inverse inertia tensor `u = I0^-1` order by order, using the recursive relation `u^(n) = a @ u^(n-1) @ u` built from the first-order inertia-tensor derivative `a` (optionally re-expressed through an extra coordinate transformation).

        :param order: the highest derivative order to compute
        :type order: int
        :param transformation: an additional forward/reverse coordinate transformation to apply to the first-order inertia-tensor derivative
        :type transformation: tuple | None
        :return: the reciprocal-inertia expansion terms, `order + 1` entries starting from `u` itself
        :rtype: list[np.ndarray]
        """
        ...

class CoriolisRotationExpansion:
    """
    Provides an expansion of the Coriolis rotation operator
    """

    def __init__(self, embedding: ModeEmbedding):
        """
        **LLM Docstring**

        Store the mode embedding and build the sub-expansions (zeta constants and reciprocal inertia) needed to compute the Coriolis rotational-coupling operator's expansion.

        :param embedding: the mode embedding describing the coordinate system
        :type embedding: ModeEmbedding
        :return: None
        :rtype: None
        """
        ...

    def get_terms(self, order, transformation=None):
        """
        **LLM Docstring**

        Compute the Coriolis rotation operator's expansion terms by combining the zeta-constant expansion (contracted with itself) with the reciprocal-inertia expansion.

        :param order: the derivative order to compute
        :type order: int
        :param transformation: an additional forward/reverse coordinate transformation forwarded to the zeta and reciprocal-inertia sub-expansions
        :type transformation: tuple | None
        :return: the Coriolis rotation expansion terms
        :rtype: list[np.ndarray]
        """
        ...