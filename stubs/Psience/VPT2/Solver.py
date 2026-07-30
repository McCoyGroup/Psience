import numpy as np, itertools, time, types
from collections import namedtuple
from McUtils.Numputils import SparseArray
import McUtils.Numputils as nput
from McUtils.Scaffolding import Logger, NullLogger, NullCheckpointer
from McUtils.Parallelizers import Parallelizer, SerialNonParallelizer
from McUtils.Data import UnitsData
from McUtils.Combinatorics import LatticePathGenerator
from ..BasisReps import Representation, BasisStateSpace, BasisMultiStateSpace, SelectionRuleStateSpace, BasisStateSpaceFilter, HarmonicOscillatorProductBasis
from .DegeneracySpecs import DegenerateMultiStateSpace, DegeneracySpec
from .Common import *
from .Corrections import *
__reload_hook__ = ['..BasisReps', '.DegeneracySpecs', '.Corrections', '.Common']
__all__ = ['PerturbationTheorySolver']

class PerturbationTheorySolver:
    """
    A solver that applies perturbation theory
    given a series of corrections and population of states.
    Supports degenerate and non-degenerate PT.
    """

    def __init__(self, perturbations, states, coupled_states=None, order=2, total_space=None, flat_total_space=None, state_space_iterations=None, state_space_terms=None, state_space_filters=None, extended_state_space_filter_generator=None, extended_state_space_postprocessor=None, target_property_rules=None, allow_sakurai_degs=False, allow_post_PT_calc=True, modify_degenerate_perturbations=False, gaussian_resonance_handling=False, ignore_odd_order_energies=False, intermediate_normalization=False, reorthogonalize_degenerate_states=None, check_overlap=True, zero_element_warning=True, degenerate_states=None, degeneracy_handlers=None, handle_strong_couplings=False, local_coupling_hamiltonian=None, local_coupling_order=None, low_frequency_mode_cutoff=0.00115, zero_order_energy_corrections=None, nondeg_hamiltonian_precision=3, memory_constrained=False, keep_hamiltonians=None, logger=None, parallelizer=None, checkpointer=None, results=None, checkpoint_keys=None, use_cached_representations=False, use_cached_basis=False):
        """

        :param perturbations:
        :type perturbations: Iterable[Representation]
        :param states:
        :type states: BasisStateSpace
        :param coupled_states:
        :type coupled_states: BasisMultiStateSpace
        :param order:
        :type order:
        :param degenerate_states:
        :type degenerate_states:
        :param degeneracy_mode:
        :type degeneracy_mode:
        :param logger:
        :type logger:
        :param parallelizer:
        :type parallelizer:
        :param checkpointer:
        :type checkpointer:
        :param results:
        :type results:
        """
        ...

    @property
    def coupled_states(self):
        """

        :return:
        :rtype:
        """
        ...

    @property
    def total_space_dim(self):
        """

        :return:
        :rtype:
        """
        ...

    @property
    def flat_total_space(self):
        """

        :return:
        :rtype:
        """
        ...

    @property
    def total_state_space(self):
        """

        :return:
        :rtype:
        """
        ...

    class PastIndexableTuple(list):

        def __getitem__(self, item):
            """
            **LLM Docstring**

            Index into the list, transparently returning `0` for any out-of-range integer index instead of raising, so higher-order perturbation terms that weren't explicitly computed are treated as identically zero.

            :param item: the index (or slice) to look up
            :type item: int | slice
            :return: the stored value, or `0` if `item` is an out-of-range integer index
            :rtype: object
            """
            ...

        def __setitem__(self, key, value):
            """
            **LLM Docstring**

            Set the value at a given index, first extending the list with zero placeholders if the index is beyond the current length. Note the extension amount (`len(self) - key`) is non-positive whenever `key >= len(self)`, so `extend` receives an empty (or negative-length, effectively empty) list and the list is never actually grown to include `key` -- this looks like an off-by-sign bug (it should extend by `key - len(self) + 1`), and the subsequent `super().__setitem__` call would raise an `IndexError` for any `key` at or beyond the original length.

            :param key: the index to set
            :type key: int
            :param value: the value to store
            :type value: object
            :return: the result of the (likely failing, for out-of-range `key`) base `list.__setitem__` call
            :rtype: object
            """
            ...

    @property
    def representations(self):
        """
        :return:
        :rtype: Iterable[SparseArray]
        """
        ...

    @representations.setter
    def representations(self, reps):
        """
        **LLM Docstring**

        The (cached) list of Hamiltonian perturbation representation matrices, computed lazily via `get_VPT_representations` and wrapped in a `PastIndexableTuple` so higher, uncomputed orders read back as zero.

        :return: the perturbation representation matrices
        :rtype: PerturbationTheorySolver.PastIndexableTuple
        """
        ...

    @property
    def degenerate_spaces(self):
        """
        :return:
        :rtype:
        """
        ...

    @classmethod
    def merge_deg_spaces(cls, new_states):
        """
        **LLM Docstring**

        Combine several independently-identified sets of degenerate state groups into one consistent set, by flattening every group (from either raw `BasisStateSpace`s or `DegenerateMultiStateSpace`s) down to its excitation vectors and merging any that share states via `DegeneracySpec.merge_state_blocks`.

        :param new_states: the separate sets of degenerate groups to merge, each either a `DegenerateMultiStateSpace` or an iterable of raw state blocks
        :type new_states: list
        :return: the merged degenerate groups
        :rtype: list[np.ndarray]
        """
        ...

    @property
    def zero_order_energies(self):
        """
        :return:
        :rtype:
        """
        ...

    def apply_VPT(self):
        """
        Applies perturbation theory to the held basis of states using the
        built representations and degenerate state spaces

        :return:
        :rtype: PerturbationTheoryCorrections
        """
        ...

    def get_VPT_representations(self):
        """
        Gets the sparse representations of the passed perturbation inside the basis of coupled states.

        :return:
        :rtype: Iterable[SparseArray]
        """
        ...

    def extend_VPT_representations(self, new_flat_space, new_states):
        """
        **LLM Docstring**

        Extend the cached Hamiltonian perturbation representation matrices to cover a newly added block of coupled states, by computing just the new matrix elements (between the new states and the full flat space) and splicing them into the existing sparse representation data rather than recomputing everything from scratch.

        :param new_flat_space: the newly added states, flattened, that the zeroth-order representation needs new diagonal elements for
        :type new_flat_space: BasisStateSpace
        :param new_states: the newly added coupled-state blocks, one per higher-order perturbation, that need new off-diagonal elements computed
        :type new_states: list[BasisStateSpace]
        :return: the extended list of representation matrices
        :rtype: list[SparseArray]
        """
        ...

    def _take_subham(self, rep, inds):
        """
        Builds a subsampled version of a representation Hamiltonian
        to allow equations to be efficiently solved in subspaces.

        :param rep: representation matrix from which to pull the subspace
        :type rep: SparseArray
        :param inds: indices for the subspace
        :type inds: np.ndarray
        :return:
        :rtype:
        """
        ...

    def _build_projector(self, inds):
        """
        Builds a projector where only inds will
        be included

        :param inds: indices for the subspace
        :type inds: np.ndarray
        :return:
        :rtype:
        """
        ...

    def _get_Pi0(self, degenerate_subspace, non_zero_cutoff=None, E0=None):
        """
        **LLM Docstring**

        Build the perturbation operator `Pi0`
        -- a diagonal matrix of `1/(E0_ref - E0_i)` for every state `i` outside the given degenerate subspace, with zero entries for states inside it --
        raising an error if any state outside the subspace is unexpectedly (near-)degenerate with the reference energy.

        :param degenerate_subspace: the indices of the states in the (single) degenerate group `Pi0` should exclude/zero out
        :type degenerate_subspace: np.ndarray
        :param non_zero_cutoff: the energy-gap magnitude below which a state is considered problematically degenerate; defaults to `Settings.non_zero_cutoff`
        :type non_zero_cutoff: float | None
        :param E0: the reference zeroth-order energy to measure gaps relative to; defaults to the average energy over `degenerate_subspace`
        :type E0: float | np.ndarray | None
        :return: the diagonal `Pi0` operator
        :rtype: SparseArray
        :raises ValueError: if more than one state outside `degenerate_subspace` is unexpectedly near-degenerate with `E0`
        """
        ...

    def load_state_spaces(self):
        """

        :return:
        :rtype:
        """
        ...

    def _generate_total_space(self):
        """
        **LLM Docstring**

        Build (if not already present) the total coupled state space and its flattened, deduplicated form spanning the target states plus every coupled-state block, computing and logging the total dimension either from scratch or by flattening an already-known total space.

        :return: None
        :rtype: None
        """
        ...

    def extend_state_spaces(self, new_targets, degenerate_states=None):
        """
        **LLM Docstring**

        Extend the solver's state spaces to additionally include a set of new target states: (re)computes the coupled-state blocks needed for the new targets (applying any configured extended-state-space filter/postprocessor), merges them into the existing per-order coupled-state spaces and flattened total space, and updates `self.states`/`self._total_space` accordingly.

        :param new_targets: the new target states to add to the solve
        :type new_targets: BasisStateSpace
        :param degenerate_states: the degenerate-state groupings driving this extension, forwarded to `extended_state_space_postprocessor` if one is configured
        :type degenerate_states: object | None
        :return: `(flat_space, new_spaces)` -- the newly added (deduplicated) flat states and the per-order newly added coupled-state blocks, or `None` if nothing new was actually found
        :rtype: tuple | None
        """
        ...

    def load_coupled_spaces(self, degenerate_spaces=None, spaces=None, wavefunction_terms=None, property_filter=None, filter_spaces=None):
        """
        Determines which states need to be coupled at which levels of correction
        to handle the PT
        :return:
        :rtype:
        """
        ...

    class StateSpaceWrapper:
        """
        Wraps a state space so that it can define stuff like __add__, __mul__, and __neg__
        """

        def __init__(self, space):
            """
            **LLM Docstring**

            Wrap a state space so it supports the symbolic `+`/`-`/unary-`-` operations used when combining state-space contributions symbolically.

            :param space: the state space to wrap
            :type space: BasisStateSpace | object
            :return: None
            :rtype: None
            """
            ...

        def __neg__(self):
            """
            **LLM Docstring**

            Negation is a no-op for state spaces (there's no notion of a "negative" set of states), so this just returns `self` unchanged.

            :return: `self`
            :rtype: PerturbationTheorySolver.StateSpaceWrapper
            """
            ...

        def simple_union(self, other):
            """
            **LLM Docstring**

            Combine this wrapped state space with another (treating `None`/`0` as "nothing to add"), unwrapping `other` first if it's also a `StateSpaceWrapper`.

            :param other: the other state space (or wrapper) to union with; `None` or `0` are treated as empty
            :type other: object | PerturbationTheorySolver.StateSpaceWrapper | None
            :return: `self` if `other` is empty, otherwise a new wrapper around the unioned space
            :rtype: PerturbationTheorySolver.StateSpaceWrapper
            """
            ...

        def __sub__(self, other):
            """
            **LLM Docstring**

            Since state spaces have no subtraction, `-` is treated the same as `+` -- both simply union the two spaces together.

            :param other: the other state space (or wrapper) to combine with
            :type other: object
            :return: the unioned wrapper
            :rtype: PerturbationTheorySolver.StateSpaceWrapper
            """
            ...

        def __rsub__(self, other):
            """
            **LLM Docstring**

            Reflected subtraction; identical to `__sub__`/`simple_union` since combining spaces is commutative and subtraction is treated as union.

            :param other: the other state space (or wrapper) to combine with
            :type other: object
            :return: the unioned wrapper
            :rtype: PerturbationTheorySolver.StateSpaceWrapper
            """
            ...

        def __add__(self, other):
            """
            **LLM Docstring**

            Combine this wrapped state space with another via union.

            :param other: the other state space (or wrapper) to combine with
            :type other: object
            :return: the unioned wrapper
            :rtype: PerturbationTheorySolver.StateSpaceWrapper
            """
            ...

        def __radd__(self, other):
            """
            **LLM Docstring**

            Reflected addition; identical to `__add__` since union is commutative.

            :param other: the other state space (or wrapper) to combine with
            :type other: object
            :return: the unioned wrapper
            :rtype: PerturbationTheorySolver.StateSpaceWrapper
            """
            ...

    class ProjectionOperatorWrapper:
        """
        Generates a symbolic form of a perturbation operator that
        either projects onto a degenerate space or projects it out
        """

        def __init__(self, space, complement=False):
            """
            **LLM Docstring**

            Set up a symbolic operator that, when applied to a state space, either restricts it to (or excludes) a fixed set of states -- used to represent projecting onto or away from a degenerate subspace.

            :param space: the fixed set of states to project onto/away from
            :type space: BasisStateSpace
            :param complement: whether this operator excludes `space` (`True`) rather than restricting to it (`False`)
            :type complement: bool
            :return: None
            :rtype: None
            """
            ...

        def get_transformed_space(self, other):
            """
            :param other:
            :type other: SelectionRuleStateSpace
            :return:
            :rtype:
            """
            ...

    class ProjectedOperator:
        """
        Generates a symbolic form of an operator where
        an operator can first be applied, then unused terms projected
        out, before returning the state space
        """

        def __init__(self, projector, operator):
            """
            **LLM Docstring**

            Set up a symbolic operator that first applies an underlying operator's state-space transformation, then filters the result through a projector -- used to compose "apply this coupling, then keep/drop degenerate states" symbolically without materializing intermediate state spaces.

            :param projector: the projection operator to apply after `operator`
            :type projector: object
            :param operator: the underlying operator whose transformed space gets projected
            :type operator: object
            :return: None
            :rtype: None
            """
            ...

        def get_transformed_space(self, other):
            """
            :param other:
            :type other: BasisStateSpace
            :return:
            :rtype:
            """
            ...

    def _could_be_a_space(self, test):
        """
        **LLM Docstring**

        Heuristically check whether an object could plausibly represent a state space: an actual `BasisStateSpace`/`BasisMultiStateSpace`, a flat list of integer indices, or a nested list of integer excitation vectors.

        :param test: the object to check
        :type test: object
        :return: whether `test` looks like it could be a state space
        :rtype: bool
        """
        ...

    def _could_be_rules(self, test):
        """
        **LLM Docstring**

        Heuristically check whether an object could plausibly represent a set of selection rules: `None` (no rules), an empty sequence, or a sequence whose first element is itself empty or starts with an integer.

        :param test: the object to check
        :type test: object
        :return: whether `test` looks like it could be a set of selection rules
        :rtype: bool
        """
        ...

    def _could_be_a_prefilter(self, test):
        """
        **LLM Docstring**

        Heuristically check whether an object could plausibly represent a single prefilter specification: `None`, a plain state space (via `_could_be_a_space`), or a `(space, rules)` pair matching `_could_be_a_space`/`_could_be_rules` respectively.

        :param test: the object to check
        :type test: object
        :return: whether `test` looks like a valid prefilter specification
        :rtype: bool
        """
        ...

    def _apply_transformation_with_filters(self, a, b, filter_space: BasisStateSpaceFilter, **opts):
        """
        **LLM Docstring**

        Apply a coupling operator's state-space transformation to `b` (interpreted relative to `a`), then run the result through any pre/post state-space filters attached to `filter_space`.

        :param a: the source/reference state space the transformation conceptually originates from
        :type a: object
        :param b: the operator/state-space object whose `get_transformed_space` should be applied
        :type b: object
        :param filter_space: the filter set (pre- and post-transformation filters) to apply to the result; if `None`, no filtering is applied
        :type filter_space: BasisStateSpaceFilter | None
        :param opts: extra options forwarded to the underlying transformation/filtering calls
        :type opts: dict
        :return: the (possibly filtered) transformed state space
        :rtype: object
        """
        ...

    def _get_new_coupled_space(self, a, b, spaces=None, ret_space=True, filter_space=None):
        """
        A symbolic version of the dot product appropriate for getting
        transformed state spaces under the operation of a on b

        :param a:
        :type a:
        :param b:
        :type b: SelectionRuleStateSpace
        :param spaces:
        :type spaces: the set of operators to which we can assign transformations
        :return:
        :rtype: StateSpaceWrapper
        """
        ...

    def _reduce_new_coupled_space(self, *terms, spaces=None, ret_space=True, filter_space=None):
        """
        Reduces through `_get_new_coupled_space` from right to left
        :param terms:
        :type terms: Iterable[SelectionRuleStateSpace]
        :param spaces:
        :type spaces: dict | None
        :return:
        :rtype:
        """
        ...

    def get_coupled_space(self, input_state_space, degenerate_space, use_second_deg, allow_PT_degs=True, wavefunction_terms=None, spaces=None, property_filter=None, filter_spaces=None):
        """
        Applies the VPT equations semi-symbolically, dispatching based on how many
        degeneracies we need to handle

        :return:
        :rtype:
        """
        ...

    def get_nondeg_coupled_space(self, input_state_space, degenerate_space=None, spaces=None, wavefunction_terms=None, property_filter=None, filter_spaces=None):
        """
        Applies the non-degenerate equations in semi-symbolic form to determine
        which states needs to be calculated.
        This will always be the initial input to a calculation and then
        certain additional states may be calculated on the fly if they are needed to handle
        truly degenerate stuff.
        The only difference there will be to add on

        :param input_state_space:
        :type input_state_space: BasisStateSpace
        :param degenerate_space:
        :type degenerate_space: BasisStateSpace
        :param spaces:
        :type spaces:
        :param wavefunction_terms: which terms to include when calculating corrections
        :type wavefunction_terms: None | Iterable[Iterable[int]]
        :param property_filter: a set of states and selection rules to allow for being targeted in state to calculate
        :type property_filter: None | Iterable[Iterable[int], Iterable[Iterable[int]]]
        :return:
        :rtype:
        """
        ...

    def get_corrections(self, non_zero_cutoff=None, handle_strong_couplings=None, check_overlap=None):
        """
        Applies the perturbation theory equations to obtain
        corrections to the wave functions and energies

        :return:
        :rtype:
        """
        ...

    @property
    def high_frequency_modes(self):
        """
        **LLM Docstring**

        The indices of the vibrational modes whose fundamental transition frequency exceeds `self.low_frequency_mode_cutoff`, used as the default set of modes considered for strong-coupling/degeneracy testing.

        :return: the high-frequency mode indices
        :rtype: list[int]
        """
        ...
    default_strong_coupling_threshold = 0.3

    def identify_strong_couplings(self, corrs):
        """
        **LLM Docstring**

        Find the strongly-coupled state pairs among the solved corrections, delegating to whichever configured `degeneracy_handlers` spec exposes a `wfc_threshold` (or, if none do, a freshly built default `'auto'` spec).

        :param corrs: the perturbation-theory corrections to search for strong couplings in
        :type corrs: PerturbationTheoryCorrections
        :return: `(strong_couplings, threshold)` -- the identified strong-coupling data and the threshold used to find it
        :rtype: tuple
        """
        ...

    def construct_strong_coupling_spaces(self, spec, sc, corrs, states, threshold):
        """
        **LLM Docstring**

        Build the degenerate state groups implied by a set of identified strong couplings, and, if the spec allows extending the state space (`spec.extend_spaces`), extend the solver's states/representations to cover any newly implicated states before returning the resulting (possibly larger) state space and representation data.

        :param spec: the degeneracy spec (typically a `StronglyCoupledDegeneracySpec`) whose group filter/extension settings govern this construction
        :type spec: DegeneracySpec
        :param sc: the raw strong-coupling data (per-state, per-order coupled-state lists) to build groups from
        :type sc: dict
        :param corrs: the perturbation-theory corrections the coupling data came from, used to collapse it into per-state coupled-state spaces
        :type corrs: PerturbationTheoryCorrections
        :param states: the current target state space
        :type states: BasisStateSpace
        :param threshold: the coupling-strength threshold to build the group filter with
        :type threshold: float
        :return: `(degenerate_states, (states, perturbations, flat_total_space, N))` -- the identified degenerate groups, and the (possibly extended) states/representations/flat space/dimension to use going forward
        :rtype: tuple
        """
        ...

    def _get_corrections(self, perturbations, states, order, flat_total_space, N, checkpointer, logger, degenerate_states, non_zero_cutoff=None, handle_strong_couplings=None):
        """
        Just exists so we can do a recursion on this...

        :param perturbations:
        :type perturbations:
        :param states:
        :type states:
        :param order:
        :type order:
        :param flat_total_space:
        :type flat_total_space:
        :param N:
        :type N:
        :param checkpointer:
        :type checkpointer:
        :param logger:
        :type logger:
        :param degenerate_states:
        :type degenerate_states:
        :param non_zero_cutoff:
        :type non_zero_cutoff:
        :param handle_strong_couplings:
        :type handle_strong_couplings:
        :return:
        :rtype:
        """
        ...
    PTResults = namedtuple('PTResults', ['corrections', 'degeneracies'])

    @staticmethod
    def _safe_dot(a, b):
        """
        **LLM Docstring**

        Generalized dot product that treats a bare `0` operand (rather than a zero-valued array) as an exact zero result, short-circuiting the multiplication -- used throughout the perturbation-theory recursion, where many terms really are the literal integer `0` rather than zero arrays.

        :param a: the left operand, or `0`
        :type a: np.ndarray | SparseArray | int | float
        :param b: the right operand, or `0`
        :type b: np.ndarray | SparseArray | int | float
        :return: `0` if either operand is a numeric zero, otherwise the dot product of `a` and `b`
        :rtype: np.ndarray | SparseArray | int
        """
        ...

    def apply_VPT_equations(self, state_index, degenerate_space_indices, degenerate_energies, zero_order_state, degenerate_subspace, degenerate_subsubspace, perturbations=None, allow_PT_degs=None, ignore_odd_orders=None, intermediate_normalization=None, non_zero_cutoff=None):
        """
        Applies VPT equations, dispatching based on how many
        degeneracies we need to handle

        :param state_index: the index of the primary state being treated using the PT
        :type state_index: int
        :param degenerate_space_indices: the indices corresponding to degeneracies with the primary state in the zero-order picture
        :type degenerate_space_indices: np.ndarray[int]
        :param degenerate_energies: the first and (possibly) second order correction to the energies
        :type degenerate_energies: Iterable[float | None]
        :param zero_order_states: the vector for the proper zero-order state corresponding ot state_index
        :type zero_order_states: np.ndarray[float]
        :param degenerate_subsubspace: the set of vectors for the zero-order states in the secondary degenerate subspace
        :type degenerate_subsubspace: tuple[np.ndarray[float], np.ndarray[int]]
        :param non_zero_cutoff: cutoff for when a term can be called zero for performance reasons
        :type non_zero_cutoff: float
        :return:
        :rtype:
        """
        ...

    def apply_VPT_nondeg_equations(self, state_index, deg_group, perturbations=None, non_zero_cutoff=None, check_overlap=None, intermediate_normalization=False, ignore_odd_orders=False):
        """
        Does the dirty work of doing the VPT iterative equations.

        :return:
        :rtype:
        """
        ...

    def apply_VPT_2k1_rules(self, existing_corrs, perturbations=None):
        """
        Apply expressions allowing for obtaining higher-order
        corrections to the energies from lower-order corrections to the
        wavefunctions

        :param existing_corrs:
        :type existing_corrs:
        :param perturbations:
        :type perturbations:
        :return:
        :rtype:
        """
        ...

    @staticmethod
    def _fmt_depert_engs(total_state_space, base_energies):
        """
        **LLM Docstring**

        Build a lazily-evaluated formatter for logging a set of deperturbed state energies: subtracting off the ground-state energy to display frequencies (rather than raw energies) whenever the ground state is present in `total_state_space`.

        :param total_state_space: the state space the energies correspond to
        :type total_state_space: BasisStateSpace
        :param base_energies: the (deperturbed) energies to format
        :type base_energies: np.ndarray
        :return: a zero-argument formatter function producing the log lines, for use with the logger's lazy message-formatting
        :rtype: callable
        """
        ...

    def apply_post_PT_variational_calc(self, degenerate_states, corrs):
        """
        Applies degenerate perturbation theory by building a representation
        for the degenerate terms in the Hamiltonian.
        This is then diagonalized, allowing the degenerate states to be expressed
        in the basis of non-degenerate states

        :param H:
        :type H: Iterable[SparseArray]
        :param corrs: the standard PerturbationTheory Corrections object that comes out of the application of non-deg PT
        :type corrs: PerturbationTheoryCorrections
        :param degenerate_states: population of degenerate states
        :type degenerate_states:
        :param logger:
        :type logger: Logger
        :return:
        :rtype:
        """
        ...

    def drop_deg_pert_els(self, perts, deg_groups):
        """

        :param perts:
        :type perts:
        :param deg_groups:
        :type deg_groups:
        :return:
        :rtype:
        """
        ...