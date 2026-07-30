import collections
import itertools
import numpy as np, enum, abc, scipy.sparse as sp
from McUtils.Combinatorics import SymmetricGroupGenerator, PermutationRelationGraph
from McUtils.Scaffolding import Logger
import McUtils.Numputils as nput
from ..BasisReps import BasisStateSpace, BasisMultiStateSpace, SelectionRuleStateSpace, BraKetSpace, HarmonicOscillatorProductBasis
__all__ = ['DegeneracySpec', 'DegenerateMultiStateSpace']
__reload_hook__ = ['..BasisReps']

class DegenerateSpaceInputFormat(enum.Enum):
    """Real access pattern: DegenerateSpaceInputFormat.<MemberName> (this is an enum with 8 members, e.g. DegenerateSpaceInputFormat.Groups == 'groups'). Collapsed into a dict below purely for compactness -- do not index it as a dict in real code:"""
    _MEMBERS = {'Groups': 'groups', 'QuantaSpecRules': 'nT', 'Polyads': 'polyads', 'Pairs': 'pairs', 'StrongCouplings': 'wfc_threshold', 'EnergyCutoff': 'energy_window', 'Callable': 'callable', 'MartinTest': 'martin_threshold'}

class DegeneracySpec(metaclass=abc.ABCMeta):
    """
    Provides a container for specifying degeneracies
    in a way that can be cleanly canonicalized
    """
    application_order = 'pre'
    group_filter = None

    def __init__(self, application_order=None, group_filter=None, energy_cutoff=0.00225, test_modes=None, max_mode_differences=None, maximize_filtered_groups=True, decoupling_overide=100, extra_groups=None, uncoupled_states=0, inconsistent_polyads=None, wavefunction_corrections=None):
        """
        **LLM Docstring**

        Store the shared configuration options used by every concrete `DegeneracySpec` subclass for identifying and filtering groups of degenerate/resonant states.

        :param application_order: when this spec's degeneracy handling should be applied relative to the main perturbation-theory solve (`'pre'`, `'mid'`, or `'post'`); defaults to the class's own `application_order`
        :type application_order: str | None
        :param group_filter: an override for the default group-filtering behavior; `None` to use the default filter, `'unfiltered'` to disable filtering entirely, or an explicit filter spec
        :type group_filter: object | None
        :param energy_cutoff: the zero-order energy difference above which two states are no longer considered candidates for the same degenerate group
        :type energy_cutoff: float
        :param test_modes: restrict degeneracy testing to only these mode indices
        :type test_modes: Iterable[int] | None
        :param max_mode_differences: per-mode caps on how much two states in the same group may differ, used to prune overly broad groups
        :type max_mode_differences: Iterable[int] | None
        :param maximize_filtered_groups: whether the group filter should try to find the largest self-consistent grouping rather than stopping at the first one found
        :type maximize_filtered_groups: bool
        :param decoupling_overide: a coupling-strength override threshold above which two states are kept coupled even if they'd otherwise be split apart by the filter
        :type decoupling_overide: float
        :param extra_groups: additional groups of states to always treat as degenerate, merged in alongside the ones this spec identifies
        :type extra_groups: Iterable | None
        :param uncoupled_states: states that should never be treated as degenerate with anything else, removed from any group they'd otherwise appear in
        :type uncoupled_states: int | Iterable | None
        :param inconsistent_polyads: accepted for interface consistency across subclasses but not used on the base class
        :type inconsistent_polyads: object | None
        :param wavefunction_corrections: perturbation-theory wavefunction corrections used (by some subclasses/filters) to assess coupling strength
        :type wavefunction_corrections: object | None
        :return: None
        :rtype: None
        """
        ...
    repr_opts = ['energy_cutoff']

    def __repr__(self):
        """
        **LLM Docstring**

        Debug string representation listing the class name and the attributes named in `repr_opts`.

        :return: string of the form `ClassName(attr1=val1, attr2=val2, ...)`
        :rtype: str
        """
        ...

    @classmethod
    def merge_state_blocks(cls, state_blocks):
        """
        **LLM Docstring**

        Merge a collection of (possibly overlapping) groups of states into a set of maximal connected groups, via `PermutationRelationGraph.make_relation_graph`.

        :param state_blocks: the groups of states to merge
        :type state_blocks: Iterable
        :return: the merged, maximal groups
        :rtype: list
        """
        ...

    def get_degenerate_group_filter(self, solver=None, evaluator=None, corrs=None, frequencies=None, zero_order_energies=None, high_frequency_modes=None, logger=None, low_frequency_mode_cutoff=0.00115, threshold=None):
        """
        **LLM Docstring**

        Resolve the group-filter specification to actually use when identifying degenerate groups, building the default filter's option dict (energy cutoff, frequencies, target modes, etc.) from this spec's settings (and any solver/evaluator context) unless an explicit override (`None`/`'unfiltered'`, or a custom filter) is already configured.

        :param solver: the perturbation-theory solver providing zero-order energies/high-frequency-mode info, if available
        :type solver: object | None
        :param evaluator: an alternate source of frequencies/energies, if `solver` isn't given
        :type evaluator: object | None
        :param corrs: wavefunction-correction data used to assess coupling strength; defaults to `self.wavefunction_corrections`
        :type corrs: object | None
        :param frequencies: the mode frequencies; inferred from `self.frequencies`/`evaluator.freqs` if not given
        :type frequencies: np.ndarray | None
        :param zero_order_energies: the zero-order state energies; inferred from `solver`/`evaluator`+`corrs` if not given
        :type zero_order_energies: np.ndarray | None
        :param high_frequency_modes: which modes count as "high frequency" for the default `test_modes`; inferred from `solver`/`frequencies` if not given
        :type high_frequency_modes: Iterable[int] | None
        :param logger: logger to pass through to the constructed filter options
        :type logger: Logger | None
        :param low_frequency_mode_cutoff: the frequency threshold used to infer `high_frequency_modes` when not given explicitly
        :type low_frequency_mode_cutoff: float
        :param threshold: an explicit coupling-strength threshold to pass through to the default filter
        :type threshold: float | None
        :return: the resolved group-filter specification (an options dict for the default filter, `None` for no filtering, or the user-supplied filter)
        :rtype: dict | None | object
        """
        ...
    format = None

    @classmethod
    def get_format_mapping(cls):
        """
        **LLM Docstring**

        The mapping from each `DegenerateSpaceInputFormat` enum value to the concrete `DegeneracySpec` subclass that implements it.

        :return: the format-to-class mapping
        :rtype: dict
        """
        ...

    @classmethod
    def get_default_spec(cls):
        """
        **LLM Docstring**

        The default degeneracy-handling spec to use when none is given explicitly (a `StronglyCoupledDegeneracySpec` with default settings).

        :return: the default spec instance
        :rtype: StronglyCoupledDegeneracySpec
        """
        ...

    @classmethod
    def from_spec(cls, spec, format=None, **kwargs):
        """
        **LLM Docstring**

        Build a concrete `DegeneracySpec` instance from a raw specification, inferring which format it represents (via `infer_format`) if not given explicitly, and treating `False`/`None` as "no degeneracy handling" and `True`/`'auto'` as the default spec.

        :param spec: the raw degeneracy specification to build from
        :type spec: object
        :param format: an explicit `DegenerateSpaceInputFormat` to use instead of inferring one
        :type format: DegenerateSpaceInputFormat | None
        :param kwargs: extra options forwarded to the resolved subclass's constructor, merged with any options inferred alongside the spec
        :type kwargs: dict
        :return: the constructed spec, or `None` if `spec` is falsy
        :rtype: DegeneracySpec | None
        """
        ...

    @staticmethod
    def _key(obj, key):
        """
        **LLM Docstring**

        Look up a named attribute or dict-style key on an object, trying attribute access first and falling back to item access, returning `None` if neither succeeds.

        :param obj: the object to look up the key on
        :type obj: object
        :param key: the attribute/key name to look up
        :type key: str
        :return: the found value, or `None` if not found by either means
        :rtype: object | None
        """
        ...

    @classmethod
    def infer_format(cls, spec):
        """
        **LLM Docstring**

        Infer which `DegenerateSpaceInputFormat` a raw specification represents: first checking for an explicit `'format'` attribute/key, then checking whether `spec` exposes a key matching one of the known format value names (pulling any accompanying `'options'`), and finally falling back to asking each registered subclass whether it can `canonicalize` the spec.

        :param spec: the raw specification to classify
        :type spec: object
        :return: `((value, options), format)` -- the extracted value/options pair and the inferred format, or `((None, None), None)` if nothing could be inferred
        :rtype: tuple
        """
        ...

    @abc.abstractmethod
    def get_groups(self, input_states, solver=None, **kwargs):
        """
        :param solver:
        :type solver:
        :param input_states:
        :type input_states:
        :return:
        :rtype:
        """
        ...

    @staticmethod
    def get_group_polyad_relation(exc):
        """
        **LLM Docstring**

        Compute the pairwise "polyad" relation vectors for a set of excitation vectors: for every pair of states, the `(positive part, negative part)` of their quantum-number difference, giving the raising/lowering pattern that would connect them.

        :param exc: the excitation vectors (one row per state) to compute pairwise relations for
        :type exc: np.ndarray
        :return: the array of `(positive, negative)` difference-vector pairs, one per state pair
        :rtype: np.ndarray
        """
        ...

    @classmethod
    def get_polyad_pairs_from_polyad_specs(cls, polyads):
        """
        **LLM Docstring**

        Flatten a list of polyad groups (each a list of transformation rules) into a flat list of all pairwise rule combinations within each polyad.

        :param polyads: the polyad groups, each a list of transformation rules
        :type polyads: list[list]
        :return: the flat list of `[rule1, rule2]` pairs
        :rtype: list[list]
        """
        ...

    def get_polyad_pairs(self, input_states=None, groups=None, solver=None, **kwargs):
        """
        :param solver:
        :type solver:
        :param input_states:
        :type input_states:
        :return:
        :rtype:
        """
        ...

    @classmethod
    @abc.abstractmethod
    def canonicalize(cls, spec):
        """
        **LLM Docstring**

        Abstract hook for validating/normalizing a raw specification into the form this subclass expects, used by `infer_format` to test whether a spec matches this subclass. Concrete subclasses must implement this.

        :param spec: the raw specification to canonicalize
        :type spec: object
        :return: never returns on the base class
        :rtype: object | None
        :raises NotImplementedError: always, on the base class
        """
        ...

    def get_group_filter(self, target_modes=None, max_mode_differences=None, maximize_groups=None, decoupling_overide=None, extra_groups=None, uncoupled_states=None, corrections=None, **kwargs):
        """
        **LLM Docstring**

        Build the group-filter callable used to prune/validate the degenerate groups this spec identifies, filling in this spec's own settings (target modes, mode-difference caps, grouping strategy, etc.) as defaults, via `DegenerateMultiStateSpace.construct_filer`.

        :param target_modes: restrict filtering to these mode indices; defaults to `self.test_modes`
        :type target_modes: Iterable[int] | None
        :param max_mode_differences: per-mode difference caps; defaults to `self.max_mode_differences`
        :type max_mode_differences: Iterable[int] | None
        :param maximize_groups: whether to maximize group size; defaults to `self.maximize_filtered_groups`
        :type maximize_groups: bool | None
        :param decoupling_overide: coupling-strength override threshold; defaults to `self.decoupling_overide`
        :type decoupling_overide: float | None
        :param extra_groups: additional always-degenerate groups; defaults to `self.extra_groups`
        :type extra_groups: Iterable | None
        :param uncoupled_states: states to always keep uncoupled; defaults to `self.uncoupled_states`
        :type uncoupled_states: int | Iterable | None
        :param corrections: wavefunction-correction data; defaults to `self.wavefunction_corrections`
        :type corrections: object | None
        :param kwargs: extra options forwarded to `DegenerateMultiStateSpace.construct_filer`
        :type kwargs: dict
        :return: the constructed group-filter callable
        :rtype: callable
        """
        ...

class EnergyCutoffDegeneracySpec(DegeneracySpec):
    """

    """
    format = DegenerateSpaceInputFormat.EnergyCutoff

    def __init__(self, cutoff, **opts):
        """
        **LLM Docstring**

        Set up a degeneracy spec that groups states purely by zero-order energy proximity.

        :param cutoff: the energy-difference cutoff below which two states are considered degenerate
        :type cutoff: float
        :param opts: extra options forwarded to the base `DegeneracySpec.__init__`
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...

    @classmethod
    def canonicalize(cls, spec):
        """
        **LLM Docstring**

        Check whether a raw specification is a bare numeric energy cutoff, matching this subclass's expected input form.

        :param spec: the raw specification to check
        :type spec: object
        :return: whether `spec` is a plain number
        :rtype: bool
        """
        ...

    @classmethod
    def _group_states_by_energy_cutoff(cls, H0, states, cutoff):
        """
        :type H: Iterable[SparseArray]
        :type states: BasisStateSpace
        :type total_state_space: BasisMultiStateSpace
        :type cutoff: float
        :rtype: Iterable[BasisStateSpace]
        """
        ...

class MartinTestDegeneracySpec(DegeneracySpec):
    application_order = 'mid'
    group_filter = 'default'

    def __init__(self, threshold=4.6e-06, test_energy_window=0.0046, convert=True, frequencies=None, **opts):
        """
        **LLM Docstring**

        Set up a degeneracy spec based on the Martin resonance test (flagging strong couplings via a `|<i|H|j>|^4 / dE^3` criterion), converting the threshold from wavenumbers to Hartrees if it looks like it was given in wavenumbers.

        :param threshold: the Martin-test threshold value that flags a resonance
        :type threshold: float
        :param test_energy_window: the zero-order energy window within which pairs of states are tested for resonance
        :type test_energy_window: float
        :param convert: whether to auto-convert `threshold` from wavenumbers to Hartrees if it appears to be given in wavenumbers (greater than `1e-2`)
        :type convert: bool
        :param frequencies: the mode frequencies used to compute zero-order energies
        :type frequencies: np.ndarray | None
        :param opts: extra options forwarded to the base `DegeneracySpec.__init__`
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...
    repr_opts = ['energy_cutoff', 'threshold']

    def prep_states(self, states: BasisStateSpace):
        """
        **LLM Docstring**

        Identify, for each input state, the nearby (within `self.window`) states reachable via a cubic (`x^3`) coupling operator, storing the resulting candidate groups for later use by `get_coupled_spaces` and returning an expanded state space including any newly discovered candidates.

        :param states: the input states to test
        :type states: BasisStateSpace
        :return: the (possibly expanded) state space including any newly identified candidate states
        :rtype: BasisStateSpace
        """
        ...

    def get_coupled_spaces(self, input_states: BasisStateSpace, solver=None):
        """
        **LLM Docstring**

        Evaluate the Martin resonance test for each input state against its candidate coupling partners (computed via `prep_states` if not already cached), using the solver's first-order Hamiltonian representation to compute the `|H_ij|^4 / dE_ij^3` test statistic, and return the groups of states that exceed `self.threshold`.

        :param input_states: the states to test
        :type input_states: BasisStateSpace
        :param solver: the perturbation-theory solver providing the first-order Hamiltonian representation and total basis
        :type solver: object | None
        :return: the list of resonant state groups (each including the originating state plus its resonant partners), or `[]` if there are no candidates
        :rtype: list[np.ndarray]
        :raises ValueError: if no `solver` is given
        """
        ...

    def get_groups(self, input_states: BasisStateSpace, solver=None, extra_groups=None, **kwargs):
        """
        **LLM Docstring**

        Build the final degenerate groups by running the Martin test (via `get_coupled_spaces`) and merging the resulting groups with any `extra_groups`.

        :param input_states: the states to identify degenerate groups for
        :type input_states: BasisStateSpace
        :param solver: the perturbation-theory solver, forwarded to `get_coupled_spaces`
        :type solver: object | None
        :param extra_groups: additional groups to merge in; defaults to `self.extra_groups`
        :type extra_groups: Iterable | None
        :param kwargs: extra options, unused
        :type kwargs: dict
        :return: the merged degenerate groups
        :rtype: list[np.ndarray]
        """
        ...

    @classmethod
    def canonicalize(cls, spec):
        """
        **LLM Docstring**

        Check whether a raw specification looks like a `MartinTestDegeneracySpec`, by checking that it has a float `threshold` attribute.

        :param spec: the raw specification to check
        :type spec: object
        :return: whether `spec.threshold` is a `float`
        :rtype: bool
        """
        ...

    def get_group_filter(self, **kwargs):
        """
        **LLM Docstring**

        Build the group-filter callable for this spec, supplying the Martin-test-specific state/basis/correction-matrix context (from the most recent `get_coupled_spaces` call) and threshold settings as defaults.

        :param kwargs: explicit overrides for any of the filter options, taking precedence over this spec's own values
        :type kwargs: dict
        :return: the constructed group-filter callable
        :rtype: callable
        """
        ...

class StronglyCoupledDegeneracySpec(DegeneracySpec):
    """

    """
    format = DegenerateSpaceInputFormat.StrongCouplings
    default_threshold = 0.3

    def __init__(self, wfc_threshold=None, state_filter=None, extend_spaces=True, iterations=None, evaluator=None, **opts):
        """
        **LLM Docstring**

        Set up a degeneracy spec based on directly inspecting the strength of perturbation-theory wavefunction-correction couplings.

        :param wfc_threshold: the coupling-strength threshold above which two states are treated as strongly coupled; `None`/`'auto'` uses `self.default_threshold`
        :type wfc_threshold: float | str | None
        :param state_filter: a custom filter function used when identifying strong couplings; defaults to a discriminating filter based on `PerturbationTheoryCorrections.default_state_filter`
        :type state_filter: callable | None
        :param extend_spaces: whether to allow the coupling search to extend beyond the originally requested state space
        :type extend_spaces: bool
        :param iterations: the number of coupling-search iterations to perform
        :type iterations: int | None
        :param evaluator: an evaluator object used to directly compute wavefunction corrections for the input states, if provided (rather than relying on precomputed `corrs` passed to `get_groups`)
        :type evaluator: object | None
        :param opts: extra options forwarded to the base `DegeneracySpec.__init__`
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...

    @property
    def application_order(self):
        """
        **LLM Docstring**

        Property getter/setter for when this spec's degeneracy handling is applied: `'pre'` if an `evaluator` is attached (so degeneracies can be identified before the main solve), otherwise `'post'`. The setter is a no-op, since the value is always derived from `self.evaluator`.

        :param ord: (setter only) ignored
        :type ord: str
        :return: (getter) `'post'` if `self.evaluator` is `None`, else `'pre'`
        :rtype: str
        """
        ...

    @application_order.setter
    def application_order(self, ord):
        """
        **LLM Docstring**

        Property getter/setter for when this spec's degeneracy handling is applied: `'pre'` if an `evaluator` is attached (so degeneracies can be identified before the main solve), otherwise `'post'`. The setter is a no-op, since the value is always derived from `self.evaluator`.

        :param ord: (setter only) ignored
        :type ord: str
        :return: (getter) `'post'` if `self.evaluator` is `None`, else `'pre'`
        :rtype: str
        """
        ...
    repr_opts = ['energy_cutoff', 'wfc_threshold']

    def prep_states(self, input_states):
        """
        **LLM Docstring**

        Return the input states unchanged; this spec doesn't need to expand the input state space ahead of time (unlike, e.g., `MartinTestDegeneracySpec`).

        :param input_states: the input states
        :type input_states: BasisStateSpace
        :return: `input_states`, unchanged
        :rtype: BasisStateSpace
        """
        ...

    def identify_strong_couplings(self, solver, corrs):
        """
        **LLM Docstring**

        Find the strongly-coupled states for every state in `corrs`, via `PerturbationTheoryCorrections.find_strong_couplings`, using either a custom `state_filter` or a default one based on energy proximity and target modes.

        :param solver: the perturbation-theory solver providing zero-order energies and high-frequency-mode info
        :type solver: object
        :param corrs: the perturbation-theory corrections to search for strong couplings in
        :type corrs: PerturbationTheoryCorrections
        :return: the strong-coupling data, as returned by `find_strong_couplings`
        :rtype: dict
        """
        ...

    def get_input_state_couplings(self, input_states):
        """
        **LLM Docstring**

        Compute (and cache) the strongly-coupled partner states for a set of input states, using `self.evaluator` to compute the necessary wavefunction corrections for any states not already cached.

        :param input_states: the states to find couplings for
        :type input_states: BasisStateSpace
        :return: a mapping from input-state index to its strongly-coupled partner states (only including entries with actual couplings)
        :rtype: dict
        """
        ...
    PTCorrectionsMatrix = collections.namedtuple('PTCorrectionsMatrix', ['initial_states', 'full_basis', 'matrices'])

    @classmethod
    def _prep_wfc_correction_space(cls, pt_corrs):
        """
        **LLM Docstring**

        Assemble a set of per-order wavefunction-correction matrices (one initial state per row, one basis state per column) from perturbation-theory correction data, building a combined basis spanning all the final states that appear across every initial state's corrections.

        :param pt_corrs: the raw correction data, with `initial_states`, `final_states` (one per initial state), and `corrections` (one array per initial state, itself one row per order)
        :type pt_corrs: object
        :return: the assembled `PTCorrectionsMatrix` namedtuple (`initial_states`, `full_basis`, `matrices`)
        :rtype: StronglyCoupledDegeneracySpec.PTCorrectionsMatrix
        """
        ...

    def get_groups(self, input_states, couplings=None, solver=None, extra_groups=None, **kwargs):
        """
        :param input_states:
        :type input_states:
        :param solver:
        :type solver:
        :return:
        :rtype:
        """
        ...

    @classmethod
    def canonicalize(cls, spec):
        """
        **LLM Docstring**

        Check whether a raw specification looks like a `StronglyCoupledDegeneracySpec`: a bare number (the coupling threshold), or a dict mapping states to coupling groups that each validate via `GroupsDegeneracySpec._validate_grp`.

        :param spec: the raw specification to check
        :type spec: object
        :return: the validated spec (unchanged number, or dict with validated group values), or `None` if it doesn't match
        :rtype: float | dict | None
        """
        ...

    @classmethod
    def get_strong_coupling_space(cls, states: BasisStateSpace, couplings: dict, extra_groups=None):
        """
        **LLM Docstring**

        Build the final degenerate groups from a mapping of state-to-coupled-states, merging each state with its coupled partners (and any `extra_groups`) via `PermutationRelationGraph.merge_groups`.

        :param states: the states to build groups for
        :type states: BasisStateSpace
        :param couplings: mapping from state index to its strongly-coupled partner states
        :type couplings: dict
        :param extra_groups: additional groups to merge in
        :type extra_groups: Iterable | None
        :return: the merged degenerate groups
        :rtype: list[np.ndarray]
        """
        ...

    def get_degenerate_group_filter(self, solver=None, evaluator=None, threshold=None, logger=None, frequencies=None, **kwargs):
        """
        **LLM Docstring**

        Build the group-filter specification for this spec, supplying the coupling threshold and an appropriate logger (from `solver`/`evaluator` if not given) before delegating to the base `get_degenerate_group_filter`.

        :param solver: the perturbation-theory solver
        :type solver: object | None
        :param evaluator: an alternate evaluator to pull a logger from if `solver` isn't given
        :type evaluator: object | None
        :param threshold: an explicit coupling threshold; defaults to `self.wfc_threshold`
        :type threshold: float | None
        :param logger: an explicit logger; inferred from `solver`/`evaluator` if not given
        :type logger: Logger | None
        :param frequencies: mode frequencies, forwarded to the base method
        :type frequencies: np.ndarray | None
        :param kwargs: extra options forwarded to the base `get_degenerate_group_filter`
        :type kwargs: dict
        :return: the resolved group-filter specification
        :rtype: dict | None | object
        """
        ...

class GroupsDegeneracySpec(DegeneracySpec):
    """

    """
    format = DegenerateSpaceInputFormat.Groups

    def __init__(self, groups, **opts):
        """
        **LLM Docstring**

        Set up a degeneracy spec from an explicit, pre-computed list of degenerate groups.

        :param groups: the explicit degenerate groups (each an array of excitation vectors)
        :type groups: Iterable[np.ndarray]
        :param opts: extra options forwarded to the base `DegeneracySpec.__init__`
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...

    @staticmethod
    def _validate_grp(grp):
        """
        **LLM Docstring**

        Validate that a candidate group looks like a proper group specification (indexable, with at least one element) and coerce it into an array.

        :param grp: the candidate group to validate
        :type grp: object
        :return: the group, coerced to an `np.ndarray`
        :rtype: np.ndarray
        :raises StopIteration: if `grp` isn't indexable (used as a sentinel for "malformatted")
        """
        ...

    @classmethod
    def canonicalize(cls, spec):
        """
        **LLM Docstring**

        Check whether a raw specification is a list of valid groups, matching this subclass's expected input form.

        :param spec: the raw specification to check
        :type spec: object
        :return: the list of validated groups, or `None` if `spec` doesn't match
        :rtype: list[np.ndarray] | None
        """
        ...

    def get_groups(self, input_states, solver=None, **kwargs):
        """
        :param solver:
        :type solver:
        :param input_states:
        :type input_states:
        :return:
        :rtype:
        """
        ...

class PolyadDegeneracySpec(DegeneracySpec):
    """

    """
    format = DegenerateSpaceInputFormat.Polyads

    def __init__(self, polyads, max_quanta=None, iterations=2, require_converged=False, extra_groups=None, extra_polyads=None, full_group_polyads=True, **opts):
        """
        **LLM Docstring**

        Set up a degeneracy spec based on polyad (raising/lowering transformation rule) relations connecting groups of resonant states.

        :param polyads: the polyad transformation rule(s) -- each a `(raise_pattern, lower_pattern)` pair, or a list of such pairs per independent polyad group
        :type polyads: Iterable
        :param max_quanta: the maximum total quanta allowed in states discovered via the polyad rules
        :type max_quanta: int | None
        :param iterations: the number of times to iteratively apply the polyad rules when building connected groups
        :type iterations: int
        :param require_converged: whether to raise an error if the iterative group-building doesn't converge within `iterations`
        :type require_converged: bool
        :param extra_groups: additional groups of states to always include
        :type extra_groups: Iterable | None
        :param extra_polyads: accepted for interface consistency (used by `TotalQuantaDegeneracySpec`) but not directly used on this class
        :type extra_polyads: Iterable | None
        :param full_group_polyads: whether `get_polyad_pairs` should derive pairwise rules from the actual discovered groups (`True`) or just from the raw `polyads` specification (`False`)
        :type full_group_polyads: bool
        :param opts: extra options forwarded to the base `DegeneracySpec.__init__`
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...

    @staticmethod
    def _validate_rule(pair):
        """
        **LLM Docstring**

        Validate that a candidate polyad rule is a `(raise_pattern, lower_pattern)` pair of equal-length vectors.

        :param pair: the candidate rule to validate
        :type pair: object
        :return: the rule, unchanged
        :rtype: tuple
        :raises StopIteration: if `pair` isn't a valid 2-element, equal-length rule (used as a sentinel for "malformatted")
        """
        ...

    @classmethod
    def _check_rule(cls, pair):
        """
        **LLM Docstring**

        Test whether a candidate polyad rule is valid, via `_validate_rule`, without raising.

        :param pair: the candidate rule to check
        :type pair: object
        :return: whether the rule is valid
        :rtype: bool
        """
        ...

    @classmethod
    def canonicalize(cls, spec):
        """
        **LLM Docstring**

        Check whether a raw specification is a list of valid polyad rules, matching this subclass's expected input form.

        :param spec: the raw specification to check
        :type spec: object
        :return: the list of validated rules, or `None` if `spec` doesn't match
        :rtype: list[tuple] | None
        """
        ...

    def get_groups(self, input_states, solver=None, **kwargs):
        """
        :param solver:
        :type solver:
        :param input_states:
        :type input_states:
        :return:
        :rtype:
        """
        ...

    @classmethod
    def get_degenerate_polyad_space(cls, states, polyadic_pairs, max_quanta=None, iterations=2, require_converged=False, extra_groups=None):
        """
        Gets degenerate spaces by using pairs of transformation rules to
        take an input state and connect it to other degenerate states

        :param states: the input states
        :type states:
        :param polyadic_pairs: the transformation rules
        :type polyadic_pairs:
        :param max_quanta: the max quanta to allow in connected states
        :type max_quanta:
        :return:
        :rtype:
        """
        ...

    def get_polyad_pairs(self, input_states=None, groups=None, solver=None, full_group_polyads=None, **kwargs):
        """
        **LLM Docstring**

        Get the pairwise polyad transformation rules used for identifying strong couplings: either derived from the actual discovered degenerate groups (via the base class's `get_polyad_pairs`), or directly from the raw `polyads` specification, depending on `full_group_polyads`.

        :param input_states: the states to identify groups for, if deriving pairs from actual groups
        :type input_states: BasisStateSpace | None
        :param groups: precomputed groups to use instead of recomputing them
        :type groups: list | None
        :param solver: the perturbation-theory solver, forwarded if deriving pairs from actual groups
        :type solver: object | None
        :param full_group_polyads: whether to derive pairs from actual groups; defaults to `self.full_group_polyads`
        :type full_group_polyads: bool | None
        :param kwargs: extra options forwarded to the base `get_polyad_pairs`
        :type kwargs: dict
        :return: the pairwise polyad rules
        :rtype: list[list]
        """
        ...

    @staticmethod
    def _is_polyad_rule(d, n_modes):
        """
        **LLM Docstring**

        Check whether a candidate object looks like a single polyad rule (a 2-element pair of length-`n_modes` vectors), as opposed to a list of such rules.

        :param d: the candidate object to check
        :type d: object
        :param n_modes: the expected vector length (number of modes)
        :type n_modes: int
        :return: whether `d` looks like a single rule
        :rtype: bool
        """
        ...

class TotalQuantaDegeneracySpec(PolyadDegeneracySpec):
    """

    """
    format = DegenerateSpaceInputFormat.QuantaSpecRules

    def __init__(self, n_T_vectors, max_quanta=3, target_modes=None, extra_polyads=None, **opts):
        """
        **LLM Docstring**

        Set up a degeneracy spec derived from one or more "nT" total-quanta vectors, automatically generating the corresponding polyad transformation rules (via `make_nt_polyad`) that connect states with the same weighted total-quanta value.

        :param n_T_vectors: one or more nT weighting vectors, each giving a per-mode weight used to compute a weighted total-quanta value
        :type n_T_vectors: Iterable[int] | Iterable[Iterable[int]]
        :param max_quanta: the maximum number of quanta to consider when constructing candidate resonant states
        :type max_quanta: int
        :param target_modes: restrict the nT-vector components used to only these modes
        :type target_modes: Iterable[int] | None
        :param extra_polyads: additional raw polyad rules to include alongside the auto-generated ones
        :type extra_polyads: Iterable | None
        :param opts: extra options forwarded to the base `PolyadDegeneracySpec.__init__`
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...

    @classmethod
    def reduce_rule(cls, a, b):
        """
        **LLM Docstring**

        Reduce a pair of raising/lowering quantum-number patterns to their simplest equivalent form: dividing out a common integer factor where possible, and canceling overlapping positive contributions shared between the two patterns.

        :param a: the first quantum-number pattern
        :type a: np.ndarray
        :param b: the second quantum-number pattern
        :type b: np.ndarray
        :return: the reduced `(a, b)` pair
        :rtype: tuple[np.ndarray, np.ndarray]
        """
        ...

    @classmethod
    def make_nt_polyad(cls, nt, target_modes=None, max_quanta=3):
        """
        **LLM Docstring**

        Build the set of pairwise polyad transformation rules connecting states that share the same weighted total-quanta value (`dot(state, nt)`), by enumerating states up to `max_quanta` total quanta, grouping them by their nT value, and forming (and reducing) the raising/lowering rule for every pair within each nonzero-nT group.

        :param nt: the nT weighting vector
        :type nt: Iterable[int]
        :param target_modes: restrict to these mode indices; defaults to the modes where `nt` is nonzero
        :type target_modes: Iterable[int] | None
        :param max_quanta: the maximum number of quanta to enumerate states up to
        :type max_quanta: int
        :return: the list of unique `(raise_pattern, lower_pattern)` rules
        :rtype: list[tuple]
        """
        ...

class CallableDegeneracySpec(DegeneracySpec):
    """

    """
    format = DegenerateSpaceInputFormat.Callable

    def __init__(self, callable, **opts):
        """
        **LLM Docstring**

        Set up a degeneracy spec that delegates entirely to a user-supplied callable for identifying degenerate groups.

        :param callable: the function used to compute degenerate groups, called as `callable(input_states, solver=solver, **kwargs)`
        :type callable: callable
        :param opts: extra options forwarded to the base `DegeneracySpec.__init__`
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...

    def get_groups(self, input_states, solver=None, **kwargs):
        """
        :param input_states:
        :type input_states:
        :param solver:
        :type solver:
        :return:
        :rtype:
        """
        ...

class DegenerateMultiStateSpace(BasisMultiStateSpace):

    @staticmethod
    def default_group_filter(group, states=None, basis=None, corrections=None, threshold=None, energy_cutoff=None, frequencies=None, energies=None, decoupling_overide=10, maximize_groups=True, target_modes=None, max_mode_differences=None, threshold_step_size=0.1, extra_groups=None, uncoupled_states=None, state_diff_filters=1, logger=None):
        """
        Excludes modes that differ in only one position, prioritizing states with fewer numbers of quanta
        (potentially add restrictions to high frequency modes...?)

        :param input_state:
        :type input_state:
        :param couplings:
        :type couplings:
        :return:
        :rtype:
        """
        ...

    @classmethod
    def construct_filer(cls, spec=None, **kwargs):
        """
        **LLM Docstring**

        Build the group-filter callable to use when identifying degenerate groups: if a `DegeneracySpec` is given, delegates to its own `get_group_filter`; otherwise wraps `default_group_filter` with the given keyword options bound in.

        :param spec: an explicit `DegeneracySpec` to build the filter from, if available
        :type spec: DegeneracySpec | None
        :param kwargs: options forwarded to `spec.get_group_filter` or bound into the default filter
        :type kwargs: dict
        :return: the constructed group-filter callable
        :rtype: callable
        """
        ...

    @classmethod
    def from_spec(cls, degenerate_states, states=None, solver=None, evaluator=None, full_basis=None, format=None, group_filter=None, log_groups=False, **kwargs):
        """
        Generates a DegenerateMultiStateSpace object from a number
        of possible specs

        :param solver: the actual applier of the perturbation theory which makes use of the degenerate states
        :type solver: PerturbationTheorySolver
        :return:
        :rtype:
        """
        ...