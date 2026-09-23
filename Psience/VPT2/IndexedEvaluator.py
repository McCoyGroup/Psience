"""Concrete indexed-block execution for analytic perturbation expressions.

This module is intentionally independent of the symbolic expression builder.
``Analytic.py`` imports it lazily after an expression has been normalized, which
keeps the symbolic DAG and concrete numerical backend as separate layers.

The first implementation establishes the backend boundary and a bounded,
immutable structural plan.  It joins all coefficient/state work items produced
for one free-index subset before performing denominator contraction and result
reduction.  Later stages can replace the coefficient-row producer without
changing the plan or contraction interfaces defined here.
"""

import collections
import dataclasses
import time
import types

import numpy as np


_NO_SELECTION = object()


class IndexedDegeneracyEvaluationContext:
    """Batch-local degeneracy masks over one shared permutation pool.

    Every concrete work item produced by one ``eval_perm`` call is routed from
    the same ``all_permutations`` array.  Evaluating a predicate against each
    private one- or two-row work-item subset defeats the array-oriented
    ``DegeneracyTestPlan`` and repeats identical tests.  This context evaluates
    each predicate, and each left/right join, once over the shared projected
    pool.  Concrete blocks recover their local masks by integer indexing.

    The cache is deliberately scoped to one ``eval_perm`` call and has a byte
    ceiling.  A large or unusual custom predicate can therefore fall back to
    recomputation without growing a runner-lifetime cache.
    """

    default_max_cache_bytes = 32 * 1024 ** 2

    def __init__(self, permutations, stats, max_cache_bytes=None):
        permutations = np.asarray(permutations, dtype=np.intp)
        if permutations.ndim != 2:
            raise ValueError("projected permutations must be two-dimensional")
        self.permutations = permutations
        self.stats = stats
        self.max_cache_bytes = (
            self.default_max_cache_bytes
            if max_cache_bytes is None else max(0, int(max_cache_bytes))
        )
        self.predicate_masks = {}
        self.selection_masks = {}
        self.block_selections = {}
        self.cache_bytes = 0
        stats.degeneracy_pool_contexts += 1
        stats.degeneracy_pool_rows += len(permutations)

    @staticmethod
    def _mask_bytes(mask):
        return mask.nbytes if isinstance(mask, np.ndarray) else 0

    def _store(self, cache, key, entry, mask=None, size=None):
        if size is None:
            size = self._mask_bytes(mask)
        if size > self.max_cache_bytes - self.cache_bytes:
            self.stats.degeneracy_pool_cache_skips += 1
            return
        cache[key] = entry
        self.cache_bytes += size
        self.stats.max_degeneracy_pool_cache_bytes = max(
            self.stats.max_degeneracy_pool_cache_bytes,
            self.cache_bytes
        )

    def predicate_mask(self, predicate, only_degenerate):
        if predicate is None:
            return not only_degenerate

        predicate_key = id(predicate)
        entry = self.predicate_masks.get(predicate_key)
        if entry is not None and entry[0] is predicate:
            self.stats.degeneracy_pool_predicate_hits += 1
            mask = entry[1]
        else:
            self.stats.degeneracy_pool_predicate_misses += 1
            tick = time.perf_counter()
            evaluate_array = getattr(predicate, 'evaluate_array', None)
            self.stats.degeneracy_states_tested += len(self.permutations)
            array_cutoff = getattr(predicate, 'scalar_state_cutoff', -1)
            if (
                    evaluate_array is None
                    or len(self.permutations) <= array_cutoff
            ):
                self.stats.degeneracy_scalar_evaluations += 1
                mask = np.fromiter(
                    (
                        bool(predicate(permutation))
                        for permutation in self.permutations
                    ),
                    dtype=bool,
                    count=len(self.permutations)
                )
            else:
                self.stats.degeneracy_array_evaluations += 1
                self.stats.degeneracy_patterns_tested += len(
                    predicate.patterns
                )
                mask = np.asarray(
                    evaluate_array(self.permutations),
                    dtype=bool
                )
            self.stats.degeneracy_seconds += time.perf_counter() - tick
            self._store(
                self.predicate_masks,
                predicate_key,
                (predicate, mask),
                mask
            )
        return mask if only_degenerate else np.logical_not(mask)

    def selection_mask(self, left_predicate, right_predicate,
                       use_left, use_right, join):
        key = (
            id(left_predicate), id(right_predicate),
            bool(use_left), bool(use_right), id(join)
        )
        entry = self.selection_masks.get(key)
        if (
                entry is not None
                and entry[0] is left_predicate
                and entry[1] is right_predicate
                and entry[2] is join
        ):
            self.stats.degeneracy_pool_selection_hits += 1
            return entry[3]

        self.stats.degeneracy_pool_selection_misses += 1
        left = self.predicate_mask(left_predicate, use_left)
        right = self.predicate_mask(right_predicate, use_right)
        selected = join(left, right)
        if np.ndim(selected) > 0:
            selected = np.asarray(selected, dtype=bool)
        self._store(
            self.selection_masks,
            key,
            (left_predicate, right_predicate, join, selected),
            selected
        )
        return selected

    def lookup_block_selection(self, selection_plan, global_indices):
        global_indices = np.asarray(global_indices, dtype=np.intp)
        index_key = global_indices.tobytes()
        key = (id(selection_plan), len(global_indices), index_key)
        entry = self.block_selections.get(key)
        if entry is not None and entry[0] is selection_plan:
            self.stats.degeneracy_block_selection_hits += 1
            return key, entry[1], entry[2]
        self.stats.degeneracy_block_selection_misses += 1
        return key, None, None

    def store_block_selection(
            self, key, selection_plan, resolved, selection_classes):
        size = len(key[2]) + sum(
            selection.nbytes
            for selection in resolved
            if isinstance(selection, np.ndarray)
        )
        before = len(self.block_selections)
        self._store(
            self.block_selections,
            key,
            (selection_plan, resolved, selection_classes),
            size=size
        )
        if len(self.block_selections) == before:
            self.stats.degeneracy_block_selection_cache_skips += 1


@dataclasses.dataclass(frozen=True)
class IndexedEnergyTerm:
    """One polynomial numerator and its denominator change product."""

    ordinal: int
    energy_changes: tuple
    polynomial: object
    contraction_shape: tuple


@dataclasses.dataclass(frozen=True)
class IndexedSubexpressionPlan:
    """Static numerical plan for one coefficient-product expression."""

    expression: object
    direct_polynomial: bool
    prefactor: object
    groups: tuple


@dataclasses.dataclass(frozen=True)
class IndexedDegeneracySelectionPlan:
    """Runtime binding from static energy terms to degeneracy predicates.

    ``degenerate_changes`` is invariant across the many concrete blocks in one
    indexed evaluator.  Binding every energy term to the same left/right
    predicate pair inside ``_select_energy_groups`` repeated millions of
    dictionary lookups in large degenerate calculations.  This plan performs
    those lookups once and retains the original term order for exact reduction
    compatibility with the legacy selector.
    """

    subexpression: object
    left_changes: object
    right_changes: object
    use_left: bool
    use_right: bool
    join: object
    selections: tuple
    shape_terms: tuple
    term_count: int


@dataclasses.dataclass(frozen=True)
class IndexedExpressionPlan:
    """Immutable structural plan for a ``PTTensorCoeffProductSum``."""

    expression: object
    subexpressions: object
    term_count: int


class ConcreteBlockTable:
    """Columnar storage for concrete state/permutation destinations.

    A full evaluation can create millions of these variable-sized records.  A
    Python object per record both consumes memory and makes record construction
    a significant part of wall time.  Keeping one list per field retains the
    irregular arrays required by the evaluator without allocating a separate
    object and attribute dictionary/slots tuple for every block.
    """

    __slots__ = (
        'state_indices', 'mask_positions', 'global_permutation_indices',
        'prefactors', 'expression_keys',
        'plans', 'permutation_subsets', 'permutation_substates',
        'tuple_states', 'permutation_frequencies', 'state_block_identities',
        'contributions'
    )

    def __init__(self):
        self.state_indices = []
        self.mask_positions = []
        self.global_permutation_indices = []
        self.prefactors = []
        self.expression_keys = []
        self.plans = []
        self.permutation_subsets = []
        self.permutation_substates = []
        self.tuple_states = []
        self.permutation_frequencies = []
        self.state_block_identities = []
        self.contributions = []

    def append(self, state_index, mask_positions, global_permutation_indices,
               prefactors, expression_key,
               plan, permutation_subsets, permutation_substates, tuple_states,
               permutation_frequencies, state_block_identity, contribution):
        index = len(self.state_indices)
        self.state_indices.append(state_index)
        self.mask_positions.append(mask_positions)
        self.global_permutation_indices.append(global_permutation_indices)
        self.prefactors.append(prefactors)
        self.expression_keys.append(expression_key)
        self.plans.append(plan)
        self.permutation_subsets.append(permutation_subsets)
        self.permutation_substates.append(permutation_substates)
        self.tuple_states.append(tuple_states)
        self.permutation_frequencies.append(permutation_frequencies)
        self.state_block_identities.append(state_block_identity)
        self.contributions.append(contribution)
        return index

    def __len__(self):
        return len(self.state_indices)


class IndexedContractionSegmentTable:
    """Compact rows for one homogeneous contraction kernel signature.

    Each segment says that all ``terms`` use the same concrete block and state
    selection.  Rows are materialized only into a bounded kernel-size buffer,
    avoiding a persistent ``(block, selection, term)`` tuple for every row.
    """

    __slots__ = ('block_indices', 'selections', 'term_runs', 'row_count')

    def __init__(self):
        self.block_indices = []
        self.selections = []
        self.term_runs = []
        self.row_count = 0

    def append(self, block_index, selection, terms):
        if len(terms) == 0:
            return
        self.block_indices.append(block_index)
        self.selections.append(selection)
        self.term_runs.append(terms)
        self.row_count += len(terms)

    def iter_blocks(self, block_size, return_segment_ids=False):
        """Yield row columns while retaining no more than ``block_size`` rows."""
        segment_index = 0
        term_offset = 0
        segment_count = len(self.block_indices)
        while segment_index < segment_count:
            block_indices = np.empty(block_size, dtype=np.intp)
            segment_ids = np.empty(block_size, dtype=np.intp)
            selections = [None] * block_size
            terms = [None] * block_size
            output_offset = 0
            while output_offset < block_size and segment_index < segment_count:
                run = self.term_runs[segment_index]
                take = min(block_size - output_offset, len(run) - term_offset)
                stop = output_offset + take
                block_indices[output_offset:stop] = self.block_indices[segment_index]
                segment_ids[output_offset:stop] = segment_index
                selections[output_offset:stop] = [
                    self.selections[segment_index]
                ] * take
                terms[output_offset:stop] = run[term_offset:term_offset + take]
                output_offset = stop
                term_offset += take
                if term_offset == len(run):
                    segment_index += 1
                    term_offset = 0
            block = (
                block_indices[:output_offset],
                selections[:output_offset],
                terms[:output_offset]
            )
            if return_segment_ids:
                block = block + (segment_ids[:output_offset],)
            yield block


@dataclasses.dataclass(frozen=True)
class IndexedPermutationRoutes:
    """Matched coefficient rows grouped by coefficient and destination state."""

    coefficient_indices: object
    state_indices: object
    mask_positions: object
    value_positions: object
    group_starts: object

    def __len__(self):
        return max(0, len(self.group_starts) - 1)

    def __iter__(self):
        for group_index in range(len(self)):
            start = self.group_starts[group_index]
            stop = self.group_starts[group_index + 1]
            yield (
                int(self.coefficient_indices[start]),
                int(self.state_indices[start]),
                self.mask_positions[start:stop],
                self.value_positions[start:stop]
            )


@dataclasses.dataclass(frozen=True)
class IndexedPermutationRoutingPlan:
    """CSR-style global-permutation to state/local-position routing.

    ``permutation_map`` can contain a global permutation more than once when
    it is needed by multiple states.  Compiling those occurrences once lets
    each coefficient row expand directly to its actual destinations instead
    of rebuilding a full global mask and splitting it for every coefficient.
    """

    global_offsets: object
    destination_positions: object
    state_offsets: object

    @classmethod
    def compile(cls, permutation_map, state_lengths, permutation_count):
        permutation_map = np.asarray(permutation_map, dtype=np.intp)
        state_lengths = np.asarray(state_lengths, dtype=np.intp)
        permutation_count = int(permutation_count)
        if np.sum(state_lengths, dtype=np.intp) != len(permutation_map):
            raise ValueError(
                "state permutation lengths do not match permutation map"
            )
        if len(permutation_map) > 0 and (
            np.min(permutation_map) < 0
            or np.max(permutation_map) >= permutation_count
        ):
            raise ValueError("permutation map contains an out-of-range index")

        state_offsets = np.empty(len(state_lengths) + 1, dtype=np.intp)
        state_offsets[0] = 0
        np.cumsum(state_lengths, out=state_offsets[1:])

        counts = np.bincount(
            permutation_map,
            minlength=permutation_count
        ).astype(np.intp, copy=False)
        global_offsets = np.empty(permutation_count + 1, dtype=np.intp)
        global_offsets[0] = 0
        np.cumsum(counts, out=global_offsets[1:])
        destination_positions = np.argsort(
            permutation_map,
            kind='stable'
        ).astype(np.intp, copy=False)
        return cls(global_offsets, destination_positions, state_offsets)

    def route(self, eval_permutations, eval_coefficients):
        """Join evaluated coefficient rows to only their real destinations."""
        eval_permutations = np.asarray(eval_permutations, dtype=np.intp)
        eval_coefficients = np.asarray(eval_coefficients, dtype=np.intp)
        if len(eval_permutations) != len(eval_coefficients):
            raise ValueError(
                "coefficient and permutation row counts do not match"
            )
        if len(eval_permutations) == 0:
            empty = np.empty(0, dtype=np.intp)
            return IndexedPermutationRoutes(
                empty, empty, empty, empty,
                np.zeros(1, dtype=np.intp)
            )

        occurrence_starts = self.global_offsets[eval_permutations]
        occurrence_counts = (
            self.global_offsets[eval_permutations + 1] - occurrence_starts
        )
        route_count = int(np.sum(occurrence_counts, dtype=np.intp))
        if route_count == 0:
            empty = np.empty(0, dtype=np.intp)
            return IndexedPermutationRoutes(
                empty, empty, empty, empty,
                np.zeros(1, dtype=np.intp)
            )

        value_positions = np.repeat(
            np.arange(len(eval_permutations), dtype=np.intp),
            occurrence_counts
        )
        output_starts = np.empty(len(occurrence_counts), dtype=np.intp)
        output_starts[0] = 0
        if len(output_starts) > 1:
            np.cumsum(
                occurrence_counts[:-1],
                out=output_starts[1:]
            )
        occurrence_positions = np.arange(route_count, dtype=np.intp)
        occurrence_positions -= np.repeat(output_starts, occurrence_counts)
        occurrence_positions += occurrence_starts[value_positions]
        destination_positions = self.destination_positions[
            occurrence_positions
        ]
        coefficient_indices = eval_coefficients[value_positions]

        # This reproduces the legacy coefficient-group order followed by
        # ascending state/local permutation position without scanning the
        # entire permutation map for each coefficient.
        sorting = np.lexsort((destination_positions, coefficient_indices))
        destination_positions = destination_positions[sorting]
        value_positions = value_positions[sorting]
        coefficient_indices = coefficient_indices[sorting]
        state_indices = np.searchsorted(
            self.state_offsets[1:],
            destination_positions,
            side='right'
        ).astype(np.intp, copy=False)
        mask_positions = (
            destination_positions - self.state_offsets[state_indices]
        )

        group_changes = np.logical_or(
            coefficient_indices[1:] != coefficient_indices[:-1],
            state_indices[1:] != state_indices[:-1]
        )
        group_starts = np.concatenate((
            np.array([0], dtype=np.intp),
            np.flatnonzero(group_changes) + 1,
            np.array([route_count], dtype=np.intp)
        ))
        return IndexedPermutationRoutes(
            coefficient_indices,
            state_indices,
            mask_positions,
            value_positions,
            group_starts
        )


@dataclasses.dataclass
class IndexedEvaluationStats:
    plan_cache_hits: int = 0
    plan_cache_misses: int = 0
    plans_compiled: int = 0
    plan_terms: int = 0
    calls: int = 0
    work_items: int = 0
    contraction_groups: int = 0
    contraction_blocks: int = 0
    contraction_rows: int = 0
    contraction_segments: int = 0
    polynomial_evaluations: int = 0
    coefficient_plans_compiled: int = 0
    coefficient_plan_hits: int = 0
    coefficient_plan_misses: int = 0
    coefficient_product_chunks: int = 0
    coefficient_tensor_groups: int = 0
    coefficient_gather_requests: int = 0
    coefficient_pattern_cache_hits: int = 0
    coefficient_pattern_cache_misses: int = 0
    coefficient_pattern_cache_evictions: int = 0
    coefficient_pattern_compilations: int = 0
    coefficient_unique_patterns: int = 0
    coefficient_rows: int = 0
    coefficient_seconds: float = 0.0
    degeneracy_mask_hits: int = 0
    degeneracy_mask_misses: int = 0
    degeneracy_mask_evictions: int = 0
    degeneracy_selection_hits: int = 0
    degeneracy_selection_misses: int = 0
    degeneracy_array_evaluations: int = 0
    degeneracy_scalar_evaluations: int = 0
    degeneracy_states_tested: int = 0
    degeneracy_patterns_tested: int = 0
    degeneracy_seconds: float = 0.0
    degeneracy_pool_contexts: int = 0
    degeneracy_pool_rows: int = 0
    degeneracy_pool_predicate_hits: int = 0
    degeneracy_pool_predicate_misses: int = 0
    degeneracy_pool_selection_hits: int = 0
    degeneracy_pool_selection_misses: int = 0
    degeneracy_pool_cache_skips: int = 0
    max_degeneracy_pool_cache_bytes: int = 0
    degeneracy_block_selection_hits: int = 0
    degeneracy_block_selection_misses: int = 0
    degeneracy_block_selection_cache_skips: int = 0
    degeneracy_selection_plan_hits: int = 0
    degeneracy_selection_plan_misses: int = 0
    degeneracy_selection_plan_evictions: int = 0
    degeneracy_selection_descriptors: int = 0
    degeneracy_group_plan_hits: int = 0
    degeneracy_group_plan_misses: int = 0
    degeneracy_group_plan_evictions: int = 0
    max_degeneracy_selection_plan_cache_items: int = 0
    max_degeneracy_group_plan_cache_items: int = 0
    energy_seconds: float = 0.0
    polynomial_seconds: float = 0.0
    reduction_seconds: float = 0.0
    reduction_scatters: int = 0
    elapsed_seconds: float = 0.0
    max_rows_per_block: int = 0
    max_pending_segments: int = 0
    max_polynomial_scratch_bytes: int = 0
    max_coefficient_pattern_cache_items: int = 0
    max_coefficient_pattern_cache_bytes: int = 0
    routing_plans_compiled: int = 0
    routing_matches: int = 0
    routing_seconds: float = 0.0

    def as_dict(self):
        return dataclasses.asdict(self)


class IndexedExpressionPlanCache:
    """Small FIFO cache retaining expressions alongside identity-keyed plans."""

    def __init__(self, max_items=128):
        self.max_items = max(0, int(max_items))
        self._plans = collections.OrderedDict()

    def get(self, expression):
        entry = self._plans.get(id(expression))
        if entry is None or entry[0] is not expression:
            return None
        return entry[1]

    def store(self, expression, plan):
        if self.max_items == 0:
            return
        key = id(expression)
        self._plans[key] = (expression, plan)
        self._plans.move_to_end(key)
        while len(self._plans) > self.max_items:
            self._plans.popitem(last=False)

    def clear(self):
        self._plans.clear()

    def __len__(self):
        return len(self._plans)


class IndexedRequestPatternCache:
    """Bounded LRU of structural unique-pattern and inverse arrays."""

    def __init__(self, max_items=256, max_bytes=256 * 1024):
        self.max_items = max(0, int(max_items))
        self.max_bytes = max(0, int(max_bytes))
        self._entries = collections.OrderedDict()
        self.bytes = 0
        self.peak_items = 0
        self.peak_bytes = 0
        self.evictions = 0

    @staticmethod
    def _entry_size(key, unique_patterns, inverse):
        return (
            unique_patterns.nbytes
            + inverse.nbytes
            + sum(len(pattern) for pattern in key) * np.dtype(np.intp).itemsize
        )

    def get(self, key):
        entry = self._entries.get(key)
        if entry is not None:
            self._entries.move_to_end(key)
        return entry

    def store(self, key, unique_patterns, inverse):
        if self.max_items == 0 or self.max_bytes == 0:
            return
        size = self._entry_size(key, unique_patterns, inverse)
        if size > self.max_bytes:
            return
        prior = self._entries.pop(key, None)
        if prior is not None:
            self.bytes -= prior[2]
        self._entries[key] = (unique_patterns, inverse, size)
        self.bytes += size
        while (
                len(self._entries) > self.max_items
                or self.bytes > self.max_bytes
        ):
            _, evicted = self._entries.popitem(last=False)
            self.bytes -= evicted[2]
            self.evictions += 1
        self.peak_items = max(self.peak_items, len(self._entries))
        self.peak_bytes = max(self.peak_bytes, self.bytes)

    def __len__(self):
        return len(self._entries)


@dataclasses.dataclass(frozen=True)
class IndexedCoefficientGatherPlan:
    """Static coefficient products paired with their concrete tensor sets."""

    coefficient_specs: object
    coefficient_tensors: object
    expansion_count: int
    max_workspace_bytes: int = 32 * 1024 ** 2
    # Structural request sets stay small and do not scale with mode count.
    # Keep the per-coefficient-plan bound conservative because one expression
    # can retain many coefficient plans at once.
    max_pattern_cache_items: int = 256
    max_pattern_cache_bytes: int = 256 * 1024
    request_pattern_cache: object = dataclasses.field(
        init=False,
        repr=False,
        compare=False
    )

    use_cached_request_patterns = True

    def __post_init__(self):
        object.__setattr__(
            self,
            'request_pattern_cache',
            IndexedRequestPatternCache(
                self.max_pattern_cache_items,
                self.max_pattern_cache_bytes
            )
        )

    @staticmethod
    def _remap_product(product, subset, num_fixed):
        return tuple(
            (
                coefficient[:2],
                tuple(
                    subset[index - num_fixed]
                    if index >= num_fixed else index
                    for index in coefficient[2:]
                )
            )
            for coefficient in product
        )

    def _active_products(self, subset, num_fixed, counts_cache, factorials):
        remapped = []
        active = []
        for product_index, product in enumerate(self.coefficient_specs):
            concrete = self._remap_product(product, subset, num_fixed)
            hits = counts_cache.get(concrete, 0)
            if hits < factorials[len(concrete)]:
                counts_cache[concrete] = hits + 1
                active.append(product_index)
                remapped.append(concrete)
        return np.asarray(active, dtype=int), remapped

    def _get_request_pattern_map(self, patterns, stats):
        """Compile structural uniqueness once, independently of mode subset."""
        # Normal plan inputs are already immutable slices of coefficient
        # specifications.  Preserve those tuples so a cache hit only allocates
        # the outer grouping key; normalize array/list inputs for direct use.
        key = tuple(
            pattern if isinstance(pattern, tuple) else tuple(pattern)
            for pattern in patterns
        )
        entry = self.request_pattern_cache.get(key)
        if entry is None:
            stats.coefficient_pattern_cache_misses += 1
            stats.coefficient_pattern_compilations += 1
            unique_patterns, inverse = np.unique(
                np.asarray(key, dtype=np.intp),
                axis=0,
                return_inverse=True
            )
            unique_patterns.setflags(write=False)
            inverse = inverse.astype(np.intp, copy=False)
            inverse.setflags(write=False)
            self.request_pattern_cache.store(
                key,
                unique_patterns,
                inverse
            )
        else:
            stats.coefficient_pattern_cache_hits += 1
            unique_patterns, inverse, _ = entry
        cache = self.request_pattern_cache
        stats.max_coefficient_pattern_cache_items = max(
            stats.max_coefficient_pattern_cache_items,
            cache.peak_items
        )
        stats.max_coefficient_pattern_cache_bytes = max(
            stats.max_coefficient_pattern_cache_bytes,
            cache.peak_bytes
        )
        stats.coefficient_unique_patterns += len(unique_patterns)
        return unique_patterns, inverse

    @staticmethod
    def _remap_request_patterns(patterns, subset, num_fixed):
        patterns = np.array(patterns, dtype=np.intp, copy=True)
        free = patterns >= num_fixed
        if np.any(free):
            subset = np.asarray(subset, dtype=np.intp)
            patterns[free] = subset[patterns[free] - num_fixed]
        return patterns

    @staticmethod
    def _gather_tensor_requests(permutations, tensor,
                                unique_patterns, inverse):
        """Evaluate precompiled distinct patterns for one coefficient tensor."""
        # permutation x request x tensor-axis
        tensor_indices = permutations[:, unique_patterns]
        gathered = tensor[tuple(
            tensor_indices[..., axis]
            for axis in range(tensor_indices.shape[-1])
        )]
        # request x permutation, restored to the possibly repeated request list
        return np.moveaxis(gathered, 0, 1)[inverse]

    def _evaluate_product_chunk(self, permutations, product_indices,
                                remapped_products, subset, num_fixed,
                                zero_cutoff, stats):
        product_count = len(product_indices)
        permutation_count = len(permutations)
        expansion_count = self.expansion_count
        products = np.ones(
            (product_count, expansion_count, permutation_count),
            dtype=float
        )
        good = np.ones(
            (product_count, expansion_count, permutation_count),
            dtype=bool
        )
        max_factors = max(
            (len(product) for product in remapped_products),
            default=0
        )

        for factor_index in range(max_factors):
            stage_products = np.asarray([
                local_index
                for local_index, product in enumerate(remapped_products)
                if factor_index < len(product)
            ], dtype=int)
            if len(stage_products) == 0:
                continue

            for expansion_index in range(expansion_count):
                request_groups = collections.OrderedDict()
                scalar_requests = []
                for local_index in stage_products:
                    global_index = int(product_indices[local_index])
                    _, indices = remapped_products[local_index][factor_index]
                    tensor = self.coefficient_tensors[
                        global_index
                    ][expansion_index][factor_index]
                    if len(indices) == 0 or np.ndim(tensor) == 0:
                        scalar_requests.append((local_index, tensor, len(indices)))
                    else:
                        tensor_key = id(tensor)
                        entry = request_groups.get(tensor_key)
                        if entry is None or entry[0] is not tensor:
                            entry = [tensor, [], [], []]
                            request_groups[tensor_key] = entry
                        entry[1].append(local_index)
                        entry[2].append(indices)
                        entry[3].append(
                            self.coefficient_specs[
                                global_index
                            ][factor_index][2:]
                        )

                for local_index, tensor, index_count in scalar_requests:
                    base_value = np.asanyarray(tensor).item()
                    # A structurally absent tensor with nonempty indices is a
                    # zero coefficient for this expansion.
                    if index_count > 0:
                        base_value = 0
                    if abs(base_value) < zero_cutoff:
                        # Match the legacy scalar branch: one rejected scalar
                        # invalidates every expansion for this product.
                        good[local_index, :, :] = False
                    else:
                        products[local_index, expansion_index, :] *= base_value

                for tensor, destinations, patterns, structural_patterns \
                        in request_groups.values():
                    stats.coefficient_tensor_groups += 1
                    stats.coefficient_gather_requests += len(patterns)
                    if self.use_cached_request_patterns:
                        old_evictions = self.request_pattern_cache.evictions
                        unique_structural, inverse = (
                            self._get_request_pattern_map(
                                structural_patterns,
                                stats
                            )
                        )
                        stats.coefficient_pattern_cache_evictions += (
                            self.request_pattern_cache.evictions
                            - old_evictions
                        )
                        unique_patterns = self._remap_request_patterns(
                            unique_structural,
                            subset,
                            num_fixed
                        )
                    else:
                        stats.coefficient_pattern_compilations += 1
                        unique_patterns, inverse = np.unique(
                            np.asarray(patterns, dtype=np.intp),
                            axis=0,
                            return_inverse=True
                        )
                        stats.coefficient_unique_patterns += len(
                            unique_patterns
                        )
                    values = self._gather_tensor_requests(
                        permutations,
                        tensor,
                        unique_patterns,
                        inverse
                    )
                    destinations = np.asarray(destinations, dtype=int)
                    active = good[destinations, expansion_index, :]
                    block = products[destinations, expansion_index, :]
                    np.multiply(block, values, out=block, where=active)
                    products[destinations, expansion_index, :] = block
                    good[destinations, expansion_index, :] = np.logical_and(
                        active,
                        np.abs(block) >= zero_cutoff
                    )

        return products, np.any(good, axis=1)

    def gather(self, permutations, subset, num_fixed,
               counts_cache, factorials, zero_cutoff, stats):
        start = time.perf_counter()
        active_products, remapped_products = self._active_products(
            subset,
            num_fixed,
            counts_cache,
            factorials
        )
        if len(active_products) == 0:
            return (
                np.empty(0, dtype=int),
                np.empty(0, dtype=int),
                np.empty((0, self.expansion_count), dtype=float)
            )

        permutation_count = len(permutations)
        # products + active mask + one gathered stage, with a conservative
        # allowance for NumPy's temporary advanced-index arrays.
        bytes_per_product = max(
            1,
            permutation_count * (17 * self.expansion_count + 24)
        )
        chunk_size = max(1, self.max_workspace_bytes // bytes_per_product)
        max_rows = len(active_products) * permutation_count
        eval_permutations = np.empty(max_rows, dtype=int)
        eval_coefficients = np.empty(max_rows, dtype=int)
        prefactors = np.empty(
            (max_rows, self.expansion_count),
            dtype=float
        )
        output_position = 0

        for chunk_start in range(0, len(active_products), chunk_size):
            chunk_stop = min(chunk_start + chunk_size, len(active_products))
            chunk_products = active_products[chunk_start:chunk_stop]
            chunk_remapped = remapped_products[chunk_start:chunk_stop]
            stats.coefficient_product_chunks += 1
            products, good = self._evaluate_product_chunk(
                permutations,
                chunk_products,
                chunk_remapped,
                subset,
                num_fixed,
                zero_cutoff,
                stats
            )
            for local_index, global_index in enumerate(chunk_products):
                active_permutations = np.flatnonzero(good[local_index])
                next_position = output_position + len(active_permutations)
                eval_permutations[output_position:next_position] = active_permutations
                eval_coefficients[output_position:next_position] = global_index
                prefactors[output_position:next_position] = products[
                    local_index, :, active_permutations
                ]
                output_position = next_position

        stats.coefficient_rows += output_position
        stats.coefficient_seconds += time.perf_counter() - start
        return (
            eval_permutations[:output_position],
            eval_coefficients[:output_position],
            prefactors[:output_position]
        )


class IndexedBlockEvaluator:
    """Serial concrete indexed-block executor.

    The current producer intentionally calls the established coefficient gather
    and state-permutation helpers.  The important boundary is that their outputs
    are converted to compact work items and joined before numerical evaluation.
    This makes coefficient gathering independently replaceable in the next
    implementation stage.
    """

    default_chunk_size = 64
    default_plan_cache_size = 128
    default_degeneracy_selection_plan_cache_size = 2048
    default_degeneracy_group_plan_cache_size = 4096
    use_compiled_permutation_routing = True
    use_shared_degeneracy_pool = True
    use_compiled_degeneracy_selection = True
    _plan_cache = IndexedExpressionPlanCache(default_plan_cache_size)

    def __init__(self, evaluator_type, expression, chunk_size=None):
        self.evaluator_type = evaluator_type
        self.expression = expression
        self.chunk_size = (
            self.default_chunk_size if chunk_size is None else max(1, int(chunk_size))
        )
        self.stats = IndexedEvaluationStats()
        self._coefficient_plans = {}
        self._degeneracy_selection_plans = collections.OrderedDict()
        self._degeneracy_group_plans = collections.OrderedDict()
        self._started = time.perf_counter()
        self.plan = self.get_plan(expression)

    @classmethod
    def clear_plan_cache(cls):
        cls._plan_cache.clear()

    def get_plan(self, expression):
        plan = self._plan_cache.get(expression)
        if plan is not None:
            self.stats.plan_cache_hits += 1
            return plan
        self.stats.plan_cache_misses += 1
        plan = self.compile_plan(expression)
        self.stats.plans_compiled += 1
        self.stats.plan_terms += plan.term_count
        self._plan_cache.store(expression, plan)
        return plan

    @staticmethod
    def _compile_subexpression(expression, analytic):
        if isinstance(expression, analytic.PTEnergyChangeProductSum):
            grouped = collections.OrderedDict()
            for ordinal, (energy_changes, polynomial) in enumerate(expression.terms.items()):
                contraction_shape = np.shape(energy_changes)
                term = IndexedEnergyTerm(
                    ordinal,
                    energy_changes,
                    polynomial,
                    contraction_shape
                )
                grouped.setdefault(contraction_shape, []).append(term)
            groups = tuple(
                (shape, tuple(terms)) for shape, terms in grouped.items()
            )
            return IndexedSubexpressionPlan(
                expression,
                False,
                expression.prefactor,
                groups
            ), sum(len(terms) for _, terms in groups)
        if isinstance(expression, (analytic.ProductPTPolynomial, analytic.ProductPTPolynomialSum)):
            return IndexedSubexpressionPlan(
                expression,
                True,
                1,
                ()
            ), 1
        raise TypeError("unsupported indexed subexpression {}".format(type(expression).__name__))

    @classmethod
    def compile_plan(cls, expression):
        # Lazy import avoids an Analytic -> IndexedEvaluator -> Analytic cycle.
        from . import Analytic as analytic

        subexpressions = {}
        term_count = 0
        for key, subexpression in expression.terms.items():
            subplan, subterms = cls._compile_subexpression(subexpression, analytic)
            subexpressions[key] = subplan
            term_count += subterms
        return IndexedExpressionPlan(
            expression,
            types.MappingProxyType(subexpressions),
            term_count
        )

    def finish(self):
        self.stats.elapsed_seconds = time.perf_counter() - self._started
        return self.stats.as_dict()

    def get_coefficient_plan(self, coefficient_specs, coefficient_tensors):
        key = (id(coefficient_specs), id(coefficient_tensors))
        entry = self._coefficient_plans.get(key)
        if (
            entry is not None
            and entry[0] is coefficient_specs
            and entry[1] is coefficient_tensors
        ):
            self.stats.coefficient_plan_hits += 1
            return entry[2]
        self.stats.coefficient_plan_misses += 1
        expansion_count = (
            0 if len(coefficient_tensors) == 0 else
            len(coefficient_tensors[0])
        )
        plan = IndexedCoefficientGatherPlan(
            coefficient_specs,
            coefficient_tensors,
            expansion_count
        )
        self.stats.coefficient_plans_compiled += 1
        self._coefficient_plans[key] = (
            coefficient_specs,
            coefficient_tensors,
            plan
        )
        return plan

    @staticmethod
    def _get_degeneracy_change_sides(expression_key, degenerate_changes):
        if degenerate_changes is None:
            return None, None
        return tuple(
            side.get(expression_key)
            for side in degenerate_changes
        )

    @staticmethod
    def _compile_degeneracy_selection_plan(
            subexpression, left_changes, right_changes,
            use_left, use_right, join):
        selection_ids = {}
        selections = []
        shape_terms = []
        term_count = 0
        for contraction_shape, terms in subexpression.groups:
            routed_terms = []
            for term in terms:
                left_predicate = (
                    None if left_changes is None else
                    left_changes.get(term.energy_changes)
                )
                right_predicate = (
                    None if right_changes is None else
                    right_changes.get(term.energy_changes)
                )
                identity = (id(left_predicate), id(right_predicate))
                entry = selection_ids.get(identity)
                if (
                        entry is None
                        or entry[0] is not left_predicate
                        or entry[1] is not right_predicate
                ):
                    selection_id = len(selections)
                    selections.append((left_predicate, right_predicate))
                    selection_ids[identity] = (
                        left_predicate, right_predicate, selection_id
                    )
                else:
                    selection_id = entry[2]
                routed_terms.append((selection_id, term))
                term_count += 1
            shape_terms.append((contraction_shape, tuple(routed_terms)))
        return IndexedDegeneracySelectionPlan(
            subexpression,
            left_changes,
            right_changes,
            bool(use_left),
            bool(use_right),
            join,
            tuple(selections),
            tuple(shape_terms),
            term_count
        )

    def _get_degeneracy_selection_plan(
            self, subexpression, expression_key, degenerate_changes,
            use_left, use_right, join):
        key = (
            id(subexpression), id(degenerate_changes),
            bool(use_left), bool(use_right), id(join)
        )
        entry = self._degeneracy_selection_plans.get(key)
        if (
                entry is not None
                and entry[0] is subexpression
                and entry[1] is degenerate_changes
                and entry[2] is join
        ):
            self.stats.degeneracy_selection_plan_hits += 1
            self._degeneracy_selection_plans.move_to_end(key)
            return entry[3]

        self.stats.degeneracy_selection_plan_misses += 1
        left_changes, right_changes = self._get_degeneracy_change_sides(
            expression_key, degenerate_changes
        )
        plan = self._compile_degeneracy_selection_plan(
            subexpression,
            left_changes,
            right_changes,
            use_left,
            use_right,
            join
        )
        self.stats.degeneracy_selection_descriptors += len(plan.selections)
        self._degeneracy_selection_plans[key] = (
            subexpression, degenerate_changes, join, plan
        )
        self._degeneracy_selection_plans.move_to_end(key)
        while (
                len(self._degeneracy_selection_plans)
                > self.default_degeneracy_selection_plan_cache_size
        ):
            self._degeneracy_selection_plans.popitem(last=False)
            self.stats.degeneracy_selection_plan_evictions += 1
        self.stats.max_degeneracy_selection_plan_cache_items = max(
            self.stats.max_degeneracy_selection_plan_cache_items,
            len(self._degeneracy_selection_plans)
        )
        return plan

    @staticmethod
    def _normalize_degeneracy_selection(joined, permutation_count):
        if np.ndim(joined) == 0:
            return None if bool(joined) else _NO_SELECTION
        active = np.flatnonzero(joined)
        if len(active) == 0:
            return _NO_SELECTION
        if len(active) == permutation_count:
            return None
        return active

    @staticmethod
    def _canonicalize_degeneracy_selections(selections):
        canonical = []
        selection_classes = []
        for selected in selections:
            if selected is _NO_SELECTION:
                selection_classes.append(-2)
            elif selected is None:
                selection_classes.append(-1)
            else:
                selection_class = None
                for candidate_id, candidate in enumerate(canonical):
                    if (
                            len(candidate) == len(selected)
                            and np.array_equal(candidate, selected)
                    ):
                        selection_class = candidate_id
                        selected = candidate
                        break
                if selection_class is None:
                    selection_class = len(canonical)
                    canonical.append(selected)
                selection_classes.append(selection_class)
        return tuple(selection_classes)

    @staticmethod
    def _compile_degeneracy_group_plan(selection_plan, selection_classes):
        groups = collections.OrderedDict()
        for contraction_shape, routed_terms in selection_plan.shape_terms:
            for selection_id, term in routed_terms:
                selection_class = selection_classes[selection_id]
                if selection_class == -2:
                    continue
                group_key = (contraction_shape, selection_class)
                group = groups.get(group_key)
                if group is None:
                    group = [selection_id, []]
                    groups[group_key] = group
                group[1].append(term)
        return tuple(
            (contraction_shape, selection_id, tuple(terms))
            for (contraction_shape, _), (selection_id, terms)
            in groups.items()
        )

    def _get_degeneracy_group_plan(self, selection_plan, selection_classes):
        key = (id(selection_plan), selection_classes)
        entry = self._degeneracy_group_plans.get(key)
        if entry is not None and entry[0] is selection_plan:
            self.stats.degeneracy_group_plan_hits += 1
            self._degeneracy_group_plans.move_to_end(key)
            return entry[1]

        self.stats.degeneracy_group_plan_misses += 1
        group_plan = self._compile_degeneracy_group_plan(
            selection_plan, selection_classes
        )
        self._degeneracy_group_plans[key] = (selection_plan, group_plan)
        self._degeneracy_group_plans.move_to_end(key)
        while (
                len(self._degeneracy_group_plans)
                > self.default_degeneracy_group_plan_cache_size
        ):
            self._degeneracy_group_plans.popitem(last=False)
            self.stats.degeneracy_group_plan_evictions += 1
        self.stats.max_degeneracy_group_plan_cache_items = max(
            self.stats.max_degeneracy_group_plan_cache_items,
            len(self._degeneracy_group_plans)
        )
        return group_plan

    def _select_energy_groups_compiled(
            self, blocks, block_index, degenerate_changes,
            only_degenerate_terms, degeneracy_context=None):
        subexpression = blocks.plans[block_index]
        use_left, use_right, join = only_degenerate_terms
        selection_plan = self._get_degeneracy_selection_plan(
            subexpression,
            blocks.expression_keys[block_index],
            degenerate_changes,
            use_left,
            use_right,
            join
        )
        self.stats.degeneracy_selection_misses += len(selection_plan.selections)
        self.stats.degeneracy_selection_hits += (
            selection_plan.term_count - len(selection_plan.selections)
        )

        permutation_subsets = blocks.permutation_subsets[block_index]
        global_indices = blocks.global_permutation_indices[block_index]
        block_selection_key = None
        resolved = selection_classes = None
        if degeneracy_context is not None:
            (
                block_selection_key,
                resolved,
                selection_classes
            ) = degeneracy_context.lookup_block_selection(
                selection_plan, global_indices
            )
        if resolved is None:
            predicate_masks = {}
            resolved = []
            for left_predicate, right_predicate in selection_plan.selections:
                if degeneracy_context is None:
                    left_matches = self._get_degeneracy_mask(
                        left_predicate,
                        permutation_subsets,
                        use_left,
                        predicate_masks
                    )
                    right_matches = self._get_degeneracy_mask(
                        right_predicate,
                        permutation_subsets,
                        use_right,
                        predicate_masks
                    )
                    joined = join(left_matches, right_matches)
                else:
                    joined = degeneracy_context.selection_mask(
                        left_predicate,
                        right_predicate,
                        use_left,
                        use_right,
                        join
                    )
                    if np.ndim(joined) > 0:
                        joined = joined[global_indices]
                resolved.append(self._normalize_degeneracy_selection(
                    joined, len(permutation_subsets)
                ))

            resolved = tuple(resolved)
            selection_classes = self._canonicalize_degeneracy_selections(
                resolved
            )
            if degeneracy_context is not None:
                degeneracy_context.store_block_selection(
                    block_selection_key,
                    selection_plan,
                    resolved,
                    selection_classes
                )
        group_plan = self._get_degeneracy_group_plan(
            selection_plan, selection_classes
        )
        return tuple(
            (contraction_shape, resolved[selection_id], terms)
            for contraction_shape, selection_id, terms in group_plan
        )

    def _select_energy_groups(self, blocks, block_index, degenerate_changes,
                              only_degenerate_terms,
                              degeneracy_context=None):
        if self.use_compiled_degeneracy_selection:
            return self._select_energy_groups_compiled(
                blocks,
                block_index,
                degenerate_changes,
                only_degenerate_terms,
                degeneracy_context=degeneracy_context
            )
        return self._select_energy_groups_legacy(
            blocks,
            block_index,
            degenerate_changes,
            only_degenerate_terms,
            degeneracy_context=degeneracy_context
        )

    def _select_energy_groups_legacy(self, blocks, block_index,
                                     degenerate_changes,
                                     only_degenerate_terms,
                                     degeneracy_context=None):
        """Attach concrete degeneracy masks to the static energy groups."""
        plan = blocks.plans[block_index]
        use_left, use_right, join = only_degenerate_terms
        if degenerate_changes is None:
            left_changes = right_changes = None
        else:
            left_changes, right_changes = [
                side.get(blocks.expression_keys[block_index], {})
                for side in degenerate_changes
            ]

        groups = collections.OrderedDict()
        # A permutation-subset object is private to this concrete block, so masks
        # cannot be reused by later work items through identity.  Keep only the
        # useful predicate reuse within this finite term scan and release all
        # masks with the block instead of churning a long-lived LRU.
        permutation_subsets = blocks.permutation_subsets[block_index]
        predicate_masks = {}
        selection_cache = {}
        canonical_selections = []
        for contraction_shape, terms in plan.groups:
            for term in terms:
                left_predicate = (
                    None if left_changes is None else
                    left_changes.get(term.energy_changes, None)
                )
                right_predicate = (
                    None if right_changes is None else
                    right_changes.get(term.energy_changes, None)
                )
                selection_key = (
                    id(left_predicate), id(right_predicate),
                    use_left, use_right
                )
                entry = selection_cache.get(selection_key)
                if (
                    entry is not None
                    and entry[0] is left_predicate
                    and entry[1] is right_predicate
                ):
                    self.stats.degeneracy_selection_hits += 1
                    selected = entry[2]
                    selection_id = entry[3]
                else:
                    self.stats.degeneracy_selection_misses += 1
                    if degeneracy_context is None:
                        left_matches = self._get_degeneracy_mask(
                            left_predicate,
                            permutation_subsets,
                            use_left,
                            predicate_masks
                        )
                        right_matches = self._get_degeneracy_mask(
                            right_predicate,
                            permutation_subsets,
                            use_right,
                            predicate_masks
                        )
                        joined = join(left_matches, right_matches)
                    else:
                        joined = degeneracy_context.selection_mask(
                            left_predicate,
                            right_predicate,
                            use_left,
                            use_right,
                            join
                        )
                        if np.ndim(joined) > 0:
                            joined = joined[
                                blocks.global_permutation_indices[block_index]
                            ]
                    if np.ndim(joined) == 0:
                        selected = None if bool(joined) else _NO_SELECTION
                    else:
                        active = np.flatnonzero(joined)
                        if len(active) == 0:
                            selected = _NO_SELECTION
                        elif len(active) == len(permutation_subsets):
                            selected = None
                        else:
                            selected = active
                    if selected is None or selected is _NO_SELECTION:
                        selection_id = selected
                    else:
                        selection_id = None
                        for candidate_id, candidate in enumerate(
                                canonical_selections):
                            if (
                                    len(candidate) == len(selected)
                                    and np.array_equal(candidate, selected)
                            ):
                                selection_id = candidate_id
                                selected = candidate
                                break
                        if selection_id is None:
                            selection_id = len(canonical_selections)
                            canonical_selections.append(selected)
                    selection_cache[selection_key] = (
                        left_predicate, right_predicate, selected, selection_id
                    )
                if selected is _NO_SELECTION:
                    continue
                group_key = (contraction_shape, selection_id)
                group = groups.get(group_key)
                if group is None:
                    group = [selected, []]
                    groups[group_key] = group
                group[1].append(term)
        return tuple(
            (contraction_shape, selected, terms)
            for (contraction_shape, _), (selected, terms) in groups.items()
        )

    def _get_degeneracy_mask(self, predicate, permutation_subsets,
                             only_degenerate, predicate_masks):
        if predicate is None:
            return not only_degenerate

        predicate_key = id(predicate)
        entry = predicate_masks.get(predicate_key)
        if entry is None or entry[0] is not predicate:
            self.stats.degeneracy_mask_misses += 1
            tick = time.perf_counter()
            evaluate = getattr(predicate, 'evaluate', None)
            self.stats.degeneracy_states_tested += len(permutation_subsets)
            array_cutoff = getattr(
                predicate,
                'scalar_state_cutoff',
                -1
            )
            if evaluate is None or len(permutation_subsets) <= array_cutoff:
                self.stats.degeneracy_scalar_evaluations += 1
                mask = np.fromiter(
                    (
                        bool(predicate(permutation))
                        for permutation in permutation_subsets
                    ),
                    dtype=bool,
                    count=len(permutation_subsets)
                )
            else:
                self.stats.degeneracy_array_evaluations += 1
                self.stats.degeneracy_patterns_tested += len(
                    predicate.patterns
                )
                mask = np.asarray(
                    evaluate(permutation_subsets),
                    dtype=bool
                )
            self.stats.degeneracy_seconds += time.perf_counter() - tick
            predicate_masks[predicate_key] = (predicate, mask)
        else:
            self.stats.degeneracy_mask_hits += 1
            mask = entry[1]
        return mask if only_degenerate else np.logical_not(mask)

    @staticmethod
    def _contract_energy_rows(energy_changes, frequency_blocks):
        """Contract row x permutation x mode frequency blocks in one call."""
        linear_factors = np.einsum(
            'rpm,rfm->rpf',
            frequency_blocks,
            energy_changes[..., 1:],
            optimize=False
        )
        return np.prod(linear_factors, axis=-1)

    def _materialize_polynomial_block(self, polynomial_factors,
                                      state_count, permutation_count):
        """Fill one bounded dense block using NumPy broadcasting assignment.

        This deliberately avoids constructing one ``broadcast_to`` view per
        row and then stacking those views.  The returned array is the only
        new numeric data buffer: scalar and array factors are assigned
        directly into their rows and are never cached in expanded form.
        """
        dtype_specs = []
        for polynomial in polynomial_factors:
            dtype = getattr(polynomial, 'dtype', None)
            if dtype is None:
                dtype = (
                    type(polynomial) if np.isscalar(polynomial) else
                    np.asanyarray(polynomial).dtype
                )
            dtype_specs.append(dtype)
        dtype = np.result_type(*dtype_specs)
        polynomial_block = np.empty(
            (
                len(polynomial_factors),
                state_count,
                permutation_count
            ),
            dtype=dtype
        )
        for row, polynomial in enumerate(polynomial_factors):
            try:
                polynomial_block[row] = polynomial
            except ValueError as exc:
                raise ValueError(
                    "polynomial factor {} with shape {} cannot broadcast to "
                    "({}, {})".format(
                        row,
                        np.shape(polynomial),
                        state_count,
                        permutation_count
                    )
                ) from exc
        self.stats.max_polynomial_scratch_bytes = max(
            self.stats.max_polynomial_scratch_bytes,
            polynomial_block.nbytes
        )
        return polynomial_block

    def _evaluate_work_items(self, blocks, output, change, baseline_shift,
                             pows, polynomial_cache, degenerate_changes,
                             only_degenerate_terms, zero_cutoff, logger,
                             evaluation_context=None,
                             degeneracy_context=None):
        evaluator = self.evaluator_type
        row_groups = collections.OrderedDict()

        for block_index in range(len(blocks)):
            plan = blocks.plans[block_index]
            if plan.direct_polynomial:
                use_left, use_right, _ = only_degenerate_terms
                if use_left or use_right:
                    raise ValueError(
                        "degenerate terms requested on {}".format(
                            plan.expression.format_expr()
                        )
                    )
                polynomial = evaluator._eval_poly(
                    polynomial_cache,
                    blocks.tuple_states[block_index],
                    blocks.permutation_substates[block_index],
                    pows,
                    plan.expression,
                    change,
                    baseline_shift,
                    False,
                    logger,
                    block_key=(
                        'full', blocks.state_block_identities[block_index]
                    )
                )
                blocks.contributions[block_index] += polynomial
                continue

            selected_groups = self._select_energy_groups(
                blocks,
                block_index,
                degenerate_changes,
                only_degenerate_terms,
                degeneracy_context=degeneracy_context
            )
            for energy_shape, selected, terms in selected_groups:
                permutation_count = (
                    len(blocks.permutation_subsets[block_index])
                    if selected is None else len(selected)
                )
                group_key = (
                    energy_shape,
                    permutation_count,
                    blocks.contributions[block_index].shape[0]
                )
                segments = row_groups.setdefault(
                    group_key,
                    IndexedContractionSegmentTable()
                )
                segments.append(block_index, selected, terms)

        self.stats.contraction_groups += len(row_groups)
        for segments in row_groups.values():
            self.stats.contraction_segments += len(segments.block_indices)
            self.stats.max_pending_segments = max(
                self.stats.max_pending_segments,
                len(segments.block_indices)
            )
            for block_indices, selected_rows, terms, segment_ids in segments.iter_blocks(
                    self.chunk_size, return_segment_ids=True):
                block_size = len(terms)
                self.stats.contraction_blocks += 1
                self.stats.contraction_rows += block_size
                self.stats.max_rows_per_block = max(
                    self.stats.max_rows_per_block,
                    block_size
                )

                selections = []
                for block_index, selected in zip(block_indices, selected_rows):
                    block_index = int(block_index)
                    block_identity = blocks.state_block_identities[block_index]
                    if selected is None:
                        selections.append((
                            blocks.permutation_frequencies[block_index],
                            blocks.permutation_substates[block_index],
                            blocks.tuple_states[block_index],
                            ('full', block_identity)
                        ))
                    else:
                        tuple_states = blocks.tuple_states[block_index]
                        selections.append((
                            blocks.permutation_frequencies[block_index][selected,],
                            blocks.permutation_substates[block_index][selected,],
                            [tuple_states[index] for index in selected],
                            ('take', block_identity, tuple(selected))
                        ))

                tick = time.perf_counter()
                energy_changes = np.asarray([
                    term.energy_changes for term in terms
                ])
                frequency_blocks = np.stack([
                    selection[0] for selection in selections
                ])
                energy_factors = self._contract_energy_rows(
                    energy_changes,
                    frequency_blocks
                )
                self.stats.energy_seconds += time.perf_counter() - tick

                tick = time.perf_counter()
                polynomial_factors = [
                    evaluator._eval_poly(
                        polynomial_cache,
                        selection[2],
                        selection[1],
                        pows,
                        term.polynomial,
                        change,
                        baseline_shift,
                        False,
                        logger,
                        block_key=selection[3]
                    )
                    for term, selection in zip(terms, selections)
                ]
                self.stats.polynomial_evaluations += len(polynomial_factors)
                self.stats.polynomial_seconds += time.perf_counter() - tick

                tick = time.perf_counter()
                state_count = blocks.contributions[int(block_indices[0])].shape[0]
                permutation_count = energy_factors.shape[1]
                polynomial_block = self._materialize_polynomial_block(
                    polynomial_factors,
                    state_count,
                    permutation_count
                )
                operator_repr = [None]

                def row_diagnostics(row_index):
                    if operator_repr[0] is None:
                        operator_repr[0] = repr(
                            self.plan.expression
                            if evaluation_context is None else evaluation_context
                        )
                    block_index = int(block_indices[row_index])
                    selection = selections[row_index]
                    term = terms[row_index]
                    return {
                        'evaluation_context': {
                            'operator': operator_repr[0],
                            'coefficient_key': blocks.expression_keys[block_index],
                            'energy_term': term.ordinal,
                            'change': None if change is None else tuple(change),
                            'baseline_shift': (
                                None if baseline_shift is None else tuple(baseline_shift)
                            )
                        },
                        'energy_changes': term.energy_changes,
                        'tuple_states': selection[2],
                        'permutation_substates': selection[1]
                    }

                scaled_rows = evaluator._divide_polynomial_block_by_energy(
                    polynomial_block,
                    energy_factors,
                    zero_cutoff,
                    row_diagnostics=row_diagnostics
                )
                segment_starts = np.concatenate((
                    np.array([0], dtype=np.intp),
                    np.flatnonzero(segment_ids[1:] != segment_ids[:-1]) + 1
                ))
                reduced_rows = np.add.reduceat(
                    scaled_rows, segment_starts, axis=0
                )
                self.stats.reduction_scatters += len(segment_starts)
                for reduced_index, row_index in enumerate(segment_starts):
                    block_index = int(block_indices[row_index])
                    selected = selected_rows[row_index]
                    scaled = reduced_rows[reduced_index]
                    if selected is None:
                        blocks.contributions[block_index] += scaled
                    else:
                        blocks.contributions[block_index][:, selected] += scaled
                self.stats.reduction_seconds += time.perf_counter() - tick

        for block_index in range(len(blocks)):
            contribution = blocks.contributions[block_index]
            plan = blocks.plans[block_index]
            if not plan.direct_polynomial:
                contribution *= plan.prefactor
            value = (
                blocks.prefactors[block_index][:, np.newaxis, :]
                * contribution[np.newaxis]
            )
            output[blocks.state_indices[block_index]][
                :, blocks.mask_positions[block_index]
            ] += value[:, 0, :]

    def _append_routed_work_item(
            self, work_items, state_index, mask_positions, value_positions,
            global_permutation_indices,
            expression_key, state_permutations, prefactors, frequencies,
            num_fixed, subset, full_set, state_permutation_cache,
            degenerate_changes):
        evaluator = self.evaluator_type
        state, permutations = state_permutations[state_index]
        state_block = state[np.newaxis]
        block_prefactors = prefactors[value_positions].T
        permutation_subsets = permutations[
            mask_positions,
        ][:, full_set]
        (
            permutation_substates,
            tuple_states,
            permutation_frequencies,
            state_block_identity
        ) = evaluator._get_state_perms(
            state_index,
            state_block,
            frequencies,
            num_fixed,
            subset,
            permutation_subsets,
            state_permutation_cache,
            mask_positions,
            full_set
        )
        if degenerate_changes is not None:
            permutation_subsets = [
                tuple(permutation)
                for permutation in permutation_subsets
            ]
        work_items.append(
            state_index,
            mask_positions,
            global_permutation_indices,
            block_prefactors,
            expression_key,
            self.plan.subexpressions[expression_key],
            permutation_subsets,
            permutation_substates,
            tuple_states,
            permutation_frequencies,
            state_block_identity,
            np.zeros((len(state_block), len(permutation_subsets)))
        )

    def _build_work_items_legacy_routing(
            self, state_permutations, all_permutations, permutation_map,
            split_spec, eval_permutations, eval_coefficients, prefactors,
            coefficient_specs, frequencies, num_fixed, subset, full_set,
            state_permutation_cache, degenerate_changes):
        """Retained mask/split implementation for controlled comparisons."""
        all_mask = np.full(len(all_permutations), False)
        subgroup_map = np.zeros(len(all_permutations), dtype=int)
        from . import Analytic as analytic
        (coefficient_groups, value_groups), _ = analytic.nput.group_by(
            np.arange(len(eval_coefficients)),
            eval_coefficients
        )
        split_permutation_indices = np.split(permutation_map, split_spec)
        work_items = ConcreteBlockTable()
        for coefficient_group, value_group in zip(
                coefficient_groups, value_groups):
            active_global_permutations = eval_permutations[value_group,]
            expression_key = coefficient_specs[coefficient_group]
            subgroup_map[active_global_permutations] = value_group
            all_mask[active_global_permutations] = True
            split_active = np.split(all_mask[permutation_map], split_spec)
            split_value_positions = np.split(
                subgroup_map[permutation_map], split_spec
            )

            for state_index, (
                    active, value_positions, _
            ) in enumerate(zip(
                split_active,
                split_value_positions,
                split_permutation_indices
            )):
                mask_positions = np.flatnonzero(active)
                if len(mask_positions) == 0:
                    continue
                self._append_routed_work_item(
                    work_items,
                    state_index,
                    mask_positions,
                    value_positions[mask_positions,],
                    eval_permutations[value_positions[mask_positions,]],
                    expression_key,
                    state_permutations,
                    prefactors,
                    frequencies,
                    num_fixed,
                    subset,
                    full_set,
                    state_permutation_cache,
                    degenerate_changes
                )
            all_mask[active_global_permutations] = False
        return work_items

    def _build_work_items_compiled_routing(
            self, state_permutations, all_permutations, permutation_map,
            eval_permutations, eval_coefficients, prefactors,
            coefficient_specs, frequencies, num_fixed, subset, full_set,
            state_permutation_cache, degenerate_changes):
        """Join coefficient rows directly to compiled state destinations."""
        routing_plan = IndexedPermutationRoutingPlan.compile(
            permutation_map,
            [len(permutations) for _, permutations in state_permutations],
            len(all_permutations)
        )
        self.stats.routing_plans_compiled += 1
        routes = routing_plan.route(eval_permutations, eval_coefficients)
        self.stats.routing_matches += len(routes.value_positions)
        work_items = ConcreteBlockTable()
        for (
                coefficient_group, state_index,
                mask_positions, value_positions
        ) in routes:
            expression_key = coefficient_specs[coefficient_group]
            self._append_routed_work_item(
                work_items,
                state_index,
                mask_positions,
                value_positions,
                eval_permutations[value_positions],
                expression_key,
                state_permutations,
                prefactors,
                frequencies,
                num_fixed,
                subset,
                full_set,
                state_permutation_cache,
                degenerate_changes
            )
        return work_items

    def eval_perm(self, expression, change, baseline_shift,
                  subset, state_permutations, all_permutations, permutation_map,
                  frequencies, coefficient_specs, coefficient_tensors,
                  num_fixed, degenerate_changes, only_degenerate_terms,
                  zero_cutoff, counts_cache, polynomial_cache,
                  state_permutation_cache, energy_cache,
                  factorials, pows, split_spec,
                  verbose, logger, log_scaled,
                  evaluation_context=None):
        """Evaluate one free-index subset through the joined block executor."""
        if verbose:
            raise NotImplementedError(
                "verbose indexed evaluation is implemented in the logging stage"
            )

        self.stats.calls += 1
        full_set = tuple(range(num_fixed)) + subset
        coefficient_plan = self.get_coefficient_plan(
            coefficient_specs,
            coefficient_tensors
        )
        eval_permutations, eval_coefficients, prefactors = coefficient_plan.gather(
            all_permutations,
            subset,
            num_fixed,
            counts_cache,
            factorials,
            zero_cutoff,
            self.stats
        )
        output = [
            np.zeros((prefactors.shape[1], len(permutations)))
            for _, permutations in state_permutations
        ]
        if prefactors.shape[0] == 0:
            return output

        tick = time.perf_counter()
        if self.use_compiled_permutation_routing:
            work_items = self._build_work_items_compiled_routing(
                state_permutations,
                all_permutations,
                permutation_map,
                eval_permutations,
                eval_coefficients,
                prefactors,
                coefficient_specs,
                frequencies,
                num_fixed,
                subset,
                full_set,
                state_permutation_cache,
                degenerate_changes
            )
        else:
            work_items = self._build_work_items_legacy_routing(
                state_permutations,
                all_permutations,
                permutation_map,
                split_spec,
                eval_permutations,
                eval_coefficients,
                prefactors,
                coefficient_specs,
                frequencies,
                num_fixed,
                subset,
                full_set,
                state_permutation_cache,
                degenerate_changes
            )
        self.stats.routing_seconds += time.perf_counter() - tick

        self.stats.work_items += len(work_items)
        degeneracy_context = (
            None if (
                degenerate_changes is None
                or not self.use_shared_degeneracy_pool
            ) else
            IndexedDegeneracyEvaluationContext(
                np.asarray(all_permutations)[:, full_set],
                self.stats
            )
        )
        self._evaluate_work_items(
            work_items,
            output,
            change,
            baseline_shift,
            pows,
            polynomial_cache,
            degenerate_changes,
            only_degenerate_terms,
            zero_cutoff,
            logger,
            evaluation_context=evaluation_context,
            degeneracy_context=degeneracy_context
        )
        return output
