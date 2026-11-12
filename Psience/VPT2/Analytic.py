from pathlib import Path
"""
Provides a symbolic approach to vibrational perturbation theory based on a Harmonic description
"""

import abc, itertools, collections, enum, math, pickle, weakref, sys, warnings
import contextlib
import copy
import contextvars
import functools
import hashlib
import inspect

import numpy as np, scipy.signal, time

from McUtils.Zachary import DensePolynomial, TensorCoefficientPoly
import McUtils.Numputils as nput
import McUtils.Devutils as dev
from McUtils.Devutils import CascadingPropertyManager, cascading_property, derived_property
import McUtils.Combinatorics as mcomb
from McUtils.Combinatorics import SymmetricGroupGenerator, IntegerPartitioner, UniquePartitions, UniquePermutations
from McUtils.Scaffolding import Logger, Checkpointer, MaxSizeCache
from McUtils.Parallelizers import (
    Parallelizer, MultiprocessingParallelizer, SharedMemoryArrayTree
)
from ..BasisReps import (
    BasisStateSpace, HarmonicOscillatorProductBasis,
    HarmonicOscillatorMatrixGenerator, HarmonicOscillatorRaisingLoweringPolyTerms
)

from .Corrections import BasicAPTCorrections
from .Common import Diagnostics, BoundedCache
# imports from merged modules
import json
from McUtils.Zachary import ElementaryFunction

__all__ = [
    'OperatorCondition',
    'SameStateCondition',
    'PolynomialInvariantCondition',
    'PolynomialRelationshipCondition',
    'EnergyWindowCondition',
    'IntegerLatticeCondition',
    'ConditionalTransformationStack',
    'BlockHamiltonianCorrection',
    'PerturbationTheoryEvaluator',
    'AnalyticPerturbationTheorySolver',
    'PolyAtom',
    'PolyAxis',
    'PolyTerm',
    'PolyPath',
    'PTTensorCoeffProductDAG',
    'PTTensorCoeffProductDAGEvaluationPlan',
    'DegeneracyTestPlan',
    'CompactDegeneracyTestPlan',
    'DegeneracyChangeIndex',
    'DegeneracyIdentificationContext',
    'SmallEnergyDenominatorWarning',
    'SurvivingEnergyDenominatorWarning',
    # 'AnalyticPerturbationTheoryDriver',
    # 'AnalyticPTCorrectionGenerator',
    # 'RaisingLoweringClasses'
]

_DEBUG_PRINT = False # statements will be moved into a proper logger later
_DEBUG_AUDIT = False
_PERMUTE_CHANGES = False # for debug purposes
_PERMUTE_FINALS = False
_TAKE_UNIQUE_CHANGES = False

# class DefaultObject:
#     def __repr__(self): return "default"
class DefaultValues(enum.Enum):
    DEFAULT = 'default'
default = DefaultValues.DEFAULT


_QUANTUM_INDEX_LABELS = "ijklmnopqrstuvwxyzabcdefgh"


def format_quantum_change(change):
    """Render a change vector as a polynomial in local mode indices.

    This is a display-only representation: the ordered integer tuple remains
    the key used by solvers, evaluators, and checkpoints.  In particular, a
    NumPy scalar is converted to an ordinary integer before formatting.
    A future solver-owned LOGGER_INDEX_CACHE can substitute D_n labels here
    without changing the calculation's keys or this algebraic fallback.
    """
    terms = []
    for index, entry in enumerate(change):
        amount = int(entry)
        if amount == 0:
            continue
        label = (
            "q_" + _QUANTUM_INDEX_LABELS[index]
            if index < len(_QUANTUM_INDEX_LABELS) else
            "q_{{{}}}".format(index)
        )
        magnitude = "" if abs(amount) == 1 else str(abs(amount))
        sign = "-" if amount < 0 else ("+" if terms else "")
        terms.append(sign + magnitude + label)
    return "".join(terms) if terms else "0"


def format_quantum_change_list(changes):
    """Render an ordered collection of change vectors for operator labels."""
    if changes is None:
        return "None"
    return "[{}]".format(
        ", ".join(format_quantum_change(change) for change in changes)
    )


class SmallEnergyDenominatorWarning(RuntimeWarning):
    """An evaluated perturbation term contains a denominator below cutoff."""


class SurvivingEnergyDenominatorWarning(SmallEnergyDenominatorWarning):
    """A small denominator retained a numerator after cutoff pruning."""


# A distinct type is required: a legacy change list can itself contain two
# entries, so its length cannot identify side-specific resonance rules.
SideDegeneracyChanges = collections.namedtuple(
    'SideDegeneracyChanges', ['left', 'right']
)


class StateDegeneracyPairs(list):
    """Actual block neighbors, ordered from one reference state to each neighbor."""


class DegeneracyTestPlan:
    """Array-oriented union of mode-pattern degeneracy predicates.

    Negative pattern entries are wildcards.  The object remains callable for
    the legacy evaluator, while indexed evaluation can test a complete
    permutation block with :meth:`evaluate` and avoid Python predicate/tree
    traversal for each state.
    """

    default_compare_workspace_bytes = 4 * 1024 ** 2
    default_scalar_state_cutoff = 32

    def __init__(self, patterns=None, checks=None, fallback_tests=(),
                 compare_workspace_bytes=None):
        if patterns is None:
            patterns = np.empty((0, 0), dtype=np.intp)
        patterns = np.asarray(patterns, dtype=np.intp)
        if patterns.ndim == 1:
            patterns = patterns[np.newaxis]
        if patterns.ndim != 2:
            raise ValueError("degeneracy patterns must be a two-dimensional array")
        if checks is None:
            checks = patterns >= 0
        checks = np.asarray(checks, dtype=bool)
        if checks.shape != patterns.shape:
            raise ValueError("degeneracy checks must match the pattern shape")
        if len(patterns) > 1:
            # Duplicate patterns arise naturally from equivalent energy-change
            # permutations.  Remove them once when the plan is compiled.
            packed = np.concatenate((patterns, checks.astype(np.intp)), axis=1)
            _, unique = np.unique(packed, axis=0, return_index=True)
            unique.sort()
            patterns = patterns[unique]
            checks = checks[unique]
        patterns.setflags(write=False)
        checks.setflags(write=False)
        self.patterns = patterns
        self.checks = checks
        self.fallback_tests = tuple(fallback_tests)
        self.scalar_state_cutoff = max(
            self.default_scalar_state_cutoff,
            2 * len(patterns)
        )
        self.compare_workspace_bytes = (
            self.default_compare_workspace_bytes
            if compare_workspace_bytes is None else
            max(1, int(compare_workspace_bytes))
        )
        # The array membership tables are compiled lazily.  Most current VPT
        # work items contain only one or two permutations.  Indexed evaluation
        # now combines those items through a shared permutation pool, while the
        # legacy evaluator still uses the scalar trie.  Build neither form
        # until its backend actually requests it.
        self.pattern_groups = None
        self.scalar_tree = None
        self.matches_all = bool(np.any(np.logical_not(np.any(checks, axis=1))))

    @property
    def cache_key(self):
        return (
            self.patterns.shape,
            self.patterns.tobytes(),
            self.checks.tobytes(),
            self.fallback_tests
        )

    def evaluate(self, states):
        if isinstance(states, np.ndarray) and states.ndim == 1:
            return self(states)
        if len(states) <= self.scalar_state_cutoff:
            return np.fromiter(
                (self(state) for state in states),
                dtype=bool,
                count=len(states)
            )
        states = np.asarray(states, dtype=np.intp)
        if states.ndim != 2:
            raise ValueError("degeneracy states must be a two-dimensional array")
        if len(self.patterns) > 0 and states.shape[1] != self.patterns.shape[1]:
            raise ValueError(
                "degeneracy state width {} does not match pattern width {}"
                .format(states.shape[1], self.patterns.shape[1])
            )

        result = np.zeros(len(states), dtype=bool)
        if len(self.patterns) > 0 and len(states) > 0:
            for positions, values, keys in self._get_pattern_groups():
                if len(positions) == 0:
                    result[:] = True
                elif len(positions) == 1:
                    result |= np.isin(
                        states[:, positions[0]],
                        values[:, 0]
                    )
                else:
                    projected = np.ascontiguousarray(states[:, positions])
                    state_keys = projected.view(keys.dtype).ravel()
                    result |= np.isin(state_keys, keys)
                if np.all(result):
                    break
        if len(self.fallback_tests) > 0 and not np.all(result):
            missing = np.flatnonzero(np.logical_not(result))
            result[missing] = np.fromiter(
                (
                    any(test(states[index]) for test in self.fallback_tests)
                    for index in missing
                ),
                dtype=bool,
                count=len(missing)
            )
        return result

    def evaluate_array(self, states):
        """Evaluate an explicitly batched state/permutation pool.

        Unlike :meth:`evaluate`, this never falls back to the scalar trie for a
        small batch.  The indexed evaluator uses it after coalescing work-item
        subsets into one shared permutation pool.
        """
        states = np.asarray(states, dtype=np.intp)
        if states.ndim != 2:
            raise ValueError("degeneracy states must be a two-dimensional array")
        if len(self.patterns) > 0 and states.shape[1] != self.patterns.shape[1]:
            raise ValueError(
                "degeneracy state width {} does not match pattern width {}"
                .format(states.shape[1], self.patterns.shape[1])
            )

        result = np.zeros(len(states), dtype=bool)
        if len(self.patterns) > 0 and len(states) > 0:
            for positions, values, keys in self._get_pattern_groups():
                if len(positions) == 0:
                    result[:] = True
                elif len(positions) == 1:
                    result |= np.isin(states[:, positions[0]], values[:, 0])
                else:
                    projected = np.ascontiguousarray(states[:, positions])
                    state_keys = projected.view(keys.dtype).ravel()
                    result |= np.isin(state_keys, keys)
                if np.all(result):
                    break
        if len(self.fallback_tests) > 0 and not np.all(result):
            missing = np.flatnonzero(np.logical_not(result))
            result[missing] = np.fromiter(
                (
                    any(test(states[index]) for test in self.fallback_tests)
                    for index in missing
                ),
                dtype=bool,
                count=len(missing)
            )
        return result

    def _get_scalar_tree(self):
        scalar_tree = self.scalar_tree
        if scalar_tree is not None:
            return scalar_tree
        scalar_tree = {}
        for pattern, check in zip(self.patterns, self.checks):
            positions = tuple(int(index) for index in np.flatnonzero(check))
            if len(positions) == 0:
                continue
            subtree = scalar_tree
            for position_index, position in enumerate(positions):
                match_values = subtree.setdefault(position, {})
                value = int(pattern[position])
                if position_index == len(positions) - 1:
                    match_values[value] = True
                else:
                    next_tree = match_values.get(value)
                    if next_tree is True:
                        break
                    if next_tree is None:
                        next_tree = {}
                        match_values[value] = next_tree
                    subtree = next_tree
        self.scalar_tree = scalar_tree
        return scalar_tree

    def _get_pattern_groups(self):
        pattern_groups = self.pattern_groups
        if pattern_groups is not None:
            return pattern_groups
        grouped = collections.OrderedDict()
        for pattern, check in zip(self.patterns, self.checks):
            positions = tuple(int(index) for index in np.flatnonzero(check))
            grouped.setdefault(positions, []).append(pattern[list(positions)])
        compiled = []
        for positions, values in grouped.items():
            positions = np.asarray(positions, dtype=np.intp)
            values = np.stack(values)
            if len(values) > 1:
                values = np.unique(values, axis=0)
            if len(positions) > 1:
                key_dtype = np.dtype((
                    np.void,
                    values.dtype.itemsize * len(positions)
                ))
                keys = np.ascontiguousarray(values).view(key_dtype).ravel()
                keys.setflags(write=False)
            else:
                keys = None
            positions.setflags(write=False)
            values.setflags(write=False)
            compiled.append((positions, values, keys))
        self.pattern_groups = tuple(compiled)
        return self.pattern_groups

    @classmethod
    def _evaluate_scalar_tree(cls, state, tree):
        for position, match_values in tree.items():
            match = match_values.get(int(state[position]))
            if isinstance(match, dict):
                match = cls._evaluate_scalar_tree(state, match)
            if match:
                return True
        return False

    def __call__(self, state):
        matched = (
            self.matches_all
            or self._evaluate_scalar_tree(state, self._get_scalar_tree())
        )
        if not matched and len(self.fallback_tests) > 0:
            matched = any(test(state) for test in self.fallback_tests)
        return bool(matched)


class CompactDegeneracyTestPlan:
    """Union of linear resonance clauses without expanded pattern rows.

For one energy-change requirement, each compiled inverse mode map projects a
concrete permutation back into the sorted mode tuple held by
``DegeneracyChangeIndex``. Membership in that small source index is exactly
the test that the expanded pattern union performs. In particular, keeping
the actual mode maps preserves the current all-equal-quanta behavior.
    """

    default_scalar_state_cutoff = 32

    def __init__(self, clauses, width, mode_index=None):
        self.clauses = tuple(clauses)
        self.width = int(width)
        self.fallback_tests = ()
        self.pattern_count = sum(clause[5] for clause in self.clauses)
        self.scalar_state_cutoff = max(
            self.default_scalar_state_cutoff,
            2 * sum(len(clause[1]) for clause in self.clauses)
        )
        # A shape-only compatibility view for callers that inspect predicate
        # arity; no Cartesian-product pattern rows are allocated.
        self.patterns = np.empty((0, self.width), dtype=np.intp)
        self.checks = np.empty((0, self.width), dtype=bool)
        self.matches_all = any(
            len(projections) > 0 and projections.shape[1] == 0
            and () in mode_set
            for _, projections, mode_set, _, _, _, _ in self.clauses
        )
        self.checked_positions = tuple(sorted({
            int(position)
            for _, projections, _, _, _, _, _ in self.clauses
            for projection in projections
            for position in projection
        }))
        # Quantum-change signatures route to the correct source rows, but two
        # different signatures can still yield the same final mode predicate.
        # Intern by the actual projected membership condition, as the expanded
        # plan did after pattern construction.
        # A projection is ordered, whereas the resulting pattern checks are
        # not: (positions 1, 0; modes 2, 7) and (positions 0, 1; modes 7, 2)
        # are the same predicate.  Canonicalize only the checked columns, not
        # a full-width pattern array, to preserve the expanded path's content
        # interning without paying for its wildcard columns.
        if self.matches_all:
            self.semantic_key = (self.width, 'all')
            self.semantic_key_bytes = 32
        else:
            groups = collections.defaultdict(list)
            for signature, projections, _, _, _, _, mode_rows in self.clauses:
                for projection in projections:
                    ordering = np.argsort(projection)
                    positions = tuple(int(position) for position in projection[ordering])
                    if mode_index is None:
                        rows = np.ascontiguousarray(mode_rows[:, ordering])
                        key_dtype = np.dtype((np.void, rows.dtype.itemsize * len(positions)))
                        encoded = np.unique(rows.view(key_dtype).ravel()).tobytes()
                    else:
                        encoded = mode_index.get_canonical_mode_rows(
                            signature, tuple(int(index) for index in ordering)
                        )
                    groups[positions].append(encoded)
            semantic_groups = []
            for positions, blocks in groups.items():
                blocks = tuple(dict.fromkeys(blocks))
                if len(blocks) == 1:
                    encoded = blocks[0]
                else:
                    key_dtype = np.dtype((np.void, np.dtype(np.intp).itemsize * len(positions)))
                    rows = np.concatenate([
                        np.frombuffer(block, dtype=key_dtype) for block in blocks
                    ])
                    encoded = np.unique(rows).tobytes()
                semantic_groups.append((positions, encoded))
            self.semantic_key = (self.width, frozenset(semantic_groups))
            self.semantic_key_bytes = sum(
                len(encoded) + 64 + 8 * len(positions)
                for positions, encoded in semantic_groups
            )

    def __call__(self, state):
        if self.matches_all:
            return True
        for _, projections, _, _, mode_trie, _, _ in self.clauses:
            for projection in projections:
                subtree = mode_trie
                for index in projection:
                    subtree = subtree.get(int(state[int(index)]))
                    if subtree is None:
                        break
                if subtree is True:
                    return True
        return False

    def evaluate(self, states):
        if isinstance(states, np.ndarray) and states.ndim == 1:
            return self(states)
        if len(states) <= self.scalar_state_cutoff:
            return np.fromiter(
                (self(state) for state in states),
                dtype=bool,
                count=len(states)
            )
        return self.evaluate_array(states)

    def evaluate_array(self, states):
        states = np.asarray(states, dtype=np.intp)
        if states.ndim != 2 or states.shape[1] != self.width:
            raise ValueError("degeneracy states do not match compact plan width")
        result = np.zeros(len(states), dtype=bool)
        if self.matches_all:
            result[:] = True
            return result
        for _, projections, _, mode_keys, _, _, _ in self.clauses:
            for projection in projections:
                if len(projection) == 0:
                    continue
                if len(projection) == 1:
                    result |= np.isin(states[:, projection[0]], mode_keys)
                else:
                    # Structured row membership avoids a state x mode-row
                    # broadcast and keeps workspace proportional to the pool.
                    dtype = np.dtype((np.void, len(projection) * np.dtype(np.intp).itemsize))
                    projected = np.ascontiguousarray(states[:, projection])
                    state_keys = projected.view(dtype).ravel()
                    result |= np.isin(state_keys, mode_keys)
                if np.all(result):
                    return result
        return result

    def evaluate_with_value_cache(self, states, cache):
        """Reuse predicate results across pools with equal checked values.

        Only the positions this plan inspects enter the key.  In particular,
        different unchecked mode assignments share one result.  Large pools
        keep the existing vectorized kernel to avoid Python row-key overhead.
        """
        states = np.asarray(states, dtype=np.intp)
        if states.ndim != 2 or states.shape[1] != self.width:
            raise ValueError("degeneracy states do not match compact plan width")
        if self.matches_all:
            return np.ones(len(states), dtype=bool)
        if len(states) > cache.max_batch_rows or cache.max_bytes == 0:
            return self.evaluate_array(states)

        if len(states) >= 4:
            projected = np.ascontiguousarray(
                states[:, self.checked_positions]
            )
            checked_values = projected.tobytes()
            found, result = cache.lookup_block(
                self, len(states), checked_values
            )
            if found:
                return result
            result = self.evaluate(states)
            cache.store_block(self, len(states), checked_values, result)
            cache.stats.degeneracy_value_evaluations += len(states)
            return result

        result = np.empty(len(states), dtype=bool)
        missing = {}
        for index, state in enumerate(states):
            values = tuple(int(state[position]) for position in self.checked_positions)
            pending = missing.get(values)
            if pending is not None:
                pending.append(index)
                continue
            found, value = cache.lookup(self, values)
            if found:
                result[index] = value
            else:
                missing[values] = [index]
        if missing:
            representatives = [indices[0] for indices in missing.values()]
            if len(representatives) > 32:
                values = self.evaluate_array(states[representatives])
            else:
                values = [self(states[index]) for index in representatives]
            for (checked_values, indices), value in zip(missing.items(), values):
                value = bool(value)
                result[indices] = value
                cache.store(self, checked_values, value)
            cache.stats.degeneracy_value_evaluations += len(representatives)
        return result


class DegeneracyChangeIndex:
    """An ordered, immutable index of resonance changes by sorted quanta."""

    canonical_mode_cache_max_bytes = 8 * 1024 ** 2

    def __init__(self, changes):
        groups = collections.OrderedDict()
        for modes, quanta in changes:
            quanta = np.asanyarray(quanta)
            sorting = np.argsort(quanta)
            signature = tuple(
                quanta[index].item()
                    if hasattr(quanta[index], 'item') else quanta[index]
                for index in sorting
            )
            groups.setdefault(signature, []).append(
                tuple(int(modes[index]) for index in sorting)
            )
        self.groups = {}
        self.mode_sets = {}
        self.mode_keys = {}
        self.mode_tries = {}
        self.canonical_mode_cache = collections.OrderedDict()
        self.canonical_mode_cache_bytes = 0
        self.canonical_mode_cache_hits = 0
        self.canonical_mode_cache_misses = 0
        self.nbytes = 0
        self.mode_index_bytes_estimate = 0
        for signature, modes in groups.items():
            block = np.asarray(modes, dtype=np.intp)
            block.setflags(write=False)
            self.groups[signature] = block
            self.nbytes += block.nbytes

    def get(self, signature):
        return self.groups.get(signature)

    def get_canonical_mode_rows(self, signature, ordering):
        """Canonical sparse row bytes, shared across checked-position maps.

        The bounded cache holds only distinct mode-tuple reorderings.  Thus
        interning equivalent predicates does not reconstruct the full
        resonance-mode x projection-map Cartesian product for every plan.
        """
        key = (signature, ordering)
        cached = self.canonical_mode_cache.get(key)
        if cached is not None:
            self.canonical_mode_cache.move_to_end(key)
            if Diagnostics.enabled:
                self.canonical_mode_cache_hits += 1
            return cached[0]
        if Diagnostics.enabled:
            self.canonical_mode_cache_misses += 1
        modes = self.groups[signature]
        rows = np.ascontiguousarray(modes[:, ordering])
        key_dtype = np.dtype((np.void, rows.dtype.itemsize * len(ordering)))
        encoded = np.unique(rows.view(key_dtype).ravel()).tobytes()
        size = len(encoded) + 64 + 8 * len(ordering)
        if size <= self.canonical_mode_cache_max_bytes:
            while (self.canonical_mode_cache and
                   self.canonical_mode_cache_bytes + size > self.canonical_mode_cache_max_bytes):
                _, (_, removed_size) = self.canonical_mode_cache.popitem(last=False)
                self.canonical_mode_cache_bytes -= removed_size
            self.canonical_mode_cache[key] = (encoded, size)
            self.canonical_mode_cache_bytes += size
        return encoded

    def get_mode_membership(self, signature):
        modes = self.groups.get(signature)
        if modes is None:
            return None
        mode_set = self.mode_sets.get(signature)
        if mode_set is None:
            mode_set = frozenset(
                tuple(int(value) for value in row) for row in modes
            )
            self.mode_sets[signature] = mode_set
            mode_trie = {}
            for mode_row in mode_set:
                subtree = mode_trie
                for value in mode_row[:-1]:
                    subtree = subtree.setdefault(value, {})
                if mode_row:
                    subtree[mode_row[-1]] = True
            self.mode_tries[signature] = mode_trie
            self.mode_index_bytes_estimate += len(mode_set) * (
                128 + 32 * modes.shape[1]
            )
            if modes.shape[1] == 1:
                self.mode_keys[signature] = modes[:, 0]
            elif modes.shape[1] > 1:
                key_dtype = np.dtype((np.void, modes.dtype.itemsize * modes.shape[1]))
                self.mode_keys[signature] = np.ascontiguousarray(modes).view(
                    key_dtype
                ).ravel()
            else:
                self.mode_keys[signature] = np.empty(0, dtype=np.intp)
        return mode_set, self.mode_keys[signature], self.mode_tries[signature]


class DegeneracyIdentificationContext:
    """Compile and bind energy-change predicates within a bounded batch cache.

    Only symbolic requirements and predicate plans are retained.  Concrete
    state/permutation masks belong to the indexed evaluator's separate batch
    cache.  Direct callers can create a short-lived context; high-level
    correction workflows share one across their expressions and DAG chunks.
    """

    default_max_items = 4096
    default_max_bytes = 32 * 1024 ** 2
    # The expanded route remains selectable for accuracy and timing controls.
    use_compact_index = True

    def __init__(self, changes, max_items=None, max_bytes=None):
        self.changes = DegeneracyChangeIndex(changes)
        self.compact_index = bool(self.use_compact_index)
        self.max_items = self.default_max_items if max_items is None else int(max_items)
        self.max_bytes = self.default_max_bytes if max_bytes is None else int(max_bytes)
        self._cache = collections.OrderedDict()
        self._cache_bytes = 0
        self._stats = collections.Counter()
        if Diagnostics.enabled:
            self._stats['change_index_bytes'] = self.changes.nbytes

    def _get(self, key, kind):
        entry = self._cache.get(key)
        if entry is None:
            if Diagnostics.enabled:
                self._stats[kind + '_misses'] += 1
            return None
        self._cache.move_to_end(key)
        if Diagnostics.enabled:
            self._stats[kind + '_hits'] += 1
        return entry[0]

    def _put(self, key, value, size):
        size = int(size) + 256  # include a conservative cache-entry allowance
        if size > self.max_bytes or self.max_items <= 0:
            if Diagnostics.enabled:
                self._stats['cache_skips'] += 1
            return value
        while self._cache and (
                len(self._cache) >= self.max_items
                or self._cache_bytes + size > self.max_bytes):
            _, (_, removed_size) = self._cache.popitem(last=False)
            self._cache_bytes -= removed_size
            if Diagnostics.enabled:
                self._stats['cache_evictions'] += 1
        self._cache[key] = (value, size)
        self._cache_bytes += size
        if Diagnostics.enabled:
            self._stats['max_cache_items'] = max(
                self._stats['max_cache_items'], len(self._cache)
            )
            self._stats['max_cache_bytes'] = max(
                self._stats['max_cache_bytes'], self._cache_bytes
            )
        return value

    def _requirement(self, echange, method):
        echange = tuple(echange)
        key = ('requirement', method, echange)
        requirement = self._get(key, 'requirement')
        if requirement is not None:
            return requirement
        nonzero_pos = [i for i, value in enumerate(echange) if value != 0]
        zero_pos = [i for i, value in enumerate(echange) if value == 0]
        sub_e = [echange[i] for i in nonzero_pos]
        sorting = np.argsort(sub_e)
        signature = tuple(sub_e[index] for index in sorting)
        inverse = np.argsort(
            [nonzero_pos[index] for index in sorting] + zero_pos
        )
        splits = np.nonzero(np.diff(signature))[0]
        blocks = (
            [list(itertools.permutations(block)) for block in
             np.split(np.arange(len(sub_e)), splits + 1)]
            if len(splits) > 0 else
            [[np.arange(len(sub_e))]]
        )
        maps = []
        for parts in itertools.product(*blocks):
            subperm = (
                np.concatenate(parts).astype(np.intp, copy=False)
                if len(parts) > 0 else np.empty(0, dtype=np.intp)
            )
            padding = (
                np.full(len(zero_pos), -1, dtype=np.intp)
                if method == 'linear' else
                np.arange(len(sub_e), len(echange), dtype=np.intp)
            )
            maps.append(np.concatenate((subperm, padding))[inverse])
        mode_maps = np.stack(maps)
        mode_maps.setflags(write=False)
        requirement = (signature, mode_maps, len(zero_pos))
        if Diagnostics.enabled:
            self._stats['energy_requirements_compiled'] += 1
        return self._put(key, requirement, mode_maps.nbytes)

    @staticmethod
    def _pattern_count(requirement, modes, nmodes, method):
        _, mode_maps, num_zero = requirement
        if method == 'linear':
            return len(modes) * len(mode_maps)
        return sum(
            math.perm(nmodes - len(set(row)), num_zero) * len(mode_maps)
            for row in modes
        )

    @classmethod
    def _fill_patterns(cls, output, requirements, nmodes, method):
        offset = 0
        for requirement, modes in requirements:
            _, maps, num_zero = requirement
            if method == 'linear':
                count = len(modes) * len(maps)
                target = output[offset:offset + count].reshape(
                    len(modes), len(maps), maps.shape[1]
                )
                # Fill by column to avoid a second M x P x W broadcast array.
                for column in range(maps.shape[1]):
                    sources = maps[:, column]
                    checked = sources >= 0
                    if np.any(checked):
                        target[:, checked, column] = modes[:, sources[checked]]
                offset += count
            else:
                for row in modes:
                    available = [
                        mode for mode in range(nmodes) if mode not in row
                    ]
                    pads = (
                        itertools.permutations(available, num_zero)
                        if num_zero > 0 else [()]
                    )
                    for pad in pads:
                        base = np.concatenate((row, np.asarray(pad, dtype=np.intp)))
                        for mode_map in maps:
                            output[offset] = base[mode_map]
                            offset += 1
        if offset != len(output):
            raise AssertionError('degeneracy pattern row count changed')

    def _compact_clause(self, echange, requirement, modes, projections=None):
        key = ('compact_clause', tuple(echange))
        cached = self._get(key, 'compact_clause')
        if cached is not None:
            return cached
        signature, maps, _ = requirement
        if projections is None:
            checked = np.flatnonzero(maps[0] >= 0)
            projections = np.stack([
                checked[np.argsort(mode_map[checked])]
                for mode_map in maps
            ])
            projections.setflags(write=False)
        mode_set, mode_keys, mode_trie = self.changes.get_mode_membership(signature)
        clause = (
            signature, projections, mode_set, mode_keys, mode_trie,
            len(modes) * len(maps), modes
        )
        if Diagnostics.enabled:
            self._stats['compact_clauses_compiled'] += 1
        return self._put(key, clause, projections.nbytes)

    def _bind_energy_key(self, ekey, nmodes, method, compiled=None):
        key = ('bound', method, int(nmodes), ekey)
        cached = self._get(key, 'bound_plan')
        if cached is not None:
            return cached
        sides = ([], [])
        changes = (PTEnergyChangeProductSum.side_change_iter(ekey) if compiled is None
                   else compiled['sides'][ekey])
        for eside, echange in changes:
            requirement = (self._requirement(echange, method) if compiled is None
                           else compiled['requirements'][method, tuple(echange)][0])
            modes = self.changes.get(requirement[0])
            if modes is not None:
                sides[eside].append((echange, requirement, modes))
                if Diagnostics.enabled:
                    self._stats['resonance_matches'] += len(modes)
        result = []
        size = 0
        for requirements in sides:
            if not requirements:
                result.append(None)
                continue
            width = len(ekey[0]) - 1  # leading element is the energy side
            count = sum(
                self._pattern_count(requirement, modes, nmodes, method)
                for _, requirement, modes in requirements
            )
            if count == 0:
                result.append(None)
                continue
            if self.compact_index and method == 'linear':
                clauses = [
                    self._compact_clause(echange, requirement, modes,
                        None if compiled is None else compiled['requirements'][method, tuple(echange)][1])
                    for echange, requirement, modes in requirements
                ]
                plan = CompactDegeneracyTestPlan(
                    clauses, width, mode_index=self.changes
                )
                result.append(plan)
                if Diagnostics.enabled:
                    self._stats['pattern_rows_avoided'] += count
                    self._stats['predicate_plans_built'] += 1
                size += sum(clause[1].nbytes + 256 for clause in clauses)
                size += plan.semantic_key_bytes
                continue
            patterns = np.full((count, width), -1, dtype=np.intp)
            self._fill_patterns(
                patterns,
                [(requirement, modes) for _, requirement, modes in requirements],
                nmodes, method
            )
            plan = DegeneracyTestPlan(patterns)
            result.append(plan)
            if Diagnostics.enabled:
                self._stats['pattern_rows_materialized'] += count
                self._stats['predicate_plans_built'] += 1
            size += patterns.nbytes + plan.patterns.nbytes + plan.checks.nbytes
            # The scalar trie is compiled lazily; reserve room for its nodes.
            size += count * (64 + 16 * width)
        bound = tuple(result)
        return self._put(key, bound, size)

    def _intern_predicate(self, plan, by_content, by_identity):
        """Reuse equal plans within one bind, as the legacy builder did.

        The indexed evaluator caches shared-pool masks by predicate identity.
        Distinct energy keys can produce the same final pattern set, so the
        bound-plan cache alone does not preserve that reuse.  Hash the existing
        read-only pattern buffer without a tuple or byte-array copy, then
        verify equality before sharing to make digest collisions harmless.
        Both local tables are capped and disappear after this bind call.
        """
        if plan.fallback_tests:
            return plan
        identity = by_identity.get(plan)
        if identity is not None and identity() is not None:
            if Diagnostics.enabled:
                self._stats['predicate_intern_identity_hits'] += 1
            return identity()

        if isinstance(plan, CompactDegeneracyTestPlan):
            key = ('compact', plan.semantic_key)
            canonical = by_content.get(key, plan)
            if canonical is plan:
                if Diagnostics.enabled:
                    self._stats['predicate_intern_misses'] += 1
                if len(by_identity) < self.max_items:
                    by_content[key] = plan
                else:
                    if Diagnostics.enabled:
                        self._stats['predicate_intern_skips'] += 1
            else:
                if Diagnostics.enabled:
                    self._stats['predicate_intern_content_hits'] += 1
            if len(by_identity) < self.max_items:
                by_identity[plan] = weakref.ref(canonical)
            return canonical

        patterns = plan.patterns
        if not patterns.flags.c_contiguous:
            patterns = np.ascontiguousarray(patterns)
        fingerprint = hashlib.blake2b(patterns.data, digest_size=16).digest()
        key = (patterns.shape, fingerprint)
        canonical = plan
        for candidate in by_content.get(key, ()):
            if (
                    candidate.scalar_state_cutoff == plan.scalar_state_cutoff
                    and candidate.compare_workspace_bytes == plan.compare_workspace_bytes
                    and np.array_equal(candidate.patterns, plan.patterns)
                    and np.array_equal(candidate.checks, plan.checks)
            ):
                canonical = candidate
                if Diagnostics.enabled:
                    self._stats['predicate_intern_content_hits'] += 1
                break
        else:
            if Diagnostics.enabled:
                self._stats['predicate_intern_misses'] += 1
            if len(by_identity) < self.max_items:
                by_content.setdefault(key, []).append(plan)
            else:
                if Diagnostics.enabled:
                    self._stats['predicate_intern_skips'] += 1

        if len(by_identity) < self.max_items:
            by_identity[plan] = weakref.ref(canonical)
        return canonical

    def bind(self, expr, nmodes, method):
        if Diagnostics.enabled:
            self._stats['identification_calls'] += 1
        all_degs = ({}, {})
        by_content = {}
        by_identity = weakref.WeakKeyDictionary()
        execution = getattr(expr, '_compiled_execution', None)
        compiled = None if execution is None else execution['degeneracy']
        terms = (compiled['terms'].items() if compiled is not None else
                 ((key, value.terms) for key, value in expr.terms.items() if isinstance(value, PTEnergyChangeProductSum)))
        for cinds, energy_keys in terms:
            for ekey in energy_keys:
                left, right = (self._bind_energy_key(ekey, nmodes, method) if compiled is None else
                               self._bind_energy_key(ekey, nmodes, method, compiled))
                if left is not None:
                    all_degs[0].setdefault(cinds, {})[ekey] = (
                        self._intern_predicate(left, by_content, by_identity)
                    )
                if right is not None:
                    all_degs[1].setdefault(cinds, {})[ekey] = (
                        self._intern_predicate(right, by_content, by_identity)
                    )
        return all_degs

    def stats(self):
        return dict(self._stats, cache_items=len(self._cache),
                    cache_bytes=self._cache_bytes,
                    mode_index_bytes_estimate=(
                        self.changes.mode_index_bytes_estimate
                    ),
                    canonical_mode_cache_bytes=(
                        self.changes.canonical_mode_cache_bytes
                    ),
                    canonical_mode_cache_hits=(
                        self.changes.canonical_mode_cache_hits
                    ),
                    canonical_mode_cache_misses=(
                        self.changes.canonical_mode_cache_misses
                    ))


class AnalyticPerturbationTheorySolver:
    """
    A re-attempt at using the recursive expressions
    to provide simpler code for getting APT expressions
    """
    def __init__(self, hamiltonian_expansion, logger=None, checkpoint=None,
                 allowed_terms=None,
                 allowed_coefficients=None,
                 disallowed_coefficients=None,
                 allowed_energy_changes=None,
                 intermediate_normalization=None,
                 polynomial_representation='path', resonance_condition=None
                 ):
        self.hamiltonian_expansion = hamiltonian_expansion
        self.logger = Logger.lookup(logger)
        self.checkpoint = Checkpointer.build_canonical(checkpoint)
        self.allowed_terms = allowed_terms
        self.allowed_coefficients = allowed_coefficients
        self.disallowed_coefficients = disallowed_coefficients
        self.allowed_energy_changes = allowed_energy_changes
        self.intermediate_normalization = intermediate_normalization
        self.resonance_condition = resonance_condition
        self._conditional_stacks = {}
        if polynomial_representation not in {'eager', 'path'}:
            raise ValueError("unknown polynomial representation {}".format(polynomial_representation))
        self.polynomial_representation = polynomial_representation
        self._op_maps = {}

    @classmethod
    def from_order(cls, order, internals=True, logger=None, checkpoint=None,
                   allowed_terms=None,
                   allowed_coefficients=None,
                   disallowed_coefficients=None,
                   allowed_energy_changes=None,
                   intermediate_normalization=None,
                   polynomial_representation='path', resonance_condition=None,
                   combine_kinetic_coefficients=True):
        logger = Logger.lookup(logger)
        if order < 2:
            raise ValueError("why")
        if internals:
            hamiltonian_terms = [
                HamiltonianExpansionTerm(
                    [
                        ['x'] * (o + 2),
                        ['p'] + ['x'] * o + ['p']
                    ] + (
                        []
                            if o < 2 else
                        [
                            ['x'] * (o - 2)
                        ]  # V' term has different shape from rest
                    ),
                    order=o,
                    identities=[0, 1] + ([] if o < 2 else [2]), # V, G, G'
                    combine_kinetic_coefficients=combine_kinetic_coefficients,
                    logger=logger,
                    allowed_coefficients=allowed_coefficients,
                    disallowed_coefficients=disallowed_coefficients,
                    polynomial_representation=polynomial_representation
                )
                for o in range(order + 1)
            ]
        else:
            # raise NotImplementedError("issues persist with Coriolis symmetries")
            hamiltonian_terms = [
                HamiltonianExpansionTerm(
                    [
                        ['x'] * (o + 2)
                    ] + (
                        [['p', 'p']] if o == 0 else []
                    )+ (
                        []
                        if o < 2 else
                        [
                            ['x', 'p'] + ['x'] * (o-1) + ['p'], # Coriolis
                            # ['p'] + ['x'] * (o-1) + ['p', 'x'],
                            ['x'] * (o - 2) # Watson
                        ]
                    ),
                    order=o,
                    identities=([0, 1] if o == 0 else [0]) + ([] if o < 2 else [3, 4]),  # V, G, G'
                    logger=logger,
                    allowed_coefficients=allowed_coefficients,
                    disallowed_coefficients=disallowed_coefficients,
                    polynomial_representation=polynomial_representation
                )
                for o in range(order + 1)
            ]

        return cls(hamiltonian_terms, logger=logger, checkpoint=checkpoint,
                   allowed_terms=allowed_terms,
                   allowed_coefficients=allowed_coefficients,
                   disallowed_coefficients=disallowed_coefficients,
                   allowed_energy_changes=allowed_energy_changes,
                   intermediate_normalization=intermediate_normalization,
                   polynomial_representation=polynomial_representation,
                   resonance_condition=resonance_condition)

    def modify_hamiltonian(self, hamiltonian_corrections):
        new_expansion = ReexpressedHamiltonian.prep_expansion(
            self.hamiltonian_expansion,
            hamiltonian_corrections,
            logger=self.logger,
            allowed_coefficients=self.allowed_coefficients,
            disallowed_coefficients=self.disallowed_coefficients
        )

        return type(self)(
            new_expansion,
            logger=self.logger, checkpoint=self.checkpoint,
            allowed_terms=self.allowed_terms,
            allowed_coefficients=self.allowed_coefficients,
            disallowed_coefficients=self.disallowed_coefficients,
            allowed_energy_changes=self.allowed_energy_changes,
            intermediate_normalization=self.intermediate_normalization,
            polynomial_representation=self.polynomial_representation,
            resonance_condition=self.resonance_condition
        )

    def get_conditional_stack(self, condition=None, *, block_recurrence=False):
        condition = self.resonance_condition if condition is None else condition
        if condition is None:
            raise ValueError('no resonance condition was supplied')
        key = (condition, 'block-recurrence') if block_recurrence else condition
        if key not in self._conditional_stacks:
            self._conditional_stacks[key] = ConditionalTransformationStack(
                self, condition, block_recurrence=block_recurrence)
        return self._conditional_stacks[key]

    def _conditional_term(self,kind,order,**options):
        if self.resonance_condition is None:
            return None
        stack = self.get_conditional_stack().extend(order)
        if kind in {'energy','hamiltonian'}:
            term = stack.hamiltonian[order]
        else:
            series = (stack.intermediate_wavefunctions(order) if self.intermediate_normalization and kind!='overlap'
                      else stack.wavefunctions)
            term = series[order]
            if kind == 'wavefunction':
                term = stack.restrict(term,retain=False)
            elif kind == 'overlap':
                term = stack.restrict(term)
        return CorrectionView(term,diagonal=kind=='energy',**options)
    filtered_block_recurrence = False
    resonance_independent = False

    _op_maps = {}
    def get_correction(self, key, cls, order, **kw):
        corr_key = (cls, key, order, bool(self.filtered_block_recurrence), bool(self.intermediate_normalization),
                    bool(self.resonance_independent))
        for option in ('operator_type','base_index'):
            if kw.get(option) is not None:corr_key += ((option,kw[option]),)
        term_filter = kw.get('allowed_terms')
        if term_filter is None:term_filter = self.allowed_terms
        if term_filter is not None:
            # Distinct filtered operators must not reuse the first correction.
            term_filter = tuple(sorted(set(tuple(t) for t in term_filter)))
            kw['allowed_terms'] = term_filter
            corr_key += (('allowed_terms',term_filter),)
        corr = self._op_maps.get(corr_key, None)
        if corr is None:
            for k,v in [
                ["logger", self.logger],
                ["checkpoint", self.checkpoint],
                ["intermediate_normalization", self.intermediate_normalization],
                ["allowed_terms", self.allowed_terms],
                ["allowed_energy_changes", self.allowed_energy_changes],
                ["allowed_coefficients", self.allowed_coefficients],
                ["disallowed_coefficients", self.disallowed_coefficients],
                ["polynomial_representation", self.polynomial_representation]
            ]:
                if kw.get(k, None) is None:
                    kw[k] = v
            corr = cls(self, order, **kw)
            self._op_maps[corr_key] = corr
        return corr

    def block_hamiltonian_correction(self, order, **kw):
        parent = self
        if not self.filtered_block_recurrence:
            parent = copy.copy(self)
            parent.filtered_block_recurrence = True
        return parent.get_correction('block-H', BlockHamiltonianCorrection, order, **kw)

    @staticmethod
    def action_expansion_powers(order, *, max_degree=None):
        """Ordered positive exponent patterns through the VPT action degree.

        A standard even correction n has degree at most n/2+1 in actions.
        Distinct permutations such as (2,1) and (1,2) represent different
        monomials on the same ordered pair of physical mode indices.
        """
        if isinstance(order,(bool,np.bool_)) or not isinstance(order,(int,np.integer)) or order < 0:
            raise ValueError('action expansion order must be a nonnegative integer')
        if max_degree is None:max_degree=int(order)//2+1
        if isinstance(max_degree,(bool,np.bool_)) or not isinstance(max_degree,(int,np.integer)) or max_degree<0:
            raise ValueError('action degree must be a nonnegative integer')
        powers = [()]
        for degree in range(1,int(max_degree)+1):
            for group in IntegerPartitioner.partitions(degree):
                for partition in group:
                    powers.extend(tuple(int(p) for p in permutation)
                                  for permutation in UniquePermutations(partition).permutations())
        return tuple(powers)

    def _action_expansion_source(self, correction, order=2, *, resonance_correction=False):
        if not isinstance(resonance_correction,(bool,np.bool_)):
            raise TypeError('resonance_correction must be boolean')
        if correction not in {'energy','hamiltonian'}:
            raise ValueError("action source must be 'energy' or 'hamiltonian'")
        if any(getattr(self,name) is not None for name in (
                'allowed_terms','allowed_coefficients','disallowed_coefficients','allowed_energy_changes')):
            raise ValueError('unrestricted action expansion requires unfiltered correction definitions')
        key=('action-source',correction,int(order),bool(resonance_correction),
             bool(self.intermediate_normalization),self.polynomial_representation,
             tuple(self.hamiltonian_expansion[:order+1]))
        cached=self._op_maps.get(key)
        source = copy.copy(self) if cached is None else cached
        if cached is None:source._op_maps,source._conditional_stacks = {},{}
        source.resonance_condition = None
        source.resonance_independent = not resonance_correction
        source.filtered_block_recurrence = resonance_correction or correction == 'hamiltonian'
        source.intermediate_normalization = self.intermediate_normalization if resonance_correction else correction == 'energy'
        # Missing higher Hamiltonian orders are zero; the lower orders can
        # still generate corrections to arbitrary perturbative order.
        source.hamiltonian_expansion = source.hamiltonian_expansion[:order+1]
        for n,term in enumerate(source.hamiltonian_expansion):
            if not isinstance(term,(OperatorExpansionTerm,ZeroTerm)):
                raise NotImplementedError('action expansion requires the native Hamiltonian expansion')
            if isinstance(term,OperatorExpansionTerm):
                if max(map(len,term.terms),default=0) > n+2:
                    raise ValueError('Hamiltonian order exceeds the VPT action degree bound')
                if any(len(t)%2 != n%2 for t in term.terms):
                    raise ValueError('Hamiltonian terms must have the native VPT order parity')
                if any(getattr(term,name) is not None for name in (
                        'allowed_terms','allowed_coefficients','disallowed_coefficients','allowed_energy_changes')):
                    raise ValueError('unrestricted action expansion requires unfiltered Hamiltonian terms')
        self._op_maps[key]=source
        return source

    def action_expansion_terms(self, order=2, correction='energy', *, resonance_correction=False):
        """Reusable symbolic action coefficients indexed by (PT order, powers).

        By default the coefficients are resonance independent. With
        resonance_correction=True they represent the exact retained-domain
        residual, including its full native normalization and F feedback.
        Both sources remain symbolic and reusable across numerical bindings.

        The zero-change block Hamiltonian uses singleton P and its full Q
        complement. No runtime resonance selection enters either source.
        Kernels retain force tensors and harmonic frequency denominators.
        """
        powers = self.action_expansion_powers(order)
        source = self._action_expansion_source(correction,order,resonance_correction=resonance_correction)
        generator = (source.energy_correction if correction == 'energy'
                     else source.block_hamiltonian_correction)
        result = {}
        for n in range(0,int(order)+1,2):
            polynomial = generator(n).get_poly_terms(())
            if nput.is_numeric(polynomial):
                polynomial = SqrtChangePoly(PTTensorCoeffProductDAG.from_sum(
                    PTTensorCoeffProductSum({():polynomial},ndim=0)),(),())
            selected_powers=powers
            if resonance_correction:
                _,polynomial=polynomial.split_change_domains()
                if nput.is_zero(polynomial):continue
                from .Evaluators import _polynomial_degree
                selected_powers=self.action_expansion_powers(order,max_degree=_polynomial_degree(polynomial.poly_obj)//2)
            for pattern in selected_powers:
                result[n,pattern] = polynomial.project_coefficients(pattern,shift=-.5,
                    preserve_domains=resonance_correction)
        return result

    def operator_action_expansion_terms(self, order, change, operator_type='transition_moment',
                                        base_index=None, *, shift=-.5, resonance_correction=False):
        """Symbolic operator coefficients keyed by (order, fixed powers, free powers).

        resonance_correction=True selects the exact resonance residual instead
        of the unrestricted operator, preserving integer change domains.

        Changed-mode slots remain external even when their powers are zero.
        Positive free powers denote additional distinct spectator-mode slots.
        Returned SqrtChangePoly kernels retain the exact transition ladder.
        This is a coefficient extraction from native equations, not a fit.
        """
        self.action_expansion_powers(order)  # shared order validation
        change = tuple(change)
        if any(int(c) != c or c == 0 for c in change) or tuple(self.hamiltonian_expansion[0].sorted_changes(change)) != change:
            raise ValueError('operator expansion requires a canonical nonzero change tuple, or ()')
        from .Evaluators import _polynomial_degree
        source = self._action_expansion_source('hamiltonian',order,resonance_correction=resonance_correction)
        terms = {}
        for n in range(order+1):
            polynomial = source.operator_correction(n,operator_type=operator_type,
                base_index=base_index).get_poly_terms(change)
            if nput.is_zero(polynomial):
                continue
            if resonance_correction:
                _,polynomial=polynomial.split_change_domains()
                if nput.is_zero(polynomial):continue
            degree = max(0,(_polynomial_degree(polynomial.poly_obj)-sum(map(abs,change)))//2)
            free_patterns = self.action_expansion_powers(order,max_degree=degree)
            for fixed_powers in itertools.product(range(degree+1),repeat=len(change)):
                for powers in free_patterns:
                    if sum(fixed_powers)+sum(powers) <= degree:
                        terms[n,fixed_powers,powers] = polynomial.project_coefficients(powers,shift,
                            fixed_powers=fixed_powers,preserve_domains=resonance_correction)
        return terms

    def block_main_wavefunction_correction(self, order):
        if order == 0:
            return IdentityTerm(logger=self.logger, polynomial_representation=self.polynomial_representation)
        return operator_sum(self.wavefunction_correction(order),
                            CorrectionView(self.overlap_correction(order), diagonal=True))

    def shifted_hamiltonian_correction(self, order, **kw):
        return self.get_correction('dH', ShiftedHamiltonianCorrection, order, **kw)

    def energy_correction(self, order, **kw):
        if self.resonance_condition is not None:
            return self._conditional_term('energy',order,**kw)
        return self.get_correction('energy', EnergyCorrection, order, **kw)

    def wavefunction_correction(self, order, **kw):
        if self.resonance_condition is not None:
            return self._conditional_term('wavefunction',order,**kw)
        return self.get_correction('proj_wfn', WavefunctionCorrection, order, **kw)

    def overlap_correction(self, order, degenerate_changes=None, **kw):
        if self.resonance_condition is not None:
            return self._conditional_term('overlap',order,**kw)
        if degenerate_changes is not None:
            key = ('ov', tuple(sorted(tuple(c) for c in degenerate_changes)))
        else:
            key = 'ov'
        return self.get_correction(key, WavefunctionOverlapCorrection, order, degenerate_changes=degenerate_changes, **kw)

    def full_wavefunction_correction(self, order, **kw):
        if self.resonance_condition is not None:
            return self._conditional_term('full',order,**kw)
        return self.get_correction('wfn', FullWavefunctionCorrection, order, **kw)

    def operator_correction(self, order, operator_type=None, **kw):
        if self.resonance_condition is not None:
            operators = self.operator_expansion_terms(order,operator_type=operator_type,
                base_index=kw.pop('base_index',None),polynomial_representation=self.polynomial_representation)
            term = self.get_conditional_stack().transform_operator(operators,order,
                                                                  terms=kw.get('allowed_terms'))[order]
            return CorrectionView(term,**kw)
        return self.get_correction('op', OperatorCorrection, order, operator_type=operator_type, **kw)
        # return OperatorCorrection(self, order, operator_type=operator_type, logger=self.logger)
    def operator_degenerate_correction(self, order, /, degenerate_changes, operator_type=None, **kw):
        if self.resonance_condition is not None:
            return {}
        key = ('op', tuple(sorted(tuple(c) for c in degenerate_changes)))
        return self.get_correction(key, OperatorDegenerateCorrection, order, operator_type=operator_type, degenerate_changes=degenerate_changes, **kw)
        # return OperatorCorrection(self, order, operator_type=operator_type, logger=self.logger)

    def reexpressed_hamiltonian(self, order, **kw):
        if self.resonance_condition is not None:
            return self._conditional_term('hamiltonian',order,**kw)
        return self.get_correction('reH', ReexpressedHamiltonian, order, **kw)
    def reexpressed_hamiltonian_degenerate_correction(self, order, /, degenerate_changes, **kw):
        if self.resonance_condition is not None:
            return {}
        key = ('reH', tuple(sorted(tuple(c) for c in degenerate_changes)))
        return self.get_correction(key, ReexpressedHamiltonianDegenerateCorrection, order, degenerate_changes=degenerate_changes, **kw)

    operator_expansion_index = 5
    @classmethod
    def operator_expansion_terms(cls, order, logger=None, base_index=None, operator_type=None,
                                 polynomial_representation='path'):
        if base_index is None:
            base_index = cls.operator_expansion_index

        base_terms = []
        padding = 0
        index_padding = 0
        if operator_type is None:
            padding = 1
        elif nput.is_numeric(operator_type):
            padding = operator_type
        elif isinstance(operator_type, str):
            if operator_type == 'transition_moment':
                padding = 1
                index_padding = 1
                base_terms.append(
                    OperatorExpansionTerm(
                        [],
                        order=0,
                        index=0,
                        identities=[base_index],
                        logger=logger,
                        polynomial_representation=polynomial_representation
                    )
                )
            else:
                raise ValueError("unknown operator type {}".format(operator_type))
        else:
            raise ValueError("unknown operator type {}".format(operator_type))

        shift_order = len(base_terms)
        expansion_terms = [
            OperatorExpansionTerm(
                [['x'] * (o + padding)],
                order=o,
                index=index_padding+o,
                identities=[base_index],
                logger=logger,
                polynomial_representation=polynomial_representation
            )
            for o in range(order+1)
        ]

        final_terms = [
            PerturbationTheoryTermSum(b,e)
            for b,e in zip(base_terms, expansion_terms)
        ] + expansion_terms[shift_order:]

        # raise Exception(final_terms, min_order, expansion_terms, shift_order)

        return final_terms

    @classmethod
    def clear_caches(cls):
        cls._op_maps.clear()
        PerturbationTheoryTermSum._sum_cache.clear()
        _clear_projected_caches()
        _clear_conditional_caches()
        ProductPTPolynomial._cache.clear()
        ProductPTPolynomial._prod_poly_cache.clear()
        PerturbationTheoryTermProduct._cache.clear()
        PerturbationTheoryTermProduct._poly_product_cache.clear()
        ShiftedEnergyBaseline._cache.clear()
        PerturbationOperator._cache.clear()
        ScaledPerturbationTheoryTerm._cache.clear()
        OperatorExpansionTerm._poly_cache.clear()
        PerturbationTheoryExpressionEvaluator._poly_cache = PerturbationTheoryExpressionEvaluator.get_cache()
        PerturbationTheoryExpressionEvaluator._ecoeff_cache = PerturbationTheoryExpressionEvaluator.get_cache()
        try:
            from .Evaluators import IndexedBlockEvaluator
        except ImportError:
            pass
        else:
            IndexedBlockEvaluator.clear_plan_cache()
        PolyPath.clear_caches()
        PTTensorCoeffProductDAG.clear_caches()

    def polynomial_cache_info(self):
        """Return backend-specific counters useful for path/eager timing comparisons."""
        return {
            'representation': self.polynomial_representation,
            'product_results': len(PerturbationTheoryTermProduct._poly_product_cache),
            **PolyPath.cache_info(),
            **PTTensorCoeffProductDAG.cache_info()
        }


class PolynomialInterface(metaclass=abc.ABCMeta):
    """
    Provides a basic interface to allow for the uniform manipulation
    of objects that dispatch down to some form of scalar multiplied by
    a sum of polynomials
    """

    def format_expr(self) -> str:
        ...

    @property
    @abc.abstractmethod
    def ndim(self) -> int:
        ...
    # @property
    # @abc.abstractmethod
    # def order(self) -> 'tuple[int]':
    #     ...
    @abc.abstractmethod
    def audit(self, target, ignore_constants=True):
        ...
    @abc.abstractmethod
    def ensure_dimension(self, ndim) -> 'Self':
        ...
    def align_dimensions(self, other) -> 'tuple[Self, Self]':
        max_dim = max(self.ndim, other.ndim)
        return self.ensure_dimension(max_dim), other.ensure_dimension(max_dim)
    # @property
    # @abc.abstractmethod
    # def order(self):
    #     ...
    @abc.abstractmethod
    def shift(self, shift)->'Self':
        ...
    @abc.abstractmethod
    def scale(self, scaling)->'Self':
        ...
    @abc.abstractmethod
    def permute(self, perm, check_perm=True, allow_padding=False)->'Self':
        ...

    @abc.abstractmethod
    def mul_simple(self, other:'PolynomialInterface')->'PolynomialInterface':
        ...
    @abc.abstractmethod
    def rmul_simple(self, other: 'PolynomialInterface') -> 'PolynomialInterface':
        ...
    @abc.abstractmethod
    def mul_along(self, other:'PolynomialInterface', inds, remainder=None, mapping=None)->'PolynomialInterface':
        ...
    @abc.abstractmethod
    def rmul_along(self, other:'PolynomialInterface', inds, remainder=None, mapping=None)->'PolynomialInterface':
        ...

    @abc.abstractmethod
    def combine(self, **kwargs) -> 'Self':
        ...
    # @abc.abstractmethod
    # def condense(self, inds=None, return_inds=False, check_inds=True):
    #     ...
    @abc.abstractmethod
    def mutate(self, *args, **kwargs) -> 'Self':
        ...

class PolyAtom:
    """An immutable, interned one-dimensional polynomial factor."""

    _cache = {}

    def __new__(cls, coeffs, shift=0):
        coeffs = np.asanyarray(coeffs)
        key = (coeffs.dtype.str, coeffs.shape, coeffs.tobytes(), int(shift))
        atom = cls._cache.get(key)
        if atom is None:
            atom = super().__new__(cls)
            atom._key = key
            atom._base_coeffs = np.array(coeffs, copy=True)
            atom._base_coeffs.flags.writeable = False
            atom.shift = int(shift)
            atom._coeffs = None
            cls._cache[key] = atom
        return atom

    def __fingerprint__(self):
        return (type(self).__module__,type(self).__qualname__,self._base_coeffs,self.shift)

    @property
    def coeffs(self):
        if self._coeffs is None:
            self._coeffs = (
                self._base_coeffs
                    if self.shift == 0 else
                DensePolynomial.compute_shifted_coeffs(self._base_coeffs, self.shift)
            )
        return self._coeffs

    @property
    def order(self):
        return len(self._base_coeffs) - 1

    def shifted(self, shift):
        return self if shift == 0 else type(self)(self._base_coeffs, self.shift + shift)

    def __hash__(self):
        return hash(self._key)

    def __eq__(self, other):
        return self is other or isinstance(other, type(self)) and self._key == other._key

    def __reduce__(self):
        return type(self), (self._base_coeffs, self.shift)


class PolyAxis:
    """A canonical product of :class:`PolyAtom` objects along one mode."""

    _cache = {}
    _materialization_count = 0

    def __new__(cls, atoms=()):
        atoms = tuple(sorted(atoms, key=lambda a: a._key))
        axis = cls._cache.get(atoms)
        if axis is None:
            axis = super().__new__(cls)
            axis.atoms = atoms
            axis._coeffs = None
            cls._cache[atoms] = axis
        return axis

    def __fingerprint__(self):
        return (type(self).__module__,type(self).__qualname__,self.atoms)

    @property
    def order(self):
        return sum(a.order for a in self.atoms)

    @property
    def coeffs(self):
        if self._coeffs is None:
            if Diagnostics.enabled:
                type(self)._materialization_count += 1
            coeffs = np.array([1])
            for atom in self.atoms:
                coeffs = scipy.signal.convolve(coeffs, atom.coeffs)
            self._coeffs = coeffs
        return self._coeffs

    def multiply(self, other):
        return type(self)(self.atoms + other.atoms)

    def shifted(self, shift):
        return self if shift == 0 else type(self)(a.shifted(shift) for a in self.atoms)

    def __hash__(self):
        return hash(self.atoms)

    def __eq__(self, other):
        return self is other or isinstance(other, type(self)) and self.atoms == other.atoms

    def __reduce__(self):
        return type(self), (self.atoms,)


class PolyTerm:
    """One separable product in a :class:`PolyPath` sum."""

    _cache = {}

    def __new__(cls, axes, steps=0):
        axes = tuple(axes)
        key = (axes, 0 if steps is None else int(steps))
        term = cls._cache.get(key)
        if term is None:
            term = super().__new__(cls)
            term.axes, term.steps = key
            cls._cache[key] = term
        return term

    @property
    def order(self):
        return tuple(a.order for a in self.axes)

    def __hash__(self):
        return hash((self.axes, self.steps))

    def __eq__(self, other):
        return self is other or (
            isinstance(other, type(self))
            and self.steps == other.steps
            and self.axes == other.axes
        )

    def __reduce__(self):
        return type(self), (self.axes, self.steps)

class TreeSerializer:
    @classmethod
    def serialize_iterable(cls, iterable, primitive_test, concat, track_shapes=True):
        if not track_shapes:
            shapes = None
        else:
            shapes = []
        flats = []
        if primitive_test(iterable[0][0]):
            if track_shapes:
                shapes = [len(iterable)] + [len(i) for i in iterable]
            flats = iterable
        else:
            if track_shapes:
                shapes.append(len(iterable))
            for obj in iterable:
                shape, flat = cls.serialize_iterable(obj, primitive_test, concat)
                if track_shapes:
                    shapes.extend(shape)
                flats.append(flat)
        flats = concat(flats)
        return shapes, flats

    @classmethod
    def build_dict_trees(cls, dict_obj):
        vals = {}
        if isinstance(next(iter(dict_obj.values())), dict):
            keys = [
                list(dict_obj.keys())
            ]  # one set of trees per tree depth
            for obj in dict_obj.values():
                subkeys, subvals = cls.build_dict_trees(obj)
                if subkeys is not None:
                    for d,k in enumerate(subkeys):
                        if len(keys) <= d+1:
                            keys.append([])
                        keys[d+1].append(k)
                for k,v in subvals.items():
                    if k not in vals: vals[k] = []
                    vals[k].append(v)
        else:
            keys = None
            vals = dict_obj
            # for k,v in dict_obj.items():
            #     if k not in vals: vals[k] = []
            #     vals[k].append(v)
        return keys, vals

    @classmethod
    def serialize_tree_dict(cls, dict_obj,
                            key_primitive_test=None, vals_primitive_test=None, concat=None):
        if key_primitive_test is None: key_primitive_test = cls.default_prim
        if vals_primitive_test is None: vals_primitive_test = cls.default_prim
        if concat is None: concat = cls.default_concat

        keys, vals = cls.build_dict_trees(dict_obj)


        key_trees = []
        for k in keys:
            shape, flat = cls.serialize_iterable(k, key_primitive_test, concat)
            key_trees.append([shape, flat])

        vals_trees = {}
        tree_shape = None
        for i,(k, v) in enumerate(vals.items()):
            shape, flat = cls.serialize_iterable(v, vals_primitive_test, concat, track_shapes=i == 0)
            if tree_shape is None: tree_shape = shape
            vals_trees[k] = flat

        return key_trees, tree_shape, vals_trees

    @staticmethod
    def default_concat(arrays):
        if nput.is_numeric(arrays[0]):
            return np.asanyarray(arrays)
        else:
            if isinstance(arrays[0][0], np.ndarray):
                max_dim = max(len(y) for x in arrays for y in x)
                arrays = [
                    [np.pad(a, [0, max_dim - len(a)]) for a in subarr]
                    for subarr in arrays
                ]
            return np.concatenate(arrays, axis=0)

    @staticmethod
    def default_prim(obj):
        return isinstance(obj, (str, np.ndarray)) or nput.is_numeric(obj)

    @classmethod
    def deserialize_subiterable(cls, shape, flat, i, pad, depth):
        if depth == 0:
            block_size = shape[i]
            if isinstance(flat, dict):
                return {key:buff[pad:pad+block_size] for key, buff in flat.items()}
            else:
                return flat[pad:pad+block_size], i+1, pad + block_size
        else:
            tree = []
            queue = collections.deque([[shape[i], depth-1, tree]])
            current = i+1
            while queue:
                rem, d, stack = queue.pop()
                if rem > 1:
                    # more blocks to process, so push back onto queue to be processed later
                    queue.append([rem-1, d, stack])

                block_size = shape[current]
                if d == 1:
                    current += 1
                    substack = []
                    for bs in shape[current:current + block_size]:
                        if isinstance(flat, dict):
                            substack.append({key:buff[pad:pad+bs] for key, buff in flat.items()})
                        else:
                            substack.append(flat[pad:pad+bs])
                        pad += bs
                        # consumed += 1
                    current += block_size
                    stack.append(substack)
                    continue
                subtree = []
                stack.append(subtree)
                queue.append([block_size, d-1, subtree])
                current += 1

            return tree, current, pad

    @classmethod
    def deserialize_iterable(cls, shape, flat, depth):
        if depth == 1:
            trees = []
            pad = 0
            current = 1
            for bs in shape[current:current+shape[0]]:
                if isinstance(flat, dict):
                    trees.append({key: buff[pad:pad + bs] for key, buff in flat.items()})
                else:
                    trees.append(flat[pad:pad + bs])
                pad += bs
        else:
            nblock = 0
            pad = 0
            # max_block = len(shape)
            # for i in range(max_block): # max possible test
            trees, nblock, pad = cls.deserialize_subiterable(shape, flat, nblock, pad, depth)
                # trees.append(subtree)
                # nblock += 1
                # if nblock >= max_block - 1:
                #     break
            # else:
            #     raise ValueError("shouldn't be able to directly exhaust this...")
        return trees
    @classmethod
    def _tuplate(cls, key):
        return tuple(
            tuple(k) if nput.is_numeric(k[0]) else cls._tuplate(k)
            for k in key
        )
    @classmethod
    def stitch_dict(cls, iterables):
        stitched = {}
        for p in zip(*iterables):
            if len(p) == 3:
                k, subkey, v = p
                new_dict = {}
                for i,sk in enumerate(subkey):
                    if isinstance(sk, (list, np.ndarray)): sk = cls._tuplate(sk)
                    new_dict[sk] = {
                        vk:vv[i]
                        for vk, vv in v.items()
                    }
                v = new_dict
            else:
                k = p[0]
                v = cls.stitch_dict(p[1:])
            # if isinstance(k, list): k = k[0]
            if isinstance(k, (list, np.ndarray)): k = cls._tuplate(k)
            stitched[k] = v
        return stitched
    @classmethod
    def rebuild_tree_dict(cls, key_trees, tree_shape, val_buffers):
        # shape = None
        total_depth = len(key_trees)
        keys = []
        for d,(shape, buffer) in enumerate(key_trees):
            st = cls.deserialize_iterable(shape, buffer, d+2)
            keys.append(st)
        vals = cls.deserialize_iterable(tree_shape, val_buffers, total_depth - 1)

        return cls.stitch_dict(keys + [vals])

class ProductPTPolynomial(PolynomialInterface):
    """
    TODO: include prefactor term so we can divide out energy changes
    """
    def __init__(self, coeffs, prefactor=1, idx=None, steps=None):
        if (
                len(coeffs) > 0 and
                (nput.is_numeric(coeffs) or nput.is_numeric(coeffs[0]))
        ):
            raise ValueError("coeffs must be a vector of vectors (not {})".format(coeffs))
        self.coeffs = [np.asanyarray(c) for c in coeffs] # coeffs along each dim independently
        if any(c.dtype == np.dtype(object) for c in self.coeffs):
            raise ValueError(self.coeffs)
        self.prefactor = prefactor
        self._order = None
        self._monics = None
        self._idx = idx
        self.steps = steps
        if self._idx is not None:
            self._cache[self._idx] = self

    _cache = {}
    @classmethod
    def lookup(cls, idx):
        return cls._cache.get(idx, None)

    def mutate(self, coeffs:'list[np.ndarray]'=default, prefactor:'Any'=default, idx=None, steps=default):
        return type(self)(
            self.coeffs if coeffs is default else coeffs,
            prefactor=self.prefactor if prefactor is default else prefactor,
            idx=idx,
            steps=self.steps if steps is default else steps
        )

    def to_state(self, serializer=None): # don't let this be directly serialized
        raise NotImplementedError(...)
    #     """
    #     Provides just the state that is needed to
    #     serialize the object
    #     :param serializer:
    #     :type serializer:
    #     :return:
    #     :rtype:
    #     """
    #     return {
    #         'coeffs': [c.tolist() for c in self.coeffs],
    #         'prefactor': self.prefactor,
    #         'steps': self.steps,
    #         # 'key': self._idx
    #     }
    # @classmethod
    # def from_state(cls, state, serializer=None):
    #     coeffs = state['coeffs']
    #     prefactor = state['prefactor']
    #     steps = state['steps']
    #     # key = state['key']
    #     # if key is not None:
    #     #     test = cls.lookup(key)
    #     # else:
    #     #     test = None
    #     # if test is None:
    #     return cls(
    #         coeffs, prefactor=prefactor, steps=steps
    #     )
    #     # else:
    #     #     return test

    @property
    def ndim(self):
        # if len(self.coeffs) == 1 and len(self.coeffs[0]) == 1:
        #     return 0
        # else:
        return len(self.coeffs)
    @property
    def order(self):
        if self._order is None:
            self._order = tuple(len(c)-1 for c in self.coeffs)
        return self._order
    def audit(self, target, ignore_constants=True):
        if len(self.coeffs) != target:
            if ignore_constants and not (
                    len(self.coeffs) == 1 and
                    len(self.coeffs[0]) == 1 and
                    target == 0
            ):  # ignore constants
                raise ValueError("{} has wrong dimension (expected {})".format(
                    self,
                    target
                ))
    def ensure_dimension(self, ndim):
        return self.pad(ndim - self.ndim) if self.ndim < ndim else self
    def pad(self, left_right_pads):
        if nput.is_numeric(left_right_pads):
            left_right_pads = [0, left_right_pads]
        l,r = left_right_pads
        new_coeffs = [
            np.array([1]) for _ in range(l)
        ] + self.coeffs + [
            np.array([1]) for _ in range(r)
        ]

        return self.mutate(new_coeffs)

    def __repr__(self):
        return "{}(<{}>)".format("Poly", ",".join(str(c) for c in self.order))

    @classmethod
    def prep_float(cls, c):
        if isinstance(c, int):
            return c
        elif c.is_integer():
            return int(c)
        else:
            return round(c, 4)
    @classmethod
    def format_simple_poly(cls, coeffs, idx):
        var = "n[{}]".format(idx)
        return "+".join(k for k in [
            "{}{}".format(
                (
                    "" if c == 1 else
                    "-" if c == -1 else
                    cls.prep_float(c)
                ) if i != 0 else cls.prep_float(c),
                "" if i == 0 else
                var if i == 1 else
                var + "^{}".format(i)
            ) for i, c in enumerate(coeffs)
            if c != 0
        ] if len(k) > 0)
    def format_expr(self):
        subpolys = [
            self.format_simple_poly(c, i)
            for i,c in enumerate(self.coeffs)
        ]
        subpolys = ["({})".format(s) for s in subpolys if len(s) > 0]
        c = self.prefactor
        poly_str = "".join(s for s in subpolys if len(s) > 0)
        if poly_str == "":
            return "" if c == 1 else str(c)
        else:
            prefac_str = "{}/s[{}]*".format(self.prep_float(c), self.steps)
            return prefac_str + poly_str

    def permute(self, new_inds, check_perm=True, allow_padding=False):
        # must be a proper permutation
        coeffs = self.coeffs
        if allow_padding:
            if len(coeffs) < len(new_inds):
                coeffs = coeffs + [
                    np.array([1.])
                    for _ in range(len(new_inds) - len(coeffs))
                ]
        if check_perm:
            if np.any(
                    np.sort(
                        np.concatenate([new_inds, np.arange(len(new_inds), len(coeffs))])
                    ) != np.arange(len(coeffs))
            ):
                raise ValueError(f"bad permutation {new_inds} from into {len(coeffs)} dimensions")
        return self.mutate(
            [coeffs[i] for i in new_inds] + coeffs[len(new_inds):]
        )

    def constant_rescale(self):
        """
        rescales so constant term is 1

        :return:
        """
        new_prefactor = self.prefactor
        new_coeffs = []
        for c in self.coeffs:
            if abs(c[0]) > 1e-8: # non-zero
                new_prefactor = self.prefactor * c[0]
                new_coeffs.append(c/c[0])
            else:
                new_coeffs.append(c)
        return self.mutate(new_coeffs, prefactor=new_prefactor)

    @property
    def monic_coeffs(self):
        if self._monics is None:
            self._monics = tuple(self._monify(c) for c in self.coeffs)
        return self._monics
    @monic_coeffs.setter
    def monic_coeffs(self, coeffs):
        self._monics = coeffs
    @staticmethod
    def _monify(coeffs, zero_thresh=1e-8):
        # scaling = abs(coeffs[-1])
        # if scaling <= zero_thresh:
        #     return np.zeros_like(coeffs), 0
        # elif scaling == 1:
        #     return coeffs, 1
        # else:
        #     return coeffs / scaling, scaling

        acoeffs = np.abs(coeffs) # np.round(np.abs(coeffs), 8)
        acmask = np.where(acoeffs > zero_thresh)
        if len(acmask) == 0 or len(acmask[0]) == 0:
            return np.zeros_like(coeffs), 0

        acmask = acmask[0]
        min_v = np.min(acoeffs[acmask])
        new_coeffs = np.zeros_like(coeffs)
        new_coeffs[acmask] = coeffs[acmask] / min_v
        return new_coeffs, min_v

    @staticmethod
    def _find_off_pos(self_coeffs, other_coeffs):
        off_pos = None
        for n,(c1,c2) in enumerate(zip(self_coeffs, other_coeffs)):
            if len(c1) != len(c2):
                if off_pos is not None:
                    return False
                off_pos = n
        return off_pos
    def combine(self, other:'ProductPTPolynomial'):
        # if other in self._checked:
        #     ...
        # self._checked.add(other)

        if any(ms < 1e-8 for mc,ms in self.monic_coeffs):
            if any(ms < 1e-8 for mc,ms in other.monic_coeffs): return True
            return other
        elif any(ms < 1e-8 for mc,ms in other.monic_coeffs):
            return self

        off_pos = self._find_off_pos(self.coeffs, other.coeffs) # if a _single_ mode differs we can still combine
        # first pass off lengths b.c. most can't be simplified...
        if off_pos is False: return off_pos

        # second pass where we actually compare the coeffs
        condensed_prefactors = [self.prefactor, other.prefactor]
        new_coeffs = []
        for n,(c1,c2,m1,m2) in enumerate(zip(self.coeffs, other.coeffs, self.monic_coeffs, other.monic_coeffs)):

            if off_pos is not None and off_pos == n: # can handle _one_ difference
                if len(c2) > len(c1):
                    c1 = np.pad(c1, [0, len(c2)-len(c1)])
                if len(c1) > len(c2):
                    c2 = np.pad(c2, [0, len(c1)-len(c2)])
                new_coeffs.append([c1, c2])
                continue


            monic_coeffs1, monic_scaling1 = m1
            monic_coeffs2, monic_scaling2 = m2
            if monic_scaling1 < 1e-8:
                if monic_scaling2 < 1e-8: return True
                new_coeffs.append(c2)
            elif monic_scaling2 < 1e-8:
                new_coeffs.append(c1)
            elif np.any(monic_coeffs1 != monic_coeffs2):
                    if off_pos is None:  # can handle _one_ difference
                        off_pos = n
                        new_coeffs.append([c1, c2])
                        continue
                    else:
                        return False
            else:
                condensed_prefactors[0] *= monic_scaling1
                condensed_prefactors[1] *= monic_scaling2

                # scaling = new_prefactor.prefactor*monic_scaling1 + other.prefactor*monic_scaling2
                # if scaling == 0:
                #     return True # special case... where they cancel perfectly

                new_coeffs.append(monic_coeffs1)

        # need to adjust 1/sqrt(2)**n prefactors
        if self.steps == other.steps:
            new_steps = self.steps
        elif self.steps < other.steps:
            new_steps = other.steps
            step_diff = other.steps - self.steps
            condensed_prefactors[0] *= np.sqrt(2) ** step_diff
        else:  # other.steps < self.steps
            new_steps = self.steps
            step_diff = self.steps - other.steps
            condensed_prefactors[1] *= np.sqrt(2) ** step_diff
        if off_pos is None: # no positions where coeffs differed
            new_prefactor = sum(condensed_prefactors)
            if new_prefactor == 0:
                return True
            new = ProductPTPolynomial(new_coeffs, prefactor=new_prefactor, steps=new_steps)
            new.monic_coeffs = [(c,1) for c in new_coeffs]
        else:
            # gotta condense the position where things differed
            new_c = condensed_prefactors[0]*new_coeffs[off_pos][0] + condensed_prefactors[1]*new_coeffs[off_pos][1]
            monic_coeffs, monic_scaling = self._monify(new_c)
            if monic_scaling == 0:
                return True
            new_coeffs[off_pos] = monic_coeffs
            new = ProductPTPolynomial(new_coeffs, prefactor=monic_scaling, steps=new_steps)
            new_monics = [
                m if k != off_pos else (monic_coeffs, 1)
                for k,m in enumerate(self.monic_coeffs)
            ]
            new.monic_coeffs = new_monics
        return new

    def condense(self, inds=None, return_inds=False, check_inds=True):
        if inds is None: inds = np.arange(len(self.coeffs))
        if check_inds:
            condensed = np.array([i for i in inds if len(self.coeffs[i]) == 1])
        else:
            condensed = inds
        if len(condensed) > 0:
            new = self.mutate(
                [c for c in self.coeffs if c not in condensed]
            )
        else:
            new = self
        if return_inds:
            return condensed, new
        else:
            return new

    def shift(self, shift):
        if len(shift) < len(self.coeffs):
            shift = list(shift) + [0] * (len(self.coeffs) - len(shift))
        if self._idx is not None:
            base_shift = self._idx[-1]
            new_shift = tuple(s + k for s,k in zip(base_shift, shift)) + base_shift[len(shift):]
            new_idx = self._idx[:-1] + (new_shift,)
        else:
            new_idx = None
        return self.mutate(
            [DensePolynomial.compute_shifted_coeffs(c, s) for c,s in zip(self.coeffs, shift)],
            idx=new_idx
        )

    def __mul__(self, other):
        if nput.is_numeric(other):
            return self.scale(other)
        else:
            raise NotImplementedError("subtle")
    def __rmul__(self, other):
        return self.__mul__(other)
    def scale(self, scalar):
        if nput.is_numeric(scalar):
            if scalar == 1:
                return self
            elif scalar == 0:
                return 0
        return self.mutate(prefactor=self.prefactor*scalar)

    def mul_simple(self, other:'PolynomialInterface'):
        if nput.is_numeric(other):
            if other == 1:
                return self
            elif other == 0:
                return 0
            else:
                raise ValueError(other)
        elif isinstance(other, ProductPTPolynomial):
            return self._poly_mul(self, other)
        else:
            return other.rmul_simple(self)
    def rmul_simple(self, other: 'PolynomialInterface') -> 'PolynomialInterface':
        return other.mul_simple(self)

    _prod_poly_cache = {}
    @classmethod
    def _poly_mul(cls, self, other):

        cacheable = other._idx is not None and self._idx is not None
        key = ((self._idx, other._idx), (0,) * len(self.coeffs)) if cacheable else None
        new = self._prod_poly_cache.get(key, None) if cacheable else None

        if new is None:
            ocs = other.coeffs
            scs = self.coeffs
            if len(ocs) < len(scs):
                ocs = ocs + [[1]]*(len(scs) - len(ocs))
            elif len(scs) < len(ocs):
                scs = scs + [[1]]*(len(ocs) - len(scs))

                # raise ValueError("not sure how to 'simply multiply' {} and {}".format(self, other))
            new = cls(
                [
                    scipy.signal.convolve(sc, oc)
                    for sc, oc in zip(scs, ocs)
                ],
                prefactor=self.prefactor * other.prefactor,
                idx=key,
                steps=self.steps + other.steps
            )
            if cacheable: cls._prod_poly_cache[key] = new

        return new

    @classmethod
    def fast_ind_remainder(cls, n, diff):
        chunks = []
        s = 0
        for d in np.sort(diff):
            if d >= n:
                chunks.append(np.arange(s, n))
                break
            if s != d:
                chunks.append(np.arange(s, d))
            s = d + 1
        else:
            if s < n:
                chunks.append(np.arange(s, n))
        if len(chunks) > 1:
            chunks = np.concatenate(chunks)
        elif len(chunks) == 1:
            chunks = chunks[0]
        else:
            chunks = np.array([])
        return chunks

    @classmethod
    def get_index_mapping(cls, dim_1, dim_2, inds, return_remainder=True):
        """
        Computes the corresponding and indices remainders for multdimensional polynomial multiplications

        :param dim_1:
        :param dim_2:
        :param inds:
        :return:
        """
        dim_1 = max([dim_1, len(inds[0])])
        dim_2 = max([dim_2, len(inds[1])])
        remainder = None
        if nput.is_numeric(inds):
            # simplest form, multiply the first indices
            k = inds
            inds = [
                list(range(k)),
                list(range(k))
            ]
            if return_remainder:
                remainder = [
                    np.arange(k, dim_1),
                    np.arange(k, dim_2)
                ]
        elif len(inds) == 0:
            # k = inds
            inds = [
                [],
                []
            ]
            if return_remainder:
                remainder = [
                    np.arange(dim_1),
                    np.arange(dim_2)
                ]
        elif nput.is_numeric(inds[0]):
            # only the 2nd poly inds are specified
            k = len(inds)
            inds = [
                list(range(k)),
                inds
            ]
            if return_remainder:
                remainder = [
                    np.arange(k, dim_1),
                    cls.fast_ind_remainder(dim_2, inds[1]),
                ]
        else:
            if return_remainder:
                remainder = [
                    cls.fast_ind_remainder(dim_1, inds[0]),
                    cls.fast_ind_remainder(dim_2, inds[1])
                ]

        # subinds = [
        #     n for n,(i1, i2) in enumerate(zip(*inds))
        #     if i1 < dim_1 and i2 < dim_2
        # ]

        new_rems = [
            [],
            [],
        ]
        new_inds = [
            [],
            []
        ]
        for i1, i2 in zip(*inds):
            if i1 < dim_1 and i2 < dim_2:
                new_inds[0].append(i1)
                new_inds[1].append(i2)
            else:
                new_rems[0].append(i1)
                new_rems[1].append(i2)


        inds = new_inds
        if return_remainder:
            remainder = [
                list(remainder[0]) + new_rems[0],
                list(remainder[1]) + new_rems[1]
            ]

        res = inds
        if return_remainder:
            res = (inds, remainder)
        return res

    def mul_along(self, other:'PolynomialInterface', inds, remainder=None, mapping=None):
        if isinstance(other, ProductPTPolynomial):
            # we take the polynomial product of indices that align
            # and direct product the remainder
            if len(inds) == 0: # pure product
                # we get to just concatenate the coeffs
                new = self.mutate(
                    #TODO: think about whether mutate is the right concept here...
                    self.coeffs+other.coeffs,
                    prefactor=self.prefactor*other.prefactor,
                    steps=self.steps + other.steps
                )
            else:
                # if remainder is None:
                (self_inds, other_inds), (self_remainder, other_remainder) = self.get_index_mapping(
                    self.ndim, other.ndim, inds, return_remainder=True
                    # len(self.coeffs), len(other.coeffs), inds, return_remainder=True
                )
                # else:
                #     (self_inds, other_inds) = self.get_index_mapping(
                #         self.ndim, other.ndim, inds, return_remainder=False
                #     )
                #     (self_remainder, other_remainder) = remainder

                # print(
                #     other_inds,
                #     self_inds,
                #     self_remainder,
                #     other_remainder
                # )
                # print(self.format_expr())
                # print(other.format_expr())

                new_coeffs = [
                                 scipy.signal.convolve(
                                     other.coeffs[ox] if len(other.coeffs) > ox else [1],
                                     self.coeffs[sx] if len(self.coeffs) > sx else [1]
                                 )
                                 for ox, sx in zip(other_inds, self_inds)
                                 # if len(other.coeffs) > ox and len(self.coeffs) > sx
                             ] + [
                                 self.coeffs[sx] if len(self.coeffs) > sx else [1]
                                 for sx in self_remainder
                                 # if len(self.coeffs) > sx
                             ] + [
                                 other.coeffs[ox] if len(other.coeffs) > ox else [1]
                                 for ox in other_remainder
                                 # if len(other.coeffs) > ox
                             ]
                if len(new_coeffs) == 0:
                    new_coeffs = [[1.]]

                new = self.mutate(
                    new_coeffs,
                    prefactor=self.prefactor*other.prefactor,
                    steps=self.steps + other.steps
                )
        else:
            new = other.rmul_along(self, inds, remainder=remainder, mapping=mapping)

        return new

    def rmul_along(self, other, inds, remainder=None, mapping=None):
        """
        We multiply every subpoly along the given indices, transposing the appropriate tensor indices

        :param other:
        :param inds:
        :param remainder:
        :return:
        """
        return other.mul_along(self, inds, remainder=remainder, mapping=mapping)

    def __add__(self, other):
        if isinstance(other, PTEnergyChangeProductSum):
            return other + self
        if isinstance(other, (int, float, np.integer, np.floating)):
            if other == 0:
                return self
            else:
                raise NotImplementedError(...)
                return self + ProductPTPolynomial([[other]])
        elif isinstance(other, ProductPTPolynomial):
            return ProductPTPolynomialSum([self, other])#, prefactor=self.prefactor)
        elif isinstance(other, ProductPTPolynomialSum):
            return other + self
        else:
            raise NotImplementedError("not sure how to add {} and {}".format(self, other))
    def __radd__(self, other):
        return self + other

class ProductPTPolynomialSum(PolynomialInterface):

    def __init__(self, polynomials, prefactor=1, ndim=None, order=None, reduced=False):
        self.polys = polynomials
        self.prefactor = prefactor
        self.reduced = reduced
        self._order = order
        self._ndim = ndim

    def prep_serialization_dict(self):
        ndim = len(self.polys[0].coeffs)
        return {
            'shared_prefactor': self.prefactor,
            'steps': self.polys[0].steps,
            'ndim': ndim,
            'prefactors': np.array([p.prefactor for p in self.polys]),
            'shapes': (
                np.concatenate([[len(c) for c in p.coeffs] for p in self.polys], axis=0)
                    if ndim > 0 else
                np.array([], dtype=int)
            ),
            'coeffs': (
                np.concatenate([np.concatenate(p.coeffs) for p in self.polys])
                    if ndim > 0 else
                np.array([], dtype=float)
            )
        }
    @classmethod
    def from_serialization_dict(cls, big_dict):
        polys = []
        coeff_vector = big_dict['coeffs']
        ndim = big_dict['ndim']
        shape_vecs = big_dict['shapes']
        steps = big_dict['steps']
        prefactors = big_dict['prefactors']
        if ndim == 0:
            return cls(
                [
                    ProductPTPolynomial([], prefactor=prefactor, steps=steps)
                    for prefactor in prefactors
                ],
                prefactor=big_dict['shared_prefactor']
            )
        shapes = shape_vecs[shape_vecs > 0].reshape(-1, ndim)
        padding = 0
        for i,shape in enumerate(shapes):
            coeffs = []
            for s in shape:
                if s > 0:
                    coeffs.append(coeff_vector[padding:padding+s])
                padding += s
            polys.append(
                ProductPTPolynomial(coeffs, prefactor=prefactors[i], steps=steps)
            )
        return cls(polys, prefactor=big_dict['shared_prefactor'])

        # {
        #     'prefactor': self.prefactor,
        #     'steps': self.polys[0].steps,
        #     'prefactors': np.array([p.prefactor for p in self.polys]),
        #     'shapes': np.array([[len(c) for c in p.coeffs] for p in self.polys]),
        #     'coeffs': np.concatenate([np.concatenate(p.coeffs) for p in self.polys])
        # }

    def to_state(self, serializer=None): # don't let this be directly serialized
        raise NotImplementedError(...)

    def mutate(self, polynomials:'list[PolynomialInterface]'=default,
               prefactor:'Any'=default,
               ndim:'Any'=default,
               order:'Any'=default,
               reduced:'Any'=default
               ):
        return type(self)(
            self.polys if polynomials is default else polynomials,
            prefactor=self.prefactor if prefactor is default else prefactor,
            ndim=self._ndim if ndim is default else ndim,
            order=self._order if order is default else order,
            reduced=self.reduced if reduced is default else reduced
        )

    @property
    def ndim(self):
        if self._ndim is None:
            self.audit()
            self._ndim = self.polys[0].ndim
        return self._ndim
    @property
    def order(self):
        if self._order is None:
            suborders = [p.order for p in self.polys]
            max_len = max(len(o) for o in suborders)
            self._order = tuple(
                max(o[i] if i < len(o) else 0 for o in suborders)
                for i in range(max_len)
            )
        return self._order
    def audit(self, target=None, ignore_constants=True):
        if target is None: target = self.polys[0].ndim
        for p in self.polys: p.audit(target, ignore_constants=ignore_constants)
    def ensure_dimension(self, ndim):
        new_ndim = self._ndim
        if new_ndim is not None:
            if new_ndim >= ndim:
                return self
            else:
                new_ndim = ndim
        new_polys = [
            k.pad(ndim - k.ndim) if k.ndim < ndim else k
            for k in self.polys
        ]
        return self.mutate(
            new_polys,
            ndim=new_ndim
        )

    def __repr__(self):
        form_counts = {}
        for p in self.polys:
            key = "<{}>".format(",".join(str(c) for c in p.order))
            form_counts[key] = form_counts.get(key, 0) + 1
        return "PSum({})".format(
            # type(self).__name__,
            "+".join("{}{}".format(k,f) for k,f in form_counts.items())
        )

    def format_expr(self):
        base_sum = "+".join(p.format_expr() for p in self.polys)
        if self.prefactor != 1:
            if "+" in base_sum: base_sum = "({})".format(base_sum)
            if self.prefactor == -1:
                base_sum = "-"+base_sum
            else:
                base_sum = "{}*{}".format(self.prefactor, base_sum)
        return base_sum

    def permute(self, new_inds, check_perm=True, allow_padding=False):
        return self.mutate(
            [p.permute(new_inds, check_perm=check_perm, allow_padding=allow_padding) for p in self.polys],
            prefactor=self.prefactor,
            reduced=self.reduced
        )

    @classmethod
    def combine_polys(cls, poly_set, cache):
        new_polys = []
        eliminated_pos = set()
        for i,p1 in enumerate(poly_set):
            if i not in eliminated_pos:
                for j,p2 in enumerate(poly_set[i+1:]):
                    j = i+j+1
                    if j not in eliminated_pos and (p1, p2) not in cache:
                        # so we can avoid trying to recombine over passes
                        newp = p1.combine(p2)
                        if newp is not False: # combined
                            eliminated_pos.add(i)
                            eliminated_pos.add(j)
                            if newp is not True: #but didn't cancel
                                new_polys.append(newp)
                            break
                        else:
                            cache.add((p1, p2))
                else:
                    new_polys.append(p1) # made it through untouched

        return new_polys
    def combine(self):
        if self.reduced: return self
        polys = self.polys
        cur_len = len(polys) + 1
        cache = set()
        while cur_len > len(polys): # hit a fixed point with this transformation
            cur_len = len(polys)
            polys = self.combine_polys(polys, cache)
        return self.mutate(
            polys,
            prefactor=self.prefactor,
            reduced=True
        )

    def condense(self, inds=None, return_inds=False, check_inds=True):
        if inds is None: inds = np.arange(self.ndim)
        if check_inds:
            o = self.order
            condensed = np.array([i for i in inds if o[i] == 0])
        else:
            condensed = inds
        if len(condensed) > 0:
            new = self.mutate(
                [p.condense for p in self.polys]
            )
        else:
            new = self
        if return_inds:
            return condensed, new
        else:
            return new

    def shift(self, shift):
        return self.mutate(
            [p.shift(shift) for p in self.polys]
        )

    def __mul__(self, other):
        if isinstance(other, (int, float, np.integer, np.floating)):
            return self.scale(other)
        else:
            raise NotImplementedError("subtle")
    def __rmul__(self, other):
        return self.__mul__(other)
    def scale(self, scalar):
        if isinstance(scalar, (int, float, np.integer, np.floating)):
            if scalar == 1:
                return self
            elif scalar == 0:
                return 0
        return self.mutate(prefactor=self.prefactor*scalar)

    def mul_simple(self, other:'ProductPTPolynomial'):
        if isinstance(other, ProductPTPolynomial):
            return self.mutate(
                [
                    poly.mul_simple(other)
                    for poly in self.polys
                ]
            )
        elif isinstance(other, ProductPTPolynomialSum):
            return self.mutate(
                [
                    p1.mul_simple(p2)
                    for p1 in self.polys
                    for p2 in other.polys
                ],
                prefactor=self.prefactor*other.prefactor
            )
        elif nput.is_numeric(other):
            if other == 1:
                return self
            elif other == 0:
                return 0
            else:
                raise ValueError(other)
        else:
            return other.rmul_simple(self)
    def rmul_simple(self, other: 'PolynomialInterface') -> 'PolynomialInterface':
        if isinstance(other, ProductPTPolynomial):
            return self.mutate(
                [
                    other.mul_simple(poly)
                    for poly in self.polys
                ]
            )
        elif nput.is_numeric(other):
            if other == 1:
                return self
            elif other == 0:
                return 0
            else:
                raise ValueError(other)
        else:
            return other.mul_simple(self)

    def mul_along(self, other:'PolynomialInterface', inds, remainder=None, mapping=None):
        if isinstance(other, ProductPTPolynomial):
            # can just distribute except at some level I have to assume
            # that every stored poly has the same dimension...?
            new = self.mutate([
                poly.mul_along(other, inds, remainder=remainder, mapping=mapping)  # how is this supposed to work?
                for poly in self.polys
            ])
        elif isinstance(other, ProductPTPolynomialSum):
            new = self.mutate(
                [
                    p1.mul_along(p2, inds, remainder=remainder, mapping=mapping)
                    for p1 in self.polys
                    for p2 in other.polys
                ],
                prefactor=self.prefactor*other.prefactor
            )
        else:
            new = other.rmul_along(self, inds, remainder=remainder, mapping=mapping)
        return new
    def rmul_along(self, other:'PolynomialInterface', inds, remainder=None, mapping=None):
        """
        We multiply every subpoly along the given indices, transposing the appropriate tensor indices

        :param other:
        :param inds:
        :param remainder:
        :return:
        """
        if isinstance(other, ProductPTPolynomial):
            return self.mutate([
                other.mul_along(poly, inds, remainder=remainder, mapping=mapping)  # how is this supposed to work?
                for poly in self.polys
            ])
        else:
            return other.mul_along(self, inds, remainder=remainder, mapping=mapping)

    def __add__(self, other):
        if isinstance(other, PTEnergyChangeProductSum):
            return other + self
        if isinstance(other, (int, float, np.integer, np.floating)):
            if other == 0:
                return self
            else:
                raise NotImplementedError(...)
                return self + ProductPTPolynomial([[other]], steps=0)
        else:
            if isinstance(other, ProductPTPolynomialSum):
                ops = other.polys
                if other.prefactor != self.prefactor:
                    ops = [p.scale(other.prefactor/self.prefactor) for p in ops]
                return self.mutate(self.polys + ops, reduced=False)
            elif isinstance(other, ProductPTPolynomial):
                if self.prefactor != 1:
                    other = other.scale(1/self.prefactor)
                    # raise NotImplementedError("need to handle prefactors")
                return self.mutate(self.polys + [other], reduced=False)
            else:
                raise TypeError("not sure what to do with {} and {}".format(type(self), type(other)))
    def __radd__(self, other):
        return self + other


class PolyPath(ProductPTPolynomialSum):
    """
    Canonical DAG-backed representation of a sum of separable polynomial products.

    Each key is an interned :class:`PolyTerm`; the corresponding value is its
    scalar coefficient.  Multiplication joins interned axis paths rather than
    convolving coefficient arrays.  ``polys`` and ``to_eager`` provide a
    compatibility boundary for the legacy implementation.
    """

    _comparison_context = contextvars.ContextVar('PolyPath.comparisons', default=None)
    _index_reduction_context = contextvars.ContextVar('PolyPath.index_reduction', default=False)

    @classmethod
    @contextlib.contextmanager
    def cached_comparisons(cls, max_items=32768):
        """Reuse structural equality of immutable paths during one evaluation.

        Identity keys avoid recursively comparing cache keys. Weak references
        prevent the cache from retaining polynomial graphs or accepting a
        recycled object id. A zero limit provides the native A/B path.
        """
        if max_items < 0:
            raise ValueError('comparison cache size must be nonnegative')
        token = cls._comparison_context.set((collections.OrderedDict(), int(max_items)))
        try:
            yield
        finally:
            cls._comparison_context.reset(token)

    _construction_count = 0
    _node_cache = weakref.WeakValueDictionary()
    _node_requests = 0
    _node_hits = 0

    @staticmethod
    def _sort_key(term):
        return (
            term.steps,
            tuple(tuple(atom._key for atom in axis.atoms) for axis in term.axes)
        )

    def __init__(self, terms, reduced=False, node=None, ndim=None, order=None):
        if isinstance(terms, dict):
            terms = terms.items()
        merged = {}
        for term, scaling in terms:
            if nput.is_zero(scaling):
                continue
            merged[term] = merged.get(term, 0) + scaling
        self._items = tuple(
            (term, scaling)
            for term, scaling in sorted(merged.items(), key=lambda x: self._sort_key(x[0]))
            if not nput.is_zero(scaling)
        )
        self.prefactor = 1
        self.reduced = reduced
        self._node = node
        self._ndim = ndim
        self._order = order
        self._hash = None
        self._combined = None
        if Diagnostics.enabled:
            type(self)._construction_count += 1

    def __fingerprint__(self):
        # Hashes, materialized axes and comparison caches are derived state.
        return (type(self).__module__,type(self).__qualname__,
                self._node,self._items,self.ndim,self.reduced)

    @classmethod
    def _from_node(cls, kind, args, ndim, order):
        if Diagnostics.enabled:
            cls._node_requests += 1
        if kind=='remap' and cls._index_reduction_context.get():
            child,pull,output_ndim=args
            # Compose pure index views. No coefficients are expanded or added.
            while child._node is not None and child._node[0]=='remap':
                grandchild,inner,_=child._node[1]
                pull=tuple(pull[i] for i in inner)
                child=grandchild
            if output_ndim==child.ndim and pull==tuple(range(child.ndim)):
                return child
            # A one-term leaf can be relabeled exactly using native axes.
            # Keeping larger leaves lazy avoids materializing polynomial sums.
            if child._node is None and len(child._items)==1:
                term,scaling=child._items[0]
                if len(term.axes)==child.ndim:
                    identity=PolyAxis((PolyAtom([1]),))
                    axes=[identity]*output_ndim
                    for old,new in enumerate(pull):axes[new]=term.axes[old]
                    return cls({PolyTerm(tuple(axes),steps=term.steps):scaling},reduced=True)
            args=(child,pull,output_ndim)
        # Addition is intentionally binary and remaps/symmetry orbits have very
        # high cardinality.  Interning those nodes costs substantially more
        # memory than their observed (<1%) hit rate at fourth order; their
        # children and canonical leaves remain interned and structurally
        # hashable.
        if kind in {'add', 'remap', 'permutation_sum'}:
            return cls((), reduced=True, node=(kind, args), ndim=ndim, order=tuple(order))
        key = (kind, args)
        path = cls._node_cache.get(key)
        if path is None:
            path = cls((), reduced=True, node=key, ndim=ndim, order=tuple(order))
            cls._node_cache[key] = path
        else:
            if Diagnostics.enabled:
                cls._node_hits += 1
        return path

    @property
    def is_zero(self):
        return self._node is None and len(self._items) == 0

    @classmethod
    def from_coeffs(cls, coeffs, prefactor=1, idx=None, steps=None):
        axes = tuple(PolyAxis((PolyAtom(c),)) for c in coeffs)
        return cls({PolyTerm(axes, steps=steps): prefactor}, reduced=True)

    @classmethod
    def from_polynomial(cls, poly):
        if isinstance(poly, cls):
            return poly
        if isinstance(poly, ProductPTPolynomial):
            return cls.from_coeffs(
                poly.coeffs,
                prefactor=poly.prefactor,
                idx=poly._idx,
                steps=poly.steps
            )
        if isinstance(poly, ProductPTPolynomialSum):
            terms = cls((), reduced=True)
            for subpoly in poly.polys:
                terms = terms + cls.from_polynomial(subpoly).scale(poly.prefactor)
            return terms
        raise TypeError("can't convert {} to PolyPath".format(type(poly)))

    @classmethod
    def clear_caches(cls):
        PolyAtom._cache.clear()
        PolyAxis._cache.clear()
        PolyTerm._cache.clear()
        cls._node_cache.clear()
        PolyAxis._materialization_count = 0
        cls._construction_count = 0
        cls._node_requests = 0
        cls._node_hits = 0

    @classmethod
    def cache_info(cls):
        return {
            'atoms': len(PolyAtom._cache),
            'axes': len(PolyAxis._cache),
            'terms': len(PolyTerm._cache),
            'nodes': len(cls._node_cache),
            'node_requests': cls._node_requests,
            'node_hits': cls._node_hits,
            'paths_created': cls._construction_count,
            'axis_materializations': PolyAxis._materialization_count
        }

    @property
    def path_terms(self):
        if self._node is not None:
            raise ValueError("DAG-backed PolyPath has no flat term table")
        return dict(self._items)

    @property
    def polys(self):
        eager = self.to_eager()
        if nput.is_numeric(eager):
            return []
        if isinstance(eager, ProductPTPolynomial):
            return [eager]
        return eager.polys

    def coefficient(self, powers, shift=0):
        """Extract [I**powers] P(I+shift) without evaluating occupations.

        Walk coefficient requests through the shared DAG. Integer path shifts,
        index maps and contact products are applied before extracting the
        requested coefficient. Only one-dimensional leaf axes materialize;
        neither a multivariate eager expansion nor state sampling is needed.
        ``shift=-.5`` changes occupation variables to actions I=n+1/2.
        """
        powers = tuple(powers)
        if any(int(p) != p or p < 0 for p in powers):
            raise ValueError('coefficient powers must be nonnegative integers')
        powers = tuple(int(p) for p in powers)
        if any(powers[self.ndim:]):
            return 0.
        powers = (powers + (0,)*self.ndim)[:self.ndim]
        if np.isscalar(shift):
            shift = (float(shift),)*self.ndim
        else:
            shift = tuple(float(s) for s in shift)
            if len(shift) != self.ndim:
                raise ValueError('formal shift dimension differs from polynomial')
        root = (self, powers, shift)
        values = {}
        pending = [(root, None)]

        def mapped(child, pull, degree, offset):
            # Native remaps are injective changes of index labels.
            if len(set(pull)) != len(pull):
                raise ValueError('coefficient extraction requires injective index maps')
            if any(d for i,d in enumerate(degree) if i not in pull):
                return None
            return (child, tuple(degree[i] for i in pull), tuple(offset[i] for i in pull))

        while pending:
            request, terms = pending.pop()
            if request in values:
                continue
            node, degree, offset = request
            if terms is not None:
                values[request] = math.fsum(
                    weight * math.prod(values[dep] for dep in deps)
                    for weight,deps in terms)
                continue
            if node._node is None:
                contributions = []
                for term, scaling in node._items:
                    weight = scaling * 2.**(-term.steps/2)
                    for i,power in enumerate(degree):
                        if i >= len(term.axes):
                            weight *= (power == 0)
                        else:
                            coeffs = term.axes[i].coeffs
                            if offset[i]:
                                coeffs = DensePolynomial.compute_shifted_coeffs(coeffs, offset[i])
                            weight *= coeffs[power] if power < len(coeffs) else 0
                        if weight == 0:
                            break
                    contributions.append(weight)
                values[request] = math.fsum(contributions)
                continue
            kind, args = node._node
            terms = []
            def append(weight, *deps):
                if weight and all(dep is not None for dep in deps):
                    terms.append((weight, deps))
            if kind == 'add':
                for child in args:
                    append(1, (child, degree, offset))
            elif kind == 'linear':
                for child,weight in args:
                    append(weight, (child, degree, offset))
            elif kind == 'scale':
                child,weight = args
                append(weight, (child, degree, offset))
            elif kind == 'shift':
                child,delta = args
                delta = tuple(delta)+(0,)*len(offset)
                append(1, (child, degree, tuple(a+b for a,b in zip(offset,delta))))
            elif kind == 'remap':
                child,pull,_ = args
                append(1, mapped(child,pull,degree,offset))
            elif kind == 'permutation_sum':
                child,_,pulls = args
                for pull in pulls:
                    append(1, mapped(child,pull,degree,offset))
            elif kind == 'project':
                child,keep,_ = args
                child_degree = [0]*child.ndim
                child_offset = [0.]*child.ndim
                for i,j in enumerate(keep):
                    child_degree[j],child_offset[j] = degree[i],offset[i]
                append(1, (child,tuple(child_degree),tuple(child_offset)))
            elif kind == 'mul':
                left,right,_,_,left_pull,right_pull = args
                for left_degree in itertools.product(*(range(p+1) for p in degree)):
                    right_degree = tuple(p-q for p,q in zip(degree,left_degree))
                    append(1, mapped(left,left_pull,left_degree,offset),
                           mapped(right,right_pull,right_degree,offset))
            else:
                raise NotImplementedError('coefficient extraction for PolyPath node '+kind)
            pending.append((request, terms))
            pending.extend((dep,None) for _,deps in terms for dep in deps if dep not in values)
        return values[root]

    def _leaf_to_eager(self):
        return [
            ProductPTPolynomial(
                [axis.coeffs for axis in term.axes],
                prefactor=scaling,
                steps=term.steps
            )
            for term, scaling in self._items
        ]

    def to_eager(self):
        if self._node is None:
            polys = self._leaf_to_eager()
            if len(polys) == 0:
                return 0
            elif len(polys) == 1:
                return polys[0]
            return ProductPTPolynomialSum(polys, reduced=self.reduced)

        kind, args = self._node
        if kind == 'add':
            eager = 0
            for child in args:
                eager = eager + child.to_eager()
            return eager
        elif kind == 'linear':
            eager = 0
            for child, scaling in args:
                contribution = child.to_eager()
                eager = eager + (
                    contribution if scaling == 1 else contribution.scale(scaling)
                )
            return eager
        elif kind == 'scale':
            child, scaling = args
            eager = child.to_eager()
            return eager if scaling == 1 else eager.scale(scaling)
        elif kind == 'shift':
            child, shift = args
            return child.to_eager().shift(shift)
        elif kind == 'remap':
            child, pull, output_ndim = args
            eager = child.to_eager()
            if nput.is_numeric(eager):
                return eager
            polys = eager.polys if isinstance(eager, ProductPTPolynomialSum) else [eager]
            remapped = []
            for poly in polys:
                coeffs = [np.array([1])] * output_ndim
                for source, target in enumerate(pull):
                    if source < len(poly.coeffs):
                        coeffs[target] = poly.coeffs[source]
                remapped.append(ProductPTPolynomial(
                    coeffs, prefactor=poly.prefactor, steps=poly.steps
                ))
            return remapped[0] if len(remapped) == 1 else ProductPTPolynomialSum(remapped)
        elif kind == 'permutation_sum':
            child, permutations, pulls = args
            eager = 0
            child_eager = child.to_eager()
            for permutation in permutations:
                eager = eager + child_eager.permute(permutation)
            return eager
        elif kind == 'project':
            child, keep, condensed = args
            eager = child.to_eager()
            if nput.is_numeric(eager):
                return eager
            polys = eager.polys if isinstance(eager, ProductPTPolynomialSum) else [eager]
            projected = []
            for poly in polys:
                scaling = poly.prefactor
                for index in condensed:
                    scaling *= poly.coeffs[index][0]
                projected.append(ProductPTPolynomial(
                    [poly.coeffs[index] for index in keep],
                    prefactor=scaling,
                    steps=poly.steps
                ))
            return (
                projected[0]
                    if len(projected) == 1 else
                ProductPTPolynomialSum(projected)
            )
        elif kind == 'mul':
            left, right, left_inds, right_inds, left_pull, right_pull = args
            left_eager, right_eager = left.to_eager(), right.to_eager()
            if left_inds is None:
                return left_eager.mul_simple(right_eager)
            return left_eager.mul_along(right_eager, [left_inds, right_inds])
        raise ValueError("unknown PolyPath node {}".format(kind))

    def prep_serialization_dict(self):
        eager = self.to_eager()
        if nput.is_numeric(eager):
            raise ValueError("cannot serialize an empty PolyPath")
        if isinstance(eager, ProductPTPolynomial):
            eager = ProductPTPolynomialSum([eager], reduced=True)
        return eager.prep_serialization_dict()

    def mutate(self, polynomials=default, prefactor=default, **kwargs):
        new = self if polynomials is default else type(self).from_polynomial(
            ProductPTPolynomialSum(polynomials)
                if not isinstance(polynomials, PolynomialInterface) else
            polynomials
        )
        if prefactor is not default:
            new = new.scale(prefactor)
        return new

    @property
    def ndim(self):
        if self._ndim is None:
            self._ndim = max([len(t.axes) for t, _ in self._items] + [0])
        return self._ndim

    @property
    def order(self):
        if self._order is None:
            self._order = tuple(
                max((t.axes[i].order if i < len(t.axes) else 0) for t, _ in self._items)
                for i in range(self.ndim)
            )
        return self._order

    def __hash__(self):
        # Deep polynomial DAGs occur at higher perturbation orders. Resolve
        # child hashes iteratively so cache lookup has no Python depth limit.
        pending = [(self,False)] if self._hash is None else []
        while pending:
            path,ready = pending.pop()
            if path._hash is not None:
                continue
            if ready or path._node is None:
                path._hash = hash(path._items if path._node is None else path._node)
                continue
            pending.append((path,True))
            arguments = [path._node]
            while arguments:
                argument = arguments.pop()
                if isinstance(argument,PolyPath):
                    if argument._hash is None:
                        pending.append((argument,False))
                elif isinstance(argument,tuple):
                    arguments.extend(argument)
        return self._hash

    def __eq__(self, other):
        if self is other:
            return True
        context = self._comparison_context.get()
        if context is None or not context[1] or not isinstance(other, PolyPath):
            return self._structurally_equal(other)
        pairs, limit = context
        key = (id(self), id(other))
        entry = pairs.get(key)
        if entry is not None and entry[0]() is self and entry[1]() is other:
            return entry[2]
        answer = self._structurally_equal(other)
        pairs[key] = (weakref.ref(self), weakref.ref(other), answer)
        if len(pairs) > limit:
            pairs.popitem(last=False)
        return answer

    def _structurally_equal(self, other):
        pending,checked = [(self,other)],set()
        while pending:
            left,right = pending.pop()
            if left is right:
                continue
            if isinstance(left,PolyPath):
                if not isinstance(right,type(left)) or left._items != right._items:
                    return False
                if left._hash is not None and right._hash is not None and left._hash != right._hash:
                    return False
                pair = (id(left),id(right))
                if pair not in checked:
                    checked.add(pair)
                    pending.append((left._node,right._node))
            elif isinstance(left,tuple):
                if not isinstance(right,tuple) or len(left) != len(right):
                    return False
                pending.extend(zip(left,right))
            elif left != right:
                return False
        return True

    def __repr__(self):
        return "PolyPath(<{}>; {} terms)".format(
            ",".join(str(o) for o in self.order),
            len(self._items) if self._node is None else "DAG"
        )

    def format_expr(self):
        eager = self.to_eager()
        return "0" if nput.is_numeric(eager) else eager.format_expr()

    def audit(self, target=None, ignore_constants=True):
        if target is None:
            target = self.ndim
        if self._node is not None:
            if self.ndim != target and not ignore_constants:
                raise ValueError("{} has wrong dimension (expected {})".format(self, target))
            return
        for term, _ in self._items:
            if len(term.axes) != target:
                if not (ignore_constants and all(axis.order == 0 for axis in term.axes)):
                    raise ValueError("{} has wrong dimension (expected {})".format(self, target))

    def ensure_dimension(self, ndim):
        if self.ndim >= ndim:
            return self
        return type(self)._from_node(
            'remap',
            (self, tuple(range(self.ndim)), ndim),
            ndim,
            self.order + (0,) * (ndim - self.ndim)
        )

    def pad(self, left_right_pads):
        if nput.is_numeric(left_right_pads):
            left_right_pads = [0, left_right_pads]
        left, right = left_right_pads
        ndim = left + self.ndim + right
        return type(self)._from_node(
            'remap',
            (self, tuple(range(left, left + self.ndim)), ndim),
            ndim,
            (0,) * left + self.order + (0,) * right
        )

    def permute(self, new_inds, check_perm=True, allow_padding=False):
        new_inds = tuple(int(i) for i in new_inds)
        ndim = self.ndim
        if allow_padding:
            ndim = max(ndim, len(new_inds))
        if check_perm:
            tail = np.arange(len(new_inds), ndim)
            if np.any(np.sort(np.concatenate([new_inds, tail])) != np.arange(ndim)):
                raise ValueError("bad permutation {} into {} dimensions".format(new_inds, ndim))
        source_order = new_inds + tuple(range(len(new_inds), ndim))
        pull = tuple(int(i) for i in np.argsort(source_order)[:self.ndim])
        padded_order = self.order + (0,) * (ndim - self.ndim)
        return type(self)._from_node(
            'remap',
            (self, pull, ndim),
            ndim,
            tuple(padded_order[i] for i in source_order)
        )

    def permutation_sum(self, permutations, check_perm=True, allow_padding=False):
        """Represent a symmetry sum without constructing each remapped child."""
        permutations = tuple(tuple(int(i) for i in perm) for perm in permutations)
        if len(permutations) == 0:
            return type(self)((), reduced=True)
        if len(permutations) == 1:
            return self.permute(
                permutations[0], check_perm=check_perm, allow_padding=allow_padding
            )

        pulls = []
        orders = []
        output_ndim = self.ndim
        for new_inds in permutations:
            ndim = self.ndim
            if allow_padding:
                ndim = max(ndim, len(new_inds))
            if check_perm:
                tail = np.arange(len(new_inds), ndim)
                if np.any(np.sort(np.concatenate([new_inds, tail])) != np.arange(ndim)):
                    raise ValueError("bad permutation {} into {} dimensions".format(new_inds, ndim))
            if ndim != output_ndim:
                raise ValueError("inconsistent permutation dimensions")
            source_order = new_inds + tuple(range(len(new_inds), ndim))
            pulls.append(tuple(int(i) for i in np.argsort(source_order)[:self.ndim]))
            padded_order = self.order + (0,) * (ndim - self.ndim)
            orders.append(tuple(padded_order[i] for i in source_order))

        order = tuple(max(o[i] for o in orders) for i in range(output_ndim))
        return type(self)._from_node(
            'permutation_sum',
            (self, permutations, tuple(pulls)),
            output_ndim,
            order
        )

    def shift(self, shift):
        shift = tuple(shift) + (0,) * max(0, self.ndim - len(shift))
        shift = shift[:self.ndim]
        if not any(shift):
            return self
        return type(self)._from_node('shift', (self, shift), self.ndim, self.order)

    def scale(self, scaling):
        if nput.is_numeric(scaling):
            if scaling == 0:
                return 0
            elif scaling == 1:
                return self
        return type(self)._from_node('scale', (self, scaling), self.ndim, self.order)

    def evaluate_polynomial(self, substates, node_cache=None, axis_cache=None):
        if node_cache is None:
            node_cache = {}
        if axis_cache is None:
            axis_cache = {}
        state_key = (substates.dtype.str, substates.shape, substates.tobytes())
        cache_key = (self, state_key)
        value = node_cache.get(cache_key)
        if value is not None:
            return value

        if self._node is None:
            value = np.zeros(substates.shape[0], dtype=float)
            for term, scaling in self._items:
                factors = []
                for state_values, axis in zip(substates.T, term.axes):
                    axis_key = (
                        axis,
                        state_values.dtype.str,
                        state_values.shape,
                        state_values.tobytes()
                    )
                    axis_value = axis_cache.get(axis_key)
                    if axis_value is None:
                        axis_value = np.polynomial.polynomial.polyval(state_values, axis.coeffs)
                        axis_cache[axis_key] = axis_value
                    factors.append(axis_value)
                poly_factor = (
                    np.prod(factors, axis=0)
                        if len(factors) > 0 else
                    np.ones(substates.shape[0])
                )
                value += scaling * poly_factor / np.sqrt(2) ** term.steps
        else:
            kind, args = self._node
            if kind == 'add':
                value = np.zeros(substates.shape[0], dtype=float)
                stack = list(args)
                while stack:
                    child = stack.pop()
                    if child._node is not None and child._node[0] == 'add':
                        stack.extend(child._node[1])
                    else:
                        value += child.evaluate_polynomial(substates, node_cache, axis_cache)
            elif kind == 'linear':
                value = np.zeros(substates.shape[0], dtype=float)
                for child, scaling in args:
                    value += scaling * child.evaluate_polynomial(
                        substates, node_cache, axis_cache
                    )
            elif kind == 'scale':
                child, scaling = args
                value = scaling * child.evaluate_polynomial(substates, node_cache, axis_cache)
            elif kind == 'shift':
                child, shift = args
                value = child.evaluate_polynomial(
                    substates[:, :child.ndim] + np.asanyarray(shift)[np.newaxis, :],
                    node_cache, axis_cache
                )
            elif kind == 'remap':
                child, pull, output_ndim = args
                value = child.evaluate_polynomial(
                    substates[:, pull], node_cache, axis_cache
                )
            elif kind == 'permutation_sum':
                child, permutations, pulls = args
                value = np.zeros(substates.shape[0], dtype=float)
                for pull in pulls:
                    value += child.evaluate_polynomial(
                        substates[:, pull], node_cache, axis_cache
                    )
            elif kind == 'project':
                child, keep, condensed = args
                child_states = np.zeros(
                    (substates.shape[0], child.ndim), dtype=substates.dtype
                )
                child_states[:, keep] = substates[:, :len(keep)]
                value = child.evaluate_polynomial(
                    child_states, node_cache, axis_cache
                )
            elif kind == 'mul':
                left, right, left_inds, right_inds, left_pull, right_pull = args
                left_value = left.evaluate_polynomial(
                    substates[:, left_pull], node_cache, axis_cache
                )
                right_value = right.evaluate_polynomial(
                    substates[:, right_pull], node_cache, axis_cache
                )
                value = left_value * right_value
            else:
                raise ValueError("unknown PolyPath node {}".format(kind))

        node_cache[cache_key] = value
        return value

    @staticmethod
    def _normalize_axis(axis):
        coeffs, scaling = ProductPTPolynomial._monify(axis.coeffs)
        if scaling == 0:
            return None, 0
        return PolyAxis((PolyAtom(coeffs),)), scaling

    def _reduce_terms(self):
        if len(self._items) < 2:
            return type(self)(self._items, reduced=True)

        max_steps = max(term.steps for term, _ in self._items)
        identity = PolyAxis((PolyAtom([1]),))
        terms = {}
        for term, scaling in self._items:
            scaling = scaling * np.sqrt(2) ** (max_steps - term.steps)
            axes = []
            for axis in term.axes:
                axis, axis_scaling = self._normalize_axis(axis)
                if axis is None:
                    scaling = 0
                    break
                axes.append(axis)
                scaling *= axis_scaling
            if nput.is_zero(scaling):
                continue
            if len(axes) < self.ndim:
                axes.extend([identity] * (self.ndim - len(axes)))
            new_term = PolyTerm(tuple(axes), max_steps)
            terms[new_term] = terms.get(new_term, 0) + scaling
        terms = {term: scaling for term, scaling in terms.items() if not nput.is_zero(scaling)}

        changed = True
        while changed and len(terms) > 1:
            changed = False
            for axis_index in range(self.ndim):
                grouped = {}
                for term, scaling in terms.items():
                    group_key = term.axes[:axis_index] + term.axes[axis_index + 1:]
                    grouped.setdefault(group_key, []).append((term.axes[axis_index], scaling))

                if not any(len(group) > 1 for group in grouped.values()):
                    continue

                new_terms = {}
                for group_key, group in grouped.items():
                    if len(group) == 1:
                        axis, scaling = group[0]
                    else:
                        max_len = max(len(axis.coeffs) for axis, _ in group)
                        coeffs = np.zeros(max_len)
                        for axis, subscaling in group:
                            axis_coeffs = axis.coeffs
                            coeffs[:len(axis_coeffs)] += subscaling * axis_coeffs
                        monic, scaling = ProductPTPolynomial._monify(coeffs)
                        if scaling == 0:
                            continue
                        axis = PolyAxis((PolyAtom(monic),))
                        changed = True

                    axes = group_key[:axis_index] + (axis,) + group_key[axis_index:]
                    new_term = PolyTerm(axes, max_steps)
                    new_terms[new_term] = new_terms.get(new_term, 0) + scaling
                terms = {
                    term: scaling
                    for term, scaling in new_terms.items()
                    if not nput.is_zero(scaling)
                }

        return type(self)(terms, reduced=True)

    def combine(self, *args, **kwargs):
        if self.is_zero:
            return self
        if self._node is not None:
            coefficients = {}
            stack = [(self, 1)]
            while stack:
                child, scaling = stack.pop()
                if child._node is not None and child._node[0] == 'add':
                    stack.extend((term, scaling) for term in child._node[1])
                elif child._node is not None and child._node[0] == 'linear':
                    stack.extend(
                        (term, scaling * subscaling)
                        for term, subscaling in child._node[1]
                    )
                elif child._node is not None and child._node[0] == 'scale':
                    term, subscaling = child._node[1]
                    stack.append((term, scaling * subscaling))
                else:
                    coefficients[child] = coefficients.get(child, 0) + scaling

            terms = tuple(
                (child, scaling)
                for child, scaling in sorted(
                    coefficients.items(), key=lambda item: (hash(item[0]), id(item[0]))
                )
                if not nput.is_zero(scaling)
            )
            if len(terms) == 0:
                return type(self)((), reduced=True)
            if len(terms) == 1:
                child, scaling = terms[0]
                return child if scaling == 1 else child.scale(scaling)
            order = tuple(max(child.order[i] for child, _ in terms) for i in range(self.ndim))
            return type(self)._from_node('linear', terms, self.ndim, order)
        if self.reduced or len(self._items) == 1:
            return self
        if self._combined is None:
            self._combined = self._reduce_terms()
        return self._combined

    def condense(self, inds=None, return_inds=False, check_inds=True):
        if inds is None:
            inds = np.arange(self.ndim)
        if check_inds:
            condensed = np.array([i for i in inds if self.order[i] == 0], dtype=int)
        else:
            condensed = np.asanyarray(inds, dtype=int)
        if len(condensed) == 0:
            new = self
        else:
            keep = tuple(int(i) for i in ProductPTPolynomial.fast_ind_remainder(
                self.ndim, condensed
            ))
            condensed_key = tuple(int(i) for i in condensed)
            new = type(self)._from_node(
                'project',
                (self, keep, condensed_key),
                len(keep),
                tuple(self.order[i] for i in keep)
            )
        return (condensed, new) if return_inds else new

    def constant_rescale(self):
        eager = self.to_eager()
        if isinstance(eager, ProductPTPolynomial):
            return type(self).from_polynomial(eager.constant_rescale())
        return self

    def mul_simple(self, other):
        if nput.is_numeric(other):
            if other == 0:
                return 0
            elif other == 1:
                return self
            raise ValueError(other)
        if not isinstance(other, (ProductPTPolynomial, ProductPTPolynomialSum)):
            return other.rmul_simple(self)
        other = type(self).from_polynomial(other)
        ndim = max(self.ndim, other.ndim)
        left_order = self.order + (0,) * (ndim - self.ndim)
        right_order = other.order + (0,) * (ndim - other.ndim)
        return type(self)._from_node(
            'mul',
            (
                self, other, None, None,
                tuple(range(self.ndim)), tuple(range(other.ndim))
            ),
            ndim,
            tuple(a + b for a, b in zip(left_order, right_order))
        )

    def rmul_simple(self, other):
        if isinstance(other, (ProductPTPolynomial, ProductPTPolynomialSum)):
            return type(self).from_polynomial(other).mul_simple(self)
        return other.mul_simple(self)

    def mul_along(self, other, inds, remainder=None, mapping=None):
        if not isinstance(other, (ProductPTPolynomial, ProductPTPolynomialSum)):
            return other.rmul_along(self, inds, remainder=remainder, mapping=mapping)
        other = type(self).from_polynomial(other)
        (left_inds, right_inds), (left_rem, right_rem) = ProductPTPolynomial.get_index_mapping(
            self.ndim, other.ndim, inds, return_remainder=True
        )
        left_pull = [None] * self.ndim
        right_pull = [None] * other.ndim
        order = []
        for output, (left, right) in enumerate(zip(left_inds, right_inds)):
            if left < self.ndim:
                left_pull[left] = output
            if right < other.ndim:
                right_pull[right] = output
            order.append(
                (self.order[left] if left < self.ndim else 0)
                + (other.order[right] if right < other.ndim else 0)
            )
        offset = len(order)
        for left in left_rem:
            if left < self.ndim:
                left_pull[left] = offset
                order.append(self.order[left])
                offset += 1
        for right in right_rem:
            if right < other.ndim:
                right_pull[right] = offset
                order.append(other.order[right])
                offset += 1
        if any(i is None for i in left_pull) or any(i is None for i in right_pull):
            raise ValueError("incomplete PolyPath multiplication mapping")
        return type(self)._from_node(
            'mul',
            (
                self, other,
                tuple(left_inds), tuple(right_inds),
                tuple(left_pull), tuple(right_pull)
            ),
            len(order),
            tuple(order)
        )

    def rmul_along(self, other, inds, remainder=None, mapping=None):
        if isinstance(other, (ProductPTPolynomial, ProductPTPolynomialSum)):
            return type(self).from_polynomial(other).mul_along(
                self, inds, remainder=remainder, mapping=mapping
            )
        return other.mul_along(self, inds, remainder=remainder, mapping=mapping)

    def __mul__(self, other):
        return self.scale(other) if nput.is_numeric(other) else self.mul_simple(other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __add__(self, other):
        if isinstance(other, PTEnergyChangeProductSum):
            return other + self
        if nput.is_numeric(other):
            if other == 0:
                return self
            raise NotImplementedError("non-zero scalar polynomial addition")
        other = type(self).from_polynomial(other)
        ndim = max(self.ndim, other.ndim)
        left = self.ensure_dimension(ndim)
        right = other.ensure_dimension(ndim)
        if left.is_zero and right.is_zero:
            return type(self)((), reduced=True)
        if left.is_zero:
            return right
        if right.is_zero:
            return left
        order = tuple(max(left.order[i], right.order[i]) for i in range(ndim))
        return type(self)._from_node('add', (left, right), ndim, order)

    def __radd__(self, other):
        return self.__add__(other)

class _PolynomialProgramTable:
    """Mapped numeric call records shared by every polynomial in one part."""
    operations = ('constant', 'leaf', 'sum', 'linear', 'scale', 'shift', 'remap',
                  'permutation_sum', 'project', 'mul')
    fields = ('opcodes', 'dimensions', 'child_offsets', 'children', 'parameter_offsets',
              'parameters', 'terms', 'axes', 'coefficients', 'coefficient_keys', 'scalars', 'signature')

    def __init__(self, *values):
        for position, (name, value) in enumerate(zip(self.fields, values)):
            # ndarray views retain the read-only mapping while avoiding the
            # memmap subclass's per-scalar indexing overhead in the interpreter.
            if position < 8: value = np.asarray(value)
            elif name == 'coefficients': value = tuple(np.asarray(array) for array in value)
            setattr(self, name, value)
        if len(values) != len(self.fields): raise ValueError('invalid polynomial table')
        self._hash = hash(('polynomial-table', self.signature))

    def __hash__(self): return self._hash
    def __eq__(self, other): return type(other) is type(self) and self.signature == other.signature

    @classmethod
    def pack(cls, roots):
        nodes, identifiers, pending = [], {}, list(roots)
        while pending:
            node = pending.pop()
            if id(node) in identifiers: continue
            identifiers[id(node)] = len(nodes)
            nodes.append(node)
            pending.extend(node.instructions[0][1])
        scalars, scalar_ids, coefficients, coefficient_keys, coefficient_ids = [], {}, [], [], {}
        def scalar(value):
            if isinstance(value, np.generic): key = type(value), value.dtype.str, value.tobytes()
            elif type(value) is float: key = float, value.hex()
            elif type(value) is complex: key = complex, value.real.hex(), value.imag.hex()
            else: key = type(value), value
            if key not in scalar_ids:
                scalar_ids[key] = len(scalars); scalars.append(value)
            return scalar_ids[key]
        opcodes, dimensions, child_offsets, children, parameter_offsets, parameters = [], [], [0], [], [0], []
        terms, axes = [], []
        for node in nodes:
            kind, calls, extra, _ = node.instructions[0]
            opcodes.append(cls.operations.index(kind)); dimensions.append(node.ndim)
            children.extend(identifiers[id(child)] for child in calls)
            child_offsets.append(len(children))
            if kind in {'constant', 'scale'}: data = [scalar(extra)]
            elif kind == 'linear': data = [scalar(value) for value in extra]
            elif kind in {'shift', 'remap', 'permutation_sum'}: data = np.asarray(extra).ravel().tolist()
            elif kind == 'mul': data = [index for pull in extra for index in pull]
            elif kind == 'project': data = [extra[0], *extra[1]]
            elif kind == 'leaf':
                first = len(terms)
                for weight, divisor, factors in extra:
                    start = len(axes)
                    for index, offset, coeffs, signature in factors:
                        if signature not in coefficient_ids:
                            coefficient_ids[signature] = len(coefficients)
                            coefficients.append(coeffs); coefficient_keys.append(signature)
                        axes.append((index, offset, coefficient_ids[signature]))
                    terms.append((scalar(weight), scalar(divisor), start, len(axes)))
                data = [first, len(terms)]
            else: data = []
            parameters.extend(data); parameter_offsets.append(len(parameters))
        arrays = [np.asarray(opcodes, dtype=np.int8), np.asarray(dimensions, dtype=np.int64),
                  np.asarray(child_offsets, dtype=np.int64), np.asarray(children, dtype=np.int64),
                  np.asarray(parameter_offsets, dtype=np.int64), np.asarray(parameters, dtype=np.int64),
                  np.asarray(terms, dtype=np.int64).reshape(-1,4), np.asarray(axes, dtype=np.int64).reshape(-1,3)]
        for array in arrays: array.flags.writeable = False
        table = cls(*arrays, tuple(coefficients), tuple(coefficient_keys), tuple(scalars),
                    dev.fingerprint(tuple(node.signature for node in nodes)))
        for root in roots:
            root.table, root.node_index = table, identifiers[id(root)]
            root.instructions, root.graph_program = (), False
        return table

    def evaluate(self, root, states, node_cache, axis_cache):
        def identity(index, data):
            return self, int(index), (data.dtype.str, data.shape, data.tobytes())
        root_key = identity(root, states)
        cached = node_cache.get(root_key)
        if cached is not None: return cached
        frames = [[root, states, root_key, 0, []]]
        while frames:
            index, data, key, cursor, values = frames[-1]
            kind = self.operations[int(self.opcodes[index])]
            children = self.children[self.child_offsets[index]:self.child_offsets[index+1]]
            extra = self.parameters[self.parameter_offsets[index]:self.parameter_offsets[index+1]]
            if cursor < len(children):
                child = int(children[cursor])
                ndim = int(self.dimensions[child])
                if kind == 'shift': argument = data[:, :ndim]+extra[None, :]
                elif kind == 'remap': argument = data[:, extra]
                elif kind == 'permutation_sum': argument = data[:, extra[cursor*ndim:(cursor+1)*ndim]]
                elif kind == 'project':
                    argument = np.zeros((len(data), int(extra[0])), dtype=data.dtype)
                    argument[:, extra[1:]] = data[:, :len(extra)-1]
                elif kind == 'mul':
                    start = 0 if cursor == 0 else int(self.dimensions[int(children[0])])
                    argument = data[:, extra[start:start+ndim]]
                else: argument = data
                frames[-1][3] += 1
                child_key = identity(child, argument)
                cached = node_cache.get(child_key)
                if cached is None: frames.append([child, argument, child_key, 0, []])
                else: values.append(cached)
                continue
            if kind == 'constant': result = np.full(len(data), self.scalars[int(extra[0])])
            elif kind == 'leaf':
                result = np.zeros(len(data))
                for weight, divisor, first, last in self.terms[int(extra[0]):int(extra[1])]:
                    factors = []
                    for column, offset, coefficient in self.axes[first:last]:
                        selected = data[:, column]
                        if offset: selected = selected+offset
                        axis_key = ('compiled-axis', self.coefficient_keys[coefficient],
                                    selected.dtype.str, selected.shape, selected.tobytes())
                        factor = axis_cache.get(axis_key)
                        if factor is None:
                            factor = np.polynomial.polynomial.polyval(selected, self.coefficients[coefficient])
                            axis_cache[axis_key] = factor
                        factors.append(factor)
                    product = np.prod(factors, axis=0) if factors else np.ones(len(data))
                    result += self.scalars[weight]*product/self.scalars[divisor]
            elif kind == 'scale': result = self.scalars[int(extra[0])]*values[0]
            elif kind == 'mul': result = values[0]*values[1]
            elif kind in {'shift', 'remap', 'project'}: result = values[0]
            else:
                result = np.zeros(len(data))
                for position, value in enumerate(values):
                    result += self.scalars[int(extra[position])]*value if kind == 'linear' else value
            node_cache[key] = result
            frames.pop()
            if frames: frames[-1][4].append(result)
        return result


class _CompiledPolynomial(PolyPath):
    """An immutable numeric program consumed by the existing polynomial backend.

    The program has no symbolic DAG to expand or polynomial factors to convolve
    at load/evaluation time. Input occupations are the only late-bound values.
    """
    program_version = 2

    def __init__(self, instructions, ndim, order, signature, zero=False, occupation_shifts=None, graph_program=False, table=None, node_index=None):
        self.instructions = tuple(instructions)
        self._ndim, self._order = int(ndim), tuple(order)
        self.signature, self.zero = signature, bool(zero)
        self.prefactor, self.reduced = 1, True
        self.occupation_shifts = occupation_shifts
        self.graph_program = bool(graph_program)
        self.table, self.node_index = table, node_index
        self._node, self._items, self._combined = None, (), None
        self._hash = hash(('compiled-polynomial', signature))

    @property
    def is_zero(self): return self.zero
    def __hash__(self): return self._hash
    def __eq__(self, other):
        return type(other) is type(self) and self.signature == other.signature
    def __fingerprint__(self): return ('compiled-polynomial', self.program_version, self.signature)
    def __repr__(self): return 'CompiledPolynomial[{}; {} instructions]'.format(self.signature[:12], len(self.instructions))
    def format_expr(self): return repr(self)

    @classmethod
    def compile(cls, polynomial, cache=None):
        """Compile shared callable nodes without unrolling coordinate contexts."""
        polynomial = PolyPath.from_polynomial(polynomial)
        cache = {} if cache is None else cache
        programs = cache.setdefault('programs', {})
        interned = cache.setdefault('interned', {})
        axes = cache.setdefault('axes', {})
        constants = cache.setdefault('constants', {})
        constant_axes = cache.setdefault('constant_axes', {})

        def describe(node):
            if not any(node.order):
                value = node.evaluate_polynomial(np.zeros((1, node.ndim), dtype=int),
                                                 node_cache=constants, axis_cache=constant_axes)
                return 'constant', (), np.asarray(value).reshape(-1)[0]
            if node._node is None:
                terms = []
                for term, weight in node._items:
                    factors = []
                    for index, axis in enumerate(term.axes):
                        identity = id(axis)
                        if identity not in axes:
                            coefficients = np.array(axis.coeffs, copy=True)
                            coefficients.flags.writeable = False
                            axes[identity] = axis, coefficients, dev.fingerprint(coefficients)
                        _, coefficients, signature = axes[identity]
                        factors.append((index, 0, coefficients, signature))
                    terms.append((weight, np.sqrt(2.)**term.steps, tuple(factors)))
                return 'leaf', (), tuple(terms)
            kind, arguments = node._node
            if kind == 'add':
                pending, children = list(arguments), []
                while pending:
                    child = pending.pop()
                    if child._node is not None and child._node[0] == 'add':
                        pending.extend(child._node[1])
                    else:
                        children.append(child)
                return 'sum', tuple(children), None
            if kind == 'linear':
                return 'linear', tuple(child for child, _ in arguments), tuple(weight for _, weight in arguments)
            if kind == 'scale': return kind, (arguments[0],), arguments[1]
            if kind in {'shift', 'remap'}:
                return kind, (arguments[0],), np.asarray(arguments[1], dtype=np.int64)
            if kind == 'permutation_sum':
                return kind, (arguments[0],)*len(arguments[2]), np.asarray(arguments[2], dtype=np.int64)
            if kind == 'project':
                return kind, (arguments[0],), (arguments[0].ndim, np.asarray(arguments[1], dtype=np.int64))
            if kind == 'mul':
                return kind, arguments[:2], tuple(np.asarray(pull, dtype=np.int64) for pull in arguments[4:6])
            raise ValueError('unsupported polynomial compilation operation: ' + kind)

        pending = [(polynomial, None)]
        while pending:
            node, description = pending.pop()
            if id(node) in programs: continue
            if description is None:
                description = describe(node)
                if description[1]:
                    pending.append((node, description))
                    pending.extend((child, None) for child in reversed(description[1]))
                    continue
            kind, children, extra = description
            children = tuple(programs[id(child)][1] for child in children)
            signature = dev.fingerprint((cls.program_version, kind,
                tuple(child.signature for child in children), extra, node.ndim, node.order))
            result = interned.get(signature)
            if result is None:
                result = cls(((kind, children, extra, ()),), node.ndim, node.order,
                             signature, node.is_zero, graph_program=True)
                interned[signature] = result
            programs[id(node)] = node, result
        return programs[id(polynomial)][1]

    @staticmethod
    def _leaf_value(terms, states, axis_cache):
        count = len(states)
        result = np.zeros(count)
        for weight, divisor, axes in terms:
            factors = []
            for index, offset, coefficients, coefficient_key in axes:
                values = np.zeros(count, dtype=states.dtype) if index < 0 else states[:, index]
                if offset: values = values+offset
                key = ('compiled-axis', coefficient_key, values.dtype.str, values.shape, values.tobytes())
                factor = axis_cache.get(key)
                if factor is None:
                    factor = np.polynomial.polynomial.polyval(values, coefficients)
                    axis_cache[key] = factor
                factors.append(factor)
            product = np.prod(factors, axis=0) if factors else np.ones(count)
            result += weight*product/divisor
        return result

    def _evaluate_graph_program(self, states, node_cache, axis_cache):
        def identity(program, data):
            return program, (data.dtype.str, data.shape, data.tobytes())
        root = identity(self, states)
        prior = node_cache.get(root)
        if prior is not None: return prior
        # A numerical call stack executes already compiled child calls and
        # coordinate bindings. It constructs no symbolic/topological plan and
        # releases completed call frames instead of retaining an unrolled DAG.
        frames = [[self, states, root, 0, []]]
        while frames:
            program, data, key, cursor, values = frames[-1]
            kind, children, extra, _ = program.instructions[0]
            if cursor < len(children):
                child = children[cursor]
                if kind == 'shift': argument = data[:, :child.ndim]+extra[None, :]
                elif kind == 'remap': argument = data[:, extra]
                elif kind == 'permutation_sum': argument = data[:, extra[cursor]]
                elif kind == 'project':
                    ndim, keep = extra
                    argument = np.zeros((len(data), ndim), dtype=data.dtype)
                    argument[:, keep] = data[:, :len(keep)]
                elif kind == 'mul': argument = data[:, extra[cursor]]
                else: argument = data
                frames[-1][3] += 1
                child_key = identity(child, argument)
                cached = node_cache.get(child_key)
                if cached is None:
                    frames.append([child, argument, child_key, 0, []])
                else:
                    values.append(cached)
                continue
            if kind == 'leaf': result = self._leaf_value(extra, data, axis_cache)
            elif kind == 'constant': result = np.full(len(data), extra)
            elif kind == 'scale': result = extra*values[0]
            elif kind == 'mul': result = values[0]*values[1]
            elif kind in {'shift', 'remap', 'project'}: result = values[0]
            elif kind in {'sum', 'linear', 'permutation_sum'}:
                result = np.zeros(len(data))
                for index, value in enumerate(values):
                    result += extra[index]*value if kind == 'linear' else value
            else: raise ValueError('unsupported compiled polynomial call: ' + kind)
            node_cache[key] = result
            frames.pop()
            if frames: frames[-1][4].append(result)
        return result

    def evaluate_polynomial(self, substates, node_cache=None, axis_cache=None):
        node_cache = {} if node_cache is None else node_cache
        axis_cache = {} if axis_cache is None else axis_cache
        if self.table is not None:
            return self.table.evaluate(self.node_index, substates, node_cache, axis_cache)
        if self.graph_program:
            return self._evaluate_graph_program(substates, node_cache, axis_cache)
        state_key = (substates.dtype.str, substates.shape, substates.tobytes())
        root_key = self, state_key
        cached = node_cache.get(root_key)
        if cached is not None:
            return cached
        values = []
        count = len(substates)
        for instruction in self.instructions:
            kind, arguments, extra, releases = instruction
            if kind == 'constant': result = np.full(count, arguments)
            elif kind == 'leaf':
                result = np.zeros(count)
                for weight, divisor, axes in arguments:
                    factors = []
                    for index, offset, coefficients, coefficient_key in axes:
                        axis_key = ('compiled-axis', coefficient_key, index, offset, state_key)
                        factor = axis_cache.get(axis_key)
                        if factor is None:
                            states = np.zeros(count, dtype=substates.dtype) if index < 0 else substates[:, index]
                            if offset: states = states + offset
                            factor = np.polynomial.polynomial.polyval(states, coefficients)
                            axis_cache[axis_key] = factor
                        factors.append(factor)
                    product = np.prod(factors, axis=0) if factors else np.ones(count)
                    result += weight * product / divisor
            elif kind == 'alias': result = values[arguments[0]]
            elif kind == 'scale': result = extra * values[arguments[0]]
            elif kind == 'mul': result = values[arguments[0]] * values[arguments[1]]
            elif kind in {'sum', 'linear'}:
                result = np.zeros(count)
                weights = (1,)*len(arguments) if kind == 'sum' else extra
                for index, weight in zip(arguments, weights):
                    result += values[index] if kind == 'sum' else weight * values[index]
            else: raise ValueError('unsupported compiled polynomial instruction: ' + kind)
            values.append(result)
            for index in releases: values[index] = None
        result = values[-1]
        node_cache[root_key] = result
        return result


class PTEnergyChangeProductSum(TensorCoefficientPoly, PolynomialInterface):
    """
    A representation of a sum of 1/energy * poly sums
    which is here so we can transpose energy change indices intelligently
    """
    def __init__(self, terms: dict, prefactor=1, canonicalize=True, reduced=False):
        if any(getattr(row[0],'_conditional_endpoint',False) for key in terms for row in key):
            terms = canonicalize_condition_paths(terms)
        if prefactor != 1:
            terms = {
                k: v.scale(prefactor)
                for k, v in terms.items()
            }
        # An empty factor product is the unit denominator. Its dimension
        # comes from the lazy numerator, not from a fictitious energy row.
        if canonicalize or () in terms:
            self._ndim = self.get_key_ndim(terms)
        else:
            self._ndim = 0 if not terms else len(next(iter(terms.keys()))[0])
        super().__init__(terms, prefactor=1, canonicalize=canonicalize)
        self.reduced = reduced
        # if canonicalize:
        #     for k in terms:
        #         if any(all(s == 0 for s in sk) for sk in k): raise ValueError("huh")
        # self.audit()

    def __fingerprint__(self):
        return (type(self).__module__,type(self).__qualname__,
                self.terms,self.prefactor,self.ndim if self.terms or self._ndim is not None else 0,self.reduced)

    def mutate(self, terms:dict=default, prefactor:'Any'=default,  canonicalize:'Any'=default, reduced:'Any'=default):
        return type(self)(
            self.terms if terms is default else terms,
            canonicalize=(terms is not default) if canonicalize is default else canonicalize,
            prefactor=self.prefactor if prefactor is default else prefactor,
            reduced=self.reduced if reduced is default else reduced
        )

    @classmethod
    def from_polynomial(cls, polynomial):
        """Embed a native polynomial as polynomial/1 without materializing it."""
        if isinstance(polynomial, cls):
            return polynomial
        if nput.is_numeric(polynomial):
            if polynomial == 0:
                return cls({})
            polynomial = PolyPath.from_coeffs([], prefactor=polynomial, steps=0)
        if not isinstance(polynomial, (ProductPTPolynomial, ProductPTPolynomialSum)):
            raise TypeError("can't promote {} to {}".format(type(polynomial), cls.__name__))
        return cls({(): polynomial})

    def __add__(self, other):
        if nput.is_numeric(other) and other == 0:
            return self
        if isinstance(other, (ProductPTPolynomial, ProductPTPolynomialSum)) or nput.is_numeric(other):
            other = type(self).from_polynomial(other)
        if not isinstance(other, PTEnergyChangeProductSum):
            return NotImplemented
        if () not in self.terms and () not in other.terms:
            return super().__add__(other)
        # Retain the lazy numerator objects and each genuine denominator key.
        new_terms = dict(self.terms)
        for key, polynomial in other.terms.items():
            value = polynomial if key not in new_terms else new_terms[key] + polynomial
            if nput.is_zero(value):
                new_terms.pop(key, None)
            else:
                new_terms[key] = value
        return type(self)(new_terms) if new_terms else 0

    def __radd__(self, other):
        return self.__add__(other)

    def flip_energy_terms(self):
        return self.mutate(
            {
                tuple(
                    (side if getattr(side, '_state_projection_factor', False) else (side + 1) % 2,) + k
                    for side, k in self.side_change_iter(key)
                ):term
                for key,term in self.terms.items()
            }
        )

    def prep_serialization_dict(self):
        new_data = {}
        for k,p in self.terms.items():
            if not isinstance(p, ProductPTPolynomialSum):
                p = ProductPTPolynomialSum([p])
            new_data[k] = p.prep_serialization_dict()
        return new_data
    @classmethod
    def from_serialization_dict(cls, big_dict):
        new_terms = {}
        for k, subdict in big_dict.items():
            new_terms[k] = ProductPTPolynomialSum.from_serialization_dict(subdict)
        return cls(new_terms)

    def to_state(self, serializer=None):# don't let this be directly serialized
        raise NotImplementedError(...)
        return {
            'term_data': self.get_serialization_data(),
            'prefactor': self.prefactor
        }
    @classmethod
    def from_state(cls, state, serializer=None):
        terms = {
            tuple(tuple(k) for k in key):p
            for key,p in serializer.deserialize(state['terms'])
        }
        prefactor = state['prefactor']
        return cls(terms, prefactor=prefactor, canonicalize=False)

    def filter(self, terms, mode='match'):
        # Energy filters leave denominator-free polynomials unchanged, as before
        # their promotion into this container.
        if any(isinstance(row[0], _ChangeDomainIndex) for key in self.terms for row in key):
            filtered = {}
            for key, poly in self.terms.items():
                divisors = tuple(row for row in key if not isinstance(row[0], _ChangeDomainIndex))
                # Integer summation support is independent of frequency
                # filtering; denominator-free equations remain unchanged.
                if not divisors or type(self)({divisors: 1}).filter(terms, mode=mode).terms:
                    filtered[key] = poly
        elif any(getattr(row[0],'_conditional_endpoint',False) for key in self.terms for row in key):
            filtered = filter_condition_paths(self.terms,terms,mode)
        else:
            filtered = dict(super().filter(terms, mode=mode).terms)
        if () in self.terms:
            filtered[()] = self.terms[()]
        return self.mutate(filtered)
    @classmethod
    def canonical_key(cls, monomial_tuple):
        l = max(len(l) for l in monomial_tuple) if len(monomial_tuple) > 0 else 0
        return tuple(
                t + (0,) * (l - len(t))
                for t in monomial_tuple
            )
        # return super().canonical_key(
        #     tuple(
        #         t + (0,) * (l - len(t))
        #         for t in monomial_tuple
        #     )
        # )

    @classmethod
    def side_change_iter(cls, key):
        return ((k[0], k[1:]) for k in key)
    @classmethod
    def format_key(self, key):
        if not key: return '1'
        if any(getattr(row[0],'_conditional_endpoint',False) for row in key):
            return format_condition_path(key)
        return "E-[{}]".format(
            "x".join(
                "[R]" if side == 1 else "[L]" + (
                    ("({})" if len(k) > 1 else "{}").format(
                        "+".join(
                            "{s}w{i}".format(s=("" if s == 1 else "-" if s == -1 else s), i=i)
                            for i,s in enumerate(k)
                        )
                    )
                )
                for side,k in self.side_change_iter(key)
                )
        )

    def __repr__(self):
        sums = []
        for k_prod,v in self.terms.items():
            sums.append("{}{}".format(v, self.format_key(k_prod)))
        return "ESum({})".format("+".join(sums) if len(sums) < 3 else "\n      +".join(sums))

    @classmethod
    def format_energy_prod_key(cls, key):
        if not key: return '1'
        if any(getattr(row[0],'_conditional_endpoint',False) for row in key):
            return format_condition_path(key)
        substr = [
            "[{}](".format("R" if side == 1 else "L")
            + "+".join(
                "{s}w[{i}]".format(s=("" if s == 1 else "-" if s == -1 else s), i=i)
                for i, s in enumerate(k)
                if s != 0
            ) + ")"
            for side,k in cls.side_change_iter(key)
        ]
        bb = "*".join(substr)
        return ("({})" if "*" in bb else "{}").format(bb)

    def format_expr(self):
        """
        Formats in a Mathematica-ingestible format
        :return:
        """
        if any(getattr(row[0],'_conditional_endpoint',False) for key in self.terms for row in key):
            return format_condition_sum(self)
        keys = []
        substrings = []
        for k,poly in self.terms.items():
            keys.append(self.format_energy_prod_key(k))
            substrings.append(poly.format_expr())

        join = " + " if sum(len(x) for x in substrings) < 50 else "\n  + "
        prefac = (
            "" if self.prefactor == 1 else
            "-" if self.prefactor == 1 else
            str(ProductPTPolynomial.prep_float(self.prefactor))
        )
        key_prods = [
            prefac + (
                "({})/{}"
                    if "+" in s else
                "{}/{}"
            ).format(s, k)
            for k,s in zip(keys, substrings)
        ]
        return join.join(key_prods)

    @classmethod
    def shift_key(cls, key, shift):
        return tuple(
            (side,) + tuple(k + (shift[i] if i < len(shift) else 0) for i,k in enumerate(ksum))
            for side,ksum in cls.side_change_iter(key)
        )
    def shift_energies(self, shift):
        new_terms = {}
        for k, p in self.terms.items():
            k2 = self.shift_key(k, shift)
            # if all(any(s != 0 for s in sk) for sk in k2):
            new_terms[k2] = p
        # if len(new_terms) == 0:
        #     raise ValueError(self.terms, shift)
        return self.mutate(new_terms)
    def shift(self, shift, shift_energies=False):
        new_terms = {}
        for k,p in self.terms.items():
            if shift_energies:
                k2 = self.shift_key(k, shift)
                # if all(any(s != 0 for s in sk) for sk in k2):
                new_terms[k2] = p.shift(shift)
            else:
                # Fixed-reference resolvents and actual-state projectors are
                # diagonal functions of the intermediate occupation. They
                # move with it, unlike relative energy-change denominators.
                shifted = tuple(
                    ((side,) + tuple(value + (shift[i] if i < len(shift) else 0)
                                     for i, value in enumerate(delta)))
                    if getattr(side, '_state_projection_factor', False) else (side,) + delta
                    for side, delta in self.side_change_iter(k)
                )
                new_terms[shifted] = p.shift(shift)
        return self.mutate(new_terms)

    def scale(self, scaling):
        if nput.is_numeric(scaling) and scaling == 1: return self
        return self.mutate({
            k: v.scale(scaling)
            for k, v in self.terms.items()
        })

    @classmethod
    def get_key_ndim(cls, terms:dict):
        return max((
            max((len(change) for _,change in cls.side_change_iter(key)), default=0)
            if key else getattr(polynomial, 'ndim', 0)
            for key,polynomial in terms.items()
        ), default=0)
    @property
    def ndim(self):
        if self._ndim is None:
            self._ndim = self.get_key_ndim(self.terms)
        return self._ndim

    def audit(self, target=None, ignore_constants=True):
        # if target is None:
        #     target = self.ndim
        for change,poly in self.terms.items():
            if any(c[0] not in {0, 1} and not getattr(c[0], '_state_projection_factor', False) for c in change):
                raise ValueError("change {} invalid, L/R index must be 0/1".format(change))
            width = len(change[0]) - 1 if change else self.ndim
            if target is not None and width != target:
                raise ValueError(
                    "({}) energy change {} doesn't have target dimension ({})".format(
                        type(self)({change: poly}),
                        change, target
                    )
                )
            if not nput.is_numeric(poly):
                try:
                    poly.audit(width, ignore_constants=ignore_constants)
                except ValueError:
                    raise ValueError("bad polynomial for changes {}: {}".format(
                        change,
                        poly
                    ))
    def ensure_dimension(self, ndim):
        nd = self.ndim
        if nd >= ndim: return self

        new_terms = {}
        for key,polys in self.terms.items():
            new_key = tuple(k + (0,) * (ndim - (len(k)-1)) for k in key)
            new_terms[new_key] = polys.ensure_dimension(ndim) if not nput.is_numeric(polys) else polys

        return self.mutate(new_terms)

    def sort(self):
        return self.mutate(
            {
                k: self.terms[k]
                for k in sorted(self.terms.keys())
            }
        )

    def _permute_changes(self, changes, new_inds):
        # rem = ProductPTPolynomial.fast_ind_remainder(len(changes), new_inds)
        new = tuple(
            changes[i] for i in new_inds
        ) + changes[len(new_inds):]
        return new
    def permute(self, new_inds, check_perm=True, allow_padding=False):
        new_terms = {}
        for energy_changes,polys in self.terms.items():
            new_ech = tuple(
                (side,)+self._permute_changes(ec, new_inds)
                for side,ec in self.side_change_iter(energy_changes)
            )
            new_terms[new_ech] = polys.permute(new_inds, check_perm=check_perm, allow_padding=allow_padding)
        return self.mutate(new_terms)

    def permutation_sum(self, permutations, check_perm=True, allow_padding=False):
        """Group a symmetry orbit by energy key before building polynomial nodes."""
        permutations = tuple(tuple(int(i) for i in perm) for perm in permutations)
        new_terms = {}
        for energy_changes, polys in self.terms.items():
            groups = {}
            for permutation in permutations:
                new_ech = tuple(
                    (side,) + self._permute_changes(ec, permutation)
                    for side, ec in self.side_change_iter(energy_changes)
                )
                groups.setdefault(new_ech, []).append(permutation)
            for new_ech, subpermutations in groups.items():
                if hasattr(polys, 'permutation_sum'):
                    contribution = polys.permutation_sum(
                        subpermutations,
                        check_perm=check_perm,
                        allow_padding=allow_padding
                    )
                else:
                    contribution = sum(
                        polys.permute(
                            permutation,
                            check_perm=check_perm,
                            allow_padding=allow_padding
                        )
                        for permutation in subpermutations
                    )
                new_terms[new_ech] = new_terms.get(new_ech, 0) + contribution
        return self.mutate(new_terms)

    @staticmethod
    def _check_neg(t1, t2):
        if len(t1) != len(t2):
            return False
        else:
            return all(
                all(ttt1 == -ttt2 for ttt1, ttt2 in zip(tt1, tt2))
                for tt1, tt2 in zip(t1, t2)
            )

    @staticmethod
    def find_term_scaling(key):
        new_keys = []
        new_scalings = []
        for k in key:
            s = np.gcd.reduce(k)
            new_scalings.append(s)
            new_keys.append([kk/s for kk in k])

        return tuple(new_keys), np.prod(new_scalings)

    # def factor_energies(self):
    #     ...


    def combine_energies(self):
        raise NotImplementedError(...)
        # terms = self.factor_energies()
        terms = self.terms

        base_keys = list(terms.keys())
        elim_pos = set()
        merges = set()
        unmerged = set()
        for i,k1 in enumerate(base_keys):
            if i not in elim_pos:
                for j,k2 in enumerate(base_keys[i+1:]):
                    j = i+1 + j
                    if j not in elim_pos:
                        if self._check_neg(k1, k2):
                            elim_pos.add(i)
                            elim_pos.add(j)
                            merges.add((k1, k2))
                            break
                else:
                    unmerged.add(k1)

        new_terms = {}
        for k1,k2 in merges:
            # sanity check
            if k1 in unmerged or k2 in unmerged:
                raise ValueError('fuck')
            new_terms[k1] = terms[k1] + terms[k2].scale(-1)
        for k in unmerged:
            new_terms[k] = terms[k]

        return new_terms

    def combine(self, combine_subterms=True, combine_energies=False):
        if self.reduced: return self

        new_terms = {}
        # raise Exception(
        #     list(self.terms.keys()),
        #     list(self.combine_energies().keys())
        # )
        base_terms = self.combine_energies() if combine_energies else self.terms
        for k,p in base_terms.items():
            if combine_subterms and isinstance(p, ProductPTPolynomialSum):
                p = p.combine()
                if isinstance(p, PolyPath):
                    if not p.is_zero:
                        new_terms[k] = p
                elif len(p.polys) > 0:
                    new_terms[k] = p
            else:
                new_terms[k] = p
        return self.mutate(new_terms, reduced=True)

    @staticmethod
    def _permute_idx(idx, inds, rem, ndim):
        """
        Builds a new change index by aligning everything that matches in inds and
        incrementing everything that doesn't in remainder
        """
        # if len(inds) == 0: return idx
        # if not nput.is_numeric(inds[0]):
        #     inds = inds[1] # we remap the right inds
        if not isinstance(inds, tuple): inds = tuple(inds)
        # if remainder is None:
        #     rem = ProductPTPolynomial.fast_ind_remainder(max_remapping, inds)
        # else:
        #     rem = remainder[1]
        side, idx = idx[0], idx[1:]
        new = tuple(
            idx[inds[i]]
                if i < len(inds) else
            idx[rem[i - len(inds)] if i - len(inds) < len(rem) else i]
            for i in range(ndim)
        )
        if all(n == 0 for n in new) and not getattr(side, '_state_projection_factor', False):
            raise ValueError(idx, inds, new)
        return (side,) + new

    @classmethod
    def _build_new_echange_key(cls, k1, k2, inds, remainder, right_state_shift=None,
                              left_ndim=None, right_ndim=None):
        l1 = (len(k1[0]) - 1 if k1 else
              len(inds[0])+len(remainder[0]) if left_ndim is None else left_ndim)
        l2 = (len(k2[0]) - 1 if k2 else
              len(inds[1])+len(remainder[1]) if right_ndim is None else right_ndim)
        i1 = len(inds[0])
        r1 = len(remainder[0])
        i2 = len(inds[1])
        r2 = len(remainder[1])
        new_key = tuple(
            (side,)
            + tuple(k[i] for i in inds[0])
            + tuple(k[j] for j in remainder[0])
            + k[i1 + r1:]
            + (0,) * (l2 - i2)
            for side,k in cls.side_change_iter(k1)
        ) + tuple(
            (side,)
            + tuple(k[i] for i in inds[1])
            + (0,) * (l1 - i1)
            + tuple(k[j] for j in remainder[1])
            + k[i2 + r2:]
            for side,k in cls.side_change_iter(k2)
        )
        if right_state_shift is not None:
            # Changes on left-only axes also shift every diagonal state
            # factor on the right, even though the right tensor has no
            # coefficient index on those axes.
            extra = ((0,) * i1
                     + tuple(right_state_shift[i] if i < len(right_state_shift) else 0
                             for i in remainder[0]))
            new_key = new_key[:len(k1)] + tuple(
                (row[0],) + tuple(value + (extra[i] if i < len(extra) else 0)
                                 for i, value in enumerate(row[1:]))
                if getattr(row[0], '_state_projection_factor', False) else row
                for row in new_key[len(k1):]
            )
        return new_key
    @classmethod
    def _pad_echange_key_right(cls, k, ndim, inds, remainder, key_ndim=None):
        ndim = max([ndim, len(inds[1])])#, len(remainder[1])])
        return cls._build_new_echange_key(k, ((0,)*(ndim+1),), inds, remainder, left_ndim=key_ndim)[:-1]
    @classmethod
    def _pad_echange_key_left(cls, k, ndim, inds, remainder):
        ndim = max([ndim, len(inds[0])])#, len(remainder[0])])
        return cls._build_new_echange_key(((0,)*(ndim+1),), k, inds, remainder)[1:]
    def mul_along(self, other:'PolynomialInterface', inds, remainder=None, mapping=None, right_state_shift=None):
        """
        We multiply every subpoly along the given indices, transposing the appropriate tensor indices

        :param other:
        :param inds:
        :param remainder:
        :return:
        """
        if isinstance(other, PTEnergyChangeProductSum):
            new_poly = self.direct_product(other,
                                       key_func=lambda k1,k2:self._build_new_echange_key(k1, k2, inds, remainder, right_state_shift,
                                           left_ndim=getattr(self.terms[k1], 'ndim', 0) if not k1 else None,
                                           right_ndim=getattr(other.terms[k2], 'ndim', 0) if not k2 else None),
                                       mul=lambda a,b,inds=inds:a.mul_along(b, inds, remainder=remainder, mapping=mapping)
                                       )
            if () in self.terms or () in other.terms:
                new_poly = new_poly.mutate(new_poly.terms)
        elif isinstance(other, (ProductPTPolynomial, ProductPTPolynomialSum)):
            new_poly = self.mutate(
                {
                    self._pad_echange_key_right(k, other.ndim, inds, remainder, key_ndim=getattr(v, 'ndim', 0)): (
                        other.scale(v)
                            if nput.is_numeric(v) else
                        v.mul_along(other, inds, remainder=remainder, mapping=mapping)
                    )
                    for k, v in self.terms.items()
                }
            )
        else:
            raise NotImplementedError(type(other))
        return new_poly

    def rmul_along(self, other, inds, remainder=None, mapping=None, right_state_shift=None):
        """
        We multiply every subpoly along the given indices, transposing the appropriate tensor indices

        :param other:
        :param inds:
        :param remainder:
        :return:
        """
        if isinstance(other, PTEnergyChangeProductSum):
            return other.mul_along(self, inds, remainder=remainder, mapping=mapping)
        elif isinstance(other, (ProductPTPolynomial, ProductPTPolynomialSum)):
            return self.mutate(
                {
                    self._build_new_echange_key(((0,)*(other.ndim+1),), k, inds, remainder, right_state_shift, right_ndim=getattr(v, 'ndim', 0))[1:]: (
                        other.scale(v)
                            if nput.is_numeric(v) else
                        other.mul_along(v, inds, remainder=remainder, mapping=mapping)
                    )
                    for k, v in self.terms.items()
                }
            )
        else:
            raise NotImplementedError(type(other))

    @staticmethod
    def _simple_product_energy_key(k1, k2, left, right):
        key = k1+k2
        if k1 and k2:
            return key
        width = max(getattr(left, 'ndim', 0), getattr(right, 'ndim', 0),
                    max((len(row)-1 for row in key), default=0))
        return tuple(row+(0,)*(width+1-len(row)) for row in key)

    def mul_simple(self, other:'ProductPTPolynomial'):
        """
        We multiply every subpoly along the given indices, transposing the appropriate tensor indices

        :param other:
        :param inds:
        :param remainder:
        :return:
        """
        if isinstance(other, (PTEnergyChangeProductSum)):
            result = self.direct_product(other,
                                       key_func=lambda k1,k2:self._simple_product_energy_key(
                                           k1,k2,self.terms[k1],other.terms[k2]),
                                       mul=lambda a,b:a.mul_simple(b)
                                       )
            return result.mutate(result.terms) if (() in self.terms or () in other.terms) else result
        elif isinstance(other, (ProductPTPolynomial, ProductPTPolynomialSum)):
            return self.mutate(
                {
                    k: (
                        other.scale(v)
                            if isinstance(v, (int, float, np.integer, np.floating)) else
                        v.mul_simple(other)
                    )
                    for k, v in self.terms.items()
                }
            )
        elif nput.is_numeric(other):
            if other == 1:
                return self
            elif other == 0:
                return 0
            else:
                raise ValueError(other)
        else:
            return other.rmul_simple(self)
    def rmul_simple(self, other: 'PolynomialInterface') -> 'PolynomialInterface':
        if isinstance(other, (ProductPTPolynomial, ProductPTPolynomialSum)):
            # pad = (0,)*other.ndim
            return type(self)(
                {
                    # tuple(
                    #     pad + k
                    #     for k in k_tup
                    # )
                    k_tup: other.scale(v) if nput.is_numeric(v) else other.mul_simple(v)
                    for k_tup, v in self.terms.items()
                },
                prefactor=self.prefactor
            )
        elif nput.is_numeric(other):
            if other == 1:
                return self
            elif other == 0:
                return 0
            else:
                raise ValueError(other)
        else:
            return other.mul_simple(self)

class _HamiltonianCoefficient(int):
    """A saved linear recipe for one normal-coordinate ladder channel.

    Numerical inputs retain their V/G conventions. Requests are prepared once
    during expression construction and serialized verbatim, including the
    source-axis permutations and rational weights; binding never derives them.
    """
    _cache = {}

    def __new__(cls, rank, positive, requests=None):
        rank, positive = int(rank), int(positive)
        building = requests is None
        if rank < 3 or not (rank+1)//2 <= positive <= rank:
            raise ValueError('invalid effective Hamiltonian channel')
        if requests is None:
            key = rank, positive
            if key in cls._cache: return cls._cache[key]
            pairs = math.comb(rank, 2)
            signs = (1,)*positive+(-1,)*(rank-positive)
            requests = [(0, tuple(range(rank)), 1, 1)]
            for i,j in itertools.combinations(range(rank), 2):
                permutation = (i,)+tuple(k for k in range(rank) if k not in (i,j))+(j,)
                requests.append((1, permutation, -signs[i]*signs[j], pairs))
            requests = tuple(requests)
        else:
            requests = tuple((int(kind), tuple(map(int, perm)), int(num), int(den))
                             for kind,perm,num,den in requests)
            if not requests or any(kind not in (0,1) or sorted(perm) != list(range(rank)) or den <= 0
                   for kind,perm,num,den in requests):
                raise ValueError('invalid effective Hamiltonian recipe')
        result = int.__new__(cls, 8+positive)
        result.rank, result.positive, result.requests = rank, positive, requests
        result._frozen = True
        if building: cls._cache[(rank,positive)] = result
        return result

    def __setattr__(self,name,value):
        if getattr(self,'_frozen',False): raise AttributeError('coefficient recipes are immutable')
        object.__setattr__(self,name,value)

    def __reduce__(self): return type(self), (self.rank,self.positive,self.requests)
    def __fingerprint__(self): return ('hamiltonian-coefficient',1,self.rank,self.positive,self.requests)

    def canonical_indices(self, indices):
        left = tuple(sorted(indices[:self.positive]))
        right = tuple(sorted(indices[self.positive:]))
        if len(left) == len(right) and right < left: left,right = right,left
        return left+right

    @staticmethod
    def operator_dependencies(operator):
        order, kind = operator
        return {(order,source) for source,_,_,_ in kind.requests} if isinstance(kind,_HamiltonianCoefficient) else {operator}

    @classmethod
    def inactive_operators(cls, operators, active):
        active = set(active)
        return {operator for operator in operators if not cls.operator_dependencies(operator)&active}

    @classmethod
    @functools.lru_cache(maxsize=32)
    def basis(cls, rank):
        """Select an exact basis from O(N**2) candidate sign patterns."""
        from fractions import Fraction
        pairs = tuple(itertools.combinations(range(rank),2))
        signs, seen = [], set()
        for count in range(3):
            for negative in itertools.combinations(range(rank),count):
                row = tuple(-1 if i in negative else 1 for i in range(rank))
                if row[0] < 0: row = tuple(-s for s in row)
                if row not in seen: seen.add(row); signs.append(row)
        matrix = [(Fraction(1),)+tuple(Fraction(-row[i]*row[j],len(pairs)) for i,j in pairs)
                  for row in signs]
        pivots = nput.rational_matrix_pivots(tuple(zip(*matrix)))
        selected = [matrix[i] for i in pivots]
        if len(selected) != 1+len(pairs): raise ValueError('incomplete ladder-channel basis')
        return tuple(signs[i] for i in pivots), nput.rational_matrix_inverse(selected)


class _BoundHamiltonianCoefficient:
    """Gather a bound saved recipe without materializing a dense tensor."""
    def __init__(self, recipe, sources):
        self.recipe, self.sources = recipe, tuple(sources)
        shapes = [source.shape for source in sources if np.ndim(source)]
        if shapes and any(shape != shapes[0] for shape in shapes):
            raise ValueError('incompatible potential/kinetic coefficient shapes')
        self.shape = () if not shapes else shapes[0]
        self.ndim = len(self.shape)
        if self.ndim and self.ndim != recipe.rank:
            raise ValueError('effective coefficient tensor rank differs from its recipe')

    def __getitem__(self, indices):
        result = 0
        for source, permutation, numerator, denominator in self.recipe.requests:
            tensor = self.sources[source]
            if np.ndim(tensor):
                result = result + (numerator/denominator)*tensor[tuple(indices[i] for i in permutation)]
        return result


class PTTensorCoeffProductSum(TensorCoefficientPoly, PolynomialInterface):
    """
    A representation for a sum of tensor coefficients * polynomial sums
    which is primarily here so we can transpose tensor coefficients intelligently
    """

    tensor_ids = ["V", "G", "U", "Z", "u", "M", "LV", "LG"] # we explicitly split Coriolis and Watson terms out

    def __init__(self, terms, prefactor=1, canonicalize=True, ndim=None, inds_map=None, reduced=False):
        # Contradictory conditional paths can reduce a coefficient to zero
        # before a later permutation. Keep zero as an absent term, rather
        # than an integer where the polynomial protocol is expected.
        terms = {key:value for key,value in terms.items() if not nput.is_zero(value)
                 and not (isinstance(value,PTEnergyChangeProductSum) and not value.terms)}
        if prefactor != 1:
            terms = {
                k: v.scale(prefactor)
                for k, v in terms.items()
            }
        super().__init__(terms, prefactor=1, canonicalize=canonicalize)
        self._inds_map = {} if inds_map is None else inds_map
        self.reduced = reduced
        self._ndim = ndim
        # if canonicalize:
        #     self.terms = self.ensure_subpoly_dim(self.terms)
        # self.audit()

    def __fingerprint__(self):
        return (type(self).__module__,type(self).__qualname__,
                self.terms,self.prefactor,self.ndim if self.terms or self._ndim is not None else 0,self.reduced)

    def mutate(self,
               terms:'Any'=default, *, prefactor:'Any'=default,
               ndim:'Any'=default, inds_map:'Any'=default,
               canonicalize:'Any'=default, reduced:'Any'=default
               ):
        return type(self)(
            self.terms if terms is default else terms,
            canonicalize=(terms is not default) if canonicalize is default else canonicalize,
            prefactor=self.prefactor if prefactor is default else prefactor,
            reduced=self.reduced if reduced is default else reduced,
            ndim=self._ndim if ndim is default else ndim,
            inds_map=self._inds_map if inds_map is default else inds_map
        )

    def prep_serialization_dict(self):
        new_data = {}
        for k, p in self.terms.items():
            if not isinstance(p, PTEnergyChangeProductSum):
                p = PTEnergyChangeProductSum({((-100,),):p})
            new_data[k] = p.prep_serialization_dict()
        return new_data
    @classmethod
    def from_serialization_dict(cls, big_dict):
        new_terms = {}
        for k,subdict in big_dict.items():
            test_key = next(iter(subdict.keys())) if len(subdict) == 1 else None
            if len(subdict) == 1 and test_key == ((-100,),):
                vals = next(iter(subdict.values()))
                poly = ProductPTPolynomialSum.from_serialization_dict(vals)
            else:
                poly = PTEnergyChangeProductSum.from_serialization_dict(subdict)
            new_terms[k] = poly
        return cls(new_terms)

    def to_state(self, serializer=None): # don't let this be directly serialized
        raise NotImplementedError(...)
        return {
            'term_data': self.get_serialization_data(),
            'prefactor': self.prefactor
        }
    @classmethod
    def from_state(cls, state, serializer=None): # don't let this be directly serialized
        raise NotImplementedError(...)
        terms = {
            tuple(tuple(k) for k in key):p
            for key,p in serializer.deserialize(state['terms'])
        }
        prefactor = state['prefactor']
        return cls(terms, prefactor=prefactor, canonicalize=False)

    def flip_energy_terms(self):
        return self.mutate(
            {
                k: t.flip_energy_terms() if hasattr(t, 'flip_energy_terms') else t
                for k, t in self.terms.items()
            }
        )

    def filter(self, terms, mode='match'):
        subterms = super().filter(terms, mode=mode)
        return self.mutate(subterms.terms) # base 'filter' doesn't
    def filter_coefficients(self, terms, mode='match'):
        return self.filter(terms, mode=mode)
    def filter_energies(self, terms, mode='match'):
        return self.mutate({
            k:v.filter(terms, mode=mode) if isinstance(v, PTEnergyChangeProductSum) else v
            for k,v in self.terms.items()
        })

    @classmethod
    def format_key(cls, key):
        return "".join([
            "{}[{}]{}".format(
                ('E{}'.format(k[1].positive) if isinstance(k[1],_HamiltonianCoefficient) else cls.tensor_ids[k[1]]),
                k[0],
                k[2:]
            )
            for k in key
        ])
    def __repr__(self):
        sums = []
        for k_prod,v in self.terms.items():
            sums.append("{}{}".format(v, self.format_key(k_prod)))
        return "TSum({})".format("+".join(sums) if len(sums) < 3 else "\n      +".join(sums))

    @classmethod
    def format_tensor_key(cls, key):
        return "".join([
            "{}[{}][{}]".format(
                ('E{}'.format(k[1].positive) if isinstance(k[1],_HamiltonianCoefficient) else cls.tensor_ids[k[1]]),
                k[0],
                ", ".join(str(i) for i in k[2:])
            )
            for k in key
        ])
    def format_expr(self):
        """
        Formats in a Mathematica-ingestible format
        :return:
        """
        keys = []
        substrings = []
        for k,esum in self.terms.items():
            keys.append(self.format_tensor_key(k))
            substrings.append(esum.format_expr())

        join = " + " if sum(len(x) for x in substrings) < 50 else "\n  + "
        key_prods = [
            (
                "{}*({})"
                    if "+" in s else
                "{}*{}"
                    if s != "" else
                "{}"
            ).format(k,s.replace("\n", "\n" + (" "*(len(k) + len(join)))))
            for k,s in zip(keys, substrings)
        ]
        return join.join(key_prods)

    @classmethod
    def coeff_product_inds(cls, key, return_counts=False):
        if len(key) == 0:
            empty = np.array([], dtype=int)
            return (empty, empty) if return_counts else empty
        return np.unique(np.concatenate([c[2:] for c in key]), return_counts=return_counts)
    def get_inds(self, key):
        if key not in self._inds_map:
            inds = self.coeff_product_inds(key)
            self._inds_map[key] = inds
        return self._inds_map[key]

    def prune_operators(self, ops):
        ops = set(ops)
        return self.mutate(
            {
                t: p
                for t, p in self.terms.items()
                if all(tt[:2] not in ops for tt in t)
            }
        )

    @property
    def ndim(self):
        if self._ndim is None:
            self._ndim = max(len(self.get_inds(coeff_inds)) for coeff_inds in self.terms.keys())
        return self._ndim
    def audit(self, target=None, required_dimension=None, ignore_constants=True):
        """
        Checks to ensure that the number of dimensions aligns with
        the number of indices in the tensor coefficients

        :return:
        """
        for coeff_indices,poly_stuff in self.terms.items():
            num_inds = len(self.get_inds(coeff_indices))
            if required_dimension is not None:
                _, counts = self.coeff_product_inds(coeff_indices, return_counts=True)
                for c, d in zip(counts, required_dimension):
                    if c < d or (c%2 != d%2): raise ValueError("bad indices {} for change order {}".format(
                        coeff_indices, required_dimension
                    ))
            if target is not None and num_inds < target: raise ValueError("too few indices ({}) for changes ({}) in {}".format(
                num_inds, target, poly_stuff
            ))
            try:
                poly_stuff.audit(target=num_inds, ignore_constants=ignore_constants)
            except ValueError:
                raise ValueError("mismatch between inds {} and poly {}".format(
                    coeff_indices,
                    poly_stuff.format_expr()
                )) from None
            # self._audit_dimension(poly_stuff, num_inds)
    def ensure_dimension(self, ndim):
        nd = self.ndim
        if nd >= ndim: return self

        return self.mutate(
            {k:p.ensure_dimension(ndim) for k,p in self.terms.items()},
            ndim=ndim,
        )
    @classmethod
    def ensure_subpoly_dim(cls, terms):
        new_terms = {}
        for k, p in terms.items():
            new_terms[k] = p.ensure_dimension(len(cls.coeff_product_inds(k)))
        return new_terms

    def sort(self):
        return self.mutate(
            {
                k: self.terms[k] for k in sorted(self.terms.keys())
            }
        )

    def permute(self, new_inds, check_perm=True, allow_padding=False):
        new_terms = {}
        inv_map = np.argsort(new_inds)
        for coeff_inds,polys in self.terms.items():
            new_key = tuple(
                ci[:2] + tuple(
                    inv_map[j] if j < len(new_inds) else j
                    for j in ci[2:]
                )
                for ci in coeff_inds
            )
            new_terms[new_key] = polys.permute(new_inds, check_perm=check_perm, allow_padding=allow_padding)
        return self.mutate(new_terms, canonicalize=True)

    def _check_equiv(self, k1, k2):
        if len(k1) != len(k2):
            return None
        if any(len(kk1) != len(kk2) or kk1[:2] != kk2[:2] for kk1, kk2 in zip(k1, k2)):
            return None
        k1_inds = self.get_inds(k1)
        k2_inds = self.get_inds(k2)
        if len(k1_inds) != len(k2_inds) or any(ki1 != ki2 for ki1, ki2 in zip(k1_inds, k2_inds)):
            return None

        for p in itertools.permutations(range(len(k1_inds))):
            swap_inds = [k2_inds[i] for i in p]
            k2_sub = tuple(ci[:2] + tuple(swap_inds[j] for j in ci[2:]) for ci in k2)
            if k1 == k2_sub:
                return p # permutation to apply to our polys

    def combine_terms(self):
        base_keys = list(self.terms.keys())
        elim_pos = set()
        merges = {}
        unmerged = set()
        for i,k1 in enumerate(base_keys):
            if i not in elim_pos:
                for j,k2 in enumerate(base_keys[i+1:]):
                    j = i+1 + j
                    if j not in elim_pos:
                        remapping = self._check_equiv(k1, k2)
                        if remapping is not None:
                            elim_pos.add(i)
                            elim_pos.add(j)
                            merges[(k1, k2)] = remapping
                            break
                else:
                    unmerged.add(k1)

        new_terms = {}
        for k1,k2 in merges:
            # sanity check
            if k1 in unmerged or k2 in unmerged:
                raise ValueError('fuck')
            new_terms[k1] = self.terms[k1] + self.terms[k2].permute(merges[(k1, k2)])
        for k in unmerged:
            new_terms[k] = self.terms[k]

        return new_terms

    @classmethod
    def reindex_terms(cls, terms):
        new_terms = {}
        for k,v in terms.items():
            ind_vars = [iii for kk in k for iii in kk[2:]]
            if len(ind_vars) == 0:
                new_terms[k] = v
                continue
            perm = np.unique(ind_vars)
            if len(perm) == perm[-1] + 1:
                new_terms[k] = v
                continue
            submap = np.argsort(perm) # just to remap things down
            mapping = dict(zip(perm, submap))
            new_key = tuple(kk[:2] + tuple(mapping[iii] for iii in kk[2:]) for kk in k)
            new_terms[new_key] = v.permute(perm)
        return new_terms
    def combine(self, combine_coeffs=False, combine_subterms=True, combine_energies=False):
        if self.reduced: return self
        new_terms = {}

        if combine_coeffs:
            raise NotImplementedError("coeff combining too subtle for current evaluator strategy")

        base_terms = self.combine_terms() if combine_coeffs else self.terms
        # base_terms = self.reindex_terms(base_terms)
        for k,p in base_terms.items():
            if combine_subterms and isinstance(p, ProductPTPolynomialSum):
                p = p.combine()
                if isinstance(p, PolyPath):
                    if not p.is_zero:
                        new_terms[k] = p
                elif len(p.polys) > 0:
                    new_terms[k] = p
            elif isinstance(p, (PTTensorCoeffProductSum, PTEnergyChangeProductSum)):
                p = p.combine(combine_subterms=combine_subterms, combine_energies=combine_energies)
                if len(p.terms) > 0:
                    new_terms[k] = p
            else:
                new_terms[k] = p
        return self.mutate(new_terms, reduced=True)

    @classmethod
    def _get_uperms(cls, p_list):
        return PerturbationTheoryTermProduct._get_uperms(p_list)
    def free_up_indices(self, start, stop):
        new_terms = {}
        for k,poly in self.terms.items():
            existing_inds = self.get_inds(k)
            perm_idx = [1] * (stop - start) + [0] * (len(existing_inds) - stop)
            perms, _ = self._get_uperms(perm_idx)
            full_perms = np.concatenate([
                np.broadcast_to(np.arange(start)[np.newaxis], (len(perms), start)),
                perms + start
            ], axis=1)
            for perm in full_perms:
                inv_map = np.argsort(perm)
                # for coeff_inds, polys in self.terms.items():
                new_key = tuple(
                    ci[:2] + tuple(
                        inv_map[j] #if j < len(new_inds) else j
                        for j in ci[2:]
                    )
                    for ci in k
                )
                new_poly = poly.permute(perm)
                new_terms[new_key] = new_terms.get(new_key, 0) + new_poly

            # _ = poly.permute(full_perms[0])
            # for p in full_perms[1:]:
            #     _ += poly.permute(p)
            # poly = _
            # new_terms[k] = poly
        return self.mutate(new_terms, canonicalize=True)
    def shift_energies(self, change):
        new_terms = {}
        for k,p in self.terms.items():
            new_terms[k] = p.shift_energies(change) if isinstance(p, PTEnergyChangeProductSum) else p
        return self.mutate(new_terms)
    def shift(self, shift):
        new_terms = {}
        for k,p in self.terms.items():
            new_terms[k] = p.shift(shift)
        return self.mutate(new_terms)
    def scale(self, scaling):
        if nput.is_numeric(scaling) and scaling == 1: return self
        return self.mutate({
            k:v.scale(scaling)
            for k,v in self.terms.items()
        })
    @staticmethod
    def _potential_symmetrizer(idx):
        return tuple(sorted(idx))
    @staticmethod
    def _kinetic_symmetrizer(idx): # first and last and everything in between
        p_sorts = list(sorted([idx[0], idx[-1]]))
        rest_sorts = list(sorted(idx[1:-1]))
        return (p_sorts[0],) + tuple(rest_sorts) + (p_sorts[1],)
    @staticmethod
    def _coriolis_symmetrizer(idx): # second and last are momenta and everything else is x
        # return idx
        # symmetry of the Coriolis term is obviously weird since it's not even Hermitian???
        if idx[1] > idx[-1] or (idx[1] == idx[-1] and idx[0] > idx[2]): # only case where we can flip stuff around...
            rest = idx[-2:1:-1] + (idx[0],)
            return (rest[0], idx[-1]) + rest[1:] + (idx[1],)
        else:
            return idx
        # p_sorts = list(sorted([idx[1], idx[-1]]))
        # rest_sorts = list(sorted((idx[0],) + idx[2:-1]))
        # new = (rest_sorts[0], p_sorts[0]) + tuple(rest_sorts[1:]) + (p_sorts[1],)
        # return new

    default_symmetrizers = None
    @classmethod
    def symmetrizers(cls):
        if cls.default_symmetrizers is None:
            return (
                cls._potential_symmetrizer,
                cls._kinetic_symmetrizer,
                cls._potential_symmetrizer,  # V'
                cls._coriolis_symmetrizer,
                cls._potential_symmetrizer,  # U
                cls._potential_symmetrizer,  # M
                cls._potential_symmetrizer,  # H (ext.) # needed for extra G-matrix terms, need to make this less hard coded
            )
        else:
            return cls.default_symmetrizers
    @classmethod
    def _symmetrize(cls, idx, symmetrizers=None):
        if isinstance(idx[1],_HamiltonianCoefficient):
            return idx[:2]+idx[1].canonical_indices(idx[2:])
        if symmetrizers is None:
            symmetrizers = cls.symmetrizers()
        return (idx[0], idx[1]) + symmetrizers[idx[1]](idx[2:])
    @classmethod
    def canonical_key(cls, monomial_tuple, symmetrizers=None):
        return super().canonical_key(tuple(cls._symmetrize(t, symmetrizers=symmetrizers) for t in monomial_tuple))

    @classmethod
    def _handle_mul(cls,
                    new_key, total_order,
                    num_fixed, mul_size,
                    poly_1, poly_2,
                    mul_inds, mul_remainder,
                    og_remainder,
                    left_baseline, baseline,
                    left_choice_x, right_perm_x,
                    left_rem_x, right_rem_x,
                    num_left, num_right,
                    logger, log_level
                    ):
        if baseline is not None:
            right_shift = [0] * num_right
            for s, i in zip(left_baseline, mul_inds[1]):
                right_shift[i] = s
            shift_poly_2 = poly_2
            if any(s != 0 for s in right_shift):
                shift_poly_2 = shift_poly_2.shift(right_shift)
                logger.log_print("baseline: {baseline}", baseline=baseline,
                                 log_level=log_level
                                 )
                logger.log_print("shift: {shift}", shift=right_shift,
                                 log_level=log_level
                                 )
                logger.log_print("{p}", p=shift_poly_2,
                                 preformatter=lambda **vars: dict(vars, p=vars['p'].format_expr()),
                                 log_level=log_level
                                 )
            if isinstance(poly_1, PTEnergyChangeProductSum):
                new_poly = poly_1.mul_along(shift_poly_2, mul_inds, mul_remainder,
                                             right_state_shift=baseline)
            elif isinstance(shift_poly_2, PTEnergyChangeProductSum):
                new_poly = shift_poly_2.rmul_along(poly_1, mul_inds, mul_remainder,
                                                 right_state_shift=baseline)
            else:
                new_poly = poly_1.mul_along(shift_poly_2, mul_inds, mul_remainder)
        else:
            new_poly = poly_1.mul_along(poly_2, mul_inds, mul_remainder)

        # finally, we need to make sure that the indices are in the right order
        # by noting that the remainder inds that were specified must be ordered
        # after the proper inds
        # so we build a permutation to apply to new_poly by constructing a proxy index
        # ordered like the indices of new_poly but with values from the target position indices
        pad_left = num_fixed
        pad_right = pad_left + len(og_remainder[0])
        pad_full = pad_right + len(og_remainder[1])
        pad_full_left = pad_right + len(og_remainder[1])
        pad_full_right = num_left + len(og_remainder[0])
        true_index = (
                [0] * num_fixed  # the first indices don't need to be reordered
                + [
                    pad_left + i if i < len(og_remainder[0]) else  # left remainder before right remainder
                    pad_right + j if j < len(og_remainder[1]) else
                    pad_full + i
                    for i, j in zip(left_choice_x, right_perm_x)
                ]
                + [
                    pad_left + i if i < len(og_remainder[0]) else
                    pad_full_left + i
                    for i in left_rem_x
                ]
                + [
                    pad_right + j if j < len(og_remainder[1]) else
                    pad_full_right + j
                    for j in right_rem_x
                ]
        )
        perm = np.argsort(true_index)  # gives the position each final index should map to
        inv_perm = np.argsort(perm)
        new_poly = new_poly.permute(perm, allow_padding=False)
        old_key = new_key
        new_key = tuple(
            ci[:2] + tuple(inv_perm[i] for i in ci[2:])
            for ci in new_key
        )
        if (
                np.any(np.not_equal(cls.coeff_product_inds(new_key), cls.coeff_product_inds(old_key)))
                or len(cls.coeff_product_inds(new_key)) != num_left + num_right - num_fixed - mul_size
        ):
            raise ValueError(old_key, new_key)
        if _DEBUG_PRINT:
            try:
                new_poly.audit(num_left + num_right - num_fixed - mul_size)
            except:
                logger.log_print([
                    "Audit failed for: {p} = {f}",
                    "Inds: {m}",
                    "Rem: {r}",
                    "Perm: {perm}"
                ],
                    p=new_poly,
                    m=mul_inds,
                    r=mul_remainder,
                    perm=perm,
                    log_level=log_level,
                    preformatter=lambda **vars: dict(vars, f=vars['p'].format_expr()),
                )
                raise
        else:
            new_poly.audit(num_left + num_right - num_fixed - mul_size)

        _, counts = cls.coeff_product_inds(new_key, return_counts=True)
        for c, d in zip(counts, total_order):
            if c < d or (c % 2 != d % 2): raise ValueError(
                "bad indices {} for expected order {}".format(
                    new_key, total_order
                ))
        for c in counts[len(total_order):]:
            if (c % 2 != 0): raise ValueError(
                "bad indices {} for expected order {}".format(
                    new_key, total_order
                ))

        return new_key, new_poly
    def _generate_direct_product_values(self,
                                        inds, remainder, index_classes, baseline,
                                        k1, k2, poly_1, poly_2
                                        ):

        left_inds = self.get_inds(k1)
        right_inds = self.get_inds(k2)
        num_left = len(left_inds)
        num_right = len(right_inds)
        num_fixed = len(inds[1])

        # we multiply indices on the left with the corresponding indices on the right
        left_fixed_inds, right_fixed_inds = [tuple(x) for x in inds]
        left_fixed_remainder, right_fixed_remainder = [tuple(y) for y in remainder]

        # there may be indices that are not specified in the overlap or remainder blocks
        # which need to be distributed over all specified indices
        left_float_remainder = tuple(range(len(left_fixed_inds) + len(left_fixed_remainder), num_left))
        right_float_remainder = tuple(range(len(right_fixed_inds) + len(right_fixed_remainder), num_right))

        left_remainder_inds = left_fixed_remainder + left_float_remainder
        right_remainder_inds = right_fixed_remainder + right_float_remainder

        base_dim = num_left + num_right - num_fixed

        num_free_left = len(left_float_remainder)
        num_free_right = len(right_float_remainder)

        num_defd_left = len(left_fixed_remainder)
        num_defd_right = len(right_fixed_remainder)
        # we take every possible combo of remaining indices from left and right
        max_mult = min(
            num_right - num_fixed,
            num_left - num_fixed
        ) # the max number of extra axes to align
        # if min_dim is not None:
        #     max_mult = min(base_dim - min_dim, max_mult)

        # def _build_remapping(fixed_inds, remainder_inds, float_inds, num_pad_rem, num_pad_float):
        #     mapping = {k:i for i,k in enumerate(fixed_inds)}
        #     rem_pad = len(fixed_inds) + rem_pad
        #     for i,k in enumerate(remainder_inds):
        #         if mapping[k] = i + len(fixed_inds) + rem_pad
        #     for i,k in enumerate(float_inds):
        #
        #     mapping = {k:padding+i for i,k in enumerate(fixed_inds)}


        # for debug purposes we construct the total order we'd expect to have
        _, left_order = self.coeff_product_inds(k1, return_counts=True)
        _, right_order = self.coeff_product_inds(k2, return_counts=True)
        total_order = [
            left_order[i] + right_order[j]
            for i,j in zip(left_fixed_inds, right_fixed_inds)
        ] + [
            left_order[i]
            for i in left_fixed_remainder
        ] + [
            right_order[j]
            for j in right_fixed_remainder
        ]

        log_level = Logger.LogLevel.Normal if _DEBUG_PRINT else Logger.LogLevel.MoreDebug
        logger = PerturbationTheoryTerm.default_logger()
        with logger.block(tag="Broadcast Mult: {inds} {remainder} {k1} {k2} | {max_mult}",
                          inds=inds, remainder=remainder, k1=k1, k2=k2, max_mult=max_mult,
                          log_level=log_level):
            logger.log_print("1: {p1}", p1=poly_1,
                             preformatter=lambda **vars: dict(vars, p1=vars['p1'].format_expr()),
                             log_level=log_level
                             )
            logger.log_print("2: {p2}", p2=poly_2,
                             preformatter=lambda **vars: dict(vars, p2=vars['p2'].format_expr()),
                             log_level=log_level
                             )
            for mul_size in range(max_mult+1):
                # we _assume_ that all permutations are handled directly by the polynomial
                # construction, and so we really just choose to overlap different numbers of terms here,
                # assuming a strict ordering, and then permuting all the remainder indices to preserve symmetry

                # computed once and cached
                # if mul_size == 0:
                #     combo_inds_left = [()]
                #     combo_rem_left = [np.arange(len(left_remainder_inds))]
                #     combo_inds_right = [()]
                #     combo_rem_right = [np.arange(len(right_remainder_inds))]
                # else:
                # We assume that all of the floating indices are identical and the symmetry is clean
                # so instead of taking all possible combinations, we consider blocks where
                # we take differing numbers of elements from the proper remainders and the floating indices
                # To do this, for now we filter out "bad" combinations from the full set of possibilities
                # by noting that for any given choice of

                uidx_left, uidx_right = index_classes
                rem_uleft = np.array([uidx_left[u] for u in left_fixed_remainder])
                rem_uright = np.array([uidx_right[u] for u in right_fixed_remainder])
                for nx_left in range(max(mul_size-num_free_left, 0), min(num_defd_left, num_free_right, mul_size)+1):
                    # m = how many terms to take from the defined remainder
                    subsize = mul_size - nx_left

                    # sample from the left remainder inds and pad with remainder inds
                    left_specd_choice = PerturbationTheoryTermProduct.get_combination_inds(num_defd_left, nx_left)
                    left_specd_rem = PerturbationTheoryTermProduct.get_combination_comp(num_defd_left, nx_left)
                    left_float_choice = np.arange(num_defd_left, num_defd_left+subsize)
                    left_float_rem = np.arange(num_defd_left+subsize, len(left_remainder_inds))

                    # we now need to reduce over the unique index classes
                    if len(uidx_left) > 0:
                        uchoice_left = nput.vector_take(rem_uleft, left_specd_choice)
                        _, uinds_left = np.unique(uchoice_left, axis=0, return_index=True)
                        left_specd_choice = left_specd_choice[uinds_left,]
                        left_specd_rem = left_specd_rem[uinds_left,]

                    combo_inds_left = np.concatenate([
                        left_specd_choice,
                        np.broadcast_to(left_float_choice[np.newaxis], (len(left_specd_choice), subsize))
                    ], axis=1)
                    combo_rem_left = np.concatenate([
                        left_specd_rem,
                        np.broadcast_to(left_float_rem[np.newaxis],
                                        (len(left_specd_rem), len(left_float_rem)))
                    ], axis=1)

                    # for the right inds, we need to pad initially by `m` remainder indices
                    # then we need to sample from the number of _right_ remainders to map onto the left
                    # choice, where again we need to choose how many of the right remainders to sample at once
                    right_match_choice = np.arange(num_defd_right, num_defd_right+nx_left)
                    # right_span = [max(mul_size-num_free_right, 0), min(num_defd_right, subsize)+1]
                    # logger.log_print("{s} {b} {d} {e}", b=subsize, d=num_defd_right, e=num_free_right, s=right_span)
                    for nx_right in range(max(mul_size-num_free_right, 0), min(num_defd_right, subsize)+1):

                        sub_subsize = subsize - nx_right # how many free indices we have in the multiplication
                        num_muld_right = num_defd_right+nx_left

                        right_specd_choice = PerturbationTheoryTermProduct.get_combination_inds(num_defd_right, nx_right)
                        right_specd_rem = PerturbationTheoryTermProduct.get_combination_comp(num_defd_right, nx_right)
                        right_float_choice = np.arange(num_muld_right, num_muld_right+sub_subsize)
                        right_float_rem = np.arange(num_muld_right+sub_subsize, len(right_remainder_inds))

                        # we now need to reduce over the unique index classes
                        if len(uidx_right) > 0:
                            uchoice_right = nput.vector_take(rem_uright, right_specd_choice)
                            _, uinds_right = np.unique(uchoice_right, axis=0, return_index=True)
                            right_specd_choice = right_specd_choice[uinds_right,]
                            right_specd_rem = right_specd_rem[uinds_right,]

                        # finally we build out these indices
                        combo_inds_right = np.concatenate([
                            np.broadcast_to(right_match_choice[np.newaxis], (len(right_specd_choice), nx_left)),
                            right_specd_choice,
                            np.broadcast_to(right_float_choice[np.newaxis], (len(right_specd_choice), subsize-nx_right))
                        ], axis=1)
                        combo_rem_right = np.concatenate([
                            right_specd_rem,
                            np.broadcast_to(right_float_rem[np.newaxis], (len(right_specd_rem), len(right_float_rem)))
                        ], axis=1)

                        # combo_inds_left = PerturbationTheoryTermProduct.get_combination_inds(len(left_remainder_inds), mul_size)
                        # combo_rem_left = PerturbationTheoryTermProduct.get_combination_comp(len(left_remainder_inds), mul_size)
                        # combo_inds_right = PerturbationTheoryTermProduct.get_combination_inds(len(right_remainder_inds), mul_size)
                        # combo_rem_right = PerturbationTheoryTermProduct.get_combination_comp(len(right_remainder_inds), mul_size)

                        # print("!", m, "l", num_defd_left, num_free_left, 'r', num_defd_right, num_free_right)
                        # print(combo_inds_left, left_specd_choice, left_float_choice)
                        # print(left_fixed_remainder, left_float_remainder)
                        # print(combo_inds_right, right_match_choice, right_specd_choice, right_float_choice)
                        # print(right_fixed_remainder, right_float_remainder)
                        for left_choice_x, left_rem_x in zip(combo_inds_left, combo_rem_left):
                            left_choice = tuple(left_remainder_inds[k] for k in left_choice_x)
                            left_mul_inds = left_fixed_inds + left_choice

                            # now we build our new key, noting that the multiplied inds will be the first _n_
                            # inds in the new set, then the left inds, then the right ones
                            left_mapping = {k: i for i, k in enumerate(left_mul_inds)}
                            left_rem_rem = [left_remainder_inds[x] for x in left_rem_x]
                            for j, k in enumerate(left_rem_rem):
                                left_mapping[k] = j + num_fixed + mul_size
                            left_key = tuple(
                                        t[:2] + tuple(left_mapping[k] for k in t[2:])
                                        for t in k1
                                    )

                            if baseline is not None:
                                # sometimes unchanged inds on the right get multiplied by changed left inds
                                left_baseline = [
                                    baseline[x] if x < len(baseline) else 0
                                    for x in left_mul_inds
                                ]
                            else:
                                left_baseline = None

                            # for right_perm in [right_choice]:
                            for right_choice_x, right_rem_x in zip(combo_inds_right, combo_rem_right):
                                # for right_perm_x in [right_choice_x]:#itertools.permutations(right_choice_x):
                                    # if any(
                                    #         i < len(remainder[0]) and j < len(remainder[1])
                                    #         for i, j in zip(left_choice_x, right_perm_x)
                                    # ): continue  # we can't multiply two remainder indices

                                right_choice = tuple(right_remainder_inds[k] for k in right_choice_x)

                                right_mul_inds = right_fixed_inds + right_choice

                                right_mapping = {k: i for i, k in enumerate(right_mul_inds)}
                                right_rem_rem = [right_remainder_inds[x] for x in right_rem_x]
                                for j,k in enumerate(right_rem_rem):
                                    right_mapping[k] = j + num_left
                                right_key = tuple(
                                    t[:2] + tuple(right_mapping[k] for k in t[2:])
                                    for t in k2
                                )

                                new_key = left_key + right_key

                                mul_inds = [
                                    left_mul_inds, # self
                                    right_mul_inds  # other
                                ]
                                mul_remainder = [
                                    left_rem_rem,
                                    right_rem_rem
                                ]

                                # if len(left_mul_inds) != len(right_mul_inds):
                                #     raise ValueError(
                                #         left_mul_inds, right_mul_inds
                                #     )
                                #
                                # if len(left_rem_rem) != len(right_rem_rem):
                                #     raise ValueError(
                                #         left_rem_rem, right_rem_rem,
                                #         nx_left, nx_right
                                #     )

                                with logger.block(tag="{m} {mul_inds} {mul_remainder}",
                                                  m=mul_size, mul_inds=mul_inds, mul_remainder=mul_remainder,
                                                  log_level=log_level):

                                    new_key, new_poly = self._handle_mul(
                                        new_key, total_order,
                                        num_fixed, mul_size,
                                        poly_1, poly_2,
                                        mul_inds, mul_remainder,
                                        remainder,
                                        left_baseline, baseline,
                                        left_choice_x, right_choice_x,
                                        left_rem_x, right_rem_x,
                                        num_left, num_right,
                                        logger, log_level
                                    )


                                    perm_blocks = []
                                    # we never permute these
                                    if num_fixed > 0:
                                        perm_blocks.append([np.arange(num_fixed)])

                                    perm_pad = num_fixed
                                    if nx_left > 0:
                                        diff_pos = np.nonzero(np.diff(rem_uleft))
                                        for sub_block_inds in np.split(left_fixed_remainder, diff_pos[0] + 1):
                                            sub_proxy = [
                                                1
                                                    if x in left_choice else
                                                0
                                                for x in sub_block_inds
                                            ]
                                            sub_perms, _ = self._get_uperms(sub_proxy)
                                            perm_blocks.append(perm_pad + sub_perms)
                                            perm_pad = perm_pad + len(sub_block_inds)
                                    else:
                                        if num_defd_left > 0:
                                            perm_blocks.append([np.arange(perm_pad, perm_pad+num_defd_left)])
                                            perm_pad = perm_pad + num_defd_left

                                    if nx_right > 0:
                                        diff_pos = np.nonzero(np.diff(rem_uright))
                                        for sub_block_inds in np.split(right_fixed_remainder, diff_pos[0] + 1):
                                            sub_proxy = [
                                                1
                                                    if x in right_choice else
                                                0
                                                for x in sub_block_inds
                                            ]
                                            sub_perms, _ = self._get_uperms(sub_proxy)
                                            perm_blocks.append(perm_pad + sub_perms)
                                            perm_pad = perm_pad + len(sub_block_inds)
                                    else:
                                        if num_defd_right > 0:
                                            perm_blocks.append([np.arange(perm_pad, perm_pad+num_defd_right)])
                                            perm_pad = perm_pad + num_defd_right

                                    # now any indices that were generated by multiplying or
                                    # ignoring the "free" remainder indices need to be symmetrized
                                    perm_idx = [  # TODO: get this faster from nx_left and nx_right
                                                   2 for i, j in zip(left_choice_x, right_choice_x)
                                                   if i >= num_defd_left and j >= num_defd_right
                                               ] + [
                                                   1 for i in left_rem_x
                                                   if i >= num_defd_left
                                               ] + [
                                                   0 for j in right_rem_x
                                                   if j >= num_defd_right
                                               ]
                                    if len(perm_idx) > 0:
                                        free_perms, _ = self._get_uperms(perm_idx)
                                        num_prev = num_fixed + num_defd_left + num_defd_right
                                        perm_blocks.append(num_prev + free_perms)
                                    group_permutations = polynomial_uses_path(new_poly)
                                    permutation_groups = {} if group_permutations else None
                                    for perm_bits in itertools.product(*perm_blocks):
                                        perm = np.concatenate(perm_bits) if perm_bits else np.array([], dtype=int)
                                        inv_map = np.argsort(perm)
                                        perm_key = tuple(
                                            ci[:2] + tuple(
                                                inv_map[j]
                                                for j in ci[2:]
                                            )
                                            for ci in new_key
                                        )
                                        perm_key = self.canonical_key(perm_key)

                                        logger.log_print("{k} [{r}]",
                                                         k=perm_key,
                                                         r=inv_map,
                                                         preformatter=lambda **vars: dict(vars, k=self.format_tensor_key(vars['k'])),
                                                         log_level=log_level
                                                         )
                                        if group_permutations:
                                            permutation_groups.setdefault(perm_key, []).append(tuple(perm))
                                        else:
                                            perm_poly = new_poly.permute(perm)
                                            logger.log_print("{p}", p=perm_poly,
                                                             preformatter=lambda **vars: dict(vars, p=vars['p'].format_expr()),
                                                             log_level=log_level
                                                             )
                                            yield perm_key, perm_poly

                                    for perm_key, permutations in (
                                            permutation_groups.items() if group_permutations else ()
                                    ):
                                        if hasattr(new_poly, 'permutation_sum'):
                                            perm_poly = new_poly.permutation_sum(permutations)
                                        else:
                                            perm_poly = sum(new_poly.permute(perm) for perm in permutations)
                                        logger.log_print("{p}", p=perm_poly,
                                                         preformatter=lambda **vars: dict(vars, p=vars['p'].format_expr()),
                                                         log_level=log_level
                                                         )

                                        yield perm_key, perm_poly#.ensure_dimension(total_dimension)

    def _adjust_key_right(self, k, other, inds, remainder):
        # TODO: make sure other isn't fucked up somehow
        mapping = np.argsort(np.concatenate([inds[0], remainder[0]]))
        return tuple(
            ci[:2] + tuple(
                mapping[i] if i < len(mapping) else i + len(remainder[1])
                for i in ci[2:]
            )
            for ci in k
        )
    def _adjust_key_left(self, k, other, inds, remainder):
        raise NotImplementedError(...)
        # TODO: make sure other isn't fucked up somehow
        mapping = np.argsort(np.concatenate([inds[1], remainder[1]]))
        return tuple(
            ci[:2] + tuple(mapping[i] if i < len(mapping) else i for i in ci[2:])
            for ci in k
        )
    def mul_along(self, other:'PolynomialInterface', inds, remainder=None, index_classes=None,
                  mapping=None, baseline=None):
        """
        We multiply every subpoly along the given indices, transposing the appropriate tensor indices

        :param other:
        :param inds:
        :param remainder:
        :return:
        """
        if isinstance(other, PTTensorCoeffProductSum):
            # we need to take the direct product of aligning (or not)
            # for all of the untouched inds
            new = self.direct_multiproduct(
                other,
                lambda *kargs: self._generate_direct_product_values(inds, remainder,
                                                                    index_classes, baseline,
                                                                    *kargs
                                                                    )
            )
        elif isinstance(other, (PTEnergyChangeProductSum, ProductPTPolynomial, ProductPTPolynomialSum)):
            new = self.mutate(
                {
                    self._adjust_key_right(k, other, inds, remainder):other.rmul_along(p, inds, remainder=remainder, mapping=mapping)
                    for k,p in self.terms.items()
                },
                canonicalize=True
            )
        else:
            raise NotImplementedError(type(other))
        return new
    def rmul_along(self, other:'ProductPTPolynomial', inds, remainder=None, mapping=None):
        """
        We multiply every subpoly along the given indices, transposing the appropriate tensor indices

        :param other:
        :param inds:
        :param remainder:
        :return:
        """
        if isinstance(other, PTTensorCoeffProductSum):
            new = other.mul_along(self, inds, remainder=remainder, mapping=mapping)
            if self._inds_map is not other._inds_map:
                self._inds_map.update(other._inds_map) # these could all be shared for even more speed up
            new._inds_map = self._inds_map
            return new
        elif isinstance(other, PTEnergyChangeProductSum):
            return self.mutate(
                {
                    self._adjust_key_left(k, other, inds, remainder):other.mul_along(p, inds, remainder=remainder, mapping=mapping)
                    for k,p in self.terms.items()
                },
                canonicalize=True
            )
        else:
            raise NotImplementedError(type(other))

    def mul_simple(self, other:'ProductPTPolynomial'):
        """
        We multiply every subpoly along the given indices, transposing the appropriate tensor indices

        :param other:
        :param inds:
        :param remainder:
        :return:
        """
        if isinstance(other, PTTensorCoeffProductSum):
            return self.direct_product(other,
                                       # key_func=lambda k1,k2:k1+tuple(
                                       #     k[:2] + tuple(i+self.ndim for i in k[2:])
                                       #     for k in k2
                                       # ),
                                       mul=lambda a,b:a.mul_simple(b)
                                       )
        elif isinstance(other, PTEnergyChangeProductSum):
            return self.mutate(
                {
                    k: v.mul_simple(other.ensure_dimension(len(self.get_inds(k))))
                    for k, v in self.terms.items()
                },
                canonicalize=False
            )
        elif isinstance(other, (ProductPTPolynomial, ProductPTPolynomialSum)):
            return self.mutate(
                {
                    k: v.mul_simple(other)
                    for k, v in self.terms.items()
                },
                canonicalize=False
            )
        else:
            return other.rmul_simple(self)
            # raise NotImplementedError(type(other))
    def rmul_simple(self, other: 'PolynomialInterface') -> 'PolynomialInterface':
        if isinstance(other, PTEnergyChangeProductSum):
            return self.mutate(
                {
                    # tuple(
                    #     k[:2] + tuple(i + other.ndim for i in k[2:])
                    #     for k in k_tup
                    # )
                    k_tup: other.ensure_dimension(len(self.get_inds(k_tup)).mul_simple(v))
                    for k_tup, v in self.terms.items()
                },
                canonicalize=False
            )
        elif isinstance(other, (ProductPTPolynomial, ProductPTPolynomialSum)):
            return self.mutate(
                {
                    # tuple(
                    #     k[:2] + tuple(i + other.ndim for i in k[2:])
                    #     for k in k_tup
                    # )
                    k_tup: other.scale(v) if nput.is_numeric(v) else other.mul_simple(v)
                    for k_tup, v in self.terms.items()
                },
                canonicalize=False
            )
        elif nput.is_numeric(other):
            if other == 1:
                return self
            elif other == 0:
                return 0
            else:
                raise ValueError(other)
        else:
            return other.mul_simple(self)


class PTTensorCoeffProductDAG(PTTensorCoeffProductSum):
    """Lazy, canonical operation DAG for tensor-coefficient expressions.

    Leaves retain the existing dictionary-backed ``PTTensorCoeffProductSum``.
    Algebra builds immutable nodes and materializes the legacy representation
    only at compatibility boundaries such as serialization or the current
    tensor evaluator.  Materialization is cached on every node, so repeated
    evaluation never replays the derivation tree.
    """

    _node_cache = weakref.WeakValueDictionary()
    _node_requests = 0
    _node_hits = 0
    _materializations = 0

    def __init__(self, node, ndim=None, reduced=False):
        self._node = node
        self._ndim = ndim
        self._inds_map = {}
        self.reduced = reduced
        self.prefactor = 1
        self._eager = None
        self._hash = None
        self._operator_keys = None

    def __fingerprint__(self):
        # Retain the symbolic DAG, never its cached operator/domain analyses.
        return (type(self).__module__,type(self).__qualname__,self._node,self._ndim,self.reduced)

    @staticmethod
    def _freeze(arg):
        if isinstance(arg, np.ndarray):
            return tuple(PTTensorCoeffProductDAG._freeze(x) for x in arg.tolist())
        if isinstance(arg, dict):
            return tuple(sorted(
                (
                    PTTensorCoeffProductDAG._freeze(key),
                    PTTensorCoeffProductDAG._freeze(value)
                )
                for key, value in arg.items()
            ))
        if isinstance(arg, (set, frozenset)):
            return tuple(sorted(PTTensorCoeffProductDAG._freeze(x) for x in arg))
        if isinstance(arg, list):
            return tuple(PTTensorCoeffProductDAG._freeze(x) for x in arg)
        if isinstance(arg, tuple):
            return tuple(PTTensorCoeffProductDAG._freeze(x) for x in arg)
        return arg

    @classmethod
    def from_sum(cls, expression):
        if isinstance(expression, cls):
            return expression
        if not isinstance(expression, PTTensorCoeffProductSum):
            raise TypeError("can't make a tensor DAG from {}".format(type(expression)))
        return cls(('leaf', expression), ndim=expression.ndim, reduced=expression.reduced)

    @classmethod
    def _from_node(cls, kind, args, ndim=None, reduced=False):
        args = cls._freeze(args)
        if Diagnostics.enabled:
            cls._node_requests += 1
        key = (kind, args, ndim, reduced)
        node = cls._node_cache.get(key)
        if node is None:
            node = cls((kind,) + args, ndim=ndim, reduced=reduced)
            cls._node_cache[key] = node
        else:
            if Diagnostics.enabled:
                cls._node_hits += 1
        return node

    @classmethod
    def clear_caches(cls):
        cls._node_cache.clear()
        cls._node_requests = 0
        cls._node_hits = 0
        cls._materializations = 0

    @classmethod
    def cache_info(cls):
        return {
            'tensor_nodes': len(cls._node_cache),
            'tensor_node_requests': cls._node_requests,
            'tensor_node_hits': cls._node_hits,
            'tensor_materializations': cls._materializations
        }

    @property
    def terms(self):
        eager = self.to_eager()
        return {} if nput.is_numeric(eager) else eager.terms

    @property
    def ndim(self):
        if self._ndim is None:
            self._ndim = self.to_eager().ndim
        return self._ndim

    @property
    def operator_keys(self):
        if self._operator_keys is None:
            kind, *args = self._node
            if kind == 'leaf':
                self._operator_keys = frozenset(
                    coefficient[:2]
                    for product in args[0].terms
                    for coefficient in product
                )
            else:
                self._operator_keys = frozenset().union(*(
                    arg.operator_keys
                    for arg in args
                    if isinstance(arg, type(self))
                ))
        return self._operator_keys

    def __hash__(self):
        if self._hash is None:
            if self._node[0] == 'leaf':
                self._hash = hash(('tensor_leaf', id(self._node[1])))
            else:
                self._hash = hash(self._node)
        return self._hash

    def __eq__(self, other):
        if self is other:
            return True
        if not isinstance(other, type(self)) or self._node[0] != other._node[0]:
            return False
        if self._node[0] == 'leaf':
            return self._node[1] is other._node[1]
        return self._node == other._node

    def __repr__(self):
        return "TensorDAG({}; ndim={})".format(self._node[0], self._ndim)

    def project_coefficients(self, powers, shift=-.5, *, fixed_powers=(), preserve_domains=False):
        options = ((tuple(fixed_powers),bool(preserve_domains))
                   if fixed_powers or preserve_domains else ())
        return type(self)._from_node('occupation_coefficient', (self,tuple(powers),float(shift))+options,
                                    ndim=self._ndim, reduced=False)

    def with_change_domain(self, change, retain):
        return type(self)._from_node('change_domain', (self, tuple(change), bool(retain)),
                                     ndim=self._ndim, reduced=self.reduced)

    def split_change_domains(self):
        """Return shared unrestricted and resonance-residual polynomial DAGs.

        This is the exact algebra P=P_res and Q=1-P_res for nonzero
        changes. Zero-change P has already become the identity in the native
        recurrence. Every residual summand therefore contains a retained
        nonzero change domain; no occupations or resonance groups are bound.
        """
        cls=type(self)
        pending=[(self,False)]
        binary={'mul_simple','rmul_simple','mul_along','rmul_along'}
        def zero(value):return nput.is_numeric(value) and value==0
        def make(node,args):
            kind=node._node[0]
            children=[a for a in args if isinstance(a,cls) or zero(a)]
            if kind=='add':
                if zero(args[0]):return args[1]
                if zero(args[1]):return args[0]
            elif kind in binary:
                if zero(args[0]) or zero(args[1]):return 0
            elif zero(args[0]):return 0
            return cls._from_node(kind,tuple(args),ndim=node._ndim,reduced=node.reduced)
        def add(left,right):
            if zero(left):return right
            if zero(right):return left
            return cls._from_node('add',(left,right),ndim=max(left._ndim,right._ndim))
        while pending:
            node,ready=pending.pop()
            if hasattr(node,'_change_domain_split'):continue
            kind,*args=node._node
            children=[arg for arg in args if isinstance(arg,cls)]
            if not ready and children:
                pending.append((node,True))
                pending.extend((child,False) for child in children if not hasattr(child,'_change_domain_split'))
                continue
            if kind=='leaf':result=node,0
            elif kind=='change_domain':
                child,delta,retain=args
                common,residual=child._change_domain_split
                if retain:
                    # P_res(common+residual), with no unrestricted part.
                    result=0,node
                else:
                    # Q_res(common+residual) = common + Q_res residual - P_res common.
                    kept=0 if zero(residual) else residual.with_change_domain(delta,False)
                    removed=0 if zero(common) else common.with_change_domain(delta,True).scale(-1)
                    result=common,add(kept,removed)
            elif kind in binary:
                left,right=args[:2]
                lu,lr=left._change_domain_split if isinstance(left,cls) else (left,0)
                ru,rr=right._change_domain_split if isinstance(right,cls) else (right,0)
                # delta(A B) = delta(A) B + A_0 delta(B), preserving order,
                # contraction maps, and shared full right-hand subgraphs.
                result=(make(node,[lu,ru]+args[2:]),
                        add(make(node,[lr,right]+args[2:]),make(node,[lu,rr]+args[2:])))
            elif kind=='add':
                left,right=args
                lu,lr=left._change_domain_split;ru,rr=right._change_domain_split
                result=add(lu,ru),add(lr,rr)
            else:
                common_args=[a._change_domain_split[0] if isinstance(a,cls) else a for a in args]
                residual_args=[a._change_domain_split[1] if isinstance(a,cls) else a for a in args]
                result=make(node,common_args),make(node,residual_args)
            node._change_domain_split=result
        return self._change_domain_split


    def change_equations(self):
        """Stream (integer index domains, coefficient key, ordinary equation)."""
        for key, polynomial in PTTensorCoeffProductDAGEvaluationPlan(self).iter_terms():
            yield from _split_change_equations(key, polynomial)

    def format_expr(self):
        if self.has_change_domains:
            return ' + '.join('ChangeSum[{}, {}]'.format(domains,
                PTTensorCoeffProductSum({key: poly}, canonicalize=False).format_expr())
                for domains, key, poly in self.change_equations()) or '0'
        eager = self.to_eager()
        return str(eager) if nput.is_numeric(eager) else eager.format_expr()

    @property
    def has_change_domains(self):
        pending = [(self, False)]
        while pending:
            node, ready = pending.pop()
            if hasattr(node, '_has_change_domains'):
                continue
            children = [child for child in node._node[1:] if isinstance(child, type(self))]
            if ready or node._node[0] == 'change_domain':
                node._has_change_domains = (node._node[0] == 'change_domain' or
                                           any(child._has_change_domains for child in children))
            else:
                pending.append((node, True))
                pending.extend((child, False) for child in children)
        return self._has_change_domains

    @staticmethod
    def _coerce_eager(expression):
        return expression.to_eager() if isinstance(expression, PTTensorCoeffProductDAG) else expression

    def to_eager(self):
        if self._eager is not None:
            return self._eager

        if Diagnostics.enabled:
            type(self)._materializations += 1
        kind, *args = self._node
        if kind == 'occupation_coefficient':
            # Explicit compatibility export streams only the projected terms;
            # never materialize the parent expression as a fallback.
            collected = {}
            for key,poly in PTTensorCoeffProductDAGEvaluationPlan(self).iter_terms():
                collected.setdefault(key,_BalancedPolynomialAccumulator()).append(poly)
            from .Evaluators import _OccupationCoefficientProjection
            collapse = _OccupationCoefficientProjection((),0.,preserve_domains=(
                bool(self._node[-1]) if len(self._node)>4 else False)).compiler()
            terms = {key:value for key,acc in collected.items()
                     if (value := collapse(key,acc.total())) is not None}
            eager = PTTensorCoeffProductSum(terms,
                                           canonicalize=False).combine()
        elif kind == 'change_domain':
            child, change, retain = args
            eager = child.to_eager()
            rows = ((_ChangeDomainIndex(retain, 0),)+(0,)*len(change),
                    (_ChangeDomainIndex(retain, 1),)+change)
            eager = eager.mul_simple(PTEnergyChangeProductSum({rows: 1}))
        elif kind == 'leaf':
            eager = args[0]
        elif kind == 'add':
            left, right = args
            left_eager, right_eager = left.to_eager(), right.to_eager()
            if nput.is_zero(left_eager):
                eager = right_eager
            elif nput.is_zero(right_eager):
                eager = left_eager
            else:
                eager = left_eager + right_eager
        elif kind == 'scale':
            child, scaling = args
            eager = child.to_eager()
            if not nput.is_zero(eager):
                eager = eager.scale(scaling)
        elif kind == 'shift':
            child, shift = args
            eager = child.to_eager()
            if not nput.is_zero(eager):
                eager = eager.shift(shift)
        elif kind == 'shift_energies':
            child, change = args
            eager = child.to_eager()
            if not nput.is_zero(eager):
                eager = eager.shift_energies(change)
        elif kind == 'permute':
            child, permutation, check_perm, allow_padding = args
            eager = child.to_eager()
            if not nput.is_zero(eager):
                eager = eager.permute(
                    permutation,
                    check_perm=check_perm,
                    allow_padding=allow_padding
                )
        elif kind == 'ensure_dimension':
            child, ndim = args
            eager = child.to_eager()
            if not nput.is_zero(eager):
                eager = eager.ensure_dimension(ndim)
        elif kind == 'free_up_indices':
            child, start, stop = args
            eager = child.to_eager()
            if not nput.is_zero(eager):
                eager = eager.free_up_indices(start, stop)
        elif kind == 'flip_energy_terms':
            eager = args[0].to_eager()
            if not nput.is_zero(eager):
                eager = eager.flip_energy_terms()
        elif kind == 'filter_coefficients':
            child, terms, mode = args
            eager = child.to_eager()
            if not nput.is_zero(eager):
                eager = eager.filter_coefficients(terms, mode=mode)
        elif kind == 'filter_energies':
            child, terms, mode = args
            eager = child.to_eager()
            if not nput.is_zero(eager):
                eager = eager.filter_energies(terms, mode=mode)
        elif kind == 'prune_operators':
            child, ops = args
            eager = child.to_eager()
            if not nput.is_zero(eager):
                eager = eager.prune_operators(ops)
        elif kind == 'combine':
            child, combine_coeffs, combine_subterms, combine_energies = args
            eager = child.to_eager()
            if not nput.is_zero(eager):
                eager = eager.combine(
                    combine_coeffs=combine_coeffs,
                    combine_subterms=combine_subterms,
                    combine_energies=combine_energies
                )
        elif kind == 'mul_along':
            (
                left, right, inds, remainder, index_classes,
                mapping, baseline
            ) = args
            left_eager = left.to_eager()
            right_eager = self._coerce_eager(right)
            if nput.is_zero(left_eager) or nput.is_zero(right_eager):
                eager = 0
            else:
                eager = left_eager.mul_along(
                    right_eager,
                    inds,
                    remainder=remainder,
                    index_classes=index_classes,
                    mapping=mapping,
                    baseline=baseline
                )
        elif kind == 'rmul_along':
            left, right, inds, remainder, mapping = args
            left_eager = self._coerce_eager(left)
            right_eager = right.to_eager()
            if nput.is_zero(left_eager) or nput.is_zero(right_eager):
                eager = 0
            else:
                eager = right_eager.rmul_along(
                    left_eager,
                    inds,
                    remainder=remainder,
                    mapping=mapping
                )
        elif kind == 'mul_simple':
            left, right = args
            left_eager = left.to_eager()
            right_eager = self._coerce_eager(right)
            eager = (
                0 if nput.is_zero(left_eager) or nput.is_zero(right_eager) else
                left_eager.mul_simple(right_eager)
            )
        elif kind == 'rmul_simple':
            left, right = args
            left_eager = self._coerce_eager(left)
            right_eager = right.to_eager()
            eager = (
                0 if nput.is_zero(left_eager) or nput.is_zero(right_eager) else
                right_eager.rmul_simple(left_eager)
            )
        else:
            raise ValueError("unknown tensor DAG node {}".format(kind))

        if isinstance(eager, PTTensorCoeffProductDAG):
            eager = eager.to_eager()
        self._eager = eager
        return eager

    def prep_serialization_dict(self):
        return self.to_eager().prep_serialization_dict()

    def mutate(self, terms=default, *, prefactor=default, ndim=default,
               inds_map=default, canonicalize=default, reduced=default):
        if terms is default:
            expression = self
        else:
            expression = type(self).from_sum(PTTensorCoeffProductSum(
                terms,
                prefactor=1 if prefactor is default else prefactor,
                canonicalize=True if canonicalize is default else canonicalize,
                ndim=None if ndim is default else ndim,
                inds_map=None if inds_map is default else inds_map,
                reduced=False if reduced is default else reduced
            ))
        if prefactor is not default and terms is default:
            expression = expression.scale(prefactor)
        return expression

    def audit(self, target=None, required_dimension=None, ignore_constants=True):
        if target is not None and self._ndim is not None and self._ndim < target:
            raise ValueError("too few tensor DAG indices for target {}".format(target))
        return self

    def get_inds(self, key):
        if key not in self._inds_map:
            self._inds_map[key] = self.coeff_product_inds(key)
        return self._inds_map[key]

    def ensure_dimension(self, ndim):
        if self._ndim is not None and self._ndim >= ndim:
            return self
        return type(self)._from_node(
            'ensure_dimension', (self, int(ndim)), ndim=ndim, reduced=self.reduced
        )

    def sort(self):
        return self

    def permute(self, new_inds, check_perm=True, allow_padding=False):
        new_inds = tuple(int(i) for i in new_inds)
        ndim = max(self.ndim, len(new_inds)) if allow_padding else self.ndim
        return type(self)._from_node(
            'permute', (self, new_inds, check_perm, allow_padding), ndim=ndim
        )

    def free_up_indices(self, start, stop):
        return type(self)._from_node(
            'free_up_indices', (self, int(start), int(stop)), ndim=self.ndim
        )

    def shift(self, shift):
        shift = tuple(shift)
        if len(shift) == 0 or not any(shift):
            return self
        return type(self)._from_node('shift', (self, shift), ndim=self.ndim)

    def shift_energies(self, change):
        change = tuple(change)
        if len(change) == 0 or not any(change):
            return self
        return type(self)._from_node('shift_energies', (self, change), ndim=self.ndim)

    def scale(self, scaling):
        if nput.is_numeric(scaling):
            if scaling == 0:
                return 0
            if scaling == 1:
                return self
        return type(self)._from_node('scale', (self, scaling), ndim=self.ndim)

    def flip_energy_terms(self):
        return type(self)._from_node('flip_energy_terms', (self,), ndim=self.ndim)

    def exclude_degenerate_changes(self, contexts, permutations, num_fixed,
                                   nmodes, statistics=None, _change_support=None):
        """Prune fixed-mode resonances before streaming this request's DAG.

        Axis sets describe permutations still possible inside a product. A
        denominator is decided here only when its transformed change is fixed
        under all those permutations. Free-axis changes remain on the normal
        assignment-mask path. Neither this DAG nor its polynomial leaves mutate.
        The returned digest identifies actual removals, independent of the
        physical row, so equal pruned expressions can share compiled entries.
        """
        memo = {}
        masks = {}
        zero_signature = b'zero'

        def signature(value):
            return hashlib.sha256(repr(value).encode('utf-8')).digest()

        def dimension(value):
            if isinstance(value, type(self)):
                if value._ndim is None:
                    raise ValueError('tensor DAG resonance pruning requires axis metadata')
                return value._ndim
            return value.ndim

        def project(delta, axes, offset):
            if offset is None:
                return None
            groups = []
            for axis, targets in enumerate(axes):
                value = delta[axis] if axis < len(delta) else 0
                if targets is None:
                    if value:
                        return None
                    continue
                targets = set(targets)
                values = [value]
                separate = []
                for existing, coefficients in groups:
                    if targets.intersection(existing):
                        targets.update(existing)
                        values.extend(coefficients)
                    else:
                        separate.append((existing, coefficients))
                groups = separate + [(targets, values)]
            if any(delta[len(axes):]):
                return None
            result = list(offset)
            for targets, values in groups:
                if len(targets) != len(values) or any(v != values[0] for v in values):
                    if any(values):
                        return None
                    continue
                for axis in targets:
                    result[axis] += values[0]
            return tuple(result)

        def energy_sum(value, frame):
            if _change_support is not None:
                return value, b''
            axes, offset, flipped = frame
            terms = {}
            for key, poly in value.terms.items():
                resonant = np.zeros(len(permutations), dtype=bool)
                for side, delta in value.side_change_iter(key):
                    if side not in (0, 1):
                        continue
                    fixed = project(delta, axes, offset)
                    if fixed is None:
                        continue
                    side = int(side) ^ flipped
                    mask_key = (side, fixed)
                    mask = masks.get(mask_key)
                    if mask is None:
                        plan = contexts[side]._bind_energy_key(
                            ((0,) + fixed,), nmodes, 'linear'
                        )[0]
                        mask = (np.zeros(len(permutations), dtype=bool)
                                if plan is None else plan.evaluate(permutations))
                        masks[mask_key] = mask
                        if statistics is not None:
                            statistics['fixed_denominators_identified'] = statistics.get('fixed_denominators_identified', 0) + 1
                    resonant |= mask
                    if np.all(resonant):
                        break
                if np.all(resonant):
                    if statistics is not None:
                        statistics['energy_keys_pruned'] = statistics.get('energy_keys_pruned', 0) + 1
                else:
                    terms[key] = poly
            if len(terms) == len(value.terms):
                return value, b''
            if not terms:
                return 0, zero_signature
            return value.mutate(terms, canonicalize=False), signature(tuple(terms))

        def at(axes, index):
            return axes[index] if index < len(axes) else None

        def product_frames(args, frame):
            left, right, inds, remainder, classes, mapping, baseline = args
            axes, offset, flipped = frame
            contact = len(inds[0])
            result = []
            for side, child in enumerate((left, right)):
                size = max(dimension(child), max((*inds[side], *remainder[side]), default=-1) + 1)
                mapped = [None] * size
                for index, source in enumerate(inds[side]):
                    mapped[source] = at(axes, index)
                start = contact + (len(remainder[0]) if side else 0)
                # The native product permutes equivalent proper remainder
                # axes when floating indices attach. Carry that entire orbit.
                labels = (classes[side] if classes is not None else ())
                groups = {}
                for index, source in enumerate(remainder[side]):
                    label = labels[source] if len(labels) > source else source
                    groups.setdefault(label, []).append((source, start + index))
                for group in groups.values():
                    targets = [at(axes, target) for _, target in group]
                    union = (None if any(t is None for t in targets) else
                             tuple(sorted({i for t in targets for i in t})))
                    for source, _ in group:
                        mapped[source] = union
                result.append((tuple(mapped), offset, flipped))
            return result

        def visit_value(value, frame):
            if isinstance(value, type(self)):
                return visit(value, frame)
            if isinstance(value, PTEnergyChangeProductSum):
                return energy_sum(value, frame)
            return value, b''

        def visit(node, frame):
            key = (node, frame)
            if key in memo:
                return memo[key]
            if statistics is not None:
                statistics['nodes_visited'] = statistics.get('nodes_visited', 0) + 1
            kind, *args = node._node
            axes, offset, flipped = frame
            if kind == 'change_domain':
                if _change_support is None:
                    raise ValueError('change-index domains require their evaluation bindings')
                child, delta, retain = args
                if retain and any(delta):
                    # Multiplication identifies aliases across children, never
                    # within one child. A retained change must therefore have
                    # the same nonzero quantum pattern as an actual block edge,
                    # even while its physical axes are still floating.
                    if 'retained_patterns' not in masks:
                        bindings=_change_support[1]
                        masks['retained_patterns']={tuple(sorted(int(x) for x in np.asarray(v)-u if x))
                            for group in bindings[0] for u in map(np.asarray,group) for v in group
                            if not np.array_equal(u,v)}
                    if tuple(sorted(int(x) for x in delta if x)) not in masks['retained_patterns']:
                        result=0,signature(('excluded_change_pattern',tuple(delta)))
                        if statistics is not None:
                            statistics['change_domains_pattern_pruned']=statistics.get('change_domains_pattern_pruned',0)+1
                        memo[key]=result
                        return result
                start = project((), axes, offset)
                end = project(delta, axes, offset)
                mask = None
                if start is not None and end is not None:
                    support_key = (start, end, retain)
                    if support_key not in masks:
                        rows, bindings, baseline = _change_support
                        labels = _filtered_block_labels(bindings[0])
                        decisions = []
                        for states, perms in rows:
                            for source in np.atleast_2d(states):
                                for perm in np.asarray(perms)[:, :num_fixed]:
                                    origin, target = source.copy(), source.copy()
                                    origin[perm] += start
                                    target[perm] += end
                                    inside = (np.array_equal(origin, target) or
                                        (tuple(origin) in labels and labels[tuple(origin)] == labels.get(tuple(target))))
                                    decisions.append(bool(np.all(origin >= 0) and np.all(target >= 0)
                                                          and (inside if retain else not inside)))
                        masks[support_key] = np.asarray(decisions, dtype=bool)
                        if statistics is not None:
                            statistics['change_domains_identified'] = statistics.get('change_domains_identified', 0) + 1
                    mask = masks[support_key]
                if mask is not None and len(mask) and not np.any(mask):
                    result = 0, signature(('excluded_change_domain', start, end, retain))
                    if statistics is not None:
                        statistics['change_domains_pruned'] = statistics.get('change_domains_pruned', 0) + 1
                else:
                    value, changed = visit(child, frame)
                    if changed and nput.is_zero(value):
                        result = 0, zero_signature
                    elif mask is not None and len(mask) and np.all(mask):
                        result = value, signature(('accepted_change_domain', start, end, retain, changed))
                        if statistics is not None:
                            statistics['change_domains_resolved'] = statistics.get('change_domains_resolved', 0) + 1
                    elif changed:
                        result = type(self)._from_node(kind, (value, delta, retain),
                            ndim=node._ndim, reduced=node.reduced), signature((kind, changed))
                    else:
                        result = node, b''
            elif kind == 'leaf':
                original = args[0]
                terms = {}
                changes = []
                for coefficient, poly in original.terms.items():
                    value, changed = visit_value(poly, frame)
                    if changed:
                        changes.append((coefficient, changed))
                    if not nput.is_zero(value):
                        terms[coefficient] = value
                if not changes:
                    result = node, b''
                elif not terms:
                    result = 0, zero_signature
                else:
                    result = type(self).from_sum(original.mutate(
                        terms, canonicalize=False, ndim=original.ndim
                    )), signature(('leaf', tuple(changes)))
            else:
                frames = [frame] * len(args)
                if kind == 'permute':
                    permutation = args[1]
                    if sorted(permutation) == list(range(len(permutation))):
                        inverse = np.argsort(permutation)
                        mapped = tuple(at(axes, int(inverse[i])) if i < len(inverse)
                                       else at(axes, i) for i in range(dimension(args[0])))
                    else:
                        mapped = (None,) * dimension(args[0])
                    frames[0] = (mapped, offset, flipped)
                elif kind == 'free_up_indices':
                    mapped = tuple(at(axes, i) if i < args[1] else None
                                   for i in range(dimension(args[0])))
                    frames[0] = (mapped, offset, flipped)
                elif kind == 'shift_energies' or kind == 'shift' and _change_support is not None:
                    frames[0] = (axes, project(args[1], axes, offset), flipped)
                elif kind == 'flip_energy_terms':
                    frames[0] = (axes, offset, flipped ^ 1)
                elif kind == 'mul_along':
                    if isinstance(args[1], type(self)):
                        frames[:2] = product_frames(args, frame)
                        # Right-state shifts on floating product axes may have
                        # multiple assignments. Leave those domains to the
                        # indexed selector rather than guessing an origin.
                        if _change_support is not None and args[-1] is not None and any(args[-1]):
                            frames[1] = (frames[1][0], None, flipped)
                    else:
                        inds, remainder = args[2:4]
                        ordered = tuple(inds[0]) + tuple(remainder[0])
                        inverse = np.argsort(ordered)
                        mapped = tuple(at(axes, int(inverse[i])) if i < len(inverse)
                                       else None
                                       for i in range(dimension(args[0])))
                        frames[0] = (mapped, offset, flipped)
                        if isinstance(args[1], PTEnergyChangeProductSum):
                            mapped = [None] * dimension(args[1])
                            for index, source in enumerate(inds[1]):
                                mapped[source] = at(axes, index)
                            frames[1] = (tuple(mapped), offset, flipped)
                elif kind == 'rmul_along':
                    frames[:2] = [((None,) * (dimension(arg) if isinstance(arg, type(self))
                                            else getattr(arg, 'ndim', 0)), offset, flipped)
                                  for arg in args[:2]]
                elif kind not in {'add', 'scale', 'shift', 'ensure_dimension',
                                  'filter_coefficients', 'filter_energies',
                                  'prune_operators', 'combine', 'mul_simple',
                                  'rmul_simple', 'occupation_coefficient'}:
                    raise NotImplementedError('resonance pruning of tensor DAG node {}'.format(kind))
                updated = list(args)
                changes = [b''] * len(args)
                # Resolve the standalone denominator before visiting a large
                # numerator subtree that it may annihilate.
                order = list(range(len(args)))
                if kind in {'mul_simple', 'rmul_simple', 'mul_along', 'rmul_along'}:
                    order = [i for i in range(2) if not isinstance(args[i], type(self))] + [
                        i for i in range(2) if isinstance(args[i], type(self))]
                eliminated = False
                for index in order:
                    updated[index], changes[index] = visit_value(args[index], frames[index])
                    if changes[index] and nput.is_zero(updated[index]) and kind != 'add':
                        eliminated = True
                        if statistics is not None:
                            statistics['zero_nodes'] = statistics.get('zero_nodes', 0) + 1
                        break
                if eliminated:
                    result = 0, zero_signature
                elif not any(changes):
                    result = node, b''
                else:
                    digest = signature((kind, tuple(changes)))
                    if kind == 'add' and nput.is_zero(updated[0]):
                        result = updated[1], digest
                    elif kind == 'add' and nput.is_zero(updated[1]):
                        result = updated[0], digest
                    else:
                        result = type(self)._from_node(kind, tuple(updated),
                            ndim=node._ndim, reduced=node.reduced), digest
            memo[key] = result
            return result

        axes = tuple((i,) for i in range(num_fixed)) + (None,) * max(0, dimension(self) - num_fixed)
        offset = (0,) * num_fixed
        if _change_support is not None and _change_support[2] is not None:
            baseline = tuple(_change_support[2])
            offset = (None if any(baseline[num_fixed:]) else
                      (baseline + (0,) * num_fixed)[:num_fixed])
        return visit(self, (axes, offset, 0))

    def select_change_domains(self, state_perms, bindings, num_fixed,
                              nmodes, baseline=None, statistics=None):
        """Bind integer summation support before expanding any polynomial.

        Only scopes whose origin and change are fixed for every assignment
        can be resolved here. Floating-index scopes retain the same integer
        domain and are selected by the native indexed evaluator.
        """
        return self.exclude_degenerate_changes(None, (), num_fixed, nmodes,
            statistics=statistics, _change_support=(state_perms, bindings, baseline))


    def _map_leaves(self, method, *args, **kwargs):
        # Map each node of the DAG once. A recursive walk visits a shared
        # subgraph once per path to it, which is exponential in the depth of
        # the sharing, and it also hits the recursion limit on deep DAGs.
        cls = type(self)
        mapped = {}
        pending = [(self, False)]
        while pending:
            node, ready = pending.pop()
            if id(node) in mapped:
                continue
            kind, *node_args = node._node
            if kind == 'leaf':
                transformed = getattr(node_args[0], method)(*args, **kwargs)
                mapped[id(node)] = cls.from_sum(transformed)
            elif not ready:
                pending.append((node, True))
                pending.extend(
                    (arg, False) for arg in node_args
                    if isinstance(arg, cls) and id(arg) not in mapped
                )
            else:
                mapped[id(node)] = cls._from_node(
                    kind,
                    tuple(mapped[id(arg)] if isinstance(arg, cls) else arg
                          for arg in node_args),
                    ndim=node._ndim, reduced=node.reduced
                )
        return mapped[id(self)]

    def filter_coefficients(self, terms, mode='match'):
        return type(self)._from_node(
            'filter_coefficients', (self, tuple(terms), mode), ndim=self.ndim
        )

    def filter_energies(self, terms, mode='match'):
        return type(self)._from_node(
            'filter_energies', (self, tuple(terms), mode), ndim=self.ndim
        )

    def prune_operators(self, ops):
        return self._map_leaves('prune_operators', tuple(ops))

    def combine(self, combine_coeffs=False, combine_subterms=True, combine_energies=False):
        if self.reduced and not combine_coeffs and combine_subterms and not combine_energies:
            return self
        return type(self)._from_node(
            'combine',
            (self, combine_coeffs, combine_subterms, combine_energies),
            ndim=self.ndim,
            reduced=True
        )

    def mul_along(self, other, inds, remainder=None, index_classes=None,
                  mapping=None, baseline=None):
        if isinstance(other, PTTensorCoeffProductSum):
            other = type(self).from_sum(other)
            ndim = self.ndim + other.ndim - len(inds[0])
        else:
            ndim = self.ndim
        return type(self)._from_node(
            'mul_along',
            (
                self, other, inds, remainder, index_classes,
                mapping, baseline
            ),
            ndim=ndim
        )

    def rmul_along(self, other, inds, remainder=None, mapping=None):
        if isinstance(other, PTTensorCoeffProductSum):
            return type(self).from_sum(other).mul_along(
                self, inds, remainder=remainder, mapping=mapping
            )
        return type(self)._from_node(
            'rmul_along', (other, self, inds, remainder, mapping), ndim=self.ndim
        )

    def mul_simple(self, other):
        if isinstance(other, PTTensorCoeffProductSum):
            other = type(self).from_sum(other)
            ndim = self.ndim + other.ndim
        else:
            ndim = self.ndim
        return type(self)._from_node('mul_simple', (self, other), ndim=ndim)

    def rmul_simple(self, other):
        if isinstance(other, PTTensorCoeffProductSum):
            return type(self).from_sum(other).mul_simple(self)
        return type(self)._from_node('rmul_simple', (other, self), ndim=self.ndim)

    def __add__(self, other):
        if nput.is_numeric(other):
            if other == 0:
                return self
            raise NotImplementedError("non-zero scalar tensor addition")
        other = type(self).from_sum(other)
        return type(self)._from_node(
            'add', (self, other), ndim=max(self.ndim, other.ndim)
        )

    def __radd__(self, other):
        return self.__add__(other)

    def __mul__(self, other):
        return self.scale(other) if nput.is_numeric(other) else self.mul_simple(other)

    def __rmul__(self, other):
        return self.__mul__(other)



class _HashableNumPyIndex:
    """Immutable, hash-once identity for a small NumPy index array.

    The compact C-order byte string preserves exact shape/dtype/content
    equality without retaining a view's potentially much larger base array.
    This is intended for integer quantum-state and permutation indices; it is
    deliberately an exact binary identity rather than a numerical comparison.
    """

    __slots__ = ('shape', 'dtype', 'data', '_hash')

    def __init__(self, values):
        values = np.asanyarray(values)
        self.shape = values.shape
        self.dtype = values.dtype.str
        self.data = values.tobytes(order='C')
        self._hash = hash((self.shape, self.dtype, self.data))

    @property
    def nbytes(self):
        return len(self.data)

    def tolist(self):
        return np.frombuffer(
            self.data,
            dtype=np.dtype(self.dtype)
        ).reshape(self.shape).tolist()

    def memory_size(self):
        return (
            sys.getsizeof(self)
            + sys.getsizeof(self.shape)
            + sys.getsizeof(self.dtype)
            + sys.getsizeof(self.data)
        )

    def __hash__(self):
        return self._hash

    def __eq__(self, other):
        return (
            self is other
            or isinstance(other, type(self))
            and self._hash == other._hash
            and self.shape == other.shape
            and self.dtype == other.dtype
            and self.data == other.data
        )


class _StatePermutationBlockIdentity:
    """Hash-once identity for the concrete quantum states in one block."""

    __slots__ = ('states', '_hash', '_size')

    def __init__(self, state_keys):
        self.states = tuple(tuple(block) for block in state_keys)
        self._hash = hash(self.states)
        self._size = sys.getsizeof(self) + sys.getsizeof(self.states) + sum(
            sys.getsizeof(block) + sum(
                key.memory_size() if isinstance(key,_HashableNumPyIndex)
                else sys.getsizeof(key) + len(key)*sys.getsizeof(0)
                for key in block)
            for block in self.states)

    def memory_size(self):
        return self._size

    def __hash__(self):
        return self._hash

    def __eq__(self, other):
        return (
            self is other
            or isinstance(other, type(self)) and self.states == other.states
        )



class _MaterializedEvaluationCache:
    """Two-level bounded cache for a materialized numerical evaluation plan.

    ``exact_cache`` holds complete polynomial/state results. A bounded exact
    cache uses flat keys so nested dictionaries cannot bypass its byte limit.
    The caller may own this cache for one batch or a complete immutable-input
    evaluation plan; concrete state keys and bounded values permit that reuse
    without accumulating unbounded arrays across successive calls.
    """

    class PathNodeView:
        """Fast root-local cache layered over reusable canonical DAG nodes."""

        _noncanonical_kinds = {'add', 'remap', 'permutation_sum'}

        def __init__(self, shared):
            self.local = {}
            self.shared = shared

        @classmethod
        def _is_shareable(cls, key):
            path = key[0]
            if isinstance(path, _PolynomialProgramTable): return True
            return (
                path._node is None
                or path._node[0] not in cls._noncanonical_kinds
            )

        def get(self, key, default=None):
            value = self.local.get(key, BoundedCache._missing)
            if value is not BoundedCache._missing:
                return value
            if self._is_shareable(key):
                value = self.shared.get(key, BoundedCache._missing)
                if value is not BoundedCache._missing:
                    self.local[key] = value
                    return value
            return default

        def __setitem__(self, key, value):
            self.local[key] = value
            if self._is_shareable(key):
                self.shared[key] = value

    def __init__(self, exact_cache, max_items=100000, max_bytes=64 * 1024 ** 2,
                 batch_cache=None):
        self.exact_cache = exact_cache
        self.batch_cache = (
            BoundedCache.sampled(max_items, max_bytes)
            if batch_cache is None else batch_cache
        )
        self.shared_path_node_cache = self.batch_cache.view('polypath_nodes')
        self.path_axis_cache = self.batch_cache.view('polypath_axes')
        self.polynomial_block_cache = self.batch_cache.view('polynomial_blocks')
        self.exact_hits = 0
        self.exact_misses = 0
        self.block_hits = 0
        self.block_misses = 0

    def lookup_polynomial_block(self, key):
        value = self.polynomial_block_cache.get(
            key, BoundedCache._missing
        )
        if Diagnostics.enabled:
            if value is BoundedCache._missing:
                self.block_misses += 1
            else:
                self.block_hits += 1
        return value

    def store_polynomial_block(self, key, value):
        self.polynomial_block_cache[key] = value

    def new_path_node_cache(self):
        return self.PathNodeView(self.shared_path_node_cache)

    def stats(self):
        return {
            **self.batch_cache.stats(),
            'exact_cache_hits': self.exact_hits,
            'exact_cache_misses': self.exact_misses,
            'block_cache_hits': self.block_hits,
            'block_cache_misses': self.block_misses
        }



class _BalancedPolynomialAccumulator:
    """Add native polynomials in a balanced tree without changing any terms.

    One block per occupied binary level bounds the number of partial sums.
    Equal-sized blocks merge as they arrive; final blocks retain stream order.
    This avoids quadratic dictionary copying in large symbolic DAG chunks.
    """
    def __init__(self):
        self.blocks = []

    def append(self, value):
        level = 0
        while level < len(self.blocks) and self.blocks[level] is not None:
            value = self.blocks[level] + value
            self.blocks[level] = None
            level += 1
        if level == len(self.blocks):
            self.blocks.append(value)
        else:
            self.blocks[level] = value

    def total(self):
        result = None
        for value in reversed(self.blocks):
            if value is not None:
                result = value if result is None else result + value
        return 0 if result is None else result


class PTTensorCoeffProductDAGEvaluationPlan:
    """Streams canonical tensor terms from a DAG without expanding it globally.

    Small shared subgraphs are replayed from a bounded, batch-local cache.  A
    node that emits more than ``max_cached_terms`` is deliberately not retained;
    this prevents a high-order expression from replacing global materialization
    with an equally unbounded evaluation cache.

    Developer note
    --------------
    This is an intentionally conservative bridge to the existing evaluator: it
    streams symbolic terms in bounded chunks, but each ``mul_along`` still uses
    the legacy symbolic product machinery.  Consequently it bounds peak storage
    without removing the dominant high-order Cartesian-product work.

    A genuinely faster evaluator should compile this DAG to numerical blocks
    whose axes describe open/fixed mode indices.  Leaf blocks should substitute
    coefficient tensors and evaluate their ``PolyPath`` values once; ``add`` and
    ``scale`` then operate directly on arrays, while ``permute``, ``shift``, and
    ``free_up_indices`` transform block metadata/views.  Most importantly,
    ``mul_along`` should become an indexed relational join/tensor contraction of
    its child blocks.  That join must preserve equality/distinctness constraints
    between mode indices, energy-denominator metadata, and final square-root
    factors; multiplying already-reduced child scalars is not equivalent.

    Keep unsupported nodes on this streaming implementation until their direct
    kernels have numerical parity tests.  The materialized evaluator must also
    remain available as the reference/timing implementation.
    """

    coalesce_combine_nodes = False
    combine_chunk_size = 32768
    progress_callback = None

    def __init__(self, expression, cache=None, max_cached_terms=256, max_mode_indices=None):
        if not isinstance(expression, PTTensorCoeffProductDAG):
            raise TypeError(type(expression))
        self.expression = expression
        self.cache = cache if cache is not None else BoundedCache.sampled()
        self.max_mode_indices = max_mode_indices
        self.term_cache = self.cache.view(('tensor_terms',max_mode_indices,self.coalesce_combine_nodes))
        self.operation_cache = self.cache.view(('tensor_operations',max_mode_indices,self.coalesce_combine_nodes))
        self.max_cached_terms = max(0, int(max_cached_terms))
        self.nodes_visited = 0
        self.node_cache_hits = 0
        self.terms_emitted = 0
        self.mode_terms_pruned = 0
        self.coalesce_input_terms = self.coalesce_output_terms = self.coalesce_batches = 0

    def _coalesced_combine(self, child, combine_coeffs, combine_subterms, combine_energies):
        pending = {}
        count = 0
        def flush():
            compact = PTTensorCoeffProductSum(
                {key:value.total() for key,value in pending.items()},
                canonicalize=False, reduced=False
            ).combine(combine_coeffs=combine_coeffs, combine_subterms=combine_subterms,
                      combine_energies=combine_energies)
            result = self._result_terms(compact)
            if Diagnostics.enabled:
                self.coalesce_output_terms += len(result)
                self.coalesce_batches += 1
            pending.clear()
            if self.progress_callback is not None:
                self.progress_callback(self)
            return result
        for key,poly in self._iter_node(child):
            accumulator = pending.get(key)
            if accumulator is None:
                accumulator = pending[key] = _BalancedPolynomialAccumulator()
            accumulator.append(poly)
            if Diagnostics.enabled:
                self.coalesce_input_terms += 1
            count += 1
            if count >= self.combine_chunk_size:
                yield from flush()
                count = 0
        if pending:
            yield from flush()

    @staticmethod
    def _singleton(key, poly):
        return PTTensorCoeffProductSum({key: poly}, canonicalize=False)

    def _result_terms(self,result):
        if nput.is_zero(result):
            return ()
        if isinstance(result, PTTensorCoeffProductDAG):
            result = result.to_eager()
        if not isinstance(result, PTTensorCoeffProductSum):
            raise TypeError("tensor DAG operation produced {}".format(type(result)))
        terms = result.terms.items()
        if self.max_mode_indices is None:
            return tuple(terms)
        kept = []
        for key,poly in terms:
            # Distinct aliases within one child cannot subsequently merge:
            # products identify aliases across children, never within one.
            count = len({index for coefficient in key for index in coefficient[2:]})
            if count <= self.max_mode_indices:
                kept.append((key,poly))
            else:
                if Diagnostics.enabled:
                    self.mode_terms_pruned += 1
        return tuple(kept)

    def _cached_operation(self, key, operation):
        result = self.operation_cache.get(key, BoundedCache._missing)
        if result is BoundedCache._missing:
            result = tuple(operation())
            if len(result) <= self.max_cached_terms:
                self.operation_cache[key] = result
        return result

    def _unary(self, method, key, poly, *args, **kwargs):
        op_key = (method, key, poly, args, tuple(sorted(kwargs.items())))
        return self._cached_operation(
            op_key,
            lambda: self._result_terms(
                getattr(self._singleton(key, poly), method)(*args, **kwargs)
            )
        )

    def _binary(self, method, left, right, *args, **kwargs):
        left_key, left_poly = left
        op_key = (method, left_key, left_poly, right, args, tuple(sorted(kwargs.items())))

        def apply():
            left_sum = self._singleton(left_key, left_poly)
            right_value = right
            if isinstance(right_value, tuple) and len(right_value) == 2 and isinstance(right_value[0], tuple):
                right_value = self._singleton(*right_value)
            return self._result_terms(getattr(left_sum, method)(right_value, *args, **kwargs))

        return self._cached_operation(op_key, apply)

    def _iter_node(self, node):
        cached = self.term_cache.get(node, BoundedCache._missing)
        if cached is not BoundedCache._missing:
            if Diagnostics.enabled:
                self.node_cache_hits += 1
            yield from cached
            return

        if Diagnostics.enabled:
            self.nodes_visited += 1
        buffered = []
        overflowed = False
        for term in self._iter_node_uncached(node):
            if not overflowed:
                if len(buffered) < self.max_cached_terms:
                    buffered.append(term)
                else:
                    overflowed = True
                    buffered = None
            yield term
        if not overflowed:
            self.term_cache[node] = tuple(buffered)

    def _iter_node_uncached(self, node):
        kind, *args = node._node
        if kind == 'occupation_coefficient':
            from .Evaluators import _OccupationCoefficientProjection
            child,powers,shift,*options = args
            transform = _OccupationCoefficientProjection(powers,shift,*options).compiler()
            for key,poly in self._iter_node(child):
                projected = transform(key,poly)
                if projected is not None:
                    yield key,projected
        elif kind == 'change_domain':
            child, change, retain = args
            rows = ((_ChangeDomainIndex(retain, 0),)+(0,)*len(change),
                    (_ChangeDomainIndex(retain, 1),)+change)
            factor = PTEnergyChangeProductSum({rows: 1})
            for key, poly in self._iter_node(child):
                yield from self._binary('mul_simple', (key, poly), factor)
        elif kind == 'leaf':
            yield from self._result_terms(args[0])
        elif kind == 'add':
            yield from self._iter_node(args[0])
            yield from self._iter_node(args[1])
        elif kind in {
            'scale', 'shift', 'shift_energies', 'ensure_dimension',
            'free_up_indices', 'flip_energy_terms', 'filter_coefficients',
            'filter_energies', 'prune_operators'
        }:
            child, *op_args = args
            for key, poly in self._iter_node(child):
                yield from self._unary(kind, key, poly, *op_args)
        elif kind == 'permute':
            child, permutation, check_perm, allow_padding = args
            for key, poly in self._iter_node(child):
                yield from self._unary(
                    'permute', key, poly, permutation,
                    check_perm=check_perm, allow_padding=allow_padding
                )
        elif kind == 'combine':
            if self.coalesce_combine_nodes:
                yield from self._coalesced_combine(*args)
            else:
                yield from self._iter_node(args[0])
        elif kind in {'mul_along', 'mul_simple'}:
            left, right, *op_args = args
            if isinstance(right, PTTensorCoeffProductDAG):
                for left_term in self._iter_node(left):
                    for right_term in self._iter_node(right):
                        yield from self._binary(kind, left_term, right_term, *op_args)
            else:
                for left_term in self._iter_node(left):
                    yield from self._binary(kind, left_term, right, *op_args)
        elif kind in {'rmul_along', 'rmul_simple'}:
            left, right, *op_args = args
            if isinstance(left, PTTensorCoeffProductDAG):
                for left_term in self._iter_node(left):
                    for right_term in self._iter_node(right):
                        yield from self._binary(
                            kind[1:], right_term, left_term, *op_args
                        )
            else:
                for key, poly in self._iter_node(right):
                    yield from self._unary(kind, key, poly, left, *op_args)
        else:
            raise NotImplementedError("direct tensor DAG evaluation of {}".format(kind))

    def iter_terms(self):
        for term in self._iter_node(self.expression):
            if Diagnostics.enabled:
                self.terms_emitted += 1
            yield term

    def stats(self):
        return {
            'dag_nodes_visited': self.nodes_visited,
            'dag_node_cache_hits': self.node_cache_hits,
            'dag_terms_emitted': self.terms_emitted,
            'dag_terms_pruned_for_modes': self.mode_terms_pruned,
            'dag_coalesce_input_terms': self.coalesce_input_terms,
            'dag_coalesce_output_terms': self.coalesce_output_terms,
            'dag_coalesce_batches': self.coalesce_batches,
            **self.cache.stats()
        }


class SqrtChangePoly(PolynomialInterface):
    def __init__(self, poly_obj:'PolynomialInterface', change, shift, canonicalize=False):
        if nput.is_numeric(poly_obj): raise ValueError("{} isn't a polynomial".format(poly_obj))
        if isinstance(poly_obj, SqrtChangePoly): raise ValueError("sqrt changes shouldn't stack...")
        self.poly_obj = poly_obj

        self.poly_change = np.asanyarray(change, dtype=int)
        # if len(self.poly_change) > 0 and np.count_nonzero(self.poly_change) != len(self.poly_change):
        #     raise ValueError("no zero changes allowed ({})".format(
        #         self.poly_change)
        #     )
        if shift is None: shift = np.zeros(len(self.poly_change), dtype=int)
        self.shift_start = np.asanyarray(shift)
        # if len(self.shift_start) != len(self.poly_change): raise ValueError(self.shift_start, self.poly_change)
        if len(self.shift_start) < len(self.poly_change):
            self.shift_start = np.pad(self.shift_start, [0, len(self.poly_change) - len(self.shift_start)])
        if canonicalize:
            self.poly_obj, self.poly_change, self.shift_start = self.canonicalize(
                self.poly_obj, self.poly_change, self.shift_start
            )

    def split_change_domains(self):
        """Reusable symbolic unrestricted and retained-resonance polynomials."""
        common,residual=PTTensorCoeffProductDAG.from_sum(self.poly_obj).split_change_domains()
        return tuple(0 if nput.is_zero(poly) else self.mutate(poly_obj=poly)
                     for poly in (common,residual))

    def project_coefficients(self, powers, shift=-.5, *, fixed_powers=(), preserve_domains=False):
        if np.any(self.shift_start) or (np.any(self.poly_change) and len(fixed_powers) != len(self.poly_change)):
            raise ValueError('coefficient projection requires an unshifted expression and powers for every changed mode')
        if any(int(p) != p or p < 0 for p in tuple(fixed_powers)+tuple(powers)):
            raise ValueError('coefficient powers must be nonnegative integers')
        dag = PTTensorCoeffProductDAG.from_sum(self.poly_obj)
        return self.mutate(poly_obj=dag.project_coefficients(powers,shift,
                           fixed_powers=fixed_powers,preserve_domains=preserve_domains))

    def mutate(self,
               poly_obj:'Any'=default,
               poly_change:'Any'=default,
               shift_start:'Any'=default,
               canonicalize=False
               ):
        return type(self)(
            self.poly_obj if poly_obj is default else poly_obj,
            self.poly_change if poly_change is default else poly_change,
            self.shift_start if shift_start is default else shift_start,
            canonicalize=canonicalize
        )

    def change_equations(self):
        """Export polynomial equations with explicit integer-change domains."""
        if isinstance(self.poly_obj, PTTensorCoeffProductDAG):
            stream = self.poly_obj.change_equations()
        else:
            stream = (row for key, poly in self.poly_obj.terms.items()
                      for row in _split_change_equations(key, poly))
        for domains, key, polynomial in stream:
            yield domains, self.mutate(poly_obj=PTTensorCoeffProductSum(
                {key: polynomial}, canonicalize=False))

    def to_state(self, serializer=None):

        is_dag = isinstance(self.poly_obj, PTTensorCoeffProductDAG)
        projected = (not is_dag and any(
            isinstance(poly, PTEnergyChangeProductSum) and (
                () in poly.terms or any(
                    getattr(row[0], '_state_projection_factor', False)
                    for key in poly.terms for row in key))
            for poly in self.poly_obj.terms.values()
        ))
        if is_dag or projected:
            # The checkpoint is already handled through the serializer's
            # pseudo-pickle protocol.  Keeping the DAG as one byte buffer also
            # preserves cross-node sharing and avoids eager tensor expansion.
            # Eager projected keys also need this protocol: flattening their
            # integer-compatible tags into a numeric tree loses the projector.
            dag_buffer = np.frombuffer(
                pickle.dumps(self.poly_obj, protocol=5), dtype=np.uint8
            ).copy()
            return {
                ('tensor_dag_v1' if is_dag else 'projected_tensor_v1'): dag_buffer,
                'change': list(self.poly_change),
                'shift': list(self.shift_start)
            }

        serial_dict = self.poly_obj.prep_serialization_dict()

        ordered_terms = []
        for k,t in serial_dict.items():
            ndim = len(PTTensorCoeffProductSum.coeff_product_inds(k))
            if len(ordered_terms) < ndim:
                ordered_terms.extend({} for _ in range(ndim - len(ordered_terms)))
            ordered_terms[ndim-1][k] = t
        # print(ordered_terms)


        poly_trees = []
        for i,d in enumerate(ordered_terms):
            if len(d) > 0:
                # print("=~!~~+~+~", i)
                poly_trees.append(TreeSerializer.serialize_tree_dict(d))
            else:
                poly_trees.append(None)

        return {
            'term_ordered_trees': poly_trees,
            'change': list(self.poly_change),
            'shift': list(self.shift_start)
        }

    @classmethod
    def from_state(cls, state, serializer=None):
        if 'tensor_dag_v1' in state or 'projected_tensor_v1' in state:
            poly_obj = pickle.loads(
                np.asanyarray(state.get('tensor_dag_v1', state.get('projected_tensor_v1')),
                              dtype=np.uint8).tobytes()
            )
            return cls(
                poly_obj,
                state['change'],
                state['shift'],
                canonicalize=False
            )
        trees = state['term_ordered_trees']
        full_dict = {}
        for subtree_data in trees:
            if subtree_data is not None:
                keys, shape, vals = subtree_data
                # print(keys)
                deserialized_dict = TreeSerializer.rebuild_tree_dict(keys, shape, vals)
                for subk, subv in deserialized_dict.items():
                    full_dict[subk] = subv
        # print(full_dict)
        poly_obj = PTTensorCoeffProductSum.from_serialization_dict(full_dict)
        # print(poly_obj.format_expr())
        return cls(
            poly_obj,
            state['change'],
            state['shift'],
            canonicalize=False
        )

    def strip(self):
        stripped = [i for i,c in enumerate(self.poly_change) if c == 0]
        if len(stripped) == 0: return self
        if not np.equal(
            stripped,
            np.arange(stripped[0], stripped[0] + len(stripped))
        ).all():
            raise ValueError("bad strip...")
        return self.mutate(
            self.poly_obj.free_up_indices(stripped[0], stripped[-1]+1),
            poly_change=[c for c in self.poly_change if c != 0],
            shift_start=[s for c,s in zip(self.poly_change, self.shift_start) if c != 0]
        )

    @property
    def ndim(self):
        return len(self.poly_change)
    def audit(self, target=None, ignore_constants=True):
        try:
            if isinstance(self.poly_obj, PTTensorCoeffProductSum):
                return self.poly_obj.audit(
                    target=self.ndim,
                    required_dimension=[abs(c) for c in self.poly_change],
                    ignore_constants=ignore_constants
                )
            else:
                return self.poly_obj.audit(
                    target=self.ndim,
                    ignore_constants=ignore_constants
                )
        except ValueError:
            raise ValueError("bad poly found in {}".format(self.short_repr()))
    def ensure_dimension(self, ndim):
        if self.ndim >= ndim: return self
        poly_change = np.pad(self.poly_change, [0, ndim - len(self.poly_change)])
        shift_start = np.pad(self.shift_start, [0, ndim - len(self.shift_start)])

        return self.mutate(
            self.poly_obj.ensure_dimension(ndim),
            poly_change,
            shift_start,
            canonicalize=False
        )

    def __repr__(self):
        return "{}({}, {})".format("SqS", self.poly_change, self.poly_obj)
    def short_repr(self):
        return "{}({}, <{}>)".format("SqS", self.poly_change, type(self.poly_obj).__name__)
    def format_sqrt_expr(self):
        subexprs = []
        for i, (s,c) in enumerate(zip(self.shift_start, self.poly_change)):
            if c == 0:
                subexprs.append("")
                continue
            if s == 0:
                subexprs.append(
                    "S[n_{idx},{delta}]".format(idx=i, delta=c, shift=s)
                )
            else:
                subexprs.append(
                    "S[n_{idx}{sign}{shift},{delta}]".format(idx=i, delta=c, sign="-" if s < 0 else "+", shift=abs(s))
                )
        return "".join(subexprs)
    def format_expr(self):
        """
        Formats in a Mathematica-ingestible format
        :return:
        """
        base_expr = self.poly_obj.format_expr()
        sqrt_exprs = self.format_sqrt_expr()
        if len(sqrt_exprs) == 0:
            return base_expr
        elif "+" in base_expr:
            return "{}*({})".format(sqrt_exprs,
                                    base_expr.replace("\n", "\n" + (" "*(len(sqrt_exprs) + 2)))
                                    )
        else:
            return "{}*{}".format(sqrt_exprs, base_expr)

    def __add__(self, other):
        if nput.is_numeric(other):
            if other == 0:
                return self
            else:
                raise NotImplementedError(...)
                return self + ProductPTPolynomial([[other]])
        elif isinstance(other, SqrtChangePoly):
            self, other = self.align_dimensions(other)
            if np.any(other.poly_change != self.poly_change):
                raise ValueError("can't add polys with different changes {} vs. {}".format(
                    self.poly_change,
                    other.poly_change
                ))
            if np.any(self.shift_start != other.shift_start):
                raise ValueError("can't add polys with different starts {} vs. {}".format(
                    self.shift_start,
                    other.shift_start
                ))
            new_poly = self.poly_obj + other.poly_obj
            if nput.is_numeric(new_poly):
                return new_poly
            else:
                return self.mutate(new_poly)
        else:
            raise NotImplementedError("not sure how to add {} and {}".format(self, other))

    def combine(self, combine_coeffs=False, combine_subterms=True, combine_energies=False):
        if nput.is_numeric(self.poly_obj): return self.poly_obj
        return self.mutate(
            self.poly_obj.combine(
                combine_coeffs=combine_coeffs, combine_subterms=combine_subterms,
                combine_energies=combine_energies
            )
        )
    def sort(self):
        if nput.is_numeric(self.poly_obj): return self.poly_obj
        return self.mutate(
            self.poly_obj.sort()
        )
    # def condense(self, inds=None, return_inds=False, check_inds=True):
    #     if inds is None: inds = np.arange(len(self.poly_change))
    #     inds = np.asanyarray(inds)
    #     self.poly_change = np.asanyarray(self.poly_change)
    #     if check_inds:
    #         inds = inds[self.poly_change[inds] == 0]
    #     condensed, new_obj = self.poly_obj.condense(inds, return_inds=True)
    #     new_change = np.delete(self.poly_change, inds)
    #     new_shift = np.delete(self.shift_start, inds)
    #     new_poly = self.mutate(
    #         new_obj,
    #         new_change,
    #         new_shift
    #     )
    #     if return_inds:
    #         return condensed, new_poly
    #     else:
    #         return new_poly

    @classmethod
    def canonicalize(cls, poly_obj, poly_change, shift_start):
        pc = np.asanyarray(poly_change)
        nonzero_pos = np.nonzero(poly_change)
        if len(nonzero_pos) == 0 or len(nonzero_pos[0]) == 0:
            return poly_obj, np.array([], dtype=int), np.array([], dtype=int)
        nonzero_pos = nonzero_pos[0]

        # lex_sort = np.lexsort([
        #     np.sign(pc[nonzero_pos]),
        #     -np.abs(pc[nonzero_pos])
        # ])
        nzp = pc[nonzero_pos]
        lex_sort = PerturbationTheoryTerm.change_sort(nzp)
        sorting = nonzero_pos[lex_sort]
        # sorting = nonzero_pos[np.argsort(pc[nonzero_pos])]
        if (
                len(nonzero_pos) == len(poly_change)
                and np.all(sorting == np.arange(len(poly_change)))
        ):
            return poly_obj, poly_change, shift_start
        else:
            return poly_obj.permute(sorting), np.take(poly_change, sorting), np.take(shift_start, sorting)

    def shift_energies(self, changes):
        return self.mutate(self.poly_obj.shift_energies(changes))
    def shift(self, shift, shift_energies=False):
        # new_shift, new_poly = self.(self.poly_change, shift)
        # if not nput.is_numeric(new_poly):
        #     new_poly = self.poly_obj.mul_simple(new_poly)
        og = self.shift_start
        if len(shift) < len(og):
            shift = list(shift) + [0] * (len(og) - len(shift))
        elif len(shift) > len(og):
            og = list(og) + [0] * (len(shift) - len(og))

        return self.mutate(
            poly_obj=self.poly_obj if not shift_energies else self.poly_obj.shift_energies(shift),
            shift_start=[i+s for i,s in zip(og, shift)]
        )
    def scale(self, scaling):
        if nput.is_numeric(scaling) and scaling == 1: return self
        return self.mutate(
            self.poly_obj.scale(scaling)
        )
    def permute(self, perm, check_perm=True, allow_padding=False):

        if check_perm:
            if np.any(np.sort(perm) != np.arange(len(perm))):
                raise ValueError("bad permutation {}")

        new = self.mutate(
            self.poly_obj.permute(perm, check_perm=check_perm, allow_padding=allow_padding),
            [self.poly_change[c] for c in perm],
            [
                self.shift_start[c] for c in perm
                if c < len(self.poly_change)
            ]
        )
        if _DEBUG_AUDIT:
            try:
                new.audit()
            except ValueError:
                raise ValueError("bad permutation {} for {}".format(
                    perm, self
                ))
        return new

    def filter_coefficients(self, terms, mode='match'):
        return self.mutate(self.poly_obj.filter_coefficients(terms, mode=mode))
    def filter_energies(self, terms, mode='match'):
        return self.mutate(self.poly_obj.filter_energies(terms, mode=mode))

    def __radd__(self, other):
        return self.__add__(other)
    def __mul__(self, other):
        raise NotImplementedError(...)
    def __rmul__(self, other):
        raise NotImplementedError(...)

    def get_change_poly(self,
                        initial_changes, extra_changes,
                        initial_shift, extra_shift
                        ):
        for ich, ish, ech, esh in zip(
                initial_changes, initial_shift,
                extra_changes, extra_shift
        ):
            if esh != ich + ish:
                raise ValueError("split shifts not supported; target shift/changes: {}/{}, initial shift/changes: {}/{}".format(
                    extra_shift, extra_changes, initial_shift, initial_changes
                ))
        new_changes = [i_sh + e_sh for i_sh, e_sh in zip(initial_changes, extra_changes)]
        sqrt_contrib = HarmonicOscillatorRaisingLoweringPolyTerms.get_direction_change_poly(
            extra_changes, initial_changes, initial_shift
        )
        if len(sqrt_contrib) == 0: sqrt_contrib = None
        if sqrt_contrib is not None:
            sqrt_contrib = make_product_polynomial(
                sqrt_contrib,
                steps=0,
                representation='path' if polynomial_uses_path(self.poly_obj) else 'eager'
            )

        return new_changes, sqrt_contrib

    def mul_simple(self, other:'PolynomialInterface'):
        if isinstance(other, SqrtChangePoly):
            return self.mul_along(other, list(range(self.ndim)))
            # return self.mutate(
            #     self.poly_obj.mul_simple(other.poly_obj),
            #     np.concatenate([self.poly_change, other.poly_change]),
            #     np.concatenate([self.shift_start, other.shift_start])
            # )
        else:
            return self.mutate(
                self.poly_obj.mul_simple(other),
                # np.pad(self.poly_change, [0, other.ndim]),
                self.poly_change,
                # np.pad(self.shift_start, [0, other.ndim])
                self.shift_start
            )
    def rmul_simple(self, other: 'PolynomialInterface') -> 'PolynomialInterface':
        if isinstance(other, SqrtChangePoly):
            other.mul_simple(self)
        else:
            return self.mutate(
                other.mul_simple(self.poly_obj),
                self.poly_change,
                # np.pad(self.poly_change, [other.ndim, 0]),
                self.shift_start
                # np.pad(self.shift_start, [other.ndim, 0])
            )
    def mul_along(self, other, inds, remainder=None, mapping=None, baseline=None):
        if isinstance(other, SqrtChangePoly):
            # self, other = self.align_dimensions(other)
            # if len(self.poly_change) != len(other.poly_change): raise ValueError(
            #     "mismatch between poly changes {} vs. {}".format(
            #         self.poly_change, other.poly_change
            #     )
            # )
            # if len(self.shift_start) != len(other.shift_start): raise ValueError(
            #     "mismatch between poly shifts {} vs. {}".format(
            #         self.shift_start, other.shift_start
            #     )
            # )

            (self_inds, other_inds) = inds
            (self_remainder, other_remainder) = remainder
            #     = ProductPTPolynomial.get_index_mapping(
            #     self.ndim,
            #     other.ndim,
            #     inds,
            #     return_remainder=True
            # )
            # baseline = self._get_mul_baseline()


            self_shift = [
                self.shift_start[c] if len(self.shift_start) > c else 0
                for c in self_inds
            ]
            other_shift = [
                other.shift_start[c] if len(other.shift_start) > c else 0
                for c in other_inds
            ]
            self_change =[
                self.poly_change[c] if len(self.poly_change) > c else 0
                for c in self_inds
            ]
            other_change = [
                other.poly_change[c] if len(other.poly_change) > c else 0
                for c in other_inds
            ]

            new_changes, sqrt_contrib = self.get_change_poly(
                self_change, other_change,
                self_shift, other_shift
            )

            final_changes = [
                                n for n in new_changes #if n != 0
                            ] + [
                                self.poly_change[r] for r in self_remainder
                                # if r < len(self.poly_change)
                            ] + [
                                other.poly_change[r] for r in other_remainder
                                # if r < len(other.poly_change)
                            ]

            final_shift = [
                              s for s, n in zip(self_shift, new_changes) #if n != 0
                          ] + [
                              self.shift_start[r] for r in self_remainder
                              # if r < len(self.shift_start)
                          ] + [
                              other.shift_start[r] for r in other_remainder
                              # if r < len(other.shift_start)
                          ]

            other_poly = other.poly_obj
            if other_shift != self_shift:
                # need to shift other to align with self
                delta_shift = [o-s for s,o in zip(self_shift, other_shift)]
                total_shift = [0] * len(other.shift_start)
                for d,c in zip(delta_shift, other_inds):
                    total_shift[c] = d
                other_poly = other_poly.shift(total_shift)
            if sqrt_contrib is not None:
                new_self = self.poly_obj.mul_along(sqrt_contrib,
                                                   [self_inds, np.arange(len(self_inds))],
                                                   remainder=[self_remainder, np.arange(len(self_inds), sqrt_contrib.ndim)]
                                                   ) #TODO: optimize out remainder calc.
                if baseline is not None:
                    baseline = [
                        baseline[x] for x in itertools.chain(self_inds, self_remainder)
                    ] + list(
                        baseline[len(self_inds) + len(self_remainder):]
                    )
                # inds and remainders must be adapted for the fact that we've moved around indices in doing this...


                new_poly = new_self.mul_along(
                    other_poly,
                    [list(range(len(self_inds))), other_inds],
                    index_classes=[
                        [
                            self.poly_change[i] for i in itertools.chain(self_inds, self_remainder)
                        ], #TODO: check if the shift needs to be included...
                        other.poly_change
                    ],
                    remainder=[
                        np.arange(len(self_inds), len(self_remainder) + len(self_inds)),
                        other_remainder
                    ],
                    baseline=baseline
                )
            else:
                # -print("md...", len(final_changes))
                new_poly = self.poly_obj.mul_along(
                    other_poly,
                    [self_inds, other_inds],
                    index_classes=[self.poly_change, other.poly_change],
                    remainder=[self_remainder, other_remainder],
                    baseline=baseline
                )

            return self.mutate(
                new_poly,
                final_changes,
                final_shift
            )
        else:
            # raise NotImplementedError(type(other))
            return self.mutate(
                self.poly_obj.mul_along(other, inds, remainder=remainder, mapping=mapping),
                self.poly_change,
                self.shift_start
            )
    def rmul_along(self, other, inds, remainder=None, mapping=None):
        if isinstance(other, SqrtChangePoly):
            return other.mul_along(self, inds, remainder=remainder, mapping=mapping)
        else:
            raise NotImplementedError(type(other))
            return type(self)(
                self.poly_obj.rmul_along(other, inds, remainder=remainder, mapping=mapping),
                self.poly_change,
                self.shift_start
            )


def make_product_polynomial(coeffs, prefactor=1, idx=None, steps=None,
                            representation='eager'):
    if representation == 'path':
        return PolyPath.from_coeffs(coeffs, prefactor=prefactor, idx=idx, steps=steps)
    return ProductPTPolynomial(coeffs, prefactor=prefactor, idx=idx, steps=steps)


def polynomial_uses_path(poly):
    if isinstance(poly, PolyPath):
        return True
    if isinstance(poly, PTTensorCoeffProductDAG):
        return True
    if isinstance(poly, SqrtChangePoly):
        return polynomial_uses_path(poly.poly_obj)
    if isinstance(poly, (PTEnergyChangeProductSum, PTTensorCoeffProductSum)):
        return any(polynomial_uses_path(p) for p in poly.terms.values())
    return False


def coerce_polynomial_representation(poly, representation):
    """Convert polynomial leaves while retaining tensor/energy/change wrappers."""
    if nput.is_numeric(poly):
        return poly
    if isinstance(poly, PTTensorCoeffProductDAG):
        if representation == 'path':
            return poly
        poly = poly.to_eager()
    if isinstance(poly, SqrtChangePoly):
        return poly.mutate(
            coerce_polynomial_representation(poly.poly_obj, representation)
        )
    if isinstance(poly, PTTensorCoeffProductSum):
        converted = poly.mutate(
            {
                key: coerce_polynomial_representation(subpoly, representation)
                for key, subpoly in poly.terms.items()
            },
            canonicalize=False
        )
        return (
            PTTensorCoeffProductDAG.from_sum(converted)
                if representation == 'path' else
            converted
        )
    if isinstance(poly, PTEnergyChangeProductSum):
        return poly.mutate(
            {
                key: coerce_polynomial_representation(subpoly, representation)
                for key, subpoly in poly.terms.items()
            },
            canonicalize=False
        )
    if representation == 'path':
        return PolyPath.from_polynomial(poly)
    elif isinstance(poly, PolyPath):
        return poly.to_eager()
    return poly

class PerturbationTheoryTerm(metaclass=abc.ABCMeta):
    """
    A generic version of one of the three terms in
    PT that will generate a correction polynomial
    """

    use_intermediate_normalization = False
    def __init__(self, logger=None, checkpoint=None,
                 allowed_terms=None,
                 allowed_energy_changes=None,
                 intermediate_normalization=None,
                 allowed_coefficients=None,
                 disallowed_coefficients=None,
                 polynomial_representation='path'):
        self._cascade_manager = CascadingPropertyManager(self)
        self._exprs = None
        self._raw_changes = {}
        self._changes = None
        self._cache = {}
        self.logger = Logger.lookup(logger)
        self.checkpoint = Checkpointer.build_canonical(checkpoint)
        self.intermediate_normalization = (
            self.use_intermediate_normalization
                if intermediate_normalization is None else
            intermediate_normalization
        )
        self.allowed_terms = allowed_terms
        self.allowed_energy_changes = allowed_energy_changes
        self.allowed_coefficients=allowed_coefficients
        self.disallowed_coefficients=disallowed_coefficients
        self.polynomial_representation = polynomial_representation

    def get_subexpressions(self) -> 'Iterable[PerturbationTheoryTerm]':
        raise NotImplementedError("just here to be overloaded")
    def __mul__(self, other):
        if isinstance(other, PerturbationTheoryTerm):
            return PerturbationTheoryTermProduct.lookup(self, other)
        else:
            return ScaledPerturbationTheoryTerm.lookup(self, other)
    def __rmul__(self, other):
        # if isinstance(other, PerturbationTheoryTerm): # can't happen
        #     return ...
        return ScaledPerturbationTheoryTerm.lookup(self, other)
    def __neg__(self):
        return ScaledPerturbationTheoryTerm.lookup(self, -1)

    @property
    def expressions(self):
        if self._exprs is None:
            self._exprs = self.get_subexpressions()
        return self._exprs

    @staticmethod
    def change_sort_key(changes):
        changes = np.asanyarray(changes)
        return -(10 * np.abs(changes) - (changes < 0))
    @classmethod
    def change_sort(cls, changes):
        return np.argsort(cls.change_sort_key(changes), axis=-1)
    @classmethod
    def sorted_changes(cls, changes):
        changes = np.asanyarray(changes)
        sorting = np.argsort(cls.change_sort_key(changes), axis=-1)
        if sorting.ndim == 1:
            return tuple(changes[i] for i in sorting)
        else:
            return nput.vector_take(changes, sorting, shared=sorting.ndim - 1)
    @abc.abstractmethod
    def get_changes(self) -> 'dict[tuple[int], Any]':
        ...
    @derived_property(depends_on=['change_rules'], settable=False)
    def changes(self):
        return self.get_changes()

    def get_serializer_key(self): # to be overridden
        return None
    @property
    def serializer_key(self):
        key = self.get_serializer_key()
        parent = getattr(self,'parent',None)
        if key is not None and any(getattr(term,'combine_kinetic_coefficients',False)
                                   for term in getattr(parent,'hamiltonian_expansion',())):
            key += '::effective-vertices-v1'
        if key is not None and getattr(getattr(self, 'parent', None), 'filtered_block_recurrence', False):
            key += '::filtered-block-v1::intermediate={}'.format(bool(self.intermediate_normalization))
        if key is not None and getattr(getattr(self, 'parent', None), 'resonance_independent', False):
            key += '::unrestricted-action-v1'
        if key is not None and self.allowed_terms is not None:
            # Keep legacy unfiltered checkpoints intact and isolate each term.
            terms = tuple(sorted(set(tuple(t) for t in self.allowed_terms)))
            key = '{}::terms={}'.format(key,terms)
        return key

    debug_logger = None
    _logger_stack = []
    @contextlib.contextmanager
    def debug_logging(self):
        try:
            PerturbationTheoryTerm._logger_stack.append(PerturbationTheoryTerm.debug_logger)
            PerturbationTheoryTerm.debug_logger = self.logger
            yield PerturbationTheoryTerm.debug_logger
        finally:
            PerturbationTheoryTerm.debug_logger = PerturbationTheoryTerm._logger_stack.pop()
    @classmethod
    def default_logger(cls) -> Logger:
        return Logger.lookup(cls.debug_logger)


    _default_filters = None
    _filter_stack = []
    @contextlib.contextmanager
    def default_filtering(self):
        try:
            PerturbationTheoryTerm._filter_stack.append(PerturbationTheoryTerm._default_filters)
            PerturbationTheoryTerm._default_filters = [
                self.allowed_terms,
                self.allowed_energy_changes,
                self.allowed_coefficients,
                self.disallowed_coefficients
            ]
            yield PerturbationTheoryTerm._default_filters
        finally:
            PerturbationTheoryTerm._default_filters = PerturbationTheoryTerm._filter_stack.pop()
    @classmethod
    def default_filters(cls) -> Logger:
        return Logger.lookup(cls.default_filters)

    def get_core_poly(self, changes, shift=None):
        terms = (term.get_poly_terms(changes, shift=shift)
                 for term in self.expressions)
        if self.polynomial_representation == 'path':
            # A left-folded sum gives the lazy tensor DAG linear depth.
            accumulator = _BalancedPolynomialAccumulator()
            for term in terms:
                accumulator.append(term)
            return accumulator.total()
        return sum(terms)

    simplify_by_default = True
    caching_enabled = True
    def get_poly_terms(self, changes, simplify=None, shift=None) -> 'SqrtChangePoly':
        if any(c == 0 for c in changes): raise ValueError(...)

        if simplify is None: simplify = self.simplify_by_default
        log_level = Logger.LogLevel.Normal if _DEBUG_PRINT else Logger.LogLevel.MoreDebug
        with self.logger.block(
                tag="{s}[{c}]", s=self, c=format_quantum_change(changes),
                log_level=log_level):
            og_changes = changes
            with self.debug_logging():
                change_perm = self.change_sort(np.asanyarray(changes))
                changes = tuple(changes[p] for p in change_perm)  # for c in changes if c != 0)
                if shift is not None:
                    pad_change_perm = list(change_perm) + list(range(len(change_perm), len(shift)))
                    shift = [shift[s] for s in pad_change_perm]
                cache_key = changes if shift is None or all(s == 0 for s in shift) else (changes, tuple(shift))
                terms = self.changes.get(cache_key, None)
                if not isinstance(terms, SqrtChangePoly):
                    if nput.is_numeric(terms): return terms
                    subterms = self.changes.get(changes, 0)
                    if nput.is_numeric(subterms): return subterms

                    start = time.time()
                    if self.caching_enabled:
                        chk = self.checkpoint
                        key = self.serializer_key
                        substr = "<"+str(changes).replace("(", "").replace(")", "").replace(",", "|")+">"
                        try:
                            if key is None:
                                terms = None
                            else:
                                terms = chk[key, substr]
                        except KeyError:
                            terms = None
                        if terms is None:
                            terms = self.get_core_poly(changes, shift=None)
                            if simplify and not nput.is_zero(terms):
                                terms = terms.combine()  # .sort()
                            self.changes[changes] = terms
                            if key is not None:
                                # try:
                                #     base = chk[key]
                                # except KeyError:
                                #     base = {}
                                # base[substr] = terms
                                chk[key, substr] = terms
                    else:
                        terms = self.get_core_poly(changes, shift=None)
                        if simplify and not nput.is_zero(terms):
                            terms = terms.combine()

                    if self.polynomial_representation == 'path' and not nput.is_numeric(terms):
                        terms = coerce_polynomial_representation(terms, 'path')
                        self.changes[changes] = terms

                    if nput.is_numeric(terms): return terms

                    if shift is not None and len(shift) > 0 and not nput.is_numeric(terms):
                        terms = terms.shift(shift)

                    if isinstance(terms, SqrtChangePoly):
                        if [c for c in terms.poly_change] != [c for c in changes]:
                            raise ValueError("[[{}]] change mismatch between expected {} and obtained {}".format(
                                self,
                                changes,
                                terms.poly_change
                            ))
                        if shift is None:
                            shift = []  # * len(changes)
                        test_shift = terms.shift_start
                        if all(s == 0 for s in test_shift): test_shift = []
                        spec_shift = shift
                        if all(s == 0 for s in spec_shift): spec_shift = []
                        if list(test_shift) != list(spec_shift):
                            raise ValueError("[[ {} ]] shift mismatch between expected {} and obtained {}".format(
                                self,
                                shift,
                                terms.shift_start
                            ))

                    if cache_key is not changes:
                        self.changes[cache_key] = terms

                    if self.disallowed_coefficients is not None:
                        terms = terms.filter_coefficients(self.disallowed_coefficients, mode='exclude')
                    if self.allowed_coefficients is not None:
                        terms = terms.filter_coefficients(self.allowed_coefficients, mode='include')
                    if self.allowed_energy_changes is not None:
                        terms = terms.filter_energies(self.allowed_energy_changes, mode='include')

                    self.logger.log_print('{s}[{c}]: {p}', p=terms, s=self,
                                          c=format_quantum_change(og_changes),
                                          preformatter=lambda **vars: dict(vars, p=vars['p'].format_expr()),
                                          log_level=log_level
                                          )
                    end = time.time()
                    # self.logger.log_print('terms: {t}', t=terms)
                    self.logger.log_print('took {e:.3f}s', e=end - start, log_level=log_level)
                terms = terms.permute(np.argsort(change_perm))

        return terms

    def __call__(self, changes, shift=None, coeffs=None, freqs=None,
                 check_sorting=True, simplify=None, return_evaluator=True):
        from .Evaluators import current_compilation
        compilation = current_compilation()
        if (return_evaluator and coeffs is None and freqs is None and compilation is not None
                and compilation[0].using_artifacts):
            # A transferred compiled library can satisfy this formal request
            # before any polynomial is constructed. Retain the ordinary builder
            # as a lazy source for a missing entry or an explicit symbolic use.
            changes = tuple(changes)
            result = PerturbationTheoryExpressionEvaluator(self, None, changes, logger=self.logger,
                expression_factory=lambda: PerturbationTheoryTerm.__call__(self, changes,
                    shift=shift, check_sorting=check_sorting, simplify=simplify, return_evaluator=False))
            result._expression_shift = None if shift is None else tuple(shift)
            return result
        with self.logger.block(tag="Building evaluator {term}[{change}]",
                term=self, change=format_quantum_change(changes)):
            with self.checkpoint: # Open once
                if coeffs is not None:
                    raise NotImplementedError(...)
                if freqs is not None:
                    raise NotImplementedError(...)
                perm = None
                if check_sorting:
                    perm = self.change_sort(np.asanyarray(changes))
                    changes = tuple(changes[p] for p in perm)
                change_key = changes if (shift is None or all(s==0 for s in shift)) else (changes, shift)
                start = time.time()
                terms = self.changes.get(change_key, None)
                if not isinstance(terms, SqrtChangePoly):
                    if not nput.is_numeric(terms):
                        terms = self.changes.get(changes, 0)
                        if not (isinstance(terms, SqrtChangePoly) or nput.is_numeric(terms)):
                            terms = self.get_poly_terms(changes, shift=None)
                            if simplify and not nput.is_zero(terms):
                                terms = terms.combine()#.sort()
                            self.changes[changes] = terms
                        if shift is not None:
                            raise NotImplementedError("reshifting not supported yet...")
                            self.changes[change_key] = terms
                if perm is not None and not nput.is_numeric(terms):
                    terms = terms.permute(np.argsort(perm))
                end = time.time()
                # self.logger.log_print('terms: {t}', t=terms)
                self.logger.log_print('took {e:.3f}s', e=end - start)
                if return_evaluator:
                    return PerturbationTheoryExpressionEvaluator(self, terms, changes,
                                                                 logger=self.logger)
                else:
                    return terms

class OperatorExpansionTerm(PerturbationTheoryTerm):
    def __init__(self, terms, order=None, identities=None, symmetrizers=None, index=None, allowed_terms=None, change_rules=None, **opts):
        super().__init__(allowed_terms=allowed_terms, **opts)

        self.terms = terms
        self.order = order
        self._change_rules = change_rules
        if index is None:
            index = order
        self.index = index
        self.identities = np.arange(len(self.terms)) if identities is None else identities
        if symmetrizers is None:
            symmetrizers = [PTTensorCoeffProductSum._potential_symmetrizer] * (max(self.identities) + 1)
        self.symmetrizers = symmetrizers

        self._grouped_terms = {}
        for i,t in enumerate(self.terms):
            self._grouped_terms[len(t)] = self._grouped_terms.get(len(t), [])
            self._grouped_terms[len(t)].append([i, tuple(t)])

        self.term_sizes = set(self._grouped_terms.keys())
        self._change_poly_cache = {}
        self._cache = {}

    def __repr__(self):
        return "M[{}]".format(self.order)

    def __fingerprint__(self):
        # Construction/evaluation populate polynomial and change caches. Their
        # hashes and object identities are not part of an operator definition.
        names = ('terms','order','index','identities','symmetrizers','_change_rules',
                 'allowed_terms','allowed_coefficients','disallowed_coefficients',
                 'allowed_energy_changes','intermediate_normalization','polynomial_representation')
        return (type(self).__module__,type(self).__qualname__,
                {name:getattr(self,name) for name in names})

    _generator_cache = {}
    @classmethod
    def _get_generator(cls, term_list, inds):
        ts = tuple(term_list[i] for i in inds)
        if ts not in cls._generator_cache:
            cls._generator_cache[ts] = HarmonicOscillatorMatrixGenerator(ts) # no reason not to cache these...
        return cls._generator_cache[ts]


    @classmethod
    def _evaluate_poly_coeffs(cls, terms, inds, delta, shift, sqrt_scale=False):
        gen = cls._get_generator(terms, inds)
        poly_contrib = gen.poly_coeffs(delta, shift=shift, sqrt_scale=sqrt_scale)
        if nput.is_numeric(poly_contrib) and poly_contrib == 0:
            return 0
        # now we get any remainder we need to
        sqrt_contrib = HarmonicOscillatorRaisingLoweringPolyTerms.get_direction_change_poly(delta, shift)
        if sqrt_contrib is not None:
            poly_contrib = scipy.signal.convolve(poly_contrib, sqrt_contrib)
        return poly_contrib

    @cascading_property
    def change_rules(self):
        return self._change_rules
    @change_rules.setter
    def change_rules(self, rules):
        self._change_rules = rules
    def get_changes(self) -> 'dict[tuple[int], Any]':
        if self.change_rules is not None:
            return {
                tuple(r): None
                for r in self.change_rules
            }
        else:
            changes = {}
            l_check = set()
            for t in self.terms:
                start = len(t) % 2
                for tl in range(start, len(t)+1, 2): # I could really just get the max_len but w/e
                    if tl == 0:
                        changes[()] = None
                        continue
                    if tl not in l_check: # just to not redo it...
                        l_check.add(tl)
                        for p in IntegerPartitioner.partitions(tl, pad=False):
                            l = p.shape[1]
                            for n in range(0, l+1):
                                base_perm = [1] * (l-n) + [-1] * n
                                signs = UniquePermutations(base_perm).permutations()
                                resigned_perms = (signs[np.newaxis] * p[:, np.newaxis]).reshape(-1, l)
                                sort_key = self.change_sort(resigned_perms)
                                sort_perms = nput.vector_take(resigned_perms, sort_key, shared=1)
                                for subp in sort_perms:
                                    changes[tuple(subp)] = None
        return changes

    def get_subexpressions(self) -> 'Iterable[PerturbationTheoryTerm]':
        raise NotImplementedError("shouldn't need this here...")

    def get_core_poly(self, changes, shift=None) -> 'SqrtChangePoly':  # TODO: CACHE THIS SHIT

        poly_contribs = {}

        abs_change = np.array([abs(c) for c in changes])
        total_change = sum(abs_change)

        og_changes = changes
        og_shift = shift
        og_abs_change = abs_change
        og_dim = np.count_nonzero(og_abs_change)

        for group_size, terms in self._grouped_terms.items():
            remainder = group_size - total_change
            if remainder < 0 or remainder % 2 != 0:
                continue # literally cannot execute this change

            # not sure what to do if group_size == 0...?
            if group_size == 0: # constant contrib
                for term_index, term_list in terms:
                    subpoly = make_product_polynomial(
                        [], steps=0,
                        representation=self.polynomial_representation
                    )
                    prefactor = (self.index, self.identities[term_index])  # type: tuple[int]
                    poly_contribs[(prefactor,)] = poly_contribs.get((prefactor,), 0) + subpoly
                continue

            total_dim = og_dim + remainder // 2

            # now we enumerate all of the ways to partition `group_size` elements across `total_dim` dimensions
            # such that each `partition - changes` is even in each dimension
            abs_change = np.pad(og_abs_change, [0, total_dim - len(og_abs_change)])
            changes = np.pad(og_changes, [0, total_dim - len(og_changes)]).astype(int)
            if shift is None:
                shift = np.zeros(len(changes), dtype=int)
            else:
                if len(shift) < total_dim:
                    shift = np.pad(shift, [0, total_dim - len(shift)]).astype(int)

            targets = SymmetricGroupGenerator(total_dim).get_terms(group_size)#.partitions(group_size, max_len=total_dim, pad=True)
            check_changes = targets - abs_change[np.newaxis, :]
            good_targets = np.all(np.logical_and(check_changes >= 0, check_changes % 2 == 0), axis=1)
            # bad_targets = np.full(len(targets), True, dtype=bool)
            # bad_pos = np.arange(len(targets))
            # for cc in itertools.permutations(abs_change):
            #     check_changes = targets[bad_targets] - np.array(cc)[np.newaxis, :]
            #     good_targets = np.all(np.logical_and(check_changes >= 0, check_changes % 2 == 0), axis=1)
            #     bad_targets[bad_pos[bad_targets][good_targets]] = False
            # good_targets = np.logical_not(bad_targets)
            targets = targets[good_targets]

            # now we drop all targets with zeros in front of non-zero terms b.c. these
            # terms should reduce out
            neg_pos = np.where(np.diff(targets, axis=1) > 0)
            zero_pos = np.where(targets[neg_pos] == 0)
            if len(zero_pos) > 0:
                zero_pos = zero_pos[0]
                if len(zero_pos) > 0:
                    bad_targs = neg_pos[0][zero_pos]
                    targets = targets[ProductPTPolynomial.fast_ind_remainder(len(targets), bad_targs)]

            # then for each of these partition sizes, we enumerate the actual partitions of the terms
            # and then take the corresponding terms to generate partitions
            # (and multiply by the appropriate tensor coefficients)

            partitioner = UniquePartitions(np.arange(group_size))
            ckey = tuple(changes)
            skey = tuple(shift)
            # Integer
            symmetric_terms = [
                (i,t) for i,t in terms if all(letter == 'x' for letter in t)
                and self.symmetrizers[self.identities[i]] is PTTensorCoeffProductSum._potential_symmetrizer
            ]
            general_terms = [(i,t) for i,t in terms if (i,t) not in symmetric_terms]
            for partition_sizes in targets:
                if symmetric_terms:
                    # Identical coordinate letters and a fully symmetric force
                    # tensor make every assignment with these sizes identical.
                    # Count them exactly instead of enumerating set partitions.
                    sizes = tuple(partition_sizes)
                    p_vec = np.split(np.arange(group_size),np.cumsum(partition_sizes)[:-1])
                    tensor_idx = tuple(np.repeat(np.arange(total_dim),partition_sizes))
                    multiplicity = UniquePermutations.count_permutations(partition_sizes)
                    for term_index,term_list in symmetric_terms:
                        base_poly = self._resolve_poly(term_list,sizes,0,ckey,skey,p_vec,changes,shift)
                        if nput.is_zero(base_poly):continue
                        if isinstance(base_poly,SqrtChangePoly):base_poly = base_poly.poly_obj
                        prefactor = (self.index,self.identities[term_index])+tensor_idx
                        poly_contribs[(prefactor,)] = (poly_contribs.get((prefactor,),0)
                                                       + base_poly.scale(multiplicity))
                if not general_terms:continue
                parts, inv_splits = partitioner.partitions(partition_sizes,
                                                          return_partitions=False, take_unique=False,
                                                          return_indices=True, return_inverse=True
                                                          )

                inv_splits = np.concatenate(inv_splits, axis=1)
                partition_sizes = tuple(partition_sizes)
                base_tup = np.array(sum(([i]*s for i,s in enumerate(partition_sizes)), []))
                for j,p_vec in enumerate(zip(*parts)):
                    tensor_idx = tuple(base_tup[inv_splits[j]])
                    for term_index,term_list in general_terms:
                        base_poly = self._resolve_poly(term_list, partition_sizes, j, ckey, skey, p_vec, changes, shift)
                        if nput.is_zero(base_poly): continue

                        term_id = self.identities[term_index]
                        if self.symmetrizers[term_id] is not None:
                            symm_idx = self.symmetrizers[term_id](tensor_idx)
                        else:
                            symm_idx = tensor_idx

                        if isinstance(base_poly, SqrtChangePoly):
                            base_poly = base_poly.poly_obj
                        prefactor = (self.index, term_id) + symm_idx #type: tuple[int]
                        poly_contribs[(prefactor,)] = poly_contribs.get((prefactor,), 0) + base_poly

        new = PTTensorCoeffProductSum(poly_contribs, canonicalize=False)
        if og_shift is None: og_shift = [0] * len(og_changes)
        if len(og_shift) < len(og_changes): og_shift = np.pad(og_shift, [0, len(og_changes) - len(og_shift)])
        # new.audit()
        return SqrtChangePoly(new, og_changes, og_shift)

    keep_ints = True
    _poly_cache = {}
    def _resolve_poly(self, term_list, partition_sizes, p_index, c_key, s_key, p_vec, changes, shift):

        key = (self.polynomial_representation, term_list, partition_sizes, p_index, c_key)
        poly = self._poly_cache.get(key, None)
        if poly is None:
            s = len(term_list)
            phase = (-1)**(sum(1 for t in term_list if t == 'p')//2)
            poly_coeffs = [
                self._evaluate_poly_coeffs(term_list, inds, delta, 0, self.keep_ints)
                for inds, delta in zip(p_vec, changes)
                if len(inds) > 0
            ]
            if any(nput.is_zero(c) for c in poly_coeffs):
                poly = 0
            else:
                # poly_coeffs = [
                #     c for c in poly_coeffs
                #     if c[0] != 1 or len(c) > 1
                # ]
                poly = make_product_polynomial(
                    poly_coeffs,
                    prefactor=phase,
                    idx=key,
                    steps=s if self.keep_ints else 0,
                    representation=self.polynomial_representation
                )
            self._poly_cache[key] = poly
        return poly

class HamiltonianExpansionTerm(OperatorExpansionTerm):
    def __init__(self, terms, order=None, identities=None, symmetrizers=None,
                 change_rules=None,
                 combine_kinetic_coefficients=False,
                 **opts):

        if identities is None:
            identities = np.arange(5) if identities is None else identities
        if symmetrizers is None:
            symmetrizers = PTTensorCoeffProductSum.symmetrizers()
        else:
            PTTensorCoeffProductSum.default_symmetrizers = symmetrizers

        if change_rules is None and nput.is_zero(order):
            change_rules = [()]

        super().__init__(terms, order=order, identities=identities, symmetrizers=symmetrizers,
                         change_rules=change_rules,
                         **opts)
        if not isinstance(combine_kinetic_coefficients,(bool,np.bool_)):
            raise TypeError('combine_kinetic_coefficients must be boolean')
        self.combine_kinetic_coefficients = bool(combine_kinetic_coefficients)

    def __fingerprint__(self):
        return super().__fingerprint__(), self.combine_kinetic_coefficients

    def get_core_poly(self, changes, shift=None):
        result = super().get_core_poly(changes,shift=shift)
        # Raw coefficient filters have their original meaning. Do not mix a
        # selected raw V/G subset into a channel which depends on both sources.
        if (not self.combine_kinetic_coefficients or not self.order
                or self.allowed_coefficients is not None or self.disallowed_coefficients is not None):
            return result
        rank = self.order+2
        if not any(identity == 1 and term == ['p']+['x']*self.order+['p']
                   for identity,term in zip(self.identities,self.terms)):
            return result
        from fractions import Fraction
        signs, inverse = _HamiltonianCoefficient.basis(rank)
        pairs = tuple(itertools.combinations(range(rank),2))
        collected, unchanged = {}, {}
        for key,poly in result.poly_obj.terms.items():
            coefficient, = key
            order,kind,*indices = coefficient
            if kind not in (0,1) or len(indices) != rank:
                unchanged[key] = poly
                continue
            ordered = tuple(sorted(indices))
            position = 0
            if kind == 1:
                position = 1+next(i for i,(a,b) in enumerate(pairs)
                    if sorted((ordered[a],ordered[b])) == sorted((indices[0],indices[-1])))
            eager = PolyPath.from_polynomial(poly).to_eager()
            leaves = eager.polys if isinstance(eager,ProductPTPolynomialSum) else [eager]
            for channel,weight in zip(signs,inverse[position]):
                if not weight: continue
                plus = tuple(i for i,s in zip(ordered,channel) if s == 1)
                minus = tuple(i for i,s in zip(ordered,channel) if s == -1)
                if len(plus)<len(minus): plus,minus = minus,plus
                recipe = _HamiltonianCoefficient(rank,len(plus))
                newkey = ((order,recipe)+recipe.canonical_indices(plus+minus),)
                monomials = collected.setdefault(newkey,collections.defaultdict(Fraction))
                for leaf in leaves:
                    if nput.is_zero(leaf): continue
                    if leaf.steps != rank: raise ValueError('unexpected primitive oscillator normalization')
                    prefactor = weight*Fraction(float(leaf.prefactor))
                    if isinstance(eager,ProductPTPolynomialSum): prefactor *= Fraction(float(eager.prefactor))
                    for powers in itertools.product(*(range(len(axis)) for axis in leaf.coeffs)):
                        value = prefactor
                        for axis,power in zip(leaf.coeffs,powers): value *= Fraction(float(axis[power]))
                        monomials[powers] += value
        terms = dict(unchanged)
        for key,monomials in collected.items():
            polynomial = 0
            for powers,value in monomials.items():
                if not value: continue
                axes = [np.eye(1,p+1,p,dtype=float)[0] for p in powers]
                polynomial += make_product_polynomial(axes,prefactor=float(value),steps=rank,
                    representation=self.polynomial_representation)
            if not nput.is_zero(polynomial): terms[key] = polynomial
        return result.mutate(PTTensorCoeffProductSum(terms,canonicalize=False))

    def __repr__(self):
        return "H[{}]".format(self.order)

class PerturbationOperator(PerturbationTheoryTerm):

    _energy_baseline = None
    def __init__(self, subterm):
        super().__init__(
            logger=subterm.logger,
            polynomial_representation=subterm.polynomial_representation
        )

        self.subterm = subterm

    _cache = {}

    @classmethod
    def lookup(cls, subterm):
        prod = cls._cache.get(subterm, None)
        if prod is None:
            prod = cls(subterm)
            cls._cache[subterm] = prod
        return prod

    def __repr__(self):
        return "{}|".format(self.subterm)

    def get_changes(self) -> 'dict[tuple[int], Any]':
        # The reduced resolvent acts in the complement of the reference-state
        # space.  In particular, an intermediate-normalized wavefunction
        # correction cannot contain the canonical zero-change component; its
        # diagonal normalization is supplied separately by
        # WavefunctionOverlapCorrection.
        return {k: None for k in self.subterm.changes if k != ()}
    def get_core_poly(self, changes, shift=None):
        base_term = self.subterm.get_core_poly(changes)
        if isinstance(base_term, PolynomialInterface):
            energy_shift = changes
            if shift is not None:
                elen = len(changes)
                blen = len(shift)
                energy_shift = tuple(
                    c + b
                    for c, b in
                    zip(
                        changes
                        if elen >= blen else
                        list(changes) + [0] * (blen - elen),
                        self._energy_baseline
                        if blen >= elen else
                        list(self._energy_baseline) + [0] * (elen - blen)
                    )
                )
            # if self._energy_baseline is not None:
            #     elen = len(changes)
            #     blen = len(self._energy_baseline)
            #     energy_shift = tuple(
            #         c - b
            #         for c, b in
            #         zip(
            #             changes
            #                 if elen >= blen else
            #             list(changes) + [0] * (blen - elen),
            #             self._energy_baseline
            #                 if blen >= elen else
            #             list(self._energy_baseline) + [0] * (elen - blen)
            #         )
            #     )
            prefactor = PTEnergyChangeProductSum.monomial((0,) + energy_shift, 1)
            base_term = base_term.mul_simple(prefactor)
        return base_term
    def get_poly_terms(self, changes, shift=None, **opts) -> 'SqrtChangePoly':
        # Keep the projector check here as well as in get_changes(): callers
        # may explicitly request a change, and a shifted nonzero change may
        # land back on the reference state.  Pad rather than zip-truncate so a
        # shift with a different dimension cannot hide a nonzero component.
        final_size = max(len(changes), 0 if shift is None else len(shift))
        final_change = np.pad(changes, (0, final_size - len(changes)))
        if shift is not None:
            final_change = final_change + np.pad(
                shift, (0, final_size - len(shift))
            )
        if not np.any(final_change):
            return 0
        return super().get_poly_terms(changes, shift=shift, **opts)

class _ShiftedEnergyBaseline(PerturbationTheoryTerm):
    def __init__(self, base_term:'PerturbationTheoryTerm'):
        super().__init__(
            logger=base_term.logger,
            polynomial_representation=base_term.polynomial_representation
        )
        self.base = base_term

    def __repr__(self):
        return "{{{}}}".format(self.base)

    def get_changes(self):
        changes = {}
        for k,v in self.base.changes.items():
            changes[k] = None
        return changes

    def get_core_poly(self, changes, shift=None) -> 'SqrtChangePoly':
        base_term = self.base.get_poly_terms(changes, shift=shift)
        if shift is None:
            final_state = changes
        else:
            if len(shift) < len(changes):
                shift = list(shift) + [0] * (len(changes) - len(shift))
            final_state = [c + s for c,s in zip(changes, shift)]
        # eshift = [-2*c for c in zip(changes, [0] * )
        return base_term.shift_energies([-c for c in final_state])

class ShiftedEnergyBaseline(PerturbationTheoryTerm):
    """
    Represents a term that will be multipled by on the left rather than the right
    for evaluating things like Y[1]M[0]Y[1], essentially changing raising operations to lowering
    """
    def __init__(self, base_term):
        super().__init__(
            logger=base_term.logger,
            polynomial_representation=base_term.polynomial_representation
        )
        self.base = base_term

    def __repr__(self):
        return "{{{}}}".format(self.base)
        # return repr(self.base)[::-1]

    # TODO: make sure changes supported are correct...
    #       NOTE: for now by doing nothing I'm assuming all operators provide symmetric changes...

    _cache = {}
    @classmethod
    def lookup(cls, base_term):
        prod = cls._cache.get(base_term, None)
        if prod is None:
            prod = cls(base_term)
            cls._cache[base_term] = prod
        return prod

    def get_changes(self):
        changes = {}
        for k,v in self.base.changes.items():
            # k = np.array([-i for i in k])
            # k = tuple(k)
            changes[k] = None
        return changes

    def get_poly_terms(self, changes, shift=None, **opts) -> 'SqrtChangePoly':
        # with self.debug_logging():
        # Operators are all Hermitian (overall) so we can just adjust the
        # energy shifts
        if shift is None:
            shift = [0] * len(changes)
        flip_changes = [-c for c in changes]

        shift_change = changes
        if len(shift_change) < len(shift):
            shift_change = list(shift_change) + [0] * (len(shift) - len(shift_change))
        flip_shift = [s+c for s,c in zip(shift, shift_change)]

        reexpressed_poly = self.base.get_poly_terms(flip_changes, flip_shift, **opts)#.shift(shift_change)
        if isinstance(reexpressed_poly, SqrtChangePoly):
            reexpressed_poly = reexpressed_poly.poly_obj.shift([c for c in shift_change])
        if nput.is_numeric(reexpressed_poly): return reexpressed_poly
        reexpressed_poly = reexpressed_poly.flip_energy_terms()
        return SqrtChangePoly(reexpressed_poly, changes, shift)

class ShiftedHamiltonianCorrection(PerturbationTheoryTerm):
    """
       Adds the wave function correction and the overlap term
       """

    def __init__(self, parent, order, allowed_terms=None, **opts):
        super().__init__(allowed_terms=allowed_terms, **opts)

        self.parent = parent
        self.order = order

    def get_serializer_key(self):  # to be overridden
        return self.__repr__()

    def __repr__(self):
        return "A[{}]".format(self.order)

    def get_changes(self):
        base_changes = {}
        for expr in self.expressions:
            for subchange in expr.changes:
                base_changes[subchange] = base_changes.get(subchange, [])
                base_changes[subchange].append(expr)  # just track which exprs generate the change
        base_changes[()] = None  # also the constant term
        return base_changes

    def get_subexpressions(self) -> 'Iterable[PerturbationTheoryTerm]':
        H = self.parent.hamiltonian_expansion
        E = self.parent.energy_correction
        k = self.order

        if k%2 == 1:
            if k >= len(H):
                return []
            else:
                return [H[k]]
        else:
            if len(H) > k:
                return [H[k], -E(k)]
            else:
                return [-E(k)]

class BlockHamiltonianCorrection(PerturbationTheoryTerm):
    """Shared F[k](delta): a change-indexed effective-Hamiltonian correction.

    Selection of block changes belongs to the evaluation plan. The equations
    expand through the same native products as W, O and Y; no evaluated block
    coefficient or reference occupation is stored in this term.
    """
    def __init__(self, parent, order, **options):
        super().__init__(**options)
        self.parent, self.order = parent, order

    def __repr__(self):
        return 'F[{}]'.format(self.order)

    get_serializer_key = __repr__

    def get_changes(self):
        return {change: None for term in self.expressions for change in term.changes}

    def get_subexpressions(self):
        H, k = self.parent.hamiltonian_expansion, self.order
        if k == 0:
            return [H[0]]
        Y = self.parent.full_wavefunction_correction
        N = self.parent.overlap_correction
        F = self.parent.block_hamiltonian_correction
        parts = [H[r] if r == k else H[r] * Y(k-r)
                 for r in range(1, min(k, len(H)-1)+1)]
        if not self.intermediate_normalization:
            parts += [-N(k-r) * F(r) for r in range(1, k)]
        base = operator_sum(*parts)
        if self.intermediate_normalization or k == 1:
            return [base]
        # Hermitian F = [H0,N] + B, with Hermitian N: F=HermitianPart(B).
        return [.5*base, .5*ShiftedEnergyBaseline.lookup(base)]

    def get_core_poly(self, changes, shift=None):
        if self.parent.resonance_independent and changes:
            return 0
        return _change_domain_polynomial(super().get_core_poly(changes), changes, True,
                                         unrestricted=self.parent.resonance_independent)


def _change_domain_polynomial(poly, change, retain, *, unrestricted=False):
    # In the unrestricted scalar reference every state is its own block:
    # P retains only zero change, Q all nonzero changes. This policy is
    # applied during construction, before any resonance domain is attached.
    if unrestricted:
        return poly if bool(change) != retain else 0
    if nput.is_numeric(poly):
        return poly
    if not change:
        return poly if retain else 0
    dag = PTTensorCoeffProductDAG.from_sum(poly.poly_obj)
    return poly.mutate(poly_obj=dag.with_change_domain(change, retain))


class _ChangeDomainIndex(int):
    """Compiler transport for an integer-change summation domain.

    These rows are generated only when lowering a domain-annotated DAG for
    evaluation. Equation export separates them into summation-domain data;
    they are never polynomial coefficients or frequency denominators.
    """
    _state_projection_factor = True
    _filtered_change_domain = True
    _cache = {}

    def __new__(cls, retain, endpoint):
        key = bool(retain), int(endpoint)
        if key not in cls._cache:
            obj = int.__new__(cls, 4096 + 2*int(retain) + endpoint)
            obj.retain, obj.endpoint = key
            cls._cache[key] = obj
        return cls._cache[key]

    def __reduce__(self):
        return type(self), (self.retain, self.endpoint)


@functools.lru_cache(maxsize=16)
def _filtered_block_labels(groups):
    return {tuple(state): index for index, group in enumerate(groups) for state in group}


def _filtered_change_factors(key, state, permutations, frequencies, baseline=None, bindings=()):
    state = np.asarray(state)[0]
    permutations = np.asarray(permutations, dtype=int)
    good = np.ones(len(permutations), dtype=bool)
    denominator = np.ones(len(permutations))
    labels = _filtered_block_labels(bindings[0]) if bindings else {}
    width = permutations.shape[-1]
    offset = np.zeros(width, dtype=int)
    if baseline is not None:
        offset[:min(width, len(baseline))] = baseline[:width]
    index = 0
    while index < len(key):
        row = key[index]
        role = row[0]
        if not isinstance(role, _ChangeDomainIndex):
            if getattr(role, '_state_projection_factor', False):
                raise TypeError('filtered change equations cannot contain foreign projector factors')
            denominator *= frequencies @ np.asarray(row[1:])
            index += 1
            continue
        following = key[index+1]
        if role.endpoint != 0 or following[0] != _ChangeDomainIndex(role.retain, 1):
            raise ValueError('change-domain index pair was not preserved')
        start = np.asarray(row[1:], dtype=int) + offset
        end = np.asarray(following[1:], dtype=int) + offset
        valid = (np.all(state[permutations]+start >= 0, axis=-1)
                 & np.all(state[permutations]+end >= 0, axis=-1))
        if np.array_equal(start, end):
            inside = np.ones(len(permutations), dtype=bool)
        else:
            if not bindings:
                raise ValueError('filtered block contraction requires resonance-change bindings')
            origins = np.broadcast_to(state, (len(permutations), len(state))).copy()
            targets = origins.copy()
            rows = np.arange(len(permutations))[:, None]
            origins[rows, permutations] += start
            targets[rows, permutations] += end
            inside = np.fromiter((
                source in labels and labels[source] == labels.get(tuple(target))
                for source, target in zip(map(tuple, origins), targets)
            ), dtype=bool, count=len(permutations))
        good &= valid & (inside if role.retain else ~inside)
        index += 2
    return good, denominator


def _split_change_equations(key, polynomial):
    """Separate ordinary polynomial equations from integer summation domains."""
    energy = PTEnergyChangeProductSum.from_polynomial(polynomial)
    for path, numerator in energy.terms.items():
        domains, divisors = [], []
        index = 0
        while index < len(path):
            row = path[index]
            if isinstance(row[0], _ChangeDomainIndex):
                end = path[index+1]
                start = tuple(row[1:])
                delta = tuple(b-a for a, b in zip(start, end[1:]))
                domains.append((row[0].retain, start, delta))
                index += 2
            else:
                if getattr(row[0], '_state_projection_factor', False):
                    raise TypeError('foreign selection factor in filtered equation')
                divisors.append(row)
                index += 1
        yield tuple(domains), key, PTEnergyChangeProductSum({tuple(divisors): numerator})


class WavefunctionCorrection(PerturbationTheoryTerm):
    def __init__(self, parent, order, allowed_terms=None, **opts):
        super().__init__(allowed_terms=allowed_terms, **opts)

        self.parent = parent
        self.order = order


    def get_serializer_key(self): # to be overridden
        return self.__repr__()
    def __repr__(self):
        return "W[{}]".format(self.order)

    def get_changes(self):
        base_changes = {}
        for expr in self.expressions:
            for subchange in expr.changes:
                base_changes[subchange] = base_changes.get(subchange, [])
                base_changes[subchange].append(expr) # just track which exprs generate the change
        # base_changes[()] = None # also the constant term
        return base_changes

    def get_core_poly(self, changes, shift=None):
        base = super().get_core_poly(changes)
        return (_change_domain_polynomial(base, changes, False,
                                           unrestricted=self.parent.resonance_independent)
                if self.parent.filtered_block_recurrence else base)

    def get_subexpressions(self):
        if self.parent.filtered_block_recurrence or self.parent.resonance_independent:
            k = self.order
            H = self.parent.hamiltonian_expansion
            Y = self.parent.full_wavefunction_correction
            X = self.parent.wavefunction_correction
            # The unrestricted scalar reference has the same normalized
            # recurrence with E in place of F. Energy feedback belongs inside
            # the resolvent and acts on the source of each correction.
            F = (self.parent.block_hamiltonian_correction if self.parent.filtered_block_recurrence
                 else self.parent.energy_correction)
            parts = [H[r] if r == k else H[r]*Y(k-r)
                     for r in range(1, min(k, len(H)-1)+1)]
            parts += [-X(k-r)*F(r) for r in range(1, k)]
            return [-PerturbationOperator.lookup(operator_sum(*parts))] if parts else []
        H = self.parent.shifted_hamiltonian_correction
        E = self.parent.energy_correction
        W = self.parent.wavefunction_correction
        k = self.order

        base_terms = [
            E(k-i) * W(i)
            for i in range(1, k)
        ] + [
            -PerturbationOperator.lookup(H(k-i) * W(i))
                if i > 0 else
            -PerturbationOperator.lookup(H(k))
            for i in range(0, k)
        ]

        return base_terms

class EnergyCorrection(PerturbationTheoryTerm):
    def __init__(self, parent, order, allowed_terms=None, **opts):
        super().__init__(allowed_terms=allowed_terms, **opts)
        self.parent = parent
        self.order = order

    def get_serializer_key(self):  # to be overridden
        return self.__repr__()
    def __repr__(self):
        return "E[{}]".format(self.order)

    def get_changes(self) -> 'dict[tuple[int], Any]':
        return {():None}

    def get_subexpressions(self) -> 'Iterable[PerturbationTheoryTerm]':
        if self.parent.filtered_block_recurrence:
            return [] if self.order % 2 else [self.parent.block_hamiltonian_correction(self.order)]
        H = self.parent.hamiltonian_expansion
        J = self.parent.shifted_hamiltonian_correction
        W = self.parent.full_wavefunction_correction
        k = self.order

        if k % 2 == 1:
            return []

        return ([H[k]] if len(H) > k else []) + [
            (J(k - i) * W(i))
            for i in range(1, k)
        ]

    def get_poly_terms(self, changes, shift=None, **opts):
        if len(changes) != 0:
            return 0
        else:
            return super().get_poly_terms(changes, shift=shift, **opts)

    # def get_poly_terms(self):
    #     H = self.parent.hamiltonian_expansion
    #     E = self.parent.energy_correction
    #     W = self.parent.wavefunction_correction
    #     k = self.order
    #
    #     # Recursive PT equations
    #     return H[k]([]) + sum(
    #         (H[k-i]*W(i))([]) - E(k-i)([])*W(i)([]) #TODO: implement caching for this kind of term...
    #         for i in range(1, k)
    #     ) # 1 to k-1

    def __call__(self, changes, shift=None, coeffs=None, freqs=None, check_sorting=None, simplify=True):
        if len(changes) > 0:
            raise ValueError("no")
        return super().__call__((), shift=shift, coeffs=coeffs, freqs=freqs, check_sorting=False, simplify=simplify)

class WavefunctionOverlapCorrection(PerturbationTheoryTerm):

    def __init__(self, parent, order, allowed_terms=None, degenerate_changes=None, **opts):
        super().__init__(allowed_terms=allowed_terms, **opts)

        self.parent = parent
        self.order = order
        self.change = degenerate_changes

    def get_serializer_key(self): # to be overridden
        if self.change is None:
            return "O[{}]".format(self.order)
        return "O[{},{}]".format(self.change, self.order)
    def __repr__(self):
        if self.change is None:
            return "O[{}]".format(self.order)
        else:
            return "O[{},{}]".format(
                format_quantum_change_list(self.change), self.order
            )

    def get_changes(self) -> 'dict[tuple[int], Any]':
        if self.parent.filtered_block_recurrence:
            changes = {change: None for term in self.expressions for change in term.changes}
            if self.change is not None:
                return {change: None for change in changes if change}
            return changes
        changes = self.change
        if changes is None:
            changes = [()]
        return {ch:None for ch in changes}

    def get_subexpressions(self) -> 'Iterable[PerturbationTheoryTerm]':
        if self.parent.filtered_block_recurrence:
            if self.order == 0:
                return [IdentityTerm(logger=self.logger, polynomial_representation=self.polynomial_representation)]
            if self.intermediate_normalization:
                return []
        W = self.parent.full_wavefunction_correction
        k = self.order
        L = lambda o:ShiftedEnergyBaseline(W(o))

        # if k == 1:
        #     return []

        return [
            -1/2 * (L(k - i) * W(i))
            for i in range(1, k)
        ]

    def get_core_poly(self, changes, shift=None):
        if self.parent.resonance_independent and changes:
            return 0
        base = super().get_core_poly(changes)
        return (_change_domain_polynomial(base, changes, True,
                                          unrestricted=self.parent.resonance_independent)
                if self.parent.filtered_block_recurrence else base)


class FullWavefunctionCorrection(PerturbationTheoryTerm):
    """
    Adds the wave function correction and the overlap term
    """

    def __init__(self, parent, order, allowed_terms=None, **opts):
        super().__init__(allowed_terms=allowed_terms, **opts)

        self.parent = parent
        self.order = order

    def get_serializer_key(self): # to be overridden
        return self.__repr__()
    def __repr__(self):
        return "Y[{}]".format(self.order)

    def get_changes(self):
        base_changes = {}
        for expr in self.expressions:
            for subchange in expr.changes:
                base_changes[subchange] = base_changes.get(subchange, [])
                base_changes[subchange].append(expr)  # just track which exprs generate the change
        base_changes[()] = None  # also the constant term
        return base_changes

    def get_subexpressions(self) -> 'Iterable[PerturbationTheoryTerm]':
        W = self.parent.wavefunction_correction
        O = self.parent.overlap_correction
        k = self.order

        if self.parent.filtered_block_recurrence:
            if k == 0:
                return [IdentityTerm(logger=self.logger, polynomial_representation=self.polynomial_representation)]
            return [W(k)] if self.intermediate_normalization else [W(k), O(k)]
        if (k % 2) == 1 or self.intermediate_normalization:
            return [W(k)]
        else:
            return [W(k), O(k)]

class OperatorCorrection(PerturbationTheoryTerm):

    def __init__(self, parent, order, operator_type=None, allowed_terms=None, wavefunction_generator=None, base_index=None, **opts):
        super().__init__(allowed_terms=allowed_terms, **opts)
        self.parent = parent
        self.order = order
        self.type = self.get_type_key(operator_type)
        self.expansion = parent.operator_expansion_terms(
            order,
            logger=self.logger,
            base_index=base_index,
            operator_type=operator_type,
            polynomial_representation=self.polynomial_representation
        )
        self._wavefunction_generator = wavefunction_generator

    @classmethod
    def get_type_key(cls, operator_type):
        if operator_type is None:
            return 1
        elif isinstance(operator_type, str):
            if operator_type == 'transition_moment':
                operator_type = 'TM'
        return operator_type

    def get_serializer_key(self):  # to be overridden
        return self.__repr__()

    repr_key = "M"
    def get_repr_key(self):
        if self.type == "TM":
            return self.repr_key
        else:
            return '{}@{}'.format(self.repr_key, self.type)
    def __repr__(self):
        return "<n|{}|m>({})".format(self.get_repr_key(), self.order)

    def get_changes(self):
        base_changes = {}
        for expr in self.expressions:
            for subchange in expr.changes:
                base_changes[subchange] = base_changes.get(subchange, [])
                base_changes[subchange].append(expr) # just track which exprs generate the change
        base_changes[()] = None # also the constant term
        return base_changes
    def default_wavefunction_generator(self, o:int):
        if self.parent.filtered_block_recurrence:
            return self.parent.block_main_wavefunction_correction(o)
        return self.parent.full_wavefunction_correction(o)
    def wavefunction_generator(self, o:int):# -> '(int)->FullWavefunctionCorrection':
        if self._wavefunction_generator is None:
            return self.default_wavefunction_generator(o)
        else:
            return self._wavefunction_generator(o)

    def get_subexpressions(self,
                          bra_wavefunction_generator=None,
                          ket_wavefunction_generator=None,
                          bounds=None,
                          const_zeros=None
                          ) -> 'Iterable[PerturbationTheoryTerm]':
        if const_zeros is None:
            const_zeros = [False, False]
        if ket_wavefunction_generator is None:
            const_zeros[1] = True
            ket_wavefunction_generator = self.wavefunction_generator
        W = ket_wavefunction_generator
        M = self.expansion
        k = self.order
        if bra_wavefunction_generator is None:
            const_zeros[0] = True
            bra_wavefunction_generator = lambda o:ShiftedEnergyBaseline(W(o))
        L = bra_wavefunction_generator

        if bounds is None:
            bounds = [0, k+1]
        n_min, n_max = bounds
        terms = [
            (i, k - i - j, j)
            for i in range(n_min, n_max)
            for j in range(n_min, n_max - i)
        ]
        if self.allowed_terms is not None:
            allowed_terms = {tuple(t) for t in self.allowed_terms}
            terms = [t for t in terms if t in allowed_terms]

        ub, uk = const_zeros
        exprs = [
            M[a]
                if (ub and uk) and i == 0 and j == 0 else
            L(j) * M[a] * W(i)
                if i > 0 and j > 0 else
            M[a] * W(i)
                if ub and i > 0 else
            L(j) * M[a]
                if uk and j > 0 else
            None
            for i, a, j in terms
        ]
        exprs = tuple(e for e in exprs if e is not None)

        return exprs

class OperatorDegenerateCorrection(OperatorCorrection):

    def __init__(self, parent, order, degenerate_changes=None, operator_type=None, allowed_terms=None, **opts):
        super().__init__(parent, order, operator_type=operator_type, allowed_terms=allowed_terms, **opts)
        self.degenerate_changes = degenerate_changes
        self.left = self.Left(parent, order, self) # hack to minimize code necessary
        self.right = self.Right(parent, order, self) # hack to minimize code necessary
        self.both = self.Both(parent, order, self) # hack to minimize code necessary

    def get_serializer_key(self):
        return "<n||{}||m>({},{})".format(
            self.get_repr_key(), self.degenerate_changes, self.order
        )

    def __repr__(self):
        return "<n||{}||m>({},{})".format(
            self.get_repr_key(),
            format_quantum_change_list(self.degenerate_changes), self.order
        )

    def default_wavefunction_generator(self, o:int):
        return self.parent.overlap_correction(o, degenerate_changes=self.degenerate_changes)

    class Left(OperatorCorrection):
        def __init__(self, parent:OperatorCorrection, order, real_parent):
            super().__init__(parent, order, logger=real_parent.logger)
            self.op = real_parent
        def get_serializer_key(self):
            return "<n||{}|m>({},{})".format(
                self.op.get_repr_key(), self.op.degenerate_changes, self.order
            )
        def __repr__(self):
            return "<n||{}|m>({},{})".format(
                self.op.get_repr_key(),
                format_quantum_change_list(self.op.degenerate_changes), self.order
            )
        def get_subexpressions(self):
            woof = self.op.get_left_degenerate_expressions()
            return woof
    class Right(OperatorCorrection):
        def __init__(self, parent, order, real_parent):
            super().__init__(parent, order, logger=real_parent.logger)
            self.op = real_parent
        def get_serializer_key(self):
            return "<n|{}||m>({},{})".format(
                self.op.get_repr_key(), self.op.degenerate_changes, self.order
            )
        def __repr__(self):
            return "<n|{}||m>({},{})".format(
                self.op.get_repr_key(),
                format_quantum_change_list(self.op.degenerate_changes), self.order
            )
        def get_subexpressions(self):
            return self.op.get_right_degenerate_expressions()
    class Both(OperatorCorrection):
        def __init__(self, parent, order, real_parent):
            super().__init__(parent, order, logger=real_parent.logger)
            self.op = real_parent
        def get_serializer_key(self):
            return "<n||{}||m>({},{})".format(
                self.op.get_repr_key(), self.op.degenerate_changes, self.order
            )
        def __repr__(self):
            return "<n||{}||m>({},{})".format(
                self.op.get_repr_key(),
                format_quantum_change_list(self.op.degenerate_changes), self.order
            )
        def get_subexpressions(self):
            return self.op.get_both_degenerate_expressions()
    def get_subexpressions(self,
                          bra_wavefunction_generator=None,
                          ket_wavefunction_generator=None
                          ) -> 'Iterable[PerturbationTheoryTerm]':
        raise NotImplementedError("ill-defined")
        W1 = lambda o: self.parent.full_wavefunction_correction(o)
        L1 = lambda o: ShiftedEnergyBaseline(W1(o))
        W2 = lambda o: self.parent.overlap_correction(o, degenerate_changes=self.degenerate_changes)
        L2 = lambda o: ShiftedEnergyBaseline(W2(o))
        return (
                super().get_subexpressions(
                    bra_wavefunction_generator=L1,
                    ket_wavefunction_generator=W2
                )
                + super().get_subexpressions(
                    bra_wavefunction_generator=L2,
                    ket_wavefunction_generator=W1
                )
                + super().get_subexpressions(
                    bra_wavefunction_generator=L2,
                    ket_wavefunction_generator=W2
                )
        )

    def get_right_degenerate_expressions(self,
                                        bra_wavefunction_generator=None,
                                        ket_wavefunction_generator=None
                                        ) -> 'Iterable[PerturbationTheoryTerm]':
        if self.parent.filtered_block_recurrence:
            main = self.parent.block_main_wavefunction_correction
            cross = lambda o: self.parent.overlap_correction(o, degenerate_changes=())
            return OperatorCorrection.get_subexpressions(self,
                bra_wavefunction_generator=lambda o: ShiftedEnergyBaseline.lookup(cross(o)),
                ket_wavefunction_generator=main, const_zeros=[False, True], bounds=None)
        # W1 = lambda o: self.parent.full_wavefunction_correction(o)
        # L1 = lambda o: ShiftedEnergyBaseline(W1(o))
        W2 = lambda o: self.parent.overlap_correction(o, degenerate_changes=self.degenerate_changes)
        L2 = lambda o: ShiftedEnergyBaseline(W2(o))
        return super().get_subexpressions(
            bra_wavefunction_generator=L2,
            ket_wavefunction_generator=None,
            # const_zeros=[True, False]
        )

    def get_left_degenerate_expressions(self,
                                         bra_wavefunction_generator=None,
                                         ket_wavefunction_generator=None
                                         ) -> 'Iterable[PerturbationTheoryTerm]':
        if self.parent.filtered_block_recurrence:
            main = self.parent.block_main_wavefunction_correction
            cross = lambda o: self.parent.overlap_correction(o, degenerate_changes=())
            return OperatorCorrection.get_subexpressions(self,
                bra_wavefunction_generator=lambda o: ShiftedEnergyBaseline.lookup(main(o)),
                ket_wavefunction_generator=cross, const_zeros=[True, False], bounds=None)
        # W1 = lambda o: self.parent.full_wavefunction_correction(o)
        # L1 = lambda o: ShiftedEnergyBaseline(W1(o))
        W2 = lambda o: self.parent.overlap_correction(o, degenerate_changes=self.degenerate_changes)
        L2 = lambda o: ShiftedEnergyBaseline(W2(o))
        return super().get_subexpressions(
            bra_wavefunction_generator=None,
            ket_wavefunction_generator=W2,
            # const_zeros=[False, True]
        )

    def get_both_degenerate_expressions(self,
                                        bra_wavefunction_generator=None,
                                        ket_wavefunction_generator=None
                                        ) -> 'Iterable[PerturbationTheoryTerm]':
        if self.parent.filtered_block_recurrence:
            main = self.parent.block_main_wavefunction_correction
            cross = lambda o: self.parent.overlap_correction(o, degenerate_changes=())
            return OperatorCorrection.get_subexpressions(self,
                bra_wavefunction_generator=lambda o: ShiftedEnergyBaseline.lookup(cross(o)),
                ket_wavefunction_generator=cross, const_zeros=[False, False], bounds=[1, self.order])
        W2 = lambda o: self.parent.overlap_correction(o, degenerate_changes=self.degenerate_changes)
        L2 = lambda o: ShiftedEnergyBaseline(W2(o))
        return super().get_subexpressions(
            bra_wavefunction_generator=L2,
            ket_wavefunction_generator=W2,
            bounds=[1, self.order],
            # const_zeros=[False, False]
        )

class DiagonalHamiltonian(OperatorExpansionTerm):
    def __init__(self):
        super().__init__(
            [['x', 'x']],
            order=0,
            identities=[5]
        )
    def __repr__(self):
        return HamiltonianExpansionTerm.__repr__(self)
    def get_changes(self):
        return {():None}
class ReexpressedHamiltonian(OperatorCorrection):
    repr_key = "H"

    def __init__(self, parent, order, allowed_terms=None, degenerate_changes=None, hamiltonian_corrections=None, **opts):
        super().__init__(parent, order, allowed_terms=allowed_terms, **opts)
        self.expansion = self.prep_expansion(
            parent.hamiltonian_expansion,
            hamiltonian_corrections,
            logger=self.logger,
            allowed_coefficients=self.allowed_coefficients,
            disallowed_coefficients=self.disallowed_coefficients
        )

    @classmethod
    def prep_expansion(self, base_expansion, corrs, *,
                       logger,
                       allowed_coefficients,
                       disallowed_coefficients,
                       base_index=5):
        if corrs is None:
            return base_expansion
        if len(corrs) < len(base_expansion):
            corrs = list(corrs) + [[]] * (len(base_expansion) - len(corrs))
        # symmetrizer = list(symmetrizers)
        new_exp = []
        for o, (b, term_list) in enumerate(zip(base_expansion, corrs)):
            if len(term_list) > 0:
                symmetrizers = list(b.symmetrizers)
                new_idx = base_index + len(term_list)
                ids = list(range(base_index, new_idx))
                if len(symmetrizers) < new_idx:
                    symmetrizers = symmetrizers + [None] * (new_idx - len(symmetrizers))
                    for i, t in zip(ids, term_list):
                        if 'p' in t:
                            symmetrizers[i] = PTTensorCoeffProductSum._kinetic_symmetrizer
                        else:
                            symmetrizers[i] = PTTensorCoeffProductSum._potential_symmetrizer

                # base_index = new_idx
                new_exp.append(
                    PerturbationTheoryTermSum(
                        b,
                        HamiltonianExpansionTerm(
                            term_list,
                            order=o,
                            identities=ids,
                            symmetrizers=tuple(symmetrizers),
                            logger=logger,
                            allowed_coefficients=allowed_coefficients,
                            disallowed_coefficients=disallowed_coefficients
                        )
                    )
                )
            else:
                new_exp.append(b)
        return new_exp

    def get_repr_key(self):
        return self.repr_key
class ReexpressedHamiltonianDegenerateCorrection(OperatorDegenerateCorrection):
    repr_key = "H"
    def __init__(self, parent, order, allowed_terms=None, hamiltonian_corrections=None, **opts):
        super().__init__(parent, order, allowed_terms=allowed_terms, **opts)
        self.expansion = ReexpressedHamiltonian.prep_expansion(
            parent.hamiltonian_expansion,
            hamiltonian_corrections,
            logger=self.logger,
            allowed_coefficients=self.allowed_coefficients,
            disallowed_coefficients=self.disallowed_coefficients
        )

    def get_repr_key(self):
        return self.repr_key

class ScaledPerturbationTheoryTerm(PerturbationTheoryTerm):
    #TODO: refactor since inheritance isn't really the right paradigm here
    def __init__(self, base_term:'PerturbationTheoryTerm', scaling):
        super().__init__(
            logger=base_term.logger,
            polynomial_representation=base_term.polynomial_representation
        )
        self.prefactor = scaling
        self.base = base_term

    _cache = {}
    @classmethod
    def lookup(cls, base_term, scaling):
        key = (base_term, scaling)
        prod = cls._cache.get(key, None)
        if prod is None:
            prod = cls(base_term, scaling)
            cls._cache[key] = prod
        return prod

    def __repr__(self):
        if nput.is_numeric(self.prefactor) and self.prefactor == -1:
            return "-{}".format(self.base)
        else:
            return "{}*{}".format(self.prefactor, self.base)

    def get_changes(self) -> 'dict[tuple[int], Any]':
        return {k:None for k in self.base.changes}
    def get_core_poly(self, changes, shift=None) -> 'SqrtChangePoly':
        # with self.debug_logging():
        base = self.base.get_poly_terms(changes, shift=shift)
        if nput.is_numeric(base): return base
        if isinstance(base, SqrtChangePoly):
            new = type(base)(self.prefactor * base.poly_obj, base.poly_change, base.shift_start)
        else:
            new = SqrtChangePoly(self.prefactor * base, changes, shift)
        return new

class PerturbationTheoryTermSum(PerturbationTheoryTerm):

    _sum_cache = weakref.WeakValueDictionary()

    @classmethod
    def lookup(cls, *terms):
        expression = cls._sum_cache.get(terms)
        if expression is None:
            expression = cls(*terms)
            cls._sum_cache[terms] = expression
        return expression

    def __init__(self, *terms):
        super().__init__(
            logger=terms[0].logger,
            polynomial_representation=terms[0].polynomial_representation
        )
        self.terms = terms
    def __repr__(self):
        return "+".join("{}".format(t) for t in self.terms)
    def get_changes(self):
        base_changes = {}
        for expr in self.expressions:
            for subchange in expr.changes:
                base_changes[subchange] = base_changes.get(subchange, [])
                base_changes[subchange].append(expr) # just track which exprs generate the change
        base_changes[()] = None # also the constant term
        return base_changes

    def get_subexpressions(self):
        return self.terms

class PerturbationTheoryTermProduct(PerturbationTheoryTerm):

    _cache = {}
    _poly_product_cache = {}
    def __init__(self, post_op, pre_op):
        super().__init__(
            logger=post_op.logger,
            polynomial_representation=post_op.polynomial_representation
        )

        self.gen1 = pre_op # we apply right-to-left
        self.gen2 = post_op

    @classmethod
    def lookup(cls, post_op, pre_op):
        key = (post_op, pre_op)
        prod = cls._cache.get(key, None)
        if prod is None:
            prod = cls(post_op, pre_op)
            cls._cache[key] = prod
        return prod

    def __repr__(self):
        return "{}*{}".format(self.gen1, self.gen2)

    _change_map = {}
    _perms_map = {}

    @classmethod
    def get_change_classes(cls, changes):
        classes = []
        for k in changes:
            if len(k) == 2 and isinstance(k[0], tuple): continue

            n = len(k)
            if n + 1 > len(classes): classes.extend([[] for _ in range(n+1 - len(classes))])
            classes[n].append(k)
        return [np.array(ch) for ch in classes]

    # change_tup = collections.namedtuple([])

    _combination_inds = {}
    _combination_comps = {}
    @classmethod
    def get_combination_inds(cls, n, r):
        if r == 0:
            return np.array([[]], dtype=int)
        elif r == n:
            return np.arange(n)[np.newaxis]
        elif (n, r) not in cls._combination_inds:
            cls._combination_inds[(n, r)] = nput.combination_indices(n, r)
        return cls._combination_inds[(n, r)]
    @classmethod
    def get_combination_comp(cls, n, r):
        if r == 0:
            return np.arange(n)[np.newaxis]
        elif r == n:
            return np.array([[]], dtype=int)
        elif (n, r) not in cls._combination_comps:
            cls._combination_comps[(n, r)] = nput.index_complement(n, cls.get_combination_inds(n, r))[0]
        return cls._combination_comps[(n, r)]

    @staticmethod
    def _fill_change_data(full_changes, n, new_changes,
                          init_changes, subsorts, contacts, rems):
        if n not in full_changes:
            full_changes[n] = {
                'changes': [],
                'sorts': [[], []],
                'init_changes': [[], []],
                'contacts': [[], []],
                'rems': [[], []]
            }

        change_data = full_changes[n]
        change_data['changes'].append(new_changes)
        for storage, contact in zip(change_data['init_changes'], init_changes):
            storage.append(contact)
        for storage, contact in zip(change_data['contacts'], contacts):
            storage.append(contact)
        for storage, contact in zip(change_data['rems'], rems):
            storage.append(contact)
        for storage, contact in zip(change_data['sorts'], [subsorts, np.argsort(subsorts, axis=1)]):
            storage.append(contact)


    # TODO: move this to a module-level cache that can be cleared
    _permutation_cache = MaxSizeCache(128, cache_type='fifo')
    @classmethod
    def _get_uperms(cls, p_list):
        key = tuple(p_list)
        if key not in cls._permutation_cache:
            cls._permutation_cache[key] = UniquePermutations(p_list).permutations(return_indices=True)
        return cls._permutation_cache[key]

    # _change_cache = MaxSizeCache(128, cache_type='fifo')
    # TODO: cache changes when the generators are the same...
    @classmethod
    def _fill_class_counts(cls, data, which, c):
        if data[which] is None:
            if len(c) == 0:
                data[which] = [np.array([], dtype=int), np.array([], dtype=int), np.array([], dtype=int)]
            else:
                classes, idx, counts = np.unique(c, return_counts=True, return_index=True)
                ord = cls.change_sort(classes)
                classes = classes[ord]
                counts = counts[ord]
                idx = idx[ord]
                mapping = np.concatenate([np.full((cnt,), i) for i, cnt in zip(idx, counts)])
                # perms = np.array([
                #     np.concatenate(p)
                #     for p in itertools.product(*[
                #         list(itertools.permutations(range(i, i+cnt)))
                #         for i, cnt in zip(idx, counts)
                #         ])
                # ])
                data[which] = [classes, counts, mapping]  # , perms]
        return data[which]

    @classmethod
    def _fill_inds_data(cls, data, which, counts, mapping, m, permute):
        if data[which][m - 1] is None:
            n1 = np.sum(counts)
            inds_1 = cls.get_combination_inds(n1, m)
            rems_1 = cls.get_combination_comp(n1, m)

            inds_mapped = mapping[inds_1]
            _, idx = np.unique(inds_mapped, axis=0, return_index=True)
            inds = inds_1[np.sort(idx),]
            rems = rems_1[np.sort(idx),]

            # inds = inds_1
            # rems = rems_1

            if permute:
                mperms = nput.permutation_indices(m, m)
                inds = nput.vector_take(inds, mperms).reshape(-1, m)
                if n1 > m:
                    rems = np.broadcast_to(rems[:, np.newaxis, :],
                                           (len(rems), len(mperms), n1 - m)).reshape(inds.shape[0], n1 - m)
                else:
                    rems = np.empty((inds.shape[0], 0), dtype=int)

            data[which][m - 1] = [inds, rems]
        return data[which][m - 1]

    @classmethod
    def _prep_product_change_data(cls,
                                  i1, ch1, inds_1, rind_1, n1,
                                  i2, ch2, inds_2, rind_2, n2,
                                  os,
                                  pad_change, initial_sort,
                                  nz, x1, x2,
                                  upairs):

        ch2_test = np.array([ch2[i2] for i2 in inds_2[x2]])
        ch2_sort = cls.change_sort(ch2_test)
        test_key = (
            tuple(ch1[inds_1[x1][s]] for s in ch2_sort),
            tuple(ch2_test[s] for s in ch2_sort)
        )
        if test_key in upairs: return None, None
        upairs.add(test_key)

        # we have to handle two combinatoric problems, now
        # firstly, we need to make sure that combinations like
        #    2  1 -1
        #    2  -  1  1  1
        # work properly, where we expect to have
        #    2  1  0  1  1 [ [0, 1, 3, 4, 2],  [0, 3, 1, 4, 2], [0, 3, 4, 1, 2] ]
        # i.e. we have [ [0] x perms([i, j, j]) x [2] ]
        # this needs to be handled because the [1, 1] block from the second poly is not equivalent
        # to the [1] from the first poly and we need to preserve symmetry
        # The second problem is in how we handle the [-1], [1] block, which could have been generated
        # by _any_ of the [1]s in the second poly, and so I think we need to scale this by
        # the number of equivalent choices?
        # [1, 1, 1] x [-1] -> [1, 1] + free index (might be handled by the later summation?)

        # we actually have a _third_ problem
        # consider [1, -1, -1] @ [-1, 1, 1] -> [0, 0, 0]
        # where we tneed to make sure that the blocks that map to the same value
        # are permuted over the different paths that get you there
        # i.e. if given two paths ([1, -1]) and ([-1, 1]) to 0 we need to make sure that we
        #      permute appropriately, meaning we take the unique permutations of [i, j, j]

        diff_pos = np.nonzero(np.diff(pad_change))
        if len(diff_pos) == 0 or len(diff_pos[0]) == 0:
            ind_blocks = [initial_sort]
            block_vals = [pad_change[0]]
        else:
            ind_blocks = np.split(initial_sort, diff_pos[0] + 1)
            block_vals = [x[0] for x in np.split(pad_change, diff_pos[0] + 1)]

        # We will handle the first problem by creating a proxy vector that we can take the unique
        # permutation of for each block in the final changes
        # To do so we take each of the ind_blocks and for the _remainder_ elements we assign
        # either 0 or 1 depending on which poly they came from, any change element by contrast is
        # treated as unique
        # And for the _product_ elements, we assign a value based on the left and right
        # indices
        inv_sort = np.argsort(initial_sort)
        perm_blocks = []
        product_inds = []
        proxy_idx = 0  # UniquePermutations expects largest to smallest...
        path_map = {}
        for left_i, right_j in zip(inds_1[x1], inds_2[x2]):
            key = (ch1[left_i], ch2[right_j])
            if key in path_map:
                val = path_map.get(key)
            else:
                val = proxy_idx
                proxy_idx += 1
                path_map[key] = val
            product_inds.append(val)
        proxy_left = len(product_inds)
        proxy_right = proxy_left + n1 - os
        origins = np.array(
            product_inds
            + [proxy_left] * (n1 - os)
            + [proxy_right] * (n2 - os)
        )
        for p, v in zip(ind_blocks, block_vals):
            # if v == 0:
            #     perm_blocks.append([p])
            # else:
            perms, _ = cls._get_uperms(origins[p])
            perm_blocks.append(nput.vector_take(p, perms))

        total_sorts = np.broadcast_to(
            initial_sort[np.newaxis],
            (np.prod([len(x) for x in perm_blocks], dtype=int), len(initial_sort))
        ).copy()
        for i, sort_prod in enumerate(itertools.product(*perm_blocks)):
            # total_perm = np.concatenate(list(sort_prod))[ch_inv]
            # total_sorts[i][rem_sort] = total_sorts[i][rem_sort][total_perm]
            total_sorts[i] = np.concatenate(sort_prod)
        # now we need to reduce the sorts so as to not dupe indices from the "origin"

        change = tuple(pad_change[:nz])

        return change, (
            ch1,
            ch2,
            inds_1[x1],
            inds_2[x2],
            rind_1[x1],
            rind_2[x2],
            total_sorts
        )

    def get_changes(self):
        #TODO: prevent _all_ changes from being enumerated, just handle the target
        #      cases by finding possible generaotrs?

        # change_keys = (
        #     # TODO: ensure these are always sorted to not have to pay the extra cost
        #     tuple(self.gen1.changes.keys()),
        #     tuple(self.gen2.changes.keys())
        # )
        #
        # if change_keys in self._change_cache:
        #     return self._change_cache[change_keys]
        # elif (change_keys[1], change_keys[0]) in self._change_cache:
        #     return self._change_cache[change_keys]

        new_changes = {}

        class_count_data_1 = [None] * len(self.gen1.changes)
        class_count_data_2 = [None] * len(self.gen2.changes)

        ind_rem_data_1 = [[None] * len(_) for _ in self.gen1.changes]
        ind_rem_data_2 = [[None] * len(_) for _ in self.gen2.changes]

        for i1,ch1 in enumerate(self.gen1.changes):
            if len(ch1) == 2 and isinstance(ch1[0], tuple): continue
            ch1 = np.array(ch1)
            if not np.equal(self.change_sort(ch1), np.arange(len(ch1))).all():
                raise ValueError(ch1)
            classes_1, counts_1, mapping_1 = self._fill_class_counts(class_count_data_1, i1, ch1)
            n1 = np.sum(counts_1)

            for i2, ch2 in enumerate(self.gen2.changes):
                if len(ch2) == 2 and isinstance(ch2[0], tuple): continue

                # we're back to one change at a time, since most of the manipulations scale only across the
                # same number of terms and the same number of unique classes
                ch2 = np.array(ch2)
                if not np.equal(self.change_sort(ch2), np.arange(len(ch2))).all():
                    raise ValueError(self, self.gen2, ch2)
                classes_2, counts_2, mapping_2 = self._fill_class_counts(class_count_data_2, i2, ch2)
                n2 = np.sum(counts_2)

                for os in range(min(n1, n2)+1): # number of contacts
                    if os == 0: # direct prod
                        inds_1 = np.array([[]], dtype=int)
                        inds_2 = np.array([[]], dtype=int)
                        rind_1 = np.array([np.arange(n1, dtype=int)])
                        rind_2 = np.array([np.arange(n2, dtype=int)])

                        if n1 > 0:
                            oval_1 = np.array([[]], dtype=int)
                            rems_1 = ch1[rind_1]
                        else:
                            oval_1 = np.array([[]], dtype=int)
                            rems_1 = np.array([[]], dtype=int)

                        if n2 > 0:
                            oval_2 = np.array([[]], dtype=int)
                            rems_2 = ch2[rind_2]
                        else:
                            oval_2 = np.array([[]], dtype=int)
                            rems_2 = np.array([[]], dtype=int)

                    else:
                        inds_1, rind_1 = self._fill_inds_data(ind_rem_data_1, i1, counts_1, mapping_1, os, False)
                        inds_2, rind_2 = self._fill_inds_data(ind_rem_data_2, i2, counts_2, mapping_2, os, True)

                        oval_1 = ch1[inds_1]
                        oval_2 = ch2[inds_2]
                        _, uinds_2 = np.unique(oval_2, axis=0, return_index=True)
                        inds_2 = inds_2[np.sort(uinds_2),]
                        rind_2 = rind_2[np.sort(uinds_2),]
                        oval_2 = ch2[inds_2]
                        rems_1 = ch1[rind_1]
                        rems_2 = ch2[rind_2]

                    cats = []
                    if os > 0:
                        cats.append(
                            (oval_1[:, np.newaxis, :] + oval_2[np.newaxis, :, :]).reshape(-1, os)
                        )
                    if os < n1:
                        pads1 = np.broadcast_to(
                            rems_1[:, np.newaxis, :],
                            rems_1.shape[:1] + rems_2.shape[:1] + (n1 - os,)
                        ).reshape(-1, n1 - os)
                        cats.append(pads1)
                    if os < n2:
                        pads2 = np.broadcast_to(
                            rems_2[np.newaxis, :, :],
                            rems_1.shape[:1] + rems_2.shape[:1] + (n2 - os,)
                        ).reshape(-1, n2 - os)
                        cats.append(pads2)

                    if len(cats) == 0:
                        # raise ValueError(ch1, ch2, os, n1, n2, counts_1, counts_2)
                        change = ()
                        if change not in new_changes:
                            new_changes[change] = []
                        new_changes[change].append((
                            ch1,
                            ch2,
                            np.array([]),
                            np.array([]),
                            np.array([]),
                            np.array([]),
                            [[]]
                        ))
                    else:
                        new_combo = np.concatenate(cats, axis=1)

                        initial_sortings = self.change_sort(new_combo)
                        new_change = nput.vector_take(new_combo, initial_sortings, shared=1)
                        nonzero_counts = np.count_nonzero(new_change, axis=1)
                        inds_1_inds, inds_2_inds = [x.flatten() for x in np.meshgrid(
                            np.arange(inds_1.shape[0], dtype=int),
                            np.arange(inds_2.shape[0], dtype=int),
                            indexing='ij'
                        )]
                        upairs = set()
                        for pad_change, initial_sort, nz, x1, x2 in zip(new_change, initial_sortings, nonzero_counts,
                                                                        inds_1_inds, inds_2_inds):

                            #TODO: should be able to cache at this depth?
                            change, change_data = self._prep_product_change_data(
                                i1, ch1, inds_1, rind_1, n1,
                                i2, ch2, inds_2, rind_2, n2,
                                os,
                                pad_change, initial_sort,
                                nz, x1, x2,
                                upairs
                            )
                            if change is None: continue

                            if change not in new_changes:
                                new_changes[change] = []
                            new_changes[change].append(change_data)

        return new_changes

    def get_expressions(self):
        raise NotImplementedError("shouldn't need this here...")

    @staticmethod
    def _freeze_product_arg(arg):
        if isinstance(arg, np.ndarray):
            return (
                'array', arg.dtype.str, tuple(arg.shape),
                tuple(arg.reshape(-1).tolist())
            )
        if isinstance(arg, (list, tuple)):
            return tuple(PerturbationTheoryTermProduct._freeze_product_arg(x) for x in arg)
        return arg

    @classmethod
    def get_poly_product_terms(cls,
                               gen1, gen2, change_1, change_2,
                               target_inds, remainder_inds, reorgs,
                               simplify=True
                               ):
        if (
                gen1.polynomial_representation != 'path'
                or gen2.polynomial_representation != 'path'
        ):
            return cls._build_poly_product_terms(
                gen1, gen2, change_1, change_2,
                target_inds, remainder_inds, reorgs,
                simplify=simplify
            )
        key = (
            gen1, gen2,
            cls._freeze_product_arg(change_1),
            cls._freeze_product_arg(change_2),
            cls._freeze_product_arg(target_inds),
            cls._freeze_product_arg(remainder_inds),
            cls._freeze_product_arg(reorgs),
            simplify,
            gen1.polynomial_representation,
            gen2.polynomial_representation
        )
        try:
            return cls._poly_product_cache[key]
        except KeyError:
            product = cls._build_poly_product_terms(
                gen1, gen2, change_1, change_2,
                target_inds, remainder_inds, reorgs,
                simplify=simplify
            )
            cls._poly_product_cache[key] = product
            return product

    @classmethod
    def _build_poly_product_terms(cls,
                                  gen1, gen2, change_1, change_2,
                                  target_inds, remainder_inds, reorgs,
                                  simplify=True
                                  ):

        log_level = Logger.LogLevel.Normal if _DEBUG_PRINT else Logger.LogLevel.MoreDebug
        logger = cls.default_logger()
        with logger.block(tag="{gen1}({ch1})[{t1}]x{gen2}({ch2})[{t2}]",
                          gen1=gen1, gen2=gen2,
                          ch1=format_quantum_change(change_1),
                          ch2=format_quantum_change(change_2),
                          t1=target_inds[0], t2=target_inds[1],
                          log_level=log_level):
            logger.log_print("remainder: {t1} {t2}", t1=remainder_inds[0], t2=remainder_inds[1], log_level=log_level)
            logger.log_print("sorting: {perm}", perm=reorgs,
                             preformatter=lambda **vars:dict(vars, perm=str(vars['perm'])),
                             log_level=log_level)

            polys_1 = gen1.get_poly_terms(change_1)
            if isinstance(polys_1, PerturbationTheoryExpressionEvaluator): polys_1 = polys_1.expr
            if nput.is_numeric(polys_1):
                if polys_1 == 0:
                    return 0
                else:
                    raise ValueError("not sure how we got a number here...")

            if _DEBUG_AUDIT:
                try:
                    polys_1.audit()
                except:
                    raise ValueError(polys_1, gen1, change_1)
            if simplify:
                polys_1 = polys_1.combine()#.sort()
                if nput.is_numeric(polys_1):
                    if polys_1 == 0:
                        return 0
                    else:
                        raise ValueError("not sure how we got a number here...")
            # we note that src_inds defines how we permute the change_2 inds so...I guess target_inds
            # tells us how we'd permute the inds of change_1?
            # if shift is not None:
            #     change_shift = change_1 + shift
            # else:
            #     change_shift = change_1
            # change_shift = [change_shift[i] for i in target_inds]

            shift_2 = [0] * max(len(change_1), len(change_2))
            for i,c in zip(*target_inds):
                # if i < len(change_1) and c < len(shift_2):
                shift_2[c] = change_1[i]

            polys_2 = gen2.get_poly_terms(change_2, shift=shift_2)
            if isinstance(polys_2, PerturbationTheoryExpressionEvaluator): polys_2 = polys_2.expr
            if nput.is_numeric(polys_2):
                if polys_2 == 0:
                    return 0
                else:
                    raise ValueError("not sure how we got a number here...")

            if _DEBUG_AUDIT:
                polys_2.audit()
            if simplify:
                polys_2 = polys_2.combine()#.sort()
            if nput.is_numeric(polys_2):
                if polys_2 == 0:
                    return 0
                else:
                    raise ValueError("not sure how we got a number here...")

            # polys_2.shift_energies(shift_2)

            logger.log_print(
                "1: {p}", p=polys_1,
                preformatter=lambda **vars: dict(vars, p=vars['p'].format_expr()),
                log_level=log_level
            )
            logger.log_print(
                "2 ([{s}]): {p}", p=polys_2, s=shift_2,
                preformatter=lambda **vars: dict(vars, p=vars['p'].format_expr()),
                log_level=log_level
            )

            # for any changes done in polys_1 that _aren't_ touched by target_inds we need to introduce a
            # baseline shift to the mul_along
            baseline = list(change_1)
            for i in target_inds[0]: baseline[i] = 0
            base_polys = polys_1.mul_along(polys_2, target_inds, remainder=remainder_inds, baseline=baseline)

            if _DEBUG_AUDIT:
                try:
                    base_polys.audit()
                except ValueError as e:
                    raise ValueError("bad polynomial found in the evaluation of {}[{}][{}]x{}[{}][{}]".format(
                        gen1, change_1, target_inds[0],
                        gen2, change_2, target_inds[1]
                    )) from e

            if reorgs is not None:
                _ = base_polys.permute(reorgs[0])
                for perm in reorgs[1:]:
                    _ += base_polys.permute(perm)
                base_polys = _

            base_polys = base_polys.strip()

            if simplify and not nput.is_numeric(base_polys):
                base_polys = base_polys.combine()#.sort()

            if nput.is_numeric(base_polys):
                if base_polys == 0:
                    return 0
                else:
                    raise ValueError("not sure how we got a number here...")

            logger.log_print(
                "{p}", p=base_polys,
                preformatter=lambda **vars: dict(vars, p=vars['p'].format_expr()),
                log_level=log_level
            )

        return base_polys

    def get_core_poly(self, changes, shift=None):

        base_changes = self.changes.get(changes, 0)
        if isinstance(base_changes, list):
            poly_changes = 0
            accumulator = (_BalancedPolynomialAccumulator()
                           if self.polynomial_representation == 'path' else None)
            for change_1, change_2, contact_1, contact_2, rem_1, rem_2, reorg in base_changes:
                # try:
                new_polys = self.get_poly_product_terms(
                    self.gen1, self.gen2, change_1, change_2,
                    [contact_1, contact_2], [rem_1, rem_2], reorg
                )
                # except:
                #     raise Exception(change_1, change_2, src_inds, target_inds, reorg)
                if nput.is_zero(new_polys): continue
                if _DEBUG_AUDIT:
                    try:
                        new_polys.audit()
                    except ValueError:
                        raise Exception("found bad polynomial in evaluating {}[{}[{}]]x{}[{}[{}]]".format(
                            self.gen1, change_1, contact_1,
                            self.gen2, change_2, contact_2
                        ))

                if accumulator is None:
                    poly_changes += new_polys
                else:
                    accumulator.append(new_polys)
            if accumulator is not None:
                poly_changes = accumulator.total()

            # if not nput.is_numeric(poly_changes):
            #
            #     log_level = Logger.LogLevel.Normal if _DEBUG_PRINT else Logger.LogLevel.MoreDebug
            #     logger = self.default_logger()
            #     logger.log_print("{p}",
            #                      p=poly_changes,
            #                      preformatter=lambda **vars: dict(vars, p=vars['p'].format_expr()),
            #                      log_level=log_level
            #                      )
        else:
            poly_changes = base_changes
        return poly_changes


class PerturbationTheoryExpressionEvaluator:
    def __init__(self, op, expr:'SqrtChangePoly', change=None, logger=None, parallelizer=None,
                 expression_factory=None):
        self.op = op
        self.expr = expr
        self._expression_factory = expression_factory
        if change is None:
            change = expr.poly_change
        self.change = change # TODO: check that expr has the same change...
        self.num_fixed = len(change)
        self.logger = Logger.lookup(logger)
        self.parallelizer = Parallelizer.lookup(parallelizer) if parallelizer is not None else None

    @property
    def expr(self):
        if '_expr' not in self.__dict__:
            self._expr = self.__dict__.pop('expr', None)
        if self._expr is None and getattr(self, '_expression_factory', None) is not None:
            self._expr = self._expression_factory()
            self._expression_factory = None
        return self._expr

    @expr.setter
    def expr(self, value):
        self._expr = value
        self._expression_factory = None
    @classmethod
    def _prune_resonant_dag(cls, expr, state_perms, changes, num_fixed, nmodes,
                           context=None):
        if (not isinstance(expr, PTTensorCoeffProductDAG) or changes is None
                or num_fixed == 0):
            return expr, b''
        from .Evaluators import _current_evaluation_plan
        plan = _current_evaluation_plan.get()
        enabled = (True if plan is None else
                   plan[1].root.options.get('prune_fixed_resonances', True))
        if not isinstance(enabled, (bool, np.bool_)):
            raise TypeError('prune_fixed_resonances must be boolean')
        if not enabled:
            return expr, b''
        permutations = np.unique(np.concatenate([
            np.asarray(perms, dtype=np.intp)[:, :num_fixed] for _, perms in state_perms
        ], axis=0), axis=0)
        if not len(permutations):
            return expr, b''
        if isinstance(changes, SideDegeneracyChanges):
            contexts = tuple(DegeneracyIdentificationContext(side or [])
                             for side in changes)
        else:
            context = context if context is not None else DegeneracyIdentificationContext(changes)
            contexts = (context, context)
        statistics = {} if Diagnostics.enabled else None
        result = expr.exclude_degenerate_changes(contexts, permutations, num_fixed,
                                                 nmodes, statistics=statistics)
        if Diagnostics.enabled:
            cls._last_dag_resonance_pruning_stats = statistics
        return result


    def __repr__(self):
        return "{}(<{}>, {})".format(type(self).__name__, self.num_fixed, self.expr)

    @staticmethod
    def _extract_coeffs(coeffs, coeff_indices, perm):
        kind = coeff_indices[1]
        if isinstance(kind,_HamiltonianCoefficient):
            row = coeffs[coeff_indices[0]] if len(coeffs)>coeff_indices[0] else ()
            coeff_tensor = _BoundHamiltonianCoefficient(kind,tuple(row[i] if i<len(row) else 0 for i in (0,1)))
        else:
            coeff_tensor = coeffs[coeff_indices[0]][kind] # order then identity
        if nput.is_zero(coeff_tensor): return 0

        idx = tuple(perm[i] for i in coeff_indices[2:])
        if len(idx) > 0:
            return coeff_tensor[idx]
        else:
            return coeff_tensor # just a number

    @classmethod
    def _eval_raw_poly(cls,
                       perm_substates, poly, change,
                       baseline_shift, pows,
                       verbose, logger,
                       path_node_cache=None, path_axis_cache=None
                       ):
        log_level = Logger.LogLevel.Normal if verbose else Logger.LogLevel.Debug
        if isinstance(poly, PolyPath):
            poly_evals = []
            for substates in perm_substates:
                if baseline_shift is not None:
                    raise NotImplementedError('sorry')
                if substates.shape[-1] == 0:
                    sqrt_factor = np.ones(substates.shape[0])
                elif isinstance(poly, _CompiledPolynomial) and poly.occupation_shifts is not None:
                    shifts_sqrts = [np.prod(n[:, np.newaxis]+offset[np.newaxis, :], axis=-1)
                                    for n, offset in zip(substates.T, poly.occupation_shifts)]
                    sqrt_factor = np.sqrt(np.prod(shifts_sqrts, axis=0))
                else:
                    shifts_sqrts = [
                        np.prod(
                            n[:, np.newaxis] + np.arange(
                                delta + 1 if delta < 0 else 1,
                                1 if delta < 0 else delta + 1
                            )[np.newaxis, :],
                            axis=-1
                        )
                        for n, delta in zip(substates.T, change)
                    ]
                    sqrt_factor = np.sqrt(np.prod(shifts_sqrts, axis=0))

                node_cache = (
                    path_node_cache.new_path_node_cache()
                    if isinstance(path_node_cache, _MaterializedEvaluationCache) else
                    path_node_cache
                )
                value = poly.evaluate_polynomial(
                    substates,
                    node_cache=node_cache,
                    axis_cache=path_axis_cache
                )
                poly_evals.append(value * sqrt_factor)
            if len(poly_evals) == 1:
                return poly_evals[0][:, np.newaxis]
            return np.moveaxis(np.array(poly_evals), 0, 1)

        if isinstance(poly, ProductPTPolynomialSum):
            subvals = [
                cls._eval_raw_poly(
                    perm_substates, p, change, baseline_shift, pows, verbose, logger,
                    path_node_cache=path_node_cache,
                    path_axis_cache=path_axis_cache
                )
                for p in poly.polys
            ]
            return poly.prefactor * np.sum(subvals, axis=0)

        if verbose:
            # logging_poly = poly.constant_rescale()
            logging_poly = poly
            # p = poly
            logger.log_print("prefactor: {p}{s}",
                             p=logging_poly.prefactor,
                             n=logging_poly.steps,
                             preformatter=lambda **vars:dict(
                                 vars,
                                 s=(
                                     ""
                                        if vars['n'] == 0 else
                                     "/Sqrt[2]"
                                        if vars['n'] == 1 else
                                     "/{}".format(2**(vars['n']//2))
                                        if vars['n']%2 == 0 else
                                     "/({}*Sqrt[2])".format(
                                         2**(vars['n']//2)
                                     )
                                 )
                             ),
                             log_level=log_level)
            logger.log_print("coeffs: {c}",
                             c=logging_poly.coeffs,
                             preformatter=lambda **vars:dict(vars, c=[np.round(c, 3) for c in vars['c']]),
                             log_level=log_level
                             )

        poly_evals = []
        for substates in perm_substates:

            if baseline_shift is not None:
                raise NotImplementedError('sorry')
                substates = substates + baseline_shift[np.newaxis]

            if substates.shape[-1] == 0:
                ct = np.ones(substates.shape[0])
            else:
                # if len(change) != len(poly.coeffs): raise ValueError(change, poly)
                poly_factor = np.prod(
                    [
                        # np.dot(c[np.newaxis, :], s[:len(c)])
                        # for s, c in zip(substate_powers, poly.coeffs)
                        np.dot(c[np.newaxis, :],
                               pows[:len(c), s]
                               # np.power(s[np.newaxis, :], np.arange(len(c))[:, np.newaxis])
                               )
                        for s, c in zip(substates.T, poly.coeffs)
                    ],
                    axis=0
                )
                if not nput.is_numeric(poly_factor): poly_factor = poly_factor[0] # we padded for the dot products
                shifts_sqrts = [
                    np.prod(
                        n[:, np.newaxis] + np.arange(
                            delta + 1 if delta < 0 else 1,
                            1 if delta < 0 else delta + 1
                        )[np.newaxis, :],
                        axis=-1
                    )
                    for n, delta in zip(substates.T, change)
                ]
                sqrt_factor = np.sqrt(np.prod(shifts_sqrts, axis=0)) / np.sqrt(2)**poly.steps

                ct = poly.prefactor * poly_factor * sqrt_factor
            poly_evals.append(ct)
        if len(poly_evals) == 1:
            return poly_evals[0][:, np.newaxis]
        else:
            return np.moveaxis(np.array(poly_evals), 0, 1)
        # return ct


    @classmethod
    def _resolve_poly_evaluation_cache(cls, cache, poly_key):
        batch_cache = cache if isinstance(cache, BoundedCache) else None
        materialized_cache = (
            cache if isinstance(cache, _MaterializedEvaluationCache) else None
        )
        if materialized_cache is not None:
            exact_cache = materialized_cache.exact_cache
            if isinstance(exact_cache, BoundedCache):
                value_cache = exact_cache.view('polynomial_values')
                batch_cache = exact_cache
            else:
                if poly_key not in exact_cache:
                    exact_cache[poly_key] = {}
                value_cache = exact_cache[poly_key]
            path_node_cache = materialized_cache
            path_axis_cache = materialized_cache.path_axis_cache
        elif batch_cache is None:
            if poly_key not in cache:
                cache[poly_key] = {}
            value_cache = cache[poly_key]
            path_node_cache = path_axis_cache = None
        else:
            value_cache = batch_cache.view('polynomial_values')
            path_node_cache = batch_cache.view('polypath_nodes')
            path_axis_cache = batch_cache.view('polypath_axes')
        return (
            value_cache, path_node_cache, path_axis_cache,
            batch_cache, materialized_cache
        )

    @classmethod
    def _lookup_poly_evaluation_cache(cls,
                                      poly_key, value_cache,
                                      tuple_states, perm_substates,
                                      batch_cache, materialized_cache):
        eval_pos = []
        eval_keys = []
        eval_states = []
        vals = np.zeros((len(perm_substates[0]), len(perm_substates)), dtype=float)
        for j, substates in enumerate(tuple_states):
            for i, state_key in enumerate(substates):
                value_key = (
                    (poly_key, state_key) if batch_cache is not None else state_key
                )
                cached_value = value_cache.get(
                    value_key, BoundedCache._missing
                )
                if cached_value is BoundedCache._missing:
                    if materialized_cache is not None:
                        if Diagnostics.enabled:
                            materialized_cache.exact_misses += 1
                    eval_pos.append((i, j))
                    eval_keys.append(value_key)
                    eval_states.append(perm_substates[j][i])
                else:
                    if materialized_cache is not None:
                        if Diagnostics.enabled:
                            materialized_cache.exact_hits += 1
                    vals[i, j] = cached_value
        return vals, eval_pos, eval_keys, eval_states

    @staticmethod
    def _stack_poly_evaluation_states(eval_states):
        if len(eval_states[0]) == 0:
            return np.empty((len(eval_states), 0), dtype=eval_states[0].dtype)
        return np.fromiter(
            eval_states,
            count=len(eval_states),
            dtype=np.dtype((eval_states[0].dtype, (len(eval_states[0]),)))
        )

    @staticmethod
    def _store_poly_evaluation_cache(value_cache, vals, eval_keys, eval_pos, evals):
        for key, value, (i, j) in zip(eval_keys, evals[:, 0], eval_pos):
            value_cache[key] = value
            vals[i, j] = value

    @classmethod
    def _eval_poly(cls, cache, tuple_states, perm_substates, pows,
                   poly, change, baseline_shift,
                   verbose, logger, block_key=None):
        # Complete (polynomial, change, baseline, state) values remain the first
        # cache level.  _MaterializedEvaluationCache adds a second, batch-local
        # level for internal PolyPath nodes and axes shared by nonidentical root
        # polynomials. Both levels can be owned by a numerical plan, with
        # entry/estimated-byte bounds and exact concrete occupation identities.
        # Expression coverage and output accumulation remain request-local.

        poly_key = (
            poly,
            None if change is None else tuple(change),
            None if baseline_shift is None else tuple(baseline_shift)
        )
        materialized_cache = (
            cache if isinstance(cache, _MaterializedEvaluationCache) else None
        )
        block_value_key = None
        if materialized_cache is not None and block_key is not None:
            block_value_key = (poly_key, block_key)
            block_value = materialized_cache.lookup_polynomial_block(
                block_value_key
            )
            if block_value is not BoundedCache._missing:
                return block_value
        (
            value_cache, path_node_cache, path_axis_cache,
            batch_cache, materialized_cache
        ) = cls._resolve_poly_evaluation_cache(cache, poly_key)
        vals, eval_pos, eval_keys, eval_states = cls._lookup_poly_evaluation_cache(
            poly_key, value_cache,
            tuple_states, perm_substates,
            batch_cache, materialized_cache
        )
        if len(eval_states) > 0:
            substates = cls._stack_poly_evaluation_states(eval_states)
            evals = cls._eval_raw_poly(
                [substates], poly, change, baseline_shift, pows, verbose, logger,
                path_node_cache=path_node_cache,
                path_axis_cache=path_axis_cache
            )
            cls._store_poly_evaluation_cache(
                value_cache, vals, eval_keys, eval_pos, evals
            )
        if block_value_key is not None:
            materialized_cache.store_polynomial_block(block_value_key, vals)
        return vals
    @classmethod
    def _compute_energy_weights(cls, energy_changes, perm_freqs):

        #TODO: find a way to avoid recomputing the same contributions repeatedly...

        if np.size(energy_changes) == 0:
            return np.ones(np.shape(perm_freqs)[:-1], dtype=np.result_type(perm_freqs, float))
        echange = np.prod(np.dot(perm_freqs, energy_changes[:, 1:].T), axis=-1)
        return echange

    @staticmethod
    def _report_energy_denominator_warning(message, category, *, logger=None):
        # The native log must retain this breakdown even with verbose=False
        # or a Python warning filter that suppresses repeated stderr messages.
        if logger is not None:
            Logger.lookup(logger).log_print(
                'WARNING: {warning_text}', warning_text=message,
                log_level=Logger.LogLevel.Warnings)
        warnings.warn(message, category, stacklevel=4)

    @staticmethod
    def _denominator_state_pairs(energy_changes, positions, context):
        """Recover exact endpoints only where native path metadata retains them.

        Ordinary relative-energy rows do not identify a unique intermediate
        occupation pair. Preserve their source state and delta without guessing.
        """
        payload = context()
        state, permutations, frequencies, baseline = payload[:4]
        projection_bindings = payload[4] if len(payload)>4 else ()
        state = np.asarray(state, dtype=int).reshape((-1, np.shape(state)[-1]))
        if len(state) != 1:
            return dict(unavailable='multiple source occupations share this contraction')
        source = state[0]
        permutations = np.asarray(permutations, dtype=int)
        frequencies = np.asarray(frequencies)
        rows = list(energy_changes)
        examples = []
        for position in positions[:4]:
            modes = permutations[int(position)]
            pairs, unresolved = [], []
            def shifted(delta):
                result = source.copy()
                result[modes[:len(delta)]] += np.asarray(delta, dtype=int)
                return result
            index = 0
            while index < len(rows):
                row = rows[index]; factor = row[0]; delta = np.asarray(row[1:], dtype=int)
                if getattr(factor, '_conditional_endpoint', False):
                    other = rows[index+1]
                    if factor.kind == 'resolvent':
                        shift = np.zeros(len(delta), dtype=int)
                        if baseline is not None:
                            shift[:min(len(shift),len(baseline))] = baseline[:len(shift)]
                        left, right = shifted(delta+shift), shifted(np.asarray(other[1:])+shift)
                        gap = float(frequencies[int(position),:len(delta)] @ (delta-np.asarray(other[1:])))
                        if np.all(left >= 0) and np.all(right >= 0):
                            pairs.append(dict(factor=index, kind='conditional resolvent',
                                left_state=tuple(int(n) for n in left), right_state=tuple(int(n) for n in right),
                                energy_difference='E(left)-E(right)', gap_hartree=gap))
                    index += 2
                    continue
                if getattr(factor, '_state_projection_factor', False):
                    if factor.kind == 'R':
                        if factor.slot >= len(projection_bindings):
                            unresolved.append(dict(factor=index,reason='numerical projection binding unavailable'))
                            index += 1
                            continue
                        _,full_freqs,reference,_ = projection_bindings[factor.slot]
                        left, right = np.asarray(reference), shifted(delta)
                        gap = float((left-right) @ np.asarray(full_freqs))
                        if np.all(right >= 0):
                            pairs.append(dict(factor=index, kind='projected resolvent',
                                left_state=tuple(int(n) for n in left), right_state=tuple(int(n) for n in right),
                                energy_difference='E(left)-E(right)', gap_hartree=gap))
                else:
                    unresolved.append(dict(factor=index, side=int(factor),
                        physical_modes_1based=tuple(int(n)+1 for n in modes[:len(delta)]),
                        energy_change=tuple(int(n) for n in delta),
                        reason='relative energy row has no unique intermediate endpoints'))
                index += 1
            examples.append(dict(position=int(position), source_state=tuple(int(n) for n in source),
                                 state_pairs=pairs, unresolved_factors=unresolved))
        return examples

    @staticmethod
    def _format_small_denominator_context(evaluation_context, positions,
                                          energy_factors, energy_changes=None,
                                          tuple_states=None,
                                          permutation_substates=None,
                                          polynomial_factor=None):
        """Build a bounded diagnostic payload for a small denominator warning."""
        details = {
            'evaluation': ({key:value for key,value in evaluation_context.items()
                            if key != '_state_pair_context'}
                           if isinstance(evaluation_context,dict) else evaluation_context),
            'positions': tuple(int(position) for position in positions[:4]),
            'denominators': tuple(
                float(energy_factors[position]) for position in positions[:4]
            )
        }
        if (energy_changes is not None and isinstance(evaluation_context, dict)
                and evaluation_context.get('_state_pair_context') is not None):
            details['denominator_state_pairs'] = (
                PerturbationTheoryExpressionEvaluator._denominator_state_pairs(
                    energy_changes, positions, evaluation_context['_state_pair_context']))
        if energy_changes is not None:
            details['energy_changes'] = tuple(
                tuple(int(value) for value in change)
                for change in np.asanyarray(energy_changes)
            )
        if tuple_states is not None:
            details['states'] = tuple(
                [
                    tuple(state.tolist())
                    if isinstance(state, _HashableNumPyIndex) else state
                    for state in tuple_states[position]
                ]
                for position in positions[:4]
            )
        if permutation_substates is not None:
            details['permutation_states'] = tuple(
                np.asanyarray(permutation_substates[position]).tolist()
                for position in positions[:4]
            )
        if polynomial_factor is not None:
            numerator = np.asanyarray(polynomial_factor)
            details['numerators'] = tuple(
                np.asanyarray(numerator[..., position]).tolist()
                for position in positions[:4]
            )
        if len(positions) > 4:
            details['additional_positions'] = len(positions) - 4
        return details

    @classmethod
    def _divide_polynomial_by_energy(cls, polynomial_factor, energy_factors,
                                     zero_cutoff, *, energy_changes=None,
                                     evaluation_context=None,
                                     tuple_states=None,
                                     permutation_substates=None, logger=None):
        """Prune evaluated numerators, then divide by numerical denominators.

        A symbolic energy change can vanish for a particular frequency vector.
        The corresponding numerator is therefore evaluated before deciding
        whether the term is harmless.  This is deliberately an evaluator-side
        decision: symbolic prefactors do not contain enough state information
        to make it safely.
        """
        numerator = np.asanyarray(polynomial_factor)
        denominator = np.asanyarray(energy_factors)
        cutoff = float(zero_cutoff)

        small_numerators = np.abs(numerator) < cutoff
        if np.any(small_numerators):
            numerator = np.where(small_numerators, 0, numerator)

        if numerator.ndim == 0:
            active_numerators = np.full(
                denominator.shape, numerator != 0, dtype=bool
            )
        else:
            active_numerators = np.any(
                numerator != 0,
                axis=tuple(range(numerator.ndim - 1))
            )

        small_denominators = np.abs(denominator) < cutoff
        if np.any(small_denominators):
            positions = np.flatnonzero(small_denominators)
            details = cls._format_small_denominator_context(
                evaluation_context, positions, denominator,
                energy_changes=energy_changes
            )
            cls._report_energy_denominator_warning(
                "energy denominator below zero_cutoff={}: {}".format(
                    cutoff, details
                ),
                SmallEnergyDenominatorWarning,
                logger=logger
            )

            surviving_positions = positions[active_numerators[positions]]
            if len(surviving_positions) > 0:
                details = cls._format_small_denominator_context(
                    evaluation_context, surviving_positions, denominator,
                    energy_changes=energy_changes,
                    tuple_states=tuple_states,
                    permutation_substates=permutation_substates,
                    polynomial_factor=numerator
                )
                cls._report_energy_denominator_warning(
                    "term survived numerator pruning at a small energy "
                    "denominator: {}".format(details),
                    SurvivingEnergyDenominatorWarning,
                    logger=logger
                )

        # A pruned numerator over an exact zero is exactly a discarded
        # contribution, not NaN.  Keep the original small denominator only
        # where a numerator survived so the diagnostic value remains visible.
        safe_denominator = np.where(
            small_denominators & np.logical_not(active_numerators),
            1,
            denominator
        )
        with np.errstate(divide='ignore', invalid='ignore'):
            scaled = numerator / safe_denominator[np.newaxis, :]
        # Ratio pruning belongs here, after both the concrete numerator and
        # denominator are known.  Doing it while constructing symbolic terms
        # can incorrectly remove a small numerator amplified by a denominator.
        return np.where(np.abs(scaled) < cutoff, 0, scaled)

    @classmethod
    def _divide_polynomial_block_by_energy(cls, polynomial_factors,
                                           energy_factors, zero_cutoff, *,
                                           row_diagnostics=None, logger=None):
        """Apply the scalar denominator policy to a homogeneous row block.

        ``polynomial_factors`` has shape ``(row, state, permutation)`` and
        ``energy_factors`` has shape ``(row, permutation)``.  Diagnostic state
        is supplied lazily because valid production calculations should never
        pay to format operator/state payloads for the warning-only path.
        """
        numerator = np.asanyarray(polynomial_factors)
        denominator = np.asanyarray(energy_factors)
        if numerator.ndim != 3 or denominator.ndim != 2:
            raise ValueError(
                "block division requires (row,state,permutation) and "
                "(row,permutation) arrays, got {} and {}".format(
                    numerator.shape, denominator.shape
                )
            )
        if numerator.shape[0] != denominator.shape[0] or numerator.shape[-1] != denominator.shape[-1]:
            raise ValueError(
                "incompatible polynomial/energy block shapes {} and {}".format(
                    numerator.shape, denominator.shape
                )
            )

        cutoff = float(zero_cutoff)
        numerator = np.where(np.abs(numerator) < cutoff, 0, numerator)
        active_numerators = np.any(numerator != 0, axis=1)
        small_denominators = np.abs(denominator) < cutoff

        if np.any(small_denominators):
            affected_rows = np.flatnonzero(np.any(small_denominators, axis=1))
            for row in affected_rows:
                diagnostics = (
                    {} if row_diagnostics is None else
                    row_diagnostics(int(row))
                    if callable(row_diagnostics) else
                    row_diagnostics[row]
                )
                positions = np.flatnonzero(small_denominators[row])
                details = cls._format_small_denominator_context(
                    diagnostics.get('evaluation_context'),
                    positions,
                    denominator[row],
                    energy_changes=diagnostics.get('energy_changes')
                )
                cls._report_energy_denominator_warning(
                    "energy denominator below zero_cutoff={}: {}".format(
                        cutoff, details
                    ),
                    SmallEnergyDenominatorWarning,
                    logger=logger
                )

                surviving_positions = positions[
                    active_numerators[row, positions]
                ]
                if len(surviving_positions) > 0:
                    details = cls._format_small_denominator_context(
                        diagnostics.get('evaluation_context'),
                        surviving_positions,
                        denominator[row],
                        energy_changes=diagnostics.get('energy_changes'),
                        tuple_states=diagnostics.get('tuple_states'),
                        permutation_substates=diagnostics.get(
                            'permutation_substates'
                        ),
                        polynomial_factor=numerator[row]
                    )
                    cls._report_energy_denominator_warning(
                        "term survived numerator pruning at a small energy "
                        "denominator: {}".format(details),
                        SurvivingEnergyDenominatorWarning,
                        logger=logger
                    )

        safe_denominator = np.where(
            small_denominators & np.logical_not(active_numerators),
            1,
            denominator
        )
        with np.errstate(divide='ignore', invalid='ignore'):
            scaled = numerator / safe_denominator[:, np.newaxis, :]
        return np.where(np.abs(scaled) < cutoff, 0, scaled)

    @staticmethod
    def _compute_coefficient_prefactor_products(perms, cinds, tensors, zero_cutoff):
        products = np.ones((len(tensors), len(perms)), dtype=float)
        good_perms = np.full((len(tensors), len(perms)), True)
        for coefficient_pos, (_, indices) in enumerate(cinds):
            for tensor_pos, tensor_list in enumerate(tensors):
                tensor = tensor_list[coefficient_pos]
                if len(indices) == 0:
                    base_value = tensor
                    if abs(base_value) < zero_cutoff:
                        good_perms[:] = False
                        break
                    products[tensor_pos] *= base_value
                else:
                    active = np.flatnonzero(good_perms[tensor_pos])
                    if len(active) == 0:
                        continue
                    active_perms = perms[active]
                    tensor_index = tuple(
                        active_perms[:, index] for index in indices
                    )
                    try:
                        base_values = tensor[tensor_index]
                    except Exception:
                        raise ValueError(tensor, active_perms, tensor_index)
                    products[tensor_pos, active] *= base_values
                    good_perms[tensor_pos, active] = (
                        np.abs(products[tensor_pos, active]) >= zero_cutoff
                    )
        return products, np.any(good_perms, axis=0)

    @staticmethod
    def _store_coefficient_prefactor_products(products, good_perms, cind_pos,
                                               eval_perms, eval_cinds, prefactors,
                                               start):
        active = np.flatnonzero(good_perms)
        next_pos = start + len(active)
        eval_perms[start:next_pos] = active
        eval_cinds[start:next_pos] = cind_pos
        prefactors[start:next_pos] = products[:, active].T
        return next_pos

    @classmethod
    def _get_prefacs(cls, perms, cinds_remapped, ctensors, counts_cache, facs, zero_cutoff):
        # compute prefactor products with fast exiting
        nevals = 0
        eval_perms = np.empty(len(perms) * len(cinds_remapped), dtype=int)
        eval_cinds = np.empty(len(perms) * len(cinds_remapped), dtype=int)
        prefactors = np.empty((len(perms) * len(cinds_remapped), len(ctensors[0])), dtype=float)
        for j, cinds in enumerate(cinds_remapped):
            cur_hits = counts_cache.get(cinds, 0)
            if cur_hits < facs[len(cinds)]:
                counts_cache[cinds] = cur_hits + 1
                products, good_perms = cls._compute_coefficient_prefactor_products(
                    perms, cinds, ctensors[j], zero_cutoff
                )
                nevals = cls._store_coefficient_prefactor_products(
                    products, good_perms, j,
                    eval_perms, eval_cinds, prefactors,
                    nevals
                )

        return eval_perms[:nevals], eval_cinds[:nevals], prefactors[:nevals]

    @classmethod
    def _eval_subcontrib(cls,
                         echanges, perm_freqs,
                         poly_cache, energy_cache,
                         tuple_states, perm_substates, pows,
                         polys, change, baseline_shift, verbose, logger,
                         block_key=None, zero_cutoff=None,
                         evaluation_context=None
                         ):
        # if echanges not in energy_cache:
        #     energy_cache[echanges] = np.array(echanges)
        # echanges = energy_cache[echanges]
        echanges = np.asanyarray(echanges)

        energy_factors = cls._compute_energy_weights(echanges, perm_freqs)
        # TODO: check size of energy factor
        poly_factor = cls._eval_poly(poly_cache,
                                     tuple_states, perm_substates, pows,
                                     polys, change, baseline_shift, verbose, logger,
                                     block_key=block_key)
        if zero_cutoff is None:
            zero_cutoff = cls.default_zero_cutoff
        # Even an entirely pruned numerator goes through the denominator check:
        # users asked to be told that the small energy term existed, not only
        # when it happened to survive this particular state evaluation.
        scaled_contrib = cls._divide_polynomial_by_energy(
            poly_factor, energy_factors, zero_cutoff,
            energy_changes=echanges,
            evaluation_context=evaluation_context,
            tuple_states=tuple_states,
            permutation_substates=perm_substates, logger=logger
        )
        return scaled_contrib, energy_factors, poly_factor
    @classmethod
    def _test_degs(cls, subchanges, echanges, perm_subsets, only_deg):
        subsubchange_test = subchanges.get(echanges, None) if subchanges is not None else None
        if subsubchange_test is None:
            base = True if not only_deg else False
            return [base] * len(perm_subsets)
        else:
            if only_deg:
                return [subsubchange_test(perm) for perm in perm_subsets]
            else:
                return [not subsubchange_test(perm) for perm in perm_subsets]
    @classmethod
    def _eval_perm_core(cls,
                        expr, state, tuple_states, perm_substates, which_perms,
                        change, baseline_shift,
                        prefacs, perm_freqs, state_block_key,
                        pows, key, perm_subsets, degenerate_changes, only_degenerate_terms,
                        poly_cache, energy_cache,
                        verbose, logger, log_level, log_scaling,
                        zero_cutoff, evaluation_context=None
                        ):
        subexpr = expr.terms[key]
        conv_subsets = None
        use_degs_left, use_degs_right, deg_join = only_degenerate_terms
        include_degs = (use_degs_left or use_degs_right)
        if isinstance(subexpr, PTEnergyChangeProductSum):
            subcontrib = np.zeros([len(state), len(perm_subsets)])

            if degenerate_changes is not None:
                left_changes, right_changes = [d.get(key, {}) for d in degenerate_changes]
                # all_subchanges = [
                #     degenerate_changes[i].get(key, None)
                #     for i in which_perms
                # ]
            else:
                left_changes = right_changes = None
                # all_subchanges = [None] * len(which_perms)

            for echanges, polys in subexpr.terms.items():
                if not echanges:
                    # A unit denominator has no resonance or small-denominator
                    # test. Evaluate exactly like the original bare numerator,
                    # including its unchanged numerical-pruning behavior.
                    if not include_degs:
                        subcontrib += cls._eval_poly(
                            poly_cache, tuple_states, perm_substates, pows,
                            polys, change, baseline_shift, verbose, logger,
                            block_key=('full', state_block_key))
                    continue
                if any(getattr(row[0], '_state_projection_factor', False) for row in echanges):
                    good, denominators = projected_path_factors(
                        echanges, state, perm_subsets, perm_freqs, baseline_shift,
                        projection_bindings=projection_inputs(evaluation_context)
                    )
                    if not np.any(good):
                        continue
                    numerator = cls._eval_poly(
                        poly_cache, tuple_states, perm_substates, pows,
                        polys, change, baseline_shift, verbose, logger,
                        block_key=('full', state_block_key)
                    )
                    numerator = np.where(good[np.newaxis, :], numerator, 0.)
                    denominators = np.where(good, denominators, 1.)
                    scaled = cls._divide_polynomial_by_energy(
                        numerator, denominators, zero_cutoff,
                        evaluation_context={
                            'operator': repr(evaluation_context),
                            '_state_pair_context': lambda: (state, perm_subsets, perm_freqs, baseline_shift, projection_inputs(evaluation_context))},
                        energy_changes=echanges, logger=logger
                    )
                    subcontrib += scaled
                    continue
                term_context = {
                    'operator': repr(expr if evaluation_context is None else evaluation_context),
                    'coefficient_key': key,
                    'change': None if change is None else tuple(change),
                    'baseline_shift': (
                        None if baseline_shift is None else tuple(baseline_shift)
                    )
                }
                # if all(any(e != 0 for e in ech) for ech in echanges):
                left_matches = cls._test_degs(left_changes, echanges, perm_subsets, use_degs_left)
                right_matches = cls._test_degs(right_changes, echanges, perm_subsets, use_degs_right)
                good_perms = np.where(deg_join(left_matches, right_matches))[0]

                # good_perms = []
                # for idx, subchanges in enumerate(all_subchanges):
                #     if subchanges is not None and echanges in subchanges:
                #         subsubchange_test = subchanges[echanges]
                #         # deg_match = subsubchange_test(perm_subsets[idx])
                #         if only_degenerate_terms:
                #             if deg_match: good_perms.append(idx)
                #         else:
                #             if not deg_match: good_perms.append(idx)
                #     elif not only_degenerate_terms:
                #         good_perms.append(idx)
                if len(good_perms) == len(perm_subsets):
                    good_perms = None
                    good_freqs = perm_freqs
                    good_perm_substates = perm_substates
                    good_tuple_states = tuple_states
                else:
                    good_perms = tuple(int(p) for p in good_perms)
                    good_freqs = perm_freqs[good_perms,]
                    good_perm_substates = perm_substates[good_perms,]
                    good_tuple_states = [tuple_states[p] for p in good_perms]

                poly_block_key = (
                    ('full', state_block_key)
                    if good_perms is None else
                    ('take', state_block_key, good_perms)
                )

                term_context['_state_pair_context'] = lambda: (
                    state, np.asarray(perm_subsets) if good_perms is None else
                    np.asarray(perm_subsets)[list(good_perms)], good_freqs, baseline_shift,
                    projection_inputs(evaluation_context))
                if len(good_perm_substates) > 0:
                    if verbose:
                        with logger.block(tag="{e} * {p}",
                                          preformatter=lambda **vars: dict(
                                              vars,
                                              e="E-[" + PTEnergyChangeProductSum.format_energy_prod_key(vars['e']) + "]",
                                              p=vars['p'].format_expr()
                                          ),
                                          pf=subexpr.prefactor,
                                          e=echanges,
                                          p=polys,
                                          log_level=log_level):
                                scaled_contrib, energy_factors, poly_factor = cls._eval_subcontrib(
                                    echanges, good_freqs,
                                    poly_cache, energy_cache,
                                    good_tuple_states, good_perm_substates, pows,
                                    polys, change, baseline_shift, verbose, logger,
                                    block_key=poly_block_key,
                                    zero_cutoff=zero_cutoff,
                                    evaluation_context=term_context
                                )
                                if not nput.is_zero(scaled_contrib):
                                    logger.log_print("engs: {ef}", ef=energy_factors.squeeze(), log_level=log_level)
                                    logger.log_print("poly: {pf}", pf=poly_factor.squeeze(), log_level=log_level)
                                    logger.log_print("sc: {pf}",
                                                     pf=(prefacs[:, np.newaxis, :] * scaled_contrib[np.newaxis, : :]).squeeze() * log_scaling,
                                             log_level=log_level)
                    else:
                        scaled_contrib, energy_factors, poly_factor = cls._eval_subcontrib(
                            echanges, good_freqs,
                            poly_cache, energy_cache,
                            good_tuple_states, good_perm_substates, pows,
                            polys, change, baseline_shift, verbose, logger,
                            block_key=poly_block_key,
                            zero_cutoff=zero_cutoff,
                            evaluation_context=term_context
                        )

                    if good_perms is None:
                        subcontrib += scaled_contrib
                    else:
                        subcontrib[:, good_perms] += scaled_contrib
            subcontrib *= subexpr.prefactor
        elif (not include_degs) and isinstance(subexpr, (ProductPTPolynomial, ProductPTPolynomialSum)):
            if verbose:
                with logger.block(tag="{p}",
                                  p=subexpr,
                                  preformatter=lambda **vars: dict(vars, p=vars['p'].format_expr()),
                                  log_level=log_level):
                    subcontrib = cls._eval_poly(poly_cache, tuple_states, perm_substates, pows,
                                                subexpr, change, baseline_shift, verbose, logger,
                                                block_key=('full', state_block_key))
            else:
                subcontrib = cls._eval_poly(poly_cache, tuple_states, perm_substates, pows,
                                            subexpr, change, baseline_shift, verbose, logger,
                                            block_key=('full', state_block_key))

                # if not nput.is_zero(subcontrib):
                #     subcontrib = subcontrib[:, np.newaxis]
        else:
            raise ValueError("degenerate terms requested on {}".format(subexpr.format_expr()))

        return subcontrib

    @staticmethod
    def _gather_state_permutations(state, sub_perms):
        if (isinstance(state,np.ndarray) and state.ndim==2
                and isinstance(sub_perms,np.ndarray) and sub_perms.ndim==2
                and sub_perms.dtype.kind in 'iu'):
            return np.moveaxis(state[:,sub_perms],0,1)
        return np.moveaxis(nput.vector_take(state, sub_perms), 0, 1)

    @staticmethod
    def _tupleize_state_permutations(perm_substates):
        return [[tuple(state) for state in state_block] for state_block in perm_substates]

    @classmethod
    def _index_state_permutations(cls, perm_substates):
        threshold = cls.hashable_numpy_index_threshold
        return [
            [
                _HashableNumPyIndex(state)
                if state.size > threshold else tuple(state)
                for state in state_block
            ]
            for state_block in perm_substates
        ]

    @staticmethod
    def _gather_frequency_permutations(freqs, sub_perms):
        if (isinstance(freqs,np.ndarray) and freqs.ndim==1
                and isinstance(sub_perms,np.ndarray) and sub_perms.ndim==2
                and sub_perms.dtype.kind in 'iu'):
            return freqs[sub_perms]
        return nput.vector_take(freqs, sub_perms)

    @classmethod
    def _get_state_perms(cls, state_idx, state, freqs, fixed, subset, sub_perms, take_cache,
                         mask_pos, full_set):
        # key = (state_idx, mask_pos, full_set)
        # if key not in take_cache:
        #     take_cache[key] = {}
        # if 'state_dat' not in take_cache[key]:
        #     arr, inds = nput.vector_take(state, sub_perms, return_spec=True)
        #     take_cache[key]['state_dat'] = (arr, inds)
        # # print(fixed, subset, len(state), inds[-1].shape, sub_perms.shape)
        # # arr, inds = nput.vector_take(state, sub_perms, return_spec=True)
        # arr, inds = take_cache[key]['state_dat']
        # inds = inds[:-1] + (np.broadcast_to(sub_perms[np.newaxis], (len(state),) + sub_perms.shape),)
        # perm_substates = np.moveaxis(arr[inds], 0, 1)
        if (
                cls.use_hashable_numpy_mask_indices
                and len(mask_pos) > cls.hashable_numpy_index_threshold
        ):
            mask_key = (
                mask_pos
                if isinstance(mask_pos, _HashableNumPyIndex) else
                _HashableNumPyIndex(mask_pos)
            )
        else:
            mask_key = (
                mask_pos if isinstance(mask_pos, tuple) else
                tuple(int(position) for position in mask_pos)
            )
        cache_key = ('state_permutations', state_idx, mask_key, full_set)
        cached = take_cache.get(cache_key, BoundedCache._missing)
        if cached is not BoundedCache._missing:
            return cached

        perm_substates = cls._gather_state_permutations(state, sub_perms)
        tuple_states = (
            cls._index_state_permutations(perm_substates)
            if cls.use_hashable_numpy_state_indices else
            cls._tupleize_state_permutations(perm_substates)
        )
        perm_freqs = cls._gather_frequency_permutations(freqs, sub_perms)
        block_identity = _StatePermutationBlockIdentity(tuple_states)
        cached = perm_substates, tuple_states, perm_freqs, block_identity
        take_cache[cache_key] = cached

        # if 'freq_dat' not in take_cache[key]:
        #     arr, inds = nput.vector_take(freqs, sub_perms, return_spec=True)
        #     take_cache[key]['freq_dat'] = (arr, inds)
        # arr, inds = take_cache[key]['freq_dat']
        # inds = inds[:-1] + (np.broadcast_to(sub_perms[np.newaxis], inds[-1].shape),)
        # perm_freqs = arr[inds]

        return cached
    @classmethod
    def _eval_perm(cls,
                   expr, change, baseline_shift,
                   subset, state_perms, all_perms, perm_map,
                   freqs, cind_specs, ctensors,
                   num_fixed, degenerate_changes, only_degenerate_terms,
                   zero_cutoff,
                   counts_cache, poly_cache, take_cache, energy_cache,
                   facs, pows, split_spec,
                   verbose, logger, log_scaled,
                   evaluation_context=None
                   ):
        log_level = Logger.LogLevel.Normal if verbose else Logger.LogLevel.Debug
        log_scaling = 219475.6 if log_scaled else 1
        #TODO: use combinatorics to skip some duplicate evaluations over permutations

        full_set = tuple(range(num_fixed)) + subset

        free = subset
        fixed = num_fixed

        cinds_remapped = [
            tuple(
                # don't need to symmetrize now that we
                # just loop over combinations (which are ordered)
                (ci[:2], tuple(free[j - fixed] if j >= fixed else j for j in ci[2:]))
                for ci in cinds
            )
            for cinds in cind_specs
        ]

        # we compute the unique prefactors and corresponding cinds
        # that can come out of the total space of perms, then
        # group these by cinds, and figure out how these reduced
        # positions map onto the total state-by-state space

        eval_perms, eval_cinds, prefactors = cls._get_prefacs(
            all_perms, cinds_remapped, ctensors, counts_cache, facs, zero_cutoff
        )
        contrib = [
            np.zeros([prefactors.shape[1], len(perms)])
            for state,perms in state_perms
        ]
        if prefactors.shape[0] == 0: return contrib

        # to figure out how the eval_perms indices map onto the state-by-state set
        # we build a mask of the size of the total unique set, then figure out how these
        # map back onto the total space

        all_perms_mask = np.full(len(all_perms), False)
        subgroup_map = np.zeros((len(all_perms),), dtype=int)

        (groups, group_inds), _ = nput.group_by(np.arange(len(eval_cinds)), eval_cinds)
        # mask = np.full(len(eval_perms), False)
        for cind_group,val_group in zip(groups, group_inds):
            # prefacs_full = prefactors[val_group,]
            pi = eval_perms[val_group,]
            g_key = cind_specs[cind_group]

            subgroup_map[pi] = val_group
            all_perms_mask[pi] = True
            split_ap_pos = np.split(all_perms_mask[perm_map], split_spec)
            split_ev_pos = np.split(subgroup_map[perm_map], split_spec)
            split_ap_inds = np.split(perm_map, split_spec)
            for state_idx, ((state, aperms), a_mask, ev_pos, a_inds) in enumerate(
                    zip(state_perms, split_ap_pos, split_ev_pos, split_ap_inds)
            ):
                # `up_pos` tells us which elements of `all_perms` correspond to `aperms`
                # `eval_perms` tells us which elements of `all_perms` to include at all
                # `val_groups` tells us which elements of `eval_perms` to include for this block
                # so in total we need to effectively take the intersection of `all_perms` and `eval_perms[val_groups]`
                # to find the `aperms` to take, and for each of those we need to track not just
                # what values of `eval_perms[val_groups]` are `True`, but also which `val_groups` they corresponded to

                mask_pos = np.where(a_mask)
                if len(mask_pos) == 0 or len(mask_pos[0]) == 0: continue
                mask_pos = mask_pos[0]

                state = state[np.newaxis] # old implementation expected state blocks
                vps = ev_pos[mask_pos,]
                prefacs = prefactors[vps,].T
                which_perms = a_inds[mask_pos]
                perm_subsets = aperms[mask_pos,][:, full_set]
                perm_substates, tuple_states, perm_freqs, state_block_key = cls._get_state_perms(
                    state_idx, state, freqs, fixed, subset, perm_subsets, take_cache, tuple(mask_pos), full_set
                )


                # subpowers = np.power(substates[np.newaxis, :, :], np.arange(max_order)[:, np.newaxis, np.newaxis])
                # perm_powers.append(np.moveaxis(subpowers, -1, 0))

                if degenerate_changes is not None:
                    perm_subsets = [tuple(p) for p in perm_subsets]

                route_context = (evaluation_context.for_route(state_idx)
                                 if hasattr(evaluation_context, 'for_route') else evaluation_context)
                if verbose:
                    with logger.block(
                            tag="{k} ({p})",
                            preformatter=lambda **vars: dict(vars, k=PTTensorCoeffProductSum.format_key(vars['k'])),
                            k=g_key,
                            p=prefacs,
                            log_level=log_level
                    ):
                        subcontrib = cls._eval_perm_core(
                            expr, state, tuple_states, perm_substates, which_perms,
                            change, baseline_shift,
                            prefacs, perm_freqs, state_block_key,
                            pows, g_key, perm_subsets, degenerate_changes, only_degenerate_terms,
                            poly_cache, energy_cache,
                            verbose, logger, log_level, log_scaling,
                            zero_cutoff, route_context
                        )
                        val = prefacs[:, np.newaxis, :] * subcontrib[np.newaxis]
                        logger.log_print("contrib: {e}", e=val,
                                         preformatter=lambda **vars: dict(vars, e=vars['e'].squeeze()*log_scaling),
                                         log_level=log_level)
                        contrib[state_idx][:, a_mask] += val[:, 0, :]

                    logger.log_print("{c}",
                                     c=val,
                                     preformatter=lambda **vars: dict(vars, c=vars['c'][:, 0, :].squeeze() * log_scaling),
                                     log_level=log_level)
                else:
                    subcontrib = cls._eval_perm_core(
                        expr, state, tuple_states, perm_substates, which_perms,
                        change, baseline_shift,
                        prefacs, perm_freqs, state_block_key,
                        pows, g_key, perm_subsets, degenerate_changes, only_degenerate_terms,
                        poly_cache, energy_cache,
                        verbose, logger, log_level, log_scaling,
                        zero_cutoff, route_context
                    )
                    val = prefacs[:, np.newaxis, :] * subcontrib[np.newaxis]
                    contrib[state_idx][:, mask_pos] += val[:, 0, :]

            # subgroup_map[pi] = 0
            all_perms_mask[pi] = False
        return contrib

    @classmethod
    def _prep_coeff_data(cls, cind_sets, coeff_lists):
        final_tensors = []
        final_cinds = []
        bound = {}
        def coefficient(coeffs, ci):
            order,kind = ci[:2]
            if not isinstance(kind,_HamiltonianCoefficient):
                return coeffs[order][kind] if len(coeffs)>order else 0
            identity = id(coeffs),order,kind
            if identity not in bound:
                row = coeffs[order] if len(coeffs)>order else ()
                sources = tuple(row[i] if i<len(row) else 0 for i in (0,1))
                bound[identity] = (_BoundHamiltonianCoefficient(kind,sources)
                                   if any(np.ndim(source) for source in sources) else 0)
            return bound[identity]
        for cinds in cind_sets:
            bad = True
            subtensors = []
            for coeffs in coeff_lists:
                try:
                    tensors = [
                        coefficient(coeffs,ci)
                        for ci in cinds
                    ]
                except:
                    print([
                        [x.shape if not nput.is_zero(x) else 0 for x in c]
                        for c in coeffs
                    ], cinds)
                    raise
                if bad and not any(nput.is_zero(t) for t in tensors): bad = False
                # indices = [ (ci[:2], ci[:2]) for ci in cinds ]
                subtensors.append(tensors)
            if not bad:
                final_tensors.append(subtensors)
                final_cinds.append(cinds)
        return final_tensors, final_cinds

    @classmethod
    def _get_max_order(cls, expr):
        if isinstance(expr, PolyPath):
            return max(expr.order + (0,)) + 1
        elif isinstance(expr, ProductPTPolynomial):
            return 1 if len(expr.coeffs) == 0 else max(len(c) for c in expr.coeffs)
        elif isinstance(expr, ProductPTPolynomialSum):
            return max(cls._get_max_order(p) for p in expr.polys)
        else:
            return max(
                [-1] +
                [
                    cls._get_max_order(e)
                    for e in expr.terms.values()
                ]
            )

    @classmethod
    def _deg_test_direct(cls, state, modes, check_pos):
        # I thought about searchsorted, but for this
        # really short state-by-state case that's probably
        # slower
        return all(
            state[c] == modes[c]
            for c in check_pos
        )
    _direct_degs = collections.namedtuple('_direct_degs', ['check_pos', 'modes'])
    @classmethod
    def _deg_test(cls, modes):
        check_pos = tuple(i for i,m in enumerate(modes) if m >= 0)
        return cls._direct_degs(check_pos, modes)

    @classmethod
    def _deg_test_tree(cls, state, tree):
        # a small amount of recursion for simplicity
        for check_pos,match_vals in tree.items(): # test all possible `check_pos`
            mode = state[check_pos]
            res = match_vals.get(mode)
            if isinstance(res, dict):
                res = cls._deg_test_tree(state, res)
            if res: return True

    @staticmethod
    def _mode_inclusion_direct(modes, tests):
        return modes in tests

    @staticmethod
    def _mode_inclusion_multi(modes, tests):
        return any(t(modes) for t in tests)
    @classmethod
    def _compile_deg_tests(cls, tests):
        if all(isinstance(t, cls._direct_degs) for t in tests):
            tree = {}
            # becomes a tree like:
            # {
            #   0: { # position in the state to check
            #        6: {  # value to match
            #              1: { # position in the state to check
            #                    4: {
            #                            2: {3:True, 5:True}
            #                       }
            #                 }
            #           },
            #        5: {  # value to match
            #              1: { # position in the state to check
            #                    2: {}
            #                 }
            #           }
            # }
            for check_pos, modes in tests:
                subtree = tree
                for n,c in enumerate(check_pos):
                    m = modes[c]
                    if c not in subtree: subtree[c] = {}
                    if n == len(check_pos) - 1: # might be better to do this with an `else`
                        subtree[c][m] = True
                    else:
                        if m not in subtree[c]: subtree[c][m] = {}
                        subtree = subtree[c][m]
                        if subtree is True: # preempted
                            break
            # import pprint
            # print(tests)
            # pprint.pprint(tree)
            test = functools.partial(cls._deg_test_tree, tree=tree)
        else:
            test = [
                t
                    if not isinstance(t, cls._direct_degs) else
                functools.partial(cls._deg_test_direct, check_pos=t.check_pos, modes=t.modes)
                for t in tests
            ]
        return test
    @classmethod
    def _make_full_deg_test(cls, tests, array_oriented=True):
        if not array_oriented:
            if isinstance(tests, set):
                return functools.partial(cls._mode_inclusion_direct, tests=tests)
            test = cls._compile_deg_tests(tests)
            if isinstance(test, list):
                test = functools.partial(cls._mode_inclusion_multi, tests=test)
            return test
        patterns = []
        checks = []
        fallback_tests = []
        for test in tests:
            if isinstance(test, cls._direct_degs):
                pattern = np.asarray(test.modes, dtype=np.intp)
                check = np.zeros(len(pattern), dtype=bool)
                check[np.asarray(test.check_pos, dtype=np.intp)] = True
                patterns.append(pattern)
                checks.append(check)
            elif isinstance(test, np.ndarray):
                pattern = np.asarray(test, dtype=np.intp)
                patterns.append(pattern)
                checks.append(pattern >= 0)
            elif isinstance(test, tuple) and all(
                    isinstance(value, (int, np.integer)) for value in test):
                pattern = np.asarray(test, dtype=np.intp)
                patterns.append(pattern)
                checks.append(np.ones(len(pattern), dtype=bool))
            else:
                fallback_tests.append(test)
        if len(patterns) == 0:
            test = cls._compile_deg_tests(fallback_tests)
            if isinstance(test, list):
                test = functools.partial(cls._mode_inclusion_multi, tests=test)
            return test
        return DegeneracyTestPlan(
            np.stack(patterns),
            np.stack(checks),
            fallback_tests=fallback_tests
        )

    default_deg_id_method = 'linear'
    deg_id_method_nmodes_switch = 50
    use_compiled_degeneracy_plans = True
    use_batched_degeneracy_pattern_construction = True
    use_compiled_degeneracy_identification = True
    @classmethod
    def _identify_possible_degeneracies_legacy(cls, all_perms, sides, expr, changes, nmodes, method=None):
        # changes is stored as `(modes, quanta)` pairs
        # representing a change between degenerate modes
        if method is None:
            method = cls.default_deg_id_method
        if method is None:
            if nmodes > cls.deg_id_method_nmodes_switch:
                method = 'linear'
            else:
                method = 'perfect'

        left_changes = []
        right_changes = []
        # initial_changes = []
        for modes,quanta in changes:
            quanta = np.asanyarray(quanta)
            sorting = np.argsort(quanta)
            modes = tuple(modes[k] for k in sorting)
            sort_quants = quanta[sorting]
            # for side in sides:
            left_changes.append([modes, 0, sort_quants])
            right_changes.append([modes, 1, sort_quants])

        all_changes = [left_changes, right_changes]
        all_degs = [{}, {}]

            # initial_changes.append([tuple(reversed(modes)), 1, np.flip(-sort_quants)])
        # the final state has similar resonance conditions, but we need to
        # adjust for the shift in the baseline state, which means we need
        # to iterate over perms
        # all_changes = [
        #     initial_changes.copy()
        #     for _ in all_perms
        # ]
        # num_changes = len(total_change)
        # for change_list, perm in zip(all_changes, all_perms):
        #     changed_pos = perm[:num_changes]
        #     for modes, quanta in changes:
        #         modes = list(modes)
        #         quanta = list(quanta) # force a copy
        #         for p,q in zip(changed_pos, total_change):
        #             try:
        #                 q_idx = modes.index(p)
        #             except ValueError:
        #                 modes.append(p)
        #                 quanta.append(q)
        #             else:
        #                 quanta[q_idx] += q
        #         quanta = np.asanyarray(quanta)
        #         sorting = np.argsort(quanta)
        #         modes = tuple(modes[k] for k in sorting)
        #         sort_quants = quanta[sorting]
        #         change_list.append([modes, 1, sort_quants])

        # all_degs = [
        #     {} for _ in all_perms
        # ]
        deg_test_cache = {}
        for cinds,t in expr.terms.items():
            if isinstance(t, PTEnergyChangeProductSum):
                # echange_tests = {}
                for ekey in t.terms.keys():
                    for eside,echange in t.side_change_iter(ekey):
                        # if known_deg: break
                        zero_pos = [
                            i for i,e in enumerate(echange)
                            if e == 0
                        ]
                        nonzero_pos = [
                            i for i, e in enumerate(echange)
                            if e != 0
                        ]
                        sub_e = [echange[i] for i in nonzero_pos]
                        echange_sorting = np.argsort(sub_e)
                        num_zero = len(zero_pos)
                        echange_inverse = np.argsort(
                                [nonzero_pos[i] for i in echange_sorting] + zero_pos
                        )
                        # need to permute over sorted echange blocks
                        echange_splits = np.nonzero(np.diff([sub_e[i] for i in echange_sorting]))
                        if len(echange_splits) > 0 and len(echange_splits[0]) > 0:
                            echange_blocks = [
                                list(itertools.permutations(p))
                                for p in np.split(np.arange(len(sub_e)), echange_splits[0] + 1)
                            ]
                        else:
                            echange_blocks = [
                                [np.arange(len(sub_e))]
                            ]
                        subpermutations = [
                            np.concatenate(parts).astype(np.intp, copy=False)
                            if len(parts) > 0 else
                            np.empty(0, dtype=np.intp)
                            for parts in itertools.product(*echange_blocks)
                        ]
                        echange_inverse = np.asarray(echange_inverse, dtype=np.intp)
                        if (
                                method == 'linear'
                                and cls.use_compiled_degeneracy_plans
                        ):
                            # A negative source is a wildcard corresponding to
                            # a zero-energy position.  These maps depend only
                            # on the energy-change term and are compiled once,
                            # outside the resonance-mode loop.
                            mode_maps = []
                            for subperm in subpermutations:
                                sources = np.concatenate((
                                    subperm,
                                    np.full(num_zero, -1, dtype=np.intp)
                                ))
                                mode_maps.append(sources[echange_inverse])
                            mode_maps = np.stack(mode_maps)
                        elif method == 'perfect':
                            mode_maps = []
                            pad_sources = np.arange(
                                len(sub_e),
                                len(sub_e) + num_zero,
                                dtype=np.intp
                            )
                            for subperm in subpermutations:
                                sources = np.concatenate((
                                    subperm,
                                    pad_sources
                                ))
                                mode_maps.append(sources[echange_inverse])
                        else:
                            mode_maps = None

                        for degs,sort_changes in zip(all_degs,all_changes):
                            for modes,side,quanta in sort_changes: # strict ordering preserved
                                if (
                                        eside == side
                                        and len(sub_e) == len(quanta)
                                        and all(sub_e[s] == q for s, q in zip(echange_sorting, quanta))
                                ):
                                    # modes = tuple(modes)
                                    if cinds not in degs: degs[cinds] = {}

                                    if method == 'perfect':
                                        if ekey not in degs[cinds]: degs[cinds][ekey] = []
                                        available_modes = [
                                            mode for mode in range(nmodes)
                                            if mode not in modes
                                        ]
                                        pad_sets = (
                                            itertools.permutations(
                                                available_modes,
                                                num_zero
                                            )
                                            if num_zero > 0 else
                                            [()]
                                        )
                                        mode_array = np.asarray(modes, dtype=np.intp)
                                        for pad_inds in pad_sets:
                                            base_modes = np.concatenate((
                                                mode_array,
                                                np.asarray(pad_inds, dtype=np.intp)
                                            ))
                                            for mode_map in mode_maps:
                                                degs[cinds][ekey].append(
                                                    base_modes[mode_map]
                                                )
                                    elif method == 'linear':
                                        if ekey not in degs[cinds]: degs[cinds][ekey] = []
                                        if (
                                                cls.use_compiled_degeneracy_plans
                                                and cls.use_batched_degeneracy_pattern_construction
                                        ):
                                            mode_array = np.asarray(modes, dtype=np.intp)
                                            patterns = np.full(
                                                mode_maps.shape,
                                                -1,
                                                dtype=np.intp
                                            )
                                            checked = mode_maps >= 0
                                            patterns[checked] = mode_array[
                                                mode_maps[checked]
                                            ]
                                            # Retain one dense chunk per
                                            # resonance match instead of one
                                            # Python/NumPy object per pattern.
                                            degs[cinds][ekey].append(patterns)
                                        elif cls.use_compiled_degeneracy_plans:
                                            mode_array = np.asarray(modes, dtype=np.intp)
                                            for mode_map in mode_maps:
                                                pattern = np.full(
                                                    len(mode_map),
                                                    -1,
                                                    dtype=np.intp
                                                )
                                                checked = mode_map >= 0
                                                pattern[checked] = mode_array[
                                                    mode_map[checked]
                                                ]
                                                degs[cinds][ekey].append(pattern)
                                        else:
                                            for subperm in subpermutations:
                                                pad_inds = (-1,) * num_zero
                                                full_mode = modes + pad_inds
                                                full_mode = tuple(
                                                    full_mode[index]
                                                    for index in subperm
                                                ) + pad_inds
                                                perm_modes = tuple(
                                                    full_mode[index]
                                                    for index in echange_inverse
                                                )
                                                degs[cinds][ekey].append(
                                                    cls._deg_test(perm_modes)
                                                )
                                    else:
                                        if ekey not in degs[cinds]: degs[cinds][ekey] = []
                                        degs[cinds][ekey].append(method(modes, echange)) # TODO: supply more info to custom methods
                    for degs in all_degs:
                        if len(degs.get(cinds, {}).get(ekey, [])) > 0:
                            raw_tests = degs[cinds][ekey]
                            array_oriented = (
                                cls.use_compiled_degeneracy_plans
                                and method in {'linear', 'perfect'}
                            )
                            if (
                                    array_oriented
                                    and all(
                                        isinstance(test, np.ndarray)
                                        for test in raw_tests
                                    )
                            ):
                                raw_patterns = (
                                    np.concatenate(raw_tests, axis=0)
                                    if raw_tests[0].ndim == 2 else
                                    np.stack(raw_tests)
                                )
                                dc_key = (
                                    'array',
                                    raw_patterns.shape,
                                    raw_patterns.tobytes()
                                )
                            else:
                                raw_patterns = None
                                dc_key = ('callable', tuple(raw_tests))
                            if dc_key not in deg_test_cache:
                                deg_test_cache[dc_key] = (
                                    DegeneracyTestPlan(raw_patterns)
                                    if raw_patterns is not None else
                                    cls._make_full_deg_test(
                                        raw_tests,
                                        array_oriented=array_oriented
                                    )
                                )
                            degs[cinds][ekey] = deg_test_cache[dc_key]
        return all_degs

    @classmethod
    def _identify_possible_degeneracies(cls, all_perms, sides, expr, changes,
                                        nmodes, method=None, context=None):
        if isinstance(changes, SideDegeneracyChanges):
            left = (
                cls._identify_possible_degeneracies(
                    all_perms, sides, expr, changes.left, nmodes,
                    method=method
                )[0]
                if changes.left else {}
            )
            right = (
                cls._identify_possible_degeneracies(
                    all_perms, sides, expr, changes.right, nmodes,
                    method=method
                )[1]
                if changes.right else {}
            )
            return [left, right]
        if method is None:
            method = cls.default_deg_id_method
        if method is None:
            method = (
                'linear' if nmodes > cls.deg_id_method_nmodes_switch
                else 'perfect'
            )
        if (
                not cls.use_compiled_degeneracy_identification
                or not cls.use_compiled_degeneracy_plans
                or method not in {'linear', 'perfect'}
        ):
            result = cls._identify_possible_degeneracies_legacy(
                all_perms, sides, expr, changes, nmodes, method=method
            )
        else:
            if context is None:
                context = DegeneracyIdentificationContext(changes)
            result = context.bind(expr, nmodes, method)
        return result


    @classmethod
    def set_cache_size(cls, max_size):
        if cls._max_cache_size != max_size:
            cls._max_cache_size = max_size
            cls._poly_cache = MaxSizeCache(cls._max_cache_size, cache_type='fifo')  # temporary hack, but these are in principle shared/reused
            cls._ecoeff_cache = MaxSizeCache(cls._max_cache_size, cache_type='fifo')
    @classmethod
    def get_cache(self):
        return MaxSizeCache(max_items=self._max_cache_size, cache_type='fifo')

    @classmethod
    def use_persistent_cache(cls, store):
        """
        Puts an on-disk layer (an `ObjectStore` or a path to one; there is no default location) under the exact-result
        polynomial cache so complete polynomial/state values survive restarts; call `flush_persistent_cache`
        (e.g. at checkpoints) to write new values. `clear_caches` and `set_cache_size` drop the layer.
        Pass ``False`` to remove it.
        """
        from .Evaluators import PersistentEvaluationCache
        if store is False:
            current = cls._poly_cache
            cls._poly_cache = current.memory if isinstance(current, PersistentEvaluationCache) else current
            return None
        if store is None or store is True:
            raise ValueError("pass an ObjectStore or a path to one (or False to remove the layer)")
        current = cls._poly_cache
        if isinstance(current, PersistentEvaluationCache):
            current = current.memory
        cls._poly_cache = PersistentEvaluationCache(store, memory=current)
        return cls._poly_cache

    @classmethod
    def flush_persistent_cache(cls):
        """Writes new values of the persistent polynomial cache (a no-op if there is none)"""
        flush = getattr(cls._poly_cache, 'flush', None)
        return 0 if flush is None else flush()

    _max_cache_size = 1e7
    _poly_cache = MaxSizeCache(_max_cache_size, cache_type='fifo') # temporary hack, but these are in principle shared/reused
    _ecoeff_cache = MaxSizeCache(_max_cache_size, cache_type='fifo')
    use_balanced_chunk_sums = True
    default_tensor_evaluation_mode = 'dag'
    default_evaluation_backend = 'indexed'
    default_dag_cache_size = 100000
    default_dag_cache_bytes = 64 * 1024 ** 2
    default_dag_chunk_size = 256
    default_dag_node_term_cache_size = 256
    default_state_permutation_cache_size = 100000
    default_state_permutation_cache_bytes = 64 * 1024 ** 2
    use_hashable_numpy_state_indices = True
    use_hashable_numpy_mask_indices = True
    # Direct construction benchmarks cross over near eight int64 elements:
    # tuples are faster below this size, while compact byte identities become
    # both faster and smaller above it.
    hashable_numpy_index_threshold = 8
    # A one-state path entry is about 180 bytes in the current representation,
    # so 350k entries lets the byte limit, rather than an unnecessarily small
    # item limit, govern the default 64 MiB cache.  Larger state vectors reach
    # the byte limit with proportionally fewer entries.
    default_materialized_path_cache_size = 350000
    default_materialized_path_cache_bytes = 64 * 1024 ** 2
    _last_dag_evaluation_stats = None
    _last_materialized_evaluation_stats = None
    _last_state_permutation_cache_stats = None
    _last_indexed_evaluation_stats = None
    _indexed_evaluation_generation = 0
    default_zero_cutoff = 1e-24

    @classmethod
    def get_last_dag_evaluation_stats(cls):
        return None if cls._last_dag_evaluation_stats is None else dict(cls._last_dag_evaluation_stats)

    @classmethod
    def get_last_materialized_evaluation_stats(cls):
        return (
            None
            if cls._last_materialized_evaluation_stats is None else
            dict(cls._last_materialized_evaluation_stats)
        )

    @classmethod
    def get_last_indexed_evaluation_stats(cls):
        return (
            None
            if cls._last_indexed_evaluation_stats is None else
            dict(cls._last_indexed_evaluation_stats)
        )

    @staticmethod
    def _compare_evaluation_backends(legacy, indexed, rtol=2e-12, atol=2e-12):
        if len(legacy) != len(indexed):
            raise AssertionError(
                "indexed backend returned {} result blocks; legacy returned {}"
                .format(len(indexed), len(legacy))
            )
        for block_index, (legacy_block, indexed_block) in enumerate(zip(legacy, indexed)):
            if np.shape(legacy_block) != np.shape(indexed_block):
                raise AssertionError(
                    "indexed result block {} has shape {}; expected {}"
                    .format(block_index, np.shape(indexed_block), np.shape(legacy_block))
                )
            np.testing.assert_allclose(
                indexed_block,
                legacy_block,
                rtol=rtol,
                atol=atol,
                equal_nan=True,
                err_msg="indexed evaluator mismatch in result block {}".format(block_index)
            )

    @classmethod
    def _new_state_permutation_cache(cls, max_items=None, max_bytes=None):
        return BoundedCache(
            cls.default_state_permutation_cache_size
            if max_items is None else max_items,
            cls.default_state_permutation_cache_bytes
            if max_bytes is None else max_bytes
        )

    _parallel_eval_main_args = None
    _cached_expansion = None
    @classmethod
    def _initialize_main_eval(cls, coeffs, *args, parallelizer=None):
        if coeffs is not None:
            cls._cached_expansion = coeffs
        cls._parallel_eval_main_args = args
        cls._cached_main_args = None
        # parallelizer.print("wat {l}", l=len(args))

    _cached_main_args = None
    _cached_stream_preparation = None
    @classmethod
    def _extract_main_args(cls):
        if cls._cached_main_args is not None:
            return cls._cached_main_args + (
                cls._cached_expansion,
                cls._parallel_eval_main_args[-1]
            )

        if cls._parallel_eval_main_args is None:
            raise ValueError("parallel eval not initialized")

        (
            expr, change, baseline_shift,
            # state_perms,  # all_perms, perm_map,
            compressed_states, perm_lens, compressed_perms,
            freqs,
            num_fixed, degenerate_changes, only_degenerate_terms,
            zero_cutoff,
            max_cache_size,
            use_materialized_path_cache,
            path_cache_size, path_cache_bytes,
            state_cache_size, state_cache_bytes,
            # counts_cache, poly_cache, take_cache, energy_cache,
            # pows,
            verbose, logger, log_scaled, log_level,
            stream_preparation_entries, ndim, evaluation_backend
        ) = cls._parallel_eval_main_args

        cls._cached_stream_preparation = None
        if stream_preparation_entries is not None:
            from .Evaluators import StreamPreparationContext
            cls._cached_stream_preparation = StreamPreparationContext(
                cls._cached_expansion, stream_preparation_entries)

        cls.set_cache_size(max_cache_size)

        decompressed_states = np.zeros((len(compressed_states), ndim), dtype=int)
        for n,(elems,pos) in enumerate(compressed_states):
            decompressed_states[n, pos] = elems
        decompressed_perms = mcomb.lehmer_decode(ndim, compressed_perms)
        split_perms = np.array_split(
            decompressed_perms,
            np.cumsum(perm_lens)
        )
        state_perms = [
            [s,p]
            for s,p in zip(decompressed_states, split_perms)
        ]

        all_perms, perm_map = np.unique(
            np.concatenate([perms for state, perms in state_perms], axis=0),
            axis=0, return_inverse=True
        )

        max_order = cls._get_max_order(expr)
        if max_order < 0:
            raise Exception(expr)
        max_state = np.max(np.concatenate([state for state, perms in state_perms]))
        pows = np.power(np.arange(max_state + 1)[np.newaxis, :], np.arange(max_order + 1)[:, np.newaxis])
        counts_cache = cls.get_cache()  # so we don't hit some coeff sets too many times
        poly_cache = (
            BoundedCache.sampled(
                path_cache_size or cls.default_materialized_path_cache_size,
                path_cache_bytes if path_cache_bytes is not None else cls.default_materialized_path_cache_bytes)
            if stream_preparation_entries is not None else
            _MaterializedEvaluationCache(
                cls._poly_cache,
                max_items=path_cache_size,
                max_bytes=path_cache_bytes
            )
            if use_materialized_path_cache else
            cls._poly_cache
        )
        take_cache = cls._new_state_permutation_cache(
            state_cache_size, state_cache_bytes
        )
        energy_cache = cls.get_cache()

        cls._cached_main_args = (
            expr, change, baseline_shift,
            state_perms, all_perms, perm_map,
            freqs,
            num_fixed, degenerate_changes, only_degenerate_terms,
            zero_cutoff,
            counts_cache, poly_cache, take_cache, energy_cache,
            pows,
            verbose, logger, log_scaled, log_level,
            ndim
        )

        return cls._cached_main_args + (cls._cached_expansion, evaluation_backend)

    @staticmethod
    def _aggregate_indexed_evaluation_stats(stats, elapsed_mode='sum'):
        """Combine independent indexed executors without hiding peak usage.

        Worker-local counters describe disjoint work and are additive.  Peak
        storage counters describe concurrent process-local bounds and retain
        their maximum; they must not be reported as if one worker allocated
        the sum.  Executors within one process group overlap in wall time,
        whereas successive free-index groups do not.
        """
        stats = [entry for entry in stats if entry is not None]
        if len(stats) == 0:
            return None
        aggregate = {}
        for key in stats[0]:
            values = [entry[key] for entry in stats]
            if key.startswith('max_'):
                aggregate[key] = max(values)
            elif key.endswith('_seconds') and elapsed_mode == 'max':
                aggregate[key] = max(values)
            else:
                aggregate[key] = sum(values)
        return aggregate
    @classmethod
    def _partition_evaluation_blocks(cls, ncomb, nproc):
        """Split combination ranks into complete, balanced half-open ranges."""
        if nproc < 1:
            raise ValueError("evaluation requires at least one process")
        block_size, remainder = divmod(ncomb, nproc)
        blocks = []
        start = 0
        for block_index in range(nproc):
            stop = start + block_size + (block_index < remainder)
            if stop > start:
                blocks.append((start, stop))
            start = stop
        return blocks

    @classmethod
    def _run_eval_block(cls,
                        block_ranges,
                        contrib_shapes,
                        free_inds,
                        cind_sets, mode_subset_selector=None
                        ):
        (
            expr, change, baseline_shift,
            state_perms, all_perms, perm_map,
            freqs,
            num_fixed, degenerate_changes, only_degenerate_terms,
            zero_cutoff,
            counts_cache, poly_cache, take_cache, energy_cache,
            pows,
            verbose, logger, log_scaled, log_level,
            ndim, coeffs, evaluation_backend
        ) = cls._extract_main_args()

        indexed_executor = None
        eval_perm = cls._eval_perm
        if evaluation_backend == 'indexed':
            from .Evaluators import IndexedBlockEvaluator
            indexed_executor = IndexedBlockEvaluator(cls, expr, preparation_context=cls._cached_stream_preparation)
            eval_perm = indexed_executor.eval_perm

        tensors, cind_specs = (cls._prep_coeff_data(cind_sets, coeffs)
            if cls._cached_stream_preparation is None else
            cls._cached_stream_preparation.prepare(cls, cind_sets))

        facs = [
            math.factorial(n)
            for n in range(max([len(cinds) for cinds in cind_specs] + [0]) + 1)
        ]
        split_spec = np.cumsum([0] + [len(perms) for state, perms in state_perms])[1:]

        contrib = [np.zeros(c) for c in contrib_shapes]
        def combinations():
            return (itertools.combinations(range(num_fixed, ndim), r=free_inds)
                    if mode_subset_selector is None else
                    mode_subset_selector.iter_subsets(ndim,num_fixed,free_inds))
        comb_iter = combinations()
        combination_position = 0
        for start, stop in sorted(block_ranges):
            if start < combination_position:
                # Keep this correct if a future scheduler assigns overlapping
                # or out-of-order ranges, without penalizing today's ordered
                # contiguous partition.
                comb_iter = combinations()
                combination_position = 0
            for subset in itertools.islice(
                    comb_iter,
                    start - combination_position,
                    stop - combination_position
            ):
                # print(n, math.comb(ndim-num_fixed, free_inds), free_inds, start, stop, subset)
                if verbose:
                    with logger.block(tag="{b} + {s}", b=tuple(range(num_fixed)), s=subset,
                                      log_level=log_level):
                        subcontrib = eval_perm(
                            expr, change, baseline_shift,
                            subset, state_perms, all_perms, perm_map,
                            freqs, cind_specs, tensors,
                            num_fixed, degenerate_changes, only_degenerate_terms,
                            zero_cutoff,
                            counts_cache, poly_cache, take_cache, energy_cache,
                            facs, pows, split_spec,
                            verbose, logger, log_scaled,
                            evaluation_context=expr
                        )
                else:
                    subcontrib = eval_perm(
                        expr, change, baseline_shift,
                        subset, state_perms, all_perms, perm_map,
                        freqs, cind_specs, tensors,
                        num_fixed, degenerate_changes, only_degenerate_terms,
                        zero_cutoff,
                        counts_cache, poly_cache, take_cache, energy_cache,
                        facs, pows, split_spec,
                        verbose, logger, log_scaled,
                        evaluation_context=expr
                    )
                for storage, corr in zip(contrib, subcontrib):
                    storage += corr
            combination_position = stop
        return contrib, (
            None if indexed_executor is None else indexed_executor.finish()
        )

    @classmethod
    def _run_eval_blocks(cls,
                         blocks,
                         contrib_shapes,
                         *args, parallelizer:Parallelizer):
        if parallelizer.on_main:
            contrib = [np.zeros(c) for c in contrib_shapes]
        else:
            contrib = None
        map_res = parallelizer.map(
            cls._run_eval_block,
            blocks,
            extra_args=(contrib_shapes,) + args,
            vectorized=True,
            aggregate=False
        )
        # print(map_res)
        if parallelizer.on_main:
            # print("!!!", np.array(subcontrib).shape, np.array(contrib).shape)
            indexed_stats = []
            for subcontrib, worker_stats in map_res:
                for storage, corr in zip(contrib, subcontrib):
                    storage += corr
                if worker_stats is not None:
                    indexed_stats.append(worker_stats)
            return (
                contrib,
                cls._aggregate_indexed_evaluation_stats(
                    indexed_stats, elapsed_mode='max'
                )
            )
        return contrib, None

    @classmethod
    def _evaluate_tensor_dag_expression(cls,
                                        state_perms, coeffs, freqs,
                                        expr, change, baseline_shift,
                                        num_fixed,
                                        op=None,
                                        logger=None,
                                        parallelizer=None,
                                        degenerate_changes=None,
                                        only_degenerate_terms=False,
                                        zero_cutoff=None,
                                        verbose=False, log_scaled=True,
                                        cache_size=None,
                                        cache_bytes=None,
                                        chunk_size=None,
                                        state_cache_size=None,
                                        state_cache_bytes=None,
                                        use_path_cache=True,
                                        evaluation_backend=None,
                                        degeneracy_identification_context=None,
                                        mode_subset_selector=None,
                                        stream_options=None):
        cache_size = cls.default_dag_cache_size if cache_size is None else cache_size
        cache_bytes = cls.default_dag_cache_bytes if cache_bytes is None else cache_bytes
        chunk_size = cls.default_dag_chunk_size if chunk_size is None else max(1, int(chunk_size))
        budget = preparation = None
        if stream_options is not None:
            from .Evaluators import StreamBudget, StreamPreparationContext
            budget = StreamBudget(stream_options)
            preparation = StreamPreparationContext(
                coeffs, budget.preparation_entries if budget.reuse_preparation else 0)
        evaluation_cache = BoundedCache.sampled(cache_size, cache_bytes)
        state_permutation_cache = cls._new_state_permutation_cache(
            state_cache_size, state_cache_bytes
        )
        # Keep tensor-plan entries and numerical PolyPath intermediates under
        # one batch-local bound, while retaining the legacy complete-result
        # cache as the first lookup level.  This gives DAG chunks the same
        # shared-node reuse as materialized evaluation without doubling the
        # configured memory ceiling or retaining intermediates across batches.
        polynomial_cache = (
            _MaterializedEvaluationCache(
                cls._poly_cache,
                batch_cache=evaluation_cache
            )
            if use_path_cache and budget is None else
            BoundedCache.sampled(cache_size, cache_bytes) if budget is not None else evaluation_cache
        )
        plan = PTTensorCoeffProductDAGEvaluationPlan(
            expr,
            cache=evaluation_cache,
            max_cached_terms=cls.default_dag_node_term_cache_size,
            max_mode_indices=state_perms[0][0].shape[-1]
        )

        result = tuple(
            np.zeros([len(coeffs), len(perms)])
            for state, perms in state_perms
        )
        pending = {}
        balanced_sums = cls.use_balanced_chunk_sums
        maximum_sum_levels = 0
        chunks = 0
        generation_seconds = accumulation_seconds = evaluation_seconds = 0.
        numerical_sums = [_BalancedPolynomialAccumulator() for _ in result]
        indexed_chunk_stats = None
        collect_stats = Diagnostics.enabled
        materializations_before = (PTTensorCoeffProductDAG.cache_info()['tensor_materializations']
                                   if collect_stats else 0)
        start = time.perf_counter()

        def evaluate_pending(reason='end'):
            nonlocal chunks, evaluation_seconds, indexed_chunk_stats, polynomial_cache
            if len(pending) == 0:
                return
            if budget is not None:
                budget.event('evaluate_start', reason, len(pending))
            tick = time.perf_counter()
            if len(pending) == 0:
                return
            chunk = PTTensorCoeffProductSum(
                {key: value.total() for key, value in pending.items()}
                    if balanced_sums else dict(pending),
                canonicalize=False, reduced=False
            )
            chunk_result = cls.evaluate_polynomial_expression(
                state_perms, coeffs, freqs,
                chunk, change, baseline_shift,
                num_fixed,
                op=op,
                logger=logger,
                parallelizer=parallelizer,
                degenerate_changes=degenerate_changes,
                only_degenerate_terms=only_degenerate_terms,
                zero_cutoff=zero_cutoff,
                verbose=verbose, log_scaled=log_scaled,
                evaluation_mode='materialized',
                evaluation_backend=evaluation_backend,
                _poly_eval_cache=polynomial_cache,
                _state_permutation_cache=state_permutation_cache,
                _degeneracy_identification_context=(
                    degeneracy_identification_context
                ),
                mode_subset_selector=mode_subset_selector,
                _stream_preparation=preparation
            )
            if collect_stats and evaluation_backend == 'indexed':
                chunk_stats = cls.get_last_indexed_evaluation_stats()
                if chunk_stats is not None:
                    indexed_chunk_stats = cls._aggregate_indexed_evaluation_stats(
                        [indexed_chunk_stats, chunk_stats])
            for index, (storage, contribution) in enumerate(zip(result, chunk_result)):
                if budget is None:
                    storage += contribution
                else:
                    numerical_sums[index].append(contribution)
            pending.clear()
            if budget is not None:
                # Numerical roots belong to this batch. Preparation and state
                # routes persist separately without keeping polynomial roots.
                polynomial_cache = BoundedCache.sampled(cache_size, cache_bytes)
                if budget.account_cache_graphs:
                    from .Evaluators import trim_to_retained_size
                    trim_to_retained_size(evaluation_cache)
            chunks += 1
            elapsed = time.perf_counter() - tick
            evaluation_seconds += elapsed
            if budget is not None:
                budget.event('evaluate_end', reason, 0, evaluation_seconds=elapsed)
                budget.flushed(reason)

        previous = time.perf_counter() if collect_stats else 0.
        for key, poly in plan.iter_terms():
            if collect_stats:
                tick = time.perf_counter()
                generation_seconds += tick - previous
            if balanced_sums:
                accumulator = pending.get(key)
                if accumulator is None:
                    accumulator = pending[key] = _BalancedPolynomialAccumulator()
                accumulator.append(poly)
                maximum_sum_levels = max(maximum_sum_levels, len(accumulator.blocks))
            else:
                pending[key] = pending.get(key, 0) + poly
            reason = budget.incoming(pending) if budget is not None else None
            if reason is None and len(pending) >= chunk_size:
                reason = 'keys'
            if collect_stats:
                accumulation_seconds += time.perf_counter() - tick
            if reason is not None:
                evaluate_pending(reason)
            if collect_stats:
                previous = time.perf_counter()
        evaluate_pending()
        if budget is not None:
            result = tuple(total.total() if total.blocks else zero
                           for total, zero in zip(numerical_sums, result))

        if not Diagnostics.enabled:
            return result
        state_cache_stats = state_permutation_cache.stats()
        cls._last_state_permutation_cache_stats = state_cache_stats
        cls._last_dag_evaluation_stats = {
            **plan.stats(),
            'dag_exact_cache_hits': getattr(polynomial_cache, 'exact_hits', 0),
            'dag_exact_cache_misses': getattr(polynomial_cache, 'exact_misses', 0),
            'dag_block_cache_hits': getattr(polynomial_cache, 'block_hits', 0),
            'dag_block_cache_misses': getattr(polynomial_cache, 'block_misses', 0),
            'dag_path_cache_enabled': use_path_cache,
            'dag_chunks': chunks,
            'dag_balanced_chunk_sums': balanced_sums,
            'dag_accumulation_maximum_levels': maximum_sum_levels,
            'dag_chunk_size': chunk_size,
            'dag_elapsed': time.perf_counter() - start,
            'dag_generation_seconds': generation_seconds,
            'dag_accumulation_seconds': accumulation_seconds,
            'dag_numerical_seconds': evaluation_seconds,
            **({} if budget is None else budget.stats()),
            **({} if preparation is None else preparation.stats()),
            'dag_materializations': (
                PTTensorCoeffProductDAG.cache_info()['tensor_materializations']
                - materializations_before
            ),
            **{
                'state_permutation_' + key: value
                for key, value in state_cache_stats.items()
            }
        }
        if indexed_chunk_stats is not None:
            cls._last_indexed_evaluation_stats = (
                indexed_chunk_stats
            )
            cls._indexed_evaluation_generation += 1
        return result

    @staticmethod
    def _normalize_coefficient_expansions(coeffs, preparation_context=None):
        """Canonicalize array-like physical tensors once for an immutable plan."""
        def coefficient(value):
            if nput.is_numeric(value):
                return value
            if isinstance(value,np.ndarray) or preparation_context is None:
                return np.asanyarray(value)
            root = preparation_context.root
            entry = root.normalized_inputs.get(id(value))
            if entry is None:
                entry = (value,np.asanyarray(value))
                root.normalized_inputs[id(value)] = entry
            return entry[1]
        return [[[coefficient(value) for value in row] for row in expansion]
                for expansion in coeffs]

    @classmethod
    def evaluate_polynomial_expression(cls,
                                       state_perms, coeffs, freqs,
                                       expr, change, baseline_shift,
                                       num_fixed,
                                       op=None,
                                       logger=None,
                                       parallelizer=None,
                                       degenerate_changes=None,
                                       only_degenerate_terms=False,
                                       zero_cutoff=None,
                                       verbose=False, log_scaled=True,
                                       evaluation_mode=None,
                                       evaluation_backend=None,
                                       dag_cache_size=None,
                                       dag_cache_bytes=None,
                                       dag_chunk_size=None,
                                       path_cache_size=None,
                                       path_cache_bytes=None,
                                       state_permutation_cache_size=None,
                                       state_permutation_cache_bytes=None,
                                       _poly_eval_cache=None,
                                       _state_permutation_cache=None,
                                       _degeneracy_identification_context=None,
                                       mode_subset_selector=None,
                                       dag_stream_options=None,
                                       _stream_preparation=None, projection_bindings=None,
                                       _max_polynomial_order=None
                                       ):

        if projection_bindings is not None:
            op = _ProjectionEvaluationContext(op,projection_bindings)

        if (
                degenerate_changes is not None
                and not isinstance(degenerate_changes, SideDegeneracyChanges)
                and _degeneracy_identification_context is None
                and cls.use_compiled_degeneracy_identification
        ):
            _degeneracy_identification_context = (
                DegeneracyIdentificationContext(degenerate_changes)
            )

        if evaluation_backend is None:
            evaluation_backend = cls.default_evaluation_backend
        if evaluation_backend not in {'legacy', 'indexed', 'compare'}:
            raise ValueError(
                "unknown analytic evaluation backend {}".format(evaluation_backend)
            )
        if evaluation_backend == 'compare':
            if parallelizer is not None:
                raise NotImplementedError(
                    "indexed/legacy comparison currently requires serial evaluation"
                )
            common_kwargs = dict(
                op=op,
                logger=logger,
                parallelizer=None,
                degenerate_changes=degenerate_changes,
                only_degenerate_terms=only_degenerate_terms,
                zero_cutoff=zero_cutoff,
                verbose=verbose,
                log_scaled=log_scaled,
                evaluation_mode=evaluation_mode,
                dag_cache_size=dag_cache_size,
                dag_cache_bytes=dag_cache_bytes,
                dag_chunk_size=dag_chunk_size,
                path_cache_size=path_cache_size,
                path_cache_bytes=path_cache_bytes,
                state_permutation_cache_size=state_permutation_cache_size,
                state_permutation_cache_bytes=state_permutation_cache_bytes,
                _degeneracy_identification_context=(
                    _degeneracy_identification_context
                ),
                mode_subset_selector=mode_subset_selector,
                dag_stream_options=dag_stream_options,
                _max_polynomial_order=_max_polynomial_order
            )
            legacy = cls.evaluate_polynomial_expression(
                state_perms, coeffs, freqs,
                expr, change, baseline_shift, num_fixed,
                evaluation_backend='legacy',
                **common_kwargs
            )
            indexed = cls.evaluate_polynomial_expression(
                state_perms, coeffs, freqs,
                expr, change, baseline_shift, num_fixed,
                evaluation_backend='indexed',
                **common_kwargs
            )
            cls._compare_evaluation_backends(legacy, indexed)
            return legacy
        use_materialized_path_cache = evaluation_mode in {
            'materialized', 'materialized_cached'
        }
        if evaluation_mode in {
            'materialized', 'materialized_cached', 'materialized_legacy'
        }:
            cls._last_materialized_evaluation_stats = None

        # ensure data is in a format where we can loop over states, perms, and degenerate
        # blocks all in one go
        # state = np.asanyarray(state)
        # smol = state.ndim == 1
        # if smol: state = state[np.newaxis]
        # ndim = state.shape[-1]

        state_perms = [
            [np.asanyarray(state), np.asanyarray(perms)]
            for state, perms in state_perms
        ]
        ndim = state_perms[0][0].shape[-1]

        all_perms, perm_map = np.unique(
            np.concatenate([perms for state, perms in state_perms], axis=0),
            axis=0, return_inverse=True
        )

        # if perms is None: perms = np.arange(ndim)
        # smol_p = perms.ndim == 1
        # if smol_p: perms = perms[np.newaxis]
        # all_p = perms.ndim == 2
        # if all_p: np.repeat(
        #         perms[np.new_axis],
        #         len(state),
        #         axis=0
        #     )

        # smol_coeffs = PerturbationTheoryEvaluator.is_single_expansion(coeffs)
        # if smol_coeffs: coeffs = [coeffs]
        if not isinstance(coeffs, SharedMemoryArrayTree):
            coeffs = cls._normalize_coefficient_expansions(coeffs)

        if only_degenerate_terms is False:
            expr, _ = cls._prune_resonant_dag(expr, state_perms, degenerate_changes,
                num_fixed, len(freqs), context=_degeneracy_identification_context)

        if isinstance(expr, PTTensorCoeffProductDAG):
            active_coefficients = {
                (order, coefficient_type)
                for expansion in coeffs
                for order, order_expansion in enumerate(expansion)
                for coefficient_type, coefficient in enumerate(order_expansion)
                if not nput.is_zero(coefficient)
            }
            inactive_coefficients = _HamiltonianCoefficient.inactive_operators(expr.operator_keys,active_coefficients)
            if len(inactive_coefficients) > 0:
                expr = expr.prune_operators(inactive_coefficients)

            if evaluation_mode is None:
                evaluation_mode = cls.default_tensor_evaluation_mode
            if evaluation_mode in {'dag', 'dag_legacy'}:
                return cls._evaluate_tensor_dag_expression(
                    state_perms, coeffs, freqs,
                    expr, change, baseline_shift,
                    num_fixed,
                    op=op,
                    logger=logger,
                    parallelizer=parallelizer,
                    degenerate_changes=degenerate_changes,
                    only_degenerate_terms=only_degenerate_terms,
                    zero_cutoff=zero_cutoff,
                    verbose=verbose, log_scaled=log_scaled,
                    cache_size=dag_cache_size,
                    cache_bytes=dag_cache_bytes,
                    chunk_size=dag_chunk_size,
                    state_cache_size=state_permutation_cache_size,
                    state_cache_bytes=state_permutation_cache_bytes,
                    use_path_cache=evaluation_mode == 'dag',
                    evaluation_backend=evaluation_backend,
                    degeneracy_identification_context=(
                        _degeneracy_identification_context
                    ),
                    mode_subset_selector=mode_subset_selector,
                    stream_options=dag_stream_options
                )
            elif evaluation_mode in {
                'materialized', 'materialized_cached', 'materialized_legacy'
            }:
                expr = expr.to_eager()
            else:
                raise ValueError("unknown tensor evaluation mode {}".format(evaluation_mode))

        # udegs, udeg_inv = np.unique(np.concatenate(degenerate_changes, axis=0), axis=0, return_inverse=True)

        if baseline_shift is not None and all(b == 0 for b in baseline_shift): baseline_shift = None

        # in general, states are mostly comprised of the same value
        # so why both recomputing the same value over and over?
        # ustate_data = [
        #     np.unique(s, return_counts=True, return_index=True)
        #     for s in state
        # ]

        # there are a few phases to the evaluation strategy
        # first we identify how many _free_ indices are required per type of term
        # then we take every combination of free indices, and for each of these
        # we want to be as efficient as possible
        # in particular, since states are mostly zeros, we can cache each type of polynomial
        # evaluation by determining what the quanta actually involved are

        logger = Logger.lookup(logger)
        log_level = Logger.LogLevel.Normal if verbose else Logger.LogLevel.Debug

        if nput.is_zero(expr) or len(expr.terms) == 0:
            res = tuple(
                np.zeros([len(coeffs), len(perms)])
                for state,perms in state_perms
            )
        else:
            execution = getattr(expr, '_compiled_execution', None)
            if execution is not None:
                _max_polynomial_order = execution['max_order']
            max_order = (cls._get_max_order(expr) if _max_polynomial_order is None
                         else _max_polynomial_order)
            if max_order < 0:
                raise Exception(expr)
            max_state = np.max(np.concatenate([state for state,perms in state_perms]))
            pows = np.power(np.arange(max_state+1)[np.newaxis, :], np.arange(max_order+1)[:, np.newaxis])

            if degenerate_changes is not None:
                # if len(degenerate_changes) != len(state_perms):
                #     raise ValueError("each state needs its own block of degeneracies")
                degenerate_changes = cls._identify_possible_degeneracies(
                    all_perms, change,
                    expr, degenerate_changes, len(freqs),
                    context=_degeneracy_identification_context
                )
                # import pprint
                # pprint.pprint([list(d.keys()) for d in degenerate_changes])
                # #     for dc in degenerate_changes
                # # ]
                if only_degenerate_terms is True:
                    only_degenerate_terms = [only_degenerate_terms, only_degenerate_terms, np.logical_or]
                elif only_degenerate_terms is False:
                    only_degenerate_terms = [only_degenerate_terms, only_degenerate_terms, np.logical_and]
                if len(only_degenerate_terms) == 2:
                    only_degenerate_terms = list(only_degenerate_terms) + [np.logical_and]
            else:
                only_degenerate_terms = [False, False, np.logical_and]

            # Group by numbers of coordinates to do the combinatorics over
            free_ind_groups = {}
            if only_degenerate_terms[0] and only_degenerate_terms[1]:
                if degenerate_changes is not None:
                    # A key with resonant denominators on both sides appears in
                    # both the left and the right dict; evaluate it once.
                    resonant_keys = dict.fromkeys(
                        coeff_indices
                        for subchanges in degenerate_changes
                        for coeff_indices in subchanges.keys()
                    )
                    for coeff_indices in resonant_keys:
                        num_inds = len(expr.get_inds(coeff_indices))
                        free_inds = num_inds - num_fixed
                        if free_inds not in free_ind_groups:
                            free_ind_groups[free_inds] = []
                        free_ind_groups[free_inds].append(coeff_indices)
            elif execution is not None:
                free_ind_groups = dict(execution['groups'])
            else:
                for coeff_indices, poly_terms in expr.terms.items():
                    num_inds = len(expr.get_inds(coeff_indices))
                    free_inds = num_inds - num_fixed
                    if free_inds not in free_ind_groups:
                        free_ind_groups[free_inds] = []
                    free_ind_groups[free_inds].append(coeff_indices)

            # # Get state counts so that we can avoid doing duplicate polynomial
            # # evaluations
            # state_vals = []
            # state_counds = []
            # for s in state:
            #     vals, counts, pos = np.unique()

            if zero_cutoff is None:
                zero_cutoff = cls.default_zero_cutoff

            indexed_executor = None
            eval_perm = cls._eval_perm
            indexed_evaluation_start = None
            if evaluation_backend == 'indexed':
                from .Evaluators import IndexedBlockEvaluator
                indexed_executor = IndexedBlockEvaluator(cls, expr, preparation_context=_stream_preparation)
                eval_perm = indexed_executor.eval_perm
                indexed_evaluation_start = time.perf_counter()
            parallel_indexed_stats = []

            if op is None: op = expr
            with logger.block(
                    tag="Evaluating {op}({ch})", op=op,
                    ch=format_quantum_change(change), log_level=log_level):

                # contrib = [np.zeros([len(p)]) for p in perms]
                contrib = tuple(
                    np.zeros([len(coeffs), len(perms)])
                    for state, perms in state_perms
                )
                counts_cache = cls.get_cache() # so we don't hit some coeff sets too many times
                if _poly_eval_cache is not None:
                    poly_cache = _poly_eval_cache
                elif use_materialized_path_cache:
                    path_cache_size = (
                        cls.default_materialized_path_cache_size
                        if path_cache_size is None else path_cache_size
                    )
                    path_cache_bytes = (
                        cls.default_materialized_path_cache_bytes
                        if path_cache_bytes is None else path_cache_bytes
                    )
                    poly_cache = _MaterializedEvaluationCache(
                        cls._poly_cache,
                        max_items=path_cache_size,
                        max_bytes=path_cache_bytes
                    )
                else:
                    poly_cache = cls._poly_cache
                # Avoid redoing complete polynomials and, in the cached
                # materialized mode, their shared PolyPath subexpressions.
                take_cache = (
                    cls._new_state_permutation_cache(
                        state_permutation_cache_size,
                        state_permutation_cache_bytes
                    )
                    if _state_permutation_cache is None else
                    _state_permutation_cache
                )
                energy_cache = cls.get_cache()
                with dev.context_wrap(parallelizer):
                    if parallelizer is not None:
                        worker_path_cache_size = (
                            cls.default_materialized_path_cache_size
                            if path_cache_size is None else path_cache_size
                        )
                        worker_path_cache_bytes = (
                            cls.default_materialized_path_cache_bytes
                            if path_cache_bytes is None else path_cache_bytes
                        )
                        # parallelizer.run(
                        #     cls._initialize_main_eval,
                        #     expr, change, baseline_shift,
                        #     state_perms, all_perms, perm_map,
                        #     freqs,
                        #     num_fixed, degenerate_changes, only_degenerate_terms,
                        #     zero_cutoff,
                        #     counts_cache, poly_cache, take_cache, energy_cache,
                        #     pows,
                        #     verbose, logger, log_scaled, log_level,
                        # )
                        nz_pos = [
                            np.nonzero(state)
                            for state, perms in state_perms
                        ]
                        compressed_states = [
                            [state[nz], nz]
                            for nz,(state, _) in zip(nz_pos, state_perms)
                        ]
                        perm_lens = np.array([
                            len(perms)
                            for state, perms in state_perms
                        ])
                        compressed_perms = mcomb.lehmer_encode(np.concatenate([
                            perms
                            for state, perms in state_perms
                        ], axis=0))
                        parallelizer.set_initializer(
                            cls._initialize_main_eval,
                            coeffs,  # Every lazy initializer generation must be self-contained.
                            expr, change, baseline_shift,
                            # state_perms, #all_perms, perm_map,
                            compressed_states, perm_lens, compressed_perms,
                            freqs,
                            num_fixed, degenerate_changes, only_degenerate_terms,
                            zero_cutoff,
                            cls._max_cache_size,
                            # The parent DAG cache is process-local. Spawned
                            # workers need their own bounded materialized/path
                            # cache rather than falling back to scalar lookup.
                            use_materialized_path_cache,
                            worker_path_cache_size, worker_path_cache_bytes,
                            state_permutation_cache_size,
                            state_permutation_cache_bytes,
                            # counts_cache, poly_cache, take_cache, energy_cache,
                            # pows,
                            verbose, logger, log_scaled, log_level,
                            (None if _stream_preparation is None else _stream_preparation.max_items),
                            ndim, evaluation_backend
                        )
                    # if cls._cached_expansion is None:
                    #     cls._cached_expansion = coeffs
                    # else:
                    #     coeffs = cls._cached_expansion
                    for free_inds, cind_sets in free_ind_groups.items():
                        # now we iterate over every subset of inds (outside of the fixed ones)
                        # excluding replacement (since that corresponds to a different coeff/poly)
                        # from McUtils.Profilers import Timer
                        # with Timer(free_inds):
                            tensors, cind_specs = (cls._prep_coeff_data(cind_sets, coeffs)
                                if _stream_preparation is None else
                                _stream_preparation.prepare(cls, cind_sets))
                            if len(tensors) > 0:

                                facs = execution['factorials'] if execution is not None else [
                                    math.factorial(n)
                                    for n in range(max([len(cinds) for cinds in cind_specs] + [0]) + 1)
                                ]
                                split_spec = np.cumsum([0] + [len(perms) for state, perms in state_perms])[1:]

                                ncomb = (math.comb(ndim-num_fixed,free_inds)
                                    if mode_subset_selector is None else
                                    mode_subset_selector.count(ndim,num_fixed,free_inds))
                                if ncomb == 0:
                                    continue
                                if ncomb > 1 and parallelizer is not None:
                                    parallelizer: Parallelizer
                                    contrib_shapes = tuple(
                                        [len(coeffs), len(perms)]
                                        for state, perms in state_perms
                                    )

                                    nproc = parallelizer.nprocs
                                    block_ranges = cls._partition_evaluation_blocks(
                                        ncomb, nproc
                                    )

                                    subcontrib, indexed_stats = parallelizer.run(
                                        cls._run_eval_blocks,
                                        block_ranges,
                                        contrib_shapes,
                                        free_inds,
                                        cind_sets,
                                        mode_subset_selector,
                                        cleanup=False
                                    )
                                    if indexed_stats is not None:
                                        parallel_indexed_stats.append(
                                            indexed_stats
                                        )
                                    for storage,corr in zip(contrib, subcontrib):
                                        storage += corr
                                else:
                                    for subset in (
                                            mode_subset_selector.iter_subsets(ndim,num_fixed,free_inds)
                                                if mode_subset_selector is not None else
                                            itertools.combinations(range(num_fixed,ndim),r=free_inds)
                                    ):
                                        if verbose:
                                            with logger.block(tag="{b} + {s}", b=tuple(range(num_fixed)), s=subset,
                                                              log_level=log_level):
                                                subcontrib = eval_perm(
                                                    expr, change, baseline_shift,
                                                    subset, state_perms, all_perms, perm_map,
                                                    freqs, cind_specs, tensors,
                                                    num_fixed, degenerate_changes, only_degenerate_terms,
                                                    zero_cutoff,
                                                    counts_cache, poly_cache, take_cache, energy_cache,
                                                    facs, pows, split_spec,
                                                    verbose, logger, log_scaled,
                                                    evaluation_context=op
                                                )
                                        else:
                                            subcontrib = eval_perm(
                                                expr, change, baseline_shift,
                                                subset, state_perms, all_perms, perm_map,
                                                freqs, cind_specs, tensors,
                                                num_fixed, degenerate_changes, only_degenerate_terms,
                                                zero_cutoff,
                                                counts_cache, poly_cache, take_cache, energy_cache,
                                                facs, pows, split_spec,
                                                verbose, logger, log_scaled,
                                                evaluation_context=op
                                            )
                                        for storage,corr in zip(contrib, subcontrib):
                                            storage += corr
            res = [expr.prefactor * c for c in contrib]

            if Diagnostics.enabled and isinstance(poly_cache, _MaterializedEvaluationCache):
                cls._last_materialized_evaluation_stats = poly_cache.stats()
            if Diagnostics.enabled and isinstance(take_cache, BoundedCache):
                cls._last_state_permutation_cache_stats = take_cache.stats()
            if Diagnostics.enabled and indexed_executor is not None:
                indexed_stats = [indexed_executor.finish()]
                indexed_stats.extend(parallel_indexed_stats)
                cls._last_indexed_evaluation_stats = (
                    cls._aggregate_indexed_evaluation_stats(indexed_stats)
                )
                cls._last_indexed_evaluation_stats['elapsed_seconds'] = (
                    time.perf_counter() - indexed_evaluation_start
                )
                cls._indexed_evaluation_generation += 1

        # if smol_coeffs: res = [r[0] for r in res]
        # if smol: res = res[0]

        return res



    def evaluate(self, state_perms, coeffs, freqs, degenerate_changes=None, only_degenerate_terms=False,
                 zero_cutoff=None, parallelizer=None,
                 verbose=False, log_scaled=True,
                 evaluation_mode=None,
                 evaluation_backend=None,
                 dag_cache_size=None,
                 dag_cache_bytes=None,
                 dag_chunk_size=None,
                 path_cache_size=None,
                 path_cache_bytes=None,
                 state_permutation_cache_size=None,
                 state_permutation_cache_bytes=None,
                 use_shared_memory=None,
                 _shared_coefficients=None,
                 _degeneracy_identification_context=None,
                 mode_subset_selector=None,
                 dag_stream_options=None, projection_bindings=None
                 ):

        from .Evaluators import _current_evaluation_mode, _validate_evaluation_mode
        evaluation_mode = _validate_evaluation_mode(
            _current_evaluation_mode.get() if evaluation_mode is None else evaluation_mode)
        use_actions = evaluation_mode == 'actions'
        if evaluation_mode in {'actions', 'direct'}:
            evaluation_mode = None  # compile selects the numerical backend

        # state = np.asanyarray(state)
        smol = nput.is_numeric(state_perms[0][0])
        if smol: state_perms = [state_perms]

        if parallelizer is not None:
            parallelizer = Parallelizer.lookup(parallelizer)


        smol_coeffs = PerturbationTheoryEvaluator.is_single_expansion(coeffs)
        if _shared_coefficients is not None:
            # The owner supplies the already-batched tree.  Keep
            # ``smol_coeffs`` from the public input so the result shape is
            # unchanged for a single expansion.
            coeffs = _shared_coefficients
        elif smol_coeffs:
            coeffs = [coeffs]

        from .Evaluators import current_compilation, _active_coefficient_operators
        compiled_context = current_compilation()
        def evaluate_compiled(compiled):
            _, _, task_options, default_parallelizer = compiled_context
            numerical = dict(projection_bindings=projection_bindings,zero_cutoff=zero_cutoff, degenerate_changes=degenerate_changes,
                only_degenerate_terms=only_degenerate_terms, verbose=verbose, log_scaled=log_scaled)
            for key,value in dict(evaluation_backend=evaluation_backend,path_cache_size=path_cache_size,
                    path_cache_bytes=path_cache_bytes,state_permutation_cache_size=state_permutation_cache_size,
                    state_permutation_cache_bytes=state_permutation_cache_bytes).items():
                if value is not None:numerical[key] = value
            numerical.update(task_options)
            numerical.setdefault('use_shared_memory', use_shared_memory)
            with _log_evaluation_time(self.logger, "evaluating compiled expression..."):
                res = compiled.evaluate(state_perms, coeffs, freqs,
                    mode_subset_selector=mode_subset_selector,
                    parallelizer=parallelizer if parallelizer is not None else default_parallelizer, **numerical)
            if smol_coeffs:res = [r[0] for r in res]
            if smol:res = res[0]
            return res

        if (compiled_context is not None and compiled_context[0].using_artifacts
                and evaluation_mode in {None, 'compiled'}):
            # File libraries contain full, system-independent programs. Execute
            # those programs directly, without controller-side physical pruning
            # or symbolic occupation-coefficient projection.
            store, compile_options, _, _ = compiled_context
            stored = self.__dict__.get('_expr', self.__dict__.get('expr'))
            change = stored.poly_change if isinstance(stored, SqrtChangePoly) else self.change
            baseline = (stored.shift_start if isinstance(stored, SqrtChangePoly)
                        else getattr(self, '_expression_shift', None))
            context = dict(operator=self.op if self.op is not None else self.expr,
                change=tuple(change), baseline=None if baseline is None else tuple(baseline))
            compiled = store.compile_factory(lambda: self, context=context, **compile_options)
            return evaluate_compiled(compiled)

        expr = self.expr
        if nput.is_zero(expr):
            res = [
                np.zeros([len(coeffs), len(p)])
                for s,p in state_perms
            ]
        else:
            if isinstance(expr, SqrtChangePoly):
                poly_obj = expr.poly_obj
                change = expr.poly_change
                shift_start = expr.shift_start
            else:
                poly_obj = expr
                change = None
                shift_start = None

            if not isinstance(poly_obj, PTTensorCoeffProductSum):
                poly_obj = PTTensorCoeffProductSum({():poly_obj})

            if change is None:
                first_echange_poly = next(iter(poly_obj.terms.values()))
                if not isinstance(first_echange_poly, PTEnergyChangeProductSum):
                    change = []
                else:
                    first_ekey = next(iter(first_echange_poly.terms.keys()))
                    change = [] if not first_ekey else np.sum([list(k) for k in first_ekey], axis=1)

            from .Evaluators import current_compilation, _active_coefficient_operators
            compiled_context = current_compilation()
            compiled_expression = self
            pruning_signature = b''
            from .Evaluators import AnalyticCorrectionCompiler
            if use_actions:
                numerical = dict(degenerate_changes=degenerate_changes,
                    only_degenerate_terms=only_degenerate_terms,projection_bindings=projection_bindings,
                    zero_cutoff=zero_cutoff,parallelizer=parallelizer,verbose=verbose,log_scaled=log_scaled,
                    evaluation_mode=evaluation_mode,evaluation_backend=evaluation_backend,
                    dag_cache_size=dag_cache_size,dag_cache_bytes=dag_cache_bytes,dag_chunk_size=dag_chunk_size,
                    path_cache_size=path_cache_size,path_cache_bytes=path_cache_bytes,
                    state_permutation_cache_size=state_permutation_cache_size,
                    state_permutation_cache_bytes=state_permutation_cache_bytes,
                    use_shared_memory=use_shared_memory,dag_stream_options=dag_stream_options,
                    _degeneracy_identification_context=_degeneracy_identification_context)
                res = AnalyticCorrectionCompiler.evaluate_action_expansion(
                    compiled_expression,state_perms,coeffs,freqs,
                    pruning_signature=pruning_signature,**numerical)
                if smol_coeffs:res = [r[0] for r in res]
                return res[0] if smol else res

            if isinstance(poly_obj, PTTensorCoeffProductDAG) and poly_obj.has_change_domains:
                if projection_bindings is None:
                    raise ValueError('filtered block equations require resonance-change bindings')
                from .Evaluators import _current_evaluation_plan
                plan = _current_evaluation_plan.get()
                enabled = True if plan is None else plan[1].root.options.get('prune_fixed_resonances', True)
                if not isinstance(enabled, (bool, np.bool_)):
                    raise TypeError('prune_fixed_resonances must be boolean')
                if enabled:
                    statistics = {} if Diagnostics.enabled else None
                    poly_obj, pruning_signature = poly_obj.select_change_domains(
                        state_perms, projection_bindings['projectors'], self.num_fixed,
                        len(freqs), baseline=shift_start, statistics=statistics)
                    if Diagnostics.enabled:
                        type(self)._last_change_domain_pruning_stats = statistics
                    if pruning_signature:
                        pruned_expr = (0 if nput.is_zero(poly_obj) else
                            expr.mutate(poly_obj=poly_obj) if isinstance(expr, SqrtChangePoly) else poly_obj)
                        compiled_expression = type(self)(self.op, pruned_expr, self.change, logger=self.logger)
                        if nput.is_zero(poly_obj):
                            result = [np.zeros((len(coeffs), len(perms))) for _,perms in state_perms]
                            if smol_coeffs: result = [row[0] for row in result]
                            return result[0] if smol else result
            if (compiled_context is not None and evaluation_mode in {None, 'compiled'}
                    and only_degenerate_terms is False and degenerate_changes is not None
                    and self.num_fixed and isinstance(poly_obj, PTTensorCoeffProductDAG)):
                with _log_evaluation_time(self.logger, "Identifying fixed-mode resonances..."):
                    pruned, pruning_signature = self._prune_resonant_dag(
                        poly_obj, state_perms, degenerate_changes, self.num_fixed,
                        len(freqs), context=_degeneracy_identification_context)
                if pruning_signature:
                    pruned_expr = (0 if nput.is_zero(pruned) else
                                   expr.mutate(poly_obj=pruned) if isinstance(expr, SqrtChangePoly) else pruned)
                    compiled_expression = type(self)(self.op, pruned_expr, self.change,
                                                     logger=self.logger)

            if compiled_context is not None and evaluation_mode in {None, 'compiled'}:
                store, compile_options, task_options, default_parallelizer = compiled_context
                compile_options = dict(compile_options)
                compile_options.setdefault('mode_dimension', len(freqs))
                active = _active_coefficient_operators(coeffs)
                if projection_bindings is not None:
                    auxiliary=projection_bindings['scalar_order']
                    active=tuple(sorted(set(active)|{(auxiliary+k,0) for k in range(projection_bindings['scalar_count'])}))
                with _log_evaluation_time(self.logger, "Preparing compiled expression..."):
                    compilation_context = dict(operator=self.op if self.op is not None else poly_obj,
                        change=tuple(change), baseline=None if shift_start is None else tuple(shift_start))
                    if pruning_signature:
                        compilation_context['excluded_resonances'] = pruning_signature.hex()
                    compiled = store.compile(compiled_expression, context=compilation_context,
                                             active_operators=active, **compile_options)
                return evaluate_compiled(compiled)

            shared_coefficients = None
            evaluator_cls = type(self)
            if isinstance(parallelizer, MultiprocessingParallelizer):
                if (
                        use_shared_memory is not False
                        and not isinstance(coeffs, SharedMemoryArrayTree)
                ):
                    # Standalone evaluator calls own one shared tree. Batch
                    # callers can pass ``_shared_coefficients`` to reuse a
                    # single arena across all correction expressions.
                    evaluator_cls._cached_expansion = None
                    coeffs = evaluator_cls._normalize_coefficient_expansions(
                        coeffs
                    )
                    shared_coefficients = SharedMemoryArrayTree.create(coeffs)
                    coeffs = shared_coefficients

            try:
                res = self.evaluate_polynomial_expression(
                    state_perms, coeffs, freqs,
                    poly_obj, change, shift_start,
                    self.num_fixed,
                    op=self.op,
                    logger=self.logger,
                    parallelizer=parallelizer,
                    zero_cutoff=zero_cutoff,
                    degenerate_changes=degenerate_changes,
                    only_degenerate_terms=only_degenerate_terms,
                    verbose=verbose, log_scaled=log_scaled,
                    evaluation_mode=evaluation_mode,
                    evaluation_backend=evaluation_backend,
                    dag_cache_size=dag_cache_size,
                    dag_cache_bytes=dag_cache_bytes,
                    dag_chunk_size=dag_chunk_size,
                    path_cache_size=path_cache_size,
                    path_cache_bytes=path_cache_bytes,
                    state_permutation_cache_size=state_permutation_cache_size,
                    state_permutation_cache_bytes=state_permutation_cache_bytes,
                    _degeneracy_identification_context=(
                        _degeneracy_identification_context
                    ),
                    mode_subset_selector=mode_subset_selector,
                    dag_stream_options=dag_stream_options, projection_bindings=projection_bindings
                )
            finally:
                if shared_coefficients is not None:
                    evaluator_cls._cached_expansion = None
                    shared_coefficients.dispose()

        if smol_coeffs: res = [r[0] for r in res]
        if smol: res = res[0]

        return res

    # def subtitute_coeffs(self, coeffs):
    #     """
    #     :param coeffs: V expansion, G expansion, U expansion, and then Coriolis expansions
    #     :return:
    #     """
    #     if not isinstance(self.expr, PTTensorCoeffProductSum):
    #         raise ValueError("coeffs have already been subbed")
    #     new_terms = []
    #     for k,p in self.expr.terms.items(): # this won't fully work since I need to map over possible pairs...
    #         c = coeffs[k[1]][k[0]]
    #         if len(k) > 2:
    #             c = c[k[2:]]
    #         if c != 0:
    #             new_terms.append(p*c)
    #     return type(self)(new_terms)
    #
    # def substitute_freqs(self, freqs):


def _compiled_evaluation_context(method):
    import functools
    signature = inspect.signature(method)
    @functools.wraps(method)
    def wrapped(self, *args, **kwargs):
        from .Evaluators import (compiled_evaluation_context, evaluation_mode_context,
                                 current_compilation, _current_evaluation_plan)
        arguments = signature.bind(self, *args, **kwargs).arguments
        mode = arguments.get('evaluation_mode', arguments.get('opts', {}).get('evaluation_mode'))
        with evaluation_mode_context(mode):
            active_plan = _current_evaluation_plan.get() is not None
            active_compilation = current_compilation() is not None
            if (self.compile_expressions and active_compilation
                    or not self.compile_expressions and active_plan and not active_compilation):
                return method(self, *args, **kwargs)
            options = dict(self.evaluation_options)
            if active_plan:
                options.pop('evaluation_mode', None)
            if mode is not None:
                options['evaluation_mode'] = mode
            with compiled_evaluation_context(
                    self.expression_store if self.compile_expressions else None,
                    self.compilation_options, options, self.parallelizer):
                return method(self, *args, **kwargs)
    return wrapped

@contextlib.contextmanager
def _log_evaluation_time(logger, tag, **options):
    """Use the analytic evaluator's normal block/elapsed-time convention."""
    with logger.block(tag=tag, **options):
        start = time.perf_counter()
        yield
        logger.log_print('took {e:.3f}s', e=time.perf_counter() - start,
                         log_level=options.get('log_level', Logger.LogLevel.Normal))


def _timed_evaluation(tag):
    """Time a complete call using its logger override and the usual block format."""
    def decorate(method):
        signature = inspect.signature(method)
        @functools.wraps(method)
        def wrapped(self, *args, **kwargs):
            arguments = signature.bind(self, *args, **kwargs).arguments
            options = arguments.get('opts', {})
            logger = arguments.get('logger', options.get('logger'))
            if logger is None:
                logger = getattr(self, 'logger', None)
                if logger is None and hasattr(self, 'solver'):
                    logger = self.solver.logger
            timing_options = {}
            if arguments.get('log_level') is not None:
                timing_options['log_level'] = arguments['log_level']
            with _log_evaluation_time(Logger.lookup(logger), tag, **timing_options):
                return method(self, *args, **kwargs)
        return wrapped
    return decorate


# ---------------------------------------------------------------------------------------------
#  Projected (formerly VPT2/Projected.py)
#  Block-projected analytic perturbation theory.
#
#  Finite model-space contractions stay in this correction layer. Sums over
#  the complementary oscillator space use the ordinary analytic term algebra.
# ---------------------------------------------------------------------------------------------

_projected_evaluators = weakref.WeakSet()


def _clear_projected_caches():
    ValidResolventTerm._cache.clear()
    ProjectedTerm._cache.clear()
    _BlockScalarTerm._cache.clear()
    _projection_reference_data.cache_clear()
    _bound_projection_row.cache_clear()
    for evaluator in list(_projected_evaluators):
        if hasattr(evaluator,'_projected_block_cache'):
            del evaluator._projected_block_cache
        if hasattr(evaluator,'_projection_runtime_inputs'):
            del evaluator._projection_runtime_inputs


def coefficient_fingerprint(expansions):
    digest = hashlib.sha256()
    for order in expansions:
        for coefficient in order:
            value = np.asarray(coefficient)
            digest.update(repr((value.shape,value.dtype.str)).encode())
            digest.update(value.tobytes() if value.dtype != object else repr(coefficient).encode())
    return digest.digest()


class ProjectionFactor(int):
    """A formal projection/reference role, containing no numerical bindings."""
    _state_projection_factor = True

    def __new__(cls, kind, slot=0):
        if kind not in {'P','Q','R','V'} or not isinstance(slot, int) or slot < 0:
            raise ValueError('projection factors require a kind and nonnegative formal slot')
        obj = int.__new__(cls, 2 + 4*slot + {'P':0,'Q':1,'R':2,'V':3}[kind])
        obj.kind, obj.slot = kind, slot
        return obj

    def __reduce__(self):
        return type(self), (self.kind,self.slot)

    def __repr__(self):
        return 'Projection[{}:{}]'.format(self.kind,self.slot)

    @property
    def label(self):
        return '{}[{}]'.format(self.kind,self.slot)


class _ProjectionEvaluationContext:
    """Numerical call data; never a term or a symbolic compilation input."""
    def __init__(self,operator,bindings):
        self.operator,self.projection_bindings = operator,bindings
    def __repr__(self):
        return repr(self.operator)
    def for_route(self, route):
        if 'route_projectors' not in self.projection_bindings:
            return self
        bindings = dict(self.projection_bindings)
        bindings['projectors'] = projection_inputs(self, route)
        bindings.pop('route_projectors')
        return type(self)(self.operator, bindings)


def projection_inputs(context, route=None):
    if not isinstance(context, _ProjectionEvaluationContext):
        return ()
    bindings = context.projection_bindings
    if route is not None and 'route_projectors' in bindings:
        # Multi-source concrete blocks contain only ordinary terms. Diagnostic
        # formatting for those blocks uses their first source as before.
        if isinstance(route, tuple):
            route = route[0]
        return bindings['route_projectors'][int(route)]
    return bindings.get('projectors', ())


class _BlockScalarTerm(PerturbationTheoryTerm):
    """Formal finite-block coefficient, supplied through numerical tensors."""
    _cache = {}
    def __init__(self,order,slot,logger=None,representation='path'):
        super().__init__(logger=logger,polynomial_representation=representation)
        self.order,self.slot = int(order),int(slot)
    @classmethod
    def lookup(cls,order,slot,logger=None,representation='path'):
        key = order,slot,representation
        if key not in cls._cache:cls._cache[key]=cls(order,slot,logger,representation)
        return cls._cache[key]
    def __repr__(self):
        return 'BlockCoefficient[{},{}]'.format(self.order,self.slot)
    def get_changes(self):return {():None}
    def get_core_poly(self,changes,shift=None):
        if changes:return 0
        poly = make_product_polynomial([],steps=0,representation=self.polynomial_representation)
        # Coefficient identity 1 denotes a kinetic tensor in native algebra;
        # use separate auxiliary rows, each with scalar identity 0.
        return SqrtChangePoly(PTTensorCoeffProductSum({((self.order+self.slot,0),):poly}),[],[])


class ValidResolventTerm(PerturbationOperator):
    """Occupation-validity marker without a system-specific model space."""
    _cache = {}
    def get_poly_terms(self, changes, shift=None, **opts):
        return PerturbationTheoryTerm.get_poly_terms(self, changes, shift=shift, **opts)
    def get_core_poly(self, changes, shift=None):
        base = super().get_core_poly(changes, shift=shift)
        if np.isscalar(base):return base
        valid = ProjectionFactor('V')
        return base.mul_simple(PTEnergyChangeProductSum.monomial((valid,)+tuple(changes),1))


class ProjectedTerm(PerturbationTheoryTerm):
    """A generic projection operator. Its role is bound only by evaluation."""
    _cache = {}
    def __init__(self,term,factor):
        super().__init__(logger=term.logger,polynomial_representation=term.polynomial_representation)
        self.term,self.factor = term,factor
    @classmethod
    def lookup(cls,term,factor):
        key=term,factor
        if key not in cls._cache:cls._cache[key]=cls(term,factor)
        return cls._cache[key]
    def __repr__(self):return '{}<{}>'.format(self.factor.label,repr(self.term))
    def get_changes(self):
        return {change:None for change in self.term.changes if change or self.factor.kind=='P'}
    def get_core_poly(self,changes,shift=None):
        base=self.term.get_poly_terms(changes,shift=shift)
        if np.isscalar(base) and base==0:return 0
        return base.mul_simple(PTEnergyChangeProductSum.monomial((self.factor,)+tuple(changes),1))


class IdentityTerm(PerturbationTheoryTerm):
    def __init__(self, logger=None, polynomial_representation='path'):
        super().__init__(logger=logger, polynomial_representation=polynomial_representation)

    def __repr__(self):
        return 'I'

    def get_changes(self):
        return {(): None}

    def get_core_poly(self, changes, shift=None):
        return SqrtChangePoly(PTTensorCoeffProductSum({(): make_product_polynomial(
            [], steps=0, representation=self.polynomial_representation
        )}), [], [])


@functools.lru_cache(maxsize=512)
def _projection_reference_data(binding, state):
    states,frequencies,reference,trivial = binding
    differences = np.asarray(states, dtype=int).reshape((-1, len(state))) - np.asarray(state)
    counts = np.count_nonzero(differences, axis=1)
    gap = (0. if reference is None else
           float(np.dot(np.asarray(reference)-state, frequencies)))
    return differences, counts, gap


def projected_path_factors(energy_key, state, permutations, perm_freqs, baseline_shift=None,
                           *, projection_bindings=()):
    """Return a projector mask and allowed denominators before any division."""
    if any(getattr(row[0],'_conditional_endpoint',False) for row in energy_key):
        return conditional_path_factors(energy_key,state,permutations,perm_freqs,baseline_shift,
            projection_bindings=projection_bindings)
    if any(isinstance(row[0], _ChangeDomainIndex) for row in energy_key):
        return _filtered_change_factors(energy_key, state, permutations, perm_freqs,
                                        baseline_shift, projection_bindings)
    state = np.asarray(state)[0]
    permutations = np.asarray(permutations, dtype=int)
    good = np.ones(len(permutations), dtype=bool)
    denominator = np.ones(len(permutations))
    for row in energy_key:
        factor, delta = row[0], np.asarray(row[1:], dtype=int)
        if not isinstance(factor, ProjectionFactor):
            denominator *= perm_freqs @ delta
            continue
        if factor.kind == 'V':
            good &= np.all(state[permutations] + delta >= 0, axis=1)
            continue
        binding = projection_bindings[factor.slot]
        if binding[3]:
            good &= np.all(state[permutations] + delta >= 0, axis=1)
            reference_gap = (0. if binding[2] is None else
                float(np.dot(np.asarray(binding[2])-state,binding[1])))
            denominator *= reference_gap - perm_freqs @ delta
            continue
        differences, counts, reference_gap = _projection_reference_data(binding, tuple(state))
        local = differences[:, permutations]
        matches = np.any(np.all(local == delta, axis=2)
                         & (np.count_nonzero(local, axis=2) == counts[:, None]), axis=0)
        good &= matches if factor.kind == 'P' else ~matches
        good &= np.all(state[permutations] + delta >= 0, axis=1)
        if factor.kind == 'R':
            denominator *= reference_gap - perm_freqs @ delta
    return good, denominator


def _compile_energy_program(key, baseline=None):
    """Validate formal factors once and lower their constant occupation data."""
    operations, position = [], 0
    ordinary = not any(getattr(row[0], '_state_projection_factor', False) for row in key)
    pure_projectors = not any(isinstance(row[0], (ConditionEndpoint, _ChangeDomainIndex)) for row in key)
    while position < len(key):
        row = key[position]
        factor, delta = row[0], np.asarray(row[1:])
        if isinstance(factor, (ConditionEndpoint, _ChangeDomainIndex)):
            if position+1 >= len(key): raise ValueError('incomplete compiled endpoint pair')
            following = key[position+1]
            end_factor = following[0]
            if (type(end_factor) is not type(factor) or factor.endpoint != 0 or end_factor.endpoint != 1
                    or factor.retain != end_factor.retain):
                raise ValueError('invalid compiled endpoint pair')
            offset = np.zeros(len(delta), dtype=int)
            if baseline is not None: offset[:min(len(offset), len(baseline))] = baseline[:len(offset)]
            start, end = delta+offset, np.asarray(following[1:])+offset
            if isinstance(factor, ConditionEndpoint):
                if factor.kind != end_factor.kind or factor.condition != end_factor.condition:
                    raise ValueError('mismatched compiled conditional endpoints')
                operations.append(('condition', factor.condition, factor.kind, factor.retain, start, end, start-end))
            else:
                operations.append(('domain', factor.retain, start, end, bool(np.array_equal(start, end))))
            position += 2
        elif isinstance(factor, ProjectionFactor):
            operations.append(('projector', factor.kind, factor.slot, delta))
            position += 1
        else:
            if getattr(factor, '_state_projection_factor', False):
                raise TypeError('unknown compiled projection factor')
            operations.append(('harmonic', delta))
            position += 1
    program = dict(version=1, operations=tuple(operations), ordinary=ordinary, pure_projectors=pure_projectors)
    if pure_projectors:
        program['deltas'] = np.asarray([row[1:] for row in key])
        program['kinds'] = np.asarray([{'R':1, 'P':2, 'Q':3, 'V':4}[row[0].kind]
                                      if isinstance(row[0], ProjectionFactor) else 0 for row in key], dtype=np.int64)
        program['slots'] = tuple(row[0].slot if isinstance(row[0], ProjectionFactor) else -1 for row in key)
        program['projection_indices'] = tuple(i for i, kind in enumerate(program['kinds']) if kind in {1, 2, 3})
        program['energy_rows'] = np.asarray(key) if ordinary else None
    return program


def _evaluate_energy_program(program, state, permutations, frequencies, bindings):
    state, permutations = np.asarray(state)[0], np.asarray(permutations, dtype=int)
    good, denominator = np.ones(len(permutations), dtype=bool), np.ones(len(permutations))
    for operation in program['operations']:
        kind, *arguments = operation
        if kind == 'harmonic': denominator *= frequencies @ arguments[0]
        elif kind == 'projector':
            role, slot, delta = arguments
            possible = np.all(state[permutations]+delta >= 0, axis=-1)
            if role == 'V':
                good &= possible
                continue
            binding = bindings[slot]
            if binding[3]:
                gap = 0. if binding[2] is None else float(np.dot(np.asarray(binding[2])-state, binding[1]))
                good &= possible
                denominator *= gap - frequencies @ delta
                continue
            differences, counts, gap = _projection_reference_data(binding, tuple(state))
            local = differences[:, permutations]
            members = np.any(np.all(local == delta, axis=-1)
                             & (np.count_nonzero(local, axis=-1) == counts[:, None]), axis=0)
            good &= possible & (members if role == 'P' else ~members)
            if role == 'R': denominator *= gap - frequencies @ delta
        elif kind == 'condition':
            condition, role, retain, start, end, gap = arguments
            possible = np.all(state[permutations]+start >= 0, axis=-1) & np.all(state[permutations]+end >= 0, axis=-1)
            members = condition.evaluate_bound(state, permutations, start, end, frequencies, bindings)
            good &= possible & (members if retain else ~members)
            if role == 'resolvent': denominator *= frequencies @ gap
        elif kind == 'domain':
            retain, start, end, identical = arguments
            possible = np.all(state[permutations]+start >= 0, axis=-1) & np.all(state[permutations]+end >= 0, axis=-1)
            if identical: members = np.ones(len(permutations), dtype=bool)
            else:
                if not bindings: raise ValueError('filtered block contraction requires resonance-change bindings')
                labels = _filtered_block_labels(bindings[0])
                origins = np.broadcast_to(state, (len(permutations), len(state))).copy()
                targets = origins.copy()
                rows = np.arange(len(permutations))[:, None]
                origins[rows, permutations] += start
                targets[rows, permutations] += end
                members = np.fromiter((source in labels and labels[source] == labels.get(tuple(target))
                    for source, target in zip(map(tuple, origins), targets)), dtype=bool, count=len(permutations))
            good &= possible & (members if retain else ~members)
        else: raise ValueError('unknown compiled denominator operation: ' + kind)
    return good, denominator


def _bind_projection_program(program, state, model_size, bindings):
    deltas, kinds = program['deltas'], program['kinds'].copy()
    gaps = np.zeros(len(deltas))
    differences = np.zeros((len(deltas), model_size, len(state)), dtype=int)
    counts = np.full((len(deltas), model_size), -1, dtype=int)
    for i in program['projection_indices']:
        binding = bindings[program['slots'][i]]
        if binding[3]:
            kinds[i] = 5
            gaps[i] = 0. if binding[2] is None else float(np.dot(np.asarray(binding[2])-state, binding[1]))
        else:
            delta, number, gaps[i] = _projection_reference_data(binding, state)
            differences[i,:len(number)], counts[i,:len(number)] = delta, number
    return deltas, kinds, gaps, differences, counts


@functools.lru_cache(maxsize=512)
def _bound_projection_row(energy_key, state, model_size, projection_bindings):
    factor_count = len(energy_key)
    deltas = np.asarray([row[1:] for row in energy_key], dtype=int)
    kinds = np.zeros(factor_count, dtype=int)
    gaps = np.zeros(factor_count)
    differences = np.zeros((factor_count, model_size, len(state)), dtype=int)
    counts = np.full((factor_count, model_size), -1, dtype=int)
    for i,row in enumerate(energy_key):
        factor = row[0]
        if isinstance(factor, ProjectionFactor):
            kinds[i] = {'R':1,'P':2,'Q':3,'V':4}[factor.kind]
            if factor.kind == 'V':continue
            binding = projection_bindings[factor.slot]
            if binding[3]:
                kinds[i] = 5
                gaps[i] = (0. if binding[2] is None else
                    float(np.dot(np.asarray(binding[2])-state,binding[1])))
                continue
            delta, number, gaps[i] = _projection_reference_data(binding,state)
            differences[i,:len(number)] = delta
            counts[i,:len(number)] = number
    return deltas,kinds,gaps,differences,counts


def projected_energy_rows(energy_keys, source_states, permutations, frequency_blocks, baseline_shift=None,
                          *, projection_bindings=(), row_projection_bindings=None, preparation_context=None,
                          programs=None):
    """Contract a bounded indexed row batch, retaining all projector tests."""
    bindings = ([projection_bindings] * len(energy_keys) if row_projection_bindings is None
                else row_projection_bindings)
    if len(bindings) != len(energy_keys):
        raise ValueError('one projector binding is required per energy row')
    if programs is not None:
        if len(programs) != len(energy_keys) or any(program is None for program in programs):
            raise ValueError('incomplete compiled denominator programs')
        if not all(program['pure_projectors'] for program in programs):
            results = [_evaluate_energy_program(program, state, perms, freqs, binding)
                       for program, state, perms, freqs, binding in zip(programs, source_states,
                           permutations, frequency_blocks, bindings)]
            return tuple(np.stack([result[k] for result in results]) for k in range(2))
    if programs is None and any(isinstance(row[0], _ChangeDomainIndex) for key in energy_keys for row in key):
        results = [_filtered_change_factors(key, state, perms, freqs, baseline_shift, binding)
                   for key, state, perms, freqs, binding in zip(
                       energy_keys, source_states, permutations, frequency_blocks, bindings)]
        return tuple(np.stack([result[k] for result in results]) for k in range(2))
    if programs is None and any(getattr(row[0],'_conditional_endpoint',False) for key in energy_keys for row in key):
        if row_projection_bindings is None:
            return conditional_energy_rows(energy_keys,source_states,permutations,frequency_blocks,baseline_shift,
                projection_bindings=projection_bindings)
        results = [conditional_energy_rows([key],[state],[perms],freqs[None],baseline_shift,
                    projection_bindings=binding)
                   for key,state,perms,freqs,binding in zip(
                       energy_keys,source_states,permutations,frequency_blocks,bindings)]
        return tuple(np.concatenate([result[k] for result in results]) for k in range(2))
    model_size = (max([1]+[len(binding[program['slots'][i]][0])
        for program, binding in zip(programs, bindings) for i in program['projection_indices']
        if not binding[program['slots'][i]][3]]) if programs is not None else max([1]+[len(binding[row[0].slot][0])
        for key,binding in zip(energy_keys,bindings) for row in key if isinstance(row[0],ProjectionFactor)
        and row[0].kind!='V' and not binding[row[0].slot][3]]))
    # Bound the projector-comparison scratch independently of the native
    # polynomial row chunk, whose bound does not include model-space size.
    factors = len(energy_keys[0])
    mode_count = len(np.asarray(source_states[0])[0])
    permutation_count, index_count = frequency_blocks.shape[1:]
    bytes_per_row = 8*factors*model_size*(mode_count+2*permutation_count*(index_count+1))
    row_limit = max(1,(16*1024**2)//max(1,bytes_per_row))
    if len(energy_keys) > row_limit:
        results = [projected_energy_rows(energy_keys[i:i+row_limit],
            source_states[i:i+row_limit],permutations[i:i+row_limit],frequency_blocks[i:i+row_limit],
            baseline_shift,projection_bindings=projection_bindings,
            row_projection_bindings=bindings[i:i+row_limit], preparation_context=preparation_context,
            programs=None if programs is None else programs[i:i+row_limit])
            for i in range(0,len(energy_keys),row_limit)]
        return tuple(np.concatenate([result[k] for result in results]) for k in range(2))
    states = np.stack([np.asarray(state)[0] for state in source_states])
    bound = ([_bind_projection_program(program, tuple(state), model_size, binding)
              for program, state, binding in zip(programs, states, bindings)] if programs is not None else
             [_bound_projection_row(key,tuple(state),model_size,binding)
              for key,state,binding in zip(energy_keys,states,bindings)])
    deltas,kinds,gaps,differences,counts = [np.stack(items) for items in zip(*bound)]
    permutations = np.stack(permutations)
    linear = (np.einsum('rpm,rfm->rpf',frequency_blocks,deltas,optimize=False)
              if preparation_context is None else
              preparation_context.energy_contractions(frequency_blocks,deltas))
    row_count,factor_count,model_count = differences.shape[:3]
    local = differences[
        np.arange(row_count)[:,None,None,None,None],
        np.arange(factor_count)[None,:,None,None,None],
        np.arange(model_count)[None,None,:,None,None],
        permutations[:,None,None,:,:]
    ]
    matches = (np.all(local == deltas[:,:,None,None,:],axis=-1)
               & (np.count_nonzero(local,axis=-1) == counts[:,:,:,None]))
    members = np.any(matches,axis=2)
    occupation = np.take_along_axis(states[:,None,:],permutations,axis=2)
    possible = np.all(occupation[:,None,:,:]+deltas[:,:,None,:] >= 0,axis=-1)
    factor_masks = np.where(kinds[:,:,None] == 0,True,
        possible & np.where(kinds[:,:,None] == 2,members,~members))
    good = np.all(factor_masks,axis=1)
    weights = np.where((kinds[:,None,:] == 1) | (kinds[:,None,:] == 5),gaps[:,None,:]-linear,
                       np.where(kinds[:,None,:] == 0,linear,1.))
    return good,np.prod(weights,axis=-1)


_runtime_roles = itertools.count(1)


class AnalyticBlockExpansion:
    """Generate canonical or intermediate-normalized block expansions.

    In the default canonical convention, the following Y/N/F definitions
    apply. With intermediate_normalization=True, Y is Omega, N[n>0]=0,
    and F is the generally non-Hermitian Bloch series P H Omega P. Its
    metric and canonical counterpart are computed lazily on request.

    Y_0=P; X_n=Q Y_n P; N_n=P Y_n P; F=Y^dagger H Y.
    Every X_n is built with a projected resolvent. No model-space energy
    difference is ever inverted. N_n and F_n are finite block contractions.

    Y[n][a] is a finite sum over model-space origins b. Its values are native
    analytic operator terms acting on |b>; the complementary oscillator
    states are summed by the existing polynomial/DAG evaluator, not enumerated
    as a Fock-state basis. The model-space overlap PY is chosen Hermitian.
    """
    def __init__(self, evaluator, states, expansions=None, freqs=None, zero_cutoff=None,
                 evaluation_backend=None, intermediate_normalization=False):
        self.intermediate_normalization = bool(intermediate_normalization)
        self._intermediate_expansion = None
        self._canonical_expansion = None
        self._metric = [np.eye(len(states))]
        self.evaluator = evaluator
        registry = getattr(evaluator,'_projection_runtime_inputs',None)
        if registry is None:
            registry = evaluator._projection_runtime_inputs = dict(projectors={},scalars={})
        self.runtime_inputs = registry
        self.states = tuple(tuple(int(n) for n in state) for state in states)
        self.freqs = np.asarray(evaluator.freqs if freqs is None else freqs)
        self._raw_expansions = expansions
        self.expansions = evaluator._prep_expansions(expansions)
        self.zero_cutoff = zero_cutoff
        self.evaluation_backend = evaluation_backend
        self.identity = IdentityTerm(logger=evaluator.solver.logger,
                                     polynomial_representation=evaluator.solver.polynomial_representation)
        self.Y = [[{state: self.identity} for state in self.states]]
        self.X = [[{} for state in self.states]]
        self.N = [np.eye(len(self.states))]
        self.F = [np.diag((np.asarray(self.states) + .5) @ self.freqs)]
        self._elements = {}
        self._coefficient_references = {}
        self._resolvents = {}
        self._inverse_normalizers = [np.eye(len(self.states))]
        self._intermediate_columns = [self.Y[0]]
        self._operator_pieces = {}
        self._generic_terms = {}
        self._hamiltonian_terms = {}

    def _scalar(self,value):
        role=next(_runtime_roles)
        self.runtime_inputs['scalars'][role]=float(value)
        return _BlockScalarTerm.lookup(len(self.expansions),role,
            self.evaluator.solver.logger,self.evaluator.solver.polynomial_representation)

    def _resolvent(self,term,state,origin):
        trivial=False
        if tuple(state)==tuple(origin):
            trivial=not any(
                tuple(other)!=tuple(state) and
                PerturbationTheoryTerm.sorted_changes((np.asarray(other)-origin)[np.asarray(other)!=origin]) in term.changes
                for other in self.states)
        role=next(_runtime_roles)
        self.runtime_inputs['projectors'][role]=(self.states,tuple(float(w) for w in self.freqs),tuple(state),trivial)
        return ProjectedTerm.lookup(term,ProjectionFactor('R',role))

    def _generic_term(self,term,scalar_order):
        key=(term,scalar_order)
        if key not in self._generic_terms:
            projectors,scalars,memo = {},{},{}
            def walk(node):
                if node in memo:return memo[node]
                if isinstance(node,ProjectedTerm):
                    role=node.factor.slot
                    slot=projectors.setdefault(role,len(projectors))
                    result=ProjectedTerm.lookup(walk(node.term),ProjectionFactor(node.factor.kind,slot))
                elif isinstance(node,_BlockScalarTerm):
                    slot=scalars.setdefault(node.slot,len(scalars))
                    result=_BlockScalarTerm.lookup(scalar_order,slot,node.logger,node.polynomial_representation)
                elif isinstance(node,PerturbationTheoryTermProduct):
                    result=PerturbationTheoryTermProduct.lookup(walk(node.gen2),walk(node.gen1))
                elif isinstance(node,PerturbationTheoryTermSum):
                    result=PerturbationTheoryTermSum.lookup(*(walk(t) for t in node.terms))
                elif isinstance(node,ScaledPerturbationTheoryTerm):
                    result=ScaledPerturbationTheoryTerm.lookup(walk(node.base),node.prefactor)
                elif isinstance(node,ShiftedEnergyBaseline):
                    result=ShiftedEnergyBaseline.lookup(walk(node.base))
                elif isinstance(node,ValidResolventTerm):
                    result=ValidResolventTerm.lookup(walk(node.subterm))
                else:
                    result=node
                memo[node]=result
                return result
            generic=walk(term)
            self._generic_terms[key]=(generic,tuple(projectors),tuple(scalars))
        generic,projectors,scalars=self._generic_terms[key]
        bindings=tuple(self.runtime_inputs['projectors'][role] for role in projectors)
        values=[self.runtime_inputs['scalars'][role] for role in scalars]
        return generic,dict(projectors=bindings,scalar_order=scalar_order,scalar_count=len(values)),values

    def _column(self,parts):
        grouped = {}
        for column, coefficient, left in parts:
            if coefficient == 0:
                continue
            # Skip exactly zero finite-block contractions at numerical assembly,
            # before requesting a generic equation. Never embed their values.
            if (isinstance(coefficient,_BlockScalarTerm)
                    and self.runtime_inputs['scalars'][coefficient.slot] == 0):
                continue
            for origin, term in column.items():
                if left is not None:
                    term = left if isinstance(term, IdentityTerm) else left * term
                if isinstance(coefficient,PerturbationTheoryTerm):
                    term = coefficient if isinstance(term,IdentityTerm) else coefficient * term
                elif coefficient != 1:
                    term = coefficient * term
                grouped.setdefault(origin, []).append(term)
        return {origin: terms[0] if len(terms) == 1 else PerturbationTheoryTermSum.lookup(*terms)
                for origin, terms in grouped.items()}

    use_batched_elements = True
    use_hermitian_symmetry = True
    max_element_batch_size = 128
    element_batch_bytes = 64 * 1024**2

    @classmethod
    def _evaluate_elements(cls, requests, logger=None, verbose=False, coefficient_batch_size=1):
        """Use native multi-route evaluation; keep pending requests bounded.

        Each request owns its final-value cache and runtime projector binding.
        Compatible requests share the existing coefficient/permutation and
        polynomial machinery, including requests from different model spaces.
        """
        iterator = iter(requests)
        while True:
            try:
                first = next(iterator)
            except StopIteration:
                return
            owner = first[0]
            subset_batch = max(1, int((getattr(owner.evaluator,'evaluation_options',None) or {}).get('batch_size', 512)))
            # Bound input permutation expansion as well as request count. This
            # is a batching target, not an RSS cap.
            limit = min(cls.max_element_batch_size,
                        max(1, cls.element_batch_bytes // max(1, 16*len(owner.freqs)*subset_batch)))
            if not cls.use_batched_elements:
                limit = 1
            elif limit >= coefficient_batch_size:
                limit = (limit//coefficient_batch_size)*coefficient_batch_size
            chunk = [first] + list(itertools.islice(iterator, limit-1))
            values = [None] * len(chunk)
            pending, groups = {}, {}
            for position, (block, tag, term, initial, final, expansions) in enumerate(chunk):
                if isinstance(term, IdentityTerm):
                    values[position] = float(tuple(initial) == tuple(final))
                    continue
                raw_expansions = expansions
                expansions = block.expansions if expansions is None else expansions
                block._coefficient_references[id(expansions)] = expansions
                key = (term, tuple(initial), tuple(final), id(expansions))
                if key in block._elements:
                    values[position] = block._elements[key]
                    continue
                identity = (block, key)
                if identity in pending:
                    pending[identity]['positions'].append(position)
                    continue
                entry = dict(block=block, key=key, positions=[position], value=0.)
                pending[identity] = entry
                log = Logger.lookup(block.evaluator.solver.logger if logger is None else logger)
                level = Logger.LogLevel.Normal if verbose else Logger.LogLevel.Debug
                changes = block.evaluator.get_diff_map([[initial, [final]]], logger=log, log_level=level)
                generic, bindings, scalars = block._generic_term(term, len(expansions))
                domain = ((id(block.evaluator), id(block._raw_expansions)) if raw_expansions is None
                          else (None, id(expansions)))
                for change, routes in changes.items():
                    group_key = (domain, id(block.freqs), generic, tuple(change), tuple(scalars),
                                 id(block.zero_cutoff), block.evaluation_backend, id(block.evaluator.parallelizer))
                    group = groups.setdefault(group_key, dict(block=block, generic=generic, change=change,
                        expansions=expansions, scalars=scalars, bindings=bindings, routes=[],
                        projectors=[], destinations=[], logger=log))
                    start = len(group['routes'])
                    group['routes'].extend(routes)
                    group['projectors'].extend([bindings['projectors']] * len(routes))
                    group['destinations'].append((entry, start))
            compatible = {}
            for group_key,group in groups.items():
                # Equality compares exact source/permutation routes and runtime
                # projectors; no Cartesian expansion of unmatched domains.
                from .Analytic import _HashableNumPyIndex
                route_key = tuple((_HashableNumPyIndex(state),_HashableNumPyIndex(perms))
                                  for state,perms in group['routes'])
                key = (group_key[1:],route_key,tuple(group['projectors']))
                compatible.setdefault(key,[]).append(group)
            from .Evaluators import current_compilation
            cohorts = compatible.values()
            # The legacy evaluator's historical multi-expansion cutoff policy
            # is coupled. Keep independent property domains separate there.
            if current_compilation() is None or any(
                    item['block'].evaluation_backend == 'legacy' for items in cohorts for item in items):
                cohorts = ([item] for items in compatible.values() for item in items)
            for cohort in cohorts:
                group = cohort[0]
                block = group['block']
                bindings = dict(group['bindings'], route_projectors=tuple(group['projectors']))
                coefficients = [list(item['expansions']) + [[value] for value in item['scalars']]
                                for item in cohort]
                generic, change = group['generic'], group['change']
                if current_compilation() is not None and current_compilation()[0].using_artifacts:
                    expression = generic(change)
                    expression.logger = group['logger']
                else:
                    expression = PerturbationTheoryExpressionEvaluator(
                        generic, generic.get_poly_terms(change), change, logger=group['logger'])
                subset_batch=max(1,int((getattr(block.evaluator,'evaluation_options',None) or {}).get('batch_size',512)))
                permutation_bytes=16*len(block.freqs)*subset_batch*len(cohort)*max(len(perms) for _,perms in group['routes'])
                # Account conservatively for the bound model-space/reference
                # metadata used while fingerprinting a numerical batch.
                projection_bytes=64*len(block.freqs)*max(
                    sum(len(binding[0])+2 for binding in route) for route in group['projectors'])
                route_limit=max(1,min(cls.max_element_batch_size,
                    cls.element_batch_bytes//max(1,permutation_bytes+projection_bytes)))
                if not cls.use_batched_elements:
                    route_limit=1
                for start in range(0,len(group['routes']),route_limit):
                    stop=min(start+route_limit,len(group['routes']))
                    route_bindings=dict(bindings,route_projectors=tuple(group['projectors'][start:stop]),
                                        independent_coefficients=len(cohort)>1)
                    with _log_evaluation_time(group['logger'], "evaluating projected element batch...",
                                              log_level=Logger.LogLevel.Normal if verbose else Logger.LogLevel.Debug):
                        result = expression.evaluate(group['routes'][start:stop], coefficients, block.freqs,
                            projection_bindings=route_bindings, only_degenerate_terms=False,
                            zero_cutoff=block.zero_cutoff, parallelizer=block.evaluator.parallelizer,
                            evaluation_backend=block.evaluation_backend, log_scaled=False, verbose=verbose)
                    for axis,item in enumerate(cohort):
                        for entry, route in item['destinations']:
                            if start<=route<stop:
                                entry['value'] += np.asarray(result[route-start][axis]).item()
            for entry in pending.values():
                entry['block']._elements[entry['key']] = entry['value']
                for position in entry['positions']:
                    values[position] = entry['value']
            for request, value in zip(chunk, values):
                yield request[1], value

    def element(self, term, initial, final, expansions=None, logger=None, verbose=False):
        return next(self._evaluate_elements(
            [(self, None, term, initial, final, expansions)], logger, verbose))[1]

    def _inner_requests(self, left, right, operator=None):
        for initial, right_term in right.items():
            if operator is not None:
                right_term = operator if isinstance(right_term, IdentityTerm) else operator * right_term
            for final, left_term in left.items():
                if isinstance(left_term, IdentityTerm):
                    expression = right_term
                elif isinstance(right_term, IdentityTerm):
                    expression = ShiftedEnergyBaseline.lookup(left_term)
                else:
                    expression = ShiftedEnergyBaseline.lookup(left_term) * right_term
                yield expression, initial, final

    def inner(self, left, right, operator=None, expansions=None, logger=None, verbose=False):
        requests = ((self, None, term, initial, final, expansions)
                    for term, initial, final in self._inner_requests(left, right, operator))
        result = 0.
        for _, value in self._evaluate_elements(requests, logger, verbose):
            result += value
        return result

    @classmethod
    def _contract_matrices(cls, specifications, logger=None, verbose=False):
        """Contract selected entries, retaining the original inner-sum order."""
        matrices = [np.zeros((len(left), len(right))) for block,left,right,op,upper in specifications]
        def requests():
            for index, (block, left, right, operator, upper) in enumerate(specifications):
                for a, column_left in enumerate(left):
                    for b in range(a if upper else 0, len(right)):
                        for term, initial, final in block._inner_requests(column_left, right[b], operator):
                            yield block, (index,a,b), term, initial, final, None
        for tag, entries in itertools.groupby(cls._evaluate_elements(requests(), logger, verbose),
                                               key=lambda item:item[0]):
            value = 0.
            for _, part in entries:
                value += part
            index,a,b = tag
            matrices[index][a,b] = value
        return matrices

    @staticmethod
    def _complete_hermitian(matrix):
        upper = np.triu(matrix)
        return upper + np.triu(matrix,1).conj().T

    def _hamiltonian(self, order):
        if order not in self._hamiltonian_terms:
            if order >= len(self.expansions) or not any(np.any(t) for t in self.expansions[order]):
                term = None
            else:
                term = self.evaluator.solver.hamiltonian_expansion[order]
            self._hamiltonian_terms[order] = term
        return self._hamiltonian_terms[order]

    def extend(self, order, logger=None, verbose=False):
        if self.intermediate_normalization:
            return self._extend_intermediate(order,logger=logger,verbose=verbose)
        self._extend_blocks([self], order, logger, verbose)
        return self

    @classmethod
    def _extend_blocks(cls, blocks, order, logger=None, verbose=False):
        """Advance independent canonical blocks together, one order at a time."""
        blocks = list(dict.fromkeys(blocks))
        ordinary = [block for block in blocks if not block.intermediate_normalization]
        for block in blocks:
            if block.intermediate_normalization:
                block._extend_intermediate(order,logger=logger,verbose=verbose)
        upper = cls.use_hermitian_symmetry
        while True:
            active = [block for block in ordinary if len(block.Y) <= order]
            if not active:
                return
            n = min(len(block.Y) for block in active)
            active = [block for block in active if len(block.Y) == n]
            log = Logger.lookup(active[0].evaluator.solver.logger if logger is None else logger)
            with _log_evaluation_time(log, "Building projected blocks at order {o}", o=n):
                specifications, locations = [], []
                for block in active:
                    for i in range(1,n):
                        locations.append((block,i))
                        specifications.append((block,block.X[i],block.X[n-i],None,upper))
                with _log_evaluation_time(log, "Contracting block normalization..."):
                    pieces = dict(zip(locations,cls._contract_matrices(specifications,log,verbose)))
                    normalizers = {}
                    for block in active:
                        size = len(block.states)
                        normalizer = np.zeros((size,size))
                        for a in range(size):
                            for b in range(a if upper else 0,size):
                                normalizer[a,b] = -.5 * sum(pieces[block,i][a,b] for i in range(1,n))
                        for i in range(1,n):
                            normalizer -= .5 * (block.N[i].conj().T @ block.N[n-i])
                        normalizers[block] = cls._complete_hermitian(normalizer) if upper else normalizer
                with _log_evaluation_time(log, "Contracting block Hamiltonian..."):
                    specifications, locations = [], []
                    identities = {block:[{state:block.identity} for state in block.states] for block in active}
                    for block in active:
                        for r in range(1,n+1):
                            h = block._hamiltonian(r)
                            if h is not None:
                                locations.append((block,r,'bare'))
                                specifications.append((block,identities[block],identities[block],h,upper))
                                locations.append((block,r,'mixed'))
                                specifications.append((block,identities[block],block.X[n-r],h,upper))
                    pieces = dict(zip(locations,cls._contract_matrices(specifications,log,verbose)))
                    hamiltonians = {}
                    for block in active:
                        energies = np.diag(block.F[0])
                        hamiltonian = (energies[:,None]-energies[None,:]) * normalizers[block]
                        for r in range(1,n+1):
                            if block._hamiltonian(r) is not None:
                                bare = pieces[block,r,'bare']
                                if upper:
                                    bare = cls._complete_hermitian(bare)
                                hamiltonian += bare @ block.N[n-r]
                                hamiltonian += pieces[block,r,'mixed']
                        for i in range(1,n):
                            hamiltonian -= block.N[i] @ block.F[n-i]
                        hamiltonians[block] = cls._complete_hermitian(hamiltonian) if upper else hamiltonian
                with _log_evaluation_time(log, "Building projected resolvents and wavefunctions..."):
                    for block in active:
                        size = len(block.states)
                        projected, full = [], []
                        for a,state in enumerate(block.states):
                            parts = [(block.X[n-r][a],1,block._hamiltonian(r))
                                     for r in range(1,n+1) if block._hamiltonian(r) is not None]
                            parts += [({origin:block.identity},block._scalar(block.N[n-r][b,a]),block._hamiltonian(r))
                                      for r in range(1,n+1) if block._hamiltonian(r) is not None
                                      for b,origin in enumerate(block.states)]
                            parts += [(block.X[i][b],block._scalar(-block.F[n-i][b,a]),None)
                                      for i in range(1,n) for b in range(size)]
                            column = {}
                            for origin,term in block._column(parts).items():
                                key = (term,state,origin)
                                if key not in block._resolvents:
                                    block._resolvents[key] = block._resolvent(term,state,origin)
                                column[origin] = block._resolvents[key]
                            projected.append(column)
                            normalization = [({origin:block.identity},block._scalar(normalizers[block][b,a]),None)
                                             for b,origin in enumerate(block.states)]
                            full.append(block._column([(column,1,None)]+normalization))
                        block.N.append(normalizers[block]);block.F.append(hamiltonians[block])
                        block.X.append(projected);block.Y.append(full)

    def _extend_intermediate(self, order, logger=None, verbose=False):
        """Bloch recurrence: P Omega P=P, B=P H Omega P.

        F holds B in this explicitly selected representation; it is generally
        non-Hermitian. No overlap contraction is needed to build it.
        """
        logger=Logger.lookup(self.evaluator.solver.logger if logger is None else logger)
        size=len(self.states)
        while len(self.Y)<=order:
            n=len(self.Y)
            with _log_evaluation_time(logger,"Building intermediate-normalized block at order {o}",o=n):
                hamiltonian=np.zeros((size,size))
                for r in range(1,n+1):
                    h=self._hamiltonian(r)
                    if h is not None:
                        for a,state in enumerate(self.states):
                            for b in range(size):
                                hamiltonian[a,b]+=self.inner(
                                    {state:self.identity},self.Y[n-r][b],h,
                                    logger=logger,verbose=verbose)
                projected=[]
                for a,state in enumerate(self.states):
                    parts=[(self.Y[n-r][a],1,self._hamiltonian(r))
                           for r in range(1,n+1) if self._hamiltonian(r) is not None]
                    parts += [(self.X[i][b],self._scalar(-self.F[n-i][b,a]),None)
                              for i in range(1,n) for b in range(size)]
                    column={}
                    for origin,term in self._column(parts).items():
                        key=(term,state,origin)
                        if key not in self._resolvents:
                            self._resolvents[key]=self._resolvent(term,state,origin)
                        column[origin]=self._resolvents[key]
                    projected.append(column)
                self.F.append(hamiltonian)
                self.N.append(np.zeros((size,size)))
                self.X.append(projected)
                self.Y.append(projected)
        return self

    def _intermediate(self):
        if self.intermediate_normalization:
            return self
        if self._intermediate_expansion is None:
            self._intermediate_expansion=type(self)(self.evaluator,self.states,
                self._raw_expansions,self.freqs,self.zero_cutoff,self.evaluation_backend,
                intermediate_normalization=True)
        return self._intermediate_expansion

    def intermediate_columns(self, order, logger=None, verbose=False):
        """Build Omega directly, without first building normalized columns."""
        block=self._intermediate().extend(order,logger=logger,verbose=verbose)
        return block.Y[order]

    def overlap_metric(self, order, logger=None, verbose=False):
        """Cache S=Omega^dagger Omega once per block and perturbation order."""
        block=self._intermediate().extend(order,logger=logger,verbose=verbose)
        size=len(block.states)
        while len(block._metric)<=order:
            n=len(block._metric)
            with _log_evaluation_time(Logger.lookup(block.evaluator.solver.logger if logger is None else logger),
                                      "Contracting deferred block metric at order {o}",o=n):
                metric=np.zeros((size,size))
                for a in range(size):
                    for b in range(size):
                        metric[a,b]=sum(block.inner(block.X[i][a],block.X[n-i][b],
                            logger=logger,verbose=verbose) for i in range(1,n))
                block._metric.append(metric)
        return block._metric[:order+1]

    def _metric_factors(self, order, logger=None, verbose=False):
        metric=self.overlap_metric(order,logger=logger,verbose=verbose)
        root=[np.eye(len(self.states))];inverse=[root[0]]
        for n in range(1,order+1):
            root.append(.5*(metric[n]-sum((root[i]@root[n-i] for i in range(1,n)),
                                          np.zeros_like(root[0]))))
            inverse.append(-sum((root[i]@inverse[n-i] for i in range(1,n+1)),
                                np.zeros_like(root[0])))
        return root,inverse

    def orthonormal_hamiltonian(self, order, logger=None, verbose=False):
        """Return S^(1/2) B S^(-1/2), consistently truncated at each order."""
        if not self.intermediate_normalization:
            return self.extend(order,logger=logger,verbose=verbose).F[:order+1]
        root,inverse=self._metric_factors(order,logger=logger,verbose=verbose)
        return [sum((root[i]@self.F[j]@inverse[n-i-j]
                     for i in range(n+1) for j in range(n+1-i)),
                    np.zeros_like(root[0])) for n in range(order+1)]

    def orthonormal_block(self, order, logger=None, verbose=False):
        """Reuse Omega and its cached metric for normalized observables.

        Y=Omega S^(-1/2). No second Q-space wavefunction recurrence is run.
        The returned block preserves the ordinary bra-H-ket contribution API.
        """
        if not self.intermediate_normalization:
            return self.extend(order,logger=logger,verbose=verbose)
        root,inverse=self._metric_factors(order,logger=logger,verbose=verbose)
        hamiltonians=self.orthonormal_hamiltonian(order,logger=logger,verbose=verbose)
        if self._canonical_expansion is None:
            self._canonical_expansion=type(self)(self.evaluator,self.states,
                self._raw_expansions,self.freqs,self.zero_cutoff,self.evaluation_backend)
            self._canonical_expansion._intermediate_expansion=self
        block=self._canonical_expansion
        while len(block.Y)<=order:
            n=len(block.Y);projected=[];full=[]
            for a,state in enumerate(self.states):
                column=self._column([(self.X[i][b],self._scalar(inverse[n-i][b,a]),None)
                    for i in range(1,n+1) for b in range(len(self.states))
                    if inverse[n-i][b,a]!=0])
                projected.append(column)
                full.append(self._column([(column,1,None)]+[
                    ({origin:self.identity},self._scalar(inverse[n][b,a]),None)
                    for b,origin in enumerate(self.states) if inverse[n][b,a]!=0]))
            block.N.append(inverse[n]);block.F.append(hamiltonians[n])
            block.X.append(projected);block.Y.append(full)
        return block

    def operator_components(self, order, column, off_diagonal=None):
        # Keep finite P coefficients outside the free-space operator products.
        # Inlining a numerical N entry into a polynomial sum would specialize
        # and duplicate the same Q-space DAG for every reference state.
        result = []
        if off_diagonal is not True and self.X[order][column]:
            result.append((self.X[order][column],1.))
        result += [({state:self.identity},self.N[order][b,column])
                   for b,state in enumerate(self.states)
                   if self.N[order][b,column] != 0 and
                   (off_diagonal is None or (b != column) == off_diagonal)]
        return result

    def _operator_inner_requests(self, other, left, right, operators, order,
                                 terms=None, branches=True, logger=None, verbose=False):
        if self.intermediate_normalization or other.intermediate_normalization:
            raise ValueError('normalized observables require orthonormal blocks or explicit metric factors')
        self.extend(order,logger=logger,verbose=verbose)
        other.extend(order,logger=logger,verbose=verbose)
        allowed = None if terms is None else {tuple(t) for t in terms}
        for i in range(order+1):
            for k in range(order+1-i):
                j = order-i-k
                if j < len(operators) and (allowed is None or (k, j, i) in allowed):
                    if branches is True:
                        for bra,bc in self.operator_components(i,left):
                            for ket,kc in other.operator_components(k,right):
                                yield bc*kc, bra, ket, operators[j]
                    else:
                        selected = {'main': True}
                        if branches is not False:
                            selected.update(branches if isinstance(branches,dict)
                                            else {key:True for key in branches})
                        for name,bra_cross,ket_cross in [('main',False,False),
                                ('left',True,False),('right',False,True),('both',True,True)]:
                            if selected.get(name,False):
                                for bra,bc in self.operator_components(i,left,bra_cross):
                                    for ket,kc in other.operator_components(k,right,ket_cross):
                                        yield bc*kc, bra, ket, operators[j]

    def operator_correction(self, other, left, right, operators, order, expansions=None,
                            terms=None, branches=True, logger=None, verbose=False):
        result = 0.
        for weight, bra, ket, operator in self._operator_inner_requests(
                other, left, right, operators, order, terms, branches, logger, verbose):
            result += weight*self.inner(bra, ket, operator, expansions, logger, verbose)
        return result


class ProjectedCorrectionEvaluator:
    """Resolve actual blocks and expose the ordinary analytic result layout."""
    def __init__(self, evaluator, blocks=None, pairs=None, expansions=None,
                 freqs=None, zero_cutoff=None, logger=None, verbose=False):
        self.evaluator = evaluator
        self.logger = Logger.lookup(evaluator.solver.logger if logger is None else logger)
        self.verbose = verbose
        self.expansions = expansions
        self.freqs = np.asarray(evaluator.freqs if freqs is None else freqs)
        self.zero_cutoff = zero_cutoff
        self.blocks = self.resolve_blocks(blocks, pairs)
        _projected_evaluators.add(evaluator)
        context_key = (coefficient_fingerprint(evaluator.expansions if expansions is None else expansions),
                       tuple(self.freqs), zero_cutoff)
        cache = getattr(evaluator, '_projected_block_cache', None)
        if cache is None or cache[0] != context_key:
            cache = evaluator._projected_block_cache = (context_key, {})
        self.cache = cache[1]

    @staticmethod
    def resolve_blocks(blocks=None, pairs=None):
        """Merge exact runtime relations without making a global change set."""
        groups = [] if blocks is None else [list(group) for group in blocks]
        if blocks is None and pairs is not None:
            if isinstance(pairs, dict):
                for state, neighbors in pairs.items():
                    groups.append([state])
                    groups.extend(list(pair) for pair in neighbors)
            else:
                groups.extend(list(pair) for pair in pairs)
        components = []
        for group in groups:
            group = {tuple(int(n) for n in state) for state in group}
            touching = [part for part in components if part & group]
            for part in touching:
                group |= part
                components.remove(part)
            components.append(group)
        return {state: tuple(sorted(group)) for group in components for state in group}

    def block(self, state, order, intermediate=False):
        state = tuple(int(n) for n in state)
        group = self.blocks.get(state, (state,))
        key=(group,'intermediate') if intermediate else group
        if not intermediate and key not in self.cache and (group,'intermediate') in self.cache:
            self.cache[key]=self.cache[(group,'intermediate')].orthonormal_block(
                order,logger=self.logger,verbose=self.verbose)
        if key not in self.cache:
            self.cache[key] = AnalyticBlockExpansion(self.evaluator, group,
                self.expansions, self.freqs, self.zero_cutoff,
                intermediate_normalization=intermediate)
        block = self.cache[key].extend(order,logger=self.logger,verbose=self.verbose)
        return block, block.states.index(state)

    def _prepare_blocks(self, states, order):
        blocks = []
        for state in states:
            state = tuple(int(n) for n in state)
            group = self.blocks.get(state,(state,))
            if group not in self.cache and (group,'intermediate') in self.cache:
                # Preserve the existing lazy conversion of a previously built
                # intermediate expansion, including its metric policy.
                self.block(state,order)
            elif group not in self.cache:
                self.cache[group] = AnalyticBlockExpansion(self.evaluator,group,
                    self.expansions,self.freqs,self.zero_cutoff)
            blocks.append(self.cache[group])
        AnalyticBlockExpansion._extend_blocks(blocks,order,self.logger,self.verbose)

    def diagonal(self, states, order, series, orders):
        self._prepare_blocks(states, order)
        corrections = []
        for n in orders:
            with self.logger.block(tag="Getting corrections at order {o}", o=n):
                with _log_evaluation_time(self.logger, "evaluating..."):
                    values = []
                    for state in states:
                        block, a = self.block(state, order)
                        values.append(getattr(block, series)[n][a,a])
                    corrections.append(np.asarray(values))
        return corrections

    def result(self, states, order, value, target_orders=None, *, prepare_canonical=True, prepare_final_blocks=False, batch_value=None):
        if target_orders is not None:
            target_orders = frozenset(target_orders)
            if not target_orders or min(target_orders) < 0 or max(target_orders) > order:
                raise ValueError('target_orders must be nonempty orders between 0 and order')
        if np.ndim(states[0][0]) == 0:
            states = [[[0] * len(self.freqs), states]]
        grouped = {}
        for initial, finals in states:
            key = tuple(int(n) for n in initial)
            row = grouped.setdefault(key, {})
            for final in finals:
                row.setdefault(tuple(int(n) for n in final), None)
        if prepare_canonical and not self.evaluator.solver.intermediate_normalization:
            needed = itertools.chain(grouped,*(finals for finals in grouped.values())) if prepare_final_blocks else grouped
            self._prepare_blocks(needed,order)
        batched_values = None if batch_value is None else iter(batch_value(
            grouped, [n for n in range(order+1) if target_orders is None or n in target_orders]))
        basis = HarmonicOscillatorProductBasis(len(self.freqs))
        initial_space = BasisStateSpace(basis, list(grouped))
        final_spaces, corrections = [], []
        for initial, finals in grouped.items():
            final_spaces.append(BasisStateSpace(basis, list(finals)))
            row = []
            for n in range(order+1):
                if target_orders is not None and n not in target_orders:
                    row.append(None)
                    continue
                with self.logger.block(tag="Getting corrections at order {o}", o=n):
                    with _log_evaluation_time(self.logger, "evaluating..."):
                        row.append(np.asarray([value(initial, final, n) if batched_values is None
                                               else next(batched_values) for final in finals]))
            corrections.append(row)
        return BasicAPTCorrections(initial_space, final_spaces, corrections)

    def wavefunctions(self, states, order, full=True, target_orders=None):
        targets = None if target_orders is None else set(target_orders)
        def value(initial, final, n):
            if targets is not None and n not in targets:
                return 0.
            block, a = self.block(initial, order,
                                  intermediate=bool(self.evaluator.solver.intermediate_normalization))
            column = block.Y[n][a] if full else block.X[n][a]
            if self.evaluator.solver.intermediate_normalization:
                column = block.intermediate_columns(n,logger=self.logger,verbose=self.verbose)[a]
                if not full and n == 0:
                    column = {}
            return -sum(block.element(term, origin, final,logger=self.logger,verbose=self.verbose) for origin, term in column.items())
        return self.result(states, order, value, target_orders=targets)

    @staticmethod
    def _operator_branches(branches):
        # Public correction rows are initial kets -> final bras. The legacy
        # "left" branch names the initial-state normalization, whereas the
        # matrix contraction's left factor is the final bra.
        if branches is True or branches is False:
            return branches
        if isinstance(branches,dict):
            return {('right' if key=='left' else 'left' if key=='right' else key):value
                    for key,value in branches.items()}
        return [('right' if key=='left' else 'left' if key=='right' else key) for key in branches]

    def hamiltonian(self, states, order, only_degenerate_terms=True, terms=None, branches=True,
                    target_orders=None, normalization='orthonormal'):
        if normalization not in {'orthonormal','intermediate'}:
            raise ValueError('unknown block normalization '+str(normalization))
        intermediate=normalization=='intermediate'
        if intermediate and (terms is not None or branches is not True):
            raise ValueError('Bloch entries are not normalized bra-H-ket contribution triples')
        def value(initial, final, n):
            block, a = self.block(initial, order, intermediate=intermediate)
            final = tuple(final)
            if final not in block.states:
                raise ValueError('Hamiltonian block endpoints must belong to the same model space')
            b = block.states.index(final)
            result = (block.F[n][b,a] if terms is None and branches is True else
                      block.operator_correction(block,b,a,self.evaluator.solver.hamiltonian_expansion,
                                                n,terms=terms,branches=self._operator_branches(branches),
                                                logger=self.logger,verbose=self.verbose))
            return result if only_degenerate_terms else -result
        return self.result(states, order, value, target_orders=target_orders, prepare_canonical=not intermediate)

    def operators(self, states, order, operators, expansions, terms=None, branches=True,
                  target_orders=None):
        return_single = self.evaluator.is_single_expansion(expansions)
        batches = [expansions] if return_single else expansions
        branches = self._operator_branches(branches)
        def values(grouped, orders):
            def requests():
                cell = 0
                for initial, finals in grouped.items():
                    right, ket = self.block(initial, order)
                    for n in orders:
                        for final in finals:
                            left, bra = self.block(final, order)
                            # A zero marker keeps empty cells in the output
                            # stream without materializing all transition work.
                            for axis, coefficients in enumerate(batches):
                                yield left, (cell,-1,axis,0.), left.identity, initial, initial, coefficients
                            for component,(weight,lhs,rhs,operator) in enumerate(left._operator_inner_requests(
                                    right,bra,ket,operators,n,terms,branches,self.logger,self.verbose)):
                                for term,origin,target in left._inner_requests(lhs,rhs,operator):
                                    for axis,coefficients in enumerate(batches):
                                        yield left,(cell,component,axis,weight),term,origin,target,coefficients
                            cell += 1
            current_cell = current_component = None
            total = np.zeros(len(batches))
            inner = np.zeros(len(batches))
            weight = 0.
            for (cell,component,axis,coefficient),value in AnalyticBlockExpansion._evaluate_elements(
                    requests(),self.logger,self.verbose,coefficient_batch_size=len(batches)):
                if (cell,component) != (current_cell,current_component):
                    if current_cell is not None:
                        total += weight*inner
                        if cell != current_cell:
                            yield -total
                            total = np.zeros(len(batches))
                    current_cell,current_component = cell,component
                    inner = np.zeros(len(batches))
                    weight = coefficient
                inner[axis] += value
            if current_cell is not None:
                total += weight*inner
                yield -total
        combined = self.result(states,order,None,target_orders=target_orders,
                               prepare_final_blocks=True,batch_value=values)
        results = [BasicAPTCorrections(combined.initial_states,combined.final_states,
                   [[None if row is None else row[:,axis] for row in orders]
                    for orders in combined.corrections])
                   for axis in range(len(batches))]
        return results[0] if return_single else results


# ---------------------------------------------------------------------------------------------
#  Conditional (formerly VPT2/Conditional.py)
#  State-independent conditional operator algebra and ordered VPT transforms.
#
#  Conditions are attached to operator endpoints, not resolved model states.
#  All sums and products are native analytic perturbation-theory terms.
# ---------------------------------------------------------------------------------------------

_stacks = weakref.WeakSet()


def _clear_conditional_caches():
    ConditionalOperatorTerm._cache.clear()
    ConditionEndpoint._cache.clear()
    operator_sector.cache_clear()
    canonical_condition_key.cache_clear()
    _resonance_block_labels.cache_clear()
    _filtered_block_labels.cache_clear()
    IntegerLatticeCondition.may_retain_change.cache_clear()
    for stack in list(_stacks):
        cache = getattr(stack.solver,'_conditional_stacks',{})
        for key,value in list(cache.items()):
            if value is stack:
                del cache[key]
        stack._transformed_operators.clear()


def _local_changes(permutations, change, modes):
    permutations, change = np.asarray(permutations), np.asarray(change)
    return np.sum((permutations[...,None] == np.asarray(modes))
                  * change[...,None,:,None], axis=-2)


def _outside_unchanged(permutations, change, modes):
    return np.all((np.asarray(change)[...,None,:] == 0)
                  | np.isin(permutations,modes), axis=-1)


def _function_value(function, points):
    # McUtils elementary functions act on a single coordinate; multivariate
    # functions perform their own index selection.
    if isinstance(function,ElementaryFunction) or getattr(function, 'elementary', False):
        if points.shape[-1] == 0:
            return function(np.zeros(points.shape[:-1]))
        return function(points[...,0 if function.idx is None else function.idx])
    return function(points)


class OperatorCondition:
    """A serializable relation on the two occupations of an operator step."""
    symmetric = False
    diagonal = None
    equivalence = False
    requires_frequencies = False

    def __setattr__(self,name,value):
        if 'key' in self.__dict__ or getattr(type(self),'key',None) is not None:
            raise AttributeError('operator conditions are immutable; construct a new condition')
        object.__setattr__(self,name,value)

    def __hash__(self):
        return hash(self.key)

    def __fingerprint__(self):
        # The same immutable definition drives equality, hashing and persistent
        # compilation. Prepared evaluators are runtime state of the functions.
        return ('operator-condition', self.key)

    def __eq__(self, other):
        if self is other:
            return True
        return isinstance(other,OperatorCondition) and self.key == other.key

    def __and__(self, other):
        return BooleanCondition('and',self,other)

    def __or__(self, other):
        return BooleanCondition('or',self,other)

    def __invert__(self):
        return BooleanCondition('not',self)

    def symmetrized(self):
        return SymmetricCondition(self)

    def evaluate_shifted(self, state, permutations, start, end, frequencies):
        raise NotImplementedError

    def evaluate_bound(self, state, permutations, start, end, frequencies, bindings=()):
        return self.evaluate_shifted(state, permutations, start, end, frequencies)

    def matches(self,initial,final,frequencies=None):
        initial,final = np.broadcast_arrays(np.asarray(initial),np.asarray(final))
        size = initial.shape[-1]
        if frequencies is None:
            if self.requires_frequencies:
                raise ValueError('this condition requires runtime frequencies')
            frequencies = np.zeros(size)
        permutations = np.broadcast_to(np.arange(size),initial.shape[:-1]+(1,size))
        values = self.evaluate_shifted(initial,permutations,np.zeros_like(initial),
            final-initial,np.broadcast_to(frequencies,initial.shape)[...,None,:])
        return values[...,0]

    __call__ = matches

    def may_retain_change(self,change):
        """False proves a change impossible; None leaves it unresolved."""
        return None


def _condition_modes(modes):
    modes = tuple(modes)
    result = tuple(int(mode) for mode in modes)
    if any(mode < 0 or mode != original for mode,original in zip(result,modes)):
        raise ValueError('condition modes must be nonnegative integer indices')
    if len(set(result)) != len(result):
        raise ValueError('condition mode indices must be distinct')
    return result


class SameStateCondition(OperatorCondition):
    symmetric = diagonal = True
    equivalence = True
    key = ('same-state',)

    def evaluate_shifted(self,state,permutations,start,end,frequencies):
        return np.broadcast_to(np.all(np.asarray(start) == end,axis=-1)[...,None],
                               np.shape(permutations)[:-1])

    def __repr__(self):
        return 'n=k'

    def may_retain_change(self,change):
        return not change


@functools.lru_cache(maxsize=16)
def _resonance_block_labels(groups):
    return {state: index for index, group in enumerate(groups) for state in group}


class ResonanceBlockCondition(OperatorCondition):
    """Formal retained-sector mask, bound to exact relations only at evaluation.

    No occupations, energies, coefficients or frequencies enter this symbolic
    object. The existing numerical evaluation context supplies mask membership.
    """
    symmetric = diagonal = equivalence = True
    key = ('runtime-resonance-blocks', 1)

    def evaluate_bound(self, state, permutations, start, end, frequencies, bindings=()):
        permutations = np.asarray(permutations, dtype=int)
        shape = permutations.shape[:-1]
        if np.array_equal(start, end):
            return np.ones(shape, dtype=bool)
        if not bindings:
            raise ValueError('symbolic block recurrence requires runtime resonance bindings')
        labels = _resonance_block_labels(bindings[0])
        permutations = permutations.reshape((-1, permutations.shape[-1]))
        left = np.broadcast_to(state, (len(permutations), len(state))).copy()
        right = left.copy()
        rows = np.arange(len(permutations))[:, None]
        left[rows, permutations] += start
        right[rows, permutations] += end
        return np.fromiter((
            (key in labels and labels[key] == labels.get(tuple(target)))
            for key, target in zip(map(tuple, left), right)
        ), dtype=bool, count=len(left)).reshape(shape)

    def may_retain_change(self, change):
        return True if not change else None

    def __repr__(self):
        return 'RuntimeResonanceBlocks'


class BooleanCondition(OperatorCondition):
    def __init__(self, operation, *conditions):
        self.operation, self.conditions = operation, tuple(conditions)
        self.symmetric = all(condition.symmetric for condition in conditions)
        self.equivalence = operation == 'and' and all(c.equivalence for c in conditions)
        self.requires_frequencies = any(c.requires_frequencies for c in conditions)
        diagonal = [condition.diagonal for condition in conditions]
        if operation == 'and':
            self.diagonal = False if False in diagonal else True if all(v is True for v in diagonal) else None
        elif operation == 'or':
            self.diagonal = True if True in diagonal else False if all(v is False for v in diagonal) else None
        elif operation == 'not':
            self.diagonal = None if diagonal[0] is None else not diagonal[0]
        else:
            raise ValueError('unknown condition composition: '+operation)
        self.key = (operation,tuple(condition.key for condition in conditions))

    def evaluate_shifted(self,*args):
        values = [condition.evaluate_shifted(*args) for condition in self.conditions]
        if self.operation == 'not':
            return ~values[0]
        return (np.logical_and.reduce(values) if self.operation == 'and'
                else np.logical_or.reduce(values))

    def evaluate_bound(self, state, permutations, start, end, frequencies, bindings=()):
        values = [condition.evaluate_bound(state, permutations, start, end, frequencies, bindings)
                  for condition in self.conditions]
        if self.operation == 'not':
            return ~values[0]
        return np.logical_and.reduce(values) if self.operation == 'and' else np.logical_or.reduce(values)

    def __repr__(self):
        return '{}({})'.format(self.operation,','.join(map(repr,self.conditions)))


class SymmetricCondition(OperatorCondition):
    symmetric = True

    def __init__(self,condition):
        self.condition, self.diagonal = condition, condition.diagonal
        self.requires_frequencies = condition.requires_frequencies
        self.key = ('symmetric',condition.key)

    def evaluate_shifted(self,state,permutations,start,end,frequencies):
        return (self.condition.evaluate_shifted(state,permutations,start,end,frequencies)
                | self.condition.evaluate_shifted(state,permutations,end,start,frequencies))

    def evaluate_bound(self, state, permutations, start, end, frequencies, bindings=()):
        return (self.condition.evaluate_bound(state, permutations, start, end, frequencies, bindings)
                | self.condition.evaluate_bound(state, permutations, end, start, frequencies, bindings))

    def __repr__(self):
        return 'Symmetric[{}]'.format(repr(self.condition))


class PolynomialInvariantCondition(OperatorCondition):
    """Retain equal McUtils polynomial/function invariants at both endpoints.

    Variables of the McUtils functions refer to the supplied ``modes``.
    ``freeze_other_modes`` retains the spectator occupations unchanged.
    No target occupation, frequency, or Hamiltonian coefficient is captured.
    """
    symmetric = diagonal = True
    equivalence = True

    def __init__(self,functions,modes,freeze_other_modes=True,tolerance=0.):
        functions = tuple(functions) if isinstance(functions,(tuple,list)) else (functions,)
        snapshot = pickle.dumps(functions,protocol=5)
        self.functions = pickle.loads(snapshot)
        self.modes = _condition_modes(modes)
        self.freeze_other_modes = bool(freeze_other_modes)
        self.tolerance = 0. if tolerance == 0 else float(tolerance)
        if self.tolerance < 0 or np.isnan(self.tolerance):
            raise ValueError('a condition tolerance must be nonnegative')
        self.equivalence = self.tolerance == 0
        self.key = ('invariants',snapshot,self.modes,self.freeze_other_modes,self.tolerance)

    def evaluate_shifted(self,state,permutations,start,end,frequencies):
        source = np.asarray(state)[...,self.modes][...,None,:]
        left = source + _local_changes(permutations,start,self.modes)
        right = source + _local_changes(permutations,end,self.modes)
        result = np.ones(np.shape(permutations)[:-1],dtype=bool)
        for function in self.functions:
            result &= np.abs(_function_value(function,left)-_function_value(function,right)) <= self.tolerance
        if self.freeze_other_modes:
            result &= _outside_unchanged(permutations,np.asarray(end)-start,self.modes)
        return result

    def __repr__(self):
        return 'Invariant[{};modes={};freeze_other_modes={};tolerance={}]'.format(
            ','.join(map(repr,self.functions)),self.modes,self.freeze_other_modes,self.tolerance)

    def __reduce__(self):
        return type(self),(pickle.loads(self.key[1]),self.modes,self.freeze_other_modes,self.tolerance)


class PolynomialRelationshipCondition(OperatorCondition):
    """An arbitrary McUtils endpoint function compared with zero.

    The arguments are the selected source occupations followed by the target
    occupations. Use ``symmetrized`` for a directional relation in a unitary
    transformation; an invariant equality is symmetric by construction.
    """
    def __init__(self,function,modes,comparison='eq',tolerance=0.,freeze_other_modes=True):
        if comparison not in {'eq','ne','lt','le','gt','ge'}:
            raise ValueError('unknown polynomial comparison: '+comparison)
        snapshot = pickle.dumps(function,protocol=5)
        self.function = pickle.loads(snapshot)
        self.modes = _condition_modes(modes)
        self.comparison = comparison
        self.tolerance = 0. if tolerance == 0 else float(tolerance)
        if self.tolerance < 0 or np.isnan(self.tolerance):
            raise ValueError('a condition tolerance must be nonnegative')
        self.freeze_other_modes = bool(freeze_other_modes)
        self.key = ('relationship',snapshot,self.modes,comparison,self.tolerance,self.freeze_other_modes)

    def evaluate_shifted(self,state,permutations,start,end,frequencies):
        source = np.asarray(state)[...,self.modes][...,None,:]
        left = source + _local_changes(permutations,start,self.modes)
        right = source + _local_changes(permutations,end,self.modes)
        value = _function_value(self.function,np.concatenate([left,right],axis=-1))
        threshold = self.tolerance
        result = {'eq':lambda:abs(value)<=threshold,'ne':lambda:abs(value)>threshold,
                  'lt':lambda:value < -threshold,'le':lambda:value <= threshold,
                  'gt':lambda:value > threshold,'ge':lambda:value >= -threshold}[self.comparison]()
        if self.freeze_other_modes:
            result &= _outside_unchanged(permutations,np.asarray(end)-start,self.modes)
        return result

    def __reduce__(self):
        return type(self),(pickle.loads(self.key[1]),self.modes,self.comparison,
                           self.tolerance,self.freeze_other_modes)

    def __repr__(self):
        return 'Relationship[{};modes={};comparison={};tolerance={};freeze_other_modes={}]'.format(
            repr(self.function),self.modes,self.comparison,self.tolerance,self.freeze_other_modes)


class EnergyWindowCondition(OperatorCondition):
    """A harmonic energy-window relation, evaluated with runtime frequencies."""
    symmetric = diagonal = True
    requires_frequencies = True

    def __init__(self,tolerance,modes=None):
        self.tolerance = 0. if tolerance == 0 else float(tolerance)
        if self.tolerance < 0 or np.isnan(self.tolerance):
            raise ValueError('a condition tolerance must be nonnegative')
        self.modes = None if modes is None else _condition_modes(modes)
        self.key = ('harmonic-window',self.tolerance,self.modes)

    def evaluate_shifted(self,state,permutations,start,end,frequencies):
        change = np.asarray(end)-start
        result = abs(np.sum(frequencies*change[...,None,:],axis=-1)) <= self.tolerance
        if self.modes is not None:
            result &= _outside_unchanged(permutations,change,self.modes)
        return result

    def __repr__(self):
        return 'EnergyWindow[tolerance={};modes={}]'.format(self.tolerance,self.modes)


class IntegerLatticeCondition(OperatorCondition):
    """A displacement belongs to the integer lattice of resonance relations.

    Hermite reduction is done once symbolically. Runtime predicates use exact
    integer divisibility and reconstruction, not a floating null-space test.
    The family is valid at every occupation with nonnegative endpoints.
    """
    symmetric = diagonal = True
    equivalence = True

    def __init__(self,relations):
        relations = np.asarray(relations,dtype=int)
        if relations.ndim != 2:
            raise ValueError('relations must be a matrix of full mode changes')
        self.modes = tuple(int(mode) for mode in np.flatnonzero(np.any(relations != 0,axis=0)))
        if not self.modes:
            raise ValueError('a resonance lattice needs at least one nonzero relation')
        active = relations[:,self.modes]
        first = active[np.flatnonzero(np.any(active != 0,axis=1))[0]]
        primitive = first//np.gcd.reduce(abs(first))
        pivot = int(np.flatnonzero(primitive)[0])
        if primitive[pivot] < 0:
            primitive = -primitive
        coefficients = active[:,pivot]//primitive[pivot]
        if (np.all(active[:,pivot] % primitive[pivot] == 0)
                and np.all(coefficients[:,None]*primitive == active)):
            vector = primitive*int(np.gcd.reduce(abs(coefficients)))
            self.basis = tuple((int(value),) for value in vector)
            self.pivots,self.inverse_numerators = (pivot,),((1,),)
            self.denominator = int(vector[pivot])
            self.key = ('integer-lattice',self.modes,self.basis)
            return
        import sympy as sp
        from sympy.matrices.normalforms import hermite_normal_form
        basis = hermite_normal_form(sp.Matrix(relations[:,self.modes].T.tolist()))
        pivots = basis.T.rref()[1]
        inverse = basis[list(pivots),:].inv()
        denominator = math.lcm(*(int(value.q) for value in inverse))
        self.basis = tuple(tuple(int(value) for value in row) for row in basis.tolist())
        self.pivots = tuple(pivots)
        self.inverse_numerators = tuple(tuple(int(value*denominator) for value in row)
                                        for row in inverse.tolist())
        self.denominator = denominator
        self.key = ('integer-lattice',self.modes,self.basis)

    def evaluate_shifted(self,state,permutations,start,end,frequencies):
        change = np.asarray(end,dtype=int)-start
        local = _local_changes(permutations,change,self.modes)
        numerators = local[...,self.pivots] @ np.asarray(self.inverse_numerators,dtype=int).T
        coordinates = numerators//self.denominator
        return (np.all(numerators % self.denominator == 0,axis=-1)
                & np.all(coordinates @ np.asarray(self.basis,dtype=int).T == local,axis=-1)
                & _outside_unchanged(permutations,change,self.modes))

    def __repr__(self):
        return 'Lattice[modes={},basis={}]'.format(self.modes,self.basis)

    @functools.lru_cache(maxsize=16384)
    def may_retain_change(self,change):
        if not change:
            return True
        size = len(self.modes)
        if len(change) > size:
            return False
        if math.perm(size,len(change)) > 512:
            return None
        permutations = list(itertools.permutations(range(size),len(change)))
        candidates = np.zeros((len(permutations),size),dtype=int)
        candidates[np.arange(len(permutations))[:,None],permutations] = change
        numerators = candidates[:,self.pivots] @ np.asarray(self.inverse_numerators,dtype=int).T
        coordinates = numerators//self.denominator
        return bool(np.any(np.all(numerators % self.denominator == 0,axis=-1)
                    & np.all(coordinates @ np.asarray(self.basis,dtype=int).T == candidates,axis=-1)))


class ConditionEndpoint(int):
    """One of two adjacent endpoint rows in a native analytic path key."""
    _state_projection_factor = True
    _conditional_endpoint = True
    _cache = {}

    def __new__(cls,condition,kind,retain,endpoint):
        if kind not in {'mask','resolvent'} or endpoint not in {0,1}:
            raise ValueError('invalid conditional path factor')
        cache_key = (condition,kind,bool(retain),endpoint)
        cached = cls._cache.get(cache_key)
        if cached is not None:
            return cached
        payload = (condition.key,kind,bool(retain),endpoint)
        number = 2 + int.from_bytes(hashlib.sha256(repr(payload).encode()).digest()[:7],'big')
        obj = int.__new__(cls,number)
        obj.condition,obj.kind,obj.retain,obj.endpoint = condition,kind,bool(retain),endpoint
        obj._frozen = True
        if len(cls._cache) >= 1024:
            cls._cache.pop(next(iter(cls._cache)))
        cls._cache[cache_key] = obj
        return obj

    def __setattr__(self,name,value):
        if getattr(self,'_frozen',False):
            raise AttributeError('conditional path factors are immutable')
        object.__setattr__(self,name,value)

    def __eq__(self,other):
        if self is other:
            return True
        return (isinstance(other,ConditionEndpoint) and self.condition == other.condition
                and (self.kind,self.retain,self.endpoint)==(other.kind,other.retain,other.endpoint))

    def __ne__(self,other):
        return not self == other

    __hash__ = int.__hash__

    def __reduce__(self):
        return type(self),(self.condition,self.kind,self.retain,self.endpoint)

    def __repr__(self):
        return '{}Endpoint[{},{}]'.format(self.kind,self.endpoint,self.condition)


class ConditionalOperatorTerm(PerturbationTheoryTerm):
    """Restrict a native operator or solve its conditional homological equation."""
    _cache = {}

    def __init__(self,term,condition,kind='mask',retain=True):
        if kind == 'resolvent' and retain:
            raise ValueError('a conditional VPT resolvent must exclude its retained relation')
        super().__init__(logger=term.logger,polynomial_representation=term.polynomial_representation)
        self.term,self.condition,self.kind,self.retain = term,condition,kind,bool(retain)

    @classmethod
    def lookup(cls,term,condition,kind='mask',retain=True):
        key = (term,condition,kind,bool(retain))
        if key not in cls._cache:
            cls._cache[key] = cls(*key)
        return cls._cache[key]

    def get_changes(self):
        return {change:None for change in self.term.changes
                if (change or not (self.kind == 'resolvent' or
                                   (not self.retain and self.condition.diagonal)))
                and not (self.retain and self.condition.may_retain_change(change) is False)}

    def get_core_poly(self,changes,shift=None):
        base = self.term.get_poly_terms(changes)
        if np.isscalar(base) and base == 0:
            return 0
        if not changes and self.kind == 'mask' and self.condition.diagonal:
            return base if self.retain else 0
        # A harmonic homological divisor is translation invariant. Keep it
        # in the native energy algebra and retain the separate endpoint mask.
        rows = ((ConditionEndpoint(self.condition,'mask',self.retain,0),)+(0,)*len(changes),
                (ConditionEndpoint(self.condition,'mask',self.retain,1),)+tuple(changes))
        if not self.retain and self.condition.may_retain_change(changes) is False:
            # This occupation-independent proof is stronger than a runtime
            # endpoint test. Keep the native divisor and omit the true mask.
            rows = ()
        if self.kind == 'resolvent':
            rows = ((0,)+tuple(changes),)+rows
        if not rows:
            return base
        result = base.mul_simple(PTEnergyChangeProductSum({rows:1}))
        return result.scale(-1) if self.kind == 'resolvent' else result

    def __repr__(self):
        return '{}[{};term@{}]'.format('R' if self.kind=='resolvent' else 'D' if self.retain else 'Q',
                                      self.condition,id(self.term))


class ZeroTerm(PerturbationTheoryTerm):
    def get_changes(self):
        return {}

    def get_core_poly(self,changes,shift=None):
        return 0

    def __repr__(self):
        return '0'


@functools.lru_cache(maxsize=8192)
def canonical_condition_key(original):
    """Cache immutable scalar-factor keys, never numerical numerators."""
    width = max((len(row)-1 for row in original),default=0)
    rows = [tuple(row)+(0,)*(width+1-len(row)) for row in original]
    legacy = any(getattr(row[0],'_state_projection_factor',False)
                 and not isinstance(row[0],ConditionEndpoint) for row in rows)
    units,masks,scale = [],{},1.
    index,possible = 0,True
    while index < len(rows):
        row = rows[index]
        factor = row[0]
        if isinstance(factor,ConditionEndpoint):
            if index+1 == len(rows):
                raise ValueError('unpaired conditional endpoint')
            other = rows[index+1]
            if (factor.endpoint != 0 or not isinstance(other[0],ConditionEndpoint)
                    or other[0].endpoint != 1 or factor.condition != other[0].condition
                    or (factor.kind,factor.retain)!=(other[0].kind,other[0].retain)):
                raise ValueError('conditional endpoint order was lost')
            start,end = row[1:],other[1:]
            if factor.kind == 'mask' and factor.condition.symmetric and end < start:
                start,end = end,start
            pair = ((factor,)+start,(other[0],)+end)
            if factor.kind == 'mask':
                key = (factor.condition,start,end)
                previous = masks.get(key)
                if previous is not None:
                    if previous != factor.retain:
                        possible = False
                        break
                    index += 2
                    continue
                masks[key] = factor.retain
            units.append(pair)
            index += 2
            continue
        if not legacy and factor in {0,1}:
            delta = tuple(int(value) for value in row[1:])
            divisor = math.gcd(*delta) if delta else 0
            if divisor:
                first = next(value for value in delta if value)
                signed = divisor if first > 0 else -divisor
                delta = tuple(value//signed for value in delta)
                scale /= signed
            row = (0,)+tuple(delta)
        units.append((row,))
        index += 1
    if not possible:
        return None,0.
    units.sort()
    return tuple(row for unit in units for row in unit),scale


def canonicalize_condition_paths(terms):
    """Combine commutative/idempotent masks and primitive harmonic divisors."""
    result = {}
    for original,polynomial in terms.items():
        key,scale = canonical_condition_key(tuple(original))
        if key is None:
            continue
        width = len(key[0])-1 if key else 0
        if scale != 1:
            polynomial = polynomial*scale if np.isscalar(polynomial) else polynomial.scale(scale)
        if np.isscalar(polynomial) and polynomial != 1:
            if polynomial == 0:
                continue
            polynomial = make_product_polynomial([np.ones(1)]*width,prefactor=polynomial,
                                                  steps=0,representation='eager')
        result[key] = result.get(key,0)+polynomial
    return result


def format_condition_path(key):
    """Render masks as multiplicative path functions, never energy divisors."""
    masks,divisors = [],[]
    def vector(values):
        return '{'+','.join(map(str,values))+'}'
    def energy(delta):
        return '('+(' + '.join('{}*w[{}]'.format(value,i)
                               for i,value in enumerate(delta) if value) or '0')+')'
    index = 0
    while index < len(key):
        row=key[index]
        if isinstance(row[0],ConditionEndpoint):
            factor,end = row[0],key[index+1]
            masks.append('PathMask[{},{},n,pi,{},{}]'.format(
                json.dumps(repr(factor.condition)),'True' if factor.retain else 'False',
                vector(row[1:]),vector(end[1:])))
            if factor.kind=='resolvent':
                divisors.append(energy(np.subtract(row[1:],end[1:])))
            index += 2
        else:
            divisors.append(energy(row[1:]))
            index += 1
    result='*'.join(masks) or '1'
    return result+('/('+'*'.join(divisors)+')' if divisors else '')


def format_condition_sum(expression):
    parts=['({})*({})'.format(poly if np.isscalar(poly) else poly.format_expr(),format_condition_path(key))
           for key,poly in expression.terms.items()]
    result=' + '.join(parts) or '0'
    return result if expression.prefactor==1 else '{}*({})'.format(expression.prefactor,result)


def filter_condition_paths(terms,test_keys,mode):
    """Filter divisor paths while retaining every occupation condition."""
    def normalize(row):
        row = tuple(row)
        delta = np.asarray(row[1:],dtype=int)
        divisor = int(np.gcd.reduce(abs(delta))) if len(delta) else 0
        if divisor:
            divisor = divisor if delta[np.flatnonzero(delta)[0]] > 0 else -divisor
            delta = delta//divisor
        return (0,)+tuple(delta)
    normalized = ([[normalize(row) for row in key] for key in test_keys] if mode=='match'
                  else [normalize(row) for row in test_keys])
    result = {}
    for key,poly in terms.items():
        divisors = tuple(row for row in key if not getattr(row[0],'_state_projection_factor',False))
        query = TensorCoefficientPoly({divisors:1},canonicalize=False)
        if query.filter(normalized,mode=mode).terms:
            result[key] = poly
    return result


def operator_sum(*terms):
    """Keep a shared native DAG and cancel equal formal summands exactly."""
    grouped = {}
    pending = [(term,1.) for term in terms]
    while pending:
        term,scale = pending.pop()
        if isinstance(term,ZeroTerm) or (np.isscalar(term) and term == 0):
            continue
        if isinstance(term,ScaledPerturbationTheoryTerm):
            pending.append((term.base,scale*term.prefactor))
        elif isinstance(term,PerturbationTheoryTermSum):
            pending.extend((child,scale) for child in term.terms)
        else:
            grouped[term] = grouped.get(term,0.)+scale
    terms = [term if scale == 1 else scale*term for term,scale in grouped.items() if scale != 0]
    if not terms:
        return ZeroTerm(logger=False)
    return terms[0] if len(terms)==1 else PerturbationTheoryTermSum.lookup(*terms)


def operator_product(left,right):
    if isinstance(left,ZeroTerm) or isinstance(right,ZeroTerm):
        return ZeroTerm(logger=False)
    if isinstance(left,IdentityTerm):
        return right
    if isinstance(right,IdentityTerm):
        return left
    return left*right


def commutator(left,right):
    if left is right or isinstance(left,IdentityTerm) or isinstance(right,IdentityTerm):
        return ZeroTerm(logger=False)
    return operator_sum(operator_product(left,right),-operator_product(right,left))


@functools.lru_cache(maxsize=32768)
def operator_sector(term,condition):
    """Prove a retained/complementary sector without inspecting any state."""
    if isinstance(term,IdentityTerm):
        return 'retained'
    if isinstance(term,OperatorExpansionTerm):
        changes = [change for change in term.changes
                   if not (len(change)==2 and isinstance(change[0],tuple))]
        if changes and all(not change for change in changes):
            return 'retained'
        if changes and all(condition.may_retain_change(change) is False for change in changes):
            return 'complement'
    if isinstance(term,ConditionalOperatorTerm) and term.condition == condition:
        return 'retained' if term.kind == 'mask' and term.retain else 'complement'
    if isinstance(term,ScaledPerturbationTheoryTerm):
        return operator_sector(term.base,condition)
    if isinstance(term,ShiftedEnergyBaseline) and condition.symmetric:
        return operator_sector(term.base,condition)
    if isinstance(term,PerturbationTheoryTermSum):
        sectors = {operator_sector(child,condition) for child in term.terms}
        return next(iter(sectors)) if len(sectors)==1 else None
    if condition.equivalence and hasattr(term,'gen1') and hasattr(term,'gen2'):
        left,right = operator_sector(term.gen1,condition),operator_sector(term.gen2,condition)
        if left == 'retained' and right == 'retained':
            return 'retained'
        if {left,right} == {'retained','complement'}:
            return 'complement'
    return None


class ConditionalTransformationStack:
    """Ordered all-order unitary transformations under a symbolic relation."""
    def __init__(self,solver,condition, *, block_recurrence=False):
        if not isinstance(condition,OperatorCondition) or not condition.symmetric:
            raise ValueError('unitary resonance handling requires a symmetric operator condition')
        self.solver = solver
        self.condition = condition if condition.diagonal else (SameStateCondition() | condition)
        self.identity = IdentityTerm(logger=solver.logger,
                                     polynomial_representation=solver.polynomial_representation)
        self.zero = ZeroTerm(logger=solver.logger,
                             polynomial_representation=solver.polynomial_representation)
        self.order = len(solver.hamiltonian_expansion)-1
        self.block_recurrence = bool(block_recurrence)
        self.hamiltonian = ([solver.hamiltonian_expansion[0]]+[self.zero]*self.order
                            if self.block_recurrence else list(solver.hamiltonian_expansion))
        self.normalizers = [self.identity]
        self.outside = [self.zero]
        self.wavefunctions = [self.identity]+[self.zero]*self.order
        self.generators = [self.zero]
        self.stages = [tuple(self.hamiltonian)]
        self._transformed_operators = {}
        _stacks.add(self)

    def restrict(self,term,retain=True,resolvent=False):
        if isinstance(term,ZeroTerm):
            return term
        sector = operator_sector(term,self.condition)
        if sector is not None:
            if (sector == 'retained') != retain:
                return self.zero
            if not resolvent:
                return term
        if isinstance(term,PerturbationTheoryTermSum):
            return operator_sum(*(self.restrict(child,retain,resolvent) for child in term.terms))
        if (resolvent and isinstance(term,ConditionalOperatorTerm)
                and term.condition == self.condition and term.kind == 'mask' and not term.retain):
            term = term.term
        return ConditionalOperatorTerm.lookup(term,self.condition,
                         kind='resolvent' if resolvent else 'mask',retain=retain)

    def _conjugate(self,series,generator,stage):
        output = list(series)
        for source,term in enumerate(series):
            nested = term
            for power in range(1,(self.order-source)//stage+1):
                nested = commutator(nested,generator)
                target = source+power*stage
                output[target] = operator_sum(output[target],(1/math.factorial(power))*nested)
        return output

    def _extend_block_recurrence(self, order, logger):
        """HY=YF with shared symbolic F, QY and full block normalization.

        All contractions are native operator DAG products. Resolvents carry
        frequency-change divisors; no reference energy or evaluated F entry
        is captured. Keep F[n] itself as a shared node in subsequent orders.
        """
        H = self.solver.hamiltonian_expansion
        while len(self.normalizers) <= order:
            n = len(self.normalizers)
            with _log_evaluation_time(logger, 'Building symbolic block recurrence at order {o}', o=n):
                overlap = operator_sum(*(operator_product(
                    ShiftedEnergyBaseline.lookup(self.wavefunctions[i]), self.wavefunctions[n-i])
                    for i in range(1, n)))
                normalizer = operator_sum(-.5*self.restrict(overlap))
                driving = operator_sum(*(operator_product(H[r], self.wavefunctions[n-r])
                                         for r in range(1, n+1)))
                inside_feedback = operator_sum(*(operator_product(
                    self.normalizers[n-r], self.hamiltonian[r]) for r in range(1, n)))
                # F=[H0,N]+B is Hermitian and N is Hermitian. Therefore
                # F=(B+B^dagger)/2; avoid expanding two cancelling H0*N
                # products with artificial extra mode-index contractions.
                inside = self.restrict(operator_sum(driving, -inside_feedback))
                F = (inside if n == 1 else operator_sum(
                    .5*inside, .5*ShiftedEnergyBaseline.lookup(inside)))
                feedback = operator_sum(*(operator_product(
                    self.outside[n-r], self.hamiltonian[r]) for r in range(1, n)))
                X = self.restrict(operator_sum(driving, -feedback), retain=False, resolvent=True)
                self.normalizers.append(normalizer)
                self.outside.append(X)
                self.hamiltonian[n] = F
                self.wavefunctions[n] = operator_sum(X, normalizer)
        return self

    def extend(self,order=None,logger=None):
        logger = Logger.lookup(self.solver.logger if logger is None else logger)
        order = self.order if order is None else order
        if order > self.order:
            raise ValueError('requested order exceeds the supplied Hamiltonian expansion')
        if self.block_recurrence:
            return self._extend_block_recurrence(order, logger)
        while len(self.generators) <= order:
            stage = len(self.generators)
            with _log_evaluation_time(logger, "Building conditional transformation at order {o}", o=stage):
                old = self.hamiltonian
                retained = self.restrict(old[stage])
                eliminated = self.restrict(old[stage],retain=False)
                generator = self.restrict(old[stage],retain=False,resolvent=True)
                # Combine the H0 and K_stage BCH contributions at every order:
                # ad_S^p(D K_stage)/p! + p ad_S^p(Q K_stage)/(p+1)!.
                # This is the exact homological identity, not a low-order rule.
                remaining = list(old)
                remaining[0] = remaining[stage] = self.zero
                updated = self._conjugate(remaining,generator,stage)
                updated[0],updated[stage] = old[0],retained
                d_nested,q_nested = retained,eliminated
                for power in range(1,(self.order-stage)//stage+1):
                    d_nested = commutator(d_nested,generator)
                    q_nested = commutator(q_nested,generator)
                    target = (power+1)*stage
                    updated[target] = operator_sum(updated[target],
                        (1/math.factorial(power))*d_nested,
                        (power/math.factorial(power+1))*q_nested)
                before = self.wavefunctions
                after = list(before)
                power_term = self.identity
                for power in range(1,self.order//stage+1):
                    power_term = operator_product(power_term,generator)
                    for source in range(self.order-power*stage+1):
                        target = source+power*stage
                        after[target] = operator_sum(after[target],
                            (1/math.factorial(power))*operator_product(before[source],power_term))
                self.wavefunctions,self.hamiltonian = after,updated
                self.generators.append(generator)
                self.stages.append(tuple(updated))
        return self

    def transform_operator(self,operators,order=None,terms=None,branches=True,logger=None):
        order = self.order if order is None else order
        logger = Logger.lookup(self.solver.logger if logger is None else logger)
        self.extend(order,logger=logger)
        allowed = None if terms is None else frozenset(tuple(term) for term in terms)
        branch_key = branches if isinstance(branches,bool) else tuple(sorted(
            branches.items() if isinstance(branches,dict) else ((name,True) for name in branches)))
        key = (tuple(operators),order,allowed,branch_key)
        if key not in self._transformed_operators:
            with _log_evaluation_time(logger, "Transforming conditional operator through order {o}", o=order):
                if allowed is None and branches is True and not self.block_recurrence:
                    series = list(operators[:self.order+1])+[self.zero]*max(0,self.order+1-len(operators))
                    for stage,generator in enumerate(self.generators[1:order+1],1):
                        series = self._conjugate(series,generator,stage)
                    self._transformed_operators[key] = series[:order+1]
                    return self._transformed_operators[key]
                corrections = []
                for n in range(order+1):
                    parts = []
                    for ket in range(n+1):
                        for bra in range(n-ket+1):
                            op = n-ket-bra
                            if op >= len(operators) or (allowed is not None and (ket,op,bra) not in allowed):
                                continue
                            if branches is True:
                                selections = [(False,False)]
                            else:
                                selected = {'main':True}
                                if branches is not False:
                                    selected.update(dict(branch_key))
                                selections = [(bc,kc) for name,bc,kc in
                                    [('main',False,False),('left',False,True),
                                     ('right',True,False),('both',True,True)] if selected.get(name,False)]
                            for bra_cross,ket_cross in selections:
                                left = self.wavefunctions[bra]
                                right = self.wavefunctions[ket]
                                if branches is not True:
                                    diagonal = SameStateCondition()
                                    def piece(term,n,cross):
                                        if n == 0:
                                            return self.zero if cross else self.identity
                                        condition = ((self.condition & ~diagonal) if cross else
                                                     (~self.condition | diagonal))
                                        return ConditionalOperatorTerm.lookup(term,condition)
                                    left,right = piece(left,bra,bra_cross),piece(right,ket,ket_cross)
                                left = self.identity if isinstance(left,IdentityTerm) else ShiftedEnergyBaseline.lookup(left)
                                parts.append(operator_product(left,operator_product(operators[op],right)))
                    corrections.append(operator_sum(*parts))
                self._transformed_operators[key] = corrections
        return self._transformed_operators[key]

    def intermediate_wavefunctions(self,order,logger=None):
        logger = Logger.lookup(self.solver.logger if logger is None else logger)
        self.extend(order,logger=logger)
        if not self.condition.equivalence:
            raise ValueError('intermediate normalization requires an equivalence relation')
        with _log_evaluation_time(logger, "Building intermediate normalization through order {o}", o=order):
            inside = [self.identity]+[self.restrict(term) for term in self.wavefunctions[1:order+1]]
            inverse = [self.identity]
            result = [self.identity]
            for n in range(1,order+1):
                inverse.append(-operator_sum(*(operator_product(inside[i],inverse[n-i])
                                                 for i in range(1,n+1))))
                result.append(operator_sum(*(operator_product(self.wavefunctions[i],inverse[n-i])
                                             for i in range(n+1))))
        return result


def conditional_path_factors(energy_key,state,permutations,perm_freqs,baseline_shift=None,
                             *,projection_bindings=()):
    """Evaluate endpoint conditions and safe denominators before division."""
    state = np.asarray(state)[0]
    permutations = np.asarray(permutations,dtype=int)
    good = np.ones(len(permutations),dtype=bool)
    denominator = np.ones(len(permutations))
    width = permutations.shape[-1]
    baseline = np.zeros(width,dtype=int)
    if baseline_shift is not None:
        baseline[:min(width,len(baseline_shift))] = baseline_shift[:width]
    index = 0
    while index < len(energy_key):
        row = energy_key[index]
        if not isinstance(row[0],ConditionEndpoint):
            if getattr(row[0],'_state_projection_factor',False):
                mask,weights = projected_path_factors((row,),[state],permutations,perm_freqs,
                    projection_bindings=projection_bindings)
                good &= mask
                denominator *= weights
            else:
                denominator *= perm_freqs @ np.asarray(row[1:])
            index += 1
            continue
        factor = row[0]
        end_row = energy_key[index+1]
        end_factor = end_row[0]
        if (factor.endpoint != 0 or not isinstance(end_factor,ConditionEndpoint)
                or end_factor.endpoint != 1 or end_factor.condition != factor.condition
                or (factor.kind,factor.retain) != (end_factor.kind,end_factor.retain)):
            raise ValueError('conditional endpoint pair was not preserved by the analytic algebra')
        start,end = np.asarray(row[1:],dtype=int)+baseline,np.asarray(end_row[1:],dtype=int)+baseline
        possible = (np.all(state[permutations]+start >= 0,axis=-1)
                    & np.all(state[permutations]+end >= 0,axis=-1))
        matches = factor.condition.evaluate_bound(state,permutations,start,end,perm_freqs,projection_bindings)
        good &= possible & (matches if factor.retain else ~matches)
        if factor.kind == 'resolvent':
            denominator *= perm_freqs @ (start-end)
        index += 2
    return good,denominator


def conditional_energy_rows(energy_keys,source_states,permutations,frequency_blocks,
                            baseline_shift=None,*,projection_bindings=()):
    values = [conditional_path_factors(key,state,perms,freqs,baseline_shift,projection_bindings=projection_bindings)
              for key,state,perms,freqs in zip(energy_keys,source_states,permutations,frequency_blocks)]
    return tuple(np.stack([value[i] for value in values]) for i in range(2))


# ---------------------------------------------------------------------------------------------
#  ConditionalEvaluator (formerly VPT2/ConditionalEvaluator.py)
#  Use the native batched evaluator for condition-based VPT expressions.
# ---------------------------------------------------------------------------------------------

class CorrectionView(PerturbationTheoryTerm):
    """A final native filtering/checkpoint view without mutating a shared DAG."""
    def __init__(self,term,diagonal=False,**options):
        super().__init__(logger=options.get('logger',term.logger),polynomial_representation=term.polynomial_representation,
            **{key:value for key,value in options.items() if key in
               {'allowed_terms','allowed_energy_changes','allowed_coefficients','disallowed_coefficients'}})
        self.term = term
        self.diagonal = bool(diagonal)

    def get_changes(self):
        return {change:None for change in self.term.changes if not self.diagonal or not change}

    def get_core_poly(self,changes,shift=None):
        return self.term.get_poly_terms(changes)

    def __repr__(self):
        return 'ConditionalCorrection[{}]'.format(self.term)


class ConditionalCorrectionEvaluator:
    """Late-bound occupations/coefficients; no finite block contractions."""
    def __init__(self,evaluator,expansions=None,freqs=None,zero_cutoff=None,parallelizer=None,
                 logger=None,verbose=False, *, blocks=None, pairs=None):
        self.evaluator = evaluator
        self.logger = Logger.lookup(evaluator.solver.logger if logger is None else logger)
        self.verbose = verbose
        resolved = ProjectedCorrectionEvaluator.resolve_blocks(blocks, pairs)
        groups = tuple(sorted(set(resolved.values())))
        self.stack = evaluator.solver.get_conditional_stack(
            ResonanceBlockCondition(), block_recurrence=evaluator.block_recurrence)
        self.projection_bindings = dict(
            projectors=(groups,), scalar_count=0, scalar_order=len(evaluator.expansions))
        self.raw_expansions = evaluator.expansions if expansions is None else expansions
        self.expansions = evaluator._prep_expansions(self.raw_expansions)
        self.freqs = evaluator.freqs if freqs is None else freqs
        self.zero_cutoff = zero_cutoff
        self.parallelizer = evaluator.parallelizer if parallelizer is None else parallelizer
        self.options = {}

    def validate_options(self,options):
        supported = {'allowed_coefficients','disallowed_coefficients','epaths',
            'only_degenerate_terms','include_degenerate_correction_terms','zero_cutoff',
            'parallelizer','target_orders','return_sorted','logger','log_scaled','terms'}
        unknown = set(options)-supported
        if unknown:
            raise TypeError('unknown conditional correction options: {}'.format(sorted(unknown)))
        self.options = {key:value for key,value in options.items() if key in
                        {'allowed_coefficients','disallowed_coefficients','epaths','return_sorted'}}

    def _diagonal_values(self,term,states,**numerical):
        from .Evaluators import current_compilation
        changes = self.evaluator.get_diff_map([[state,[state]] for state in states], logger=self.logger)
        values = []
        for change,permutations in changes.items():
            with _log_evaluation_time(self.logger, "Building evaluator {s}[{c}]", s=term, c=change):
                if current_compilation() is not None and current_compilation()[0].using_artifacts:
                    expression = term(change)
                    expression.logger = self.logger
                else:
                    expression = PerturbationTheoryExpressionEvaluator(
                        term,term.get_poly_terms(change),change,logger=self.logger)
            with _log_evaluation_time(self.logger, "evaluating..."):
                values.extend(expression.evaluate(permutations,self.expansions,self.freqs,
                    zero_cutoff=self.zero_cutoff,parallelizer=self.parallelizer,log_scaled=False,
                    verbose=self.verbose,projection_bindings=self.projection_bindings,**numerical))
        return np.asarray(values)[:,0]

    def energies(self,states,order,**numerical):
        self.stack.extend(order,logger=self.logger)
        values = []
        for n in range(0,order+1,2):
            with self.logger.block(tag="Getting corrections at order {o}", o=n):
                values.append(self._diagonal_values(self.stack.hamiltonian[n],states,**numerical))
        return values

    def overlaps(self,states,order,**numerical):
        self.stack.extend(order,logger=self.logger)
        values = []
        for n in range(order+1):
            with self.logger.block(tag="Getting overlap corrections at order {o}", o=n):
                values.append(self._diagonal_values(self.stack.wavefunctions[n],states,**numerical))
        return np.stack(values,axis=-1)

    def _evaluate(self,series,states,order,expansions=None,target_orders=None,sign=1.):
        def generator(n,**options):
            term = series[n] if sign == 1 else sign*series[n]
            return CorrectionView(term,logger=self.logger,**options)
        return self.evaluator.get_state_by_state_corrections(generator,states,order=order,
            expansions=self.raw_expansions if expansions is None else expansions,freqs=self.freqs,
            zero_cutoff=self.zero_cutoff,parallelizer=self.parallelizer,
            only_degenerate_terms=False,include_degenerate_correction_terms=False,
            target_orders=target_orders,logger=self.logger,verbose=self.verbose,
            projection_bindings=self.projection_bindings,**self.options)

    def wavefunctions(self,states,order,full=True,target_orders=None):
        self.stack.extend(order,logger=self.logger)
        if self.evaluator.solver.intermediate_normalization:
            series = self.stack.intermediate_wavefunctions(order,logger=self.logger)
        else:
            series = self.stack.wavefunctions
        if not full:
            series = [self.stack.restrict(term,retain=False) for term in series]
        return self._evaluate(series,states,order,target_orders=target_orders)

    def hamiltonian(self,states,order,only_degenerate_terms=True,terms=None,branches=True,
                    target_orders=None, normalization='orthonormal'):
        self.stack.extend(order,logger=self.logger)
        if normalization == 'intermediate':
            if terms is not None or branches is not True:
                raise ValueError('Bloch entries are not normalized bra-H-ket contribution triples')
            omega = self.stack.intermediate_wavefunctions(order, logger=self.logger)
            H = self.evaluator.solver.hamiltonian_expansion
            series = [H[0]]+[self.stack.restrict(operator_sum(*(
                operator_product(H[r], omega[n-r]) for r in range(1,n+1))))
                for n in range(1,order+1)]
            return self._evaluate(series,states,order,target_orders=target_orders,
                                  sign=-1. if only_degenerate_terms else 1.)
        series = (self.stack.hamiltonian if terms is None and branches is True else
                  self.stack.transform_operator(self.evaluator.solver.hamiltonian_expansion,
                                                order,terms,branches,logger=self.logger))
        return self._evaluate(series,states,order,target_orders=target_orders,
                              sign=-1. if only_degenerate_terms else 1.)

    def operators(self,states,order,operators,expansions,terms=None,branches=True,
                  target_orders=None):
        series = self.stack.transform_operator(operators,order,terms,branches,logger=self.logger)
        return self._evaluate(series,states,order,expansions,target_orders)


class PerturbationTheoryEvaluator:

    degeneracy_handling_modes = ('projected', 'filtered', 'conditional')
    # Filtered effective Hamiltonian, off-diagonal elements: True evaluates only the
    # resonant terms of the main <n|H|m> expression and negates them (see
    # _compute_corr) instead of evaluating the non-resonant complement. Both give
    # the same numbers; on HNO3/criegee the resonant set is the larger one (about
    # two thirds of the (term, assignment) pairs), so the complement is the default.
    use_resonant_complement = False

    @classmethod
    def resolve_degeneracy_handling(cls, mode, evaluation_mode=None):
        """
        With no explicit method, actions selects 'filtered'; direct selects 'projected'.
        'projected': degenerate blocks are treated by block projection
        (`ProjectedCorrectionEvaluator`).
        'conditional': resonance membership is inferred from the same resolved
        degeneracy-spec blocks/pairs as the other modes and bound at evaluation.
        No independent resonance_condition is supplied to the runner. block_recurrence=True
        (default) selects the shared symbolic F/Q/N recurrence; False selects
        the existing ordered Conditional transformation. Neither flag value
        changes 'projected' dispatch.
        'filtered': the native change-indexed recurrence. With resolved blocks,
        block_recurrence=True (default) includes shared symbolic F[k] feedback
        and complete block normalization. False keeps the original scalar
        filtered recurrence. Integer change domains bind at evaluation time;
        neither choice uses the Conditional stack or numerical projection.
        """
        if mode is None:
            mode = 'filtered' if evaluation_mode == 'actions' else cls.degeneracy_handling_modes[0]
        if mode not in cls.degeneracy_handling_modes:
            raise ValueError("degeneracy_handling must be one of {}, got {!r}".format(
                cls.degeneracy_handling_modes, mode))
        return mode

    def _projected_corrections(self, blocks, pairs, expansions=None, freqs=None, zero_cutoff=None,
                               parallelizer=None, logger=None, verbose=False):
        from .Evaluators import _current_evaluation_mode
        if (_current_evaluation_mode.get() == 'actions'
                and (self.degeneracy_handling == 'conditional'
                     or (self.degeneracy_handling == 'projected' and (blocks is not None or pairs is not None)))):
            raise ValueError('actions evaluation requires the native filtered recurrence')
        self._prepare_filtered_recurrence(blocks, pairs)
        if self.degeneracy_handling == 'conditional':
            return ConditionalCorrectionEvaluator(self, expansions, freqs, zero_cutoff, parallelizer,
                logger=logger, verbose=verbose, blocks=blocks, pairs=pairs)
        if blocks is None and pairs is None:
            return None
        if self.degeneracy_handling == 'filtered':
            return None
        evaluator = (self if parallelizer is None else
                     type(self)(self.solver, self.expansions, self.freqs, parallelizer=parallelizer,
                        expressions_folder=self.expressions_folder,compilation_options=self.compilation_options,
                        evaluation_options=self.evaluation_options,compile=self.compile_expressions,
                        degeneracy_handling=self.degeneracy_handling,block_recurrence=self.block_recurrence))
        return ProjectedCorrectionEvaluator(evaluator, blocks=blocks, pairs=pairs,
            expansions=expansions, freqs=freqs, zero_cutoff=zero_cutoff,
            logger=logger, verbose=verbose)

    def _prepare_filtered_recurrence(self, blocks, pairs):
        groups = ()
        if (self.degeneracy_handling == 'filtered' and self.block_recurrence
                and self.solver.resonance_condition is None):
            resolved = ProjectedCorrectionEvaluator.resolve_blocks(blocks, pairs)
            groups = tuple(sorted({group for group in resolved.values() if len(group)>1}))
        enabled = bool(groups)
        if self.solver.filtered_block_recurrence != enabled:
            self.solver = copy.copy(self.solver)
            self.solver.filtered_block_recurrence = enabled
        self._filtered_change_bindings = (None if not enabled else dict(
            projectors=(groups,), scalar_count=0, scalar_order=len(self.expansions)))

    @staticmethod
    def _validate_projected_options(opts):
        # Complete corrections and individual (bra,operator,ket) orders are
        # supported. A coefficient/denominator filter after contracting a
        # finite block would lose the provenance of its numerical entries.
        for key in ('allowed_coefficients', 'disallowed_coefficients', 'epaths'):
            if opts.get(key) is not None:
                raise NotImplementedError('runtime {} filtering of projected block diagrams is not supported'.format(key))
        unknown = set(opts) - {'allowed_coefficients','disallowed_coefficients','epaths',
            'only_degenerate_terms','include_degenerate_correction_terms','zero_cutoff',
            'parallelizer','target_orders','return_sorted','logger','log_scaled','terms'}
        if unknown:
            raise TypeError('unknown projected correction options: {}'.format(sorted(unknown)))

    @staticmethod
    def _state_degeneracy_pairs(degenerate_states, state):
        return degenerate_states.get(tuple(int(x) for x in state), [])

    def __init__(self, solver:AnalyticPerturbationTheorySolver, expansion, freqs=None, parallelizer=None,
                 expressions_folder=None, compilation_options=None, evaluation_options=None, compile=True,
                 degeneracy_handling=None, block_recurrence=True, evaluation_mode=None):
        if not isinstance(block_recurrence, (bool, np.bool_)):
            raise TypeError("block_recurrence must be boolean")
        self.block_recurrence = bool(block_recurrence)
        self.expressions_folder = expressions_folder
        self.compile_expressions = compile is not False
        self._expression_store = None
        self.compilation_options = dict(compilation_options or {})
        self.evaluation_options = dict(evaluation_options or {})
        if evaluation_mode is not None:
            self.evaluation_options['evaluation_mode'] = evaluation_mode
        self.solver = solver
        self.expansions = expansion
        if freqs is None:
            # Match the harmonic diagonal used by the classic solver.
            potential, kinetic = self.expansions[0][:2]
            kinetic_diagonal = np.diag(kinetic) if np.ndim(kinetic) == 2 else kinetic
            freqs = np.diag(potential) + kinetic_diagonal
        self.freqs = freqs
        from .Evaluators import _validate_evaluation_mode
        self.evaluation_options['evaluation_mode'] = _validate_evaluation_mode(
            self.evaluation_options.get('evaluation_mode'), nmodes=len(freqs))
        self.degeneracy_handling = self.resolve_degeneracy_handling(
            degeneracy_handling, self.evaluation_options['evaluation_mode'])
        self.parallelizer = Parallelizer.lookup(parallelizer) if parallelizer is not None else None

    @property
    def expression_store(self):
        """Lazily reuse a bounded RAM store, or the explicitly configured disk store."""
        from .Evaluators import AnalyticExpressionStore
        folder=None if self.expressions_folder is None else Path(self.expressions_folder).expanduser().resolve()
        budget=self.compilation_options.get('cache_bytes',512*1024**2)
        if (self._expression_store is None or self._expression_store.folder != folder
                or folder is None and self._expression_store.memory_cache_bytes != budget):
            self._expression_store=AnalyticExpressionStore(folder,
                memory_cache_bytes=self.compilation_options.get('cache_bytes',512*1024**2))
        return self._expression_store

    def compilation_statistics(self):
        folder=None if self.expressions_folder is None else Path(self.expressions_folder).expanduser().resolve()
        budget=self.compilation_options.get('cache_bytes',512*1024**2)
        if (self._expression_store is None or self._expression_store.folder != folder
                or folder is None and self._expression_store.memory_cache_bytes != budget):return dict(requests=0,compilations=0,cache_hits=0,failures=0,
            compilation_seconds=0.,cache_lookup_seconds=0.,preparation_seconds=0.,evaluation_seconds=0.)
        return self._expression_store.statistics()

    def precompile(self, correction='energy', order=None, **options):
        from .Evaluators import AnalyticCorrectionCompiler
        if not self.compile_expressions:raise ValueError('compilation disabled by compile=False')
        return AnalyticCorrectionCompiler(self).precompile(correction=correction, order=order, **options)

    def precompile_action_expansion(self, order=2, correction='energy', **options):
        from .Evaluators import AnalyticCorrectionCompiler
        if not self.compile_expressions:
            raise ValueError('compilation disabled by compile=False')
        return AnalyticCorrectionCompiler(self).precompile_action_expansion(order,correction,**options)

    @_compiled_evaluation_context
    def get_action_expansion(self, order=2, correction='energy', **options):
        """Return ([C, w, X, Y, ...], symbolic kernels) in actions I=n+1/2.

        Tensor r stores each monomial coefficient once at nondecreasing mode
        indices i1<=...<=ir, with no permutation multiplicity factor. The
        maximum rank is order//2+1. Odd perturbative orders vanish under the
        native VPT parity; all lower even corrections are included. Both
        source choices are independent of runtime resonance groups.
        """
        from .Evaluators import AnalyticCorrectionCompiler
        return AnalyticCorrectionCompiler(self).action_expansion(order,correction,**options)

    def precompile_frequency_corrections(self, initial_state, final_state, order=None, **options):
        from .Evaluators import AnalyticCorrectionCompiler
        if not self.compile_expressions:raise ValueError('compilation disabled by compile=False')
        return AnalyticCorrectionCompiler(self).precompile_frequencies(initial_state, final_state, order=order, **options)

    def get_frequency_corrections(self, initial_state, final_state, order=None, **options):
        from .Evaluators import AnalyticCorrectionCompiler
        if not self.compile_expressions:raise ValueError('compilation disabled by compile=False')
        return AnalyticCorrectionCompiler(self).frequencies(initial_state, final_state, order=order, **options)

    def modify_hamiltonian(self, hamiltonian_corrections):
        if hamiltonian_corrections is None:
            return self

        hamiltonian_corrections, expansion_corrections = self._prep_ham_corrs(hamiltonian_corrections)
        exps = self._prep_operator_expansion(self.expansions, expansion_corrections)

        return type(self)(
            self.solver.modify_hamiltonian(hamiltonian_corrections),
            exps, parallelizer=self.parallelizer,
            expressions_folder=self.expressions_folder,compilation_options=self.compilation_options,
            evaluation_options=self.evaluation_options,compile=self.compile_expressions,
            degeneracy_handling=self.degeneracy_handling, block_recurrence=self.block_recurrence
        )

    @_timed_evaluation('Energy corrections (total)')
    @_compiled_evaluation_context
    def get_energy_corrections(self, states, order=None, expansions=None, freqs=None,
                               zero_cutoff=None, degenerate_states=None, verbose=False,
                               logger=None,
                               parallelizer=None,
                               evaluation_mode=None,
                               dag_cache_size=None,
                               dag_cache_bytes=None,
                               dag_chunk_size=None, degenerate_blocks=None):
        projected = self._projected_corrections(degenerate_blocks, degenerate_states,
                                               expansions, freqs, zero_cutoff, parallelizer,
                                               logger=logger, verbose=verbose)
        if projected is not None:
            if order is None: order = len(self.expansions) - 1
            if hasattr(projected,'energies'):
                return projected.energies(states,order,evaluation_mode=evaluation_mode,
                    dag_cache_size=dag_cache_size,dag_cache_bytes=dag_cache_bytes,
                    dag_chunk_size=dag_chunk_size)
            return projected.diagonal(states, order, 'F', range(0, order+1, 2))
        if isinstance(degenerate_states, dict) and not self.solver.filtered_block_recurrence:
            results = [
                self.get_energy_corrections(
                    [state], order=order, expansions=expansions, freqs=freqs,
                    zero_cutoff=zero_cutoff,
                    degenerate_states=self._state_degeneracy_pairs(degenerate_states, state),
                    verbose=verbose, logger=logger, parallelizer=parallelizer,
                    evaluation_mode=evaluation_mode,
                    dag_cache_size=dag_cache_size, dag_cache_bytes=dag_cache_bytes,
                    dag_chunk_size=dag_chunk_size
                )
                for state in states
            ]
            if not results:
                raise ValueError("energy corrections require at least one state")
            return [np.concatenate([result[i] for result in results]) for i in range(len(results[0]))]
        expansions = self._prep_expansions(expansions)
        if freqs is None: freqs = self.freqs
        if order is None: order = len(self.expansions) - 1

        degenerate_changes = (None if self.solver.filtered_block_recurrence else
                              self.get_degenerate_changes(degenerate_states))
        degeneracy_identification_context = (
            None if degenerate_changes is None else
            DegeneracyIdentificationContext(degenerate_changes)
        )

        if logger is None: logger = self.solver.logger
        logger = Logger.lookup(logger)
        if parallelizer is None: parallelizer = self.parallelizer
        parallelizer = Parallelizer.lookup(parallelizer) if parallelizer is not None else None

        PerturbationTheoryExpressionEvaluator._cached_expansion = None # to make paralellism work right
        from .Evaluators import _current_evaluation_mode
        generator=self.solver.energy_correction
        if (_current_evaluation_mode.get()=='actions' and self.block_recurrence
                and not self.solver.filtered_block_recurrence):
            generator=self.solver._action_expansion_source('energy',order).energy_correction
        corrs = []
        for o in range(0, order + 1, 2):
            with logger.block(tag="Getting corrections at order {o}", o=o):
                evaluator = generator(o)([])
                with logger.block(tag="evaluating..."):
                    start = time.time()
                    subcorrs = evaluator.evaluate(
                        [[s, [np.arange(len(s))]] for s in states],
                        expansions, freqs, verbose=verbose,
                        degenerate_changes=degenerate_changes,
                        projection_bindings=getattr(self, '_filtered_change_bindings', None), zero_cutoff=zero_cutoff, log_scaled=True,
                        parallelizer=parallelizer,
                        evaluation_mode=evaluation_mode,
                        dag_cache_size=dag_cache_size,
                        dag_cache_bytes=dag_cache_bytes,
                        dag_chunk_size=dag_chunk_size,
                        _degeneracy_identification_context=(
                            degeneracy_identification_context
                        )
                    )

                    end = time.time()
                    logger.log_print('took {e:.3f}s', e=end - start)
                corrs.append(np.concatenate(subcorrs))

        return corrs

    @staticmethod
    def is_single_expansion(expansion, min_order=None):
        if min_order is None:
            return (
                    nput.is_numeric(expansion[0][0])
                    or isinstance(expansion[0][0], np.ndarray)
                    or (not nput.is_numeric(expansion[0][0][0]) and nput.is_numeric(expansion[0][0][0][0]))
            )
        else:
            def _extract(a, m): return _extract(a[0], m-1) if m > 0 else a
            return (
                nput.is_numeric(expansion[0])
                or (
                        not nput.is_numeric(expansion[0][0])
                        and nput.is_numeric(_extract(expansion[0], min_order))
                )
            )


    def _prep_expansions(self, expansions):
        #TODO: check symmetry of G-matrix before manipulating
        if expansions is None: expansions = self.expansions
        return [
            [
                tensor
                    if i != 1 else
                ( np.moveaxis(tensor, -1, 0) if not nput.is_numeric(tensor) else tensor )
                for i,tensor in enumerate(order_expansion)
            ]
            for order_expansion in expansions
        ]
    @_timed_evaluation('Overlap corrections (total)')
    @_compiled_evaluation_context
    def get_overlap_corrections(self, states, order=None, expansions=None,
                                degenerate_states=None, freqs=None,
                                zero_cutoff=None, verbose=False,
                                parallelizer=None,
                                evaluation_mode=None,
                                dag_cache_size=None,
                                dag_cache_bytes=None,
                                dag_chunk_size=None, degenerate_blocks=None):
        projected = self._projected_corrections(degenerate_blocks, degenerate_states,
                                               expansions, freqs, zero_cutoff, parallelizer,
                                               verbose=verbose)
        if projected is not None:
            if order is None: order = len(self.expansions) - 1
            if hasattr(projected,'overlaps'):
                return projected.overlaps(states,order,evaluation_mode=evaluation_mode,
                    dag_cache_size=dag_cache_size,dag_cache_bytes=dag_cache_bytes,
                    dag_chunk_size=dag_chunk_size)
            return np.stack(projected.diagonal(states, order, 'N', range(order+1)), axis=-1)
        if isinstance(degenerate_states, dict) and not self.solver.filtered_block_recurrence:
            results = [
                self.get_overlap_corrections(
                    [state], order=order, expansions=expansions,
                    degenerate_states=self._state_degeneracy_pairs(degenerate_states, state),
                    freqs=freqs, zero_cutoff=zero_cutoff, verbose=verbose,
                    parallelizer=parallelizer, evaluation_mode=evaluation_mode,
                    dag_cache_size=dag_cache_size, dag_cache_bytes=dag_cache_bytes,
                    dag_chunk_size=dag_chunk_size
                )
                for state in states
            ]
            if not results:
                raise ValueError("overlap corrections require at least one state")
            return np.concatenate(results, axis=0)
        expansions = self._prep_expansions(expansions)
        if freqs is None: freqs = self.freqs
        if order is None: order = len(expansions) - 1

        overlap_evaluators = [
            self.solver.overlap_correction(o)([]) for o in range(2, order + 1)
        ]

        # corrs = [
        #     evaluator.evaluate(states, expansions, freqs, verbose=verbose, log_scaled=False)
        #     for evaluator in overlap_evaluators
        # ]

        degenerate_changes = (None if self.solver.filtered_block_recurrence else
                              self.get_degenerate_changes(degenerate_states))
        degeneracy_identification_context = (
            None if degenerate_changes is None else
            DegeneracyIdentificationContext(degenerate_changes)
        )
        if parallelizer is None: parallelizer = self.parallelizer
        parallelizer = Parallelizer.lookup(parallelizer) if parallelizer is not None else None

        corrs = [
            np.concatenate(
                evaluator.evaluate(
                    [[s, [np.arange(len(s))]] for s in states],
                    expansions, freqs, verbose=verbose,
                    degenerate_changes=degenerate_changes,
                        projection_bindings=getattr(self, '_filtered_change_bindings', None), zero_cutoff=zero_cutoff, log_scaled=True,
                    parallelizer=parallelizer,
                    evaluation_mode=evaluation_mode,
                    dag_cache_size=dag_cache_size,
                    dag_cache_bytes=dag_cache_bytes,
                    dag_chunk_size=dag_chunk_size,
                    _degeneracy_identification_context=(
                        degeneracy_identification_context
                    )
                )
            )[:, np.newaxis]
            for evaluator in overlap_evaluators
        ]

        return np.concatenate([np.ones((len(states), 1)), np.zeros((len(states), 1))] + corrs, axis=1)

    @_timed_evaluation('Grouping state changes and permutations')
    def get_diff_map(self, state_map, logger=None, log_level=None):
        """Group endpoints by change; logging options control only the timer."""
        # states = []
        change_perms = {}
        for initial, finals in state_map:
            initial = np.asanyarray(initial)
            # states.append(initial)
            finals = np.asanyarray(finals)
            diffs = finals - initial[np.newaxis, :]
            perms = PerturbationTheoryTerm.change_sort(diffs)
            sort_diffs = nput.vector_take(diffs, perms, shared=1)
            # raise Exception(sort_diffs, perms)
            for ch, perm_blocks in zip(*nput.group_by(perms, sort_diffs)[0]):
                ch = tuple(ch[ch != 0])
                change_perms[ch] = change_perms.get(ch, [])
                change_perms[ch].append([initial, perm_blocks])

        # for ch, (state_block, pblock) in change_perms.items():
        #     pblock = np.unique(np.concatenate(pblock, axis=0), axis=0)
        #     change_perms[ch] = [np.array(state_block), pblock]

        return change_perms

    @staticmethod
    def get_finals(initial, change, perms):
        base_change = np.pad(change, [0, len(initial) - len(change)])
        changes = nput.vector_take(base_change, np.argsort(perms, axis=1))
        return np.asanyarray(initial)[np.newaxis] + changes

    def _reformat_corrections(self, order, corrs, change_map, num_expansions):
        # reshape to allow for better manipulation

        ndim = next(iter(change_map.values()))[0][0].shape[0]
        basis = HarmonicOscillatorProductBasis(ndim)
        change_states = {
            k:BasisStateSpace(basis, [s for s,p in state_perms])
            for k,state_perms in change_map.items()
        }
        change_final_states = {
            change: [self.get_finals(s, change, p) for s, p in state_perms]
            for change, state_perms in change_map.items()
        }
        total_basis = None
        for state_space in change_states.values():
            if total_basis is None:
                total_basis = state_space
            else:
                total_basis = total_basis.union(state_space)

        finals_map = [
            None for _ in range(len(total_basis))
        ]
        corrs_map = [
            [None]*(order+1)
            for _ in range(len(total_basis))
        ]


        for change, subcorrs in corrs.items():
            initial_states = change_states[change]
            final_states = change_final_states[change]
            basis_idx = total_basis.find(initial_states)
            for o in range(order+1):
                if subcorrs[o] is None:
                    # A targeted evaluation did not calculate this order. Keep
                    # the usual rectangular correction layout for callers.
                    continue
                for state_pos, finals, exp_corr_vals in zip(basis_idx, final_states, subcorrs[o]):
                    if corrs_map[state_pos][o] is None:
                        if num_expansions is None:
                            corrs_map[state_pos][o] = exp_corr_vals
                        else:
                            corrs_map[state_pos][o] = [cc for cc in exp_corr_vals]
                    else:
                        if num_expansions is None:
                            corrs_map[state_pos][o] = np.concatenate([corrs_map[state_pos][o], exp_corr_vals])
                        else:
                            corrs_map[state_pos][o] = [
                                np.concatenate([c1, corr])
                                for c1, corr in zip(corrs_map[state_pos][o], exp_corr_vals)
                            ]

            for state_pos, finals in zip(basis_idx, final_states):
                if finals_map[state_pos] is None:
                    finals_map[state_pos] = finals.astype(int)
                else:
                    finals_map[state_pos] = np.concatenate([finals_map[state_pos], finals.astype(int)], axis=0)

        if num_expansions is None:
            for finals, order_block in zip(finals_map, corrs_map):
                for o, vals in enumerate(order_block):
                    if vals is None:
                        order_block[o] = np.zeros(len(finals))

        # all_corrs = [
        #
        # ]
        # final_states = {
        #     state_space:self.get_finals(s, change, p)
        # }
        #
        # for exp_corrs in final_corrs:
        #     for n, cb in enumerate(exp_corrs):
        #         for o, c in enumerate(cb):
        #             cb[o] = np.concatenate(c, axis=1)
        #
        # for n, s in enumerate(final_states):
        #     final_states[n] = np.concatenate(s, axis=0)
        #

        if num_expansions is None:
            all_corrs = BasicAPTCorrections(total_basis, [BasisStateSpace(basis, f) for f in finals_map], corrs_map)
        else:
            # currently is initial x order x operator
            all_corrs = [
                BasicAPTCorrections(total_basis, [BasisStateSpace(basis, f) for f in finals_map],
                                    [[c[i] for c in ord_block] for ord_block in corrs_map])
                for i in range(num_expansions)
            ]

        return all_corrs
    @classmethod
    def _compute_corr(cls,  gen_corrs, only_degs, key, gen, expr, *args, only_degenerate_terms=None,
                      resonant_complement=False, **kwargs):
        if kwargs.pop('_filtered_block_recurrence', False):
            kwargs['degenerate_changes'] = None
            new_corrs = [-value for value in expr.evaluate(*args, only_degenerate_terms=False, **kwargs)]
            return new_corrs if gen_corrs is None else [a+b for a,b in zip(gen_corrs,new_corrs)]
        if resonant_complement:
            # For a != b the full (non-degenerate) main expression is
            # <psi_a|H|psi_b> = 0 at every order, so the sum with the resonant
            # terms excluded equals minus the sum over the resonant terms only.
            # Only valid for the Hamiltonian itself, off the diagonal, with the
            # same H in the operator and in the wavefunction corrections.
            new_corrs = [-c for c in expr.evaluate(*args, only_degenerate_terms=True, **kwargs)]
        else:
            new_corrs = expr.evaluate(*args, **kwargs)
        if only_degenerate_terms is True or only_degenerate_terms is False:
            only_degenerate_terms = [only_degenerate_terms, only_degenerate_terms]
        other_degs = only_degenerate_terms[0] is only_degenerate_terms[1]
        if other_degs and not only_degs:
            new_corrs = [-c for c in new_corrs]
        # print(f"|{key}", gen, "->", new_corrs)
        if gen_corrs is None:
            gen_corrs = new_corrs
        else:
            gen_corrs = [
                g + d
                for g, d in zip(gen_corrs, new_corrs)
            ]
        return gen_corrs
    @classmethod
    def _build_corrections(cls, corr_gen, degenerate_corr_gen, expansions, order,
                           terms, allowed_coefficients, disallowed_coefficients,
                           epaths,
                           change_map, degenerate_changes, only_degenerate_terms,
                           include_degenerate_correction_terms,
                           freqs, verbose, logger, zero_cutoff, log_scaled,
                           parallelizer, shared_coefficients=None,
                           degeneracy_identification_context=None,
                           target_orders=None,
                           resonant_complement=False, projection_bindings=None,
                           filtered_block_recurrence=False):
        corrs = {}
        if degenerate_changes is not None:
            changes = (
                list(degenerate_changes.left) + list(degenerate_changes.right)
                if isinstance(degenerate_changes, SideDegeneracyChanges) else
                degenerate_changes
            )
            allowed_degenerate_changes = list({PerturbationTheoryTerm.sorted_changes(c) for _,c in changes})
        else:
            allowed_degenerate_changes = None
        if filtered_block_recurrence:
            allowed_degenerate_changes = ()
        for o in range(order + 1):
            if target_orders is not None and o not in target_orders:
                for change in change_map:
                    corrs.setdefault(change, []).append(None)
                continue
            generator = corr_gen(o,
                                 allowed_terms=terms,
                                 allowed_energy_changes=epaths,
                                 allowed_coefficients=allowed_coefficients,
                                 disallowed_coefficients=disallowed_coefficients
                                 )
            if (degenerate_changes is not None or filtered_block_recurrence) and degenerate_corr_gen is not None:
                deg_gen = degenerate_corr_gen(
                    o,
                    allowed_terms=terms,
                    allowed_energy_changes=epaths,
                    allowed_coefficients=allowed_coefficients,
                    degenerate_changes=allowed_degenerate_changes,
                    disallowed_coefficients=disallowed_coefficients
                )
                if not isinstance(deg_gen, dict):
                    if hasattr(deg_gen, 'both'): # TODO: figure out why the corrections really use the neither form...
                        deg_gen = {
                            'left':deg_gen.left,
                            'right':deg_gen.right,
                            'both':deg_gen.both
                        }
                    else:
                        deg_gen = {
                            'both':deg_gen
                        }
                deg_left_gen = deg_gen.get('left')
                deg_right_gen = deg_gen.get('right')
                deg_both_gen = deg_gen.get('both')
                deg_neither_gen = deg_gen.get('neither')
            else:
                deg_left_gen = deg_right_gen = deg_both_gen = deg_neither_gen = None
            if include_degenerate_correction_terms is True:
                include_degenerate_correction_terms = ['left', 'right', 'both', 'neither']
            elif include_degenerate_correction_terms is False:
                include_degenerate_correction_terms = []
            if not isinstance(include_degenerate_correction_terms, dict):
                include_degenerate_correction_terms = {k:True for k in include_degenerate_correction_terms}
            include_degenerate_correction_terms = dict(
                {
                    'main':True,
                    'left':False,
                    'right':False,
                    'both':False,
                    'neither':False
                },
                **include_degenerate_correction_terms
            )
            if not include_degenerate_correction_terms['main']:
                main_gen = None
            else:
                main_gen = generator
            if not include_degenerate_correction_terms['left']:
                deg_left_gen = None
            if not include_degenerate_correction_terms['right']:
                deg_right_gen = None
            if not include_degenerate_correction_terms['both']:
                deg_both_gen = None
            if not include_degenerate_correction_terms['neither']:
                deg_neither_gen = None
            for change, state_perms in change_map.items():
                corr_block = corrs.get(change, [])
                corrs[change] = corr_block
                with logger.block(tag="Getting corrections at order {o}", o=o):
                    expr = None if main_gen is None else main_gen(change)
                    left_deg_expr = None if deg_left_gen is None else deg_left_gen(change)
                    right_deg_expr = None if deg_right_gen is None else deg_right_gen(change)
                    both_deg_expr = None if deg_both_gen is None else deg_both_gen(change)
                    neither_deg_expr = None if deg_neither_gen is None else deg_neither_gen(change)
                    gen_corrs = None
                    with logger.block(tag="evaluating..."):
                        start = time.time()
                        eval_args = (
                            state_perms,
                            expansions,
                            freqs
                        )
                        eval_kwargs = dict(
                            degenerate_changes=degenerate_changes,
                            projection_bindings=projection_bindings,
                            _filtered_block_recurrence=filtered_block_recurrence,
                            verbose=verbose,
                            zero_cutoff=zero_cutoff,
                            log_scaled=log_scaled,
                            parallelizer=parallelizer,
                            _shared_coefficients=shared_coefficients,
                            _degeneracy_identification_context=(
                                degeneracy_identification_context
                            )
                        )
                        if expr is not None:
                            gen_corrs = cls._compute_corr(
                                gen_corrs, only_degenerate_terms,
                                "!", main_gen, expr,
                                *eval_args,
                                only_degenerate_terms=only_degenerate_terms,
                                resonant_complement=(
                                    resonant_complement
                                    and degenerate_changes is not None
                                    and len(change) > 0
                                ),
                                **eval_kwargs
                            )
                        if left_deg_expr is not None:
                            gen_corrs = cls._compute_corr(
                                gen_corrs, only_degenerate_terms,
                                "L", deg_left_gen, left_deg_expr,
                                *eval_args,
                                only_degenerate_terms=[True, False],
                                **eval_kwargs
                            )
                        if right_deg_expr is not None:
                            gen_corrs = cls._compute_corr(
                                gen_corrs, only_degenerate_terms,
                                "R", deg_right_gen, right_deg_expr,
                                *eval_args,
                                only_degenerate_terms=[False, True],
                                **eval_kwargs
                            )
                        if both_deg_expr is not None:
                            gen_corrs = cls._compute_corr(
                                gen_corrs, only_degenerate_terms,
                                "B", deg_both_gen, both_deg_expr,
                                *eval_args,
                                only_degenerate_terms=[True, True],
                                **eval_kwargs
                            )
                        if neither_deg_expr is not None:
                            gen_corrs = cls._compute_corr(
                                gen_corrs, only_degenerate_terms,
                                "N", deg_neither_gen, neither_deg_expr,
                                *eval_args,
                                only_degenerate_terms=[False, False],
                                **eval_kwargs
                            )
                        end = time.time()
                        logger.log_print('took {e:.3f}s', e=end - start)
                corr_block.append(gen_corrs)
        return corrs

    @classmethod
    def _build_state_specific_corrections(cls, corr_gen, degenerate_corr_gen,
                                          expansions, order, terms,
                                          allowed_coefficients, disallowed_coefficients,
                                          epaths, change_map, degenerate_states,
                                          only_degenerate_terms,
                                          include_degenerate_correction_terms,
                                          freqs, verbose, logger, zero_cutoff,
                                          log_scaled, parallelizer,
                                          shared_coefficients=None,
                                          target_orders=None,
                                          resonant_complement=False, projection_bindings=None):
        """Evaluate one state/permutation row at a time and restore the old shape."""
        result = {
            change: [
                None if target_orders is not None and o not in target_orders
                else [None] * len(entries)
                for o in range(order + 1)
            ]
            for change, entries in change_map.items()
        }
        change_cache = {}
        for change, entries in change_map.items():
            for entry_index, (left_state, permutations) in enumerate(entries):
                row_values = [[] for _ in range(order + 1)]
                left_pairs = cls._state_degeneracy_pairs(degenerate_states, left_state)
                for permutation in permutations:
                    one_perm = np.asarray(permutation)[np.newaxis, :]
                    right_state = cls.get_finals(left_state, change, one_perm)[0]
                    right_pairs = cls._state_degeneracy_pairs(degenerate_states, right_state)
                    left_key = tuple(int(x) for x in left_state)
                    right_key = tuple(int(x) for x in right_state)
                    if left_key not in change_cache:
                        change_cache[left_key] = cls.get_degenerate_changes(left_pairs) or []
                    if right_key not in change_cache:
                        change_cache[right_key] = cls.get_degenerate_changes(right_pairs) or []
                    left_changes = change_cache[left_key]
                    right_changes = change_cache[right_key]
                    side_changes = (
                        SideDegeneracyChanges(left_changes, right_changes)
                        if left_changes or right_changes else None
                    )
                    one_result = cls._build_corrections(
                        corr_gen, degenerate_corr_gen, expansions, order,
                        terms, allowed_coefficients, disallowed_coefficients,
                        epaths, {change: [[left_state, one_perm]]},
                        side_changes, only_degenerate_terms,
                        include_degenerate_correction_terms,
                        freqs, verbose, logger, zero_cutoff, log_scaled,
                        parallelizer,
                        shared_coefficients=shared_coefficients,
                        target_orders=target_orders,
                        resonant_complement=resonant_complement, projection_bindings=projection_bindings
                    )[change]
                    for o in range(order + 1):
                        if result[change][o] is not None:
                            row_values[o].append(one_result[o][0])
                for o in range(order + 1):
                    if result[change][o] is not None:
                        result[change][o][entry_index] = np.concatenate(
                            row_values[o], axis=-1
                        )
        return result

    @classmethod
    def get_degenerate_changes(cls, degenerate_pairs):
        if degenerate_pairs is None: return None

        finals = []
        directed = isinstance(degenerate_pairs, StateDegeneracyPairs)
        for start,end in degenerate_pairs:
            orientations = [[end, start]] if directed else [[start, end], [end, start]]
            for e,s in orientations:
                diff = np.subtract(e, s)
                nzp = np.nonzero(diff)
                if len(nzp) > 0: nzp = nzp[0]
                if len(nzp) > 0:
                    finals.append([nzp, diff[nzp,]])
        if len(finals) == 0: finals = None

        return finals

    def _sort_corrections(self, corrs:BasicAPTCorrections, state_pairs):
        basis = corrs.initial_states.basis
        initial_basis = BasisStateSpace(basis, [s[0] for s in state_pairs])
        corr_sorting = np.argsort(initial_basis.find(corrs.initial_states))
        corr_data = list(zip(corrs.final_states, corrs.corrections))
        new_finals = []
        new_corrs = []
        for n,(finals, corr_block) in enumerate(corr_data[i] for i in corr_sorting):
            final_basis = BasisStateSpace(basis, state_pairs[n][1])
            subcorr_sorting = np.argsort(final_basis.find(finals, missing_val=len(final_basis)+1))
            corr_block = [c[subcorr_sorting,] for c in corr_block]
            finals = final_basis.take_subspace(subcorr_sorting)
            new_corrs.append(corr_block)
            new_finals.append(finals)

        return BasicAPTCorrections(initial_basis, new_finals, new_corrs)

    @_timed_evaluation('State-by-state corrections (total)')
    @_compiled_evaluation_context
    def get_state_by_state_corrections(self, generator, states, order=None,
                                       terms=None, epaths=None,
                                       expansions=None, freqs=None, verbose=False,
                                       allowed_coefficients=None, disallowed_coefficients=None,
                                       degenerate_states=None, only_degenerate_terms=False,
                                       degenerate_correction_generator=None,
                                       include_degenerate_correction_terms=True,
                                       log_scaled=False,
                                       zero_cutoff=None,
                                       return_sorted=False,
                                       logger=None,
                                       parallelizer=None,
                                       target_orders=None,
                                       resonant_complement=False, projection_bindings=None, evaluation_mode=None):
        PerturbationTheoryExpressionEvaluator._cached_expansion = None # to make paralellism work right

        spex = expansions is None or self.is_single_expansion(expansions)
        if spex:
            expansions = self._prep_expansions(expansions)
            if order is None: order = len(expansions) - 1
        else:
            expansions = [self._prep_expansions(e) for e in expansions]
            if order is None: order = len(expansions[0]) - 1
        if freqs is None: freqs = self.freqs
        if target_orders is not None:
            target_orders = frozenset(target_orders)
            if not target_orders or min(target_orders) < 0 or max(target_orders) > order:
                raise ValueError("target_orders must be nonempty orders between 0 and order")
            if not spex:
                raise ValueError("target_orders requires a single coefficient expansion")
        if logger is None: logger=self.solver.logger
        logger = Logger.lookup(logger)
        if parallelizer is None: parallelizer=self.parallelizer
        parallelizer = Parallelizer.lookup(parallelizer) if parallelizer is not None else None


        all_gs = nput.is_numeric(states[0][0])
        if all_gs:
            states = [
                [[0] * len(states[0]), states]
            ]
        change_map = self.get_diff_map(states, logger=logger)

        filtered_block = bool(self.solver.filtered_block_recurrence)
        if filtered_block and projection_bindings is None:
            projection_bindings = self._filtered_change_bindings
        state_specific = isinstance(degenerate_states, dict) and not filtered_block
        degenerate_changes = (
            None if state_specific or filtered_block else self.get_degenerate_changes(degenerate_states)
        )

        degeneracy_identification_context = (
            None if degenerate_changes is None else
            DegeneracyIdentificationContext(degenerate_changes)
        )

        shared_coefficients = None
        if isinstance(parallelizer, MultiprocessingParallelizer):
            coefficient_batch = expansions if not spex else [expansions]
            coefficient_batch = (
                PerturbationTheoryExpressionEvaluator
                ._normalize_coefficient_expansions(coefficient_batch)
            )
            shared_coefficients = SharedMemoryArrayTree.create(
                coefficient_batch
            )
            PerturbationTheoryExpressionEvaluator._cached_expansion = None
        try:
            build_corrections = (
                self._build_state_specific_corrections if state_specific else
                self._build_corrections
            )
            corrs = build_corrections(
                generator, degenerate_correction_generator, expansions, order,
                terms, allowed_coefficients, disallowed_coefficients,
                epaths, change_map,
                degenerate_states if state_specific else degenerate_changes,
                only_degenerate_terms,
                include_degenerate_correction_terms,
                freqs, verbose, logger, zero_cutoff, log_scaled, parallelizer,
                shared_coefficients=shared_coefficients,
                target_orders=target_orders,
                resonant_complement=resonant_complement,
                projection_bindings=projection_bindings,
                **({} if state_specific else {
                    'degeneracy_identification_context': degeneracy_identification_context,
                    'filtered_block_recurrence': filtered_block
                })
            )
        finally:
            if shared_coefficients is not None:
                PerturbationTheoryExpressionEvaluator._cached_expansion = None
                shared_coefficients.dispose()
        with _log_evaluation_time(logger, "Assembling and sorting correction results"):
            new_corrs = self._reformat_corrections(order, corrs, change_map, None if spex else len(expansions))
            if return_sorted:
                if isinstance(new_corrs, BasicAPTCorrections):
                    new_corrs = self._sort_corrections(new_corrs, states)
                else:
                    new_corrs = [self._sort_corrections(nc, states) for nc in new_corrs]
        return new_corrs

    @_timed_evaluation('Matrix corrections (total)')
    def get_matrix_corrections(self, states, order=None, expansions=None, freqs=None,
                               zero_cutoff=None, verbose=False):
        gen = lambda o, **kw: self.solver.hamiltonian_expansion[o]
        return self.get_state_by_state_corrections(gen, states, order=order, expansions=expansions,
                                                   freqs=freqs, zero_cutoff=zero_cutoff, verbose=verbose)
    @_timed_evaluation('Full wavefunction corrections (total)')
    def get_full_wavefunction_corrections(self, states, order=None, expansions=None, freqs=None,
                                          zero_cutoff=None, degenerate_states=None, verbose=False,
                                          degenerate_blocks=None):
        projected = self._projected_corrections(degenerate_blocks, degenerate_states,
                                               expansions, freqs, zero_cutoff, verbose=verbose)
        if projected is not None:
            return projected.wavefunctions(states, len(self.expansions)-1 if order is None else order)
        return self.get_state_by_state_corrections(self.solver.full_wavefunction_correction, states, degenerate_states=degenerate_states,
                                                   order=order, expansions=expansions, freqs=freqs, zero_cutoff=zero_cutoff,
                                                   degenerate_correction_generator=lambda *args,**kwargs:{
                                                       **({} if self.solver.filtered_block_recurrence else {'neither':self.solver.overlap_correction(*args,**kwargs)})
                                                   },
                                                   verbose=verbose
                                                   )
    @_timed_evaluation('Wavefunction corrections (total)')
    def get_wavefunction_corrections(self, states, order=None, expansions=None, freqs=None,
                                          zero_cutoff=None, degenerate_states=None, verbose=False,
                                          target_orders=None, degenerate_blocks=None):
        projected = self._projected_corrections(degenerate_blocks, degenerate_states,
                                               expansions, freqs, zero_cutoff, verbose=verbose)
        if projected is not None:
            return projected.wavefunctions(states, len(self.expansions)-1 if order is None else order,
                                            full=False, target_orders=target_orders)
        return self.get_state_by_state_corrections(self.solver.wavefunction_correction, states, degenerate_states=degenerate_states,
                                                   order=order, expansions=expansions, freqs=freqs, zero_cutoff=zero_cutoff,
                                                   verbose=verbose, target_orders=target_orders)

    @classmethod
    def _prep_ham_corrs(cls, base_corrs):
        full_corrs = []
        exps = []
        for term_list_exps in base_corrs:
            if nput.is_zero(term_list_exps):
                full_corrs.append([])
                exps.append([])
            else:
                full_corrs.append(
                    [term_list for term_list, exp in term_list_exps]
                )
                exps.append(
                    [exp for term_list, exp in term_list_exps]
                )
        return full_corrs, exps

    @_timed_evaluation('Hamiltonian operator corrections (optional)')
    @_compiled_evaluation_context
    def get_hamiltonian_operator_corrections(self, states, order=None, expansions=None,
            freqs=None, degenerate_states=None, degenerate_blocks=None, terms=None,
            verbose=False, **opts):
        """Full ordered Hamiltonian matrix elements in the native operator-correction sign.

        Unlike the effective-Hamiltonian triangular assembler, this route evaluates
        both directions independently, including across distinct resolved blocks.
        No model-space grouping or numerical cutoff is changed.
        """
        if opts.pop('only_degenerate_terms',False):
            raise ValueError('Hamiltonian operator contributions require the full operator')
        if not opts.pop('include_diagonal',True):
            raise ValueError('Hamiltonian operator contributions include diagonal elements')
        if order is None:order=len(self.expansions)-1
        projected=self._projected_corrections(degenerate_blocks,degenerate_states,
            expansions,freqs,opts.get('zero_cutoff'),opts.get('parallelizer'),
            logger=opts.get('logger'),verbose=verbose)
        if projected is not None:
            if hasattr(projected,'validate_options'):projected.validate_options(opts)
            else:self._validate_projected_options(opts)
            # ConditionalCorrectionEvaluator prepares raw coefficients itself.
            # A second preparation would transpose internal kinetic tensors twice.
            coefficients=expansions if hasattr(projected,'energies') else self._prep_expansions(expansions)
            return projected.operators(states,order,self.solver.hamiltonian_expansion,
                coefficients,terms=terms,
                branches=opts.get('include_degenerate_correction_terms',True),
                target_orders=opts.get('target_orders'))
        return self.get_reexpressed_hamiltonian(states,order=order,expansions=expansions,
            freqs=freqs,degenerate_states=degenerate_states,degenerate_blocks=degenerate_blocks,
            only_degenerate_terms=False,include_diagonal=True,terms=terms,
            verbose=verbose,**opts)

    @_timed_evaluation('Effective Hamiltonian corrections (total)')
    @_compiled_evaluation_context
    def get_reexpressed_hamiltonian(self, states, order=None, expansions=None, freqs=None,
                                    degenerate_states=None, only_degenerate_terms=False,
                                    verbose=False, include_diagonal=False,
                                    hamiltonian_corrections=None,
                                    degenerate_blocks=None,
                                    evaluation_mode=None, **opts):
        """

        :param states:
        :param order:
        :param expansions:
        :param freqs:
        :param degenerate_states:
        :param only_degenerate_terms:
        :param verbose:
        :param include_diagonal:
        :param hamiltonian_corrections:  `[[(order, terms), expansion], ...]`
        :param opts:
        :return:
        """
        normalization=opts.pop('normalization','orthonormal')
        if normalization not in {'orthonormal','intermediate'}:
            raise ValueError('unknown block normalization '+str(normalization))
        projected = self._projected_corrections(degenerate_blocks, degenerate_states,
                                               expansions, freqs, opts.get('zero_cutoff'), opts.get('parallelizer'),
                                               logger=opts.get('logger'), verbose=verbose)
        if projected is not None:
            if hasattr(projected,'validate_options'):
                projected.validate_options(opts)
            else:
                self._validate_projected_options(opts)
            if hamiltonian_corrections is not None:
                projected = self.modify_hamiltonian(hamiltonian_corrections)._projected_corrections(
                    degenerate_blocks, degenerate_states, expansions, freqs, opts.get('zero_cutoff'),opts.get('parallelizer'),
                    logger=opts.get('logger'), verbose=verbose)
                if hasattr(projected,'validate_options'):
                    projected.validate_options(opts)
            if order is None: order = len(self.expansions) - 1
            if nput.is_numeric(states[0][0]):
                states = [[states[i], states[i+(0 if include_diagonal else 1):]]
                          for i in range(len(states)-(0 if include_diagonal else 1))]
            symbolic_block = (isinstance(projected, ConditionalCorrectionEvaluator)
                              and projected.stack.block_recurrence)
            if normalization=='intermediate' and hasattr(projected,'validate_options') and not symbolic_block:
                raise NotImplementedError('direct Bloch blocks require explicit resolved model-space blocks')
            block_options=({'normalization':normalization} if symbolic_block or not hasattr(projected,'validate_options') else {})
            return projected.hamiltonian(states, order, only_degenerate_terms,
                terms=opts.get('terms'), branches=opts.get('include_degenerate_correction_terms',True),
                target_orders=opts.get('target_orders'),**block_options)
        if self.solver.filtered_block_recurrence and opts.get('terms') is None and opts.get('include_degenerate_correction_terms', True) is True:
            if hamiltonian_corrections is not None:
                return self.modify_hamiltonian(hamiltonian_corrections).get_reexpressed_hamiltonian(
                    states, order=order, expansions=expansions, freqs=freqs,
                    degenerate_states=degenerate_states, degenerate_blocks=degenerate_blocks,
                    only_degenerate_terms=only_degenerate_terms, include_diagonal=include_diagonal,
                    normalization=normalization, verbose=verbose, **opts)
            intermediate = normalization == 'intermediate'
            if bool(self.solver.intermediate_normalization) != intermediate:
                evaluator = copy.copy(self)
                evaluator.solver = copy.copy(self.solver)
                evaluator.solver.intermediate_normalization = intermediate
                evaluator.solver._op_maps = {}
                return evaluator.get_reexpressed_hamiltonian(states, order=order,
                    expansions=expansions, freqs=freqs, degenerate_states=degenerate_states,
                    degenerate_blocks=degenerate_blocks, only_degenerate_terms=only_degenerate_terms,
                    include_diagonal=include_diagonal, normalization=normalization,
                    verbose=verbose, **opts)
            if order is None: order = len(self.expansions)-1
            if nput.is_numeric(states[0][0]):
                states = [[states[i], states[i+(0 if include_diagonal else 1):]]
                          for i in range(len(states)-(0 if include_diagonal else 1))]
            sign = -1 if only_degenerate_terms else 1
            return self.get_state_by_state_corrections(
                lambda n, **options: CorrectionView(sign*self.solver.block_hamiltonian_correction(n), **options),
                states, order=order, expansions=expansions, freqs=freqs,
                projection_bindings=self._filtered_change_bindings,
                verbose=verbose, **{k:v for k,v in opts.items() if k!='include_degenerate_correction_terms'})
        if normalization=='intermediate':
            raise ValueError('direct Bloch blocks require explicit degenerate_blocks')
        if freqs is None: freqs = self.freqs
        if hamiltonian_corrections is not None:
            hamiltonian_corrections, expansion_corrections = self._prep_ham_corrs(hamiltonian_corrections)
            extra_expansion = expansion_corrections
        else:
            extra_expansion = [np.diag(freqs)]

        if expansions is None:
            exps = self._prep_operator_expansion(expansions, extra_expansion)
        elif self.is_single_expansion(expansions):
            exps = self._prep_operator_expansion(expansions, extra_expansion)
            if order is None: order = len(expansions) - 1
        else:
            if order is None: order = len(expansions[0]) - 1
            exps = [
                self._prep_operator_expansion(exp, extra_expansion)
                for exp in expansions
            ]

        utri_block = nput.is_numeric(states[0][0])
        if utri_block:
            states = [
                [states[i], states[i+(0 if include_diagonal else 1):]]
                for i in range(len(states) - (0 if include_diagonal else 1))
            ]

        if hamiltonian_corrections is None:
            corr_gen=self.solver.reexpressed_hamiltonian
            deg_corr_gen=self.solver.reexpressed_hamiltonian_degenerate_correction
        else:
            corr_gen=lambda *a,**kw:self.solver.reexpressed_hamiltonian(
                *a, hamiltonian_corrections=hamiltonian_corrections, **kw
            )
            deg_corr_gen= lambda order,**kw:self.solver.reexpressed_hamiltonian_degenerate_correction(
                order, hamiltonian_corrections=hamiltonian_corrections, **kw
            )

        cur_symm = PTTensorCoeffProductSum.default_symmetrizers
        try:
            return self.get_state_by_state_corrections(corr_gen, states,
                                                       degenerate_states=degenerate_states,
                                                       only_degenerate_terms=only_degenerate_terms,
                                                       resonant_complement=(
                                                           self.use_resonant_complement
                                                           and only_degenerate_terms is True
                                                           and hamiltonian_corrections is None
                                                       ),
                                                       degenerate_correction_generator=deg_corr_gen,
                                                       order=order, expansions=exps, freqs=freqs, verbose=verbose,
                                                       log_scaled=True,
                                                       **opts)
        finally:
            PTTensorCoeffProductSum.default_symmetrizers = cur_symm

    def _prep_operator_expansion(self, expansions, operator_expansion):
        if expansions is None: expansions = self.expansions
        # expansions = self._prep_expansions(expansions)
        # if order is None: order = len(operator_expansion) - 1

        max_len = max(len(e) for e in expansions)

        exps = []
        for base, op in zip(expansions, operator_expansion):
            # we always make the operator be the 5th element
            if nput.is_numeric_array_like(op):
                op = np.asanyarray(op)
                if len(np.unique(op.shape)) < 2:
                    op = [op]
                else:
                    op = list(op)
            elif nput.is_numeric(op):
                op = [op]
            exps.append(
                base + [0] * (max_len - len(base)) + op
            )
        if len(operator_expansion) > len(expansions):
            for op in operator_expansion[len(expansions):]:
                if nput.is_numeric_array_like(op):
                    op = np.asanyarray(op)
                    if len(np.unique(op.shape)) < 2:
                        op = [op]
                    else:
                        op = list(op)
                elif nput.is_numeric(op):
                    op = [op]
                exps.append(
                    [0] * (max_len) + op
                )
        elif len(expansions) > len(operator_expansion):
            for base in expansions[len(operator_expansion):]:
                exps.append(base + [0] * (max_len+1 - len(base)))

        max_len = max(len(e) for e in exps)
        exps = [
            e + [0] * (max_len - len(e))
            for e in exps
        ]

        return exps

    @_timed_evaluation('Operator corrections (total)')
    @_compiled_evaluation_context
    def get_operator_corrections(self, operator_expansion, states,
                                 order=None, expansions=None, freqs=None,
                                 degenerate_states=None,
                                 operator_type=None,
                                 check_single=True,
                                 terms=None, min_order=1, verbose=False,
                                 degenerate_blocks=None, operator_action_expansion=False,
                                 evaluation_mode=None, **opts):
        logger = Logger.lookup(self.solver.logger if opts.get('logger') is None else opts['logger'])
        with _log_evaluation_time(logger, "Preparing operator coefficients"):
            if check_single and self.is_single_expansion(operator_expansion, min_order=min_order):
                exps = self._prep_operator_expansion(expansions, operator_expansion)
                if order is None: order = len(operator_expansion) - 1
                base_index = len(exps[-1]) - 1
            else:
                if order is None: order = len(operator_expansion[0]) - 1
                exps = [
                    self._prep_operator_expansion(expansions, o)
                    for o in operator_expansion
                ]
                base_index = len(exps[0][-1]) - 1

        if not isinstance(operator_action_expansion,(bool,np.bool_)):
            raise TypeError('operator_action_expansion must be boolean')
        from .Evaluators import _current_evaluation_mode
        action_token = _current_evaluation_mode.set(
            'actions' if operator_action_expansion else _current_evaluation_mode.get())
        cur_symm = PTTensorCoeffProductSum.default_symmetrizers
        symms = PTTensorCoeffProductSum.symmetrizers()
        if len(symms) < base_index + 1:
            symms = symms + (PTTensorCoeffProductSum._potential_symmetrizer,) * (base_index + 1 - len(symms))
            PTTensorCoeffProductSum.default_symmetrizers = symms
        try:
            if operator_type is None:
                operator_type = min_order
            projected = self._projected_corrections(degenerate_blocks, degenerate_states,
                                                   expansions, freqs, opts.get('zero_cutoff'),opts.get('parallelizer'),
                                                   logger=opts.get('logger'), verbose=verbose)
            if projected is not None:
                if operator_action_expansion:
                    raise ValueError('operator action expansion requires the native filtered recurrence')
                if hasattr(projected,'validate_options'):
                    projected.validate_options(opts)
                else:
                    self._validate_projected_options(opts)
                if opts.get('only_degenerate_terms',False):
                    raise ValueError('physical block properties exclude all model-space resolvent denominators')
                operators = self.solver.operator_expansion_terms(order,
                    operator_type=operator_type, base_index=base_index,
                    polynomial_representation=self.solver.polynomial_representation)
                single = self.is_single_expansion(exps)
                batches = [exps] if single else exps
                if hasattr(projected,'energies'):
                    return projected.operators(states,order,operators,exps,terms=terms,
                        branches=opts.get('include_degenerate_correction_terms',True),
                        target_orders=opts.get('target_orders'))
                prepared = [self._prep_expansions(coefficients) for coefficients in batches]
                return projected.operators(states, order, operators,
                    prepared[0] if single else prepared,
                    terms=terms, branches=opts.get('include_degenerate_correction_terms',True),
                    target_orders=opts.get('target_orders'))
            return self.get_state_by_state_corrections(
                lambda *a, **kw: self.solver.operator_correction(*a, operator_type=operator_type, base_index=base_index, **kw),
                states, order=order,
                expansions=exps, freqs=freqs, terms=terms, degenerate_states=degenerate_states,
                degenerate_correction_generator=(
                    lambda *a, **kw: self.solver.operator_degenerate_correction(*a, operator_type=operator_type, base_index=base_index, **kw)
                ),
                verbose=verbose, **opts)
        finally:
            _current_evaluation_mode.reset(action_token)
            PTTensorCoeffProductSum.default_symmetrizers = cur_symm

    @_timed_evaluation('Expression corrections (total)')
    def evaluate_expressions(self, states, exprs, expansions=None, operator_expansions=None,
                             degenerate_states=None, zero_cutoff=None,
                             verbose=False):
        if self.solver.resonance_condition is not None:
            # Conditional expressions already carry their endpoint restrictions.
            degenerate_states = None
        order = len(exprs) - 1
        expansions = (
            self._prep_operator_expansion(expansions, operator_expansions)
                if operator_expansions is not None else
            expansions #self._prep_expansions(expansions)
        )
        return self.get_state_by_state_corrections(
            lambda i, **kw: exprs[i],
            states,
            expansions=expansions,
            degenerate_states=degenerate_states,
            order=order,
            verbose=verbose,
            zero_cutoff=zero_cutoff
        )
