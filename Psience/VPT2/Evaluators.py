"""
Numerical evaluation backends for the analytic VPT expressions in `Analytic`.
"""
# imports from merged modules
import collections
import contextlib
import os
import sys
import time
import numpy as np
import contextvars
import dataclasses
import functools
import hashlib
import inspect
import itertools
import json
import marshal
import math
from pathlib import Path
import pickle
import uuid
import weakref
from McUtils.Parallelizers import Parallelizer, MultiprocessingParallelizer
from McUtils.Devutils import atomic_write, atomic_write_json, file_hash, freeze, fingerprint, implementation_signature, EntryStore
from .Analytic import ConditionalTransformationStack, SameStateCondition, operator_sum, operator_product, ZeroTerm, projection_inputs
from .Analytic import ShiftedEnergyBaseline
from .Common import Diagnostics, BoundedCache
import types

__all__ = ['AnalyticExpressionStore', 'CompiledAnalyticExpression', 'CompiledAnalyticEvaluation', 'AnalyticCorrectionCompiler']


# ---------------------------------------------------------------------------------------------
#  BoundedStream (formerly VPT2/BoundedStream.py)
#  Opt-in storage policy for the native analytic tensor-DAG evaluator.
#
#  This module does not change any symbolic products, conditions, or denominators.
#  Byte thresholds are checked periodically and describe the reachable Python
#  payload, not a hard bound on process RSS or a kernel's transient workspace.
# ---------------------------------------------------------------------------------------------

def retained_size(root, stop_after=None):
    """Count shared reachable payload once; do not visit classes or weak caches.

    Native symbolic objects expose their instance dictionaries. Container
    storage and ndarray owning buffers are charged without double counting
    shared objects. A view retains its base, which is therefore also visited.
    The accounting walk itself is temporary and does not retain old batches.
    """
    todo = [root]
    # Reflection is invariant for a type throughout this walk. In large
    # PolyPath DAGs, rediscovering every class's slots dominates bookkeeping.
    layouts = {}
    getsizeof = sys.getsizeof
    seen = set()
    size = 0
    while todo:
        obj = todo.pop()
        identity = id(obj)
        if identity in seen:
            continue
        seen.add(identity)
        size += getsizeof(obj)
        if stop_after is not None and size > stop_after:
            return size
        if isinstance(obj, np.ndarray):
            if obj.base is not None:
                todo.append(obj.base)
            if obj.dtype.hasobject:
                todo.extend(obj.flat)
        elif isinstance(obj, dict):
            todo.extend(obj.keys())
            todo.extend(obj.values())
        elif isinstance(obj, (tuple, list, set, frozenset)):
            todo.extend(obj)
        else:
            cls = type(obj)
            if cls not in layouts:
                slots = []
                if cls.__module__.startswith('Psience.VPT2.'):
                    for base in cls.__mro__:
                        names = getattr(base, '__slots__', ())
                        if isinstance(names, str):
                            names = (names,)
                        slots.extend(name for name in names
                                     if name not in {'__weakref__', '__dict__'})
                    layouts[cls] = tuple(slots)
                else:
                    layouts[cls] = None
            slots = layouts[cls]
            if slots is not None:
                attrs = getattr(obj, '__dict__', None)
                if attrs is not None:
                    todo.append(attrs)
                for name in slots:
                    value = getattr(obj, name, None)
                    if value is not None:
                        todo.append(value)
    return size


def trim_to_retained_size(cache):
    """Measure a cache's reachable payload once (at a batch boundary) and evict the
    oldest half until it fits `cache.max_bytes`."""
    measured = retained_size(cache._data)
    while cache.max_bytes is not None and measured > cache.max_bytes and len(cache):
        for _ in range(max(1, len(cache) // 2)):
            cache.popitem()
        measured = retained_size(cache._data)
    return measured


class StreamPreparationContext:
    """Bound immutable preparation to an evaluation plan, with bounded values.

    Child domains separate coefficient sets, including request-specific scalar
    values. Pure polynomial values, tensor gathers and frequency contractions
    can be reused across those domains. Coverage counts and output accumulators
    belong to each expression evaluation and are deliberately not stored here.
    """
    def __init__(self, coefficients, max_items=32, *, immutable_expressions=False,
                 evaluation_plan=None, options=None):
        self.coefficients = coefficients
        self.immutable_expressions = bool(immutable_expressions)
        self.max_items = max(0, int(max_items))
        self.data = collections.OrderedDict()
        self.plans = collections.OrderedDict()
        self.data_hits = self.data_misses = 0
        self.plan_hits = self.plan_misses = 0
        self.evaluation_plan = evaluation_plan
        self.options = {} if options is None else dict(options)
        self.pid = os.getpid()
        self.domains = collections.OrderedDict()
        self.bindings = {}
        self.copied_inputs = {}
        self.normalized_inputs = {}
        self.borrowed_handles = {}
        self.parallelizers = {}
        self.shared_arrays = {}
        self.numeric = {}
        self.closed = False

    @property
    def root(self):
        return self if self.evaluation_plan is None else self.evaluation_plan

    def bind(self, value):
        """Bind an input array once; identity cannot be recycled during the plan."""
        if isinstance(value, np.ndarray) and value.ndim:
            root = self.root
            identity = id(value)
            entry = root.bindings.get(identity)
            if entry is None:
                # Bind layout once. Strong input ownership prevents identity or
                # buffer-address recycling; equivalent views still get equal keys.
                key = ('array',(value.__array_interface__['data'][0],value.dtype.str,
                                value.shape,value.strides))
                entry = (value,key)
                root.bindings[identity] = entry
            return entry[1]
        if isinstance(value, (tuple, list)):
            return tuple(self.bind(part) for part in value)
        if isinstance(value, np.ndarray):
            return ('scalar', value.dtype.str, value.tobytes())
        return ('scalar', value)

    def domain(self, coefficients):
        root = self.root
        key = root.bind(coefficients)
        entry = root.domains.get(key)
        if entry is None:
            entry = type(self)(coefficients, self.max_items,
                               immutable_expressions=True, evaluation_plan=root)
            root._store(root.domains, key, entry)
        return entry

    def cache(self, name, max_bytes):
        from .Analytic import PerturbationTheoryExpressionEvaluator as Expression
        root = self.root
        if name not in root.numeric:
            items = root.options.get('path_cache_size',10000)
            if items is None:
                items = Expression.default_materialized_path_cache_size
            root.numeric[name] = BoundedCache.sampled(items,max_bytes)
        return root.numeric[name]

    def polynomial_cache(self):
        from .Analytic import (_MaterializedEvaluationCache,
                               PerturbationTheoryExpressionEvaluator as Expression)
        root = self.root
        if 'polynomials' not in root.numeric:
            byte_limit = root.options.get('path_cache_bytes')
            paths = root.cache('polynomial_paths',64*1024**2 if byte_limit is None else byte_limit)
            # Preserve the native exact-value cache and its existing capacity /
            # persistence controls. path_cache_* governs DAG/path work, not the
            # exact polynomial cache; putting both under that smaller budget
            # evicts useful exact values during discovery.
            root.numeric['polynomials'] = _MaterializedEvaluationCache(
                Expression._poly_cache, batch_cache=paths)
        else:
            root.numeric['polynomials'].exact_cache = Expression._poly_cache
        return root.numeric['polynomials']

    def energy_contractions(self, frequencies, deltas):
        """Cache only immutable frequency/delta sums, never projector masks."""
        from .Analytic import _HashableNumPyIndex
        root = self.root
        cache = root.cache('energy_sums', 32*1024**2)
        keys = [(_HashableNumPyIndex(f), _HashableNumPyIndex(d))
                for f,d in zip(frequencies,deltas)]
        ready, missing = {}, {}
        for index,key in enumerate(keys):
            value = cache.get(key)
            if value is not None:
                ready[key] = value
            elif key not in missing:
                missing[key] = index
        if missing:
            indices = list(missing.values())
            values = np.einsum('rpm,rfm->rpf', frequencies[indices], deltas[indices], optimize=False)
            for key,value in zip(missing, values):
                value = value.copy()
                value.setflags(write=False)
                ready[key] = value
                cache[key] = value
        return np.stack([ready[key] for key in keys])

    def share_coefficients(self, coefficients):
        """One native shared arena per physical tensor, reused across requests."""
        from McUtils.Parallelizers import SharedMemoryArrayTree
        root = self.root
        def share(value):
            if not value.ndim:
                return value
            key = root.bind(value)
            if key not in root.shared_arrays:
                root.shared_arrays[key] = SharedMemoryArrayTree.create(value)
            return root.shared_arrays[key]
        return SharedMemoryArrayTree._map_tree(coefficients, share)

    def close(self):
        root = self.root
        for domain in root.domains.values():
            domain.data.clear()
            domain.plans.clear()
            domain.coefficients = None
            domain.closed = True
        root.domains.clear()
        root.data.clear()
        root.plans.clear()
        root.numeric.clear()
        root.copied_inputs.clear()
        root.normalized_inputs.clear()
        root.bindings.clear()
        root.borrowed_handles.clear()
        try:
            with contextlib.ExitStack() as cleanup:
                for owner in root.shared_arrays.values():
                    cleanup.callback(owner.dispose)
                root.shared_arrays.clear()
        finally:
            root.parallelizers.clear()
            root.closed = True

    def _store(self, cache, key, value):
        if self.max_items:
            cache[key] = value
            while len(cache) > self.max_items:
                cache.popitem(last=False)
        return value

    def prepare(self, evaluator, coefficient_sets):
        key = tuple(coefficient_sets)
        if key in self.data:
            if Diagnostics.enabled:
                self.data_hits += 1
            return self.data[key]
        if Diagnostics.enabled:
            self.data_misses += 1
        return self._store(self.data, key,
                           evaluator._prep_coeff_data(coefficient_sets, self.coefficients))

    def coefficient_plan(self, specs, tensors, create):
        key = (id(specs), id(tensors))
        entry = self.plans.get(key)
        if entry is not None:
            if Diagnostics.enabled:
                self.plan_hits += 1
            return entry[2], True
        if Diagnostics.enabled:
            self.plan_misses += 1
        plan = create()
        self._store(self.plans, key, (specs, tensors, plan))
        return plan, False

    def stats(self):
        return dict(stream_preparation_hits=self.data_hits,
                    stream_preparation_misses=self.data_misses,
                    stream_coefficient_plan_hits=self.plan_hits,
                    stream_coefficient_plan_misses=self.plan_misses,
                    stream_preparation_entries=len(self.data),
                    stream_coefficient_plan_entries=len(self.plans))


class StreamBudget:
    """Independent input-count and periodically checked reachable-byte limits."""
    def __init__(self, options):
        options = dict(options)
        self.max_terms = max(1, int(options.pop('max_terms', 65536)))
        byte_limit = options.pop('max_bytes', None)
        self.max_bytes = None if byte_limit is None else max(1, int(byte_limit))
        self.check_interval = max(1, int(options.pop('check_interval', 128)))
        self.callback = options.pop('callback', None)
        self.bounded_byte_accounting = bool(options.pop('bounded_byte_accounting', True))
        self.preparation_entries = max(0, int(options.pop('preparation_entries', 32)))
        self.reuse_preparation = bool(options.pop('reuse_preparation', True))
        self.account_cache_graphs = bool(options.pop('account_cache_graphs', True))
        if options:
            raise TypeError('unknown stream options: {}'.format(tuple(options)))
        self.terms = self.total_terms = self.maximum_terms = 0
        self.bytes = self.maximum_bytes = 0
        self.scan_seconds = 0.
        self.flush_reasons = collections.Counter()
        self.batch = 0

    def incoming(self, pending):
        self.terms += 1
        self.total_terms += 1
        self.maximum_terms = max(self.maximum_terms, self.terms)
        if self.max_bytes is not None and (
                self.terms == 1 or self.terms % self.check_interval == 0
                or self.terms >= self.max_terms):
            started = time.perf_counter()
            # A flush decision only needs to know whether the limit was
            # crossed. Preserve the same checkpoints and flush decisions.
            self.bytes = retained_size(pending, stop_after=(
                self.max_bytes if self.bounded_byte_accounting else None))
            self.scan_seconds += time.perf_counter() - started
            self.maximum_bytes = max(self.maximum_bytes, self.bytes)
        if self.max_bytes is not None and self.bytes >= self.max_bytes:
            return 'bytes'
        if self.terms >= self.max_terms:
            return 'terms'
        return None

    def event(self, phase, reason, keys, **extra):
        data = dict(phase=phase, batch=self.batch, reason=reason,
                    input_terms=self.terms, total_input_terms=self.total_terms,
                    coefficient_keys=keys, measured_pending_bytes=self.bytes,
                    pending_bytes_capped=bool(self.bounded_byte_accounting and
                                             self.max_bytes is not None and self.bytes > self.max_bytes),
                    **extra)
        if self.callback is not None:
            self.callback(data)

    def flushed(self, reason):
        self.flush_reasons[reason] += 1
        self.batch += 1
        self.terms = self.bytes = 0

    def stats(self):
        return dict(stream_max_input_terms=self.max_terms,
                    stream_max_pending_bytes=self.max_bytes,
                    stream_byte_check_interval=self.check_interval,
                    stream_input_terms=self.total_terms,
                    stream_bounded_byte_accounting=self.bounded_byte_accounting,
                    stream_peak_bytes_is_lower_bound=bool(self.bounded_byte_accounting and
                        self.max_bytes is not None and self.maximum_bytes > self.max_bytes),
                    stream_peak_batch_terms=self.maximum_terms,
                    stream_peak_measured_pending_bytes=self.maximum_bytes,
                    stream_accounting_seconds=self.scan_seconds,
                    stream_flush_reasons=dict(self.flush_reasons))


# ---------------------------------------------------------------------------------------------
#  PersistentCaches (formerly VPT2/PersistentCaches.py)
#  An optional on-disk layer under the VPT2 exact-result cache (``PerturbationTheoryExpressionEvaluator._poly_cache``).
#
#  That cache maps ``(polynomial, change, baseline_shift)`` to a dict of ``{state index: value}``. The
#  persistent layer keys each of those dicts by the same expansion fingerprint `CompiledExpressions` uses
#  (definition + implementation signature), keeps the in-memory cache in front, loads a polynomial's values
#  the first time it is looked up, and writes the dicts that grew back to an `ObjectStore` on `flush()`.
#  Polynomials that can't be fingerprinted, and value dicts whose keys aren't uniform integer state indices,
#  are simply not persisted.
# ---------------------------------------------------------------------------------------------

_NAMESPACE = 'vpt2-poly-values'


class PersistentEvaluationCache(collections.abc.MutableMapping):

    def __init__(self, store, memory=None, fingerprint=None, signature=None):
        """
        :param store: an `ObjectStore` (or a path to one)
        :param memory: the in-memory mapping that stays in front (default: a dict)
        :param fingerprint: ``fingerprint(poly_key) -> str`` (default: the `CompiledExpressions` one)
        :param signature: the implementation signature folded into every key (default: McUtils + Psience)
        """
        from McUtils.Scaffolding import ObjectStore
        if not isinstance(store, ObjectStore):
            store = ObjectStore(store)
        self.store = store
        self.memory = {} if memory is None else memory
        if fingerprint is None or signature is None:
            fingerprint = fingerprint or _fingerprint
            signature = signature if signature is not None else _implementation_signature()
        self._fingerprint = fingerprint
        self._signature = signature
        self._digests = {}
        self._persisted = {}  # poly_key -> number of entries already on disk
        self._touched = {}    # poly_keys handed out or set (the in-memory cache need not be iterable)
        self.disk_hits = 0
        self.disk_writes = 0

    # ---- identity ---- #
    def _store_key(self, poly_key):
        try:
            digest = self._digests[poly_key]
        except KeyError:
            try:
                digest = self._fingerprint(poly_key)
            except Exception:  # not a definition we can identify: memory only
                digest = None
            self._digests[poly_key] = digest
        return None if digest is None else self.store.key_for(_NAMESPACE, self._signature, digest)

    # ---- (de)serialization of one {state: value} dict ---- #
    @staticmethod
    def _encode(values):
        from .Analytic import _HashableNumPyIndex
        if not values or not all(isinstance(k, _HashableNumPyIndex) for k in values):
            return None
        first = next(iter(values))
        if any(k.shape != first.shape or k.dtype != first.dtype for k in values):
            return None
        states = np.frombuffer(b''.join(k.data for k in values), dtype=np.dtype(first.dtype))
        return {
            'states': states.reshape((len(values),) + tuple(first.shape)),
            'values': np.fromiter(values.values(), dtype=float, count=len(values))
        }

    @staticmethod
    def _decode(record):
        from .Analytic import _HashableNumPyIndex
        return {_HashableNumPyIndex(s): float(v) for s, v in zip(record['states'], record['values'])}

    # ---- mapping ---- #
    def _load(self, poly_key):
        key = self._store_key(poly_key)
        if key is None:
            return None
        record = self.store.get(key, default=None)
        if record is None:
            return None
        self.disk_hits += 1
        values = self._decode(record)
        self._persisted[poly_key] = len(values)
        return values

    def __contains__(self, poly_key):
        if poly_key in self.memory:
            return True
        values = self._load(poly_key)
        if values is None:
            return False
        self.memory[poly_key] = values
        self._touched[poly_key] = None
        return True

    def __getitem__(self, poly_key):
        if poly_key not in self:
            raise KeyError(poly_key)
        self._touched[poly_key] = None  # callers add states to the returned dict in place
        return self.memory[poly_key]

    def __setitem__(self, poly_key, value):
        self.memory[poly_key] = value
        self._touched[poly_key] = None

    def __delitem__(self, poly_key):
        del self.memory[poly_key]
        self._touched.pop(poly_key, None)

    def __iter__(self):
        return iter(self._touched)

    def __len__(self):
        return len(self._touched)

    def flush(self):
        """Writes every polynomial whose value dict grew since it was loaded or last written; returns how many."""
        written = 0
        for poly_key in list(self._touched):
            values = self.memory.get(poly_key, None)
            if not isinstance(values, dict) or len(values) <= self._persisted.get(poly_key, 0):
                continue
            key = self._store_key(poly_key)
            record = None if key is None else self._encode(values)
            if record is None:
                continue
            self.store.put(record, key=key, replace=True, meta={'kind': _NAMESPACE, 'entries': len(values)})
            self._persisted[poly_key] = len(values)
            written += 1
        self.disk_writes += written
        return written


# ---------------------------------------------------------------------------------------------
#  CompiledExpressions (formerly VPT2/CompiledExpressions.py)
#  Persistent, bounded analytic compilation and native parallel contractions.
#
#  Only final expressions are compiled. Coefficients, frequencies, conditions,
#  state requests and numerical cutoffs remain inputs to the native evaluator.
#  The optional occupation projection is applied after every symbolic product.
# ---------------------------------------------------------------------------------------------

_current_evaluation_plan = contextvars.ContextVar('analytic_evaluation_plan', default=None)
_worker_evaluation_plans = {}

_current_compilation = contextvars.ContextVar('analytic_compilation', default=None)
_current_free_index_reduction = contextvars.ContextVar('analytic_free_index_reduction', default=False)
_current_coefficient_fingerprints = contextvars.ContextVar('analytic_coefficient_fingerprints', default=None)


@contextlib.contextmanager
def compiled_evaluation_context(store, compilation_options=None, evaluation_options=None,
                                parallelizer=None):
    """Bind immutable coefficient tensors for one public evaluation.

    Array content is fingerprinted on first binding and reused within this
    context; changing a bound tensor requires a new context. Scalar coefficient
    values, occupations and projector bindings remain request-specific. Nested
    public calls reuse the existing context through _compiled_evaluation_context.
    Every actually used live parallelizer receives worker-cache cleanup, including
    per-request overrides. Pool ownership remains with the native contexts.
    """
    current = _current_evaluation_plan.get()
    active = _current_compilation.get()
    owner = current is None or (active is not None and active[0] is not store)
    plan = (uuid.uuid4().hex, StreamPreparationContext(
        None, immutable_expressions=True, options=evaluation_options)) if owner else current
    cache = _current_coefficient_fingerprints.get()
    if owner:
        cache = {}
    coefficient_token = _current_coefficient_fingerprints.set(cache)
    plan_token = _current_evaluation_plan.set(plan)
    token = _current_compilation.set(None if store is None else
        (store, compilation_options or {}, evaluation_options or {}, parallelizer))
    manager = parallelizer if owner and parallelizer is not None else contextlib.nullcontext()
    try:
        with manager:
            try:
                yield plan[1]
            finally:
                if owner:
                    try:
                        with contextlib.ExitStack() as workers:
                            for backend in tuple(plan[1].parallelizers.values()):
                                if backend.active:
                                    workers.callback(backend.apply,_release_evaluation_plan,
                                                     plan[0],cleanup=False)
                    finally:
                        try:
                            _release_evaluation_plan(plan[0])
                        finally:
                            plan[1].close()
    finally:
        _current_compilation.reset(token)
        _current_evaluation_plan.reset(plan_token)
        _current_coefficient_fingerprints.reset(coefficient_token)
        if owner:
            cache.clear()


def _release_evaluation_plan(identity, parallelizer=None):
    from .Analytic import PerturbationTheoryExpressionEvaluator
    plan = _worker_evaluation_plans.pop(identity, None)
    if plan is not None:
        plan.close()
    IndexedBlockEvaluator.clear_plan_cache()
    PerturbationTheoryExpressionEvaluator._cached_expansion = None
    PerturbationTheoryExpressionEvaluator._cached_stream_preparation = None


def _evaluation_plan(context):
    identity = context.get('plan_identity')
    if identity is None:
        return StreamPreparationContext(None, immutable_expressions=True,
                                        options=context['numerical_options'])
    current = _current_evaluation_plan.get()
    if current is not None and current[0] == identity and current[1].pid == os.getpid():
        return current[1]
    if identity not in _worker_evaluation_plans:
        _worker_evaluation_plans[identity] = StreamPreparationContext(
            None, immutable_expressions=True, options=context['numerical_options'])
    return _worker_evaluation_plans[identity]



def _coefficient_binding(value):
    # Only physical tensors are immutable within a public evaluation.
    cache = _current_coefficient_fingerprints.get()
    if cache is None or not isinstance(value, np.ndarray) or value.ndim == 0 or value.dtype.hasobject:
        return None
    key = id(value)
    entry = cache.get(key)
    if entry is None or entry[0]() is not value:
        def dropped(ref, key=key):
            old = cache.get(key)
            if old is not None and old[0] is ref:
                cache.pop(key, None)
        entry = (weakref.ref(value, dropped), _freeze(value), {})
        cache[key] = entry
    return entry


def _active_coefficient_operators(coefficients):
    active = set()
    for expansion in coefficients:
        for order, row in enumerate(expansion):
            for kind, coefficient in enumerate(row):
                entry = _coefficient_binding(coefficient)
                if entry is None:
                    nonzero = bool(np.any(coefficient))
                else:
                    metadata = entry[2]
                    if 'nonzero' not in metadata:
                        metadata['nonzero'] = bool(np.any(coefficient))
                    nonzero = metadata['nonzero']
                if nonzero:
                    active.add((order, kind))
    return tuple(sorted(active))

def current_compilation():
    return _current_compilation.get()


def _atomic_json(path, value):
    atomic_write_json(path, value, sort_keys=True)


def _atomic_pickle(path, value):
    # Preserve shared native PolyPath nodes. JSON/NPZ object conversion expands
    # those graphs; ordinary pickle is the storage format, not a new algebra.
    atomic_write(path, lambda stream: pickle.dump(value, stream, protocol=5), fsync=False, keep_mode=False)


_file_hash = file_hash

# Evaluation caches and handles that are not part of an expression's definition.
_RUNTIME_ATTRS = ('logger', 'checkpoint', '_cascade_manager', '_changes', '_exprs', '_combined', '_eager',
                  '_op_maps', '_conditional_stacks', '_cache', '_grouped_terms')


def _freeze_solver(value, kind, rec):
    names = ['hamiltonian_expansion', 'allowed_terms', 'allowed_coefficients',
             'disallowed_coefficients', 'allowed_energy_changes',
             'intermediate_normalization', 'polynomial_representation', 'resonance_condition']
    return [kind, {name: rec(getattr(value, name)) for name in names}]


def _reject_parallelizer(value, kind, rec):
    raise TypeError('a parallelizer is runtime state, not an expression definition')


def _freeze_options():
    from McUtils.Parallelizers import SharedMemoryArrayTree
    from .Analytic import AnalyticPerturbationTheorySolver
    return dict(
        pre_handlers=[(lambda v: isinstance(v, SharedMemoryArrayTree), lambda v, rec: rec(v.tree))],
        handlers=[(lambda v: isinstance(v, AnalyticPerturbationTheorySolver), _freeze_solver),
                  (lambda v: isinstance(v, Parallelizer), _reject_parallelizer)],
        exclude_attrs=_RUNTIME_ATTRS
    )


def _freeze(value, seen=None):
    """Fingerprint definitions without including mutable evaluation caches."""
    return freeze(value, seen=seen, **_freeze_options())


def _fingerprint(value):
    return fingerprint(value, **_freeze_options())


def _implementation_signature():
    return implementation_signature('McUtils', 'Psience')


def _polynomial_degree(expression):
    """Bound operator letters from the actual DAG, independent of VPT order."""
    from .Analytic import PTTensorCoeffProductDAG, PTTensorCoeffProductSum
    def leaf(value):
        if isinstance(value, PTTensorCoeffProductSum):
            return max((sum(len(c)-2 for c in key) for key in value.terms), default=0)
        return 0
    if not isinstance(expression, PTTensorCoeffProductDAG):
        return leaf(expression)
    cache = {};pending = [(expression, False)]
    while pending:
        node, ready = pending.pop()
        if node in cache:continue
        kind, *args = node._node
        children = [arg for arg in args if isinstance(arg, PTTensorCoeffProductDAG)]
        if not ready:
            pending.append((node, True));pending.extend((child, False) for child in children)
        elif kind == 'leaf':cache[node] = leaf(args[0])
        elif kind == 'add':cache[node] = max(cache[child] for child in children)
        elif kind in {'mul_along','mul_simple','rmul_along','rmul_simple'}:
            cache[node] = sum(cache[arg] if isinstance(arg, PTTensorCoeffProductDAG) else leaf(arg)
                              for arg in args[:2])
        else:cache[node] = max((cache[child] for child in children), default=0)
    return cache[expression]


def _canonical_free_index_key(key, num_fixed, *, max_permutations=720,
                              symmetrizers=None, cache=None):
    """Return (representative, new_to_old_axes), or None for bounded fallback.

    External axes [0,num_fixed) remain pointwise fixed. Equality classes and
    coefficient identities are preserved. Tensor-slot symmetries come only
    from PTTensorCoeffProductSum.canonical_key; G never acquires V symmetry.
    The representative must be evaluated on new_to_old_axes, and numerator,
    denominator and projector deltas must be permuted with the same binding.
    In particular, do not replace a source term by this key on an unchanged
    sorted mode subset.
    """
    from Psience.VPT2.Analytic import PTTensorCoeffProductSum
    key = tuple(tuple(int(i) for i in coefficient) for coefficient in key)
    if num_fixed < 0:
        raise ValueError("num_fixed must be nonnegative")
    indices = sorted({i for coefficient in key for i in coefficient[2:]})
    rank = len(indices)
    if indices != list(range(rank)) or num_fixed > rank:
        return None
    if math.factorial(rank-num_fixed) > max_permutations:
        return None
    cache_key=(num_fixed,key,tuple(PTTensorCoeffProductSum.symmetrizers()
               if symmetrizers is None else symmetrizers))
    if cache is not None:
        cached=cache.get(cache_key,())
        if cached:return cached
    fixed = tuple(range(num_fixed))
    representative = None
    binding = None
    for tail in itertools.permutations(range(num_fixed, rank)):
        old_to_new = fixed + tail
        candidate = PTTensorCoeffProductSum.canonical_key(
            tuple(coefficient[:2] + tuple(old_to_new[i] for i in coefficient[2:])
                  for coefficient in key), symmetrizers=symmetrizers)
        inverse = [0] * rank
        for old, new in enumerate(old_to_new):
            inverse[new] = old
        inverse = tuple(inverse)
        if representative is None or (candidate, inverse) < (representative, binding):
            representative, binding = candidate, inverse
    result=(representative,binding)
    if cache is not None:cache[cache_key]=result
    return result


def _free_index_value_key(value):
    """Exact structural key; includes every denominator/projector factor."""
    from .Analytic import PTEnergyChangeProductSum, PolyPath
    if isinstance(value,PTEnergyChangeProductSum):
        return (PTEnergyChangeProductSum,value.prefactor,frozenset(
            (energy,PolyPath.from_polynomial(polynomial))
            for energy,polynomial in value.terms.items()))
    return (PolyPath,PolyPath.from_polynomial(value))


def _normalize_index_paths(value,cache):
    """Rebuild existing lazy paths through exact native remap simplification.

    The walk is iterative and memoized per immutable part. Original paths stay
    untouched; no conversion to eager polynomial sums is used.
    """
    from .Analytic import PolyPath, PTEnergyChangeProductSum
    def path(root):
        root=PolyPath.from_polynomial(root)
        pending=[(root,False)]
        def children(arg):
            if isinstance(arg,PolyPath):yield arg
            elif isinstance(arg,tuple):
                for item in arg:yield from children(item)
        def remap(arg):
            if isinstance(arg,PolyPath):return cache[id(arg)][1]
            if isinstance(arg,tuple):return tuple(remap(item) for item in arg)
            return arg
        while pending:
            node,ready=pending.pop()
            entry=cache.get(id(node))
            if entry is not None and entry[0] is node:continue
            if node._node is None:
                cache[id(node)]=(node,node);continue
            if not ready:
                pending.append((node,True))
                pending.extend((child,False) for child in children(node._node[1]))
            else:
                kind,args=node._node
                reduced=type(node)._from_node(kind,remap(args),node.ndim,node.order)
                cache[id(node)]=(node,reduced)
        return cache[id(root)][1]
    if isinstance(value,PTEnergyChangeProductSum):
        return value.mutate({energy:path(poly) for energy,poly in value.terms.items()})
    return path(value)


def _prepare_free_index_plan(expression,num_fixed,cache=None,*,rank=None):
    """Lower equivalent numerator/denominator plans with explicit bindings.

    Original coefficient keys and values stay intact for native gathering,
    numerical cutoffs, external-degeneracy fallback and reference evaluation.
    These annotations share preparation, not ordering sectors or scalar sums.
    """
    # Compiler parts already have a known homogeneous rank. Return before
    # scanning values, canonicalizing coefficient keys, walking lazy graphs,
    # or hashing numerators when no dummy-index permutation exists.
    if rank is not None and rank-num_fixed<2:
        expression.__dict__.pop('_free_index_plan_data',None)
        return None
    bindings={};values={};interned={};representatives=set();fallbacks=0;paths={}
    passthrough=0
    for key,value in expression.terms.items():
        key_rank=len({i for coefficient in key for i in coefficient[2:]})
        if key_rank-num_fixed<2:
            bindings[key]=None;values[key]=value;passthrough+=1
            continue
        orbit=_canonical_free_index_key(key,num_fixed,cache=cache)
        if orbit is None:
            bindings[key]=None;values[key]=value;fallbacks+=1;passthrough+=1
            continue
        else:
            representative,binding=orbit;representatives.add(representative)
            if binding==tuple(range(len(binding))):
                canonical=value
            else:
                canonical=value.ensure_dimension(len(binding)).permute(binding,allow_padding=True)
        canonical=_normalize_index_paths(canonical,paths)
        identity=_free_index_value_key(canonical)
        canonical=interned.setdefault(identity,canonical)
        bindings[key]=binding;values[key]=canonical
    if passthrough==len(expression.terms):
        expression.__dict__.pop('_free_index_plan_data',None)
        return None
    plan_count=len({id(value) for value in values.values()})
    expression._free_index_plan_data=dict(format=1,num_fixed=num_fixed,
        bindings=bindings,values=values,stats=dict(
            source_keys=len(expression.terms),coefficient_patterns=len(representatives),
            canonical_plans=plan_count,shared_plans=len(values)-plan_count,
            fallback_keys=fallbacks,native_passthrough_keys=passthrough))
    return expression._free_index_plan_data



@dataclasses.dataclass(frozen=True)
class _OccupationProjection:
    before: int
    after: int
    background: int = 0

    def compiler(self):
        from .Analytic import PolyPath, PTEnergyChangeProductSum
        from .Analytic import canonical_condition_key
        projection = self
        class Transform:
            def __init__(self):
                self.cache = BoundedCache.sampled(100000, 128*1024**2)
                self.nodes = self.cache.view('projection_nodes');self.axes = self.cache.view('projection_axes')
                self.memo = {};self.constants = {}
            def difference(self, value, rank):
                if np.isscalar(value):return np.longdouble(0)
                key = (value, rank)
                if key not in self.memo:
                    path = PolyPath.from_polynomial(value)
                    states = np.full((2, max(rank, path.ndim)), projection.background, dtype=int)
                    states[:,0] = [projection.before, projection.after]
                    values = np.asarray(path.evaluate_polynomial(states,node_cache=self.nodes,axis_cache=self.axes),dtype=np.longdouble)
                    if len(self.memo)>=10000:self.memo.pop(next(iter(self.memo)))
                    self.memo[key] = values[1]-values[0]
                return self.memo[key]
            def constant(self,value,rank):
                key = (float(value),rank)
                if key not in self.constants:
                    self.constants[key] = PolyPath.from_coeffs([np.ones(1)]*rank,prefactor=key[0],steps=0)
                return self.constants[key]
            def __call__(self,key,value):
                rank = len({i for c in key for i in c[2:]})
                if rank == 0:return None
                if isinstance(value,PTEnergyChangeProductSum):
                    rows = {}
                    for energy,numerator in value.terms.items():
                        if any(getattr(row[0],'_state_projection_factor',False) or row[0] not in {0,1} for row in energy):
                            raise ValueError('occupation projection requires occupation-independent harmonic denominators')
                        canonical,scale = canonical_condition_key(energy)
                        if canonical is None:raise ValueError('invalid harmonic denominator')
                        rows[canonical] = rows.get(canonical,np.longdouble(0))+scale*self.difference(numerator,rank)
                    rows = {energy:self.constant(weight,rank) for energy,weight in rows.items() if weight!=0}
                    return value.mutate(rows,reduced=True) if rows else None
                weight = self.difference(value,rank)
                return self.constant(weight,rank) if weight else None
        return Transform()


class _CompiledEntries(EntryStore):
    # compiled expressions live at folder/<signature>/manifest.json, aliases at folder/aliases/<signature>.json
    meta_name = 'manifest.json'
    shard = 0
    alias_field = 'compiled'
    stage_prefix = '.compile-'


class AnalyticExpressionStore:
    format_version = 1
    def __init__(self, folder=None, *, memory_cache_size=128, memory_cache_bytes=512*1024**2):
        self.folder = None if folder is None else Path(folder).expanduser().resolve()
        self.entries = None
        if self.folder is not None:
            self.folder.mkdir(parents=True,exist_ok=True)
            self.entries = _CompiledEntries(self.folder)
        self.memory_cache_size = max(0, int(memory_cache_size))
        self.memory_cache_bytes = None if memory_cache_bytes is None else max(0, int(memory_cache_bytes))
        self._memory = collections.OrderedDict()
        self._aliases = collections.OrderedDict()
        self.memory_bytes = 0
        self.timings = dict(requests=0, compilations=0, cache_hits=0, failures=0,
            compilation_seconds=0., cache_lookup_seconds=0., preparation_seconds=0., evaluation_seconds=0.)

    def statistics(self):
        return dict(self.timings, storage='memory' if self.folder is None else 'disk',
                    memory_entries=len(self._memory), memory_payload_bytes=self.memory_bytes)

    def clear(self):
        """Release reusable in-memory entries; never remove persistent cache files."""
        self._memory.clear(); self._aliases.clear(); self.memory_bytes=0

    def _remember(self, compiled):
        size = sum(part['bytes'] for part in compiled.manifest['parts'])
        key = compiled.manifest['signature']
        if self.memory_cache_size and (self.memory_cache_bytes is None or size <= self.memory_cache_bytes):
            previous = self._memory.pop(key, None)
            if previous is not None:self.memory_bytes -= sum(p['bytes'] for p in previous.manifest['parts'])
            self._memory[key]=compiled; self.memory_bytes+=size
            while len(self._memory)>self.memory_cache_size or (self.memory_cache_bytes is not None and self.memory_bytes>self.memory_cache_bytes):
                removed, value=self._memory.popitem(last=False)
                self.memory_bytes-=sum(p['bytes'] for p in value.manifest['parts'])
                for alias,target in list(self._aliases.items()):
                    if target==removed:del self._aliases[alias]
        return compiled

    def _get_memory(self, signature):
        value=self._memory.get(signature)
        if value is not None:self._memory.move_to_end(signature)
        return value

    def _signature(self,context,active_operators,projection,options):
        return _fingerprint(dict(format=self.format_version,implementation=_implementation_signature(),
            definition=context,active_operators=active_operators,projection=projection,options=options))

    def _find_alias(self,signature):
        if self.folder is None:
            target=self._aliases.get(signature)
            return None if target is None else self._get_memory(target)
        target=self.entries.resolve_alias(signature)
        return None if target is None else CompiledAnalyticExpression.load(self.folder/target)

    def compile(self, expression, *, context=None, active_operators=None, projection=None,
                max_terms=32768, max_bytes=64*1024**2, cache_size=100000,
                cache_bytes=512*1024**2, node_cache_size=8192, mode_dimension=None, progress=None,
                coalesce_combine_nodes=True,bounded_byte_accounting=True,free_index_reduction=True):
        """Compile native expressions; early DAG combine-node coalescing is on by default.

        Set coalesce_combine_nodes=False (also accepted through evaluator
        compilation_options) to disable that optimization for diagnostics or
        benchmarking. Bounded final-part accumulation remains active.
        """
        return self.compile_factory(lambda:expression, context=expression if context is None else context,
            active_operators=active_operators,projection=projection,max_terms=max_terms,max_bytes=max_bytes,
            cache_size=cache_size,cache_bytes=cache_bytes,node_cache_size=node_cache_size,
            mode_dimension=mode_dimension,progress=progress,coalesce_combine_nodes=coalesce_combine_nodes,
            bounded_byte_accounting=bounded_byte_accounting,free_index_reduction=free_index_reduction)

    def compile_factory(self, factory, **options):
        """Compile or reuse an entry, accounting separately for lookup and actual construction."""
        started=time.perf_counter(); building=[]; completed=False
        from .Analytic import PolyPath
        reduction=options.get('free_index_reduction',True)
        if not isinstance(reduction,(bool,np.bool_)):
            raise TypeError('free_index_reduction must be boolean')
        index_token=PolyPath._index_reduction_context.set(bool(reduction))
        self.timings['requests']+=1
        def build():
            building.append(time.perf_counter())
            return factory()
        try:
            result=self._compile_factory(build, **options)
            result.store=self; completed=True
            return result
        finally:
            PolyPath._index_reduction_context.reset(index_token)
            elapsed=time.perf_counter()-started
            construction=0. if not building else time.perf_counter()-building[0]
            self.timings['compilation_seconds']+=construction
            self.timings['cache_lookup_seconds']+=max(0.,elapsed-construction)
            self.timings['compilations' if building else 'cache_hits']+=int(completed)
            self.timings['failures']+=int(not completed)

    def _compile_factory(self, factory, *, context, active_operators=None, projection=None,
                         max_terms=32768,max_bytes=64*1024**2,cache_size=100000,
                          cache_bytes=512*1024**2,node_cache_size=8192,mode_dimension=None,progress=None,
                          coalesce_combine_nodes=True,bounded_byte_accounting=True,free_index_reduction=True):
        from .Analytic import (SqrtChangePoly,PTTensorCoeffProductDAG,PTTensorCoeffProductSum,
            PTTensorCoeffProductDAGEvaluationPlan,_BalancedPolynomialAccumulator)
        if not isinstance(coalesce_combine_nodes, (bool, np.bool_)):
            raise TypeError('coalesce_combine_nodes must be boolean')
        coalescing_enabled = bool(coalesce_combine_nodes)
        if not isinstance(free_index_reduction,(bool,np.bool_)):
            raise TypeError('free_index_reduction must be boolean')
        free_index_reduction=bool(free_index_reduction)
        options = dict(max_terms=max_terms,max_bytes=max_bytes,cache_size=cache_size,
                       cache_bytes=cache_bytes,node_cache_size=node_cache_size,mode_dimension=mode_dimension,
                       coalesce_combine_nodes=coalescing_enabled,bounded_byte_accounting=bool(bounded_byte_accounting),
                       free_index_reduction=free_index_reduction)
        if max_terms<1 or (max_bytes is not None and max_bytes<1):raise ValueError('invalid compilation budget')
        signature = self._signature(context,active_operators,projection,options)
        directory = None if self.folder is None else self.folder/signature
        if self.folder is None:
            cached=self._get_memory(signature)
            if cached is not None:return cached
        elif self.entries.has(signature):return CompiledAnalyticExpression.load(directory)
        alias=self._find_alias(signature)
        if alias is None and active_operators is not None:
            # A cache compiled without pruning is valid for a narrower active
            # coefficient set. The reverse substitution would lose terms.
            alias=self._find_alias(self._signature(context,None,projection,options))
        if alias is not None:return alias
        started = time.perf_counter()
        wrapper = factory()
        expr = wrapper.expr
        if np.isscalar(expr):
            if expr!=0:raise ValueError('nonzero scalar expression needs a native identity polynomial')
            change = tuple(wrapper.change) if wrapper.change is not None else ()
            baseline = None;poly = PTTensorCoeffProductSum({})
        elif isinstance(expr,SqrtChangePoly):
            change = tuple(int(i) for i in expr.poly_change)
            baseline = tuple(int(i) for i in expr.shift_start);poly = expr.poly_obj
        else:
            change = tuple(wrapper.change) if wrapper.change is not None else ()
            baseline = None;poly = expr
        if not isinstance(poly,PTTensorCoeffProductSum):poly = PTTensorCoeffProductSum({():poly})
        if projection is not None and (change or baseline and any(baseline)):
            raise ValueError('final occupation projection requires a diagonal expression')
        if isinstance(poly,PTTensorCoeffProductDAG) and active_operators is not None:
            poly = poly.prune_operators(poly.operator_keys-set(active_operators))
        degree = _polynomial_degree(poly)
        maximum = len(change)+max(0,(degree-sum(abs(i) for i in change))//2)
        if mode_dimension is not None:
            if mode_dimension<1:raise ValueError('mode_dimension must be positive')
            maximum=min(maximum,mode_dimension)
        if isinstance(poly,PTTensorCoeffProductDAG):
            class CompilationPlan(PTTensorCoeffProductDAGEvaluationPlan):
                coalesce_combine_nodes = coalescing_enabled
            plan = CompilationPlan(poly,
                cache=BoundedCache.sampled(cache_size,cache_bytes),max_cached_terms=node_cache_size,
                max_mode_indices=maximum)
            iterator = plan.iter_terms()
            if progress is not None:
                plan.progress_callback=lambda p:progress(dict(stage='generating',**p.stats()))
        else:
            plan = None;iterator = poly.terms.items()
        transform = projection.compiler() if projection is not None else None
        budget = StreamBudget(dict(max_terms=max_terms,max_bytes=max_bytes,
            bounded_byte_accounting=bounded_byte_accounting))
        staging = contextlib.nullcontext(None) if self.folder is None else self.entries.stage(signature)
        memory_parts={}
        if free_index_reduction and not hasattr(self,'_free_index_orbits'):
            self._free_index_orbits=BoundedCache.sampled(32768,32*1024**2)
        can_reduce=free_index_reduction and not (
            baseline and any(baseline[wrapper.num_fixed:]))
        with staging as staged:
            stage = None if staged is None else Path(staged.path)
            pending={};parts=[]
            def flush(reason):
                if not pending:return
                compact = PTTensorCoeffProductSum({key:acc.total() for key,acc in pending.items()},
                    canonicalize=False,reduced=False).combine()
                groups = {}
                for key,value in compact.terms.items():
                    rank = len({i for coefficient in key for i in coefficient[2:]})
                    groups.setdefault(rank,{})[key] = value
                for rank,terms in sorted(groups.items()):
                    name = 'part-{:06d}-r{}.pkl'.format(len(parts),rank)
                    part = PTTensorCoeffProductSum(terms,canonicalize=False,reduced=True)
                    if can_reduce:
                        _prepare_free_index_plan(part,wrapper.num_fixed,self._free_index_orbits,rank=rank)
                    if stage is None:
                        payload=pickle.dumps(part,protocol=5)
                        part_hash=hashlib.sha256(payload).hexdigest();part_bytes=len(payload)
                        memory_parts[name]=part
                    else:
                        _atomic_pickle(stage/name,part)
                        part_hash=_file_hash(stage/name);part_bytes=(stage/name).stat().st_size
                    parts.append(dict(file=name,rank=rank,sha256=part_hash,bytes=part_bytes,
                        coefficient_keys=len(terms),energy_rows=sum(len(p.terms) if hasattr(p,'terms') else 1 for p in terms.values()),
                        free_index_plan_stats=(getattr(part,'_free_index_plan_data',{}).get('stats'))))
                pending.clear();budget.flushed(reason)
            for key,value in iterator:
                if transform is not None:value = transform(key,value)
                if value is None:continue
                pending.setdefault(key,_BalancedPolynomialAccumulator()).append(value)
                reason = budget.incoming(pending)
                if reason:flush(reason)
                if progress is not None:progress(dict(stage='compiling',parts=len(parts),**budget.stats()))
            flush('end')
            manifest = dict(format=self.format_version,complete=True,signature=signature,parts=parts,
                change=change,baseline_shift=baseline,num_fixed=wrapper.num_fixed,
                projection=None if projection is None else dataclasses.asdict(projection),
                active_operators=None if active_operators is None else sorted(active_operators),
                operator_letter_bound=degree,maximum_mode_indices=maximum,options=options,
                seconds=time.perf_counter()-started,stream_stats=budget.stats(),
                plan_stats=None if plan is None else plan.stats(),implementation=_implementation_signature())
            if stage is not None:_atomic_json(stage/'manifest.json',manifest)
        compiled=(self._remember(CompiledAnalyticExpression(None,manifest,parts=memory_parts))
                  if directory is None else CompiledAnalyticExpression.load(directory))
        if projection is None:
            native=[wrapper.op if wrapper.op is not None else poly]
            from .Analytic import CorrectionView
            if isinstance(wrapper.op,CorrectionView) and not any(getattr(wrapper.op,name) is not None for name in
                    ['allowed_terms','allowed_coefficients','disallowed_coefficients','allowed_energy_changes']):
                native.append(wrapper.op.term)
            for op in native:
                key=self._signature(dict(operator=op,change=change,baseline=baseline),active_operators,None,options)
                if self.folder is None:
                    if signature in self._memory:
                        self._aliases[key]=signature
                        while len(self._aliases)>4*self.memory_cache_size:self._aliases.popitem(last=False)
                else:self.entries.set_alias(key,signature)
        return compiled


@dataclasses.dataclass(frozen=True)
class _TemplateSubset:
    rank: int
    def count(self,ndim,num_fixed,free_count):return int(free_count==self.rank-num_fixed)
    def iter_subsets(self,ndim,num_fixed,free_count):
        if free_count==self.rank-num_fixed:yield tuple(range(num_fixed,self.rank))


class CompiledAnalyticExpression:
    def __init__(self,folder,manifest,parts=None):
        self.folder=None if folder is None else Path(folder)
        self.manifest=manifest;self.parts=parts;self.store=None
    @classmethod
    def load(cls,folder):
        folder=Path(folder).resolve();manifest=json.loads((folder/'manifest.json').read_text())
        if not manifest.get('complete') or manifest.get('format')!=AnalyticExpressionStore.format_version:
            raise ValueError('incomplete or unsupported compiled expression')
        if manifest['implementation']!=_implementation_signature():raise ValueError('compiled expression source version changed')
        for part in manifest['parts']:
            if Path(part['file']).name!=part['file']:raise ValueError('invalid expression part path')
            if _file_hash(folder/part['file'])!=part['sha256']:raise ValueError('compiled expression checksum mismatch')
        return cls(folder,manifest)
    def prepare_evaluation(self,state_perms,coeffs,freqs,*,anchored_modes=(),parallelizer=None,
                           batch_size=512,task_size=32768,resume=True,use_shared_memory=None,**numerical_options):
        started=time.perf_counter()
        try:
            return CompiledAnalyticEvaluation(self,state_perms,coeffs,freqs,anchored_modes=anchored_modes,
                parallelizer=parallelizer,batch_size=batch_size,task_size=task_size,resume=resume,
                numerical_options=numerical_options,use_shared_memory=use_shared_memory)
        finally:
            if self.store is not None:self.store.timings['preparation_seconds']+=time.perf_counter()-started
    def evaluate(self,state_perms,coeffs,freqs,**options):
        return self.prepare_evaluation(state_perms,coeffs,freqs,**options).run()


class CompiledAnalyticEvaluation:
    """A reusable numerical task plan; the parallelizer never sees a symbolic DAG."""
    def __init__(self,compiled,state_perms,coeffs,freqs,*,anchored_modes=(),parallelizer=None,
                 batch_size=512,task_size=32768,resume=True,numerical_options=None,use_shared_memory=None):
        self.compiled=compiled
        self.state_perms=[(np.asarray(state,dtype=int),np.asarray(perms,dtype=int)) for state,perms in state_perms]
        from McUtils.Parallelizers import SharedMemoryArrayTree
        from .Analytic import PerturbationTheoryExpressionEvaluator
        active = _current_evaluation_plan.get()
        self.coeffs = (coeffs if isinstance(coeffs,SharedMemoryArrayTree) else
            PerturbationTheoryExpressionEvaluator._normalize_coefficient_expansions(
                coeffs,preparation_context=None if active is None else active[1]))
        self.freqs=np.asarray(freqs);self.anchored_modes=tuple(anchored_modes)
        cap=compiled.manifest['options'].get('mode_dimension')
        if cap is not None and len(self.freqs)>cap:
            raise ValueError('compiled expression was pruned for fewer modes; recompile for the larger basis')
        self.parallelizer=Parallelizer.lookup(parallelizer) if parallelizer is not None else None
        self.batch_size=int(batch_size);self.task_size=int(task_size);self.resume=bool(resume)
        if min(self.batch_size,self.task_size)<1:raise ValueError('batch/task sizes must be positive')
        self.numerical_options=dict(numerical_options or {})
        forbidden={'parallelizer','mode_subset_selector','evaluation_mode'} & self.numerical_options.keys()
        if forbidden:raise TypeError('compiled numerical options override scheduling: {}'.format(sorted(forbidden)))
        self.use_shared_memory=use_shared_memory is not False
        projection=compiled.manifest['projection']
        if projection is not None:
            if len(self.anchored_modes)!=1:raise ValueError('projected transition requires its one anchored physical mode')
            anchor=self.anchored_modes[0]
            for state,perms in self.state_perms:
                expected=np.full(len(self.freqs),projection['background'],dtype=int);expected[anchor]=projection['after']
                if not np.array_equal(state,expected) or not np.array_equal(perms,np.arange(len(self.freqs))[None]):
                    raise ValueError('occupation-specialized expression used outside its target-state domain')
        elif self.anchored_modes:raise ValueError('anchoring requires a matching final occupation projection')
        signature=self._fingerprint_context(dict(expression=compiled.manifest['signature'],states=self.state_perms,
            coefficients=self.coeffs,frequencies=self.freqs,anchors=self.anchored_modes,numerical=self.numerical_options,
            batch_size=self.batch_size,task_size=self.task_size))
        self.folder=None if compiled.folder is None else compiled.folder/'evaluations'/signature
        if self.folder is not None:self.folder.mkdir(parents=True,exist_ok=True)
        tasks=[];ndim=len(self.freqs);num_fixed=compiled.manifest['num_fixed']
        for part in compiled.manifest['parts']:
            rank=part['rank'];free=rank-num_fixed
            if free<0:raise ValueError('coefficient aliases do not cover the fixed change')
            count=math.comb(ndim-num_fixed-len(self.anchored_modes),free-len(self.anchored_modes)) if free>=len(self.anchored_modes) and free<=ndim-num_fixed else 0
            for start in range(0,count,self.task_size):
                stop=min(start+self.task_size,count);index=len(tasks)
                tasks.append(dict(index=index,part=part,start=start,stop=stop,
                                  file='task-{:08d}.pkl'.format(index)))
        self.tasks=tasks;self.signature=signature
        if self.folder is not None:_atomic_json(self.folder/'plan.json',dict(signature=signature,expression_signature=compiled.manifest['signature'],
            kind='numerical-task-plan',task_count=len(tasks),subset_part_count=sum(t['stop']-t['start'] for t in tasks),
            batch_size=self.batch_size,task_size=self.task_size,tasks=tasks))


    @staticmethod
    def _fingerprint_context(context):
        # Physical coefficient tensors are bound as immutable for the active
        # evaluation context. Cache their frozen content, not runtime inputs.
        cache = _current_coefficient_fingerprints.get()
        if cache is None:
            return _fingerprint(context)
        from McUtils.Parallelizers import SharedMemoryArrayTree
        arrays = {}
        def bind(value):
            if isinstance(value, np.ndarray):
                # Per-element scalar coefficients can change between requests.
                if value.ndim == 0 or value.dtype.hasobject:
                    return
                entry = _coefficient_binding(value)
                arrays[id(value)] = entry[1]
            elif isinstance(value, SharedMemoryArrayTree):
                bind(value.tree)
            elif isinstance(value, (list, tuple)):
                for part in value:
                    bind(part)
            elif isinstance(value, dict):
                for part in value.values():
                    bind(part)
        bind(context['coefficients'])
        options = _freeze_options()
        options['pre_handlers'] = [
            (lambda value: isinstance(value, np.ndarray) and id(value) in arrays,
             lambda value, rec: arrays[id(value)])
        ] + options['pre_handlers']
        return fingerprint(context, **options)

    def _context(self,coeffs=None):
        plan = _current_evaluation_plan.get()
        from McUtils.Parallelizers import SharedMemoryArrayTree
        inputs = self.coeffs.tree if isinstance(self.coeffs,SharedMemoryArrayTree) else self.coeffs
        bindings = None if plan is None else plan[1].bind(inputs)
        return dict(plan_identity=None if plan is None else plan[0],
            coefficient_bindings=bindings,coefficient_owner_pid=os.getpid(),expression_folder=None if self.compiled.folder is None else str(self.compiled.folder),
            expression_parts=self.compiled.parts,manifest=self.compiled.manifest,
            evaluation_folder=None if self.folder is None else str(self.folder),signature=self.signature,state_perms=self.state_perms,
            coeffs=self.coeffs if coeffs is None else coeffs,freqs=self.freqs,anchors=self.anchored_modes,
            batch_size=self.batch_size,resume=self.resume,numerical_options=self.numerical_options)

    def run(self):
        started=time.perf_counter()
        try:
            if _current_evaluation_plan.get() is None:
                manager = compiled_evaluation_context(self.compiled.store,
                    evaluation_options=self.numerical_options,parallelizer=self.parallelizer)
            else:
                manager = contextlib.nullcontext() if self.parallelizer is None else self.parallelizer
            with manager:return self._run()
        finally:
            active = _current_evaluation_plan.get()
            if active is not None and self.parallelizer is not None and not self.parallelizer.active:
                # Per-request pools keep their native lifetime/default policy.
                # Dead pools have no worker state left for the plan to release.
                active[1].parallelizers.pop(id(self.parallelizer),None)
            if self.compiled.store is not None:self.compiled.store.timings['evaluation_seconds']+=time.perf_counter()-started

    def _run(self):
        backend=self.parallelizer
        nprocs=1 if backend is None else backend.nprocs
        bins=[[] for _ in range(nprocs)];loads=[0]*nprocs
        for task in sorted(self.tasks,key=lambda t:(t['stop']-t['start'])*max(1,t['part']['energy_rows']),reverse=True):
            index=min(range(nprocs),key=lambda i:loads[i]);bins[index].append(task)
            loads[index]+=(task['stop']-task['start'])*max(1,task['part']['energy_rows'])
        owner=None
        active = _current_evaluation_plan.get()
        if active is not None and backend is not None:
            active[1].parallelizers[id(backend)] = backend
        try:
            from McUtils.Parallelizers import SharedMemoryArrayTree
            coefficients = self.coeffs
            if isinstance(backend,MultiprocessingParallelizer) and self.use_shared_memory and not isinstance(self.coeffs,SharedMemoryArrayTree):
                normalized = self.coeffs
                if active is not None:
                    coefficients = active[1].share_coefficients(normalized)
                else:
                    owner = SharedMemoryArrayTree.create(normalized)
                    coefficients = owner
            context=self._context(coefficients)
            if backend is None:
                records=_run_compiled_tasks(self.tasks,context)
            else:
                records=backend.apply(_map_compiled_bins,bins,context,cleanup=False)
            records=sorted(records,key=lambda r:r['index'])
            if [r['index'] for r in records]!=list(range(len(self.tasks))):raise ValueError('incomplete or duplicated compiled task coverage')
            outputs=[np.zeros((len(self.coeffs),len(perms)),dtype=np.result_type(
                *[r['values'][position].dtype for r in records],float))
                for position,(_,perms) in enumerate(self.state_perms)]
            for position,out in enumerate(outputs):
                for index in np.ndindex(out.shape):
                    values=[r['values'][position][index] for r in records]
                    if any(np.iscomplexobj(v) for v in values):out[index]=sum(values)
                    else:out[index]=math.fsum(values)
            if self.folder is not None:_atomic_json(self.folder/'result.json',dict(complete=True,signature=self.signature,task_count=len(records),
                subset_part_count=sum(t['stop']-t['start'] for t in self.tasks),
                resumed_tasks=sum(int(r.get('resumed',False)) for r in records)))
            return outputs
        finally:
            if owner is not None:
                from .Analytic import PerturbationTheoryExpressionEvaluator
                IndexedBlockEvaluator.clear_plan_cache()
                PerturbationTheoryExpressionEvaluator._cached_expansion=None
                owner.dispose()


def _map_compiled_bins(bins,context,parallelizer=None):
    return parallelizer.map(_run_assigned_compiled_bins,bins,extra_args=(context,),vectorized=True,aggregate=True)


def _run_assigned_compiled_bins(bins,context):
    return _run_compiled_tasks([task for group in bins for task in group],context)


def _run_compiled_tasks(tasks,context):
    from .Analytic import PolyPath
    limit = context['numerical_options'].get('path_equality_cache_size', 32768)
    enabled=context['numerical_options'].get('free_index_reduction',
        context['manifest']['options'].get('free_index_reduction',False))
    if not isinstance(enabled,(bool,np.bool_)):
        raise TypeError('free_index_reduction must be boolean')
    token=_current_free_index_reduction.set(bool(enabled))
    index_token=PolyPath._index_reduction_context.set(bool(enabled))
    try:
        with PolyPath.cached_comparisons(limit):
            return _run_compiled_tasks_inner(tasks,context)
    finally:
        PolyPath._index_reduction_context.reset(index_token)
        _current_free_index_reduction.reset(token)


def _run_compiled_tasks_inner(tasks,context):
    from .Analytic import PerturbationTheoryExpressionEvaluator as Expression
    from McUtils.Parallelizers import SharedMemoryArrayTree
    evaluation_plan = _evaluation_plan(context)
    def restore(value, binding=None):
        if isinstance(value, SharedMemoryArrayTree):
            if value.owner or context.get('coefficient_owner_pid') == os.getpid():
                # A caller's local handle is borrowed, including attach() handles.
                # Only transport-created worker attachments belong to this plan.
                # Retain a borrowed handle while prepared views refer to it.
                evaluation_plan.borrowed_handles[id(value)] = value
                return value.tree
            # Unpickled descriptors attach on each dispatch. Intern the native
            # attachment so tensor identities and gathers survive that dispatch.
            key = ('shared', value.descriptor.name)
            prior = evaluation_plan.shared_arrays.get(key)
            if prior is None:
                evaluation_plan.shared_arrays[key] = value
                prior = value
            elif prior is not value:
                value.dispose()
            return prior.tree
        if isinstance(value,np.ndarray) and value.ndim and binding is not None:
            # Non-shared dispatch pickles new ndarray objects each time. Intern
            # them by the parent's immutable plan binding, not their new worker
            # addresses, so old transport copies do not accumulate in a plan.
            return evaluation_plan.copied_inputs.setdefault(binding,value)
        if isinstance(value, (list, tuple)):
            return type(value)(restore(part,None if binding is None else binding[index])
                               for index,part in enumerate(value))
        return value
    coefficients = Expression._normalize_coefficient_expansions(
        restore(context['coeffs'],context.get('coefficient_bindings')),
        preparation_context=evaluation_plan)
    preparation = evaluation_plan.domain(coefficients)
    polynomial_cache = evaluation_plan.polynomial_cache()
    manifest=context['manifest'];num_fixed=manifest['num_fixed'];ndim=len(context['freqs'])
    anchors=context['anchors'];available=tuple(i for i in range(num_fixed,ndim) if i not in anchors)
    cache=collections.OrderedDict();records=[]
    for task in tasks:
        target=None if context['evaluation_folder'] is None else Path(context['evaluation_folder'])/task['file']
        if target is not None and context['resume'] and target.exists():
            with target.open('rb') as stream:record=pickle.load(stream)
            if record['signature']!=context['signature'] or record['index']!=task['index']:raise ValueError('compiled result context mismatch')
            if record['start']!=task['start'] or not task['start']<=record['stop']<=task['stop']:
                raise ValueError('compiled result range mismatch')
            if record['complete']:
                if record['stop']!=task['stop']:raise ValueError('completed compiled task has a missing range')
                record=dict(record,resumed=True);records.append(record);continue
        else:record=None
        part=task['part'];key=part['file']
        if key not in cache:
            if context.get('expression_parts') is not None:
                cache[key]=context['expression_parts'][key]
            else:
                path=Path(context['expression_folder'])/key
                if _file_hash(path)!=part['sha256']:raise ValueError('compiled part changed in worker')
                with path.open('rb') as stream:cache[key]=pickle.load(stream)
            # Compiled parts are immutable for the lifetime of this worker.
            # Compute invariant polynomial degree once per load, not per batch.
            cache[key]=(cache[key],Expression._get_max_order(cache[key])
                        if context['numerical_options'].get('cache_polynomial_order',True) else None)
            while len(cache)>2:cache.popitem(last=False)
        poly,max_order=cache[key];rank=part['rank'];free=rank-num_fixed
        batches=[];position=task['start']
        journal=None if target is None else target.with_suffix('.batches')
        written_batches=0
        if record is not None:
            position=record['stop']
            if 'batch_count' in record:
                for ordinal in range(record['batch_count']):
                    with (journal/('batch-{:08d}.pkl'.format(ordinal))).open('rb') as stream:
                        batch_record=pickle.load(stream)
                    if batch_record['signature']!=context['signature'] or batch_record['index']!=task['index']:
                        raise ValueError('compiled batch context mismatch')
                    batches.append(batch_record['batch'])
                written_batches=len(batches)
            else:
                # Read old checkpoints without changing their summation order.
                batches=record['batches']
            cursor=task['start']
            for batch in batches:
                if batch['start']!=cursor or not cursor<batch['stop']<=task['stop']:
                    raise ValueError('invalid partial compiled task coverage')
                cursor=batch['stop']
            if cursor!=position:raise ValueError('partial task stop does not match its batches')
        iterator=itertools.islice(itertools.combinations(available,free-len(anchors)),position,task['stop'])
        def total_values():
            arrays=[np.zeros((len(coefficients),len(perms))) for _,perms in context['state_perms']]
            for route,array in enumerate(arrays):
                for index in np.ndindex(array.shape):
                    array[index]=math.fsum(batch['values'][route][index] for batch in batches)
            return arrays
        options=dict(evaluation_mode='materialized',evaluation_backend='indexed',logger=False,log_scaled=False,
            path_cache_size=10000,path_cache_bytes=64*1024**2,
            state_permutation_cache_size=10000,state_permutation_cache_bytes=16*1024**2)
        options.update(context['numerical_options'])
        for name in ('path_equality_cache_size','cache_polynomial_order','free_index_reduction'):
            options.pop(name,None)
        options['_max_polynomial_order']=max_order
        options['_poly_eval_cache']=polynomial_cache
        options['_stream_preparation']=preparation
        for extra_batch in iter(lambda:list(itertools.islice(iterator,context['batch_size'])),[]):
            subsets=[anchors+tuple(extra) for extra in extra_batch]
            requests=[]
            for state,perms in context['state_perms']:
                expanded=[]
                for subset in subsets:
                    order=tuple(range(num_fixed))+subset+tuple(i for i in range(num_fixed,ndim) if i not in subset)
                    expanded.extend(perms[:,order])
                requests.append([state,np.asarray(expanded,dtype=int)])
            result=Expression.evaluate_polynomial_expression(requests,coefficients,context['freqs'],poly,
                manifest['change'],manifest['baseline_shift'],num_fixed,
                mode_subset_selector=_TemplateSubset(rank),**options)
            batch_values=[]
            for array,(_,perms) in zip(result,context['state_perms']):
                array=np.asarray(array).reshape(len(coefficients),len(subsets),len(perms))
                value=np.zeros((len(coefficients),len(perms)))
                for index in np.ndindex(value.shape):value[index]=math.fsum(array[index[0],:,index[1]])
                if not np.isfinite(value).all():raise FloatingPointError('nonfinite compiled contraction')
                batch_values.append(value)
            position+=len(subsets);batches.append(dict(start=position-len(subsets),stop=position,values=batch_values))
            if target is not None and position < task['stop']:
                journal.mkdir(exist_ok=True)
                # Publish only new batch payloads, then a constant-size cursor.
                # An interrupted unpublished tail is overwritten on resume.
                for ordinal in range(written_batches,len(batches)):
                    _atomic_pickle(journal/('batch-{:08d}.pkl'.format(ordinal)),
                        dict(signature=context['signature'],index=task['index'],batch=batches[ordinal]))
                written_batches=len(batches)
                _atomic_pickle(target,dict(complete=False,signature=context['signature'],index=task['index'],
                    start=task['start'],stop=position,batch_count=written_batches))
        if position!=task['stop']:raise ValueError('compiled subset range incomplete')
        values=total_values()
        record=dict(complete=True,signature=context['signature'],index=task['index'],start=task['start'],stop=position,
                    values=values,batches=batches,resumed=record is not None)
        if target is not None:
            _atomic_pickle(target,record)
            if journal.exists():
                for ordinal in range(written_batches):
                    (journal/('batch-{:08d}.pkl'.format(ordinal))).unlink(missing_ok=True)
                if not any(journal.iterdir()):
                    journal.rmdir()
        records.append(record)
    if context.get('plan_identity') is None:
        evaluation_plan.close()
        IndexedBlockEvaluator.clear_plan_cache()
    return records


class AnalyticCorrectionCompiler:
    """Correction selection and final projection belong to the evaluator layer."""
    def __init__(self,evaluator,folder=None):
        self.evaluator=evaluator
        self.store=evaluator.expression_store if folder is None else AnalyticExpressionStore(folder)

    def _generator(self,correction,order,energy_method):
        from .Analytic import SameStateCondition
        if correction=='energy' and energy_method=='auto':
            energy_method='bloch' if self.evaluator.solver.resonance_condition is None or type(self.evaluator.solver.resonance_condition) is SameStateCondition else 'native'
        if correction=='energy' and energy_method=='bloch':
            if self.evaluator.solver.resonance_condition is not None and type(self.evaluator.solver.resonance_condition) is not SameStateCondition:
                raise ValueError('scalar Bloch energies must not replace a nontrivial resonance condition')
            from .Analytic import AnalyticPerturbationTheorySolver
            source=self.evaluator.solver
            prepared=self.evaluator._prep_expansions(self.evaluator.expansions)
            # Remove identically zero Hamiltonian orders before enumerating
            # operator paths, rather than only pruning the resulting tensor DAG.
            expansion=[term if n==0 or n<len(prepared) and any(np.any(c) for c in prepared[n])
                       else ZeroTerm(logger=False) for n,term in enumerate(source.hamiltonian_expansion)]
            solver=AnalyticPerturbationTheorySolver(expansion,logger=source.logger,
                polynomial_representation=source.polynomial_representation,
                resonance_condition=source.resonance_condition,
                allowed_terms=source.allowed_terms,allowed_coefficients=source.allowed_coefficients,
                disallowed_coefficients=source.disallowed_coefficients,
                allowed_energy_changes=source.allowed_energy_changes,
                intermediate_normalization=source.intermediate_normalization)
            terms=AnalyticBlochSeries(solver).hermitian_hamiltonian(order)
            return lambda n:terms[n]
        names={'energy':'energy_correction','overlap':'overlap_correction',
               'wavefunction':'wavefunction_correction','full_wavefunction':'full_wavefunction_correction',
               'operator':'operator_correction'}
        if callable(correction):return correction
        if correction not in names:raise ValueError('unknown analytic correction '+str(correction))
        return getattr(self.evaluator.solver,names[correction])

    def precompile(self,correction='energy',order=None,orders=None,changes=((),),energy_method='native',projection=None,
                   compilation_options=None):
        from .Analytic import PerturbationTheoryExpressionEvaluator
        evaluator=self.evaluator
        if order is None:order=len(evaluator.expansions)-1
        if orders is None:orders=range(0,order+1,2 if correction=='energy' else 1)
        prepared=evaluator._prep_expansions(evaluator.expansions)
        active=tuple((n,k) for n,row in enumerate(prepared) for k,c in enumerate(row) if np.any(c))
        # An operator correction can introduce coefficients absent from the
        # Hamiltonian payload. Their absence here is not evidence of a zero.
        if correction not in {'energy','overlap','wavefunction','full_wavefunction'}:
            active=None
        options=dict(evaluator.compilation_options);options.update(compilation_options or {})
        options.setdefault('mode_dimension',len(evaluator.freqs))
        source=dict(solver=evaluator.solver,correction=correction,energy_method=energy_method)
        generator=[None]
        def expression(n,change):
            if generator[0] is None:generator[0]=self._generator(correction,order,energy_method)
            term=generator[0](n)
            return PerturbationTheoryExpressionEvaluator(term,term.get_poly_terms(change),change)
        return {(n,tuple(change)):self.store.compile_factory(lambda n=n,c=tuple(change):expression(n,c),
                    context=dict(source,order=n,change=tuple(change)),active_operators=active,projection=projection,**options)
                for n in orders for change in changes}

    def _frequency_request(self,initial,final):
        from .Analytic import SameStateCondition
        initial=np.asarray(initial,dtype=int);final=np.asarray(final,dtype=int)
        if initial.shape!=final.shape or initial.ndim!=1 or np.any(initial<0) or np.any(final<0):
            raise ValueError('frequency states must be nonnegative vectors of the same dimension')
        if len(initial)!=len(self.evaluator.freqs):raise ValueError('state dimension differs from expansion')
        changed=np.flatnonzero(initial!=final);projection=None;anchors=()
        if len(changed)==1 and (self.evaluator.solver.resonance_condition is None or type(self.evaluator.solver.resonance_condition) is SameStateCondition):
            if len(changed)==1:
                anchor=int(changed[0]);rest=np.delete(initial,anchor)
                if len(rest)==0 or np.all(rest==rest[0]):
                    background=int(rest[0]) if len(rest) else 0
                    projection=_OccupationProjection(int(initial[anchor]),int(final[anchor]),background)
                    anchors=(anchor,)
        return initial,final,projection,anchors

    def precompile_frequencies(self,initial,final,order=None,energy_method='auto',specialize=True,**options):
        _,_,projection,_=self._frequency_request(initial,final)
        if not specialize:projection=None
        return self.precompile('energy',order=order,energy_method=energy_method,projection=projection,**options)

    def frequencies(self,initial,final,order=None,energy_method='auto',parallelizer=None,specialize=True,**options):
        initial,final,projection,anchors=self._frequency_request(initial,final)
        numerical=dict(self.evaluator.evaluation_options);numerical.update(options)
        if not specialize or numerical.get('degenerate_changes') is not None or numerical.get('only_degenerate_terms'):
            projection=None;anchors=()
        compiled=self.precompile('energy',order=order,energy_method=energy_method,projection=projection)
        prepared=self.evaluator._prep_expansions(self.evaluator.expansions)
        identity=np.arange(len(initial))[None]
        state_perms=[(final,identity)] if projection is not None else [(initial,identity),(final,identity)]
        backend=self.evaluator.parallelizer if parallelizer is None else parallelizer
        results=[]
        for key,expression in compiled.items():
            values=expression.evaluate(state_perms,[prepared],self.evaluator.freqs,anchored_modes=anchors,
                parallelizer=backend,**numerical)
            results.append(float(values[0][0,0]) if projection is not None else float(values[1][0,0]-values[0][0,0]))
        return np.asarray(results)


# --------------------------------------------------------------------------------------------- #
#  Bloch series (formerly AnalyticBloch.py)
# Condition-based analytic Bloch wave operator, using native projected terms.
#
# Omega_0=I, D Omega_n=0 for n>0. With V=H-H0:
#   T_n = sum_{r=1}^n H_r Omega_{n-r}
#   K_n = D T_n
#   Omega_n = R_Q (T_n - sum_{r=1}^{n-1} Omega_{n-r} K_r).
# R_Q is the native conditional homological inverse; it cannot propagate through
# the retained relation. No perturbation order is singled out.
# --------------------------------------------------------------------------------------------- #

class AnalyticBlochSeries:
    def __init__(self,solver,condition=None):
        condition=SameStateCondition() if condition is None else condition
        if not condition.equivalence:
            raise ValueError('Bloch model spaces require an equivalence relation')
        self.projector=ConditionalTransformationStack(solver,condition)
        self.condition=self.projector.condition
        self.h0=solver.hamiltonian_expansion[0]
        self.expansion=solver.hamiltonian_expansion
        self.identity,self.zero=self.projector.identity,self.projector.zero
        self.wave_operator=[self.identity]
        self.effective_hamiltonian=[self.h0]

    def extend(self,order):
        if order>=len(self.expansion):raise ValueError('Hamiltonian expansion too short')
        while len(self.wave_operator)<=order:
            n=len(self.wave_operator)
            driving=operator_sum(*(operator_product(self.expansion[r],self.wave_operator[n-r])
                                  for r in range(1,n+1)))
            correction=self.projector.restrict(driving)
            feedback=operator_sum(*(operator_product(self.wave_operator[n-r],self.effective_hamiltonian[r])
                                   for r in range(1,n)))
            self.effective_hamiltonian.append(correction)
            self.wave_operator.append(self.projector.restrict(operator_sum(driving,-feedback),
                                                            retain=False,resolvent=True))
        return self

    def metric(self,order):
        self.extend(order)
        result=[self.identity]
        for n in range(1,order+1):
            result.append(self.projector.restrict(operator_sum(*(
                operator_product(ShiftedEnergyBaseline.lookup(self.wave_operator[i]),self.wave_operator[n-i])
                for i in range(n+1)))))
        return result

    def _product_series(self,left,right,order):
        return [operator_sum(*(operator_product(left[i],right[n-i])
                    for i in range(n+1) if i<len(left) and n-i<len(right))) for n in range(order+1)]

    def _metric_power(self,order,power):
        delta=self.metric(order);delta[0]=self.zero
        result=[self.identity]+[self.zero]*order
        product=[self.identity]+[self.zero]*order;coefficient=1.
        for exponent in range(1,order+1):
            product=self._product_series(product,delta,order)
            coefficient*= (power-exponent+1)/exponent
            for n in range(1,order+1):result[n]=operator_sum(result[n],coefficient*product[n])
        return result

    def hermitian_hamiltonian(self,order):
        self.extend(order)
        if isinstance(self.condition,SameStateCondition):
            # Each retained block is a scalar, so the metric similarity cancels.
            return self.effective_hamiltonian[:order+1]
        return self._product_series(self._metric_power(order,.5),
            self._product_series(self.effective_hamiltonian,self._metric_power(order,-.5),order),order)


# ---------------------------------------------------------------------------------------------
#  IndexedEvaluator (formerly VPT2/IndexedEvaluator.py)
#  Concrete indexed-block execution for analytic perturbation expressions.
#
#  This module is intentionally independent of the symbolic expression builder.
#  ``Analytic.py`` imports it lazily after an expression has been normalized, which
#  keeps the symbolic DAG and concrete numerical backend as separate layers.
#
#  The first implementation establishes the backend boundary and a bounded,
#  immutable structural plan.  It joins all coefficient/state work items produced
#  for one free-index subset before performing denominator contraction and result
#  reduction.  Later stages can replace the coefficient-row producer without
#  changing the plan or contraction interfaces defined here.
# ---------------------------------------------------------------------------------------------

_NO_SELECTION = object()


# --------------------------------------------------------------------------------------------- #
#  optional compiled ordinary harmonic denominator contraction (formerly _IndexedKernels.py)
# --------------------------------------------------------------------------------------------- #

def _energy_rows(energy_changes,frequencies):
    rows,permutations,modes=frequencies.shape
    result=np.ones((rows,permutations))
    for row in range(rows):
        for permutation in range(permutations):
            product=1.
            for factor in range(energy_changes.shape[1]):
                value=0.
                for mode in range(modes):
                    value+=frequencies[row,permutation,mode]*energy_changes[row,factor,mode+1]
                product*=value
            result[row,permutation]=product
    return result

@functools.lru_cache(maxsize=1)
def get_energy_row_kernel():
    try:
        from numba import njit
    except ImportError:
        return None
    return njit(cache=True,fastmath=False)(_energy_rows)


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
        if stats is not None:
            stats.degeneracy_pool_contexts += 1
            stats.degeneracy_pool_rows += len(permutations)

    @staticmethod
    def _mask_bytes(mask):
        return mask.nbytes if isinstance(mask, np.ndarray) else 0

    def _store(self, cache, key, entry, mask=None, size=None):
        if size is None:
            size = self._mask_bytes(mask)
        if size > self.max_cache_bytes - self.cache_bytes:
            if self.stats is not None:
                self.stats.degeneracy_pool_cache_skips += 1
            return
        cache[key] = entry
        self.cache_bytes += size
        if self.stats is not None:
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
            if self.stats is not None:
                self.stats.degeneracy_pool_predicate_hits += 1
            mask = entry[1]
        else:
            if self.stats is not None:
                self.stats.degeneracy_pool_predicate_misses += 1
            tick = time.perf_counter()
            evaluate_array = getattr(predicate, 'evaluate_array', None)
            if self.stats is not None:
                self.stats.degeneracy_states_tested += len(self.permutations)
            array_cutoff = getattr(predicate, 'scalar_state_cutoff', -1)
            if (
                    evaluate_array is None
                    or len(self.permutations) <= array_cutoff
            ):
                if self.stats is not None:
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
                if self.stats is not None:
                    self.stats.degeneracy_array_evaluations += 1
                    self.stats.degeneracy_patterns_tested += len(
                        predicate.patterns
                    )
                mask = np.asarray(
                    evaluate_array(self.permutations),
                    dtype=bool
                )
            if self.stats is not None:
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
            if self.stats is not None:
                self.stats.degeneracy_pool_selection_hits += 1
            return entry[3]

        if self.stats is not None:
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
            if self.stats is not None:
                self.stats.degeneracy_block_selection_hits += 1
            return key, entry[1], entry[2]
        if self.stats is not None:
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
            if self.stats is not None:
                self.stats.degeneracy_block_selection_cache_skips += 1


@dataclasses.dataclass(frozen=True)
class IndexedEnergyTerm:
    """One polynomial numerator and its denominator change product."""

    ordinal: int
    energy_changes: tuple
    polynomial: object
    contraction_shape: tuple
    projected: bool = False
    constant_numerator: object = None


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
    index_bindings: object = None


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
        'contributions', 'reference_states'
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
        self.reference_states = []

    def append(self, state_index, mask_positions, global_permutation_indices,
               prefactors, expression_key,
               plan, permutation_subsets, permutation_substates, tuple_states,
               permutation_frequencies, state_block_identity, contribution,
               reference_state=None):
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
        self.reference_states.append(reference_state)
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
    batched_state_blocks: int = 0
    batched_source_states: int = 0
    unmasked_selection_fast_paths: int = 0
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
    coefficient_gather_cache_hits: int = 0
    coefficient_gather_cache_misses: int = 0
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
    free_index_source_plans: int = 0
    free_index_canonical_plans: int = 0
    free_index_shared_plans: int = 0
    free_index_bound_blocks: int = 0
    free_index_fallback_calls: int = 0
    shared_value_input_blocks: int = 0
    shared_value_output_blocks: int = 0
    shared_value_blocks_eliminated: int = 0
    shared_value_input_points: int = 0
    shared_value_unique_points: int = 0
    shared_value_duplicate_points: int = 0
    shared_value_disjoint_blocks: int = 0
    shared_value_pool_hits: int = 0
    shared_value_pool_misses: int = 0
    shared_value_sort_calls: int = 0
    shared_value_seconds: float = 0.0
    max_shared_value_permutations: int = 0
    energy_scalar_lanes: int = 0
    energy_computed_scalar_lanes: int = 0
    energy_rows_computed: int = 0
    energy_rows_reused: int = 0
    energy_row_cache_evictions: int = 0
    max_energy_row_cache_bytes: int = 0

    def as_dict(self):
        return dataclasses.asdict(self)




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

    preparation_context: object = dataclasses.field(default=None, repr=False, compare=False)
    _tensor_metadata: object = dataclasses.field(init=False,repr=False,compare=False)
    _gather_cache: object = dataclasses.field(init=False,repr=False,compare=False)

    use_cached_request_patterns = True

    def __post_init__(self):
        object.__setattr__(
            self,
            'request_pattern_cache',
            BoundedCache(
                self.max_pattern_cache_items,
                self.max_pattern_cache_bytes,
                sizer=self._request_pattern_size,
                lru=True
            )
        )

        context = self.preparation_context
        metadata = {}
        for product in self.coefficient_tensors:
            for expansion in product:
                for tensor in expansion:
                    identity = id(tensor)
                    if identity not in metadata:
                        scalar = np.ndim(tensor)==0
                        metadata[identity] = (scalar,None if scalar or context is None else context.bind(tensor))
        cache = None if context is None else context.cache('coefficient_gathers',32*1024**2)
        if cache is not None and (cache.max_items==0 or cache.max_bytes==0):
            cache = None
        object.__setattr__(self,'_tensor_metadata',types.MappingProxyType(metadata))
        object.__setattr__(self,'_gather_cache',cache)

    @staticmethod
    def _request_pattern_size(key, value):
        unique_patterns, inverse = value
        return (unique_patterns.nbytes + inverse.nbytes
                + sum(len(pattern) for pattern in key) * np.dtype(np.intp).itemsize)

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
            if stats is not None:
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
            self.request_pattern_cache[key] = (unique_patterns, inverse)
        else:
            if stats is not None:
                stats.coefficient_pattern_cache_hits += 1
            unique_patterns, inverse = entry
        cache = self.request_pattern_cache
        if stats is not None:
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
    def _gather_tensor_requests(permutations, tensor, unique_patterns, inverse,
                                *, cache=None, tensor_binding=None, stats=None):
        """Gather distinct tensor factors, optionally reusing the plan binding."""
        tensor_indices = permutations[:,unique_patterns]
        key = gathered = None
        if cache is not None:
            from .Analytic import _HashableNumPyIndex
            if tensor_binding is None:
                raise ValueError('cached tensor gathering requires a plan-bound tensor')
            key = (tensor_binding,_HashableNumPyIndex(tensor_indices))
            gathered = cache.get(key)
        if gathered is None:
            if stats is not None and cache is not None:
                stats.coefficient_gather_cache_misses += 1
            gathered = tensor[tuple(tensor_indices[...,axis]
                                    for axis in range(tensor_indices.shape[-1]))]
            if cache is not None:
                gathered.setflags(write=False)
                cache[key] = gathered
        elif stats is not None:
            stats.coefficient_gather_cache_hits += 1
        # Advanced indexing owns the returned request plane; callers cannot
        # mutate the cached factor while accumulating a different product.
        return np.moveaxis(gathered,0,1)[inverse]

    def _evaluate_product_chunk(self, permutations, product_indices,
                                remapped_products, subset, num_fixed,
                                zero_cutoff, stats, independent_coefficients=False):
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
                    if len(indices) == 0 or self._tensor_metadata[id(tensor)][0]:
                        scalar_requests.append((local_index, tensor, len(indices)))
                    else:
                        tensor_key = id(tensor)
                        entry = request_groups.get(tensor_key)
                        if entry is None:
                            entry = [tensor, [], [], [], self._tensor_metadata[tensor_key][1]]
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
                        # Preserve legacy coupled planes by default. Batched
                        # independent properties retain each scalar call's cutoff.
                        if independent_coefficients:
                            good[local_index,expansion_index,:] = False
                        else:
                            good[local_index,:,:] = False
                    else:
                        products[local_index, expansion_index, :] *= base_value

                for tensor, destinations, patterns, structural_patterns, tensor_binding \
                        in request_groups.values():
                    if stats is not None:
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
                        if stats is not None:
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
                        if stats is not None:
                            stats.coefficient_pattern_compilations += 1
                        unique_patterns, inverse = np.unique(
                            np.asarray(patterns, dtype=np.intp),
                            axis=0,
                            return_inverse=True
                        )
                        if stats is not None:
                            stats.coefficient_unique_patterns += len(
                                unique_patterns
                            )
                    values = self._gather_tensor_requests(
                        permutations,tensor,unique_patterns,inverse,
                        cache=self._gather_cache,tensor_binding=tensor_binding,stats=stats)
                    destinations = np.asarray(destinations, dtype=int)
                    active = good[destinations, expansion_index, :]
                    block = products[destinations, expansion_index, :]
                    np.multiply(block, values, out=block, where=active)
                    products[destinations, expansion_index, :] = block
                    good[destinations, expansion_index, :] = np.logical_and(
                        active,
                        np.abs(block) >= zero_cutoff
                    )

        if independent_coefficients:
            products = np.where(good,products,0.)
        return products, np.any(good, axis=1)

    def gather(self, permutations, subset, num_fixed,
               counts_cache, factorials, zero_cutoff, stats, *, independent_coefficients=False):
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
            if stats is not None:
                stats.coefficient_product_chunks += 1
            products, good = self._evaluate_product_chunk(
                permutations,
                chunk_products,
                chunk_remapped,
                subset,
                num_fixed,
                zero_cutoff,
                stats,
                independent_coefficients=independent_coefficients
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

        if stats is not None:
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
    use_batched_state_routing = True
    max_batched_states = 128
    use_unmasked_selection = True
    use_jit_energy_contractions = False
    use_constant_numerators = True
    use_compiled_permutation_routing = True
    use_shared_degeneracy_pool = True
    use_compiled_degeneracy_selection = True
    # keyed by expression identity; entries keep the expression alive so ids are never reused
    _plan_cache = BoundedCache(default_plan_cache_size)

    def __init__(self, evaluator_type, expression, chunk_size=None, preparation_context=None):
        self.preparation_context = preparation_context
        self.evaluator_type = evaluator_type
        self.expression = expression
        self.chunk_size = (
            self.default_chunk_size if chunk_size is None else max(1, int(chunk_size))
        )
        self.stats = IndexedEvaluationStats() if Diagnostics.enabled else None
        self._coefficient_plans = {}
        self._degeneracy_selection_plans = collections.OrderedDict()
        self._degeneracy_group_plans = collections.OrderedDict()
        self._started = time.perf_counter()
        self._free_index_plan = None
        self._native_plan = None
        data=getattr(expression,'_free_index_plan_data',None)
        if _current_free_index_reduction.get() and data is not None:
            if data.get('format')!=1:raise ValueError('unsupported free-index plan format')
            self._free_index_plan=self.get_plan(expression,reduced=True)
            if self.stats is not None:
                stats=data['stats']
                self.stats.free_index_source_plans=stats['source_keys']
                self.stats.free_index_canonical_plans=stats['canonical_plans']
                self.stats.free_index_shared_plans=stats['shared_plans']
        self.plan=self._free_index_plan
        if self.plan is None:
            self._native_plan=self.plan=self.get_plan(expression)

    @classmethod
    def clear_plan_cache(cls):
        cls._plan_cache.clear()

    def get_plan(self, expression, reduced=False):
        namespace='free-index-v1' if reduced else None
        cache_plans = (self.preparation_context is None
                       or getattr(self.preparation_context, 'immutable_expressions', False))
        entry = self._plan_cache.get((id(expression), namespace)) if cache_plans else None
        plan = entry[1] if entry is not None and entry[0] is expression else None
        if plan is not None:
            if self.stats is not None:
                self.stats.plan_cache_hits += 1
            return plan
        if self.stats is not None:
            self.stats.plan_cache_misses += 1
        plan = self.compile_plan(expression,reduced=reduced)
        if self.stats is not None:
            self.stats.plans_compiled += 1
            self.stats.plan_terms += plan.term_count
        if cache_plans:
            self._plan_cache[(id(expression), namespace)] = (expression, plan)
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
                    contraction_shape,
                    any(getattr(row[0], '_state_projection_factor', False) for row in energy_changes),
                    (np.asarray(polynomial.evaluate_polynomial(
                        np.zeros((1,polynomial.ndim),dtype=int))).reshape(-1)[0]
                     if isinstance(polynomial,analytic.PolyPath) and not any(polynomial.order)
                     else None)
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
    def compile_plan(cls, expression, reduced=False):
        # Lazy import avoids an Analytic -> IndexedEvaluator -> Analytic cycle.
        from . import Analytic as analytic

        subexpressions = {}
        term_count = 0
        data=getattr(expression,'_free_index_plan_data',None) if reduced else None
        terms=expression.terms if data is None else data['values']
        interned={}
        for key, subexpression in terms.items():
            # Compiler-interned canonical values share the exact native plan.
            identity=id(subexpression)
            if identity not in interned:
                interned[identity]=cls._compile_subexpression(subexpression, analytic)
            subplan, subterms = interned[identity]
            subexpressions[key] = subplan
            term_count += subterms
        return IndexedExpressionPlan(
            expression,
            types.MappingProxyType(subexpressions),
            term_count,
            None if data is None else types.MappingProxyType(data['bindings'])
        )

    def finish(self):
        if self.stats is not None:
            self.stats.elapsed_seconds = time.perf_counter() - self._started
            return self.stats.as_dict()
        return None

    def get_coefficient_plan(self, coefficient_specs, coefficient_tensors):
        if self.preparation_context is not None:
            def create():
                return IndexedCoefficientGatherPlan(coefficient_specs, coefficient_tensors,
                    0 if not coefficient_tensors else len(coefficient_tensors[0]),
                    preparation_context=self.preparation_context)
            plan, hit = self.preparation_context.coefficient_plan(
                coefficient_specs, coefficient_tensors, create)
            if hit:
                if self.stats is not None:
                    self.stats.coefficient_plan_hits += 1
            else:
                if self.stats is not None:
                    self.stats.coefficient_plan_misses += 1
                    self.stats.coefficient_plans_compiled += 1
            return plan
        key = (id(coefficient_specs), id(coefficient_tensors))
        entry = self._coefficient_plans.get(key)
        if (
            entry is not None
            and entry[0] is coefficient_specs
            and entry[1] is coefficient_tensors
        ):
            if self.stats is not None:
                self.stats.coefficient_plan_hits += 1
            return entry[2]
        if self.stats is not None:
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
        if self.stats is not None:
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
            if self.stats is not None:
                self.stats.degeneracy_selection_plan_hits += 1
            self._degeneracy_selection_plans.move_to_end(key)
            return entry[3]

        if self.stats is not None:
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
        if self.stats is not None:
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
            if self.stats is not None:
                self.stats.degeneracy_selection_plan_evictions += 1
        if self.stats is not None:
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
            if self.stats is not None:
                self.stats.degeneracy_group_plan_hits += 1
            self._degeneracy_group_plans.move_to_end(key)
            return entry[1]

        if self.stats is not None:
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
            if self.stats is not None:
                self.stats.degeneracy_group_plan_evictions += 1
        if self.stats is not None:
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
        if self.stats is not None:
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
        if self.use_unmasked_selection and degenerate_changes is None:
            use_left,use_right,join=only_degenerate_terms
            selected=join(not use_left,not use_right)
            if isinstance(selected,(bool,np.bool_)):
                if self.stats is not None:
                    self.stats.unmasked_selection_fast_paths += 1
                return (tuple((shape,None,terms) for shape,terms in blocks.plans[block_index].groups)
                        if selected else ())
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
                    if self.stats is not None:
                        self.stats.degeneracy_selection_hits += 1
                    selected = entry[2]
                    selection_id = entry[3]
                else:
                    if self.stats is not None:
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
            if self.stats is not None:
                self.stats.degeneracy_mask_misses += 1
            tick = time.perf_counter()
            evaluate = getattr(predicate, 'evaluate', None)
            if self.stats is not None:
                self.stats.degeneracy_states_tested += len(permutation_subsets)
            array_cutoff = getattr(
                predicate,
                'scalar_state_cutoff',
                -1
            )
            if evaluate is None or len(permutation_subsets) <= array_cutoff:
                if self.stats is not None:
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
                if self.stats is not None:
                    self.stats.degeneracy_array_evaluations += 1
                    self.stats.degeneracy_patterns_tested += len(
                        predicate.patterns
                    )
                mask = np.asarray(
                    evaluate(permutation_subsets),
                    dtype=bool
                )
            if self.stats is not None:
                self.stats.degeneracy_seconds += time.perf_counter() - tick
            predicate_masks[predicate_key] = (predicate, mask)
        else:
            if self.stats is not None:
                self.stats.degeneracy_mask_hits += 1
            mask = entry[1]
        return mask if only_degenerate else np.logical_not(mask)

    @classmethod
    def _contract_energy_rows(cls, energy_changes, frequency_blocks):
        """Contract ordinary harmonic denominator rows; projected masks stay separate."""
        if np.size(energy_changes) == 0:
            return np.ones(np.shape(frequency_blocks)[:-1], dtype=np.result_type(frequency_blocks, float))
        if (cls.use_jit_energy_contractions
                and frequency_blocks.dtype == np.dtype('float64')
                and energy_changes.dtype.kind in 'biuf'):
            kernel=get_energy_row_kernel()
            if kernel is not None:
                return kernel(np.ascontiguousarray(energy_changes,dtype=float),
                              np.ascontiguousarray(frequency_blocks))
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
        if self.stats is not None:
            self.stats.max_polynomial_scratch_bytes = max(
                self.stats.max_polynomial_scratch_bytes,
                polynomial_block.nbytes
            )
        return polynomial_block

    # Bound the joined permutation axis as well as the existing row chunk.
    # This remains a diagnostic native control, independent of compiler coalescing.
    use_shared_canonical_values = True
    use_shared_energy_rows = True
    max_energy_row_cache_bytes = 32 * 1024**2
    max_shared_value_permutations = 4096

    def _share_canonical_value_blocks(self, blocks):
        """Evaluate only genuine overlapping physical inputs once.

        Distinct ordered summation sectors must retain their native blocks.
        Joining them only widens vectors without removing scalar work and can
        destroy useful numerator-cache locality. Exact overlap is checked once
        per cached physical pool; every original coefficient/scatter survives.
        """
        from .Analytic import _StatePermutationBlockIdentity
        tick = time.perf_counter()

        def component(index):
            return (id(blocks.permutation_frequencies[index]),
                    id(blocks.permutation_substates[index]),
                    id(blocks.state_block_identities[index]))

        groups = collections.OrderedDict()
        for index, plan in enumerate(blocks.plans):
            groups.setdefault((id(plan), blocks.state_indices[index]), []).append(index)
        cap = max(1, int(self.max_shared_value_permutations))
        chunks = []
        for indices in groups.values():
            pending, size, components = [], 0, set()
            for index in indices:
                key = component(index)
                count = 0 if key in components else len(blocks.permutation_subsets[index])
                if pending and size + count > cap:
                    chunks.append(pending)
                    pending, size, components = [], 0, set()
                    count = len(blocks.permutation_subsets[index])
                pending.append(index)
                size += count
                components.add(key)
            if pending:
                chunks.append(pending)
        if self.stats is not None:
            self.stats.shared_value_input_blocks += len(blocks)
        if len(chunks) == len(blocks):
            if self.stats is not None:
                self.stats.shared_value_output_blocks += len(blocks)
                self.stats.shared_value_seconds += time.perf_counter() - tick
            return blocks, None, None

        decisions, pools = [], {}
        for indices in chunks:
            first = indices[0]
            components = {component(index): index for index in indices}
            keys = tuple(sorted(components))
            pool_key = (blocks.state_indices[first], keys)
            pooled = pools.get(pool_key)
            if pooled is None:
                if self.stats is not None:
                    self.stats.shared_value_pool_misses += 1
                if len(keys) == 1:
                    source = components[keys[0]]
                    pooled = (blocks.permutation_subsets[source], blocks.permutation_substates[source],
                              blocks.tuple_states[source], blocks.permutation_frequencies[source],
                              blocks.state_block_identities[source], {keys[0]: None})
                else:
                    ordered = [components[key] for key in keys]
                    physical = np.concatenate([blocks.permutation_subsets[index] for index in ordered], axis=0)
                    unique, take, inverse = np.unique(physical, axis=0, return_index=True, return_inverse=True)
                    if self.stats is not None:
                        self.stats.shared_value_sort_calls += 1
                    substates = np.concatenate([blocks.permutation_substates[index] for index in ordered], axis=0)[take]
                    frequencies = np.concatenate([blocks.permutation_frequencies[index] for index in ordered], axis=0)[take]
                    tuples = [entry for index in ordered for entry in blocks.tuple_states[index]]
                    tuples = [tuples[int(index)] for index in take]
                    inverse_routes, cursor = {}, 0
                    for key, index in zip(keys, ordered):
                        count = len(blocks.permutation_subsets[index])
                        inverse_routes[key] = inverse[cursor:cursor+count]
                        cursor += count
                    pooled = (unique, substates, tuples, frequencies,
                              _StatePermutationBlockIdentity(tuples), inverse_routes)
                pools[pool_key] = pooled
            else:
                if self.stats is not None:
                    self.stats.shared_value_pool_hits += 1
            incoming = sum(len(blocks.permutation_subsets[index]) for index in indices)
            unique_count = len(pooled[0])
            if self.stats is not None:
                self.stats.shared_value_input_points += incoming
                self.stats.shared_value_unique_points += unique_count
                self.stats.shared_value_duplicate_points += incoming - unique_count
            if unique_count == incoming:
                if self.stats is not None:
                    self.stats.shared_value_disjoint_blocks += len(indices)
                decisions.extend(([index], None) for index in indices)
            else:
                decisions.append((indices, pooled))
        if self.stats is not None:
            self.stats.shared_value_output_blocks += len(decisions)
            self.stats.shared_value_blocks_eliminated += len(blocks) - len(decisions)
        if len(decisions) == len(blocks):
            if self.stats is not None:
                self.stats.shared_value_seconds += time.perf_counter() - tick
            return blocks, None, None

        joined = ConcreteBlockTable()
        routes, aliases = [None] * len(blocks), []
        for indices, pooled in decisions:
            first = indices[0]
            if pooled is None:
                pooled = (blocks.permutation_subsets[first], blocks.permutation_substates[first],
                          blocks.tuple_states[first], blocks.permutation_frequencies[first],
                          blocks.state_block_identities[first], {component(first): None})
            physical, substates, tuples, frequencies, identity, inverse_routes = pooled
            position = np.arange(len(physical), dtype=np.intp)
            group = joined.append(blocks.state_indices[first], position, position, None,
                                  blocks.expression_keys[first], blocks.plans[first], physical,
                                  substates, tuples, frequencies, identity,
                                  np.zeros((blocks.contributions[first].shape[0], len(physical))),
                                  blocks.reference_states[first])
            aliases.append(tuple(blocks.expression_keys[index] for index in indices))
            for index in indices:
                routes[index] = (group, inverse_routes[component(index)])
            if self.stats is not None:
                self.stats.max_shared_value_permutations = max(self.stats.max_shared_value_permutations, len(physical))
        if self.stats is not None:
            self.stats.shared_value_seconds += time.perf_counter() - tick
        return joined, routes, aliases

    def _evaluate_work_items(self, blocks, output, change, baseline_shift,
                             pows, polynomial_cache, degenerate_changes,
                             only_degenerate_terms, zero_cutoff, logger,
                             evaluation_context=None,
                             degeneracy_context=None, uniform_output=None):
        evaluator = self.evaluator_type
        original_blocks=blocks;shared_routes=shared_aliases=None
        if (self.use_shared_canonical_values and self.plan.index_bindings is not None
                and degenerate_changes is None and not any(only_degenerate_terms[:2])):
            blocks,shared_routes,shared_aliases=self._share_canonical_value_blocks(blocks)
        row_groups = collections.OrderedDict()
        energy_row_cache = collections.OrderedDict()
        energy_row_cache_bytes = 0
        frequency_keys = collections.OrderedDict()
        frequency_key_bytes = 0
        projector_domains = {}

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
                if terms and not terms[0].energy_changes:
                    # A unit factor is a bare numerator: it has no energy
                    # contraction, resonance selection, or extra cutoff.
                    if not any(only_degenerate_terms[:2]):
                        tick = time.perf_counter()
                        for term in terms:
                            polynomial = evaluator._eval_poly(
                                polynomial_cache,
                                blocks.tuple_states[block_index],
                                blocks.permutation_substates[block_index],
                                pows, term.polynomial, change, baseline_shift,
                                False, logger,
                                block_key=('full', blocks.state_block_identities[block_index]))
                            blocks.contributions[block_index] += polynomial
                            if self.stats is not None:
                                self.stats.polynomial_evaluations += 1
                        if self.stats is not None:
                            self.stats.polynomial_seconds += time.perf_counter() - tick
                    continue
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

        if self.stats is not None:
            self.stats.contraction_groups += len(row_groups)
        for segments in row_groups.values():
            if self.stats is not None:
                self.stats.contraction_segments += len(segments.block_indices)
                self.stats.max_pending_segments = max(
                    self.stats.max_pending_segments,
                    len(segments.block_indices)
                )
            for block_indices, selected_rows, terms, segment_ids in segments.iter_blocks(
                    self.chunk_size, return_segment_ids=True):
                block_size = len(terms)
                if self.stats is not None:
                    self.stats.contraction_blocks += 1
                    self.stats.contraction_rows += block_size
                    self.stats.energy_scalar_lanes += block_size*(
                        len(blocks.permutation_subsets[int(block_indices[0])])
                        if selected_rows[0] is None else len(selected_rows[0]))
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
                projected_call = any(term.projected for term in terms)

                def contract(indices):
                    if projected_call:
                        from .Analytic import projected_energy_rows
                        permutation_rows, source_rows = [], []
                        for index in indices:
                            block_index = int(block_indices[index])
                            permutations = blocks.permutation_subsets[block_index]
                            selected = selected_rows[index]
                            if selected is not None:
                                permutations = np.asarray(permutations)[selected]
                            permutation_rows.append(permutations)
                            source_rows.append(blocks.reference_states[block_index])
                        mask, factors = projected_energy_rows(
                            tuple(terms[index].energy_changes for index in indices),
                            source_rows, permutation_rows,
                            np.stack([selections[index][0] for index in indices]), baseline_shift,
                            row_projection_bindings=[
                                projection_inputs(evaluation_context, blocks.state_indices[int(block_indices[index])])
                                for index in indices],
                            preparation_context=self.preparation_context)
                        return mask, np.where(mask, factors, 1.)
                    changes = np.asarray([terms[index].energy_changes for index in indices])
                    frequencies = np.stack([selections[index][0] for index in indices])
                    if self.preparation_context is None:
                        return None, self._contract_energy_rows(changes, frequencies)
                    from .Analytic import _HashableNumPyIndex
                    cache = self.preparation_context.cache('ordinary_energy',32*1024**2)
                    keys = [(_HashableNumPyIndex(c),_HashableNumPyIndex(f)) for c,f in zip(changes,frequencies)]
                    ready, missing = {}, {}
                    for index,key in enumerate(keys):
                        value = cache.get(key)
                        if value is not None:
                            ready[key] = value
                        elif key not in missing:
                            missing[key] = index
                    if missing:
                        selected = list(missing.values())
                        values = self._contract_energy_rows(changes[selected],frequencies[selected])
                        for key,value in zip(missing,values):
                            value = value.copy()
                            value.setflags(write=False)
                            ready[key] = value
                            cache[key] = value
                    return None, np.stack([ready[key] for key in keys])

                share_energy = (self.use_shared_energy_rows and self.plan.index_bindings is not None
                    and degenerate_changes is None and not any(only_degenerate_terms[:2])
                    and not any(getattr(row[0], '_conditional_endpoint', False)
                                for term in terms for row in term.energy_changes))
                if not share_energy:
                    projected_mask, energy_factors = contract(range(block_size))
                    if self.stats is not None:
                        self.stats.energy_rows_computed += block_size
                        self.stats.energy_computed_scalar_lanes += energy_factors.size
                else:
                    # A batch-local cache: the formal row, physical frequency/
                    # occupation binding, source state and projection route must
                    # all agree. Numerators and coefficient products are free to
                    # differ. Runtime inputs are fixed for this call and no
                    # array values are rescanned or globally cached.
                    keys, ready, missing = [], {}, collections.OrderedDict()
                    for index, term in enumerate(terms):
                        block_index = int(block_indices[index])
                        # Frequency bindings describe physical mode indices,
                        # not the source-state row that happened to gather them.
                        frequency_key = frequency_keys.get(block_index)
                        if frequency_key is None:
                            physical = np.asarray(blocks.permutation_subsets[block_index])
                            frequency_key = (physical.shape,physical.dtype.str,physical.tobytes())
                            capacity = max(0,int(self.max_energy_row_cache_bytes))
                            size = len(frequency_key[2])
                            if size <= capacity:
                                while frequency_keys and (len(frequency_keys)>=128 or frequency_key_bytes+size>capacity):
                                    _,removed = frequency_keys.popitem(last=False)
                                    frequency_key_bytes -= len(removed[2])
                                frequency_keys[block_index] = frequency_key
                                frequency_key_bytes += size
                        else:
                            frequency_keys.move_to_end(block_index)
                        source = blocks.state_indices[block_index] if term.projected else None
                        if term.projected:
                            if source not in projector_domains:
                                bindings = projection_inputs(evaluation_context,source)
                                state = np.asarray(blocks.reference_states[block_index])[0]
                                projector_domains[source] = all(binding[3] and
                                    (binding[2] is None or np.array_equal(binding[2],state))
                                    for binding in bindings)
                            if projector_domains[source]:
                                # Equal local occupations still remain in the
                                # selection key: valid-state masks can differ.
                                source = None
                        key = (term.energy_changes, projected_call, source,
                               selections[index][3], frequency_key)
                        keys.append(key)
                        entry = energy_row_cache.get(key)
                        if entry is not None:
                            energy_row_cache.move_to_end(key)
                            ready[key] = entry
                        elif key not in missing:
                            missing[key] = index
                    if missing:
                        indices = list(missing.values())
                        masks, factors = contract(indices)
                        if self.stats is not None:
                            self.stats.energy_rows_computed += len(indices)
                            self.stats.energy_computed_scalar_lanes += factors.size
                        capacity = max(0, int(self.max_energy_row_cache_bytes))
                        for position, key in enumerate(missing):
                            # Own bounded row arrays; a view would retain an
                            # entire parent chunk while accounting only one row.
                            entry = (None if masks is None else masks[position].copy(),
                                     factors[position].copy())
                            ready[key] = entry
                            size = entry[1].nbytes + (0 if entry[0] is None else entry[0].nbytes) + len(key[-1][2])
                            if size <= capacity:
                                while energy_row_cache and energy_row_cache_bytes + size > capacity:
                                    removed_key, removed = energy_row_cache.popitem(last=False)
                                    energy_row_cache_bytes -= (removed[1].nbytes +
                                        (0 if removed[0] is None else removed[0].nbytes) + len(removed_key[-1][2]))
                                    if self.stats is not None:
                                        self.stats.energy_row_cache_evictions += 1
                                energy_row_cache[key] = entry
                                energy_row_cache_bytes += size
                                if self.stats is not None:
                                    self.stats.max_energy_row_cache_bytes = max(
                                        self.stats.max_energy_row_cache_bytes, energy_row_cache_bytes)
                    if self.stats is not None:
                        self.stats.energy_rows_reused += block_size - len(missing)
                    energy_factors = np.stack([ready[key][1] for key in keys])
                    projected_mask = np.stack([ready[key][0] for key in keys]) if projected_call else None
                if self.stats is not None:
                    self.stats.energy_seconds += time.perf_counter() - tick

                tick = time.perf_counter()
                polynomial_factors = [
                    term.constant_numerator
                    if (self.use_constant_numerators and term.constant_numerator is not None
                        and change is not None and len(change)==0 and baseline_shift is None) else
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
                if self.stats is not None:
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
                if projected_mask is not None:
                    polynomial_block = np.where(projected_mask[:, None, :], polynomial_block, 0.)
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
                            '_state_pair_context': lambda: (
                                blocks.reference_states[block_index],
                                blocks.permutation_subsets[block_index] if selected_rows[row_index] is None else
                                np.asarray(blocks.permutation_subsets[block_index])[selected_rows[row_index]],
                                selection[0], baseline_shift,
                                projection_inputs(
                                    evaluation_context, blocks.state_indices[block_index])),
                            'coefficient_key': blocks.expression_keys[block_index],
                            'shared_coefficient_keys': (None if shared_aliases is None else shared_aliases[block_index]),
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
                    row_diagnostics=row_diagnostics, logger=logger
                )
                segment_starts = np.concatenate((
                    np.array([0], dtype=np.intp),
                    np.flatnonzero(segment_ids[1:] != segment_ids[:-1]) + 1
                ))
                reduced_rows = np.add.reduceat(
                    scaled_rows, segment_starts, axis=0
                )
                if self.stats is not None:
                    self.stats.reduction_scatters += len(segment_starts)
                for reduced_index, row_index in enumerate(segment_starts):
                    block_index = int(block_indices[row_index])
                    selected = selected_rows[row_index]
                    scaled = reduced_rows[reduced_index]
                    if selected is None:
                        blocks.contributions[block_index] += scaled
                    else:
                        blocks.contributions[block_index][:, selected] += scaled
                if self.stats is not None:
                    self.stats.reduction_seconds += time.perf_counter() - tick

        if shared_routes is not None:
            # Apply the native prefactor once before routing one shared value
            # back to each original coefficient block in its original order.
            for index,plan in enumerate(blocks.plans):
                if not plan.direct_polynomial:blocks.contributions[index]*=plan.prefactor
        for block_index in range(len(original_blocks)):
            if shared_routes is None:
                contribution=blocks.contributions[block_index]
                plan=blocks.plans[block_index]
                if not plan.direct_polynomial:contribution*=plan.prefactor
            else:
                group,inverse=shared_routes[block_index]
                contribution=blocks.contributions[group]
                if inverse is not None:contribution=contribution[:,inverse]
            value = (
                original_blocks.prefactors[block_index][:, np.newaxis, :]
                * contribution[np.newaxis]
            )
            state_indices=original_blocks.state_indices[block_index]
            if isinstance(state_indices,tuple):
                if uniform_output is not None:
                    # A routed batch has one shared permutation and distinct
                    # target states. Preserve coefficient/term accumulation order.
                    uniform_output[np.asarray(state_indices),:,0] += value[:,:,0].T
                else:
                    for row,state_index in enumerate(state_indices):
                        output[state_index][:,original_blocks.mask_positions[block_index]] += value[:,row,:]
            else:
                output[state_indices][:,original_blocks.mask_positions[block_index]] += value[:,0,:]

    def _append_routed_work_item(
            self, work_items, state_index, mask_positions, value_positions,
            global_permutation_indices,
            expression_key, state_permutations, prefactors, frequencies,
            num_fixed, subset, full_set, state_permutation_cache,
            degenerate_changes):
        evaluator = self.evaluator_type
        if isinstance(state_index,tuple):
            state_block=np.asarray([state_permutations[i][0] for i in state_index])
            permutations=state_permutations[state_index[0]][1]
        else:
            state,permutations=state_permutations[state_index]
            state_block=state[np.newaxis]
        block_prefactors = prefactors[value_positions].T
        binding=(None if self.plan.index_bindings is None else
                 self.plan.index_bindings.get(expression_key))
        if binding is not None and binding!=tuple(range(len(binding))):
            # Gather original physical modes in canonical polynomial order.
            # Coefficient prefactors and output destinations remain original.
            full_set=tuple(full_set[i] for i in binding)
            if self.stats is not None:
                self.stats.free_index_bound_blocks += 1
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
            np.zeros((len(state_block), len(permutation_subsets))),
            reference_state=state_block
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
        if (self.use_batched_state_routing and degenerate_changes is None
                and len(state_permutations)>1 and len(all_permutations)==1
                and all(len(permutations)==1 for _,permutations in state_permutations)
                and np.all(np.asarray(permutation_map)==0)
                and np.all(np.asarray(eval_permutations)==0)
                and np.all(np.diff(eval_coefficients)>0)):
            plans=[self.plan.subexpressions[coefficient_specs[i]] for i in eval_coefficients]
            # Conditional/finite-state endpoint masks currently bind one source
            # row. Keep their original evaluator unless all terms prove ordinary.
            if all(not any(term.projected for _,terms in plan.groups for term in terms)
                   for plan in plans):
                work_items=ConcreteBlockTable()
                single=np.array([0],dtype=np.intp)
                maximum=max(1,int(self.max_batched_states))
                for position,coefficient_group in enumerate(eval_coefficients):
                    for first in range(0,len(state_permutations),maximum):
                        indices=tuple(range(first,min(first+maximum,len(state_permutations))))
                        self._append_routed_work_item(work_items,indices,single,
                            np.array([position],dtype=np.intp),single,
                            coefficient_specs[coefficient_group],state_permutations,
                            prefactors,frequencies,num_fixed,subset,full_set,
                            state_permutation_cache,degenerate_changes)
                        if self.stats is not None:
                            self.stats.batched_state_blocks += 1
                            self.stats.batched_source_states += len(indices)
                if self.stats is not None:
                    self.stats.routing_matches += len(eval_permutations)*len(state_permutations)
                return work_items
        routing_plan = IndexedPermutationRoutingPlan.compile(
            permutation_map,
            [len(permutations) for _, permutations in state_permutations],
            len(all_permutations)
        )
        if self.stats is not None:
            self.stats.routing_plans_compiled += 1
        routes = routing_plan.route(eval_permutations, eval_coefficients)
        if self.stats is not None:
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
            logger.log_print("Gathering indexed coefficients for subset {s}...", s=subset)
            coefficient_start = time.perf_counter()

        if self.stats is not None:
            self.stats.calls += 1
        if self._free_index_plan is not None:
            data=expression._free_index_plan_data
            if degenerate_changes is None and num_fixed==data['num_fixed']:
                self.plan=self._free_index_plan
            else:
                # External legacy predicates bind original energy-row axes.
                # Keep the entire native path until their transport is proven.
                if self.stats is not None:
                    self.stats.free_index_fallback_calls += 1
                if self._native_plan is None:self._native_plan=self.get_plan(expression)
                self.plan=self._native_plan
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
            self.stats,
            independent_coefficients=getattr(evaluation_context,'projection_bindings',{}).get(
                'independent_coefficients',False)
        )
        if verbose:
            logger.log_print("Gathered {n} coefficient rows; took {e:.3f}s",
                             n=len(eval_permutations), e=time.perf_counter()-coefficient_start)
        uniform_output=None
        if (self.use_batched_state_routing and len(state_permutations)>1
                and all(len(permutations)==1 for _,permutations in state_permutations)):
            uniform_output=np.zeros((len(state_permutations),prefactors.shape[1],1))
            output=[uniform_output[i] for i in range(len(state_permutations))]
        else:
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
        routing_seconds = time.perf_counter() - tick
        if self.stats is not None:
            self.stats.routing_seconds += routing_seconds
        if verbose:
            logger.log_print("Routed {n} indexed blocks; took {e:.3f}s",
                             n=len(work_items), e=routing_seconds)

        if self.stats is not None:
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
        if verbose:
            logger.log_print("Evaluating {n} indexed blocks...", n=len(work_items))
            evaluation_start = time.perf_counter()
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
            degeneracy_context=degeneracy_context,
            uniform_output=uniform_output
        )
        if verbose:
            logger.log_print("Evaluated indexed blocks; took {e:.3f}s",
                             e=time.perf_counter()-evaluation_start)
        return output
