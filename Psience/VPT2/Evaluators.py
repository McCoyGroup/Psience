"""
Numerical evaluation backends for the analytic VPT expressions in `Analytic`.
"""
# imports from merged modules
import collections
from collections.abc import Mapping
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

_current_evaluation_mode = contextvars.ContextVar('analytic_evaluation_mode', default=None)


def _validate_evaluation_mode(mode, *, nmodes=None):
    # Resolve the automatic policy once, when the active mode count is known
    # during construction. None in an evaluation context means inheritance.
    if nmodes is not None:
        if mode is None:
            return 'actions' if nmodes > 6 else 'direct'
        if mode not in {'direct', 'actions'}:
            raise ValueError("evaluation_mode must be 'direct', 'actions', or None during construction")
    # Retain the existing low-level backend controls for expression evaluation.
    if mode not in {None, 'direct', 'actions', 'compiled', 'dag', 'dag_legacy',
                    'materialized', 'materialized_cached', 'materialized_legacy'}:
        raise ValueError("unknown analytic evaluation mode {!r}".format(mode))
    return mode


@contextlib.contextmanager
def evaluation_mode_context(mode):
    """Select a controller evaluation policy; None inherits the current call."""
    mode = _validate_evaluation_mode(mode)
    token = _current_evaluation_mode.set(_current_evaluation_mode.get() if mode is None else mode)
    try:
        yield
    finally:
        _current_evaluation_mode.reset(token)

__all__ = ['AnalyticExpressionStore', 'AnalyticExpressionIR', 'CompiledAnalyticExpression', 'CompiledAnalyticEvaluation', 'AnalyticCorrectionCompiler']


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
    # This selects symbolic coefficient evaluation in the controller, rather
    # than a numerical worker option (workers contract the projected kernels).
    compilation_options = dict(compilation_options or {})
    artifact_options = {name: compilation_options.pop(name) for name in
                        ('ir_file', 'compiled_file', 'storage_backend', 'reuse_only')
                        if name in compilation_options}
    evaluation_options = dict(evaluation_options or {})
    mode = _validate_evaluation_mode(evaluation_options.pop('evaluation_mode', None))
    mode_token = _current_evaluation_mode.set(_current_evaluation_mode.get() if mode is None else mode)
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
    artifacts = (store.artifact_session(**artifact_options,
                    evaluation_folder=evaluation_options.get('evaluation_folder'))
                 if owner and store is not None else contextlib.nullcontext())
    try:
        with artifacts, manager:
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
        _current_evaluation_mode.reset(mode_token)
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
                  '_op_maps', '_conditional_stacks', '_cache', '_grouped_terms', '_change_domain_split')


def _freeze_solver(value, kind, rec):
    names = ['hamiltonian_expansion', 'allowed_terms', 'allowed_coefficients',
             'disallowed_coefficients', 'allowed_energy_changes',
             'intermediate_normalization', 'polynomial_representation', 'resonance_condition',
             'filtered_block_recurrence', 'resonance_independent']
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


@dataclasses.dataclass(frozen=True)
class _OccupationCoefficientProjection:
    powers: tuple
    shift: float = -.5
    fixed_powers: tuple = ()
    preserve_domains: bool = False
    kind: str = dataclasses.field(default='coefficient', init=False)

    def compiler(self):
        from .Analytic import PolyPath, PTEnergyChangeProductSum, canonical_condition_key
        projection = self
        class Transform:
            def __init__(self):
                self.memo = {};self.constants = {}
            def coefficient(self, value, rank):
                powers = projection.fixed_powers + projection.powers
                if rank < len(powers):
                    return 0.
                if np.isscalar(value):
                    return value if not any(powers) else 0.
                key = value,rank
                if key not in self.memo:
                    if len(self.memo) >= 10000:
                        self.memo.pop(next(iter(self.memo)))
                    self.memo[key] = PolyPath.from_polynomial(value).coefficient(
                        powers, shift=projection.shift)
                return self.memo[key]
            def constant(self, value, rank):
                key = float(value),rank
                if key not in self.constants:
                    if len(self.constants) >= 10000:
                        self.constants.pop(next(iter(self.constants)))
                    self.constants[key] = PolyPath.from_coeffs([np.ones(1)]*rank,
                                                             prefactor=key[0], steps=0)
                return self.constants[key]
            def __call__(self, key, value):
                rank = len({i for c in key for i in c[2:]})
                if rank < len(projection.fixed_powers) + len(projection.powers):
                    return None
                if isinstance(value,PTEnergyChangeProductSum):
                    rows = {}
                    for energy,numerator in value.terms.items():
                        domains, harmonic = [], []
                        for row in energy:
                            if getattr(row[0],'_filtered_change_domain',False) and projection.preserve_domains:
                                domains.append(row)
                            elif getattr(row[0],'_state_projection_factor',False) or row[0] not in {0,1}:
                                raise ValueError('coefficient projection requires harmonic divisors and native change domains')
                            else:
                                harmonic.append(row)
                        canonical,scale = canonical_condition_key(tuple(harmonic))
                        if canonical is None:
                            raise ValueError('invalid harmonic denominator')
                        rows.setdefault(canonical+tuple(domains),[]).append(scale*self.coefficient(numerator,rank))
                    rows = {energy:self.constant(weight,rank) for energy,values in rows.items()
                            if (weight := math.fsum(values)) != 0}
                    return value.mutate(rows,reduced=True) if rows else None
                weight = self.coefficient(value,rank)
                return self.constant(weight,rank) if weight else None
        return Transform()


@dataclasses.dataclass(frozen=True)
class _AnchoredModeSubset:
    anchors: tuple
    alternatives: tuple | None = None

    def __post_init__(self):
        if self.alternatives is None:return
        anchors=tuple(self.anchors)
        supports=sorted({frozenset(s)-set(anchors) for s in self.alternatives},key=lambda s:(len(s),tuple(sorted(s))))
        minimal=[]
        for support in supports:
            if not any(prior<=support for prior in minimal):minimal.append(support)
        common=set.intersection(*(set(s) for s in minimal)) if minimal else set()
        object.__setattr__(self,'anchors',anchors+tuple(sorted(common)))
        object.__setattr__(self,'alternatives',tuple(tuple(sorted(s-common)) for s in minimal))

    @staticmethod
    @functools.lru_cache(maxsize=32768)
    def _count(available, supports, k):
        n=available.bit_count()
        if k<0 or k>n or not supports:return 0
        if 0 in supports:return math.comb(n,k)
        supports=tuple(s for s in supports if s.bit_count()<=k and s&available==s)
        if not supports:return 0
        bit=available&-available
        rest=available^bit
        included=tuple(sorted({s&~bit for s in supports}))
        excluded=tuple(s for s in supports if not s&bit)
        return _AnchoredModeSubset._count(rest,included,k-1)+_AnchoredModeSubset._count(rest,excluded,k)

    def _inputs(self,ndim,num_fixed,free_count):
        available=sum(1<<i for i in range(num_fixed,ndim) if i not in self.anchors)
        supports=(0,) if self.alternatives is None else tuple(sum(1<<i for i in s) for s in self.alternatives)
        return available,supports,free_count-len(self.anchors)

    def count(self,ndim,num_fixed,free_count):
        return self._count(*self._inputs(ndim,num_fixed,free_count))

    def iter_subsets(self,ndim,num_fixed,free_count):
        available,supports,k=self._inputs(ndim,num_fixed,free_count)
        if supports==(0,):
            modes=(i for i in range(num_fixed,ndim) if i not in self.anchors)
            if 0<=k<=available.bit_count():
                for rest in itertools.combinations(modes,k):yield self.anchors+rest
            return
        def visit(available,supports,k,prefix):
            if not self._count(available,supports,k):return
            if k==0:
                yield self.anchors+prefix;return
            bit=available&-available;mode=bit.bit_length()-1;rest=available^bit
            included=tuple(sorted({s&~bit for s in supports}))
            excluded=tuple(s for s in supports if not s&bit)
            yield from visit(rest,included,k-1,prefix+(mode,))
            yield from visit(rest,excluded,k,prefix)
        yield from visit(available,supports,k,())



class _CompiledEntries(EntryStore):
    # compiled expressions live at folder/<signature>/manifest.json, aliases at folder/aliases/<signature>.json
    meta_name = 'manifest.json'
    shard = 0
    alias_field = 'compiled'
    stage_prefix = '.compile-'


class _AnalyticIRCodec:
    """Versioned semantic records, deliberately independent of Python pickling.

    References name nodes in a per-part graph; numerical buffers are ordinary
    non-object .npy files. Only explicitly supported mathematical types are
    accepted. Caches, compiled annotations and arbitrary object dictionaries
    never enter the record. Graph traversal is iterative, including decoding.
    """
    version = 2
    record_format = 'psience-analytic-expression'
    path_operations = frozenset({'add', 'linear', 'scale', 'shift', 'remap',
                                 'permutation_sum', 'project', 'mul'})

    def __init__(self, directory, prefix):
        self.directory = Path(directory)
        self.prefix = prefix
        self.nodes = []
        self.pending = []
        self.identities = {}
        self.arrays = []
        self.array_blocks = {}
        self.array_values = {}

    def _value(self, value):
        if value is None or type(value) in (str, bool):
            return value
        if type(value) is int:
            # Ordinary JSON readers preserve these integers exactly. Keep an
            # explicit decimal tag beyond the binary64-safe integer range.
            return value if abs(value) <= 2**53-1 else {'integer': str(value)}
        if type(value) is float:
            return {'float': value.hex()}
        if type(value) is complex:
            return {'complex': [value.real.hex(), value.imag.hex()]}
        if type(value) is bytes:
            return {'bytes': value.hex()}
        if isinstance(value, np.generic):
            dtype = value.dtype
            if (dtype.kind not in 'biufcSU' or dtype.kind == 'f' and dtype.itemsize > 8
                    or dtype.kind == 'c' and dtype.itemsize > 16):
                raise TypeError('unsupported portable IR scalar precision/type')
            return {'scalar': [dtype.newbyteorder('<').str, self._value(value.item())]}
        if type(value) in (tuple, list):
            items = [self._value(v) for v in value]
            return items if type(value) is tuple else {'list': items}
        if type(value) is dict:
            return {'mapping': [[self._value(k), self._value(v)] for k, v in value.items()]}
        identity = id(value)
        if identity not in self.identities:
            if isinstance(value, np.ndarray) and not value.dtype.hasobject and value.dtype.fields is None:
                array = value.astype(value.dtype.newbyteorder('<'), copy=False)
                content = (array.dtype.str, array.shape, hashlib.sha256(array.tobytes(order='C')).digest())
                previous = self.array_values.get(content)
                if previous is not None:
                    # Published numerical arrays are immutable. Equal factors
                    # can share one mapped node even when the compiler created
                    # separate ndarray objects for different numeric programs.
                    return {'ref': previous}
                self.array_values[content] = len(self.nodes)
            self.identities[identity] = len(self.nodes)
            self.nodes.append(None)
            self.pending.append(value)
        return {'ref': self.identities[identity]}

    def _record(self, value):
        from . import Analytic as a
        from McUtils.Zachary.Symbolic import ElementaryFunctions as ef
        kind = type(value)
        if isinstance(value, (np.ndarray, np.generic)):
            array = np.asarray(value)
            if array.dtype.hasobject or array.dtype.fields is not None or array.dtype.kind not in 'biufcSU':
                raise TypeError('IR requires primitive, non-object numerical/string arrays')
            if (array.dtype.kind == 'f' and array.dtype.itemsize > 8
                    or array.dtype.kind == 'c' and array.dtype.itemsize > 16):
                raise TypeError('portable IR extended-precision encoding is not implemented')
            # Explicit byte order, never a native-size alias such as intp.
            array = array.astype(array.dtype.newbyteorder('<'), copy=False)
            dtype = array.dtype.str
            if dtype not in self.array_blocks:
                name = '{}-buffer-{:03d}.npy'.format(self.prefix, len(self.array_blocks))
                self.array_blocks[dtype] = dict(file=name, size=0, arrays=[])
                self.arrays.append(name)
            block = self.array_blocks[dtype]
            offset = block['size']
            block['arrays'].append(array)
            block['size'] += array.size
            return dict(kind='array', file=block['file'], dtype=dtype,
                        shape=list(array.shape), offset=offset, count=array.size)
        if isinstance(value, ef.Functionlike):
            if kind is ef.Variable: function, fields = 'variable', (value.name, value.idx)
            elif kind is ef.Scalar: function, fields = 'scalar', (value.scalar, value.idx)
            elif kind is ef.Power: function, fields = 'power', (value.power, value.idx)
            elif kind in {ef.Exponent, ef.Logarithm}:
                function, fields = ('exponent' if kind is ef.Exponent else 'logarithm'), (value.base, value.idx)
            elif kind in {ef.Identity, ef.Exp, ef.Ln, ef.Sin, ef.Cos}:
                function, fields = {ef.Identity:'identity', ef.Exp:'exp', ef.Ln:'ln', ef.Sin:'sin', ef.Cos:'cos'}[kind], (value.idx,)
            elif kind in {ef.ElementarySummation, ef.ElementaryProduct, ef.ElementaryComposition,
                          ef.Summation, ef.Product, ef.Composition}:
                function = {ef.ElementarySummation:'elementary-sum', ef.ElementaryProduct:'elementary-product',
                            ef.ElementaryComposition:'elementary-composition', ef.Summation:'sum',
                            ef.Product:'product', ef.Composition:'composition'}[kind]
                fields = value.functions, (value.idx if isinstance(value, ef.ElementaryFunction) else value._inds)
            else: raise TypeError('unsupported portable elementary function: ' + kind.__name__)
            return dict(kind='elementary-function', fields=self._value((function, fields)))
        if kind is a.PolyAtom:
            fields = (value._base_coeffs, value.shift)
            tag = 'polynomial-factor'
        elif kind is a.PolyAxis:
            tag, fields = 'polynomial-axis', value.atoms
        elif kind is a.PolyTerm:
            tag, fields = 'polynomial-term', (value.axes, value.steps)
        elif kind is a.PolyPath:
            if value._node is not None and value._node[0] not in self.path_operations:
                raise TypeError('unsupported IR polynomial operation: ' + value._node[0])
            tag = 'polynomial-path'
            fields = (value._items, value._node, value.ndim, value.order, value.reduced)
        elif kind is a.ProductPTPolynomial:
            tag, fields = 'product-polynomial', (value.coeffs, value.prefactor, value.steps)
        elif kind is a.ProductPTPolynomialSum:
            tag, fields = 'polynomial-sum', (value.polys, value.prefactor, value.reduced)
        elif kind in (a.PTEnergyChangeProductSum, a.PTTensorCoeffProductSum):
            tag = 'energy-sum' if kind is a.PTEnergyChangeProductSum else 'tensor-sum'
            fields = (value.terms, value.prefactor, value.reduced, value.ndim)
        elif kind is a.SqrtChangePoly:
            tag, fields = 'sqrt-change', (value.poly_obj, value.poly_change, value.shift_start)
        elif kind is a.ProjectionFactor:
            tag, fields = 'projection', (value.kind, value.slot)
        elif kind is a._ChangeDomainIndex:
            tag, fields = 'change-domain', (value.retain, value.endpoint)
        elif kind is a.ConditionEndpoint:
            tag, fields = 'condition-endpoint', (value.condition, value.kind, value.retain, value.endpoint)
        elif kind is a.SameStateCondition:
            tag, fields = 'same-state', ()
        elif kind is a.ResonanceBlockCondition:
            tag, fields = 'resonance-blocks', ()
        elif kind is a.BooleanCondition:
            tag, fields = 'boolean-condition', (value.operation, value.conditions)
        elif kind is a.SymmetricCondition:
            tag, fields = 'symmetric-condition', (value.condition,)
        elif kind is a.EnergyWindowCondition:
            tag, fields = 'energy-window', (value.tolerance, value.modes)
        elif kind is a.IntegerLatticeCondition:
            tag = 'integer-lattice'
            fields = (value.modes, value.basis)
        elif kind is a.PolynomialInvariantCondition:
            tag, fields = 'polynomial-invariant', (value.functions, value.modes, value.freeze_other_modes, value.tolerance)
        elif kind is a.PolynomialRelationshipCondition:
            tag, fields = 'polynomial-relationship', (value.function, value.modes, value.comparison,
                                                    value.tolerance, value.freeze_other_modes)
        else:
            raise TypeError('unsupported analytic IR value: ' + kind.__name__)
        return dict(kind=tag, fields=self._value(fields))

    @staticmethod
    def storage_backend(backend):
        if type(backend) is not str or backend not in {'portable', 'pickle'}:
            raise ValueError("storage_backend must be 'portable' or 'pickle'")
        return backend

    @classmethod
    def record_suffix(cls, backend):
        return '.json' if cls.storage_backend(backend) == 'portable' else '.pkl'

    @classmethod
    def storage_version(cls, backend):
        return 2 if cls.storage_backend(backend) == 'pickle' else 1

    @staticmethod
    def content_hash(record):
        """Identify the same mathematical/program data across storage backends."""
        digest = hashlib.sha256()
        encoder = json.JSONEncoder(sort_keys=True, separators=(',', ':'), allow_nan=False)
        for chunk in encoder.iterencode(record): digest.update(chunk.encode('utf-8'))
        return digest.hexdigest()

    @classmethod
    def write_record(cls, path, record, backend):
        if cls.storage_backend(backend) == 'portable':
            atomic_write_json(path, record, indent=None, sort_keys=True)
        else:
            _atomic_pickle(path, record)

    def write(self, expression, *, storage_backend='portable'):
        suffix = self.record_suffix(storage_backend)
        root = self._value(expression)
        position = 0
        while position < len(self.pending):
            self.nodes[position] = self._record(self.pending[position])
            position += 1
        hashes = {}
        for dtype, block in self.array_blocks.items():
            destination = self.directory / block['file']
            with destination.open('wb') as stream:
                np.lib.format.write_array_header_2_0(stream, dict(descr=dtype,
                    fortran_order=False, shape=(block['size'],)))
                for array in block['arrays']:
                    array.ravel(order='C').tofile(stream)
            hashes[block['file']] = _file_hash(destination)
        for node in self.nodes:
            if node['kind'] == 'array':
                node['sha256'] = hashes[node['file']]
        record = dict(format=self.record_format, version=self.version,
                      root=root, nodes=self.nodes)
        self.content_sha256 = self.content_hash(record)
        if storage_backend == 'portable':
            self.write_record(self.directory / (self.prefix + suffix), record, storage_backend)
            return [self.prefix + suffix] + self.arrays
        envelope = dict(format=self.record_format, version=self.version,
                        representation='python', payload=expression)
        try:
            atomic_write(self.directory/(self.prefix+suffix), lambda stream:
                _AnalyticPickler(stream, type(self)).dump(envelope), fsync=False, keep_mode=False)
        except RecursionError:
            # Python's object pickler is recursive. Keep exceptionally deep
            # expressions in the same validated, flat reference graph instead.
            self.write_record(self.directory/(self.prefix+suffix),
                dict(format=self.record_format, version=self.version,
                     representation='graph', payload=record), storage_backend)
            return [self.prefix+suffix]+self.arrays
        # These pools were used to compute the backend-neutral content digest.
        # Pickle stores its immutable numerical data inline in the part instead.
        for name in self.arrays: (self.directory/name).unlink()
        return [self.prefix+suffix]

    @classmethod
    def read(cls, record, array_loader, *, consume=False):
        from . import Analytic as a
        if record.get('format') != cls.record_format or record.get('version') not in (1, cls.version):
            raise ValueError('unsupported analytic IR record')
        records = record['nodes']
        values, active = {}, set()

        def references(value):
            todo = [value]
            while todo:
                item = todo.pop()
                if isinstance(item, dict):
                    if set(item) == {'ref'}:
                        index = item['ref']
                        if type(index) is not int or not 0 <= index < len(records):
                            raise ValueError('invalid analytic IR reference')
                        yield index
                    else:
                        todo.extend(item.values())
                elif isinstance(item, list):
                    todo.extend(item)

        def decode(value):
            if value is None or type(value) in (str, bool):
                return value
            if type(value) is int and abs(value) <= 2**53-1:
                return value
            if type(value) is list:
                return tuple(decode(v) for v in value)
            if not isinstance(value, dict) or len(value) != 1:
                raise ValueError('invalid analytic IR value')
            tag, content = next(iter(value.items()))
            if tag == 'ref': return values[content]
            if tag == 'integer': return int(content)
            if tag == 'float': return float.fromhex(content)
            if tag == 'complex': return complex(*(float.fromhex(v) for v in content))
            if tag == 'bytes': return bytes.fromhex(content)
            if tag == 'scalar':
                dtype = np.dtype(content[0])
                if dtype.hasobject or dtype.fields is not None or dtype.kind not in 'biufcSU':
                    raise ValueError('invalid analytic IR scalar type')
                return np.asarray(decode(content[1]), dtype=dtype)[()]
            if tag == 'tuple': return tuple(decode(v) for v in content)
            if tag == 'list': return [decode(v) for v in content]
            if tag == 'mapping': return {decode(k): decode(v) for k, v in content}
            raise ValueError('unknown analytic IR value tag: ' + tag)

        def build(node):
            kind = node['kind']
            if kind == 'array':
                buffer = array_loader(node)
                offset, count = node['offset'], node['count']
                shape = node['shape']
                if (buffer.dtype.str != node['dtype'] or buffer.ndim != 1
                        or type(offset) is not int or type(count) is not int
                        or offset < 0 or count < 0 or offset + count > len(buffer)
                        or any(type(d) is not int or d < 0 for d in shape)
                        or math.prod(shape) != count):
                    raise ValueError('analytic IR array descriptor mismatch')
                return buffer[offset:offset + count].reshape(shape)
            fields = decode(node['fields'])
            return cls.restore_node(kind, fields)

        stack = [(index, False) for index in references(record['root'])]
        while stack:
            index, ready = stack.pop()
            if index in values:
                continue
            if ready:
                values[index] = build(records[index])
                if consume: records[index] = None
                active.remove(index)
            else:
                if index in active:
                    raise ValueError('cyclic analytic IR graph')
                active.add(index)
                stack.append((index, True))
                stack.extend((child, False) for child in references(records[index].get('fields'))
                             if child not in values)
        return decode(record['root'])

    @classmethod
    def restore_node(cls, kind, fields):
        from . import Analytic as a
        custom = cls._restore_custom(kind, fields)
        if custom is not NotImplemented:
            return custom
        if kind == 'elementary-function':
            from McUtils.Zachary.Symbolic import ElementaryFunctions as ef
            name, arguments = fields
            simple = {'variable':ef.Variable, 'scalar':ef.Scalar, 'power':ef.Power,
                      'exponent':ef.Exponent, 'logarithm':ef.Logarithm, 'identity':ef.Identity,
                      'exp':ef.Exp, 'ln':ef.Ln, 'sin':ef.Sin, 'cos':ef.Cos}
            if name == 'variable': return ef.Variable(*arguments)
            if name in simple:
                return simple[name](*arguments[:-1], idx=arguments[-1])
            univariate = {'elementary-sum':ef.ElementarySummation, 'elementary-product':ef.ElementaryProduct,
                          'elementary-composition':ef.ElementaryComposition}
            multivariate = {'sum':ef.Summation, 'product':ef.Product, 'composition':ef.Composition}
            if name in univariate: return univariate[name](*arguments[0], idx=arguments[1])
            if name in multivariate: return multivariate[name](*arguments[0], indices=arguments[1])
            raise ValueError('unknown portable function operation: ' + name)
        constructors = {
            'polynomial-factor': lambda: a.PolyAtom(*fields),
            'polynomial-axis': lambda: a.PolyAxis(fields),
            'polynomial-term': lambda: a.PolyTerm(*fields),
            'product-polynomial': lambda: a.ProductPTPolynomial(fields[0], prefactor=fields[1], steps=fields[2]),
            'polynomial-sum': lambda: a.ProductPTPolynomialSum(fields[0], prefactor=fields[1], reduced=fields[2]),
            'sqrt-change': lambda: a.SqrtChangePoly(*fields),
            'projection': lambda: a.ProjectionFactor(*fields),
            'change-domain': lambda: a._ChangeDomainIndex(*fields),
            'condition-endpoint': lambda: a.ConditionEndpoint(*fields),
            'same-state': a.SameStateCondition,
            'resonance-blocks': a.ResonanceBlockCondition,
            'boolean-condition': lambda: a.BooleanCondition(fields[0], *fields[1]),
            'symmetric-condition': lambda: a.SymmetricCondition(*fields),
            'energy-window': lambda: a.EnergyWindowCondition(*fields)
            ,'polynomial-invariant': lambda: a.PolynomialInvariantCondition(*fields)
            ,'polynomial-relationship': lambda: a.PolynomialRelationshipCondition(*fields)
        }
        if kind == 'polynomial-path':
            items, operation, ndim, order, reduced = fields
            if operation is not None and operation[0] not in cls.path_operations:
                raise ValueError('unknown analytic IR polynomial operation')
            return a.PolyPath(items, reduced=reduced, node=operation, ndim=ndim, order=order)
        if kind in ('energy-sum', 'tensor-sum'):
            terms, prefactor, reduced, ndim = fields
            constructor = a.PTEnergyChangeProductSum if kind == 'energy-sum' else a.PTTensorCoeffProductSum
            result = constructor(terms, prefactor=prefactor, canonicalize=False, reduced=reduced)
            result._ndim = ndim
            return result
        if kind == 'integer-lattice':
            # Only the mathematical lattice is in IR. Its evaluation program
            # will belong to the compiled artifact, not this semantic record.
            modes, basis = fields
            relations = np.zeros((len(basis[0]), max(modes) + 1), dtype=int)
            relations[:, modes] = np.asarray(basis, dtype=int).T
            return a.IntegerLatticeCondition(relations)
        if kind not in constructors:
            raise ValueError('unknown analytic IR node kind: ' + kind)
        return constructors[kind]()


    @classmethod
    def _restore_custom(cls, kind, fields):
        return NotImplemented


class _CompiledEnergyProgram(dict):
    """A shared denominator program, retained as a graph node on disk."""


class _CompiledAnalyticCodec(_AnalyticIRCodec):
    record_format = 'psience-compiled-expression'

    def _record(self, value):
        from . import Analytic as a
        from McUtils.Zachary.Symbolic import ElementaryFunctions as ef
        kind = type(value)
        if kind is a._PolynomialProgramTable:
            tag, fields = 'polynomial-table', tuple(getattr(value, name) for name in value.fields)
        elif kind is _CompiledEnergyProgram:
            tag, fields = 'energy-program', tuple(value.items())
        elif kind is ef._FunctionProgram:
            tag, fields = 'function-program', (value.instructions, value.signature, value.elementary, value.idx)
        elif isinstance(value, a.OperatorCondition):
            tag = 'condition-program'
            name, fields = _condition_program_fields(value)
            fields = (name, fields, tuple(getattr(value, field) for field in
                       ('symmetric', 'diagonal', 'equivalence', 'requires_frequencies')), value.key)
        elif kind is a.ConditionEndpoint:
            tag, fields = 'compiled-endpoint', (value.condition, value.kind, value.retain, value.endpoint, int(value))
        elif kind is a._CompiledPolynomial and value.table is not None:
            tag, fields = 'table-polynomial', (value.table, value.node_index, value.ndim, value.order,
                                              value.signature, value.zero, value.occupation_shifts)
        elif kind is a._CompiledPolynomial:
            tag = 'polynomial-program'
            fields = (value.instructions, value.ndim, value.order, value.signature, value.zero,
                      value.occupation_shifts, value.graph_program)
        elif kind is a.PTEnergyChangeProductSum and hasattr(value, '_compiled_terms'):
            tag, fields = 'compiled-energy-plan', (value._compiled_terms, value.prefactor, value.reduced, value.ndim)
        elif kind in (a.PTTensorCoeffProductSum, a.PTEnergyChangeProductSum):
            tag = 'compiled-tensor' if kind is a.PTTensorCoeffProductSum else 'compiled-energy'
            fields = (value.terms, value.prefactor, value.reduced, value.ndim,
                      value._inds_map if kind is a.PTTensorCoeffProductSum else None)
        elif kind is IndexedSubexpressionPlan:
            tag = 'subexpression-plan'
            fields = (value.expression, value.direct_polynomial, value.prefactor, value.groups)
        elif kind is IndexedEnergyTerm:
            tag = 'energy-term-plan'
            fields = tuple(getattr(value, name) for name in
                           ('ordinal', 'energy_changes', 'polynomial', 'contraction_shape',
                            'projected', 'constant_numerator', 'denominator_program'))
        else:
            return super()._record(value)
        return dict(kind=tag, fields=self._value(fields))

    @classmethod
    def _restore_custom(cls, kind, fields):
        from . import Analytic as a
        if kind == 'polynomial-table': return a._PolynomialProgramTable(*fields)
        if kind == 'table-polynomial':
            table, index, ndim, order, signature, zero, occupations = fields
            return a._CompiledPolynomial((), ndim, order, signature, zero, occupations, table=table, node_index=index)
        if kind == 'energy-program': return _CompiledEnergyProgram(fields)
        if kind == 'compiled-energy-plan':
            terms, prefactor, reduced, ndim = fields
            result = object.__new__(a.PTEnergyChangeProductSum)
            result.terms = {term.energy_changes: term.polynomial for term in terms}
            result.prefactor, result.reduced, result._ndim = prefactor, reduced, ndim
            result._shape, result._compiled_terms = None, terms
            return result
        if kind == 'function-program':
            from McUtils.Zachary.Symbolic.ElementaryFunctions import _FunctionProgram
            return _FunctionProgram(*fields)
        if kind == 'condition-program':
            name, arguments, flags, key = fields
            result = _restore_condition_program(name, arguments)
            for field, value in zip(('symmetric', 'diagonal', 'equivalence', 'requires_frequencies'), flags):
                object.__setattr__(result, field, value)
            object.__setattr__(result, 'key', key)
            return result
        if kind == 'compiled-endpoint':
            condition, role, retain, endpoint, number = fields
            result = int.__new__(a.ConditionEndpoint, number)
            for name, value in (('condition', condition), ('kind', role), ('retain', retain),
                                ('endpoint', endpoint), ('_frozen', True)):
                object.__setattr__(result, name, value)
            return result
        if kind == 'polynomial-program':
            return a._CompiledPolynomial(*fields)
        if kind in {'compiled-tensor', 'compiled-energy'}:
            constructor = a.PTTensorCoeffProductSum if kind == 'compiled-tensor' else a.PTEnergyChangeProductSum
            result = object.__new__(constructor)
            result.terms, result.prefactor, result.reduced, result._ndim, indices = fields
            result._shape = None
            if indices is not None: result._inds_map = indices
            return result
        if kind == 'subexpression-plan': return IndexedSubexpressionPlan(*fields)
        if kind == 'energy-term-plan': return IndexedEnergyTerm(*fields)
        return NotImplemented


def _restore_analytic_array(dtype, shape, data, scalar):
    array = np.frombuffer(data, dtype=np.dtype(dtype)).reshape(shape)
    return array[()] if scalar else array


def _restore_analytic_pickle_node(compiled, kind, fields):
    codec = _CompiledAnalyticCodec if compiled else _AnalyticIRCodec
    return codec.restore_node(kind, fields)


class _AnalyticPickler(pickle.Pickler):
    """Use the same whitelisted node fields while preserving Python memoization."""
    def __init__(self, stream, codec):
        super().__init__(stream, protocol=5)
        self.codec = codec(Path('.'), '')
        self.compiled = codec.record_format == 'psience-compiled-expression'
        # _record describes exactly the same fields for both storage modes.
        # Native pickle traverses them itself instead of making a tagged graph.
        self.codec._value = lambda value: value

    def reducer_override(self, value):
        if isinstance(value, (np.ndarray, np.generic)):
            array = np.asarray(value)
            dtype = array.dtype
            if (dtype.hasobject or dtype.fields is not None or dtype.kind not in 'biufcSU'
                    or dtype.kind == 'f' and dtype.itemsize > 8
                    or dtype.kind == 'c' and dtype.itemsize > 16):
                raise TypeError('unsupported analytic pickle numerical precision/type')
            array = array.astype(dtype.newbyteorder('<'), copy=False)
            return _restore_analytic_array, (array.dtype.str, array.shape, array.tobytes(order='C'),
                                             isinstance(value, np.generic))
        if isinstance(value, (type, types.FunctionType, types.BuiltinFunctionType)):
            return NotImplemented
        record = self.codec._record(value)
        return _restore_analytic_pickle_node, (self.compiled, record['kind'], record['fields'])


def _condition_program_layouts():
    from . import Analytic as a
    return {
        'same-state': (a.SameStateCondition, ()), 'resonance-blocks': (a.ResonanceBlockCondition, ()),
        'boolean': (a.BooleanCondition, ('operation', 'conditions')),
        'symmetric': (a.SymmetricCondition, ('condition',)),
        'energy-window': (a.EnergyWindowCondition, ('tolerance', 'modes')),
        'integer-lattice': (a.IntegerLatticeCondition, ('modes', 'basis', 'pivots', 'inverse_numerators', 'denominator')),
        'invariant': (a.PolynomialInvariantCondition, ('functions', 'modes', 'freeze_other_modes', 'tolerance')),
        'relationship': (a.PolynomialRelationshipCondition, ('function', 'modes', 'comparison', 'tolerance', 'freeze_other_modes'))
    }


def _condition_program_fields(condition):
    for name, (constructor, attributes) in _condition_program_layouts().items():
        if type(condition) is constructor:
            return name, tuple(getattr(condition, attribute) for attribute in attributes)
    raise TypeError('unsupported compiled condition: ' + type(condition).__name__)


def _restore_condition_program(name, fields):
    layouts = _condition_program_layouts()
    if name not in layouts: raise ValueError('unknown compiled condition program: ' + name)
    constructor, attributes = layouts[name]
    if len(fields) != len(attributes): raise ValueError('invalid compiled condition fields')
    result = object.__new__(constructor)
    for attribute, value in zip(attributes, fields): object.__setattr__(result, attribute, value)
    return result


class AnalyticExpressionIR:
    """A lazy, portable simplified expression; no numerical system is bound."""
    format_version = 1
    artifact_kind = 'analytic-expression-ir'
    version_field = 'ir_version'
    codec = _AnalyticIRCodec

    def __init__(self, source, manifest, archive=None, prefix='', expression=None):
        self.source = Path(source)
        self.manifest = manifest
        self.archive = archive
        self.prefix, self.expression = prefix, expression
        self.storage_backend = self.codec.storage_backend(manifest.get('storage_backend', 'portable'))

    @classmethod
    def load(cls, source, *, expression=None):
        from McUtils.Scaffolding import IndexedTarStore
        source = Path(source).expanduser().resolve()
        archive = None if source.is_dir() else IndexedTarStore(source, read_only=True)
        try:
            catalog = cls._catalog(source, archive)
            if catalog is not None:
                if expression is None:
                    if len(catalog['entries']) != 1:
                        raise ValueError('select an expression from the archive: ' + ', '.join(catalog['entries']))
                    expression = next(iter(catalog['entries']))
                if expression not in catalog['entries']: raise KeyError(expression)
                prefix = catalog['entries'][expression]
                if IndexedTarStore.normalize_name(prefix) != prefix:
                    raise ValueError('invalid expression catalog prefix')
            else:
                if expression is not None: raise KeyError(expression)
                prefix = '' if archive is None else 'expression'
            member = prefix+'/manifest.json' if prefix else 'manifest.json'
            data = (source/member).read_bytes() if archive is None else archive.read(member)
            manifest = json.loads(data)
            if (manifest.get('kind') != cls.artifact_kind
                    or manifest.get(cls.version_field) != cls.format_version or not manifest.get('complete')):
                raise ValueError('incomplete or unsupported analytic expression IR')
            storage_version = manifest.get('storage_version', 1)
            if (storage_version not in {1, 2} or
                    storage_version == 2 and manifest.get('storage_backend') != 'pickle'):
                raise ValueError('unsupported analytic storage version')
            return cls(source, manifest, archive, prefix, expression)
        except BaseException:
            if archive is not None:
                archive.close()
            raise

    @classmethod
    def _catalog(cls, source, archive=None):
        if archive is None:
            path = Path(source)/'manifest.json'
            data = path.read_bytes() if path.is_file() else None
        else:
            try: data = archive.read('manifest.json')
            except FileNotFoundError: data = None
        if data is None: return None
        manifest = json.loads(data)
        if manifest.get('kind') == cls.artifact_kind+'-collection':
            if (manifest.get('catalog_version') != 1 or not manifest.get('complete')
                    or not isinstance(manifest.get('entries'), dict)):
                raise ValueError('unsupported or incomplete analytic expression collection')
            return manifest
        return None

    @classmethod
    def expression_names(cls, source):
        from McUtils.Scaffolding import IndexedTarStore
        source = Path(source).expanduser().resolve()
        manager = contextlib.nullcontext(None) if source.is_dir() else IndexedTarStore(source, read_only=True)
        with manager as archive:
            catalog = cls._catalog(source, archive)
            return () if catalog is None else tuple(catalog['entries'])

    @staticmethod
    def _name(name):
        if type(name) is not str or not name or Path(name).name != name or name in {'.', '..'}:
            raise ValueError('invalid analytic IR member path')
        return name

    def _read(self, name):
        with self._open(name) as stream:
            return stream.read()

    def _member(self, name):
        name = self._name(name)
        return self.prefix+'/'+name if self.prefix else name

    def _open(self, name):
        member = self._member(name)
        return (self.source/member).open('rb') if self.archive is None else self.archive.open(member)

    def iter_parts(self):
        for part in self.manifest['parts']:
            yield self.read_part(part)

    def _read_record(self, part):
        if Path(part['file']).suffix != self.codec.record_suffix(self.storage_backend):
            raise ValueError('analytic part does not match its storage backend')
        data = self._read(part['file'])
        if hashlib.sha256(data).hexdigest() != part['sha256']:
            raise ValueError('analytic IR part checksum mismatch')
        record = json.loads(data) if self.storage_backend == 'portable' else pickle.loads(data)
        if (not isinstance(record, dict) or record.get('format') != self.codec.record_format
                or record.get('version') not in (1, self.codec.version)):
            raise ValueError('unsupported analytic expression record')
        return record

    def read_part(self, part):
        record = self._read_record(part)
        return self._decode_record(record)

    def _decode_record(self, record):
        if record.get('representation') == 'python':
            if self.storage_backend != 'pickle' or self.manifest.get('storage_version', 1) != 2:
                raise ValueError('native analytic data requires pickle storage version 2')
            return record['payload']
        if record.get('representation') == 'graph':
            if self.storage_backend != 'pickle' or self.manifest.get('storage_version', 1) != 2:
                raise ValueError('flattened analytic data requires pickle storage version 2')
            record = record['payload']
        buffers = {}
        def array_loader(descriptor):
            name = self._name(descriptor['file'])
            if name in buffers:
                expected, buffer = buffers[name]
                if expected != descriptor['sha256']:
                    raise ValueError('inconsistent analytic IR buffer checksum')
                return buffer
            if self.archive is None:
                path = self.source/self._member(name)
                if _file_hash(path) != descriptor['sha256']:
                    raise ValueError('analytic IR array checksum mismatch')
                buffer = np.load(path, mmap_mode='r', allow_pickle=False)
            else:
                digest = hashlib.sha256()
                with self._open(name) as stream:
                    for block in iter(lambda: stream.read(1024 * 1024), b''):
                        digest.update(block)
                if digest.hexdigest() != descriptor['sha256']:
                    raise ValueError('analytic IR array checksum mismatch')
                buffer = self.archive.memmap(self._member(name))
            buffers[name] = descriptor['sha256'], buffer
            return buffer
        return self.codec.read(record, array_loader, consume=True)

    def save(self, file, *, storage_backend=None):
        """Repack or change storage without simplifying or compiling expressions."""
        target = AnalyticExpressionStore._pack_readers(
            {self.expression or 'expression': self}, file, storage_backend=storage_backend,
            single=self.expression is None)
        return type(self).load(target)

    def close(self):
        if self.archive is not None:
            self.archive.close()

    def __enter__(self): return self
    def __exit__(self, *exc): self.close()


class _CompiledAnalyticArchive(AnalyticExpressionIR):
    artifact_kind = 'compiled-analytic-expression'
    version_field = 'plan_version'
    codec = _CompiledAnalyticCodec

    @classmethod
    def load(cls, source, *, expression=None):
        result = super().load(source, expression=expression)
        versions = result.manifest.get('program_versions', {})
        if (versions.get('polynomial') not in {1, 2, 3} or versions.get('denominator') not in {1, 2} or
                {key: value for key, value in versions.items() if key not in {'polynomial', 'denominator'}} !=
                {'function':1, 'degeneracy':1}):
            result.close()
            raise ValueError('unsupported compiled analytic program versions')
        return result

    def read_part(self, part):
        payload = super().read_part(part)
        expression = payload['expression']
        execution = payload['execution']
        if execution.get('version') != 1:
            raise ValueError('unsupported compiled expression execution version')
        execution['plans'] = {name: IndexedExpressionPlan(expression,
            types.MappingProxyType(plans), count,
            None if bindings is None else types.MappingProxyType(bindings))
            for name, (plans, count, bindings) in payload['plans'].items()}
        expression._compiled_execution = execution
        if payload['free_index'] is not None:
            expression._free_index_plan_data = payload['free_index']
        return expression


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
        self._artifacts = None
        self.timings = dict(requests=0, compilations=0, cache_hits=0, failures=0,
            compilation_seconds=0., cache_lookup_seconds=0., preparation_seconds=0., evaluation_seconds=0.,
            operator_coefficient_contractions=0,operator_coefficient_cache_hits=0,operator_coefficient_batches=0)

    def statistics(self):
        return dict(self.timings, storage='memory' if self.folder is None else 'disk',
                    memory_entries=len(self._memory), memory_payload_bytes=self.memory_bytes)

    def clear(self):
        """Release reusable in-memory entries; never remove persistent cache files."""
        self._memory.clear(); self._aliases.clear(); self.memory_bytes=0

    @property
    def using_artifacts(self):
        return self._artifacts is not None

    @contextlib.contextmanager
    def artifact_session(self, *, ir_file=None, compiled_file=None, storage_backend=None,
                         reuse_only=False, evaluation_folder=None):
        """Reuse/build whole libraries during one public runner calculation.

        Entries are identified by versioned formal requests, without producer
        source hashes or numerical system inputs. On normal exit, publish each
        library atomically from temporary staging.
        """
        if ir_file is None and compiled_file is None:
            if storage_backend is not None or reuse_only:
                raise ValueError('storage_backend/reuse_only requires compiled_file')
            yield
            return
        if compiled_file is None:
            raise ValueError('ir_file requires a separate compiled_file')
        if self.using_artifacts:
            raise RuntimeError('an expression artifact session is already active')
        if not isinstance(reuse_only, (bool, np.bool_)):
            raise TypeError('reuse_only must be boolean')
        if storage_backend is not None: _AnalyticIRCodec.storage_backend(storage_backend)
        paths = dict(compiled=Path(compiled_file).expanduser().resolve(),
                     ir=None if ir_file is None else Path(ir_file).expanduser().resolve())
        if paths['ir'] == paths['compiled']:
            raise ValueError('compiled_file must be separate from ir_file')
        if evaluation_folder is not None:
            destination = Path(evaluation_folder).expanduser().resolve()
            for path in paths.values():
                if path is not None and (destination == path.parent or path.parent in destination.parents):
                    raise ValueError('evaluation_folder must be outside the expression artifact folders')
        import shutil
        import tempfile
        with tempfile.TemporaryDirectory(prefix='analytic-run-') as folder:
            with contextlib.ExitStack() as readers:
                state = dict(folder=Path(folder), paths=paths, readers=readers,
                             sources=dict(ir={}, compiled={}), dirty=set(), results=[],
                             backend=storage_backend, reuse_only=bool(reuse_only))
                for kind, path in paths.items():
                    if path is None or not path.exists(): continue
                    reader_type = AnalyticExpressionIR if kind == 'ir' else _CompiledAnalyticArchive
                    names = reader_type.expression_names(path)
                    if not names:
                        raise ValueError('runner expression files must be named libraries')
                    for name in names:
                        reader = readers.enter_context(reader_type.load(path, expression=name))
                        state['sources'][kind][name] = reader
                        if storage_backend is not None and reader.storage_backend != storage_backend:
                            state['dirty'].add(kind)
                self._artifacts = state
                candidates = {}
                try:
                    yield
                    for kind in state['dirty']:
                        if paths[kind] is None: continue
                        candidate = state['folder']/(kind+'.tar')
                        self._pack_readers(state['sources'][kind], candidate,
                                           storage_backend=storage_backend)
                        candidates[kind] = candidate
                finally:
                    self._artifacts = None
            # Close original immutable archive handles before replacing a file.
            for kind, candidate in candidates.items():
                target = paths[kind]
                target.parent.mkdir(parents=True, exist_ok=True)
                def publish(stream, candidate=candidate):
                    with candidate.open('rb') as source:
                        shutil.copyfileobj(source, stream, length=1024*1024)
                atomic_write(target, publish, keep_mode=False)
            for name, result in state['results']:
                final = CompiledAnalyticExpression.load(paths['compiled'], expression=name)
                result.folder, result.artifact = final.folder, final.artifact
                result.expression, result.manifest = final.expression, final.manifest

    def _compile_artifact(self, factory, *, context, active_operators=None, projection=None,
                          mode_dimension=None, free_index_reduction=True, progress=None, **options):
        state = self._artifacts
        unknown = set(options) - {'max_terms', 'max_bytes', 'cache_size', 'cache_bytes',
                                  'node_cache_size', 'coalesce_combine_nodes', 'bounded_byte_accounting'}
        if unknown: raise TypeError('unknown artifact compilation options: ' + ', '.join(sorted(unknown)))
        if not isinstance(free_index_reduction, (bool, np.bool_)):
            raise TypeError('free_index_reduction must be boolean')
        request = _fingerprint(dict(request_version=1, definition=context, projection=projection,
                                    active_operators=active_operators, mode_dimension=mode_dimension))
        def source_ir():
            ir = state['sources']['ir'].get(request)
            if ir is None:
                ir = self.export_ir(factory(), staging_folder=state['folder']/request/'ir',
                    context=context, active_operators=active_operators, projection=projection,
                    mode_dimension=mode_dimension, progress=progress,
                    storage_backend=state['backend'] or 'portable', **options)
                state['readers'].enter_context(ir)
                ir.manifest.update(request_version=1, request_signature=request)
                state['sources']['ir'][request] = ir
                state['dirty'].add('ir')
            elif ir.manifest.get('request_version') != 1 or ir.manifest.get('request_signature') != request:
                raise ValueError('unsupported runner expression IR request metadata')
            return ir

        source = state['sources']['compiled'].get(request)
        if source is not None:
            if source.manifest.get('request_version') != 1 or source.manifest.get('request_signature') != request:
                raise ValueError('unsupported runner expression request metadata')
            if source.manifest['options']['free_index_reduction'] == bool(free_index_reduction):
                if state['paths']['ir'] is not None and request not in state['sources']['ir']:
                    if state['reuse_only']:
                        raise ValueError('ir_file lacks the mathematical source for request ' + request)
                    source_ir()
                self.timings['requests'] += 1
                self.timings['cache_hits'] += 1
                result = CompiledAnalyticExpression(None, source.manifest,
                    artifact=source.source, expression=source.expression)
                result.store = self
                state['results'].append((request, result))
                return result
        if state['reuse_only']:
            raise ValueError('compiled_file has no compatible expression for request ' + request)
        ir = source_ir()
        result = self.compile_ir(ir, staging_folder=state['folder']/request/'compiled',
            free_index_reduction=free_index_reduction, storage_backend=state['backend'] or 'portable', progress=progress)
        compiled = state['readers'].enter_context(_CompiledAnalyticArchive.load(result.artifact))
        compiled.manifest.update(request_version=1, request_signature=request)
        state['sources']['compiled'][request] = compiled
        state['dirty'].add('compiled')
        state['results'].append((request, result))
        return result

    def export_ir(self, expression, file=None, *, staging_folder=None, storage_backend='portable', **options):
        """Simplify/stream an expression to portable IR without execution plans.

        A directory may be retained with ``file=None, staging_folder=...``.
        Otherwise only the final self-contained uncompressed tar is published.
        Physical pruning is opt-in through active_operators/mode_dimension and
        is recorded in the manifest. No coefficients or frequencies are used.
        """
        import tempfile
        from McUtils.Scaffolding import IndexedTarStore
        storage_backend = _AnalyticIRCodec.storage_backend(storage_backend)
        if file is None and staging_folder is None:
            raise ValueError('IR export needs a file or a staging folder')
        if isinstance(expression, (str, os.PathLike)):
            # Data migration deliberately discards the old execution annotations.
            # Validate its complete manifest and part checksums, then accept only
            # the mathematical types supported by the explicit IR codec.
            expression = CompiledAnalyticExpression.load(expression, validate_source=False)
        if isinstance(expression, CompiledAnalyticExpression):
            if options: raise TypeError('a simplified cached expression cannot be specialized again during export')
            return self._export_cached_ir(expression, file, staging_folder=staging_folder,
                                          storage_backend=storage_backend)
        if options.pop('free_index_reduction', False):
            raise ValueError('execution free-index plans belong to IR compilation')
        manager = (tempfile.TemporaryDirectory(prefix='analytic-ir-') if staging_folder is None
                   else contextlib.nullcontext(str(Path(staging_folder).expanduser().resolve())))
        with manager as folder:
            store = type(self)(Path(folder) / 'expressions')
            from .Analytic import PolyPath
            token = PolyPath._index_reduction_context.set(False)
            try:
                source_context = options.pop('context', None)
                directory = store._compile_factory(lambda: expression,
                    context=expression if source_context is None else source_context, free_index_reduction=False,
                    part_format='ir', storage_backend=storage_backend, **options)
            finally:
                PolyPath._index_reduction_context.reset(token)
            if file is None:
                return AnalyticExpressionIR.load(directory)
            target = Path(file).expanduser().resolve()
            target.parent.mkdir(parents=True, exist_ok=True)
            archive = IndexedTarStore(Path(folder) / 'build.tar')
            archive.write_trees([('expression', directory, {})])
            # Build locking belongs to the staging store. Publish through the
            # existing atomic writer without leaving a target-side lock file.
            atomic_write(target, archive.save, keep_mode=False)
        return AnalyticExpressionIR.load(target)

    def _export_cached_ir(self, expression, file, *, staging_folder=None, storage_backend='portable'):
        import tempfile
        from McUtils.Scaffolding import IndexedTarStore
        if expression.artifact is not None:
            raise ValueError('retain the mathematical source IR when compiling portable programs')
        if expression.manifest.get('active_operators') is not None and expression.manifest.get('operator_keys') is None:
            raise ValueError('specialized legacy cache lacks its operator domain; re-export from the original expression')
        manager = (tempfile.TemporaryDirectory(prefix='analytic-cache-ir-') if staging_folder is None else
                   contextlib.nullcontext(str(Path(staging_folder).expanduser().resolve())))
        with manager as folder:
            entries = _CompiledEntries(Path(folder)/'expressions')
            key = _fingerprint(('cached-ir', expression.manifest['signature'], storage_backend))
            directory = Path(entries.entry_path(key))
            operators = set()
            with entries.stage(key) as staged:
                stage, parts = Path(staged.path), []
                for original, part in zip(expression.manifest['parts'], expression.iter_parts()):
                    name = 'part-{:06d}-r{}'.format(len(parts), original['rank'])
                    codec = _AnalyticIRCodec(stage, name)
                    files = codec.write(part, storage_backend=storage_backend)
                    operators.update(coefficient[:2] for product in part.terms for coefficient in product)
                    parts.append(dict(original, file=files[0], buffers=files[1:], sha256=_file_hash(stage/files[0]),
                                      content_sha256=codec.content_sha256,
                                      bytes=sum((stage/member).stat().st_size for member in files), free_index_plan_stats=None))
                manifest = dict(expression.manifest, parts=parts, kind=AnalyticExpressionIR.artifact_kind,
                                ir_version=AnalyticExpressionIR.format_version,
                                storage_backend=storage_backend, storage_version=_AnalyticIRCodec.storage_version(storage_backend),
                                operator_keys=expression.manifest.get('operator_keys', sorted(operators)),
                                specializations=dict(active_operators=expression.manifest.get('active_operators'),
                                    mode_dimension=expression.manifest['options'].get('mode_dimension')))
                manifest['signature'] = _fingerprint((key, [part['sha256'] for part in parts]))
                manifest['source_legacy_signature'] = expression.manifest['signature']
                manifest['options'] = dict(manifest['options'], free_index_reduction=False)
                _atomic_json(stage/'manifest.json', manifest)
            if file is None: return AnalyticExpressionIR.load(directory)
            target = Path(file).expanduser().resolve()
            target.parent.mkdir(parents=True, exist_ok=True)
            archive = IndexedTarStore(Path(folder)/'build.tar')
            archive.write_trees([('expression', directory, {})])
            atomic_write(target, archive.save, keep_mode=False)
        return AnalyticExpressionIR.load(target)

    @staticmethod
    def _pack_readers(readers, file, *, staging_folder=None, storage_backend=None, single=False):
        import shutil
        import tempfile
        from McUtils.Scaffolding import IndexedTarStore
        if storage_backend is not None: _AnalyticIRCodec.storage_backend(storage_backend)
        if not readers or any(type(name) is not str or not name for name in readers):
            raise ValueError('an expression collection needs nonempty string names')
        if single and len(readers) != 1: raise ValueError('single-expression output needs one reader')
        kinds = {reader.artifact_kind for reader in readers.values()}
        if len(kinds) != 1: raise ValueError('cannot mix expression IR and compiled programs in a collection')
        if staging_folder is not None: Path(staging_folder).mkdir(parents=True, exist_ok=True)
        if file is None:
            if staging_folder is None: raise ValueError('a collection needs an output file or staging folder')
            manager = contextlib.nullcontext(tempfile.mkdtemp(prefix='analytic-library-', dir=staging_folder))
        else:
            manager = tempfile.TemporaryDirectory(prefix='analytic-library-', dir=staging_folder)
        with manager as temporary:
            root = Path(temporary)
            entries = {}
            for index, (name, reader) in enumerate(readers.items()):
                prefix = 'expression' if single else 'entries/{:06d}'.format(index)
                destination = root/prefix
                destination.mkdir(parents=True)
                backend = reader.storage_backend if storage_backend is None else storage_backend
                version = (reader.manifest.get('storage_version', 1) if storage_backend is None else
                           reader.codec.storage_version(backend))
                parts = []
                for part in reader.manifest['parts']:
                    record = None
                    if backend == reader.storage_backend and version == reader.manifest.get('storage_version', 1):
                        buffers = part.get('buffers')
                        if buffers is None:
                            record = reader._read_record(part)
                            buffers = sorted({node['file'] for node in record['nodes'] if node['kind'] == 'array'})
                        for member in (part['file'], *buffers):
                            with reader._open(member) as source, (destination/reader._name(member)).open('wb') as target:
                                shutil.copyfileobj(source, target, length=1024*1024)
                        if _file_hash(destination/part['file']) != part['sha256']:
                            raise ValueError('analytic part checksum mismatch while copying')
                        digest = part.get('content_sha256')
                        if digest is None:
                            if record is None: record = reader._read_record(part)
                            digest = reader.codec.content_hash(record)
                        parts.append(dict(part, buffers=list(buffers), content_sha256=digest))
                    else:
                        payload = reader._decode_record(reader._read_record(part))
                        codec = reader.codec(destination, Path(part['file']).stem)
                        members = codec.write(payload, storage_backend=backend)
                        parts.append(dict(part, file=members[0], buffers=members[1:],
                            sha256=_file_hash(destination/members[0]), content_sha256=codec.content_sha256,
                            bytes=sum((destination/member).stat().st_size for member in members)))
                _atomic_json(destination/'manifest.json', dict(reader.manifest, parts=parts,
                             storage_backend=backend, storage_version=version))
                entries[name] = prefix
            if not single:
                catalog = dict(kind=next(iter(kinds))+'-collection', catalog_version=1, complete=True, entries=entries)
                _atomic_json(root/'manifest.json', catalog)
            if file is None: return root/'expression' if single else root
            target = Path(file).expanduser().resolve()
            if any(target == reader.source for reader in readers.values()):
                raise ValueError('collection output must be separate from its inputs')
            target.parent.mkdir(parents=True, exist_ok=True)
            archive = IndexedTarStore(root/'build.tar')
            archive.write_trees([('expression', root/'expression', {})] if single else
                               [('entries', root/'entries', {}),
                                ('manifest.json', None, {'': (root/'manifest.json').read_bytes()})])
            atomic_write(target, archive.save, keep_mode=False)
        return target

    def pack_ir(self, sources, file, *, staging_folder=None, storage_backend=None):
        """Pack named, already simplified expressions into one indexed IR file."""
        with contextlib.ExitStack() as stack:
            readers = {name: (source if isinstance(source, AnalyticExpressionIR) else
                              stack.enter_context(AnalyticExpressionIR.load(source)))
                       for name, source in sources.items()}
            if any(reader.artifact_kind != AnalyticExpressionIR.artifact_kind for reader in readers.values()):
                raise ValueError('pack_ir requires simplified mathematical expressions')
            return self._pack_readers(readers, file, staging_folder=staging_folder, storage_backend=storage_backend)

    def _compile_ir_collection(self, sources, file, *, staging_folder=None, **options):
        import tempfile
        if file is not None:
            target = Path(file).expanduser().resolve()
            if any(target == (source.source if isinstance(source, AnalyticExpressionIR) else
                               Path(source).expanduser().resolve()) for source in sources.values()):
                raise ValueError('compiled output must be separate from its source IR')
        with tempfile.TemporaryDirectory(prefix='analytic-library-compile-') as temporary, contextlib.ExitStack() as stack:
            readers = {}
            for index, (name, source) in enumerate(sources.items()):
                compiled = self.compile_ir(source, staging_folder=Path(temporary)/str(index), **options)
                readers[name] = stack.enter_context(_CompiledAnalyticArchive.load(compiled.artifact))
            started = time.perf_counter()
            target = self._pack_readers(readers, file, staging_folder=staging_folder)
            self.timings['compilation_seconds'] += time.perf_counter()-started
        results = {name: CompiledAnalyticExpression.load(target, expression=name) for name in sources}
        for value in results.values(): value.store = self
        return results

    @staticmethod
    def _lower_ir_part(part, manifest, *, free_index_reduction=True, rank=None, orbit_cache=None):
        from . import Analytic as a
        num_fixed = manifest['num_fixed']
        baseline = manifest['baseline_shift']
        if free_index_reduction and not (baseline and any(baseline[num_fixed:])):
            _prepare_free_index_plan(part, num_fixed, orbit_cache, rank=rank)
        source_plans = {'native': IndexedBlockEvaluator.compile_plan(part)}
        if hasattr(part, '_free_index_plan_data'):
            source_plans['reduced'] = IndexedBlockEvaluator.compile_plan(part, reduced=True)
        max_order = a.PerturbationTheoryExpressionEvaluator._get_max_order(part)
        memo, programs, polynomial_cache = {}, {}, {}
        conditions, energy_keys, energy_programs, lowered_terms, lowered_subplans = {}, {}, {}, {}, {}
        occupation_shifts = tuple(np.arange(delta+1 if delta < 0 else 1,
                                             1 if delta < 0 else delta+1, dtype=np.int64)
                                  for delta in manifest['change'])

        def condition(value):
            identity = id(value)
            if identity not in conditions:
                name, arguments = _condition_program_fields(value)
                def prepare(argument):
                    from McUtils.Zachary.Symbolic.ElementaryFunctions import Functionlike
                    if isinstance(argument, a.OperatorCondition): return condition(argument)
                    if isinstance(argument, Functionlike): return argument.compile(mode='program')
                    if isinstance(argument, tuple): return tuple(prepare(item) for item in argument)
                    return argument
                arguments = tuple(prepare(argument) for argument in arguments)
                result = _restore_condition_program(name, arguments)
                for field in ('symmetric', 'diagonal', 'equivalence', 'requires_frequencies'):
                    object.__setattr__(result, field, getattr(value, field))
                object.__setattr__(result, 'key', ('compiled-condition', _fingerprint((name, arguments))))
                conditions[identity] = result
            return conditions[identity]

        def energy_key(key):
            if key not in energy_keys:
                energy_keys[key] = tuple((a.ConditionEndpoint(condition(row[0].condition), row[0].kind,
                                                   row[0].retain, row[0].endpoint),)+row[1:]
                             if isinstance(row[0], a.ConditionEndpoint) else row for row in key)
            return energy_keys[key]

        def lower_term(term):
            identity = id(term)
            if identity not in lowered_terms:
                key = energy_key(term.energy_changes)
                if key not in energy_programs:
                    energy_programs[key] = _CompiledEnergyProgram(a._compile_energy_program(key, baseline))
                lowered_terms[identity] = dataclasses.replace(term, polynomial=lower(term.polynomial),
                    energy_changes=key, denominator_program=energy_programs[key])
            return lowered_terms[identity]

        def lower(value):
            identity = id(value)
            if identity not in memo:
                if isinstance(value, a.PTEnergyChangeProductSum):
                    result = object.__new__(a.PTEnergyChangeProductSum)
                    result.terms = {energy_key(key): lower(poly) for key, poly in value.terms.items()}
                    result.prefactor, result.reduced, result._ndim = value.prefactor, value.reduced, value.ndim
                    result._shape = None
                else:
                    result = a._CompiledPolynomial.compile(value, cache=polynomial_cache)
                    result = programs.setdefault(result.signature, result)
                    result.occupation_shifts = occupation_shifts
                memo[identity] = result
            return memo[identity]

        indices = {key: tuple(part.get_inds(key)) for key in part.terms}
        free_groups, patterns = {}, {}
        for key, positions in indices.items():
            free_groups.setdefault(len(positions)-num_fixed, []).append(key)
            for coefficient in key:
                pattern = tuple(int(i) for i in coefficient[2:])
                if pattern:
                    table = patterns.setdefault(len(pattern), {})
                    table.setdefault(pattern, len(table))
        pattern_tables = {arity: (np.asarray(tuple(table), dtype=np.int64), table)
                          for arity, table in patterns.items()}
        maximum_factors = max(map(len, part.terms), default=0)
        execution = dict(version=1, max_order=max_order,
            groups=tuple((rank, tuple(keys)) for rank, keys in free_groups.items()),
            factorials=tuple(math.factorial(n) for n in range(maximum_factors+1)),
            coefficient_patterns=pattern_tables)
        plans = {}
        for name, source in source_plans.items():
            subexpressions = {}
            for key, subplan in source.subexpressions.items():
                if id(subplan) not in lowered_subplans:
                    groups = tuple((shape, tuple(lower_term(term) for term in terms))
                                   for shape, terms in subplan.groups)
                    lowered = lower(subplan.expression)
                    if isinstance(lowered, a.PTEnergyChangeProductSum):
                        lowered._compiled_terms = tuple(sorted((term for _, terms in groups for term in terms),
                                                               key=lambda term: term.ordinal))
                    lowered_subplans[id(subplan)] = IndexedSubexpressionPlan(lowered,
                        subplan.direct_polynomial, subplan.prefactor, groups)
                subexpressions[key] = lowered_subplans[id(subplan)]
            plans[name] = subexpressions, source.term_count, (
                None if source.index_bindings is None else dict(source.index_bindings))
        free_index = getattr(part, '_free_index_plan_data', None)
        if free_index is not None:
            free_index = dict(free_index, values={key: lower(value) for key, value in free_index['values'].items()})
        expression = object.__new__(a.PTTensorCoeffProductSum)
        expression.terms = {key: lower(value) for key, value in part.terms.items()}
        expression.prefactor, expression.reduced, expression._ndim = part.prefactor, part.reduced, part.ndim
        expression._shape, expression._inds_map = None, indices
        requirements, sides, energy_terms = {}, {}, {}
        identification = a.DegeneracyIdentificationContext([])
        for value in expression.terms.values():
            if not isinstance(value, a.PTEnergyChangeProductSum): continue
            for key in value.terms:
                if key in sides: continue
                sides[key] = tuple(a.PTEnergyChangeProductSum.side_change_iter(key))
                for _, change in sides[key]:
                    for method in ('linear', 'perfect'):
                        identity = method, tuple(change)
                        if identity in requirements: continue
                        requirement = identification._requirement(change, method)
                        maps = requirement[1]
                        projections = None
                        if method == 'linear':
                            checked = np.flatnonzero(maps[0] >= 0)
                            projections = np.stack([checked[np.argsort(mapping[checked])] for mapping in maps])
                        requirements[identity] = requirement, projections
        energy_terms = {key: tuple(value.terms) for key, value in expression.terms.items()
                        if isinstance(value, a.PTEnergyChangeProductSum)}
        execution['degeneracy'] = dict(sides=sides, requirements=requirements, terms=energy_terms)
        a._PolynomialProgramTable.pack(tuple(programs.values()))
        return dict(expression=expression, execution=execution, plans=plans, free_index=free_index)

    def compile_ir(self, source, file=None, *, staging_folder=None, free_index_reduction=True,
                   storage_backend='portable', progress=None):
        """Compile portable IR without a solver or any physical system inputs.

        The output is an independent self-contained archive (or a retained
        staging directory). All expression planning occurs before publication.
        """
        import tempfile
        from McUtils.Scaffolding import IndexedTarStore
        storage_backend = _CompiledAnalyticCodec.storage_backend(storage_backend)
        if file is None and staging_folder is None:
            raise ValueError('IR compilation needs an output file or staging folder')
        if not isinstance(free_index_reduction, (bool, np.bool_)):
            raise TypeError('free_index_reduction must be boolean')
        if isinstance(source, Mapping):
            return self._compile_ir_collection(source, file, staging_folder=staging_folder,
                free_index_reduction=free_index_reduction, storage_backend=storage_backend, progress=progress)
        if not isinstance(source, AnalyticExpressionIR):
            if file is not None and Path(source).expanduser().resolve() == Path(file).expanduser().resolve():
                raise ValueError('compiled output must be separate from its source IR')
            names = AnalyticExpressionIR.expression_names(source)
            if names:
                with contextlib.ExitStack() as stack:
                    sources = {name: stack.enter_context(AnalyticExpressionIR.load(source, expression=name)) for name in names}
                    return self._compile_ir_collection(sources, file, staging_folder=staging_folder,
                        free_index_reduction=free_index_reduction, storage_backend=storage_backend, progress=progress)
        owner = not isinstance(source, AnalyticExpressionIR)
        reader = AnalyticExpressionIR.load(source) if owner else source
        started = time.perf_counter()
        self.timings['requests'] += 1
        completed = False
        try:
            semantic = {key: reader.manifest[key] for key in
                        ('change', 'baseline_shift', 'num_fixed', 'projection', 'active_operators',
                         'operator_letter_bound', 'maximum_mode_indices', 'specializations')}
            semantic['operator_keys'] = reader.manifest.get('operator_keys')
            if semantic['active_operators'] is not None and semantic['operator_keys'] is None:
                raise ValueError('specialized IR needs operator-domain metadata; re-export it')
            signature = _fingerprint(dict(version=_CompiledAnalyticArchive.format_version,
                parts=[part.get('content_sha256', part['sha256']) for part in reader.manifest['parts']], semantic=semantic,
                free_index_reduction=bool(free_index_reduction)))
            manager = (tempfile.TemporaryDirectory(prefix='analytic-compile-') if staging_folder is None
                       else contextlib.nullcontext(str(Path(staging_folder).expanduser().resolve())))
            with manager as folder:
                entries = _CompiledEntries(Path(folder) / 'compiled')
                storage_key = _fingerprint((signature, storage_backend))
                directory = Path(entries.entry_path(storage_key))
                with entries.stage(storage_key) as staged:
                    stage = Path(staged.path)
                    parts = []
                    orbit_cache = BoundedCache.sampled(32768, 32*1024**2) if free_index_reduction else None
                    for original, expression in zip(reader.manifest['parts'], reader.iter_parts()):
                        payload = self._lower_ir_part(expression, reader.manifest,
                                                      free_index_reduction=free_index_reduction,
                                                      rank=original['rank'], orbit_cache=orbit_cache)
                        name = 'part-{:06d}-r{}'.format(len(parts), original['rank'])
                        codec = _CompiledAnalyticCodec(stage, name)
                        files = codec.write(payload, storage_backend=storage_backend)
                        parts.append(dict(original, file=files[0], buffers=files[1:], sha256=_file_hash(stage/files[0]),
                                          content_sha256=codec.content_sha256,
                                          bytes=sum((stage/member).stat().st_size for member in files)))
                        if progress is not None:
                            progress(dict(stage='lowering', parts=len(parts), source_parts=len(reader.manifest['parts'])))
                    manifest = dict(semantic, format=self.format_version, kind=_CompiledAnalyticArchive.artifact_kind,
                        plan_version=_CompiledAnalyticArchive.format_version, complete=True, signature=signature,
                        storage_backend=storage_backend, storage_version=_AnalyticIRCodec.storage_version(storage_backend),
                        program_versions=dict(polynomial=3, function=1, denominator=2, degeneracy=1),
                        parts=parts, source_ir_signature=reader.manifest['signature'],
                        options=dict(free_index_reduction=bool(free_index_reduction),
                                     mode_dimension=semantic['specializations']['mode_dimension']),
                        seconds=time.perf_counter()-started)
                    for key in ('request_version', 'request_signature'):
                        if key in reader.manifest: manifest[key] = reader.manifest[key]
                    _atomic_json(stage/'manifest.json', manifest)
                if file is None:
                    result = CompiledAnalyticExpression.load(directory)
                else:
                    target = Path(file).expanduser().resolve()
                    if target == reader.source:
                        raise ValueError('compiled output must be separate from its source IR')
                    target.parent.mkdir(parents=True, exist_ok=True)
                    archive = IndexedTarStore(Path(folder)/'build.tar')
                    archive.write_trees([('expression', directory, {})])
                    atomic_write(target, archive.save, keep_mode=False)
                    result = CompiledAnalyticExpression.load(target)
                result.store = self
                completed = True
                return result
        finally:
            if owner: reader.close()
            self.timings['compilation_seconds'] += time.perf_counter()-started
            self.timings['compilations' if completed else 'failures'] += 1

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
                coalesce_combine_nodes=True,bounded_byte_accounting=True,free_index_reduction=True,
                file=None,staging_folder=None,storage_backend=None):
        """Compile native expressions; early DAG combine-node coalescing is on by default.

        Pass an IR path/reader plus ``file=...`` for standalone portable
        compilation. A native expression plus ``file=...`` runs the same two
        explicit stages as a convenience. Omitting a destination preserves the
        existing in-memory/legacy-cache interface.

        Set coalesce_combine_nodes=False (also accepted through evaluator
        compilation_options) to disable that optimization for diagnostics or
        benchmarking. Bounded final-part accumulation remains active.
        """
        if self.using_artifacts and file is None and staging_folder is None:
            if storage_backend is not None:
                raise ValueError('select storage_backend on the artifact session')
            return self._compile_artifact(lambda: expression,
                context=expression if context is None else context, active_operators=active_operators,
                projection=projection, mode_dimension=mode_dimension, free_index_reduction=free_index_reduction,
                max_terms=max_terms, max_bytes=max_bytes, cache_size=cache_size, cache_bytes=cache_bytes,
                node_cache_size=node_cache_size, progress=progress,
                coalesce_combine_nodes=coalesce_combine_nodes, bounded_byte_accounting=bounded_byte_accounting)
        if storage_backend is not None:
            _AnalyticIRCodec.storage_backend(storage_backend)
            if file is None and staging_folder is None:
                raise ValueError('storage_backend needs an output file or staging folder')
        backend = 'portable' if storage_backend is None else storage_backend
        if isinstance(expression, (AnalyticExpressionIR, str, os.PathLike)):
            if active_operators is not None or projection is not None or mode_dimension is not None:
                raise ValueError('IR specializations are fixed at export')
            return self.compile_ir(expression, file, staging_folder=staging_folder,
                free_index_reduction=free_index_reduction, storage_backend=backend, progress=progress)
        if file is not None or staging_folder is not None:
            import tempfile
            with tempfile.TemporaryDirectory(prefix='analytic-ir-compile-') as temporary:
                with self.export_ir(expression, staging_folder=temporary, context=context, storage_backend=backend,
                    active_operators=active_operators, projection=projection, mode_dimension=mode_dimension,
                    max_terms=max_terms, max_bytes=max_bytes, cache_size=cache_size, cache_bytes=cache_bytes,
                    node_cache_size=node_cache_size, coalesce_combine_nodes=coalesce_combine_nodes,
                    bounded_byte_accounting=bounded_byte_accounting, progress=progress) as ir:
                    return self.compile_ir(ir, file, staging_folder=staging_folder,
                        free_index_reduction=free_index_reduction, storage_backend=backend, progress=progress)
        return self.compile_factory(lambda:expression, context=expression if context is None else context,
            active_operators=active_operators,projection=projection,max_terms=max_terms,max_bytes=max_bytes,
            cache_size=cache_size,cache_bytes=cache_bytes,node_cache_size=node_cache_size,
            mode_dimension=mode_dimension,progress=progress,coalesce_combine_nodes=coalesce_combine_nodes,
            bounded_byte_accounting=bounded_byte_accounting,free_index_reduction=free_index_reduction)

    def compile_factory(self, factory, **options):
        """Compile or reuse an entry, accounting separately for lookup and actual construction."""
        if self.using_artifacts:
            return self._compile_artifact(factory, **options)
        backend = options.pop('storage_backend', None)
        if backend is not None:
            _AnalyticIRCodec.storage_backend(backend)
            raise ValueError('storage_backend needs compile(..., file=...) or export_ir/compile_ir')
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
                          coalesce_combine_nodes=True,bounded_byte_accounting=True,free_index_reduction=True,
                          part_format='pickle',storage_backend='portable'):
        from .Analytic import (SqrtChangePoly,PTTensorCoeffProductDAG,PTTensorCoeffProductSum,
            PTTensorCoeffProductDAGEvaluationPlan,_BalancedPolynomialAccumulator)
        if not isinstance(coalesce_combine_nodes, (bool, np.bool_)):
            raise TypeError('coalesce_combine_nodes must be boolean')
        coalescing_enabled = bool(coalesce_combine_nodes)
        if not isinstance(free_index_reduction,(bool,np.bool_)):
            raise TypeError('free_index_reduction must be boolean')
        free_index_reduction=bool(free_index_reduction)
        if part_format not in {'pickle', 'ir'}:
            raise ValueError('unknown analytic part format')
        if part_format == 'ir' and (self.folder is None or free_index_reduction):
            raise ValueError('IR building requires staging and no execution plans')
        if part_format == 'ir': _AnalyticIRCodec.storage_backend(storage_backend)
        options = dict(max_terms=max_terms,max_bytes=max_bytes,cache_size=cache_size,
                       cache_bytes=cache_bytes,node_cache_size=node_cache_size,mode_dimension=mode_dimension,
                       coalesce_combine_nodes=coalescing_enabled,bounded_byte_accounting=bool(bounded_byte_accounting),
                       free_index_reduction=free_index_reduction)
        if part_format == 'ir': options.update(part_format='ir', storage_backend=storage_backend)
        if max_terms<1 or (max_bytes is not None and max_bytes<1):raise ValueError('invalid compilation budget')
        signature = self._signature(context,active_operators,projection,options)
        directory = None if self.folder is None else self.folder/signature
        if self.folder is None:
            cached=self._get_memory(signature)
            if cached is not None:return cached
        elif self.entries.has(signature):return directory if part_format == 'ir' else CompiledAnalyticExpression.load(directory)
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
            change = tuple(int(i) for i in wrapper.change) if wrapper.change is not None else ()
            baseline = None;poly = PTTensorCoeffProductSum({})
        elif isinstance(expr,SqrtChangePoly):
            change = tuple(int(i) for i in expr.poly_change)
            baseline = tuple(int(i) for i in expr.shift_start);poly = expr.poly_obj
        else:
            change = tuple(int(i) for i in wrapper.change) if wrapper.change is not None else ()
            baseline = None;poly = expr
        if not isinstance(poly,PTTensorCoeffProductSum):poly = PTTensorCoeffProductSum({():poly})
        operator_keys = (poly.operator_keys if isinstance(poly, PTTensorCoeffProductDAG) else
                         {coefficient[:2] for key in poly.terms for coefficient in key})
        if projection is not None and (any(change) or baseline and any(baseline)):
            raise ValueError('final occupation projection requires a zero-change numerator')
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
        projected = (isinstance(poly,PTTensorCoeffProductDAG)
                     and poly._node[0] == 'occupation_coefficient'
                     and getattr(projection,'kind',None) == 'coefficient'
                     and poly._node[2:] == ((projection.powers,projection.shift) if not (
                         projection.fixed_powers or projection.preserve_domains) else
                         (projection.powers,projection.shift,projection.fixed_powers,projection.preserve_domains)))
        transform = projection.compiler() if projection is not None and not projected else None
        constant_transform = (_OccupationCoefficientProjection((),0.,preserve_domains=projection.preserve_domains).compiler()
                              if getattr(projection,'kind',None) == 'coefficient' else None)
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
                totals = {key:acc.total() for key,acc in pending.items()}
                if constant_transform is not None:
                    totals = {key:value for key,total in totals.items()
                              if (value := constant_transform(key,total)) is not None}
                compact = PTTensorCoeffProductSum(totals,
                    canonicalize=False,reduced=False).combine()
                groups = {}
                for key,value in compact.terms.items():
                    rank = len({i for coefficient in key for i in coefficient[2:]})
                    groups.setdefault(rank,{})[key] = value
                for rank,terms in sorted(groups.items()):
                    name = 'part-{:06d}-r{}.{}'.format(len(parts),rank,'json' if part_format == 'ir' else 'pkl')
                    part = PTTensorCoeffProductSum(terms,canonicalize=False,reduced=True)
                    if can_reduce:
                        _prepare_free_index_plan(part,wrapper.num_fixed,self._free_index_orbits,rank=rank)
                    if part_format == 'ir':
                        codec = _AnalyticIRCodec(stage, Path(name).stem)
                        files = codec.write(part, storage_backend=storage_backend)
                        name = files[0]
                        part_hash = _file_hash(stage / name)
                        part_bytes = sum((stage / member).stat().st_size for member in files)
                    elif stage is None:
                        payload=pickle.dumps(part,protocol=5)
                        part_hash=hashlib.sha256(payload).hexdigest();part_bytes=len(payload)
                        memory_parts[name]=part
                    else:
                        _atomic_pickle(stage/name,part)
                        part_hash=_file_hash(stage/name);part_bytes=(stage/name).stat().st_size
                    parts.append(dict(file=name,rank=rank,sha256=part_hash,bytes=part_bytes,
                        coefficient_keys=len(terms),energy_rows=sum(len(p.terms) if hasattr(p,'terms') else 1 for p in terms.values()),
                        free_index_plan_stats=(getattr(part,'_free_index_plan_data',{}).get('stats'))))
                    if part_format == 'ir': parts[-1].update(content_sha256=codec.content_sha256, buffers=files[1:])
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
                operator_keys=sorted(operator_keys),
                operator_letter_bound=degree,maximum_mode_indices=maximum,options=options,
                seconds=time.perf_counter()-started,stream_stats=budget.stats(),
                plan_stats=None if plan is None else plan.stats(),implementation=_implementation_signature())
            if part_format == 'ir':
                manifest.update(kind='analytic-expression-ir', ir_version=AnalyticExpressionIR.format_version,
                                storage_backend=storage_backend, storage_version=_AnalyticIRCodec.storage_version(storage_backend),
                                specializations=dict(active_operators=manifest['active_operators'],
                                                     mode_dimension=mode_dimension))
            if stage is not None:_atomic_json(stage/'manifest.json',manifest)
        if part_format == 'ir': return directory
        compiled=(self._remember(CompiledAnalyticExpression(None,manifest,parts=memory_parts))
                  if directory is None else CompiledAnalyticExpression.load(directory))
        if projection is None:
            native=[wrapper.op if wrapper.op is not None else poly]
            from .Analytic import CorrectionView
            if isinstance(wrapper.op,CorrectionView) and not any(getattr(wrapper.op,name) is not None for name in
                    ['allowed_terms','allowed_coefficients','disallowed_coefficients','allowed_energy_changes']):
                native.append(wrapper.op.term)
            for op in native:
                native_context = dict(operator=op,change=change,baseline=baseline)
                if isinstance(context,dict) and 'excluded_resonances' in context:
                    native_context['excluded_resonances'] = context['excluded_resonances']
                key=self._signature(native_context,active_operators,None,options)
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
    def __init__(self,folder,manifest,parts=None,artifact=None,expression=None):
        self.folder=None if folder is None else Path(folder)
        self.manifest=manifest;self.parts=parts;self.store=None
        self.artifact = None if artifact is None else Path(artifact)
        self.expression = expression

    @classmethod
    def expression_names(cls, source):
        return _CompiledAnalyticArchive.expression_names(source)

    def iter_parts(self):
        """Read bounded parts from a compatible legacy cache or portable file."""
        if self.artifact is not None:
            with _CompiledAnalyticArchive.load(self.artifact, expression=self.expression) as reader:
                yield from reader.iter_parts()
            return
        for part in self.manifest['parts']:
            if self.parts is not None:
                yield self.parts[part['file']]
            else:
                path = self.folder/part['file']
                if _file_hash(path) != part['sha256']: raise ValueError('compiled part checksum mismatch')
                with path.open('rb') as stream: yield pickle.load(stream)

    def save(self, file, *, free_index_reduction=True, storage_backend=None):
        """Publish one self-contained compiled file, including legacy migration."""
        if self.artifact is not None:
            with _CompiledAnalyticArchive.load(self.artifact, expression=self.expression) as reader:
                target = AnalyticExpressionStore._pack_readers({self.expression or 'expression':reader}, file,
                    storage_backend=storage_backend, single=self.expression is None)
            return type(self).load(target)
        import tempfile
        store = AnalyticExpressionStore()
        with tempfile.TemporaryDirectory(prefix='analytic-cache-migration-') as temporary:
            backend = 'portable' if storage_backend is None else storage_backend
            with store.export_ir(self, staging_folder=temporary, storage_backend=backend) as ir:
                return store.compile_ir(ir, file, free_index_reduction=free_index_reduction,
                                        storage_backend=backend)
    @classmethod
    def load(cls,folder,*,expression=None,validate_source=True):
        folder=Path(folder).resolve()
        if folder.is_file():
            with _CompiledAnalyticArchive.load(folder, expression=expression) as reader:
                return cls(None, reader.manifest, artifact=folder, expression=reader.expression)
        manifest=json.loads((folder/'manifest.json').read_text())
        if manifest.get('kind') in {_CompiledAnalyticArchive.artifact_kind, _CompiledAnalyticArchive.artifact_kind+'-collection'}:
            with _CompiledAnalyticArchive.load(folder, expression=expression) as reader:
                return cls(folder, reader.manifest, artifact=folder, expression=reader.expression)
        if manifest.get('kind') == 'analytic-expression-ir':
            raise ValueError('simplified expression IR must be compiled before evaluation')
        if not manifest.get('complete') or manifest.get('format')!=AnalyticExpressionStore.format_version:
            raise ValueError('incomplete or unsupported compiled expression')
        if validate_source and manifest['implementation']!=_implementation_signature():
            raise ValueError('compiled expression source version changed')
        for part in manifest['parts']:
            if Path(part['file']).name!=part['file']:raise ValueError('invalid expression part path')
            if _file_hash(folder/part['file'])!=part['sha256']:raise ValueError('compiled expression checksum mismatch')
        return cls(folder,manifest)
    def evaluate_coefficient(self, coeffs, freqs, *, indices=(), **options):
        """Contract a symbolic action coefficient with physical tensors.

        External mode indices bind the requested action powers. The constant
        numerator uses an internal zero index carrier, with no reference-state
        binding in the expression or its symbolic cache.
        """
        projection = self.manifest['projection']
        if projection is None or projection.get('kind') != 'coefficient':
            raise ValueError('expression is not an action coefficient kernel')
        indices = tuple(indices)
        if any(not isinstance(i,(int,np.integer)) for i in indices):
            raise ValueError('coefficient mode indices must be integers')
        carrier = [(np.zeros(len(freqs),dtype=int),np.arange(len(freqs))[None])]
        return self.evaluate(carrier,coeffs,freqs,anchored_modes=indices,**options)[0][:,0]

    def prepare_evaluation(self,state_perms,coeffs,freqs,*,anchored_modes=(),mode_subset_selector=None,parallelizer=None,
                           batch_size=512,task_size=32768,resume=True,use_shared_memory=None,
                           evaluation_folder=None,**numerical_options):
        started=time.perf_counter()
        try:
            return CompiledAnalyticEvaluation(self,state_perms,coeffs,freqs,anchored_modes=anchored_modes,mode_subset_selector=mode_subset_selector,
                parallelizer=parallelizer,batch_size=batch_size,task_size=task_size,resume=resume,
                numerical_options=numerical_options,use_shared_memory=use_shared_memory,
                evaluation_folder=evaluation_folder)
        finally:
            if self.store is not None:self.store.timings['preparation_seconds']+=time.perf_counter()-started
    def evaluate(self,state_perms,coeffs,freqs,**options):
        return self.prepare_evaluation(state_perms,coeffs,freqs,**options).run()


class CompiledAnalyticEvaluation:
    """A reusable numerical task plan; the parallelizer never sees a symbolic DAG."""
    def __init__(self,compiled,state_perms,coeffs,freqs,*,anchored_modes=(),mode_subset_selector=None,parallelizer=None,
                 batch_size=512,task_size=32768,resume=True,numerical_options=None,use_shared_memory=None,
                 evaluation_folder=None):
        self.compiled=compiled
        self.state_perms=[(np.asarray(state,dtype=int),np.asarray(perms,dtype=int)) for state,perms in state_perms]
        from McUtils.Parallelizers import SharedMemoryArrayTree
        from .Analytic import PerturbationTheoryExpressionEvaluator
        active = _current_evaluation_plan.get()
        self.coeffs = (coeffs if isinstance(coeffs,SharedMemoryArrayTree) else
            PerturbationTheoryExpressionEvaluator._normalize_coefficient_expansions(
                coeffs,preparation_context=None if active is None else active[1]))
        self.freqs=np.asarray(freqs);self.anchored_modes=tuple(anchored_modes)
        active_operators = compiled.manifest.get('active_operators')
        if compiled.artifact is not None and active_operators is not None:
            relevant = compiled.manifest.get('operator_keys')
            if relevant is None: raise ValueError('compiled specialization lacks its operator domain')
            coefficients = self.coeffs.tree if isinstance(self.coeffs, SharedMemoryArrayTree) else self.coeffs
            unexpected = (set(_active_coefficient_operators(coefficients)) & set(map(tuple, relevant))
                          - set(map(tuple, active_operators)))
            if unexpected:
                raise ValueError('compiled specialization excludes active coefficient operators: {}'.format(sorted(unexpected)))
        self.mode_subset_selector=mode_subset_selector or _AnchoredModeSubset(self.anchored_modes)
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
        if projection is not None and projection.get('kind') == 'coefficient':
            if (len(self.anchored_modes) != len(projection['powers'])
                    or len(set(self.anchored_modes)) != len(self.anchored_modes)
                    or any(i < compiled.manifest['num_fixed'] or i >= len(self.freqs) for i in self.anchored_modes)):
                raise ValueError('coefficient kernel requires distinct free mode positions for its powers')
        elif projection is not None:
            if len(self.anchored_modes)!=1:raise ValueError('projected transition requires its one anchored physical mode')
            anchor=self.anchored_modes[0]
            for state,perms in self.state_perms:
                expected=np.full(len(self.freqs),projection['background'],dtype=int);expected[anchor]=projection['after']
                if not np.array_equal(state,expected) or not np.array_equal(perms,np.arange(len(self.freqs))[None]):
                    raise ValueError('occupation-specialized expression used outside its target-state domain')
        elif self.anchored_modes:raise ValueError('anchoring requires a matching final occupation projection')
        signature=self._fingerprint_context(dict(expression=compiled.manifest['signature'],states=self.state_perms,
            coefficients=self.coeffs,frequencies=self.freqs,anchors=self.anchored_modes,subsets=self.mode_subset_selector,numerical=self.numerical_options,
            batch_size=self.batch_size,task_size=self.task_size))
        if evaluation_folder is None and active is not None:
            evaluation_folder = active[1].root.options.get('evaluation_folder')
        destination = None if evaluation_folder is None else Path(evaluation_folder).expanduser().resolve()
        if destination is not None and compiled.folder is not None:
            root = compiled.folder.resolve() if compiled.store is None or compiled.store.folder is None else compiled.store.folder
            if destination == root or root in destination.parents:
                raise ValueError('evaluation_folder must be outside the compiled expression store')
        self.folder=None if destination is None else destination/signature
        if self.folder is not None:self.folder.mkdir(parents=True,exist_ok=True)
        tasks=[];ndim=len(self.freqs);num_fixed=compiled.manifest['num_fixed']
        for part in compiled.manifest['parts']:
            rank=part['rank'];free=rank-num_fixed
            if free<0:raise ValueError('coefficient aliases do not cover the fixed change')
            count=self.mode_subset_selector.count(ndim,num_fixed,free)
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
            expression_archive=None if self.compiled.artifact is None else str(self.compiled.artifact),
            expression_name=self.compiled.expression,
            evaluation_folder=None if self.folder is None else str(self.folder),signature=self.signature,state_perms=self.state_perms,
            coeffs=self.coeffs if coeffs is None else coeffs,freqs=self.freqs,anchors=self.anchored_modes,subset_selector=self.mode_subset_selector,
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
    source = context.get('expression_archive')
    manager = (contextlib.nullcontext(None) if source is None else
               _CompiledAnalyticArchive.load(source, expression=context.get('expression_name')))
    with manager as archive:
        return _run_compiled_tasks_bound(tasks, context, archive)


def _run_compiled_tasks_bound(tasks,context,archive):
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
            if archive is not None:
                cache[key] = archive.read_part(part)
            elif context.get('expression_parts') is not None:
                cache[key]=context['expression_parts'][key]
            else:
                path=Path(context['expression_folder'])/key
                if _file_hash(path)!=part['sha256']:raise ValueError('compiled part changed in worker')
                with path.open('rb') as stream:cache[key]=pickle.load(stream)
            # Compiled parts are immutable for the lifetime of this worker.
            # Compute invariant polynomial degree once per load, not per batch.
            execution = getattr(cache[key], '_compiled_execution', None)
            cache[key]=(cache[key], execution['max_order'] if execution is not None else
                        Expression._get_max_order(cache[key])
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
        iterator=itertools.islice(context['subset_selector'].iter_subsets(ndim,num_fixed,free),position,task['stop'])
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
        for name in ('path_equality_cache_size','cache_polynomial_order','free_index_reduction','prune_fixed_resonances'):
            options.pop(name,None)
        options['_max_polynomial_order']=max_order
        options['_poly_eval_cache']=polynomial_cache
        options['_stream_preparation']=preparation
        for extra_batch in iter(lambda:list(itertools.islice(iterator,context['batch_size'])),[]):
            subsets=[tuple(extra) for extra in extra_batch]
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
        if not self.store.using_artifacts:
            options.setdefault('mode_dimension',len(evaluator.freqs))
        else:
            active = None
        for key in ('ir_file', 'compiled_file', 'storage_backend', 'reuse_only'): options.pop(key, None)
        source=dict(solver=evaluator.solver,correction=correction,energy_method=energy_method)
        generator=[None]
        def expression(n,change):
            if generator[0] is None:generator[0]=self._generator(correction,order,energy_method)
            term=generator[0](n)
            return PerturbationTheoryExpressionEvaluator(term,term.get_poly_terms(change),change)
        return {(n,tuple(change)):self.store.compile_factory(lambda n=n,c=tuple(change):expression(n,c),
                    context=dict(source,order=n,change=tuple(change)),active_operators=active,projection=projection,**options)
                for n in orders for change in changes}

    def _action_terms(self, order, correction):
        return self.evaluator.solver.action_expansion_terms(order,correction)

    def precompile_action_expansion(self, order=2, correction='energy', *, compilation_options=None):
        from .Analytic import PerturbationTheoryExpressionEvaluator
        terms = self._action_terms(order,correction)
        source = self.evaluator.solver._action_expansion_source(correction,order)
        options = dict(self.evaluator.compilation_options)
        options.update(compilation_options or {})
        if not self.store.using_artifacts:
            options.setdefault('mode_dimension',len(self.evaluator.freqs))
        for key in ('ir_file', 'compiled_file', 'storage_backend', 'reuse_only'): options.pop(key, None)
        return {key:self.store.compile(
                    PerturbationTheoryExpressionEvaluator(None,kernel,(),logger=False),
                    context=dict(solver=source,order=key[0],powers=key[1],source=correction),
                    projection=_OccupationCoefficientProjection(key[1]),active_operators=None,**options)
                for key,kernel in terms.items()}

    def action_expansion(self, order=2, correction='energy', *, parallelizer=None, **options):
        from .Analytic import PerturbationTheoryExpressionEvaluator
        evaluator = self.evaluator
        kernels = self._action_terms(order,correction)
        numerical = dict(evaluator.evaluation_options);numerical.update(options)
        numerical.pop('evaluation_mode',None)
        if any(numerical.get(k) is not None and numerical.get(k) is not False for k in (
                'degenerate_changes','projection_bindings','only_degenerate_terms','allowed_coefficients',
                'disallowed_coefficients','allowed_energy_changes','epaths')):
            raise ValueError('resonance-independent action expansion cannot use runtime diagram filters')
        numerical.pop('only_degenerate_terms',None)
        numerical.setdefault('zero_cutoff',0.)
        prepared = evaluator._prep_expansions(evaluator.expansions)
        freqs = evaluator.freqs
        ndim = len(freqs)
        backend = evaluator.parallelizer if parallelizer is None else parallelizer
        maximum_degree = max(sum(powers) for _,powers in kernels)
        tensors = [0.] + [np.zeros((ndim,)*degree) for degree in range(1,maximum_degree+1)]
        compiled = (self.precompile_action_expansion(order,correction)
                    if evaluator.compile_expressions else None)
        for key,kernel in kernels.items():
            _,powers = key
            indices = ([()] if not powers else [(i,) for i in range(ndim)] if len(powers)==1
                       else itertools.combinations(range(ndim),len(powers)))
            for anchors in indices:
                if compiled is not None:
                    value = compiled[key].evaluate_coefficient([prepared],freqs,indices=anchors,
                                                              parallelizer=backend,**numerical)[0]
                else:
                    # The projected kernels are constant in occupations. This
                    # index carrier only reuses the native tensor contraction;
                    # no reference energy or state was used to derive them.
                    carrier = [(np.zeros(ndim,dtype=int),np.arange(ndim)[None])]
                    wrapper = PerturbationTheoryExpressionEvaluator(None,kernel,(),logger=False)
                    value = wrapper.evaluate(carrier,[prepared],freqs,
                        mode_subset_selector=_AnchoredModeSubset(tuple(anchors)),
                        parallelizer=backend,evaluation_mode='direct',**numerical)[0][0,0]
                if not powers:
                    tensors[0] += value
                else:
                    indices = tuple(mode for mode,power in zip(anchors,powers) for _ in range(power))
                    tensors[sum(powers)][indices] += value
        return tensors,kernels

    @staticmethod
    def _resonance_subsets(wrapper,state,permutation,groups):
        """Necessary integer mode support of every retained nonzero path.

        A PT vertex of order r has at most r+2 <= 3r operator letters.
        Distances along source->u->v->target therefore give a conservative
        order bound. One property vertex can supply its base rank for free.
        For diagonal VPT2 this reduces exactly to the previous H1 R(F1)
        source-neighbor schedule, including its odd one/three-quanta test.
        """
        term=wrapper.op
        while getattr(term,'order',None) is None and hasattr(term,'term'):term=term.term
        order=getattr(term,'order',None)
        observable=getattr(term,'type',0)
        if not isinstance(observable,(int,np.integer)):observable=None
        state=np.asarray(state,dtype=int);target=state.copy()
        target[np.asarray(permutation[:wrapper.num_fixed],dtype=int)]+=wrapper.expr.poly_change
        def cost(distance,rank=0):
            value=max(0,(int(distance)-rank+2)//3)
            return value+((int(distance)-rank-value)%2)
        inverse=np.argsort(permutation);fixed=set(permutation[:wrapper.num_fixed]);supports=set()
        for group in groups:
            for u in group:
                u=np.asarray(u,dtype=int)
                for v in group:
                    if np.array_equal(u,v):continue
                    v=np.asarray(v,dtype=int)
                    distances=(np.sum(abs(state-u)),np.sum(abs(u-v)),np.sum(abs(v-target)))
                    if order is not None and observable is not None:
                        costs=[cost(d) for d in distances]
                        bound=min(sum(costs)-costs[i]+cost(distances[i],observable) for i in range(3))
                        if bound>order:continue
                    physical=set(np.flatnonzero((state!=u)|(state!=v)))-fixed
                    supports.add(tuple(sorted(int(inverse[i]) for i in physical)))
        return _AnchoredModeSubset((),tuple(sorted(supports)))

    @staticmethod
    def evaluate_action_expansion(wrapper,state_perms,coeffs,freqs,*,
                                  pruning_signature=b'',**numerical):
        from .Analytic import (PTTensorCoeffProductDAG,SqrtChangePoly,
                               PerturbationTheoryExpressionEvaluator as Expression)
        expression=wrapper.expr
        if (not isinstance(expression,SqrtChangePoly)
                or not isinstance(expression.poly_obj,PTTensorCoeffProductDAG)
                or not expression.poly_obj.has_change_domains):
            return AnalyticCorrectionCompiler._evaluate_action_coefficients(
                wrapper,state_perms,coeffs,freqs,pruning_signature=pruning_signature,**numerical)
        bindings=numerical.get('projection_bindings')
        if bindings is None:raise ValueError('filtered action polynomials require resonance-change bindings')
        common,residual=expression.split_change_domains()
        term=wrapper.op
        while getattr(term,'order',None) is None and hasattr(term,'term'):term=term.term
        if (type(term).__name__ in {'EnergyCorrection','BlockHamiltonianCorrection'}
                and not wrapper.change):
            # Scalar eigenvalues are independent of wavefunction normalization.
            # Use the same native scalar recurrence in intermediate normalization
            # for the unrestricted part, avoiding block-normalization expressions
            # that cancel there. The exact normalized residual above is retained.
            common=term.parent._action_expansion_source('energy',term.order).energy_correction(
                term.order).get_poly_terms(())
        output=[np.zeros((len(coeffs),len(perms))) for _,perms in state_perms]
        if not isinstance(common,(int,float)):
            base=Expression(None,common,wrapper.change,logger=wrapper.logger)
            options=dict(numerical,projection_bindings=None,degenerate_changes=None)
            bound=(getattr(term,'order',0)//2+1
                   if type(term).__name__ in {'EnergyCorrection','BlockHamiltonianCorrection'} else None)
            output=AnalyticCorrectionCompiler._evaluate_action_coefficients(
                base,state_perms,coeffs,freqs,max_degree=bound,**options)
        if isinstance(residual,(int,float)):return output
        groups=bindings['projectors'][0]
        if not groups:return output
        delta=Expression(None,residual,wrapper.change,logger=wrapper.logger)
        numerical=dict(numerical)
        # Select sparse integer supports before numerical subset dispatch.
        # Kernels are not reprojected and the full filtered expression is never
        # contracted as an alternative implementation.
        numerical['evaluation_mode']='direct'
        compilation=current_compilation()
        for row,(state,permutations) in enumerate(state_perms):
            for column,permutation in enumerate(permutations):
                selector=AnalyticCorrectionCompiler._resonance_subsets(wrapper,state,permutation,groups)
                if not selector.alternatives:
                    if compilation is not None:
                        stats=compilation[0].timings;stats['resonance_rows_skipped']=stats.get('resonance_rows_skipped',0)+1
                    continue
                values=delta.evaluate([(state,np.asarray([permutation]))],coeffs,freqs,
                    mode_subset_selector=selector,**numerical)
                output[row][:,column]+=np.asarray(values[0]).reshape(len(coeffs))
                if compilation is not None:
                    stats=compilation[0].timings;stats['resonance_rows']=stats.get('resonance_rows',0)+1
        return output

    @staticmethod
    def _evaluate_action_coefficients(wrapper, state_perms, coeffs, freqs, *,
                                    shift=0., pruning_signature=b'', max_degree=None, **numerical):
        """Evaluate exact change-indexed occupation coefficients.

        Symbolic projection never binds occupations. At evaluation time only
        monomials nonzero on the requested occupations need to be contracted.
        Unresolved native change domains remain in each coefficient kernel.
        """
        from .Analytic import (PTTensorCoeffProductDAG, SqrtChangePoly,
                               PerturbationTheoryExpressionEvaluator as Expression)
        expression = wrapper.expr
        if not isinstance(expression,SqrtChangePoly) or np.any(expression.shift_start):
            raise ValueError('action coefficient evaluation requires an unshifted SqrtChangePoly')
        change = tuple(int(x) for x in expression.poly_change)
        fixed = wrapper.num_fixed
        if len(change) != fixed:
            raise ValueError('action coefficient change must match its external mode indices')
        degree = max(0,(_polynomial_degree(expression.poly_obj)-sum(map(abs,change)))//2)
        if max_degree is not None:degree=min(degree,max_degree)
        has_domains = (isinstance(expression.poly_obj,PTTensorCoeffProductDAG)
                       and expression.poly_obj.has_change_domains)
        active = _current_evaluation_plan.get()
        root = None if active is None else active[1].root
        cache = BoundedCache.sampled(10000,64*1024**2) if root is None else root.cache(
            'operator_action_coefficients',64*1024**2)
        force_tree = coeffs.tree if hasattr(coeffs,'tree') else coeffs
        force_key = None if root is None else root.bind(force_tree)
        frequency_key = None if root is None else root.bind(np.asarray(freqs))
        filter_key = _fingerprint({name:numerical.get(name) for name in (
            'degenerate_changes','only_degenerate_terms','projection_bindings','zero_cutoff')})
        compilation = current_compilation()
        compiled = compilation is not None and numerical.get('evaluation_mode') in {None,'compiled'}
        groups = collections.OrderedDict()
        output = [np.zeros((len(coeffs),len(perms))) for _,perms in state_perms]
        token = _current_evaluation_mode.set(None)
        try:
            for row,(state,permutations) in enumerate(state_perms):
                state = np.asarray(state)
                for column,permutation in enumerate(permutations):
                    permutation = np.asarray(permutation,dtype=int)
                    occupations = state[permutation]
                    actions = occupations-shift
                    occupied = np.flatnonzero(actions)
                    ladder = 1.
                    for occupation,delta in zip(occupations,change):
                        offsets = np.arange(delta+1,1) if delta < 0 else np.arange(1,delta+1)
                        ladder *= np.prod(occupation+offsets,dtype=float)
                    if ladder == 0:
                        continue
                    ladder = np.sqrt(ladder)
                    for power in range(degree+1):
                        for monomial in itertools.combinations_with_replacement(occupied,power):
                            counts = collections.Counter(monomial)
                            fixed_powers = tuple(counts.get(i,0) for i in range(fixed))
                            anchors = tuple(sorted(i for i in counts if i >= fixed))
                            powers = tuple(counts[i] for i in anchors)
                            pattern = fixed_powers,powers
                            cache_key = (expression.poly_obj,change,pattern,tuple(permutation),anchors,
                                         force_key,frequency_key,filter_key,
                                         tuple(state) if has_domains else None)
                            value = cache.get(cache_key,BoundedCache._missing)
                            weight = ladder*math.prod(actions[i]**p for i,p in counts.items())
                            if value is not BoundedCache._missing:
                                output[row][:,column] += value*weight
                                if compiled:compilation[0].timings['operator_coefficient_cache_hits'] += 1
                                continue
                            pending = groups.setdefault((pattern,anchors),collections.OrderedDict())
                            if cache_key not in pending:
                                carrier = state if has_domains else np.zeros(len(state),dtype=int)
                                pending[cache_key] = (carrier,permutation,[])
                            elif compiled:
                                compilation[0].timings['operator_coefficient_cache_hits'] += 1
                            pending[cache_key][2].append((row,column,weight))

            # Contract a whole physical coefficient table per symbolic pattern.
            # Batching avoids a numerical task-plan dispatch for every cell.
            kernels = {}
            for (pattern,anchors),pending in groups.items():
                fixed_powers,powers = pattern
                if pattern not in kernels:
                    projected = expression.project_coefficients(powers,shift,
                        fixed_powers=fixed_powers,preserve_domains=True)
                    coefficient = projected.mutate(poly_change=(0,)*fixed)
                    kernel = Expression(wrapper.op,coefficient,wrapper.change,logger=wrapper.logger)
                    if compiled:
                        store,compile_options,_,_ = compilation
                        compile_options = dict(compile_options)
                        if not store.using_artifacts:
                            compile_options.setdefault('mode_dimension',len(freqs))
                        projection = _OccupationCoefficientProjection(powers,shift,fixed_powers,True)
                        kernel = store.compile(kernel,projection=projection,
                            context=dict(operator=wrapper.op if wrapper.op is not None else expression.poly_obj,
                                change=change,fixed_powers=fixed_powers,powers=powers,shift=shift,
                                excluded_resonances=pruning_signature.hex()),
                            active_operators=None if store.using_artifacts else _active_coefficient_operators(coeffs),**compile_options)
                    kernels[pattern] = kernel
                carriers = collections.OrderedDict()
                for key,(carrier,permutation,targets) in pending.items():
                    carriers.setdefault(tuple(carrier),[]).append((key,permutation,targets))
                requests = [(np.asarray(state,dtype=int),np.asarray([row[1] for row in rows]))
                            for state,rows in carriers.items()]
                if compiled:
                    _,_,task_options,default_parallelizer = compilation
                    options = {k:v for k,v in numerical.items() if v is not None and k not in {
                        'parallelizer','evaluation_mode','dag_stream_options','use_shared_memory'}}
                    options.update(task_options)
                    values = kernels[pattern].evaluate(requests,coeffs,freqs,anchored_modes=anchors,
                        parallelizer=numerical.get('parallelizer') or default_parallelizer,
                        use_shared_memory=numerical.get('use_shared_memory'),**options)
                    compilation[0].timings['operator_coefficient_contractions'] += len(pending)
                    compilation[0].timings['operator_coefficient_batches'] += 1
                else:
                    options = dict(numerical)
                    options['evaluation_mode'] = options.get('evaluation_mode') or 'dag'
                    values = kernels[pattern].evaluate(requests,coeffs,freqs,
                        mode_subset_selector=_AnchoredModeSubset(anchors),**options)
                for rows,block in zip(carriers.values(),values):
                    block = np.asarray(block).reshape(len(coeffs),len(rows))
                    for index,(key,_,targets) in enumerate(rows):
                        value = block[:,index].copy()
                        value.setflags(write=False)
                        cache[key] = value
                        for row,column,weight in targets:
                            output[row][:,column] += value*weight
        finally:
            _current_evaluation_mode.reset(token)
        return output

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
        numerical.pop('evaluation_mode',None)
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
    denominator_program: object = None


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
    structural_patterns: object = dataclasses.field(default=None, repr=False, compare=False)
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

    def _bind_request_pattern_map(self, patterns):
        """Select a compiled request table for the concrete active tensor group.

        Which requests survive depends on actual tensors, aliases and cutoffs.
        The abstract patterns and their IDs were already compiled; this step
        instantiates only the selected physical group, without pattern analysis.
        """
        pool, identifiers = self.structural_patterns[len(patterns[0])]
        selected, inverse, positions = [], [], {}
        for pattern in patterns:
            identifier = identifiers[tuple(pattern)]
            if identifier not in positions:
                positions[identifier] = len(selected)
                selected.append(identifier)
            inverse.append(positions[identifier])
        return pool[np.asarray(selected, dtype=np.intp)], np.asarray(inverse, dtype=np.intp)

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
                    if self.structural_patterns is not None:
                        unique_structural, inverse = self._bind_request_pattern_map(structural_patterns)
                        unique_patterns = self._remap_request_patterns(unique_structural, subset, num_fixed)
                    elif self.use_cached_request_patterns:
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
        execution = getattr(expression, '_compiled_execution', None)
        if execution is not None:
            name = 'reduced' if reduced else 'native'
            if name not in execution['plans']:
                raise ValueError('compiled expression lacks requested structural plan')
            if self.stats is not None: self.stats.plan_cache_hits += 1
            return execution['plans'][name]
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
        execution = getattr(getattr(getattr(self, 'plan', None), 'expression', None), '_compiled_execution', None)
        patterns = None if execution is None else execution['coefficient_patterns']
        if self.preparation_context is not None:
            def create():
                return IndexedCoefficientGatherPlan(coefficient_specs, coefficient_tensors,
                    0 if not coefficient_tensors else len(coefficient_tensors[0]),
                    preparation_context=self.preparation_context, structural_patterns=patterns)
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
            expansion_count,
            structural_patterns=patterns
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
                            preparation_context=self.preparation_context,
                            programs=([terms[index].denominator_program for index in indices]
                                      if terms[indices[0]].denominator_program is not None else None))
                        return mask, np.where(mask, factors, 1.)
                    changes = np.stack([terms[index].denominator_program['energy_rows']
                        if terms[index].denominator_program is not None else np.asarray(terms[index].energy_changes)
                        for index in indices])
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
                    and all(term.denominator_program['pure_projectors']
                            if term.denominator_program is not None else
                            not any(getattr(row[0], '_conditional_endpoint', False)
                                    or getattr(row[0], '_filtered_change_domain', False)
                                    for row in term.energy_changes) for term in terms))
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
