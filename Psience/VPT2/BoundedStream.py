"""Opt-in storage policy for the native analytic tensor-DAG evaluator.

This module does not change any symbolic products, conditions, or denominators.
Byte thresholds are checked periodically and describe the reachable Python
payload, not a hard bound on process RSS or a kernel's transient workspace.
"""
import collections
import contextlib
import os
import sys
import time
import numpy as np


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
        from .Analytic import (_BoundedEvaluationCache,
                               PerturbationTheoryExpressionEvaluator as Expression)
        root = self.root
        if name not in root.numeric:
            items = root.options.get('path_cache_size',10000)
            if items is None:
                items = Expression.default_materialized_path_cache_size
            root.numeric[name] = _BoundedEvaluationCache(items,max_bytes)
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
            self.data_hits += 1
            return self.data[key]
        self.data_misses += 1
        return self._store(self.data, key,
                           evaluator._prep_coeff_data(coefficient_sets, self.coefficients))

    def coefficient_plan(self, specs, tensors, create):
        key = (id(specs), id(tensors))
        entry = self.plans.get(key)
        if entry is not None:
            self.plan_hits += 1
            return entry[2], True
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
