"""Opt-in storage policy for the native analytic tensor-DAG evaluator.

This module does not change any symbolic products, conditions, or denominators.
Byte thresholds are checked periodically and describe the reachable Python
payload, not a hard bound on process RSS or a kernel's transient workspace.
"""
import collections
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
    seen = set()
    size = 0
    while todo:
        obj = todo.pop()
        identity = id(obj)
        if identity in seen:
            continue
        seen.add(identity)
        size += sys.getsizeof(obj)
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
        elif type(obj).__module__.startswith('Psience.VPT2.'):
            attrs = getattr(obj, '__dict__', None)
            if attrs is not None:
                todo.append(attrs)
            # Some immutable condition records use slots rather than __dict__.
            for base in type(obj).__mro__:
                slots = getattr(base, '__slots__', ())
                if isinstance(slots, str):
                    slots = (slots,)
                for name in slots:
                    if name not in {'__weakref__', '__dict__'}:
                        value = getattr(obj, name, None)
                        if value is not None:
                            todo.append(value)
    return size


class StreamPreparationContext:
    """Reuse coefficient preparation within one call, without retaining batches.

    The coefficient owner fixes this context's lifetime. It is never reused
    across calls or coefficient updates. Workers receive their own context.
    Expression plans remain batch-local; coefficient routing plans may persist.
    """
    def __init__(self, coefficients, max_items=32):
        self.coefficients = coefficients
        self.max_items = max(0, int(max_items))
        self.data = collections.OrderedDict()
        self.plans = collections.OrderedDict()
        self.data_hits = self.data_misses = 0
        self.plan_hits = self.plan_misses = 0

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
        if entry is not None and entry[0] is specs and entry[1] is tensors:
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
            self.bytes = retained_size(pending)
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
                    stream_peak_batch_terms=self.maximum_terms,
                    stream_peak_measured_pending_bytes=self.maximum_bytes,
                    stream_accounting_seconds=self.scan_seconds,
                    stream_flush_reasons=dict(self.flush_reasons))
