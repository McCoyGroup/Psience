
import sys, collections, contextlib
import numpy as np
from McUtils.Numputils import SparseArray

__all__ = [
    "PerturbationTheoryException",
    "Settings",
    "Diagnostics",
    "BoundedCache"
]
class PerturbationTheoryException(Exception):
    pass

class Settings:
    non_zero_cutoff = 1.0e-14

class Diagnostics:
    """
    Switch for the evaluation statistics the analytic evaluators can record
    (indexed-evaluator counters, DAG/cache statistics, compilation summaries).
    Off by default: normal runs do no bookkeeping and the `get_last_*_stats`
    accessors return `None`. Enable with `Diagnostics.enabled = True` or
    `with Diagnostics.collect(): ...` when profiling or testing.
    """
    enabled = False

    @classmethod
    @contextlib.contextmanager
    def collect(cls, enabled=True):
        prev = cls.enabled
        cls.enabled = enabled
        try:
            yield cls
        finally:
            cls.enabled = prev

def _safe_dot(a, b):
    # generalizes the dot product so that we can use 0 as a special value...
    if (
            isinstance(a, (int, np.integer, float, np.floating)) and a == 0
            or isinstance(b, (int, np.integer, float, np.floating)) and b == 0
    ):
        return 0

    if isinstance(a, np.ndarray):
        doots = np.dot(a, b)
    else:
        doots = a.dot(b)

    if isinstance(b, np.ndarray) and isinstance(doots, SparseArray):
        doots = doots.asarray()

    return doots

class BoundedCache:
    """
    Item- and (optionally) byte-bounded mapping with FIFO (default) or LRU eviction.
    This is the one cache type behind the analytic evaluators' batch, state-permutation,
    plan and request-pattern caches.

    Entry sizes come from `sizer(key, value)`; the default charges only the
    ndarray payload of the value (directly or one level inside a tuple/list)
    plus a fixed per-entry overhead, so accounting never walks symbolic keys.
    Hit/miss/eviction counts are kept only when `Diagnostics.enabled`.
    """
    _missing = object()
    entry_overhead = 128

    class View:
        """A namespaced window onto a shared cache."""
        def __init__(self, parent, namespace):
            self.parent = parent
            self.namespace = namespace
        def get(self, key, default=None):
            return self.parent.get((self.namespace, key), default)
        def __contains__(self, key):
            return (self.namespace, key) in self.parent
        def __getitem__(self, key):
            return self.parent[(self.namespace, key)]
        def __setitem__(self, key, value):
            self.parent[(self.namespace, key)] = value

    def __init__(self, max_items=None, max_bytes=None, sizer=None, lru=False):
        self.max_items = None if max_items is None else max(0, int(max_items))
        self.max_bytes = None if max_bytes is None else max(0, int(max_bytes))
        self.sizer = self.payload_size if sizer is None else sizer
        self.lru = lru
        self._data = collections.OrderedDict()
        self._sizes = {}
        self.nbytes = 0
        self.hits = self.misses = self.evictions = 0
        self.peak_items = self.peak_bytes = 0

    @classmethod
    def payload_size(cls, key, value):
        if isinstance(value, np.ndarray):
            return cls.entry_overhead + value.nbytes
        if isinstance(value, (tuple, list)):
            return cls.entry_overhead + sum(v.nbytes for v in value if isinstance(v, np.ndarray))
        return cls.entry_overhead

    @classmethod
    def sampled_size(cls, key, value):
        """Sizer that also charges keys: arrays/bytes in full, containers by a
        bounded sample (16 items, depth 4), objects with `memory_size()` by that."""
        return cls._sampled_size(key) + cls._sampled_size(value)

    @classmethod
    def _sampled_size(cls, value, depth=0):
        memory_size = getattr(value, 'memory_size', None)
        if memory_size is not None and not isinstance(value, type):
            return memory_size()
        if isinstance(value, np.ndarray):
            return sys.getsizeof(value) + value.nbytes
        if isinstance(value, (bytes, bytearray, memoryview)):
            return sys.getsizeof(value)
        if isinstance(value, (tuple, list)) and depth < 4:
            return sys.getsizeof(value) + sum(
                max(0, cls._sampled_size(item, depth=depth + 1) - sys.getsizeof(item))
                for item in value[:16]
            )
        return sys.getsizeof(value)

    @classmethod
    def sampled(cls, max_items=100000, max_bytes=64 * 1024 ** 2):
        """The analytic evaluators' numeric cache: FIFO, keys and values charged by `sampled_size`."""
        return cls(max_items, max_bytes, sizer=cls.sampled_size)

    def view(self, namespace):
        return self.View(self, namespace)

    def __len__(self):
        return len(self._data)
    def __contains__(self, key):
        return key in self._data
    def keys(self):
        return self._data.keys()

    def get(self, key, default=None):
        value = self._data.get(key, self._missing)
        if value is self._missing:
            if Diagnostics.enabled: self.misses += 1
            return default
        if self.lru: self._data.move_to_end(key)
        if Diagnostics.enabled: self.hits += 1
        return value

    def __getitem__(self, key):
        value = self.get(key, self._missing)
        if value is self._missing:
            raise KeyError(key)
        return value

    def __setitem__(self, key, value):
        if self.max_items == 0 or self.max_bytes == 0:
            return
        if key in self._data:
            self.pop(key)
        size = self.sizer(key, value) if self.max_bytes is not None else 0
        if self.max_bytes is not None and size > self.max_bytes:
            return
        self._data[key] = value
        self._sizes[key] = size
        self.nbytes += size
        while ((self.max_items is not None and len(self._data) > self.max_items)
               or (self.max_bytes is not None and self.nbytes > self.max_bytes)):
            self.popitem()
            if Diagnostics.enabled: self.evictions += 1
        if Diagnostics.enabled:
            self.peak_items = max(self.peak_items, len(self._data))
            self.peak_bytes = max(self.peak_bytes, self.nbytes)

    def pop(self, key, default=None):
        value = self._data.pop(key, self._missing)
        if value is self._missing:
            return default
        self.nbytes -= self._sizes.pop(key)
        return value

    def popitem(self):
        """Remove and return the oldest (FIFO) / least recently used (LRU) entry."""
        key, value = self._data.popitem(last=False)
        self.nbytes -= self._sizes.pop(key)
        return key, value

    def clear(self):
        self._data.clear(); self._sizes.clear(); self.nbytes = 0

    def stats(self):
        return {
            'cache_items': len(self._data), 'cache_bytes': self.nbytes,
            'cache_peak_items': self.peak_items, 'cache_peak_bytes': self.peak_bytes,
            'cache_hits': self.hits, 'cache_misses': self.misses, 'cache_evictions': self.evictions,
            'cache_max_items': self.max_items, 'cache_max_bytes': self.max_bytes
        }
