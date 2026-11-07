"""
An optional on-disk layer under the VPT2 exact-result cache (``PerturbationTheoryExpressionEvaluator._poly_cache``).

That cache maps ``(polynomial, change, baseline_shift)`` to a dict of ``{state index: value}``. The
persistent layer keys each of those dicts by the same expansion fingerprint `CompiledExpressions` uses
(definition + implementation signature), keeps the in-memory cache in front, loads a polynomial's values
the first time it is looked up, and writes the dicts that grew back to an `ObjectStore` on `flush()`.
Polynomials that can't be fingerprinted, and value dicts whose keys aren't uniform integer state indices,
are simply not persisted.
"""
import collections.abc

import numpy as np

__all__ = ["PersistentEvaluationCache"]

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
            from . import CompiledExpressions as ce
            fingerprint = fingerprint or ce._fingerprint
            signature = signature if signature is not None else ce._implementation_signature()
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
