
import contextlib
import numpy as np
from McUtils.Numputils import SparseArray

__all__ = [
    "PerturbationTheoryException",
    "Settings",
    "Diagnostics"
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