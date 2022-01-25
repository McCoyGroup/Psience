
import numpy as np

__all__ = [
    "PerturbationTheoryException",
    "Settings"
]
class PerturbationTheoryException(Exception):
    pass

class Settings:
    non_zero_cutoff = 1.0e-14

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