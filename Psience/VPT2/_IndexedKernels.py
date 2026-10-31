"""Optional compiled ordinary harmonic denominator contraction."""
import functools
import numpy as np

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
