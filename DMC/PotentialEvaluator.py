import numpy as np, multiprocessing as mp
from McUtils.CPotentialLib import CPotential

__all__ = ["PotentialEvaluator"]

##################################################################################################################
##
##                                                  PotentialEvaluator
##
##

class PotentialEvaluator:

    def __init__(self, func, mode="single"):

        self.f = func
        if type(func).__name__ == "PyCapsule": #poor man's type check...
            self.f = CPotential(func, mode = mode)
        self.mode = mode
        self._pool = None

    def _call_serial(self, atoms, coords):
        return np.array([ self.f(atoms, coord) for coord in coords ])
    def _call_parallel(self, atoms, coords, map_switch_bytes = 50000, chunk_size = 100):
        from sys import getsizeof as size
        if self._pool is None:
            self._pool = mp.Pool()
        if size(coords[0])<map_switch_bytes:
            res = self._pool.map(self.f, coords, chunk_size)
        else:
            res = self._pool.imap(self.f, coords, chunk_size)
        return np.array(res)

    def call(self, atoms, coords):
        if self.mode == "single":
            return self._call_serial(atoms, coords)
        elif self.mode == "parallel":
            return self._call_parallel(atoms, coords)
        else:
            return self.f(atoms, coords)

    def __call__(self, *args, **kwargs):
        return self.call(*args, **kwargs)