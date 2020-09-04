"""
Provides an abstract Hamiltonian object that can be used when building representations
"""

__all__ = [
    "TermComputer"
]

import numpy as np, itertools as ip
from .Operators import Operator

#TODO: add in some level of support for caching
class TermComputer:
    """
    A TermComputer provides a simple interface to compute only some elements of high-dimensional tensors.
    It takes a tensor shape and a function to compute tensor elements.
    The `compute` function should be able to take a block of indices and return all the matrix elements.
    """
    def __init__(self, compute, n_quanta):
        if isinstance(compute, Operator):
            operator = compute
            compute = lambda inds, c=compute: c[inds]
        else:
            operator = None
        self.operator = operator
        self.compute = compute
        self.dims = n_quanta
    @property
    def diag(self):
        ndims = int(np.prod(self.dims))
        return self[np.arange(ndims), np.arange(ndims)]
    def get_element(self, n, m):
        """
        Computes term elements.
        Determines first whether it needs to pull single elements or blocks of them.

        :param n:
        :type n:
        :param m:
        :type m:
        :return:
        :rtype:
        """

        dims = self.dims
        ndims = int(np.prod(dims))
        idx = (n, m)

        # There are two possible modes for this pulling individual elements or pulling blocks
        # the block pulling is quite a bit faster, so we try to detect first off if we want to do that
        pull_elements = True
        if isinstance(n, np.ndarray) and isinstance(m, np.ndarray):
            if len(n.shape) > 1 and len(m.shape) > 1:
                pull_elements = False
        if pull_elements:
            pull_elements = all(isinstance(x, (int, np.integer)) for x in idx)
            if not pull_elements:
                pull_elements = all(not isinstance(x, (int, np.integer, slice)) for x in idx)
                if pull_elements:
                    e1 = len(idx[0])
                    pull_elements = all(len(x) == e1 for x in idx)
        # We figure out the row spec
        if not isinstance(n, int):
            if isinstance(n, np.ndarray):
                n = n.flatten()
            if not isinstance(n, slice):
                n = np.array(n)
            n = np.arange(ndims)[n]
        else:
            n = [n]
        # Then the column spec
        if not isinstance(m, int):
            if isinstance(m, np.ndarray):
                m = m.flatten()
            if not isinstance(m, slice):
                m = np.array(m)
            m = np.arange(ndims)[m]
        else:
            m = [m]

        if pull_elements:
            # If we're just pulling elements we need only unravel those indices
            n = np.unravel_index(n, dims)
            m = np.unravel_index(m, dims)
        else:
            # If we're pulling blocks we need to compute the product of the row
            #  and column indices to get the total index spec
            blocks = np.array(list(ip.product(n, m)))
            n = np.unravel_index(blocks[:, 0], dims)
            m = np.unravel_index(blocks[:, 1], dims)

        # we define a temporary helper to pad repeat the indices if necessary
        def pad_lens(a, b):
            if isinstance(a, (int, np.integer)) and not isinstance(b, (int, np.integer)):
                a = np.full((len(b),), a)
            if isinstance(b, (int, np.integer)) and not isinstance(a, (int, np.integer)):
                b = np.full((len(a),), b)
            return a, b
        i = tuple(pad_lens(a, b) for a, b in zip(n, m))
        els = self.compute(i)
        if not pull_elements:
            shp = (len(np.unique(blocks[:, 0])), len(np.unique(blocks[:, 1])))
            # for sparse arrays this something happens in-place :|
            els = els.reshape(shp).squeeze()
        return els

    def __getitem__(self, item):
        if not isinstance(item, tuple):
            item = (item,)
        if len(item) == 1:
            item = item + (slice(None, None, None),)
        if len(item) > 2:
            raise Exception("index spec '{}' must be of dimension 2".format(item))
        return self.get_element(*item)