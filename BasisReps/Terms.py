"""
Provides an abstract Hamiltonian object that can be used when building representations
"""

__all__ = [
    "TermComputer",
    "ExpansionTerm"
]

import numpy as np, itertools as ip, scipy.sparse as sp

from McUtils.Numputils import SparseArray

from .Operators import Operator

#TODO: add in some level of support for caching
class TermComputer:
    """
    A TermComputer provides a simple interface to compute only some elements of high-dimensional tensors.
    It takes a tensor shape and a function to compute tensor elements.
    The `compute` function should be able to take a block of indices and return all the matrix elements.
    """

    def __init__(self, compute, n_quanta):
        """
        :param compute: the function that turns indices into values
        :type compute: callable | Operator
        :param n_quanta: the total quanta used in the representations (necessary for shape reasons)
        :type n_quanta: tuple[int]
        """
        if isinstance(compute, Operator):
            operator = compute
            compute = self._compute_op_els
        else:
            operator = None
        self.operator = operator
        self.compute = compute
        self.dims = n_quanta
    def _compute_op_els(self, inds):
        return self.operator[inds] #compute: c[inds]

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
        # the block pulling can be quite a bit faster, so we try to detect first off if we want to do that
        pull_elements = True
        # first we check if we've got something like from `np.ix_(n, m)`
        if isinstance(n, np.ndarray) and isinstance(m, np.ndarray):
            if len(n.shape) > 1 and len(m.shape) > 1:
                pull_elements = False
        if pull_elements:
            # next we check to see if idx is really just a single element
            pull_elements = all(isinstance(x, (int, np.integer)) for x in idx)
            if not pull_elements:
                # if not a single element, we make sure there are no slices
                pull_elements = all(not isinstance(x, (int, np.integer, slice)) for x in idx)
                # print(">?>", pull_elements, idx)
                if pull_elements:
                    e1 = len(idx[0])
                    pull_elements = all(len(x) == e1 for x in idx)
                    # print(">??", pull_elements)

        # We figure out the row spec
        if not isinstance(n, int):
            if isinstance(n, np.ndarray):
                n = n.flatten()
            if not isinstance(n, slice):
                n = np.array(n, dtype=int)
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
        if isinstance(els, int) and els == 0:
            # short-circuited :|
            return els
        elif not pull_elements:

            shp = (len(np.unique(blocks[:, 0])), len(np.unique(blocks[:, 1])))
            extra_shp = els.shape[:-1] # some terms will return higher-dimensional results?
            # for sparse arrays this happens in-place :|
            els = els.reshape(extra_shp + shp).squeeze()

        return els

    def __getitem__(self, item):
        if not isinstance(item, tuple):
            item = (item,)
        if len(item) == 1:
            item = item + (slice(None, None, None),)
        if len(item) > 2:
            raise Exception("index spec '{}' must be of dimension 2".format(item))
        return self.get_element(*item)

    def __rmul__(self, other):
        if isinstance(other, (int, float, np.integer, np.floating)):
            return ExpansionTerm([other], [self], self.dims)
        else:
            raise TypeError("operator * not defined for objects of type {0} and {1} (only numbers are supported with {0})".format(
                type(self).__name__,
                type(other).__name__
            ))
    def __mul__(self, other):
        if isinstance(other, (int, float, np.integer, np.floating)):
            return ExpansionTerm([other], [self], self.dims)
        else:
            raise TypeError("operator * not defined for objects of type {0} and {1} (only numbers are supported with {0})".format(
                type(self).__name__,
                type(other).__name__
            ))
    def __add__(self, other):
        if isinstance(other, TermComputer):
            if other.dims != self.dims:
                if isinstance(other, ExpansionTerm):
                    return other + self
                raise ValueError("Can't combine TermComputer objects with dim {} and {}".format(
                    self.dims,
                    other.dims
                ))
            return ExpansionTerm([1, 1], [self, other], self.dims)
        else:
            raise TypeError("operator * not defined for objects of type {0} and {1} (only subclasses of TermComputer are supported with {0})".format(
                type(self).__name__,
                type(other).__name__
            ))

    def __repr__(self):
        return "{}(<{}>, {})".format(
            type(self).__name__,
            ", ".join(str(s) for s in self.dims) if hasattr(self, 'dims') else '???',
            self.operator if self.operator is not None else self.compute
        )

class ExpansionTerm(TermComputer):
    """
    Provides support for terms that look like `1/2 pGp + 1/2 dV/dQdQ QQ` by computing each term on its own
    """
    def __init__(self, coeffs, computers, n_quanta):
        """
        :param coeffs: The expansion coefficients
        :type coeffs: Iterable[float]
        :param compute: the functions that turns indices into values
        :type compute: Iterable[callable | Operator]
        :param n_quanta: the total quanta used in the representations (necessary for shape reasons)
        :type n_quanta: tuple[int]
        """
        self.coeffs = np.array(coeffs)
        self.computers = [TermComputer(c, n_quanta) if not isinstance(c, TermComputer) else c for c in computers]
        super().__init__(None, n_quanta)

    def __rmul__(self, other):
        if isinstance(other, (int, float, np.integer, np.floating)):
            return type(self)(self.coeffs * other, self.computers, self.dims)
        else:
            raise TypeError(
                "operator * not defined for objects of type {0} and {1} (only numbers are supported with {0})".format(
                    type(self).__name__,
                    type(other).__name__
                ))
    def __mul__(self, other):
        if isinstance(other, (int, float, np.integer, np.floating)):
            return type(self)(self.coeffs*other, self.computers, self.dims)
        else:
            raise TypeError(
                "operator * not defined for objects of type {0} and {1} (only numbers are supported with {0})".format(
                    type(self).__name__,
                    type(other).__name__
                ))
    def __add__(self, other):
        if isinstance(other, TermComputer):
            if other.dims != self.dims:
                raise ValueError("Can't combine TermComputer objects with dim {} and {}".format(
                    self.dims,
                    other.dims
                ))
            if isinstance(other, ExpansionTerm):
                return type(self)(
                    np.concatenate([self.coeffs, other.coeffs]),
                    self.computers + other.computers,
                    self.dims
                )
            else:
                return type(self)(
                    np.concatenate([self.coeffs, [1]]),
                    self.computers + [other],
                    self.dims
                )
        else:
            raise TypeError(
                "operator * not defined for objects of type {0} and {1} (only subclasses of TermComputer are supported with {0})".format(
                    type(self).__name__,
                    type(other).__name__
                ))
    def get_element(self, n, m):

        # try:
        els = None
        for c, t in zip(self.coeffs, self.computers):
            if not (isinstance(c, (int, float, np.integer, np.floating)) and c == 0):
                bits = t.get_element(n, m)
                scaled = bits * c
                if isinstance(scaled, (int, float, np.integer, np.floating)) and scaled != 0:
                    # raise Exception(bits, c, scaled, n,m, t)
                    raise ValueError(" ".join([
                        "Adding a constant ({}) to a sparse operator ({}) would cast to dense.",
                        "Error likely occurred in getting elements for {}.",
                        "Explicitly subclass {} if you truly need the constant shift.",
                    ]).format(
                        scaled,
                        els,
                        t,
                        type(self).__name__
                    ))
                else:
                    if els is None:
                        els = scaled
                    else:
                        if isinstance(scaled, (SparseArray, sp.spmatrix)):
                            scaled = scaled.toarray()
                            # import McUtils.Plots as plt
                            #
                            # plt.ArrayPlot(scaled).show()
                        # print(scaled.shape, els.shape)
                        els += scaled

        return els

    def __repr__(self):
        return "{}(<{}>, ({}), ({}))".format(
            type(self).__name__,
            ", ".join(str(s) for s in self.dims) if hasattr(self, 'dims') else '???',
            self.coeffs,
            self.computers
        )