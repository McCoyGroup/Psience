
__all__ = [
    "LegendreBasis",
    "ChebyshevBasis"
]

import scipy.sparse as sp, numpy as np
from .Bases import *

class LegendreBasis(RepresentationBasis):
    name='Legendre'
    def __init__(self, n_quanta):
        super().__init__(None, n_quanta)
    def __eq__(self, other):
        return isinstance(other, type(self)) and other.dimensions == self.dimensions

    selection_rules_mapping = {
        'x':[-1, 1],
        'p':[-1, 1],
        "I":[0]
    }

    def p(self, n):
        raise NotImplementedError("...")

    def p2(self, n):
        bands = [
            [(1 - np.arange(1, n+1)) * np.arange(1, n+1), 0]
        ]
        return sp.csr_matrix(sp.diags([b[0] for b in bands], [b[1] for b in bands]))

    def x(self, n):
        ar = np.arange(1, n)
        bands = [
            [np.sqrt(ar ** 2 / ((2 * ar - 1) * (2 * ar + 1))),  1],
            [np.sqrt(ar ** 2 / ((2 * ar - 1) * (2 * ar + 1))), -1]
        ]
        return sp.csr_matrix(sp.diags([b[0] for b in bands], [b[1] for b in bands]))

    def __repr__(self):
        return "LegendreBasis({})".format(
            ",".join(str(x) for x in self.dimensions)
        )

class ChebyshevBasis(RepresentationBasis):
    name='Chebyshev'
    def __init__(self, n_quanta):
        super().__init__(None, n_quanta)
    def __eq__(self, other):
        return isinstance(other, type(self)) and other.dimensions == self.dimensions

    selection_rules_mapping = {
        'x':[-1, 1],
        'p':[-1, 1],
        "I":[0]
    }

    def p(self, n):
        raise NotImplementedError("...")

    def p2(self, n):
        """
        0  -> (Which[# \[Equal] 1, 3, # \[Equal] 2, 1/2,
     True, (-1 - 8 (# - 1) - 4 (# - 1)^2)]*1/8 &),
    2  -> (If[# \[Equal] 1,
    15/(8*Sqrt[2]), (15 + 16 (# - 1) + 4 (# - 1)^2)*1/16] &),
    -2 -> (If[#2 \[Equal] 1,
    15/(8*Sqrt[2]), (15 + 16 (#2 - 1) + 4 (#2 - 1)^2)*1/16] &)
        :param n:
        :type n:
        :return:
        :rtype:
        """
        ar = np.arange(2, n)
        a = np.concatenate([[3, 1/2], (-1 - 8*ar - 4*ar**2)/8])
        ar = np.arange(1, n)
        b = np.concatenate([[15/(8*np.sqrt(2))], (15 + 16*ar + 4*ar**2)/16])
        bands = [
            [b, 2],
            [a, 0],
            [b, -2]
        ]
        return sp.csr_matrix(sp.diags([b[0] for b in bands], [b[1] for b in bands]))

    def x(self, n):
        ar = np.concatenate([[1/np.sqrt(2)], np.full(n-1, 1/2)])
        bands = [
            [ar,  1],
            [ar, -1]
        ]
        return sp.csr_matrix(sp.diags([b[0] for b in bands], [b[1] for b in bands]))

    def __repr__(self):
        return "ChebyshevBasis({})".format(
            ",".join(str(x) for x in self.dimensions)
        )

class LaguerreBasis(RepresentationBasis):
    name='Laguerre'
    def __init__(self, n_quanta):
        super().__init__(None, n_quanta)
    def __eq__(self, other):
        return isinstance(other, type(self)) and other.dimensions == self.dimensions

    selection_rules_mapping = {
        'x':[-1, 0, 1],
        'p':[-1, 0, 1],
        "I":[0]
    }

    def p(self, n):
        raise NotImplementedError("...")

    def x(self, n):
        bands = [
            [-np.arange(1, n),    1],
            [2*np.arange(n)+3,   0],
            [-np.arange(1, n),   -1]
        ]
        return sp.csr_matrix(sp.diags([b[0] for b in bands], [b[1] for b in bands]))

    def __repr__(self):
        return "LaguerreBasis({})".format(
            ",".join(str(x) for x in self.dimensions)
        )