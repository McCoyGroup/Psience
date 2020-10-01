"""
Provides representations based on starting from a harmonic oscillator basis
"""

__all__ = [
    "HarmonicOscillatorBasis"
]

import scipy.sparse as sp, numpy as np, functools
from .Bases import *
from McUtils.Data import WavefunctionData

class HarmonicOscillatorBasis(RepresentationBasis):
    """
    Provides a concrete implementation of RepresentationBasis using the H.O.
    Need to make it handle units a bit better.
    Currently 1D, need to make multidimensional in the future.
    """
    name='HarmonicOscillator'
    def __init__(self, n_quanta, m=None, re=None, dimensionless=True):
        self.dimensionless = dimensionless
        if not dimensionless and any(x is None for x in (m, re)):
            raise ValueError("if not dimensionless, parameters cannot be 'None'")
        self.data = WavefunctionData['HarmonicOscillator']
        super().__init__(self.data['Wavefunction'], n_quanta)

    @functools.lru_cache(maxsize=128)
    def p(self, n):
        mat = self.pmatrix_ho(n)
        if not self.dimensionless:
            raise NotImplementedError
        return mat
    @staticmethod
    def pmatrix_ho(n):
        """
        There's one big subtlety to what we're doing here, which is that
          for efficiency reasons we return an entirely real matrix
        The reason for that is we assumed people would mostly use it in the context
          of stuff like pp, pQp, or pQQp, in which case the imaginary part pulls out
          and becomes a negative sign
        We actually use this assumption across _all_ of our representations
        :param n:
        :type n:
        :return:
        :rtype: sp.csr_matrix
        """
        # the imaginary terms pull out and just become a negative sign
        ar = 1 / np.sqrt(2) * np.sqrt(np.arange(1, n))
        bands = [
            [ ar,  1],
            [-ar, -1]
        ]
        return sp.csr_matrix(sp.diags([b[0] for b in bands], [b[1] for b in bands]))

    @functools.lru_cache(maxsize=128)
    def x(self, n):
        return self.qmatrix_ho(n)
    @staticmethod
    def qmatrix_ho(n):
        """

        :param n:
        :type n:
        :return:
        :rtype: sp.csr_matrix
        """

        ar = 1 / np.sqrt(2) * np.sqrt(np.arange(1, n))
        bands = [
            [ar,  1],
            [ar, -1]
        ]
        return sp.csr_matrix(sp.diags([b[0] for b in bands], [b[1] for b in bands]))