from Peeves.TestUtils import *
from unittest import TestCase
from Psience.BasisReps import *
import sys, os, numpy as np

class BasisSetTests(TestCase):

    @debugTest
    def test_HOBasis(self):

        n = 7
        basis = HarmonicOscillatorBasis(n)

        xx = basis.representation('x', 'x')[:, :].todense()
        targ =np.zeros((n, n))
        targ[np.arange(n), np.arange(n)] = np.arange(n) + 1/2
        targ[np.arange(n-2), np.arange(2, n)] = np.sqrt(np.arange(1, n-1)*(np.arange(1, n-1)+1)/4)
        targ[np.arange(2, n), np.arange(n-2)] = np.sqrt(np.arange(1, n-1)*(np.arange(1, n-1)+1)/4)

        # raise Exception([
        #     targ**2,
        #     xx**2
        #     ])

        self.assertLess(np.average(np.abs(xx - targ)), 1e-14)