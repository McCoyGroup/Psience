"""Extracted from AnalyticModelsTests.test_RedundantG via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest AnalyticModelsTests.test_RedundantG"""

from Peeves.TestUtils import *
from Psience.AnalyticModels import *
from McUtils.Plots import *
import McUtils.Numputils as nput
from unittest import TestCase
import sys, h5py, math, numpy as np

class AnalyticModelsTests(TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        np.set_printoptions(linewidth=int(100000000.0))
    maxDiff = None

    def check_expr(self, test, real, raise_error=True):
        if not isinstance(test, int):
            test = test.expand().simplify().expand()
        if not isinstance(real, int):
            real = real.expand().simplify().expand()
        diff = test - real
        if not isinstance(diff, int):
            diff = diff.expand().simplify().expand()
        msg = '\nTest: {}\nReal: {}\nDiff: {}'.format(test, real, diff)
        if not raise_error:
            return msg
        else:
            self.assertEquals(diff, 0, msg=msg)

    @debugTest
    def test_RedundantG(self):
        from Psience.Molecools import Molecule
        from Psience.Psience.AnalyticModels.AnalyticJacobianDotCalculator import InternalJacobianDisplacements
        from Psience.AnalyticModels import AnalyticKineticEnergyConstructor, GeometricFunction
        file_name = 'nh3.fchk'
        test_internals = [[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1], [3, 0, 1, 2]]
        mol = Molecule.from_file(TestManager.test_data(file_name), internals=test_internals)
        coords = [(0, 1), (2, 0), (0, 3), (2, 0, 1), (3, 0, 1), (2, 0, 3), (2, 1, 0), (3, 1, 0), (3, 1, 2), (1, 2), (1, 3), (2, 3), (3, 0, 1, 2), (1, 0, 2, 3), (2, 0, 3, 1), (0, 1, 3, 2)]
        B = np.array([nput.dist_vec(mol.coords, *inds) if len(inds) == 2 else nput.angle_vec(mol.coords, inds[1], inds[0], inds[2]) if len(inds) == 3 else nput.dihed_vec(mol.coords, *inds) if len(inds) == 4 else nput.book_vec(mol.coords, *inds) if len(inds) == 5 and inds[-1] == -1 else None for inds in coords])
        M = np.diag(np.repeat(1 / np.sqrt(mol.atomic_masses)[:, np.newaxis], 3))
        raise Exception(AnalyticKineticEnergyConstructor.g((1, 3), (1, 2, 3), method='direct'))
        B = B @ M
        g_direct = B @ B.T
        g_red = AnalyticKineticEnergyConstructor.g_matrix(coords, return_function=True, method='direct')(mol.atomic_masses, mol.coords)
        print((g_direct,))
        print((g_red,))
        print(np.round(g_red - g_direct, 8))
        print(np.linalg.eigvalsh(g_direct))
        print(np.linalg.eigvalsh(g_red))
        raise Exception(...)
