
from Peeves.TestUtils import *
from Psience.Data import *
from McUtils.Coordinerds import cartesian_to_zmatrix
from McUtils.Plots import *
from unittest import TestCase
import sys, h5py, math, numpy as np

class DataTests(TestCase):

    maxDiff = None
    @validationTest
    def test_FChkFileDipoleSurface(self):
        fchk = TestManager.test_data("HOD_freq.fchk")
        surf = DipoleSurface.from_fchk_file(fchk)
        surf_center = surf.surfs[0].base.data['center']
        self.assertIsInstance(surf_center, np.ndarray)
        self.assertTrue(
            np.allclose(surf(surf_center) - np.array([s.base.data['ref'] for s in surf.surfs]), 0.)
        )
        self.assertEquals(surf([[0, 0, 0], [1, 0, 0], [0, 1, 0]]).shape, (1, 3))
        self.assertEquals(surf([
            [[0, 0, 0], [1, 0, 0], [0, 1, 0]],
            [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
        ]).shape, (2, 3))

    @validationTest
    def test_LogFileDipoleSurface(self):
        log = TestManager.test_data("water_OH_scan.log")
        conv = lambda x: cartesian_to_zmatrix(
            x, ordering=[[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]]
        ).coords[:, 0, 0]
        surf = DipoleSurface.from_log_file(log, conv)
        dips = surf(np.arange(.5, 2, .1))
        self.assertEquals(dips.shape, ((2-.5)/.1, 3))

        # surf_center = surf.surfs[0].base.data['center']
        # self.assertIsInstance(surf_center, np.ndarray)
        # self.assertTrue(
        #     np.allclose(surf(surf_center) - np.array([s.base.data['ref'] for s in surf.surfs]), 0.)
        # )
        # self.assertEquals(surf([[0, 0, 0], [1, 0, 0], [0, 1, 0]]).shape, (1, 3))
        # self.assertEquals(surf([
        #     [[0, 0, 0], [1, 0, 0], [0, 1, 0]],
        #     [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
        # ]).shape, (2, 3))

    @validationTest
    def test_LogFilePotentialSurface(self):
        log = TestManager.test_data("water_OH_scan.log")
        conv = lambda x: np.linalg.norm(x[:, 0] - x[:, 1], axis=1)
        surf = PotentialSurface.from_log_file(log, conv)
        pots = surf(np.arange(.5, 2, .1))
        self.assertEquals(pots.shape, ((2-.5)/.1,))

    def check_expr(self, test, real, raise_error=True):
        if not isinstance(test, int):
            test = test.expand().simplify().expand()
        if not isinstance(real, int):
            real = real.expand().simplify().expand()
        diff = (test - real)
        if not isinstance(diff, int):
            diff = diff.expand().simplify().expand()
        msg = "\nTest: {}\nReal: {}\nDiff: {}".format(test, real, diff)
        if not raise_error:
            return msg
        else:
            self.assertEquals(diff, 0, msg=msg)

    @validationTest
    def test_GmatrixElements(self):
        import sympy as sym
        class SymbolicCaller:
            def __init__(self, sym):
                self.sym = sym
            def __getitem__(self, item):
                if isinstance(item, int):
                    return self.sym(item)
                else:
                    return self.sym(*item)
        m = SymbolicCaller(AnalyticGMatrixConstructor.symbolic_m)
        r = SymbolicCaller(AnalyticGMatrixConstructor.symbolic_r)
        a = SymbolicCaller(AnalyticGMatrixConstructor.symbolic_a)
        t = SymbolicCaller(AnalyticGMatrixConstructor.symbolic_t)
        y = SymbolicCaller(AnalyticGMatrixConstructor.symbolic_y)
        sin = sym.sin; cos = sym.cos; cot = sym.cot; tan = sym.tan
        L = SymbolicCaller(AnalyticGMatrixConstructor.lam)

        self.check_expr(
            AnalyticGMatrixConstructor.g([1, 2], [1, 2]),
            1 / m[1] + 1 / m[2]
        )
        self.check_expr(
            AnalyticGMatrixConstructor.g([1, 2], [1, 3]),
            cos(a[2,1,3])/m[1]
        )
        self.check_expr(
            AnalyticGMatrixConstructor.g([1, 2], [3, 4]),
            0
        )
        self.check_expr(
            AnalyticGMatrixConstructor.g([1, 2], [1, 2, 3]),
            -sin(a[1,2,3])/(m[2]*r[2,3])
        )
        self.check_expr(
            AnalyticGMatrixConstructor.g([1, 2], [1, 3, 4]),
            sin(a[2,1,3])*cos(t[2,1,3,4])/(m[1]*r[1,3])
        )
        self.check_expr(
            # This requires some non-automatic simplifications...
            AnalyticGMatrixConstructor.g([1, 2], [3, 1, 4]),
            -sin(a[3, 1, 4])/m[1]*(
                              1/r[1, 4]*cos(a[2, 1, 3])*sin(a[3, 1, 4])
                              + L[3, 1, 4]*sin(a[2, 1, 3])*cos(y[3, 1, 2, 4])
                          )
        )
        self.check_expr(
            AnalyticGMatrixConstructor.g([1, 2], [1, 2, 3, 4]),
            -sin(a[1, 2, 3])*sin(t[1, 2, 3, 4])*cot(a[2, 3, 4])/(m[2]*r[2,3])
        )
        self.check_expr(
            AnalyticGMatrixConstructor.g([1, 2], [3, 2, 1, 4]),
            0
        )
        self.check_expr(
            AnalyticGMatrixConstructor.g([1, 2], [1, 3, 4, 5]),
            -sin(a[2, 1, 3])*sin(t[2, 1, 3, 4])/(m[1]*r[1,3]*sin(a[1,3,4]))
        )
        self.check_expr(
            AnalyticGMatrixConstructor.g([1, 2], [4, 3, 1, 5]),
            -sin(a[2, 1, 3]) / m[1] * (
                    cot(a[1, 3, 4]) / r[1, 3] * sin(t[4, 3, 1, 2])
                     + L[5, 1, 3] * sin(t[4, 3, 1, 5] - t[4, 3, 1, 2])
            )
        )
        self.check_expr(
            AnalyticGMatrixConstructor.g([1, 2, 3], [1, 2, 3]),
            1/(m[1]*r[1,2]**2) + 1/(m[3]*r[2,3]**2)
            + 1/m[2]*(1/r[1,2]**2 + 1/r[2,3]**2 - 2*cos(a[1,2,3])/(r[1,2]*r[2,3]))
        )
        self.check_expr(
            AnalyticGMatrixConstructor.g([1, 2, 3], [2, 1, 4]),
            -1 / r[1, 2] * cos(t[3, 2, 1, 4]) * (
                L[2, 1, 4] / m[1] * sin(a[2, 1, 4]) + L[1, 2, 3] / m[2] * sin(a[1, 2, 3])
            )
        )
        self.check_expr(
            AnalyticGMatrixConstructor.g([1, 2, 3], [1, 2, 4]),
            1/m[2]*sin(a[1, 2, 3])*sin(a[1, 2, 4])*(
                L[1, 2, 3]*L[1, 2, 4]*cos(y[1, 2, 3, 4]) + 1/(r[2, 3]*r[2, 4])
            ) + 1/(m[1] * r[1, 2]**2)*cos(y[1, 2, 3, 4])
        )
        self.check_expr(
            AnalyticGMatrixConstructor.g([1, 2, 3], [1, 2, 3, 4]),
            1/r[2, 3]*sin(t[1, 2, 3, 4])*(
                L[3, 2, 1]/m[2]*sin(a[1, 2, 3])*cot(a[2, 3, 4]) - L[4, 3, 2]/m[3]
            )
        )
