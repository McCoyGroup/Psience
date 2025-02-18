
from Peeves.TestUtils import *
from Psience.AnalyticModels import *
from McUtils.Plots import *
import McUtils.Numputils as nput
from unittest import TestCase
import sys, h5py, math, numpy as np

class AnalyticModelsTests(TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        np.set_printoptions(linewidth=int(1e8))

    maxDiff = None

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

    @debugTest
    def test_RedundantG(self):
        from Psience.Molecools import Molecule
        from Psience.Psience.AnalyticModels.AnalyticJacobianDotCalculator import InternalJacobianDisplacements
        from Psience.AnalyticModels import AnalyticKineticEnergyConstructor, GeometricFunction

        file_name = "nh3.fchk"
        test_internals = [[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1], [3, 0, 1, 2]]
        mol = Molecule.from_file(TestManager.test_data(file_name), internals=test_internals)

        coords = [
            (0, 1)
            , (2, 0)
            , (0, 3)
            , (2, 0, 1)
            , (3, 0, 1)
            , (2, 0, 3)
            , (2, 1, 0)
            , (3, 1, 0)
            , (3, 1, 2)
            , (1, 2)
            , (1, 3)
            , (2, 3)
            , (3, 0, 1, 2)
            , (1, 0, 2, 3)
            , (2, 0, 3, 1)
            , (0, 1, 3, 2)
        ]

        B = np.array([
            (
                nput.dist_vec(mol.coords, *inds)
                if len(inds) == 2 else
                nput.angle_vec(mol.coords, inds[1], inds[0], inds[2])
                if len(inds) == 3 else
                nput.dihed_vec(mol.coords, *inds)
                if len(inds) == 4 else
                nput.book_vec(mol.coords, *inds)
                if len(inds) == 5 and inds[-1] == -1 else
                None
            )
            for inds in coords
        ])

        M = np.diag(np.repeat(1 / np.sqrt(mol.atomic_masses)[:, np.newaxis], 3))

        # 2.41247616e-05
        # raise Exception(
        #     (lambda b: b @ M @ M @ b.T)(np.array([
        #             nput.angle_vec(mol.coords, 1, 2, 0),
        #             nput.dihed_vec(mol.coords, 2, 0, 3, 1)
        #         ])),
        #     # AnalyticKineticEnergyConstructor.g(
        #     #     (2, 0), (3, 0, 1, 2),
        #     #     return_function=True,
        #     #     method='direct'
        #     # )(mol.atomic_masses, mol.coords),
        #     AnalyticKineticEnergyConstructor.g(
        #         (2, 1, 0), (2, 0, 3, 1),
        #         return_function=True,
        #         method='direct'
        #     )(mol.atomic_masses, mol.coords)
        # )

        raise Exception(
            AnalyticKineticEnergyConstructor.g(
                    (1, 3), (1, 2, 3),
                    method='direct'
                )
        )

        # raise Exception(
        #     (lambda b:b@M@M@b.T)(np.array([
        #         nput.angle_vec(mol.coords, 0, 3, 1),
        #         nput.dihed_vec(mol.coords, 3, 0, 1, 2)
        #     ])),
        #     # AnalyticKineticEnergyConstructor.g(
        #     #     (3, 0), (3, 0),
        #     #     return_function=True,
        #     #     method='direct'
        #     # )(mol.atomic_masses, mol.coords),
        #     # AnalyticKineticEnergyConstructor.g(
        #     #     (3, 0, 1, 2), (3, 0, 1, 2),
        #     #     return_function=True,
        #     #     method='direct'
        #     # )(mol.atomic_masses, mol.coords),
        #     AnalyticKineticEnergyConstructor.g(
        #         (2, 3), (3, 0, 1, 2),
        #         return_function=True,
        #         method='direct'
        #     )(mol.atomic_masses, mol.coords),
        #     # AnalyticKineticEnergyConstructor.g(
        #     #     (2, 3), (3, 0, 1, 2),
        #     #     method='direct'
        #     # ),
        #     # AnalyticKineticEnergyConstructor.g(
        #     #     (3, 0, 1), (3, 0, 1, 2),
        #     #     return_function=True,
        #     #     method='direct'
        #     # )(mol.atomic_masses, mol.coords),
        #     # AnalyticKineticEnergyConstructor.g(
        #     #     (3, 0, 1), (3, 0, 1, 2),
        #     #     method='direct'
        #     # )
        # )

        # print(B)
        # for c in coords:
        #     print(
        #         [
        #             GeometricFunction.from_expr(e.as_expr())(
        #                 mol.atomic_masses, mol.coords
        #             )
        #             for e in InternalJacobianDisplacements.displacement_vectors(c)
        #         ]
        #     )
        # raise Exception(...)

        B = B @ M

        g_direct = B @ B.T

        g_red = AnalyticKineticEnergyConstructor.g_matrix(
            coords,
            return_function=True,
            method='direct'
        )(mol.atomic_masses, mol.coords)

        # print(mol.g_matrix)
        print((g_direct,))
        print((g_red,))
        print(np.round(g_red - g_direct, 8))

        print(np.linalg.eigvalsh(g_direct))
        print(np.linalg.eigvalsh(g_red))

        raise Exception(...)
        #     np.linalg.eigh(g_red)[0] * 219475.6
        # )

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
        m = SymbolicCaller(AnalyticKineticEnergyConstructor.symbolic_m)
        r = SymbolicCaller(AnalyticKineticEnergyConstructor.symbolic_r)
        a = SymbolicCaller(AnalyticKineticEnergyConstructor.symbolic_a)
        t = SymbolicCaller(AnalyticKineticEnergyConstructor.symbolic_t)
        y = SymbolicCaller(AnalyticKineticEnergyConstructor.symbolic_y)
        sin = sym.sin; cos = sym.cos; cot = sym.cot; tan = sym.tan
        L = SymbolicCaller(AnalyticKineticEnergyConstructor.lam)

        self.check_expr(
            AnalyticKineticEnergyConstructor.g([1, 2], [1, 2]),
            1 / m[1] + 1 / m[2]
        )
        self.check_expr(
            AnalyticKineticEnergyConstructor.g([1, 2], [1, 3]),
            cos(a[2,1,3])/m[1]
        )
        self.check_expr(
            AnalyticKineticEnergyConstructor.g([1, 2], [3, 4]),
            0
        )
        self.check_expr(
            AnalyticKineticEnergyConstructor.g([1, 2], [1, 2, 3]),
            -sin(a[1,2,3])/(m[2]*r[2,3])
        )
        self.check_expr(
            AnalyticKineticEnergyConstructor.g([1, 2], [1, 3, 4]),
            sin(a[2,1,3])*cos(t[2,1,3,4])/(m[1]*r[1,3])
        )
        self.check_expr(
            # This requires some non-automatic simplifications...
            AnalyticKineticEnergyConstructor.g([1, 2], [3, 1, 4]),
            -sin(a[3, 1, 4])/m[1]*(
                              1/r[1, 4]*cos(a[2, 1, 3])*sin(a[3, 1, 4])
                              + L[3, 1, 4]*sin(a[2, 1, 3])*cos(y[3, 1, 2, 4])
                          )
        )
        self.check_expr(
            AnalyticKineticEnergyConstructor.g([1, 2], [1, 2, 3, 4]),
            -sin(a[1, 2, 3])*sin(t[1, 2, 3, 4])*cot(a[2, 3, 4])/(m[2]*r[2,3])
        )
        self.check_expr(
            AnalyticKineticEnergyConstructor.g([1, 2], [3, 2, 1, 4]),
            0
        )
        self.check_expr(
            AnalyticKineticEnergyConstructor.g([1, 2], [1, 3, 4, 5]),
            -sin(a[2, 1, 3])*sin(t[2, 1, 3, 4])/(m[1]*r[1,3]*sin(a[1,3,4]))
        )
        self.check_expr(
            AnalyticKineticEnergyConstructor.g([1, 2], [4, 3, 1, 5]),
            -sin(a[2, 1, 3]) / m[1] * (
                    cot(a[1, 3, 4]) / r[1, 3] * sin(t[4, 3, 1, 2])
                     + L[5, 1, 3] * sin(t[4, 3, 1, 5] - t[4, 3, 1, 2])
            )
        )
        self.check_expr(
            AnalyticKineticEnergyConstructor.g([1, 2, 3], [1, 2, 3]),
            1/(m[1]*r[1,2]**2) + 1/(m[3]*r[2,3]**2)
            + 1/m[2]*(1/r[1,2]**2 + 1/r[2,3]**2 - 2*cos(a[1,2,3])/(r[1,2]*r[2,3]))
        )
        self.check_expr(
            AnalyticKineticEnergyConstructor.g([1, 2, 3], [2, 1, 4]),
            -1 / r[1, 2] * cos(t[3, 2, 1, 4]) * (
                L[2, 1, 4] / m[1] * sin(a[2, 1, 4]) + L[1, 2, 3] / m[2] * sin(a[1, 2, 3])
            )
        )
        self.check_expr(
            AnalyticKineticEnergyConstructor.g([1, 2, 3], [1, 2, 4]),
            1/m[2]*sin(a[1, 2, 3])*sin(a[1, 2, 4])*(
                L[1, 2, 3]*L[1, 2, 4]*cos(y[1, 2, 3, 4]) + 1/(r[2, 3]*r[2, 4])
            ) + 1/(m[1] * r[1, 2]**2)*cos(y[1, 2, 3, 4])
        )
        self.check_expr(
            AnalyticKineticEnergyConstructor.g([1, 2, 3], [1, 2, 3, 4]),
            1/r[2, 3]*sin(t[1, 2, 3, 4])*(
                L[3, 2, 1]/m[2]*sin(a[1, 2, 3])*cot(a[2, 3, 4]) - L[4, 3, 2]/m[3]
            )
        )
