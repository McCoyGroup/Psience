import collections
import gc
import itertools
import math
import os

import scipy.linalg

import McUtils.Zachary
from Peeves.TestUtils import *
from unittest import TestCase

from McUtils.Data import UnitsData, PotentialData, AtomData
from McUtils.Zachary import Interpolator, FiniteDifferenceDerivative
import McUtils.Plots as plt
import McUtils.Numputils as nput
from McUtils.GaussianInterface import GaussianLogReader
from McUtils.Extensions import ModuleLoader

from McUtils.Scaffolding import Checkpointer

from Psience.Modes import *
from Psience.Molecools import *

import numpy as np

class LocalTests(TestCase):

    def setUp(self) -> None:
        np.set_printoptions(linewidth=1e8)

    @validationTest
    def test_Water(self):
        mol = Molecule.from_file(
            TestManager.test_data('HOD_freq.fchk'),
            internals=[
                [0, -1, -1, -1],
                [1,  0, -1, -1],
                [2,  0,  1, -1]
            ]
        )

        ortho = ObliqueModeGenerator.from_molecule(mol, frequency_scaled=True)

        f, g, u, ui = ortho.run()

        import scipy

        g12 = scipy.linalg.fractional_matrix_power(ortho.g, 1/2)
        gi12 = scipy.linalg.fractional_matrix_power(ortho.g, -1/2)
        # freq2, L = scipy.linalg.eigh(ortho.f, ortho.g, type=3)
        # raise Exception(...)
        # print(ui)
        # gfg = g12 @ ortho.f @ g12
        gvals, gvecs = np.linalg.eigh(ortho.g)
        gvals = gvals[(2, 1, 0),]
        gvecs = gvecs[:, (2, 1, 0)]

        # print(g12)
        # print(gvecs@np.diag(np.power(gvals, 1/2))@gvecs.T)
        # raise Exception(...)

        fvals, fvecs = np.linalg.eigh(ortho.f)
        fvecs = fvecs[:, (2, 1, 0)]
        fvals = fvals[(2, 1, 0),]

        freq2, Q = np.linalg.eigh(g12 @ ortho.f @ g12)
        freq2 = freq2[(2, 1, 0),]
        Q = Q[:, (2, 1, 0)]

        F = np.diag(1 / np.power(freq2, 1/4))
        FF = np.diag(1 / np.power(freq2, 1/8))
        Fi = np.diag(np.power(freq2, 1/4))
        L = g12 @ Q @ F
        l12 = np.diag(np.power(gvals, 1/2))
        print("="*50)
        LU, LS, LV = np.linalg.svd(g12 @ Q)
        print(np.power(gvals, 1/2))
        print(l12@F)
        # print(
        #     np.linalg.svd(l12@np.diag([1, -1, 1])@F)[1]
        # )
        print('-'*50)

        print(ui)

        print('-'*50)

        U, S, V = np.linalg.svd(L)
        print(S)
        print(
            np.linalg.svd(l12@gvecs.T@Q@F)[1]
        )
        print('-' * 50)
        print(U)
        print(V)
        print(U.T @ V)
        print(gvecs.T @ U)
        print('-' * 50)
        print(LV)
        print(gvecs)
        print(fvecs)
        raise Exception(...)
        # LU, LS, LV = np.linalg.svd(np.diag(gvals) @ LV @)
        print('-'*50)
        print(U)
        print(gvecs)
        print(fvecs)
        print('-'*50)
        print(ui)
        print(U @ np.diag(S) @ U.T)
        raise Exception(...)
        # print(LU)
        # print("-"*50)
        # print(gvecs)
        # print(LS, np.sqrt(gvals))
        print("="*50)
        # print(LV)
        # print("-"*50)
        # print(gvecs.T @ Q)
        # print("-"*50)
        print(ui)
        print("-"*50)
        print(g12)
        # print("-"*50)
        # print(gvecs @ F@np.diag(np.power(gvals, 1/2)) @ gvecs.T)
        print("-"*50)
        print(Q)
        print("-"*50)
        print(gvecs)
        print("-"*50)
        print(np.linalg.eigh(ui)[0])
        print(np.linalg.eigh(ui)[1])
        print("-"*50)
        print(gvals / fvals)

        raise Exception(...)
        L = g12 @ Q @ F
        Li = Fi @ Q.T @ gi12
        print("-"*50)
        print(L.T @ ortho.f @ L)
        print("-"*50)
        print(Li @ ortho.g @ Li.T)
        print("="*50)
        print(ortho.rotation)
        print("-"*50)
        print(Q)
        print("-"*50)
        print(Q.T @ ortho.rotation)
        print("-"*50)
        V = gvecs.T @ Q
        print("-"*50)
        print(V.T @ np.diag(np.power(gvals, 1/2)) @ V)
        print("-"*50)
        print(ui)
        # print(freq3, freq2)
        raise Exception(...)
        print(Qgfg.T @ ortho.f @ Qgfg + Qgfg.T @ ortho.g @ Qgfg)

        F = np.diag([7, 5, 3])
        print(np.linalg.svd(ortho.f @ F)[1])
        print(np.linalg.svd(ortho.f)[1])

        raise Exception(...)



        # print(L @ np.diag(1/np.power(freq2, 1/4)))
        # print(L / np.power(freq2, 1/4)[np.newaxis, :])
        F = np.diag(1 / np.power(freq2, 1/4))
        rU, S, rV = np.linalg.svd(L)
        print(...)
        gvals, gvecs = np.linalg.eigh(ortho.g)
        gvals = gvals[(2, 1, 0),]
        gvecs = gvecs[:, (2, 1, 0)]
        fvals, fvecs = np.linalg.eigh(ortho.f)
        fvecs = fvecs[:, (2, 1, 0)]
        fvals = fvals[(2, 1, 0),]

        fvals, fvecs = np.linalg.eigh(ortho.f)

        # print(rU)
        print(rV)
        print(Qgfg)

        raise Exception(...)
        # R1 = rU @ rV

        # rU3, S3, rV3 = np.linalg.svd(np.diag(S) @ rV @ F)
        # print(S3)
        # print(F @ np.diag(S))

        # print(R1)
        LF = L @ F
        rU2, S2, rV2 = np.linalg.svd(LF)
        # print(rU)
        # print(rU2)
        print(np.diag(S) @ rV @ F)
        raise Exception(...)
        freq12 = np.diag(np.power(freq2, 1/8))
        print(g12)
        print(ui)
        raise Exception(...)

        s, l = np.linalg.eigh(f)
        print(s)
        print(l)
        print(ortho.rotation)
        print(ui)
        print("-"*20)
        H = ortho.modes @ ortho.modes.T
        sh, qh = np.linalg.eigh(H)
        P = qh @ np.diag(np.sqrt(sh)) @ qh.T
        print(P)
        s, l = np.linalg.eigh(P@ortho.f@P)
        print(s)
        print(l)


        # print("-"*50)
        # print(u)
        #
        # raise Exception(...)


        # freq2, modes = scipy.linalg.eigh(ortho.f, ortho.g, type=3)
        # freqs = np.sqrt(freq2)
        # modes = modes / np.sqrt(freqs)[np.newaxis, :]
        # print(np.linalg.inv(modes))

        # print(np.linalg.norm((f - g).flatten()))
        # print("-"*50)
        # print(f)
        # print("-"*50)
        # print(ortho.f)
        # print("-"*50)
        # freq2, modes = scipy.linalg.eigh(ortho.f, ortho.g, type=3)
        # print(np.sqrt(freq2))
        # print("-"*50)
        # print(u)
        # print("-"*50)
        # print(ui)

        self.assertLess(abs(np.linalg.norm((f - g).flatten())), 1e-15)

    @validationTest
    def test_OCHH(self):
        zmatrix = [
            [0, -1, -1, -1],
            [1,  0, -1, -1],
            [2,  1,  0, -1],
            [3,  1,  0,  2]
        ]

        mol = Molecule.from_file(
            TestManager.test_data('OCHH_freq.fchk'),
            internals=zmatrix
        )

        iterative = False

        ortho = ObliqueModeGenerator.from_molecule(mol, sel=[0, 1, 3, 2, 4])
        f, g, u, ui = ortho.run()
        # print(np.linalg.norm((f - g).flatten()))
        # print("-"*50)
        # print(u)
        # print("-"*50)
        # print(ui)

        self.assertLess(abs(np.linalg.norm((f - g).flatten())), 1e-15)

        from McUtils.Misc import TeX

        def format_mat(lhs, m, label=None, digits=2):
            f_og = m.astype(object)
            f_og[np.tril_indices_from(f_og, -1)] = ""
            fsym = TeX.bold(lhs).as_expr()
            TeX.Writer.real_digits = digits
            fexpr = fsym.Eq(TeX.Matrix(f_og))
            return TeX.Equation(fexpr, label=label).format_tex()

        sys = 'ochh'
        print(
            format_mat('f', ortho.f * UnitsData.convert("Hartrees", "Wavenumbers"), label='f_'+sys)
        )
        print(
            format_mat('g', ortho.g * UnitsData.convert("Hartrees", "Wavenumbers"), label='g_'+sys)
        )
        print(
            format_mat('F', f * UnitsData.convert("Hartrees", "Wavenumbers"), label='F_'+sys)
        )
        print(
            format_mat('P', u, label='P_'+sys, digits=3)
        )

        # g_test = u@ortho.g@u
        # print(
        #     format_mat('G_scaled', g_test * UnitsData.convert("Hartrees", "Wavenumbers"), label='G_ochh')
        # )

    @validationTest
    def test_HOONO(self):
        zmatrix  = [
            [1, -1, -1, -1],
            [2,  1, -1, -1],
            [3,  2,  1, -1],
            [0,  1,  2,  3],
            [4,  3,  2,  1]
        ]

        mol = Molecule.from_file(
            TestManager.test_data('HOONO_freq.fchk'),
            internals=zmatrix
        )

        ortho = ObliqueModeGenerator.from_molecule(mol)
        f, g, u, ui = ortho.run()
        # print(np.linalg.norm((f - g).flatten()))
        # print("-"*50)
        # print(u)
        # print("-"*50)
        # print(ui)


        self.assertLess(abs(np.linalg.norm((f - g).flatten())), 1e-15)

    @validationTest
    def test_NH3(self):
        tag = 'NH3 Internals'
        file_name = "nh3.fchk"

        _ = -1
        zmatrix = [
            [0, _, _, _],
            [1, 0, _, _],
            [2, 0, 1, _],
            [3, 0, 1, 2]
        ]

        def conv(r, t, f, **kwargs):
            cp1 = np.cos(f[..., 3])  # skip three for embedding
            ct1 = np.cos(t[..., 2])  # skip two for embedding
            ct2 = np.cos(t[..., 3])
            st1 = np.sin(t[..., 2])
            st2 = np.sin(t[..., 3])
            f[..., 3] = np.arccos(st1 * st2 * cp1 + ct1 * ct2)
            return np.array([r, t, f])

        def inv(r, t, f, **kwargs):
            cp1 = np.cos(f[..., 3])
            ct1 = np.cos(t[..., 2])
            ct2 = np.cos(t[..., 3])
            st1 = np.sin(t[..., 2])
            st2 = np.sin(t[..., 3])
            f[..., 3] = np.arccos((cp1 - ct1 * ct2) / (st1 * st2))
            return np.array([r, t, f])

        # def conv(crds, **kwargs):
        #     return crds
        # def inv(crds, **kwargs):
        #     return crds


        internals = {
            'zmatrix':zmatrix,
            'conversion':conv,
            'inverse':inv,
            # 'converter_options':{
            #     'pointwise':False,
            #     # 'jacobian_prep':ZMatrixCoordinateSystem.jacobian_prep_coordinates
            # }
        }

        mol = Molecule.from_file(
            TestManager.test_data('nh3.fchk'),
            internals=zmatrix
        )
        ortho = ObliqueModeGenerator.from_molecule(mol)
        f, g, u, ui = ortho.run()

        self.assertLess(abs(np.linalg.norm((f - g).flatten())), 1e-15)

        mol = Molecule.from_file(
            TestManager.test_data('nh3.fchk'),
            internals=zmatrix
        )
        ortho = ObliqueModeGenerator.from_molecule(mol, sel=[0, 1, 3, 2, 4, 5])
        f, g, u, ui = ortho.run()


        from McUtils.Misc import TeX

        def format_mat(lhs, m, label=None, digits=1):
            f_og = m.astype(object)
            f_og[np.tril_indices_from(f_og, -1)] = ""
            fsym = TeX.bold(lhs).as_expr()
            TeX.Writer.real_digits = digits
            fexpr = fsym.Eq(TeX.Matrix(f_og))
            return TeX.Equation(fexpr, label=label).format_tex()

        sys = 'nh3_nosymm'
        print(
            format_mat('F', ortho.f * UnitsData.convert("Hartrees", "Wavenumbers"), label='f_' + sys)
        )
        print(
            format_mat('G', ortho.g * UnitsData.convert("Hartrees", "Wavenumbers"), label='g_' + sys)
        )
        print(
            format_mat('F', f * UnitsData.convert("Hartrees", "Wavenumbers"), label='F_' + sys)
        )
        print(
            format_mat('A', u, label='P_' + sys, digits=3)
        )

        mol = Molecule.from_file(
            TestManager.test_data('nh3.fchk'),
            internals=internals
        )

        ortho = ObliqueModeGenerator.from_molecule(mol, sel=[0, 1, 3, 2, 4, 5])
        f, g, u, ui = ortho.run()
        # print(np.linalg.norm((f - g).flatten()))
        # print("-"*50)
        # print(u)
        # print("-"*50)
        # print(ui)

        self.assertLess(abs(np.linalg.norm((f - g).flatten())), 1e-15)

        from McUtils.Misc import TeX

        def format_mat(lhs, m, label=None, digits=1):
            f_og = m.astype(object)
            f_og[np.tril_indices_from(f_og, -1)] = ""
            fsym = TeX.bold(lhs).as_expr()
            TeX.Writer.real_digits = digits
            fexpr = fsym.Eq(TeX.Matrix(f_og))
            return TeX.Equation(fexpr, label=label).format_tex()

        print("="*50)
        sys = 'nh3_symm'
        print(
            format_mat('F', ortho.f * UnitsData.convert("Hartrees", "Wavenumbers"), label='f_' + sys)
        )
        print(
            format_mat('G', ortho.g * UnitsData.convert("Hartrees", "Wavenumbers"), label='g_' + sys)
        )
        print(
            format_mat("F", f * UnitsData.convert("Hartrees", "Wavenumbers"), label='F_' + sys)
        )
        print(
            format_mat('A', u, label='P_' + sys, digits=3)
        )

    @validationTest
    def test_Dimer(self):
        COM = -3
        A = -2
        C = -1
        X = 1000
        LHF = 0
        LO = 1
        SH = 2
        RO = 3
        RH1 = 4
        RH2 = 5

        internals = [
            [LHF,   X,   X,   X],
            [ LO, LHF,   X,   X],
            [ SH,  LO, LHF,   X],
            [ RO,  LO, LHF,  SH],
            [RH1,  RO,  LO, LHF],
            [RH2,  RO, RH1,  LO]
        ]
        sel = [0, 1, 6, 9, 3, #stretches
               5,  #SH OOP
               2,  #SH bend
               4,  #HF_OO bend
               10, #HOH bend
               7,  #H1_OO bend
               8,  #HOO-HF OOP
               11  #HOH-HF OOP
               ]#[0, 1, 9, 6, 2, 10, 3, 4, 7, 8, 5, 11]

        dimer = Molecule.from_file(
            TestManager.test_data('water_dimer_freq.fchk'),
            internals=internals
        )

        ics = dimer.internal_coordinates * np.array([
            UnitsData.convert("BohrRadius", "Angstroms"),
            180/np.pi,
            180/np.pi
        ])[np.newaxis]
        ics = ics.flatten()
        bad_coords = dimer._get_embedding_coords()
        good_coords = np.setdiff1d(np.arange(len(ics)), bad_coords)
        ics_1 = ics[good_coords][sel]


        pd = dimer.potential_derivatives
        dimer.potential_derivatives = [
            np.zeros_like(pd[0]),
            pd[1]
        ]

        ortho = ObliqueModeGenerator.from_molecule(dimer, sel=sel)
        f, g, u, ui = ortho.run()
        # print("-"*50)
        # print(f)
        # print("-"*50)
        # print(u)
        # print("-"*50)
        # print(ui)

        print(
            "Oblique Monomer 1",
            np.linalg.eigvalsh(f[:, [0, 1, 6]][[0, 1, 6], :]) * 219475.6
        )

        self.assertLess(abs(np.linalg.norm((f - g).flatten())), 1e-15)

        from Psience.VPT2 import VPTRunner
        # VPTRunner.run_simple(
        #     TestManager.test_data('water_dimer_freq.fchk'),
        #     2,
        #     internals=internals
        # )
        """
0 0 0 0 0 0 0 0 0 0 0 1    3935.27376     97.33660      3753.87931     57.72733
0 0 0 0 0 0 0 0 0 0 1 0    3915.15761    114.99756      3741.00783     82.38512
0 0 0 0 0 0 0 0 0 1 0 0    3813.91883     11.14760      3647.24782      6.15614
0 0 0 0 0 0 0 0 1 0 0 0    3718.74059    296.76286      3582.53802    185.65100
0 0 0 0 0 0 0 1 0 0 0 0    1650.27150     36.30532      1596.35380     39.01877
0 0 0 0 0 0 1 0 0 0 0 0    1629.31503     87.42836      1585.55269     67.22019
0 0 0 0 0 1 0 0 0 0 0 0     630.33014     91.54623       539.18888     91.11127
0 0 0 0 1 0 0 0 0 0 0 0     360.39517     51.44266       305.47782      0.73307
0 0 0 1 0 0 0 0 0 0 0 0     184.18194    152.66165       155.93803    185.48737
0 0 1 0 0 0 0 0 0 0 0 0     154.94037    128.25909       137.66338     65.62913
0 1 0 0 0 0 0 0 0 0 0 0     147.11043     64.66132       152.23110     11.18758
1 0 0 0 0 0 0 0 0 0 0 0     127.00706    103.95165       110.15448    153.49999
0 0 0 0 0 2 0 0 0 0 0 0    1260.66028      0.00000      1013.58568     15.00090
0 0 0 0 2 0 0 0 0 0 0 0     720.79034      0.00000       555.63944      6.08444
0 0 0 2 0 0 0 0 0 0 0 0     368.36387      0.00000       314.48326     41.36952
0 0 2 0 0 0 0 0 0 0 0 0     309.88075      0.00000       265.29283     20.96075
0 2 0 0 0 0 0 0 0 0 0 0     294.22085      0.00000       307.62441      3.29242
0 0 0 1 0 0 0 0 1 0 0 0    3902.92253      0.00000      3752.77368     22.29680
1 1 0 0 0 0 0 0 0 0 0 0     274.11749      0.00000       255.48766     13.27823
"""

        from McUtils.Misc import TeX

        def format_mat(lhs, m, label=None, digits=1):
            f_og = m.astype(object)
            f_og[np.tril_indices_from(f_og, -1)] = ""
            fsym = TeX.bold(lhs).as_expr()
            TeX.Writer.real_digits = digits
            fexpr = fsym.Eq(TeX.Matrix(f_og))
            return TeX.Equation(fexpr, label=label).format_tex()

        print("=" * 50)
        print(np.sqrt(scipy.linalg.eigvalsh(ortho.f, ortho.g, type=3)) * 219475.6)

        sys = 'dimer_1'
        print(
            format_mat('f', ortho.f[5:, :][:, 5:] * UnitsData.convert("Hartrees", "Wavenumbers"), label='f_' + sys)
        )
        # print(
        #     format_mat('g', ortho.g * UnitsData.convert("Hartrees", "Wavenumbers"), label='g_' + sys)
        # )
        print(
            format_mat('F', f[5:, :][:, 5:] * UnitsData.convert("Hartrees", "Wavenumbers"), label='F_' + sys)
        )
        # print(
        #     format_mat('P', u, label='P_' + sys, digits=3)
        # )

        print(np.linalg.eigvalsh(f[5:, :][:, 5:]) * UnitsData.convert("Hartrees", "Wavenumbers"))

        print("="*50)

        internals = [
            [LHF,   X,   X,   X],
            [ LO, LHF,   X,   X],
            [ SH,  LO, LHF,   X],
            [ RO,  LO, LHF,  SH],
            # [ RO,  LO, LHF,   X],
            # [ SH,  LO,  RO, LHF],
            [RH1,  RO,  LO, LHF],
            [RH2,  RO,  LO, LHF]
        ]
        sel = [
            0, 3, 6, 9, 1, #stretches
            5, #SH_OOP
            2, #
            4, #
            7, #H1_OO bend
            10, #H2_OO bend
            8, #H1-OOHF OOP
            11 #H2-OOHF OOP
        ]

        # from Psience.VPT2 import VPTRunner
        # VPTRunner.run_simple(
        #     TestManager.test_data('water_dimer_freq.fchk'),
        #     2,
        #     internals=internals
        # )
        """
0 0 0 0 0 0 0 0 0 0 0 1    3935.27376     97.33660      3753.89114     57.73399
0 0 0 0 0 0 0 0 0 0 1 0    3915.15761    114.99756      3741.01251     82.41541
0 0 0 0 0 0 0 0 0 1 0 0    3813.91883     11.14760      3647.25710      6.15108
0 0 0 0 0 0 0 0 1 0 0 0    3718.74059    296.76285      3582.53252    185.63408
0 0 0 0 0 0 0 1 0 0 0 0    1650.27151     36.30532      1596.60966     39.22153
0 0 0 0 0 0 1 0 0 0 0 0    1629.31503     87.42836      1585.58546     67.16024
0 0 0 0 0 1 0 0 0 0 0 0     630.33868     91.54864       504.79327     85.17639
0 0 0 0 1 0 0 0 0 0 0 0     360.39519     51.44266       305.26515      0.67423
0 0 0 1 0 0 0 0 0 0 0 0     184.18194    152.66166       155.66394    182.99421
0 0 1 0 0 0 0 0 0 0 0 0     154.94037    128.25909       137.41059     67.44037
0 1 0 0 0 0 0 0 0 0 0 0     147.11065     64.66109       145.20484     34.43950
1 0 0 0 0 0 0 0 0 0 0 0     127.00715    103.95200       106.57198    111.69437
0 0 0 0 0 2 0 0 0 0 0 0    1260.67737      0.00000       919.65889     13.60571
0 0 0 0 2 0 0 0 0 0 0 0     720.79037      0.00000       555.23478      6.08001
0 0 0 2 0 0 0 0 0 0 0 0     368.36387      0.00000       313.89588     41.29244
0 0 2 0 0 0 0 0 0 0 0 0     309.88075      0.00000       264.75383     20.91816
0 2 0 0 0 0 0 0 0 0 0 0     294.22130      0.00000       293.22913      3.13868
0 0 0 1 0 0 0 0 1 0 0 0    3902.92253      0.00000      3752.46724     22.29499
1 1 0 0 0 0 0 0 0 0 0 0     274.11780      0.00000       244.55840     12.71015
"""

        dimer = Molecule.from_file(
            TestManager.test_data('water_dimer_freq.fchk'),
            internals=internals
        )

        ics = dimer.internal_coordinates * np.array([
            UnitsData.convert("BohrRadius", "Angstroms"),
            180 / np.pi,
            180 / np.pi
        ])[np.newaxis]
        ics = ics.flatten()
        bad_coords = dimer._get_embedding_coords()
        good_coords = np.setdiff1d(np.arange(len(ics)), bad_coords)
        ics_2 = ics[good_coords][sel]

        TeX.Writer.real_digits = 3
        # raise Exception(
        #     TeX.Array(
        #     np.array(
        #         [
        #             np.concatenate([ics_1, ics_2]),
        #             np.concatenate([ics_1, ics_2])
        #             ]
        #     ).T
        # ).format_tex()
        # )

        pd = dimer.potential_derivatives
        dimer.potential_derivatives = [
            np.zeros_like(pd[0]),
            pd[1]
        ]

        ortho = ObliqueModeGenerator.from_molecule(dimer, sel=sel)
        print(np.sqrt(scipy.linalg.eigvalsh(ortho.f, ortho.g, type=3)) * 219475.6)
        f, g, u, ui = ortho.run()

        print(
            "Oblique Monomer 2",
            np.linalg.eigvalsh(f[:, [0, 4, 6]][[0, 4, 6], :]) * 219475.6
        )

        sys = 'dimer_2'
        print(
            format_mat('f', ortho.f[5:, :][:, 5:] * UnitsData.convert("Hartrees", "Wavenumbers"), label='f_' + sys)
        )
        # print(
        #     format_mat('g', ortho.g * UnitsData.convert("Hartrees", "Wavenumbers"), label='g_' + sys)
        # )
        # sel = [0, ]
        print(
            format_mat('F', f[5:, :][:, 5:] * UnitsData.convert("Hartrees", "Wavenumbers"), label='F_' + sys)
        )
        # print(
        #     format_mat('P', u, label='P_' + sys, digits=3)
        # )

        print(
            np.sqrt(
                scipy.linalg.eigvalsh(ortho.f, ortho.g, type=3)
            ) * UnitsData.convert("Hartrees", "Wavenumbers")
        )
        print(
            np.sqrt(
                scipy.linalg.eigvalsh(ortho.f[:5, :][:, :5], ortho.g[:5, :][:, :5], type=3)
            ) * UnitsData.convert("Hartrees", "Wavenumbers")
        )
        print(np.linalg.eigvalsh(f[:5, :][:, :5]) * UnitsData.convert("Hartrees", "Wavenumbers"))
        print(
            np.sqrt(
                scipy.linalg.eigvalsh(ortho.f[5:, :][:, 5:], ortho.g[5:, :][:, 5:], type=3)
            ) * UnitsData.convert("Hartrees", "Wavenumbers")
        )
        print(np.linalg.eigvalsh(f[5:, :][:, 5:]) * UnitsData.convert("Hartrees", "Wavenumbers"))

        """
        [ 127.38278305  147.37823763  167.8307369   360.273703    630.61290521 1634.62791975 1654.82208119]
        [ 127.04607029  147.21749535  166.72211701  360.22808449  630.47924308 1636.68361974 1656.14710237]
        """

        # subdimer = ObliqueModeGenerator(
        #     ortho.f[:6, :6],
        #     ortho.g[:6, :6]
        # )
        #
        # fs, gs, us, uis = subdimer.run()
        # print("-" * 50)
        # print(f[:6, :6] - fs)
        # # print("-" * 50)
        # # print(ui)

    @validationTest
    def test_TwoState(self):

        w = 1500 * UnitsData.convert("Wavenumbers", "Hartrees")
        l = 100 * UnitsData.convert("Wavenumbers", "Hartrees")
        f = np.array([
            [w, l],
            [l, w]
        ])
        g = np.array([
            [w, 0],
            [0, w]
        ])

        """
        \textbf{F} = \left(\begin{tabular}{rrrrr}
3892.4 & -29.4 &   -0.0 &   -0.0 &  -59.9 \\
       & 176.2 &  -20.3 &  -20.3 &   78.0 \\
       &       & 3869.7 &  -65.3 &   -5.6 \\
       &       &        & 3869.7 &   -5.6 \\
       &       &        &        & 3735.3
       """

        ortho = ObliqueModeGenerator(f, g)
        f, g, u, ui = ortho.run()
        raise Exception(
            np.sqrt(scipy.linalg.eigvalsh(ortho.f, ortho.g, type=3)) * UnitsData.convert("Hartrees", "Wavenumbers"),
            np.linalg.eigvalsh(f) * UnitsData.convert("Hartrees", "Wavenumbers")
        )

    @validationTest
    def test_LocalModeComplement(self):

        mol = Molecule.from_file(
            TestManager.test_data('OCHH_freq.fchk')
        )

        nms = mol.get_normal_modes()
        locs = nms.localize(internals=[(0, 1), (1, 2), (1, 3)])
        comp = locs.get_complement()
        full_loc = nms.apply_transformation(
            np.concatenate(
                [
                    locs.localizing_transformation[0],
                    comp.localizing_transformation[0]
                ],
                axis=1
            )
        )

        import McUtils.Formatters as mfmt

        print(locs.local_freqs * UnitsData.hartrees_to_wavenumbers)
        print(comp.local_freqs * UnitsData.hartrees_to_wavenumbers)
        print(
            mfmt.TableFormatter("{:.0f}").format(
                full_loc.local_hessian * UnitsData.hartrees_to_wavenumbers
            )
        )

    @validationTest
    def test_OCHHObliqueVPT(self):
        ochh = Molecule.from_file(
            TestManager.test_data('OCHH_freq.fchk')
        )

        nms = NormalModes.from_molecule(ochh, dimensionless=False)
        # print(np.round(nms.matrix, 8))
        # print(np.round(ochh.normal_modes.modes.basis.matrix, 8))
        # print(phases)
        # raise Exception(...)

        make_oblique = True
        (tf, tf_inv), loc_nms = nms.get_localized_modes(
            [
                [1, 0],
                [1, 2],
                [1, 3],
                # [1, 2, 3],
                [1, 0, 2],
                [1, 0, 3],
                [2, 1, 0, 3]
            ],
            rediagonalize=True,
            make_oblique=make_oblique,
            mode_selection=[0, 1, 2],
            symmetrizer=lambda crds:np.array([
                crds[:, 0],
                1/np.sqrt(2)*(crds[:, 1] + crds[:, 2]),
                1/np.sqrt(2)*(crds[:, 1] - crds[:, 2]),
                1/np.sqrt(2)*(crds[:, 3] + crds[:, 4]),
                1/np.sqrt(2)*(crds[:, 3] - crds[:, 4]),
                crds[:, 5]
            ]).T
        )

        pot_expansion = list(ochh.potential_derivatives)
        nm_tf = tf.T
        pot_expansion[2] = np.tensordot(nm_tf, pot_expansion[2], axes=[1, 0])
        pot_expansion[3] = np.tensordot(nm_tf,
                                        np.tensordot(nm_tf, pot_expansion[3], axes=[1, 1]),
                                        axes=[1, 1]
                                        )

        from Psience.VPT2 import VPTRunner
        full_runner = VPTRunner.construct(ochh, 2,
                            degeneracy_specs='auto',
                            logger=False
                            )[0]
        subrunner = VPTRunner.construct(ochh, 2,
                            degeneracy_specs='auto',
                            mode_selection=[3, 4, 5],
                            logger=False
                            )[0]
        new_runner = VPTRunner.construct(
            [ochh.atoms, ochh.coords],
            2,
            degeneracy_specs='auto',
            modes={'freqs':loc_nms.freqs, 'matrix':loc_nms.matrix, 'inverse':loc_nms.inverse},
            logger=False,
            potential_derivatives=pot_expansion
        )[0]

        full_runner.print_tables(print_intensities=False)
        """
        >>------------------------- Degenerate Energies -------------------------
        State           Harmonic   Anharmonic     Harmonic   Anharmonic
                             ZPE          ZPE    Frequency    Frequency
        0 0 0 0 0 0   5866.87172   5796.34178            -            - 
        0 0 0 0 0 1            -            -   3061.70145   2987.29160 
        0 0 0 0 1 0            -            -   2977.64049   2843.35728 
        0 0 0 1 0 0            -            -   1727.08266   1693.95557 
        """
        # subrunner.print_tables(print_intensities=False)
        """
        >>------------------------- Degenerate Energies -------------------------
        State     Harmonic   Anharmonic     Harmonic   Anharmonic
                       ZPE          ZPE    Frequency    Frequency
        0 0 0   3883.21230   3848.71494            -            - 
        0 0 1            -            -   3061.70145   2924.82377 
        0 1 0            -            -   2977.64049   2847.67134 
        1 0 0            -            -   1727.08266   1711.71453 
        0 0 2            -            -   6123.40290   5804.53642 
        0 2 0            -            -   5955.28097   5604.60379 
        2 0 0            -            -   3454.16533   3408.09426 
        0 1 1            -            -   6039.34194   5641.53164 
        1 0 1            -            -   4788.78411   4639.61866 
        1 1 0            -            -   4704.72315   4556.23881 
        >>--------------------------------------------------<<
        """
        new_runner.print_tables(print_intensities=False)
        """
        State     Harmonic   Anharmonic     Harmonic   Anharmonic
               ZPE          ZPE    Frequency    Frequency
        0 0 0   3861.83994   3821.09271            -            - 
        0 0 1            -            -   3055.66055   2912.14354 
        0 1 0            -            -   2959.78211   2827.49408 
        1 0 0            -            -   1708.23723   1674.16514 
        0 0 2            -            -   6111.32109   5778.49325 
        0 2 0            -            -   5919.56421   5565.55617 
        2 0 0            -            -   3416.47445   3322.56151 
        0 1 1            -            -   6015.44265   5607.36174 
        1 0 1            -            -   4763.89777   4579.59631 
        1 1 0            -            -   4668.01933   4491.76492 
        >>--------------------------------------------------<<
        """

    @validationTest
    def test_WaterDimerObliqueVPT(self):
        dimer = Molecule.from_file(
            TestManager.test_data('water_dimer_freq.fchk')
        )

        nms = dimer.normal_modes.modes.basis.to_new_modes()#
        # nms = NormalModes.from_molecule(dimer, dimensionless=False)
        # from Psience.VPT2 import VPTRunner
        # VPTRunner.run_simple(dimer, 2,
        #                     degeneracy_specs='auto'
        #                     )

        make_oblique = True
        (tf, tf_inv), loc_nms = nms.get_localized_modes(
            (
                self.setupDefaultInternals(2, remapping=[1, 0, 2, 3, 4, 5])
                + [
                    [0, 3, 4],
                    [5, 3, 1]
                ]
            ),
            # fixed_coords=[1, 3],
            rediagonalize=True,
            make_oblique=make_oblique,
            mode_selection=[0, 1, 2]
            # symmetrizer=lambda crds: np.array([
            #     crds[:, 0],
            #     1 / np.sqrt(2) * (crds[:, 1] + crds[:, 2]),
            #     1 / np.sqrt(2) * (crds[:, 1] - crds[:, 2]),
            #     1 / np.sqrt(2) * (crds[:, 3] + crds[:, 4]),
            #     1 / np.sqrt(2) * (crds[:, 3] - crds[:, 4]),
            #     crds[:, 5]
            # ]).T
        )

        pot_expansion = list(dimer.potential_derivatives)
        nm_tf = tf.T
        # print(np.round(tf, 2))
        # print(np.round(tf_inv, 2))
        # print(np.round(nms.matrix - loc_nms.matrix @ tf_inv, 2))
        # print(np.round(nms.inverse - tf @ loc_nms.inverse, 2))
        # raise Exception(loc_nms.freqs, nms.freqs)
        pot_expansion[2] = np.tensordot(nm_tf, pot_expansion[2], axes=[1, 0])
        pot_expansion[3] = np.tensordot(nm_tf,
                                        np.tensordot(nm_tf, pot_expansion[3], axes=[1, 1]),
                                        axes=[1, 1]
                                        )

        from Psience.VPT2 import VPTRunner
        full_runner = VPTRunner.construct(dimer, 2,
                                          degeneracy_specs='auto',
                                          logger=False
                                          )[0]
        hf_runner = VPTRunner.construct(dimer, 2,
                                        degeneracy_specs='auto',
                                        mode_selection=[i for i,f in enumerate(nms.freqs) if f > 1000 / 219475.6],
                                        logger=False
                                        )[0]

        hoh_runner = VPTRunner.construct(dimer, 2,
                                        degeneracy_specs='auto',
                                        mode_selection=[-6 ,-4, -1],
                                        logger=False
                                        )[0]
        new_runner = VPTRunner.construct(
            [dimer.atoms, dimer.coords],
            2,
            degeneracy_specs='auto',
            # mode_selection=[i for i,f in enumerate(loc_nms.freqs) if f > 1000 / 219475.6],
            modes={'freqs': loc_nms.freqs, 'matrix': loc_nms.matrix, 'inverse': loc_nms.inverse},
            logger=False,
            potential_derivatives=pot_expansion
        )[0]
        # raise Exception(
        #     hf_runner.hamiltonian.V_terms[1][0] * 219475,
        #     new_runner.hamiltonian.V_terms[1][0] * 219475
        # )

        # full_runner.print_tables(print_intensities=False)
        """
        >>------------------------- Deperturbed Energies -------------------------
        :: State                       Harmonic   Anharmonic     Harmonic   Anharmonic
                                         ZPE          ZPE    Frequency    Frequency
        0 0 0 0 0 0 0 0 0 0 0 0  10133.32562   9928.34735            -            - 
        0 0 0 0 0 0 0 0 0 0 0 1            -            -   3935.27376   3779.83820 
        0 0 0 0 0 0 0 0 0 0 1 0            -            -   3915.15761   3741.22006 
        0 0 0 0 0 0 0 0 0 1 0 0            -            -   3813.91883   3647.33470 
        0 0 0 0 0 0 0 0 1 0 0 0            -            -   3718.74059   3582.50325 
        0 0 0 0 0 0 0 1 0 0 0 0            -            -   1650.27150   1615.63293 
        0 0 0 0 0 0 1 0 0 0 0 0            -            -   1629.31503   1585.61078 
        0 0 0 0 0 1 0 0 0 0 0 0            -            -    630.33864    570.09157 
        0 0 0 0 1 0 0 0 0 0 0 0            -            -    360.39517    331.79023 
        0 0 0 1 0 0 0 0 0 0 0 0            -            -    184.18194    218.02583 
        0 0 1 0 0 0 0 0 0 0 0 0            -            -    154.94037    205.30373 
        0 1 0 0 0 0 0 0 0 0 0 0            -            -    147.11065    169.95357 
        1 0 0 0 0 0 0 0 0 0 0 0            -            -    127.00715    167.93124 
        >>------------------------- Degenerate Energies -------------------------
        State                       Harmonic   Anharmonic     Harmonic   Anharmonic
                                         ZPE          ZPE    Frequency    Frequency
        0 0 0 0 0 0 0 0 0 0 0 0  10133.32562   9928.33507            -            - 
        0 0 0 0 0 0 0 0 0 0 0 1            -            -   3935.27376   3895.91045 
        0 0 0 0 0 0 0 0 0 0 1 0            -            -   3915.15761   3751.34005 
        0 0 0 0 0 0 0 0 0 1 0 0            -            -   3813.91883   3762.88143 
        0 0 0 0 0 0 0 0 1 0 0 0            -            -   3718.74059   3582.51554 
        0 0 0 0 0 0 0 1 0 0 0 0            -            -   1650.27150   1584.31960 
        0 0 0 0 0 0 1 0 0 0 0 0            -            -   1629.31503   1579.03776 
        0 0 0 0 0 1 0 0 0 0 0 0            -            -    630.33864    579.75883 
        0 0 0 0 1 0 0 0 0 0 0 0            -            -    360.39517    329.74797 
        0 0 0 1 0 0 0 0 0 0 0 0            -            -    184.18194    220.03342 
        0 0 1 0 0 0 0 0 0 0 0 0            -            -    154.94037    105.79500 
        0 1 0 0 0 0 0 0 0 0 0 0            -            -    147.11065    215.94945 
        1 0 0 0 0 0 0 0 0 0 0 0            -            -    127.00715    201.30696 
        """
        # hf_runner.print_tables(print_intensities=False)
        """
        >>------------------------- Deperturbed Energies -------------------------
        State           Harmonic   Anharmonic     Harmonic   Anharmonic
                             ZPE          ZPE    Frequency    Frequency
        0 0 0 0 0 0   9331.33866   9203.97083            -            - 
        0 0 0 0 0 1            -            -   3935.27375   3748.96957 
        0 0 0 0 1 0            -            -   3915.15760   3727.63027 
        0 0 0 1 0 0            -            -   3813.91883   3646.25947 
        0 0 1 0 0 0            -            -   3718.74059   3522.16962 
        0 1 0 0 0 0            -            -   1650.27151   1586.97884 
        1 0 0 0 0 0            -            -   1629.31502   1571.37950 
        0 0 0 0 0 2            -            -   7870.54750   7407.27914 
        0 0 0 0 2 0            -            -   7830.31521   7320.69281 
        0 0 0 2 0 0            -            -   7627.83767   7213.15140 
        0 0 2 0 0 0            -            -   7437.48119   6896.34749 
        0 2 0 0 0 0            -            -   3300.54303   3144.83427 
        2 0 0 0 0 0            -            -   3258.63005   3115.28710 
        0 0 0 0 1 1            -            -   7850.43135   7476.50529 
        0 0 0 1 0 1            -            -   7749.19258   7242.46255 
        0 0 1 0 0 1            -            -   7654.01434   7268.32213 
        0 1 0 0 0 1            -            -   5585.54526   5333.66574 
        1 0 0 0 0 1            -            -   5564.58877   5287.02147 
        0 0 0 1 1 0            -            -   7729.07644   7372.61204 
        0 0 1 0 1 0            -            -   7633.89820   7183.92249 
        0 1 0 0 1 0            -            -   5565.42912   5278.58862 
        1 0 0 0 1 0            -            -   5544.47263   5296.36071 
        0 0 1 1 0 0            -            -   7532.65943   7163.68204 
        0 1 0 1 0 0            -            -   5464.19035   5231.80428 
        1 0 0 1 0 0            -            -   5443.23386   5201.28059 
        0 1 1 0 0 0            -            -   5369.01211   5087.28544 
        1 0 1 0 0 0            -            -   5348.05562   5091.69518 
        1 1 0 0 0 0            -            -   3279.58654   3151.62003 
        >>--------------------------------------------------<<
        >>------------------------- Degenerate Energies -------------------------
        State           Harmonic   Anharmonic     Harmonic   Anharmonic
                             ZPE          ZPE    Frequency    Frequency
        0 0 0 0 0 0   9331.33866   9203.97083            -            - 
        0 0 0 0 0 1            -            -   3935.27375   3748.96957 
        0 0 0 0 1 0            -            -   3915.15760   3727.63027 
        0 0 0 1 0 0            -            -   3813.91883   3646.25947 
        0 0 1 0 0 0            -            -   3718.74059   3522.16962 
        0 1 0 0 0 0            -            -   1650.27151   1586.97884 
        1 0 0 0 0 0            -            -   1629.31502   1571.37950 
        0 0 0 0 0 2            -            -   7870.54750   7431.94954 
        0 0 0 0 2 0            -            -   7830.31521   7320.69281 
        0 0 0 2 0 0            -            -   7627.83767   7188.48099 
        0 0 2 0 0 0            -            -   7437.48119   6870.94252 
        0 2 0 0 0 0            -            -   3300.54303   3136.67763 
        2 0 0 0 0 0            -            -   3258.63005   3114.03466 
        0 0 0 0 1 1            -            -   7850.43135   7476.50529 
        0 0 0 1 0 1            -            -   7749.19258   7242.46255 
        0 0 1 0 0 1            -            -   7654.01434   7268.32213 
        0 1 0 0 0 1            -            -   5585.54526   5334.52274 
        1 0 0 0 0 1            -            -   5564.58877   5286.16447 
        0 0 0 1 1 0            -            -   7729.07644   7372.61204 
        0 0 1 0 1 0            -            -   7633.89820   7209.32746 
        0 1 0 0 1 0            -            -   5565.42912   5273.91429 
        1 0 0 0 1 0            -            -   5544.47263   5301.03504 
        0 0 1 1 0 0            -            -   7532.65943   7163.68204 
        0 1 0 1 0 0            -            -   5464.19035   5231.80428 
        1 0 0 1 0 0            -            -   5443.23386   5201.28059 
        0 1 1 0 0 0            -            -   5369.01211   5081.97346 
        1 0 1 0 0 0            -            -   5348.05562   5097.00716 
        1 1 0 0 0 0            -            -   3279.58654   3161.02911 
        >>--------------------------------------------------<<
        """
        # hoh_runner.print_tables(print_intensities=False)
        """
        0 0 0   4641.66468   4620.60458            -            - 
        0 0 1            -            -   3935.27375   3973.35109 
        0 1 0            -            -   3718.74059   3573.65580 
        1 0 0            -            -   1629.31502   1584.61633 
        0 0 2            -            -   7870.54750   8030.90518 
        0 2 0            -            -   7437.48119   7014.34155 
        2 0 0            -            -   3258.63005   3158.09074 
        0 1 1            -            -   7654.01434   7522.32324 
        1 0 1            -            -   5564.58877   5490.39975 
        1 1 0            -            -   5348.05562   5158.72630 
        """
        """[1634.66514245 3718.35800754 3914.2578522 ]"""
        new_runner.print_tables(print_intensities=False)
        """
        Oblique - Fixed Oxygens
        0 0 0   4605.90346   4542.81409            -            - 
        0 0 1            -            -   3894.26511   3720.46361 
        0 1 0            -            -   3678.41124   3478.55113 
        1 0 0            -            -   1639.13057   1584.63395 
        0 0 2            -            -   7788.53022   7276.36721 
        0 2 0            -            -   7356.82247   6760.76489 
        2 0 0            -            -   3278.26113   3137.96333 
        0 1 1            -            -   7572.67635   7209.44253 
        1 0 1            -            -   5533.39568   5287.85701 
        1 1 0            -            -   5317.54180   5034.04152 
        Unfixed Oxygens
        0 0 0   4540.91903   4482.78604            -            - 
        0 0 1            -            -   3869.09698   3706.50784 
        0 1 0            -            -   3568.31191   3370.52545 
        1 0 0            -            -   1644.42916   1594.51342 
        0 0 2            -            -   7738.19396   7219.57674 
        0 2 0            -            -   7136.62383   6524.71866 
        2 0 0            -            -   3288.85833   3155.78147 
        0 1 1            -            -   7437.40889   7143.09928 
        1 0 1            -            -   5513.52614   5296.65491 
        1 1 0            -            -   5212.74108   4936.06448 
        """
        """
        Local - Fixed Oxygens
        0 0 0   4606.28678   4543.28703            -            - 
        0 0 1            -            -   3894.44212   3720.47509 
        0 1 0            -            -   3678.64903   3479.17402 
        1 0 0            -            -   1639.48241   1585.19859 
        0 0 2            -            -   7788.88424   7276.44092 
        0 2 0            -            -   7357.29805   6762.33886 
        2 0 0            -            -   3278.96482   3139.18788 
        0 1 1            -            -   7573.09115   7209.80003 
        1 0 1            -            -   5533.92453   5288.27712 
        1 1 0            -            -   5318.13143   5035.62013  
        Unfixed Oxygens
        0 0 0   4586.98380   4525.18051            -            - 
        0 0 1            -            -   3883.65169   3713.29168 
        0 1 0            -            -   3645.78138   3448.02059 
        1 0 0            -            -   1644.53454   1591.67466 
        0 0 2            -            -   7767.30339   7246.93650 
        0 2 0            -            -   7291.56276   6700.87118 
        2 0 0            -            -   3289.06907   3150.98654 
        0 1 1            -            -   7529.43307   7188.50544 
        1 0 1            -            -   5528.18623   5296.34689 
        1 1 0            -            -   5290.31592   5007.32050 
        """

    @validationTest
    def test_WaterTrimerObliqueVPT(self):
        trimer = Molecule.from_file(
            TestManager.test_data('water_trimer_freq.fchk')
        )

        nms = trimer.normal_modes.modes.basis.to_new_modes()  #

        # print(nms.origin)
        # print(nms.matrix)
        # raise Exception(...)

        # nms = NormalModes.from_molecule(dimer, dimensionless=False)
        # from Psience.VPT2 import VPTRunner
        # VPTRunner.run_simple(dimer, 2,
        #                     degeneracy_specs='auto'
        #                     )

        make_oblique = True
        (tf, tf_inv), loc_nms = nms.get_localized_modes(
            [
                [0, 3],
                [0, 4],
                [0, 3, 4],

                [1, 5],
                [1, 6],
                [1, 5, 6],

                [2, 7],
                [2, 8],
                [2, 7, 8],

                [0, 1],
                [0, 2],
                [1, 2],

                [3, 0, 1, 2],
                [4, 0, 1, 2],
                [5, 1, 2, 0],
                [6, 1, 2, 0],
                [7, 2, 0, 1],
                [8, 2, 0, 1]
            ],
            rediagonalize=True,
            make_oblique=make_oblique,
            mode_selection=[0, 1, 2]
        )

        # print(loc_nms.matrix)
        # raise Exception(...)

        pot_expansion = list(trimer.potential_derivatives)
        nm_tf = tf.T
        pot_expansion[2] = np.tensordot(nm_tf, pot_expansion[2], axes=[1, 0])
        pot_expansion[3] = np.tensordot(nm_tf,
                                        np.tensordot(nm_tf, pot_expansion[3], axes=[1, 1]),
                                        axes=[1, 1]
                                        )

        from Psience.VPT2 import VPTRunner
        # full_runner = VPTRunner.construct(dimer, 2,
        #                                   degeneracy_specs='auto',
        #                                   logger=False
        #                                   )[0]
        # hf_runner = VPTRunner.construct(dimer, 2,
        #                                 degeneracy_specs='auto',
        #                                 mode_selection=[i for i, f in enumerate(nms.freqs) if f > 1000 / 219475.6],
        #                                 logger=False
        #                                 )[0]

        # hoh_runner = VPTRunner.construct(dimer, 2,
        #                                  degeneracy_specs='auto',
        #                                  mode_selection=[-6, -4, -1],
        #                                  logger=False
        #                                  )[0]
        new_runner = VPTRunner.construct(
            [trimer.atoms, trimer.coords],
            2,
            degeneracy_specs='auto',
            # mode_selection=[i for i,f in enumerate(loc_nms.freqs) if f > 1000 / 219475.6],
            modes={'freqs': loc_nms.freqs, 'matrix': loc_nms.matrix, 'inverse': loc_nms.inverse},
            logger=False,
            potential_derivatives=pot_expansion
        )[0]

        new_runner.print_tables(print_intensities=False)
        """
        Local
        >>------------------------- States Energies -------------------------
        State     Harmonic   Anharmonic     Harmonic   Anharmonic
                       ZPE          ZPE    Frequency    Frequency
        0 0 0   4491.20675   4420.53846            -            - 
        0 0 1            -            -   3830.55603   3665.56065 
        0 1 0            -            -   3557.98659   3266.90514 
        1 0 0            -            -   1593.87087   1575.27283 
        0 0 2            -            -   7661.11207   7150.36935 
        0 2 0            -            -   7115.97319   6252.61703 
        2 0 0            -            -   3187.74174   3123.52550 
        0 1 1            -            -   7388.54263   6929.91203 
        1 0 1            -            -   5424.42691   5274.90038 
        1 1 0            -            -   5151.85747   4824.95528 
        >>--------------------------------------------------<<
        """
        """
        Oblique
        >>------------------------- States Energies -------------------------
        State     Harmonic   Anharmonic     Harmonic   Anharmonic
                       ZPE          ZPE    Frequency    Frequency
        0 0 0   4421.44047   4354.49515            -            - 
        0 0 1            -            -   3783.69417   3627.53420 
        0 1 0            -            -   3492.24781   3198.11735 
        1 0 0            -            -   1566.93896   1551.67921 
        0 0 2            -            -   7567.38833   7058.66029 
        0 2 0            -            -   6984.49562   6098.99773 
        2 0 0            -            -   3133.87792   3071.62672 
        0 1 1            -            -   7275.94198   6852.53424 
        1 0 1            -            -   5350.63312   5232.82697 
        1 1 0            -            -   5059.18677   4729.12687 
        >>--------------------------------------------------<<
        """

    def runMBPolModel(self,
                      model_name,
                      base_coords,
                      internals,
                      coord_sel=None,
                      mode_selection=None,
                      monomer=0,
                      states=2,
                      logger=False,
                      degeneracy_specs='auto',
                      fix_oxygens=False,
                      mass_weight=True,
                      unitary=True,
                      make_oblique=True,
                      reoptimize=False,
                      print_normal_modes=False,
                      print_local_modes=False,
                      print_target_modes=False,
                      reexpand=False,
                      return_runner=False
                      ):
        from Psience.VPT2 import VPTRunner

        loader = ModuleLoader(TestManager.current_manager().test_data_dir)
        mbpol = loader.load("LegacyMBPol").potential

        base_coords = np.asanyarray(base_coords)
        nwaters = base_coords.shape[0] // 3

        expansion_file = os.path.expanduser(f'~/Documents/Postdoc/local_VPT/{model_name}.hdf5')
        with Checkpointer.from_file(expansion_file) as chk:
            try:
                if reoptimize: raise Exception("reoptimizing")
                coords = chk['coords']
            except Exception as e:
                print("Reoptimizing because:", e)
                import scipy.optimize as opt

                min_opt = opt.minimize(lambda x: mbpol(x, nwaters=nwaters)[0], base_coords.flatten(),
                                       tol=1e-8
                                       )
                print(min_opt.x.reshape(3*nwaters, 3))
                print(min_opt.fun * 219475.6)
                print(min_opt.message)
                val_base, grad_base = mbpol(base_coords, deriv_order=1, nwaters=nwaters)
                print(val_base * 219475.6)
                print(grad_base * 219475.6)
                val_final, grad_opt = mbpol(min_opt.x.reshape(3*nwaters, 3), deriv_order=1, nwaters=nwaters)
                print(grad_opt * 219475.6)

                coords = min_opt.x
                chk['coords'] = coords
        coords = np.asanyarray(coords).reshape((3*nwaters, 3))

        try:
            if reexpand: raise Exception("reexpand")
            expansion = chk['force_field']
        except Exception as e:
            print("Reexpanding because:", e)
            expansion = mbpol(coords, nwaters=nwaters, deriv_order=4)
            chk['force_field'] = expansion
            # raise Exception(expansion[2].shape)
        expansion = [np.array(x) for x in expansion]

        mol = Molecule(
            ["O", "H", "H"] * nwaters,
            coords,
            potential_derivatives=expansion[1:]
        )

        nms = NormalModes.from_molecule(mol)
        if print_normal_modes:
            with np.printoptions(threshold=1e8):
                print(nms.origin)
                print(nms.matrix)

        if coord_sel is None and monomer is not None:
            if isinstance(monomer, (int, np.integer)):
                coord_sel = [0 + 3 * monomer, 1 + 3 * monomer, 2 + 3 * monomer]
            else:
                monomers = list(sorted(monomer)) # type:list
                coord_sel = sum([
                    [0 + 3 * m, 1 + 3 * m, 2 + 3 * m]
                    for m in monomers
                    ],
                    []
                )
                # include the OO coordinates
                o_pos = nwaters * 3
                for i,m1 in enumerate(monomers):
                    shift_1 = sum(nwaters-j-1 for j in range(m1))
                    for m2 in monomers[i+1:]:
                        oo_mode = o_pos + shift_1 + (m2 - m1) - 1
                        coord_sel.append(oo_mode)
                        # print(internals[oo_mode])

        if internals is not None:
            (tf, tf_inv), loc_nms, target_modes = nms.get_localized_modes(
                internals,
                rediagonalize=True,
                fixed_coords=[3*i for i in range(nwaters)] if fix_oxygens else None,
                mass_weight=mass_weight,
                unitary=unitary,
                make_oblique=make_oblique,
                mode_selection=coord_sel,
                return_target_modes=True
            )

            if print_target_modes:
                with np.printoptions(threshold=1e8):
                    print(target_modes)
        else:
            loc_nms = nms

        if print_local_modes:
            with np.printoptions(threshold=1e8):
                # if mass_weight:
                #     print(loc_nms.make_mass_weighted(loc_nms.masses / 1832).matrix)
                # else:
                print(loc_nms.matrix)

        pot_expansion = list(mol.potential_derivatives)
        pot_expansion = [x.reshape((nwaters * 9,) * (n + 1)) for n, x in enumerate(pot_expansion)]
        # raise Exception([p.shape for p in tet.potential_derivatives])
        # nm_tf = tf.T
        # raise Exception(nm_tf.shape, [x.shape for x in pot_expansion])
        # pot_expansion[2] = np.tensordot(nm_tf, pot_expansion[2], axes=[1, 0])
        # pot_expansion[3] = np.tensordot(nm_tf,
        #                                 np.tensordot(nm_tf, pot_expansion[3], axes=[1, 1]),
        #                                 axes=[1, 1]
        #                                 )

        if isinstance(mode_selection, (int, float, np.integer, np.floating)):
             mode_selection = [i for i,f in enumerate(loc_nms.freqs) if f > mode_selection]

        new_runner = VPTRunner.construct(
            [mol.atoms, mol.coords],
            states,
            degeneracy_specs=degeneracy_specs,
            mode_selection=mode_selection,
            modes={'freqs': loc_nms.freqs, 'matrix': loc_nms.matrix, 'inverse': loc_nms.inverse},
            logger=logger,
            potential_derivatives=pot_expansion
        )[0]
        if return_runner:
            return new_runner
        else:
            new_runner.print_tables(print_intensities=False)
    def runMBPolMonomers(self,
                         model_name,
                         base_coords,
                         internals,
                         print_local_modes=True,
                         **opts
                         ):
        nwaters = len(base_coords) // 3
        for i in range(nwaters):
            print("=" * 25, "Monomer", i, "=" * 25)
            self.runMBPolModel(
                model_name,
                base_coords,
                internals,
                **opts,
                print_local_modes=print_local_modes,
                monomer=i
            )
    def runMBPolDimers(self,
                       model_name,
                       base_coords,
                       internals,
                       print_local_modes=True,
                       **opts
                       ):
        nwaters = len(base_coords) // 3
        for i in itertools.combinations(range(nwaters), 2):
            print("=" * 25, "Dimer", i, "=" * 25)
            self.runMBPolModel(
                model_name,
                base_coords,
                internals,
                print_local_modes=print_local_modes,
                **opts,
                monomer=list(i)
            )

    def runMBPolNMers(self,
                      model_name,
                      base_coords,
                      internals,
                      k,
                      print_local_modes=True,
                      **opts
                      ):
        nwaters = len(base_coords) // 3
        for i in itertools.combinations(range(nwaters), k):
            print("=" * 25, "NMer", i, "=" * 25)
            self.runMBPolModel(
                model_name,
                base_coords,
                internals,
                print_local_modes=print_local_modes,
                **opts,
                monomer=list(i)
            )

    def runMBE(self,
               model_name,
               base_coords,
               internals,
               order,
               run_non_oblique=True,
               **opts
               ):

        if isinstance(order, (int, np.integer)):
            order = range(1, order+1)

        for i,k in enumerate(order):
            if i > 0:
                print("~"*100)
            if k == 1:
                runner = self.runMBPolMonomers
            elif k == 2:
                runner = self.runMBPolDimers
            else:
                runner = lambda *a,_o=k,**kw:self.runMBPolNMers(*a,_o,**kw)

            runner(
                model_name,
                base_coords,
                internals,
                **opts
            )
            if run_non_oblique:
                runner(
                    model_name,
                    base_coords,
                    internals,
                    make_oblique=False,
                    **opts
                )

    def setupDefaultInternals(self, nwaters, remapping=None,
                              include_hoh_rocks=False,
                              include_hoh_oop=True,
                              include_hydrogen_oop=False,
                              include_oxygen_oop=False):
        internals = []
        for i in range(nwaters):
            # monomer internals
            monomer = [
                np.array([0, 1]) + 3*i,
                np.array([0, 2]) + 3*i,
                np.array([0, 1, 2]) + 3*i
            ]
            internals.extend([x.tolist() for x in monomer])

        # OO stretches
        waters_pos = np.arange(nwaters) * 3
        combos = itertools.combinations(waters_pos, 2)
        internals.extend(combos)

        if include_hydrogen_oop:
            if nwaters > 2:
                nats = nwaters * 3
                for i in range(nwaters):
                    diheds = np.array([
                        [1, 0, 3, 6],
                        [2, 0, 3, 6]
                    ]) + 3*i
                    internals.extend((diheds % nats).tolist())
            else:
                internals.extend([
                    [1, 0, 3, 2],
                    [4, 3, 0, 5]
                ])

        if include_hoh_rocks:
            for i in range(nwaters):
                # monomer internals
                internals.append(
                    {"rock":np.array([0, 1, 2]) + 3 * i}
                )
        if include_hoh_oop:
            for i in range(nwaters):
                # monomer internals
                internals.extend([
                    {"oop":np.array([1, 0, 2]) + 3 * i},
                    {"oop":np.array([2, 0, 1]) + 3 * i}
                ])

        if include_oxygen_oop:
            if nwaters > 3:
                nats = nwaters * 3
                for i in range(nwaters):
                    diheds = np.array([
                        [0, 3, 6, 9]
                    ]) + 3 * i
                    internals.extend((diheds % nats).tolist())

        if remapping is not None:
            internals = [ [remapping[i] for i in idx] for idx in internals ]

        # for i in internals: print(i)
        # raise Exception(len(internals))
        return internals

    @validationTest
    def test_WaterDimerMBPol(self):
        dimer_coords = [
             [ 2.86771112e+00, -3.17080506e-05, -2.23401093e-01],
             [ 3.63213828e+00, -3.73743134e-04,  1.41524327e+00],
             [ 1.07451468e+00,  3.17123593e-05,  9.76384644e-02],
             [-2.63391749e+00, -1.95872838e-04,  2.00041963e-01],
             [-3.31628543e+00,  1.43506723e+00, -6.69046768e-01],
             [-3.31598121e+00, -1.43556339e+00, -6.69113632e-01]
        ]
        dimer_internals = self.setupDefaultInternals(2,
                                                     include_hoh_rocks=True,
                                                     # include_hydrogen_oop=True,
                                                     # include_oxygen_oop=True
                                                     ) + [
            # [0, 2, 3]
            # {'rock':[0, 2, 1]}
            # [3, 4, 0],
            # [3, 5, 0],
            # [2, 0, 3, 1],
            # [5, 3, 0, 4]
        ]
        # for i in dimer_internals: print(i)
        # coord_sel = list(range(len(dimer_internals)))
        coord_sel = None
        unitary = True

        # self.runMBPolModel(
        #     'dimer',
        #     dimer_coords,
        #     None,
        #     fix_oxygens=False,
        #     mass_weight=False,
        #     make_oblique=False,
        #     subsel=None,
        #     print_normal_modes=True
        # )
        """
        0 0 0 0 0 0 0 0 0 0 0 0  10144.77584   9923.11324            -            - 
        0 0 0 0 0 0 0 0 0 0 0 1            -            -   3934.80518   3741.99045 
        0 0 0 0 0 0 0 0 0 0 1 0            -            -   3915.37781   3728.90061 
        0 0 0 0 0 0 0 0 0 1 0 0            -            -   3832.01218   3647.22314 
        0 0 0 0 0 0 0 0 1 0 0 0            -            -   3751.90382   3596.50725 
        0 0 0 0 0 0 0 1 0 0 0 0            -            -   1664.90489   1610.51239 
        0 0 0 0 0 0 1 0 0 0 0 0            -            -   1650.35778   1600.56668 
        0 0 0 0 0 1 0 0 0 0 0 0            -            -    602.06023    493.96275 
        0 0 0 0 1 0 0 0 0 0 0 0            -            -    353.61723    338.86730 
        0 0 0 1 0 0 0 0 0 0 0 0            -            -    183.09187    150.34762 
        0 0 1 0 0 0 0 0 0 0 0 0            -            -    146.65736    135.63871 
        0 1 0 0 0 0 0 0 0 0 0 0            -            -    139.43585    122.53374 
        1 0 0 0 0 0 0 0 0 0 0 0            -            -    115.32749     98.06471 
        """
        # dimer_internals = dimer_internals[:7]

        self.runMBPolDimers(
            'dimer',
            dimer_coords,
            dimer_internals,
            make_oblique=True,
            # print_target_modes=True,
            print_local_modes=True,
            states=1,
            unitary=unitary,
            coord_sel=coord_sel
        )
        self.runMBPolDimers(
            'dimer',
            dimer_coords,
            dimer_internals,
            make_oblique=False,
            print_local_modes=True,
            states=1,
            unitary=unitary,
            # logger=True,
            # degeneracy_specs=None
            coord_sel=coord_sel
        )
        raise Exception(...)
        """
          0 0 0 0 0 0 0 0 0 0 0 0   5853.10002  -6074.76262
          0 0 0 0 0 0 0 0 0 0 0 1   4855.74021  -5270.21753
          0 0 0 0 0 0 0 0 0 0 1 0   4892.10490  -5300.24470
          0 0 0 0 0 0 0 0 0 1 0 0   4885.68015  -5292.13179
          0 0 0 0 0 0 0 0 1 0 0 0   5266.83748  -5643.89664
          0 0 0 0 0 0 0 1 0 0 0 0   5971.60496  -6247.66005
          0 0 0 0 0 0 1 0 0 0 0 0   5954.60280  -6226.05650
          0 0 0 0 0 1 0 0 0 0 0 0   6844.30755  -7174.06762
          0 0 0 0 1 0 0 0 0 0 0 0   8054.93504  -8291.34756
          0 0 0 1 0 0 0 0 0 0 0 0   7741.35679  -7995.76363
          0 0 1 0 0 0 0 0 0 0 0 0  10720.29560 -10952.97684
          0 1 0 0 0 0 0 0 0 0 0 0  12571.55265 -12810.11736
          1 0 0 0 0 0 0 0 0 0 0 0  15945.16875 -16184.09412
          """
        # print("?"*100)
        # print(
        #     #43259.11917984811 -8466.495366325204
        #     #75789.83266522287 -7139.485937973492
        #     np.max(runner.hamiltonian.V_terms[2] * 219475.6),
        #     np.min(runner.hamiltonian.V_terms[2] * 219475.6)
        # )
        # print("=" * 25, "Monomer 01", "=" * 25)
        # self.runMBPolModel(
        #     'dimer',
        #     dimer_coords,
        #     dimer_internals,
        #     make_oblique=False,
        #     print_local_modes=True,
        #     unitary=False,
        #     # coord_sel=[i for i,x in enumerate(dimer_internals) if any(j in x for j in [0, 1, 2])],
        #     monomer=[0, 1],
        #     mode_selection=1000 / UnitsData.hartrees_to_wavenumbers#
        # )
        raise Exception(...)
        for i in range(2):
            print("="*25, "Monomer", i, "="*25)
            self.runMBPolModel(
                'dimer',
                dimer_coords,
                dimer_internals,
                fix_oxygens=False,
                mass_weight=True,
                print_local_modes=True,
                print_target_modes=False,
                monomer=i
            )
        for i in range(2):
            print("="*25, "Monomer", i, "="*25)
            self.runMBPolModel(
                'dimer',
                dimer_coords,
                dimer_internals,
                make_oblique=False,
                print_local_modes=True,
                fix_oxygens=False,
                mass_weight=True,
                monomer=i
            )

    @validationTest
    def test_WaterTrimerMBPol(self):
        init_coords = np.array(
            [[-2.84907277, -1.09595083, -0.14633201],
             [-1.26117558, -2.00130798, -0.04933503],
             [-3.93857322, -1.87593528,  1.06738756],
             [ 2.3830576 , -1.90840511, -0.15181051],
             [ 2.32690962, -0.08075627, -0.04173429],
             [ 3.63762974, -2.43398881,  1.0390116 ],
             [ 0.46194949,  3.00918878,  0.19210532],
             [-1.09663158,  2.05078279,  0.10653446],
             [ 0.34950138,  4.23739783, -1.13107505]]
        )

        nwaters = init_coords.shape[0] // 3
        internals = self.setupDefaultInternals(nwaters)

        # self.runMBPolModel(
        #     'trimer',
        #     init_coords,
        #     internals,
        #     coord_sel=None,
        #     monomer=None
        # )
        # raise Exception(...)

        # self.runMBPolNMers(
        #     'trimer',
        #     init_coords,
        #     internals,
        #     3
        # )
        # self.runMBPolNMers(
        #     'trimer',
        #     init_coords,
        #     internals,
        #     3,
        #     make_oblique=False,
        # )
        # raise Exception(...)

        # self.runMBPolDimers(
        #     'trimer',
        #     init_coords,
        #     internals,
        # )
        # self.runMBPolDimers(
        #     'trimer',
        #     init_coords,
        #     internals,
        #     make_oblique=False,
        # )
        # raise Exception(...)

        self.runMBPolMonomers(
            'trimer',
            init_coords,
            internals,
        )
        self.runMBPolMonomers(
            'trimer',
            init_coords,
            internals,
            make_oblique=False
        )

    @validationTest
    def test_WaterTetramerMBPol(self):

        from Psience.VPT2 import VPTRunner

        loader = ModuleLoader(TestManager.current_manager().test_data_dir)
        mbpol = loader.load("LegacyMBPol").potential

        optimize_struct = False
        if optimize_struct:
            import scipy.optimize as opt
            # opt_ints = opt.minimize(pot_i, ints, method='nelder-mead', tol=1e-8)
            # raise Exception(opt_ints.fun, opt_ints.message)
            base_atoms = ["O", "H", "H", "O", "H", "H", "O", "H", "H"]
            base_coords = np.array([
                 [-1.92362674, -0.17446612, -0.00722898 ],
                 [-2.49650201, -0.36265251, -0.75533246 ],
                 [-1.26537422, -0.90289372,  0.00571245 ],
                 [ 0.17446612, -1.92362674,  0.00722898 ],
                 [ 0.90289372, -1.26537422, -0.00571245 ],
                 [ 0.36265251, -2.49650201,  0.75533246 ],
                 [-0.17446612,  1.92362674,  0.00722898 ],
                 [-0.36265251,  2.49650201,  0.75533246 ],
                 [-0.90289372,  1.26537422, -0.00571245 ],
                 [ 1.92362674,  0.17446612, -0.00722898 ],
                 [ 1.26537422,  0.90289372,  0.00571245 ],
                 [ 2.49650201,  0.36265251, -0.75533246 ]
            ]) * UnitsData.convert("Angstroms", "BohrRadius")

            print(base_coords)
            min_opt = opt.minimize(lambda x:mbpol(x, nwaters=4)[0], base_coords.flatten(),
                                   tol=1e-8
                                   )
            print(min_opt.x.reshape(12, 3))
            print(min_opt.fun * 219475.6)
            print(min_opt.message)
            val_base, grad_base = mbpol(base_coords, deriv_order=1, nwaters=4)
            print(val_base * 219475.6)
            print(grad_base * 219475.6)
            val_final, grad_opt = mbpol(min_opt.x.reshape(12, 3), deriv_order=1, nwaters=4)
            print(grad_opt * 219475.6)

            raise Exception(...)

        base_struct = np.array(
            [[-3.65820826, -0.33541594, -0.08484091],
             [-4.71931702, -0.6970077,  -1.5036615],
             [-2.43608816, -1.71114672, -0.01353123],
             [ 0.33625233, -3.65983196,  0.08037927],
             [ 1.71188582, -2.43757664,  0.00961531],
             [ 0.69821123, -4.72183449,  1.49844432],
             [-0.33391615,  3.65905136,  0.08479249],
             [-0.697188,    4.72024229,  1.50312788],
             [-1.7095794,   2.43694241,  0.01187531],
             [ 3.66093091,  0.33481905, -0.08199767],
             [ 2.43885954,  1.71052585, -0.00966898],
             [ 4.72268458,  0.69798994, -1.49994101]]
        )

        reexpand = False
        expansion_file = os.path.expanduser('~/Documents/Postdoc/local_VPT/tetramer.json')
        with Checkpointer.from_file(expansion_file) as chk:
            try:
                if reexpand: raise Exception("reexpand")
                expansion = chk['force_field']
            except Exception as e:
                print("Reexpanding because:", e)
                expansion = mbpol(base_struct, nwaters=4, deriv_order=4)
                chk['force_field'] = expansion
                # raise Exception(expansion[2].shape)
            expansion = [np.array(x) for x in expansion]

        tet = Molecule(
            ["O", "H", "H"] * 4,
            base_struct,
            potential_derivatives=expansion[1:]
        )

        nms = NormalModes.from_molecule(tet)
        # with np.printoptions(threshold=1e8):
        #     print(nms.origin)

        make_oblique = True
        (tf, tf_inv), loc_nms = nms.get_localized_modes(
            [
                [0, 1],
                [0, 2],
                [1, 0, 2],

                [3, 4],
                [3, 5],
                [4, 3, 5],

                [6, 7],
                [6, 8],
                [7, 6, 8],

                [9, 10],
                [9, 11],
                [10, 9, 11],

                [0, 3],
                [0, 6],
                [0, 9],
                [3, 6],
                [3, 9],
                [6, 9],

                [1, 0, 3, 6],
                [2, 0, 3, 6],
                [4, 3, 6, 9],
                [5, 3, 6, 9],
                [7, 6, 9, 0],
                [8, 6, 9, 0],
                [10, 9, 0, 6],
                [11, 9, 0, 6]

            ],
            rediagonalize=True,
            make_oblique=make_oblique,
            mode_selection=[0, 1, 2]
        )

        # print(loc_nms.matrix)

        pot_expansion = list(tet.potential_derivatives)
        pot_expansion = [x.reshape((36,) * (n+1)) for n,x in enumerate(pot_expansion)]
        # raise Exception([p.shape for p in tet.potential_derivatives])
        # nm_tf = tf.T
        # raise Exception(nm_tf.shape, [x.shape for x in pot_expansion])
        # pot_expansion[2] = np.tensordot(nm_tf, pot_expansion[2], axes=[1, 0])
        # pot_expansion[3] = np.tensordot(nm_tf,
        #                                 np.tensordot(nm_tf, pot_expansion[3], axes=[1, 1]),
        #                                 axes=[1, 1]
        #                                 )

        new_runner = VPTRunner.construct(
            [tet.atoms, tet.coords],
            2,
            degeneracy_specs='auto',
            # mode_selection=[i for i,f in enumerate(loc_nms.freqs) if f > 1000 / 219475.6],
            modes={'freqs': loc_nms.freqs, 'matrix': loc_nms.matrix, 'inverse': loc_nms.inverse},
            logger=False,
            potential_derivatives=pot_expansion
        )[0]
        new_runner.print_tables(print_intensities=False)

        """
        >>------------------------- States Energies -------------------------
        State     Harmonic   Anharmonic     Harmonic   Anharmonic
                       ZPE          ZPE    Frequency    Frequency
        0 0 0   4347.63550   4278.19960            -            - 
        0 0 1            -            -   3817.43749   3639.22905 
        0 1 0            -            -   3409.54283   3192.21494 
        1 0 0            -            -   1468.29068   1409.40140 
        0 0 2            -            -   7634.87498   7107.97125 
        0 2 0            -            -   6819.08566   6188.52147 
        2 0 0            -            -   2936.58135   2772.00850 
        0 1 1            -            -   7226.98032   6814.39790 
        1 0 1            -            -   5285.72816   5050.23334 
        1 1 0            -            -   4877.83350   4575.82347 
        >>--------------------------------------------------<<
        """

    @validationTest
    def test_WaterTetramerAutoMBPol(self):

        tetramer_coords = np.array(
            [[-3.65820826, -0.33541594, -0.08484091],
             [-4.71931702, -0.6970077,  -1.5036615],
             [-2.43608816, -1.71114672, -0.01353123],
             [ 0.33625233, -3.65983196,  0.08037927],
             [ 1.71188582, -2.43757664,  0.00961531],
             [ 0.69821123, -4.72183449,  1.49844432],
             [-0.33391615,  3.65905136,  0.08479249],
             [-0.697188,    4.72024229,  1.50312788],
             [-1.7095794,   2.43694241,  0.01187531],
             [ 3.66093091,  0.33481905, -0.08199767],
             [ 2.43885954,  1.71052585, -0.00966898],
             [ 4.72268458,  0.69798994, -1.49994101]]
        )
        tetramer_internals = self.setupDefaultInternals(4)

        for i in range(4):
            print("="*25, "Monomer", i, "="*25)
            self.runMBPolModel(
                'tetramer',
                tetramer_coords,
                tetramer_internals,
                # print_normal_modes=True,
                fix_oxygens=True,
                print_local_modes=True,
                make_oblique=True,
                monomer=i
            )
        """
        Fixed Oxygens
        State     Harmonic   Anharmonic     Harmonic   Anharmonic
               ZPE          ZPE    Frequency    Frequency
        0 0 0   4416.12196   4341.36653            -            - 
        0 0 1            -            -   3867.44051   3689.29570 
        0 1 0            -            -   3499.56458   3253.93961 
        1 0 0            -            -   1465.23883   1400.27635 
        0 0 2            -            -   7734.88102   7210.54022 
        0 2 0            -            -   6999.12917   6280.18982 
        2 0 0            -            -   2930.47765   2758.24158 
        0 1 1            -            -   7367.00510   6937.85746 
        1 0 1            -            -   5332.67934   5074.76260 
        1 1 0            -            -   4964.80341   4623.72268 
        :: State     Harmonic   Anharmonic     Harmonic   Anharmonic
                       ZPE          ZPE    Frequency    Frequency
        0 0 0   4409.23264   4337.66754            -            - 
        0 0 1            -            -   3865.64674   3685.81833 
        0 1 0            -            -   3500.74444   3278.54904 
        1 0 0            -            -   1452.07409   1391.12498 
        0 0 2            -            -   7731.29349   7203.38967 
        0 2 0            -            -   7001.48887   6329.54676 
        2 0 0            -            -   2904.14819   2722.70830 
        0 1 1            -            -   7366.39118   6959.54929 
        1 0 1            -            -   5317.72084   5058.59851 
        1 1 0            -            -   4952.81853   4685.20391 
        >>------------------------- States Energies -------------------------
        :: State     Harmonic   Anharmonic     Harmonic   Anharmonic
                       ZPE          ZPE    Frequency    Frequency
        0 0 0   4396.53055   4330.67292            -            - 
        0 0 1            -            -   3867.86616   3689.21376 
        0 1 0            -            -   3500.12522   3297.30679 
        1 0 0            -            -   1425.06972   1381.99528 
        0 0 2            -            -   7735.73233   7210.54449 
        0 2 0            -            -   7000.25043   6366.05456 
        2 0 0            -            -   2850.13944   2700.71669 
        0 1 1            -            -   7367.99138   6981.29233 
        1 0 1            -            -   5292.93588   5054.89854 
        1 1 0            -            -   4925.19493   4736.01146 
        >>------------------------- Degenerate Energies -------------------------
        :: State     Harmonic   Anharmonic     Harmonic   Anharmonic
                       ZPE          ZPE    Frequency    Frequency
        0 0 0   4388.11394   4313.99718            -            - 
        0 0 1            -            -   3865.47856   3685.09229 
        0 1 0            -            -   3501.14913   3268.84032 
        1 0 0            -            -   1409.60020   1343.27414 
        0 0 2            -            -   7730.95713   7201.03235 
        0 2 0            -            -   7002.29826   6310.55717 
        2 0 0            -            -   2819.20040   2631.89138 
        0 1 1            -            -   7366.62769   6949.18239 
        1 0 1            -            -   5275.07876   5010.64858 
        1 1 0            -            -   4910.74933   4606.49400 
        3 1 0            -            -   7729.94972   7117.83065 
        """
        """
        Unfixed Oxygens
        >>------------------------- States Energies -------------------------
        State     Harmonic   Anharmonic     Harmonic   Anharmonic
                       ZPE          ZPE    Frequency    Frequency
        0 0 0   4347.82130   4278.37603            -            - 
        0 0 1            -            -   3817.40794   3639.20924 
        0 1 0            -            -   3409.60121   3192.25666 
        1 0 0            -            -   1468.63345   1409.71425 
        0 0 2            -            -   7634.81588   7107.93419 
        0 2 0            -            -   6819.20243   6188.60057 
        2 0 0            -            -   2937.26690   2772.60713 
        0 1 1            -            -   7227.00915   6814.41753 
        1 0 1            -            -   5286.04139   5050.54305 
        1 1 0            -            -   4878.23467   4576.15568 
        >>--------------------------------------------------<<
        >>------------------------- States Energies -------------------------
        State     Harmonic   Anharmonic     Harmonic   Anharmonic
                       ZPE          ZPE    Frequency    Frequency
        0 0 0   4337.19153   4267.02868            -            - 
        0 0 1            -            -   3813.15665   3633.08785 
        0 1 0            -            -   3432.51692   3230.87923 
        1 0 0            -            -   1428.70949   1360.85540 
        0 0 2            -            -   7626.31331   7096.00616 
        0 2 0            -            -   6865.03384   6263.07440 
        2 0 0            -            -   2857.41898   2650.16046 
        0 1 1            -            -   7245.67358   6847.41795 
        1 0 1            -            -   5241.86614   4990.69388 
        1 1 0            -            -   4861.22641   4602.37649 
        >>--------------------------------------------------<<
        >>------------------------- States Energies -------------------------
        :: State     Harmonic   Anharmonic     Harmonic   Anharmonic
                       ZPE          ZPE    Frequency    Frequency
        0 0 0   4333.97011   4272.77122            -            - 
        0 0 1            -            -   3813.57114   3639.40713 
        0 1 0            -            -   3425.73183   3255.57140 
        1 0 0            -            -   1428.63727   1385.66943 
        0 0 2            -            -   7627.14227   7115.04417 
        0 2 0            -            -   6851.46365   6312.99322 
        2 0 0            -            -   2857.27453   2701.40701 
        0 1 1            -            -   7239.30296   6885.60976 
        1 0 1            -            -   5242.20840   5013.65749 
        1 1 0            -            -   4854.36909   4706.58793 
        >>--------------------------------------------------<<
        >>------------------------- States Energies -------------------------
        :: State     Harmonic   Anharmonic     Harmonic   Anharmonic
                       ZPE          ZPE    Frequency    Frequency
        0 0 0   4316.77780   4243.77467            -            - 
        0 0 1            -            -   3814.38088   3634.01963 
        0 1 0            -            -   3423.40714   3209.60377 
        1 0 0            -            -   1395.76759   1324.39463 
        0 0 2            -            -   7628.76176   7094.73003 
        0 2 0            -            -   6846.81428   6221.43853 
        2 0 0            -            -   2791.53518   2581.70824 
        0 1 1            -            -   7237.78802   6824.82898 
        1 0 1            -            -   5210.14847   4963.10465 
        1 1 0            -            -   4819.17473   4520.72412 
        >>--------------------------------------------------<<
        """
        for i in range(4):
            print("="*25, "Monomer", i, "="*25)
            self.runMBPolModel(
                'tetramer',
                tetramer_coords,
                tetramer_internals,
                # print_normal_modes=True,
                fix_oxygens=True,
                # print_local_modes=True,
                make_oblique=False,
                monomer=i
            )
        """
        Fixed Oxygens
        >>------------------------- States Energies -------------------------
        State     Harmonic   Anharmonic     Harmonic   Anharmonic
                       ZPE          ZPE    Frequency    Frequency
        0 0 0   4446.74294   4372.96863            -            - 
        0 0 1            -            -   3867.89856   3683.14800 
        0 1 0            -            -   3500.63979   3253.62322 
        1 0 0            -            -   1524.94754   1460.37533 
        0 0 2            -            -   7735.79712   7198.37237 
        0 2 0            -            -   7001.27957   6280.17439 
        2 0 0            -            -   3049.89507   2887.24797 
        0 1 1            -            -   7368.53835   6931.06927 
        1 0 1            -            -   5392.84610   5115.57141 
        1 1 0            -            -   5025.58732   4679.81147 
        >>------------------------- State Energies -------------------------
        State     Harmonic   Anharmonic     Harmonic   Anharmonic
                       ZPE          ZPE    Frequency    Frequency
        0 0 0   4446.74932   4376.22866            -            - 
        0 0 1            -            -   3866.84095   3684.38243 
        0 1 0            -            -   3501.36902   3265.49258 
        1 0 0            -            -   1525.28868   1467.88434 
        0 0 2            -            -   7733.68190   7201.44552 
        0 2 0            -            -   7002.73804   6365.03519 
        2 0 0            -            -   3050.57735   2896.93899 
        0 1 1            -            -   7368.20997   6944.31790 
        1 0 1            -            -   5392.12963   5127.54553 
        1 1 0            -            -   5026.65769   4732.47364 
        2 1 0            -            -   6551.94637   6076.64163 
        3 0 0            -            -   4575.86603   4275.63920 
        >>------------------------- Degenerate Energies -------------------------
        State     Harmonic   Anharmonic     Harmonic   Anharmonic
                       ZPE          ZPE    Frequency    Frequency
        0 0 0   4439.87334   4371.95775            -            - 
        0 0 1            -            -   3868.12762   3687.29421 
        0 1 0            -            -   3501.06103   3271.52553 
        1 0 0            -            -   1510.55803   1461.46560 
        0 0 2            -            -   7736.25524   7206.76006 
        0 2 0            -            -   7002.12206   6405.87638 
        2 0 0            -            -   3021.11607   2883.59044 
        0 1 1            -            -   7369.18865   6953.13493 
        1 0 1            -            -   5378.68565   5128.43452 
        1 1 0            -            -   5011.61907   4752.20873 
        2 1 0            -            -   6522.17710   6066.83046 
        3 0 0            -            -   4531.67410   4247.97885 
        >>------------------------- States Energies -------------------------
        State     Harmonic   Anharmonic     Harmonic   Anharmonic
                       ZPE          ZPE    Frequency    Frequency
        0 0 0   4430.99934   4359.39333            -            - 
        0 0 1            -            -   3866.74336   3682.57054 
        0 1 0            -            -   3501.63488   3260.34331 
        1 0 0            -            -   1493.62044   1434.40626 
        0 0 2            -            -   7733.48673   7197.39727 
        0 2 0            -            -   7003.26977   6294.07736 
        2 0 0            -            -   2987.24088   2835.26956 
        0 1 1            -            -   7368.37825   6937.47374 
        1 0 1            -            -   5360.36380   5089.55887 
        1 1 0            -            -   4995.25532   4670.82504 
        """

        """
        Unfixed Oxygens
        >>------------------------- States Energies -------------------------
        State     Harmonic   Anharmonic     Harmonic   Anharmonic
                       ZPE          ZPE    Frequency    Frequency
        0 0 0   4435.62685   4363.51521            -            - 
        0 0 1            -            -   3857.62035   3672.69228 
        0 1 0            -            -   3480.17927   3242.78225 
        1 0 0            -            -   1533.45408   1469.69736 
        0 0 2            -            -   7715.24069   7176.42403 
        0 2 0            -            -   6960.35854   6268.41369 
        2 0 0            -            -   3066.90815   2904.20984 
        0 1 1            -            -   7337.79962   6907.83262 
        1 0 1            -            -   5391.07442   5118.09650 
        1 1 0            -            -   5013.63335   4679.62910 
        >>------------------------- Degenerate Energies -------------------------
        State     Harmonic   Anharmonic     Harmonic   Anharmonic
                       ZPE          ZPE    Frequency    Frequency
        0 0 0   4433.48851   4364.18840            -            - 
        0 0 1            -            -   3855.79219   3672.29672 
        0 1 0            -            -   3492.65590   3264.46934 
        1 0 0            -            -   1518.52894   1460.86671 
        0 0 2            -            -   7711.58437   7177.11374 
        0 2 0            -            -   6985.31180   6363.60260 
        2 0 0            -            -   3037.05787   2882.53023 
        0 1 1            -            -   7348.44809   6929.03099 
        1 0 1            -            -   5374.32112   5108.86700 
        1 1 0            -            -   5011.18484   4722.21170 
        2 1 0            -            -   6529.71377   6069.08419 
        3 0 0            -            -   4555.58681   4255.49326 
        >>------------------------- Degenerate Energies -------------------------
        State     Harmonic   Anharmonic     Harmonic   Anharmonic
                       ZPE          ZPE    Frequency    Frequency
        0 0 0   4432.04218   4365.41115            -            - 
        0 0 1            -            -   3855.73960   3674.51604 
        0 1 0            -            -   3488.49961   3269.16821 
        1 0 0            -            -   1519.84516   1470.52287 
        0 0 2            -            -   7711.47919   7181.53183 
        0 2 0            -            -   6976.99922   6418.06951 
        2 0 0            -            -   3039.69031   2899.32243 
        0 1 1            -            -   7344.23921   6936.15192 
        1 0 1            -            -   5375.58475   5125.12464 
        1 1 0            -            -   5008.34477   4761.61316 
        2 1 0            -            -   6528.18993   6080.26683 
        3 0 0            -            -   4559.53547   4269.19292 
        >>------------------------- States Energies -------------------------
        State     Harmonic   Anharmonic     Harmonic   Anharmonic
                       ZPE          ZPE    Frequency    Frequency
        0 0 0   4419.44530   4348.96610            -            - 
        0 0 1            -            -   3856.37422   3670.90921 
        0 1 0            -            -   3487.85125   3254.42126 
        1 0 0            -            -   1494.66514   1435.09334 
        0 0 2            -            -   7712.74844   7173.45020 
        0 2 0            -            -   6975.70249   6291.17293 
        2 0 0            -            -   2989.33028   2836.00424 
        0 1 1            -            -   7344.22547   6917.86264 
        1 0 1            -            -   5351.03936   5079.27678 
        1 1 0            -            -   4982.51639   4665.46163 
        """

    @validationTest#(log_file='~/Documents/Postdoc/local_VPT/hexamer_mbe_rocks.txt')
    def test_WaterHexamerMBPol(self):
        hexamer_base_coords = np.array([
            [-0.35921654,  0.66032868,  1.70403349],
            [ 0.48387738,  0.24478184,  1.46606149],
            [-1.00096128,  0.17370798,  1.16214939],
            [ 1.77202874, -0.59659445,  0.12135318],
            [ 1.35538225, -0.06078792, -0.60578423],
            [ 2.72228751, -0.49953465,  0.01340080],
            [ 0.04288025,  2.88312174,  0.12951739],
            [-0.54206536,  3.63113699,  0.27293631],
            [-0.17328769,  2.23111583,  0.82893297],
            [ 0.31403336,  0.76377541, -1.64191878],
            [ 0.18199259,  1.63209002, -1.21518887],
            [-0.51692626,  0.29334277, -1.46887566],
            [ 0.13037902, -2.79310762,  0.00956184],
            [ 0.16772631, -3.31339836,  0.81700174],
            [ 0.86584103, -2.15058446,  0.08135296],
            [-1.72933868, -0.85118903, -0.37267206],
            [-1.20686995, -1.67534742, -0.25579395],
            [-2.62103595, -1.13064966, -0.59766950]
        ]) * UnitsData.convert("Angstroms", "BohrRadius")
        hexamer_internals = self.setupDefaultInternals(6,
                                                       include_hoh_rocks=True,
                                                       include_hoh_oop=True)
        self.runMBPolMonomers(
            'hexamer',
            hexamer_base_coords,
            hexamer_internals,
            print_local_modes=False
        )
        raise Exception(...)
        # for i in [3, 4, 5, 9, 10, 11, 24, 34, 36]:
        #     print(hexamer_internals[i])
        # raise Exception(...)

        self.runMBPolModel(
            'hexamer',
            hexamer_base_coords,
            hexamer_internals,
            coord_sel=[3, 4, 5, 34],
            monomer=None,
            # monomer=[1, 3],
            unitary=False,
            make_oblique=True,
            # return_runner=True
        )
        self.runMBPolModel(
            'hexamer',
            hexamer_base_coords,
            hexamer_internals,
            coord_sel=[9, 10, 11, 36],
            monomer=None,
            # monomer=[1, 3],
            unitary=False,
            make_oblique=True,
            # return_runner=True
        )
        self.runMBPolModel(
            'hexamer',
            hexamer_base_coords,
            hexamer_internals,
            coord_sel=[3, 4, 5, 34],
            monomer=None,
            # monomer=[1, 3],
            unitary=False,
            make_oblique=False,
            # return_runner=True
        )
        self.runMBPolModel(
            'hexamer',
            hexamer_base_coords,
            hexamer_internals,
            coord_sel=[9, 10, 11, 36],
            monomer=None,
            # monomer=[1, 3],
            unitary=False,
            make_oblique=False,
            # return_runner=True
        )
        self.runMBPolModel(
            'hexamer',
            hexamer_base_coords,
            hexamer_internals,
            coord_sel=[3, 4, 5, 9, 10, 11, 24, 34, 36],
            monomer=None,
            # monomer=[1, 3],
            unitary=False,
            make_oblique=True,
            # return_runner=True
        )
        self.runMBPolModel(
            'hexamer',
            hexamer_base_coords,
            hexamer_internals,
            coord_sel=[3, 4, 5, 9, 10, 11, 24, 34, 36],
            monomer=None,
            # monomer=[1, 3],
            unitary=False,
            make_oblique=False,
            # return_runner=True
        )
        # t1 = wtf_runner.hamiltonian.V_terms[1]
        # wtf_runner = self.runMBPolModel(
        #     'hexamer',
        #     hexamer_base_coords,
        #     hexamer_internals,
        #     monomer=[0, 4],
        #     unitary=False,
        #     make_oblique=True,
        #     return_runner=True
        # )
        # t2 = wtf_runner.hamiltonian.V_terms[1]
        #
        # print(...)
        # print(t1[-1]*219475.6)
        # print(t2[-1]*219475.6)
        # raise Exception(...)

        # self.runMBE(
        #     'hexamer',
        #     hexamer_base_coords,
        #     hexamer_internals,
        #     2,
        #     unitary=True
        # )

    @validationTest
    def test_WaterOctamerMBPol(self):
        init_coords = np.array([[1.36660032,1.36660032,1.318721],
                                [1.51174383,0.42459088,1.51490506],
                                [0.42459088,1.51174383,1.51490506],
                                [-1.46840116,-1.46840116,-1.34546438],
                                [-2.10311105,-2.10311105,-1.68929956],
                                [-1.52349463,-1.52349463,-0.35604814],
                                [1.36660032,-1.36660032,-1.318721],
                                [1.51174383,-0.42459088,-1.51490506],
                                [0.42459088,-1.51174383,-1.51490506],
                                [-1.46840116,1.46840116,1.34546438],
                                [-1.52349463,1.52349463,0.35604814],
                                [-2.10311105,2.10311105,1.68929956],
                                [-1.36660032,-1.36660032,1.318721],
                                [-0.42459088,-1.51174383,1.51490506],
                                [-1.51174383,-0.42459088,1.51490506],
                                [-1.36660032,1.36660032,-1.318721],
                                [-1.51174383,0.42459088,-1.51490506],
                                [-0.42459088,1.51174383,-1.51490506],
                                [1.46840116,1.46840116,-1.34546438],
                                [1.52349463,1.52349463,-0.35604814],
                                [2.10311105,2.10311105,-1.68929956],
                                [1.46840116,-1.46840116,1.34546438],
                                [1.52349463,-1.52349463,0.35604814],
                                [2.10311105,-2.10311105,1.68929956]]) * UnitsData.convert("Angstroms", "BohrRadius")

        nwaters = init_coords.shape[0] // 3
        internals = self.setupDefaultInternals(nwaters)

        for i in range(nwaters):
            print("=" * 25, "Monomer", i, "=" * 25)
            self.runMBPolModel(
                'water_{}'.format(nwaters),
                init_coords,
                internals,
                # print_normal_modes=i==0,
                print_local_modes=True,
                make_oblique=True,
                monomer=i
            )
        for i in range(nwaters):
            print("=" * 25, "Monomer", i, "=" * 25)
            self.runMBPolModel(
                'water_{}'.format(nwaters),
                init_coords,
                internals,
                # print_normal_modes=i==0,
                print_local_modes=True,
                make_oblique=False,
                monomer=i
            )

    @debugTest
    def test_MassLocalization(self):

        mol = Molecule.from_file(
            TestManager.test_data('proplybenz.hess'),
        )

        loc0 = mol.get_normal_modes().localize(
            atoms=[18]
        )#.make_oblique()
        print(loc0.local_freqs * UnitsData.hartrees_to_wavenumbers)
        loc1 = mol.get_normal_modes().localize(
            mass_scaling=[1/12, 1/12],
            atoms=[8, 18]
        )#.make_oblique()
        print(loc1.local_freqs * UnitsData.hartrees_to_wavenumbers)
        # print(loc0.make_mass_weighted().coords_by_modes @ loc0.make_mass_weighted().modes_by_coords)
        print(
            np.round(
                loc0.make_mass_weighted().coords_by_modes @
                loc1.make_mass_weighted().modes_by_coords,
                2
            )
        )

        return

    @validationTest
    def test_Localizers(self):

        mol = Molecule.from_file(
            TestManager.test_data('HOONO_freq.fchk'),
        )

        fig = mol.plot(backend='x3d')
        fig.show()

        raise Exception(...)

        # trms = mol.translation_rotation_modes[1][0]
        # base_F = mol.potential_derivatives[1]
        # base_G = mol.g_matrix
        # g12 = scipy.linalg.fractional_matrix_power(base_G, 1/2)
        # g12i = scipy.linalg.fractional_matrix_power(base_G, -1/2)
        # proj = (np.eye(trms.shape[0]) - trms @ trms.T) @ g12
        # # base_modes
        # # print(trms.shape)
        # proj_F = proj @ base_F @ proj.T
        # print(proj_F)
        # # raise Exception(
        # #     mol.normal_modes.modes.freqs**2,
        # #     np.linalg.eigh(proj_F)[0][6:]
        # # )
        # # print(base_G)
        # # print(base_F)
        # print("="*50)

        base_modes = mol.normal_modes.modes.basis.to_new_modes()
        new_modes = LocalizedModes.localize_by_masses(base_modes, [0])
        final_freqs = np.array(
            # [base_modes.freqs*219475.6]
            [n.freqs*219475.6 for n in new_modes]
        )
        print(final_freqs.T.tolist())