import collections
import gc
import itertools
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

        iterative = False
        if iterative:
            ortho = ObliqueModeGeneratorIterative.from_molecule(mol, logger=True)
            f, g, u, its, _ = ortho.run()

            # with np.printoptions(linewidth=1e8):
            #     print('Initial:')
            #     print("F")
            #     print(ortho.f)
            #     print("-" * 50)
            #     print("G")
            #     print(ortho.g)
            #     print("="*50)
            #     print("Iterations:", its)
            #     print("F")
            #     print(f)
            #     print("-" * 50)
            #     print("G")
            #     print(g)

            rl, a, rr = np.linalg.svd(u)
            r = (rl @ rr)
            # with np.printoptions(linewidth=1e8, suppress=True):
            #     print(np.round(r, 8))

            u2 = u @ r.T
            u2i = np.linalg.inv(u2)
            # print(u2)
            # print("-"*50)
            #
            # print(u.T@ ortho.g @u)
            # print(ui@ ortho.f @ui.T)
            # print("-"*50)
            gp = u2.T @ ortho.g @ u2
            fp = u2i @ ortho.f @ u2i.T
            # print(gp)
            # print(fp)

            print(u2)
            print(
                "Asymmetry Norms",
                "U", np.linalg.norm((u - u.T).flatten()),
                "S", np.linalg.norm((u2 - u2.T).flatten())
            )
            print(
                "Off-Diag Norms  ",
                "U:", ortho.off_diag_norm(f, g),
                "S:", ortho.off_diag_norm(fp, gp),
            )
            print(
                "Diag Norms  ",
                "U:", np.linalg.norm(np.diag(f - g)),
                "S:", np.linalg.norm(np.diag(fp - gp)),
            )

            self.assertLess(its, 250)

        else:

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

        iterative = False
        if iterative:
            ortho = ObliqueModeGeneratorIterative.from_molecule(dimer,
                                                             max_iterations=5000,
                                                             scaling=1,
                                                             damping=1.001,
                                                             tolerance=1e-14
                                                             )
            f, g, u, its, conv = ortho.run()

            rl, a, rr = np.linalg.svd(u)
            r = (rl @ rr)
            u2 = u @ r.T
            u2i = np.linalg.inv(u2)
            # print(u2)
            # print("-" * 50)

            ui = np.linalg.inv(u)
            # print(u.T @ ortho.g @ u)
            # print(ui @ ortho.f @ ui.T)
            # print("-" * 50)
            gp = u2.T @ ortho.g @ u2
            fp = u2i @ ortho.f @ u2i.T


            print(u2)
            print(
                "Asymmetry Norms",
                "U", np.linalg.norm((u - u.T).flatten()),
                "S", np.linalg.norm((u2 - u2.T).flatten())
            )
            print(
                "Off-Diag Norms  ",
                "U:", ortho.off_diag_norm(f, g),
                "S:", ortho.off_diag_norm(fp, gp),
            )
            print(
                "Diag Norms  ",
                "U:",np.linalg.norm(np.diag(f - g)),
                "S:",np.linalg.norm(np.diag(fp - gp)),
            )

            # plt.Plot(np.arange(len(conv)), np.array(conv),
            #          linewidth=1
            #          ).show()
            #
            # with np.printoptions(linewidth=1e8):
            #     print(f)
            #     print("-" * 50)
            #     print(g)

            self.assertLess(its, 5000)

        else:
            ortho = ObliqueModeGenerator.from_molecule(dimer, sel=sel)
            f, g, u, ui = ortho.run()
            # print("-"*50)
            # print(f)
            # print("-"*50)
            # print(u)
            # print("-"*50)
            # print(ui)

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

    @debugTest
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
                [0, 1, 2],
                [0, 1, 3],
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

        make_oblique = False
        (tf, tf_inv), loc_nms = nms.get_localized_modes(
            [
                [1, 2], # LHS SP OH
                [0, 1], # free OH
                [0, 1, 2], # Free bend
                [0, 1, 3, 2], # OOP SP

                [1, 3], # OO
                [4, 3],
                [5, 3],
                [4, 3, 1],
                [5, 3, 1],
                [5, 3, 1, 4],
                [4, 3, 1, 2],
                [4, 3, 1, 0]
            ],
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

        # print(dimer.coords)
        # print(nms.matrix)
        # print(loc_nms.matrix)
        # raise Exception(...)

        # raise Exception(loc_nms.freqs * 219475.6)
        # masses = dimer.atomic_masses
        # undim_modes = nms.make_dimensionless(masses)
        # undim_modes = loc_nms.make_dimensionless(masses)
        # trip_mass = np.broadcast_to(
        #     np.sqrt(masses)[:, np.newaxis],
        #     (len(masses), undim_modes.matrix.shape[0] // len(masses))
        # ).flatten()
        # mw_der = dimer.potential_derivatives[1] / (trip_mass[:, np.newaxis] * trip_mass[np.newaxis, :])
        # raise Exception(
        #     undim_modes.freqs * 219475.6,
        #     np.diag(undim_modes.matrix.T @ undim_modes.matrix) * 219475.6,
        #     # undim_modes.matrix.T @ undim_modes.matrix * 219475.6,
        #     np.diag(undim_modes.inverse @ mw_der @ undim_modes.inverse.T) * 219475.6,
        # )

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
        hoh_runner.print_tables(print_intensities=False)
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
        new_runner.print_tables(print_intensities=False)
        """
        Oblique
        0 0 0   4543.47464   4485.06828            -            - 
        0 0 1            -            -   3865.65444   3698.22906 
        0 1 0            -            -   3574.53344   3376.02145 
        1 0 0            -            -   1646.76140   1599.13137 
        0 0 2            -            -   7731.30887   7201.12091 
        0 2 0            -            -   7149.06688   6534.81736 
        2 0 0            -            -   3293.52280   3164.80781 
        0 1 1            -            -   7440.18788   7135.05098 
        1 0 1            -            -   5512.41584   5292.38363 
        1 1 0            -            -   5221.29484   4951.77943 
        """
        """
        Local
        0 0 1            -            -   3881.93936   3706.82985 
        0 1 0            -            -   3649.28873   3454.91359 
        1 0 0            -            -   1647.81095   1595.20414 
        0 0 2            -            -   7763.87872   7234.00043 
        0 2 0            -            -   7298.57745   6714.99623 
        2 0 0            -            -   3295.62190   3158.06919 
        0 1 1            -            -   7531.22809   7187.01672 
        1 0 1            -            -   5529.75031   5285.86021 
        1 1 0            -            -   5297.09968   5025.75608 
        >>--------------------------------------------------<<
        """

