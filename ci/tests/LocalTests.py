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

from Psience.LocalModes import *
from Psience.Molecools import *

import numpy as np

class LocalTests(TestCase):

    def setUp(self) -> None:
        np.set_printoptions(linewidth=1e8)

    @debugTest
    def test_Water(self):
        mol = Molecule.from_file(
            TestManager.test_data('HOD_freq.fchk'),
            internals=[
                [0, -1, -1, -1],
                [1,  0, -1, -1],
                [2,  0,  1, -1]
            ]
        )

        iterative = False
        if iterative:
            ortho = BlockLocalFGOrthogonalizer.from_molecule(mol, logger=True)
            f, g, u, its, _ = ortho.run()
            ui = np.linalg.inv(u)

            rl, a, rr = np.linalg.svd(u)
            r = (rl @ rr)
            # with np.printoptions(linewidth=1e8, suppress=True):
            #     print(np.round(r, 8))

            u2 = u @ r.T
            u2i = np.linalg.inv(u2)
            print(np.round(u2, 8))
            print(np.round(u2i, 8))
            # print(u2)
            # print("-"*50)
            #
            # print(u.T@ ortho.g @u)
            # print(ui@ ortho.f @ui.T)
            # print("-"*50)
            gp = u2.T@ ortho.g @u2
            fp = u2i@ ortho.f @u2i.T
            # print(gp)
            # print(fp)

            freq2, modes = scipy.linalg.eigh(ortho.f, ortho.g, type=3)
            freqs = np.sqrt(freq2)
            modes = modes / np.sqrt(freqs)[np.newaxis, :]
            rU, a, rM = np.linalg.svd(modes)
            R = (rU @ rM)
            # with np.printoptions(linewidth=1e8, suppress=True):
            #     print(np.round(r, 8))

            no_rot = modes @ R.T
            print("-"*100)
            print(no_rot)
            print("-"*100)
            print(modes.T @ ortho.f @ modes)
            print("-"*100)
            print(u2i)

            raise Exception(
                freqs,
                # modes.T @ ortho.f @ modes,
                modes @ R.T
            )


            # print(u2)
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

            # with np.printoptions(linewidth=1e8):
            #     print("Iterations:", its)
            #     print(f)
            #     print("-" * 50)
            #     print(g)

            self.assertLess(its, 250)

        else:
            ortho = BlockLocalFGOrthogonalizer.from_molecule(mol, frequency_scaled=True)

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
        if iterative:
            ortho = BlockLocalFGOrthogonalizerIterative.from_molecule(mol, logger=True)
            f, g, u, its, _ = ortho.run()

            # with np.printoptions(linewidth=1e8):
            #     print("Iterations:", its)
            #     print(f)
            #     print("-" * 50)
            #     print(g)

            rl, a, rr = np.linalg.svd(u)
            r = (rl @ rr)
            # with np.printoptions(linewidth=1e8, suppress=True):
            #     print(np.round(r, 8))

            u2 = u @ r.T
            u2i = np.linalg.inv(u2)
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
                "U:", np.linalg.norm(np.diag(f - g)),
                "S:", np.linalg.norm(np.diag(fp - gp)),
            )

            self.assertLess(its, 250)

        else:
            ortho = BlockLocalFGOrthogonalizer.from_molecule(mol, sel=[0, 1, 3, 2, 4])
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

        iterative = False
        if iterative:
            ortho = BlockLocalFGOrthogonalizerIterative.from_molecule(mol, logger=True)
            f, g, u, its, _ = ortho.run()

            # with np.printoptions(linewidth=1e8):
            #     print("Iterations:", its)
            #     print(f)
            #     print("-" * 50)
            #     print(g)

            rl, a, rr = np.linalg.svd(u)
            r = (rl @ rr)
            # with np.printoptions(linewidth=1e8, suppress=True):
            #     print(np.round(r, 8))

            u2 = u @ r.T
            u2i = np.linalg.inv(u2)
            # print(u2)
            # print("-" * 50)

            # ui = np.linalg.inv(u)
            # print(u.T @ ortho.g @ u)
            # print(ui @ ortho.f @ ui.T)
            # print("-" * 50)
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
                "U:",np.linalg.norm(np.diag(f - g)),
                "S:",np.linalg.norm(np.diag(fp - gp)),
            )

            self.assertLess(its, 250)

        else:
            ortho = BlockLocalFGOrthogonalizer.from_molecule(mol)
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
        ortho = BlockLocalFGOrthogonalizer.from_molecule(mol)
        f, g, u, ui = ortho.run()

        self.assertLess(abs(np.linalg.norm((f - g).flatten())), 1e-15)

        iterative = False
        if iterative:
            ortho = BlockLocalFGOrthogonalizerIterative.from_molecule(mol, logger=True)
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
            ortho = BlockLocalFGOrthogonalizer.from_molecule(mol, sel=[0, 1, 3, 2, 4, 5])
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

            ortho = BlockLocalFGOrthogonalizer.from_molecule(mol, sel=[0, 1, 3, 2, 4, 5])
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
            ortho = BlockLocalFGOrthogonalizerIterative.from_molecule(dimer,
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
            ortho = BlockLocalFGOrthogonalizer.from_molecule(dimer, sel=sel)
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

            ortho = BlockLocalFGOrthogonalizer.from_molecule(dimer, sel=sel)
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

            # subdimer = BlockLocalFGOrthogonalizer(
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

        ortho = BlockLocalFGOrthogonalizer(f, g)
        f, g, u, ui = ortho.run()
        raise Exception(
            np.sqrt(scipy.linalg.eigvalsh(ortho.f, ortho.g, type=3)) * UnitsData.convert("Hartrees", "Wavenumbers"),
            np.linalg.eigvalsh(f) * UnitsData.convert("Hartrees", "Wavenumbers")
        )





