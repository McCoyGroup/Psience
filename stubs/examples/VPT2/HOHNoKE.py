"""Extracted from VPT2Tests.test_HOHNoKE via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VPT2Tests.test_HOHNoKE"""

import itertools
import tempfile
try:
    from Peeves.TestUtils import *
    from Peeves import BlockProfiler
except:
    pass
from unittest import TestCase
from Psience.VPT2 import *
from Psience.Molecools import Molecule
from Psience.BasisReps import HarmonicOscillatorProductBasis, BasisStateSpace
from McUtils.Data import UnitsData
import McUtils.Plots as plt
import McUtils.Numputils as nput
from McUtils.Scaffolding import *
from McUtils.Profilers import Timer
from McUtils.Parallelizers import Parallelizer, SerialNonParallelizer, MultiprocessingParallelizer
from McUtils.Zachary import FiniteDifferenceDerivative
import sys, os, numpy as np, itertools as ip

class VPT2Tests(TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        np.set_printoptions(linewidth=int(100000000.0))

    @inactiveTest
    def test_HOHNoKE(self):
        print()
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
        dimer_internals = [[LHF, X, X, X], [LO, LHF, X, X], [SH, LO, LHF, X], [RH2, SH, LO, LHF], [RO, LO, RH2, LHF], [RH1, RO, RH2, LHF]]
        hoono_internals = [[1, -1, -1, -1], [2, 1, -1, -1], [3, 2, 1, -1], [0, 1, 2, 3], [4, 3, 2, 1]]
        test_internals = None
        file_name = 'HOH_freq.fchk'
        runner1, _ = VPTRunner.construct(TestManager.test_data(file_name), 2, internals=test_internals, logger=False)
        G_new = runner1.hamiltonian.G_terms.base_terms.optimize_coordinates()
        raise Exception([G.shape for G in G_new])
        runner2, _ = VPTRunner.construct(TestManager.test_data(file_name), 2, internals=[[1, -1, -1, -1], [2, 1, -1, -1], [0, 1, 2, -1]], logger=False)
        from McUtils.Zachary import TensorDerivativeConverter

        def symm_terms(v_terms):
            test_V = [0, 0, 0, 0]
            for i in [3, 2, 1]:
                test_V[i] = v_terms[i - 1]
            test_V[0] = np.zeros(len(test_V[1]))
            n = len(test_V[0])
            for i in range(n):
                for j in range(n):
                    for k in range(n):
                        perms = list(itertools.permutations([i, j, k]))
                        v3 = sum((test_V[2][p] for p in perms)) / len(perms)
                        for p in perms:
                            test_V[2][p] = v3
                        for l in range(n):
                            perms = list(itertools.permutations([i, j, k, l]))
                            v4 = sum((test_V[3][p] for p in perms)) / len(perms)
                            for p in perms:
                                test_V[3][p] = v4
            return test_V
        test_V1 = symm_terms(runner1.hamiltonian.V_terms)
        n = len(test_V1[0])
        test_G = [np.zeros((n,) * (i + 2)) if isinstance(g, int) and g == 0 else g for i, g in enumerate(reversed([runner1.hamiltonian.G_terms[n] for n in [2, 1, 0]]))]
        V3 = test_V1[2]
        freqs = np.diag(test_V1[1])
        V3 = np.tensordot(V3, np.tensordot(V3, np.diag(1 / freqs), axes=[2, -1]), axes=[2, -1])
        V3 = V3 + np.moveaxis(V3, 2, 1) + np.moveaxis(V3, 2, 0)
        ders = [np.eye(len(test_V1[1])), -test_V1[2] / (3 * freqs[np.newaxis, np.newaxis, :]), -(test_V1[3] - 5 / 9 * V3) / (4 * freqs[np.newaxis, np.newaxis, np.newaxis, :]), 0]
        convs = TensorDerivativeConverter.convert_fast(ders, test_V1, order=len(ders))
        R2 = test_V1[2] / (3 * freqs[np.newaxis, np.newaxis, :])
        R3 = (test_V1[3] - 1 / 9 * V3) / (4 * freqs[np.newaxis, np.newaxis, np.newaxis, :])
        dersR = [np.eye(len(test_V1[1])), R2, R3, 0]
        G0 = np.diag(freqs)
        Q1 = ders[0]
        Q2 = ders[1]
        Q3 = ders[2]
        R1 = dersR[0]
        G1_Q = np.tensordot(R2, G0, axes=[1, 0])
        W = freqs[:, np.newaxis] / freqs[np.newaxis, :]
        V1 = test_V1[2]
        VV = test_V1[2][:, :, :, np.newaxis, np.newaxis] * test_V1[2][np.newaxis, np.newaxis, :, :, :]
        G1 = test_G[1] + test_V1[2] / 3 * (W + W.T)[np.newaxis]
        G2_QR = 2 * np.tensordot(ders[1], test_V1[2] / 3, axes=[2, 0])
        G2_RR = np.tensordot(R2, R2 * freqs[np.newaxis, :, np.newaxis], axes=[1, 1])
        G2_RR = sum((G2_RR.transpose(*p) for p in [[0, 2, 1, 3], [0, 2, 3, 1]]))
        G2_R3 = R3 * freqs[np.newaxis, np.newaxis, :, np.newaxis]
        G2_R3 = G2_R3 + G2_R3.transpose(0, 1, 3, 2)
        G2_RG = np.tensordot(R2, test_G[1], axes=[1, 1])
        G2_RG = sum((G2_RG.transpose(*p) for p in [[0, 2, 1, 3], [2, 0, 1, 3], [0, 2, 3, 1], [2, 0, 3, 1]]))
        G2_QG = np.tensordot(Q2, test_G[1], axes=[2, 0])
        G2 = G2_QR + G2_RR + G2_RG + G2_QG + G2_R3 + test_G[2]
        ggg = [G0, G1, G2]
        ggg_test = runner1.hamiltonian.reexpress_G([Q1, Q2, Q3], [R1, R2, R3])
        (Q, R), _ = runner1.hamiltonian.get_potential_optimized_coordinates()
        rrr, _ = VPTRunner.construct(TestManager.test_data(file_name), 2, internals=test_internals, kinetic_terms=[G0, test_G[1], test_G[2]], potential_terms=[G0, test_V1[2], test_V1[3]], include_pseudopotential=False, include_coriolis_coupling=True, logger=False)
        rrr.print_tables(print_intensities=False)
        '\n        :: State    <0|dH(2)|0>  <0|dH(1)|1> \n              0 0 0    100.74075   -156.70145\n              0 0 1    296.33507   -545.08549\n              0 1 0    301.26122   -538.53511\n              1 0 0    108.13556   -213.67502\n              0 0 2    587.47448  -1127.09543\n              0 2 0    597.90621  -1104.59990\n              2 0 0    109.32398   -292.47319\n              0 1 1    687.69956  -1284.13087\n              1 0 1    312.98434   -634.36364\n              1 1 0    326.60393   -633.50201\n        '
        '\n              0 0 0    105.49844   -157.99275\n              0 0 1    335.26400   -580.55272\n              0 1 0    279.87458   -513.67920\n              1 0 0    105.42497   -207.50600\n              0 0 2    703.29509  -1239.45907\n              0 2 0    554.24618  -1057.46833\n              2 0 0    108.34454   -288.04967\n              0 1 1    678.58740  -1271.55370\n              1 0 1    349.24337   -667.17858\n              1 1 0    259.61523   -563.04545\n              '
        '\n            0 0 0     10.85263   -144.38370\n            0 0 1    103.73295   -428.81417\n            0 1 0    102.30966   -416.98116\n            1 0 0    -35.81457   -148.46241\n            0 0 2    249.91709   -862.39340\n            0 2 0    252.94485   -836.86328\n            2 0 0    -83.05123   -180.01127\n            0 1 1    313.14314   -985.39758\n            1 0 1     18.26562   -419.47037\n            1 1 0      2.24708   -387.69952\n            '
        rrr, _ = VPTRunner.construct(TestManager.test_data(file_name), 2, internals=test_internals, kinetic_terms=ggg_test, potential_terms=[G0, 0, 0], include_pseudopotential=False, include_coriolis_coupling=True, logger=False, zero_element_warning=False)
        rrr.print_tables(print_intensities=False)
        '\n          0 0 0    176.43652   -163.08293\n          0 0 1    289.53377   -468.96990\n          0 1 0    298.51621   -466.47580\n          1 0 0    170.71509   -206.94025\n          '
        raise Exception(...)
        '\n              0 0 0    176.43652   -163.08293\n              0 0 1    289.53377   -468.96990\n              0 1 0    298.51621   -466.47580\n              1 0 0    170.71509   -206.94025\n              0 0 2    455.91832   -926.22498\n              0 2 0    479.53208   -916.91148\n              2 0 0    153.91367   -267.74859\n              0 1 1    529.20837  -1056.32539\n              1 0 1    285.83732   -537.90233\n              1 1 0    301.48690   -539.07068\n        '
        '\n          0 0 0     88.73554   -172.09145\n          0 0 1    239.96002   -516.11032\n          0 1 0    191.24982   -455.91603\n          1 0 0    108.69731   -241.63993\n          0 0 2    461.67666  -1028.70222\n          0 2 0    353.91552   -887.99927\n          2 0 0    138.27476   -348.84148\n          0 1 1    455.22336  -1079.05125\n          1 0 1    308.63734   -657.43414\n          1 1 0    183.18822   -517.48003\n          '
        raise Exception(...)
        raise Exception(...)
        VPTRunner.run_simple(TestManager.test_data(file_name), 2, internals=[[1, -1, -1, -1], [2, 1, -1, -1], [0, 1, 2, -1]], include_pseudopotential=False, potential_derivatives=[np.zeros_like(m) if i != 1 else m for i, m in enumerate(mol.potential_derivatives)], zero_element_warning=False, calculate_intensities=False)
        VPTRunner.run_simple(TestManager.test_data(file_name), 2, include_pseudopotential=False, potential_derivatives=[np.zeros_like(m) if i != 1 else m for i, m in enumerate(mol.potential_derivatives)], internals=[[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]], zero_element_warning=False, calculate_intensities=False)
