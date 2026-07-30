"""Extracted from VPT2Tests.test_HOHVPTAnneManip via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VPT2Tests.test_HOHVPTAnneManip"""

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

    @validationTest
    def test_HOHVPTAnneManip(self):
        runner, _ = VPTRunner.helpers.run_anne_job(TestManager.test_data('vpt2_helpers_api/hod/r'), return_runner=True, order=2, expansion_order=2)
        H = runner.hamiltonian
        V = H.V_terms
        freqs = np.diag(V[0])
        G = H.G_terms
        U = H.pseudopotential_term
        D = H.expansion_options['dipole_terms']
        frequency_shift = np.array([-1, 0, 0]) * UnitsData.convert('Wavenumbers', 'Hartrees')
        new_freqs = freqs + frequency_shift
        scaling_factor = np.sqrt(new_freqs) / np.sqrt(freqs)
        G[2]
        v_expansion = [V[0] * (scaling_factor[:, np.newaxis] * scaling_factor[np.newaxis, :]), V[1] * (scaling_factor[:, np.newaxis, np.newaxis] * scaling_factor[np.newaxis, :, np.newaxis] * scaling_factor[np.newaxis, np.newaxis, :]), V[2] * (scaling_factor[:, np.newaxis, np.newaxis, np.newaxis] * scaling_factor[np.newaxis, :, np.newaxis, np.newaxis] * scaling_factor[np.newaxis, np.newaxis, :, np.newaxis] * scaling_factor[np.newaxis, np.newaxis, np.newaxis, :])]
        g_expansion = [G[0] * (scaling_factor[:, np.newaxis] * scaling_factor[np.newaxis, :]), G[1] * (scaling_factor[:, np.newaxis, np.newaxis] / scaling_factor[np.newaxis, :, np.newaxis] / scaling_factor[np.newaxis, np.newaxis, :]), G[2] * (scaling_factor[:, np.newaxis, np.newaxis, np.newaxis] * scaling_factor[np.newaxis, :, np.newaxis, np.newaxis] / scaling_factor[np.newaxis, np.newaxis, :, np.newaxis] / scaling_factor[np.newaxis, np.newaxis, np.newaxis, :])]
        u_expansion = [U[0]]
        d_expansion = [[D[a][0], D[a][1] * scaling_factor[:], D[a][2] * (scaling_factor[:, np.newaxis] * scaling_factor[np.newaxis, :]), D[a][3] * (scaling_factor[:, np.newaxis, np.newaxis] * scaling_factor[np.newaxis, :, np.newaxis] * scaling_factor[np.newaxis, np.newaxis, :])] for a in range(3)]
        new_runner, _ = VPTRunner.construct(runner.system.mol, runner.states.state_list, potential_terms=v_expansion, kinetic_terms=g_expansion, pseudopotential_terms=u_expansion, dipole_terms=d_expansion)
        runner.print_tables()
        new_runner.print_tables()
