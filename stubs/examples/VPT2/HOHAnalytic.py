"""Extracted from VPT2Tests.test_HOHAnalytic via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VPT2Tests.test_HOHAnalytic"""

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
    def test_HOHAnalytic(self):
        file_name = 'HOH_freq.fchk'
        runner, states = AnalyticVPTRunner.construct(TestManager.test_data(file_name), 2, degeneracy_specs={'polyads': [[[0, 0, 1], [0, 1, 0]]]})
        classic, _ = runner.construct_classic_runner(states)
        runner.run_VPT(states, operator_expansions={'x1': [0, np.eye(9)[:, 0], np.zeros((9, 9))], 'y1': [0, np.eye(9)[:, 1], np.zeros((9, 9))], 'z1': [0, np.eye(9)[:, 2], np.zeros((9, 9))]}, calculate_intensities=True)
        runner, states = AnalyticVPTRunner.construct(TestManager.test_data(file_name), 2, internals=[[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]])
        classic2, _ = runner.construct_classic_runner(states, zero_element_warning=False)
        runner.run_VPT(states)
