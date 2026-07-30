"""Extracted from VPT2Tests.test_AnalyticOCHHMultiple via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VPT2Tests.test_AnalyticOCHHMultiple"""

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
    def test_AnalyticOCHHMultiple(self):
        file_name = 'OCHH_freq.fchk'
        file_name = 'OCHH_freq.fchk'
        runner = AnalyticVPTRunner.run_simple(TestManager.test_data(file_name), [[0, [[0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 1, 0], [0, 0, 0, 2, 0, 0]]], [[[0, 0, 0, 0, 0, 1]], [[0, 0, 0, 0, 1, 0], [0, 0, 0, 2, 0, 0]]]], mixed_derivative_handling_mode='analytical', order=2, expansion_order=2, degeneracy_specs={'polyads': [[[0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 1, 0]], [[0, 0, 0, 0, 0, 1], [0, 1, 0, 1, 0, 0]]]})
        raise Exception(...)
        classic, _ = runner.construct_classic_runner(states, zero_element_warning=False)
        runner.run_VPT(states)
        raise Exception(...)
