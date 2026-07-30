"""Extracted from VPT2Tests.test_AnalyticHOONO via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VPT2Tests.test_AnalyticHOONO"""

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
    def test_AnalyticHOONO(self):
        file_name = 'HOONO_freq.fchk'
        state = VPTStateMaker(9)
        from McUtils.Parallelizers import Parallelizer
        from McUtils.Profilers import Timer, BlockProfiler
        AnalyticVPTRunner.run_simple(TestManager.test_data(file_name), 2, degeneracy_specs=[[state(1), state(8, 7)], [state(3), state([6, 2])], [state(6), state([9, 2])]])
        return
        runner, states = AnalyticVPTRunner.construct(TestManager.test_data(file_name), 2, logger=True)
        og, _ = runner.construct_classic_runner(TestManager.test_data(file_name), states, mode_selection=np.arange(len(states[0])), zero_element_warning=False, logger=False)
