"""Extracted from VPT2Tests.test_TrimerAnalytic via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VPT2Tests.test_TrimerAnalytic"""

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
    def test_TrimerAnalytic(self):
        file_name = 'water_trimer_freq.fchk'
        runner, states = AnalyticVPTRunner.construct(TestManager.test_data(file_name), 2, logger=True)
        og, _ = runner.construct_classic_runner(TestManager.test_data(file_name), states, mode_selection=np.arange(len(states[0])), logger=False)
        with BlockProfiler(print_options={'show_all': True}):
            spec = runner.get_spectrum(states, verbose=False)
        print(np.array(spec).T)
