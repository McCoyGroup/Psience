"""Extracted from VPT2Tests.test_PyreneAnalytic via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VPT2Tests.test_PyreneAnalytic"""

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
    def test_PyreneAnalytic(self):
        file_name = os.path.expanduser('~/Documents/Postdoc/Projects/VPT/anne_tests/pyrene_try3.fchk')
        state = VPTStateMaker(72)
        degs = [[state(8), state(9)], [state(15), state(69)], [state(16), state(39)], [state(19), state(20)], [state(19), state(15)], [state(20), state(16)], [state(21), state(70)], [state(22), state(18)], [state(23), state(24)], [state(24), state(55)], [state(27), state(28)], [state(35), state(66)], [state(36), state(37)], [state(37), state(15)], [state(38), state(69)], [state(39), state(40)], [state(41), state(71)], [state(42), state(43)], [state(67), state(68)], [state(67), state(37)], [state(68), state(69)], [state(71), state(18)], [state(71), state(22)]]
        par = Parallelizer.lookup(('multiprocessing', {'initialization_timeout': 25}))
        with par:
            with Timer():
                corrs = AnalyticVPTRunner.run_simple(file_name, 1, degeneracy_specs=degs, logger=True, parallelizer=par)
