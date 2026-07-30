"""Extracted from VPT2Tests.test_PartialRebuild via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VPT2Tests.test_PartialRebuild"""

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
    def test_PartialRebuild(self):
        state = VPTStateMaker(7)
        corrs = AnalyticVPTRunner.run_simple(mol, [[[state()], [state(1), state(2), state(3), state([1, 2]), state([2, 2]), state([3, 2]), state(1, 2), state(1, 3), state(2, 3)]]], full_surface_mode_selection=[108 - 108, 108 - 107, 108 - 106, 108 - 105, 108 - 21, 108 - 20, 108 - 19, 108 - 1], mode_selection=[108 - 108, 108 - 107, 108 - 106, 108 - 105, 108 - 21, 108 - 20, 108 - 19], expressions_file=os.path.expanduser('exprs.hdf5'))
