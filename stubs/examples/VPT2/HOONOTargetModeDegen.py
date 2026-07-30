"""Extracted from VPT2Tests.test_HOONOTargetModeDegen via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VPT2Tests.test_HOONOTargetModeDegen"""

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
    def test_HOONOTargetModeDegen(self):
        file_name = 'HOONO_freq.fchk'
        VPTRunner.run_simple(TestManager.test_data(file_name), 2, logger=True, degeneracy_specs={'wfc_threshold': 0.3, 'max_mode_differences': [0 if i in {1, 2, 3, 4} else -1 for i in range(9)]}, target_property='frequencies', calculate_intensities=False)
