"""Extracted from VPT2Tests.test_BlockLabels via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VPT2Tests.test_BlockLabels"""

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
    def test_BlockLabels(self):
        VPTRunner.run_simple(TestManager.test_data('i_hoh_opt.fchk'), VPTStateSpace.get_state_list_from_quanta(4, 6) + [[0, 1, 2, 2, 0, 0]], initial_states=[[0, 0, 0, 0, 0, 0], [0, 0, 0, 2, 0, 0], [0, 1, 0, 2, 0, 0], [0, 0, 0, 0, 1, 0]], degeneracy_specs={'wfc_threshold': 0.3}, logger=True, plot_spectrum=False)
