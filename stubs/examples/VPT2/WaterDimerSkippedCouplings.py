"""Extracted from VPT2Tests.test_WaterDimerSkippedCouplings via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VPT2Tests.test_WaterDimerSkippedCouplings"""

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
    def test_WaterDimerSkippedCouplings(self):
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
        internals = [[LHF, X, X, X], [LO, LHF, X, X], [SH, LO, LHF, X], [RH2, SH, LO, LHF], [RO, LO, RH2, LHF], [RH1, RO, RH2, LHF]]
        VPTRunner.run_simple(TestManager.test_data('water_dimer_freq.fchk'), 1, degeneracy_specs='auto', mode_selection=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
