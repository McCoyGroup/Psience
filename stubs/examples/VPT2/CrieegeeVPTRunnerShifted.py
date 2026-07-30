"""Extracted from VPT2Tests.test_CrieegeeVPTRunnerShifted via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VPT2Tests.test_CrieegeeVPTRunnerShifted"""

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
    def test_CrieegeeVPTRunnerShifted(self):
        freqs = VPTSystem('criegee_eq_anh.fchk').mol.normal_modes.modes.freqs
        freqs = freqs.copy()
        freqs[1] += 10 / UnitsData.convert('Hartrees', 'Wavenumbers')
        VPTRunner.run_simple(2, logger=True, degeneracy_specs='auto', corrected_fundamental_frequencies=freqs)
