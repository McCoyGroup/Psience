"""Extracted from VPT2Tests.test_OCHHFasterDegenSubspace via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VPT2Tests.test_OCHHFasterDegenSubspace"""

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
    def test_OCHHFasterDegenSubspace(self):
        file_name = 'OCHH_freq.fchk'
        VPTRunner.run_simple(TestManager.test_data(file_name), [[0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 1]], logger=True, degeneracy_specs='auto', extended_space_target_property='frequencies')
        '\n::> IR Data\n  > Initial State: 0 0 0 0 0 0 \n                         Harmonic                  Anharmonic\nState             Frequency    Intensity       Frequency    Intensity\n  0 0 1 0 0 1    4588.74227      0.00000      4446.79734      0.00134\n  0 1 2 0 0 0    4306.24556      0.00000      4167.25991      0.29964\n  0 1 1 1 0 0    4506.28742      0.00000      4345.62126      0.48923\n  '
        VPTRunner.run_simple(TestManager.test_data(file_name), [[0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 1]], logger=True, degeneracy_specs='auto')
        '\n::> IR Data\n  > Initial State: 0 0 0 0 0 0 \n                         Harmonic                  Anharmonic\nState             Frequency    Intensity       Frequency    Intensity\n  0 0 1 0 0 1    4588.74227      0.00000      4446.79734      0.00000\n  0 1 2 0 0 0    4306.24556      0.00000      4167.25991      0.26680\n  0 1 1 1 0 0    4506.28742      0.00000      4345.62126      0.47084\n  '
