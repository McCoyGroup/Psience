"""Extracted from VPT2Tests.test_HOONOFasterDegenSubspace via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VPT2Tests.test_HOONOFasterDegenSubspace"""

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
    def test_HOONOFasterDegenSubspace(self):
        file_name = 'HOONO_freq.fchk'
        VPTRunner.run_simple(TestManager.test_data(file_name), [[0, 0, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 1, 0, 0], [1, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 1, 0, 0], [0, 0, 1, 0, 0, 0, 0, 1, 0]], logger=True, degeneracy_specs='auto', extended_space_target_property='frequencies')
        '\nState                   Frequency    Intensity       Frequency    Intensity\n  1 0 0 0 0 0 1 0 0    1789.34288      0.00000      1669.11650      0.00850\n  1 0 0 0 0 0 0 1 0    1924.23563      0.00000      1825.06905      0.30108\n  0 0 1 0 0 0 1 0 0    1957.70875      0.00000      1886.51516      0.05466\n  0 0 1 0 0 0 0 1 0    2092.60150      0.00000      2045.13170      0.02941\n  1 0 0 2 0 0 0 0 0    1787.50010      0.00000      2413.94904     80.96353\n  0 0 1 2 0 0 0 0 0    1955.86597      0.00000      1971.94764      2.37340\n  1 0 0 1 1 0 0 0 0    1908.01157      0.00000      2131.62528      0.72881\n  0 0 1 1 1 0 0 0 0    2076.37744      0.00000      2080.76165      2.72838\n  '
        VPTRunner.run_simple(TestManager.test_data(file_name), [[0, 0, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 1, 0, 0], [1, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 1, 0, 0], [0, 0, 1, 0, 0, 0, 0, 1, 0]], logger=True, degeneracy_specs='auto')
        '\nState                   Frequency    Intensity       Frequency    Intensity\n  1 0 0 0 0 0 1 0 0    1789.34288      0.00000      1669.11650      0.04924\n  1 0 0 0 0 0 0 1 0    1924.23563      0.00000      1825.06905      0.07528\n  0 0 1 0 0 0 1 0 0    1957.70875      0.00000      1886.51516      0.00766\n  0 0 1 0 0 0 0 1 0    2092.60150      0.00000      2045.13170      0.01768\n  1 0 0 2 0 0 0 0 0    1787.50010      0.00000      2413.94904      0.08504\n  0 0 1 2 0 0 0 0 0    1955.86597      0.00000      1971.94764      0.00159\n  1 0 0 1 1 0 0 0 0    1908.01157      0.00000      2131.62528      0.00035\n  0 0 1 1 1 0 0 0 0    2076.37744      0.00000      2080.76165      0.00332\n        '
