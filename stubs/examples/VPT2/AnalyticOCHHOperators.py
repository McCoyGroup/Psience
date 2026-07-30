"""Extracted from VPT2Tests.test_AnalyticOCHHOperators via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VPT2Tests.test_AnalyticOCHHOperators"""

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
    def test_AnalyticOCHHOperators(self):
        """
        Run OCHH, add in the single degeneracy by hand
        :return:
        """
        file_name = 'OCHH_freq.fchk'
        AnalyticVPTRunner.run_simple(TestManager.test_data(file_name), [[0, [[0, 0, 0, 0, 0, 1], [0, 1, 0, 1, 0, 0], [0, 0, 0, 1, 1, 0], [0, 0, 0, 0, 1, 0]]], [[0, 0, 0, 0, 1, 0], [[0, 0, 0, 0, 1, 1], [0, 1, 0, 1, 1, 0]]]], calculate_intensities=True, degeneracy_specs=[[[0, 0, 0, 0, 0, 1], [0, 1, 0, 1, 0, 0]]], operator_terms={'q1': [0, [0, 0, 0, 0, 0, 1]], 'q2': [0, [0, 0, 0, 0, 1, 0]]})
        '\n::> Operator Corrections:\n  > =========================================================== () ===========================================================\n  >            |         Operator          |          Order 0          |          Order 1          |          Order 2         \n  >   States   |     q1      |     q2      |     q1      |     q2      |     q1      |     q2      |     q1      |     q2     \n  > --------------------------------------------------------------------------------------------------------------------------\n  >       1(1) | -0.51944545 |  0.00000107 | -0.70710678 | -0.00000000 | -0.00000000 | -0.00000000 |  0.18766133 |  0.00000107\n  >   3(1)5(1) | -0.50347989 | -0.00005836 | -0.00000000 | -0.00000000 | -0.50347989 | -0.00005836 | -0.00000000 | -0.00000000\n  >   2(1)3(1) | -0.00000000 | -0.00169639 | -0.00000000 | -0.00000000 | -0.00000000 | -0.00169639 | -0.00000000 | -0.00000000\n  >       2(1) |  0.00004036 | -0.68109787 | -0.00000000 | -0.70710678 | -0.00000000 | -0.00000000 |  0.00004036 |  0.02600891\n  > ============================================================ 2(1) ============================================================\n  >                |         Operator          |          Order 0          |          Order 1          |          Order 2         \n  >     States     |     q1      |     q2      |     q1      |     q2      |     q1      |     q2      |     q1      |     q2     \n  > ------------------------------------------------------------------------------------------------------------------------------\n  >       1(1)2(1) | -0.52457444 |  0.00000107 | -0.70710678 | -0.00000000 | -0.00000000 | -0.00000000 |  0.18253234 |  0.00000107\n  >   2(1)3(1)5(1) | -0.50347989 | -0.00005836 | -0.00000000 | -0.00000000 | -0.50347989 | -0.00005836 | -0.00000000 | -0.00000000\n<::\n    '
