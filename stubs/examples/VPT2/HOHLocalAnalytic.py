"""Extracted from VPT2Tests.test_HOHLocalAnalytic via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VPT2Tests.test_HOHLocalAnalytic"""

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
    def test_HOHLocalAnalytic(self):
        file_name = 'HOH_freq.fchk'
        OHH = Molecule.from_file(TestManager.test_data(file_name), internals=[[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]])
        XR = OHH.get_internals_by_cartesians(1, strip_embedding=True)[0]
        RX = OHH.get_cartesians_by_internals(1, strip_embedding=True)[0]
        AnalyticVPTRunner.run_simple(TestManager.test_data(file_name), 1, local_modes={'matrix': XR, 'inverse': RX, 'sort_freqs': True}, local_mode_coupling_order=2, mixed_derivative_handling_mode='analytical', degeneracy_specs={'polyads': [[[0, 0, 1], [0, 1, 0]], [[0, 0, 1], [2, 0, 0]]]}, calculate_intensities=True)
        '\n2nd order, Primary\n::        |     Harmonic      |    Anharmonic     |    Deperturbed   \n:: States |  Freq.   |  Int.  |  Freq.   |  Int.  |  Freq.   |  Int. \n:: ------------------------------------------------------------------\n::     () |    0.000 |  0.000 |    0.000 |  0.000 |    0.000 |  0.000\n::   1(1) | 3873.508 | 35.315 | 3636.877 |  4.658 | 3688.685 | 34.090\n::   2(1) | 3873.508 | 35.315 | 3752.804 | 64.983 | 3688.980 | 33.922\n::   3(1) | 1641.379 | 65.147 | 1597.321 | 69.553 | 1597.321 | 69.553\n::   3(2) | 3282.757 |  0.000 | 3171.933 |  0.413 | 3183.950 |  0.922\n\n2nd order, Degenerate\n::        |     Harmonic      |    Anharmonic     |    Deperturbed   \n:: States |  Freq.   |  Int.  |  Freq.   |  Int.  |  Freq.   |  Int. \n:: ------------------------------------------------------------------\n::     () |    0.000 |  0.000 |    0.000 |  0.000 |    0.000 |  0.000\n::   1(1) | 3873.508 | 35.315 | 3636.877 |  4.658 | 3688.685 | 34.090\n::   2(1) | 3873.508 | 35.315 | 3752.804 | 64.983 | 3688.980 | 33.922\n::   3(1) | 1641.379 | 65.147 | 1597.321 | 69.553 | 1597.321 | 69.553\n::   3(2) | 3282.757 |  0.000 | 3171.933 |  0.413 | 3183.950 |  0.922\n\n1st order, Degenerate\n::        |     Harmonic      |    Anharmonic     |    Deperturbed   \n:: States |  Freq.   |  Int.  |  Freq.   |  Int.  |  Freq.   |  Int. \n:: ------------------------------------------------------------------\n::     () |    0.000 |  0.000 |    0.000 |  0.000 |    0.000 |  0.000\n::   1(1) | 3873.508 | 35.315 | 3635.945 |  4.898 | 3688.685 | 34.297\n::   2(1) | 3873.508 | 35.315 | 3752.804 | 65.133 | 3688.980 | 34.113\n::   3(1) | 1641.379 | 65.147 | 1597.321 | 67.286 | 1597.321 | 67.286\n::   3(2) | 3282.757 |  0.000 | 3172.866 |  0.418 | 3183.950 |  0.922\n\n1st order, Primary\n::     () |    0.000 |  0.000 |    0.000 |  0.000 |    0.000 |  0.000\n::   1(1) | 3873.508 | 35.315 | 3629.147 |  4.667 | 3685.563 | 33.922\n::   2(1) | 3873.508 | 35.315 | 3752.794 | 64.662 | 3685.858 | 33.751\n::   3(1) | 1641.379 | 65.147 | 1578.730 | 69.539 | 1578.730 | 69.539\n::   3(2) | 3282.757 |  0.000 | 3136.249 |  0.505 | 3146.768 |  0.996\n\n1st order, Primary, Internals\n::        |     Harmonic      |    Anharmonic     |    Deperturbed   \n:: States |  Freq.   |  Int.  |  Freq.   |  Int.  |  Freq.   |  Int. \n:: ------------------------------------------------------------------\n::     () |    0.000 |  0.000 |    0.000 |  0.000 |    0.000 |  0.000\n::   1(1) | 3873.508 | 35.315 | 3633.060 |  3.796 | 3692.440 | 33.437\n::   2(1) | 3873.508 | 35.315 | 3753.652 | 63.777 | 3692.745 | 32.936\n::   3(1) | 1641.379 | 65.147 | 1570.743 | 69.153 | 1570.743 | 69.153\n::   3(2) | 3282.757 |  0.000 | 3108.288 |  0.879 | 3109.816 |  1.069\n>>--------------------------------------------------<<\n'
