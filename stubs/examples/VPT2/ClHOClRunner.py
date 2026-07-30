"""Extracted from VPT2Tests.test_ClHOClRunner via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VPT2Tests.test_ClHOClRunner"""

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
    def test_ClHOClRunner(self):
        file_name = 'cl_hocl.fchk'
        state = VPTStateMaker(6)
        COM = -3
        A = -2
        C = -1
        _ = 1000
        O = 0
        H = 1
        Cl = 2
        X = 3
        VPTRunner.run_simple(TestManager.test_data(file_name), [state(), state([1, 1]), state([1, 2]), state([1, 3]), state([1, 2], [5, 1]), state([1, 1], [2, 2])], logger=True, handle_strong_couplings=True, strong_coupling_test_modes=list(range(3, 6)))
        '\n        > [[0 0 0 0 0 2]\n        >  [1 0 0 0 0 2]\n        >  [0 1 0 0 0 2]\n        >  [0 0 0 0 2 1]\n        >  [0 0 0 0 4 0]\n        >  [0 2 0 0 0 2]\n        >  [0 1 0 0 2 1]]\n                                 Harmonic                  Anharmonic\n        State             Frequency    Intensity       Frequency    Intensity\n          0 0 0 0 0 1    2709.16096   2782.25433      2285.38768   2611.96281\n          0 0 0 0 0 2    5418.32192      0.00000      4140.45935     19.13726\n          0 0 0 0 0 3    8127.48288      0.00000      5353.97008      0.16004\n          0 1 0 0 0 2    5699.88024      0.00000      4605.01290      3.93389\n          0 0 0 0 2 1    5592.48466      0.00000      5023.09956      7.99053\n        '
        VPTRunner.run_simple(TestManager.test_data(file_name), [state(), state([1, 1]), state([1, 2]), state([1, 3]), state([1, 2], [5, 1]), state([1, 1], [2, 2])], logger=True, handle_strong_couplings=True)
        '        \n        > [[0 0 0 0 0 2]\n        >  [0 0 0 0 2 1]\n        >  [0 0 0 0 4 0]]\n                                 Harmonic                  Anharmonic\n        State             Frequency    Intensity       Frequency    Intensity\n          0 0 0 0 0 1    2709.16096   2782.25433      2285.38768   2611.96281\n          0 0 0 0 0 2    5418.32192      0.00000      4096.10015     21.49976\n          0 0 0 0 0 3    8127.48288      0.00000      5353.97008      0.16004\n          0 1 0 0 0 2    5699.88024      0.00000      4374.63719      2.66632\n          0 0 0 0 2 1    5592.48466      0.00000      4964.16010      6.86648\n        '
        VPTRunner.run_simple(TestManager.test_data(file_name), [state(), state([1, 1]), state([1, 2]), state([1, 3]), state([1, 2], [5, 1]), state([1, 1], [2, 2])], logger=True, handle_strong_couplings=True, strong_coupling_test_modes=list(range(3, 6)), internals=[[Cl, _, _, _], [O, Cl, _, _], [X, O, Cl, _], [H, O, Cl, X]])
        '        \n        > [[0 0 0 0 0 2]\n        >  [1 0 0 0 0 2]\n        >  [0 1 0 0 0 2]\n        >  [0 0 0 0 2 1]\n        >  [0 0 0 0 4 0]\n        >  [0 2 0 0 0 2]\n        >  [0 1 0 0 2 1]]\n                                 Harmonic                  Anharmonic\n        State             Frequency    Intensity       Frequency    Intensity\n          0 0 0 0 0 1    2709.16096   2782.25433      2305.91708   2613.72238\n          0 0 0 0 0 2    5418.32192      0.00000      4174.53466     17.82679\n          0 0 0 0 0 3    8127.48288      0.00000      5353.97246      0.16004\n          0 1 0 0 0 2    5699.88024      0.00000      4645.21985      6.56429\n          0 0 0 0 2 1    5592.48466      0.00000      5114.27886      9.19226\n        '
        VPTRunner.run_simple(TestManager.test_data(file_name), [state(), state([1, 1]), state([1, 2]), state([1, 3]), state([1, 2], [5, 1]), state([1, 1], [2, 2])], logger=True, handle_strong_couplings=True, internals=[[Cl, _, _, _], [O, Cl, _, _], [X, O, Cl, _], [H, O, Cl, X]])
        '\n        > [[0 0 0 0 0 2]\n        >  [0 0 0 0 2 1]\n        >  [0 0 0 0 4 0]]\n                                 Harmonic                  Anharmonic\n        State             Frequency    Intensity       Frequency    Intensity\n          0 0 0 0 0 1    2709.16096   2782.25433      2305.91708   2613.72238\n          0 0 0 0 0 2    5418.32192      0.00000      4130.15930     22.20869\n          0 0 0 0 0 3    8127.48288      0.00000      5353.97246      0.16004\n          0 1 0 0 0 2    5699.88024      0.00000      4374.63638      2.66634\n          0 0 0 0 2 1    5592.48466      0.00000      5053.08138      7.79496\n        '
