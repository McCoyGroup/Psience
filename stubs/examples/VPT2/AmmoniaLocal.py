"""Extracted from VPT2Tests.test_AmmoniaLocal via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VPT2Tests.test_AmmoniaLocal"""

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
    def test_AmmoniaLocal(self):
        file_name = 'nh3.fchk'
        nh3 = Molecule.from_file(TestManager.test_data(file_name))
        XR, RX = nput.internal_coordinate_tensors(nh3.coords, [(0, 1), (0, 2), (0, 3), (1, 0, 2), (1, 0, 3), (2, 0, 3)], return_inverse=True, order=1)
        VPTRunner.run_simple(TestManager.test_data(file_name), 1, local_modes={'matrix': XR[1], 'inverse': RX[0], 'sort_freqs': True}, mixed_derivative_handling_mode='analytical', local_mode_couplings=True, degeneracy_specs={'polyads': [[[0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 1, 0]], [[0, 0, 0, 0, 0, 1], [0, 0, 0, 1, 0, 0]], [[0, 0, 1, 0, 0, 0], [0, 1, 0, 0, 0, 0]], [[0, 0, 1, 0, 0, 0], [1, 0, 0, 0, 0, 0]]]}, calculate_intensities=True)
        return
