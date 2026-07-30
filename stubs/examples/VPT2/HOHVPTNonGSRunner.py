"""Extracted from VPT2Tests.test_HOHVPTNonGSRunner via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VPT2Tests.test_HOHVPTNonGSRunner"""

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
    def test_HOHVPTNonGSRunner(self):
        file_name = 'HOH_freq.fchk'
        mol = Molecule.from_file(TestManager.test_data(file_name), internals=[[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]])
        ics = mol.internal_coordinates
        derivs = mol.coords.jacobian(ics.system, [1, 2], all_numerical=True, converter_options={'reembed': True})
        derivs = [x.reshape((9,) * (i + 1) + (3, 3)) for i, x in enumerate(derivs)]
        VPTRunner.run_simple(TestManager.test_data(file_name), 2, initial_states=1, operators={'OH1': [ics[1, 0], derivs[0][:, 1, 0], derivs[1][:, :, 1, 0]], 'OH2': [ics[2, 0], derivs[0][:, 2, 0], derivs[1][:, :, 2, 0]], 'HOH': [ics[2, 1], derivs[0][:, 2, 1], derivs[1][:, :, 2, 1]]}, logger=True)
