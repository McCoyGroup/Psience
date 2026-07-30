"""Extracted from VPT2Tests.test_HOHCorrectedDegeneracies via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VPT2Tests.test_HOHCorrectedDegeneracies"""

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
    def test_HOHCorrectedDegeneracies(self):
        VPTRunner.run_simple(TestManager.test_data('HOH_freq.fchk'), 2, zero_order_energy_corrections=[[(0, 1, 0), (4681.56364 + 3800) * UnitsData.convert('Wavenumbers', 'Hartrees')], [(0, 2, 0), (4681.56364 + 7800) * UnitsData.convert('Wavenumbers', 'Hartrees')], [(0, 0, 2), (4681.56364 + 7801) * UnitsData.convert('Wavenumbers', 'Hartrees')], [(0, 3, 1), (4681.56364 + 3801) * UnitsData.convert('Wavenumbers', 'Hartrees')]], degeneracy_specs={'wfc_threshold': 0.3, 'extra_groups': [[[0, 2, 0], [0, 1, 1]], [[0, 0, 2], [0, 3, 1]]]})
