"""Extracted from VPT2Tests.test_HOHCorrectedPostfilters via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VPT2Tests.test_HOHCorrectedPostfilters"""

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
    def test_HOHCorrectedPostfilters(self):
        VPTAnalyzer.run_VPT(TestManager.test_data('HOH_freq.fchk'), 2, zero_order_energy_corrections=[[(0, 1, 0), (4681.56364 + 3800) * UnitsData.convert('Wavenumbers', 'Hartrees')], [(0, 2, 0), (4681.56364 + 7800) * UnitsData.convert('Wavenumbers', 'Hartrees')], [(0, 0, 2), (4681.56364 + 7801) * UnitsData.convert('Wavenumbers', 'Hartrees')]], degeneracy_specs='auto', basis_filters=[{'max_quanta': [2, -1, -1]}, {'excluded_transitions': [[1, 0, 0], [0, 1, 0], [0, 0, 1]]}]).print_output_tables()
