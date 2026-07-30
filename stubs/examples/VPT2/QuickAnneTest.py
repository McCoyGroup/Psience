"""Extracted from VPT2Tests.test_QuickAnneTest via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VPT2Tests.test_QuickAnneTest"""

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
    def test_QuickAnneTest(self):
        AnalyticVPTRunner.run_simple(os.path.expanduser('~/Documents/Postdoc/Projects/Anne_Misc/formate/Dformate_OH+ezWeD2z_anh_tz.fchk'), states=2, order=2, operator_chunk_size=100000, degeneracy_specs={'wfc_threshold': 0.3}, logger=os.path.expanduser('~/Documents/Postdoc/Projects/Anne_Misc/formate/formate.out'), internals=VPTRunner.helpers.parse_zmatrix(os.path.expanduser('~/Documents/Postdoc/Projects/Anne_Misc/formate/z_mat_WeD2z.dat')))
