"""Extracted from VPT2Tests.test_SelAnharmAnalytic via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VPT2Tests.test_SelAnharmAnalytic"""

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
    def test_SelAnharmAnalytic(self):
        file_name = 'freq_anion.fchk'
        state = VPTStateMaker(108)
        corrs = AnalyticVPTRunner.run_simple(file_name, [state(), state(21), state(20), state(19)], full_surface_mode_selection=[108 - 21, 108 - 20, 108 - 19, 108 - 1], logger=True, expressions_file=os.path.expanduser('exprs.hdf5'))
        print(corrs)
        raise Exception(...)
