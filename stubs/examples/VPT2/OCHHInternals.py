"""Extracted from VPT2Tests.test_OCHHInternals via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VPT2Tests.test_OCHHInternals"""

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
    def test_OCHHInternals(self):
        tag = 'OCHH Internals'
        file_name = 'OCHH_freq.fchk'
        zmatrix = [[0, -1, -1, -1], [1, 0, -1, -1], [2, 1, 0, -1], [3, 1, 0, 2]]

        def conv(r, t, f, **kwargs):
            return np.array([r ** 2 + 1, t, f])

        def inv(r2, t, f, **kwargs):
            return np.array([np.sqrt(r2 - 1), t, f])
        internals = {'zmatrix': zmatrix, 'conversion': conv, 'inverse': inv}
        print('Cart:')
        print('ZMat:')
        print('Custom:')
        VPTAnalyzer.run_VPT(TestManager.test_data(file_name), 2, internals=internals).print_output_tables()
