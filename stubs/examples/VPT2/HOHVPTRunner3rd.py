"""Extracted from VPT2Tests.test_HOHVPTRunner3rd via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VPT2Tests.test_HOHVPTRunner3rd"""

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
    def test_HOHVPTRunner3rd(self):
        """
        test that runner works for 3rd order PT, too

        :return:
        :rtype:
        """
        file_name = 'HOH_freq.fchk'
        handling_mode = 'unhandled'
        logger = Logger()
        with logger.block(tag='Internals 2nd Order triad'):
            VPTRunner.run_simple(TestManager.test_data(file_name), 3, logger=logger, order=2, internals=[[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]], expansion_order=2, degeneracy_specs=[[[0, 0, 1], [2, 0, 0]], [[0, 1, 0], [2, 0, 0]]])
        with logger.block(tag='Internals 3rd Order triad'):
            VPTRunner.run_simple(TestManager.test_data(file_name), 3, logger=logger, order=3, internals=[[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]], expansion_order=2, degeneracy_specs=[[[0, 0, 1], [2, 0, 0]], [[0, 1, 0], [2, 0, 0]]])
        with logger.block(tag='Internals 2nd Order dyad'):
            VPTRunner.run_simple(TestManager.test_data(file_name), 3, logger=logger, order=2, internals=[[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]], expansion_order=2, degeneracy_specs=[[[0, 1, 0], [2, 0, 0]]])
        with logger.block(tag='Internals 3rd Order dyad'):
            VPTRunner.run_simple(TestManager.test_data(file_name), 3, logger=logger, order=3, internals=[[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]], expansion_order=2, degeneracy_specs=[[[0, 1, 0], [2, 0, 0]]])
        with logger.block(tag='Cartesians 2nd Order triad'):
            VPTRunner.run_simple(TestManager.test_data(file_name), 3, logger=logger, order=2, expansion_order=2, degeneracy_specs=[[[0, 0, 1], [2, 0, 0]], [[0, 1, 0], [2, 0, 0]]])
        with logger.block(tag='Cartesians 3rd Order triad'):
            VPTRunner.run_simple(TestManager.test_data(file_name), 3, logger=logger, order=3, expansion_order=2, mixed_derivative_handling_mode=handling_mode, degeneracy_specs=[[[0, 0, 1], [2, 0, 0]], [[0, 1, 0], [2, 0, 0]]])
        with logger.block(tag='Cartesians 2nd Order dyad'):
            VPTRunner.run_simple(TestManager.test_data(file_name), 3, logger=logger, order=2, expansion_order=2, degeneracy_specs=[[[0, 1, 0], [2, 0, 0]]])
        with logger.block(tag='Cartesians 3rd Order dyad'):
            VPTRunner.run_simple(TestManager.test_data(file_name), 3, logger=logger, order=3, expansion_order=2, degeneracy_specs=[[[0, 1, 0], [2, 0, 0]]])
