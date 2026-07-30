"""Extracted from VPT2Tests.test_ResultsFileAnalysis via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VPT2Tests.test_ResultsFileAnalysis"""

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
    def test_ResultsFileAnalysis(self):
        temp_file = os.path.expanduser('~/Desktop/test_results.hdf5')
        log_file = os.path.expanduser('~/Desktop/test_results.txt')
        if not os.path.exists(temp_file):
            VPTRunner.run_simple(TestManager.test_data('i_hoh_opt.fchk'), 2, degeneracy_specs={'polyads': [[[0, 0, 0, 0, 1, 0], [0, 0, 0, 2, 0, 0]], [[0, 0, 0, 1, 0, 0], [0, 0, 2, 0, 0, 0]]], 'extra_groups': [[[0, 0, 0, 0, 1, 0], [0, 1, 0, 0, 1, 0], [1, 0, 0, 0, 1, 0], [0, 0, 0, 2, 0, 0], [0, 1, 0, 2, 0, 0], [1, 0, 0, 2, 0, 0], [0, 0, 2, 1, 0, 0], [0, 1, 2, 1, 0, 0], [1, 0, 2, 1, 0, 0], [0, 0, 4, 0, 0, 0], [0, 1, 4, 0, 0, 0], [1, 0, 4, 0, 0, 0]]]}, results=temp_file, logger=log_file, plot_spectrum=False)
        analyzer = VPTAnalyzer(temp_file)
        shifted_spec = analyzer.shifted_transformed_spectrum(analyzer.degenerate_states[4], analyzer.deperturbed_hamiltonians[4], [0, -50 / UnitsData.hartrees_to_wavenumbers])
        shifted_spec.plot()
        print(shifted_spec.frequencies, shifted_spec.intensities)
        with analyzer.log_parser as parser:
            for i, block in enumerate(parser.get_blocks()):
                for subblock in block.lines:
                    print(subblock.tag)
        from McUtils.Scaffolding import LogParser
        with LogParser(log_file) as parser:
            for i, block in enumerate(parser.get_blocks()):
                for subblock in block.lines:
                    print(subblock.tag)
