"""Extracted from VPT2Tests.test_AnalyticWFC via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VPT2Tests.test_AnalyticWFC"""

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

    @debugTest
    def test_AnalyticWFC(self):
        import pickle
        file_name = 'OCHH_freq.fchk'
        corr = AnalyticVPTRunner.run_simple(TestManager.test_data(file_name), [[0, [[0, 0, 0, 0, 0, 1], [0, 1, 0, 1, 0, 0], [0, 0, 0, 1, 1, 0], [0, 0, 0, 0, 1, 0]]], [[0, 0, 0, 0, 1, 0], [[0, 0, 0, 0, 1, 1], [0, 1, 0, 1, 1, 0]]]], degeneracy_specs={'wfc_threshold': 0.3})
        with tempfile.NamedTemporaryFile() as tf:
            with open(tf.name, 'wb') as f:
                pickle.dump(corr, f)
            with open(tf.name, 'rb') as f:
                corr2 = pickle.load(f)
                print(corr2)
