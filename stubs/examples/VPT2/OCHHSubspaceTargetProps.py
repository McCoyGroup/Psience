"""Extracted from VPT2Tests.test_OCHHSubspaceTargetProps via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VPT2Tests.test_OCHHSubspaceTargetProps"""

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
    def test_OCHHSubspaceTargetProps(self):
        file_name = 'OCHH_freq.fchk'
        VPTRunner.run_simple(TestManager.test_data(file_name), [[0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 1], [0, 1, 2, 0, 0, 0]], logger=True, target_property='wavefunctions', degeneracy_specs='auto')
        '\nState           Harmonic   Anharmonic     Harmonic   Anharmonic\n                     ZPE          ZPE    Frequency    Frequency\n0 0 0 0 0 0   5866.87156   5785.95861            -            - \n0 0 1 0 0 1            -            -   4588.74227   4415.63305 \n0 1 2 0 0 0            -            -   4306.24556   4191.96329 \n'
        '\n::> IR Data\n  > Initial State: 0 0 0 0 0 0 \n                         Harmonic                  Anharmonic\nState             Frequency    Intensity       Frequency    Intensity\n  0 0 1 0 0 1    4588.74227      0.00000      4446.79734      0.00164\n  0 1 2 0 0 0    4306.24556      0.00000      4167.25991      0.30489\n  0 1 1 1 0 0    4506.28742      0.00000      4345.62126      0.48928\n  '
        VPTRunner.run_simple(TestManager.test_data(file_name), [[0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 1], [0, 1, 2, 0, 0, 0]], logger=True, target_property='intensities', degeneracy_specs='auto')
        '\n::> IR Data\n  > Initial State: 0 0 0 0 0 0 \n                         Harmonic                  Anharmonic\nState             Frequency    Intensity       Frequency    Intensity\n  0 0 1 0 0 1    4588.74227      0.00000      4446.79734      0.00164\n  0 1 2 0 0 0    4306.24556      0.00000      4167.25991      0.30489\n  0 1 1 1 0 0    4506.28742      0.00000      4345.62126      0.48928\n  '
        VPTRunner.run_simple(TestManager.test_data(file_name), [[0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 1], [0, 1, 2, 0, 0, 0]], logger=True, target_property='frequencies', degeneracy_specs='auto')
        '\n::> IR Data\n  > Initial State: 0 0 0 0 0 0 \n                         Harmonic                  Anharmonic\nState             Frequency    Intensity       Frequency    Intensity\n  0 0 1 0 0 1    4588.74227      0.00000      4458.35353      0.00402\n  0 1 2 0 0 0    4306.24556      0.00000      4170.69772      0.23905\n  0 1 1 1 0 0    4506.28742      0.00000      4330.62726      0.55570\n  '
