"""Extracted from VPT2Tests.test_HOHLocal via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VPT2Tests.test_HOHLocal"""

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
    def test_HOHLocal(self):
        file_name = 'HOH_freq.fchk'
        OHH = Molecule.from_file(TestManager.test_data(file_name), internals=[[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]])
        XR = OHH.get_internals_by_cartesians(1, strip_embedding=True)[0]
        RX = OHH.get_cartesians_by_internals(1, strip_embedding=True)[0]
        VPTRunner.run_simple(TestManager.test_data(file_name), 1, local_modes={'matrix': XR, 'inverse': RX, 'sort_freqs': True}, mixed_derivative_handling_mode='analytical', local_mode_couplings=True, degeneracy_specs={'polyads': [[[0, 0, 1], [0, 1, 0]], [[0, 0, 1], [2, 0, 0]]]}, calculate_intensities=True)
        '\n        ::> Deperturbed IR Data\n          > Initial State: 0 0 0 \n                           Harmonic                  Anharmonic\n        State       Frequency    Intensity       Frequency    Intensity\n          0 0 1    3873.50774     35.31527      3688.68476     35.44348\n          0 1 0    3873.50774     35.31527      3688.97960     35.28916\n          1 0 0    1641.37852     65.14715      1597.32063     67.22674\n          2 0 0    3282.75704      0.00000      3183.94992      0.15040\n        ::> IR Data\n          > Initial State: 0 0 0 \n                           Harmonic                  Anharmonic\n        State       Frequency    Intensity       Frequency    Intensity\n          0 0 1    3873.50774     35.31527      3636.87729      4.18690\n          0 1 0    3873.50774     35.31527      3752.80408     67.80845\n          1 0 0    1641.37852     65.14715      1597.32063     67.22674\n          2 0 0    3282.75704      0.00000      3171.93290      0.00656\n        '
        return
        VPTRunner.run_simple(TestManager.test_data(file_name), 1, mode_selection=[0, 2, 1], mixed_derivative_handling_mode='analytical', degeneracy_specs={'polyads': [[[0, 0, 1], [0, 1, 0]], [[0, 0, 1], [2, 0, 0]]]}, calculate_intensities=True)
        '\n        ::> Deperturbed IR Data\n          > Initial State: 0 0 0 \n                           Harmonic                  Anharmonic\n        State       Frequency    Intensity       Frequency    Intensity\n          0 0 1    3803.29959      4.14283      3612.43457      3.45008\n          0 1 0    3937.52466     67.02051      3744.73518     64.43622\n          1 0 0    1622.30304     67.45626      1572.70760     68.11278\n          2 0 0    3244.60608      0.00000      3126.93384      1.18410\n        ::> IR Data\n          > Initial State: 0 0 0 \n                           Harmonic                  Anharmonic\n        State       Frequency    Intensity       Frequency    Intensity\n          0 0 1    3803.29959      4.14283      3623.17785      2.78848\n          0 1 0    3937.52466     67.02051      3744.73518     64.43622\n          1 0 0    1622.30304     67.45626      1572.70760     68.11278\n          2 0 0    3244.60608      0.00000      3116.19056      1.75788\n        '
        raise Exception(...)
        runner, states = VPTRunner.construct(OHH, 2, degeneracy_specs={'polyads': [[[0, 0, 1], [0, 1, 0]]]})
        classic, _ = runner.construct_classic_runner(states)
        extra_coupling_ham = np.diag(runner.eval.freqs) + np.array([[0, 100, 100], [100, 0, 100], [100, 100, 0]]) * UnitsData.convert('Wavenumbers', 'Hartrees')
        runner.run_VPT(states, calculate_intensities=False, degenerate_zero_order_hamiltonian=extra_coupling_ham)
        raise Exception(...)
