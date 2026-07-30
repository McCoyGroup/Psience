"""Extracted from VPT2Tests.test_HOHPartialQuartic via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VPT2Tests.test_HOHPartialQuartic"""

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
    def test_HOHPartialQuartic(self):
        VPTRunner.helpers.run_anne_job(TestManager.test_data('vpt2_helpers_api/hod/x'), states=[[0, 0, 0], [0, 0, 1], [0, 1, 0]], order=2, expansion_order=2, calculate_intensities=False, zero_element_warning=False, include_only_mode_couplings=[1, 2], include_coriolis_coupling=False)
        VPTRunner.helpers.run_anne_job(TestManager.test_data('vpt2_helpers_api/hod/x_no_bend'), states=[[0, 0, 0], [0, 0, 1], [0, 1, 0]], order=2, expansion_order=2, calculate_intensities=False, zero_element_warning=False, include_coriolis_coupling=True)
        raise Exception(...)
        VPTRunner.helpers.run_anne_job(TestManager.test_data('vpt2_helpers_api/hod/x_decoupled'), states=1, order=2, expansion_order=2, calculate_intensities=False, zero_element_warning=False, include_coriolis_coupling=False)
        '\n::> States Energies\n  > State     Harmonic   Anharmonic     Harmonic   Anharmonic\n               ZPE          ZPE    Frequency    Frequency\n0 0 0   4052.91097   4001.04707            -            - \n0 0 1            -            -   3873.84521   3688.94993 \n0 1 0            -            -   2810.03028   2723.42011 \n1 0 0            -            -   1421.94645   1383.13454 \n'
        VPTRunner.helpers.run_anne_job(TestManager.test_data('vpt2_helpers_api/hod/x'), states=1, order=2, expansion_order=2, mode_selection=[0, 2], calculate_intensities=False, zero_element_warning=False, include_coriolis_coupling=False)
        '\n0 0 0   4052.91097   3994.84632            -            - \n0 0 1            -            -   3873.84521   3685.79215 \n0 1 0            -            -   2810.03028   2706.14630 \n1 0 0            -            -   1421.94645   1383.40761 \n\n0 0 0   4052.91097   4033.89584            -            - \n0 0 1            -            -   3873.84521   3697.97351 \n0 1 0            -            -   2810.03028   2902.77195 \n1 0 0            -            -   1421.94645   1380.70462\n'
        '\n        array([1383.40761359, 2706.14629421, 3685.79215281])\n        array([1380.70462067, 2902.77194309, 3697.97351711]))\n        '
        raise Exception(...)
        raise Exception(...)
        '\n        0   1973.85245   1992.46835            -            - \n        1            -            -   3947.70491   4042.24072\n        '
        runner1 = VPTRunner.helpers.run_anne_job(os.path.expanduser('~/Desktop/r_as'), states=1, order=2, expansion_order=2, calculate_intensities=False, operator_coefficient_threshold=-1)
        '\n  State       Frequency    Intensity       Frequency    Intensity\n  0 0 1    3947.69802     75.46200      3761.26977     70.72838\n  0 1 0    3821.87392      5.56098      3652.31667      4.82837\n  1 0 0    1628.37574      0.00000      1590.97610      0.00483\n  State       Frequency    Intensity       Frequency    Intensity\n  0 0 1    3947.69802     75.46200      3874.95435     71.78714\n  0 1 0    3821.87392      0.00000      3864.33291      0.00111\n  1 0 0    1628.37574      0.00000      1620.53058      0.00003\n'
        '\n        ::> States Energies\n  > State     Harmonic   Anharmonic     Harmonic   Anharmonic\n               ZPE          ZPE    Frequency    Frequency\n0 0 0   4698.97384   4628.17891            -            - \n0 0 1            -            -   3947.69802   3761.26977 \n0 1 0            -            -   3821.87392   3652.31667 \n1 0 0            -            -   1628.37574   1590.97610 \n'
        runner2 = VPTRunner.helpers.run_anne_job(os.path.expanduser('~/Desktop/r_a'), states=1, order=2, expansion_order=2, calculate_intensities=False, operator_coefficient_threshold=-1)
        raise Exception(...)
        raise Exception(runner1[0].ham_opts.opts['potential_terms'][1], runner2[0].ham_opts.opts['potential_terms'][1])
