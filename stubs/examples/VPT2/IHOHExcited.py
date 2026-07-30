"""Extracted from VPT2Tests.test_IHOHExcited via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VPT2Tests.test_IHOHExcited"""

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
    def test_IHOHExcited(self):
        wfns = VPTRunner.run_simple(TestManager.test_data('i_hoh_opt.fchk'), VPTStateSpace.get_state_list_from_quanta(4, 6) + [[0, 1, 2, 2, 0, 0]], initial_states=[[0, 0, 0, 0, 0, 0], [0, 0, 0, 2, 0, 0], [0, 1, 0, 2, 0, 0], [0, 0, 0, 0, 1, 0]], degeneracy_specs={'polyads': [[[0, 0, 0, 0, 1, 0], [0, 0, 0, 2, 0, 0]], [[0, 0, 0, 1, 0, 0], [0, 0, 2, 0, 0, 0]]], 'extra_groups': [[[0, 0, 0, 0, 1, 0], [0, 1, 0, 0, 1, 0], [1, 0, 0, 0, 1, 0], [0, 0, 0, 2, 0, 0], [0, 1, 0, 2, 0, 0], [1, 0, 0, 2, 0, 0], [0, 0, 2, 1, 0, 0], [0, 1, 2, 1, 0, 0], [1, 0, 2, 1, 0, 0], [0, 0, 4, 0, 0, 0], [0, 1, 4, 0, 0, 0], [1, 0, 4, 0, 0, 0]]]}, target_property='wavefunctions', logger=os.path.expanduser('~/Desktop/specks/run_wfns.txt'), plot_spectrum=False)
        multispec = wfns.get_spectrum().frequency_filter(600, 4400)
        multispec.plot().savefig(os.path.expanduser(f'~/Desktop/specks/full.pdf'), transparent=True)
        for state, spec in zip(wfns.initial_states, multispec):
            s = ''.join((str(s) for s in state))
            spec.plot(plot_range=[[600, 4400], [0, 500]], padding=[[0, 0], [0, 0]], image_size=[0.75 * 20, 5.48 * 20]).savefig(os.path.expanduser(f'~/Desktop/specks/state_{s}.pdf'), transparent=True)
        multispec = wfns.get_deperturbed_spectrum().frequency_filter(600, 4400)
        for state, spec in zip(wfns.initial_states, multispec):
            s = ''.join((str(s) for s in state))
            spec.plot(plot_range=[[600, 4400], [0, 500]], padding=[[0, 0], [0, 0]], image_size=[0.75 * 20, 5.48 * 20]).savefig(os.path.expanduser(f'~/Desktop/specks/state_depert_{s}.pdf'), transparent=True)
