"""Extracted from VPT2Tests.test_TermRephasing via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VPT2Tests.test_TermRephasing"""

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
    def test_TermRephasing(self):
        molg = Molecule.from_file(TestManager.test_data('ts_taup_anh_b2_qz_old.log'), 'gspec')
        no_dummy_pos = [i for i, a in enumerate(molg.atoms) if a != 'X']
        mol2 = molg.modify(atoms=np.array(molg.atoms)[no_dummy_pos,], coords=molg.coords[no_dummy_pos, :], potential_derivatives=molg.potential_derivatives, dipole_derivatives=molg.dipole_derivatives, internals=None)
        runner, _ = mol2.setup_VPT(states=1, degeneracy_specs='auto', logger=True, mode_selection=list(range(2, 12)))
        runner.print_tables(print_intensities=False)
