"""Extracted from MolecoolsTests.test_AIMNetExpansions via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest MolecoolsTests.test_AIMNetExpansions"""

import itertools
import os.path
import pprint
from Peeves.TestUtils import *
from unittest import TestCase
from Peeves import BlockProfiler
from Psience.Molecools import Molecule, MolecularNormalModes
from Psience.Data import DipoleSurface
from McUtils.GaussianInterface import GaussianFChkReader, GaussianLogReader
from McUtils.Plots import *
import McUtils.Plots as plt
from McUtils.Coordinerds import cartesian_to_zmatrix
from McUtils.Data import UnitsData
import numpy as np, scipy
import McUtils.Numputils as nput
import McUtils.Profilers as prof
from McUtils.Formatters import TableFormatter
import McUtils.Formatters as mfmt

class MolecoolsTests(TestCase):

    def setUp(self):
        self.test_log_water = TestManager.test_data('water_OH_scan.log')
        self.test_log_freq = TestManager.test_data('water_freq.log')
        self.test_HOD = TestManager.test_data('HOD_freq.fchk')
        self.test_fchk = TestManager.test_data('water_freq.fchk')
        self.test_log_h2 = TestManager.test_data('outer_H2_scan_new.log')

    def tearDown(self):
        ...

    @classmethod
    def setup_OCHH(cls, optimize=True):
        from McUtils.Extensions import ModuleLoader
        loader = ModuleLoader(os.path.expanduser('~/Documents/Postdoc/Projects/DGB'))
        h2co_mod = loader.load('H2COPot')

        def internal_pot(coords, order=None):
            coords = coords[..., (0, 1, 3, 2, 4, 5)]
            vals = h2co_mod.InternalsPotential.get_pot(coords)
            return vals
        ochh = Molecule.from_file(TestManager.test_data('OCHH_freq.fchk'), energy_evaluator={'potential_function': internal_pot, 'distance_units': 'Angstroms', 'energy_units': 'Wavenumbers', 'strip_embedding': True}, internals=[[0, -1, -1, -1], [1, 0, -1, -1], [2, 1, 0, -1], [3, 1, 2, 0]])
        if optimize:
            base_dip = ochh.dipole_derivatives
            ochh = ochh.optimize(method='conjugate-gradient', max_iterations=50, stencil=3, prevent_oscillations=3, restart_interval=15).modify(dipole_derivatives=base_dip)
        return ochh

    @validationTest
    def test_AIMNetExpansions(self):
        methanol = Molecule(['C', 'O', 'H', 'H', 'H', 'H'], [[-0.71174571, 0.0161939, 0.02050266], [1.71884591, -1.07310118, -0.2778059], [-1.30426891, 0.02589585, 1.99632677], [-0.77962613, 1.94036941, -0.7197672], [-2.02413643, -1.14525287, -1.05166036], [2.91548382, -0.08353621, 0.65084457]], energy_evaluator='aimnet2')
        expansion = methanol.calculate_energy(order=3)
