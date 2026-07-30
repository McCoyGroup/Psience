"""Extracted from MolecoolsTests.test_ZMatOpt via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest MolecoolsTests.test_ZMatOpt"""

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
    def test_ZMatOpt(self):
        import McUtils.Coordinerds as coordops
        import warnings
        warnings.filterwarnings('ignore', category=RuntimeWarning)
        z1 = Molecule.from_string('C(=O)O', 'smi', energy_evaluator='aimnet2', conf_get_options={'random_seed': 12321}).optimize(max_iterations=80)
        z1.internals = z1.get_bond_zmatrix()
        disp_coords = z1.get_displaced_coordinates([[np.pi / 6]], which=[-1], use_internals='reembed', strip_embedding=True)[0]
        z1 = z1.modify(coords=disp_coords)
        opt = z1.optimize(max_iterations=25)
        opt2 = z1.modify(internals=None).optimize(max_iterations=15)
        print(opt.calculate_energy() - z1.calculate_energy())
        print(opt.calculate_energy() - opt2.calculate_energy())
        print(coordops.set_zmatrix_embedding(np.round(z1.internal_coordinates - z1.modify(coords=opt.coords).internal_coordinates, 2)))
        print(coordops.set_zmatrix_embedding(np.round(z1.internal_coordinates - z1.modify(coords=opt2.coords).internal_coordinates, 2)))
