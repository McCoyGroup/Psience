"""Extracted from MolecoolsTests.test_LocalModeCHModel via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest MolecoolsTests.test_LocalModeCHModel"""

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
    def test_LocalModeCHModel(self):
        from Psience.BasisReps import modify_internal_hamiltonian
        methanol = Molecule.from_file(TestManager.test_data('methanol_vpt_3.fchk'))
        base_nms = methanol.get_normal_modes()
        internals = {(0, 1): 'OH', (1, 2): 'CO', (2, 3): 'CH3_stretch', (2, 4): 'CH3_stretch', (2, 5): 'CH3_stretch', (1, 2, 3): 'CH3_bend', (1, 2, 4): 'CH3_bend', (1, 2, 5): 'CH3_bend'}
        loc_modes = base_nms.localize(internals=internals)
        ob_modes = loc_modes.make_oblique()
        base_hess = ob_modes.local_hessian
        new_hess = modify_internal_hamiltonian(base_hess, {(0, 1): 'OH', (1, 2): 'CO', (2, 3): 'CH3_stretch', (2, 4): 'CH3_stretch', (2, 5): 'CH3_stretch', (1, 2, 3): 'CH3_bend', (1, 2, 4): 'CH3_bend', (1, 2, 5): 'CH3_bend'}, scaling_types={'OH': 0.93, 'CH3_stretch': 0.96}, coupling_types={(1, 'CH3_stretch', 'CH3_stretch'): -22 * UnitsData.convert('Wavenumbers', 'Hartrees')})
        g_mat = loc_modes.local_gmatrix

        def print_arr(header, array=None):
            if array is None:
                array = header
                header = []
            elif isinstance(header, str):
                header = [header]
            if array.ndim == 1:
                array = array[np.newaxis]
            print(*header, TableFormatter('{:.0f}').format(array * UnitsData.hartrees_to_wavenumbers))
        print()
        print('=' * 25, 'Unscaled', '=' * 25)
        freqs_old, _ = scipy.linalg.eigh(base_hess, g_mat, type=3)
        print_arr('Freqs:', np.sqrt(freqs_old))
        print('Hessian:')
        print_arr(base_hess)
        print('=' * 25, 'Scaled', '=' * 25)
        freqs_new, _ = scipy.linalg.eigh(new_hess, g_mat, type=3)
        print_arr('Freqs:', np.sqrt(freqs_new))
        print('Hessian:')
        print_arr(new_hess)
