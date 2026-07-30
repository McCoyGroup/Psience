"""Extracted from MolecoolsTests.test_1DPotentialReps via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest MolecoolsTests.test_1DPotentialReps"""

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
    def test_1DPotentialReps(self):
        ochh = self.setup_OCHH(optimize=True)
        int_ochh = ochh.modify(internals=[[0, -1, -1, -1], [1, 0, -1, -1], [2, 1, 0, -1], [3, 1, 2, 0]])
        scan_disps = [-1.0, 1.0, 51]
        scan_angles = int_ochh.get_scan_coordinates([scan_disps], which=[[3, 1]], internals='reembed')
        pot_vals = int_ochh.calculate_energy(scan_angles)
        woof_pots = ochh.get_1d_potentials([(0, 1), (1, 2), (1, 3), (2, 1, 3), (0, 1, 2, 3)])
        woof_pots_4 = ochh.get_1d_potentials([(0, 1), (1, 2), (1, 3), (2, 1, 3), (0, 1, 2, 3)], poly_expansion_order=4)
        woof_pots_morse = ochh.get_1d_potentials([(0, 1), (1, 2), (1, 3), (2, 1, 3), (0, 1, 2, 3)], quartic_potential_cutoff=0)
        angs, pot_vals_angs_appx = woof_pots[3](scan_angles)
        angs, pot_vals_angs_appx_4 = woof_pots_4[3](scan_angles)
        angs, pot_vals_angs_appx_morse = woof_pots_morse[3](scan_angles)
        disp_vals = np.linspace(*scan_disps)
        ploots = plt.Plot(disp_vals, pot_vals * 219475.6, color='black', plot_range=[None, [None, 35000]])
        shang_vals = ochh.calculate_energy() + pot_vals_angs_appx[0]
        shang_vals_4 = ochh.calculate_energy() + pot_vals_angs_appx_4[0]
        shang_vals_morse = ochh.calculate_energy() + pot_vals_angs_appx_morse[0]
        plt.Plot(np.linspace(*scan_disps), shang_vals * 219475.6, figure=ploots, linestyle='dashed')
        plt.Plot(np.linspace(*scan_disps), shang_vals_4 * 219475.6, figure=ploots, linestyle='dashed')
        plt.Plot(np.linspace(*scan_disps), shang_vals_morse * 219475.6, figure=ploots, linestyle='dashed')
        ploots.show()
        return
