"""Extracted from MolecoolsTests.test_PartialQuartic via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest MolecoolsTests.test_PartialQuartic"""

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
    def test_PartialQuartic(self):
        import McUtils.Coordinerds as coordops
        methanol_zmatrix = coordops.reindex_zmatrix(coordops.functionalized_zmatrix(3, {(2, 1, 0): [[0, -1, -2, -3], [1, -1, 0, -2], [2, -1, 0, 1]]}), [5, 1, 0, 2, 3, 4])
        crds_aimnet_opt = np.array([[-0.71174571, 0.0161939, 0.02050266], [1.71884591, -1.07310118, -0.2778059], [-1.30426891, 0.02589585, 1.99632677], [-0.77962613, 1.94036941, -0.7197672], [-2.02413643, -1.14525287, -1.05166036], [2.91548382, -0.08353621, 0.65084457]])
        crds_bad = np.array([[[[0.99603661, -0.03076131, 0.317282], [2.28225031, -0.60719144, 0.1594239], [0.68719378, -0.0280038, 1.36425206], [0.95774141, 0.98826296, -0.07215449], [0.29937143, -0.64460867, -0.24823604], [2.91548382, -0.08353621, 0.65084457]]], [[[0.99603661, -0.03076131, 0.317282], [2.28225031, -0.60719144, 0.1594239], [0.68248684, -0.02562726, 1.36284309], [0.96011584, 0.98746852, -0.07445194], [0.30154935, -0.64537247, -0.25008224], [2.91548382, -0.08353621, 0.65084457]]], [[[0.99603661, -0.03076131, 0.317282], [2.28225031, -0.60719144, 0.1594239], [0.67778774, -0.02325085, 1.36140798], [0.9625001, 0.98666785, -0.07673702], [0.3037351, -0.64614144, -0.25191702], [2.91548382, -0.08353621, 0.65084457]]]]) / UnitsData.bohr_to_angstroms
        methanol = Molecule(['C', 'O', 'H', 'H', 'H', 'H'], crds_aimnet_opt, energy_evaluator='rdkit').optimize()
        woof = methanol.partial_force_field(order=3, modes=methanol.get_normal_modes(project_transrot=False))
        print(woof[-1][0, 0])
        print(woof[-1][0, 1])
        return
        proj_dir = os.path.expanduser('~/Documents/Postdoc/Projects/CoordinatePaper/ml_fd_tests/')
        os.makedirs(proj_dir, exist_ok=True)
        from McUtils.Formatters import TableFormatter
        modes = methanol.get_normal_modes()
        for step_size in [1]:
            for stencil in [5]:
                print(f'Mesh spacing: {step_size}')
                print(f'Stencil: {stencil}')
                print(f'Freqs: {modes.freqs * UnitsData.hartrees_to_wavenumbers}')
                expansion = methanol.partial_force_field(order=3, modes=modes, mesh_spacing=step_size, stencil=stencil)
                terms = expansion[-1]
                for term in terms:
                    print('=' * 100)
                    print(TableFormatter('{:.3f}').format(term * UnitsData.hartrees_to_wavenumbers))
