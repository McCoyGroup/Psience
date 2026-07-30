"""Extracted from MolecoolsTests.test_DihedOptRPNMScan via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest MolecoolsTests.test_DihedOptRPNMScan"""

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
    def test_DihedOptRPNMScan(self):
        import McUtils.Coordinerds as coordops
        from Psience.Modes import ReactionPathModes
        methanol = Molecule.from_string('methanol', energy_evaluator='aimnet2').optimize(max_iterations=250)
        methanol_zmatrix = coordops.reindex_zmatrix(coordops.functionalized_zmatrix(3, {(2, 1, 0): [[0, -1, -2, -3], [1, -1, 0, -2], [2, -1, 0, 1]]}), [5, 1, 0, 2, 3, 4])
        me_int = methanol.modify(internals=methanol_zmatrix)
        scan_coord = [(2, 0, 1, 5)]
        scan_spec = [-np.deg2rad(180), np.deg2rad(0), 6]
        me_coords = me_int.get_scan_coordinates([scan_spec], which=coordops.zmatrix_indices(methanol_zmatrix, scan_coord), internals='reembed', strip_embedding=True)
        opt_mols = [methanol.modify(coords=c).optimize(coordinate_constraints=scan_coord) for c in me_coords]
        opt_coords = np.array([o.coords for o in opt_mols])
        me_aimnet2 = methanol.get_energy_function('aimnet2')
        me_expansions = [me_aimnet2(c, order=2) for c in opt_coords]
        me_grads = [exp[1] for exp in me_expansions]
        me_hesss = [exp[2] for exp in me_expansions]
        me_proj = nput.translation_rotation_projector(opt_coords, methanol.atomic_masses, mass_weighted=True)
        rpnms_opt = ReactionPathModes.get_rp_modes(me_grads, me_hesss, methanol.atomic_masses, projector=me_proj)
