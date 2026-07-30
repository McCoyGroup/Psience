"""Extracted from MolecoolsTests.test_InternalProjectedModes via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest MolecoolsTests.test_InternalProjectedModes"""

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
    def test_InternalProjectedModes(self):
        import McUtils.Coordinerds as crops
        methanol_zmatrix = crops.functionalized_zmatrix(3, {(2, 1, 0): [[0, -1, -2, -3], [1, -1, 0, -2], [2, -1, 0, 1]]})
        methanol_zmatrix = crops.set_zmatrix_embedding(methanol_zmatrix)
        me_ints = Molecule.from_file(TestManager.test_data('methanol_vpt_1.fchk'), internals=methanol_zmatrix)
        nms = me_ints.get_normal_modes(project_transrot=False)
        locs = nms.localize(coordinate_constraints=crops.zmatrix_indices(methanol_zmatrix, [(3, 2, 1, 0)]))
        print()
        from McUtils.Formatters import TableFormatter
        print(TableFormatter('{:.0f}').format(nms.freqs[np.newaxis] * UnitsData.hartrees_to_wavenumbers))
        print(TableFormatter('{:.0f}').format(locs.local_hessian * UnitsData.hartrees_to_wavenumbers))
        me_carts = Molecule.from_file(TestManager.test_data('methanol_vpt_1.fchk'))
        nms_carts = me_carts.get_normal_modes(project_transrot=False, use_internals=False)
        loc_2 = nms_carts.apply_transformation(locs.localizing_transformation).make_dimensionless()
        f_nmw = me_carts.potential_derivatives[1]
        g12 = nput.fractional_power(me_carts.g_matrix, 1 / 2)
        f_mw = g12 @ f_nmw @ g12
        f_cart = nput.tensor_reexpand([loc_2.coords_by_modes], [0, f_mw])[-1]
        print(TableFormatter('{:.3f}').format(f_cart * UnitsData.hartrees_to_wavenumbers))
        return
        runner, _ = me_ints.setup_VPT(states=2, degeneracy_specs='auto', cartesian_analytic_deriv_order=-1, internal_by_cartesian_derivative_method='fast', modes=nms_carts, mode_transformation=locs.localizing_transformation)
        runner.print_tables()
