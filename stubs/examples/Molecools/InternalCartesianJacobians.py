"""Extracted from MolecoolsTests.test_InternalCartesianJacobians via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest MolecoolsTests.test_InternalCartesianJacobians"""

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

    @inactiveTest
    def test_InternalCartesianJacobians(self):
        import McUtils.Plots as plt
        m = Molecule.from_file(TestManager.test_data('HOH_freq.fchk'), zmatrix=[[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]])
        intcds = m.internal_coordinates
        carts = m.coords
        ijacsnum, ijacs2num = intcds.jacobian(carts.system, [1, 2], all_numerical=True, mesh_spacing=0.01)
        ijacs, ijacs2 = intcds.jacobian(carts.system, [1, 2], analytic_deriv_order=1, converter_options=dict(reembed=False))
        jacs, jacs2 = carts.jacobian(intcds.system, [1, 2], mesh_spacing=1e-05)
        meh1 = ijacs.squeeze().reshape(9, 9)
        meh0 = ijacsnum.squeeze().reshape(9, 9)
        meh2 = jacs.squeeze().reshape(9, 9)
        itest = np.dot(meh1, meh2)
        itest2 = np.dot(meh2, meh1)
        self.assertTrue(np.allclose(np.eye(9), itest))
        good_sel = (...,) + np.ix_((3, 5, 6), (3, 5, 6))
        meh12 = ijacs2.squeeze().reshape(9, 9, 9)
        meh12 = meh12.transpose(2, 0, 1).reshape(3, 3, 9, 9)
        meh22 = ijacs2num.squeeze().reshape(9, 9, 9)
        meh22 = meh22.transpose(2, 0, 1).reshape(3, 3, 9, 9)
        meh12 = meh12[good_sel]
        meh22 = meh22[good_sel]
        ps = dict(vmin=-0.05, vmax=0.05)
        plt.TensorPlot(meh12, plot_style=ps)
        plt.TensorPlot(meh22, plot_style=ps).show()
        self.assertAlmostEquals(meh22[1, 1, 0, 0], 0.009235, places=6)
        self.assertTrue(np.allclose(meh12, meh22))
