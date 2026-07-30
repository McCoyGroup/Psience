"""Extracted from MolecoolsTests.test_AddDummyAtomJacobians via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest MolecoolsTests.test_AddDummyAtomJacobians"""

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
    def test_AddDummyAtomJacobians(self):
        file_name = TestManager.test_data('HOONO_freq.fchk')
        mol = Molecule.from_file(file_name)
        n_pos = mol.atom_positions['N']
        o_pos = mol.atom_positions['O']
        normal = nput.vec_crosses(mol.coords[o_pos[0]] - mol.coords[o_pos[1]], mol.coords[n_pos[0]] - mol.coords[o_pos[1]], normalize=True)
        mol2 = mol.insert_atoms('X', mol.coords[o_pos[1]] + 5 * normal, 3, handle_properties=False)
        mol2.zmatrix = [[1, -1, -1, -1], [2, 1, -1, -1], [3, 2, 1, -1], [5, 2, 1, 3], [0, 1, 2, 5], [4, 3, 2, 5]]
        jacobians_no_dummy = mol2.coords.jacobian(mol2.internal_coordinates.system, [1, 2], stencil=3, all_numerical=True, converter_options=dict(strip_dummies=True))
        self.assertEquals(jacobians_no_dummy[0].shape, (5, 3, 5, 3))
        self.assertEquals(jacobians_no_dummy[1].shape, (5, 3, 5, 3, 5, 3))
        jacobians = mol2.coords.jacobian(mol2.internal_coordinates.system, [1, 2, 3], stencil=5, all_numerical=True, converter_options=dict(strip_dummies=False))
        self.assertEquals(jacobians[0].shape, (6, 3, 6, 3))
        self.assertEquals(jacobians[1].shape, (6, 3, 6, 3, 6, 3))
        self.assertEquals(jacobians[2].shape, (6, 3, 6, 3, 6, 3, 6, 3))
        jacobians_analytic = mol2.coords.jacobian(mol2.internal_coordinates.system, [1, 2], stencil=5, analytic_deriv_order=1, converter_options=dict(strip_dummies=False))
        self.assertEquals(jacobians_analytic[0].shape, (6, 3, 6, 3))
        self.assertEquals(jacobians_analytic[1].shape, (6, 3, 6, 3, 6, 3))
        jacobians_no_dummy_analytic = mol2.coords.jacobian(mol2.internal_coordinates.system, [1, 2], stencil=3, analytic_deriv_order=1, converter_options=dict(strip_dummies=True))
        self.assertEquals(jacobians_no_dummy_analytic[0].shape, (5, 3, 5, 3))
        self.assertEquals(jacobians_no_dummy_analytic[1].shape, (5, 3, 5, 3, 5, 3))
        self.assertTrue(np.allclose(jacobians[0][0, 0][:2], jacobians_no_dummy[0][0, 0][:2]))
        self.assertTrue(np.allclose(jacobians[1][0, 0, 0, 0][:2], jacobians_no_dummy[1][0, 0, 0, 0][:2]))
        jacobians_no_dummy = mol2.internal_coordinates.jacobian(mol2.coords.system, [1, 2], stencil=3, all_numerical=True, converter_options=dict(strip_dummies=True))
        self.assertEquals(jacobians_no_dummy[0].shape, (5, 3, 5, 3))
        self.assertEquals(jacobians_no_dummy[1].shape, (5, 3, 5, 3, 5, 3))
        jacobians = mol2.internal_coordinates.jacobian(mol2.coords.system, [1, 2], stencil=3, all_numerical=True, converter_options=dict(strip_dummies=False))
        self.assertEquals(jacobians[0].shape, (6, 3, 6, 3))
        self.assertEquals(jacobians[1].shape, (6, 3, 6, 3, 6, 3))
