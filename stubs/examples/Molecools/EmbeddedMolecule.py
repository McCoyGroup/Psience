"""Extracted from MolecoolsTests.test_EmbeddedMolecule via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest MolecoolsTests.test_EmbeddedMolecule"""

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
    def test_EmbeddedMolecule(self):
        file_name = TestManager.test_data('HOH_freq.fchk')
        mol1 = Molecule.from_file(file_name)
        mol = mol1.get_embedded_molecule()
        init_mat = mol1.normal_modes.modes.basis.matrix
        self.assertTrue(np.allclose(mol.moments_of_inertia, mol1.moments_of_inertia), msg='(HOH) Moments of inertia changed post-rotation: {} to {}'.format(mol1.moments_of_inertia, mol.moments_of_inertia))
        self.assertTrue(np.allclose(mol.inertial_axes, np.eye(3)), msg='(HOH) Principle axes are not identity matrix post-rotation: {}'.format(mol.inertial_axes))
        norms_1 = np.linalg.norm(mol.normal_modes.modes.basis.matrix, axis=0)
        norms_2 = np.linalg.norm(init_mat, axis=0)
        self.assertTrue(np.allclose(norms_1, norms_2), msg='(HOH) Normal modes renomalized:{} different from {}'.format(norms_1, norms_2))
        file_name = TestManager.test_data('tbhp_180.fchk')
        mol1 = Molecule.from_file(file_name)
        mol = mol1.get_embedded_molecule()
        init_mat = mol1.normal_modes.modes.basis.matrix
        self.assertTrue(np.allclose(mol.moments_of_inertia, mol1.moments_of_inertia), msg='(TBHP) Moments of inertia changed post-rotation: {} to {}'.format(mol1.moments_of_inertia, mol.moments_of_inertia))
        self.assertTrue(np.allclose(mol.inertial_axes, np.eye(3)), msg='(TBHP) Principle axes are not identity matrix post-rotation: {}'.format(mol.inertial_axes))
        norms_1 = np.linalg.norm(mol.normal_modes.modes.basis.matrix, axis=0)
        norms_2 = np.linalg.norm(init_mat, axis=0)
        self.assertTrue(np.allclose(norms_1, norms_2), msg='(TBHP) Normal modes renomalized: {} different from {}'.format(norms_1, norms_2))
        file_name = TestManager.test_data('HOONO_freq.fchk')
        mol1 = Molecule.from_file(file_name)
        mol = mol1.get_embedded_molecule()
        init_mat = mol1.normal_modes.modes.basis.matrix
        self.assertTrue(np.allclose(mol.moments_of_inertia, mol1.moments_of_inertia), msg='(HOONO) Moments of inertia changed post-rotation: {} to {}'.format(mol1.moments_of_inertia, mol.moments_of_inertia))
        self.assertTrue(np.allclose(mol.inertial_axes, np.eye(3)), msg='(HOONO) Principle axes are not identity matrix post-rotation: {}'.format(mol.inertial_axes))
        norms_1 = np.linalg.norm(mol.normal_modes.modes.basis.matrix, axis=0)
        norms_2 = np.linalg.norm(init_mat, axis=0)
        self.assertTrue(np.allclose(norms_1, norms_2), msg='(HOONO) Normal modes renomalized: {} different from {}'.format(norms_1, norms_2))
