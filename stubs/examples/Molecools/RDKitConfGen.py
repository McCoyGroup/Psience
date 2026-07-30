"""Extracted from MolecoolsTests.test_RDKitConfGen via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest MolecoolsTests.test_RDKitConfGen"""

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
    def test_RDKitConfGen(self):
        from Psience.Molecools import Molecule
        mols = Molecule.from_string('Cc1ccc(OC[C:6]([NH:5]/[N:4]=[CH:3]/[c:2]2cc(=O)[nH]c(=O)[nH:1]2)=[O:7])c([N+](=O)[O-])c1', 'smi', confgen_opts={'distance_constraints': {(4, 5): [1.2, 1.4]}, 'random_seed': 100})
        mols = Molecule.from_string('Cc1ccc(OC[C:6]([NH:5]/[N:4]=[CH:3]/[c:2]2cc(=O)[nH]c(=O)[nH:1]2)=[O:7])c([N+](=O)[O-])c1', 'smi', num_confs=10, confgen_opts={'distance_constraints': {(4, 5): [1.2, 1.4]}, 'random_seed': 100})
        raise Exception(mols)
        print(Molecule(['O', 'H', 'H'], [[0, 0, 0], [0, 1, 0], [1, 0, 0]]).bonds)
        mol = Molecule.from_string('Cc1ccc(OC[C:6]([NH:5]/[N:4]=[CH:3]/[c:2]2cc(=O)[nH]c(=O)[nH:1]2)=[O:7])c([N+](=O)[O-])c1', 'smi', confgen_opts={'distance_constraints': {(4, 5): [1.2, 1.4]}, 'random_seed': 100})
        print(Molecule(mol.atoms, mol.coords).bonds[:3])
        print(np.linalg.norm(mol.coords[4] - mol.coords[5]) * UnitsData.bohr_to_angstroms)
        Molecule.from_string('NC(N)=NC(=O)c1[nH]c(C(=O)O)c2cc3c(cc12)c1c2cc4c(C(=O)O)[nH]c(C(=O)N=C(N)N)c4cc2c2c4cc5c(C(=O)O)[nH:1][c:2]([C:3](=O)[N:4]=[C:5](N)[NH2:6])c5cc4c4c5cc6c(C(=O)O)[nH]c(C(=O)N=C(N)N)c6cc5c5c6cc7c(C(=O)O)[nH]c(C(=O)N=C(N)N)c7cc6c3c3c1c2c4c53', 'smi')
