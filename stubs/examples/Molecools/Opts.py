"""Extracted from MolecoolsTests.test_Opts via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest MolecoolsTests.test_Opts"""

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
    def test_Opts(self):
        import McUtils.Coordinerds as coordops
        import McUtils.Formatters as mfmt
        import numpy as np
        np.seterr(all='raise')
        mol = Molecule.from_file(TestManager.test_data('OCHH_freq.fchk'), energy_evaluator='aimnet2')
        fig = mol.plot(use_default_bonds=False)
        mol.plot(use_default_bonds=False, figure=fig)
        return
        mol = Molecule.from_file(TestManager.test_data('react_samp.xyz'), energy_evaluator='rdkit')
        coordops.validate_zmatrix(mol.get_bond_zmatrix())
        zmcs = mol.get_bond_zmatrix()
        n = 4
        constraints = [z for z in coordops.extract_zmatrix_internals(zmcs) if len(z) == n][:3]
        opt = mol.optimize(coordinate_constraints=constraints)
        opt2 = mol.optimize()
        print(constraints)
        idx = coordops.zmatrix_indices(zmcs, constraints)
        pop_vals = lambda m: coordops.extract_zmatrix_values(m.modify(internals=zmcs).internal_coordinates, idx)
        print(pop_vals(mol))
        print(pop_vals(opt))
        print(pop_vals(opt2))
        return
        zmat = mol.get_bond_zmatrix(validate=True)
        int_mol = mol.modify(internals=zmat)
        print(mol.coords[:5])
        exp = int_mol.get_cartesians_by_internals(method='classic', use_direct_expansions=True, orthogonalize_derivatives=False, allow_fd=False, order=1, strip_embedding=False)
        exp = mol.modify(internals=zmat).get_cartesians_by_internals(allow_fd=False, order=1, strip_embedding=True)
        return
        mol = Molecule.from_file(TestManager.test_data('tbhp_180.fchk'), energy_evaluator='aimnet2')
        ugh = mol.modify(internals=zmat).optimize(coordinate_constraints=[c for c in coordops.extract_zmatrix_internals(zmat) if len(c) == 4], track_best=True, max_iterations=50)
        print(np.concatenate([ugh.modify(internals=zmat).internal_coordinates, mol.modify(internals=zmat).internal_coordinates], axis=-1))
        print((ugh.calculate_energy() - mol.calculate_energy()) * UnitsData.convert('Hartrees', 'Kilocalories/Mole'))
