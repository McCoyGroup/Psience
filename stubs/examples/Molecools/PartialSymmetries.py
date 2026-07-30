"""Extracted from MolecoolsTests.test_PartialSymmetries via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest MolecoolsTests.test_PartialSymmetries"""

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
    def test_PartialSymmetries(self):
        import McUtils.Coordinerds as coordops
        nh3 = Molecule.from_file(TestManager.test_data('nh3.fchk'), internals=coordops.spoke_zmatrix(3)).get_embedded_molecule()
        coeffs, basis, expansions = nh3.get_symmetrized_internals(return_expansions=True)
        cpmo3m = Molecule.from_file(TestManager.test_data('cpmo3m_opt.xyz'), units='Angstroms')
        zmat = coordops.reindex_zmatrix(coordops.functionalized_zmatrix(1, {(0, -1, -2): coordops.chain_zmatrix(2), (0, 1, 2): coordops.chain_zmatrix(2), (0, 3, 4): coordops.chain_zmatrix(2), (0, 5, 6): coordops.chain_zmatrix(2), (0, 7, 8): coordops.chain_zmatrix(2), (0, 9, 10): coordops.chain_zmatrix(2), (0, 11, 12): coordops.chain_zmatrix(2), (0, 13, 14): coordops.chain_zmatrix(2)}), [10, 11, 12, 13, 14, 15, 16, 0, 9, 1, 8, 2, 7, 3, 6, 4, 5])
        cpmo3m_int = cpmo3m.modify(internals=zmat)
        coeffs, intenrals, expansions = cpmo3m_int.get_symmetrized_internals(atom_selection=list(range(11)), tol=0.8, permutation_tol=0.9, return_expansions=True, return_point_group=False, extra_internals=[(0, 1)])
        a1_modes, a1_inv = expansions[0]
        cpmo3m_int.animate_coordinate(0, coordinate_expansion=[a1_inv[0][(15,),]]).show()
