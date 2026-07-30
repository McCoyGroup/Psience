"""Extracted from MolecoolsTests.test_InternalConv via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest MolecoolsTests.test_InternalConv"""

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
    def test_InternalConv(self):
        test_inverse = False
        test_vpt = False
        test_direct_expansions = False
        test_numerical_expansions = ['all']
        nh3 = Molecule.from_file(TestManager.test_data('nh3.fchk'), internals=[[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1], [3, 0, 1, 2]])
        nh3_2 = nh3.copy()
        test_coords = nh3.coords
        if test_direct_expansions:
            with np.printoptions(linewidth=100000000.0, suppress=True):
                e1 = nput.angle_vec(test_coords, 0, 1, 2, method='old', angle_ordering='jik', order=2)
                e2 = nput.angle_vec(test_coords, 0, 1, 2, method='expansion', angle_ordering='jik', order=2)
                for i1, i2 in zip(e1, e2):
                    print(':::' * 50)
                    print(':::' * 50)
                    print(i1.shape)
                    w = np.where(np.abs(i1 - i2) > 1e-06)
                    print(w)
                    if len(w[0]) > 0:
                        print('=' * 50)
                        print(np.round(i1, 6))
                        print('-' * 50)
                        print(np.round(i2, 6))
                        print('=' * 50)
            e1 = nput.dihed_vec(test_coords, 3, 0, 1, 2, method='old', order=2)
            e2 = nput.dihed_vec(test_coords, 3, 0, 1, 2, method='expansion', order=2)
            with np.printoptions(linewidth=100000000.0, suppress=True):
                for i1, i2 in zip(e1, e2):
                    print(':::' * 50)
                    print(':::' * 50)
                    print(i1.shape)
                    w = np.where(np.abs(i1 - i2) > 1e-06)
                    print(w)
                    if len(w[0]) > 0:
                        print('=' * 50)
                        print(np.round(i1, 6))
                        print('-' * 50)
                        print(np.round(i2, 6))
                        print('=' * 50)
        if test_numerical_expansions:
            if isinstance(test_numerical_expansions, bool):
                test_numerical_expansions = ['dihedrals', 'all']
            exp0 = nh3.get_internals_by_cartesians(4, strip_embedding=True, analytic_derivative_order=0, all_numerical=True)
            if 'dihedrals' in test_numerical_expansions:
                with np.printoptions(linewidth=100000000.0, suppress=True):
                    e1 = [e[..., -1] for e in exp0]
                    e2 = nput.dihed_vec(test_coords, *nh3.internals['zmatrix'][-1], method='old', order=2)[1:]
                    for i1, i2 in zip(e1, e2):
                        print(':::' * 50)
                        print(':::' * 50)
                        print(i1.shape)
                        i1[np.abs(i1) < 1e-06] = 0
                        i2[np.abs(i2) < 1e-06] = 0
                        d1 = -i1.copy()
                        d1[np.abs(d1) < 1e-06] = 1
                        r = i2 / d1
                        r[np.logical_and(np.abs(i1) < 1e-06, np.abs(i2) < 1e-06)] = 1
                        w = np.where(np.abs(r - 1) > 0.01)
                        r[np.logical_and(np.abs(i1) < 1e-06, np.abs(i2) < 1e-06)] = 0
                        print(w)
                        if len(w[0]) > 0:
                            print('=' * 50)
                            print(np.round(i1, 6)[2])
                            print('.' * 10)
                            print(np.round(i2, 6)[2])
                            print('=' * 50)
                            print(np.round(r, 6))
            if 'all' in test_numerical_expansions:
                exp1 = nh3_2.get_internals_by_cartesians(4, strip_embedding=True, analytic_derivative_order=-1, all_numerical=False)
                for x in range(6):
                    print('##' * 50)
                    print('TESTING:', x)
                    print('##' * 50)
                    e1 = [e[..., x] for e in exp0]
                    e2 = [e[..., x] for e in exp1]
                    for i1, i2 in zip(e1, e2):
                        print(':::' * 50)
                        print(':::' * 50)
                        print(i1.shape)
                        i1[np.abs(i1) < 1e-06] = 0
                        i2[np.abs(i2) < 1e-06] = 0
                        d1 = (-1 if x == 5 else 1) * i1.copy()
                        d1[np.abs(d1) < 1e-06] = 1
                        r = i2 / d1
                        r[np.logical_and(np.abs(i1) < 1e-06, np.abs(i2) < 1e-06)] = 1
                        w = np.where(np.abs(r - 1) > 0.01)
                        r[np.logical_and(np.abs(i1) < 1e-06, np.abs(i2) < 1e-06)] = 0
                        print(w)
                        if len(w[0]) > 0:
                            print('=' * 50)
                            print(np.round(i1, 6)[2])
                            print('.' * 10)
                            print(np.round(i2, 6)[2])
                            print('=' * 50)
                            print(np.round(r, 6))
        if test_inverse:
            with np.printoptions(linewidth=100000000.0, suppress=True):
                exp0 = nh3.get_internals_by_cartesians(3, strip_embedding=True)
                inv0 = nh3.get_cartesians_by_internals(3, strip_embedding=True, reembed=True, method='fast')
                for t in nput.tensor_reexpand(inv0, exp0):
                    print(np.round(t, 5))
        if test_vpt:
            nh3.setup_VPT(states=1, logger=False, cartesian_analytic_deriv_order=-1, cartesian_by_internal_derivative_method='fast')[0].print_tables(print_intensities=False, print_energy_corrections=False)
            nh3.setup_VPT(states=1, logger=False, cartesian_analytic_deriv_order=0, cartesian_by_internal_derivative_method='old')[0].print_tables(print_intensities=False, print_energy_corrections=False)
            nh3 = Molecule.from_file(nh3.source_file)
            nh3.setup_VPT(states=1, logger=False)[0].print_tables(print_intensities=False, print_energy_corrections=False)
