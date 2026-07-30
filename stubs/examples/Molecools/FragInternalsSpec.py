"""Extracted from MolecoolsTests.test_FragInternalsSpec via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest MolecoolsTests.test_FragInternalsSpec"""

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
    def test_FragInternalsSpec(self):
        import McUtils.Coordinerds as coordops
        ints3 = [[2.59092749, 0.0, 0.0], [2.89495928, 1.93458804, 0.0], [2.89963503, 1.8239319, 3.37828743], [2.83053209, 1.83815421, 5.06413309], [2.79207623, 2.00065824, 3.29813193], [2.52079324, 2.15594405, 4.10504516], [2.82036909, 2.05789597, 0.9634783], [2.82837755, 1.91974668, 1.11715569], [2.8622616, 1.93798376, 5.26578887], [2.09534423, 1.97234113, 4.22162045], [2.0556097, 1.98211883, 5.54531358], [2.12037191, 1.97468473, 4.16514767], [2.10238981, 1.86667929, 0.80257191], [2.07437367, 1.97295648, 4.17603043], [2.12921129, 1.8901408, 5.36779409], [2.04513453, 2.1175231, 3.14158509], [2.05011789, 2.05580992, 3.1416388], [2.06074897, 2.453391, 0.31696486], [2.09040482, 1.89999553, 2.33908385], [2.11914653, 1.84889244, 1.97498365], [2.13523637, 1.74785729, 5.16708695], [2.08252609, 1.9355175, 4.49096479], [2.07995171, 1.95503155, 5.23240844], [2.09065558, 1.93768086, 2.16018763]]
        zm3 = [[0, -1, -2, -3], [1, 0, -1, -2], [2, 1, 0, -1], [3, 2, 1, 0], [4, 3, 2, 1], [5, 4, 3, 2], [6, 5, 4, 3], [7, 5, 4, 3], [8, 4, 3, 2], [9, 8, 4, 3], [10, 1, 0, 2], [11, 2, 1, 0], [12, 2, 11, 1], [13, 3, 2, 1], [14, 3, 13, 2], [15, 4, 3, 2], [16, 6, 5, 4], [17, 6, 16, 5], [18, 7, 6, 5], [19, 7, 18, 6], [20, 7, 18, 19], [21, 8, 7, 6], [22, 8, 21, 7], [23, 9, 8, 7], [24, 9, 23, 8]]
        bits = Molecule.from_string('FC1CCC(C(=C)C)CC1.C1OCC(C(OC)O)C1C(OC)O', 'smi', confgen_opts=dict(verbose=True, random_seed=12321))
        zm = bits.get_bond_zmatrix(validate=True, connect_fragments=False)
        spec = coordops.InternalSpec.from_zmatrix(*(np.array(z).tolist() for z in zm))
        ints = spec.cartesians_to_internals(bits.coords)
        zm2, _ = spec.get_zmat_conv()
        carts, exp = spec.internals_to_cartesians(ints, reference_cartesians=bits.coords, order=1)
        uuuh = spec.get_direct_inverses(carts, order=1)
        print(uuuh[1][0].shape, exp[0][1].shape)
        print(np.round(uuuh[1] @ exp[0][1], 8)[:6, :6])
        return
