"""Extracted from MolecoolsTests.test_RedundantConversion via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest MolecoolsTests.test_RedundantConversion"""

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
    def test_RedundantConversion(self):
        gggg = Molecule(['C', 'C', 'H', 'H', 'H', 'C', 'C', 'H', 'H', 'H', 'C', 'C', 'H', 'H', 'H', 'C', 'H', 'H', 'H'], [[0.0, 0.0, 0.0], [2.5221177, 0.0, 0.0], [-1.074814, 1.74867225, 0.0], [-1.06173943, -1.75324789, -0.02644452], [3.50660496, -1.80459232, -0.09993741], [4.12940149, 2.25082363, 0.07118096], [3.53603628, 4.43884306, 1.1781188], [5.96989258, 2.06662578, -0.83065536], [4.82076346, 6.03519987, 1.12058693], [1.77236274, 4.68977978, 2.19737026], [0.05431618, 1.81919915, 6.66571583], [0.10361848, 4.11457325, 7.68790077], [1.66160335, 1.19828653, 5.53987418], [-1.46076767, 4.82280651, 8.81630128], [1.71219257, 5.36239542, 7.43802933], [-2.06783807, -0.02189412, 6.91272006], [-2.80613635, -0.54136474, 5.04934629], [-1.42168595, -1.77751975, 7.80094268], [-3.61859877, 0.74452369, 8.04209746]], internals={'primitives': [(0, 1), (0, 2), (0, 3), (15, 10, 11, 13)], 'untransformed_coordinates': [(15, 10, 11, 13)]})
        print(gggg.internal_coordinates.converter_options['redundant_transformation'])
