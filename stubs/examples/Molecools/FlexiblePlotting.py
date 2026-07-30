"""Extracted from MolecoolsTests.test_FlexiblePlotting via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest MolecoolsTests.test_FlexiblePlotting"""

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
    def test_FlexiblePlotting(self):
        from Psience.Molecools import Molecule
        mol = Molecule.from_string('O=[C:2]([NH:1]c1cccc(C(F)(F)F)c1)[O:3]/[N:4]=[C:5](\\C1COc2ccccc2O1)[NH2:6]', 'smi')
        fig1 = mol.plot(highlight_bonds=[(0, 1), (2, 3), (4, 5)], bond_style={(0, 1): {'color': 'blue'}}, bond_radius=5, use_default_radii=False, atom_style={i: {'color': '#FF00FF'} for i in range(5)}, image_size=[800, 500], background='blue', atom_labels={i: {} for i in range(27)}, draw_coords={(0, 2): {'label': 'r<msub>1</msub>'}, (0, 1, 2): {'scaling': 0.5, 'styles': {'color': 'red'}}, (3, 4, 5): {'scaling': 1, 'label': {'text': 'a', 'offset': [3.1, 0]}}, (5, 6, 7): {'color': None, 'label': {'text': '6'}}}, label_style={'color': 'gray', 'font_size': 7}, backend='2d')
