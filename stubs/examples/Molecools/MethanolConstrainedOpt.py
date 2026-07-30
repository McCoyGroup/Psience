"""Extracted from MolecoolsTests.test_MethanolConstrainedOpt via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest MolecoolsTests.test_MethanolConstrainedOpt"""

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
    def test_MethanolConstrainedOpt(self):
        import McUtils.Coordinerds as coordops
        methanol = Molecule(['C', 'O', 'H', 'H', 'H', 'H'], [[-0.6988896, 0.00487717, 0.00528378], [1.69694605, -1.08628154, -0.22505594], [-1.27384807, -0.22559494, 1.98568702], [-0.59371792, 2.01534286, -0.59617633], [-2.04665278, -0.99128091, -1.21504054], [2.91616232, 0.28293736, 0.04530201]])
        eng0 = methanol.calculate_energy()
        methanol_zmatrix = coordops.reindex_zmatrix(coordops.functionalized_zmatrix(3, {(2, 1, 0): [[0, -1, -2, -3], [1, -1, 0, -2], [2, -1, 0, 1]]}), [5, 1, 0, 2, 3, 4])
        meth_int = methanol.modify(internals=methanol_zmatrix)
        ints0 = meth_int.internal_coordinates
        meth_int = meth_int.optimize(max_displacement=0.5, max_iterations=500, coordinate_constraints=[(0, 1)])
        eng1 = meth_int.calculate_energy()
        ints1 = meth_int.internal_coordinates
