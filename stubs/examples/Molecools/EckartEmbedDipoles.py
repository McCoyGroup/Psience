"""Extracted from MolecoolsTests.test_EckartEmbedDipoles via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest MolecoolsTests.test_EckartEmbedDipoles"""

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

    @inactiveTest
    def test_EckartEmbedDipoles(self):
        scan_file = TestManager.test_data('tbhp_030.log')
        ref_file = TestManager.test_data('tbhp_180.fchk')
        scan = Molecule.from_file(scan_file)
        ref = Molecule.from_file(ref_file)
        sel = np.where(ref.masses > 3)[0]
        pax_rot = ref.principle_axis_frame(sel=sel, inverse=True)
        rot_ref = pax_rot.apply(ref)
        transf = scan.eckart_frame(rot_ref, sel=sel)
        carts, dips = DipoleSurface.get_log_values(scan_file, keys=('StandardCartesianCoordinates', 'OptimizedDipoleMoments'))
        rot_dips = np.array([np.dot(t.transformation_function.transform, d) for t, d in zip(transf, dips)])
        self.assertTrue(np.allclose(np.linalg.norm(dips, axis=1) - np.linalg.norm(rot_dips, axis=1), 0.0))
