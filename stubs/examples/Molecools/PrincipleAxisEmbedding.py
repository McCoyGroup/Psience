"""Extracted from MolecoolsTests.test_PrincipleAxisEmbedding via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest MolecoolsTests.test_PrincipleAxisEmbedding"""

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
    def test_PrincipleAxisEmbedding(self):
        ref_file = TestManager.test_data('tbhp_180.fchk')
        ref = Molecule.from_file(ref_file)
        self.assertEquals(ref.center_of_mass.tolist(), [-0.10886336323443993, -7.292720327263524e-05, -0.04764041570644441])
        ref_inerts = [[-0.9998646051394727, -1.6944914059526497e-05, -0.016455123887595957], [-0.016455007408638932, 0.004930578772682442, 0.9998524501765987], [-6.419087070136426e-05, -0.9999878444790397, 0.004930190026343585]]
        inerts = ref.inertial_axes
        test_inerts = (inerts * np.array([-1, 1, 1])).T
        self.assertTrue(np.allclose(test_inerts, ref_inerts), msg="principle axes {} and {} don't align".format(test_inerts, ref_inerts))
        pax_rot = ref.principle_axis_frame()
        self.assertTrue(np.allclose(pax_rot.transformation_function.transform, inerts.T))
        rot_ref = pax_rot.apply(ref)
        self.assertTrue(np.allclose(rot_ref.center_of_mass, [0.0, 0.0, 0.0]), msg='COM: {} was {}'.format(rot_ref.center_of_mass, ref.center_of_mass))
        test_coords = np.matmul((ref.coords - ref.center_of_mass[np.newaxis])[:, np.newaxis, :], inerts).squeeze()
        self.assertTrue(np.allclose(rot_ref.coords, test_coords))
        mathematica_coords = np.array([[0.0002094928525160645, 0.0885212868882308, -0.8400509406910139], [2.389396575506497, 1.697491740062459, -0.8428256390972853], [2.435043833038253, 2.934952064361808, 0.7950074811481486], [4.0560845074996985, 0.4921123166233054, -0.8003781737352631], [-0.004983484850171475, -1.5885626031388058, 1.2992229461755922], [-0.00017490151872158886, -0.0018815600632167903, 3.5774728125123842], [-0.004314406779968471, -1.3424852433777361, 4.810480604689872], [-0.004312429484356625, -1.7659250558813848, -3.0429810385290326], [-1.6805757842711242, -2.9559004963767235, -2.984461679814903], [1.663962078887355, -2.9669237481136603, -2.9820756778710344], [0.0004171884239172418, -0.7242576512048614, -4.816727043081511], [-2.3797319162701913, 1.7110998385574014, -0.8442221100234485], [-4.053502667206945, 0.5153958278660512, -0.8051208327551433], [-2.439171179603177, 2.871593767591361, -2.543401568931165], [-2.419963556488472, 2.947396453869957, 0.7945604672548087], [2.4576648430627377, 2.8566629998551765, -2.5425989365331256]])
        self.assertTrue(np.allclose(rot_ref.coords, -mathematica_coords[:, (2, 1, 0)]), msg='{} but mathematica {}'.format(rot_ref.coords, -mathematica_coords[:, (2, 1, 0)]))
