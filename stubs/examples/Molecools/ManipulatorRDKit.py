"""Extracted from MolecoolsTests.test_ManipulatorRDKit via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest MolecoolsTests.test_ManipulatorRDKit"""

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
    def test_ManipulatorRDKit(self):
        from Psience.Molecools import Molecule
        import McUtils.Numputils as nput
        import numpy as np
        bits = Molecule.from_string('C1CCOCC1=O', 'smi').get_embedded_molecule()

        def prep_up_vector(view_angle):
            if view_angle is None or isinstance(view_angle, str):
                view_angle = 0
            return nput.rotation_matrix([0, 0, 1], view_angle) @ np.array([0, 1, 0])

        def plot_mol(view_angle, view_dist, view_x=0, view_z=0, **etc):
            uuuh = bits.plot(backend='rdkit', image_size=[500, 500], highlight_atoms=[0, 1, 2], draw_coords={(0, 4): 'r', (0, 1, 2): {'label': 't', 'label_style': {'font_size': 22, 'font_family': 'Arial', 'color': 'red'}}}, postdraw=[{'pattern': 'labels_1', 'classes': True, 'replacement': {'text': 'θ', 'mode': 'svg'}}], view_settings={'view_distance': 5 if view_dist is None or isinstance(view_dist, str) else view_dist, 'up_vector': prep_up_vector(view_angle), 'view_vector': [0, 0, 1]})
            return uuuh.to_widget()
