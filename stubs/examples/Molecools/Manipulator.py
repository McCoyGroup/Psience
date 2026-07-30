"""Extracted from MolecoolsTests.test_Manipulator via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest MolecoolsTests.test_Manipulator"""

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
    def test_Manipulator(self):
        from Psience.Molecools import Molecule
        import McUtils.Numputils as nput
        import numpy as np
        bits = Molecule.from_string('C1CCOCC1=O', 'smi', confgen_opts=dict(random_seed=21232)).get_embedded_molecule()
        vv = [0, 0, 1]
        uv = [1, 0, 0]

        def prep_up_vector(view_angle):
            if view_angle is None or isinstance(view_angle, str):
                view_angle = 0
            return nput.rotation_matrix(vv, view_angle) @ np.array(uv)

        def plot_mol(view_angle, view_dist, view_x=0, view_z=0, **etc):
            uuuh = bits.plot(backend='matplotlib', image_size=[500, 500], highlight_atoms=[0, 1, 2], draw_coords={(0, 4): {'label': 'r', 'line_color': 'pink', 'label_style': {'font_size': 20, 'color': 'red'}}, (0, 1, 2): {'label': 'a', 'line_color': 'pink', 'label_style': {'font_size': 20, 'color': 'blue'}}}, plot_range=[[-3, 3], [-3, 3], [-3, 3]], view_settings={'view_distance': 5 if view_dist is None or isinstance(view_dist, str) else view_dist, 'up_vector': prep_up_vector(view_angle), 'view_vector': vv}, principle_axes=True, include_save_buttons=True)
            return uuuh.to_widget()
        import McUtils.Jupyter as interactive
        woof = interactive.Manipulator(plot_mol, ['view_angle', {'range': [-np.pi, np.pi], 'value': 0, 'continuous_update': True}], ['view_dist', {'range': [1, 20], 'value': 10}])
        woof
