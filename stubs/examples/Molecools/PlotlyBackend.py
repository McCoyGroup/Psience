"""Extracted from MolecoolsTests.test_PlotlyBackend via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest MolecoolsTests.test_PlotlyBackend"""

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
    def test_PlotlyBackend(self):
        from Psience.Molecools import Molecule
        tests = [[[2.38332071, -0.49643036, 5.55111512e-17], [0.00658769136, 0.98876683, 3.33066907e-16], [-2.38990841, -0.49233647, -4.99600361e-16], [-2.08988725, -3.00745588, 0.741033629], [-0.523349285, -4.08233895, -1.09291015], [2.14118284, -3.26708995, -0.417890874], [3.66973327, 0.255947571, -1.50495176], [3.43121803, -0.192532829, 1.82150928], [-0.0608841279, 2.25686365, 1.71104297], [-0.114242408, 2.33347006, -1.6298735], [-3.20094923, -0.604268579, -1.95106407], [-3.84381689, 0.336703601, 1.30762207], [-0.878968145, -3.32719027, -3.00837486], [-0.611932834, -6.15747003, -1.09757998], [3.41377283, -4.02279921, -1.90838061], [2.60551826, -4.35184035, 1.33902027]], [[2.26805368, 0.592936157, 0.0], [0.00257086576, -1.18386057, -2.22044605e-16], [-2.27062454, 0.590924417, 1.11022302e-16], [-2.30791921, 1.88012754, -2.29124159], [-0.254755984, 3.46313899, -2.50768436], [2.15182481, 1.90602459, -2.56714591], [4.06416242, -0.386408559, 0.31366895], [1.8467326, 2.01434259, 1.49524714], [-0.0223771964, -2.45799578, 1.63765593], [-0.0837319822, -2.19336565, -1.85026755], [-3.97754998, -0.575479502, 0.212293297], [-2.0568066, 1.91328775, 1.59105223], [-0.349214054, 4.5290863, -4.28173613], [-0.130106083, 4.71540063, -0.838457033], [1.91592871, 0.349950245, -3.95763078], [3.83541515, 3.03601575, -2.92206855]], [[-2.37496538, -0.516913078, -2.22044605e-16], [-0.00151564096, 1.03283778, 4.4408921e-16], [2.37648102, -0.515924698, 0.0], [2.36112435, -2.49907484, 1.67326377], [0.267785361, -4.03453305, 1.54551068], [-2.20303996, -2.64979149, 1.86935105], [-3.92449341, 0.752133206, 0.681237472], [-2.88710909, -1.17905511, -1.90284281], [-0.00253225012, 2.11878872, -1.81438423], [0.0379870511, 2.45152654, 1.54550417], [3.91483037, 0.797233384, 0.625195636], [2.89487215, -1.07093185, -1.97610252], [0.133465406, -5.19598512, -0.168889241], [0.392752782, -5.37624976, 3.17017578], [-2.57208138, -2.04623188, 3.82601529], [-3.70060088, -4.03037183, 1.3634215]], [[-2.25547156, 0.600809059, -8.8817842e-16], [0.000222534122, -1.20179599, 1.33226763e-15], [2.25524902, 0.600986929, -4.4408921e-16], [1.75066984, 2.75425742, -1.34530337], [0.183205907, 2.48701259, -3.38252635], [-2.40251269, 1.4253382, -2.75825525], [-3.98808226, -0.309181461, 0.67301456], [-1.82228951, 2.31487316, 1.16334253], [-0.0407583757, -2.46686804, 1.6377395], [0.0428858789, -2.32777829, -1.78495144], [4.01707996, -0.390110707, -0.551681364], [2.41274739, 1.21229589, 2.03414022], [-0.164379751, 4.40744415, -4.13914932], [1.14647365, 1.44677024, -4.91357472], [-3.96825875, 2.74648401, -3.11420515], [-2.73132791, -0.393039213, -3.79281507]]]
        uuh2 = Molecule(['C', 'C', 'C', 'O', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'], tests[3])
        view_vector, right_vector, up_vector = nput.view_matrix([1, 0, 0], [0, 1, 0], output_order=['x', 'y', 'z']).T
        fig = uuh2.plot(backend='x3d', image_size=[500, 500], highlight_atoms=[0, 1, 2], draw_coords={(0, 4): {'label': 'r', 'line_color': 'pink', 'label_style': {'font_size': 22, 'color': 'red'}}, (0, 1, 2): {'label': 'O', 'line_color': 'pink', 'label_style': {'font_size': 40, 'color': 'blue'}}}, view_settings={'up_vector': up_vector, 'view_vector': view_vector, 'view_distance': 9})
        fig.show()
        return
        uuh3 = Molecule(['C', 'C', 'C', 'O', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'], [[-2.25547156, 0.600809059, -8.8817842e-16], [0.000222534122, -1.20179599, 1.33226763e-15], [2.25524902, 0.600986929, -4.4408921e-16], [1.75066984, 2.75425742, -1.34530337], [0.183205907, 2.48701259, -3.38252635], [-2.40251269, 1.4253382, -2.75825525], [-3.98808226, -0.309181461, 0.67301456], [-1.82228951, 2.31487316, 1.16334253], [-0.0407583757, -2.46686804, 1.6377395], [0.0428858789, -2.32777829, -1.78495144], [4.01707996, -0.390110707, -0.551681364], [2.41274739, 1.21229589, 2.03414022], [-0.164379751, 4.40744415, -4.13914932], [1.14647365, 1.44677024, -4.91357472], [-3.96825875, 2.74648401, -3.11420515], [-2.73132791, -0.393039213, -3.79281507]])
        ploot = uuh3.plot(backend='matplotlib3D', highlight_atoms=[0, 1, 2], draw_coords={(0, 4): {'label': 'r', 'label_style': {'color': 'blue', 'font_size': 20}}, (0, 1, 2): {'label': 'q', 'label_style': {'color': 'red', 'font_size': 32, 'font_family': 'serif'}}}, view_settings={'view_vector': [0, 0, 1]})
        ploot.write('/Users/Mark/Desktop/why.html')
