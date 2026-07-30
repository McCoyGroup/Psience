"""Extracted from MolecoolsTests.test_SufaceTriangulation via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest MolecoolsTests.test_SufaceTriangulation"""

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
    def test_SufaceTriangulation(self):
        from McUtils.Plots import ColorPalette
        propylbenzene = Molecule.from_file(TestManager.test_data('proplybenz.hess'))
        mesh = propylbenzene.get_surface_mesh()
        mol_plot = propylbenzene.plot(backend='x3d', include_save_buttons=True)
        mesh.plot(function=lambda pts: pts[:, 2], transparency=0.2, figure=mol_plot)
        mol_plot.show()
        return
        return
        pts = surf.sampling_points
        dm = nput.distance_matrix(pts)
        np.fill_diagonal(dm, 100)
        print(np.min(dm))
        print(np.max(np.min(dm, axis=1)))
        pts2 = SphereUnionSurface.adjust_point_cloud_density(pts, centers=surf.centers, radii=surf.radii, min_component=0.6, max_iterations=250)
        dm = nput.distance_matrix(pts2)
        np.fill_diagonal(dm, 100)
        print(np.min(dm))
        print(np.max(np.min(dm, axis=1)))
        pts3 = SphereUnionSurface.point_cloud_repulsion(pts, surf.centers, surf.radii, max_iterations=25)
        dm = nput.distance_matrix(pts3)
        np.fill_diagonal(dm, 100)
        print(np.min(dm))
        print(np.max(np.min(dm, axis=1)))
        mol_plot = propylbenzene.plot(backend='x3d', image_size=[950, 700], include_save_buttons=True)
        plt.Sphere(pts * UnitsData.convert('BohrRadius', 'Angstroms'), 0.1, color='teal').plot(mol_plot)
        mol_plot.show()
        mol_plot = propylbenzene.plot(backend='x3d', image_size=[950, 700], include_save_buttons=True)
        plt.Sphere(pts2 * UnitsData.convert('BohrRadius', 'Angstroms'), 0.1, color='purple').plot(mol_plot)
        mol_plot.show()
        mol_plot = propylbenzene.plot(backend='x3d', image_size=[950, 700], include_save_buttons=True)
        plt.Sphere(pts3 * UnitsData.convert('BohrRadius', 'Angstroms'), 0.1, color='red').plot(mol_plot)
        mol_plot.show()
