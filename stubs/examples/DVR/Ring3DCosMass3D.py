"""Extracted from DVRTests.test_Ring3DCosMass3D via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest DVRTests.test_Ring3DCosMass3D"""

from Peeves.TestUtils import *
from unittest import TestCase
from McUtils.Data import UnitsData, PotentialData
from McUtils.Zachary import Interpolator
from McUtils.Extensions import ModuleLoader
import McUtils.Plots as plt
from Psience.DVR import *
from Psience.Molecools import Molecule
import numpy as np

class DVRTests(TestCase):

    def ho(self, grid, k=1):
        return k / 2 * np.power(grid, 2)

    def ho_2D(self, grid, k1=1, k2=1):
        return k1 / 2 * np.power(grid[:, 0], 2) + k2 / 2 * np.power(grid[:, 1], 2)

    def ho_3D(self, grid, k1=1, k2=1, k3=1):
        return k1 / 2 * np.power(grid[:, 0], 2) + k2 / 2 * np.power(grid[:, 1], 2) + k3 / 2 * np.power(grid[:, 2], 2)

    def cos2D(self, grid):
        return np.cos(grid[..., 0]) * np.cos(grid[..., 1])

    def cos3D(self, grid):
        return np.cos(grid[..., 0]) * np.cos(grid[..., 1]) * np.cos(grid[..., 2])

    def cos_sin_pot(self, grid):
        return UnitsData.convert('Wavenumbers', 'Hartrees') * 2500 / 8 * ((2 + np.cos(grid[..., :, 0])) * (2 + np.sin(grid[..., :, 1])) - 1)

    def setupMBPolModel(self):
        ...

    @validationTest
    def test_Ring3DCosMass3D(self):
        dvr_3D = RingNDDVR((15,) * 3)
        res_basic = dvr_3D.run(potential_function=self.cos3D, mass=1, domain=((0, 2 * np.pi),) * 3, divs=(15,) * 3, flavor='[0,2pi]')
        g_el = lambda vals: np.full(len(vals), 1 / 2)
        gd_el = lambda vals: np.zeros(len(vals))
        res = dvr_3D.run(potential_function=self.cos3D, domain=((0, 2 * np.pi),) * 3, divs=(15,) * 3, g=[[g_el, 0, 0], [0, g_el, 0], [0, 0, g_el]], g_deriv=[gd_el, gd_el, gd_el], flavor='[0,2pi]', num_wavefunctions=2)
        self.assertTrue(np.allclose(res.wavefunctions.energies, res_basic.wavefunctions.energies))
        g_el = lambda vals: 2 + np.cos(vals[..., 0])
        gd_el = lambda vals: -np.cos(vals) / 10
        res = dvr_3D.run(potential_function=self.cos3D, domain=((0, 2 * np.pi),) * 3, divs=(15,) * 3, g=[[g_el, 0, 0], [0, g_el, 0], [0, 0, g_el]], g_deriv=[gd_el, gd_el, gd_el], flavor='[0,2pi]', num_wavefunctions=2)
        self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)
