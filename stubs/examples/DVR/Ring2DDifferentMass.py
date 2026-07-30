"""Extracted from DVRTests.test_Ring2DDifferentMass via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest DVRTests.test_Ring2DDifferentMass"""

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
    def test_Ring2DDifferentMass(self):
        dvr_2D = RingNDDVR((45, 45))
        g_tt = lambda vals: np.full(len(vals), 2)
        gd_tt = lambda vals: np.zeros(len(vals))
        g_HH = lambda vals: np.full(len(vals), 3)
        gd_HH = lambda vals: np.zeros(len(vals))
        zero_pot = lambda v: np.zeros(len(v))
        res = dvr_2D.run(potential_function=zero_pot, g=[[g_tt, 0], [0, g_HH]], g_deriv=[gd_tt, gd_HH])
        self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)
        f1 = res.wavefunctions.energies
        ke1 = res.kinetic_energy
        g_tH = lambda vals: np.full(len(vals), 1)
        res = dvr_2D.run(potential_function=zero_pot, g=[[g_tt, g_tH], [g_tH, g_HH]], g_deriv=[gd_tt, gd_HH])
        self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)
        f2 = res.wavefunctions.energies
        ke2 = res.kinetic_energy
        self.assertLess(np.max(np.abs(f2 - f1)), 5)
        dvr_3D = RingNDDVR((15, 15, 15))
        g_tH = lambda vals: 2 + np.cos(vals[..., 0]) * np.cos(2 * vals[..., 1])
        res = dvr_3D.run(potential_function=self.cos3D, g=[[g_tt, g_tH, 0], [g_tH, g_HH, 0], [0, 0, g_HH]], g_deriv=[gd_tt, gd_HH, gd_HH], flavor='[0,2pi]', diag_mode='dense')
        self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)
