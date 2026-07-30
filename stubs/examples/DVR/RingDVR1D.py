"""Extracted from DVRTests.test_RingDVR1D via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest DVRTests.test_RingDVR1D"""

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
    def test_RingDVR1D(self):
        dvr_1D = RingDVR()
        npts = 5
        n = (5 - 1) // 2
        res = dvr_1D.run(potential_function=np.sin, domain=(0, 2 * np.pi), divs=npts, mass=1, result='grid')
        self.assertTrue(np.allclose(res.grid, 2 * np.pi * np.arange(1, npts + 1) / npts))
        res = dvr_1D.run(potential_function=np.sin, domain=(0, 2 * np.pi), divs=5, mass=1, result='kinetic_energy')
        self.assertTrue(np.allclose(res.kinetic_energy, [[1.0, -0.5854101966249685, 0.08541019662496847, 0.0854101966249685, -0.5854101966249681], [-0.5854101966249685, 1.0, -0.5854101966249685, 0.08541019662496847, 0.0854101966249685], [0.08541019662496847, -0.5854101966249685, 1.0, -0.5854101966249685, 0.08541019662496847], [0.0854101966249685, 0.08541019662496847, -0.5854101966249685, 1.0, -0.5854101966249685], [-0.5854101966249681, 0.0854101966249685, 0.08541019662496847, -0.5854101966249685, 1.0]]))
        res = dvr_1D.run(potential_function=np.sin, domain=(0, 2 * np.pi), divs=5, mass=1, result='potential_energy')
        self.assertTrue(np.allclose(np.diag(res.potential_energy), np.sin(res.grid)))
