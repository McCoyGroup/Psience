"""Extracted from DVRTests.test_RingDVR2DExplicitMass via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest DVRTests.test_RingDVR2DExplicitMass"""

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
    def test_RingDVR2DExplicitMass(self):
        dvr_2D = RingNDDVR((25, 25))
        res = dvr_2D.run(potential_function=self.cos_sin_pot, domain=((0, 2 * np.pi),) * 2, mass=[1 / (2 * 0.000197), 1 / (2 * 0.000197)], divs=(25, 25), flavor='[0,2pi]', diag_mode='dense')
        print(UnitsData.convert('Hartrees', 'Wavenumbers') * (res.wavefunctions[:5].energies[1:] - res.wavefunctions[:5].energies[0]))
        self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)
