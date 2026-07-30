"""Extracted from DVRTests.test_MBPolDVR via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest DVRTests.test_MBPolDVR"""

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

    @debugTest
    def test_MBPolDVR(self):
        loader = ModuleLoader(TestManager.current_manager().test_data_dir)
        mbpol = loader.load('LegacyMBPol').potential
        mol = Molecule.from_file(TestManager.test_data('water_freq.fchk'))
        mol.zmatrix = [[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]]
        disps_r = np.linspace(-0.2, 0.5, 25) / UnitsData.bohr_to_angstroms
        disps_grid = np.array(np.meshgrid(disps_r, disps_r))
        dr_coords = mol.get_displaced_coordinates(np.moveaxis(disps_grid, 0, -1).reshape(-1, disps_grid.shape[0]), [0, 1], use_internals='convert', strip_embedding=True)
        pot = mbpol(dr_coords)[0]
        plt.ContourPlot(*disps_grid, pot.reshape(disps_grid[0].shape) * 219475.6).show()
