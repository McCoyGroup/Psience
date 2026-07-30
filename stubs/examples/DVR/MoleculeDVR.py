"""Extracted from DVRTests.test_MoleculeDVR via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest DVRTests.test_MoleculeDVR"""

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
    def test_MoleculeDVR(self):
        scan_coords = Molecule.from_file(TestManager.test_data('water_HOH_scan.log'))
        scan_coords.zmatrix = [[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]]
        g = scan_coords.g_matrix
        a_vals = scan_coords.bond_angle(1, 0, 2)
        g_func = Interpolator(a_vals, g[:, 2, 2])
        g_deriv = g_func.derivative(2)
        pot = scan_coords.potential_surface.unshare(coordinates=((1, 0, 2),))
        min_pos = np.argmin(pot.base.interp_data[1])
        g_eq = g[min_pos][2, 2]
        carts = CartesianDVR(domain=(np.min(a_vals), np.max(a_vals)), divs=251, mass=1 / g_eq, potential_function=pot, nodeless_ground_state=True)
        res_const = carts.run()
        carts = CartesianDVR(domain=(np.min(a_vals), np.max(a_vals)), divs=251, g=g_func, g_deriv=g_deriv, potential_function=pot, nodeless_ground_state=True)
        res = carts.run()
        print((res.wavefunctions.energies[1] - res.wavefunctions.energies[0]) * UnitsData.convert('Hartrees', 'Wavenumbers'), (res_const.wavefunctions.energies[1] - res_const.wavefunctions.energies[0]) * UnitsData.convert('Hartrees', 'Wavenumbers'))
        grid = plt.GraphicsGrid(nrows=1, ncols=3, subimage_size=(400, 400), spacings=[70, 0], padding=[[50, 0], [50, 50]], figure_label='Water HOH DVR')
        res.plot_potential(figure=grid[0, 0], zero_shift=True, plot_units='wavenumbers')
        grid[0, 0].plot_label = 'HOH Potential'
        res.wavefunctions[(0, 3, 7),].plot(figure=grid[0, 1])
        grid[0, 1].plot_label = 'HOH Wavefunctions'
        res_const.wavefunctions[(0, 3, 7),].plot(figure=grid[0, 2])
        grid[0, 2].plot_label = 'Constant G'
        grid.show()
