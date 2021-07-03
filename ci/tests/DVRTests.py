
from Peeves.TestUtils import *
from unittest import TestCase

from McUtils.Data import UnitsData, PotentialData
from McUtils.Zachary import Interpolator
import McUtils.Plots as plt

from Psience.DVR import *
from Psience.Molecools import Molecule
import numpy as np

class DVRTests(TestCase):

    #region setup
    def ho(self, grid, k=1):
        return k/2*np.power(grid, 2)
    def ho_2D(self, grid, k1=1, k2=1):
        return k1/2*np.power(grid[:, 0], 2) + k2/2*np.power(grid[:, 1], 2)
    def ho_3D(self, grid, k1=1, k2=1, k3=1):
        return k1/2*np.power(grid[:, 0], 2) + k2/2*np.power(grid[:, 1], 2) + k3/2*np.power(grid[:, 2], 2)

    def cos2D(self, grid):
        return np.cos(grid[..., 0]) * np.cos(grid[..., 1])
    def cos3D(self, grid):
        return np.cos(grid[..., 0]) * np.cos(grid[..., 1]) * np.cos(grid[..., 2])

    def cos_sin_pot(self, grid):
        return UnitsData.convert("Wavenumbers", "Hartrees")* 2500 / 8 * ((2 + np.cos(grid[..., :, 0])) * (2 + np.sin(grid[..., :, 1])) - 1)
    #endregion

    @validationTest
    def test_1D(self):
        dvr_1D = CartesianDVR(domain=(-5, 5), divs=250)
        pot = dvr_1D.run(potential_function=self.ho, result='potential_energy')
        self.assertIsInstance(pot.potential_energy, np.ndarray)

    @validationTest
    def test_energies_1D(self):
        dvr_1D = CartesianDVR()
        res = dvr_1D.run(potential_function=self.ho, domain=(-5, 5), divs=250, mass=1)
        # print(e[:5], file=sys.stderr)
        self.assertIsInstance(res.wavefunctions.energies, np.ndarray)
        self.assertTrue(np.allclose(res.wavefunctions.energies[:5].tolist(), [1/2, 3/2, 5/2, 7/2, 9/2]))

    @validationTest
    def test_energies_2D(self):
        dvr_2D = CartesianNDDVR(((-5, 5, 25), (-5, 5, 25)))
        res = dvr_2D.run(potential_function=self.ho_2D, mass=1)
        self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)

    @validationTest
    def test_energies_3D(self):
        dvr_3D = CartesianNDDVR(((-5, 5, 25), (-5, 5, 25), (-5, 5, 25)))
        res = dvr_3D.run(potential_function=self.ho_3D, mass=1)
        # print(res[0][:5], file=sys.stderr)
        self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)

    @validationTest
    def test_RingDVR1D(self):
        dvr_1D = RingDVR()
        npts = 5
        n = (5 - 1) // 2
        res = dvr_1D.run(potential_function=np.sin,
                         domain=(0, 2 * np.pi),
                         divs=npts,
                         mass=1,
                         result='grid'
                         )
        self.assertTrue(np.allclose(
            res.grid,
            (2 * np.pi) * np.arange(1, npts + 1) / npts
        ))

        res = dvr_1D.run(potential_function=np.sin,
                         domain=(0, 2 * np.pi),
                         divs=5,
                         mass=1,
                         result='kinetic_energy'
                         )
        self.assertTrue(np.allclose(res.kinetic_energy,
                                    [
                                        [1.0, -0.5854101966249685, 8.541019662496847e-2, 8.54101966249685e-2,
                                         -0.5854101966249681],
                                        [-0.5854101966249685, 1.0, -0.5854101966249685, 8.541019662496847e-2,
                                         8.54101966249685e-2],
                                        [8.541019662496847e-2, -0.5854101966249685, 1.0, -0.5854101966249685,
                                         8.541019662496847e-2],
                                        [8.54101966249685e-2, 8.541019662496847e-2, -0.5854101966249685, 1.0,
                                         -0.5854101966249685],
                                        [-0.5854101966249681, 8.54101966249685e-2, 8.541019662496847e-2,
                                         -0.5854101966249685, 1.0]
                                    ]
                                    ))

        res = dvr_1D.run(potential_function=np.sin,
                         domain=(0, 2 * np.pi),
                         divs=5,
                         mass=1,
                         result='potential_energy'
                         )
        self.assertTrue(np.allclose(np.diag(res.potential_energy), np.sin(res.grid)))

    @validationTest
    def test_RingDVR1DExplicitMass(self):

        dvr_1D = RingDVR()
        res = dvr_1D.run(potential_function=np.sin,
                         domain=(0, 2 * np.pi),
                         mass=1/(2*0.000197),
                         divs=251,
                         flavor='[0,2pi]'
                         )

        print(
            UnitsData.convert("Hartrees", "Wavenumbers")*(
                    res.wavefunctions[:5].energies[1:] -
                    res.wavefunctions[:5].energies[0]
            )
            )
        # self.assertTrue(np.allclose(
        #     res.wavefunctions[:5].energies.tolist(), [-0.536281, 0.341958, 0.854909, 2.05781, 2.08047],
        #     atol=.03  # different eigensolvers?
        # ))
        # self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)

    @validationTest
    def test_RingDVR2DExplicitMass(self):

        dvr_2D = RingNDDVR((25, 25))
        res = dvr_2D.run(potential_function=self.cos_sin_pot,
                         domain=((0, 2 * np.pi),)*2,
                         mass=[1/(2*0.000197), 1/(2*0.000197)],
                         divs=(25, 25),
                         flavor='[0,2pi]',
                         diag_mode='dense'
                         )

        print(
            UnitsData.convert("Hartrees", "Wavenumbers") * (
                    res.wavefunctions[:5].energies[1:] -
                    res.wavefunctions[:5].energies[0]
            )
        )
        # self.assertTrue(np.allclose(
        #     res.wavefunctions[:5].energies.tolist(), [-0.536281, 0.341958, 0.854909, 2.05781, 2.08047],
        #     atol=.03  # different eigensolvers?
        # ))
        self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)

    @validationTest
    def test_RingDVR1DCosMass(self):
        dvr_1D = RingDVR()
        res = dvr_1D.run(potential_function=np.sin,
                         g=np.cos,
                         g_deriv=lambda g:-np.cos(g),
                         domain=(0, 2*np.pi),
                         divs=251
                         )
        self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)

    @validationTest
    def test_Ring3D(self):
        dvr_3D = RingNDDVR((15,) * 3)
        res = dvr_3D.run(mass=1, potential_function=self.cos3D)
        # print(res[0][:5], file=sys.stderr)
        self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)

    @validationTest
    def test_Ring3DCosMass3D(self):
        dvr_3D = RingNDDVR((15,) * 3)

        res_basic = dvr_3D.run(potential_function=self.cos3D,
                               mass=1,
                               domain=((0, 2 * np.pi),) * 3,
                               divs=(15,) * 3,
                               flavor='[0,2pi]'
                               )

        g_el = lambda vals: np.full(len(vals), 1/2)
        gd_el = lambda vals: np.zeros(len(vals))

        res = dvr_3D.run(potential_function=self.cos3D,
                         domain=((0, 2 * np.pi),) * 3,
                         divs=(15,) * 3,
                         g=[
                             [g_el, 0, 0],
                             [0, g_el, 0],
                             [0, 0, g_el]
                         ],
                         g_deriv=[gd_el, gd_el, gd_el],
                         flavor='[0,2pi]',
                         num_wavefunctions=2
                         )

        self.assertTrue(np.allclose(res.wavefunctions.energies, res_basic.wavefunctions.energies))

        g_el = lambda vals: (2 + np.cos(vals[..., 0])) # plausible in size????
        gd_el = lambda vals: -np.cos(vals)/10
        res = dvr_3D.run(potential_function=self.cos3D,
                         domain=((0, 2*np.pi),) * 3,
                         divs=(15,) * 3,
                         g=[
                             [g_el, 0, 0],
                             [0, g_el, 0],
                             [0, 0, g_el]
                         ],
                         g_deriv=[gd_el, gd_el, gd_el],
                         flavor='[0,2pi]',
                         num_wavefunctions=2
                         )
        # print(res[0][:5], file=sys.stderr)
        self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)

    @validationTest
    def test_Ring2DDifferentMass(self):

        dvr_2D = RingNDDVR((15, 15))

        g_tt = lambda vals: (2 + np.cos(vals[..., 0])); gd_tt = lambda vals: -np.cos(vals[..., 0])
        g_HH = lambda vals: (2 + np.cos(2*vals[..., 1])); gd_HH = lambda vals: -np.sin(vals[..., 1])

        res = dvr_2D.run(potential_function=self.cos2D,
                         g=[
                             [g_tt, 0],
                             [0, g_HH]
                         ],
                         g_deriv=[gd_tt, gd_HH],
                         flavor='[0,2pi]'
                         )
        # print(res[0][:5], file=sys.stderr)
        self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)

        dvr_3D = RingNDDVR((15, 15, 15))
        g_tH = lambda vals: (2 + np.cos(vals[..., 0])*np.cos(2*vals[..., 1]))
        res = dvr_3D.run(potential_function=self.cos3D,
                         g=[
                             [g_tt, g_tH,    0],
                             [g_tH, g_HH,    0],
                             [0,       0, g_HH]
                         ],
                         g_deriv=[gd_tt, gd_HH, gd_HH],
                         flavor='[0,2pi]',
                         diag_mode='dense'
                         )
        # print(res[0][:5], file=sys.stderr)
        self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)

        g_tH = lambda vals: (2 + np.cos(vals[..., 0])*np.sin(vals[..., 1]))

        res = dvr_2D.run(potential_function=self.cos2D,
                         domain=((0, 2 * np.pi),) * 2,
                         divs=(15,) * 2,
                         g=[
                             [g_tt, g_tH],
                             [g_tH, g_HH]
                         ],
                         g_deriv=[gd_tt, gd_HH],
                         flavor='[0,2pi]'
                         )
        # print(res[0][:5], file=sys.stderr)
        self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)

    @debugTest
    def test_MoleculeDVR(self):

        scan_coords = Molecule.from_file(TestManager.test_data("water_HOH_scan.log"))
        scan_coords.zmatrix = [
            [0, -1, -1, -1],
            [1,  0, -1, -1],
            [2,  0,  1, -1]
        ]
        g = scan_coords.g_matrix

        a_vals = scan_coords.bond_angle(1, 0, 2)

        g_func = Interpolator(a_vals, g[:, 2, 2])
        g_deriv = g_func.derivative(2)

        # Get potential
        pot = scan_coords.potential_surface.load(coordinates=((1, 0, 2),))
        min_pos = np.argmin(pot.base.interp_data[1])
        g_eq = g[min_pos][2, 2]

        carts = CartesianDVR(domain=(np.min(a_vals), np.max(a_vals)), divs=251,
                             mass=1/g_eq,
                             potential_function=pot,
                             nodeless_ground_state=True
                             )
        res_const = carts.run()

        carts = CartesianDVR(domain=(np.min(a_vals), np.max(a_vals)),
                             divs=251,
                             g=g_func,
                             g_deriv=g_deriv,
                             potential_function=pot,
                             nodeless_ground_state=True
                             )
        res = carts.run()

        print(
            (res.wavefunctions.energies[1] - res.wavefunctions.energies[0])*UnitsData.convert("Hartrees", "Wavenumbers"),
            (res_const.wavefunctions.energies[1] - res_const.wavefunctions.energies[0])*UnitsData.convert("Hartrees", "Wavenumbers")
        )

        grid = plt.GraphicsGrid(nrows=1, ncols=3,
                                subimage_size=(400, 400),
                                spacings=[70, 0],
                                padding=[[50, 0], [50, 50]],
                                figure_label='Water HOH DVR'
                                )

        res.plot_potential(figure=grid[0, 0], zero_shift=True, plot_units='wavenumbers'); grid[0, 0].plot_label = 'HOH Potential'
        res.wavefunctions[(0, 3, 7),].plot(figure=grid[0, 1]); grid[0, 1].plot_label = 'HOH Wavefunctions'
        res_const.wavefunctions[(0, 3, 7),].plot(figure=grid[0, 2]); grid[0, 2].plot_label = 'Constant G'
        # wf_ploot = res.wavefunctions[(0, 3, 7),]
        # wf_ploot.plot(figure=grid[0, 2]); grid[0, 2].plot_label ='HOH G-Matrix'

        grid.show()





