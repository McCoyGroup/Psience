
from Peeves.TestUtils import *
from unittest import TestCase

from McUtils.Data import UnitsData
from Psience.DVR import *
import numpy as np

class DVRTests(TestCase):

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

    @validationTest
    def test_1D(self):
        dvr_1D = DVR("ColbertMiller1D")
        pot = dvr_1D.run(potential_function=self.ho, result='potential_energy')
        self.assertIsInstance(pot.potential_energy, np.ndarray)

    @validationTest
    def test_energies_1D(self):
        dvr_1D = DVR("ColbertMiller1D")
        res = dvr_1D.run(potential_function=self.ho, domain=(-5, 5), divs=250)
        # print(e[:5], file=sys.stderr)
        self.assertIsInstance(res.wavefunctions.energies, np.ndarray)
        self.assertTrue(np.allclose(res.wavefunctions.energies[:5].tolist(), [1/2, 3/2, 5/2, 7/2, 9/2]))

    @validationTest
    def test_energies_2D(self):
        dvr_2D = DVR("ColbertMillerND")
        res = dvr_2D.run(potential_function=self.ho_2D, domain=((-5, 5), (-5, 5)), divs=(25, 25))
        # print(res[0][:5], file=sys.stderr)
        self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)

    @validationTest
    def test_energies_3D(self):
        dvr_3D = DVR("ColbertMillerND")
        res = dvr_3D.run(potential_function=self.ho_3D, domain=((-5, 5),)*3, divs=(15,)*3)
        # print(res[0][:5], file=sys.stderr)
        self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)

    @debugTest
    def test_RingDVR1D(self):
        dvr_1D = DVR("ColbertMiller1D")
        npts = 5
        n = (5-1)//2
        res = dvr_1D.run(potential_function=np.sin,
                         domain=(0, 2 * np.pi),
                         divs=npts,
                         flavor='[0,2pi]',
                         result='grid'
                         )
        self.assertTrue(np.allclose(
            res.grid,
            (2*np.pi) * np.arange(1, npts+1)/npts
        ))

        res = dvr_1D.run(potential_function=np.sin,
                         domain=(0, 2 * np.pi),
                         divs=5,
                         flavor='[0,2pi]',
                         result='kinetic_energy'
                         )
        self.assertTrue(np.allclose(res.kinetic_energy,
                                    [
                                        [1.0, -0.5854101966249685, 8.541019662496847e-2, 8.54101966249685e-2, -0.5854101966249681],
                                        [-0.5854101966249685, 1.0, -0.5854101966249685, 8.541019662496847e-2, 8.54101966249685e-2],
                                        [8.541019662496847e-2, -0.5854101966249685, 1.0, -0.5854101966249685, 8.541019662496847e-2],
                                        [8.54101966249685e-2, 8.541019662496847e-2, -0.5854101966249685, 1.0, -0.5854101966249685],
                                        [-0.5854101966249681, 8.54101966249685e-2, 8.541019662496847e-2, -0.5854101966249685, 1.0]
                                    ]
                                    ))

        res = dvr_1D.run(potential_function=np.sin,
                         domain=(0, 2 * np.pi),
                         divs=5,
                         flavor='[0,2pi]',
                         result='potential_energy'
                         )
        self.assertTrue(np.allclose(np.diag(res.potential_energy), np.sin(res.grid)))

        res = dvr_1D.run(potential_function=np.sin,
                         domain=(0, 2 * np.pi),
                         divs=251,
                         flavor='[0,2pi]'
                         )

        self.assertTrue(np.allclose(
            res.wavefunctions[:5].energies.tolist(), [-0.536281, 0.341958, 0.854909, 2.05781, 2.08047],
            atol=.03 # different eigensolvers?
        ))
        self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)

    @debugTest
    def test_RingDVR1D(self):
        dvr_1D = DVR("ColbertMiller1D")
        npts = 5
        n = (5 - 1) // 2
        res = dvr_1D.run(potential_function=np.sin,
                         domain=(0, 2 * np.pi),
                         divs=npts,
                         flavor='[0,2pi]',
                         result='grid'
                         )
        self.assertTrue(np.allclose(
            res.grid,
            (2 * np.pi) * np.arange(1, npts + 1) / npts
        ))

        res = dvr_1D.run(potential_function=np.sin,
                         domain=(0, 2 * np.pi),
                         divs=5,
                         flavor='[0,2pi]',
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
                         flavor='[0,2pi]',
                         result='potential_energy'
                         )
        self.assertTrue(np.allclose(np.diag(res.potential_energy), np.sin(res.grid)))

    @debugTest
    def test_RingDVR1DExplicitMass(self):

        dvr_1D = DVR("ColbertMiller1D")
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

    @debugTest
    def test_RingDVR2DExplicitMass(self):

        dvr_2D = DVR("ColbertMillerND")
        res = dvr_2D.run(potential_function=self.cos2D,
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
        self.assertTrue(np.allclose(
            res.wavefunctions[:5].energies.tolist(), [-0.536281, 0.341958, 0.854909, 2.05781, 2.08047],
            atol=.03  # different eigensolvers?
        ))
        self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)

    @validationTest
    def test_RingDVR1DCosMass(self):
        dvr_1D = DVR("ColbertMiller1D")
        res = dvr_1D.run(potential_function=np.sin,
                         g=np.cos,
                         g_deriv=lambda g:-np.cos(g),
                         domain=(0, 2*np.pi),
                         divs=15,
                         flavor='[0,2pi]'
                         )
        # print(res[0][:5], file=sys.stderr)
        self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)


    @validationTest
    def test_Ring3D(self):
        dvr_3D = DVR("ColbertMillerND")
        res = dvr_3D.run(potential_function=self.cos3D,
                         domain=((0, 2*np.pi),) * 3,
                         divs=(15,) * 3,
                         flavor='[0,2pi]'
                         )
        # print(res[0][:5], file=sys.stderr)
        self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)

    @validationTest
    def test_Ring3DCosMass3D(self):
        dvr_3D = DVR("ColbertMillerND")

        res_basic = dvr_3D.run(potential_function=self.cos3D,
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

        g_el = lambda vals: (2 + np.cos(vals)) # plausible in size????
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

        dvr_2D = DVR("ColbertMillerND")

        g_tt = lambda vals: (2 + np.cos(vals)); gd_tt = lambda vals: -np.cos(vals)
        g_HH = lambda vals: (2 + np.cos(2*vals)); gd_HH = lambda vals: -np.sin(vals)

        res = dvr_2D.run(potential_function=self.cos2D,
                         domain=((0, 2 * np.pi),) * 2,
                         divs=(15,) * 2,
                         g=[
                             [g_tt, 0],
                             [0, g_HH]
                         ],
                         g_deriv=[gd_tt, gd_HH],
                         flavor='[0,2pi]'
                         )
        # print(res[0][:5], file=sys.stderr)
        self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)

        g_tH = lambda vals: (2 + np.cos(vals[..., 0])*np.cos(2*vals[..., 1]))
        res = dvr_2D.run(potential_function=self.cos3D,
                         domain=((0, 2 * np.pi),) * 3,
                         divs=(15,) * 3,
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


