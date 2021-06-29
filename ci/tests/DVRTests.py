
from Peeves.TestUtils import *
from unittest import TestCase
from Psience.DVR import *
import numpy as np

class DVRTests(TestCase):

    def ho(self, grid, k=1):
        return k/2*np.power(grid, 2)
    def ho_2D(self, grid, k1=1, k2=1):
        return k1/2*np.power(grid[:, 0], 2) + k2/2*np.power(grid[:, 1], 2)
    def ho_3D(self, grid, k1=1, k2=1, k3=1):
        return k1/2*np.power(grid[:, 0], 2) + k2/2*np.power(grid[:, 1], 2) + k3/2*np.power(grid[:, 2], 2)

    def cos3D(self, grid):
        return np.cos(grid[..., 0]) * np.cos(grid[..., 1]) * np.cos(grid[..., 2])

    @validationTest
    def test_1D(self):
        dvr_1D = DVR("ColbertMiller1D")
        pot = dvr_1D.run(potential_function=self.ho, result='potential_energy')
        self.assertIsInstance(pot.potential_energy, np.ndarray)

    @debugTest
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
        res = dvr_1D.run(potential_function=np.sin,
                         domain=(0, 2 * np.pi),
                         divs=15,
                         flavor='[0,2pi]'
                         )
        # print(res[0][:5], file=sys.stderr)
        self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)

    @debugTest
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


    @debugTest
    def test_Ring3D(self):
        dvr_3D = DVR("ColbertMillerND")
        res = dvr_3D.run(potential_function=self.cos3D,
                         domain=((0, 2*np.pi),) * 3,
                         divs=(15,) * 3,
                         flavor='[0,2pi]'
                         )
        # print(res[0][:5], file=sys.stderr)
        self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)

    @debugTest
    def test_Ring3DCosMass3D(self):
        dvr_3D = DVR("ColbertMillerND")

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


