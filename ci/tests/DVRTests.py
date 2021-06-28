
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

    @validationTest
    def test_energies_1D(self):
        dvr_1D = DVR("ColbertMiller1D")
        res = dvr_1D.run(potential_function=self.ho, divs=150)
        # print(e[:5], file=sys.stderr)
        self.assertIsInstance(res.wavefunctions.energies, np.ndarray)

    @validationTest
    def test_energies_2D(self):
        dvr_2D = DVR("ColbertMillerND")
        res = dvr_2D.run(potential_function=self.ho_2D, divs=(25, 25))
        # print(res[0][:5], file=sys.stderr)
        self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)

    @debugTest
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
    def test_Ring3DCosMass(self):
        dvr_3D = DVR("ColbertMillerND")
        res = dvr_3D.run(potential_function=self.cos3D,
                         domain=((0, 2*np.pi),) * 3,
                         divs=(15,) * 3,
                         flavor='[0,2pi]'
                         )
        # print(res[0][:5], file=sys.stderr)
        self.assertIsInstance(res.wavefunctions[0].data, np.ndarray)


