
try:
    from Peeves.TestUtils import *
    from Peeves import BlockProfiler
except:
    pass
from unittest import TestCase

from McUtils.Coordinerds import CartesianCoordinates3D
from Psience.Vibronic import *
from Psience.Molecools import Molecule
from Psience.Modes import *

import numpy as np


class VibronicTests(TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        np.set_printoptions(linewidth=int(1e8))

    @classmethod
    def fake_mode_data(cls, n):
        fake_f = np.random.rand(3*n, 3*n)
        fake_f = fake_f @ fake_f.T

        origin = np.random.rand(n, 3)
        fake_mol = Molecule(['H']*n, origin, masses=np.ones(n))
        fake_mol = fake_mol.get_embedded_molecule(embed_properties=False)
        _, tr_modes = fake_mol.translation_rotation_modes
        tr_modes = tr_modes[0]
        proj = np.eye(3*n) - tr_modes @ tr_modes.T
        fake_f = proj @ fake_f @ proj
        _, modes = np.linalg.eigh(fake_f)

        modes = modes[:, 6:]
        inv = modes.T
        freqs = np.ones(3*n-6)
        return NormalModes(
            CartesianCoordinates3D,
            modes,
            inverse=inv,
            freqs=freqs,
            origin=fake_mol.coords,
            masses=fake_mol.masses
        )
    @classmethod
    def shift_modes(cls, modes:NormalModes, shift):
        new_origin = modes.origin + (modes.matrix @ shift).reshape(modes.origin.shape)

        return NormalModes(
            modes.basis,
            modes.matrix,
            inverse=modes.inverse,
            freqs=modes.freqs,
            origin=new_origin,
            masses=modes.masses
        )

    @debugTest
    def test_FCFs1D(self):
        np.random.seed(123123123)
        gs = self.fake_mode_data(3)
        es = self.shift_modes(gs, [1, 0, 0])

        # self.assertTrue(np.allclose(
        #     FranckCondonModel.get_fcfs(
        #         gs,
        #         es,
        #         [
        #             [0, 0, 0],
        #             [1, 0, 0],
        #             [2, 0, 0],
        #             [0, 1, 0],
        #             [0, 0, 1]
        #         ],
        #         embed=True,
        #         mass_weight=False
        #     ),
        #     [0.7788007830714044, -0.5506953149031834, 0.27534765745159184, 0, 0]
        #     )
        # )
        #
        # gs.freqs = np.array([1, 2, 2])
        # es.freqs = np.array([1, 2, 2])
        # self.assertTrue(np.allclose(
        #     FranckCondonModel.get_fcfs(
        #         gs[[0, 1]],
        #         es[[0, 1]],
        #         [
        #             [0, 0],
        #             [1, 0],
        #             [2, 0],
        #             [0, 1]
        #         ],
        #         embed=False,
        #         mass_weight=False
        #     ),
        #     [0.7788007830714044, -0.5506953149031834, 0.27534765745159184, 0]
        # ))
        #
        # gs.freqs = np.array([2, 2, 2])
        # es = self.shift_modes(gs, [1, 0, 0])
        # gs.freqs = np.array([1, 2, 2])
        # # self.assertTrue(np.allclose(
        #
        # self.assertTrue(np.allclose(
        #     FranckCondonModel.get_fcfs(
        #         gs[[0, 1]],
        #         es[[0, 1]],
        #         [
        #             [0, 0],
        #             [1, 0],
        #             [2, 0],
        #             [0, 1]
        #         ],
        #         embed=False,
        #         mass_weight=False
        #     ),
        #     [0.6957401109084786, -0.46382674060565227, 0.3826375391742289, 0]
        # ))


        gs.freqs = np.array([1, 2, 2])
        es = self.shift_modes(gs, [0, 0, 0])
        a = np.pi/6
        c = np.cos(a); s = np.sin(a)
        rot = np.array([
            [c, -s, 0],
            [s,  c, 0],
            [0,  0, 1]
        ])
        es.matrix = es.matrix @ rot
        es.inverse = rot.T @ es.inverse
        # self.assertTrue(np.allclose(

        print(
        FranckCondonModel.get_fcfs(
            gs[[0, 1]],
            es[[0, 1]],
            [
                [0, 0],
                [1, 0],
                [2, 0],
                [0, 1],
                [0, 2],
                [1, 1],
                [3, 1]
            ],
            embed=False,
            mass_weight=False
        )
        )

    @inactiveTest
    def test_FCFsNH3(self):

        gs = Molecule.from_file('/Users/Mark/Documents/Postdoc/FCFs/nh3_s0.fchk')
        es = Molecule.from_file('/Users/Mark/Documents/Postdoc/FCFs/nh3_s1.fchk')

        print()
        print(
            FranckCondonModel.get_fcfs(
                gs.normal_modes.modes.basis.to_new_modes(),
                es.normal_modes.modes.basis.to_new_modes(),
                [
                    [0, 0, 0, 0, 0, 0],
                    [1, 0, 0, 0, 0, 0],
                    [0, 1, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 0, 1],
                ],
                embed=True,
                mass_weight=True
            )
        )
