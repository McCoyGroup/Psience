
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

    def test_FCFs(self):
        np.random.seed(123123123)
        gs = self.fake_mode_data(3)
        es = self.shift_modes(gs, [1, 0, 0])

        print(
            FranckCondonModel.get_fcfs(
                gs,
                es,
                [
                    [0, 0, 0],
                    [1, 0, 0],
                    [2, 0, 0],
                    # [0, 1, 0],
                    # [0, 0, 1]
                ],
                embed=False
            )
        )
