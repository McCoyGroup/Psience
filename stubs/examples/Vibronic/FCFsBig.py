"""Extracted from VibronicTests.test_FCFsBig via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VibronicTests.test_FCFsBig"""

import scipy.linalg
try:
    from Peeves.TestUtils import *
    from Peeves import BlockProfiler
except:
    pass
from unittest import TestCase
import os, sys, numpy as np
from McUtils.Coordinerds import CartesianCoordinates3D
from McUtils.Data import UnitsData
from Psience.Vibronic import *
from Psience.Molecools import Molecule
from Psience.Modes import *

class VibronicTests(TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        np.set_printoptions(linewidth=int(100000000.0))

    @classmethod
    def fake_mode_data(cls, n):
        fake_f = np.random.rand(3 * n, 3 * n)
        fake_f = fake_f @ fake_f.T
        origin = np.random.rand(n, 3)
        fake_mol = Molecule(['H'] * n, origin, masses=np.ones(n))
        fake_mol = fake_mol.get_embedded_molecule(embed_properties=False)
        _, tr_modes = fake_mol.translation_rotation_modes
        tr_modes = tr_modes[0]
        proj = np.eye(3 * n) - tr_modes @ tr_modes.T
        fake_f = proj @ fake_f @ proj
        _, modes = np.linalg.eigh(fake_f)
        modes = modes[:, 6:]
        inv = modes.T
        freqs = np.ones(3 * n - 6)
        return NormalModes(CartesianCoordinates3D, modes, inverse=inv, freqs=freqs, origin=fake_mol.coords, masses=fake_mol.masses)

    @classmethod
    def shift_modes(cls, modes: NormalModes, shift):
        new_origin = modes.origin + (modes.matrix @ shift).reshape(modes.origin.shape)
        return NormalModes(modes.basis, modes.matrix, inverse=modes.inverse, freqs=modes.freqs, origin=new_origin, masses=modes.masses)

    @classmethod
    def load_log_mol(cls, log_file, ref_file):
        klow_exc = Molecule.from_file(ref_file)
        from McUtils.GaussianInterface import GaussianLogReader
        woof = log_file
        with GaussianLogReader(woof) as reader:
            nm_data = reader.parse(['InputCartesianCoordinates', 'NormalModes'])
            carts = nm_data['InputCartesianCoordinates'][1][-1]
            freqs, masses, matrix = nm_data['NormalModes'][0]
            renorm = np.diag(1 / np.sqrt(masses))
            matrix = matrix @ renorm
            gvals, gvecs = np.linalg.eigh(matrix.T @ matrix)
            g12 = gvecs @ np.diag(1 / np.sqrt(gvals)) @ gvecs.T
            matrix = matrix @ g12
            cart_mass = np.sqrt(klow_exc.atomic_masses)
            matrix = np.reshape(matrix.reshape(-1, 3, matrix.shape[-1]) / cart_mass[:, np.newaxis, np.newaxis], matrix.shape)
            klow = Molecule(klow_exc.atoms, carts * UnitsData.convert('Angstroms', 'BohrRadius'), masses=klow_exc.masses, normal_modes={'freqs': freqs * UnitsData.convert('Wavenumbers', 'Hartrees'), 'matrix': matrix})
        return (klow, klow_exc)

    @inactiveTest
    def test_FCFsBig(self):
        root = os.path.expanduser('~/Documents/Postdoc/Projects/Tim-Ned-FCFs/FCFs')
        path = lambda *p: os.path.join(root, *p)
        for s in [2, 4, 6]:
            s = f'Na{s}'
            fc_model = FranckCondonModel.from_files(path(f'{s}_m0_s0.fchk'), path(f'{s}_m0_s1.fchk'), logger=False, mode_selection=np.arange(80, dtype=int))
            main_modes = [1, 4, 7]
            aux_modes = [0, 2, 5, 8, 9, 10, 12, 13, 15, 20, 24, 25, 27, 28, 33, 37, 41, 44, 46, 47, 51, 55, 56, 59, 60]
            q = 3
            t = 1000
            oof, (_, excitations) = fc_model.get_spectrum({'threshold': t / UnitsData.hartrees_to_wavenumbers, 'max_quanta': q}, return_states=True)
            np.savez(path(f'{s}_m0_{t}_{q}_spec.npz'), freqs=oof.frequencies, fcfs=oof.intensities, excitations=excitations)
            with open(path(f'Na2_m0_{t}_{q}_quanta.txt'), 'w+') as out:
                print('Freqs:', oof.frequencies, file=out)
                print('FCFs:', oof.intensities, file=out)
            oof.plot().savefig(path(f'{s}_m0_{t}_{q}_quanta.png'))
            with np.printoptions(linewidth=100000000.0, threshold=None):
                print(oof.intensities[:1500])
