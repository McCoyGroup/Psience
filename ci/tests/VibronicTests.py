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

    @validationTest
    def test_FCFsAnalytic(self):
        # print()

        np.random.seed(123123123)
        gs = self.fake_mode_data(3)
        es = self.shift_modes(gs, [1, 0, 0])

        self.assertTrue(np.allclose(
            FranckCondonModel.get_fcfs(
                gs,
                es,
                [
                    [0, 0, 0],
                    [1, 0, 0],
                    [2, 0, 0],
                    [0, 1, 0],
                    [0, 0, 1]
                ],
                embed=True,
                mass_weight=False
            ),
            [0.7788007830714044, 0.5506953149031834, 0.27534765745159184, 0, 0]
        ))

        gs.freqs = np.array([1, 2, 2])
        es.freqs = np.array([1, 2, 2])
        self.assertTrue(np.allclose(
            FranckCondonModel.get_fcfs(
                gs[[0, 1]],
                es[[0, 1]],
                [
                    [0, 0],
                    [1, 0],
                    [2, 0],
                    [0, 1]
                ],
                embed=False,
                mass_weight=False
            ),
            [0.7788007830714044, 0.5506953149031834, 0.27534765745159184, 0]
        ))

        gs.freqs = np.array([2, 2, 2])
        es = self.shift_modes(gs, [1, 0, 0])
        gs.freqs = np.array([1, 2, 2])
        self.assertTrue(np.allclose(
            FranckCondonModel.get_fcfs(
                gs[[0, 1]],
                es[[0, 1]],
                [
                    [0, 0],
                    [1, 0],
                    [2, 0],
                    [0, 1]
                ],
                embed=False,
                mass_weight=False
            ),
            [0.6957401109084786, 0.46382674060565227, 0.3826375391742289, 0]
        ))

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
        es = self.shift_modes(es, [1, 0, 0])
        self.assertTrue(np.allclose(
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
                    [3, 1],
                    [4, 0]
                ],
                embed=False,
                mass_weight=False
            ),
            [0.7496767974532862, 0.5782925969208352, 0.26724127584978014, 0.07869565469361299, 0.05403238910624021, -0.0505874828034262, -0.06072786570450098, 0.008309307609240951]
        ))

        gs.freqs = np.array([1, 1, 3])
        es = self.shift_modes(gs, [0, 0, 0])
        es.freqs = np.array([2, 5, 7])
        a = np.pi / 6
        c = np.cos(a); s = np.sin(a)
        rot = (
                np.eye(3)
                @ np.array([
                    [c, -s, 0],
                    [s,  c, 0],
                    [0,  0, 1]
                ])
                @ np.array([
                    [c,  0, -s],
                    [0,  1,  0],
                    [s,  0,  c]
                ])
                # @ np.array([
                #     [1,  0,  0],
                #     [0,  c, -s],
                #     [0,  s,  c]
                # ])
        )
        # rot = rot.T
        es.matrix = es.matrix @ rot
        es.inverse = rot.T @ es.inverse
        es = self.shift_modes(es, [.5, 0, 0])
        self.assertTrue(np.allclose(
            FranckCondonModel.get_fcfs(
                gs,
                es,
                [
                    # [0, 0, 0],
                    # [1, 0, 0],
                    [2, 0, 0],
                    [1, 1, 0],
                    [1, 0, 1],
                    # [1, 2, 1]
                ],
                embed=False,
                mass_weight=False
            ),
            [0.16796089572281053, -1.0693197339113112e-16, -0.10980392299108514]
        ))

    @validationTest
    def test_FCFsNH3(self):
        # print()

        fc_model = FranckCondonModel.from_files(
            TestManager.test_data('nh3_s0.fchk'),
            TestManager.test_data('nh3_s1.fchk'),
            logger=True
        )

        # od = fc_model.get_overlap_data()
        #
        # print(fc_model.gs_nms.freqs)
        # print(fc_model.es_nms.freqs)
        # uuugh = np.power(
        #     fc_model.get_overlaps(
        #         [
        #             [0, 0, 0, 0, 0, 0],
        #             [1, 0, 0, 0, 0, 0],
        #             [3, 0, 0, 0, 0, 0],
        #             [0, 0, 1, 0, 0, 0],
        #             [0, 0, 0, 1, 0, 0],
        #             [0, 0, 0, 0, 1, 0],
        #             [0, 0, 0, 0, 0, 1]
        #         ]
        #     ),
        #     2
        # )
        # print(uuugh)
        # raise Exception(...)

        uuugh = np.power(
                fc_model.get_overlaps(
                    [
                        [0, 0, 0, 0, 0, 0],
                        [1, 0, 0, 0, 0, 0],
                        [2, 0, 0, 0, 0, 0],
                        [3, 0, 0, 0, 0, 0],
                        [0, 1, 0, 0, 0, 0],
                        [0, 0, 1, 0, 0, 0],
                        [0, 0, 0, 1, 0, 0],
                        [0, 0, 0, 0, 1, 0],
                        [0, 0, 0, 0, 0, 1]
                    ]
                ),
                2
            )
        print(uuugh)
        # print(np.sum(uuugh))
        # self.assertTrue(np.allclose(
        #     uuugh,
        #     [2.79672434e-02, 2.43447521e-02, 3.66375168e-05, 1.45628441e-02, 5.93123060e-08, 1.35127810e-08, 6.98421424e-08,
        #      3.91535362e-09, 6.62467294e-02]
        # ))
        # print(uuugh)
        # print(uuugh / np.max(uuugh))

        spec = fc_model.get_spectrum({'threshold':4000 / UnitsData.hartrees_to_wavenumbers, 'max_quanta':6})
        # spec.plot().show()
        # print(spec.intensities)

    @classmethod
    def load_log_mol(cls, log_file, ref_file):

        klow_exc = Molecule.from_file(ref_file)
        from McUtils.GaussianInterface import GaussianLogReader
        woof = log_file
        with GaussianLogReader(woof) as reader:
            nm_data = reader.parse(['InputCartesianCoordinates', 'NormalModes'])
            carts = nm_data['InputCartesianCoordinates'][1][-1]
            freqs, matrix = nm_data['NormalModes'][0]
            gvals, gvecs = np.linalg.eigh(matrix.T @ matrix)
            g12 = gvecs @ np.diag(1/np.sqrt(gvals)) @ gvecs.T
            matrix = matrix @ g12
            cart_mass = np.sqrt(klow_exc.atomic_masses)
            matrix = np.reshape(
                matrix.reshape(-1, 3, matrix.shape[-1]) / cart_mass[:, np.newaxis, np.newaxis],
                matrix.shape
            )
            # altmat = klow_exc.normal_modes.modes.basis.matrix

            klow = Molecule(
                klow_exc.atoms,
                carts * UnitsData.convert("Angstroms", "BohrRadius"),
                masses=klow_exc.masses,
                normal_modes={
                    'freqs':freqs * UnitsData.convert("Wavenumbers", "Hartrees"),
                    'matrix':matrix
                }
            )
            # klow = klow.get_embedded_molecule()
            # klow_modes = klow.normal_modes.modes.basis.to_new_modes().make_mass_weighted()
            # klow_exc = klow_exc.get_embedded_molecule()
            # new_modes = klow_exc.normal_modes.modes.basis.to_new_modes().make_mass_weighted()

        return klow, klow_exc
        # raise Exception(...)
    @debugTest
    def test_FCFsBig(self):

        root = os.path.expanduser('~/Documents/Postdoc/FCFs/')
        path = lambda *p:os.path.join(root, *p)

        s = 6
        s = f"Na{s}"
        # s = "K_low"
        # s = "K_min"

        # fc_model = FranckCondonModel.from_mols(
        #     *reversed(
        #         self.load_log_mol(
        #         path(f'{s}_m0_ex_freq.log'),
        #         path(f'{s}_m0_s0.fchk'))
        #     ),
        #     logger=True
        # )

        # gs, exc = self.load_log_mol(path(f'{s}_m0_freq.log'), path(f'{s}_m0_s1.fchk'))
        # fc_model = FranckCondonModel.from_mols(
        #     gs, exc,
        #     logger=True
        # )

        fc_model = FranckCondonModel.from_files(
            path(f'{s}_m0_s0.fchk'),
            path(f'{s}_m0_s1.fchk'),
            logger=True
        )

        main_modes = [1, 4, 7]
        aux_modes = [0, 2, 5, 8,  9, 10, 12, 13, 15, 20, 24, 25, 27, 28,
                     33, 37, 41, 44, 46, 47, 51, 55, 56, 59, 60]
        # nz_modes = [ 0,  1,  2,  4,  5,  7,  8,  9, 10, 12, 13, 15, 20, 24, 25, 27, 28,
        #              33, 37, 41, 44, 46, 47, 51, 55, 56, 59, 60 ]


        q = 3
        # t0 = 0000
        t = 2000
        # with BlockProfiler('FCFs'):
        oof, (_, excitations) = fc_model.get_spectrum(
            {
                'threshold': t / UnitsData.hartrees_to_wavenumbers,
                # 'min_quanta': q,
                # 'min_freq': t0 / UnitsData.hartrees_to_wavenumbers,
                'max_quanta': q,
                # 'max_state': [
                #     (
                #         q
                #             if s in main_modes else
                #         2
                #             if s in aux_modes else
                #         0
                #     )
                #     for s in range(len(fc_model.es_nms.freqs))
                # ]
            },
            return_states=True
        )
        oof.plot().savefig(path(f'{s}_m0_{t}_{q}_quanta.png'))
        np.savez(path(f'{s}_m0_{t}_{q}_spec.npz'),
                 freqs=oof.frequencies,
                 fcfs=oof.intensities,
                 excitations=excitations)
        # with open(path(f'Na2_m0_{t}_{q}_quanta.txt'), 'w+') as out:
        #     print('Freqs:', oof.frequencies, file=out)
        #     print('FCFs:', oof.intensities, file=out)
