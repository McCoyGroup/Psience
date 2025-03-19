
from Peeves.TestUtils import *
from Psience.Data import *
from McUtils.Coordinerds import cartesian_to_zmatrix
from McUtils.Plots import *
import McUtils.Numputils as nput
from unittest import TestCase
import sys, h5py, math, numpy as np

class DataTests(TestCase):

    maxDiff = None
    @validationTest
    def test_FChkFileDipoleSurface(self):
        fchk = TestManager.test_data("HOD_freq.fchk")
        surf = DipoleSurface.from_fchk_file(fchk)
        surf_center = surf.surfs[0].base.data['center']
        self.assertIsInstance(surf_center, np.ndarray)
        self.assertTrue(
            np.allclose(surf(surf_center) - np.array([s.base.data['ref'] for s in surf.surfs]), 0.)
        )
        self.assertEquals(surf([[0, 0, 0], [1, 0, 0], [0, 1, 0]]).shape, (1, 3))
        self.assertEquals(surf([
            [[0, 0, 0], [1, 0, 0], [0, 1, 0]],
            [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
        ]).shape, (2, 3))

    @validationTest
    def test_LogFileDipoleSurface(self):
        log = TestManager.test_data("water_OH_scan.log")
        conv = lambda x: cartesian_to_zmatrix(
            x, ordering=[[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]]
        ).coords[:, 0, 0]
        surf = DipoleSurface.from_log_file(log, conv)
        dips = surf(np.arange(.5, 2, .1))
        self.assertEquals(dips.shape, ((2-.5)/.1, 3))

        # surf_center = surf.surfs[0].base.data['center']
        # self.assertIsInstance(surf_center, np.ndarray)
        # self.assertTrue(
        #     np.allclose(surf(surf_center) - np.array([s.base.data['ref'] for s in surf.surfs]), 0.)
        # )
        # self.assertEquals(surf([[0, 0, 0], [1, 0, 0], [0, 1, 0]]).shape, (1, 3))
        # self.assertEquals(surf([
        #     [[0, 0, 0], [1, 0, 0], [0, 1, 0]],
        #     [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
        # ]).shape, (2, 3))

    @validationTest
    def test_LogFilePotentialSurface(self):
        log = TestManager.test_data("water_OH_scan.log")
        conv = lambda x: np.linalg.norm(x[:, 0] - x[:, 1], axis=1)
        surf = PotentialSurface.from_log_file(log, conv)
        pots = surf(np.arange(.5, 2, .1))
        self.assertEquals(pots.shape, ((2-.5)/.1,))


    @debugTest
    def test_PotentialRegistry(self):

        from Psience.Molecools import Molecule

        h2co_mod = PotentialRegistryAPI().get_potential('H2COPot')
        def cart_pot(coords, order=None):
            return h2co_mod.Potential.get_pot(coords)
        def internal_pot(coords, order=None):
            # internals = np.moveaxis(
            #     np.array([rOC, rCH1, rCH2, aOCH1, aOCH2, dOCHH]),
            #     0, -1
            # )
            coords = coords[..., (0, 1, 3, 2, 4, 5)]
            vals = h2co_mod.InternalsPotential.get_pot(coords, threading_mode='serial')
            # if coords.ndim > 1:
            #     vv = vals.reshape(-1)
            #     cc = coords.reshape((-1, 6))
            # else:
            #     vv = [vals]
            #     cc = [coords]
            # for c, v in zip(cc, vv):
            #     print(c, "==>", v)
            return vals

        ochh_base = Molecule.from_file(
            TestManager.test_data('OCHH_freq.fchk'),
            energy_evaluator={
                'potential_function':cart_pot,
                'permutation':[2, 3, 1, 0],
                "distance_units": "Angstroms",
                "energy_units": "Wavenumbers"
            }
        )

        cart_eng = ochh_base.calculate_energy()

        #
        # opt_ochh = ochh_base.optimize(max_displacement=.1, max_iterations=50)
        # opt_plot = opt_ochh.plot(backend='x3d')
        # base_plot = ochh_base.plot(backend='x3d', highlight_atoms=[0, 1, 2, 3], figure=opt_plot)
        # base_plot.show()
        # return
        # print(ochh_base.calculate_energy(), opt_ochh.calculate_energy())
        # print(opt_ochh.coords - ochh_base.coords)

        ochh_base = Molecule.from_file(
            TestManager.test_data('OCHH_freq.fchk'),
            energy_evaluator={
                'potential_function': internal_pot,
                # 'permutation': [2, 3, 1, 0],
                "distance_units": "Angstroms",
                "energy_units": "Wavenumbers",
                "strip_embedding": True,
                # "flatten_internals": True
            },
            internals=[
                [0, -1, -1, -1],
                [1,  0, -1, -1],
                [2,  1,  0, -1],
                [3,  1,  0,  2],
            ]
        )
        # raise Exception(ochh_base.internal_coordinates)

        opt_ochh = ochh_base.optimize(
            # method='quasi-newton'
            method='conjugate-gradient'
            # method='gradient-descent'
            , max_iterations=100
            , stencil=3
            # , logger=True
            # , max_displacement=.01
            , prevent_oscillations=3
            , restart_interval=15
            # , mesh_spacing=1e-2
        )
        opt_ochh = opt_ochh.modify(
            coords=nput.eckart_embedding(
                ochh_base.coords, opt_ochh.coords, masses=opt_ochh.masses
            ).coordinates
        )
        print("...")
        b1 = ochh_base.calculate_energy()
        b2 = opt_ochh.calculate_energy()
        print(b1, b2)
        print(opt_ochh.coords - ochh_base.coords)
        print(opt_ochh.internal_coordinates - ochh_base.internal_coordinates)
        # print(opt_ochh.internal_coordinates - ochh_base.internal_coordinates)

        # opt_ochh = ochh_base.optimize(max_displacement=.1, max_iterations=50)
        # opt_plot = opt_ochh.plot(backend='x3d')
        # base_plot = ochh_base.plot(backend='x3d', highlight_atoms=[0, 1, 2, 3], figure=opt_plot)
        # base_plot.show()