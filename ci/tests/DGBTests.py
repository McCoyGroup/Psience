import collections
import gc
import itertools
import os

import scipy.linalg

import McUtils.Zachary
from Peeves.TestUtils import *
from Peeves import BlockProfiler
from unittest import TestCase

from McUtils.Data import UnitsData, PotentialData, AtomData
from McUtils.Zachary import Interpolator, FiniteDifferenceDerivative
import McUtils.Plots as plt
import McUtils.Numputils as nput
from McUtils.GaussianInterface import GaussianFChkReader, GaussianLogReader
from McUtils.Extensions import ModuleLoader

from McUtils.Scaffolding import Checkpointer

from Psience.DGB import *
from Psience.Molecools import Molecule
from Psience.AIMD import AIMDSimulator

import numpy as np

class DGBTests(TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        np.set_printoptions(linewidth=int(1e8))

    @validationTest
    def test_Harmonic(self):
        ndivs = [10]*3
        domain = [[-1, 1], [-1, 1], [-1, 1]] #, [-1, 1], [-1, 1]] # unit cube
        ndim = len(domain)

        pts = np.array(
            np.meshgrid(*(np.linspace(d[0], d[1], n) for d,n in zip(domain, ndivs)))
        ).T.reshape(-1, ndim)

        np.random.seed(0)
        npts = int(np.prod(ndivs))
        pts = np.random.uniform(low=[d[0] for d in domain], high=[d[1] for d in domain], size=(npts, ndim) )

        pot = lambda c:np.sum(c**2, axis=-1)/2

        wfns = DGB.run(pts, pot, alphas=1, clustering_radius=.05, min_singular_value=1e-4)
        # wfns = ham.get_wavefunctions()

        # wfns[0].project(1).plot().show()

        test_es = np.sort(np.sum(list(itertools.product(*[np.arange(8)+1/2]*ndim)), axis=-1))
        self.assertLess(
            np.linalg.norm(
                wfns.energies[:8] - test_es[:8]
            ),
            .05 # not sure why it's further off but it's still giving qualitatively correct results for now
        )

    @validationTest
    def test_Morse(self):
        d = 2
        ndivs = [15]*d
        domain = [[-.45, 2]]*d # , [-1, 1], [-1, 1]] # unit cube
        ndim = len(domain)

        pts = np.array(
            np.meshgrid(*(np.linspace(d[0], d[1], n) for d, n in zip(domain, ndivs)))
        ).T.reshape(-1, ndim)

        w = 2.05; wx=.1; mu=1
        de = (w ** 2) / (4 * wx)
        a = np.sqrt(2 * mu * wx)

        pot = lambda c,de=de,a=a,mu=mu: np.sum(de*(1-np.exp(-a*c))**2, axis=-1)

        # ham = DGB(pts, pot, alphas=1, optimize_centers=False, clustering_radius=.01, min_singular_value=1e-4)
        #
        # raise Exception(ham.V)

        wfns = DGB.run(pts, pot, alphas=2, optimize_centers=False, clustering_radius=.01, min_singular_value=1e-4)
        e = wfns.energies

        base_energies = [w*(np.arange(3) + 1 / 2) - wx*(np.arange(3) + 1 / 2)**2] * ndim

        test_es = np.sort(np.sum(list(itertools.product(*base_energies)), axis=-1))

        self.assertLess(
            np.linalg.norm(
                e[:3] - test_es[:3]
            ),
            .0025,
            msg="{} != {}".format(e[:3], test_es[:3])
        )

    @validationTest
    def test_SampleRotated(self):
        ndim = d = 2
        reduced_mass = (
                               1 / (16 * UnitsData.convert("AtomicMassUnits", "ElectronMass"))
                               + 1 / (1.01 * UnitsData.convert("AtomicMassUnits", "ElectronMass"))
                       ) ** (-1)
        masses = [reduced_mass] * ndim
        re = np.array([0.957836 * UnitsData.convert("Angstroms", "BohrRadius"), ] * ndim)

        w = 1#1750.47 * UnitsData.convert("Wavenumbers", "Hartrees")
        wx = .1#2 * 84.11 * UnitsData.convert("Wavenumbers", "Hartrees")
        mu = 1#reduced_mass
        de = (w ** 2) / (4 * wx)
        a = .5#np.sqrt(2 * mu * wx)

        ang = np.deg2rad(45)
        rot_mat = np.array([
            [np.cos(ang), np.sin(ang)],
            [-np.sin(ang), np.cos(ang)],
        ])

        def simple_morse(c, de=de, a=a, deriv_order=None):
            base = c.shape[:-2]
            ndim = c.shape[-1]

            c = c.reshape(-1, ndim)
            if deriv_order is None:
                res = np.sum(de * (1 - np.exp(-a * c)) ** 2, axis=-1)
            else:
                n = deriv_order
                m = ((-1) ** (n + 1) * 2 * a ** n * de) * np.exp(-2 * a * c) * (np.exp(a * c) - (2 ** (n - 1)))
                if n == 1:
                    return m
                res = np.zeros(c.shape[:-1] + (ndim,) * deriv_order)
                for k in range(ndim):
                    idx = (...,) + (k,) * n
                    res[idx] = m[..., k]
                # for j in range(deriv_order):
                #     res = np.tensordot(res, rot_mat, axes=[1, 1])
                # res = res

            res = res.reshape(base + res.shape[1:])

            return res

        test_ham = DGB.construct(
            np.array([
                [-1,  1],
                [ 1, .5]
            ]),
            simple_morse,
            alphas=[[1, .1], [1, .1]],
            masses=[1],
            expansion_degree=2,
            transformations=np.array([
                np.eye(2),

                np.array([
                    [ 1 / np.sqrt(2), -1 / np.sqrt(2)],
                    [ 1 / np.sqrt(2),  1 / np.sqrt(2)]
                ]).T
            ])
        )

        self.assertTrue(
            np.allclose(
                test_ham.gaussians.coords.centers,
                np.array([
                    [-1, 1],
                    [ 1,.5]
                ])
            )
        )
        self.assertTrue(
            np.allclose(
                test_ham.gaussians.alphas,
                np.array([
                    [1, .1],
                    [1, .1]
                ])
            )
        )
        self.assertTrue(
            np.allclose(
                test_ham.S,
                [[1, 0.2859716],
                 [0.2859716, 1]]
            )
        )
        self.assertTrue(
            np.allclose(
                test_ham.T,
                [[ 0.55,      -0.03266702],
                 [-0.03266702, 0.55      ]]
            )
        )
        self.assertTrue(
            np.allclose(
                test_ham.V,
                [[2.23291299,  0.66597918],
                 [0.66597918,  0.99361619]]
            )
        )

        # A = np.array([
        #     [-1,  0, 0],
        #     [ 1, -3, 0],
        #     [ 0,  0, 1]
        # ])
        # _, tf = np.linalg.eigh(A.T@A)
        # # need to resort the columns of tf to
        # # match my analytic test
        # tf = tf[:, (2, 1, 0)]
        # if np.linalg.det(tf) < 0:
        #     tf[:, 2] *= -1
        # test_ham = DGB(
        #     np.array([
        #         [-1,  1, 1],
        #         [ 1, -2, 0]
        #     ]),
        #     simple_morse,
        #     optimize_centers=False,
        #     alphas=[[1, .1, .5], [1, .5, .5]],
        #     clustering_radius=-1,
        #     min_singular_value=-10,  # 0.0001,
        #     num_svd_vectors=10000,
        #     expansion_degree=2,
        #     transformations=np.array([
        #         np.eye(3),
        #         # np.array([
        #         #     [0.85065080835204, -0.5257311121191337, 0.0],
        #         #     [0.0, 0.0, 1.0 ],
        #         #     [0.5257311121191336, 0.8506508083520398, 0.0]
        #         # ])
        #         tf.T
        #         # np.eye(3)
        #     ])
        # )
        #
        # self.assertTrue(
        #     np.allclose(
        #         test_ham.centers,
        #         np.array([
        #             [-1, 1, 1],
        #             [ 1, -2, 0]
        #         ])
        #     ) and
        #     np.allclose(
        #         test_ham.alphas,
        #         np.array([
        #             [ 1, .1, .5],
        #             [ 1, .5, .5]
        #         ])
        #     ) and
        #     np.allclose(
        #         test_ham.S,
        #         [[1, 0.05927224],
        #          [0.05927224, 1]]
        #     ) and
        #     np.allclose(
        #         test_ham.T,
        #         [[0.8, -0.04597853],
        #          [-0.04597853, 1.0]]
        #     ) and
        #     np.allclose(
        #         test_ham.V,
        #         [[2.66034212,  0.42746656],
        #          [0.42746656, 11.74998476]]
        #     )
        # )

    w2h = UnitsData.convert("Wavenumbers", "Hartrees")
    @staticmethod
    def multiply_model_functions(
            func1: 'dict[int|tuple, dict|tuple]',
            func2: 'dict[int|tuple, dict|tuple]'
    ):
        new_func = {}
        for k, f1 in func1.items():
            for k2, f2 in func2.items():
                if isinstance(k, int):
                    k = (k,)
                    f1 = (f1,)
                if isinstance(k2, int):
                    k2 = (k2,)
                    f2 = (f2,)
                new_func[k + k2] = f1 + f2
        return new_func
    @classmethod
    def buildWaterModel(cls,*,
                        oh_model=False,
                        atoms=None,
                        w=3869.47 * w2h,
                        wx=84 * w2h,
                        w2=3869.47 * w2h,
                        wx2=84 * w2h,
                        ka=1600 ** 2 / 150 * w2h,
                        dudr1=None, dudr2=None, duda=None,
                        dipole=None,
                        dipole_magnitude=None,
                        dipole_direction='auto'
                        ):
        base_water = Molecule.from_file(
            TestManager.test_data('water_freq.fchk'),
            internals=[
                [0, -1, -1, -1],
                [1,  0, -1, -1],
                [2,  0, 1, -1],
            ]
        )

        if oh_model:
            w2 = None
            wx2 = None
            ka = None
            mol = Molecule(
                base_water.atoms[:2],
                # ['O', 'O'],
                base_water.coords[:2],
                internals=[[0, -1, -1, -1], [1, 0, -1, -1]]
            ).get_embedded_molecule(load_properties=False)
        elif atoms is not None:
            mol = Molecule(
                atoms,
                # ['O', 'O'],
                base_water.coords,
                internals=[
                    [0, -1, -1, -1],
                    [1,  0, -1, -1],
                    [2,  0,  1, -1],
                ]
            ).get_embedded_molecule(load_properties=False)
        else:
            mol = base_water.get_embedded_molecule(load_properties=False)

        r1 = 0
        r2 = 1
        a12 = 2
        potential_params={}
        if wx is not None:
            potential_params[r1]={'morse':{'w':w, 'wx':wx}}
        elif w is not None:
            potential_params[r1]={'harmonic':{'k':w}}
        if wx2 is not None:
            potential_params[r2]={'morse':{'w':w2, 'wx':wx2}}
        elif w2 is not None:
            potential_params[r2]={'harmonic':{'k':w2}}
        elif w2 is not None:
            potential_params[r2]={'harmonic':{'k':w2}}
        if ka is not None:
            potential_params[a12]={'harmonic':{'k':ka}}

        if dipole is None and dudr1 is not None or dudr2 is not None or duda is not None:
            if oh_model:
                dipole = [
                    {r1: {'linear': {'eq': 0, 'scaling': dudr1}}},
                    0,
                    0
                ]
            else:
                dipole = [
                    {
                        (r1, a12): ({'linear': {'eq': 0, 'scaling': dudr1}}, {'sin': {'eq': 0}}),
                        (r2, a12): ({'linear': {'eq': 0, 'scaling': dudr2}}, {'sin': {'eq': 0, 'scaling':-1}})
                    },
                    {
                        (r1, a12): ({'linear': {'eq': 0, 'scaling': dudr1}}, {'cos': {'eq': 0, 'scaling':1/2}}),
                        (r2, a12): ({'linear': {'eq': 0, 'scaling': dudr2}}, {'cos': {'eq': 0, 'scaling':1/2}})
                    },
                    0
                ]

        # if dipole is not None:
        #     if isinstance(dipole, str) and dipole == 'auto':
        #         dipole_magnitude = 'auto'
        # elif dudr1 is not None or dudr2 is not None or duda is not None:
        #
        #     dipole_magnitude = {}
        #     if dudr1 is not None:
        #         dipole_magnitude[r1] = {'linear': {'eq': 0, 'scaling': dudr1}}
        #     if dudr2 is not None:
        #         dipole_magnitude[r2] = {'linear': {'eq': 0, 'scaling': dudr2}}
        #     if duda is not None:
        #         dipole_magnitude[a12] = {'linear': {'eq': 0, 'scaling': duda}}
        # if dipole_magnitude is not None:
        #     if isinstance(dipole_magnitude, str) and dipole_magnitude == 'auto':
        #         dipole_magnitude = {
        #             r1: {'linear': {'eq': 0, 'scaling': 1 / 5.5}},
        #             r2: {'linear': {'eq': 0, 'scaling': 1 / 5.5}}
        #         }
        #     if isinstance(dipole_direction, str) and dipole_direction == 'auto':
        #         if oh_model:
        #             dipole_direction = [
        #                 1,
        #                 0,
        #                 0
        #             ]
        #         else:
        #             dipole_direction = [
        #                 {
        #                     a12:{'sin': {'eq': 0}}
        #                 },
        #                 {
        #                     a12: {'cos': {'eq': 0}}
        #                 },
        #                 0
        #             ]
        #     dipole = [
        #         0
        #             if isinstance(d, int) and d == 0 else
        #         dipole_magnitude
        #             if isinstance(d, int) and d == 1 else
        #         cls.multiply_model_functions(dipole_magnitude, d)
        #         for d in dipole_direction
        #     ]

        return mol, mol.get_model(
            potential_params,
            dipole=dipole
        )
    @classmethod
    def buildTrajectory(cls,
                        mol,
                        cart_pot_func,
                        steps,
                        timestep=.5,
                        initial_energies=None,
                        initial_displacements=None,
                        displaced_coords=None,
                        seed=0
                        ):

        if initial_displacements is not None:
            init_pos = mol.get_displaced_coordinates(
                initial_displacements,
                which=displaced_coords,
                internals='reembed'
            )

            sim = AIMDSimulator(
                mol.masses,
                init_pos,
                lambda c: -cart_pot_func(c, deriv_order=1)[1].reshape(c.shape),
                timestep=timestep,
                track_kinetic_energy=True
            )
        else:
            mol.potential_derivatives = cart_pot_func(mol.coords, deriv_order=2)[1:]
            nms = mol.normal_modes.modes.basis
            sim = AIMDSimulator(
                mol.atomic_masses,
                [mol.coords] * len(initial_energies),
                lambda c: -cart_pot_func(c, deriv_order=1)[1].reshape(c.shape),
                velocities=AIMDSimulator.mode_energies_to_velocities(nms.inverse.T, mol.atomic_masses, initial_energies, inverse=nms.matrix.T),
                timestep=timestep,
                track_kinetic_energy=True
            )

        sim.propagate(steps)
        raise Exception(np.array(sim.trajectory).shape)
        coords = np.array(sim.trajectory).reshape((-1,) + mol.coords.shape)
        coords = mol.embed_coords(coords)

        return coords, sim

    @staticmethod
    def buildRunDGB(
            coords,
            pot,
            dipole,
            *,
            logger=True,
            plot_wavefunctions=True,
            plot_spectrum=True,
            **opts
    ):

        dgb = DGB.construct(
            np.round(coords, 8),  # this ends up really mattering to keep optimize_centers stable
            pot,
            logger=logger,
            **opts
        )

        logger = dgb.logger
        with logger.block(tag="Running DGB"):
            logger.log_print("num coords: {c}", c=len(dgb.gaussians.coords.centers))
            with logger.block(tag="S"):
                logger.log_print(logger.prep_array(dgb.S[:5, :5]))
            with logger.block(tag="T"):
                logger.log_print(logger.prep_array(dgb.T[:5, :5]))
            with logger.block(tag="V"):
                logger.log_print(logger.prep_array(dgb.V[:5, :5]))

            wfns_cart = dgb.get_wavefunctions()
            with logger.block(tag="Energies"):
                logger.log_print(
                    logger.prep_array(wfns_cart.energies[:5] * UnitsData.convert("Hartrees", "Wavenumbers"))
                )
            with logger.block(tag="Frequencies"):
                logger.log_print(
                    logger.prep_array(wfns_cart.frequencies()[:5] * UnitsData.convert("Hartrees", "Wavenumbers"))
                )

            if plot_wavefunctions:
                for i in range(4):
                    wfns_cart[i].plot().show()

            spec = wfns_cart[:4].get_spectrum(dipole)
            with logger.block(tag="Intensities"):
                logger.log_print(
                    logger.prep_array(spec.intensities)
                )
            if plot_spectrum:
                spec.plot().show()

    @classmethod
    def setupCartesianModelDGB(cls):
        ...

    @validationTest
    def test_ModelPotentialAIMD1DOG(self):

        mol, model = self.buildWaterModel(
            oh_model=True,
            dudr1=1/2
        )

        check_freqs = False
        if check_freqs:
            freqs = model.normal_modes()[0]
            raise Exception(freqs * UnitsData.convert("Hartrees", "Wavenumbers"))

        check_anh = False
        if check_anh:
            model.run_VPT(order=2, states=5)
            """
            ============================================= IR Data ==============================================
            Initial State: 0 
                           Harmonic                  Anharmonic
            State   Frequency    Intensity       Frequency    Intensity
              1    3869.47000    257.06585      3701.47000    251.27204
              2    7738.94000      0.00000      7234.94000      5.21706
              3   11608.41000      0.00000     10600.41000      0.22125
              4   15477.88000      0.00000     13797.88000      0.00000
              5   19347.35000      0.00000     16827.35000      0.00000
            ====================================================================================================
            """
            raise Exception(...)

        cart_pot_func = model.potential
        cart_dipole_func = lambda coords,deriv_order=None: (
            model.dipole(coords)
                if deriv_order is None else
            [np.moveaxis(d, -1, 1) for d in model.dipole(coords, deriv_order=deriv_order)]
        )

        pot_derivs = model.v(order=8, evaluate='constants', lambdify=True)
        def r_pot_func(rs, deriv_order=None):
            rs = rs.reshape(-1, 1)
            if deriv_order is None:
                return pot_derivs[0](rs)
            elif deriv_order > 8:
                raise ValueError("only imp'd up to 8th order")
            else:
                return [p(rs) for p in pot_derivs[:deriv_order+1]]
        mu_derivs = model.mu(order=8, evaluate='constants', lambdify=True)
        def r_dipole_func(rs, deriv_order=None):
            rs = rs.reshape(-1, 1)
            if deriv_order is None:
                u = np.moveaxis(
                        np.array([m[0](rs) for m in mu_derivs]),
                        0, 1
                    )
                return u
            elif deriv_order > 8:
                raise ValueError("only imp'd up to 8th order")
            else:
                u = [
                    np.moveaxis(
                        np.array([m[i](rs) for m in mu_derivs]),
                        0, 1
                    )
                    for i in range(deriv_order+1)
                ]
                return u

        # plot_points = False
        # if plot_points:
        #     ke_list = np.array(sim.kinetic_energies).flatten()
        #     r_vals = Molecule(mol.atoms, coords, internals=mol.internals).internal_coordinates[:, 1, 0]
        #
        #     ke_plot = plt.Plot(
        #         r_vals,
        #         ke_list
        #     )
        #
        #     pe_list = cart_pot_func(coords)
        #     plt.Plot(r_vals, pe_list, figure=ke_plot)
        #     plt.ScatterPlot(r_vals, pe_list, figure=ke_plot)
        #     plt.Plot(r_vals, pe_list + ke_list, figure=ke_plot).show()


        # coords = mol.get_displaced_coordinates(
        #     np.linspace(-.65, .65, 51).reshape(-1, 1),
        #     which=[[1, 0]],
        #     internals='reembed'
        # )

        coords, sim = self.buildTrajectory(mol, cart_pot_func,
                                           50,
                                           timestep=5,
                                           initial_energies=[[15000*self.w2h], [-15000*self.w2h]]
                                           )


        """
        TS = 5, 10000 KE
        25:
        Energies: [ 1913.73487519  5615.20480647  9148.67494111 12514.14640237 15711.62241258]
        1D Freqs: [ 3701.46993128  7234.94006592 10600.41152718 13797.88753739 16827.3755888 ]
        50:
        Energies: [ 1913.73487989  5615.20487993  9148.67555268 12514.14981935 15711.63643239]
        1D Freqs: [ 3701.47000004  7234.94067279 10600.41493946 13797.9015525  16827.4240311 ]
        
        TS = 1, 10000 KE
        25:
        Energies: [ 1913.76077586  5616.47422902  9162.75780088 12701.98968031 16420.3295596 ]
        1D Freqs: [ 3702.71345316  7248.99702502 10788.22890445 14506.56878375 19188.42196799]
        50:
        Energies: [ 1913.73489262  5615.20504874  9148.67934699 12514.20896499 15714.09594936]
        1D Freqs: [ 3701.47015612  7234.94445437 10600.47407237 13800.36105674 16859.6813473 ]
"""

        plot_points = False
        if plot_points:
            r_vals = Molecule(mol.atoms, coords, internals=mol.internals).internal_coordinates[:, 1, 0]
            # raise Exception(np.min(r_vals), np.max(r_vals))
            pe_list = cart_pot_func(coords) * UnitsData.hartrees_to_wavenumbers
            r_sort = np.argsort(r_vals)
            pe_plot = plt.Plot(r_vals[r_sort], pe_list[r_sort])
            plt.ScatterPlot(r_vals, pe_list, figure=pe_plot).show()
            raise Exception(...)

        sort_points = True
        if sort_points:
            # r_vals = Molecule(mol.atoms, coords, internals=mol.internals).internal_coordinates[:, 1, 0]
            # raise Exception(np.min(r_vals), np.max(r_vals))
            pe_list = cart_pot_func(coords) * UnitsData.hartrees_to_wavenumbers
            coords = coords[np.argsort(pe_list)]

        from McUtils.Data import PotentialData
        freq = 3869.47 * self.w2h
        anh = 84 * self.w2h
        De = (freq ** 2) / (4 * anh)
        muv = (1/model.vals[model.m(0)] + 1/model.vals[model.m(1)])
        a = np.sqrt(2 * anh / muv)
        re = model.vals[model.r(0, 1)]

        def morse_basic(r,
                        re=re,
                        alpha=a,
                        De=De,
                        deriv_order=None,
                        _morse=PotentialData["MorsePotential"]
                        ):
            return _morse(r, re=re, alpha=alpha, De=De, deriv_order=deriv_order)

        pairwise_potential_functions = {
            (0, 1):morse_basic
        }

        mass_vec = np.array(
            [
                mol.masses[0] * UnitsData.convert("AtomicMassUnits", "ElectronMass"),
                mol.masses[1] * UnitsData.convert("AtomicMassUnits", "ElectronMass")
            ]
        )

        r_vals = np.abs(coords.reshape(-1, 6)[:, 3] - coords.reshape(-1, 6)[:, 0])
        red_mass = 1/(1/mass_vec[0] + 1/mass_vec[1])

        r_crd = r_vals.view(np.ndarray).reshape(-1, 1, 1)
        self.buildRunDGB(
            r_crd,
            r_pot_func,
            r_dipole_func,
            alphas='virial'
            # np.sqrt(np.abs(red_mass*r_pot_func(r_crd, deriv_order=2)[2].reshape(-1, 1)))
            # alphas=100,
            , optimize_centers=True
            , masses=[red_mass]
            , expansion_degree=2
            # , quadrature_degree=4
            # quadrature_degree=3,
            # min_singular_value=1e-8,
            , plot_wavefunctions=True
            , plot_spectrum=False
        )
        """
        >>------------------------- Running DGB -------------------------
:: num coords: 16
::> S
  > [[1.         0.9925669  0.97109643 0.87837011 0.84262766]
  >  [0.9925669  1.         0.93587884 0.9276011  0.77948889]
  >  [0.97109643 0.93587884 1.         0.7553171  0.94236045]
  >  [0.87837011 0.9276011  0.7553171  1.         0.55475199]
  >  [0.84262766 0.77948889 0.94236045 0.55475199 1.        ]]
<::
::> T
  > [[ 0.0088153   0.0087815   0.00775012  0.00618133  0.00442927]
  >  [ 0.0087815   0.00915033  0.00701388  0.00762487  0.00319938]
  >  [ 0.00775012  0.00701388  0.00815306  0.00304996  0.00636977]
  >  [ 0.00618133  0.00762487  0.00304996  0.01022226 -0.00076633]
  >  [ 0.00442927  0.00319938  0.00636977 -0.00076633  0.00721211]]
<::
::> V
  > evauating integrals with 2-degree expansions
  > expanding about 136 points...
  > adding up all derivative contributions...
  > [[0.00220383 0.00226341 0.00218052 0.00267003 0.00224616]
  >  [0.00226341 0.00241846 0.00205612 0.00312139 0.00189822]
  >  [0.00218052 0.00205612 0.002543   0.00194964 0.0031494 ]
  >  [0.00267003 0.00312139 0.00194964 0.00489011 0.00126777]
  >  [0.00224616 0.00189822 0.0031494  0.00126777 0.00471625]]
<::
:: solving with min_singular_value=0.0001
:: solving with subspace size 15
::> Energies
  > [ 1903.77688585  5604.13185152  9135.56244159 11562.10564082 12497.8951657 ]
<::
::> Frequencies
  > [ 3700.35496567  7231.78555574  9658.32875497 10594.11827985 13784.63473029]
<::
:: evauating integrals with 2-degree expansions
:: expanding about 136 points...
:: adding up all derivative contributions...
::> Intensities
  > [2.51317857e+02 5.53642925e+00 1.41867406e-01]
<::
"""

        self.buildRunDGB(
            coords,
            cart_pot_func,
            cart_dipole_func,
            masses=mass_vec,
            modes='normal'
            # , transformations='reaction_path',
            , optimize_centers=True
            , alphas='virial'
            # , quadrature_degree=3
            , expansion_degree=2
            , pairwise_potential_functions=pairwise_potential_functions
            , plot_wavefunctions=False
            , plot_spectrum=False
        )
        """
        >>------------------------- Running DGB -------------------------
:: num coords: 16
::> S
  > [[1.         0.9925669  0.97109642 0.87837011 0.84262766]
  >  [0.9925669  1.         0.93587882 0.9276011  0.77948889]
  >  [0.97109642 0.93587882 1.         0.75531707 0.94236046]
  >  [0.87837011 0.9276011  0.75531707 1.         0.55475199]
  >  [0.84262766 0.77948889 0.94236046 0.55475199 1.        ]]
<::
::> T
  > [[ 0.0088153   0.0087815   0.00775012  0.00618133  0.00442927]
  >  [ 0.0087815   0.00915033  0.00701387  0.00762487  0.00319938]
  >  [ 0.00775012  0.00701387  0.00815306  0.00304996  0.00636977]
  >  [ 0.00618133  0.00762487  0.00304996  0.01022226 -0.00076633]
  >  [ 0.00442927  0.00319938  0.00636977 -0.00076633  0.00721211]]
<::
::> V
  > evauating integrals with 2-degree expansions
  > expanding about 136 points...
  > adding up all derivative contributions...
  > [[0.00224589 0.00230473 0.00222236 0.00270563 0.00228409]
  >  [0.00230473 0.00245964 0.00209605 0.00315856 0.00193301]
  >  [0.00222236 0.00209605 0.00258708 0.00198109 0.00319271]
  >  [0.00270563 0.00315856 0.00198109 0.00492896 0.00129199]
  >  [0.00228409 0.00193301 0.00319271 0.00129199 0.00476405]]
<::
:: solving with min_singular_value=0.0001
:: solving with subspace size 15
::> Energies
  > [ 1913.69991002  5615.15838793  9148.63078537 12514.18923649 15711.99693727]
<::
::> Frequencies
  > [ 3701.45847791  7234.93087535 10600.48932646 13798.29702724 16828.69542188]
<::
:: evauating integrals with 2-degree expansions
:: expanding about 136 points...
:: adding up all derivative contributions...
::> Intensities
  > [251.23202594   5.56303597   0.25175567]
<::
>>--------------------------------------------------<<
"""

        self.buildRunDGB(
            np.round(coords, 8), # this ends up really mattering to keep optimize_centers stable
            cart_pot_func,
            cart_dipole_func,
            masses=mass_vec
            , cartesians=[0]
            # , transformations='reaction_path',
            , optimize_centers=True
            # , alphas=[[160, 10]]
            , alphas='masses'
            # , alphas={'method':'virial', 'scaling':1}
            # , quadrature_degree=9
            , expansion_degree=2 # order 2 is less precise annoyingly...
            , pairwise_potential_functions={
                'functions':pairwise_potential_functions,
                'quadrature_degree':9
            }
            , plot_wavefunctions=False
            # , alphas=[.05]
        )

        """
>>------------------------- Running DGB -------------------------
:: num coords: 85
::> S
  > [[1.         0.99776758 0.99757665 0.99074271 0.9796392 ]
  >  [0.99776758 1.         0.9907248  0.97955825 0.99079762]
  >  [0.99757665 0.9907248  1.         0.99777641 0.96355397]
  >  [0.99074271 0.97955825 0.99777641 1.         0.94408894]
  >  [0.9796392  0.99079762 0.96355397 0.94408894 1.        ]]
<::
::> T
  > [[0.00544321 0.00541892 0.00541684 0.00534266 0.00522269]
  >  [0.00541892 0.00544321 0.00534247 0.00522181 0.00534326]
  >  [0.00541684 0.00534247 0.00544321 0.00541901 0.0050501 ]
  >  [0.00534266 0.00522181 0.00541901 0.00544321 0.00484321]
  >  [0.00522269 0.00534326 0.0050501  0.00484321 0.00544321]]
<::
::> V
  > evauating integrals with 2-degree expansions
  > expanding about 3655 points...
  > adding up all derivative contributions...
  > [[0.00758899 0.00788831 0.00731749 0.00709651 0.00860499]
  >  [0.00788831 0.00829758 0.00750701 0.00718572 0.00924987]
  >  [0.00731749 0.00750701 0.00715677 0.00703804 0.00799182]
  >  [0.00709651 0.00718572 0.00703804 0.00701469 0.0074619 ]
  >  [0.00860499 0.00924987 0.00799182 0.0074619  0.01071821]]
<::
:: solving with min_singular_value=0.0001
:: solving with subspace size 13
::> Energies
  > [ 2511.15650725  6213.70399517  9750.28681508 13138.44162533 16437.52203889]
<::
::> Frequencies
  > [ 3702.54748792  7239.13030783 10627.28511808 13926.36553164 17110.59884166]
<::
:: evauating integrals with 2-degree expansions
:: expanding about 3655 points...
:: adding up all derivative contributions...
::> Intensities
  > [251.31774692   5.6350833    0.26484567]
<::
"""

        # crd = nm_dgb.gaussians.coords.centers.flatten()
        # # plt.Plot(r_vals, crd.flatten()).show()
        # r_sort = np.argsort(r_vals)
        # slope = (np.diff(crd[r_sort]) / np.diff(r_vals[r_sort]))[0]

        # raise ValueError(
        #     (slope**2) * nm_dgb.pot.potential_function(nm_dgb.gaussians.coords.centers[:5], deriv_order=2)[2],
        #     r_pot_func(r_vals.reshape(-1, 1)[:5], deriv_order=2)[2]
        # )

        # for i in range(2):
        #     wfns_nm[i].plot().show()


    @classmethod
    def plot_dgb_potential(cls,
                           dgb, mol, potential,
                           coordinate_sel=None,
                           domain=None, domain_padding=1,
                           potential_cutoff=17000,
                           plot_cartesians=None,
                           plot_atoms=True,
                           cmap=None,
                           plot_points=100,
                           levels=24,
                           **plot_styles
                           ):
        def cutoff_pot(points, cutoff=potential_cutoff / UnitsData.hartrees_to_wavenumbers):
            values = potential(points)
            values[values > cutoff] = cutoff
            return values

        if isinstance(dgb, DGBCoords):
            coords = dgb
        else:
            coords = dgb.gaussians.coords

        if plot_cartesians is None:
            plot_cartesians = isinstance(coords, DGBCartesians)
            if coordinate_sel is None:
                coordinate_sel = [0, 1]
        if plot_cartesians:
            figure = mol.plot_molecule_function(
                cutoff_pot,
                axes=coordinate_sel,
                domain=domain,
                domain_padding=domain_padding,
                cmap=cmap,
                plot_points=plot_points,
                levels=levels,
                plot_atoms=plot_atoms,
                **plot_styles
            )
        else:
            if domain is None:
                from McUtils.Zachary import Mesh
                domain = Mesh(coords.centers).bounding_box
            points = DGBWavefunction.prep_plot_grid(
                domain,
                domain_padding=domain_padding,
            )
            if points.shape[-1] == 1:
                figure = plt.Plot(
                    *np.moveaxis(points, -1, 0),
                    potential(points),
                    **plot_styles
                )
            else:
                figure = plt.TriContourPlot(
                    *np.moveaxis(points, -1, 0),
                    cutoff_pot(points),
                    cmap=cmap,
                    levels=levels,
                    **plot_styles
                )

        return figure

    @classmethod
    def plot_gaussians(cls,
                       dgb, mol,
                       *,
                       domain=None,
                       domain_padding=1,
                       cmap='RdBu',
                       plot_dir=None,
                       plot_name='gaussian_{i}.pdf',
                       **plot_options
                       ):
        n = len(dgb.gaussians.coords.centers)
        wfns = DGBWavefunctions(
            np.zeros(n),
            np.eye(n),
            dgb
        )
        cls.plot_wavefunctions(
            wfns, dgb, mol,
            cmap=cmap,
            plot_name=plot_name,
            plot_dir=plot_dir,
            domain=domain,
            domain_padding=domain_padding,
            potential_styles=dict(
                domain=domain,
                domain_padding=domain_padding
            ),
            **plot_options
        )

    default_num_plot_wfns = 5
    @classmethod
    def plot_wavefunctions(cls,
                           wfns, dgb, mol,
                           which=True,
                           coordinate_sel=None,
                           cartesians=None,
                           plot_dir=None,
                           plot_name='wfn_{i}.pdf',
                           plot_potential=True,
                           plot_atoms=None,
                           plot_centers=True,
                           potential_styles=None,
                           **plot_options
                           ):

        figure = None
        if cartesians:
            wfns = wfns.as_cartesian_wavefunction()
            dgb = wfns.hamiltonian

        if coordinate_sel is None:
            coordinate_sel = list(range(dgb.gaussians.alphas.shape[-1]))

        figs = []
        if which is True:
            which = cls.default_num_plot_wfns
        if isinstance(which, int):
            which = list(range(which))
        if isinstance(which, list) and len(which) > 0:
            pot = dgb.pot.potential_function
            if plot_potential:
                if plot_atoms is None:
                    plot_atoms = bool(plot_centers)
                if potential_styles is None:
                    potential_styles = {}
                pot_figure = cls.plot_dgb_potential(
                    dgb, mol, pot,
                    plot_atoms=plot_atoms,
                    **potential_styles
                )
            else:
                pot_figure = None

            if plot_dir is not None:
                os.makedirs(plot_dir, exist_ok=True)
            if (
                    coordinate_sel is None and wfns.gaussians.alphas.shape[-1] == 1
                    or len(coordinate_sel) == 1
            ):
                for k in ['cmap', 'levels', 'plotter', 'contour_levels']:
                    if k in plot_options: del plot_options[k]
                    if k in potential_styles: del potential_styles[k]
            for i in which:
                if i < len(wfns):
                    if pot_figure is not None:
                        figure = pot_figure.copy()

                    if isinstance(dgb.gaussians.coords, DGBCartesians):
                        figs.append(
                            wfns[i].plot_cartesians(
                                coordinate_sel,
                                plot_centers=plot_centers,
                                figure=figure,
                                **plot_options
                            )
                        )
                    else:
                        if wfns.gaussians.alphas.shape[-1] > 1:
                            wfn = wfns[i]
                            if coordinate_sel is not None:
                                wfn = wfn.project(coordinate_sel)
                            figs.append(
                                wfn.plot(
                                    figure=figure,
                                    plot_centers=plot_centers,
                                    **plot_options
                                )
                            )
                        else:
                            figs.append(
                                wfns[i].plot(
                                    figure=figure,
                                    plot_centers=plot_centers,
                                    **plot_options
                                )
                            )

                    if plot_dir is not None:
                        fig = figs.pop()
                        fig.savefig(os.path.join(plot_dir, plot_name.format(i=i)))
                        fig.close()

            if plot_dir is None:
                figs[0].show()
    @classmethod
    def runDGB(cls,
               dgb: DGB,
               mol,
               plot_centers=True,
               plot_wavefunctions=True,
               plot_spectrum=False,
               pot_cmap='viridis',
               wfn_cmap='RdBu',
               wfn_points=100,
               wfn_contours=12,
               plot_dir=None,
               plot_potential=True,
               pot_points=100,
               domain=None,
               domain_padding=1,
               potential_cutoff=15000,
               mode=None,
               nodeless_ground_state=None,
               min_singular_value=None,
               subspace_size=None,
               plot_similarity=False,
               similarity_cutoff=None,
               similarity_chunk_size=None,
               similar_det_cutoff=None,
               **plot_options
               ):

        print("--->", len(dgb.gaussians.coords.centers))
        print(dgb.S[:5, :5])
        print(dgb.T[:5, :5])
        print(dgb.V[:5, :5])

        try:
            wfns, spec = dgb.run(
                calculate_spectrum=plot_spectrum,
                mode=mode,
                nodeless_ground_state=nodeless_ground_state,
                subspace_size=subspace_size,
                min_singular_value=min_singular_value,
                similarity_cutoff=similarity_cutoff,
                similarity_chunk_size=similarity_chunk_size,
                similar_det_cutoff=similar_det_cutoff
            )
        except Exception as e:

            if plot_wavefunctions is not False:
                print(e)

                if isinstance(plot_wavefunctions, str) and plot_wavefunctions == 'cartesians':
                    plot_wavefunctions = {'cartesians':None}
                cartesian_plot_axes = None
                if isinstance(plot_wavefunctions, dict):
                    if 'cartesians' in plot_wavefunctions:
                        cartesian_plot_axes=plot_wavefunctions['cartesians']
                        dgb = dgb.as_cartesian_dgb()
                        raise Exception(...)
                    else:
                        raise ValueError(plot_wavefunctions)

                pot = dgb.pot.potential_function
                figure = cls.plot_dgb_potential(
                    dgb, mol, pot,
                    domain=domain, domain_padding=domain_padding,
                    potential_cutoff=potential_cutoff,
                )

                dgb.gaussians.plot_centers(
                    figure,
                    xyz_sel=cartesian_plot_axes
                )

                figure.show()

            raise e
        else:
            if plot_similarity:
                plt.ArrayPlot(dgb.get_similarity_matrix()).show()

            if plot_spectrum:
                spec[:5].plot().show()

            use_cartesians = False
            if isinstance(plot_wavefunctions, str) and plot_wavefunctions == 'cartesians':
                plot_wavefunctions = {'cartesians': None}
            coordinate_sel = None
            if isinstance(plot_wavefunctions, dict):
                if 'cartesians' in plot_wavefunctions:
                    use_cartesians = True
                    coordinate_sel = plot_wavefunctions['cartesians']
                    plot_wavefunctions = plot_wavefunctions.get('num', True)
                elif 'modes' in plot_wavefunctions:
                    coordinate_sel = plot_wavefunctions['modes']
                    plot_wavefunctions = plot_wavefunctions.get('num', True)
                    plot_potential = False
                else:
                    raise ValueError(plot_wavefunctions)

            if plot_wavefunctions:
                cls.plot_wavefunctions(
                    wfns,
                    dgb,
                    mol,
                    cartesians=use_cartesians,
                    coordinate_sel=coordinate_sel,
                    plot_dir=plot_dir,
                    contour_levels=wfn_contours,
                    cmap=wfn_cmap,
                    plot_points=wfn_points,
                    plot_centers={'color': 'red'} if plot_centers else False,
                    domain=domain,
                    domain_padding=domain_padding,
                    plot_potential=plot_potential,
                    scaling=.2,
                    plotter=plt.TriContourLinesPlot,
                    potential_styles=dict(
                        domain=domain,
                        domain_padding=domain_padding,
                        cmap=pot_cmap,
                        plot_points=pot_points
                    ),
                    **plot_options
                )



            return wfns

    @classmethod
    def getMorseParameters(cls, w=None, wx=None, m1=None, m2=None, re=None):
        if w is None:
            w = 3869.47 * cls.w2h
        freq = w
        if wx is None:
            wx = 84 * cls.w2h
        anh = wx
        De = (freq ** 2) / (4 * anh)
        if m1 is None:
            m1 = AtomData["O", "Mass"] * UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
        if m2 is None:
            m2 = AtomData["H", "Mass"] * UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
        muv = (1 / m1 + 1 / m2)
        a = np.sqrt(2 * anh / muv)

        if re is None:
            re = 1.82534

        return (De, a, re)

    @classmethod
    def setupMorseFunction(cls, model, i, j, w=None, wx=None):

        from McUtils.Data import PotentialData

        if isinstance(model, float):
            m1 = model
            m2 = i
            re = j
        else:
            m1 = model.vals[model.m(i)]
            m2 = model.vals[model.m(j)]
            re = model.vals[model.r(i, j)]

        De, a, re = cls.getMorseParameters(w=w, wx=wx, m1=m1, m2=m2, re=re)

        def morse_basic(r,
                        re=re,
                        alpha=a,
                        De=De,
                        deriv_order=None,
                        _morse=PotentialData["MorsePotential"]
                        ):
            return _morse(r, re=re, alpha=alpha, De=De, deriv_order=deriv_order)

        return morse_basic

    @classmethod
    def plot_interpolation_error(cls, dgb, pot):
        sel = slice(None)  # slice(15,30)

        centers = dgb.gaussians.overlap_data['centers'][sel]
        embpot = dgb.gaussians.coords.embed_function(pot)
        realpots, _, real_hess = [d * 219475 for d in embpot(centers, deriv_order=2)]
        interpots, _, inter_hess = [d * 219475 for d in dgb.pot.potential_function(centers, deriv_order=2)]

        ords = np.argsort(realpots)
        devs = interpots[ords] - realpots[ords]
        rows, cols = np.triu_indices_from(dgb.S)
        utris = dgb.S[rows, cols]
        unscaled_devs = devs
        devs = devs * utris
        max_dev_pos = np.flip(np.argsort(np.abs(devs)))[:5]

        inter_trace = np.sum(np.sum(np.abs(inter_hess[ords]), axis=-1), axis=-1)
        real_trace = np.sum(np.sum(np.abs(real_hess[ords]), axis=-1), axis=-1)
        dev_trace = inter_trace - real_trace
        dev_hess = inter_hess - real_hess
        hess_plot = plt.ScatterPlot(realpots[ords], inter_trace - real_trace)

        print("Mean Absolute Error:", np.mean(np.abs(unscaled_devs)), "Std:", np.std(unscaled_devs))
        print("Mean Scaled Error:", np.mean(np.abs(devs)), "Std:", np.std(devs))
        print("Mean Hessian Error:", np.mean(np.abs(dev_hess.flatten())),
              "Std:", np.std(dev_hess.flatten()),
              "Max:", np.max(np.abs(dev_hess.flatten()))
              )
        print("Mean Summed Hessian Error:", np.mean(np.abs(dev_trace.flatten())),
              "Std:", np.std(dev_trace.flatten()),
              "Max:", np.max(np.abs(dev_trace.flatten())),
              )
        print("Maximum (Scaled) Error:", devs[max_dev_pos])
        print("Maximum (Scaled) Interpolation Error:")
        for l, r, c, tt, ii, ov in zip(
                dgb.gaussians.coords.centers[rows[sel][ords[max_dev_pos]]],
                dgb.gaussians.coords.centers[cols[sel][ords[max_dev_pos]]],
                dgb.gaussians.overlap_data['centers'][sel][ords[max_dev_pos]],
                realpots[ords[max_dev_pos]],
                interpots[ords[max_dev_pos]],
                utris[sel][ords[max_dev_pos]]
        ):
            print(f"Centers: {c} ({ov}) <- {l} {r}")
            print(f"  Error: {ii - tt} <- {tt} {ii}")

        bad_bad = np.abs(devs) > 50

        # center_plot=plt.ScatterPlot(
        #     centers[:, 0],
        #     centers[:, 1],
        #     c=unscaled_devs
        # )
        # center_plot = plt.ScatterPlot(
        #     centers[:, 0][bad_bad],
        #     centers[:, 1][bad_bad],
        #     c='blue'
        # )
        # plt.ScatterPlot(
        #     dgb.gaussians.coords.centers[:, 0],
        #     dgb.gaussians.coords.centers[:, 1],
        #     c='red',
        #     figure=center_plot
        # )
        # woof = plt.ScatterPlot(realpots[ords], unscaled_devs / realpots[ords])
        scaled_dev_plot = plt.ScatterPlot(realpots[ords], unscaled_devs)
        dev_plot = plt.ScatterPlot(realpots[ords], devs
                                   # plot_range=[None, [-100, 100]]
                                   )
        dev_plot.show()

    @validationTest
    def test_ModelPotentialAIMD1D(self):
        mol, model = self.buildWaterModel(
            w2=None, wx2=None,
            ka=None,
            oh_model=True,
            dudr1=1 / 5.5,
            dudr2=None#1 / 5.5
            # dipole_direction=[1, 0, 0]
        )

        check_freqs = False
        if check_freqs:
            freqs = model.normal_modes()[0]
            raise Exception(freqs * UnitsData.convert("Hartrees", "Wavenumbers"))

        check_anh = False
        if check_anh:
            model.run_VPT(order=2, states=5, degeneracy_specs='auto', logger=True)
            """
            ZPE: 1934.73500   1913.73500
            ============================================= IR Data ==============================================
            Initial State: 0 
                           Harmonic                  Anharmonic
            State   Frequency    Intensity       Frequency    Intensity
              1    3869.47000     33.99218      3701.47000     33.22606
              2    7738.94000      0.00000      7234.94000      0.68986
              3   11608.41000      0.00000     10600.41000      0.02926
              4   15477.88000      0.00000     13797.88000      0.00000
              5   19347.35000      0.00000     16827.35000      0.00000
            ====================================================================================================
            """
            raise Exception(...)

        # mol.potential_derivatives = model.potential(mol.coords, deriv_order=2)[1:]
        # raise Exception(mol.coords, mol.normal_modes.modes)

        sim = model.setup_AIMD(
            initial_energies=[
                [ 3800 * 3 * self.w2h],
                [-3800 * 3 * self.w2h]
            ],
            timestep=10
        )
        sim.propagate(35)
        coords = sim.extract_trajectory(flatten=True, embed=mol.coords)

        cartesians = False
        with BlockProfiler(inactive=True):

            dgb = model.setup_DGB(
                np.round(coords, 8),
                optimize_centers=1e-8,
                # optimize_centers=False,
                alphas={'method':'virial', 'scaling':2},
                modes=None if cartesians else 'normal',
                cartesians=[0, 1] if cartesians else None,
                quadrature_degree=3,
                # expansion_degree=2,
                # pairwise_potential_functions={
                #     (0, 1): self.setupMorseFunction(model, 0, 1)
                #     # (0, 2): self.setupMorseFunction(model, 0, 2)
                # }
            )

            self.runDGB(dgb, mol,
                        domain_padding=10,
                        plot_spectrum=True,
                        plot_wavefunctions=False#{'cartesians':[0, 1]} if not cartesians else True
                        )

    @validationTest
    def test_ModelPotentialPhasedAIMD1D(self):
        mol, model = self.buildWaterModel(
            w2=None, wx2=None,
            ka=None,
            oh_model=True,
            dudr1=1 / 5.5,
            dudr2=None  # 1 / 5.5
            # dipole_direction=[1, 0, 0]
        )

        check_freqs = False
        if check_freqs:
            freqs = model.normal_modes()[0]
            raise Exception(freqs * UnitsData.convert("Hartrees", "Wavenumbers"))

        # mol.potential_derivatives = model.potential(mol.coords, deriv_order=2)[1:]
        # raise Exception(
        #     self.getMorseParameters(),
        #     mol.coords,
        #     mol.normal_modes.modes.basis.matrix
        # )

        check_anh = False
        if check_anh:
            model.run_VPT(order=2, states=8, degeneracy_specs=None, logger=True)
            """
            ZPE: 1934.73500   1913.73500
            ============================================= IR Data ==============================================
            Initial State: 0 
                           Harmonic                  Anharmonic
            State   Frequency    Intensity       Frequency    Intensity
              1    3869.47000     33.99218      3701.47000     33.22606
              2    7738.94000      0.00000      7234.94000      0.68986
              3   11608.41000      0.00000     10600.41000      0.02926
              4   15477.88000      0.00000     13797.88000      0.00000
              5   19347.35000      0.00000     16827.35000      0.00000
              6   23216.82000      0.00000     19688.82000      0.00000
              7   27086.29000      0.00000     22382.29000      0.00000
              8   30955.76000      0.00000     24907.76000      0.00000
            ====================================================================================================
            """
            raise Exception(...)

        # mol.potential_derivatives = model.potential(mol.coords, deriv_order=2)[1:]
        # raise Exception(mol.coords, mol.normal_modes.modes)

        sim = model.setup_AIMD(
            initial_energies=[
                [ 2*3701 * self.w2h],
                [-2*3701 * self.w2h]
            ],
            timestep=10,
            track_velocities=True
        )
        sim.propagate(15)
        coords, velocities = sim.extract_trajectory(flatten=True, embed=mol.coords)
        momenta = velocities * mol.masses[np.newaxis, :, np.newaxis]
        all_pots = model.potential(coords)
        pots_scaling = 1 + (all_pots - np.min(all_pots)) / (np.max(all_pots) - np.min(all_pots))
        # momenta = -model.potential(coords, deriv_order=1)[1]
        # momenta = np.full(velocities.shape, .0)

        cartesians = False
        with BlockProfiler(inactive=True):

            dgb = model.setup_DGB(
                np.round(coords[1:], 8),
                # optimize_centers=1e-8,
                optimize_centers=False,
                modes=None if cartesians else 'normal',
                # alphas={'method':'virial', 'scaling':pots_scaling[1:, np.newaxis]},
                cartesians=[0, 1] if cartesians else None,
                quadrature_degree=3,
                expansion_degree=-1,
                pairwise_potential_functions={
                    (0, 1): self.setupMorseFunction(model, 0, 1)
                    # (0, 2): self.setupMorseFunction(model, 0, 2)
                },
                # momenta=15*(100*momenta[1:]),
                # momenta=150*momenta[1:],
                # momenta=np.ones_like(150*momenta[1:]) / 10,
                # momenta=np.zeros_like(momenta[1:])
            )

            plot_gaussians = False
            if plot_gaussians:
                self.plot_gaussians(dgb, mol,
                                    which=[0, 1, 2, 5],
                                    domain_padding=10,
                                    scaling=.2
                                    )
                raise Exception(...)

            """
            >>------------------------- Running distributed Gaussian basis calculation -------------------------
            :: diagonalizing in the space of 15 S functions
            :: ZPE: 1874.6654445484646
            :: Frequencies: [ 3701.65896883  7231.8802379  10591.834868   13773.9302937  16763.72207843 19442.90514644 21366.89341503 23107.37699882 25564.19893302 28107.929736   30768.02886696 34146.63799842 38497.7450025  48064.99243391]
            >>--------------------------------------------------<<
            """

            # print(np.linalg.eigvalsh(dgb.S))
            # print(dgb.gaussians.coords.centers[:5])
            # print(dgb.gaussians.momenta[:5])
            # print(dgb.gaussians.alphas[:5])

            self.runDGB(dgb, mol,
                        domain_padding=10,
                        plot_spectrum=False,
                        plot_wavefunctions=True  # {'cartesians':[0, 1]} if not cartesians else True
                        )

    @validationTest
    def test_ModelPotentialAIMD2D(self):
        mol, model = self.buildWaterModel(
            # w2=None, wx2=None,
            ka=None,
            dudr1=1/5.5,
            dudr2=1/5.5
            # dipole_direction=[1, 0, 0]
        )

        check_freqs = False
        if check_freqs:
            freqs = model.normal_modes()[0]
            raise Exception(freqs * UnitsData.convert("Hartrees", "Wavenumbers"))

        check_anh = False
        if check_anh:
            from Psience.VPT2 import VPTRunner

            VPTRunner.run_simple(
                [mol.atoms, mol.coords],
                potential_derivatives=model.potential(mol.coords, deriv_order=4)[1:],
                dipole_derivatives=model.dipole(mol.coords, deriv_order=3),
                order=2, states=3,
                logger=True,
                degeneracy_specs='auto',
                calculate_intensities=True
            )

            # model.run_VPT(order=2, states=2, degeneracy_specs='auto', logger=True)
            """
            ZPE: 3869.37229   3814.75070 
            State     Frequency    Intensity       Frequency    Intensity
              0 1    3896.87028     64.98650      3726.72077     63.58462
              1 0    3841.87433      0.26781      3675.96177      0.25708
              0 2    7793.74057      0.00000      7415.60644      0.00317
              2 0    7683.74866      0.00000      7220.11984      0.00861
              1 1    7738.74461      0.00000      7236.25926      1.32009
              0 3   11690.61085      0.00000     10986.10686      0.00048
              3 0   11525.62299      0.00000     10895.67849      0.00008
              1 2   11635.61490      0.00000     10593.72292      0.00056
              2 1   11580.61894      0.00000     10596.33867      0.05435
            """
            raise Exception(...)

        check_dvr = False
        if check_dvr:
            raise Exception("do the proper 2D thing...")
            print("Running DVR...")
            dvr = model.setup_DVR(
                domain=[[1, 4], [1, 4], [np.deg2rad(60), np.deg2rad(160)]],
                divs=[800, 800, 800], po_divs=[25, 25, 25],
                potential_optimize=True,
                logger=True
            )
            po_data = dvr.run()
            """
            :: PotentialOptimizedDVR([WavefunctionBasisDVR(None, pts=20, pot=None), WavefunctionBasisDVR(None, pts=20, pot=None), WavefunctionBasisDVR(None, pts=20, pot=None)], pot=SympyExpr(0.203038951525208*(1 - exp(-0.813301570558368*sqrt(2)*(r[1,2] - 1.8253409520594)))**2 + 0.203038951525208*(1 - exp(-0.813301570558368*sqrt(2)*(r[2,3] - 1.82534095205941)))**2 + 0.255579575354735*(0.551593470847119*a[1,2,3] - 1)**2))
            :: g: [[SympyExpr(0.0005786177281533848), SympyExpr(3.42971451934982e-5*cos(a[1,2,3])), SympyExpr(-3.42971451934982e-5*sin(a[1,2,3])/r[2,3])], [SympyExpr(3.42971451934982e-5*cos(a[1,2,3])), SympyExpr(0.0005786177281533848), SympyExpr(0.0)], [SympyExpr(-3.42971451934982e-5*sin(a[1,2,3])/r[2,3]), SympyExpr(0.0), SympyExpr(0.000578617728153385/r[2,3]**2 - 6.85942903869965e-5*cos(a[1,2,3])/(r[1,2]*r[2,3]) + 0.000578617728153385/r[1,2]**2)]]
            :: mass: [None, None, None]
            :: g_deriv: [SympyExpr(0.0), SympyExpr(0.0), SympyExpr(6.85942903869965e-5*cos(a[1,2,3])/(r[1,2]*r[2,3]))]
            :: domain: [[1, 4], [1, 4], [1.0471975511965976, 2.792526803190927]]
            :: divs: [800, 800, 800]
            :: potential_function: SympyExpr(0.203038951525208*(1 - exp(-0.813301570558368*sqrt(2)*(r[1,2] - 1.8253409520594)))**2 + 0.203038951525208*(1 - exp(-0.813301570558368*sqrt(2)*(r[2,3] - 1.82534095205941)))**2 + 0.255579575354735*(0.551593470847119*a[1,2,3] - 1)**2)
            ::> constructing grid
             """
            raise Exception(
                po_data.wavefunctions.energies[1] * UnitsData.hartrees_to_wavenumbers,
                po_data.wavefunctions.frequencies() * UnitsData.hartrees_to_wavenumbers
            )
        # mol.potential_derivatives = model.potential(mol.coords, deriv_order=2)[1:]
        # raise Exception(mol.coords, mol.normal_modes.modes)

        symm, asymm = model.normal_modes()[0]
        sim = model.setup_AIMD(
            initial_energies=[
                [ symm,  asymm],
                [-symm,  asymm],
                [-symm, -asymm],
                [ symm, -asymm],
                # [2000 * self.w2h, 0],
                # [0, 2000 * self.w2h],
                # [-2000 * self.w2h, 0],
                # [0, -2000 * self.w2h],
                # [-15000 * self.w2h, -15000 * self.w2h],
                # [15000 * self.w2h, -15000 * self.w2h],
                # [10000 * self.w2h, 0],
                # [0, 10000 * self.w2h],
                # [-10000 * self.w2h, 0],
                # [0, -10000 * self.w2h]
            ],
            timestep=25,
            track_velocities = True
        )
        sim.propagate(25)
        coords, velocities = sim.extract_trajectory(flatten=True, embed=mol.coords)
        momenta = 250 * velocities * mol.masses[np.newaxis, :, np.newaxis]
        # momenta = None

        cartesians=False
        with BlockProfiler(inactive=True):

            dgb = model.setup_DGB(
                coords,
                optimize_centers=1e-14,
                # optimize_centers=False,
                modes=None if cartesians else 'normal',
                cartesians=[0, 1] if cartesians else None,
                quadrature_degree=3,
                expansion_degree=2,
                # pairwise_potential_functions={
                #     (0, 1):self.setupMorseFunction(model, 0, 1),
                #     (0, 2):self.setupMorseFunction(model, 0, 2)
                # },
                # transformations='diag'
                momenta=momenta
            )

            plot_gaussians = True
            if plot_gaussians:
                self.plot_gaussians(dgb, mol)
            raise Exception(...)




            # print(dgb.gaussians.coords.centers[:3])
            # print(dgb.gaussians.overlap_data['init_covariances'][:3])
            # print(dgb.gaussians.overlap_data['initial_momenta'][:3])

            """
            >>------------------------- Running distributed Gaussian basis calculation -------------------------
            :: diagonalizing in the space of 59 S functions
            :: ZPE: 3815.0254181251166
            :: Frequencies: [ 3676.14207931  3726.99361301  7221.08377964  7236.78275062  7418.52358929 10595.91399038 10597.43855691 10901.20911454 11001.80376921 13797.34317942 13799.82687628 14293.18090478 14329.25062911 14640.80903924 16841.23176496 16849.01710935 17531.49955265 17572.95249918 17891.12454688 18319.17520319]
            >>--------------------------------------------------<<
            """

            """
            >>------------------------- Running distributed Gaussian basis calculation -------------------------
            :: diagonalizing in the space of 41 S functions
            :: ZPE: 3815.1110786450845
            :: Frequencies: [ 3676.77263767  3727.70483878  7227.82029744  7237.90664201  7440.76152723 10601.58120583 10618.10732722 10920.95275887 11104.66728657 13808.26706889 13827.90646153 14416.22589799 14471.29425293 15204.74679781 16858.16199974 16902.38702413 17668.82463041 17864.53996371 18809.9243302  18939.68118366]
            >>--------------------------------------------------<<
            """

            self.runDGB(dgb, mol,
                        similarity_cutoff=.95,
                        plot_spectrum=False,
                        plot_wavefunctions=False
                        # mode='classic'
                        )

    @validationTest
    def test_ModelPotentialAIMD2D_HOD(self):
        mol, model = self.buildWaterModel(
            # w2=None, wx2=None,
            ka=None,
            atoms=['O', 'H', 'D'],
            w2=3869.47 * self.w2h / np.sqrt(2),
            wx2=84 * self.w2h / np.sqrt(2),
            dudr1=1 / 5.5,
            dudr2=1 / 5.5
            # dipole_direction=[1, 0, 0]
        )

        check_freqs = False
        if check_freqs:
            freqs = model.normal_modes()[0]
            raise Exception(freqs * UnitsData.convert("Hartrees", "Wavenumbers"))

        check_anh = False
        if check_anh:
            model.run_VPT(order=2, states=2, degeneracy_specs='auto', logger=True)
            """
            ZPE:     3302.64651                   3220.08468
                             Harmonic                  Anharmonic
            State     Frequency    Intensity       Frequency    Intensity
              0 1    3870.20673     34.00202      3702.08096     33.19339
              1 0    2735.08628     16.06164      2616.40747     15.75924
              0 2    7740.41347      0.00000      7236.30645      0.65964
              2 0    5470.17256      0.00000      5114.40644      0.37370
              1 1    6605.29301      0.00000      6317.94784      0.00436
            """
            raise Exception(...)

        # raise Exception(model.v(0, evaluate='constants'))
        # mol.potential_derivatives = model.potential(mol.coords, deriv_order=2)[1:]
        # raise Exception(mol.coords, mol.normal_modes.modes)

        sim = model.setup_AIMD(
            initial_energies=[
                [5000 * self.w2h, 5000 * self.w2h],
                [-5000 * self.w2h, 5000 * self.w2h],
                [-5000 * self.w2h, -5000 * self.w2h],
                [5000 * self.w2h, -5000 * self.w2h],
                # [2000 * self.w2h, 0],
                # [0, 2000 * self.w2h],
                # [-2000 * self.w2h, 0],
                # [0, -2000 * self.w2h],
                # [-15000 * self.w2h, -15000 * self.w2h],
                # [15000 * self.w2h, -15000 * self.w2h],
                # [10000 * self.w2h, 0],
                # [0, 10000 * self.w2h],
                # [-10000 * self.w2h, 0],
                # [0, -10000 * self.w2h]
            ],
            timestep=15
        )
        sim.propagate(25)
        coords = sim.extract_trajectory(flatten=True, embed=mol.coords)

        cartesians = False
        with BlockProfiler(inactive=True):

            dgb = model.setup_DGB(
                np.round(coords, 8),
                optimize_centers=1e-8,
                # optimize_centers=False,
                modes=None if cartesians else 'normal',
                cartesians=[0, 1] if cartesians else None,
                # quadrature_degree=3,
                expansion_degree=2,
                pairwise_potential_functions={
                    (0, 1): self.setupMorseFunction(model, 0, 1),
                    (0, 2): self.setupMorseFunction(model, 0, 2,
                                                    w=3869.47 * self.w2h / np.sqrt(2),
                                                    wx=84 * self.w2h / np.sqrt(2),
                                                    )
                }
            )

            self.runDGB(dgb, mol,
                        # vmin=-.05,
                        # vmax=.05,
                        # domain=[[-20, 20], [-20, 20]],
                        plot_wavefunctions={'cartesians':[0, 1]} if not cartesians else True,
                        plot_spectrum=False
                        # mode='classic'
                        )

    @validationTest
    def test_ModelPotentialAIMD2D_sym_bend(self):
        mol, model = self.buildWaterModel(
            # w2=None, wx2=None,
            # ka=None,
            dudr1=1 / 5.5,
            dudr2=1 / 5.5
            # dipole_direction=[1, 0, 0]
        )

        check_freqs = False
        if check_freqs:
            freqs = model.normal_modes()[0]
            raise Exception(freqs * UnitsData.convert("Hartrees", "Wavenumbers"))

        check_anh = False
        if check_anh:
            model.run_VPT(order=2, states=2,
                          logger=True,
                          mode_selection=[0],
                          degeneracy_specs='auto'
                          )
            """
            > State    Harmonic   Anharmonic     Harmonic   Anharmonic
                         ZPE          ZPE    Frequency    Frequency
            0 0   2732.22799   2656.66754            -            - 
            0 1            -            -   3843.25802   3760.42618 
            1 0            -            -   1621.19795   1608.98591 
            0 2            -            -   7686.51604   7438.96418 
            2 0            -            -   3242.39590   3207.43724 
            1 1            -            -   5464.45597   5368.24375
            """
            raise Exception(...)

        # mol.potential_derivatives = model.potential(mol.coords, deriv_order=2)[1:]
        # raise Exception(mol.coords, mol.normal_modes.modes)

        bend, symm, asymm = model.normal_modes()[0]
        sim = model.setup_AIMD(
            initial_energies=np.array([
                [ bend,  symm, 0],
                [ bend, -symm, 0],
                [-bend,  symm, 0],
                [-bend, -symm, 0],
                [ bend, 0,  asymm],
                [ bend, 0, -asymm],
                [-bend, 0,  asymm],
                [-bend, 0, -asymm]
            ]),
            timestep=15
        )
        sim.propagate(5)
        coords = sim.extract_trajectory(flatten=True, embed=mol.coords)

        cartesians = False
        with BlockProfiler(inactive=True):

            dgb = model.setup_DGB(
                np.round(coords, 8),
                optimize_centers=1e-6,
                # optimize_centers=False,
                modes=None if cartesians else 'normal',
                coordinate_selection=[0, 1],
                cartesians=[0, 1] if cartesians else None,
                quadrature_degree=3,
                # expansion_degree=2,
                # pairwise_potential_functions={
                #     (0, 1): self.setupMorseFunction(model, 0, 1),
                #     (0, 2): self.setupMorseFunction(model, 0, 2)
                # }
            )

            # submol, submodel_stretch = self.buildWaterModel(
            #     # w2=None, wx2=None,
            #     ka=None
            # )
            # submol, submodel_bend = self.buildWaterModel(
            #     w=0, wx=None,
            #     w2=0, wx2=None,
            #     # ka=None,
            # )

            # print(
            #     dgb.evaluate_multiplicative_operator(
            #         submodel_bend.potential,
            #         quadrature_degree=4
            #     )[:5, :5]
            # )
            # print(
            #     dgb.evaluate_multiplicative_operator(
            #         submodel_bend.potential,
            #         expansion_degree=2
            #     )[:5, :5]
            # )
            #     dgb.evaluate_multiplicative_operator(
            #         submodel_stretch.potential,
            #         expansion_degree=2,
            #         pairwise_functions={
            #             (0, 1): self.setupMorseFunction(model, 0, 1),
            #             (0, 2): self.setupMorseFunction(model, 0, 2)
            #         }
            #     )[:5, :5]

            # print(
            #     dgb.evaluate_multiplicative_operator(
            #         submodel_stretch.potential,
            #         quadrature_degree=3
            #     )[:5, :5] -
            #     dgb.evaluate_multiplicative_operator(
            #         submodel_stretch.potential,
            #         expansion_degree=2,
            #         pairwise_functions={
            #             (0, 1): self.setupMorseFunction(model, 0, 1),
            #             (0, 2): self.setupMorseFunction(model, 0, 2)
            #         }
            #     )[:5, :5]
            # )
            raise Exception(...)

            """
>>------------------------- Running distributed Gaussian basis calculation -------------------------
:: diagonalizing in the space of 27 S functions
:: ZPE: 2703.0532865370874
:: Frequencies: [ 1611.61510053  3215.4525035   3761.65133913  4822.94422724  5383.49202567  6528.60083868  7121.5055917   7461.58631351  8451.7611391   8802.65305029  9126.15952187 10068.00946755 10655.33557818 11035.91654539 11686.79214732 12368.50320681 12618.40116626 13452.85838251 13927.47046538 14512.77618706]
>>--------------------------------------------------<<
[[0.00623843 0.00730029 0.00521558 0.00503987 0.00518823]
 [0.00730029 0.03210892 0.00082444 0.0032473  0.01149502]
 [0.00521558 0.00082444 0.02563066 0.00815058 0.00171581]
 [0.00503987 0.0032473  0.00815058 0.01458135 0.00088887]
 [0.00518823 0.01149502 0.00171581 0.00088887 0.01424601]]
 
            evauating integrals with 3-order quadrature
[[0.00624626 0.00731284 0.00521744 0.00504326 0.00519344]
 [0.00731284 0.03218069 0.00082591 0.00325166 0.0115091 ]
 [0.00521744 0.00082591 0.0256305  0.00815148 0.00171665]
 [0.00504326 0.00325166 0.00815148 0.01458469 0.00088981]
 [0.00519344 0.0115091  0.00171665 0.00088981 0.01425408]]
 
            
            """

            # type(self).default_num_plot_wfns = 5
            type(self).default_num_plot_wfns = 1
            self.runDGB(dgb, mol,
                        # vmin=-.05,
                        # vmax=.05,
                        # domain=[[-2, 2], [-2, 2]],
                        # plot_wavefunctions=False,
                        plot_wavefunctions={'cartesians': [0, 1]} if not cartesians else True
                        )

    def symmetrizePoints(self, coords, equivalent_atoms, unique_tol=7):
        coords = np.asanyarray(coords)
        all_perms = [
            np.unique(list(itertools.permutations(a)), axis=0)
            for a in equivalent_atoms
        ]
        stacks = [coords]
        for swaps in itertools.product(*all_perms):
            perm = np.arange(coords.shape[1])
            for og,sw in zip(equivalent_atoms, swaps):
                perm[og] = perm[sw]
            stacks.append(coords[:, perm])
        coords = np.concatenate(stacks, axis=0)
        if unique_tol is not None:
            _, upos = np.unique(np.round(coords, unique_tol), return_index=True, axis=0)
            coords = coords[np.sort(upos)]
        return coords
    @validationTest
    def test_ModelPotentialAIMD3D(self):

        # raise Exception(
        #     self.symmetrizePoints(
        #         [[0, 1, -1], [0, 1, 1]],
        #         [[1, 2]]
        #     )
        # )

        mol, model = self.buildWaterModel(
            # w2=None, wx2=None,
            # ka=None,
            dudr1=1 / 5.5,
            dudr2=1 / 5.5
            # dipole_direction=[1, 0, 0]
        )

        check_freqs = False
        if check_freqs:
            freqs = model.normal_modes()[0]
            raise Exception(freqs * UnitsData.convert("Hartrees", "Wavenumbers"))

        check_anh = False
        if check_anh:
            from Psience.VPT2 import VPTRunner

            VPTRunner.run_simple(
                [mol.atoms, mol.coords],
                potential_derivatives=model.potential(mol.coords, deriv_order=4)[1:],
                dipole_derivatives=model.dipole(mol.coords, deriv_order=3),
                order=2, states=3,
                logger=True,
                degeneracy_specs='auto',
                calculate_intensities=True,
                include_coriolis_coupling=True
            )

            """
              0 0 0    4680.66314                   4611.38521
                              Harmonic                  Anharmonic
            State       Frequency    Intensity       Frequency    Intensity
              0 0 1    3896.87028     64.98650      3719.85792     62.04864
              0 1 0    3843.25804      0.17386      3676.15022      0.13738
              1 0 0    1621.19796     64.86522      1603.43661     64.14564
              0 0 2    7793.74057      0.00000      7405.54785      0.00156
              0 2 0    7686.51607      0.00000      7216.32913      0.00897
              2 0 0    3242.39591      0.00000      3197.00895      0.07403
              0 1 1    7740.12832      0.00000      7229.92435      1.31011
              1 0 1    5518.06824      0.00000      5308.87368      0.09319
              1 1 0    5464.45599      0.00000      5278.21350      0.07390
              0 0 3   11690.61085      0.00000     10968.42588      0.00035
              0 3 0   11529.77411      0.00000     10889.98015      0.00002
              3 0 0    4863.59387      0.00000      4780.71702      0.02966
              0 1 2   11636.99860      0.00000     10585.27120      0.00057
              1 0 2    9414.93853      0.00000      8988.16177      0.00045
              0 2 1   11583.38636      0.00000     10587.47926      0.05388
              2 0 1    7139.26620      0.00000      6888.02518      0.00327
              1 2 0    9307.71403      0.00000      8809.00009      0.00197
              2 1 0    7085.65395      0.00000      6870.41251      0.00035
              1 1 1    9361.32628      0.00000      8817.56679      0.00646
              """
            raise Exception(...)

        check_dvr = False
        if check_dvr:
            raise ValueError("you don't want to rerun this")
            dvr = model.setup_DVR(
                domain=[[1, 4], [1, 4], [np.deg2rad(60), np.deg2rad(160)]],
                divs=[800, 800, 800], po_divs=[25, 25, 25],
                potential_optimize=True,
                logger=True
            )
            po_data = dvr.run()
            """
            :: PotentialOptimizedDVR([WavefunctionBasisDVR(None, pts=20, pot=None), WavefunctionBasisDVR(None, pts=20, pot=None), WavefunctionBasisDVR(None, pts=20, pot=None)], pot=SympyExpr(0.203038951525208*(1 - exp(-0.813301570558368*sqrt(2)*(r[1,2] - 1.8253409520594)))**2 + 0.203038951525208*(1 - exp(-0.813301570558368*sqrt(2)*(r[2,3] - 1.82534095205941)))**2 + 0.255579575354735*(0.551593470847119*a[1,2,3] - 1)**2))
            :: g: [[SympyExpr(0.0005786177281533848), SympyExpr(3.42971451934982e-5*cos(a[1,2,3])), SympyExpr(-3.42971451934982e-5*sin(a[1,2,3])/r[2,3])], [SympyExpr(3.42971451934982e-5*cos(a[1,2,3])), SympyExpr(0.0005786177281533848), SympyExpr(0.0)], [SympyExpr(-3.42971451934982e-5*sin(a[1,2,3])/r[2,3]), SympyExpr(0.0), SympyExpr(0.000578617728153385/r[2,3]**2 - 6.85942903869965e-5*cos(a[1,2,3])/(r[1,2]*r[2,3]) + 0.000578617728153385/r[1,2]**2)]]
            :: mass: [None, None, None]
            :: g_deriv: [SympyExpr(0.0), SympyExpr(0.0), SympyExpr(6.85942903869965e-5*cos(a[1,2,3])/(r[1,2]*r[2,3]))]
            :: domain: [[1, 4], [1, 4], [1.0471975511965976, 2.792526803190927]]
            :: divs: [800, 800, 800]
            :: potential_function: SympyExpr(0.203038951525208*(1 - exp(-0.813301570558368*sqrt(2)*(r[1,2] - 1.8253409520594)))**2 + 0.203038951525208*(1 - exp(-0.813301570558368*sqrt(2)*(r[2,3] - 1.82534095205941)))**2 + 0.255579575354735*(0.551593470847119*a[1,2,3] - 1)**2)
            ::> constructing grid
            <::
            ::> constructing potential matrix
              > evaluating potential function over grid
            <::
            ::> constructing kinetic matrix
              > evaluating kinetic coupling
            <::
            ::> building Hamiltonian
            <::
            ::> evaluating wavefunctions
              ::> diagonalizing Hamiltonian
                > dimension=(8000, 8000)
                > density=9.750%
                > mode=None
              <::
            <::
            >>--------------------------------------------------<<
            ERROR

            ======================================================================
            ERROR: test_ModelPotentialAIMD (tests.DGBTests.DGBTests)
            ----------------------------------------------------------------------
            Traceback (most recent call last):
              File "/Users/Mark/Documents/UW/Research/Development/Peeves/Peeves/TestUtils.py", line 501, in Debug
                return fn(*args, **kwargs)
              File "/Users/Mark/Documents/UW/Research/Development/Psience/ci/tests/DGBTests.py", line 1040, in test_ModelPotentialAIMD
                raise Exception(po_data.wavefunctions.frequencies()*UnitsData.hartrees_to_wavenumbers)
            Exception: (
                6234.520135705428, 
                array([ 
                    1603.52651473,  3197.92783543,  
                    3677.14362091,  3719.74235567,  
                    4791.35107014,  5281.20835861,  5308.52076799,  
                    6402.00507514,  6876.02718339,  6885.85546794,  
                    7220.93261924,  7232.7182684 ,  7405.52232462,  
                    8049.27247436,  8459.98906896,  8479.50398517,  
                    8821.70119285,  8825.8627132 ,  8988.5644022 ,  
                    9731.71367219, 10046.96812294, 10127.39430738, 
                    10403.06169691, 10407.06227839
                    ]))
             """
            raise Exception(
                po_data.wavefunctions.energies[1] * UnitsData.hartrees_to_wavenumbers,
                po_data.wavefunctions.frequencies() * UnitsData.hartrees_to_wavenumbers
            )

        # mol.potential_derivatives = model.potential(mol.coords, deriv_order=2)[1:]
        # raise Exception(mol.coords, mol.normal_modes.modes)

        bend, symm, asymm = model.normal_modes()[0]
        init_e = np.array([
                [-bend,     0,       0 ],
                [ bend,     0,       0 ],
                [    0,  -symm,       0 ],
                # [    0,   symm,       0 ],
                [    0,     0,   asymm ],
                # [ bend,  symm,       0 ],
                # [-bend,  symm,       0 ],
                # [ bend,     0,   asymm ],
                # [-bend,     0,   asymm ],
                # [    0,  symm,   asymm ],
                # [ bend,  symm,   asymm ],
                # [-bend,  symm,   asymm ],
            ])
        sim = model.setup_AIMD(
            initial_energies=init_e * 1,
            timestep=25,
            track_velocities=True
        )
        sim.propagate(100)


        # coords = sim.extract_trajectory(flatten=True, embed=mol.coords)
        coords, velocities = sim.extract_trajectory(flatten=True, embed=mol.coords)
        coords = coords[len(init_e) - 1:]
        velocities = velocities[len(init_e) - 1:]
        momenta = 150 * velocities * mol.masses[np.newaxis, :, np.newaxis]
        virial_scaling = 1.2
        # all_pots = model.potential(coords)
        # virial_scaling = 1 + (all_pots - np.min(all_pots)) / (np.max(all_pots) - np.min(all_pots))
        # virial_scaling = virial_scaling[:, np.newaxis]
        pruning_energy = 3600 / UnitsData.hartrees_to_wavenumbers

        mol.potential_derivatives = model.potential(mol.coords, deriv_order=2)[1:]
        modes = mol.normal_modes.modes.basis.to_new_modes()
        emb_coords = modes.embed_coords(coords)
        new_emb = emb_coords * np.array([1, 1, -1])[np.newaxis, :]
        new_coords = modes.unembed_coords(new_emb)

        cartesians = False
        use_interpolation = True
        plot_interpolation_error = True
        use_quadrature = False
        permute_coords = True
        if use_interpolation:
            if permute_coords:
                interp_coords = np.concatenate([coords, new_coords], axis=0)
            else:
                interp_coords = coords
            potential_data = {
                'centers':interp_coords,
                'values':model.potential(interp_coords, deriv_order=2)
            }
        else:
            potential_data = None
        use_momenta = True
        use_pairwise = True

        if use_quadrature and use_interpolation:
            raise ValueError("don't use interpolation with quadrature...")

        # plot_wfns = {'modes':[2, 1]}
        plot_wfns = {'cartesians':[0, 1], 'num':10}

        with BlockProfiler(inactive=True):
            """
            >>------------------------- Running distributed Gaussian basis calculation -------------------------
            :: diagonalizing in the space of 15 S functions
            :: ZPE: 4569.059696337776
            :: Frequencies: [ 1608.80754098  3216.94449611  3703.47869608  3762.41780866  4835.82080583  5420.2895992   5433.19405061  7250.01177644  7323.31350762  7501.31168876  7617.50706198  7893.4242054   9746.29740443 10228.83473985]
            >>--------------------------------------------------<<
            """

            """
            >>------------------------- Running distributed Gaussian basis calculation -------------------------
            :: diagonalizing in the space of 56 S functions
            :: ZPE: 4581.861159506666
            :: Frequencies: [1593.42440048 3179.02517045 3679.80951784 3711.92403488 4757.05374131 5285.49196186 5302.43012877 6321.90439255 6872.37334386 6893.83413192 7232.93199651 7311.54294011 7406.72275565 7890.02544108 8457.85459221 8855.08195512 8879.86908611 8979.08988067 9243.53666465 9448.79285329]
            >>--------------------------------------------------<<
            """

            """
            :: diagonalizing in the space of 6 S functions
            :: ZPE: 4621.103254514144
            :: Frequencies: [1626.13229446 1854.91547825 3542.66011746 3814.91511096 5645.33283672]
            """

            """
            :: ZPE: 4597.450569703509
            :: Frequencies: [ 1532.0238538   3303.70601926  3886.06002078  3964.43747624  4728.98993909  5346.90145145  5693.6426873   6405.30581875  7154.85197981  7398.30642509  8923.81178447  9713.71335386 10610.13098991 11231.01060612]
            :: ZPE: 4592.909813102955
            :: Frequencies: [ 1638.0134135   3340.09636088  3581.51512274  3825.31857303  5495.38492396  5602.10440597  5937.30478358  6935.0265671   7342.87842492  7409.26799591  7525.80528395  7636.42783572  8023.87015034  8660.93003546  9119.65216576  9505.35854727  9725.2442756  10216.54968368 10446.88673994 10589.9953759 ]
            """

            """
            :: ZPE: 4611.094983585687
            :: Frequencies: [1603.44863295 3198.13903952 3679.81865526 3725.47652423 4781.54053    5286.23312573 5321.31083824 6355.83780495 6882.2687276  6957.89008281 7235.86055485 7302.76906113 7416.36229189 7920.22057333 8467.07742527 8570.37531735 8836.64381987 8881.92981176 9006.90005472 9495.1200398 ]
            
            :: ZPE: 4571.150876676482
            :: Frequencies: [1604.93242976 3201.2762706  3679.48073856 3723.99293599 4788.08329001 5286.6499235  5319.70366332 6366.13763555 6890.23518753 6959.79446469 7233.18188159 7302.72415051 7413.80482603 7978.23252103 8491.44159641 8579.49532908 8837.49145321 8877.54093539 9004.81997599 9603.16119709]
            """

            dgb = model.setup_DGB(
                coords,
                potential_function=potential_data,
                # optimize_centers=False,
                optimize_centers=[
                    {
                        'method': 'energy-cutoff',
                        'cutoff': pruning_energy
                    } if pruning_energy is not None else None,
                    1e-14
                ],
                alphas={'method':'virial', 'scaling':virial_scaling},
                modes=None if cartesians else 'normal',
                cartesians=[0, 1] if cartesians else None,
                quadrature_degree=3
                , expansion_degree=2 if not use_quadrature else None
                , pairwise_potential_functions={
                    (0, 1):self.setupMorseFunction(model, 0, 1),
                    (0, 2):self.setupMorseFunction(model, 0, 2)
                } if use_pairwise and not use_quadrature else None
                , transformations='diag'
                , momenta=momenta if use_momenta else None
            )

            """
            >>------------------------- Running distributed Gaussian basis calculation -------------------------
            :: diagonalizing in the space of 29 S functions
            :: ZPE: 4621.070203025893
            :: Frequencies: [ 1612.5368789   3208.56213954  3696.97053809  3754.6146347   4812.80012189  5312.65709999  6411.84633085  6680.27922037  7257.77195695  7277.87948641  7377.19723447  7518.10407279  8797.66839764  9233.26327849  9925.21216708 10414.99363543 10625.47426277 10919.38280277 11015.37832137 12314.26015755]
            >>--------------------------------------------------<<
            """

            """
            >>------------------------- Running distributed Gaussian basis calculation -------------------------
            :: diagonalizing in the space of 28 S functions
            :: ZPE: 4575.775602600208
            :: Frequencies: [ 1613.07878096  3213.70373479  3696.06725792  3746.82351792  4812.61991409  5334.04205722  5380.7963598   6594.09044668  7153.21525645  7199.20742615  7359.9240031   7446.44199528  7466.5521983   8173.17755511  8946.95608661  9098.22416268  9688.77529614 10509.60876202 10650.76512812 11113.93952066]
            >>--------------------------------------------------<<
            """

            """
            :: diagonalizing in the space of 5 S functions
            :: ZPE: 4631.159993716219
            :: Frequencies: [1654.5274286  3399.10904144 3768.42654735 4351.43562365]
            
            >>------------------------- Running distributed Gaussian basis calculation -------------------------
            :: diagonalizing in the space of 50 S functions
            :: ZPE: 4581.899925357468
            :: Frequencies: [1594.40031422 3180.8893923  3680.26800032 3714.95577085 4757.77891927 5285.57205191 5306.77728956 6322.88067822 6885.56680522 6974.93397215 7236.77106831 7311.54476953 7407.40616614 7929.08822608 8497.51104151 8878.33658414 8942.86126327 9000.15199805 9520.71309612 9993.134724  ]
            >>--------------------------------------------------<<

            """

            """
            Quad (diag on top)
            :: Frequencies: [ 1608.15422901  3210.43539753  3687.80981187  3723.55153359  4895.08134347  5324.51974615  5349.55404715  6537.97195165  6996.55458136  7041.99024524  7290.16851734  7323.84521184  7463.58005119  8398.91185649  8786.60637452  8913.35832829  9112.61849297  9495.17548506  9503.36436839 10199.46614382]
            :: Frequencies: [1602.92862805 3200.09079743 3680.20944591 3722.43343526 4788.04990973 5290.4664062  5320.35773164 6431.87241655 6921.59675425 6925.25078426 7236.71298577 7242.86318064 7423.00166202 8029.54060322 8546.01311377 8551.19640146 8885.17965835 8933.58490114 9033.7236885  9748.41422037]
            
            Expansions (diag on top)
            :: Frequencies: [1602.72784175 3194.17457145 3678.12473949 3720.4027005  4777.70463931 5282.38796837 5309.00508594 6350.68799335 6873.27560286 6887.92343644 7224.47734662 7236.35889162 7408.21856099 7923.43823514 8458.76134409 8478.22811052 8827.51438024 8836.75664778 8995.88311051 9522.79437823]
            :: Frequencies: [1602.94342352 3200.22645253 3680.44751438 3721.98040484 4788.27296366 5290.77563873 5320.00898146 6431.92197999 6922.26482632 6925.11439887 7237.19784537 7242.6045962  7422.55836083 8029.53907426 8546.6466649  8551.0311153  8886.09236692 8933.92493092 9033.73227324 9749.00043004]
           """

            plot_gaussians = False
            if plot_gaussians:
                self.plot_gaussians(dgb, mol,
                                    cartesians=True,
                                    coordinate_sel=[0, 1]
                                    )
                raise Exception(...)

            if use_interpolation and plot_interpolation_error:
                self.plot_interpolation_error(dgb, model.potential)
                raise Exception(...)

            # print(dgb.gaussians.alphas[:3])
            # print(dgb.gaussians.coords.centers[:3])
            # print(dgb.gaussians.transformations[0][:3])
            # print(dgb.S[0, :3])
            # print(dgb.T[0, :3])
            # raise Exception(...)

            # type(self).default_num_plot_wfns = 5
            type(self).default_num_plot_wfns = 5
            self.runDGB(dgb, mol,
                        similarity_cutoff=.9,
                        # similarity_chunk_size=5,
                        # vmin=-.05,
                        # vmax=.05,
                        # domain=[[-2, 2], [-2, 2]],
                        # plot_wavefunctions=False,
                        plot_centers=True,
                        plot_spectrum=False,
                        plot_wavefunctions=plot_wfns,
                        # plot_wavefunctions={'cartesians':[0, 1]} if not cartesians else True
                        )
            plt.Graphics().show()

    @validationTest
    def test_ModelPotentialAIMD3DHOD(self):
        mol, model = self.buildWaterModel(
            # w2=None, wx2=None,
            # ka=None,
            w2=3869.47 * self.w2h / np.sqrt(2),
            wx2=84 * self.w2h / np.sqrt(2),
            atoms=['O', 'H', 'D'],
            dudr1=1 / 5.5,
            dudr2=1 / 5.5
            # dipole_direction=[1, 0, 0]
        )

        check_freqs = False
        if check_freqs:
            freqs = model.normal_modes()[0]
            raise Exception(freqs * UnitsData.convert("Hartrees", "Wavenumbers"))

        # from Psience.VPT2 import VPTRunner
        #
        # runner, _ = VPTRunner.construct(
        #     [mol.atoms, mol.coords],
        #     potential_derivatives=model.potential(mol.coords, deriv_order=4)[1:],
        #     order=2, states=3,
        #     logger=True,
        #     degeneracy_specs='auto'
        # )
        # raise Exception(
        #     runner.hamiltonian.modes.basis.matrix,
        #     runner.hamiltonian.coriolis_terms.base_terms.modes
        # )
        check_anh = False
        if check_anh:
            from Psience.VPT2 import VPTRunner

            VPTRunner.run_simple(
                [mol.atoms, mol.coords],
                potential_derivatives=model.potential(mol.coords, deriv_order=4)[1:],
                order=2, states=3,
                logger=True,
                degeneracy_specs='auto',
                calculate_intensities=False,
                include_coriolis_coupling=True
            )

            # model.run_VPT(order=2, states=3,
            #               logger=True,
            #               degeneracy_specs='auto'
            #               )
            """
            ZPE:       4013.07238                   3923.87672
            ============================================= IR Data ==============================================
            Initial State: 0 0 0 
                               Harmonic                  Anharmonic
            State       Frequency    Intensity       Frequency    Intensity
              0 0 1    3870.81584     33.92454      3695.37083     32.29979
              0 1 0    2737.23347     16.11050      2608.55203     14.12031
              1 0 0    1418.09545     49.64519      1402.64554     49.13414
              0 0 2    7741.63168      0.00000      7222.94025      0.65500
              0 2 0    5474.46694      0.00000      5106.15439      0.34140
              2 0 0    2836.19091      0.00000      2805.76221      1.55285
              0 1 1    6608.04931      0.00000      6301.78208      0.00245
              1 0 1    5288.91129      0.00000      5084.86993      0.06913
              1 1 0    4155.32893      0.00000      4001.03724      0.08449
              2 1 0    5573.42438      0.00000      5369.60935      0.04358
              3 0 0    4254.28636      0.00000      4200.81468      0.00838
            ====================================================================================================
            """
            raise Exception(...)

        # mol.potential_derivatives = model.potential(mol.coords, deriv_order=2)[1:]
        # raise Exception(mol.coords, mol.normal_modes.modes)

        sim = model.setup_AIMD(
            initial_energies=np.array([
                [5000, 7000, 0],
                [5000, -7000, 0],
                [-5000, 7000, 0],
                [-5000, -7000, 0],
                [5000, 0, 7000],
                [5000, 0, -7000],
                [-5000, 0, 7000],
                [-5000, 0, -7000],
                # [2000 * self.w2h, 0],
                # [0, 2000 * self.w2h],
                # [-2000 * self.w2h, 0],
                # [0, -2000 * self.w2h],
                # [-15000 * self.w2h, -15000 * self.w2h],
                # [15000 * self.w2h, -15000 * self.w2h],
                # [10000 * self.w2h, 0],
                # [0, 10000 * self.w2h],
                # [-10000 * self.w2h, 0],
                # [0, -10000 * self.w2h]
            ]) * self.w2h * .8,
            timestep=10
        )
        sim.propagate(3)
        coords = sim.extract_trajectory(flatten=True, embed=mol.coords)

        cartesians = False
        with BlockProfiler(inactive=True):

            dgb = model.setup_DGB(
                np.round(coords, 8),
                # optimize_centers=False,
                optimize_centers={
                    'method': 'gram-schmidt',
                    'norm_cutoff': 1e-14,
                    'allow_pivoting': True
                },
                # alphas=[1, 2, 3],
                # optimize_centers=False,
                modes=None if cartesians else 'normal',
                cartesians=[0, 1] if cartesians else None,
                quadrature_degree=3,
                # expansion_degree=2,
                # pairwise_potential_functions={
                #     (0, 1): self.setupMorseFunction(model, 0, 1),
                #     (0, 2): self.setupMorseFunction(model, 0, 2,
                #                                     w=3869.47 * self.w2h / np.sqrt(2),
                #                                     wx=84 * self.w2h / np.sqrt(2),
                #                                     )
                # }
            )

            # print(dgb.gaussians.coords.centers[:3])
            # print(dgb.gaussians.alphas[:3])

            type(self).default_num_plot_wfns = 5
            self.runDGB(dgb, mol,
                        # similarity_chunk_size=5,
                        # vmin=-.05,
                        # vmax=.05,
                        # domain=[[-2, 2], [-2, 2]],
                        # plot_wavefunctions=False,
                        # mode='classic',
                        # subspace_size=15,
                        plot_wavefunctions={'cartesians': [0, 1]} if not cartesians else True
                        )

    @classmethod
    def declusterPoints(cls, points, radius):
        pivots = np.arange(len(points))
        dec_pts = points.reshape(points.shape[0], -1)
        for i in range(len(points)):
            cur_pos = pivots[i]
            dists = np.linalg.norm(
                dec_pts[cur_pos][np.newaxis, :] - dec_pts[pivots[i+1:], :],
                axis=1
            )
            good_pos = np.where(dists > radius)
            if len(good_pos) == 0 or len(good_pos[0]) == 0:
                break
            pivots = np.concatenate([pivots[:i+1], pivots[i+1:][good_pos]])
        return points[pivots]
    @classmethod
    def getMBPolModel(cls, atoms=None, ref=None, embed=True):
        loader = ModuleLoader(TestManager.current_manager().test_data_dir)
        mbpol = loader.load("LegacyMBPol").MBPol

        b2a = UnitsData.convert("BohrRadius", "Angstroms")

        def potential(coords, deriv_order=0, chunk_size=int(5e5)):

            coords = coords.reshape(-1, 9)

            just_vals = deriv_order is None
            if just_vals:
                deriv_order = 0

            chunks = [[] for _ in range(deriv_order + 1)]
            num_chunks = int(len(coords) / chunk_size) + 1

            for coords in np.array_split(coords, num_chunks):
                # interp.logger.log_print("evaluating energies")
                energies = mbpol.get_pot(coords=coords.reshape(-1, 3, 3) * b2a, nwaters=1,
                                         threading_vars=['energy', 'coords'], threading_mode='omp')
                if deriv_order > 0:
                    derivs = []
                    grads = lambda c: mbpol.get_pot_grad(
                        nwaters=1, coords=c.reshape(-1, 3, 3) * b2a, threading_vars=['energy', 'grad', 'coords'],
                        threading_mode='omp'
                    )['grad'].reshape(c.shape) * b2a
                    # interp.logger.log_print("evaluating forces")
                    derivs.append(grads(coords))
                    if deriv_order > 1:
                        # interp.logger.log_print("evaluating Hessians")
                        hess_fun = FiniteDifferenceDerivative(
                            grads,
                            function_shape=(9, 9)
                        )

                        # chunks = []
                        # num_chunks = int(len(coords)/1000)
                        # for a in np.array_split(coords, num_chunks):
                        #     chunks.append(hess_fun.derivatives(a).compute_derivatives(1))
                        # hess = np.concatenate(chunks, axis=0)
                        new_derivs = hess_fun.derivatives(coords).derivative_tensor(list(range(1, deriv_order)))
                        # print([d.shape for d in new_derivs])
                        derivs.extend(
                            np.moveaxis(d, -2, 0)
                            for i, d in enumerate(new_derivs)
                        )
                    # interp.logger.log_print("done")
                    for i, d in enumerate([energies] + derivs):
                        chunks[i].append(d)
                else:
                    chunks[0].append(energies)

            for i, c in enumerate(chunks):
                chunks[i] = np.concatenate(c, axis=0)

            if just_vals:
                chunks = chunks[0]

            return chunks

        if ref is None:
            ref = np.array([
                [0.00000000e+00, 6.56215885e-02, 0.00000000e+00],
                [7.57391014e-01, -5.20731105e-01, 0.00000000e+00],
                [-7.57391014e-01, -5.20731105e-01, 0.00000000e+00]
            ]) * UnitsData.convert("Angstroms", "BohrRadius")

        if atoms is None:
            atoms = ["O", "H", "H"]

        ref_mol = Molecule(
            atoms,
            ref
        )
        if embed:
            ref_mol = ref_mol.get_embedded_molecule(load_properties=False)

        return potential, ref_mol
    @debugTest
    def test_WaterAIMD(self):
        pot, mol = self.getMBPolModel()

        get_anharmonicity = False
        if get_anharmonicity:
            mbpol_deriv = pot(mol.coords, deriv_order=4)
            mol.internals = [[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]]
            conv = mol.get_cartesians_by_internals(4, strip_embedding=True)
            from McUtils.Zachary import TensorDerivativeConverter
            internals = TensorDerivativeConverter(conv, [d.squeeze() for d in mbpol_deriv[1:]]).convert()
            g = (1 / mol.atomic_masses[0] + 1 / mol.atomic_masses[1])
            w = np.sqrt(internals[1][0, 0] * g)
            wx = (g/(4*w)) ** 2 * (internals[3][0, 0, 0, 0] - 5 / 3 * internals[2][0, 0, 0]/internals[1][0, 0])
            raise Exception(
                w * UnitsData.hartrees_to_wavenumbers, # 3891.634745263701
                wx * UnitsData.hartrees_to_wavenumbers # 185.31310314514627
            )

        check_freqs = False
        if check_freqs:
            mol.potential_derivatives = pot(mol.coords, deriv_order=2)[1:]
            freqs = mol.normal_modes.modes.freqs
            raise Exception(freqs * UnitsData.convert("Hartrees", "Wavenumbers"))

        # mol.potential_derivatives = pot(mol.coords, deriv_order=2)[1:]
        # mol.internals = [[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]]
        # grid = np.linspace(-.5, .5, 15)
        # modes = mol.normal_modes.modes.basis
        # disps = grid[:, np.newaxis, np.newaxis] * modes.matrix[:, 0].reshape(mol.coords.shape)[np.newaxis]
        # coords = mol.coords[np.newaxis] + disps
        # woof = (coords - modes.origin[np.newaxis]).reshape(len(coords), -1)  @ modes.inverse.T
        # raise Exception(woof)
        # import McUtils.Plots as plt
        # plt.Plot(
        #     np.array(pot(coords)[0]) - np.array(pot(mol.coords)[0]),
        #     (mol.normal_modes.modes.freqs[0]**2)*(grid**2) / 2
        # ).show()
        # raise Exception(...)

        check_anh = False
        if check_anh:
            from Psience.VPT2 import VPTRunner

            VPTRunner.run_simple(
                [mol.atoms, mol.coords],
                potential_derivatives=[x.reshape((9,)*(i+1)) for i,x in enumerate(pot(mol.coords, deriv_order=4)[1:])],
                order=2, states=2,
                logger=True,
                degeneracy_specs='auto',
                calculate_intensities=False,
                include_coriolis_coupling=True
            )
            """
            State     Harmonic   Anharmonic     Harmonic   Anharmonic
                           ZPE          ZPE    Frequency    Frequency
            0 0 0   4713.07975   4636.90793            -            - 
            0 0 1            -            -   3944.32593   3753.07778 
            0 1 0            -            -   3832.76457   3654.53255 
            1 0 0            -            -   1649.06900   1594.42231 
            0 0 2            -            -   7888.65187   7438.71221 
            0 2 0            -            -   7665.52914   7192.33980 
            2 0 0            -            -   3298.13800   3152.62843 
            0 1 1            -            -   7777.09050   7240.72933 
            1 0 1            -            -   5593.39493   5326.66828 
            1 1 0            -            -   5481.83357   5232.92568 
            """
            raise Exception(...)

        mol.potential_derivatives = pot(mol.coords, deriv_order=2)[1:]
        bend, symm, asym = mol.normal_modes.modes.freqs

        base_dir = os.path.expanduser('~/Documents/Postdoc/AIMD-Spec/')
        for steps in [250]:#[10, 25, 50, 100, 150]:
            os.makedirs(base_dir, exist_ok=True)
            timestep = 15
            ntraj = 50

            energy_scaling = .8
            seed = 12213123

            svd_contrib_cutoff = None #1.5*1e-2
            svd_min_value = None #1e-4
            max_overlap_cutoff = 1 - 1e-4
            gs_cutoff = 14
            declustering_radius = None# 1e-4

            sim_cutoff = .9
            plot_similarity = False
            # pruning_energy = 1000
            pruning_energy = None
            pruning_probabilities = [[250, None], [500, .8], [700, .5], [1000, .25]]
            if pruning_probabilities is not None:
                pruning_probabilities = [
                    [e / UnitsData.hartrees_to_wavenumbers, p]
                    for e,p in pruning_probabilities
                ]

            use_interpolation = True
            plot_interpolation_error = False

            use_momenta = True
            momentum_scaling = 1000

            use_pairwise_potentials = False
            morse_w = 3891.634745263701 * self.w2h
            morse_wx = 185.31310314514627 * self.w2h
            expansion_degree = 2
            virial_scaling = 1
            transformation_method = 'diag' #4642.62343373526, 1624.36701583 3196.27137197 3474.1645502  3740.31164635
            # transformation_method = 'rpath' #4597.724422233999, 1613.48466441 3197.30480264 3559.70550684 3858.89799206
            # transformation_method = None #1698.2354715  3295.17558598 3667.70526921 3929.12204053

            run_profiler = True
            plot_gaussians = False
            plot_wfns = not run_profiler
            # plot_wfns = {'modes':[2, 1]}
            # plot_wfns = {'cartesians':[0, 1], 'num':15}

            plot_dir = os.path.join(base_dir, 'Figures', 'mbpol')
            if not plot_wfns:
                plot_dir = None
            if plot_dir is not None:
                all_params = (
                    timestep,
                    ntraj,
                    energy_scaling,

                    seed,

                    max_overlap_cutoff ,
                    gs_cutoff,
                    declustering_radius,
                    sim_cutoff,
                    pruning_energy,
                    None if not use_pairwise_potentials else morse_wx,
                    None if not use_pairwise_potentials else morse_w
                    )
                param_hash = str(hash(all_params))[:25]
                plot_dir = os.path.join(plot_dir,
                                        '{}_E{}_s{}_nt{}_ts{}_sim{}_h{}'.format(
                                            'interp' if use_interpolation else 'exact',
                                            pruning_energy,
                                            steps,
                                            ntraj,
                                            timestep,
                                            int(sim_cutoff * 100),
                                            param_hash
                                        )
                                        )

            if pruning_energy is not None:
                pruning_energy = pruning_energy / UnitsData.hartrees_to_wavenumbers

            pot_cmap = 'Greys_r'
            pot_points = 250
            wfn_cmap = 'coolwarm'
            wnf_points = 100
            wfn_contours = 10
            plot_centers = False
            plot_styles = {
                'frame': True,
                # 'ticks_style': (False, False),  # , dict(left=False, right=False))
                'ticks': ([-2, -1, 0, 1, 2], [-1, -.5, 0]),
                'axes_labels': ['x (bohr)', 'y (bohr)'],
                'padding': ([60, 0], [40, 0]),
                'aspect_ratio': 'auto'
            }
            plot_styles = dict(
                plot_styles,
                plot_range=[[-2.5, 2.5], [-1.5, 0.5]],
                image_size=[500, 500 * 2/5]
            )

            # initial_energies = np.array([
            #                 [ 3.0, 0.0, 0.0],
            #                 [ 1.0, 0.0, 0.0],
            #                 [ 0.0, 3.0, 0.0],
            #                 [ 0.0, 0.0, -3.0],
            #                 [ 2.0, 0.0, -2.0],
            #                 [ 0.0, 2.0, -2.0],
            #                 [-2.0, 2.0, 2.0],
            #                 # [ 0.7, 0.7, 0.0],
            #                 # [-0.7, 0.7, 0.0],
            #             ]) * np.array([bend, symm, asym])[np.newaxis, :]
            np.random.seed(seed)
            dirs = np.random.normal(0, 1, size=(ntraj, 3))
            dirs = dirs / np.linalg.norm(dirs, axis=1)[:, np.newaxis]
            dirs = np.concatenate([
                [
                    [1, 0, 0],
                    [0, 1, 0],
                    [0, 0, 1]
                ],
                dirs
            ])
            initial_energies = energy_scaling * dirs * np.array([bend, symm, asym])[np.newaxis, :]
            # initial_energies = initial_energies / np.linalg.norm(initial_energies, axis=1)[:, np.newaxis]
            # initial_energies = initial_energies * np.array([bend, symm, asym])[np.newaxis, :]
            # ninit = len(initial_energies)

            sim = mol.setup_AIMD(
                pot,
                initial_energies=np.array(initial_energies) * energy_scaling,
                timestep=timestep  # DON'T MESS WITH THIS
                , track_velocities=True
            )
            sim.propagate(steps)
            coords, velocities = sim.extract_trajectory(flatten=True, embed=mol.coords)
            momenta = momentum_scaling * velocities * mol.masses[np.newaxis, :, np.newaxis]

            # plt.HistogramPlot(pot(coords)[0]*219475.6).show()
            # raise Exception(...)

            # plops = pot(coords)[0] * 219475
            # coords = coords[plops < 7500]


            # os.makedirs(run_dir, exist_ok=True)
            # if save_plots:
            #     plot_dir=run_dir
            #     os.makedirs(plot_dir, exist_ok=True)

            print("="*25, "Steps:", steps, "Points:", len(coords), "="*25)

            cartesians = False
            if use_interpolation:
                augment = False
                if augment:
                    with Checkpointer.from_file(os.path.join(base_dir, interp_traj_file)) as chk:
                        traj = np.concatenate([
                            chk['coords'],
                            coords
                        ])
                        pot_vals = pot(traj, deriv_order=2)
                        interp_data = {'centers':traj, 'values':pot_vals}
                else:
                    interp_data = {'centers':coords, 'values':pot(coords, deriv_order=2)}
            else:
                interp_data = None

            with BlockProfiler(inactive=not run_profiler):
                # crd = np.round(coords[len(initial_energies)-1:], 8)
                # if declustering_radius is not None:
                #     crd = self.declusterPoints(crd, declustering_radius)
                print("input points:", len(coords))

                dgb = DGB.construct(
                    coords,
                    pot if not use_interpolation else interp_data,
                    masses=mol.atomic_masses,
                    alphas={'method':'virial', 'scaling':virial_scaling},
                    transformations=transformation_method,
                    optimize_centers=[
                        x for x in
                        [
                            {
                                'method': 'energy-cutoff',
                                'probabilities':pruning_probabilities,
                                'cutoff': pruning_energy
                            } if (pruning_energy is not None or pruning_probabilities is not None) else None,
                            {
                                'method': 'decluster',
                                'cluster_radius': declustering_radius
                            } if declustering_radius is not None else None,
                            {
                                'method': 'svd',
                                'contrib_cutoff': svd_contrib_cutoff,
                                'min_value': svd_min_value
                            } if svd_contrib_cutoff is not None else None,
                            {
                                'method': 'gram-schmidt',
                                'norm_cutoff': 10 ** (-gs_cutoff),
                                'max_overlap_cutoff': max_overlap_cutoff
                            } if gs_cutoff is not None else None,
                        ]
                        if x is not None
                    ],
                    modes=None if cartesians else 'normal',
                    cartesians=[0, 1] if cartesians else None,
                    quadrature_degree=3,
                    expansion_degree=expansion_degree,
                    pairwise_potential_functions={
                        (0, 1): self.setupMorseFunction(
                            mol.atomic_masses[0],
                            mol.atomic_masses[1],
                            np.linalg.norm(mol.coords[0] - mol.coords[1]),
                            w=morse_w,
                            wx=morse_wx
                            ),
                        (0, 2): self.setupMorseFunction(
                            mol.atomic_masses[0],
                            mol.atomic_masses[2],
                            np.linalg.norm(mol.coords[0] - mol.coords[2]),
                            w=morse_w,
                            wx=morse_wx
                        )
                    } if use_pairwise_potentials else None,
                    momenta=momenta if use_momenta else None,
                    logger=True
                )


                # plt.Graphics().show()
                # raise Exception(...)

                if use_interpolation and plot_interpolation_error:
                    self.plot_interpolation_error(dgb, pot)
                    raise Exception(...)

                if plot_gaussians:
                    self.plot_gaussians(dgb, mol, 10, plot_dir=plot_dir)
                    raise Exception(...)

                type(self).default_num_plot_wfns = 5
                if plot_wfns is True:
                    plot_wfns = {'cartesians': [0, 1]} if not cartesians else True
                wfns = self.runDGB(dgb, mol,
                                   mode='similarity',
                                   similarity_cutoff=sim_cutoff,
                                   plot_similarity=plot_similarity,
                                   plot_wavefunctions=plot_wfns,
                                   plot_spectrum=False,
                                   pot_cmap=pot_cmap,
                                   pot_points=pot_points,
                                   wfn_cmap=wfn_cmap,
                                   wfn_points=wnf_points,
                                   wfn_contours=wfn_contours,
                                   plot_centers=plot_centers,
                                   plot_dir=plot_dir,
                                   **plot_styles
                                   )

                if plot_dir is not None:
                    with open(os.path.join(plot_dir, 'freqs.txt'), 'w+') as woof:
                        print("ZPE:", wfns.energies[0] * UnitsData.hartrees_to_wavenumbers, file=woof)
                        print("Freqs:", wfns.frequencies() * UnitsData.hartrees_to_wavenumbers, file=woof)


            if plot_dir is not None or plot_wfns is False:
                plt.Graphics().show()

    @validationTest
    def test_WaterFromGauss(self):

        check_anh = False
        if check_anh:
            from Psience.VPT2 import VPTRunner
            VPTRunner.run_simple(TestManager.test_data('h2o_aimd_opt.fchk'),
                                 states=2, order=2,
                                 degeneracy_specs='auto'
                                 )
            """
            0 0 0      4671.11065                   4598.17667
            State       Frequency    Intensity       Frequency    Intensity
              0 0 1    3922.23915     56.86267      3734.52710     50.58253
              0 1 0    3817.17658      9.21138      3640.75276      6.13469
              1 0 0    1602.80557     66.70198      1554.84636     69.96507
              0 0 2    7844.47830      0.00000      7403.75133      0.00786
              0 2 0    7634.35317      0.00000      7164.50137      0.87675
              2 0 0    3205.61114      0.00000      3078.10542      2.38482
              0 1 1    7739.41573      0.00000      7209.82293      3.76134
              1 0 1    5525.04472      0.00000      5271.68806      4.38346
              1 1 0    5419.98215      0.00000      5180.54068      0.06344
              """
            raise Exception(...)

        mol = Molecule.from_file(TestManager.test_data('h2o_aimd_opt.fchk'))
        mol.potential_derivatives; mol.dipole_derivatives; # kludge
        mol = mol.get_embedded_molecule(embed_properties=True)

        get_initial_velocities = False
        if get_initial_velocities:
            nms = mol.normal_modes.modes.basis
            bend, sym, asym = nms.freqs
            initial_velocities = AIMDSimulator.mode_energies_to_velocities(
                nms.inverse.T,
                mol.atomic_masses,
                [
                    [bend, sym, asym]
                ],
                inverse=nms.matrix.T
            )
            with np.printoptions(suppress=True):
                print()
                print(initial_velocities)
            raise Exception(...)


        with GaussianFChkReader(TestManager.test_data('h2o_aimd_opt.fchk')) as parser:
            ref_eng = parser.parse(['Total Energy'])['Total Energy']

        coords = None
        grads = None
        hess = None
        engs = None
        # traj_file = TestManager.test_data('h2o_aimd.log')
        for traj_file in [
            os.path.expanduser('~/Documents/Postdoc/AIMD-Spec/bomd_h2o_rand.log'),
            os.path.expanduser('~/Documents/Postdoc/AIMD-Spec/bomd_h2o_rand_5.log'),
            # os.path.expanduser('~/Documents/Postdoc/AIMD-Spec/bomd_h2o_dupes.log'),
            os.path.expanduser('~/Documents/Postdoc/AIMD-Spec/bomd_h2o_long.log'),
            # os.path.expanduser('~/Documents/Postdoc/AIMD-Spec/bomd_h2o_init.log'),
            # os.path.expanduser('~/Documents/Postdoc/AIMD-Spec/bomd_h2o_rand_1.log')
        ]:
            with GaussianLogReader(traj_file) as parser:
                traj_data = parser.parse(['AIMDTrajectory'])['AIMDTrajectory']

            c = traj_data.coords
            g = traj_data.vals.gradients#[:3]
            h = traj_data.vals.hessians#[:3]
            e = traj_data.vals.energies - ref_eng


            # # raise Exception(e[:5]*219475, e[100:105]*219475)
            #
            # dupe_pos = np.where(np.diff(e) == 0)[0] + 1
            # print(dupe_pos)
            # good_pos = np.setdiff1d(np.arange(len(e)), dupe_pos)
            # e = e[good_pos]
            # g = g[good_pos]
            # h = h[good_pos]
            #
            # print(c.shape, g.shape, h.shape, e.shape)
            #
            # c = c[:len(e)]

            coords = c if coords is None else np.concatenate([coords, c])
            engs = e if engs is None else np.concatenate([engs, e])
            grads = g if grads is None else np.concatenate([grads, g])
            hess = h if hess is None else np.concatenate([hess, h])


        # nvals = min(len(coords), len(grads))
        # coords = coords[:nvals]
        # grads = grads[:nvals]
        # hess = hess[:nvals]

        pot, _ = self.getMBPolModel()
        new_engs, new_grads, new_hess = pot(coords, deriv_order=2)
        # engs, grads, hess = new_engs, new_grads, new_hess

        # raise Exception(
        #     # len(engs), # 135
        #     np.max(np.abs(engs - new_engs) * 219475),
        #     engs[:10] * 219475,
        #     new_engs[:10] * 219475,
        #     np.where(np.abs(engs - new_engs) * 219475 > 1000)
        # )

        embedding = mol.get_embedding_data(coords)

        rots = embedding.rotations
        # need to reembed derivatives
        coords = nput.vec_tensordot( #...?
            embedding.coord_data.coords,
            rots,
            shared=1,
            axes=[-1, 2]
        )
        grads = nput.vec_tensordot(
            grads.reshape(grads.shape[0], -1, 3),
            embedding.rotations,
            shared=1,
            axes=[-1, 1] # dv/dx has inverse embedding
        ).reshape(grads.shape)
        hess = nput.vec_tensordot(
            hess.reshape((hess.shape[0], hess.shape[1], -1, 3)),
            embedding.rotations,
            shared=1,
            axes=[-1, 1]
        ).reshape(hess.shape)
        hess = np.moveaxis(hess, -2, -1)
        hess = nput.vec_tensordot(
            hess.reshape((hess.shape[0], hess.shape[1], -1, 3)),
            embedding.rotations,
            shared=1,
            axes=[-1, 1]
        ).reshape(hess.shape)

        coords = np.concatenate([[mol.coords], coords], axis=0)
        engs = np.concatenate([[0], engs], axis=0)
        grads = np.concatenate([[mol.potential_derivatives[0]], grads], axis=0)
        hess = np.concatenate([[mol.potential_derivatives[1]], hess], axis=0)
        values = [engs, grads, hess]

        # mw = np.diag(1/np.sqrt(np.broadcast_to(mol.atomic_masses[:, np.newaxis], (3, 3)).flatten()))
        # ugh = np.linalg.eigvalsh(mw@hess@mw)
        # raise Exception(np.sign(ugh) * np.sqrt(np.abs(ugh)) * 219475)

        interp_data = {'centers':coords, 'values':values}

        augment_trajectory = False
        if augment_trajectory:
            pot, _ = self.getMBPolModel()

            nms = mol.normal_modes.modes.basis
            bend, sym, asym = nms.freqs
            sim = mol.setup_AIMD(
                pot,
                initial_energies=np.array([
                    [ bend,  0, asym],
                    [-bend,  0, asym],
                    [ bend,  0,-asym],
                    [-bend,  0,-asym],
                    [ bend, sym,   0],
                    [-bend, sym,   0],
                    [ bend, sym,   0],
                    [-bend, sym,   0]
                ]) * .8,
                timestep=60
            )
            sim.propagate(20)
            sub_coords = sim.extract_trajectory(flatten=True, embed=mol.coords)
            coords = np.concatenate([coords, sub_coords])

        print("="*25, len(coords), "Points Sampled", "="*25)
        dgb = DGB.construct(
            coords,
            interp_data,
            masses=mol.atomic_masses,
            alphas='auto',
            transformations='diag',
            modes='normal',
            expansion_degree=2,
            optimize_centers={
                'method':'energy-cutoff',
                'cutoff':1000 / UnitsData.hartrees_to_wavenumbers
            },
            # pairwise_potential_functions={
            #     (0, 1): self.setupMorseFunction(
            #         mol.atomic_masses[0],
            #         mol.atomic_masses[1],
            #         np.linalg.norm(mol.coords[0] - mol.coords[1])
            #     ),
            #     (0, 2): self.setupMorseFunction(
            #         mol.atomic_masses[0],
            #         mol.atomic_masses[2],
            #         np.linalg.norm(mol.coords[0] - mol.coords[2])
            #     )
            # },
            logger=True
        )

        self.runDGB(dgb, mol,
                    similarity_cutoff=.9,
                    plot_centers=True,
                    # mode='classic',
                    # subspace_size=15,
                    plot_spectrum=False,
                    plot_wavefunctions={'cartesians': [0, 1]}
                    )