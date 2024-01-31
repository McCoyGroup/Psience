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
        if plot_cartesians:
            figure = mol.plot_molecule_function(
                cutoff_pot,
                axes=[0, 1],
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
    def plot_gaussians(cls, dgb, mol, nmax=10, start_at=0, sel=None):
        n = len(dgb.gaussians.coords.centers)
        wfns = DGBWavefunctions(
            np.zeros(n),
            np.eye(n),
            dgb
        )
        if sel is not None:
            subdgb_coords = dgb.gaussians.coords[:, sel]
            subpot = subdgb_coords.embed_function(dgb.pot.potential_function.og_fn)
        else:
            subdgb_coords = dgb.gaussians.coords  # [:, [0, 1]]
            subpot = dgb.pot.potential_function  # subdgb_coords.embed_function(dgb.pot.potential_function.og_fn)
        for i in range(start_at, min(n, nmax)):
            figure = cls.plot_dgb_potential(
                subdgb_coords, mol, subpot
            )
            wf = wfns[i]
            if sel is not None:
                wf = wf.project(sel)
            plot = wf.plot(
                plotter=plt.TriContourLinesPlot,
                domain_padding=2,
                cmap='RdBu',
                figure=figure,
                plot_centers={'color': 'red'}
            )
            plot.show()

        # wfns = wfns.as_cartesian_wavefunction()
        # wfns[12].plot_cartesians(
        #     [0, 1],
        #     contour_levels=16,
        #     # cmap='RdBu',
        #     plot_centers={'color': 'red'},
        #     domain_padding=.5,
        # ).show()

    default_num_plot_wfns = 5
    @classmethod
    def runDGB(cls,
               dgb: DGB,
               mol,
               plot_centers=True,
               plot_atoms=True,
               plot_potential=True,
               plot_wavefunctions=True,
               plot_spectrum=False,
               pot_cmap='viridis',
               pot_points=100,
               wfn_cmap='RdBu',
               wfn_points=100,
               wfn_contours=12,
               domain=None,
               domain_padding=1,
               potential_cutoff=15000,
               mode=None,
               nodeless_ground_state=None,
               min_singular_value=None,
               subspace_size=None,
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

            if plot_spectrum:
                spec[:5].plot().show()

            figure = None
            if isinstance(plot_wavefunctions, str) and plot_wavefunctions=='cartesians':
                plot_wavefunctions = {'cartesians':None}
            cartesian_plot_axes = None
            coordinate_sel = None
            if isinstance(plot_wavefunctions, dict):
                if 'cartesians' in plot_wavefunctions:
                    cartesian_plot_axes=plot_wavefunctions['cartesians']
                    plot_wavefunctions = plot_wavefunctions.get('num', True)
                    wfns = wfns.as_cartesian_wavefunction()
                    dgb = wfns.hamiltonian
                elif 'modes' in plot_wavefunctions:
                    coordinate_sel = plot_wavefunctions['modes']
                    plot_wavefunctions = plot_wavefunctions.get('num', True)
                    plot_potential = False
                else:
                    raise ValueError(plot_wavefunctions)
            if plot_wavefunctions is True:
                plot_wavefunctions = cls.default_num_plot_wfns
            figs = []
            if isinstance(plot_wavefunctions, int):
                plot_wavefunctions = list(range(plot_wavefunctions))
            if isinstance(plot_wavefunctions, list):

                pot = dgb.pot.potential_function
                if plot_potential:
                    pot_figure = cls.plot_dgb_potential(
                        dgb, mol, pot,
                        cmap=pot_cmap,
                        plot_points=pot_points,
                        domain=domain, domain_padding=domain_padding,
                        potential_cutoff=potential_cutoff,
                        plot_atoms=plot_atoms
                    )
                else:
                    pot_figure = None
                for i in plot_wavefunctions:
                    if i < len(wfns):

                        if pot_figure is not None:
                            figure = pot_figure.copy()

                        if isinstance(dgb.gaussians.coords, DGBCartesians):
                            figs.append(
                                wfns[i].plot_cartesians(
                                    cartesian_plot_axes,
                                    contour_levels=wfn_contours,
                                    cmap=wfn_cmap,
                                    plot_points=wfn_points,
                                    figure=figure,
                                    plot_centers={'color':'red'} if plot_centers else False,
                                    domain=domain,
                                    domain_padding=.5,
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
                                        contour_levels=wfn_contours,
                                        cmap=wfn_cmap,
                                        plot_points=wfn_points,
                                        plotter=plt.TriContourLinesPlot ,
                                        plot_centers={'color':'red'} if plot_centers else False,
                                        domain=domain,
                                        domain_padding=domain_padding,
                                        figure=figure,
                                        **plot_options
                                    )
                                )
                            else:
                                figs.append(
                                    wfns[i].plot(
                                        plot_centers={'color': 'red'} if plot_centers else False,
                                        domain=domain,
                                        cmap=wfn_cmap,
                                        plot_points=wfn_points,
                                        domain_padding=domain_padding,
                                        figure=figure,
                                        scaling=-.1,
                                        shift=wfns.energies[i],
                                        **plot_options
                                    )
                                )

                figs[0].show()

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
            m1 = AtomData["O", "Mass"] * UnitsData.convert("AtomicMassUnits", "AtomicUnitsOfMass")
        if m2 is None:
            m2 = AtomData["O", "Mass"] * UnitsData.convert("AtomicMassUnits", "AtomicUnitsOfMass")
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
                [5000 * self.w2h],
                [-5000 * self.w2h]
            ],
            timestep=10
        )
        sim.propagate(10)
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
                    (0, 1): self.setupMorseFunction(model, 0, 1)
                    # (0, 2): self.setupMorseFunction(model, 0, 2)
                }
            )

            self.runDGB(dgb, mol,
                        domain_padding=10,
                        plot_spectrum=True,
                        plot_wavefunctions=False#{'cartesians':[0, 1]} if not cartesians else True
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
            timestep=25
        )
        sim.propagate(35)
        coords = sim.extract_trajectory(flatten=True, embed=mol.coords)

        cartesians=False
        with BlockProfiler(inactive=True):

            dgb = model.setup_DGB(
                coords,
                optimize_centers=1e-14,
                # optimize_centers=False,
                modes=None if cartesians else 'normal',
                cartesians=[0, 1] if cartesians else None,
                # quadrature_degree=3,
                expansion_degree=2,
                pairwise_potential_functions={
                    (0, 1):self.setupMorseFunction(model, 0, 1),
                    (0, 2):self.setupMorseFunction(model, 0, 2)
                }
            )

            self.runDGB(dgb, mol,
                        similarity_cutoff=.95,
                        plot_spectrum=True,
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

    @validationTest
    def test_ModelPotentialAIMD3D(self):
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
        sim = model.setup_AIMD(
            initial_energies=np.array([
                [ bend,  symm,       0 ],
                [ bend, -symm,       0 ],
                [-bend,  symm,       0 ],
                [-bend, -symm,       0 ],
                [ bend,     0,   asymm ],
                [ bend,     0,  -asymm ],
                [-bend,     0,   asymm ],
                [-bend,     0,  -asymm ],
            ]) * 1,
            timestep=25
        )
        sim.propagate(50)
        coords = sim.extract_trajectory(flatten=True, embed=mol.coords)

        cartesians = False
        with BlockProfiler(inactive=True):

            dgb = model.setup_DGB(
                coords,
                # optimize_centers=False,
                optimize_centers=1e-14,
                # optimize_centers=False,
                modes=None if cartesians else 'normal',
                cartesians=[0, 1] if cartesians else None,
                quadrature_degree=3
                , expansion_degree=2
                , pairwise_potential_functions={
                    (0, 1):self.setupMorseFunction(model, 0, 1),
                    (0, 2):self.setupMorseFunction(model, 0, 2)
                }
                , transformations='diag'
            )
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
                self.plot_gaussians(dgb, mol, 10, sel=[0, 1])
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
                        # similarity_chunk_size=5,
                        # vmin=-.05,
                        # vmax=.05,
                        # domain=[[-2, 2], [-2, 2]],
                        # plot_wavefunctions=False,
                        plot_centers=False,
                        plot_spectrum=True,
                        # plot_wavefunctions=False,
                        plot_wavefunctions={'cartesians':[0, 1]} if not cartesians else True
                        )

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
                    'overlap_cutoff': 1e-14,
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
    @validationTest
    def test_WaterAIMD(self):
        pot, mol = self.getMBPolModel()

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

        plot_dir = None
        save_plots = False
        for steps in [100]:#[10, 25, 50, 100, 150]:
            base_dir = os.path.expanduser('~/Documents/Postdoc/AIMD-Spec/water_new/')
            os.makedirs(base_dir, exist_ok=True)
            timestep = 15
            ntraj = 50

            energy_scaling = 1
            gs_cutoff = 14
            pruning_energy = 1000 / UnitsData.hartrees_to_wavenumbers
            # pruning_energy = None
            # bend = bend / energy_scaling
            use_interpolation = True
            plot_interpolation_error = True

            """
            >>------------------------- Running distributed Gaussian basis calculation -------------------------
            :: ZPE: 4579.793892345433
            :: Frequencies: [1595.35033389 3155.18176385 3656.9448941  3754.82924308 4672.58609003 5236.81240866 5332.21066675 6142.14536336 6780.07394829 6874.65167784 7200.61819466 7247.70256692 7442.44249515 7555.52756485 8283.55505483 8377.95782511 8765.05498695 8808.5239685  8897.08668035 9001.93290534]
            >>--------------------------------------------------<<
            
            :: diagonalizing in the space of 12 S functions
            :: ZPE: 4601.454181757252
            :: Frequencies: [1621.09465946 3175.19070097 3687.94757724 3790.43474694 5174.30285754 5435.00507871 5722.86447725 7112.86489917 7538.33457983 7735.17129795 7937.37086245]
            
            >>------------------------- Running distributed Gaussian basis calculation -------------------------
            :: diagonalizing in the space of 12 S functions
            :: ZPE: 4602.475402196309
            :: Frequencies: [1617.49559128 3145.73144444 3683.56942544 3766.45797359 5113.10993611 5333.96662647 5911.9881934  7073.52852616 7521.64746949 7962.9483711  8091.36424309]
            >>--------------------------------------------------<<
            
            No Pruning
            >>------------------------- Running distributed Gaussian basis calculation -------------------------
            :: diagonalizing in the space of 12 S functions
            :: ZPE: 4602.475401562745
            :: Frequencies: [1617.49559106 3145.731444   3683.56942494 3766.45797307 5113.1099354  5333.96662574 5911.98819259 7073.52852519 7521.64746846 7962.94837    8091.36424197]
            >>--------------------------------------------------<<
            """

            sim_cutoff = .8
            use_pairwise_potentials = False
            expansion_degree = 2
            virial_scaling = 1/2

            transformation_method = 'diag' #4642.62343373526, 1624.36701583 3196.27137197 3474.1645502  3740.31164635
            # transformation_method = 'rpath' #4597.724422233999, 1613.48466441 3197.30480264 3559.70550684 3858.89799206
            # transformation_method = None #1698.2354715  3295.17558598 3667.70526921 3929.12204053

            plot_gaussians = False
            plot_wfns = False
            # plot_wfns = {'modes':[2, 1]}
            plot_wfns = {'cartesians':[0, 1], 'num':5}

            pot_cmap = 'Greys_r'
            pot_points = 250
            wfn_cmap = 'coolwarm'
            wnf_points = 100
            wfn_contours = 10
            plot_centers = False

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
            np.random.seed(12213123)
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
            ninit = len(initial_energies)


            interp_traj_file='traj_50000_10_1_12.hdf5'
            augment = False
            if steps > 1000:
                with Checkpointer.from_file(os.path.join(base_dir, f'traj_{steps}_{timestep}_{energy_scaling}_{ninit}.hdf5')) as chk:
                    try:
                        coords = chk['coords']
                    except Exception as e:
                        sim = mol.setup_AIMD(
                            pot,
                            initial_energies=np.array(initial_energies) * energy_scaling,
                            timestep=timestep#DON'T MESS WITH THIS
                        )
                        sim.propagate(steps)
                        coords = sim.extract_trajectory(flatten=True, embed=mol.coords)
                        chk['coords'] = coords

                plt.Plot(np.arange(10), np.arange(10)).show()
                raise Exception(...)
            else:
                sim = mol.setup_AIMD(
                    pot,
                    initial_energies=np.array(initial_energies) * energy_scaling,
                    timestep=timestep  # DON'T MESS WITH THIS
                )
                sim.propagate(steps)
                coords = sim.extract_trajectory(flatten=True, embed=mol.coords)

            plops = pot(coords)[0] * 219475
            coords = coords[plops < 7500]

            run_dir = os.path.join(base_dir, f'steps_{steps}/')

            os.makedirs(run_dir, exist_ok=True)
            if save_plots:
                plot_dir=run_dir
                os.makedirs(plot_dir, exist_ok=True)
            print("="*25, "Steps:", steps, "="*25)

            cartesians = False
            if use_interpolation:
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

            with BlockProfiler(inactive=True):
                print("="*25, steps, "="*25)
                crd = np.round(coords[len(initial_energies)-1:], 8)

                dgb = DGB.construct(
                    crd,
                    pot if not use_interpolation else interp_data,
                    masses=mol.atomic_masses,
                    alphas={'method':'virial', 'scaling':virial_scaling},
                    transformations=transformation_method,
                    optimize_centers=[
                        x for x in
                        [
                            {
                                'method': 'energy-cutoff',
                                'cutoff': pruning_energy
                            } if pruning_energy is not None else None,
                            10 ** (-gs_cutoff) if gs_cutoff is not None else None,
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
                            np.linalg.norm(mol.coords[0] - mol.coords[1])
                            ),
                        (0, 2): self.setupMorseFunction(
                            mol.atomic_masses[0],
                            mol.atomic_masses[2],
                            np.linalg.norm(mol.coords[0] - mol.coords[2])
                        )
                    } if use_pairwise_potentials else None,
                    logger=True
                )

                if use_interpolation and plot_interpolation_error:
                    sel = slice(None)#slice(15,30)

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

                    inter_trace = np.trace(np.abs(inter_hess[ords]), axis1=1, axis2=2)
                    real_trace = np.trace(np.abs(real_hess[ords]), axis1=1, axis2=2)
                    dev_trace = inter_trace - real_trace
                    dev_hess = inter_hess - real_hess
                    hess_plot = plt.ScatterPlot(realpots[ords], inter_trace - real_trace)

                    print("Mean Absolute Error:", np.mean(np.abs(unscaled_devs)), "Std:", np.std(unscaled_devs))
                    print("Mean Scaled Error:", np.mean(np.abs(devs)), "Std:", np.std(devs))
                    print("Mean Hessian Error:", np.mean(np.abs(dev_hess.flatten())), "Std:", np.std(dev_hess.flatten()))
                    print("Maximum (Scaled) Error:", devs[max_dev_pos])
                    print("Maximum (Scaled) Interpolation Error:")
                    for l,r,c, tt, ii, ov in zip(
                            dgb.gaussians.coords.centers[rows[sel][ords[max_dev_pos]]],
                            dgb.gaussians.coords.centers[cols[sel][ords[max_dev_pos]]],
                            dgb.gaussians.overlap_data['centers'][sel][ords[max_dev_pos]],
                            realpots[ords[max_dev_pos]],
                            interpots[ords[max_dev_pos]],
                            utris[sel][ords[max_dev_pos]]
                    ):
                        print(f"Centers: {c} ({ov}) <- {l} {r}")
                        print(f"  Error: {ii-tt} <- {tt} {ii}")

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
                    dev_plot = plt.ScatterPlot(realpots[ords], devs,
                                               plot_range=[None, [-100, 100]]
                                               )
                    dev_plot.show()
                    raise Exception(...)

                if plot_gaussians:
                    self.plot_gaussians(dgb, mol, 10)

                    raise Exception(...)

                type(self).default_num_plot_wfns = 5
                if plot_wfns is True:
                    plot_wfns = {'cartesians': [0, 1]} if not cartesians else True
                self.runDGB(dgb, mol,
                            mode='similarity',
                            similarity_cutoff=sim_cutoff,
                            plot_wavefunctions=plot_wfns,
                            plot_spectrum=False,
                            pot_cmap=pot_cmap,
                            pot_points=pot_points,
                            wfn_cmap=wfn_cmap,
                            wfn_points=wnf_points,
                            wfn_contours=wfn_contours,
                            plot_centers=plot_centers
                            )

    @debugTest
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

        with GaussianLogReader(TestManager.test_data('h2o_aimd.log')) as parser:
            traj_data = parser.parse(['AIMDCoordinates', 'AIMDValues'])

        coords = traj_data['AIMDCoordinates']#[:3]
        grads = traj_data['AIMDValues'].gradients#[:3]
        hess = traj_data['AIMDValues'].hessians#[:3]
        engs = traj_data['AIMDValues'].energies - ref_eng
        dupe_pos = np.where(np.diff(engs) == 0)[0]
        good_pos = np.setdiff1d(np.arange(len(engs)), dupe_pos)
        engs = engs[good_pos]
        grads = grads[good_pos]
        hess = hess[good_pos]

        nvals = min(len(coords), len(grads))
        coords = coords[:nvals]
        grads = grads[:nvals]
        hess = hess[:nvals]

        # pot, _ = self.getMBPolModel()
        # new_engs, new_grads, new_hess = pot(coords, deriv_order=2)
        # engs = new_engs

        # raise Exception(
        #     np.abs(engs - new_engs) * 219475
        # )

        embedding = mol.get_embedding_data(coords)

        rots = embedding.rotations[:nvals]
        # need to reembed derivatives
        coords = nput.vec_tensordot(
            embedding.coord_data.coords,
            rots,
            shared=1,
            axes=[-1, 2]
        )
        grads = nput.vec_tensordot(
            grads.reshape(grads.shape[0], -1, 3),
            embedding.rotations,
            shared=1,
            axes=[-1, 2]
        ).reshape(grads.shape)
        hess = nput.vec_tensordot(
            hess.reshape((hess.shape[0], hess.shape[1], -1, 3)),
            embedding.rotations,
            shared=1,
            axes=[-1, 2]
        ).reshape(hess.shape)
        hess = np.moveaxis(hess, -2, -1)
        hess = nput.vec_tensordot(
            hess.reshape((hess.shape[0], hess.shape[1], -1, 3)),
            embedding.rotations,
            shared=1,
            axes=[-1, 2]
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
                    similarity_cutoff=.5,
                    plot_centers=True,
                    # mode='classic',
                    # subspace_size=15,
                    plot_spectrum=False,
                    plot_wavefunctions={'cartesians': [0, 1]}
                    )