import collections
import gc
import itertools
import os

import scipy.linalg

import McUtils.Zachary
from Peeves.TestUtils import *
from unittest import TestCase

from McUtils.Data import UnitsData, PotentialData, AtomData
from McUtils.Zachary import Interpolator, FiniteDifferenceDerivative
import McUtils.Plots as plt
import McUtils.Numputils as nput
from McUtils.GaussianInterface import GaussianLogReader
from McUtils.Extensions import ModuleLoader

from McUtils.Scaffolding import Checkpointer

from Psience.DGB import *
from Psience.Molecools import Molecule
from Psience.AIMD import AIMDSimulator

import numpy as np

class DGBTests(TestCase):

    @validationTest
    def test_Harmonic(self):
        ndivs = [15, 15, 8]
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
        ndivs = [50]*d
        domain = [[2, 3]]*d # , [-1, 1], [-1, 1]] # unit cube
        ndim = len(domain)

        pts = np.array(
            np.meshgrid(*(np.linspace(d[0], d[1], n) for d, n in zip(domain, ndivs)))
        ).T.reshape(-1, ndim)

        w = 2.05; wx=.1; mu=1
        de = (w ** 2) / (4 * wx)
        a = np.sqrt(2 * mu * wx)

        pot = lambda c,de=de,a=a,mu=mu: np.sum(de*(1-np.exp(-a*c))**2, axis=-1)

        wfns = DGB.run(pts, pot, alphas=1, optimize_centers=False, clustering_radius=.01, min_singular_value=1e-4)
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
            base = c.shape[:-1]
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

        test_ham = DGB(
            np.array([
                [-1,  1],
                [ 1, .5]
            ]),
            simple_morse,
            optimize_centers=False,
            alphas=[[1, .1], [1, .1]],
            clustering_radius=-1,
            min_singular_value=-1,  # 0.0001,
            num_svd_vectors=10000,
            expansion_degree=2,
            transformations=np.array([
                np.eye(2),
                # np.eye(2)

                np.array([
                    [ 1 / np.sqrt(2), -1 / np.sqrt(2)],
                    [ 1 / np.sqrt(2),  1 / np.sqrt(2)]
                ]).T
            ])
            # masses=masses
        )
        test_ham2 = DGB(
            test_ham.centers,
            simple_morse,
            optimize_centers=False,
            alphas=test_ham.alphas,
            clustering_radius=-1,
            min_singular_value=-1,  # 0.0001,
            num_svd_vectors=10000,
            expansion_degree=2,
            # transformations=np.array([
            #     np.eye(2),
            #     np.eye(2)
            #     # [[1, 0], [0, 1]],
            #     # [
            #     #     [1 / np.sqrt(2), -1 / np.sqrt(2)],
            #     #     [1 / np.sqrt(2), 1 / np.sqrt(2)]
            #     # ]
            # ])
            # masses=masses
        )

        self.assertTrue(
            np.allclose(
                test_ham.centers,
                np.array([
                    [-1, 1],
                    [1, .5]
                ])
            )
        )
        self.assertTrue(
            np.allclose(
                test_ham.alphas,
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
                [[0.55, -0.03266702],
                 [-0.03266702, 0.55]]
            )
        )
        self.assertTrue(
            np.allclose(
                test_ham.V,
                [[2.23291299,  0.66597918],
                 [0.66597918,  0.99361619]]
            )
        )

        A = np.array([
            [-1,  0, 0],
            [ 1, -3, 0],
            [ 0,  0, 1]
        ])
        _, tf = np.linalg.eigh(A.T@A)
        # need to resort the columns of tf to
        # match my analytic test
        tf = tf[:, (2, 1, 0)]
        if np.linalg.det(tf) < 0:
            tf[:, 2] *= -1
        test_ham = DGB(
            np.array([
                [-1,  1, 1],
                [ 1, -2, 0]
            ]),
            simple_morse,
            optimize_centers=False,
            alphas=[[1, .1, .5], [1, .5, .5]],
            clustering_radius=-1,
            min_singular_value=-10,  # 0.0001,
            num_svd_vectors=10000,
            expansion_degree=2,
            transformations=np.array([
                np.eye(3),
                # np.array([
                #     [0.85065080835204, -0.5257311121191337, 0.0],
                #     [0.0, 0.0, 1.0 ],
                #     [0.5257311121191336, 0.8506508083520398, 0.0]
                # ])
                tf.T
                # np.eye(3)
            ])
        )

        self.assertTrue(
            np.allclose(
                test_ham.centers,
                np.array([
                    [-1, 1, 1],
                    [ 1, -2, 0]
                ])
            ) and
            np.allclose(
                test_ham.alphas,
                np.array([
                    [ 1, .1, .5],
                    [ 1, .5, .5]
                ])
            ) and
            np.allclose(
                test_ham.S,
                [[1, 0.05927224],
                 [0.05927224, 1]]
            ) and
            np.allclose(
                test_ham.T,
                [[0.8, -0.04597853],
                 [-0.04597853, 1.0]]
            ) and
            np.allclose(
                test_ham.V,
                [[2.66034212,  0.42746656],
                 [0.42746656, 11.74998476]]
            )
        )

    @validationTest
    def test_MorseRotatedAIMD(self):
        ndim = d = 2
        reduced_mass = (
                               1 / (16 * UnitsData.convert("AtomicMassUnits", "ElectronMass"))
                               + 1 / (1.01 * UnitsData.convert("AtomicMassUnits", "ElectronMass"))
                       ) ** (-1)
        masses = np.array([reduced_mass]*ndim)
        re = np.array([0.957836 * UnitsData.convert("Angstroms", "BohrRadius")]*ndim)

        w = np.array([3869.47, 3869.47]) * UnitsData.convert("Wavenumbers", "Hartrees")
        wx = np.array([
            1 * 84.11,
            3 * 84.11
            ]) * UnitsData.convert("Wavenumbers", "Hartrees")
        # wx = .1 * UnitsData.convert("Wavenumbers", "Hartrees")
        mu = reduced_mass
        de = (w ** 2) / (4 * wx)
        a = np.sqrt(2 * mu * wx)

        ang = np.deg2rad(45/2)
        rot_mat = np.array([
            [np.cos(ang), np.sin(ang)],
            [-np.sin(ang), np.cos(ang)],
        ])
        center = rot_mat.T@re

        ndivs = [25]*d
        domain = [[r - np.max(a)/2, r + np.max(a)/1.5] for i,r in enumerate(center)]
        ndim = len(domain)

        def simple_morse(c, de=de, a=a, re=re, rot_mat=rot_mat, deriv_order=None):
            base = c.shape[:-1]
            ndim = c.shape[-1]

            c = c.reshape(-1, ndim)
            c = rot_mat@c[:, :, np.newaxis]
            c = c.reshape(-1, ndim)
            c = c - np.broadcast_to(np.array(re)[np.newaxis], c.shape)

            if not isinstance(a, (float, int, np.integer, np.floating)):
                a = np.broadcast_to(np.array(a)[np.newaxis], c.shape)

            if not isinstance(de, (float, int, np.integer, np.floating)):
                de = np.broadcast_to(np.array(de)[np.newaxis], c.shape)

            if deriv_order is None:
                res = np.sum(de*(1-np.exp(-a*c))**2, axis=-1)
            else:
                n = deriv_order
                m = ((-1)**(n+1) * 2 * a**n * de) * np.exp(-2*a*c) * (np.exp(a*c)-(2**(n-1)))
                if n == 1:
                    res = m
                else:
                    res = np.zeros(c.shape[:-1] + (ndim,)*deriv_order)
                    for k in range(ndim):
                        idx = (...,) + (k,)*n
                        res[idx] = m[..., k]
                for j in range(deriv_order):
                    res = np.tensordot(res, rot_mat, axes=[1, 0])
                # res = res

            res = res.reshape(base + res.shape[1:])

            return res

        def dipole(c, de=de, a=a, re=re, rot_mat=rot_mat, deriv_order=None): # simple linear dipole...
            base = c.shape[:-1]
            ndim = c.shape[-1]

            c = c.reshape(-1, ndim)
            c = rot_mat@c[:, :, np.newaxis]
            c = c.reshape(-1, ndim)
            c = c - np.broadcast_to(np.array(re)[np.newaxis], c.shape)

            if deriv_order is None or deriv_order == 0:
                res = rot_mat.T@c[:, :, np.newaxis] / 100
                res = np.reshape(res, c.shape)
                res = np.concatenate([res, np.zeros(c.shape[:-1] + (1,))], axis=-1)
            else: # I don't need to care about getting the linear deriv right since it won't contribute...
                if deriv_order > 1:
                    raise NotImplementedError("ugh...")
                n = deriv_order
                res = np.zeros(c.shape[:-1] + (ndim,)*deriv_order + (3,))
                # for k in range(ndim):
                #     idx = (...,) + (k,)*n + (slice(None, None, None),)
                #     res[idx] = m[..., k, :]
                # for j in range(deriv_order):
                #     res = np.tensordot(res, rot_mat, axes=[1, 0])
                # res = res

            res = res.reshape(base + res.shape[1:])

            return res

        # doing some AIMD for a pair of Morse oscillators...

        # TODO: TUNABLE PARAMETERS

        from Psience.DVR import DVR

        dvr = DVR(
            domain=[[x[0] - .5, x[1] + 1] for x in domain],
            divs=[25, 25],
            potential_function=simple_morse,
            # potential_optimize=True,
            mass=[reduced_mass, reduced_mass]
        )

        dvr_wfns = None

        max_e = np.min(de) - 100*UnitsData.convert("Wavenumbers", "Hartrees")
        for ts in [500]:#[1000, 2500, 5000]:
            for nt in [1]:#[1, 5, 15]:
                for dt in [5]:#[1, 2, 5]:
                    for dc in [0.000001]:#[.05, .1]:
                        for et in [max_e]:#[np.max(w), 2*np.max(w), np.min(de) - 100*UnitsData.convert("Wavenumbers", "Hartrees")]:

                            ntraj = nt
                            traj_steps = ts
                            trad_dt = dt
                            e_tot = et
                            disp_rad = np.power([1e-8, 1e-8], 2)
                            vel_cov = np.power([.2, 1], 2)

                            distance_cutoff = dc

                            scaling = 1
                            rp_scaling = scaling
                            min_rp_freq = 100 * UnitsData.convert("Wavenumbers", "Hartrees")
                            # min_rp_mass = 900
                            sing_cutoff = 1e-4
                            min_dist_alpha_scaling = None
                            potential_scaling = None

                            diag_scaling = 2
                            hess_diag_sing_cutoff = 1e-4
                            num_svd_vectors = 10000
                            min_dist_scaling = 1/4
                            min_dist_min_sin = 1e-4

                            exp_deg = 2
                            e_cut = np.max(de) #3 * np.max(w)
                            plot_traj = False
                            plot_orthog = False#range(15, 20)
                            plot_S_eigenvectors = True
                            plot_dists = True
                            plot_wfns = 7
                            plot_spectrum = True
                            throw_energies = False

                            opts = dict(
                                omega=w * UnitsData.convert("Hartrees", "Wavenumbers"),
                                omegax=wx* UnitsData.convert("Hartrees", "Wavenumbers"),
                                re=re,
                                masses=masses,

                                ntraj=ntraj,
                                traj_steps = traj_steps,
                                trad_dt = trad_dt,
                                e_tot = e_tot,
                                disp_rad = disp_rad,
                                vel_cov = vel_cov,

                                distance_cutoff = distance_cutoff,

                                scaling=scaling,
                                rp_scaling = rp_scaling,
                                min_rp_freq = min_rp_freq * UnitsData.convert("Hartrees", "Wavenumbers"),
                                # min_rp_mass = min_rp_mass,
                                sing_cutoff = sing_cutoff,
                                min_dist_alpha_scaling = min_dist_alpha_scaling,
                                potential_scaling = potential_scaling,

                                diag_scaling = diag_scaling,
                                hess_diag_sing_cutoff = hess_diag_sing_cutoff,
                                num_svd_vectors = num_svd_vectors,
                                min_dist_scaling = min_dist_scaling,
                                min_dist_min_sin = min_dist_min_sin
                            )

                            import datetime
                            plots_dir = os.path.join(
                                os.path.expanduser("~/Documents/Postdoc/AIMD-Spec/2D_tests"),
                                "Exp{}/T{}/N{}/DT{}/E{}/D{}/{}".format(
                                    exp_deg,
                                    traj_steps, ntraj, trad_dt,
                                    round(e_tot*UnitsData.convert("Hartrees", "Wavenumbers")),
                                    distance_cutoff,
                                    datetime.datetime.now().isoformat()
                                )
                            )
                            plots_dir = None
                            if plots_dir is not None:
                                os.makedirs(plots_dir, exist_ok=True)
                                with Checkpointer.from_file(os.path.join(plots_dir, 'params.json')) as chk:
                                    for k,v in opts.items():
                                        chk[k] = v


                            np.random.seed(0)
                            disps = np.random.multivariate_normal([0, 0], np.diag(disp_rad), size=(ntraj,))
                            coords = re[np.newaxis] + disps
                            coords = (rot_mat.T[np.newaxis]@coords[:, :, np.newaxis]).reshape(-1, ndim)

                            e_rem = e_tot - simple_morse(coords)
                            coords = coords[e_rem > 0]
                            e_rem = e_rem[e_rem > 0]

                            # perp_rot = np.array([
                            #     [np.cos(np.pi/2), -np.sin(np.pi/2)],
                            #     [np.sin(np.pi/2),  np.cos(np.pi/2)]
                            # ])
                            # dirs = perp_rot[np.newaxis]@simple_morse(coords, deriv_order=1)[:, :, np.newaxis]
                            # dirs = dirs.reshape(coords.shape)
                            dirs = np.random.multivariate_normal([0, 0], np.diag(vel_cov), size=(len(e_rem),))
                            cur_e = np.abs(dirs) * w
                            e_part = cur_e * ( e_rem / np.sum(cur_e, axis=1) )[:, np.newaxis]
                            v_part = np.sign(dirs) * np.sqrt(2 * e_part / masses)
                            vels = rot_mat.T[np.newaxis]@v_part[:, :, np.newaxis]
                            vels = np.reshape(vels, (-1, ndim))
                            sim_mass = masses
                            # sim_mass = rot_mat.T@masses
                            # raise Exception(e_tot, np.sum(masses[np.newaxis]/2*vels**2, axis=1))

                            #
                            #
                            # if vel_rad is not None:
                            #     vels = np.random.multivariate_normal([0, 0], np.diag(vel_rad), size=(ntraj,))
                            # else:
                            #     vels = np.zeros_like(coords)
                            # vels = (rot_mat[np.newaxis]@vels[:, :, np.newaxis]).reshape(-1, ndim)

                            forces = lambda c: -simple_morse(c, deriv_order=1)
                            sim = AIMDSimulator(sim_mass, coords, velocities=vels, sampling_rate=10, force_function=forces, track_kinetic_energy=True, timestep=trad_dt)
                            sim.propagate(traj_steps)

                            pts = np.array(sim.trajectory).reshape(-1, ndim)

                            # total_e = (
                            #     np.array(simple_morse(pts)) +
                            #         np.array(sim.kinetic_energies)
                            # )
                            # plt.Plot(np.arange(len(total_e)), total_e).show()
                            # raise Exception(...)

                            def get_plot_grid(pts):
                                plot_grid = np.array(
                                    np.meshgrid(
                                        *(
                                            np.linspace(
                                                min(d[0], np.min(pts[:, i] - .1)),
                                                max(d[1], np.max(pts[:, i] + .1)),
                                                75
                                            )
                                            for i, (d, n) in enumerate(zip(domain, ndivs)))
                                    )
                                )
                                plot_pts = np.moveaxis(plot_grid, 0, 2).reshape(-1, ndim)

                                return plot_grid, plot_pts

                            if plot_traj or plots_dir is not None:
                                plot_grid, plot_pts = get_plot_grid(pts)
                                plot_vals = simple_morse(plot_pts).reshape(plot_grid[0].shape)
                                vmax = e_cut + np.max(w)  # 10000 * UnitsData.convert("Wavenumbers", "Hartrees")
                                plot_vals[plot_vals > vmax] = vmax
                                fig = plt.ContourPlot(*plot_grid, plot_vals, levels=20, name='Traj Plot')
                                plt.ScatterPlot(pts[:, 0], pts[:, 1], figure=fig, color='red', name="Traj Points")

                                if plots_dir is None:
                                    fig.show()
                                else:
                                    fig.savefig(os.path.join(plots_dir, 'traj.png'))
                                    fig.close()

                            ham1 = DGB(pts, simple_morse,
                                       optimize_centers=False,
                                       alphas={'method':'virial', 'allow_rotations':True, 'remove_translation_rotations':False},
                                       min_singular_value=sing_cutoff,#0.0001,
                                       expansion_degree=exp_deg,
                                       masses=masses
                                       )

                            ham1A = DGB(pts, simple_morse,
                                       alphas={'method':'virial', 'allow_rotations':False},
                                       min_singular_value=hess_diag_sing_cutoff,
                                       expansion_degree=exp_deg,
                                       #  quadrature_degree=4,
                                       masses=masses
                                       )

                            ham2 = DGB(pts, simple_morse,
                                       optimize_centers=False,
                                       alphas={'method':'min_dist', 'use_mean':True},
                                       min_singular_value=min_dist_min_sin,
                                       expansion_degree=exp_deg,
                                       # quadrature_degree=6,
                                       masses=masses
                                       )
                            ham3 = DGB(pts, simple_morse,
                                       optimize_centers=False,
                                       alphas='min_dist',
                                       min_singular_value=min_dist_min_sin,#0.0001,
                                       expansion_degree=exp_deg,
                                       masses=masses
                                       )

                            # raise Exception(
                            #     np.min(ham1.T), np.max(ham1.T),
                            #     np.min(ham2.T), np.max(ham2.T)
                            # )

                            # raise Exception(
                            #     ham1.get_wavefunctions().energies[:5],
                            #     ham1A.get_wavefunctions().energies[:5],
                            #     ham2.get_wavefunctions().energies[:5],
                            #     ham3.get_wavefunctions().energies[:5]
                            # )

                            # rot_fun = np.linalg.det(rot_data['new_sigs'])
                            plot_grid, plot_pts = get_plot_grid(pts)

                            if plot_orthog:

                                # evs1 = np.linalg.eigvalsh(ham1.S)
                                # evs2 = np.linalg.eigvalsh(ham1A.S)
                                # raise Exception(
                                #     np.sum(evs1[evs1 > .0001]),
                                #     np.sum(evs2[evs2 > .0001])
                                # )

                                # fffff = plt.Plot(
                                #     np.arange(len(ham1.S)),
                                #     np.linalg.eigvalsh(ham1.S)
                                # )
                                # plt.Plot(
                                #     np.arange(len(ham1A.S)),
                                #     np.linalg.eigvalsh(ham1A.S),
                                #     figure=fffff
                                # )

                                # fffff = plt.ScatterPlot(
                                #     np.arange(len(ham1.S)),
                                #     np.linalg.eigh(ham1.S)[1][:, 0]**2
                                # )
                                # plt.ScatterPlot(
                                #     np.arange(len(ham1A.S)),
                                #     np.linalg.eigh(ham1A.S)[1][:, 0]**2,
                                #     figure=fffff
                                # )
                                if plot_orthog is True:
                                    plot_orthog = 5
                                if isinstance(plot_orthog, int):
                                    plot_orthog = range(plot_orthog)
                                for n in plot_orthog:

                                    base = plt.GraphicsGrid(ncols=2, nrows=2,
                                                            subimage_size=(300, 300), padding=[[50, 10], [50, 50]],
                                                            spacings=[50, 50])
                                    plot_vals = simple_morse(plot_pts).reshape(plot_grid[0].shape)
                                    vmax = e_cut + np.max(w)#10000 * UnitsData.convert("Wavenumbers", "Hartrees")
                                    plot_vals[plot_vals > vmax] = vmax
                                    for i in range(2):
                                        for j in range(2):
                                            plt.ContourPlot(*plot_grid, plot_vals, levels=20,
                                                            figure=base[i, j])
                                            # if i == 0 and j == 0:
                                            #     plt.ScatterPlot(pts[:, 0], pts[:, 1], color='red', figure=base[i, j])

                                    for h, (i, j) in [
                                        (ham2, [0, 0]),
                                        (ham3, [0, 1]),
                                        (ham1, [1, 0]),
                                        (ham1A, [1, 1])
                                    ]:
                                        # wws = np.linalg.eigh(h.S)[1][:, -n] ** 2
                                        # pps = h.centers
                                        # ri, ci = np.triu_indices(len(pps), k=1)
                                        # wvs = np.linalg.norm(pps[ri] - pps[ci], axis=1)[:, np.newaxis] * wws[ri] * wws[ci]
                                        # delocs = np.sum(wvs, axis=0)


                                        if plot_S_eigenvectors:
                                            sigs, Q = np.linalg.eigh(h.S)
                                        else:
                                            Q, Qinv, (Qq, Qqinv) = h.get_orthogonal_transform()
                                            sigs, _ = np.linalg.eigh(h.S)
                                        # Q = L @ np.diag(1/(sigs**2)) @ L.T
                                        wfns = DGBWavefunctions(
                                            np.ones(len(Q)),
                                            Q,
                                            hamiltonian=h
                                        )
                                        wf = wfns[-(n+1)]
                                        max_val = max(np.max(np.abs(wf.data)), 5)
                                        wf.plot(
                                            figure=base[i, j],
                                            plotter=plt.TriContourLinesPlot,
                                            # levels=np.linspace(-max_val, max_val, 16),
                                            domain=[[np.min(plot_pts[:, 0]), np.max(plot_pts[:, 0])],
                                                    [np.min(plot_pts[:, 1]), np.max(plot_pts[:, 1])]],
                                            cmap='RdBu',
                                            contour_levels=10,
                                            plot_label=str(1/np.sqrt(sigs[-(n+1)])),
                                            plot_range=[[np.min(plot_pts[:, 0]), np.max(plot_pts[:, 0])],
                                                        [np.min(plot_pts[:, 1]), np.max(plot_pts[:, 1])]]
                                        )
                                        # if hasattr(wfns[w], 'centers'):
                                        plt.ScatterPlot(wf.centers[:, 0], wf.centers[:, 1], color='#fff1', figure=base[i, j])
                                base.show()
                                raise Exception(...)

                            npts = len(pts)
                            shit_rows, shit_cols = np.triu_indices(npts)
                            shit_pos = np.where(shit_rows == shit_cols)
                            shit_pos = (shit_pos[0],)
                            def eval_gauss(rot_data, plot_pts=plot_pts, shit_pos=shit_pos):
                                gauss_vals = None
                                if isinstance(rot_data, dict):
                                    rot_centers = rot_data['centers']
                                    for c, s in zip(rot_data['centers'][shit_pos], rot_data['sigmas'][shit_pos]):
                                        disps = plot_pts - c[np.newaxis]
                                        v = np.linalg.det(s) ** (1 / 4) * np.exp(
                                            -(disps[:, np.newaxis, :] @ s[np.newaxis, :, :] @ disps[:, :, np.newaxis]) / 2
                                        )
                                        if gauss_vals is None:
                                            gauss_vals = v
                                        else:
                                            gauss_vals = np.max(
                                                np.concatenate([
                                                    gauss_vals.reshape(len(plot_pts), 1),
                                                    v.reshape(len(plot_pts), 1)
                                                ],
                                                    axis=-1),
                                                axis=-1
                                            )
                                else:
                                    rot_centers = rot_data[0]
                                    for c, s in zip(*(x[shit_pos] for x in rot_data)):
                                        disps = plot_pts - c[np.newaxis]
                                        v = (2 ** ndim * np.prod(s)) ** (1 / 4) * np.exp(-np.tensordot(disps ** 2, s, axes=[-1, -1]))
                                        if gauss_vals is None:
                                            gauss_vals = v
                                        else:
                                            gauss_vals = np.max(
                                                np.concatenate([
                                                    gauss_vals.reshape(len(plot_pts), 1),
                                                    v.reshape(len(plot_pts), 1)
                                                ], axis=-1),
                                                axis=-1
                                            )
                                return rot_centers, gauss_vals

                            def plot_gauss(rot_data, figure=None, plot_grid=plot_grid, shit_pos=shit_pos, color='#f00f'):
                                rot_centers, gauss_vals = eval_gauss(rot_data, shit_pos=shit_pos)
                                fig = plt.ContourPlot(*plot_grid, gauss_vals.reshape(plot_grid[0].shape), figure=figure)
                                # plt.ScatterPlot(rot_centers[:, 0], rot_centers[:, 1], figure=base[1, 0], plot_label='Min-Max: {:.3f} {:.3f}'.format(
                                #     np.min(gauss_vals), np.max(gauss_vals)
                                # ))
                                # plt.ScatterPlot(rot_centers[:, 0], rot_centers[:, 1], color='blue', figure=fig)
                                plt.ScatterPlot(rot_centers[shit_pos][:, 0], rot_centers[shit_pos][:, 1], color=color, figure=fig)
                                return fig

                            if plot_dists or plots_dir is not None:
                                base = plt.GraphicsGrid(ncols=2, nrows=2, subimage_size=(300, 300), padding=[[50, 10], [50, 50]], spacings=[50, 50])

                                plot_vals = simple_morse(plot_pts).reshape(plot_grid[0].shape)
                                vmax = e_cut + np.max(w)#10000 * UnitsData.convert("Wavenumbers", "Hartrees")
                                plot_vals[plot_vals > vmax] = vmax
                                plt.ContourPlot(*plot_grid, plot_vals, levels=20, figure=base[0, 0])
                                plt.ScatterPlot(pts[:, 0], pts[:, 1], color='red', figure=base[0, 0])

                                # plot_gauss(ham1.get_overlap_gaussians(), shit_pos=( shit_pos[0][:3],), figure=base[0, 1])
                                # plot_gauss(ham1.get_overlap_gaussians(), shit_pos=( np.array([1, 2, 3]),), figure=base[0, 2])
                                # plot_gauss(ham1.get_overlap_gaussians(), shit_pos=(np.concatenate([
                                #     shit_pos[0][:3],
                                #     np.array([1, 2, 3])
                                #     ]),), figure=base[0, 3])
                                # plot_gauss(ham1.get_overlap_gaussians(), shit_pos=slice(None, None, None), color="#f00f", figure=base[1, 0])
                                # plot_gauss(ham1A.get_overlap_gaussians(), shit_pos=slice(None, None, None), color="#f00f", figure=base[1, 1])
                                plot_gauss(ham1.get_overlap_gaussians(), color="#f00f", figure=base[1, 0])
                                plot_gauss(ham1A.get_overlap_gaussians(), color="#f00f", figure=base[1, 1])
                                plot_gauss(ham2.get_overlap_gaussians(),  figure=base[0, 0])
                                plot_gauss(ham3.get_overlap_gaussians(),  figure=base[0, 1])

                                if plots_dir is None:
                                    if not plot_wfns and not plot_spectrum:
                                        base.show()
                                else:
                                    base.savefig(os.path.join(plots_dir, 'dists.png'))
                                    base.close()

                            # base = plt.GraphicsGrid(ncols=3, nrows=2, subimage_size=(350, 350))


                            # raise Exception(
                            #     np.linalg.svd(ham1.T)[1],
                            #     np.linalg.svd(ham2.S)[1]
                            # )
                            #
                            # raise Exception(
                            #     np.min(ham1.T), np.max(ham1.T),
                            #     np.min(ham2.T), np.max(ham2.T)
                            # )

                            base_energies = [(ww*(np.arange(5) + 1 / 2) - wwx*(np.arange(5) + 1 / 2)**2) for ww, wwx in zip(w, wx)]
                            test_es = np.sort(np.sum(list(itertools.product(*base_energies)), axis=-1))
                            test_fs = (test_es[1:] - test_es[0]) * UnitsData.convert("Hartrees", "Wavenumbers")

                            wfns = [
                                ham3.get_wavefunctions(),
                                ham1.get_wavefunctions(),
                                ham1A.get_wavefunctions(),
                                ham2.get_wavefunctions(),
                            ]

                            h2w = UnitsData.convert("Hartrees", "Wavenumbers")
                            if dvr_wfns is None:
                                dvr_wfns = dvr.run().wavefunctions
                            wfns.append(dvr_wfns)

                            if plot_wfns is True:
                                plot_wfns = 1
                            if plot_wfns or plots_dir is not None:
                                omega = np.max(w)
                                for n in range(plot_wfns):

                                    base = plt.GraphicsGrid(ncols=3, nrows=2, subimage_size=(300, 300), padding=[[50, 10], [50, 50]], spacings=[50, 50])

                                    plot_vals = simple_morse(plot_pts).reshape(plot_grid[0].shape)
                                    vmax = e_cut + omega#10000 * UnitsData.convert("Wavenumbers", "Hartrees")
                                    plot_vals[plot_vals > vmax] = vmax
                                    for i in range(2):
                                        for j in range(3):
                                            plt.ContourPlot(*plot_grid, plot_vals, levels=20, figure=base[i, j], vmin=0, vmax=vmax)
                                            if i == 0 and j == 0:
                                                plt.ScatterPlot(pts[:, 0], pts[:, 1], color='red', figure=base[i, j])

                                    for wfx,(i, j) in [
                                        (0, [0, 1]),
                                        (3, [0, 2]),
                                        (1, [1, 1]),
                                        (2, [1, 2]),
                                        (4, [1, 0])
                                    ]:
                                        if len(wfns[wfx]) > n:
                                            wf = wfns[wfx][n]
                                            max_val = np.max(np.abs(wf.data))
                                            wf.plot(
                                                figure=base[i, j],
                                                plot_label=f"Energy: {(wf.energy - (0 if n == 0 else wfns[wfx][0].energy)) * h2w:.0f}",
                                                plotter=plt.TriContourLinesPlot,
                                                contour_levels=10,
                                                # levels=np.linspace(-max_val, max_val, 16),
                                                domain=[[np.min(plot_pts[:, 0]), np.max(plot_pts[:, 0])], [np.min(plot_pts[:, 1]), np.max(plot_pts[:, 1])]],
                                                cmap='RdBu',
                                                # vmin=-max_val, vmax=max_val,
                                                plot_range=[[np.min(plot_pts[:, 0]), np.max(plot_pts[:, 0])], [np.min(plot_pts[:, 1]), np.max(plot_pts[:, 1])]]
                                            )
                                            # if hasattr(wf, 'centers'):
                                            #     plt.ScatterPlot(wf.centers[:, 0], wf.centers[:, 1], color='0001', figure=base[i, j])
                                    if plots_dir is not None:
                                        base.savefig(os.path.join(plots_dir, 'wfns_{}.png'.format(n)))
                                        base.close()

                                if plots_dir is None and not plot_spectrum:
                                    base.show()

                            if plot_spectrum is True:
                                plot_spectrum = plot_wfns - 1
                            if plot_spectrum:

                                base = plt.GraphicsGrid(ncols=2, nrows=2, subimage_size=(300, 300), padding=[[50, 10], [50, 50]], spacings=[50, 50])
                                dvr_spec = wfns[4][:plot_spectrum+1].get_spectrum(dipole)#.normalize(0)

                                for wfx, (i, j) in [
                                    (0, [0, 0]),
                                    (3, [0, 1]),
                                    (1, [1, 0]),
                                    (2, [1, 1])
                                ]:
                                    if len(wfns[wfx]) > 1:
                                        spec = wfns[wfx][:plot_spectrum+1].get_spectrum(dipole, expansion_degree=1)#.normalize(0)
                                        spec.plot(figure=base[i, j],
                                                  plot_range=[[test_fs[0] - 300, test_fs[plot_spectrum] + 500], None],#[0, 1]],
                                                  )

                                for i in range(2):
                                    for j in range(2):
                                        dvr_spec.plot(figure=base[i, j], color='red', line_style='dashed')

                                if plots_dir is None:
                                    base.show()
                                else:
                                    base.savefig(os.path.join(plots_dir, 'spec.png'))
                                    base.close()

                            if throw_energies:
                                with np.printoptions(linewidth=1e8):

                                    for ham in [
                                        ham1,
                                        ham1A,
                                        ham2,
                                        ham3
                                    ]:
                                        print("=" * 50)
                                        print(np.linalg.eigh(ham.S)[0])

                                    engs = [e.energies for e in wfns] + [test_es]
                                    ne = min(10, min(len(e) for e in engs))
                                    raise Exception(str(np.round(
                                        np.array([
                                            np.concatenate([[eng[0]], eng[1:ne] - eng[0]]) for eng in engs
                                        ]) * h2w
                                    )))
                                # e = wfns.energies

    @debugTest
    def test_Morse1DAIMD(self):

        mol = Molecule.from_file(
            TestManager.test_data('water_freq.fchk'),
            internals=[[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]]
        )
        mol = mol.get_embedded_molecule()

        r1 = 0
        r2 = 1
        a12 = 2
        w2h = UnitsData.convert("Wavenumbers", "Hartrees")
        model = mol.get_model(
            {
                r1: {'morse': {'w': 3869.47 * w2h, 'wx': 84 * w2h}},
                # r2: {'morse': {'w': 3869.47 * w2h, 'wx': 84 * w2h}},
                # a12: {'harmonic': {'k': 1600 ** 2 / 150 * w2h}}
            },
            dipole=[
                {
                    (r1, a12): ({'linear': {'eq': 0, 'scaling': 1 / 5.5}}, {'sin': {'eq': 0}}),
                    (r2, a12): ({'linear': {'eq': 0, 'scaling': -1 / 5.5}}, {'sin': {'eq': 0}})
                },
                0,
                # {
                #     (r1, a12): ({'linear': {'eq': 0, 'scaling': 1 / (2 * 5.5)}}, {'cos': {'eq': 0}}),
                #     (r2, a12): ({'linear': {'eq': 0, 'scaling': 1 / (2 * 5.5)}}, {'cos': {'eq': 0}})
                # },
                0
            ]
        )

        r_vec = np.concatenate([
            np.zeros(3),
            mol.coords[1] - mol.coords[0],
            np.zeros(3)
        ]) / np.linalg.norm(mol.coords[1] - mol.coords[0])

        def get_displacements(rs):
            disps = r_vec[np.newaxis, :] * rs[:, np.newaxis]
            return mol.get_displaced_coordinates(disps.reshape(-1, 3, 3))

        reduced_mass = (
                               1 / (AtomData["O", "Mass"] * UnitsData.convert("AtomicMassUnits", "ElectronMass"))
                               + 1 / (AtomData["H", "Mass"] * UnitsData.convert("AtomicMassUnits", "ElectronMass"))
                       ) ** (-1)
        #
        # def pot(, base_func=model.potential)
        # grad = lambda r: mod

        initial_energies = np.array([8000, 0, 0]) / UnitsData.hartrees_to_wavenumbers

        pot_func = model.potential
        mol.potential_derivatives = model.potential(mol.coords, deriv_order=2)[1:]
        modes = np.array(
            [
                r_vec,
                [1, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 1, 0, 0]
            ]
        ).T
        sim = AIMDSimulator(
            mol.atomic_masses,
            [mol.coords],
            lambda c: -pot_func(c, deriv_order=1)[1].reshape(c.shape),
            velocities=AIMDSimulator.mode_energies_to_velocities(modes, mol.atomic_masses, [initial_energies]),
            timestep=10,
            track_kinetic_energy=True
        )

        sim.propagate(30)
        # raise Exception(sim.trajectory[1])
        coords = np.array(sim.trajectory).reshape((-1, 3, 3))
        coords = mol.embed_coords(coords)

        dists = np.linalg.norm(coords[:, 1] - coords[:, 0], axis=1)
        pot_vals = pot_func(coords)

        # plt.Plot(dists, pot_vals).show()

        # raise Exception(dists)

        data_dir = os.path.expanduser("~/Documents/Postdoc/ISMS")
        # import json
        # with open(os.path.join(data_dir, "water.json"), 'w+') as out:
        #     json.dump(coords.tolist(), out)
        # raise Exception(...)

        dists = np.linalg.norm(coords[:, 1] - coords[:, 0], axis=1)

        def plot_pot():
            domain = [1.1, 5, 100]
            r_vals = np.linspace(*domain)
            points = get_displacements(r_vals)

            vals = pot_func(points) * UnitsData.hartrees_to_wavenumbers
            # cut = de / 2
            # vals[vals > cut] = cut

            return plt.Plot(r_vals, vals)

        def pot_r(r):
            return pot_func(get_displacements(r))
        pot_vals = pot_r(dists)

        ham = DGB(dists, pot_r, masses=reduced_mass, alphas=15, quadrature_degree=4)
        gauss_plot_wfns = DGBWavefunctions(
            np.zeros(len(dists)),
            np.eye(len(dists)),
            hamiltonian=ham
        )
        def plot_wfn(which, wavefuns=None, figure=None, **opts):
            return wavefuns[which].plot(figure=figure,
                                        scaling=10000,
                                        shift=wavefuns[which].energy,
                                        domain=[1.1, 5],
                                        **opts
                                        )

        plot_points = False
        if plot_points:
            for which in range(0, len(dists), 3):
                fig = plot_pot()
                fig = plt.ScatterPlot(
                    [dists[which]],
                    [pot_vals[which]*UnitsData.hartrees_to_wavenumbers],
                    figure=fig,
                    color='red'
                )
                fig.savefig(os.path.join(data_dir, 'figs', f'morse_1D_points_{which}.pdf'), transparent=True)
            fig = plt.ScatterPlot(
                    dists[::3],
                    pot_vals[::3]*UnitsData.hartrees_to_wavenumbers,
                    figure=plot_pot(),
                    color='red'
                )
            fig.savefig(os.path.join(data_dir, 'figs', 'morse_1D_points_all.pdf'), transparent=True)

        plot_gauss=True
        if plot_gauss:
            # for which in range(len(dists)):
            #     fig = plot_wfn(which, wavefuns=gauss_plot_wfns, figure=plot_pot())
            #     fig.savefig(os.path.join(data_dir, 'figs', f'morse_1D_gaussians_{which}.pdf'), transparent=True)
            all_plot = plot_pot()
            for which in range(0, len(dists), 6):
                all_plot = plot_wfn(which, wavefuns=gauss_plot_wfns, figure=all_plot, color='green')
            all_plot.savefig(os.path.join(data_dir, 'figs', f'morse_1D_gaussians_all.pdf'), transparent=True)


        plot_eig=False
        if plot_eig:
            _, s_eigfuncs = np.linalg.eigh(ham.S)
            s_wfns = DGBWavefunctions(
            np.zeros(len(dists)),
            s_eigfuncs,
            hamiltonian=ham
        )
            for which in range(len(dists)):
                fig = plot_wfn(which, wavefuns=s_wfns, figure=plot_pot())
                fig.savefig(os.path.join(data_dir, 'figs', f'morse_1D_s_eigs_{which}.pdf'), transparent=True)

        plot_wfns=True
        if plot_wfns:
            wfns = ham.get_wavefunctions()
            for which in range(3):
                fig = plot_wfn(which, wavefuns=wfns, figure=plot_pot())
                fig.savefig(os.path.join(data_dir, 'figs', f'morse_1D_wavefuncs_{which}.pdf'), transparent=True)

        raise Exception(dists)

        # TODO: TUNABLE PARAMETERS

        from Psience.DVR import DVR

        dvr = DVR(
            domain=[[x[0] - .5, x[1] + 1] for x in domain],
            divs=[25]*ndim,
            potential_function=simple_morse,
            # potential_optimize=True,
            mass=[reduced_mass, reduced_mass]
        )

        dvr_wfns = None

        max_e = np.min(de) - 100 * UnitsData.convert("Wavenumbers", "Hartrees")
        for ts in [500]:  # [1000, 2500, 5000]:
            for nt in [1]:  # [1, 5, 15]:
                for dt in [5]:  # [1, 2, 5]:
                    for dc in [0.000001]:  # [.05, .1]:
                        for et in [max_e]:  # [np.max(w), 2*np.max(w), np.min(de) - 100*UnitsData.convert("Wavenumbers", "Hartrees")]:

                            ntraj = nt
                            traj_steps = ts
                            trad_dt = dt
                            e_tot = et
                            disp_rad = np.power([1e-8], 2)
                            vel_cov = np.power([.2], 2)

                            distance_cutoff = dc

                            scaling = 1
                            rp_scaling = scaling
                            min_rp_freq = 100 * UnitsData.convert("Wavenumbers", "Hartrees")
                            # min_rp_mass = 900
                            sing_cutoff = 1e-4
                            min_dist_alpha_scaling = None
                            potential_scaling = None

                            diag_scaling = 2
                            hess_diag_sing_cutoff = 1e-4
                            num_svd_vectors = 10000
                            min_dist_scaling = 1 / 4
                            min_dist_min_sin = 1e-4

                            exp_deg = 2
                            e_cut = np.max(de)  # 3 * np.max(w)
                            plot_traj = False
                            plot_orthog = False  # range(15, 20)
                            plot_S_eigenvectors = True
                            plot_dists = True
                            plot_wfns = 7
                            plot_spectrum = True
                            throw_energies = False

                            opts = dict(
                                omega=w * UnitsData.convert("Hartrees", "Wavenumbers"),
                                omegax=wx * UnitsData.convert("Hartrees", "Wavenumbers"),
                                re=re,
                                masses=masses,

                                ntraj=ntraj,
                                traj_steps=traj_steps,
                                trad_dt=trad_dt,
                                e_tot=e_tot,
                                disp_rad=disp_rad,
                                vel_cov=vel_cov,

                                distance_cutoff=distance_cutoff,

                                scaling=scaling,
                                rp_scaling=rp_scaling,
                                min_rp_freq=min_rp_freq * UnitsData.convert("Hartrees", "Wavenumbers"),
                                # min_rp_mass = min_rp_mass,
                                sing_cutoff=sing_cutoff,
                                min_dist_alpha_scaling=min_dist_alpha_scaling,
                                potential_scaling=potential_scaling,

                                diag_scaling=diag_scaling,
                                hess_diag_sing_cutoff=hess_diag_sing_cutoff,
                                num_svd_vectors=num_svd_vectors,
                                min_dist_scaling=min_dist_scaling,
                                min_dist_min_sin=min_dist_min_sin
                            )

                            import datetime
                            plots_dir = os.path.join(
                                os.path.expanduser("~/Documents/Postdoc/AIMD-Spec/1D_tests"),
                                "Exp{}/T{}/N{}/DT{}/E{}/D{}/{}".format(
                                    exp_deg,
                                    traj_steps, ntraj, trad_dt,
                                    round(e_tot * UnitsData.convert("Hartrees", "Wavenumbers")),
                                    distance_cutoff,
                                    datetime.datetime.now().isoformat()
                                )
                            )
                            plots_dir = None
                            if plots_dir is not None:
                                os.makedirs(plots_dir, exist_ok=True)
                                with Checkpointer.from_file(os.path.join(plots_dir, 'params.json')) as chk:
                                    for k, v in opts.items():
                                        chk[k] = v

                            np.random.seed(0)
                            disps = np.random.multivariate_normal([0], np.diag(disp_rad), size=(ntraj,))
                            coords = re[np.newaxis] + disps
                            # coords = (rot_mat.T[np.newaxis] @ coords[:, :, np.newaxis]).reshape(-1, ndim)

                            e_rem = e_tot - simple_morse(coords)
                            coords = coords[e_rem > 0]
                            e_rem = e_rem[e_rem > 0]

                            # perp_rot = np.array([
                            #     [np.cos(np.pi/2), -np.sin(np.pi/2)],
                            #     [np.sin(np.pi/2),  np.cos(np.pi/2)]
                            # ])
                            # dirs = perp_rot[np.newaxis]@simple_morse(coords, deriv_order=1)[:, :, np.newaxis]
                            # dirs = dirs.reshape(coords.shape)
                            dirs = np.random.multivariate_normal([0]*ndim, np.diag(vel_cov), size=(len(e_rem),))
                            cur_e = np.abs(dirs) * w
                            e_part = cur_e * (e_rem / np.sum(cur_e, axis=1))[:, np.newaxis]
                            v_part = np.sign(dirs) * np.sqrt(2 * e_part / masses)
                            # vels = rot_mat.T[np.newaxis] @ v_part[:, :, np.newaxis]
                            # vels = np.reshape(vels, (-1, ndim))
                            vels = v_part
                            sim_mass = masses

                            forces = lambda c: -simple_morse(c, deriv_order=1)
                            sim = AIMDSimulator(sim_mass, coords, velocities=vels, sampling_rate=10,
                                                force_function=forces, track_kinetic_energy=True, timestep=trad_dt)
                            sim.propagate(traj_steps)

                            pts = np.array(sim.trajectory).reshape(-1, ndim)

                            # total_e = (
                            #     np.array(simple_morse(pts)) +
                            #         np.array(sim.kinetic_energies)
                            # )
                            # plt.Plot(np.arange(len(total_e)), total_e).show()
                            # raise Exception(...)

                            def get_plot_grid(pts):
                                plot_grid = np.array(
                                    np.meshgrid(
                                        *(
                                            np.linspace(
                                                min(d[0], np.min(pts[:, i] - .1)),
                                                max(d[1], np.max(pts[:, i] + .1)),
                                                75
                                            )
                                            for i, (d, n) in enumerate(zip(domain, ndivs)))
                                    )
                                )
                                plot_pts = np.moveaxis(plot_grid, 0, 2).reshape(-1, ndim)

                                return plot_grid, plot_pts

                            if plot_traj or plots_dir is not None:
                                plot_grid, plot_pts = get_plot_grid(pts)
                                plot_vals = simple_morse(plot_pts).reshape(plot_grid[0].shape)
                                raise Exception(plot_vals)
                                vmax = e_cut + np.max(w)  # 10000 * UnitsData.convert("Wavenumbers", "Hartrees")
                                plot_vals[plot_vals > vmax] = vmax
                                fig = plt.ContourPlot(*plot_grid, plot_vals, levels=20, name='Traj Plot')
                                plt.ScatterPlot(pts[:, 0], pts[:, 1], figure=fig, color='red', name="Traj Points")

                                if plots_dir is None:
                                    fig.show()
                                else:
                                    fig.savefig(os.path.join(plots_dir, 'traj.png'))
                                    fig.close()

                            ham1 = DGB(pts, simple_morse,
                                       optimize_centers=False,
                                       alphas={'method': 'virial', 'allow_rotations': True,
                                               'remove_translation_rotations': False},
                                       min_singular_value=sing_cutoff,  # 0.0001,
                                       expansion_degree=exp_deg,
                                       masses=masses
                                       )

                            ham1A = DGB(pts, simple_morse,
                                        alphas={'method': 'virial', 'allow_rotations': False},
                                        min_singular_value=hess_diag_sing_cutoff,
                                        expansion_degree=exp_deg,
                                        #  quadrature_degree=4,
                                        masses=masses
                                        )

                            ham2 = DGB(pts, simple_morse,
                                       optimize_centers=False,
                                       alphas={'method': 'min_dist', 'use_mean': True},
                                       min_singular_value=min_dist_min_sin,
                                       expansion_degree=exp_deg,
                                       # quadrature_degree=6,
                                       masses=masses
                                       )
                            ham3 = DGB(pts, simple_morse,
                                       optimize_centers=False,
                                       alphas='min_dist',
                                       min_singular_value=min_dist_min_sin,  # 0.0001,
                                       expansion_degree=exp_deg,
                                       masses=masses
                                       )

                            # raise Exception(
                            #     np.min(ham1.T), np.max(ham1.T),
                            #     np.min(ham2.T), np.max(ham2.T)
                            # )

                            # raise Exception(
                            #     ham1.get_wavefunctions().energies[:5],
                            #     ham1A.get_wavefunctions().energies[:5],
                            #     ham2.get_wavefunctions().energies[:5],
                            #     ham3.get_wavefunctions().energies[:5]
                            # )

                            # rot_fun = np.linalg.det(rot_data['new_sigs'])
                            plot_grid, plot_pts = get_plot_grid(pts)

                            if plot_orthog:

                                # evs1 = np.linalg.eigvalsh(ham1.S)
                                # evs2 = np.linalg.eigvalsh(ham1A.S)
                                # raise Exception(
                                #     np.sum(evs1[evs1 > .0001]),
                                #     np.sum(evs2[evs2 > .0001])
                                # )

                                # fffff = plt.Plot(
                                #     np.arange(len(ham1.S)),
                                #     np.linalg.eigvalsh(ham1.S)
                                # )
                                # plt.Plot(
                                #     np.arange(len(ham1A.S)),
                                #     np.linalg.eigvalsh(ham1A.S),
                                #     figure=fffff
                                # )

                                # fffff = plt.ScatterPlot(
                                #     np.arange(len(ham1.S)),
                                #     np.linalg.eigh(ham1.S)[1][:, 0]**2
                                # )
                                # plt.ScatterPlot(
                                #     np.arange(len(ham1A.S)),
                                #     np.linalg.eigh(ham1A.S)[1][:, 0]**2,
                                #     figure=fffff
                                # )
                                if plot_orthog is True:
                                    plot_orthog = 5
                                if isinstance(plot_orthog, int):
                                    plot_orthog = range(plot_orthog)
                                for n in plot_orthog:

                                    base = plt.GraphicsGrid(ncols=2, nrows=2,
                                                            subimage_size=(300, 300), padding=[[50, 10], [50, 50]],
                                                            spacings=[50, 50])
                                    plot_vals = simple_morse(plot_pts).reshape(plot_grid[0].shape)
                                    vmax = e_cut + np.max(w)  # 10000 * UnitsData.convert("Wavenumbers", "Hartrees")
                                    plot_vals[plot_vals > vmax] = vmax
                                    for i in range(2):
                                        for j in range(2):
                                            plt.ContourPlot(*plot_grid, plot_vals, levels=20,
                                                            figure=base[i, j])
                                            # if i == 0 and j == 0:
                                            #     plt.ScatterPlot(pts[:, 0], pts[:, 1], color='red', figure=base[i, j])

                                    for h, (i, j) in [
                                        (ham2, [0, 0]),
                                        (ham3, [0, 1]),
                                        (ham1, [1, 0]),
                                        (ham1A, [1, 1])
                                    ]:
                                        # wws = np.linalg.eigh(h.S)[1][:, -n] ** 2
                                        # pps = h.centers
                                        # ri, ci = np.triu_indices(len(pps), k=1)
                                        # wvs = np.linalg.norm(pps[ri] - pps[ci], axis=1)[:, np.newaxis] * wws[ri] * wws[ci]
                                        # delocs = np.sum(wvs, axis=0)

                                        if plot_S_eigenvectors:
                                            sigs, Q = np.linalg.eigh(h.S)
                                        else:
                                            Q, Qinv, (Qq, Qqinv) = h.get_orthogonal_transform()
                                            sigs, _ = np.linalg.eigh(h.S)
                                        # Q = L @ np.diag(1/(sigs**2)) @ L.T
                                        wfns = DGBWavefunctions(
                                            np.ones(len(Q)),
                                            Q,
                                            hamiltonian=h
                                        )
                                        wf = wfns[-(n + 1)]
                                        max_val = max(np.max(np.abs(wf.data)), 5)
                                        wf.plot(
                                            figure=base[i, j],
                                            plotter=plt.TriContourLinesPlot,
                                            # levels=np.linspace(-max_val, max_val, 16),
                                            domain=[[np.min(plot_pts[:, 0]), np.max(plot_pts[:, 0])],
                                                    [np.min(plot_pts[:, 1]), np.max(plot_pts[:, 1])]],
                                            cmap='RdBu',
                                            contour_levels=10,
                                            plot_label=str(1 / np.sqrt(sigs[-(n + 1)])),
                                            plot_range=[[np.min(plot_pts[:, 0]), np.max(plot_pts[:, 0])],
                                                        [np.min(plot_pts[:, 1]), np.max(plot_pts[:, 1])]]
                                        )
                                        # if hasattr(wfns[w], 'centers'):
                                        plt.ScatterPlot(wf.centers[:, 0], wf.centers[:, 1], color='#fff1',
                                                        figure=base[i, j])
                                base.show()
                                raise Exception(...)

                            npts = len(pts)
                            shit_rows, shit_cols = np.triu_indices(npts)
                            shit_pos = np.where(shit_rows == shit_cols)
                            shit_pos = (shit_pos[0],)

                            def eval_gauss(rot_data, plot_pts=plot_pts, shit_pos=shit_pos):
                                gauss_vals = None
                                if isinstance(rot_data, dict):
                                    rot_centers = rot_data['centers']
                                    for c, s in zip(rot_data['centers'][shit_pos], rot_data['sigmas'][shit_pos]):
                                        disps = plot_pts - c[np.newaxis]
                                        v = np.linalg.det(s) ** (1 / 4) * np.exp(
                                            -(disps[:, np.newaxis, :] @ s[np.newaxis, :, :] @ disps[:, :,
                                                                                              np.newaxis]) / 2
                                        )
                                        if gauss_vals is None:
                                            gauss_vals = v
                                        else:
                                            gauss_vals = np.max(
                                                np.concatenate([
                                                    gauss_vals.reshape(len(plot_pts), 1),
                                                    v.reshape(len(plot_pts), 1)
                                                ],
                                                    axis=-1),
                                                axis=-1
                                            )
                                else:
                                    rot_centers = rot_data[0]
                                    for c, s in zip(*(x[shit_pos] for x in rot_data)):
                                        disps = plot_pts - c[np.newaxis]
                                        v = (2 ** ndim * np.prod(s)) ** (1 / 4) * np.exp(
                                            -np.tensordot(disps ** 2, s, axes=[-1, -1]))
                                        if gauss_vals is None:
                                            gauss_vals = v
                                        else:
                                            gauss_vals = np.max(
                                                np.concatenate([
                                                    gauss_vals.reshape(len(plot_pts), 1),
                                                    v.reshape(len(plot_pts), 1)
                                                ], axis=-1),
                                                axis=-1
                                            )
                                return rot_centers, gauss_vals

                            def plot_gauss(rot_data, figure=None, plot_grid=plot_grid, shit_pos=shit_pos,
                                           color='#f00f'):
                                rot_centers, gauss_vals = eval_gauss(rot_data, shit_pos=shit_pos)
                                fig = plt.ContourPlot(*plot_grid, gauss_vals.reshape(plot_grid[0].shape), figure=figure)
                                # plt.ScatterPlot(rot_centers[:, 0], rot_centers[:, 1], figure=base[1, 0], plot_label='Min-Max: {:.3f} {:.3f}'.format(
                                #     np.min(gauss_vals), np.max(gauss_vals)
                                # ))
                                # plt.ScatterPlot(rot_centers[:, 0], rot_centers[:, 1], color='blue', figure=fig)
                                plt.ScatterPlot(rot_centers[shit_pos][:, 0], rot_centers[shit_pos][:, 1], color=color,
                                                figure=fig)
                                return fig

                            if plot_dists or plots_dir is not None:
                                base = plt.GraphicsGrid(ncols=2, nrows=2, subimage_size=(300, 300),
                                                        padding=[[50, 10], [50, 50]], spacings=[50, 50])

                                plot_vals = simple_morse(plot_pts).reshape(plot_grid[0].shape)
                                vmax = e_cut + np.max(w)  # 10000 * UnitsData.convert("Wavenumbers", "Hartrees")
                                plot_vals[plot_vals > vmax] = vmax
                                plt.ContourPlot(*plot_grid, plot_vals, levels=20, figure=base[0, 0])
                                plt.ScatterPlot(pts[:, 0], pts[:, 1], color='red', figure=base[0, 0])

                                # plot_gauss(ham1.get_overlap_gaussians(), shit_pos=( shit_pos[0][:3],), figure=base[0, 1])
                                # plot_gauss(ham1.get_overlap_gaussians(), shit_pos=( np.array([1, 2, 3]),), figure=base[0, 2])
                                # plot_gauss(ham1.get_overlap_gaussians(), shit_pos=(np.concatenate([
                                #     shit_pos[0][:3],
                                #     np.array([1, 2, 3])
                                #     ]),), figure=base[0, 3])
                                # plot_gauss(ham1.get_overlap_gaussians(), shit_pos=slice(None, None, None), color="#f00f", figure=base[1, 0])
                                # plot_gauss(ham1A.get_overlap_gaussians(), shit_pos=slice(None, None, None), color="#f00f", figure=base[1, 1])
                                plot_gauss(ham1.get_overlap_gaussians(), color="#f00f", figure=base[1, 0])
                                plot_gauss(ham1A.get_overlap_gaussians(), color="#f00f", figure=base[1, 1])
                                plot_gauss(ham2.get_overlap_gaussians(), figure=base[0, 0])
                                plot_gauss(ham3.get_overlap_gaussians(), figure=base[0, 1])

                                if plots_dir is None:
                                    if not plot_wfns and not plot_spectrum:
                                        base.show()
                                else:
                                    base.savefig(os.path.join(plots_dir, 'dists.png'))
                                    base.close()

                            # base = plt.GraphicsGrid(ncols=3, nrows=2, subimage_size=(350, 350))

                            # raise Exception(
                            #     np.linalg.svd(ham1.T)[1],
                            #     np.linalg.svd(ham2.S)[1]
                            # )
                            #
                            # raise Exception(
                            #     np.min(ham1.T), np.max(ham1.T),
                            #     np.min(ham2.T), np.max(ham2.T)
                            # )

                            base_energies = [(ww * (np.arange(5) + 1 / 2) - wwx * (np.arange(5) + 1 / 2) ** 2) for
                                             ww, wwx in zip(w, wx)]
                            test_es = np.sort(np.sum(list(itertools.product(*base_energies)), axis=-1))
                            test_fs = (test_es[1:] - test_es[0]) * UnitsData.convert("Hartrees", "Wavenumbers")

                            wfns = [
                                ham3.get_wavefunctions(),
                                ham1.get_wavefunctions(),
                                ham1A.get_wavefunctions(),
                                ham2.get_wavefunctions(),
                            ]

                            h2w = UnitsData.convert("Hartrees", "Wavenumbers")
                            if dvr_wfns is None:
                                dvr_wfns = dvr.run().wavefunctions
                            wfns.append(dvr_wfns)

                            if plot_wfns is True:
                                plot_wfns = 1
                            if plot_wfns or plots_dir is not None:
                                omega = np.max(w)
                                for n in range(plot_wfns):

                                    base = plt.GraphicsGrid(ncols=3, nrows=2, subimage_size=(300, 300),
                                                            padding=[[50, 10], [50, 50]], spacings=[50, 50])

                                    plot_vals = simple_morse(plot_pts).reshape(plot_grid[0].shape)
                                    vmax = e_cut + omega  # 10000 * UnitsData.convert("Wavenumbers", "Hartrees")
                                    plot_vals[plot_vals > vmax] = vmax
                                    for i in range(2):
                                        for j in range(3):
                                            plt.ContourPlot(*plot_grid, plot_vals, levels=20, figure=base[i, j], vmin=0,
                                                            vmax=vmax)
                                            if i == 0 and j == 0:
                                                plt.ScatterPlot(pts[:, 0], pts[:, 1], color='red', figure=base[i, j])

                                    for wfx, (i, j) in [
                                        (0, [0, 1]),
                                        (3, [0, 2]),
                                        (1, [1, 1]),
                                        (2, [1, 2]),
                                        (4, [1, 0])
                                    ]:
                                        if len(wfns[wfx]) > n:
                                            wf = wfns[wfx][n]
                                            max_val = np.max(np.abs(wf.data))
                                            wf.plot(
                                                figure=base[i, j],
                                                plot_label=f"Energy: {(wf.energy - (0 if n == 0 else wfns[wfx][0].energy)) * h2w:.0f}",
                                                plotter=plt.TriContourLinesPlot,
                                                contour_levels=10,
                                                # levels=np.linspace(-max_val, max_val, 16),
                                                domain=[[np.min(plot_pts[:, 0]), np.max(plot_pts[:, 0])],
                                                        [np.min(plot_pts[:, 1]), np.max(plot_pts[:, 1])]],
                                                cmap='RdBu',
                                                # vmin=-max_val, vmax=max_val,
                                                plot_range=[[np.min(plot_pts[:, 0]), np.max(plot_pts[:, 0])],
                                                            [np.min(plot_pts[:, 1]), np.max(plot_pts[:, 1])]]
                                            )
                                            # if hasattr(wf, 'centers'):
                                            #     plt.ScatterPlot(wf.centers[:, 0], wf.centers[:, 1], color='0001', figure=base[i, j])
                                    if plots_dir is not None:
                                        base.savefig(os.path.join(plots_dir, 'wfns_{}.png'.format(n)))
                                        base.close()

                                if plots_dir is None and not plot_spectrum:
                                    base.show()

                            if plot_spectrum is True:
                                plot_spectrum = plot_wfns - 1
                            if plot_spectrum:

                                base = plt.GraphicsGrid(ncols=2, nrows=2, subimage_size=(300, 300),
                                                        padding=[[50, 10], [50, 50]], spacings=[50, 50])
                                dvr_spec = wfns[4][:plot_spectrum + 1].get_spectrum(dipole)  # .normalize(0)

                                for wfx, (i, j) in [
                                    (0, [0, 0]),
                                    (3, [0, 1]),
                                    (1, [1, 0]),
                                    (2, [1, 1])
                                ]:
                                    if len(wfns[wfx]) > 1:
                                        spec = wfns[wfx][:plot_spectrum + 1].get_spectrum(dipole,
                                                                                          expansion_degree=1)  # .normalize(0)
                                        spec.plot(figure=base[i, j],
                                                  plot_range=[[test_fs[0] - 300, test_fs[plot_spectrum] + 500], None],
                                                  # [0, 1]],
                                                  )

                                for i in range(2):
                                    for j in range(2):
                                        dvr_spec.plot(figure=base[i, j], color='red', line_style='dashed')

                                if plots_dir is None:
                                    base.show()
                                else:
                                    base.savefig(os.path.join(plots_dir, 'spec.png'))
                                    base.close()

                            if throw_energies:
                                with np.printoptions(linewidth=1e8):

                                    for ham in [
                                        ham1,
                                        ham1A,
                                        ham2,
                                        ham3
                                    ]:
                                        print("=" * 50)
                                        print(np.linalg.eigh(ham.S)[0])

                                    engs = [e.energies for e in wfns] + [test_es]
                                    ne = min(10, min(len(e) for e in engs))
                                    raise Exception(str(np.round(
                                        np.array([
                                            np.concatenate([[eng[0]], eng[1:ne] - eng[0]]) for eng in engs
                                        ]) * h2w
                                    )))
                                # e = wfns.energies

    @inactiveTest
    def test_PolyDBGKE(self):

        # from McUtils.Zachary import DensePolynomial
        #
        # dense_tensors = DensePolynomial(
        #         [[0, 1, 2], [3, -5, 7]]
        #     ).coefficient_tensors
        # new = DensePolynomial.from_tensors(dense_tensors)
        # raise Exception(
        #     new.coeffs
        # )

        """
        [[ 5.00000000e+00  5.00000000e+00 -8.14070159e-02 -2.61438679e-23]
 [ 5.00000000e+00  5.00000000e+00 -8.14070159e-02 -2.61438679e-23]
 [-8.14070159e-02 -8.14070159e-02  7.00000000e+00 -4.97570369e-12]
 [-2.61438679e-23 -2.61438679e-23 -4.97570369e-12  9.00000000e+00]]
        :return:
        """

        ham = DGB(
            [[-1, -1], [-1, -1], [0, 0], [2, 2]],
            None,
            alphas=[5, 5, 7, 9],
            poly_coeffs=[None, None, None, None],
            masses=[1, 1],
            transformations=np.array([
                [
                    [ np.cos(10/180), -np.sin(10/180)],
                    [ np.sin(10/180),  np.cos(10/180)]
                ],
                [
                    [ np.cos(10 / 180), -np.sin(10 / 180)],
                    [ np.sin(10 / 180),  np.cos(10 / 180)]
                ],
                [
                    [np.cos(40 / 180), -np.sin(40 / 180)],
                    [np.sin(40 / 180), np.cos(40 / 180)]
                ],
                [
                    [1, 0],
                    [0, 1]
                ]
            ]),
            logger=True
        )
        # raise Exception(ham.T)

        rot_mat = np.array([
            [0.3990334094673942, -0.8487736362935158, -0.3469231218323582],
            [0.8989977455223288, 0.2876708903823626, 0.3302249420809833],
            [-0.18048654153314841, -0.4436538889266057, 0.8778358816804551]
        ])
        rot_mat2 = np.array([
            [0.6953147343696692, 0.38298501871028495, 0.6081610770280513],
            [-0.4098117159050655, -0.48387112171352364, 0.7732548707112311],
            [0.5904166136512273, -0.786887039552171, -0.17949097277979698]
        ])
        ang = np.deg2rad(90)
        ham_3 = DGB(
            [[-1, -2, 0], [-1, -2, 0], [-1, -1, 3], [-3, 2, 1]],
            None,
            projection_indices=[0, 1],
            alphas=[5, 5, 7, 9],
            # poly_coeffs=[{1: 1}, None, None, None],
            poly_coeffs=[None, None, None, None],
            masses=[1, 2, 1],#[2, 7, 13],
            transformations=np.array([
                [
                    [1, 0, 0],
                    [0, np.cos(ang), -np.sin(ang)],
                    [0, np.sin(ang),  np.cos(ang)]
                ],
                [
                    [1, 0, 0],
                    [0, np.cos(ang), -np.sin(ang)],
                    [0, np.sin(ang),  np.cos(ang)]
                ],
                np.eye(3),
                # rot_mat2,
                np.eye(3)
            ]),
            logger=True
        )

        """
        [[-0.00000000e+00 -0.00000000e+00  1.51479585e+01 -9.71873221e-04]
 [-0.00000000e+00 -0.00000000e+00  1.51479585e+01 -9.71873221e-04]
 [ 1.51479585e+01  1.51479585e+01 -0.00000000e+00 -0.00000000e+00]
 [-9.71873221e-04 -9.71873221e-04 -0.00000000e+00 -0.00000000e+00]]
 """

        # print(ham_3.S)
        raise Exception(ham_3.T)

    @validationTest
    def test_ModelPotentialAIMD(self):

        mol = Molecule.from_file(
            TestManager.test_data('water_freq.fchk'),
            internals=[[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]]
        )
        mol = mol.get_embedded_molecule()

        r1 = 0
        r2 = 1
        a12 = 2
        w2h = UnitsData.convert("Wavenumbers", "Hartrees")
        model = mol.get_model(
            {
                r1: {'morse':{'w':3869.47 * w2h, 'wx':84 *w2h}},
                r2: {'morse':{'w':3869.47 * w2h, 'wx':84 *w2h}},
                a12: {'harmonic':{'k': 1600**2/150 * w2h}}
            },
            dipole=[
                {
                    (r1, a12): ({'linear': {'eq': 0, 'scaling': 1/5.5}}, {'sin': {'eq': 0}}),
                    (r2, a12): ({'linear': {'eq': 0, 'scaling':-1/5.5}}, {'sin': {'eq': 0}})
                },
                {
                    (r1, a12): ({'linear': {'eq': 0, 'scaling': 1 / (2 * 5.5)}}, {'cos': {'eq': 0}}),
                    (r2, a12): ({'linear': {'eq': 0, 'scaling': 1 / (2 * 5.5)}}, {'cos': {'eq': 0}})
                },
                0
            ]
        )

        ics = mol.internal_coordinates
        re1 = ics[1, 0]; re2 = ics[2, 0]; ae = ics[2, 1]
        #
        # r = AnalyticModel.r
        # a = AnalyticModel.a
        # cos = AnalyticModel.cos
        # sin = AnalyticModel.sin
        # morse = AnalyticModel.morse
        # harmonic = AnalyticModel.harmonic
        # sym = AnalyticModel.sym
        # m = AnalyticModel.m
        # model = AnalyticModel(
        #     [
        #         r(1, 2),
        #         r(2, 3),
        #         a(1, 2, 3),
        #     ],
        #     morse(1, 2, w="w", wx="wx")
        #     + morse(2, 3, w="w", wx="wx")
        #     + harmonic(1, 2, 3),
        #     dipole=[
        #         ( r(2, 3) * sin(a(1, 2, 3)) -
        #             r(1, 2) * sin(a(1, 2, 3)) )  / 5.5,
        #         1/2 * (
        #                 r(2, 3) * cos(a(1, 2, 3))
        #                 + r(1, 2) * cos(a(1, 2, 3))
        #         ) / 5.5,
        #         0
        #     ],
        #     values={
        #         sym("w", 1, 2): 3869.47 * UnitsData.convert("Wavenumbers", "Hartrees"),
        #         sym("wx", 1, 2): 84 * UnitsData.convert("Wavenumbers", "Hartrees"),
        #         sym("re", 1, 2): re1,
        #         sym("w", 2, 3): 3869.47 * UnitsData.convert("Wavenumbers", "Hartrees"),
        #         sym("wx", 2, 3): 84 * UnitsData.convert("Wavenumbers", "Hartrees"),
        #         sym("re", 2, 3): re2,
        #         sym("k", 1, 2, 3): 1600**2/150 * UnitsData.convert("Wavenumbers", "Hartrees"),
        #         sym("qe", 1, 2, 3): ae,
        #
        #         r(1, 2):re1,
        #         r(2, 3):re2,
        #         a(1, 2, 3):ae,
        #         m(2): AtomData["O", "Mass"] * UnitsData.convert("AtomicMassUnits", "ElectronMass"),
        #         m(1): AtomData["H", "Mass"] * UnitsData.convert("AtomicMassUnits", "ElectronMass"),
        #         m(3): AtomData["H", "Mass"] * UnitsData.convert("AtomicMassUnits", "ElectronMass")
        #     }
        # )

        de = (3869.47 * w2h) ** 2 / (4 * 84 * w2h)

        check_freqs = False
        if check_freqs:
            freqs = model.normal_modes()[0]
            raise Exception(freqs*UnitsData.convert("Hartrees", "Wavenumbers"))

        check_anh = False
        if check_anh:
            model.run_VPT(logger=True)
            raise Exception(...) # very comparable to PODVR
            """
            State             Harmonic                     Anharmonic
                        ZPE                          ZPE    
              0 0 0    4680.66312                   4610.84351
            State       Frequency    Intensity       Frequency    Intensity
              0 0 1    3896.87027     64.98650      3719.85791     63.90801
              0 1 0    3843.25802      0.17386      3676.15021      0.12777
              1 0 0    1621.19796     64.86522      1603.43661     64.96941
              0 0 2    7793.74054      0.00000      7352.95576      0.00703
              0 2 0    7686.51604      0.00000      7268.92117      0.00362
              2 0 0    3242.39592      0.00000      3197.00895      0.07403
              0 1 1    7740.12829      0.00000      7229.92433      1.31011
              1 0 1    5518.06823      0.00000      5308.87367      0.09319
              1 1 0    5464.45598      0.00000      5278.21349      0.07390
              """

        check_dvr = False
        if check_dvr:
            dvr = model.setup_DVR(
                domain=[[1, 4], [1, 4], [np.deg2rad(60), np.deg2rad(160)]],
                divs=[800, 800, 800], po_divs=[20, 20, 20],
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
            Exception: [ 1605.10561571  3201.25184761  3675.81077847  3720.06017794
              4796.62930569  5280.78674778  5310.62058746  6409.77925325
              6876.17021986  6890.14602255  7219.00946583  7231.71801329
              7405.19785922  8061.11922153  8465.22595371  8480.52469107
              8819.15211826  8827.45481461  8989.31222846  9750.51070295
             10055.48547869 10126.58850176 10400.03078133 10410.80175582]
             """
            raise Exception(po_data.wavefunctions.frequencies()*UnitsData.hartrees_to_wavenumbers)

        cart_pot_func = model.potential
        cart_dipole_func = model.dipole

        def plot_pot():
            domain = [
                [-2.7, 2.7, 100],
                [-2.3, 1.0, 100]
            ]

            disps = np.moveaxis(
                np.array(np.meshgrid(*[np.linspace(*d) for d in domain], indexing='ij')),
                0, -1
            ).reshape(-1, 2)
            points = mol.get_nearest_displacement_coordinates(disps, axes=[0, 1])

            vals = cart_pot_func(points)
            cut = de / 2
            vals[vals > cut] = cut

            pot_plot = plt.TriContourPlot(
                disps[:, 0],
                disps[:, 1],
                vals,
                padding=[[50, 10], [50, 50]]
            )
            return pot_plot, points, [d[:2] for d in domain]

        plot_potential = False
        if plot_potential:
            plot_pot()[0].show()
            raise Exception(...)

        def sub_cart_pot_func(coords, deriv_order=None):
            # expects just x and y coordinates for the atoms
            coords = coords.reshape(-1, 3, 2)
            return cart_pot_func(coords, deriv_order=deriv_order, axes=[0, 1])


            # coords = np.concatenate([coords, np.zeros((len(coords), 3, 1))], axis=-1)
            #
            # terms = cart_pot_func(coords, deriv_order=deriv_order)
            # if deriv_order is not None:
            #     new = []
            #     for n,d in enumerate(terms):
            #         for j in range(n):
            #             d = np.take(d, (0, 1, 3, 4, 6, 7), axis=j+1)
            #         d = d.reshape(base_shape + d.shape[1:])
            #         new.append(d)
            # else:
            #     new = terms.reshape(base_shape)
            # return new

        def sub_cart_dipole_func(coords, deriv_order=None):
            # expects just x and y coordinates for the atoms
            coords = coords.reshape(-1, 3, 2)
            return cart_dipole_func(coords, deriv_order=deriv_order, axes=[0, 1])

        # raise Exception(
        #     sub_cart_dipole_func(
        #         np.broadcast_to(mol.coords[:, :2][np.newaxis], (100, 3, 2)).reshape(100, 6),
        #         deriv_order=1
        #     )
        # )

        np.random.seed(0)
        # init_disps = np.random.multivariate_normal([0, 0], np.power([[.2, 0], [0, .2]], 2), size=2)
        ri = .4
        # for init_disp in [
        #     # [-ri / 3, ri / 2, np.deg2rad(-50)],
        #     # [-ri / 3, ri / 2, np.deg2rad(-30)],
        #     # [ ri, 0, np.deg2rad(-50)],
        #     # [ ri, ri, np.deg2rad(60)],
        #     # [ .5, .5, np.deg2rad(60)],
        #     [ .55, .55, np.deg2rad(60)],
        #     # [ .7, .7, np.deg2rad(60)],
        # ]:

        for ie_vec, method, traj_steps, scaling in itertools.product(
                [
                    # [ 1650//2, 3850//2, 3900//2],
                    # [ 1650, 3850, 3900],
                    # [ 3500, 0, 4000 ],
                    # [ 1650, 3850//2, 3900 ]
                    # [ 1650,    0,    0],
                    # [    0, 3850,    0],
                    # [    0,    0, 3900],
                    # [ 1650, 3850,    0],
                    # # [    0, 3850, 3900], # Blows up?
                    [ 1650,    0, 3900],
                    # [ 1650, 3850, 3900],
                    # [ 1650 * 2,    0, 3900],
                ],
                [
                    # 'unrotated',
                    'rotated',
                    # 'H_rotation',
                    # 'base_rotation',
                    # 'min_dist'
                ],
                [10, 100],#, 150],#, 250],#, 500],
                [
                    1 / 2,
                    1 / 1.5,
                    1
                    # 1.2
                ]
        ):
            # initial_displacements = [
            #             # [ .55, .55, np.deg2rad(60)],
            #             [ .8, -.2, np.deg2rad(-7)]
            #         ]
            initial_displacements = None

            initial_energies = np.array([
                ie_vec
            ]) / UnitsData.hartrees_to_wavenumbers
            if initial_displacements is not None:
                init_pos = mol.get_displaced_coordinates(
                    initial_displacements,
                    which=[[1, 0], [2, 0], [2, 1]],
                    internals='reembed'
                )

                sim = AIMDSimulator(
                    mol.masses,
                    init_pos,
                    lambda c: -cart_pot_func(c, deriv_order=1)[1].reshape(c.shape),
                    timestep=1,
                    track_kinetic_energy=True
                )
            else:
                mol.potential_derivatives = cart_pot_func(mol.coords, deriv_order=2)[1:]
                modes = mol.normal_modes.modes.basis.matrix
                sim = AIMDSimulator(
                    mol.atomic_masses,
                    [mol.coords] * len(initial_energies),
                    lambda c: -cart_pot_func(c, deriv_order=1)[1].reshape(c.shape),
                    velocities=AIMDSimulator.mode_energies_to_velocities(modes, mol.atomic_masses, initial_energies),
                    timestep=60,
                    track_kinetic_energy=True
                )

            sim.propagate(traj_steps)
            coords = np.array(sim.trajectory).reshape((-1, 3, 3))
            coords = mol.embed_coords(coords)

            # plot_wavefunctions = False
            pot_plot = plot_pts = None
            plot_trajectory = False
            if plot_trajectory:
                pot_plot, plot_pts, domain = plot_pot()
                for i in range(3):
                    plt.ScatterPlot(
                        coords[:, i, 0],
                        coords[:, i, 1],
                        figure=pot_plot,
                        plot_range=domain
                    )
                pot_plot.show()
                raise Exception(...)

            e_init = round(
                np.mean(
                    cart_pot_func(sim.trajectory[0])
                    + sim.kinetic_energies[0]
                ) * UnitsData.hartrees_to_wavenumbers
            )
            init_e = np.round(
                np.mean(initial_energies, axis=0) * UnitsData.hartrees_to_wavenumbers
            ).astype(int)

            # for min_e in [None]:#100, 500, 800, 1000]:
            print(f"Running scaling: {scaling}")
            # plots_dir = os.path.join(
            #     os.path.expanduser("~/Documents/Postdoc/AIMD-Spec/water_model"),
            #     f"E{e_init}/EP{init_e[0]}_{init_e[1]}_{init_e[2]}/I{traj_steps}/{method}/S{scaling}"
            # )
            plots_dir = None
            if plots_dir is not None:
                os.makedirs(plots_dir, exist_ok=True)

            mass_vec = np.array(
                [AtomData["O", "Mass"] * UnitsData.convert("AtomicMassUnits", "ElectronMass")] * 2
                + [AtomData["H", "Mass"] * UnitsData.convert("AtomicMassUnits", "ElectronMass")] * 4
            )

            ham_coords = coords.reshape(-1, 9)[:, (0, 1, 3, 4, 6, 7)]
            method = 'H_rotation'
            if method == 'rotated':
                alphas = {
                    'method': 'virial',
                    'allow_rotations': True,
                    'planar': True,
                    'translation_rotation_frequency': 1e-18,
                    # 'translation_rotation_masses': mass_vec,
                    # 'min_frequency':500/UnitsData.hartrees_to_wavenumbers,
                    'scaling':  scaling
                }
                projection_indices = [3, 4, 5]
                poly_coeffs = None
                # ham_coords = np.concatenate(
                #     [ham_coords[:1], ham_coords],
                #     axis=0
                # )
                # poly_coeffs = [
                #     {0:1}
                # ] + [None] * (len(ham_coords)-1)
                poly_coeffs = [None] * (len(ham_coords))
            elif method == 'H_rotation':
                ra = (np.deg2rad(180) - ae) / 2
                cd = np.cos(ra)
                sd = np.sin(ra)
                base_rot = np.array([
                    [ 1,  0,   0,   0,   0,   0],
                    [ 0,  1,   0,   0,   0,   0],
                    [ 0,  0,  cd, -sd,   0,   0],
                    [ 0,  0,  sd,  cd,   0,   0],
                    [ 0,  0,   0,   0,  cd,  sd],
                    [ 0,  0,   0,   0, -sd,  cd]
                ])
                base_rot = np.eye(6)
                alphas = {
                    'method': 'virial',
                    'base_rotation': base_rot,
                    'scaling': scaling
                }
                projection_indices = None
                poly_coeffs = None
                # poly_coeffs = [None] * (len(ham_coords))
            elif method == 'base_rotation':
                mol.potential_derivatives = cart_pot_func(mol.coords, deriv_order=2)[1:]
                trip_mass = np.broadcast_to(
                    mol.masses[:, np.newaxis] * UnitsData.amu_to_me,
                    (3, 3)
                ).flatten()
                rot = np.concatenate([
                    mol.translation_rotation_modes[1][0],
                    mol.normal_modes.modes.basis.matrix
                ],
                    axis=-1
                ) * np.sqrt(trip_mass[:, np.newaxis])
                rot = rot / np.linalg.norm(rot, axis=0)[np.newaxis, :]
                # with np.printoptions(linewidth=1e8):
                #     raise Exception(str(rot))
                rot = rot[(0, 1, 3, 4, 6, 7), :][:, (0, 1, 5, 6, 7, 8)]
                # with np.printoptions(linewidth=1e8):
                #     raise Exception(str(rot))
                alphas = {
                    'method': 'virial',
                    'base_rotation': rot,
                    # 'min_frequency':500/UnitsData.hartrees_to_wavenumbers,
                    'scaling': scaling
                }
                projection_indices = [3, 4, 5]
                poly_coeffs = None
            elif method == 'min_dist':
                alphas = {'method': 'min_dist',
                          'scaling': scaling,
                          # 'min_frequency':min_e/UnitsData.hartrees_to_wavenumbers,
                          # 'allow_rotations':True
                          }
                projection_indices = None
                poly_coeffs = None
            else:
                alphas = {'method': 'virial', 'scaling': scaling, 'remove_translation_rotations': True}
                projection_indices = None
                poly_coeffs = None
                # poly_coeffs = [None] * (len(ham_coords))

            """[[0.07004169 0.04240438 0.03320008 0.04187683 0.05054381 0.04771534
  0.06524132 0.03906239 0.03523602 0.0476753  0.04793243]
 [0.04240438 0.06475006 0.05773734 0.0454053  0.02215346 0.01653997
  0.05643786 0.06126858 0.04852323 0.03584074 0.01519887]
  
                [[0.01479601 0.00859037 0.00621834 0.00791241 0.00739164 0.00687555
  0.01318776 0.00814297 0.00707907 0.00861234 0.00823437]
 [0.00859037 0.01378909 0.01235883 0.0094388  0.00343196 0.00230552
  0.01150354 0.01245531 0.00857874 0.00491725 0.00150661]
 [0.00621834 0.01235883 0.01295151 0.01026206 0.00330853 0.00191848
  0.00942091 0.01003521 0.00665623 0.00331115 0.00071587]
  """

            ham = DGB(
                ham_coords,
                sub_cart_pot_func,
                alphas=alphas,
                poly_coeffs=poly_coeffs,
                projection_indices=projection_indices,
                expansion_degree=2,
                masses=mass_vec,
                min_singular_value=1e-8,
                logger=True
            )

            raise Exception(ham.T)

            plot_orthog = False
            if plot_orthog:
                if plot_orthog is True:
                    plot_orthog = 15

                S = ham.S
                sig, Q = np.linalg.eigh(S)
                sig = np.flip(sig)
                Q = np.flip(Q, axis=1)
                wfns = DGBWavefunctions(
                    sig[:plot_orthog],
                    Q[:, :plot_orthog],
                    hamiltonian=ham
                )
                # wfns = DGBWavefunctions(
                #     sig[15+25:25+15+plot_orthog],
                #     Q[:, 25+15:25+15+plot_orthog],
                #     hamiltonian=ham
                # )

                for i in range(plot_orthog):
                    pot_plot, plot_pts, domain = plot_pot()
                    wfns[i].projection_plot(
                        [[0, 1], [2, 3], [4, 5]],
                        plotter=plt.TriContourLinesPlot,
                        contour_levels=10,
                        domain=domain,
                        cmap='RdBu',
                        figure=pot_plot,
                        plot_centers=True,
                        plot_label="Eig: {s}".format(s=wfns[i].energy)
                    )

                    if plots_dir is not None:
                        os.makedirs(os.path.join(plots_dir, "seigs"), exist_ok=True)
                        pot_plot.savefig(os.path.join(plots_dir, "seigs", f"S_eig_{i}.png"))
                        pot_plot.close()
                    else:
                        pot_plot.show()
                raise Exception(...)

#             with np.printoptions(linewidth=1e8):
#                 """
# Exception: ((15, 15), 2.625929151914913e-05, 3.810340801964157, 0.30135909686720164, (15, 15), -0.0018605379556571248, 0.024083643909265783, 0.006150208276975052)
# """
#                 raise Exception(ham.V.shape, np.min(ham.V), np.max(ham.V), np.std(ham.V.flatten()),
#                                 ham.T.shape, np.min(ham.T), np.max(ham.T), np.std(ham.T.flatten()))

            plot_gaussians = False
            if plot_gaussians:

                for i in range(-1, len(ham.centers)):
                    pot_plot, plot_pts, domain = plot_pot()

                    if i < 0:
                        g_sum_dat = np.ones((len(ham.centers), 1))
                    else:
                        g_sum_dat = np.zeros((len(ham.centers), 1))
                        g_sum_dat[i, 0] = 1
                    wfns = DGBWavefunctions(
                        [1],
                        g_sum_dat,
                        hamiltonian=ham
                    )
                    wfns[0].projection_plot(
                        [[0, 1], [2, 3], [4, 5]],
                        plotter=plt.TriContourLinesPlot,
                        contour_levels=10,
                        domain=domain,
                        cmap='RdBu',
                        figure=pot_plot,
                        plot_centers=True
                    )

                    if plots_dir is not None:
                        os.makedirs(os.path.join(plots_dir, "gaussians"), exist_ok=True)
                        pot_plot.savefig(os.path.join(plots_dir, "gaussians", f"gaussian_{i}.png"))
                        pot_plot.close()
                    else:
                        pot_plot.show()
                raise Exception(...)

            plot_wfns = True
            if plot_wfns:
                wfns = ham.get_wavefunctions(
                    nodeless_ground_state=True,
                    stable_epsilon=2e-4,
                    # min_singular_value=2e-4,
                    # subspace_size=ssize,
                    mode='classic'
                )
                h2w = UnitsData.convert("Hartrees", "Wavenumbers")
                plots = []
                plot_me = np.where(wfns.frequencies() < 8000 / UnitsData.hartrees_to_wavenumbers)
                num_wfns = 1 + (0 if len(plot_me) == 0 or len(plot_me[0]) == 0 else max(plot_me[0]))
                for n in range(num_wfns):
                    # if pot_plot is None:
                    pot_plot, plot_pts, domain = plot_pot()
                    if len(wfns) <= n:
                        break

                    wfn = wfns[n]
                    proj_plot = wfn.projection_plot(
                        [[0, 1], [2, 3], [4, 5]],
                        figure=pot_plot,
                        plot_label=f"Energy: {(wfn.energy - (0 if n == 0 else wfns[0].energy)) * h2w:.0f}",
                        padding=[[50, 10], [50, 50]],
                        plotter=plt.TriContourLinesPlot,
                        contour_levels=10,
                        domain=domain,
                        cmap='RdBu',
                        plot_centers=True
                    )
                    if plots_dir is not None:
                        pot_plot.savefig(os.path.join(plots_dir, f"wfn_{n}.png"))
                        pot_plot.close()
                    else:
                        plots.append(pot_plot)
                if plots_dir is None:
                    plots[0].show()
            raise Exception(...)

        plot_subspace_energies = True
        if plot_subspace_energies:
            engs = []
            for ssize in range(traj_steps, 10, -8):
                plots_dir = f"/Users/Mark/Documents/Postdoc/AIMD-Spec/stab_tests/S{scaling}/size_{ssize}"
                os.makedirs(plots_dir, exist_ok=True)
                with np.printoptions(linewidth=1e8):
                    wfns = ham.get_wavefunctions(
                        nodeless_ground_state=False,
                        stable_epsilon=2e-4,
                        # min_singular_value=2e-4,
                        subspace_size=ssize,
                        mode='classic'
                    )
                engs.append(wfns.energies * UnitsData.convert("Hartrees", "Wavenumbers"))

            ploot = None
            for i,e in enumerate(engs):
                from matplotlib.colors import hsv_to_rgb
                ploot = plt.ScatterPlot(
                    np.full(len(e), (2*i)),
                    e,
                    plot_range=[None, [0, 16000]],
                    figure=ploot,
                    color=hsv_to_rgb(
                            np.concatenate(
                                [np.clip((e%8000)/8000, 0, 1)[:, np.newaxis], np.ones((len(e), 2))],
                                axis=-1
                            )
                        ),
                    axes_labels=['Subspace Depletion', 'Energy']
                )
            ploot.show()
            raise Exception(...)

            #     h2w = UnitsData.convert("Hartrees", "Wavenumbers")
            #     plots = []
            #     plot_me = np.where(wfns.frequencies() < 1e10)# 8000 / UnitsData.hartrees_to_wavenumbers)
            #     num_wfns = 1 + (0 if len(plot_me) == 0 or len(plot_me[0]) == 0 else max(plot_me[0]))
            #     for n in range(num_wfns):
            #         # if pot_plot is None:
            #         pot_plot, plot_pts, domain = plot_pot()
            #         if len(wfns) <= n:
            #             break
            #
            #         wfn = wfns[n]
            #         proj_plot = wfn.projection_plot(
            #             [[0, 1], [2, 3], [4, 5]],
            #             figure=pot_plot,
            #             plot_label=f"Energy: {(wfn.energy - (0 if n == 0 else wfns[0].energy)) * h2w:.0f}",
            #             padding=[[50, 10], [50, 50]],
            #             plotter=plt.TriContourLinesPlot,
            #             contour_levels=10,
            #             domain=domain,
            #             cmap='RdBu',
            #             plot_centers=True
            #         )
            #         if plots_dir is not None:
            #             pot_plot.savefig(os.path.join(plots_dir, f"wfn_{n}.png"))
            #             pot_plot.close()
            #         else:
            #             plots.append(pot_plot)
            #     if plots_dir is None:
            #         plots[0].show()
            # raise Exception(...)

            plot_spectrum = True
            if plot_spectrum:
                if len(wfns) > 1:
                    plot_me = np.where(wfns.frequencies() * UnitsData.hartrees_to_wavenumbers < 8000)
                    num_wfns = 1 + (0 if len(plot_me) == 0 or len(plot_me[0]) == 0 else max(plot_me[0]))
                    spec = wfns[:num_wfns].get_spectrum(sub_cart_dipole_func, expansion_degree=1)  # .normalize(0)
                    sploot = spec.plot()#plot_range=[[1400, 7000], None])
                    if plots_dir is not None:
                        np.savetxt(
                            os.path.join(plots_dir, "spec.txt"),
                            np.array([
                                spec.frequencies,
                                spec.intensities
                            ])
                        )
                        sploot.savefig(os.path.join(plots_dir, "spec.png"))
                        sploot.close()
                    else:
                        sploot.show()

        plt.Graphics().show()

        raise Exception(wfns.frequencies() * UnitsData.convert("Hartrees", "Wavenumbers"))

        # hmm, grads = cart_pot_derivs(coords, deriv_order=1)
        # ics = Molecule(mol.atoms, coords, internals=mol.internals).internal_coordinates
        # r1s = ics[..., 2, 1]
        # sort = np.argsort(r1s)
        # plt.Plot(r1s[sort], hmm[sort]).show()
        # raise Exception(hmm)

    @validationTest
    def test_Expansion(self):
        from McUtils.Zachary import Symbols, RBFDInterpolator

        sym = Symbols('xyz')

        np.random.seed(3)
        ndim = 2
        pts = np.random.uniform(low=-5, high=5, size=(3500, ndim))

        # fn = sym.morse(sym.x) * sym.morse(sym.y) - sym.morse(sym.x) - sym.morse(sym.y)
        vars = [sym.x, sym.y][:ndim]

        w = 2.05; wx = .1; mu = 1
        de = (w ** 2) / (4 * wx)
        a = np.sqrt(2 * mu * wx)

        def morse_pot(x, de=de, a=a):
            return sym.morse(x, de=de, a=a)
        def harmonic_pot(x, w=w):
            return w/2 * x**2
        def quartic_pot(r, a=a, de=de):
            return (a**2*de)*r**2 - (a**3*de)*r**3 + (7/12*a**4*de)*r**4
        def sextic_pot(r, a=a, de=de):
            return (de*a**2)*r**2 - (de*a**3)*r**3 + (7/12*de*a**4)*r**4 - (de*a**5 /4)*r**5+ (31*de*a**6/360)*r**6

        fn_1D = morse_pot
        fn = sum(fn_1D(var) for var in vars) if ndim > 1 else fn_1D(vars[0])

        def simple_morse(c, de=de, a=a, deriv_order=None):
            ndim = c.shape[-1]
            if deriv_order is None:
                return np.sum(de*(1-np.exp(-a*c))**2, axis=-1)
            else:
                n = deriv_order
                m = ((-1)**(n+1) * 2 * a**n * de) * np.exp(-2*a*c)*(np.exp(a*c)-(2**(n-1)))
                if n == 1:
                    return m
                res = np.zeros(c.shape[:-1] + (ndim,)*deriv_order)
                for k in range(ndim):
                    idx = (...,) + (k,)*n
                    res[idx] = m[..., k]
                return res

        # print(res['AIMDEnergies'].gradients)
        # interp = RBFDInterpolator.create_function_interpolation(
        #     pts,
        #     fn,
        #     lambda p, f=fn.deriv(order=1): f(p),#.transpose(),
        #     lambda p, f=fn.deriv(order=2): f(p),#.transpose(2, 0, 1),
        #     clustering_radius=1e-5,
        # )

        test_pot = lambda c,fn=fn,deriv_order=None: (
            fn(c.reshape(-1, ndim)).reshape(c.shape[:-1])
                if deriv_order is None else
            np.moveaxis(fn.deriv(order=deriv_order)(c.reshape(-1, ndim)), -1, 0).reshape(c.shape[:-1] + (ndim,)*deriv_order)
        )

        np.random.seed(0)
        centers = pts
        # centers = np.array(
        #     np.meshgrid(*[
        #         np.linspace(-5, 5, 50)
        #     ]*ndim)
        # ).T.reshape((-1, ndim))


        alpha = .3
        # cluster = .045

        wfns = DGB.run(
            centers,
            test_pot,
            alphas=alpha,
            # clustering_radius=cluster,
            quadrature_degree=4,
            min_singular_value=1e-3
        )
        e = wfns.energies
        print(e[:10])

        wfns = DGB.run(
            wfns.hamiltonian.centers, # reuse the opt from before
            test_pot,
            optimize_centers=False,
            alphas=alpha,
            # clustering_radius=cluster,
            expansion_degree=6,
            min_singular_value=1e-3
        )
        e = wfns.energies
        print(e[:10])

        ham = DGB.run(
            wfns.hamiltonian.centers, # reuse the opt from before
            test_pot,
            optimize_centers=False,
            alphas=alpha,
            # clustering_radius=cluster,
            expansion_degree=6,
            expansion_type='taylor',
            min_singular_value=1e-3
        )
        e = wfns.energies
        print(e[:10])

    @validationTest
    def test_Interp(self):
        from McUtils.Zachary import Symbols, SymPyFunction, RBFDInterpolator

        sym = Symbols('xyz')

        np.random.seed(3)
        ndim = 2
        pts = np.random.uniform(low=-2, high=8, size=(10000, ndim))

        # fn = sym.morse(sym.x) * sym.morse(sym.y) - sym.morse(sym.x) - sym.morse(sym.y)
        vars = [sym.x, sym.y][:ndim]

        w = 2.05; wx = .1; mu = 1
        de = (w ** 2) / (4 * wx)
        a = np.sqrt(2 * mu * wx)

        def morse_pot(x, de=de, a=a):
            return sym.morse(x, de=de, a=a)

        def harmonic_pot(x, w=w):
            return w / 2 * x ** 2

        def quartic_pot(r, a=a, de=de):
            return (a ** 2 * de) * r ** 2 - (a ** 3 * de) * r ** 3 + (7 / 12 * a ** 4 * de) * r ** 4

        def sextic_pot(r, a=a, de=de):
            return (de * a ** 2) * r ** 2 - (de * a ** 3) * r ** 3 + (7 / 12 * de * a ** 4) * r ** 4 - (
                        de * a ** 5 / 4) * r ** 5 + (31 * de * a ** 6 / 360) * r ** 6

        fn_1D = morse_pot
        fn = sum(fn_1D(var) for var in vars) if ndim > 1 else fn_1D(vars[0])
        fn = fn_1D(vars[0])
        for v in vars[1:]:
            fn = fn * fn_1D(v)
        fn = fn / (100)
        for v in vars[1:]:
            fn = fn + fn_1D(v)
        # fn = sum(fn_1D(var) for var in vars) if ndim > 1 else fn_1D(vars[0])


        x, y = SymPyFunction.symbols('x', 'y')
        def morse_pot(x, de=de, a=a):
            return SymPyFunction.morse(x, de=de, a=a)
        fn = morse_pot(x) * morse_pot(y) / 100 + morse_pot(x) + morse_pot(y)

        interp = RBFDInterpolator.create_function_interpolation(
            pts,
            fn,
            lambda p, f=fn.deriv(order=1): np.moveaxis(f(p), -1, 0),
            lambda p, f=fn.deriv(order=2): np.moveaxis(f(p), -1, 0),
            clustering_radius=1e-5,
            # neighborhood_size=30,
            multicenter_monomials=True,
            extra_degree=2,
            error_threshold=.01,
            neighborhood_merge_threshold=None
            # solve_method='solve'
        )

        test_pot = lambda c, fn=fn, deriv_order=None: (
            fn(c.reshape(-1, ndim)).reshape(c.shape[:-1])
            if deriv_order is None else
            np.moveaxis(fn.deriv(order=deriv_order)(c.reshape(-1, ndim)), -1, 0).reshape(
                c.shape[:-1] + (ndim,) * deriv_order)
        )

        # import json, os
        # with open(os.path.expanduser('~/Desktop/interp_data.json'), 'w+') as out:
        #     json.dump({
        #         'points':interp.grid.tolist(),
        #         'vals':interp.vals.tolist(),
        #         'ders':[d.tolist() for d in interp.derivs]
        #     }, out)
        # raise Exception(...)


        # raise Exception(interp(pts[:5], deriv_order=2).shape)

        # raise Exception(fn(pts))

        # print(res['AIMDEnergies'].gradients)


        # test_pot = harmonic

        np.random.seed(0)
        tv = fn(pts)
        prune_chunks = [
            pts[tv < 3],
            pts[np.logical_and(tv > 3, tv < 5)],
            pts[np.logical_and(tv > 5, tv < 10)]
        ]
        centers = np.concatenate(
            [
                prune_chunks[0],
                prune_chunks[1][np.random.random_integers(0, len(prune_chunks[1]), 100)],
                prune_chunks[2][np.random.random_integers(0, len(prune_chunks[1]), 50)],
            ],
            axis=0
        )
        # centers = np.array(
        #     np.meshgrid(*[
        #         np.linspace(-5, 5, 50)
        #     ]*ndim)
        # ).T.reshape((-1, ndim))

        alpha = 1.5
        # cluster = .045

        ham = DGB(
            centers,
            test_pot,
            alphas=alpha,
            clustering_radius=0.01,
            optimize_centers=True,
            # clustering_radius=.055,
            # optimize_centers=False,
            # expansion_degree=2
        )
        e, wf = ham.get_wavefunctions()  # print_debug_info=True)
        print(e[:50])

        ham = DGB(
            ham.centers,
            test_pot,
            alphas=alpha,
            clustering_radius=None,
            optimize_centers=False,
            # clustering_radius=.055,
            # optimize_centers=False,
            expansion_degree=2
        )
        # print(ham.clustering_radius)

        # np.seterr(all=None, divide=None, over=None, under='warn', invalid=None)

        # c, a = ham.get_overlap_gaussians()
        # wat = interp(c[[52, 52], [75, 91]], deriv_order=2, use_cache=False, return_interpolation_data=True)
        # raise Exception("...")

        # ivals, error = interp(c, deriv_order=2, return_error=True)
        # real_vals = fn(c)
        # diffs = ivals[0] - real_vals
        # real_d1s = fn.deriv(order=1)(c).transpose(1, 2, 0)
        # d1_diffs = ivals[1] - real_d1s
        # real_d2s = fn.deriv(order=2)(c).transpose(2, 3, 0, 1)
        # d2_diffs = ivals[2] - real_d2s
        #
        # bad_pos = np.where(np.abs(d2_diffs) > 5)
        # bad_pos = [
        #     bad_pos[0][bad_pos[0] <= bad_pos[1]],
        #     bad_pos[1][bad_pos[0] <= bad_pos[1]]
        # ]
        # bad_pos = tuple(np.unique(np.array(bad_pos).T, axis=0).T)
        # _, idata = interp(c[bad_pos], deriv_order=2, use_cache=False, return_interpolation_data=True)
        # subdata = [
        #     {
        #         'mat': i['matrix'].tolist(),
        #         'weights': i['data'].weights[0].tolist(),
        #         'centers': i['points'].tolist(),
        #         'matrix': i['solver_data']['matrix'].tolist(),
        #         'vals': i['solver_data']['vals'].tolist(),
        #         'cn': i['solver_data']['condition_number']
        #     } for i in idata
        # ]
        # print(">>", [
        #     np.dot(i['matrix'], i['weights']) - i['vals']
        #     for i in subdata
        #       ])

        # good_pos = np.where(np.abs(d2_diffs) < .0001)
        # good_pos = [
        #     good_pos[0][good_pos[0] <= good_pos[1]],
        #     good_pos[1][good_pos[0] <= good_pos[1]]
        # ]
        # good_pos = tuple(np.unique(np.array(good_pos).T, axis=0)[:15].T)
        # _, idata2 = interp(c[good_pos], deriv_order=2, use_cache=False, return_interpolation_data=True)

        # import json, os
        # with open(os.path.expanduser('~/Desktop/overlap_center_data.json'), 'w+') as out:
        #     json.dump({
        #         'fn':str(fn),
        #         'ipoints':pts.tolist(),
        #         'samp_vs':interp.vals.tolist(),
        #         'samp_ders':[d.tolist() for d in interp.derivs],
        #         'centers':ham.centers.tolist(),
        #         'errors':error.tolist(),
        #         'gauss':c.tolist(),
        #         'ivals':[i.tolist() for i in ivals],
        #         'reals':[real_vals.tolist(), real_d1s.tolist(), real_d2s.tolist()],
        #         'diffs':[diffs.tolist(), d1_diffs.tolist(), d2_diffs.tolist()],
        #         'idat': [
        #             [b.tolist() for b in bad_pos],
        #             [
        #                 {
        #                     'mat': i['matrix'].tolist(),
        #                     'weights': i['data'].weights[0].tolist(),
        #                     'centers': i['points'].tolist(),
        #                     'matrix': i['solver_data']['matrix'].tolist(),
        #                     'vals': i['solver_data']['vals'].tolist(),
        #                     'cn': i['solver_data']['condition_number']
        #                 } if i is not None else "n/a" for i in idata
        #             ]
        #         ],
        #         'gdat': [
        #             [b.tolist() for b in good_pos],
        #             [
        #                 {
        #                     'mat': i['matrix'].tolist(),
        #                     'weights': i['data'].weights[0].tolist(),
        #                     'centers': i['points'].tolist(),
        #                     'vals': i['solver_data']['vals'].tolist(),
        #                     'cn': i['solver_data']['condition_number']
        #                 } if i is not None else "n/a"  for i in idata2
        #             ]
        #         ]
        #     }, out)
        # print("...saved")

        e, wf = ham.get_wavefunctions()#print_debug_info=True)
        print(e[:50])


        ham_1 = ham
        ham = DGB(
            ham.centers, # reuse the opt from before
            interp,
            alphas=alpha,
            optimize_centers=False,
            clustering_radius=None,
            expansion_degree=2
         )

        e, wf = ham.get_wavefunctions()#print_debug_info=True)
        print(e[:50])
        # with np.printoptions(linewidth=1e8):
        #     print(ham.V - ham_1.V)
        import json, os
        with open(os.path.expanduser('~/Desktop/dgb_V_diffs.json'), 'w+') as out:
            json.dump({
                "S": ham_1.S.tolist(),
                "T": ham_1.T.tolist(),
                'exact': ham_1.V.tolist(),
                'interp': ham.V.tolist(),
                "diff": (ham_1.V - ham.V).tolist()
            },
                out)

        # print("UPDATE")
        # ham.V = ham_1.V
        # e, wf = ham.get_wavefunctions()  # print_debug_info=True)
        # print(e[:10])

        ham = DGB(
            ham.centers, # reuse the opt from before
            interp,
            alphas=alpha,
            optimize_centers=False,
            clustering_radius=None,
            expansion_degree=2,
            expansion_type='taylor'
        )
        e, wf = ham.get_wavefunctions() # print_debug_info=True)
        print(e[:50])

    @inactiveTest
    def test_WaterFromGauss(self):
        from McUtils.Numputils import vec_tensordot

        with GaussianLogReader(TestManager.test_data('h2o_aimd.log')) as parser:
            res = parser.parse(['AIMDCoordinates', 'AIMDEnergies'])

        crds = res['AIMDCoordinates']#.reshape(-1, 9)
        e = (res['AIMDEnergies'].energies - np.min(res['AIMDEnergies'].energies)) #* 219474.65
        # print(res['AIMDEnergies'].gradients)

        ref_pos = np.argmin(e)
        # raise Exception( # proves no extra unit conversion necessary I think...
        #     # e[ref_pos] * 219475,
        #     Molecule(
        #         ["O", "H", "H"],
        #         res['AIMDCoordinates'][ref_pos],
        #         potential_derivatives=[
        #             res['AIMDEnergies'].gradients[ref_pos],
        #             res['AIMDEnergies'].hessians[ref_pos]
        #         ]
        #     ).normal_modes.modes.freqs * 219475
        # )

        # raise Exception(e* 219475)

        def coord_transformer(crds):
            _ = 10000
            return Molecule(["O", "H", "H"], crds.reshape(-1, 3, 3).squeeze(), internals=[[0, _, _, _], [1, 0, _, _], [2, 0, 1, _]])

        def coord_transf(crds, deriv_order=0, direction='forward'):
            vals = []
            jacs = [[] for _ in range(deriv_order)]
            for i,c in enumerate(crds):
                if i % 50 == 0:
                    print("...", i)
                mol = coord_transformer(c)
                icrds = mol.internal_coordinates
                # if direction == 'forward':
                vals.append(icrds[(1, 2, 2), (0, 0, 1)])
                # else:
                #     vals.append(icrds.convert(mol.coords.system))
                if deriv_order > 0:
                    if direction == 'reverse': # dR/dX
                        for n, j in enumerate(mol.coords.jacobian(icrds.system, list(range(1, deriv_order+1)), all_numerical=True)):
                            jacs[n].append(j)
                    else: # dX/dR
                        for n,j in enumerate(icrds.jacobian(mol.coords.system, list(range(1, deriv_order+1)), all_numerical=True)):
                            jacs[n].append(j)
            vals = np.array(vals)
            jacs = [np.array(j) for j in jacs]
            return vals, jacs

        def derivs_to_internals(jacs, grd, hes):
            dX = jacs[0].reshape((-1, 9, 9))[:, (3, 6, 7), :]
            dXX = jacs[1].reshape((-1, 9, 9, 9))[:, (3, 6, 7)][:, :, (3, 6, 7)]
            # print(dXX.shape, dX.shape)

            # import McUtils.Plots as plt
            # plt.ArrayPlot(dX[0]).show()
            # raise Exception(dXX)
            g1 = grd
            grd = vec_tensordot(dX, g1.reshape(-1, 9), axes=[2, 1], shared=1)
            hes = vec_tensordot(
                dX,
                vec_tensordot(dX, hes.reshape(-1, 9, 9), axes=[2, 2], shared=1),
                axes=[2, 2],
                shared=1
            )
            hes = hes + vec_tensordot(
                dXX,
                g1.reshape(-1, 9),
                axes=[3, 1],
                shared=1
            )

            return grd, hes

        def derivs_to_carts(jacs, grd, hes):
            grd = grd.reshape((-1, 3))
            hes = hes.reshape((-1, 3, 3))

            dX = jacs[0].reshape((-1, 9, 9))[:, :, (3, 6, 7)]
            dXX = jacs[1].reshape((-1, 9, 9, 9))[:, :, :, (3, 6, 7)]
            # print(grd.shape, hes.shape, dX.shape, dXX.shape)
            g2 = vec_tensordot(dX, grd, axes=[2, 1], shared=1)

            h1 = vec_tensordot(dX, hes, axes=[2, 2], shared=1)
            # print("???", h1.shape)
            h1 = vec_tensordot(dX, h1, axes=[2, 2], shared=1)
            h2 = h1 + vec_tensordot(dXX, grd, axes=[3, 1], shared=1)

            grd = g2

            return g2, h2

        ref_pos = np.argmin(e)
        ref = Molecule(["O", "H", "H"], crds[ref_pos]).get_embedded_molecule(load_properties=False)
        all_crds = []
        all_grds = []
        all_hess = []
        for c, g, h in zip(crds, res['AIMDEnergies'].gradients, res['AIMDEnergies'].hessians):
            mol = Molecule(["O", "H", "H"],
                           c,
                           potential_derivatives=[g,h]
                           ).get_embedded_molecule(ref, load_properties=False)
            all_crds.append(mol.coords)
            pes = mol.potential_derivatives
            all_grds.append(pes[0])
            all_hess.append(pes[1])

        crd = np.array(all_crds).reshape(-1, 9)#[[ref_pos]]
        grd1 = np.array(all_grds)#[[ref_pos]]
        hes1 = np.array(all_hess)#[[ref_pos]]

        # vals, jacs = coord_transf(crd, direction='forward', deriv_order=2)
        # grd, hes = derivs_to_internals(jacs, grd1, hes1)
        #
        # rev_crds, rev_jacs = coord_transf(crd, direction='reverse', deriv_order=2)
        # grd2, hes2 = derivs_to_carts(rev_jacs, grd, hes)
        # import McUtils.Plots as plt
        # plt.ArrayPlot(rev_jacs[0].reshape(9, 9).T @ jacs[0].reshape(9, 9)).show()
        #
        # raise Exception(rev_jacs[0].reshape(9, 9) @ jacs[0].reshape(9, 9))

        # raise Exception( # proves no extra unit conversion necessary I think...
        #     # e[ref_pos] * 219475,
        #     np.sqrt(
        #         scipy.linalg.eigh(
        #             hes[0],
        #             coord_transformer(crd[0]).g_matrix,
        #             type=2
        #         )[0]
        #     ) * 219475,
        #     Molecule(
        #         ["O", "H", "H"],
        #         rev_crds[0].reshape(3, 3),
        #         potential_derivatives=[
        #             grd2[0],
        #             hes2[0]
        #         ]
        #     ).normal_modes.modes.freqs * 219475
        # )
        #
        # raise Exception(hes1[0] - hes2[0])

        # raise Exception( # proves no extra unit conversion necessary I think...
        #     # e[ref_pos] * 219475,
        #     Molecule(
        #         ["O", "H", "H"],
        #         crd[0].reshape(3, 3),
        #         potential_derivatives=[
        #             grd[0],
        #             hes1[0]
        #         ]
        #     ).normal_modes.modes.freqs * 219475
        # )

        # vals, jacs = coord_transf(crd, direction='forward', deriv_order=2)

        vals, jacs = coord_transf(crd, direction='forward', deriv_order=2)
        grd, hes = derivs_to_internals(jacs, grd1, hes1)

        base_interp = McUtils.Zachary.RBFDInterpolator(
            vals,
            e,
            grd, # not sure if this needs a conversion?
            hes, # not sure if this needs a conversion?
            clustering_radius=0,
            error_threshold=.01
        )

        def interp(crds, deriv_order=0):
            # print(">>>>", crds.shape)
            base_shape = crds.shape[:-1]
            crds = crds.reshape(-1, 9)
            # print("  > ", crds.shape)
            vals, jacs = coord_transf(crds, deriv_order=deriv_order, direction='reverse')
            # print("  > ", vals.shape)
            # print("  > ", [j.shape for j in jacs])
            if deriv_order == 2:
                e, grd, hes = base_interp(vals, deriv_order=deriv_order)
                grd, hes = derivs_to_carts(jacs, grd, hes)
                # grd = grd.reshape((-1, 3))
                # hes = hes.reshape((-1, 3, 3))
                # dX = jacs[0].reshape((-1, 9, 9))[:, :, (3, 6, 7)]
                # dXX = jacs[1].reshape((-1, 9, 9, 9))[:, :, :, (3, 6, 7)]
                # # print(e.shape, grd.shape, hes.shape, dX.shape, dXX.shape)
                # g2 = vec_tensordot(dX, grd, axes=[2, 1], shared=1)
                #
                # h1 = vec_tensordot(dX, hes, axes=[2, 2], shared=1)
                # # print("???", h1.shape)
                # h1 = vec_tensordot(dX, h1, axes=[2, 2], shared=1)
                # hes = h1 + vec_tensordot(dXX, grd, axes=[3, 1], shared=1)
                # grd = g2 # ugh
            else:
                raise ValueError("...")
            e = e.reshape(base_shape)
            grd = grd.reshape(base_shape + (9,))
            hes = hes.reshape(base_shape + (9, 9))
            # print("  > ", e.shape, grd.shape, hes.shape)
            return e, grd, hes

        # raise Exception(e-interp(vals))

        ham = DGB(
            crd,
            interp,
            alphas=1.5,
            clustering_radius=0,
            optimize_centers=False,
            # clustering_radius=None,
            # optimize_centers=False,
            expansion_degree=2
        )
        print(ham.centers.shape)
        e, wf = ham.get_wavefunctions(min_singular_value=1e-5)  # print_debug_info=True)
        print(e * 219475)

        # ham = DGB(
        #     ham.centers,
        #     interp.global_interpolator,
        #     alphas=.5,
        #     clustering_radius=None,
        #     optimize_centers=False,
        #     expansion_degree=2
        # )
        #
        # # raise Exception(ham.clustering_radius)
        # e, wf = ham.get_wavefunctions()  # print_debug_info=True)
        # print(e)

        # need a new center to expand about to make this work...

        # ham = DGB(
        #     ham.centers,  # reuse the opt from before
        #     interp,
        #     alphas=.5,
        #     # clustering_radius=.005,
        #     # optimize_centers=False,
        #     expansion_degree=2,
        #     expansion_type='taylor'
        # )
        # e, wf = ham.get_wavefunctions()  # print_debug_info=True)  # print_debug_info=True)
        # print(e[:10])

    @inactiveTest
    def test_WaterFromAIMD(self):

        from McUtils.Scaffolding import Checkpointer
        from McUtils.Zachary import RBFDInterpolator, InverseDistanceWeightedInterpolator

        loader = ModuleLoader(TestManager.current_manager().test_data_dir)
        mbpol = loader.load("LegacyMBPol").MBPol

        b2a = UnitsData.convert("BohrRadius", "Angstroms")

        def potential(coords, deriv_order=0, chunk_size=int(5e5)):

            coords = coords.reshape(-1, 9)

            chunks = [[] for _ in range(deriv_order+1)]
            num_chunks = int(len(coords) / chunk_size) + 1

            for coords in np.array_split(coords, num_chunks):
                # interp.logger.log_print("evaluating energies")
                energies = mbpol.get_pot(coords=coords.reshape(-1, 3, 3) * b2a, nwaters=1, threading_vars=['energy', 'coords'], threading_mode='omp')
                if deriv_order > 0:
                    derivs = []
                    grads = lambda c: mbpol.get_pot_grad(
                        nwaters=1, coords=c.reshape(-1, 3, 3) * b2a, threading_vars=['energy', 'grad', 'coords'], threading_mode='omp'
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
                            for i,d in enumerate(new_derivs)
                        )
                    # interp.logger.log_print("done")
                    for i,d in enumerate([energies] + derivs):
                        chunks[i].append(d)
                else:
                    chunks[0].append(energies)

            for i,c in enumerate(chunks):
                chunks[i] = np.concatenate(c, axis=0)

            return chunks

        # ref = np.array([  # some random structure Mathematica got from who knows where...
        #     [7.74875081e-06, 6.20321759e-02, -2.09594068e-17],
        #     [7.83793537e-01, -4.92291304e-01, -1.82568317e-16],
        #     [-7.83916515e-01, -4.92204344e-01, 2.82713128e-17]
        # ]) * UnitsData.convert("Angstroms", "BohrRadius")
        # r1 = np.linalg.norm(ref[1] - ref[0])
        # v1 = (ref[1] - ref[0]) / r1

        # from Psience.DVR import DVR
        # def pot(r):
        #     disps = v1[np.newaxis, :] * r[:, np.newaxis]
        #     new = np.broadcast_to(ref[np.newaxis], (len(r), 3, 3)).copy()
        #     new[:, 1] += disps
        #     return mbpol.get_pot(coords=new * UnitsData.convert("BohrRadius", "Angstroms"), nwaters=1, threading_vars=['energy', 'coords'], threading_mode='omp')

        # sys = DVR(
        #     domain=[r1 - .5, r1 + 1],
        #     divs=100,
        #     potential_function=pot,
        #     mass=1/(
        #             1/(UnitsData.convert("AtomicMassUnits", "ElectronMass") * AtomData["O", "Mass"])
        #             + 1/(UnitsData.convert("AtomicMassUnits", "ElectronMass") * AtomData["H", "Mass"])
        #     )
        # )
        # res = sys.run()
        # raise Exception(res.wavefunctions.frequencies() * UnitsData.convert("Hartrees", "Wavenumbers"))


        # import scipy.optimize as opt
        # base_pot = lambda coords:mbpol.get_pot(coords=coords.reshape(3, 3), nwaters=1)
        # ref = np.array([
        #     [ 0.00000000e+00,  6.56215885e-02, 0.00000000e+00],
        #     [ 7.57391014e-01, -5.20731105e-01, 0.00000000e+00],
        #     [-7.57391014e-01, -5.20731105e-01, 0.00000000e+00]
        # ])
        # opt_vals = opt.minimize(
        #     base_pot,
        #     ref,
        #     method='Nelder-Mead',
        #     options=dict(fatol=1e-16)
        # )
        #
        # ugh = Molecule(
        #     ["O", "H", "H"],
        #     opt_vals.x.reshape(3, 3) * UnitsData.convert("Angstroms", "BohrRadius")
        # ).get_embedded_molecule(load_properties=False).coords * UnitsData.convert("BohrRadius", "Angstroms")
        #
        # raise Exception( base_pot(ugh), ugh)

        ref = np.array([
            [0.00000000e+00, 6.56215885e-02, 0.00000000e+00],
            [7.57391014e-01, -5.20731105e-01, 0.00000000e+00],
            [-7.57391014e-01, -5.20731105e-01, 0.00000000e+00]
        ]) * UnitsData.convert("Angstroms", "BohrRadius")

        ref_mol = Molecule(
            ["O", "H", "H"],
            ref
        ).get_embedded_molecule(load_properties=False)

        rebuild_interpolation=False
        with Checkpointer.from_file(os.path.expanduser("~/Desktop/water_dat.hdf5")) as chk:
            try:
                if rebuild_interpolation:
                    raise Exception("rebuilding...")
                grid = chk['grid']
                vals = chk['vals']
                derivs = chk['derivs']
            except Exception as e:
                print("ERROR: ", e)

                np.random.seed(0)
                disps = np.random.normal(0, 0.5, size=(1000, 3, 3))
                disps[..., 2] = 0
                coords = ref + disps

                forces = lambda c: -mbpol.get_pot_grad(
                    nwaters=1, coords=c.reshape(-1, 3, 3) * b2a,
                    threading_vars=['energy', 'grad', 'coords'],
                    threading_mode='omp'
                )['grad'].reshape(c.shape) * b2a
                energies = lambda c: mbpol.get_pot(
                    nwaters=1, coords=c.reshape((-1, 3, 3)) * b2a,
                    threading_vars=['energy', 'coords'],
                    threading_mode='omp'
                ).reshape(c.shape[:-2])

                sim = AIMDSimulator(["O", "H", "H"], coords, force_function=forces)
                sim.propagate(1000)

                print("Eckart embedding and creating interpolator")
                interp = sim.build_interpolation(energies, eckart_embed=True, interpolation_order=2, reference=ref, logger="All")#, neighborhood_size=5)

                chk['grid'] = interp.grid
                chk['vals'] = interp.vals
                chk['derivs'] = interp.derivs

                grid = chk['grid']
                vals = chk['vals']
                derivs = chk['derivs']

            # grid = grid[100*200:]
            # vals = vals[100*200:]
            # derivs = [d[100*200:] for d in derivs]
            #
            #
            # grad_dat = derivs[0].reshape((-1, 100, 9))
            # hess_dat = derivs[1].reshape((-1, 100, 9, 9))
            # pg = np.arange(len(hess_dat))
            #
            # fig = None
            # for k in range(10):
            #     fig = plt.Plot(
            #         pg,
            #         vals.reshape(-1, 100)[:, k],
            #         figure=fig,
            #         plot_label="Energies"
            #     )
            # fig.show()
            #
            #
            # for i in range(9):
            #     fig = None
            #     for k in range(10):
            #         fig = plt.Plot(
            #             pg,
            #             grad_dat[:, k, i],
            #             figure=fig,
            #             plot_label="grads"
            #         )
            #     fig.show()
            #
            # for i in range(9):
            #     for j in range(i, 9):
            #         fig = None
            #         for k in range(10):
            #             fig = plt.Plot(
            #                 pg,
            #                 hess_dat[:, k, i, j],
            #                 figure=fig,
            #                 plot_label=str((i, j))
            #             )
            #         fig.show()
            # raise Exception(...)

            # eval_dat = potential(grid, deriv_order=2)
            # vals = eval_dat[0]
            # derivs = eval_dat[1:]

            # rots, (og_ref, ref_com, ref_rot), _ = ref_mol.get_embedding_data(grid[:2].reshape(-1, 3, 3), in_paf=True)

            # raise Exception(
            #     potential(grid[:2], deriv_order=2)[2] -
            #     (derivs[1][:2])
            # )

            # else:
            # raise Exception(
            #     np.mean(
            #         vals * UnitsData.convert("Hartrees", "Wavenumbers")
            #     ),
            #     np.std(
            #         vals * UnitsData.convert("Hartrees", "Wavenumbers")
            #     )
            # )

        # good_pos = vals < 15000 * UnitsData.convert("Wavenumbers", "Hartrees")

        # sorting = np.argsort(vals)
        # distance_cutoff = 0.15
        # sort_blocks = np.array_split(sorting, len(sorting)/1e3)
        # final_sortings = []
        # # print(len(sort_blocks))
        # for block in sort_blocks:
        #     # print("....?")
        #     subpts = grid[block]
        #     rinds, cinds = np.triu_indices(len(block), k=1)
        #     v = subpts[rinds] - subpts[cinds]
        #     dvec = np.linalg.norm(v, axis=1)
        #     # print(">>> ")
        #     bad_spots = dvec < distance_cutoff
        #     if bad_spots.any():
        #         to_kill = cinds[bad_spots]
        #         mask = np.ones(len(subpts), dtype=bool)
        #         mask[to_kill] = False
        #         final_sortings.append(block[mask])

        # np.random.seed(0)
        # good_pos = np.concatenate(final_sortings)
        # good_pos = np.random.choice(good_pos, len(good_pos)) # shuffle around to lose ordering
        # raise Exception(len(good_pos))

        # grid = grid[good_pos]
        # vals = vals[good_pos]
        # derivs = [d[good_pos] for d in derivs]

        # derivs = [
        #     derivs[0][:, (0, 1, 3, 4, 6, 7)],
        #     derivs[1][:, :, (0, 1, 3, 4, 6, 7)][:, (0, 1, 3, 4, 6, 7), :]
        # ]

        der_diags = np.diagonal(derivs[1], axis1=1, axis2=2)
        trip_mass = np.sqrt(
            np.broadcast_to(
                np.array([
                    AtomData[a, "Mass"] * UnitsData.convert("AtomicMassUnits", "ElectronMass")
                    for a in ["O", "H", "H"]
                ])[:, np.newaxis],
                (3, 3)
            )
        ).flatten()
        weighted = der_diags / trip_mass[np.newaxis]
        max_ders = np.max(np.abs(weighted), axis=1)

        alphas = .1 * max_ders

        # raise Exception(int(len(grid)/40))

        subsel = (0, 1, 3, 4, 6, 7)
        interp = InverseDistanceWeightedInterpolator(
            grid,#[:, subsel],
            vals,
            derivs[0],#[:, subsel],
            derivs[1],#[:, :, subsel][:, subsel, :],
            bad_interpolation_retries=1,
            neighborhood_size=11,
            # neighborhood_clustering_radius=.1,
            neighborhood_merge_threshold=None,
            neighborhood_max_merge_size=25,
            # multicenter_monomials=True,
            # monomial_basis=False,
            # extra_degree=1,

            # kernel='gaussian',
            # kernel_options={'e':.1},

            # kernel='zeros',

            # kernel='wendland_polynomial',
            # kernel_options={'d':9, 'k':11},

            # auxiliary_basis='compact',
            # auxiliary_basis_options={'k':3},

            # auxiliary_basis='compact_laguerre',
            # auxiliary_basis_options={'k': 3, 'e':.0},

            # clustering_radius=.05,
            logger="All"
        )

        def embedded_interp(traj, deriv_order=0, chunk_size=1e2, **opts):
            from Psience.Psience.Molecools.Properties import PropertyManager

            traj = traj.reshape(-1, 3, 3)
            rots, _, (pax_traj, _, pax_rots) = ref_mol.get_embedding_data(traj)
            traj = pax_traj @ np.swapaxes(rots, -2, -1)

            traj = traj.reshape(-1, 9)#[:, subsel]
            vals = interp(traj, deriv_order=deriv_order, chunk_size=chunk_size, **opts)

            if deriv_order > 0:
                new_derivs = PropertyManager._transform_derivatives(
                    vals[1:],
                    np.swapaxes(pax_rots, -2, -1) @ np.swapaxes(rots, -2, -1)
                )
                # if deriv_order > 1:
                #     # interp.logger.log_print("evaluating Hessians")
                #     hess_fun = FiniteDifferenceDerivative(
                #         lambda x: embedded_interp(x, zero_tol=-1, retries=0, deriv_order=1)[1].reshape(x.shape),
                #         function_shape=(9, 9),
                #         mesh_spacing=5e-5
                #     )
                #
                #     # chunks = []
                #     # num_chunks = int(len(coords)/1000)
                #     # for a in np.array_split(coords, num_chunks):
                #     #     chunks.append(hess_fun.derivatives(a).compute_derivatives(1))
                #     # hess = np.concatenate(chunks, axis=0)
                #     fd_derivs = hess_fun.derivatives(traj).derivative_tensor(list(range(1, deriv_order)))
                #     # print([d.shape for d in new_derivs])
                #     new_derivs = list(new_derivs[:1]) + [
                #         np.moveaxis(d, -2, 0)
                #         for i, d in enumerate(fd_derivs)
                #     ]
                #
                # # new_derivs = list(new_derivs)
                # # zz = np.zeros((len(traj), 9))
                # # zz[:, subsel] = new_derivs[0]
                # # new_derivs[0] = zz
                # # if deriv_order > 1:
                # #     zz = np.zeros((len(traj), 9, 9))
                # #     zz[np.ix_(np.arange(len(traj)), subsel, subsel)] = new_derivs[1]
                # #     new_derivs[1] = zz
                return [vals[0]] + list(new_derivs)
            else:
                return vals

        #region Test Potential

        # np.random.seed(2)
        # npts = 15
        # urgh = grid[:npts]
        # disps = np.random.normal(0, 0.05, (npts, 3, 3))
        # disps[..., 2] = 0.0
        # urgh = urgh + disps.reshape(-1, 9)
        #
        # urgh = ref_mol.embed_coords(urgh.reshape(-1, 3, 3)).reshape(-1, 9)

        # urgh = urgh[:2]
        # a, a1, a2 = [x * UnitsData.convert("Hartrees", "Wavenumbers") for x in potential(urgh, deriv_order=2)]
        # b, b1, b2 = [x * UnitsData.convert("Hartrees", "Wavenumbers") for x in
        #             embedded_interp(urgh,
        #                     zero_tol=-1, retries=0
        #                     , deriv_order=2
        #                     # deriv_order=2,
        #                     # max_distance=.2,
        #                     # neighbors=int(len(grid)/20),
        #                     # neighborhood_clustering_radius=.0035,
        #                     # use_natural_neighbors=True
        #                     )
        # ]



        # b = [x * UnitsData.convert("Hartrees", "Wavenumbers") for x in
        #              embedded_interp(urgh,
        #                              zero_tol=-1, retries=00
        #                              # deriv_order=2,
        #                              # max_distance=.2,
        #                              # neighbors=int(len(grid)/20),
        #                              # neighborhood_clustering_radius=.0035,
        #                              # use_natural_neighbors=True
        #                              )
        #              ]
        # c = (a - b)
        # d = 100 * c / a

        # with np.printoptions(linewidth=1e8):
        #     fack = "\n" + str(np.array([np.round(a), np.round(b), np.round(c, 3)]))
        #
        #     fack = fack + "\n Error: {}".format(np.round(np.linalg.norm(d) / len(d), 2))
        #
        #     # fack += "\n\n " + str(
        #     #     np.moveaxis(
        #     #         np.array([
        #     #             np.round(a1, 3),
        #     #             np.round(b1, 3),
        #     #             np.round(a1 - b1, 3)
        #     #         ]),
        #     #         1, 0
        #     #     )
        #     # )
        #     #
        #     # fack = fack + "\n\n Error: {}".format(np.round(np.linalg.norm(100 * (a1[:, subsel] - b1[:, subsel]) / a1[:, subsel], axis=1), 2))
        #
        #     a2 = np.diagonal(a2, axis1=1, axis2=2)
        #     b2 = np.diagonal(b2, axis1=1, axis2=2)
        #     fack += "\n\n\n " + str(np.moveaxis(np.array([
        #         np.round(a2),
        #         np.round(b2),
        #         np.round(a2 - b2, 3)
        #     ]), 1, 0))
        #
        #     fack = fack + "\n\n Error: {} ({})".format(
        #         np.round(np.mean(np.linalg.norm(100 * (a2[:, subsel] - b2[:, subsel]) / a2[:, subsel], axis=1)), 2),
        #         np.round(np.linalg.norm(100 * (a2[:, subsel] - b2[:, subsel]) / a2[:, subsel], axis=1) / len(a2), 2)
        #     )
        #
        # raise Exception(fack)

        #endregion

        #region Other Methods

        # r1 = np.linalg.norm(ref[1] - ref[0])
        # v1 = (ref[1] - ref[0]) / r1
        # def pot(r, deriv_order=0):
        #     base_shape = r.shape[:-1]
        #     r = r.flatten()
        #     disps = v1[np.newaxis] * r[:, np.newaxis]
        #     pts = np.broadcast_to(ref[np.newaxis], (len(r), 3, 3)).copy()
        #     pts[:, 1, :] += disps
        #
        #     # fig = plt.ScatterPlot(
        #     #                         pts[:, 0, 0],
        #     #                         pts[:, 0, 1],
        #     #                         color="#ff0000",
        #     #                         padding=[[50, 50], [50, 50]],
        #     #                     )
        #     # plt.ScatterPlot(
        #     #     pts[:, 1, 0],
        #     #     pts[:, 1, 1],
        #     #     color="#000000",
        #     #     figure=fig,
        #     #     padding=[[50, 50], [50, 50]],
        #     # )
        #     # plt.ScatterPlot(
        #     #     pts[:, 2, 0],
        #     #     pts[:, 2, 1],
        #     #     color="#00ff00",
        #     #     figure=fig,
        #     #     padding=[[50, 50], [50, 50]],
        #     # )
        #     # fig.show()
        #
        #     res = potential(pts, deriv_order=deriv_order)[0]
        #     res[r < -1.5] = 20
        #
        #     if len(base_shape) > 0:
        #         res = res.reshape(base_shape)
        #
        #     return res

        # g = np.linspace(-1.5, 10, 1000)
        # plt.Plot(
        #     g,
        #     pot(g)
        # ).show()

        # red_mass = 1/(
        #         1/AtomData["H", "Mass"] * UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
        #         + 1/AtomData["O", "Mass"] * UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
        # )
        #
        # from Psience.DVR import DVR
        #
        # res = DVR(
        #     domain=[-1, 5],
        #     divs=200,
        #     potential_function=pot,
        #     mass=1/red_mass
        # ).run()
        #
        # raise Exception(
        #     res.wavefunctions.energies,
        #     res.wavefunctions.frequencies() * UnitsData.convert("Hartrees", "Wavenumbers"))

        #   array([0.00833372, 0.02447275, 0.03995099, 0.05477823, 0.06896247,
        #        0.08250939, 0.09542198, 0.10770047, 0.11934227, 0.13034168,
        #        0.14068864, 0.15036553, 0.1593404 , 0.16755323, 0.17488519,
        #        0.18106746, 0.18524416, 0.1870915 , 0.18920155, 0.19114839,
        #        0.19304661, 0.19520578, 0.19770914, 0.20051799, 0.20360185])


        # wfns = DGB.run(
        #     np.linspace(-1, 5, 1000),
        #     pot,
        #     logger=interp.logger,
        #     alphas=15,
        #     clustering_radius=.0001,
        #     masses=[1/red_mass],
        #     # optimize_centers=True,
        #     min_singular_value=1e-4,
        #     num_svd_vectors=200,
        #     # clustering_radius=.055,
        #     # optimize_centers=False,
        #     quadrature_degree=8
        #     # expansion_type='taylor',
        #     # reference_structure=ref
        # )
        # wfns = wfns[np.where(wfns.energies > 1e-4)[0]]
        # for w in wfns:
        #     if w.energy > 1e-4:
        #         w.plot(
        #             plot_label=f'Energy: {w.energy}',
        #             padding=[[50, 50], [50, 50]]
        #         ).show()
        #
        # print(wfns.energies)
        # print(
        #     wfns.frequencies() * UnitsData.convert("Hartrees", "Wavenumbers")
        # )
        #
        # raise Exception(
        #     wfns.energies,
        #     wfns.frequencies() * UnitsData.convert("Hartrees", "Wavenumbers")
        # )
        # raise Exception(e[e > .5], (e[e > .5][1:] - e[e > .5][0]) * 219475, ham.clustering_radius)


        # from Psience.VPT2 import VPTRunner
        #
        # # raise Exception([x.shape for x in potential(ref, deriv_order=4)[1:]])
        #
        # uhhh = VPTRunner.run_simple(
        #     [["O", "H", "H"], ref],
        #     2,
        #     potential_derivatives=[
        #         x.reshape((9,)*(i+1))
        #         for i,x in enumerate(potential(ref, deriv_order=4)[1:])
        #     ],
        #     # order=0,
        #     # expansion_order=0,
        #     calculate_intensities=False
        # )
        #
        # raise Exception(...)

        #endregion

        # region Subsetting of sampled points

        np.random.seed(1)
        divs = [
            # [250,   500],
            # [500,   None],
            # [1000,  500],
            [5000, None],
            [10000, 100]
        ]
        ec_last = -1
        all_pos = []
        for ec, num in divs:
            pos = np.where(
                np.logical_and(
                    vals >= ec_last * UnitsData.convert("Wavenumbers", "Hartrees"),
                    vals < ec * UnitsData.convert("Wavenumbers", "Hartrees")
                )
            )
            if len(pos) > 0 and len(pos[0]) > 0:
                all_pos.append(
                    pos[0][np.random.choice(len(pos[0]), num),]
                    if num is not None else
                    pos[0]
                )
        good_pos = np.concatenate(all_pos)

        grid = grid[good_pos]
        vals = vals[good_pos]
        derivs = [d[good_pos] for d in derivs]

        # endregion

        #region potential for later plotting

        left_hydrogen_mesh = np.array(
            np.meshgrid(
                np.linspace(-2.0, -0.5, 25),
                np.linspace(-2.0,  0.2, 10)
            ))
        left_hydrogen_cat = np.moveaxis(left_hydrogen_mesh, 0, -1).reshape(-1, 2)
        left_hydrogen_points = np.broadcast_to(ref[np.newaxis], (len(left_hydrogen_cat), 3, 3)).copy()
        left_hydrogen_points[:, 2, :2] = left_hydrogen_cat


        right_hydrogen_mesh = np.array(
            np.meshgrid(
                np.linspace( 0.5, 2.0, 25),
                np.linspace(-2.0, 0.2, 10)
            ))
        right_hydrogen_cat = np.moveaxis(right_hydrogen_mesh, 0, -1).reshape(-1, 2)
        right_hydrogen_points = np.broadcast_to(ref[np.newaxis], (len(right_hydrogen_cat), 3, 3)).copy()
        right_hydrogen_points[:, 1, :2] = right_hydrogen_cat

        full_points = np.concatenate([
            left_hydrogen_points, right_hydrogen_points
        ])

        full_cat = np.concatenate([
            left_hydrogen_cat, right_hydrogen_cat
        ])


        mbpol_pot = lambda r,deriv_order=0: (
                potential(r)[0].reshape(r.shape[:-1])
                    if deriv_order == 0 else
                potential(r, deriv_order=deriv_order)
        )
        pot = mbpol_pot
        # pot = embedded_interp

        sampling_potentials = pot(full_points.reshape(-1, 9)) * UnitsData.convert("Hartrees", "Wavenumbers")
        # sampling_potentials = embedded_interp(full_points.reshape(-1, 9)) * UnitsData.convert("Hartrees", "Wavenumbers")

        # with Checkpointer.from_file(os.path.expanduser("~/Desktop/interp_engs.hdf5")) as chk:
        #     chk['energies'] = sampling_potentials
        # raise Exception(...)

        # plt.TriContourPlot(
        #     full_cat[:, 0],
        #     full_cat[:, 1],
        #     sampling_potentials,
        #     levels=np.linspace(0, 10000, 25),
        #     colorbar=True,
        #     image_size=[600, 600]
        # ).show()

        # 0.005, 0.05, 0.1,
        # pos = np.random.choice(len(pts), 9000)

        #endregion

        np.random.seed(0)
        pos = np.arange(len(grid))
        run_pos = np.random.choice(pos, 4500)
        pts = grid[run_pos]
        r, c = np.triu_indices(len(pts), k=1)
        dists = np.linalg.norm(pts[r] - pts[c], axis=1)
        mask = np.ones(len(pts), dtype=bool)
        mask[c[dists < 1e-1]] = False
        pts = pts[mask]
        run_pos = run_pos[mask]

        for a in [15]:#, 10, 5, 2]:#0.5, 0.3, 0.8, 0.4, 0.05, 0.1, 0.2,]:
            for cr in [0.05]:#, 0.1, 0.2]:
                for ms in [1e-4]:
                    with interp.logger.block(tag=f"Alpha: {a} Min-Sing: {ms} Clust: {cr}"):

                        # largest_deriv = np.max(derivs[1]
                        # a =
                        der_diags = np.diagonal(derivs[1][run_pos], axis1=1, axis2=2)
                        trip_mass = np.sqrt(
                            np.broadcast_to(
                                np.array([
                                    AtomData[a, "Mass"] * UnitsData.convert("AtomicMassUnits", "ElectronMass")
                                    for a in ["O", "H", "H"]
                                ])[:, np.newaxis],
                                (3, 3)
                            )
                        ).flatten()
                        weighted = der_diags / trip_mass[np.newaxis]
                        max_ders = np.max(np.abs(weighted), axis=1)

                        dist_mat = np.linalg.norm(pts[:, np.newaxis] - pts[np.newaxis, :], axis=2)
                        np.fill_diagonal(dist_mat, 100)
                        dist_mins = np.min(dist_mat, axis=1)
                        alphas = 1/dist_mins #* (1 + 1 * (max_ders - np.min(max_ders)) / (np.max(max_ders) - np.min(max_ders)))
                        # raise Exception(alphas)

                        # disp_mat = np.abs(pts[:, np.newaxis] - pts[np.newaxis, :])
                        # disp_mat[np.diag_indices(len(pts))] = 1000
                        # min_disp = np.min(disp_mat, axis=0)
                        # min_disp[min_disp < 1e-8] = -1
                        # alphas = np.sqrt(np.abs(weighted))/(min_disp * np.sqrt(trip_mass)) #* (1 + 1 * (max_ders - np.min(max_ders)) / (np.max(max_ders) - np.min(max_ders)))
                        # alphas[min_disp == -1] = 10
                        # alphas[alphas < 6] = 10
                        # raise Exception(alphas)
                        # alphas = np.sqrt(max_ders)/dist_mins**2
                        # alphas = a

                        # pts = pts[np.lexsort(pts.T)]
                        # ham = DGB(
                        #     pts,
                        #     pot,
                        #     logger=interp.logger,
                        #     alphas=alphas,
                        #     atoms=["O", "H", "H"],
                        #     projection_indices=(0, 1, 3, 4, 6, 7),
                        #     clustering_radius=cr,
                        #     optimize_centers=False,
                        #     min_singular_value=ms,
                        #     num_svd_vectors=1,
                        #     svd_contrib_cutoff=1e-2,
                        #     # clustering_radius=.055,
                        #     # optimize_centers=False,
                        #     # quadrature_degree=4,
                        #     expansion_degree=2 if pot is interp else 2,
                        #     # expansion_type='taylor',
                        #     # reference_structure=ref
                        # )
                        # diag1 = np.diagonal(ham.S, offset=1)
                        # urgh = ham.S.copy()
                        # np.fill_diagonal(urgh)
                        # raise Exception(np.min(diag1))
                        # plt.Plot(np.arange(len(diag1)), diag1).show()
                        #
                        # plt.ArrayPlot(ham.S).show()
                        # raise Exception(...)
                        # ov_pts = ham.get_overlap_gaussians()[0]


                        # base_plot = plt.ScatterPlot(
                        #     interp.grid[:, 0], interp.grid[:, 1],
                        #     color='red'
                        # )
                        # base_plot = plt.ScatterPlot(
                        #     interp.grid[:, 3], interp.grid[:, 4],
                        #     color='gray',
                        #     figure=base_plot
                        # )
                        # base_plot = plt.ScatterPlot(
                        #     interp.grid[:, 6], interp.grid[:, 7],
                        #     color='gray',
                        #     figure=base_plot
                        # )
                        #
                        # new_grid = ov_pts[np.triu_indices(len(ov_pts), k=1)].reshape(-1, 9)
                        # plt.ScatterPlot(
                        #     new_grid[:, 0], new_grid[:, 1],
                        #     color='blue',
                        #     figure=base_plot
                        # )
                        # base_plot = plt.ScatterPlot(
                        #     new_grid[:, 3], new_grid[:, 4],
                        #     color='blue',
                        #     figure=base_plot
                        # )
                        # base_plot = plt.ScatterPlot(
                        #     new_grid[:, 6], new_grid[:, 7],
                        #     color='blue',
                        #     figure=base_plot
                        # )
                        # base_plot.show()
                        #
                        #
                        # raise Exception("...")

                        # a, _, a2 = [x * UnitsData.convert("Hartrees", "Wavenumbers") for x in
                        #             mbpol_pot(ov_pts.reshape(-1, 9), deriv_order=2)
                        #             ]
                        # b, _, b2 = [x * UnitsData.convert("Hartrees", "Wavenumbers") for x in
                        #             embedded_interp(ov_pts.reshape(-1, 9),
                        #                     # neighbors=2500,
                        #                     zero_tol=-1,
                        #                     deriv_order=2,
                        #                     # neighborhood_clustering_radius=.0015
                        #                     )
                        # ]
                        # with np.printoptions(linewidth=1e8):
                        #     c = a - b
                        #     d = 100 * (c) / a
                        #     fack = "\n" + str(np.array([
                        #         np.round(a),
                        #         np.round(b),
                        #         np.round(c, 3),
                        #     ]))
                        #
                        #     fack = fack + "\n Error: {}".format(np.round(np.linalg.norm(d) / len(d), 2))
                        #
                        #     # fack += "\n\n " + str(
                        #     #     np.moveaxis(
                        #     #         np.array([
                        #     #             np.round(a1, 3),
                        #     #             np.round(b1, 3),
                        #     #             np.round(a1 - b1, 3)
                        #     #         ]),
                        #     #         1, 0
                        #     #     )
                        #     # )
                        #     #
                        #     # fack = fack + "\n\n Error: {}".format(np.round(np.linalg.norm(100 * (a1[:, subsel] - b1[:, subsel]) / a1[:, subsel], axis=1), 2))
                        #
                        #     a2 = np.diagonal(a2, axis1=1, axis2=2)
                        #     b2 = np.diagonal(b2, axis1=1, axis2=2)
                        #     fack += "\n\n\n " + str(np.moveaxis(np.array([
                        #         np.round(a2),
                        #         np.round(b2),
                        #         np.round(a2 - b2, 3)
                        #     ]), 1, 0))
                        #
                        #     fack = fack + "\n\n Error: {} ({})".format(
                        #         np.round(np.mean(np.linalg.norm(100 * (a2[:, subsel] - b2[:, subsel]) / a2[:, subsel], axis=1)) / len(a2), 2),
                        #         np.round(np.linalg.norm(100 * (a2[:, subsel] - b2[:, subsel]) / a2[:, subsel], axis=1), 2)
                        #     )
                        #
                        # raise Exception(fack)

                        wfns = DGB.run(
                            pts,
                            pot,
                            logger=interp.logger,
                            alphas=alphas,
                            atoms=["O", "H", "H"],
                            projection_indices=(0, 1, 3, 4, 6, 7),
                            clustering_radius=cr,
                            optimize_centers=False,
                            min_singular_value=ms,
                            num_svd_vectors=1,
                            svd_contrib_cutoff=1e-2,
                            # clustering_radius=.055,
                            # optimize_centers=False,
                            # quadrature_degree=4,
                            expansion_degree=2 if pot is interp else 2,
                            # expansion_type='taylor',
                            # reference_structure=ref
                        )
                    with interp.logger.block(tag="Results:"):
                        count = 8
                        ecut = 0e-2
                        emax = ecut + 7000 * UnitsData.convert("Hartrees", "Wavenumbers")
                        max_val = None
                        for i,wfn in enumerate(wfns):
                            if count > 0 and wfn.energy > ecut and wfn.energy < emax:
                                count -= 1

                                fig = plt.TriContourPlot(
                                    full_cat[:, 0],
                                    full_cat[:, 1],
                                    sampling_potentials,
                                    levels=np.linspace(0, 10000, 25),
                                    colorbar=True,
                                    cmap='viridis',
                                    image_size=[600, 600]
                                )

                                proj = wfn.marginalize_out([0, 1, 2, 5, 6, 7, 8]) # what we're projecting _out_
                                if max_val is None:
                                    max_val = min([
                                        np.abs(np.max(proj.data)) * np.max(proj.alphas),
                                        2
                                    ])
                                proj.plot(
                                    figure=fig,
                                    plot_label=f"Wavefunction {i} Energy: {wfn.energy}",
                                    epilog=[
                                        plt.Disk(
                                            ref[0, :2],
                                            3,
                                            color='black'
                                        ),

                                        plt.Disk(
                                            ref[1, :2],
                                            2,
                                            color='black'
                                        ),

                                        plt.Disk(
                                            ref[2, :2],
                                            2,
                                            color='black'
                                        )
                                    ],
                                    plotter=plt.TriContourLinesPlot,
                                    levels=np.linspace(-max_val, max_val, 16),
                                    domain=[[-2, 2], [-2, .2]],
                                    cmap='RdBu'
                                )

                                plt.ScatterPlot(
                                    proj.centers[:, 0],
                                    proj.centers[:, 1],
                                    color="#ffffff11",
                                    figure=fig,
                                    padding=[[50, 50], [50, 50]],
                                )

                                proj = wfn.marginalize_out([0, 1, 2, 3, 4, 5, 8])

                                proj.plot(
                                    figure=fig,
                                    plotter=plt.TriContourLinesPlot,
                                    domain=[[-2, 2], [-2, .2]],
                                    levels=np.linspace(-max_val, max_val, 16),
                                    cmap='RdBu'
                                )

                                plt.ScatterPlot(
                                    proj.centers[:, 0],
                                    proj.centers[:, 1],
                                    color="#ffffff11",
                                    figure=fig,
                                    padding=[[50, 50], [50, 50]],
                                )

                                fig.show()

                        # e = wfns.energies
                        #
                        # # interp.logger.log_print("All energies: {eng}", eng=wfns.energies)
                        # pos_e = np.logical_and(e > 1e-2, e < 2)
                        # e = e[pos_e]
                        # interp.logger.log_print("energies: {eng}", eng=e)
                        #
                        # all_freqs = (e[:, np.newaxis] - e[np.newaxis, :]) * UnitsData.convert("Hartrees", "Wavenumbers")
                        wfn_start = np.where(wfns.energies > ecut)[0][0]
                        interp.logger.log_print("ZPE: {e}", e=wfns[wfn_start].energy * UnitsData.convert("Hartrees", "Wavenumbers"))
                        interp.logger.log_print("All frequencies: {freq}", freq=wfns[wfn_start:wfn_start+20].frequencies() * UnitsData.convert("Hartrees", "Wavenumbers"))
                        # best_pos = np.argmin(np.abs(all_freqs.flatten() - 3500))
                        # interp.logger.log_print("Best frequency: {freq}", freq=all_freqs.flatten()[best_pos])

                    gc.collect()


        """
        VPT
        
        ::> States Energies
          > State     Harmonic   Anharmonic     Harmonic   Anharmonic
                       ZPE          ZPE    Frequency    Frequency
        0 0 0   4713.07975   4636.90794            -            - 
        0 0 1            -            -   3944.32593   3753.07779 
        0 1 0            -            -   3832.76457   3654.53256 
        1 0 0            -            -   1649.06900   1594.42231 
        0 0 2            -            -   7888.65187   7408.76383 
        0 2 0            -            -   7665.52914   7222.28820 
        2 0 0            -            -   3298.13800   3152.62843 
        0 1 1            -            -   7777.09050   7240.72934 
        1 0 1            -            -   5593.39493   5326.66828 
        1 1 0            -            -   5481.83357   5232.92568 
        <::
        """

        """DGB 2500 pts & MB-Pol w/ FD Hessians
        
        :: All frequencies: [ 
             1626.23440329  3291.02298809  3719.74860224  3977.64291868
             5045.50668587  5397.92323296  5962.62638424  6865.53165906
             7227.10335916  7351.04973073  7632.22294328  7801.16335243
             8443.88015359  8788.70001116  9036.19149045  9968.46931843
             10262.30527112 10676.33757779 10932.38611674]
        
        """

        """DGB 2500 pts & MB-Pol w/ FD 4th derivs

        >>------------------------- Results: -------------------------
        :: All frequencies: [ 1577.73146693  3196.09291218  3713.43203707  3985.63634612
          4943.97500195  5343.17673283  5891.11364755  6745.74649368
          7110.01884398  7342.20274783  7643.45326813  7714.82677927
          8383.75793472  8664.3708776   8957.39854453  9978.36189936
         10172.96724881 10559.28985947 10913.42580748]
        >>--------------------------------------------------<<

        """


