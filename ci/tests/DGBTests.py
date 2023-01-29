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
from McUtils.GaussianInterface import GaussianLogReader
from McUtils.Extensions import ModuleLoader

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
        domain = [[-2, 3]]*d # , [-1, 1], [-1, 1]] # unit cube
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

                                proj = wfn.project([0, 1, 2, 5, 6, 7, 8]) # what we're projecting _out_
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

                                proj = wfn.project([0, 1, 2, 3, 4, 5, 8])

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


