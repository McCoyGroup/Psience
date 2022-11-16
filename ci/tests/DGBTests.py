import itertools

import scipy.linalg

import McUtils.Zachary
from Peeves.TestUtils import *
from unittest import TestCase

from McUtils.Data import UnitsData, PotentialData
from McUtils.Zachary import Interpolator
import McUtils.Plots as plt
from McUtils.GaussianInterface import GaussianLogReader

from Psience.DGB import *
from Psience.Molecools import Molecule
import numpy as np

class DGBTests(TestCase):

    @validationTest
    def test_Harmonic(self):
        ndivs = [8, 8, 8]
        domain = [[-1, 1], [-1, 1], [-1, 1]] #, [-1, 1], [-1, 1]] # unit cube
        ndim = len(domain)

        pts = np.array(
            np.meshgrid(*(np.linspace(d[0], d[1], n) for d,n in zip(domain, ndivs)))
        ).T.reshape(-1, ndim)

        np.random.seed(0)
        npts = int(np.prod(ndivs))
        pts = np.random.uniform(low=[d[0] for d in domain], high=[d[1] for d in domain], size=(npts, ndim) )

        pot = lambda c:np.sum(c**2, axis=-1)/2

        ham = DGB(pts, pot, alphas=1, clustering_radius=.05)
        e, wf = ham.get_wavefunctions()

        test_es = np.sort(np.sum(list(itertools.product(*[np.arange(8)+1/2]*ndim)), axis=-1))
        self.assertLess(
            np.linalg.norm(
                e[:8] - test_es[:8]
            ),
            .0025
        )

    @validationTest
    def test_Morse(self):
        d = 2
        ndivs = [100]*d
        domain = [[-2, 3]]*d # , [-1, 1], [-1, 1]] # unit cube
        ndim = len(domain)

        pts = np.array(
            np.meshgrid(*(np.linspace(d[0], d[1], n) for d, n in zip(domain, ndivs)))
        ).T.reshape(-1, ndim)

        w = 2.05; wx=.1; mu=1
        de = (w ** 2) / (4 * wx)
        a = np.sqrt(2 * mu * wx)

        pot = lambda c,de=de,a=a,mu=mu: np.sum(de*(1-np.exp(-a*c))**2, axis=-1)

        ham = DGB(pts, pot, alphas=1, clustering_radius=.1)
        e, wf = ham.get_wavefunctions()

        base_energies = [w*(np.arange(3) + 1 / 2) - wx*(np.arange(3) + 1 / 2)**2] * ndim

        test_es = np.sort(np.sum(list(itertools.product(*base_energies)), axis=-1))

        self.assertLess(
            np.linalg.norm(
                e[:3] - test_es[:3]
            ),
            .0025
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

        ham = DGB(
            centers,
            test_pot,
            alphas=alpha,
            # clustering_radius=cluster,
            quadrature_degree=6
        )
        e, wf = ham.get_wavefunctions()#print_debug_info=True)
        print(e[:10])

        ham = DGB(
            ham.centers, # reuse the opt from before
            test_pot,
            optimize_centers=False,
            alphas=alpha,
            # clustering_radius=cluster,
            expansion_degree=6
        )
        e, wf = ham.get_wavefunctions()#print_debug_info=True)
        print(e[:10])

        ham = DGB(
            ham.centers, # reuse the opt from before
            test_pot,
            optimize_centers=False,
            alphas=alpha,
            # clustering_radius=cluster,
            expansion_degree=6,
            expansion_type='taylor'
        )
        e, wf = ham.get_wavefunctions()#print_debug_info=True)
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

    @debugTest
    def test_Water(self):
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

