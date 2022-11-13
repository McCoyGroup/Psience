import itertools

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

    @debugTest
    def test_Expansion(self):
        from McUtils.Zachary import Symbols, RBFDInterpolator

        sym = Symbols('xyz')

        # np.random.seed(3)
        # fn = sym.morse(sym.x)
        # pts = np.random.uniform(low=-.5, high=1.2, size=30)

        # test_1D(sym.morse(sym.x), pts)
        # test_1D(sym.sin(sym.x), pts)
        # test_1D(sym.morse(sym.x)*sym.sin(sym.x), pts)
        #
        # raise Exception(...)

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
        # test_pot = simple_morse

        # with np.printoptions(linewidth=1e8):
        #     print("...")
        #     print(simple_morse(pts[:3], deriv_order=3))
        #     print("...")
        #     print(test_pot(pts[:3], deriv_order=3))
        #     print("...")
        #     print(simple_morse(pts[:3], deriv_order=3) - test_pot(pts[:3], deriv_order=3))
        # raise Exception(...)

        # interp = RBFDInterpolator.create_function_interpolation(
        #     pts,
        #     fn,
        #     lambda p, f=fn.deriv(order=1): f(p),#.transpose(),
        #     lambda p, f=fn.deriv(order=2): f(p),#.transpose(2, 0, 1),
        #     clustering_radius=1e-5,
        # )

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
            centers,
            test_pot,
            alphas=alpha,
            # clustering_radius=cluster,
            expansion_degree=6
        )
        e, wf = ham.get_wavefunctions()#print_debug_info=True)
        print(e[:10])

        ham = DGB(
            centers,
            test_pot,
            alphas=alpha,
            # clustering_radius=cluster,
            expansion_degree=6,
            expansion_type='taylor'
        )
        e, wf = ham.get_wavefunctions()#print_debug_info=True)
        print(e[:10])

    @inactiveTest
    def test_Interp(self):
        from McUtils.Zachary import Symbols, RBFDInterpolator

        sym = Symbols('xyz')

        np.random.seed(3)
        ndim = 2
        pts = np.random.uniform(low=-5, high=5, size=(1000, ndim))

        w = 2.05; wx = .1; mu = 1
        de = (w ** 2) / (4 * wx)
        a = np.sqrt(2 * mu * wx)
        fn = sym.morse(sym.x, de=de, a=a) + (0 if ndim == 1 else sym.morse(sym.y, de=de, a=a))

        # fn = 1/2*sym.x**2 + (0 if ndim == 1 else 1/2*sym.y**2)
        test_pot = lambda c, fn=fn, deriv_order=None: (
            fn(c.reshape(-1, ndim)).reshape(c.shape[:-1])
            if deriv_order is None else
            np.moveaxis(fn.deriv(order=deriv_order)(c.reshape(-1, ndim)), -1, 0).reshape(
                c.shape[:-1] + (ndim,) * deriv_order)
        )

        interp = RBFDInterpolator.create_function_interpolation(
            pts,
            test_pot,
            lambda p, f=fn.deriv(order=1): f(p),#.transpose(),
            lambda p, f=fn.deriv(order=2): f(p),#.transpose(2, 0, 1),
            clustering_radius=1e-5,
        )

        # raise Exception(fn(pts))

        # print(res['AIMDEnergies'].gradients)


        # test_pot = harmonic

        np.random.seed(0)
        centers = pts  # [np.unique(np.random.random_integers(low=0, high=len(pts)-1, size=100))]
        # centers = np.array(
        #     np.meshgrid(*[
        #         np.linspace(-5, 5, 50)
        #     ]*ndim)
        # ).T.reshape((-1, ndim))

        alpha = .3
        cluster = .05

        ham = DGB(
            centers,
            interp,
            alphas=alpha,
            clustering_radius=cluster,
            quadrature_degree=4
        )
        e, wf = ham.get_wavefunctions()
        print(e[:10])

        ham = DGB(
            centers,
            test_pot,
            expansion_degree=4,
            alphas=alpha,
            clustering_radius=cluster
        )
        e, wf = ham.get_wavefunctions()  # print_debug_info=True)
        print(e[:10])

        # ham = DGB(
        #     centers,
        #     test_pot,
        #     expansion_degree=8,
        #     expansion_type='taylor',
        #     alphas=alpha,
        #     clustering_radius=cluster
        # )
        # e, wf = ham.get_wavefunctions()  # print_debug_info=True)
        # print(e[:10])

    @inactiveTest
    def test_Water(self):
        with GaussianLogReader(TestManager.test_data('h2o_aimd.log')) as parser:
            res = parser.parse(['AIMDCoordinates', 'AIMDEnergies'])

        crds = res['AIMDCoordinates'].reshape(-1, 9)
        e = (res['AIMDEnergies'].energies - np.min(res['AIMDEnergies'].energies)) #* 219474.65
        # print(res['AIMDEnergies'].gradients)
        interp = McUtils.Zachary.RBFDInterpolator(
            crds,
            e,
            res['AIMDEnergies'].gradients,
            res['AIMDEnergies'].hessians,
            clustering_radius=0
        )

        ham = DGB(crds, interp, alphas=1, clustering_radius=.1)
        e, wf = ham.get_wavefunctions()

        print(e)


