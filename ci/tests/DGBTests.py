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

        ham = DGB(pts, pot, clustering_radius=.05)
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
    def test_Interp(self):
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
        ndim = 1
        pts = np.random.uniform(low=-.5, high=1.2, size=(1000, ndim))

        # fn = sym.morse(sym.x) * sym.morse(sym.y) - sym.morse(sym.x) - sym.morse(sym.y)

        w = 2.05; wx = .1; mu = 1
        de = (w ** 2) / (4 * wx)
        a = np.sqrt(2 * mu * wx)
        fn = sym.morse(sym.x, de=de, a=a) + sym.morse(sym.y, de=de, a=a)
        fn = (1/2*sym.x**2)# + 1/2*sym.y**2)
        # raise Exception(fn(pts))

        # print(res['AIMDEnergies'].gradients)
        # interp = RBFDInterpolator.create_function_interpolation(
        #     pts,
        #     fn,
        #     lambda p, f=fn.deriv(order=1): f(p),#.transpose(),
        #     lambda p, f=fn.deriv(order=2): f(p),#.transpose(2, 0, 1),
        #     clustering_radius=1e-5,
        # )

        test_pot = lambda c,fn=fn,deriv_order=None: (
            fn(c).reshape(c.shape[:-1])
                if deriv_order is None else
            fn.deriv(order=deriv_order)(c).reshape(c.shape[:-1] + (ndim,)*deriv_order)
        )


        np.random.seed(0)
        centers = pts[np.unique(np.random.random_integers(low=0, high=len(pts)-1, size=100))]
        centers = np.linspace(-1, 1, 50).reshape(-1, 1)
        ham = DGB(
            centers,
            test_pot,
            alphas=1,
            clustering_radius=0.1
        )
        e, wf = ham.get_wavefunctions()
        print(e)

        # for k in range(1, 6):
        #     print(DGB.polyint_1D_coeffs(k))
        # raise Exception(DGB.polyint_1D_coeffs(4))

        ham = DGB(
            centers,
            test_pot,
            expansion_degree=6,
            alphas=1
        )
        # ham = DGB(
        #     centers,
        #     test_pot,
        #     expansion_degree=2, expansion_type='taylor', alphas=1, clustering_radius=.1)
        e, wf = ham.get_wavefunctions()
        print(e)

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


