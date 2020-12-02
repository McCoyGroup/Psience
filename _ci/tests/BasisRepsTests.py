from Peeves.TestUtils import *
from unittest import TestCase
from Psience.BasisReps import *
import sys, os, numpy as np

class BasisSetTests(TestCase):


    def get_states(self, n_quanta, n_modes, max_quanta = None):
        import itertools as ip

        if max_quanta is None:
            max_quanta = n_quanta
        return tuple(sorted(
            [p for p in ip.product(*(range(n_quanta+1) for i in range(n_modes))) if all(x<=max_quanta for x in p) and sum(p) <= n_quanta],
            key=lambda p: (
                    sum(p)
                    + sum(1 for v in p if v != 0) * n_quanta ** (-1)
                    + sum(v * n_quanta ** (-i - 2) for i, v in enumerate(p))
            )
        ))

    #region 1D Basis Tests
    @debugTest
    def test_HOBasis1DX(self):
        n = 7
        basis = HarmonicOscillatorBasis(n)

        term = ['x']
        iphase = (-1) ** (term.count("p") // 2)
        rep1 = basis.representation(*term)
        rep2 = Representation(super(HarmonicOscillatorBasis, basis).operator(*term), basis)
        xx = rep1[:, :].todense()
        x2 = iphase * rep2[:, :].todense()

        self.assertLess(np.average(np.abs(xx - x2)), 1e-14)

    @debugTest
    def test_HOBasis1DXX(self):

        n = 7
        basis = HarmonicOscillatorBasis(n)

        rep1 = basis.representation('x', 'x')
        rep2 = Representation(super(HarmonicOscillatorBasis, basis).operator('x', 'x'), basis)
        xx = rep1[:, :].todense()
        x2 = rep2[:, :].todense()
        targ =np.zeros((n, n))
        targ[np.arange(n), np.arange(n)] = np.arange(n) + 1/2
        targ[np.arange(n-2), np.arange(2, n)] = np.sqrt(np.arange(1, n-1)*(np.arange(1, n-1)+1)/4)
        targ[np.arange(2, n), np.arange(n-2)] = np.sqrt(np.arange(1, n-1)*(np.arange(1, n-1)+1)/4)

        # raise Exception([
        #     targ**2,
        #     xx**2
        #     ])

        self.assertLess(np.average(np.abs(x2 - targ)), 1e-14)
        self.assertLess(np.average(np.abs(xx - targ)), 1e-14)

    @debugTest
    def test_HOBasis1DPXP(self):
        n = 7
        basis = HarmonicOscillatorBasis(n)

        term = ['p', 'x', 'p']
        iphase = (-1) ** (term.count("p") // 2)
        rep1 = basis.representation(*term)
        rep2 = Representation(super(HarmonicOscillatorBasis, basis).operator(*term), basis)
        xx = rep1[:, :].todense()
        x2 = iphase * rep2[:, :].todense()

        self.assertLess(np.average(np.abs(xx - x2)), 1e-14)

    @debugTest
    def test_HOBasis1DPP(self):
        n = 7
        basis = HarmonicOscillatorBasis(n)

        rep1 = basis.representation('p', 'p')
        rep2 = Representation(super(HarmonicOscillatorBasis, basis).operator('p', 'p'), basis)
        xx = rep1[:, :].todense()
        x2 = -rep2[:, :].todense()
        # targ = np.zeros((n, n))
        # targ[np.arange(n), np.arange(n)] = np.arange(n) + 1 / 2
        # targ[np.arange(n - 2), np.arange(2, n)] = np.sqrt(np.arange(1, n - 1) * (np.arange(1, n - 1) + 1) / 4)
        # targ[np.arange(2, n), np.arange(n - 2)] = np.sqrt(np.arange(1, n - 1) * (np.arange(1, n - 1) + 1) / 4)

        # raise Exception([
        #     targ**2,
        #     xx**2
        #     ])

        self.assertLess(np.average(np.abs(xx - x2)), 1e-14)

    @debugTest
    def test_HOBasis1DXXX(self):
        n = 7
        basis = HarmonicOscillatorBasis(n)

        term = ['x', 'x', 'x']
        iphase = (-1) ** (term.count("p") // 2)
        rep1 = basis.representation(*term)
        rep2 = Representation(super(HarmonicOscillatorBasis, basis).operator(*term), basis)
        xx = rep1[:, :].todense()
        x2 = iphase * rep2[:, :].todense()

        self.assertLess(np.average(np.abs(xx - x2)), 1e-14)

    @debugTest
    def test_HOBasis1DPPXX(self):
        n = 7
        basis = HarmonicOscillatorBasis(n)

        term = ['p', 'p', 'x', 'x']
        iphase = (-1) ** (term.count("p") // 2)
        rep1 = basis.representation(*term)
        rep2 = Representation(super(HarmonicOscillatorBasis, basis).operator(*term), basis)
        xx = rep1[:, :].todense()
        x2 = iphase * rep2[:, :].todense()
        # targ = np.zeros((n, n))
        # targ[np.arange(n), np.arange(n)] = np.arange(n) + 1 / 2
        # targ[np.arange(n - 2), np.arange(2, n)] = np.sqrt(np.arange(1, n - 1) * (np.arange(1, n - 1) + 1) / 4)
        # targ[np.arange(2, n), np.arange(n - 2)] = np.sqrt(np.arange(1, n - 1) * (np.arange(1, n - 1) + 1) / 4)

        # raise Exception([
        #     targ**2,
        #     xx**2
        #     ])

        self.assertLess(np.average(np.abs(xx - x2)), 1e-14)

    #endregion

    @validationTest
    def test_HOBasis2DXX(self):
        from Peeves import Timer, BlockProfiler

        n = 7
        m = 10
        oppo = HarmonicOscillatorProductBasis((n,) * m)
        oppo2 = SimpleProductBasis(HarmonicOscillatorBasis, (n,) * m)

        term = ['x', 'x']
        iphase = (-1) ** (term.count("p") // 2)
        n_terms = len(term)
        xxpp1 = oppo.representation(*term)
        xxpp2 = oppo2.representation(*term)

        usr = os.path.expanduser('~')
        job_is_dumb = [
            os.path.join(usr, "Documents/Python/config/python3.7/lib/python3.7/"),
            os.path.join(usr, "Documents/UW/Research/Development")
        ]

        states = (
            (0, 0, 0, 0, 0),
            (0, 1, 2, 3, 4)
        )
        # with Timer("New style"):
            # with BlockProfiler("New Style", strip_dirs=job_is_dumb, num_lines=10, sort_by='tottime', filter="Psience"):
        vals1 = xxpp1[states]
        self.assertEquals(vals1.shape, (m,) * n_terms + (len(states[0]),))

        # with Timer("Old style"):
            # with BlockProfiler("Old Style", strip_dirs=job_is_dumb, num_lines=10, sort_by='tottime', filter="Psience"):
        vals2 = xxpp2[states]
        self.assertEquals(vals2.shape, (m,) * n_terms + (len(states[0]),))

        wat = np.roll(np.arange(n_terms + 1), 1)
        # print(wat)
        v1 = vals1.toarray().transpose(wat)
        v2 = iphase * vals2.toarray().transpose(wat)

        # print([np.max(v) for v in v1])
        # print([np.max(v) for v in v2])
        # print([np.max(v) for v in np.abs(v1 - v2)])
        # print(v1[0], v2[0])
        # print(v1[1], v2[1])
        # print(vals1.toarray()[:, :, -1] - vals2.toarray()))

        self.assertLess(np.max(np.abs(v1 - v2)), 1.0e-14)

    @validationTest
    def test_HOBasis2DPP(self):
        from Peeves import Timer, BlockProfiler

        n = 7
        m = 10
        oppo = HarmonicOscillatorProductBasis((n,) * m)
        oppo2 = SimpleProductBasis(HarmonicOscillatorBasis, (n,) * m)

        term = ['p', 'p']
        iphase = (-1) ** (term.count("p") // 2)
        n_terms = len(term)
        xxpp1 = oppo.representation(*term)
        xxpp2 = oppo2.representation(*term)

        usr = os.path.expanduser('~')
        job_is_dumb = [
            os.path.join(usr, "Documents/Python/config/python3.7/lib/python3.7/"),
            os.path.join(usr, "Documents/UW/Research/Development")
        ]

        states = (
            (0, 0, 0, 0, 0),
            (0, 1, 2, 3, 4)
        )
        # with Timer("New style"):
            # with BlockProfiler("New Style", strip_dirs=job_is_dumb, num_lines=10, sort_by='tottime', filter="Psience"):
        vals1 = xxpp1[states]
        self.assertEquals(vals1.shape, (m,) * n_terms + (len(states[0]),))

        # with Timer("Old style"):
            # with BlockProfiler("Old Style", strip_dirs=job_is_dumb, num_lines=10, sort_by='tottime', filter="Psience"):
        vals2 = xxpp2[states]
        self.assertEquals(vals2.shape, (m,) * n_terms + (len(states[0]),))

        wat = np.roll(np.arange(n_terms + 1), 1)
        # print(wat)
        v1 = vals1.toarray().transpose(wat)
        v2 = iphase * vals2.toarray().transpose(wat)

        # print([np.max(v) for v in v1])
        # print([np.max(v) for v in v2])
        # print([np.max(v) for v in np.abs(v1 - v2)])
        # print(v1[0], v2[0])
        # print(v1[1], v2[1])
        # print(vals1.toarray()[:, :, -1] - vals2.toarray()))

        self.assertLess(np.max(np.abs(v1 - v2)), 1.0e-14)

    @validationTest
    def test_HOBasis3DXXX(self):
        from Peeves import Timer, BlockProfiler

        n = 15
        m = 5
        oppo = HarmonicOscillatorProductBasis((n,)*m)
        oppo2 = SimpleProductBasis(HarmonicOscillatorBasis, (n,)*m)

        term = ['x', 'x', 'x']
        iphase = (-1)**(term.count("p") // 2)
        n_terms = len(term)
        xxpp1 = oppo.representation(*term)
        xxpp2 = oppo2.representation(*term)

        usr = os.path.expanduser('~')
        job_is_dumb = [
            os.path.join(usr, "Documents/Python/config/python3.7/lib/python3.7/"),
            os.path.join(usr, "Documents/UW/Research/Development")
        ]

        quant_states = self.get_states(1, m)
        states = oppo.ravel_state_inds(quant_states)
        with Timer("New style"):
            # with BlockProfiler("New Style", strip_dirs=job_is_dumb, num_lines=10, sort_by='tottime', filter="Psience"):
                vals1 = xxpp1[np.ix_(states, states)]
        self.assertEquals(vals1.shape, (m,)*n_terms + (len(states), len(states)))

        with Timer("Old style"):
            # with BlockProfiler("Old Style", strip_dirs=job_is_dumb, num_lines=10, sort_by='tottime', filter="Psience"):
                vals2 = xxpp2[np.ix_(states, states)]
        self.assertEquals(vals2.shape, (m,)*n_terms + (len(states), len(states)))

        # print(wat)
        v1 = vals1.toarray()
        v2 = iphase * vals2.toarray()

        # print([np.max(v) for v in v1])
        # print([np.max(v) for v in v2])
        # print([np.max(v) for v in np.abs(v1 - v2)])

        # for i in range(m):
        #     for j in range(m):
        #         for k in range(m):
        #             r = np.round(v1[i, j, k] - v2[j, i, k], 3)
        #             if np.max(np.abs(r)) > 0:
        #                 print((i, j, k), v1[i, j, k])
        #                 print((i, j, k), v2[i, j, k])
                        # print((i, j, k), r)
            # print(v2[i, i, i])
        # print([np.where(v != 0.) for v in np.round(np.abs(v1 - v2), 5)])
        # print(v1[0], v2[0])
        # print(v1[1], v2[1])
        # print(vals1.toarray()[:, :, -1] - vals2.toarray()))

        self.assertLess(np.max(np.abs(v1 - v2)), 1.0e-14)

        # eval = oppo[0, 1, 0, 1] # same as [1, 0, 1, 0] but not [0, 0, 1, 1]
        # self.assertEquals(eval.terms, (("x", "p"), ("x", "p")))

    @debugTest
    def test_HOBasis3DPXP(self):
        from Peeves import Timer, BlockProfiler

        n = 10
        m = 3
        oppo = HarmonicOscillatorProductBasis((n,) * m)
        oppo2 = SimpleProductBasis(HarmonicOscillatorBasis, (n,) * m)

        term = ['p', 'x', 'p']
        iphase = (-1) ** (term.count("p") // 2)
        n_terms = len(term)

        g1 = np.array([
            [[-1.81146079e-04,  3.97836803e-05,  2.91649691e-05],
             [ 3.97836803e-05,  2.63572358e-05,  2.37597837e-04],
             [ 2.91649691e-05,  2.37597837e-04, -3.38457268e-05]],

            [[-4.36589189e-04,  2.79004059e-05, -1.50059967e-05],
             [ 2.79004059e-05, -1.44188965e-06,  3.49657651e-06],
             [-1.50059967e-05,  3.49657652e-06,  3.11501367e-06]],

            [[-8.10821036e-04,  6.31615150e-06,  5.50255712e-05],
             [ 6.31615151e-06,  4.05569426e-06,  3.51303496e-08],
             [ 5.50255712e-05,  3.51303696e-08, -3.80070492e-06]]])
        xxpp1 = 2*oppo.representation(*term, coeffs=g1, axes=[[0, 1, 2], [1, 0, 2]])
        xxpp1 = xxpp1 + xxpp1
        xxpp2 = 2*oppo2.representation(*term, coeffs=g1,  axes=[[0, 1, 2], [1, 0, 2]])
        xxpp2 = xxpp2 + xxpp2

        usr = os.path.expanduser('~')
        job_is_dumb = [
            os.path.join(usr, "Documents/Python/config/python3.7/lib/python3.7/"),
            os.path.join(usr, "Documents/UW/Research/Development")
        ]

        quant_states = self.get_states(9, 3, max_quanta=10)
        states = oppo.ravel_state_inds(quant_states)
        # print(quant_states)
        import itertools as ip
        wat = np.array(list(ip.product(states, states))).T
        # with Timer("New style"):
        vals1 = xxpp1[wat[0], wat[1]]
        # with Timer("New style Broadcasting"):
        #     vals2 = xxpp1[np.ix_(states, states)]

        # with Timer("Old style"):
        vals2 = xxpp2[wat[0], wat[1]]
        # with Timer("Old style Broadcasting"):
        #     vals2 = xxpp2[np.ix_(states, states)]

        # import McUtils.Plots as plt
        # n = len(states)
        # plt.ArrayPlot(vals1.reshape((n, n)))
        # plt.ArrayPlot(iphase * vals2.reshape((n, n))).show()

        v1 = vals1
        v2 = iphase * vals2

        self.assertLess(np.max(np.abs(v1 - v2)), 1.0e-14)

    @validationTest
    def test_HOBasis4DPXXP(self):
        from Peeves import Timer, BlockProfiler

        n = 15
        m = 5
        oppo = HarmonicOscillatorProductBasis((n,) * m)
        oppo2 = SimpleProductBasis(HarmonicOscillatorBasis, (n,) * m)

        term = ['p', 'x', 'x', 'p']
        iphase = (-1) ** (term.count("p") // 2)

        xxpp1 = oppo.representation(*term)
        xxpp2 = oppo2.representation(*term)

        usr = os.path.expanduser('~')
        job_is_dumb = [
            os.path.join(usr, "Documents/Python/config/python3.7/lib/python3.7/"),
            os.path.join(usr, "Documents/UW/Research/Development")
        ]

        quant_states = self.get_states(4, m)
        states = oppo.ravel_state_inds(quant_states)
        # print(quant_states)
        import itertools as ip
        wat = np.array(list(ip.product(states, states))).T
        with Timer("New style"):
            vals1 = xxpp1[wat[0], wat[1]]

        with Timer("Old style"):
            vals2 = xxpp2[wat[0], wat[1]]

        v1 = vals1.toarray()
        v2 = iphase * vals2.toarray()

        self.assertLess(np.max(np.abs(v1 - v2)), 1.0e-14)

    @inactiveTest
    def test_HOSelRuleTerms(self):
        from Peeves import Timer, BlockProfiler

        n = 15
        m = 6
        basis = HarmonicOscillatorProductBasis((n,) * m)

        states = BasisStateSpace(
            basis,
            self.get_states(2, m)
        )

        transitions_h1 = [
            [-1],
            [1],
            [-3],
            [3],
            [-1, -1, -1],
            [-1, -1, 1],
            [-1, 1, 1],
            [1, 1, 1],
            [1, 2],
            [-1, 2],
            [1, -2],
            [-1, -2]
        ]

        with BlockProfiler("Selection Rules"):
            h1_space = states.apply_selection_rules(
                transitions_h1,
                1
            )

    # @debugTest
    # def test_HOPXP(self):
    #
    #     new_basi = HarmonicOscillatorProductBasis((10, 10, 10))
    #     old_prod = SimpleProductBasis(HarmonicOscillatorBasis, (10, 10, 10))
    #
    #     new_rep = new_basis.representation('p', 'x', 'p')
    #     old_rep = old_basis.representation('p', 'x', 'p')



