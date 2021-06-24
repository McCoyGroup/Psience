from Peeves import Timer, BlockProfiler

from McUtils.Scaffolding import *
import McUtils.Plots as plt
from McUtils.Combinatorics import CompleteSymmetricGroupSpace
from Peeves.TestUtils import *
from unittest import TestCase
from Psience.BasisReps import *
import sys, os, numpy as np

class BasisSetTests(TestCase):

    def get_states(self, n_quanta, n_modes, max_quanta=None):
        return [np.flip(x) for x in BasisStateSpace.from_quanta(
            HarmonicOscillatorProductBasis(n_modes),
            range(n_quanta)
        ).excitations]

    #region 1D Basis Tests
    @validationTest
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

    @validationTest
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

    @validationTest
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

    @validationTest
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

    @validationTest
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

    @validationTest
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

        n = 10
        m = 2
        oppo = HarmonicOscillatorProductBasis((n,) * m)
        oppo2 = SimpleProductBasis(HarmonicOscillatorBasis, (n,) * m)

        term = ['p', 'p']
        iphase = (-1) ** (term.count("p") // 2)
        n_terms = len(term)

        g1 = np.array(
            [[-1.81146079e-04, 3.97836803e-05],
             [3.97836803e-05, 2.63572358e-05]])
        xxpp1 = 2 * oppo.representation(*term, coeffs=g1, axes=[[0, 1], [1, 0]])
        xxpp1 = xxpp1 + xxpp1
        xxpp2 = 2 * oppo2.representation(*term, coeffs=g1, axes=[[0, 1], [1, 0]])
        xxpp2 = xxpp2 + xxpp2

        usr = os.path.expanduser('~')
        job_is_dumb = [
            os.path.join(usr, "Documents/Python/config/python3.7/lib/python3.7/"),
            os.path.join(usr, "Documents/UW/Research/Development")
        ]

        quant_states = BasisStateSpace(
            oppo,
            self.get_states(9, m, max_quanta=10)
        )
        brakets = quant_states.get_representation_brakets()

        # with Timer("New style"):
            # with BlockProfiler("New Style", strip_dirs=job_is_dumb, num_lines=10, sort_by='tottime', filter="Psience"):
        vals1 = xxpp1[brakets]

        # with Timer("Old style"):
            # with BlockProfiler("Old Style", strip_dirs=job_is_dumb, num_lines=10, sort_by='tottime', filter="Psience"):
        vals2 = xxpp2[brakets]

        v1 = vals1
        v2 = iphase * vals2

        # n = len(quant_states)
        # plt.ArrayPlot(v1.reshape((n, n)))
        # plt.ArrayPlot(v2.reshape((n, n))).show()

        self.assertLess(np.max(np.abs(v1 - v2)), 1.0e-14)

    @validationTest
    def test_HarmHam(self):

        n = 10
        m = 3
        basis = HarmonicOscillatorProductBasis((n,) * m)
        G, V = [
            np.array([[6.47886479e-03, 5.17641431e-12, -1.12922679e-12],
                      [5.17641431e-12, 1.28034398e-02, -3.15629792e-12],
                      [-1.12922679e-12, -3.15629792e-12, 1.76505371e-02]]),
            np.array([[6.47886478e-03, -8.45595180e-13, -1.01327126e-11],
                      [-8.45595549e-13, 1.28034398e-02, -4.72136245e-12],
                      [-1.01327124e-11, -4.72136255e-12, 1.76505372e-02]])]

        mommy = (1 / 2) * basis.representation('p', 'p', coeffs=G)
        possy = (1 / 2) * basis.representation('x', 'x', coeffs=V)
        H0 = ( mommy + possy )

        states = BasisStateSpace(basis, self.get_states(2, 3, max_quanta=10), mode='excitations')

        diag_inds = BraKetSpace(states, states)

        # raise Exception(diag_inds.state_pairs)

        diags = H0[diag_inds]

        self.assertEquals(np.average(diags), 0.036932841734999985)

    @validationTest
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
        xxpp1 = oppo.representation(*term, coeffs=g1, axes=[[0, 1, 2], [1, 0, 2]])
        # xxpp1 = xxpp1 + xxpp1
        xxpp2 = oppo2.representation(*term, coeffs=g1,  axes=[[0, 1, 2], [1, 0, 2]])
        # xxpp2 = xxpp2 + xxpp2

        usr = os.path.expanduser('~')
        job_is_dumb = [
            os.path.join(usr, "Documents/Python/config/python3.7/lib/python3.7/"),
            os.path.join(usr, "Documents/UW/Research/Development")
        ]

        quant_states = BasisStateSpace(
            oppo,
            self.get_states(3, 3, max_quanta=10)
        )
        inds = quant_states.get_representation_brakets()

        #
        # # raise Exception(inds.bras.indices)
        #
        # quant_states = BasisStateSpace(
        #     oppo,
        #     self.get_states(3, 3, max_quanta=10)
        # )
        # new_stuff = quant_states.apply_selection_rules([[-1, 1]])
        # inds2 = new_stuff.get_representation_brakets()
        #
        # plt.ArrayPlot(inds2.adjacency_matrix().toarray()).show()


        # inds = BasisStateSpace(
        #     oppo,
        #     (
        #         [0, 0, 0],
        #         [1, 0, 0],
        #     )
        # ).get_representation_brakets()


        # with Timer("New style"):
        vals1 = xxpp1[inds]
        # with Timer("Old style"):
        vals2 = xxpp2[inds]

        v1 = vals1
        v2 = iphase * vals2

        # n = len(quant_states)
        # plt.ArrayPlot(v1.reshape((n, n)))
        # plt.ArrayPlot(v2.reshape((n, n)))
        # plt.ArrayPlot(v1.reshape((n, n)) - v1.reshape((n, n)).T,
        #               plot_style=dict(vmin=-1.0e-5, vmax=1.0e-5))
        # plt.ArrayPlot(v2.reshape((n, n)) - v2.reshape((n, n)).T,
        #               plot_style=dict(vmin=-1.0e-5, vmax=1.0e-5))
        # plt.ArrayPlot((v1 - v2).reshape((n, n)),
        #               plot_style=dict(vmin=-1.0e-5, vmax=1.0e-5)).show()

        self.assertTrue(
            np.allclose(
                v1[:15],
                [0.00000000e+00, -2.86578374e-04, 0.00000000e+00, 3.29150701e-06,
                 -1.53766049e-04, 0.00000000e+00, -1.59263719e-06, 0.00000000e+00,
                 -5.52442364e-06, 1.24871307e-06, -6.66923918e-05, 0.00000000e+00,
                 -3.81027078e-05, 0.00000000e+00, -1.61862393e-04],
                atol=1.0e-5
            ))

        self.assertLess(np.max(np.abs(v1 - v2)), 1.0e-14)

    @validationTest
    def test_HOBasis3DXXX(self):

        n = 15
        m = 5
        oppo = HarmonicOscillatorProductBasis((n,) * m)
        oppo2 = SimpleProductBasis(HarmonicOscillatorBasis, (n,) * m)

        term = ['x', 'x', 'x']
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
        # with Timer("New style"):
        vals1 = xxpp1[wat[0], wat[1]]

        # with Timer("Old style"):
        vals2 = xxpp2[wat[0], wat[1]]

        v1 = vals1.asarray()
        v2 = iphase * vals2.asarray()

        # with JSONCheckpointer(os.path.expanduser("~/Desktop/test_terms.json")) as chk:
        #     chk['XXX_3D_new'] = v1
        #     chk['XXX_3D_old'] = v2

        self.assertLess(np.max(np.abs(v1 - v2)), 1.0e-14)

    @validationTest
    def test_HOBasis3DXXX2D(self):

        n = 15
        m = 2
        oppo = HarmonicOscillatorProductBasis((n,) * m)
        oppo2 = SimpleProductBasis(HarmonicOscillatorBasis, (n,) * m)

        term = ['x', 'x', 'x']
        iphase = (-1) ** (term.count("p") // 2)

        xxpp1 = oppo.representation(*term)
        xxpp2 = oppo2.representation(*term)

        # usr = os.path.expanduser('~')
        # job_is_dumb = [
        #     os.path.join(usr, "Documents/Python/config/python3.7/lib/python3.7/"),
        #     os.path.join(usr, "Documents/UW/Research/Development")
        # ]

        states = BasisStateSpace.from_quanta(oppo, range(10))
        brakets = states.get_representation_brakets()
        vals1 = xxpp1[brakets]
        vals2 = xxpp2[brakets]

        v1 = vals1.asarray()
        v2 = iphase * vals2.asarray()

        # with JSONCheckpointer(os.path.expanduser("~/Desktop/test_terms.json")) as chk:
        #     chk['XXX_exc'] = states.excitations
        #     chk['XXX_3D_new'] = v1
        #     chk['XXX_3D_old'] = v2

        self.assertLess(np.max(np.abs(v1 - v2)), 2.0e-14)

    @validationTest
    def test_HOBasis3DXXX2DContracted(self):
        n = 15
        m = 2
        oppo = HarmonicOscillatorProductBasis((n,) * m)
        oppo2 = SimpleProductBasis(HarmonicOscillatorBasis, (n,) * m)

        term = ['x', 'x', 'x']
        iphase = (-1) ** (term.count("p") // 2)

        xxpp1 = oppo.representation(*term, coeffs=np.ones((m, m, m)))
        xxpp2 = oppo2.representation(*term, coeffs=np.ones((m, m, m)))

        # usr = os.path.expanduser('~')
        # job_is_dumb = [
        #     os.path.join(usr, "Documents/Python/config/python3.7/lib/python3.7/"),
        #     os.path.join(usr, "Documents/UW/Research/Development")
        # ]

        states = BasisStateSpace.from_quanta(oppo, range(10))
        brakets = states.get_representation_brakets()
        vals1 = xxpp1[brakets]
        vals2 = xxpp2[brakets]

        v1 = vals1
        v2 = iphase * vals2

        with JSONCheckpointer(os.path.expanduser("~/Desktop/test_terms.json")) as chk:
            chk['XXX_exc'] = states.excitations
            chk['XXX_3D_new'] = v1
            chk['XXX_3D_old'] = v2

        self.assertLess(np.max(np.abs(v1 - v2)), 2.0e-14)

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
        # with Timer("New style"):
        vals1 = xxpp1[wat[0], wat[1]]

        # with Timer("Old style"):
        vals2 = xxpp2[wat[0], wat[1]]

        v1 = vals1.toarray()
        v2 = iphase * vals2.toarray()

        self.assertLess(np.max(np.abs(v1 - v2)), 1.0e-14)

    @inactiveTest
    def test_HOSelRuleTerms(self):
        """
        Profiler to see how quickly terms can be generated

        :return:
        :rtype:
        """

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

    @validationTest
    def test_GenerateSelectionRuleSpace(self):
        """
        Tests (and profiles) the generation of a state
        space from a set of selection rules and initial states.
        Mostly here to more easily speed up state space generation
        for use in VPT2.

        :return:
        :rtype:
        """

        basis = HarmonicOscillatorProductBasis(8)
        rules = basis.selection_rules("x", "x", "x", "x")

        states = BasisStateSpace.from_quanta(basis, 3)

        # with BlockProfiler(""):
        h2_space = states.apply_selection_rules(rules, iterations=1)

        self.assertEquals(h2_space.nstates, 120)

    @validationTest
    def test_GenerateFilteredSelectionRuleSpace(self):
        """
        Tests (and profiles) the generation of a state
        space from a set of selection rules and initial states.
        Mostly here to more easily speed up state space generation
        for use in VPT2.

        :return:
        :rtype:
        """

        basis = HarmonicOscillatorProductBasis(8)
        rules = basis.selection_rules("x", "x", "x", "x")

        states = BasisStateSpace.from_quanta(basis, 3)

        h2_space = states.apply_selection_rules(rules, iterations=1)

        sub_h2_space = h2_space.spaces[0].take_subspace(np.arange(10))

        h2_space2 = states.apply_selection_rules(rules,
                                                 filter_space=sub_h2_space,
                                                 iterations=1
                                                 )

        uinds, counts = np.unique(h2_space2.indices, return_counts=True)
        sorting = np.argsort(h2_space2.indices)
        ind_tag = (hash(tuple(uinds)), hash(tuple(counts)), hash(tuple(sorting)))
        # raise Exception(ind_tag)
        self.assertEquals(ind_tag, (320425735722628681, 4044592283957769633))#, 1029650262075525554))
        # raise Exception(ind_tag)

    @validationTest
    def test_StateIndexing(self):
        """
        Tests indexing state specs through a more
        intelligent lexicographic order
        :return:
        :rtype:
        """

        ndim = 6
        indexer = PermutationStateIndexer(ndim)

        states = BasisStateSpace.from_quanta(HarmonicOscillatorProductBasis(ndim), range(10)).excitations
        # print(states)
        inds = indexer.to_indices(states)

        # print(states[44:])

        checks = inds != np.arange(len(states))
        self.assertFalse(
            checks.any()
            , msg="{} no good ({} out of {})".format(states[checks], inds[checks], inds)
        )

        # np.random.seed(0)
        # some_sel = np.arange(len(states))
        some_sel = np.unique(np.random.choice(np.arange(len(states)), 100))
        rev = indexer.from_indices(inds[some_sel,])
        self.assertTrue((states[some_sel,] == rev).all(),
                        msg="{} != {}".format(states[some_sel,], rev))

    # @debugTest
    # def test_HOPXP(self):
    #
    #     new_basi = HarmonicOscillatorProductBasis((10, 10, 10))
    #     old_prod = SimpleProductBasis(HarmonicOscillatorBasis, (10, 10, 10))
    #
    #     new_rep = new_basis.representation('p', 'x', 'p')
    #     old_rep = old_basis.representation('p', 'x', 'p')

    @validationTest
    def test_FindIndices(self):
        ndim = 6
        states = BasisStateSpace.from_quanta(HarmonicOscillatorProductBasis(ndim), range(5))
        test_1 = states.find(states)
        ntest = np.arange(len(test_1))
        self.assertEquals(tuple(test_1), tuple(ntest))

        sel = np.random.choice(ntest, 15)
        _, upos = np.unique(sel, return_index=True)
        sel = sel[np.sort(upos)]
        states2 = states.take_subspace(sel)
        test_2 = states2.find(states2)
        self.assertEquals(tuple(test_2), tuple(np.arange(len(sel))))

    @inactiveTest
    def test_PermIndices(self):
        ndim = 3
        indexer = PermutationStateIndexer(ndim)

        states = [[0, 0, 0],
         [0, 0, 1],
         [0, 0, 2],
         [0, 0, 3],
         [0, 0, 4],
         [0, 0, 5],
         [0, 0, 6],
         [0, 1, 0],
         [0, 1, 1],
         [0, 1, 2],
         [0, 2, 3],
         [0, 2, 4],
         [3, 1, 0],
         [3, 1, 1],
         [3, 1, 2],
         [3, 2, 0],
         [3, 2, 1],
         [3, 3, 0],
         [4, 0, 0],
         [4, 0, 1],
         [4, 0, 2],
         [4, 1, 0],
         [4, 1, 1],
         [4, 2, 0],
         [5, 0, 0],
         [5, 0, 1],
         [5, 1, 0],
         [6, 0, 0]]
        raise Exception(indexer.to_indices(states))

    @inactiveTest
    def test_OperatorAdjacencyGraph(self):
        """
        Tests building an adjacency graph for an operator
        under an initial set of states

        :return:
        :rtype:
        """

        from McUtils.Numputils import SparseArray

        ndim = 4
        basis = HarmonicOscillatorProductBasis(ndim)
        oppo = basis.representation("x", "x", "x", "x", coeffs=np.ones((ndim,)*4))
        rules = basis.selection_rules("x", "x", "x", "x")

        states = BasisStateSpace.from_quanta(basis, 3)
        h2_space = states.apply_selection_rules(rules, iterations=1)
        bk = h2_space.get_representation_brakets()
        flat_total_space = h2_space.to_single().take_unique()

        # pull from get_vpt2_reps or whatever
        # sub = oppo[bk]
        # flat_total_space = h2_space.to_single().take_unique()
        # N = len(flat_total_space)
        #
        # row_inds = flat_total_space.find(bk.bras)
        # col_inds = flat_total_space.find(bk.kets)
        #
        # up_tri = np.array([row_inds, col_inds]).T
        # low_tri = np.array([col_inds, row_inds]).T
        # # but now we need to remove the duplicates, because many sparse matrix implementations
        # # will sum up any repeated elements
        # full_inds = np.concatenate([up_tri, low_tri])
        # full_dat = np.concatenate([sub, sub])
        #
        # _, idx = np.unique(full_inds, axis=0, return_index=True)
        # sidx = np.sort(idx)
        # full_inds = full_inds[sidx]
        # full_dat = full_dat[sidx]
        # adj_mat = SparseArray((full_dat, full_inds.T), shape=(N, N))

        adj_arr = bk.adjacency_matrix(total_space=flat_total_space).toarray()

        # plt.ArrayPlot(adj_arr).show()

        # np.savetxt(os.path.expanduser("~/Desktop/bleh.dat"), adj_arr)

        # adj_dat = adj_mat.data.reshape(N, N)
        # import networkx
        # graph = networkx.from_scipy_sparse_matrix(adj_dat)

    @inactiveTest
    def test_PermIndexingChange(self):
        import json
        ndim = 5
        basis = HarmonicOscillatorProductBasis(ndim, indexer=PermutationStateIndexer(ndim))
        rules = basis.selection_rules("x", "x", "x", "x")

        full_states = BasisStateSpace.from_quanta(basis, 1)
        for x in range(4):
            states = full_states.take_subspace([x])
            # if len(states) == 0:
            #     raise ValueError(states)

            # with BlockProfiler(""):
            h2_space = states.apply_selection_rules(rules, iterations=1)

            print(h2_space)

            print(states.indices.tolist(), np.sort(h2_space.indices).tolist())

    @debugTest
    def test_NewOrthogonalityCalcs(self):

        n = 15
        m = 4
        oppo = HarmonicOscillatorProductBasis((n,) * m)

        states = BasisStateSpace.from_quanta(oppo, range(10))
        brakets = states.get_representation_brakets()
        orthog_1 = brakets.get_non_orthog([0, 0, 1])

        brakets2 = states.get_representation_brakets()
        brakets2.preindex_trie_enabled=False
        brakets2.aggressive_caching_enabled = False
        orthog_2 = brakets2.get_non_orthog([0, 0, 1])

        self.assertTrue( (orthog_1==orthog_2).all() )

        m = 2
        oppo = HarmonicOscillatorProductBasis((n,) * m)

        states = BasisStateSpace.from_quanta(oppo, range(10))

        brakets = states.get_representation_brakets()
        orthog_1 = brakets.get_non_orthog([0, 0, 1])

        brakets2 = states.get_representation_brakets()
        brakets2.preindex_trie_enabled = False
        brakets2.aggressive_caching_enabled = False
        orthog_2 = brakets2.get_non_orthog([0, 0, 1])

        # raise Exception(orthog_1, orthog_2)

        self.assertTrue((orthog_1 == orthog_2).all())

    @inactiveTest
    def test_BasisRepMatrixOps(self):

        n = 15 # totally meaningless these days
        m = 4
        basis = HarmonicOscillatorProductBasis((n,) * m)

        mat = StateSpaceMatrix(basis)

        self.assertEquals(mat.array.shape[0], 0)

        states = BasisStateSpace.from_quanta(basis, range(10))
        brakets = BraKetSpace(states, states)
        vals = np.ones(len(brakets))
        mat_2 = StateSpaceMatrix(brakets, vals)

        def wat(state_space):
            return np.ones(len(state_space))
        sub_brakets = BasisStateSpace.from_quanta(basis, range(4)).get_representation_brakets()
        mat2_vals = mat_2.compute_values(wat, sub_brakets)

        self.assertEquals(mat2_vals.tolist(), mat_2[sub_brakets].tolist())

    @inactiveTest
    def test_ImprovedRepresentations(self):
        n = 15  # totally meaningless these days
        m = 4
        basis = HarmonicOscillatorProductBasis((n,) * m)

        x_rep = basis.representation('x', coeffs=np.ones((m,)))
        initial_space = basis.get_state_space(range(4))
        initial_states = StateSpaceMatrix.identity_from_space(initial_space)

        x = x_rep.apply(initial_states)
        x2 = x_rep.apply(x)
        # plt.ArrayPlot(x.array.asarray(), aspect_ratio=x.shape[0]/x.shape[1], image_size=300)
        # plt.ArrayPlot(x2.array.asarray(), aspect_ratio=x2.shape[0]/x2.shape[1], image_size=300).show()

        x2_rep = basis.representation('x', 'x', coeffs=np.ones((m, m)))
        x22 = x2_rep.apply(initial_states)

        plt.ArrayPlot(x2.array.asarray(), aspect_ratio=x2.shape[0]/x2.shape[1], image_size=200)
        plt.ArrayPlot(x22.array.asarray(), aspect_ratio=x22.shape[0]/x22.shape[1], image_size=200).show()

        self.assertTrue(np.allclose(x2.array.asarray(), x22.array.asarray()))
        # raise Exception(mat_2.array.asarray())

    @validationTest
    def test_StateConnections(self):

        n = 15  # totally meaningless these days
        m = 12
        basis = HarmonicOscillatorProductBasis((n,) * m)

        init_state = basis.get_state_space(2)

        x_rep = basis.representation('x', 'x', 'x', coeffs=np.ones((m, m, m)))
        with BlockProfiler():
            for i in range(15):
                tf = x_rep.get_transformed_space(init_state)
            # tf = x_rep.get_transformed_space(init_state)



        # np.set_printoptions(threshold=100000)
        # raise Exception(tf.excitations)

    @validationTest
    def test_StateSpaceTakeProfile(self):

        n = 15  # totally meaningless these days
        m = 12
        basis = HarmonicOscillatorProductBasis((n,) * m)

        x_rep = basis.representation('x', 'x', coeffs=np.ones((m,)))
        init_space = basis.get_state_space(2)
        init_space.full_basis = CompleteSymmetricGroupSpace(m)
        tf = x_rep.get_transformed_space(init_space)#, logger=Logger())

        # with BlockProfiler():
        #     for i in range(500):
        #         exc = tf.excitations

        # raise Exception(len(exc))

    @validationTest
    def test_PermutationallyReducedStateSpace(self):
        n = 15  # totally meaningless these days
        m = 4
        basis = HarmonicOscillatorProductBasis((n,) * m)

        x_rep = basis.representation('x', 'x', coeffs=np.ones((m,)))
        init_space = basis.get_state_space(2)
        tf = x_rep.get_transformed_space(init_space)

        tf_1 = tf.to_single().take_unique()
        reduced = tf_1.permutationally_reduce()
        expanded = reduced.permutationally_expand()

        self.assertEquals(
            np.sort(tf_1.indices).tolist(),
            np.sort(expanded.indices).tolist()
        )

        sub = reduced.take_subspace([1])
        perms = np.array([[0, 2, 3, 1], [2, 3, 1, 0]])
        direct_prod = sub.permutation_direct_product(perms)
        # raise Exception(
        #     sub.permutationally_expand().excitations,
        #     direct_prod.permutationally_expand().excitations,
        # )

        tf_2 =  x_rep.get_transformed_space(reduced)
        tf_22 = x_rep.get_transformed_space(tf_1)

        new_rep = tf_2.get_space(1).take_permutations(0).permutationally_expand()
        old_rep = tf_22.get_space(0)
        self.assertEquals(
            np.sort(new_rep.indices).tolist(),
            np.sort(old_rep.indices).tolist()
        )

    @validationTest
    def test_TransformedReduced(self):
        n = 15  # totally meaningless these days
        m = 10
        basis = HarmonicOscillatorProductBasis((n,) * m)

        x_rep = basis.representation('x', 'x', 'x', coeffs=np.ones((m, m, m)) )
        init_space = basis.get_state_space(2)

        with BlockProfiler("standard method"):
            tf_space = x_rep.operator.get_transformed_space(init_space).get_representation_brakets()
            tf_els = x_rep.operator.get_elements(tf_space)
        # raise Exception(tf_els)

        with BlockProfiler("new method"):
            tf, bk = x_rep.operator.apply_reduced(init_space)
        self.assertEquals(np.sort(tf).tolist(), np.sort(tf_els).tolist())
