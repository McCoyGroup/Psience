from Peeves.TestUtils import *
from Peeves import Timer, BlockProfiler
from unittest import TestCase
from Psience.VPT2 import *
from McUtils.Data import UnitsData
import sys, os, numpy as np, itertools as ip
import cProfile, pstats, io

class VPTTests(TestCase):

    def setUp(self):
        self.h2w = UnitsData.convert("Hartrees", "Wavenumbers")

    def save_wfns(self, file, wfns):
        """
        We save the corrections so that we can reload them later

        :return:
        :rtype:
        """
        corrs = wfns.corrs
        corrs.savez(file)

    def load_wfns(self, mol, basis, file):
        """
        We save the corrections so that we can reload them later

        :return:
        :rtype:
        """
        from Psience.VPT2.Hamiltonian import PerturbationTheoryCorrections
        return PerturbationTheoryWavefunctions(
            mol,
            basis,
            PerturbationTheoryCorrections.loadz(file)
        )

    wfn_file_dir = os.path.expanduser("~/Desktop/")
    usr = os.path.expanduser('~')
    job_is_dumb = [
        os.path.join(usr, "Documents/Python/config/python3.7/lib/python3.7/"),
        os.path.join(usr, "Documents/UW/Research/Development")
    ]

    def get_VPT2_wfns_and_ham(self,
                              fchk,
                              internals,
                              states,
                              save_coeffs=False,
                              save_wfns=False,
                              regenerate=True,
                              coupled_states=None,
                              mode_selection=None,
                              v2=None,
                              t2=None,
                              v3=None,
                              t3=None,
                              t4=None,
                              v4=None,
                              coriolis=None,
                              watson=None,
                              degeneracies=None,
                              log=False,
                              get_breakdown=False
                              ):

        hammer = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data(fchk),
            internals=internals,
            mode_selection=mode_selection,
            log = log
        )

        if get_breakdown:
            bd = hammer.get_breakdown(states=states, degeneracies=degeneracies, coupled_states=coupled_states)
            return bd, hammer

        wfn_file = os.path.join(self.wfn_file_dir, fchk.replace("fchk", "npz"))
        if regenerate or not os.path.exists(wfn_file):

            if save_coeffs:
                coeffs_file = os.path.join(self.wfn_file_dir, fchk.replace(".fchk", "_coeffs.npz"))
                np.savez(coeffs_file,
                         G=hammer.H0.computers[0].operator.coeffs,
                         F=hammer.H0.computers[1].operator.coeffs,
                         dGdQ=hammer.H1.computers[0].operator.coeffs,
                         dFdQ=hammer.H1.computers[1].operator.coeffs,
                         dGdQQ=hammer.H2.computers[0].operator.coeffs,
                         dFdQQ=hammer.H2.computers[1].operator.coeffs,
                         coriolis=hammer.H2.computers[2].operator.coeffs,
                         watson=hammer.H2.computers[3].operator.coeffs
                         )

            if t2 is not None:
                hammer.H0.computers[0].operator.coeffs = t2
            if v2 is not None:
                hammer.H0.computers[1].operator.coeffs = v2
            if t3 is not None:
                hammer.H1.computers[0].operator.coeffs = t3
            if v3 is not None:
                hammer.H1.computers[1].operator.coeffs = v3
            if t4 is not None:
                hammer.H2.computers[0].operator.coeffs = t4
            if v4 is not None:
                hammer.H2.computers[1].operator.coeffs = v4
            if coriolis is not None:
                hammer.H2.computers[2].operator.coeffs = coriolis
                if watson is None and isinstance(coriolis, (int, float)) and coriolis == 0:
                    watson = 0
            if watson is not None:
                hammer.H2.computers[3].operator.coeffs = watson

            wfns = hammer.get_wavefunctions(states, coupled_states=coupled_states, degeneracies=degeneracies)
            if save_wfns:
                self.save_wfns(wfn_file, wfns)
        else:
            wfns = self.load_wfns(hammer.molecule, hammer.basis, wfn_file)

        return wfns, hammer

    def get_VPT2_wfns(self, fchk,
                      internals,
                      states,
                      save_coeffs=False,
                      regenerate=False,
                      coupled_states=None,
                      mode_selection=None,
                      v2=None,
                      t2=None,
                      v3=None,
                      t3=None,
                      t4=None,
                      v4=None,
                      coriolis=None,
                      watson=None,
                      degeneracies=None,
                      log=False,
                      get_breakdown=False
                      ):
        return self.get_VPT2_wfns_and_ham(
            fchk,
            internals,
            states,
            regenerate=regenerate,
            save_coeffs=save_coeffs,
            coupled_states=coupled_states,
            mode_selection=mode_selection,
            v2=v2,
            t2=t2,
            v3=v3,
            t3=t3,
            t4=t4,
            v4=v4,
            coriolis=coriolis,
            watson=watson,
            degeneracies=degeneracies,
            log=log,
            get_breakdown=get_breakdown
        )[0]
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
    def profile_block(self, block):

        pr = cProfile.Profile()
        pr.enable()
        exc = None
        res = None
        try:
            res = block()
        except Exception as e:
            # we're gonna just throw errors to trace how the program is performing
            # we want to catch this error for later re-throwing, though
            exc = e
        finally:
            pr.disable()
        s = io.StringIO()
        sortby = 'cumulative'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats(50)
        stat_block = s.getvalue()
        usr = os.path.expanduser('~')
        stat_block = stat_block.replace(
            os.path.join(usr, "Documents/Python/config/python3.7/lib/python3.7/"),
            ""
        )
        stat_block = stat_block.replace(
            os.path.join(usr, "Documents/UW/Research/Development"), ""
        )
        return exc, stat_block, res

    #region Test Representations

    @validationTest
    def test_RepresentQQQ(self):
        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("HOD_freq.fchk"),
            internals=None # doesn't matter for this
        )

        QQQ = ham.H1.computers[1].operator
        QQQ.coeffs = 1

        bras = (
            (0, 0, 0),
        )
        kets = (
            (0, 0, 1),
        )

        legit = np.array([
            [[0.        , 0.        , 0.35355339],
             [0.        , 0.        , 0.        ],
             [0.35355339, 0.        , 0.        ]],

            [[0.        , 0.        , 0.        ],
             [0.        , 0.        , 0.35355339],
             [0.        , 0.35355339, 0.        ]],

            [[0.35355339, 0.        , 0.        ],
             [0.        , 0.35355339, 0.        ],
             [0.        , 0.        , 1.06066017]]
        ])

        appx = QQQ[bras, kets].toarray().squeeze()

        self.assertTrue(np.allclose(legit, appx))

    @validationTest
    def test_RepresentPQP(self):
        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("HOD_freq.fchk"),
            internals=None  # doesn't matter for this
        )

        pQp = ham.H1.computers[0].operator
        pQp.coeffs = 1

        bras = (
            (0, 0, 0),
            # (0, 0, 0),
            # (0, 0, 0),
            # (0, 0, 0),
            # (0, 0, 0)
        )
        kets = (
            (0, 0, 1),
            # (0, 0, 2),
            # (0, 1, 1),
            # (1, 0, 1),
            # (1, 1, 0)
        )

        appx = pQp[bras, kets].toarray().transpose(3, 0, 1, 2)

        legit = np.array(
            [[
                [[0.,          0.,        -0.35355339],
                 [0.,          0.,         0.        ],
                 [-0.35355339, 0.,         0.        ]],

                [[0.,          0.,         0.        ],
                 [0.,          0.,        -0.35355339],
                 [0.,         -0.35355339, 0.        ]],

                [[ 0.35355339, 0.,         0.        ],
                 [ 0.,         0.35355339, 0.        ],
                 [ 0.,         0.,        -1.06066   ]]
            ]]
        )

        self.assertTrue(np.allclose(legit, appx))

    @validationTest
    def test_RepresentQQQQ(self):
        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("HOD_freq.fchk"),
            internals=None  # doesn't matter for this
        )

        QQQQ = ham.H2.computers[1].operator
        QQQQ.coeffs = 1

        bras = (
            (0, 0, 0),
        )
        kets = (
            (0, 0, 0),
        )

        legit = np.array([
            [
                [[0.75, 0., 0.],
                 [0., 0.25, 0.],
                 [0., 0., 0.25]],

                [[0., 0.25, 0.],
                 [0.25, 0., 0.],
                 [0., 0., 0.]],

                [[0., 0., 0.25],
                 [0., 0., 0.],
                 [0.25, 0., 0.]]
            ],
            [
                [[0., 0.25, 0.],
                 [0.25, 0., 0.],
                 [0., 0., 0.]],

                [[0.25, 0., 0.],
                 [0., 0.75, 0.],
                 [0., 0., 0.25]],

                [[0., 0., 0.],
                 [0., 0., 0.25],
                 [0., 0.25, 0.]]
            ],
            [
                [[0.,   0., 0.25],
                 [0.,   0., 0.  ],
                 [0.25, 0., 0.  ]],

                [[0., 0., 0.],
                 [0., 0., 0.25],
                 [0., 0.25, 0.]],

                [[0.25, 0., 0.],
                 [0., 0.25, 0.],
                 [0., 0., 0.75]]
            ]
        ])

        # import McUtils.Plots as plt
        # plt.ArrayPlot(
        #     legit.reshape(3*3, 3*3)
        # ).show()

        appx = QQQQ[bras, kets].toarray().squeeze()
        # raise Exception(appx)

        self.assertTrue(np.allclose(legit, appx))

    @validationTest
    def test_RepresentPQQP2(self):
        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("OCHD_freq.fchk"),
            internals=None  # doesn't matter for this
        )

        pQQp = ham.H2.computers[0].operator
        pQQp.coeffs = 1

        bras = (
            (0, 0, 0, 0, 0, 0),
        )
        kets = (
            (0, 0, 0, 0, 0, 0),
        )

        appx = pQQp[bras, kets].toarray().squeeze()

        # import McUtils.Plots as plt
        # plt.ArrayPlot(
        #     appx.reshape(6 * 6, 6 * 6)
        # ).show()

        self.assertAlmostEquals(np.min(appx), -0.75)
        self.assertAlmostEquals(np.max(appx), 0.0)
        self.assertAlmostEquals(appx[0, 0, 0, 0], -0.75)
        self.assertAlmostEquals(appx[0, 0, 1, 1], -0.25)
        self.assertAlmostEquals(appx[0, 1, 0, 1], -0.25)
        self.assertEquals(appx[0, 1, 0, 0], 0.0)

    @validationTest
    def test_RepresentQQQQ2(self):
        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("OCHD_freq.fchk"),
            internals=None  # doesn't matter for this
        )

        QQQQ = ham.H2.computers[1].operator
        QQQQ.coeffs = 1

        bras = (
            (0, 0, 0, 0, 0, 0),
        )
        kets = (
            (0, 0, 0, 0, 0, 0),
        )

        appx = QQQQ[bras, kets].toarray().squeeze()

        self.assertAlmostEquals(np.min(appx), 0.0)
        self.assertAlmostEquals(np.max(appx),  .75)
        self.assertAlmostEquals(appx[0, 0, 0, 0],  .75)
        self.assertAlmostEquals(appx[0, 0, 1, 1],  .25)
        self.assertAlmostEquals(appx[0, 1, 0, 1],  .25)
        self.assertEquals(appx[0, 1, 0, 0], 0)

        # import McUtils.Plots as plt
        # plt.ArrayPlot(
        #     appx.reshape(6 * 6, 6 * 6)
        # ).show()

    @validationTest
    def test_RepresentPQP(self):
        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("HOD_freq.fchk"),
            internals=None  # doesn't matter for this
        )

        pQp = ham.H1.computers[0].operator
        pQp.coeffs = 1

        bras = (
            (0, 0, 0),
            # (0, 0, 0),
            # (0, 0, 0),
            # (0, 0, 0),
            # (0, 0, 0)
        )
        kets = (
            (0, 0, 1),
            # (0, 0, 2),
            # (0, 1, 1),
            # (1, 0, 1),
            # (1, 1, 0)
        )

        appx = pQp[bras, kets].toarray().transpose(3, 0, 1, 2)

        legit = np.array(
            [[[[0., 0., 0.35355339],
               [0., 0., 0.],
               [-0.35355339, 0., 0.]],

              [[0., 0., 0.],
               [0., 0., 0.35355339],
               [0., -0.35355339, 0.]],

              [[-0.35355339, 0., 0.],
               [0., -0.35355339, 0.],
               [0., 0., -0.35355339]]]]
        )

        self.assertTrue(np.allclose(legit, appx))

    @validationTest
    def test_RepresentFullQQQQ(self):

        n_modes = 6
        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("OCHD_freq.fchk"),
            internals=None  # doesn't matter for this
        )

        pp = ham.H0.computers[0].operator
        pp.coeffs = 0

        QQ = ham.H0.computers[1].operator
        QQ.coeffs = 0

        pQp = ham.H1.computers[0].operator
        pQp.coeffs = 0

        QQQ = ham.H2.computers[1].operator
        QQQ.coeffs = 0

        pQQp = ham.H2.computers[0].operator
        pQQp.coeffs = 0

        QQQQ = ham.H2.computers[1].operator
        QQQQ.coeffs = 1

        coupled_states = self.get_states(4, n_modes, max_quanta=4)
        h0, h1, h2 = ham.get_representations(coupled_states)

        sub = 24*h2[0, 0, :, :, 0, 0].toarray()
        expected = np.array([
            [0.75, 0.0,  0.0,  0.0,  0.0,  0.0 ],
            [0.0,  0.25, 0.0,  0.0,  0.0,  0.0 ],
            [0.0,  0.0,  0.25, 0.0,  0.0,  0.0 ],
            [0.0,  0.0,  0.0,  0.25, 0.0,  0.0 ],
            [0.0,  0.0,  0.0,  0.0,  0.25, 0.0 ],
            [0.0,  0.0,  0.0,  0.0,  0.0,  0.25]
        ])

        self.assertTrue(np.allclose(sub, expected))

    @validationTest
    def test_RepresentFullPQQP(self):

        n_modes = 6
        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("OCHD_freq.fchk"),
            internals=None  # doesn't matter for this
        )

        pp = ham.H0.computers[0].operator
        pp.coeffs = 0

        QQ = ham.H0.computers[1].operator
        QQ.coeffs = 0

        pQp = ham.H1.computers[0].operator
        pQp.coeffs = 0

        QQQ = ham.H2.computers[1].operator
        QQQ.coeffs = 0

        pQQp = ham.H2.computers[0].operator
        pQQp.coeffs = 1

        QQQQ = ham.H2.computers[1].operator
        QQQQ.coeffs = 0

        coupled_states = self.get_states(4, n_modes, max_quanta=4)
        h0, h1, h2 = ham.get_representations(coupled_states)

        sub = 4 * h2[0, 0, :, :, 0, 0].toarray()
        expected = np.array([
            [0.75, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.25, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.25, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.25, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.25, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.25]
        ])

        self.assertTrue(np.allclose(sub, expected))

    #endregion

    #region Test Inputs

    @validationTest
    def test_TestQuarticsCartesians(self):
        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("HOD_freq.fchk"),
            internals=None
        )

        v4_Gaussian = [
            [1, 1, 1, 1, 1517.96213],
            [2, 1, 1, 1, 98.59961],
            [2, 2, 1, 1, 8.99887],
            [2, 2, 2, 1, -50.96655],
            [2, 2, 2, 2, 804.29611],
            [3, 1, 1, 1, 142.08091],
            [3, 2, 1, 1, -18.73606],
            [3, 2, 2, 1, -22.35470],
            [3, 2, 2, 2, 71.81011],
            [3, 3, 1, 1, -523.38920],
            [3, 3, 2, 1, -4.05652],
            [3, 3, 2, 2, -95.43623],
            [3, 3, 3, 1, -145.84374],
            [3, 3, 3, 2, -41.06991],
            [3, 3, 3, 3, 83.41603]
        ]

        legit = np.zeros((3, 3, 3, 3))
        mode_mapping = [
            2, 1, 0
        ]
        for i, j, k, l, v in v4_Gaussian:
            i = mode_mapping[i - 1];
            j = mode_mapping[j - 1]
            k = mode_mapping[k - 1];
            l = mode_mapping[l - 1]
            for perm in ip.permutations((i, j, k, l)):
                legit[perm] = v

        v4 = self.h2w * ham.V_terms[2]

        print_errors = True
        if print_errors:
            if not np.allclose(legit, v4, atol=.001):
                diff = legit - v4
                bad_pos = np.array(np.where(np.abs(diff) > .001)).T
                print("Gaussian/This Disagreements:\n" + "\n".join(
                    "{:>.0f} {:>.0f} {:>.0f} {:>.0f} {:>5.3f} (Actual: {:>8.3f} This: {:>8.3f})".format(
                        i, j, k, l, diff[i, j, k, l], legit[i, j, k, l], v4[i, j, k, l]
                    ) for i, j, k, l in bad_pos
                ))

        self.assertTrue(np.allclose(legit, v4, atol=.1))  # testing to within .001 wavenumbers

    @inactiveTest
    def test_TestCubicsInternals(self):
        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("HOD_freq.fchk"),
            internals=[
                [0, -1, -1, -1],
                [1, 0, -1, -1],
                [2, 0, 1, -1]
            ]
        )

        # it turns out Anne and I disagree pretty fucking dramatically on these...
        # but is it an issue? If I use her cubic force constants the energies are way, way off...
        # RESOLVED: this is from an old bad version of Anne's code
        v4_Anne = [
            [1, 1, 1,   198.63477267],
            [2, 1, 1,    41.05987944],
            [3, 1, 1,   429.45742955],
            [1, 2, 1,    41.05987944],
            [2, 2, 1,   -90.66588863],
            [3, 2, 1,    82.31540784],
            [1, 3, 1,   429.45742955],
            [2, 3, 1,    82.31540784],
            [3, 3, 1,  -153.50694953],
            [1, 1, 2,    41.16039749],
            [2, 1, 2,   -90.73683527],
            [3, 1, 2,    82.32690754],
            [1, 2, 2,   -90.73683527],
            [2, 2, 2, -1588.89967595],
            [3, 2, 2, 91.00951548],
            [1, 3, 2, 82.32690754],
            [2, 3, 2, 91.00951548],
            [3, 3, 2, -172.34377091],
            [1, 1, 3, 430.44067645],
            [2, 1, 3, 82.33125106],
            [3, 1, 3, -153.78923868],
            [1, 2, 3, 82.33125106],
            [2, 2, 3, 90.96564514],
            [3, 2, 3, -172.45333790],
            [1, 3, 3, -153.78923868],
            [2, 3, 3, -172.45333790],
            [3, 3, 3, -2558.24567375]
        ]

        legit = np.zeros((3, 3, 3))
        for k, j, i, v in v4_Anne:
            i = i - 1;
            j = j - 1;
            k = k - 1
            legit[i, j, k] = v

        v3 = self.h2w * ham.V_terms[1]

        # for unknown reasons, Anne and I disagree on the order of like 5 cm^-1,
        # but it's unclear which set of derivs. is right since my & hers both
        # yield a max deviation from the Gaussian result of ~.1 cm^-1
        print_errors = True
        if print_errors:
            if not np.allclose(legit, v3, atol=.1):
                diff = legit - v3
                bad_pos = np.array(np.where(np.abs(diff) > .1)).T
                print("Anne/This Disagreements:\n" + "\n".join(
                    "{:>.0f} {:>.0f} {:>.0f} {:>8.3f} (Anne: {:>10.3f} This: {:>10.3f})".format(
                        i, j, k, diff[i, j, k], legit[i, j, k], v3[i, j, k]
                    ) for i, j, k in bad_pos
                ))

        self.assertTrue(np.allclose(legit, v3, atol=10))

    @validationTest
    def test_TestQuarticsInternals(self):
        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("HOD_freq.fchk"),
            internals=[
                [0, -1, -1, -1],
                [1,  0, -1, -1],
                [2,  0,  1, -1]
            ]
        )

        v4_Anne = [
            [1, 1, 1,  -37.03937000],
            [2, 1, 1,  -32.30391126],
            [3, 1, 1,  -33.08215609],
            [1, 2, 1,  -32.30391126],
            [2, 2, 1,    3.57147725],
            [3, 2, 1,    9.77124742],
            [1, 3, 1,  -33.08215609],
            [2, 3, 1,    9.77124742],
            [3, 3, 1,    3.08396862],
            [1, 1, 2,    3.53204514],
            [2, 1, 2,   66.35374213],
            [3, 1, 2,   -8.46713126],
            [1, 2, 2,   66.35374213],
            [2, 2, 2,  804.47871323],
            [3, 2, 2,  -51.44004640],
            [1, 3, 2,   -8.46713126],
            [2, 3, 2,  -51.44004640],
            [3, 3, 2,   10.60086681],
            [1, 1, 3,    2.67361974],
            [2, 1, 3,    3.14497676],
            [3, 1, 3,  111.80682105],
            [1, 2, 3,    3.14497676],
            [2, 2, 3,   10.60153758],
            [3, 2, 3,   97.05643377],
            [1, 3, 3,  111.80682105],
            [2, 3, 3,   97.05643377],
            [3, 3, 3, 1519.00602277]
        ]

        legit = np.zeros((3, 3, 3, 3))
        for k, j, i, v in v4_Anne:
            i = i - 1;
            j = j - 1;
            k = k - 1
            for perm in ip.permutations((i, i, j, k)):
                legit[perm] = v

        v4 = self.h2w * ham.V_terms[2]

        print_errors = True
        if print_errors:
            if not np.allclose(legit, v4, atol=0):
                diff = legit - v4
                bad_pos = np.array(np.where(np.abs(diff) > 0)).T
                print("Anne/This Disagreements:\n" + "\n".join(
                    "{:>.0f} {:>.0f} {:>.0f} {:>.0f} {:>8.3f} (Anne: {:>10.3f} This: {:>10.3f})".format(
                        i, j, k, l, diff[i, j, k, l], legit[i, j, k, l], v4[i, j, k, l]
                    ) for i, j, k, l in bad_pos
                ))

        self.assertTrue(np.allclose(legit, v4, atol=1))

    @validationTest
    def test_TestCubicsInternalsHOH(self):
        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("HOH_freq.fchk"),
            internals=[
                [0, -1, -1, -1],
                [1, 0, -1, -1],
                [2, 0, 1, -1]
            ]
        )

        v4_Anne = [
        [1,   1,   1,    323.8035817067],
        [2,   1,   1,   -132.0532052153],
        [1,   2,   1,   -132.0532052153],
        [1,   1,   2,   -131.6177927698],

        [2,   2,   1,   -117.6853050409],
        [2,   1,   2,   -117.8181841885],
        [1,   2,   2,   -117.8181841885],

        [3,   3,   1,   -214.1124120699],
        [3,   1,   3,   -214.2976999111],
        [1,   3,   3,   -214.2976999111],

        [2,   2,   2,  -1829.4315686310],

        [3,   3,   2,  -1812.9480307020],
        [3,   2,   3,  -1812.9726116114],
        [2,   3,   3,  -1812.9726116114]
        ]

        legit = np.zeros((3, 3, 3))
        for k, j, i, v in v4_Anne:
            i = i - 1;
            j = j - 1;
            k = k - 1
            legit[i, j, k] = v

        v3 = self.h2w * ham.V_terms[1]

        print_errors = False
        if print_errors:
            if not np.allclose(legit, v3, atol=.1):
                diff = legit - v3
                bad_pos = np.array(np.where(np.abs(diff) > .1)).T
                print("Anne/This Disagreements:\n" + "\n".join(
                    "{:>.0f} {:>.0f} {:>.0f} {:>8.3f} (Anne: {:>10.3f} This: {:>10.3f})".format(
                        i, j, k, diff[i, j, k], legit[i, j, k], v3[i, j, k]
                    ) for i, j, k in bad_pos
                ))

        self.assertTrue(np.allclose(legit, v3, atol=1))

    @validationTest
    def test_TestQuarticsInternalsHOH(self):
        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("HOH_freq.fchk"),
            internals=[
                [0, -1, -1, -1],
                [1,  0, -1, -1],
                [2,  0,  1, -1]
            ]
        )

        v4_Anne = [
            [1, 1, 1, 1,  -49.5572026693],

            [2, 2, 1, 1,   19.9982906121],
            [1, 1, 2, 2,   19.8598962862],

            [3, 3, 1, 1,   -7.6481006128],
            [1, 1, 3, 3,   -7.8037499292],

            [2, 2, 2, 2,  768.5651528907],

            [3, 3, 2, 2,  763.1648005165],
            [2, 2, 3, 3,  763.1728458966],

            [3, 3, 3, 3,  764.3609369845]
        ]

        legit = np.zeros((3, 3, 3, 3))
        for i, j, k, l, v in v4_Anne:
            i = i - 1; j = j - 1; k = k - 1; l = l - 1
            for perm in ip.permutations((i, j, k, l)):
                legit[perm] = v

        v4 = self.h2w * ham.V_terms[2]
        for i, j, k, l in ip.product(range(3), range(3), range(3), range(3)):
            if i == j:
                if k != l:
                    v4[i, j, k, l] = 0.
            elif i == k:
                if j != l:
                    v4[i, j, k, l] = 0.
            elif i == l:
                if j != k:
                    v4[i, j, k, l] = 0.
            else:
                v4[i, j, k, l] = 0.


        print_errors = True
        if print_errors:
            if not np.allclose(legit, v4, atol=1):
                diff = legit - v4
                bad_pos = np.array(np.where(np.abs(diff) > 1)).T
                print("Anne/This Disagreements:\n" + "\n".join(
                    "{:>.0f} {:>.0f} {:>.0f} {:>.0f} {:>8.3f} (Anne: {:>10.3f} This: {:>10.3f})".format(
                        i, j, k, l, diff[i, j, k, l], legit[i, j, k, l], v4[i, j, k, l]
                    ) for i, j, k, l in bad_pos
                ))

        self.assertTrue(np.allclose(legit, v4, atol=1))

    @validationTest
    def test_TestGQInternalsHOH(self):
        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("HOH_freq.fchk"),
            internals=[
                [0, -1, -1, -1],
                [1,  0, -1, -1],
                [2,  0,  1, -1]
            ]
        )

        v3_Anne = [
        [1,   1,   1,   -42.7571572306],

        [2,   1,   1,    11.2415784419],
        [1,   2,   1,    11.2415784419],
        [1,   1,   2,  -236.4538948819],

        [2,   2,   1,    45.2894929663],
        [2,   1,   2,    14.0668758162],
        [1,   2,   2,    14.0668758162],

        [3,   3,   1,   -47.2340749349],
        [3,   1,   3,    10.2421465679],
        [1,   3,   3,    10.2421465679],

        [2,   2,   2,    -0.6940443704],

        [3,   3,   2,     0.3058655908],
        [3,   2,   3,    -1.0551886296],
        [2,   3,   3,    -1.0551886296]
        ]

        legit = np.zeros((3, 3, 3))
        for k, j, i, v in v3_Anne:
            i = i - 1
            j = j - 1
            k = k - 1
            legit[i, j, k] = v

        v3 = self.h2w * ham.G_terms[1]

        print_errors = False
        if print_errors:
            if not np.allclose(legit, v3, atol=1):
                diff = legit - v3
                bad_pos = np.array(np.where(np.abs(diff) > 1)).T
                print("Anne/This Disagreements:\n" + "\n".join(
                    "{:>.0f} {:>.0f} {:>.0f} {:>8.3f} (Anne: {:>10.3f} This: {:>10.3f})".format(
                        i, j, k, diff[i, j, k], legit[i, j, k], v3[i, j, k]
                    ) for i, j, k in bad_pos
                ))

        self.assertTrue(np.allclose(legit, v3, atol=10))

    @validationTest
    def test_TestGQQInternalsHOH(self):
        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("HOH_freq.fchk"),
            internals=[
                [0, -1, -1, -1],
                [1,  0, -1, -1],
                [2,  0,  1, -1]
            ]
        )

        # it turns out Anne's ij/ij numbers are questionable
        v4_Anne = [
            [1, 1, 1, 1, -0.0813253351],
            [2, 2, 1, 1,  3.9823166834],
            # [2, 1, 2, 1, 27.7264838055],
            # [1, 2, 1, 2, 27.7264838055],
            [1, 1, 2, 2, 49.3985875508],
            [3, 3, 1, 1, -2.5877826647],
            # [3, 1, 3, 1, 20.1052292927],
            # [1, 3, 1, 3, 20.1052292927],
            [1, 1, 3, 3, 48.8826338651],
            [2, 2, 2, 2,  0.2270703177],
            [3, 3, 2, 2, -0.0001049353],
            # [3, 2, 3, 2, -1.6145992219],
            # [2, 3, 2, 3, -1.6145992219],
            [2, 2, 3, 3,  0.2233861641],
            [3, 3, 3, 3, -0.0000015626]
        ]

        legit = np.zeros((3, 3, 3, 3))
        for l, k, j, i, v in v4_Anne:
            i = i - 1
            j = j - 1
            k = k - 1
            l = l - 1
            legit[i, j, k, l] = v

        v4 = self.h2w * ham.G_terms[2]
        import itertools as ip
        for i, j, k, l in ip.product(range(3), range(3), range(3), range(3)):
            s = list(sorted([i, j, k, l]))
            if s[0] != s[1] or s[2] != s[3]: # Anne doesn't have numbers for these values
                # print(s)
                v4[i, j, k, l] = 0.

        # import McUtils.Plots as plt
        #
        # plt.TensorPlot(legit, plot_style={'vmin':-30, 'vmax':30})
        # plt.TensorPlot(v4, plot_style={'vmin':-30, 'vmax':30}).show()

        print_errors = True
        if print_errors:
            if not np.allclose(legit, v4, atol=5):
                diff = legit - v4
                bad_pos = np.array(np.where(np.abs(diff) > .5)).T
                print("Anne/This Disagreements:\n" + "\n".join(
                    "{:>.0f} {:>.0f} {:>.0f} {:>.0f} {:>8.3f} (Anne: {:>10.3f} This: {:>10.3f})".format(
                        i, j, k, l, diff[i, j, k, l], legit[i, j, k, l], v4[i, j, k, l]
                    ) for i, j, k, l in bad_pos
                ))

        self.assertTrue(np.allclose(legit, v4, atol=5))

    @validationTest
    def test_TestCubicsCartesians2(self):
        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("OCHD_freq.fchk"),
            internals=None
        )

        v3_Gaussian = [
                [2,  1,  1,     363.57094],
                [2,  2,  2,   -1975.93951],
                [3,  1,  1,     -96.95870],
                [3,  2,  2,     121.92201],
                [3,  3,  2,     108.24195],
                [3,  3,  3,    1222.60350],
                [4,  1,  1,     -18.34694],
                [4,  2,  2,      55.71496],
                [4,  3,  2,      37.03147],
                [4,  3,  3,    -119.05324],
                [4,  4,  2,      38.87584],
                [4,  4,  3,    -107.82518],
                [4,  4,  4,    -551.81346],
                [5,  1,  1,      36.64257],
                [5,  2,  2,     -66.26172],
                [5,  3,  2,     -66.45267],
                [5,  3,  3,     -19.01573],
                [5,  4,  2,     -88.00773],
                [5,  4,  3,     -16.16525],
                [5,  4,  4,     -84.27197],
                [5,  5,  2,     235.75713],
                [5,  5,  3,      12.71311],
                [5,  5,  4,     -77.86211],
                [5,  5,  5,     -34.43596],
                [6,  1,  1,      26.75955],
                [6,  2,  2,     -31.06476],
                [6,  3,  2,     -18.05063],
                [6,  3,  3,     -41.63776],
                [6,  4,  2,      58.71180],
                [6,  4,  3,      59.11324],
                [6,  4,  4,     -38.22535],
                [6,  5,  2,     -96.15434],
                [6,  5,  3,      -9.13930],
                [6,  5,  4,     -14.76207],
                [6,  5,  5,      78.97040],
                [6,  6,  2,      -0.55703],
                [6,  6,  3,    -153.39899],
                [6,  6,  4,     -15.21350],
                [6,  6,  5,     -18.10616],
                [6,  6,  6,     -58.56657]
        ]

        legit = np.zeros((6, 6, 6))
        mode_mapping = [1, 5, 4, 3, 2, 0]
        for i, j, k, v in v3_Gaussian:
            i = mode_mapping[i - 1]
            j = mode_mapping[j - 1]
            k = mode_mapping[k - 1]
            for perm in ip.permutations((i, j, k)):
                legit[perm] = v

        v3 = self.h2w * ham.V_terms[1]

        print_errors = True
        if print_errors:
            if not np.allclose(legit, v3, atol=1):
                diff = legit - v3
                bad_pos = np.array(np.where(np.abs(diff) > 1)).T
                print("Gaussian/This Disagreements:\n" + "\n".join(
                    "{:>.0f} {:>.0f} {:>.0f} {:>.0f} {:>5.3f} (Actual: {:>8.3f} This: {:>8.3f})".format(
                        i, j, k, diff[i, j, k], legit[i, j, k], v3[i, j, k]
                    ) for i, j, k in bad_pos
                ))

        self.assertTrue(np.allclose(legit, v3, atol=1))  # testing to within a wavenumber

    @validationTest
    def test_TestQuarticsCartesians2(self):
        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("OCHD_freq.fchk"),
            internals=None
        )

        v4_Gaussian = [
            [1,  1,  1,  1,     164.58095],
            [2,  2,  1,  1,    -417.46602],
            [2,  2,  2,  2,    1103.84367],
            [3,  2,  1,  1,       3.04182],
            [3,  2,  2,  2,     -60.94479],
            [3,  3,  1,  1,    -106.04643],
            [3,  3,  2,  2,       3.61805],
            [3,  3,  3,  2,      50.31773],
            [3,  3,  3,  3,     578.95057],
            [4,  2,  1,  1,      20.36969],
            [4,  2,  2,  2,     -29.02635],
            [4,  3,  1,  1,       3.47924],
            [4,  3,  2,  2,     -13.46981],
            [4,  3,  3,  2,      -0.91433],
            [4,  3,  3,  3,     -52.25043],
            [4,  4,  1,  1,     -13.93558],
            [4,  4,  2,  2,     -39.88255],
            [4,  4,  3,  2,      -8.24727],
            [4,  4,  3,  3,      -7.24933],
            [4,  4,  4,  2,      -2.23878],
            [4,  4,  4,  3,      37.30299],
            [4,  4,  4,  4,     160.67711],
            [5,  2,  1,  1,     -13.89043],
            [5,  2,  2,  2,      62.03997],
            [5,  3,  1,  1,      -1.50577],
            [5,  3,  2,  2,      43.97280],
            [5,  3,  3,  2,      -7.25862],
            [5,  3,  3,  3,      -6.99604],
            [5,  4,  1,  1,     -13.33270],
            [5,  4,  2,  2,     114.53101],
            [5,  4,  3,  3,      14.40884],
            [5,  4,  4,  2,     -11.44982],
            [5,  4,  4,  3,       5.42367],
            [5,  4,  4,  4,      29.13016],
            [5,  5,  1,  1,      31.34711],
            [5,  5,  2,  2,    -361.04782],
            [5,  5,  3,  2,       8.83851],
            [5,  5,  3,  3,      -9.46934],
            [5,  5,  4,  2,      28.19922],
            [5,  5,  4,  3,       5.08660],
            [5,  5,  4,  4,      13.43015],
            [5,  5,  5,  2,     -14.95264],
            [5,  5,  5,  3,     -20.24561],
            [5,  5,  5,  4,     -19.14598],
            [5,  5,  5,  5,      65.90101],
            [6,  2,  1,  1,      -6.68214],
            [6,  2,  2,  2,       2.95096],
            [6,  3,  1,  1,       6.40218],
            [6,  3,  2,  2,     -16.45996],
            [6,  3,  3,  2,     -16.44845],
            [6,  3,  3,  3,     -28.58871],
            [6,  4,  1,  1,       4.50738],
            [6,  4,  2,  2,     -28.51506],
            [6,  4,  3,  3,      53.07429],
            [6,  4,  4,  2,      -3.57254],
            [6,  4,  4,  3,      -9.09429],
            [6,  4,  4,  4,      11.54065],
            [6,  5,  1,  1,     -17.18941],
            [6,  5,  2,  2,      97.35767],
            [6,  5,  3,  3,     -37.37982],
            [6,  5,  4,  4,       4.02941],
            [6,  5,  5,  2,     -15.86701],
            [6,  5,  5,  3,      12.25290],
            [6,  5,  5,  4,       7.18948],
            [6,  5,  5,  5,     -30.32976],
            [6,  6,  1,  1,       1.57620],
            [6,  6,  2,  2,     -32.47267],
            [6,  6,  3,  2,     -14.65910],
            [6,  6,  3,  3,    -181.29487],
            [6,  6,  4,  2,      -0.62290],
            [6,  6,  4,  3,       0.80724],
            [6,  6,  4,  4,      -5.26763],
            [6,  6,  5,  2,      12.21672],
            [6,  6,  5,  3,       4.09931],
            [6,  6,  5,  4,      -3.58662],
            [6,  6,  5,  5,       8.33588],
            [6,  6,  6,  2,       9.67099],
            [6,  6,  6,  3,       4.36194],
            [6,  6,  6,  4,      -5.35107],
            [6,  6,  6,  5,      -5.93064],
            [6,  6,  6,  6,      42.05633]
        ]

        legit = np.zeros((6, 6, 6, 6))
        mode_mapping = [1, 5, 4, 3, 2, 0]
        for i, j, k, l, v in v4_Gaussian:
            i = mode_mapping[i - 1]
            j = mode_mapping[j - 1]
            k = mode_mapping[k - 1]
            l = mode_mapping[l - 1]
            for perm in ip.permutations((i, j, k, l)):
                legit[perm] = v

        v4 = self.h2w * ham.V_terms[2]

        print_errors = True
        if print_errors:
            if not np.allclose(legit, v4, atol=1):
                diff = legit - v4
                bad_pos = np.array(np.where(np.abs(diff) > 1)).T
                print("Gaussian/This Disagreements:\n" + "\n".join(
                    "{:>.0f} {:>.0f} {:>.0f} {:>.0f} {:>5.3f} (Actual: {:>8.3f} This: {:>8.3f})".format(
                        i, j, k, l, diff[i, j, k, l], legit[i, j, k, l], v4[i, j, k, l]
                    ) for i, j, k, l in bad_pos
                ))

        self.assertTrue(np.allclose(legit, v4, atol=1))  # testing to within a wavenumber

    @validationTest
    def test_TestCubicsCartesians3(self):
        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("CH2DT_freq.fchk"),
            internals=None
        )

        v3_Gaussian = [
            [1,  1,  1,       0.06815],
            [2,  1,  1,   -1395.70040],
            [2,  2,  2,   -1333.19164],
            [3,  1,  1,     -61.21371],
            [3,  2,  2,     -72.75916],
            [3,  3,  2,     111.94447],
            [3,  3,  3,   -1181.91933],
            [4,  1,  1,      47.93552],
            [4,  2,  2,      64.19649],
            [4,  3,  2,      10.38925],
            [4,  3,  3,     199.69522],
            [4,  4,  2,      74.18204],
            [4,  4,  3,     152.48473],
            [4,  4,  4,     909.14866],
            [5,  1,  1,    -103.49761],
            [5,  2,  2,      30.68413],
            [5,  3,  2,     -24.92683],
            [5,  3,  3,       8.92570],
            [5,  4,  2,      33.14072],
            [5,  4,  3,      -1.79770],
            [5,  4,  4,      -9.43933],
            [5,  5,  2,     195.91273],
            [5,  5,  3,      -2.37845],
            [5,  5,  4,       4.81627],
            [5,  5,  5,      67.87295],
            [6,  2,  1,     -15.11274],
            [6,  3,  1,      67.34301],
            [6,  4,  1,      43.83666],
            [6,  5,  1,     -34.74422],
            [6,  6,  2,     179.13158],
            [6,  6,  3,       3.08917],
            [6,  6,  4,      17.43233],
            [6,  6,  5,     -33.45784],
            [7,  1,  1,     -11.74252],
            [7,  2,  2,      -5.74352],
            [7,  3,  2,      65.00439],
            [7,  3,  3,     -22.17977],
            [7,  4,  2,      56.39770],
            [7,  4,  3,      11.40643],
            [7,  4,  4,     -11.28175],
            [7,  5,  2,       2.60202],
            [7,  5,  3,      25.24120],
            [7,  5,  4,      18.82210],
            [7,  5,  5,      -1.05280],
            [7,  6,  1,     244.67887],
            [7,  6,  6,      13.00576],
            [7,  7,  2,     338.17722],
            [7,  7,  3,      17.93616],
            [7,  7,  4,      15.22551],
            [7,  7,  5,      55.22135],
            [7,  7,  7,      32.37393],
            [8,  2,  1,      76.86835],
            [8,  3,  1,     -23.32698],
            [8,  4,  1,      65.32973],
            [8,  5,  1,     247.14911],
            [8,  6,  2,      34.76624],
            [8,  6,  3,     -33.94786],
            [8,  6,  4,      13.86767],
            [8,  6,  5,       8.61294],
            [8,  7,  1,      81.98867],
            [8,  7,  6,     -28.48111],
            [8,  8,  2,     182.23674],
            [8,  8,  3,      88.65287],
            [8,  8,  4,     -77.54773],
            [8,  8,  5,     -76.75741],
            [8,  8,  7,     -12.69855],
            [9,  1,  1,     -54.45594],
            [9,  2,  2,     -30.24115],
            [9,  3,  2,       1.30117],
            [9,  3,  3,      33.58217],
            [9,  4,  2,     -19.51444],
            [9,  4,  3,      53.68501],
            [9,  4,  4,     -18.44456],
            [9,  5,  2,      78.51601],
            [9,  5,  3,      -1.11472],
            [9,  5,  4,     -18.35336],
            [9,  5,  5,      13.22045],
            [9,  6,  1,     -57.39319],
            [9,  6,  6,      28.70793],
            [9,  7,  2,     -56.58162],
            [9,  7,  3,      66.81247],
            [9,  7,  4,      21.64467],
            [9,  7,  5,     -11.02663],
            [9,  7,  7,      73.63799],
            [9,  8,  1,      82.30035],
            [9,  8,  6,      21.26192],
            [9,  8,  8,     -49.88854],
            [9,  9,  2,      10.73372],
            [9,  9,  3,     122.77338],
            [9,  9,  4,     -77.62983],
            [9,  9,  5,      15.46243],
            [9,  9,  7,     -14.66681],
            [9,  9,  9,     -23.60486]
        ]

        legit = np.zeros((9, 9, 9))
        mode_mapping = list(reversed(range(9)))
        for i, j, k, v in v3_Gaussian:
            i = mode_mapping[i - 1]
            j = mode_mapping[j - 1]
            k = mode_mapping[k - 1]
            for perm in ip.permutations((i, j, k)):
                legit[perm] = v

        v3 = self.h2w * ham.V_terms[1]

        print_errors = True
        if print_errors:
            if not np.allclose(legit, v3, atol=1):
                diff = legit - v3
                bad_pos = np.array(np.where(np.abs(diff) > 1)).T
                print("Gaussian/This Disagreements:\n" + "\n".join(
                    "{:>.0f} {:>.0f} {:>.0f} {:>8.3f} (Gaussian: {:>8.3f} This: {:>8.3f})".format(
                        i, j, k, diff[i, j, k], legit[i, j, k], v3[i, j, k]
                    ) for i, j, k in bad_pos
                ))

        self.assertTrue(np.allclose(legit, v3, atol=1))  # testing to within a wavenumber

    @validationTest
    def test_TestQuarticsCartesians3(self):
        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("CH2DT_freq.fchk"),
            internals=None
        )

        v4_Gaussian = [
             [1,  1,  1,  1,     559.91879],
             [2,  1,  1,  1,      -0.01425],
             [2,  2,  1,  1,     534.37271],
             [2,  2,  2,  1,       0.01872],
             [2,  2,  2,  2,     505.94646],
             [3,  2,  1,  1,      23.63452],
             [3,  2,  2,  2,      23.39780],
             [3,  3,  1,  1,      -3.02667],
             [3,  3,  2,  2,       2.70186],
             [3,  3,  3,  2,     -50.84416],
             [3,  3,  3,  3,     549.73449],
             [4,  2,  1,  1,     -19.19644],
             [4,  2,  2,  2,     -21.57758],
             [4,  3,  1,  1,      -3.55084],
             [4,  3,  2,  2,      -2.52690],
             [4,  3,  3,  2,       7.10251],
             [4,  3,  3,  3,     -77.07600],
             [4,  4,  1,  1,      -5.84458],
             [4,  4,  2,  2,      -2.86652],
             [4,  4,  3,  2,       5.88020],
             [4,  4,  3,  3,      19.57459],
             [4,  4,  4,  2,      29.29115],
             [4,  4,  4,  3,      73.27334],
             [4,  4,  4,  4,     383.60251],
             [5,  2,  1,  1,      24.88197],
             [5,  2,  2,  2,     -21.19389],
             [5,  3,  1,  1,      12.79430],
             [5,  3,  2,  2,       9.96353],
             [5,  3,  3,  2,       1.86175],
             [5,  3,  3,  3,      -3.77074],
             [5,  4,  1,  1,     -23.42343],
             [5,  4,  2,  2,     -19.40889],
             [5,  4,  3,  3,       2.51083],
             [5,  4,  4,  2,       0.59586],
             [5,  4,  4,  3,      -1.33782],
             [5,  4,  4,  4,      -4.17980],
             [5,  5,  1,  1,    -250.51096],
             [5,  5,  2,  2,    -202.73815],
             [5,  5,  3,  2,      -6.26073],
             [5,  5,  3,  3,      -0.71817],
             [5,  5,  4,  2,       4.41076],
             [5,  5,  4,  3,       0.74368],
             [5,  5,  4,  4,       0.70606],
             [5,  5,  5,  2,     -21.51808],
             [5,  5,  5,  3,      -5.28183],
             [5,  5,  5,  4,       6.20833],
             [5,  5,  5,  5,      26.72854],
             [6,  1,  1,  1,       2.33145],
             [6,  2,  2,  1,       6.58357],
             [6,  3,  3,  1,      -8.68835],
             [6,  4,  4,  1,       2.71077],
             [6,  5,  5,  1,      -0.37850],
             [6,  6,  1,  1,    -170.82753],
             [6,  6,  2,  2,    -167.87589],
             [6,  6,  3,  2,       4.53916],
             [6,  6,  3,  3,     -39.54289],
             [6,  6,  4,  2,      -1.53365],
             [6,  6,  4,  3,       4.63850],
             [6,  6,  4,  4,      -4.88103],
             [6,  6,  5,  2,       6.10598],
             [6,  6,  5,  3,      -0.78072],
             [6,  6,  5,  4,       1.41560],
             [6,  6,  5,  5,       7.68429],
             [6,  6,  6,  1,      -6.26849],
             [6,  6,  6,  6,      31.63034],
             [7,  2,  1,  1,       3.96811],
             [7,  2,  2,  2,       1.07101],
             [7,  3,  1,  1,     -29.28447],
             [7,  3,  2,  2,     -24.89279],
             [7,  3,  3,  2,      -6.33488],
             [7,  3,  3,  3,       7.24231],
             [7,  4,  1,  1,     -34.33637],
             [7,  4,  2,  2,     -31.16510],
             [7,  4,  3,  3,     -12.56004],
             [7,  4,  4,  2,       2.26610],
             [7,  4,  4,  3,      -0.96095],
             [7,  4,  4,  4,      -3.85367],
             [7,  5,  1,  1,      -6.15658],
             [7,  5,  2,  2,      -3.77331],
             [7,  5,  3,  3,       4.83778],
             [7,  5,  4,  4,       0.56545],
             [7,  5,  5,  2,      -0.49611],
             [7,  5,  5,  3,       3.81880],
             [7,  5,  5,  4,       3.66655],
             [7,  5,  5,  5,       0.68441],
             [7,  6,  6,  2,      -6.16266],
             [7,  6,  6,  3,       9.92742],
             [7,  6,  6,  4,       9.97438],
             [7,  6,  6,  5,      -6.53707],
             [7,  7,  1,  1,    -266.13407],
             [7,  7,  2,  2,    -252.95567],
             [7,  7,  3,  2,      -5.95357],
             [7,  7,  3,  3,     -23.82614],
             [7,  7,  4,  2,       3.22457],
             [7,  7,  4,  3,       2.88829],
             [7,  7,  4,  4,      -2.95230],
             [7,  7,  5,  2,     -20.08846],
             [7,  7,  5,  3,      -0.50168],
             [7,  7,  5,  4,       2.45986],
             [7,  7,  5,  5,      16.39147],
             [7,  7,  6,  1,      -7.80515],
             [7,  7,  6,  6,      58.98337],
             [7,  7,  7,  2,      -8.29049],
             [7,  7,  7,  3,      10.64617],
             [7,  7,  7,  4,      17.31625],
             [7,  7,  7,  5,       0.22639],
             [7,  7,  7,  7,     120.94216],
             [8,  1,  1,  1,      -0.19388],
             [8,  2,  2,  1,     -37.14979],
             [8,  3,  3,  1,      10.82796],
             [8,  4,  4,  1,       9.21324],
             [8,  6,  1,  1,     -21.26666],
             [8,  6,  2,  2,     -20.63124],
             [8,  6,  3,  3,      55.19087],
             [8,  6,  4,  4,     -15.05070],
             [8,  6,  5,  5,      -4.05642],
             [8,  6,  6,  1,      12.88020],
             [8,  6,  6,  6,       8.97881],
             [8,  7,  7,  1,      12.60408],
             [8,  7,  7,  6,      18.63533],
             [8,  8,  1,  1,    -158.16189],
             [8,  8,  2,  2,    -153.45228],
             [8,  8,  3,  2,       0.59511],
             [8,  8,  3,  3,     -82.61984],
             [8,  8,  4,  2,      -1.29526],
             [8,  8,  4,  3,       1.38179],
             [8,  8,  4,  4,     -59.10940],
             [8,  8,  5,  2,      33.63349],
             [8,  8,  5,  3,       0.54180],
             [8,  8,  5,  4,       8.20902],
             [8,  8,  5,  5,      51.76688],
             [8,  8,  6,  1,       0.65592],
             [8,  8,  6,  6,      10.13160],
             [8,  8,  7,  2,       8.74278],
             [8,  8,  7,  3,      -0.95384],
             [8,  8,  7,  4,       5.30499],
             [8,  8,  7,  5,      16.84072],
             [8,  8,  7,  7,      26.63937],
             [8,  8,  8,  1,      37.98053],
             [8,  8,  8,  6,       4.21083],
             [8,  8,  8,  8,      31.80347],
             [9,  2,  1,  1,      17.71928],
             [9,  2,  2,  2,       5.45841],
             [9,  3,  1,  1,       7.75018],
             [9,  3,  2,  2,       8.95414],
             [9,  3,  3,  2,      -6.71352],
             [9,  3,  3,  3,     -15.27981],
             [9,  4,  1,  1,      -0.70145],
             [9,  4,  2,  2,      -2.85353],
             [9,  4,  3,  3,     -18.34517],
             [9,  4,  4,  2,      -7.16397],
             [9,  4,  4,  3,      16.28174],
             [9,  4,  4,  4,      -8.86546],
             [9,  5,  1,  1,     -61.29902],
             [9,  5,  2,  2,     -54.31218],
             [9,  5,  3,  3,      13.88669],
             [9,  5,  4,  4,       0.81207],
             [9,  5,  5,  2,     -10.61006],
             [9,  5,  5,  3,      -0.91776],
             [9,  5,  5,  4,       0.53076],
             [9,  5,  5,  5,      13.38706],
             [9,  6,  6,  2,      -8.28108],
             [9,  6,  6,  3,      -2.67799],
             [9,  6,  6,  4,      -1.08179],
             [9,  6,  6,  5,       4.63593],
             [9,  7,  1,  1,      38.55319],
             [9,  7,  2,  2,      35.70565],
             [9,  7,  3,  3,     -61.72840],
             [9,  7,  4,  4,      21.12623],
             [9,  7,  5,  5,      -3.30041],
             [9,  7,  6,  6,     -13.65712],
             [9,  7,  7,  2,     -14.84869],
             [9,  7,  7,  3,     -11.46409],
             [9,  7,  7,  4,       1.19955],
             [9,  7,  7,  5,       7.60685],
             [9,  7,  7,  7,     -23.96468],
             [9,  8,  8,  2,       5.16779],
             [9,  8,  8,  3,       5.41853],
             [9,  8,  8,  4,      -1.26488],
             [9,  8,  8,  5,      21.91611],
             [9,  8,  8,  7,       5.10241],
             [9,  9,  1,  1,     -20.40716],
             [9,  9,  2,  2,     -22.29169],
             [9,  9,  3,  2,      11.19169],
             [9,  9,  3,  3,    -151.24397],
             [9,  9,  4,  2,      -6.63351],
             [9,  9,  4,  3,      -1.72724],
             [9,  9,  4,  4,     -77.25698],
             [9,  9,  5,  2,      -6.89194],
             [9,  9,  5,  3,       0.16269],
             [9,  9,  5,  4,       0.17722],
             [9,  9,  5,  5,       3.20777],
             [9,  9,  6,  1,       5.91328],
             [9,  9,  6,  6,      -0.99577],
             [9,  9,  7,  2,       5.96663],
             [9,  9,  7,  3,      -4.33510],
             [9,  9,  7,  4,       0.08697],
             [9,  9,  7,  5,      -0.71670],
             [9,  9,  7,  7,       2.64062],
             [9,  9,  8,  1,      -4.73727],
             [9,  9,  8,  6,      -5.01838],
             [9,  9,  8,  8,      11.21814],
             [9,  9,  9,  2,      -0.98925],
             [9,  9,  9,  3,      10.83010],
             [9,  9,  9,  4,       4.44835],
             [9,  9,  9,  5,       2.63948],
             [9,  9,  9,  7,      11.12104],
             [9,  9,  9,  9,      25.22733]
        ]

        legit = np.zeros((9, 9, 9, 9))
        mode_mapping = list(reversed(range(9)))
        for i, j, k, l, v in v4_Gaussian:
            i = mode_mapping[i - 1]
            j = mode_mapping[j - 1]
            k = mode_mapping[k - 1]
            l = mode_mapping[l - 1]
            for perm in ip.permutations((i, j, k, l)):
                legit[perm] = v

        v4 = self.h2w * ham.V_terms[2]

        print_errors = True
        if print_errors:
            if not np.allclose(legit, v4, atol=1):
                diff = legit - v4
                bad_pos = np.array(np.where(np.abs(diff) > 1)).T
                print("Gaussian/This Disagreements:\n" + "\n".join(
                    "{:>.0f} {:>.0f} {:>.0f} {:>.0f} {:>8.3f} (Actual: {:>8.3f} This: {:>8.3f})".format(
                        i, j, k, l, diff[i, j, k, l], legit[i, j, k, l], v4[i, j, k, l]
                    ) for i, j, k, l in bad_pos
                ))

        self.assertTrue(np.allclose(legit, v4, atol=1))  # testing to within a wavenumber

    @validationTest
    def test_OCHTCoriolisCouplings(self):
        # for unclear reasons this isn't working...?
        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("OCHT_freq.fchk"),
            internals=None
        )

        x, y, z = ham.coriolis_terms.get_zetas()

        x_els, y_els, z_els = [[[2, 1,     -0.40065],
                                [3, 1,      0.53649],
                                [4, 1,     -0.37833],
                                [5, 1,     -0.58017],
                                [6, 1,      0.26819]],
                               [[3, 2,      0.23121],
                                [4, 2,     -0.21924],
                                [4, 3,     -0.34391],
                                [5, 2,     -0.91672],
                                [5, 3,     -0.00152],
                                [5, 4,     -0.14017],
                                [6, 2,      0.06115],
                                [6, 3,     -0.90544],
                                [6, 4,      0.03576],
                                [6, 5,     -0.30877]],
                               [[2, 1,      0.72244],
                                [3, 1,      0.46646],
                                [4, 1,      0.16840],
                                [5, 1,     -0.33669],
                                [6, 1,     -0.34465]]]
        mapping = [2, 6, 5, 4, 3, 1]

        gaussian_x = np.zeros((6, 6))
        for j, i, v in x_els:
            i = mapping[i-1]-1
            j = mapping[j-1]-1
            gaussian_x[i, j] =  v
            gaussian_x[j, i] = -v
        gaussian_y = np.zeros((6, 6))
        for i, j, v in y_els:
            i = mapping[i-1]-1
            j = mapping[j-1]-1
            gaussian_y[i, j] =  v
            gaussian_y[j, i] = -v
        gaussian_z = np.zeros((6, 6))
        for i, j, v in z_els:
            i = mapping[i-1]-1
            j = mapping[j-1]-1
            gaussian_z[i, j] = v
            gaussian_z[j, i] = -v

        print_report = True
        if print_report:
            print("-"*50, "x", "-"*50)
            for a, b in zip(x, gaussian_z):
                print(
                    ("[" + "{:>8.3f}"*6 + "]").format(*a)
                    + ("[" + "{:>8.3f}"*6+ "]").format(*b)
                )
            print("-" * 50, "y", "-" * 50)
            for a, b in zip(y, gaussian_x):
                print(
                    ("[" + "{:>8.3f}" * 6 + "]").format(*a)
                    + ("[" + "{:>8.3f}" * 6 + "]").format(*b)
                )
            print("-" * 50, "z", "-" * 50)
            for a, b in zip(z, gaussian_y):
                print(
                    ("[" + "{:>8.3f}" * 6 + "]").format(*a)
                    + ("[" + "{:>8.3f}" * 6 + "]").format(*b)
                )

        sum_diff = np.abs(sum([x, y, z])) - np.abs(sum([gaussian_x, gaussian_y, gaussian_z]))
        # print(np.round(sum_diff, 3))
        self.assertLess(np.max(sum_diff), .001)

    @validationTest
    def test_HODCoriolisCouplings(self):

        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("HOD_freq.fchk"),
            internals=None
        )

        x, y, z = ham.coriolis_terms.get_zetas()

        x_els, y_els, z_els = [[

                                ],
                               [[2, 1,  0.01299],
                                [3, 1,  0.84340],
                                [3, 2, -0.53712]],
                               [

                               ]]
        mapping = [3, 2, 1] # remap modes to go from highest to lowest frequency

        gaussian_x = np.zeros((3, 3))
        for i, j, v in x_els:
            i = mapping[i - 1] - 1
            j = mapping[j - 1] - 1
            gaussian_x[i, j] = v
            gaussian_x[j, i] = -v
        gaussian_y = np.zeros((3, 3))
        for i, j, v in y_els:
            i = mapping[i - 1] - 1
            j = mapping[j - 1] - 1
            gaussian_y[i, j] = v
            gaussian_y[j, i] = -v
        gaussian_z = np.zeros((3, 3))
        for i, j, v in z_els:
            i = mapping[i - 1] - 1
            j = mapping[j - 1] - 1
            gaussian_z[i, j] = v
            gaussian_z[j, i] = -v

        print_report = False
        if print_report:
            print("-" * 50, "x", "-" * 50)
            for a, b in zip(x, gaussian_z):
                print(
                    ("[" + "{:>8.3f}" * 3 + "]").format(*a)
                    + ("[" + "{:>8.3f}" * 3 + "]").format(*b)
                )
            print("-" * 50, "y", "-" * 50)
            for a, b in zip(y, gaussian_x):
                print(
                    ("[" + "{:>8.3f}" * 3 + "]").format(*a)
                    + ("[" + "{:>8.3f}" * 3 + "]").format(*b)
                )
            print("-" * 50, "z", "-" * 50)
            for a, b in zip(z, gaussian_y):
                print(
                    ("[" + "{:>8.3f}" * 3 + "]").format(*a)
                    + ("[" + "{:>8.3f}" * 3 + "]").format(*b)
                )

        sum_diff = np.abs(sum([x, y, z])) - np.abs(sum([gaussian_x, gaussian_y, gaussian_z]))
        # print(np.round(sum_diff, 3))
        self.assertLess(np.max(sum_diff), .001)

    #endregion

    #region Test Algorithm

    @validationTest
    def test_SecondOrderEnergyCorrection(self):

        n_modes = 6
        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("HOD_freq.fchk"),
            internals=None  # doesn't matter for this
        )

        # states = self.get_states(3, 3)
        # h0, h1, h2 = ham.get_representations(states)

        # we supply hard-coded versions of H0, H1, and H2 to separate errors in those from errors in
        # the perturbation theory code
        h0 = np.array([
            [0.01846642088155693, 0.0, 0.0, 0.0, 4.13557989936697e-11, 2.0055884904224275e-11, -6.95971058561895e-13, -7.825322884799421e-13, -4.501742831265434e-12, -3.01100483487076e-12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.03611695804433793, -3.938830212179731e-12, -5.630969623303555e-12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.163034401735935e-11, 0.0, 0.0, -1.106667775363189e-12, 2.0055884904224275e-11, -6.366425766291432e-12, 0.0, -6.95971058561895e-13, 0.0, -3.01100483487076e-12],
            [0.0, -3.938830212179731e-12, 0.031269860697454216, 2.1654094707438123e-12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.47378133203069e-11, 0.0, 4.13557989936697e-11, -1.106667775363189e-12, 0.0, -4.258203873845191e-12, 0.0, -6.95971058561895e-13, -4.501742831265434e-12],
            [0.0, -5.630969623303555e-12, 2.1654094707438123e-12, 0.024945285665992512, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.2054576087328073e-12, 0.0, 0.0, 4.13557989936697e-11, 2.0055884904224275e-11, -6.366425766291432e-12, -4.258203873845191e-12, -7.825322884799421e-13],
            [4.13557989936697e-11, 0.0, 0.0, 0.0, 0.053767495207118925, 0.0, 0.0, -5.570347105949472e-12, -7.963393610586807e-12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [2.0055884904224275e-11, 0.0, 0.0, 0.0, 0.0, 0.044073300513351496, 0.0, -5.570347105949472e-12, 0.0, 3.0623514416170453e-12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [-6.95971058561895e-13, 0.0, 0.0, 0.0, 0.0, 0.0, 0.031424150450428096, 0.0, -7.963393610586807e-12, 3.0623514416170453e-12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [-7.825322884799421e-13, 0.0, 0.0, 0.0, -5.570347105949472e-12, -5.570347105949472e-12, 0.0, 0.04892039786023521, 2.1654094707438123e-12, -5.630969623303555e-12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [-4.501742831265434e-12, 0.0, 0.0, 0.0, -7.963393610586807e-12, 0.0, -7.963393610586807e-12, 2.1654094707438123e-12, 0.042595822828773514, -3.938830212179731e-12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [-3.01100483487076e-12, 0.0, 0.0, 0.0, 0.0, 3.0623514416170453e-12, 3.0623514416170453e-12, -5.630969623303555e-12, -3.938830212179731e-12, 0.03774872548188979, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 7.163034401735935e-11, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.07141803236989991, 0.0, 0.0, -6.8222540498825955e-12, 0.0, -9.753125483438739e-12, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 3.47378133203069e-11, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05687674032924876, 0.0, 0.0, -6.8222540498825955e-12, 0.0, 3.750599222519115e-12, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, -1.2054576087328073e-12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.03790301523486367, 0.0, 0.0, 0.0, 0.0, -9.753125483438739e-12, 3.750599222519115e-12, 0.0],
            [0.0, -1.106667775363189e-12, 4.13557989936697e-11, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -6.8222540498825955e-12, 0.0, 0.0, 0.06657093502301621, -7.877660424359464e-12, 2.1654094707438123e-12, 0.0, 0.0, 0.0, -7.963393610586807e-12],
            [0.0, 2.0055884904224275e-11, -1.106667775363189e-12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -6.8222540498825955e-12, 0.0, -7.877660424359464e-12, 0.06172383767613249, 0.0, -5.630969623303555e-12, 0.0, 0.0, 3.0623514416170453e-12],
            [0.0, -6.366425766291432e-12, 0.0, 4.13557989936697e-11, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -9.753125483438739e-12, 0.0, 0.0, 2.1654094707438123e-12, 0.0, 0.060246359991554504, 0.0, -1.1261939246607113e-11, 0.0, -5.570347105949472e-12],
            [0.0, 0.0, -4.258203873845191e-12, 2.0055884904224275e-11, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.750599222519115e-12, 0.0, 0.0, -5.630969623303555e-12, 0.0, 0.05055216529778707, 0.0, 4.330818941487625e-12, -5.570347105949472e-12],
            [0.0, -6.95971058561895e-13, 0.0, -6.366425766291432e-12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -9.753125483438739e-12, 0.0, 0.0, -1.1261939246607113e-11, 0.0, 0.04907468761320909, -3.938830212179731e-12, 3.0623514416170453e-12],
            [0.0, 0.0, -6.95971058561895e-13, -4.258203873845191e-12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.750599222519115e-12, 0.0, 0.0, 0.0, 4.330818941487625e-12, -3.938830212179731e-12, 0.044227590266325376, -7.963393610586807e-12],
            [0.0, -3.01100483487076e-12, -4.501742831265434e-12, -7.825322884799421e-13, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -7.963393610586807e-12, 3.0623514416170453e-12, -5.570347105949472e-12, -5.570347105949472e-12, 3.0623514416170453e-12, -7.963393610586807e-12, 0.055399262644670794]
        ])
        h1 = np.array([
            [0.0, -0.0016443729078422903, -0.0013877537983702287, -3.1781411870495364e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0016823474007175237, -0.0010448715868643585, 0.00013686781455426965, -0.00019641231752939666, 0.00010350166682222738, -0.00017691574963922763, -0.0001050919627883023, 0.00048491503933251114, 4.360275047362864e-05, 0.00013002287392060706],
            [-0.0016443729078422903, 0.0, 0.0, 0.0, -0.00523940564189364, 0.00010350166682222738, 0.00048491503933251114,-0.001665522761637432, -0.00028197806440769413, 0.00013002287392060706, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [-0.0013877537983702287, 0.0, 0.0, 0.0, -0.00019641231752939666, -0.003772350918724143, 4.360275047362864e-05, -0.0014979994468940748, 0.00013002287392060706, -0.00018040389094212117, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [-3.1781411870495364e-05, 0.0, 0.0, 0.0, -0.00017691574963922763, -0.0001050919627883023, 0.0001921163050302903, 0.00013002287392060706, -0.0009585994826195703, -0.0013260901972936532, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.00523940564189364, -0.00019641231752939666, -0.00017691574963922763, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.009985732955126163, 0.0, 0.0, -0.001943291724904636, 0.00014637346094821538, -0.0005321747169448928, 0.0, 0.0006857734252227202, 0.0, 0.0001838801117172495],
            [0.0, 0.00010350166682222738, -0.003772350918724143, -0.0001050919627883023, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.006836674794419553, 0.0, -0.0002777689632672036, -0.0013516259859458594, 0.0, -0.00032902637001374697, 0.0, 6.166360107657553e-05, 0.0001838801117172495],
            [0.0, 0.00048491503933251114, 4.360275047362864e-05, 0.0001921163050302903, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0005256339386890709, 0.0, 0.0, -0.00025019665253719876, -0.00014862247907162577, -0.00027282605739684974, -0.0012644265962170774, 0.0001838801117172495],
            [0.0, -0.001665522761637432, -0.0014979994468940748, 0.00013002287392060706, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.00034019611319326625, 0.0001792701456041638, 0.0, -0.005032402308249186, -0.004165175553782936, 0.0001838801117172495, 0.0001838801117172495, 4.360275047362864e-05, 0.00048491503933251114, -0.00043060054347931993],
            [0.0, -0.00028197806440769413, 0.00013002287392060706, -0.0009585994826195703, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.00030642706703427756, 0.0, 0.00083989748547817, 0.0001838801117172495, -0.0001050919627883023, -0.004269575563228618, 0.00010350166682222738, -0.00016171519424816507, 0.0001838801117172495, -0.0016038591605608565],
            [0.0, 0.00013002287392060706, -0.00018040389094212117, -0.0013260901972936532, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.00018202461901647743, 7.552217917007273e-05, -0.00017691574963922763, 0.0001838801117172495, -0.00019641231752939666, -0.0036851454177768855, 0.0001838801117172495, -1.8067620546314358e-05, -0.0008122260216713549],
            [-0.0016823474007175237, 0.0, 0.0, 0.0, -0.00998573295512616, 0.0, 0.0, -0.00034019611319326625, -0.00030642706703427756, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [-0.0010448715868643585, 0.0, 0.0, 0.0, 0.0, -0.00683667479441955, 0.0, 0.0001792701456041638, 0.0, -0.00018202461901647743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.00013686781455426965, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0005256339386890707, 0.0, 0.00083989748547817, 7.552217917007273e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [-0.00019641231752939666, 0.0, 0.0, 0.0, -0.001943291724904636, -0.0002777689632672036, 0.0, -0.005032402308249186, 0.0001838801117172495, -0.00017691574963922763, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.00010350166682222738, 0.0, 0.0, 0.0, 0.00014637346094821538, -0.0013516259859458594, 0.0, -0.004165175553782936, -0.0001050919627883023, 0.0001838801117172495, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [-0.00017691574963922763, 0.0, 0.0, 0.0, -0.0005321747169448928, 0.0, -0.00025019665253719876, 0.0001838801117172495, -0.004269575563228618, -0.00019641231752939666, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [-0.0001050919627883023, 0.0, 0.0, 0.0, 0.0, -0.00032902637001374697, -0.00014862247907162577, 0.0001838801117172495, 0.00010350166682222738, -0.0036851454177768855, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.00048491503933251114, 0.0, 0.0, 0.0, 0.0006857734252227202, 0.0, -0.00027282605739684974, 4.360275047362864e-05, -0.00016171519424816507, 0.0001838801117172495, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [4.360275047362864e-05, 0.0, 0.0, 0.0, 0.0, 6.166360107657553e-05, -0.0012644265962170774, 0.00048491503933251114, 0.0001838801117172495, -1.8067620546314358e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.00013002287392060706, 0.0, 0.0, 0.0, 0.0001838801117172495, 0.0001838801117172495, 0.0001838801117172495, -0.00043060054347931993, -0.0016038591605608565, -0.0008122260216713549, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        ])
        h2 = np.array([
            [0.00016893728625542488, 0.0, 0.0, 0.0, 0.00040424868443598213, 0.0002891097654542219, -0.00021553164706434878, 2.4818662079841186e-05, -1.4875021568903252e-05, 6.836862596819799e-06, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0007406312583562726, 2.4818662079841186e-05, -1.4875021568903252e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0014060757014638186, -2.3700926027540236e-05, -6.782156904121915e-05, 0.00011451632427152226, 0.0002963577824377973, 9.340316626644939e-05, -1.800567732468782e-05, -0.0006369307960913957, -3.267339178858443e-06, -1.4505080007608414e-05],
            [0.0, 2.4818662079841186e-05, 0.0005778002375752901, 6.836862596819799e-06, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.585167791337891e-05, 0.000874773833110212, -1.909869789865215e-05, 0.00041149670141955755, -5.952319552865031e-06, -1.5091032339285272e-05, 6.750851507516751e-05, -3.267339178858443e-06, -0.0002923835608554791, -4.033889464119042e-05],
            [0.0, -1.4875021568903252e-05, 6.836862596819799e-06, -0.00013587049214358845, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.607174743933356e-05, 3.339378441666386e-05, -0.0003345208870153343, -1.5091032339285272e-05, -1.800567732468782e-05, -1.7150464591064757e-05, 0.0002122578516630916, -0.00013850686067176722, -2.3411131310370607e-05, 2.0197946700226607e-05],
            [0.00040424868443598224, 0.0, 0.0, 0.0, 0.0021768682764619924, 5.125121959241341e-06, -0.00029797419586326533, 0.00011451632427152226, 9.340316626644939e-05, -1.5091032339285272e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.00028910976545422194, 0.0, 0.0, 0.0, 5.125121959241341e-06, 0.0014447435276446511, -5.434250938887216e-05, -5.952319552865031e-06, -1.800567732468782e-05, 6.750851507516751e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [-0.00021553164706434878, 0.0, 0.0, 0.0, -0.00029797419586326533, -5.434250938887216e-05, -0.0003931693436894778, -3.267339178858443e-06, -0.00013850686067176722, -2.3411131310370607e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [2.4818662079841186e-05, 0.0, 0.0, 0.0, 0.00011451632427152226, -5.952319552865031e-06, -3.267339178858443e-06, 0.0011597444535946205, -1.4505080007608414e-05, -4.033889464119042e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [-1.4875021568903252e-05, 0.0, 0.0, 0.0, 9.340316626644939e-05, -1.800567732468782e-05, -0.00013850686067176722, -1.4505080007608414e-05, -0.00016012491176927153, 2.0197946700226607e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [6.836862596819799e-06, 0.0, 0.0, 0.0, -1.5091032339285272e-05, 6.750851507516751e-05, -2.3411131310370607e-05, -4.033889464119042e-05, 2.0197946700226607e-05, 0.0001643074403985325, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0014060757014638184, 4.585167791337891e-05, 6.607174743933356e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.004477648340572581, 0.0, 0.0, 0.00023751937798615797, 8.876971628392951e-06, 0.0002545543908341627, 0.0, -0.0005161064465796555, 0.0, -2.6138434750307097e-05],
            [0.0, -2.3700926027540236e-05, 0.0008747738331102118, 3.339378441666386e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.002769767156463507, 0.0, 8.876971628392951e-06, -5.756732938867862e-05, 0.0, 0.00015351962184508721, 0.0, -9.412398727231534e-05, -3.118674795105016e-05],
            [0.0, -6.782156904121915e-05, -1.909869789865215e-05, -0.0003345208870153343, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0006029592683822425, 0.0, 0.0, -0.0005161064465796555, -9.412398727231534e-05, -0.00031350684139956923, -6.918711939376203e-05, -5.659197463343198e-06],
            [0.0, 0.00011451632427152226, 0.0004114967014195576, -1.5091032339285272e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00023751937798615808, 8.876971628392951e-06, 0.0, 0.0026062317156188225, 0.00010389536369897936, -3.5847022612036623e-05, -2.1341942604428205e-05, -4.620715379614582e-06, -0.00029797419586326533, 5.739181161707375e-05],
            [0.0, 0.00029635778243779736, -5.952319552865031e-06, -1.800567732468782e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.876971628392951e-06, -5.7567329388678686e-05, 0.0, 0.00010389536369897936, 0.0020369379875824645, -2.5463873072287226e-05, -6.580276771347766e-05, -5.434250938887216e-05, -4.620715379614582e-06, 3.732645039659695e-05],
            [0.0, 9.340316626644939e-05, -1.5091032339285272e-05, -1.7150464591064686e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0002545543908341628, 0.0, -0.0005161064465796555, -3.5847022612036623e-05, -2.5463873072287226e-05, 0.0006801637146099174, 5.125121959241341e-06, -3.403621320332538e-05, -2.1341942604428205e-05, 0.00010798164591380537],
            [0.0, -1.800567732468782e-05, 6.750851507516751e-05, 0.00021225785166309165, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00015351962184508727, -9.412398727231534e-05, -2.1341942604428205e-05, -6.580276771347766e-05, 5.125121959241341e-06, 0.0009225657116901491, -2.5463873072287226e-05, 4.868939299170736e-05, -1.2486997910581915e-05],
            [0.0, -0.0006369307960913957, -3.267339178858443e-06, -0.00013850686067176722, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0005161064465796555, 0.0, -0.00031350684139956945, -4.620715379614582e-06, -5.434250938887216e-05, -3.403621320332538e-05, -2.5463873072287226e-05, -0.0010133721550416915, 1.5577231320612022e-05, -5.359319598894117e-05],
            [0.0, -3.267339178858443e-06, -0.0002923835608554791, -2.3411131310370607e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -9.412398727231534e-05, -6.918711939376206e-05, -0.00029797419586326533, -4.620715379614582e-06, -2.1341942604428205e-05, 4.868939299170736e-05, 1.5577231320612022e-05, -0.00020167642992510108, -0.0001745182153211429],
            [0.0, -1.4505080007608414e-05, -4.033889464119042e-05, 2.0197946700226607e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.6138434750307097e-05, -3.118674795105016e-05, -5.659197463343198e-06, 5.739181161707375e-05, 3.732645039659695e-05, 0.00010798164591380537, -1.2486997910581915e-05, -5.359319598894117e-05, -0.0001745182153211429, 0.0001503032646913323]
        ])

        class FakePert:
            def __init__(self, h):
                self.h=h
            def __getitem__(self, item):
                return self.h
        corrs = ham._get_VPT_corrections(
            [FakePert(h) for h in [h0, h1, h2]],
            [0],
            np.arange(len(h0)),
            len(h0),
            2
        )

        e_corr = (corrs.energy_corrs[0] * self.h2w)

        self.assertTrue(np.allclose(
            e_corr,
            [4052.91093, 0., -50.05456], # from Gaussian
            .01 # within .01 wavenumbers
        ))

    @inactiveTest
    def test_HOHNielsenEnergies(self):
        # need to figure out the units on this shit if I want it to work...
        states = (
            (0, 0, 0),
            (0, 0, 1), (0, 1, 0), (1, 0, 0),
            (0, 0, 2), (0, 2, 0), (2, 0, 0),
            (0, 1, 1), (1, 0, 1), (1, 1, 0)
        )

        hammer = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("HOD_freq.fchk"),
            internals=None
        )

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        e_harm, e_corr = hammer.get_Nielsen_energies(states)

        # wfns = hammer.get_wavefunctions(states, coupled_states=None)
        energies = h2w * (e_harm + e_corr)
        zero_ord = h2w * e_harm

        # print(wfns.corrs.coupled_states)
        # print([
        #     np.max(np.abs(wfns.corrs.wfn_corrections[i, 1]))
        #     for i in range(len(wfns.corrs.wfn_corrections))
        # ])
        # print([
        #     np.max(np.abs(wfns.corrs.wfn_corrections[i, 2]))
        #     for i in range(len(wfns.corrs.wfn_corrections))
        # ])

        gaussian_states = [(0, 0, 1), (0, 1, 0), (1, 0, 0),
                           (0, 0, 2), (0, 2, 0), (2, 0, 0),
                           (0, 1, 1), (1, 0, 1), (1, 1, 0)]
        gaussian_energies = np.array([4681.564, 4605.953])
        gaussian_freqs = np.array([
            [3937.525, 3744.734],
            [3803.300, 3621.994],
            [1622.303, 1572.707],

            [7875.049, 7391.391],
            [7606.599, 7155.881],
            [3244.606, 3117.366],

            [7740.824, 7200.364],
            [5559.828, 5294.379],
            [5425.603, 5174.665]
        ])

        my_energies = np.array([zero_ord[0], energies[0]])
        my_freqs = np.column_stack([
            zero_ord[1:] - zero_ord[0],
            energies[1:] - energies[0]
        ])

        print_report = False
        if print_report:
            print("Gaussian:\n",
                  "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*gaussian_energies, "-", "-"),
                  *(
                      "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                      zip(gaussian_states, gaussian_freqs)
                  )
                  )
            print("State Energies:\n",
                  "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*my_energies, "-", "-"),
                  *(
                      "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                      zip(states[1:], my_freqs)
                  )
                  )

        print_diffs = True
        if print_diffs:
            print("Difference Energies:\n",
                  "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*(my_energies - gaussian_energies), "-", "-"),
                  *(
                      "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                      zip(states[1:], my_freqs - gaussian_freqs[:len(my_freqs)])
                  )
                  )
        self.assertLess(np.max(np.abs(my_freqs - gaussian_freqs[:len(my_freqs)])), 1.5)

    #endregion

    #region Test Systems

    @validationTest
    def test_DODVPTCartesians(self):

        internals = None
        # [
        #     [0, -1, -1, -1],
        #     [1, 0, -1, -1],
        #     [2, 0, 1, -1]
        # ]
        states = (
            (0, 0, 0),
            (0, 0, 1), (0, 1, 0), (1, 0, 0),
            (0, 0, 2), (0, 2, 0), (2, 0, 0),
            (0, 1, 1), (1, 0, 1), (1, 1, 0)
        )

        with BlockProfiler("DOD Cartesians",
                           print_res=False,
                           strip_dirs=self.job_is_dumb,
                           sort_by='ncalls'
                           ):
            wfns = self.get_VPT2_wfns(
                "DOD_freq.fchk",
                internals,
                states,
                save_coeffs=True,
                regenerate=True
                # , coupled_states=self.get_states(5, 3)
                , log=True
                # , watson=0.
                # , v3=0
                # , v4=0
            )

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        # wfns = hammer.get_wavefunctions(states, coupled_states=None)
        energies = h2w * wfns.energies
        zero_ord = h2w * wfns.zero_order_energies

        # print(wfns.corrs.coupled_states)
        # print([
        #     np.max(np.abs(wfns.corrs.wfn_corrections[i, 1]))
        #     for i in range(len(wfns.corrs.wfn_corrections))
        # ])
        # print([
        #     np.max(np.abs(wfns.corrs.wfn_corrections[i, 2]))
        #     for i in range(len(wfns.corrs.wfn_corrections))
        # ])

        gaussian_states = [(0, 0, 1), (0, 1, 0), (1, 0, 0),
                           (0, 0, 2), (0, 2, 0), (2, 0, 0),
                           (0, 1, 1), (1, 0, 1), (1, 1, 0)]
        gaussian_energies = np.array([3406.854, 3366.588])
        gaussian_freqs = np.array([
            [2884.314, 2779.931],
            [2742.310, 2648.051],
            [1187.085, 1160.758],

            [5768.628, 5504.706],
            [5484.619, 5250.472],
            [2374.171, 2306.554],

            [5626.624, 5341.492],
            [4071.399, 3928.727],
            [3929.395, 3798.041]
        ])

        my_energies = np.array([zero_ord[0], energies[0]])
        my_freqs = np.column_stack([
            zero_ord[1:] - zero_ord[0],
            energies[1:] - energies[0]
        ])

        print_report = False
        if print_report:
            print("Gaussian:\n",
                  "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*gaussian_energies, "-", "-"),
                  *(
                      "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                      zip(gaussian_states, gaussian_freqs)
                  )
                  )
            print("State Energies:\n",
                  "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*my_energies, "-", "-"),
                  *(
                      "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                      zip(states[1:], my_freqs)
                  )
                  )

        print_diffs = True
        if print_diffs:
            print("Difference Energies:\n",
                  "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*(my_energies - gaussian_energies), "-", "-"),
                  *(
                      "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                      zip(states[1:], my_freqs - gaussian_freqs[:len(my_freqs)])
                  )
                  )
        self.assertLess(np.max(np.abs(my_freqs - gaussian_freqs[:len(my_freqs)])), 1.5)

    @validationTest
    def test_DODVPTInternals(self):

        internals = [
            [0, -1, -1, -1],
            [1,  0, -1, -1],
            [2,  0,  1, -1]
        ]
        states = (
            (0, 0, 0),
            (0, 0, 1), (0, 1, 0), (1, 0, 0),
            (0, 0, 2), (0, 2, 0), (2, 0, 0),
            (0, 1, 1), (1, 0, 1), (1, 1, 0)
        )

        wfns = self.get_VPT2_wfns(
            "DOD_freq.fchk",
            internals,
            states,
            save_coeffs=True,
            regenerate=True
            , coupled_states=self.get_states(5, 3)
            # , watson=0.
            # , v3=0
            # , v4=0
        )

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        # wfns = hammer.get_wavefunctions(states, coupled_states=None)
        energies = h2w * wfns.energies
        zero_ord = h2w * wfns.zero_order_energies

        gaussian_states = [(0, 0, 1), (0, 1, 0), (1, 0, 0),
                           (0, 0, 2), (0, 2, 0), (2, 0, 0),
                           (0, 1, 1), (1, 0, 1), (1, 1, 0)]
        gaussian_energies = np.array([3406.854, 3366.588])
        gaussian_freqs = np.array([
            [2884.314, 2779.931],
            [2742.310, 2648.051],
            [1187.085, 1160.758],

            [5768.628, 5504.706],
            [5484.619, 5250.472],
            [2374.171, 2306.554],

            [5626.624, 5341.492],
            [4071.399, 3928.727],
            [3929.395, 3798.041]
        ])

        my_energies = np.array([zero_ord[0], energies[0]])
        my_freqs = np.column_stack([zero_ord[1:] - zero_ord[0], energies[1:] - energies[0]])

        print_report = False
        if print_report:
            print("Gaussian:\n",
                  "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*gaussian_energies, "-", "-"),
                  *(
                      "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                      zip(gaussian_states, gaussian_freqs)
                  )
                  )
            print("State Energies:\n",
                  "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*my_energies, "-", "-"),
                  *(
                      "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                      zip(states[1:], my_freqs)
                  )
                  )

        print_diffs = True
        if print_diffs:
            print("Difference Energies:\n",
                  "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*(my_energies - gaussian_energies), "-", "-"),
                  *(
                      "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                      zip(states[1:], my_freqs - gaussian_freqs[:len(my_freqs)])
                  )
                  )
        self.assertLess(np.max(np.abs(my_freqs - gaussian_freqs[:len(my_freqs)])), 1.5)

    @validationTest
    def test_HOHVPTInternals(self):

        internals = [
            [0, -1, -1, -1],
            [1,  0, -1, -1],
            [2,  0,  1, -1]
        ]
        states = (
            (0, 0, 0),
            (0, 0, 1), (0, 1, 0), (1, 0, 0),
            (0, 0, 2), (0, 2, 0), (2, 0, 0),
            (0, 1, 1), (1, 0, 1), (1, 1, 0)
        )

        wfns = self.get_VPT2_wfns(
            "HOH_freq.fchk",
            internals,
            states,
            save_coeffs=True,
            regenerate=True
            , coupled_states = self.get_states(5, 3)
            # , watson=0.
            # , v3=0
            # , v4=0
        )

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        # wfns = hammer.get_wavefunctions(states, coupled_states=None)
        energies = h2w * wfns.energies
        zero_ord = h2w * wfns.zero_order_energies

        # print(wfns.corrs.coupled_states)
        # print([
        #     np.max(np.abs(wfns.corrs.wfn_corrections[i, 1]))
        #     for i in range(len(wfns.corrs.wfn_corrections))
        # ])
        # print([
        #     np.max(np.abs(wfns.corrs.wfn_corrections[i, 2]))
        #     for i in range(len(wfns.corrs.wfn_corrections))
        # ])

        gaussian_states = [(0, 0, 1), (0, 1, 0), (1, 0, 0),
                           (0, 0, 2), (0, 2, 0), (2, 0, 0),
                           (0, 1, 1), (1, 0, 1), (1, 1, 0)]
        gaussian_energies = np.array([4681.564, 4605.953])
        gaussian_freqs = np.array([
            [3937.525, 3744.734],
            [3803.300, 3621.994],
            [1622.303, 1572.707],

            [7875.049, 7391.391],
            [7606.599, 7155.881],
            [3244.606, 3117.366],

            [7740.824, 7200.364],
            [5559.828, 5294.379],
            [5425.603, 5174.665]
        ])

        my_energies = np.array([zero_ord[0], energies[0]])
        my_freqs = np.column_stack([
            zero_ord[1:] - zero_ord[0],
            energies[1:] - energies[0]
        ])

        print_report = False
        if print_report:
            print("Gaussian:\n",
                  "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*gaussian_energies, "-", "-"),
                  *(
                      "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                      zip(gaussian_states, gaussian_freqs)
                  )
                  )
            print("State Energies:\n",
                  "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*my_energies, "-", "-"),
                  *(
                      "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                      zip(states[1:], my_freqs)
                  )
                  )

        print_diffs = True
        if print_diffs:
            print("Difference Energies:\n",
                  "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*(my_energies - gaussian_energies), "-", "-"),
                  *(
                      "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                      zip(states[1:], my_freqs - gaussian_freqs[:len(my_freqs)])
                  )
                  )
        self.assertLess(np.max(np.abs(my_freqs - gaussian_freqs[:len(my_freqs)])), 1.5)

    @validationTest
    def test_HOHVPTCartesians(self):

        internals = None
        # [
        #     [0, -1, -1, -1],
        #     [1, 0, -1, -1],
        #     [2, 0, 1, -1]
        # ]
        states = (
            (0, 0, 0),
            (0, 0, 1), (0, 1, 0), (1, 0, 0),
            (0, 0, 2), (0, 2, 0), (2, 0, 0),
            (0, 1, 1), (1, 0, 1), (1, 1, 0)
        )

        wfns = self.get_VPT2_wfns(
            "HOH_freq.fchk",
            internals,
            states,
            regenerate=True
            # , coupled_states = self.get_states(5, 3)
            # , v3=0
            # , v4=0
        )

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        # wfns = hammer.get_wavefunctions(states, coupled_states=None)
        energies = h2w * wfns.energies
        zero_ord = h2w * wfns.zero_order_energies

        # print(len(self.get_states(5, 3)), len(wfns.corrs.coupled_states))
        # print([
        #     np.max(np.abs(wfns.corrs.wfn_corrections[i, 1]))
        #     for i in range(len(wfns.corrs.wfn_corrections))
        # ])
        # print([
        #     np.max(np.abs(wfns.corrs.wfn_corrections[i, 2]))
        #     for i in range(len(wfns.corrs.wfn_corrections))
        # ])

        gaussian_states = [(0, 0, 1), (0, 1, 0), (1, 0, 0),
                           (0, 0, 2), (0, 2, 0), (2, 0, 0),
                           (0, 1, 1), (1, 0, 1), (1, 1, 0)]
        gaussian_energies = np.array([4681.564, 4605.953])
        gaussian_freqs = np.array([
            [3937.525, 3744.734],
            [3803.300, 3621.994],
            [1622.303, 1572.707],

            [7875.049, 7391.391],
            [7606.599, 7155.881],
            [3244.606, 3117.366],

            [7740.824, 7200.364],
            [5559.828, 5294.379],
            [5425.603, 5174.665]
        ])

        my_energies = np.array([zero_ord[0], energies[0]])
        my_freqs = np.column_stack([
            zero_ord[1:] - zero_ord[0],
            energies[1:] - energies[0]
        ])

        print_report = False
        if print_report:
            print("Gaussian:\n",
                  "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*gaussian_energies, "-", "-"),
                  *(
                      "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                      zip(gaussian_states, gaussian_freqs)
                  )
                  )
            print("State Energies:\n",
                  "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*my_energies, "-", "-"),
                  *(
                      "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                      zip(states[1:], my_freqs)
                  )
                  )

        print_diffs = True
        if print_diffs:
            print("Difference Energies:\n",
                  "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*(my_energies - gaussian_energies), "-", "-"),
                  *(
                      "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                      zip(states[1:], my_freqs - gaussian_freqs[:len(my_freqs)])
                  )
                  )
        self.assertLess(np.max(np.abs(my_freqs - gaussian_freqs[:len(my_freqs)])), 1.5)

    @validationTest
    def test_HODVPTInternals(self):

        internals = [
            [0, -1, -1, -1],
            [1,  0, -1, -1],
            [2,  0,  1, -1]
        ]
        states = (
            (0, 0, 0),
            (0, 0, 1), (0, 1, 0), (1, 0, 0),
            (0, 0, 2), (0, 2, 0), (2, 0, 0),
            (0, 1, 1), (1, 0, 1), (1, 1, 0)
        )

        coupled_states = self.get_states(5, 3)
        wfns = self.get_VPT2_wfns(
            "HOD_freq.fchk",
            internals,
            states,
            save_coeffs=True,
            regenerate=True,
            coupled_states=coupled_states
            # , v3=legit/self.h2w
            # , v4=0
        )

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        # wfns = hammer.get_wavefunctions(states, coupled_states=None)
        energies = h2w * wfns.energies
        zero_ord = h2w * wfns.zero_order_energies
        # raise Exception([energies, zero_ord])

        # print([
        #     np.max(np.abs(wfns.corrs.wfn_corrections[i, 1]))
        #     for i in range(len(wfns.corrs.wfn_corrections))
        # ])
        # print([
        #     np.max(np.abs(wfns.corrs.wfn_corrections[i, 2]))
        #     for i in range(len(wfns.corrs.wfn_corrections))
        # ])

        gaussian_states = [(0, 0, 1), (0, 1, 0), (1, 0, 0),
                           (0, 0, 2), (0, 2, 0), (2, 0, 0),
                           (0, 1, 1), (1, 0, 1), (1, 1, 0)]
        gaussian_energies = np.array([4052.912, 3994.844])
        gaussian_freqs = np.array([
            [3873.846, 3685.815],
            [2810.031, 2706.132],
            [1421.946, 1383.391],
            [7747.692, 7202.835],
            [5620.062, 5323.917],
            [2843.893, 2749.027],
            [6683.877, 6377.958],
            [5295.792, 5044.721],
            [4231.977, 4072.407]
        ])

        my_energies = np.array([zero_ord[0], energies[0]])
        my_freqs = np.column_stack([
            zero_ord[1:] - zero_ord[0],
            energies[1:] - energies[0]
        ])

        print_report = False
        if print_report:
            print("Gaussian:\n",
                  "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*gaussian_energies, "-", "-"),
                  *(
                      "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                      zip(gaussian_states, gaussian_freqs)
                  )
                  )
            print("State Energies:\n",
                  "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*my_energies, "-", "-"),
                  *(
                      "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                      zip(states[1:], my_freqs)
                  )
                  )

        print_diffs = True
        if print_diffs:
            print("Difference Energies:\n",
                  "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*(my_energies-gaussian_energies), "-", "-"),
                  *(
                      "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                      zip(states[1:], my_freqs - gaussian_freqs[:len(my_freqs)])
                  )
                  )
        self.assertLess(np.max(np.abs(my_freqs - gaussian_freqs[:len(my_freqs)])), 1)

    @validationTest
    def test_HODVPTCartesians(self):

        internals = None

        states = (
            (0, 0, 0),
            (0, 0, 1), (0, 1, 0), (1, 0, 0),
            (0, 0, 2), (0, 2, 0), (2, 0, 0),
            (0, 1, 1), (1, 0, 1), (1, 1, 0)
        )

        coupled_states = self.get_states(5, 3)
        wfns = self.get_VPT2_wfns(
            "HOD_freq.fchk",
            internals,
            states,
            regenerate=True,
            coupled_states=coupled_states
            # , v3=0
            # , v4=0
        )

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        # wfns = hammer.get_wavefunctions(states, coupled_states=None)
        energies = h2w * wfns.energies
        zero_ord = h2w * wfns.zero_order_energies
        # raise Exception([energies, zero_ord])

        print([
            np.max(np.abs(wfns.corrs.wfn_corrections[i, 1]))
            for i in range(len(wfns.corrs.wfn_corrections))
        ])
        # print([
        #     np.max(np.abs(wfns.corrs.wfn_corrections[i, 2]))
        #     for i in range(len(wfns.corrs.wfn_corrections))
        # ])

        gaussian_states = [(0, 0, 1), (0, 1, 0), (1, 0, 0),
                           (0, 0, 2), (0, 2, 0), (2, 0, 0),
                           (0, 1, 1), (1, 0, 1), (1, 1, 0)]
        gaussian_energies = np.array([4052.912, 3994.844])
        gaussian_freqs = np.array([
            [3873.846, 3685.815],
            [2810.031, 2706.132],
            [1421.946, 1383.391],
            [7747.692, 7202.835],
            [5620.062, 5323.917],
            [2843.893, 2749.027],
            [6683.877, 6377.958],
            [5295.792, 5044.721],
            [4231.977, 4072.407]
        ])

        my_energies = np.array([zero_ord[0], energies[0]])
        my_freqs = np.column_stack([
            zero_ord[1:] - zero_ord[0],
            energies[1:] - energies[0]
        ])

        print_report = False
        if print_report:
            print("Gaussian:\n",
                  "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*gaussian_energies, "-", "-"),
                  *(
                      "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                      zip(gaussian_states, gaussian_freqs)
                  )
                  )
            print("State Energies:\n",
                  "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*my_energies, "-", "-"),
                  *(
                      "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                      zip(states[1:], my_freqs)
                  )
                  )

        print_diffs = True
        if print_diffs:
            print("Difference Energies:\n",
                  "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*(my_energies - gaussian_energies), "-", "-"),
                  *(
                      "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                      zip(states[1:], my_freqs - gaussian_freqs[:len(my_freqs)])
                  )
                  )
        self.assertLess(np.max(np.abs(my_freqs - gaussian_freqs[:len(my_freqs)])), 1)

    @validationTest
    def test_HOTVPTInternals(self):

        internals = [
            [0, -1, -1, -1],
            [1,  0, -1, -1],
            [2,  0,  1, -1]
        ]
        states = (
            (0, 0, 0),
            (0, 0, 1), (0, 1, 0), (1, 0, 0),
            (0, 0, 2), (0, 2, 0), (2, 0, 0),
            (0, 1, 1), (1, 0, 1), (1, 1, 0)
        )
        n_modes = 3
        coupled_states = self.get_states(5, n_modes, max_quanta=5)
        def block(self=self, internals=internals, states=states, coupled_states=coupled_states):
            return self.get_VPT2_wfns(
                "HOT_freq.fchk",
                internals,
                states,
                regenerate=True,
                coupled_states=coupled_states
            )
        wfns = block()
        exc, stat_block, wfns = self.profile_block(block)

        do_profile = False
        if exc is not None:
            try:
                raise exc
            except:
                raise Exception(stat_block)
        elif do_profile:
            raise Exception(stat_block)

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")
        engs = h2w * wfns.energies
        # wfns = hammer.get_wavefunctions(states, coupled_states=None)
        energies = h2w * wfns.energies
        zero_ord = h2w * wfns.zero_order_energies
        # raise Exception([energies, zero_ord])

        gaussian_states = [(0, 0, 1), (0, 1, 0), (1, 0, 0),
                           (0, 0, 2), (0, 2, 0), (2, 0, 0),
                           (0, 1, 1), (1, 0, 1), (1, 1, 0)]
        gaussian_energies = np.array([3789.117, 3736.855])
        gaussian_freqs = np.array([
            [3872.325, 3692.047],
            [2357.491, 2285.849],
            [1348.418, 1313.165],
            [7744.651, 7214.810],
            [4714.982, 4509.256],
            [2696.835, 2607.130],
            [6229.816, 5973.755],
            [5220.743, 4987.363],
            [3705.909, 3584.756]
        ])

        my_energies = np.array([zero_ord[0], energies[0]])
        my_freqs = np.column_stack([
            zero_ord[1:] - zero_ord[0],
            energies[1:] - energies[0]
        ])

        print_report = False
        if print_report:
            print("Gaussian:\n",
                  "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*gaussian_energies, "-", "-"),
                  *(
                      "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                      zip(gaussian_states, gaussian_freqs)
                  )
                  )
            print("State Energies:\n",
                  "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*my_energies, "-", "-"),
                  *(
                      "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                      zip(states[1:], my_freqs)
                  )
                  )
        print_diffs = True
        if print_diffs:
            print("Difference Energies:\n",
                  "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*(my_energies - gaussian_energies), "-", "-"),
                  *(
                      "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                      zip(states[1:], my_freqs - gaussian_freqs[:len(my_freqs)])
                  )
                  )

        self.assertLess(np.max(np.abs(my_freqs - gaussian_freqs[:len(my_freqs)])), 1)

    @validationTest
    def test_HOTVPTCartesians(self):

        internals = None
        states = (
            (0, 0, 0),
            (0, 0, 1), (0, 1, 0), (1, 0, 0),
            (0, 0, 2), (0, 2, 0), (2, 0, 0),
            (0, 1, 1), (1, 0, 1), (1, 1, 0)
        )
        n_modes = 3
        coupled_states = self.get_states(5, n_modes, max_quanta=5)

        def block(self=self, internals=internals, states=states, coupled_states=coupled_states):
            return self.get_VPT2_wfns(
                "HOT_freq.fchk",
                internals,
                states,
                regenerate=True,
                coupled_states=coupled_states
                # , coriolis=0
            )

        # wfns = block()
        exc, stat_block, wfns = self.profile_block(block)

        do_profile = False
        if exc is not None:
            try:
                raise exc
            except:
                raise Exception(stat_block)
        elif do_profile:
            raise Exception(stat_block)

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")
        engs = h2w * wfns.energies
        # wfns = hammer.get_wavefunctions(states, coupled_states=None)
        energies = h2w * wfns.energies
        zero_ord = h2w * wfns.zero_order_energies
        # raise Exception([energies, zero_ord])

        gaussian_states = [(0, 0, 1), (0, 1, 0), (1, 0, 0),
                           (0, 0, 2), (0, 2, 0), (2, 0, 0),
                           (0, 1, 1), (1, 0, 1), (1, 1, 0)]
        gaussian_energies = np.array([3789.117, 3736.855])
        gaussian_freqs = np.array([
            [3872.325, 3692.047],
            [2357.491, 2285.849],
            [1348.418, 1313.165],
            [7744.651, 7214.810],
            [4714.982, 4509.256],
            [2696.835, 2607.130],
            [6229.816, 5973.755],
            [5220.743, 4987.363],
            [3705.909, 3584.756]
        ])

        my_energies = np.array([zero_ord[0], energies[0]])
        my_freqs = np.column_stack([
            zero_ord[1:] - zero_ord[0],
            energies[1:] - energies[0]
        ])

        print_report = False
        if print_report:
            print("Gaussian:\n",
                  "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*gaussian_energies, "-", "-"),
                  *(
                      "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                      zip(gaussian_states, gaussian_freqs)
                  )
                  )
            print("State Energies:\n",
                  "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*my_energies, "-", "-"),
                  *(
                      "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                      zip(states[1:], my_freqs)
                  )
                  )
        print_difference = True
        if print_difference:
            print("Difference Energies:\n",
                  "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*(my_energies - gaussian_energies), "-", "-"),
                  *(
                      "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                      zip(states[1:], my_freqs - gaussian_freqs[:len(my_freqs)])
                  )
                  )
        self.assertLess(np.max(np.abs(my_freqs - gaussian_freqs[:len(my_freqs)])), 1)

    @inactiveTest
    def test_OCHHVPTInternals(self):

        internals = [
            [0, -1, -1, -1],
            [1,  0, -1, -1],
            [2,  1,  0, -1],
            [3,  1,  0,  2]
        ]

        n_modes = 3 * 4 - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)

        states = self.get_states(2, n_modes)
        coupled_states = self.get_states(5, n_modes, max_quanta=5)

        def block(self=self,
                  internals=internals,
                  states=states,
                  coupled_states=coupled_states,
                  mode_selection=mode_selection
                  ):
            return self.get_VPT2_wfns(
                "OCHH_freq.fchk",
                internals,
                states,
                regenerate=True,
                coupled_states=coupled_states,
                mode_selection=mode_selection
            )

        # block()
        exc, stat_block, wfns = self.profile_block(block)

        do_profile = False
        if do_profile:
            if exc is not None:
                try:
                    raise exc
                except:
                    raise Exception(stat_block)
            else:
                raise Exception(stat_block)
        elif exc is not None:
            raise exc

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")
        engs = h2w * wfns.energies
        freqs = engs[1:] - engs[0]
        harm_engs = h2w * wfns.zero_order_energies
        harm_freq = harm_engs[1:] - harm_engs[0]

        gaussian_engs = [5866.872, 5785.953]
        gaussian_freqs = np.array([
            [3061.701,    2849.454],
            [2977.640,    2820.563],
            [1727.083,    1695.163],
            [1527.041,    1493.139],
            [1252.164,    1231.359],
            [1188.114,    1166.697],
            # 2 quanta
            [6123.403,    5719.555],
            [5955.281,    5578.213],
            [3454.165,    3372.197],
            [3054.082,    2988.763],
            [2504.328,    2459.646],
            [2376.227,    2327.621],
            # mixed states
            [6039.342,    5586.305],
            [4788.784,    4589.165],
            [4704.723,    4510.744],
            [4588.742,    4374.903],
            [4504.681,    4278.834],
            [3254.123,    3181.610],
            [4313.865,    4114.701],
            [4229.804,    4043.264],
            [2979.247,    2967.726],
            [2779.205,    2710.416],
            [4249.815,    4043.462],
            [4165.754,    3978.440],
            [2915.196,    2854.766],
            [2715.155,    2657.682],
            [2440.278,    2404.805]
        ])

        print_report = False
        if print_report:
            if n_modes == 6:
                print("Gaussian Energies:\n",
                      ('0 ' * n_modes + "{:>8.3f} {:>8.3f} {:>8} {:>8}\n").format(*gaussian_engs, "-", "-"),
                      *(
                          ('{:<1.0f} ' * n_modes + "{:>8} {:>8} {:>8.3f} {:>8.3f}\n").format(*s, "-", "-", *e) for s, e
                      in
                          zip(states[1:], gaussian_freqs)
                      )
                      )
            print("State Energies:\n",
                  ('0 ' * n_modes + "{:>8.3f} {:>8.3f} {:>8} {:>8}\n").format(harm_engs[0], engs[0], "-", "-"),
                  *(
                      ('{:<1.0f} ' * n_modes + "{:>8} {:>8} {:>8.3f} {:>8.3f}\n").format(*s, "-", "-", e1, e2) for
                  s, e1, e2 in
                      zip(states[1:], harm_freq, freqs)
                  )
                  )
        ns = len(states) - 1
        print_difference = True
        if print_difference:
            print("Difference Energies:\n",
                  ('0 ' * n_modes + "{:>8.3f} {:>8.3f} {:>8} {:>8}\n").format(harm_engs[0] - gaussian_engs[0],
                                                                              engs[0] - gaussian_engs[1], "-", "-"),
                  *(
                      ('{:<1.0f} ' * n_modes + "{:>8} {:>8} {:>8.3f} {:>8.3f}\n").format(*s, "-", "-", e1, e2) for
                      s, e1, e2 in
                      zip(states[1:], harm_freq[:ns] - gaussian_freqs[:ns, 0], freqs[:ns] - gaussian_freqs[:ns, 1])
                  )
                  )

        self.assertLess(
            np.max(np.abs(freqs[:ns] - gaussian_freqs[:ns, 1])),
            1)

    @inactiveTest
    def test_OCHHVPTCartesians(self):

        internals = None

        n_modes = 3 * 4 - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)

        states = self.get_states(2, n_modes)
        coupled_states = self.get_states(5, n_modes, max_quanta=5)

        def block(self=self,
                  internals=internals,
                  states=states,
                  coupled_states=coupled_states,
                  mode_selection=mode_selection
                  ):
            return self.get_VPT2_wfns(
                "OCHH_freq.fchk",
                internals,
                states,
                regenerate=True,
                coupled_states=coupled_states,
                mode_selection=mode_selection
                # , degeneracies=(
                #     coupled_states.index((0, 0, 0, 0, 0, 1)),
                #     coupled_states.index((0, 0, 2, 0, 0, 0))
                # )
            )

        # block()
        exc, stat_block, wfns = self.profile_block(block)

        do_profile = False
        if do_profile:
            if exc is not None:
                try:
                    raise exc
                except:
                    raise Exception(stat_block)
            else:
                raise Exception(stat_block)
        elif exc is not None:
            raise exc

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")
        engs = h2w * wfns.energies
        freqs = engs[1:] - engs[0]
        harm_engs = h2w * wfns.zero_order_energies
        harm_freq = harm_engs[1:] - harm_engs[0]

        gaussian_engs = [5866.872, 5785.953]
        gaussian_freqs = np.array([
            [3061.701,    2849.454],
            [2977.640,    2820.563],
            [1727.083,    1695.163],
            [1527.041,    1493.139],
            [1252.164,    1231.359],
            [1188.114,    1166.697],
            # 2 quanta
            [6123.403,    5719.555],
            [5955.281,    5578.213],
            [3454.165,    3372.197],
            [3054.082,    2988.763],
            [2504.328,    2459.646],
            [2376.227,    2327.621],
            # mixed states
            [6039.342,    5586.305],
            [4788.784,    4589.165],
            [4704.723,    4510.744],
            [4588.742,    4374.903],
            [4504.681,    4278.834],
            [3254.123,    3181.610],
            [4313.865,    4114.701],
            [4229.804,    4043.264],
            [2979.247,    2967.726],
            [2779.205,    2710.416],
            [4249.815,    4043.462],
            [4165.754,    3978.440],
            [2915.196,    2854.766],
            [2715.155,    2657.682],
            [2440.278,    2404.805]
        ])

        print_report = False
        if print_report:
            if n_modes == 6:
                print("Gaussian Energies:\n",
                      ('0 ' * n_modes + "{:>8.3f} {:>8.3f} {:>8} {:>8}\n").format(*gaussian_engs, "-", "-"),
                      *(
                          ('{:<1.0f} ' * n_modes + "{:>8} {:>8} {:>8.3f} {:>8.3f}\n").format(*s, "-", "-", *e) for s, e
                          in
                          zip(states[1:], gaussian_freqs)
                      )
                      )
            print("State Energies:\n",
                  ('0 ' * n_modes + "{:>8.3f} {:>8.3f} {:>8} {:>8}\n").format(harm_engs[0], engs[0], "-", "-"),
                  *(
                      ('{:<1.0f} ' * n_modes + "{:>8} {:>8} {:>8.3f} {:>8.3f}\n").format(*s, "-", "-", e1, e2) for
                      s, e1, e2 in
                      zip(states[1:], harm_freq, freqs)
                  )
                  )
        ns = len(states) - 1
        print_difference = True
        if print_difference:
            print("Difference Energies:\n",
                  ('0 ' * n_modes + "{:>8.3f} {:>8.3f} {:>8} {:>8}\n").format(harm_engs[0] - gaussian_engs[0],
                                                                              engs[0] - gaussian_engs[1], "-", "-"),
                  *(
                      ('{:<1.0f} ' * n_modes + "{:>8} {:>8} {:>8.3f} {:>8.3f}\n").format(*s, "-", "-", e1, e2) for
                      s, e1, e2 in
                      zip(states[1:], harm_freq[:ns] - gaussian_freqs[:ns, 0], freqs[:ns] - gaussian_freqs[:ns, 1])
                  )
                  )

        self.assertLess(
            np.max(np.abs(freqs[:ns] - gaussian_freqs[:ns, 1])),
            1)

    @inactiveTest
    def test_OCHDVPTInternals(self):

        internals = [
            [0, -1, -1, -1],
            [1,  0, -1, -1],
            [2,  1,  0, -1],
            [3,  1,  0,  2]
        ]

        n_modes = 3*4 - 6
        mode_selection = None#[5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)

        states = self.get_states(2, n_modes)
        coupled_states = self.get_states(5, n_modes)
        def block(self=self,
                  internals=internals,
                  states=states,
                  coupled_states=coupled_states,
                  mode_selection=mode_selection
                  ):
            return self.get_VPT2_wfns(
                "OCHD_freq.fchk",
                internals,
                states,
                regenerate=True,
                coupled_states=coupled_states,
                mode_selection=mode_selection
                # , degeneracies=(
                #     coupled_states.index((0, 0, 0, 0, 0, 1)),
                #     coupled_states.index((0, 0, 2, 0, 0, 0))
                # )
            )

        # block()
        exc, stat_block, wfns = self.profile_block(block)

        do_profile = False
        if do_profile:
            if exc is not None:
                try:
                    raise exc
                except:
                    raise Exception(stat_block)
            else:
                raise Exception(stat_block)
        elif exc is not None:
            raise exc

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")
        engs = h2w * wfns.energies
        freqs = engs[1:] - engs[0]
        harm_engs = h2w * wfns.zero_order_energies
        harm_freq = harm_engs[1:] - harm_engs[0]

        gaussian_engs = [5235.162,    5171.477]
        gaussian_freqs = np.array([
            [3022.813, 2881.940],
            [2221.268, 2152.210],
            [1701.548, 1673.311],
            [1417.539, 1388.280],
            [1076.474, 1058.211],
            [1030.681, 1015.305],
            # 2 quanta
            [6045.626, 5574.482],
            [4442.536, 4203.337],
            [3403.097, 3327.254],
            [2835.078, 2744.748],
            [2152.949, 2098.305],
            [2061.363, 2022.956],
            # mixed states
            [5244.081, 4990.144],
            [4724.361, 4529.189],
            [3922.817, 3809.397],
            [4440.352, 4192.841],
            [3638.807, 3518.724],
            [3119.087, 3060.953],
            [4099.287, 3896.226],
            [3297.743, 3183.261],
            [2778.023, 2725.597],
            [2494.013, 2448.626],
            [4053.494, 3866.190],
            [3251.950, 3146.116],
            [2732.230, 2681.778],
            [2448.220, 2403.196],
            [2107.156, 2074.925]
        ])

        print_report = False
        if print_report:
            if n_modes == 6:
                print("Gaussian Energies:\n",
                      ('0 ' * n_modes + "{:>8.3f} {:>8.3f} {:>8} {:>8}\n").format(*gaussian_engs, "-", "-"),
                      *(
                          ('{:<1.0f} ' * n_modes + "{:>8} {:>8} {:>8.3f} {:>8.3f}\n").format(*s, "-", "-", *e) for s, e in
                            zip(states[1:], gaussian_freqs)
                          )
                      )
            print("State Energies:\n",
                  ('0 '*n_modes + "{:>8.3f} {:>8.3f} {:>8} {:>8}\n").format(harm_engs[0], engs[0], "-", "-"),
                  *(
                      ('{:<1.0f} '*n_modes + "{:>8} {:>8} {:>8.3f} {:>8.3f}\n").format(*s, "-", "-", e1, e2) for s, e1, e2 in
                        zip(states[1:], harm_freq, freqs)
                  )
            )
        ns = len(states) - 1
        print_difference = True
        if print_difference:
            print("Difference Energies:\n",
                  ('0 ' * n_modes + "{:>8.3f} {:>8.3f} {:>8} {:>8}\n").format(harm_engs[0] - gaussian_engs[0],
                                                                              engs[0] - gaussian_engs[1], "-", "-"),
                  *(
                      ('{:<1.0f} ' * n_modes + "{:>8} {:>8} {:>8.3f} {:>8.3f}\n").format(*s, "-", "-", e1, e2) for
                      s, e1, e2 in
                      zip(states[1:], harm_freq[:ns] - gaussian_freqs[:ns, 0], freqs[:ns] - gaussian_freqs[:ns, 1])
                  )
                  )

        self.assertLess(
            np.max(np.abs(freqs[:ns] - gaussian_freqs[:ns, 1])),
            1)

    @inactiveTest
    def test_OCHDVPTCartesians(self):

        internals = None

        n_modes = 3 * 4 - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)

        states = self.get_states(2, n_modes)
        coupled_states = None#self.get_states(5, n_modes)

        # def block(self=self,
        #           internals=internals,
        #           states=states,
        #           coupled_states=coupled_states,
        #           mode_selection=mode_selection
        #           ):
        #     return

        with BlockProfiler("OCHD Cartesians",
                           print_res=False,#'raise',
                           strip_dirs=self.job_is_dumb,
                           sort_by='cumtime'
                           ):
            wfns = self.get_VPT2_wfns(
                "OCHD_freq.fchk",
                internals,
                states,
                regenerate=True,
                coupled_states=coupled_states,
                mode_selection=mode_selection,
                log=True
            )

        # block()
        # exc, stat_block, wfns = self.profile_block(block)
        #
        # do_profile = False
        # if do_profile:
        #     if exc is not None:
        #         try:
        #             raise exc
        #         except:
        #             raise Exception(stat_block)
        #     else:
        #         raise Exception(stat_block)
        # elif exc is not None:
        #     raise exc

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")
        engs = h2w * wfns.energies
        freqs = engs[1:] - engs[0]
        harm_engs = h2w * wfns.zero_order_energies
        harm_freq = harm_engs[1:] - harm_engs[0]

        gaussian_engs = [5235.162, 5171.477]
        gaussian_freqs = np.array([
            [3022.813, 2881.940],
            [2221.268, 2152.210],
            [1701.548, 1673.311],
            [1417.539, 1388.280],
            [1076.474, 1058.211],
            [1030.681, 1015.305],
            # 2 quanta
            [6045.626, 5574.482],
            [4442.536, 4203.337],
            [3403.097, 3327.254],
            [2835.078, 2744.748],
            [2152.949, 2098.305],
            [2061.363, 2022.956],
            # mixed states
            [5244.081, 4990.144],
            [4724.361, 4529.189],
            [3922.817, 3809.397],
            [4440.352, 4192.841],
            [3638.807, 3518.724],
            [3119.087, 3060.953],
            [4099.287, 3896.226],
            [3297.743, 3183.261],
            [2778.023, 2725.597],
            [2494.013, 2448.626],
            [4053.494, 3866.190],
            [3251.950, 3146.116],
            [2732.230, 2681.778],
            [2448.220, 2403.196],
            [2107.156, 2074.925]
        ])


        print_report = False
        if print_report:
            if n_modes == 6:
                print("Gaussian Energies:\n",
                      ('0 ' * n_modes + "{:>8.3f} {:>8.3f} {:>8} {:>8}\n").format(*gaussian_engs, "-", "-"),
                      *(
                          ('{:<1.0f} ' * n_modes + "{:>8} {:>8} {:>8.3f} {:>8.3f}\n").format(*s, "-", "-", *e) for s, e
                      in zip(states[1:], gaussian_freqs)
                      )
                      )
            print("State Energies:\n",
                  ('0 ' * n_modes + "{:>8.3f} {:>8.3f} {:>8} {:>8}\n").format(harm_engs[0], engs[0], "-", "-"),
                  *(
                      ('{:<1.0f} ' * n_modes + "{:>8} {:>8} {:>8.3f} {:>8.3f}\n").format(*s, "-", "-", e1, e2) for
                  s, e1, e2 in
                      zip(states[1:], harm_freq, freqs)
                  )
                  )
        ns = len(states) - 1
        print_difference = True
        if print_difference:
            print("Difference Energies:\n",
                  ('0 ' * n_modes + "{:>8.3f} {:>8.3f} {:>8} {:>8}\n").format(harm_engs[0] - gaussian_engs[0],
                                                                              engs[0] - gaussian_engs[1], "-", "-"),
                  *(
                      ('{:<1.0f} ' * n_modes + "{:>8} {:>8} {:>8.3f} {:>8.3f}\n").format(*s, "-", "-", e1, e2) for
                      s, e1, e2 in
                      zip(states[1:], harm_freq[:ns] - gaussian_freqs[:ns, 0], freqs[:ns] - gaussian_freqs[:ns, 1])
                  )
                  )

        self.assertLess(
            np.max(np.abs(freqs[:ns] - gaussian_freqs[:ns, 1])),
            1)

    @inactiveTest
    def test_OCHTVPTInternals(self):

        internals = [
            [0, -1, -1, -1],
            [1,  0, -1, -1],
            [2,  1,  0, -1],
            [3,  1,  0,  2]
        ]

        n_modes = 3 * 4 - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)

        states = self.get_states(2, n_modes)
        coupled_states = self.get_states(5, n_modes, max_quanta=5)

        def block(self=self,
                  internals=internals,
                  states=states,
                  coupled_states=coupled_states,
                  mode_selection=mode_selection
                  ):
            return self.get_VPT2_wfns(
                "OCHT_freq.fchk",
                internals,
                states,
                regenerate=True,
                coupled_states=coupled_states,
                mode_selection=mode_selection
                # , t2=0
                # , watson=0
                # , degeneracies=1/self.h2w
                # , degeneracies=(
                #     coupled_states.index((0, 0, 0, 0, 0, 1)),
                #     coupled_states.index((0, 0, 2, 0, 0, 0))
                # )
            )

        # block()
        exc, stat_block, wfns = self.profile_block(block)

        do_profile = False
        if do_profile:
            if exc is not None:
                try:
                    raise exc
                except:
                    raise Exception(stat_block)
            else:
                raise Exception(stat_block)
        elif exc is not None:
            raise exc

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")
        engs = h2w * wfns.energies
        freqs = engs[1:] - engs[0]
        harm_engs = h2w * wfns.zero_order_energies
        harm_freq = harm_engs[1:] - harm_engs[0]

        gaussian_engs = [4970.730, 4912.805]
        gaussian_freqs = np.array([
            [3022.143,  2864.375],
            [1918.480,  1859.357],
            [1660.249,  1630.227],
            [1399.294,  1370.966],
            [1036.668,  1019.245],
            [ 904.625,   892.062],

            [6044.287,  5601.333],
            [3836.959,  3677.621],
            [3320.498,  3246.038],
            [2798.588,  2724.050],
            [2073.336,  2035.811],
            [1809.250,  1778.733],

            [4940.623,  4725.335],
            [4682.392,  4497.932],
            [3578.729,  3474.240],
            [4421.438,  4220.692],
            [3317.774,  3227.798],
            [3059.543,  3000.454],
            [4058.812,  3870.413],
            [2955.148,  2867.701],
            [2696.917,  2644.455],
            [2435.962,  2393.534],
            [3926.768,  3756.997],
            [2823.105,  2746.008],
            [2564.874,  2515.823],
            [2303.919,  2263.695],
            [1941.293,  1911.106]
        ])

        print_report = False
        if print_report:
            if n_modes == 6:
                print("Gaussian Energies:\n",
                      ('0 ' * n_modes + "{:>8.3f} {:>8.3f} {:>8} {:>8}\n").format(*gaussian_engs, "-", "-"),
                      *(
                          ('{:<1.0f} ' * n_modes + "{:>8} {:>8} {:>8.3f} {:>8.3f}\n").format(*s, "-", "-", *e) for s, e
                      in
                          zip(states[1:], gaussian_freqs)
                      )
                      )
            print("State Energies:\n",
                  ('0 ' * n_modes + "{:>8.3f} {:>8.3f} {:>8} {:>8}\n").format(harm_engs[0], engs[0], "-", "-"),
                  *(
                      ('{:<1.0f} ' * n_modes + "{:>8} {:>8} {:>8.3f} {:>8.3f}\n").format(*s, "-", "-", e1, e2) for
                  s, e1, e2 in
                      zip(states[1:], harm_freq, freqs)
                  )
                  )

        print_difference = True
        if print_difference:
            print("Difference Energies:\n",
                  ('0 ' * n_modes + "{:>8.3f} {:>8.3f} {:>8} {:>8}\n").format(harm_engs[0] - gaussian_engs[0],
                                                                              engs[0] - gaussian_engs[1], "-", "-"),
                  *(
                      ('{:<1.0f} ' * n_modes + "{:>8} {:>8} {:>8.3f} {:>8.3f}\n").format(*s, "-", "-", e1, e2) for
                      s, e1, e2 in
                      zip(states[1:], harm_freq - gaussian_freqs[:, 0], freqs - gaussian_freqs[:, 1])
                  )
                  )

        self.assertLess(
            np.max(np.abs(freqs-gaussian_freqs[:, 1])[:len(states)-1]),
            1)

    @inactiveTest
    def test_OCHTVPTCartesians(self):

        internals = None

        n_modes = 3 * 4 - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)

        states = self.get_states(2, n_modes)

        def block(self=self,
                  internals=internals,
                  states=states,
                  mode_selection=mode_selection
                  ):
            return self.get_VPT2_wfns(
                "OCHT_freq.fchk",
                internals,
                states,
                regenerate=True,
                # coupled_states=self.get_states(9, n_modes),
                # freq_threshold=10000/self.h2w, #only couple stuff within 5000 wavenumbers
                # coupled_states=coupled_states,
                mode_selection=mode_selection
                # , degeneracies = (
                #     (
                #         (0, 0, 0, 0, 0, 1),
                #     ),
                #     (
                #         (0, 0, 1, 1, 0, 0),
                #     )
                # )
                # , degeneracies=50/self.h2w
                , log=True
            )

        # block()
        exc, stat_block, wfns = self.profile_block(block)
        do_profile = False
        if do_profile:
            if exc is not None:
                try:
                    raise exc
                except:
                    raise Exception(stat_block)
            else:
                raise Exception(stat_block)
        elif exc is not None:
            raise exc


        # print(len(coupled_states), len(wfns.corrs.coupled_states))
        # print(wfns.corrs.degenerate_states)

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")
        engs = h2w * wfns.energies
        # print([
        #     np.max(np.abs(wfns.corrs.wfn_corrections[i, 1]))
        #     for i in range(len(wfns.corrs.wfn_corrections))
        #     ])
        # print([
        #     np.max(np.abs(wfns.corrs.wfn_corrections[i, 2]))
        #     for i in range(len(wfns.corrs.wfn_corrections))
        # ])
        freqs = engs[1:] - engs[0]
        harm_engs = h2w * wfns.zero_order_energies
        harm_freq = harm_engs[1:] - harm_engs[0]

        gaussian_engs = [4970.730, 4912.805]
        gaussian_freqs = np.array([
            [3022.143, 2864.375],
            [1918.480, 1859.357],
            [1660.249, 1630.227],
            [1399.294, 1370.966],
            [1036.668, 1019.245],
            [ 904.625,  892.062],

            [6044.287, 5601.333],
            [3836.959, 3677.621],
            [3320.498, 3246.038],
            [2798.588, 2724.050],
            [2073.336, 2035.811],
            [1809.250, 1778.733],

            [4940.623, 4725.335],
            [4682.392, 4497.932],
            [3578.729, 3474.240],
            [4421.438, 4220.692],
            [3317.774, 3227.798],
            [3059.543, 3000.454],
            [4058.812, 3870.413],
            [2955.148, 2867.701],
            [2696.917, 2644.455],
            [2435.962, 2393.534],
            [3926.768, 3756.997],
            [2823.105, 2746.008],
            [2564.874, 2515.823],
            [2303.919, 2263.695],
            [1941.293, 1911.106]
        ])

        ns = len(states)-1
        print_report = True
        if print_report:
            if n_modes == 6:
                print("Gaussian Energies:\n",
                      ('0 ' * n_modes + "{:>8.3f} {:>8.3f} {:>8} {:>8}\n").format(*gaussian_engs, "-", "-"),
                      *(
                          ('{:<1.0f} ' * n_modes + "{:>8} {:>8} {:>8.3f} {:>8.3f}\n").format(*s, "-", "-", *e) for s, e
                          in
                          zip(states[1:], gaussian_freqs[:ns])
                      )
                      )
            print("State Energies:\n",
                  ('0 ' * n_modes + "{:>8.3f} {:>8.3f} {:>8} {:>8}\n").format(harm_engs[0], engs[0], "-", "-"),
                  *(
                      ('{:<1.0f} ' * n_modes + "{:>8} {:>8} {:>8.3f} {:>8.3f}\n").format(*s, "-", "-", e1, e2) for
                      s, e1, e2 in
                      zip(states[1:], harm_freq[:ns], freqs[:ns])
                  )
                  )
        print_difference = True
        if print_difference:
            print("Difference Energies:\n",
                  ('0 ' * n_modes + "{:>8.3f} {:>8.3f} {:>8} {:>8}\n").format(harm_engs[0] - gaussian_engs[0],
                                                                              engs[0] - gaussian_engs[1], "-", "-"),
                  *(
                      ('{:<1.0f} ' * n_modes + "{:>8} {:>8} {:>8.3f} {:>8.3f}\n").format(*s, "-", "-", e1, e2) for
                      s, e1, e2 in
                      zip(states[1:], harm_freq[:ns] - gaussian_freqs[:ns, 0], freqs[:ns] - gaussian_freqs[:ns, 1])
                  )
                  )

        self.assertLess(
            np.max(np.abs(freqs[:ns] - gaussian_freqs[:ns, 1])),
            1)

    @inactiveTest
    def test_CH2DTVPTCartesians(self):

        internals = None

        n_modes = 5 * 3 - 6
        mode_selection = None
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)

        states = self.get_states(2, n_modes)[:5]
        # coupled_states = self.get_states(5, n_modes, max_quanta=5)
        #

        def block(self=self,
                  internals=internals,
                  states=states,
                  coupled_states=None,
                  mode_selection=mode_selection
                  ):
            return self.get_VPT2_wfns(
                "CH2DT_freq.fchk",
                internals,
                states,
                regenerate=True,
                coupled_states=coupled_states,
                mode_selection=mode_selection,
                log=True
            )

        # block()
        exc, stat_block, wfns = self.profile_block(block)

        do_profile = False
        if do_profile:
            if exc is not None:
                try:
                    raise exc
                except:
                    raise Exception(stat_block)
            else:
                raise Exception(stat_block)
        elif exc is not None:
            raise exc

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")
        engs = h2w * wfns.energies
        freqs = engs[1:] - engs[0]
        harm_engs = h2w * wfns.zero_order_energies
        harm_freq = harm_engs[1:] - harm_engs[0]

        gaussian_engs = [8311.227, 8205.865]
        gaussian_freqs = np.array([
             [3206.212, 3056.426],
             [3141.846, 3001.995],
             [2328.234, 2246.906],
             [1951.131, 1907.512],
             [1445.225, 1413.150],
             [1311.183, 1285.762],
             [1230.305, 1205.049],
             [1048.980, 1028.622],
             [ 959.340,  944.038],

             [6412.425, 6049.075],
             [6283.691, 5956.951],
             [4656.467, 4431.987],
             [3902.262, 3737.670],
             [2890.450, 2815.938],
             [2622.365, 2568.324],
             [2460.609, 2402.604],
             [2097.961, 2053.655],
             [1918.679, 1868.498],

             [6348.058, 5946.349],
             [5534.446, 5304.546],
             [5470.079, 5252.170],
             [5157.343, 4946.234],
             [5092.977, 4895.552],
             [4279.364, 4131.666],
             [4651.438, 4449.136],
             [4587.071, 4407.847],
             [3773.459, 3657.329],
             [3396.356, 3302.019],
             [4517.395, 4325.685],
             [4453.028, 4280.602],
             [3639.416, 3518.558],
             [3262.313, 3175.698],
             [2756.408, 2698.844],
             [4436.517, 4250.290],
             [4372.150, 4204.590],
             [3558.538, 3448.813],
             [3181.436, 3094.860],
             [2675.530, 2618.649],
             [2541.487, 2483.423],
             [4255.193, 4077.953],
             [4190.826, 4027.168],
             [3377.214, 3269.456],
             [3000.111, 2904.013],
             [2494.205, 2435.125],
             [2360.163, 2319.432],
             [2279.285, 2236.240],
             [4165.552, 3999.186],
             [4101.185, 3948.975],
             [3287.573, 3183.220],
             [2910.471, 2822.890],
             [2404.565, 2356.088],
             [2270.522, 2230.076],
             [2189.644, 2143.874],
             [2008.320, 1973.845]
        ])

        print_report = False
        if print_report:
            print("Gaussian Energies:\n",
                  ('0 ' * n_modes + "{:>8.3f} {:>8.3f} {:>8} {:>8}\n").format(*gaussian_engs, "-", "-"),
                  *(
                      ('{:<1.0f} ' * n_modes + "{:>8} {:>8} {:>8.3f} {:>8.3f}\n").format(*s, "-", "-", *e) for s, e
                      in
                      zip(states[1:], gaussian_freqs)
                  )
                )
            print("State Energies:\n",
                  ('0 ' * n_modes + "{:>8.3f} {:>8.3f} {:>8} {:>8}\n").format(harm_engs[0], engs[0], "-", "-"),
                  *(
                      ('{:<1.0f} ' * n_modes + "{:>8} {:>8} {:>8.3f} {:>8.3f}\n").format(*s, "-", "-", e1, e2) for
                      s, e1, e2 in
                      zip(states[1:], harm_freq, freqs)
                  )
                  )
        ns = len(states) - 1
        print_difference = True
        if print_difference:
            print("Difference Energies:\n",
                  ('0 ' * n_modes + "{:>8.3f} {:>8.3f} {:>8} {:>8}\n").format(harm_engs[0] - gaussian_engs[0],
                                                                              engs[0] - gaussian_engs[1], "-", "-"),
                  *(
                      ('{:<1.0f} ' * n_modes + "{:>8} {:>8} {:>8.3f} {:>8.3f}\n").format(*s, "-", "-", e1, e2) for
                      s, e1, e2 in
                      zip(states[1:], harm_freq[:ns] - gaussian_freqs[:ns, 0], freqs[:ns] - gaussian_freqs[:ns, 1])
                  )
                  )

        self.assertLess(
            np.max(np.abs(freqs[:ns] - gaussian_freqs[:ns, 1])),
            1)

    @inactiveTest
    def test_WaterDimerVPTCartesians(self):

        internals = None

        n_modes = 6 * 3 - 6
        mode_selection = None
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)

        states = self.get_states(2, n_modes)[:24]

        # coupled_states = self.get_states(5, n_modes, max_quanta=5)
        #

        def block(self=self,
                  internals=internals,
                  states=states,
                  coupled_states=None,
                  mode_selection=mode_selection
                  ):
            return self.get_VPT2_wfns(
                "water_dimer_freq.fchk",
                internals,
                states,
                regenerate=True,
                coupled_states=coupled_states,
                mode_selection=mode_selection,
                log=True
            )

        # block()
        exc, stat_block, wfns = self.profile_block(block)

        do_profile = False
        if do_profile:
            if exc is not None:
                try:
                    raise exc
                except:
                    raise Exception(stat_block)
            else:
                raise Exception(stat_block)
        elif exc is not None:
            raise exc

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")
        engs = h2w * wfns.energies
        freqs = engs[1:] - engs[0]
        harm_engs = h2w * wfns.zero_order_energies
        harm_freq = harm_engs[1:] - harm_engs[0]

        gaussian_engs = [10133.860,    9909.756]
        gaussian_freqs = np.array([
             [3935.490,    3742.918],
             [3914.939,    3752.151],
             [3814.079,    3652.414],
             [3718.192,    3584.139],
             [1650.023,    1592.653],
             [1629.210,    1585.962],
             [ 631.340,     505.605],
             [ 362.703,     295.834],
             [ 183.777,     141.372],
             [ 154.306,     110.995],
             [ 146.544,     150.517],
             [ 127.117,      69.163],

             [7870.980,    7393.560],
             [7829.879,    7368.493],
             [7628.159,    7224.882],
             [7436.384,    7016.025],
             [3300.045,    3152.473],
             [3258.421,    3144.157],
             [1262.679,     921.053],
             [ 725.405,     488.907],
             [ 367.554,     268.882],
             [ 308.612,     207.465],
             [ 293.089,     299.766],
             [ 254.234,     114.677],

             [7729.019,    7402.976],
             [7633.131,    7271.264],
             [7532.271,    7230.663],
             [5564.962,    5328.224],
             [5464.102,    5244.056],
             [5368.215,    5164.597],
             [5544.150,    5337.031],
             [5443.290,    5222.111],
             [5347.402,    5168.407],
             [3279.233,    3172.374],
             [4277.642,    4024.523],
             [4176.782,    3928.515],
             [4080.894,    3889.457],
             [2012.725,    1852.952],
             [1991.913,    1862.320],
             [4098.716,    3895.805],
             [3997.856,    3794.279],
             [3901.969,    3739.502],
             [1833.800,    1732.294],
             [1812.987,    1729.354],
             [ 546.479,     389.065],
             [4069.245,    3864.621],
             [3968.385,    3763.445],
             [3872.498,    3704.835],
             [1804.329,    1699.128],
             [1783.516,    1700.178],
             [ 517.009,     374.506],
             [ 338.083,     235.655],
             [7850.429,    7494.650],
             [7749.569,    7239.308],
             [7653.682,    7322.974],
             [5585.513,    5334.869],
             [5564.700,    5314.396],
             [4298.193,    4016.063],
             [4119.267,    3875.578],
             [4089.796,    3839.370],
             [4546.279,    4261.388],
             [4445.419,    4159.774],
             [4349.531,    4139.077],
             [2281.362,    2107.393],
             [2260.550,    2094.511],
             [ 994.042,     745.791],
             [ 815.116,     620.620],
             [ 785.646,     595.362],
             [4566.830,    4249.695],
             [4061.484,    3903.055],
             [3960.624,    3844.234],
             [3864.736,    3750.792],
             [1796.567,    1745.999],
             [1775.755,    1736.874],
             [ 509.247,     405.832],
             [ 330.321,     287.882],
             [ 300.850,     261.693],
             [4082.035,    3892.356],
             [ 777.884,     646.479],
             [4042.056,    3844.889],
             [3941.197,    3692.207],
             [3845.309,    3657.100],
             [1777.140,    1656.439],
             [1756.328,    1651.982],
             [ 489.820,     336.402],
             [ 310.894,     207.357],
             [ 281.423,     169.590],
             [4062.607,    3810.202],
             [ 758.457,     525.680],
             [ 273.662,     201.333]
        ])

        print_report = True
        if print_report:
            print("Gaussian Energies:\n",
                  ('0 ' * n_modes + "{:>8.3f} {:>8.3f} {:>8} {:>8}\n").format(*gaussian_engs, "-", "-"),
                  *(
                      ('{:<1.0f} ' * n_modes + "{:>8} {:>8} {:>8.3f} {:>8.3f}\n").format(*s, "-", "-", *e) for s, e
                      in
                      zip(states[1:], gaussian_freqs)
                  )
                  )
            print("State Energies:\n",
                  ('0 ' * n_modes + "{:>8.3f} {:>8.3f} {:>8} {:>8}\n").format(harm_engs[0], engs[0], "-", "-"),
                  *(
                      ('{:<1.0f} ' * n_modes + "{:>8} {:>8} {:>8.3f} {:>8.3f}\n").format(*s, "-", "-", e1, e2) for
                      s, e1, e2 in
                      zip(states[1:], harm_freq, freqs)
                  )
                  )
        ns = len(states) - 1
        print_difference = True
        if print_difference:
            print("Difference Energies:\n",
                  ('0 ' * n_modes + "{:>8.3f} {:>8.3f} {:>8} {:>8}\n").format(harm_engs[0] - gaussian_engs[0],
                                                                              engs[0] - gaussian_engs[1], "-", "-"),
                  *(
                      ('{:<1.0f} ' * n_modes + "{:>8} {:>8} {:>8.3f} {:>8.3f}\n").format(*s, "-", "-", e1, e2) for
                      s, e1, e2 in
                      zip(states[1:], harm_freq[:ns] - gaussian_freqs[:ns, 0], freqs[:ns] - gaussian_freqs[:ns, 1])
                  )
                  )

        self.assertLess(
            np.max(np.abs(freqs[:ns] - gaussian_freqs[:ns, 1])),
            1)

    #endregion Test Systems

    #region Test Action Expansions

    @debugTest
    def test_HOHCartesianActionExpansion(self):

        internals = None

        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("HOH_freq.fchk"),
            internals=internals
        )

        expansion, wfns = ham.get_action_expansion()

        g0, w, x = expansion


        raise Exception([
            self.h2w * x,
            self.h2w * w,
            self.h2w * g0
        ])

    @inactiveTest
    def test_OCHHCartesianActionExpansion(self):

        internals = None

        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("OCHH_freq.fchk"),
            internals=internals
        )

        expansion, wfns = ham.get_action_expansion()

        g0, w, x = expansion

        wax = ham.get_Nielsen_xmatrix()
        wat = np.sum(wax, axis=0)
        wat_vib = np.sum(wax[:2], axis=0)

        raise Exception([
            self.h2w * x,
            wat * self.h2w,
            wax[0] * self.h2w,
            wax[1] * self.h2w,
            wax[2] * self.h2w
            # self.h2w * w,
            # self.h2w * g0
        ])

    @inactiveTest
    def test_WaterDimerCartesianActionExpansion(self):

        internals = None

        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("water_dimer_freq.fchk"),
            internals=internals
        )

        expansion, wfns = ham.get_action_expansion()

        g0, w, x = expansion

        wat = np.sum(ham.get_Nielsen_xmatrix(), axis=0)

        raise Exception([
            self.h2w * x,
            wat * self.h2w
            # self.h2w * w,
            # self.h2w * g0
        ])

        # raise Exception([
        #     self.h2w * x,
        #     self.h2w * w,
        #     self.h2w * g0
        # ])

    #endregion Test Action Expansions

    #region Test Terms
    @validationTest
    def test_HODVPTCartesianPotential(self):

        internals = None

        n_modes = 3 * 3 - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)

        max_quanta = 4
        states = (
            (0, 0, 0),
        )  # we're only interested in the correction of the ground state w/o perturbation or accounting for the energy or any shit
        coupled_states = self.get_states(5, n_modes, max_quanta=5)

        def block(self=self,
                  internals=internals,
                  states=states,
                  coupled_states=coupled_states,
                  mode_selection=mode_selection
                  ):
            return self.get_VPT2_wfns(
                "HOD_freq.fchk",
                internals,
                states,
                regenerate=True,
                coupled_states=coupled_states,
                mode_selection=mode_selection
                # , v3=0
                # , v4=0
                , coriolis=0
                # , degeneracies=(
                #     coupled_states.index((0, 0, 0, 0, 0, 1)),
                #     coupled_states.index((0, 0, 2, 0, 0, 0))
                # )
            )

        # block()
        exc, stat_block, wfns = self.profile_block(block)

        do_profile = False
        if do_profile:
            if exc is not None:
                try:
                    raise exc
                except:
                    raise Exception(stat_block)
            else:
                raise Exception(stat_block)
        elif exc is not None:
            raise exc

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")
        eng_corrs = h2w * wfns.corrs.energy_corrs[0]
        gaussian_corrs = np.array([4052.91093, 0., -50.05456])

        self.assertLess(np.max(np.abs(eng_corrs - gaussian_corrs)), .5)

    @validationTest
    def test_HODVPTCartesianCoriolis(self):

        internals = None

        n_modes = 3 * 3 - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)

        max_quanta = 4
        states = (
            (0, 0, 0),
        )  # we're only interested in the correction of the ground state w/o perturbation or accounting for the energy or any shit
        coupled_states = self.get_states(4, n_modes, max_quanta=max_quanta)

        def block(self=self,
                  internals=internals,
                  states=states,
                  coupled_states=coupled_states,
                  mode_selection=mode_selection
                  ):
            return self.get_VPT2_wfns(
                "HOD_freq.fchk",
                internals,
                states,
                regenerate=True,
                coupled_states=coupled_states,
                mode_selection=mode_selection
                , v3=0
                , v4=0
                # , coriolis=np.ones((3, 3, 3, 3))
            )

        # block()
        exc, stat_block, wfns = self.profile_block(block)

        do_profile = False
        if do_profile:
            if exc is not None:
                try:
                    raise exc
                except:
                    raise Exception(stat_block)
            else:
                raise Exception(stat_block)
        elif exc is not None:
            raise exc

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")
        eng_corrs = h2w * wfns.corrs.energy_corrs[0]
        gaussian_corrs = np.array([4052.91093, 0., -8.01373])

        self.assertLess(np.max(np.abs(eng_corrs - gaussian_corrs)), .5)

    @validationTest
    def test_HOHCartesianPotential(self):

        internals = None

        n_modes = 3 * 3 - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)

        max_quanta = 4
        states = (
            (0, 0, 0),
        )  # we're only interested in the correction of the ground state w/o perturbation or accounting for the energy or any shit
        coupled_states = self.get_states(4, n_modes, max_quanta=max_quanta)

        def block(self=self,
                  internals=internals,
                  states=states,
                  coupled_states=coupled_states,
                  mode_selection=mode_selection
                  ):
            return self.get_VPT2_wfns(
                "HOH_freq.fchk",
                internals,
                states,
                regenerate=True,
                coupled_states=coupled_states,
                mode_selection=mode_selection
                # , v3=0
                # , v4=0
                , coriolis=0
                # , degeneracies=(
                #     coupled_states.index((0, 0, 0, 0, 0, 1)),
                #     coupled_states.index((0, 0, 2, 0, 0, 0))
                # )
            )

        # block()
        exc, stat_block, wfns = self.profile_block(block)

        do_profile = False
        if do_profile:
            if exc is not None:
                try:
                    raise exc
                except:
                    raise Exception(stat_block)
            else:
                raise Exception(stat_block)
        elif exc is not None:
            raise exc

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")
        eng_corrs = h2w * wfns.corrs.energy_corrs[0]
        gaussian_corrs = np.array([4681.56363, 0., -64.98016])

        self.assertLess(np.max(np.abs(eng_corrs - gaussian_corrs)), .5)

    @validationTest
    def test_HOHCartesianCoriolis(self):

        internals = None

        n_modes = 3 * 3 - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)

        max_quanta = 4
        states = (
            (0, 0, 0),
        )  # we're only interested in the correction of the ground state w/o perturbation or accounting for the energy or any shit
        coupled_states = self.get_states(4, n_modes, max_quanta=max_quanta)

        def block(self=self,
                  internals=internals,
                  states=states,
                  coupled_states=coupled_states,
                  mode_selection=mode_selection
                  ):
            return self.get_VPT2_wfns(
                "HOH_freq.fchk",
                internals,
                states,
                regenerate=True,
                coupled_states=coupled_states,
                mode_selection=mode_selection
                , v3=0
                , v4=0
                # , coriolis=np.ones((3, 3, 3, 3))
            )

        # block()
        exc, stat_block, wfns = self.profile_block(block)

        do_profile = False
        if do_profile:
            if exc is not None:
                try:
                    raise exc
                except:
                    raise Exception(stat_block)
            else:
                raise Exception(stat_block)
        elif exc is not None:
            raise exc

        h2w = self.h2w
        eng_corrs = h2w * wfns.corrs.energy_corrs[0]
        gauss_corrs = np.array([4681.56363, 0., -10.63019]) # rounding _just_ barely fucks up on thos

        self.assertLess(np.max(np.abs(eng_corrs - gauss_corrs)), .5) # less than .5 wavenumbers off

    @validationTest
    def test_OCHDVPTCartesiansPotential(self):

        internals = None

        n_modes = 3 * 4 - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)

        max_quanta = 4
        states = (
            (0, 0, 0, 0, 0, 0),
            (0, 0, 0, 0, 0, 1)
        )  # we're only interested in the correction of the ground state w/o perturbation or accounting for the energy or any shit
        coupled_states = [
            x for x in self.get_states(4, n_modes, max_quanta=max_quanta)
            if all(x[i] < 3 for i in range(5))
        ]

        def block(self=self,
                  internals=internals,
                  states=states,
                  coupled_states=coupled_states,
                  mode_selection=mode_selection
                  ):
            return self.get_VPT2_wfns(
                "OCHD_freq.fchk",
                internals,
                states,
                regenerate=True,
                coupled_states=coupled_states,
                mode_selection=mode_selection
                , coriolis=0
            )

        # block()
        exc, stat_block, wfns = self.profile_block(block)

        do_profile = False
        if do_profile:
            if exc is not None:
                try:
                    raise exc
                except:
                    raise Exception(stat_block)
            else:
                raise Exception(stat_block)
        elif exc is not None:
            raise exc

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")
        eng_corrs = h2w * wfns.corrs.energy_corrs[0]
        gaussian_corrs = np.array([5235.16208, 0., -63.06083])

        self.assertLess(np.max(np.abs(eng_corrs - gaussian_corrs)), .5)

    @validationTest
    def test_OCHDVPTCartesiansCoriolis(self):

        internals = None

        n_modes = 3 * 4 - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)

        max_quanta = 4
        states = (
            (0, 0, 0, 0, 0, 0),
            (0, 0, 0, 0, 0, 1)
        )  # we're only interested in the correction of the ground state w/o perturbation or accounting for the energy or any shit
        coupled_states = [
            x for x in self.get_states(4, n_modes, max_quanta=max_quanta)
            if all(x[i] < 3 for i in range(5))
        ]

        def block(self=self,
                  internals=internals,
                  states=states,
                  coupled_states=coupled_states,
                  mode_selection=mode_selection
                  ):
            return self.get_VPT2_wfns(
                "OCHD_freq.fchk",
                internals,
                states,
                regenerate=True,
                coupled_states=coupled_states,
                mode_selection=mode_selection
                , v3=0
                , v4=0
            )

        # block()
        exc, stat_block, wfns = self.profile_block(block)

        do_profile = False
        if do_profile:
            if exc is not None:
                try:
                    raise exc
                except:
                    raise Exception(stat_block)
            else:
                raise Exception(stat_block)
        elif exc is not None:
            raise exc

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")
        eng_corrs = h2w * wfns.corrs.energy_corrs[0]
        gaussian_corrs = np.array([5235.16208, 0., -0.62420])

        self.assertLess(np.max(np.abs(eng_corrs - gaussian_corrs)), .5)

    @validationTest
    def test_OCHTCartesianPotential(self):

        # This is _expected_ to fail

        internals = None

        n_modes = 3 * 4 - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)

        max_quanta = 4
        states = (
            (0, 0, 0, 0, 0, 0),
            # (0, 0, 0, 0, 0, 1)
        )  # we're only interested in the correction of the ground state w/o perturbation or accounting for the energy or any shit
        coupled_states = self.get_states(4, n_modes, max_quanta=max_quanta)

        def block(self=self,
                  internals=internals,
                  states=states,
                  coupled_states=coupled_states,
                  mode_selection=mode_selection
                  ):
            return self.get_VPT2_wfns(
                "OCHT_freq.fchk",
                internals,
                states,
                regenerate=True,
                coupled_states=coupled_states,
                mode_selection=mode_selection
                , coriolis=0
            )

        # block()
        exc, stat_block, wfns = self.profile_block(block)

        do_profile = False
        if do_profile:
            if exc is not None:
                try:
                    raise exc
                except:
                    raise Exception(stat_block)
            else:
                raise Exception(stat_block)
        elif exc is not None:
            raise exc

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")
        eng_corrs = h2w * wfns.corrs.energy_corrs[0]
        gaussian_corrs = np.array([4970.72973, 0., -57.45822])

        self.assertLess(np.max(np.abs(eng_corrs - gaussian_corrs)), .5)

    @validationTest
    def test_OCHTCartesianCoriolis(self):

        internals = None

        n_modes = 3 * 4 - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)

        max_quanta = 4
        states = (
            (0, 0, 0, 0, 0, 0),
            # (0, 0, 0, 0, 0, 1)
        )  # we're only interested in the correction of the ground state w/o perturbation or accounting for the energy or any shit
        coupled_states = self.get_states(4, n_modes, max_quanta=max_quanta)

        def block(self=self,
                  internals=internals,
                  states=states,
                  coupled_states=coupled_states,
                  mode_selection=mode_selection
                  ):
            return self.get_VPT2_wfns(
                "OCHT_freq.fchk",
                internals,
                states,
                regenerate=True,
                coupled_states=coupled_states,
                mode_selection=mode_selection
                , v3=0
                , v4=0
            )

        # block()
        exc, stat_block, wfns = self.profile_block(block)

        do_profile = False
        if do_profile:
            if exc is not None:
                try:
                    raise exc
                except:
                    raise Exception(stat_block)
            else:
                raise Exception(stat_block)
        elif exc is not None:
            raise exc

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")
        eng_corrs = h2w * wfns.corrs.energy_corrs[0]
        gaussian_corrs = np.array([4970.72973, 0., -0.46674])

        self.assertLess(np.max(np.abs(eng_corrs - gaussian_corrs)), .5)

    @validationTest
    def test_OCHHVPTCartesianPotential(self):

        # This is _expected_ to fail

        internals = None

        n_modes = 3 * 4 - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)

        max_quanta = 4
        states = (
            (0, 0, 0, 0, 0, 0),
            # (0, 0, 0, 0, 0, 1)
        )  # we're only interested in the correction of the ground state w/o perturbation or accounting for the energy or any shit
        coupled_states = self.get_states(4, n_modes, max_quanta=max_quanta)

        def block(self=self,
                  internals=internals,
                  states=states,
                  coupled_states=coupled_states,
                  mode_selection=mode_selection
                  ):
            return self.get_VPT2_wfns(
                "OCHH_freq.fchk",
                internals,
                states,
                regenerate=True,
                coupled_states=coupled_states,
                mode_selection=mode_selection
                , coriolis=0
            )

        # block()
        exc, stat_block, wfns = self.profile_block(block)

        do_profile = False
        if do_profile:
            if exc is not None:
                try:
                    raise exc
                except:
                    raise Exception(stat_block)
            else:
                raise Exception(stat_block)
        elif exc is not None:
            raise exc

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")
        eng_corrs = h2w * wfns.corrs.energy_corrs[0]
        gaussian_corrs = np.array([5866.87157, 0., -80.04786])

        self.assertLess(np.max(np.abs(eng_corrs - gaussian_corrs)), .5)

    @validationTest
    def test_OCHHVPTCartesianCoriolis(self):

        internals = None

        n_modes = 3 * 4 - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)

        max_quanta = 4
        states = (
            (0, 0, 0, 0, 0, 0),
            # (0, 0, 0, 0, 0, 1)
        )  # we're only interested in the correction of the ground state w/o perturbation or accounting for the energy or any shit
        coupled_states = self.get_states(4, n_modes, max_quanta=max_quanta)

        def block(self=self,
                  internals=internals,
                  states=states,
                  coupled_states=coupled_states,
                  mode_selection=mode_selection
                  ):
            return self.get_VPT2_wfns(
                "OCHH_freq.fchk",
                internals,
                states,
                regenerate=True,
                coupled_states=coupled_states,
                mode_selection=mode_selection
                , v3=0
                , v4=0
            )

        # block()
        exc, stat_block, wfns = self.profile_block(block)

        do_profile = False
        if do_profile:
            if exc is not None:
                try:
                    raise exc
                except:
                    raise Exception(stat_block)
            else:
                raise Exception(stat_block)
        elif exc is not None:
            raise exc

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")
        eng_corrs = h2w * wfns.corrs.energy_corrs[0]
        gaussian_corrs = np.array([5866.87157, 0., -0.87104])

        self.assertLess(np.max(np.abs(eng_corrs - gaussian_corrs)), .5)

    #endregion Test Components

    #region Test Intensities

    @validationTest
    def test_HODIntensities(self):

        internals = [
            [0, -1, -1, -1],
            [1,  0, -1, -1],
            [2,  0,  1, -1]
        ]
        states = (
            (0, 0, 0),
            (0, 0, 1), (0, 1, 0), (1, 0, 0),
            (0, 0, 2), (0, 2, 0), (2, 0, 0),
            (0, 1, 1), (1, 0, 1), (1, 1, 0)
        )
        coupled_states = self.get_states(5, 3, max_quanta=5)

        wfns = self.get_VPT2_wfns(
            "HOD_freq.fchk",
            internals,
            states,
            regenerate=True,
            coupled_states=coupled_states
        )

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        # trying to turn off the orthogonality condition
        # states = wfns.corrs.states
        # for i,s in enumerate(states):
        #     wfns.corrs.wfn_corrections[i, 2, s] = 0 # turn off second correction
        engs = h2w * wfns.energies
        freqs = engs - engs[0]
        ints = wfns.intensities

        harm_engs = h2w * wfns.zero_order_energies
        harm_freqs = harm_engs - harm_engs[0]
        harm_ints = wfns.zero_order_intensities

        plot_specs = False
        if plot_specs:
            import McUtils.Plots as plt
            s = plt.StickPlot(freqs, ints,
                          aspect_ratio=.5,
                          plot_legend=("Anharmonic", "Harmonic"),
                              image_size=500
                          )
            plt.StickPlot(harm_freqs, harm_ints, figure=s, plot_style=dict(linefmt='red'),
                          aspect_ratio = .5,
                          axes_labels = ["Frequency (cm$^{-1}$)", "Intensity (km mol$^{-1}$)"],
                          plot_legend=("Anharmonic", "Harmonic")
                          )
            s.show()

        # plt.TensorPlot(np.array(tms)).show()
        # print(tms[0])
        gaussian_harm_freqs = [
            0.,
            3873.846,
            2810.031,
            1421.946
        ]
        gaussian_harm_ints = [
            0.,
            40.02416012,
            17.23521259,
            57.65661496
        ]
        gaussian_freqs = [
            0.,
            3685.815,
            2706.132,
            1383.391,
            7202.835,
            5323.917,
            2749.027,
            6377.958,
            5044.721,
            4072.407
        ]
        gaussian_ints = [
            0.,
            36.93616327,
            14.01377595,
            58.17413771,
            1.16678234,
            0.55376888,
            3.30245090,
            0.19253704,
            1.57560144,
            1.82673531
        ]

        print_specs = True
        if print_specs:
            report = "Harmonic:   State     Freq.   Int.    Gaussian: Freq.   Int.\n" + "\n".join(
                " " * 12 + "{} {:>7.2f} {:>7.4f}           {:>7.2f} {:>7.4f}".format(s, f, i, gf, g)
                for s, f, i, gf, g in zip(states, harm_freqs, harm_ints, gaussian_harm_freqs, gaussian_harm_ints)
            )
            print(report)
            report = "Anharmonic: State     Freq.   Int.    Gaussian: Freq.   Int.\n" + "\n".join(
                " " * 12 + "{} {:>7.2f} {:>7.4f}           {:>7.2f} {:>7.4f}".format(s, f, i, gf, g)
                for s, f, i, gf, g in zip(states, freqs, ints, gaussian_freqs, gaussian_ints)
            )
            print(report)

        self.assertEquals(
            [round(x, 2) for x in gaussian_harm_ints[1:]],
            list(np.round(harm_ints[1:4], 2))
        )
        self.assertEquals(
            [round(x, 2) for x in gaussian_ints[1:]],
            list(np.round(ints[1:10], 2))
        )

    @validationTest
    def test_HODIntensitiesCartesian(self):

        internals = None
        states = (
            (0, 0, 0),
            (0, 0, 1), (0, 1, 0), (1, 0, 0),
            (0, 0, 2), (0, 2, 0), (2, 0, 0),
            (0, 1, 1), (1, 0, 1), (1, 1, 0)
        )
        coupled_states = self.get_states(5, 3, max_quanta=5)
        wfns = self.get_VPT2_wfns(
            "HOD_freq.fchk",
            internals,
            states,
            regenerate=True,
            coupled_states=coupled_states
        )

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        # trying to turn off the orthogonality condition
        # states = wfns.corrs.states
        # for i,s in enumerate(states):
        #     wfns.corrs.wfn_corrections[i, 2, s] = 0 # turn off second correction
        engs = h2w * wfns.energies
        freqs = engs - engs[0]
        ints = wfns.intensities

        harm_engs = h2w * wfns.zero_order_energies
        harm_freqs = harm_engs - harm_engs[0]
        harm_ints = wfns.zero_order_intensities

        plot_specs = False
        if plot_specs:
            import McUtils.Plots as plt
            s = plt.StickPlot(freqs, ints,
                              aspect_ratio=.5,
                              plot_legend=("Anharmonic", "Harmonic"),
                              image_size=500
                              )
            plt.StickPlot(harm_freqs, harm_ints, figure=s, plot_style=dict(linefmt='red'),
                          aspect_ratio=.5,
                          axes_labels=["Frequency (cm$^{-1}$)", "Intensity (km mol$^{-1}$)"],
                          plot_legend=("Anharmonic", "Harmonic")
                          )
            s.show()

        # plt.TensorPlot(np.array(tms)).show()
        # print(tms[0])
        gaussian_harm_freqs = [
            0.,
            3873.846,
            2810.031,
            1421.946
        ]
        gaussian_harm_ints = [
            0.,
            40.02416012,
            17.23521259,
            57.65661496
        ]
        gaussian_freqs = [
            0.,
            3685.815,
            2706.132,
            1383.391,
            7202.835,
            5323.917,
            2749.027,
            6377.958,
            5044.721,
            4072.407
        ]
        gaussian_ints = [
            0.,
            36.93616327,
            14.01377595,
            58.17413771,
            1.16678234,
            0.55376888,
            3.30245090,
            0.19253704,
            1.57560144,
            1.82673531
        ]

        print_specs = True
        if print_specs:
            report = "Harmonic:   State     Freq.   Int.    Gaussian: Freq.   Int.\n" + "\n".join(
                " " * 12 + "{} {:>7.2f} {:>7.4f}           {:>7.2f} {:>7.4f}".format(s, f, i, gf, g)
                for s, f, i, gf, g in zip(states, harm_freqs, harm_ints, gaussian_harm_freqs, gaussian_harm_ints)
            )
            print(report)
            report = "Anharmonic: State     Freq.   Int.    Gaussian: Freq.   Int.\n" + "\n".join(
                " " * 12 + "{} {:>7.2f} {:>7.4f}           {:>7.2f} {:>7.4f}".format(s, f, i, gf, g)
                for s, f, i, gf, g in zip(states, freqs, ints, gaussian_freqs, gaussian_ints)
            )
            print(report)

        self.assertEquals(
            [round(x, 2) for x in gaussian_harm_ints[1:]],
            list(np.round(harm_ints[1:4], 2))
        )
        self.assertEquals(
            [round(x, 2) for x in gaussian_ints[1:]],
            list(np.round(ints[1:10], 2))
        )

    #endregion Test Intensities

    #region Intensity Breakdowns

    @validationTest
    def test_HODCartesianIntensityBreakdown(self):

        internals = [
            [0, -1, -1, -1],
            [1,  0, -1, -1],
            [2,  0,  1, -1]
        ]
        states = (
            (0, 0, 0),
            (0, 0, 1), (0, 1, 0), (1, 0, 0),
            (0, 0, 2), (0, 2, 0), (2, 0, 0),
            (0, 1, 1), (1, 0, 1), (1, 1, 0)
        )
        coupled_states = self.get_states(5, 3, max_quanta=5)

        all_wfns = self.get_VPT2_wfns(
            "HOD_freq.fchk",
            internals,
            states,
            regenerate=True,
            coupled_states=coupled_states,
            get_breakdown=True
        )

        import io, csv
        s = io.StringIO()
        # with open('/Users/Mark/Desktop/HOD_freq.csv', "w+") as s:
        writer = csv.writer(s)
        for k, wnf in all_wfns.items():
            writer.writerow([k + " Wavefunctions"])
            wnf.write_CSV_breakdown(s, wnf.generate_intensity_breakdown(), padding=[" "] * 2)
        print(s.getvalue())

    @inactiveTest
    def test_HODIntensityBreakdown(self):

        internals = [
            [0, -1, -1, -1],
            [1,  0, -1, -1],
            [2,  0,  1, -1]
        ]
        states = (
            (0, 0, 0),
            (0, 0, 1), (0, 1, 0), (1, 0, 0),
            (0, 0, 2), (0, 2, 0), (2, 0, 0),
            (0, 1, 1), (1, 0, 1), (1, 1, 0)
        )
        coupled_states = self.get_states(5, 3, max_quanta=5)

        all_wfns = self.get_VPT2_wfns(
            "HOD_freq.fchk",
            internals,
            states,
            regenerate=True,
            coupled_states=coupled_states,
            get_breakdown=True
        )

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        # trying to turn off the orthogonality condition
        # states = wfns.corrs.states
        # for i,s in enumerate(states):
        #     wfns.corrs.wfn_corrections[i, 2, s] = 0 # turn off second correction

        for mode in ['standard', 'intuitive']:

            terms = []
            plot_specs = False
            if plot_specs:
                import McUtils.Plots as plt
                g = plt.GraphicsGrid(nrows=len(all_wfns), ncols=4, subimage_size=(200, 100))
            for i, wfns in enumerate(all_wfns):
                dts = wfns.dipole_terms
                wfns.dipole_partitioning=mode
                wfn_terms=[]

                dipole_breakdowns = [
                    [[d[0], d[1], np.zeros(d[2].shape), np.zeros(d[3].shape)] for d in dts],
                    [[d[0], d[1], d[2], np.zeros(d[3].shape)] for d in dts],
                    [[d[0], d[1], np.zeros(d[2].shape), d[3]] for d in dts],
                    dts
                ]

                engs = h2w * wfns.energies
                freqs = engs - engs[0]

                harm_engs = h2w * wfns.zero_order_energies
                harm_freqs = harm_engs - harm_engs[0]

                wfn_terms.append(freqs.tolist())
                for j, dt in enumerate(dipole_breakdowns):
                    wfns.dipole_terms = dt
                    ints = wfns.intensities
                    harm_ints = wfns.zero_order_intensities

                    wfn_terms.append(ints.tolist())

                    if plot_specs:
                        s = plt.StickPlot(freqs, ints,
                                          aspect_ratio=.5,
                                          plot_legend=("Anharmonic", "Harmonic"),
                                          image_size=500,
                                          figure=g[i, j]
                                          )
                        plt.StickPlot(harm_freqs, harm_ints, figure=s, plot_style=dict(linefmt='red'),
                                      aspect_ratio=.5,
                                      axes_labels=["$\nu$ (cm$^{-1}$)", "I (km mol$^{-1}$)"],
                                      plot_legend=("Anharmonic", "Harmonic")
                                      )
                    corrs = wfns.transition_moment_corrections
                    wfn_terms.append([[[[list(a) for a in z] for z in y] for y in x] for x in corrs])

                wfn_terms.append(wfns.corrs.wfn_corrections[:, :, wfns.corrs.coupled_states].tolist())

                terms.append(wfn_terms)
            if plot_specs:
                g.show()


            import json

            if wfns.dipole_partitioning == 'intuitive':
                dump_file = "/Users/Mark/Documents/UW/Research/WaterPT/anne_HOD_intuitive.json"
            else:
                dump_file = "/Users/Mark/Documents/UW/Research/WaterPT/anne_HOD.json"
            with open(dump_file, "w+") as f:
                json.dump(terms, f)

    @validationTest
    def test_HODCartesianIntensityBreakdown(self):

        internals = None
        states = (
            (0, 0, 0),
            (0, 0, 1), (0, 1, 0), (1, 0, 0),
            (0, 0, 2), (0, 2, 0), (2, 0, 0),
            (0, 1, 1), (1, 0, 1), (1, 1, 0)
        )
        coupled_states = self.get_states(5, 3, max_quanta=5)

        all_wfns = self.get_VPT2_wfns(
            "HOD_freq.fchk",
            internals,
            states,
            regenerate=True,
            coupled_states=coupled_states,
            get_breakdown=True
        )

        import io, csv
        s = io.StringIO()
        # with open('/Users/Mark/Desktop/HOD_freq.csv', "w+") as s:
        writer = csv.writer(s)
        for k,wnf in all_wfns.items():
            writer.writerow([k + " Wavefunctions"])
            wnf.write_CSV_breakdown(s, wnf.generate_intensity_breakdown(), padding=[" "]*2)
        print(s.getvalue())

        # h2w = UnitsData.convert("Hartrees", "Wavenumbers")
        #
        # # trying to turn off the orthogonality condition
        # # states = wfns.corrs.states
        # # for i,s in enumerate(states):
        # #     wfns.corrs.wfn_corrections[i, 2, s] = 0 # turn off second correction
        #
        #
        # for mode in ['standard', 'intuitive']:
        #
        #     terms = []
        #     plot_specs = False
        #     if plot_specs:
        #         import McUtils.Plots as plt
        #         g = plt.GraphicsGrid(nrows=len(all_wfns), ncols=4, subimage_size=(200, 100))
        #
        #     for i, wfns in enumerate(all_wfns):
        #         dts = wfns.dipole_terms
        #         wfns.dipole_partitioning=mode
        #         wfn_terms = []
        #
        #         dipole_breakdowns = [
        #             [[d[0], d[1], np.zeros(d[2].shape), np.zeros(d[3].shape)] for d in dts],
        #             [[d[0], d[1], d[2], np.zeros(d[3].shape)] for d in dts],
        #             [[d[0], d[1], np.zeros(d[2].shape), d[3]] for d in dts],
        #             dts
        #         ]
        #
        #         engs = h2w * wfns.energies
        #         freqs = engs - engs[0]
        #
        #         harm_engs = h2w * wfns.zero_order_energies
        #         harm_freqs = harm_engs - harm_engs[0]
        #
        #         wfn_terms.append(freqs.tolist())
        #         for j, dt in enumerate(dipole_breakdowns):
        #             wfns.dipole_terms = dt
        #             ints = wfns.intensities
        #             harm_ints = wfns.zero_order_intensities
        #
        #             wfn_terms.append(ints.tolist())
        #
        #             if plot_specs:
        #                 s = plt.StickPlot(freqs, ints,
        #                                   aspect_ratio=.5,
        #                                   plot_legend=("Anharmonic", "Harmonic"),
        #                                   image_size=500,
        #                                   figure=g[i, j]
        #                                   )
        #                 plt.StickPlot(harm_freqs, harm_ints, figure=s, plot_style=dict(linefmt='red'),
        #                               aspect_ratio=.5,
        #                               axes_labels=["$\nu$ (cm$^{-1}$)", "I (km mol$^{-1}$)"],
        #                               plot_legend=("Anharmonic", "Harmonic")
        #                               )
        #
        #             corrs = wfns.transition_moment_corrections
        #             wfn_terms.append([[[[list(a) for a in z] for z in y] for y in x] for x in corrs])
        #
        #         wfn_terms.append(wfns.corrs.wfn_corrections[:, :, wfns.corrs.coupled_states].tolist())
        #
        #         terms.append(wfn_terms)
        #     if plot_specs:
        #         g.show()
        #
        #     import json
        #
        #     if wfns.dipole_partitioning == 'intuitive':
        #         dump_file = "/Users/Mark/Documents/UW/Research/WaterPT/anne_HOD_cart_intuitive.json"
        #     else:
        #         dump_file = "/Users/Mark/Documents/UW/Research/WaterPT/anne_HOD_cart.json"
        #     with open(dump_file, "w+") as f:
        #         json.dump(terms, f)

    @inactiveTest
    def test_HOHIntensityBreakdown(self):

        internals = [
            [0, -1, -1, -1],
            [1,  0, -1, -1],
            [2,  0,  1, -1]
        ]
        states = (
            (0, 0, 0),
            (0, 0, 1), (0, 1, 0), (1, 0, 0),
            (0, 0, 2), (0, 2, 0), (2, 0, 0),
            (0, 1, 1), (1, 0, 1), (1, 1, 0)
        )
        coupled_states = self.get_states(5, 3, max_quanta=5)

        all_wfns = self.get_VPT2_wfns(
            "HOH_freq.fchk",
            internals,
            states,
            regenerate=True,
            coupled_states=coupled_states,
            get_breakdown=True
        )

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        # trying to turn off the orthogonality condition
        # states = wfns.corrs.states
        # for i,s in enumerate(states):
        #     wfns.corrs.wfn_corrections[i, 2, s] = 0 # turn off second correction

        for mode in ['standard', 'intuitive']:

            terms = []
            plot_specs = False
            if plot_specs:
                import McUtils.Plots as plt
                g = plt.GraphicsGrid(nrows=len(all_wfns), ncols=4, subimage_size=(200, 100))

            for i, wfns in enumerate(all_wfns):
                dts = wfns.dipole_terms
                wfns.dipole_partitioning = mode
                wfn_terms = []

                dipole_breakdowns = [
                    [[d[0], d[1], np.zeros(d[2].shape), np.zeros(d[3].shape)] for d in dts],
                    [[d[0], d[1], d[2], np.zeros(d[3].shape)] for d in dts],
                    [[d[0], d[1], np.zeros(d[2].shape), d[3]] for d in dts],
                    dts
                ]

                engs = h2w * wfns.energies
                freqs = engs - engs[0]

                harm_engs = h2w * wfns.zero_order_energies
                harm_freqs = harm_engs - harm_engs[0]

                wfn_terms.append(freqs.tolist())
                for j, dt in enumerate(dipole_breakdowns):
                    wfns.dipole_terms = dt
                    ints = wfns.intensities
                    harm_ints = wfns.zero_order_intensities

                    wfn_terms.append(ints.tolist())

                    if plot_specs:
                        s = plt.StickPlot(freqs, ints,
                                          aspect_ratio=.5,
                                          plot_legend=("Anharmonic", "Harmonic"),
                                          image_size=500,
                                          figure=g[i, j]
                                          )
                        plt.StickPlot(harm_freqs, harm_ints, figure=s, plot_style=dict(linefmt='red'),
                                      aspect_ratio=.5,
                                      axes_labels=["$\nu$ (cm$^{-1}$)", "I (km mol$^{-1}$)"],
                                      plot_legend=("Anharmonic", "Harmonic")
                                      )
                    corrs = wfns.transition_moment_corrections
                    wfn_terms.append([[[[list(a) for a in z] for z in y] for y in x] for x in corrs])

                wfn_terms.append(wfns.corrs.wfn_corrections[:, :, wfns.corrs.coupled_states].tolist())

                terms.append(wfn_terms)
            if plot_specs:
                g.show()

            import json

            if wfns.dipole_partitioning == 'intuitive':
                dump_file = "/Users/Mark/Documents/UW/Research/WaterPT/anne_HOH_intuitive.json"
            else:
                dump_file = "/Users/Mark/Documents/UW/Research/WaterPT/anne_HOH.json"
            with open(dump_file, "w+") as f:
                json.dump(terms, f)

    @inactiveTest
    def test_HOHCartesianIntensityBreakdown(self):

        internals = None
        states = (
            (0, 0, 0),
            (0, 0, 1), (0, 1, 0), (1, 0, 0),
            (0, 0, 2), (0, 2, 0), (2, 0, 0),
            (0, 1, 1), (1, 0, 1), (1, 1, 0)
        )
        coupled_states = self.get_states(5, 3, max_quanta=5)

        all_wfns = self.get_VPT2_wfns(
            "HOH_freq.fchk",
            internals,
            states,
            regenerate=True,
            coupled_states=coupled_states,
            get_breakdown=True
        )

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        # trying to turn off the orthogonality condition
        # states = wfns.corrs.states
        # for i,s in enumerate(states):
        #     wfns.corrs.wfn_corrections[i, 2, s] = 0 # turn off second correction

        for mode in ['standard', 'intuitive']:

            terms = []
            plot_specs = False
            if plot_specs:
                import McUtils.Plots as plt
                g = plt.GraphicsGrid(nrows=len(all_wfns), ncols=4, subimage_size=(200, 100))

            for i, wfns in enumerate(all_wfns):
                dts = wfns.dipole_terms
                wfns.dipole_partitioning = mode
                wfn_terms = []

                dipole_breakdowns = [
                    [[d[0], d[1], np.zeros(d[2].shape), np.zeros(d[3].shape)] for d in dts],
                    [[d[0], d[1], d[2], np.zeros(d[3].shape)] for d in dts],
                    [[d[0], d[1], np.zeros(d[2].shape), d[3]] for d in dts],
                    dts
                ]

                engs = h2w * wfns.energies
                freqs = engs - engs[0]

                harm_engs = h2w * wfns.zero_order_energies
                harm_freqs = harm_engs - harm_engs[0]

                wfn_terms.append(freqs.tolist())
                for j, dt in enumerate(dipole_breakdowns):
                    wfns.dipole_terms = dt
                    ints = wfns.intensities
                    harm_ints = wfns.zero_order_intensities

                    wfn_terms.append(ints.tolist())

                    if plot_specs:
                        s = plt.StickPlot(freqs, ints,
                                          aspect_ratio=.5,
                                          plot_legend=("Anharmonic", "Harmonic"),
                                          image_size=500,
                                          figure=g[i, j]
                                          )
                        plt.StickPlot(harm_freqs, harm_ints, figure=s, plot_style=dict(linefmt='red'),
                                      aspect_ratio=.5,
                                      axes_labels=["$\nu$ (cm$^{-1}$)", "I (km mol$^{-1}$)"],
                                      plot_legend=("Anharmonic", "Harmonic")
                                      )
                    corrs = wfns.transition_moment_corrections
                    wfn_terms.append([[[[list(a) for a in z] for z in y] for y in x] for x in corrs])

                wfn_terms.append(wfns.corrs.wfn_corrections[:, :, wfns.corrs.coupled_states].tolist())

                terms.append(wfn_terms)

            if plot_specs:
                g.show()

            import json

            if wfns.dipole_partitioning == 'intuitive':
                dump_file = "/Users/Mark/Documents/UW/Research/WaterPT/anne_HOH_cart_intuitive.json"
            else:
                dump_file = "/Users/Mark/Documents/UW/Research/WaterPT/anne_HOH_cart.json"
            with open(dump_file, "w+") as f:
                json.dump(terms, f)

    #endregion Intensity Breakdowns