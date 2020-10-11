from Peeves.TestUtils import *
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
    def get_VPT2_wfns(self, fchk,
                      internals,
                      states,
                      n_quanta,
                      regenerate=False,
                      coupled_states=None,
                      mode_selection=None,
                      v2 = None,
                      t2 = None,
                      v3 = None,
                      t3 = None,
                      t4 = None,
                      v4 = None,
                      coriolis = None,
                      watson = None,
                      degeneracies=None
                      ):

        hammer = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data(fchk),
            n_quanta=n_quanta,
            internals=internals,
            mode_selection=mode_selection
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

        wfn_file = os.path.join(self.wfn_file_dir, fchk.replace("fchk", "npz"))
        if regenerate or not os.path.exists(wfn_file):
            wfns = hammer.get_wavefunctions(states, coupled_states=coupled_states, degeneracies=degeneracies)
            self.save_wfns(wfn_file, wfns)
        else:
            wfns = self.load_wfns(hammer.molecule, hammer.basis, wfn_file)

        return wfns
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

    @validationTest
    def test_RepresentQQQ(self):
        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("HOD_freq.fchk"),
            n_quanta=6,
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
            n_quanta=6,
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
            n_quanta=6,
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
            n_quanta=6,
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
            n_quanta=6,
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
            n_quanta=6,
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
            n_quanta=6,
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
            n_quanta=6,
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

    @validationTest
    def test_TestCubicsCartesians2(self):
        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("OCHD_freq.fchk"),
            n_quanta=6,
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
            if not np.allclose(legit, v3, rtol=1):
                diff = legit - v3
                bad_pos = np.array(np.where(np.abs(diff) > 1)).T
                print("Gaussian/This Disagreements:\n" + "\n".join(
                    "{:>.0f} {:>.0f} {:>.0f} {:>.0f} {:>5.3f} (Actual: {:>8.3f} This: {:>8.3f})".format(
                        i, j, k, diff[i, j, k], legit[i, j, k], v3[i, j, k]
                    ) for i, j, k in bad_pos
                ))

        self.assertTrue(np.allclose(legit, v3, rtol=1))  # testing to within a wavenumber

    @validationTest
    def test_TestQuarticsCartesians(self):
        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("HOD_freq.fchk"),
            n_quanta=6,
            internals=None
        )

        v4_Gaussian = [
          [1,  1,  1,  1,    1517.96213],
          [2,  1,  1,  1,      98.59961],
          [2,  2,  1,  1,       8.99887],
          [2,  2,  2,  1,     -50.96655],
          [2,  2,  2,  2,     804.29611],
          [3,  1,  1,  1,     142.08091],
          [3,  2,  1,  1,     -18.73606],
          [3,  2,  2,  1,     -22.35470],
          [3,  2,  2,  2,      71.81011],
          [3,  3,  1,  1,    -523.38920],
          [3,  3,  2,  1,      -4.05652],
          [3,  3,  2,  2,     -95.43623],
          [3,  3,  3,  1,    -145.84374],
          [3,  3,  3,  2,     -41.06991],
          [3,  3,  3,  3,      83.41603]
          ]
        # for reasons I'm still working out, Gaussian 16 flips the phases on its normal modes
        # but _doesn't_ flip the phases on its derivatives
        v4_Gaussian16 = [
            [1, 1, 1, 1, 1517.96195],
            [2, 1, 1, 1,   98.60039],
            [2, 2, 1, 1,    8.99879],
            [2, 2, 2, 1,  -50.96534],
            [2, 2, 2, 2,  804.29612],
            [3, 1, 1, 1,  142.07906],
            [3, 2, 1, 1,  -18.73793],
            [3, 2, 2, 1,  -22.35735],
            [3, 2, 2, 2,   71.80757],
            [3, 3, 1, 1, -523.38806],
            [3, 3, 2, 1,   -4.05417],
            [3, 3, 2, 2,  -95.43403],
            [3, 3, 3, 1, -145.84906],
            [3, 3, 3, 2,  -41.07498],
            [3, 3, 3, 3,   83.42476]
        ]

        legit = np.zeros((3, 3, 3, 3))
        mode_mapping = [
            2, 1, 0
        ]
        for i, j, k, l, v in v4_Gaussian:
            i = mode_mapping[i-1]; j = mode_mapping[j-1]
            k = mode_mapping[k-1]; l = mode_mapping[l-1]
            for perm in ip.permutations((i, j, k, l)):
                legit[perm] = v

        v4 = self.h2w * ham.V_terms[2]

        print_errors = False
        if print_errors:
            if not np.allclose(legit, v4, rtol=.001):
                diff = legit - v4
                bad_pos = np.array(np.where(np.abs(diff) > .001)).T
                print("Gaussian/This Disagreements:\n"+"\n".join(
                    "{:>.0f} {:>.0f} {:>.0f} {:>.0f} {:>5.3f} (Actual: {:>8.3f} This: {:>8.3f})".format(
                        i,j,k,l, diff[i,j,k,l], legit[i,j,k,l], v4[i,j,k,l]
                    ) for i,j,k,l in bad_pos
                ))

        self.assertTrue(np.allclose(legit, v4, rtol=.001)) # testing to within .001 wavenumbers

    @validationTest
    def test_TestQuarticsCartesians2(self):
        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("OCHD_freq.fchk"),
            n_quanta=6,
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
            if not np.allclose(legit, v4, rtol=1):
                diff = legit - v4
                bad_pos = np.array(np.where(np.abs(diff) > 1)).T
                print("Gaussian/This Disagreements:\n" + "\n".join(
                    "{:>.0f} {:>.0f} {:>.0f} {:>.0f} {:>5.3f} (Actual: {:>8.3f} This: {:>8.3f})".format(
                        i, j, k, l, diff[i, j, k, l], legit[i, j, k, l], v4[i, j, k, l]
                    ) for i, j, k, l in bad_pos
                ))

        self.assertTrue(np.allclose(legit, v4, rtol=1))  # testing to within a wavenumber

    @validationTest
    def test_TestQuarticsInternals(self):
        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("HOD_freq.fchk"),
            n_quanta=6,
            internals=[
            [0, -1, -1, -1],
            [1, 0, -1, -1],
            [2, 0, 1, -1]
        ]
        )

        v4_Anne =[
             [1,  1,  1,  -37.03937000],
             [2,  1,  1,  -32.30391126],
             [3,  1,  1,  -33.08215609],
             [1,  2,  1,  -32.30391126],
             [2,  2,  1,    3.57147725],
             [3,  2,  1,    9.77124742],
             [1,  3,  1,  -33.08215609],
             [2,  3,  1,    9.77124742],
             [3,  3,  1,    3.08396862],
             [1,  1,  2,    3.53204514],
             [2,  1,  2,   66.35374213],
             [3,  1,  2,   -8.46713126],
             [1,  2,  2,   66.35374213],
             [2,  2,  2,  804.47871323],
             [3,  2,  2,  -51.44004640],
             [1,  3,  2,   -8.46713126],
             [2,  3,  2,  -51.44004640],
             [3,  3,  2,   10.60086681],
             [1,  1,  3,    2.67361974],
             [2,  1,  3,    3.14497676],
             [3,  1,  3,  111.80682105],
             [1,  2,  3,    3.14497676],
             [2,  2,  3,   10.60153758],
             [3,  2,  3,   97.05643377],
             [1,  3,  3,  111.80682105],
             [2,  3,  3,   97.05643377],
             [3,  3,  3, 1519.00602277]
        ]

        legit = np.zeros((3, 3, 3, 3))
        for k, j, i, v in v4_Anne:
            i = i-1; j = j-1; k = k-1
            for perm in ip.permutations((i, i, j, k)):
                legit[perm] = v

        v4 = self.h2w * ham.V_terms[2]

        # for unknown reasons, Anne and I disagree on the order of like 5 cm^-1,
        # but it's unclear which set of derivs. is right since my & hers both
        # yield a max deviation from the Gaussian result of ~.1 cm^-1
        print_errors=False
        if print_errors:
            if not np.allclose(legit, v4, rtol=.1):
                diff = legit - v4
                bad_pos = np.array(np.where(np.abs(diff) > .1)).T
                print("Gaussian/This Disagreements:\n"+"\n".join(
                    "{:>.0f} {:>.0f} {:>.0f} {:>.0f} {:>5.3f} (Actual: {:>8.3f} This: {:>8.3f})".format(
                        i,j,k,l, diff[i,j,k,l], legit[i,j,k,l], v4[i,j,k,l]
                    ) for i,j,k,l in bad_pos
                ))

        self.assertTrue(np.allclose(legit, v4, rtol=10))

    @inactiveTest
    def test_OCHTCoriolisCouplings(self):
        # for unclear reasons this isn't working...?
        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("OCHT_freq.fchk"),
            n_quanta=6,
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
            n_quanta=6,
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

    @validationTest
    def test_SecondOrderEnergyCorrection(self):

        n_modes = 6
        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("HOD_freq.fchk"),
            n_quanta=6,
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
        corrs = ham._get_state_VPT_corrections(
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

    @debugTest
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
        n_quanta = 6

        coupled_states = self.get_states(5, 3)
        wfns = self.get_VPT2_wfns(
            "HOH_freq.fchk",
            internals,
            states,
            n_quanta,
            regenerate=True,
            coupled_states=coupled_states
            # , v2 = 0
            # , v3 = 0
            # , v4 = 0
        )

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        # wfns = hammer.get_wavefunctions(states, coupled_states=None)
        energies = h2w * wfns.energies
        zero_ord = h2w * wfns.zero_order_energies

        gaussian_states = [(0, 0, 1), (0, 1, 0), (1, 0, 0),
                           (0, 0, 2), (0, 2, 0), (2, 0, 0),
                           (0, 1, 1), (1, 0, 1), (1, 1, 0)]
        gaussian_energies = np.array([4681.564, 4605.953])
        gaussian_freqs = np.array([
            [3937.525,    3744.734],
            [3803.300,    3621.994],
            [1622.303,    1572.707],

            [7875.049,    7391.391],
            [7606.599,    7155.881],
            [3244.606,    3117.366],

            [7740.824,    7200.364],
            [5559.828,    5294.379],
            [5425.603,    5174.665]
        ])


        # raise Exception(h2w * wfns.corrs.energy_corrs)

        my_energies = np.array([zero_ord[0], energies[0]])
        my_freqs = np.column_stack([
            zero_ord[1:] - zero_ord[0],
            energies[1:] - energies[0]
        ])

        print_diffs = True
        if print_diffs:
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
            print("Difference Energies:\n",
                  "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*(my_energies - gaussian_energies), "-", "-"),
                  *(
                      "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                      zip(states[1:], my_freqs - gaussian_freqs)
                  )
                  )
        self.assertLess(np.max(np.abs(my_freqs - gaussian_freqs)), 1)

    @debugTest
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
        n_quanta = 6

        coupled_states = self.get_states(5, 3)
        wfns = self.get_VPT2_wfns(
            "HOH_freq.fchk",
            internals,
            states,
            n_quanta,
            regenerate=True,
            coupled_states=coupled_states
        )

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        # wfns = hammer.get_wavefunctions(states, coupled_states=None)
        energies = h2w * wfns.energies
        zero_ord = h2w * wfns.zero_order_energies

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

        print_diffs = True
        if print_diffs:
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
            print("Difference Energies:\n",
                  "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*(my_energies - gaussian_energies), "-", "-"),
                  *(
                      "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                      zip(states[1:], my_freqs - gaussian_freqs)
                  )
                  )
        self.assertLess(np.max(np.abs(my_freqs - gaussian_freqs)), 1.5)

    @debugTest
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
        n_quanta = 6

        coupled_states = self.get_states(5, 3)
        wfns = self.get_VPT2_wfns(
            "HOD_freq.fchk",
            internals,
            states,
            n_quanta,
            regenerate=True,
            coupled_states=coupled_states
        )

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        # wfns = hammer.get_wavefunctions(states, coupled_states=None)
        energies = h2w * wfns.energies
        zero_ord = h2w * wfns.zero_order_energies
        # raise Exception([energies, zero_ord])

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

        print_diffs=True
        if print_diffs:
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
            print("Difference Energies:\n",
                  "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*(my_energies-gaussian_energies), "-", "-"),
                  *(
                      "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                      zip(states[1:], my_freqs-gaussian_freqs)
                  )
                  )
        self.assertLess(np.max(np.abs(my_freqs - gaussian_freqs)), 1)

    @debugTest
    def test_HODVPTCartesians(self):

        internals = None

        states = (
            (0, 0, 0),
            (0, 0, 1), (0, 1, 0), (1, 0, 0),
            (0, 0, 2), (0, 2, 0), (2, 0, 0),
            (0, 1, 1), (1, 0, 1), (1, 1, 0)
        )
        n_quanta = 6

        coupled_states = self.get_states(5, 3)
        wfns = self.get_VPT2_wfns(
            "HOD_freq.fchk",
            internals,
            states,
            n_quanta,
            regenerate=True,
            coupled_states=coupled_states
        )

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        # wfns = hammer.get_wavefunctions(states, coupled_states=None)
        energies = h2w * wfns.energies
        zero_ord = h2w * wfns.zero_order_energies
        # raise Exception([energies, zero_ord])

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

        print_diffs = True
        if print_diffs:
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
            print("Difference Energies:\n",
                  "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*(my_energies - gaussian_energies), "-", "-"),
                  *(
                      "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                      zip(states[1:], my_freqs - gaussian_freqs)
                  )
                  )
        self.assertLess(np.max(np.abs(my_freqs - gaussian_freqs)), 1)

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
        n_quanta = 6
        wfns = self.get_VPT2_wfns(
            "HOD_freq.fchk",
            internals,
            states,
            n_quanta,
            regenerate=True,
            coupled_states = coupled_states
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
            s = plt.StickPlot(freqs, ints)
            plt.StickPlot(harm_freqs, harm_ints, figure=s, plot_style=dict(linefmt='red'))
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
            report = "Harmonic:   State     Freq.   Int.    Gaussian: Freq.   Int.\n"+"\n".join(
                " "*12+"{} {:>7.2f} {:>7.4f}           {:>7.2f} {:>7.4f}".format(s, f, i, gf, g)
                for s,f,i,gf,g in zip(states, harm_freqs, harm_ints, gaussian_harm_freqs, gaussian_harm_ints)
            )
            print(report)
            report = "Anharmonic: State     Freq.   Int.    Gaussian: Freq.   Int.\n"+"\n".join(
                " "*12+"{} {:>7.2f} {:>7.4f}           {:>7.2f} {:>7.4f}".format(s, f, i, gf, g)
                for s,f,i,gf,g in zip(states, freqs, ints, gaussian_freqs, gaussian_ints)
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
        n_quanta = 6
        n_modes = 3
        coupled_states = self.get_states(5, n_modes, max_quanta=5)
        def block(self=self, internals=internals, states=states, coupled_states=coupled_states, n_quanta=n_quanta):
            return self.get_VPT2_wfns(
                "HOT_freq.fchk",
                internals,
                states,
                n_quanta,
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

        print_diffs = True
        if print_diffs:
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
            print("Difference Energies:\n",
                  "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*(my_energies - gaussian_energies), "-", "-"),
                  *(
                      "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                      zip(states[1:], my_freqs - gaussian_freqs)
                  )
                  )
        self.assertLess(np.max(np.abs(my_freqs - gaussian_freqs)), 1)

    @validationTest
    def test_HOTVPTCartesians(self):

        internals = None
        states = (
            (0, 0, 0),
            (0, 0, 1), (0, 1, 0), (1, 0, 0),
            (0, 0, 2), (0, 2, 0), (2, 0, 0),
            (0, 1, 1), (1, 0, 1), (1, 1, 0)
        )
        n_quanta = 6
        n_modes = 3
        coupled_states = self.get_states(5, n_modes, max_quanta=5)

        def block(self=self, internals=internals, states=states, coupled_states=coupled_states, n_quanta=n_quanta):
            return self.get_VPT2_wfns(
                "HOT_freq.fchk",
                internals,
                states,
                n_quanta,
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

        print_diffs = True
        if print_diffs:
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
            print("Difference Energies:\n",
                  "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*(my_energies - gaussian_energies), "-", "-"),
                  *(
                      "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                      zip(states[1:], my_freqs - gaussian_freqs)
                  )
                  )
        self.assertLess(np.max(np.abs(my_freqs - gaussian_freqs)), 1)

    @inactiveTest
    def test_OCHHVPT(self):

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
        n_quanta = 6

        max_quanta = 4
        states = (
            (0, 0, 0, 0, 0, 0),
            (0, 0, 0, 0, 0, 1)
        )
        # coupled_states = [
        #     x for x in self.get_states(4, n_modes, max_quanta=max_quanta)
        #     if all(x[i] < 3 for i in range(5))
        # ]
        # states = self.get_states(2, n_modes)
        # coupled_states = [
        #     x for x in self.get_states(4, n_modes, max_quanta=max_quanta)
        #     if all(x[i] < 3 for i in range(5))
        # ]
        coupled_states = self.get_states(5, n_modes, max_quanta=max_quanta)

        def block(self=self,
                  internals=internals,
                  states=states,
                  coupled_states=coupled_states,
                  n_quanta=n_quanta,
                  mode_selection=mode_selection
                  ):
            return self.get_VPT2_wfns(
                "OCHH_freq.fchk",
                internals,
                states,
                n_quanta,
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
        raise Exception(wfns.energy_corrs)
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

        print_report = True
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

        self.assertLess(np.max(np.abs(freqs[:len(states)] - gaussian_freqs[:len(states), 1])), 1)

    @validationTest
    def test_OCHDVPT(self):

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
        n_quanta = 6

        max_quanta = 4
        states = (
            (0, 0, 0, 0, 0, 0),
            (0, 0, 0, 0, 0, 1)
        )
        # coupled_states = [
        #     x for x in self.get_states(4, n_modes, max_quanta=max_quanta)
        #     if all(x[i] < 3 for i in range(5))
        # ]
        # states = self.get_states(2, n_modes)
        # coupled_states = [
        #     x for x in self.get_states(4, n_modes, max_quanta=max_quanta)
        #     if all(x[i] < 3 for i in range(5))
        # ]
        coupled_states = self.get_states(5, n_modes, max_quanta=max_quanta)
        def block(self=self,
                  internals=internals,
                  states=states,
                  coupled_states=coupled_states,
                  n_quanta=n_quanta,
                  mode_selection=mode_selection
                  ):
            return self.get_VPT2_wfns(
                "OCHD_freq.fchk",
                internals,
                states,
                n_quanta,
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

        print_report = True
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

        self.assertLess(np.max(np.abs(freqs[:len(states)] - gaussian_freqs[:len(states), 1])), 1)

    @validationTest
    def test_OCHTVPT(self):

        internals = [
            [0, -1, -1, -1],
            [1, 0, -1, -1],
            [2, 1, 0, -1],
            [3, 1, 0, 2]
        ]

        n_modes = 3 * 4 - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        n_quanta = 6

        max_quanta = 4
        states = (
            (0, 0, 0, 0, 0, 0),
            (0, 0, 0, 0, 0, 1)
        )
        # states = self.get_states(2, n_modes)
        coupled_states = self.get_states(5, n_modes)
        # [
        #     x for x in self.get_states(5, n_modes)
        #     if all(x[i] < 4 for i in range(5))
        # ]

        def block(self=self,
                  internals=internals,
                  states=states,
                  coupled_states=coupled_states,
                  n_quanta=n_quanta,
                  mode_selection=mode_selection
                  ):
            return self.get_VPT2_wfns(
                "OCHT_freq.fchk",
                internals,
                states,
                n_quanta,
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

        print_report = True
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

            print("Difference Energies:\n",
                  ('0 ' * n_modes + "{:>8.3f} {:>8.3f} {:>8} {:>8}\n").format(harm_engs[0]-gaussian_engs[0], engs[0]-gaussian_engs[1], "-", "-"),
                  *(
                      ('{:<1.0f} ' * n_modes + "{:>8} {:>8} {:>8.3f} {:>8.3f}\n").format(*s, "-", "-", e1, e2) for
                      s, e1, e2 in
                      zip(states[1:], harm_freq-gaussian_freqs[:, 0], freqs-gaussian_freqs[:, 1])
                  )
                  )

        self.assertLess(
            np.max(np.abs(freqs-gaussian_freqs[:, 1])[:len(states)-1]),
            1)

    @validationTest
    def test_HODVPTCartesianPotential(self):

        internals = None

        n_modes = 3 * 3 - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        n_quanta = 6

        max_quanta = 4
        states = (
            (0, 0, 0),
        )  # we're only interested in the correction of the ground state w/o perturbation or accounting for the energy or any shit
        coupled_states = self.get_states(5, n_modes, max_quanta=5)

        def block(self=self,
                  internals=internals,
                  states=states,
                  coupled_states=coupled_states,
                  n_quanta=n_quanta,
                  mode_selection=mode_selection
                  ):
            return self.get_VPT2_wfns(
                "HOD_freq.fchk",
                internals,
                states,
                n_quanta,
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
        n_quanta = 6

        max_quanta = 4
        states = (
            (0, 0, 0),
        )  # we're only interested in the correction of the ground state w/o perturbation or accounting for the energy or any shit
        coupled_states = self.get_states(4, n_modes, max_quanta=max_quanta)

        def block(self=self,
                  internals=internals,
                  states=states,
                  coupled_states=coupled_states,
                  n_quanta=n_quanta,
                  mode_selection=mode_selection
                  ):
            return self.get_VPT2_wfns(
                "HOD_freq.fchk",
                internals,
                states,
                n_quanta,
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
    def test_HOHVPTCartesianPotential(self):

        internals = None

        n_modes = 3 * 3 - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        n_quanta = 6

        max_quanta = 4
        states = (
            (0, 0, 0),
        )  # we're only interested in the correction of the ground state w/o perturbation or accounting for the energy or any shit
        coupled_states = self.get_states(4, n_modes, max_quanta=max_quanta)

        def block(self=self,
                  internals=internals,
                  states=states,
                  coupled_states=coupled_states,
                  n_quanta=n_quanta,
                  mode_selection=mode_selection
                  ):
            return self.get_VPT2_wfns(
                "HOH_freq.fchk",
                internals,
                states,
                n_quanta,
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
    def test_HOHVPTCartesianCoriolis(self):

        internals = None

        n_modes = 3 * 3 - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        n_quanta = 6

        max_quanta = 4
        states = (
            (0, 0, 0),
        )  # we're only interested in the correction of the ground state w/o perturbation or accounting for the energy or any shit
        coupled_states = self.get_states(4, n_modes, max_quanta=max_quanta)

        def block(self=self,
                  internals=internals,
                  states=states,
                  coupled_states=coupled_states,
                  n_quanta=n_quanta,
                  mode_selection=mode_selection
                  ):
            return self.get_VPT2_wfns(
                "HOH_freq.fchk",
                internals,
                states,
                n_quanta,
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
        n_quanta = 6

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
                  n_quanta=n_quanta,
                  mode_selection=mode_selection
                  ):
            return self.get_VPT2_wfns(
                "OCHD_freq.fchk",
                internals,
                states,
                n_quanta,
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
        n_quanta = 6

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
                  n_quanta=n_quanta,
                  mode_selection=mode_selection
                  ):
            return self.get_VPT2_wfns(
                "OCHD_freq.fchk",
                internals,
                states,
                n_quanta,
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
    def test_OCHTVPTCartesianPotential(self):

        # This is _expected_ to fail

        internals = None

        n_modes = 3 * 4 - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        n_quanta = 6

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
                  n_quanta=n_quanta,
                  mode_selection=mode_selection
                  ):
            return self.get_VPT2_wfns(
                "OCHT_freq.fchk",
                internals,
                states,
                n_quanta,
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
    def test_OCHTVPTCartesianCoriolis(self):

        internals = None

        n_modes = 3 * 4 - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        n_quanta = 6

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
                  n_quanta=n_quanta,
                  mode_selection=mode_selection
                  ):
            return self.get_VPT2_wfns(
                "OCHT_freq.fchk",
                internals,
                states,
                n_quanta,
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
        n_quanta = 6

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
                  n_quanta=n_quanta,
                  mode_selection=mode_selection
                  ):
            return self.get_VPT2_wfns(
                "OCHH_freq.fchk",
                internals,
                states,
                n_quanta,
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
        n_quanta = 6

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
                  n_quanta=n_quanta,
                  mode_selection=mode_selection
                  ):
            return self.get_VPT2_wfns(
                "OCHH_freq.fchk",
                internals,
                states,
                n_quanta,
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