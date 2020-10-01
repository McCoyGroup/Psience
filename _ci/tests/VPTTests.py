from Peeves.TestUtils import *
from unittest import TestCase
from Psience.VPT2 import *
from McUtils.Data import UnitsData
import sys, os, numpy as np
import cProfile, pstats, io

class VPTTests(TestCase):

    def save_wfns(self, file, wfns):
        """
        We save the corrections so that we can reload them later

        :return:
        :rtype:
        """
        np.savez(file,
                 **{k:v for k,v in zip(wfns.corrs._fields, wfns.corrs)}
                 )
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
            PerturbationTheoryCorrections(**np.load(file))
        )

    wfn_file_dir = os.path.expanduser("~/Desktop/")
    def get_VPT2_wfns(self, fchk, internals, states, n_quanta, regenerate=False, coupled_states=None):

        hammer = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data(fchk),
            n_quanta=n_quanta,
            internals=internals
        )

        wfn_file = os.path.join(self.wfn_file_dir, fchk.replace("fchk", "npz"))
        if regenerate or not os.path.exists(wfn_file):
            wfns = hammer.get_wavefunctions(states, coupled_states=coupled_states)
            self.save_wfns(wfn_file, wfns)
        else:
            wfns = self.load_wfns(hammer.molecule, hammer.basis, wfn_file)

        return wfns
    def get_states(self, n_quanta, n_modes):
        import itertools as ip

        return tuple(sorted(
            [p for p in ip.product(*(range(n_quanta+1) for i in range(n_modes))) if sum(p) <= n_quanta],
            key=lambda p: sum(p) + sum(1 for v in p if v != 0) * n_quanta ** (-1) + sum(
                v * n_quanta ** (-i - 2) for i, v in enumerate(p))
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
    def test_HODVPTCartesians(self):
        # this isn't expected to give fully accurate results
        internals = None

        hammer = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("HOD_freq.fchk"),
            n_quanta=6,
            internals=internals
        )
        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        states = (
            (0, 0, 0),
            (0, 0, 1), (0, 1, 0), (1, 0, 0),
            (0, 0, 2), (0, 2, 0), (2, 0, 0),
            (0, 1, 1), (1, 0, 1), (1, 1, 0)
        )
        wfns = hammer.get_wavefunctions(states, coupled_states=None)
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

        print_diffs = False
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
        self.assertLess(np.max(np.abs(my_freqs - gaussian_freqs)), 35)

    @validationTest
    def test_HODVPTInternals(self):

        internals = [
            [0, -1, -1, -1],
            [1, 0, -1, -1],
            [2, 0, 1, -1]
        ]
        states = (
            (0, 0, 0),
            (0, 0, 1), (0, 1, 0), (1, 0, 0),
            (0, 0, 2), (0, 2, 0), (2, 0, 0),
            (0, 1, 1), (1, 0, 1), (1, 1, 0)
        )
        n_quanta = 6
        wfns = self.get_VPT2_wfns(
            "HOD_freq.fchk",
            internals,
            states,
            n_quanta,
            regenerate=False
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

        print_diffs=False
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
        n_quanta = 6
        wfns = self.get_VPT2_wfns(
            "HOD_freq.fchk",
            internals,
            states,
            n_quanta,
            regenerate=False
        )

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")
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
            [1, 0, -1, -1],
            [2, 0, 1, -1]
        ]
        states = (
            (0, 0, 0),
            (0, 0, 1), (0, 1, 0), (1, 0, 0),
            (0, 0, 2), (0, 2, 0), (2, 0, 0),
            (0, 1, 1), (1, 0, 1), (1, 1, 0)
        )
        n_quanta = 6

        wfns = self.get_VPT2_wfns(
            "HOT_freq.fchk",
            internals,
            states,
            n_quanta,
            regenerate=False
        )

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

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

        print_diffs = False
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
    def test_WaterVPT(self):

        hammer = PerturbationTheoryHamiltonian.from_fchk(TestManager.test_data("HOD_freq.fchk"), n_quanta=5)
        h2w = UnitsData.convert("Hartrees", "Wavenumbers")
        get_ind = hammer.get_state_indices
        get_qn = hammer.get_state_quantum_numbers
        np.set_printoptions(suppress=True)

        # np.savetxt("/Users/Mark/Desktop/H0.txt", hammer.H0[:, :])
        # np.savetxt("/Users/Mark/Desktop/H1.txt", hammer.H1[:, :])
        # np.savetxt("/Users/Mark/Desktop/H2.txt", hammer.H2[:, :])
        # coeffs, corrs = hammer.get_corrections([0, 1])

        states = ((0, 0, 0), (0, 0, 1), (0, 1, 0), (1, 0, 0))
        "100 300 010 210 120 030 201 001 111 021 102 012 003"
        c_states = (
            (0, 0, 0),
            (0, 0, 1),
            (0, 0, 3),
            (0, 1, 0),
            (0, 1, 2),
            (0, 2, 1),
            (0, 3, 0),
            (1, 0, 2),
            (1, 0, 0),
            (1, 1, 1),
            (1, 2, 0),
            (2, 0, 1),
            (2, 1, 0),
            (3, 0, 0)
            )
        # print(hammer.H1[get_ind([c_states[0]]), get_ind([c_states[4]])])
        # cf, cr = hammer.get_corrections(c_states[:2], coupled_states=c_states)
        # print(h2w*cr[1])

        # hammer.H1[get_ind([0, 1, 0]), get_ind([0, 4, 0])]

        thresh = 300 / h2w
        coeffs, corrs = hammer.get_corrections(states, coupled_states=None)
        energies = h2w*sum(corrs)
        corrs = [h2w*x for x in corrs]

        print("State Energies:\n",
              *(
                  "{:<1.0f} {:<1.0f} {:<1.0f} {:>5.0f} {:>5.0f} {:>5.0f} {:>5.0f}\n".format(*x) for x in np.round(
                  np.column_stack((states, corrs[0], energies, corrs[0] - corrs[0][0], (energies - energies[0])))
              ))
              )
        print("Target:\n",
              """
 0 0 0  4053  3995     0     0
 0 0 1     -     -  3874  3685
 0 1 0     -     -  2810  2706
 1 0 0     -     -  1422  1383
 """
              )

        # print(h2w*hammer.martin_test())

    @debugTest
    def test_FormaldehydeVPT(self):

        internals = [
            [0, -1, -1, -1],
            [1,  0, -1, -1],
            [2,  1,  0, -1],
            [3,  1,  0,  2]
        ]
        n_modes = 3*4 - 6
        states = self.get_states(2, n_modes)
        n_quanta = 6
        coupled_states = self.get_states(2, n_modes)
        def block(self=self, internals=internals, states=states, coupled_states=coupled_states, n_quanta=n_quanta):
            return self.get_VPT2_wfns(
                "OCHD_freq.fchk",
                internals,
                states,
                n_quanta,
                regenerate=False,
                coupled_states=coupled_states
            )

        block()
        exc, stat_block, wfns = self.profile_block(block)

        if exc is None:
            h2w = UnitsData.convert("Hartrees", "Wavenumbers")
            engs = h2w * wfns.energies
            print(engs[1:] - engs[0])
            Exception(stat_block)
        else:
            try:
                raise exc
            except:
                raise Exception(stat_block)

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")
        engs = h2w * wfns.energies
        print(engs[1:] - engs[0])