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
    def get_VPT2_wfns(self, fchk, internals, states, n_quanta,
                      regenerate=False,
                      coupled_states=None,
                      mode_selection=None
                      ):

        hammer = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data(fchk),
            n_quanta=n_quanta,
            internals=internals,
            mode_selection=mode_selection
        )

        wfn_file = os.path.join(self.wfn_file_dir, fchk.replace("fchk", "npz"))
        if regenerate or not os.path.exists(wfn_file):
            wfns = hammer.get_wavefunctions(states, coupled_states=coupled_states)
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

    @debugTest
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
        n_modes = 3
        coupled_states = None#self.get_states(5, n_modes)
        def block(self=self, internals=internals, states=states, coupled_states=coupled_states, n_quanta=n_quanta):
            return self.get_VPT2_wfns(
                "HOT_freq.fchk",
                internals,
                states,
                n_quanta,
                regenerate=True,
                coupled_states=coupled_states
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
        n_quanta = 10
        max_quanta = 5

        states = self.get_states(1, n_modes, max_quanta=max_quanta)
        coupled_states=self.get_states(4, n_modes, max_quanta = max_quanta)

        states = ((0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 1))
        coupled_states = self.get_states(5, n_modes, max_quanta=max_quanta)
        # raise Exception(len(coupled_states))

        # hammer = PerturbationTheoryHamiltonian.from_fchk(
        #     TestManager.test_data("OCHD_freq.fchk"),
        #     n_quanta=n_quanta,
        #     internals=None,
        #     mode_selection=mode_selection
        # )
        #
        # h2w = UnitsData.convert("Hartrees", "Wavenumbers")
        # v4 = hammer.V_terms[2]*h2w
        # inds = np.indices(v4.shape).transpose((1, 2, 3, 4, 0))
        # bleh = v4.flatten()
        # blinds = inds.reshape((len(bleh), 4))
        # agh = "\n".join(
        #     "{0[0]} {0[1]} {0[2]} {0[3]} {1:>7.2f}".format(i, v) for i,v in zip(blinds, bleh)
        # )
        # raise Exception(agh)

        def block(self=self,
                  internals=internals, states=states, coupled_states=coupled_states,
                  n_quanta=n_quanta, mode_selection=mode_selection
                  ):
            return self.get_VPT2_wfns(
                "OCHD_freq.fchk",
                internals,
                states,
                n_quanta,
                regenerate=True,
                coupled_states=coupled_states,
                mode_selection=mode_selection
            )
        # block()
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
        freqs = engs[1:] - engs[0]
        harm_engs = h2w * wfns.zero_order_energies
        harm_freq = harm_engs[1:] - harm_engs[0]

        gaussian_engs = [5235.162,    5171.477]
        gaussian_freqs = [
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
        ]

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