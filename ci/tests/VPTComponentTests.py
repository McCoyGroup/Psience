
try:
    from Peeves.TestUtils import *
    from Peeves import BlockProfiler
except:
    pass
from unittest import TestCase

from Psience.VPT2 import *
from Psience.Molecools import Molecule
from Psience.BasisReps import HarmonicOscillatorProductBasis, BasisStateSpace

from McUtils.Data import UnitsData
import McUtils.Plots as plt
import McUtils.Numputils as nput
from McUtils.Scaffolding import *
from McUtils.Parallelizers import SerialNonParallelizer, MultiprocessingParallelizer
from McUtils.Zachary import FiniteDifferenceDerivative

import sys, os, numpy as np, itertools as ip

class VPT2Tests(TestCase):

    #region setup

    gaussian_data = {} # storage to allow me to duplicate Gaussain data less often...
    analytic_data = {} # same idea

    def setUp(self):
        import warnings
        np.seterr(all='raise')
        warnings.filterwarnings('error', category=np.VisibleDeprecationWarning)

        self.h2w = UnitsData.convert("Hartrees", "Wavenumbers")

    def __getstate__(self):
        return {}
    def __setstate__(self, state):
        pass

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
                              mol_spec,
                              internals,
                              states,
                              save_coeffs=False,
                              save_wfns=False,
                              regenerate=True,
                              coupled_states=None,
                              mode_selection=None,
                              potential_terms=None,
                              kinetic_terms=None,
                              coriolis_terms=None,
                              pseudopotential_terms=None,
                              v2=None,
                              t2=None,
                              v3=None,
                              t3=None,
                              t4=None,
                              v4=None,
                              coriolis=None,
                              watson=None,
                              degeneracies=None,
                              selection_rules=None,
                              log=False,
                              verbose=False,
                              pre_run_script=None,
                              get_breakdown=False,
                              parallelized=False,
                              initialization_timeout=1,
                              direct_sum_chunk_size=None,
                              processes=None,
                              checkpoint=None,
                              chunk_size=None,
                              order=2,
                              expansion_order=None
                              , memory_constrained=False
                              , allow_sakurai_degs=False
                              , allow_post_PT_calc=True
                              , modify_degenerate_perturbations=False
                              , gaussian_resonance_handling=False
                              , ignore_odd_order=True
                              , intermediate_normalization=False
                              , state_space_iterations=None
                              , zero_element_warning = False
                              , **solver_opts
                              ):
        if parallelized:
            pverb = verbose == 'all'
            parallelizer = MultiprocessingParallelizer(verbose=pverb,
                                                       processes=processes,
                                                       initialization_timeout=initialization_timeout
                                                       )
        else:
            parallelizer = SerialNonParallelizer()

        if direct_sum_chunk_size is not None:
            BasisStateSpace.direct_sum_chunk_size = direct_sum_chunk_size

        with parallelizer:
            if isinstance(mol_spec, str):
                hammer = PerturbationTheoryHamiltonian.from_fchk(
                    TestManager.test_data(mol_spec),
                    internals=internals,
                    mode_selection=mode_selection,
                    log=log,
                    parallelizer=parallelizer,
                    checkpoint=checkpoint,
                    operator_chunk_size=chunk_size,
                    selection_rules=selection_rules,
                    potential_terms=potential_terms,
                    kinetic_terms=kinetic_terms,
                    coriolis_terms=coriolis_terms,
                    pseudopotential_terms=pseudopotential_terms,
                    include_pseudopotential=False if watson is False else True,
                    coriolis_coupling=False if coriolis is False else True
                )
            else:
                hammer = PerturbationTheoryHamiltonian(
                    mol_spec,
                    mode_selection=mode_selection,
                    logger=log,
                    parallelizer=parallelizer,
                    checkpoint=checkpoint,
                    operator_chunk_size=chunk_size,
                    selection_rules=selection_rules,
                    potential_terms=potential_terms,
                    kinetic_terms=kinetic_terms,
                    coriolis_terms=coriolis_terms,
                    pseudopotential_terms=pseudopotential_terms,
                    include_pseudopotential=False if watson is False else True,
                    include_coriolis_coupling=False if coriolis is False else True
                )

            if pre_run_script is not None:
                pre_run_script(hammer, states)

            if get_breakdown:
                bd = hammer.get_breakdown(states=states, degeneracies=degeneracies, coupled_states=coupled_states)
                return bd, hammer

            # wfn_file = os.path.join(self.wfn_file_dir, fchk.replace("fchk", "npz"))
            # if regenerate or not os.path.exists(wfn_file):

                # if save_coeffs:
                #     coeffs_file = os.path.join(self.wfn_file_dir, fchk.replace(".fchk", "_coeffs.npz"))
                #     np.savez(coeffs_file,
                #              G=hammer.H0.computers[0].operator.coeffs,
                #              F=hammer.H0.computers[1].operator.coeffs,
                #              dGdQ=hammer.H1.computers[0].operator.coeffs,
                #              dFdQ=hammer.H1.computers[1].operator.coeffs,
                #              dGdQQ=hammer.H2.computers[0].operator.coeffs,
                #              dFdQQ=hammer.H2.computers[1].operator.coeffs,
                #              coriolis=hammer.H2.computers[2].operator.coeffs,
                #              watson=hammer.H2.computers[3].operator.coeffs
                #              )

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
                if internals is None and isinstance(t4, (int, float)) and t4 == 0:
                    if coriolis is None:
                        coriolis = 0
            if v4 is not None:
                hammer.H2.computers[1].operator.coeffs = v4
            if coriolis is not None:
                hammer.H2.computers[2].operator.coeffs = coriolis
                if watson is None and isinstance(coriolis, (int, float)) and coriolis == 0:
                    watson = 0
            if isinstance(watson, (int, float, np.integer, np.floating, np.ndarray)) and watson is not False and watson is not True:
                hammer.H2.computers[3].operator.coeffs = watson

            wfns = hammer.get_wavefunctions(states, coupled_states=coupled_states, degeneracies=degeneracies,
                                            order=order,
                                            expansion_order=expansion_order
                                            , memory_constrained=memory_constrained
                                            , allow_sakurai_degs=allow_sakurai_degs
                                            , allow_post_PT_calc=allow_post_PT_calc
                                            , modify_degenerate_perturbations=modify_degenerate_perturbations
                                            , gaussian_resonance_handling=gaussian_resonance_handling
                                            , ignore_odd_order_energies=ignore_odd_order
                                            , zero_element_warning = zero_element_warning
                                            , intermediate_normalization=intermediate_normalization
                                            , state_space_iterations=state_space_iterations
                                            , verbose=verbose
                                            , **solver_opts
                                            )

        return wfns, hammer

    def get_VPT2_wfns(self,
                      mol_spec,
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
                      selection_rules=None,
                      log=False,
                      verbose=False,
                      pre_run_script=None,
                      parallelized=False,
                      initialization_timeout=1,
                      checkpoint=None,
                      get_breakdown=False,
                      chunk_size=None,
                      order=2,
                      expansion_order=None
                      , memory_constrained=False
                      , allow_sakurai_degs=False
                      , allow_post_PT_calc=True
                      , ignore_odd_order=True
                      , gaussian_resonance_handling=False
                      , modify_degenerate_perturbations=False
                      , intermediate_normalization=False
                      , state_space_iterations=None
                      , zero_element_warning=False
                      ):
        return self.get_VPT2_wfns_and_ham(
            mol_spec,
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
            selection_rules=selection_rules,
            log=log,
            verbose=verbose,
            parallelized=parallelized,
            pre_run_script=pre_run_script,
            checkpoint=checkpoint,
            get_breakdown=get_breakdown,
            chunk_size=chunk_size,
            order=order,
            expansion_order=expansion_order
            , memory_constrained=memory_constrained
            , allow_sakurai_degs=allow_sakurai_degs
            , allow_post_PT_calc=allow_post_PT_calc
            , gaussian_resonance_handling=gaussian_resonance_handling
            , ignore_odd_order=ignore_odd_order
            , modify_degenerate_perturbations=modify_degenerate_perturbations
            , intermediate_normalization=intermediate_normalization
            , state_space_iterations=state_space_iterations
            , zero_element_warning=zero_element_warning
        )[0]

    @staticmethod
    def get_states(n_quanta, n_modes, max_quanta=None, target_modes=None):
        whee = [np.flip(x) for x in BasisStateSpace.from_quanta(
            HarmonicOscillatorProductBasis(n_modes),
            range(n_quanta)
        ).excitations]
        if target_modes is not None:
            whee = [
                p for p in whee if sum(p) == 0 or any(p[i] > 0 for i in target_modes)
            ]
        return whee

    @staticmethod
    def get_degenerate_polyad_space(states, polyadic_pairs):
        # we build a graph of connected states by the polyadic rules
        polyadic_pairs = np.array(polyadic_pairs)
        states = np.array(states)
        poss_degs = [[] for _ in states]
        check_list = states.tolist()
        for n,s in enumerate(states): # build list-of-lists structure
            for i, nt_spec in enumerate(polyadic_pairs):
                if np.all(s - nt_spec[0] >= 0):
                    new = (s - nt_spec[0] + nt_spec[1]).tolist()
                    if new not in check_list:
                        check_list.append(new)
                        poss_degs.append([])
                    # else:
                    #     poss_degs[idx].append(slist)
                    poss_degs[n].append(new)

        # from the populated lists build the real connection graph
        groups = [[] for _ in check_list]
        new_checks = []
        for i,s1 in enumerate(check_list):
            if s1 not in new_checks:
                new_checks.append(s1)
                groups[i].append(s1)
                groups[i].extend(poss_degs[i])
                for s2,p in zip(check_list[i+1:], poss_degs[i+1:]):
                    if s2 not in new_checks:
                        if s2 in poss_degs[i]:
                            new_checks.append(s2)
                            groups[i].extend(p)

        return [g for g in groups if len(g) > 1]

        # degeneracies = {}
        # for s in states:
        #     for i, nt_spec in enumerate(nt_specs):
        #         nt = nt_spec * s
        #         if np.all(nt % 0 == 0):
        #             if i not in degeneracies:
        #                 degeneracies[i] = {}
        #             key = tuple(nt.astype(int))
        #             if key in degeneracies[i]:
        #                 degeneracies[i].append(s)
        #             else:
        #                 degeneracies[i] = [s]

    def print_energy_block(self, tag, wfns, states, zpe, freqs, real_fmt='{:>10.3f}', dash_fmt='{:>10}'):

        if wfns is not None:
            print(
                tag,
                wfns.format_energies_table(
                    states=states, zpe=zpe, freqs=freqs,
                    real_fmt=real_fmt, dash_fmt=dash_fmt).replace(
                    "\n", "\n  "
                ),
                sep="\n  "
            )
        else:
            print(
                tag,
                PerturbationTheoryWavefunctions._format_energies_table(
                    states, zpe, freqs,
                    real_fmt=real_fmt, dash_fmt=dash_fmt
                ).replace(
                    "\n", "\n  "
                ),
                sep="\n  "
            )

    def run_PT_test(self,
                    tag, mol_spec, internals, mode_selection,
                    states, gaussian_energies, gaussian_freqs,
                    print_report=False,
                    print_diffs=True,
                    log=False,
                    energy_order=None,
                    nielsen_tolerance=1,
                    gaussian_tolerance=1,
                    print_profile=False,
                    profile_filter=None,
                    profiling_mode=None,
                    pre_wfns_script=None,
                    post_wfns_script=None,
                    print_specs=True,
                    calculate_intensities=False,
                    print_x=False,
                    invert_x=False,
                    **opts
                    ):
        with BlockProfiler(tag, print_res=print_profile, mode=profiling_mode, inactive=not print_profile):#, filter=profile_filter):
            wfns, hammer = self.get_VPT2_wfns_and_ham(
                mol_spec,
                internals,
                states,
                regenerate=True,
                pre_run_script=pre_wfns_script,
                mode_selection=mode_selection
                , log=log,
                **opts
            )
        if post_wfns_script is not None:
           post_wfns_script(wfns, hammer)

        # Pure PT energies
        h2w = UnitsData.convert("Hartrees", "Wavenumbers")
        if energy_order is None:
            engs = h2w * wfns.energies
        else:
            engs = h2w * wfns.energies_to_order(energy_order)

        # raise Exception(h2w*wfns.corrs.energy_corrs)
        freqs = engs[1:] - engs[0]
        harm_engs = h2w * wfns.zero_order_energies
        harm_freq = harm_engs[1:] - harm_engs[0]

        my_energies = np.array([harm_engs[0], engs[0]])
        my_freqs = np.column_stack([harm_freq, freqs])

        if print_report:
            self.print_energy_block("State Energies:", wfns, states, my_energies, my_freqs)

        if internals is None and nielsen_tolerance is not None:
            # Energies from Nielsen expressions
            e_harm, e_corr, x = hammer.get_Nielsen_energies(states, return_split=True)
            if print_x:
                with np.printoptions(linewidth=10000000):  # infinite line width basically...
                    print("="*25 + "X-Matrix:" + "="*25,
                          repr(x * h2w).strip("array()").replace("       ", " "),
                          "=" * 50,
                          sep="\n"
                          )

            e_corr = np.sum(e_corr, axis=0)
            energies = h2w * (e_harm + e_corr)
            zero_ord = h2w * e_harm

            nielsen_engs = np.array([zero_ord[0], energies[0]])
            nielsen_freqs = np.column_stack([zero_ord[1:] - zero_ord[0], energies[1:] - energies[0]])

            if print_report:
                self.print_energy_block("Nielsen Energies:", None, states, nielsen_engs, nielsen_freqs)

            if print_diffs:
                self.print_energy_block("Nielsen Difference Energies:", None, states, my_energies - nielsen_engs,
                                        my_freqs - nielsen_freqs)

        if gaussian_freqs is not None:
            ns = min(len(freqs), len(gaussian_freqs))
            if print_report:
                self.print_energy_block("Reference Energies:", None, states, gaussian_energies, gaussian_freqs[:ns])

            if print_diffs:
                self.print_energy_block("Reference Difference Energies:", None, states, my_energies - gaussian_energies,
                                        my_freqs[:ns] - gaussian_freqs[:ns])

        if calculate_intensities:
            engs = h2w * wfns.energies
            freqs = engs - engs[0]
            ints = wfns.intensities

            harm_engs = h2w * wfns.zero_order_energies
            harm_freqs = harm_engs - harm_engs[0]
            harm_ints = wfns.zero_order_intensities

            print_specs = True
            if print_specs:
                print(wfns.format_intensities_table())
                # n_modes= wfns.corrs.total_basis.ndim
                # padding = np.max([len(str("0 " * n_modes)), 1]) + 1
                # bar = "Frequency    Intensity"
                # spacer = "    "
                # report = (
                #         " " * (padding+2) + " " * (len("Frequency")-2) + "Harmonic" + " " * (len("    Intensity")+2 - len("Harmonic")) + spacer + " " * (len("Frequency")-2) + "Anharmonic\n"
                #         "State" + " " * (padding - 3) + bar + spacer + bar + "\n"
                #          )+ "\n".join(
                #     (
                #         (
                #                 "  " + ("{:<1.0f} " * n_modes) + " "
                #                 + "{:>9.2f} {:>12.4f}"
                #                 + spacer
                #                 + "{:>9.2f} {:>12.4f}"
                #         ).format(*s, hf, hi, f, i)
                #         for s, hf, hi, f, i in zip(states[1:], harm_freqs[1:], harm_ints[1:], freqs[1:], ints[1:])
                #     ))
                # print(report)

        # if invert_x:
        #
        #     hammer.get_

        if internals is None and nielsen_tolerance is not None:
            self.assertLess(np.max(np.abs(my_freqs - nielsen_freqs)), nielsen_tolerance)

        if gaussian_freqs is not None:
            self.assertLess(np.max(np.abs(my_freqs[:ns] - gaussian_freqs[:ns])), gaussian_tolerance)

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

    def write_intensity_breakdown(self, s, all_wfns, plot_spec, write_wavefunctions=True):
        """

        :param s:
        :type s:
        :param all_wfns:
        :type all_wfns: dict
        :param plot_spec:
        :type plot_spec:
        :return:
        :rtype:
        """

        import io, csv
        from Psience.Spectra import DiscreteSpectrum

        writer = csv.writer(s)
        for k, wnf in all_wfns.items():
            wnf = wnf #type: PerturbationTheoryWavefunctions
            writer.writerow([k + " Wavefunctions"])
            for dp in [
                wnf.DipolePartitioningMethod.Intuitive,
                wnf.DipolePartitioningMethod.Standard
            ]:
                writer.writerow([" ", dp.value.capitalize() + " Dipole Partitioning"])
                wnf.dipole_partitioning = dp
                bd = wnf.generate_intensity_breakdown(include_wavefunctions=write_wavefunctions)
                wnf.write_CSV_breakdown(s, bd, padding=[" "] * 2)
                if plot_spec:
                    freqs = self.h2w * bd["frequencies"]
                    ploot = None
                    colors = ["red", "blue", "green", "black"]
                    for i, vals in enumerate(bd["breakdowns"].values()):
                        spec = DiscreteSpectrum(freqs, vals["intensities"])
                        if ploot is None:
                            ploot = spec.plot()
                        else:
                            ploot = spec.plot(figure=ploot, plot_style=dict(linefmt=colors[i]))
                    ploot.show()
            # break

    #endregion

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
            [2, 2, 1,  -90.66588863],
            [3, 2, 1,    82.31540784],
            [1, 3, 1,   429.45742955],
            [2, 3, 1,    82.31540784],
            [3, 3, 1,  -153.50694953],
            [1, 1, 2,    41.16039749],
            [2, 1, 2,   -90.73683527],
            [3, 1, 2,    82.32690754],
            [1, 2, 2,   -90.73683527],
            [2, 2, 2, -1588.89967595],
            [3, 2, 2,    91.00951548],
            [1, 3, 2,    82.32690754],
            [2, 3, 2,    91.00951548],
            [3, 3, 2,  -172.34377091],
            [1, 1, 3,   430.44067645],
            [2, 1, 3,    82.33125106],
            [3, 1, 3,  -153.78923868],
            [1, 2, 3,    82.33125106],
            [2, 2, 3,    90.96564514],
            [3, 2, 3,  -172.45333790],
            [1, 3, 3,  -153.78923868],
            [2, 3, 3,  -172.45333790],
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

        print_report = True
        if print_report:
            print("-" * 50, "x", "-" * 50)
            for a, b in zip(x, gaussian_x):
                print(
                    ("[" + "{:>8.3f}" * 3 + "]").format(*a)
                    + ("[" + "{:>8.3f}" * 3 + "]").format(*b)
                )
            print("-" * 50, "y", "-" * 50)
            for a, b in zip(y, gaussian_y):
                print(
                    ("[" + "{:>8.3f}" * 3 + "]").format(*a)
                    + ("[" + "{:>8.3f}" * 3 + "]").format(*b)
                )
            print("-" * 50, "z", "-" * 50)
            for a, b in zip(z, gaussian_z):
                print(
                    ("[" + "{:>8.3f}" * 3 + "]").format(*a)
                    + ("[" + "{:>8.3f}" * 3 + "]").format(*b)
                )

        sum_diff = np.abs(sum([x, y, z])) - np.abs(sum([gaussian_x, gaussian_y, gaussian_z]))
        # print(np.round(sum_diff, 3))
        self.assertLess(np.max(sum_diff), .001)

    @validationTest
    def test_HOHCoriolisCouplings(self):

        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("HOH_freq.fchk"),
            internals=None
        )

        x, y, z = ham.coriolis_terms.get_zetas()

        x_els, y_els, z_els = [
            [

            ],
            [
                [2, 1, -0.960228e-2],
                [3, 1,  0.999954]
            ],
            [

            ]
        ]
        mapping = [3, 2, 1]  # remap modes to go from highest to lowest frequency

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

        print_report = True
        if print_report:
            print("-" * 50, "x", "-" * 50)
            fmt = "{:>8.5f}"
            for a, b in zip(x, gaussian_x):
                print(
                    ("[" +  fmt * 3 + "]").format(*a)
                    + ("[" + fmt * 3 + "]").format(*b)
                )
            print("-" * 50, "y", "-" * 50)
            for a, b in zip(y, gaussian_y):
                print(
                    ("[" + fmt * 3 + "]").format(*a)
                    + ("[" + fmt * 3 + "]").format(*b)
                )
            print("-" * 50, "z", "-" * 50)
            for a, b in zip(z, gaussian_z):
                print(
                    ("[" + fmt * 3 + "]").format(*a)
                    + ("[" + fmt * 3 + "]").format(*b)
                )

        sum_diff = np.abs(sum([x, y, z])) - np.abs(sum([gaussian_x, gaussian_y, gaussian_z]))
        # print(np.round(sum_diff, 3))
        self.assertLess(np.max(sum_diff), .001)

    @validationTest
    def test_HOONOCoriolisCouplings(self):

        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("HOONO_freq.fchk"),
            internals=None
        )

        x, y, z = ham.coriolis_terms.get_zetas()

        x_els, y_els, z_els = [
            [
                [8,     1,        0.36276],
                [8,     2,       -0.49613],
                [8,     3,        0.68484],
                [8,     4,       -0.33710],
                [8,     5,        0.14662],
                [8,     6,       -0.11881],
                [8,     7,       -0.06335],
                [9,     1,        0.57992],
                [9,     2,        0.45045],
                [9,     3,        0.33040],
                [9,     4,        0.52671],
                [9,     5,       -0.13830],
                [9,     6,        0.19700],
                [9,     7,        0.12751]
            ],
            [
                [8,      1,      -0.36436],
                [8,      2,       0.16770],
                [8,      3,       0.27720],
                [8,      4,      -0.31168],
                [8,      5,      -0.14792],
                [8,      6,       0.69773],
                [8,      7,      -0.39550],
                [9,      1,      -0.61150],
                [9,      2,       0.17607],
                [9,      3,       0.55239],
                [9,      4,       0.33293],
                [9,      5,      -0.03255],
                [9,      6,      -0.42064],
                [9,      7,      -0.03290]

            ],
            [
              [2,      1,             0.30785],
              [3,      1,             0.88920],
              [3,      2,            -0.04935],
              [4,      1,             0.25299],
              [4,      2,             0.40061],
              [4,      3,            -0.25043],
              [5,      1,            -0.02792],
              [5,      2,             0.57161],
              [5,      3,             0.07597],
              [5,      4,            -0.34347],
              [6,      1,             0.10760],
              [6,      2,            -0.39341],
              [6,      3,             0.31003],
              [6,      4,            -0.69746],
              [6,      5,             0.16906],
              [7,      1,             0.16256],
              [7,      2,            -0.27993],
              [7,      3,             0.17555],
              [7,      4,             0.01092],
              [7,      5,             0.66449],
              [7,      6,             0.13739]
            ]
        ]
        mapping = [9, 8, 7, 6, 5, 4, 2, 3, 1]  # remap modes to go from highest to lowest frequency

        gaussian_x = np.zeros((9, 9))
        for i, j, v in x_els:
            i = mapping[i - 1] - 1
            j = mapping[j - 1] - 1
            gaussian_x[i, j] = v
            gaussian_x[j, i] = -v
        gaussian_y = np.zeros((9, 9))
        for i, j, v in y_els:
            i = mapping[i - 1] - 1
            j = mapping[j - 1] - 1
            gaussian_y[i, j] = -v
            gaussian_y[j, i] = v
        gaussian_z = np.zeros((9, 9))
        for i, j, v in z_els:
            i = mapping[i - 1] - 1
            j = mapping[j - 1] - 1
            gaussian_z[i, j] = v
            gaussian_z[j, i] = -v

        print_report = True
        if print_report:
            print("-" * 50, "x", "-" * 50)
            fmt = "{:>8.5f}"
            for a, b in zip(x, gaussian_x):
                print(
                    ("[" + fmt * 3 + "]").format(*a)
                    + ("[" + fmt * 3 + "]").format(*b)
                )
            print("-" * 50, "y", "-" * 50)
            for a, b in zip(y, gaussian_y):
                print(
                    ("[" + fmt * 3 + "]").format(*a)
                    + ("[" + fmt * 3 + "]").format(*b)
                )
            print("-" * 50, "z", "-" * 50)
            for a, b in zip(z, gaussian_z):
                print(
                    ("[" + fmt * 3 + "]").format(*a)
                    + ("[" + fmt * 3 + "]").format(*b)
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

    @validationTest
    def test_HOHNielsenEnergies(self):
        states = (
            (0, 0, 0),
            (0, 0, 1), (0, 1, 0), (1, 0, 0),
            (0, 0, 2), (0, 2, 0), (2, 0, 0),
            (0, 1, 1), (1, 0, 1), (1, 1, 0)
        )

        hammer = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("HOH_freq.fchk"),
            internals=None,
            log=False
        )

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        gaussian_x = np.array([
            [
                [-0.965669e+02, -0.357168e+03,  0.428508e+02],
                [-0.357168e+03, -0.921241e+02,  0.553701e+02],
                [ 0.428508e+02,  0.553701e+02, -0.122624e+02],
            ],
            [
                [ 0.475281e+02,  0.190802e+03, -0.924902e+02],
                [ 0.190802e+03,  0.480711e+02, -0.754058e+02],
                [-0.924902e+02, -0.754058e+02, -0.176138e+01],
            ],
            [
                [0.000000e+00, 0.172747e-02, 0.265776e+02],
                [0.172747e-02, 0.000000e+00, 0.000000e+00],
                [0.265776e+02, 0.000000e+00, 0.000000e+00],
            ]
        ])

        gaussian_x = np.array([x[np.ix_((2, 1, 0), (2, 1, 0))] for x in gaussian_x])

        e_harm, e_corr, x = hammer.get_Nielsen_energies(states, return_split=True)
        e_harm, e_corr_gaussian = hammer._get_Nielsen_energies_from_x(states, hammer.modes.freqs, gaussian_x, return_split=True)

        # raise Exception(np.round(gaussian_x[2], 4), np.round(x[2] * h2w, 4))

        self.assertTrue(
            np.allclose(e_corr*h2w, e_corr_gaussian, atol=1),
            msg="diff: \n {}".format(
                np.round(e_corr_gaussian - e_corr*h2w)
            )
        )

        x = x * h2w
        self.assertTrue(
            np.allclose(gaussian_x, x, atol=5),
            msg="wat \n {} \n {}".format(
                gaussian_x, x
            )
        )

        # wfns = hammer.get_wavefunctions(states, coupled_states=None)

        e_corr = np.sum(e_corr, axis=0)
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
            self.print_energy_block("Gaussian Energies:", states, gaussian_energies, gaussian_freqs)
            self.print_energy_block("Nielsen Energies:", states, my_energies, my_freqs)

        print_diffs = True
        if print_diffs:
            self.print_energy_block("Difference Energies:", states, my_energies - gaussian_energies, my_freqs - gaussian_freqs[:len(my_freqs)])

        self.assertLess(np.max(np.abs(my_freqs - gaussian_freqs[:len(my_freqs)])), 1.5)

    #endregion

    #region Test Systems

    #region Analytic Models

    @inactiveTest
    def test_TwoMorseCartesiansDegenerate(self):

        # internals = [
        #     [0, -1,  -1, -1],
        #     [1,  0,  -1, -1],
        #     [2,  0,   1, -1],
        #     [3,  2,   1,  0]
        # ]
        internals = None

        re_1 = 1.0
        re_2 = 1.0

        mol0 = Molecule(
            ["O", "H", "O", "H"],
            np.array([
                [0.000000, 0.000000, 0.000000],
                [re_1,     0.000000, 0.000000],
                [0.000000, 1.000000, 0.000000],
                [re_2,     1.000000, 0.000000]
                ]),
            zmatrix=internals,
            guess_bonds=False,
        )

        # masses = mol0.masses
        # anchor_pos = np.average(mol0.coords[:2], weights=masses[:2]/np.sum(masses[:2]), axis=0)
        # mol0.coords[2] += anchor_pos

        # put in PAF
        mol = mol0
        # mol = mol0.principle_axis_frame()(mol0)
        # raise Exception(mol.coords)

        from McUtils.Zachary import FiniteDifferenceDerivative

        masses = mol0.masses * UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")

        w2h = UnitsData.convert("Wavenumbers", "Hartrees")
        w1 = 3500*w2h; wx1 = 100*w2h
        mu_1 = (1 / masses[0] + 1 / masses[1]) ** (-1)
        De_1 = (w1 ** 2) / (4*wx1)
        a_1 = np.sqrt(2 * mu_1 * wx1)

        w2 = 3500*w2h; wx2 = 100*w2h
        mu_2 = (1 / masses[2] + 1 / masses[3]) ** (-1)
        De_2 = (w2 ** 2) / (4 * wx2)
        a_2 = np.sqrt(2 * mu_2 * wx2)

        # mass_weights = masses[:2] / np.sum(masses[:2])
        def two_morse(carts, *,
                      De_1 = De_1, a_1=a_1, re_1=re_1,
                      De_2 = De_2, a_2=a_2, re_2=re_2
                      ):

            # anchor_pos = np.average(carts[..., :2, :], weights=mass_weights, axis=-2)
            r1 = nput.vec_norms(carts[..., 1, :] - carts[..., 0, :]) - re_1
            r2 = nput.vec_norms(carts[..., 3, :] - carts[..., 2, :]) - re_2

            return (
                De_1 * (1 - np.exp(-a_1 * r1))**2
                + De_2 * (1 - np.exp(-a_2 * r2))**2
            )

        deriv_gen = FiniteDifferenceDerivative(two_morse,
                                               function_shape=((None, None), 0),
                                               mesh_spacing=1e-3,
                                               stencil=5
                                               ).derivatives(mol.coords)

        # v1, v2 = deriv_gen.derivative_tensor([1, 2])
        v1, v2, v3, v4 = deriv_gen.derivative_tensor([1, 2, 3, 4])

        mol.potential_derivatives = [v1, v2, v3, v4]

        modes = mol.normal_modes[-1, -2]
        stretches = modes.basis.matrix.T
        symm = (stretches[1] + stretches[0])*np.sqrt(1/2)
        asym = (stretches[1] - stretches[0])*np.sqrt(1/2)
        new_mat = np.array([symm, asym]).T
        modes.basis.matrix = new_mat
        modes.basis.inverse = new_mat.T
        ham = PerturbationTheoryHamiltonian(mol, modes=modes, n_quanta=50, include_coriolis_coupling=False, include_pseudo_potential=False)

        # mode_selection=[-1, -2]
        # ham = PerturbationTheoryHamiltonian(mol, mode_selection=mode_selection, n_quanta=50)

        # import json
        # raise Exception(json.dumps(ham.V_terms[1].tolist()))

        states = ham.basis.get_state_space(range(5))
        wfns = ham.get_wavefunctions(
            states
            , degeneracies=(
                ([1, 0], [0, 1]),
                ([2, 0], [1, 1], [0, 2]),
                ([3, 0], [2, 1], [1, 2], [0, 3]),
                ([4, 0], [3, 1], [2, 2], [1, 3], [0, 4])
            )
            , allow_sakurai_degs=False
            , allow_post_PT_calc=True
            , order=2
            , coupled_states = [self.get_states(11, 2), self.get_states(11, 2)]
        )

        with JSONCheckpointer(os.path.expanduser('~/Desktop/test_dat.json')) as woof:
            woof['h1'] = wfns.corrs.hams[1].asarray()
            woof['h2'] = wfns.corrs.hams[2].asarray()
            woof['states'] = wfns.corrs.total_basis.excitations
        # raise Exception(np.diag(wfns.corrs.hams[0].asarray()))
        states = states.excitations

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        # wfns = hammer.get_wavefunctions(states, coupled_states=None)
        energies = h2w * wfns.energies
        zero_ord = h2w * wfns.zero_order_energies

        print_report = True
        n_modes = len(states[0])
        harm_engs = zero_ord
        engs = energies
        harm_freq = zero_ord[1:] - zero_ord[0]
        freqs = energies[1:] - energies[0]
        if print_report:
            print("State Energies:\n",
                  ('0 ' * n_modes + "{:>8.3f} {:>8.3f} {:>8} {:>8}\n").format(harm_engs[0], engs[0], "-", "-"),
                  *(
                      ('{:<1.0f} ' * n_modes + "{:>8} {:>8} {:>8.3f} {:>8.3f}\n").format(*s, "-", "-", e1, e2) for
                      s, e1, e2 in
                      zip(states[1:], harm_freq, freqs)
                  )
                  )


    class single_morse:
        # mass_weights = masses[:2] / np.sum(masses[:2])
        def __init__(self, De_1, a_1, re_1):  # , De_2, a_2, re_2, kb):
            self.De_1 = De_1
            self.a_1 = a_1
            self.re_1 = re_1

        def __call__(self, carts):
            # anchor_pos = np.average(carts[..., :2, :], weights=mass_weights, axis=-2)
            r1_vecs = carts[..., 1, :] - carts[..., 0, :]
            r1 = nput.vec_norms(r1_vecs) - self.re_1

            return self.De_1 * (1 - np.exp(-self.a_1 * r1)) ** 2
    analytic_data['Morse3500/50'] = {
        'zpe': np.array([1750.000, 1737.5000000]),
        'freqs': np.array([
            [ 3500.000, 3400.0000000],
            [ 7000.000, 6700.0000000],
            [10500.000, 9900.0000000],

            [14000.000, 13000.000000],
            [17500.000, 16000.000000]
        ])
    }
    @validationTest
    #paper
    def test_OneMorseCartesiansNonDeg(self):

        import warnings
        np.seterr(all='raise')
        warnings.filterwarnings('error', category=np.VisibleDeprecationWarning)

        tag="TwoMorseCartesiansNonDeg"

        # Set up system

        re_1 = 1.0
        re_2 = 1.0
        internals = [
            [0, -1, -1, -1],
            [1,  0, -1, -1],
            [2,  0,  1, -1]
        ]
        # internals = None
        mol = Molecule(
            ["O", "H", "H"],
            np.array([
                [0.000000, 0.000000, 0.000000],
                [re_1, 0.000000, 0.000000],
                [0.000000, re_2, 0.000000],
            ]),
            zmatrix=internals,
            guess_bonds=False,
        )

        masses = mol.masses * UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")

        w2h = UnitsData.convert("Wavenumbers", "Hartrees")
        w1 = 3500 * w2h; wx1 = 50 * w2h
        mu_1 = (1 / masses[0] + 1 / masses[1]) ** (-1)
        De_1 = (w1 ** 2) / (4 * wx1)
        a_1 = np.sqrt(2 * mu_1 * wx1)

        morse = self.single_morse(De_1, a_1, re_1)

        deriv_gen = FiniteDifferenceDerivative(morse,
                                               function_shape=((None, None), 0),
                                               mesh_spacing=1e-3,
                                               stencil=7,
                                               parallelizer=MultiprocessingParallelizer()#verbose=True)
                                               ).derivatives(mol.coords)
        # v1, v2, v3, v4 = deriv_gen.derivative_tensor([1, 2, 3, 4, 5, 6])

        # raise Exception()

        import time
        start = time.time()
        mol.potential_derivatives = deriv_gen.derivative_tensor([1, 2])#, 3, 4, 5, 6])
        end = time.time()
        print('getting input derivs took {:.3f}s'.format(end-start))
        # raise Exception([x.shape for x in mol.potential_derivatives])

        n_atoms = 3
        n_modes = 3 * n_atoms - 6
        mode_selection = [-1]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(5, n_modes)

        print_report = True
        reference_energies = self.analytic_data['Morse3500/50']['zpe']
        reference_freqs = self.analytic_data['Morse3500/50']['freqs']

        # def post_wfns_script(wfns, ham):
        #     plt.ArrayPlot(wfns.corrs.hams[2].asarray()).show()
        #     raise Exception(wfns.corrs.hams[2].asarray())

        post_wfns_script = None

        self.run_PT_test(
            tag,
            mol,
            internals,
            mode_selection,
            states,
            reference_energies,
            reference_freqs,
            log=True,
            verbose=True,
            post_wfns_script=post_wfns_script,
            print_report=print_report,
            # ignore_odd_order=False,
            intermediate_normalization=False,
            order=6,
            expansion_order=6
            , potential_terms=[
                np.array([[0.01594717]]),
                np.array([[[-0.00808669]]]), #np.array([[[0.]]]),
                np.array([[[[0.00318943]]]]), #np.array([[[[0.]]]])
                np.array([[[[[-0.00115532]]]]]),
                np.array([[[[[[0.000403561]]]]]]),
                np.array([[[[[[[-0.000138629]]]]]]]),
                np.array([[[[[[[[0.0000472371]]]]]]]]) # for some reason my FD was giving me np.array([[[[[[0.00038929]]]]]])
            ]
            , kinetic_terms=[
                np.array([[0.01594717]]),
                0,
                0,
                0,
                0,
                0,
                0
            ]
            , watson=False
            # pseudopotential_terms = [
            #     None,
            #     None,
            #     None,
            #     None,
            #     None
            # ]
            # , parallelized=True
        )

    @validationTest
    def test_OneMorseCartesiansNonDegMostlyNumerical(self):

        import warnings
        np.seterr(all='raise')
        warnings.filterwarnings('error', category=np.VisibleDeprecationWarning)

        tag = "TwoMorseCartesiansNonDeg"

        # Set up system

        re_1 = 1.0
        re_2 = 1.0
        internals = [
            [0, -1, -1, -1],
            [1, 0, -1, -1],
            [2, 0, 1, -1]
        ]
        # internals = None
        mol = Molecule(
            ["O", "H", "H"],
            np.array([
                [0.000000, 0.000000, 0.000000],
                [re_1, 0.000000, 0.000000],
                [0.000000, re_2, 0.000000],
            ]),
            zmatrix=internals,
            guess_bonds=False,
        )

        masses = mol.masses * UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")

        w2h = UnitsData.convert("Wavenumbers", "Hartrees")
        w1 = 3500 * w2h;
        wx1 = 50 * w2h
        mu_1 = (1 / masses[0] + 1 / masses[1]) ** (-1)
        De_1 = (w1 ** 2) / (4 * wx1)
        a_1 = np.sqrt(2 * mu_1 * wx1)

        morse = self.single_morse(De_1, a_1, re_1)

        deriv_gen = FiniteDifferenceDerivative(morse,
                                               function_shape=((None, None), 0),
                                               mesh_spacing=1e-3,
                                               stencil=7,
                                               parallelizer=MultiprocessingParallelizer()  # verbose=True)
                                               ).derivatives(mol.coords)
        # v1, v2, v3, v4 = deriv_gen.derivative_tensor([1, 2, 3, 4, 5, 6])

        # raise Exception()

        import time
        start = time.time()
        mol.potential_derivatives = deriv_gen.derivative_tensor([1, 2, 3, 4, 5, 6])
        end = time.time()
        print('getting input derivs took {:.3f}s'.format(end - start))
        # raise Exception([x.shape for x in mol.potential_derivatives])

        n_atoms = 3
        n_modes = 3 * n_atoms - 6
        mode_selection = [-1]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(5, n_modes)

        print_report = False
        reference_energies = self.analytic_data['Morse3500/50']['zpe']
        reference_freqs = self.analytic_data['Morse3500/50']['freqs']

        # def post_wfns_script(wfns, ham):
        #     plt.ArrayPlot(wfns.corrs.hams[2].asarray()).show()
        #     raise Exception(wfns.corrs.hams[2].asarray())

        post_wfns_script = None

        self.run_PT_test(
            tag,
            mol,
            internals,
            mode_selection,
            states,
            reference_energies,
            reference_freqs,
            log=True,
            verbose=True,
            post_wfns_script=post_wfns_script,
            print_report=print_report,
            gaussian_tolerance=2,
            # ignore_odd_order=False,
            intermediate_normalization=False,
            order=4,
            expansion_order=4
            # , potential_terms=[
            #     np.array([[0.01594717]]),
            #     np.array([[[-0.00808669]]]), #np.array([[[0.]]]),
            #     np.array([[[[0.00318943]]]]), #np.array([[[[0.]]]])
            #     np.array([[[[[-0.00115532]]]]]),
            #     np.array([[[[[[0.000403561]]]]]]),
            #     np.array([[[[[[[-0.000138629]]]]]]]),
            #     np.array([[[[[[[[0.0000472371]]]]]]]]) # for some reason my FD was giving me np.array([[[[[[0.00038929]]]]]])
            # ]
            , kinetic_terms=[
                None,#np.array([[0.01594717]]),
                None,
                None,
                0,
                0,
                0,
                0
            ]
            , watson=False
            # pseudopotential_terms = [
            #     None,
            #     None,
            #     None,
            #     None,
            #     None
            # ]
            , parallelized=True
        )

    class harmonically_coupled_morse:
        # mass_weights = masses[:2] / np.sum(masses[:2])
        def __init__(self,
                     De_1, a_1, re_1,
                     De_2, a_2, re_2,
                     kb, b_e
                     ):
            self.De_1 = De_1
            self.a_1 = a_1
            self.re_1 = re_1
            self.De_2 = De_2
            self.a_2 = a_2
            self.re_2 = re_2
            self.kb = kb
            self.b_e = b_e

        def __call__(self, carts):
            v1 = carts[..., 1, :] - carts[..., 0, :]
            v2 = carts[..., 2, :] - carts[..., 0, :]
            r1 = nput.vec_norms(v1) - self.re_1
            r2 = nput.vec_norms(v2) - self.re_2
            bend, _ = nput.vec_angles(v1, v2)
            bend = bend - self.b_e

            return (
                    self.De_1 * (1 - np.exp(-self.a_1 * r1)) ** 2
                    + self.De_2 * (1 - np.exp(-self.a_2 * r2)) ** 2
                    + self.kb * bend**2
            )
    analytic_data['Morse3500/50+3000/100'] = {
        'zpe': np.array([1750 + 1500, 3212.500]),
        'freqs': np.array([
            [ 3500.000,  3400.000],
            [ 3000.000,  2800.000],
            [ 7000.000,  6700.000],
            [ 6000.000,  5400.000],
            [ 6500.000,  6200.000],
            [10500.000,  9900.000],
            [ 9000.000,  7800.000],
            [10000.000,  9500.000],
            [ 9500.000,  8800.000],
            [14000.000, 13000.000],
            [12000.000, 10000.000],
            [13500.000, 12700.000],
            [12500.000, 11200.000],
            [13000.000, 12100.000]
        ])
    }
    @inactiveTest
    def test_TwoMorseCartesiansNonDeg(self):

        import warnings
        np.seterr(all='raise')
        warnings.filterwarnings('error', category=np.VisibleDeprecationWarning)

        tag = "TwoMorseCartesiansNonDeg"

        # Set up system

        re_1 = 1.0
        re_2 = 1.0
        b_e = np.deg2rad(90.)#(104.5)
        internals = [
            [0, -1, -1, -1],
            [1,  0, -1, -1],
            [2,  0,  1, -1]
        ]

        # raise Exception(
        #                 np.array([
        #         [0.000000, 0.000000, 0.000000],
        #         [re_1, 0.000000, 0.000000],
        #         np.dot(
        #             nput.rotation_matrix([0, 0, 1], b_e),
        #             [re_2, 0.000000, 0.000000]
        #         )
        #     ])
        # )

        # internals = None
        mol = Molecule(
            ["O", "H", "H"],
            np.array([
                [0.000000, 0.000000, 0.000000],
                [re_1, 0.000000, 0.000000],
                np.dot(
                    nput.rotation_matrix([0, 0, 1], b_e),
                    [re_2, 0.000000, 0.000000]
                )
            ]),
            zmatrix=internals,
            guess_bonds=False,
        )

        masses = mol.masses * UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")

        w2h = UnitsData.convert("Wavenumbers", "Hartrees")
        w1 = 3500 * w2h; wx1 = 50 * w2h
        mu_1 = (1 / masses[0] + 1 / masses[1]) ** (-1)
        De_1 = (w1 ** 2) / (4 * wx1)
        a_1 = np.sqrt(2 * mu_1 * wx1)

        w2 = 3000 * w2h; wx2 = 100 * w2h
        mu_2 = (1 / masses[0] + 1 / masses[2]) ** (-1)
        De_2 = (w2 ** 2) / (4 * wx2)
        a_2 = np.sqrt(2 * mu_2 * wx2)

        kb = 0. #1500 * w2h

        morse = self.harmonically_coupled_morse(
            De_1, a_1, re_1,
            De_2, a_2, re_2,
            kb, b_e
        )

        deriv_gen = FiniteDifferenceDerivative(morse,
                                               function_shape=((None, None), 0),
                                               mesh_spacing=1e-3,
                                               stencil=7,
                                               parallelizer=MultiprocessingParallelizer()  # verbose=True)
                                               ).derivatives(mol.coords)

        # raise Exception()

        import time
        start = time.time()
        mol.potential_derivatives = deriv_gen.derivative_tensor([1, 2, 3, 4])#, 5, 6])
        end = time.time()
        print('getting input derivs took {:.3f}s'.format(end - start))
        # raise Exception([x.shape for x in mol.potential_derivatives])

        n_atoms = 3
        n_modes = 3 * n_atoms - 6
        mode_selection = [-2, -1]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(5, n_modes)

        print_report = True
        reference_energies = self.analytic_data['Morse3500/50+3000/100']['zpe']
        reference_freqs = self.analytic_data['Morse3500/50+3000/100']['freqs']

        # def post_wfns_script(wfns, ham):
        #     plt.ArrayPlot(wfns.corrs.hams[2].asarray()).show()
        #     raise Exception(wfns.corrs.hams[2].asarray())

        post_wfns_script = None

        self.run_PT_test(
            tag,
            mol,
            internals,
            mode_selection,
            states,
            reference_energies,
            reference_freqs,
            log=True,
            verbose=True,
            post_wfns_script=post_wfns_script,
            print_report=print_report,
            # ignore_odd_order=False,
            intermediate_normalization=False,
            order=2,
            expansion_order=2
            # , potential_terms=[
            #     np.array([[0.01594717]]),
            #     np.array([[[-0.00808669]]]), #np.array([[[0.]]]),
            #     np.array([[[[0.00318943]]]]), #np.array([[[[0.]]]])
            #     np.array([[[[[-0.00115532]]]]]),
            #     np.array([[[[[[0.000403561]]]]]]),
            #     np.array([[[[[[[-0.000138629]]]]]]]),
            #     np.array([[[[[[[[0.0000472371]]]]]]]]) # for some reason my FD was giving me np.array([[[[[[0.00038929]]]]]])
            # ]
            # , kinetic_terms=[
            #     np.array([[0.01594717]]),
            #     0,
            #     0,
            #     0,
            #     0,
            #     0,
            #     0
            # ]
            , watson=False
            # pseudopotential_terms = [
            #     None,
            #     None,
            #     None,
            #     None,
            #     None
            # ]
            # , parallelized=True
        )

    analytic_data['NedModel'] = {
        'zpe': np.array([1750 + 1500, 3212.500]),
        'freqs': np.array([
            [3500.000, 3400.000],
            [3000.000, 2800.000],
            [7000.000, 6700.000],
            [6000.000, 5400.000],
            [6500.000, 6200.000],
            [10500.000, 9900.000],
            [9000.000, 7800.000],
            [10000.000, 9500.000],
            [9500.000, 8800.000],
            [14000.000, 13000.000],
            [12000.000, 10000.000],
            [13500.000, 12700.000],
            [12500.000, 11200.000],
            [13000.000, 12100.000]
        ])
    }
    @validationTest # was on
    def test_TwoMorseCartesiansAndBendNonDeg(self):

        import warnings
        np.seterr(all='raise')
        warnings.filterwarnings('error', category=np.VisibleDeprecationWarning)

        tag = "TwoMorseCartesiansNed"

        # Set up system

        re_1 = 0.9575
        re_2 = 0.9575
        b_e = np.deg2rad(104.5)
        internals = [
            [0, -1, -1, -1],
            [1,  0, -1, -1],
            [2,  0,  1, -1]
        ]
        # internals = None
        mol = Molecule(
            ["O", "H", "H"],
            np.array([
                [0.000000, 0.000000, 0.000000],
                [re_1, 0.000000, 0.000000],
                np.dot(
                    nput.rotation_matrix([0, 0, 1], b_e),
                    [re_2, 0.000000, 0.000000]
                )
            ]),
            zmatrix=internals,
            guess_bonds=False,
        )

        masses = mol.masses * UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")

        erg2h = UnitsData.convert("Ergs", "Hartrees")
        cm2borh = UnitsData.convert("InverseAngstroms", "InverseBohrRadius")
        De_1 = 8.84e-12 * erg2h
        a_1 = 2.175 * cm2borh

        De_2 = 8.84e-12 * erg2h
        a_2 = 2.175 * cm2borh

        hz2h = UnitsData.convert("Hertz", "Hartrees")
        kb = 7.2916e14/(2 * np.pi) * hz2h

        morse = self.harmonically_coupled_morse(
            De_1, a_1, re_1,
            De_2, a_2, re_2,
            kb, b_e
        )

        deriv_gen = FiniteDifferenceDerivative(morse,
                                               function_shape=((None, None), 0),
                                               mesh_spacing=1e-3,
                                               stencil=7,
                                               parallelizer=MultiprocessingParallelizer()  # verbose=True)
                                               ).derivatives(mol.coords)

        # raise Exception()

        import time
        start = time.time()
        mol.potential_derivatives = deriv_gen.derivative_tensor([1, 2 , 3, 4, 5, 6])
        end = time.time()
        print('getting input derivs took {:.3f}s'.format(end - start))
        # raise Exception([x.shape for x in mol.potential_derivatives])

        # we rigorously zero out small terms for
        # numerical stability
        mat = mol.normal_modes.modes.basis.matrix
        bad_spots = np.where(np.abs(mat) < 1e-14)
        mat[bad_spots] = 0.

        n_atoms = 3
        n_modes = 3 * n_atoms - 6
        mode_selection = [-2, -1]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)

        # states = self.get_states(5, n_modes)
        states = [
            [0, 0],
            [0, 1],
            [1, 0],
            [0, 5],
            [5, 0],
            [1, 4],
            [4, 1],
            [2, 3],
            [3, 2],
        ]

        print_report = True
        # reference_energies = self.analytic_data['NedModel']['zpe']
        # reference_freqs = self.analytic_data['NedModel']['freqs']

        # def post_wfns_script(wfns, ham):
        #     plt.ArrayPlot(wfns.corrs.hams[2].asarray()).show()
        #     raise Exception(wfns.corrs.hams[2].asarray())

        """
        State Energies:
  0 0 3869.908 3827.940        -        - 
  0 1        -        - 3898.075 3727.463 
  1 0        -        - 3841.741 3676.347 
  0 5        -        - 19490.376 18406.154 
  5 0        -        - 19208.703 16919.042 
  1 4        -        - 19434.041 18548.837 
  4 1        -        - 19265.037 16121.717 
  2 3        -        - 19377.707 17704.042 
  3 2        -        - 19321.372 16639.193 
  """

        post_wfns_script = None

        self.run_PT_test(
            tag,
            mol,
            internals,
            mode_selection,
            states,
            None,
            None,
            log=True,
            verbose=True,
            post_wfns_script=post_wfns_script,
            print_report=print_report,
            # ignore_odd_order=False,
            intermediate_normalization=False,
            order=4,
            expansion_order=4
            , kinetic_terms=[
                None,
                None,
                None,
                0,
                0
            ]
            , watson=False
            # pseudopotential_terms = [
            #     None,
            #     None,
            #     None,
            #     None,
            #     None
            # ]
            # , parallelized=True
        )

    #endregion

    #region Water Analogs

    gaussian_data['DOD'] = {
        'zpe': np.array([3406.854, 3366.588]),
        'freqs': np.array([
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
    }
    @validationTest
    def test_DODVPTCartesians(self):

        tag = 'DOD Internals'
        file_name = "DOD_freq.fchk"

        internals = None

        n_atoms = 3
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        print_report = False

        gaussian_energies = self.gaussian_data['DOD']['zpe']
        gaussian_freqs = self.gaussian_data['DOD']['freqs']

        self.run_PT_test(
            tag,
            file_name,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            print_report=print_report
        )

    @validationTest
    def test_DODVPTInternals(self):

        tag = 'DOD Internals'
        file_name = "DOD_freq.fchk"

        internals = [
            [0, -1, -1, -1],
            [1, 0, -1, -1],
            [2, 0, 1, -1]
        ]

        n_atoms = 3
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        print_report = False

        gaussian_energies = self.gaussian_data['DOD']['zpe']
        gaussian_freqs = self.gaussian_data['DOD']['freqs']

        self.run_PT_test(
            tag,
            file_name,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            print_report=print_report
        )

    gaussian_data['HOH'] = {
        'zpe': np.array([4681.564, 4605.953]),
        'freqs': np.array([
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
    }
    @debugTest
    def test_HOHVPTInternals(self):

        tag = 'HOH Internals'
        file_name = "HOH_freq.fchk"

        internals = [
            [0, -1, -1, -1],
            [1,  0, -1, -1],
            [2,  0,  1, -1]
        ]

        n_atoms = 3
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        print_report = True

        gaussian_energies = self.gaussian_data['HOH']['zpe']
        gaussian_freqs = self.gaussian_data['HOH']['freqs']

        """
        State Energies:
          0 0 0 4681.555 4605.955        -        - 
          0 0 1        -        - 3937.525 3744.737 
          0 1 0        -        - 3803.300 3621.995 
          1 0 0        -        - 1622.286 1572.725 
          0 0 2        -        - 7875.049 7391.393 
          0 2 0        -        - 7606.599 7155.885 
          2 0 0        -        - 3244.571 3117.430 
          0 1 1        -        - 7740.824 7200.366 
          1 0 1        -        - 5559.810 5294.413 
          1 1 0        -        - 5425.585 5174.691
        """

        self.run_PT_test(
            tag,
            file_name,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            log=False,
            verbose=False,
            print_report=print_report,
            calculate_intensities=True
        )

    @validationTest
    def test_HOHVPTInternalsEmbedded(self):

        tag = 'HOH Cartesians'
        file_name = "HOH_freq.fchk"

        mol = Molecule.from_file(TestManager.test_data(file_name)).get_embedded_molecule()
        mol.zmatrix  = [
            [0, -1, -1, -1],
            [1,  0, -1, -1],
            [2,  0,  1, -1]
        ]
        internals = []

        n_atoms = 3
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        print_report = False

        gaussian_energies = self.gaussian_data['HOH']['zpe']
        gaussian_freqs = self.gaussian_data['HOH']['freqs']

        self.run_PT_test(
            tag,
            mol,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            log=False,
            verbose=True,
            print_report=print_report,
            calculate_intensities=False
        )

    @validationTest
    def test_HOHVPTCartesians(self):

        tag = 'HOH Cartesians'
        file_name = "HOH_freq.fchk"

        internals = None

        n_atoms = 3
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        print_report = False

        gaussian_energies = self.gaussian_data['HOH']['zpe']
        gaussian_freqs = self.gaussian_data['HOH']['freqs']

        # import McUtils.Misc as mcmisc
        #
        # with mcmisc.without_numba():
        # os.remove(os.path.expanduser('~/hoh.hdf5'))
        self.run_PT_test(
            tag,
            file_name,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            log=True,
            verbose=True,
            print_report=print_report,
            calculate_intensities=True,
            print_x=True
            # , chunk_size=200
            # zero_order_energy_corrections = [
            #     [(0, 1, 0), 5500 * UnitsData.convert("Wavenumbers", "Hartrees")]
            # ],
            # , memory_constrained=True
            # , state_space_terms=((1, 0), (2, 0))
            # , checkpoint=os.path.expanduser('~/hoh.hdf5'),
            # , watson=False
        )

    @validationTest
    def test_HOHVPTCartesiansFiltered(self):

        import warnings
        np.seterr(all='raise')
        warnings.filterwarnings('error', category=np.VisibleDeprecationWarning)

        tag = 'HOH Cartesians'
        file_name = "HOH_freq.fchk"

        internals = None

        n_atoms = 3
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        print_report = False

        gaussian_energies = self.gaussian_data['HOH']['zpe']
        gaussian_freqs = self.gaussian_data['HOH']['freqs']

        # import McUtils.Misc as mcmisc
        #
        # with mcmisc.without_numba():
        # os.remove(os.path.expanduser('~/hoh.hdf5'))
        self.run_PT_test(
            tag,
            file_name,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            log=True,
            verbose=True,
            print_report=print_report,
            calculate_intensities=True,
            state_space_filters={
                # (2, 0): BasisStateSpace(
                #     HarmonicOscillatorProductBasis(n_modes),
                #     states[:1],
                # ).apply_selection_rules([[]]), # get energies right
                (1, 1): BasisStateSpace(
                    HarmonicOscillatorProductBasis(n_modes),
                    states[:1]
                ).apply_selection_rules([[-1], [], [1]])

            }
            # target_property_rules=([0], [
            #     [-1,  0, 1],
            #     [-2,  0, 2],
            #     [-3, -1, 1, 3],
            # ])
            # zero_order_energy_corrections = [
            #     [(0, 1, 0), 5500 * UnitsData.convert("Wavenumbers", "Hartrees")]
            # ],
            # , memory_constrained=True
            # , state_space_terms=((1, 0), (2, 0))
            # , checkpoint=os.path.expanduser('~/hoh.hdf5'),
            # , watson=False
        )

    @validationTest
    def test_HOHVPTCartesiansEmbdeddd(self):

        tag = 'HOH Cartesians'
        file_name = "HOH_freq.fchk"

        mol = Molecule.from_file(TestManager.test_data(file_name)).get_embedded_molecule()
        # mol.source_file = None
        # mol.zmatrix = [
        #     [0, -1, -1, -1],
        #     [1,  0, -1, -1],
        #     [2,  0, 1, -1]
        # ]
        internals = None

        n_atoms = 3
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        print_report = False

        gaussian_energies = self.gaussian_data['HOH']['zpe']
        gaussian_freqs = self.gaussian_data['HOH']['freqs']

        # def floop(ham, states):
        #     raise Exception(*args)
        self.run_PT_test(
            tag,
            mol,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            print_x=True,
            print_report=print_report,
            log=True
        )

    @validationTest
    def test_HOHVPTCartesians4thOrder(self):

        import warnings
        np.seterr(all='raise')
        warnings.filterwarnings('error', category=np.VisibleDeprecationWarning)

        tag = 'HOH Cartesians'
        file_name = "HOH_freq.fchk"

        internals = None

        n_atoms = 3
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        print_report = True

        gaussian_energies = self.gaussian_data['HOH']['zpe']
        gaussian_freqs = self.gaussian_data['HOH']['freqs']

        # import McUtils.Misc as mcmisc
        #
        # with mcmisc.without_numba():
        # os.remove(os.path.expanduser('~/hoh.hdf5'))
        self.run_PT_test(
            tag,
            file_name,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            log=True,
            verbose=True,
            order=4,
            expansion_order=2,
            print_report=print_report,
            calculate_intensities=True
            # , checkpoint=os.path.expanduser('~/hoh.hdf5'),
            # , watson=False
        )

    @validationTest
    def test_HOHVPTInternals4thOrder(self):

        tag = 'HOH Internals'
        file_name = "HOH_freq.fchk"

        internals = [
            [0, -1, -1, -1],
            [1, 0, -1, -1],
            [2, 0, 1, -1]
        ]

        n_atoms = 3
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        print_report = True

        gaussian_energies = self.gaussian_data['HOH']['zpe']
        gaussian_freqs = self.gaussian_data['HOH']['freqs']

        self.run_PT_test(
            tag,
            file_name,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            log=True,
            verbose=True,
            order=4,
            expansion_order=2,
            print_report=print_report,
            calculate_intensities=True
        )

    @inactiveTest
    def test_HOTVPTCartesians4thOrder(self):

        tag = 'HOT Cartesians'
        file_name = "HOT_freq.fchk"

        internals = [
            [0, -1, -1, -1],
            [1, 0, -1, -1],
            [2, 0, 1, -1]
        ]
        internals = None

        n_atoms = 3
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        basis = HarmonicOscillatorProductBasis(n_modes)
        coupled_states = [
            # BasisStateSpace(basis, self.get_states(20, n_modes)).apply_selection_rules(
            #     basis.selection_rules("x", "x", "x")
            # ),
            self.get_states(15, n_modes),
            self.get_states(15, n_modes)
            # BasisStateSpace(basis, self.get_states(20, n_modes)).apply_selection_rules(
            #     basis.selection_rules("x", "x", "x", "x")
            # )
        ]
        # coupled_states=None

        print_report = True
        log = True
        nielsen_tolerance = None
        order = 4

        gaussian_energies = self.gaussian_data['HOT']['zpe']
        gaussian_freqs = self.gaussian_data['HOT']['freqs']

        self.run_PT_test(
            tag,
            file_name,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            order=order,
            nielsen_tolerance=nielsen_tolerance,
            log=log,
            print_report=print_report
            # , energy_order=2
            , coupled_states=coupled_states
            , intermediate_normalization=True
        )

    @inactiveTest
    def test_HOHVPTCartesiansDegenerate(self):

        internals = None

        basis = HarmonicOscillatorProductBasis(3)

        states = self.get_states(6, 3)

        wfns = self.get_VPT2_wfns(
            "HOH_freq.fchk",
            internals,
            states,
            regenerate=True,
            log=True
            # , degeneracies = {
            #     'NT':(1, 2, 2),
            #     'MartinTest': False
            # }
            # , degeneracies=.0005 # anything within 1 Hartree (i.e. everything)
            # , degeneracies=(1, 2, 2) # use N_T = 1 bend + 2 stretch
            , degeneracies=(
                # [[0, 0, 1], [0, 1, 0], [2, 0, 0]],
                [[0, 0, 2], [0, 2, 0], [0, 1, 1]],
            )
            , allow_sakurai_degs=True
            , allow_post_PT_calc=True
            # , order=4
            # degeneracies=100/self.h2w # any pair of states within 100 wavenumbers can be treated as degenerate
            # , coupled_states = [
            #     BasisStateSpace(basis, self.get_states(12, 3)).apply_selection_rules(
            #         basis.selection_rules("x", "x", "x")
            #     ),
            #     BasisStateSpace(basis, self.get_states(12, 3)).apply_selection_rules(
            #         basis.selection_rules("x", "x", "x", "x")
            #     )
            # ]
            # , coupled_states=[
            #     BasisStateSpace(basis, self.get_states(12, 3)),
            #     BasisStateSpace(basis, self.get_states(12, 3))
            # ]
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

        ns = min(len(my_freqs), len(gaussian_freqs))
        print_report = False
        if print_report:
            self.print_energy_block("Gaussian Energies:", states, gaussian_energies, gaussian_freqs[:ns])
            self.print_energy_block("State Energies:", states, my_energies, my_freqs[:ns])

        print_diffs = True
        if print_diffs:
            self.print_energy_block("Difference Energies:", states, my_energies - gaussian_energies,
                                    my_freqs[:ns] - gaussian_freqs[:ns])

        self.assertLess(np.max(np.abs(my_freqs[:ns] - gaussian_freqs[:ns])), 1.5)

    @inactiveTest
    def test_HOHVPTInternalsDegenerate(self):

        internals = [
            [0, -1, -1, -1],
            [1,  0, -1, -1],
            [2,  0,  1, -1]
        ]

        states = self.get_states(3, 3)

        wfns = self.get_VPT2_wfns(
            "HOH_freq.fchk",
            internals,
            states,
            regenerate=True,
            log=True
            # , degeneracies = {
            #     'NT':(1, 2, 2),
            #     'MartinTest': False
            # }
            # , degeneracies=1  # anything within 1 Hartree (i.e. everything)
            , degeneracies=(1, 2, 2) # use N_T = 1 bend + 2 stretch
            # degeneracies=100/self.h2w # any pair of states within 100 wavenumbers can be treated as degenerate
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

        print_report = True
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
        self.assertLess(np.max(np.abs(my_freqs - gaussian_freqs[:len(my_freqs)])), 100)

    gaussian_data['HOD'] = {
        'zpe': np.array([4052.912, 3994.844]),
        'freqs': np.array([
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
    }
    @validationTest
    def test_HODVPTInternals(self):

        tag = 'HOD Internals'
        file_name = "HOD_freq.fchk"

        internals = [
            [0, -1, -1, -1],
            [1,  0, -1, -1],
            [2,  0,  1, -1]
        ]

        n_atoms = 3
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        print_report = False

        gaussian_energies = self.gaussian_data['HOD']['zpe']
        gaussian_freqs = self.gaussian_data['HOD']['freqs']

        self.run_PT_test(
            tag,
            file_name,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            print_report=print_report
        )

    @validationTest
    def test_HODVPTCartesians(self):

        tag = 'HOD Cartesians'
        file_name = "HOD_freq.fchk"

        internals = None

        n_atoms = 3
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        print_report = False

        gaussian_energies = self.gaussian_data['HOD']['zpe']
        gaussian_freqs = self.gaussian_data['HOD']['freqs']

        self.run_PT_test(
            tag,
            file_name,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            print_report=print_report
        )

    gaussian_data['HOT'] = {
        'zpe': np.array([3789.117, 3736.855]),
        'freqs': np.array([
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
    }
    @validationTest
    def test_HOTVPTInternals(self):

        tag = 'HOT Cartesians'
        file_name = "HOT_freq.fchk"

        internals = [
            [0, -1, -1, -1],
            [1, 0, -1, -1],
            [2, 0, 1, -1]
        ]

        n_atoms = 3
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        print_report = False

        gaussian_energies = self.gaussian_data['HOT']['zpe']
        gaussian_freqs = self.gaussian_data['HOT']['freqs']

        self.run_PT_test(
            tag,
            file_name,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            print_report=print_report
        )

    @validationTest
    def test_HOTVPTCartesians(self):

        tag = 'HOT Cartesians'
        file_name = "HOT_freq.fchk"

        internals = None

        n_atoms = 3
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        print_report = False

        gaussian_energies = self.gaussian_data['HOT']['zpe']
        gaussian_freqs = self.gaussian_data['HOT']['freqs']

        self.run_PT_test(
            tag,
            file_name,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            print_report=print_report
        )

    #endregion Water Analogs

    #region Formaldehyde Analogs

    gaussian_data['OCHH'] = {
        'zpe': [5866.872, 5785.953],
        'freqs': np.array([
            [3061.701, 2849.454],
            [2977.640, 2820.563],
            [1727.083, 1695.163],
            [1527.041, 1493.139],
            [1252.164, 1231.359],
            [1188.114, 1166.697],
            # 2 quanta
            [6123.403, 5719.555],
            [5955.281, 5578.213],
            [3454.165, 3372.197],
            [3054.082, 2988.763],
            [2504.328, 2459.646],
            [2376.227, 2327.621],
            # mixed states
            [6039.342, 5586.305],
            [4788.784, 4589.165],
            [4588.742, 4374.903],
            [4313.865, 4114.701],
            [4249.815, 4043.462],
            [4704.723, 4510.744],
            [4504.681, 4278.834],
            [4229.804, 4043.264],
            [4165.754, 3978.440],
            [3254.123, 3181.610],
            [2979.247, 2967.726],
            [2915.196, 2854.766],
            [2779.205, 2710.416],
            [2715.155, 2657.682],
            [2440.278, 2404.805]
        ])
    }
    @validationTest
    def test_OCHHVPTInternals(self):

        tag = 'OCHH Internals'
        file_name = "OCHH_freq.fchk"

        internals = [
            [0, -1, -1, -1],
            [1,  0, -1, -1],
            [2,  1,  0, -1],
            [3,  1,  0,  2]
        ]

        n_atoms = 4
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        gaussian_energies = self.gaussian_data['OCHH']['zpe']
        gaussian_freqs = self.gaussian_data['OCHH']['freqs']

        # ExpansionTerms.cartesian_analytic_deriv_order = 0
        # ExpansionTerms.cartesian_fd_mesh_spacing = 1.0e-3
        # ExpansionTerms.cartesian_fd_stencil = 9

        print_report = False
        nielsen_tolerance = 1
        gaussian_tolerance = 100
        self.run_PT_test(
            tag,
            file_name,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            log=False,
            verbose=False,
            calculate_intensities=True,
            print_report=print_report,
            nielsen_tolerance=nielsen_tolerance,
            gaussian_tolerance=gaussian_tolerance
        )

    @validationTest
    def test_OCHHVPTCartesians(self):

        tag = 'OCHH Cartesians'
        file_name = "OCHH_freq.fchk"

        internals = None

        n_atoms = 4
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        gaussian_energies = self.gaussian_data['OCHH']['zpe']
        gaussian_freqs = self.gaussian_data['OCHH']['freqs']

        print_report = False
        nielsen_tolerance = 1
        gaussian_tolerance = 100
        self.run_PT_test(
            tag,
            file_name,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            log=False,
            verbose=False,
            calculate_intensities=True,
            print_report=print_report,
            nielsen_tolerance=nielsen_tolerance,
            gaussian_tolerance=gaussian_tolerance
        )

    @validationTest
    def test_OCHHVPTCartesiansFiltered(self):

        tag = 'OCHH Cartesians'
        file_name = "OCHH_freq.fchk"

        internals = None

        n_atoms = 4
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        gaussian_energies = self.gaussian_data['OCHH']['zpe']
        gaussian_freqs = self.gaussian_data['OCHH']['freqs']

        print_report = False
        nielsen_tolerance = 1
        gaussian_tolerance = 100
        self.run_PT_test(
            tag,
            file_name,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            log=True,
            verbose=True,
            calculate_intensities=True,
            print_report=print_report,
            nielsen_tolerance=nielsen_tolerance,
            gaussian_tolerance=gaussian_tolerance,
            state_space_filters={
                # (2, 0): BasisStateSpace(
                #     HarmonicOscillatorProductBasis(n_modes),
                #     states[:1],
                # ).apply_selection_rules([[]]), # get energies right
                (1, 1): BasisStateSpace(
                    HarmonicOscillatorProductBasis(n_modes),
                    states[:1]
                ).apply_selection_rules([[-1], [], [1]])

            }
        )

    @validationTest
    def test_OCHHVPTCartesiansDegenerate(self):

        tag = 'OCHH Cartesians'
        file_name = "OCHH_freq.fchk"

        internals = None

        n_atoms = 4
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        # nt_spec = np.array([
        #     [0, 0, 0, 0, 0, 1],
        #     [0, 1, 0, 1, 0, 0]
        # ])
        # degeneracies = []
        # for s in states:
        #     if np.dot(nt_spec[0], s) > 0:
        #         degeneracies.append([
        #             s,
        #             s - nt_spec[0] + nt_spec[1]
        #         ])
        # # raise Exception(degeneracies)

        from McUtils.GaussianInterface import GaussianLogReader
        with GaussianLogReader(TestManager.test_data('OCHH_freq_16.log')) as reader:
            sort_spec = np.flip([4, 0, 1, 2, 5, 3])
            x = reader.parse("XMatrix")["XMatrix"][np.ix_(sort_spec, sort_spec)]
            # raise Exception(x)
        def pre_wfns_script(hammer, states):
            harm, corrs = hammer.get_Nielsen_energies(states, x_mat=x)
            harm = harm * self.h2w
            anharms = harm + corrs
            # print(states)
            print(np.column_stack([
                harm[1:] - harm[0],
                anharms[1:] - anharms[0]
            ]))
        pre_wfns_script = None

        degeneracies = self.get_degenerate_polyad_space(
            states,
            [
                [[0, 0, 0, 0, 0, 1], [0, 1, 0, 1, 0, 0]]
            ]
        )
        degeneracies = [np.array(x).tolist() for x in degeneracies]
        states = np.array(states).tolist()
        for pair in degeneracies:
            for p in pair:
                if p not in states:
                    states.append(p)

        # degeneracies = None

        gaussian_energies = self.gaussian_data['OCHH']['zpe']
        gaussian_freqs = self.gaussian_data['OCHH']['freqs']

        print_report = False
        nielsen_tolerance = None
        gaussian_tolerance = 1
        self.run_PT_test(
            tag,
            file_name,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            print_report=print_report,
            nielsen_tolerance=nielsen_tolerance,
            pre_wfns_script=pre_wfns_script,
            log=True,
            verbose=True,
            calculate_intensities=True,
            gaussian_tolerance=gaussian_tolerance,
            degeneracies=degeneracies
            # , allow_post_PT_calc=False
            # , invert_x=True
            # , modify_degenerate_perturbations=True
        )

    @validationTest
    def test_OCHHVPTInternalsDegenerate(self):

        tag = 'OCHH Cartesians'
        file_name = "OCHH_freq.fchk"

        internals = [
            [0, -1, -1, -1],
            [1,  0, -1, -1],
            [2,  1,  0, -1],
            [3,  1,  0,  2]
        ]

        n_atoms = 4
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        # nt_spec = np.array([
        #     [0, 0, 0, 0, 0, 1],
        #     [0, 1, 0, 1, 0, 0]
        # ])
        # degeneracies = []
        # for s in states:
        #     if np.dot(nt_spec[0], s) > 0:
        #         degeneracies.append([
        #             s,
        #             s - nt_spec[0] + nt_spec[1]
        #         ])
        # # raise Exception(degeneracies)

        from McUtils.GaussianInterface import GaussianLogReader
        with GaussianLogReader(TestManager.test_data('OCHH_freq_16.log')) as reader:
            sort_spec = np.flip([4, 0, 1, 2, 5, 3])
            x = reader.parse("XMatrix")["XMatrix"][np.ix_(sort_spec, sort_spec)]
            # raise Exception(x)

        def pre_wfns_script(hammer, states):
            harm, corrs = hammer.get_Nielsen_energies(states, x_mat=x)
            harm = harm * self.h2w
            anharms = harm + corrs
            # print(states)
            print(np.column_stack([
                harm[1:] - harm[0],
                anharms[1:] - anharms[0]
            ]))

        pre_wfns_script = None

        degeneracies = self.get_degenerate_polyad_space(
            states,
            [
                [[0, 0, 0, 0, 0, 1], [0, 1, 0, 1, 0, 0]]
            ]
        )
        degeneracies = [np.array(x).tolist() for x in degeneracies]
        states = np.array(states).tolist()
        for pair in degeneracies:
            for p in pair:
                if p not in states:
                    states.append(p)

        # degeneracies = None

        gaussian_energies = self.gaussian_data['OCHH']['zpe']
        gaussian_freqs = self.gaussian_data['OCHH']['freqs']

        print_report = True
        nielsen_tolerance = None
        gaussian_tolerance = 1
        self.run_PT_test(
            tag,
            file_name,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            print_report=print_report,
            nielsen_tolerance=nielsen_tolerance,
            pre_wfns_script=pre_wfns_script,
            log=False,
            verbose=False,
            calculate_intensities=True,
            gaussian_tolerance=gaussian_tolerance,
            degeneracies=degeneracies
            , allow_post_PT_calc=False
            # , invert_x=True
            # , modify_degenerate_perturbations=True
        )

    gaussian_data['OCHD'] = {
        'zpe': [5235.162,  5171.477],
        'freqs': np.array([
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
            [4440.352, 4192.841],
            [4099.287, 3896.226],
            [4053.494, 3866.190],
            [3922.817, 3809.397],
            [3638.807, 3518.724],
            [3297.743, 3183.261],
            [3251.950, 3146.116],
            [3119.087, 3060.953],
            [2778.023, 2725.597],
            [2732.230, 2681.778],
            [2494.013, 2448.626],
            [2448.220, 2403.196],
            [2107.156, 2074.925]
        ])
    }
    @validationTest
    def test_OCHDVPTInternals(self):

        tag = 'OCHD Internals'
        file_name = "OCHD_freq.fchk"

        internals = [
            [0, -1, -1, -1],
            [1, 0, -1, -1],
            [2, 1, 0, -1],
            [3, 1, 0, 2]
        ]

        n_atoms = 4
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        gaussian_energies = self.gaussian_data['OCHD']['zpe']
        gaussian_freqs = self.gaussian_data['OCHD']['freqs']

        print_report = False
        nielsen_tolerance = 1
        gaussian_tolerance = 50
        self.run_PT_test(
            tag,
            file_name,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            print_report=print_report,
            nielsen_tolerance=nielsen_tolerance,
            gaussian_tolerance=gaussian_tolerance
        )

    @validationTest
    def test_OCHDVPTCartesians(self):

        tag = 'OCHD Cartesians'
        file_name = "OCHD_freq.fchk"

        internals = None

        n_atoms = 4
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        gaussian_energies = self.gaussian_data['OCHD']['zpe']
        gaussian_freqs = self.gaussian_data['OCHD']['freqs']

        print_report = False
        nielsen_tolerance = 1
        gaussian_tolerance = 50
        self.run_PT_test(
            tag,
            file_name,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            print_report=print_report,
            nielsen_tolerance=nielsen_tolerance,
            gaussian_tolerance=gaussian_tolerance
        )

    gaussian_data['OCHT'] = {
        'zpe': [4970.730, 4912.805],
        'freqs': np.array([
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
            [4421.438,  4220.692],
            [4058.812,  3870.413],
            [3926.768,  3756.997],
            [3578.729,  3474.240],
            [3317.774,  3227.798],
            [2955.148,  2867.701],
            [2823.105,  2746.008],
            [3059.543,  3000.454],
            [2696.917,  2644.455],
            [2564.874,  2515.823],
            [2435.962,  2393.534],
            [2303.919,  2263.695],
            [1941.293,  1911.106]
        ])
    }
    @validationTest
    def test_OCHTVPTInternals(self):

        tag = 'OCHT Internals'
        file_name = "OCHT_freq.fchk"

        internals = [
            [0, -1, -1, -1],
            [1,  0, -1, -1],
            [2,  1,  0, -1],
            [3,  1,  0,  2]
        ]

        n_atoms = 4
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        gaussian_energies = self.gaussian_data['OCHT']['zpe']
        gaussian_freqs = self.gaussian_data['OCHT']['freqs']

        print_report = False
        nielsen_tolerance = 1
        gaussian_tolerance = 50
        self.run_PT_test(
            tag,
            file_name,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            print_report=print_report,
            nielsen_tolerance=nielsen_tolerance,
            gaussian_tolerance=gaussian_tolerance
        )

    @validationTest
    def test_OCHTVPTCartesians(self):

        tag = 'OCHT Cartesians'
        file_name = "OCHT_freq.fchk"

        internals = None

        n_atoms = 4
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        gaussian_energies = self.gaussian_data['OCHT']['zpe']
        gaussian_freqs = self.gaussian_data['OCHT']['freqs']

        print_report = False
        nielsen_tolerance=1
        gaussian_tolerance=50
        self.run_PT_test(
            tag,
            file_name,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            print_report=print_report,
            nielsen_tolerance=nielsen_tolerance,
            gaussian_tolerance=gaussian_tolerance
            , log=True
            # , print_profile=True
        )

    @inactiveTest
    def test_OCHTVPTCartesiansDegenerate(self):

        tag = 'OCHT Cartesians'
        file_name = "OCHT_freq.fchk"

        internals = None

        n_atoms = 4
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        gaussian_energies = self.gaussian_data['OCHT']['zpe']
        gaussian_freqs = self.gaussian_data['OCHT']['freqs']

        print_report = False
        nielsen_tolerance=1
        gaussian_tolerance=50
        self.run_PT_test(
            tag,
            file_name,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            print_report=print_report,
            nielsen_tolerance=nielsen_tolerance,
            gaussian_tolerance=gaussian_tolerance
            , log=True
            , degeneracies=degeneracies
            # , print_profile=True
        )

    #endregion Formaldehyde Analogs

    #region X-HOCL



    #Paper
    @validationTest
    def test_ClHOClVPTInternals(self):

        tag = 'ClHOCl Internals'
        file_name = TestManager.test_data("cl_hocl.fchk")

        n_atoms = 4
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        # internals = [
        #     [0, -1, -1, -1],  # H
        #     [1,  0, -1, -1],  # O
        #     [2,  1,  0, -1],  # O
        #     [3,  2,  1,  0],  # N
        #     [4,  3,  2,  0]   # O
        # ]
        internals = [
            [1, -1, -1, -1],
            [2,  1, -1, -1],
            [3,  2,  1, -1],
            [0,  1,  2,  3]
        ]

        gaussian_energies = None#self.gaussian_data['HOONO']['zpe']
        gaussian_freqs = None#self.gaussian_data['HOONO']['freqs']

        print_report = True
        nielsen_tolerance = 10
        gaussian_tolerance = 10
        # from Psience.VPT2 import PotentialTerms
        # PotentialTerms.hessian_tolerance = None
        self.run_PT_test(
            tag,
            file_name,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            log=True,
            verbose=True,
            calculate_intensities=True,
            print_report=print_report,
            nielsen_tolerance=nielsen_tolerance,
            gaussian_tolerance=gaussian_tolerance
        )
    @inactiveTest
    def test_ClHOClVPTInternalsDummy(self):

        tag = 'Cl-HOCl Internals'
        file_name = TestManager.test_data("cl_hocl.fchk")

        mol = Molecule.from_file(file_name).get_embedded_molecule()
        n_pos = mol.atom_positions["N"]
        o_pos = mol.atom_positions["O"]

        normal = -nput.vec_crosses(
            mol.coords[o_pos[0]] - mol.coords[o_pos[1]],
            mol.coords[n_pos[0]] - mol.coords[o_pos[1]],
            normalize=True
        )

        mol = mol.insert_atoms("X", mol.coords[o_pos[1]] + 5 * normal, 3, handle_properties=False)

        # internals = [
        #     [0, -1, -1, -1],  # H
        #     [1,  0, -1, -1],  # O
        #     [2,  1,  0, -1],  # O
        #     [3,  2,  1,  0],  # N
        #     [4,  3,  2,  0]  # O
        # ]
        # mol.zmatrix = [
        #     [0, -1, -1, -1],  # H
        #     [1,  0, -1, -1],  # O
        #     [2,  1,  0, -1],  # O
        #     [3,  2,  1,  0],  # X
        #     [4,  2,  1,  3],  # N
        #     [5,  4,  2,  3]   # O
        # ]
        mol.zmatrix = [
            [1, -1, -1, -1],  # O
            [2,  1, -1, -1],  # O
            [4,  2,  1, -1],  # N
            [3,  2,  1,  4],  # X
            [0,  1,  2,  3],  # H
            [5,  4,  2,  3]   # O
        ]
        # mol.zmatrix = [
        #     [1, -1, -1, -1],  # O
        #     [2,  1, -1, -1],  # O
        #     [4,  2,  1, -1],  # N
        #     [3,  2,  1,  4],  # X
        #     [0,  1,  2,  3],  # H
        #     [5,  4,  2,  3]   # O
        # ]

        # raise Exception(mol.internal_coordinates)
        # mol.zmatrix = [
        #     [1, -1, -1, -1],  # O
        #     [2,  1, -1, -1],  # O
        #     [4,  2,  1, -1],  # N
        #     [3,  2,  1,  4],  # X
        #     [0,  1,  2,  3],  # H
        #     [5,  4,  2,  3]  # O
        # ]

        n_atoms = 5
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        gaussian_energies = self.gaussian_data['HOONO']['zpe']
        gaussian_freqs = self.gaussian_data['HOONO']['freqs']

        print_report = True
        nielsen_tolerance = 10
        gaussian_tolerance = 10

        PotentialTerms.hessian_tolerance = None
        internals = None
        self.run_PT_test(
            tag,
            mol,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            log=True,
            verbose=True,
            print_report=print_report,
            nielsen_tolerance=nielsen_tolerance,
            gaussian_tolerance=gaussian_tolerance
        )
    @validationTest
    def test_ClHOClVPTCartesians(self):

        tag = 'Cl-HOCl Cartesians'
        file_name = "cl_hocl.fchk"

        internals = None

        n_atoms = 4
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        gaussian_energies = None#self.gaussian_data['HOONO']['zpe']
        gaussian_freqs = None#self.gaussian_data['HOONO']['freqs']
        pre_wfns_script = None

        degeneracies = None
        print_report = False
        nielsen_tolerance = 10 if degeneracies is None else 500
        gaussian_tolerance = 10 if degeneracies is not None else 50
        self.run_PT_test(
            tag,
            file_name,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            log=True,
            verbose=True,
            calculate_intensities=True,
            degeneracies=degeneracies,
            print_report=print_report,
            nielsen_tolerance=nielsen_tolerance,
            gaussian_tolerance=gaussian_tolerance,
            pre_wfns_script=pre_wfns_script,
            gaussian_resonance_handling=True
            # zero_element_warning=False

        )

    #endregion

    #region HOONO

    gaussian_data['HOONO'] = {
        'zpe': [5144.56105, 5052.05930],
        'freqs': np.array([
            [3486.855,  3299.146],
            [1568.502,  1558.272],
            [1433.609,  1395.836],
            [ 970.877,   947.362],
            [ 836.395,   814.545],
            [ 715.883,   799.170],
            [ 524.099,   482.999],
            [ 397.168,   371.575],
            [ 355.733,   270.102],
            # 2 quanta
            [6973.711,  6378.861],
            [3137.004,  3097.108],
            [2867.219,  2783.602],
            [1941.754,  1889.021],
            [1672.790,  1622.255],
            [1431.767,  1543.473],
            [1048.199,   949.341],
            [ 794.335,   736.075],
            [ 711.467,   309.620],
            # mixed states
            [5055.357,  4854.609],
            [4920.465,  4656.617],
            [4457.732,  4250.173],
            [4323.250,  4118.829],
            [4202.739,  4110.476],
            [4010.955,  3804.739],
            [3884.023,  3672.896],
            [3842.589,  3629.964],
            [3002.112,  2946.860],
            [2539.379,  2502.223],
            [2404.897,  2371.654],
            [2284.385,  2386.161],
            [2092.601,  2042.056],
            [1965.670,  1930.208],
            [1924.236,  1833.376],
            [2404.486,  2340.072],
            [2270.004,  2203.572],
            [2149.493,  2194.956],
            [1957.709,  1884.793],
            [1830.777,  1768.335],
            [1789.343,  1680.838],
            [1807.272,  1752.785],
            [1686.760,  1737.147],
            [1494.976,  1428.035],
            [1368.044,  1314.186],
            [1326.610,  1212.119],
            [1552.278,  1609.052],
            [1360.494,  1292.706],
            [1233.562,  1182.927],
            [1192.128,  1081.094],
            [1239.983,  1265.933],
            [1113.051,  1161.080],
            [1071.617,  1350.184],
            [ 921.267,   848.123],
            [ 879.833,   706.552],
            [ 752.901,   627.196],
        ])
    }
    #Paper
    @validationTest
    def test_HOONOVPTInternals(self):

        tag = 'HOONO Internals'
        file_name = TestManager.test_data("HOONO_freq.fchk")

        n_atoms = 5
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        # internals = [
        #     [0, -1, -1, -1],  # H
        #     [1,  0, -1, -1],  # O
        #     [2,  1,  0, -1],  # O
        #     [3,  2,  1,  0],  # N
        #     [4,  3,  2,  0]   # O
        # ]
        internals = [
            [1, -1, -1, -1],
            [2,  1, -1, -1],
            [3,  2,  1, -1],
            [0,  1,  2,  3],
            [4,  3,  2,  1]
        ]

        gaussian_energies = self.gaussian_data['HOONO']['zpe']
        gaussian_freqs = self.gaussian_data['HOONO']['freqs']

        print_report = True
        nielsen_tolerance = 10
        gaussian_tolerance = 10
        # from Psience.VPT2 import PotentialTerms
        # PotentialTerms.hessian_tolerance = None
        ExpansionTerms.cartesian_analytic_deriv_order = 0
        self.run_PT_test(
            tag,
            file_name,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            log=True,
            verbose=True,
            calculate_intensities=True,
            print_report=print_report,
            nielsen_tolerance=nielsen_tolerance,
            gaussian_tolerance=gaussian_tolerance
        )
    @inactiveTest
    def test_HOONOVPTInternalsDummy(self):

        tag = 'HOONO Internals'
        file_name = TestManager.test_data("HOONO_freq.fchk")

        mol = Molecule.from_file(file_name)#.get_embedded_molecule()
        n_pos = mol.atom_positions["N"]
        o_pos = mol.atom_positions["O"]

        normal = -nput.vec_crosses(
            mol.coords[o_pos[0]] - mol.coords[o_pos[1]],
            mol.coords[n_pos[0]] - mol.coords[o_pos[1]],
            normalize=True
        )

        mol = mol.insert_atoms("X", mol.coords[o_pos[1]] + 5 * normal, 3, handle_properties=False)

        # internals = [
        #     [0, -1, -1, -1],  # H
        #     [1,  0, -1, -1],  # O
        #     [2,  1,  0, -1],  # O
        #     [3,  2,  1,  0],  # N
        #     [4,  3,  2,  0]  # O
        # ]
        # mol.zmatrix = [
        #     [0, -1, -1, -1],  # H
        #     [1,  0, -1, -1],  # O
        #     [2,  1,  0, -1],  # O
        #     [3,  2,  1,  0],  # X
        #     [4,  2,  1,  3],  # N
        #     [5,  4,  2,  3]   # O
        # ]
        mol.zmatrix = [
            [1, -1, -1, -1],  # O
            [2,  1, -1, -1],  # O
            [4,  2,  1, -1],  # N
            [3,  2,  1,  4],  # X
            [0,  1,  2,  3],  # H
            [5,  4,  2,  3]   # O
        ]
        # mol.zmatrix = [
        #     [1, -1, -1, -1],  # O
        #     [2,  1, -1, -1],  # O
        #     [4,  2,  1, -1],  # N
        #     [3,  2,  1,  4],  # X
        #     [0,  1,  2,  3],  # H
        #     [5,  4,  2,  3]   # O
        # ]

        # raise Exception(mol.internal_coordinates)
        # mol.zmatrix = [
        #     [1, -1, -1, -1],  # O
        #     [2,  1, -1, -1],  # O
        #     [4,  2,  1, -1],  # N
        #     [3,  2,  1,  4],  # X
        #     [0,  1,  2,  3],  # H
        #     [5,  4,  2,  3]  # O
        # ]

        # raise Exception(mol.normal_modes.modes.freqs*UnitsData.convert("Hartrees", "Wavenumbers"))

        n_atoms = 5
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        gaussian_energies = self.gaussian_data['HOONO']['zpe']
        gaussian_freqs = self.gaussian_data['HOONO']['freqs']

        print_report = True
        nielsen_tolerance = 10
        gaussian_tolerance = 10

        ExpansionTerms.cartesian_analytic_deriv_order = 0
        internals = None
        self.run_PT_test(
            tag,
            mol,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            log=True,
            verbose=True,
            print_report=print_report,
            nielsen_tolerance=nielsen_tolerance,
            gaussian_tolerance=gaussian_tolerance
        )
    @validationTest
    def test_HOONOVPTInternalsEmbed(self):

        tag = 'HOONO Internals'
        file_name = TestManager.test_data("HOONO_freq.fchk")

        mol = Molecule.from_file(TestManager.test_data(file_name)).get_embedded_molecule()
        mol.zmatrix = [
            [1, -1, -1, -1],
            [2,  1, -1, -1],
            [3,  2,  1, -1],
            [0,  1,  2,  3],
            [4,  3,  2,  1]
        ]
        internals = []

        n_atoms = 5
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        gaussian_energies = self.gaussian_data['HOONO']['zpe']
        gaussian_freqs = self.gaussian_data['HOONO']['freqs']

        print_report = True
        nielsen_tolerance = 10
        gaussian_tolerance = 10
        # from Psience.VPT2 import PotentialTerms
        # PotentialTerms.hessian_tolerance = None
        self.run_PT_test(
            tag,
            mol,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            log=True,
            verbose=True,
            print_report=print_report,
            nielsen_tolerance=nielsen_tolerance,
            gaussian_tolerance=gaussian_tolerance
        )
    @validationTest
    def test_HOONOVPTCartesians(self):

        tag = 'HOONO Cartesians'
        file_name = "HOONO_freq.fchk"

        internals = None

        n_atoms = 5
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        gaussian_energies = self.gaussian_data['HOONO']['zpe']
        gaussian_freqs = self.gaussian_data['HOONO']['freqs']

        # """
        # 0.337366D+01  0.708754D+00  0.685684D+00  0.704381D+00  0.386143D-01 0.189447D+00  0.238868D-01  0.000000D+00  0.000000D+00
        # """

        # def pre_wfns_script(hammer, states):
        #     harm, corrs, x = hammer.get_Nielsen_energies(states, return_split=True)
        #     x_flips = np.argsort([8, 6, 7, 5, 4, 3, 2, 1, 0])
        #     _ = []
        #     for y in x:
        #         _.append(y[np.ix_(x_flips, x_flips)])
        #     x = np.array(_)
        #
        #     import json
        #     raise Exception(
        #         json.dumps((x[2]*self.h2w).tolist())
        #         )
        pre_wfns_script = None

        # nt_spec = np.array([
        #         [0, 0, 0, 0, 0, 0, 1, 0, 0],
        #         [0, 0, 0, 2, 0, 0, 0, 0, 0]
        #     ])
        # degeneracies = [
        #     [
        #         s,
        #         (s - nt_spec[0] + nt_spec[1])
        #     ]for s in states if np.dot(s, nt_spec[0]) > 0
        # ]
        # degeneracies = [np.array(x).tolist() for x in degeneracies]
        # states = np.array(states).tolist()
        # for pair in degeneracies:
        #     for p in pair:
        #         if p not in states:
        #             states.append(p)

        degeneracies = None
        print_report = False
        nielsen_tolerance = 10 if degeneracies is None else 500
        gaussian_tolerance = 10 if degeneracies is not None else 50
        self.run_PT_test(
            tag,
            file_name,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            log=False,
            verbose=False,
            print_x=True,
            calculate_intensities=True,
            degeneracies=degeneracies,
            print_report=print_report,
            nielsen_tolerance=nielsen_tolerance,
            gaussian_tolerance=gaussian_tolerance,
            pre_wfns_script=pre_wfns_script
            # zero_element_warning=False

        )

    @validationTest
    def test_HOONOVPTCartesiansEmbdedded(self):

        tag = 'HOONO Cartesians'
        file_name = "HOONO_freq.fchk"

        mol = Molecule.from_file(TestManager.test_data(file_name)).get_embedded_molecule()  # embed_properties=False)

        # mol.zmatrix = [
        #     [0, -1, -1, -1],
        #     [1,  0, -1, -1],
        #     [2,  0, 1, -1]
        # ]
        internals = None

        n_atoms = 5
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        gaussian_energies = self.gaussian_data['HOONO']['zpe']
        gaussian_freqs = self.gaussian_data['HOONO']['freqs']

        degeneracies = None
        print_report = False
        nielsen_tolerance = 10 if degeneracies is None else 500
        gaussian_tolerance = 10 if degeneracies is not None else 50
        self.run_PT_test(
            tag,
            mol,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            log=False,
            verbose=False,
            calculate_intensities=True,
            degeneracies=degeneracies,
            print_report=print_report,
            nielsen_tolerance=nielsen_tolerance,
            gaussian_tolerance=gaussian_tolerance,
            print_x=True
            # zero_element_warning=False

        )

    @validationTest
    def test_HOONOVPTCartesiansDegenerate(self):

        tag = 'HOONO Cartesians'
        file_name = "HOONO_freq.fchk"

        internals = None

        n_atoms = 5
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        gaussian_energies = self.gaussian_data['HOONO']['zpe']
        gaussian_freqs = self.gaussian_data['HOONO']['freqs']

        # """
        # 0.337366D+01  0.708754D+00  0.685684D+00  0.704381D+00  0.386143D-01 0.189447D+00  0.238868D-01  0.000000D+00  0.000000D+00
        # """

        # def pre_wfns_script(hammer, states):
        #     harm, corrs, x = hammer.get_Nielsen_energies(states, return_split=True)
        #     x_flips = np.argsort([8, 6, 7, 5, 4, 3, 2, 1, 0])
        #     _ = []
        #     for y in x:
        #         _.append(y[np.ix_(x_flips, x_flips)])
        #     x = np.array(_)
        #
        #     import json
        #     raise Exception(
        #         json.dumps((x[2]*self.h2w).tolist())
        #         )
        pre_wfns_script = None

        nt_spec = np.array([
            [0, 0, 0, 0, 0, 0, 1, 0, 0],
            [0, 0, 0, 2, 0, 0, 0, 0, 0]
        ])
        degeneracies = [
            [
                s,
                (s - nt_spec[0] + nt_spec[1])
            ] for s in states if np.dot(s, nt_spec[0]) > 0
        ]
        degeneracies = [np.array(x).tolist() for x in degeneracies]
        states = np.array(states).tolist()
        for pair in degeneracies:
            for p in pair:
                if p not in states:
                    states.append(p)

        print_report = False
        nielsen_tolerance = 10 if degeneracies is None else 500
        gaussian_tolerance = 10 if degeneracies is not None else 50
        self.run_PT_test(
            tag,
            file_name,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            log=True,
            verbose=True,
            calculate_intensities=True,
            degeneracies=degeneracies,
            print_report=print_report,
            nielsen_tolerance=nielsen_tolerance,
            gaussian_tolerance=gaussian_tolerance,
            pre_wfns_script=pre_wfns_script,
            gaussian_resonance_handling=True
            # zero_element_warning=False

        )

    #endregion

    #region Methane

    gaussian_data['CH2DT'] = {
        'zpe': [8311.227, 8205.865],
        'freqs': np.array([
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
             [5157.343, 4946.234],
             [4651.438, 4449.136],
             [4517.395, 4325.685],
             [4436.517, 4250.290],
             [4255.193, 4077.953],
             [4165.552, 3999.186],
             [5470.079, 5252.170],
             [5092.977, 4895.552],
             [4587.071, 4407.847],
             [4453.028, 4280.602],
             [4372.150, 4204.590],
             [4190.826, 4027.168],
             [4101.185, 3948.975],
             [4279.364, 4131.666],
             [3773.459, 3657.329],
             [3639.416, 3518.558],
             [3558.538, 3448.813],
             [3377.214, 3269.456],
             [3287.573, 3183.220],
             [3396.356, 3302.019],
             [3262.313, 3175.698],
             [3181.436, 3094.860],
             [3000.111, 2904.013],
             [2910.471, 2822.890],
             [2756.408, 2698.844],
             [2675.530, 2618.649],
             [2494.205, 2435.125],
             [2404.565, 2356.088],
             [2541.487, 2483.423],
             [2360.163, 2319.432],
             [2270.522, 2230.076],
             [2279.285, 2236.240],
             [2189.644, 2143.874],
             [2008.320, 1973.845]
        ])
    }
    @validationTest
    def test_CH2DTVPTCartesians(self):

        tag = 'CH2DT Cartesians'
        file_name = "CH2DT_freq.fchk"

        internals = None

        n_atoms = 5
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        gaussian_energies = self.gaussian_data['CH2DT']['zpe']
        gaussian_freqs = self.gaussian_data['CH2DT']['freqs']

        print_report = True
        nielsen_tolerance = 1
        gaussian_tolerance = 50
        self.run_PT_test(
            tag,
            file_name,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            print_report=print_report,
            nielsen_tolerance=nielsen_tolerance,
            gaussian_tolerance=gaussian_tolerance
            , log=True
        )

    #endregion Methane

    #region Water Clusters

    gaussian_data['WaterDimer'] = {
        'zpe': [10133.860, 9909.756],
        'freqs': np.array([
         [ 3935.49,  3742.92],
         [ 3914.94,  3752.15],
         [ 3814.08,  3652.41],
         [ 3718.19,  3584.14],
         [ 1650.02,  1592.65],
         [ 1629.21,  1585.96],
         [  631.34,  505.605],
         [ 362.703,  295.834],
         [ 183.777,  141.372],
         [ 154.306,  110.995],
         [ 146.544,  150.517],
         [ 127.117,   69.163],
         [ 7870.98,  7393.56],
         [ 7829.88,  7368.49],
         [ 7628.16,  7224.88],
         [ 7436.38,  7016.02],
         [ 3300.05,  3152.47],
         [ 3258.42,  3144.16],
         [ 1262.68,  921.053],
         [ 725.405,  488.907],
         [ 367.554,  268.882],
         [ 308.612,  207.465],
         [ 293.089,  299.766],
         [ 254.234,  114.677],
         [ 7850.43,  7494.65],
         [ 7749.57,  7239.31],
         [ 7653.68,  7322.97],
         [ 5585.51,  5334.87],
         [  5564.7,   5314.4],
         [ 4566.83,   4249.7],
         [ 4298.19,  4016.06],
         [ 4119.27,  3875.58],
         [  4089.8,  3839.37],
         [ 4082.03,  3892.36],
         [ 4062.61,   3810.2],
         [ 7729.02,  7402.98],
         [ 7633.13,  7271.26],
         [ 5564.96,  5328.22],
         [ 5544.15,  5337.03],
         [ 4546.28,  4261.39],
         [ 4277.64,  4024.52],
         [ 4098.72,  3895.81],
         [ 4069.25,  3864.62],
         [ 4061.48,  3903.06],
         [ 4042.06,  3844.89],
         [ 7532.27,  7230.66],
         [  5464.1,  5244.06],
         [ 5443.29,  5222.11],
         [ 4445.42,  4159.77],
         [ 4176.78,  3928.52],
         [ 3997.86,  3794.28],
         [ 3968.39,  3763.45],
         [ 3960.62,  3844.23],
         [  3941.2,  3692.21],
         [ 5368.22,   5164.6],
         [  5347.4,  5168.41],
         [ 4349.53,  4139.08],
         [ 4080.89,  3889.46],
         [ 3901.97,   3739.5],
         [  3872.5,  3704.84],
         [ 3864.74,  3750.79],
         [ 3845.31,   3657.1],
         [ 3279.23,  3172.37],
         [ 2281.36,  2107.39],
         [ 2012.73,  1852.95],
         [  1833.8,  1732.29],
         [ 1804.33,  1699.13],
         [ 1796.57,    1746.],
         [ 1777.14,  1656.44],
         [ 2260.55,  2094.51],
         [ 1991.91,  1862.32],
         [ 1812.99,  1729.35],
         [ 1783.52,  1700.18],
         [ 1775.76,  1736.87],
         [ 1756.33,  1651.98],
         [ 994.042,  745.791],
         [ 815.116,   620.62],
         [ 785.646,  595.362],
         [ 777.884,  646.479],
         [ 758.457,   525.68],
         [ 546.479,  389.065],
         [ 517.009,  374.506],
         [ 509.247,  405.832],
         [  489.82,  336.402],
         [ 338.083,  235.655],
         [ 330.321,  287.882],
         [ 310.894,  207.357],
         [  300.85,  261.693],
         [ 281.423,   169.59],
         [ 273.662,  201.333]
    ])
    }
    #Paper
    @validationTest
    def test_WaterDimerVPTCartesians(self):
        # the high-frequency stuff agrees with Gaussian, but not the low-freq

        tag = 'Water Dimer Cartesians'
        file_name = "water_dimer_freq.fchk"

        internals = None

        n_atoms = 6
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)#[:6]

        gaussian_energies = self.gaussian_data['WaterDimer']['zpe']
        gaussian_freqs = self.gaussian_data['WaterDimer']['freqs']

        chk = os.path.expanduser('~/Desktop/dimer_chk2.hdf5')
        print_report = False
        nielsen_tolerance = 50
        gaussian_tolerance = 50
        self.run_PT_test(
            tag,
            file_name,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            log=True,
            verbose=True,
            print_profile=True,
            # profile_filter='Combinatorics/Permutations',
            print_report=print_report,
            nielsen_tolerance=nielsen_tolerance,
            gaussian_tolerance=gaussian_tolerance,
            calculate_intensities=True
            , checkpoint=chk
            , use_cached_representations=False
            , state_space_filters={
                # (2, 0): BasisStateSpace(
                #     HarmonicOscillatorProductBasis(n_modes),
                #     states[:1],
                # ).apply_selection_rules([[]]), # get energies right
                (1, 1): BasisStateSpace(
                    HarmonicOscillatorProductBasis(n_modes),
                    states[:1]
                ).apply_selection_rules([[-1], [], [1]])

            }
            # , parallelized=True
        )

    @validationTest
    def test_WaterDimerVPTCartesiansParallel(self):
        # the high-frequency stuff agrees with Gaussian, but not the low-freq

        tag = 'Water Dimer Cartesians'
        file_name = "water_dimer_freq.fchk"

        internals = None

        n_atoms = 6
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        gaussian_energies = self.gaussian_data['WaterDimer']['zpe']
        gaussian_freqs = self.gaussian_data['WaterDimer']['freqs']

        print_report = False
        nielsen_tolerance = 50
        gaussian_tolerance = 50
        self.run_PT_test(
            tag,
            file_name,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            log=True,
            verbose=True,
            print_report=print_report,
            nielsen_tolerance=nielsen_tolerance,
            gaussian_tolerance=gaussian_tolerance
            , parallelized=True
        )

    @validationTest
    def test_WaterDimerVPTCartesiansPararalelMemConstrained(self):
        # the high-frequency stuff agrees with Gaussian, but not the low-freq

        tag = 'Water Dimer Cartesians'
        file_name = "water_dimer_freq.fchk"

        internals = None

        n_atoms = 6
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        gaussian_energies = self.gaussian_data['WaterDimer']['zpe']
        gaussian_freqs = self.gaussian_data['WaterDimer']['freqs']

        print_report = False
        nielsen_tolerance = 1
        gaussian_tolerance = 50
        self.run_PT_test(
            tag,
            file_name,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            print_report=print_report,
            nielsen_tolerance=nielsen_tolerance,
            gaussian_tolerance=gaussian_tolerance
            , parallelized=True
            , chunk_size=int(2e4)
        )

    @validationTest
    def test_WaterDimerVPTCartesiansSubspace(self):
        # the high-frequency stuff agrees with Gaussian, but not the low-freq

        internals = None

        n_modes = 6 * 3 - 6
        mode_selection = [-3, -2, -1]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)

        states = self.get_states(2, n_modes)

        # coupled_states = self.get_states(5, n_modes, max_quanta=5)

        with BlockProfiler("WaterDimer", print_res=False):
            wfns = self.get_VPT2_wfns(
                "water_dimer_freq.fchk",
                internals,
                states,
                regenerate=True,
                mode_selection=mode_selection,
                log=True
                # , parallelized=True
            )

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")
        engs = h2w * wfns.energies
        freqs = engs[1:] - engs[0]
        harm_engs = h2w * wfns.zero_order_energies
        harm_freq = harm_engs[1:] - harm_engs[0]

        gaussian_engs = [10133.860, 9909.756]
        gaussian_freqs = np.array([
            [3935.490, 3742.918],
            [3914.939, 3752.151],
            [3814.079, 3652.414],
            [3718.192, 3584.139],
            [1650.023, 1592.653],
            [1629.210, 1585.962],
            [631.340, 505.605],
            [362.703, 295.834],
            [183.777, 141.372],
            [154.306, 110.995],
            [146.544, 150.517],
            [127.117, 69.163],

            [7870.980, 7393.560],
            [7829.879, 7368.493],
            [7628.159, 7224.882],
            [7436.384, 7016.025],
            [3300.045, 3152.473],
            [3258.421, 3144.157],
            [1262.679, 921.053],
            [725.405, 488.907],
            [367.554, 268.882],
            [308.612, 207.465],
            [293.089, 299.766],
            [254.234, 114.677],

            [7850.429, 7494.650],
            [7749.569, 7239.308],
            [7729.019, 7402.976],
            [7653.682, 7322.974],
            [7633.131, 7271.264],
            [7532.271, 7230.663],
            [5585.513, 5334.869],
            [5564.962, 5328.224],
            [5464.102, 5244.056],
            [5368.215, 5164.597],
            [5564.700, 5314.396],
            [5544.150, 5337.031],
            [5443.290, 5222.111],
            [5347.402, 5168.407],
            [3279.233, 3172.374],
            [4566.830, 4249.695],
            [4546.279, 4261.388],
            [4445.419, 4159.774],
            [4349.531, 4139.077],
            [2281.362, 2107.393],
            [2260.550, 2094.511],
            [4298.193, 4016.063],
            [4277.642, 4024.523],
            [4176.782, 3928.515],
            [4080.894, 3889.457],
            [2012.725, 1852.952],
            [1991.913, 1862.320],
            [994.042, 745.791],
            [4119.267, 3875.578],
            [4098.716, 3895.805],
            [3997.856, 3794.279],
            [3901.969, 3739.502],
            [1833.800, 1732.294],
            [1812.987, 1729.354],
            [815.116, 620.620],
            [546.479, 389.065],
            [4089.796, 3839.370],
            [4069.245, 3864.621],
            [3968.385, 3763.445],
            [3872.498, 3704.835],
            [1804.329, 1699.128],
            [1783.516, 1700.178],
            [785.646, 595.362],
            [517.009, 374.506],
            [338.083, 235.655],
            [4082.035, 3892.356],
            [4061.484, 3903.055],
            [3960.624, 3844.234],
            [3864.736, 3750.792],
            [1796.567, 1745.999],
            [1775.755, 1736.874],
            [777.884, 646.479],
            [509.247, 405.832],
            [330.321, 287.882],
            [300.850, 261.693],
            [4062.607, 3810.202],
            [4042.056, 3844.889],
            [3941.197, 3692.207],
            [3845.309, 3657.100],
            [1777.140, 1656.439],
            [1756.328, 1651.982],
            [758.457, 525.680],
            [489.820, 336.402],
            [310.894, 207.357],
            [281.423, 169.590],
            [273.662, 201.333]
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

        # self.assertLess(
        #     np.max(np.abs(freqs[:ns] - gaussian_freqs[:ns, 1])),
        #     1)

    @validationTest
    def test_WaterDimerVPTCartesiansHarmonic(self):
        # the high-frequency stuff agrees with Gaussian, but not the low-freq
        internals = None

        n_modes = 6 * 3 - 6
        mode_selection = None#[-3, -2, -1]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)

        states = self.get_states(2, n_modes)

        # coupled_states = self.get_states(5, n_modes, max_quanta=5)

        with BlockProfiler("WaterDimer", print_res=False):
            wfns = self.get_VPT2_wfns(
                "water_dimer_freq.fchk",
                internals,
                states,
                regenerate=True,
                mode_selection=mode_selection,
                log=True,
                order=0
                # , parallelized=True
            )

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")
        engs = h2w * wfns.energies
        freqs = engs[1:] - engs[0]
        harm_engs = h2w * wfns.zero_order_energies
        harm_freq = harm_engs[1:] - harm_engs[0]

        gaussian_engs = [10133.860, 9909.756]
        gaussian_freqs = np.array([
            [3935.490, 3742.918],
            [3914.939, 3752.151],
            [3814.079, 3652.414],
            [3718.192, 3584.139],
            [1650.023, 1592.653],
            [1629.210, 1585.962],
            [631.340, 505.605],
            [362.703, 295.834],
            [183.777, 141.372],
            [154.306, 110.995],
            [146.544, 150.517],
            [127.117, 69.163],

            [7870.980, 7393.560],
            [7829.879, 7368.493],
            [7628.159, 7224.882],
            [7436.384, 7016.025],
            [3300.045, 3152.473],
            [3258.421, 3144.157],
            [1262.679, 921.053],
            [725.405, 488.907],
            [367.554, 268.882],
            [308.612, 207.465],
            [293.089, 299.766],
            [254.234, 114.677],

            [7850.429, 7494.650],
            [7749.569, 7239.308],
            [7729.019, 7402.976],
            [7653.682, 7322.974],
            [7633.131, 7271.264],
            [7532.271, 7230.663],
            [5585.513, 5334.869],
            [5564.962, 5328.224],
            [5464.102, 5244.056],
            [5368.215, 5164.597],
            [5564.700, 5314.396],
            [5544.150, 5337.031],
            [5443.290, 5222.111],
            [5347.402, 5168.407],
            [3279.233, 3172.374],
            [4566.830, 4249.695],
            [4546.279, 4261.388],
            [4445.419, 4159.774],
            [4349.531, 4139.077],
            [2281.362, 2107.393],
            [2260.550, 2094.511],
            [4298.193, 4016.063],
            [4277.642, 4024.523],
            [4176.782, 3928.515],
            [4080.894, 3889.457],
            [2012.725, 1852.952],
            [1991.913, 1862.320],
            [994.042, 745.791],
            [4119.267, 3875.578],
            [4098.716, 3895.805],
            [3997.856, 3794.279],
            [3901.969, 3739.502],
            [1833.800, 1732.294],
            [1812.987, 1729.354],
            [815.116, 620.620],
            [546.479, 389.065],
            [4089.796, 3839.370],
            [4069.245, 3864.621],
            [3968.385, 3763.445],
            [3872.498, 3704.835],
            [1804.329, 1699.128],
            [1783.516, 1700.178],
            [785.646, 595.362],
            [517.009, 374.506],
            [338.083, 235.655],
            [4082.035, 3892.356],
            [4061.484, 3903.055],
            [3960.624, 3844.234],
            [3864.736, 3750.792],
            [1796.567, 1745.999],
            [1775.755, 1736.874],
            [777.884, 646.479],
            [509.247, 405.832],
            [330.321, 287.882],
            [300.850, 261.693],
            [4062.607, 3810.202],
            [4042.056, 3844.889],
            [3941.197, 3692.207],
            [3845.309, 3657.100],
            [1777.140, 1656.439],
            [1756.328, 1651.982],
            [758.457, 525.680],
            [489.820, 336.402],
            [310.894, 207.357],
            [281.423, 169.590],
            [273.662, 201.333]
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

        # self.assertLess(
        #     np.max(np.abs(freqs[:ns] - gaussian_freqs[:ns, 1])),
        #     1)

    @debugTest
    def test_WaterDimerVPTInternals(self):

        """
      1          1           0        0.000000    0.000000    0.000000
      2          8           0        0.000000    0.000000    0.960485
      3          1           0        0.938795    0.000000    1.199654
      4          8           0        2.745004    0.000000    1.930528
      5          1           0        2.858217    0.761260    2.508248
      6          1           0        2.858217   -0.761260    2.508248
      """

        tag = 'Water Dimer Internals'
        file_name = "water_dimer_freq.fchk"

        COM = -3
        A = -2
        C = -1
        X = 1000
        LHF = 0
        LO = 1
        SH = 2
        RO = 3
        RH1 = 4
        RH2 = 5
        # internals = [
        #     [LHF,  X,   X,   X],
        #     [LO, LHF,   X,   X],
        #     [SH,  LO, LHF,   X],
        #     [RH2, SH,  LO, LHF], # get out of plane
        #     [RO,  LO, RH2, LHF],
        #     [RH1, RO, RH2, LHF]
        # ]

        # internals = [
        #     [LHF,   X,   X,   X],
        #     [LO,  LHF,   X,   X],
        #     [SH,   LO, LHF,   X],
        #     [RO,   LO,  SH, LHF],
        #     [RH1,  RO,  LO, LHF],
        #     [RH2,  RO, RH1,  LO]
        # ]

        internals = [
            [LHF,   X,   X,    X],
            [LO,  LHF,   X,    X],
            [SH,   LO,  LHF,   X],
            [RO,   LO,  LHF,   C],
            [RH1,  RO,   SH, LHF],
            [RH2,  RO,  RH1, LHF]
        ]

        n_atoms = 6
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)  # [:6]

        gaussian_energies = self.gaussian_data['WaterDimer']['zpe']
        gaussian_freqs = self.gaussian_data['WaterDimer']['freqs']

        # chk = os.path.expanduser('~/Desktop/dimer_chk2.hdf5')
        print_report = False
        nielsen_tolerance = 50
        gaussian_tolerance = 50
        self.run_PT_test(
            tag,
            file_name,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            log=False,
            verbose=False,
            print_profile=False,
            # profile_filter='Combinatorics/Permutations',
            print_report=print_report,
            nielsen_tolerance=nielsen_tolerance,
            gaussian_tolerance=gaussian_tolerance,
            calculate_intensities=True
            # , checkpoint=chk
            , use_cached_representations=False
            , state_space_filters={
                # (2, 0): BasisStateSpace(
                #     HarmonicOscillatorProductBasis(n_modes),
                #     states[:1],
                # ).apply_selection_rules([[]]), # get energies right
                (1, 1): BasisStateSpace(
                    HarmonicOscillatorProductBasis(n_modes),
                    states[:1]
                ).apply_selection_rules([[-1], [], [1]])

            }
            # , parallelized=True
        )

    @validationTest
    def test_WaterTrimerVPTCartesians(self):
        tag = 'Water Trimer Cartesians'
        file_name = "water_trimer_freq.fchk"

        internals = None

        n_atoms = 9
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes, target_modes=[-1, -2, -3, -4, -5, -6]) # [:6]
        # raise Exception(states)

        gaussian_energies = None#self.gaussian_data['WaterDimer']['zpe']
        gaussian_freqs = None#self.gaussian_data['WaterDimer']['freqs']

        print_report = True
        nielsen_tolerance = 50
        gaussian_tolerance = 50
        self.run_PT_test(
            tag,
            file_name,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            log=True,
            verbose=True,
            print_profile=True,
            # profiling_mode='deterministic',
            # profile_filter='Combinatorics/Permutations',
            print_report=print_report,
            nielsen_tolerance=nielsen_tolerance,
            gaussian_tolerance=gaussian_tolerance
            # , checkpoint=chk
            # , parallelized=True
            , initialization_timeout=2
            , chunk_size=int(5e6)
            # , memory_constrained=True
            , state_space_terms=((1, 0), (2, 0))
            , calculate_intensities=True
            # , direct_sum_chunk_size=int(1e3)
        )

    #endregion Water Clusters

    #endregion Test Systems

    #region Test Action Expansions

    @validationTest
    def test_HOHCartesianActionExpansion(self):

        internals = None

        ham = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("HOH_freq.fchk"),
            internals=internals,
            log=False
        )

        expansion, wfns = ham.get_action_expansion()

        g0, w, x = [self.h2w * y for y in expansion]

        gaussian_x = np.array([
            [-0.490388e+02,  0.0000000000,  0.0000000000],
            [-0.166364e+03, -0.440530e+02,  0.0000000000],
            [-0.230618e+02, -0.200356e+02, -0.140238e+02]
        ])[np.ix_((2, 1, 0), (2, 1, 0))]

        self.assertTrue(
            np.allclose(gaussian_x, x, atol=.5),
            msg="wat \n {} \n {}".format(
                gaussian_x, x
            )
        )

        # raise Exception([
        #     self.h2w * x,
        #     self.h2w * w,
        #     self.h2w * g0
        # ])

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
    def test_HODCartesianPotential(self):

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
    def test_HODCartesianCoriolis(self):

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
    def test_OCHDCartesiansPotential(self):

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
    def test_OCHDCartesiansCoriolis(self):

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
        # coupled_states = self.get_states(5, 3, max_quanta=5)

        wfns = self.get_VPT2_wfns(
            "HOD_freq.fchk",
            internals,
            states,
            regenerate=True,
            # coupled_states=coupled_states,
            log=True
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
    def test_HODIntensitiesHarmonic(self):

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
        # coupled_states = self.get_states(5, 3, max_quanta=5)

        wfns = self.get_VPT2_wfns(
            "HOD_freq.fchk",
            internals,
            states,
            regenerate=True,
            # coupled_states=coupled_states,
            log=True,
            order=0
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

        # self.assertEquals(
        #     [round(x, 2) for x in gaussian_harm_ints[1:]],
        #     list(np.round(harm_ints[1:4], 2))
        # )
        # self.assertEquals(
        #     [round(x, 2) for x in gaussian_ints[1:]],
        #     list(np.round(ints[1:10], 2))
        # )

        wfns.dipole_paritioning="intuitive"
        ints = wfns.intensities

        harm_ints = wfns.zero_order_intensities

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

    @validationTest
    def test_HODIntensitiesCartesian(self):

        internals = None
        states = (
            (0, 0, 0),
            (0, 0, 1), (0, 1, 0), (1, 0, 0),
            (0, 0, 2), (0, 2, 0), (2, 0, 0),
            (0, 1, 1), (1, 0, 1), (1, 1, 0)
        )
        coupled_states = None #self.get_states(5, 3, max_quanta=5)
        wfns = self.get_VPT2_wfns(
            "HOD_freq.fchk",
            internals,
            states,
            regenerate=True,
            coupled_states=coupled_states,
            log=False
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
        self.assertTrue(
            np.allclose(
                gaussian_ints[1:],
                ints[1:10],
                atol=3
            )
        )

    @validationTest
    def test_HOHIntensitiesCartesian4thOrder(self):

        internals = None
        states = (
            (0, 0, 0),
            (0, 0, 1), (0, 1, 0), (1, 0, 0),
            (0, 0, 2), (0, 2, 0), (2, 0, 0),
            (0, 1, 1), (1, 0, 1), (1, 1, 0)
        )
        # coupled_states = self.get_states(5, 3, max_quanta=5)

        wfns = self.get_VPT2_wfns(
            "HOH_freq.fchk",
            internals,
            states,
            regenerate=True
            # coupled_states=coupled_states,
            , log=True
            , verbose=True
            , order=4
            , expansion_order=2
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

        # breakdown = wfns.generate_intensity_breakdown(include_wavefunctions=False)
        #
        # import json
        # raise Exception(
        #     json.dumps(
        #         [
        #             (np.array(breakdown['frequencies']) * self.h2w).tolist(),
        #             breakdown['breakdowns']['Full']['intensities'].tolist(),
        #             breakdown['breakdowns']['Linear']['intensities'].tolist(),
        #             breakdown['breakdowns']['Quadratic']['intensities'].tolist(),
        #             breakdown['breakdowns']['Cubic']['intensities'].tolist(),
        #             breakdown['breakdowns']['Constant']['intensities'].tolist()
        #         ]
        #     )
        # )

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
            3937.525,
            3803.300,
            1622.303
        ]
        gaussian_harm_ints = [
            0.,
            67.02031831,
            4.14281519,
            67.45606550
        ]
        gaussian_freqs = [
            0.,
            3744.734,
            3621.994,
            1572.707,
            7391.391,
            7155.881,
            3117.366,
            7200.364,
            5294.379,
            5174.665
        ]
        gaussian_ints = [
            0.,
            64.17034278,
            3.11359907,
            68.32307006,
            0.01482583,
            0.31497697,
            0.55501343,
            2.20970792,
            3.76327875,
            0.06230526
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

        # self.assertEquals(
        #     [round(x, 2) for x in gaussian_harm_ints[1:]],
        #     list(np.round(harm_ints[1:4], 2))
        # )
        # self.assertEquals(
        #     [round(x, 2) for x in gaussian_ints[1:]],
        #     list(np.round(ints[1:10], 2))
        # )

    @validationTest
    def test_HOHIntensities(self):

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
        # coupled_states = self.get_states(5, 3, max_quanta=5)

        # internals=None
        with BlockProfiler("HOH Intenstities", print_res=False):
            wfns = self.get_VPT2_wfns(
                "HOH_freq.fchk",
                internals,
                states,
                regenerate=True
                # coupled_states=coupled_states,
                , log=True
                , order=2
            )

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        # trying to turn off the orthogonality condition
        # states = wfns.corrs.states
        # for i,s in enumerate(states):
        #     wfns.corrs.wfn_corrections[i, 2, s] = 0 # turn off second correction

        wfns.dipole_partitioning = 'standard'#'intuitive'
        engs = h2w * wfns.energies
        freqs = engs - engs[0]
        ints = wfns.intensities

        harm_engs = h2w * wfns.zero_order_energies
        harm_freqs = harm_engs - harm_engs[0]
        harm_ints = wfns.zero_order_intensities

        # breakdown = wfns.generate_intensity_breakdown(include_wavefunctions=False)

        # import json
        # raise Exception(
        #     json.dumps(
        #         [
        #             (np.array(breakdown['frequencies']) * self.h2w).tolist(),
        #             breakdown['breakdowns']['Full']['intensities'].tolist(),
        #             breakdown['breakdowns']['Linear']['intensities'].tolist(),
        #             breakdown['breakdowns']['Quadratic']['intensities'].tolist(),
        #             breakdown['breakdowns']['Cubic']['intensities'].tolist(),
        #             breakdown['breakdowns']['Constant']['intensities'].tolist()
        #             ]
        #     )
        # )

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
            3937.525,
            3803.300,
            1622.303
        ]
        gaussian_harm_ints = [
            0.,
            67.02031831,
             4.14281519,
            67.45606550
        ]
        gaussian_freqs = [
            0.,
            3744.734,
            3621.994,
            1572.707,
            7391.391,
            7155.881,
            3117.366,
            7200.364,
            5294.379,
            5174.665
        ]
        gaussian_ints = [
            0.,
            64.17034278,
            3.11359907,
            68.32307006,
            0.01482583,
            0.31497697,
            0.55501343,
            2.20970792,
            3.76327875,
            0.06230526
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
    def test_HOHIntensitiesInternalCartesianEquality(self):
        # Add test to check that internal/cartesian intensities
        # should be the same

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
        # coupled_states = [
        #             BasisStateSpace.from_quanta(
        #                 HarmonicOscillatorProductBasis(3),
        #                 range(12)
        #             )
        #                  ]*2

        with BlockProfiler("HOH Intenstities", print_res=False):
            wfns_internal, ham = self.get_VPT2_wfns_and_ham(
                "HOH_freq.fchk",
                internals,
                states,
                regenerate=True
                # , coupled_states=coupled_states
                , log=False
                , order=2
                # , v3 = 0
                # , v4 = 0
                # , t3 = 0
                # , t4 = 0
            )

        internals = None
        with BlockProfiler("HOH Intenstities", print_res=False):
            wfns_cartesian, ham = self.get_VPT2_wfns_and_ham(
                "HOH_freq.fchk",
                internals,
                states,
                regenerate=True
                # , coupled_states=coupled_states
                , log=False
                , order=2
                # , v3 = 0
                # , v4 = 0
                # , t3 = 0
                # , t4 = 0
            )


        int_internal = wfns_internal.intensities
        int_cartesin = wfns_cartesian.intensities

        self.assertTrue(np.allclose(int_internal, int_cartesin, atol=.3),
            msg='internals and cartesians differ; difference: {}'.format(int_internal-int_cartesin)
        )


    @inactiveTest
    def test_HOHIntensitiesSandbox(self):

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
        # coupled_states = [
        #             BasisStateSpace.from_quanta(
        #                 HarmonicOscillatorProductBasis(3),
        #                 range(12)
        #             )
        #                  ]*2

        internals=None
        with BlockProfiler("HOH Intenstities", print_res=False):
            wfns, ham = self.get_VPT2_wfns_and_ham(
                "HOH_freq.fchk",
                internals,
                states,
                regenerate=True
                # , coupled_states=coupled_states
                , log=False
                , order=6
                # , v3 = 0
                # , v4 = 0
                # , t3 = 0
                # , t4 = 0
            )

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        # trying to turn off the orthogonality condition
        # states = wfns.corrs.states
        # for i,s in enumerate(states):
        #     wfns.corrs.wfn_corrections[i, 2, s] = 0 # turn off second correction

        # wfns.dipole_partitioning = 'intuitive'
        engs = h2w * wfns.energies
        freqs = engs - engs[0]
        ints = wfns.intensities

        harm_engs = h2w * wfns.zero_order_energies
        harm_freqs = harm_engs - harm_engs[0]
        harm_ints = wfns.zero_order_intensities

        import json
        breakdown = wfns.generate_intensity_breakdown(include_wavefunctions=False)
        raise Exception(
            json.dumps(
                [
                    (np.array(breakdown['frequencies']) * self.h2w).tolist(),
                    # ints.tolist(),
                    breakdown['breakdowns']['Full']['intensities'].tolist(),
                    breakdown['breakdowns']['Linear']['intensities'].tolist(),
                    breakdown['breakdowns']['Quadratic']['intensities'].tolist(),
                    breakdown['breakdowns']['Cubic']['intensities'].tolist(),
                    breakdown['breakdowns']['Constant']['intensities'].tolist(),
                    breakdown['breakdowns']['Order0']['intensities'].tolist(),
                    breakdown['breakdowns']['Order1']['intensities'].tolist()
                    ]
            )
        )

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

        print_specs = True
        if print_specs:
            report = "Harmonic:   State     Freq.   Int.  \n" + "\n".join(
                " " * 12 + "{} {:>7.2f} {:>7.4f}       ".format(s, f, i)
                for s, f, i in zip(states, harm_freqs, harm_ints)
            )
            print(report)
            report = "Anharmonic: State     Freq.   Int.\n" + "\n".join(
                " " * 12 + "{} {:>7.2f} {:>7.4f}  ".format(s, f, i)
                for s, f, i in zip(states, freqs, ints)
            )
            print(report)

    @validationTest
    def test_OCHHIntensitiesSanbox(self):

        tag = "OCHH Intenstities"
        file_name = "OCHH_freq.fchk"

        internals = [
            [0, -1, -1, -1],
            [1,  0, -1, -1],
            [2,  1,  0, -1],
            [3,  1,  0,  2]
        ]

        # internals = None

        n_atoms = 4
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = self.get_states(3, n_modes)

        degeneracies = (
            [
                [0, 0, 0, 0, 0, 1],
                [0, 1, 0, 1, 0, 0],
            ],
            [
                [0, 0, 0, 0, 0, 2],
                [0, 1, 0, 1, 0, 1],
            ],
            [
                [0, 0, 0, 0, 1, 1],
                [0, 1, 0, 1, 1, 0],
            ],
            [
                [0, 0, 0, 1, 0, 1],
                [0, 1, 0, 2, 0, 0],
            ],
            [
                [0, 0, 1, 0, 0, 1],
                [0, 1, 1, 1, 0, 0],
            ],
            [
                [0, 1, 0, 0, 0, 1],
                [0, 2, 0, 1, 0, 0],
            ],
            [
                [1, 0, 0, 0, 0, 1],
                [1, 1, 0, 1, 0, 0],
            ]

            #
            # [0, 0, 0, 0, 1, 1],
            # [0, 1, 0, 1, 1, 0],
            #
            # [0, 0, 0, 1, 0, 1],
            # [0, 1, 0, 2, 0, 0],
            #
            # [0, 0, 1, 0, 0, 1],
            # [0, 1, 1, 1, 0, 0],
            #
            # [0, 1, 0, 0, 0, 1],
            # [0, 2, 0, 1, 0, 0],
            #
            # [1, 0, 0, 0, 0, 1],
            # [1, 1, 0, 1, 0, 0]
            # [
            #     [0, 0, 0, 0, 0, 2],
            #     [0, 1, 0, 0, 0, 1],
            #     [0, 1, 0, 1, 0, 1],
            # ]
            # [
            #     [0, 0, 0, 0, 1, 1],
            #     [0, 1, 0, 1, 1, 0]
            # ],
            # [
            #     [1, 0, 0, 0, 0, 1],
            #     [1, 1, 0, 1, 0, 0]
            # ],
            # [
            #     [0, 0, 0, 1, 0, 1],
            #     [0, 1, 0, 2, 0, 0]
            # ],
            # [
            #     [0, 0, 0, 0, 0, 2],
            #     [0, 1, 0, 1, 0, 1],
            #     [0, 2, 0, 2, 0, 0]
            # ]
        )

        degeneracies = [np.array(x).tolist() for x in degeneracies]
        states = np.array(states).tolist()
        for pair in degeneracies:
            for p in pair:
                if p not in states:
                    states.append(p)

        basis = HarmonicOscillatorProductBasis(n_modes)

        with BlockProfiler(tag, inactive=True):
            wfns = self.get_VPT2_wfns(
                file_name,
                internals,
                states,
                regenerate=True
                # coupled_states=coupled_states,
                , log=True
                , order=2
                , degeneracies=degeneracies
                # , v3 = 0
                # , t3 = 0
                # , v4 = 0
                # , t4 = 0
                # , selection_rules=[
                #     [r for r in basis.selection_rules("x", "x", "x") if len(r) <= 2],
                #     [r for r in basis.selection_rules("x", "x", "x", "x") if len(r) <= 3]
                # ]
                # , state_space_iterations=5
            )

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        # with BlockProfiler(tag, print_res=True):
        #     raise Exception(wfns.intensities)

        # import json
        # with BlockProfiler(tag, print_res=True):
        #     breakdown = wfns.generate_intensity_breakdown(include_wavefunctions=False)
        # raise Exception(
        #     json.dumps(
        #         [
        #             (np.array(breakdown['frequencies']) * self.h2w).tolist(),
        #             breakdown['breakdowns']['Full']['intensities'].tolist(),
        #             breakdown['breakdowns']['Linear']['intensities'].tolist(),
        #             breakdown['breakdowns']['Quadratic']['intensities'].tolist(),
        #             breakdown['breakdowns']['Cubic']['intensities'].tolist(),
        #             breakdown['breakdowns']['Constant']['intensities'].tolist(),
        #             breakdown['breakdowns']['Order0']['intensities'].tolist(),
        #             breakdown['breakdowns']['Order1']['intensities'].tolist()
        #         ]
        #     )
        # )

        # trying to turn off the orthogonality condition
        # states = wfns.corrs.states
        # for i,s in enumerate(states):
        #     wfns.corrs.wfn_corrections[i, 2, s] = 0 # turn off second correction
        engs = h2w * wfns.energies
        freqs = engs - engs[0]

        with BlockProfiler("Intensities", print_res=True):
            ints = wfns.intensities

        harm_engs = h2w * wfns.zero_order_energies
        harm_freqs = harm_engs - harm_engs[0]
        harm_ints = wfns.zero_order_intensities

        print_specs = True
        if print_specs:
            report = "Harmonic:   State     Freq.   Int.  \n" + "\n".join(
                " " * 12 + "{} {:>7.2f} {:>7.4f}       ".format(s, f, i)
                for s, f, i in zip(states, harm_freqs, harm_ints)
            )
            print(report)
            report = "Anharmonic: State     Freq.   Int.\n" + "\n".join(
                " " * 12 + "{} {:>7.2f} {:>7.4f}  ".format(s, f, i)
                for s, f, i in zip(states, freqs, ints)
            )
            print(report)


    @validationTest
    def test_WaterDimerIntensitiesCartesian(self):

        internals = None

        n_modes = 6 * 3 - 6
        mode_selection = None
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)

        states = self.get_states(3, n_modes)


        with BlockProfiler('Water dimer wfns', print_res=True):
            wfns = self.get_VPT2_wfns(
                "water_dimer_freq.fchk",
                internals,
                states,
                regenerate=True,
                log=True,
                parallelized=True
            )

        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        # trying to turn off the orthogonality condition
        # states = wfns.corrs.states
        # for i,s in enumerate(states):
        #     wfns.corrs.wfn_corrections[i, 2, s] = 0 # turn off second correction
        engs = h2w * wfns.energies
        freqs = engs - engs[0]


        with BlockProfiler('Water dimer intensities', print_res=True):
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

        print_specs = True
        if print_specs:
            report = "Harmonic:   State     Freq.   Int.  \n" + "\n".join(
                " " * 12 + "{} {:>7.2f} {:>7.4f}         ".format(s, f, i)
                for s, f, i in zip(states, harm_freqs, harm_ints)
            )
            print(report)
            report = "Anharmonic: State     Freq.   Int. \n" + "\n".join(
                " " * 12 + "{} {:>7.2f} {:>7.4f} ".format(s, f, i)
                for s, f, i in zip(states, freqs, ints)
            )
            print(report)

    #endregion Test Intensities

    #region Intensity Breakdowns

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

        all_wfns = self.get_VPT2_wfns(
            "HOD_freq.fchk",
            internals,
            states,
            regenerate=True,
            get_breakdown=True
        )

        plot_spec = False
        with open('/Users/Mark/Desktop/HOD_freq.csv', "w+") as s:
            # with io.StringIO() as s:
            self.write_intensity_breakdown(s, all_wfns, plot_spec)

    @inactiveTest
    def test_HODCartesianIntensityBreakdown(self):
        internals = None
        states = (
            (0, 0, 0),
            (0, 0, 1), (0, 1, 0), (1, 0, 0),
            (0, 0, 2), (0, 2, 0), (2, 0, 0),
            (0, 1, 1), (1, 0, 1), (1, 1, 0)
        )

        all_wfns = self.get_VPT2_wfns(
            "HOD_freq.fchk",
            internals,
            states,
            regenerate=True,
            get_breakdown=True
        )

        plot_spec = False
        with open(os.path.expanduser('~/Desktop/HOD_freq_carts.csv'), "w+") as s:
            # with io.StringIO() as s:
            self.write_intensity_breakdown(s, all_wfns, plot_spec)

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

        all_wfns = self.get_VPT2_wfns(
            "HOH_freq.fchk",
            internals,
            states,
            regenerate=True,
            get_breakdown=True
            # , log=True
        )

        plot_spec = False
        with open(os.path.expanduser('~/Desktop/HOH_freq.csv'), "w+") as s:
            # with io.StringIO() as s:
            self.write_intensity_breakdown(s, all_wfns, plot_spec)

    @inactiveTest
    def test_HOHCartesianIntensityBreakdown(self):
        internals = None
        states = (
            (0, 0, 0),
            (0, 0, 1), (0, 1, 0), (1, 0, 0),
            (0, 0, 2), (0, 2, 0), (2, 0, 0),
            (0, 1, 1), (1, 0, 1), (1, 1, 0)
        )

        all_wfns = self.get_VPT2_wfns(
            "HOH_freq.fchk",
            internals,
            states,
            regenerate=True,
            get_breakdown=True
            # , log=True
        )

        plot_spec=False
        with open(os.path.expanduser('~/Desktop/HOH_freq_carts.csv'), "w+") as s:
        # with io.StringIO() as s:
            self.write_intensity_breakdown(s, all_wfns, plot_spec)

    @inactiveTest
    def test_WaterDimerCartesianIntensityBreakdown(self):

        internals = None

        n_modes = 6 * 3 - 6
        mode_selection = None
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)

        states = self.get_states(2, n_modes)

        all_wfns = self.get_VPT2_wfns(
            "water_dimer_freq.fchk",
            internals,
            states,
            regenerate=True,
            log=True,
            get_breakdown=True
        )

        plot_spec = False

        import io

        with open('/Users/Mark/Desktop/dimer_freq_cart.csv', "w+") as s:
        # with io.StringIO() as s:
            self.write_intensity_breakdown(s, all_wfns, plot_spec, write_wavefunctions=False)

    @inactiveTest
    def test_WaterDimerInternalIntensityBreakdown(self):

        """
         H                  ! 1: Shared H
         O 1 r1             ! 2: Right-side O
         H 2 r2 1 a1        ! 4: Right-side H
         O 2 r3 3 a2 1 t1 0 ! 5: LH O
         H 4 r4 3 a3 1 t2 0 ! 6: LH H
         H 4 r5 3 a4 1 t3 0 ! 7: LH H
        :return:
        :rtype:
        """
        COM = -3
        A = -2
        C = -1
        X = 1000
        LHF = 0
        LO = 1
        SH = 2
        RO = 3
        RH1 = 4
        RH2 = 5
        internals = [
            [ LHF,   X,   X,   X],
            [  LO, LHF,   X,   X],
            [  SH,  LO, LHF,   X],
            [  RO,  LO, LHF,   C],
            [ RH1,  RO,  SH, LHF],
            [ RH2,  RO, RH1, LHF]
        ]

        n_modes = 6 * 3 - 6
        mode_selection = None
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)

        states = self.get_states(2, n_modes)

        all_wfns = self.get_VPT2_wfns(
            "water_dimer_freq.fchk",
            internals,
            states,
            regenerate=True,
            log=True,
            get_breakdown=True
        )

        plot_spec = False

        import io

        with open('/Users/Mark/Desktop/dimer_freq.csv', "w+") as s:
            # with io.StringIO() as s:
            self.write_intensity_breakdown(s, all_wfns, plot_spec, write_wavefunctions=False)

    #endregion Intensity Breakdowns