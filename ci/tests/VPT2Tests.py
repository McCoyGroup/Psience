
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
                              , hamiltonian_options=None
                              , **solver_opts
                              ):
        if log is True:
            if verbose:
                log = Logger(log_level=LogLevel.All)
            else:
                log = Logger()
        if parallelized:
            parallelizer = MultiprocessingParallelizer(
                       logger=log,
                       processes=processes,
                       initialization_timeout=initialization_timeout
                       )
            # with parallelizer:
            #     raise Exception(parallelizer.comm, len(parallelizer.comm.locations))
        else:
            parallelizer = SerialNonParallelizer()

        if direct_sum_chunk_size is not None:
            BasisStateSpace.direct_sum_chunk_size = direct_sum_chunk_size

        if hamiltonian_options is None:
            hamiltonian_options = {}

        with parallelizer:
            if isinstance(mol_spec, str):
                hammer = PerturbationTheoryHamiltonian.from_fchk(
                    TestManager.test_data(mol_spec),
                    internals=internals,
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
                    include_coriolis_coupling=False if coriolis is False else True,
                    **hamiltonian_options
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
                    include_coriolis_coupling=False if coriolis is False else True,
                    **hamiltonian_options
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

    def print_energy_block(self, tag, wfns, states, zpe, freqs, real_fmt='{:>12.5f}'):

        if wfns is not None:
            print(
                tag,
                wfns.format_energies_table(
                    states=states, zpe=zpe, freqs=freqs,
                    real_fmt=real_fmt).replace(
                    "\n", "\n  "
                ),
                sep="\n  "
            )
        else:
            print(
                tag,
                PerturbationTheoryWavefunctions._format_energies_table(
                    states, zpe, freqs,
                    real_fmt=real_fmt
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
                    verbose=False,
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
                    hamiltonian_options=None,
                    **opts
                    ):
        if log and verbose:
            log = Logger(log_level=LogLevel.All)
        with BlockProfiler(tag, print_res=print_profile, mode=profiling_mode, inactive=not print_profile):#, filter=profile_filter):
            wfns, hammer = self.get_VPT2_wfns_and_ham(
                mol_spec,
                internals,
                states,
                regenerate=True,
                pre_run_script=pre_wfns_script,
                mode_selection=mode_selection
                , log=log,
                hamiltonian_options=hamiltonian_options,
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
            print("Energy Corrections:")
            print(wfns.format_energy_corrections_table())
            self.print_energy_block("State Energies:", wfns, states, my_energies, my_freqs)

        if print_diffs:
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
            ints = wfns.intensities
            print_specs = True
            if print_specs:
                if wfns.degenerate_transformation is not None:
                    print("Deperturbed Results")
                    if print_report:
                        for a, m in zip(["X", "Y", "Z"], wfns.format_deperturbed_dipole_contribs_tables()):
                            print("{} Dipole Contributions".format(a))
                            print(m)
                    print(wfns.format_deperturbed_intensities_table())
                    print("Degenerate Results")

                if print_report:
                    for a, m in zip(["X", "Y", "Z"], wfns.format_dipole_contribs_tables()):
                        print("{} Dipole Contributions".format(a))
                        print(m)
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

    #region Test Systems

    #region Analytic Models

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
    @validationTest #paper
    def test_OneMorseCartesiansNonDeg(self):

        import warnings
        np.seterr(all='raise')
        warnings.filterwarnings('error', category=np.VisibleDeprecationWarning)

        tag="OneMorseCartesiansNonDeg"

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
                                               stencil=7
                                               # parallelizer=MultiprocessingParallelizer()#verbose=True)
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
        states = VPTStateSpace.get_state_list_from_quanta(4, n_modes)

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
                [    re_1, 0.000000, 0.000000],
                [0.000000,     re_2, 0.000000],
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
        states = VPTStateSpace.get_state_list_from_quanta(4, n_modes)

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
        states = VPTStateSpace.get_state_list_from_quanta(4, n_modes)

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
    @validationTest
    def test_TwoMorseCartesiansAndBendNonDeg(self):

        import warnings
        np.seterr(all='raise')
        warnings.filterwarnings('error', category=np.VisibleDeprecationWarning)

        tag = "TwoMorseCartesiansNed"

        # Set up system


        cm2borr = UnitsData.convert("Angstroms", "BohrRadius")
        re_1 = 0.9575 * cm2borr
        re_2 = 0.9575 * cm2borr
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
                                               # parallelizer=MultiprocessingParallelizer()  # verbose=True)
                                               ).derivatives(mol.coords)


        # with BlockProfiler():
        mol.potential_derivatives = deriv_gen.derivative_tensor([1, 2, 3, 4])#, 5, 6])
            # raise Exception("wheee")
        # raise Exception([x.shape for x in mol.potential_derivatives])

        # we rigorously zero out small terms for
        # numerical stability
        # mat = mol.normal_modes.modes.basis.matrix
        # bad_spots = np.where(np.abs(mat) < 1e-14)
        # mat[bad_spots] = 0.

        n_atoms = 3
        n_modes = 3 * n_atoms - 6
        mode_selection = [-2, -1]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)

        # states = VPTStateSpace.get_state_list_from_quanta(5, n_modes)
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
            order=2,
            expansion_order=2
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
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

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
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

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
    @validationTest
    def test_HOHVPTInternals(self):

        tag = 'HOH Internals'
        file_name = "HOH_freq.fchk"

        # raise Exception(
        #     Molecule.from_file(TestManager.test_data("HOH_freq.fchk")).moments_of_inertia
        # )

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
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

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
            print_report=print_report,
            calculate_intensities=True
        )

    @validationTest
    def test_HOHVPTInternalsBackprop(self):

        tag = 'HOH Internals'
        file_name = "HOH_freq.fchk"

        # raise Exception(
        #     Molecule.from_file(TestManager.test_data("HOH_freq.fchk")).moments_of_inertia
        # )

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
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

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
            hamiltonian_options = {"backpropagate_internals":True},
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
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

        print_report = True

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
            verbose=False,
            print_report=print_report,
            calculate_intensities=True
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
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

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
            verbose=False,
            print_report=print_report,
            calculate_intensities=True,
            print_profile=False,
            print_x=False
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
    def test_HOHVPTCartesiansOnly2Quanta(self):

        tag = 'HOH Cartesians'
        file_name = "HOH_freq.fchk"

        internals = None

        n_atoms = 3
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = VPTStateSpace.get_state_list_from_quanta([0, 2], n_modes)

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
            verbose=False,
            print_report=print_report,
            calculate_intensities=True,
            print_profile=False,
            print_x=False
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
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

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
            log=False,
            verbose=False,
            print_report=print_report,
            calculate_intensities=True,
            state_space_filters=VPTStateSpace.get_state_space_filter(states, n_modes, target='intensities')
        )

    @validationTest
    def test_HOHVPTCartesiansSubsample(self):

        tag = 'HOH Cartesians'
        file_name = "HOH_freq.fchk"

        internals = None

        n_atoms = 3
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)[:3]

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
            log=False,
            verbose=False,
            print_report=print_report,
            calculate_intensities=True,
            print_profile=False,
            print_x=False
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
    def test_HOHVPTCartesiansSubsampleFiltered(self):

        tag = 'HOH Cartesians'
        file_name = "HOH_freq.fchk"

        internals = None

        n_atoms = 3
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)[:3]

        print_report = False

        gaussian_energies = self.gaussian_data['HOH']['zpe']
        gaussian_freqs = self.gaussian_data['HOH']['freqs']

        # import McUtils.Misc as mcmisc
        #
        # with mcmisc.without_numba():
        # os.remove(os.path.expanduser('~/hoh.hdf5'))
        # zero_to_ones = BasisStateSpace(
        #     HarmonicOscillatorProductBasis(n_modes),
        #     states[:1]

        # raise Exception(filt_space, one_qs)
        self.run_PT_test(
            tag,
            file_name,
            internals,
            mode_selection,
            states,
            gaussian_energies,
            gaussian_freqs,
            log=True,
            verbose=False,
            print_report=print_report,
            calculate_intensities=True,
            state_space_filters = VPTStateSpace.get_state_space_filter(states, n_modes, target='intensities')
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
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

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
            print_x=False,
            order=2,
            print_report=print_report,
            log=True,
            verbose=True,
            calculate_intensities=True
        )

    @validationTest
    def test_HOHVPTCartesiansSubspace(self):

        tag = 'HOH Cartesians'
        file_name = "HOH_freq.fchk"

        internals = None

        n_atoms = 3
        n_modes = 3 * n_atoms - 6
        mode_selection = [1, 2]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

        print_report = True

        gaussian_energies = None# self.gaussian_data['HOH']['zpe']
        gaussian_freqs = None# self.gaussian_data['HOH']['freqs']

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
            calculate_intensities=True,
            print_x=False
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
    def test_HOHVPTCartesians4thOrder(self):

        tag = 'HOH Cartesians'
        file_name = "HOH_freq.fchk"

        internals = None

        n_atoms = 3
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

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
            log=False,
            verbose=False,
            order=4,
            expansion_order=2,
            print_report=print_report
            , hamiltonian_options=dict(
                allow_higher_dipole_terms=True,
                # coordinate_derivatives={
                #     "int":[None, None, None, 0, 0],
                #     "cart":[None, None, None, 0, 0]
                # }
            )
            , calculate_intensities=True
            # , checkpoint=os.path.expanduser('~/hoh.hdf5'),
            # , watson=False
        )

    @validationTest
    def test_HOHVPTCartesians4thOrderKE(self):

        tag = 'HOH Cartesians'
        file_name = "HOH_freq.fchk"

        internals = None

        n_atoms = 3
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

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
            expansion_order=4,
            hamiltonian_options=dict(
                allow_higher_potential_terms=True,
                allow_higher_dipole_terms=True
            ),
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
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

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
            log=False,
            verbose=False,
            order=4,
            expansion_order=2,
            print_report=print_report
            , hamiltonian_options=dict(
                allow_higher_dipole_terms=True,
                # coordinate_derivatives={
                #     "int":[None, None, None, 0, 0],
                #     "cart":[None, None, None, 0, 0]
                # }
            )
            , calculate_intensities=True
        )

    @validationTest
    def test_HOHVPTInternals4thOrderKE(self):

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
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

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
            expansion_order=4,
            hamiltonian_options=dict(
                allow_higher_potential_terms=True,
                allow_higher_dipole_terms=True,
                # coordinate_derivatives={
                #     "int":[None, None, None, 0, 0],
                #     "cart":[None, None, None, 0, 0]
                # }
            ),
            print_report=print_report,
            calculate_intensities=True
        )

    @inactiveTest
    def test_HOHVPTCartesiansDegenerate(self):

        internals = None

        basis = HarmonicOscillatorProductBasis(3)

        states = VPTStateSpace.get_state_list_from_quanta(6, 3)

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
            #     BasisStateSpace(basis, VPTStateSpace.get_state_list_from_quanta(12, 3)).apply_selection_rules(
            #         basis.selection_rules("x", "x", "x")
            #     ),
            #     BasisStateSpace(basis, VPTStateSpace.get_state_list_from_quanta(12, 3)).apply_selection_rules(
            #         basis.selection_rules("x", "x", "x", "x")
            #     )
            # ]
            # , coupled_states=[
            #     BasisStateSpace(basis, VPTStateSpace.get_state_list_from_quanta(12, 3)),
            #     BasisStateSpace(basis, VPTStateSpace.get_state_list_from_quanta(12, 3))
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
            self.print_energy_block("Difference Energies:", states,
                                    my_energies - gaussian_energies,
                                    my_freqs[:ns] - gaussian_freqs[:ns]
                                    )

        self.assertLess(np.max(np.abs(my_freqs[:ns] - gaussian_freqs[:ns])), 1.5)

    @validationTest
    def test_HOHVPTRunner(self):

        file_name = "HOH_freq.fchk"
        VPTRunner.run_simple(
            TestManager.test_data(file_name),
            3,
            logger=True
        )

        system = VPTSystem(TestManager.test_data(file_name))
        states = VPTStateSpace.from_system_and_quanta(system, 3)
        pt_opts = VPTSolverOptions(state_space_filters=states.get_filter("intensities"))
        run_opts = VPTRuntimeOptions(logger=True)
        runner = VPTRunner(system, states, runtime_options=run_opts, solver_options=pt_opts)
        runner.print_tables()

    @validationTest
    def test_HOHVPTRunnerShifted(self):

        file_name = "HOH_freq.fchk"
        VPTRunner.run_simple(
            TestManager.test_data(file_name),
            3,
            logger=True,
            corrected_fundamental_frequencies=np.array([1600, 3775, 3880])/self.h2w
        )

    @validationTest
    def test_GetDegenerateSpaces(self):

        base_states = [
            [0, 0, 1],
            [0, 1, 0],
            [0, 2, 1],
            [0, 4, 0]
        ]

        degenerate_states = VPTStateSpace.get_degenerate_polyad_space(
            base_states,
            [
                [
                    [0, 2, 0],
                    [0, 0, 1]
                ]
            ],
        )

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
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

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
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

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
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

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
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

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
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

        gaussian_energies = self.gaussian_data['OCHH']['zpe']
        gaussian_freqs = self.gaussian_data['OCHH']['freqs']

        # ExpansionTerms.cartesian_analytic_deriv_order = 0
        # ExpansionTerms.cartesian_fd_mesh_spacing = 1.0e-3
        # ExpansionTerms.cartesian_fd_stencil = 9

        print_report = True
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
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

        gaussian_energies = self.gaussian_data['OCHH']['zpe']
        gaussian_freqs = self.gaussian_data['OCHH']['freqs']

        print_report = True
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
            gaussian_tolerance=gaussian_tolerance
        )

    @validationTest
    def test_OCHHVPTInternalsDummy(self):

        tag = 'OCHH Internals'
        file_name = TestManager.test_data("OCHH_freq.fchk")

        mol = Molecule.from_file(file_name)  # .get_embedded_molecule()
        o_pos = 0
        c_pos = 1
        h_pos = 2

        normal = -nput.vec_crosses(
            mol.coords[o_pos] - mol.coords[c_pos],
            mol.coords[h_pos] - mol.coords[c_pos],
            normalize=True
        )

        mol = mol.insert_atoms("X", mol.coords[c_pos] + 5 * normal, c_pos+1, handle_properties=True)

        # Dummy somehow broken...
        mol.zmatrix = [
            [0, -1, -1, -1],  # O
            [1,  0, -1, -1],  # C
            [2,  1,  0, -1],  # X
            [3,  1,  0,  2],  # H
            [4,  1,  0,  2],  # H
        ]

        internals = None

        n_atoms = 4
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

        gaussian_energies = self.gaussian_data['OCHH']['zpe']
        gaussian_freqs = self.gaussian_data['OCHH']['freqs']

        # ExpansionTerms.cartesian_analytic_deriv_order = 0
        # ExpansionTerms.cartesian_fd_mesh_spacing = 1.0e-3
        # ExpansionTerms.cartesian_fd_stencil = 9

        print_report = True
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
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

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
            state_space_filters = VPTStateSpace.get_state_space_filter(states, n_modes, target='intensities')
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
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

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

        # from McUtils.GaussianInterface import GaussianLogReader
        # with GaussianLogReader(TestManager.test_data('OCHH_freq_16.log')) as reader:
        #     sort_spec = np.flip([4, 0, 1, 2, 5, 3])
        #     x = reader.parse("XMatrix")["XMatrix"][np.ix_(sort_spec, sort_spec)]
        #     # raise Exception(x)
        # def pre_wfns_script(hammer, states):
        #     harm, corrs = hammer.get_Nielsen_energies(states, x_mat=x)
        #     harm = harm * self.h2w
        #     anharms = harm + corrs
        #     # print(states)
        #     print(np.column_stack([
        #         harm[1:] - harm[0],
        #         anharms[1:] - anharms[0]
        #     ]))
        # pre_wfns_script = None

        degeneracies = VPTStateSpace.get_degenerate_polyad_space(
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
            # pre_wfns_script=pre_wfns_script,
            log=True,
            verbose='all',
            calculate_intensities=True,
            gaussian_tolerance=gaussian_tolerance,
            degeneracies=degeneracies,
            gaussian_resonance_handling=True
            , state_space_filters = VPTStateSpace.get_state_space_filter(states, n_modes, target='intensities')
        )

    @debugTest
    def test_OCHHVPTCartesiansNonDegenerateSubsample2(self):

        tag = 'OCHH Cartesians'
        file_name = "OCHH_freq.fchk"

        internals = None

        n_atoms = 4
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = VPTStateSpace.get_state_list_from_quanta([0, 2], n_modes)[:5]

        degeneracies = VPTStateSpace.get_degenerate_polyad_space(
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

        degeneracies = None

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
            # pre_wfns_script=pre_wfns_script,
            log=True,
            verbose='all',
            calculate_intensities=True,
            gaussian_tolerance=gaussian_tolerance,
            degeneracies=degeneracies,
            gaussian_resonance_handling=True
            , state_space_filters=VPTStateSpace.get_state_space_filter(states, n_modes, target='intensities')
        )

    @debugTest
    def test_OCHHVPTCartesiansDegenerateSubsample2(self):

        tag = 'OCHH Cartesians'
        file_name = "OCHH_freq.fchk"

        internals = None

        n_atoms = 4
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = VPTStateSpace.get_state_list_from_quanta([0, 2], n_modes)[:5]

        degeneracies = VPTStateSpace.get_degenerate_polyad_space(
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
            # pre_wfns_script=pre_wfns_script,
            log=True,
            verbose='all',
            calculate_intensities=True,
            gaussian_tolerance=gaussian_tolerance,
            degeneracies=degeneracies,
            gaussian_resonance_handling=True
            , state_space_filters=VPTStateSpace.get_state_space_filter(states, n_modes, target='intensities')
        )

    @validationTest
    def test_OCHHVPTCartesiansDegenerateSubsampleFiltered(self):

        tag = 'OCHH Cartesians'
        file_name = "OCHH_freq.fchk"

        internals = None

        n_atoms = 4
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = VPTStateSpace.get_state_list_from_quanta(4, n_modes)
        np.random.seed(0)
        subsel = np.unique(np.random.randint(len(VPTStateSpace.get_state_list_from_quanta(3, n_modes)), len(states), 10))
        states = states[:1] + [states[s] for s in subsel] + VPTStateSpace.get_state_list_from_quanta(2, n_modes)[1:]

        degeneracies = VPTStateSpace.get_degenerate_polyad_space(
            states,
            [
                # [[0, 0, 0, 0, 0, 1], [0, 1, 0, 1, 0, 0]],
                [[0, 0, 1, 0, 0, 0], [0, 0, 0, 1, 0, 0]],
            ]
        )
        degeneracies = [np.array(x).tolist() for x in degeneracies]
        states = np.array(states).tolist()
        flat_degs = []
        for pair in degeneracies:
            for p in pair:
                if p not in states:
                    states.append(p)
                if p not in flat_degs:
                    flat_degs.append(p)

        # degeneracies = None

        gaussian_energies = self.gaussian_data['OCHH']['zpe']
        gaussian_freqs = self.gaussian_data['OCHH']['freqs']

        # assert [0, 1, 1, 0, 1, 0] in filter_space
        # assert [0, 1, 0, 1, 1, 0] in filter_space
        # sorter = np.argsort(HarmonicOscillatorProductBasis(n_modes).ravel_state_inds(filter_space))
        # filter_space = [filter_space[s] for s in sorter]

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
            # pre_wfns_script=pre_wfns_script,
            log=True,
            verbose=True,
            calculate_intensities=True,
            gaussian_tolerance=gaussian_tolerance,
            degeneracies=degeneracies,
            gaussian_resonance_handling=False
            # these filters work fine for _three_ quantum resonances it seems?
            # as we push out of that space though I think I need more stuff?
            , state_space_filters = VPTStateSpace.get_state_space_filter(states, n_modes, target='intensities')
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
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

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

        degeneracies = VPTStateSpace.get_degenerate_polyad_space(
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

    @validationTest
    def test_OCHHVPTRunner(self):

        file_name = "OCHH_freq.fchk"
        VPTRunner.run_simple(
            TestManager.test_data(file_name),
            3,
            logger=True,
            degeneracy_specs=[
                [[0, 0, 0, 0, 0, 1], [0, 1, 0, 1, 0, 0]]
            ]
        )

    @validationTest
    def test_OCHHVPTRunnerPolyadExtended(self):

        file_name = "OCHH_freq.fchk"
        VPTRunner.run_simple(
            TestManager.test_data(file_name),
            3,
            logger=True,
            degeneracy_specs=[
                    [[0, 0, 0, 0, 0, 1], [0, 1, 0, 1, 0, 0]]
                ]
        )

        VPTRunner.run_simple(
            TestManager.test_data(file_name),
            3,
            logger=Logger(log_level=Logger.LogLevel.All),
            degeneracy_specs={
                "polyads": [
                    [[0, 0, 0, 0, 0, 1], [0, 1, 0, 1, 0, 0]]
                ],
                "extra_groups":[
                    [[0, 2, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1]]
                ]
            }
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
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

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
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

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
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

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
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

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
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

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
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

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
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

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
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

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
        # ExpansionTerms.cartesian_fd_mesh_spacing = 1.0e-7
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
            hamiltonian_options=dict(cartesian_analytic_deriv_order=0)
        )
    @validationTest
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

        mol = mol.insert_atoms("X", mol.coords[o_pos[1]] + 5 * normal, 5, handle_properties=True)

        # Dummy doing nothing
        # mol.zmatrix = [
        #     [1, -1, -1, -1],
        #     [2,  1, -1, -1],
        #     [3,  2,  1, -1],
        #     [0,  1,  2,  3],
        #     [4,  3,  2,  1],
        #     [5,  2,  1,  3]
        # ]

        # #Dummy somehow broken...
        # mol.zmatrix = [
        #     [1, -1, -1, -1],  # O
        #     [2,  1, -1, -1],  # O
        #     [3,  2,  1, -1],  # N
        #     [5,  2,  1,  3],  # X
        #     [0,  1,  2,  5],  # H
        #     [4,  3,  2,  5]   # O
        # ]
        # ExpansionTerms.cartesian_analytic_deriv_order = 0

        # # With only ref to H somehow clean...
        # mol.zmatrix = [
        #     [1, -1, -1, -1],  # O
        #     [2,  1, -1, -1],  # O
        #     [3,  2,  1, -1],  # N
        #     [5,  2,  1,  3],  # X
        #     [0,  1,  2,  5],  # H
        #     [4,  3,  2,  1]   # O
        # ]
        # ExpansionTerms.cartesian_analytic_deriv_order = 0

        # # Broken when O refs dummy...?
        # mol.zmatrix = [
        #     [1, -1, -1, -1],  # O
        #     [2,  1, -1, -1],  # O
        #     [3,  2,  1, -1],  # N
        #     [5,  2,  1,  3],  # X
        #     [0,  1,  2,  3],  # H
        #     [4,  3,  2,  5]   # O
        # ]
        # ExpansionTerms.cartesian_analytic_deriv_order = 0

        ## Fine if O defined against original OO pair...?
        mol.zmatrix = [
            [1, -1, -1, -1],  # O
            [2,  1, -1, -1],  # O
            [5,  2,  1, -1],  # X
            [3,  2,  1,  5],  # N
            [4,  1,  2,  5],  # O
            [0,  1,  2,  5],  # H
        ]

        n_atoms = 5
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

        gaussian_energies = self.gaussian_data['HOONO']['zpe']
        gaussian_freqs = self.gaussian_data['HOONO']['freqs']

        print_report = True
        nielsen_tolerance = None
        gaussian_tolerance = 10

        internals = None
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
            print_report=print_report,
            calculate_intensities=True,
            nielsen_tolerance=nielsen_tolerance,
            gaussian_tolerance=gaussian_tolerance
            , state_space_filters = VPTStateSpace.get_state_space_filter(states, n_modes, target='intensities')
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
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

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
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

        gaussian_energies = self.gaussian_data['HOONO']['zpe']
        gaussian_freqs = self.gaussian_data['HOONO']['freqs']

        pre_wfns_script = None

        degeneracies = None
        print_report = True
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
            print_x=False,
            calculate_intensities=True,
            degeneracies=degeneracies,
            print_report=print_report,
            nielsen_tolerance=nielsen_tolerance,
            gaussian_tolerance=gaussian_tolerance,
            pre_wfns_script=pre_wfns_script
            , state_space_filters = VPTStateSpace.get_state_space_filter(states, n_modes, target='intensities')
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
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

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
            # print_x=True
            # zero_element_warning=False

        )

    @validationTest
    def test_HOONOVPTInternalsDegenerate(self):

        tag = 'HOONO Internals'
        file_name = TestManager.test_data("HOONO_freq.fchk")

        n_atoms = 5
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

        # internals = [
        #     [0, -1, -1, -1],  # H
        #     [1,  0, -1, -1],  # O
        #     [2,  1,  0, -1],  # O
        #     [3,  2,  1,  0],  # N
        #     [4,  3,  2,  0]   # O
        # ]
        internals = [
            [1, -1, -1, -1],
            [2, 1, -1, -1],
            [3, 2, 1, -1],
            [0, 1, 2, 3],
            [4, 3, 2, 1]
        ]

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

        nt_spec = np.array([
            [0, 0, 0, 1, 0, 0, 0, 0, 0],
            [2, 0, 0, 0, 0, 0, 0, 0, 0]
        ])
        degeneracies2 = [
            [
                s,
                (s - nt_spec[0] + nt_spec[1])
            ] for s in states if np.dot(s, nt_spec[0]) > 0
        ]
        degeneracies += [np.array(x).tolist() for x in degeneracies2]
        states = np.array(states).tolist()
        for pair in degeneracies:
            for p in pair:
                if p not in states:
                    states.append(p)

        gaussian_energies = self.gaussian_data['HOONO']['zpe']
        gaussian_freqs = self.gaussian_data['HOONO']['freqs']

        print_report = True
        nielsen_tolerance = 10
        gaussian_tolerance = 10
        # from Psience.VPT2 import PotentialTerms
        # PotentialTerms.hessian_tolerance = None
        ExpansionTerms.cartesian_analytic_deriv_order = 0
        # ExpansionTerms.cartesian_fd_mesh_spacing = 1.0e-7
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
            state_space_filters = VPTStateSpace.get_state_space_filter(states, n_modes, target='intensities')

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
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

        gaussian_energies = self.gaussian_data['HOONO']['zpe']
        gaussian_freqs = self.gaussian_data['HOONO']['freqs']

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
        nt_spec = np.array([
            [0, 0, 0, 1, 0, 0, 0, 0, 0],
            [2, 0, 0, 0, 0, 0, 0, 0, 0]
        ])
        degeneracies2 = [
            [
                s,
                (s - nt_spec[0] + nt_spec[1])
            ] for s in states if np.dot(s, nt_spec[0]) > 0
        ]
        degeneracies += [np.array(x).tolist() for x in degeneracies2]
        states = np.array(states).tolist()
        for pair in degeneracies:
            for p in pair:
                if p not in states:
                    states.append(p)

        print_report = True
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
            gaussian_resonance_handling=False,
            state_space_filters = VPTStateSpace.get_state_space_filter(states, n_modes, target='intensities')
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
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)

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
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)  # [:6]

        gaussian_energies = self.gaussian_data['WaterDimer']['zpe']
        gaussian_freqs = self.gaussian_data['WaterDimer']['freqs']

        # chk = os.path.expanduser('~/Desktop/dimer_chk2.hdf5')
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
            print_profile=False,
            # profile_filter='Combinatorics/Permutations',
            print_report=print_report,
            nielsen_tolerance=nielsen_tolerance,
            gaussian_tolerance=gaussian_tolerance,
            calculate_intensities=True
            # , checkpoint=chk
            , use_cached_representations=False
            , state_space_filters = VPTStateSpace.get_state_space_filter(states, n_modes, target='intensities')
            # , parallelized=True
        )

    @validationTest
    def test_WaterDimerVPTInternalsBackprop(self):

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
            [LHF, X, X, X],
            [LO, LHF, X, X],
            [SH, LO, LHF, X],
            [RO, LO, LHF, C],
            [RH1, RO, SH, LHF],
            [RH2, RO, RH1, LHF]
        ]

        n_atoms = 6
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)  # [:6]

        gaussian_energies = self.gaussian_data['WaterDimer']['zpe']
        gaussian_freqs = self.gaussian_data['WaterDimer']['freqs']

        # chk = os.path.expanduser('~/Desktop/dimer_chk2.hdf5')
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
            verbose=False,
            print_profile=False,
            hamiltonian_options = {"backpropagate_internals":True},
            # profile_filter='Combinatorics/Permutations',
            print_report=print_report,
            nielsen_tolerance=nielsen_tolerance,
            gaussian_tolerance=gaussian_tolerance,
            calculate_intensities=True
            # , checkpoint=chk
            , use_cached_representations=False
            , state_space_filters = VPTStateSpace.get_state_space_filter(states, n_modes, target='intensities')
            # , parallelized=True
        )

    @validationTest
    def test_WaterDimerVPTInternalsDummy(self):
        # Borked up somehow...

        tag = 'Water Dimer Internals'
        file_name = TestManager.test_data("water_dimer_freq.fchk")

        COM = -3
        A = -2
        C = -1
        X = 1000
        LHF = 0
        LO = 1
        SH = 2

        mol = Molecule.from_file(file_name)  # .get_embedded_molecule()

        normal = -nput.vec_crosses(
            mol.coords[LO] - mol.coords[SH],
            mol.coords[LHF] - mol.coords[SH],
            normalize=True
        )

        mol = mol.insert_atoms("X", mol.coords[SH] + 5 * normal, SH + 1, handle_properties=True)

        Z = 3
        RO = 3 + 1
        RH1 = 4 + 1
        RH2 = 5 + 1

        # Dummy somehow broken...
        mol.zmatrix = [
            [LHF,   X,   X,    X],
            [LO,  LHF,   X,    X],
            [SH,   LO,  LHF,   X],
            [ Z,   SH,  LHF,  LO],
            [RO,   LO,    Z, LHF],
            [RH1,  RO,   SH, LHF],
            [RH2,  RO,  RH1, LHF]
        ]

        # mol.plot()[0].show()

        n_atoms = 6
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)  # [:6]

        gaussian_energies = self.gaussian_data['WaterDimer']['zpe']
        gaussian_freqs = self.gaussian_data['WaterDimer']['freqs']

        # chk = os.path.expanduser('~/Desktop/dimer_chk2.hdf5')
        print_report = True
        nielsen_tolerance = None
        gaussian_tolerance = 50
        self.run_PT_test(
            tag,
            mol,
            None,
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
            , state_space_filters = VPTStateSpace.get_state_space_filter(states, n_modes, target='intensities')
            # , parallelized=True
        )

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
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)#[:6]

        gaussian_energies = self.gaussian_data['WaterDimer']['zpe']
        gaussian_freqs = self.gaussian_data['WaterDimer']['freqs']

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
            verbose=False,
            print_profile=False,
            # profile_filter='Combinatorics/Permutations',
            print_report=print_report,
            nielsen_tolerance=nielsen_tolerance,
            gaussian_tolerance=gaussian_tolerance,
            calculate_intensities=True
            # , checkpoint=chk
            , use_cached_representations=False
            , state_space_filters = VPTStateSpace.get_state_space_filter(states, n_modes, target='intensities')
        )

    @validationTest
    def test_WaterDimerVPTCartesiansDegenerate(self):

        tag = 'Water Dimer Cartesians'
        file_name = "water_dimer_freq.fchk"

        internals = None

        n_atoms = 6
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)  # [:6]

        np.random.seed(0)
        subsel = np.unique(np.random.randint(len(VPTStateSpace.get_state_list_from_quanta(2, n_modes)), len(states), 20))
        states = VPTStateSpace.get_state_list_from_quanta(1, n_modes) + [states[s] for s in subsel]

        degeneracies = VPTStateSpace.get_degenerate_polyad_space(
            states,
            [
                # [[0, 0, 0, 0, 0, 1], [0, 1, 0, 1, 0, 0]],
                [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0]],
            ]
        )
        degeneracies = [np.array(x).tolist() for x in degeneracies]
        states = np.array(states).tolist()
        flat_degs = []
        for pair in degeneracies:
            for p in pair:
                if p not in states:
                    states.append(p)
                if p not in flat_degs:
                    flat_degs.append(p)

        gaussian_energies = self.gaussian_data['WaterDimer']['zpe']
        gaussian_freqs = self.gaussian_data['WaterDimer']['freqs']

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
            degeneracies=degeneracies,
            log=True,
            verbose=True,
            print_profile=False,
            # profile_filter='Combinatorics/Permutations',
            print_report=print_report,
            nielsen_tolerance=nielsen_tolerance,
            gaussian_tolerance=gaussian_tolerance,
            calculate_intensities=True
            # , checkpoint=chk
            , use_cached_representations=False
            , state_space_filters=VPTStateSpace.get_state_space_filter(states, n_modes, target='intensities')
        )

    @validationTest
    def test_WaterDimerVPTCartesiansSubterms(self):
        # the high-frequency stuff agrees with Gaussian, but not the low-freq

        tag = 'Water Dimer Cartesians'
        file_name = "water_dimer_freq.fchk"

        internals = None

        n_atoms = 6
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)  # [:6]

        gaussian_energies = self.gaussian_data['WaterDimer']['zpe']
        gaussian_freqs = self.gaussian_data['WaterDimer']['freqs']

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
            verbose=False,
            print_profile=False,
            # profile_filter='Combinatorics/Permutations',
            print_report=print_report,
            nielsen_tolerance=nielsen_tolerance,
            gaussian_tolerance=gaussian_tolerance,
            calculate_intensities=True
            # , checkpoint=chk
            , use_cached_representations=False
            , state_space_terms=((1, 0), (2, 0))
        )

    @validationTest
    def test_WaterDimerVPTCartesiansRestrictedFilter(self):
        # the high-frequency stuff agrees with Gaussian, but not the low-freq

        tag = 'Water Dimer Cartesians'
        file_name = "water_dimer_freq.fchk"

        internals = None

        n_atoms = 6
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = VPTStateSpace.get_state_list_from_quanta(2, n_modes)[:5]

        gaussian_energies = self.gaussian_data['WaterDimer']['zpe']
        gaussian_freqs = self.gaussian_data['WaterDimer']['freqs']

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
            verbose=False,
            print_profile=False,
            # profile_filter='Combinatorics/Permutations',
            print_report=print_report,
            nielsen_tolerance=nielsen_tolerance,
            gaussian_tolerance=gaussian_tolerance,
            calculate_intensities=True
            # , checkpoint=chk
            , use_cached_representations=False
            , state_space_filters = VPTStateSpace.get_state_space_filter(states, n_modes, 'intensities')
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
        states = VPTStateSpace.get_state_list_from_quanta(2, n_modes)

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
            verbose=False,
            print_profile=True,
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
        states = VPTStateSpace.get_state_list_from_quanta(2, n_modes)

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

        tag = 'Water Dimer Cartesians'
        file_name = "water_dimer_freq.fchk"

        internals = None

        n_atoms = 6
        n_modes = 3 * n_atoms - 6
        mode_selection = [-3, -2, -1]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = VPTStateSpace.get_state_list_from_quanta(2, n_modes)  # [:6]

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
            log=True,
            verbose=True,
            print_profile=True,
            # profile_filter='Combinatorics/Permutations',
            print_report=print_report,
            nielsen_tolerance=nielsen_tolerance,
            gaussian_tolerance=gaussian_tolerance,
            calculate_intensities=True
            # , checkpoint=chk
            , use_cached_representations=False
            , state_space_filters = VPTStateSpace.get_state_space_filter(states, n_modes, target='intensities')
            # , parallelized=True
        )

    @validationTest
    def test_WaterDimerVPTCartesiansHarmonic(self):
        # the high-frequency stuff agrees with Gaussian, but not the low-freq

        tag = 'Water Dimer Cartesians'
        file_name = "water_dimer_freq.fchk"

        internals = None

        n_atoms = 6
        n_modes = 3 * n_atoms - 6
        mode_selection = None  # [5, 4, 3]
        if mode_selection is not None and len(mode_selection) < n_modes:
            n_modes = len(mode_selection)
        states = VPTStateSpace.get_state_list_from_quanta(2, n_modes)#[:6]

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
            order=0,
            log=True,
            verbose=True,
            print_profile=True,
            # profile_filter='Combinatorics/Permutations',
            print_report=print_report,
            nielsen_tolerance=nielsen_tolerance,
            gaussian_tolerance=gaussian_tolerance,
            calculate_intensities=True
            # , checkpoint=chk
            , use_cached_representations=False
            , state_space_filters=VPTStateSpace.get_state_space_filter(states, n_modes, 'intensities')
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

        states = VPTStateSpace.get_state_list_from_quanta(3, n_modes)  # [:6]

        # np.random.seed(0)
        # subsel = np.unique(np.random.randint(len(VPTStateSpace.get_state_list_from_quanta(2, n_modes)), len(states), 20))
        # subsel2 = np.unique(np.random.randint(len(VPTStateSpace.get_state_list_from_quanta(1, n_modes)), len(VPTStateSpace.get_state_list_from_quanta(2, n_modes))-1, 20))
        # states = VPTStateSpace.get_state_list_from_quanta(1, n_modes) + [states[s] for s in subsel2]  + [states[s] for s in subsel]
        # raise Exception(states)

        gaussian_energies = None#self.gaussian_data['WaterDimer']['zpe']
        gaussian_freqs = None#self.gaussian_data['WaterDimer']['freqs']

        print_report = False
        print_diffs = False
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
            verbose=False,
            print_profile=False,
            # profile_filter='Combinatorics/Permutations',
            print_report=print_report,
            print_diffs=print_diffs,
            nielsen_tolerance=nielsen_tolerance,
            gaussian_tolerance=gaussian_tolerance,
            calculate_intensities=True
            # , checkpoint=chk
            , use_cached_representations=False
            , state_space_filters= VPTStateSpace.get_state_space_filter(states, n_modes, 'intensities')
            # , parallelized=True
            # , processes=5
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