"""
A little package of utilities for setting up/running VPT jobs
"""

import numpy as np

from ..BasisReps import BasisStateSpace, HarmonicOscillatorProductBasis
from ..Molecools import Molecule

from .Hamiltonian import PerturbationTheoryHamiltonian

__all__ = [
    "VPTRunner"
]

class VPTRunner:
    """
    A helper class to make it easier to run jobs by making the inputs/options
    clear and making it easier to customize run options
    """

    class InputSystem:
        """
        Provides a little helper for setting up the input
        system for a VPT job
        """

        def __init__(self,
                     mol,
                     internals=None,
                     modes=None,
                     mode_selection=None,
                     potential_derivatives=None,
                     dipole_derivatives=None
                     ):
            """
            :param mol: the molecule or system specification to use (doesn't really even need to be a molecule)
            :type mol: str | Molecule
            :param internals: the Z-matrix for the internal coordinates (in the future will support a general function for this too)
            :type internals:
            :param modes: the normal modes to use if not already supplied by the Molecule
            :type modes:
            :param mode_selection: the subset of normal modes to do perturbation theory on
            :type mode_selection:
            :param potential_derivatives: the derivatives of the potential to use for expansions
            :type potential_derivatives: Iterable[np.ndarray]
            :param dipole_derivatives: the set of dipole derivatives to use for expansions
            :type dipole_derivatives: Iterable[np.ndarray]
            """
            if not isinstance(mol, Molecule):
                mol = self.load_molecule_from_spec(mol)
            mol.zmatrix = internals
            self.mol = mol
            if modes is not None:
                self.mol.normal_modes = modes
            if mode_selection is not None:
                self.mol.normal_modes.modes = self.mol.normal_modes.modes[mode_selection]
            if potential_derivatives is not None:
                self.potential_derivatives = potential_derivatives
            if dipole_derivatives is not None:
                self.dipole_derivatives = dipole_derivatives

        @classmethod
        def load_molecule_from_spec(cls, spec):
            if isinstance(spec, str):
                mol = Molecule.from_file(spec)
            else:
                raise NotImplementedError("don't have molecule loading from spec {}".format(spec))

            return mol

        @classmethod
        def from_harmonic_scan(cls, scan_array):
            raise NotImplementedError("loading from a Harmonic scan needs to be implemented")

    class HamiltonianOptions:
        """
        Provides a helper to keep track of the levers available for
        setting up the Hamiltonian
        """

        def __init__(self,
                     potential_derivatives=None,
                     coriolis_coupling=None,
                     include_pseudopotential=None,
                     potential_terms=None,
                     kinetic_terms=None,
                     coriolis_terms=None,
                     pseudopotential_terms=None,
                     undimensionalize_normal_modes=None,
                     use_numerical_jacobians=None,
                     eckart_embed_derivatives=None,
                     strip_dummy_atoms=None,
                     strip_embedding_coordinates=None,
                     mixed_derivative_handling_mode=None,
                     backpropagate_internals=None,
                     zero_mass_term=None,
                     internal_fd_mesh_spacing=None,
                     internal_fd_stencil=None,
                     cartesian_fd_mesh_spacing=None,
                     cartesian_fd_stencil=None,
                     cartesian_analytic_deriv_order=None,
                     internal_by_cartesian_order=None,
                     cartesian_by_internal_order=None,
                     jacobian_warning_threshold=None,
                     check_input_force_constants=None,
                     hessian_tolerance=None,
                     grad_tolerance=None,
                     freq_tolerance=None,
                     g_derivative_threshold=None
                     ):
            """
            :param potential_derivatives: the potential derivatives to use when building the expansions
            :type potential_derivatives: Iterable[np.ndarray]
            :param coriolis_coupling: whether or not to include Coriolis coupling in Cartesian normal mode calculation
            :type coriolis_coupling: bool
            :param include_pseudopotential: whether or not to include the pseudopotential/Watson term
            :type include_pseudopotential: bool
            :param potential_terms: explicit values for the potential terms (e.g. from analytic models)
            :type potential_terms: Iterable[np.ndarray]
            :param kinetic_terms: explicit values for the kinetic terms (e.g. from analytic models)
            :type kinetic_terms: Iterable[np.ndarray]
            :param coriolis_terms: explicit values for the Coriolis terms
            :type coriolis_terms: Iterable[np.ndarray]
            :param pseudopotential_terms: explicit values for the psuedopotential terms
            :type pseudopotential_terms: Iterable[np.ndarray]
            :param undimensionalize_normal_modes: whether or not to convert normal modes into dimensional coordinates
            :type undimensionalize_normal_modes: bool
            :param use_numerical_jacobians: whether or not to use numerical differentiation when getting coordinate transformations
            :type use_numerical_jacobians: bool
            :param eckart_embed_derivatives: whether or not to use Eckart embedding when getting Cartesian to internal transformations (needed for proper results)
            :type eckart_embed_derivatives: bool
            :param strip_dummy_atoms: whether or not to strip off dummy atoms when doing transformations
            :type strip_dummy_atoms: bool
            :param strip_embedding_coordinates: whether or not to strip off translation/rotation embedding coordinates when doing transformations
            :type strip_embedding_coordinates: bool
            :param mixed_derivative_handling_mode: how to handle differences between numerical/analytical mixed derivatives of potential/dipole terms
            :type mixed_derivative_handling_mode: bool
            :param backpropagate_internals: whether or not to do Cartesian coordinate calculations with values backpropagated from internals
            :type backpropagate_internals: bool
            :param zero_mass_term: a placeholder value for dummy atom masses
            :type zero_mass_term: float
            :param internal_fd_mesh_spacing: mesh spacing for finite difference of Cartesian coordinates with internals
            :type internal_fd_mesh_spacing: float
            :param internal_fd_stencil: stencil for finite difference of Cartesian coordinates with internals
            :type internal_fd_stencil: int
            :param cartesian_fd_mesh_spacing: mesh spacing for finite difference of internal coordinates with respect to Cartesians
            :type cartesian_fd_mesh_spacing: float
            :param cartesian_fd_stencil: stencil for finite difference of internal coordinates with respect to Cartesians
            :type cartesian_fd_stencil: int
            :param cartesian_analytic_deriv_order: order of analytic derivatives to use for derivatives of internal coordinates with respect to Cartesians (supports `0` or `1`)
            :type cartesian_analytic_deriv_order: int
            :param jacobian_warning_threshold: the value at which to warn that the Jacobian is ill-conditions
            :type jacobian_warning_threshold: float
            :param check_input_force_constants: whether or not to check that the input force constants match the input frequencies
            :type check_input_force_constants: bool
            :param hessian_tolerance: the deviation to allow when transforming from Cartesian to internal Hessian
            :type hessian_tolerance: float
            :param grad_tolerance: the size of the norm of the gradient above which to print a warning
            :type grad_tolerance: float
            :param freq_tolerance: the deviation from the input frequencies to allow when transforming from Cartesians to internals
            :type freq_tolerance: float
            :param g_derivative_threshold: the size of the norm of any G-matrix derivative above which to print a warning
            :type g_derivative_threshold: float
            """
            all_opts = dict(
                potential_derivatives=potential_derivatives,
                coriolis_coupling=coriolis_coupling,
                include_pseudopotential=include_pseudopotential,
                potential_terms=potential_terms,
                kinetic_terms=kinetic_terms,
                coriolis_terms=coriolis_terms,
                pseudopotential_terms=pseudopotential_terms,
                undimensionalize=undimensionalize_normal_modes,
                numerical_jacobians=use_numerical_jacobians,
                eckart_embed_derivatives=eckart_embed_derivatives,
                strip_dummies=strip_dummy_atoms,
                strip_embedding=strip_embedding_coordinates,
                mixed_derivative_handling_mode=mixed_derivative_handling_mode,
                backpropagate_internals=backpropagate_internals,
                zero_mass_term=zero_mass_term,
                internal_fd_mesh_spacing=internal_fd_mesh_spacing,
                internal_fd_stencil=internal_fd_stencil,
                cartesian_fd_mesh_spacing=cartesian_fd_mesh_spacing,
                cartesian_fd_stencil=cartesian_fd_stencil,
                cartesian_analytic_deriv_order=cartesian_analytic_deriv_order,
                internal_by_cartesian_order=internal_by_cartesian_order,
                cartesian_by_internal_order=cartesian_by_internal_order,
                jacobian_warning_threshold=jacobian_warning_threshold,
                check_input_force_constants=check_input_force_constants,
                hessian_tolerance=hessian_tolerance,
                grad_tolerance=grad_tolerance,
                freq_tolerance=freq_tolerance,
                g_derivative_threshold=g_derivative_threshold
            )

            real_opts = {}
            for o, v in all_opts.items():
                if v is not None:
                    real_opts[o] = v

            self.opts = real_opts

    class RuntimeOptions:
        """
        Provides a helper to keep track of the options available
        for configuring the way the code runs
        """

        def __init__(self,
                     operator_chunk_size=None,
                     logger=None,
                     verbose=None,
                     checkpoint=None,
                     parallelizer=None,
                     memory_constrained=None,
                     checkpoint_keys=None,
                     use_cached_representations=None,
                     use_cached_basis=None
                     ):
            """
            :param operator_chunk_size: the number of representation matrix elements to calculate at once
            :type operator_chunk_size: int
            :param logger: the `Logger` object to use when logging the status of the calculation
            :type logger: Logger
            :param verbose: whether or not to be verbose in log output
            :type verbose: bool
            :param checkpoint: the checkpoint file or `Checkpointer` object to use
            :type checkpoint: str
            :param parallelizer: the `Parallelizer` object to use when parallelizing pieces of the calculation
            :type parallelizer: Parallelizer
            :param memory_constrained: whether or not to attempt memory optimizations
            :type memory_constrained: bool
            :param checkpoint_keys: the keys to write to the checkpoint file
            :type checkpoint_keys: Iterable[str]
            :param use_cached_representations: whether or not to try to load representation matrices from the checkpoint
            :type use_cached_representations: bool
            :param use_cached_basis: whether or not to try to load the bases to use from the checkpoint
            :type use_cached_basis: bool
            """
            ham_run_opts = dict(
                operator_chunk_size=operator_chunk_size,
                logger=logger,
                checkpoint=checkpoint,
                parallelizer=parallelizer
            )
            real_ham_opts = {}
            for o, v in ham_run_opts.items():
                if v is not None:
                    real_ham_opts[o] = v
            self.ham_opts = real_ham_opts

            solver_run_opts = dict(
                operator_chunk_size=operator_chunk_size,
                memory_constrained=memory_constrained,
                verbose=verbose,
                checkpoint_keys=checkpoint_keys,
                use_cached_representations=use_cached_representations,
                use_cached_basis=use_cached_basis
            )
            real_solver_run_opts = {}
            for o, v in solver_run_opts.items():
                if v is not None:
                    real_solver_run_opts[o] = v
            self.solver_opts = real_solver_run_opts

    class PerturbationTheoryOptions:
        """
        Provides a helper to keep track of the options available
        for configuring the way the perturbation theory is applied
        """

        def __init__(self,
                     coupled_states=None,
                     order=2,
                     expansion_order=None,
                     total_space=None,
                     flat_total_space=None,
                     state_space_iterations=None,
                     state_space_terms=None,
                     state_space_filters=None,
                     allow_post_PT_calc=None,
                     modify_degenerate_perturbations=False,
                     gaussian_resonance_handling=False,
                     ignore_odd_order_energies=False,
                     intermediate_normalization=False,
                     zero_element_warning=None,
                     degenerate_states=None,
                     zero_order_energy_corrections=None,
                     ):
            """
            :param order: the order of perturbation theory to apply
            :type order: int
            :param expansion_order: the order to go to in the expansions of the perturbations
            :type expansion_order: int
            :param degenerate_states: the set of degeneracies to handle
            :type degenerate_states: Iterable[BasisStateSpace]
            :param coupled_states: explicit bases of states to use at each order in the perturbation theory
            :type coupled_states: Iterable[SelectionRuleStateSpace]
            :param total_space: the total state spaces at each order in the perturbation theory
            :type total_space: Iterable[BasisStateSpace]
            :param flat_total_space: the union of all of the total state spaces
            :type flat_total_space: BasisStateSpace
            :param state_space_iterations: the order to go to when getting the `coupled_states`
            :type state_space_iterations: int
            :param state_space_terms: the explicit set of terms to include, as a tuple `(i, j)` which indicates `(H(i), |n(j)>)`
            :type state_space_terms: Iterable[(int, int)]
            :param state_space_filters: filters that can be used to cut down on the size of bases (see `VPTRunner.get_state_space_filter`)
            :type state_space_filters: dict
            :param allow_post_PT_calc: whether to do the post-perturbation theory variational calculation for degeneracy handling
            :type allow_post_PT_calc: bool
            :param modify_degenerate_perturbations: whether to modify the perturbation representation matrices themselves when doing degeneracy handling
            :type modify_degenerate_perturbations: bool
            :param gaussian_resonance_handling: whether or not to skip the post-PT variational calculation for states with more than two quanta of excitation
            :type gaussian_resonance_handling: bool
            :param ignore_odd_order_energies: whether or not to skip actually calculating the energy corrections for odd orders
            :type ignore_odd_order_energies: bool
            :param intermediate_normalization: whether or not to use 'intermediate normalization' in the wavefunctions
            :type intermediate_normalization: bool
            :param zero_element_warning: whether or not to warn if an element of the representations evaluated to zero (i.e. we wasted effort)
            :type zero_element_warning: bool
            :param zero_order_energy_corrections: energies to use for the zero-order states instead of the diagonal of `H(0)`
            :type zero_order_energy_corrections: dict
            """
            all_opts = dict(
                order=order,
                expansion_order=expansion_order,
                coupled_states=coupled_states,
                total_space=total_space,
                flat_total_space=flat_total_space,
                degenerate_states=degenerate_states,
                state_space_iterations=state_space_iterations,
                state_space_terms=state_space_terms,
                state_space_filters=state_space_filters,
                allow_post_PT_calc=allow_post_PT_calc,
                modify_degenerate_perturbations=modify_degenerate_perturbations,
                gaussian_resonance_handling=gaussian_resonance_handling,
                ignore_odd_order_energies=ignore_odd_order_energies,
                intermediate_normalization=intermediate_normalization,
                zero_element_warning=zero_element_warning,
                zero_order_energy_corrections=zero_order_energy_corrections
            )

            real_opts = {}
            for o, v in all_opts.items():
                if v is not None:
                    real_opts[o] = v

            self.opts = real_opts

    def __init__(self,
                 system,
                 hamiltonian_options=None,
                 perturbation_theory_options=None,
                 runtime_options=None
                 ):
        """
        :param system: the system to run perturbation theory on
        :type system: InputSystem
        :param hamiltonian_options: options to configure the Hamiltonian
        :type hamiltonian_options: HamiltonianOptions
        :param perturbation_theory_options: options to configure the way the perturbation theory is applied
        :type perturbation_theory_options: PerturbationTheoryOptions
        :param runtime_options: options to configure the way the code runs
        :type runtime_options: RuntimeOptions
        """

        if not isinstance(system, self.InputSystem):
            raise NotImplementedError("{} needs to be fed a {} as its `{}` argument".format(
                type(self).__name__,
                self.InputSystem.__name__,
                "system"
            ))

        self.system = system

        self._ham = None
        self._wfns = None
        self.ham_opts = self.HamiltonianOptions() if hamiltonian_options is None else hamiltonian_options
        self.runtime_opts = self.RuntimeOptions() if runtime_options is None else runtime_options
        self.pt_opts = self.PerturbationTheoryOptions() if perturbation_theory_options is None else perturbation_theory_options

    @classmethod
    def get_states(cls, n_quanta, n_modes, target_modes=None, only_target_modes=False):
        """
        Gets states up to `n_quanta` over `n_modes`

        :param n_quanta: the number of quanta to provide excitations for
        :type n_quanta: int | Iterable[int]
        :param n_modes: the number of modes in the system
        :type n_modes: int
        :param target_modes: modes that must be excited
        :type target_modes: Iterable[int]
        :param only_target_modes: whether or not to _only_ support excitations in the `target_modes`
        :type only_target_modes: bool
        :return:
        :rtype:
        """
        if isinstance(n_quanta, int):
            n_quanta = range(n_quanta+1)
        whee = [np.flip(x) for x in BasisStateSpace.from_quanta(
            HarmonicOscillatorProductBasis(n_modes),
            n_quanta
        ).excitations]
        if target_modes is not None:
            whee = [
                p for p in whee if sum(p) == 0 or any(p[i] > 0 for i in target_modes)
            ]
            if only_target_modes:
                whee = [p for p in whee if all(j in target_modes or x == 0 for j,x in enumerate(p))]
        return whee

    @classmethod
    def get_degenerate_polyad_space(cls, states, polyadic_pairs, max_quanta=None):
        """
        Gets degenerate spaces by using pairs of transformation rules to
        take an input state and connect it to other degenerate states

        :param states: the input states
        :type states:
        :param polyadic_pairs: the transformation rules
        :type polyadic_pairs:
        :param max_quanta: the max quanta to allow in connected states
        :type max_quanta:
        :return:
        :rtype:
        """
        # we build a graph of connected states by the polyadic rules
        polyadic_pairs = np.array(polyadic_pairs)
        states = np.array(states)
        poss_degs = [[] for _ in states]
        check_list = states.tolist()
        for n, s in enumerate(states):  # build list-of-lists structure
            for i, nt_spec in enumerate(polyadic_pairs):
                if np.all(s - nt_spec[0] >= 0):
                    new = (s - nt_spec[0] + nt_spec[1]).tolist()
                    if max_quanta is not None and sum(new) > max_quanta:
                        continue
                    if new not in check_list:
                        check_list.append(new)
                        poss_degs.append([])
                    # else:
                    #     poss_degs[idx].append(slist)
                    poss_degs[n].append(new)

        # from the populated lists build the real connection graph
        groups = [[] for _ in check_list]
        new_checks = []
        for i, s1 in enumerate(check_list):
            if s1 not in new_checks:
                new_checks.append(s1)
                groups[i].append(s1)
                groups[i].extend(poss_degs[i])
                for s2, p in zip(check_list[i + 1:], poss_degs[i + 1:]):
                    if s2 not in new_checks:
                        if s2 in poss_degs[i]:
                            new_checks.append(s2)
                            groups[i].extend(p)

        return [g for g in groups if len(g) > 1]

    @classmethod
    def get_state_space_filter(cls, states, n_modes=None, target='wavefunctions'):
        """
        Gets `state_space_filters` for the input `states` targeting some property

        :param states: the input states
        :type states:
        :param n_modes:
        :type n_modes: int
        :param target: the property to target, one of `('frequencies', 'intensities', 'wavefunctions')`
        :type target: str
        :return:
        :rtype:
        """

        if n_modes is None:
            if isinstance(states, BasisStateSpace):
                n_modes = states.ndim
            else:
                n_modes = len(states[0])

        if target == 'wavefunctions':
            return None
        elif target == 'intensities':
            return {
                    (1, 1): (
                        cls.get_states(1, n_modes),
                        (
                            cls.get_states([2, 3], n_modes),
                            [
                                x for x in HarmonicOscillatorProductBasis(n_modes).selection_rules("x", "x", "x")
                                if sum(x) in [-1, 1]
                            ]
                        )
                    ) if any(sum(s) == 3 for s in states) else cls.get_states(2, n_modes),
                    (2, 0): (
                        cls.get_states(1, n_modes),
                        (
                            cls.get_states([3], n_modes),
                            [
                                x for x in HarmonicOscillatorProductBasis(n_modes).selection_rules("x", "x", "x", "x")
                                if sum(x) in [-2, 0]#, 2]
                            ]
                        ),
                        (None, [[]])  # selection rules to apply to remainder
                    ) if any(sum(s) == 3 for s in states) else (
                        cls.get_states(1, n_modes),
                        (None, [[]])  # selection rules to apply to remainder
                    )
            }
        elif target == 'frequencies':
            return {
                (1, 1): ([],),
                (2, 0): (None, [[]])
            }

    def get_Hamiltonian(self):
        return PerturbationTheoryHamiltonian(
            self.system.mol,
            **self.ham_opts.opts,
            **self.runtime_opts.ham_opts
        )

    def get_wavefunctions(self):
        return self.get_Hamiltonian().get_wavefunctions(
            **self.pt_opts.opts,
            **self.runtime_opts.solver_opts
        )

    def print_tables(self, wfns=None):
        """
        Prints a bunch of formatted output data from a PT run

        :param wfns:
        :type wfns:
        :return:
        :rtype:
        """

        if wfns is None:
            wfns = self.get_wavefunctions()

        print("Energy Corrections:")
        print(wfns.format_energy_corrections_table())
        if wfns.degenerate_transformation is not None:
            print("Deperturbed Energies:")
            print(wfns.format_deperturbed_energies_table())
            print("Degenerate Energies:")
            print(wfns.format_energies_table())
        else:
            print("States Energies:")
            print(wfns.format_energies_table())

        ints = wfns.intensities #
        if wfns.degenerate_transformation is not None:
            print("Deperturbed IR Data:")
            for a, m in zip(["X", "Y", "Z"], wfns.format_deperturbed_dipole_contribs_tables()):
                print("{} Dipole Contributions".format(a))
                print(m)
            print(wfns.format_deperturbed_intensities_table())
            print("Degenerate IR Data:")
            for a, m in zip(["X", "Y", "Z"], wfns.format_dipole_contribs_tables()):
                print("{} Dipole Contributions".format(a))
                print(m)
            print(wfns.format_intensities_table())
        else:
            print("IR Data:")
            for a, m in zip(["X", "Y", "Z"], wfns.format_dipole_contribs_tables()):
                print("{} Dipole Contributions".format(a))
                print(m)
            print(wfns.format_intensities_table())