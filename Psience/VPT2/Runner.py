"""
A little package of utilities for setting up/running VPT jobs
"""

import numpy as np

from ..BasisReps import BasisStateSpace, HarmonicOscillatorProductBasis
from ..Molecools import Molecule

from .Hamiltonian import PerturbationTheoryHamiltonian

__all__ = [
    "VPTRunner",
    "InputSystem",
    "HamiltonianOptions",
    "RuntimeOptions",
    "PerturbationTheoryOptions"
]

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
        :param potential_derivatives:
        :type potential_derivatives:
        :param coriolis_coupling:
        :type coriolis_coupling:
        :param include_pseudopotential:
        :type include_pseudopotential:
        :param potential_terms:
        :type potential_terms:
        :param kinetic_terms:
        :type kinetic_terms:
        :param coriolis_terms:
        :type coriolis_terms:
        :param pseudopotential_terms:
        :type pseudopotential_terms:
        :param undimensionalize_normal_modes:
        :type undimensionalize_normal_modes:
        :param use_numerical_jacobians:
        :type use_numerical_jacobians:
        :param eckart_embed_derivatives:
        :type eckart_embed_derivatives:
        :param strip_dummy_atoms:
        :type strip_dummy_atoms:
        :param strip_embedding_coordinates:
        :type strip_embedding_coordinates:
        :param mixed_derivative_handling_mode:
        :type mixed_derivative_handling_mode:
        :param backpropagate_internals:
        :type backpropagate_internals:
        :param zero_mass_term:
        :type zero_mass_term:
        :param internal_fd_mesh_spacing:
        :type internal_fd_mesh_spacing:
        :param internal_fd_stencil:
        :type internal_fd_stencil:
        :param cartesian_fd_mesh_spacing:
        :type cartesian_fd_mesh_spacing:
        :param cartesian_fd_stencil:
        :type cartesian_fd_stencil:
        :param cartesian_analytic_deriv_order:
        :type cartesian_analytic_deriv_order:
        :param internal_by_cartesian_order:
        :type internal_by_cartesian_order:
        :param cartesian_by_internal_order:
        :type cartesian_by_internal_order:
        :param jacobian_warning_threshold:
        :type jacobian_warning_threshold:
        :param check_input_force_constants:
        :type check_input_force_constants:
        :param hessian_tolerance:
        :type hessian_tolerance:
        :param grad_tolerance:
        :type grad_tolerance:
        :param freq_tolerance:
        :type freq_tolerance:
        :param g_derivative_threshold:
        :type g_derivative_threshold:
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
        for o,v in all_opts.items():
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
                 total_space=None,
                 flat_total_space=None,
                 state_space_iterations=None,
                 state_space_terms=None,
                 state_space_filters=None,
                 target_property_rules=None,
                 allow_sakurai_degs=None,
                 allow_post_PT_calc=None,
                 modify_degenerate_perturbations=False,
                 gaussian_resonance_handling=False,
                 ignore_odd_order_energies=False,
                 intermediate_normalization=False,
                 zero_element_warning=None,
                 degenerate_states=None,
                 zero_order_energy_corrections=None,
                 ):
        all_opts = dict(
            order=order,
            coupled_states=coupled_states,
            total_space=total_space,
            flat_total_space=flat_total_space,
            degenerate_states=degenerate_states,
            state_space_iterations=state_space_iterations,
            state_space_terms=state_space_terms,
            state_space_filters=state_space_filters,
            target_property_rules=target_property_rules,
            allow_sakurai_degs=allow_sakurai_degs,
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

class VPTRunner:
    """
    A helper class to make it easier to run jobs
    """

    def __init__(self,
                 system,
                 hamiltonian_options=None,
                 runtime_options=None,
                 perturbation_theory_options=None
                 ):

        if not isinstance(system, InputSystem):
            raise NotImplementedError("{} needs to be fed a {} as its `{}` argument".format(
                type(self).__name__,
                InputSystem.__name__,
                "system"
            ))

        self.system = system

        self._ham = None
        self._wfns = None
        self.ham_opts = HamiltonianOptions() if hamiltonian_options is None else hamiltonian_options
        self.runtime_opts = RuntimeOptions() if runtime_options is None else runtime_options
        self.pt_opts = PerturbationTheoryOptions() if perturbation_theory_options is None else perturbation_theory_options

    @classmethod
    def get_states(cls, n_quanta, n_modes, target_modes=None, only_target_modes=False):
        if isinstance(n_quanta, int):
            n_quanta = range(n_quanta)
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
    def get_state_space_filter(cls, n_modes, target='wavefunctions'):
        if target == 'wavefunctions':
            return None
        elif target == 'intensities':
            return {
                    (1, 1): (
                        cls.get_states(3, n_modes),
                        (
                            cls.get_states(4, n_modes)[len(cls.get_states(3, n_modes)):],
                            [
                                x for x in HarmonicOscillatorProductBasis(n_modes).selection_rules("x", "x", "x")
                                if sum(x) in [-1, 1]
                            ]
                        )
                    ),
                    (2, 0): (
                        cls.get_states(2, n_modes),
                        (
                            cls.get_states(4, n_modes)[len(cls.get_states(3, n_modes)):],
                            [
                                x for x in HarmonicOscillatorProductBasis(n_modes).selection_rules("x", "x", "x", "x")
                                if sum(x) in [-2, 0, 2]
                            ]
                        ),
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
        if wfns is None:
            wfns = self.get_wavefunctions()

        print("Energy Corrections:")
        print(wfns.format_energy_corrections_table())

        print("ZPE and Frequencies:")
        print(wfns.format_energy_corrections_table())

        print("Energy Corrections:")
        print(wfns.format_energy_corrections_table())

        print("ZPE and Frequencies:")
        print(wfns.format_energy_corrections_table())

        for axis, mat in zip(["X", "Y", "Z"], wfns.format_dipole_contribs_tables()):
            print("{} Dipole Contributions".format(axis))
            print(mat)

        print(wfns.format_intensities_table())