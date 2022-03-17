"""
A little package of utilities for setting up/running VPT jobs
"""

import numpy as np, sys

from McUtils.Scaffolding import ParameterManager, Checkpointer
from McUtils.Zachary import FiniteDifferenceDerivative

from ..BasisReps import BasisStateSpace, HarmonicOscillatorProductBasis
from ..Molecools import Molecule

from .DegeneracySpecs import DegeneracySpec
from .Hamiltonian import PerturbationTheoryHamiltonian
from .StateFilters import PerturbationTheoryStateSpaceFilter

__all__ = [
    "VPTRunner",
    "VPTSystem",
    "VPTStateSpace",
    "VPTStateMaker",
    "VPTHamiltonianOptions",
    "VPTRuntimeOptions",
    "VPTSolverOptions"
]

class VPTSystem:
    """
    Provides a little helper for setting up the input
    system for a VPT job
    """

    __props__ = (
        "internals",
        "modes",
        "mode_selection",
        "potential_derivatives",
        "potential_function",
        "order",
        "dipole_derivatives",
        "dummy_atoms",
        "eckart_embed"
    )
    def __init__(self,
                 mol,
                 internals=None,
                 dummy_atoms=None,
                 modes=None,
                 mode_selection=None,
                 potential_derivatives=None,
                 potential_function=None,
                 order=2,
                 dipole_derivatives=None,
                 eckart_embed=False
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
            mol = Molecule.from_spec(mol)
        if internals is not None:
            if dummy_atoms is not None:
                dummy_specs = []
                for where, coord_spec in dummy_atoms:
                    if any(isinstance(x, (float, np.floating)) for x in coord_spec):
                        pos = coord_spec
                    else:
                        ref_1, ref_2, ref_3 = coord_spec
                        a = mol.coords[ref_1] - mol.coords[ref_2]
                        b = mol.coords[ref_3] - mol.coords[ref_2]
                        c = np.cross(a, b)
                        pos = mol.coords[ref_2] + 2 * c / np.linalg.norm(c)
                    dummy_specs.append([where, pos])
                mol = mol.insert_atoms(["X"]*len(dummy_atoms), [d[1] for d in dummy_specs], [d[0] for d in dummy_specs])
            mol.zmatrix = internals
        self.mol = mol
        if modes is not None:
            self.mol.normal_modes.modes = modes
        if potential_derivatives is not None:
            self.mol.potential_derivatives = potential_derivatives
        elif potential_function is not None:
            self.get_potential_derivatives(potential_function, order=order)
        if dipole_derivatives is not None:
            self.mol.dipole_derivatives = dipole_derivatives

        if mode_selection is not None:
            self.mol.normal_modes.modes = self.mol.normal_modes.modes[mode_selection]

        if eckart_embed:
            self.mol = self.mol.get_embedded_molecule()

    @property
    def nmodes(self):
        """
        Provides the number of modes in the system

        :return:
        :rtype:
        """
        return len(self.mol.normal_modes.modes.freqs)

    def get_potential_derivatives(self, potential_function, order=2, **fd_opts):
        """
        Computes potential derivatives for the given function through finite difference

        :param potential_function:
        :type potential_function:
        :param order:
        :type order:
        :param fd_opts:
        :type fd_opts:
        :return:
        :rtype:
        """
        deriv_gen = FiniteDifferenceDerivative(potential_function,
                                               function_shape=((None, None), 0),
                                               stencil=5 + order,
                                               mesh_spacing=1e-3,
                                               **fd_opts
                                               ).derivatives(self.mol.coords)
        self.mol.potential_derivatives = deriv_gen.derivative_tensor(list(range(1, order+3)))

    @classmethod
    def from_harmonic_scan(cls, scan_array):
        raise NotImplementedError("loading from a Harmonic scan needs to be implemented")

class VPTStateSpace:
    """
    Provides a helper to make it easier to set up the input
    state spaces/degenerate spaces to run the perturbation theory
    """

    __props__ = (
        "degeneracy_specs",
    )
    def __init__(self,
                 states,
                 degeneracy_specs=None
                 ):
        self.state_list = states
        self.degenerate_states = self.build_degenerate_state_spaces(degeneracy_specs)
        if self.degenerate_states is not None:
            self.degenerate_states = [np.array(x).tolist() for x in self.degenerate_states]
            states = np.array(self.state_list).tolist()
            for pair in self.degenerate_states:
                for p in pair:
                    if p not in states:
                        states.append(p)
            self.state_list = states

    @classmethod
    def from_system_and_quanta(cls, system, quanta, target_modes=None, only_target_modes=False, **opts):
        """
        Takes a system and a number of quanta and constructs a state space
        based on that

        :param system:
        :type system:
        :param quanta:
        :type quanta:
        :param opts: any of the options supported by
        :type opts:
        :return:
        :rtype:
        """
        n_modes = system.nmodes
        states = cls.get_state_list_from_quanta(
            quanta,
            n_modes,
            target_modes=target_modes,
            only_target_modes=only_target_modes
        )

        return cls(
            states,
            **opts
        )

    @classmethod
    def get_state_list_from_quanta(cls, n_quanta, n_modes, target_modes=None, only_target_modes=False):
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

    def build_degenerate_state_spaces(self, degeneracy_specs):
        """
        :param degeneracy_specs:
        :type degeneracy_specs:
        :return:
        :rtype:
        """

        spec = DegeneracySpec.from_spec(degeneracy_specs)
        if spec is None:
            return None
        else:
            return spec.get_groups(self.state_list)

        # elif isinstance(degeneracy_specs, dict):
        #     # dispatch on mode
        #     degeneracy_specs = degeneracy_specs.copy()
        #     if 'polyads' in degeneracy_specs:
        #         polyads = degeneracy_specs['polyads']
        #         del degeneracy_specs['polyads']
        #         return DegenerateMultiStateSpace.get_degenerate_polyad_space(
        #             self.state_list,
        #             polyads,
        #             **degeneracy_specs
        #         )
        #     else:
        #         NotImplementedError(
        #             "couldn't infer degenerate space construction mode from spec {}".format(degeneracy_specs)
        #         )
        # elif all(DegenerateMultiStateSpace._is_polyad_rule(d, n_modes) for d in degeneracy_specs):
        #     return DegenerateMultiStateSpace.get_degenerate_polyad_space(
        #         self.state_list,
        #         degeneracy_specs
        #     )
        # else:
        #     raise NotImplementedError("don't know what to do with degeneracy spec {}".format(degeneracy_specs))

    def get_filter(self, target_property, order=2):
        """
        Obtains a state space filter for the given target property
        using the states we want to get corrections for

        :param target_property:
        :type target_property:
        :param order:
        :type order:
        :return:
        :rtype:
        """
        return self.get_state_space_filter(self.state_list,
                                           target=target_property,
                                           order=order
                                           )
    @classmethod
    def get_state_space_filter(cls, states, n_modes=None, order=2, target='wavefunctions'):
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

        if order != 2:
            raise ValueError("state space filters currently only implemented at second order")

        if target == 'wavefunctions':
            return None
        elif target == 'intensities':
            return PerturbationTheoryStateSpaceFilter.from_property_rules(
                cls.get_state_list_from_quanta(0, n_modes),
                states,
                [
                    HarmonicOscillatorProductBasis(n_modes).selection_rules("x", "x", "x"),
                    HarmonicOscillatorProductBasis(n_modes).selection_rules("x", "x", "x", "x")
                ],
                [
                    # (),
                    # (),
                    # ()
                    HarmonicOscillatorProductBasis(n_modes).selection_rules("x"),
                    HarmonicOscillatorProductBasis(n_modes).selection_rules("x", "x"),
                    HarmonicOscillatorProductBasis(n_modes).selection_rules("x", "x", "x")
                ]
            )

            return {
                    (1, 1): (
                        cls.get_state_list_from_quanta(1, n_modes),
                        (
                            cls.get_state_list_from_quanta([2, 3], n_modes),
                            [
                                x for x in HarmonicOscillatorProductBasis(n_modes).selection_rules("x", "x", "x")
                                if sum(x) in [-1, 1]
                            ]
                        )
                    ) if any(sum(s) == 3 for s in states) else (
                        cls.get_state_list_from_quanta(1, n_modes),
                        (
                            cls.get_state_list_from_quanta([2], n_modes),
                            [
                                x for x in HarmonicOscillatorProductBasis(n_modes).selection_rules("x", "x", "x")
                                if sum(x) in [-1, 1]
                            ]
                        )
                    ),
                    (2, 0): (
                        cls.get_state_list_from_quanta(1, n_modes),
                        (
                            cls.get_state_list_from_quanta([3], n_modes),
                            [
                                x for x in HarmonicOscillatorProductBasis(n_modes).selection_rules("x", "x", "x", "x")
                                if sum(x) in [-2, 0]#, 2]
                            ]
                        ),
                        (None, [[]])  # selection rules to apply to remainder
                    ) if any(sum(s) == 3 for s in states) else (
                        cls.get_state_list_from_quanta(1, n_modes),
                        (None, [[]])  # selection rules to apply to remainder
                    )
            }
        elif target == 'frequencies':
            return {
                (1, 1): ([],),
                (2, 0): (None, [[]])
            }

class VPTHamiltonianOptions:
    """
    Provides a helper to keep track of the levers available for
    setting up the Hamiltonian
    """

    __props__ = (
         "include_coriolis_coupling",
         "include_pseudopotential",
         "potential_terms",
         "kinetic_terms",
         "coriolis_terms",
         "pseudopotential_terms",
         "undimensionalize_normal_modes",
         "use_numerical_jacobians",
         "eckart_embed_derivatives",
         "eckart_embed_planar_ref_tolerance",
         "strip_dummy_atoms",
         "strip_embedding_coordinates",
         "mixed_derivative_handling_mode",
         "backpropagate_internals",
         "direct_propagate_cartesians",
         "zero_mass_term",
         "internal_fd_mesh_spacing",
         "internal_fd_stencil",
         "cartesian_fd_mesh_spacing",
         "cartesian_fd_stencil",
         "cartesian_analytic_deriv_order",
         "internal_by_cartesian_order",
         "cartesian_by_internal_order",
         "jacobian_warning_threshold",
         "check_input_force_constants",
         "hessian_tolerance",
         "grad_tolerance",
         "freq_tolerance",
         "g_derivative_threshold"
    )

    def __init__(self,
                 include_coriolis_coupling=None,
                 include_pseudopotential=None,
                 potential_terms=None,
                 kinetic_terms=None,
                 coriolis_terms=None,
                 pseudopotential_terms=None,
                 undimensionalize_normal_modes=None,
                 use_numerical_jacobians=None,
                 eckart_embed_derivatives=None,
                 eckart_embed_planar_ref_tolerance=None,
                 strip_dummy_atoms=None,
                 strip_embedding_coordinates=None,
                 mixed_derivative_handling_mode=None,
                 backpropagate_internals=None,
                 direct_propagate_cartesians=None,
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
        :param include_coriolis_coupling: whether or not to include Coriolis coupling in Cartesian normal mode calculation
        :type include_coriolis_coupling: bool
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
            include_coriolis_coupling=include_coriolis_coupling,
            include_pseudopotential=include_pseudopotential,
            potential_terms=potential_terms,
            kinetic_terms=kinetic_terms,
            coriolis_terms=coriolis_terms,
            pseudopotential_terms=pseudopotential_terms,
            undimensionalize=undimensionalize_normal_modes,
            numerical_jacobians=use_numerical_jacobians,
            eckart_embed_derivatives=eckart_embed_derivatives,
            eckart_embed_planar_ref_tolerance=eckart_embed_planar_ref_tolerance,
            strip_dummies=strip_dummy_atoms,
            strip_embedding=strip_embedding_coordinates,
            mixed_derivative_handling_mode=mixed_derivative_handling_mode,
            backpropagate_internals=backpropagate_internals,
            direct_propagate_cartesians=direct_propagate_cartesians,
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

class VPTRuntimeOptions:
    """
    Provides a helper to keep track of the options available
    for configuring the way the code runs
    """

    __props__ = (
        "operator_chunk_size",
        "logger",
        "verbose",
        "checkpoint",
        "results",
        "memory_constrained",
        "checkpoint_keys",
        "use_cached_representations",
        "use_cached_basis"
    )
    def __init__(self,
                 operator_chunk_size=None,
                 logger=None,
                 verbose=None,
                 checkpoint=None,
                 results=None,
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
            results=results,
            use_cached_representations=use_cached_representations,
            use_cached_basis=use_cached_basis
        )
        real_solver_run_opts = {}
        for o, v in solver_run_opts.items():
            if v is not None:
                real_solver_run_opts[o] = v
        self.solver_opts = real_solver_run_opts

class VPTSolverOptions:
    """
    Provides a helper to keep track of the options available
    for configuring the way the perturbation theory is applied
    """

    __props__ = (
        "order",
        "expansion_order",
        "coupled_states",
        "total_space",
        "flat_total_space",
        "state_space_iterations",
        "state_space_terms",
        "state_space_filters",
        "allow_post_PT_calc",
        "modify_degenerate_perturbations",
        "gaussian_resonance_handling",
        "ignore_odd_order_energies",
        "intermediate_normalization",
        "zero_element_warning",
        "degenerate_states",
        "zero_order_energy_corrections",
        "handle_strong_couplings",
        "strong_coupling_test_modes",
        "strong_couplings_state_filter",
        "strongly_coupled_group_filter"
    )
    def __init__(self,
                 order=2,
                 expansion_order=None,
                 coupled_states=None,
                 total_space=None,
                 flat_total_space=None,
                 state_space_iterations=None,
                 state_space_terms=None,
                 state_space_filters=None,
                 allow_post_PT_calc=None,
                 modify_degenerate_perturbations=None,
                 gaussian_resonance_handling=None,
                 ignore_odd_order_energies=None,
                 intermediate_normalization=None,
                 zero_element_warning=None,
                 degenerate_states=None,
                 handle_strong_couplings=None,
                 strong_coupling_test_modes=None,
                 strong_couplings_state_filter=None,
                 strongly_coupled_group_filter=None,
                 zero_order_energy_corrections=None

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
            handle_strong_couplings=handle_strong_couplings,
            strong_coupling_test_modes=strong_coupling_test_modes,
            strong_couplings_state_filter=strong_couplings_state_filter,
            strongly_coupled_group_filter=strongly_coupled_group_filter,
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

    @staticmethod
    def get_zero_order_energies(corrected_fundamental_freqs, states):
        """

        :param corrected_fundamental_freqs:
        :type corrected_fundamental_freqs:
        :param states:
        :type states:
        :return:
        :rtype:
        """
        corrected_fundamental_freqs = np.asanyarray(corrected_fundamental_freqs)

        return [
            (s, np.dot(np.array(s) + 1 / 2, corrected_fundamental_freqs))
            for s in states
        ]

class VPTRunner:
    """
    A helper class to make it easier to run jobs by making the inputs/options
    clear and making it easier to customize run options
    """

    def __init__(self,
                 system,
                 states,
                 hamiltonian_options=None,
                 solver_options=None,
                 runtime_options=None
                 ):
        """
        :param system: the system to run perturbation theory on
        :type system: VPTSystem
        :param hamiltonian_options: options to configure the Hamiltonian
        :type hamiltonian_options: VPTHamiltonianOptions
        :param solver_options: options to configure the way the perturbation theory is applied
        :type solver_options: VPTSolverOptions
        :param runtime_options: options to configure the way the code runs
        :type runtime_options: VPTRuntimeOptions
        """

        if not isinstance(system, VPTSystem):
            raise NotImplementedError("{} needs to be fed a {} as its `{}` argument".format(
                type(self).__name__,
                VPTSystem.__name__,
                "system"
            ))

        if not isinstance(states, VPTStateSpace):
            raise NotImplementedError("{} needs to be fed a {} as its `{}` argument".format(
                type(self).__name__,
                VPTStateSpace.__name__,
                "states"
            ))

        self.system = system
        self.states = states

        self._ham = None
        self._wfns = None
        self.ham_opts = VPTHamiltonianOptions() if hamiltonian_options is None else hamiltonian_options
        self.runtime_opts = VPTRuntimeOptions() if runtime_options is None else runtime_options
        self.pt_opts = VPTSolverOptions() if solver_options is None else solver_options

    def get_Hamiltonian(self):
        return PerturbationTheoryHamiltonian(
            self.system.mol,
            **self.ham_opts.opts,
            **self.runtime_opts.ham_opts
        )

    @property
    def hamiltonian(self):
        if self._ham is None:
            self._ham = self.get_Hamiltonian()
        return self._ham

    def get_wavefunctions(self):
        pt_opts = self.pt_opts.opts.copy()
        if 'degenerate_states' not in pt_opts:
            pt_opts['degenerate_states'] = self.states.degenerate_states
        return self.hamiltonian.get_wavefunctions(
            self.states.state_list,
            **pt_opts,
            **self.runtime_opts.solver_opts
        )

    @classmethod
    def print_output_tables(cls, wfns=None, file=None, print_intensities=True, sep_char="=", sep_len=100):
        """
        Prints a bunch of formatted output data from a PT run

        :param wfns:
        :type wfns:
        :return:
        :rtype:
        """

        if wfns.logger is not None:
            def print_block(label, *args, **kwargs):
                with wfns.logger.block(tag=label):
                    wfns.logger.log_print(" ".join("{}".format(x) for x in args), **kwargs)
        else:
            if file is None:
                file = sys.stdout

            def print_label(label, file=file, **opts):
                lablen = len(label) + 2
                split_l = int(np.floor((sep_len - lablen) / 2))
                split_r = int(np.ceil((sep_len - lablen) / 2))
                print(sep_char * split_l, label, sep_char * split_r, **opts, file=file)

            def print_footer(label=None, file=file, **opts):
                print(sep_char * sep_len, **opts, file=file)

            def print_block(label, *args, file=file, **kwargs):
                print_label(label, file=file, **kwargs)
                print(*args, file=file, **kwargs)
                print_footer(file=file, **kwargs)

        print_block("Energy Corrections", wfns.format_energy_corrections_table())
        if wfns.degenerate_transformation is not None:
            print_block("Deperturbed Energies",
                        wfns.format_deperturbed_energies_table()
                        )
            print_block(
                "Degenerate Energies",
                wfns.format_energies_table()
            )
        else:
            print_block("States Energies",
                wfns.format_energies_table()
            )

        if print_intensities:
            ints = wfns.intensities # to make sure they're computed before printing starts
            if wfns.degenerate_transformation is not None:
                for a, m in zip(["X", "Y", "Z"], wfns.format_deperturbed_dipole_contribs_tables()):
                    print_block("{} Deperturbed Dipole Contributions".format(a), m)

                print_block("Deperturbed IR Data",
                            wfns.format_deperturbed_intensities_table()
                            )

            for a, m in zip(["X", "Y", "Z"], wfns.format_dipole_contribs_tables()):
                print_block("{} Dipole Contributions".format(a), m)
            print_block("IR Data", wfns.format_intensities_table())

    def print_tables(self, wfns=None, file=None, print_intensities=True, sep_char="=", sep_len=100):
        """
        Prints a bunch of formatted output data from a PT run

        :param wfns:
        :type wfns:
        :return:
        :rtype:
        """

        if wfns is None:
            wfns = self.get_wavefunctions()

        self.print_output_tables(wfns=wfns, file=file,
                                 print_intensities=print_intensities, sep_char=sep_char, sep_len=sep_len)

        return wfns

    @classmethod
    def run_simple(cls,
                   system,
                   states,
                   target_property=None,
                   corrected_fundamental_frequencies=None,
                   **opts
                   ):

        full_opts = (
                VPTSystem.__props__
                + VPTStateSpace.__props__
                + VPTSolverOptions.__props__
                + VPTHamiltonianOptions.__props__
                + VPTRuntimeOptions.__props__
                + VPTSolverOptions.__props__
        )

        misses = set(opts.keys()).difference(set(full_opts))
        if len(misses) > 0:
            raise ValueError("{}: options {} not valid, full listing is {}".format(
                cls.__name__,
                misses,
                "\n  ".join(("",) + full_opts)
            ))

        par = ParameterManager(**opts)
        sys = VPTSystem(system, **par.filter(VPTSystem))
        if isinstance(states, int) or isinstance(states[0], int):
            states = VPTStateSpace.from_system_and_quanta(
                sys,
                states,
                **par.filter(VPTStateSpace)
            )
        else:
            states = VPTStateSpace(
                states,
                **par.filter(VPTStateSpace)
            )

        order = 2 if 'order' not in opts.keys() else opts['order']
        if target_property is None and order == 2:
            target_property = 'intensities'
        if target_property is not None and 'state_space_filters' not in opts:
            par.ops['state_space_filters'] = states.get_filter(target_property, order=order)

            # print(par.ops['state_space_filters'])

        if corrected_fundamental_frequencies is not None and (
            'zero_order_energy_corrections' not in opts
            or opts['zero_order_energy_corrections'] is None
        ):
            par.ops['zero_order_energy_corrections'] = VPTSolverOptions.get_zero_order_energies(
                corrected_fundamental_frequencies,
                states.state_list
            )

        runner = cls(
            sys,
            states,
            hamiltonian_options=VPTHamiltonianOptions(**par.filter(VPTHamiltonianOptions)),
            runtime_options=VPTRuntimeOptions(**par.filter(VPTRuntimeOptions)),
            solver_options=VPTSolverOptions(**par.filter(VPTSolverOptions))
        )

        logger = runner.hamiltonian.logger
        with logger.block(tag="Starting Perturbation Theory Runner"):
            with logger.block(tag="Hamiltonian Options"):
                for k,v in runner.ham_opts.opts.items():
                    logger.log_print("{k}: {v:<100.100}", k=k, v=v, preformatter=lambda *a,k=k,v=v,**kw:dict({'k':k, 'v':str(v)}, **kw))
            with logger.block(tag="Solver Options"):
                opts = runner.pt_opts.opts
                for k,v in opts.items():
                    logger.log_print("{k}: {v:<100.100}", k=k, v=v, preformatter=lambda *a,k=k,v=v,**kw:dict({'k':k, 'v':str(v)}, **kw))
            with logger.block(tag="Runtime Options"):
                opts = dict(runner.runtime_opts.ham_opts, **runner.runtime_opts.solver_opts)
                for k,v in opts.items():
                    logger.log_print("{k}: {v:<100.100}", k=k, v=v, preformatter=lambda *a,k=k,v=v,**kw:dict({'k':k, 'v':str(v)}, **kw))

            with logger.block(tag="States"):
                logger.log_print(
                    states.state_list,
                    message_prepper=lambda a:str(np.array(a)).splitlines()
                )
            with logger.block(tag="Degeneracies"):
                if states.degenerate_states is None:
                    logger.log_print("None")
                else:
                    ds = states.degenerate_states
                    if isinstance(ds, list):
                        for x in states.degenerate_states:
                            logger.log_print(
                                x,
                                message_prepper=lambda a:str(np.array(a)).splitlines()
                            )
                    else:
                        logger.log_print(str(ds))

            return runner.print_tables()

class VPTStateMaker:
    """
    A tiny but useful class to make states based on their quanta
    of excitation
    """

    def __init__(self, ndim, mode='low-high'):
        self.ndim = ndim
        self.mode = mode

    def make_state(self, *specs, mode=None):

        if mode is None:
            mode = self.mode

        state = [0] * self.ndim
        for s in specs:
            if isinstance(s, (int, np.integer)):
                i = s
                q = 1
            else:
                i,q = s
            if mode == 'low-high':
                state[-i] = q
            elif mode == 'high-low':
                state[i-1] = q
            elif mode == 'normal':
                state[i] = q
            else:
                raise ValueError("don't know what to do with filling mode '{}'".format(mode))
        return state

    def __call__(self, *specs, mode=None):
        return self.make_state(*specs, mode=mode)