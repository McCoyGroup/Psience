"""
A little package of utilities for setting up/running VPT jobs
"""

import numpy as np, sys, os, itertools, scipy

from McUtils.Data import UnitsData, AtomData
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

__reload_hook__ = ["..BasisRepss", "..Molecools", ".DegeneracySpecs", ".Hamiltonian", ".StateFilters"]

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

        if eckart_embed is True:
            self.mol = self.mol.get_embedded_molecule()
        else:
            if isinstance(eckart_embed, np.ndarray):
                ref = self.mol.copy()
                ref.coords = eckart_embed
                eckart_embed = ref
            if isinstance(eckart_embed, Molecule):
                self.mol = self.mol.get_embedded_molecule(eckart_embed)


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
                 degeneracy_specs=None,
                 system=None
                 ):
        """
        :param states: A list of states or a number of quanta to target
        :type states: list | int
        :param degeneracy_specs: A specification of degeneracies, either as polyads or explicit groups of states
        :type degeneracy_specs: list | dict
        """
        if not isinstance(states, BasisStateSpace):
            states = BasisStateSpace(
                HarmonicOscillatorProductBasis(len(states[0])),
                states
            )
        basis = states.basis
        self.state_list = states.excitations.tolist()
        self.degenerate_states = self.build_degenerate_state_spaces(degeneracy_specs, states, system=system)
        if self.degenerate_states is not None:
            if isinstance(self.degenerate_states, tuple) and isinstance(self.degenerate_states[0], DegeneracySpec):
                self.degenerate_states, new_states = self.degenerate_states
                if isinstance(new_states, BasisStateSpace):
                    new_states = new_states.excitations.tolist()
                states = np.asanyarray(self.state_list).tolist()
                for state in new_states:
                    if state not in states:
                        states.append(state)
                self.states = BasisStateSpace(basis, states)
                self.state_list = states
            else:
                self.degenerate_states = [np.array(x).tolist() for x in self.degenerate_states]
                states = np.asanyarray(self.state_list).tolist()
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
            system=system,
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

    def build_degenerate_state_spaces(self, degeneracy_specs, states, system=None):
        """
        :param degeneracy_specs:
        :type degeneracy_specs:
        :return:
        :rtype:
        """

        spec = DegeneracySpec.from_spec(degeneracy_specs)
        if hasattr(spec, 'frequencies') and spec.frequencies is None:
            spec.frequencies = system.mol.normal_modes.modes.freqs
        # raise Exception(spec)
        if spec is None:
            return None
        elif hasattr(spec, 'prep_states'):
            return (spec, spec.prep_states(states))
        else:
            return spec.get_groups(states)

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

    def filter_generator(self, target_property, order=2):
        def filter(states):
            return self.get_state_space_filter(states, target=target_property, order=order)
        return filter
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

        # if order != 2:
        #     raise ValueError("state space filters currently only implemented at second order")
        if target == 'wavefunctions':
            return None
        elif target == 'intensities':
            return PerturbationTheoryStateSpaceFilter.from_property_rules(
                cls.get_state_list_from_quanta(0, n_modes),
                states,
                [
                    HarmonicOscillatorProductBasis(n_modes).selection_rules(*["x"]*i) for i in range(3, order+3)
                ],
                [
                    HarmonicOscillatorProductBasis(n_modes).selection_rules(*["x"]*i) for i in range(1, order+2)
                ],
                order=order
            )
        elif target == 'frequencies':
            # return {
            #     (1, 1): ([],),
            #     (2, 0): (None, [[0]])
            # }
            return PerturbationTheoryStateSpaceFilter.from_property_rules(
                cls.get_state_list_from_quanta(0, n_modes),
                states,
                [
                    HarmonicOscillatorProductBasis(n_modes).selection_rules(*["x"] * i) for i in range(3, order + 3)
                ],
                None,
                order=order
                # [
                #     # (),
                #     # (),
                #     # ()
                #     # HarmonicOscillatorProductBasis(n_modes).selection_rules("x"),
                #     # HarmonicOscillatorProductBasis(n_modes).selection_rules("x", "x"),
                #     # HarmonicOscillatorProductBasis(n_modes).selection_rules("x", "x", "x")
                # ]
            )

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
         "dipole_terms",
         "dipole_derivatives",
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
         "use_internal_modes",
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
         "g_derivative_threshold",
         "gmatrix_tolerance",
        'use_cartesian_kinetic_energy'
    )

    def __init__(self,
                 include_coriolis_coupling=None,
                 include_pseudopotential=None,
                 potential_terms=None,
                 kinetic_terms=None,
                 coriolis_terms=None,
                 pseudopotential_terms=None,
                 dipole_terms=None,
                 dipole_derivatives=None,
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
                 g_derivative_threshold=None,
                 gmatrix_tolerance=None,
                 use_internal_modes=None,
                 use_cartesian_kinetic_energy=None
                 ):
        """
        :param include_coriolis_coupling: whether or not to include Coriolis coupling in Cartesian normal mode calculation
        :type include_coriolis_coupling: bool
        :param include_pseudopotential: whether or not to include the pseudopotential/Watson term
        :type include_pseudopotential: bool
        :param potential_terms: explicit values for the potential terms (e.g. from analytic models), should be a list of tensors starting with the Hessian with each axis of length `nmodes`
        :type potential_terms: Iterable[np.ndarray]
        :param kinetic_terms: explicit values for the kinetic terms (e.g. from analytic models), same format as for the potential
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
            dipole_terms=dipole_terms,
            dipole_derivatives=dipole_derivatives,
            undimensionalize=undimensionalize_normal_modes,
            numerical_jacobians=use_numerical_jacobians,
            use_internal_modes=use_internal_modes,
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
            g_derivative_threshold=g_derivative_threshold,
            gmatrix_tolerance=gmatrix_tolerance,
            use_cartesian_kinetic_energy=use_cartesian_kinetic_energy
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
            results=results,
            parallelizer=parallelizer
        )
        real_ham_opts = {}
        for o, v in ham_run_opts.items():
            if v is not None:
                real_ham_opts[o] = v
        self.ham_opts = real_ham_opts

        solver_run_opts = dict(
            # operator_chunk_size=operator_chunk_size,
            memory_constrained=memory_constrained,
            verbose=verbose,
            checkpoint_keys=checkpoint_keys,
            # results=results,
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
        # "strong_coupling_test_modes",
        # "strong_couplings_state_filter",
        # "strongly_coupled_group_filter",
        # "extend_strong_coupling_spaces",
        # "strong_coupling_zero_order_energy_cutoff",
        "low_frequency_mode_cutoff"
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
                 extend_strong_coupling_spaces=None,
                 strong_coupling_zero_order_energy_cutoff=None,
                 low_frequency_mode_cutoff=None,
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
            extend_strong_coupling_spaces=extend_strong_coupling_spaces,
            strong_coupling_zero_order_energy_cutoff=strong_coupling_zero_order_energy_cutoff,
            low_frequency_mode_cutoff=low_frequency_mode_cutoff,
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
    def print_output_tables(cls, wfns=None, file=None,
                            print_intensities=True,
                            print_energies=True,
                            print_energy_corrections=True,
                            print_transition_moments=True,
                            logger=None, sep_char="=", sep_len=100):
        """
        Prints a bunch of formatted output data from a PT run

        :param wfns:
        :type wfns:
        :return:
        :rtype:
        """

        if logger is None:
            logger = wfns.logger
        if logger is not None:
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

        if print_energy_corrections:
            print_block("Energy Corrections", wfns.format_energy_corrections_table())
        if print_energies:
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
            if print_transition_moments:
                if wfns.degenerate_transformation is not None:
                    for a, m in zip(["X", "Y", "Z"], wfns.format_deperturbed_dipole_contribs_tables()):
                        print_block("{} Deperturbed Dipole Contributions".format(a), m)

                    print_block("Deperturbed IR Data",
                                wfns.format_deperturbed_intensities_table()
                                )

                for a, m in zip(["X", "Y", "Z"], wfns.format_dipole_contribs_tables()):
                    print_block("{} Dipole Contributions".format(a), m)
            print_block("IR Data", wfns.format_intensities_table())

    def print_tables(self,
                     wfns=None, file=None,
                     print_intensities=True,
                     print_energy_corrections=True,
                     print_transition_moments=True,
                     sep_char="=", sep_len=100):
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
                                 print_intensities=print_intensities,
                                 print_energy_corrections=print_energy_corrections,
                                 print_transition_moments=print_transition_moments,
                                 sep_char=sep_char, sep_len=sep_len)

        return wfns

    @classmethod
    def construct(cls,
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
                system=sys,
                **par.filter(VPTStateSpace)
            )

        order = 2 if 'order' not in opts.keys() else opts['order']
        if target_property is None and order == 2:
            target_property = 'intensities'
        if target_property is not None and 'state_space_filters' not in opts:
            if 'expansion_order' in opts.keys():
                expansion_order = opts['expansion_order']
            else:
                expansion_order = order
            par.ops['state_space_filters'] = states.filter_generator(target_property, order=order)

            # print(par.ops['state_space_filters'])

        if corrected_fundamental_frequencies is not None and (
                'zero_order_energy_corrections' not in opts
                or opts['zero_order_energy_corrections'] is None
        ):
            par.ops['zero_order_energy_corrections'] = VPTSolverOptions.get_zero_order_energies(
                corrected_fundamental_frequencies,
                states.state_list
            )

        hops = VPTHamiltonianOptions(**par.filter(VPTHamiltonianOptions))
        rops = VPTRuntimeOptions(**par.filter(VPTRuntimeOptions))
        sops = VPTSolverOptions(**par.filter(VPTSolverOptions))

        runner = cls(
            sys,
            states,
            hamiltonian_options=hops,
            runtime_options=rops,
            solver_options=sops
        )
        return runner, (sys, states, hops, rops, sops)
    @classmethod
    def run_simple(cls,
                   system,
                   states,
                   target_property=None,
                   corrected_fundamental_frequencies=None,
                   calculate_intensities=True,
                   **opts
                   ):

        runner, (system, states, hops, rops, sops) = cls.construct(
            system,
            states,
            target_property=target_property,
            corrected_fundamental_frequencies=corrected_fundamental_frequencies,
            **opts
        )

        logger = runner.hamiltonian.logger
        with np.printoptions(linewidth=1e8, threshold=1e6):
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

                return runner.print_tables(print_intensities=calculate_intensities)

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

class AnneInputHelpers:

    @classmethod
    def _check_file(cls, no_file):
        if len(no_file.splitlines()) < 2 and '.' in no_file:
            raise FileNotFoundError("{} is not a file".format(no_file))

    @classmethod
    def get_tensor_idx(cls, line, inds, m, start_at=0):
        bits = line.split()
        idx = tuple(int(i) - 1 for i in bits if '.' not in i)
        m = max(m, max(idx[start_at:]))
        val = float(bits[len(idx)])
        inds[idx] = val
        return m

    @classmethod
    def parse_tensor(cls, block, dims=None):
        inds = {}
        m = 0
        if os.path.isfile(block):
            with open(block) as f:
                for line in f:
                    line = line.strip()
                    if len(line) > 0:
                        m = cls.get_tensor_idx(line, inds, m)
        else:
            cls._check_file(block)
            for line in block.splitlines():
                line = line.strip()
                if len(line) > 0:
                    m = cls.get_tensor_idx(line, inds, m)
        if dims is None:
            m = m + 1
            n = len(next(iter(inds)))
            dims = (m,) * n
        a = np.zeros(dims)
        for i, v in inds.items():
            for p in itertools.permutations(i):
                a[p] = v
        return a
    @classmethod
    def parse_dipole_tensor(cls, block, dims=None):
        inds = {}
        m = 0
        if os.path.isfile(block):
            with open(block) as f:
                for line in f:
                    line = line.strip()
                    if len(line) > 0:
                        m = cls.get_tensor_idx(line, inds, m, start_at=1)
        else:
            cls._check_file(block)
            for line in block.splitlines():
                line = line.strip()
                if len(line) > 0:
                    m = cls.get_tensor_idx(line, inds, m, start_at=1)
        if dims is None:
            m = m + 1
            n = len(next(iter(inds)))
            dims = (3,) + (m,) * (n-1)
        a = np.zeros(dims)
        for i, v in inds.items():
            u = i[0]
            i = i[1:]
            for p in itertools.permutations(i):
                p = (u,) + p
                a[p] = v
        # print(">>", a.shape)
        # a = np.transpose(a, list(range(1, a.ndim)) + [0]) # needs to be innermost I guess...
        # print(">", a.shape)
        return a

    @classmethod
    def parse_freqs_line(cls, line):
        data = [float(x) for x in line.split()]
        if len(data) > 0:
            return np.array(data)
        else:
            return None

    @classmethod
    def parse_modes_line(cls, line, nmodes):
        data = [float(x) for x in line.split()]
        if len(data) > 0:
            l = len(data)
            n = int(l/nmodes)
            return np.array(data).reshape(nmodes, n).T
        else:
            return None

    @classmethod
    def parse_modes(cls, block):
        inds = {}
        m = 0
        freqs = None
        L = None
        Linv = None
        if os.path.isfile(block):
            with open(block) as f:
                for line in f:
                    line = line.strip()
                    if len(line) > 0:
                        if freqs is None:
                            freqs = cls.parse_freqs_line(line)
                        elif L is None:
                            L = cls.parse_modes_line(line, len(freqs))
                        elif Linv is None:
                            Linv = cls.parse_modes_line(line, L.shape[1])
                            break
        else:
            cls._check_file(block)
            for line in block.splitlines():
                line = line.strip()
                if len(line) > 0:
                    if freqs is None:
                        freqs = cls.parse_freqs_line(line)
                    elif L is None:
                        L = cls.parse_modes_line(line, len(freqs))
                    elif Linv is None:
                        Linv = cls.parse_modes_line(line, L.shape[1])
                        break
        return freqs, L, Linv

    @classmethod
    def parse_coords(cls, block):
        coords = []
        if os.path.isfile(block):
            with open(block) as f:
                for line in f:
                    line = line.strip()
                    if len(line) > 0:
                        coords.append([float(x) for x in line.split()])
        else:
            cls._check_file(block)
            for line in block.splitlines():
                line = line.strip()
                if len(line) > 0:
                    coords.append([float(x) for x in line.split()])
        return np.array(coords)

    @classmethod
    def parse_atoms(cls, block):
        coords = []
        if os.path.isfile(block):
            with open(block) as f:
                for line in f:
                    line = line.strip()
                    if len(line) > 0:
                        try:
                            coords.extend([int(x) for x in line.split()])
                        except ValueError:
                            pass
        else:
            cls._check_file(block)
            for line in block.splitlines():
                line = line.strip()
                if len(line) > 0:
                    try:
                        coords.extend([int(x) for x in line.split()])
                    except ValueError:
                        pass
        return [AtomData[x]["Symbol"] for x in coords]

    @classmethod
    def parse_masses(cls, block):
        coords = []
        if os.path.isfile(block):
            with open(block) as f:
                for line in f:
                    line = line.strip()
                    if len(line) > 0:
                        try:
                            coords.extend([float(x) for x in line.split()])
                        except ValueError:
                            pass
        else:
            cls._check_file(block)
            for line in block.splitlines():
                line = line.strip()
                if len(line) > 0:
                    try:
                        coords.extend([float(x) for x in line.split()])
                    except ValueError:
                        pass
        return coords

    @classmethod
    def parse_zmatrix(cls, block):
        _ = 10000
        zmat = [[0, _, _, _]]
        n = 1
        if os.path.isfile(block):
            with open(block) as f:
                for line in f:
                    line = line.strip()
                    if len(line) > 0:
                        split = line.split()
                        if len(split) == 4:
                            split = split[:3]
                        zmat.append([n] + [int(x)-1 for x in split] + [_]*(3-len(split)))
                        n += 1
                    elif len(zmat) > 1:
                        break
        else:
            cls._check_file(block)
            for line in block.splitlines():
                line = line.strip()
                if len(line) > 0:
                    split = line.split()
                    if len(split) == 4:
                        split = split[:3]
                    zmat.append([n] + [int(x) - 1 for x in split] + [_] * (3 - len(split)))
                    n += 1
                elif len(zmat) > 1:
                    break
        return np.array(zmat)

    @classmethod
    def standard_sorting(cls, zmat):
        """
        converts from [r1, r2, r3, ..., a1, a2, ..., t1, t2, ...] coords
        to standard zmat coords
        :param zmat:
        :type zmat:
        :return:
        :rtype:
        """
        nats = len(zmat)
        ncoords = 3*nats - 6
        if nats < 4:
            return None
        else:
            r_coords = [0, 1, 3]
            a_coords = [2, 4]
            t_coords = [5]
            if nats > 4:
                extra = np.arange(6, ncoords+1)
                r_coords += extra[::4].tolist()
                a_coords += extra[1::4].tolist()
                t_coords += extra[2::4].tolist()
            return np.argsort(np.concatenate([r_coords, a_coords, t_coords]))

    @classmethod
    def get_internal_FG(cls, freqs, modes, inv, sorting=None):
        modes = modes.T * np.sqrt(freqs)[:, np.newaxis]
        inv = inv.T / np.sqrt(freqs)[np.newaxis, :]
        G = np.dot(modes.T, modes)
        # print(np.dot(modes, modes.T))
        # print(np.dot(modes.T, modes))
        F = np.dot(np.dot(inv, np.diag(freqs ** 2)), inv.T)
        if sorting is not None:
            G = G[np.ix_(sorting, sorting)]
            # print(G)
            # print(np.round(F*219465))
            F = F[np.ix_(sorting, sorting)]
            # print(F)
        return F, G

    @classmethod
    def renormalize_modes(cls, freqs, modes, inv, sorting=None, type=2):
        if type == 0:
            if sorting is not None:
                modes = modes[sorting, :]
                inv = inv[:, sorting]
            modes, inv = inv.T, modes.T
            freq = freqs
        else:
            F, G = cls.get_internal_FG(freqs, modes, inv, sorting=sorting)
            if type==1:
                G = np.linalg.inv(G)
            freq, modes = scipy.linalg.eigh(F, G, type=type)
            if type == 2:
                modes = modes * np.sqrt(freqs)[np.newaxis, :]
                inv = np.linalg.inv(modes)
            else:
                inv = (modes.T / np.sqrt(freqs)[:, np.newaxis])
                modes = np.linalg.inv(inv)
        return np.sqrt(freq), modes, inv

    @classmethod
    def rerotate_force_field(cls, old_inv, new_modes, old_field, dim_skips=0, sorting=None):
        mid_field = []
        final = -1 - dim_skips
        for f in old_field:
            for i in range(f.ndim - dim_skips):
                f = np.tensordot(old_inv, f, axes=[1, final])
                # print(np.round(f*219465))
            if sorting is not None:
                f = f[np.ix_(*(sorting,)*(f.ndim - dim_skips) )]
                # print(np.round(f*219465))
            mid_field.append(f)
        new_field = []
        for f in mid_field:
            for i in range(f.ndim - dim_skips):
                f = np.tensordot(new_modes, f, axes=[1, final])
            new_field.append(f)
        return new_field, mid_field

    @classmethod
    def reexpress_normal_modes(cls, base_modes, old_field, dipole, sorting=None, type=2):
        freq, matrix, inv = cls.renormalize_modes(*base_modes, sorting=sorting, type=type)
        potential_terms = cls.rerotate_force_field(
            base_modes[1],
            matrix.T,
            old_field,
            sorting=sorting
        )[0]
        if dipole is not None:
            # print(dipole[0].shape, base_modes[1].shape, matrix.shape)
            dipole = [
                    [0] + cls.rerotate_force_field(
                        base_modes[1],
                        matrix.T,
                        # /np.sqrt(freq)[np.newaxis, :],  # we divide by the sqrt of the frequencies to get the units to work
                        [d[a] for d in dipole],
                        sorting=sorting,
                        dim_skips=1
                        # [
                        #     np.diag(freq),  # in the file I was given the frequencies were in Hartree
                        #     cubics * helpers.convert("Wavenumbers", "Hartrees"),
                        #     quartics * helpers.convert("Wavenumbers", "Hartrees")
                        # ]
                    )[0]
            for a in range(3)
            ]
            # dipole = [np.zeros(3)] + dipole
        return (freq, matrix, inv), potential_terms, dipole

    # @classmethod
    # def rephase_normal_modes(cls, ...):

    convert = UnitsData.convert
    @staticmethod
    def mass(atom):
        return AtomData[atom]["Mass"] * UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")

    @classmethod
    def extract_term_lists(cls, checkpoint, terms, skip_dimensions=0, threshold=0, aggregator=None):
        from McUtils.Scaffolding import Checkpointer

        data = Checkpointer.from_file(checkpoint)
        terms = data[terms]
        if aggregator is not None:
            terms = aggregator(terms)

        term_files = []
        for t in terms:
            term_list = []
            ndim = t.ndim - skip_dimensions
            use_t = threshold<=0
            for idx in np.ndindex(*t.shape):
                if (use_t or np.abs(t[idx]) > threshold) and all(idx[i] <= idx[i+1] for i in range(skip_dimensions, ndim-1)):
                    term_list.append(tuple(i+1 for i in idx) + (t[idx],))
            term_files.append(term_list)
        return term_files

    @classmethod
    def write_term_lists(cls, terms, file_template=None, int_fmt="{:>3.0f}", float_fmt="{:>16.8e}", index_function=None):
        import io

        res = []
        if index_function is None:
            index_function = lambda i:i
        for i, t in enumerate(terms):
            if file_template is None:
                file = io.StringIO()
            else:
                file = file_template.format(index_function(i))
            with file if file_template is None else open(file, 'w+') as out:
                out.writelines(
                    " ".join((int_fmt if isinstance(x, int) else float_fmt).format(x) for x in row)+"\n"
                    for row in t
                )
            res.append(file)
        return res

    @classmethod
    def extract_terms(cls, chk, out, terms, default_output='output.hdf5', aggregator=None, index_function=None, skip_dimensions=0):
        if os.path.isdir(chk):
            woof = os.getcwd()
            try:
                os.chdir(chk)
                return cls.write_term_lists(
                    cls.extract_term_lists(default_output, terms, aggregator=aggregator, skip_dimensions=skip_dimensions),
                    file_template=out,
                    index_function=index_function
                )
            finally:
                os.chdir(woof)
        else:
            return cls.write_term_lists(
                cls.extract_term_lists(chk, terms, aggregator=aggregator, skip_dimensions=skip_dimensions),
                file_template=out,
                index_function=index_function
            )
    @classmethod
    def extract_potential(cls, chk, out='potential_expansion_{}.dat'):
        return cls.extract_terms(chk, out, 'potential_terms', index_function=lambda i:i+2)
    @classmethod
    def extract_gmatrix(cls, chk, out='gmatrix_expansion_{}.dat'):
        return cls.extract_terms(chk, out, 'gmatrix_terms')
    @classmethod
    def extract_dipole_expansion(cls, chk, out='dipole_expansion_{}.dat'):
        def agg(terms):
            return [
                np.array([terms[a][i] for a in ['x', 'y', 'z']])
                for i in range(len(terms['x']))
            ]
        return cls.extract_terms(chk, out, 'dipole_terms', aggregator=agg, skip_dimensions=1)

    @classmethod
    def run_anne_job(cls, base_dir,
                     states=2,
                     calculate_intensities=None,
                     return_analyzer=False,
                     return_runner=False,
                     modes_file='nm_int.dat',
                     atoms_file='atom.dat',
                     masses_file='mass.dat',
                     coords_file='cart_ref.dat',
                     zmat_file='z_mat.dat',
                     potential_files=('cub.dat', 'quart.dat', 'quintic.dat', 'sextic.dat'),
                     dipole_files=('lin_dip.dat', 'quad_dip.dat', "cub_dip.dat", "quart_dip.dat", 'quintic_dip.dat'),
                     results_file=None,#'output.hdf5',
                     order=None,
                     expansion_order=None,
                     energy_units=None,
                     type=0,
                     **opts
                     ):
        from .Analyzer import VPTAnalyzer

        og_dir = os.getcwd()
        try:
            if base_dir is not None:
                os.chdir(base_dir)

            base_freqs, base_mat, base_inv = cls.parse_modes(modes_file)
            if energy_units is None:
                if np.max(base_freqs) > 1:
                    conv = cls.convert("Wavenumbers", "Hartrees")
                else:
                    conv = 1
            else:
                conv = cls.convert(energy_units, "Hartrees")
            base_freqs *= conv
            potential = [cls.parse_tensor(f) for f in potential_files if os.path.isfile(f)]
            if energy_units is None:
                if np.max(np.abs(potential[0])) > 1:
                    conv = cls.convert("Wavenumbers", "Hartrees")
                else:
                    conv = 1
            else:
                conv = cls.convert(energy_units, "Hartrees")
            potential = [t * conv for t in potential]
            if os.path.exists(masses_file):
                masses = cls.parse_masses(masses_file)
            else:
                masses = None
            if masses is not None:
                masses = np.asanyarray(masses)
                if np.min(np.abs(masses)) < 100:
                    masses = cls.convert('AtomicMassUnits', 'AtomicUnitOfMass')*masses
            if masses is None or os.path.isfile(atoms_file):
                atoms = cls.parse_atoms(atoms_file)
            else:
                atoms = ['H']*len(masses)
            coords = cls.parse_coords(coords_file)
            if zmat_file is None:
                zmat = None
            else:
                if os.path.isfile(zmat_file):
                    zmat = cls.parse_zmatrix(zmat_file)
                else:
                    zmat = None
            sorting = cls.standard_sorting(zmat)  # we need to re-sort our internal coordinates
            if os.path.exists(dipole_files[0]):
                dipole_terms = [cls.parse_dipole_tensor(f) for f in dipole_files if os.path.isfile(f)]
                # if energy_units is None:
                #     if np.max(np.abs(potential[0])) > 1:
                #         conv = cls.convert("Wavenumbers", "Hartrees")
                #     else:
                #         conv = 1
                # else:
                #     conv = cls.convert(energy_units, "Hartrees")
            else:
                dipole_terms = None
            # raise Exception(sorting)
            (freq, matrix, inv), potential_terms, dipole_terms = cls.reexpress_normal_modes(
                (base_freqs, base_mat, base_inv),
                [np.diag(base_freqs)] + potential,
                dipole_terms,
                sorting=sorting,
                type=type
            )
            # if type == 0:
            #     # (freq, matrix, inv) = (base_freqs, base_mat, base_inv)
            #     potential_terms = [np.diag(base_freqs)] + potential
            # if dipole_terms is not None:
            #     dipole_terms = [
            #         [0] + [d[a] for d in dipole_terms]
            #         for a in range(3)
            #     ]

            if calculate_intensities is None:
                calculate_intensities = dipole_terms is not None

            # raise Exception(np.diag(potential_terms[0]) * UnitsData.convert("Hartrees", "Wavenumbers"), freq*UnitsData.convert("Hartrees", "Wavenumbers"))
            if return_analyzer:
                runner = lambda *a, **kw:VPTAnalyzer.run_VPT(*a, calculate_intensities=calculate_intensities, **kw)
            elif return_runner:
                runner = VPTRunner.construct
            else:
                runner = lambda *a, **kw:VPTRunner.run_simple(*a, calculate_intensities=calculate_intensities, **kw)

            if expansion_order is None:
                expansion_order = {}
            if isinstance(expansion_order, dict):
                if 'potential' not in expansion_order:
                    expansion_order['potential'] = len(potential_terms) - 1
                if 'kinetic' not in expansion_order:
                    expansion_order['kinetic'] = 2
                if dipole_terms is not None and 'dipole' not in expansion_order:
                    expansion_order['dipole'] = len(dipole_terms) - 1
            if order is None:
                if isinstance(expansion_order, int):
                    if expansion_order > 2:
                        order = expansion_order
                    else:
                        order = 2
                else:
                    order = expansion_order['potential']
            # raise Exception([a.shape for a in dipole_terms])

            res = runner(
                [atoms, coords, dict(masses=masses)],
                states,
                modes={
                    "freqs": freq,
                    "matrix": matrix,
                    "inverse": inv
                },
                potential_terms=potential_terms,
                dipole_terms=dipole_terms,
                internals=zmat,
                order=order,
                expansion_order=expansion_order,
                results=results_file,
                **opts
            )
            if return_analyzer:
                res.print_output_tables(
                    print_intensities=calculate_intensities,
                    print_energies=not calculate_intensities
                )

        finally:
            os.chdir(og_dir)

        return res
    @classmethod
    def run_fchk_job(cls,
                     base_dir,
                     states=2,
                     calculate_intensities=None,
                     return_analyzer=False,
                     return_runner=False,
                     # modes_file='nm_int.dat',
                     # atoms_file='atom.dat',
                     # masses_file='mass.dat',
                     # coords_file='cart_ref.dat',
                     zmat_file='z_mat.dat',
                     # potential_files=('cub.dat', 'quart.dat', 'quintic.dat', 'sextic.dat'),
                     # dipole_files=('lin_dip.dat', 'quad_dip.dat', "cub_dip.dat", "quart_dip.dat", 'quintic_dip.dat'),
                     fchk_file='fchk.fchk',
                     results_file='output.hdf5',
                     **opts
                     ):
        from .Analyzer import VPTAnalyzer

        curdir = os.getcwd()
        try:
            os.chdir(base_dir)

            if return_analyzer:
                runner = lambda *a, **kw:VPTAnalyzer.run_VPT(*a, calculate_intensities=calculate_intensities, **kw)
            elif return_runner:
                runner = VPTRunner.construct
            else:
                runner = lambda *a, **kw:VPTRunner.run_simple(*a, calculate_intensities=calculate_intensities, **kw)
            # raise Exception([a.shape for a in dipole_terms])

            if zmat_file is None:
                zmat = None
            else:
                if os.path.isfile(zmat_file):
                    zmat = cls.parse_zmatrix(zmat_file)
                else:
                    zmat = None

            res = runner(
                fchk_file,
                states,
                internals=zmat,
                results=results_file,
                **opts
            )
            if return_analyzer:
                res.print_output_tables(
                    print_intensities=calculate_intensities,
                    print_energies=not calculate_intensities
                )

        finally:
            os.chdir(curdir)

    @classmethod
    def get_internal_expansion(cls, fchk, internals, states=2, **opts):
        test_runner, _ = VPTRunner.construct(
            fchk,
            states,
            internals=internals,
            logger=False
        )
        test_freqs = test_runner.hamiltonian.modes.freqs
        test_V = test_runner.hamiltonian.V_terms
        test_modes = test_V.internal_L_matrix
        test_inver = test_V.internal_L_inverse
        return {
            "runner":test_runner,
            "freqs":test_freqs,
            "molecule":test_runner.hamiltonian.molecule,
            "kinetic":test_runner.hamiltonian.G_terms,
            "potential":test_V,
            "modes":[test_modes, test_inver],
            "states":states,
            "fchk":fchk,
            "zmatrix":internals
        }

    @classmethod
    def run_internal_expansion(cls, expansion_data, calculate_intensities=False, **opts):
        test_V = expansion_data['potential']
        return VPTRunner.run_simple(
            expansion_data['fchk'],
            expansion_data['states'],
            modes={
                "freqs": expansion_data['freqs'],
                "matrix": expansion_data['modes'][0],
                "inverse": expansion_data['modes'][1]
            },
            internals=expansion_data['zmatrix'],
            potential_terms=[test_V[0], test_V[1], test_V[2]],
            undimensionalize_normal_modes=False,
            calculate_intensities=calculate_intensities
        )

    # @classmethod
    # def make_nt_polyad(cls, nt):  #this should be of the form num of quanta give same E
    #     results=[]
    #     for i in range(len(nt)):
    #         for j in range(i+1,len(nt)):
    #             a=nt[i]
    #             b=nt[j]
    #             if a <= 0 or b <= 0:
    #                 continue  # this skips to the next step in the loop
    #             if a>b:
    #                 a, b = b, a
    #                 i, j = j, i
    #             if b%a==0:
    #                 b=b//a
    #                 a=1
    #             baselist=[0]*len(nt)
    #             baselist[i]=a
    #             twolist=[0]*len(nt)
    #             twolist[j]=b
    #             results.append([baselist,twolist])
    #     return results

VPTRunner.helpers = AnneInputHelpers