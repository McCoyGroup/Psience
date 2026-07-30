"""
A little package of utilities for setting up/running VPT jobs
"""
import numpy as np, sys, os, itertools, scipy, traceback as tb, math
from McUtils.Data import UnitsData, AtomData
from McUtils.Scaffolding import ParameterManager, Logger
from McUtils.Zachary import FiniteDifferenceDerivative
from McUtils.Extensions import ModuleLoader
import McUtils.Numputils as nput
import McUtils.Iterators as itut
from McUtils.Formatters import TableFormatter
import McUtils.Devutils as dev
import McUtils.Formatters as mfmt
from ..BasisReps import BasisStateSpace, HarmonicOscillatorProductBasis, BasisStateSpaceFilter, SelectionRuleStateSpace, StateMaker
from ..Molecools import Molecule
from ..Spectra import DiscreteSpectrum
from .DegeneracySpecs import DegeneracySpec, DegenerateMultiStateSpace
from .Hamiltonian import PerturbationTheoryHamiltonian
from .Analytic import PerturbationTheoryEvaluator, AnalyticPerturbationTheorySolver
from .Corrections import AnalyticPerturbationTheoryCorrections
VPTStateMaker = StateMaker
__all__ = ['VPTRunner', 'VPTSystem', 'VPTStateSpace', 'VPTStateMaker', 'VPTHamiltonianOptions', 'VPTRuntimeOptions', 'VPTSolverOptions', 'AnalyticVPTRunner']
__reload_hook__ = ['..BasisReps', '..Molecools', '.DegeneracySpecs', '.Hamiltonian', '.StateFilters', '.Analytic']

class VPTSystem:
    """
    Provides a little helper for setting up the input
    system for a VPT job

    :details:
    When using functions of internal (Z-matrix/polyspherical) coordinates, a sample form of the conversion function is
    ```python
    def conv(r, t, f, **kwargs):
        '''
        Takes the bond lengths (`r`), angles `(t)` and dihedrals `(f)`,
        and returns new coordinates that are functions of these coordinates
        '''
        ... # convert the coordinates
        return np.array([r, t, f])
    ```
    and then the inverse function will take the output of `conv` and return the original Z-matrix/polyspherical coordinates.
    """
    __props__ = ('internals', 'modes', 'local_modes', 'mode_selection', 'mode_transformation', 'full_surface_mode_selection', 'potential_derivatives', 'potential_function', 'order', 'dipole_derivatives', 'dummy_atoms', 'eckart_embed')

    def __init__(self, mol, internals=None, dummy_atoms=None, modes=None, local_modes=None, mode_selection=None, mode_transformation=None, full_surface_mode_selection=None, potential_derivatives=None, potential_function=None, order=2, dipole_derivatives=None, eckart_embed=False, copy_mol=False):
        """
        :param mol: the molecule or system specification to use (doesn't really even need to be a molecule)
        :type mol: str | list | Molecule
        :param internals: the Z-matrix for the internal coordinates optionally with a specification of a `conversion` and `inverse`
        To supply a conversion function, provide a `dict` like so
        ```python
        {
            'zmatrix': [[atom1, bond1, angle1, dihed1], [atom2, bond2, angle2, dihed2], ...] or None,
            'conversion': 'a function to convert from Z-matrix coordinates to desired coordinates',
            'inverse': 'the inverse conversion'
        }
        ```
        :type internals: list | dict
        :param modes: the normal modes to use if not already supplied by the Molecule
        :type modes: MolecularVibrations|dict
        :param potential_derivatives: the derivatives of the potential to use for expansions
        :type potential_derivatives: Iterable[np.ndarray]
        :param dipole_derivatives: the set of dipole derivatives to use for expansions
        :type dipole_derivatives: Iterable[np.ndarray]
        """
        ...

    def prep_local_modes(self, dRdX, dXdR=None, sort_freqs=False):
        """
        **LLM Docstring**

        Build a set of "local mode" normal-mode data (frequencies, mode matrix, and its inverse) from a set of Cartesian-to-internal Jacobians, by rescaling the force-constant and G-matrix diagonals into a locally-diagonal (Duschinsky-like) basis.

        :param dRdX: the internals-by-Cartesians Jacobian; used to derive `dXdR` if not given
        :type dRdX: np.ndarray | None
        :param dXdR: the Cartesians-by-internals Jacobian; derived from `dRdX` via pseudo-inverse if not given
        :type dXdR: np.ndarray | None
        :param sort_freqs: whether to sort the resulting modes by ascending frequency
        :type sort_freqs: bool
        :return: a dict with `'matrix'`, `'inverse'`, and `'freqs'` keys describing the local-mode basis
        :rtype: dict
        """
        ...

    @property
    def nmodes(self):
        """
        Provides the number of modes in the system

        :return:
        :rtype:
        """
        ...

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
        ...

    @classmethod
    def from_harmonic_scan(cls, scan_array):
        """
        **LLM Docstring**

        Intended to build a `VPTSystem` from a harmonic potential-energy scan array. Not implemented.

        :param scan_array: the scan data to build the system from
        :type scan_array: object
        :return: never returns
        :rtype: VPTSystem
        :raises NotImplementedError: always
        """
        ...

class VPTStateSpace:
    """
    Provides a helper to make it easier to set up the input
    state spaces/degenerate spaces to run the perturbation theory

    :details:

    There are multiple possible values for the `degeneracy_specs`.
    The simplest is to use the automatic approach, in which we supply a numeric type (`int`, `float`, etc.) to use as the `WFC` threshold.
    The next simplest is to explicitly supply the groups we want, like

    ```python
    [
        [ # the first resonant space
            state_1,
            state_2,
            state_3
        ],
        [ # the second
            state_5, state_11, ...
        ],
        ...
    ]
    ```

    We can also supply pairs of relations for determining resonances, like

    ```python
    [
        [state_1,  state_2], # A first relation
        [state_3,  state_4],  # another relation
        ...
    ]
    ```

    To allow for extra options, you can also supply a `dict`. If you wanted to have a different `wfc_threshold` and you wanted to do the secondary resonant space splitting step with a very large threshold, you could do that by supplying

    ```python
    {
        'wfc_threshold':.1,
        'energy_cutoff':1.0 # in Hartree
    }
    ```

    or you can explicitly add extra groups to the pairs of polyad rules by saying

    ```python
    {
        'polyads':[
                [state_1,  state_2], # A first relation
                [state_3,  state_4],  # another relation
                ...
            ],
        'extra_groups': [
            [ # the first resonant space
                state_a,
                state_b,
                state_c
            ],
            [ # the second
                state_d, state_e, ...
            ],
            ...
        ]
    }
    ```

    This also allows us to define more resonance handling strategies.

    The Martin Test is supported,
    ```python
    {
        'martin_threshold':.1/219465, #in Hartree
    }
    ```

    As are total quanta vectors/polyads
    ```python
    {
        'nT': [1, 1, 1, 0, 2, 2, 0] # e.g.
    }
    ```
    """
    __props__ = ('degeneracy_specs',)

    def __init__(self, states, degeneracy_specs=None, system=None, frequencies=None, evaluator=None):
        """
        :param states: A list of states or a number of quanta to target
        :type states: list | int
        :param degeneracy_specs: A specification of degeneracies, either as polyads, explicit groups of states, or parameters to a method. (see Details for more info)
        :type degeneracy_specs: 'auto' | list | dict
        """
        ...

    @classmethod
    def from_system_and_spec(cls, system, spec, **opts):
        """
        **LLM Docstring**

        Build a `VPTStateSpace` from a system and a flexible state specification: passes an already-built `VPTStateSpace` through unchanged, dispatches a bare integer to `from_system_and_quanta` (states up to that many quanta), and otherwise treats `spec` as an explicit state list.

        :param system: the system (molecule/mode context) the states are defined relative to
        :type system: VPTSystem
        :param spec: the state specification: an existing `VPTStateSpace`, an integer quanta cutoff, or an explicit list of states
        :type spec: VPTStateSpace | int | Iterable
        :param opts: extra options forwarded to the constructor used
        :type opts: dict
        :return: the resolved state space
        :rtype: VPTStateSpace
        """
        ...

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
        ...

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
        ...

    def build_degenerate_state_spaces(self, degeneracy_specs, states, system=None, evaluator=None, freqs=None) -> '(None|DegeneracySpec, None|list[np.ndarray])':
        """
        :param degeneracy_specs:
        :type degeneracy_specs:
        :return:
        :rtype:
        """
        ...

    def filter_generator(self, target_property, order=2, initial_states=None, postfilters=None):
        """
        **LLM Docstring**

        Build a callable that produces a state-space filter for a given target property, by binding `target_property`/`order`/`initial_states`/`postfilters` and delegating each call to `get_state_space_filter`; used to supply `state_space_filter_generator`-style hooks that need to be re-evaluated against different candidate state sets.

        :param target_property: the property (e.g. `'intensities'`) the filter should be built for
        :type target_property: str
        :param order: the perturbation-theory order to build the filter for
        :type order: int
        :param initial_states: the initial states the filtered property is computed relative to
        :type initial_states: Iterable | None
        :param postfilters: extra post-transformation filters to apply
        :type postfilters: object | None
        :return: the constructed filter-generating function
        :rtype: callable
        """
        ...

    def get_filter(self, target_property, order=2, initial_states=None, postfilters=None):
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
        ...

    @classmethod
    def get_state_space_filter(cls, states, initial_states=None, n_modes=None, order=2, target='wavefunctions', postfilters=None, **opts):
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
        ...

class VPTHamiltonianOptions:
    """
    Provides a helper to keep track of the levers available for
    setting up the Hamiltonian
    """
    __props__ = ('mode_selection', 'mode_transformation', 'local_mode_couplings', 'local_mode_coupling_order', 'full_surface_mode_selection', 'include_potential', 'include_gmatrix', 'include_coriolis_coupling', 'include_pseudopotential', 'include_only_mode_couplings', 'potential_terms', 'kinetic_terms', 'coriolis_terms', 'pseudopotential_terms', 'include_dipole', 'dipole_terms', 'dipole_derivatives', 'undimensionalize_normal_modes', 'use_numerical_jacobians', 'eckart_embed_derivatives', 'eckart_embed_planar_ref_tolerance', 'strip_dummy_atoms', 'strip_embedding_coordinates', 'mixed_derivative_handling_mode', 'mixed_derivative_warning_threshold', 'mixed_derivative_handle_zeros', 'rephase_modes', 'backpropagate_internals', 'direct_propagate_cartesians', 'zero_mass_term', 'use_internal_modes', 'internal_fd_mesh_spacing', 'internal_fd_stencil', 'cartesian_fd_mesh_spacing', 'cartesian_fd_stencil', 'cartesian_analytic_deriv_order', 'cartesian_by_internal_derivative_method', 'internal_by_cartesian_order', 'cartesian_by_internal_order', 'expansion_handling_mode', 'jacobian_warning_threshold', 'check_input_force_constants', 'hessian_tolerance', 'grad_tolerance', 'freq_tolerance', 'g_derivative_threshold', 'gmatrix_tolerance', 'use_cartesian_kinetic_energy', 'operator_coefficient_threshold', 'imaginary_frequency_handling_mode')

    def __init__(self, mode_selection=None, mode_transformation=None, local_mode_couplings=None, local_mode_coupling_order=None, full_surface_mode_selection=None, include_potential=None, include_gmatrix=None, include_coriolis_coupling=None, include_pseudopotential=None, include_only_mode_couplings=None, potential_terms=None, kinetic_terms=None, coriolis_terms=None, pseudopotential_terms=None, include_dipole=None, dipole_terms=None, dipole_derivatives=None, undimensionalize_normal_modes=None, use_numerical_jacobians=None, eckart_embed_derivatives=None, eckart_embed_planar_ref_tolerance=None, strip_dummy_atoms=None, strip_embedding_coordinates=None, mixed_derivative_handling_mode=None, mixed_derivative_warning_threshold=None, mixed_derivative_handle_zeros=None, rephase_modes=None, backpropagate_internals=None, direct_propagate_cartesians=None, zero_mass_term=None, internal_fd_mesh_spacing=None, internal_fd_stencil=None, cartesian_fd_mesh_spacing=None, cartesian_fd_stencil=None, cartesian_analytic_deriv_order=None, cartesian_by_internal_derivative_method=None, internal_by_cartesian_order=None, cartesian_by_internal_order=None, expansion_handling_mode=None, jacobian_warning_threshold=None, check_input_force_constants=None, hessian_tolerance=None, grad_tolerance=None, freq_tolerance=None, g_derivative_threshold=None, gmatrix_tolerance=None, use_internal_modes=None, use_cartesian_kinetic_energy=None, operator_coefficient_threshold=None, imaginary_frequency_handling_mode=None):
        """
        :param mode_selection: the set of the supplied normal modes to do perturbation theory on
        :type mode_selection: Iterable[int]|None
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
        :param operator_coefficient_threshold: the minimum size of a coefficient to keep when evaluating representation terms
        :type operator_coefficient_threshold: float|None
        """
        ...

class VPTRuntimeOptions:
    """
    Provides a helper to keep track of the options available
    for configuring the way the code runs
    """
    __props__ = ('operator_chunk_size', 'matrix_element_threshold', 'logger', 'verbose', 'checkpoint', 'results', 'memory_constrained', 'checkpoint_keys', 'use_cached_representations', 'use_cached_basis', 'nondeg_hamiltonian_precision')

    def __init__(self, operator_chunk_size=None, matrix_element_threshold=None, nondeg_hamiltonian_precision=None, logger=None, verbose=None, checkpoint=None, results=None, parallelizer=None, memory_constrained=None, checkpoint_keys=None, use_cached_representations=None, use_cached_basis=None):
        """
        :param operator_chunk_size: the number of representation matrix elements to calculate in at one time
        :type operator_chunk_size: int|None default:None
        :param matrix_element_threshold: the minimum size of matrix element to keep
        :type matrix_element_threshold: float|None default:None
        :param nondeg_hamiltonian_precision: the precision with which to print out elements in the degenerate coupling Hamiltonians in the log file
        :type nondeg_hamiltonian_precision: int
        :param logger: the `Logger` object to use when logging the status of the calculation (`True` means log normally)
        :type logger: str|Logger|bool|None default:None
        :param results: the `Checkpointer` to write corrections out to
        :type results: str|Checkpointer|None default:None
        :param parallelizer: the `Parallelizer` to use for parallelizing the evaluation of matrix elements
        :type parallelizer: Parallelizer|None default:None
        :param memory_constrained: whether or not to attempt memory optimizations (`None` means attempt for >20D problems)
        :type memory_constrained: bool|None
        :param checkpoint: the `Checkpointer` to write Hamiltonians and other bits out to
        :type checkpoint: str|Checkpointer|None default:None
        :param checkpoint_keys: which keys to save in the checkpoint
        :type checkpoint_keys: Iterable[str]|None
        :param use_cached_representations: whether other not to use Hamiltonian reps from the checkpoint
        :type use_cached_representations: bool
        :param use_cached_basis: whether other not to use bases from the checkpoint
        :type use_cached_basis: bool
        """
        ...

class VPTSolverOptions:
    """
    Provides a helper to keep track of the options available
    for configuring the way the perturbation theory is applied

    :details:
    The `basis_postfilters` have multiple possible values.
    Here are the currently supported cases

    ```python
    {
        'max_quanta': [2, -1, 1, -1, ...] # the max number of quanta allowed in a given mode in the basis (-1 means infinity)
    }
    ```

    - for excluding transitions

    ```python
    {
        'excluded_transitions': [[0, 0, 1, 0, ...], [1, 0, 0, 0, ...], ...] # a set of transitions that are forbidden on the input states
    }
    ```

    - for excluding based on a test

    ```python
    {
        'test': func # a function that takes the basis and tests if states should be allowed
    }
    ```
    """
    __props__ = ('order', 'expansion_order', 'coupled_states', 'total_space', 'flat_total_space', 'state_space_iterations', 'state_space_terms', 'state_space_filters', 'extended_state_space_filter_generator', 'extended_state_space_postprocessor', 'allow_post_PT_calc', 'modify_degenerate_perturbations', 'gaussian_resonance_handling', 'ignore_odd_order_energies', 'intermediate_normalization', 'zero_element_warning', 'degenerate_states', 'zero_order_energy_corrections', 'handle_strong_couplings', 'low_frequency_mode_cutoff', 'check_overlap')

    def __init__(self, order=2, expansion_order=None, coupled_states=None, total_space=None, flat_total_space=None, state_space_iterations=None, state_space_terms=None, state_space_filters=None, extended_state_space_filter_generator=None, extended_state_space_postprocessor=None, allow_post_PT_calc=None, modify_degenerate_perturbations=None, gaussian_resonance_handling=None, ignore_odd_order_energies=None, intermediate_normalization=None, zero_element_warning=None, degenerate_states=None, handle_strong_couplings=None, strong_coupling_test_modes=None, strong_couplings_state_filter=None, strongly_coupled_group_filter=None, extend_strong_coupling_spaces=None, strong_coupling_zero_order_energy_cutoff=None, low_frequency_mode_cutoff=None, zero_order_energy_corrections=None, check_overlap=None):
        """
        :param order: the order of perturbation theory to apply
        :type order: int
        :type order: int
        :param expansion_order: the order to go to in the expansions of the perturbations, this can be supplied for different properties independently, like
        ```python
        expansion_order = {
            'potential':some_int,
            'kinetic':some_int,
            'dipole':some_int
        }
        ```
        :type expansion_order: int | dict
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
        :param low_frequency_mode_cutoff: the energy below which to consider a mode to be "low frequency"
        :type low_frequency_mode_cutoff: float (default:500 cm-1)
        :param zero_order_energy_corrections: energies to use for the zero-order states instead of the diagonal of `H(0)`
        :type zero_order_energy_corrections: dict
        :param check_overlap: whether or not to ensure states are normalized in the VPT
        :type check_overlap: bool default:True
        """
        ...

    @staticmethod
    def _harmonic_energies(corrected_fundamental_freqs, states):
        """

        :param corrected_fundamental_freqs:
        :type corrected_fundamental_freqs:
        :param states:
        :type states:
        :return:
        :rtype:
        """
        ...

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
        ...

class VPTRunner:
    """
    A helper class to make it easier to run jobs by making the inputs/options
    clear and making it easier to customize run options
    """

    def __init__(self, system, states, initial_states=None, hamiltonian_options=None, solver_options=None, runtime_options=None):
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
        ...

    def get_Hamiltonian(self):
        """
        **LLM Docstring**

        Build a `PerturbationTheoryHamiltonian` for this runner's system, combining the configured Hamiltonian and runtime options.

        :return: the constructed Hamiltonian
        :rtype: PerturbationTheoryHamiltonian
        """
        ...

    @property
    def hamiltonian(self):
        """
        **LLM Docstring**

        The (cached) `PerturbationTheoryHamiltonian` for this runner, built lazily via `get_Hamiltonian` the first time it's needed.

        :return: the Hamiltonian
        :rtype: PerturbationTheoryHamiltonian
        """
        ...

    def get_wavefunctions(self, **opts):
        """
        **LLM Docstring**

        Run the full VPT calculation and return the resulting wavefunctions, combining the solver options, degenerate-state/degeneracy-handler settings, and runtime options configured on this runner (with any explicitly passed `opts` taking precedence).

        :param opts: extra options overriding the runner's configured solver/runtime options
        :type opts: dict
        :return: the computed perturbation-theory wavefunctions
        :rtype: PerturbationTheoryWavefunctions
        """
        ...

    def get_solver(self):
        """
        **LLM Docstring**

        Build a `PerturbationTheorySolver` for this runner's target states, without running the full wavefunction calculation, combining the configured solver and runtime options.

        :return: the constructed solver
        :rtype: PerturbationTheorySolver
        """
        ...

    @classmethod
    def print_output_tables(cls, wfns=None, file=None, print_intensities=True, print_energies=True, print_energy_corrections=True, print_transition_moments=True, operators=None, logger=None, sep_char='=', sep_len=100):
        """
        Prints a bunch of formatted output data from a PT run

        :param wfns:
        :type wfns:
        :return:
        :rtype:
        """
        ...

    def print_tables(self, wfns=None, file=None, print_intensities=True, print_energy_corrections=True, print_transition_moments=True, operators=None, logger=None, sep_char='=', sep_len=100):
        """
        Prints a bunch of formatted output data from a PT run

        :param wfns:
        :type wfns:
        :return:
        :rtype:
        """
        ...

    def get_Nielsen_energies(self, return_split=False, return_X=False, **potential_params):
        """
        **LLM Docstring**

        Compute harmonic and anharmonic (Nielsen-formula) vibrational energies for the target states, via the Hamiltonian's own `get_Nielsen_energies`.

        :param return_split: whether to also return the anharmonicity split out separately
        :type return_split: bool
        :param return_X: whether to also return the underlying X anharmonicity-constant matrix
        :type return_X: bool
        :param potential_params: extra options forwarded to `PerturbationTheoryHamiltonian.get_Nielsen_energies`
        :type potential_params: dict
        :return: `(harmonic, total)` energies, or `(harmonic, total, x_matrix)` if `return_split`/`return_X` is set
        :rtype: tuple
        """
        ...

    def print_Nielsen_frequencies(self, logger=None, state_formatting='vector', **potential_params):
        """
        **LLM Docstring**

        Compute the Nielsen-formula harmonic and anharmonic transition frequencies (relative to the ground state) and print them as a formatted state/frequency table.

        :param logger: `None` to print directly, `True` to use the Hamiltonian's own logger, or an explicit logger to print into a block
        :type logger: Logger | bool | None
        :param state_formatting: `'vector'` to display states as raw excitation vectors, or another format understood by `VPTStateMaker.parse_state`
        :type state_formatting: str
        :param potential_params: extra options forwarded to `get_Nielsen_energies`
        :type potential_params: dict
        :return: None
        :rtype: None
        """
        ...

    @classmethod
    def _quanta_coupled_state_spaces(cls, states, order, track=False):
        """
        **LLM Docstring**

        Build, per perturbative order, the set of states coupled to each input state purely by total-quantum-number-change parity/magnitude rules (states within `order` quanta of each other, alternating between odd/even changes by order), used as a coupling space to extend for degenerate treatment.

        :param states: the reference states to build coupled spaces for
        :type states: BasisStateSpace
        :param order: the highest perturbative order to build coupled spaces for
        :type order: int
        :param track: whether to precompute change-index tracking on the resulting spaces
        :type track: bool
        :return: the per-order list of `SelectionRuleStateSpace` coupled spaces
        :rtype: list[SelectionRuleStateSpace]
        """
        ...

    @classmethod
    def _matrix_element_filler(cls, new_states, new_spaces, degeneracy_spec):
        """
        **LLM Docstring**

        Extend a set of coupled-state spaces to additionally include quantum-number-based coupling states (via `_quanta_coupled_state_spaces`), either across a single degeneracy-spec-free set of states or per degenerate group. Note the `degeneracy_spec is None` branch contains an unconditional `raise Exception(qspaces, new_spaces)` placed before the union logic that would otherwise complete it, so that branch currently always raises rather than returning extended spaces -- this looks like leftover debugging code.

        :param new_states: the states to compute quantum-number-based couplings for
        :type new_states: BasisStateSpace
        :param new_spaces: the existing per-order coupled-state spaces to extend
        :type new_spaces: list[SelectionRuleStateSpace]
        :param degeneracy_spec: the degenerate groups to compute couplings within, or `None` to compute them across all of `new_states` at once
        :type degeneracy_spec: Iterable[BasisStateSpace] | None
        :return: the extended per-order coupled-state spaces
        :rtype: list[SelectionRuleStateSpace]
        :raises Exception: unconditionally, if `degeneracy_spec` is `None` (leftover debug code)
        """
        ...

    @classmethod
    def construct(cls, system, states, target_property=None, extended_space_target_property=None, basis_filters=None, initial_states=None, corrected_fundamental_frequencies=None, **opts):
        """
        **LLM Docstring**

        Top-level constructor that assembles a fully configured `VPTRunner` (and, if target states are given, its resolved `VPTStateSpace`) from a molecule/system specification, target states, and a large set of options split automatically across `VPTSystem`, `VPTStateSpace`, `VPTHamiltonianOptions`, `VPTRuntimeOptions`, and `VPTSolverOptions`.

        :param system: the molecule or system specification to run VPT on
        :type system: str | list | Molecule | VPTSystem
        :param states: the target states, given as a states object, quanta cutoff, or explicit list
        :type states: VPTStateSpace | int | Iterable
        :param target_property: the property (e.g. `'intensities'`) used to build a default state-space filter if `state_space_filters` isn't explicitly given
        :type target_property: str | None
        :param extended_space_target_property: accepted for interface consistency but not directly used in this method's body
        :type extended_space_target_property: object | None
        :param basis_filters: extra post-transformation filters merged into the generated state-space filter
        :type basis_filters: object | None
        :param initial_states: the initial states used both for filter construction and as the runner's own initial states
        :type initial_states: VPTStateSpace | int | Iterable | None
        :param corrected_fundamental_frequencies: replacement fundamental frequencies to use for the underlying molecule's normal modes
        :type corrected_fundamental_frequencies: np.ndarray | None
        :param opts: the full set of system/state-space/Hamiltonian/runtime/solver options, validated against the union of all their `__props__`
        :type opts: dict
        :return: the constructed runner, or `(runner, states)` if target states were given
        :rtype: VPTRunner | tuple[VPTRunner, VPTStateSpace]
        :raises ValueError: if `opts` contains any key not recognized by any of the option classes
        """
        ...

    @classmethod
    def run_simple(cls, system, states, target_property=None, corrected_fundamental_frequencies=None, calculate_intensities=True, plot_spectrum=False, operators=None, **opts):
        """
        The standard runner for VPT.
        Makes a runner using the `construct` method and then calls that
        runner's `print_tables` method after printing out run info.

        :param system: the system spec, either as a `Molecule`, molecule spec (atoms, coords, opts) or a file to construct a `Molecule`
        :type system: list|str|Molecule
        :param states: the states to get corrections for either an `int` (up to that many quanta) or an explicit state list
        :type states: int|list
        :param target_property: the target property to get corrections for (one of 'frequencies', 'intensities', 'wavefunctions')
        :type target_property: str
        :param corrected_fundamental_frequencies: a set of fundamental frequencies to use to get new zero-order energies
        :type corrected_fundamental_frequencies: Iterable[float]|None
        :param calculate_intensities: whether or not to calculate energies
        :type calculate_intensities: bool default:True
        :param opts: options that work for a `VPTSystem`, `VPTStateSpace`, `VPTRuntimeOptions`, `VPTSolverOptions`, or `VPTHamiltonianOptions` object which will be filtered automatically
        :type opts:
        """
        ...

    class helpers:
        """
        A stub to be replaced with the AnneInputHelpers interface
        """

        @classmethod
        def run_anne_job(cls, base_dir, states=2, calculate_intensities=None, return_analyzer=False, return_runner=False, modes_file=('nm_int.dat', 'modes.dat'), atoms_file='atom.dat', masses_file='mass.dat', coords_file='cart_ref.dat', zmat_file='z_mat.dat', potential_files=('cub.dat', 'quart.dat', 'quintic.dat', 'sextic.dat'), dipole_files=('lin_dip.dat', 'quad_dip.dat', 'cub_dip.dat', 'quart_dip.dat', 'quintic_dip.dat'), coordinate_transformation=None, coordinate_transformation_file='coordinate_transformation.py', results_file=None, order=None, expansion_order=None, energy_units=None, normalization_type=0, **opts):
            """
            **LLM Docstring**

            Stub placeholder documenting the intended `AnneInputHelpers.run_anne_job` interface for running a VPT job from Anne-format input files; the actual implementation lives on `AnneInputHelpers` and is bound onto `VPTRunner.helpers` at import time, so this method's body (`...`) is never executed.

            :param base_dir: the directory containing the Anne-format input files
            :type base_dir: str
            :param states: the target states (or quanta cutoff) to compute
            :type states: int | Iterable
            :param calculate_intensities: whether to compute IR intensities
            :type calculate_intensities: bool | None
            :param return_analyzer: whether to return a `VPTAnalyzer` instead of raw results
            :type return_analyzer: bool
            :param return_runner: whether to return the constructed `VPTRunner` instead of raw results
            :type return_runner: bool
            :param modes_file: the normal-modes input file(s)
            :type modes_file: str | tuple
            :param atoms_file: the atoms input file
            :type atoms_file: str
            :param masses_file: the masses input file
            :type masses_file: str
            :param coords_file: the Cartesian-coordinates input file
            :type coords_file: str
            :param zmat_file: the Z-matrix input file
            :type zmat_file: str
            :param potential_files: the potential-derivative input files
            :type potential_files: tuple
            :param dipole_files: the dipole-derivative input files
            :type dipole_files: tuple
            :param coordinate_transformation: an explicit `[conversion, inverse]` coordinate-transformation pair
            :type coordinate_transformation: list | None
            :param coordinate_transformation_file: a Python module file defining `conversion`/`inverse` functions
            :type coordinate_transformation_file: str
            :param results_file: where to store the results checkpoint
            :type results_file: str | None
            :param order: the perturbation-theory order
            :type order: int | None
            :param expansion_order: the per-term expansion orders
            :type expansion_order: int | dict | None
            :param energy_units: the units the input energies are given in
            :type energy_units: str | None
            :param normalization_type: which mode-renormalization convention to use
            :type normalization_type: int
            :param opts: extra options forwarded to the underlying runner
            :type opts: dict
            :return: never actually executed (stub body is `...`)
            :rtype: object
            """
            ...
        convert = UnitsData.convert

class AnneInputHelpers:

    @classmethod
    def _check_file(cls, no_file):
        """
        **LLM Docstring**

        Sanity-check that a supposed file path actually looks like a file (rather than raw file content passed by mistake): raises if the string has fewer than 2 lines and contains a `.`, which would be typical of a bare (nonexistent) filename rather than real multi-line data.

        :param no_file: the string that was expected to be an existing file path but wasn't found on disk
        :type no_file: str
        :return: None
        :rtype: None
        :raises FileNotFoundError: if `no_file` looks like a filename rather than raw data
        """
        ...

    @classmethod
    def get_tensor_idx(cls, line, inds, m, start_at=0):
        """
        **LLM Docstring**

        Parse one line of a flattened force-constant/tensor data file into its (1-indexed-to-0-indexed) tuple of indices and its value, inserting the value into the running `inds` dict and updating the running maximum index seen.

        :param line: the raw data line to parse, whitespace-separated indices followed by a value
        :type line: str
        :param inds: the running dict mapping index tuples to values, updated in place
        :type inds: dict
        :param m: the running maximum index seen so far
        :type m: int
        :param start_at: which index position to start considering for the running maximum (used to skip a leading Cartesian-component index for dipole tensors)
        :type start_at: int
        :return: the updated maximum index
        :rtype: int
        """
        ...

    @classmethod
    def parse_tensor(cls, block, dims=None):
        """
        **LLM Docstring**

        Parse a symmetric force-constant tensor from an Anne-format data file (or its raw string content), filling in every index permutation from the (upper-triangular-only) parsed entries.

        :param block: the file path, or the raw file content as a string
        :type block: str
        :param dims: the tensor's shape; inferred from the largest parsed index and the number of index columns if not given
        :type dims: tuple | None
        :return: the assembled, fully-symmetrized tensor
        :rtype: np.ndarray
        """
        ...

    @classmethod
    def parse_dipole_tensor(cls, block, dims=None):
        """
        **LLM Docstring**

        Parse a dipole-derivative tensor from an Anne-format data file (or its raw string content), treating the first index column as the Cartesian (x/y/z) component and symmetrizing over the remaining (mode) indices.

        :param block: the file path, or the raw file content as a string
        :type block: str
        :param dims: the tensor's shape (Cartesian-component axis first); inferred from the largest parsed mode index if not given
        :type dims: tuple | None
        :return: the assembled dipole-derivative tensor
        :rtype: np.ndarray
        """
        ...

    @classmethod
    def parse_freqs_line(cls, line):
        """
        **LLM Docstring**

        Parse a whitespace-separated line of frequency values into a NumPy array.

        :param line: the raw line to parse
        :type line: str
        :return: the parsed frequencies, or `None` if the line was empty
        :rtype: np.ndarray | None
        """
        ...

    @classmethod
    def parse_modes_line(cls, line, nmodes):
        """
        **LLM Docstring**

        Parse a whitespace-separated, flattened block of mode-matrix values into a `(nmodes, ncols)` array, inferring the number of columns from the total value count and the known number of modes.

        :param line: the raw (possibly multi-line, joined) block of values to parse
        :type line: str
        :param nmodes: the number of modes (rows, before transposing) the flattened data should reshape into
        :type nmodes: int
        :return: the parsed, transposed mode matrix, or `None` if the line was empty
        :rtype: np.ndarray | None
        """
        ...

    @classmethod
    def _parse_modes(cls, line_iter, energy_units=None):
        """
        **LLM Docstring**

        Parse the frequencies, mode matrix, and inverse mode matrix out of a sequence of file lines (separated into blank-line-delimited blocks), deriving the inverse from the mode matrix's transpose (with an appropriate unit conversion) if the file doesn't explicitly provide a third block.

        :param line_iter: an iterator over the raw lines to parse
        :type line_iter: Iterable[str]
        :param energy_units: the units the parsed frequencies are given in, used to derive the fallback inverse-matrix scaling; if `None`, units are inferred from the frequency magnitudes
        :type energy_units: str | None
        :return: `(freqs, L, Linv)` -- the frequencies, mode matrix, and its inverse
        :rtype: tuple[np.ndarray, np.ndarray, np.ndarray]
        """
        ...

    @classmethod
    def parse_modes(cls, block, energy_units=None):
        """
        **LLM Docstring**

        Parse the frequencies/mode matrix/inverse from an Anne-format modes file, or the first file that exists among a list of candidate filenames.

        :param block: the file path (or raw file content), or an iterable of candidate file paths to try in order
        :type block: str | Iterable[str]
        :param energy_units: the units the frequencies are given in
        :type energy_units: str | None
        :return: `(freqs, L, Linv)`
        :rtype: tuple[np.ndarray, np.ndarray, np.ndarray]
        """
        ...

    @classmethod
    def parse_coords(cls, block):
        """
        **LLM Docstring**

        Parse a Cartesian-coordinates data file (or its raw string content) into an `(natoms, 3)` array.

        :param block: the file path, or the raw file content as a string
        :type block: str
        :return: the parsed coordinates
        :rtype: np.ndarray
        """
        ...

    @classmethod
    def parse_atoms(cls, block):
        """
        **LLM Docstring**

        Parse an atomic-number data file (or its raw string content) into a list of element symbols.

        :param block: the file path, or the raw file content as a string
        :type block: str
        :return: the parsed element symbols
        :rtype: list[str]
        """
        ...

    @classmethod
    def parse_masses(cls, block):
        """
        **LLM Docstring**

        Parse a masses data file (or its raw string content) into a list of floating-point mass values.

        :param block: the file path, or the raw file content as a string
        :type block: str
        :return: the parsed masses
        :rtype: list[float]
        """
        ...

    @classmethod
    def parse_zmatrix(cls, block):
        """
        **LLM Docstring**

        Parse a Z-matrix connectivity data file (or its raw string content) into a standard `(natoms, 4)` Z-matrix ordering array, prepending a placeholder dummy-origin row and padding each atom's row out to 4 columns (dropping a redundant 4th input column, if present).

        :param block: the file path, or the raw file content as a string
        :type block: str
        :return: the parsed Z-matrix ordering array
        :rtype: np.ndarray
        """
        ...

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
        ...

    @classmethod
    def get_internal_FG(cls, freqs, modes, inv, sorting=None):
        """
        **LLM Docstring**

        Build the internal-coordinate force-constant (`F`) and kinetic-energy (`G`) matrices from a set of normal-mode frequencies and (mass-weighted, dimensionless) mode/inverse matrices, optionally reordering both to a standard coordinate ordering.

        :param freqs: the mode frequencies
        :type freqs: np.ndarray
        :param modes: the mode matrix (frequency-undimensionalized inside this method)
        :type modes: np.ndarray
        :param inv: the inverse mode matrix
        :type inv: np.ndarray
        :param sorting: a permutation to reorder the resulting `F`/`G` matrices' rows/columns
        :type sorting: np.ndarray | None
        :return: `(F, G)` -- the internal-coordinate force-constant and kinetic-energy matrices
        :rtype: tuple[np.ndarray, np.ndarray]
        """
        ...

    @classmethod
    def renormalize_modes(cls, freqs, modes, inv, sorting=None, type=2):
        """
        **LLM Docstring**

        Re-derive a set of normal modes (with a chosen dimensionality convention) from the internal-coordinate `F`/`G` matrices, via a generalized eigenvalue solve, optionally reordering to a standard coordinate ordering first.

        :param freqs: the original mode frequencies
        :type freqs: np.ndarray
        :param modes: the original mode matrix
        :type modes: np.ndarray
        :param inv: the original inverse mode matrix
        :type inv: np.ndarray
        :param sorting: a permutation to reorder the coordinates before re-solving
        :type sorting: np.ndarray | None
        :param type: the mode-dimensionality convention: `0` to just reorder/transpose the given matrices directly, `1`/`2` to re-solve the generalized eigenproblem `F v = freq^2 G v` (in the `G^-1`-inverted or standard form, respectively)
        :type type: int
        :return: `(freq, modes, inv)` -- the (square-rooted) frequencies and the renormalized mode/inverse matrices
        :rtype: tuple[np.ndarray, np.ndarray, np.ndarray]
        """
        ...

    @classmethod
    def rerotate_force_field(cls, old_inv, new_modes, old_field, dim_skips=0, sorting=None):
        """
        **LLM Docstring**

        Re-express a force-field expansion (a list of derivative tensors) from one mode basis to another, by first rotating every non-skipped axis back through the old inverse-mode matrix (into internal coordinates) and optionally reordering, then rotating forward again through the new mode matrix.

        :param old_inv: the inverse mode matrix to rotate the old basis's derivative axes back through
        :type old_inv: np.ndarray
        :param new_modes: the (transposed) new mode matrix to rotate the intermediate (internal-coordinate) tensors forward through
        :type new_modes: np.ndarray
        :param old_field: the force-field expansion (list of derivative tensors) in the old mode basis
        :type old_field: list[np.ndarray]
        :param dim_skips: number of trailing axes to leave untouched by the rotation (e.g. a leading Cartesian-component axis)
        :type dim_skips: int
        :param sorting: a permutation to reorder the intermediate (internal-coordinate) tensor axes
        :type sorting: np.ndarray | None
        :return: `(new_field, mid_field)` -- the re-expressed force field in the new mode basis, and the intermediate internal-coordinate-basis force field
        :rtype: tuple[list[np.ndarray], list[np.ndarray]]
        """
        ...

    @classmethod
    def reexpress_normal_modes(cls, base_modes, old_field, dipole, sorting=None, type=2):
        """
        **LLM Docstring**

        Re-express a potential (and, optionally, dipole) expansion from one normal-mode basis into a renormalized one, via `renormalize_modes` followed by `rerotate_force_field`.

        :param base_modes: the `(freqs, modes, inv)` triple describing the original mode basis
        :type base_modes: tuple
        :param old_field: the potential-energy expansion (force constants) in the original mode basis
        :type old_field: list[np.ndarray]
        :param dipole: the dipole-derivative expansion (one list per Cartesian component) in the original mode basis, or `None`
        :type dipole: list[list[np.ndarray]] | None
        :param sorting: a permutation to reorder the intermediate internal-coordinate axes
        :type sorting: np.ndarray | None
        :param type: the mode-dimensionality convention forwarded to `renormalize_modes`
        :type type: int
        :return: `((freq, matrix, inv), potential_terms, dipole)` -- the renormalized mode data, the re-expressed potential terms, and the re-expressed dipole terms (or `None`)
        :rtype: tuple
        """
        ...
    convert = UnitsData.convert

    @staticmethod
    def mass(atom):
        """
        **LLM Docstring**

        Look up an atom's mass (in atomic units) from its element symbol.

        :param atom: the element symbol
        :type atom: str
        :return: the atom's mass, converted from atomic mass units to atomic units
        :rtype: float
        """
        ...

    @classmethod
    def extract_term_lists(cls, checkpoint, terms, skip_dimensions=0, threshold=0, aggregator=None):
        """
        **LLM Docstring**

        Read a set of stored expansion terms out of a checkpoint file and flatten each one into a list of `(index..., value)` rows -- restricted to non-redundant (sorted) index combinations and, optionally, thresholded to only significant values -- for later text-file export.

        :param checkpoint: the checkpoint file path to read from
        :type checkpoint: str
        :param terms: the checkpoint key naming the expansion to extract
        :type terms: str
        :param skip_dimensions: number of leading tensor axes to exclude from the symmetry/threshold filtering (e.g. a leading Cartesian-component axis)
        :type skip_dimensions: int
        :param threshold: minimum absolute value required for an entry to be included; `0` includes everything
        :type threshold: float
        :param aggregator: an optional function to transform the raw list of stored term tensors before flattening (e.g. to reassemble dipole components)
        :type aggregator: callable | None
        :return: a list, one entry per expansion order, of flattened `(idx1, idx2, ..., value)` row lists (1-indexed)
        :rtype: list[list[tuple]]
        """
        ...

    @classmethod
    def write_term_lists(cls, terms, file_template=None, int_fmt='{:>3.0f}', float_fmt='{:>16.8e}', index_function=None):
        """
        **LLM Docstring**

        Write a set of flattened term-list data (as produced by `extract_term_lists`) out to text files (or in-memory string buffers, if no file template is given), one file per expansion order.

        :param terms: the per-order flattened term-list data to write
        :type terms: list[list[tuple]]
        :param file_template: a filename template (with a single `{}` placeholder for the order index) to write each order's file to; if `None`, an in-memory `io.StringIO` is used instead
        :type file_template: str | None
        :param int_fmt: the format string used for integer-valued columns (indices)
        :type int_fmt: str
        :param float_fmt: the format string used for floating-point-valued columns
        :type float_fmt: str
        :param index_function: a function mapping the order index to the value substituted into `file_template`; identity by default
        :type index_function: callable | None
        :return: the list of written file paths (or in-memory `StringIO` buffers), one per order
        :rtype: list
        """
        ...

    @classmethod
    def extract_terms(cls, chk, out, terms, default_output='output.hdf5', aggregator=None, index_function=None, skip_dimensions=0):
        """
        **LLM Docstring**

        Extract and write out a named expansion's terms from a checkpoint (or a directory containing a default-named checkpoint file), via `extract_term_lists` and `write_term_lists`.

        :param chk: the checkpoint file path, or a directory containing `default_output`
        :type chk: str
        :param out: the output filename template
        :type out: str
        :param terms: the checkpoint key naming the expansion to extract
        :type terms: str
        :param default_output: the default checkpoint filename to look for within `chk` if it's a directory
        :type default_output: str
        :param aggregator: an optional function to transform the raw term data before flattening
        :type aggregator: callable | None
        :param index_function: a function mapping the order index to the value substituted into the output filename template
        :type index_function: callable | None
        :param skip_dimensions: number of leading tensor axes to exclude from symmetry filtering
        :type skip_dimensions: int
        :return: the list of written output files/buffers
        :rtype: list
        """
        ...

    @classmethod
    def extract_potential(cls, chk, out='potential_expansion_{}.dat'):
        """
        **LLM Docstring**

        Extract and write out the potential-energy expansion terms from a checkpoint, via `extract_terms`.

        :param chk: the checkpoint file path or directory
        :type chk: str
        :param out: the output filename template
        :type out: str
        :return: the list of written output files/buffers
        :rtype: list
        """
        ...

    @classmethod
    def extract_gmatrix(cls, chk, out='gmatrix_expansion_{}.dat'):
        """
        **LLM Docstring**

        Extract and write out the G-matrix expansion terms from a checkpoint, via `extract_terms`.

        :param chk: the checkpoint file path or directory
        :type chk: str
        :param out: the output filename template
        :type out: str
        :return: the list of written output files/buffers
        :rtype: list
        """
        ...

    @classmethod
    def extract_dipole_expansion(cls, chk, out='dipole_expansion_{}.dat'):
        """
        **LLM Docstring**

        Extract and write out the dipole-derivative expansion terms from a checkpoint, via `extract_terms`, first regrouping the per-axis (`x`/`y`/`z`) stored terms into single per-order tensors with a leading Cartesian-component axis.

        :param chk: the checkpoint file path or directory
        :type chk: str
        :param out: the output filename template
        :type out: str
        :return: the list of written output files/buffers
        :rtype: list
        """
        ...

    @classmethod
    def run_anne_job(cls, base_dir, states=2, initial_states=0, calculate_intensities=None, return_analyzer=False, return_runner=False, modes_file=('nm_int.dat', 'modes.dat'), atoms_file='atom.dat', masses_file='mass.dat', coords_file='cart_ref.dat', zmat_file='z_mat.dat', potential_files=('cub.dat', 'quart.dat', 'quintic.dat', 'sextic.dat'), dipole_files=('lin_dip.dat', 'quad_dip.dat', 'cub_dip.dat', 'quart_dip.dat', 'quintic_dip.dat'), coordinate_transformation=None, coordinate_transformation_file='coordinate_transformation.py', results_file=None, order=None, expansion_order=None, energy_units=None, normalization_type=0, **opts):
        """
        **LLM Docstring**

        Run a full VPT calculation from a directory of Anne-format input files (normal modes, atoms, masses, coordinates, Z-matrix, potential/dipole expansions, and an optional custom coordinate transformation), parsing and re-expressing the inputs as needed before dispatching to `VPTRunner.construct`/`run_simple` or `VPTAnalyzer.run_VPT`.

        :param base_dir: the directory containing the Anne-format input files (or `None`/`"."` to use the current directory)
        :type base_dir: str | None
        :param states: the target states (or quanta cutoff) to compute
        :type states: int | Iterable
        :param initial_states: the initial states to compute transitions from
        :type initial_states: int | Iterable
        :param calculate_intensities: whether to compute IR intensities; defaults to whether dipole data is available
        :type calculate_intensities: bool | None
        :param return_analyzer: whether to run via `VPTAnalyzer.run_VPT` and return an analyzer
        :type return_analyzer: bool
        :param return_runner: whether to run via `VPTRunner.construct` and return the constructed runner (without running)
        :type return_runner: bool
        :param modes_file: the normal-modes input file(s), tried in order
        :type modes_file: str | tuple
        :param atoms_file: the atoms input file
        :type atoms_file: str
        :param masses_file: the masses input file
        :type masses_file: str
        :param coords_file: the Cartesian-coordinates input file
        :type coords_file: str
        :param zmat_file: the Z-matrix input file
        :type zmat_file: str | None
        :param potential_files: the potential-derivative input files
        :type potential_files: tuple
        :param dipole_files: the dipole-derivative input files
        :type dipole_files: tuple
        :param coordinate_transformation: an explicit `[conversion, inverse]` coordinate-transformation function pair
        :type coordinate_transformation: list | None
        :param coordinate_transformation_file: a Python module file defining `conversion`/`inverse` functions, used if `coordinate_transformation` isn't given
        :type coordinate_transformation_file: str
        :param results_file: where to store the results checkpoint
        :type results_file: str | None
        :param order: the perturbation-theory order
        :type order: int | None
        :param expansion_order: the per-term expansion orders
        :type expansion_order: int | dict | None
        :param energy_units: the units the input energies/frequencies are given in; inferred from magnitude if not given
        :type energy_units: str | None
        :param normalization_type: which mode-renormalization convention to use when re-expressing into internal coordinates
        :type normalization_type: int
        :param opts: extra options forwarded to the underlying runner/analyzer
        :type opts: dict
        :return: the run results (wavefunctions, analyzer, or runner, depending on `return_analyzer`/`return_runner`)
        :rtype: object
        :raises ValueError: if `coordinate_transformation` is given without both a conversion and inverse function
        """
        ...

    @classmethod
    def run_fchk_job(cls, base_dir, states=2, calculate_intensities=None, return_analyzer=False, return_runner=False, zmat_file='z_mat.dat', fchk_file='fchk.fchk', results_file='output.hdf5', **opts):
        """
        **LLM Docstring**

        Run a VPT calculation using a Gaussian FChk file for the molecule/potential/dipole data, together with an optional Z-matrix file for the internal-coordinate specification, dispatching to `VPTRunner.construct`/`run_simple` or `VPTAnalyzer.run_VPT`.

        :param base_dir: the directory containing the FChk (and optional Z-matrix) file
        :type base_dir: str
        :param states: the target states (or quanta cutoff) to compute
        :type states: int | Iterable
        :param calculate_intensities: whether to compute IR intensities
        :type calculate_intensities: bool | None
        :param return_analyzer: whether to run via `VPTAnalyzer.run_VPT` and return an analyzer
        :type return_analyzer: bool
        :param return_runner: whether to run via `VPTRunner.construct` and return the constructed runner (without running)
        :type return_runner: bool
        :param zmat_file: the Z-matrix input file; if missing, Cartesian coordinates are used directly
        :type zmat_file: str | None
        :param fchk_file: the Gaussian FChk file to load the molecule/potential/dipole data from
        :type fchk_file: str
        :param results_file: where to store the results checkpoint
        :type results_file: str
        :param opts: extra options forwarded to the underlying runner/analyzer
        :type opts: dict
        :return: the run results (wavefunctions, analyzer, or runner, depending on `return_analyzer`/`return_runner`)
        :rtype: object
        """
        ...

    @classmethod
    def get_internal_expansion(cls, fchk, internals, states=2, **opts):
        """
        **LLM Docstring**

        Build a throwaway `VPTRunner` for a given FChk file and internal-coordinate specification, purely to extract the resulting internal-coordinate potential/kinetic expansion and mode transformation data (without running the full perturbation theory).

        :param fchk: the Gaussian FChk file to load
        :type fchk: str
        :param internals: the internal-coordinate (Z-matrix) specification to use
        :type internals: object
        :param states: the (nominal) target states used to construct the throwaway runner
        :type states: int | Iterable
        :param opts: extra options forwarded to `VPTRunner.construct`
        :type opts: dict
        :return: a dict with `'runner'`, `'freqs'`, `'molecule'`, `'kinetic'`, `'potential'`, `'modes'`, `'states'`, `'fchk'`, and `'zmatrix'` entries describing the extracted internal-coordinate expansion
        :rtype: dict
        """
        ...

    @classmethod
    def run_internal_expansion(cls, expansion_data, calculate_intensities=False, **opts):
        """
        **LLM Docstring**

        Run a VPT calculation directly from a previously extracted internal-coordinate expansion (as produced by `get_internal_expansion`), reusing its modes and truncated (harmonic/cubic/quartic) potential expansion.

        :param expansion_data: the expansion data dict, as returned by `get_internal_expansion`
        :type expansion_data: dict
        :param calculate_intensities: whether to compute IR intensities
        :type calculate_intensities: bool
        :param opts: accepted for interface consistency but not used in this method's body
        :type opts: dict
        :return: the computed VPT wavefunctions
        :rtype: PerturbationTheoryWavefunctions
        """
        ...
VPTRunner.helpers = AnneInputHelpers

class MultiVPTStateSpace:
    """
    Generalizes a VPTStateSpace to pairs of initial and final spaces
    """

    def __init__(self, state_space_pairs, system=None, degeneracy_specs=None, evaluator=None, **opts):
        """
        **LLM Docstring**

        Build a collection of `(initial, target)` `VPTStateSpace` pairs sharing a common system, resolving each pair's spec, merging any per-pair degenerate-state groupings that overlap across pairs into a single consistent set (re-including the merged groups' states into every affected pair), and computing a flattened union of all states/degenerate-state pairs across the whole collection.

        :param state_space_pairs: the `(initial, target)` state-space specs, or a bare single target-state spec (paired with the ground state `0`)
        :type state_space_pairs: object | list[list]
        :param system: the system (molecule/mode context) the states are defined relative to
        :type system: VPTSystem | None
        :param degeneracy_specs: the degeneracy specification(s) to apply when resolving each `VPTStateSpace`
        :type degeneracy_specs: object | None
        :param evaluator: an evaluator to forward to `VPTStateSpace.from_system_and_spec` for degeneracy handling
        :type evaluator: object | None
        :param opts: extra options forwarded to `VPTStateSpace.from_system_and_spec`
        :type opts: dict
        :return: None
        :rtype: None
        """
        ...

    @property
    def state_list_pairs(self):
        """
        **LLM Docstring**

        The `[initial_state_list, target_state_list]` pairs for every `(initial, target)` space pair in this collection.

        :return: the list of state-list pairs
        :rtype: list[list]
        """
        ...

class AnalyticVPTRunner:

    def __init__(self, expansions, order=None, expansion_order=None, freqs=None, internals=True, logger=None, hamiltonian=None, checkpoint=None, dipole_expansion=None, allowed_terms=None, allowed_coefficients=None, disallowed_coefficients=None, allowed_energy_changes=None, intermediate_normalization=None, local_mode_couplings=None, local_mode_coupling_order=None, parallelizer=None):
        """
        **LLM Docstring**

        Set up an analytic (symbolic) VPT evaluator, either wrapping an already-built `PerturbationTheoryEvaluator` directly or constructing one from a raw set of Taylor-expansion coefficients and an `AnalyticPerturbationTheorySolver` built for the requested order, and resolving the dipole expansion (from the attached classic Hamiltonian, if any and not given directly) needed for transition-moment/intensity calculations.

        :param expansions: the raw per-order `[V, G, ...]` expansion coefficients, or an already-built `PerturbationTheoryEvaluator`
        :type expansions: list | PerturbationTheoryEvaluator
        :param order: the perturbation-theory order to build the solver for
        :type order: int | None
        :param expansion_order: the per-term expansion orders used to resolve the dipole expansion's default order
        :type expansion_order: dict | None
        :param freqs: the mode frequencies
        :type freqs: np.ndarray | None
        :param internals: whether the expansion is expressed in internal coordinates
        :type internals: bool
        :param logger: logger for diagnostics
        :type logger: Logger | None
        :param hamiltonian: an associated classic `PerturbationTheoryHamiltonian`, used as a fallback source for the dipole expansion and other molecule-specific data
        :type hamiltonian: PerturbationTheoryHamiltonian | None
        :param checkpoint: checkpoint file/object for caching intermediate symbolic expressions
        :type checkpoint: str | Checkpointer | None
        :param dipole_expansion: explicit dipole-derivative expansion terms to use instead of deriving them from `hamiltonian`
        :type dipole_expansion: list[np.ndarray] | None
        :param allowed_terms: restrict the analytic solver to only these perturbation term types
        :type allowed_terms: object | None
        :param allowed_coefficients: restrict the analytic solver to only these expansion coefficients
        :type allowed_coefficients: object | None
        :param disallowed_coefficients: exclude these expansion coefficients from the analytic solver
        :type disallowed_coefficients: object | None
        :param allowed_energy_changes: restrict the analytic solver to only these energy-change patterns
        :type allowed_energy_changes: object | None
        :param intermediate_normalization: whether to use intermediate (rather than full) normalization in the analytic solver
        :type intermediate_normalization: bool | None
        :param local_mode_couplings: extra local-mode coupling terms to inject into the Hamiltonian when running VPT
        :type local_mode_couplings: list | None
        :param local_mode_coupling_order: the perturbative order local-mode couplings should be injected at
        :type local_mode_coupling_order: int | None
        :param parallelizer: parallelization backend for the evaluator
        :type parallelizer: Parallelizer | None
        :return: None
        :rtype: None
        """
        ...

    @classmethod
    def from_hamiltonian(cls, ham, order, expansion_order=None, logger=None, checkpoint=None, parallelizer=None, allowed_terms=None, allowed_coefficients=None, disallowed_coefficients=None, allowed_energy_changes=None, take_diagonal_v4_terms=True, intermediate_normalization=None, corrected_fundamental_frequencies=None, **opts):
        """
        A driver powered by a classic PerturbationTheoryHamiltonian object

        :param ham:
        :return:
        """
        ...

    @classmethod
    def construct(cls, system, states=None, *, order=2, expressions_file=None, allowed_terms=None, allowed_coefficients=None, disallowed_coefficients=None, allowed_energy_changes=None, mixed_derivative_handling_mode='analytical', degeneracy_specs=None, corrected_fundamental_frequencies=None, parallelizer=None, **settings) -> '(AnalyticVPTRunner, VPTMultiStateSpace)':
        """
            **LLM Docstring**

            Build an `AnalyticVPTRunner` (and, if target states are given, its resolved `MultiVPTStateSpace`) directly from a molecule/system specification: builds a throwaway single-state `VPTRunner` to assemble the classic Hamiltonian, then derives the analytic evaluator from it via `from_hamiltonian`.

            :param system: the molecule or system specification to run VPT on
            :type system: str | list | Molecule | VPTSystem
            :param states: the target states (possibly as `(initial, target)` pairs)
            :type states: object | None
            :param order: the perturbation-theory order
            :type order: int
            :param expressions_file: checkpoint file for caching the analytic expressions
            :type expressions_file: str | None
            :param allowed_terms: restrict the analytic solver to only these perturbation term types
            :type allowed_terms: object | None
            :param allowed_coefficients: restrict the analytic solver to only these expansion coefficients
            :type allowed_coefficients: object | None
            :param disallowed_coefficients: exclude these expansion coefficients from the analytic solver
            :type disallowed_coefficients: object | None
            :param allowed_energy_changes: restrict the analytic solver to only these energy-change patterns
            :type allowed_energy_changes: object | None
            :param mixed_derivative_handling_mode: the mixed-derivative handling mode to use when building the classic Hamiltonian
            :type mixed_derivative_handling_mode: str
            :param degeneracy_specs: the degeneracy specification(s) to build the `MultiVPTStateSpace` with
            :type degeneracy_specs: object | None
            :param corrected_fundamental_frequencies: replacement fundamental frequencies to use
            :type corrected_fundamental_frequencies: np.ndarray | None
            :param parallelizer: parallelization backend for the evaluator
            :type parallelizer: Parallelizer | None
            :param settings: extra options forwarded to `VPTRunner.construct`
            :type settings: dict
            :return: the constructed runner, or `(runner, states)` if `states` was given
            :rtype: AnalyticVPTRunner | tuple[AnalyticVPTRunner, MultiVPTStateSpace]
            """
        ...

    @classmethod
    def from_file(cls, file_name, order=2, allowed_terms=None, allowed_coefficients=None, disallowed_coefficients=None, allowed_energy_changes=None, expressions_file=None, **settings):
        """
        **LLM Docstring**

        Build an `AnalyticVPTRunner` from a molecule file (e.g. an FChk), via a throwaway single-state `VPTRunner.construct` followed by `from_hamiltonian`.

        :param file_name: the molecule file to load
        :type file_name: str
        :param order: the perturbation-theory order
        :type order: int
        :param allowed_terms: restrict the analytic solver to only these perturbation term types
        :type allowed_terms: object | None
        :param allowed_coefficients: restrict the analytic solver to only these expansion coefficients
        :type allowed_coefficients: object | None
        :param disallowed_coefficients: exclude these expansion coefficients from the analytic solver
        :type disallowed_coefficients: object | None
        :param allowed_energy_changes: restrict the analytic solver to only these energy-change patterns
        :type allowed_energy_changes: object | None
        :param expressions_file: checkpoint file for caching the analytic expressions
        :type expressions_file: str | None
        :param settings: extra options forwarded to `VPTRunner.construct`
        :type settings: dict
        :return: the constructed runner
        :rtype: AnalyticVPTRunner
        """
        ...

    def construct_classic_runner(self, states, system=None, logger=None, corrected_fundamental_frequencies=None, potential_terms=None, kinetic_terms=None, coriolis_terms=None, pseudopotential_terms=None, dipole_terms=None, initial_states=None, **opts):
        """
        **LLM Docstring**

        Build a classic `VPTRunner` reproducing this analytic evaluator's expansion data (rescaling the stored symbolic-coefficient conventions back into ordinary Taylor-expansion derivative tensors for potential/kinetic/Coriolis/pseudopotential/dipole terms), for cross-validation against the classic (matrix) VPT solver, or for use where a classic runner interface is required.

        :param states: the target states, or a `MultiVPTStateSpace` (in which case its flattened state space and per-pair initial states are used)
        :type states: object | MultiVPTStateSpace
        :param system: the molecule/system to attach; derived from the associated Hamiltonian, or built as a dummy placeholder molecule, if not given
        :type system: Molecule | None
        :param logger: logger for the constructed runner; defaults to the associated Hamiltonian's logger
        :type logger: Logger | None
        :param corrected_fundamental_frequencies: replacement fundamental frequencies; defaults to this evaluator's own frequencies
        :type corrected_fundamental_frequencies: np.ndarray | None
        :param potential_terms: explicit potential-expansion tensors, bypassing the default rescaling
        :type potential_terms: list[np.ndarray] | None
        :param kinetic_terms: explicit kinetic-expansion tensors, bypassing the default rescaling
        :type kinetic_terms: list[np.ndarray] | None
        :param coriolis_terms: explicit Coriolis-expansion tensors, bypassing the default rescaling
        :type coriolis_terms: list[np.ndarray] | None
        :param pseudopotential_terms: explicit pseudopotential-expansion tensors, bypassing the default rescaling
        :type pseudopotential_terms: list[np.ndarray] | None
        :param dipole_terms: explicit dipole-expansion tensors; defaults to this evaluator's own dipole expansion
        :type dipole_terms: list[np.ndarray] | None
        :param initial_states: the initial states for the classic runner; derived from `states` if it's a `MultiVPTStateSpace`
        :type initial_states: object | None
        :param opts: extra options forwarded to `VPTRunner.construct`
        :type opts: dict
        :return: `(runner, states)`, the constructed classic runner and its resolved state space
        :rtype: tuple
        """
        ...

    @classmethod
    def clear_caches(cls):
        """
        **LLM Docstring**

        Clear the global caches used by the underlying `AnalyticPerturbationTheorySolver` (e.g. cached symbolic expressions shared across instances).

        :return: None
        :rtype: None
        """
        ...

    @classmethod
    def prep_multispace(self, states, freqs, system=None, degeneracy_specs=None):
        """
        **LLM Docstring**

        Coerce a raw state specification into a `MultiVPTStateSpace`, passing an already-built one through unchanged and normalizing a bare `[[raise, lower], ...]` polyad-pair specification into the `{'polyads': ...}` dict form expected by `DegeneracySpec`.

        :param states: the raw state specification, or an already-built `MultiVPTStateSpace`
        :type states: object | MultiVPTStateSpace
        :param freqs: the mode frequencies used to build the underlying `VPTStateSpace`(s)
        :type freqs: np.ndarray
        :param system: the system (molecule/mode context) to associate with the state space
        :type system: object | None
        :param degeneracy_specs: the degeneracy specification(s) to apply
        :type degeneracy_specs: object | None
        :return: the resolved `MultiVPTStateSpace`
        :rtype: MultiVPTStateSpace
        """
        ...

    class _dummy_system:

        def __init__(self, runner):
            """
            **LLM Docstring**

            Wrap this `AnalyticVPTRunner` so it exposes just enough of the `VPTSystem` interface (specifically, `nmodes`) to be usable in places that expect a full system object but only actually need the mode count.

            :param runner: the `AnalyticVPTRunner` to wrap
            :type runner: AnalyticVPTRunner
            :return: None
            :rtype: None
            """
            ...

    def prep_states(self, states, degeneracy_specs=None):
        """
        **LLM Docstring**

        Coerce a raw state specification into a `MultiVPTStateSpace` using this evaluator's own frequencies and a lightweight dummy system object, via `prep_multispace`.

        :param states: the raw state specification
        :type states: object
        :param degeneracy_specs: the degeneracy specification(s) to apply
        :type degeneracy_specs: object | None
        :return: the resolved `MultiVPTStateSpace`
        :rtype: MultiVPTStateSpace
        """
        ...

    def evaluate_expressions(self, states, exprs, zero_cutoff=None, operator_expansions=None, degeneracy_specs=None, verbose=False):
        """
        **LLM Docstring**

        Evaluate a set of arbitrary symbolic perturbation-theory expressions numerically for the given target states, via the underlying `PerturbationTheoryEvaluator`.

        :param states: the target states to evaluate the expressions for
        :type states: object
        :param exprs: the symbolic expression(s) to evaluate
        :type exprs: object
        :param zero_cutoff: the magnitude below which a term is treated as exactly zero
        :type zero_cutoff: float | None
        :param operator_expansions: extra operator expansions the expressions may reference
        :type operator_expansions: object | None
        :param degeneracy_specs: the degeneracy specification(s) to apply when resolving the state space
        :type degeneracy_specs: object | None
        :param verbose: whether to log detailed evaluation progress
        :type verbose: bool
        :return: the evaluated expression results
        :rtype: object
        """
        ...

    def get_matrix_corrections(self, states, order=None, degeneracy_specs=None, zero_cutoff=None, verbose=False):
        """
        **LLM Docstring**

        Compute the perturbative matrix-element corrections for the given target states, via the underlying `PerturbationTheoryEvaluator`.

        :param states: the target states to compute corrections for
        :type states: object
        :param order: the perturbation-theory order to compute to
        :type order: int | None
        :param degeneracy_specs: the degeneracy specification(s) to apply
        :type degeneracy_specs: object | None
        :param zero_cutoff: the magnitude below which a term is treated as exactly zero
        :type zero_cutoff: float | None
        :param verbose: whether to log detailed evaluation progress
        :type verbose: bool
        :return: the matrix-element corrections
        :rtype: object
        """
        ...

    def get_energy_corrections(self, states, order=None, degeneracy_specs=None, zero_cutoff=None, verbose=False):
        """
        **LLM Docstring**

        Compute the perturbative energy corrections for the given target states, via the underlying `PerturbationTheoryEvaluator`.

        :param states: the target states to compute energy corrections for
        :type states: object
        :param order: the perturbation-theory order to compute to
        :type order: int | None
        :param degeneracy_specs: the degeneracy specification(s) to apply
        :type degeneracy_specs: object | None
        :param zero_cutoff: the magnitude below which a term is treated as exactly zero
        :type zero_cutoff: float | None
        :param verbose: whether to log detailed evaluation progress
        :type verbose: bool
        :return: the per-order energy corrections
        :rtype: np.ndarray
        """
        ...

    def get_overlap_corrections(self, states, order=None, degeneracy_specs=None, zero_cutoff=None, verbose=False):
        """
        **LLM Docstring**

        Compute the perturbative wavefunction-overlap corrections for the given target states, via the underlying `PerturbationTheoryEvaluator`.

        :param states: the target states to compute overlap corrections for
        :type states: object
        :param order: the perturbation-theory order to compute to
        :type order: int | None
        :param degeneracy_specs: the degeneracy specification(s) to apply
        :type degeneracy_specs: object | None
        :param zero_cutoff: the magnitude below which a term is treated as exactly zero
        :type zero_cutoff: float | None
        :param verbose: whether to log detailed evaluation progress
        :type verbose: bool
        :return: the overlap corrections
        :rtype: object
        """
        ...

    @classmethod
    def prep_eval_state_pairs(cls, states):
        """
        **LLM Docstring**

        Flatten a `MultiVPTStateSpace`'s `(initial, final)` state-list pairs into a flat list of `[single_initial_state, final_states]` pairs, one per individual initial state (rather than grouped by pair-block), for direct consumption by the underlying evaluator.

        :param states: the multi-state space to flatten
        :type states: MultiVPTStateSpace
        :return: the flattened `[initial_state, final_states]` pairs
        :rtype: list[list]
        """
        ...

    def get_full_wavefunction_corrections(self, states, order=None, degeneracy_specs=None, zero_cutoff=None, verbose=False):
        """
        **LLM Docstring**

        Compute the full (all-component) perturbative wavefunction corrections for the given target states, via the underlying `PerturbationTheoryEvaluator`.

        :param states: the target states to compute wavefunction corrections for
        :type states: object
        :param order: the perturbation-theory order to compute to
        :type order: int | None
        :param degeneracy_specs: the degeneracy specification(s) to apply
        :type degeneracy_specs: object | None
        :param zero_cutoff: the magnitude below which a term is treated as exactly zero
        :type zero_cutoff: float | None
        :param verbose: whether to log detailed evaluation progress
        :type verbose: bool
        :return: the full wavefunction corrections
        :rtype: object
        """
        ...

    def get_wavefunction_corrections(self, states, order=None, degeneracy_specs=None, zero_cutoff=None, verbose=False):
        """
        **LLM Docstring**

        Compute the perturbative wavefunction corrections for the given target states, via the underlying `PerturbationTheoryEvaluator`.

        :param states: the target states to compute wavefunction corrections for
        :type states: object
        :param order: the perturbation-theory order to compute to
        :type order: int | None
        :param degeneracy_specs: the degeneracy specification(s) to apply
        :type degeneracy_specs: object | None
        :param zero_cutoff: the magnitude below which a term is treated as exactly zero
        :type zero_cutoff: float | None
        :param verbose: whether to log detailed evaluation progress
        :type verbose: bool
        :return: the wavefunction corrections
        :rtype: object
        """
        ...

    @classmethod
    def unflatten_corr(cls, states, corrs):
        """
        **LLM Docstring**

        Regroup a flat correction result (expressed over the combined initial/final state spaces) back into the per-`(initial, final)`-block structure implied by a `MultiVPTStateSpace`'s state-list pairs.

        :param states: the multi-state space whose block structure the flat corrections should be regrouped into
        :type states: MultiVPTStateSpace
        :param corrs: the flat correction data (with `initial_states`/`final_states`/`corrections` attributes) to regroup
        :type corrs: object
        :return: the per-block regrouped correction arrays
        :rtype: list[list[np.ndarray]]
        """
        ...

    def get_operator_corrections(self, operator_expansion, states, order=None, terms=None, degeneracy_specs=None, verbose=False, operator_type=None, check_single=True, **opts):
        """
        **LLM Docstring**

        Compute the perturbative corrections to one or more arbitrary operator expansions for the given target states, via the underlying `PerturbationTheoryEvaluator`, then regroup the flat results back into per-block structure via `unflatten_corr`.

        :param operator_expansion: the raw operator expansion(s) to evaluate, each a list of derivative tensors
        :type operator_expansion: list
        :param states: the target states to compute corrections for
        :type states: object
        :param order: the perturbation-theory order to compute to
        :type order: int | None
        :param terms: restrict the evaluation to specific term contributions
        :type terms: object | None
        :param degeneracy_specs: the degeneracy specification(s) to apply
        :type degeneracy_specs: object | None
        :param verbose: whether to log detailed evaluation progress
        :type verbose: bool
        :param operator_type: the operator's symmetry/rank type; inferred from the expansion tensor shapes if not given
        :type operator_type: str | int | None
        :param check_single: whether to auto-detect and wrap a single (rather than multiple) operator expansion
        :type check_single: bool
        :param opts: extra options forwarded to the underlying evaluator
        :type opts: dict
        :return: the per-block regrouped operator corrections
        :rtype: list
        """
        ...

    def construct_corrections_vectors(self, states, corrs):
        """
        **LLM Docstring**

        Assemble a set of flat per-order correction matrices spanning the full flat state space, from one or more raw per-block correction results, indexing rows by the union of every block's initial states.

        :param states: the state space whose flat state list defines the column indexing
        :type states: object
        :param corrs: one (or a list of) raw correction result(s), each with `initial_states`/`final_states`/`corrections`
        :type corrs: object | list
        :return: `(init_states, op_corr_mats)` -- the union of initial states used for row indexing, and the assembled per-order correction matrices (single set if `corrs` was a single result)
        :rtype: tuple
        """
        ...

    def construct_corrections_matrix(self, group, corrs):
        """
        **LLM Docstring**

        Assemble a set of square per-order correction matrices restricted to a single group of states, from one or more raw per-block correction results.

        :param group: the group of states (typically a degenerate block) to build the matrix over
        :type group: object
        :param corrs: one (or a list of) raw correction result(s), each with `initial_states`/`final_states`/`corrections`
        :type corrs: object | list
        :return: the assembled per-order correction matrices, restricted to `group` (single set if `corrs` was a single result)
        :rtype: list
        """
        ...

    def get_transition_moment_corrections(self, states, dipole_expansion=None, order=None, degeneracy_specs=None, axes=None, **opts):
        """
        **LLM Docstring**

        Compute the perturbative transition-dipole-moment corrections for the given target states, using this evaluator's dipole expansion (or an explicitly given one), scaled by the appropriate factorial normalization, via `get_operator_corrections`.

        :param states: the target states to compute transition moments for
        :type states: object
        :param dipole_expansion: an explicit dipole-derivative expansion to use instead of `self.dipole_expansion`
        :type dipole_expansion: list[np.ndarray] | None
        :param order: the perturbation-theory order to compute to; inferred from `self.expansion_order`/the dipole expansion length if not given
        :type order: int | None
        :param degeneracy_specs: the degeneracy specification(s) to apply
        :type degeneracy_specs: object | None
        :param axes: which Cartesian axes to compute transition moments for; defaults to all three
        :type axes: Iterable[int] | None
        :param opts: extra options forwarded to `get_operator_corrections`
        :type opts: dict
        :return: the per-block regrouped transition-moment corrections
        :rtype: list
        """
        ...

    def get_freqs(self, states, order=None, degeneracy_specs=None, return_corrections=False, verbose=False):
        """
        **LLM Docstring**

        Compute the vibrational transition frequencies (in wavenumbers) for the given target states, relative to the ground/reference state.

        :param states: the target states to compute frequencies for
        :type states: object
        :param order: the perturbation-theory order to compute to
        :type order: int | None
        :param degeneracy_specs: the degeneracy specification(s) to apply
        :type degeneracy_specs: object | None
        :param return_corrections: whether to also return the raw per-order energy corrections
        :type return_corrections: bool
        :param verbose: whether to log detailed evaluation progress
        :type verbose: bool
        :return: `(zpe, freqs)`, or `((zpe, freqs), corrections)` if `return_corrections` is set
        :rtype: tuple
        """
        ...

    def get_reexpressed_hamiltonian(self, states, order=None, degeneracy_specs=None, only_degenerate_terms=True, verbose=False, hamiltonian_corrections=None, **opts):
        """
        **LLM Docstring**

        Build the deperturbed (degenerate-block) Hamiltonian matrices for each degenerate group of states, re-expressing the analytic Hamiltonian corrections in the basis of each group and summing them into square correction matrices.

        :param states: the target states (used to resolve degenerate groupings)
        :type states: object
        :param order: the perturbation-theory order to compute to
        :type order: int | None
        :param degeneracy_specs: the degeneracy specification(s) to apply
        :type degeneracy_specs: object | None
        :param only_degenerate_terms: whether to include only the strictly degenerate-coupling terms in the reexpressed Hamiltonian
        :type only_degenerate_terms: bool
        :param verbose: whether to log detailed evaluation progress
        :type verbose: bool
        :param hamiltonian_corrections: extra Hamiltonian corrections to fold in
        :type hamiltonian_corrections: object | None
        :param opts: extra options forwarded to the underlying evaluator
        :type opts: dict
        :return: `(all_mats, all_corrs)` -- the summed Hamiltonian matrix and the individual per-order correction matrices, one entry per degenerate group; `None` if there are no degenerate states
        :rtype: tuple | None
        """
        ...

    def get_wfc_test_states(self, input_states: BasisStateSpace, energy_window):
        """
        **LLM Docstring**

        Identify the candidate states that could plausibly be strongly coupled (via wavefunction-correction magnitude) to a given set of input states, based purely on energy proximity (within `energy_window`) and the maximum possible change in quantum numbers implied by the Hamiltonian expansion's term structure.

        :param input_states: the states to find coupling candidates for
        :type input_states: BasisStateSpace
        :param energy_window: the energy window (around each input state's own energy) to search within
        :type energy_window: float
        :return: the candidate coupling states for each input state
        :rtype: list[BasisStateSpace]
        """
        ...

    def get_test_wfn_corrs(self, input_states: BasisStateSpace, energy_window):
        """
        We take the expansions and frequencies that we have and at find the possible terms
        that could possibly lead to a correction greater than the specified threshold
        To do this, we first determine from the expansions what magnitude of energy difference
        could possible lead to terms above this threshold
        """
        ...

    def format_energies_table(self, states, energies, energy_corrections, zpe_pos, number_format='.3f'):
        """
        **LLM Docstring**

        Format a table of state energies/frequencies alongside their per-order corrections, converting to wavenumbers and displaying every non-zero-point-energy state as a frequency shift relative to the ZPE.

        :param states: the states the energies correspond to
        :type states: MultiVPTStateSpace
        :param energies: the total energies (in Hartrees), one per state
        :type energies: np.ndarray
        :param energy_corrections: the per-order energy corrections (in Hartrees)
        :type energy_corrections: np.ndarray
        :param zpe_pos: the index of the zero-point-energy (reference) state
        :type zpe_pos: int
        :param number_format: the format spec used for each numeric column
        :type number_format: str
        :return: the formatted table
        :rtype: str
        """
        ...

    def format_degenerate_energies_table(self, states, energies, deperturbed_energies, zpe_pos, number_format='.3f'):
        """
        **LLM Docstring**

        Format a table comparing each state's degenerate-perturbation-theory-corrected energy/frequency against its deperturbed counterpart, both relative to the zero-point energy.

        :param states: the states the energies correspond to
        :type states: MultiVPTStateSpace
        :param energies: the degenerate-corrected total energies (in Hartrees)
        :type energies: np.ndarray
        :param deperturbed_energies: the deperturbed total energies (in Hartrees)
        :type deperturbed_energies: np.ndarray
        :param zpe_pos: the index of the zero-point-energy (reference) state
        :type zpe_pos: int
        :param number_format: the format spec used for each numeric column
        :type number_format: str
        :return: the formatted table
        :rtype: str
        """
        ...

    def format_transition_moment_table(self, states, transition_moments, transition_moment_corrections, number_format='.8f'):
        """
        **LLM Docstring**

        Format a table of transition-dipole moments (and their per-order corrections) for each initial-state block in a `MultiVPTStateSpace`, with a labeled header separating each initial state's sub-table when there's more than one.

        :param states: the initial/final state pairs the transition moments correspond to
        :type states: MultiVPTStateSpace
        :param transition_moments: the per-axis, per-block final transition moments
        :type transition_moments: list
        :param transition_moment_corrections: the per-axis, per-block, per-order transition-moment corrections
        :type transition_moment_corrections: list
        :param number_format: the format spec used for each numeric column
        :type number_format: str
        :return: the formatted table(s), concatenated with initial-state header separators
        :rtype: str
        """
        ...

    def format_operators_table(self, states, keys, operator_values, operator_corrections, number_format='.8f'):
        """
        **LLM Docstring**

        Format a table of arbitrary operator expectation values (and their per-order corrections) for each initial-state block in a `MultiVPTStateSpace`, with a labeled header separating each initial state's sub-table when there's more than one.

        :param states: the initial/final state pairs the operator values correspond to
        :type states: MultiVPTStateSpace
        :param keys: the operator names/labels, used as column headers
        :type keys: Iterable
        :param operator_values: the per-operator, per-block final operator values
        :type operator_values: list
        :param operator_corrections: the per-operator, per-block, per-order operator corrections
        :type operator_corrections: list
        :param number_format: the format spec used for each numeric column
        :type number_format: str
        :return: the formatted table(s), concatenated with initial-state header separators
        :rtype: str
        """
        ...

    def format_spectrum_table(self, states, harmonic_spectra, spectra, deperturbed_spectra=None, number_format='.3f'):
        """
        **LLM Docstring**

        Format a table of harmonic, anharmonic, and (optionally) deperturbed IR spectra for each initial-state block in a `MultiVPTStateSpace`, with a labeled header separating each initial state's sub-table when there's more than one.

        :param states: the initial/final state pairs the spectra correspond to
        :type states: MultiVPTStateSpace
        :param harmonic_spectra: the per-block harmonic `DiscreteSpectrum` objects
        :type harmonic_spectra: list
        :param spectra: the per-block anharmonic `DiscreteSpectrum` objects
        :type spectra: list
        :param deperturbed_spectra: the per-block deperturbed `DiscreteSpectrum` objects, if degenerate treatment was applied
        :type deperturbed_spectra: list | None
        :param number_format: the format spec used for each numeric column
        :type number_format: str
        :return: the formatted table(s), concatenated with initial-state header separators
        :rtype: str
        """
        ...

    def prep_operators(self, operator_expansions, operator_terms, order=None):
        """
        **LLM Docstring**

        Normalize a user-supplied operator specification (raw expansion coefficients, in either list or named-dict form) into the fully-expanded per-mode-basis operator terms needed for correction calculations, using `self.ham.prep_operator_terms` to expand raw finite coefficients if pre-expanded `operator_terms` weren't already given directly.

        :param operator_expansions: raw operator expansion coefficients (a single operator's coefficient list, a list of such lists, or a name-keyed dict of them), used if `operator_terms` isn't given
        :type operator_expansions: object | None
        :param operator_terms: already-expanded per-mode-basis operator terms (a single list, a list of lists, or a name-keyed dict), used directly if given
        :type operator_terms: object | None
        :param order: the expansion order to expand `operator_expansions` to; inferred from the coefficient list length if not given
        :type order: int | None
        :return: `(keys, operator_terms)` -- the resolved operator names/indices and their expanded per-mode-basis terms
        :rtype: tuple
        """
        ...
    matrix_formatting_options = dict(linewidth=100000000.0, threshold=100000000.0, suppress=True, precision=3)

    def format_matrix(self, ham):
        """
        **LLM Docstring**

        Format a matrix as a plain-text string using this class's standard print options (fixed precision, no truncation, suppressed scientific notation), stripping the outer bracket characters for a cleaner look.

        :param ham: the matrix to format
        :type ham: np.ndarray | object
        :return: the formatted matrix string
        :rtype: str
        """
        ...

    def modify_hamiltonian(self, hamiltonian_corrections):
        """
        **LLM Docstring**

        Build a new `AnalyticVPTRunner` whose underlying evaluator has extra Hamiltonian corrections applied, via `PerturbationTheoryEvaluator.modify_hamiltonian`, preserving this runner's order/expansion-order/dipole-expansion/logger settings.

        :param hamiltonian_corrections: the Hamiltonian corrections to apply (e.g. local-mode coupling terms)
        :type hamiltonian_corrections: object
        :return: the new runner with the modified Hamiltonian
        :rtype: AnalyticVPTRunner
        """
        ...

    @classmethod
    def _prep_deg_pair_msg(cls, pairs, max_pairs=1000, fmt='{l} <-> {r}'):
        """
        **LLM Docstring**

        Format a list of degenerate state pairs as `"left <-> right"` log lines, truncating the middle of very long lists (showing only the first and last halves, separated by an ellipsis) to keep the log readable.

        :param pairs: the `(left, right)` state pairs to format
        :type pairs: list[tuple]
        :param max_pairs: the maximum number of pairs to display before truncating
        :type max_pairs: int
        :param fmt: the format string used for each pair, given `l`/`r` keyword values
        :type fmt: str
        :return: the formatted log lines
        :rtype: list[str]
        """
        ...
    hamiltonian_correction_modification_type = 'degenerate'

    def run_VPT(self, states, calculate_intensities=True, operator_expansions=None, operator_terms=None, operator_type=None, order=None, verbose=False, degeneracy_specs=None, handle_degeneracies=True, zero_cutoff=None, transition_moment_terms=None, hamiltonian_corrections=None, clear_caches=True, hamiltonian_correction_type=None, only_degenerate_terms=True, force_return_on_crash=True):
        """
        **LLM Docstring**

        Top-level entry point for running a full analytic VPT calculation on a set of target states: optionally injects local-mode coupling terms (and/or applies other Hamiltonian corrections, either up front via `modify_hamiltonian` or as post-hoc corrections during the solve), resolves the degenerate-state structure, computes energies/wavefunction corrections/transition moments/operator values (with strong-coupling-based degenerate-state extension if requested), and assembles the results into formatted output tables and spectra.

        :param states: the target states (and, optionally, initial states) to compute
        :type states: object
        :param calculate_intensities: whether to compute transition moments/IR intensities
        :type calculate_intensities: bool
        :param operator_expansions: extra operator expansions to additionally evaluate corrections for
        :type operator_expansions: object | None
        :param operator_terms: pre-expanded operator terms to additionally evaluate corrections for
        :type operator_terms: object | None
        :param operator_type: the operator symmetry/rank type, forwarded to the correction evaluator
        :type operator_type: str | int | None
        :param order: the perturbation-theory order to compute to
        :type order: int | None
        :param verbose: whether to log detailed evaluation progress
        :type verbose: bool
        :param degeneracy_specs: the degeneracy specification(s) to apply
        :type degeneracy_specs: object | None
        :param handle_degeneracies: whether to apply degenerate-perturbation-theory handling at all
        :type handle_degeneracies: bool
        :param zero_cutoff: the magnitude below which a term is treated as exactly zero
        :type zero_cutoff: float | None
        :param transition_moment_terms: restrict the transition-moment evaluation to specific term contributions
        :type transition_moment_terms: object | None
        :param hamiltonian_corrections: extra Hamiltonian corrections to apply, either up front (`'primary'`) or as post-hoc adjustments (`'degenerate'`), per `hamiltonian_correction_type`
        :type hamiltonian_corrections: object | None
        :param clear_caches: whether to clear the global analytic-solver caches before running
        :type clear_caches: bool
        :param hamiltonian_correction_type: `'primary'` to fold `hamiltonian_corrections` into the Hamiltonian itself before solving, or `'degenerate'`/other to apply them only within degenerate blocks; defaults to `self.hamiltonian_correction_modification_type`
        :type hamiltonian_correction_type: str | None
        :param only_degenerate_terms: whether the reexpressed Hamiltonian should include only strictly degenerate-coupling terms
        :type only_degenerate_terms: bool
        :param force_return_on_crash: whether to catch exceptions during the run and still return whatever partial results were computed, rather than propagating the error
        :type force_return_on_crash: bool
        :return: the computed VPT results (energies, wavefunction data, spectra, and formatted tables)
        :rtype: object
        """
        ...

    @classmethod
    def run_simple(cls, system, states, calculate_intensities=True, operator_expansions=None, operator_terms=None, operator_type=None, verbose=False, return_runner=False, degeneracy_specs=None, degeneracy_states=None, handle_degeneracies=True, zero_cutoff=None, clear_caches=True, hamiltonian_correction_type=None, hamiltonian_corrections=None, only_degenerate_terms=True, force_return_on_crash=True, **opts):
        """
        **LLM Docstring**

        Convenience one-shot entry point: builds an `AnalyticVPTRunner` for the given system (via `construct`) and immediately runs the full VPT calculation on it (via `run_VPT`).

        :param system: the molecule or system specification to run VPT on
        :type system: str | list | Molecule | VPTSystem
        :param states: the target states (and, optionally, initial states) to compute
        :type states: object
        :param calculate_intensities: whether to compute transition moments/IR intensities
        :type calculate_intensities: bool
        :param operator_expansions: extra operator expansions to additionally evaluate corrections for
        :type operator_expansions: object | None
        :param operator_terms: pre-expanded operator terms to additionally evaluate corrections for
        :type operator_terms: object | None
        :param operator_type: the operator symmetry/rank type
        :type operator_type: str | int | None
        :param verbose: whether to log detailed evaluation progress
        :type verbose: bool
        :param return_runner: whether to also return the constructed runner alongside the results
        :type return_runner: bool
        :param degeneracy_specs: the degeneracy specification(s) to apply
        :type degeneracy_specs: object | None
        :param handle_degeneracies: whether to apply degenerate-perturbation-theory handling at all
        :type handle_degeneracies: bool
        :param zero_cutoff: the magnitude below which a term is treated as exactly zero
        :type zero_cutoff: float | None
        :param clear_caches: whether to clear the global analytic-solver caches before running
        :type clear_caches: bool
        :param hamiltonian_correction_type: how to apply `hamiltonian_corrections`, forwarded to `run_VPT`
        :type hamiltonian_correction_type: str | None
        :param hamiltonian_corrections: extra Hamiltonian corrections to apply
        :type hamiltonian_corrections: object | None
        :param only_degenerate_terms: whether the reexpressed Hamiltonian should include only strictly degenerate-coupling terms
        :type only_degenerate_terms: bool
        :param force_return_on_crash: whether to catch exceptions during the run and still return partial results
        :type force_return_on_crash: bool
        :param opts: extra options forwarded to `construct`
        :type opts: dict
        :return: the computed VPT results, or `(runner, results)` if `return_runner` is set
        :rtype: object | tuple
        """
        ...