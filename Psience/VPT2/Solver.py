
import numpy as np, itertools, time, types
import scipy.sparse as sp
from collections import namedtuple
from functools import partial

from McUtils.Numputils import SparseArray
import McUtils.Numputils as nput
from McUtils.Scaffolding import Logger, NullLogger, NullCheckpointer
from McUtils.Parallelizers import Parallelizer, SerialNonParallelizer
from McUtils.Data import UnitsData
from McUtils.Combinatorics import LatticePathGenerator

from ..BasisReps import (
    Representation, BasisStateSpace, BasisMultiStateSpace,
    SelectionRuleStateSpace, BasisStateSpaceFilter, HarmonicOscillatorProductBasis, StateMaker
)

from .DegeneracySpecs import DegenerateMultiStateSpace, DegeneracySpec
from .Common import *
from .Corrections import *

__reload_hook__ = [ "..BasisReps", ".DegeneracySpecs", ".Corrections", ".Common" ]

__all__ = [
    "PerturbationTheorySolver"
]

class PerturbationTheorySolver:
    """
    A solver that applies perturbation theory
    given a series of corrections and population of states.
    Supports degenerate and non-degenerate PT.
    """

    evaluation_mode = None  # Old restored solvers retain their existing backend.
    default_wavefunction_workspace_bytes = 64*1024**2
    state_formatting = 'quanta'

    def __init__(self,
                 perturbations, states,
                 coupled_states=None,
                 order=2,
                 total_space=None,
                 flat_total_space=None,
                 state_space_iterations=None,
                 state_space_terms=None,
                 state_space_filters=None,
                 extended_state_space_filter_generator=None,
                 extended_state_space_postprocessor=None,
                 target_property_rules=None,
                 allow_sakurai_degs=False,
                 allow_post_PT_calc=True,
                 modify_degenerate_perturbations=False,
                 gaussian_resonance_handling=False,
                 ignore_odd_order_energies=False,
                 intermediate_normalization=False,
                 reorthogonalize_degenerate_states=None,
                 check_overlap=True,
                 zero_element_warning=True,
                 degenerate_states=None,
                 degeneracy_handlers=None,
                 handle_strong_couplings=False,
                 local_coupling_hamiltonian=None,
                 local_coupling_order=None,
                 # strong_coupling_test_modes=None,
                 # strong_couplings_state_filter=None,
                 # strongly_coupled_group_filter=None,
                 # strong_coupling_zero_order_energy_cutoff=4.5e-3,
                 low_frequency_mode_cutoff=1.15e-3,
                 # extend_strong_coupling_spaces=True,
                 zero_order_energy_corrections=None,
                 nondeg_hamiltonian_precision=3,
                 memory_constrained=False,
                 wavefunction_workspace_bytes=None,
                 sparse_low_order_corrections=False,
                 state_space_storage="auto",
                 state_space_memory_limit=5_000_000_000,
                 state_space_block_size=65536,
                 state_space_resident_bytes=None,
                 state_space_scratch=None,
                 correction_storage="auto",
                 correction_memory_limit=None,
                 correction_block_size=None,
                 correction_block_bytes=64*1024**2,
                 correction_spill_directory=None,
                 correction_resident_bytes=None,
                 correction_spill_bytes=None,
                 correction_spill_segment_bytes=None,
                 resonance_memory_constrained=True,
                 resonance_block_size=None,
                 hamiltonian_storage=None,
                 keep_hamiltonians=None,
                 logger=None,
                 parallelizer=None,
                 checkpointer=None,
                 results=None,
                 checkpoint_keys=None,
                 use_cached_representations=False,
                 use_cached_basis=False,
                 use_block_recurrence=True,
                 strong_coupling_property_filter=False,
                 strong_coupling_correction_filter=None,
                 state_formatting='quanta',
                 state_space_flush_buffer=None,
                 evaluation_mode=None
                 ):
        """

        :param perturbations:
        :type perturbations: Iterable[Representation]
        :param states:
        :type states: BasisStateSpace
        :param coupled_states:
        :type coupled_states: BasisMultiStateSpace
        :param order:
        :type order:
        :param degenerate_states:
        :type degenerate_states:
        :param degeneracy_mode:
        :type degeneracy_mode:
        :param strong_coupling_property_filter: enable property-use screening in strong-coupling discovery
        :param strong_coupling_correction_filter: callable(source, endpoints, correction_order),
            returning the structurally useful endpoint space. Runner supplies the property planner.
        :param use_block_recurrence: use the matrix recurrence for multi-state resonance blocks;
            False retains the original scalar recurrence and reorthogonalization path
        :type use_block_recurrence: bool
        :param logger:
        :type logger:
        :param state_formatting: 'quanta' for compact occupation labels (default),
            or 'vector' for the original excitation-vector diagnostics
        :param parallelizer:
        :type parallelizer:
        :param checkpointer:
        :type checkpointer:
        :param results:
        :type results:
        """

        # if memory_constrained:
        #     raise NotImplementedError('memory constraint handling currently broken')

        self.perts = perturbations
        self.state_formatting = state_formatting
        if state_space_storage not in {'auto', 'memory', 'external'}:
            raise ValueError('state_space_storage must be auto, memory or external')
        if state_space_memory_limit is not None and (not np.isfinite(state_space_memory_limit) or state_space_memory_limit < 0):
            raise ValueError('state_space_memory_limit must be nonnegative bytes or None')
        self.requested_state_space_storage = state_space_storage
        self.state_space_memory_limit = state_space_memory_limit
        self.state_space_block_size = int(state_space_block_size)
        if self.state_space_block_size < 1:
            raise ValueError('state_space_block_size must be positive')
        estimator = getattr(state_space_filters, 'estimate_storage', None)
        if estimator is not None:
            self.state_space_estimate = dict(estimator(states),target_property=getattr(state_space_filters,'target_property',None))
        else:
            rules = [getattr(p, 'selection_rules', None) or
                     states.basis.selection_rules(*['x']*(h+2))
                     for h,p in enumerate(perturbations[1:],1)]
            if not rules:
                rules = [states.basis.selection_rules('x','x','x')]
            self.state_space_estimate = BasisStateSpaceFilter.estimate_property_storage(
                states, states, rules, None, order=order, filtered=False)
        if state_space_storage == 'auto':
            state_space_storage = ('external' if state_space_memory_limit is not None and
                self.state_space_estimate['estimated_bytes'] > state_space_memory_limit else 'memory')
        self.state_space_storage = state_space_storage
        self._state_space_store = None
        self.state_space_flush_buffer = state_space_flush_buffer
        if state_space_storage == 'external':
            from ..BasisReps.ExternalStateSpaces import StateSpaceDiskStore
            self._state_space_store = StateSpaceDiskStore(state_space_scratch,state_space_block_size,
                memory_limit=state_space_memory_limit, resident_bytes=state_space_resident_bytes,
                flush_buffer=state_space_flush_buffer)
            self.state_space_flush_buffer = self._state_space_store.flush_buffer
        if evaluation_mode is not None and not (isinstance(evaluation_mode, str) and evaluation_mode == 'stored'):
            modes = {'stored', 'positions', 'runtime', 'diagonal'}
            if isinstance(evaluation_mode, str):
                if evaluation_mode not in modes-{'diagonal'}:
                    raise ValueError('Unknown evaluation_mode: '+evaluation_mode)
                evaluation_mode = {0: 'diagonal', **{
                    k: ('stored' if k == 1 else evaluation_mode)
                    for k in range(1, len(perturbations))}}
            if not isinstance(evaluation_mode, dict):
                raise TypeError('evaluation_mode must be a mode string or an order-to-mode dictionary')
            if any(not isinstance(k, (int, np.integer)) or k < 0 or v not in modes
                   or (k != 0 and v == 'diagonal') for k, v in evaluation_mode.items()):
                raise ValueError('Invalid perturbation-order evaluation policy')
            if any(v != 'stored' for v in evaluation_mode.values()):
                self.evaluation_mode = dict(evaluation_mode)
            implicit = any(v in ('positions', 'runtime') for v in evaluation_mode.values())
            if implicit and use_cached_representations:
                raise ValueError('Runtime/positions bindings must be rebuilt; disable use_cached_representations')
            if implicit and modify_degenerate_perturbations:
                raise ValueError('In-place modification of implicit perturbations is unsupported; use stored evaluation')
        self._reps = None
        self.order = order
        self.state_space_iterations = state_space_iterations
        self.state_space_terms = state_space_terms

        if isinstance(state_space_filters, (types.FunctionType, types.MethodType, types.LambdaType)):
            self.state_space_filter_generator = state_space_filters
            self.state_space_filters = BasisStateSpaceFilter.from_data(states, self.state_space_filter_generator(
                states, **({'compact':True} if getattr(state_space_filters,'estimate_storage',None) is not None and self._state_space_store is not None else {})))
        else:
            self.state_space_filter_generator = None
            self.state_space_filters = BasisStateSpaceFilter.from_data(states, state_space_filters)

        if isinstance(extended_state_space_filter_generator, (types.FunctionType, types.MethodType, types.LambdaType)):
            self.extended_state_space_filter_generator = extended_state_space_filter_generator
        elif extended_state_space_filter_generator is not None:
            raise ValueError("bad filter generator {}".format(extended_state_space_filter_generator))
        else:
            self.extended_state_space_filter_generator = None
        if isinstance(extended_state_space_postprocessor, (types.FunctionType, types.MethodType, types.LambdaType)):
            self.extended_state_space_postprocessor = extended_state_space_postprocessor
        elif extended_state_space_postprocessor is not None:
            raise ValueError("bad filter generator {}".format(extended_state_space_postprocessor))
        else:
            self.extended_state_space_postprocessor = None

        self.strong_coupling_property_filter = strong_coupling_property_filter
        self.strong_coupling_correction_filter = strong_coupling_correction_filter
        if self.strong_coupling_correction_filter is None:
            self.strong_coupling_correction_filter = getattr(
                state_space_filters, 'correction_filter', None)
        self.target_property_rules=target_property_rules

        self.logger = logger
        self.parallelizer = parallelizer
        self.checkpointer = checkpointer
        self.results = results
        self.checkpoint_keys = checkpoint_keys
        self.use_cached_representations=use_cached_representations
        self.use_cached_basis=use_cached_basis

        self.states = states
        self.full_basis = states.full_basis

        if degeneracy_handlers is None:
            degeneracy_handlers = [degeneracy_handlers]
        elif nput.is_atomic(degeneracy_handlers) or isinstance(degeneracy_handlers, dict):
            degeneracy_handlers = [degeneracy_handlers]
        elif not all(nput.is_atomic(d) or isinstance(d, (DegeneracySpec, dict)) for d in degeneracy_handlers):
            degeneracy_handlers = [degeneracy_handlers]
        elif isinstance(degeneracy_handlers, DegeneracySpec):
            degeneracy_handlers = [degeneracy_handlers]
        self.degeneracy_handlers = [
            DegeneracySpec.from_spec(d) if not (d is None or isinstance(d, DegeneracySpec)) else d
            for d in degeneracy_handlers
        ]

        self.degenerate_states = degenerate_states
        self._deg_states = None
        if len(self.degeneracy_handlers) == 1 and self.degeneracy_handlers[0] is None:
            if handle_strong_couplings is None:
                handle_strong_couplings = True
            if handle_strong_couplings is not False:
                self.degeneracy_handlers = [
                    DegeneracySpec.from_spec({
                        'wfc_threshold': 'auto' if handle_strong_couplings is True else handle_strong_couplings
                    })
                ]

        self.handle_strong_couplings = False
        for dspec in self.degeneracy_handlers:
            if hasattr(dspec, 'wfc_threshold'):
                self.handle_strong_couplings = True
                if dspec.iterations is None: dspec.iterations = order//2
        self.nondeg_hamiltonian_precision=nondeg_hamiltonian_precision

        self.local_coupling_hamiltonian = local_coupling_hamiltonian
        self.local_coupling_order = local_coupling_order

        # self.extend_strong_coupling_spaces = extend_strong_coupling_spaces
        # self.strong_coupling_test_modes = strong_coupling_test_modes
        # self.strong_couplings_state_filter = strong_couplings_state_filter
        # self.strongly_coupled_group_filter = strongly_coupled_group_filter
        # self.strong_coupling_zero_order_energy_cutoff = strong_coupling_zero_order_energy_cutoff
        self.low_frequency_mode_cutoff = low_frequency_mode_cutoff
        # self.degeneracy_mode = degeneracy_mode
        self.allow_sakurai_degs = allow_sakurai_degs
        self.allow_post_PT_calc = allow_post_PT_calc
        self.zero_order_energy_corrections = zero_order_energy_corrections
        self.ignore_odd_orders = ignore_odd_order_energies
        self.drop_perturbation_degs = modify_degenerate_perturbations
        self.intermediate_normalization = intermediate_normalization
        if reorthogonalize_degenerate_states is None:
            reorthogonalize_degenerate_states = not intermediate_normalization
        self.reorthogonalize_degenerate_states = reorthogonalize_degenerate_states
        self.use_block_recurrence = use_block_recurrence
        self.check_overlap = check_overlap
        self.gaussian_resonance_handling = gaussian_resonance_handling
        self.zero_element_warning = zero_element_warning

        self.state_space_storage = state_space_storage
        self.state_space_scratch = state_space_scratch
        self.correction_storage = correction_storage
        self.correction_block_bytes = correction_block_bytes
        self.correction_spill_directory = correction_spill_directory
        self.correction_resident_bytes = correction_resident_bytes
        self.correction_spill_bytes = correction_spill_bytes
        self.correction_spill_segment_bytes = correction_spill_segment_bytes
        if correction_spill_bytes is not None and correction_storage != 'compact':
            raise ValueError('Correction spill requires explicit correction_storage="compact"')
        self.correction_memory_limit = correction_memory_limit
        self.correction_block_size = self.state_space_block_size if correction_block_size is None else correction_block_size
        self.resonance_memory_constrained = resonance_memory_constrained
        self.resonance_block_size = self.state_space_block_size if resonance_block_size is None else resonance_block_size
        if not isinstance(self.resonance_block_size,(int,np.integer)) or self.resonance_block_size < 1:
            raise ValueError('resonance_block_size must be a positive integer')
        self._correction_default_memory_limit = state_space_memory_limit
        PerturbationTheoryCorrections.correction_storage_policy(correction_storage, correction_memory_limit,
                                  self.correction_block_size, state_space_memory_limit, (0, 0, 0))
        self.hamiltonian_storage = hamiltonian_storage
        self.memory_constrained=memory_constrained
        if wavefunction_workspace_bytes is None:
            wavefunction_workspace_bytes = self.default_wavefunction_workspace_bytes
        if (not isinstance(wavefunction_workspace_bytes, (int, np.integer))
                or wavefunction_workspace_bytes < 1):
            raise ValueError('wavefunction_workspace_bytes must be positive integer bytes')
        self.wavefunction_workspace_bytes = int(wavefunction_workspace_bytes)
        if sparse_low_order_corrections and (not memory_constrained or
                correction_storage not in {'sparse', 'compact'}):
            raise ValueError('sparse_low_order_corrections requires memory_constrained and sparse/compact storage')
        self.sparse_low_order_corrections = bool(sparse_low_order_corrections)
        self.keep_hamiltonians=keep_hamiltonians

        self._coupled_states = coupled_states
        self._total_space = total_space
        self._flat_space = flat_total_space
        self._total_dim = None

        self._zo_engs = None
    @property
    def coupled_states(self):
        """

        :return:
        :rtype:
        """
        if self._coupled_states is None:
            self.load_state_spaces()
        return self._coupled_states
    @property
    def total_space_dim(self):
        """

        :return:
        :rtype:
        """
        if self._total_dim is None:
            self.load_state_spaces()
        return self._total_dim
    @property
    def flat_total_space(self):
        """

        :return:
        :rtype:
        """
        if self._flat_space is None:
            self.load_state_spaces()
        return self._flat_space
    @property
    def total_state_space(self):
        """

        :return:
        :rtype:
        """
        if self._total_space is None:
            self.load_state_spaces()
        return self._total_space

    class PastIndexableTuple(list):
        def __getitem__(self, item):
            """
            **LLM Docstring**

            Index into the list, transparently returning `0` for any out-of-range integer index instead of raising, so higher-order perturbation terms that weren't explicitly computed are treated as identically zero.

            :param item: the index (or slice) to look up
            :type item: int | slice
            :return: the stored value, or `0` if `item` is an out-of-range integer index
            :rtype: object
            """
            if isinstance(item, (int, np.integer)) and item >= len(self):
                return 0
            else:
                return super().__getitem__(item)
        def __setitem__(self, key, value):
            """
            **LLM Docstring**

            Set the value at a given index, first extending the list with zero placeholders if the index is beyond the current length. Note the extension amount (`len(self) - key`) is non-positive whenever `key >= len(self)`, so `extend` receives an empty (or negative-length, effectively empty) list and the list is never actually grown to include `key` -- this looks like an off-by-sign bug (it should extend by `key - len(self) + 1`), and the subsequent `super().__setitem__` call would raise an `IndexError` for any `key` at or beyond the original length.

            :param key: the index to set
            :type key: int
            :param value: the value to store
            :type value: object
            :return: the result of the (likely failing, for out-of-range `key`) base `list.__setitem__` call
            :rtype: object
            """
            if isinstance(key, (int, np.integer)) and key >= len(self):
                self.extend([0]*(len(self)-key))
            return super().__setitem__(key, value)

    @property
    def representations(self):
        """
        :return:
        :rtype: Iterable[SparseArray]
        """
        if self._reps is None:
            self._reps = self.PastIndexableTuple(self.get_VPT_representations())
        return self._reps
    @representations.setter
    def representations(self, reps):
        """
        **LLM Docstring**

        The (cached) list of Hamiltonian perturbation representation matrices, computed lazily via `get_VPT_representations` and wrapped in a `PastIndexableTuple` so higher, uncomputed orders read back as zero.

        :return: the perturbation representation matrices
        :rtype: PerturbationTheorySolver.PastIndexableTuple
        """
        self._reps = reps

    @property
    def degenerate_spaces(self):
        """
        :return:
        :rtype:
        """
        if self._deg_states is None:
            new_states = []
            if self.degenerate_states is not None:
                new_states.append(
                    DegenerateMultiStateSpace.from_spec(
                        self.degenerate_states,
                        solver=self,
                        full_basis=self.full_basis,
                        group_filter=None,
                        log_groups=True
                    )
                )
            for spec in self.degeneracy_handlers:
                if spec is not None and not isinstance(spec, DegeneracySpec):
                    raise ValueError("expected `{}` object".format(DegeneracySpec.__name__))
                if spec is not None and spec.application_order != 'pre':
                    spec = None

                if spec is None:
                    group_filter = None
                else:
                    if spec.group_filter is not None:
                        group_filter = spec.get_degenerate_group_filter(self)
                    else:
                        group_filter = None
                new_states.append(DegenerateMultiStateSpace.from_spec(
                    spec,
                    solver=self,
                    full_basis=self.full_basis,
                    group_filter=group_filter,
                    default_group_filter=(partial(self.filter_degenerate_group, block_size=self.resonance_block_size, state_formatting=self.state_formatting)
                        if self.resonance_memory_constrained else None),
                    log_groups=True
                ))
            if len(new_states) > 0:
                if len(new_states) == 1:
                    self._deg_states = new_states[0]
                else:
                    self._deg_states = self.merge_deg_spaces(new_states)
            else:
                self._deg_states = None
        return self._deg_states

    @classmethod
    def merge_deg_spaces(cls, new_states):
        """
        **LLM Docstring**

        Combine several independently-identified sets of degenerate state groups into one consistent set, by flattening every group (from either raw `BasisStateSpace`s or `DegenerateMultiStateSpace`s) down to its excitation vectors and merging any that share states via `DegeneracySpec.merge_state_blocks`.

        :param new_states: the separate sets of degenerate groups to merge, each either a `DegenerateMultiStateSpace` or an iterable of raw state blocks
        :type new_states: list
        :return: the merged degenerate groups
        :rtype: list[np.ndarray]
        """
        all_blocks = [
            b.excitations if isinstance(b, BasisStateSpace) else b
            for deg_spaces in new_states
            for b in (
                deg_spaces.flat if isinstance(deg_spaces, DegenerateMultiStateSpace) else deg_spaces
            )
        ]
        merged = DegeneracySpec.merge_state_blocks(all_blocks)
        basis = None
        for b in new_states:
            if isinstance(b, DegenerateMultiStateSpace):
                basis = b.basis
                break
        else:
            basis = HarmonicOscillatorProductBasis(len(all_blocks[0][0]))
        return DegenerateMultiStateSpace(
            [BasisStateSpace(basis, b) for b in merged]
        )


    @property
    def zero_order_energies(self):
        """
        :return:
        :rtype:
        """
        if self._zo_engs is None:
            H0 = self.representations[0]
            e_vec_full = (H0._mapped_diagonal if hasattr(H0,'_mapped_diagonal') else
                          np.diag(H0) if isinstance(H0, np.ndarray) else H0.diag)
            if isinstance(e_vec_full, SparseArray):
                e_vec_full = e_vec_full.asarray()
            self._zo_engs = e_vec_full
            if self.zero_order_energy_corrections is not None:
                if callable(self.zero_order_energy_corrections):
                    if self._state_space_store is None:
                        energies = self.zero_order_energy_corrections(self.flat_total_space.excitations)
                    else:
                        energies = self._state_space_store.array('zero_order',float,(len(self.flat_total_space),))
                        offset = 0
                        for block in self.flat_total_space.iter_excitations():
                            energies[offset:offset+len(block)] = self.zero_order_energy_corrections(block)
                            offset += len(block)
                    self._zo_engs = energies
                else:
                    states=np.asarray([s[0] for s in self.zero_order_energy_corrections])
                    energies=np.asarray([s[1] for s in self.zero_order_energy_corrections])
                    inds = self.flat_total_space.find(states) # slow but w/e
                    biggies = np.where(inds >= len(self._zo_engs))
                    if len(biggies) > 0:
                        if len(biggies[0]) > 0:
                            for k,v in zip(states[biggies], energies[biggies]):
                                self.logger.log_print(
                                    "WARNING: zero-order correction {} for state {} not used as state is not in basis".format(
                                        v * UnitsData.convert("Hartrees", "Wavenumbers"),
                                        StateMaker.format_state(k, self.state_formatting)
                                    )
                                )
                    # for k,v in self.zero_order_energy_corrections:
                    #     if not isinstance(k, int):
                    #         k = self.flat_total_space.find([k]) # slow but w/e
                    #     if k >= len(self._zo_engs):
                    #         self.logger.log_print(
                    #             "WARNING: zero-order correction {} for state {} not used as state is not in basis".format(
                    #                 v * UnitsData.convert("Hartrees", "Wavenumbers"),
                    #                 k
                    #             )
                    #         )
                    self._zo_engs[inds] = energies
                if self.evaluation_mode is not None and self.evaluation_mode.get(0, 'stored') != 'stored':
                    self.representations[0] = nput.EvaluatedSparseMatrix.from_diagonal_data(self._zo_engs)
                else:
                    self.representations[0] = SparseArray.from_diag(self._zo_engs)
        return self._zo_engs

    def apply_VPT(self):
        """
        Applies perturbation theory to the held basis of states using the
        built representations and degenerate state spaces

        :return:
        :rtype: PerturbationTheoryCorrections
        """

        corrs = self.get_corrections()
        degenerate_states = corrs.degeneracies
        corrs = corrs.corrections
        corrs.state_formatting = self.state_formatting

        # degeneracy_mode = self.degeneracy_mode
        # degenerate_states = self.degenerate_spaces
        if (
                self.allow_post_PT_calc
                and degenerate_states is not None
                and any(len(x) > 1 for x in degenerate_states)
        ):
            with self.logger.block(tag="Applying post-PT variational calc."):
                if self.gaussian_resonance_handling:
                    self.logger.log_print('WARNING: Doing Gaussian resonance handling and not doing variational calculation involving states with more than 2 quanta of excitation')
                deg_engs, deg_transf, deg_hams = self.apply_post_PT_variational_calc(degenerate_states, corrs)
                corrs.degenerate_energies = deg_engs
                corrs.degenerate_transf = deg_transf
                corrs.degenerate_hamiltonians = deg_hams

        return corrs

    #region Get Matrix Inputs
    def get_VPT_representations(self):
        """
        Gets the sparse representations of the passed perturbation inside the basis of coupled states.

        :return:
        :rtype: Iterable[SparseArray]
        """

        logger = self.logger
        if logger is None:
            logger = NullLogger()

        with self.checkpointer as checkpointer:
            with self.logger.block(tag='getting representations'):
                self.logger.log_print('trying to load from checkpoint...')
                try:
                    if self.use_cached_representations:
                        H = checkpointer['representations']
                    else:
                        H = None
                except (OSError, KeyError):
                    H = None
                if (H is not None and self.evaluation_mode is not None
                        and self.evaluation_mode.get(0, 'stored') == 'diagonal'):
                    H = list(H)
                    if not isinstance(H[0], nput.EvaluatedSparseMatrix):
                        old = H[0]
                        values = old._mapped_diagonal if hasattr(old, '_mapped_diagonal') else old.diag
                        H[0] = nput.EvaluatedSparseMatrix.from_diagonal_data(values)
                if H is None:
                    self.logger.log_print('failed to load, building instead...')

                    par = Parallelizer.lookup(self.parallelizer)
                    with par:  # we put an outermost block here to just make sure everything is clean

                        # diag_inds = BraKetSpace(self.flat_total_space, self.flat_total_space)
                        # N = len(self.flat_total_space)

                        n_spaces = len(self.total_state_space.spaces)
                        # raise Exception(len(self.total_state_space.spaces))
                        H = [np.zeros(1)] * min(len(self.perts), n_spaces)
                        with logger.block(tag="building {}".format(self.perts[0])):
                            start = time.time()
                            H[0] = self.perts[0].get_representation_matrix(self.flat_total_space, self.flat_total_space,
                                                                           zero_element_warning=self.zero_element_warning,
                                                                           diagonal=True,
                                                                           **({} if self.evaluation_mode is None else dict(
                                                                               evaluation_mode=self.evaluation_mode.get(0, 'stored'),
                                                                               evaluation_working_memory=self.wavefunction_workspace_bytes))
                                                                           )
                            end = time.time()
                            logger.log_print("took {t:.3f}s", t=end - start)
                            self.perts[0].clear_cache()
                            # if self.memory_constrained and not isinstance(self.checkpointer, NullCheckpointer):
                            #     self.checkpointer['H[0]'] = H[0]
                            #     H[0] = None
                            #     gc.collect()

                        # import tracemalloc
                        for i, h in enumerate(self.perts[1:]):
                            # calculate matrix elements in the coupled subspace
                            if n_spaces > i + 1:
                                cs = self.total_state_space[i + 1]
                                with logger.block(tag="building {}".format(h)):
                                    start = time.time()
                                    # tracemalloc.start()
                                    # pre = tracemalloc.take_snapshot()
                                    H[i + 1] = h.get_representation_matrix(cs, self.flat_total_space,
                                                                           zero_element_warning=self.zero_element_warning,
                                                                           diagonal=False,
                                                                           memory_constrained=self.memory_constrained,
                                                                           deduplication_method=self.state_space_storage,
                                                                           deduplication_directory=self.state_space_scratch,
                                                                           deduplication_flush_buffer=self.state_space_flush_buffer,
                                                                           **({} if self.evaluation_mode is None else dict(
                                                                               evaluation_mode=self.evaluation_mode.get(i+1, 'stored'),
                                                                               evaluation_working_memory=self.wavefunction_workspace_bytes))
                                                                           )
                                    h.clear_cache()
                                    end = time.time()
                                    logger.log_print("took {t:.3f}s", t=end - start)
                                    # raise Exception("woop")

                    # raise Exception("...")

                    if self.use_cached_representations:
                        try:
                            checkpointer['representations'] = H
                        except KeyError:
                            pass

            if self.hamiltonian_storage is not None and self._state_space_store is None:
                from pathlib import Path
                from McUtils.Numputils import save_sparse_matrix, load_sparse_matrix
                directory = Path(self.hamiltonian_storage)
                directory.mkdir(parents=True, exist_ok=True)
                for order, matrix in enumerate(H):
                    if isinstance(matrix, SparseArray) and not isinstance(matrix, nput.EvaluatedSparseMatrix):
                        path = directory / "H{}.npz".format(order)
                        save_sparse_matrix(path, matrix)
                        H[order] = load_sparse_matrix(path)

            return H
    def extend_VPT_representations(self, new_flat_space, new_states):
        """
        **LLM Docstring**

        Extend the cached Hamiltonian perturbation representation matrices to cover a newly added block of coupled states, by computing just the new matrix elements (between the new states and the full flat space) and splicing them into the existing sparse representation data rather than recomputing everything from scratch.

        :param new_flat_space: the newly added states, flattened, that the zeroth-order representation needs new diagonal elements for
        :type new_flat_space: BasisStateSpace
        :param new_states: the newly added coupled-state blocks, one per higher-order perturbation, that need new off-diagonal elements computed
        :type new_states: list[BasisStateSpace]
        :return: the extended list of representation matrices
        :rtype: list[SparseArray]
        """

        if self._state_space_store is not None or (self.evaluation_mode is not None and
                any(mode != 'stored' for mode in self.evaluation_mode.values())):
            # Rebuild from the merged recipes into mapped CSR. The resident
            # splice below materializes all old/new COO entries for np.unique.
            self._zo_engs = None
            cached = self.use_cached_representations
            self.use_cached_representations = False
            try:
                H = self.get_VPT_representations()
            finally:
                self.use_cached_representations = cached
            if cached:
                with self.checkpointer as checkpointer:
                    try:
                        checkpointer['representations'] = H
                    except KeyError:
                        pass
            return H

        logger = self.logger
        if logger is None:
            logger = NullLogger()

        with self.checkpointer as checkpointer:
            with self.logger.block(tag='extending representations'):
                par = Parallelizer.lookup(self.parallelizer)
                with par:
                    n_spaces = len(new_states) + 1
                    N = len(self.flat_total_space)
                    H = [np.zeros(1)] * min(len(self.perts), n_spaces)
                    with logger.block(tag="building {}".format(self.perts[0])):
                        start = time.time()
                        subh = self.perts[0].get_representation_matrix(new_flat_space,
                                                                       self.flat_total_space,
                                                                       zero_element_warning=self.zero_element_warning,
                                                                       diagonal=True
                                                                       )
                        # old_N = self.representations[0].shape[0]
                        old_vals, old_inds = self.representations[0].block_data
                        new_vals, new_inds = subh.block_data
                        full_dat = np.concatenate([old_vals, new_vals])
                        full_inds = tuple(
                            np.concatenate([o.astype(int), i.astype(int)]) for o,i in zip(old_inds, new_inds)
                        )
                        H[0] = SparseArray.from_data(
                            (
                                full_dat,
                                full_inds
                            ),
                            shape=(N, N)
                        )
                        end = time.time()
                        logger.log_print("took {t:.3f}s", t=end - start)
                        self.perts[0].clear_cache()
                        self._zo_engs = None

                    # dm = np.abs(
                    #     self.zero_order_energies[np.newaxis, :] - self.zero_order_energies[:, np.newaxis]
                    # )
                    # np.fill_diagonal(dm, 10)
                    #
                    # print(np.min(dm))
                    #
                    # self._get_Pi0([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2]])
                    # raise Exception(...)

                    for i, h in enumerate(self.perts[1:]):
                        if n_spaces > i + 1:
                            cs = new_states[i]
                            with logger.block(tag="building {}".format(h)):
                                start = time.time()
                                if len(cs) > 0:
                                    #TODO: add a way to avoid calculating terms
                                    #      that I've already calculated
                                    subh = h.get_representation_matrix(cs,
                                                                       self.flat_total_space,
                                                                       zero_element_warning=self.zero_element_warning,
                                                                       diagonal=False,
                                                                       memory_constrained=self.memory_constrained
                                                                       )
                                    old_vals, old_inds = self.representations[i+1].block_data
                                    new_vals, new_inds = subh.block_data

                                    #TODO: I wouldn't need this if I could be assured of no dupes
                                    full_dat = np.concatenate([old_vals, new_vals])
                                    full_inds = tuple(np.concatenate([o, i]) for o, i in zip(old_inds, new_inds))
                                    uinds, usort = np.unique(
                                        np.array(full_inds).T,
                                        return_index=True,
                                        axis=0
                                    )
                                    full_dat = full_dat[usort]
                                    full_inds = uinds.T
                                else:
                                    full_dat, full_inds = self.representations[i+1].block_data
                                H[i + 1] = SparseArray.from_data(
                                    (
                                        full_dat,
                                        full_inds
                                    ),
                                    shape=(N, N)

                                )
                                h.clear_cache()
                                end = time.time()
                                logger.log_print("took {t:.3f}s", t=end - start)

                # raise Exception("...")

                if self.use_cached_representations:
                    try:
                        checkpointer['representations'] = H
                    except KeyError:
                        pass

            return H

    def _take_subham(self, rep, inds):
        """
        Builds a subsampled version of a representation Hamiltonian
        to allow equations to be efficiently solved in subspaces.

        :param rep: representation matrix from which to pull the subspace
        :type rep: SparseArray
        :param inds: indices for the subspace
        :type inds: np.ndarray
        :return:
        :rtype:
        """
        ind_pairs = np.array(list(itertools.product(inds, inds))).T
        return np.reshape(rep[ind_pairs], (len(inds), len(inds)))
    def _build_projector(self, inds):
        """
        Builds a projector where only inds will
        be included

        :param inds: indices for the subspace
        :type inds: np.ndarray
        :return:
        :rtype:
        """
        shp = self.representations[0].shape
        return SparseArray.from_data(
                (
                    np.ones(len(inds)),
                    (inds, inds)
                ),
                shape=shp
            )
    def _get_Pi0(self, degenerate_subspace, non_zero_cutoff=None, E0=None):
        """
        **LLM Docstring**

        Build the perturbation operator `Pi0`
        -- a diagonal matrix of `1/(E0_ref - E0_i)` for every state `i` outside the given degenerate subspace, with zero entries for states inside it --
        raising an error if any state outside the subspace is unexpectedly (near-)degenerate with the reference energy.

        :param degenerate_subspace: the indices of the states in the (single) degenerate group `Pi0` should exclude/zero out
        :type degenerate_subspace: np.ndarray
        :param non_zero_cutoff: the energy-gap magnitude below which a state is considered problematically degenerate; defaults to `Settings.non_zero_cutoff`
        :type non_zero_cutoff: float | None
        :param E0: the reference zeroth-order energy to measure gaps relative to; defaults to the average energy over `degenerate_subspace`
        :type E0: float | np.ndarray | None
        :return: the diagonal `Pi0` operator
        :rtype: SparseArray
        :raises ValueError: if more than one state outside `degenerate_subspace` is unexpectedly near-degenerate with `E0`
        """
        # generate the perturbation operator
        e_vec_full = self.zero_order_energies
        if E0 is None:
            E0 = np.average(e_vec_full[degenerate_subspace]) # better to use the first or the average? Not clear...not even sure this is a useful codepath
        e_vec = e_vec_full - E0
        e_vec[degenerate_subspace] = 1
        if non_zero_cutoff is None:
            non_zero_cutoff = Settings.non_zero_cutoff
        zero_checks = np.where(np.abs(e_vec) < non_zero_cutoff)[0]
        if len(zero_checks) > 0:
            if isinstance(E0, (int, float, np.integer, np.floating)):
                Et = [E0]
            else:
                Et = E0
            bad_vec = np.concatenate([Et, e_vec_full[zero_checks]])
            if len(zero_checks) > 10:
                #TODO: add better message to specify which members of degenerate subspace
                raise ValueError(
                    "degeneracies encountered: states {} and {} other states are degenerate (average energy: {} stddev: {})".format(
                        degenerate_subspace,
                        len(zero_checks),
                        np.average(bad_vec),
                        np.std(bad_vec)
                    ))
            else:
                #TODO: add better message
                raise ValueError(
                    "degeneracies encountered: states {} and {} are degenerate (average energy: {} stddev: {})".format(
                        StateMaker.format_states(self.flat_total_space.take_subspace(degenerate_subspace), self.state_formatting, join=True),
                        StateMaker.format_states(self.flat_total_space.take_subspace(zero_checks), self.state_formatting, join=True),
                        np.average(bad_vec),
                        np.std(bad_vec)
                    ))
        pi = 1 / e_vec
        pi[degenerate_subspace] = 0
        return SparseArray.from_diag(pi)
    #endregion

    #region Get Coupled Spaces
    def load_state_spaces(self):
        """

        :return:
        :rtype:
        """

        logger = self.logger
        logger.log_print('state-space storage: {storage}; estimated state-space RAM {cost} bytes',
                         storage=self.state_space_storage, cost=self.state_space_estimate['estimated_bytes'])
        with self.logger.block(tag='getting basis'):
            if self._coupled_states is None:
                self.logger.log_print('trying to load from checkpoint...')
                with self.checkpointer:
                    try:
                        if self.use_cached_basis:
                            self._coupled_states = self.checkpointer['coupled_states']
                    except (OSError, KeyError):
                        self._coupled_states = None
                    if self._coupled_states is None:
                        self.logger.log_print('fail to load, building instead...')
                        parallelizer = Parallelizer.lookup(self.parallelizer)
                        if parallelizer.nprocs > 1:
                            parallelizer.printer = self.logger.log_print
                            self.logger.log_print('parallelizing over {nproc} processors',
                                                  nproc=parallelizer.nprocs
                                                  )
                        with parallelizer:  # we put an outermost block here to just make sure everything is clean

                            start = time.time()
                            self._coupled_states = self.load_coupled_spaces()
                            # _ = [len(s) for s in [self.states] + list(self._coupled_states)]
                            end = time.time()

                            self.logger.log_print(
                                ['H({i}): {s}'.format(i=i, s=s) for i,s in enumerate([self.states] + list(self._coupled_states))]
                                + ["took {t}s..."],
                                t=round(end - start, 3)
                            )
                        # raise Exception('break')
                        if self.use_cached_basis:
                            try:
                                self.checkpointer['coupled_states'] = self._coupled_states
                            except KeyError:
                                pass

                    # raise Exception('break')

            elif len(self._coupled_states) != len(self.perts) - 1:
                raise ValueError("coupled states must be specified for all perturbations (got {}, expected {})".format(
                    len(self._coupled_states),
                    len(self.perts) - 1
                ))
            elif any(not isinstance(cs, (BasisStateSpace, SelectionRuleStateSpace)) for cs in self._coupled_states):
                self._coupled_states = [
                    BasisStateSpace(self.states.basis, cs, full_basis=self.full_basis)
                        if not isinstance(cs, (BasisStateSpace, SelectionRuleStateSpace)) else
                    cs
                    for cs in self._coupled_states
                ]

            with logger.block(tag="precomputing coupled space indices"):
                # otherwise they get calculated twice
                start = time.time()
                for s in self._coupled_states:
                    if s is not None and self._state_space_store is None:
                        logger.log_print('generating indices for {s}', s=s)
                        new = s.indices

                # inds = [s.indices for s in self._coupled_states]
                end = time.time()
                logger.log_print(
                    [
                        "took {t:.3f}s"
                    ],
                    t=end - start
                )

            self._generate_total_space()
    def _generate_total_space(self):
        """
        **LLM Docstring**

        Build (if not already present) the total coupled state space and its flattened, deduplicated form spanning the target states plus every coupled-state block, computing and logging the total dimension either from scratch or by flattening an already-known total space.

        :return: None
        :rtype: None
        """
        if self._total_space is None:
            with self.logger.block(tag="generating total space"):

                start = time.time()

                space_list = [self.states] + [s for s in self._coupled_states if s is not None]
                self._total_space = BasisMultiStateSpace(np.array(space_list, dtype=object))
                if self._state_space_store is not None:
                    def rank_blocks():
                        yield self.states.indices
                        for space in self._coupled_states:
                            if space is not None:
                                yield from space.iter_index_blocks()
                    self._flat_space = self._state_space_store.frontier(self.states.basis,rank_blocks())
                    flat_space = self._flat_space
                else:
                    flat_space = self.states.take_unique().to_single(track_excitations=False)
                    for s in self._coupled_states:
                        if s is not None:
                            flat_space = flat_space.union(s.to_single(track_excitations=False), track_excitations=False)
                    # flat_space = self._total_space.to_single()
                    self._flat_space = flat_space.take_unique(track_excitations=False)
                # raise Exception(
                #     self._flat_space.find(
                #         space_list[1].get_representation_brakets(other=None).bras
                #     )
                # )
                self._total_dim = len(self.flat_total_space)

                end = time.time()
                self.logger.log_print(
                    [
                        "total coupled space dimension: {d} (contracted from {f})",
                        "took {t:.3f}s"
                    ],
                    d=self.total_space_dim,
                    f=len(flat_space),
                    t=end - start
                )
                # raise Exception("break")
        else:
            if self._flat_space is None:
                with self.logger.block(tag="generating total space"):
                    start = time.time()
                    self._flat_space = self._total_space.to_single().take_unique()
                    self._total_dim = len(self._flat_space)

                    end = time.time()
                    self.logger.log_print(
                        [
                            "total coupled space dimension: {d}",
                            "took {t:.3f}s"
                        ],
                        d=self.total_space_dim,
                        t=end - start
                    )
            else:
                self._total_dim = len(self._flat_space)

    def extend_state_spaces(self, new_targets, degenerate_states=None):
        """
        **LLM Docstring**

        Extend the solver's state spaces to additionally include a set of new target states: (re)computes the coupled-state blocks needed for the new targets (applying any configured extended-state-space filter/postprocessor), merges them into the existing per-order coupled-state spaces and flattened total space, and updates `self.states`/`self._total_space` accordingly.

        :param new_targets: the new target states to add to the solve
        :type new_targets: BasisStateSpace
        :param degenerate_states: the degenerate-state groupings driving this extension, forwarded to `extended_state_space_postprocessor` if one is configured
        :type degenerate_states: object | None
        :return: `(flat_space, new_spaces)` -- the newly added (deduplicated) flat states and the per-order newly added coupled-state blocks, or `None` if nothing new was actually found
        :rtype: tuple | None
        """
        with self.logger.block(tag='extending basis'):
            parallelizer = Parallelizer.lookup(self.parallelizer)
            if parallelizer.nprocs > 1:
                parallelizer.printer = self.logger.log_print
                self.logger.log_print('parallelizing over {nproc} processors',
                                      nproc=parallelizer.nprocs
                                      )
            with parallelizer:  # we put an outermost block here to just make sure everything is clean
                start = time.time()
                existing_spaces = {self.perts[0]:None}
                for p,cs in zip(self.perts[1:], self.coupled_states):
                    existing_spaces[p] = ({None:cs}, cs)
                if self.extended_state_space_filter_generator is not None:
                    filters = BasisStateSpaceFilter.from_data(new_targets,
                                                              self.extended_state_space_filter_generator(new_targets, check_subspaces=False, **({'compact':True} if self._state_space_store is not None and getattr(self.extended_state_space_filter_generator,'estimate_storage',None) is not None else {}))
                                                              )
                elif self.state_space_filter_generator is not None:
                    filters = BasisStateSpaceFilter.from_data(new_targets,
                                                              self.state_space_filter_generator(new_targets, check_subspaces=False, **({'compact':True} if self._state_space_store is not None and getattr(self.state_space_filter_generator,'estimate_storage',None) is not None else {}))
                                                              )
                else:
                    filters = None #TODO: is this really the right choice...?
                new_spaces = self.load_coupled_spaces([new_targets],
                                                      filter_spaces=filters
                                                      )#, spaces=existing_spaces)
                if self.extended_state_space_postprocessor is not None:
                    new_spaces = self.extended_state_space_postprocessor(new_targets, new_spaces, degenerate_states)
                self.states = self.states.union(new_targets)

                if all(n is None for n in new_spaces):
                    return None
                for i,(old,new) in enumerate(zip(self._coupled_states, new_spaces)):
                    if new is not None:
                        _ = new.indices # precomputing
                        self._coupled_states[i] = old.union(new)
                        # print(new)
                        # new_spaces[i] = new.difference(old)
                        # print(new)
                # _ = [len(s) for s in [self.states] + list(self._coupled_states)]
                end = time.time()

                self.logger.log_print(
                    ['H({i}): {s}'.format(i=i, s=s) for i, s in
                     enumerate([self.states] + list(self._coupled_states))]
                    + ["took {t}s..."],
                    t=round(end - start, 3)
                )

            flat_space = None
            for s in new_spaces:
                if s is not None:
                    if flat_space is None:
                        flat_space = s.to_single(track_excitations=False)
                    else:
                        flat_space = flat_space.union(s.to_single(track_excitations=False), track_excitations=False)
            flat_space = flat_space.take_unique(track_excitations=False)
            flat_space = flat_space.difference(self.flat_total_space)

            new_spaces = [s for s in new_spaces if s is not None]

            self._flat_space = self.flat_total_space.union(flat_space, track_excitations=False)
            self._total_dim = len(self.flat_total_space)

            space_list = [self.states] + [s for s in self._coupled_states if s is not None]
            self._total_space = BasisMultiStateSpace(np.array(space_list, dtype=object))

        return flat_space, new_spaces

    def load_coupled_spaces(self,
                            degenerate_spaces=None,
                            spaces=None,
                            wavefunction_terms=None,
                            property_filter=None,
                            filter_spaces=None
                            ):
        """
        Determines which states need to be coupled at which levels of correction
        to handle the PT
        :return:
        :rtype:
        """

        total_state_spaces = []
        # loop over the degenerate sets and build a full
        # set of connected states
        simple_spaces = []
        if degenerate_spaces is None:
            degenerate_spaces = self.degenerate_spaces
        for deg_group in degenerate_spaces:
            # self.logger.log_print('loading {g} coupled states...', g=deg_group.indices)
            if len(deg_group) > 1 and self.allow_sakurai_degs:
                raise NotImplementedError("True degeneracy handling needs some patches")
                second_deg = self._use_second_deg_PT(deg_group)
                deg_space = deg_group
                spaces = self.get_coupled_space(None, deg_space, second_deg,
                                                      allow_PT_degs=True, spaces=spaces)
                total_state_spaces.append(spaces)
            else:
                deg_space = None if self.allow_sakurai_degs else deg_group
                simple_spaces.append(deg_space)

        if len(simple_spaces) > 0:
            space = simple_spaces[0]
            for s in simple_spaces[1:]:
                space = space.union(s)

            # raise Exception(simple_spaces[0].full_basis)

            if wavefunction_terms is None:
                wavefunction_terms = self.state_space_terms
            if property_filter is None:
                property_filter = self.target_property_rules
            if filter_spaces is None:
                filter_spaces = self.state_space_filters

            spaces = self.get_coupled_space(
                space,
                None, False,
                allow_PT_degs=self.allow_sakurai_degs,
                spaces=spaces,
                wavefunction_terms=wavefunction_terms,
                property_filter=property_filter,
                filter_spaces=filter_spaces
                )
            total_state_spaces.append(spaces)

        coupled_states = [spaces[h][1] if spaces[h] is not None else None for h in self.perts]

        return coupled_states[1:]

    class StateSpaceWrapper:
        """
        Wraps a state space so that it can define stuff like __add__, __mul__, and __neg__
        """
        def __init__(self, space):
            """
            **LLM Docstring**

            Wrap a state space so it supports the symbolic `+`/`-`/unary-`-` operations used when combining state-space contributions symbolically.

            :param space: the state space to wrap
            :type space: BasisStateSpace | object
            :return: None
            :rtype: None
            """
            self.space = space
        def __neg__(self):
            """
            **LLM Docstring**

            Negation is a no-op for state spaces (there's no notion of a "negative" set of states), so this just returns `self` unchanged.

            :return: `self`
            :rtype: PerturbationTheorySolver.StateSpaceWrapper
            """
            return self
        def simple_union(self, other):
            """
            **LLM Docstring**

            Combine this wrapped state space with another (treating `None`/`0` as "nothing to add"), unwrapping `other` first if it's also a `StateSpaceWrapper`.

            :param other: the other state space (or wrapper) to union with; `None` or `0` are treated as empty
            :type other: object | PerturbationTheorySolver.StateSpaceWrapper | None
            :return: `self` if `other` is empty, otherwise a new wrapper around the unioned space
            :rtype: PerturbationTheorySolver.StateSpaceWrapper
            """
            if (
                other is None
                or isinstance(other, (int, np.integer, float, np.floating)) and other == 0
            ):
                return self
            if isinstance(other, type(self)):
                other = other.space
            if other is None:
                return self
            else:
                return type(self)(self.space.union(other))
        def __sub__(self, other):
            """
            **LLM Docstring**

            Since state spaces have no subtraction, `-` is treated the same as `+` -- both simply union the two spaces together.

            :param other: the other state space (or wrapper) to combine with
            :type other: object
            :return: the unioned wrapper
            :rtype: PerturbationTheorySolver.StateSpaceWrapper
            """
            return self.simple_union(other)
        def __rsub__(self, other):
            """
            **LLM Docstring**

            Reflected subtraction; identical to `__sub__`/`simple_union` since combining spaces is commutative and subtraction is treated as union.

            :param other: the other state space (or wrapper) to combine with
            :type other: object
            :return: the unioned wrapper
            :rtype: PerturbationTheorySolver.StateSpaceWrapper
            """
            return self.simple_union(other)
        def __add__(self, other):
            """
            **LLM Docstring**

            Combine this wrapped state space with another via union.

            :param other: the other state space (or wrapper) to combine with
            :type other: object
            :return: the unioned wrapper
            :rtype: PerturbationTheorySolver.StateSpaceWrapper
            """
            return self.simple_union(other)
        def __radd__(self, other):
            """
            **LLM Docstring**

            Reflected addition; identical to `__add__` since union is commutative.

            :param other: the other state space (or wrapper) to combine with
            :type other: object
            :return: the unioned wrapper
            :rtype: PerturbationTheorySolver.StateSpaceWrapper
            """
            return self.simple_union(other)
        # def __mul__(self, other):
        #     return type(self)(self.space.intersection(other))
    class ProjectionOperatorWrapper:
        """
        Generates a symbolic form of a perturbation operator that
        either projects onto a degenerate space or projects it out
        """
        def __init__(self, space, complement=False):
            """
            **LLM Docstring**

            Set up a symbolic operator that, when applied to a state space, either restricts it to (or excludes) a fixed set of states -- used to represent projecting onto or away from a degenerate subspace.

            :param space: the fixed set of states to project onto/away from
            :type space: BasisStateSpace
            :param complement: whether this operator excludes `space` (`True`) rather than restricting to it (`False`)
            :type complement: bool
            :return: None
            :rtype: None
            """
            self.space = space
            self.complement = complement
        def get_transformed_space(self, other):
            """
            :param other:
            :type other: SelectionRuleStateSpace
            :return:
            :rtype:
            """
            if self.complement:
                return other.drop_states(self.space)
            else:
                return other.take_states(self.space)
    class ProjectedOperator:
        """
        Generates a symbolic form of an operator where
        an operator can first be applied, then unused terms projected
        out, before returning the state space
        """
        def __init__(self, projector, operator):
            """
            **LLM Docstring**

            Set up a symbolic operator that first applies an underlying operator's state-space transformation, then filters the result through a projector -- used to compose "apply this coupling, then keep/drop degenerate states" symbolically without materializing intermediate state spaces.

            :param projector: the projection operator to apply after `operator`
            :type projector: object
            :param operator: the underlying operator whose transformed space gets projected
            :type operator: object
            :return: None
            :rtype: None
            """
            self.proj = projector
            self.op = operator

        def get_transformed_space(self, other):
            """
            :param other:
            :type other: BasisStateSpace
            :return:
            :rtype:
            """
            a = self.op
            if (
                    isinstance(a, (int, np.integer, float, np.floating)) and a == 0
            ):
                return None

            wtf1 = self.op.get_transformed_space(other)
            contracted = self.proj.get_transformed_space(wtf1)
            return contracted

    def _could_be_a_space(self, test): # nasty checking code that I don't want to redupe all the time
        """
        **LLM Docstring**

        Heuristically check whether an object could plausibly represent a state space: an actual `BasisStateSpace`/`BasisMultiStateSpace`, a flat list of integer indices, or a nested list of integer excitation vectors.

        :param test: the object to check
        :type test: object
        :return: whether `test` looks like it could be a state space
        :rtype: bool
        """
        if isinstance(test, (BasisStateSpace, BasisMultiStateSpace)):
            return True
        else:
            try:
                lt = len(test)
            except TypeError:
                return False

            try:
                if lt > 0 and isinstance(test[0], (int, np.integer)):
                    return True
                else:
                    lt = len(test[0])
                    if lt > 0 and isinstance(test[0][0], (int, np.integer)):
                        return True
                    else:
                        return False
            except TypeError:
                return False
    def _could_be_rules(self, test):
        """
        **LLM Docstring**

        Heuristically check whether an object could plausibly represent a set of selection rules: `None` (no rules), an empty sequence, or a sequence whose first element is itself empty or starts with an integer.

        :param test: the object to check
        :type test: object
        :return: whether `test` looks like it could be a set of selection rules
        :rtype: bool
        """
        # selection rule options
        if test is None:
            return True

        try:
            lt = len(test)
        except TypeError:
            return False

        if lt == 0:
            return True

        try:
            lt = len(test[0])
            res = (
                    lt == 0
                    or isinstance(test[0][0], (int, np.integer))
            )
        except TypeError:
            return False

        return res
    def _could_be_a_prefilter(self, test):
        """
        **LLM Docstring**

        Heuristically check whether an object could plausibly represent a single prefilter specification: `None`, a plain state space (via `_could_be_a_space`), or a `(space, rules)` pair matching `_could_be_a_space`/`_could_be_rules` respectively.

        :param test: the object to check
        :type test: object
        :return: whether `test` looks like a valid prefilter specification
        :rtype: bool
        """
        return (
            test is None
            or self._could_be_a_space(test)
            or (
                    len(test) == 2 and
                    test[0] is None or self._could_be_a_space(test[0])
                    and self._could_be_rules(test[1])
            )
        )
    def _apply_transformation_with_filters(self, a, b, filter_space:BasisStateSpaceFilter, **opts):
        """
        **LLM Docstring**

        Apply a coupling operator's state-space transformation to `b` (interpreted relative to `a`), then run the result through any pre/post state-space filters attached to `filter_space`.

        :param a: the source/reference state space the transformation conceptually originates from
        :type a: object
        :param b: the operator/state-space object whose `get_transformed_space` should be applied
        :type b: object
        :param filter_space: the filter set (pre- and post-transformation filters) to apply to the result; if `None`, no filtering is applied
        :type filter_space: BasisStateSpaceFilter | None
        :param opts: extra options forwarded to the underlying transformation/filtering calls
        :type opts: dict
        :return: the (possibly filtered) transformed state space
        :rtype: object
        """

        if self._state_space_store is not None:
            b = self._state_space_store.coerce(b)
        if filter_space is not None:
            prefilters = filter_space.prefilters
            postfilters = filter_space.postfilters
        else:
            prefilters = None
            postfilters = None

        if prefilters is not None:
            # this means we are able to filter _before_ we apply the selection rules
            # which we do by iteratively breaking the space to transform up into
            # chunks that match the prefilters and applying the associated selection
            # rules

            new = None # the new space we will slowly build
            b_remainder = b #.as_sorted()

            for n, (filter_space, filter_rules) in enumerate(prefilters):

                # take the intersection/remainder
                if filter_space is None:
                    b = b_remainder
                    b_remainder = None
                elif len(b_remainder) > 0:
                    b = b_remainder.intersection(filter_space)
                    if n < len(prefilters): # just a cheap opt...
                        b_remainder = b_remainder.difference(filter_space)
                else:
                    b = b_remainder

                if len(b) > 0:
                    # if filter_rules is not None and len(filter_rules) == 0:
                    #     new = b.to_single().take_unique()
                    # else:
                    if new is None:
                        new = a.get_transformed_space(b, rules=filter_rules, **opts)
                    else:
                        new = new.union(
                            a.get_transformed_space(b, rules=filter_rules, **opts)
                        )

        else:  # no prefilters
            new = a.get_transformed_space(b, **opts)

        if new is not None: # possible to have None if we had to do no work
            if postfilters is not None:
                for filter in postfilters:
                    if isinstance(new, BasisMultiStateSpace):
                        new = new.map(filter.apply)
                    else:
                        new = filter.apply(new)

        return new

    def _get_new_coupled_space(self, a, b, spaces=None, ret_space=True, filter_space=None):
        """
        A symbolic version of the dot product appropriate for getting
        transformed state spaces under the operation of a on b

        :param a:
        :type a:
        :param b:
        :type b: SelectionRuleStateSpace
        :param spaces:
        :type spaces: the set of operators to which we can assign transformations
        :return:
        :rtype: StateSpaceWrapper
        """
        if spaces is None:
            raise ValueError("...spaces shouldn't be None")

        if isinstance(b, self.StateSpaceWrapper):
            b = b.space

        if isinstance(a, self.ProjectedOperator):
            op = a.op
        else:
            op = a
        if (
                isinstance(op, (int, np.integer, float, np.floating)) and op == 0
                or isinstance(b, (int, np.integer, float, np.floating)) and b == 0
                or op is None
                or b is None
        ):
            return 0

        logger = self.logger
        if isinstance(a, self.StateSpaceWrapper):
            raise NotImplementedError("we shouldn't be here")
            new = a * b
        elif isinstance(a, self.ProjectionOperatorWrapper):
            raise NotImplementedError("we shouldn't be here")
            # uhhh...
            new = a.get_transformed_space(b, parallelizer=self.parallelizer, logger=logger)
        elif isinstance(a, (self.ProjectedOperator, Representation)):
            cur = spaces[op] #type: SelectionRuleStateSpace
            proj = None if not isinstance(a, self.ProjectedOperator) else a.proj
            if cur is None:

                new = self._apply_transformation_with_filters(
                    a, b, filter_space,
                    track_excitations= not self.memory_constrained,
                    parallelizer=self.parallelizer, logger=logger
                )

                # we track not only the output SelectionRuleStateSpace
                # but also which projection operators have been applied
                # so that we can make sure we calculate any pieces that
                # need to be calculated

                spaces[op] = (
                    {proj:b},
                    new
                )
                # reduce to a single space to feed to the next round
                # of terms

                if new is not None:
                    new = new.to_single().take_unique()
            else:
                projections, cur = cur
                # figure out what stuff we've already calculated
                rep_space = projections[None] if None in projections else None
                if proj is not None:
                    sub_rep = projections[proj] if proj in projections else None
                    if sub_rep is not None:
                        if rep_space is not None:
                            rep_space = rep_space.union(sub_rep, track_excitations=False)
                        else:
                            rep_space = sub_rep

                if rep_space is None:
                    # means we can't determine which parts we have and have not calculated
                    # so we calculate everything and associate it to proj

                    new = self._apply_transformation_with_filters(
                        a, b, filter_space,
                        track_excitations=not self.memory_constrained,
                        parallelizer=self.parallelizer, logger=logger
                    )

                    cur = cur.union(new)
                    projections[proj] = b
                    spaces[op] = (projections, cur)
                    # reduce to a single space to feed to the next round
                    # of terms
                    new = new.to_single(track_excitations=not self.memory_constrained).take_unique()
                else:
                    # means we've potentially calculated some of this already,
                    # so we figure out what parts we've already calculated in this
                    # projected space (rep_space is the current space of the representations)
                    diffs = b.difference(rep_space)
                    # if diffs.full_basis is None:
                    #     raise ValueError(diffs.full_basis, b.full_basis)
                    if len(diffs) > 0:
                        # raise Exception(projections, rep_space, diffs)
                        # we have an initial space we've already transformed, so we
                        # make sure not to recompute that
                        b_sels = SelectionRuleStateSpace(b, [], ignore_shapes=True, # just some type fuckery
                                                         changes=[] if SelectionRuleStateSpace.track_change_positions else None
                                                         )
                        existing = cur.intersection(b_sels, handle_subspaces=False)
                        # and now we do extra transformations where we need to


                        if len(diffs) > 0:
                            new_new = self._apply_transformation_with_filters(
                                a, diffs, filter_space,
                                track_excitations=not self.memory_constrained,
                                parallelizer=self.parallelizer, logger=logger
                            )

                            if new_new is None:
                                new = b.to_single().take_unique()
                            else:
                                # next we add the new stuff to the cache
                                cur = cur.union(new_new)
                                projections[proj] = rep_space.union(diffs)
                                spaces[op] = (projections, cur)

                                # TODO: find a way to make this not cause memory spikes...
                                if ret_space:
                                    new = existing.union(new_new).to_single(track_excitations=not self.memory_constrained).take_unique()
                                else:
                                    new = b.to_single().take_unique()
                        else:
                            new = b.to_single().take_unique()

                    else:
                        # means we already calculated everything
                        # so we don't need to worry about this
                        if cur is not None:
                            b_sels = SelectionRuleStateSpace(b, [],
                                                             changes=[] if SelectionRuleStateSpace.track_change_positions else None,
                                                             ignore_shapes=True) # just some type fuckery
                            new = cur.intersection(b_sels, handle_subspaces=False)
                            new = new.to_single().take_unique()
                        else:
                            new = None

        else:
            raise TypeError("don't know what to do with {} and {}".format(a, b))

        return self.StateSpaceWrapper(new)
    def _reduce_new_coupled_space(self, *terms,
                                  spaces=None,
                                  ret_space=True,
                                  filter_space=None
                                  ):
        """
        Reduces through `_get_new_coupled_space` from right to left
        :param terms:
        :type terms: Iterable[SelectionRuleStateSpace]
        :param spaces:
        :type spaces: dict | None
        :return:
        :rtype:
        """
        import functools
        return functools.reduce(
            lambda a, b:self._get_new_coupled_space(b, a, spaces, ret_space=ret_space, filter_space=filter_space),
            reversed(terms[:-1]),
            terms[-1]
        )
    def get_coupled_space(self,
                            input_state_space,
                            degenerate_space,
                            use_second_deg,
                            allow_PT_degs=True,
                            wavefunction_terms=None,
                            spaces=None,
                            property_filter=None,
                            filter_spaces=None
                            ):
        """
        Applies the VPT equations semi-symbolically, dispatching based on how many
        degeneracies we need to handle

        :return:
        :rtype:
        """

        if not allow_PT_degs:
            spaces = self.get_nondeg_coupled_space(
                input_state_space,
                degenerate_space,
                spaces=spaces,
                wavefunction_terms=wavefunction_terms,
                property_filter=property_filter,
                filter_spaces=filter_spaces
            )
        else:
            raise NotImplementedError("True degeneracy handling was scrubbed for lack of utility")

        # raise Exception(spaces)
        return spaces

    def get_nondeg_coupled_space(self,
                                 input_state_space,
                                 degenerate_space=None,
                                 spaces=None,
                                 wavefunction_terms=None,
                                 property_filter=None,
                                 filter_spaces=None
                                 ):
        """
        Applies the non-degenerate equations in semi-symbolic form to determine
        which states needs to be calculated.
        This will always be the initial input to a calculation and then
        certain additional states may be calculated on the fly if they are needed to handle
        truly degenerate stuff.
        The only difference there will be to add on

        :param input_state_space:
        :type input_state_space: BasisStateSpace
        :param degenerate_space:
        :type degenerate_space: BasisStateSpace
        :param spaces:
        :type spaces:
        :param wavefunction_terms: which terms to include when calculating corrections
        :type wavefunction_terms: None | Iterable[Iterable[int]]
        :param property_filter: a set of states and selection rules to allow for being targeted in state to calculate
        :type property_filter: None | Iterable[Iterable[int], Iterable[Iterable[int]]]
        :return:
        :rtype:
        """

        # holder for perts to map to their stored states
        spaces = {h:None for h in self.perts} if spaces is None else spaces
        # final state spaces for each energy, corr, overlap, etc.
        input_wrapper = self.StateSpaceWrapper(input_state_space)

        order = self.state_space_iterations if self.state_space_iterations is not None else self.order
        order = order+1
        E = [None]*order
        E[0] = input_wrapper
        corrs = [None]*order
        corrs[0] = input_wrapper

        D = degenerate_space
        if D is None:
            D = input_state_space
        # pi = self.ProjectionOperatorWrapper(D, complement=True)
        # piD = self.ProjectionOperatorWrapper(D, complement=False)

        dot = lambda *terms, spaces=spaces,ret_space=True,filter_space=None: self._reduce_new_coupled_space(*terms,
                                                                                          spaces=spaces,
                                                                                          ret_space=ret_space,
                                                                                          filter_space=filter_space
                                                                                          )

        H = self.PastIndexableTuple(self.perts)

        if filter_spaces is None:
            filter_spaces = {}
        if property_filter is None:
            # This is intentionally written to parallel the non-degenerate VPT equations
            for k in range(1, order):
                ####################################################################################################
                #                                       *IMPLEMENTATION NOTE*
                # The states for the energy end up being totally subsumed in the states
                # for the corrections, so we just leave this part off
                ####################################################################################################
                #         En^(k) = <n^(0)|H^(k)|n^(0)> + sum(<n^(0)|H^(k-i)|n^(i)> - E^(k-i)<n^(0)|n^(i)>, i=1...k-1)
                # but we drop the shift part because that doesn't affect the state space at all
                # E[k] = sum(
                #             dot(self.ProjectedOperator(piD, H[k - i]), corrs[i])
                #             for i in range(0, k)
                #             )

                #   <n^(0)|n^(k)> = -1/2 sum(<n^(i)|n^(k-i)>, i=1...k-1)
                #         |n^(k)> = sum(Pi_n (En^(k-i) - H^(k-i)) |n^(i)>, i=0...k-1) + <n^(0)|n^(k)> |n^(0)>
                # but we drop the energy and overlap parts of this because they don't affect the overall state space

                self.logger.log_print(
                    'getting states for ' +
                        '+'.join('H({})|n({})>'.format(k-i, i)
                                 for i in range(0, k)
                                 if (
                                         not isinstance(H[k - i], (int, np.integer))
                                         and (wavefunction_terms is None or (k-i, i) in wavefunction_terms)
                                 )
                                 )
                    )

                with self.logger.block(tag='getting states for order {k}'.format(k=k)):
                    corrs[k] = sum(corrs[i] for i in range(0, k)) # this all in here from energies
                    for i in range(0, k):
                        if not isinstance(H[k - i], (int, np.integer)):
                            if wavefunction_terms is None or (k-i, i) in wavefunction_terms:
                                self.logger.log_print('H({a})|n({b})>', a=k - i, b=i)
                                if k < order-1:
                                    corrs[k] += dot(H[k - i], corrs[i],
                                                    filter_space=filter_spaces.get((k-i, i), None)
                                                    )
                                else:
                                    dot(H[k - i], corrs[i], ret_space=False,
                                                    filter_space=filter_spaces.get((k-i, i), None)
                                                    )
        else:
            raise NotImplementedError("property filters not here yet")

        return spaces

    #endregion

    #region Apply Equations

    def get_corrections(self,
                        non_zero_cutoff=None,
                        handle_strong_couplings=None,
                        check_overlap=None
                        ):
        """
        Applies the perturbation theory equations to obtain
        corrections to the wave functions and energies

        :return:
        :rtype:
        """
        # We use the iterative equations
        #            En^(k) = <n^(0)|H^(k)|n^(0)> + sum(<n^(0)|H^(k-i)|n^(i)> - E^(k-i)<n^(0)|n^(i)>, i=1...k-1)
        #     <n^(0)|n^(k)> = -1/2 sum(<n^(i)|n^(k-i)>, i=1...k-1)
        #           |n^(k)> = sum(Pi_n (En^(k-i) - H^(k-i)) |n^(i)>, i=1...k-1) + <n^(0)|n^(k)> |n^(0)>
        #  where Pi_n is the perturbation operator [1/(E_m-E_n) for m!=n]

        perturbations = self.representations
        states = self.states
        order = self.order

        flat_total_space = self.flat_total_space
        N = self.total_space_dim

        checkpointer = self.checkpointer
        logger = self.logger

        degenerate_states = self.degenerate_spaces
        if handle_strong_couplings is None:
            handle_strong_couplings = self.handle_strong_couplings

        if non_zero_cutoff is None:
            non_zero_cutoff = Settings.non_zero_cutoff

        return self._get_corrections(
            perturbations,
            states,
            order,
            flat_total_space,
            N,
            checkpointer,
            logger,
            degenerate_states,
            handle_strong_couplings=handle_strong_couplings,
            non_zero_cutoff=non_zero_cutoff
        )

    @property
    def high_frequency_modes(self):
        """
        **LLM Docstring**

        The indices of the vibrational modes whose fundamental transition frequency exceeds `self.low_frequency_mode_cutoff`, used as the default set of modes considered for strong-coupling/degeneracy testing.

        :return: the high-frequency mode indices
        :rtype: list[int]
        """
        fundamentals = BasisStateSpace.from_quanta(self.total_state_space.basis, [0, 1])
        where = self.flat_total_space.find(fundamentals)
        test_freqs = self.zero_order_energies[where[1:]] - self.zero_order_energies[where[0]]
        test_modes = np.where(test_freqs > self.low_frequency_mode_cutoff)
        if len(test_modes) > 0:
            test_modes = list(test_modes[0])
        return test_modes

    default_strong_coupling_threshold = .3
    @staticmethod
    def _group_matrix_values(matrix, rows, columns):
        value = matrix[rows, columns]
        if sp.issparse(value):value = value.toarray()
        return np.asarray(value).reshape(-1)

    @staticmethod
    def _group_edge_positions(matrix, predicate):
        """Row-major positions selected from stored CSR entries, not a dense mask."""
        keep = np.flatnonzero(predicate(matrix.data))
        rows = np.searchsorted(matrix.indptr, keep, side='right') - 1
        return rows, matrix.indices[keep]

    @staticmethod
    def _threshold_group_graph(matrix, threshold):
        graph = matrix.copy()
        graph.data = (graph.data > threshold).astype(int)
        graph.eliminate_zeros()
        return graph

    @staticmethod
    def _set_group_graph_value(matrix, row, column, value):
        # Retain the legacy integer adjacency assignment semantics.
        value = np.asarray(value,dtype=matrix.dtype).item()
        start,stop = matrix.indptr[row:row+2]
        offset = np.searchsorted(matrix.indices[start:stop],column)
        if start+offset < stop and matrix.indices[start+offset] == column:
            matrix.data[start+offset] = value
            return matrix
        if value == 0:return matrix
        return matrix + sp.csr_matrix(([value],([row],[column])),shape=matrix.shape,dtype=matrix.dtype)

    class _GroupCorrectionView:
        """Compute group weights without making (order, group, all targets) arrays."""
        def __init__(self, matrices, group_columns, state_rows, block_size):
            self.matrices = matrices
            self.columns = np.asarray(group_columns)
            self.state_rows = np.asarray(state_rows)
            self.block_size = block_size

        @staticmethod
        def _take_row(matrix, row, columns):
            if isinstance(matrix, SparseArray):
                data = getattr(matrix,'data',None)
                if sp.issparse(data) and data.format == 'csr' and data.has_canonical_format:
                    start, stop = data.indptr[row:row+2]
                    indices, values = data.indices[start:stop], data.data[start:stop]
                    positions = np.searchsorted(indices, columns)
                    valid = positions < len(indices)
                    result = np.zeros(len(columns), dtype=values.dtype)
                    found = np.flatnonzero(valid)
                    found = found[indices[positions[found]] == columns[found]]
                    result[found] = values[positions[found]]
                    return result
                # Other supported sparse formats use bounded stored-entry iteration.
                result = np.zeros(len(columns), dtype=float)
                sorting = np.argsort(columns)
                ordered = columns[sorting]
                for values, (rows, cols) in matrix.iter_nonzero_blocks():
                    selected = rows == row
                    vals, cols = values[selected], cols[selected]
                    positions = np.searchsorted(ordered, cols)
                    valid = positions < len(ordered)
                    found = np.flatnonzero(valid)
                    found = found[ordered[positions[found]] == cols[found]]
                    np.add.at(result, sorting[positions[found]], vals[found])
                return result
            return np.asarray(matrix[row, columns]).reshape(-1)

        def maxima(self):
            """Legacy sorting statistic, needed only when threshold is None."""
            maximum = np.zeros(len(self.columns))
            for matrix in self.matrices:
                for row in range(matrix.shape[0]):
                    for start in range(0, len(self.columns), self.block_size):
                        stop = min(start+self.block_size, len(self.columns))
                        values = np.abs(self._take_row(matrix, row, self.columns[start:stop]))
                        values[self.state_rows[start:stop] == row] = 0
                        np.maximum(maximum[start:stop], values, out=maximum[start:stop])
            return maximum

        def symmetric_strengths(self, found_pos):
            """Same ordered sum of absolute coefficients and reciprocal directions."""
            count = len(self.columns)
            values, columns, indptr, nnz = [], [], [0], 0
            for local, source in enumerate(self.state_rows):
                if found_pos[local]:
                    for start in range(0, count, self.block_size):
                        stop = min(start+self.block_size, count)
                        total = None
                        for matrix in self.matrices:
                            part = np.abs(self._take_row(matrix, int(source), self.columns[start:stop]))
                            total = part if total is None else total + part
                        if total is None:total = np.zeros(stop-start)
                        if start <= local < stop:total[local-start] = 0
                        keep = np.flatnonzero(total != 0)
                        if len(keep):
                            values.append(total[keep]);columns.append(keep+start)
                            nnz += len(keep)
                indptr.append(nnz)
            data = np.concatenate(values) if values else np.empty(0)
            indices = np.concatenate(columns) if columns else np.empty(0,dtype=np.intp)
            matrix = sp.csr_matrix((data,indices,np.asarray(indptr)),shape=(count,count))
            matrix = (matrix.T + matrix).tocsr()
            matrix.sum_duplicates();matrix.sort_indices();matrix.eliminate_zeros()
            return matrix

    class _GroupPairRules:
        """Implicit incompatible-pair relation; no group-by-group-by-mode tensor."""
        def __init__(self, excitations, energies, energy_cutoff, block_size):
            self.excitations = excitations
            self.energies = energies
            self.cutoff = energy_cutoff
            self.block_size = block_size

        def bad(self, rows, columns):
            rows, columns = np.asarray(rows), np.asarray(columns)
            if self.cutoff is not None and self.energies is not None:
                return np.abs(self.energies[rows] - self.energies[columns]) > self.cutoff
            left, right = self.excitations[rows], self.excitations[columns]
            return (np.sum(np.abs(left-right),axis=-1) == 1) & (np.sum(np.abs(left+right),axis=-1) > 0)

        def row(self, row, difference=None, zeroed=False):
            pieces = []
            for start in range(0,len(self.excitations),self.block_size):
                columns = np.arange(start,min(start+self.block_size,len(self.excitations)))
                if difference is None:mask = self.bad(np.full(len(columns),row),columns)
                else:
                    distances = (np.zeros(len(columns),dtype=int) if zeroed else
                        np.sum(np.abs(self.excitations[row]-self.excitations[columns]),axis=-1))
                    mask = distances == difference
                pieces.append(columns[mask])
            return np.concatenate(pieces) if pieces else np.empty(0,dtype=int)

        def any_bad(self):
            if self.cutoff is not None and self.energies is not None:
                values = self.energies[~np.isnan(self.energies)]
                return bool(len(values) and np.max(values)-np.min(values) > self.cutoff)
            return any(len(self.row(row)) for row in range(len(self.excitations)))

        def remove_bad_edges(self, matrix):
            matrix = matrix.copy()
            for start in range(0,len(matrix.data),self.block_size):
                stop = min(start+self.block_size,len(matrix.data))
                positions = np.arange(start,stop)
                rows = np.searchsorted(matrix.indptr,positions,side='right')-1
                bad = self.bad(rows,matrix.indices[start:stop])
                matrix.data[start:stop][bad] = 0
            matrix.eliminate_zeros()
            return matrix

        def components_conflict(self, labels):
            if self.cutoff is not None and self.energies is not None:
                count = int(np.max(labels))+1 if len(labels) else 0
                minimum, maximum = np.full(count,np.inf),np.full(count,-np.inf)
                valid = ~np.isnan(self.energies)
                np.minimum.at(minimum,labels[valid],self.energies[valid])
                np.maximum.at(maximum,labels[valid],self.energies[valid])
                with np.errstate(invalid='ignore'):
                    return bool(np.any(maximum-minimum > self.cutoff))
            for row in range(len(labels)):
                bad = self.row(row)
                if np.any(labels[bad] == labels[row]):return True
            return False

    @staticmethod
    def _get_sparse_groups(matrix):
        # Explicit stored zeros are not edges in the legacy dense graph.
        matrix.eliminate_zeros()
        _, labels = sp.csgraph.connected_components(matrix,return_labels=True)
        from McUtils import Numputils as nput
        (_, groups), _ = nput.group_by(np.arange(matrix.shape[0]),labels)
        return groups, labels

    @staticmethod
    def _filter_group_vector_rows(excitations, filters, block_size):
        """Bounded equivalent of the legacy vector-rule broadcast/where branch."""
        filters = np.asarray(filters)
        if filters.ndim == 1:filters = filters[None]
        result = []
        for row in range(len(excitations)):
            found = False
            for start in range(0,len(excitations),block_size):
                difference = excitations[row]-excitations[start:start+block_size]
                for rule in filters:
                    if np.any(np.all(difference == rule,axis=-1) | np.all(difference == -rule,axis=-1)):
                        found = True;break
                if found:break
            if found:result.append(row)
        return np.asarray(result,dtype=int)

    @staticmethod
    def _filter_group_mode_differences(group, limits, block_size):
        """Preserve the original row-major greedy removal order in bounded blocks."""
        excitations = group.excitations
        limits = np.asarray(limits)
        tested = np.flatnonzero(limits > -1)
        dead = set()
        if len(tested):
            for row in range(len(excitations)):
                for start in range(row,len(excitations),block_size):
                    columns = np.arange(start,min(start+block_size,len(excitations)))
                    differences = np.abs(excitations[row,tested]-excitations[columns][:,tested])
                    for column in columns[np.any(differences > limits[tested],axis=-1)]:
                        if row in dead or column in dead:continue
                        if np.sum(excitations[row,tested]) > np.sum(excitations[column,tested]):dead.add(row)
                        else:dead.add(int(column))
        keep = np.setdiff1d(np.arange(len(excitations)),list(dead))
        return group.take_subspace(keep) if len(keep) else None

    @classmethod
    def filter_degenerate_group(cls, group, states=None, basis=None, corrections=None, threshold=None, energy_cutoff=None, frequencies=None, energies=None, decoupling_overide=10, maximize_groups=True, target_modes=None, max_mode_differences=None, threshold_step_size=0.1, extra_groups=None, uncoupled_states=None, state_diff_filters=1, logger=None, block_size=65536, state_formatting="quanta"):
        if not isinstance(block_size,(int,np.integer)) or block_size < 1:raise ValueError('block_size must be a positive integer')
        if len(group) < 2:
            return group
        if uncoupled_states is not None:
            if corrections is not None:
                if hasattr(corrections, 'total_basis'):
                    states = corrections.states
                    basis = corrections.total_basis
                elif hasattr(corrections, 'matrices'):
                    basis = corrections.full_basis
                    states = corrections.initial_states
            smol = isinstance(group, BasisStateSpace)
            if smol:
                group = [group]
            new_groups = []
            if nput.is_int(uncoupled_states):
                uncoupled_states = BasisStateSpace.from_quanta(basis.basis if basis is not None else group[0].basis if len(group) > 0 else HarmonicOscillatorProductBasis(len(states[0])), uncoupled_states)
            for g in group:
                uncoupled_pos = g.find(uncoupled_states, dtype=int, missing_val=-1)
                found_pos = np.where(uncoupled_pos > -1)[0]
                if len(found_pos) > 0:
                    if logger is not None:
                        logger.log_print('WARNING: states {g0} were in degenerate group {g2}, coupling has been removed', g0=uncoupled_states.excitations, g2=g.excitations,
                            preformatter=None if state_formatting == "vector" else logger.preformat_keys({key: lambda s: StateMaker.format_states(s, state_formatting, join=True) for key in ['g0', 'g2']}))
                    mask = np.setdiff1d(np.arange(len(g)), found_pos)
                    g = g.take_subspace(mask)
                new_groups.append(g)
            if smol:
                new_groups = new_groups[0]
            group = new_groups
        if target_modes is not None:
            sub = group.take_subdimensions(target_modes)
            inds = sub.indices
            exc = sub.excitations
        else:
            inds = group.indices
            exc = group.excitations
        if corrections is not None:
            if isinstance(corrections, np.ndarray):
                state_inds = basis.find(group)
                group_inds = state_inds
                found_pos = state_inds < len(states)
            elif hasattr(corrections, 'total_basis'):
                states = corrections.states
                basis = corrections.total_basis
                group_inds = basis.find(group)
                state_inds = states.find(group, dtype=int, missing_val=len(states))
                found_pos = state_inds < len(states)
                corrections = cls._GroupCorrectionView(corrections.wfn_corrections[1:], group_inds, state_inds, block_size)
                sorting = (np.arange(len(group_inds)) if threshold is not None else
                    np.argsort(-corrections.maxima() - 1*found_pos))
            elif hasattr(corrections, 'matrices'):
                basis = corrections.full_basis
                states = corrections.initial_states
                group_inds = basis.find(group)
                state_inds = states.find(group, dtype=int, missing_val=len(states))
                found_pos = state_inds < len(states)
                corrections = cls._GroupCorrectionView(corrections.matrices[1:], group_inds, state_inds, block_size)
                sorting = (np.arange(len(group_inds)) if threshold is not None else
                    np.argsort(-corrections.maxima() - 1*found_pos))
            elif basis is not None and states is not None:
                group_inds = basis.find(group)
                state_inds = states.find(group, dtype=int, missing_val=len(states))
                found_pos = state_inds < len(states)
                corrections = cls._GroupCorrectionView([corrections], group_inds, state_inds, block_size)
                sorting = (np.arange(len(group_inds)) if threshold is not None else
                    np.argsort(-corrections.maxima() - 1*found_pos))
            else:
                corrections = None
        if corrections is None:
            sorting = np.argsort(inds)
            corrections = None
            if basis is not None:
                group_inds = basis.find(group)
            else:
                group_inds = np.arange(len(group))
        local_energies = None if energies is None else np.asarray(energies[group_inds])
        if local_energies is None and frequencies is not None and basis is not None:
            # Evaluate only this group's energies, never the full basis table.
            local_energies = np.tensordot(group.excitations + 1/2, frequencies, axes=[-1,0])
        exc = exc[sorting]
        rules = cls._GroupPairRules(exc, local_energies, energy_cutoff, block_size)
        elims = []
        if threshold is not None and corrections is not None:
            e_vec = None if energy_cutoff is None else local_energies
            if rules.any_bad():
                N = len(group_inds)
                if isinstance(corrections, cls._GroupCorrectionView):
                    base_mat = corrections.symmetric_strengths(found_pos)
                else:
                    # Legacy ndarray input is already a group-by-target view.
                    selected = sp.csr_matrix(corrections[:, state_inds[found_pos]])
                    coordinates = selected.tocoo()
                    base_mat = sp.coo_matrix((coordinates.data,(coordinates.row,np.flatnonzero(found_pos)[coordinates.col])),shape=(N,N)).tocsr()
                    base_mat = base_mat.T + base_mat
                    base_mat.setdiag(0);base_mat.eliminate_zeros()
                base_mat = rules.remove_bad_edges(base_mat)
                get_groups = cls._get_sparse_groups
                corr_mat = cls._threshold_group_graph(base_mat, threshold)
                groups, A = get_groups(corr_mat)
                target_threshold = threshold
                step_size = threshold_step_size
                while rules.components_conflict(A) and decoupling_overide > threshold:
                    threshold = threshold + step_size
                    corr_mat = cls._threshold_group_graph(base_mat, threshold)
                    groups, A = get_groups(corr_mat)
                if maximize_groups or (decoupling_overide > 0 and decoupling_overide < threshold):
                    base_thresh = target_threshold
                    if not maximize_groups:
                        base_thresh = decoupling_overide
                    dropped_pos = cls._group_edge_positions(base_mat, lambda values: (values > base_thresh) & (values < threshold))
                    if len(dropped_pos) > 0:
                        sorting = np.argsort(-cls._group_matrix_values(base_mat, *dropped_pos))
                        dropped_pos = tuple((d[sorting] for d in dropped_pos))
                        for i, j in zip(*dropped_pos):
                            if e_vec is not None:
                                g_i = None
                                g_j = None
                                for n, g in enumerate(groups):
                                    ix, jx = np.searchsorted(g, [i, j])
                                    if ix < len(g) and g[ix] == i:
                                        g_i = [n, g]
                                    if jx < len(g) and g[jx] == j:
                                        g_j = [n, g]
                                    if g_i is not None and g_j is not None:
                                        break
                                else:
                                    if g_j is not None:
                                        g_i, g_j = (g_j, g_i)
                                        i, j = (j, i)
                                    if g_i is not None:
                                        e_g = e_vec[g_i[1]]
                                        w = max(abs(e_vec[i] - np.min(e_g)), abs(e_vec[i] - np.max(e_g)))
                                        if w < energy_cutoff:
                                            groups[g_i[0]] = np.sort(np.concatenate([g_i[1], [j]]))
                                        continue
                            corr_mat = cls._set_group_graph_value(corr_mat, i, j, base_mat[i, j])
                            _, A = get_groups(corr_mat)
                            if rules.components_conflict(A):
                                corr_mat = cls._set_group_graph_value(corr_mat, i, j, 0.0)
                            else:
                                groups = _
                    if decoupling_overide > 0 and decoupling_overide < threshold or extra_groups is not None:
                        if (decoupling_overide > 0 and decoupling_overide < threshold) is not None:
                            need_pos = cls._group_edge_positions(base_mat, lambda values: (values > decoupling_overide) & (values < threshold))
                        else:
                            need_pos = ()
                        if extra_groups is not None:
                            for g in extra_groups:
                                try:
                                    new_inds = group.find(g)
                                except IndexError:
                                    pass
                                else:
                                    new_inds = list(itertools.combinations(new_inds, 2))
                                    new_inds = ([x[0] for x in new_inds], [x[1] for x in new_inds])
                                    if len(need_pos) > 0:
                                        need_pos = (np.concatenate([need_pos[0], new_inds[0]]), np.concatenate([need_pos[1], new_inds[1]]))
                                    else:
                                        need_pos = new_inds
                        if len(need_pos) > 0:
                            sorting = np.argsort(-cls._group_matrix_values(base_mat, *need_pos))
                            need_pos = tuple((d[sorting] for d in need_pos))
                            for i, j in zip(*need_pos):
                                g_i = None
                                g_j = None
                                for n, g in enumerate(groups):
                                    ix, jx = np.searchsorted(g, [i, j])
                                    if ix < len(g) and g[ix] == i:
                                        g_i = [n, g]
                                    if jx < len(g) and g[jx] == j:
                                        g_j = [n, g]
                                    if g_i is not None and g_j is not None:
                                        if g_i[0] == g_j[0]:
                                            break
                                        i_couplings = cls._group_matrix_values(base_mat, i, g_i[1])
                                        j_couplings = cls._group_matrix_values(base_mat, j, g_j[1])
                                        if np.linalg.norm(j_couplings) > np.linalg.norm(i_couplings):
                                            g_i, g_j = (g_j, g_i)
                                            i, j = (j, i)
                                        groups[g_i[0]] = np.sort(np.concatenate([g_i[1], [j]]))
                                        groups[g_j[0]] = np.setdiff1d(g_j[1], [j])
                                        break
                                else:
                                    if g_j is not None:
                                        g_i, g_j = (g_j, g_i)
                                        i, j = (j, i)
                                    if g_i is not None:
                                        groups[g_i[0]] = np.sort(np.concatenate([g_i[1], [j]]))
                                        continue
                                    else:
                                        groups.append(np.array([i, j]))
                group = [group.take_subspace(g) for g in groups if len(g) > 1]
        else:
            zeroed_rows = np.zeros(len(exc),dtype=bool)
            global_bad_rows = None
            if energy_cutoff is not None and local_energies is not None:
                # The legacy branch uses the first coordinate of global where().
                global_bad_rows = np.asarray([i for i in range(len(exc)) if len(rules.row(i))],dtype=int)
            elif state_diff_filters is not None and not nput.is_int(state_diff_filters):
                global_bad_rows = cls._filter_group_vector_rows(exc,state_diff_filters,block_size)
            for n in np.arange(len(exc)):
                if global_bad_rows is not None:bad_pos = global_bad_rows
                elif state_diff_filters is None:bad_pos = np.empty(0,dtype=int)
                else:bad_pos = rules.row(n,difference=state_diff_filters,zeroed=zeroed_rows[n])
                if len(bad_pos):
                    kills = bad_pos[bad_pos > n]
                    elims.extend(sorting[kills])
                    zeroed_rows[bad_pos] = True
            if len(elims) > 0:
                mask = np.setdiff1d(np.arange(len(exc)), elims)
                group = group.take_subspace(mask)
        if max_mode_differences is not None:
            smol = isinstance(group, BasisStateSpace)
            if smol:
                group = [group]
            new_groups = []
            for g in group:
                g = cls._filter_group_mode_differences(g,max_mode_differences,block_size)
                if g is None:return []
                new_groups.append(g)
            if smol:
                new_groups = new_groups[0]
            group = new_groups
        if isinstance(group, BasisStateSpace):
            return group
        else:
            return [g for g in group if len(g.excitations) > 1]

    def identify_strong_couplings(self, corrs):
        """
        **LLM Docstring**

        Find the strongly-coupled state pairs among the solved corrections, delegating to whichever configured `degeneracy_handlers` spec exposes a `wfc_threshold` (or, if none do, a freshly built default `'auto'` spec).

        :param corrs: the perturbation-theory corrections to search for strong couplings in
        :type corrs: PerturbationTheoryCorrections
        :return: `(strong_couplings, threshold)` -- the identified strong-coupling data and the threshold used to find it
        :rtype: tuple
        """
        for spec in self.degeneracy_handlers:
            if hasattr(spec, 'wfc_threshold'):
                return spec.identify_strong_couplings(
                    self,
                    corrs
                ), spec.wfc_threshold
        else:
            spec = DegeneracySpec.from_spec('auto')
            return spec.identify_strong_couplings(
                self,
                corrs
            ), spec.wfc_threshold

    def construct_strong_coupling_spaces(self, spec, sc, corrs, states, threshold):
        """
        **LLM Docstring**

        Build the degenerate state groups implied by a set of identified strong couplings, and, if the spec allows extending the state space (`spec.extend_spaces`), extend the solver's states/representations to cover any newly implicated states before returning the resulting (possibly larger) state space and representation data.

        :param spec: the degeneracy spec (typically a `StronglyCoupledDegeneracySpec`) whose group filter/extension settings govern this construction
        :type spec: DegeneracySpec
        :param sc: the raw strong-coupling data (per-state, per-order coupled-state lists) to build groups from
        :type sc: dict
        :param corrs: the perturbation-theory corrections the coupling data came from, used to collapse it into per-state coupled-state spaces
        :type corrs: PerturbationTheoryCorrections
        :param states: the current target state space
        :type states: BasisStateSpace
        :param threshold: the coupling-strength threshold to build the group filter with
        :type threshold: float
        :return: `(degenerate_states, (states, perturbations, flat_total_space, N))` -- the identified degenerate groups, and the (possibly extended) states/representations/flat space/dimension to use going forward
        :rtype: tuple
        """
        groups = self._strong_coupling_groups(spec, sc, corrs, threshold)
        return groups, self._extend_strong_coupling_spaces(spec, groups, states)

    def _strong_coupling_groups(self, spec, sc, corrs, threshold):
        """Finish native group filtering while discovery coefficients are live."""
        sc = corrs.collapse_strong_couplings(sc)
        # self.logger.log_print(
        #     corrs.states.excitations,
        #     message_prepper=lambda a: str(np.array(a)).splitlines()
        # )
        with self.logger.block(tag='unpruned spaces states', log_level=self.logger.LogLevel.Never):
            for s in sc.values():
                self.logger.log_print(
                    s,
                    message_prepper=lambda a: StateMaker.format_states(a, self.state_formatting),
                    log_level=self.logger.LogLevel.Never
                )

        group_filter = spec.get_degenerate_group_filter(self, corrs=corrs, threshold=threshold)
        degenerate_states = DegenerateMultiStateSpace.from_spec(
            spec,
            couplings=sc,
            solver=self,
            group_filter=group_filter,
            default_group_filter=(partial(self.filter_degenerate_group, block_size=self.resonance_block_size, state_formatting=self.state_formatting)
                if self.resonance_memory_constrained else None),
            log_groups=True
        )

        return degenerate_states

    def _extend_strong_coupling_spaces(self, spec, degenerate_states, states):
        """Extend representations after finalized groups no longer need amplitudes."""
        if spec.extend_spaces:
            missing = degenerate_states.to_single().difference(states)
            if len(missing) > 0:
                with self.logger.block(tag='new states'):
                    self.logger.log_print(
                        missing.excitations,
                        message_prepper=lambda a: StateMaker.format_states(a, self.state_formatting)
                    )
                extended_spaces = self.extend_state_spaces(missing, degenerate_states=degenerate_states)
                if extended_spaces is not None:
                    new_perts = self.extend_VPT_representations(*extended_spaces)
                    # new_perts2 = self.get_VPT_representations()
                    # diffs = [x.asarray()-y.asarray() for x,y in zip(new_perts, new_perts2)]
                    # raise Exception(
                    #     [[np.min(d), np.max(d)] for d in diffs]
                    # )
                    self._zo_engs = None
                    self.representations = perturbations = self.PastIndexableTuple(new_perts)
                    flat_total_space = self.flat_total_space
                    N = self.total_space_dim
                else:
                    flat_total_space = self.flat_total_space
                    N = self.total_space_dim
                    perturbations = self.representations
                states = states.union(missing)
            else:
                flat_total_space = self.flat_total_space
                N = self.total_space_dim
                perturbations = self.representations
        else:
            flat_total_space = self.flat_total_space
            N = self.total_space_dim
            perturbations = self.representations

        return states, perturbations, flat_total_space, N

    def _get_corrections(self,
                         perturbations,
                         states,
                         order,
                         flat_total_space,
                         N,
                         checkpointer,
                         logger,
                         degenerate_states,
                         non_zero_cutoff=None,
                         handle_strong_couplings=None
                         ):
        """
        Just exists so we can do a recursion on this...

        :param perturbations:
        :type perturbations:
        :param states:
        :type states:
        :param order:
        :type order:
        :param flat_total_space:
        :type flat_total_space:
        :param N:
        :type N:
        :param checkpointer:
        :type checkpointer:
        :param logger:
        :type logger:
        :param degenerate_states:
        :type degenerate_states:
        :param non_zero_cutoff:
        :type non_zero_cutoff:
        :param handle_strong_couplings:
        :type handle_strong_couplings:
        :return:
        :rtype:
        """

        # print("????", states.excitations)

        with checkpointer:

            all_energies = np.zeros((len(states), order + 1))
            all_overlaps = np.zeros((len(states), order + 1))
            mode, limit, block_size = PerturbationTheoryCorrections.correction_storage_policy(
                self.correction_storage, self.correction_memory_limit,
                self.correction_block_size, self._correction_default_memory_limit,
                (len(states), order + 1, N))
            completed_row_cutoff = (
                Settings.non_zero_cutoff if non_zero_cutoff is None else non_zero_cutoff
            ) if self.use_block_recurrence and mode in {'sparse', 'compact'} else None
            self.correction_storage_report = dict(
                mode=mode, memory_limit=limit, block_size=block_size,
                memory_limit_role='dense/sparse selection threshold, not sparse retention cap',
                completed_row_cutoff=completed_row_cutoff,
                dense_bytes=len(states) * (order + 1) * N * np.dtype(float).itemsize)
            all_corrs = PerturbationTheoryCorrections.create_correction_array(
                mode, len(states), order, N, block_size,
                cutoff=completed_row_cutoff, block_bytes=self.correction_block_bytes,
                spill_directory=self.correction_spill_directory or self.state_space_scratch,
                resident_bytes=self.correction_resident_bytes,
                spill_bytes=self.correction_spill_bytes,
                spill_segment_bytes=self.correction_spill_segment_bytes)
            logger.log_print('Correction storage: {s}', s=self.correction_storage_report)
            all_energy_corrs = np.full((len(states), order + 1), None, dtype=object)

            with logger.block(tag="Applying Perturbation Theory"):
                logger.log_print(
                    [
                        "order: {o}",
                        "states: {n}",
                        "deg. spaces: {d}",
                        'deg. handling: {dm}',
                    ],
                    o=order,
                    n=len(states.indices),
                    d=0 if degenerate_states is None else len([1 for x in degenerate_states if len(x) > 1]),
                    dm=(
                        'Sakurai' if self.allow_sakurai_degs else
                        'mod. H' if self.drop_perturbation_degs else
                        'standard'
                    )
                )
                start = time.time()

                _ = []
                if degenerate_states is None:
                    degenerate_states = [[i] for i in states.indices]
                for deg_group in degenerate_states:
                    if not hasattr(deg_group, 'indices'):
                        deg_group = BasisStateSpace(flat_total_space.basis, deg_group, full_basis=self.full_basis)
                    deg_group = deg_group.as_sorted()
                    deg_group.deg_find_inds = None
                    _.append(deg_group)
                degenerate_states = _

                if self.drop_perturbation_degs:
                    dropped_els, perturbations = self.drop_deg_pert_els(perturbations, degenerate_states)

                    for deg_group in degenerate_states:
                        for n in deg_group.indices:
                            d2 = deg_group.take_states([n])
                            d2.deg_find_inds = None
                            energies, overlaps, corrs, ecorrs = self.apply_VPT_equations(n, deg_group,
                                                                                 None, None, None, None,
                                                                                 allow_PT_degs=False,
                                                                                 non_zero_cutoff=non_zero_cutoff,
                                                                                 perturbations=perturbations
                                                                                 )
                            res_index = states.find(n)
                            # need to explicitly update zero order energy in case we supplied a correction
                            if d2.deg_find_inds is not None:
                                energies[0] = self.zero_order_energies[d2.deg_find_inds[0]]
                            all_energies[res_index] = energies
                            all_energy_corrs[res_index] = ecorrs
                            all_corrs[res_index] = corrs
                            all_overlaps[res_index] = overlaps
                            del corrs
                else:
                    # loop over the degenerate sets
                    for deg_group in degenerate_states:
                        # we use this to build a pertubation operator that removes
                        # then entire set of degenerate states
                        deg_inds = flat_total_space.find(deg_group)
                        if len(deg_group) > 1:
                            if self.allow_sakurai_degs:
                                raise NotImplementedError("True degeneracy handling was purged")
                                deg_engs, zero_order_states, main_zero_states, subspaces = self._get_deg_eq_inputs(deg_inds)
                                main_subspace = (np.array([d[0] for d in deg_engs]), main_zero_states)
                            else:
                                deg_engs = zero_order_states = subspaces = main_subspace = [None]*len(deg_group)
                        else:
                            deg_engs = zero_order_states = subspaces = main_subspace = [None]

                        res_inds = states.find(deg_group.indices, missing_val=-1)
                        if self.use_block_recurrence and len(deg_group) > 1:
                            energies, overlaps, corrs, ecorrs = self.apply_VPT_block_equations(
                                deg_group, perturbations=perturbations,
                                non_zero_cutoff=non_zero_cutoff
                            )
                            present = res_inds > -1
                            all_energies[res_inds[present]] = energies[present]
                            all_energy_corrs[res_inds[present]] = ecorrs[present]
                            if mode != 'dense':
                                all_corrs.set_rows(res_inds, corrs)
                            else:
                                all_corrs[res_inds[present]] = corrs[present]
                            all_overlaps[res_inds[present]] = overlaps[present]
                            # Packing owns its rows; drop the dense block before
                            # evaluating the next solve on the RHS of an assignment.
                            del corrs
                            continue
                        for n, res_index, de, zo, s in zip(deg_group.indices, res_inds, deg_engs, zero_order_states, subspaces):
                            if res_index > -1:
                                energies, overlaps, corrs, ecorrs = self.apply_VPT_equations(
                                    n, deg_group, de, zo, main_subspace, s,
                                    allow_PT_degs=self.allow_sakurai_degs,
                                    non_zero_cutoff=non_zero_cutoff,
                                    perturbations=perturbations
                                )
                                all_energies[res_index] = energies
                                all_energy_corrs[res_index] = ecorrs
                                all_corrs[res_index] = corrs
                                all_overlaps[res_index] = overlaps
                                del corrs

                        # now we reorthogonalize degenerate states
                        if self.reorthogonalize_degenerate_states:
                            sparse_rows = mode != 'dense'
                            for k in range(1, order+1):
                                for i, (n, x) in enumerate(zip(res_inds, deg_inds)):
                                    for m, y in zip(res_inds[i + 1:], deg_inds[i+1:]):
                                        if sparse_rows:
                                            onm = -1 / 2 * sum(all_corrs.dot(n, j, m, k-j) for j in range(1, k))
                                            onn = -1 / 2 * sum(all_corrs.dot(n, j, n, k-j) for j in range(1, k))
                                            omm = -1 / 2 * sum(all_corrs.dot(m, j, m, k-j) for j in range(1, k))
                                            all_corrs.set_element(n, k, y, onm)
                                            all_corrs.set_element(m, k, x, onm)
                                            all_corrs.set_element(n, k, x, onn)
                                            all_corrs.set_element(m, k, y, omm)
                                        else:
                                            onm = -1 / 2 * sum(np.dot(all_corrs[n][j], all_corrs[m][k-j]) for j in range(1, k))
                                            onn = -1 / 2 * sum(np.dot(all_corrs[n][j], all_corrs[n][k-j]) for j in range(1, k))
                                            omm = -1 / 2 * sum(np.dot(all_corrs[m][j], all_corrs[m][k-j]) for j in range(1, k))
                                            all_corrs[n][k][y] = onm
                                            all_corrs[m][k][x] = onm
                                            all_corrs[n][k][x] = onn
                                            all_corrs[m][k][y] = omm

                end = time.time()
                logger.log_print(
                    "took {t}s",
                    t=round(end - start, 3)
                )

            # raise Exception(
            #     (np.sum(all_energies, axis=1) * UnitsData.convert("Hartrees", "Wavenumbers") - 4605.5)
            # )

            # now we recompute reduced state spaces for use in results processing
            # and we also convert the correction vectors to sparse representations
            tci = flat_total_space.indices
            N = len(tci)
            nstates = len(all_corrs)

            corr_mats, corr_inds = PerturbationTheoryCorrections.create_coupling_matrix(
                all_corrs,
                states,
                flat_total_space,
                nstates,
                order,
                non_zero_cutoff=non_zero_cutoff,
                # filters=self.state_space_filters
                filters=None,
                logger=logger,
                state_formatting=self.state_formatting
            )

            cs_states = SelectionRuleStateSpace(states, corr_inds, None)
            total_states = self.flat_total_space
            corrs = PerturbationTheoryCorrections.from_dicts(
                {
                    "states": states,
                    "coupled_states": cs_states,
                    "total_states": total_states,
                    "degenerate_states": degenerate_states
                },
                {
                    "energies": all_energies,
                    'energy_corrections':all_energy_corrs,
                    "wavefunctions": corr_mats,
                    "degenerate_transformation": None,
                    "degenerate_energies": None
                }
                , logger=self.logger
                , nondeg_hamiltonian_precision=self.nondeg_hamiltonian_precision
                , state_formatting=self.state_formatting
            )

            # if (
            #         degenerate_states is None
            #         or all(len(x) == 1 for x in degenerate_states)
            # ):
            sc, degenerate_correction_threshold = self.identify_strong_couplings(corrs)
            if len(sc) > 0:
                with self.logger.block(tag="Strongly coupled states (threshold={})".format(degenerate_correction_threshold)):
                    self.logger.log_print(
                        corrs.format_strong_couplings_report(sc, join=False)
                    )

                    if handle_strong_couplings and (sc is not None and len(sc) > 0):
                        for spec in self.degeneracy_handlers:
                            if hasattr(spec, 'wfc_threshold'):
                                construct = self.construct_strong_coupling_spaces
                                if (self.memory_constrained and mode == 'compact' and
                                        getattr(construct, '__func__', None) is
                                        PerturbationTheorySolver.construct_strong_coupling_spaces):
                                    degenerate_states = self._strong_coupling_groups(
                                        spec, sc, corrs, degenerate_correction_threshold)
                                    # The group-filter closure has returned. No amplitudes
                                    # from this discarded pass are needed by extension.
                                    corrs = corr_mats = corr_inds = cs_states = None
                                    meta = self._extend_strong_coupling_spaces(spec, degenerate_states, states)
                                else:
                                    # Preserve existing subclass/instance overrides.
                                    degenerate_states, meta = construct(
                                        spec, sc, corrs, states, degenerate_correction_threshold)
                                states, perturbations, flat_total_space, N = meta
                                if self.correction_storage == 'compact':
                                    # This branch returns the recomputed result. Its old
                                    # correction matrices and support views are obsolete;
                                    # release them before allocating the next pass.
                                    corrs = corr_mats = corr_inds = cs_states = None
                                if spec.iterations > 1:
                                    spec.iterations -= 1
                                    degenerate_states = None
                                    corrs = self._get_corrections(
                                        perturbations,
                                        states,
                                        order,
                                        flat_total_space,
                                        N,
                                        checkpointer,
                                        logger,
                                        degenerate_states,
                                        handle_strong_couplings=True,
                                        non_zero_cutoff=non_zero_cutoff
                                    )
                                else:
                                    with self.logger.block(tag="Redoing PT with strong couplings handled"):
                                        with self.logger.block(tag="Degenerate groups:"):
                                            for g in degenerate_states.flat:
                                                if len(g) > 1:
                                                    self.logger.log_print(g, message_prepper=lambda a: StateMaker.format_states(a, self.state_formatting))
                                    if self.degenerate_states is not None:
                                        degenerate_states = self.merge_deg_spaces(
                                            [degenerate_states, self.degenerate_states]
                                        )
                                    corrs = self._get_corrections(
                                        perturbations,
                                        states,
                                        order,
                                        flat_total_space,
                                        N,
                                        checkpointer,
                                        logger,
                                        degenerate_states,
                                        handle_strong_couplings=False,
                                        non_zero_cutoff=non_zero_cutoff
                                    )
                                return corrs
            else:
                sc = None


            # raise Exception(sc)

            if self.results is None or isinstance(self.results, NullCheckpointer):
                try:
                    checkpointer['corrections'] = {
                        "states": states.excitations,
                        "total_states": self._archive_state_vectors(total_states),
                        'energies': all_energies,
                        'wavefunctions': corr_mats
                    }
                except KeyError:
                    pass
            else:
                with self.results:
                    try:
                        self.results['corrections'] = {
                            "states": states.excitations,
                            "total_states": self._archive_state_vectors(total_states),
                            'energies': all_energies,
                            'wavefunctions': corr_mats
                        }
                    except KeyError:
                        pass

            # with self.logger.block(tag="overlap matrix"):
            #     self.logger.log_print(str(np.sum(corrs.get_overlap_matrices(), axis=0)).splitlines())

        return self.PTResults(corrs, degenerate_states)

    def _archive_state_vectors(self, space):
        # Null output sinks must not force materialization of the global basis.
        if isinstance(self.results, NullCheckpointer):
            return None
        if not hasattr(space, 'store'):
            return space.excitations
        dtype = BasisStateSpace.excitations_dtype
        result = space.store.array('archived_states',dtype,(len(space),space.ndim))
        offset = 0
        for block in space.iter_excitations():
            result[offset:offset+len(block)] = block
            offset += len(block)
        return result

    PTResults = namedtuple("PTResults", ["corrections", "degeneracies"])

    @staticmethod
    def _safe_dot(a, b):
        """
        **LLM Docstring**

        Generalized dot product that treats a bare `0` operand (rather than a zero-valued array) as an exact zero result, short-circuiting the multiplication -- used throughout the perturbation-theory recursion, where many terms really are the literal integer `0` rather than zero arrays.

        :param a: the left operand, or `0`
        :type a: np.ndarray | SparseArray | int | float
        :param b: the right operand, or `0`
        :type b: np.ndarray | SparseArray | int | float
        :return: `0` if either operand is a numeric zero, otherwise the dot product of `a` and `b`
        :rtype: np.ndarray | SparseArray | int
        """
        # generalizes the dot product so that we can use 0 as a special value...
        if (
                isinstance(a, (int, np.integer, float, np.floating)) and a == 0
                or isinstance(b, (int, np.integer, float, np.floating)) and b == 0
        ):
            return 0

        if isinstance(a, np.ndarray):
            doots = np.dot(a, b)
        else:
            doots = a.dot(b)

        if isinstance(b, np.ndarray) and isinstance(doots, SparseArray):
            doots = doots.asarray()

        return doots
    def apply_VPT_equations(self,
                            state_index,
                            degenerate_space_indices,
                            degenerate_energies,
                            zero_order_state,
                            degenerate_subspace,
                            degenerate_subsubspace,
                            perturbations=None,
                            allow_PT_degs=None,
                            ignore_odd_orders=None,
                            intermediate_normalization=None,
                            non_zero_cutoff=None
                            ):
        """
        Applies VPT equations, dispatching based on how many
        degeneracies we need to handle

        :param state_index: the index of the primary state being treated using the PT
        :type state_index: int
        :param degenerate_space_indices: the indices corresponding to degeneracies with the primary state in the zero-order picture
        :type degenerate_space_indices: np.ndarray[int]
        :param degenerate_energies: the first and (possibly) second order correction to the energies
        :type degenerate_energies: Iterable[float | None]
        :param zero_order_states: the vector for the proper zero-order state corresponding ot state_index
        :type zero_order_states: np.ndarray[float]
        :param degenerate_subsubspace: the set of vectors for the zero-order states in the secondary degenerate subspace
        :type degenerate_subsubspace: tuple[np.ndarray[float], np.ndarray[int]]
        :param non_zero_cutoff: cutoff for when a term can be called zero for performance reasons
        :type non_zero_cutoff: float
        :return:
        :rtype:
        """
        if non_zero_cutoff is None:
            non_zero_cutoff = Settings.non_zero_cutoff

        if ignore_odd_orders is None:
            ignore_odd_orders=self.ignore_odd_orders
        if allow_PT_degs is None:
            allow_PT_degs = self.allow_sakurai_degs
        if intermediate_normalization is None:
            intermediate_normalization = self.intermediate_normalization
        return self.apply_VPT_nondeg_equations(state_index, degenerate_space_indices, non_zero_cutoff=non_zero_cutoff,
                                               ignore_odd_orders=ignore_odd_orders,
                                               intermediate_normalization=intermediate_normalization,
                                               perturbations=perturbations
                                               )

    def _apply_VPT_block_equations_in_place(self, group, perturbations=None,
                                           non_zero_cutoff=None):
        """Evaluate the same block recurrence using C[k] as its RHS storage."""
        if non_zero_cutoff is None:
            non_zero_cutoff = Settings.non_zero_cutoff
        H = self.representations if perturbations is None else perturbations
        inds = self.flat_total_space.find(group)
        nstates, nbasis = len(inds), len(self.flat_total_space)
        order = self.order
        e0 = self.zero_order_energies[inds]
        sparse_low = self.sparse_low_order_corrections
        if sparse_low:
            # Existing CSR selectors contain only one entry per reference.
            selector = sp.csr_matrix((np.ones(nstates), inds, np.arange(nstates+1)),
                                     shape=(nstates, nbasis))
            C = [selector, None]+[np.zeros((nstates, nbasis)) for _ in range(2, order+1)]
        else:
            C = np.zeros((order+1, nstates, nbasis), dtype=float)
            C[0, np.arange(nstates), inds] = 1
        F = np.zeros((order+1, nstates, nstates), dtype=float)
        F[0] = np.diag(e0)
        energies = np.zeros((nstates, order+1), dtype=float)
        energies[:, 0] = e0
        overlaps = np.zeros_like(energies)
        overlaps[:, 0] = 1
        energy_terms = np.full((nstates, order+1), None, dtype=object)
        budget = self.wavefunction_workspace_bytes
        p_order = np.argsort(inds)
        p_sorted = inds[p_order]
        # One denominator buffer is shared by all orders and reference states.
        chunk_size = max(1, min(nbasis, budget//(4*np.dtype(float).itemsize)))
        denominator = np.empty(chunk_size, dtype=float)
        basis_energies = self.zero_order_energies

        def project_product(rows, columns, product):
            lo = np.searchsorted(p_sorted, rows.start)
            hi = np.searchsorted(p_sorted, rows.stop)
            if hi > lo:
                projection[p_order[lo:hi], columns] = product[p_sorted[lo:hi]-rows.start]

        for k in range(1, order+1):
            if self.intermediate_normalization or k == 1:
                normalizer = np.zeros((nstates, nstates))
            else:
                # C is real by construction; conjugating it would copy a full
                # B-by-N correction slice for no numerical change.
                products = (C[i] @ C[k-i].T for i in range(1, k))
                normalizer = -.5*sum(x.toarray() if sp.issparse(x) else x for x in products)
                if not self.reorthogonalize_degenerate_states:
                    normalizer = np.diag(np.diag(normalizer))
            F[k] = e0[:, None]*normalizer-normalizer*e0[None, :]
            if sparse_low and k == 1:
                h1 = H[1].data if isinstance(H[1], SparseArray) else H[1]
                if sp.issparse(h1):
                    # Select before transpose: never copy the full Hamiltonian.
                    rhs = -h1[:, inds].T.tocsr()
                elif isinstance(h1, (int, float, np.integer, np.floating)) and h1 == 0:
                    rhs = sp.csr_matrix((nstates, nbasis))
                else:
                    rhs = sp.csr_matrix(-np.asarray(h1)[:, inds].T)
                if not np.can_cast(rhs.dtype, float, casting='same_kind'):
                    raise TypeError('Block recurrence requires real Hamiltonian coefficients')
                # The dense recurrence accumulates in float64 even when an
                # input operator uses narrower floating-point/integer storage.
                rhs = rhs.astype(float, copy=False)
                rhs.sum_duplicates()
                rhs.sort_indices()
                C[k] = rhs
                projection = -rhs[:, inds].T.toarray()
                F[k] += projection
                terms = [np.diag(projection)]
            else:
                rhs = C[k]  # Lower-order C[t] is never overwritten.
                terms = []
            for s in (() if sparse_low and k == 1 else range(k, 0, -1)):
                t = k-s
                projection = np.zeros((nstates, nstates))
                if not isinstance(H[s], (int, float, np.integer, np.floating)):
                    nput.matrix_dot_into(H[s], C[t].T, rhs.T, alpha=-1, beta=1,
                        working_memory=budget, block_callback=project_product)
                elif H[s] != 0:
                    raise TypeError('Expected a matrix or zero Hamiltonian term')
                term = projection
                if t > 0:
                    nput.matrix_dot_into(F[s].T, C[t], rhs, alpha=1, beta=1,
                        working_memory=budget)
                    term = term-C[t][:, inds].T @ F[s]
                F[k] += term
                terms.append(np.diag(term))
            energies[:, k] = np.diag(F[k])
            if self.ignore_odd_orders and k % 2:
                energies[:, k] = 0
            for a, energy in enumerate(e0):
                for start in range(0, nbasis, chunk_size):
                    stop = min(nbasis, start+chunk_size)
                    gap = denominator[:stop-start]
                    np.subtract(basis_energies[start:stop], energy, out=gap)
                    lo = np.searchsorted(p_sorted, start)
                    hi = np.searchsorted(p_sorted, stop)
                    local_p = p_sorted[lo:hi]-start
                    gap[local_p] = 1
                    if k == 1 and np.any(np.abs(gap) < non_zero_cutoff):
                        # Keep the established full diagnostic and exception.
                        # This error-only path never builds a resolvent CSR.
                        self._get_Pi0(inds, E0=energy, non_zero_cutoff=non_zero_cutoff)
                        raise ValueError('Unexpected unprojected near-zero denominator')
                    np.reciprocal(gap, out=gap)
                    gap[local_p] = 0
                    if sp.issparse(rhs):
                        lo_row, hi_row = rhs.indptr[a:a+2]
                        columns = rhs.indices[lo_row:hi_row]
                        lo_col, hi_col = np.searchsorted(columns, [start, stop])
                        values = rhs.data[lo_row+lo_col:lo_row+hi_col]
                        np.multiply(values, gap[columns[lo_col:hi_col]-start], out=values)
                    else:
                        np.multiply(rhs[a, start:stop], gap, out=rhs[a, start:stop])
                energy_terms[a, k] = (None if self.ignore_odd_orders and k % 2 else
                                      np.asarray([term[a] for term in terms]))
            if sp.issparse(C[k]):
                # C1 has zero P projection and zero first-order normalizer.
                C[k].eliminate_zeros()
            else:
                C[k][:, inds] = normalizer.T
            overlaps[:, k] = np.diag(normalizer)
        # Only the explicit sparse-low-order option changes the block container.
        corrections = tuple(C) if sparse_low else C.transpose(1, 0, 2)
        return energies, overlaps, corrections, energy_terms

    def apply_VPT_block_equations(self, group, perturbations=None, non_zero_cutoff=None):
        """Apply the same recurrence to all columns of one projected block.

        Corrections are stored as (order, reference state, basis state), so
        right multiplication of ket columns by F is F.T @ C in this layout.
        With sparse_low_order_corrections=True, corrections are a tuple of
        (reference, basis) matrices by order: CSR selectors/C1 and dense C2+.
        The default retains the existing (reference, order, basis) ndarray.
        F is constructed from HY = YF; higher-order normalization is imposed
        before advancing the recurrence. No denominator is inverted inside P.
        """
        if self.memory_constrained and self.order > 0:
            return self._apply_VPT_block_equations_in_place(
                group, perturbations=perturbations, non_zero_cutoff=non_zero_cutoff)
        if non_zero_cutoff is None:
            non_zero_cutoff = Settings.non_zero_cutoff
        H = self.representations if perturbations is None else perturbations
        inds = self.flat_total_space.find(group)
        nstates, nbasis = len(inds), len(self.flat_total_space)
        order = self.order
        e0 = self.zero_order_energies[inds]
        C = np.zeros((order + 1, nstates, nbasis), dtype=float)
        C[0, np.arange(nstates), inds] = 1
        F = np.zeros((order + 1, nstates, nstates), dtype=float)
        F[0] = np.diag(e0)
        energies = np.zeros((nstates, order + 1), dtype=float)
        energies[:, 0] = e0
        overlaps = np.zeros_like(energies)
        overlaps[:, 0] = 1
        energy_terms = np.full((nstates, order + 1), None, dtype=object)
        # Reuse the existing denominator checks and sign convention for Pi0.
        resolvents = [self._get_Pi0(inds, E0=e, non_zero_cutoff=non_zero_cutoff) for e in e0]
        dot = self._safe_dot
        for k in range(1, order + 1):
            if self.intermediate_normalization or k == 1:
                normalizer = np.zeros((nstates, nstates))
            else:
                normalizer = -.5 * sum(C[i].conj() @ C[k-i].T for i in range(1, k))
                if not self.reorthogonalize_degenerate_states:
                    normalizer = np.diag(np.diag(normalizer))
            # PH0 Yk - P Yk F0, including unequal unperturbed energies in P.
            F[k] = e0[:, None] * normalizer - normalizer * e0[None, :]
            rhs = np.zeros((nstates, nbasis))
            terms = []
            for s in range(k, 0, -1):
                t = k - s
                hc = dot(H[s], C[t].T)
                if isinstance(hc, (int, float, np.integer, np.floating)):
                    hc = np.zeros((nstates, nbasis))
                else:
                    hc = np.asarray(hc).T
                rhs -= hc
                term = hc[:, inds].T
                if t > 0:
                    rhs += F[s].T @ C[t]
                    term = term - C[t][:, inds].T @ F[s]
                F[k] += term
                terms.append(np.diag(term))
            # F[k] is a matrix even when its odd-order diagonal vanishes.
            energies[:, k] = np.diag(F[k])
            if self.ignore_odd_orders and k % 2:
                energies[:, k] = 0
            for a in range(nstates):
                C[k, a] = dot(resolvents[a], rhs[a])
                energy_terms[a, k] = (None if self.ignore_odd_orders and k % 2 else
                                      np.asarray([term[a] for term in terms]))
            C[k][:, inds] = normalizer.T
            overlaps[:, k] = np.diag(normalizer)
        return energies, overlaps, C.transpose(1, 0, 2), energy_terms

    def apply_VPT_nondeg_equations(self,
                                   state_index,
                                   deg_group,
                                   perturbations=None,
                                   non_zero_cutoff=None,
                                   check_overlap=None,
                                   intermediate_normalization=False,
                                   ignore_odd_orders=False
                                   ):
        """
        Does the dirty work of doing the VPT iterative equations.

        :return:
        :rtype:
        """

        if non_zero_cutoff is None:
            non_zero_cutoff = Settings.non_zero_cutoff

        if intermediate_normalization:
            check_overlap=False
        elif check_overlap is None:
            check_overlap = self.check_overlap

        n = state_index
        e_vec_full = self.zero_order_energies

        order = self.order
        total_state_space = self.flat_total_space

        energies = np.zeros((order + 1,), dtype=float)
        overlaps = np.zeros((order + 1,), dtype=float)
        energy_corrs = [None] * (order + 1)
        corrs = np.zeros((order + 1, len(total_state_space)), dtype=float)  # can I make this less expensive in general?

        # find the state index in the coupled subspace
        n_ind = total_state_space.find(n)

        logger = NullLogger() if self.logger is None else self.logger

        block_tag_formatter = lambda:'getting corrections for state {}/{}'.format(
            n,
            StateMaker.format_state(total_state_space.take_subspace(n_ind).excitations[0], self.state_formatting)
        )

        with logger.block(
            tag=block_tag_formatter,
            log_level=logger.LogLevel.Debug
        ):
            D = deg_group
            deg_inds = (n_ind,)
            if D is not None:
                if D.deg_find_inds is None:
                    D.deg_find_inds = total_state_space.find(D)
                deg_inds = D.deg_find_inds
                if len(D) > 1:
                    logger.log_print('Degenerate space: {D}',
                                          D=D,
                                          preformatter=logger.preformat_keys({"D":lambda D:D.indices}),
                                          log_level=logger.LogLevel.Debug
                                          )
            E0 = e_vec_full[n_ind]
            pi = self._get_Pi0(deg_inds, E0=E0, non_zero_cutoff=non_zero_cutoff)

            energies[0] = E0
            logger.log_print('Zero-order energy: {e}',
                                  e=E0 * UnitsData.convert("Hartrees", "Wavenumbers"),
                                  log_level=logger.LogLevel.Debug
                                  )
            # self.logger.log_print("{n}: E0={E}", n=n_ind, E=E0)
            overlaps[0] = 1
            corrs[0, n_ind] = 1
            H = self.representations if perturbations is None else perturbations


            dot = self._safe_dot
            takeDiag = lambda h, n_ind: h[n_ind, n_ind] if not isinstance(h, (int, np.integer, float, np.floating)) else 0.
            take = lambda h, el: h[el] if not isinstance(h, (int, np.integer, float, np.floating)) else 0.
            for k in range(1, order + 1):  # to actually go up to target order
                #         En^(k) = <n^(0)|H^(k)|n^(0)> + sum(<n^(0)|H^(k-i)|n^(i)> - E^(k-i)<n^(0)|n^(i)>, i=1...k-1)
                if ignore_odd_orders and k % 2 == 1:
                    logger.log_print('Skipping order {k} for the energy (assumed to be 0)', k=k, log_level=logger.LogLevel.Debug)
                    energy_terms = None
                    Ek = 0
                # elif ignore_odd_orders: # Tried to get the 2n + 1 trick working but...it doesn't work?
                #     Ek = (
                #             takeDiag(H[k], n_ind)
                #             + sum(dot(take(H[k - i], n_ind), corrs[i]) - energies[k - i] * overlaps[i]
                #                   for i in range(1, (k + 1) // 2))
                #     )
                else:
                    energy_terms = [takeDiag(H[k], n_ind)] + [
                            dot(take(H[k - i], n_ind), corrs[i]) - energies[k - i] * overlaps[i]
                            for i in range(1, k)
                        ]
                    energy_terms = np.array([
                        x.flatten()[0] if not isinstance(x, (int, float, np.integer, np.floating)) else x
                        for x in energy_terms
                    ])
                    logger.log_print(
                        energy_terms,
                        message_prepper=lambda energy_terms:['Energy terms at order {k} in cm^-1:'] + [
                            '{} = {}'.format(s, e) for s, e in
                            zip(
                                ["<n(0)|H({})|n(0)>".format(k)] + [
                                    "<n(0)|H({0})-E({0})|n({1})>".format(
                                        k - i, i
                                    ) for i in range(1, k)
                                ],
                                energy_terms * UnitsData.convert("Hartrees", "Wavenumbers")
                            )
                        ],
                        k=k,
                        log_level=logger.LogLevel.Debug
                    )
                    Ek = np.sum(energy_terms)
                energy_corrs[k] = energy_terms
                energies[k] = Ek
                #   <n^(0)|n^(k)> = -1/2 sum(<n^(i)|n^(k-i)>, i=1...k-1)
                #         |n^(k)> = sum(Pi_n (En^(k-i) - H^(k-i)) |n^(i)>, i=0...k-1) + <n^(0)|n^(k)> |n^(0)>
                corrs[k] = sum(
                    dot(pi, energies[k - i] * corrs[i] - dot(H[k - i], corrs[i]))
                        if abs(energies[k - i]) > non_zero_cutoff else
                     -dot(pi, dot(H[k - i], corrs[i])) # just cut out a potentially unnecessary dense cast
                    for i in range(0, k)
                )

                if check_overlap:
                    should_be_zero = corrs[k][deg_inds]
                    if (np.abs(should_be_zero) > 0).any():
                        raise ValueError("Perturbation operator should have made overlap of state {} with {} zero...got {} instead".format(
                            n, D, should_be_zero
                        ))

                if intermediate_normalization:
                    ok = 0.0
                else:
                    ok = -1 / 2 * sum(dot(corrs[i], corrs[k - i]) for i in range(1, k))
                    logger.log_print([
                        'Overlap at order {k}:'
                        '<n(0)|n({k})> = {ok}'
                        ], k=k, ok=ok, log_level=logger.LogLevel.Debug)
                overlaps[k] = ok
                corrs[k][n_ind] = ok  # pi (the perturbation operator) ensures it's zero before this

            if check_overlap:
                # full_wfn = np.sum(corrs, axis=0)
                # ov = np.dot(full_wfn, full_wfn)
                ov_parts = [ [dot(corrs[k-i], corrs[i]) for i in range(k+1)] for k in range(order+1)]
                ov = sum(np.sum(v) for v in ov_parts)
                if abs(ov - 1) > .005:
                    raise ValueError(
                        "state {} isn't normalized (overlap = {}, bits {})".format(
                            state_index, ov, ov_parts
                        ))

        return energies, overlaps, corrs, energy_corrs

    def apply_VPT_2k1_rules(self,
                           existing_corrs,
                           perturbations=None
                           ):
        """
        Apply expressions allowing for obtaining higher-order
        corrections to the energies from lower-order corrections to the
        wavefunctions

        :param existing_corrs:
        :type existing_corrs:
        :param perturbations:
        :type perturbations:
        :return:
        :rtype:
        """
        ...

    #endregion

    #region Handle Post-PT Variational Stuff
    @staticmethod
    def _fmt_depert_engs(total_state_space, base_energies, state_formatting='quanta'):
        """
        **LLM Docstring**

        Build a lazily-evaluated formatter for logging a set of deperturbed state energies: subtracting off the ground-state energy to display frequencies (rather than raw energies) whenever the ground state is present in `total_state_space`.

        :param total_state_space: the state space the energies correspond to
        :type total_state_space: BasisStateSpace
        :param base_energies: the (deperturbed) energies to format
        :type base_energies: np.ndarray
        :return: a zero-argument formatter function producing the log lines, for use with the logger's lazy message-formatting
        :rtype: callable
        """

        def fmt(*a, base_energies=base_energies):
            """
            **LLM Docstring**

            Format each state's padded occupation label (or legacy vector) beside its energy/frequency.

            :param a: unused positional arguments, accepted for logger message-formatter compatibility
            :type a: tuple
            :param base_energies: the energies to format, captured from the enclosing scope
            :type base_energies: np.ndarray
            :return: the formatted log lines, with a descriptive label first
            :rtype: list[str]
            """
            states = total_state_space.excitations
            state_sums = np.sum(total_state_space.excitations, axis=1)
            gs_pos = np.where(state_sums == 0)
            label = "Deperturbed States/Energies:"
            if len(gs_pos) > 0 and len(gs_pos[0]) > 0:
                label = "Deperturbed States/Frequencies:"
                gs_pos = gs_pos[0]
                if len(gs_pos) > 0:
                    gs_pos = gs_pos[0]
                gs_eng = base_energies[gs_pos]
                base_energies = base_energies - gs_eng
                base_energies[gs_pos] = gs_eng

            if state_formatting == 'vector':
                labels = [('{:.0f} ' * len(s)).format(*s) for s in states]
            else:
                labels = [StateMaker.parse_state(s) for s in states]
                width = max([len('State')] + [len(s) for s in labels]) + 1
                labels = [s.ljust(width) for s in labels]

            return [label] + [
                s + " {:12.4f}".format(e)
                for s, e in zip(
                    labels,
                    UnitsData.convert("Hartrees", "Wavenumbers") * base_energies
                )
            ]

        return fmt
    def apply_post_PT_variational_calc(self, degenerate_states, corrs):
        """
        Applies degenerate perturbation theory by building a representation
        for the degenerate terms in the Hamiltonian.
        This is then diagonalized, allowing the degenerate states to be expressed
        in the basis of non-degenerate states

        :param H:
        :type H: Iterable[SparseArray]
        :param corrs: the standard PerturbationTheory Corrections object that comes out of the application of non-deg PT
        :type corrs: PerturbationTheoryCorrections
        :param degenerate_states: population of degenerate states
        :type degenerate_states:
        :param logger:
        :type logger: Logger
        :return:
        :rtype:
        """

        # we pull the states and total states from the corrections object
        total_state_space = corrs.states  # type: BasisStateSpace

        # set up space to store the degenerate energies and rotations coming from the
        # sets of diagonalizations
        energies = np.zeros(len(total_state_space))
        base_energies = corrs.energies  # for when we're not rotating

        logger = NullLogger() if self.logger is None else self.logger
        logger.log_print(
            None,
            message_prepper=self._fmt_depert_engs(total_state_space, base_energies, self.state_formatting)
        )

        rotations = []
        rotation_vals = []
        rotation_row_inds = []
        rotation_col_inds = []

        ndeg_ham_corrs = []
        for group_num,group in enumerate(degenerate_states):
            deg_inds, H_nd, deg_rot, deg_engs = corrs.get_degenerate_transformation(
                group,
                self.representations,
                label="Block {group_num}".format(group_num=group_num),
                gaussian_resonance_handling=self.gaussian_resonance_handling,
                local_coupling_hamiltonian=self.local_coupling_hamiltonian,
                local_coupling_order=self.local_coupling_order
            )
            if H_nd is not None:
                ndeg_ham_corrs.append(H_nd)
                rotations.append(deg_rot)
                energies[deg_inds] = deg_engs
                rotation_vals.append(deg_rot.flatten())
                deg_rows, deg_cols = np.array([p for p in itertools.product(deg_inds, deg_inds)]).T
                rotation_row_inds.append(deg_rows)
                rotation_col_inds.append(deg_cols)
            else:
                for i in deg_inds:
                    # we'll be a little bit inefficient for now and speed up later
                    energies[i] = base_energies[i]
                    rotation_vals.append([1.])
                    rotation_row_inds.append([i])
                    rotation_col_inds.append([i])

        if self.results is None or isinstance(self.results, NullCheckpointer):
            try:
                self.checkpointer["degenerate_data"] = {
                    "states": [d.excitations for d in degenerate_states],
                    "energies": energies,
                    "hamiltonians": ndeg_ham_corrs,
                    "rotations": rotations
                }
            except KeyError:
                pass
        else:
            with self.results:
                try:
                    self.results["degenerate_data"] = {
                        "states": [d.excitations for d in degenerate_states],
                        "energies": energies,
                        "hamiltonians": ndeg_ham_corrs,
                        "rotations": rotations
                    }
                except KeyError:
                    pass

        rotation_vals = np.concatenate(rotation_vals)
        rotation_row_inds = np.concatenate(rotation_row_inds)
        rotation_col_inds = np.concatenate(rotation_col_inds)

        rotations = SparseArray.from_data(
            (
                rotation_vals,
                (
                    rotation_row_inds,
                    rotation_col_inds
                )
            ),
            shape=(len(energies), len(energies))
        )

        return energies, rotations, ndeg_ham_corrs
    def drop_deg_pert_els(self, perts, deg_groups):
        """

        :param perts:
        :type perts:
        :param deg_groups:
        :type deg_groups:
        :return:
        :rtype:
        """
        deg_grop_inds = []
        for g in deg_groups:
            if g.deg_find_inds is None:
                g.deg_find_inds = self.flat_total_space.find(g)
            deg_grop_inds.append(g.deg_find_inds)
        pert_blocks = []
        perts = self.PastIndexableTuple([perts[0]] + [p.copy() for p in perts[1:]])

        block_logger = NullLogger() if self.logger is None else self.logger

        with block_logger.block(tag='modifying perturbations', log_level=block_logger.LogLevel.Debug):
            for d,g in zip(deg_grop_inds, deg_groups):
                if len(d) > 1:
                    block_logger.log_print(
                        None,
                        message_prepper=lambda *a: ["dropping elements coupling degenerate space:"] + StateMaker.format_states(g, self.state_formatting),
                        log_level=block_logger.LogLevel.Debug
                    )
                    idx = tuple(np.array([x for x in itertools.product(d, d) if x[0] != x[1]]).T)
                    els = []
                    for p in perts[1:]:
                        els.append(p[idx].flatten())
                        p[idx] = 0.
                    pert_blocks.append([idx, els])

                    triu = np.where(idx[0] > idx[1])
                    def pad_els(el, triu=triu):
                        """
                        **LLM Docstring**

                        Reconstruct the full symmetric matrix of zeroed-out coupling-element magnitudes (in wavenumbers, rounded to integers) for one degenerate group from its flattened upper-triangular values, for logging.

                        :param el: the flattened coupling-element values (for every off-diagonal pair in the group) to reassemble
                        :type el: np.ndarray
                        :param triu: the upper-triangular index mask (relative to `el`'s ordering) identifying which flattened entries belong in the upper triangle, captured from the enclosing scope
                        :type triu: tuple[np.ndarray, np.ndarray]
                        :return: the reconstructed, symmetric, wavenumber-converted coupling matrix for the group
                        :rtype: np.ndarray
                        """
                        e = np.zeros((len(d), len(d)))
                        e[np.triu_indices_from(e, k=1)] = el[triu]
                        e = np.round(e * UnitsData.convert("Hartrees", "Wavenumbers")).astype(int)
                        e[np.tril_indices_from(e, k=-1)] = e[np.triu_indices_from(e, k=1)]
                        return e
                    block_logger.log_print(
                        None,
                        message_prepper = lambda *a: ["zeroed out coupling elements:"] + sum(
                            (str(pad_els(e)).splitlines() for e in els),
                            []
                        ),
                        log_level=block_logger.LogLevel.Debug
                    )

        return pert_blocks, perts

    #endregion
