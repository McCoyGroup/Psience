
import numpy as np, itertools, time, types
from collections import namedtuple

from McUtils.Numputils import SparseArray
from McUtils.Scaffolding import Logger, NullLogger, NullCheckpointer
from McUtils.Parallelizers import Parallelizer, SerialNonParallelizer
from McUtils.Data import UnitsData
from McUtils.Combinatorics import LatticePathGenerator

from ..BasisReps import Representation, BasisStateSpace, BasisMultiStateSpace, SelectionRuleStateSpace, BraKetSpace

from .StateFilters import PerturbationTheoryStateSpaceFilter
from .DegeneracySpecs import DegenerateMultiStateSpace, DegeneracySpec
from .Common import *
from .Corrections import *

__reload_hook__ = [ "..BasisReps", ".DegeneracySpecs", ".Corrections", ".StateFilters", ".Common" ]

__all__ = [
    "PerturbationTheorySolver"
]

class PerturbationTheorySolver:
    """
    A solver that applies perturbation theory
    given a series of corrections and population of states.
    Supports degenerate and non-degenerate PT.
    """

    def __init__(self,
                 perturbations, states,
                 coupled_states=None,
                 order=2,
                 total_space=None,
                 flat_total_space=None,
                 state_space_iterations=None,
                 state_space_terms=None,
                 state_space_filters=None,
                 target_property_rules=None,
                 allow_sakurai_degs=False,
                 allow_post_PT_calc=True,
                 modify_degenerate_perturbations=False,
                 gaussian_resonance_handling=False,
                 ignore_odd_order_energies=False,
                 intermediate_normalization=False,
                 zero_element_warning=True,
                 degenerate_states=None,
                 handle_strong_couplings=False,
                 # strong_coupling_test_modes=None,
                 # strong_couplings_state_filter=None,
                 # strongly_coupled_group_filter=None,
                 # strong_coupling_zero_order_energy_cutoff=4.5e-3,
                 low_frequency_mode_cutoff=1.15e-3,
                 # extend_strong_coupling_spaces=True,
                 zero_order_energy_corrections=None,
                 memory_constrained=False,
                 keep_hamiltonians=None,
                 logger=None,
                 parallelizer=None,
                 checkpointer=None,
                 results=None,
                 checkpoint_keys=None,
                 use_cached_representations=True,
                 use_cached_basis=True
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
        :param logger:
        :type logger:
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
        self._reps = None
        self.order = order
        self.state_space_iterations = state_space_iterations
        self.state_space_terms = state_space_terms

        if isinstance(state_space_filters, (types.FunctionType, types.MethodType, types.LambdaType)):
            self.state_space_filter_generator = state_space_filters
            self.state_space_filters = PerturbationTheoryStateSpaceFilter.from_data(states, self.state_space_filter_generator(states))
        else:
            self.state_space_filter_generator = None
            self.state_space_filters = PerturbationTheoryStateSpaceFilter.from_data(states, state_space_filters)
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

        self.degeneracy_spec = degenerate_states
        self._deg_states = None
        if self.degeneracy_spec is None:
            if handle_strong_couplings is None:
                handle_strong_couplings = True
            if handle_strong_couplings is not False:
                self.degeneracy_spec = DegeneracySpec.from_spec({
                    'wfc_threshold':'auto' if handle_strong_couplings is True else handle_strong_couplings
                })
        self.handle_strong_couplings = hasattr(self.degeneracy_spec, 'wfc_threshold')

        # self.handle_strong_couplings = handle_strong_couplings
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
        self.gaussian_resonance_handling = gaussian_resonance_handling
        self.zero_element_warning = zero_element_warning

        self.memory_constrained=memory_constrained
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

    class PastIndexableTuple(tuple):
        def __getitem__(self, item):
            if isinstance(item, (int, np.integer)) and item >= len(self):
                return 0
            else:
                return super().__getitem__(item)

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
        self._reps = reps

    @property
    def degenerate_spaces(self):
        """
        :return:
        :rtype:
        """
        if self._deg_states is None:
            spec = self.degeneracy_spec
            if isinstance(spec, DegeneracySpec) and spec.application_order != 'pre':
                spec = None
            self._deg_states = DegenerateMultiStateSpace.from_spec(
                spec,
                solver=self,
                full_basis=self.full_basis
            )
        return self._deg_states

    @property
    def zero_order_energies(self):
        """
        :return:
        :rtype:
        """
        if self._zo_engs is None:
            H0 = self.representations[0]
            e_vec_full = np.diag(H0) if isinstance(H0, np.ndarray) else H0.diag
            if isinstance(e_vec_full, SparseArray):
                e_vec_full = e_vec_full.asarray()
            self._zo_engs = e_vec_full
            if self.zero_order_energy_corrections is not None:
                for k,v in self.zero_order_energy_corrections:
                    if not isinstance(k, int):
                        k = self.flat_total_space.find([k]) # slow but w/e
                    if k >= len(self._zo_engs):
                        self.logger.log_print(
                            "WARNING: zero-order correction {} for state {} not used as state is not in basis".format(
                                v * UnitsData.convert("Hartrees", "Wavenumbers"),
                                k
                            )
                        )
                    self._zo_engs[k] = v
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
                                                                           diagonal=True
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
                                                                           memory_constrained=self.memory_constrained
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

            return H
    def extend_VPT_representations(self, new_flat_space, new_states):

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
                raise ValueError(
                    "degeneracies encountered: states {} and {} other states are degenerate (average energy: {} stddev: {})".format(
                        degenerate_subspace,
                        len(zero_checks),
                        np.average(bad_vec),
                        np.std(bad_vec)
                    ))
            else:
                raise ValueError(
                    "degeneracies encountered: states {} and {} are degenerate (average energy: {} stddev: {})".format(
                        degenerate_subspace,
                        zero_checks,
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
                    if not isinstance(cs, (BasisStateSpace, SelectionRuleStateSpace))
                    else cs
                    for cs in self._coupled_states
                ]

            with logger.block(tag="precomputing coupled space indices"):
                # otherwise they get calculated twice
                start = time.time()
                for s in self._coupled_states:
                    if s is not None:
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
        if self._total_space is None:
            with self.logger.block(tag="generating total space"):

                start = time.time()

                space_list = [self.states] + [s for s in self._coupled_states if s is not None]
                self._total_space = BasisMultiStateSpace(np.array(space_list, dtype=object))
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

    def extend_state_spaces(self, new_targets):
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
                    flat = cs.to_single().take_unique()
                    existing_spaces[p] = ({None:cs}, cs)
                if self.state_space_filter_generator is not None:
                    filters = PerturbationTheoryStateSpaceFilter.from_data(new_targets, self.state_space_filter_generator(new_targets))
                    # raise Exception(filters[(1, 1)].prefilters)
                else:
                    filters = None
                new_spaces = self.load_coupled_spaces([new_targets],
                                                      filter_spaces=filters
                                                      )#, spaces=existing_spaces)
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
            self.space = space
        def __neg__(self):
            return self
        def simple_union(self, other):
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
            return self.simple_union(other)
        def __rsub__(self, other):
            return self.simple_union(other)
        def __add__(self, other):
            return self.simple_union(other)
        def __radd__(self, other):
            return self.simple_union(other)
        # def __mul__(self, other):
        #     return type(self)(self.space.intersection(other))
    class ProjectionOperatorWrapper:
        """
        Generates a symbolic form of a perturbation operator that
        either projects onto a degenerate space or projects it out
        """
        def __init__(self, space, complement=False):
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
        return (
            test is None
            or self._could_be_a_space(test)
            or (
                    len(test) == 2 and
                    test[0] is None or self._could_be_a_space(test[0])
                    and self._could_be_rules(test[1])
            )
        )
    def _apply_transformation_with_filters(self, a, b, filter_space:PerturbationTheoryStateSpaceFilter, **opts):

        if filter_space is not None:
            prefilters = filter_space.prefilters
            postfilter = filter_space.postfilter
        else:
            prefilters = None
            postfilter = None

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
            if postfilter is not None:
                if not isinstance(postfilter, (BasisStateSpace, BasisMultiStateSpace)):
                    postfilter = BasisStateSpace(b.basis, postfilter)
                if isinstance(postfilter, SelectionRuleStateSpace):
                    postfilter = postfilter.to_single().take_unique()
                new = new.take_states(postfilter)

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
                        b_sels = SelectionRuleStateSpace(b, [], ignore_shapes=True)  # just some type fuckery
                        existing = cur.intersection(b_sels, handle_subspaces=False)
                        # and now we do extra transformations where we need to


                        if len(diffs) > 0:
                            new_new = self._apply_transformation_with_filters(
                                a, diffs, filter_space,
                                track_excitations=not self.memory_constrained,
                                parallelizer=self.parallelizer, logger=logger
                            )

                            if new_new is None:
                                new = b
                            else:
                                # next we add the new stuff to the cache
                                cur = cur.union(new_new)
                                projections[proj] = rep_space.union(diffs)
                                spaces[op] = (projections, cur)

                                # TODO: find a way to make this not cause memory spikes...
                                if ret_space:
                                    new = existing.union(new_new).to_single(track_excitations=not self.memory_constrained).take_unique()
                                else:
                                    new = b
                        else:
                            new = b

                    else:
                        # means we already calculated everything
                        # so we don't need to worry about this
                        if cur is not None:
                            b_sels = SelectionRuleStateSpace(b, [], ignore_shapes=True) # just some type fuckery
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
                                                    filter_space=filter_spaces[(k-i, i)] if filter_spaces is not None and (k-i, i) in filter_spaces else None
                                                    )
                                else:
                                    dot(H[k - i], corrs[i], ret_space=False,
                                                    filter_space=filter_spaces[(k-i, i)] if filter_spaces is not None and (k-i, i) in filter_spaces else None
                                                    )
        else:
            raise NotImplementedError("property filters not here yet")

        return spaces

    #endregion

    #region Apply Equations

    def get_corrections(self,
                        non_zero_cutoff=None,
                        handle_strong_couplings=None,
                        check_overlap=True
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

        if isinstance(self.degeneracy_spec, DegeneracySpec):
            if self.degeneracy_spec.group_filter is not None:
                group_filter = self.get_degenerate_group_filter(group_filter=self.degeneracy_spec.group_filter)
            else:
                group_filter = None
            degenerate_states = DegenerateMultiStateSpace.from_spec(
                None if handle_strong_couplings else self.degeneracy_spec,
                solver=self,
                full_basis=self.full_basis,
                group_filter=group_filter,
                log_groups=True
            )#self.degeneracy_spec.get_groups(states, solver=self)
        # degenerate_states = None,
        # degeneracy_mode = None,
        # logger = None,
        # checkpointer = None,

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
        fundamentals = BasisStateSpace.from_quanta(self.total_state_space.basis, [0, 1])
        where = self.flat_total_space.find(fundamentals)
        test_freqs = self.zero_order_energies[where[1:]] - self.zero_order_energies[where[0]]
        test_modes = np.where(test_freqs > self.low_frequency_mode_cutoff)
        if len(test_modes) > 0:
            test_modes = list(test_modes[0])
        return test_modes

    default_strong_coupling_threshold = .3
    def identify_strong_couplings(self, corrs):
        spec = self.degeneracy_spec
        if spec is None:
            spec = DegeneracySpec.from_spec('auto')
        return spec.identify_strong_couplings(
            self,
            corrs
        ), spec.wfc_threshold

    def get_degenerate_group_filter(self, corrs=None, threshold=None, group_filter=None):
        return self.degeneracy_spec.get_degenerate_group_filter(
            self,
            corrs=corrs,
            threshold=threshold
        )

    def construct_strong_coupling_spaces(self, sc, corrs, states, threshold):
        sc = corrs.collapse_strong_couplings(sc)
        with self.logger.block(tag='unpruned spaces states', log_level=self.logger.LogLevel.Never):
            for s in sc.values():
                self.logger.log_print(
                    s,
                    message_prepper=lambda a: str(a.excitations).splitlines(),
                    log_level=self.logger.LogLevel.Never
                )

        group_filter = self.get_degenerate_group_filter(corrs=corrs, threshold=threshold)
        degenerate_states = DegenerateMultiStateSpace.from_spec(
            self.degeneracy_spec,
            couplings=sc,
            solver=self,
            group_filter=group_filter,
            log_groups=True
        )

        if self.degeneracy_spec.extend_spaces:
            missing = degenerate_states.to_single().difference(states)
            if len(missing) > 0:
                with self.logger.block(tag='new states'):
                    self.logger.log_print(
                        missing.excitations,
                        message_prepper=lambda a: str(np.array(a)).splitlines()
                    )
                extended_spaces = self.extend_state_spaces(missing)
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

        return degenerate_states, (states, perturbations, flat_total_space, N)

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

        with checkpointer:

            all_energies = np.zeros((len(states), order + 1))
            all_overlaps = np.zeros((len(states), order + 1))
            all_corrs = np.zeros((len(states), order + 1, N))
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
                    d=len([1 for x in degenerate_states if len(x) > 1]),
                    dm=(
                        'Sakurai' if self.allow_sakurai_degs else
                        'mod. H' if self.drop_perturbation_degs else
                        'standard'
                    )
                )
                start = time.time()

                _ = []
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
                            all_energies[res_index] = energies
                            all_energy_corrs[res_index] = ecorrs
                            all_corrs[res_index] = corrs
                            all_overlaps[res_index] = overlaps
                else:
                    # loop over the degenerate sets
                    for deg_group in degenerate_states:
                        # logger.log_print(str(deg_group.excitations))
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

                        # now we reorthogonalize degenerate states
                        if not self.intermediate_normalization:
                            for k in range(1, order+1):
                                for i, (n, x) in enumerate(zip(res_inds, deg_inds)):
                                    for m, y in zip(res_inds[i + 1:], deg_inds[i+1:]):
                                        onn = -1 / 2 * np.sum(np.dot(all_corrs[n][i], all_corrs[n][k - i]) for i in range(1, k))
                                        onm = -1 / 2 * np.sum(np.dot(all_corrs[n][i], all_corrs[m][k - i]) for i in range(1, k))
                                        omm = -1 / 2 * np.sum(np.dot(all_corrs[m][i], all_corrs[m][k - i]) for i in range(1, k))
                                        all_corrs[n][k][x] = onn
                                        all_corrs[n][k][y] = onm
                                        all_corrs[m][k][x] = onm
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
                logger=logger
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
            )

            if (
                    degenerate_states is None
                    or all(len(x) == 1 for x in degenerate_states)
            ):
                sc, degenerate_correction_threshold = self.identify_strong_couplings(corrs)
                if len(sc) > 0:
                    with self.logger.block(tag="Strongly coupled states (threshold={})".format(degenerate_correction_threshold)):
                        self.logger.log_print(
                            corrs.format_strong_couplings_report(sc, join=False)
                        )

                        if handle_strong_couplings and (sc is not None and len(sc) > 0):
                            degenerate_states, meta = self.construct_strong_coupling_spaces(sc, corrs, states, degenerate_correction_threshold)
                            states, perturbations, flat_total_space, N = meta
                            with self.logger.block(tag="Redoing PT with strong couplings handled"):
                                with self.logger.block(tag="Degenerate groups:"):
                                    for g in degenerate_states.flat:
                                        if len(g) > 1:
                                            self.logger.log_print(str(g.excitations).splitlines())

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
                        "total_states": total_states.excitations,
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
                            "total_states": total_states.excitations,
                            'energies': all_energies,
                            'wavefunctions': corr_mats
                        }
                    except KeyError:
                        pass

            # with self.logger.block(tag="overlap matrix"):
            #     self.logger.log_print(str(np.sum(corrs.get_overlap_matrices(), axis=0)).splitlines())

        return self.PTResults(corrs, degenerate_states)

    PTResults = namedtuple("PTResults", ["corrections", "degeneracies"])

    @staticmethod
    def _safe_dot(a, b):
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

    def apply_VPT_nondeg_equations(self,
                                   state_index,
                                   deg_group,
                                   perturbations=None,
                                   non_zero_cutoff=None,
                                   check_overlap=True,
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
            total_state_space.take_subspace([n_ind]).excitations[0]
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
                                  e=E0[0] * UnitsData.convert("Hartrees", "Wavenumbers"),
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
                    if (should_be_zero > 0).any():
                        raise ValueError("Perturbation operator should have made overlap of state {} with {} zero...got {} instead".format(
                            n, D, should_be_zero
                        ))

                if intermediate_normalization:
                    ok = 0.0
                else:
                    ok = -1 / 2 * np.sum(dot(corrs[i], corrs[k - i]) for i in range(1, k))
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
    def _fmt_depert_engs(total_state_space, base_energies):

        def fmt(*a, base_energies=base_energies):
            states = total_state_space.excitations
            state_sums = np.sum(total_state_space.excitations, axis=1)
            gs_pos = np.where(state_sums == 0)
            label = "Deperturbed States/Energies:"
            if len(gs_pos) > 0:
                label = "Deperturbed States/Frequencies:"
                gs_pos = gs_pos[0]
                if len(gs_pos) > 0:
                    gs_pos = gs_pos[0]
                gs_eng = base_energies[gs_pos]
                base_energies = base_energies - gs_eng
                base_energies[gs_pos] = gs_eng

            ndim = len(states[0])
            fmt_str = "{:.0f} " * ndim + " {:12.4f}"

            return [label] + [
                fmt_str.format(*s, e)
                for s, e in zip(
                    total_state_space.excitations,
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
            message_prepper=self._fmt_depert_engs(total_state_space, base_energies)
        )

        rotations = []
        rotation_vals = []
        rotation_row_inds = []
        rotation_col_inds = []

        ndeg_ham_corrs = []
        for group in degenerate_states:
            deg_inds, H_nd, deg_rot, deg_engs = corrs.get_degenerate_transformation(group, self.representations, gaussian_resonance_handling=self.gaussian_resonance_handling)
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
                        lambda *a:["dropping elements coupling degenerate space:"] + str(g.excitations).splitlines(),
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

    # def _martin_test(cls, h_reps, states, threshold, total_coupled_space):
    #     """
    #     Applies the Martin Test to a set of states and perturbations to determine which resonances need to be
    #     treated variationally. Everything is done within the set of indices for the representations.
    #
    #     :param h_reps: The representation matrices of the perturbations we're applying.
    #     :type h_reps: Iterable[np.ndarray | SparseArray]
    #     :param states: The indices of the states to which we're going apply to the Martin test.
    #     :type states: np.ndarray
    #     :param threshold: The threshold for what should be treated variationally (in the same energy units as the Hamiltonians)
    #     :type threshold: float
    #     :return: Pairs of coupled states
    #     :rtype: tuple[BasisStateSpace, BasisStateSpace]
    #     """
    #
    #     raise NotImplementedError("This is fucked up :weep:; need to do full non-degenerate calc per pair of states")
    #
    #     H0 = h_reps[0]
    #     H1 = h_reps[1]
    #     energies = np.diag(H0) if isinstance(H0, np.ndarray) else H0.diag
    #
    #     # the 'states' should already be indices within the space over which we do the H1 calculation
    #     # basically whichever states we need to treat as degenerate for
    #     state_energies = energies[states]
    #     diffs = state_energies[:, np.newaxis] - energies[np.newaxis, :]
    #     for n, s in enumerate(states):
    #         diffs[n, s] = 1
    #
    #     deg_states = []
    #     for s in states:
    #         # pull the blocks out of H1 that correspond to each the `states` we fed in...
    #         H1_block = H1[s, :]
    #         if isinstance(H1_block, SparseArray):
    #             nzvals = H1_block.block_vals
    #             nzinds, _ = H1_block.block_inds
    #             H1_block = nzvals
    #             diffs = energies[s] - energies[nzinds]  # do I need an abs ?
    #         else:
    #             # compute the energy differences
    #             diffs = energies[s] - energies  # do I need an abs ?
    #             nzinds = np.arange(len(energies))
    #
    #         s_pos = np.where(nzinds == s)[0]
    #         H1_block[s_pos] = 0
    #         diffs[s_pos] = 1
    #
    #         anh_eff = (np.abs(H1_block) ** 4) / (diffs ** 3)
    #         big = np.where(np.abs(anh_eff) > threshold)[0]
    #         if len(big) > 0:
    #             deg_states.extend((s, nzinds[d]) for d in big)
    #
    #     if len(deg_states) == 0:
    #         return None
    #     else:
    #         new_degs = np.array(deg_states).T
    #
    #         # raise Exception(new_degs)
    #
    #         # we now have indices inside the space of coupled states...
    #         # so now we need to broadcast these back into their indices in the overall basis of states
    #         tc_inds = total_coupled_space.indices
    #         basis = total_coupled_space.basis
    #
    #         degs = (
    #             BasisStateSpace(basis, tc_inds[new_degs[0]], mode='indices'),
    #             BasisStateSpace(basis, tc_inds[new_degs[1]], mode='indices')
    #         )
    #
    #         return degs

    # def _prep_degeneracies_spec(self, degeneracies):
    #     if (
    #             degeneracies is not None
    #             and not isinstance(degeneracies, (int, float, np.integer, np.floating))
    #     ):
    #         if isinstance(degeneracies[0], (int, np.integer)):
    #             degs = BasisStateSpace(self.basis, degeneracies, mode="indices")
    #             degeneracies = (degs, degs)
    #         elif isinstance(degeneracies[0][0], (int, np.integer)):
    #             degs = BasisStateSpace(self.basis, degeneracies, mode="excitations")
    #             degeneracies = (degs, degs)
    #         else:
    #             degeneracies = (
    #                 BasisStateSpace(self.basis, degeneracies[0]),
    #                 BasisStateSpace(self.basis, degeneracies[1])
    #             )
    #
    #     return degeneracies