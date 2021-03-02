
import numpy as np, itertools, time

from McUtils.Numputils import SparseArray
from McUtils.Scaffolding import Logger, NullLogger
from McUtils.Parallelizers import Parallelizer
from McUtils.Data import UnitsData

from ..BasisReps import BasisStateSpace, BasisMultiStateSpace, SelectionRuleStateSpace, BraKetSpace
from .Common import PerturbationTheoryException


__reload_hook__ = [ "..BasisReps" ]

class PerturbationTheoryCorrections:
    """
    Represents a set of corrections from perturbation theory.
    Can be used to correct other operators in the basis of the original calculation.

    """
    def __init__(self,
                 hamiltonians,
                 states,
                 coupled_states,
                 total_basis,
                 energy_corrs,
                 wfn_corrections,
                 degenerate_states=None,
                 degenerate_transformation=None,
                 degenerate_energies=None
                 ):
        """
        :param hamiltonians:
        :type hamiltonians: Iterable[SparseArray]
        :param states:
        :type states: BasisStateSpace
        :param coupled_states:
        :type coupled_states: BasisMultiStateSpace
        :param total_basis:
        :type total_basis: BasisMultiStateSpace
        :param energy_corrs:
        :type energy_corrs: np.ndarray
        :param wfn_corrections:
        :type wfn_corrections: Iterable[SparseArray]
        :param degenerate_states:
        :type degenerate_states: None | np.ndarray
        :param degenerate_transformation:
        :type degenerate_transformation: None | np.ndarray
        :param degenerate_energies:
        :type degenerate_energies: None | np.ndarray
        """
        self.hams = hamiltonians
        self.states = states
        self.coupled_states = coupled_states
        self.total_basis = total_basis
        self.energy_corrs = energy_corrs
        self.wfn_corrections = wfn_corrections
        self.degenerate_states = degenerate_states
        self.degenerate_transf = degenerate_transformation
        self.degenerate_energies = degenerate_energies

    @classmethod
    def from_dicts(cls,
                 states,
                 corrections,
                 hamiltonians
                 ):
        """
        :param states: a dict with the states described by the corrections, the set of states coupled, and the size of the overall basis
        :type states: dict
        :param corrections: the corrections generated, including the corrections for the energies, wavefunctions, and a transformation from degenerate PT
        :type corrections: dict
        :param hamiltonians: the set of Hamiltonian matrices used as an expansion
        :type hamiltonians: Iterable[np.ndarray]
        """
        state_space = states['states']
        coupled_states = states['coupled_states']
        total_basis = states['total_states']
        energy_corrs = corrections['energies']
        wfn_corrections = corrections['wavefunctions']
        if 'degenerate_states' in states:
            degenerate_states = states['degenerate_states']
        else:
            degenerate_states = None

        if 'degenerate_transformation' in corrections:
            degenerate_transf = corrections['degenerate_transformation']
        else:
            degenerate_transf = None

        if 'degenerate_energies' in corrections:
            degenerate_energies = corrections['degenerate_energies']
        else:
            degenerate_energies = None

        return cls(
            hamiltonians,
            state_space,
            coupled_states,
            total_basis,
            energy_corrs,
            wfn_corrections,
            degenerate_states=degenerate_states,
            degenerate_transformation=degenerate_transf,
            degenerate_energies=degenerate_energies
        )

    @property
    def degenerate(self):
        return self.degenerate_transf is not None

    @property
    def energies(self):
        if self.degenerate:
            return self.degenerate_energies
        else:
            return np.sum(self.energy_corrs, axis=1)

    @property
    def order(self):
        return len(self.energy_corrs[0])

    def take_subspace(self, space):
        """
        Takes only those elements that are in space
        :param space:
        :type space:
        :return:
        :rtype:
        """

        new_states = self.states.find(space)
        # raise Exception(new_states)
        return type(self)(
            self.hams,
            self.states.take_subspace(new_states),
            self.coupled_states.take_subspace(new_states),
            self.total_basis,
            self.energy_corrs[new_states],
            [w[new_states, :] for w in self.wfn_corrections],
            # not sure what to do with all this...
            degenerate_states=self.degenerate_states,
            degenerate_transformation=self.degenerate_transf,
            degenerate_energies=self.degenerate_energies
        )



    def operator_representation(self, operator_expansion, order=None, subspace=None):
        """
        Generates the representation of the operator in the basis of stored states

        :param operator_expansion: the expansion of the operator
        :type operator_expansion: Iterable[float] | Iterable[np.ndarray]
        :param order: the order of correction to go up to
        :type order: Iterable[float] | Iterable[np.ndarray]
        :param subspace: the subspace of terms in which the operator expansion is defined
        :type subspace: None | BasisStateSpace
        :return: the set of representation matrices for this operator
        :rtype: Iterable[np.ndarray]
        """

        mordor = self.order - 1
        if order is None:
            order = mordor
        if order > mordor:
            raise PerturbationTheoryException("{}: can't correct up to order {} when zero-order states were only corrected up to order {}".format(
                type(self).__name__,
                order,
                mordor
            ))
        order = order + 1 # so that we actually do get up to the request order after accounting for the zeros...
        if len(operator_expansion) < order:
            operator_expansion = list(operator_expansion) + [0]*(order - len(operator_expansion))

        # we stopped supporting indexing based on the total set of inds...
        if subspace is None:
            wfn_corrs = self.wfn_corrections[:order]
        else:
            # need to make the subspace good with the subspace in which the corrections are defined...
            subspace_sel = self.total_basis.find(subspace, check=True)
            wfn_corrs = []
            # raise Exception(subspace_sel)
            for k in range(order):
                wfn_corrs.append(self.wfn_corrections[k][:, subspace_sel])
        # raise Exception(wfn_corrs)

        # generalizes the dot product so that we can use 0 as a special value...
        dot = PerturbationTheorySolver._safe_dot

        # does the dirty work of acutally applying the rep...
        reps = [np.zeros(1)] * order
        for k in range(order):
            op = None
            # apply each thing up to requested order...
            for a in range(k+1):
                for b in range(k-a+1):
                    c = k - (a + b)
                    rop = operator_expansion[c]
                    if isinstance(rop, (int, float, np.integer, np.floating)): # constant reps...
                        if rop != 0: # cheap easy check
                            subrep = rop * dot(wfn_corrs[a], wfn_corrs[b].T)
                            if op is None:
                                op = subrep
                            else:
                                op += subrep
                    else:
                        subrep = dot(dot(wfn_corrs[a], rop), wfn_corrs[b].T)
                        if op is None:
                            op = subrep
                        else:
                            op += subrep
            reps[k] = op

        return reps

    def savez(self, file):
        keys = dict(
            states=self.states,
            coupled_states=self.coupled_states,
            total_states=self.total_basis,
            energies=self.energy_corrs,
            wavefunctions=self.wfn_corrections,
            hamiltonians=self.hams
        )
        if self.degenerate_states is not None:
            keys['degenerate_states'] = self.degenerate_states
        if self.degenerate_transf is not None:
            keys['degenerate_transformation'] = self.degenerate_transf
        if self.degenerate_energies is not None:
            keys['degenerate_energies'] = self.degenerate_energies
        np.savez(file, **keys)
    @classmethod
    def loadz(cls, file):
        keys = np.load(file)
        return cls(
            {
                "states":keys['states'],
                 "coupled_states":keys['coupled_states'],
                 "total_states":keys['total_states'],
                 "degenerate_states":keys['degenerate_states'] if 'degenerate_states' in keys else None
             },
            {
                "energies":keys['energies'],
                "wavefunctions":keys['wavefunctions'],
                "degenerate_transformation": keys['degenerate_transformation'] if 'degenerate_transformation' in keys else None,
                "degenerate_energies": keys['degenerate_energies'] if 'degenerate_energies' in keys else None
            },
            keys['hamiltonians']
        )

class DegenerateMultiStateSpace(BasisMultiStateSpace):

    @classmethod
    def from_spec(cls,
                  solver,
                  degenerate_states
                  ):
        """
        Generates a DegenerateMultiStateSpace object from a number
        of possible specs

        :param solver: the actual applier of the perturbation theory which makes use of the degenerate states
        :type solver: PerturbationTheorySolver
        :return:
        :rtype:
        """

        total_state_space = solver.flat_total_space
        logger = solver.logger
        H = solver.representations
        states = solver.states

        if degenerate_states is not None:

            with logger.block(tag="getting degeneracies"):
                if isinstance(degenerate_states, dict):
                    if 'MartinTest' in degenerate_states:
                        martin_test = degenerate_states['MartinTest']
                    if 'states' in degenerate_states:
                        degenerate_states = degenerate_states['states']
                    elif 'NT' in degenerate_states:
                        logger.log_print(
                            "NT vector: {s}",
                            s=degenerate_states
                        )
                        degenerate_states = cls._group_states_by_nt_spec(H, states, total_state_space,
                                                                         degenerate_states['NT'])
                    elif 'energy_cutoff' in degenerate_states:
                        logger.log_print(
                            "energy cutoff: {s}",
                            s=degenerate_states
                        )
                        degenerate_states = cls._group_states_by_energy_cutoff(H, states, total_state_space,
                                                                               degenerate_states['energy_cutoff'])
                    elif 'generator' in degenerate_states:
                        logger.log_print(
                            "callable: {s}",
                            s=degenerate_states
                        )

                        # we assume we have some kind of callable
                        try:
                            degenerate_states = degenerate_states['generator'](H, states)
                        except (TypeError, ValueError):
                            pass

                        if not isinstance(degenerate_states[0], (BasisStateSpace, BasisMultiStateSpace)):
                            raise NotImplementedError("can't deal with non-BasisStateSpace specs for degeneracies")

                    else:
                        raise NotImplementedError("unsure what to do with degeneracy spec {}".format(degenerate_states))

                else:

                    def _is_degenerate_NT_spec(spec):
                        test1 = isinstance(spec, np.ndarray) and spec.dtype == np.dtype(int)
                        if test1:
                            return test1
                        else:
                            try:
                                it = iter(spec)
                            except TypeError:
                                return False
                            else:
                                return all(isinstance(x, int) for x in it)

                    def _is_degenerate_state_spec(spec):
                        test1 = all(isinstance(x, (BasisStateSpace, BasisMultiStateSpace)) for x in spec)
                        if test1:
                            return test1
                        else:
                            spec = np.asarray(spec)
                            return (
                                    spec.dtype == int
                                    and (spec.ndim == 2 or spec.ndim == 3)
                            )

                    # we dispatch on the types of degeneracy specs we support
                    if isinstance(degenerate_states, (int, np.integer, np.floating, float)):
                        logger.log_print(
                            "energy cutoff: {s}",
                            s=degenerate_states
                        )
                        degenerate_states = cls._group_states_by_energy_cutoff(H, states, total_state_space,
                                                                               degenerate_states)
                    elif _is_degenerate_NT_spec(degenerate_states):
                        logger.log_print(
                            "N_T vector: {s}",
                            s=degenerate_states
                        )
                        degenerate_states = cls._group_states_by_nt_spec(H, states, total_state_space,
                                                                         degenerate_states)
                    elif _is_degenerate_state_spec(degenerate_states):
                        if not isinstance(degenerate_states[0], (BasisStateSpace, BasisMultiStateSpace)):
                            degenerate_states = [BasisStateSpace(states.basis, x) for x in degenerate_states]
                    elif degenerate_states is not None:
                        logger.log_print(
                            "callable: {s}",
                            s=degenerate_states
                        )

                        # we assume we have some kind of callable
                        try:
                            degenerate_states = degenerate_states(H, states)
                        except (TypeError, ValueError):
                            pass

                        if not isinstance(degenerate_states[0], (BasisStateSpace, BasisMultiStateSpace)):
                            raise NotImplementedError("can't deal with non-BasisStateSpace specs for degeneracies")

                logger.log_print(
                    "degenerate state sets found {n}",
                    n=len([x for x in degenerate_states if len(x) > 1])
                )

                # if martin_test:
                #     # need to apply Martin test to every pair of states to figure out if they are truly
                #     # going to be significantly affected by the near degeneracy
                #     raise NotImplementedError("Don't have Martin test applying cleanly yet")
                #     logger.log_print(
                #         "applying Martin test with threshold {}",
                #         thresh
                #     )
                #     degenerate_states = cls._martin_test(
                #         H,
                #         state_inds,  # state indices in the coupled_states
                #         thresh,
                #         total_state_space
                #     )
                # else:
                #     logger.log_print(
                #         "skipping Martin test"
                #     )

        # build groups of degenerate states for use later
        if degenerate_states is None:
            groups = [[x] for x in states.indices]  # we're gonna loop through this later so why not destructure now...
        else:
            groups = [[]] * len(degenerate_states)
            deg_sets = [set(d.indices) for d in degenerate_states]
            for x in states.indices:
                for i, d in enumerate(deg_sets):
                    if x in d:
                        if len(groups[i]) == 0:
                            groups[i] = []
                        groups[i].append(x)
                        break
                else:
                    groups.append([x])

        # now turn these into proper BasisStateSpace objects so we can work with them more easily
        groups = [total_state_space.take_states(g) for g in groups]

        return cls(np.array(groups))

    @classmethod
    def _group_states_by_energy_cutoff(cls, H, states, total_state_space, cutoff):
        """
        :type H: Iterable[SparseArray]
        :type states: BasisStateSpace
        :type total_state_space: BasisMultiStateSpace
        :type cutoff: float
        :rtype: Iterable[BasisStateSpace]
        """
        # we look for states with energies within a range...
        # so initially we pull the sets of energies
        energies = (H[0].diag if isinstance(H[0], SparseArray) else np.diag(H[0]))
        degenerate_groups = []
        # then we look through the input states
        for n in total_state_space.find(states):
            # we only want to apply this once per degenerate group
            # NOTE: this is a path to subtlety, since
            #   if state a is within 50 cm^-1 of state b, and state b is within of c,
            #   you might argue a and c are degenerate
            #   we are wagering that states are distinct _enough_ such that this is not
            #   an issue, but if it is a different strategy will be required
            if all(n not in d for d in degenerate_groups):
                e_diffs = np.abs(energies - energies[n])
                inds = np.where(e_diffs < cutoff)[0]
                degenerate_groups.append(set(inds))
        # raise Exception(degenerate_groups)
        degenerate_groups = [total_state_space.take_subspace(np.array(list(d), dtype=int)) for d in degenerate_groups]
        return degenerate_groups

    @classmethod
    def _group_states_by_nt_spec(cls, H, states, total_state_space, q_vec):
        """
        :type H: Iterable[SparseArray]
        :type states: BasisStateSpace
        :type total_state_space: BasisMultiStateSpace
        :type cutoff: Iterable[int]
        :rtype: Iterable[BasisStateSpace]
        """
        # we build the total N_t to compare against once...
        tot_n_t = np.dot(total_state_space.excitations, q_vec)
        degenerate_groups = {}
        # then we look through the input states
        for vec in states.excitations:
            # base n_t
            n_t = np.dot(q_vec, vec)
            if n_t not in degenerate_groups:
                degenerate_groups[n_t] = np.where(tot_n_t == n_t)[0]
        degenerate_groups = [total_state_space.take_subspace(np.array(d)) for d in degenerate_groups.values()]
        return degenerate_groups

class PerturbationTheorySolver:
    """
    A solver that applies perturbation theory
    given a series of corrections and population of states.
    Supports degenerate and non-degenerate PT.
    """

    def __init__(self, perturbations, states, coupled_states, order=2,
                 allow_sakurai_degs=True, allow_post_PT_calc=True,
                 degeneracy_spec=None, degeneracy_mode="basis",
                 logger=None, parallelizer=None, checkpointer=None):
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
        """
        self.perts = perturbations
        self._reps = None
        self.order = order

        self.logger = logger
        self.parallelizer = parallelizer
        self.checkpointer = checkpointer

        self.states = states

        self.degeneracy_spec = degeneracy_spec
        self._deg_states = None
        self.degeneracy_mode = degeneracy_mode
        self.allow_sakurai_degs = allow_sakurai_degs
        self.allow_post_PT_calc = allow_post_PT_calc

        if len(coupled_states) != len(perturbations) - 1:
            raise ValueError("coupled states must be specified for all perturbations (got {}, expected {})".format(
                len(coupled_states),
                len(perturbations) - 1
            ))

        space_list = [states] + list(coupled_states)
        self.total_state_space = BasisMultiStateSpace(np.array(space_list, dtype=object))
        self.flat_total_space = self.total_state_space.to_single().take_unique()
        self.total_space_dim = len(self.flat_total_space)

        self._zo_engs = None

    class PastIndexableTuple(tuple):
        def __getitem__(self, item):
            if item >= len(self):
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

    @property
    def degenerate_states(self):
        if self._deg_states is None:
            self._deg_states = DegenerateMultiStateSpace.from_spec(self, self.degeneracy_spec)
        return self._deg_states

    @property
    def zero_order_energies(self):
        if self._zo_engs is None:
            H0 = self.representations[0]
            e_vec_full = np.diag(H0) if isinstance(H0, np.ndarray) else H0.diag
            if isinstance(e_vec_full, SparseArray):
                e_vec_full = e_vec_full.toarray()
            self._zo_engs = e_vec_full
        return self._zo_engs

    def get_VPT_representations(self):
        """
        Gets the sparse representations of the passed perturbation inside the basis of coupled states.

        :return:
        :rtype: Iterable[SparseArray]
        """

        logger = self.logger
        if logger is None:
            logger = NullLogger()

        par = Parallelizer.lookup(self.parallelizer)
        with par:  # we put an outermost block here to just make sure everything is clean


            with logger.block(tag="generating total space"):
                start = time.time()

                diag_inds = BraKetSpace(self.flat_total_space, self.flat_total_space)
                N = len(self.flat_total_space)

                end = time.time()
                logger.log_print(
                    [
                        "total coupled space dimensions: {d}",
                        "took {t:.3f}s"
                    ],
                    d=N,
                    t=end - start
                )

            H = [np.zeros(1)] * len(self.perts)
            with logger.block(tag="getting H0"):
                start = time.time()
                logger.log_print(["calculating diagonal elements"])
                diag = self.perts[0][diag_inds]
                logger.log_print(["constructing sparse representation"])
                H[0] = SparseArray.from_diag(diag)
                end = time.time()
                logger.log_print("took {t:.3f}s", t=end - start)

            # print(flat_total_space.indices)
            for i, h in enumerate(self.perts[1:]):
                # calculate matrix elements in the coupled subspace
                cs = self.total_state_space[i + 1]
                with logger.block(tag="getting H" + str(i + 1)):
                    start = time.time()
                    H[i + 1] = self._build_representation_matrix(h, cs)
                    end = time.time()
                    logger.log_print("took {t:.3f}s", t=end - start)


            return H

    def _build_representation_matrix(self, h, cs):
        """
        Actively constructs a perturbation theory Hamiltonian representation

        :param h:
        :type h:
        :param cs:
        :type cs:
        :return:
        :rtype:
        """
        logger = self.logger
        m_pairs = cs.get_representation_brakets()# something for the future... freq_threshold=self.freq_threshold

        if len(m_pairs) > 0:
            logger.log_print(["coupled space dimension {d}"], d=len(m_pairs))
            sub = h[m_pairs]
            SparseArray.clear_caches()
        else:
            logger.log_print('no states to couple!')
            sub = 0

        logger.log_print("constructing sparse representation...")

        N = self.total_space_dim
        if isinstance(sub, (int, np.integer, np.floating, float)):
            if sub == 0:
                sub = SparseArray.empty((N, N), dtype=float)
            else:
                raise ValueError("Using a constant shift of {} will force Hamiltonians to be dense...".format(sub))
                sub = np.full((N, N), sub)
        else:
            # figure out the appropriate inds for this data in the sparse representation
            row_inds = self.flat_total_space.find(m_pairs.bras)
            col_inds = self.flat_total_space.find(m_pairs.kets)

            # upper triangle of indices
            up_tri = np.array([row_inds, col_inds]).T
            # lower triangle is made by transposition
            low_tri = np.array([col_inds, row_inds]).T
            # but now we need to remove the duplicates, because many sparse matrix implementations
            # will sum up any repeated elements
            full_inds = np.concatenate([up_tri, low_tri])
            full_dat = np.concatenate([sub, sub])

            _, idx = np.unique(full_inds, axis=0, return_index=True)
            sidx = np.sort(idx)
            full_inds = full_inds[sidx]
            full_dat = full_dat[sidx]
            sub = SparseArray.from_data((full_dat, full_inds.T), shape=(N, N))

        return sub

    def apply_VPT(self):
        """
        Applies perturbation theory to the held basis of states using the
        built representations and degenerate state spaces

        :return:
        :rtype: PerturbationTheoryCorrections
        """

        corrs = self.get_corrections()

        # import McUtils.Plots as plt

        # wat = []
        # for k in range(2+1):
        #     for i in range(k+1):
        #         c1 = corrs.wfn_corrections[i].asarray()
        #         c2 = corrs.wfn_corrections[k-i].asarray()
        #         wat.append(np.dot(c1, c2.T))
        # woof = np.sum(wat, axis=0)
        # plt.ArrayPlot(woof, plot_style=dict(vmin=-1e-2, vmax=1e-2)).show()

        # watttt = self.get_transformed_Hamiltonian(corrs)
        # plt.ArrayPlot(watttt,
        #               plot_style=dict(
        #                   vmin=-1e-6,
        #                   vmax=1e-6
        #               )
        #               ).show()

        # degeneracy_mode = self.degeneracy_mode
        degenerate_states = self.degenerate_states
        if (
                self.allow_post_PT_calc
                and degenerate_states is not None
                and any(len(x) > 1 for x in degenerate_states)
        ):
            with self.logger.block(tag="Applying post-PT variational calc."):
                deg_engs, deg_transf = self.apply_post_PT_variational_calc(degenerate_states, corrs)
                corrs.degenerate_energies = deg_engs
                corrs.degenerate_transf = deg_transf

        # tf = deg_transf.asarray()
        # tf_wat = np.dot(np.dot(tf, watttt), tf.T)
        #
        # plt.ArrayPlot(tf_wat).show()

        return corrs

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

        # this will be built from a series of block-diagonal matrices
        # so we store the relevant values and indices to compose the SparseArray
        rotation_vals = []
        rotation_row_inds = []
        rotation_col_inds = []

        for group in degenerate_states:
            # we apply the degenerate PT on a group-by-group basis
            # by transforming the H reps into the non-degenerate basis
            deg_inds = total_state_space.find(group)
            if len(deg_inds) > 1:
                deg_engs, deg_rot = self.get_degenerate_rotation(group, corrs)
                energies[deg_inds] = deg_engs
                rotation_vals.append(deg_rot.flatten())
                deg_rows, deg_cols = np.array([p for p in itertools.product(deg_inds, deg_inds)]).T
                rotation_row_inds.append(deg_rows)
                rotation_col_inds.append(deg_cols)
            else:
                # we'll be a little bit inefficient for now and speed up later
                energies[deg_inds[0]] = base_energies[deg_inds[0]]
                rotation_vals.append([1.])
                rotation_row_inds.append(deg_inds)
                rotation_col_inds.append(deg_inds)

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

        return energies, rotations

    def get_transformed_Hamiltonian(self, corrs, deg_group=None):
        if deg_group is not None:
            subcorrs = corrs.take_subspace(deg_group)
            inds = self.flat_total_space.find(deg_group)
            subhams = [SparseArray.from_data(self._take_subham(H, inds)) for H in self.representations]
            # H_nd =[
            #     x.asarray() if isinstance(x, SparseArray) else x
            #     for x in subcorrs.operator_representation(subhams, subspace=deg_group)
            # ]
            H_nd = np.sum([
                x.asarray() if isinstance(x, SparseArray) else x
                for x in subcorrs.operator_representation(subhams, subspace=deg_group)
            ],
                axis=0
            )
        else:
            subhams = self.representations
            H_nd = np.sum([
                x.asarray() if isinstance(x, SparseArray) else x
                for x in corrs.operator_representation(subhams)
            ],
                axis=0
            )
        return H_nd

    def get_degenerate_rotation(self, deg_group, corrs):

        logger = self.logger

        with logger.block(tag="states"):
            logger.log_print(
                str(
                    corrs.states.take_states(deg_group).excitations
                ).splitlines()
            )

        H_nd = self.get_transformed_Hamiltonian(corrs.take_subspace(deg_group), None)
        # H_nd = self.get_transformed_Hamiltonian(corrs, deg_group)

        with logger.block(tag="non-degenerate Hamiltonian"):
            logger.log_print(
                str(
                    np.round(H_nd * UnitsData.convert("Hartrees", "Wavenumbers"))
                ).splitlines()
            )

        deg_engs, deg_transf = np.linalg.eigh(H_nd)

        for i in range(len(deg_transf)):
            max_ov = np.max(deg_transf[:, i] ** 2)
            ov_thresh = .5
            if max_ov < ov_thresh:  # there must be a single mode that has more than 50% of the initial state character?
                logger.log_print(
                    "    state {i} is more than 50% mixed",
                    i=i
                )
            #     raise PerturbationTheoryException("mode {} is has no contribution of greater than {}".format(
            #         i, ov_thresh
            #     ))

        # we pick the terms with the max contribution from each input state
        # and zero out the contributions so that two states can't map
        # to the same input state
        sort_transf = np.abs(deg_transf.copy())
        sorting = [-1] * len(deg_transf)
        for i in range(len(deg_transf)):
            o = np.argmax(sort_transf[i, :])
            sorting[i] = o
            sort_transf[:, o] = 0.  # np.zeros(len(sort_transf))

        with logger.block(tag='contributions'):
            logger.log_print(
                str(np.round(100 * (deg_transf ** 2)).astype(int)).splitlines()
            )

        logger.log_print('sorting: {s}', s=sorting)

        # sorting = np.argsort(sorting)

        #
        # print(sorting)
        # # if len(sorting) != len(np.unique(sorting)):
        # #     raise PerturbationTheoryException("After diagonalizing can't distinguish modes...")
        deg_engs = deg_engs[sorting,]

        self.logger.log_print("degenerate energies {e}",
                         e=np.round(deg_engs * UnitsData.convert("Hartrees", "Wavenumbers")))

        deg_transf = deg_transf[:, sorting]

        return deg_engs, deg_transf

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
    def _build_subprojector(self, states, inds):
        """
        Builds a subspace projector where only inds will
        be included but also the projection will be onto states

        :param states: state vectors in the inds subspace
        :type states: np.ndarray
        :param inds: indices for the subspace
        :type inds: np.ndarray
        :return:
        :rtype: SparseArray
        """

        raise NotImplementedError("I'm not sure I need this?")

        shp = self.representations[0].shape
        ind_pairs = np.array(list(itertools.product(inds, inds))).T
        vals = np.dot(states.T, states)
        return SparseArray.from_data(
                (
                    np.ones(len(inds)),
                    (inds, inds)
                ),
                shape=shp
            )

    def _get_Pi0(self, degenerate_subspace, non_zero_cutoff=1.0e-14, E0=None):
        # generate the perturbation operator
        e_vec_full = self.zero_order_energies
        if E0 is None:
            E0 = np.average(e_vec_full[degenerate_subspace])
            # h2w = UnitsData.convert("Hartrees", "Wavenumbers")
            # gs = self.zero_order_energies[0] * h2w
            # raise Exception(E0*h2w - gs, e_vec_full[degenerate_subspace]*h2w - gs)
        e_vec = e_vec_full - E0
        e_vec[degenerate_subspace] = 1
        zero_checks = np.where(np.abs(e_vec) < non_zero_cutoff)[0]
        if len(zero_checks) > 0:
            bad_vec = np.concatenate([[E0], e_vec_full[zero_checks]])
            raise ValueError(
                "degeneracies encountered: states {} and {} other states are degenerate (average energy: {} stddev: {})".format(
                    degenerate_subspace,
                    len(zero_checks),
                    np.average(bad_vec),
                    np.std(bad_vec)
                ))
        pi = 1 / e_vec
        pi[degenerate_subspace] = 0
        return SparseArray.from_diag(pi)
    def _get_Pi1(self, D, G, E1, non_zero_cutoff=1.0e-14, singular_check=1e10):
        """
        Returns the first-order perturbation operator in the D_n subspace

        :param D: degenerate subspace
        :type D: np.ndarray[int]
        :param G: degenerate subspace within D corresponding to E1
        :type G: np.ndarray
        :param E1:
        :type E1: float
        :return:
        :rtype:
        """

        H1D = self._take_subham(self.representations[1], D)
        if np.max(np.abs(H1D)) < non_zero_cutoff:
            return np.zeros((len(D), len(D)))
        Pn = np.dot(G.T, G)
        Rn = np.eye(len(Pn)) - Pn
        dH = Rn@(H1D - E1 * np.eye(len(Pn)))@Rn

        pi = np.linalg.inv(dH)

        if np.max(np.abs(pi)) > singular_check:
            raise ValueError("singular perturbation operator in degenerate subspace {}".format(D))

        return pi
        # raise Exception(pi, Rn@pi@Rn)
    def _get_Pi2(self, D, G, n, E2, singular_check=1e10):
        """
        Returns the first-order perturbation operator in the D_n subspace

        :param D: degenerate subspace
        :type D: np.ndarray[int]
        :param G: degenerate subsubspace corresponding to E1
        :type G: np.ndarray
        :param E1:
        :type E1: float
        :return:
        :rtype:
        """

        H2D = self._take_subham(self.representations[2], D)
        Pg = np.dot(G.T, G)
        dH = Pg@(H2D - E2 * np.eye(len(Pg)))@Pg
        Rn = np.eye(len(Pg)) - np.dot(n[:, np.newaxis], n[np.newaxis, :])

        pi = np.linalg.inv(dH)
        # pi = np.dot(np.dot(Rn, pi), Rn)

        if np.max(np.abs(pi)) > singular_check:
            raise ValueError("singular perturbation operator in degenerate subspace {}".format(G))

        return pi
    def _get_V_projector(self, D, G):
        Pn = np.dot(G.T, G)
        Rn = np.eye(len(Pn)) - Pn

        ind_pairs = np.array(list(itertools.product(D, D))).T
        vals = Rn.flatten()

        shp = self.representations[0].shape
        return SparseArray.from_data(
            (
                vals,
                ind_pairs
            ),
            shape=shp
        )

    def _get_secondary_degenerate_inputs(self, deg_inds, deg_transf, subspace):
        """
        Gets the eigenvalues of Ps (H2  - H1 Pi_U H1) Ps where E_n is defined to be the
        average zero-order energy in s.
        The closer to degenerate the proper zero-order states are, the better this
        approximation (it is exact when the degeneracy is perfect)

        :param deg_inds:
        :type deg_inds:
        :param deg_transf:
        :type deg_transf:
        :param subspace: secondary degenerate subspace (I called it G in my OG stuff...)
        :type subspace:
        :return:
        :rtype:
        """

        PiU = self._get_Pi0(deg_inds)
        subtf = deg_transf[subspace, :]
        Ps = subtf.T@subtf

        H2 = self._take_subham(self.representations[2], deg_inds)
        H1UH1 = self._safe_dot(self._safe_dot(self.representations[1], PiU), self.representations[1])
        if not isinstance(H1UH1, (int, float, np.integer, np.floating)):
            H1UH1 = self._take_subham(H1UH1, deg_inds)

        subham = np.dot(np.dot(Ps, H2 - H1UH1), Ps)

        # import McUtils.Plots as plt
        # h1 = self._take_subham(self.representations[1], deg_inds)
        # plt.ArrayPlot(h1)
        # plt.ArrayPlot(H2)
        # plt.ArrayPlot(H1UH1)
        # plt.ArrayPlot(Ps)
        # plt.ArrayPlot(subham).show()

        eng2, deg_transf2 = np.linalg.eigh(subham)

        # raise Exception(eng2*UnitsData.convert("Hartrees", "Wavenumbers"))

        # now we need to get the appropriate sorting to match up the
        # secondary degenerate transformation and the OG terms
        overlaps = np.dot(deg_transf2.T, subtf.T)

        sort_transf = np.abs(overlaps)
        sorting = [-1] * len(deg_transf2)
        for i in range(len(deg_transf2)):
            o = np.argmax(sort_transf[i, :])
            sorting[i] = o
            sort_transf[:, o] = 0.
            # print(sort_transf)

        new_eng = eng2[sorting]
        deg_transf2 = deg_transf2.T[sorting]

        # raise Exception(sorting)

        # raise Exception(deg_transf2, new_eng, sorting)

        return new_eng, deg_transf2

    def _get_deg_eq_inputs(self, deg_inds, degeneracy_cutoff=.00001): # within a few wavenumbers or so
        """
        Diagonalizes the perturbations in the degenerate subspace to
        get a cleaner basis in which to do the perturbation theory.
        This comes from Sakurai.

        :param deg_inds:
        :type deg_inds:
        :return:
        :rtype: tuple[Iterable[SparseArray], SparseArray, SparseArray]
        """

        H1 = self.representations[1]
        subham = self._take_subham(H1, deg_inds)

        new_eng1, deg_transf = np.linalg.eigh(subham)

        # now we split this into degenerate subspaces by grouping up
        # runs of degenerate states
        deg_spaces = []
        deg_set = set()
        for i,a in enumerate(new_eng1):
            if i+1 < len(new_eng1):
                b_ind = i+1
            else:
                b_ind = i
            b = new_eng1[b_ind]
            if abs(a - b) < degeneracy_cutoff:
                deg_set.add(i)
                deg_set.add(b_ind)
            else:
                deg_spaces.append(deg_set)
                deg_set.add(i)
        deg_spaces.append(deg_set)
        deg_spaces = [np.sort(np.array(list(s))) for s in deg_spaces]
        # now we handle secondary degeneracies...

        subspaces = [ ]
        new_engs = [ ]
        for s in deg_spaces:
            if len(s) > 1:
                new_eng2, subs_transf = self._get_secondary_degenerate_inputs(deg_inds, deg_transf, s)
                deg_transf[:, s] = subs_transf.T

            else:
                new_eng2 = [None]
                subs_transf = None
            new_engs.extend(zip(new_eng1[s], new_eng2))
            subspaces.extend([subs_transf] * len(s))
        new_eng = np.full(len(new_engs), None)
        new_subspace = np.full(len(subspaces), None)
        for i in range(len(new_engs)):
            new_eng[i] = new_engs[i]
            new_subspace[i] = subspaces[i]

        # we sort now to get the "best" mapping back onto the OG states
        sort_transf = np.abs(deg_transf.copy())
        sorting = [-1] * len(deg_transf)
        for i in range(len(deg_transf)):
            o = np.argmax(sort_transf[i, :])
            sorting[i] = o
            sort_transf[:, o] = 0.  # np.zeros(len(sort_transf))

        new_eng = new_eng[sorting]
        new_subspace = new_subspace[sorting]
        deg_transf = deg_transf.T[sorting]

        return new_eng, deg_transf, new_subspace

    def get_corrections(self, non_zero_cutoff=1.0e-14):
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

        degenerate_states = self.degenerate_states
        # degenerate_states = None,
        # degeneracy_mode = None,
        # logger = None,
        # checkpointer = None,

        checkpointer['indices'] = self.total_state_space
        checkpointer['representations'] = perturbations

        all_energies = np.zeros((len(states), order + 1))
        all_overlaps = np.zeros((len(states), order + 1))
        all_corrs = np.zeros((len(states), order + 1, N))

        with logger.block(tag="applying non-degenerate PT"):
            logger.log_print(
                [
                    "order: {o}",
                    "states: {n}"
                ],
                o=order,
                n=len(states.indices)
            )
            start = time.time()

            # loop over the degenerate sets
            for deg_group in degenerate_states:
                # we use this to build a pertubation operator that removes
                # then entire set of degenerate states
                deg_inds = flat_total_space.find(deg_group)
                if hasattr(deg_group, 'indices'):
                    deg_group = deg_group.indices

                if len(deg_group) > 1:
                    deg_engs, zero_order_states, subspaces = self._get_deg_eq_inputs(deg_inds)
                else:
                    deg_engs = zero_order_states = subspaces = [None]
                for n, de, zo, s in zip(deg_group, deg_engs, zero_order_states, subspaces):
                    energies, overlaps, corrs = self.apply_VPT_equations(n, deg_inds, de, zo, s,
                                                                         allow_PT_degs=self.allow_sakurai_degs,
                                                                         non_zero_cutoff=non_zero_cutoff)

                    res_index = states.find(n)
                    all_energies[res_index] = energies
                    all_corrs[res_index] = corrs
                    all_overlaps[res_index] = overlaps

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

        corr_inds = [[] for i in range(nstates)]
        corr_mats = [None] * (order + 1)

        for o in range(order + 1):
            non_zeros = []
            for i, corr in enumerate(all_corrs):
                # we find the non-zero elements within the o level of correction for the ith state
                nonzi = np.where(np.abs(corr[o]) > non_zero_cutoff)[0]
                # print(nonzi)
                # then we pull these out
                vals = corr[o][nonzi,]
                # and we add the values and indices to the list
                non_zeros.append(
                    (
                        vals,
                        np.column_stack([
                            np.full(len(nonzi), i),
                            nonzi
                        ])
                    )
                )

                # and then we add the appropriate basis indices to the list of basis data
                wat = tci[nonzi,]
                corr_inds[i].append(wat)

            # now we build the full mat rep for this level of correction
            vals = np.concatenate([x[0] for x in non_zeros])
            inds = np.concatenate([x[1] for x in non_zeros], axis=0).T
            # print(inds, N)
            corr_mats[o] = SparseArray.from_data(
                (
                    vals,
                    inds
                ),
                shape=(nstates, N)
            )

        # now we build state reps from corr_inds
        for i, dat in enumerate(corr_inds):
            cat = np.concatenate(dat)
            _, upos = np.unique(cat, return_index=True)
            full_dat = cat[np.sort(upos)]
            corr_inds[i] = flat_total_space.take_states(
                full_dat)  # BasisStateSpace(states.basis, full_dat, mode="indices")

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
                "wavefunctions": corr_mats,
                "degenerate_transformation": None,
                "degenerate_energies": None
            },
            perturbations # we probably want to ditch this for memory reasons...
        )

        return corrs

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
                            degenerate_subsubspace,
                            allow_PT_degs=True,
                            non_zero_cutoff=1.0e-14
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
        if not allow_PT_degs:
            return self.apply_VPT_nondeg_equations(state_index, degenerate_space_indices, non_zero_cutoff=non_zero_cutoff)
        if len(degenerate_space_indices) == 1:
            return self.apply_VPT_nondeg_equations(state_index, None, non_zero_cutoff=non_zero_cutoff)
        elif len(degenerate_subsubspace[0]) == 1:
            return self.apply_VPT_deg_equations(state_index, degenerate_space_indices, degenerate_energies[0],
                                                zero_order_state, non_zero_cutoff=non_zero_cutoff)
        else:
            return self.apply_VPT_second_deg_equations(state_index, degenerate_space_indices, degenerate_energies,
                                                zero_order_state, degenerate_subsubspace, non_zero_cutoff=non_zero_cutoff)
    def apply_VPT_nondeg_equations(self,
                                   state_index,
                                   degenerate_space_indices,
                                   non_zero_cutoff=1.0e-14,
                                   check_overlap=True
                                   ):
        """
        Does the dirty work of doing the VPT iterative equations.
        Needs to be adapted to include the two types of degeneracies that can
        be introduced in Sakurai's approach.

        :return:
        :rtype:
        """

        n = state_index
        e_vec_full = self.zero_order_energies

        order = self.order
        total_state_space = self.flat_total_space

        energies = np.zeros((order + 1,), dtype=float)
        overlaps = np.zeros((order + 1,), dtype=float)
        corrs = np.zeros((order + 1, len(total_state_space)), dtype=float)  # can I make this less expensive in general?

        # find the state index in the coupled subspace
        n_ind = total_state_space.find(n)
        D = degenerate_space_indices
        if D is None:
            D = (n_ind,)
        E0 = e_vec_full[n_ind]
        pi = self._get_Pi0(D, non_zero_cutoff=non_zero_cutoff)
        energies[0] = E0
        overlaps[0] = 1
        corrs[0, n_ind] = 1
        H = self.representations
        dot = self._safe_dot
        take = lambda h, *els: h[els] if not isinstance(h, (int, np.integer, float, np.floating)) else 0.
        for k in range(1, order + 1):  # to actually go up to k
            #         En^(k) = <n^(0)|H^(k)|n^(0)> + sum(<n^(0)|H^(k-i)|n^(i)> - E^(k-i)<n^(0)|n^(i)>, i=1...k-1)
            Ek = (
                    take(H[k], n_ind, n_ind)
                    + sum(dot(take(H[k - i], n_ind), corrs[i]) - energies[k - i] * overlaps[i]
                for i in range(1, k)
            ))
            energies[k] = Ek
            #   <n^(0)|n^(k)> = -1/2 sum(<n^(i)|n^(k-i)>, i=1...k-1)
            #         |n^(k)> = sum(Pi_n (En^(k-i) - H^(k-i)) |n^(i)>, i=0...k-1) + <n^(0)|n^(k)> |n^(0)>
            corrs[k] = sum(
                dot(pi, energies[k - i] * corrs[i] - dot(H[k - i], corrs[i]))
                for i in range(0, k)
            )
            ok = -1 / 2 * np.sum(dot(corrs[i], corrs[k - i]) for i in range(1, k))
            overlaps[k] = ok
            corrs[k][n_ind] = ok  # pi (the perturbation operator) e nsures it's zero before this

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

        return energies, overlaps, corrs
    def apply_VPT_deg_equations(self,
                                state_index,
                                degenerate_space_indices,
                                degenerate_energy,
                                zero_order_state,
                                non_zero_cutoff=1.0e-14,
                                check_overlap=True
                                ):
        """
        Does the dirty work of doing the VPT iterative equations.
        Needs to be adapted to include the two types of degeneracies that can
        be introduced in Sakurai's approach.

        :return:
        :rtype:
        """

        n = state_index
        e_vec_full = self.zero_order_energies

        order = self.order
        total_state_space = self.flat_total_space

        energies = np.zeros((order + 2,), dtype=float)
        overlaps = np.zeros((order + 1,), dtype=float)
        corrs = np.zeros((order + 1, len(total_state_space)), dtype=float)  # can I make this less expensive in general?

        # find the state index in the coupled subspace
        n_ind = total_state_space.find(n)
        E0 = e_vec_full[n_ind]
        E1 = degenerate_energy

        D = degenerate_space_indices
        piU = self._get_Pi0(D, E0=E0, non_zero_cutoff=non_zero_cutoff)
        piDn = self._get_Pi1(D, zero_order_state[np.newaxis, :], E1)

        energies[0] = E0
        energies[1] = E1
        corrs[0, D] = zero_order_state

        H = self.representations
        dot = self._safe_dot
        for k in range(1, order + 1):  # to actually go up to k
            # Pu |n^(k)>
            corrs[k] = sum(
                dot(piU, energies[k - i] * corrs[i] - dot(H[k - i], corrs[i]))
                for i in range(1, k) # goes up to k-1
            ) - dot(piU, dot(H[k], corrs[0]))
            # E^(k+1)
            Ekp1 = (
                    dot(corrs[0], dot(H[k+1], corrs[0]))
                    + dot(corrs[0], dot(H[1], corrs[k]))
                    + sum(
                          dot(corrs[0], dot(H[k - i], corrs[i])) - energies[k - i] * overlaps[i]
                          for i in range(1, k)
                          )
            )
            energies[k+1] = Ekp1
            # PDn|n^(k)>
            dHDnk = sum(
                energies[k + 1 - i] * corrs[i] - dot(H[k + 1 - i], corrs[i])
                for i in range(0, k)
            ) - dot(H[1], corrs[k])
            corrs[k][D] = dot(piDn, dHDnk[D])

            # #   <n^(0)|n^(k)> = -1/2 sum(<n^(i)|n^(k-i)>, i=1...k-1)
            ok = -1 / 2 * sum(dot(corrs[i], corrs[k - i]) for i in range(1, k))
            overlaps[k] = ok
            corrs[k][D] += zero_order_state * ok

        if check_overlap:
            full_wfn = np.sum(corrs, axis=0)
            ov = np.dot(full_wfn, full_wfn)
            if abs(ov - 1) > .2:
                raise ValueError("state {} isn't normalized (overlap = {})".format(state_index, ov))


        return energies, overlaps, corrs

    def apply_VPT_second_deg_equations(self,
                                       state_index,
                                       degenerate_space_indices,
                                       degenerate_energies,
                                       zero_order_state,
                                       degenerate_subsubspace,
                                       non_zero_cutoff=1.0e-14,
                                       check_overlap=True
                                       ):
        """
        Does the dirty work of doing the VPT iterative equations.
        Needs to be adapted to include the two types of degeneracies that can
        be introduced in Sakurai's approach.

        :return:
        :rtype:
        """
        n = state_index
        e_vec_full = self.zero_order_energies

        order = self.order
        total_state_space = self.flat_total_space

        energies = np.zeros((order + 3,), dtype=float)
        overlaps = np.zeros((order + 1,), dtype=float)
        corrs = np.zeros((order + 2, len(total_state_space)), dtype=float)  # can I make this less expensive in general?

        # find the state index in the coupled subspace
        n_ind = total_state_space.find(n)
        E0 = e_vec_full[n_ind]
        E1 = degenerate_energies[0]
        E2 = degenerate_energies[1]

        D = degenerate_space_indices
        piU = self._get_Pi0(D, E0=E0, non_zero_cutoff=non_zero_cutoff)
        piV = self._get_Pi1(D, degenerate_subsubspace, E1)
        piG = self._get_Pi2(D, degenerate_subsubspace, zero_order_state, E2)
        # projV = self._get_V_projector(D, degenerate_subsubspace)

        energies[0] = E0
        energies[1] = E1
        energies[2] = E2
        corrs[0, D] = zero_order_state

        overlaps[0] = 1

        H = self.representations
        dot = self._safe_dot
        E = energies

        def eval_wf_bits(k):
            corrs[k] = sum(
                dot(piU, E[k - i] * corrs[i] - dot(H[k - i], corrs[i]))
                for i in range(1, k)  # goes up to k-1
            ) - dot(piU, dot(H[k], corrs[0]))
            # Pv |n^(k)>
            dHV = sum(
                energies[k + 1 - i] * corrs[i] - dot(H[k + 1 - i], corrs[i])
                for i in range(0, k)
            ) - dot(H[1], corrs[k])
            corrs[k][D] = dot(piV, dHV[D])
            watt = dot(corrs[0], corrs[k])
            if np.abs(watt) > 0.005:
                raise ValueError(
                    "state {} isn't orthogonal to V at order {} (contrib {})".format(
                        state_index, order, watt
                    )
                )

        eval_wf_bits(1)
        # eval_wf_bits(2)
        # raise Exception(piV)#dot(projV, corrs[1])[D])
        # c1 = corrs[1].copy()
        # c1[D] = 0.
        # raise Exception(dot(c1, c1))

        for k in range(1, order + 1):  # to actually go up to k
            # Pu |n^(k)>
            eval_wf_bits(k+1)

            # E^(k+2)
            Ecorrkp2 = (
                dot(corrs[0], dot(H[2], corrs[k]))
                - dot(corrs[0],
                      dot(dot(dot(H[1], piU), H[1]), corrs[k]))
                + dot(corrs[0],
                      dot(H[1], dot(piU,
                              (
                                  sum(
                                      E[k+1-i] * corrs[i] - dot(H[k+1-i], corrs[i])
                                      for i in range(1, k)
                                  )
                                  - dot(H[k+1], corrs[0])
                                  + E1*corrs[k]
                              )
                              )))
            )
            Ekp2 = (
                    dot(corrs[0], dot(H[k + 2], corrs[0]))
                    + Ecorrkp2
                    + sum(
                        dot(corrs[0], dot(H[k + 2 - i], corrs[i]))
                        - energies[k + 2 - i] * overlaps[i]
                        for i in range(1, k)
                    )
            )
            energies[k + 2] = Ekp2

            # PDn|n^(k)>
            dHG1_bits = [
                energies[k+2 - i] * corrs[i] - dot(H[k+2 - i], corrs[i])
                for i in range(0, k)
            ]
            dHG1 = sum(dHG1_bits)

            dHG2 = sum(
                dot(H[j], corrs[k + 2 - j])
                for j in range(1, 3)
            )

            # testFleh = sum(
            #     dot(piG, dot(projV, corrs[k + 2 - j])[D])
            #     for j in range(1, 3)
            # )

            dHG = dHG1 - dHG2
            ketG = dot(piG, dHG[D])
            # ketG1 = dot(piG, dHG1[D])
            # ketG1_bits = [dot(piG, x[D]) for x in dHG1_bits]
            # ketG2 = dot(piG, dHG2[D])
            corrs[k][D] += ketG

            # <n^(0)|n^(k)> = -1/2 sum(<n^(i)|n^(k-i)>, i=1...k-1)
            ok = -1 / 2 * sum(dot(corrs[i], corrs[k - i]) for i in range(1, k))
            corrs[k][D] += corrs[0][D]*ok - corrs[0][D]*dot(corrs[0][D], corrs[k][D])

            check_overlap = False
            if check_overlap:
                true_ov = dot(corrs[0], corrs[k])
                overlaps[k] = true_ov
                if abs(true_ov - ok) > .005:
                    # raise ValueError(
                    #     dot(corrs[0][D], ketG1),
                    #     dot(corrs[0][D], ketG2),
                    #     [dot(corrs[0][D], x[D]) for x in dHG1_bits],
                    #     [dot(corrs[0][D], b) for b in ketG1_bits],
                    #     # dot(projV, corrs[k ])[D],
                    #     # dot(projV, corrs[k + 1])[D],
                    #     ok
                    # )
                    raise ValueError(
                        "state {} fails overlap relationship at order {} (contrib {}, expected {})".format(
                            state_index, order, true_ov, ok
                        ))

            # overlaps[k] = ok
            # corrs[k][D] += zero_order_state * ok

        if check_overlap:
            # full_wfn = np.sum(corrs, axis=0)
            # ov = np.dot(full_wfn, full_wfn)
            ov_parts = [[dot(corrs[k - i], corrs[i]) for i in range(k + 1)] for k in range(order + 1)]
            ov = sum(np.sum(v) for v in ov_parts)
            if abs(ov - 1) > .005:
                raise ValueError(
                    "state {} isn't normalized (overlap = {}, bits {})".format(
                        state_index, ov, ov_parts
                    ))

        return energies[:order+1], overlaps, corrs[:order+1]

    def _martin_test(cls, h_reps, states, threshold, total_coupled_space):
        """
        Applies the Martin Test to a set of states and perturbations to determine which resonances need to be
        treated variationally. Everything is done within the set of indices for the representations.

        :param h_reps: The representation matrices of the perturbations we're applying.
        :type h_reps: Iterable[np.ndarray | SparseArray]
        :param states: The indices of the states to which we're going apply to the Martin test.
        :type states: np.ndarray
        :param threshold: The threshold for what should be treated variationally (in the same energy units as the Hamiltonians)
        :type threshold: float
        :return: Pairs of coupled states
        :rtype: tuple[BasisStateSpace, BasisStateSpace]
        """

        raise NotImplementedError("This is fucked up :weep:; need to do full non-degenerate calc per pair of states")

        H0 = h_reps[0]
        H1 = h_reps[1]
        energies = np.diag(H0) if isinstance(H0, np.ndarray) else H0.diag

        # the 'states' should already be indices within the space over which we do the H1 calculation
        # basically whichever states we need to treat as degenerate for
        state_energies = energies[states]
        diffs = state_energies[:, np.newaxis] - energies[np.newaxis, :]
        for n, s in enumerate(states):
            diffs[n, s] = 1

        deg_states = []
        for s in states:
            # pull the blocks out of H1 that correspond to each the `states` we fed in...
            H1_block = H1[s, :]
            if isinstance(H1_block, SparseArray):
                nzvals = H1_block.block_vals
                nzinds, _ = H1_block.block_inds
                H1_block = nzvals
                diffs = energies[s] - energies[nzinds]  # do I need an abs ?
            else:
                # compute the energy differences
                diffs = energies[s] - energies  # do I need an abs ?
                nzinds = np.arange(len(energies))

            s_pos = np.where(nzinds == s)[0]
            H1_block[s_pos] = 0
            diffs[s_pos] = 1

            anh_eff = (np.abs(H1_block) ** 4) / (diffs ** 3)
            big = np.where(np.abs(anh_eff) > threshold)[0]
            if len(big) > 0:
                deg_states.extend((s, nzinds[d]) for d in big)

        if len(deg_states) == 0:
            return None
        else:
            new_degs = np.array(deg_states).T

            # raise Exception(new_degs)

            # we now have indices inside the space of coupled states...
            # so now we need to broadcast these back into their indices in the overall basis of states
            tc_inds = total_coupled_space.indices
            basis = total_coupled_space.basis

            degs = (
                BasisStateSpace(basis, tc_inds[new_degs[0]], mode='indices'),
                BasisStateSpace(basis, tc_inds[new_degs[1]], mode='indices')
            )

            return degs

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