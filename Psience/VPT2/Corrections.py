
import numpy as np, itertools

from McUtils.Numputils import SparseArray
import McUtils.Numputils as nput
from McUtils.Data import UnitsData
from McUtils.Scaffolding import NullLogger, Checkpointer

from ..BasisReps import BasisStateSpace, BasisMultiStateSpace, SelectionRuleStateSpace
from .Common import PerturbationTheoryException, _safe_dot

__all__ = [
    "PerturbationTheoryCorrections"
]



class PerturbationTheoryCorrections:
    """
    Represents a set of corrections from perturbation theory.
    Can be used to correct other operators in the basis of the original calculation.

    """
    def __init__(self,
                 states,
                 coupled_states,
                 total_basis,
                 energy_corrs,
                 wfn_corrections,
                 all_energy_corrections=None,
                 degenerate_states=None,
                 degenerate_transformation=None,
                 degenerate_energies=None,
                 degenerate_hamiltonians=None,
                 logger=None
                 ):
        """
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
        self.states = states
        self.coupled_states = coupled_states
        self.total_basis = total_basis
        self.energy_corrs = energy_corrs
        self.all_energy_corrs = all_energy_corrections
        self.wfn_corrections = wfn_corrections
        self.degenerate_states = degenerate_states
        self.degenerate_transf = degenerate_transformation
        self.degenerate_energies = degenerate_energies
        self.degenerate_hamiltonians = degenerate_hamiltonians
        self.logger = logger

    @classmethod
    def from_dicts(cls,
                   states,
                   corrections,
                   **opts
                   ):
        """
        :param states: a dict with the states described by the corrections, the set of states coupled, and the size of the overall basis
        :type states: dict
        :param corrections: the corrections generated, including the corrections for the energies, wavefunctions, and a transformation from degenerate PT
        :type corrections: dict
        """
        state_space = states['states']
        coupled_states = states['coupled_states']
        total_basis = states['total_states']
        energy_corrs = corrections['energies']
        all_energy_corrs = corrections['energy_corrections'] if 'energy_corrections' in corrections else None
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
            state_space,
            coupled_states,
            total_basis,
            energy_corrs,
            wfn_corrections,
            all_energy_corrections=all_energy_corrs,
            degenerate_states=degenerate_states,
            degenerate_transformation=degenerate_transf,
            degenerate_energies=degenerate_energies,
            **opts
        )

    @property
    def degenerate(self):
        """

        :return:
        :rtype:
        """
        return self.degenerate_transf is not None

    @property
    def energies(self):
        """

        :return:
        :rtype:
        """
        if self.degenerate:
            return self.degenerate_energies
        else:
            return np.sum(self.energy_corrs, axis=1)

    @property
    def order(self):
        """

        :return:
        :rtype:
        """
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
        # print("? =", new_states)
        return type(self)(
            self.states.take_subspace(new_states),
            self.coupled_states.take_states(space),
            self.total_basis,
            self.energy_corrs[new_states],
            [w[new_states, :] for w in self.wfn_corrections],
            # not sure what to do with all this...
            degenerate_states=self.degenerate_states,
            degenerate_transformation=self.degenerate_transf,
            degenerate_energies=self.degenerate_energies,
            logger=self.logger
        )

    @classmethod
    def create_coupling_matrix(cls, corrs,
                               states:BasisStateSpace, flat_total_space:BasisStateSpace,
                               nstates, order,
                               filters=None,
                               non_zero_cutoff=1.0e-14,
                               logger=None
                               ):
        """

        :param corrs:
        :type corrs:
        :param states:
        :type states:
        :param flat_total_space:
        :type flat_total_space:
        :param nstates:
        :type nstates:
        :param order:
        :type order:
        :param filters:
        :type filters:
        :param non_zero_cutoff:
        :type non_zero_cutoff:
        :return:
        :rtype:
        """

        # now we recompute reduced state spaces for use in results processing
        # and we also convert the correction vectors to sparse representations
        tci = flat_total_space.indices
        N = len(tci)

        is_transp = len(corrs) == order and len(corrs[0]) == nstates
        is_sparse = isinstance(corrs[0], SparseArray) if is_transp else isinstance(corrs, SparseArray)
        # nstates = len(all_corrs)

        corr_mats = [None] * (order + 1)
        corr_inds = [[] for i in range(nstates)]
        for o in range(order + 1):
            if filters is not None:
                # we use this to check that we're only
                # allowing transitions that our filters support
                nquanta_rules = [
                    v for k,v in filters.items() if sum(k) == o
                ]
            else:
                nquanta_rules = None
            non_zeros = []
            if is_sparse:
                if is_transp:
                    sp_vals, sp_inds = corrs[o].block_data
                    nonzi = np.where(np.abs(sp_vals) > non_zero_cutoff)[0]
                    sp_vals = sp_vals[nonzi,]
                    sp_inds = tuple(s[nonzi] for s in sp_inds)
                    corr_mats[o] = SparseArray.from_data(
                        (
                            sp_vals,
                            sp_inds
                        ),
                        shape=(nstates, N),
                        cache_block_data=False
                    )
                    ind_keys, ind_vals = nput.group_by(sp_inds[1], sp_inds[0])[0]
                    for i,v in zip(ind_keys, ind_vals):
                        corr_inds[i].append(v)
                else:
                    raise NotImplementedError("constructing final coupling matrix from (order, nstates, N) `SparseArray` not supported")
            else:
                initial_quanta = np.sum(states.excitations, axis=1)
                for i in range(nstates):
                    if is_transp:
                        nonzi = np.where(np.abs(corrs[o, i]) > non_zero_cutoff)[0]
                        vals = corrs[o, i][nonzi,]
                    else:
                        nonzi = np.where(np.abs(corrs[i, o]) > non_zero_cutoff)[0]
                        vals = corrs[i, o][nonzi,]

                    if len(nonzi) > 0:
                        # we attempt to filter out things that can't touch based on our filter rules
                        target_quanta = np.sum(flat_total_space.take_subspace(nonzi).excitations, axis=1)
                        if nquanta_rules is not None and len(nquanta_rules) > 0:
                            # from .Solver import PerturbationTheoryStateSpaceFilter
                            mask = None
                            for f in nquanta_rules:
                                # f:PerturbationTheoryStateSpaceFilter
                                for (filter_space, filter_rules) in f.prefilters:
                                    is_in = states.take_subspace([i]).intersection(filter_space)
                                    if len(is_in) > 0:
                                        q_diffs = target_quanta - initial_quanta[i]
                                        poss_diffs = np.unique([sum(x) for x in filter_rules])
                                        if mask is None:
                                            mask = np.isin(q_diffs, poss_diffs)
                                        else:
                                            mask = np.logical_or(mask, np.isin(q_diffs, poss_diffs))
                            if mask is not None:
                                nonzi = nonzi[mask]
                                vals = vals[mask]
                    else:
                        if logger is not None:
                            logger.log_print("No corrections for state {s} at order {o}",
                                s=states.excitations[i],
                                o=o
                            )

                    # and then we add the appropriate basis indices to the list of basis data
                    non_zeros.append(
                        (
                            vals,
                            np.column_stack([
                                np.full(len(nonzi), i),
                                nonzi
                            ])
                        )
                    )

                    corr_inds[i].append(tci[nonzi,])

                # now we build the full mat rep for this level of correction
                vals = np.concatenate([x[0] for x in non_zeros])
                inds = np.concatenate([x[1] for x in non_zeros], axis=0).T
                corr_mats[o] = SparseArray.from_data(
                    (
                        vals,
                        inds
                    ),
                    shape=(nstates, N),
                    cache_block_data=False
                )

        # now we build state reps from corr_inds
        for i, dat in enumerate(corr_inds): #TODO: this might break with pruning...I can't really be sure at this point
            spaces = []
            for substates in dat:
                _, upos = np.unique(substates, return_index=True)
                usubs = substates[np.sort(upos)]
                spaces.append(flat_total_space.take_states(usubs))
            corr_inds[i] = BasisMultiStateSpace(np.array(spaces, dtype=object))

        return corr_mats, corr_inds

    def prune(self, threshold=.1, in_place=False):
        """
        Returns corrections with couplings less than the given cutoff set to zero

        :param threshold:
        :type threshold:
        :return:
        :rtype:
        """
        if not in_place:
            import copy
            new = copy.copy(self)
            new.prune(threshold=threshold, in_place=False)
            # might need to work harder here...
            return new

        for o in range(1, self.order):
            sp_vals, sp_inds = self.wfn_corrections[o].block_data
            mask = abs(sp_vals) > threshold
            sp_vals = sp_vals[mask]
            sp_inds = tuple(s[mask] for s in sp_inds)
            self.wfn_corrections[o].block_vals = sp_vals
            self.wfn_corrections[o].block_inds = sp_inds
        return self

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
    def get_transformed_Hamiltonians(self, hams, deg_group=None):
        """

        :param corrs:
        :type corrs:
        :param deg_group:
        :type deg_group:
        :return:
        :rtype:
        """
        if deg_group is not None:
            subcorrs = self.take_subspace(deg_group)
            inds = self.total_basis.find(deg_group)
            subhams = [SparseArray.from_data(self._take_subham(H, inds)) for H in hams]
            # H_nd =[
            #     x.asarray() if isinstance(x, SparseArray) else x
            #     for x in subcorrs.operator_representation(subhams, subspace=deg_group)
            # ]
            H_nd = [
                x.asarray() if isinstance(x, SparseArray) else x
                for x in subcorrs.operator_representation(subhams, subspace=deg_group, logger_symbol="H")
            ]

        else:
            subhams = hams
            H_nd = [
                x.asarray() if isinstance(x, SparseArray) else x
                for x in self.operator_representation(subhams, logger_symbol="H", logger_conversion=UnitsData.convert("Hartrees", "Wavenumbers"))
            ]
        return H_nd
    def get_degenerate_rotation(self, deg_group, hams):
        """

        :param deg_group:
        :type deg_group:
        :param corrs:
        :type corrs:
        :return:
        :rtype:
        """

        logger = self.logger
        # raise Exception(corrs.states.excitations, deg_group)
        with logger.block(tag="states"):
            logger.log_print(
                str(
                    self.states.take_states(deg_group).excitations
                ).splitlines()
            )
        subdegs = self.take_subspace(deg_group)

        # from McUtils.Scaffolding import JSONSerializer
        # import os
        # with open(os.path.expanduser("~/Desktop/wat6.json"), "w+") as woof:
        #     JSONSerializer().serialize(woof, subdegs)

        # H_nd = self.get_transformed_Hamiltonians(corrs, deg_group)
        # for h in H_nd[1:]:
        #     np.fill_diagonal(h, 0.)
        H_nd_corrs = subdegs.get_transformed_Hamiltonians(hams, deg_group=None)
        # import McUtils.Plots as plt
        # plt.TensorPlot(np.array(H_nd)).show()
        H_nd = np.sum(H_nd_corrs, axis=0)
        if np.sum(H_nd) == 0:
            raise Exception(subdegs.wfn_corrections)
        #     raise Exception(deg_group.excitations,
        #                     self.states.take_states(deg_group).excitations,
        #                     # self.coupled_states.take_states(deg_group).excitations
        #                     )
        # overlaps = np.sum(subdegs.get_overlap_matrices(), axis=0)

        with logger.block(tag="non-degenerate Hamiltonian"):
            logger.log_print(
                str(
                    np.round(H_nd * UnitsData.convert("Hartrees", "Wavenumbers")).astype(int)
                ).splitlines()
            )

        deg_engs, deg_transf = np.linalg.eigh(H_nd)

        ov_thresh = .5
        for i in range(len(deg_transf)):
            max_ov = np.max(deg_transf[:, i] ** 2)
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
        # # if len(sorting) != len(np.unique(sorting)):
        # #     raise PerturbationTheoryException("After diagonalizing can't distinguish modes...")
        deg_engs = deg_engs[sorting,]

        self.logger.log_print("degenerate energies {e}",
                              e=np.round(deg_engs * UnitsData.convert("Hartrees", "Wavenumbers")))

        deg_transf = deg_transf[:, sorting]

        return H_nd_corrs, deg_engs, deg_transf

    def get_degenerate_transformation(self, group, hams, gaussian_resonance_handling=False):
        # this will be built from a series of block-diagonal matrices
        # so we store the relevant values and indices to compose the SparseArray

        # we apply the degenerate PT on a group-by-group basis
        # by transforming the H reps into the non-degenerate basis
        # print(">", group.excitations)
        deg_inds = self.states.find(group, missing_val=-1)
        mask = deg_inds > -1
        deg_inds = deg_inds[mask]
        if not mask.all():
            bad = group.take_subspace(np.where(np.logical_not(mask))[0])
            self.logger.log_print(
                "WARNING: got degeneracy spec including states {} that are not in space of corrected states".format(
                    bad.excitations
                ))
        group = group.take_subspace(np.where(mask)[0])
        # print(">..", deg_inds, self.states.find(group, missing_val=-1))

        if len(deg_inds) == 1 or (
                gaussian_resonance_handling and np.max(np.sum(group.excitations, axis=1)) > 2):
            H_nd = deg_engs = deg_rot = None
        elif len(deg_inds) > 1:
            H_nd, deg_engs, deg_rot = self.get_degenerate_rotation(group, hams)
        else:
            H_nd = deg_engs = deg_rot = None
            # raise NotImplementedError("Not sure what to do when no states in degeneracy spec are in total space")

        return deg_inds, H_nd, deg_rot, deg_engs

    @staticmethod
    def default_state_filter(state, couplings, energy_cutoff=None, energies=None, basis=None, target_modes=None):
        """
        Excludes modes that differ in only one position, prioritizing states with fewer numbers of quanta
        (potentially add restrictions to high frequency modes...?)

        :param input_state:
        :type input_state:
        :param couplings:
        :type couplings:
        :return:
        :rtype:
        """

        if target_modes is None:
            target_modes = np.arange(len(state.excitations[0]))

        if energy_cutoff is not None:
            state_ind = basis.find(state)
            coupling_inds = basis.find(couplings)
            diff_mask = np.abs(energies[coupling_inds] - energies[state_ind]) > energy_cutoff
        else:
            exc_1 = state.excitations[0, target_modes]
            exc_2 = couplings.excitations[:, target_modes]
            diffs = exc_2 - exc_1[np.newaxis, :]
            diff_sums = np.sum(diffs != 0, axis=1)
            diff_mask = np.logical_or(
                np.sum(couplings.excitations, axis=1) == 0, # drop ground state
                diff_sums == 1 # find where changes are only in one position
            )
        # now drop these modes
        # print(couplings.excitations)
        if diff_mask.any():
            if diff_mask.all():
                couplings = None
            else:
                # print(">>>", state.excitations[0])
                # print(couplings.excitations)
                couplings = couplings.take_subspace(np.where(np.logical_not(diff_mask))[0])
                # print("===")
                # print(couplings.excitations)
        # print(couplings.excitations)
        return couplings

    def find_strong_couplings(self, threshold=.1, state_filter=None):
        """
        Finds positions in the expansion matrices where the couplings are too large

        :param threshold:
        :type threshold:
        :return:
        :rtype:
        """

        if state_filter is None:
            state_filter = self.default_state_filter

        order = self.order
        strong_couplings = {}
        for o in range(1, order):
            sp_vals, sp_inds = self.wfn_corrections[o].block_data
            nonzi = np.where(np.abs(sp_vals) > threshold)[0]
            sp_inds = tuple(s[nonzi] for s in sp_inds)
            ind_keys, ind_vals = nput.group_by(sp_inds[1], self.states.indices[sp_inds[0]])[0]
            # print(self.states.indices[sp_inds[0]])
            # print(tci[sp_inds[1]])
            for i, v in zip(ind_keys, ind_vals):
                _, upos = np.unique(v, return_index=True)
                usubs = v[np.sort(upos)]
                states = state_filter(self.states.take_states([i]), self.total_basis.take_subspace(usubs))
                if states is not None and len(states) > 0:
                    if i not in strong_couplings:
                        strong_couplings[i] = [None for _ in range(order)]
                    strong_couplings[i][o] = states

        return strong_couplings

    def format_strong_couplings_report(self, couplings=None, threshold=.1, int_fmt="{:>3.0f}", padding="{:<8}", join=True, use_excitations=True):
        if couplings is None:
            couplings = self.find_strong_couplings(threshold=threshold)

        list_fmt = " ".join(int_fmt for _ in range(self.total_basis.ndim))
        coupling_statements = []
        for i,v in sorted(couplings.items(), key=lambda k:k[0]):
            if use_excitations:
                i = self.states.take_states([i]).excitations[0]
                coupling_statements.append(padding.format("state:") + list_fmt.format(*i))
                # print(v)
                for n,l in enumerate(v):
                    if l is not None and len(l) > 0:
                        coupling_statements.extend(padding.format(" order {}".format(n) if j == 0 else "")+list_fmt.format(*e) for j,e in enumerate(l.excitations))
            else:
                coupling_statements.append(list_fmt.format(i))
                coupling_statements.append(padding + str(v.indices))

        return coupling_statements if not join else "\n".join(coupling_statements)

    def collapse_strong_couplings(self, sc:dict):
        """

        :param sc:
        :type sc:
        :return:
        :rtype:
        """
        new = {}
        for k,v in sc.items():
            s = None
            for s2 in v:
                if s2 is not None:
                    if s is None:
                        s = s2
                    else:
                        s = s.concatenate(s2)
            new[k] = s
        return new

    # def apply_martin_test(self, hams):
    #     zoos =
    #     (np.abs(H1_block) ** 4) / (diffs ** 3)

    @staticmethod
    def _fmt_operator_rep(full_ops, operator_symbol, conversion, real_fmt="{:>.8e}", padding_fmt='{:>16}'):
        tag_line = None
        rep_lines = None

        op_dim = None

        for (a, b, c), subrep in full_ops:
            if isinstance(subrep, SparseArray):
                subrep = subrep.asarray()
            elif isinstance(subrep, (int, float, np.integer, np.floating)):
                if subrep == 0:
                    if op_dim is None:
                        raise ValueError("was lazy and haven't filled operator dim yet...")
                    subrep = np.zeros(op_dim)
                else:
                    raise ValueError("don't know what to do with representation '{}'".format(subrep))

            if op_dim is None:
                op_dim = subrep.shape

            if conversion is not None:
                subrep = subrep * conversion

            subrep_lines = [
                " ".join(padding_fmt.format(real_fmt.format(e)) for e in line)
                for line in subrep
            ]
            line_len = len(subrep_lines[0])

            if rep_lines is None:
                rep_lines = subrep_lines
            else:
                rep_lines = [x + " " + y for x,y in zip(rep_lines, subrep_lines)]

            tag_fmt = "{:<" + str(line_len) + "}"
            base_tag=tag_fmt.format("<{a}|{A}({c})|{b}>".format(A=operator_symbol, a=a, c=c, b=b))
            if tag_line is None:
                tag_line = base_tag
            else:
                tag_line += " " + base_tag

        # we want to return a line list so the logger can add any prefixes it needs
        rep_lines.insert(0, tag_line)
        return rep_lines

    def operator_representation(self, operator_expansion, order=None, subspace=None, contract=True,
                                logger_symbol="A",
                                logger_conversion=None
                                ):
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
            for k in range(order):
                wfn_corrs.append(self.wfn_corrections[k][:, subspace_sel])

        # generalizes the dot product so that we can use 0 as a special value...
        dot = _safe_dot
        logger = self.logger
        logger = None if logger is None or isinstance(logger, NullLogger) else logger

        # does the dirty work of acutally applying the rep...
        reps = [[] for _ in range(order)]
        full_ops = []
        for k in range(order):
            tags = []
            op = []
            # apply each thing up to requested order...
            for a in range(k+1): # if k == 2: a=0, a=1, a=2
                for b in range(k-a+1): # if k==2, a==0: b=0, b=1, b=2; a==1: b=0, b=1
                    c = k - (a + b) # a + b + c == k
                    rop = operator_expansion[c]
                    if isinstance(rop, (int, float, np.integer, np.floating)): # constant reps...
                        if rop != 0: # cheap easy check
                            subrep = rop * dot(wfn_corrs[a], wfn_corrs[b].T)
                            op.append(subrep)
                        else:
                            subrep = 0
                            op.append(0)
                    else:
                        subrep = dot(dot(wfn_corrs[a], rop), wfn_corrs[b].T)
                        op.append(subrep)

                    full_ops.append([
                        (a, b, c),
                        subrep
                    ])

            if contract:
                op = sum(op)
            reps[k] = op

        if logger is not None:
            logger.log_print(full_ops, logger_symbol, logger_conversion, message_prepper=self._fmt_operator_rep)

        return reps

    def get_overlap_matrices(self):
        """
        Returns the overlap matrices for the set of corrections
        at each order of correction

        :return:
        :rtype:
        """

        wat = []
        for k in range(2 + 1):
            ov = None
            for i in range(k + 1):
                c1 = self.wfn_corrections[i].asarray()
                c2 = self.wfn_corrections[k - i].asarray()
                if ov is None:
                    ov = np.dot(c1, c2.T)
                else:
                    ov += np.dot(c1, c2.T)
            wat.append(ov)

        return wat

    # def checkpoint_save(self, checkpoint:Checkpointer):
    #     """
    #     Writes correction arrays to checkpoint
    #
    #     :param checkpoint:
    #     :type checkpoint:
    #     :return:
    #     :rtype:
    #     """
    #     import gc, sys
    #
    #     # self.states = states
    #     # self.coupled_states = coupled_states
    #     # self.total_basis = total_basis
    #     # self.energy_corrs = energy_corrs
    #     # self.all_energy_corrs = all_energy_corrections
    #     # self.wfn_corrections = wfn_corrections
    #     # self.degenerate_states = degenerate_states
    #     # self.degenerate_transf = degenerate_transformation
    #     # self.degenerate_energies = degenerate_energies
    #     # self.logger = logger
    #     checkpoint["wfn_corrections"] = self.wfn_corrections
    #     print(">>>>", sys.getrefcount(self.wfn_corrections))
    #     self.wfn_corrections = None
    #     gc.collect()
    #
    # def checkpoint_reload(self, checkpoint:Checkpointer):
    #     self.wfn_corrections = checkpoint['wfn_corrections']
    #
    # def disk_backed(self):
    #     return self.chk_backer(self)
    #
    # class chk_backer:
    #     def __init__(self, parent):
    #         self.parent = parent
    #         self.chk = None
    #     def load_chk(self):
    #         import tempfile as tf
    #         target = tf.NamedTemporaryFile(suffix=".hdf5").name
    #         self.chk = Checkpointer.from_file(target)
    #     def unload_chk(self):
    #         import os
    #         os.remove(self.chk.checkpoint_file)
    #     def __enter__(self):
    #         self.load_chk()
    #         self.chk.__enter__()
    #         self.parent.checkpoint_save(self.chk)
    #     def __exit__(self, exc_type, exc_val, exc_tb):
    #         self.parent.checkpoint_reload(self.chk)
    #         self.chk.__exit__(exc_type, exc_val, exc_tb)
    #         self.unload_chk()

    def savez(self, file):
        raise NotImplementedError("old and wrong now")
        keys = dict(
            states=self.states,
            coupled_states=self.coupled_states,
            total_states=self.total_basis,
            energies=self.energy_corrs,
            wavefunctions=self.wfn_corrections
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
        return cls.from_dicts(
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

    def to_state(self, serializer=None):
        keys = dict(
            states=self.states,
            coupled_states=self.coupled_states,
            total_states=self.total_basis,
            energies=self.energy_corrs,
            wavefunctions=self.wfn_corrections,
            degenerate_states=self.degenerate_states,
            degenerate_transformations=self.degenerate_transf,
            degenerate_energies=self.degenerate_energies
        )
        return keys
    @classmethod
    def from_state(cls, data, serializer=None):
        return cls.from_dicts(
            {
                "states": serializer.deserialize(data['states']),
                "coupled_states": serializer.deserialize(data['coupled_states']),
                "total_states": serializer.deserialize(data['coupled_states']),
                "degenerate_states": serializer.deserialize(data['degenerate_states']),
            },
            {
                "energies": data['energies'],
                "wavefunctions": data['wavefunctions'],
                "degenerate_transformation": data['degenerate_transformations'],
                "degenerate_energies": data['degenerate_energies']
            },
            data['hamiltonians'] # we probably want to ditch this for memory reasons...
        )