
import numpy as np, itertools, collections, dataclasses

from McUtils.Numputils import SparseArray
import McUtils.Numputils as nput
from McUtils.Data import UnitsData
from McUtils.Scaffolding import Logger, NullLogger, Checkpointer

from ..Spectra import DiscreteSpectrum
from ..BasisReps import (
    BasisStateSpace, BasisMultiStateSpace, SelectionRuleStateSpace, HarmonicOscillatorProductBasis,
    StateMaker
)
from .Common import PerturbationTheoryException, _safe_dot

__all__ = [
    "PerturbationTheoryCorrections",
    "AnalyticPerturbationTheoryCorrections",
    "BasicAPTCorrections"
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
                 nondeg_hamiltonian_precision=3,
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
        self.nondeg_hamiltonian_precision = nondeg_hamiltonian_precision

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
    def get_degenerate_rotation(self, deg_group, hams, label=None, zero_point_energy=None,
                                local_coupling_hamiltonian=None,
                                local_coupling_order=None):
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
        with logger.block(tag="states" if label is None else label):
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
        if local_coupling_order is None:
            local_coupling_order = -1
        if local_coupling_order < 0:
            local_coupling_order = len(hams) + local_coupling_order
        if local_coupling_order < len(hams) - 1:
            #TODO: calculate only the required terms at this order
            raise NotImplementedError("local mode couplings only provided at highest order")
        if local_coupling_hamiltonian is not None:
            hc = local_coupling_hamiltonian.get_representation_matrix(deg_group, deg_group,
                                                                      zero_element_warning=False,
                                                                      diagonal=False,
                                                                      # memory_constrained=self.memory_constrained
                                                                      )
        else:
            hc = None

        H_nd_corrs = subdegs.get_transformed_Hamiltonians(hams, deg_group=None)
        group_inds = self.states.find(deg_group)
        if hc is not None:
            H_nd_corrs[local_coupling_order] += hc.asarray()
        # zero_order_engs = self.energy_corrs[group_inds, 0]
        # raise Exception(
        #     (H_nd_corrs[0]-np.diag(np.full(len(group_inds), self.energy_corrs[0, 0])))*UnitsData.convert("Hartrees", "Wavenumbers"),
        #     (zero_order_engs-self.energy_corrs[0, 0])*UnitsData.convert("Hartrees", "Wavenumbers")
        # )
        # np.fill_diagonal(H_nd_corrs[0], zero_order_engs)

        # import McUtils.Plots as plt
        # plt.TensorPlot(np.array(H_nd)).show()
        H_nd = np.sum(H_nd_corrs, axis=0)
        if np.sum(H_nd) == 0:
            raise ValueError("No corrections from ", subdegs.wfn_corrections)
        #     raise Exception(deg_group.excitations,
        #                     self.states.take_states(deg_group).excitations,
        #                     # self.coupled_states.take_states(deg_group).excitations
        #                     )
        # overlaps = np.sum(subdegs.get_overlap_matrices(), axis=0)

        with logger.block(tag="non-degenerate Hamiltonian"):
            if zero_point_energy is None:
                zero_point_energy = self.energies[0]
            with np.printoptions(precision=self.nondeg_hamiltonian_precision, suppress=True):
                logger.log_print(
                    str(
                        (H_nd - np.diag(np.full(len(group_inds), zero_point_energy))) * UnitsData.convert("Hartrees", "Wavenumbers")
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

    def get_degenerate_transformation(self, group, hams, gaussian_resonance_handling=False, label=None,
                                      zero_point_energy=None,
                                      local_coupling_hamiltonian=None,
                                      local_coupling_order=None
                                      ):
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
            H_nd, deg_engs, deg_rot = self.get_degenerate_rotation(group, hams, label=label,
                                                                   zero_point_energy=zero_point_energy,
                                                                   local_coupling_hamiltonian=local_coupling_hamiltonian,
                                                                   local_coupling_order=local_coupling_order
                                                                   )
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

    def operator_representation(self, operator_expansion,
                                order=None, subspace=None, contract=True,
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
        raise NotImplementedError("old and wrong now")
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


BasicAPTCorrections = collections.namedtuple("PTCorrections",
                                       ['initial_states', 'final_states', 'corrections']
                                       )

@dataclasses.dataclass
class AnalyticPerturbationTheoryCorrections:
    states: BasisStateSpace
    state_lists: 'list[tuple[np.ndarray, np.ndarray]]'
    _energies: np.ndarray = None
    _transition_moments: 'Iterable[np.ndarray]' = None
    _spectra: 'Iterable[DiscreteSpectrum]' = None
    _deperturbed_energies: np.ndarray = None
    _deperturbed_transition_moments: 'Iterable[np.ndarray]' = None
    _deperturbed_spectra: DiscreteSpectrum = None
    degenerate_states: 'Iterable[BasisStateSpace]' = None
    only_degenerate_terms: 'bool' = True
    _degenerate_hamiltonians: 'Iterable[np.ndarray]' = None
    _degenerate_coefficients: 'Iterable[np.ndarray]' = None
    _degenerate_state_list_transformations: 'Iterable[list[np.ndarray, np.ndarray]]' = None
    energy_corrections: BasicAPTCorrections = None
    transition_moment_corrections: 'Iterable[BasicAPTCorrections]' = None
    degenerate_hamiltonian_corrections: 'Iterable[BasicAPTCorrections]' = None
    operator_corrections: 'Iterable[BasicAPTCorrections]' = None
    _deperturbed_operator_values: 'Iterable[np.ndarray]' = None
    _operator_values: 'Iterable[np.ndarray]' = None
    operator_keys: 'Iterable[Any]' = None
    logger: 'Logger' = None

    _zpe_pos: int = None

    def get_zpe_pos(self) -> int:
        if self._zpe_pos is None:
            basis = self.states
            try:
                zpe_pos = basis.find([np.zeros(basis.ndim, dtype=int)])
            except IndexError:
                zpe_pos = 0  # do this better in the future...
            else:
                zpe_pos = zpe_pos[0]
            self._zpe_pos = zpe_pos
        return self._zpe_pos

    @property
    def energies(self) -> np.ndarray:
        if self._energies is None:
            if self.degenerate_states is None:
                self._energies = np.sum(self.energy_corrections, axis=0)
            else:
                self._deperturbed_energies = np.sum(self.energy_corrections, axis=0)
                self._energies, (self._degenerate_hamiltonians, self._degenerate_coefficients) = self.get_degenerate_transformations(
                    self.states, self._deperturbed_energies
                )
        return self._energies

    @property
    def deperturbed_energies(self) -> np.ndarray:
        if self._deperturbed_energies is None:
            self._deperturbed_energies = np.sum(self.energy_corrections, axis=0)
        return self._deperturbed_energies

    @classmethod
    def handle_degenerate_transformation(cls, degenerate_ham):
        deg_engs, deg_transf = np.linalg.eigh(degenerate_ham)

        # ov_thresh = .5
        # for i in range(len(deg_transf)):
        #     max_ov = np.max(deg_transf[:, i] ** 2)
        #     if max_ov < ov_thresh:  # there must be a single mode that has more than 50% of the initial state character?
        #         if logger is not None:
        #             logger.log_print(
        #                 "! state {i} is more than 50% mixed",
        #                 i=i
        #             )

        # we pick the terms with the max contribution from each input state
        # and zero out the contributions so that two states can't map
        # to the same input state
        sort_transf = np.abs(deg_transf.copy())
        sorting = [-1] * len(deg_transf)
        for i in range(len(deg_transf)):
            o = np.argmax(sort_transf[i, :])
            sorting[i] = o
            sort_transf[:, o] = 0.  # np.zeros(len(sort_transf))
        # print(">>>>", sorting)

        # self.logger.log_print('sorting: {s}', s=sorting)

        deg_engs = deg_engs[sorting,]
        deg_transf = deg_transf[:, sorting]
        return deg_engs, deg_transf

    def get_degenerate_transformations(self, basis, energies):
        energies = energies.copy()

        hams = []
        transf = []
        for block_num, (group, hcorr) in enumerate(
                zip(self.degenerate_states, self.degenerate_hamiltonian_corrections)
        ):
            ham = np.sum(hcorr, axis=0)
            # with self.logger.block(tag="Degenerate Block {n}", n=block_num + 1):
            #     self.logger.log_print("{s}",
            #                           s=group,
            #                           preformatter=lambda **kw: dict(
            #                               kw,
            #                               s=self.format_matrix(kw['s'])
            #                           ))
            e_pos = basis.find(group)
            if not self.only_degenerate_terms:
                ham = -ham  # only deg terms
            engs = energies[e_pos,]
            ham[np.diag_indices_from(ham)] = engs
            deg_engs, deg_mixing = self.handle_degenerate_transformation(ham)

            # with self.logger.block(tag='contributions'):
            #     self.logger.log_print(
            #         "{contribs}",
            #         contribs=deg_mixing,
            #         preformatter=lambda **kw: dict(
            #             kw,
            #             contribs=self.format_matrix(np.round(100 * (deg_mixing ** 2)))
            #         )
            #     )

            energies[e_pos] = deg_engs
            transf.append(deg_mixing)
            hams.append(ham)

        return energies, (hams, transf)

    @property
    def degenerate_hamiltonians(self):
        if self._degenerate_hamiltonians is None:
            e_base = self.energies # TODO: this is bad practice...
        return self._degenerate_hamiltonians
    @property
    def degenerate_coefficients(self):
        if self._degenerate_coefficients is None:
            e_base = self.energies # TODO: this is bad practice...
        return self._degenerate_coefficients

    def get_freqs(self):
        return self.energies - self.energies[self.get_zpe_pos()]

    def get_deperturbed_freqs(self):
        if self.degenerate_states is not None:
            return self.deperturbed_energies - self.deperturbed_energies[self.get_zpe_pos()]
        else:
            return self.get_freqs()

    @property
    def degenerate_transformation_pairs(self):
        if self._degenerate_state_list_transformations is None:
            self._degenerate_state_list_transformations = self._get_degenerate_tfs_mats()
        return self._degenerate_state_list_transformations
    def _get_degenerate_tfs_mats(self, logger=None):
        #TODO: add checks to ensure that our blocks are complete and we have proper unitary tfs at the end
        if logger is None:
            logger = self.logger
        logger = Logger.lookup(logger)
        all_degs = BasisStateSpace(
            HarmonicOscillatorProductBasis(len(self.state_lists[0][0][0])),
            np.concatenate(self.degenerate_states, axis=0)
        )
        deg_map_row = np.concatenate([
            [i] * len(g) for i, g in enumerate(self.degenerate_states)
        ])
        deg_map_col = np.concatenate([
            np.arange(len(g)) for i, g in enumerate(self.degenerate_states)
        ])
        for i, block in enumerate(self.degenerate_states):
            with logger.block(tag="Degenerate block {i}", i=i + 1):
                logger.log_print(
                    "{blocks}",
                    blocks=block,
                    preformatter=lambda **kw: dict(
                        kw,
                        blocks="\n".join(StateMaker.parse_state(e) for e in kw['blocks'])
                    )
                )

        tfs = []
        for block_idx, (init_states, final_states) in enumerate(self.state_lists):
            initial_space = BasisStateSpace(
                HarmonicOscillatorProductBasis(len(init_states[0])),
                init_states
            )
            final_space = BasisStateSpace(
                HarmonicOscillatorProductBasis(len(init_states[0])),
                final_states
            )
            init_pos = all_degs.find(initial_space, missing_val=-1)
            final_pos = all_degs.find(final_space, missing_val=-1)

            # subcorr is a n_init x n_final object, but we need to figure out the transformation to apply to each axis
            # to do so we find the appropriate transformation and insert it
            row_tf = np.eye(len(init_states))
            for n,i in enumerate(init_pos):
                if i != -1:
                    col_pos = deg_map_row[i]
                    row_pos = deg_map_col[i]
                    deg_block = self.degenerate_coefficients[col_pos]
                    block_idx = initial_space.find(self.degenerate_states[col_pos])
                    # print(row_pos, col_pos, col_tf.shape)
                    row_tf[block_idx, n] = deg_block[:, row_pos]

            nz_init_pos = [i for i in init_pos if i > -1]
            nz_init = [s for i, s in zip(init_pos, init_states) if i > -1]
            if len(nz_init_pos) > 0:
                with logger.block(tag="Initial states ({ix})",
                                  ix=[(deg_map_row[i],  deg_map_col[i]) for i in nz_init_pos]):
                    logger.log_print(
                        "{initials}",
                        initials=nz_init,
                        preformatter=lambda **kw: dict(
                            kw,
                            initials="\n".join(StateMaker.parse_state(e) for e in kw['initials'])
                        )
                    )
                    logger.log_print(
                        "{tf}",
                        tf=row_tf,
                        preformatter=lambda **kw:dict(kw, tf="\n".join(logger.prep_array(kw['tf'])))
                    )
            if np.sum(np.abs((row_tf @ row_tf.T) - np.eye(len(init_states))).flatten()) > 1e-3:
                raise ValueError("Non-unitary row tf, something wrong with initial state degs")

            col_tf = np.eye(len(final_states))
            for n,i in enumerate(final_pos):
                if i != -1:
                    col_pos = deg_map_row[i]
                    deg_block = self.degenerate_coefficients[col_pos]
                    block_idx = final_space.find(self.degenerate_states[col_pos])
                    row_pos = deg_map_col[i]
                    col_tf[block_idx, n] = deg_block[:, row_pos]
            nz_final_pos = [i for i in final_pos if i > -1]
            nz_final = [s for i, s in zip(final_pos, final_states) if i > -1]
            if len(nz_final_pos) > 0:
                with logger.block(tag="Final states ({ix})",
                                  ix=[(deg_map_row[i],  deg_map_col[i]) for i in nz_final_pos]):
                    logger.log_print(
                        "{finals}",
                        finals=nz_final,
                        preformatter=lambda **kw: dict(
                            kw,
                            finals="\n".join(StateMaker.parse_state(e) for e in kw['finals'])
                        )
                    )
                    logger.log_print(
                        "{tf}",
                        tf=col_tf,
                        preformatter=lambda **kw:dict(kw, tf="\n".join(logger.prep_array(kw['tf'])))
                    )
            if np.sum(np.abs((row_tf @ row_tf.T) - np.eye(len(init_states))).flatten()) > 1e-3:
                raise ValueError("Non-unitary col tf, something wrong with final state degs")

            tfs.append([row_tf, col_tf])
        return tfs
    def _apply_degs_to_corrs(self, corrs, logger=None):
        if logger is None:
            logger = self.logger
        logger = Logger.lookup(logger)


        has_subcorrs = not (isinstance(corrs[0], np.ndarray) and corrs[0].ndim == 2)
        if has_subcorrs:
            all_tms = [[] for _ in corrs]  # num axes
        else:
            all_tms = []
        tfs = self.degenerate_transformation_pairs

        for block_idx, (row_tf, col_tf) in enumerate(tfs):

            if has_subcorrs:
                for a, (storage, axis_moms) in enumerate(zip(all_tms, corrs)):
                    tf_corr = row_tf.T @ axis_moms[block_idx] @ col_tf
                    storage.append(tf_corr)
            else:
                all_tms.append(row_tf.T @ corrs[block_idx] @ col_tf)

        return all_tms

    @property
    def transition_moments(self):
        if self._transition_moments is None:
            logger = Logger.lookup(self.logger)
            null = NullLogger()
            if self.degenerate_states is None:
                self._transition_moments = self.deperturbed_transition_moments
            else:
                tmoms = self.deperturbed_transition_moments

                self._transition_moments = self._apply_degs_to_corrs(
                    tmoms,
                    logger = self.logger
                )

        return self._transition_moments

    @property
    def harmonic_transition_moments(self):
        return [
            [corr_block[0] for corr_block in axis]
            for axis in self.transition_moment_corrections
            # np.sum([np.sum(corr.corrections[n], axis=0) ** 2 for corr in self.transition_moment_corrections], axis=0)
            # for n in range(len(self.transition_moment_corrections[0]))
        ]

    @property
    def deperturbed_transition_moments(self):
        if self._deperturbed_transition_moments is None:
            self._deperturbed_transition_moments = [
                [np.sum(corr_block, axis=0) for corr_block in axis]
                for axis in self.transition_moment_corrections
                # np.sum([np.sum(corr.corrections[n], axis=0) ** 2 for corr in self.transition_moment_corrections], axis=0)
                # for n in range(len(self.transition_moment_corrections[0]))
            ]
        return self._deperturbed_transition_moments

    def get_spectra(self, energies, transition_moments):
        spectra = []
        for block_idx, (init_states, final_states) in enumerate(self.state_lists):
            block_specs = []
            init_engs = energies[self.states.find(init_states),]
            final_engs = energies[self.states.find(final_states),]
            for i, init in enumerate(init_states):
                submoms = np.array([
                    axis_moms[block_idx][i]
                    for axis_moms in transition_moments
                ]).T
                freqs = final_engs - init_engs[i]
                block_specs.append(DiscreteSpectrum.from_transition_moments(freqs, submoms))
            spectra.append(block_specs)
        return spectra

    @property
    def harmonic_spectra(self):
        return self.get_spectra(
            self.energy_corrections[0],
            self.harmonic_transition_moments
        )

    @property
    def deperturbed_spectra(self):
        if self._deperturbed_spectra is None:
            self._deperturbed_spectra = self.get_spectra(
                self.deperturbed_energies,
                self.deperturbed_transition_moments
            )
        return self._deperturbed_spectra

    @property
    def spectra(self):
        if self._spectra is None:
            if self.degenerate_states is None:
                self._spectra = self.deperturbed_spectra
            else:
                self._spectra = self.get_spectra(
                    self.energies,
                    self.transition_moments
                )
        return self._spectra

    @property
    def deperturbed_operator_values(self):
        if self._deperturbed_operator_values is None:
            self._deperturbed_operator_values = [
                [np.sum(corr_block, axis=0) for corr_block in op]
                for op in self.operator_corrections
            ]
        return self._deperturbed_operator_values

    @property
    def operator_values(self):
        if self._operator_values is None:
            if self.degenerate_states is None:
                self._operator_values = self.deperturbed_operator_values
            else:
                tmoms = self.deperturbed_operator_values

                self._operator_values = self._apply_degs_to_corrs(
                    tmoms,
                    logger=self.logger
                )

                # ops = self.deperturbed_operator_values
                # all_tms = [[] for _ in ops]  # num operators
                # for block_idx, (init_states, final_states) in enumerate(self.state_lists):
                #     for storage, axis_moms in zip(all_tms, ops):
                #         storage.append(
                #             self.apply_degenerate_transformations(
                #                 init_states, final_states,
                #                 axis_moms[block_idx]
                #             )
                #         )
                # self._operator_values = all_tms

        return self._operator_values
