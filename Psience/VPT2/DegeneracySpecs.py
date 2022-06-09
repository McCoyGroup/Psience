

import numpy as np, enum, abc, scipy.sparse as sp
from McUtils.Combinatorics import SymmetricGroupGenerator, PermutationRelationGraph
import McUtils.Numputils as nput
from ..BasisReps import BasisStateSpace, BasisMultiStateSpace, SelectionRuleStateSpace, BraKetSpace, HarmonicOscillatorProductBasis

__all__ = [
    "DegeneracySpec",
    "DegenerateMultiStateSpace"
]

__reload_hook__ = ["..BasisReps"]

class DegenerateSpaceInputFormat(enum.Enum):
    Groups = "groups"
    QuantaSpecRules = "nT"
    Polyads = "polyads"
    StrongCouplings = "wfc_threshold"
    EnergyCutoff = "energy_window"
    Callable = "callable"
    MartinTest = "martin_threshold"

class DegeneracySpec(metaclass=abc.ABCMeta):
    """
    Provides a container for specifying degeneracies
    in a way that can be cleanly canonicalized
    """

    application_order = 'pre'
    group_filter = None
    def __init__(self, application_order=None, group_filter=None, energy_cutoff=2.25e-3, test_modes=None, maximize_filtered_groups=True):
        if application_order is None:
            application_order = self.application_order
        self.application_order = application_order
        self.energy_cutoff = energy_cutoff
        if group_filter is None:
            group_filter = self.group_filter
        self.group_filter = group_filter
        self.test_modes = test_modes
        self.maximize_filtered_groups = maximize_filtered_groups

    repr_opts = ['energy_cutoff']
    def __repr__(self):
        return "{}({})".format(type(self).__name__, ", ".join('{}={}'.format(k, getattr(self, k)) for k in self.repr_opts))

    def get_degenerate_group_filter(self, solver, corrs=None, threshold=None):
        group_filter = self.group_filter
        if isinstance(group_filter, str) and group_filter == 'default':
            group_filter = None
        if group_filter is None:
            zero_order_energy_cutoff = self.energy_cutoff
            test_modes = self.test_modes
            if test_modes is None:
                test_modes = solver.high_frequency_modes
            group_filter = dict(
                corrections=corrs,
                energy_cutoff=zero_order_energy_cutoff,
                energies=solver.zero_order_energies,
                threshold=threshold,
                target_modes=test_modes,
                maximize_groups=self.maximize_filtered_groups
            )
        elif isinstance(group_filter, str) and group_filter == 'unfiltered':
            group_filter = None
        return group_filter

    format = None
    @classmethod
    def get_format_mapping(cls):
        return {
            DegenerateSpaceInputFormat.Groups: GroupsDegeneracySpec,
            DegenerateSpaceInputFormat.QuantaSpecRules: TotalQuantaDegeneracySpec,
            DegenerateSpaceInputFormat.Polyads: PolyadDegeneracySpec,
            DegenerateSpaceInputFormat.StrongCouplings: StronglyCoupledDegeneracySpec,
            DegenerateSpaceInputFormat.EnergyCutoff: EnergyCutoffDegeneracySpec,
            DegenerateSpaceInputFormat.Callable: CallableDegeneracySpec,
            DegenerateSpaceInputFormat.MartinTest: MartinTestDegeneracySpec,
        }

    @classmethod
    def get_default_spec(cls):
        return StronglyCoupledDegeneracySpec()
    @classmethod
    def from_spec(cls, spec, format=None, **kwargs):
        if spec is None:
            return None
        elif (
                isinstance(spec, str) and spec == 'auto'
                or spec is True
        ):
            return cls.get_default_spec()

        if format is None:
            (spec, opts), format = cls.infer_format(spec)
            if opts is not None:
                kwargs = dict(opts, **kwargs)
        res_cls = cls.get_format_mapping()[format]

        return res_cls(spec, **kwargs)

    @staticmethod
    def _key(obj, key):

        try:
            val = getattr(obj, key)
        except AttributeError:
            val = None
        if val is None:
            try:
                val = obj[key]
            except (TypeError, KeyError, AttributeError, IndexError):
                val = None
        return val


    @classmethod
    def infer_format(cls, spec):

        if spec is None:
            return (None, None), None

        fmt = cls._key(spec, 'format')
        if fmt is not None:
            return spec, fmt

        for fmt in DegenerateSpaceInputFormat:
            key = fmt.value
        #     ['states', DegenerateSpaceInputFormat.Groups],
        #     ['nT', DegenerateSpaceInputFormat.QuantaSpecRules],
        #     ['energy_cutoff', DegenerateSpaceInputFormat.EnergyCutoff],
        #     ['polyads', DegenerateSpaceInputFormat.Polyads],
        #     ['couplings', DegenerateSpaceInputFormat.StrongCouplings],
        #     ['callable', DegenerateSpaceInputFormat.Callable]
        # ]:
            val = cls._key(spec, key)
            if val is not None:
                opts = cls._key(spec, 'options')
                if opts is None and isinstance(spec, dict):
                    opts = spec.copy()
                    del opts[key]
                return (val, opts), fmt
        else:
            for fmt,kls in cls.get_format_mapping().items():
                new_spec = kls.canonicalize(spec)
                if new_spec is not None:
                    return (new_spec, None), fmt
            else:
                return (None, None), None

    @abc.abstractmethod
    def get_groups(self, input_states, solver=None, **kwargs):
        """
        :param solver:
        :type solver:
        :param input_states:
        :type input_states:
        :return:
        :rtype:
        """

        raise NotImplementedError("abstract interface")

    @classmethod
    @abc.abstractmethod
    def canonicalize(cls, spec):
        raise NotImplementedError("abstract interface")

    def get_group_filter(self, **kwargs):
        return DegenerateMultiStateSpace.construct_filer(**kwargs)

class EnergyCutoffDegeneracySpec(DegeneracySpec):
    """

    """

    format = DegenerateSpaceInputFormat.EnergyCutoff
    def __init__(self, cutoff, **opts):
        super().__init__(**opts)
        self.cutoff = cutoff

    @classmethod
    def canonicalize(cls, spec):
        return isinstance(spec, (int, float, np.integer, np.floating))

    @classmethod
    def _group_states_by_energy_cutoff(cls, H0, states, cutoff):
        """
        :type H: Iterable[SparseArray]
        :type states: BasisStateSpace
        :type total_state_space: BasisMultiStateSpace
        :type cutoff: float
        :rtype: Iterable[BasisStateSpace]
        """
        # we look for states with energies within a range...
        # so initially we pull the sets of energies

        diag_inds = BraKetSpace(states, states)
        energies = H0[diag_inds, diag_inds]
        degenerate_groups = []
        # then we look through the input states
        for n, e in enumerate(energies):
            # we only want to apply this once per degenerate group
            # NOTE: this is a path to subtlety, since
            #   if state a is within 50 cm^-1 of state b, and state b is within of c,
            #   you might argue a and c are degenerate
            #   we are wagering that states are distinct _enough_ such that this is not
            #   an issue, but if it is a different strategy will be required
            if all(n not in d for d in degenerate_groups):
                e_diffs = np.abs(energies - e)
                inds = np.where(e_diffs < cutoff)[0]
                degenerate_groups.append(set(inds))
        # raise Exception(degenerate_groups)
        degenerate_groups = [states.take_subspace(np.array(list(d), dtype=int)) for d in degenerate_groups]
        return degenerate_groups

class MartinTestDegeneracySpec(DegeneracySpec):
    application_order = 'mid'
    group_filter = 'default'
    def __init__(self, threshold=4.6e-6, test_energy_window=4.6e-3, convert=True, frequencies=None, **opts):
        super().__init__(**opts)
        if convert and threshold > 1e-2: # help coerce it into Hartrees
            threshold = threshold / 219465 # only need a rough conversion...
        self.threshold = threshold
        self.window = test_energy_window
        self.frequencies = frequencies
        self._test_groups = None
        self._states = None
        self._basis = None
        self._matrix = None

    repr_opts = ['energy_cutoff', 'threshold']
    def prep_states(self, states:BasisStateSpace):
        state_list = states.excitations
        zo_engs = np.dot(state_list, self.frequencies)
        extended_states = states.basis.operator('x', 'x', 'x').get_transformed_space(states)
        key_inds = []
        exc_groups = []
        for i, (s, e0, exc) in enumerate(zip(state_list, zo_engs, extended_states)):
            e1 = np.dot(exc.excitations, self.frequencies)
            eng_diffs = e1 - e0
            pos = np.where(np.abs(eng_diffs < self.window))
            if len(pos) > 0 and len(pos[0]) > 0:
                key_inds.append(i)
                exc_groups.append(exc.take_subspace(pos[0]))
        if len(key_inds) > 0:
            self._test_groups = SelectionRuleStateSpace(
                states.take_subspace(key_inds),
                exc_groups
            )
            new = states.union(self._test_groups.to_single().take_unique())
        else:
            new = states
            self._test_groups = []
        return new
    def get_coupled_spaces(self, input_states:BasisStateSpace, solver=None):
        # raise NotImplementedError("...")
        if self._test_groups is None:
            self.prep_states(input_states)
        elif len(self._test_groups) == 0:
            return []
        if solver is None:
            raise ValueError("need a solver with fully evaluated perturbations first...")
        total_basis = solver.flat_total_space
        # engs = solver.zero_order_energies
        h1 = solver.representations[1]
        groups = self._test_groups #type:SelectionRuleStateSpace
        keys = groups.representative_space
        key_inds = total_basis.find(keys)
        spaces = []
        martin_inds = [[], []]
        martin_vals = []
        for i, (k, subspace) in enumerate(zip(key_inds, groups.flat)):
            e0 = np.dot(keys.excitations[i], self.frequencies)
            e1 = np.dot(subspace.excitations, self.frequencies)
            p = total_basis.find(subspace)
            e_diffs = np.abs(e1 - e0)
            repr_elems = h1[k, p].asarray()
            # diffs = subspace.excitations - keys.excitations[np.newaxis, k]
            # diff_sums = np.sum(np.abs(diffs), axis=1)
            # raise Exception(diffs, keys.excitations[np.newaxis, k], subspace.excitations)
            # weights = 1
            test_val = (repr_elems**4)/(e_diffs**3)
            # if (keys.excitations[i] == np.array([0, 0, 0, 1, 0, 0, 0, 0, 0])).all():
            #     raise Exception(test_val*219465, repr_elems, e_diffs, subspace.excitations)
            # if (keys.excitations[i] == np.array([0, 0, 0, 0, 0, 0, 1, 0, 0])).all():
            #     raise Exception(test_val*219465, repr_elems, e_diffs, subspace.excitations)
            martin_inds[0].extend([i]*len(subspace))
            martin_inds[1].extend(p)
            martin_vals.extend(test_val)
            pos = np.where(test_val > self.threshold)
            if len(pos) > 0 and len(pos[0]) > 0:
                # print(pos[0], subspace.excitations)
                new = np.concatenate([
                        keys.excitations[i:i+1],
                        subspace.excitations[pos[0],]#.take_subspace(pos[0]).excitations
                    ], axis=0)
                # print(keys.excitations[i:i+1], new)
                new_sums = np.sum(new, axis=1)
                new = new[new_sums > 0]
                spaces.append(new)

        self._states = keys
        self._basis = solver.flat_total_space
        self._matrix = nput.SparseArray.from_data(
            (
                np.array(martin_vals),
                tuple(np.array(x) for x in martin_inds)
            ),
            shape=(len(self._states), len(self._basis))
        )
        # raise Exception("...")
        return spaces
    # @classmethod
    def get_groups(self, input_states:BasisStateSpace, solver=None, **kwargs):
        groups = self.get_coupled_spaces(input_states, solver=solver)
        indexer = SymmetricGroupGenerator(input_states.ndim)
        groups = PermutationRelationGraph.merge_groups([(indexer.to_indices(g), g) for g in groups])
        return [g[1] for g in groups]

    @classmethod
    def canonicalize(cls, spec):
        return isinstance(spec.threshold, float)

    def get_group_filter(self, threshold=None, states=None, basis=None, corrections=None, **kwargs):
        return DegenerateMultiStateSpace.construct_filer(
            states=self._states,
            basis=self._basis,
            corrections=self._matrix,
            threshold=self.threshold,
            threshold_step_size=1/219465,
            **kwargs
        )

class StronglyCoupledDegeneracySpec(DegeneracySpec):
    """

    """
    application_order = 'post'
    format = DegenerateSpaceInputFormat.StrongCouplings
    default_threshold=.3
    def __init__(self, wfc_threshold=None, state_filter=None, extend_spaces=True, **opts):
        super().__init__(**opts)
        if wfc_threshold is None or isinstance(wfc_threshold, str) and wfc_threshold == 'auto':
            wfc_threshold = self.default_threshold
        self.wfc_threshold = wfc_threshold
        self.state_filter = state_filter
        self.extend_spaces=extend_spaces

    repr_opts = ['energy_cutoff', 'wfc_threshold']
    def prep_states(self, input_states):
        return input_states

    def identify_strong_couplings(self, solver, corrs):
        degenerate_correction_threshold = self.wfc_threshold
        test_modes = self.test_modes
        if test_modes is None:
            test_modes = solver.high_frequency_modes

            # raise Exception(test_modes, test_freqs * 219465, fundamentals.excitations, where)

        state_filter = self.state_filter
        zero_order_energy_cutoff = self.energy_cutoff
        if state_filter is None:
            state_filter = lambda state, couplings: corrs.default_state_filter(state, couplings,
                                                                               energy_cutoff=zero_order_energy_cutoff,
                                                                               energies=solver.zero_order_energies,
                                                                               basis=solver.flat_total_space,
                                                                               target_modes=test_modes
                                                                               )
        sc = corrs.find_strong_couplings(threshold=degenerate_correction_threshold, state_filter=state_filter)
        return sc

    def get_groups(self, input_states, couplings=None, solver=None, **kwargs):
        """
        :param input_states:
        :type input_states:
        :param solver:
        :type solver:
        :return:
        :rtype:
        """
        if couplings is None:
            raise ValueError("need couplings")
        return self.get_strong_coupling_space(input_states, couplings)

    @classmethod
    def canonicalize(cls, spec):
        if isinstance(spec, (int, float, np.integer, np.floating)):
            return spec
        else:
            try:
                return {k:GroupsDegeneracySpec._validate_grp(v) for k,v in spec.items()}
            except (np.VisibleDeprecationWarning, AttributeError, StopIteration):
                return None

    @classmethod
    def get_strong_coupling_space(cls, states: BasisStateSpace, couplings: dict):
        indexer = SymmetricGroupGenerator(states.ndim)
        groups = [None] * len(states.indices)
        for n, i in enumerate(states.indices):
            if i in couplings:
                groups[n] = [i] + list(couplings[i].indices)
            else:
                groups[n] = [i]

        groups = PermutationRelationGraph.merge_groups([(g, indexer.from_indices(g)) for g in groups])

        return [g[1] for g in groups]

class GroupsDegeneracySpec(DegeneracySpec):
    """

    """

    format = DegenerateSpaceInputFormat.Groups
    def __init__(self, groups, **opts):
        super().__init__(**opts)
        self.groups = groups

    @staticmethod
    def _validate_grp(grp):
        try:
            b = grp[0]
        except (IndexError, TypeError):
            raise StopIteration("malformatted")
        else:
            return np.asanyarray(grp)
    @classmethod
    def canonicalize(cls, spec):
        try:
            return [cls._validate_grp(s) for s in spec]
        except (np.VisibleDeprecationWarning, AttributeError, StopIteration):
            return None

    def get_groups(self, input_states, solver=None, **kwargs):
        """
        :param solver:
        :type solver:
        :param input_states:
        :type input_states:
        :return:
        :rtype:
        """
        # do I want to canonicalize if group size is too small?
        return self.groups

class PolyadDegeneracySpec(DegeneracySpec):
    """

    """

    format = DegenerateSpaceInputFormat.Polyads
    def __init__(self, polyads,
                 max_quanta=None, max_iterations=2,
                 require_converged=False, extra_groups=None,
                 **opts):
        super().__init__(**opts)
        self.polyads = polyads
        self.max_quanta = max_quanta
        self.max_iterations = max_iterations
        self.require_converged= require_converged
        self.extra_groups = extra_groups

    @staticmethod
    def _validate_rule(pair):
        try:
            if len(pair) != 2:
                raise StopIteration("malformatted")
            b = pair[0]
        except (IndexError, TypeError):
            raise StopIteration("malformatted")
        else:
            if len(pair[0]) == len(pair[1]):
                return pair
            else:
                raise StopIteration("malformatted")
    @classmethod
    def canonicalize(cls, spec):
        try:
            return [cls._validate_rule(s) for s in spec]
        except (np.VisibleDeprecationWarning, AttributeError, StopIteration):
            return None

    def get_groups(self, input_states, solver=None, **kwargs):
        """
        :param solver:
        :type solver:
        :param input_states:
        :type input_states:
        :return:
        :rtype:
        """
        if isinstance(input_states, BasisStateSpace):
            input_states = input_states.excitations
        return self.get_degenerate_polyad_space(
            input_states,
            self.polyads,
            max_quanta=self.max_quanta,
            max_iterations=self.max_iterations,
            require_converged=self.require_converged,
            extra_groups=self.extra_groups
        )

    @classmethod
    def get_degenerate_polyad_space(cls, states, polyadic_pairs,
                                    max_quanta=None, max_iterations=2,
                                    require_converged=False, extra_groups=None):
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

        graph = PermutationRelationGraph(polyadic_pairs)
        groups = graph.build_state_graph(states,
                                         extra_groups=extra_groups,
                                         max_sum=max_quanta,
                                         max_iterations=max_iterations,
                                         raise_iteration_error=require_converged
                                         )

        grp = [g for g in groups if len(g) > 1]
        return grp

    @staticmethod
    def _is_polyad_rule(d, n_modes):
        try:
            return (
                    len(d) == 2
                    and len(d[0]) == n_modes
                    and len(d[1]) == n_modes
            )
        except TypeError:
            return False

class TotalQuantaDegeneracySpec(PolyadDegeneracySpec):
    """

    """

    format = DegenerateSpaceInputFormat.QuantaSpecRules
    def __init__(self, n_T_vectors, max_quanta=3, **opts):
        self.max_quanta = max_quanta
        self.nt_vecs = [n_T_vectors] if isinstance(n_T_vectors[0], int) else n_T_vectors
        rules = sum((self.make_nt_polyad(v, max_quanta=max_quanta) for v in self.nt_vecs), [])
        # raise Exception(rules)
        super().__init__(rules, **opts)

    @classmethod
    def reduce_rule(cls, a, b):
        a = a.copy()
        b = b.copy()

        m1 = np.min(a[a>0])
        m2 = np.min(b[b>0])
        if m2 > m1:
            m1, m2 = m2, m1
        if m1 % m2 == 0:
            a = a//m2
            b = b//m2
        reduce_pos = np.where(np.logical_and(a > 0, b > 0))
        if len(reduce_pos) > 0:
            reduce_pos = reduce_pos[0]
            if len(reduce_pos) > 0:
                pairs = np.concatenate([a[reduce_pos][:, np.newaxis], b[reduce_pos][:, np.newaxis]], axis=1)
                min_vals = np.min(pairs, axis=1)
                a[reduce_pos] = a[reduce_pos] - min_vals
                b[reduce_pos] = b[reduce_pos] - min_vals
        return a,b

    @classmethod
    def make_nt_polyad(cls, nt, max_quanta=3):  # this should be of the form num of quanta give same E
        # we want to choose different subsets of the nT vector to construct rules
        # we'll do this in a kinda dumb way where we force a max number of quanta for these polyad
        # rules and explicitly test states of a given form
        states = BasisStateSpace.from_quanta(
            HarmonicOscillatorProductBasis(len(nt)),
            list(range(1, max_quanta+1))
        ).excitations

        nts = np.dot(states, nt)
        (keys, groups), _ = nput.group_by(states, nts)
        rules = set()

        for k,g in zip(keys, groups):
            if k > 0:
                for n,i in enumerate(g):
                    for j in g[n+1:]:
                        a, b = cls.reduce_rule(i, j)
                        rules.add((tuple(a), tuple(b)))
        return list(rules)

class CallableDegeneracySpec(DegeneracySpec):
    """

    """

    format = DegenerateSpaceInputFormat.Callable
    def __init__(self, callable, **opts):
        super().__init__(**opts)
        self.callable = callable

    def get_groups(self, input_states, solver=None, **kwargs):
        """
        :param input_states:
        :type input_states:
        :param solver:
        :type solver:
        :return:
        :rtype:
        """
        return self.callable(input_states, solver=solver, **kwargs)

class DegenerateMultiStateSpace(BasisMultiStateSpace):

    @staticmethod
    def default_group_filter(group,
                             states=None,
                             basis=None,
                             corrections=None,
                             threshold=None,
                             energy_cutoff=None,
                             energies=None,
                             decoupling_overide=-1,
                             maximize_groups=True,
                             target_modes=None,
                             threshold_step_size=.1
                             ):
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

        if len(group) < 2:
            return group

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
            elif hasattr(corrections, 'total_basis'): # PT corrections
                states = corrections.states
                basis = corrections.total_basis
                group_inds = corrections.total_basis.find(group)
                state_inds = corrections.states.find(group, dtype=int, missing_val=len(states))
                found_pos = state_inds < len(states)
                kill_spots = (np.arange(len(group_inds))[found_pos], state_inds[found_pos])
                corr_mat = [np.abs(x[:, group_inds].asarray().T) for x in corrections.wfn_corrections[1:]]
                max_corr = None
                for c in corr_mat:
                    c[kill_spots] = 0
                    if max_corr is None:
                        max_corr = np.max(c, axis=1)
                    else:
                        submax = np.max(c, axis=1)
                        p = submax > max_corr
                        max_corr[p] = submax[p]
                corrections = sum(corr_mat)
                if threshold is not None:
                    sorting = np.arange(len(max_corr))
                else:
                    sorting = np.argsort(-max_corr - 1 * (found_pos))  # prioritize states inside requested packet for marginal cases
            elif basis is not None and states is not None:
                group_inds = basis.find(group)
                state_inds = states.find(group, dtype=int, missing_val=len(states))
                found_pos = state_inds < len(states)
                found_states = state_inds[found_pos]
                kill_spots = (np.arange(len(group_inds))[found_pos], found_states)
                corrections = np.abs(corrections[:, group_inds].asarray().T)
                max_corr = None
                corrections[kill_spots] = 0
                # print(group.excitations)
                # print(states.excitations)
                # raise Exception(corrections)
                if max_corr is None:
                    max_corr = np.max(corrections, axis=1)
                else:
                    submax = np.max(corrections, axis=1)
                    p = submax > max_corr
                    max_corr[p] = submax[p]
                if threshold is not None:
                    sorting = np.arange(len(max_corr))
                else:
                    sorting = np.argsort(-max_corr - 1*found_pos)  # prioritize states inside requested packet for marginal cases

        else:
            sorting = np.argsort(inds)
            corrections = None

        exc = exc[sorting]
        diffs = exc[:, np.newaxis] - exc[np.newaxis, :]  # difference matrix (s, s, m)
        tots = exc[:, np.newaxis] + exc[np.newaxis, :]  # difference matrix (s, s, m)
        diff_sums = np.sum(np.abs(diffs), axis=2)  # (s, s)
        tots_sums = np.sum(np.abs(tots), axis=2)  # (s, s)
        elims = []
        if threshold is not None and corrections is not None:
            if energy_cutoff is None:
                bad_pos = np.where(np.logical_and(
                    diff_sums == 1,
                    tots_sums > 0
                ))  # np.logical_or(diff_sums == 1, diff_sums == 0))
            else:
                e_vec = energies[group_inds]
                bad_pos = np.where(np.abs(np.subtract.outer(e_vec, e_vec)) > energy_cutoff)
                # print(">>>", group.excitations)
                # print(bad_pos)
            if len(bad_pos) > 0 and len(bad_pos[0]) > 0:
                N = len(group_inds)
                base_corr = corrections
                base_vals = base_corr[:, state_inds[found_pos]]
                base_mat = np.zeros((N, N), dtype=float)
                base_mat[:, np.where(found_pos)[0]] = base_vals
                base_mat = (base_mat.T + base_mat) - 2*np.diag(np.diag(base_mat))
                if decoupling_overide > 0:
                    truly_bad = base_mat[bad_pos] < decoupling_overide
                else:
                    truly_bad = np.full(len(base_mat[bad_pos]), True)
                # introduce an override so things that are super strongly coupled can remain coupled
                bad_pos = tuple(b[truly_bad] for b in bad_pos)
                base_mat[bad_pos] = 0  # refuse to couple problem states
                base_mat[bad_pos[1], bad_pos[0]] = 0  # refuse to couple problem states

                def get_groups(cm):
                    _, labs = sp.csgraph.connected_components(cm, return_labels=True)
                    (grp_keys, groups), _ = nput.group_by(np.arange(N), labs)
                    A = np.zeros((N, N))
                    for g in groups:
                        A[np.ix_(g, g)] = 1
                    return groups, A

                corr_mat = (base_mat > threshold).astype(int)
                groups, A = get_groups(corr_mat)

                target_threshold = threshold
                step_size = threshold_step_size
                while (A[bad_pos] > 0).any():
                    threshold = threshold + step_size
                    corr_mat = (base_mat > threshold).astype(int)
                    groups, A = get_groups(corr_mat)

                if maximize_groups:
                    dropped_pos = np.where(np.logical_and(base_mat > target_threshold, base_mat < threshold))
                    if len(dropped_pos) > 0:
                        # figure out if any can be added without causing breakdowns
                        sorting = np.argsort(base_mat[dropped_pos])
                        dropped_pos = tuple(d[sorting] for d in dropped_pos)
                        for i,j in zip(*dropped_pos): # this could be very expensive...
                            corr_mat[i, j] = base_mat[i, j]
                            _, A = get_groups(corr_mat)
                            if (A[bad_pos] > 0).any():
                                corr_mat[i, j] = 0.
                            else:
                                groups = _

                group = [group.take_subspace(g) for g in groups if len(g) > 1]

        else:
            for n in np.arange(len(exc)):
                diff_vec = diff_sums[n]
                bad_pos = np.where(diff_vec == 1)
                if len(bad_pos) > 0:
                    bad_pos = bad_pos[0]
                    if len(bad_pos) > 0:
                        kills = bad_pos[bad_pos > n]
                        elims.extend(sorting[kills])
                        diff_sums[bad_pos, :] = 0

            if len(elims) > 0:
                mask = np.setdiff1d(np.arange(len(exc)), elims)
                group = group.take_subspace(mask)

        return group

    @classmethod
    def construct_filer(cls, spec=None, **kwargs):
        if spec is not None:
            return spec.get_group_filter(**kwargs)
        else:
            return lambda group:cls.default_group_filter(group, **kwargs)

    @classmethod
    def from_spec(cls,
                  degenerate_states,
                  solver=None,
                  full_basis=None,
                  format=None,
                  group_filter=None,
                  log_groups=False,
                  **kwargs
                  ):
        """
        Generates a DegenerateMultiStateSpace object from a number
        of possible specs

        :param solver: the actual applier of the perturbation theory which makes use of the degenerate states
        :type solver: PerturbationTheorySolver
        :return:
        :rtype:
        """

        logger = solver.logger
        states = solver.states #type: BasisStateSpace

        if degenerate_states is not None:
            if not isinstance(degenerate_states, DegeneracySpec):
                degenerate_states = DegeneracySpec.from_spec(degenerate_states, format=format)
            with logger.block(tag="getting degeneracies"):
                deg_states = degenerate_states.get_groups(states, solver=solver, **kwargs)
                if hasattr(group_filter, 'items'):
                    group_filter = cls.construct_filer(**dict(group_filter, spec=degenerate_states))
                new = GroupsDegeneracySpec.canonicalize(deg_states)
                if new is None:
                    raise ValueError("{} returned invalid degenerate subspace spec {}".format(
                        degenerate_states,
                        deg_states
                    ))
                else:
                    degenerate_states = new
                    del deg_states

                degenerate_states = [
                    BasisStateSpace(states.basis, s)#.as_sorted(track_indices=False, track_excitations=False)
                    if not isinstance(s, BasisStateSpace) else s
                    for s in degenerate_states
                ]

                if log_groups:
                    with logger.block(tag="initial degenerate groups:"):
                        for g in degenerate_states:
                            if len(g) > 1:
                                logger.log_print(str(g.excitations).splitlines())
                else:
                    logger.log_print(
                        "{n} degenerate state sets found",
                        n=len([x for x in degenerate_states if len(x) > 1]),
                        # s=[x for x in degenerate_states if len(x) > 1]
                    )

        # build groups of degenerate states for use later
        if degenerate_states is None:
            groups = states.split(1) #[[x] for x in states.indices]  # we're gonna loop through this later so why not destructure now...
        else:
            groups = [None] * len(degenerate_states)
            deg_sets = [(set(d.indices),d) for d in degenerate_states]
            for x in states.indices:
                for i, (d,bs) in enumerate(deg_sets):
                    # if i == 11:
                    #     raise Exception(bs.excitations, bs.indices)
                    if x in d:
                        if groups[i] is None:
                            groups[i] = bs.indices
                        # groups[i].append(x)
                        break
                else:
                    groups.append([x])

        # now turn these into proper BasisStateSpace objects so we can work with them more easily
        ugh = []
        for i,g in enumerate(groups):
            if g is None:
                continue #???
            # g = np.sort(np.array(g))
            if not isinstance(g, BasisStateSpace):
                # if np.any(np.asanyarray(g) < 0):
                #     raise Exception(i, g)
                g = BasisStateSpace(states.basis, np.array(g), mode=BasisStateSpace.StateSpaceSpec.Indices, full_basis=full_basis)
            if group_filter is not None:
                g = group_filter(g)
                if isinstance(g, BasisStateSpace):
                    ugh.append(g)
                else:
                    # for gg in g:
                    #     print(gg.excitations, gg)
                    ugh.extend(g)
            else:
                ugh.append(g)

        for n,x in enumerate(states.indices):
            for g in ugh:
                if x in g.indices:
                    break
            else:
                ugh.append(states.take_subspace([n]))

        # now make a numpy array for initialization
        arrs = np.full(len(ugh), None)
        for i, g in enumerate(ugh):
            arrs[i] = g

        # raise Exception(arrs)

        return cls(arrs)





