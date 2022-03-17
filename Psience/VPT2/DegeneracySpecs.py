

import numpy as np, enum, abc
from McUtils.Combinatorics import SymmetricGroupGenerator, PermutationRelationGraph
from ..BasisReps import BasisStateSpace, BasisMultiStateSpace, BraKetSpace

__all__ = [
    "DegeneracySpec",
    "DegenerateMultiStateSpace"
]

__reload_hook__ = ["..BasisReps"]

class DegenerateSpaceInputFormat(enum.Enum):
    Groups = "groups"
    QuantaSpecRules = "n_T"
    Polyads = "polyads"
    StrongCouplings = "strong_couplings"
    EnergyCutoff = "energy_cutoff"
    Callable = "callable"

class DegeneracySpec(metaclass=abc.ABCMeta):
    """
    Provides a container for specifying degeneracies
    in a way that can be cleanly canonicalized
    """

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
        }

    @classmethod
    def from_spec(cls, spec, format=None, **kwargs):
        if spec is None:
            return None

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

        for key,fmt in [
            ['states', DegenerateSpaceInputFormat.Groups],
            ['nT', DegenerateSpaceInputFormat.QuantaSpecRules],
            ['energy', DegenerateSpaceInputFormat.EnergyCutoff],
            ['polyads', DegenerateSpaceInputFormat.Polyads],
            ['couplings', DegenerateSpaceInputFormat.StrongCouplings],
            ['callable', DegenerateSpaceInputFormat.Callable]
        ]:
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
    def get_groups(self, input_states, solver=None):
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


class TotalQuantaDegeneracySpec(DegeneracySpec):
    """

    """

    format = DegenerateSpaceInputFormat.QuantaSpecRules
    def __init__(self, n_T_vector):
        self.n_T = n_T_vector

    def get_groups(self, input_states, solver=None):
        return self._group_states_by_nt_spec(input_states, self.n_T)
    @classmethod
    def canonicalize(cls, n_T_vector):
        try:
            test_el = n_T_vector[0]
        except TypeError:
            test_el = None
        if isinstance(test_el, (int, np.integer)):
            return np.asanyarray(n_T_vector)
        else:
            return None

    @classmethod
    def _group_states_by_nt_spec(cls, states, q_vec):
        """
        :type H: Iterable[SparseArray]
        :type states: BasisStateSpace
        :type total_state_space: BasisMultiStateSpace
        :type cutoff: Iterable[int]
        :rtype: Iterable[BasisStateSpace]
        """
        # we build the total N_t to compare against once...
        tot_n_t = np.dot(states.excitations, q_vec)
        degenerate_groups = {}
        # then we look through the input states
        for vec in states.excitations:
            # base n_t
            n_t = np.dot(q_vec, vec)
            if n_t not in degenerate_groups:
                degenerate_groups[n_t] = np.where(tot_n_t == n_t)[0]
        degenerate_groups = [states.take_subspace(np.array(d)) for d in degenerate_groups.values()]
        return degenerate_groups

class EnergyCutoffDegeneracySpec(DegeneracySpec):
    """

    """

    format = DegenerateSpaceInputFormat.EnergyCutoff
    def __init__(self, cutoff):
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

class GroupsDegeneracySpec(DegeneracySpec):
    """

    """

    format = DegenerateSpaceInputFormat.Groups
    def __init__(self, groups):
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

    def get_groups(self, input_states, solver=None):
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
                 require_converged=False, extra_groups=None
                 ):
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

    def get_groups(self, input_states, solver=None):
        """
        :param solver:
        :type solver:
        :param input_states:
        :type input_states:
        :return:
        :rtype:
        """
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

        return [g for g in groups if len(g) > 1]

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

class StronglyCoupledDegeneracySpec(DegeneracySpec):
    """

    """

    format = DegenerateSpaceInputFormat.StrongCouplings
    def __init__(self, couplings):
        self.couplings = couplings

    def get_groups(self, input_states, solver=None):
        """
        :param input_states:
        :type input_states:
        :param solver:
        :type solver:
        :return:
        :rtype:
        """
        return self.get_strong_coupling_space(input_states, self.couplings)

    @classmethod
    def canonicalize(cls, spec):
        try:
            return {k:GroupsDegeneracySpec._validate_grp(g) for k,v in spec.items()}
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

class CallableDegeneracySpec(DegeneracySpec):
    """

    """

    format = DegenerateSpaceInputFormat.Callable
    def __init__(self, callable):
        self.callable = callable

    def get_groups(self, input_states, solver=None):
        """
        :param input_states:
        :type input_states:
        :param solver:
        :type solver:
        :return:
        :rtype:
        """
        return self.callable(input_states, solver=solver)

class DegenerateMultiStateSpace(BasisMultiStateSpace):

    @staticmethod
    def default_group_filter(group, target_modes=None):
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
        if len(group) == 1:
            return group

        if target_modes is not None:
            sub = group.take_subdimensions(target_modes)
            inds = sub.indices
            exc = sub.excitations
        else:
            inds = group.indices
            exc = group.excitations
        sorting = np.argsort(inds)
        exc = exc[sorting]
        diffs = exc[:, np.newaxis] - exc[np.newaxis, :]  # difference matrix (s, s, m)
        diff_sums = np.sum(diffs != 0, axis=2)  # (s, s)
        bad_pos = np.where(diff_sums == 1) # np.logical_or(diff_sums == 1, diff_sums == 0))
        if len(bad_pos) > 0 and len(bad_pos[0]) > 0:
            kills = np.unique(bad_pos[1][bad_pos[1] > bad_pos[0]]) # upper triangle
            # prek = kills
            kills = sorting[kills] # OG indices
            # raise Exception(prek, kills, exc[prek], group.excitations[kills])

            # now drop all of these from the total space
            mask = np.setdiff1d(np.arange(len(exc)), kills)
            group = group.take_subspace(mask)
            # raise Exception(group.excitations, kills)#, np.array(bad_pos).T)

        return group

    @classmethod
    def from_spec(cls,
                  degenerate_states,
                  solver=None,
                  full_basis=None,
                  format=None,
                  group_filter=None
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
                deg_states = degenerate_states.get_groups(states, solver=solver)
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
                    BasisStateSpace(states.basis, s) if not isinstance(s, BasisStateSpace) else s
                    for s in degenerate_states
                ]

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
                    if x in d:
                        if groups[i] is None:
                            groups[i] = bs.indices
                        # groups[i].append(x)
                        break
                else:
                    groups.append([x])

        # now turn these into proper BasisStateSpace objects so we can work with them more easily
        ugh = [None]*len(groups)
        for i,g in enumerate(groups):
            # g = np.sort(np.array(g))
            if not isinstance(g, BasisStateSpace):
                g = BasisStateSpace(states.basis, np.array(g), mode=BasisStateSpace.StateSpaceSpec.Indices, full_basis=full_basis)
            if group_filter is not None:
                g = group_filter(g)
            ugh[i] = g
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
            # if len(g) > 1:
            #     raise Exception(ugh[i].indices, g, ugh[i].excitations,
            #             states.basis.unravel_state_inds(np.arange(10)))

        return cls(arrs)





