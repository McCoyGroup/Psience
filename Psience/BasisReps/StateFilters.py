
import numpy as np, abc, math, collections
from McUtils.Combinatorics import UniquePermutations
from McUtils.Data import UnitsData
from .StateSpaces import BasisStateSpace, BasisMultiStateSpace, SelectionRuleStateSpace

__all__ = [
    "BasisStateSpaceFilter", "PropertyCorrectionFilter"
]

class Postfilter(metaclass=abc.ABCMeta):
    native_keys = frozenset({'partitions','parity','min_quanta','max_quanta',
                             'energy_windows','energy_units','frequencies'})

    def get_permutation_constraints(self):
        """Optional numeric constraint description; None preserves apply()."""
        return None

    def __init__(self, input_space=None):
        self.input_space = input_space
    @abc.abstractmethod
    def apply(self, state_space:BasisStateSpace)->BasisStateSpace:
        pass
    @classmethod
    def get_spec_map(cls):
        return {
            'test':ExclusionTestFilter,
            'max_quanta':MaxQuantaFilter,
            'intersected_space':IntersectionSpaceFilter,
            'excluded_space':DifferenceSpaceFilter,
            'excluded_transitions':ExcludedTransitionFilter,
        }
    @classmethod
    def construct(cls, input_space, spec):
        if isinstance(spec,Postfilter):
            return spec
        if isinstance(spec, BasisStateSpace):
            spec = {'intersected_space':spec}
        if any(key in cls.native_keys for key in spec):
            if set(spec) == {'max_quanta'}:
                return MaxQuantaFilter(input_space=input_space,**spec)
            unknown = set(spec)-cls.native_keys-{'test'}
            if unknown:
                raise ValueError('combine native and set/callback filters as a list; unknown keys: {}'.format(sorted(unknown)))
            return ExclusionTestFilter(input_space=input_space,**spec)
        for k,c in cls.get_spec_map().items():
            if k in spec:
                return c(input_space=input_space, **spec)
        else:
            raise ValueError("Don't know how to build filter from spec {}".format(spec))

class ExclusionTestFilter(Postfilter):
    def __init__(self, test=None, input_space=None, *, partitions=None, parity=None,
                 min_quanta=None, max_quanta=None, energy_windows=None,
                 energy_units='Hartrees', frequencies=None, **opts):
        super().__init__(input_space=input_space, **opts)
        self.test = test
        self._constraints = None
        if any(value is not None for value in [partitions,parity,min_quanta,max_quanta,energy_windows]):
            if input_space is None:
                raise ValueError('structured basis filters require an input space')
            constraints = dict(partitions=partitions,parity=parity,
                               minimum=min_quanta,maximum=max_quanta)
            if energy_windows is not None:
                if frequencies is None:
                    raise ValueError('energy_windows require harmonic frequencies; run_simple binds these automatically')
                constraints.update(weights=np.asarray(frequencies,dtype=float),
                    windows=np.asarray(energy_windows,dtype=float)*UnitsData.convert(energy_units,'Hartrees'))
            elif frequencies is not None or energy_units != 'Hartrees':
                raise ValueError('energy_units/frequencies require energy_windows')
            self._constraints = UniquePermutations.prepare_constraints(constraints,input_space.ndim)
        elif test is None:
            raise ValueError('a filter requires a test or structured constraints')

    def get_permutation_constraints(self):
        return self._constraints if self.test is None else None

    def apply(self, state_space:BasisStateSpace) ->BasisStateSpace:
        if self._constraints is not None:
            if self.test is None and getattr(state_space,'store',None) is not None:
                group = getattr(state_space.basis.indexer,'symmetric_group',None)
                if group is not None:
                    return state_space.store.frontier(state_space.basis,
                        group.iter_filtered_indices(state_space.indices,self._constraints,
                                                    state_space.store.block_size))
            vals = UniquePermutations.test_constraints(state_space.excitations,self._constraints)
            if self.test is not None:
                vals &= self.test(state_space)
        else:
            vals = self.test(state_space)
        return state_space.take_subspace(np.flatnonzero(vals))


class MaxQuantaFilter(ExclusionTestFilter):
    def __init__(self, max_quanta=None, input_space:BasisStateSpace=None, **opts):
        if isinstance(max_quanta, (int, np.integer)):
            max_quanta = [max_quanta] * input_space.ndim
        self.max_quanta = np.asanyarray(max_quanta)
        super().__init__(self.check_mode_quanta, input_space=input_space, **opts)
    def get_permutation_constraints(self):
        return UniquePermutations.prepare_constraints(
            {'maximum':self.max_quanta},self.input_space.ndim)

    def check_mode_quanta(self, state_space):
        if len(state_space) == 0:
            return np.array([], dtype=bool)
        exc = state_space.excitations
        qs = self.max_quanta.copy()
        qs[qs < 0] = np.max(exc)+1
        return np.all(exc <= qs[np.newaxis], axis=1)
class IntersectionSpaceFilter(Postfilter):
    def __init__(self, intersected_space=None, **opts):
        super().__init__(**opts)
        self.space = intersected_space
    def apply(self, state_space:BasisStateSpace) ->BasisStateSpace:
        return state_space.intersection(self.space)
class DifferenceSpaceFilter(Postfilter):
    def __init__(self, excluded_space=None, **opts):
        super().__init__(**opts)
        self.space = excluded_space
    def apply(self, state_space:BasisStateSpace) ->BasisStateSpace:
        return state_space.difference(self.space)
class ExcludedTransitionFilter(DifferenceSpaceFilter):
    def __init__(self, excluded_transitions=None, input_space=None, **opts):
        bad_transitions = np.asanyarray(excluded_transitions)
        perms = (
                input_space.excitations[:, np.newaxis, :]
                + bad_transitions[np.newaxis, :, :]
        ).reshape(-1, bad_transitions.shape[-1])
        perms = perms[np.all(perms > 0, axis=1)]
        excluded_space = BasisStateSpace(input_space.basis, perms)
        super().__init__(excluded_space=excluded_space, input_space=input_space, **opts)
    def apply(self, state_space:BasisStateSpace) ->BasisStateSpace:
        return state_space.difference(self.space)

class BasisStateSpaceFilter:
    """
    Provides an easier constructor for the VPT state space filters
    """

    def __init__(self, input_space, prefilters, postfilters):
        """
        :param input_space:
        :type input_space: BasisStateSpace
        :param prefilters:
        :type prefilters:
        :param postfilters:
        :type postfilters:
        """
        self.input_space = input_space
        self._prefilters = None
        self._raw_prefilters = prefilters
        self._postfilters = None
        self._raw_postfilters = postfilters

    @classmethod
    def bind_postfilters(cls, filters, *, modes=None, mode_selection=None, mode_transformation=None):
        """Bind automatic harmonic weights without mutating caller filter specs."""
        if filters is None or isinstance(filters,BasisStateSpace) or isinstance(filters,Postfilter):
            return filters
        if isinstance(filters,dict):
            stage_map = any(isinstance(k,tuple) and len(k)==2 for k in filters)
            if stage_map:
                return {key:cls.bind_postfilters(value,modes=modes,mode_selection=mode_selection,
                                                mode_transformation=mode_transformation)
                        for key,value in filters.items()}
            result = dict(filters)
            if result.get('energy_windows') is not None and result.get('frequencies') is None:
                if mode_transformation is not None:
                    raise ValueError('energy_windows with transformed modes require explicit frequencies in Hartrees')
                if modes is None:
                    raise ValueError('automatic energy-window frequencies require molecular modes')
                frequencies = np.asarray(modes.freqs,dtype=float)
                selection = mode_selection
                if isinstance(selection,dict):
                    if selection.get('modes') is not None:
                        frequencies = frequencies[np.asarray(selection['modes'],dtype=int)]
                    selection = selection.get('derivatives')
                if selection is not None:
                    frequencies = frequencies[selection if isinstance(selection,slice) else np.asarray(selection,dtype=int)]
                result['frequencies'] = np.array(frequencies,copy=True)
            return result
        return [cls.bind_postfilters(item,modes=modes,mode_selection=mode_selection,
                                     mode_transformation=mode_transformation) for item in filters]

    @classmethod
    def from_postfilters(cls, input_space, postfilters, order=2):
        """Apply explicit postfilters without property-quanta restrictions."""
        if not isinstance(input_space,BasisStateSpace):
            from .HarmonicOscillator import HarmonicOscillatorProductBasis
            input_space = BasisStateSpace(HarmonicOscillatorProductBasis(len(input_space[0])),input_space)
        stages = isinstance(postfilters,dict) and any(isinstance(k,tuple) and len(k)==2 for k in postfilters)
        rules = []
        for h in range(1,order+1):
            for w in range(order-h+1):
                filters = postfilters.get((h,w)) if stages else postfilters
                if filters is not None:
                    rules.append(dict(order=h+w,expansion_term=h,prefilters=None,postfilters=filters))
        return cls.from_rules(input_space,*rules)


    @classmethod
    def from_data(cls, input_space, data):
        """
        Works to canonicalize inputs and initialize appropriately from there

        :param data:
        :type data:
        :return:
        :rtype:
        """

        if data is None:
            return None
        elif isinstance(data, dict):
            if not all(isinstance(v, cls) for v in data.values()):
                return cls._from_old_style_rules(input_space, data)
            else:
                return data
        else:
            return cls.from_rules(input_space, *data)

    @classmethod
    def _from_old_style_rules(cls, input_space, rules):
        """
        Builds a set of filter spaces from a dict keyed by Hamiltonian/wfn pairs
        and giving prefilters/postfilters

        :param rules:
        :type rules:
        :return:
        :rtype:
        """
        new_rules = []
        for k,v in rules.items():
            if not isinstance(v, dict):
                v = {'prefilters':v}
            new_rules.append(dict({'order':k[0]+k[1], 'expansion_term':k[0]}, **v))

        return cls.from_rules(input_space, *new_rules)


    @classmethod
    def from_rules(cls, input_space, *rules):
        """
        Builds a set of filter spaces from dicts of rules

        :param rules:
        :type rules:
        :return:
        :rtype:
        """

        filters = {}
        for r in rules:
            if not isinstance(r, dict):
                raise ValueError("filter rule must be a dict with `order`, `expansion_term`"
                                 ", and one or both of `prefilters` and `postfilter` for keys")
            for k in ["order", "expansion_term"]:
                if k not in r:
                    raise ValueError("filter rule needs `order`, `expansion_term`"
                                     ", and one or both of `prefilters` and `postfilter` for keys")
            if 'prefilters' not in r and 'postfilters' not in r and 'postfilter' not in r:
                raise ValueError("filter rule needs one or both of `prefilters` and `postfilter` for keys")

            h_term = r['expansion_term']
            y_term = r['order'] - h_term

            filters[(h_term, y_term)] = cls(
                input_space,
                r['prefilters'] if 'prefilters' in r else None,
                r.get('postfilters',r.get('postfilter'))
            )

        return filters

    @property
    def prefilters(self):
        """

        :return:
        :rtype:
        """
        if self._prefilters is None:
            self._prefilters = self.canonicalize_prefilters(self.input_space.basis,
                                                            self._raw_prefilters)
        return self._prefilters

    @property
    def postfilters(self):
        """

        :return:
        :rtype:
        """
        if self._postfilters is None:
            self._postfilters = self.canonicalize_postfilters(self.input_space,
                                                            self._raw_postfilters)
        return self._postfilters
    @classmethod
    def canonicalize_postfilters(self, input_space, filters):
        if filters is None:
            return None
        else:
            if hasattr(filters, 'items') or isinstance(filters, (BasisStateSpace,Postfilter)):
                filters = [filters]
            return [
                Postfilter.construct(input_space, f) for f in filters
            ]


    def _could_be_a_space(self, test): # nasty checking code that I don't want to redupe all the time
        if isinstance(test, (BasisStateSpace, BasisMultiStateSpace)):
            return True
        else:
            try:
                lt = len(test)
            except TypeError:
                return False
            else:
                if lt == 0:
                    return True

            try:
                if lt > 0 and isinstance(test[0], (int, np.integer)):
                    return True
                else:
                    lt = len(test[0])
                    if lt > 0 and isinstance(test[0][0], (int, np.integer)):
                        return True
                    else:
                        return False
            except (TypeError, IndexError):
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
                    len(test) == 2
                    and (test[0] is None or self._could_be_a_space(test[0]))
                    and (self._could_be_rules(test[1]))
            )
        )

    def canonicalize_prefilters(self, basis, prefilters):
        """
        Puts the prefilters in canonical form...

        :param basis:
        :type basis:
        :param prefilters:
        :type prefilters:
        :return:
        :rtype:
        """

        if prefilters is None:
            return None

        if self._could_be_a_prefilter(prefilters):
            prefilters = (prefilters,)

        final_filters = []
        for n, filter_space in enumerate(prefilters):
            if filter_space is None or self._could_be_a_space(filter_space):
                filter_rules = None
            else:
                filter_space, filter_rules = filter_space

            if filter_space is not None:
                if isinstance(filter_space, (int, np.integer)):
                    filter_space = BasisStateSpace.from_quanta(
                        self.input_space.basis,
                        (filter_space,)
                    )
                elif isinstance(filter_space, (list, tuple)) and all(
                        isinstance(f, (int, np.integer)) for f in filter_space
                ):
                    filter_space = BasisStateSpace.from_quanta(
                        self.input_space.basis,
                        filter_space
                    )
                if not isinstance(filter_space, (BasisStateSpace, BasisMultiStateSpace)):
                    filter_space = BasisStateSpace(basis, filter_space)
                if isinstance(filter_space, SelectionRuleStateSpace):
                    filter_space = filter_space.to_single().take_unique()

            final_filters.append([filter_space, filter_rules])

        return final_filters

    @classmethod
    def estimate_property_storage(cls, initial_space, target_space, perturbation_rules,
                                  property_rules, order=2, filtered=True, postfilters=None,
                                  compact=False, **opts):
        """Conservative state/adjacency RAM estimate, without generating states.

        The filtered case uses exactly the target property's rule tree. Counts
        of complete quanta manifolds and signed rule permutations are upper
        bounds; explicit filters and numerical zeros may reduce the actual work.
        """
        initial = initial_space.excitations if isinstance(initial_space, BasisStateSpace) else np.asarray(initial_space)
        targets = target_space.excitations if isinstance(target_space, BasisStateSpace) else np.asarray(target_space)
        ndim = targets.shape[-1]
        actual = collections.Counter(int(q) for q in np.sum(targets,axis=1))
        for q in np.sum(initial,axis=1):
            actual[int(q)] = max(actual[int(q)], 1)
        if filtered:
            rules = cls._create_prop_rules_tree(initial, targets, perturbation_rules,
                                                property_rules, order=order, **opts)
        else:
            reachable = [set(actual)]
            rules = {}
            for k in range(1,order+1):
                next_q = set().union(*reachable)
                for h in range(1,min(k,len(perturbation_rules))+1):
                    by_q = {}
                    for q in reachable[k-h]:
                        allowed = [r for r in perturbation_rules[h-1]
                                   if q >= -sum(x for x in r if x < 0)]
                        by_q[q] = allowed
                        next_q.update(q+sum(r) for r in allowed)
                    rules[(h,k-h)] = by_q
                reachable.append(next_q)
        all_quanta = set(actual)
        pairs = 0
        largest_step = 0
        for (h,y),by_q in rules.items():
            step = 0
            for q,changes in by_q.items():
                sources = actual.get(int(q),0) if y == 0 else math.comb(ndim+int(q)-1,int(q))
                permutations = 0
                signatures = {tuple(sorted(int(x) for x in r if x)) for r in changes}
                for signature in signatures:
                    if len(signature) > ndim:
                        continue
                    all_quanta.add(int(q+sum(signature)))
                    counts = collections.Counter(signature)
                    count = math.factorial(ndim)//math.factorial(ndim-len(signature))
                    for c in counts.values():
                        count //= math.factorial(c)
                    permutations += count
                step += sources*permutations
            largest_step = max(largest_step,step)
            pairs += step
        nstates = sum(math.comb(ndim+q-1,q) for q in all_quanta if q >= 0)
        # Occupation rows, ranks, change ids, masks, and direct-sum count/offset
        # arrays, plus canonical frontier data. This estimates state-space RAM,
        # not operator tensors, Hamiltonian values or wavefunction amplitudes.
        per_pair = ndim*np.dtype(BasisStateSpace.excitations_dtype).itemsize + 8*4 + 1
        cost = pairs*per_pair + nstates*(ndim+8)
        return dict(estimated_states=int(nstates), estimated_pairs=int(pairs),
                    largest_step_pairs=int(largest_step), estimated_bytes=int(cost),
                    quanta=sorted(int(q) for q in all_quanta if q >= 0),
                    estimate_kind='property-rule upper bound' if filtered else 'unfiltered-rule upper bound')

    @classmethod
    def from_property_rules(cls,
                            initial_space,
                            target_space,
                            perturbation_rules,
                            property_rules,
                            order=2,
                            postfilters=None
                            , compact=False
                            , **opts
                            ):
        """
        :param initial_space:
        :type initial_space:
        :param target_space:
        :type target_space:
        :param perturbation_rules:
        :type perturbation_rules:
        :param property_rules:
        :type property_rules:
        :param order:
        :type order:
        :return:
        :rtype:
        """
        int_rules = cls._create_prop_rules_tree(
            initial_space,
            target_space,
            perturbation_rules,
            property_rules,
            order=order,
            **opts
        )

        if isinstance(initial_space, BasisStateSpace):
            basis = initial_space.basis
            ndim = initial_space.ndim
        else:
            from .HarmonicOscillator import HarmonicOscillatorProductBasis
            ndim = len(initial_space[0])
            basis = HarmonicOscillatorProductBasis(ndim)

        if postfilters is None:
            postfilters = {}
            single_post = False
        else:
            single_post = not hasattr(postfilters, 'items')
            if not single_post:
                single_post = not any(
                    len(k) == 2 and isinstance(k[0], (int, np.integer)) and isinstance(k[1], (int, np.integer))
                    for k in postfilters
                )

        if compact:
            from .ExternalStateSpaces import QuantaStateSpace
        make_quanta = (lambda q: QuantaStateSpace(basis,q)) if compact else (lambda q: BasisStateSpace.from_quanta(basis,[q]))
        full_rules = {
            k:{
                'prefilters':tuple((make_quanta(i), r) for i,r in v.items()),
                'postfilters':(postfilters if single_post else (None if k not in postfilters else postfilters[k]))
            } for k,v in int_rules.items()
        }

        return full_rules

    @classmethod
    def property_wavefunction_quanta(cls, initial_space, target_space,
                                     perturbation_rules, property_rules, order=2):
        """Select WFC endpoints that may occur in a property expansion.

        This is a selection-rule screen, not a numerical VPT calculation. Reuse
        the property planner's path trees, retaining the cumulative perturbative
        order at every prefix. Prefixes include lower-order wavefunctions used
        inside higher-order corrections. Closed energy/normalization paths keep
        their endpoint manifolds too: equal total quanta does not imply equal
        states. The result maps (source quanta, WFC order) to endpoint quanta.
        No Hamiltonian matrices or correction magnitudes are constructed.
        """
        rules = list(perturbation_rules)
        rules += [()] * max(0, order-len(rules))
        _, initial = cls._build_nquanta_trees(initial_space, rules, order)
        _, target = cls._build_nquanta_trees(target_space, rules, order)
        result = {}

        def mark(final_quanta, paths):
            for path in paths:
                source_quanta = int(final_quanta-sum(d for h,d in path))
                quanta = source_quanta
                correction_order = 0
                for h,change in path:
                    quanta += int(change)
                    correction_order += int(h)
                    result.setdefault((source_quanta,correction_order),set()).add(quanta)

        for tree in (initial,target):
            for endpoints in tree:
                for q,paths in endpoints.items():
                    closed = [p for p in paths if sum(d for h,d in p)==0]
                    mark(q,closed)
        if property_rules is None:
            return result
        property_quanta = [set(sum(r) for r in changes) for changes in property_rules]
        for total in range(order+1):
            for a in range(total+1):
                for b in range(total-a+1):
                    c = total-a-b
                    if c >= len(property_quanta):
                        continue
                    for qi,pi in initial[a].items():
                        for qf,pf in target[b].items():
                            if (qf-qi in property_quanta[c]
                                    or qi-qf in property_quanta[c]):
                                mark(qi,pi)
                                mark(qf,pf)
        return result

    @classmethod
    def _build_nquanta_trees(cls, initial_space, perturbation_rules, order):
        """
        Build the tree of paths from the initial/target spaces to the
        final numbers of quanta based on the perturbation rules

        :param initial_space:
        :type initial_space:
        :param target_space:
        :type target_space:
        :param perturbation_rules:
        :type perturbation_rules:
        :param order:
        :type order:
        :return:
        :rtype:
        """

        if isinstance(initial_space, BasisStateSpace):
            initial_space = initial_space.excitations
        init_nquanta = np.unique(np.sum(initial_space, axis=1))
        pert_rules_nquanta = [np.unique([sum(x) for x in p]) for p in perturbation_rules]
        init_places = [{k: [()] for k in init_nquanta}]
        for o in range(order):
            # we determine what rules apply at each order
            # and where this can take the initial and final spaces
            r1_dict = {}
            for j,r1 in enumerate(tuple(init_places)):
                p = pert_rules_nquanta[o-j]
                for r,v in r1.items():
                    for k in p:
                        key = (1+o-j, k)
                        n = r + k
                        if n >= 0:
                            if n not in r1_dict:
                                r1_dict[n] = [t+(key,) for t in v]
                            else:
                                r1_dict[n].extend(t+(key,) for t in v)
            init_places.append(r1_dict)

        return init_nquanta, init_places

    @classmethod
    def _create_prop_rules_tree(cls,
                                initial_space,
                                target_space,
                                perturbation_rules,
                                property_rules,
                                order=2,
                                eliminate_negatives=True,
                                subspace=True,
                                check_subspaces=True
                                ):

        # we really want to do this one at a time for every initial quanta
        # and target quanta pair

        init_nquanta, init_places = cls._build_nquanta_trees(initial_space, perturbation_rules, order)
        targ_nquanta, targ_places = cls._build_nquanta_trees(target_space, perturbation_rules, order)
        # print(init_places)

        # with those sets built we can now figure out which property rules can connect
        # the terms at which orders
        rules = []
        def check_zero_rule(x):
            return (
                sum(y[1] for y in x) == 0 # doesn't change number of quanta
                # and # and the term will be iterated enough this rule can apply
            )
        for k in range(order+1):
            for p1, v1 in init_places[k].items():
                clean_v = [x for x in v1 if check_zero_rule(x)]
                if len(clean_v)>0:
                    rules.append(
                        ((p1, clean_v), None, (None, None))
                    )
            for p2, v2 in targ_places[k].items():
                clean_v = [x for x in v2 if check_zero_rule(x)]
                if len(clean_v)>0:
                    rules.append(
                        ((None, None), None, (p2, clean_v))
                    )
            if property_rules is not None:
                prop_rules_nquanta = [np.unique([sum(x) for x in p]) for p in property_rules]
                for i in range(k+1):
                    for j in range(k-i+1):
                        matches = [
                            ((p1, v1), x, (p2, v2))
                            for p1,v1 in init_places[i].items()
                            for p2,v2 in targ_places[j].items()
                            for x in prop_rules_nquanta[k-i-j]
                            if p2+x==p1 or p1+x==p2
                        ]
                        rules.extend(matches)

        # now we need to take these matches and use them to determine which terms are truly needed
        # for both the initial and target spaces
        keys = {}
        for (n,v), _, _ in rules:
            if n is not None:
                if n not in keys:
                    keys[n] = set(v)
                else:
                    keys[n].update(v)
        for _, _, (n,v) in rules:
            if n is not None:
                if n not in keys:
                    keys[n] = set(v)
                else:
                    keys[n].update(v)
        full_rules = cls._construct_rule_table_from_paths(perturbation_rules, keys)

        if eliminate_negatives:
            full_rules = {
                # we can eliminate rules that would require dropping below zero
                ord:{
                    k:[r for r in rule_vec if -sum(x for x in r if x < 0) <= k]
                    for k,rule_vec in rules.items()
                }
                for ord,rules in full_rules.items()
            }

        if isinstance(target_space, BasisStateSpace):
            target_space = target_space.excitations
        space_counts = dict(zip(*np.unique(np.sum(target_space, axis=1), return_counts=True)))
        quants = list(space_counts.keys())
        if check_subspaces:
            is_subspace = {
                k:space_counts[k] < v
                for k,v in zip(
                    quants,
                    BasisStateSpace.num_states_with_quanta(quants, len(target_space[0]))
                )
            }
        else:
            is_subspace = {}

        tree_ordering = list(sorted(full_rules.keys()))
        tree_groups = [[full_rules[x] for x in tree_ordering if x[0] == i+1] for i in range(order)]
        new_rules = {
            k:cls._prune_transition_tree(
                tree_groups[k[0]-1][:k[1]],
                t,
                is_subspace
            )
            for k,t in full_rules.items()
        }
        if property_rules is None:
            new_rules[(order, 0)] = {k:((),) for k in new_rules[(order, 0)]}
        # print(new_rules, order)

        # # Do one final prune based on rules
        # # this is necessary because when working with
        # # subspaces we can't make use of symmetry as effectively
        # # as we otherwise could
        # if property_rules is not None:
        #     prop_rules_nquanta = [np.unique([sum(x) for x in p]) for p in property_rules]
        #     for ord in range(order+1):
        #         for k in range(ord+1):
        #             for j in range(ord+1 - k): # all pairs of rules
        #                a = ord - k - j
        #                start_quanta = [
        #
        #                ]
        #                possible_qfs = [
        #                    [qi + p for qi in start_quanta for p in dqs if qi + p >= 0]
        #                    for dqs in prop_rules_nquanta
        #                ]
        #                for init_ham in range(1, k+1):
        #                    for final_ham in range()
        #                # for qi in init_nquanta:
        #                #     for qf in targ_nquanta:
        #                #         for k in
        #
        #
        #                print(k, a, j, prop_rules_nquanta[a])

        if not subspace: # TODO: infer subspace from number of states provided...
            # we're not transforming full spaces (index by quanta) and so we need to
            # include back in at higher order some of the transformations from lower order
            # that we excluded via symmetry
            for o in range(1, order+1):
                for j in range(1, order+1 - o): # all the higher order terms
                    prev = (o, j-1)
                    cur = (o, j)
                    new_terms = {}
                    max_quant = max(sum(t) + oo for oo,v in new_rules[prev].items() for t in v)
                    # print(max_quant, cur, prev)
                    for k, v in new_rules[cur].items():
                        base_rules = set(v)
                        added_terms = [
                            t for t in
                                new_rules[prev].get(k, ())
                            if 0 <= sum(t) + k and sum(t) + k <= max_quant
                        ]
                        # print(k, added_terms)
                        base_rules.update(added_terms)
                        new_terms[k] = tuple(base_rules)
                    new_rules[cur] = new_terms

        """
        {3: ((-2, 1), (-1,), (1, -1, -1)), 1: ((1, 1, -1), (1,), (2, -1)), 0: ((1,),)}
        {3: ((-2, 1), (-1,), (1, -1, -1)), 1: ((1, 1, -1), (1,), (2, -1)), 4: ((-1, -1, -1), (-3,), (-2, -1)), 0: ((1,),), 2: ((-2, 1), (-1,), (1, -1, -1))}
"""
        # del new_rules[(1, 1)][2]
        # del new_rules[(1, 1)][4]
        # new_rules[(1, 1)][1] = ((1, 1, -1), (1,), (2, -1), (3,), (2, 1))

        # raise Exception(new_rules[(1, 1)])

        return new_rules

    @classmethod
    def _prune_full_rules(self, init_space, targ_space, full_rules, property_rules, order):
        # we now check which rules truly can get us from the init space
        # to the targ space in under the target order
        raise NotImplementedError("???")
        for key, states in full_rules.items():
            ...


    @staticmethod
    def _construct_rule_table_from_paths(perturbation_rules, keys):
        init_full_rules = {}
        for t,r in keys.items():
            for k in r:
                cur_d = t - sum(d for h,d in k)  # this will track which numbers of quanta to filter on
                for n, (h, d) in enumerate(k):
                    # h is which Hamiltonian term we're working with
                    # and d is which displacement is allowed
                    key = (h, n)
                    subrules = set(s for s in perturbation_rules[h - 1] if sum(s) == d)
                    if key not in init_full_rules:
                        init_full_rules[key] = {cur_d:subrules}
                    elif cur_d not in init_full_rules[key]:
                        init_full_rules[key][cur_d] = subrules
                    else:
                        init_full_rules[key][cur_d].update(subrules)
                    cur_d += d
                # print(init_full_rules[key])

        return init_full_rules

    @staticmethod
    def _prune_transition_tree(prev_trees, ttree, subspace_specs):
        # remove terms by symmetry
        new_tree = {}
        for k,paths in ttree.items():
            # we go through the proposed changes for the current
            # tree and check if any of them have a symmetric
            # pair starting from a state with fewer quanta
            #
            # We also need to be careful about subspaces, in particular
            # if we are transforming from a subspace, we can miss terms
            # notably when a state is accessed by an intermediate
            subpaths = {}
            for v in paths:
                shift = sum(v)
                new = shift + k
                duped = False
                if new < k:
                    w = tuple(sorted((-x for x in v), key=lambda x:-100*abs(x)-x))
                    for t in prev_trees:
                        if new in t:
                            if w in t[new]:
                                if subspace_specs.get(new, False):
                                    if new not in subpaths:
                                        subpaths[new] = []
                                    subpaths[new].append(w)
                                duped = True
                                break
                if not duped:
                    if k not in subpaths:
                        subpaths[k] = []
                    subpaths[k].append(v)

            if len(subpaths) > 0:
                for k,sp in subpaths.items():
                    if k in new_tree:
                        new_tree[k] = new_tree[k] + sp
                    else:
                       new_tree[k] = sp

        return {k:tuple(v) for k,v in new_tree.items()}


    @classmethod
    def generate_nquanta_filter(cls, initials, rules, finals):
        """
        Takes the initial number of quanta, a set of possible rules, and
        a set of final numbers of quanta and determines which rules apply

        :param initial:
        :type initial:
        :param rules:
        :type rules:
        :return:
        :rtype:
        """

        return {
            k: [r for r in rules if k + sum(r)] in finals for k in initials
        }



class PropertyCorrectionFilter:
    """Order-aware structural screen for wavefunction-correction endpoints.

    Uses the same quanta-path planner as property-space generation. Closed
    energy and normalization paths are retained as well as operator paths.
    A failed test proves the coefficient unused under these selection rules;
    a passed test does not assert a nonzero numerical operator coefficient.
    """
    def __init__(self, initial_space, target_space, perturbation_rules,
                 property_rules, order=2, unfiltered=False):
        self.initial_space = initial_space
        self.target_space = target_space
        self.perturbation_rules = perturbation_rules
        self.property_rules = property_rules
        self.order = order
        self.unfiltered = unfiltered
        self._uses = None

    @classmethod
    def from_property(cls, target_space, initial_space=None, order=2,
                      target_property='intensities', property_rules=None,
                      perturbation_rules=None):
        if target_property not in {'intensities', 'frequencies', 'wavefunctions'}:
            raise ValueError('unknown correction-filter target: '+str(target_property))
        basis = target_space.basis
        if initial_space is None:
            initial_space = np.zeros((1, target_space.ndim), dtype=int)
        if perturbation_rules is None:
            perturbation_rules = [basis.selection_rules(*['x']*(h+2))
                                  for h in range(1, order+1)]
        if property_rules is None and target_property == 'intensities':
            property_rules = [basis.selection_rules(*['x']*(c+1))
                              for c in range(order+1)]
        return cls(initial_space, target_space, perturbation_rules,
                   property_rules, order=order,
                   unfiltered=target_property == 'wavefunctions')

    @staticmethod
    def operator_rules(basis, operators, order=2, include_intensities=True):
        """Union rules for Taylor operators O_c of degree c and dipole mu_c.

        Mirrors classic operator expansion ordering. Cartesian transformations
        preserve derivative degree. Representation matrices have no polynomial
        degree metadata: return None so callers can conservatively retain all.
        No numerical coefficient is inspected or pruned.
        """
        if not isinstance(operators, dict) or 'expansion' in operators:
            operators = {'Operator': operators}
        rules = [set(basis.selection_rules(*['x']*(c+1)))
                 if include_intensities else set() for c in range(order+1)]
        for expansion in operators.values():
            if isinstance(expansion, dict):
                expansion = expansion['expansion']
            if not isinstance(expansion, (list, tuple)):
                return None  # Unknown polynomial support: conservatively retain all.
            for c in range(min(len(expansion), order+1)):
                rules[c].update([()] if c == 0 else basis.selection_rules(*['x']*c))
        return [tuple(sorted(r)) for r in rules]

    @property
    def uses(self):
        if self._uses is None:
            self._uses = BasisStateSpaceFilter.property_wavefunction_quanta(
                self.initial_space, self.target_space, self.perturbation_rules,
                self.property_rules, order=self.order)
        return self._uses

    def __call__(self, source, endpoints, order):
        if self.unfiltered:
            return endpoints
        if not 1 <= order <= self.order:
            raise ValueError('correction order is outside the planned property order')
        if len(source) != 1:
            raise ValueError('correction filtering requires one source state')
        quanta = int(np.sum(source.excitations[0]))
        allowed = self.uses.get((quanta, order), ())
        mask = np.isin(np.sum(endpoints.excitations, axis=-1), list(allowed))
        return endpoints.take_subspace(np.flatnonzero(mask))
