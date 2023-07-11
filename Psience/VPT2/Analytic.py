import itertools

import numpy as np, scipy.signal

from McUtils.Zachary import DensePolynomial, TensorCoefficientPoly, PureMonicPolynomial
from McUtils.Combinatorics import SymmetricGroupGenerator, IntegerPartitioner, UniquePartitions, IntegerPartitionPermutations
from ..BasisReps import HarmonicOscillatorMatrixGenerator, HarmonicOscillatorRaisingLoweringPolyTerms

__all__ = [
    'AnalyticPerturbationTheoryDriver',
    'AnalyticPTCorrectionGenerator',
    'RaisingLoweringClasses'
]

class RaisingLoweringClasses:
    """
    A tree providing the classes of raising and lowering operations to get
    to a given set of changes (for multidimensional cases) over a given number of terms
    """
    def __init__(self, nterms, changes, check=True):
        self.nterms = nterms
        self.changes = changes
        if check:
            t_change = sum(abs(x) for x in changes)
            self.empty = t_change > nterms or (t_change % 2) != (nterms % 2)
        else:
            self.empty = False

    @classmethod
    def _enum_changes_rec(cls, nterms, changes, change_sums): #TODO: Optimize this
        if len(changes) == 1:
            k = changes[0]
            s = nterms
            yield [((k + s) // 2, (k - s) // 2)]
        else:
            # we need to take at least as many terms as are required for the first change
            # and we can take _up to_ the total number of terms minus the number required
            # to make all of the remaining changes
            k = changes[0]
            for s in range(k, nterms - change_sums[0]+1):
                for t in cls._enum_changes_rec(nterms - s, changes[1:], change_sums[1:]):
                    yield [((k + s) // 2, (k - s) // 2)] + t

    def __iter__(self):
        if not self.empty:
            change_sums = np.flip(np.cumsum(np.flip(self.changes)))[1:]
            for t in self._enum_changes_rec(self.nterms, self.changes, change_sums): yield t

class AnalyticPTPolynomialGenerators:
    """
    Provides a generic wrapper for an operator term
    given the operator identities (strings of p/q terms)
    and the breakpoint for switching from `Pi_n` to `Pi_m`

    Makes use of `HarmonicMatrixGenerator` to get polynomial coefficients in the
    different dimensions
    """

    def __init__(self, matrix_generators):
        self.generator_lists = matrix_generators # a 2D array of generators, as dimension x steps

    @classmethod
    def get_1D_path_contrib(cls, generators, path):
        """
        For a given set of energy changes in 1D provides the polynomial
        contribution and associated states


        :param generators:
        :param path:
        :return:
        """
        shifts = [0]*len(path)
        k = 0
        poly = None
        for i,(generator_list,(a,b)) in enumerate(zip(generators, path)):
            d = a - b
            poly_contribs = [
                g.poly_coeffs(d, shift=k)
                for g in generator_list
                if g is not None
            ]
            if len(poly_contribs) > 0:
                poly_lens = [len(x) for x in poly_contribs if isinstance(x, np.ndarray)]
                max_contrib_len = max(poly_lens) if len(poly_lens) > 0 else 1
                poly_contrib = sum(
                    np.pad(pc, [0, max_contrib_len - len(pc)])
                        if isinstance(pc, np.ndarray) else
                    np.zeros(max_contrib_len)
                    for pc in poly_contribs
                )

                if d != 0:
                    dir_change = np.sign(d) == -np.sign(k)  # avoid the zero case
                    if dir_change:
                        sqrt_contrib = HarmonicOscillatorRaisingLoweringPolyTerms.get_sqrt_remainder_coeffs(d, k)
                        # print(f" s({a},{b})({k})>", sqrt_contrib)
                        poly_contrib = scipy.signal.convolve(poly_contrib, sqrt_contrib)
                if poly is None:
                    poly = poly_contrib
                else:
                    poly = scipy.signal.convolve(poly, poly_contrib)
                k = k + d
                shifts[i] = k

        return poly, shifts

    def get_correction_poly_generator(self, changes):
        """
        Provides the polynomial corresponding to changing quanta
        by the proscribed integer partition

        :param changes:
        :type changes:
        :return:
        :rtype:
        """

        paths_1D = [
            [
                list(reversed(path)) for path in
                HarmonicOscillatorMatrixGenerator.get_paths(
                    [
                        # b.c. of way we store terms tedious to get lens
                        ([ggg.size for ggg in g if ggg is not None] + [0])[0]
                        for g in gl
                    ],
                    delta
                )
            ]
            for gl,delta in zip(self.generator_lists, changes)
        ]

        poly_shifts = [
            [self.get_1D_path_contrib(gl, p) for p in paths]
            for gl,paths in zip(self.generator_lists, paths_1D)
        ]
        polys = [
            [p[0] for p in pi]
            for pi in poly_shifts
        ]
        shifts = tuple(
            tuple(tuple(p[1]) for p in pi)
            for pi in poly_shifts
        )

        # full set of polynomial terms are now direct product of these 1D terms...
        return AnalyticPTProductContrib(polys, shifts, changes)

class AnalyticPTProductContrib:
    def __init__(self, polys_1D, shifts_1D, changes):
        self.polys = polys_1D
        self.shifts = shifts_1D
        self.changes = changes
        self.ndim = len(changes)

    def __repr__(self):
        return "{}({})".format(type(self).__name__, self.polys)

    def __radd__(self, other):
        return self.__add__(other)
    def __add__(self, other):
        if isinstance(other, AnalyticPTProductContrib):
            return AnalyticPTProductContribSum([self, other])
        elif isinstance(other, AnalyticPTProductContribSum):
            return other + self
        elif isinstance(other, (int, float, np.integer, np.floating)) and other == 0:
            return self
        else:
            raise TypeError("can't add {} and {}".format(self, other))


    def _substitute_coeffs(self, tensor_poly, coeff_tensors):
        contribs = [
            scaling * np.prod([coeff_tensors[k[0]][k[1]][k[2:]] for k in index_tuple_tuple])
            for index_tuple_tuple, scaling in tensor_poly.terms.items()
        ]
        # raise Exception(tensor_poly, tensor_poly.scaling*np.sum(contribs))
        return tensor_poly.scaling*np.sum(contribs)

    def get_numerators(self, coefficient_tensors):
        return [
            [
                np.sum([# I should probably have added these up earlier...
                    self._substitute_coeffs(tcp, coefficient_tensors)
                        if isinstance(tcp, TensorCoefficientPoly) else
                    tcp
                    for tcp in poly_set
                ])
                    if poly_set is not None else
                np.array([1])
                for poly_set in polys_1D
            ]
            for polys_1D in self.polys
        ]

    def _eval_shift(self, path, change, crossover):
        # we just go up normally from the left (starting at 0)
        # but from the right we need to effectively start from the final state
        num_right = len(path) - crossover
        num_left = crossover
        full_path = path - np.pad(np.full(num_right, change), [num_left, 0])
        if num_right > 0:
            # means we don't end on a Pi_n term so we drop the final element in the path
            # (based on the way we're stepping)
            full_path = full_path[:-1]
        if num_left == 0:
            # means we _only_ have Pi_m terms and so we need to
            # add on the fact that we end at -change (from the context of the end state)
            full_path = np.concatenate([[-change], full_path])
        #TODO: what do I do with terms that have a (Pi_n)(Pi_m) in them???
        return full_path

    def get_quanta_shifts(self, crossover):
        return [
            np.array([
                self._eval_shift(path, delta, crossover)
                for path in np.array(paths_1D)
            ])
            for delta,paths_1D in zip(self.changes, self.shifts)
        ]
    def get_denominators(self, freqs, crossover):
        shifts = self.get_quanta_shifts(crossover)
        return [
            [
                np.array([
                    freq * path
                    # self._eval_energy(freq, path)
                    for path in path_set
                ])
                for path_set in paths_1D
            ]
            for freq,paths_1D in zip(freqs, shifts)
        ]

    class direct_product_evaluator:
        def __init__(self, polynomial_coeffs, paths, freqs):
            self.pcs = polynomial_coeffs
            self.paths = paths
            self.freqs = freqs
        def evaluate(self, initial_states, path_exponents, crossover):
            # take direct product of numerators
            # and divide by direct sum product of denominators

            raise NotImplementedError("need to eval direct product")

    def get_poly_evaluator(self, coefficient_tensors, freqs):
        numerators = self.get_numerators(coefficient_tensors)

        return self.direct_product_evaluator(numerators, self.shifts, freqs)

        # need to turn changes + freqs + crossover into energy denominators
        # and need to substitute coefficient tensors into polynomials

    # def eval(self, initial_states, coefficient_tensors, freqs, crossover):
    #     """
    #
    #     :param initial_states: sets of arrays for initial quanta for each
    #     :return:
    #     """
    #     base_eval = self.get_poly_evaluator(coefficient_tensors, freqs, crossover)
    #     raise Exception(...)
    #     for n, poly in zip(initial_states, self.polys):
    #         ...

    def expand(self):
        """
        Expands the 1D terms as a direct product polynomial
        :return:
        """
        raise NotImplementedError("whoops")

class AnalyticPTProductContribSum:
    def __init__(self, contribs):
        self.contribs = tuple(contribs)
    def __repr__(self):
        return "{}({})".format(type(self).__name__, self.contribs)

    def __radd__(self, other):
        return self.__add__(other)
    def __add__(self, other):
        if isinstance(other, AnalyticPTProductContribSum):
            return type(self)(self.contribs+other.contribs)
        elif isinstance(other, AnalyticPTProductContrib):
            return type(self.contribs+(other,))
        elif isinstance(other, (int, float, np.integer, np.floating)) and other == 0:
            return self
        else:
            raise TypeError("can't add {} and {}".format(self, other))

    class direct_product_evaluator_sum:
        def __init__(self, evaluators):
            self.evaluators = evaluators


    def get_poly_evaluator(self, coefficient_tensors, freqs):
        return self.direct_product_evaluator_sum([
            contrib.get_poly_evaluator(coefficient_tensors, freqs)
            for contrib in self.contribs
        ])

    def eval(self, quanta, coefficient_tensors, freqs, crossover):
        ...

class AnalyticPTCorrectionSum(PureMonicPolynomial):
    """
    Represents a sum of polynomial corrections over energy path products
    """

    @classmethod
    def canonical_key(cls, path):
        return path
        # # we need a way to sort indices, which we do by grouping by key length,
        # # doing a standard sort for each length, and reconcatenating
        # s_groups = {}
        # for index_tuple in monomial_tuple:
        #     l = len(index_tuple)
        #     grp = s_groups.get(l, [])
        #     s_groups[l] = grp
        #     grp.append(index_tuple)
        # t = tuple(
        #     grp
        #     for l in sorted(s_groups.keys())
        #     for grp in sorted(s_groups[l])
        # )
        # return t

    def evaluate(self, quanta, coefficient_tensors, freqs, crossover):
        for path, contrib_gen in self.terms.items():
            contrib = contrib_gen.eval(quanta, coefficient_tensors, freqs, crossover)
        ...


class PTTensorCoeffProductSum(TensorCoefficientPoly):
    """
    Only here for the better formatting
    """
    def __repr__(self):
        sums = []
        for k_prod,v in self.terms.items():
            ks = [
                "{}[{}]{}".format(
                    ["V", "G", "V'"][k[1]],
                    k[0],
                    k[2:]
                )
                for k in k_prod
            ]
            v = round(v, 5)
            sums.append("{}{}".format(v, "".join(ks)))
        return "+".join(sums)
class AnalyticPTCorrectionGenerator:
    """
    Takes strings of `x`/`p` terms plus tensor identities and uses them to
    generate a correction
    """

    def __init__(self, term_lists, operator_indices=None):
        if operator_indices is None and isinstance(term_lists[0][0], dict): # operator inds encoded in terms
            operator_indices = [
                [next(iter(t.keys())) for t in tt]
                for tt in term_lists
            ]
            term_lists = [
                [next(iter(t.values())) for t in tt]
                for tt in term_lists
            ]
        self.term_lists = term_lists
        self.sizes = [len(t[0]) for t in term_lists]
        self.total_terms = sum(self.sizes)
        if operator_indices is None:
            operator_indices = [
                [(i, j) for j in range(len(t))]
                for i,t in enumerate(self.term_lists)
            ]
        self.operator_indices = [
            [tuple(t) for t in tt]
            for tt in operator_indices
        ]

    def _build_operator_split_trees(self, boxes, steps, generator, cache):
        # TODO: huge amount of waste that could be pruned here...
        #       by making use of the tree structure of the problem
        #       and turning the sets of boxes into a trie so invalid combos
        #       aren't tried repeatedly, also we have a bunch of permutation symmetry
        #       we aren't making use of that could reduce work across the entrie calculation
        s = steps[0]
        if s not in cache:
            cache[s] = generator.get_terms(s) # could probably be faster...

        parts = cache[s]
        new_boxes = boxes[np.newaxis] - parts # this could be done smarter I think
        good_places = np.all(new_boxes >= 0, axis=1)

        if len(steps) > 1:
            new_trees = []
            for subbox, part in zip(new_boxes[good_places,], parts[good_places,]):
                subtree = self._build_operator_split_trees(subbox, steps[1:], generator, cache)
                new_trees.append(
                    np.concatenate(
                        [
                            np.broadcast_to(part[np.newaxis, np.newaxis], (subtree.shape[0], 1, subtree.shape[2])),
                            subtree  # num_terms x num_steps x num_terms
                        ],
                        axis=1
                    )
                )
            tree = np.concatenate(new_trees, axis=0)
        else:
            tree = parts[good_places,][:, np.newaxis, :]

        return tree

    def enumerate_operator_partitions(self, changes):
        # we have a set of operators with sizes {s_i}
        # which must be distributed over D dimensions where
        # D depends on how many terms are left over after subtracting
        # off the requisite terms to cause `changes`
        #
        # To do this we enumerate the partitions of the total number of steps
        # that distribute enough terms into each dimension and then for each of
        # these partitions we find the sets of partitions of the sizes that lead
        # to a valid solutions

        changes = np.asanyarray(changes, dtype=int)
        changes = np.abs(changes).astype(int)
        total_change = np.sum(changes)
        remainder_terms = self.total_terms - total_change
        if remainder_terms % 2 != 0:
            return [], 0
            # raise ValueError("need even remainder by symmetry")
        total_dim = len([c for c in changes if c != 0]) + remainder_terms // 2

        changes = np.pad(changes, [0, total_dim - len(changes)])

        targets = IntegerPartitioner.partitions(self.total_terms, max_len=total_dim, pad=True)
        check_changes = targets - changes[np.newaxis, :]
        good_targets = np.all(check_changes >= 0, axis=1)
        targets = targets[good_targets]

        generator = SymmetricGroupGenerator(total_dim)
        cache = {}

        base_tree = np.concatenate(
            [
                self._build_operator_split_trees(t, self.sizes, generator, cache)
                for t in targets
            ], axis=0
        )

        even_check = np.sum(base_tree, axis=1) - changes[np.newaxis, :]
        symmetric_trees = np.all(even_check % 2 == 0, axis=1)

        return base_tree[symmetric_trees,], total_dim

    class PTMatrixGenerator:
        def __init__(self, generator_sets, prefactor_sets):
            self.generator_lists = generator_sets
            self.prefactor_lists = prefactor_sets
            self.size = len(generator_sets[0].terms)
            # self.terms = generator_sets[0].terms
        def poly_coeffs(self, change, shift=0):
            base_lists = []
            for gen, prefactor in zip(self.generator_lists, self.prefactor_lists):
                base = gen.poly_coeffs(change, shift=shift)
                if not (isinstance(base, (int, float, np.integer, np.floating)) and base == 0):
                    if not isinstance(prefactor, (int, float, np.integer, np.floating)):
                        base = np.array([prefactor*b if b != 0 else 0 for b in base], dtype=object)
                    else:
                        base = prefactor * base
                base_lists.append(base)
            return sum(base_lists)
        def __repr__(self):
            return "{}({})".format(type(self).__name__, [g.terms for g in self.generator_lists])

    def _build_matrix_generator(self, step, index_list, cache, prefactor):
        if len(index_list) == 0:
            return None

        # print(index_list, self.term_lists)
        gens = []
        for terms in self.term_lists[step]:
            k = tuple([terms[i] for i in index_list])
            if k not in cache:
                # if k[1] > 2:
                #     raise ValueError("didn't expect to have more than one p?")
                # elif k[1] == 2:
                #     term_list = ['x'] * k[0] + ['p']
                # elif k[1] == 1:
                #     term_list = ['x'] * k[0] + ['p']
                # else:
                #     term_list = ['x'] * k[0]
                cache[k] = HarmonicOscillatorMatrixGenerator(k)
            gens.append(cache[k])
        mg = self.PTMatrixGenerator(gens, prefactor)
        return mg

    def _build_matgen_prefactor_list(self, operator_splits, operator_index, op_cache):
        def _prefactor(op_idx, num_subops, inv):
            base_list = [i for i,v in enumerate(inv) for _ in range(len(v))]
            # print("prefactor:", (k,) + tuple(base_list[j] for j in np.concatenate(inv)), inv)
            base_tup = tuple(base_list[j] for j in np.concatenate(inv))
            return [
                PTTensorCoeffProductSum.monomial(self.operator_indices[op_idx][subop_idx] + base_tup)
                for subop_idx in range(num_subops)
            ]

        index_sets = operator_splits[0]
        inverse_sets = operator_splits[1]
        prefactors = [
            _prefactor(operator_index, len(self.term_lists[operator_index]), i)
            for i in zip(*inverse_sets)
        ] # orderd by partition

        # print(">>", operator_index, index_sets)
        return [
            [
                self._build_matrix_generator(operator_index, p, op_cache, f)
                for p,f in zip(pset, prefactors)
            ]
            for pset in index_sets
        ]

    def build_matrix_generators(self, operator_sets):


        # subbies = [
        #     UniqueSubsets(t)
        #     for t in term_encodings
        # ]
        op_cache = {}
        parts = [UniquePartitions(np.arange(len(t[0]))) for t in self.term_lists]
        op_subs = [
            [
                # [
                #     self._build_matrix_generator(subby.get_subsets(s), op_cache)
                #     # this actually provides an over-estimate, instead we need to take all unique disjoint subsets
                #     # of
                #     for subby, s in zip(subbies, dim_subs)
                # ]
                # for dim_subs in op_list
                p.partitions(ol,
                             take_unique=False,
                             return_partitions=False,
                             return_indices=True,
                             return_inverse=True
                             ) # I actually don't want to do the unique terms...? (loses scaling)
                for p, ol in zip(parts, op_list)
            ]
            for op_list in operator_sets
        ]
        # raise Exception(
        #     parts[1].partitions(
        #         operator_sets[1][1],
        #         take_unique=False,
        #         return_partitions=False,
        #         return_indices=True,
        #         return_inverse=True
        #     )[0]
        # )

        # now convert these encoded lists into proper matrix generators
        mat_gens = [
            [
                self._build_matgen_prefactor_list(term_partition_lists, operator_index, op_cache)
                for operator_index, term_partition_lists in enumerate(sub_list) # iterate across operators themselves
            ]
            for sub_list in op_subs # iterate over ways to divide operators across coords
        ]

        return mat_gens

    def get_correction(self, changes):
            # Step 1: determine the distinct ways to split the terms up
            #         that can generate the change
            # Step 2: create a set of matrix generators for each of these sets
            # Step 3: set up the 1D poly/energy change lists that will eventually
            #         be direct product-ed
            # Step 4: return these + the associated tensor coefficient product to be
            #         collated into a single sparse polynomial

            if self.total_terms == 0: # just a product of tensor coefficients assuming orthogonality applies
                cs = np.sum(np.abs(changes))
                if cs > 0:
                    return 0
                else:
                    polys = [
                        [
                            TensorCoefficientPoly.monomial(self.operator_indices[i][j])
                            for j,_ in enumerate(t)
                        ]
                        for i, t in enumerate(self.term_lists)
                    ]
                    shifts = tuple(
                        (0,)*len(t)
                        for t in self.term_lists
                    )
                    return AnalyticPTProductContrib(polys, shifts, [])

            operator_sets, dimension = self.enumerate_operator_partitions(changes)
            if len(operator_sets) == 0: # _not_ an empty operator, just no valid paths
                return 0
            changes = np.pad(changes, [0, dimension - len(changes)]).astype(int)
            mat_gens = self.build_matrix_generators(operator_sets)

            polys = [
                AnalyticPTPolynomialGenerators(
                    [
                        [tm[i] for tm in term_mat_list]
                        for i in range(dimension)
                    ]
                ).get_correction_poly_generator(changes)
                for term_mat_list in mat_gens
            ]
            # for p,gl in zip(polys, mat_gens):
            #     for i,(pl,sh) in enumerate(zip(p.polys, p.shifts)):
            #         print("Generators({i})".format(i=i), [g[i] for g in reversed(gl)])
            #         print("Shifts({i})".format(i=i), sh)
            #         print("Poly-Lens({i})".format(i=i), [len(py) if py is not None else 0 for py in pl])
            # raise Exception(...)

            return AnalyticPTCorrectionSum({
                g.shifts:g for g in polys
            })

class AnalyticPerturbationTheoryDriver:
    """

    """
    def __init__(self,
                 hamiltonian_terms
                 ):
        self.terms = hamiltonian_terms
        self._cache = {}

    @classmethod
    def from_order(cls, order, internals=True):
        if internals:
            hamiltonian_terms = [
                [
                    [
                        {(o,0):['x']*(o+2)},
                        {(o,1):['p'] + ['x']*o + ['p']}
                    ]
                ] + (
                    []
                        if o < 2 else
                    [
                        [
                            {(o,2):['x'] * (o-2)}
                        ] # V' term has different shape from rest
                    ]
                )
                for o in range(order+1)
            ]
        else:
            raise NotImplementedError("Cartesians coming soon?")

        return cls(hamiltonian_terms)

    def _take_term_direct_prod(self, terms):
        prod_inds = itertools.product(*[range(len(t)) for t in terms])
        return [
            [t[i] for t,i in zip(terms,p)]
            for p in prod_inds
        ]
    def energy_terms(self, order, return_indices=False):
        term_indices= set(itertools.chain(*( # few enough I can be a bit inefficient
                list(itertools.permutations(part))
                for order_parts in IntegerPartitioner.partitions(order, pad=False)
                for part in order_parts
        )))

        terms = [
            self._take_term_direct_prod([self.terms[i] for i in term_index_list])
            for term_index_list in term_indices
        ]

        # raise Exception(terms)

        if return_indices:
            return terms, term_indices
        else:
            return terms

    def energy_correction_driver(self, order):
        return AnalyticPTEnergyDriver(
            [
                self.energy_terms(o, return_indices=True)
                for o in range(1, order+1)
            ]
        )

class AnalyticPTOverlapDriver:
    """
    Provides a driver to get the wfn contribution for a state from the state itself
    at any given order (needed assuming we don't use intermediate normalization)
    """

class AnalyticPTEnergyDriver:
    def __init__(self, terms_inds):
        self.term_inds = terms_inds
        self.order = len(terms_inds)
        self._base_terms = None

    @property
    def generic_corrections(self):
        if self._base_terms is None:
            self._base_terms = self._evaluate_base_corrs()
        return self._base_terms
    def _evaluate_base_corrs(self): # this could be shared across orders...
        term_map = {}
        for order_list, order_inds in self.term_inds:
            for terms,inds in zip(order_list, order_inds):
                term_map[inds] = sum(AnalyticPTCorrectionGenerator(t).get_correction([]) for t in terms)
        return term_map

    def _enumerate_ov_e_pi_prods(self, index_tuple):
        #TODO: figure out parities

        remainder = self.order - sum(index_tuple)
        nterms = len(index_tuple)
        if remainder == 0: # just put a Pi between every term
            return [[]], [[]], [[0] + [1]*(nterms - 1) + [0]]
        elif len(index_tuple) == 1: # only one term, so can only multiply by the overlap
            return [[remainder]], [[]], [[0, 0]]
        else: # need to figure out how to split up terms across the energies and overlaps
            raise NotImplementedError("still working out how to split integer partition permutations over E and Ov")


    def get_poly_evaluator(self, tensor_coefficients, freqs):
        poly_terms = self.generic_corrections
        for index_tuple, poly_sum in poly_terms:
            # evaluate sets of overlap * energy_correction * Pi's to weight polynomials
            overlap, ecorrs, exponents = self._enumerate_ov_e_pi_prods(index_tuple)










