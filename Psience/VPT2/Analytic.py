import numpy as np, scipy.signal

from McUtils.Zachary import DensePolynomial, TensorCoefficientPoly
from McUtils.Combinatorics import SymmetricGroupGenerator, IntegerPartitioner
from ..BasisReps import HarmonicOscillatorMatrixGenerator, HarmonicOscillatorRaisingLoweringPolyTerms

__all__ = [
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
        for i,(g,(a,b)) in enumerate(zip(generators, path)):
            d = a - b
            if g is not None:
                poly_contrib = g.get_poly_coeffs(d, shift=k)
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
                shifts[i] = k + d

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
            HarmonicOscillatorMatrixGenerator.get_paths(
                [g.size for g in gl],
                delta
            )
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
        shifts = [
            [p[1] for p in pi]
            for pi in poly_shifts
        ]
        # full set of polynomial terms are now direct product of these 1D terms...
        return self.DirectProductContribPoly(polys, shifts)

class AnalyticPTCorrectionGenerator:
    """
    Takes strings of `x`/`p` terms plus tensor identities and uses them to
    generate a correction
    """

    def __init__(self, terms):
        self.terms = terms
        self.sizes = [len(t) for t in terms]
        self.total_terms = sum(self.sizes)

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

        changes = np.abs(changes)
        total_change = np.sum(changes)
        remainder_terms = self.total_terms - total_change
        if remainder_terms % 2 != 0:
            raise ValueError("need even remainder by symmetry")
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

        return base_tree[symmetric_trees,]

    def get_correction(self, changes):
            # Step 1: determine the distinct ways to split the terms up
            #         that can generate the change
            # Step 2: create a set of matrix generators for each of these sets
            # Step 3: set up the 1D poly/energy change lists that will eventually
            #         be direct product-ed
            # Step 4: return these + the associated tensor coefficient product to be
            #         collated into a single sparse polynomial
            operator_sets = self.enumerate_operator_partitions(changes)

            raise Exception(operator_sets)





