import abc
import itertools

import numpy as np, scipy.signal

from McUtils.Zachary import DensePolynomial, TensorCoefficientPoly, PureMonicPolynomial
from McUtils.Combinatorics import SymmetricGroupGenerator, IntegerPartitioner, UniquePartitions, IntegerPartitionPermutations
from ..BasisReps import HarmonicOscillatorMatrixGenerator, HarmonicOscillatorRaisingLoweringPolyTerms

__all__ = [
    'AnalyticPerturbationTheorySolver',
    'AnalyticPerturbationTheoryDriver',
    'AnalyticPTCorrectionGenerator',
    'RaisingLoweringClasses'
]

class AnalyticPerturbationTheorySolver:
    """
    A re-attempt at using the recursive expressions
    to provide simpler code for getting APT expressions
    """
    def __init__(self, hamiltonian_expansion):
        self.hamiltonian_expansion = hamiltonian_expansion

    @classmethod
    def from_order(cls, order, internals=True):
        if order < 2:
            raise ValueError("why")
        if internals:
            hamiltonian_terms = [
                HamiltonianExpansionTerm(
                    [
                        ['x'] * (o + 2),
                        ['p'] + ['x'] * o + ['p']
                    ] + (
                        []
                            if o < 2 else
                        [
                            ['x'] * (o - 2)
                        ]  # V' term has different shape from rest
                    ),
                    order=o,
                    identities=[0, 1] + ([] if o < 2 else [2]) # V, G, G'
                )
                for o in range(order + 1)
            ]
        else:
            hamiltonian_terms = [
                HamiltonianExpansionTerm(
                    [
                        ['x'] * (o + 2)
                    ] + (
                        [['p', 'p']] if o == 0 else []
                    )+ (
                        []
                        if o < 2 else
                        [
                            ['x', 'p'] + ['x'] * (o-1) + ['p'], # Coriolis
                            ['x'] * (o - 2) # Watson
                        ]
                    ),
                    order=o,
                    identities=([0, 1] if o == 0 else [0]) + ([] if o < 2 else [3, 4])  # V, G, G'
                )
                for o in range(order + 1)
            ]

        return cls(hamiltonian_terms)

    def energy_correction(self, order):
        return EnergyCorrection(self, order)

    def wavefunction_correction(self, order):
        return WavefunctionCorrection(self, order)

    def overlap_correction(self, order):
        return WavefunctionOverlapCorrection(self, order)

class ProductPTPolynomial:
    """
    TODO: include prefactor term so we can divide out energy changes
    """
    def __init__(self, coeffs, prefactor=1):
        if (
                isinstance(coeffs, (int, float, np.integer, np.floating))
                or isinstance(coeffs[0], (int, float, np.integer, np.floating))
        ):
            raise ValueError("coeffs must be a vector of vectors (not {})".format(coeffs))
        self.coeffs = [np.asanyarray(c) for c in coeffs] # coeffs along each dim independently
        self.prefactor = prefactor
        self._order = None

    @property
    def order(self):
        if self._order is None:
            self._order = tuple(len(c)-1 for c in self.coeffs)
        return self._order
    def __repr__(self):
        return "{}(<{}>)".format(type(self).__name__, ",".join(str(c) for c in self.order))

    @staticmethod
    def _monify(coeffs):
        acoeffs = np.abs(coeffs)
        mcoeffs = np.max(acoeffs)
        if mcoeffs == 0:
            raise ValueError(coeffs)
        min_pos = np.argmin(acoeffs + (1+mcoeffs) * (acoeffs == 0))
        min_v = coeffs[min_pos]
        if min_v == 0:
            raise ValueError(coeffs, min_pos, acoeffs + (1+mcoeffs) * (acoeffs == 0))
        return coeffs / min_v, min_v
    def combine(self, other:'ProductPTPolynomial'):
        off_pos = None # if a _single_ mode differs we can still combine
        # first pass off lengths b.c. most can't be simplified...
        for n,(c1,c2) in enumerate(zip(self.coeffs, other.coeffs)):
            if len(c1) != len(c2):
                if off_pos is not None:
                    return False
                off_pos = n
        # second pass where we actually compare the coeffs
        new_coeffs = []
        for n,(c1,c2) in enumerate(zip(self.coeffs, other.coeffs)):
            if off_pos is not None and off_pos == n: # can handle _one_ difference
                if len(c2) > len(c1):
                    c1 = np.pad(c1, [0, len(c2)-len(c1)])
                if len(c1) > len(c2):
                    c2 = np.pad(c2, [0, len(c1)-len(c2)])
                new_coeffs.append(c1+c2)
                continue

            monic_coeffs1, monic_scaling1 = self._monify(c1)
            monic_coeffs2, monic_scaling2 = self._monify(c2)
            if np.any(monic_coeffs1 != monic_coeffs2):
                if off_pos is None:  # can handle _one_ difference
                    off_pos = n
                    new_coeffs.append(c1 + c2)
                    continue
                else:
                    return False

            scaling = self.prefactor*monic_scaling1 + other.prefactor*monic_scaling2
            if scaling == 0:
                return True # special case... where they cancel perfectly

            new_coeffs.append(monic_coeffs1 * scaling)

        return ProductPTPolynomial(new_coeffs)

    def mul_simple(self, other:'ProductPTPolynomial'):
        ocs = other.coeffs
        scs = self.coeffs
        if len(ocs) < len(scs):
            ocs = ocs + [[1]]*(len(scs) - len(ocs))
        elif len(scs) < len(ocs):
            scs = scs + [[1]]*(len(scs) - len(scs))

            # raise ValueError("not sure how to 'simply multiply' {} and {}".format(self, other))
        return type(self)(
            [
                scipy.signal.convolve(sc, oc)
                for sc, oc in zip(scs, ocs)
            ],
            prefactor=self.prefactor * other.prefactor
        )
    def mul_along(self, other:'ProductPTPolynomial', inds, remainder=None):
        if isinstance(other, ProductPTPolynomialSum):
            return type(other)(
                [
                    self.mul_along(o, inds, remainder=remainder)
                    for o in other.polys
                ],
                prefactor=self.prefactor * other.prefactor
            )

        if len(inds) == 0:
            # we get to just concatenate the coeffs
            return type(self)(self.coeffs+other.coeffs, prefactor=self.prefactor*other.prefactor)

        if isinstance(inds[0], (int, np.integer)):
            # we assume the inds line up with the inds for the other coeffs
            inds = [
                np.arange(len(inds)),
                inds
            ]

            # remainder is just an optimization since we get it for free
            if remainder is None:
                remainder = np.setdiff1d(np.arange(len(self.coeffs)), inds[1])
            if len(remainder) == 0 or isinstance(remainder[0], (int, np.integer)):
                remainder = [
                    np.arange(len(inds[1]), len(other.coeffs)),
                    remainder
                    ]

        other_inds, self_inds = inds # axes to convolve
        if remainder is None:
            remainder = [
                np.setdiff1d(np.arange(len(other.coeffs)), other_inds),
                np.setdiff1d(np.arange(len(self.coeffs)), self_inds)
            ]
        other_remainder, self_remainder = remainder

        # self_remapping = np.argsort(np.concatenate([self_inds, self_remainder])) # literally just the inverse...(I could pass this in)
        # other_remapping = np.argsort(np.concatenate([self_inds, self_remainder])) # should avoid actually doing this sort (when `remainder` is/was `None)

        # print(other_inds, self_inds, len(other.coeffs), len(self.coeffs))
        new_coeffs = [
            scipy.signal.convolve(other.coeffs[ox], self.coeffs[sx])
            #     self._remap_tensor_coefficients(other.coeffs[ox], other_remapping),
            #     self._remap_tensor_coefficients(self.coeffs[sx], self_remapping)
            # )
            for ox, sx in zip(other_inds, self_inds)
        ] + [
            self.coeffs[sx] for sx in self_remainder
        ] + [
            other.coeffs[ox] for ox in other_remainder
        ]

        return type(self)(new_coeffs, prefactor=self.prefactor*other.prefactor)

    # def __mul__(self, other):
    #     if isinstance(other, ProductPTPolynomial):
    #         return self.mul_simple(other)
    #     elif isinstance(other, (ProductPTPolynomialSum, PTTensorCoeffProductSum)):
    #         return other.mul_simple(self)
    #     else:
    #         return type(self)(self.coeffs, prefactor=self.prefactor*other)

    # @classmethod
    # def _remap_tcp_tuple(cls, tcp_tuple, perm):
    #     return (
    #             tcp_tuple[:2] # order + subtensor identity
    #             + tuple(perm[p] for p in tcp_tuple[2:])
    #     )
    #
    # @classmethod
    # def _remap_tcp_inds(cls, p: 'TensorCoefficientPoly', mapping):
    #     return type(p)(
    #         {
    #             (cls._remap_tcp_tuple(k, mapping) for k in tensor_coeff_product): v
    #             for tensor_coeff_product, v in p.terms.items()
    #         },
    #         prefactor=p.scaling
    #     )


    # @classmethod
    # def _remap_tensor_coefficients(cls, coefficient_list, perm):
    #     # should be straight-forward if we can just permute the
    #     # indices of the appropriate tensor coefficient, but that's not assured?
    #     if isinstance(coefficient_list, np.ndarray) and coefficient_list.dtype != np.dtype(object):
    #         return coefficient_list
    #     else:
    #         return np.array(
    #             [
    #                 cls._remap_tcp_inds(p, perm) if isinstance(p, TensorCoefficientPoly) else p
    #                 for p in coefficient_list
    #             ],
    #             dtype=object
    #         )

    # def transpose(self, perm):
    #     """
    #     Reorders coeffs & attempts to reorder any tensor coefficients if necessary
    #
    #     :param perm:
    #     :return:
    #     """
    #
    #     remapping = np.argsort(perm)
    #
    #     return [
    #         self._remap_tensor_coefficients(self.coeffs[i], remapping)
    #         for i in perm
    #     ]

    # def __mul__(self, other):
    #     raise NotImplementedError(...)

    def __add__(self, other):
        if isinstance(other, (int, float, np.integer, np.floating)):
            if other == 0:
                return self
            else:
                return self + ProductPTPolynomial([[other]])
        elif isinstance(other, ProductPTPolynomial):
            if self.prefactor != 1 or other.prefactor != 1:
                raise NotImplementedError("need to include prefactors")
            return ProductPTPolynomialSum([self, other])
        elif isinstance(other, ProductPTPolynomialSum):
            return other + self
        else:
            raise NotImplementedError("not sure how to add {} and {}".format(self, other))
    def __radd__(self, other):
        return self + other

class ProductPTPolynomialSum:

    def __init__(self, polynomials, prefactor=1):
        self.polys = polynomials
        self.prefactor = prefactor
        self._order = None

    def __repr__(self):
        form_counts = {}
        for p in self.polys:
            key = "<{}>".format(",".join(str(c) for c in p.order))
            form_counts[key] = form_counts.get(key, 0) + 1
        return "{}({})".format(
            type(self).__name__,
            "+".join("{}{}".format(k,f) for k,f in form_counts.items())
        )

    @property
    def order(self):
        if self._order is None:
            suborders = [p.order for p in self.polys]
            max_len = max(len(o) for o in suborders)
            self._order = tuple(
                max(o[i] if i < len(o) else 0 for o in suborders)
                for i in range(max_len)
            )
        return self._order


    @classmethod
    def combine_polys(cls, poly_set, cache):
        new_polys = []
        eliminated_pos = set()
        for i,p1 in enumerate(poly_set):
            if i not in eliminated_pos:
                for j,p2 in enumerate(poly_set[i+1:]):
                    j = i+j+1
                    if j not in eliminated_pos and (p1, p2) not in cache:
                        # so we can avoid trying to recombine over passes
                        newp = p1.combine(p2)
                        if newp is not False: # combined
                            eliminated_pos.add(i)
                            eliminated_pos.add(j)
                            if newp is not True: #but didn't cancel
                                new_polys.append(newp)
                            break
                        else:
                            cache.add((p1, p2))
                else:
                    new_polys.append(p1) # made it through untouched

        return new_polys
    def combine(self):
        polys = self.polys
        cur_len = len(polys) + 1
        cache = set()
        while cur_len > len(polys): # hit a fixed point with this transformation
            cur_len = len(polys)
            polys = self.combine_polys(polys, cache)
        return type(self)(polys, prefactor=self.prefactor)

    def mul_simple(self, other:'ProductPTPolynomial'):
        if isinstance(other, ProductPTPolynomial):
            return type(self)(
                [
                    poly.mul_simple(other)
                    for poly in self.polys
                ],
                prefactor=self.prefactor
            )
        elif isinstance(other, ProductPTPolynomialSum):
            return type(self)(
                [
                    p1.mul_simple(p2)
                    for p1 in self.polys
                    for p2 in other.polys
                ],
                prefactor=self.prefactor*other.prefactor
            )
        else:
            raise NotImplementedError(type(other))

    def mul_along(self, other:'ProductPTPolynomial', inds, remainder=None):
        if isinstance(other, ProductPTPolynomial):
            # can just distribute except at some level I have to assume
            # that every stored poly has the same dimension...?
            return type(self)(
                [
                    poly.mul_along(other, inds, remainder=remainder) # how is this supposed to work?
                    for poly in self.polys # the fuck is inds supposed to mean here...????
                ],
                prefactor=self.prefactor
            )
        elif isinstance(other, ProductPTPolynomialSum):
            return type(self)(
                [
                    p1.mul_along(p2, inds, remainder=remainder)
                    for p1 in self.polys
                    for p2 in other.polys
                ],
                prefactor=self.prefactor*other.prefactor
            )
        else:
            raise NotImplementedError(type(other))


    # def __mul__(self, other):
    #     if isinstance(other, (ProductPTPolynomial, ProductPTPolynomialSum)):
    #         return self.mul_simple(other)
    #     elif isinstance(other, PTTensorCoeffProductSum):
    #         return other.mul_simple(self)
    #     else:
    #         return type(self)(self.polys, prefactor=self.prefactor*other)

    def __add__(self, other):
        if isinstance(other, (int, float, np.integer, np.floating)):
            if other == 0:
                return self
            else:
                return self + ProductPTPolynomial([[other]])
        else:
            if isinstance(other, ProductPTPolynomialSum):
                if other.prefactor != 1 or self.prefactor != 1:
                    print(other.prefactor, self.prefactor)
                    raise NotImplementedError("need to handle prefactors")
                return type(self)(self.polys + other.polys)
            elif isinstance(other, ProductPTPolynomial):
                if other.prefactor != 1 or self.prefactor != 1:
                    raise NotImplementedError("need to handle prefactors")
                return type(self)(self.polys + [other])
            else:
                raise TypeError("not sure what to do with {} and {}".format(type(self), type(other)))
    def __radd__(self, other):
        return self + other

class PTEnergyChangeProductSum(TensorCoefficientPoly):
    """
    A representation of a sum of 1/energy * poly sums
    which is here so we can transpose energy change indices intelligently
    """
    def __repr__(self):
        sums = []
        for k_prod,v in self.terms.items():
            ks = [
                "E-[{}]".format(k)
                for k in k_prod
            ]
            sums.append("{}{}".format(v, "".join(ks)))
        return "ESum({})".format("+".join(sums) if len(sums) < 3 else "\n      +".join(sums))

    @staticmethod
    def _permute_idx(idx, remapping):
        return idx[:2] + tuple(remapping[k] if k < len(remapping) else k for k in idx[2:])
    def mul_along(self, other:'ProductPTPolynomial', inds, remainder=None):
        """
        We multiply every subpoly along the given indices, transposing the appropriate tensor indices

        :param other:
        :param inds:
        :param remainder:
        :return:
        """
        if isinstance(other, PTEnergyChangeProductSum):
            return self.direct_product(other,
                                       key_func=lambda k1,k2:k1+tuple(self._permute_idx(idx, inds) for idx in k2),
                                       mul=lambda a,b:a.mul_along(b, inds, remainder=remainder)
                                       )
        else:
            raise NotImplementedError(type(other))

    def mul_simple(self, other:'ProductPTPolynomial'):
        """
        We multiply every subpoly along the given indices, transposing the appropriate tensor indices

        :param other:
        :param inds:
        :param remainder:
        :return:
        """
        if isinstance(other, PTTensorCoeffProductSum):
            return self.direct_product(other,
                                       # key_func=lambda k1,k2:k1+tuple(k2[i] for i in inds), #inds says how k2 permutes?
                                       mul=lambda a,b:a.mul_simple(b)
                                       )
        elif isinstance(other, (ProductPTPolynomial, ProductPTPolynomialSum)):
            return type(self)(
                {
                    k: v.mul_simple(other)
                    for k, v in self.terms.items()
                },
                prefactor=self.prefactor
            )

class PTTensorCoeffProductSum(TensorCoefficientPoly):
    """
    A representation for a sum of tensor coefficients * polynomial sums
    which is primarily here so we can transpose tensor coefficients intelligently
    """

    def __repr__(self):
        sums = []
        for k_prod,v in self.terms.items():
            ks = [
                "{}[{}]{}".format(
                    ["V", "G", "V'", "Z", "U"][k[1]], # we explicitly split Coriolis and Watson terms out
                    k[0],
                    k[2:]
                )
                for k in k_prod
            ]
            sums.append("{}{}".format(v, "".join(ks)))
        return "PTSum({})".format("+".join(sums) if len(sums) < 3 else "\n      +".join(sums))

    def combine(self):
        new_terms = {}
        for k,p in self.terms.items():
            if isinstance(p, ProductPTPolynomialSum):
                p = p.combine()
                if len(p.polys) > 0:
                    new_terms[k] = p
            else:
                new_terms[k] = p
        return type(self)(new_terms, prefactor=self.prefactor)

    @staticmethod
    def _potential_symmetrizer(idx):
        return tuple(sorted(idx))
    @staticmethod
    def _kinetic_symmetrizer(idx): # first and last and everything in between
        p_sorts = list(sorted([idx[0], idx[-1]]))
        rest_sorts = list(sorted(idx[1:-1]))
        return (p_sorts[0],) + tuple(rest_sorts) + (p_sorts[1],)
    @staticmethod
    def _coriolis_symmetrizer(idx): # second and last are momenta and everything else is x
        p_sorts = list(sorted([idx[1], idx[-1]]))
        rest_sorts = list(sorted((idx[0],) + idx[2:-1]))
        return (rest_sorts[0], p_sorts[0]) + tuple(rest_sorts[1:]) + (p_sorts[1],)

    @classmethod
    def symmetrizers(cls):
        return (
            cls._potential_symmetrizer,
            cls._kinetic_symmetrizer,
            cls._potential_symmetrizer,  # V'
            cls._coriolis_symmetrizer,
            cls._potential_symmetrizer  # "U
        )
    @classmethod
    def _symmetrize(cls, idx):
        return (idx[0], idx[1]) + cls.symmetrizers()[idx[1]](idx[2:])
    @classmethod
    def canonical_key(cls, monomial_tuple):
        return super().canonical_key(tuple(cls._symmetrize(t) for t in monomial_tuple))

    @staticmethod
    def _permute_idx(idx, remapping):
        return idx[:2] + tuple(remapping[k] if k < len(remapping) else k for k in idx[2:])
    def mul_along(self, other:'ProductPTPolynomial', inds, remainder=None):
        """
        We multiply every subpoly along the given indices, transposing the appropriate tensor indices

        :param other:
        :param inds:
        :param remainder:
        :return:
        """
        if isinstance(other, PTTensorCoeffProductSum):
            return self.direct_product(other,
                                       key_func=lambda k1,k2:k1+tuple(self._permute_idx(idx, inds) for idx in k2),
                                       mul=lambda a,b:a.mul_along(b, inds, remainder=remainder)
                                       )
        else:
            raise NotImplementedError(type(other))
    def mul_simple(self, other:'ProductPTPolynomial'):
        """
        We multiply every subpoly along the given indices, transposing the appropriate tensor indices

        :param other:
        :param inds:
        :param remainder:
        :return:
        """
        if isinstance(other, PTTensorCoeffProductSum):
            return self.direct_product(other,
                                       # key_func=lambda k1,k2:k1+tuple(k2[i] for i in inds), #inds says how k2 permutes?
                                       mul=lambda a,b:a.mul_simple(b)
                                       )
        elif isinstance(other, (ProductPTPolynomial, ProductPTPolynomialSum)):
            return type(self)(
                {
                    k: v.mul_simple(other)
                    for k, v in self.terms.items()
                },
                prefactor=self.prefactor
            )

# class InverseEnergyFactor:
#     def __init__(self, changes, prefactor=1):
#         self.changes = changes
#         self.prefactor = prefactor
#
#     def __repr__(self):
#         return "DE[{}]".format(tuple(self.changes))
#
#     def __mul__(self, other):
#         if isinstance(other, (int, float, np.integer, np.floating)):
#             if other == 0:
#                 return 0
#             elif other == 1:
#                 return self
#         if isinstance(other, (ProductPTPolynomial, ProductPTPolynomialSum, PTTensorCoeffProductSum)):
#             return other.__mul__(self)
#         else:
#             return type(self)(self.changes, other*self.prefactor)
#     def __rmul__(self, other):
#         return self * other

class FullChangePathPolyTerm:
    """
    Built by taking products of wave function corrections with the Hamiltonian
    expansion elements
    """
    def __init__(self, terms:"dict[tuple[int],ProductPTPolynomialSum]", term_sizes=None):
        self.terms = terms
        if term_sizes is None:
            term_sizes = {sum(abs(x) for x in k) for k in terms.keys()}
        self.term_sizes = term_sizes

    # def _build_poly(self, poly_1, poly_2, poly_1_inds, rem_inds, transpose):
    #     """
    #     We multiply poly_1 and poly_2 along the specified `poly_1_inds`, just convolving
    #     the polynomial coefficients for each axis and copying in the remaining inds of poly_2
    #
    #     :param poly_1:
    #     :param poly_2:
    #     :param poly_1_inds:
    #     :param rem_inds:
    #     :param transpose:
    #     :return:
    #     """

    def op_mult(self, op):

        new_terms = {}
        for term_1, poly_1 in self.terms.items():
            for term_2, poly_2 in op.terms.items():
                l1 = len(term_1)
                l2 = len(term_2)
                # min_poly_len = max(len(term_1), len(term_2)) # everything lines up
                # max_poly_len = len(term_1) + len(term_2) # nothing lines up

                # we'll take partitions of term_2 for every possible overlap size
                # then permute the overlapping part and add that
                # this _should_ give the minimal number of additions
                # the partitions themselves are maybe allowed to be unique
                # although in the cases of reducing over unique partitions I guess
                # it's possible that we need to worry about the different permutations
                # of tensor coefficient terms???

                #TODO: handle overlap == 0 term

                perm_inds = np.array(list(itertools.permutations(range(l1))))
                inv_perm = np.argsort(perm_inds)
                for overlaps in range(1, min(l1, l2)+1):
                    #TODO: handle overlap == len(term_2)

                    (ov_parts, rem_part), (ov_inds, rem_inds) = UniquePartitions(term_2).partitions(
                        [overlaps, l2 - overlaps],
                        return_indices=True
                    ) #!need to track indices
                    if l1 == overlaps:
                        pad_ov = ov_parts
                    else:
                        pad_ov = np.pad(ov_parts, [[0, 0], [0, l1-overlaps]])
                    for p, r, io, ir in zip(pad_ov, rem_part, ov_inds, rem_inds):
                        perm_ov = p[perm_inds] # I think I can actually take the unique perms...?
                        perm_ov, uperm_inds = np.unique(perm_ov, return_index=True, axis=0)

                        # it might be necessary to apply the transpose here

                        newts = np.concatenate([
                            np.array(term_1)[np.newaxis] + perm_ov,
                            r[np.newaxis]
                        ],
                            axis=1
                        )

                        newt_transposes = np.argsort(newts, axis=1)

                        # build new polys
                        target_axes = inv_perm[uperm_inds][:, :overlaps] # axes in term_1 that the selected parts of
                                                                         # term_2 are multiplied by
                        for full_changes, transp, targ_inds in zip(newts[newt_transposes], newt_transposes, target_axes):
                            self._build_new_poly(poly_1, poly_2, targ_inds)



                # for plen in range(min_poly_len, max_poly_len+1):
                #     # we want to take all permutations of the padded version of term_2
                #     # that will create distinct polynomials...
                #     ...



        # take direct sum of terms, accounting for changes
        # and permute so that they remain ordered --> permuting tensor coefficients too
        raise NotImplementedError(...)

class PerturbationTheoryTerm(metaclass=abc.ABCMeta):
    """
    A generic version of one of the three terms in
    PT that will generate a correction polynomial
    """
    def __init__(self):
        self._exprs = None
        self._changes = None

    def get_subexpresions(self) -> 'Iterable[PerturbationTheoryTerm]':
        raise NotImplementedError("just here to be overloaded")
    def __mul__(self, other):
        if isinstance(other, PerturbationTheoryTerm):
            return PerturbationTheoryTermProduct(self, other)
        else:
            return ScaledPerturbationTheoryTerm(self, other)
    def __rmul__(self, other):
        # if isinstance(other, PerturbationTheoryTerm): # can't happen
        #     return ...
        return ScaledPerturbationTheoryTerm(self, other)
    def __neg__(self):
        return ScaledPerturbationTheoryTerm(self, -1)

    @property
    def expressions(self):
        if self._exprs is None:
            self._exprs = self.get_subexpresions()
        return self._exprs

    @abc.abstractmethod
    def get_changes(self) -> 'dict[tuple[int], Any]':
        ...
    @property
    def changes(self):
        if self._changes is None:
            self._changes = self.get_changes()
        return self._changes


    def get_poly_terms(self, changes, shift=None) -> 'PTTensorCoeffProductSum':
        return sum(
            term.get_poly_terms(changes, shift=shift)
            for term in self.expressions
        )

    def __call__(self, changes, shift=None, coeffs=None, freqs=None, check_sorting=True):
        if coeffs is not None:
            raise NotImplementedError(...)
        if freqs is not None:
            raise NotImplementedError(...)
        if check_sorting:
            t = np.asanyarray([c for c in changes if c != 0])
            sort_key = -(np.abs(t) * 10 - (t < 0))
            changes = tuple(t[np.argsort(sort_key)])
        terms = self.changes.get(changes, None)
        if not isinstance(terms, (PTTensorCoeffProductSum, int, np.integer)):
            self.changes[changes] = self.get_poly_terms(changes, shift=None).combine()
        terms = self.changes[changes]
        if shift is not None:
            raise NotImplementedError("reshifting not supported yet...")
        return PerturbationTheoryEvaluator(terms, len(changes))

class HamiltonianExpansionTerm(PerturbationTheoryTerm):
    def __init__(self, terms, order=None, identities=None, symmetrizers=None):
        super().__init__()

        self.terms = terms
        self.order = order
        self.identities = np.arange(len(terms)) if identities is None else identities
        if symmetrizers is None:
            symmetrizers = PTTensorCoeffProductSum.symmetrizers()
        self.symmetrizers = symmetrizers

        self._grouped_terms = {}
        for t in self.terms:
            self._grouped_terms[len(t)] = self._grouped_terms.get(len(t), [])
            self._grouped_terms[len(t)].append(t)

        self.term_sizes = set(self._grouped_terms.keys())
        #TODO: need to include parity effects for all integer partitions"
        # self.changes =

        self._change_poly_cache = {}

    def __repr__(self):
        return "H[{}]".format(self.order)

    _generator_cache = {}
    @classmethod
    def _get_generator(cls, term_list, inds):
        ts = tuple(term_list[i] for i in inds)
        if ts not in cls._generator_cache:
            cls._generator_cache[ts] = HarmonicOscillatorMatrixGenerator(ts) # no reason not to cache these...
        return cls._generator_cache[ts]


    @classmethod
    def _evaluate_poly_coeffs(cls, terms, inds, delta, shift):
        gen = cls._get_generator(terms, inds)
        poly_contrib = gen.poly_coeffs(delta, shift=shift)
        if isinstance(poly_contrib, (int, float, np.integer, np.floating)) and poly_contrib == 0:
            return 0
        # now we get any remainder we need to
        sqrt_contrib = HarmonicOscillatorRaisingLoweringPolyTerms.get_direction_change_poly(delta, shift)
        if sqrt_contrib is not None:
            poly_contrib = scipy.signal.convolve(poly_contrib, sqrt_contrib)
        return poly_contrib

    def get_changes(self) -> 'dict[tuple[int], Any]':
        changes = {}
        l_check = set()
        for t in self.terms:
            start = len(t) % 2
            for tl in range(start, len(t)+1, 2): # I could really just get the max_len but w/e
                if tl == 0:
                    changes[()] = None
                    continue
                if tl not in l_check: # just to not redo it...
                    l_check.add(tl)
                    for p in IntegerPartitioner.partitions(tl, pad=False):
                        l = p.shape[1]
                        for n in range(0, l+1):
                            base_perm = [1] * (l-n) + [-1] * n
                            signs = np.array(list(itertools.permutations(base_perm)))
                            resigned_perms = (signs[np.newaxis] * p[:, np.newaxis]).reshape(-1, l)
                            sort_key = -(10*p[:, np.newaxis] - (signs < 0)[np.newaxis]).reshape(-1, l)
                            _, uinds = np.unique(np.sort(sort_key), axis=0, return_index=True)
                            resigned_perms = resigned_perms[uinds]
                            for subp in resigned_perms:
                                changes[tuple(subp)] = None
        return changes

    def get_subexpresions(self) -> 'Iterable[PerturbationTheoryTerm]':
        raise NotImplementedError("shouldn't need this here...")

    def get_poly_terms(self, changes, shift=None) -> 'PTTensorCoeffProductSum': #TODO: CACHE THIS SHIT
        #NOTE: this explicitly excludes the sqrt contribution to the poly

        poly_contribs = {}

        abs_change = np.array([abs(c) for c in changes])
        total_change = sum(abs_change)
        for group_size, terms in self._grouped_terms.items():
            remainder = group_size - total_change
            if remainder < 0 or remainder % 2 != 0:
                continue # literally cannot execute this change

            # not sure what to do if group_size == 0...?
            if group_size == 0: # constant contrib
                for term_index, term_list in enumerate(terms):
                    subpolys = [
                        ProductPTPolynomial([np.array([1])])
                    ]
                    prefactor = (self.order, self.identities[term_index])  # type: tuple[int]
                    poly_contribs[(prefactor,)] = poly_contribs.get((prefactor,), 0) + ProductPTPolynomialSum(subpolys)
                continue

            total_dim = np.count_nonzero(abs_change) + remainder // 2

            # now we enumerate all of the ways to partition `group_size` elements across `total_dim` dimensions
            # such that each `partition - changes` is even in each dimension
            abs_change = np.pad(abs_change, [0, total_dim - len(changes)])
            changes = np.pad(changes, [0, total_dim - len(changes)]).astype(int)
            if shift is None:
                shift = np.zeros(len(changes), dtype=int)
            else:
                shift = np.pad(shift, [0, total_dim - len(shift)]).astype(int)

            targets = SymmetricGroupGenerator(total_dim).get_terms(group_size)#.partitions(group_size, max_len=total_dim, pad=True)
            check_changes = targets - abs_change[np.newaxis, :]
            good_targets = np.all(np.logical_and(check_changes >= 0, check_changes % 2 == 0), axis=1)
            # bad_targets = np.full(len(targets), True, dtype=bool)
            # bad_pos = np.arange(len(targets))
            # for cc in itertools.permutations(abs_change):
            #     check_changes = targets[bad_targets] - np.array(cc)[np.newaxis, :]
            #     good_targets = np.all(np.logical_and(check_changes >= 0, check_changes % 2 == 0), axis=1)
            #     bad_targets[bad_pos[bad_targets][good_targets]] = False
            # good_targets = np.logical_not(bad_targets)
            targets = targets[good_targets]

            # then for each of these partition sizes, we enumerate the actual partitions of the terms
            # and then take the corresponding terms to generate partitions
            # (and multiply by the appropriate tensor coefficients)

            partitioner = UniquePartitions(np.arange(group_size))
            for partition_sizes in targets:
                parts, inv_splits = partitioner.partitions(partition_sizes,
                                                          return_partitions=False, take_unique=False,
                                                          return_indices=True, return_inverse=True
                                                          )

                inv_splits = np.concatenate(inv_splits, axis=1)
                base_tup = np.array(sum(([i]*s for i,s in enumerate(partition_sizes)), []))
                for j,p_vec in enumerate(zip(*parts)):
                    tensor_idx = tuple(base_tup[inv_splits[j]])
                    for term_index,term_list in enumerate(terms):
                        poly_coeffs = [
                            self._evaluate_poly_coeffs(term_list, inds, delta, s)
                            for inds, delta, s in zip(p_vec, changes, shift)
                        ]
                        if any(isinstance(c, (int, float, np.integer, np.floating)) and c == 0 for c in poly_coeffs):
                            # total thing is zero
                            continue
                        if self.symmetrizers[term_index] is not None:
                            symm_idx = self.symmetrizers[term_index](tensor_idx)
                        else:
                            symm_idx = tensor_idx
                        prefactor = (self.order, self.identities[term_index]) + symm_idx #type: tuple[int]
                        poly_contribs[(prefactor,)] = poly_contribs.get((prefactor,), 0) + ProductPTPolynomial(poly_coeffs)

        return PTTensorCoeffProductSum(poly_contribs)

    # def get_correction_poly(self, changes, shifts=None):
    #     #TODO: check to make sure changes is sorted??? Or not???
    #     t = tuple(changes)
    #     if t not in self._change_poly_cache:
    #         self._change_poly_cache[t] = FullChangePathPolyTerm(
    #             self.get_poly_terms(changes, shifts)
    #         )
    #     return self._change_poly_cache[t]

    # def __call__(self, changes, shifts=None):
    #     return self.get_correction_poly(changes, shifts)

class PerturbationOperator(PerturbationTheoryTerm):
    def __init__(self, subterm):
        super().__init__()

        self.subterm = subterm

    def __repr__(self):
        return "P{}".format(self.subterm)

    def get_changes(self) -> 'dict[tuple[int], Any]':
        return {k:None for k in self.subterm.changes if k != ()}

    def get_poly_terms(self, changes, shift=None) -> 'PTTensorCoeffProductSum':
        if len(changes) == 0:
            return 0

        base_term = self.subterm.get_poly_terms(changes, shift=shift)
        if isinstance(base_term, PTTensorCoeffProductSum):
            if shift is not None:
                changes = tuple(c + s for c,s in zip(changes, shift))
            prefactor = InverseEnergyFactor(changes)
            base_term = base_term * prefactor

        return base_term

class WavefunctionCorrection(PerturbationTheoryTerm):
    def __init__(self, parent, order):
        super().__init__()

        self.parent = parent
        self.order = order

    def __repr__(self):
        return "Y[{}]".format(self.order)

    # def get_poly_terms(self, changes, shift=None) -> 'PTTensorCoeffProductSum':
    #     if self.order == 0:
    #         if len(changes) == 0:
    #             return PTTensorCoeffProductSum({():ProductPTPolynomial([[1]])}) # constant contrib?
    #         else:
    #             return 0
    #     else:
    #         return super().get_poly_terms(changes, shift=shift)

    def get_changes(self):
        base_changes = {}
        for expr in self.expressions:
            for subchange in expr.changes:
                base_changes[subchange] = base_changes.get(subchange, [])
                base_changes[subchange].append(expr) # just track which exprs generate the change
        return base_changes

    def get_subexpresions(self):
        H = self.parent.hamiltonian_expansion
        E = self.parent.energy_correction
        W = self.parent.wavefunction_correction
        k = self.order

        base_terms = [
            E(k-i) * W(i)
            for i in range(0, k)
            if i > 0
        ] + [
            -PerturbationOperator(H[k-i] * W(i))
                if i > 0 else
            -PerturbationOperator(H[k])
            for i in range(0, k)
        ]

        #TODO: include energy weighting

        return base_terms

class WavefunctionOverlapCorrection(PerturbationTheoryTerm):
    """
    Provides a slight optimization on the base `WavefunctionCorrection`
    """

    def __init__(self, parent, order):
        super().__init__()

        self.parent = parent
        self.order = order

    def __repr__(self):
        return "O[{}]".format(self.order)

    def get_changes(self) -> 'dict[tuple[int], Any]':
        return {():None}

    def get_subexpresions(self) -> 'Iterable[PerturbationTheoryTerm]':
        W = self.parent.wavefunction_correction
        k = self.order

        if k == 1:
            return []

        return [
            -1/2 * (W[k - i] * W(i)) for i in range(1, k)
        ]

class EnergyCorrection(PerturbationTheoryTerm):
    def __init__(self, parent, order):
        super().__init__()

        self.parent = parent
        self.order = order

    def __repr__(self):
        return "E[{}]".format(self.order)

    def get_changes(self) -> 'dict[tuple[int], Any]':
        return {():None}

    def get_subexpresions(self) -> 'Iterable[PerturbationTheoryTerm]':
        H = self.parent.hamiltonian_expansion
        E = self.parent.energy_correction
        W = self.parent.wavefunction_correction
        O = self.parent.overlap_correction
        k = self.order

        if k % 2 == 1:
            return []

        return [H[k]] + [
            (H[k - i] * W(i)) for i in range(1, k)
        ] + [
            -E(k-i) * O(i) for i in range(1, k)
        ]

    def get_poly_terms(self, changes, shift=None):
        if len(changes) != 0:
            return 0
        else:
            return super().get_poly_terms(changes, shift=shift)

    # def get_poly_terms(self):
    #     H = self.parent.hamiltonian_expansion
    #     E = self.parent.energy_correction
    #     W = self.parent.wavefunction_correction
    #     k = self.order
    #
    #     # Recursive PT equations
    #     return H[k]([]) + sum(
    #         (H[k-i]*W(i))([]) - E(k-i)([])*W(i)([]) #TODO: implement caching for this kind of term...
    #         for i in range(1, k)
    #     ) # 1 to k-1

    def __call__(self, changes, shift=None, coeffs=None, freqs=None, check_sorting=None):
        if len(changes) > 0:
            raise ValueError("no")
        return super().__call__((), shift=shift, coeffs=coeffs, freqs=freqs, check_sorting=False)

class ScaledPerturbationTheoryTerm(PerturbationTheoryTerm):
    #TODO: refactor since inheritance isn't really the right paradigm here
    def __init__(self, base_term:'PerturbationTheoryTerm', scaling):
        super().__init__()
        self.prefactor = scaling
        self.base = base_term


    def __repr__(self):
        if isinstance(self.prefactor, (int, np.integer)) and self.prefactor == -1:
            return "-{}".format(self.base)
        else:
            return "{}*{}".format(self.prefactor, self.base)

    def get_changes(self) -> 'dict[tuple[int], Any]':
        return {k:None for k in self.base.changes}
    def get_poly_terms(self, changes, shift=None) -> 'PTTensorCoeffProductSum':
        base = self.base.get_poly_terms(changes, shift=shift)
        return self.prefactor * base

class PerturbationTheoryTermProduct(PerturbationTheoryTerm):
    def __init__(self, gen1, gen2):
        super().__init__()

        self.gen1 = gen1
        self.gen2 = gen2


    def __repr__(self):
        return "{}*{}".format(self.gen1, self.gen2)

    @classmethod
    def change_direct_sum_generator(self, change_1, change_2):
        """
        Enumerates the direct sum + indices for terms taken from `change_1` and `change_2`

        :param change_1:
        :param change_2:
        :return:
        """
        l1 = len(change_1)
        l2 = len(change_2)
        # min_poly_len = max(len(term_1), len(term_2)) # everything lines up
        # max_poly_len = len(term_1) + len(term_2) # nothing lines up

        # we'll take partitions of term_2 for every possible overlap size
        # then permute the overlapping part and add that
        # this _should_ give the minimal number of additions
        # the partitions themselves are maybe allowed to be unique
        # although in the cases of reducing over unique partitions I guess
        # it's possible that we need to worry about the different permutations
        # of tensor coefficient terms???

        # TODO: handle overlap == 0 term

        perm_inds = np.array(list(itertools.permutations(range(l1))))
        inv_perm = np.argsort(perm_inds)
        for overlaps in range(1, min(l1, l2) + 1):
            # TODO: handle overlap == len(term_2)

            (ov_parts, rem_part), (ov_inds, rem_inds) = UniquePartitions(change_2).partitions(
                [overlaps, l2 - overlaps],
                return_indices=True
            )  # !need to track indices
            if l1 == overlaps:
                pad_ov = ov_parts
            else:
                pad_ov = np.pad(ov_parts, [[0, 0], [0, l1 - overlaps]])
            for p, r, io, ir in zip(pad_ov, rem_part, ov_inds, rem_inds):
                perm_ov = p[perm_inds]  # I think I can actually take the unique perms...?
                perm_ov, uperm_inds = np.unique(perm_ov, return_index=True, axis=0)

                # it might be necessary to apply the transpose here

                if len(r) > 0:
                    newts = np.concatenate([
                        np.array(change_1)[np.newaxis] + perm_ov,
                        np.broadcast_to(r[np.newaxis], (len(perm_ov), len(r)))
                    ],
                        axis=1
                    )
                else:
                    newts = np.array(change_1)[np.newaxis] + perm_ov

                newt_transposes = np.argsort(newts, axis=1)

                # build new polys
                target_axes = inv_perm[uperm_inds][:, :overlaps]  # axes in term_1 that the selected parts of
                from_inds = perm_inds[uperm_inds][:, :overlaps]  # axes in term_1 that the selected parts of
                # term_2 are multiplied by
                for base_change, transp, targ_inds, src_inds in zip(newts, newt_transposes, target_axes, from_inds):
                    yield base_change[transp], transp, targ_inds, src_inds

    def get_changes(self):
        changes = {}
        for change_1 in self.gen1.changes:
            for change_2 in self.gen2.changes:
                for change_sum, transpose, target_indices, src_inds in self.change_direct_sum_generator(change_1, change_2):
                    sorting = np.argsort(-(np.abs(change_sum) * 10 - (change_sum < 0)))
                    k = tuple(c for c in change_sum[sorting] if c != 0)
                    clist = changes.get(k, [])
                    changes[k] = clist
                    clist.append(
                        [change_1, change_2, transpose, target_indices, src_inds] # stuff to build the term later if needed
                    )
        return changes
    # @property
    # def changes(self):
    #     if self._changes is None:
    #         self._changes = self.build_changes()
    #     return self._changes

    def get_expressions(self):
        raise NotImplementedError("shouldn't need this here...")

    def get_poly_terms(self, changes, shift=None): # need to align this with
        t = tuple(changes)
        # check sorting?
        base_changes = self.changes.get(t, 0)
        if isinstance(base_changes, list):
            poly_changes = 0
            for change_1, change_2, reorg, target_inds, src_inds in base_changes:
                polys_1 = self.gen1.get_poly_terms(change_1, shift=shift)
                # we note that src_inds defines how we permute the change_2 inds so...I guess target_inds
                # tells us how we'd permute the inds of change_1?
                if shift is not None:
                    change_shift = change_1 + shift
                else:
                    change_shift = change_1
                change_shift = [change_shift[i] for i in target_inds]
                polys_2 = self.gen2.get_poly_terms(change_2, shift=change_shift)
                base_polys = polys_1.mul_along(polys_2, target_inds, remainder=None)
                # we need to also account for both possible direction changes...
                if shift is not None:
                    sqrt_contrib_1 = HarmonicOscillatorRaisingLoweringPolyTerms.get_direction_change_poly(change_1, shift)
                else:
                    sqrt_contrib_1 = None
                if sqrt_contrib_1 is not None:
                    # for every term in base_polys we need to convolve with this change...
                    base_polys = base_polys.mul_simple(ProductPTPolynomial(sqrt_contrib_1))

                sqrt_contrib_2 = HarmonicOscillatorRaisingLoweringPolyTerms.get_direction_change_poly(change_2, change_shift)
                if sqrt_contrib_2 is not None:
                    # need to permute
                    sqrt_contrib_2 = [sqrt_contrib_2[i] for i in target_inds]
                    base_polys = base_polys.mul_simple(ProductPTPolynomial(sqrt_contrib_2))
                poly_changes += base_polys

            self.changes[t] = poly_changes
            base_changes = poly_changes

        return base_changes

class PerturbationTheoryEvaluator:
    def __init__(self, expr:'PTTensorCoeffProductSum', num_fixed):
        self.expr = expr
        self.num_fixed = num_fixed
    def __repr__(self):
        return "{}(<{}>, {})".format(type(self).__name__, self.num_fixed, self.expr)

    @staticmethod
    def _extract_coeffs(coeffs, coeff_indices, fixed, free):
        coeff_tensor = coeffs[coeff_indices[0]][coeff_indices[1]] # order then identity
        if isinstance(coeff_tensor, (int, float, np.integer, np.floating)) and coeff_tensor == 0:
            return 0

        idx = tuple(free[j-fixed] if j > fixed else j for j in coeff_indices[2:])
        if len(idx) > 0:
            return coeff_tensor[idx]
        else:
            return coeff_tensor # just a number
    def evaluate(self, state, coeffs, freqs, zero_cutoff=None):
        # we do state-by-state evaluation for now although we will at some
        # point need to do this in batches
        ndim = len(state)
        # max_order = None
        free_ind_groups = {}
        for coeff_indices,poly_terms in self.expr.terms.items():
            num_inds = np.unique(np.concatenate(coeff_indices))
            free_inds = len(num_inds) - self.num_fixed
            if free_inds not in free_ind_groups:
                free_ind_groups[free_inds] = []
            free_ind_groups[free_inds].append(coeff_indices)
        #     max_order = max(max_order, poly_terms.order) if max_order is not None else poly_terms.order
        #
        # state_polys = np.power(state, np.arange(max_order+1))

        if zero_cutoff is None:
            zero_cutoff = 0

        # TODO: numba-ify this part
        for free_inds,cind_sets in free_ind_groups.items():
            # now we iterate over every subset of inds (outside of the fixed ones)
            # excluding replacement (since that corresponds to a different coeff/poly)
            for subset in itertools.combinations(range(self.num_fixed, ndim), r=free_inds):
                prefactors = np.array([
                    np.prod([self._extract_coeffs(coeffs, ci, self.num_fixed, subset) for ci in cinds])
                    for cinds in cind_sets
                ])
                good_pref = np.where(np.abs(prefactors) > zero_cutoff) # don't bother to include useless terms
                if len(good_pref) == 0:
                    continue
                good_pref = good_pref[0]

                for g in good_pref:
                    subexpr = self.expr.terms[cind_sets[g]]
                    base_prefactor = subexpr.prefactor
                    if isinstance(subexpr, ProductPTPolynomialSum):
                        subprefactors = [p.prefactor for p in subexpr.polys]
                        poly_coeffs = [p.coeffs for p in subexpr.polys]
                    else:
                        subprefactors = []
                        poly_coeffs = subexpr.coeffs
                    print(poly_coeffs)

                # raise Exception(cind_sets)

    # def subtitute_coeffs(self, coeffs):
    #     """
    #     :param coeffs: V expansion, G expansion, U expansion, and then Coriolis expansions
    #     :return:
    #     """
    #     if not isinstance(self.expr, PTTensorCoeffProductSum):
    #         raise ValueError("coeffs have already been subbed")
    #     new_terms = []
    #     for k,p in self.expr.terms.items(): # this won't fully work since I need to map over possible pairs...
    #         c = coeffs[k[1]][k[0]]
    #         if len(k) > 2:
    #             c = c[k[2:]]
    #         if c != 0:
    #             new_terms.append(p*c)
    #     return type(self)(new_terms)
    #
    # def substitute_freqs(self, freqs):




#region Older Attempt (More Effficient?)

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
        """
        Holder for correction polys, indexed by dimension & path with paths given by `shifts_1D`

        :param polys_1D:
        :param shifts_1D:
        :param changes:
        """
        if (
                isinstance(polys_1D, (int, float, np.integer, np.floating, TensorCoefficientPoly))
                or isinstance(polys_1D[0], (int, float, np.integer, np.floating, TensorCoefficientPoly))
                or isinstance(polys_1D[0][0], (int, float, np.integer, np.floating, TensorCoefficientPoly))
        ): # means we fucked up
            raise TypeError("needed a list of lists of arrays of 1D poly coeffs")
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


    @staticmethod
    def _extract_tensor_coeff(coeff_tensors, k):
        base = coeff_tensors[k[0]][k[1]]
        return base if not isinstance(base, np.ndarray) else base[k[2:]]
    def _substitute_coeffs(self, tensor_poly, coeff_tensors):

        contribs = [
            scaling * np.prod([self._extract_tensor_coeff(coeff_tensors, k) for k in index_tuple_tuple])
            for index_tuple_tuple, scaling in tensor_poly.terms.items()
        ]
        return tensor_poly.scaling*np.sum(contribs)

    def _aggregate_subbed_coeff_lists(self, coefficient_tensors, poly_set):
        if poly_set is None:
            return np.array([1])
        else:
            sublists =  [
                    self._substitute_coeffs(tcp, coefficient_tensors)
                        if isinstance(tcp, TensorCoefficientPoly) else
                    tcp
                    for tcp in poly_set
                ]
            return sublists #np.sum(sublists, axis=0)

    def substitute_coefficients(self, coefficient_tensors):
        return type(self)(
            [
                [
                    self._aggregate_subbed_coeff_lists(coefficient_tensors, path_poly_set)
                    for path_poly_set in polys_1D
                ]
                for polys_1D in self.polys
            ],
            self.shifts,
            self.changes
        )

    def _eval_shift(self, path, change, crossover):
        #TODO: what do I do with terms that have a (Pi_n)(Pi_m) in them???

        path = np.concatenate([[0], path]) # eval strategy drops this...
        if crossover == -1 or crossover == len(path):
            return path
        elif crossover == 0:
            return np.concatenate([[0], path]) - path
        else:
            # we just go up normally from the left (starting at 0)
            # but from the right we need to effectively start from the final state
            num_right = len(path) - crossover
            num_left = crossover
            return path - np.pad(np.full(num_right, change), [num_left, 0])

    def get_quanta_shifts(self, crossover):
        return [
            np.array([
                self._eval_shift(path, delta, crossover)
                for path in np.array(paths_1D)
            ])
            for delta,paths_1D in zip(self.changes, self.shifts)
        ]
    def substitute_paths(self, freqs, exponents, crossover):
        shifts = self.get_quanta_shifts(crossover)
        return [
            [
                np.array([
                    np.power(freq * path, exponents)
                    # self._eval_energy(freq, path)
                    for path in path_set
                ])
                for path_set in paths_1D
            ]
            for freq,paths_1D in zip(freqs, shifts)
        ]

    # class direct_product_evaluator:
    #     def __init__(self, polynomial_coeffs, paths, freqs):
    #         self.pcs = polynomial_coeffs
    #         self.paths = paths
    #         self.freqs = freqs
    #     def evaluate(self, initial_states, path_exponents, crossover):
    #         # take direct product of numerators
    #         # and divide by direct sum product of denominators
    #
    #         raise NotImplementedError("need to eval direct product")

    # def get_poly_evaluator(self, coefficient_tensors):
    #     numerators = self.get_numerators(coefficient_tensors)
    #
    #     return self.direct_product_evaluator(numerators, self.shifts, freqs)

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

    # class direct_product_evaluator_sum:
    #     def __init__(self, evaluators):
    #         self.evaluators = evaluators
    #
    #
    # def get_poly_evaluator(self, coefficient_tensors, freqs):
    #     return self.direct_product_evaluator_sum([
    #         contrib.get_poly_evaluator(coefficient_tensors, freqs)
    #         for contrib in self.contribs
    #     ])
    #
    # def eval(self, quanta, coefficient_tensors, freqs, crossover):
    #     ...

    def substitute_coefficients(self, coefficient_tensors):
        return type(self)(tuple(
            contrib.substitute_coefficients(coefficient_tensors)
            for contrib in self.contribs
        ))

    def substitute_paths(self, freqs, exponents, crossover):
        return type(self)(tuple(
            contrib.substitute_paths(freqs, exponents, crossover)
            for contrib in self.contribs
        ))

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

    # def evaluate(self, quanta, coefficient_tensors, freqs, crossover):
    #     for path, contrib_gen in self.terms.items():
    #         contrib = contrib_gen.eval(quanta, coefficient_tensors, freqs, crossover)
    #     ...

    def substitute_coefficients(self, coeffs):
        return type(self)({
            path:sub.substitute_coefficients(coeffs)
            for path,sub in self.terms.items()
        })

    def substitute_paths(self, freqs, exponents, crossover):
        return type(self)({
            path:sub.substitute_paths(freqs, exponents, crossover)
            for path,sub in self.terms.items()
        })

# class PTTensorCoeffProductSum(TensorCoefficientPoly):
#     """
#     Only here for the better formatting
#     """
#     def __repr__(self):
#         sums = []
#         for k_prod,v in self.terms.items():
#             ks = [
#                 "{}[{}]{}".format(
#                     ["V", "G", "V'"][k[1]],
#                     k[0],
#                     k[2:]
#                 )
#                 for k in k_prod
#             ]
#             # v = round(v, 5)
#             sums.append("{}{}".format(v, "".join(ks)))
#         return "+".join(sums)

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
                        np.array([TensorCoefficientPoly.monomial(self.operator_indices[i][j])], dtype=object)
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

        return AnalyticPTCorrectionSum({ # Not even sure this makes sense to do
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

    @classmethod
    def get_integer_partition_permutation_tuples(cls, order, unique=True, max_part=None):
        base_parts = IntegerPartitioner.partitions(order, pad=False)
        if max_part is not None:
            base_parts = [b[np.max(b) < max_part] for b in base_parts]
        perm_iterator = itertools.chain(*(  # few enough I can be a bit inefficient
            list(itertools.permutations(part))
            for order_parts in base_parts
            for part in order_parts
        ))
        if unique:
            return set(perm_iterator)
        else:
            return list(perm_iterator)
    def energy_terms(self, order, return_indices=False):
        term_indices = self.get_integer_partition_permutation_tuples(order)

        terms = [
            self._take_term_direct_prod([self.terms[i] for i in term_index_list])
            for term_index_list in term_indices
        ]

        if return_indices:
            return terms, term_indices
        else:
            return terms

    def energy_correction_driver(self, order):
        return AnalyticPTEnergyDriver(
            [
                self.energy_terms(o, return_indices=True)
                for o in range(1, order+1)
            ],
            self
        )

    def overlap_terms(self, order, return_indices=False):
        term_indices = self.get_integer_partition_permutation_tuples(order, max_part=order-1)

        terms = [
            self._take_term_direct_prod([self.terms[i] for i in term_index_list])
            for term_index_list in term_indices
        ]

        if return_indices:
            return terms, term_indices
        else:
            return terms

    def overlap_correction_driver(self, order):
        return AnalyticPTOverlapDriver(
            [
                self.overlap_terms(o, return_indices=True)
                for o in range(1, order+1)
            ],
            self
        )

class AnalyticPTOperatorElementDriver(metaclass=abc.ABCMeta):
    def __init__(self, terms_inds, changes, parent:AnalyticPerturbationTheoryDriver):
        self.parent = parent
        self.changes = changes
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
                term_map[inds] = sum(AnalyticPTCorrectionGenerator(t).get_correction(self.changes) for t in terms)
        return term_map

    @abc.abstractmethod
    def enumerate_overlap_energy_exponent_tuples(self, index_tuple):
        #TODO: figure out parities

        remainder = self.order - sum(index_tuple)
        nterms = len(index_tuple)
        if remainder == 0: # just put a Pi between every term
            return [0], [[]], [[0] + [1]*(nterms - 1) + [0]]
        elif len(index_tuple) == 1: # only one term, so can only multiply by the overlap
            return [remainder], [[]], [[0, 0]]
        else: # need to figure out how to split up terms across the energies and overlaps
            # somehow we only ever have one overlap term so we'll just loop through those
            for ov in range(remainder+1):
                e_terms = []
                # enumerate all possible
                e_parts = IntegerPartitioner.partitions(remainder - ov)
                for part in e_parts:
                    engs, counts = np.unique(part, return_counts=True)
            raise NotImplementedError("still working out how to split integer partition permutations over E and Ov")


    def get_poly_evaluator(self, coefficient_tensors, freqs):
        poly_terms = self.generic_corrections
        polys = []
        for index_tuple, poly_sum in poly_terms.items():
            # evaluate sets of overlap * energy_correction * Pi's to weight polynomials
            if isinstance(poly_sum, (int, float, np.integer, np.floating)) and poly_sum == 0:
                continue

            base_poly = poly_sum.substitute_coefficients(coefficient_tensors)
            overlap, ecorrs, exponents = self.enumerate_overlap_energy_exponent_tuples(index_tuple)
            # poly_scaling = 0
            for ov, ec, exp in zip(overlap, ecorrs, exponents): # gotta dd up all these products...
                if ov == 1 or any(e%2 == 1 for e in ec): # ov[1] always zero as well as odd-order energy corrs
                    continue

                subpoly = base_poly.substitute_paths(freqs, exponents, -1)

                sub_scaling = 1
                if ov > 0: # ov[0] always one...
                    sub_scaling *= self.parent.overlap_correction_driver(ov).get_poly_evaluator(coefficient_tensors, freqs)
                for ec in ec:
                    # this should be cached...
                    sub_scaling *= self.parent.energy_correction_driver(ec).get_poly_evaluator(coefficient_tensors, freqs)

                # poly_scaling = poly_scaling + sub_scaling

                subpoly = subpoly * sub_scaling
                polys.append(subpoly) # this could be made more efficient by not having every exponent set act differently...

        raise Exception(polys)
        return self._poly_evaluator(polys)

class AnalyticCorrectionPolyEvaluator:
    def __init__(self, polys):
        self.polys = polys
    def eval(self, quanta):
        # this will need to be optimized and adapted
        # for resonance handling
        return sum(
            subp.evaluate(quanta)
            for subp in self.polys
        )

class AnalyticPTEnergyDriver(AnalyticPTOperatorElementDriver):

    def __init__(self, terms_inds, parent:AnalyticPerturbationTheoryDriver):
        super().__init__(terms_inds, [], parent)

    def enumerate_overlap_energy_exponent_tuples(self, index_tuple):
        #TODO: figure out parities

        remainder = self.order - sum(index_tuple)
        nterms = len(index_tuple)
        if remainder == 0: # just put a Pi between every term
            return [0], [[]], [[0] + [1]*(nterms - 1) + [0]]
        elif len(index_tuple) == 1: # only one term, so can only multiply by the overlap
            return [remainder], [[]], [[0, 0]]
        else: # need to figure out how to split up terms across the energies and overlaps
            # somehow we only ever have one overlap term so we'll just loop through those
            for ov in range(remainder+1):
                e_terms = []
                # enumerate all possible
                e_parts = IntegerPartitioner.partitions(remainder - ov)
                for part in e_parts:
                    engs, counts = np.unique(part, return_counts=True)
            raise NotImplementedError("still working out how to split integer partition permutations over E and Ov")


class AnalyticPTOverlapDriver(AnalyticPTOperatorElementDriver):

    def __init__(self, terms_inds, parent:AnalyticPerturbationTheoryDriver):
        super().__init__(terms_inds, [], parent)

    def enumerate_overlap_energy_exponent_tuples(self, index_tuple):
        #TODO: figure out parities
        raise Exception(...)

        remainder = self.order - sum(index_tuple)
        nterms = len(index_tuple)
        if remainder == 0: # just put a Pi between every term
            return [0], [[]], [[0] + [1]*(nterms - 1) + [0]]
        elif len(index_tuple) == 1: # only one term, so can only multiply by the overlap
            return [remainder], [[]], [[0, 0]]
        else: # need to figure out how to split up terms across the energies and overlaps
            # somehow we only ever have one overlap term so we'll just loop through those
            for ov in range(remainder+1):
                e_terms = []
                # enumerate all possible
                e_parts = IntegerPartitioner.partitions(remainder - ov)
                for part in e_parts:
                    engs, counts = np.unique(part, return_counts=True)
            raise NotImplementedError("still working out how to split integer partition permutations over E and Ov")

#endregion