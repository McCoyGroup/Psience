import abc
import itertools

import numpy as np, scipy.signal, time

from McUtils.Zachary import DensePolynomial, TensorCoefficientPoly
import McUtils.Numputils as nput
from McUtils.Combinatorics import SymmetricGroupGenerator, IntegerPartitioner, UniquePartitions, UniquePermutations
from McUtils.Scaffolding import Logger
from ..BasisReps import HarmonicOscillatorMatrixGenerator, HarmonicOscillatorRaisingLoweringPolyTerms

__all__ = [
    'AnalyticPerturbationTheorySolver',
    # 'AnalyticPerturbationTheoryDriver',
    # 'AnalyticPTCorrectionGenerator',
    # 'RaisingLoweringClasses'
]

class AnalyticPerturbationTheorySolver:
    """
    A re-attempt at using the recursive expressions
    to provide simpler code for getting APT expressions
    """
    def __init__(self, hamiltonian_expansion, logger=None):
        self.hamiltonian_expansion = hamiltonian_expansion
        self.logger = Logger.lookup(logger)

    @classmethod
    def from_order(cls, order, internals=True, logger=None):
        logger = Logger.lookup(logger)
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
                    identities=[0, 1] + ([] if o < 2 else [2]), # V, G, G'
                    logger=logger
                )
                for o in range(order + 1)
            ]
        else:
            # raise NotImplementedError("issues persist with Coriolis symmetries")
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
                            # ['p'] + ['x'] * (o-1) + ['p', 'x'],
                            ['x'] * (o - 2) # Watson
                        ]
                    ),
                    order=o,
                    identities=([0, 1] if o == 0 else [0]) + ([] if o < 2 else [3, 4]),  # V, G, G'
                    logger=logger
                )
                for o in range(order + 1)
            ]

        return cls(hamiltonian_terms, logger=logger)

    _op_maps = {}
    def get_correction(self, key, cls, order, **kw):
        k = (key, order)
        corr = self._op_maps.get(k, None)
        if corr is None:
            corr = cls(self, order, logger=self.logger, **kw)
            self._op_maps[k] = corr
        return corr

    def energy_correction(self, order, **kw):
        return self.get_correction('energy', EnergyCorrection, order, **kw)

    def wavefunction_correction(self, order, **kw):
        return self.get_correction('wfn', WavefunctionCorrection, order, **kw)

    def overlap_correction(self, order, **kw):
        return self.get_correction('ov', WavefunctionOverlapCorrection, order, **kw)

    def operator_correction(self, order, operator_type=None, **kw):
        return self.get_correction('op', OperatorCorrection, order, operator_type=operator_type, **kw)
        # return OperatorCorrection(self, order, operator_type=operator_type, logger=self.logger)

    @classmethod
    def operator_expansion_terms(cls, order, logger=None, operator_type=None):

        if operator_type is not None:
            raise NotImplementedError("only one type of expansion currently supported")

        return [
            OperatorExpansionTerm(
                [['x'] * (o + 1)],
                order=o,
                identities=[5],
                logger=logger
            )
            for o in range(order+1)
        ]

class ProductPTPolynomial:
    """
    TODO: include prefactor term so we can divide out energy changes
    """
    def __init__(self, coeffs, prefactor=1, idx=None):
        if (
                isinstance(coeffs, (int, float, np.integer, np.floating))
                or isinstance(coeffs[0], (int, float, np.integer, np.floating))
        ):
            raise ValueError("coeffs must be a vector of vectors (not {})".format(coeffs))
        self.coeffs = [np.asanyarray(c) for c in coeffs] # coeffs along each dim independently
        if any(c.dtype == np.dtype(object) for c in self.coeffs):
            raise ValueError(self.coeffs)
        self.prefactor = prefactor
        self._order = None
        self._monics = None
        self._idx = idx

    @property
    def order(self):
        if self._order is None:
            self._order = tuple(len(c)-1 for c in self.coeffs)
        return self._order
    def __repr__(self):
        return "{}(<{}>)".format("Poly", ",".join(str(c) for c in self.order))

    def permute(self, new_inds):
        return type(self)([self.coeffs[i] for i in new_inds], prefactor=self.prefactor)

    def constant_rescale(self):
        """
        rescales so constant term is 1

        :return:
        """
        new_prefactor = self.prefactor
        new_coeffs = []
        for c in self.coeffs:
            if abs(c[0]) > 1e-8: # non-zero
                new_prefactor *= c[0]
                new_coeffs.append(c/c[0])
            else:
                new_coeffs.append(c)
        return type(self)(new_coeffs, prefactor=new_prefactor)

    @property
    def monic_coeffs(self):
        if self._monics is None:
            self._monics = tuple(self._monify(c) for c in self.coeffs)
        return self._monics
    @monic_coeffs.setter
    def monic_coeffs(self, coeffs):
        self._monics = coeffs
    @staticmethod
    def _monify(coeffs, zero_thresh=1e-8):
        # print(coeffs)
        scaling = abs(coeffs[-1])
        if scaling <= zero_thresh:
            return np.zeros_like(coeffs), 0
        elif scaling == 1:
            return coeffs, 1
        else:
            return coeffs / scaling, scaling

        acoeffs = np.abs(coeffs) # np.round(np.abs(coeffs), 8)
        acmask = np.where(acoeffs > zero_thresh)
        if len(acmask) == 0 or len(acmask[0]) == 0:
            return np.zeros_like(coeffs), 0

        acmask = acmask[0]
        min_v = np.min(acoeffs[acmask])
        new_coeffs = np.zeros_like(coeffs)
        new_coeffs[acmask] = coeffs[acmask] / min_v
        # print(coeffs, new_coeffs)
        return new_coeffs, min_v

    @staticmethod
    def _find_off_pos(self_coeffs, other_coeffs):
        off_pos = None
        for n,(c1,c2) in enumerate(zip(self_coeffs, other_coeffs)):
            if len(c1) != len(c2):
                if off_pos is not None:
                    return False
                off_pos = n
        return off_pos
    def combine(self, other:'ProductPTPolynomial'):
        # if other in self._checked:
        #     ...
        # self._checked.add(other)

        if any(ms < 1e-8 for mc,ms in self.monic_coeffs):
            if any(ms < 1e-8 for mc,ms in other.monic_coeffs): return True
            return other
        elif any(ms < 1e-8 for mc,ms in other.monic_coeffs):
            return self

        off_pos = self._find_off_pos(self.coeffs, other.coeffs) # if a _single_ mode differs we can still combine
        # first pass off lengths b.c. most can't be simplified...
        if off_pos is False: return off_pos

        # second pass where we actually compare the coeffs
        condensed_prefactors = [self.prefactor, other.prefactor]
        new_coeffs = []
        for n,(c1,c2,m1,m2) in enumerate(zip(self.coeffs, other.coeffs, self.monic_coeffs, other.monic_coeffs)):

            if off_pos is not None and off_pos == n: # can handle _one_ difference
                if len(c2) > len(c1):
                    c1 = np.pad(c1, [0, len(c2)-len(c1)])
                if len(c1) > len(c2):
                    c2 = np.pad(c2, [0, len(c1)-len(c2)])
                new_coeffs.append([c1, c2])
                continue


            monic_coeffs1, monic_scaling1 = m1
            monic_coeffs2, monic_scaling2 = m2
            if monic_scaling1 < 1e-8:
                if monic_scaling2 < 1e-8: return True
                new_coeffs.append(c2)
            elif monic_scaling2 < 1e-8:
                new_coeffs.append(c1)
            elif np.any(monic_coeffs1 != monic_coeffs2):
                    if off_pos is None:  # can handle _one_ difference
                        off_pos = n
                        new_coeffs.append([c1, c2])
                        continue
                    else:
                        return False
            else:
                condensed_prefactors[0] *= monic_scaling1
                condensed_prefactors[1] *= monic_scaling2

                # scaling = new_prefactor.prefactor*monic_scaling1 + other.prefactor*monic_scaling2
                # if scaling == 0:
                #     return True # special case... where they cancel perfectly

                new_coeffs.append(monic_coeffs1)

        if off_pos is None: # no positions where coeffs differed
            new_prefactor = sum(condensed_prefactors)
            if new_prefactor == 0:
                return True
            new = ProductPTPolynomial(new_coeffs, prefactor=new_prefactor)
            new.monic_coeffs = self.monic_coeffs
        else:
            # gotta condense the position where things differed
            new_c = condensed_prefactors[0]*new_coeffs[off_pos][0] + condensed_prefactors[1]*new_coeffs[off_pos][1]
            monic_coeffs, monic_scaling = self._monify(new_c)
            if monic_scaling == 0:
                return True
            new_coeffs[off_pos] = monic_coeffs
            new = ProductPTPolynomial(new_coeffs, prefactor=monic_scaling)
            new_monics = [
                m if k != off_pos else (monic_coeffs, 1)
                for k,m in enumerate(self.monic_coeffs)
            ]
            new.monic_coeffs = new_monics
        return new

    def shift(self, shift):
        if self._idx is not None:
            base_shift = self._idx[-1]
            new_shift = tuple(s + k for s,k in zip(base_shift, shift)) + base_shift[len(shift):]
            new_idx = self._idx[:-1] + (new_shift,)
        else:
            new_idx = None
        return type(self)(
            [DensePolynomial.compute_shifted_coeffs(c, s) for c,s in zip(self.coeffs, shift)],
            prefactor=self.prefactor,
            idx=new_idx
        )

    def __mul__(self, other):
        if isinstance(other, (int, float, np.integer, np.floating)):
            return self.scale(other)
        else:
            raise NotImplementedError("subtle")
    def __rmul__(self, other):
        return self.__mul__(other)
    def scale(self, scalar):
        if nput.is_numeric(scalar):
            if scalar == 1:
                return self
            elif scalar == 0:
                return 0
        return type(self)(self.coeffs, prefactor=self.prefactor*scalar)

    def mul_simple(self, other:'ProductPTPolynomial'):
        if nput.is_numeric(other):
            if other == 1:
                return self
            elif other == 0:
                return 0
            else:
                raise ValueError(other)

        return self._poly_mul(self, other)


    _prod_poly_cache = {}
    @classmethod
    def _poly_mul(cls, self, other):

        bad_keys = other._idx is not None or self._idx is not None
        key = ((self._idx, other._idx), (0,) * len(self.coeffs)) if not bad_keys else None
        if bad_keys:
            new = None
        else:
            new = self._prod_poly_cache.get(key, None)

        if new is None:
            ocs = other.coeffs
            scs = self.coeffs
            if len(ocs) < len(scs):
                ocs = ocs + [[1]]*(len(scs) - len(ocs))
            elif len(scs) < len(ocs):
                scs = scs + [[1]]*(len(scs) - len(scs))

                # raise ValueError("not sure how to 'simply multiply' {} and {}".format(self, other))
            # print(scs, ocs)
            new = type(self)(
                [
                    scipy.signal.convolve(sc, oc)
                    for sc, oc in zip(scs, ocs)
                ],
                prefactor=self.prefactor * other.prefactor,
                idx=key
            )
            if not bad_keys: cls._prod_poly_cache[key] = new

        return new

    @classmethod
    def fast_ind_remainder(cls, n, diff):
        chunks = []
        s = 0
        for d in np.sort(diff):
            if d >= n:
                chunks.append(np.arange(s, n))
                break
            if s != d:
                chunks.append(np.arange(s, d))
            s = d + 1
        else:
            if s < n:
                chunks.append(np.arange(s, n))
        if len(chunks) > 1:
            chunks = np.concatenate(chunks)
        elif len(chunks) == 1:
            chunks = chunks[0]
        else:
            # print("???", n, diff)
            chunks = np.array([])
        # print("...", n, diff, np.setdiff1d(np.arange(n), diff), "->", chunks)
        return chunks
    def mul_along(self, other:'ProductPTPolynomial', inds, remainder=None, mapping=None):
        # print(">>>", 1)
        if isinstance(other, ProductPTPolynomialSum):
            new = type(other)(
                [
                    self.mul_along(o, inds, remainder=remainder, mapping=mapping)
                    for o in other.polys
                ],
                prefactor=self.prefactor * other.prefactor
            )
        elif isinstance(other, (PTTensorCoeffProductSum, PTEnergyChangeProductSum)):
            new = other.rmul_along(self, inds, remainder=None, mapping=mapping)
        else:
            if len(inds) == 0:
                # we get to just concatenate the coeffs
                new = type(self)(self.coeffs+other.coeffs, prefactor=self.prefactor*other.prefactor)
            else:
                if isinstance(inds[0], (int, np.integer)):
                    # we assume the inds line up with the inds for the other coeffs
                    inds = [np.arange(len(inds)), inds]

                    # remainder is just an optimization since we get it for free
                    if remainder is None:
                        remainder = self.fast_ind_remainder(len(self.coeffs), inds[1])
                    if len(remainder) == 0 or isinstance(remainder[0], (int, np.integer)):
                        remainder = [
                            np.arange(len(inds[1]), len(other.coeffs)),
                            remainder
                            ]

                other_inds, self_inds = inds # axes to convolve
                if remainder is None:
                    remainder = [
                        self.fast_ind_remainder(len(other.coeffs), other_inds),
                        self.fast_ind_remainder(len(self.coeffs), self_inds)
                    ]
                other_remainder, self_remainder = remainder

                # self_remapping = np.argsort(np.concatenate([self_inds, self_remainder])) # literally just the inverse...(I could pass this in)
                # other_remapping = np.argsort(np.concatenate([self_inds, self_remainder])) # should avoid actually doing this sort (when `remainder` is/was `None)

                new_coeffs = [
                                 scipy.signal.convolve(other.coeffs[ox], self.coeffs[sx])
                                 #     self._remap_tensor_coefficients(other.coeffs[ox], other_remapping),
                                 #     self._remap_tensor_coefficients(self.coeffs[sx], self_remapping)
                                 # )
                                 for ox, sx in zip(other_inds, self_inds)
                                 if len(other.coeffs) > ox and len(self.coeffs) > sx
                             ] + [
                                 self.coeffs[sx] for sx in self_remainder
                                 if len(self.coeffs) > sx
                             ] + [
                                 other.coeffs[ox] for ox in other_remainder
                                 if len(other.coeffs) > ox
                             ]


                new = type(self)(new_coeffs, prefactor=self.prefactor*other.prefactor)

        # print("<<<", 1)
        return new

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
                raise NotImplementedError(...)
                return self + ProductPTPolynomial([[other]])
        elif isinstance(other, ProductPTPolynomial):
            # if other.prefactor != self.prefactor:
            #     print("????", other.prefactor, self.prefactor)
            #     other = other.scale(1/self.prefactor)
            #     # raise NotImplementedError("need to include prefactors")
            return ProductPTPolynomialSum([self, other])#, prefactor=self.prefactor)
        elif isinstance(other, ProductPTPolynomialSum):
            return other + self
        else:
            raise NotImplementedError("not sure how to add {} and {}".format(self, other))
    def __radd__(self, other):
        return self + other

class ProductPTPolynomialSum:

    def __init__(self, polynomials, prefactor=1, reduced=False):
        self.polys = polynomials
        self.prefactor = prefactor
        self.reduced = reduced
        self._order = None
        for p in self.polys:
            if p.order == 0: raise Exception(p.coeffs)

    def __repr__(self):
        form_counts = {}
        for p in self.polys:
            key = "<{}>".format(",".join(str(c) for c in p.order))
            form_counts[key] = form_counts.get(key, 0) + 1
        return "PSum({})".format(
            # type(self).__name__,
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

    def permute(self, new_inds):
        return type(self)([p.permute(new_inds) for p in self.polys],
                          prefactor=self.prefactor, reduced=self.reduced)

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
        if self.reduced: return self
        polys = self.polys
        cur_len = len(polys) + 1
        cache = set()
        while cur_len > len(polys): # hit a fixed point with this transformation
            cur_len = len(polys)
            polys = self.combine_polys(polys, cache)
        return type(self)(polys, prefactor=self.prefactor, reduced=True)

    def shift(self, shift):
        return type(self)([p.shift(shift) for p in self.polys],
                          prefactor=self.prefactor, reduced=self.reduced)

    def __mul__(self, other):
        if isinstance(other, (int, float, np.integer, np.floating)):
            return self.scale(other)
        else:
            raise NotImplementedError("subtle")
    def __rmul__(self, other):
        return self.__mul__(other)
    def scale(self, scalar):
        if isinstance(scalar, (int, float, np.integer, np.floating)):
            if scalar == 1:
                return self
            elif scalar == 0:
                return 0
        return type(self)(self.polys, prefactor=self.prefactor*scalar)

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
        elif nput.is_numeric(other):
            if other == 1:
                return self
            elif other == 0:
                return 0
            else:
                raise ValueError(other)
        else:
            raise NotImplementedError(type(other))

    def mul_along(self, other:'ProductPTPolynomial', inds, remainder=None, mapping=None):
        # print(">>>", 2)
        if isinstance(other, ProductPTPolynomial):
            # can just distribute except at some level I have to assume
            # that every stored poly has the same dimension...?
            new = type(self)(
                [
                    poly.mul_along(other, inds, remainder=remainder, mapping=mapping) # how is this supposed to work?
                    for poly in self.polys # the fuck is inds supposed to mean here...????
                ],
                prefactor=self.prefactor
            )
        elif isinstance(other, ProductPTPolynomialSum):
            new = type(self)(
                [
                    p1.mul_along(p2, inds, remainder=remainder, mapping=mapping)
                    for p1 in self.polys
                    for p2 in other.polys
                ],
                prefactor=self.prefactor*other.prefactor
            )
        elif isinstance(other, (PTTensorCoeffProductSum, PTEnergyChangeProductSum)):
            new = other.rmul_along(self, inds, remainder=remainder, mapping=mapping)
        else:
            raise NotImplementedError(type(other))
        # print("<<<", 2)
        return new


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
                ops = other.polys
                if other.prefactor != self.prefactor:
                    ops = [p.scale(other.prefactor/self.prefactor) for p in ops]
                return type(self)(self.polys + ops, prefactor=self.prefactor)
            elif isinstance(other, ProductPTPolynomial):
                if self.prefactor != 1:
                    other = other.scale(1/self.prefactor)
                    # raise NotImplementedError("need to handle prefactors")
                return type(self)(self.polys + [other], prefactor=self.prefactor)
            else:
                raise TypeError("not sure what to do with {} and {}".format(type(self), type(other)))
    def __radd__(self, other):
        return self + other

class PTEnergyChangeProductSum(TensorCoefficientPoly):
    """
    A representation of a sum of 1/energy * poly sums
    which is here so we can transpose energy change indices intelligently
    """
    def __init__(self, terms: dict, prefactor=1, canonicalize=True, reduced=False):
        super().__init__(terms, prefactor=prefactor, canonicalize=canonicalize)
        self.reduced = reduced

    @classmethod
    def format_key(self, key):
        return "".join([
                "E-[{}]".format(k)
                for k in key
            ])
    def __repr__(self):
        sums = []
        for k_prod,v in self.terms.items():
            sums.append("{}{}".format(v, self.format_key(k_prod)))
        return "ESum({})".format("+".join(sums) if len(sums) < 3 else "\n      +".join(sums))

    def shift(self, shift):
        return type(self)(
            {
                k: p.shift(shift) for k,p in self.terms.items()
            },
            prefactor=self.prefactor
        )

    def sort(self):
        return type(self)(
            {
                k: self.terms[k] for k in sorted(self.terms.keys())
            },
            prefactor=self.prefactor
        )

    def _permute_changes(self, changes, new_inds):
        new = tuple(
            changes[i]
                if i < len(changes) else
            0
            for i in new_inds
        )
        for j in range(1, len(new)+1):
            if new[-j] != 0:
                break
        else:
            j = 1 # just to shut the IDE up...
        new = new[:len(new) - (j-1)]
        return new
    def permute(self, new_inds):
        new_terms = {}
        for energy_changes,polys in self.terms.items():
            new_ech = tuple(
                self._permute_changes(ec, new_inds)
                for ec in energy_changes
            )
            new_terms[new_ech] = polys.permute(new_inds)
        return type(self)(new_terms, prefactor=self.prefactor)

    @staticmethod
    def _check_neg(t1, t2):
        if len(t1) != len(t2):
            return False
        else:
            return all(
                all(ttt1 == -ttt2 for ttt1, ttt2 in zip(tt1, tt2))
                for tt1, tt2 in zip(t1, t2)
            )

    def combine_energies(self):
        base_keys = list(self.terms.keys())
        elim_pos = set()
        merges = set()
        unmerged = set()
        for i,k1 in enumerate(base_keys):
            if i not in elim_pos:
                for j,k2 in enumerate(base_keys[i+1:]):
                    j = i+1 + j
                    if j not in elim_pos:
                        if self._check_neg(k1, k2):
                            elim_pos.add(i)
                            elim_pos.add(j)
                            merges.add((k1, k2))
                            break
                else:
                    unmerged.add(k1)

        new_terms = {}
        for k1,k2 in merges:
            # sanity check
            if k1 in unmerged or k2 in unmerged:
                raise ValueError('fuck')
            # print(self.terms[k2].prefactor, self.terms[k2].scale(-1))
            new_terms[k1] = self.terms[k1] + self.terms[k2].scale(-1)
        for k in unmerged:
            new_terms[k] = self.terms[k]

        return new_terms

    def combine(self, combine_subterms=True, combine_energies=True):
        if self.reduced: return self

        new_terms = {}
        # raise Exception(
        #     list(self.terms.keys()),
        #     list(self.combine_energies().keys())
        # )
        base_terms = self.combine_energies() if combine_energies else self.terms
        for k,p in base_terms.items():
            if combine_subterms and isinstance(p, ProductPTPolynomialSum):
                p = p.combine()
                if len(p.polys) > 0:
                    new_terms[k] = p
            else:
                new_terms[k] = p
        return type(self)(new_terms, prefactor=self.prefactor, reduced=True)

    @staticmethod
    def _permute_idx(idx, inds, mapping=None):
        if mapping is not None:
            return tuple(mapping.get(k, k) for k in idx)
        else:
            # print(inds, idx, len(inds))
            return tuple(idx[k] if k < len(idx) else k for k in inds)
    def mul_along(self, other:'ProductPTPolynomial', inds, remainder=None, mapping=None):
        """
        We multiply every subpoly along the given indices, transposing the appropriate tensor indices

        :param other:
        :param inds:
        :param remainder:
        :return:
        """
        # print(">>>", 3)
        if isinstance(other, PTEnergyChangeProductSum):
            # print("...", 1, list(self.terms.keys()))
            # print("...", 2, list(other.terms.keys()))
            # print("idx", inds)
            # if not nput.is_numeric(inds[0]):
            #     raise ValueError(...)
            new_poly = self.direct_product(other,
                                       key_func=lambda k1,k2,inds=inds:k1+tuple(self._permute_idx(idx, inds, mapping=mapping) for idx in k2),
                                       mul=lambda a,b,inds=inds:a.mul_along(b, inds, remainder=remainder, mapping=mapping)
                                       )
        # elif isinstance(other, PTTensorCoeffProductSum):
        #     return other.mul_along(self, inds, remainder=remainder)
        elif isinstance(other, (ProductPTPolynomial, ProductPTPolynomialSum)):
            new_poly = type(self)(
                {
                    k: (
                        other.scale(v)
                            if isinstance(v, (int, float, np.integer, np.floating)) else
                        v.mul_along(other, inds, remainder=remainder, mapping=mapping)
                    )
                    for k, v in self.terms.items()
                },
                prefactor=self.prefactor
            )
        else:
            raise NotImplementedError(type(other))
        # print("<<<", 3)
        return new_poly

    def rmul_along(self, other, inds, remainder=None, mapping=None):
        """
        We multiply every subpoly along the given indices, transposing the appropriate tensor indices

        :param other:
        :param inds:
        :param remainder:
        :return:
        """
        if isinstance(other, PTEnergyChangeProductSum):
            return other.mul_along(self, inds, remainder=remainder, mapping=mapping)
        # elif isinstance(other, PTTensorCoeffProductSum):
        #     return other.mul_along(self, inds, remainder=remainder)
        elif isinstance(other, (ProductPTPolynomial, ProductPTPolynomialSum)):
            return type(self)(
                {
                    k: (
                        other.scale(v)
                            if isinstance(v, (int, float, np.integer, np.floating)) else
                        other.mul_along(v, inds, remainder=remainder, mapping=mapping)
                    )
                    for k, v in self.terms.items()
                },
                prefactor=self.prefactor
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
        if isinstance(other, (PTTensorCoeffProductSum, PTEnergyChangeProductSum)):
            return self.direct_product(other,
                                       # key_func=lambda k1,k2:k1+tuple(k2[i] for i in inds), #inds says how k2 permutes?
                                       mul=lambda a,b:a.mul_simple(b)
                                       )
        elif isinstance(other, (ProductPTPolynomial, ProductPTPolynomialSum)):
            return type(self)(
                {
                    k: other.scale(v) if isinstance(v, (int, float, np.integer, np.floating)) else v.mul_simple(other)
                    for k, v in self.terms.items()
                },
                prefactor=self.prefactor
            )
        elif nput.is_numeric(other):
            if other == 1:
                return self
            elif other == 0:
                return 0
            else:
                raise ValueError(other)
        else:
            raise NotImplementedError(type(other))

class PTTensorCoeffProductSum(TensorCoefficientPoly):
    """
    A representation for a sum of tensor coefficients * polynomial sums
    which is primarily here so we can transpose tensor coefficients intelligently
    """

    def __init__(self, terms, prefactor=1, canonicalize=True, inds_map=None, reduced=False):
        super().__init__(terms, prefactor=prefactor, canonicalize=canonicalize)
        self._inds_map = {} if inds_map is None else inds_map
        self.reduced = reduced

    @classmethod
    def format_key(self, key):
        return "".join([
            "{}[{}]{}".format(
                ["V", "G", "V'", "Z", "U", "M"][k[1]],  # we explicitly split Coriolis and Watson terms out
                k[0],
                k[2:]
            )
            for k in key
        ])
    def __repr__(self):
        sums = []
        for k_prod,v in self.terms.items():
            sums.append("{}{}".format(v, self.format_key(k_prod)))
        return "TSum({})".format("+".join(sums) if len(sums) < 3 else "\n      +".join(sums))

    def print_tree(self):
        for k,esum in self.terms.items():
            print("="*20, self.format_key(k), self.prefactor, "="*20)
            if isinstance(esum, PTEnergyChangeProductSum):
                for e,pp in esum.terms.items():
                    print("  ", "-"*10, esum.format_key(e), esum.prefactor, "-"*10)
                    if isinstance(pp, ProductPTPolynomialSum):
                        print("  ::", pp.prefactor)
                        for p in pp.polys:
                            p = p.constant_rescale()
                            print("   >", p.prefactor, p.coeffs)
            elif isinstance(esum, ProductPTPolynomialSum):
                print("  ::", esum.prefactor)
                for p in esum.polys:
                    p = p.constant_rescale()
                    print("   >", p.prefactor, p.coeffs)
            else:
                p = esum.constant_rescale()
                print("   >", p.prefactor, p.coeffs)

    @classmethod
    def coeff_product_inds(cls, key, return_counts=False):
        return np.unique(np.concatenate([c[2:] for c in key]), return_counts=return_counts)
    def get_inds(self, key):
        if key not in self._inds_map:
            inds = self.coeff_product_inds(key)
            self._inds_map[key] = inds
        return self._inds_map[key]

    def prune_operators(self, ops):
        ops = set(ops)
        return type(self)(
            {
                t: p
                for t, p in self.terms.items()
                if all(tt[:2] not in ops for tt in t)
            },
            prefactor=self.prefactor
        )

    def _audit_dimension(self, poly, num_inds):
        if isinstance(poly, ProductPTPolynomial):
            if len(poly.coeffs) != num_inds:
                if not (
                        len(poly.coeffs) == 1 and
                        len(poly.coeffs[0]) == 1 and
                        num_inds == 0
                ):  # ignore constants
                    raise ValueError("{} is fucked (expected {})".format(
                        poly,
                        num_inds
                    ))
        elif isinstance(poly, ProductPTPolynomialSum):
            for p in poly.polys:
                self._audit_dimension(p, num_inds)
        elif isinstance(poly, PTEnergyChangeProductSum):
            for p in poly.terms.values():
                self._audit_dimension(p, num_inds)
        else:
            raise ValueError('wat')
    def audit(self):
        """
        Checks to ensure that the number of dimensions aligns with
        the number of indices in the tensor coefficients

        :return:
        """
        for coeff_indices,poly_stuff in self.terms.items():
            num_inds = len(self.get_inds(coeff_indices))
            self._audit_dimension(poly_stuff, num_inds)

    def sort(self):
        return type(self)(
            {
                k: self.terms[k] for k in sorted(self.terms.keys())
            },
            prefactor=self.prefactor
        )

    def permute(self, new_inds):
        new_terms = {}
        for coeff_inds,polys in self.terms.items():
            new_inds = tuple(ci[:2] + tuple(new_inds[j] for j in ci[2:]) for ci in coeff_inds)
            new_terms[new_inds] = polys.permute(new_inds)
        return type(self)(new_terms, prefactor=self.prefactor, inds_map=self._inds_map)

    def _check_equiv(self, k1, k2):
        if len(k1) != len(k2):
            return None
        if any(len(kk1) != len(kk2) or kk1[:2] != kk2[:2] for kk1, kk2 in zip(k1, k2)):
            return None
        k1_inds = self.get_inds(k1)
        k2_inds = self.get_inds(k2)
        if len(k1_inds) != len(k2_inds) or any(ki1 != ki2 for ki1, ki2 in zip(k1_inds, k2_inds)):
            return None

        for p in itertools.permutations(range(len(k1_inds))):
            swap_inds = [k2_inds[i] for i in p]
            k2_sub = tuple(ci[:2] + tuple(swap_inds[j] for j in ci[2:]) for ci in k2)
            if k1 == k2_sub:
                return p # permutation to apply to our polys

    def combine_terms(self):
        base_keys = list(self.terms.keys())
        elim_pos = set()
        merges = {}
        unmerged = set()
        for i,k1 in enumerate(base_keys):
            if i not in elim_pos:
                for j,k2 in enumerate(base_keys[i+1:]):
                    j = i+1 + j
                    if j not in elim_pos:
                        remapping = self._check_equiv(k1, k2)
                        if remapping is not None:
                            elim_pos.add(i)
                            elim_pos.add(j)
                            merges[(k1, k2)] = remapping
                            break
                else:
                    unmerged.add(k1)

        new_terms = {}
        for k1,k2 in merges:
            # sanity check
            if k1 in unmerged or k2 in unmerged:
                raise ValueError('fuck')
            # print(self.terms[k2].prefactor, self.terms[k2].scale(-1))
            new_terms[k1] = self.terms[k1] + self.terms[k2].permute(merges[(k1, k2)])
        for k in unmerged:
            new_terms[k] = self.terms[k]

        return new_terms

    def combine(self, combine_coeffs=False, combine_subterms=True, combine_energies=True):
        if self.reduced: return self
        new_terms = {}

        if combine_coeffs:
            raise NotImplementedError("coeff combining too subtle for current evaluator strategy")

        base_terms = self.combine_terms() if combine_coeffs else self.terms
        for k,p in base_terms.items():
            if combine_subterms and isinstance(p, ProductPTPolynomialSum):
                p = p.combine()
                if len(p.polys) > 0:
                    new_terms[k] = p
            elif isinstance(p, (PTTensorCoeffProductSum, PTEnergyChangeProductSum)):
                p = p.combine(combine_subterms=combine_subterms, combine_energies=combine_energies)
                if len(p.terms) > 0:
                    new_terms[k] = p
            else:
                new_terms[k] = p
        return type(self)(new_terms, prefactor=self.prefactor, reduced=True)

    def shift(self, shift):
        new_terms = {}
        for k,p in self.terms.items():
            new_terms[k] = p.shift(shift)
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
        # symmetry of the Coriolis term is obviously weird since it's not even Hermitian???
        if idx[1] > idx[-1] or (idx[1] == idx[-1] and idx[0] > idx[2]): # only case where we can flip stuff around...
            rest = idx[-2:1:-1] + (idx[0],)
            return (rest[0], idx[-1]) + rest[1:] + (idx[1],)
        else:
            return idx
        # p_sorts = list(sorted([idx[1], idx[-1]]))
        # rest_sorts = list(sorted((idx[0],) + idx[2:-1]))
        # new = (rest_sorts[0], p_sorts[0]) + tuple(rest_sorts[1:]) + (p_sorts[1],)
        # return new

    @classmethod
    def symmetrizers(cls):
        return (
            cls._potential_symmetrizer,
            cls._kinetic_symmetrizer,
            cls._potential_symmetrizer,  # V'
            cls._coriolis_symmetrizer,
            cls._potential_symmetrizer,  # U
            cls._potential_symmetrizer,  # M
        )
    @classmethod
    def _symmetrize(cls, idx):
        return (idx[0], idx[1]) + cls.symmetrizers()[idx[1]](idx[2:])
    @classmethod
    def canonical_key(cls, monomial_tuple):
        return super().canonical_key(tuple(cls._symmetrize(t) for t in monomial_tuple))

    def _generate_direct_product_values(self, inds, k1, k2, poly_1, poly_2):
        # print("="*19)
        inds = tuple(inds)

        left_inds = self.get_inds(k1)#np.unique(np.concatenate([c[2:] for c in k1]))
        right_inds = self.get_inds(k2)#np.unique(np.concatenate([c[2:] for c in k2]))
        # num_inds = len(np.unique(np.concatenate([left_inds, right_inds]))) # how many inds do we have total
        num_left = len(left_inds)
        num_right = len(right_inds)
        num_fixed = len(inds)

        left_fixed_inds = tuple(np.arange(num_fixed))
        right_fixed_inds = inds
        left_remainder_inds = np.arange(num_fixed, num_left)
        right_remainder_inds = ProductPTPolynomial.fast_ind_remainder(num_right, inds)

        max_mult = min(num_right - num_fixed, num_left - num_fixed) # the max number of extra axes to align
        for mul_size in range(max_mult+1):
            if mul_size == 0:
                new_poly = poly_1.mul_along(poly_2, inds)
                left_mapping = {k: i for i, k in enumerate(left_fixed_inds)}
                left_rem_rem = left_remainder_inds
                for j, k in enumerate(left_rem_rem):
                    left_mapping[k] = j + num_fixed + mul_size
                right_mapping = {k: i for i, k in enumerate(right_fixed_inds)}
                right_rem_rem = right_remainder_inds
                for j, k in enumerate(right_rem_rem):
                    right_mapping[k] = j + num_left

                new_key = tuple(
                    t[:2] + tuple(left_mapping[k] for k in t[2:])
                    for t in k1
                ) + tuple(
                    t[:2] + tuple(right_mapping[k] for k in t[2:])
                    for t in k2
                )
                # print(">>>", k1, k2, inds, new_key)
                yield new_key, new_poly
            else:
                for left_choice_x in itertools.combinations(range(len(left_remainder_inds)), r=mul_size):
                # for left_choice in itertools.combinations(left_remainder_inds, r=mul_size):
                    left_choice = tuple(left_remainder_inds[k] for k in left_choice_x)
                    for right_choice_x in itertools.combinations(range(len(right_remainder_inds)), r=mul_size):
                        # right_choice = right_remainder_inds[right_choice_x,]
                        right_choice = tuple(right_remainder_inds[k] for k in right_choice_x)
                        for right_perm in itertools.permutations(right_choice):
                            mul_inds = [
                                left_fixed_inds + left_choice, # left
                                right_fixed_inds + right_perm  # right
                            ]
                            new_poly = poly_1.mul_along(poly_2, mul_inds, mapping=right_mapping)
                            # now we build our new key, noting that the multiplied inds will be the first _n_
                            # inds in the new set, then the left inds, then the right ones
                            left_mapping = {k:i for i,k in enumerate(mul_inds[0])}
                            left_rem_rem = np.delete(left_remainder_inds, left_choice_x)
                            # left_rem_rem = np.setdiff1d(left_remainder_inds, left_choice)
                            for j,k in enumerate(left_rem_rem):
                                left_mapping[k] = j + num_fixed + mul_size
                            right_mapping = {k:i for i,k in enumerate(mul_inds[1])}
                            right_rem_rem = np.delete(right_remainder_inds, right_choice_x)
                            for j,k in enumerate(right_rem_rem):
                                right_mapping[k] = j + num_left

                            new_key = tuple(
                                t[:2] + tuple(left_mapping[k] for k in t[2:])
                                for t in k1
                            ) + tuple(
                                t[:2] + tuple(right_mapping[k] for k in t[2:])
                                for t in k2
                            )
                            # print(">2>", k1, k2, inds, new_key)
                            yield new_key, new_poly

    def mul_along(self, other:'ProductPTPolynomial', inds, remainder=None, mapping=None):
        """
        We multiply every subpoly along the given indices, transposing the appropriate tensor indices

        :param other:
        :param inds:
        :param remainder:
        :return:
        """
        # print(">>>", inds)
        if isinstance(other, PTTensorCoeffProductSum):
            # we need to take the direct product of aligning (or not)
            # for all of the untouched inds
            new = self.direct_multiproduct(other, lambda *kargs: self._generate_direct_product_values(inds, *kargs))
            # new.audit()
            if self._inds_map is not other._inds_map:
                self._inds_map.update(other._inds_map) # these could all be shared for even more speed up
            new._inds_map = self._inds_map
            # try:
            # new.audit()
            # except ValueError:
            #     raise ValueError("direct product of tensors of dimensions {} and {} gave invalid result".format(
            #         self, other
            #     ))
            # return new
        elif isinstance(other, PTEnergyChangeProductSum):
            new = type(self)(
                {k:other.rmul_along(p, inds, remainder=remainder, mapping=mapping) for k,p in self.terms.items()},
                prefactor=self.prefactor,
                inds_map=self._inds_map
            )
        else:
            raise NotImplementedError(type(other))
        # print("<<<", 4)
        return new
    def rmul_along(self, other:'ProductPTPolynomial', inds, remainder=None, mapping=None):
        """
        We multiply every subpoly along the given indices, transposing the appropriate tensor indices

        :param other:
        :param inds:
        :param remainder:
        :return:
        """
        if isinstance(other, PTTensorCoeffProductSum):
            new = other.mul_along(self, inds, remainder=remainder, mapping=mapping)
            if self._inds_map is not other._inds_map:
                self._inds_map.update(other._inds_map) # these could all be shared for even more speed up
            new._inds_map = self._inds_map
            return new
        elif isinstance(other, PTEnergyChangeProductSum):
            return type(self)(
                {k:other.mul_along(p, inds, remainder=remainder, mapping=mapping) for k,p in self.terms.items()},
                prefactor=self.prefactor,
                inds_map = self._inds_map
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
        elif isinstance(other, PTEnergyChangeProductSum):
            return type(self)(
                {
                    k: other.mul_simple(v)
                    for k, v in self.terms.items()
                },
                prefactor=self.prefactor
            )
        elif isinstance(other, (ProductPTPolynomial, ProductPTPolynomialSum)):
            return type(self)(
                {
                    k: v.mul_simple(other)
                    for k, v in self.terms.items()
                },
                prefactor=self.prefactor
            )

class PerturbationTheoryTerm(metaclass=abc.ABCMeta):
    """
    A generic version of one of the three terms in
    PT that will generate a correction polynomial
    """
    def __init__(self, logger=None):
        self._exprs = None
        self._raw_changes = {}
        self._changes = None
        self._cache = {}
        self.logger = Logger.lookup(logger)

    def get_subexpresions(self) -> 'Iterable[PerturbationTheoryTerm]':
        raise NotImplementedError("just here to be overloaded")
    def __mul__(self, other):
        if isinstance(other, PerturbationTheoryTerm):
            return PerturbationTheoryTermProduct.lookup(self, other)
        else:
            return ScaledPerturbationTheoryTerm.lookup(self, other)
    def __rmul__(self, other):
        # if isinstance(other, PerturbationTheoryTerm): # can't happen
        #     return ...
        return ScaledPerturbationTheoryTerm.lookup(self, other)
    def __neg__(self):
        return ScaledPerturbationTheoryTerm.lookup(self, -1)

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


    def get_poly_terms(self, changes, simplify=True, shift=None) -> 'PTTensorCoeffProductSum':
        with self.logger.block(tag="Terms for {}[{}]".format(self, changes)):
            start = time.time()
            changes = tuple(changes)
            terms = self._raw_changes.get(changes, None)
            if not isinstance(terms, (PTTensorCoeffProductSum, int, np.integer)):
                terms = sum(
                    term.get_poly_terms(changes, shift=None)
                    for term in self.expressions
                )
                self._raw_changes[changes] = terms
            if simplify and not nput.is_zero(terms):
                # print("...simplifying...", end = None)
                terms = terms.combine().sort()
                # print("done")
            if simplify:
                self.changes[changes] = terms
            if shift is not None and len(shift) > 0 and not nput.is_numeric(terms):
                # if len(shift) == 0:
                #     print("fuck")
                #     raise Exception(...)
                terms = terms.shift(shift)
                # print("  .", [type(t) for t in terms.terms.keys()])
                # print(terms.terms)
                # shift = tuple(shift)
                # t = self._cache.get((changes, shift), None)
                # if t is None:
                #     t = terms.shift(shift)
                #     self._cache[(changes, shift)] = t
                # terms = t
            end = time.time()
            self.logger.log_print('took {e:.3f}s', e=end-start)
        return terms

    def __call__(self, changes, shift=None, coeffs=None, freqs=None, check_sorting=True, simplify=True, return_evaluator=True):
        with self.logger.block(tag="Building evaluator {}[{}]".format(self, changes)):
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
                terms = self.get_poly_terms(changes, shift=None)
                if simplify and not nput.is_zero(terms):
                    # print("...simplifying...", end = None)
                    terms = terms.combine().sort()
                    # print("done")
                if simplify:
                    self.changes[changes] = terms
            # terms = self.changes[changes]
            if shift is not None:
                raise NotImplementedError("reshifting not supported yet...")
            if return_evaluator:
                return PerturbationTheoryEvaluator(terms, changes)
            else:
                return terms

class OperatorExpansionTerm(PerturbationTheoryTerm):
    def __init__(self, terms, order=None, identities=None, symmetrizers=None, logger=None):
        super().__init__(logger=logger)

        self.terms = terms
        self.order = order
        self.identities = np.arange(len(self.terms)) if identities is None else identities
        if symmetrizers is None:
            symmetrizers = [PTTensorCoeffProductSum._potential_symmetrizer] * (max(self.identities) + 1)
        self.symmetrizers = symmetrizers

        self._grouped_terms = {}
        for i,t in enumerate(self.terms):
            self._grouped_terms[len(t)] = self._grouped_terms.get(len(t), [])
            self._grouped_terms[len(t)].append([i, tuple(t)])

        self.term_sizes = set(self._grouped_terms.keys())
        self._change_poly_cache = {}
        self._cache = {}

    def __repr__(self):
        return "M[{}]".format(self.order)

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
        with self.logger.block(tag="Terms for {}[{}]".format(self, changes)):
            start = time.time()
            key = (tuple(changes), tuple(shift) if shift is not None else None)
            terms = self._raw_changes.get(key, None)
            if terms is None:
                terms = self._get_poly_terms(changes, shift=shift)
                self._raw_changes[key] = terms
            end = time.time()
            self.logger.log_print("took {e:.3f}s", e=end-start)
        return terms

    def _get_poly_terms(self, changes, shift=None) -> 'PTTensorCoeffProductSum':  # TODO: CACHE THIS SHIT

        poly_contribs = {}

        abs_change = np.array([abs(c) for c in changes])
        total_change = sum(abs_change)
        for group_size, terms in self._grouped_terms.items():
            remainder = group_size - total_change
            if remainder < 0 or remainder % 2 != 0:
                continue # literally cannot execute this change

            # not sure what to do if group_size == 0...?
            if group_size == 0: # constant contrib
                for term_index, term_list in terms:
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
                if len(shift) < total_dim:
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

            # now we drop all targets with zeros in front of non-zero terms b.c. these
            # terms should reduce out
            neg_pos = np.where(np.diff(targets, axis=1) > 0)
            zero_pos = np.where(targets[neg_pos] == 0)
            if len(zero_pos) > 0:
                zero_pos = zero_pos[0]
                if len(zero_pos) > 0:
                    bad_targs = neg_pos[0][zero_pos]
                    targets = targets[ProductPTPolynomial.fast_ind_remainder(len(targets), bad_targs)]
            # print(zero_pos, np.diff(targets, axis=1))

            # then for each of these partition sizes, we enumerate the actual partitions of the terms
            # and then take the corresponding terms to generate partitions
            # (and multiply by the appropriate tensor coefficients)

            partitioner = UniquePartitions(np.arange(group_size))
            ckey = tuple(changes)
            skey = tuple(shift)
            # Integer
            for partition_sizes in targets:
                parts, inv_splits = partitioner.partitions(partition_sizes,
                                                          return_partitions=False, take_unique=False,
                                                          return_indices=True, return_inverse=True
                                                          )

                inv_splits = np.concatenate(inv_splits, axis=1)
                partition_sizes = tuple(partition_sizes)
                base_tup = np.array(sum(([i]*s for i,s in enumerate(partition_sizes)), []))
                for j,p_vec in enumerate(zip(*parts)):
                    tensor_idx = tuple(base_tup[inv_splits[j]])
                    for term_index,term_list in terms:
                        base_poly = self._resolve_poly(term_list, partition_sizes, j, ckey, skey, p_vec, changes, shift)
                        if nput.is_zero(base_poly): continue

                        term_id = self.identities[term_index]
                        if self.symmetrizers[term_id] is not None:
                            symm_idx = self.symmetrizers[term_id](tensor_idx)
                        else:
                            symm_idx = tensor_idx
                        prefactor = (self.order, term_id) + symm_idx #type: tuple[int]
                        poly_contribs[(prefactor,)] = poly_contribs.get((prefactor,), 0) + base_poly

        new = PTTensorCoeffProductSum(poly_contribs, canonicalize=False)
        new.audit()
        return new

    _poly_cache = {}
    @classmethod
    def _resolve_poly(cls, term_list, partition_sizes, p_index, c_key, s_key, p_vec, changes, shift):
        # print(term_list, p_index, changes, shift, p_index)

        key = (term_list, partition_sizes, p_index, c_key, s_key)
        poly = cls._poly_cache.get(key, None)
        if poly is None:
            poly_coeffs = [
                cls._evaluate_poly_coeffs(term_list, inds, delta, s)
                for inds, delta, s in zip(p_vec, changes, shift)
            ]
            if any(nput.is_zero(c) for c in poly_coeffs):
                poly = 0
            else:
                poly_coeffs = [
                    c for c in poly_coeffs
                    if c[0] != 1 or len(c) > 1
                ]
                poly = ProductPTPolynomial(poly_coeffs, idx=key)
            cls._poly_cache[key] = poly
        return poly


    # def get_correction_poly(self, changes, shifts=None):
    #     #TODO: check to make sure changes is sorted??? Or not???
    #     t = tuple(changes)
    #     if t not in self._change_poly_cache:
    #         self._change_poly_cache[t] = FullChangePathPolyTerm(
    #             self.get_poly_terms(changes, shifts)
    #         )
    #     return self._change_poly_cache[t]

class HamiltonianExpansionTerm(OperatorExpansionTerm):
    def __init__(self, terms, order=None, identities=None, symmetrizers=None, logger=None):

        if identities is None:
            identities = np.arange(5) if identities is None else identities
        if symmetrizers is None:
            symmetrizers = PTTensorCoeffProductSum.symmetrizers()

        super().__init__(terms, order=order, identities=identities, symmetrizers=symmetrizers, logger=logger)

    def __repr__(self):
        return "H[{}]".format(self.order)

class PerturbationOperator(PerturbationTheoryTerm):
    def __init__(self, subterm):
        super().__init__(logger=subterm.logger)

        self.subterm = subterm

    _cache = {}

    @classmethod
    def lookup(cls, subterm):
        prod = cls._cache.get(subterm, None)
        if prod is None:
            prod = cls(subterm)
            cls._cache[subterm] = prod
        return prod

    def __repr__(self):
        return "P{}".format(self.subterm)

    def get_changes(self) -> 'dict[tuple[int], Any]':
        return {k:None for k in self.subterm.changes if k != ()}

    def get_poly_terms(self, changes, shift=None) -> 'PTTensorCoeffProductSum':
        if shift is not None:
            final_change = tuple(c + s for c, s in zip(changes, shift))
        else:
            final_change = tuple(changes)
        # print(changes, shift)
        final_change = tuple(c for c in final_change if c != 0)
        if len(final_change) == 0:
            return 0

        base_term = self.subterm.get_poly_terms(changes, shift=None)
        if isinstance(base_term, PTTensorCoeffProductSum):
            prefactor = PTEnergyChangeProductSum.monomial(changes, 1)
            base_term = base_term.mul_simple(prefactor)

        if shift is not None:
            base_term = base_term.shift(shift)

        return base_term

class WavefunctionCorrection(PerturbationTheoryTerm):
    def __init__(self, parent, order, logger=None):
        super().__init__(logger=logger)

        self.parent = parent
        self.order = order

    def __repr__(self):
        return "Y[{}]".format(self.order)

    def get_poly_terms(self, changes, shift=None) -> 'PTTensorCoeffProductSum':
        if len(changes) == 0:
            return self.parent.overlap_correction(self.order).get_poly_terms([], shift=shift)
        else:
            return super().get_poly_terms(changes, shift=shift)

    def get_changes(self):
        base_changes = {}
        for expr in self.expressions:
            for subchange in expr.changes:
                base_changes[subchange] = base_changes.get(subchange, [])
                base_changes[subchange].append(expr) # just track which exprs generate the change
        base_changes[()] = None # also the constant term
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
            -PerturbationOperator.lookup(H[k-i] * W(i))
                if i > 0 else
            -PerturbationOperator.lookup(H[k])
            for i in range(0, k)
            if len(H) > k - i
        ]

        return base_terms

class WavefunctionOverlapCorrection(PerturbationTheoryTerm):
    """
    Provides a slight optimization on the base `WavefunctionCorrection`
    """

    def __init__(self, parent, order, logger=None):
        super().__init__(logger=logger)

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
            -1/2 * (W(k - i) * W(i)) for i in range(1, k)
        ]

class EnergyCorrection(PerturbationTheoryTerm):
    def __init__(self, parent, order, intermediate_normalization=False, logger=None):
        super().__init__(logger=logger)

        self.parent = parent
        self.order = order
        self.int_norm = intermediate_normalization

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

        return ([H[k]] if len(H) > k else []) + [
            (H[k - i] * W(i)) for i in range(1, k)
            if len(H) > k - i
        ] + (
            [
                -E(k - i) * O(i) for i in range(1, k)
                if (k - i) % 2 == 0
            ]
            if self.int_norm else
            []
        )

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

    def __call__(self, changes, shift=None, coeffs=None, freqs=None, check_sorting=None, simplify=True):
        if len(changes) > 0:
            raise ValueError("no")
        return super().__call__((), shift=shift, coeffs=coeffs, freqs=freqs, check_sorting=False, simplify=simplify)

class OperatorCorrection(PerturbationTheoryTerm):

    def __init__(self, parent, order, operator_type=None, logger=None):
        super().__init__(logger=logger)

        self.parent = parent
        self.order = order
        self.expansion = parent.operator_expansion_terms(order, logger=self.logger, operator_type=operator_type)

    def __repr__(self):
        return "<n|M{}|m>".format(self.order)

    def get_changes(self):
        base_changes = {}
        for expr in self.expressions:
            for subchange in expr.changes:
                base_changes[subchange] = base_changes.get(subchange, [])
                base_changes[subchange].append(expr) # just track which exprs generate the change
        base_changes[()] = None # also the constant term
        return base_changes

    def get_subexpresions(self) -> 'Iterable[PerturbationTheoryTerm]':
        H = self.parent.hamiltonian_expansion
        E = self.parent.energy_correction
        W = self.parent.wavefunction_correction
        M = self.expansion
        k = self.order
        L = lambda o:FlippedPerturbationTheoryTerm.lookup(W(o))

        exprs = (
                [M[k]]
        )+ [
            L(k-i) * M[i]
            for i in range(0, k)
        ] + [
            M[i] * W(k-i)
            for i in range(0, k)
        ] + [
            L(i) * M[k-i-j] * W(j)
            for i in range(1, k)
            for j in range(1, k-i+1)
        ]

        return exprs

class ScaledPerturbationTheoryTerm(PerturbationTheoryTerm):
    #TODO: refactor since inheritance isn't really the right paradigm here
    def __init__(self, base_term:'PerturbationTheoryTerm', scaling):
        super().__init__(logger=base_term.logger)
        self.prefactor = scaling
        self.base = base_term

    _cache = {}
    @classmethod
    def lookup(cls, base_term, scaling):
        key = (base_term, scaling)
        prod = cls._cache.get(key, None)
        if prod is None:
            prod = cls(base_term, scaling)
            cls._cache[key] = prod
        return prod

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

    _cache = {}
    def __init__(self, post_op, pre_op):
        super().__init__(logger=post_op.logger)

        self.gen1 = pre_op # we apply right-to-left
        self.gen2 = post_op

    @classmethod
    def lookup(cls, post_op, pre_op):
        key = (post_op, pre_op)
        prod = cls._cache.get(key, None)
        if prod is None:
            prod = cls(post_op, pre_op)
            cls._cache[key] = prod
        return prod

    def __repr__(self):
        return "{}*{}".format(self.gen1, self.gen2)

    _change_map = {}
    _perms_map = {}
    @classmethod
    def change_direct_sum_generator(cls, change_1, change_2):
        """
        Enumerates the direct sum + indices for terms taken from `change_1` and `change_2`

        :param change_1:
        :param change_2:
        :return:
        """

        # changes = cls.

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

        perm_inds, inv_perm = cls._perms_map.get(l1, (None, None))
        if perm_inds is None:
            perm_inds = np.array(list(itertools.permutations(range(l1))))
            inv_perm = np.argsort(perm_inds)
            cls._perms_map[l1] = (perm_inds, inv_perm)

        # perm_inds = np.asanyarray(list(itertools.permutations(range(l1))))
        # inv_perm = np.argsort(perm_inds)
        # print("::>", change_1, change_2)
        perm_sorting = None
        arr_1 = np.asanyarray(change_1)
        arr_2 = np.asanyarray(change_2)

        # print("?????", arr_2)
        change_map = cls._change_map
        for overlaps in range(1, min(l1, l2) + 1):
            # TODO: handle overlap == len(term_2) -> not sure if I've done this or not?
            (ov_parts, rem_part), (ov_inds, rem_inds) = UniquePartitions(arr_2).partitions(
                [overlaps, l2 - overlaps],
                return_indices=True
            )  # !need to track indices
            if l1 == overlaps:
                pad_ov = ov_parts
            else:
                pad_ov = np.pad(ov_parts, [[0, 0], [0, l1 - overlaps]])

            change_counts = tuple(np.unique(arr_2, return_counts=True)[1])
            for p, r, io, ir in zip(pad_ov, rem_part, ov_inds, rem_inds):
                # print("-"*50)
                # print("???", legit_inds[:10])
                # print("...", perm_ov_1[:10])
                # print(",,,", perm_inds[:10])
                # print("_"*50)


                p_counts = tuple(np.unique(p, return_counts=True)[1])
                uperm_inds = change_map.get((p_counts, change_counts), None)
                if uperm_inds is None:
                    legit_inds = change_map.get(p_counts, None)
                    if legit_inds is None:
                        legit_inds, perm_ov_1 = UniquePermutations(p).permutations(return_indices=True)
                        change_map[p_counts] = legit_inds
                    perm_inds2, sortings, _, ulegit_inds, uperm_inds = nput.intersection(legit_inds, perm_inds,
                                                                               sortings=(None, perm_sorting),
                                                                               assume_unique=True, return_indices=True
                                                                               )
                    change_map[(p_counts, change_counts)] = uperm_inds
                perm_inds2 = perm_inds[uperm_inds]

                # print("--->", p, arr_2, np.unique(p, return_counts=True)[1], np.unique(arr_2, return_counts=True)[1])
                # print(uperm_inds)
                #
                #     _, perm_sorting = sortings
                # ar1_indices, ar2_indices

                perm_ov = p[perm_inds2]  # I think I can actually take the unique perms...?
                # print("???", p, change_1, change_2, len(perm_ov), end=" ")
                # print(perm_inds)
                # po = len(perm_ov)

                # perm_ov = p[perm_inds]
                # perm_ov, uperm_inds = np.unique(perm_ov, return_index=True, axis=0)

                # print("?"*50)
                # print(perm_ov1)
                # print(perm_ov)
                # print("_"*50)
                # if po == 120:
                #     print(perm_ov)
                #     raise Exception(...)
                # print("->", len(perm_ov))

                # it might be necessary to apply the transpose here

                if len(r) > 0:
                    newts = np.concatenate([
                        arr_1[np.newaxis] + perm_ov,
                        np.broadcast_to(r[np.newaxis], (len(perm_ov), len(r)))
                    ],
                        axis=1
                    )
                else:
                    newts = arr_1[np.newaxis] + perm_ov

                newt_transposes = np.argsort(newts, axis=1)

                # build new polys
                target_axes = inv_perm[uperm_inds][:, :overlaps]  # axes in term_1 that the selected parts of
                from_inds = perm_inds[uperm_inds][:, :overlaps]  # axes in term_1 that the selected parts of
                                                                 # term_2 are multiplied by
                for base_change, transp, targ_inds, src_inds in zip(newts, newt_transposes, target_axes, from_inds):
                    yield base_change[transp], transp, targ_inds, src_inds
        # print("<::")
        # handle overlap == 0 term
        base_change = np.concatenate([change_1, change_2])
        transp = np.argsort(base_change)
        targ_inds = ()
        src_inds = ()
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

    def get_poly_terms(self, changes, allowed_paths=None, shift=None, simplify=True):
        t = tuple(changes)
        shift = tuple(shift) if shift is not None else None
        # check sorting?
        base_changes = self.changes.get(t, 0)
        if allowed_paths is not None:
            allowed_paths = set(allowed_paths)
        if isinstance(base_changes, list):
            poly_changes = 0
            for change_1, change_2, reorg, target_inds, src_inds in base_changes:
                if allowed_paths is not None and (change_1, change_2) not in allowed_paths:
                    continue

                polys_1 = self.gen1.get_poly_terms(change_1, shift=shift)
                if isinstance(polys_1, PerturbationTheoryEvaluator): polys_1 = polys_1.expr
                if nput.is_numeric(polys_1):
                    if polys_1 == 0:
                        continue
                    else:
                        raise ValueError("not sure how we got a number here...")
                if simplify:
                    polys_1 = polys_1.combine()
                # we note that src_inds defines how we permute the change_2 inds so...I guess target_inds
                # tells us how we'd permute the inds of change_1?
                if shift is not None:
                    change_shift = change_1 + shift
                else:
                    change_shift = change_1
                change_shift = [change_shift[i] for i in target_inds]

                polys_2 = self.gen2.get_poly_terms(change_2, shift=change_shift)
                if isinstance(polys_2, PerturbationTheoryEvaluator): polys_2 = polys_2.expr
                if nput.is_numeric(polys_2):
                    if polys_2 == 0:
                        continue
                    else:
                        raise ValueError("not sure how we got a number here...")
                if simplify:
                    polys_2 = polys_2.combine()

                base_polys = polys_1.mul_along(polys_2, target_inds, remainder=None, mapping=None)

                # print(">>>>", change_1, change_2)
                # print("polys_1")
                # polys_1.prune_operators([(1, 1)]).print_tree()
                # print("polys_2")
                # polys_2.prune_operators([(1, 1)]).print_tree()
                # print('final')
                # base_polys.prune_operators([(1, 1)]).print_tree()
                # print("<<<<")

                # we need to also account for both possible direction changes...
                # EDIT: NOT ACTUALLY THIS ALREADY IN THERE
                # if shift is not None:
                #     sqrt_contrib_1 = HarmonicOscillatorRaisingLoweringPolyTerms.get_direction_change_poly(change_1, shift)
                # else:
                #     sqrt_contrib_1 = None
                # if sqrt_contrib_1 is not None:
                #     # for every term in base_polys we need to convolve with this change...
                #     base_polys = base_polys.mul_simple(ProductPTPolynomial(sqrt_contrib_1))

                # # this is already in polys_2
                # sqrt_contrib_2 = HarmonicOscillatorRaisingLoweringPolyTerms.get_direction_change_poly(change_2, change_shift)
                #
                # if sqrt_contrib_2 is not None:
                #     # need to permute
                #     sqrttoooo = [[1]] * len(change_1) # number of terms in base_polys????
                #     for i,c in zip(target_inds, sqrt_contrib_2):
                #         sqrttoooo[i] = c
                #     base_polys = base_polys.mul_simple(ProductPTPolynomial(sqrttoooo))

                poly_changes += base_polys

            self.changes[t] = poly_changes
            base_changes = poly_changes

        return base_changes

class FlippedPerturbationTheoryTerm(PerturbationTheoryTerm):
    """
    Represents a term that will be multipled by on the left rather than the right
    for evaluating things like Y[1]M[0]Y[1]
    """
    def __init__(self, base_term):
        super().__init__(logger=base_term.logger)
        self.base = base_term

    def __repr__(self):
        return "/{}".format(self.base)

    # TODO: make sure changes supported are correct...
    #       NOTE: for now by doing nothing I'm assuming all operators provide symmetric changes...

    _cache = {}
    @classmethod
    def lookup(cls, base_term):
        prod = cls._cache.get(base_term, None)
        if prod is None:
            prod = cls(base_term)
            cls._cache[base_term] = prod
        return prod

    def get_changes(self):
        return self.base.changes
        # change_map = {}
        # for k,v in self.base.changes.items():
        #     print(k)
        #     k = tuple(-k)
        #     change_map[tuple(re)]



    def get_poly_terms(self, changes, shift=None) -> 'PTTensorCoeffProductSum':
        # reading right to left is equivalent to reading
        # left to right with a shift and negative changes
        if shift is None:
            shift = [0] * len(changes)
        flip_changes = [-c for c in changes]
        # if len(shift) != len(changes):
        #     raise ValueError(shift, changes)
        if len(changes) < len(shift):
            changes = list(changes) + [0]*(len(changes) - len(shift))
        shift = [s+c for s,c in zip(shift, changes)]
        return self.base.get_poly_terms(flip_changes, shift)

class PerturbationTheoryEvaluator:
    def __init__(self, expr:'PTTensorCoeffProductSum', change):
        self.expr = expr
        self.change = change
        self.num_fixed = len(change)
    def __repr__(self):
        return "{}(<{}>, {})".format(type(self).__name__, self.num_fixed, self.expr)

    @staticmethod
    def _extract_coeffs(coeffs, coeff_indices):
        coeff_tensor = coeffs[coeff_indices[0]][coeff_indices[1]] # order then identity
        if nput.is_zero(coeff_tensor): return 0

        # print(free, fixed, coeff_indices)
        idx = coeff_indices[2:]
        if len(idx) > 0:
            return coeff_tensor[idx]
        else:
            return coeff_tensor # just a number

    @classmethod
    def _eval_poly(cls, state, fixed, inds, poly, change, verbose):
        if isinstance(poly, ProductPTPolynomialSum):
            subvals = [
                cls._eval_poly(state, fixed, inds, p, change, verbose)
                for p in poly.polys
            ]
            return poly.prefactor * np.sum(subvals)

        p = poly.constant_rescale()
        if verbose:
            print("     |", round(p.prefactor, 8), [np.round(c, 8) for c in p.coeffs])

        if fixed > 0 and len(inds) > 0:
            substates = np.concatenate([state[:fixed], state[inds,]])
        elif fixed > 0:
            substates = state[:fixed]
        else:
            substates = state[inds,]

        # print(substates, inds, fixed)

        poly_factor = np.prod([
            np.dot(c, np.power(s, np.arange(len(c))))
            for s,c in zip(substates, poly.coeffs)
        ])
        sqrt_factor = np.sqrt(
            np.prod(
                [
                    n + i
                    for n, delta in zip(substates, change)
                    for i in range(delta + 1 if delta < 0 else 1, 1 if delta < 0 else delta + 1)
                ]
            )
        )
        return poly.prefactor * poly_factor * sqrt_factor

    def _eval_perm(self, subset, state, coeffs, freqs, cind_sets, zero_cutoff, eval_cache, verbose):
        #TODO: use combinatorics to skip some duplicate evaluations over permutations

        # print("?", cind_sets)
        full_set = tuple(range(self.num_fixed)) + subset

        free = subset
        fixed = self.num_fixed
        cinds_remapped = [
            tuple(
                # we need to avoid doing the exact same term twice, _but_ order matters
                # so we can't sort the tuples
                PTTensorCoeffProductSum._symmetrize(
                    ci[:2] + tuple(free[j - fixed] if j >= fixed else j for j in ci[2:])
                )
                for ci in cinds
            )
            for cinds in cind_sets
        ]
        prefactors = np.array([
            np.prod([
                self._extract_coeffs(coeffs, ci) for ci in cinds])
                    if eval_cache.get(cinds, 0) < np.math.factorial(len(cinds)) else # should cache this
                0
            for cinds in cinds_remapped
        ])
        for cinds in cinds_remapped:
            eval_cache[cinds] = eval_cache.get(cinds, 0) + 1
        good_pref = np.where(np.abs(prefactors) > zero_cutoff)  # don't bother to include useless terms
        if len(good_pref) == 0:
            return 0
        good_pref = good_pref[0]

        contrib = 0
        for g in good_pref:
            if verbose:
                print("","-"*10, PTTensorCoeffProductSum.format_key(cind_sets[g]), prefactors[g], "-"*10)
            subexpr = self.expr.terms[cind_sets[g]]
            if isinstance(subexpr, PTEnergyChangeProductSum):
                subcontrib = 0
                for echanges, polys in subexpr.terms.items():
                    energy_factor = np.prod([
                        np.sum([c * freqs[s] for c, s in zip(eecc, full_set)])
                        for eecc in echanges
                    ])
                    # TODO: check size of energy factor
                    if hasattr(polys, 'polys'):
                        if verbose:
                            print("  ::", polys.prefactor)
                        # for p in polys.polys:
                        #     print("->", p.prefactor, p.coeffs)
                    poly_factor = self._eval_poly(state, self.num_fixed, subset, polys, self.change, verbose)
                    if verbose:
                        print("  ..", PTEnergyChangeProductSum.format_key(echanges), poly_factor, energy_factor, poly_factor*prefactors[g]/energy_factor * 219475)
                    subcontrib += poly_factor / energy_factor
                subcontrib *= subexpr.prefactor
            elif isinstance(subexpr, (ProductPTPolynomial, ProductPTPolynomialSum)):
                subcontrib = self._eval_poly(state, self.num_fixed, subset, subexpr, self.change, verbose)
            else:
                raise ValueError("how the hell did we end up with {}".format(subexpr))

            if verbose:
                print(" ->", prefactors[g]*subcontrib*219475)
            contrib += prefactors[g] * subcontrib
            # base_prefactor = subexpr.prefactor
            # if isinstance(subexpr, ProductPTPolynomialSum):
            #     subprefactors = [p.prefactor for p in subexpr.polys]
            #     poly_coeffs = [p.coeffs for p in subexpr.polys]
            # else:
            #     subprefactors = []
            #     poly_coeffs = subexpr.coeffs
            # print(poly_coeffs)

        if verbose:
            print("::>", contrib * 219475)

        return contrib

    def evaluate(self, state, coeffs, freqs, zero_cutoff=None, verbose=True):
        # we do state-by-state evaluation for now although we will at some
        # point need to do this in batches

        state = np.asanyarray(state)
        ndim = len(state)
        # max_order = None
        free_ind_groups = {}
        for coeff_indices,poly_terms in self.expr.terms.items():
            num_inds = len(self.expr.get_inds(coeff_indices))
            free_inds = num_inds - self.num_fixed
            if free_inds not in free_ind_groups:
                free_ind_groups[free_inds] = []
            free_ind_groups[free_inds].append(coeff_indices)
        #     max_order = max(max_order, poly_terms.order) if max_order is not None else poly_terms.order
        #
        # state_polys = np.power(state, np.arange(max_order+1))

        if zero_cutoff is None:
            zero_cutoff = 1e-18

        # TODO: numba-ify this part
        contrib = 0
        eval_cache = {} # so we don't hit some coeff sets too many times
        for free_inds,cind_sets in free_ind_groups.items():
            # now we iterate over every subset of inds (outside of the fixed ones)
            # excluding replacement (since that corresponds to a different coeff/poly)
            if free_inds == 0:
                if verbose:
                    print("="*7, (), "="*7)
                contrib += self._eval_perm((), state, coeffs, freqs, cind_sets, zero_cutoff, eval_cache, verbose)
            else:
                for subset in itertools.combinations(range(self.num_fixed, ndim), r=free_inds):
                    for subperm in itertools.permutations(subset):
                        if verbose:
                            print("="*7, subperm, "="*7)
                        contrib += self._eval_perm(subperm, state, coeffs, freqs, cind_sets, zero_cutoff, eval_cache, verbose)




                # raise Exception(cind_sets)

        return self.expr.prefactor * contrib

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
