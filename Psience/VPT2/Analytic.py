"""
Provides a symbolic approach to vibrational perturbation theory based on a Harmonic description
"""

import abc, itertools, collections, enum
import numpy as np, scipy.signal, time

from McUtils.Zachary import DensePolynomial, TensorCoefficientPoly
import McUtils.Numputils as nput, McUtils.Combinatorics as comb
from McUtils.Combinatorics import SymmetricGroupGenerator, IntegerPartitioner, UniquePartitions, UniquePermutations
from McUtils.Scaffolding import Logger, Checkpointer
from ..BasisReps import (
    BasisStateSpace, HarmonicOscillatorProductBasis,
    HarmonicOscillatorMatrixGenerator, HarmonicOscillatorRaisingLoweringPolyTerms
)

__all__ = [
    'PerturbationTheoryEvaluator',
    'AnalyticPerturbationTheorySolver',
    # 'AnalyticPerturbationTheoryDriver',
    # 'AnalyticPTCorrectionGenerator',
    # 'RaisingLoweringClasses'
]


_DEBUG_PRINT = True # statements will be moved into a proper logger later

# class DefaultObject:
#     def __repr__(self): return "default"
class DefaultValues(enum.Enum):
    DEFAULT = 'default'
default = DefaultValues.DEFAULT

class AnalyticPerturbationTheorySolver:
    """
    A re-attempt at using the recursive expressions
    to provide simpler code for getting APT expressions
    """
    def __init__(self, hamiltonian_expansion, logger=None, checkpoint=None, intermediate_normalization=False):
        self.hamiltonian_expansion = hamiltonian_expansion
        self.logger = Logger.lookup(logger)
        self.checkpoint = Checkpointer.build_canonical(checkpoint)
        self.intermediate_normalization = intermediate_normalization

    @classmethod
    def from_order(cls, order, internals=True, logger=None, checkpoint=None, intermediate_normalization=False):
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

        return cls(hamiltonian_terms, logger=logger, checkpoint=checkpoint, intermediate_normalization=intermediate_normalization)

    _op_maps = {}
    def get_correction(self, key, cls, order, allowed_terms=None, allowed_changes=None, **kw):
        k = (key, order)
        corr = self._op_maps.get(k, None)
        if corr is None:
            corr = cls(self, order,
                       logger=self.logger, checkpoint=self.checkpoint,
                       intermediate_normalization=self.intermediate_normalization,
                       allowed_terms=allowed_terms,
                       allowed_changes=allowed_changes,
                       **kw
                       )
            self._op_maps[k] = corr
        return corr

    def shifted_hamiltonian_correction(self, order, **kw):
        return self.get_correction('dH', ShiftedHamiltonianCorrection, order, **kw)

    def energy_correction(self, order, **kw):
        return self.get_correction('energy', EnergyCorrection, order, **kw)

    def wavefunction_correction(self, order, **kw):
        return self.get_correction('proj_wfn', WavefunctionCorrection, order, **kw)

    def overlap_correction(self, order, **kw):
        return self.get_correction('ov', WavefunctionOverlapCorrection, order, **kw)

    def full_wavefunction_correction(self, order, **kw):
        return self.get_correction('wfn', FullWavefunctionCorrection, order, **kw)

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

class PolynomialInterface(metaclass=abc.ABCMeta):
    """
    Provides a basic interface to allow for the uniform manipulation
    of objects that dispatch down to some form of scalar multiplied by
    a sum of polynomials
    """

    def format_expr(self) -> str:
        ...

    @property
    @abc.abstractmethod
    def ndim(self) -> int:
        ...
    # @property
    # @abc.abstractmethod
    # def order(self) -> 'tuple[int]':
    #     ...
    @abc.abstractmethod
    def audit(self, target, ignore_constants=True):
        ...
    @abc.abstractmethod
    def ensure_dimension(self, ndim) -> 'Self':
        ...
    def align_dimensions(self, other) -> 'tuple[Self, Self]':
        max_dim = max(self.ndim, other.ndim)
        return self.ensure_dimension(max_dim), other.ensure_dimension(max_dim)
    # @property
    # @abc.abstractmethod
    # def order(self):
    #     ...
    @abc.abstractmethod
    def shift(self, shift)->'Self':
        ...
    @abc.abstractmethod
    def scale(self, scaling)->'Self':
        ...
    @abc.abstractmethod
    def permute(self, perm, check_perm=True)->'Self':
        ...

    @abc.abstractmethod
    def mul_simple(self, other:'PolynomialInterface')->'PolynomialInterface':
        ...
    @abc.abstractmethod
    def rmul_simple(self, other: 'PolynomialInterface') -> 'PolynomialInterface':
        ...
    @abc.abstractmethod
    def mul_along(self, other:'PolynomialInterface', inds, remainder=None, mapping=None)->'PolynomialInterface':
        ...
    @abc.abstractmethod
    def rmul_along(self, other:'PolynomialInterface', inds, remainder=None, mapping=None)->'PolynomialInterface':
        ...

    @abc.abstractmethod
    def combine(self, **kwargs) -> 'Self':
        ...
    # @abc.abstractmethod
    # def condense(self, inds=None, return_inds=False, check_inds=True):
    #     ...
    @abc.abstractmethod
    def mutate(self, *args, **kwargs) -> 'Self':
        ...

class ProductPTPolynomial(PolynomialInterface):
    """
    TODO: include prefactor term so we can divide out energy changes
    """
    def __init__(self, coeffs, prefactor=1, idx=None, sqrt_scaling=None):
        if (
                nput.is_numeric(coeffs)
                or nput.is_numeric(coeffs[0])
        ):
            raise ValueError("coeffs must be a vector of vectors (not {})".format(coeffs))
        self.coeffs = [np.asanyarray(c) for c in coeffs] # coeffs along each dim independently
        if any(c.dtype == np.dtype(object) for c in self.coeffs):
            raise ValueError(self.coeffs)
        self.prefactor = prefactor
        self._order = None
        self._monics = None
        self._idx = idx

    def mutate(self, coeffs:'list[np.ndarray]'=default, prefactor:'Any'=default, idx=None):
        return type(self)(
            self.coeffs if coeffs is default else coeffs,
            prefactor=self.prefactor if prefactor is default else prefactor,
            idx=idx
        )

    @property
    def ndim(self):
        return len(self.coeffs)
    @property
    def order(self):
        if self._order is None:
            self._order = tuple(len(c)-1 for c in self.coeffs)
        return self._order
    def audit(self, target, ignore_constants=True):
        if len(self.coeffs) != target:
            if ignore_constants and not (
                    len(self.coeffs) == 1 and
                    len(self.coeffs[0]) == 1 and
                    target == 0
            ):  # ignore constants
                raise ValueError("{} has wrong dimension (expected {})".format(
                    self,
                    target
                ))
    def ensure_dimension(self, ndim):
        return self.pad(ndim - self.ndim) if self.ndim < ndim else self
    def pad(self, left_right_pads):
        if nput.is_numeric(left_right_pads):
            left_right_pads = [0, left_right_pads]
        l,r = left_right_pads
        new_coeffs = [
            np.array([1]) for _ in range(l)
        ] + self.coeffs + [
            np.array([1]) for _ in range(r)
        ]

        return self.mutate(new_coeffs)

    def __repr__(self):
        return "{}(<{}>)".format("Poly", ",".join(str(c) for c in self.order))

    @classmethod
    def format_simple_poly(cls, coeffs, idx):
        var = "n[{}]".format(idx)
        return "+".join(k for k in [
            "{}{}".format(
                (
                    "" if c == 1 else
                    "-" if c == -1 else
                    round(c, 4)
                ) if i != 0 else round(c, 4),
                "" if i == 0 else
                var if i == 1 else
                var + "^{}".format(i)
            ) for i, c in enumerate(coeffs)
            if c != 0
        ] if len(k) > 0)
    def format_expr(self):
        subpolys = [
            self.format_simple_poly(c, i)
            for i,c in enumerate(self.coeffs)
        ]
        subpolys = ["({})".format(s) for s in subpolys if len(s) > 0]
        c = self.prefactor
        poly_str = "".join(s for s in subpolys if len(s) > 0)
        if poly_str == "":
            return "" if c == 1 else str(c)
        else:
            prefac_str = "" if c == 1 else "-" if c == -1 else "{}*".format(round(c, 4))
            return prefac_str + poly_str

    def permute(self, new_inds, check_perm=True):
        # must be a proper permutation
        if check_perm:
            if np.any(
                    np.sort(np.concatenate([new_inds, np.arange(len(new_inds), len(self.coeffs))]))
                        != np.arange(len(self.coeffs))
            ):
                raise ValueError("bad permutation {}")
        return self.mutate(
            [self.coeffs[i] for i in new_inds] + self.coeffs[len(new_inds):]
        )

    def constant_rescale(self):
        """
        rescales so constant term is 1

        :return:
        """
        new_prefactor = self.prefactor
        new_coeffs = []
        for c in self.coeffs:
            if abs(c[0]) > 1e-8: # non-zero
                new_prefactor = self.prefactor * c[0]
                new_coeffs.append(c/c[0])
            else:
                new_coeffs.append(c)
        return self.mutate(new_coeffs, prefactor=new_prefactor)

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
        # scaling = abs(coeffs[-1])
        # if scaling <= zero_thresh:
        #     return np.zeros_like(coeffs), 0
        # elif scaling == 1:
        #     return coeffs, 1
        # else:
        #     return coeffs / scaling, scaling

        acoeffs = np.abs(coeffs) # np.round(np.abs(coeffs), 8)
        acmask = np.where(acoeffs > zero_thresh)
        if len(acmask) == 0 or len(acmask[0]) == 0:
            return np.zeros_like(coeffs), 0

        acmask = acmask[0]
        min_v = np.min(acoeffs[acmask])
        new_coeffs = np.zeros_like(coeffs)
        new_coeffs[acmask] = coeffs[acmask] / min_v
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
            new.monic_coeffs = [(c,1) for c in new_coeffs]
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

    def condense(self, inds=None, return_inds=False, check_inds=True):
        if inds is None: inds = np.arange(len(self.coeffs))
        if check_inds:
            condensed = np.array([i for i in inds if len(self.coeffs[i]) == 1])
        else:
            condensed = inds
        if len(condensed) > 0:
            new = self.mutate(
                [c for c in self.coeffs if c not in condensed]
            )
        else:
            new = self
        if return_inds:
            return condensed, new
        else:
            return new


    def shift(self, shift):
        if len(shift) < len(self.coeffs):
            shift = list(shift) + [0] * (len(self.coeffs) - len(shift))
        if self._idx is not None:
            base_shift = self._idx[-1]
            new_shift = tuple(s + k for s,k in zip(base_shift, shift)) + base_shift[len(shift):]
            new_idx = self._idx[:-1] + (new_shift,)
        else:
            new_idx = None
        return self.mutate(
            [DensePolynomial.compute_shifted_coeffs(c, s) for c,s in zip(self.coeffs, shift)],
            idx=new_idx
        )

    def __mul__(self, other):
        if nput.is_numeric(other):
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
        return self.mutate(prefactor=self.prefactor*scalar)

    def mul_simple(self, other:'PolynomialInterface'):
        if nput.is_numeric(other):
            if other == 1:
                return self
            elif other == 0:
                return 0
            else:
                raise ValueError(other)
        elif isinstance(other, ProductPTPolynomial):
            return self._poly_mul(self, other)
        else:
            return other.rmul_simple(self)
    def rmul_simple(self, other: 'PolynomialInterface') -> 'PolynomialInterface':
        return other.mul_simple(self)

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
            new = cls(
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
            chunks = np.array([])
        return chunks

    @classmethod
    def get_index_mapping(cls, dim_1, dim_2, inds, return_remainder=True):
        """
        Computes the corresponding and indices remainders for multdimensional polynomial multiplications

        :param dim_1:
        :param dim_2:
        :param inds:
        :return:
        """
        remainder = None
        if nput.is_numeric(inds):
            # simplest form, multiply the first indices
            k = inds
            inds = [
                list(range(k)),
                list(range(k))
            ]
            if return_remainder:
                remainder = [
                    np.arange(k, dim_1),
                    np.arange(k, dim_2)
                ]
        elif len(inds) == 0:
            # k = inds
            inds = [
                [],
                []
            ]
            if return_remainder:
                remainder = [
                    np.arange(dim_1),
                    np.arange(dim_2)
                ]
        elif nput.is_numeric(inds[0]):
            # only the 2nd poly inds are specified
            k = len(inds)
            inds = [
                list(range(k)),
                inds
            ]
            if return_remainder:
                remainder = [
                    np.arange(k, dim_1),
                    cls.fast_ind_remainder(dim_2, inds[1]),
                ]
        else:
            if return_remainder:
                remainder = [
                    cls.fast_ind_remainder(dim_1, inds[0]),
                    cls.fast_ind_remainder(dim_2, inds[1])
                ]

        # subinds = [
        #     n for n,(i1, i2) in enumerate(zip(*inds))
        #     if i1 < dim_1 and i2 < dim_2
        # ]

        new_rems = [
            [],
            [],
        ]
        new_inds = [
            [],
            []
        ]
        for i1, i2 in zip(*inds):
            if i1 < dim_1 and i2 < dim_2:
                new_inds[0].append(i1)
                new_inds[1].append(i2)
            else:
                new_rems[0].append(i1)
                new_rems[1].append(i2)


        inds = new_inds
        remainder = [
            list(remainder[0]) + new_rems[0],
            list(remainder[1]) + new_rems[1]
        ]

        res = inds
        if return_remainder:
            res = (inds, remainder)
        return res

    def mul_along(self, other:'PolynomialInterface', inds, remainder=None, mapping=None):
        if isinstance(other, ProductPTPolynomial):
            # we take the polynomial product of indices that align
            # and direct product the remainder
            if len(inds) == 0: # pure product
                # we get to just concatenate the coeffs
                new = self.mutate(
                    #TODO: think about whether mutate is the right concept here...
                    self.coeffs+other.coeffs,
                    prefactor=self.prefactor*other.prefactor
                )
            else:
                (self_inds, other_inds), (self_remainder, other_remainder) = self.get_index_mapping(
                    len(self.coeffs), len(other.coeffs), inds, return_remainder=True
                )

                new_coeffs = [
                                 scipy.signal.convolve(
                                     other.coeffs[ox] if len(other.coeffs) > ox else [1],
                                     self.coeffs[sx] if len(self.coeffs) > sx else [1]
                                 )
                                 for ox, sx in zip(other_inds, self_inds)
                                 # if len(other.coeffs) > ox and len(self.coeffs) > sx
                             ] + [
                                 self.coeffs[sx] if len(self.coeffs) > sx else [1]
                                 for sx in self_remainder
                                 # if len(self.coeffs) > sx
                             ] + [
                                 other.coeffs[ox] if len(other.coeffs) > ox else [1]
                                 for ox in other_remainder
                                 # if len(other.coeffs) > ox
                             ]

                new = self.mutate(
                    new_coeffs,
                    prefactor=self.prefactor*other.prefactor
                )
        else:
            new = other.rmul_along(self, inds, remainder=remainder, mapping=mapping)

        return new

    def rmul_along(self, other, inds, remainder=None, mapping=None):
        """
        We multiply every subpoly along the given indices, transposing the appropriate tensor indices

        :param other:
        :param inds:
        :param remainder:
        :return:
        """
        return other.mul_along(self, inds, remainder=remainder, mapping=mapping)

    def __add__(self, other):
        if isinstance(other, (int, float, np.integer, np.floating)):
            if other == 0:
                return self
            else:
                raise NotImplementedError(...)
                return self + ProductPTPolynomial([[other]])
        elif isinstance(other, ProductPTPolynomial):
            return ProductPTPolynomialSum([self, other])#, prefactor=self.prefactor)
        elif isinstance(other, ProductPTPolynomialSum):
            return other + self
        else:
            raise NotImplementedError("not sure how to add {} and {}".format(self, other))
    def __radd__(self, other):
        return self + other

class ProductPTPolynomialSum(PolynomialInterface):

    def __init__(self, polynomials, prefactor=1, ndim=None, order=None, reduced=False):
        self.polys = polynomials
        self.prefactor = prefactor
        self.reduced = reduced
        self._order = order
        self._ndim = ndim

    def mutate(self, polynomials:'list[PolynomialInterface]'=default,
               prefactor:'Any'=default,
               ndim:'Any'=default,
               order:'Any'=default,
               reduced:'Any'=default
               ):
        return type(self)(
            self.polys if polynomials is default else polynomials,
            prefactor=self.prefactor if prefactor is default else prefactor,
            ndim=self._ndim if ndim is default else ndim,
            order=self._order if order is default else order,
            reduced=self.reduced if reduced is default else reduced
        )

    @property
    def ndim(self):
        if self._ndim is None:
            self.audit()
            self._ndim = self.polys[0].ndim
        return self._ndim
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
    def audit(self, target=None, ignore_constants=True):
        if target is None: target = self.polys[0].ndim
        for p in self.polys: p.audit(target, ignore_constants=ignore_constants)
    def ensure_dimension(self, ndim):
        new_ndim = self._ndim
        if new_ndim is not None:
            if new_ndim >= ndim:
                return self
            else:
                new_ndim = ndim
        new_polys = [
            k.pad(ndim - k.ndim) if k.ndim < ndim else k
            for k in self.polys
        ]
        return self.mutate(
            new_polys,
            ndim=new_ndim
        )

    def __repr__(self):
        form_counts = {}
        for p in self.polys:
            key = "<{}>".format(",".join(str(c) for c in p.order))
            form_counts[key] = form_counts.get(key, 0) + 1
        return "PSum({})".format(
            # type(self).__name__,
            "+".join("{}{}".format(k,f) for k,f in form_counts.items())
        )

    def format_expr(self):
        base_sum = "+".join(p.format_expr() for p in self.polys)
        if self.prefactor != 1:
            if "+" in base_sum: base_sum = "({})".format(base_sum)
            if self.prefactor == -1:
                base_sum = "-"+base_sum
            else:
                base_sum = "{}*{}".format(self.prefactor, base_sum)
        return base_sum

    def permute(self, new_inds, check_perm=True):
        return self.mutate(
            [p.permute(new_inds, check_perm=check_perm) for p in self.polys],
            prefactor=self.prefactor,
            reduced=self.reduced
        )

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
        return self.mutate(
            polys,
            prefactor=self.prefactor,
            reduced=True
        )

    def condense(self, inds=None, return_inds=False, check_inds=True):
        if inds is None: inds = np.arange(self.ndim)
        if check_inds:
            o = self.order
            condensed = np.array([i for i in inds if o[i] == 0])
        else:
            condensed = inds
        if len(condensed) > 0:
            new = self.mutate(
                [p.condense for p in self.polys]
            )
        else:
            new = self
        if return_inds:
            return condensed, new
        else:
            return new

    def shift(self, shift):
        return self.mutate(
            [p.shift(shift) for p in self.polys]
        )

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
        return self.mutate(prefactor=self.prefactor*scalar)

    def mul_simple(self, other:'ProductPTPolynomial'):
        if isinstance(other, ProductPTPolynomial):
            return self.mutate(
                [
                    poly.mul_simple(other)
                    for poly in self.polys
                ]
            )
        elif isinstance(other, ProductPTPolynomialSum):
            return self.mutate(
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
            return other.rmul_simple(self)
    def rmul_simple(self, other: 'PolynomialInterface') -> 'PolynomialInterface':
        if isinstance(other, ProductPTPolynomial):
            return self.mutate(
                [
                    other.mul_simple(poly)
                    for poly in self.polys
                ]
            )
        elif nput.is_numeric(other):
            if other == 1:
                return self
            elif other == 0:
                return 0
            else:
                raise ValueError(other)
        else:
            return other.mul_simple(self)

    def mul_along(self, other:'PolynomialInterface', inds, remainder=None, mapping=None):
        if isinstance(other, ProductPTPolynomial):
            # can just distribute except at some level I have to assume
            # that every stored poly has the same dimension...?
            new = self.mutate([
                poly.mul_along(other, inds, remainder=remainder, mapping=mapping)  # how is this supposed to work?
                for poly in self.polys
            ])
        elif isinstance(other, ProductPTPolynomialSum):
            new = self.mutate(
                [
                    p1.mul_along(p2, inds, remainder=remainder, mapping=mapping)
                    for p1 in self.polys
                    for p2 in other.polys
                ],
                prefactor=self.prefactor*other.prefactor
            )
        else:
            new = other.rmul_along(self, inds, remainder=remainder, mapping=mapping)
        return new
    def rmul_along(self, other:'PolynomialInterface', inds, remainder=None, mapping=None):
        """
        We multiply every subpoly along the given indices, transposing the appropriate tensor indices

        :param other:
        :param inds:
        :param remainder:
        :return:
        """
        if isinstance(other, ProductPTPolynomial):
            return self.mutate([
                other.mul_along(poly, inds, remainder=remainder, mapping=mapping)  # how is this supposed to work?
                for poly in self.polys
            ])
        else:
            return other.mul_along(self, inds, remainder=remainder, mapping=mapping)

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
                return self.mutate(self.polys + ops)
            elif isinstance(other, ProductPTPolynomial):
                if self.prefactor != 1:
                    other = other.scale(1/self.prefactor)
                    # raise NotImplementedError("need to handle prefactors")
                return self.mutate(self.polys + [other])
            else:
                raise TypeError("not sure what to do with {} and {}".format(type(self), type(other)))
    def __radd__(self, other):
        return self + other

class PTEnergyChangeProductSum(TensorCoefficientPoly, PolynomialInterface):
    """
    A representation of a sum of 1/energy * poly sums
    which is here so we can transpose energy change indices intelligently
    """
    def __init__(self, terms: dict, prefactor=1, canonicalize=True, reduced=False):
        if canonicalize:
            self._ndim = self.get_key_ndim(terms)
        else:
            self._ndim = len(next(iter(terms.keys()))[0]) # we assume all of these have the same length
        super().__init__(terms, prefactor=prefactor, canonicalize=canonicalize)
        self.reduced = reduced
        if canonicalize:
            for k in terms:
                if any(all(s == 0 for s in sk) for sk in k): raise ValueError("huh")
        # self.audit()

    def mutate(self, terms:dict=default, prefactor:'Any'=default,  canonicalize:'Any'=default, reduced:'Any'=default):
        return type(self)(
            self.terms if terms is default else terms,
            canonicalize=(terms is not default) if canonicalize is default else canonicalize,
            prefactor=self.prefactor if prefactor is default else prefactor,
            reduced=self.reduced if reduced is default else reduced
        )

    @classmethod
    def canonical_key(cls, monomial_tuple):
        l = max(len(l) for l in monomial_tuple) if len(monomial_tuple) > 0 else 0
        return super().canonical_key(
            tuple(
                t + (0,) * (l - len(t))
                for t in monomial_tuple
            )
        )

    @classmethod
    def format_key(self, key):
        return "E-[{}]".format(
            "x".join(
                ("({})" if len(k) > 1 else "{}").format(
                    "+".join(
                        "{s}w{i}".format(s=("" if s == 1 else "-" if s == -1 else s), i=i)
                        for i,s in enumerate(k)
                    )
                )
                for k in key
                )
        )
    def __repr__(self):
        sums = []
        for k_prod,v in self.terms.items():
            sums.append("{}{}".format(v, self.format_key(k_prod)))
        return "ESum({})".format("+".join(sums) if len(sums) < 3 else "\n      +".join(sums))

    @classmethod
    def format_energy_prod_key(cls, key):
        substr = [
            "+".join(
                "{s}w[{i}]".format(s=("" if s == 1 else "-" if s == -1 else s), i=i)
                for i, s in enumerate(k)
            )
            for k in key
        ]
        bb = "*".join(
            ("({})" if len(s) > 3 else "{}").format(s)
            for s in substr
        )
        return ("({})" if "*" in bb else "{}").format(bb)

    def format_expr(self):
        """
        Formats in a Mathematica-ingestible format
        :return:
        """
        keys = []
        substrings = []
        for k,poly in self.terms.items():
            keys.append(self.format_energy_prod_key(k))
            substrings.append(poly.format_expr())

        join = " + " if sum(len(x) for x in substrings) < 80 else "\n  + "
        key_prods = [
            (
                "({})/{}"
                    if "+" in s else
                "{}/{}"
            ).format(s, k)
            for k,s in zip(keys, substrings)
        ]
        return join.join(key_prods)

    @staticmethod
    def shift_key(key, shift):
        return tuple(
            tuple(k + (shift[i] if i < len(shift) else 0) for i,k in enumerate(ksum))
            for ksum in key
        )
    def shift(self, shift, shift_energies=False):
        new_terms = {}
        for k,p in self.terms.items():
            if shift_energies:
                k2 = self.shift_key(k, shift)
                if all(any(s != 0 for s in sk) for sk in k2):
                    new_terms[k2] = p.shift(shift)
            else:
                new_terms[k] = p.shift(shift)
        return self.mutate(new_terms)

    def scale(self, scaling):
        if nput.is_numeric(scaling) and scaling == 1: return self
        return self.mutate(
            prefactor=self.scaling * scaling
        )

    @classmethod
    def get_key_ndim(cls, terms:dict):
        if len(terms) == 0: return 0
        return max(
            max(len(k) for k in key_tup)
            for key_tup in terms.keys()
        )
    @property
    def ndim(self):
        if self._ndim is None:
            self._ndim = self.get_key_ndim(self.terms)
        return self._ndim

    def audit(self, target=None, ignore_constants=True):
        # if target is None:
        #     target = self.ndim
        for change,poly in self.terms.items():
            if target is not None and len(change[0]) != target:
                raise ValueError(
                    "({}) energy change {} doesn't have target dimension ({})".format(
                        type(self)({change: poly}),
                        change, target
                    )
                )
            if not nput.is_numeric(poly):
                try:
                    poly.audit(len(change[0]), ignore_constants=ignore_constants)
                except ValueError:
                    raise ValueError("bad polynomial for changes {}: {}".format(
                        change,
                        poly
                    ))
    def ensure_dimension(self, ndim):
        nd = self.ndim
        if nd >= ndim: return self

        new_terms = {}
        for key,polys in self.terms.items():
            new_key = tuple(k + (0,) * (ndim - len(k)) for k in key)
            new_terms[new_key] = polys.ensure_dimension(ndim) if not nput.is_numeric(polys) else polys

        return self.mutate(new_terms)

    def sort(self):
        return self.mutate(
            {
                k: self.terms[k] for k in sorted(self.terms.keys())
            }
        )

    def _permute_changes(self, changes, new_inds):
        # rem = ProductPTPolynomial.fast_ind_remainder(len(changes), new_inds)
        new = tuple(
            changes[i] for i in new_inds
        ) + changes[len(new_inds):]
        return new
    def permute(self, new_inds, check_perm=True):
        new_terms = {}
        for energy_changes,polys in self.terms.items():
            new_ech = tuple(
                self._permute_changes(ec, new_inds)
                for ec in energy_changes
            )
            new_terms[new_ech] = polys.permute(new_inds, check_perm=check_perm)
        return self.mutate(new_terms)

    @staticmethod
    def _check_neg(t1, t2):
        if len(t1) != len(t2):
            return False
        else:
            return all(
                all(ttt1 == -ttt2 for ttt1, ttt2 in zip(tt1, tt2))
                for tt1, tt2 in zip(t1, t2)
            )

    @staticmethod
    def find_term_scaling(key):
        new_keys = []
        new_scalings = []
        for k in key:
            s = np.gcd.reduce(k)
            new_scalings.append(s)
            new_keys.append([kk/s for kk in k])

        return tuple(new_keys), np.prod(new_scalings)

    # def factor_energies(self):
    #     ...


    def combine_energies(self):
        # terms = self.factor_energies()
        terms = self.terms

        base_keys = list(terms.keys())
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
            new_terms[k1] = terms[k1] + terms[k2].scale(-1)
        for k in unmerged:
            new_terms[k] = terms[k]

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
        return self.mutate(new_terms, reduced=True)

    @staticmethod
    def _permute_idx(idx, inds, rem, ndim):
        """
        Builds a new change index by aligning everything that matches in inds and
        incrementing everything that doesn't in remainder
        """
        # if len(inds) == 0: return idx
        # if not nput.is_numeric(inds[0]):
        #     inds = inds[1] # we remap the right inds
        if not isinstance(inds, tuple): inds = tuple(inds)
        # if remainder is None:
        #     rem = ProductPTPolynomial.fast_ind_remainder(max_remapping, inds)
        # else:
        #     rem = remainder[1]
        new = tuple(
            idx[inds[i]]
                if i < len(inds) else
            idx[rem[i - len(inds)] if i - len(inds) < len(rem) else i]
            for i in range(ndim)
        )
        if all(n == 0 for n in new): raise ValueError(idx, inds, new)
        return new
    @classmethod
    def _build_new_echange_key(cls, k1, k2, inds, remainder):
        l1 = len(k1[0])
        l2 = len(k2[0]) #if len(k2) > 0 else 0
        # r1 = l1 - len(inds[0]) - len(remainder[0])
        # r2 = l2 - len(inds[1]) - len(remainder[1])
        # ndim = l1 + l2 - len(inds[0])
        i1 = len(inds[0])
        r1 = len(remainder[0])
        i2 = len(inds[1])
        r2 = len(remainder[1])
        new_key = tuple(
            tuple(k[i] for i in inds[0])
            + tuple(k[j] for j in remainder[0])
            + k[i1 + r1:]
            + (0,) * (l2 - i2)
            for k in k1
        ) + tuple(
            tuple(k[i] for i in inds[1])
            + (0,) * (l1 - i1)
            + tuple(k[j] for j in remainder[1])
            + k[i2 + r2:]
            for k in k2
        )
        return new_key
    @classmethod
    def _pad_echange_key_right(cls, k, ndim, inds, remainder):
        return cls._build_new_echange_key(k, ((0,)*ndim,), inds, remainder)[:-1]
    @classmethod
    def _pad_echange_key_left(cls, k, ndim, inds, remainder):
        return cls._build_new_echange_key(((0,)*ndim,), k, inds, remainder)[1:]
    def mul_along(self, other:'PolynomialInterface', inds, remainder=None, mapping=None):
        """
        We multiply every subpoly along the given indices, transposing the appropriate tensor indices

        :param other:
        :param inds:
        :param remainder:
        :return:
        """
        if isinstance(other, PTEnergyChangeProductSum):
            new_poly = self.direct_product(other,
                                       key_func=lambda k1,k2:self._build_new_echange_key(k1, k2, inds, remainder),
                                       mul=lambda a,b,inds=inds:a.mul_along(b, inds, remainder=remainder, mapping=mapping)
                                       )
        elif isinstance(other, (ProductPTPolynomial, ProductPTPolynomialSum)):
            new_poly = self.mutate(
                {
                    self._pad_echange_key_right(k, other.ndim, inds, remainder): (
                        other.scale(v)
                            if nput.is_numeric(v) else
                        v.mul_along(other, inds, remainder=remainder, mapping=mapping)
                    )
                    for k, v in self.terms.items()
                }
            )
        else:
            raise NotImplementedError(type(other))
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
        elif isinstance(other, (ProductPTPolynomial, ProductPTPolynomialSum)):
            return self.mutate(
                {
                    self._pad_echange_key_left(k, other.ndim, inds, remainder): (
                        other.scale(v)
                            if nput.is_numeric(v) else
                        other.mul_along(v, inds, remainder=remainder, mapping=mapping)
                    )
                    for k, v in self.terms.items()
                }
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
        if isinstance(other, (PTEnergyChangeProductSum)):
            return self.direct_product(other,
                                       # key_func=lambda k1,k2:k1+tuple(k2[i] for i in inds), #inds says how k2 permutes?
                                       mul=lambda a,b:a.mul_simple(b)
                                       )
        elif isinstance(other, (ProductPTPolynomial, ProductPTPolynomialSum)):
            return self.mutate(
                {
                    k: (
                        other.scale(v)
                            if isinstance(v, (int, float, np.integer, np.floating)) else
                        v.mul_simple(other)
                    )
                    for k, v in self.terms.items()
                }
            )
        elif nput.is_numeric(other):
            if other == 1:
                return self
            elif other == 0:
                return 0
            else:
                raise ValueError(other)
        else:
            return other.rmul_simple(self)
    def rmul_simple(self, other: 'PolynomialInterface') -> 'PolynomialInterface':
        if isinstance(other, (ProductPTPolynomial, ProductPTPolynomialSum)):
            # pad = (0,)*other.ndim
            return type(self)(
                {
                    # tuple(
                    #     pad + k
                    #     for k in k_tup
                    # )
                    k_tup: other.scale(v) if nput.is_numeric(v) else other.mul_simple(v)
                    for k_tup, v in self.terms.items()
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
            return other.mul_simple(self)

class PTTensorCoeffProductSum(TensorCoefficientPoly, PolynomialInterface):
    """
    A representation for a sum of tensor coefficients * polynomial sums
    which is primarily here so we can transpose tensor coefficients intelligently
    """

    tensor_ids = ["V", "G", "U", "Z", "u", "M"] # we explicitly split Coriolis and Watson terms out

    def __init__(self, terms, prefactor=1, canonicalize=True, ndim=None, inds_map=None, reduced=False):
        super().__init__(terms, prefactor=prefactor, canonicalize=canonicalize)
        self._inds_map = {} if inds_map is None else inds_map
        self.reduced = reduced
        self._ndim = ndim
        # if canonicalize:
        #     self.terms = self.ensure_subpoly_dim(self.terms)
        # self.audit()

    def mutate(self,
               terms:'Any'=default, *, prefactor:'Any'=default,
               ndim:'Any'=default, inds_map:'Any'=default,
               canonicalize:'Any'=default, reduced:'Any'=default
               ):
        return type(self)(
            self.terms if terms is default else terms,
            canonicalize=(terms is not default) if canonicalize is default else canonicalize,
            prefactor=self.prefactor if prefactor is default else prefactor,
            reduced=self.reduced if reduced is default else reduced,
            ndim=self._ndim if ndim is default else ndim,
            inds_map=self._inds_map if inds_map is default else inds_map
        )

    @classmethod
    def format_key(cls, key):
        return "".join([
            "{}[{}]{}".format(
                cls.tensor_ids[k[1]],
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

    @classmethod
    def format_tensor_key(cls, key):
        return "".join([
            "{}[{}][{}]".format(
                cls.tensor_ids[k[1]],  # we explicitly split Coriolis and Watson terms out
                k[0],
                ", ".join(str(i) for i in k[2:])
            )
            for k in key
        ])
    def format_expr(self):
        """
        Formats in a Mathematica-ingestible format
        :return:
        """
        keys = []
        substrings = []
        for k,esum in self.terms.items():
            keys.append(self.format_tensor_key(k))
            substrings.append(esum.format_expr())

        join = " + " if sum(len(x) for x in substrings) < 80 else "\n  + "
        key_prods = [
            (
                "{}*({})"
                    if "+" in s else
                "{}*{}"
                    if s != "" else
                "{}"
            ).format(k,s)
            for k,s in zip(keys, substrings)
        ]
        return join.join(key_prods)

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
        return self.mutate(
            {
                t: p
                for t, p in self.terms.items()
                if all(tt[:2] not in ops for tt in t)
            }
        )

    @property
    def ndim(self):
        if self._ndim is None:
            self._ndim = max(len(self.get_inds(coeff_inds)) for coeff_inds in self.terms.keys())
        return self._ndim
    def audit(self, target=None, ignore_constants=True):
        """
        Checks to ensure that the number of dimensions aligns with
        the number of indices in the tensor coefficients

        :return:
        """
        for coeff_indices,poly_stuff in self.terms.items():
            num_inds = len(self.get_inds(coeff_indices))
            if target is not None and num_inds < target: raise ValueError("too few indices ({}) for changes ({}) in {}".format(
                num_inds, target, poly_stuff
            ))
            try:
                poly_stuff.audit(target=num_inds, ignore_constants=ignore_constants)
            except ValueError:
                raise ValueError("mismatch between inds {} and poly {}".format(
                    coeff_indices,
                    poly_stuff
                )) from None
            # self._audit_dimension(poly_stuff, num_inds)
    def ensure_dimension(self, ndim):
        nd = self.ndim
        if nd >= ndim: return self

        return self.mutate(
            {k:p.ensure_dimension(ndim) for k,p in self.terms.items()},
            ndim=ndim,
        )
    @classmethod
    def ensure_subpoly_dim(cls, terms):
        new_terms = {}
        for k, p in terms.items():
            new_terms[k] = p.ensure_dimension(len(cls.coeff_product_inds(k)))
        return new_terms

    def sort(self):
        return self.mutate(
            {
                k: self.terms[k] for k in sorted(self.terms.keys())
            }
        )

    def permute(self, new_inds, check_perm=True):
        new_terms = {}
        for coeff_inds,polys in self.terms.items():
            new_key = tuple(
                ci[:2] + tuple(new_inds[j] if j < len(new_inds) else j for j in ci[2:])
                for ci in coeff_inds
            )
            new_terms[new_key] = polys.permute(new_inds, check_perm=check_perm)
        return self.mutate(new_terms)

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
            new_terms[k1] = self.terms[k1] + self.terms[k2].permute(merges[(k1, k2)])
        for k in unmerged:
            new_terms[k] = self.terms[k]

        return new_terms

    @classmethod
    def reindex_terms(cls, terms):
        new_terms = {}
        for k,v in terms.items():
            ind_vars = [iii for kk in k for iii in kk[2:]]
            if len(ind_vars) == 0:
                new_terms[k] = v
                continue
            perm = np.unique(ind_vars)
            if len(perm) == perm[-1] + 1:
                new_terms[k] = v
                continue
            submap = np.argsort(perm) # just to remap things down
            mapping = dict(zip(perm, submap))
            new_key = tuple(kk[:2] + tuple(mapping[iii] for iii in kk[2:]) for kk in k)
            new_terms[new_key] = v.permute(perm)
        return new_terms
    def combine(self, combine_coeffs=False, combine_subterms=True, combine_energies=True):
        if self.reduced: return self
        new_terms = {}

        if combine_coeffs:
            raise NotImplementedError("coeff combining too subtle for current evaluator strategy")

        base_terms = self.combine_terms() if combine_coeffs else self.terms
        base_terms = self.reindex_terms(base_terms)
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
        return self.mutate(new_terms, reduced=True)

    def shift(self, shift):
        new_terms = {}
        for k,p in self.terms.items():
            new_terms[k] = p.shift(shift)
        return self.mutate(new_terms)
    def scale(self, scaling):
        if nput.is_numeric(scaling) and scaling == 1: return self
        return self.mutate(
            prefactor=self.scaling * scaling
        )
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

    def _generate_direct_product_values(self,
                                        inds, remainder, min_dim, baseline,
                                        k1, k2, poly_1, poly_2
                                        ):

        left_inds = self.get_inds(k1)
        right_inds = self.get_inds(k2)
        num_left = len(left_inds)
        num_right = len(right_inds)
        num_fixed = len(inds[1])

        # we multiply indices on the left with the corresponding indices on the right
        left_fixed_inds, right_fixed_inds = [tuple(x) for x in inds]
        left_remainder_inds, right_remainder_inds = [tuple(y) for y in remainder]
        test_num_left = len(left_remainder_inds) + len(left_fixed_inds)
        if test_num_left < num_left:
            left_remainder_inds = left_remainder_inds + tuple(range(test_num_left, num_left))
        test_num_right = len(right_remainder_inds) + len(right_fixed_inds)
        if test_num_right < num_right:
            right_remainder_inds = right_remainder_inds + tuple(range(test_num_right, num_right))

        base_dim = num_left + num_right - num_fixed

        # we take every possible combo of remaining indices from left and right
        max_mult = min(num_right - num_fixed, num_left - num_fixed) # the max number of extra axes to align
        if min_dim is not None:
            max_mult = min(base_dim - min_dim, max_mult)
        if _DEBUG_PRINT:
            print("||"*20, base_dim, inds, k1, k2, "|", max_mult)
        for mul_size in range(max_mult+1):
            total_dimension = base_dim - mul_size
            if mul_size == 0: # basically just a direct product, every remaining index gets duped

                if _DEBUG_PRINT:
                    print("||||||0", inds, remainder)
                    print("|1", poly_1.format_expr())
                    print("|2", poly_2.format_expr())
                new_poly = poly_1.mul_along(poly_2, inds, remainder)
                if _DEBUG_PRINT:
                    print("|", new_poly.format_expr())
                    print("|||")

                left_mapping = {k: i for i, k in enumerate(left_fixed_inds)}
                left_rem_rem = left_remainder_inds
                for j, k in enumerate(left_rem_rem):
                    left_mapping[k] = j + num_fixed + mul_size
                right_mapping = {k: i for i, k in enumerate(right_fixed_inds)}
                right_rem_rem = right_remainder_inds
                for j, k in enumerate(right_rem_rem):
                    right_mapping[k] = j + num_left

                new_key = tuple(
                    t[:2] + tuple(left_mapping.get(k, k) for k in t[2:])
                    for t in k1
                ) + tuple(
                    t[:2] + tuple(right_mapping.get(k, k) for k in t[2:])
                    for t in k2
                )

                true_index = (
                        [0] * num_fixed  # the first indices don't need to be reordered
                        + [1] * len(remainder[0])
                        + [3] * (num_left - num_fixed - len(remainder[0]))
                        + [2] * len(remainder[1])
                        + [3] * (num_right - num_fixed - len(remainder[1]))
                )
                perm = np.argsort(true_index)  # gives the position each final index should map to
                # print("="*30)
                # print(inds, remainder, num_left, num_right)
                # print(poly_1)
                # print(poly_2)
                # print(new_poly)
                # print("...=..."*10)
                new_poly = new_poly.permute(perm)
                old_key = new_key
                new_key = tuple(
                    ci[:2] + tuple(perm[i] for i in ci[2:])
                    for ci in new_key
                )
                if (
                        np.any(np.not_equal(self.get_inds(new_key), self.get_inds(old_key)))
                        or len(self.get_inds(new_key)) != num_left + num_right - num_fixed
                ):
                    raise ValueError(old_key, new_key)
                new_poly.audit(num_left + num_right - num_fixed)

                yield new_key, new_poly#.ensure_dimension(total_dimension)
            else:
                # computed once and cached
                combo_inds_left = PerturbationTheoryTermProduct.get_combination_inds(len(left_remainder_inds), mul_size)
                combo_rem_left = PerturbationTheoryTermProduct.get_combination_comp(len(left_remainder_inds), mul_size)
                combo_inds_right = PerturbationTheoryTermProduct.get_combination_inds(len(right_remainder_inds), mul_size)
                combo_rem_right = PerturbationTheoryTermProduct.get_combination_comp(len(right_remainder_inds), mul_size)

                for left_choice_x, left_rem_x in zip(combo_inds_left, combo_rem_left):
                    left_choice = tuple(left_remainder_inds[k] for k in left_choice_x)
                    left_mul_inds = left_fixed_inds + left_choice

                    # now we build our new key, noting that the multiplied inds will be the first _n_
                    # inds in the new set, then the left inds, then the right ones
                    left_mapping = {k: i for i, k in enumerate(left_mul_inds)}
                    left_rem_rem = [left_remainder_inds[x] for x in left_rem_x]
                    for j, k in enumerate(left_rem_rem):
                        left_mapping[k] = j + num_fixed + mul_size
                    left_key = tuple(
                                t[:2] + tuple(left_mapping[k] for k in t[2:])
                                for t in k1
                            )

                    if baseline is not None:
                        # sometimes unchanged inds on the right get multiplied by changed left inds
                        left_baseline = [
                            baseline[x] if x < len(baseline) else 0
                            for x in left_mul_inds
                        ]
                    else:
                        left_baseline = None

                    for right_choice_x, right_rem_x in zip(combo_inds_right, combo_rem_right):
                        if any(
                            i < len(remainder[0]) and j < len(remainder[1])
                            for i,j in zip(left_choice_x, right_choice_x)
                        ): continue # we can't multiply two remainder indices

                        right_choice = tuple(right_remainder_inds[k] for k in right_choice_x)
                        for right_perm in itertools.permutations(right_choice):
                            right_mul_inds = right_fixed_inds + right_perm

                            right_mapping = {k: i for i, k in enumerate(right_mul_inds)}
                            right_rem_rem = [right_remainder_inds[x] for x in right_rem_x]
                            for j,k in enumerate(right_rem_rem):
                                right_mapping[k] = j + num_left
                            right_key = tuple(
                                t[:2] + tuple(right_mapping[k] for k in t[2:])
                                for t in k2
                            )

                            new_key = left_key + right_key

                            mul_inds = [
                                left_mul_inds, # self
                                right_mul_inds  # other
                            ]
                            mul_remainder = [
                                left_rem_rem,
                                right_rem_rem
                            ]

                            if _DEBUG_PRINT:
                                print("||||||m", mul_inds, mul_size)
                                print("|1", poly_1.format_expr())
                            if baseline is not None:
                                right_shift = [0] * num_right
                                for s,i in zip(left_baseline, right_mul_inds):
                                    right_shift[i] = s
                                shift_poly_2 = poly_2
                                if any(s != 0 for s in right_shift):
                                    shift_poly_2 = shift_poly_2.shift(right_shift)
                                if _DEBUG_PRINT:
                                    print("|2", poly_2.format_expr())
                                    print("|b", baseline, right_shift)
                                    print("|2b", shift_poly_2.format_expr())
                                new_poly = poly_1.mul_along(shift_poly_2, mul_inds, mul_remainder)
                            else:
                                if _DEBUG_PRINT:
                                    print("|2", poly_2.format_expr())
                                new_poly = poly_1.mul_along(poly_2, mul_inds, mul_remainder)
                            if _DEBUG_PRINT:
                                print("||", new_poly.format_expr())

                            # finally, we need to make sure that the indices are in the right order
                            # by noting that the remainder inds that were specified must be ordered
                            # after the proper inds
                            # so we build a permutation to apply to new_poly by constructing a proxy index
                            # ordered like the indices of new_poly but with values from the target position indices
                            true_index = (
                                [0] * num_fixed # the first indices don't need to be reordered
                                + [
                                    1 if i < len(remainder[0]) else # left remainder before right remainder
                                    2 if j < len(remainder[1]) else
                                    3
                                    for i,j in zip(left_choice_x, right_choice_x)
                                ]
                                + [
                                    1 if i < len(remainder[0]) else
                                    3
                                    for i in left_rem_x
                                ]
                                + [
                                    2 if j < len(remainder[1]) else
                                    3
                                    for i in right_rem_x
                                ]
                            )
                            perm = np.argsort(true_index) # gives the position each final index should map to
                            new_poly = new_poly.permute(perm)
                            old_key = new_key
                            new_key = tuple(
                                ci[:2] + tuple(perm[i] for i in ci[2:])
                                for ci in new_key
                            )
                            if (
                                    np.any(np.not_equal(self.get_inds(new_key), self.get_inds(old_key)))
                                    or len(self.get_inds(new_key)) != num_left + num_right - num_fixed - mul_size
                            ):
                                raise ValueError(old_key, new_key)
                            new_poly.audit(num_left + num_right - num_fixed - mul_size)

                            yield new_key, new_poly#.ensure_dimension(total_dimension)

    def _adjust_key_right(self, k, other, inds, remainder):
        # TODO: make sure other isn't fucked up somehow
        mapping = np.argsort(np.concatenate([inds[0], remainder[0]]))
        return tuple(
            ci[:2] + tuple(
                mapping[i] if i < len(mapping) else i + len(remainder[1])
                for i in ci[2:]
            )
            for ci in k
        )
    def _adjust_key_left(self, k, other, inds, remainder):
        raise NotImplementedError(...)
        # TODO: make sure other isn't fucked up somehow
        mapping = np.argsort(np.concatenate([inds[1], remainder[1]]))
        return tuple(
            ci[:2] + tuple(mapping[i] if i < len(mapping) else i for i in ci[2:])
            for ci in k
        )
    def mul_along(self, other:'PolynomialInterface', inds, remainder=None, min_dim=None, mapping=None,
                  baseline=None):
        """
        We multiply every subpoly along the given indices, transposing the appropriate tensor indices

        :param other:
        :param inds:
        :param remainder:
        :return:
        """
        if isinstance(other, PTTensorCoeffProductSum):
            # we need to take the direct product of aligning (or not)
            # for all of the untouched inds
            new = self.direct_multiproduct(other,
                                           lambda *kargs: self._generate_direct_product_values(inds, remainder,
                                                                                               min_dim, baseline,
                                                                                               *kargs
                                                                                               )
                                           )
        elif isinstance(other, (PTEnergyChangeProductSum, ProductPTPolynomial)):
            new = self.mutate(
                {
                    self._adjust_key_right(k, other, inds, remainder):other.rmul_along(p, inds, remainder=remainder, mapping=mapping)
                    for k,p in self.terms.items()
                },
                canonicalize=True
            )
        else:
            raise NotImplementedError(type(other))
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
            return self.mutate(
                {
                    self._adjust_key_left(k, other, inds, remainder):other.mul_along(p, inds, remainder=remainder, mapping=mapping)
                    for k,p in self.terms.items()
                },
                canonicalize=True
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
                                       # key_func=lambda k1,k2:k1+tuple(
                                       #     k[:2] + tuple(i+self.ndim for i in k[2:])
                                       #     for k in k2
                                       # ),
                                       mul=lambda a,b:a.mul_simple(b)
                                       )
        elif isinstance(other, PTEnergyChangeProductSum):
            return self.mutate(
                {
                    k: v.mul_simple(other.ensure_dimension(len(self.get_inds(k))))
                    for k, v in self.terms.items()
                },
                canonicalize=False
            )
        elif isinstance(other, (ProductPTPolynomial, ProductPTPolynomialSum)):
            return self.mutate(
                {
                    k: v.mul_simple(other)
                    for k, v in self.terms.items()
                },
                canonicalize=False
            )
        else:
            return other.rmul_simple(self)
            # raise NotImplementedError(type(other))
    def rmul_simple(self, other: 'PolynomialInterface') -> 'PolynomialInterface':
        if isinstance(other, PTEnergyChangeProductSum):
            return self.mutate(
                {
                    # tuple(
                    #     k[:2] + tuple(i + other.ndim for i in k[2:])
                    #     for k in k_tup
                    # )
                    k_tup: other.ensure_dimension(len(self.get_inds(k_tup)).mul_simple(v))
                    for k_tup, v in self.terms.items()
                },
                canonicalize=False
            )
        elif isinstance(other, (ProductPTPolynomial, ProductPTPolynomialSum)):
            return self.mutate(
                {
                    # tuple(
                    #     k[:2] + tuple(i + other.ndim for i in k[2:])
                    #     for k in k_tup
                    # )
                    k_tup: other.scale(v) if nput.is_numeric(v) else other.mul_simple(v)
                    for k_tup, v in self.terms.items()
                },
                canonicalize=False
            )
        elif nput.is_numeric(other):
            if other == 1:
                return self
            elif other == 0:
                return 0
            else:
                raise ValueError(other)
        else:
            return other.mul_simple(self)

class SqrtChangePoly(PolynomialInterface):
    def __init__(self, poly_obj:'PolynomialInterface', change, shift, canonicalize=False):
        if nput.is_numeric(poly_obj): raise ValueError("{} isn't a polynomial".format(poly_obj))
        if isinstance(poly_obj, SqrtChangePoly): raise ValueError("sqrt changes shouldn't stack...")
        self.poly_obj = poly_obj

        self.poly_change = np.asanyarray(change, dtype=int)
        # if len(self.poly_change) > 0 and np.count_nonzero(self.poly_change) != len(self.poly_change):
        #     raise ValueError("no zero changes allowed ({})".format(
        #         self.poly_change)
        #     )
        if shift is None: shift = np.zeros(len(self.poly_change), dtype=int)
        self.shift_start = np.asanyarray(shift)
        # if len(self.shift_start) != len(self.poly_change): raise ValueError(self.shift_start, self.poly_change)
        if len(self.shift_start) < len(self.poly_change):
            self.shift_start = np.pad(self.shift_start, [0, len(self.poly_change) - len(self.shift_start)])
        if canonicalize:
            self.poly_obj, self.poly_change, self.shift_start = self.canonicalize(
                self.poly_obj, self.poly_change, self.shift_start
            )

    def mutate(self,
               poly_obj:'Any'=default,
               poly_change:'Any'=default,
               shift_start:'Any'=default,
               canonicalize=False
               ):
        return type(self)(
            self.poly_obj if poly_obj is default else poly_obj,
            self.poly_change if poly_change is default else poly_change,
            self.shift_start if shift_start is default else shift_start,
            canonicalize=canonicalize
        )

    def strip(self):
        return self.mutate(
            poly_change=[c for c in self.poly_change if c != 0],
            shift_start=[s for c,s in zip(self.poly_change, self.shift_start) if c != 0]
        )

    @property
    def ndim(self):
        return len(self.poly_change)
    def audit(self, target=None, ignore_constants=True):
        try:
            return self.poly_obj.audit(target=self.ndim, ignore_constants=ignore_constants)
        except ValueError:
            raise ValueError("bad poly found in {}".format(self.short_repr()))
    def ensure_dimension(self, ndim):
        if self.ndim >= ndim: return self
        poly_change = np.pad(self.poly_change, [0, ndim - len(self.poly_change)])
        shift_start = np.pad(self.shift_start, [0, ndim - len(self.shift_start)])

        return self.mutate(
            self.poly_obj.ensure_dimension(ndim),
            poly_change,
            shift_start,
            canonicalize=False
        )

    def __repr__(self):
        return "{}({}, {})".format("SqS", self.poly_change, self.poly_obj)
    def short_repr(self):
        return "{}({}, <{}>)".format("SqS", self.poly_change, type(self.poly_obj).__name__)
    def format_sqrt_expr(self):
        subexprs = []
        for i, (s,c) in enumerate(zip(self.shift_start, self.poly_change)):
            if c == 0:
                subexprs.append("")
                continue
            if s == 0:
                subexprs.append(
                    "S[n_{idx},{delta}]".format(idx=i, delta=c, shift=s)
                )
            else:
                subexprs.append(
                    "S[n_{idx}{sign}{shift},{delta}]".format(idx=i, delta=c, sign="-" if s < 0 else "+", shift=abs(s))
                )
        return "".join(subexprs)
    def format_expr(self):
        """
        Formats in a Mathematica-ingestible format
        :return:
        """
        base_expr = self.poly_obj.format_expr()
        sqrt_exprs = self.format_sqrt_expr()
        if len(sqrt_exprs) == 0:
            return base_expr
        elif "+" in base_expr:
            return "{}*({})".format(sqrt_exprs, base_expr)
        else:
            return "{}*{}".format(sqrt_exprs, base_expr)

    def __add__(self, other):
        if nput.is_numeric(other):
            if other == 0:
                return self
            else:
                raise NotImplementedError(...)
                return self + ProductPTPolynomial([[other]])
        elif isinstance(other, SqrtChangePoly):
            self, other = self.align_dimensions(other)
            if np.any(other.poly_change != self.poly_change):
                raise ValueError("can't add polys with different changes {} vs. {}".format(
                    self.poly_change,
                    other.poly_change
                ))
            if np.any(self.shift_start != other.shift_start):
                raise ValueError("can't add polys with different starts {} vs. {}".format(
                    self.shift_start,
                    other.shift_start
                ))
            new_poly = self.poly_obj + other.poly_obj
            if nput.is_numeric(new_poly):
                return new_poly
            else:
                return self.mutate(new_poly)
        else:
            raise NotImplementedError("not sure how to add {} and {}".format(self, other))

    def combine(self, combine_coeffs=False, combine_subterms=True, combine_energies=False):
        if nput.is_numeric(self.poly_obj): return self.poly_obj
        return self.mutate(
            self.poly_obj.combine(
                combine_coeffs=combine_coeffs, combine_subterms=combine_subterms,
                combine_energies=combine_energies
            )
        )
    def sort(self):
        if nput.is_numeric(self.poly_obj): return self.poly_obj
        return self.mutate(
            self.poly_obj.sort()
        )
    # def condense(self, inds=None, return_inds=False, check_inds=True):
    #     if inds is None: inds = np.arange(len(self.poly_change))
    #     inds = np.asanyarray(inds)
    #     self.poly_change = np.asanyarray(self.poly_change)
    #     if check_inds:
    #         inds = inds[self.poly_change[inds] == 0]
    #     condensed, new_obj = self.poly_obj.condense(inds, return_inds=True)
    #     new_change = np.delete(self.poly_change, inds)
    #     new_shift = np.delete(self.shift_start, inds)
    #     new_poly = self.mutate(
    #         new_obj,
    #         new_change,
    #         new_shift
    #     )
    #     if return_inds:
    #         return condensed, new_poly
    #     else:
    #         return new_poly

    @classmethod
    def canonicalize(cls, poly_obj, poly_change, shift_start):
        pc = np.asanyarray(poly_change)
        nonzero_pos = np.nonzero(poly_change)
        if len(nonzero_pos) == 0 or len(nonzero_pos[0]) == 0:
            return poly_obj, np.array([], dtype=int), np.array([], dtype=int)
        nonzero_pos = nonzero_pos[0]

        # lex_sort = np.lexsort([
        #     np.sign(pc[nonzero_pos]),
        #     -np.abs(pc[nonzero_pos])
        # ])
        nzp = pc[nonzero_pos]
        lex_sort = PerturbationTheoryTerm.change_sort(nzp)
        sorting = nonzero_pos[lex_sort]
        # sorting = nonzero_pos[np.argsort(pc[nonzero_pos])]
        if (
                len(nonzero_pos) == len(poly_change)
                and np.all(sorting == np.arange(len(poly_change)))
        ):
            return poly_obj, poly_change, shift_start
        else:
            return poly_obj.permute(sorting), np.take(poly_change, sorting), np.take(shift_start, sorting)

    def shift(self, shift):
        # new_shift, new_poly = self.(self.poly_change, shift)
        # if not nput.is_numeric(new_poly):
        #     new_poly = self.poly_obj.mul_simple(new_poly)
        og = self.shift_start
        if len(shift) < len(og):
            shift = list(shift) + [0] * (len(og) - len(shift))
        elif len(shift) > len(og):
            og = list(og) + [0] * (len(shift) - len(og))

        return self.mutate(
            shift_start=[i+s for i,s in zip(og, shift)]
        )
    def scale(self, scaling):
        if nput.is_numeric(scaling) and scaling == 1: return self
        return self.mutate(
            self.poly_obj.scale(scaling)
        )
    def permute(self, perm, check_perm=True):

        if check_perm:
            if np.any(np.sort(perm) != np.arange(len(perm))):
                raise ValueError("bad permutation {}")

        return self.mutate(
            self.poly_obj.permute(perm, check_perm=check_perm),
            [self.poly_change[c] for c in perm],
            [
                self.shift_start[c] for c in perm
                if c < len(self.poly_change)
            ]
        )

    def __radd__(self, other):
        return self.__add__(other)
    def __mul__(self, other):
        raise NotImplementedError(...)
    def __rmul__(self, other):
        raise NotImplementedError(...)

    def get_change_poly(self,
                        initial_changes, extra_changes,
                        initial_shift, extra_shift
                        ):
        for ich, ish, ech, esh in zip(
                initial_changes, initial_shift,
                extra_changes, extra_shift
        ):
            if esh != ich + ish:
                raise ValueError("split shifts not supported; target shift/changes: {}/{}, initial shift/changes: {}/{}".format(
                    extra_shift, extra_changes, initial_shift, initial_changes
                ))
        new_changes = [i_sh + e_sh for i_sh, e_sh in zip(initial_changes, extra_changes)]
        sqrt_contrib = HarmonicOscillatorRaisingLoweringPolyTerms.get_direction_change_poly(
            extra_changes, initial_changes, initial_shift
        )
        if len(sqrt_contrib) == 0: sqrt_contrib = None
        if sqrt_contrib is not None:
            sqrt_contrib = ProductPTPolynomial(sqrt_contrib)

        return new_changes, sqrt_contrib

    def mul_simple(self, other:'PolynomialInterface'):
        if isinstance(other, SqrtChangePoly):
            return self.mul_along(other, list(range(self.ndim)))
            # return self.mutate(
            #     self.poly_obj.mul_simple(other.poly_obj),
            #     np.concatenate([self.poly_change, other.poly_change]),
            #     np.concatenate([self.shift_start, other.shift_start])
            # )
        else:
            return self.mutate(
                self.poly_obj.mul_simple(other),
                # np.pad(self.poly_change, [0, other.ndim]),
                self.poly_change,
                # np.pad(self.shift_start, [0, other.ndim])
                self.shift_start
            )
    def rmul_simple(self, other: 'PolynomialInterface') -> 'PolynomialInterface':
        if isinstance(other, SqrtChangePoly):
            other.mul_simple(self)
        else:
            return self.mutate(
                other.mul_simple(self.poly_obj),
                self.poly_change,
                # np.pad(self.poly_change, [other.ndim, 0]),
                self.shift_start
                # np.pad(self.shift_start, [other.ndim, 0])
            )
    def mul_along(self, other, inds, remainder=None, mapping=None, baseline=None):
        if isinstance(other, SqrtChangePoly):
            # self, other = self.align_dimensions(other)
            # if len(self.poly_change) != len(other.poly_change): raise ValueError(
            #     "mismatch between poly changes {} vs. {}".format(
            #         self.poly_change, other.poly_change
            #     )
            # )
            # if len(self.shift_start) != len(other.shift_start): raise ValueError(
            #     "mismatch between poly shifts {} vs. {}".format(
            #         self.shift_start, other.shift_start
            #     )
            # )

            (self_inds, other_inds) = inds
            (self_remainder, other_remainder) = remainder
            #     = ProductPTPolynomial.get_index_mapping(
            #     self.ndim,
            #     other.ndim,
            #     inds,
            #     return_remainder=True
            # )
            # baseline = self._get_mul_baseline()


            self_shift = [
                self.shift_start[c] if len(self.shift_start) > c else 0
                for c in self_inds
            ]
            other_shift = [
                other.shift_start[c] if len(other.shift_start) > c else 0
                for c in other_inds
            ]
            self_change =[
                self.poly_change[c] if len(self.poly_change) > c else 0
                for c in self_inds
            ]
            other_change = [
                other.poly_change[c] if len(other.poly_change) > c else 0
                for c in other_inds
            ]

            new_changes, sqrt_contrib = self.get_change_poly(
                self_change, other_change,
                self_shift, other_shift
            )

            final_changes = [
                                n for n in new_changes #if n != 0
                            ] + [
                                self.poly_change[r] for r in self_remainder
                                if r < len(self.poly_change)
                            ] + [
                                other.poly_change[r] for r in other_remainder
                                if r < len(other.poly_change)
                            ]

            final_shift = [
                              s for s, n in zip(self_shift, new_changes) #if n != 0
                          ] + [
                              self.shift_start[r] for r in self_remainder
                              if r < len(self.shift_start)
                          ] + [
                              other.shift_start[r] for r in other_remainder
                              if r < len(other.shift_start)
                          ]

            other_poly = other.poly_obj
            if other_shift != self_shift:
                # need to shift other to align with self
                delta_shift = [o-s for s,o in zip(self_shift, other_shift)]
                total_shift = [0] * len(other.shift_start)
                for d,c in zip(delta_shift, other_inds):
                    total_shift[c] = d
                other_poly = other_poly.shift(total_shift)
            if sqrt_contrib is not None:
                new_self = self.poly_obj.mul_along(sqrt_contrib,
                                                   [self_inds, np.arange(len(self_inds))],
                                                   remainder=[self_remainder, np.arange(len(self_inds), sqrt_contrib.ndim)]
                                                   ) #TODO: optimize out remainder calc.
                if baseline is not None:
                    baseline = [
                        baseline[x] for x in itertools.chain(self_inds, self_remainder)
                    ] + list(
                        baseline[len(self_inds) + len(self_remainder):]
                    )
                # inds and remainders must be adapted for the fact that we've moved around indices in doing this...
                new_poly = new_self.mul_along(
                    other_poly,
                    [list(range(len(self_inds))), other_inds],
                    min_dim=len(final_changes),
                    remainder=[
                        np.arange(len(self_inds), len(self_remainder) + len(self_inds)),
                        other_remainder
                    ],
                    baseline=baseline
                )
            else:
                # -print("md...", len(final_changes))
                new_poly = self.poly_obj.mul_along(
                    other_poly,
                    [self_inds, other_inds],
                    remainder=[self_remainder, other_remainder],
                    min_dim=len(final_changes),
                    baseline=baseline
                )

            return self.mutate(
                new_poly,
                final_changes,
                final_shift
            )
        else:
            # raise NotImplementedError(type(other))
            return self.mutate(
                self.poly_obj.mul_along(other, inds, remainder=remainder, mapping=mapping),
                self.poly_change,
                self.shift_start
            )
    def rmul_along(self, other, inds, remainder=None, mapping=None):
        if isinstance(other, SqrtChangePoly):
            return other.mul_along(self, inds, remainder=remainder, mapping=mapping)
        else:
            raise NotImplementedError(type(other))
            return type(self)(
                self.poly_obj.rmul_along(other, inds, remainder=remainder, mapping=mapping),
                self.poly_change,
                self.shift_start
            )

class PerturbationTheoryTerm(metaclass=abc.ABCMeta):
    """
    A generic version of one of the three terms in
    PT that will generate a correction polynomial
    """
    def __init__(self, logger=None, checkpoint=None,
                 allowed_terms=None,
                 allowed_changes=None,
                 intermediate_normalization=False):
        self._exprs = None
        self._raw_changes = {}
        self._changes = None
        self._cache = {}
        self.logger = Logger.lookup(logger)
        self.checkpoint = Checkpointer.build_canonical(checkpoint)
        self.intermediate_normalization = intermediate_normalization
        self.allowed_terms = allowed_terms
        self.allowed_changes = allowed_changes

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

    @staticmethod
    def change_sort_key(changes):
        return -(10 * np.abs(changes) - (changes < 0))
    @classmethod
    def change_sort(cls, changes):
        return np.argsort(cls.change_sort_key(changes), axis=-1)
    @abc.abstractmethod
    def get_changes(self) -> 'dict[tuple[int], Any]':
        ...
    @property
    def changes(self):
        if self._changes is None:
            self._changes = self.get_changes()
            if self.allowed_changes is not None:
                acs = set(tuple(k) for k in self.allowed_changes)
                self._changes = {
                    k:v for k,v in self._changes.items()
                    if k in acs
                }
        return self._changes

    def get_serializer_key(self): # to be overridden
        return None
    @property
    def serializer_key(self):
        return self.get_serializer_key()

    def get_poly_terms(self, changes, simplify=True, shift=None) -> 'SqrtChangePoly':
        with self.logger.block(tag="Terms for {}[{}]".format(self, changes)):
            start = time.time()
            changes = tuple(c for c in changes if c != 0)
            cache_key = changes if shift is None or all(s == 0 for s in shift) else (changes, tuple(shift))
            terms = self.changes.get(cache_key, None)
            if not isinstance(terms, SqrtChangePoly):
                if nput.is_numeric(terms): return terms

                with self.checkpoint as chk:
                    key = self.serializer_key
                    try:
                        if key is None:
                            terms = None
                        else:
                            terms = chk[key][changes]
                    except KeyError:
                        terms = None
                    if terms is None:
                        # terms = self._raw_changes.get(changes, None)
                        # if not isinstance(terms, (PTTensorCoeffProductSum, int, np.integer)):
                        terms = sum(
                            term.get_poly_terms(changes, shift=None)
                            for term in self.expressions
                        )
                        # self._raw_changes[changes] = terms
                        if simplify and not nput.is_zero(terms):
                            # print("...simplifying...", end = None)
                            terms = terms.combine()#.sort()
                            # print("done")
                        # if simplify:
                        self.changes[changes] = terms
                        if key is not None:
                            try:
                                base = chk[key]
                            except KeyError:
                                base = {}
                            base[changes] = terms
                            chk[key] = base
                if shift is not None and len(shift) > 0 and not nput.is_numeric(terms):
                    if not isinstance(terms, SqrtChangePoly): raise ValueError(...)
                    terms = terms.shift(shift)
                    self.changes[cache_key] = terms
            end = time.time()
            # self.logger.log_print('terms: {t}', t=terms)
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
            change_key = changes if (shift is None or all(s==0 for s in shift)) else (changes, shift)
            terms = self.changes.get(change_key, None)
            if not isinstance(terms, SqrtChangePoly):
                if nput.is_numeric(terms): return terms

                terms = self.changes.get(changes, 0)
                if not (isinstance(terms, SqrtChangePoly) or nput.is_numeric(terms)):
                    terms = self.get_poly_terms(changes, shift=None)
                    if simplify and not nput.is_zero(terms):
                        terms = terms.combine().sort()
                    self.changes[changes] = terms
                if shift is not None:
                    raise NotImplementedError("reshifting not supported yet...")
                    self.changes[change_key] = terms
            if return_evaluator:
                return PerturbationTheoryExpressionEvaluator(self, terms, changes, logger=self.logger)
            else:
                return terms

class OperatorExpansionTerm(PerturbationTheoryTerm):
    def __init__(self, terms, order=None, identities=None, symmetrizers=None, logger=None, checkpoint=None,
                 allowed_terms=None, allowed_changes=None, intermediate_normalization=False):
        super().__init__(logger=logger, checkpoint=checkpoint,
                         allowed_terms=allowed_terms,
                         allowed_changes=allowed_changes,
                         intermediate_normalization=intermediate_normalization)

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
        if nput.is_numeric(poly_contrib) and poly_contrib == 0:
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
                            sorting = np.argsort(sort_key, axis=-1)
                            # raise Exception(
                            #     sort_key.shape,
                            #     sorting.shape,
                            #     sort_key,
                            #     sorting,
                            #     sort_key[sorting].shape,
                            #     np.sort(sort_key).shape
                            # )
                            # _, uinds = np.unique(sort_key[sorting], axis=0, return_index=True)
                            # resigned_perms = resigned_perms[uinds]
                            # sorting = sorting[uinds]
                            # # resigned_perms = resigned_perms[np.argsort(sort_key, axis=1)]
                            # for subp, sort in zip(resigned_perms, sorting):
                            #     changes[tuple(subp[sort])] = None
                            #
                            # _, uinds = np.unique(np.sort(sort_key), axis=0, return_index=True)
                            # resigned_perms = resigned_perms[uinds]
                            for subp, ss in zip(resigned_perms, sorting):
                                changes[tuple(subp[ss])] = None
        return changes

    def get_subexpresions(self) -> 'Iterable[PerturbationTheoryTerm]':
        raise NotImplementedError("shouldn't need this here...")

    def get_poly_terms(self, changes, shift=None) -> 'SqrtChangePoly': #TODO: CACHE THIS SHIT
        #NOTE: this explicitly excludes the sqrt contribution to the poly
        with self.logger.block(tag="Terms for {}[{}]".format(self, changes)):
            start = time.time()
            changes = tuple(changes)
            shift = None if shift is None or all(s == 0 for s in shift) else tuple(shift)
            key = changes if shift is None else (changes, shift)
            terms = self.changes.get(key, None)
            if not isinstance(terms, SqrtChangePoly):
                if not nput.is_numeric(terms):
                    terms = self._get_poly_terms(changes, shift=shift)
                    terms.audit()
                    self.changes[key] = terms
            end = time.time()
            self.logger.log_print("took {e:.3f}s", e=end-start)
        return terms

    def _get_poly_terms(self, changes, shift=None) -> 'SqrtChangePoly':  # TODO: CACHE THIS SHIT

        poly_contribs = {}

        abs_change = np.array([abs(c) for c in changes])
        total_change = sum(abs_change)

        og_changes = changes
        og_shift = shift
        og_abs_change = abs_change
        og_dim = np.count_nonzero(og_abs_change)

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

            total_dim = og_dim + remainder // 2

            # now we enumerate all of the ways to partition `group_size` elements across `total_dim` dimensions
            # such that each `partition - changes` is even in each dimension
            abs_change = np.pad(og_abs_change, [0, total_dim - len(og_abs_change)])
            changes = np.pad(og_changes, [0, total_dim - len(og_changes)]).astype(int)
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

                        if isinstance(base_poly, SqrtChangePoly):
                            base_poly = base_poly.poly_obj
                        prefactor = (self.order, term_id) + symm_idx #type: tuple[int]
                        poly_contribs[(prefactor,)] = poly_contribs.get((prefactor,), 0) + base_poly

        new = PTTensorCoeffProductSum(poly_contribs, canonicalize=False)
        if og_shift is None: og_shift = [0] * len(og_changes)
        if len(og_shift) < len(og_changes): og_shift = np.pad(og_shift, [0, len(og_changes) - len(og_shift)])
        # new.audit()
        return SqrtChangePoly(new, og_changes, og_shift)

    _poly_cache = {}
    @classmethod
    def _resolve_poly(cls, term_list, partition_sizes, p_index, c_key, s_key, p_vec, changes, shift):
        # print(term_list, p_index, changes, shift, p_index)

        key = (term_list, partition_sizes, p_index, c_key)
        poly = cls._poly_cache.get(key, None)
        if poly is None:
            phase = (-1)**(sum(1 for t in term_list if t == 'p')//2)
            poly_coeffs = [
                cls._evaluate_poly_coeffs(term_list, inds, delta, 0)
                for inds, delta in zip(p_vec, changes)
            ]
            if any(nput.is_zero(c) for c in poly_coeffs):
                poly = 0
            else:
                poly_coeffs = [
                    c for c in poly_coeffs
                    if c[0] != 1 or len(c) > 1
                ]
                poly = ProductPTPolynomial(poly_coeffs, prefactor=phase, idx=key)
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
        final_change = tuple(c for c in final_change if c != 0)
        if len(final_change) == 0:
            return 0

        base_term = self.subterm.get_poly_terms(changes, shift=shift)
        if isinstance(base_term, PolynomialInterface):
            prefactor = PTEnergyChangeProductSum.monomial(changes, 1)
            base_term = base_term.mul_simple(prefactor)
            # base_term.audit()

        # if shift is not None and not nput.is_numeric(base_term):
        #     base_term = base_term.shift(shift)

        return base_term

class ShiftedHamiltonianCorrection(PerturbationTheoryTerm):
    """
       Adds the wave function correction and the overlap term
       """

    def __init__(self, parent, order, logger=None, checkpoint=None,
                       allowed_terms=None, allowed_changes=None, intermediate_normalization=False):
        super().__init__(logger=logger, checkpoint=checkpoint,
                         allowed_terms=allowed_terms, allowed_changes=allowed_changes,
                         intermediate_normalization=intermediate_normalization)

        self.parent = parent
        self.order = order

    def get_serializer_key(self):  # to be overridden
        return self.__repr__()

    def __repr__(self):
        return "DH[{}]".format(self.order)

    def get_changes(self):
        base_changes = {}
        for expr in self.expressions:
            for subchange in expr.changes:
                base_changes[subchange] = base_changes.get(subchange, [])
                base_changes[subchange].append(expr)  # just track which exprs generate the change
        base_changes[()] = None  # also the constant term
        return base_changes

    def get_subexpresions(self) -> 'Iterable[PerturbationTheoryTerm]':
        H = self.parent.hamiltonian_expansion
        E = self.parent.energy_correction
        k = self.order

        if k == 1:
            if len(H) == 0:
                return []
            else:
                return [H[k]]
        else:
            if len(H) > k:
                return [H[k], -E(k)]
            else:
                return [-E(k)]

class WavefunctionCorrection(PerturbationTheoryTerm):
    def __init__(self, parent, order, logger=None, checkpoint=None,
                 allowed_terms=None, allowed_changes=None, intermediate_normalization=False):
        super().__init__(logger=logger, checkpoint=checkpoint,
                 allowed_terms=allowed_terms, allowed_changes=allowed_changes,
                         intermediate_normalization=intermediate_normalization)

        self.parent = parent
        self.order = order


    def get_serializer_key(self): # to be overridden
        return self.__repr__()
    def __repr__(self):
        return "W[{}]".format(self.order)

    # def get_poly_terms(self, changes, shift=None) -> 'PTTensorCoeffProductSum':
    #     if len(changes) == 0:
    #         return self.parent.overlap_correction(self.order).get_poly_terms([], shift=shift)
    #     else:
    #         return super().get_poly_terms(changes, shift=shift)

    def get_changes(self):
        base_changes = {}
        for expr in self.expressions:
            for subchange in expr.changes:
                base_changes[subchange] = base_changes.get(subchange, [])
                base_changes[subchange].append(expr) # just track which exprs generate the change
        # base_changes[()] = None # also the constant term
        return base_changes

    def get_subexpresions(self):
        H = self.parent.shifted_hamiltonian_correction
        E = self.parent.energy_correction
        W = self.parent.wavefunction_correction
        k = self.order

        base_terms = [
            E(k-i) * W(i)
            for i in range(1, k)
        ] + [
            -PerturbationOperator.lookup(H(k-i) * W(i))
                if i > 0 else
            -PerturbationOperator.lookup(H(k))
            for i in range(0, k)
        ]

        return base_terms

class EnergyCorrection(PerturbationTheoryTerm):
    def __init__(self, parent, order, logger=None, checkpoint=None,
                 allowed_terms=None, allowed_changes=None,intermediate_normalization=False):
        super().__init__(logger=logger, checkpoint=checkpoint,
                         allowed_terms=allowed_terms, allowed_changes=allowed_changes,
                         intermediate_normalization=intermediate_normalization)

        self.parent = parent
        self.order = order

    def get_serializer_key(self):  # to be overridden
        return self.__repr__()
    def __repr__(self):
        return "E[{}]".format(self.order)

    def get_changes(self) -> 'dict[tuple[int], Any]':
        return {():None}

    def get_subexpresions(self) -> 'Iterable[PerturbationTheoryTerm]':
        H = self.parent.hamiltonian_expansion
        J = self.parent.shifted_hamiltonian_correction
        W = self.parent.full_wavefunction_correction
        k = self.order

        if k % 2 == 1:
            return []

        return ([H[k]] if len(H) > k else []) + [
            (J(k - i) * W(i))
            for i in range(1, k)
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

    def __call__(self, changes, shift=None, coeffs=None, freqs=None, check_sorting=None, simplify=True):
        if len(changes) > 0:
            raise ValueError("no")
        return super().__call__((), shift=shift, coeffs=coeffs, freqs=freqs, check_sorting=False, simplify=simplify)

class WavefunctionOverlapCorrection(PerturbationTheoryTerm):
    """
    Provides a slight optimization on the base `WavefunctionCorrection`
    """

    def __init__(self, parent, order, logger=None, checkpoint=None,
                 allowed_terms=None, allowed_changes=None, intermediate_normalization=False):
        super().__init__(logger=logger, checkpoint=checkpoint,
                         allowed_terms=allowed_terms, allowed_changes=allowed_changes,
                         intermediate_normalization=intermediate_normalization)

        self.parent = parent
        self.order = order

    def get_serializer_key(self): # to be overridden
        return self.__repr__()
    def __repr__(self):
        return "O[{}]".format(self.order)

    def get_changes(self) -> 'dict[tuple[int], Any]':
        return {():None}

    def get_subexpresions(self) -> 'Iterable[PerturbationTheoryTerm]':
        W = self.parent.full_wavefunction_correction
        k = self.order
        L = lambda o:FlippedPerturbationTheoryTerm.lookup(W(o))

        if k == 1:
            return []

        return [
            -1/2 * (L(k - i) * W(i)) for i in range(1, k)
        ]

class FullWavefunctionCorrection(PerturbationTheoryTerm):
    """
    Adds the wave function correction and the overlap term
    """

    def __init__(self, parent, order, logger=None, checkpoint=None,
                 allowed_terms=None, allowed_changes=None, intermediate_normalization=False):
        super().__init__(logger=logger, checkpoint=checkpoint,
                         allowed_terms=allowed_terms, allowed_changes=allowed_changes,
                         intermediate_normalization=intermediate_normalization)

        self.parent = parent
        self.order = order

    def get_serializer_key(self): # to be overridden
        return self.__repr__()
    def __repr__(self):
        return "Y[{}]".format(self.order)

    def get_changes(self):
        base_changes = {}
        for expr in self.expressions:
            for subchange in expr.changes:
                base_changes[subchange] = base_changes.get(subchange, [])
                base_changes[subchange].append(expr)  # just track which exprs generate the change
        base_changes[()] = None  # also the constant term
        return base_changes

    def get_subexpresions(self) -> 'Iterable[PerturbationTheoryTerm]':
        W = self.parent.wavefunction_correction
        O = self.parent.overlap_correction
        k = self.order

        if k == 1 or self.intermediate_normalization:
            return [W(k)]
        else:
            return [W(k), O(k)]

class OperatorCorrection(PerturbationTheoryTerm):

    def __init__(self, parent, order, operator_type=None, logger=None, checkpoint=None,
                 allowed_terms=None, allowed_changes=None, intermediate_normalization=False):
        super().__init__(logger=logger, checkpoint=checkpoint,
                         allowed_terms=allowed_terms, allowed_changes=allowed_changes,
                         intermediate_normalization=intermediate_normalization)

        self.parent = parent
        self.order = order
        self.expansion = parent.operator_expansion_terms(order, logger=self.logger, operator_type=operator_type)

    def get_serializer_key(self):  # to be overridden
        return self.__repr__()
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
        W = self.parent.full_wavefunction_correction
        M = self.expansion
        k = self.order
        L = lambda o:FlippedPerturbationTheoryTerm.lookup(W(o))

        terms = [
            (i, k - i - j, j)
            for i in range(0, k + 1)
            for j in range(0, k + 1 - i)
        ]
        if self.allowed_terms is not None:
            allowed_terms = {tuple(t) for t in self.allowed_terms}
            terms = [t for t in terms if t in allowed_terms]

        exprs = tuple(
            M[a]
                if i == 0 and j == 0 else
            L(j) * M[a] * W(i)
                if i > 0 and j > 0 else
            M[a] * W(i)
                if i > 0 else
            L(j) * M[a]

            for i, a, j in terms
        )

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
        if nput.is_numeric(self.prefactor) and self.prefactor == -1:
            return "-{}".format(self.base)
        else:
            return "{}*{}".format(self.prefactor, self.base)

    def get_changes(self) -> 'dict[tuple[int], Any]':
        return {k:None for k in self.base.changes}
    def get_poly_terms(self, changes, shift=None) -> 'SqrtChangePoly':
        base = self.base.get_poly_terms(changes, shift=shift)
        if nput.is_numeric(base): return base
        if isinstance(base, SqrtChangePoly):
            new = type(base)(self.prefactor * base.poly_obj, base.poly_change, base.shift_start)
        else:
            new = SqrtChangePoly(self.prefactor * base, changes, shift)
        return new

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

        # We can work through the many classes of changes just by seeing what generates [1, -1] from H2/H1.H1
        # 1. [1], [-1], [], [0, 1], [:]
        # 2. [-1], [1], [], [1, 0], [:]
        # 3. [2], [-1, -1], [0], [0, 1], [:]
        # 4. [-3], [2, 1], [0], [1, 0], [:]
        # 5. [2, 1], [-3], [0], [1, 0], [:]
        # 5. [1], [1, -1, -1], [1], [1, 2, 0], [:2]
        # 6. [1, 1, -1], [1, -1, -1], [[0, 2], [1, 0]], [2, 3, 0, 1]


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

                perm_ov = p[perm_inds2]  # I think I can actually take the unique perms...?

                if len(r) > 0:
                    newts = np.concatenate([
                        arr_1[np.newaxis] + perm_ov,
                        np.broadcast_to(r[np.newaxis], (len(perm_ov), len(r)))
                    ],
                        axis=1
                    )
                else:
                    newts = arr_1[np.newaxis] + perm_ov

                newt_transposes = cls.change_sort(newts)

                # build new polys
                from_axes = inv_perm[uperm_inds]  # axes in change_1 that are added to the corresponding
                                                  # parts of change_2 to build the whole change
                if len(io) < len(p):
                    io = np.concatenate([io, np.arange(len(change_2), len(change_2) + len(p))])
                target_axes = io[perm_inds[uperm_inds]]  # axes in change_2 that are added to the corresponding
                                                         # parts of change_1[io] to build the whole change
                # print(newts, change_1, change_2)
                for base_change, transp, targ_inds, src_inds in zip(newts, newt_transposes, target_axes, from_axes):
                    expected_change = tuple(
                        change_1[n] + (change_2[i] if i < len(change_2) else 0)
                        for n,i in enumerate(targ_inds)
                    ) + tuple(
                        change_1[len(targ_inds):]
                    ) + tuple(
                        change_2[i]
                        for i in range(len(change_2)) if i not in targ_inds
                    )
                    if tuple(base_change) != expected_change:
                        raise ValueError(
                            base_change, change_1, change_2, targ_inds, src_inds, p, r, io
                        )
                    yield base_change[transp], transp, targ_inds, src_inds
        # print("<::")
        # handle overlap == 0 term
        base_change = np.concatenate([change_1, change_2])
        transp = np.argsort(base_change)
        targ_inds = ()
        src_inds = ()
        yield base_change[transp], transp, targ_inds, src_inds

    @classmethod
    def get_change_classes(cls, changes):
        classes = []
        for k in changes:
            if len(k) == 2 and isinstance(k[0], tuple): continue

            n = len(k)
            if n + 1 > len(classes): classes.extend([[] for _ in range(n+1 - len(classes))])
            classes[n].append(k)
        return [np.array(ch) for ch in classes]

    # change_tup = collections.namedtuple([])

    _combination_inds = {}
    _combination_comps = {}
    @classmethod
    def get_combination_inds(cls, n, r):
        if r == 0:
            return np.array([[]], dtype=int)
        elif r == n:
            return np.arange(n)[np.newaxis]
        elif (n, r) not in cls._combination_inds:
            cls._combination_inds[(n, r)] = nput.combination_indices(n, r)
        return cls._combination_inds[(n, r)]
    @classmethod
    def get_combination_comp(cls, n, r):
        if r == 0:
            return np.arange(n)[np.newaxis]
        elif r == n:
            return np.array([[]], dtype=int)
        elif (n, r) not in cls._combination_comps:
            cls._combination_comps[(n, r)] = nput.index_complement(n, cls.get_combination_inds(n, r))[0]
        return cls._combination_comps[(n, r)]

    @staticmethod
    def _fill_change_data(full_changes, n, new_changes,
                          init_changes, subsorts, contacts, rems):
        if n not in full_changes:
            full_changes[n] = {
                'changes': [],
                'sorts': [],
                'init_changes': [[], []],
                'contacts': [[], []],
                'rems': [[], []]
            }

        change_data = full_changes[n]
        change_data['changes'].append(new_changes)
        change_data['sorts'].append(subsorts)
        for storage, contact in zip(change_data['init_changes'], init_changes):
            storage.append(contact)
        for storage, contact in zip(change_data['contacts'], contacts):
            storage.append(contact)
        for storage, contact in zip(change_data['rems'], rems):
            storage.append(contact)

    def get_changes(self):

        changes_1 = self.get_change_classes(self.gen1.changes)
        changes_2 = self.get_change_classes(self.gen2.changes)

        full_changes = {} # indexed by number of modes in final change

        for n1,ch1 in enumerate(changes_1):
            if len(ch1) == 0: continue
            for n2,ch2 in enumerate(changes_2):
                if len(ch2) == 0: continue

                if n1 == 0:
                    k = len(ch2)
                    self._fill_change_data(full_changes, n2, ch2,
                                           [np.broadcast_to(np.array([[]]), (k, 0)), ch2],
                                           np.broadcast_to(np.arange(n2)[np.newaxis], (k, n2)),
                                           [np.broadcast_to(np.array([[]]), (k, 0)), np.broadcast_to(np.array([[]]), (k, 0))],
                                           [np.broadcast_to(np.array([[]]), (k, 0)), np.broadcast_to(np.arange(n2)[np.newaxis], (k, n2))]
                                           )
                elif n2 == 0:
                    k = len(ch1)
                    self._fill_change_data(full_changes, n1, ch1,
                                           [ch1, np.broadcast_to(np.array([[]]), (k, 0))],
                                           np.broadcast_to(np.arange(n1)[np.newaxis], (k, n1)),
                                           [np.broadcast_to(np.array([[]]), (k, 0)), np.broadcast_to(np.array([[]]), (k, 0))],
                                           [np.broadcast_to(np.arange(n1)[np.newaxis], (k, n1)), np.broadcast_to(np.array([[]]), (k, 0))]
                                           )
                else:
                    for m in range(min(n1, n2)+1):
                        if m == 0:
                            inds_1 = np.array([[]])
                            inds_2 = np.array([[]])
                            rems_1 = np.arange(n1)[np.newaxis]
                            rems_2 = np.arange(n2)[np.newaxis]

                            uinds = None

                            # pure direct products, concatenate final combinations
                            new_combo = np.concatenate([
                                np.broadcast_to(ch1[:, np.newaxis, :], (ch1.shape[0], ch2.shape[0], n1)).reshape(-1, n1),
                                np.broadcast_to(ch2[np.newaxis, :, :], (ch1.shape[0], ch2.shape[0], n2)).reshape(-1, n2),
                            ], axis=1)
                        else:
                            inds_1 = self.get_combination_inds(n1, m)
                            inds_2 = self.get_combination_inds(n2, m)
                            rems_1 = self.get_combination_comp(n1, m)
                            rems_2 = self.get_combination_comp(n2, m)

                            # ch1 is a c1 x n1 array (c1 changes of len n1)
                            # we need to broadcast every change pair, so we get a
                            # c1 x k1 x n1 array where k is n1 choose m

                            # add overlaps, concatenate the remainder
                            ov_1 = nput.vector_take(ch1, inds_1) # c1 x k1 x m
                            ov_2 = nput.vector_take(ch2, inds_2) # c2 x k2 x m
                            re_1 = nput.vector_take(ch1, rems_1) # c1 x k1 x n1 - m
                            re_2 = nput.vector_take(ch2, rems_2) # c2 x k2 x n2 - m

                            # remove duplicate changes --> handled later on
                            uinds_1 = []
                            k1 = len(inds_1)
                            for n,block in enumerate(ov_1): # k1 x m blocks
                                _, uu = np.unique(block, axis=0, return_index=True)
                                uu = np.sort(uu)
                                uinds_1.append(k1*n + uu ) # need to add on the block offset
                            uinds_2 = []
                            k2 = len(inds_2)
                            for n,block in enumerate(ov_2): # k2 x m blocks
                                _, uu = np.unique(block, axis=0, return_index=True)
                                uu = np.sort(uu)
                                uinds_2.append(k2*n + uu ) # need to add on the block offset
                            uinds = [np.concatenate(uinds_1), np.concatenate(uinds_2)]

                            combos = ov_1[:, :, np.newaxis, np.newaxis] + ov_2[np.newaxis, np.newaxis, :, :]
                            cats = [combos.reshape(-1, m)]
                            if m < n1:
                                pads1 = np.broadcast_to(
                                    re_1[:, :, np.newaxis, np.newaxis],
                                    re_1.shape[:2] + re_2.shape[:2] + re_1.shape[-1:]
                                ).reshape(-1, n1 - m)
                                cats.append(pads1)
                            if m < n2:
                                pads2 = np.broadcast_to(re_2[np.newaxis, np.newaxis, :, :],
                                                        re_1.shape[:2] + re_2.shape[:2] + re_2.shape[-1:]
                                                        ).reshape(-1, n2 - m)
                                cats.append(pads2)

                            new_combo = np.concatenate(
                                cats,
                                axis=1
                            )


                        # print(new_combo.shape,
                        #       (ch1.shape, inds_1.shape, ch2.shape, inds_2.shape))
                        # print(new_combo)
                        if uinds is not None:
                            real_inds = np.ravel_multi_index(
                                [t.flatten() for t in np.meshgrid(*uinds, indexing='ij')],
                                (len(ch1) * len(inds_1), len(ch2) * len(inds_2))
                            )
                        else:
                            real_inds = np.arange(len(new_combo))
                        unwrap_inds = np.unravel_index(
                            real_inds,
                            (len(ch1), len(inds_1), len(ch2), len(inds_2))
                        )

                        new_combo = new_combo[real_inds,]
                        sorts = self.change_sort(new_combo)
                        nzs = np.count_nonzero(new_combo, axis=1)

                        # we split the array up by the numbers of zeros to make further ops faster since those are dropped
                        (numzs, ind_blocks), _ = nput.group_by(np.arange(len(nzs)), nzs)
                        for n, ix in zip(numzs, ind_blocks):

                            subsorts = sorts[ix]
                            new_changes = nput.vector_take(new_combo[ix], subsorts[:, :n], shared=1)
                            full_ch1 = ch1[unwrap_inds[0][ix]]
                            full_ch2 = ch2[unwrap_inds[2][ix]]
                            contact_1 = inds_1[unwrap_inds[1][ix]]
                            contact_2 = inds_2[unwrap_inds[3][ix]]
                            uncontact_1 = rems_1[unwrap_inds[1][ix]]
                            uncontact_2 = rems_2[unwrap_inds[3][ix]]

                            if _DEBUG_PRINT:
                                combo_terms = []
                                if contact_1.shape[1] > 0 and contact_2.shape[1] > 0:
                                    combo_terms.append(
                                        nput.vector_take(full_ch1, contact_1, shared=1)
                                        + nput.vector_take(full_ch2, contact_2, shared=1)
                                    )
                                elif contact_1.shape[1] > 0:
                                    combo_terms.append(
                                        nput.vector_take(full_ch1, contact_1, shared=1)
                                    )
                                elif contact_2.shape[1] > 0:
                                    combo_terms.append(
                                        nput.vector_take(full_ch2, contact_2, shared=1)
                                    )
                                if uncontact_1.shape[1] > 0:
                                    combo_terms.append(
                                        nput.vector_take(full_ch1, uncontact_1, shared=1)
                                    )
                                if uncontact_2.shape[1] > 0:
                                    combo_terms.append(
                                        nput.vector_take(full_ch2, uncontact_2, shared=1)
                                    )
                                # we do a quick check
                                test_changes = nput.vector_take(
                                    np.concatenate(combo_terms, axis=1),
                                    subsorts[:, :n],
                                    shared=1
                                )
                                if not np.allclose(new_changes, test_changes):
                                    print(";_;"*25)
                                    # print(new_combo)
                                    # print(new_combo[real_inds,])
                                    # print(ch1)
                                    # print(ch2)
                                    print(new_changes)
                                    print(test_changes)
                                    # print(combo_terms)
                                    raise ValueError(...)

                            self._fill_change_data(
                                full_changes, n, new_changes,
                                [full_ch1, full_ch2],
                                subsorts,
                                [contact_1, contact_2],
                                [uncontact_1, uncontact_2]
                            )

                    # new_changes = nput.vector_take(new_combo, sorts)

        changes = {}
        for n,ndat in full_changes.items():
            # we can now build the final changes by taking the _unique_ changes and splitting up the
            # stored values on a per-change basis
            # the final remaining difficulty is that we need a way to split the extra data
            # by the block it came from--as these objects are of arbitrary size and cannot be
            # packed into indvidual numpy arrays
            # to that point, we construct a block-column point array that we will use to find the actual elements
            # we are sampling
            block_lens = [len(x) for x in ndat['sorts']]
            # block_starts = [0] + list(np.cumsum(block_lens, dtype=int)[:-1])
            col_ptrs = np.concatenate([np.arange(x) for x in block_lens])
            row_ptrs = np.concatenate([np.full(x, n) for n,x in enumerate(block_lens)])

            concat_change = np.concatenate(ndat['changes'], axis=0)
            (uchanges, uinds), sorting = nput.group_by(
                np.arange(len(concat_change)),
                concat_change
            )

            for change, inds in zip(uchanges, uinds):
                change = tuple(change)
                changes[change] = [
                    (
                        ndat['init_changes'][0][r][c],
                        ndat['init_changes'][1][r][c],
                        ndat['contacts'][0][r][c],
                        ndat['contacts'][1][r][c],
                        ndat['rems'][0][r][c],
                        ndat['rems'][1][r][c],
                        ndat['sorts'][r][c]
                    )
                    for r, c in zip(row_ptrs[inds], col_ptrs[inds])
                ]

        return changes
    # @property
    # def changes(self):
    #     if self._changes is None:
    #         self._changes = self.build_changes()
    #     return self._changes

    def get_expressions(self):
        raise NotImplementedError("shouldn't need this here...")

    @classmethod
    def get_poly_product_terms(cls,
                               gen1, gen2, change_1, change_2,
                               target_inds, remainder_inds, reorg,
                               simplify=True
                               ):
        polys_1 = gen1.get_poly_terms(change_1)
        if isinstance(polys_1, PerturbationTheoryExpressionEvaluator): polys_1 = polys_1.expr
        if nput.is_numeric(polys_1):
            if polys_1 == 0:
                return 0
            else:
                raise ValueError("not sure how we got a number here...")
        try:
            polys_1.audit()
        except:
            raise ValueError(polys_1, gen1, change_1)
        if simplify:
            polys_1 = polys_1.combine()
            if nput.is_numeric(polys_1):
                if polys_1 == 0:
                    return 0
                else:
                    raise ValueError("not sure how we got a number here...")
        # we note that src_inds defines how we permute the change_2 inds so...I guess target_inds
        # tells us how we'd permute the inds of change_1?
        # if shift is not None:
        #     change_shift = change_1 + shift
        # else:
        #     change_shift = change_1
        # change_shift = [change_shift[i] for i in target_inds]

        shift_2 = [0] * max(len(change_1), len(change_2))
        for i,c in zip(*target_inds):
            # if i < len(change_1) and c < len(shift_2):
            shift_2[c] = change_1[i]

        if _DEBUG_PRINT:
            print("="*50)
            print(gen1, gen2)
            print("ch", change_1, change_2, target_inds, reorg, shift_2)

        polys_2 = gen2.get_poly_terms(change_2, shift=shift_2)
        if isinstance(polys_2, PerturbationTheoryExpressionEvaluator): polys_2 = polys_2.expr
        if nput.is_numeric(polys_2):
            if polys_2 == 0:
                return 0
            else:
                raise ValueError("not sure how we got a number here...")

        polys_2.audit()
        if simplify:
            polys_2 = polys_2.combine()
        if nput.is_numeric(polys_2):
            if polys_2 == 0:
                return 0
            else:
                raise ValueError("not sure how we got a number here...")
        if _DEBUG_PRINT:
            print("p1", polys_1.format_expr())
            print("p2", polys_2.format_expr())

        # for any changes done in polys_1 that _aren't_ touched by target_inds we need to introduce a
        # baseline shift to the mul_along
        baseline = list(change_1)
        for i in target_inds[0]: baseline[i] = 0
        base_polys = polys_1.mul_along(polys_2, target_inds, remainder=remainder_inds, baseline=baseline)
        try:
            base_polys.audit()
        except ValueError as e:
            raise ValueError("bad polynomial found in the evaluation of {}[{}][{}]x{}[{}][{}]".format(
                gen1, change_1, target_inds[0],
                gen2, change_2, target_inds[1]
            )) from e

        if _DEBUG_PRINT:
            print("3", base_polys.format_expr())
            print("_"*50)

        if reorg is not None:
            base_polys = base_polys.permute(reorg)
        base_polys = base_polys.strip()

        return base_polys

    def get_poly_terms(self, changes, allowed_paths=None, shift=None, simplify=True) -> SqrtChangePoly:
        t = tuple(changes)
        shift = tuple(shift) if shift is not None else None
        # check sorting?
        key = t if shift is None or all(s == 0 for s in shift) else (t, shift)
        base_changes = self.changes.get(key, None)
        if allowed_paths is not None:
            allowed_paths = set(allowed_paths)
        if not isinstance(base_changes, SqrtChangePoly):
            if nput.is_numeric(base_changes): return base_changes

            base_changes = self.changes.get(t, 0)
            if isinstance(base_changes, list):
                poly_changes = 0
                for change_1, change_2, contact_1, contact_2, rem_1, rem_2, reorg in base_changes:
                    if allowed_paths is not None and (change_1, change_2) not in allowed_paths:
                        continue
                    # try:
                    new_polys = self.get_poly_product_terms(
                        self.gen1, self.gen2, change_1, change_2,
                        [contact_1, contact_2], [rem_1, rem_2], reorg
                    )
                    # except:
                    #     raise Exception(change_1, change_2, src_inds, target_inds, reorg)
                    if nput.is_zero(new_polys): continue
                    try:
                        new_polys.audit()
                    except ValueError:
                        raise Exception("found bad polynomial in evaluating {}[{}[{}]]x{}[{}[{}]]".format(
                            self.gen1, change_1, contact_1,
                            self.gen2, change_2, contact_2
                            ))

                    # if isinstance(new_polys, SqrtChangePoly):
                    #     new_polys = new_polys.canonical_sort()

                    # -print("||||||>",
                    #       change_1, change_2, [contact_1, contact_2],
                    #       reorg, new_polys.poly_change
                    #       )

                    poly_changes += new_polys

                if not nput.is_numeric(poly_changes):
                    try:
                        poly_changes.audit()
                    except ValueError:
                        raise Exception("found bad polynomial in evaluating {}[{}]".format(self, changes))
                self.changes[t] = poly_changes
                base_changes = poly_changes
            if not nput.is_numeric(base_changes):
                if shift is not None:
                    if not isinstance(base_changes, SqrtChangePoly): raise ValueError(...)
                    base_changes = base_changes.shift(shift)
                    self.changes[key] = base_changes

        if nput.is_numeric(base_changes):
            return base_changes
        elif isinstance(base_changes, SqrtChangePoly):
            if [c for c in base_changes.poly_change] != [c for c in changes]:
                raise ValueError("[[{}]] change mismatch between expected {} and obtained {}".format(
                    self,
                    changes,
                    base_changes.poly_change
                ))
            if shift is None:
                shift = [] #* len(changes)
            test_shift = base_changes.shift_start
            if all(s == 0 for s in test_shift): test_shift = []
            spec_shift = shift
            if all(s == 0 for s in spec_shift): spec_shift = []
            if list(test_shift) != list(spec_shift):
                raise ValueError("[[ {} ]] shift mismatch between expected {} and obtained {}".format(
                    self,
                    shift,
                    base_changes.shift_start
                ))
            new = base_changes
        elif isinstance(base_changes, PolynomialInterface):
            new = SqrtChangePoly(base_changes, changes, shift)
        else:
            raise ValueError("got bad value {}".format(base_changes))

        new.audit()

        return new

class FlippedPerturbationTheoryTerm(PerturbationTheoryTerm):
    """
    Represents a term that will be multipled by on the left rather than the right
    for evaluating things like Y[1]M[0]Y[1], essentially changing raising operations to lowering
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
        changes = {}
        for k,v in self.base.changes.items():
            k = tuple(-i for i in k)
            changes[k] = v
        return changes

    def get_poly_terms(self, changes, shift=None) -> 'SqrtChangePoly':
        # Operators are all Hermitian (overall) so we can just adjust the
        # energy shifts
        if shift is None:
            shift = [0] * len(changes)
        flip_changes = [-c for c in changes]

        shift_change = changes
        if len(shift_change) < len(shift):
            shift_change = list(shift_change) + [0] * (len(shift) - len(shift_change))
        flip_shift = [s+c for s,c in zip(shift, shift_change)]
        # raise Exception(flip_changes, shift)
        # if len(shift) != len(changes):
        #     raise ValueError(shift, changes)
        # lc = len(changes)
        # ls = len(shift)
        # if ls < lc:
        #     changes = list(changes) + [0]*(lc - ls)
        # elif ls < lc:
        #     raise ValueError(shift, changes)
        #     shift = list(shift) + [0]*(lc - ls)
        # shift = [s+c for s,c in zip(shift, changes)]
        reexpressed_poly = self.base.get_poly_terms(flip_changes, flip_shift)#.shift(shift_change)
        if isinstance(reexpressed_poly, SqrtChangePoly):
            reexpressed_poly = reexpressed_poly.poly_obj.shift([c for c in shift_change])
        if nput.is_numeric(reexpressed_poly): return reexpressed_poly
        return SqrtChangePoly(reexpressed_poly, changes, shift)


class PerturbationTheoryExpressionEvaluator:
    def __init__(self, op, expr:'SqrtChangePoly', change=None, logger=None):
        self.op = op
        self.expr = expr
        if change is None:
            change = expr.poly_change
        self.change = change # TODO: check that expr has the same change...
        self.num_fixed = len(change)
        self.logger = Logger.lookup(logger)
    def __repr__(self):
        return "{}(<{}>, {})".format(type(self).__name__, self.num_fixed, self.expr)

    @staticmethod
    def _extract_coeffs(coeffs, coeff_indices, perm):

        coeff_tensor = coeffs[coeff_indices[0]][coeff_indices[1]] # order then identity
        if nput.is_zero(coeff_tensor): return 0

        idx = tuple(perm[i] for i in coeff_indices[2:])
        if len(idx) > 0:
            return coeff_tensor[idx]
        else:
            return coeff_tensor # just a number

    @classmethod
    def _eval_poly(cls, state, fixed, inds,
                   poly, change, baseline_shift,
                   verbose, logger):
        log_level = Logger.LogLevel.Normal if verbose else Logger.LogLevel.Debug
        if isinstance(poly, ProductPTPolynomialSum):
            subvals = [
                cls._eval_poly(state, fixed, inds, p, change, baseline_shift, verbose, logger)
                for p in poly.polys
            ]
            return poly.prefactor * np.sum(subvals, axis=0)

        logging_poly = poly.constant_rescale()
        # p = poly
        logger.log_print("prefactor: {p:.3f}", p=logging_poly.prefactor, log_level=log_level)
        logger.log_print("coeffs: {c}",
                         c=logging_poly.coeffs,
                         preformatter=lambda **vars:dict(vars, c=[np.round(c, 3) for c in vars['c']]),
                         log_level=log_level
                         )

        if fixed > 0 and len(inds) > 0:
            substates1 = state[:, :fixed]
            substates2 = state[:, inds,]
            substates = np.concatenate([substates1, substates2], axis=-1)
        elif fixed > 0:
            substates = state[:, :fixed]
        else:
            substates = state[:, inds,]
        if baseline_shift is not None:
            substates = substates + baseline_shift[np.newaxis]

        if substates.shape[-1] == 0:
            return np.ones(substates.shape[0])
        else:
            # if len(change) != len(poly.coeffs): raise ValueError(change, poly)
            poly_factor = np.prod(
                [
                    np.dot(c[np.newaxis, :], np.power(s[np.newaxis, :], np.arange(len(c))[:, np.newaxis]))
                    for s, c in zip(substates.T, poly.coeffs)
                ],
                axis=0
            )
            if not nput.is_numeric(poly_factor): poly_factor = poly_factor[0] # we padded for the dot products
            shifts_sqrts = [
                np.prod(
                    n[:, np.newaxis] + np.arange(
                        delta + 1 if delta < 0 else 1,
                        1 if delta < 0 else delta + 1
                    )[np.newaxis, :],
                    axis=-1
                )
                for n, delta in zip(substates.T, change)
            ]
            sqrt_factor = np.sqrt(np.prod(shifts_sqrts, axis=0))

            ct = poly.prefactor * poly_factor * sqrt_factor
            return ct

    @classmethod
    def _compute_energy_weights(cls, energy_changes, freqs, full_inds, perms):

        ec_pad = max(len(e) for e in energy_changes)
        if len(full_inds) < ec_pad:
            full_inds = list(full_inds) + list(ProductPTPolynomial.fast_ind_remainder(ec_pad, full_inds))
        ec = np.array([list(e) + [0]*(ec_pad - len(e)) for e in energy_changes])
        full_inds = full_inds[:ec.shape[1]]
        freq_pulls = np.array([freqs[p[full_inds,],] for p in perms])
        return np.prod(np.dot(freq_pulls, ec.T), axis=-1)

    @classmethod
    def _eval_perm(cls, expr, change, baseline_shift,
                   subset, state, perms, coeffs, freqs, cind_sets, num_fixed,
                   zero_cutoff, eval_cache, verbose, logger, log_scaled):
        log_level = Logger.LogLevel.Normal if verbose else Logger.LogLevel.Debug
        log_scaling = 219475.6 if log_scaled else 1
        #TODO: use combinatorics to skip some duplicate evaluations over permutations

        full_set = tuple(range(num_fixed)) + subset

        free = subset
        fixed = num_fixed
        cinds_remapped = [
            tuple(
                # we need to avoid doing the exact same term twice, _but_ order matters
                # so we can't sort the tuples
                PTTensorCoeffProductSum._symmetrize(
                    ci[:2] + tuple(free[j - fixed] if j >= fixed and len(free) > j - fixed else j for j in ci[2:])
                )
                for ci in cinds
            )
            for cinds in cind_sets
        ]
        prefactors = np.array([
            [
                np.prod([cls._extract_coeffs(coeffs, ci, perm) for ci in cinds])
                    if eval_cache.get(cinds, 0) < np.math.factorial(len(cinds)) else # should cache this
                0
                for cinds in cinds_remapped
            ]
            for perm in perms
        ])

        contrib = np.zeros([len(state), len(perms)])
        for cinds in cinds_remapped:
            eval_cache[cinds] = eval_cache.get(cinds, 0) + 1
        good_pref = np.where(np.abs(prefactors) > zero_cutoff) # don't bother to include useless terms
        if len(good_pref) == 0 or len(good_pref[0]) == 0: return contrib
        prefac_groups = nput.group_by(good_pref[0], good_pref[1])[0] # group perms by corresponding cinds

        # with logger.block(tag='fs: {fs}', fs=full_set):
        for g,pi in zip(*prefac_groups):
            prefacs = prefactors[pi,][:, g]
            g_key = cind_sets[g]
            with logger.block(
                    tag="{k} ({p})",
                    k=PTTensorCoeffProductSum.format_key(g_key),
                    p=prefacs,
                    log_level=log_level
            ):
                subexpr = expr.terms[g_key]
                if isinstance(subexpr, PTEnergyChangeProductSum):
                    subcontrib = np.zeros([len(state), len(perms)])
                    for echanges, polys in subexpr.terms.items():
                        with logger.block(tag="{e} * {p}",
                                          e=PTEnergyChangeProductSum.format_key(echanges),
                                          p=polys.format_expr(),
                                          log_level=log_level):
                            energy_factors = cls._compute_energy_weights(echanges, freqs, full_set, perms[pi,])
                            # TODO: check size of energy factor
                            poly_factor = cls._eval_poly(state, num_fixed, subset, polys, change, baseline_shift, verbose, logger)
                            if not nput.is_zero(poly_factor):
                                scaled_contrib = poly_factor[:, np.newaxis] / energy_factors[np.newaxis, :]
                                logger.log_print("engs: {ef}", ef=energy_factors.squeeze(), log_level=log_level)
                                logger.log_print("poly: {pf}", pf=poly_factor.squeeze(), log_level=log_level)
                                logger.log_print("sc: {pf}",
                                                 pf=(prefacs[np.newaxis, :]*scaled_contrib).squeeze()*log_scaling,
                                                 log_level=log_level)

                        subcontrib += scaled_contrib
                    subcontrib *= subexpr.prefactor
                elif isinstance(subexpr, (ProductPTPolynomial, ProductPTPolynomialSum)):
                    subcontrib = cls._eval_poly(state, num_fixed, subset, subexpr, change, baseline_shift, verbose, logger)
                    if not nput.is_zero(subcontrib):
                        subcontrib = subcontrib[:, np.newaxis]
                else:
                    raise ValueError("how the hell did we end up with {}".format(subexpr))

                val = prefacs[np.newaxis, :] * subcontrib
                logger.log_print("contrib: {e}", e=val.squeeze()*log_scaling, log_level=log_level)
                contrib += val

            logger.log_print("{c}", c=contrib.squeeze() * log_scaling, log_level=log_level)
        return contrib

    @classmethod
    def evaluate_polynomial_expression(cls,
                                       state, coeffs, freqs,
                                       expr, change, baseline_shift,
                                       num_fixed,
                                       op=None,
                                       logger=None,
                                       perms=None, zero_cutoff=None, verbose=False, log_scaled=True
                                       ):
        state = np.asanyarray(state)
        smol = state.ndim == 1
        if smol: state = state[np.newaxis]
        ndim = state.shape[-1]
        if perms is None: perms = np.arange(ndim)
        smol_p = perms.ndim == 1
        if smol_p: perms = perms[np.newaxis]

        if baseline_shift is not None and all(b == 0 for b in baseline_shift): baseline_shift = None

        logger = Logger.lookup(logger)

        if nput.is_zero(expr):
            res = np.zeros([len(state), len(perms)])
        else:

            # max_order = None
            free_ind_groups = {}
            for coeff_indices, poly_terms in expr.terms.items():
                num_inds = len(expr.get_inds(coeff_indices))
                free_inds = num_inds - num_fixed
                if free_inds not in free_ind_groups:
                    free_ind_groups[free_inds] = []
                free_ind_groups[free_inds].append(coeff_indices)
            #     max_order = max(max_order, poly_terms.order) if max_order is not None else poly_terms.order
            #
            # state_polys = np.power(state, np.arange(max_order+1))

            if zero_cutoff is None:
                zero_cutoff = 1e-18

            log_level = Logger.LogLevel.Normal if verbose else Logger.LogLevel.Debug
            if op is None: op = expr
            with logger.block(tag="Evaluating {op}({ch})", op=op, ch=change, log_level=log_level):

                # TODO: numba-ify this part
                contrib = np.zeros([len(state), len(perms)])
                eval_cache = {}  # so we don't hit some coeff sets too many times
                for free_inds, cind_sets in free_ind_groups.items():
                    # now we iterate over every subset of inds (outside of the fixed ones)
                    # excluding replacement (since that corresponds to a different coeff/poly)
                    if free_inds == 0:
                        with logger.block(tag="()", log_level=log_level):
                            contrib += cls._eval_perm(
                                expr, change, baseline_shift,
                                (), state, perms, coeffs, freqs,
                                cind_sets, num_fixed,
                                zero_cutoff,
                                eval_cache, verbose, logger, log_scaled
                            )
                    else:
                        for subset in itertools.combinations(range(num_fixed, ndim), r=free_inds):
                            for subperm in itertools.permutations(subset):
                                with logger.block(tag="{s}", s=subperm, log_level=log_level):
                                    contrib += cls._eval_perm(
                                        expr, change, baseline_shift,
                                        subperm, state, perms, coeffs, freqs,
                                        cind_sets, num_fixed,
                                        zero_cutoff, eval_cache,
                                        verbose, logger, log_scaled
                                    )
                    # raise Exception(cind_sets)
            res = expr.prefactor * contrib

        if smol_p: res = res[:, 0]
        if smol: res = res[0]
        return res



    def evaluate(self, state, coeffs, freqs, perms=None, zero_cutoff=None, verbose=False, log_scaled=True):
        # we do state-by-state evaluation for now although we will at some
        # point need to do this in batches

        state = np.asanyarray(state)
        smol = state.ndim == 1
        if smol: state = state[np.newaxis]
        ndim = state.shape[-1]
        if perms is None: perms = np.arange(ndim)
        smol_p = perms.ndim == 1
        if smol_p: perms = perms[np.newaxis]


        expr = self.expr
        if nput.is_zero(expr):
            return np.zeros([len(state), len(perms)])
        else:
            if isinstance(expr, SqrtChangePoly):
                poly_obj = expr.poly_obj
                change = expr.poly_change
                shift_start = expr.shift_start
            else:
                poly_obj = expr
                change = None
                shift_start = None

            if not isinstance(poly_obj, PTTensorCoeffProductSum):
                poly_obj = PTTensorCoeffProductSum({():poly_obj})

            if change is None:
                first_echange_poly = next(iter(poly_obj.terms.values()))
                if not isinstance(first_echange_poly, PTEnergyChangeProductSum):
                    change = []
                else:
                    first_ekey = next(iter(first_echange_poly.terms.keys()))
                    change = np.sum([list(k) for k in first_ekey], axis=1)


            return self.evaluate_polynomial_expression(
                state, coeffs, freqs,
                poly_obj, change, shift_start,
                self.num_fixed,
                op=self.op,
                logger=self.logger,
                perms=perms, zero_cutoff=zero_cutoff,
                verbose=verbose, log_scaled=log_scaled
            )


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

class PerturbationTheoryEvaluator:

    def __init__(self, solver:AnalyticPerturbationTheorySolver, expansion, freqs=None):
        self.solver = solver
        self.expansions = expansion
        if freqs is None: freqs = 2*np.diag(self.expansions[0][0])
        self.freqs = freqs

    def get_energy_corrections(self, states, order=None, expansions=None, freqs=None, verbose=False):
        if expansions is None: expansions = self.expansions
        if freqs is None: freqs = self.freqs
        if order is None: order = len(self.expansions) - 1

        energy_evaluators = [
            self.solver.energy_correction(o)([]) for o in range(order + 1) if o % 2 == 0
        ]

        corrs = [
            evaluator.evaluate(states, expansions, freqs, verbose=verbose, log_scaled=True)
            for evaluator in energy_evaluators
        ]

        return corrs

    def get_overlap_corrections(self, states, order=None, expansions=None, freqs=None, verbose=False):
        if expansions is None: expansions = self.expansions
        if freqs is None: freqs = self.freqs
        if order is None: order = len(self.expansions) - 1

        overlap_evaluators = [
            self.solver.overlap_correction(o)([]) for o in range(2, order + 1)
        ]

        corrs = [
            evaluator.evaluate(states, expansions, freqs, verbose=verbose, log_scaled=False)
            for evaluator in overlap_evaluators
        ]

        return np.concatenate([np.ones((len(states), 1)), np.zeros((len(states), 1))] + corrs, axis=1)

    def get_diff_map(self, state_map):
        states = []
        change_perms = {}
        for initial, finals in state_map:
            initial = np.asanyarray(initial)
            states.append(initial)
            finals = np.asanyarray(finals)
            diffs = finals - initial[np.newaxis, :]
            sort_keys = -(10*np.abs(diffs) - np.sign(diffs))
            perms = np.argsort(sort_keys, axis=1)
            sort_diffs = nput.vector_take(diffs, perms, shared=1)
            for ch, perm_blocks in zip(*nput.group_by(np.argsort(perms, axis=1), sort_diffs)[0]):
                ch = tuple(ch[ch != 0])
                change_perms[ch] = change_perms.get(ch, [])
                change_perms[ch].append(perm_blocks)

        for ch, pblock in change_perms.items():
            pblock = np.unique(np.concatenate(pblock, axis=0), axis=0)
            change_perms[ch] = pblock

        return np.array(states), change_perms

    @staticmethod
    def get_finals(initial, change, perms):
        base_change = np.pad(change, [0, len(initial) - len(change)])
        changes = nput.vector_take(base_change, perms)
        return np.asanyarray(initial)[np.newaxis] + changes

    PTCorrections = collections.namedtuple("PTCorrections",
                                           ['initial_states', 'final_states', 'corrections']
                                           )
    def _reformat_corrections(self, states, order, corrs, change_map):
        # reshape to allow for better manipulation
        final_corrs = [
            [[] for _ in range(order + 1)]
            for _ in range(len(states))
        ]
        final_states = [
            []
            for _ in range(len(states))
        ]
        for change, corr_block in corrs.items():
            for n, s in enumerate(states):
                for o, c in enumerate(corr_block):
                    final_corrs[n][o].append(c)
                finals = self.get_finals(s, change, change_map[change])
                final_states[n].append(finals)

        for n, cb in enumerate(final_corrs):
            for o, c in enumerate(cb):
                cb[o] = np.concatenate(c, axis=1)

        for n, s in enumerate(final_states):
            final_states[n] = np.concatenate(s, axis=0)

        return self.PTCorrections(states, final_states, [np.concatenate(c, axis=0).T for c in final_corrs])

    def _build_corrections(self, corr_gen, states, expansions, order, terms, paths, change_map, freqs, verbose):
        corrs = {}
        for o in range(order + 1):
            generator = corr_gen(o, allowed_terms=terms, allowed_changes=paths)
            for change, perms in change_map.items():
                corr_block = corrs.get(change, [])
                corrs[change] = corr_block
                expr = generator(change)
                gen_corrs = expr.evaluate(
                    states,
                    expansions,
                    freqs,
                    perms=perms,
                    verbose=verbose,
                    log_scaled=False
                )
                corr_block.append(gen_corrs)
        return corrs

    def get_state_by_state_corrections(self, generator, states, order=None,
                                       terms=None, paths=None,
                                       expansions=None, freqs=None, verbose=False):
        if expansions is None: expansions = self.expansions
        if freqs is None: freqs = self.freqs
        if order is None: order = len(self.expansions) - 1

        all_gs = nput.is_numeric(states[0][0])
        if all_gs:
            states = [
                [[0] * len(states[0]), states]
            ]
        states, change_map = self.get_diff_map(states)

        corrs = self._build_corrections(generator, states, expansions, order, terms, paths, change_map, freqs, verbose)
        return self._reformat_corrections(states, order, corrs, change_map)

    def get_matrix_corrections(self, states, order=None, expansions=None, freqs=None, verbose=False):
        gen = lambda o, **kw: self.solver.hamiltonian_expansion[o]
        return self.get_state_by_state_corrections(gen, states,
                                                   order=order, expansions=expansions, freqs=freqs, verbose=verbose)
    def get_wavefunction_corrections(self, states, order=None, expansions=None, freqs=None, verbose=False):
        return self.get_state_by_state_corrections(self.solver.wavefunction_correction, states,
                                                   order=order, expansions=expansions, freqs=freqs, verbose=verbose)

    def _prep_operator_expansion(self, expansions, operator_expansion):
        if expansions is None: expansions = self.expansions
        # if order is None: order = len(operator_expansion) - 1

        exps = []
        for base, op in zip(expansions, operator_expansion):
            # we always make the operator be the 5th element
            exps.append(
                base + [0] * (5 - len(base)) + [op]
            )

        return exps
    def get_operator_corrections(self, operator_expansion, states, order=None, expansions=None, freqs=None,
                                 terms=None, verbose=False):
        if expansions is None: expansions = self.expansions
        if order is None: order = len(operator_expansion) - 1

        exps = self._prep_operator_expansion(expansions, operator_expansion)

        return self.get_state_by_state_corrections(self.solver.operator_correction, states,
                                                   order=order, terms=terms, expansions=exps, freqs=freqs, verbose=verbose)

    def evaluate_expressions(self, states, exprs, expansions=None, operator_expansions=None, verbose=False):
        order = len(exprs) - 1
        expansions = (
            self._prep_operator_expansion(expansions, operator_expansions)
                if operator_expansions is not None else
            expansions
        )
        return self.get_state_by_state_corrections(
            lambda i, **kw: exprs[i],
            states,
            expansions=expansions,
            order=order,
            verbose=verbose
        )