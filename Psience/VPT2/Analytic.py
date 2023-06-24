import numpy as np

from McUtils.Combinatorics import StirlingS1, Binomial

__all__ = [
    'PTPoly',
    'TensorCoeffPoly',
    'AnalyticPTOperator',
    'RaisingLoweringClasses'
]

class PTPoly:
    """
    The core data structure, a polynomial induced by some set of raising and lowering operations,
    with switching between 1D and nD cases depending on the shape of the input coefficients
    """
    def __init__(self,
                 coeffs,
                 fourier_coeffs=None,
                 prefactor=None,
                 shift=None,
                 fourier_shape=None,
                 poly_shape=None
                 ):
        self._scaling = prefactor
        self._shift = shift
        self._poly_coeffs = np.asanyarray(coeffs) if coeffs is not None else None
        self._poly_shape = poly_shape
        self._fourier_coeffs = np.asanyarray(fourier_coeffs) if fourier_coeffs is not None else None
        self._fourier_shape = fourier_shape
    def __repr__(self):
        return "{}({}, {})".format(type(self).__name__, self.coeffs, self.scaling)

    @property
    def shape(self):
        if self._poly_coeffs is not None:
            return self._poly_coeffs.shape
        else:
            return self._poly_shape
    @property
    def fourier_shape(self):
        if self._fourier_coeffs is not None:
            return self._fourier_coeffs.shape
        else:
            return self._fourier_shape


    @property
    def scaling(self):
        return 1 if self._scaling is None else self._scaling
    @scaling.setter
    def scaling(self, s):
        self._scaling = s

    @property
    def num_fourier(self):
        return 2*len(self._poly_coeffs) if self._fourier_coeffs is None else len(self._fourier_coeffs)
    @scaling.setter
    def scaling(self, s):
        self._scaling = s

    @classmethod
    def _unpadded_ifft(cls, coeffs, shape):
        if shape is None:
            if not all(c%2 == 0 for c in coeffs.shape):
                raise ValueError("can't invert fourier coeffs with shape {} without known target shape".format(
                    coeffs.shape
                ))
            return cls._unpadded_ifft(
                coeffs,
                [c//2 for c in coeffs.shape]
            )
        if (
                len(coeffs.shape) != len(shape)
                or any(s > c for c, s in zip(coeffs.shape, shape))
        ):
            raise ValueError("can't coerce fourier coeffs of shape {} into shape {}".format(
                coeffs.shape,
                shape
            ))
        poly_coeffs = cls.ifft(coeffs)
        for s in reversed(shape):
            poly_coeffs = np.moveaxis(poly_coeffs, -1, 0)
            poly_coeffs = poly_coeffs[:s]
        return poly_coeffs
    @property
    def coeffs(self):
        if self._poly_coeffs is None:
            # base_coeffs = self.rifft(self._fourier_coeffs)
            # int_coeffs = np.round(np.real(base_coeffs))
            # if np.linalg.norm(
            #     int_coeffs - base_coeffs
            # ) > 1e-10:
            #     raise ValueError('fourier coeffs not inverted cleanly')
            self._poly_coeffs = self._unpadded_ifft(self._fourier_coeffs, self._poly_shape)
        if self._shift is not None and (
                not isinstance(self._shift, (int, float, np.integer, np.floating))
                or self._shift != 0
        ):
            self._poly_coeffs = self._compute_shifted_coeffs(
                self._poly_coeffs,
                self._shift
            )
            self._fourier_coeffs = None
            self._shift = None
        return self._poly_coeffs
    @coeffs.setter
    def coeffs(self, cs):
        self._poly_coeffs = cs

    @classmethod
    def _padded_fft(cls, coeffs, shape):
        if shape is None:
            return cls._padded_fft(
                coeffs,
                [2*c for c in coeffs.shape]
            )
        elif coeffs.shape != shape:
            if (
                    len(coeffs.shape) != len(shape)
                or any(s < c for c,s in zip(coeffs.shape, shape))
            ):
                raise ValueError("can't coerce coeffs of shape {} into shape {}".format(
                    coeffs.shape,
                    shape
                ))
            coeffs = np.pad(
                coeffs,
                [
                    [0, s - c]
                    for c,s in zip(coeffs.shape, shape)
                ]
            )
        return cls.fft(coeffs)
    @property
    def fourier_coeffs(self):
        if self._shift is not None and (
                not isinstance(self._shift, (int, float, np.integer, np.floating))
                or self._shift != 0
        ):
            self._poly_coeffs = self._compute_shifted_coeffs(
                self.coeffs,
                self._shift
            )
            self._shift = None
            self._fourier_coeffs = None
        if self._fourier_coeffs is None:
            self._fourier_coeffs = self._padded_fft(self._poly_coeffs, self._fourier_shape)
        return self._fourier_coeffs
    @fourier_coeffs.setter
    def fourier_coeffs(self, cs):
        self._fourier_coeffs = cs

    def __mul__(self, other):
        if isinstance(other, PTPoly):
            if self.fourier_shape != other.fourier_shape: # need to cast to consistent shape
                if len(self.fourier_shape) != len(other.fourier_shape):
                    raise ValueError("can't force shape {} into shape {}".format(
                        self.shape,
                        other.shape
                    ))
                consistent_shape = [
                    max(s1, s2)
                    for s1, s2 in zip(
                        self.fourier_shape,
                        other.fourier_shape
                    )
                ]
                fcs = self._padded_fft(self.coeffs, consistent_shape)
                ocs = other._padded_fft(other.coeffs, consistent_shape)
            else:
                fcs = self.fourier_coeffs
                ocs = other.fourier_coeffs
            return PTPoly(
                None,
                fourier_coeffs=fcs*ocs,
                poly_shape=[
                    s+o-1
                    for s,o in zip(self.shape, other.shape)
                ],
                shift=None,
                prefactor=self.scaling * other.scaling
            )
        else:
            return PTPoly(
                self._poly_coeffs,
                fourier_coeffs=self._fourier_coeffs,
                shift=self.shift,
                prefactor=self.scaling * other
            )

    def __add__(self, other):
        if isinstance(other, PTPoly):
            if self.shape != other.shape:  # need to cast to consistent shape
                if len(self.fourier_shape) != len(other.fourier_shape):
                    raise ValueError("can't force shape {} into shape {}".format(
                        self.shape,
                        other.shape
                    ))
                consistent_shape = [
                    max(s1, s2)
                    for s1, s2 in zip(
                        self.fourier_shape,
                        other.fourier_shape
                    )
                ]
                fcs = np.pad(self.coeffs,
                             [[0, c-s] for c,s in zip(consistent_shape, self.shape)]
                             )
                ocs = np.pad(other.coeffs,
                             [[0, c-s] for c,s in zip(consistent_shape, other.shape)]
                             )
            else:
                fcs = self.coeffs
                ocs = other.coeffs
            return PTPoly(
                fcs + ocs,
                fourier_coeffs=None,
                shift=None,
                prefactor=self.scaling * other.scaling
            )
        else:
            if self.scaling != 1:
                other = other / self.scaling
            new = self.coeffs.copy().flatten()
            new[0] += other
            new = new.reshape(self.coeffs.shape)
            return PTPoly(
                new,
                fourier_coeffs=None,
                fourier_shape=self.fourier_shape,
                shift=self.shift,
                prefactor=self.scaling
            )

    def shift(self, shift):
        if not isinstance(shift, (int, np.integer, float, np.floating)):
            shift = np.asanyarray(shift)
        return PTPoly(
            self._poly_coeffs,
            fourier_coeffs=self._fourier_coeffs,
            poly_shape=self._poly_shape,
            fourier_shape=self._fourier_shape,
            shift=(0 if self._shift is None else self._shift) + shift,
            prefactor=self._scaling
        )
    @classmethod
    def _compute_shifted_coeffs(cls, poly_coeffs, shift):
        # if fourier_coeffs is None:
        #     raise ValueError("need fourier coeffs for shifted coeff calc")
        if poly_coeffs.ndim == 1:
            shift = [shift]
        factorial_terms = np.array([np.math.factorial(x) for x in range(poly_coeffs.shape[0])])
        for s in poly_coeffs.shape[1:]:
            factorial_terms = np.expand_dims(factorial_terms, -1) * np.reshape(
                np.array([np.math.factorial(x) for x in range(s)]),
                [1]*factorial_terms.ndim + [s]
            )

        shift_terms = np.power(shift[0], np.arange(poly_coeffs.shape[0]))
        for f,s in zip(shift[1:], poly_coeffs.shape[1:]):
            shift_terms = np.expand_dims(shift_terms, -1) * np.reshape(
                np.power(f, np.arange(s)),
                [1]*shift_terms.ndim + [s]
            )

        rev_fac = factorial_terms
        for i in range(factorial_terms.ndim):
            rev_fac = np.flip(rev_fac, axis=i)
            shift_terms = np.flip(shift_terms, axis=i)

        poly_terms = poly_coeffs * factorial_terms
        shift_terms = shift_terms / rev_fac

        padded_shape = [
            s+o-1
            for s, o in zip(
                poly_terms.shape,
                shift_terms.shape
            )
        ]

        fcs = cls._padded_fft(poly_terms, padded_shape)
        ocs = cls._padded_fft(shift_terms, padded_shape)

        new = cls.ifft(fcs * ocs)
        for s in reversed(poly_terms.shape):
            new = np.moveaxis(new, -1, 0)
            new = new[-s:]
        new = new / factorial_terms

        return new

    @classmethod
    def _cook_nd_args(cls, a, invreal=0):  # pulled from numpy
        shapeless = 1
        s = list(a.shape)
        axes = list(range(-len(s), 0))
        # if invreal and shapeless:
        #     s[-1] = (a.shape[axes[-1]] - 1) * 2
        return s, axes

    @classmethod
    def _execute_dft(cls, x, sign=-1):
        # This is slow compared to Numpy, but fast enough for convenience
        # especially as we expect to be able to
        N = len(x)
        n = np.arange(N)
        k = n.reshape((N, 1))
        e = np.exp(sign*2j * np.pi * k * n / N)

        return np.dot(e, x)
        # return [sum(eei*xi for eei,xi in zip(ee, x)) for ee in e]

    # @classmethod
    # def _execute_real_dft(cls, x):
    #     # need a simple real-only implementation of Cooley-Tukey DFT

    @classmethod
    def _base_fft(cls, a, n=None, inv_norm=None, axis=-1, sign=None):  # adapted from numpy.fft
        a = np.asanyarray(a)
        if n is None:
            n = a.shape[axis]
        # output = _raw_fft(a, n, axis, False, True, inv_norm)

        if axis < 0:
            axis = a.ndim + axis
        if n is None:
            n = a.shape[axis]

        if a.shape[axis] != n:
            s = list(a.shape)
            index = [slice(None)] * len(s)
            if s[axis] > n:
                index[axis] = slice(0, n)
                a = a[tuple(index)]
            else:
                index[axis] = slice(0, s[axis])
                s[axis] = n
                z = np.zeros(s, a.dtype.char)
                z[tuple(index)] = a
                a = z

        fct = 1 / inv_norm(n)
        if axis == a.ndim - 1:
            r = fct * cls._execute_dft(a, sign=sign)
        else:
            a = np.swapaxes(a, axis, -1)
            r = fct * cls._execute_dft(a, sign=sign)
            r = np.swapaxes(r, axis, -1)

        return r

    @classmethod
    def symbolic_fft(cls, a, n=None, axis=-1):
        return cls._base_fft(a, n=n, axis=axis, inv_norm=lambda n:1, sign=-1)
    @classmethod
    def symbolic_ifft(cls, a, n=None, axis=-1):
        return cls._base_fft(a, n=n, axis=axis, inv_norm=lambda n:n, sign=1)

    @classmethod
    def symbolic_fftn(cls, a):
        a = np.asanyarray(a)
        s, axes = cls._cook_nd_args(a)
        itl = list(range(len(axes)))
        itl.reverse()
        for ii in itl:
            a = cls.symbolic_fft(a, n=s[ii], axis=axes[ii])
        return a
    @classmethod
    def symbolic_ifftn(cls, a):
        a = np.asanyarray(a)
        s, axes = cls._cook_nd_args(a)
        itl = list(range(len(axes)))
        itl.reverse()
        for ii in itl:
            a = cls.symbolic_ifft(a, n=s[ii], axis=axes[ii])
        return a

    @classmethod
    def fft(cls, a, clip=True):
        a = np.asanyarray(a)
        if a.dtype == object:
            if a.ndim > 1:
                base = cls.symbolic_fftn(a)
                if clip:
                    base = np.apply_along_axis(
                        lambda x: x.clip() if isinstance(x, TensorCoeffPoly) else x,
                        -1,
                        base
                    )
            else:
                base = cls.symbolic_fft(a)
                if clip:
                    base = np.array([
                        x.clip() if isinstance(x, TensorCoeffPoly) else x
                        for x in base
                    ], dtype=base.dtype)
            return base
        else:
            if a.ndim > 1:
                return np.fft.fftn(a)
            else:
                return np.fft.fft(a)

    @classmethod
    def ifft(cls, a, clip=True):
        a = np.asanyarray(a)
        if a.dtype == object:
            if a.ndim > 1:
                base = cls.symbolic_ifftn(a)
                if clip:
                    base = np.apply_along_axis(
                        lambda x: x.real().clip() if isinstance(x, TensorCoeffPoly) else np.real(x),
                        -1,
                        base
                    )
            else:
                base = cls.symbolic_ifft(a)
                if clip:
                    base = np.array([
                        x.real().clip() if isinstance(x, TensorCoeffPoly) else np.real(x)
                        for x in base
                    ], dtype=base.dtype)
            return base
        else:
            if a.ndim > 1:
                base = np.fft.ifftn(a)
            else:
                base = np.fft.ifft(a)
            if clip:
                base = np.real(base)
            return base

class RaisingLoweringPolynomial(PTPoly):
    """
    The polynomial induced by _a_ raising operations and _b_ lowering ops
    """
    _stirlings = None
    _binomials = None
    def __init__(self, a, b, prefactor, shift=None):
        if b > a:
            shift = b-a + (0 if shift is None else shift)
            a, b = b, a
        prefactor = self.binom(a+b, a)*(1 if prefactor is None else prefactor)
        super().__init__(
            self.get_reduced_raising_lowering_coeffs(a, b),
            prefactor=prefactor,
            shift=shift
        )
    @classmethod
    def s1(cls, i, j):
        if cls._stirlings is None or cls._stirlings.shape[0] <= i or cls._stirlings.shape[0] <= j:
            cls._stirlings = StirlingS1(2**np.ceil(np.log2(max([64, i, j]))))
        return cls._stirlings[i, j]
    @classmethod
    def binom(cls, i, j):
        if cls._binomials is None or cls._binomials.shape[0] <= i or cls._binomials.shape[0] <= j:
            cls._binomials = Binomial(2 ** np.ceil(np.log2(max([64, i, j]))))
        return cls._binomials[i, j]
    @classmethod
    def get_reduced_raising_lowering_coeffs(cls, a, b):
        return np.array([
            sum(
                (cls.s1(b-j, w)*cls.binom(a, w)*cls.binom(b, w)*np.math.factorial(w)/2**w)
                for w in range(0, b-j+1)
            )
        for j in range(0, b+1)
        ])

# class PTProductPoly:
#     """
#     Provides a nicer way to handle multidimensional product polynomials
#     """
#     def __init__(self, poly_terms):
#         self.terms = poly_terms # an iterable of 1D polys
#     def

class TensorCoeffPoly:
    """
    A semi-symbolic representation of a polynomial of tensor
    coefficients
    """

    def __init__(self, terms:dict, prefactor=1):
        self.terms = terms
        self.prefactor = prefactor
    def expand(self):
        if self.prefactor == 1:
            return self
        else:
            return TensorCoeffPoly({k:self.prefactor*v for k,v in self.terms.items()}, prefactor=1)
    @classmethod
    def monomial(cls, idx, value=1):
        return cls({(idx,):value})
    def __repr__(self):
        return "{}({},{})".format(type(self).__name__, self.terms,self.prefactor)

    @classmethod
    def idx_hash(cls, idx_tuple):
        # hashes tuples of indices to give a fast check that the tuples
        # are the same independent of ordering
        return sum(hash(x) for x in idx_tuple)

    @classmethod
    def _canonical_idx(cls, idx):
        # we need a way to sort indices, which we do by grouping by key length,
        # doing a standard sort for each length, and reconcatenating
        s_groups = {}
        for i in idx:
            l = len(i)
            grp = s_groups.get(l, [])
            s_groups[l] = grp
            grp.append(i)
        t = tuple(
            i
            for k in sorted(s_groups.keys())
            for i in sorted(s_groups[k])
        )
        return t

    def __mul__(self, other):
        if isinstance(other, (int, float, np.integer, np.floating)):
            if other == 1:
                return self
            elif other == 0:
                return TensorCoeffPoly({})
        if isinstance(other, TensorCoeffPoly):
            new_terms = {}
            term_hashes = {}
            for k,v in other.terms.items():
                for k2,v2 in self.terms.items():
                    new_hash = self.idx_hash(k) + self.idx_hash(k2)
                    if new_hash in term_hashes:
                        t = term_hashes[new_hash]
                    else:
                        t = self._canonical_idx(k + k2)
                        term_hashes[new_hash] = t
                    new_terms[t] = new_terms.get(t, 0) + v * v2
                    if new_terms[t] == 0:
                        del new_terms[t]
            return TensorCoeffPoly(new_terms, self.prefactor*other.prefactor)
        else:
            return TensorCoeffPoly(self.terms, self.prefactor*other)
            # new_terms = {}
            # for k,v in self.terms.items():
            #     new_terms[k] = self.prefactor*other*v
            # return TensorCoeffPoly(new_terms)
    def __rmul__(self, other):
        return self * other

    def __add__(self, other):
        if isinstance(other, (int, float, np.integer, np.floating)):
            if other == 0:
                return self
        self = self.expand()
        new_terms = {}
        if isinstance(other, TensorCoeffPoly):
            other = other.expand()
            for s in other.terms.keys() & self.terms.keys():
                v = other.terms[s]
                v2 = self.terms[s]
                new = v + v2
                if new != 0:
                    new_terms[s] = new
            for o in other.terms.keys() - self.terms.keys():
                new_terms[o] = other.terms[o]
            for o in self.terms.keys() - other.terms.keys():
                new_terms[o] = self.terms[o]
        else:
            for k2, v2 in self.terms.items():
                new = v2 + other
                if new != 0:
                    new_terms[k2] = new
        if len(new_terms) == 0:
            return 0
        return TensorCoeffPoly(new_terms)
    def __radd__(self, other):
        return self + other
    def __truediv__(self, other):
        return self * (1/other)
    # def __rtruediv__(self, other):
    #     return other * (1/self)

    def real(self):
        return TensorCoeffPoly({i:np.real(v) for i,v in self.terms.items()}, self.prefactor)
    def clip(self, decimals=14):
        threshold = 10**(-decimals)
        new_terms = {
            i: np.round(v, decimals)
            for i, v in self.terms.items()
            if np.abs(self.prefactor*v) > threshold
        }
        if len(new_terms) == 0:
            return 0
        return TensorCoeffPoly(new_terms, self.prefactor)


class PossiblePathTree:
    """
    A tree representing the different paths to change quanta by _k_ given
    `block_sizes` encoding how the different numbers raising/lowering operations
    for the different operators
    """
    def __init__(self, k, block_sizes):
        self.delta = k
        self.blocks = block_sizes


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

class AnalyticPTOperator:
    """
    Provides a concrete representation for an operator
    that induces a polynomial for a given set of quanta changes
    """
    def __init__(self, terms, coeff_generator):
        self.spec = terms
        self.nterms = len(terms)
        self.coeffs = coeff_generator
    def get_poly(self, changes, shifts):
        t_change = sum(abs(x) for x in changes)
        if t_change > self.nterms or (t_change%2) != (self.nterms%2):
            return 0
        # enumerate all paths over the given number of terms
        # that can lead to these changes

        for c in RaisingLoweringClasses(self.nterms, changes):
            # enumerate all paths of length

            for (a,b),s in zip(c, shifts):
                p = RaisingLoweringPolynomial(a, b, s)

class OperatorProductTerm:
    """
    Provides a generic wrapper for an operator term
    given the operator identities (strings of p/q terms)
    and the breakpoint for switching from `Pi_n` to `Pi_m`
    """

    def __init__(self, pt_operators):
        self.operators = pt_operators

    def get_poly(self, changes):
        """

        :param changes:
        :type changes:
        :return:
        :rtype:
        """



