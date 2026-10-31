"""Selection-rule-driven harmonic product contractions.

The oscillator polynomials come from Psience's native matrix generators.
This layer enumerates only index multisets with the parity and degree needed
to connect each requested bra/ket, and reduces directly into the result.
"""
import itertools
import math
import numpy as np
from numba import njit
from McUtils.Combinatorics import UniquePermutations

from .FactorizedContractions import FactorizedCoefficients, FactorizedHarmonicContraction

__all__ = ['HarmonicDirectContraction', 'QuadraticProductCoefficients', 'FactorizedCoefficients']


class QuadraticProductCoefficients(FactorizedCoefficients):
    """Exact sum_r weights[r] factors[r,i,j] factors[r,k,l]."""
    _quadratic_product = True
    def __init__(self, factors, weights):
        factors=np.asanyarray(factors)
        if factors.ndim!=3 or factors.shape[1]!=factors.shape[2]:
            raise ValueError('quadratic factors must have shape (rank, modes, modes)')
        super().__init__((factors,factors),((0,1),(2,3)),weights)
        # Preserve the public leading-Coriolis adapter's factor attribute.
        self.factors=self.factor_tensors[0]


@njit(cache=True)
def _canonical_coefficients(coefficients, nmodes, degree):
    result = np.zeros_like(coefficients)
    slots = np.empty(degree, dtype=np.int64)
    for index in range(len(coefficients)):
        code = index
        for j in range(degree-1, -1, -1):
            slots[j] = code % nmodes
            code //= nmodes
        slots.sort()
        target = 0
        for j in range(degree):
            target = target * nmodes + slots[j]
        for j in range(coefficients.shape[1]):
            result[target, j] += coefficients[index, j]
    return result


@njit(cache=True)
def _contract(bras, kets, coefficients, tables, permutations,
              permutation_counts, degree, phase, threshold, symmetric):
    size, nmodes = bras.shape
    result = np.zeros((size, coefficients.shape[1]))
    minimum = np.empty(degree, dtype=np.int64)
    slots = np.empty(degree, dtype=np.int64)
    ordered = np.empty(degree, dtype=np.int64)
    choices = np.empty(degree//2, dtype=np.int64)
    unique_modes = np.empty(degree, dtype=np.int64)
    masks = np.empty(degree, dtype=np.int64)
    for pair in range(size):
        used = 0
        possible = True
        for mode in range(nmodes):
            delta = abs(int(bras[pair, mode])-int(kets[pair, mode]))
            if used + delta > degree:
                possible = False
                break
            for j in range(delta):
                minimum[used] = mode
                used += 1
        if not possible or (degree-used) % 2:
            continue
        extra = (degree-used)//2
        for j in range(extra):
            choices[j] = 0
        finished = False
        while not finished:
            for j in range(used):
                slots[j] = minimum[j]
            for j in range(extra):
                slots[used+2*j] = slots[used+2*j+1] = choices[j]
            slots.sort()
            equality = 0
            for j in range(degree-1):
                if slots[j] == slots[j+1]:
                    equality |= 1 << j
            nperms = 1 if symmetric else permutation_counts[equality]
            for perm in range(nperms):
                code = 0
                for j in range(degree):
                    mode = slots[j] if symmetric else slots[permutations[equality, perm, j]]
                    ordered[j] = mode
                    code = code * nmodes + mode
                nonzero = False
                for c in range(coefficients.shape[1]):
                    nonzero |= coefficients[code,c] != 0
                if not nonzero:
                    continue
                count = 0
                for j in range(degree):
                    mode = ordered[j]
                    found = -1
                    for u in range(count):
                        if unique_modes[u] == mode:
                            found = u
                            break
                    if found < 0:
                        found = count
                        unique_modes[count] = mode
                        masks[count] = 0
                        count += 1
                    masks[found] |= 1 << j
                value = phase
                for u in range(count):
                    mode = unique_modes[u]
                    value *= tables[masks[u], bras[pair,mode], kets[pair,mode]]
                if abs(value) >= threshold:
                    for c in range(coefficients.shape[1]):
                        result[pair,c] += value*coefficients[code,c]
            if extra == 0:
                finished = True
            else:
                j = extra-1
                while j >= 0 and choices[j] == nmodes-1:
                    j -= 1
                if j < 0:
                    finished = True
                else:
                    start = choices[j]+1
                    for k in range(j, extra):
                        choices[k] = start
    return result


class HarmonicDirectContraction:
    """A per-call numerical plan; coefficient arrays are never cached across calls."""
    def __init__(self, operator):
        self.factor_plan=None
        if isinstance(operator.coeffs,FactorizedCoefficients):
            self.factor_plan=FactorizedHarmonicContraction(operator)
            self.component_shape=self.factor_plan.component_shape
            return
        from .HarmonicOscillator import HarmonicProductOperatorTermEvaluator
        if not isinstance(operator.funcs, HarmonicProductOperatorTermEvaluator):
            raise TypeError('direct contraction requires the native harmonic term evaluator')
        self.word = tuple(operator.funcs.terms)
        self.degree = len(self.word)
        self.nmodes = operator.mode_n
        if self.degree == 0:
            c = np.asarray(operator.coeffs, dtype=float)
            self.component_shape = c.shape
            self.coefficients = c.reshape(1,-1)
            return
        if self.degree > 6:
            self.factor_plan=FactorizedHarmonicContraction.from_dense_operator(operator)
            self.component_shape=self.factor_plan.component_shape
            return
        self.nmodes = operator.mode_n
        self.threshold = operator.zero_threshold
        self.phase = operator.funcs.i_phase
        self.component_shape = ()
        self.symmetric = len(set(self.word)) == 1
        self.tables = None
        self.table_size = 0
        c = np.asanyarray(operator.coeffs)
        if np.iscomplexobj(c):
            raise TypeError('direct contraction requires real coefficients')
        axes = operator.axes
        if axes is None:
            axes = (tuple(range(c.ndim)),)*2
        left, right = tuple(axes[0]), tuple(axes[1])
        if len(left) != self.degree or sorted(left) != list(range(self.degree)):
            raise ValueError('direct contractions must contract every operator axis')
        ordered_axes = [right[left.index(i)] % c.ndim for i in range(self.degree)]
        remaining = [i for i in range(c.ndim) if i not in ordered_axes]
        c = np.transpose(c, ordered_axes+remaining)
        if c.shape[:self.degree] != (self.nmodes,)*self.degree:
            raise ValueError('coefficient mode axes do not match the oscillator basis')
        self.component_shape = c.shape[self.degree:]
        c = np.ascontiguousarray(c.reshape((self.nmodes**self.degree,-1)), dtype=float)
        # Respect explicitly requested coefficient exclusions before grouping.
        if operator.skipped_indices is not None:
            c = c.copy()
            for pos in operator.skipped_indices:
                c[np.ravel_multi_index(pos,(self.nmodes,)*self.degree)] = 0
        self.coefficients = (_canonical_coefficients(c,self.nmodes,self.degree)
                             if self.symmetric else c)
        patterns = 1 << (self.degree-1)
        nperms = math.factorial(self.degree)
        self.permutations = np.zeros((patterns,nperms,self.degree),dtype=np.int64)
        self.permutation_counts = np.zeros(patterns,dtype=np.int64)
        for equality in range(patterns):
            representative = np.zeros(self.degree,dtype=int)
            for j in range(1,self.degree):
                representative[j] = representative[j-1] + (not (equality & (1 << (j-1))))
            permutations,_ = UniquePermutations(representative).permutations(return_indices=True)
            self.permutations[equality,:len(permutations)] = permutations
            self.permutation_counts[equality] = len(permutations)

    def evaluate(self, pairs):
        if self.factor_plan is not None:
            return self.factor_plan.evaluate(pairs)
        from .HarmonicOscillator import HarmonicOscillatorMatrixGenerator
        bras,kets = pairs.state_pairs
        bras,kets = np.ascontiguousarray(bras.T),np.ascontiguousarray(kets.T)
        if bras.shape[1] != self.nmodes:
            raise ValueError('bra/ket mode count does not match the contraction coefficients')
        if self.degree == 0:
            return (np.all(bras == kets, axis=1)[:,None] * self.coefficients).reshape(
                (len(bras),)+self.component_shape)
        if not len(bras):
            return np.empty((0,)+self.component_shape)
        extent = int(max(np.max(bras),np.max(kets)))+1
        if extent > self.table_size:
            tables = np.zeros((1 << self.degree,extent,extent))
            i,j = np.indices((extent,extent))
            for mask in range(1,1 << self.degree):
                word = tuple(self.word[k] for k in range(self.degree) if mask & (1 << k))
                gen = HarmonicOscillatorMatrixGenerator.load_cached(word)
                tables[mask] = gen.evaluate_state_terms((i.ravel(),j.ravel())).reshape(extent,extent)
            self.tables,self.table_size = tables,extent
        result = _contract(bras,kets,self.coefficients,self.tables,
                           self.permutations,self.permutation_counts,self.degree,self.phase,
                           self.threshold,self.symmetric)
        return result.reshape((len(bras),)+self.component_shape)
