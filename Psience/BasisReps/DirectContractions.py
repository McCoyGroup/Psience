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

__all__ = ['HarmonicDirectContraction', 'QuadraticProductCoefficients']


class QuadraticProductCoefficients:
    """Exact sum_r weights[r] factors[r,i,j] factors[r,k,l]."""
    _quadratic_product = True
    def __init__(self, factors, weights):
        self.factors = np.ascontiguousarray(factors, dtype=float)
        self.weights = np.ascontiguousarray(weights, dtype=float)
        if self.factors.ndim != 3 or self.factors.shape[1] != self.factors.shape[2]:
            raise ValueError('quadratic factors must have shape (rank, modes, modes)')
        if len(self.weights) != len(self.factors):
            raise ValueError('one weight is required per factor')
        self.shape = (self.factors.shape[1],) * 4
        self.ndim = 4

    def asarray(self):
        return np.einsum('r,rij,rkl->ijkl', self.weights, self.factors, self.factors)


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
def _contract(bras, kets, coefficients, factors, weights, tables, permutations,
              permutation_counts, degree, phase, threshold, symmetric):
    size, nmodes = bras.shape
    result = np.zeros((size, coefficients.shape[1]))
    minimum = np.empty(degree, dtype=np.int64)
    slots = np.empty(degree, dtype=np.int64)
    ordered = np.empty(degree, dtype=np.int64)
    choices = np.empty(degree//2, dtype=np.int64)
    unique_modes = np.empty(degree, dtype=np.int64)
    masks = np.empty(degree, dtype=np.int64)
    factored = len(weights) > 0
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
                if factored:
                    coefficient = 0.
                    for r in range(len(weights)):
                        coefficient += (weights[r]*factors[r,ordered[0],ordered[1]]
                                        *factors[r,ordered[2],ordered[3]])
                    if coefficient == 0:
                        continue
                else:
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
                    if factored:
                        result[pair,0] += value*coefficient
                    else:
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
        if not 1 <= self.degree <= 6:
            raise ValueError('direct harmonic contractions support product degree 1 through 6')
        self.nmodes = operator.mode_n
        self.threshold = operator.zero_threshold
        self.phase = operator.funcs.i_phase
        self.component_shape = ()
        self.symmetric = len(set(self.word)) == 1
        self.tables = None
        self.table_size = 0
        if isinstance(operator.coeffs, QuadraticProductCoefficients):
            if operator.skipped_indices is not None or operator.coefficient_threshold is not None:
                raise ValueError('factorized coefficients require exact unfiltered support')
            if self.word != ('x','p','x','p') or operator.axes is not None:
                raise ValueError('quadratic products require ordered x,p,x,p and identity axes')
            self.factors, self.weights = operator.coeffs.factors, operator.coeffs.weights
            self.coefficients = np.zeros((1,1))
            self.symmetric = False
        else:
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
            self.factors, self.weights = np.zeros((0,1,1)), np.zeros(0)
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
        result = _contract(bras,kets,self.coefficients,self.factors,self.weights,self.tables,
                           self.permutations,self.permutation_counts,self.degree,self.phase,
                           self.threshold,self.symmetric)
        return result.reshape((len(bras),)+self.component_shape)
