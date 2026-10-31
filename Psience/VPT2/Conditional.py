"""State-independent conditional operator algebra and ordered VPT transforms.

Conditions are attached to operator endpoints, not resolved model states.
All sums and products are native analytic perturbation-theory terms.
"""
import hashlib
import functools
import math
import pickle
import itertools
import json
import weakref
import numpy as np

from McUtils.Zachary.Symbolic.ElementaryFunctions import ElementaryFunction
from .Analytic import (PerturbationTheoryTerm, PerturbationTheoryTermSum,
    ScaledPerturbationTheoryTerm, PTEnergyChangeProductSum,
    ShiftedEnergyBaseline, PerturbationTheoryExpressionEvaluator, make_product_polynomial,
    OperatorExpansionTerm)
from .Projected import IdentityTerm

_stacks = weakref.WeakSet()


def clear_caches():
    ConditionalOperatorTerm._cache.clear()
    ConditionEndpoint._cache.clear()
    operator_sector.cache_clear()
    canonical_condition_key.cache_clear()
    IntegerLatticeCondition.may_retain_change.cache_clear()
    for stack in list(_stacks):
        cache = getattr(stack.solver,'_conditional_stacks',{})
        for key,value in list(cache.items()):
            if value is stack:
                del cache[key]
        stack._transformed_operators.clear()


def _local_changes(permutations, change, modes):
    permutations, change = np.asarray(permutations), np.asarray(change)
    return np.sum((permutations[...,None] == np.asarray(modes))
                  * change[...,None,:,None], axis=-2)


def _outside_unchanged(permutations, change, modes):
    return np.all((np.asarray(change)[...,None,:] == 0)
                  | np.isin(permutations,modes), axis=-1)


def _function_value(function, points):
    # McUtils elementary functions act on a single coordinate; multivariate
    # functions perform their own index selection.
    if isinstance(function,ElementaryFunction):
        if points.shape[-1] == 0:
            return function(np.zeros(points.shape[:-1]))
        return function(points[...,0 if function.idx is None else function.idx])
    return function(points)


class OperatorCondition:
    """A serializable relation on the two occupations of an operator step."""
    symmetric = False
    diagonal = None
    equivalence = False
    requires_frequencies = False

    def __setattr__(self,name,value):
        if 'key' in self.__dict__ or getattr(type(self),'key',None) is not None:
            raise AttributeError('operator conditions are immutable; construct a new condition')
        object.__setattr__(self,name,value)

    def __hash__(self):
        return hash(self.key)

    def __eq__(self, other):
        if self is other:
            return True
        return isinstance(other,OperatorCondition) and self.key == other.key

    def __and__(self, other):
        return BooleanCondition('and',self,other)

    def __or__(self, other):
        return BooleanCondition('or',self,other)

    def __invert__(self):
        return BooleanCondition('not',self)

    def symmetrized(self):
        return SymmetricCondition(self)

    def evaluate_shifted(self, state, permutations, start, end, frequencies):
        raise NotImplementedError

    def matches(self,initial,final,frequencies=None):
        initial,final = np.broadcast_arrays(np.asarray(initial),np.asarray(final))
        size = initial.shape[-1]
        if frequencies is None:
            if self.requires_frequencies:
                raise ValueError('this condition requires runtime frequencies')
            frequencies = np.zeros(size)
        permutations = np.broadcast_to(np.arange(size),initial.shape[:-1]+(1,size))
        values = self.evaluate_shifted(initial,permutations,np.zeros_like(initial),
            final-initial,np.broadcast_to(frequencies,initial.shape)[...,None,:])
        return values[...,0]

    __call__ = matches

    def may_retain_change(self,change):
        """False proves a change impossible; None leaves it unresolved."""
        return None


def _condition_modes(modes):
    modes = tuple(modes)
    result = tuple(int(mode) for mode in modes)
    if any(mode < 0 or mode != original for mode,original in zip(result,modes)):
        raise ValueError('condition modes must be nonnegative integer indices')
    if len(set(result)) != len(result):
        raise ValueError('condition mode indices must be distinct')
    return result


class SameStateCondition(OperatorCondition):
    symmetric = diagonal = True
    equivalence = True
    key = ('same-state',)

    def evaluate_shifted(self,state,permutations,start,end,frequencies):
        return np.broadcast_to(np.all(np.asarray(start) == end,axis=-1)[...,None],
                               np.shape(permutations)[:-1])

    def __repr__(self):
        return 'n=k'

    def may_retain_change(self,change):
        return not change


class BooleanCondition(OperatorCondition):
    def __init__(self, operation, *conditions):
        self.operation, self.conditions = operation, tuple(conditions)
        self.symmetric = all(condition.symmetric for condition in conditions)
        self.equivalence = operation == 'and' and all(c.equivalence for c in conditions)
        self.requires_frequencies = any(c.requires_frequencies for c in conditions)
        diagonal = [condition.diagonal for condition in conditions]
        if operation == 'and':
            self.diagonal = False if False in diagonal else True if all(v is True for v in diagonal) else None
        elif operation == 'or':
            self.diagonal = True if True in diagonal else False if all(v is False for v in diagonal) else None
        elif operation == 'not':
            self.diagonal = None if diagonal[0] is None else not diagonal[0]
        else:
            raise ValueError('unknown condition composition: '+operation)
        self.key = (operation,tuple(condition.key for condition in conditions))

    def evaluate_shifted(self,*args):
        values = [condition.evaluate_shifted(*args) for condition in self.conditions]
        if self.operation == 'not':
            return ~values[0]
        return (np.logical_and.reduce(values) if self.operation == 'and'
                else np.logical_or.reduce(values))

    def __repr__(self):
        return '{}({})'.format(self.operation,','.join(map(repr,self.conditions)))


class SymmetricCondition(OperatorCondition):
    symmetric = True

    def __init__(self,condition):
        self.condition, self.diagonal = condition, condition.diagonal
        self.requires_frequencies = condition.requires_frequencies
        self.key = ('symmetric',condition.key)

    def evaluate_shifted(self,state,permutations,start,end,frequencies):
        return (self.condition.evaluate_shifted(state,permutations,start,end,frequencies)
                | self.condition.evaluate_shifted(state,permutations,end,start,frequencies))

    def __repr__(self):
        return 'Symmetric[{}]'.format(repr(self.condition))


class PolynomialInvariantCondition(OperatorCondition):
    """Retain equal McUtils polynomial/function invariants at both endpoints.

    Variables of the McUtils functions refer to the supplied ``modes``.
    ``freeze_other_modes`` retains the spectator occupations unchanged.
    No target occupation, frequency, or Hamiltonian coefficient is captured.
    """
    symmetric = diagonal = True
    equivalence = True

    def __init__(self,functions,modes,freeze_other_modes=True,tolerance=0.):
        functions = tuple(functions) if isinstance(functions,(tuple,list)) else (functions,)
        snapshot = pickle.dumps(functions,protocol=5)
        self.functions = pickle.loads(snapshot)
        self.modes = _condition_modes(modes)
        self.freeze_other_modes = bool(freeze_other_modes)
        self.tolerance = 0. if tolerance == 0 else float(tolerance)
        if self.tolerance < 0 or np.isnan(self.tolerance):
            raise ValueError('a condition tolerance must be nonnegative')
        self.equivalence = self.tolerance == 0
        self.key = ('invariants',snapshot,self.modes,self.freeze_other_modes,self.tolerance)

    def evaluate_shifted(self,state,permutations,start,end,frequencies):
        source = np.asarray(state)[...,self.modes][...,None,:]
        left = source + _local_changes(permutations,start,self.modes)
        right = source + _local_changes(permutations,end,self.modes)
        result = np.ones(np.shape(permutations)[:-1],dtype=bool)
        for function in self.functions:
            result &= np.abs(_function_value(function,left)-_function_value(function,right)) <= self.tolerance
        if self.freeze_other_modes:
            result &= _outside_unchanged(permutations,np.asarray(end)-start,self.modes)
        return result

    def __repr__(self):
        return 'Invariant[{};modes={};freeze_other_modes={};tolerance={}]'.format(
            ','.join(map(repr,self.functions)),self.modes,self.freeze_other_modes,self.tolerance)

    def __reduce__(self):
        return type(self),(pickle.loads(self.key[1]),self.modes,self.freeze_other_modes,self.tolerance)


class PolynomialRelationshipCondition(OperatorCondition):
    """An arbitrary McUtils endpoint function compared with zero.

    The arguments are the selected source occupations followed by the target
    occupations. Use ``symmetrized`` for a directional relation in a unitary
    transformation; an invariant equality is symmetric by construction.
    """
    def __init__(self,function,modes,comparison='eq',tolerance=0.,freeze_other_modes=True):
        if comparison not in {'eq','ne','lt','le','gt','ge'}:
            raise ValueError('unknown polynomial comparison: '+comparison)
        snapshot = pickle.dumps(function,protocol=5)
        self.function = pickle.loads(snapshot)
        self.modes = _condition_modes(modes)
        self.comparison = comparison
        self.tolerance = 0. if tolerance == 0 else float(tolerance)
        if self.tolerance < 0 or np.isnan(self.tolerance):
            raise ValueError('a condition tolerance must be nonnegative')
        self.freeze_other_modes = bool(freeze_other_modes)
        self.key = ('relationship',snapshot,self.modes,comparison,self.tolerance,self.freeze_other_modes)

    def evaluate_shifted(self,state,permutations,start,end,frequencies):
        source = np.asarray(state)[...,self.modes][...,None,:]
        left = source + _local_changes(permutations,start,self.modes)
        right = source + _local_changes(permutations,end,self.modes)
        value = _function_value(self.function,np.concatenate([left,right],axis=-1))
        threshold = self.tolerance
        result = {'eq':lambda:abs(value)<=threshold,'ne':lambda:abs(value)>threshold,
                  'lt':lambda:value < -threshold,'le':lambda:value <= threshold,
                  'gt':lambda:value > threshold,'ge':lambda:value >= -threshold}[self.comparison]()
        if self.freeze_other_modes:
            result &= _outside_unchanged(permutations,np.asarray(end)-start,self.modes)
        return result

    def __reduce__(self):
        return type(self),(pickle.loads(self.key[1]),self.modes,self.comparison,
                           self.tolerance,self.freeze_other_modes)

    def __repr__(self):
        return 'Relationship[{};modes={};comparison={};tolerance={};freeze_other_modes={}]'.format(
            repr(self.function),self.modes,self.comparison,self.tolerance,self.freeze_other_modes)


class EnergyWindowCondition(OperatorCondition):
    """A harmonic energy-window relation, evaluated with runtime frequencies."""
    symmetric = diagonal = True
    requires_frequencies = True

    def __init__(self,tolerance,modes=None):
        self.tolerance = 0. if tolerance == 0 else float(tolerance)
        if self.tolerance < 0 or np.isnan(self.tolerance):
            raise ValueError('a condition tolerance must be nonnegative')
        self.modes = None if modes is None else _condition_modes(modes)
        self.key = ('harmonic-window',self.tolerance,self.modes)

    def evaluate_shifted(self,state,permutations,start,end,frequencies):
        change = np.asarray(end)-start
        result = abs(np.sum(frequencies*change[...,None,:],axis=-1)) <= self.tolerance
        if self.modes is not None:
            result &= _outside_unchanged(permutations,change,self.modes)
        return result

    def __repr__(self):
        return 'EnergyWindow[tolerance={};modes={}]'.format(self.tolerance,self.modes)


class IntegerLatticeCondition(OperatorCondition):
    """A displacement belongs to the integer lattice of resonance relations.

    Hermite reduction is done once symbolically. Runtime predicates use exact
    integer divisibility and reconstruction, not a floating null-space test.
    The family is valid at every occupation with nonnegative endpoints.
    """
    symmetric = diagonal = True
    equivalence = True

    def __init__(self,relations):
        relations = np.asarray(relations,dtype=int)
        if relations.ndim != 2:
            raise ValueError('relations must be a matrix of full mode changes')
        self.modes = tuple(int(mode) for mode in np.flatnonzero(np.any(relations != 0,axis=0)))
        if not self.modes:
            raise ValueError('a resonance lattice needs at least one nonzero relation')
        active = relations[:,self.modes]
        first = active[np.flatnonzero(np.any(active != 0,axis=1))[0]]
        primitive = first//np.gcd.reduce(abs(first))
        pivot = int(np.flatnonzero(primitive)[0])
        if primitive[pivot] < 0:
            primitive = -primitive
        coefficients = active[:,pivot]//primitive[pivot]
        if (np.all(active[:,pivot] % primitive[pivot] == 0)
                and np.all(coefficients[:,None]*primitive == active)):
            vector = primitive*int(np.gcd.reduce(abs(coefficients)))
            self.basis = tuple((int(value),) for value in vector)
            self.pivots,self.inverse_numerators = (pivot,),((1,),)
            self.denominator = int(vector[pivot])
            self.key = ('integer-lattice',self.modes,self.basis)
            return
        import sympy as sp
        from sympy.matrices.normalforms import hermite_normal_form
        basis = hermite_normal_form(sp.Matrix(relations[:,self.modes].T.tolist()))
        pivots = basis.T.rref()[1]
        inverse = basis[list(pivots),:].inv()
        denominator = math.lcm(*(int(value.q) for value in inverse))
        self.basis = tuple(tuple(int(value) for value in row) for row in basis.tolist())
        self.pivots = tuple(pivots)
        self.inverse_numerators = tuple(tuple(int(value*denominator) for value in row)
                                        for row in inverse.tolist())
        self.denominator = denominator
        self.key = ('integer-lattice',self.modes,self.basis)

    @classmethod
    def from_pairs(cls,pairs):
        if isinstance(pairs,dict):
            pairs = [pair for neighbors in pairs.values() for pair in neighbors]
        return cls([np.subtract(end,start) for start,end in pairs])

    def evaluate_shifted(self,state,permutations,start,end,frequencies):
        change = np.asarray(end,dtype=int)-start
        local = _local_changes(permutations,change,self.modes)
        numerators = local[...,self.pivots] @ np.asarray(self.inverse_numerators,dtype=int).T
        coordinates = numerators//self.denominator
        return (np.all(numerators % self.denominator == 0,axis=-1)
                & np.all(coordinates @ np.asarray(self.basis,dtype=int).T == local,axis=-1)
                & _outside_unchanged(permutations,change,self.modes))

    def __repr__(self):
        return 'Lattice[modes={},basis={}]'.format(self.modes,self.basis)

    @functools.lru_cache(maxsize=16384)
    def may_retain_change(self,change):
        if not change:
            return True
        size = len(self.modes)
        if len(change) > size:
            return False
        if math.perm(size,len(change)) > 512:
            return None
        permutations = list(itertools.permutations(range(size),len(change)))
        candidates = np.zeros((len(permutations),size),dtype=int)
        candidates[np.arange(len(permutations))[:,None],permutations] = change
        numerators = candidates[:,self.pivots] @ np.asarray(self.inverse_numerators,dtype=int).T
        coordinates = numerators//self.denominator
        return bool(np.any(np.all(numerators % self.denominator == 0,axis=-1)
                    & np.all(coordinates @ np.asarray(self.basis,dtype=int).T == candidates,axis=-1)))


class ConditionEndpoint(int):
    """One of two adjacent endpoint rows in a native analytic path key."""
    _state_projection_factor = True
    _conditional_endpoint = True
    _cache = {}

    def __new__(cls,condition,kind,retain,endpoint):
        if kind not in {'mask','resolvent'} or endpoint not in {0,1}:
            raise ValueError('invalid conditional path factor')
        cache_key = (condition,kind,bool(retain),endpoint)
        cached = cls._cache.get(cache_key)
        if cached is not None:
            return cached
        payload = (condition.key,kind,bool(retain),endpoint)
        number = 2 + int.from_bytes(hashlib.sha256(repr(payload).encode()).digest()[:7],'big')
        obj = int.__new__(cls,number)
        obj.condition,obj.kind,obj.retain,obj.endpoint = condition,kind,bool(retain),endpoint
        obj._frozen = True
        if len(cls._cache) >= 1024:
            cls._cache.pop(next(iter(cls._cache)))
        cls._cache[cache_key] = obj
        return obj

    def __setattr__(self,name,value):
        if getattr(self,'_frozen',False):
            raise AttributeError('conditional path factors are immutable')
        object.__setattr__(self,name,value)

    def __eq__(self,other):
        if self is other:
            return True
        return (isinstance(other,ConditionEndpoint) and self.condition == other.condition
                and (self.kind,self.retain,self.endpoint)==(other.kind,other.retain,other.endpoint))

    def __ne__(self,other):
        return not self == other

    __hash__ = int.__hash__

    def __reduce__(self):
        return type(self),(self.condition,self.kind,self.retain,self.endpoint)

    def __repr__(self):
        return '{}Endpoint[{},{}]'.format(self.kind,self.endpoint,self.condition)


class ConditionalOperatorTerm(PerturbationTheoryTerm):
    """Restrict a native operator or solve its conditional homological equation."""
    _cache = {}

    def __init__(self,term,condition,kind='mask',retain=True):
        if kind == 'resolvent' and retain:
            raise ValueError('a conditional VPT resolvent must exclude its retained relation')
        super().__init__(logger=term.logger,polynomial_representation=term.polynomial_representation)
        self.term,self.condition,self.kind,self.retain = term,condition,kind,bool(retain)

    @classmethod
    def lookup(cls,term,condition,kind='mask',retain=True):
        key = (term,condition,kind,bool(retain))
        if key not in cls._cache:
            cls._cache[key] = cls(*key)
        return cls._cache[key]

    def get_changes(self):
        return {change:None for change in self.term.changes
                if (change or not (self.kind == 'resolvent' or
                                   (not self.retain and self.condition.diagonal)))
                and not (self.retain and self.condition.may_retain_change(change) is False)}

    def get_core_poly(self,changes,shift=None):
        base = self.term.get_poly_terms(changes)
        if np.isscalar(base) and base == 0:
            return 0
        if not changes and self.kind == 'mask' and self.condition.diagonal:
            return base if self.retain else 0
        # A harmonic homological divisor is translation invariant. Keep it
        # in the native energy algebra and retain the separate endpoint mask.
        rows = ((ConditionEndpoint(self.condition,'mask',self.retain,0),)+(0,)*len(changes),
                (ConditionEndpoint(self.condition,'mask',self.retain,1),)+tuple(changes))
        if not self.retain and self.condition.may_retain_change(changes) is False:
            # This occupation-independent proof is stronger than a runtime
            # endpoint test. Keep the native divisor and omit the true mask.
            rows = ()
        if self.kind == 'resolvent':
            rows = ((0,)+tuple(changes),)+rows
        if not rows:
            return base
        result = base.mul_simple(PTEnergyChangeProductSum({rows:1}))
        return result.scale(-1) if self.kind == 'resolvent' else result

    def __repr__(self):
        return '{}[{};term@{}]'.format('R' if self.kind=='resolvent' else 'D' if self.retain else 'Q',
                                      self.condition,id(self.term))


class ZeroTerm(PerturbationTheoryTerm):
    def get_changes(self):
        return {}

    def get_core_poly(self,changes,shift=None):
        return 0

    def __repr__(self):
        return '0'


@functools.lru_cache(maxsize=8192)
def canonical_condition_key(original):
    """Cache immutable scalar-factor keys, never numerical numerators."""
    width = max((len(row)-1 for row in original),default=0)
    rows = [tuple(row)+(0,)*(width+1-len(row)) for row in original]
    legacy = any(getattr(row[0],'_state_projection_factor',False)
                 and not isinstance(row[0],ConditionEndpoint) for row in rows)
    units,masks,scale = [],{},1.
    index,possible = 0,True
    while index < len(rows):
        row = rows[index]
        factor = row[0]
        if isinstance(factor,ConditionEndpoint):
            if index+1 == len(rows):
                raise ValueError('unpaired conditional endpoint')
            other = rows[index+1]
            if (factor.endpoint != 0 or not isinstance(other[0],ConditionEndpoint)
                    or other[0].endpoint != 1 or factor.condition != other[0].condition
                    or (factor.kind,factor.retain)!=(other[0].kind,other[0].retain)):
                raise ValueError('conditional endpoint order was lost')
            start,end = row[1:],other[1:]
            if factor.kind == 'mask' and factor.condition.symmetric and end < start:
                start,end = end,start
            pair = ((factor,)+start,(other[0],)+end)
            if factor.kind == 'mask':
                key = (factor.condition,start,end)
                previous = masks.get(key)
                if previous is not None:
                    if previous != factor.retain:
                        possible = False
                        break
                    index += 2
                    continue
                masks[key] = factor.retain
            units.append(pair)
            index += 2
            continue
        if not legacy and factor in {0,1}:
            delta = tuple(int(value) for value in row[1:])
            divisor = math.gcd(*delta) if delta else 0
            if divisor:
                first = next(value for value in delta if value)
                signed = divisor if first > 0 else -divisor
                delta = tuple(value//signed for value in delta)
                scale /= signed
            row = (0,)+tuple(delta)
        units.append((row,))
        index += 1
    if not possible:
        return None,0.
    units.sort()
    return tuple(row for unit in units for row in unit),scale


def canonicalize_condition_paths(terms):
    """Combine commutative/idempotent masks and primitive harmonic divisors."""
    result = {}
    for original,polynomial in terms.items():
        key,scale = canonical_condition_key(tuple(original))
        if key is None:
            continue
        width = len(key[0])-1 if key else 0
        if scale != 1:
            polynomial = polynomial*scale if np.isscalar(polynomial) else polynomial.scale(scale)
        if np.isscalar(polynomial) and polynomial != 1:
            if polynomial == 0:
                continue
            polynomial = make_product_polynomial([np.ones(1)]*width,prefactor=polynomial,
                                                  steps=0,representation='eager')
        result[key] = result.get(key,0)+polynomial
    return result


def format_condition_path(key):
    """Render masks as multiplicative path functions, never energy divisors."""
    masks,divisors = [],[]
    def vector(values):
        return '{'+','.join(map(str,values))+'}'
    def energy(delta):
        return '('+(' + '.join('{}*w[{}]'.format(value,i)
                               for i,value in enumerate(delta) if value) or '0')+')'
    index = 0
    while index < len(key):
        row=key[index]
        if isinstance(row[0],ConditionEndpoint):
            factor,end = row[0],key[index+1]
            masks.append('PathMask[{},{},n,pi,{},{}]'.format(
                json.dumps(repr(factor.condition)),'True' if factor.retain else 'False',
                vector(row[1:]),vector(end[1:])))
            if factor.kind=='resolvent':
                divisors.append(energy(np.subtract(row[1:],end[1:])))
            index += 2
        else:
            divisors.append(energy(row[1:]))
            index += 1
    result='*'.join(masks) or '1'
    return result+('/('+'*'.join(divisors)+')' if divisors else '')


def format_condition_sum(expression):
    parts=['({})*({})'.format(poly if np.isscalar(poly) else poly.format_expr(),format_condition_path(key))
           for key,poly in expression.terms.items()]
    result=' + '.join(parts) or '0'
    return result if expression.prefactor==1 else '{}*({})'.format(expression.prefactor,result)


def filter_condition_paths(terms,test_keys,mode):
    """Filter divisor paths while retaining every occupation condition."""
    from McUtils.Zachary.Polynomials import TensorCoefficientPoly
    def normalize(row):
        row = tuple(row)
        delta = np.asarray(row[1:],dtype=int)
        divisor = int(np.gcd.reduce(abs(delta))) if len(delta) else 0
        if divisor:
            divisor = divisor if delta[np.flatnonzero(delta)[0]] > 0 else -divisor
            delta = delta//divisor
        return (0,)+tuple(delta)
    normalized = ([[normalize(row) for row in key] for key in test_keys] if mode=='match'
                  else [normalize(row) for row in test_keys])
    result = {}
    for key,poly in terms.items():
        divisors = tuple(row for row in key if not getattr(row[0],'_state_projection_factor',False))
        query = TensorCoefficientPoly({divisors:1},canonicalize=False)
        if query.filter(normalized,mode=mode).terms:
            result[key] = poly
    return result


def operator_sum(*terms):
    """Keep a shared native DAG and cancel equal formal summands exactly."""
    grouped = {}
    pending = [(term,1.) for term in terms]
    while pending:
        term,scale = pending.pop()
        if isinstance(term,ZeroTerm) or (np.isscalar(term) and term == 0):
            continue
        if isinstance(term,ScaledPerturbationTheoryTerm):
            pending.append((term.base,scale*term.prefactor))
        elif isinstance(term,PerturbationTheoryTermSum):
            pending.extend((child,scale) for child in term.terms)
        else:
            grouped[term] = grouped.get(term,0.)+scale
    terms = [term if scale == 1 else scale*term for term,scale in grouped.items() if scale != 0]
    if not terms:
        return ZeroTerm(logger=False)
    return terms[0] if len(terms)==1 else PerturbationTheoryTermSum.lookup(*terms)


def operator_product(left,right):
    if isinstance(left,ZeroTerm) or isinstance(right,ZeroTerm):
        return ZeroTerm(logger=False)
    if isinstance(left,IdentityTerm):
        return right
    if isinstance(right,IdentityTerm):
        return left
    return left*right


def commutator(left,right):
    if left is right or isinstance(left,IdentityTerm) or isinstance(right,IdentityTerm):
        return ZeroTerm(logger=False)
    return operator_sum(operator_product(left,right),-operator_product(right,left))


@functools.lru_cache(maxsize=32768)
def operator_sector(term,condition):
    """Prove a retained/complementary sector without inspecting any state."""
    if isinstance(term,IdentityTerm):
        return 'retained'
    if isinstance(term,OperatorExpansionTerm):
        changes = [change for change in term.changes
                   if not (len(change)==2 and isinstance(change[0],tuple))]
        if changes and all(not change for change in changes):
            return 'retained'
        if changes and all(condition.may_retain_change(change) is False for change in changes):
            return 'complement'
    if isinstance(term,ConditionalOperatorTerm) and term.condition == condition:
        return 'retained' if term.kind == 'mask' and term.retain else 'complement'
    if isinstance(term,ScaledPerturbationTheoryTerm):
        return operator_sector(term.base,condition)
    if isinstance(term,ShiftedEnergyBaseline) and condition.symmetric:
        return operator_sector(term.base,condition)
    if isinstance(term,PerturbationTheoryTermSum):
        sectors = {operator_sector(child,condition) for child in term.terms}
        return next(iter(sectors)) if len(sectors)==1 else None
    if condition.equivalence and hasattr(term,'gen1') and hasattr(term,'gen2'):
        left,right = operator_sector(term.gen1,condition),operator_sector(term.gen2,condition)
        if left == 'retained' and right == 'retained':
            return 'retained'
        if {left,right} == {'retained','complement'}:
            return 'complement'
    return None


class ConditionalTransformationStack:
    """Ordered all-order unitary transformations under a symbolic relation."""
    def __init__(self,solver,condition):
        if not isinstance(condition,OperatorCondition) or not condition.symmetric:
            raise ValueError('unitary resonance handling requires a symmetric operator condition')
        self.solver = solver
        self.condition = condition if condition.diagonal else (SameStateCondition() | condition)
        self.identity = IdentityTerm(logger=solver.logger,
                                     polynomial_representation=solver.polynomial_representation)
        self.zero = ZeroTerm(logger=solver.logger,
                             polynomial_representation=solver.polynomial_representation)
        self.order = len(solver.hamiltonian_expansion)-1
        self.hamiltonian = list(solver.hamiltonian_expansion)
        self.wavefunctions = [self.identity]+[self.zero]*self.order
        self.generators = [self.zero]
        self.stages = [tuple(self.hamiltonian)]
        self._transformed_operators = {}
        _stacks.add(self)

    def restrict(self,term,retain=True,resolvent=False):
        if isinstance(term,ZeroTerm):
            return term
        sector = operator_sector(term,self.condition)
        if sector is not None:
            if (sector == 'retained') != retain:
                return self.zero
            if not resolvent:
                return term
        if isinstance(term,PerturbationTheoryTermSum):
            return operator_sum(*(self.restrict(child,retain,resolvent) for child in term.terms))
        if (resolvent and isinstance(term,ConditionalOperatorTerm)
                and term.condition == self.condition and term.kind == 'mask' and not term.retain):
            term = term.term
        return ConditionalOperatorTerm.lookup(term,self.condition,
                         kind='resolvent' if resolvent else 'mask',retain=retain)

    def _conjugate(self,series,generator,stage):
        output = list(series)
        for source,term in enumerate(series):
            nested = term
            for power in range(1,(self.order-source)//stage+1):
                nested = commutator(nested,generator)
                target = source+power*stage
                output[target] = operator_sum(output[target],(1/math.factorial(power))*nested)
        return output

    def extend(self,order=None):
        order = self.order if order is None else order
        if order > self.order:
            raise ValueError('requested order exceeds the supplied Hamiltonian expansion')
        while len(self.generators) <= order:
            stage = len(self.generators)
            old = self.hamiltonian
            retained = self.restrict(old[stage])
            eliminated = self.restrict(old[stage],retain=False)
            generator = self.restrict(old[stage],retain=False,resolvent=True)
            # Combine the H0 and K_stage BCH contributions at every order:
            # ad_S^p(D K_stage)/p! + p ad_S^p(Q K_stage)/(p+1)!.
            # This is the exact homological identity, not a low-order rule.
            remaining = list(old)
            remaining[0] = remaining[stage] = self.zero
            updated = self._conjugate(remaining,generator,stage)
            updated[0],updated[stage] = old[0],retained
            d_nested,q_nested = retained,eliminated
            for power in range(1,(self.order-stage)//stage+1):
                d_nested = commutator(d_nested,generator)
                q_nested = commutator(q_nested,generator)
                target = (power+1)*stage
                updated[target] = operator_sum(updated[target],
                    (1/math.factorial(power))*d_nested,
                    (power/math.factorial(power+1))*q_nested)
            before = self.wavefunctions
            after = list(before)
            power_term = self.identity
            for power in range(1,self.order//stage+1):
                power_term = operator_product(power_term,generator)
                for source in range(self.order-power*stage+1):
                    target = source+power*stage
                    after[target] = operator_sum(after[target],
                        (1/math.factorial(power))*operator_product(before[source],power_term))
            self.wavefunctions,self.hamiltonian = after,updated
            self.generators.append(generator)
            self.stages.append(tuple(updated))
        return self

    def transform_operator(self,operators,order=None,terms=None,branches=True):
        order = self.order if order is None else order
        self.extend(order)
        allowed = None if terms is None else frozenset(tuple(term) for term in terms)
        branch_key = branches if isinstance(branches,bool) else tuple(sorted(
            branches.items() if isinstance(branches,dict) else ((name,True) for name in branches)))
        key = (tuple(operators),order,allowed,branch_key)
        if key not in self._transformed_operators:
            if allowed is None and branches is True:
                series = list(operators[:self.order+1])+[self.zero]*max(0,self.order+1-len(operators))
                for stage,generator in enumerate(self.generators[1:order+1],1):
                    series = self._conjugate(series,generator,stage)
                self._transformed_operators[key] = series[:order+1]
                return self._transformed_operators[key]
            corrections = []
            for n in range(order+1):
                parts = []
                for ket in range(n+1):
                    for bra in range(n-ket+1):
                        op = n-ket-bra
                        if op >= len(operators) or (allowed is not None and (ket,op,bra) not in allowed):
                            continue
                        if branches is True:
                            selections = [(False,False)]
                        else:
                            selected = {'main':True}
                            if branches is not False:
                                selected.update(dict(branch_key))
                            selections = [(bc,kc) for name,bc,kc in
                                [('main',False,False),('left',False,True),
                                 ('right',True,False),('both',True,True)] if selected.get(name,False)]
                        for bra_cross,ket_cross in selections:
                            left = self.wavefunctions[bra]
                            right = self.wavefunctions[ket]
                            if branches is not True:
                                diagonal = SameStateCondition()
                                def piece(term,n,cross):
                                    if n == 0:
                                        return self.zero if cross else self.identity
                                    condition = ((self.condition & ~diagonal) if cross else
                                                 (~self.condition | diagonal))
                                    return ConditionalOperatorTerm.lookup(term,condition)
                                left,right = piece(left,bra,bra_cross),piece(right,ket,ket_cross)
                            left = self.identity if isinstance(left,IdentityTerm) else ShiftedEnergyBaseline.lookup(left)
                            parts.append(operator_product(left,operator_product(operators[op],right)))
                corrections.append(operator_sum(*parts))
            self._transformed_operators[key] = corrections
        return self._transformed_operators[key]

    def intermediate_wavefunctions(self,order):
        self.extend(order)
        if not self.condition.equivalence:
            raise ValueError('intermediate normalization requires an equivalence relation')
        inside = [self.identity]+[self.restrict(term) for term in self.wavefunctions[1:order+1]]
        inverse = [self.identity]
        result = [self.identity]
        for n in range(1,order+1):
            inverse.append(-operator_sum(*(operator_product(inside[i],inverse[n-i])
                                             for i in range(1,n+1))))
            result.append(operator_sum(*(operator_product(self.wavefunctions[i],inverse[n-i])
                                         for i in range(n+1))))
        return result


def conditional_path_factors(energy_key,state,permutations,perm_freqs,baseline_shift=None):
    """Evaluate endpoint conditions and safe denominators before division."""
    state = np.asarray(state)[0]
    permutations = np.asarray(permutations,dtype=int)
    good = np.ones(len(permutations),dtype=bool)
    denominator = np.ones(len(permutations))
    width = permutations.shape[-1]
    baseline = np.zeros(width,dtype=int)
    if baseline_shift is not None:
        baseline[:min(width,len(baseline_shift))] = baseline_shift[:width]
    index = 0
    while index < len(energy_key):
        row = energy_key[index]
        if not isinstance(row[0],ConditionEndpoint):
            if getattr(row[0],'_state_projection_factor',False):
                from .Projected import projected_path_factors
                mask,weights = projected_path_factors((row,),[state],permutations,perm_freqs)
                good &= mask
                denominator *= weights
            else:
                denominator *= perm_freqs @ np.asarray(row[1:])
            index += 1
            continue
        factor = row[0]
        end_row = energy_key[index+1]
        end_factor = end_row[0]
        if (factor.endpoint != 0 or not isinstance(end_factor,ConditionEndpoint)
                or end_factor.endpoint != 1 or end_factor.condition != factor.condition
                or (factor.kind,factor.retain) != (end_factor.kind,end_factor.retain)):
            raise ValueError('conditional endpoint pair was not preserved by the analytic algebra')
        start,end = np.asarray(row[1:],dtype=int)+baseline,np.asarray(end_row[1:],dtype=int)+baseline
        possible = (np.all(state[permutations]+start >= 0,axis=-1)
                    & np.all(state[permutations]+end >= 0,axis=-1))
        matches = factor.condition.evaluate_shifted(state,permutations,start,end,perm_freqs)
        good &= possible & (matches if factor.retain else ~matches)
        if factor.kind == 'resolvent':
            denominator *= perm_freqs @ (start-end)
        index += 2
    return good,denominator


def conditional_energy_rows(energy_keys,source_states,permutations,frequency_blocks,
                            baseline_shift=None):
    values = [conditional_path_factors(key,state,perms,freqs,baseline_shift)
              for key,state,perms,freqs in zip(energy_keys,source_states,permutations,frequency_blocks)]
    return tuple(np.stack([value[i] for value in values]) for i in range(2))
