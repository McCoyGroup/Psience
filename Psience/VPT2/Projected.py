"""Block-projected analytic perturbation theory.

Finite model-space contractions stay in this correction layer. Sums over
the complementary oscillator space use the ordinary analytic term algebra.
"""
import hashlib
import functools
import weakref
import itertools
import numpy as np
from McUtils.Scaffolding import Logger

from .Analytic import (
    PerturbationTheoryTerm, PerturbationTheoryTermSum,
    PerturbationTheoryTermProduct, ScaledPerturbationTheoryTerm,
    PTTensorCoeffProductSum, PTEnergyChangeProductSum, SqrtChangePoly,
    ShiftedEnergyBaseline, make_product_polynomial,
    PerturbationTheoryExpressionEvaluator, PerturbationOperator, _log_evaluation_time,
)
from ..BasisReps import BasisStateSpace, HarmonicOscillatorProductBasis
from .Corrections import BasicAPTCorrections

_projected_evaluators = weakref.WeakSet()


def clear_caches():
    ValidResolventTerm._cache.clear()
    ProjectedTerm._cache.clear()
    _BlockScalarTerm._cache.clear()
    _projection_reference_data.cache_clear()
    _bound_projection_row.cache_clear()
    for evaluator in list(_projected_evaluators):
        if hasattr(evaluator,'_projected_block_cache'):
            del evaluator._projected_block_cache
        if hasattr(evaluator,'_projection_runtime_inputs'):
            del evaluator._projection_runtime_inputs


def coefficient_fingerprint(expansions):
    digest = hashlib.sha256()
    for order in expansions:
        for coefficient in order:
            value = np.asarray(coefficient)
            digest.update(repr((value.shape,value.dtype.str)).encode())
            digest.update(value.tobytes() if value.dtype != object else repr(coefficient).encode())
    return digest.digest()


class ProjectionFactor(int):
    """A formal projection/reference role, containing no numerical bindings."""
    _state_projection_factor = True

    def __new__(cls, kind, slot=0):
        if kind not in {'P','Q','R','V'} or not isinstance(slot, int) or slot < 0:
            raise ValueError('projection factors require a kind and nonnegative formal slot')
        obj = int.__new__(cls, 2 + 4*slot + {'P':0,'Q':1,'R':2,'V':3}[kind])
        obj.kind, obj.slot = kind, slot
        return obj

    def __reduce__(self):
        return type(self), (self.kind,self.slot)

    def __repr__(self):
        return 'Projection[{}:{}]'.format(self.kind,self.slot)

    @property
    def label(self):
        return '{}[{}]'.format(self.kind,self.slot)


class _ProjectionEvaluationContext:
    """Numerical call data; never a term or a symbolic compilation input."""
    def __init__(self,operator,bindings):
        self.operator,self.projection_bindings = operator,bindings
    def __repr__(self):
        return repr(self.operator)


def projection_inputs(context):
    return (context.projection_bindings.get('projectors',())
            if isinstance(context,_ProjectionEvaluationContext) else ())


class _BlockScalarTerm(PerturbationTheoryTerm):
    """Formal finite-block coefficient, supplied through numerical tensors."""
    _cache = {}
    def __init__(self,order,slot,logger=None,representation='path'):
        super().__init__(logger=logger,polynomial_representation=representation)
        self.order,self.slot = int(order),int(slot)
    @classmethod
    def lookup(cls,order,slot,logger=None,representation='path'):
        key = order,slot,representation
        if key not in cls._cache:cls._cache[key]=cls(order,slot,logger,representation)
        return cls._cache[key]
    def __repr__(self):
        return 'BlockCoefficient[{},{}]'.format(self.order,self.slot)
    def get_changes(self):return {():None}
    def get_core_poly(self,changes,shift=None):
        if changes:return 0
        poly = make_product_polynomial([],steps=0,representation=self.polynomial_representation)
        # Coefficient identity 1 denotes a kinetic tensor in native algebra;
        # use separate auxiliary rows, each with scalar identity 0.
        return SqrtChangePoly(PTTensorCoeffProductSum({((self.order+self.slot,0),):poly}),[],[])


class ValidResolventTerm(PerturbationOperator):
    """Occupation-validity marker without a system-specific model space."""
    _cache = {}
    def get_poly_terms(self, changes, shift=None, **opts):
        return PerturbationTheoryTerm.get_poly_terms(self, changes, shift=shift, **opts)
    def get_core_poly(self, changes, shift=None):
        base = super().get_core_poly(changes, shift=shift)
        if np.isscalar(base):return base
        valid = ProjectionFactor('V')
        return base.mul_simple(PTEnergyChangeProductSum.monomial((valid,)+tuple(changes),1))


class ProjectedTerm(PerturbationTheoryTerm):
    """A generic projection operator. Its role is bound only by evaluation."""
    _cache = {}
    def __init__(self,term,factor):
        super().__init__(logger=term.logger,polynomial_representation=term.polynomial_representation)
        self.term,self.factor = term,factor
    @classmethod
    def lookup(cls,term,factor):
        key=term,factor
        if key not in cls._cache:cls._cache[key]=cls(term,factor)
        return cls._cache[key]
    def __repr__(self):return '{}<{}>'.format(self.factor.label,repr(self.term))
    def get_changes(self):
        return {change:None for change in self.term.changes if change or self.factor.kind=='P'}
    def get_core_poly(self,changes,shift=None):
        base=self.term.get_poly_terms(changes,shift=shift)
        if np.isscalar(base) and base==0:return 0
        return base.mul_simple(PTEnergyChangeProductSum.monomial((self.factor,)+tuple(changes),1))


class IdentityTerm(PerturbationTheoryTerm):
    def __init__(self, logger=None, polynomial_representation='path'):
        super().__init__(logger=logger, polynomial_representation=polynomial_representation)

    def __repr__(self):
        return 'I'

    def get_changes(self):
        return {(): None}

    def get_core_poly(self, changes, shift=None):
        return SqrtChangePoly(PTTensorCoeffProductSum({(): make_product_polynomial(
            [], steps=0, representation=self.polynomial_representation
        )}), [], [])


@functools.lru_cache(maxsize=512)
def _projection_reference_data(binding, state):
    states,frequencies,reference,trivial = binding
    differences = np.asarray(states, dtype=int).reshape((-1, len(state))) - np.asarray(state)
    counts = np.count_nonzero(differences, axis=1)
    gap = (0. if reference is None else
           float(np.dot(np.asarray(reference)-state, frequencies)))
    return differences, counts, gap


def projected_path_factors(energy_key, state, permutations, perm_freqs, baseline_shift=None,
                           *, projection_bindings=()):
    """Return a projector mask and allowed denominators before any division."""
    if any(getattr(row[0],'_conditional_endpoint',False) for row in energy_key):
        from .Conditional import conditional_path_factors
        return conditional_path_factors(energy_key,state,permutations,perm_freqs,baseline_shift,
            projection_bindings=projection_bindings)
    state = np.asarray(state)[0]
    permutations = np.asarray(permutations, dtype=int)
    good = np.ones(len(permutations), dtype=bool)
    denominator = np.ones(len(permutations))
    for row in energy_key:
        factor, delta = row[0], np.asarray(row[1:], dtype=int)
        if not isinstance(factor, ProjectionFactor):
            denominator *= perm_freqs @ delta
            continue
        if factor.kind == 'V':
            good &= np.all(state[permutations] + delta >= 0, axis=1)
            continue
        binding = projection_bindings[factor.slot]
        if binding[3]:
            good &= np.all(state[permutations] + delta >= 0, axis=1)
            reference_gap = (0. if binding[2] is None else
                float(np.dot(np.asarray(binding[2])-state,binding[1])))
            denominator *= reference_gap - perm_freqs @ delta
            continue
        differences, counts, reference_gap = _projection_reference_data(binding, tuple(state))
        local = differences[:, permutations]
        matches = np.any(np.all(local == delta, axis=2)
                         & (np.count_nonzero(local, axis=2) == counts[:, None]), axis=0)
        good &= matches if factor.kind == 'P' else ~matches
        good &= np.all(state[permutations] + delta >= 0, axis=1)
        if factor.kind == 'R':
            denominator *= reference_gap - perm_freqs @ delta
    return good, denominator


@functools.lru_cache(maxsize=512)
def _bound_projection_row(energy_key, state, model_size, projection_bindings):
    factor_count = len(energy_key)
    deltas = np.asarray([row[1:] for row in energy_key], dtype=int)
    kinds = np.zeros(factor_count, dtype=int)
    gaps = np.zeros(factor_count)
    differences = np.zeros((factor_count, model_size, len(state)), dtype=int)
    counts = np.full((factor_count, model_size), -1, dtype=int)
    for i,row in enumerate(energy_key):
        factor = row[0]
        if isinstance(factor, ProjectionFactor):
            kinds[i] = {'R':1,'P':2,'Q':3,'V':4}[factor.kind]
            if factor.kind == 'V':continue
            binding = projection_bindings[factor.slot]
            if binding[3]:
                kinds[i] = 5
                gaps[i] = (0. if binding[2] is None else
                    float(np.dot(np.asarray(binding[2])-state,binding[1])))
                continue
            delta, number, gaps[i] = _projection_reference_data(binding,state)
            differences[i,:len(number)] = delta
            counts[i,:len(number)] = number
    return deltas,kinds,gaps,differences,counts


def projected_energy_rows(energy_keys, source_states, permutations, frequency_blocks, baseline_shift=None,
                          *, projection_bindings=()):
    """Contract a bounded indexed row batch, retaining all projector tests."""
    if any(getattr(row[0],'_conditional_endpoint',False) for key in energy_keys for row in key):
        from .Conditional import conditional_energy_rows
        return conditional_energy_rows(energy_keys,source_states,permutations,frequency_blocks,baseline_shift,
            projection_bindings=projection_bindings)
    model_size = max([1]+[len(projection_bindings[row[0].slot][0])
        for key in energy_keys for row in key if isinstance(row[0],ProjectionFactor)
        and row[0].kind!='V' and not projection_bindings[row[0].slot][3]])
    # Bound the projector-comparison scratch independently of the native
    # polynomial row chunk, whose bound does not include model-space size.
    factors = len(energy_keys[0])
    mode_count = len(np.asarray(source_states[0])[0])
    permutation_count, index_count = frequency_blocks.shape[1:]
    bytes_per_row = 8*factors*model_size*(mode_count+2*permutation_count*(index_count+1))
    row_limit = max(1,(16*1024**2)//max(1,bytes_per_row))
    if len(energy_keys) > row_limit:
        results = [projected_energy_rows(energy_keys[i:i+row_limit],
            source_states[i:i+row_limit],permutations[i:i+row_limit],frequency_blocks[i:i+row_limit],
            baseline_shift,projection_bindings=projection_bindings)
            for i in range(0,len(energy_keys),row_limit)]
        return tuple(np.concatenate([result[k] for result in results]) for k in range(2))
    states = np.stack([np.asarray(state)[0] for state in source_states])
    bound = [_bound_projection_row(key,tuple(state),model_size,projection_bindings)
             for key,state in zip(energy_keys,states)]
    deltas,kinds,gaps,differences,counts = [np.stack(items) for items in zip(*bound)]
    permutations = np.stack(permutations)
    linear = np.einsum('rpm,rfm->rpf',frequency_blocks,deltas,optimize=False)
    row_count,factor_count,model_count = differences.shape[:3]
    local = differences[
        np.arange(row_count)[:,None,None,None,None],
        np.arange(factor_count)[None,:,None,None,None],
        np.arange(model_count)[None,None,:,None,None],
        permutations[:,None,None,:,:]
    ]
    matches = (np.all(local == deltas[:,:,None,None,:],axis=-1)
               & (np.count_nonzero(local,axis=-1) == counts[:,:,:,None]))
    members = np.any(matches,axis=2)
    occupation = np.take_along_axis(states[:,None,:],permutations,axis=2)
    possible = np.all(occupation[:,None,:,:]+deltas[:,:,None,:] >= 0,axis=-1)
    factor_masks = np.where(kinds[:,:,None] == 0,True,
        possible & np.where(kinds[:,:,None] == 2,members,~members))
    good = np.all(factor_masks,axis=1)
    weights = np.where((kinds[:,None,:] == 1) | (kinds[:,None,:] == 5),gaps[:,None,:]-linear,
                       np.where(kinds[:,None,:] == 0,linear,1.))
    return good,np.prod(weights,axis=-1)


_runtime_roles = itertools.count(1)


class AnalyticBlockExpansion:
    """Generate canonical or intermediate-normalized block expansions.

    In the default canonical convention, the following Y/N/F definitions
    apply. With intermediate_normalization=True, Y is Omega, N[n>0]=0,
    and F is the generally non-Hermitian Bloch series P H Omega P. Its
    metric and canonical counterpart are computed lazily on request.

    Y_0=P; X_n=Q Y_n P; N_n=P Y_n P; F=Y^dagger H Y.
    Every X_n is built with a projected resolvent. No model-space energy
    difference is ever inverted. N_n and F_n are finite block contractions.

    Y[n][a] is a finite sum over model-space origins b. Its values are native
    analytic operator terms acting on |b>; the complementary oscillator
    states are summed by the existing polynomial/DAG evaluator, not enumerated
    as a Fock-state basis. The model-space overlap PY is chosen Hermitian.
    """
    def __init__(self, evaluator, states, expansions=None, freqs=None, zero_cutoff=None,
                 evaluation_backend=None, intermediate_normalization=False):
        self.intermediate_normalization = bool(intermediate_normalization)
        self._intermediate_expansion = None
        self._canonical_expansion = None
        self._metric = [np.eye(len(states))]
        self.evaluator = evaluator
        registry = getattr(evaluator,'_projection_runtime_inputs',None)
        if registry is None:
            registry = evaluator._projection_runtime_inputs = dict(projectors={},scalars={})
        self.runtime_inputs = registry
        self.states = tuple(tuple(int(n) for n in state) for state in states)
        self.freqs = np.asarray(evaluator.freqs if freqs is None else freqs)
        self._raw_expansions = expansions
        self.expansions = evaluator._prep_expansions(expansions)
        self.zero_cutoff = zero_cutoff
        self.evaluation_backend = evaluation_backend
        self.identity = IdentityTerm(logger=evaluator.solver.logger,
                                     polynomial_representation=evaluator.solver.polynomial_representation)
        self.Y = [[{state: self.identity} for state in self.states]]
        self.X = [[{} for state in self.states]]
        self.N = [np.eye(len(self.states))]
        self.F = [np.diag((np.asarray(self.states) + .5) @ self.freqs)]
        self._elements = {}
        self._coefficient_references = {}
        self._resolvents = {}
        self._inverse_normalizers = [np.eye(len(self.states))]
        self._intermediate_columns = [self.Y[0]]
        self._operator_pieces = {}

    def _scalar(self,value):
        role=next(_runtime_roles)
        self.runtime_inputs['scalars'][role]=float(value)
        return _BlockScalarTerm.lookup(len(self.expansions),role,
            self.evaluator.solver.logger,self.evaluator.solver.polynomial_representation)

    def _resolvent(self,term,state,origin):
        trivial=False
        if tuple(state)==tuple(origin):
            trivial=not any(
                tuple(other)!=tuple(state) and
                PerturbationTheoryTerm.sorted_changes((np.asarray(other)-origin)[np.asarray(other)!=origin]) in term.changes
                for other in self.states)
        role=next(_runtime_roles)
        self.runtime_inputs['projectors'][role]=(self.states,tuple(float(w) for w in self.freqs),tuple(state),trivial)
        return ProjectedTerm.lookup(term,ProjectionFactor('R',role))

    def _generic_term(self,term,scalar_order):
        projectors,scalars,memo = {},{},{}
        def walk(node):
            if node in memo:return memo[node]
            if isinstance(node,ProjectedTerm):
                role=node.factor.slot
                slot=projectors.setdefault(role,len(projectors))
                result=ProjectedTerm.lookup(walk(node.term),ProjectionFactor(node.factor.kind,slot))
            elif isinstance(node,_BlockScalarTerm):
                slot=scalars.setdefault(node.slot,len(scalars))
                result=_BlockScalarTerm.lookup(scalar_order,slot,node.logger,node.polynomial_representation)
            elif isinstance(node,PerturbationTheoryTermProduct):
                result=PerturbationTheoryTermProduct.lookup(walk(node.gen2),walk(node.gen1))
            elif isinstance(node,PerturbationTheoryTermSum):
                result=PerturbationTheoryTermSum.lookup(*(walk(t) for t in node.terms))
            elif isinstance(node,ScaledPerturbationTheoryTerm):
                result=ScaledPerturbationTheoryTerm.lookup(walk(node.base),node.prefactor)
            elif isinstance(node,ShiftedEnergyBaseline):
                result=ShiftedEnergyBaseline.lookup(walk(node.base))
            elif isinstance(node,ValidResolventTerm):
                result=ValidResolventTerm.lookup(walk(node.subterm))
            else:
                result=node
            memo[node]=result
            return result
        generic=walk(term)
        bindings=tuple(self.runtime_inputs['projectors'][role] for role in projectors)
        values=[self.runtime_inputs['scalars'][role] for role in scalars]
        return generic,dict(projectors=bindings,scalar_order=scalar_order,scalar_count=len(values)),values

    def _column(self,parts):
        grouped = {}
        for column, coefficient, left in parts:
            if coefficient == 0:
                continue
            # Skip exactly zero finite-block contractions at numerical assembly,
            # before requesting a generic equation. Never embed their values.
            if (isinstance(coefficient,_BlockScalarTerm)
                    and self.runtime_inputs['scalars'][coefficient.slot] == 0):
                continue
            for origin, term in column.items():
                if left is not None:
                    term = left if isinstance(term, IdentityTerm) else left * term
                if isinstance(coefficient,PerturbationTheoryTerm):
                    term = coefficient if isinstance(term,IdentityTerm) else coefficient * term
                elif coefficient != 1:
                    term = coefficient * term
                grouped.setdefault(origin, []).append(term)
        return {origin: terms[0] if len(terms) == 1 else PerturbationTheoryTermSum.lookup(*terms)
                for origin, terms in grouped.items()}

    def element(self, term, initial, final, expansions=None, logger=None, verbose=False):
        if isinstance(term, IdentityTerm):
            return float(tuple(initial) == tuple(final))
        expansions = self.expansions if expansions is None else expansions
        self._coefficient_references[id(expansions)] = expansions
        key = (term, tuple(initial), tuple(final), id(expansions))
        if key not in self._elements:
            logger = Logger.lookup(self.evaluator.solver.logger if logger is None else logger)
            level = Logger.LogLevel.Normal if verbose else Logger.LogLevel.Debug
            changes = self.evaluator.get_diff_map([[initial, [final]]], logger=logger, log_level=level)
            value = 0.
            generic,projection_bindings,scalar_values=self._generic_term(term,len(expansions))
            numerical_expansions=list(expansions)+[[value] for value in scalar_values]
            for change, state_perms in changes.items():
                with _log_evaluation_time(logger, "Building projected evaluator {s}[{c}]",
                                          s=term, c=change, log_level=level):
                    expression = PerturbationTheoryExpressionEvaluator(
                        generic, generic.get_poly_terms(change), change, logger=logger)
                with _log_evaluation_time(logger, "evaluating projected element...", log_level=level):
                    result = expression.evaluate(state_perms, numerical_expansions, self.freqs,
                        projection_bindings=projection_bindings,
                        only_degenerate_terms=False, zero_cutoff=self.zero_cutoff,
                        parallelizer=self.evaluator.parallelizer,
                        evaluation_backend=self.evaluation_backend,
                        log_scaled=False, verbose=verbose)
                value += float(result[0][0])
            self._elements[key] = value
        return self._elements[key]

    def inner(self, left, right, operator=None, expansions=None, logger=None, verbose=False):
        result = 0.
        for initial, right_term in right.items():
            if operator is not None:
                right_term = operator if isinstance(right_term, IdentityTerm) else operator * right_term
            for final, left_term in left.items():
                if isinstance(left_term, IdentityTerm):
                    expression = right_term
                elif isinstance(right_term, IdentityTerm):
                    expression = ShiftedEnergyBaseline.lookup(left_term)
                else:
                    expression = ShiftedEnergyBaseline.lookup(left_term) * right_term
                result += self.element(expression, initial, final, expansions, logger=logger, verbose=verbose)
        return result


    def _hamiltonian(self, order):
        if order >= len(self.expansions):
            return None
        if all(np.all(np.asarray(t) == 0) for t in self.expansions[order]):
            return None
        return self.evaluator.solver.hamiltonian_expansion[order]

    def extend(self, order, logger=None, verbose=False):
        if self.intermediate_normalization:
            return self._extend_intermediate(order,logger=logger,verbose=verbose)
        logger = Logger.lookup(self.evaluator.solver.logger if logger is None else logger)
        size = len(self.states)
        while len(self.Y) <= order:
            n = len(self.Y)
            with _log_evaluation_time(logger, "Building projected block ({size} states) at order {o}",
                                      size=size, o=n):
                with _log_evaluation_time(logger, "Contracting block normalization..."):
                    normalizer = np.zeros((size, size))
                    for a in range(size):
                        for b in range(size):
                            normalizer[a, b] = -.5 * sum(
                                self.inner(self.X[i][a], self.X[n-i][b], logger=logger, verbose=verbose) for i in range(1, n)
                            )
                    for i in range(1,n):
                        normalizer -= .5 * (self.N[i].T @ self.N[n-i])
                with _log_evaluation_time(logger, "Contracting block Hamiltonian..."):
                    energies = np.diag(self.F[0])
                    hamiltonian = (energies[:, None] - energies[None, :]) * normalizer
                    for r in range(1, n+1):
                        h = self._hamiltonian(r)
                        if h is not None:
                            bare = np.asarray([[self.element(h,initial,final,logger=logger,verbose=verbose) for initial in self.states]
                                               for final in self.states])
                            hamiltonian += bare @ self.N[n-r]
                            for a, state in enumerate(self.states):
                                for b in range(size):
                                    hamiltonian[a, b] += self.inner(
                                        {state: self.identity}, self.X[n-r][b], h, logger=logger, verbose=verbose)
                    for i in range(1, n):
                        hamiltonian -= self.N[i] @ self.F[n-i]
                with _log_evaluation_time(logger, "Building projected resolvents and wavefunctions..."):
                    projected, full = [], []
                    for a, state in enumerate(self.states):
                        parts = [(self.X[n-r][a], 1, self._hamiltonian(r))
                                 for r in range(1, n+1) if self._hamiltonian(r) is not None]
                        parts += [({origin:self.identity},self._scalar(self.N[n-r][b,a]),self._hamiltonian(r))
                                  for r in range(1,n+1) if self._hamiltonian(r) is not None
                                  for b,origin in enumerate(self.states)]
                        parts += [(self.X[i][b], self._scalar(-self.F[n-i][b, a]), None)
                                  for i in range(1, n) for b in range(size)]
                        source = self._column(parts)
                        column = {}
                        for origin, term in source.items():
                            key = (term, state, origin)
                            if key not in self._resolvents:
                                self._resolvents[key] = self._resolvent(term,state,origin)
                            column[origin] = self._resolvents[key]
                        projected.append(column)
                        normalization = [({origin: self.identity}, self._scalar(normalizer[b, a]), None)
                                         for b, origin in enumerate(self.states)]
                        full.append(self._column([(column, 1, None)] + normalization))
                self.N.append(normalizer)
                self.F.append(hamiltonian)
                self.X.append(projected)
                self.Y.append(full)
        return self

    def _extend_intermediate(self, order, logger=None, verbose=False):
        """Bloch recurrence: P Omega P=P, B=P H Omega P.

        F holds B in this explicitly selected representation; it is generally
        non-Hermitian. No overlap contraction is needed to build it.
        """
        logger=Logger.lookup(self.evaluator.solver.logger if logger is None else logger)
        size=len(self.states)
        while len(self.Y)<=order:
            n=len(self.Y)
            with _log_evaluation_time(logger,"Building intermediate-normalized block at order {o}",o=n):
                hamiltonian=np.zeros((size,size))
                for r in range(1,n+1):
                    h=self._hamiltonian(r)
                    if h is not None:
                        for a,state in enumerate(self.states):
                            for b in range(size):
                                hamiltonian[a,b]+=self.inner(
                                    {state:self.identity},self.Y[n-r][b],h,
                                    logger=logger,verbose=verbose)
                projected=[]
                for a,state in enumerate(self.states):
                    parts=[(self.Y[n-r][a],1,self._hamiltonian(r))
                           for r in range(1,n+1) if self._hamiltonian(r) is not None]
                    parts += [(self.X[i][b],self._scalar(-self.F[n-i][b,a]),None)
                              for i in range(1,n) for b in range(size)]
                    column={}
                    for origin,term in self._column(parts).items():
                        key=(term,state,origin)
                        if key not in self._resolvents:
                            self._resolvents[key]=self._resolvent(term,state,origin)
                        column[origin]=self._resolvents[key]
                    projected.append(column)
                self.F.append(hamiltonian)
                self.N.append(np.zeros((size,size)))
                self.X.append(projected)
                self.Y.append(projected)
        return self

    def _intermediate(self):
        if self.intermediate_normalization:
            return self
        if self._intermediate_expansion is None:
            self._intermediate_expansion=type(self)(self.evaluator,self.states,
                self._raw_expansions,self.freqs,self.zero_cutoff,self.evaluation_backend,
                intermediate_normalization=True)
        return self._intermediate_expansion

    def intermediate_columns(self, order, logger=None, verbose=False):
        """Build Omega directly, without first building normalized columns."""
        block=self._intermediate().extend(order,logger=logger,verbose=verbose)
        return block.Y[order]

    def intermediate_hamiltonian(self, order, logger=None, verbose=False):
        """Return the Bloch series B; use a general, not symmetric eigensolver.

        Individual entries differ from the orthonormal representation. Exact
        similarity preserves eigenvalues; truncated series can differ beyond
        their retained perturbation order.
        """
        return self._intermediate().extend(order,logger=logger,verbose=verbose).F[:order+1]

    def overlap_metric(self, order, logger=None, verbose=False):
        """Cache S=Omega^dagger Omega once per block and perturbation order."""
        block=self._intermediate().extend(order,logger=logger,verbose=verbose)
        size=len(block.states)
        while len(block._metric)<=order:
            n=len(block._metric)
            with _log_evaluation_time(Logger.lookup(block.evaluator.solver.logger if logger is None else logger),
                                      "Contracting deferred block metric at order {o}",o=n):
                metric=np.zeros((size,size))
                for a in range(size):
                    for b in range(size):
                        metric[a,b]=sum(block.inner(block.X[i][a],block.X[n-i][b],
                            logger=logger,verbose=verbose) for i in range(1,n))
                block._metric.append(metric)
        return block._metric[:order+1]

    def _metric_factors(self, order, logger=None, verbose=False):
        metric=self.overlap_metric(order,logger=logger,verbose=verbose)
        root=[np.eye(len(self.states))];inverse=[root[0]]
        for n in range(1,order+1):
            root.append(.5*(metric[n]-sum((root[i]@root[n-i] for i in range(1,n)),
                                          np.zeros_like(root[0]))))
            inverse.append(-sum((root[i]@inverse[n-i] for i in range(1,n+1)),
                                np.zeros_like(root[0])))
        return root,inverse

    def orthonormal_hamiltonian(self, order, logger=None, verbose=False):
        """Return S^(1/2) B S^(-1/2), consistently truncated at each order."""
        if not self.intermediate_normalization:
            return self.extend(order,logger=logger,verbose=verbose).F[:order+1]
        root,inverse=self._metric_factors(order,logger=logger,verbose=verbose)
        return [sum((root[i]@self.F[j]@inverse[n-i-j]
                     for i in range(n+1) for j in range(n+1-i)),
                    np.zeros_like(root[0])) for n in range(order+1)]

    def orthonormal_block(self, order, logger=None, verbose=False):
        """Reuse Omega and its cached metric for normalized observables.

        Y=Omega S^(-1/2). No second Q-space wavefunction recurrence is run.
        The returned block preserves the ordinary bra-H-ket contribution API.
        """
        if not self.intermediate_normalization:
            return self.extend(order,logger=logger,verbose=verbose)
        root,inverse=self._metric_factors(order,logger=logger,verbose=verbose)
        hamiltonians=self.orthonormal_hamiltonian(order,logger=logger,verbose=verbose)
        if self._canonical_expansion is None:
            self._canonical_expansion=type(self)(self.evaluator,self.states,
                self._raw_expansions,self.freqs,self.zero_cutoff,self.evaluation_backend)
            self._canonical_expansion._intermediate_expansion=self
        block=self._canonical_expansion
        while len(block.Y)<=order:
            n=len(block.Y);projected=[];full=[]
            for a,state in enumerate(self.states):
                column=self._column([(self.X[i][b],self._scalar(inverse[n-i][b,a]),None)
                    for i in range(1,n+1) for b in range(len(self.states))
                    if inverse[n-i][b,a]!=0])
                projected.append(column)
                full.append(self._column([(column,1,None)]+[
                    ({origin:self.identity},self._scalar(inverse[n][b,a]),None)
                    for b,origin in enumerate(self.states) if inverse[n][b,a]!=0]))
            block.N.append(inverse[n]);block.F.append(hamiltonians[n])
            block.X.append(projected);block.Y.append(full)
        return block

    def operator_piece(self, order, column, off_diagonal=False):
        key = (order, column, off_diagonal)
        if key not in self._operator_pieces:
            if order == 0:
                result = {} if off_diagonal else self.Y[0][column]
            elif off_diagonal:
                result = self._column([({state: self.identity}, self._scalar(self.N[order][b,column]), None)
                                       for b,state in enumerate(self.states) if b != column])
            else:
                result = self._column([(self.X[order][column],1,None),
                    ({self.states[column]:self.identity},self._scalar(self.N[order][column,column]),None)])
            self._operator_pieces[key] = result
        return self._operator_pieces[key]

    def operator_components(self, order, column, off_diagonal=None):
        # Keep finite P coefficients outside the free-space operator products.
        # Inlining a numerical N entry into a polynomial sum would specialize
        # and duplicate the same Q-space DAG for every reference state.
        result = []
        if off_diagonal is not True and self.X[order][column]:
            result.append((self.X[order][column],1.))
        result += [({state:self.identity},self.N[order][b,column])
                   for b,state in enumerate(self.states)
                   if self.N[order][b,column] != 0 and
                   (off_diagonal is None or (b != column) == off_diagonal)]
        return result

    def operator_correction(self, other, left, right, operators, order, expansions=None,
                            terms=None, branches=True, logger=None, verbose=False):
        if self.intermediate_normalization or other.intermediate_normalization:
            raise ValueError('normalized observables require orthonormal blocks or explicit metric factors')
        self.extend(order,logger=logger,verbose=verbose)
        other.extend(order,logger=logger,verbose=verbose)
        allowed = None if terms is None else {tuple(t) for t in terms}
        result = 0.
        for i in range(order+1):
            for k in range(order+1-i):
                j = order-i-k
                if j < len(operators) and (allowed is None or (k, j, i) in allowed):
                    if branches is True:
                        for bra,bc in self.operator_components(i,left):
                            for ket,kc in other.operator_components(k,right):
                                result += bc*kc*self.inner(bra,ket,operators[j],expansions,logger=logger,verbose=verbose)
                    else:
                        selected = {'main': True}
                        if branches is not False:
                            selected.update(branches if isinstance(branches,dict)
                                            else {key:True for key in branches})
                        for name,bra_cross,ket_cross in [('main',False,False),
                                ('left',True,False),('right',False,True),('both',True,True)]:
                            if selected.get(name,False):
                                for bra,bc in self.operator_components(i,left,bra_cross):
                                    for ket,kc in other.operator_components(k,right,ket_cross):
                                        result += bc*kc*self.inner(bra,ket,operators[j],expansions,logger=logger,verbose=verbose)
        return result


class ProjectedCorrectionEvaluator:
    """Resolve actual blocks and expose the ordinary analytic result layout."""
    def __init__(self, evaluator, blocks=None, pairs=None, expansions=None,
                 freqs=None, zero_cutoff=None, logger=None, verbose=False):
        self.evaluator = evaluator
        self.logger = Logger.lookup(evaluator.solver.logger if logger is None else logger)
        self.verbose = verbose
        self.expansions = expansions
        self.freqs = np.asarray(evaluator.freqs if freqs is None else freqs)
        self.zero_cutoff = zero_cutoff
        groups = [] if blocks is None else [list(group) for group in blocks]
        if blocks is None and pairs is not None:
            if isinstance(pairs, dict):
                for state, neighbors in pairs.items():
                    groups.append([state])
                    groups.extend(list(pair) for pair in neighbors)
            else:
                groups.extend(list(pair) for pair in pairs)
        components = []
        for group in groups:
            group = {tuple(int(n) for n in state) for state in group}
            touching = [part for part in components if part & group]
            for part in touching:
                group |= part
                components.remove(part)
            components.append(group)
        self.blocks = {state: tuple(sorted(group)) for group in components for state in group}
        _projected_evaluators.add(evaluator)
        context_key = (coefficient_fingerprint(evaluator.expansions if expansions is None else expansions),
                       tuple(self.freqs), zero_cutoff)
        cache = getattr(evaluator, '_projected_block_cache', None)
        if cache is None or cache[0] != context_key:
            cache = evaluator._projected_block_cache = (context_key, {})
        self.cache = cache[1]

    def block(self, state, order, intermediate=False):
        state = tuple(int(n) for n in state)
        group = self.blocks.get(state, (state,))
        key=(group,'intermediate') if intermediate else group
        if not intermediate and key not in self.cache and (group,'intermediate') in self.cache:
            self.cache[key]=self.cache[(group,'intermediate')].orthonormal_block(
                order,logger=self.logger,verbose=self.verbose)
        if key not in self.cache:
            self.cache[key] = AnalyticBlockExpansion(self.evaluator, group,
                self.expansions, self.freqs, self.zero_cutoff,
                intermediate_normalization=intermediate)
        block = self.cache[key].extend(order,logger=self.logger,verbose=self.verbose)
        return block, block.states.index(state)

    def diagonal(self, states, order, series, orders):
        corrections = []
        for n in orders:
            with self.logger.block(tag="Getting corrections at order {o}", o=n):
                with _log_evaluation_time(self.logger, "evaluating..."):
                    values = []
                    for state in states:
                        block, a = self.block(state, order)
                        values.append(getattr(block, series)[n][a,a])
                    corrections.append(np.asarray(values))
        return corrections

    def result(self, states, order, value, target_orders=None):
        if target_orders is not None:
            target_orders = frozenset(target_orders)
            if not target_orders or min(target_orders) < 0 or max(target_orders) > order:
                raise ValueError('target_orders must be nonempty orders between 0 and order')
        if np.ndim(states[0][0]) == 0:
            states = [[[0] * len(self.freqs), states]]
        grouped = {}
        for initial, finals in states:
            key = tuple(int(n) for n in initial)
            row = grouped.setdefault(key, {})
            for final in finals:
                row.setdefault(tuple(int(n) for n in final), None)
        basis = HarmonicOscillatorProductBasis(len(self.freqs))
        initial_space = BasisStateSpace(basis, list(grouped))
        final_spaces, corrections = [], []
        for initial, finals in grouped.items():
            final_spaces.append(BasisStateSpace(basis, list(finals)))
            row = []
            for n in range(order+1):
                if target_orders is not None and n not in target_orders:
                    row.append(None)
                    continue
                with self.logger.block(tag="Getting corrections at order {o}", o=n):
                    with _log_evaluation_time(self.logger, "evaluating..."):
                        row.append(np.asarray([value(initial, final, n) for final in finals]))
            corrections.append(row)
        return BasicAPTCorrections(initial_space, final_spaces, corrections)

    def wavefunctions(self, states, order, full=True, target_orders=None):
        targets = None if target_orders is None else set(target_orders)
        def value(initial, final, n):
            if targets is not None and n not in targets:
                return 0.
            block, a = self.block(initial, order,
                                  intermediate=bool(self.evaluator.solver.intermediate_normalization))
            column = block.Y[n][a] if full else block.X[n][a]
            if self.evaluator.solver.intermediate_normalization:
                column = block.intermediate_columns(n,logger=self.logger,verbose=self.verbose)[a]
                if not full and n == 0:
                    column = {}
            return -sum(block.element(term, origin, final,logger=self.logger,verbose=self.verbose) for origin, term in column.items())
        return self.result(states, order, value, target_orders=targets)

    @staticmethod
    def _operator_branches(branches):
        # Public correction rows are initial kets -> final bras. The legacy
        # "left" branch names the initial-state normalization, whereas the
        # matrix contraction's left factor is the final bra.
        if branches is True or branches is False:
            return branches
        if isinstance(branches,dict):
            return {('right' if key=='left' else 'left' if key=='right' else key):value
                    for key,value in branches.items()}
        return [('right' if key=='left' else 'left' if key=='right' else key) for key in branches]

    def hamiltonian(self, states, order, only_degenerate_terms=True, terms=None, branches=True,
                    target_orders=None, normalization='orthonormal'):
        if normalization not in {'orthonormal','intermediate'}:
            raise ValueError('unknown block normalization '+str(normalization))
        intermediate=normalization=='intermediate'
        if intermediate and (terms is not None or branches is not True):
            raise ValueError('Bloch entries are not normalized bra-H-ket contribution triples')
        def value(initial, final, n):
            block, a = self.block(initial, order, intermediate=intermediate)
            final = tuple(final)
            if final not in block.states:
                raise ValueError('Hamiltonian block endpoints must belong to the same model space')
            b = block.states.index(final)
            result = (block.F[n][b,a] if terms is None and branches is True else
                      block.operator_correction(block,b,a,self.evaluator.solver.hamiltonian_expansion,
                                                n,terms=terms,branches=self._operator_branches(branches),
                                                logger=self.logger,verbose=self.verbose))
            return result if only_degenerate_terms else -result
        return self.result(states, order, value, target_orders=target_orders)

    def operators(self, states, order, operators, expansions, terms=None, branches=True,
                  target_orders=None):
        def value(initial, final, n):
            right, ket = self.block(initial, order)
            left, bra = self.block(final, order)
            return -left.operator_correction(right, bra, ket, operators, n,
                    expansions=expansions, terms=terms,branches=self._operator_branches(branches),
                    logger=self.logger,verbose=self.verbose)
        return self.result(states, order, value, target_orders=target_orders)
