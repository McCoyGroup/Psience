"""Block-projected analytic perturbation theory.

Finite model-space contractions stay in this correction layer. Sums over
the complementary oscillator space use the ordinary analytic term algebra.
"""
import hashlib
import functools
import weakref
import numpy as np

from .Analytic import (
    PerturbationTheoryTerm, PerturbationTheoryTermSum,
    PTTensorCoeffProductSum, PTEnergyChangeProductSum, SqrtChangePoly,
    ShiftedEnergyBaseline, make_product_polynomial,
    PerturbationTheoryExpressionEvaluator, PerturbationOperator,
)
from ..BasisReps import BasisStateSpace, HarmonicOscillatorProductBasis
from .Corrections import BasicAPTCorrections

_projected_evaluators = weakref.WeakSet()


def clear_caches():
    ValidResolventTerm._cache.clear()
    _projection_reference_data.cache_clear()
    _bound_projection_row.cache_clear()
    for evaluator in list(_projected_evaluators):
        if hasattr(evaluator,'_projected_block_cache'):
            del evaluator._projected_block_cache


def coefficient_fingerprint(expansions):
    digest = hashlib.sha256()
    for order in expansions:
        for coefficient in order:
            value = np.asarray(coefficient)
            digest.update(repr((value.shape,value.dtype.str)).encode())
            digest.update(value.tobytes() if value.dtype != object else repr(coefficient).encode())
    return digest.digest()


class StateProjectionFactor(int):
    """Serializable diagonal projector/resolvent factor in an analytic path.

    An ordinary energy-change row records a relative denominator. This factor
    additionally retains its actual model space and fixed reference energy.
    Its occupation shift is transformed by the existing polynomial algebra.
    """
    _state_projection_factor = True

    def __new__(cls, kind, states, freqs, reference=None):
        states = tuple(tuple(int(n) for n in state) for state in states)
        freqs = tuple(float(w) for w in freqs)
        reference = None if reference is None else tuple(int(n) for n in reference)
        if kind not in {'P', 'Q', 'R', 'V'}:
            raise ValueError('unknown state projection factor: {}'.format(kind))
        if any(len(state) != len(freqs) for state in states):
            raise ValueError('model states and frequency vector have different dimensions')
        if kind == 'R' and reference not in states:
            raise ValueError('a projected resolvent reference must belong to its model space')
        payload = (kind, states, freqs, reference)
        number = 2 + int.from_bytes(hashlib.sha256(repr(payload).encode()).digest()[:7], 'big')
        obj = int.__new__(cls, number)
        obj.kind, obj.states, obj.freqs, obj.reference = payload
        return obj

    def __reduce__(self):
        return type(self), (self.kind, self.states, self.freqs, self.reference)

    @property
    def label(self):
        return self.kind if self.kind != 'R' else 'R' + str(self.reference)


class ValidResolventTerm(PerturbationOperator):
    """Reference-relative resolvent when the block projector is redundant."""
    _cache = {}

    def get_poly_terms(self, changes, shift=None, **opts):
        # A surrounding product can shift this column's reference occupation.
        # Its own zero-change component is excluded by get_changes; comparing
        # to the outer product's initial state would exclude additional paths.
        return PerturbationTheoryTerm.get_poly_terms(self, changes, shift=shift, **opts)

    def get_core_poly(self, changes, shift=None):
        base = super().get_core_poly(changes, shift=shift)
        if np.isscalar(base):
            return base
        valid = StateProjectionFactor('V', (), ())
        return base.mul_simple(PTEnergyChangeProductSum.monomial((valid,)+tuple(changes),1))


class ProjectedTerm(PerturbationTheoryTerm):
    """Apply P, Q, or a Q-projected fixed-reference resolvent to a term."""
    def __init__(self, term, states, freqs, reference=None, kind='R'):
        super().__init__(logger=term.logger,
                         polynomial_representation=term.polynomial_representation)
        self.term = term
        self.factor = StateProjectionFactor(kind, states, freqs, reference)

    def __repr__(self):
        return '{}<{}>'.format(self.factor.label, id(self.term))

    @classmethod
    def resolve(cls, term, states, freqs, reference, origin):
        # This is an order-independent support test, not a low-order shortcut.
        # If no nonzero change can reach another P state, Q excludes only the
        # reference, which the ordinary relative resolvent already removes.
        if tuple(reference) == tuple(origin):
            reachable = False
            for state in states:
                delta = np.subtract(state,origin)
                change = tuple(PerturbationTheoryTerm.sorted_changes(delta[delta != 0]))
                if change and change in term.changes:
                    reachable = True
                    break
            if not reachable:
                return -ValidResolventTerm.lookup(term)
        return cls(term,states,freqs,reference)

    def get_changes(self):
        return {change: None for change in self.term.changes
                if change or self.factor.kind == 'P'}

    def get_core_poly(self, changes, shift=None):
        base = self.term.get_poly_terms(changes, shift=shift)
        if np.isscalar(base) and base == 0:
            return 0
        delta = tuple(changes)
        factor = PTEnergyChangeProductSum.monomial((self.factor,) + delta, 1)
        return base.mul_simple(factor)


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
def _projection_reference_data(factor, state):
    differences = np.asarray(factor.states, dtype=int).reshape((-1, len(state))) - np.asarray(state)
    counts = np.count_nonzero(differences, axis=1)
    gap = (0. if factor.reference is None else
           float(np.dot(np.asarray(factor.reference)-state, factor.freqs)))
    return differences, counts, gap


def projected_path_factors(energy_key, state, permutations, perm_freqs):
    """Return a projector mask and allowed denominators before any division."""
    state = np.asarray(state)[0]
    permutations = np.asarray(permutations, dtype=int)
    good = np.ones(len(permutations), dtype=bool)
    denominator = np.ones(len(permutations))
    for row in energy_key:
        factor, delta = row[0], np.asarray(row[1:], dtype=int)
        if not isinstance(factor, StateProjectionFactor):
            denominator *= perm_freqs @ delta
            continue
        if factor.kind == 'V':
            good &= np.all(state[permutations] + delta >= 0, axis=1)
            continue
        differences, counts, reference_gap = _projection_reference_data(factor, tuple(state))
        local = differences[:, permutations]
        matches = np.any(np.all(local == delta, axis=2)
                         & (np.count_nonzero(local, axis=2) == counts[:, None]), axis=0)
        good &= matches if factor.kind == 'P' else ~matches
        good &= np.all(state[permutations] + delta >= 0, axis=1)
        if factor.kind == 'R':
            denominator *= reference_gap - perm_freqs @ delta
    return good, denominator


@functools.lru_cache(maxsize=512)
def _bound_projection_row(energy_key, state, model_size):
    factor_count = len(energy_key)
    deltas = np.asarray([row[1:] for row in energy_key], dtype=int)
    kinds = np.zeros(factor_count, dtype=int)
    gaps = np.zeros(factor_count)
    differences = np.zeros((factor_count, model_size, len(state)), dtype=int)
    counts = np.full((factor_count, model_size), -1, dtype=int)
    for i,row in enumerate(energy_key):
        factor = row[0]
        if isinstance(factor, StateProjectionFactor):
            kinds[i] = {'R':1,'P':2,'Q':3,'V':4}[factor.kind]
            delta, number, gaps[i] = _projection_reference_data(factor,state)
            differences[i,:len(number)] = delta
            counts[i,:len(number)] = number
    return deltas,kinds,gaps,differences,counts


def projected_energy_rows(energy_keys, source_states, permutations, frequency_blocks):
    """Contract a bounded indexed row batch, retaining all projector tests."""
    model_size = max(1,max(len(row[0].states) for key in energy_keys for row in key
                          if isinstance(row[0],StateProjectionFactor)))
    # Bound the projector-comparison scratch independently of the native
    # polynomial row chunk, whose bound does not include model-space size.
    factors = len(energy_keys[0])
    mode_count = len(np.asarray(source_states[0])[0])
    permutation_count, index_count = frequency_blocks.shape[1:]
    bytes_per_row = 8*factors*model_size*(mode_count+2*permutation_count*(index_count+1))
    row_limit = max(1,(16*1024**2)//max(1,bytes_per_row))
    if len(energy_keys) > row_limit:
        results = [projected_energy_rows(energy_keys[i:i+row_limit],
            source_states[i:i+row_limit],permutations[i:i+row_limit],frequency_blocks[i:i+row_limit])
            for i in range(0,len(energy_keys),row_limit)]
        return tuple(np.concatenate([result[k] for result in results]) for k in range(2))
    states = np.stack([np.asarray(state)[0] for state in source_states])
    bound = [_bound_projection_row(key,tuple(state),model_size)
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
    weights = np.where(kinds[:,None,:] == 1,gaps[:,None,:]-linear,
                       np.where(kinds[:,None,:] == 0,linear,1.))
    return good,np.prod(weights,axis=-1)


class AnalyticBlockExpansion:
    """Generate normalized block wavefunctions and Hamiltonians at any order.

    Y_0=P; X_n=Q Y_n P; N_n=P Y_n P; F=Y^dagger H Y.
    Every X_n is built with a projected resolvent. No model-space energy
    difference is ever inverted. N_n and F_n are finite block contractions.

    Y[n][a] is a finite sum over model-space origins b. Its values are native
    analytic operator terms acting on |b>; the complementary oscillator
    states are summed by the existing polynomial/DAG evaluator, not enumerated
    as a Fock-state basis. The model-space overlap PY is chosen Hermitian.
    """
    def __init__(self, evaluator, states, expansions=None, freqs=None, zero_cutoff=None,
                 evaluation_backend=None):
        self.evaluator = evaluator
        self.states = tuple(tuple(int(n) for n in state) for state in states)
        self.freqs = np.asarray(evaluator.freqs if freqs is None else freqs)
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

    @staticmethod
    def _column(parts):
        grouped = {}
        for column, coefficient, left in parts:
            if coefficient == 0:
                continue
            for origin, term in column.items():
                if left is not None:
                    term = left if isinstance(term, IdentityTerm) else left * term
                if coefficient != 1:
                    term = coefficient * term
                grouped.setdefault(origin, []).append(term)
        return {origin: terms[0] if len(terms) == 1 else PerturbationTheoryTermSum.lookup(*terms)
                for origin, terms in grouped.items()}

    def element(self, term, initial, final, expansions=None):
        if isinstance(term, IdentityTerm):
            return float(tuple(initial) == tuple(final))
        expansions = self.expansions if expansions is None else expansions
        self._coefficient_references[id(expansions)] = expansions
        key = (term, tuple(initial), tuple(final), id(expansions))
        if key not in self._elements:
            changes = self.evaluator.get_diff_map([[initial, [final]]])
            value = 0.
            for change, state_perms in changes.items():
                expression = PerturbationTheoryExpressionEvaluator(
                    term, term.get_poly_terms(change), change, logger=term.logger)
                result = expression.evaluate(state_perms, expansions, self.freqs,
                    only_degenerate_terms=False, zero_cutoff=self.zero_cutoff,
                    parallelizer=self.evaluator.parallelizer,
                    evaluation_backend=self.evaluation_backend,
                    log_scaled=False)
                value += float(result[0][0])
            self._elements[key] = value
        return self._elements[key]

    def inner(self, left, right, operator=None, expansions=None):
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
                result += self.element(expression, initial, final, expansions)
        return result


    def _hamiltonian(self, order):
        if order >= len(self.expansions):
            return None
        if all(np.all(np.asarray(t) == 0) for t in self.expansions[order]):
            return None
        return self.evaluator.solver.hamiltonian_expansion[order]

    def extend(self, order):
        size = len(self.states)
        while len(self.Y) <= order:
            n = len(self.Y)
            normalizer = np.zeros((size, size))
            for a in range(size):
                for b in range(size):
                    normalizer[a, b] = -.5 * sum(
                        self.inner(self.X[i][a], self.X[n-i][b]) for i in range(1, n)
                    )
            for i in range(1,n):
                normalizer -= .5 * (self.N[i].T @ self.N[n-i])
            energies = np.diag(self.F[0])
            hamiltonian = (energies[:, None] - energies[None, :]) * normalizer
            for r in range(1, n+1):
                h = self._hamiltonian(r)
                if h is not None:
                    bare = np.asarray([[self.element(h,initial,final) for initial in self.states]
                                       for final in self.states])
                    hamiltonian += bare @ self.N[n-r]
                    for a, state in enumerate(self.states):
                        for b in range(size):
                            hamiltonian[a, b] += self.inner(
                                {state: self.identity}, self.X[n-r][b], h)
            for i in range(1, n):
                hamiltonian -= self.N[i] @ self.F[n-i]
            projected, full = [], []
            for a, state in enumerate(self.states):
                parts = [(self.X[n-r][a], 1, self._hamiltonian(r))
                         for r in range(1, n+1) if self._hamiltonian(r) is not None]
                parts += [({origin:self.identity},self.N[n-r][b,a],self._hamiltonian(r))
                          for r in range(1,n+1) if self._hamiltonian(r) is not None
                          for b,origin in enumerate(self.states)]
                parts += [(self.X[i][b], -self.F[n-i][b, a], None)
                          for i in range(1, n) for b in range(size)]
                source = self._column(parts)
                column = {}
                for origin, term in source.items():
                    key = (term, state, origin)
                    if key not in self._resolvents:
                        self._resolvents[key] = ProjectedTerm.resolve(term,self.states,self.freqs,state,origin)
                    column[origin] = self._resolvents[key]
                projected.append(column)
                normalization = [({origin: self.identity}, normalizer[b, a], None)
                                 for b, origin in enumerate(self.states)]
                full.append(self._column([(column, 1, None)] + normalization))
            self.N.append(normalizer)
            self.F.append(hamiltonian)
            self.X.append(projected)
            self.Y.append(full)
        return self

    def intermediate_columns(self, order):
        """Return Omega=Y(PY)^-1; P Omega P=P at every retained order."""
        self.extend(order)
        while len(self._intermediate_columns) <= order:
            n = len(self._intermediate_columns)
            inverse = -sum(self.N[i] @ self._inverse_normalizers[n-i]
                           for i in range(1,n+1))
            self._inverse_normalizers.append(inverse)
            columns = []
            for a, state in enumerate(self.states):
                parts = [(self.X[i][b], self._inverse_normalizers[n-i][b,a], None)
                         for i in range(1,n+1) for b in range(len(self.states))]
                columns.append(self._column(parts))
            self._intermediate_columns.append(columns)
        return self._intermediate_columns[order]

    def operator_piece(self, order, column, off_diagonal=False):
        key = (order, column, off_diagonal)
        if key not in self._operator_pieces:
            if order == 0:
                result = {} if off_diagonal else self.Y[0][column]
            elif off_diagonal:
                result = self._column([({state: self.identity}, self.N[order][b,column], None)
                                       for b,state in enumerate(self.states) if b != column])
            else:
                result = self._column([(self.X[order][column],1,None),
                    ({self.states[column]:self.identity},self.N[order][column,column],None)])
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
                            terms=None, branches=True):
        self.extend(order); other.extend(order)
        allowed = None if terms is None else {tuple(t) for t in terms}
        result = 0.
        for i in range(order+1):
            for k in range(order+1-i):
                j = order-i-k
                if j < len(operators) and (allowed is None or (k, j, i) in allowed):
                    if branches is True:
                        for bra,bc in self.operator_components(i,left):
                            for ket,kc in other.operator_components(k,right):
                                result += bc*kc*self.inner(bra,ket,operators[j],expansions)
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
                                        result += bc*kc*self.inner(bra,ket,operators[j],expansions)
        return result


class ProjectedCorrectionEvaluator:
    """Resolve actual blocks and expose the ordinary analytic result layout."""
    def __init__(self, evaluator, blocks=None, pairs=None, expansions=None,
                 freqs=None, zero_cutoff=None):
        self.evaluator = evaluator
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

    def block(self, state, order):
        state = tuple(int(n) for n in state)
        group = self.blocks.get(state, (state,))
        if group not in self.cache:
            self.cache[group] = AnalyticBlockExpansion(self.evaluator, group,
                self.expansions, self.freqs, self.zero_cutoff)
        block = self.cache[group].extend(order)
        return block, block.states.index(state)

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
            corrections.append([None if target_orders is not None and n not in target_orders
                                else np.asarray([value(initial, final, n) for final in finals])
                                for n in range(order+1)])
        return BasicAPTCorrections(initial_space, final_spaces, corrections)

    def wavefunctions(self, states, order, full=True, target_orders=None):
        targets = None if target_orders is None else set(target_orders)
        def value(initial, final, n):
            if targets is not None and n not in targets:
                return 0.
            block, a = self.block(initial, order)
            column = block.Y[n][a] if full else block.X[n][a]
            if self.evaluator.solver.intermediate_normalization:
                column = block.intermediate_columns(n)[a]
                if not full and n == 0:
                    column = {}
            return -sum(block.element(term, origin, final) for origin, term in column.items())
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
                    target_orders=None):
        def value(initial, final, n):
            block, a = self.block(initial, order)
            final = tuple(final)
            if final not in block.states:
                raise ValueError('Hamiltonian block endpoints must belong to the same model space')
            b = block.states.index(final)
            result = (block.F[n][b,a] if terms is None and branches is True else
                      block.operator_correction(block,b,a,self.evaluator.solver.hamiltonian_expansion,
                                                n,terms=terms,branches=self._operator_branches(branches)))
            return result if only_degenerate_terms else -result
        return self.result(states, order, value, target_orders=target_orders)

    def operators(self, states, order, operators, expansions, terms=None, branches=True,
                  target_orders=None):
        def value(initial, final, n):
            right, ket = self.block(initial, order)
            left, bra = self.block(final, order)
            return -left.operator_correction(right, bra, ket, operators, n,
                    expansions=expansions, terms=terms,branches=self._operator_branches(branches))
        return self.result(states, order, value, target_orders=target_orders)
