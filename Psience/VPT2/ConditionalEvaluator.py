"""Use the native batched evaluator for condition-based VPT expressions."""
import numpy as np
from .Analytic import PerturbationTheoryTerm, PerturbationTheoryExpressionEvaluator


class CorrectionView(PerturbationTheoryTerm):
    """A final native filtering/checkpoint view without mutating a shared DAG."""
    def __init__(self,term,diagonal=False,**options):
        super().__init__(logger=term.logger,polynomial_representation=term.polynomial_representation,
            **{key:value for key,value in options.items() if key in
               {'allowed_terms','allowed_energy_changes','allowed_coefficients','disallowed_coefficients'}})
        self.term = term
        self.diagonal = bool(diagonal)

    def get_changes(self):
        return {change:None for change in self.term.changes if not self.diagonal or not change}

    def get_core_poly(self,changes,shift=None):
        return self.term.get_poly_terms(changes)

    def __repr__(self):
        return 'ConditionalCorrection[{}]'.format(self.term)


class ConditionalCorrectionEvaluator:
    """Late-bound occupations/coefficients; no finite block contractions."""
    def __init__(self,evaluator,expansions=None,freqs=None,zero_cutoff=None,parallelizer=None):
        self.evaluator = evaluator
        self.stack = evaluator.solver.get_conditional_stack()
        self.raw_expansions = evaluator.expansions if expansions is None else expansions
        self.expansions = evaluator._prep_expansions(self.raw_expansions)
        self.freqs = evaluator.freqs if freqs is None else freqs
        self.zero_cutoff = zero_cutoff
        self.parallelizer = evaluator.parallelizer if parallelizer is None else parallelizer
        self.options = {}

    def validate_options(self,options):
        supported = {'allowed_coefficients','disallowed_coefficients','epaths',
            'only_degenerate_terms','include_degenerate_correction_terms','zero_cutoff',
            'parallelizer','target_orders','return_sorted','logger','log_scaled','terms'}
        unknown = set(options)-supported
        if unknown:
            raise TypeError('unknown conditional correction options: {}'.format(sorted(unknown)))
        self.options = {key:value for key,value in options.items() if key in
                        {'allowed_coefficients','disallowed_coefficients','epaths','return_sorted'}}

    def _diagonal_values(self,term,states,**numerical):
        changes = self.evaluator.get_diff_map([[state,[state]] for state in states])
        values = []
        for change,permutations in changes.items():
            expression = PerturbationTheoryExpressionEvaluator(term,term.get_poly_terms(change),change)
            values.extend(expression.evaluate(permutations,self.expansions,self.freqs,
                zero_cutoff=self.zero_cutoff,parallelizer=self.parallelizer,log_scaled=False,
                **numerical))
        return np.asarray(values)[:,0]

    def energies(self,states,order,**numerical):
        self.stack.extend(order)
        return [self._diagonal_values(self.stack.hamiltonian[n],states,**numerical)
                for n in range(0,order+1,2)]

    def overlaps(self,states,order,**numerical):
        self.stack.extend(order)
        return np.stack([self._diagonal_values(self.stack.wavefunctions[n],states,**numerical)
                         for n in range(order+1)],axis=-1)

    def _evaluate(self,series,states,order,expansions=None,target_orders=None,sign=1.):
        def generator(n,**options):
            term = series[n] if sign == 1 else sign*series[n]
            return CorrectionView(term,**options)
        return self.evaluator.get_state_by_state_corrections(generator,states,order=order,
            expansions=self.raw_expansions if expansions is None else expansions,freqs=self.freqs,
            zero_cutoff=self.zero_cutoff,parallelizer=self.parallelizer,
            only_degenerate_terms=False,include_degenerate_correction_terms=False,
            target_orders=target_orders,**self.options)

    def wavefunctions(self,states,order,full=True,target_orders=None):
        self.stack.extend(order)
        if self.evaluator.solver.intermediate_normalization:
            series = self.stack.intermediate_wavefunctions(order)
        else:
            series = self.stack.wavefunctions
        if not full:
            series = [self.stack.restrict(term,retain=False) for term in series]
        return self._evaluate(series,states,order,target_orders=target_orders)

    def hamiltonian(self,states,order,only_degenerate_terms=True,terms=None,branches=True,
                    target_orders=None):
        self.stack.extend(order)
        series = (self.stack.hamiltonian if terms is None and branches is True else
                  self.stack.transform_operator(self.evaluator.solver.hamiltonian_expansion,
                                                order,terms,branches))
        return self._evaluate(series,states,order,target_orders=target_orders,
                              sign=-1. if only_degenerate_terms else 1.)

    def operators(self,states,order,operators,expansions,terms=None,branches=True,
                  target_orders=None):
        series = self.stack.transform_operator(operators,order,terms,branches)
        return self._evaluate(series,states,order,expansions,target_orders)
