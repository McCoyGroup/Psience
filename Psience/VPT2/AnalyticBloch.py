"""Condition-based analytic Bloch wave operator, using native projected terms.

Omega_0=I, D Omega_n=0 for n>0. With V=H-H0:
  T_n = sum_{r=1}^n H_r Omega_{n-r}
  K_n = D T_n
  Omega_n = R_Q (T_n - sum_{r=1}^{n-1} Omega_{n-r} K_r).
R_Q is the native conditional homological inverse; it cannot propagate through
the retained relation. No perturbation order is singled out.
"""
from .Conditional import (ConditionalTransformationStack,SameStateCondition,
    operator_sum,operator_product,ZeroTerm)
from .Analytic import ShiftedEnergyBaseline

class AnalyticBlochSeries:
    def __init__(self,solver,condition=None):
        condition=SameStateCondition() if condition is None else condition
        if not condition.equivalence:
            raise ValueError('Bloch model spaces require an equivalence relation')
        self.projector=ConditionalTransformationStack(solver,condition)
        self.condition=self.projector.condition
        self.h0=solver.hamiltonian_expansion[0]
        self.expansion=solver.hamiltonian_expansion
        self.identity,self.zero=self.projector.identity,self.projector.zero
        self.wave_operator=[self.identity]
        self.effective_hamiltonian=[self.h0]

    def extend(self,order):
        if order>=len(self.expansion):raise ValueError('Hamiltonian expansion too short')
        while len(self.wave_operator)<=order:
            n=len(self.wave_operator)
            driving=operator_sum(*(operator_product(self.expansion[r],self.wave_operator[n-r])
                                  for r in range(1,n+1)))
            correction=self.projector.restrict(driving)
            feedback=operator_sum(*(operator_product(self.wave_operator[n-r],self.effective_hamiltonian[r])
                                   for r in range(1,n)))
            self.effective_hamiltonian.append(correction)
            self.wave_operator.append(self.projector.restrict(operator_sum(driving,-feedback),
                                                            retain=False,resolvent=True))
        return self

    def metric(self,order):
        self.extend(order)
        result=[self.identity]
        for n in range(1,order+1):
            result.append(self.projector.restrict(operator_sum(*(
                operator_product(ShiftedEnergyBaseline.lookup(self.wave_operator[i]),self.wave_operator[n-i])
                for i in range(n+1)))))
        return result

    def _product_series(self,left,right,order):
        return [operator_sum(*(operator_product(left[i],right[n-i])
                    for i in range(n+1) if i<len(left) and n-i<len(right))) for n in range(order+1)]

    def _metric_power(self,order,power):
        delta=self.metric(order);delta[0]=self.zero
        result=[self.identity]+[self.zero]*order
        product=[self.identity]+[self.zero]*order;coefficient=1.
        for exponent in range(1,order+1):
            product=self._product_series(product,delta,order)
            coefficient*= (power-exponent+1)/exponent
            for n in range(1,order+1):result[n]=operator_sum(result[n],coefficient*product[n])
        return result

    def normalized_wavefunctions(self,order):
        self.extend(order)
        return self._product_series(self.wave_operator,self._metric_power(order,-.5),order)

    def hermitian_hamiltonian(self,order):
        self.extend(order)
        if isinstance(self.condition,SameStateCondition):
            # Each retained block is a scalar, so the metric similarity cancels.
            return self.effective_hamiltonian[:order+1]
        return self._product_series(self._metric_power(order,.5),
            self._product_series(self.effective_hamiltonian,self._metric_power(order,-.5),order),order)
