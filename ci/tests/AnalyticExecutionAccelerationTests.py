"""Differential tests for execution-only analytic accelerations."""
import unittest,itertools,warnings
from unittest.mock import patch
import numpy as np
from Psience.VPT2 import AnalyticVPTRunner
from Psience.VPT2.Conditional import SameStateCondition,PolynomialInvariantCondition,IntegerLatticeCondition
from Psience.VPT2.IndexedEvaluator import IndexedBlockEvaluator as Indexed
from Psience.VPT2.Analytic import (PerturbationTheoryExpressionEvaluator as Expression,
    AnalyticPerturbationTheorySolver,HamiltonianExpansionTerm,PerturbationTheoryEvaluator)
from Psience.VPT2.Conditional import ZeroTerm
from McUtils.Zachary import Symbols

class AccelerationTests(unittest.TestCase):
    def run_energy(self,nmodes=4,condition=None,batch=True,jit=True,max_states=128):
        rng=np.random.default_rng(102)
        w=np.sqrt(np.arange(2,nmodes+2))*.017
        c3=rng.normal(size=(nmodes,)*3)*.00006
        c4=rng.normal(size=(nmodes,)*4)*.000004
        expansions=[[np.diag(w)/2,np.diag(w)/2],[c3,0],[c4,0,0,0,0]]
        r=AnalyticVPTRunner(expansions,order=2,freqs=w,internals=False,
            logger=False,resonance_condition=condition)
        states=np.vstack([np.zeros((1,nmodes),int),np.eye(nmodes,dtype=int),
            np.full((1,nmodes),2,dtype=int)])
        with patch.object(Indexed,'use_batched_state_routing',batch),\
             patch.object(Indexed,'use_unmasked_selection',batch),\
             patch.object(Indexed,'use_jit_energy_contractions',jit),\
             patch.object(Indexed,'max_batched_states',max_states):
            result=np.asarray(r.eval.get_energy_corrections(states,order=2,evaluation_mode='dag'))
            stats=Expression.get_last_indexed_evaluation_stats()
        return result,stats

    def test_batched_energy_matches_unbatched_and_respects_state_cap(self):
        reference,_=self.run_energy(condition=SameStateCondition(),batch=False,jit=False)
        for cap in [2,128]:
            actual,stats=self.run_energy(condition=SameStateCondition(),batch=True,max_states=cap)
            np.testing.assert_allclose(actual,reference,atol=2e-14,rtol=2e-12)
            self.assertGreater(stats['batched_state_blocks'],0)
            self.assertGreater(stats['batched_source_states'],stats['batched_state_blocks'])

    def test_ordinary_analytic_recurrence_matches(self):
        a,_=self.run_energy(batch=False,jit=False)
        b,_=self.run_energy(batch=True)
        np.testing.assert_allclose(a,b,atol=2e-14,rtol=2e-12)

    def test_occupation_condition_retains_projected_fallback(self):
        q=Symbols('q').q
        condition=PolynomialInvariantCondition(q*(q-1),[0])
        a,_=self.run_energy(condition=condition,batch=False,jit=False)
        b,_=self.run_energy(condition=condition,batch=True)
        np.testing.assert_allclose(a,b,atol=2e-14,rtol=2e-12)

    def test_jit_energy_kernel_shapes_and_empty_products(self):
        rng=np.random.default_rng(71)
        for rows,perms,modes,factors in [(1,1,0,0),(23,1,3,1),(17,7,4,3),(2,4,3,0)]:
            changes=rng.integers(-3,4,size=(rows,factors,modes+1))
            frequencies=rng.uniform(.01,.1,size=(rows,perms,modes))
            with patch.object(Indexed,'use_jit_energy_contractions',False):a=Indexed._contract_energy_rows(changes,frequencies)
            with patch.object(Indexed,'use_jit_energy_contractions',True):b=Indexed._contract_energy_rows(changes,frequencies)
            np.testing.assert_allclose(a,b,atol=2e-18,rtol=3e-14)

    def test_no_numba_uses_native_energy_fallback(self):
        rng=np.random.default_rng(3)
        f=rng.random((3,4,2));c=rng.integers(-2,3,size=(3,2,3))
        with patch.object(Indexed,'use_jit_energy_contractions',False):a=Indexed._contract_energy_rows(c,f)
        with patch.object(Indexed,'use_jit_energy_contractions',True),\
             patch('Psience.VPT2._IndexedKernels.get_energy_row_kernel',return_value=None):b=Indexed._contract_energy_rows(c,f)
        np.testing.assert_array_equal(a,b)

    def test_fourth_order_energies_match_with_batching(self):
        rng=np.random.default_rng(66);w=np.array([.037,.0613])
        h0=AnalyticPerturbationTheorySolver.from_order(2,internals=False,logger=False).hamiltonian_expansion[0]
        zero=ZeroTerm(logger=False)
        h2=HamiltonianExpansionTerm([['x','x']],order=2,identities=[0],logger=False)
        hamiltonian=[h0,zero,h2,zero,zero]
        expansions=[[np.diag(w)/2,np.diag(w)/2],[0],[rng.normal(size=(2,2))*1e-5],[0],[0]]
        states=np.array([[0,0],[1,0],[0,1],[2,0]])
        results=[]
        for batching in [False,True]:
            solver=AnalyticPerturbationTheorySolver(hamiltonian,logger=False,resonance_condition=SameStateCondition())
            r=AnalyticVPTRunner(PerturbationTheoryEvaluator(solver,expansions,freqs=w),order=4)
            with patch.object(Indexed,'use_batched_state_routing',batching),\
                 patch.object(Indexed,'use_unmasked_selection',batching),\
                 patch.object(Indexed,'use_jit_energy_contractions',False):
                results.append(np.asarray(r.eval.get_energy_corrections(states,order=4,evaluation_mode='dag')))
        np.testing.assert_allclose(*results,atol=3e-13,rtol=3e-11)
        self.assertEqual(results[1].shape,(3,4))
        self.assertGreater(np.max(np.abs(results[1][-1])),1e-15)

    def test_batched_multiple_coefficient_expansions_and_irregular_permutations(self):
        solver=AnalyticPerturbationTheorySolver.from_order(2,internals=False,logger=False)
        evaluator=solver.energy_correction(2)([])
        keys={c for product in evaluator.expr.poly_obj.to_eager().terms for c in product}
        rng=np.random.default_rng(720)
        tree=[]
        for kind in range(max(k[0] for k in keys)+1):
            orders=[k[1] for k in keys if k[0]==kind]
            values=[]
            for order in range(max(orders,default=0)+1):
                ranks=[len(k)-2 for k in keys if k[:2]==(kind,order)]
                values.append(0 if not ranks else rng.normal(size=(3,)*max(ranks))*.0002)
            tree.append(values)
        coefficients=[tree,[[np.asarray(v)*.7 for v in values] for values in tree]]
        w=np.array([.041,.066,.097]);identity=np.arange(3)[None,:]
        for irregular in [False,True]:
            states=[[np.array([0,0,0]),identity],[np.array([1,0,0]),identity],
                    [np.array([0,1,0]),np.array([[1,0,2]]) if irregular else identity]]
            reference=evaluator.evaluate(states,coefficients,w,evaluation_mode='dag',evaluation_backend='legacy')
            actual=evaluator.evaluate(states,coefficients,w,evaluation_mode='dag',evaluation_backend='indexed')
            for a,b in zip(actual,reference):np.testing.assert_allclose(a,b,atol=2e-14,rtol=2e-12)
            stats=Expression.get_last_indexed_evaluation_stats()
            self.assertEqual(stats['batched_state_blocks']>0,not irregular)

    def test_batch_limit_keeps_a_single_state_tail_and_input_order(self):
        w=np.array([.04])
        expansions=[[np.diag(w)/2,np.diag(w)/2],[np.full((1,1,1),1e-5),0],
                    [np.full((1,1,1,1),1e-6),0,0,0,0]]
        states=np.random.default_rng(812).permutation(129)[:,None]
        results=[]
        for batching in [False,True]:
            r=AnalyticVPTRunner(expansions,order=2,freqs=w,internals=False,
                resonance_condition=SameStateCondition(),logger=False)
            with patch.object(Indexed,'use_batched_state_routing',batching),\
                 patch.object(Indexed,'max_batched_states',128),\
                 patch.object(Indexed,'use_jit_energy_contractions',False):
                results.append(np.asarray(r.eval.get_energy_corrections(states,order=2,evaluation_mode='dag')))
        np.testing.assert_allclose(*results,atol=2e-13,rtol=2e-12)
        self.assertEqual(results[1].shape,(2,129))

if __name__=='__main__':unittest.main()
