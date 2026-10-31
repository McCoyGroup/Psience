import itertools
import unittest
import warnings

import numpy as np
from McUtils.Numputils import SparseArray
from McUtils.Scaffolding import NullLogger,NullCheckpointer
from McUtils.Data import UnitsData
from Psience.BasisReps import HarmonicOscillatorProductBasis,BasisStateSpace
from Psience.VPT2 import AnalyticVPTRunner
from Psience.VPT2.Analytic import SmallEnergyDenominatorWarning
from Psience.VPT2.Corrections import AnalyticPerturbationTheoryCorrections
from Psience.VPT2.Solver import PerturbationTheorySolver

def tensor(n,rank,terms):
    a=np.zeros((n,)*rank)
    for coefficient,indices in terms:
        perms=set(itertools.permutations(indices))
        for p in perms:a[p]=coefficient/len(perms)
    return a

def analytic(w,v3,v4,dipole=None):
    w=np.array(w)
    return AnalyticVPTRunner([[np.diag(w)/2,np.diag(w)/2],[v3,0],[v4,0,0,0,0]],
        order=2,freqs=w,internals=False,logger=False,dipole_expansion=dipole)

class VPT2ResonanceProjectionTests(unittest.TestCase):
    def test_separable_and_spectator_invariance(self):
        a,b=[1,0,0],[0,1,0]
        v3=tensor(3,3,[(.002,(0,0,0)),(.002,(1,1,1))])
        v4=tensor(3,4,[(.0002,(0,0,0,0)),(.0002,(1,1,1,1))])
        for spectator in [.002,.2]:
            r=analytic([.0172,.0173,spectator],v3,v4)
            space=r.prep_states([[0,0,0],a,b],degeneracy_specs={'groups':[[a,b]]})
            for specific in [False,True]:
                for only in [False,True]:
                    with self.subTest(spectator=spectator,specific=specific,only=only):
                        _,corrections=r.get_reexpressed_hamiltonian(space,order=2,
                            state_specific_degeneracies=specific,only_degenerate_terms=only)
                        self.assertLess(np.max(abs(np.array(corrections[0])[:,0,1])),1e-12)

    def test_exact_fermi_never_divides_by_forbidden_denominator(self):
        a,b=[1,0,0],[0,2,0]
        g=.0006
        r=analytic([.0344,.0172,.2],tensor(3,3,[(g,(0,1,1))]),np.zeros((3,)*4))
        space=r.prep_states([[0,0,0],a,b],degeneracy_specs={'groups':[[a,b]]})
        for specific in [False,True]:
            degs=r._eval_degeneracies(space,specific)
            with self.subTest(specific=specific),warnings.catch_warnings():
                warnings.simplefilter('error',SmallEnergyDenominatorWarning)
                for only in [False,True]:
                    _,blocks=r.get_reexpressed_hamiltonian(space,order=2,
                        state_specific_degeneracies=specific,only_degenerate_terms=only)
                    self.assertTrue(np.isfinite(blocks).all())
                    physical=AnalyticPerturbationTheoryCorrections(
                        space.flat_space.state_space,space.state_list_pairs,logger=NullLogger())
                    physical.degenerate_states=space.flat_space.degenerate_states
                    physical.degenerate_hamiltonian_corrections=blocks
                    physical.only_degenerate_terms=only
                    energies=(space.flat_space.state_space.excitations+.5)@r.eval.freqs
                    _,(hams,_)=physical.get_degenerate_transformations(
                        space.flat_space.state_space,energies)
                    np.testing.assert_allclose(hams[0][0,1],g/2,atol=1e-13)
                wf=r.eval.get_wavefunction_corrections([[a,[b]],[b,[a]]],order=1,
                    degenerate_states=degs,target_orders=[1])
                for block in wf.corrections:np.testing.assert_array_equal(block[1],0.)

    def test_coupled_operator_and_hamiltonian_match_classic(self):
        w=np.array([.0201,.0217,.0313]);g=[0,0,0];a=[1,0,0];b=[0,1,0];c=[0,0,1]
        states=[g,a,b,c]
        v3=tensor(3,3,[(.0003,(0,0,0)),(.0004,(1,1,1)),(.00025,(2,2,2)),
            (.00008,(0,1,1)),(.0001,(0,0,1)),(.00007,(0,1,2))])
        v4=tensor(3,4,[(.00001,(0,0,0,0)),(.00002,(0,0,1,1))])
        basis=HarmonicOscillatorProductBasis(3)
        def rep(rank,coeff):
            return basis.representation(*(['x']*rank),coeffs=coeff,
                axes=[list(range(rank)),list(range(rank))])
        h0=rep(2,np.diag(w)/2)+basis.representation('p','p',coeffs=np.diag(w)/2,axes=[[0,1],[0,1]])
        solver=PerturbationTheorySolver([h0,rep(3,v3),rep(4,v4)],BasisStateSpace(basis,states),
            order=2,degenerate_states=[np.array([a,b])],handle_strong_couplings=False,
            logger=NullLogger(),checkpointer=NullCheckpointer(None),results=NullCheckpointer(None),
            zero_element_warning=False)
        corr=solver.apply_VPT()
        mu1=np.array([.2,.1,.3]);mu2=np.diag([.02,-.01,.01]);mu2[0,2]=mu2[2,0]=.003
        mu3=tensor(3,3,[(.0004,(0,0,0)),(-.0002,(1,1,1))])
        operators=[rep(1,mu1),rep(2,mu2/2),rep(3,mu3/6)]
        matrices=[op.get_representation_matrix(corr.total_basis,corr.total_basis,
            zero_element_warning=False) for op in operators]
        reference=[v.asarray() if isinstance(v,SparseArray) else v for v in corr.operator_representation(matrices)]
        native_initial=int(corr.states.find([g])[0])
        dipole=[[0.,mu1,mu2,mu3]]+[[0.,np.zeros(3),np.zeros((3,3)),np.zeros((3,3,3))] for _ in range(2)]
        for specific in [False,True]:
            r=analytic(w,v3,v4,dipole)
            space=r.prep_states([[states,states]],degeneracy_specs={'groups':[[a,b]]})
            _,h=r.get_reexpressed_hamiltonian(space,order=2,state_specific_degeneracies=specific)
            np.testing.assert_allclose(np.array(h[0])[:,0,1],np.array(corr.degenerate_hamiltonians[0])[:,0,1],atol=2e-12)
            moments=r.get_transition_moment_corrections(space,order=2,state_specific_degeneracies=specific)
            finals=space.state_list_pairs[0][1]
            positions=corr.states.find(finals)
            for initial_index,initial in enumerate(space.state_list_pairs[0][0]):
                native_row=int(corr.states.find([initial])[0])
                for order in range(3):
                    np.testing.assert_allclose(moments[0][0][order][initial_index],-np.array(reference[order])[native_row,positions],atol=2e-11)
            # A quadratic zeroth-order operator exercises normalization on
            # both legs within the same resonance block, not just g -> P.
            quadratic=[mu2/2,mu3/6,np.zeros((3,)*4)]
            qm=[rep(k+2,v).get_representation_matrix(corr.total_basis,corr.total_basis,
                zero_element_warning=False) for k,v in enumerate(quadratic)]
            qref=[v.asarray() if isinstance(v,SparseArray) else v for v in corr.operator_representation(qm)]
            qr=analytic(w,v3,v4)
            qs=qr.prep_states([[states,states]],degeneracy_specs={'groups':[[a,b]]})
            qcorr=qr.get_operator_corrections([quadratic],qs,order=2,
                operator_type=2,check_single=False,state_specific_degeneracies=specific)
            for initial_index,initial in enumerate(qs.state_list_pairs[0][0]):
                native_row=int(corr.states.find([initial])[0])
                for order in range(3):
                    np.testing.assert_allclose(qcorr[0][0][order][initial_index],-np.array(qref[order])[native_row,positions],atol=2e-11)



"""Additional portable regressions appended by build_generic_patch.py."""
import pickle
import json
from Psience.VPT2.Analytic import (
    AnalyticPerturbationTheorySolver, PerturbationTheoryEvaluator,
    HamiltonianExpansionTerm, SqrtChangePoly, PerturbationTheoryExpressionEvaluator,
)
from Psience.VPT2.Projected import AnalyticBlockExpansion, StateProjectionFactor


class AnalyticProjectedBlockTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.w = np.array([.04, .091, .137])
        cls.force = .004
        parent = AnalyticPerturbationTheorySolver.from_order(4, internals=False, logger=False)
        terms = list(parent.hamiltonian_expansion)
        terms[1] = HamiltonianExpansionTerm([['x']], order=1, identities=[0], logger=False)
        solver = AnalyticPerturbationTheorySolver(terms, logger=False)
        coefficients = [[np.diag(cls.w)/2, np.diag(cls.w)/2],
                        [np.array([cls.force, 0., 0.])], [0.], [0.], [0.]]
        cls.evaluator = PerturbationTheoryEvaluator(solver, coefficients, freqs=cls.w)
        cls.block = AnalyticBlockExpansion(cls.evaluator, [(0,0,0), (1,0,0)],
                                          zero_cutoff=1e-30).extend(4)
        cls.position = (1/cls.force) * terms[1]

    def test_block_feedback_in_second_order_fermi_wavefunction(self):
        a, b, k = [1,0,0], [0,2,0], [0,1,0]
        cubic = tensor(3, 3, [(.0006, (0,1,1)), (.0004, (1,1,1))])
        runner = analytic([.04,.02,.13], cubic, np.zeros((3,)*4))
        result = runner.eval.get_wavefunction_corrections([[a,[k]]], order=2,
                                                          degenerate_blocks=[[a,b]])
        # The generic evaluator has the established overall minus sign.
        np.testing.assert_allclose(result.corrections[0][2], [.0016125], atol=1e-13)

    def test_fourth_order_orthogonality_including_odd_block_normalization(self):
        block = self.block
        self.assertGreater(abs(block.N[3][0,1]), 0.)
        for n in range(5):
            gram = np.array([[sum(block.inner(block.Y[i][a], block.Y[n-i][b])
                                  for i in range(n+1)) for b in range(2)] for a in range(2)])
            np.testing.assert_allclose(gram, np.eye(2) if n == 0 else 0., atol=2e-12)
            np.testing.assert_allclose(block.F[n], block.F[n].T, atol=2e-12)

    def test_fourth_order_exact_oscillator_energies_and_dipoles(self):
        block = self.block
        moments = [np.array([[block.operator_correction(block, a, b, [self.position], n)
                              for b in range(2)] for a in range(2)]) for n in range(5)]
        energy_errors, dipole_errors = [], []
        for strength in [1., .5]:
            h = sum(value * strength**n for n, value in enumerate(block.F))
            energies, rotation = np.linalg.eigh(h)
            exact = np.diag(block.F[0]) - strength**2 * self.force**2/(2*self.w[0])
            energy_errors.append(np.max(abs(energies-exact)))
            mu = rotation.T @ sum(value * strength**n for n, value in enumerate(moments)) @ rotation
            dipole_errors.append(max(np.max(abs(np.diag(mu)+strength*self.force/self.w[0])),
                                     abs(abs(mu[0,1])-1/np.sqrt(2))))
        self.assertLess(energy_errors[1], energy_errors[0]/20)
        self.assertLess(dipole_errors[1], dipole_errors[0]/20)

    def test_intermediate_normalization_has_no_higher_order_P_components(self):
        block = self.block
        for n in range(1,5):
            for column in block.intermediate_columns(n):
                for final in block.states:
                    value = sum(block.element(term, origin, final) for origin, term in column.items())
                    self.assertEqual(value, 0.)

    def test_nonzero_hamiltonian_terms_at_every_order(self):
        forces = [0., .004, .0003, -.0005, .0002]
        parent = AnalyticPerturbationTheorySolver.from_order(4, internals=False, logger=False)
        terms = list(parent.hamiltonian_expansion)
        for n in range(1,5):
            terms[n] = HamiltonianExpansionTerm([['x']], order=n, identities=[0], logger=False)
        coefficients = [[np.diag(self.w)/2, np.diag(self.w)/2]] + [
            [np.array([force,0.,0.])] for force in forces[1:]]
        evaluator = PerturbationTheoryEvaluator(
            AnalyticPerturbationTheorySolver(terms,logger=False), coefficients, freqs=self.w)
        block = AnalyticBlockExpansion(evaluator,self.block.states,zero_cutoff=1e-30).extend(4)
        position = (1/forces[1])*terms[1]
        moments = [np.array([[block.operator_correction(block,a,b,[position],n)
                              for b in range(2)] for a in range(2)]) for n in range(5)]
        energy_errors, moment_errors = [], []
        for strength in [1.,.5,.25]:
            total_force = sum(force*strength**n for n,force in enumerate(forces))
            energies, rotation = np.linalg.eigh(sum(f*strength**n for n,f in enumerate(block.F)))
            exact = np.diag(block.F[0])-total_force**2/(2*self.w[0])
            energy_errors.append(np.max(abs(energies-exact)))
            mu = rotation.T @ sum(m*strength**n for n,m in enumerate(moments)) @ rotation
            moment_errors.append(max(np.max(abs(np.diag(mu)+total_force/self.w[0])),
                                     abs(abs(mu[0,1])-1/np.sqrt(2))))
        for errors in [energy_errors,moment_errors]:
            self.assertLess(errors[1],errors[0]/20)
            self.assertLess(errors[2],errors[1]/20)

    def test_projector_metadata_survives_expression_serialization(self):
        block = self.block
        initial, final = (1,0,0), (2,0,0)
        term = block.X[1][1][initial]
        original = term.get_poly_terms((1,))
        eager = original.mutate(poly_obj=original.poly_obj.to_eager())
        # Cross a real text-serialization boundary; a numeric tree must not
        # silently convert a projector tag to a plain integer.
        restore = lambda poly: SqrtChangePoly.from_state(json.loads(json.dumps(
            poly.to_state(),default=lambda value: value.tolist())))
        restored = [restore(poly) for poly in [original,eager]]
        for change, permutations in self.evaluator.get_diff_map([[initial,[final]]]).items():
            for polynomial in [original] + restored:
                evaluator = PerturbationTheoryExpressionEvaluator(term, polynomial, change)
                for backend in ['legacy', 'indexed', 'compare']:
                    result = evaluator.evaluate(permutations, block.expansions, self.w,
                                                evaluation_backend=backend)
                    np.testing.assert_allclose(result[0], [-self.force/self.w[0]], atol=1e-13)
        factor = StateProjectionFactor('R', block.states, self.w, initial)
        restored = pickle.loads(pickle.dumps(factor))
        self.assertIsInstance(restored, StateProjectionFactor)
        self.assertEqual(restored.states, block.states)
        self.assertEqual(restored.reference, initial)

    def test_projected_block_recursion_matches_both_numerical_backends(self):
        blocks = [AnalyticBlockExpansion(self.evaluator, self.block.states,
                   zero_cutoff=1e-30, evaluation_backend=backend).extend(3)
                  for backend in ['legacy', 'indexed']]
        np.testing.assert_allclose(blocks[0].F, blocks[1].F, atol=2e-14)
        np.testing.assert_allclose(blocks[0].N, blocks[1].N, atol=2e-14)

    def test_operator_branch_decomposition_at_arbitrary_order(self):
        block = self.block
        cross = []
        for n in [2,3,4]:
            full = block.operator_correction(block,0,0,[self.position],n)
            parts = []
            for name in ['main','left','right','both']:
                selection = {'main':False, name:True}
                parts.append(block.operator_correction(block,0,0,[self.position],n,branches=selection))
            self.assertAlmostEqual(full,sum(parts),places=13)
            cross.extend(parts[1:])
        self.assertGreater(max(abs(v) for v in cross),0.)

    def test_public_operator_term_orders_use_initial_ket_and_final_bra(self):
        quadratic = np.diag([1.,0.,0.])
        coefficients = [[quadratic,np.zeros((3,)*3)]]
        pairs = [[self.block.states[0],[self.block.states[1]]]]
        def correction(term):
            result = self.evaluator.get_operator_corrections(coefficients,pairs,order=1,
                operator_type=2,check_single=False,terms=[term],degenerate_blocks=[self.block.states])
            return result[0].corrections[0][1][0]
        self.assertAlmostEqual(correction((1,0,0)),0.,places=13)
        self.assertAlmostEqual(correction((0,0,1)),self.force/(np.sqrt(2)*self.w[0]),places=13)

    def test_block_cache_tracks_coefficient_changes_and_clear_caches(self):
        a,b=[1,0,0],[0,2,0]
        cubic=tensor(3,3,[(.0006,(0,1,1))])
        runner=analytic([.04,.02,.13],cubic,np.zeros((3,)*4))
        def coupling():
            value=runner.eval.get_reexpressed_hamiltonian([a,b],order=1,
                degenerate_blocks=[[a,b]],only_degenerate_terms=True)
            return value.corrections[0][1][0]
        self.assertAlmostEqual(coupling(),.0003,places=13)
        runner.eval.expansions[1][0] *= 2
        self.assertAlmostEqual(coupling(),.0006,places=13)
        AnalyticPerturbationTheorySolver.clear_caches()
        self.assertFalse(hasattr(runner.eval,'_projected_block_cache'))
        self.assertAlmostEqual(coupling(),.0006,places=13)


if __name__ == "__main__":
    unittest.main()
