import json
import pickle
import unittest
import warnings
import numpy as np

from McUtils.Zachary import Symbols
from Psience.BasisReps import HarmonicOscillatorProductBasis, BasisStateSpace
from Psience.VPT2 import AnalyticVPTRunner
from Psience.VPT2.Analytic import (
    AnalyticPerturbationTheorySolver, PerturbationTheoryEvaluator,
    HamiltonianExpansionTerm, PerturbationTheoryExpressionEvaluator,
    ShiftedEnergyBaseline, ScaledPerturbationTheoryTerm,
    PerturbationTheoryTermSum, PerturbationTheoryTermProduct,
    PTTensorCoeffProductSum, PTTensorCoeffProductDAG,
    PTTensorCoeffProductDAGEvaluationPlan, SqrtChangePoly,
    make_product_polynomial, SmallEnergyDenominatorWarning, PolyPath,
)
from Psience.VPT2.Conditional import (
    PolynomialInvariantCondition, PolynomialRelationshipCondition,
    IntegerLatticeCondition, SameStateCondition, EnergyWindowCondition,
    ConditionalOperatorTerm, ConditionalTransformationStack, ConditionEndpoint,
    ZeroTerm, operator_sector,
)
from Psience.VPT2.Projected import IdentityTerm


class MatrixInterpreter:
    """Test-only interpretation of native operator terms in a Psience basis."""
    def __init__(self,primitive,states,frequencies):
        self.primitive,self.states,self.frequencies = primitive,states,frequencies
        self.identity = np.eye(len(states))
        self.zero = np.zeros_like(self.identity)
        self.energies = (np.asarray(states)+.5) @ frequencies
        self.memo,self.masks = {},{}

    def __call__(self,term):
        if term in self.memo:
            return self.memo[term]
        if term in self.primitive:
            result = self.primitive[term]
        elif isinstance(term,IdentityTerm):
            result = self.identity
        elif isinstance(term,ZeroTerm):
            result = self.zero
        elif isinstance(term,ScaledPerturbationTheoryTerm):
            result = term.prefactor*self(term.base)
        elif isinstance(term,PerturbationTheoryTermSum):
            result = sum((self(child) for child in term.terms),self.zero)
        elif isinstance(term,PerturbationTheoryTermProduct):
            result = self(term.gen2) @ self(term.gen1)
        elif isinstance(term,ShiftedEnergyBaseline):
            result = self(term.base).T
        elif isinstance(term,ConditionalOperatorTerm):
            if term.condition not in self.masks:
                self.masks[term.condition] = term.condition.matches(
                    self.states[None,:,:],self.states[:,None,:],self.frequencies)
            mask = self.masks[term.condition]
            selected = mask if term.retain else ~mask
            result = np.where(selected,self(term.term),0.)
            if term.kind == 'resolvent':
                gap = self.energies[None,:]-self.energies[:,None]
                result = result/np.where(selected,gap,1.)
        else:
            raise TypeError(type(term))
        self.memo[term] = result
        return result


class ConditionalProjectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.w,cls.force = np.array([.04]),.004
        q = Symbols('q').q
        cls.condition = PolynomialInvariantCondition(q*(q-1),[0])
        h0 = AnalyticPerturbationTheorySolver.from_order(4,internals=False,logger=False).hamiltonian_expansion[0]
        h1 = HamiltonianExpansionTerm([['x']],order=1,identities=[0],logger=False)
        cls.hs = [h0,h1]+[ZeroTerm(logger=False) for _ in range(3)]
        cls.solver = AnalyticPerturbationTheorySolver(cls.hs,logger=False,resonance_condition=cls.condition)
        cls.coefficients = [[np.diag(cls.w)/2,np.diag(cls.w)/2],[np.array([cls.force])]]+[[0.] for _ in range(3)]
        cls.evaluator = PerturbationTheoryEvaluator(cls.solver,cls.coefficients,freqs=cls.w)
        cls.position = (1/cls.force)*h1

    def element(self,term,initial,final,backend='indexed',evaluator=None,coefficients=None,freqs=None):
        evaluator = self.evaluator if evaluator is None else evaluator
        initial,final = np.atleast_1d(initial),np.atleast_1d(final)
        coefficients = evaluator._prep_expansions(coefficients)
        freqs = evaluator.freqs if freqs is None else freqs
        value = 0.
        for change,permutations in evaluator.get_diff_map([[initial,[final]]]).items():
            expression = PerturbationTheoryExpressionEvaluator(term,term.get_poly_terms(change),change)
            value += float(expression.evaluate(permutations,coefficients,freqs,
                evaluation_backend=backend,zero_cutoff=1e-30)[0][0])
        return value

    def test_occupation_dependent_polynomial_condition(self):
        self.assertTrue(self.condition.matches([0],[1]))
        self.assertFalse(self.condition.matches([1],[2]))
        selected = ConditionalOperatorTerm.lookup(self.position,self.condition)
        self.assertAlmostEqual(self.element(selected,0,1),1/np.sqrt(2))
        self.assertEqual(self.element(selected,1,2),0.)
        a,b = Symbols('a b').vars
        pair = PolynomialRelationshipCondition(a+b-1,[0]).symmetrized()
        self.assertTrue(pair.matches([0],[1]))
        self.assertFalse(pair.matches([1],[2]))
        self.assertTrue(EnergyWindowCondition(.05).matches([0],[1],[.04]))

    def test_products_shift_the_condition_reference(self):
        selected = ConditionalOperatorTerm.lookup(self.position,self.condition)
        for backend in ['legacy','indexed']:
            self.assertAlmostEqual(self.element(selected*self.position,0,2,backend),0.)
            self.assertAlmostEqual(self.element(self.position*selected,0,2,backend),1/np.sqrt(2))

    def test_disjoint_axes_shift_the_condition_reference(self):
        w = np.array([.04,.07]);a,b = Symbols('a b').vars
        condition = PolynomialInvariantCondition(a*b,[0,1])
        hs = [AnalyticPerturbationTheorySolver.from_order(2,internals=False,logger=False).hamiltonian_expansion[0]]
        hs += [HamiltonianExpansionTerm([['x']],order=n,identities=[0],logger=False) for n in [1,2]]
        evaluator = PerturbationTheoryEvaluator(AnalyticPerturbationTheorySolver(hs,logger=False),
            [[np.diag(w)/2,np.diag(w)/2],[np.array([.004,0.])],[np.array([0.,.003])]],freqs=w)
        selected = ConditionalOperatorTerm.lookup(hs[1],condition)
        self.assertAlmostEqual(self.element(selected*hs[2],[0,0],[1,1],evaluator=evaluator),0.)
        self.assertAlmostEqual(self.element(hs[2]*selected,[0,0],[1,1],evaluator=evaluator),.004*.003/2)

    def test_adjoint_transposes_directional_conditions(self):
        a,b = Symbols('a b').vars
        condition = PolynomialRelationshipCondition(a,[0])
        term = ConditionalOperatorTerm.lookup(self.position,condition)
        self.assertAlmostEqual(self.element(term,0,1),1/np.sqrt(2))
        self.assertEqual(self.element(term,1,0),0.)
        self.assertAlmostEqual(self.element(ShiftedEnergyBaseline.lookup(term),1,0),1/np.sqrt(2))
        self.assertEqual(self.element(ShiftedEnergyBaseline.lookup(term),0,1),0.)

    def test_mask_algebra_and_primitive_divisors(self):
        q = ConditionalOperatorTerm.lookup(self.position,self.condition,retain=False)
        contradiction = ConditionalOperatorTerm.lookup(q,self.condition)
        self.assertEqual(self.element(contradiction,1,2),0.)
        inverse = ConditionalOperatorTerm.lookup(self.position*self.position,self.condition,
                                                 kind='resolvent',retain=False)
        self.assertAlmostEqual(self.element(inverse,2,4),-np.sqrt(12)/(4*self.w[0]))
        bad = PTTensorCoeffProductSum({((1,0,0),):0})
        self.assertFalse(bad.permute([0]).terms)

    def test_immutable_self_contained_serialization(self):
        term = ConditionalOperatorTerm.lookup(self.position,self.condition,kind='resolvent',retain=False)
        polynomial = term.get_poly_terms((1,))
        eager = polynomial.mutate(poly_obj=polynomial.poly_obj.to_eager())
        for original in [polynomial,eager]:
            state = json.loads(json.dumps(original.to_state(),default=lambda value:value.tolist()))
            restored = SqrtChangePoly.from_state(state)
            for change,permutations in self.evaluator.get_diff_map([[[1],[[2]]]]).items():
                expression = PerturbationTheoryExpressionEvaluator(term,restored,change)
                value = expression.evaluate(permutations,self.evaluator._prep_expansions(None),self.w)[0][0]
                self.assertAlmostEqual(value,-1/self.w[0])
        clone = pickle.loads(pickle.dumps(self.condition))
        left = ConditionEndpoint(self.condition & self.condition,'mask',True,0)
        right = ConditionEndpoint(self.condition & clone,'mask',True,0)
        self.assertEqual(left,right)
        self.assertEqual(hash(left),hash(right))
        with self.assertRaises(AttributeError):
            self.condition.tolerance = .01
        with self.assertRaises(AttributeError):
            left.retain = False
        self.assertIn('PathMask',polynomial.poly_obj.to_eager().format_expr())

    def test_divisor_filter_keeps_conditions(self):
        term = ConditionalOperatorTerm.lookup(self.position,self.condition,kind='resolvent',retain=False)
        polynomial = term.get_poly_terms((1,)).filter_energies([[(0,1)]])
        for initial,expected in [(0,0.),(1,-1/self.w[0])]:
            for change,permutations in self.evaluator.get_diff_map([[[initial],[[initial+1]]]]).items():
                expression = PerturbationTheoryExpressionEvaluator(term,polynomial,change)
                actual = expression.evaluate(permutations,self.evaluator._prep_expansions(None),self.w)[0][0]
                self.assertAlmostEqual(actual,expected)

    def test_same_formula_reuses_runtime_coefficients(self):
        term = ConditionalOperatorTerm.lookup(self.hs[1],self.condition,kind='resolvent',retain=False)
        first = self.element(term,1,2)
        coefficients = [[np.diag(self.w)/2,np.diag(self.w)/2],[np.array([2*self.force])]]+[[0.] for _ in range(3)]
        self.assertAlmostEqual(self.element(term,1,2,coefficients=coefficients),2*first)
        self.assertAlmostEqual(self.element(term,1,2,coefficients=coefficients,freqs=self.w/2),4*first)

    def test_public_batch_and_coefficient_filtering(self):
        states = [[[0],[[0],[1],[2]]]]
        operators = [[np.array([scale]),np.zeros((1,1)),np.zeros((1,1,1))] for scale in [1.,2.,-3.]]
        values = self.evaluator.get_operator_corrections(operators,states,order=2,operator_type=1,check_single=False)
        for component,scale in enumerate([1.,2.,-3.]):
            for n in range(3):
                np.testing.assert_allclose(values[component].corrections[0][n],
                    scale*values[0].corrections[0][n],atol=2e-12)
        filtered = self.evaluator.get_operator_corrections(operators[:1],states,order=2,
            operator_type=1,check_single=False,disallowed_coefficients=[(1,0,0)])
        np.testing.assert_allclose(filtered[0].corrections[0][0],values[0].corrections[0][0])
        for n in [1,2]:
            np.testing.assert_allclose(filtered[0].corrections[0][n],0.,atol=2e-12)

    def test_exact_fermi_feedback_and_forbidden_denominators(self):
        w = np.array([.04,.02,.13]);cubic = np.zeros((3,)*3)
        cubic[0,1,1]=cubic[1,0,1]=cubic[1,1,0]=.0006/3
        cubic[1,1,1]=.0004
        condition = IntegerLatticeCondition([[1,-2,0]])
        runner = AnalyticVPTRunner([[np.diag(w)/2,np.diag(w)/2],[cubic,0],[np.zeros((3,)*4),0,0,0,0]],
            order=2,freqs=w,internals=False,logger=False,resonance_condition=condition)
        a,b,k = [1,0,0],[0,2,0],[0,1,0]
        with warnings.catch_warnings():
            warnings.simplefilter('error',SmallEnergyDenominatorWarning)
            value = runner.eval.get_full_wavefunction_corrections([[a,[k]]],order=2)
            self.assertAlmostEqual(-value.corrections[0][2][0],-.0016125)
            h = runner.eval.get_reexpressed_hamiltonian([a,b],order=2,only_degenerate_terms=True)
            self.assertAlmostEqual(h.corrections[0][1][0],.0003)

    def test_intermediate_normalization_and_cache_clear(self):
        stack = self.solver.get_conditional_stack().extend(2)
        for term in stack.intermediate_wavefunctions(2)[1:]:
            for initial in [0,1]:
                for final in [0,1]:
                    self.assertAlmostEqual(self.element(term,initial,final),0.,places=12)
        AnalyticPerturbationTheorySolver.clear_caches()
        self.assertFalse(self.solver._conditional_stacks)
        self.assertIsNot(stack,self.solver.get_conditional_stack())

    def test_sector_proof_uses_relations_not_resolved_states(self):
        cubic = HamiltonianExpansionTerm([['x']*3],order=1,identities=[0],logger=False)
        oh = IntegerLatticeCondition([[1,-1,0]])
        fermi = IntegerLatticeCondition([[1,-2,0]])
        self.assertEqual(operator_sector(cubic,oh),'complement')
        self.assertIsNone(operator_sector(cubic,fermi))

    def test_integer_lattice_preserves_congruences_and_rank(self):
        even = IntegerLatticeCondition([[2,-2,0]])
        self.assertFalse(even.matches([0,2,0],[1,1,0]))
        self.assertTrue(even.matches([0,2,0],[2,0,0]))
        rank_two = IntegerLatticeCondition([[2,0,0],[0,3,0]])
        self.assertTrue(rank_two.matches([0,0,1],[2,3,1]))
        self.assertFalse(rank_two.matches([0,0,1],[1,3,1]))
        self.assertFalse(rank_two.matches([0,0,1],[2,3,2]))

    def test_function_serialization_after_runtime_compilation(self):
        from McUtils.Zachary.Symbolic.ElementaryFunctions import SymPyFunction
        condition = PolynomialInvariantCondition(SymPyFunction('q*(q-1)'),[0])
        self.assertTrue(condition.matches([0],[1]))
        restored = pickle.loads(pickle.dumps(condition))
        self.assertEqual(restored,condition)
        self.assertFalse(restored.matches([1],[2]))
        for modes in [[-1],[0,0],[.5]]:
            with self.assertRaises(ValueError):
                PolynomialInvariantCondition(Symbols('q').q,modes)

    def test_mode_pruning_matches_materialized_evaluation(self):
        polynomial=make_product_polynomial([np.ones(1),np.ones(1)],steps=0,representation='path')
        leaf=PTTensorCoeffProductDAG.from_sum(PTTensorCoeffProductSum(
            {((0,0,0),(0,0,1)):polynomial}))
        dag=leaf.mul_along(leaf,([0,1],[0,1]),remainder=([],[]),
                           index_classes=([0,0],[0,0]))
        plan=PTTensorCoeffProductDAGEvaluationPlan(dag,max_mode_indices=1)
        self.assertFalse(list(plan.iter_terms()))
        self.assertGreater(plan.mode_terms_pruned,0)
        values=[]
        for modes in [1,2]:
            expression=PerturbationTheoryExpressionEvaluator(self.hs[1],SqrtChangePoly(dag,[],[]),(),logger=False)
            permutations=[[np.zeros(modes,dtype=int),np.arange(modes)[None,:]]]
            coefficients=[[np.arange(1,modes+1,dtype=float)]]
            materialized=expression.evaluate(permutations,coefficients,np.ones(modes),evaluation_mode='materialized')
            streamed=expression.evaluate(permutations,coefficients,np.ones(modes),evaluation_mode='dag')
            np.testing.assert_allclose(streamed,materialized,atol=1e-13)
            values.append(float(np.asarray(streamed).ravel()[0]))
        self.assertEqual(values[0],0.)
        self.assertGreater(values[1],0.)

    def test_deep_polynomial_cache_keys_do_not_recurse(self):
        left = PolyPath.from_coeffs([np.array([1.,2.])],steps=0)
        right = PolyPath.from_coeffs([np.array([1.,2.])],steps=0)
        for _ in range(1500):
            left = PolyPath((),node=('shift',(left,(0,))),ndim=1,order=(1,))
            right = PolyPath((),node=('shift',(right,(0,))),ndim=1,order=(1,))
        self.assertEqual(hash(left),hash(right))
        self.assertEqual(left,right)
        self.assertNotEqual(left,PolyPath((),node=('shift',(right,(1,))),ndim=1,order=(1,)))

    def test_fourth_order_polynomial_dipole_matches_matrix_oracle(self):
        basis=HarmonicOscillatorProductBasis(1)
        states=BasisStateSpace(basis,np.arange(10)[:,None])
        x=basis.representation('x',coeffs=np.ones(1),axes=[[0],[0]]).get_representation_matrix(
            states,states,zero_element_warning=False)
        x=x.asarray() if hasattr(x,'asarray') else np.asarray(x)
        bare={self.hs[0]:np.diag((np.arange(10)+.5)*self.w[0]),self.hs[1]:self.force*x}
        matrix=MatrixInterpreter(bare,states.excitations,self.w)
        stack=self.solver.get_conditional_stack().extend(4)
        term=stack.transform_operator([self.position],4)[4]
        np.testing.assert_allclose(self.element(term,0,1),matrix(term)[1,0],atol=2e-12)

    def test_eighth_order_unitary_stack_with_higher_H_terms(self):
        order=8;w=.04;forces=[0.,.012,.001,-.002,.0008]+[0.]*4
        basis=HarmonicOscillatorProductBasis(1)
        states=BasisStateSpace(basis,np.arange(order+6)[:,None])
        x=basis.representation('x',coeffs=np.ones(1),axes=[[0],[0]]).get_representation_matrix(
            states,states,zero_element_warning=False)
        x=x.asarray() if hasattr(x,'asarray') else np.asarray(x)
        hs=[self.hs[0]];bare={hs[0]:np.diag((np.arange(len(states))+.5)*w)}
        for n in range(1,order+1):
            term=HamiltonianExpansionTerm([['x']],order=n,identities=[0],logger=False) if forces[n] else ZeroTerm(logger=False)
            hs.append(term)
            if forces[n]:bare[term]=forces[n]*x
        solver=AnalyticPerturbationTheorySolver(hs,logger=False,resonance_condition=self.condition)
        stack=solver.get_conditional_stack().extend(order)
        matrix=MatrixInterpreter(bare,states.excitations,np.array([w]))
        u=[matrix(term) for term in stack.wavefunctions]
        h=[matrix(term) for term in stack.hamiltonian]
        moments=[matrix(term) for term in stack.transform_operator([(1/forces[1])*hs[1]],order)]
        for n in range(order+1):
            gram=sum((u[i].T@u[n-i] for i in range(n+1)),matrix.zero)
            direct_h=sum((u[i].T@bare.get(hs[j],matrix.zero)@u[n-i-j]
                          for i in range(n+1) for j in range(n-i+1)),matrix.zero)
            direct_m=sum((u[i].T@x@u[n-i] for i in range(n+1)),matrix.zero)
            np.testing.assert_allclose(gram,matrix.identity if n==0 else 0.,atol=2e-12)
            np.testing.assert_allclose(direct_h,h[n],atol=2e-12)
            np.testing.assert_allclose(direct_m,moments[n],atol=2e-12)
        energy_errors,moment_errors=[],[]
        for strength in [1.,.5]:
            energies,rotation=np.linalg.eigh(sum(value[:2,:2]*strength**n for n,value in enumerate(h)))
            force=sum(f*strength**n for n,f in enumerate(forces))
            energy_errors.append(np.max(abs(energies-((np.arange(2)+.5)*w-force**2/(2*w)))))
            mu=rotation.T@sum(value[:2,:2]*strength**n for n,value in enumerate(moments))@rotation
            moment_errors.append(max(np.max(abs(np.diag(mu)+force/w)),abs(abs(mu[0,1])-1/np.sqrt(2))))
        self.assertLess(energy_errors[1],energy_errors[0]/200)
        self.assertLess(moment_errors[1],moment_errors[0]/200)


if __name__=='__main__':
    unittest.main()
