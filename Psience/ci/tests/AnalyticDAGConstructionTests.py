import types
import unittest
import tempfile
from unittest import mock

import numpy as np
import Psience.VPT2.Analytic as Analytic


class AnalyticDAGConstructionTests(unittest.TestCase):

    def test_long_path_term_sums_remain_shallow_and_evaluable(self):
        primitive = Analytic.OperatorExpansionTerm(
            [['x']], order=0, identities=[0], logger=False
        )
        term = Analytic.PerturbationTheoryTermSum(*([primitive] * 1500))
        poly = term.get_core_poly((1,))
        self._check_long_path_sum(primitive, poly)

    def test_long_path_product_sums_remain_shallow_and_evaluable(self):
        primitive = Analytic.OperatorExpansionTerm(
            [['x']], order=0, identities=[0], logger=False
        )
        poly = primitive.get_poly_terms((1,))
        product = types.SimpleNamespace(
            changes={(1,): [((), (), [], [], [], [], [])] * 1500},
            polynomial_representation='path', gen1=primitive, gen2=primitive,
            get_poly_product_terms=lambda *args: poly
        )
        result = Analytic.PerturbationTheoryTermProduct.get_core_poly(product, (1,))
        self._check_long_path_sum(primitive, result)

    def test_prune_operators_maps_shared_nodes_once(self):
        # x -> x + x forty times: 41 distinct nodes but 2**40 root-to-leaf paths
        primitive = Analytic.OperatorExpansionTerm(
            [['x']], order=0, identities=[0], logger=False
        )
        dag = primitive.get_poly_terms((1,)).poly_obj
        self.assertIsInstance(dag, Analytic.PTTensorCoeffProductDAG)
        for _ in range(40):
            dag = Analytic.PTTensorCoeffProductDAG._from_node('add', (dag, dag), ndim=dag.ndim)
        calls = []
        original = Analytic.PTTensorCoeffProductDAG._from_node.__func__
        def counting(cls, *args, **kwargs):
            calls.append(args[0])
            return original(cls, *args, **kwargs)
        with mock.patch.object(Analytic.PTTensorCoeffProductDAG, '_from_node', classmethod(counting)):
            pruned = dag.prune_operators([(9, 9)])
        self.assertLessEqual(len(calls), 41)
        self.assertEqual(pruned.operator_keys, frozenset({(0, 0)}))

    def _check_long_path_sum(self, primitive, poly):
        # Inspect without recursively walking the very tree under test.
        depths = {}
        pending = [(poly.poly_obj, False)]
        while pending:
            node, ready = pending.pop()
            if id(node) in depths:
                continue
            children = [arg for arg in node._node[1:]
                        if isinstance(arg, Analytic.PTTensorCoeffProductDAG)]
            if ready:
                depths[id(node)] = 1 + max(
                    (depths[id(child)] for child in children), default=0
                )
            else:
                pending.append((node, True))
                pending.extend((child, False) for child in children)
        self.assertEqual(poly.poly_obj.operator_keys, frozenset({(0, 0)}))
        self.assertLess(depths[id(poly.poly_obj)], 30)
        # Pruning and streaming must remain safe too; key lookup alone is
        # insufficient to protect a linearly deep tensor DAG.
        pruned = poly.poly_obj.prune_operators([(9, 9)])
        self.assertEqual(pruned.operator_keys, poly.poly_obj.operator_keys)
        evaluator = Analytic.PerturbationTheoryExpressionEvaluator(
            primitive, poly, (1,)
        )
        result = evaluator.evaluate(
            [[np.array([0]), np.array([[0]])]],
            [[np.array([.1])]], np.array([.04]),
            evaluation_mode='dag', evaluation_backend='indexed',
            zero_cutoff=1e-30
        )[0][0]
        self.assertAlmostEqual(result, 1500 * .1 / np.sqrt(2), places=10)


class DAGResonancePruningTests(unittest.TestCase):
    @staticmethod
    def leaf(energy_keys, dimension=2):
        polynomial = Analytic.make_product_polynomial([np.ones(1)] * dimension,
                                               steps=0, representation='path')
        energy = Analytic.PTEnergyChangeProductSum({key: polynomial for key in energy_keys})
        coefficients = ((0, 0) + tuple(range(dimension)),)
        return Analytic.PTTensorCoeffProductDAG.from_sum(Analytic.PTTensorCoeffProductSum(
            {coefficients: energy}, ndim=dimension, canonicalize=False))

    @staticmethod
    def prune(expression, changes, fixed=1, permutations=None):
        if permutations is None:
            permutations = np.array([[0, 1]])
        return Analytic.PerturbationTheoryExpressionEvaluator._prune_resonant_dag(
            expression, [[np.array([1, 0]), permutations]], changes, fixed, 2)

    def test_resonant_factor_prunes_numerator_before_visiting_it(self):
        source = self.leaf([()], 2)
        factor = Analytic.PTEnergyChangeProductSum.monomial((0, 1, 0), 1)
        expression = source.mul_simple(factor)
        changes = [([0], [1])]
        with mock.patch.object(Analytic.PTTensorCoeffProductDAGEvaluationPlan, 'iter_terms',
                               side_effect=AssertionError('DAG expanded during pruning')), \
             mock.patch.object(Analytic.PTTensorCoeffProductDAG, 'to_eager',
                               side_effect=AssertionError('eager fallback during pruning')):
            pruned, key = self.prune(expression, changes)
        self.assertEqual(pruned, 0)
        self.assertTrue(key)
        self.assertIsNone(source._eager)

    def test_missing_metadata_raises_without_materialization(self):
        leaf = self.leaf([((0, 1, 0),)])
        missing = Analytic.PTTensorCoeffProductDAG(leaf._node)
        with mock.patch.object(Analytic.PTTensorCoeffProductDAG, 'to_eager',
                               side_effect=AssertionError('eager fallback')):
            with self.assertRaisesRegex(ValueError, 'axis metadata'):
                self.prune(missing, [([0], [1])])

    def test_partial_energy_sum_keeps_free_denominator(self):
        expression = self.leaf([((0, 1, 0),), ((0, 0, 1),)])
        pruned, key = self.prune(expression, [([0], [1])])
        self.assertTrue(key)
        remaining = next(iter(pruned.to_eager().terms.values()))
        self.assertEqual(set(remaining.terms), {((0, 0, 1),)})
        self.assertIsNone(expression._eager)

    def test_different_rows_are_not_globally_pruned(self):
        expression = self.leaf([((0, 1, 0),)])
        pruned, key = self.prune(expression, [([0], [1])],
                                permutations=np.array([[0, 1], [1, 0]]))
        self.assertIs(pruned, expression)
        self.assertFalse(key)

    def test_side_specific_and_flipped_changes(self):
        expression = self.leaf([((0, 1, 0),)])
        sides = Analytic.SideDegeneracyChanges([], [([0], [1])])
        self.assertIs(self.prune(expression, sides)[0], expression)
        self.assertEqual(self.prune(expression.flip_energy_terms(), sides)[0], 0)

    def test_permutation_and_energy_shift(self):
        expression = self.leaf([((0, 0, 1),)])
        self.assertEqual(self.prune(expression.permute((1, 0)), [([0], [1])])[0], 0)
        expression = self.leaf([((0, 2, 0),)])
        self.assertEqual(self.prune(expression.shift_energies((-1, 0)), [([0], [1])])[0], 0)

    def test_shared_dag_can_be_reused_with_different_changes(self):
        expression = self.leaf([((0, 1, 0),)])
        self.assertEqual(self.prune(expression, [([0], [1])])[0], 0)
        unpruned, key = self.prune(expression, [([1], [1])])
        self.assertIs(unpruned, expression)
        self.assertFalse(key)

    def test_identical_removals_share_a_signature_across_physical_modes(self):
        expression = self.leaf([((0, 1, 0),), ((0, 0, 1),)])
        _, first = self.prune(expression, [([0], [1])])
        _, second = self.prune(expression, [([1], [1])], permutations=np.array([[1, 0]]))
        self.assertEqual(first, second)

    def test_compiled_aliases_do_not_leak_pruning_between_requests(self):
        from Psience.VPT2.Evaluators import AnalyticExpressionStore, compiled_evaluation_context
        polynomial = Analytic.make_product_polynomial([np.ones(1)] * 2,
                                               steps=0, representation='path')
        eager = Analytic.PTTensorCoeffProductSum({((0, 0, 0, 0, 1, 1),):
            Analytic.PTEnergyChangeProductSum({((0, 1, 0),): polynomial,
                                       ((0, 0, 1),): polynomial})})
        dag = Analytic.PTTensorCoeffProductDAG.from_sum(eager)
        wrapper = Analytic.PerturbationTheoryExpressionEvaluator(
            None, Analytic.SqrtChangePoly(dag, [1], [0]), (1,), logger=False)
        reference = Analytic.PerturbationTheoryExpressionEvaluator(
            None, Analytic.SqrtChangePoly(eager, [1], [0]), (1,), logger=False)
        arguments = ([[np.array([1, 1]), np.array([[0, 1]])]],
                     [[np.ones((2,) * 4)]], np.array([.1, .17]))
        changes = [[([0], [1])], None, [([1], [1])], [([0], [1])]]
        expected = [reference.evaluate(*arguments, degenerate_changes=change,
                                       evaluation_mode='materialized') for change in changes]
        with tempfile.TemporaryDirectory() as folder:
            for location in [None, folder]:
                for _ in range(2 if location is not None else 1):
                    store = AnalyticExpressionStore(location)
                    with compiled_evaluation_context(store):
                        for change, value in zip(changes, expected):
                            actual = wrapper.evaluate(*arguments, degenerate_changes=change)
                            np.testing.assert_allclose(actual, value, rtol=1e-13, atol=1e-13)

    def test_product_fixed_contact_and_remainder_orbits(self):
        polynomial = Analytic.make_product_polynomial([np.ones(1)] * 2,
                                               steps=0, representation='path')
        bare = Analytic.PTTensorCoeffProductDAG.from_sum(Analytic.PTTensorCoeffProductSum(
            {((0, 0, 0, 0, 1, 1),): polynomial}))
        factor = Analytic.PTEnergyChangeProductSum.monomial((0, 1, 0), 1)
        left = bare.mul_simple(factor)
        product = left.mul_along(bare, ((0,), (0,)), remainder=((1,), (1,)),
                                 index_classes=((0, 0), (0, 0)))
        pruned, signature = self.prune(product, [([0], [1])])
        self.assertEqual(pruned, 0)
        self.assertTrue(signature)

    def test_uniform_remainder_orbit_is_fixed_but_one_member_is_not(self):
        polynomial = Analytic.make_product_polynomial([np.ones(1)] * 2,
                                               steps=0, representation='path')
        bare = Analytic.PTTensorCoeffProductDAG.from_sum(Analytic.PTTensorCoeffProductSum(
            {((0, 0, 0, 0, 1, 1),): polynomial}))
        for delta, changes, eliminated in [((1, 1), [([0, 1], [1, 1])], True),
                                            ((1, 0), [([0], [1])], False)]:
            left = bare.mul_simple(Analytic.PTEnergyChangeProductSum.monomial((0,) + delta, 1))
            product = left.mul_along(bare, ((), ()), remainder=((0, 1), ()),
                                     index_classes=((0, 0), (0, 0)))
            result, key = self.prune(product, changes, fixed=2)
            if eliminated:
                self.assertEqual(result, 0)
                self.assertTrue(key)
            else:
                self.assertIs(result, product)
                self.assertFalse(key)

    def test_compiled_all_zero_pruned_request_preserves_result_shape(self):
        from Psience.VPT2.Evaluators import AnalyticExpressionStore, compiled_evaluation_context
        polynomial = Analytic.make_product_polynomial([np.ones(1)] * 2,
                                               steps=0, representation='path')
        bare = Analytic.PTTensorCoeffProductDAG.from_sum(Analytic.PTTensorCoeffProductSum(
            {((0, 0, 0, 0, 1, 1),): polynomial}))
        dag = bare.mul_simple(Analytic.PTEnergyChangeProductSum.monomial((0, 1, 0), 1))
        wrapper = Analytic.PerturbationTheoryExpressionEvaluator(None,
            Analytic.SqrtChangePoly(dag, [1], [0]), (1,), logger=False)
        with compiled_evaluation_context(AnalyticExpressionStore()):
            value = wrapper.evaluate([[np.array([1, 1]), np.array([[0, 1]])]],
                [[np.ones((2,) * 4)]], np.array([.1, .17]),
                degenerate_changes=[([0], [1])])
        np.testing.assert_array_equal(value, [[0.]])

    def test_explicit_opt_out_skips_prepass_in_compiled_and_uncompiled_calls(self):
        from Psience.VPT2.Evaluators import AnalyticExpressionStore, compiled_evaluation_context
        expression = self.leaf([((0, 1, 0),), ((0, 0, 1),)])
        leaf = expression._node[1]
        expression = Analytic.PTTensorCoeffProductDAG.from_sum(leaf.mutate(
            {((0, 0, 0, 0, 1, 1),): next(iter(leaf.terms.values()))}, canonicalize=False))
        wrapper = Analytic.PerturbationTheoryExpressionEvaluator(
            None, Analytic.SqrtChangePoly(expression, [1], [0]), (1,), logger=False)
        arguments = ([[np.array([1, 1]), np.array([[0, 1]])]],
                     [[np.ones((2,) * 4)]], np.array([.1, .17]))
        changes = [([0], [1])]
        for store in [None, AnalyticExpressionStore()]:
            with compiled_evaluation_context(store, evaluation_options={'prune_fixed_resonances': False}), \
                 mock.patch.object(Analytic.PTTensorCoeffProductDAG, 'exclude_degenerate_changes',
                                   side_effect=AssertionError('disabled prepass was called')):
                actual = wrapper.evaluate(*arguments, degenerate_changes=changes)
            with compiled_evaluation_context(store, evaluation_options={'prune_fixed_resonances': True}):
                expected = wrapper.evaluate(*arguments, degenerate_changes=changes)
            np.testing.assert_allclose(actual, expected, rtol=1e-13, atol=1e-13)
            self.assertGreater(float(np.asarray(actual).max()), 0.)

    def test_opt_out_rejects_non_boolean_values(self):
        from Psience.VPT2.Evaluators import compiled_evaluation_context
        expression = self.leaf([((0, 1, 0),)])
        with compiled_evaluation_context(None, evaluation_options={'prune_fixed_resonances': 'off'}):
            with self.assertRaisesRegex(TypeError, 'must be boolean'):
                self.prune(expression, [([0], [1])])

    def test_uncompiled_runner_public_method_honors_opt_out(self):
        from Psience.VPT2 import AnalyticVPTRunner
        frequencies = np.array([.0201, .0217])
        cubic = np.zeros((2,) * 3)
        cubic[0, 1, 1] = .00008
        runner = AnalyticVPTRunner(
            [[np.diag(frequencies)/2, np.diag(frequencies)/2], [cubic, 0],
             [np.zeros((2,) * 4), 0, 0, 0, 0]], order=2, freqs=frequencies,
            internals=False, logger=False, compile=False,
            evaluation_options={'prune_fixed_resonances': False})
        states = runner.prep_states([[0, 0], [1, 0], [0, 1]],
                                    degeneracy_specs={'groups': [[[1, 0], [0, 1]]]})
        with mock.patch.object(Analytic.PTTensorCoeffProductDAG, 'exclude_degenerate_changes',
                               side_effect=AssertionError('disabled runner prepass')):
            result = runner.get_reexpressed_hamiltonian(states, order=2,
                                                        state_specific_degeneracies=True)
        self.assertIsNotNone(result)

    def test_streamed_pruning_matches_native_materialized_mask(self):
        polynomial = Analytic.make_product_polynomial([np.ones(1)] * 2,
                                               steps=0, representation='path')
        leaf = Analytic.PTTensorCoeffProductDAG.from_sum(Analytic.PTTensorCoeffProductSum(
            {((0, 0, 0, 0, 1, 1),): polynomial}))
        factor = Analytic.PTEnergyChangeProductSum({
            ((0, 1, 0),): polynomial, ((1, 0, 1),): polynomial})
        base = leaf.mul_simple(factor)
        expressions = [base, base.flip_energy_terms(), base.permute((1, 0)),
                       base.shift_energies((1, 0)),
                       base.mul_along(leaf, ((0,), (0,)), remainder=((1,), (1,)),
                                      index_classes=((0, 0), (0, 0)))]
        coefficients = [[np.arange(1, 82, dtype=float).reshape((3,) * 4) / 100]]
        frequencies = np.array([.071, .113, .193])
        for expression in expressions:
            eager = expression.to_eager()
            for permutations in [np.array([[0, 1, 2]]), np.array([[1, 0, 2]]),
                                 np.array([[0, 1, 2], [1, 0, 2]])]:
                for changes in [[([0], [1])], [([1], [1])],
                                Analytic.SideDegeneracyChanges([([0], [1])], [([1], [1])])]:
                    arguments = [[np.array([1, 2, 1]), permutations]], coefficients, frequencies
                    reference = Analytic.PerturbationTheoryExpressionEvaluator(
                        None, Analytic.SqrtChangePoly(eager, [1], [0]), (1,), logger=False)
                    candidate = Analytic.PerturbationTheoryExpressionEvaluator(
                        None, Analytic.SqrtChangePoly(expression, [1], [0]), (1,), logger=False)
                    expected = reference.evaluate(*arguments, degenerate_changes=changes,
                                                  evaluation_mode='materialized')
                    actual = candidate.evaluate(*arguments, degenerate_changes=changes,
                                                 evaluation_mode='dag')
                    np.testing.assert_allclose(actual, expected, rtol=2e-12, atol=2e-12)


import inspect,itertools,tempfile,unittest
from unittest import mock
import numpy as np
import Psience.VPT2.Analytic as A
from Psience.VPT2 import AnalyticVPTRunner


class SymbolicBlockRecurrenceTests(unittest.TestCase):
    groups={'groups':[[[1,0,0],[0,2,0]]]}
    states=[[0,0,0],[1,0,0],[0,2,0]]

    def runner(self,block=True,compile=True,scale=1.,freq_scale=1.,folder=None,order=2,
               degeneracy_handling=None):
        w=np.array([.041,.0201,.0313])*freq_scale
        v3=np.zeros((3,)*3);v4=np.zeros((3,)*4)
        for indices,value in [((0,0,0),.0003),((0,1,1),.0002),((0,0,1),.0001)]:
            for p in set(itertools.permutations(indices)):v3[p]=value*scale
        # Keep only cubic perturbations, so order-two amplitudes are quadratic.
        dipole=[[0.,np.array([.2,.1,.3]),np.zeros((3,3)),np.zeros((3,)*3)]]
        dipole += [[0.,np.zeros(3),np.zeros((3,3)),np.zeros((3,)*3)] for _ in range(2)]
        expansions=[[np.diag(w)/2,np.diag(w)/2],[v3,0],[v4,0,0,0,0]]
        if order>2:
            expansions += [[0]*5 for _ in range(order-2)]
            for axis in dipole:axis += [np.zeros((3,)*(n+1)) for n in range(3,order+1)]
        return AnalyticVPTRunner(expansions,order=order,freqs=w,internals=False,
            dipole_expansion=dipole,logger=False,compile=compile,
            expressions_folder=folder,block_recurrence=block,
            degeneracy_handling=('conditional' if block else 'projected')
            if degeneracy_handling is None else degeneracy_handling)

    def values(self,runner,groups=None,order=2):
        states=runner.prep_states(self.states,degeneracy_specs=self.groups if groups is None else groups)
        return (runner.get_energy_corrections(states,order=order),
                runner.get_transition_moment_corrections(states,order=order),
                runner.get_reexpressed_hamiltonian(states,order=order))

    def flat(self,value):
        if isinstance(value,A.BasicAPTCorrections):
            yield value.initial_states.excitations
            for states in value.final_states:yield states.excitations
            yield from self.flat(value.corrections)
        elif isinstance(value,(list,tuple)):
            for part in value:yield from self.flat(part)
        elif value is not None:yield np.asarray(value)

    def assertValues(self,left,right):
        a,b=list(self.flat(left)),list(self.flat(right))
        self.assertEqual(len(a),len(b))
        for x,y in zip(a,b):np.testing.assert_allclose(x,y,atol=1e-11,rtol=1e-10)

    def test_default_flag_and_exact_runtime_relation(self):
        self.assertIs(inspect.signature(AnalyticVPTRunner.run_simple).parameters['block_recurrence'].default,True)
        self.assertEqual(A.PerturbationTheoryEvaluator.resolve_degeneracy_handling(None),'projected')
        self.assertEqual(A.PerturbationTheoryEvaluator.resolve_degeneracy_handling('conditional'),'conditional')
        condition=A.ResonanceBlockCondition()
        self.assertEqual(vars(condition),{})
        groups=(((1,0),(0,2)),)
        p=np.array([[0,1]])
        self.assertTrue(condition.evaluate_bound(np.array([0,0]),p,np.array([1,0]),np.array([0,2]),None,(groups,))[0])
        self.assertFalse(condition.evaluate_bound(np.array([1,0]),p,np.array([1,0]),np.array([0,2]),None,(groups,))[0])
        with self.assertRaisesRegex(ValueError,'bindings'):
            condition.evaluate_bound(np.array([0,0]),p,np.array([1,0]),np.array([0,2]),None)

    def test_feedback_is_a_shared_symbolic_subgraph(self):
        runner=self.runner()
        stack=runner.eval.solver.get_conditional_stack(A.ResonanceBlockCondition(),block_recurrence=True).extend(2)
        seen=set();pending=[stack.outside[2]]
        while pending:
            term=pending.pop()
            if id(term) in seen:continue
            seen.add(id(term))
            self.assertNotIsInstance(term,(A._BlockScalarTerm,A.ProjectedTerm))
            for name in ['term','base','subterm','gen1','gen2','terms']:
                value=getattr(term,name,())
                pending.extend([value] if isinstance(value,A.PerturbationTheoryTerm) else
                               [x for x in value if isinstance(x,A.PerturbationTheoryTerm)] if isinstance(value,(list,tuple)) else [])
        self.assertIn(id(stack.hamiltonian[1]),seen)
        poly=stack.outside[2].get_poly_terms((1,-1))
        self.assertIsInstance(poly.poly_obj,A.PTTensorCoeffProductDAG)
        self.assertFalse(any(hasattr(stack,name) for name in ['states','freqs','energies','expansions']))

    def test_compiled_and_uncompiled_match_canonical_blocks_without_numeric_F(self):
        expected=self.values(self.runner(False))
        for compile in [True,False]:
            with mock.patch.object(A.AnalyticBlockExpansion,'_scalar',side_effect=AssertionError('detached numeric F')):
                actual=self.values(self.runner(compile=compile))
            self.assertValues(actual,expected)

    def test_explicit_opt_out_retains_legacy_dispatch(self):
        runner=self.runner(False)
        with mock.patch.object(runner.eval.solver,'get_conditional_stack',side_effect=AssertionError('block path used')):
            self.assertIsInstance(runner.eval._projected_corrections([self.groups['groups'][0]],None),
                                  A.ProjectedCorrectionEvaluator)
        with self.assertRaisesRegex(TypeError,'boolean'):
            AnalyticVPTRunner(runner.eval,order=2,block_recurrence='yes')

    def test_conditional_mode_selects_requested_recurrence(self):
        for block in [True,False]:
            runner=self.runner(block,degeneracy_handling='conditional')
            evaluator=runner.eval._projected_corrections([self.groups['groups'][0]],None)
            self.assertIsInstance(evaluator,A.ConditionalCorrectionEvaluator)
            self.assertIs(evaluator.stack.block_recurrence,block)

    def test_projected_keeps_native_dispatch_with_either_block_flag(self):
        for compile in [True,False]:
            expected=self.values(self.runner(False,compile=compile))
            runner=self.runner(True,compile=compile,degeneracy_handling='projected')
            with mock.patch.object(runner.eval.solver,'get_conditional_stack',
                                   side_effect=AssertionError('projected rerouted to conditional')):
                self.assertIsInstance(runner.eval._projected_corrections([self.groups['groups'][0]],None),
                                      A.ProjectedCorrectionEvaluator)
                actual=self.values(runner)
            self.assertValues(actual,expected)

    def test_filtered_keeps_its_recurrences_with_either_block_flag(self):
        for compile in [True,False]:
            expected=self.values(self.runner(False,compile=compile,degeneracy_handling='filtered'))
            runner=self.runner(True,compile=compile,degeneracy_handling='filtered')
            with mock.patch.object(runner.eval.solver,'get_conditional_stack',
                                   side_effect=AssertionError('filtered rerouted to conditional')):
                self.assertIsNone(runner.eval._projected_corrections([self.groups['groups'][0]],None))
                actual=self.values(runner)
            self.assertValues(actual,expected)

    def test_cloning_preserves_explicit_opt_out(self):
        runner=self.runner(False)
        self.assertIs(runner.eval.modify_hamiltonian(None).block_recurrence,False)
        self.assertIs(runner.modify_hamiltonian(None).eval.block_recurrence,False)
        self.assertIs(inspect.signature(AnalyticVPTRunner.from_file).parameters['block_recurrence'].default,True)

    def test_rebinding_parameters_preserves_both_feedback_dependencies(self):
        runner=self.runner()
        space=runner.prep_states(self.states,degeneracy_specs=self.groups)
        first=runner.get_transition_moment_corrections(space,order=2)
        runner.eval.expansions[1][0] *= 2
        second=runner.get_transition_moment_corrections(space,order=2)
        reference=self.runner(False,scale=2.)
        self.assertValues(second,reference.get_transition_moment_corrections(
            reference.prep_states(self.states,degeneracy_specs=self.groups),order=2))
        a,b=np.asarray(first[0][0][2]),np.asarray(second[0][0][2])
        self.assertGreater(np.max(np.abs(a)),1e-8)
        np.testing.assert_allclose(b,4*a,atol=1e-12,rtol=1e-10)

    def test_runtime_blocks_are_not_frozen_in_memory_or_disk_compilation(self):
        alternate={'groups':[[[1,0,0],[0,1,0]]]}
        with tempfile.TemporaryDirectory() as folder:
            for location in [None,folder]:
                runner=self.runner(folder=location)
                for groups in [self.groups,alternate,self.groups]:
                    self.assertValues(self.values(runner,groups),self.values(self.runner(False),groups))

    def test_frequency_parameters_are_not_captured_by_the_stack(self):
        with tempfile.TemporaryDirectory() as folder:
            for scale in [1.,1.03,1.]:
                self.assertValues(self.values(self.runner(freq_scale=scale,folder=folder)),
                                  self.values(self.runner(False,freq_scale=scale)))

    def test_intermediate_normalization_matches_native_block_evaluator(self):
        results=[]
        for block in [False,True]:
            runner=self.runner(block)
            states=runner.prep_states(self.states,degeneracy_specs=self.groups)
            results.append(runner.get_reexpressed_hamiltonian(states,order=2,normalization='intermediate'))
        self.assertValues(*results)

    def test_third_order_feedback_and_normalization_match_native_blocks(self):
        # A third-order outside-block column exercises F2 feedback directly.
        # Keep a selected correction request, rather than building a full
        # third-order dipole surface in an otherwise second-order CI module.
        values=[]
        for block in [False,True]:
            runner=self.runner(block,order=3)
            space=runner.prep_states([[[[1,0,0]],[[0,1,0],[2,1,0]]]],
                                     degeneracy_specs=self.groups)
            values.append(runner.get_full_wavefunction_corrections(space,order=3))
        self.assertValues(*values)
        self.assertGreater(max(np.max(np.abs(row[3])) for row in values[1].corrections),1e-10)


    def test_conditional_inference_uses_existing_group_and_polyad_specs(self):
        specs = [
            dict(self.groups, group_filter='unfiltered'),
            {'polyads': self.groups['groups'], 'group_filter': 'unfiltered'},
            {'nT': [2, 1, 10], 'group_filter': 'unfiltered'},
            [dict(self.groups, group_filter='unfiltered'),
             {'polyads': self.groups['groups'], 'group_filter': 'unfiltered'}]
        ]
        for spec in specs:
            for compile in [True, False]:
                runner = self.runner(compile=compile)
                self.assertIsNone(runner.eval.solver.resonance_condition)
                with mock.patch.object(A.AnalyticBlockExpansion, '_scalar',
                                       side_effect=AssertionError('numerical projected fallback')):
                    actual = self.values(runner, groups=spec)
                self.assertValues(actual, self.values(self.runner(False), groups=spec))

    def test_conditional_no_resonances_do_not_reuse_previous_blocks(self):
        runner = self.runner()
        for spec in [self.groups, None, False, self.groups]:
            results = []
            for current in [runner, self.runner(False)]:
                space = current.prep_states(self.states, degeneracy_specs=spec)
                results.append((current.get_energy_corrections(space, order=2),
                                current.get_transition_moment_corrections(space, order=2)))
                if current is runner:
                    evaluator = current.eval._projected_corrections(
                        space.flat_space.degenerate_states, None)
                    self.assertIsInstance(evaluator, A.ConditionalCorrectionEvaluator)
                    self.assertIsInstance(evaluator.stack.condition, A.ResonanceBlockCondition)
                    self.assertEqual(vars(evaluator.stack.condition), {})
            self.assertValues(*results)

    def test_conditional_inference_preserves_existing_group_filter(self):
        spec = {'groups': [[[1, 0, 0], [0, 2, 0], [0, 0, 1]]]}
        runner = self.runner()
        space = runner.prep_states(self.states, degeneracy_specs=spec)
        evaluator = runner.eval._projected_corrections(space.flat_space.degenerate_states, None)
        resolved = A.ProjectedCorrectionEvaluator.resolve_blocks(space.flat_space.degenerate_states)
        self.assertEqual(evaluator.projection_bindings['projectors'],
                         (tuple(sorted(set(resolved.values()))),))
        self.assertValues(self.values(runner, groups=spec),
                          self.values(self.runner(False), groups=spec))

    def test_runner_resonances_have_no_independent_condition_input(self):
        for constructor in [AnalyticVPTRunner, AnalyticVPTRunner.construct,
                            AnalyticVPTRunner.from_hamiltonian, AnalyticVPTRunner.from_file]:
            self.assertNotIn('resonance_condition', inspect.signature(constructor).parameters)
        condition = A.SameStateCondition()
        with self.assertRaises(TypeError):
            AnalyticVPTRunner(self.runner().eval, order=2, resonance_condition=condition)
        for constructor, source in [(AnalyticVPTRunner.construct, 'not-a-molecule-file'),
                                    (AnalyticVPTRunner.from_file, 'not-a-molecule-file'),
                                    (AnalyticVPTRunner.from_hamiltonian, None)]:
            with self.assertRaisesRegex(TypeError, 'degeneracy_specs'):
                constructor(source, order=2, resonance_condition=condition)
        evaluator = self.runner(False, compile=False).eval
        evaluator.solver.resonance_condition = condition
        for mode in ['filtered', 'projected', 'conditional']:
            with self.assertRaisesRegex(TypeError, 'degeneracy_specs'):
                AnalyticVPTRunner(evaluator, order=2, degeneracy_handling=mode)
            self.assertIs(evaluator.block_recurrence, False)
            self.assertIs(evaluator.compile_expressions, False)

    def test_conditional_evaluator_infers_its_own_stack_from_resolved_groups(self):
        signature = inspect.signature(A.ConditionalCorrectionEvaluator).parameters
        self.assertNotIn('stack', signature)
        self.assertNotIn('resonance_groups', signature)
        self.assertIn('blocks', signature)
        self.assertIn('pairs', signature)
        runner = self.runner()
        evaluator = runner.eval._projected_corrections(self.groups['groups'], None)
        self.assertIsInstance(evaluator.stack.condition, A.ResonanceBlockCondition)
        self.assertEqual(evaluator.projection_bindings['projectors'],
                         ((((0, 2, 0), (1, 0, 0)),),))


if __name__ == '__main__':
    unittest.main()
