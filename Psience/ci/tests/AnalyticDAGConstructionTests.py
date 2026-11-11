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
        for block in [False,True]:
            for compile in [True,False]:
                expected=self.values(self.runner(False,compile=compile,
                    degeneracy_handling='projected' if block else 'filtered'))
                runner=self.runner(block,compile=compile,degeneracy_handling='filtered')
                with mock.patch.object(runner.eval.solver,'get_conditional_stack',
                                       side_effect=AssertionError('filtered rerouted to conditional')):
                    self.assertIsNone(runner.eval._projected_corrections([self.groups['groups'][0]],None))
                    actual=self.values(runner)
                self.assertIs(runner.eval.solver.filtered_block_recurrence,block)
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

class FilteredBlockRecurrenceTests(unittest.TestCase):
    # The fixture is deliberately shared with the existing symbolic-block tests.
    # Comparisons use the independent numerical projected recurrence.
    groups=SymbolicBlockRecurrenceTests.groups
    states=SymbolicBlockRecurrenceTests.states
    flat=SymbolicBlockRecurrenceTests.flat
    values=SymbolicBlockRecurrenceTests.values
    assertValues=SymbolicBlockRecurrenceTests.assertValues

    def runner(self, block=True, **options):
        options.setdefault('degeneracy_handling', 'filtered' if block else 'projected')
        return SymbolicBlockRecurrenceTests.runner(self, block, **options)

    def test_native_filtered_recurrence_and_normalization_match_matrix_blocks(self):
        for groups in [self.groups, {'groups':[[[1,0,0],[0,1,0]]]},
                       {'groups':[[[1,0,0],[0,1,0],[0,2,0]]]}]:
            for compile in [True,False]:
                expected=self.values(self.runner(False,compile=compile),groups)
                with mock.patch.object(A.ConditionalTransformationStack,'extend',side_effect=AssertionError('Conditional used')), \
                     mock.patch.object(A.ProjectedCorrectionEvaluator,'__init__',side_effect=AssertionError('projected used')), \
                     mock.patch.object(A.AnalyticBlockExpansion,'_scalar',side_effect=AssertionError('numeric F used')):
                    actual=self.values(self.runner(compile=compile),groups)
                self.assertValues(actual,expected)

    def test_energy_filters_preserve_integer_domains_and_only_filter_divisors(self):
        term=A.HamiltonianExpansionTerm([['x']],order=1,identities=[0],logger=False)
        for base in [term.get_poly_terms((1,)), A.PerturbationOperator.lookup(term).get_poly_terms((1,))]:
            scoped=A._change_domain_polynomial(base,(1,),True)
            for mode, terms in [('include',[(0,1)]), ('exclude',[(0,1)]), ('match',[((0,1),)])]:
                bare=list(base.filter_energies(terms,mode=mode).change_equations())
                selected=list(scoped.filter_energies(terms,mode=mode).change_equations())
                self.assertEqual(len(bare),len(selected))
                for (_,expected),(domains,actual) in zip(bare,selected):
                    self.assertEqual(domains,((True,(0,),(1,)),))
                    self.assertEqual(expected.format_expr(),actual.format_expr())

    def test_F_is_shared_and_exports_ordinary_change_equations(self):
        runner=self.runner();runner.eval._projected_corrections(self.groups['groups'],None)
        solver=runner.eval.solver
        f=solver.block_hamiltonian_correction(1)
        self.assertIs(f,solver.block_hamiltonian_correction(1))
        seen=set();pending=[solver.wavefunction_correction(2)]
        while pending:
            term=pending.pop()
            if term in seen:continue
            seen.add(term)
            self.assertNotIsInstance(term,(A.ConditionalOperatorTerm,A._BlockScalarTerm,A.ProjectedTerm))
            if isinstance(term,A.PerturbationTheoryTermProduct):pending += [term.gen1,term.gen2]
            elif isinstance(term,(A.ScaledPerturbationTheoryTerm,A.ShiftedEnergyBaseline)):pending += [term.base]
            elif isinstance(term,A.PerturbationOperator):pending += [term.subterm]
            elif isinstance(term,A.OperatorExpansionTerm):pass
            else:pending += list(term.expressions)
        self.assertIn(f,seen)
        poly=solver.wavefunction_correction(2).get_poly_terms((1,-1))
        equations=list(poly.change_equations())
        self.assertTrue(equations)
        self.assertTrue(any(domains for domains,_ in equations))
        for domains,equation in equations:
            for retain,origin,change in domains:
                self.assertIsInstance(retain,bool)
                self.assertTrue(all(type(v) in (int,np.int64,np.int32) for v in origin+change))
            for value in equation.poly_obj.terms.values():
                self.assertIsInstance(value,A.PTEnergyChangeProductSum)
                for path in value.terms:
                    self.assertTrue(all(type(row[0]) is int for row in path))
        self.assertNotIn('PathMask',poly.format_expr())
        self.assertFalse(any(hasattr(f,name) for name in ['states','freqs','energies','expansions']))

    def test_pruning_rebinding_and_parameter_dependence(self):
        alternate={'groups':[[[1,0,0],[0,1,0]]]}
        with tempfile.TemporaryDirectory() as folder:
            runner=self.runner(folder=folder)
            for groups in [self.groups,alternate,self.groups]:
                self.assertValues(self.values(runner,groups),self.values(self.runner(False),groups))
            runner.eval.expansions[1][0] *= 2
            self.assertValues(self.values(runner),self.values(self.runner(False,scale=2)))
            for scale in [1.03,1.]:
                self.assertValues(self.values(self.runner(freq_scale=scale,folder=folder)),
                                  self.values(self.runner(False,freq_scale=scale)))

    def test_pruning_opt_out_preserves_results_and_cache_reuse(self):
        runner=self.runner()
        expected=self.values(self.runner(False))
        for enabled in [True,False,True]:
            runner.eval.evaluation_options['prune_fixed_resonances']=enabled
            self.assertValues(self.values(runner),expected)

    def test_legacy_and_block_policies_do_not_alias_in_one_disk_cache(self):
        from Psience.VPT2.Evaluators import _fingerprint
        with tempfile.TemporaryDirectory() as folder:
            hashes={}
            for block in [False,True,False,True]:
                runner=self.runner(block,folder=folder,degeneracy_handling='filtered')
                runner.eval.evaluation_options['prune_fixed_resonances']=False
                actual=self.values(runner)
                reference=self.runner(False,degeneracy_handling='projected' if block else 'filtered')
                reference.eval.evaluation_options['prune_fixed_resonances']=False
                self.assertValues(actual,self.values(reference))
                hashes[block]=_fingerprint(runner.eval.solver)
            self.assertNotEqual(hashes[False],hashes[True])

    def test_no_resonance_and_opt_out_preserve_legacy_filtered(self):
        for groups in [None,[],[[[1,0,0]]]]:
            a=self.runner();b=self.runner(False,degeneracy_handling='filtered')
            args=dict(degeneracy_specs=None if groups is None else {'groups':groups})
            sa=a.prep_states(self.states,**args);sb=b.prep_states(self.states,**args)
            self.assertValues(a.get_energy_corrections(sa),b.get_energy_corrections(sb))
            self.assertFalse(a.eval.solver.filtered_block_recurrence)
        runner=self.runner(False,degeneracy_handling='filtered')
        self.values(runner)
        self.assertFalse(runner.eval.solver.filtered_block_recurrence)
        for mode in ['conditional','projected']:
            runner=self.runner(degeneracy_handling=mode)
            selected=runner.eval._projected_corrections(self.groups['groups'],None)
            self.assertIsInstance(selected,A.ConditionalCorrectionEvaluator if mode=='conditional' else A.ProjectedCorrectionEvaluator)
            self.assertFalse(runner.eval.solver.filtered_block_recurrence)

    def test_intermediate_H_blocks_and_third_order_feedback(self):
        values=[]
        for block in [False,True]:
            runner=self.runner(block)
            space=runner.prep_states(self.states,degeneracy_specs=self.groups)
            values.append(runner.get_reexpressed_hamiltonian(space,order=2,normalization='intermediate'))
        self.assertValues(*values)
        values=[]
        for block in [False,True]:
            runner=self.runner(block,order=3)
            space=runner.prep_states([[[[1,0,0]],[[0,1,0],[2,1,0]]]],degeneracy_specs=self.groups)
            values.append(runner.get_full_wavefunction_corrections(space,order=3))
        self.assertValues(*values)

    def test_change_domain_prepass_skips_rejected_subtree_before_expansion(self):
        poly=A.make_product_polynomial([np.ones(1)],steps=0,representation='path')
        leaf=A.PTTensorCoeffProductDAG.from_sum(A.PTTensorCoeffProductSum({((1,0,0),):poly}))
        dag=leaf.with_change_domain((1,),True)
        rows=[[np.array([2]),np.array([[0]])]]
        groups=(((0,),(1,)),)
        stats={}
        with mock.patch.object(A.PTTensorCoeffProductDAGEvaluationPlan,'iter_terms',side_effect=AssertionError('expanded before pruning')):
            selected,digest=dag.select_change_domains(rows,(groups,),1,1,statistics=stats)
        self.assertEqual(selected,0);self.assertTrue(digest)
        self.assertEqual(stats['change_domains_pruned'],1)
        self.assertIs(dag._node[1],leaf)
        selected,_=dag.select_change_domains([[np.array([0]),np.array([[0]])]],(groups,),1,1)
        self.assertIs(selected,leaf)

    def test_fourth_order_polynomial_recurrence_matches_finite_basis(self):
        from Psience.BasisReps import HarmonicOscillatorProductBasis,BasisStateSpace
        order=4;w=np.array([.04]);force=.004
        hs=[A.AnalyticPerturbationTheorySolver.from_order(order,internals=False,logger=False).hamiltonian_expansion[0],
            A.HamiltonianExpansionTerm([['x']],order=1,identities=[0],logger=False)]
        hs += [A.ZeroTerm(logger=False) for _ in range(order-1)]
        solver=A.AnalyticPerturbationTheorySolver(hs,logger=False)
        solver.filtered_block_recurrence=True
        coefficients=[[np.diag(w)/2,np.diag(w)/2],[np.array([force])]]+[[0.] for _ in range(order-1)]
        evaluator=A.PerturbationTheoryEvaluator(solver,coefficients,freqs=w,compile=False)
        basis=HarmonicOscillatorProductBasis(1);states=BasisStateSpace(basis,np.arange(14)[:,None])
        x=basis.representation('x',coeffs=np.ones(1),axes=[[0],[0]]).get_representation_matrix(states,states,zero_element_warning=False)
        x=x.asarray() if hasattr(x,'asarray') else np.asarray(x)
        energies=(states.excitations[:,0]+.5)*w[0]
        P=np.eye(len(states),dtype=bool);P[:2,:2]=True
        gap=energies[:,None]-energies[None,:]
        Y=[np.eye(len(states))];X=[np.zeros_like(x)];N=[np.eye(len(states))];F=[np.diag(energies)]
        bindings=dict(projectors=((((0,),(1,)),),),scalar_count=0,scalar_order=len(coefficients))
        def element(term,a,b):
            total=0.
            for change,rows in evaluator.get_diff_map([[np.array([a]),[np.array([b])]]]).items():
                expression=A.PerturbationTheoryExpressionEvaluator(term,term.get_poly_terms(change),change,logger=False)
                value = float(np.asarray(expression.evaluate(rows,coefficients,w,projection_bindings=bindings,zero_cutoff=1e-30)).ravel()[0])
                exported = 0.
                if not A.nput.is_numeric(expression.expr):
                    for domains, equation in expression.expr.change_equations():
                        valid = True
                        for retain, origin, delta in domains:
                            left = a + (origin[0] if origin else 0)
                            right = left + (delta[0] if delta else 0)
                            inside = left == right or left in (0,1) and right in (0,1)
                            valid &= left >= 0 and right >= 0 and (inside if retain else not inside)
                        if valid:
                            bare = A.PerturbationTheoryExpressionEvaluator(term,equation,change,logger=False)
                            exported += float(np.asarray(bare.evaluate(rows,coefficients,w,zero_cutoff=1e-30)).ravel()[0])
                self.assertAlmostEqual(value,exported,places=12)
                total += value
            return total
        for n in range(1,order+1):
            norm=-.5*sum((Y[j].T@Y[n-j] for j in range(1,n)),np.zeros_like(x))
            N.append(np.where(P,norm,0))
            raw=force*x@Y[n-1]
            outside=raw-sum((X[n-r]@F[r] for r in range(1,n)),np.zeros_like(x))
            X.append(np.where(P,0,-outside/np.where(P,1,gap)))
            B=raw-sum((N[n-r]@F[r] for r in range(1,n)),np.zeros_like(x))
            F.append(np.where(P,(B+B.T)/2,0));Y.append(X[n]+N[n])
            for source,target in [(0,0),(0,1),(0,2),(1,0),(1,2)]:
                self.assertAlmostEqual(element(solver.full_wavefunction_correction(n),source,target),Y[n][target,source],places=12)
                self.assertAlmostEqual(element(solver.block_hamiltonian_correction(n),source,target),F[n][target,source],places=12)


import itertools
import math
import tempfile
import unittest
from unittest import mock

import numpy as np
import Psience.VPT2.Analytic as A
from Psience.VPT2 import AnalyticVPTRunner, PerturbationTheoryHamiltonian
from Psience.VPT2.Evaluators import _OccupationCoefficientProjection


class SymbolicActionExpansionTests(unittest.TestCase):
    @staticmethod
    def model(ndim=3, internals=False, scale=1., freq_scale=1., compile=True, folder=None):
        rng = np.random.default_rng(410)
        w = np.array([.013,.029,.047][:ndim])*freq_scale
        def symmetric(rank, magnitude):
            x = rng.normal(size=(ndim,)*rank)*magnitude*scale
            return sum(x.transpose(p) for p in itertools.permutations(range(rank)))/math.factorial(rank)
        v3,v4 = symmetric(3,.0003),symmetric(4,.00002)
        if internals:
            # G derivatives are symmetric within coordinate and momentum
            # slots, and need not have potential-tensor permutation symmetry.
            g1 = rng.normal(size=(ndim,)*3)*.0001*scale
            g1 = (g1+g1.swapaxes(-1,-2))/2
            g2 = rng.normal(size=(ndim,)*4)*.00001*scale
            g2 = (g2+g2.swapaxes(0,1))/2
            g2 = (g2+g2.swapaxes(-1,-2))/2
            coefficients = [[np.diag(w)/2,np.diag(w)/2],[v3,g1],[v4,g2,.000007*scale]]
        else:
            coefficients = [[np.diag(w)/2,np.diag(w)/2],[v3,0],[v4,0,0,0,0]]
        solver = A.AnalyticPerturbationTheorySolver.from_order(2,internals=internals,logger=False)
        return A.PerturbationTheoryEvaluator(solver,coefficients,w,compile=compile,expressions_folder=folder)

    def assertExpansion(self, left, right, atol=1e-12):
        for a,b in zip(left,right):
            np.testing.assert_allclose(a,b,rtol=1e-12,atol=atol)

    @staticmethod
    def native_energies(evaluator, states):
        # Oracle only: the implementation never samples these states.
        solver = A.AnalyticPerturbationTheorySolver(evaluator.solver.hamiltonian_expansion,logger=False)
        prepared = evaluator._prep_expansions(evaluator.expansions)
        identity = np.arange(len(evaluator.freqs))[None]
        total = np.zeros(len(states))
        for n in (0,2):
            term = solver.energy_correction(n)
            wrapper = A.PerturbationTheoryExpressionEvaluator(term,term.get_poly_terms(()),(),logger=False)
            total += np.asarray(wrapper.evaluate([(s,identity) for s in states],prepared,
                                                  evaluator.freqs,zero_cutoff=0)).reshape(-1)
        return total

    def test_polypath_all_operations_without_sampling_or_eager_expansion(self):
        p = A.PolyPath.from_coeffs([[2,3,4],[1,-2]],prefactor=3,steps=2)
        q = A.PolyPath.from_coeffs([[3,-1],[2,0,1]],steps=0)
        projected = A.PolyPath._from_node('project',(p,(1,),(0,)),1,(1,))
        linear = A.PolyPath._from_node('linear',((p,2.),(q,-3.)),2,(2,2))
        paths = [p,p+q,linear,p.scale(2.5),p.shift((2,-1)),p.permute((1,0)),
                 p.permutation_sum(((0,1),(1,0))),p.mul_simple(q),
                 p.mul_along(q,[(0,),(1,)]),projected,p.permute((2,0,1),allow_padding=True)]
        def oracle(path, powers):
            eager = path.to_eager()
            polynomials = eager.polys if isinstance(eager,A.ProductPTPolynomialSum) else [eager]
            value = 0.
            for poly in polynomials:
                weight = poly.prefactor*2**(-poly.steps/2)
                for i,power in enumerate(powers):
                    axis = poly.coeffs[i]
                    weight *= sum(axis[k]*math.comb(k,power)*(-.5)**(k-power)
                                  for k in range(power,len(axis)))
                value += weight
            return value
        for path in paths:
            for powers in itertools.product(range(3),repeat=path.ndim):
                expected = oracle(path,powers)
                with (mock.patch.object(A.PolyPath,'to_eager',side_effect=AssertionError('eager expansion')),
                      mock.patch.object(A.PolyPath,'evaluate_polynomial',side_effect=AssertionError('state sampling'))):
                    self.assertAlmostEqual(path.coefficient(powers,shift=-.5),expected,places=11)
        with self.assertRaises(ValueError):p.coefficient((-1,))
        with self.assertRaises(ValueError):p.coefficient((1.5,))
        self.assertEqual(p.coefficient((0,0,1)),0)

    def test_deep_and_shared_coefficient_requests_are_iterative(self):
        path = A.PolyPath.from_coeffs([[1,2,3]],steps=0)
        for _ in range(40):
            path = A.PolyPath._from_node('add',(path,path),1,(2,))
        for i in range(1800):
            path = path.shift((1 if i%2 == 0 else -1,))
        self.assertEqual(path.coefficient((2,),shift=-.5),3*2.**40)
        self.assertEqual(path.coefficient((1,),shift=-.5),-2.**40)

    def test_one_mode_closed_form_and_both_sources(self):
        w,a,b = .043,.0003,.00002
        solver = A.AnalyticPerturbationTheorySolver.from_order(2,internals=False,logger=False)
        coefficients = [[np.array([[w/2]]),np.array([[w/2]])],
                        [np.full((1,)*3,a),0],[np.full((1,)*4,b),0,0,0,0]]
        expected = [3*b/8-7*a*a/(16*w),np.array([w]),np.array([[1.5*b-15*a*a/(4*w)]])]
        for compiled in (False,True):
            evaluator = A.PerturbationTheoryEvaluator(solver,coefficients,np.array([w]),compile=compiled)
            for source in ('energy','hamiltonian'):
                expansion,kernels = evaluator.get_action_expansion(correction=source)
                self.assertExpansion(expansion,expected,atol=1e-15)
                self.assertEqual(set(kernels),{(n,p) for n in (0,2) for p in ((),(1,),(2,),(1,1))})

    def test_coupled_modes_reconstruct_energies_and_match_nielsen(self):
        evaluator = self.model()
        states = np.array(list(itertools.product(range(3),repeat=3)))
        expected = self.native_energies(evaluator,states)
        v3,v4 = evaluator.expansions[1][0],evaluator.expansions[2][0]
        nielsen = PerturbationTheoryHamiltonian._get_Nielsen_xmat(evaluator.freqs,6*v3,24*v4,None,None).sum(axis=0)
        reference = None
        for compiled in (False,True):
            evaluator.compile_expressions = compiled
            for source in ('energy','hamiltonian'):
                result,_ = evaluator.get_action_expansion(correction=source)
                C,linear,X = result
                I = states+.5
                actual = C+I@linear+np.einsum('si,ij,sj->s',I,X,I)
                np.testing.assert_allclose(actual,expected,atol=1e-13,rtol=1e-13)
                np.testing.assert_allclose(X,np.triu(nielsen),atol=1e-13,rtol=1e-12)
                if reference is None:reference = result
                self.assertExpansion(result,reference)

    def test_internal_kinetic_and_scalar_terms_preserve_constant(self):
        evaluator = self.model(internals=True)
        states = np.array(list(itertools.product(range(3),repeat=3)))
        expected = self.native_energies(evaluator,states)
        for source in ('energy','hamiltonian'):
            (C,linear,X),_ = evaluator.get_action_expansion(correction=source)
            I = states+.5
            np.testing.assert_allclose(C+I@linear+np.einsum('si,ij,sj->s',I,X,I),expected,atol=1e-13,rtol=1e-13)
            self.assertNotEqual(C,0)

    def test_lower_degree_second_order_terms_preserve_linear_actions(self):
        w = np.array([.023,.041]);shift = np.array([.0002,.0003]);constant = .0007
        base = A.AnalyticPerturbationTheorySolver.from_order(2,internals=False,logger=False)
        perturbation = A.HamiltonianExpansionTerm([['x','x'],[]],order=2,identities=[0,2],logger=False)
        solver = A.AnalyticPerturbationTheorySolver(
            [base.hamiltonian_expansion[0],A.ZeroTerm(logger=False),perturbation],logger=False)
        coefficients = [[np.diag(w)/2]*2,[0,0],[np.diag(shift),0,constant,0,0]]
        evaluator = A.PerturbationTheoryEvaluator(solver,coefficients,w)
        for source in ('energy','hamiltonian'):
            result,_ = evaluator.get_action_expansion(correction=source)
            self.assertExpansion(result,[constant,w+shift,np.zeros((2,2))],atol=1e-15)

    def test_cubic_occupation_coefficients_cancel_symbolically(self):
        for internals in (False,True):
            evaluator = self.model(internals=internals)
            for source in ('energy','hamiltonian'):
                kernels = evaluator.solver.action_expansion_terms(correction=source)
                raw = kernels[2,()].poly_obj._node[1]
                for powers in ((3,),(2,1),(1,2),(1,1,1)):
                    # Export explicitly collects identical force products and
                    # canonical frequency rows; no numerical tensors are bound.
                    expression = raw.project_coefficients(powers).to_eager()
                    for value in expression.terms.values():
                        rows = value.terms.values() if isinstance(value,A.PTEnergyChangeProductSum) else [value]
                        for numerator in rows:
                            self.assertAlmostEqual(A.PolyPath.from_polynomial(numerator).coefficient(()),0.,places=13)

    def test_resonance_groups_flags_and_conditional_stack_do_not_enter(self):
        evaluator = self.model()
        expected,kernels = evaluator.get_action_expansion()
        evaluator.solver.resonance_condition = A.ResonanceBlockCondition()
        evaluator.solver.filtered_block_recurrence = True
        for handling in ('conditional','filtered','projected'):
            evaluator.degeneracy_handling = handling
            evaluator._filtered_change_bindings = dict(projectors=((((1,0,0),(0,2,0)),),))
            with (mock.patch.object(A.ConditionalTransformationStack,'extend',side_effect=AssertionError('Conditional route')),
                 mock.patch.object(A.ProjectedCorrectionEvaluator,'__init__',side_effect=AssertionError('projected route')),
                 mock.patch.object(A.PolyPath,'evaluate_polynomial',side_effect=AssertionError('occupation sampling')),
                 mock.patch.object(PerturbationTheoryHamiltonian,'_invert_action_expansion_tensors',side_effect=AssertionError('state fit'))):
                # Only construction is under the polynomial-evaluation guard:
                # numerical contraction legitimately evaluates constant kernels.
                symbolic = evaluator.solver.action_expansion_terms()
                self.assertEqual(len(symbolic),8)
            for source in ('energy','hamiltonian'):
                actual,_ = evaluator.get_action_expansion(correction=source)
                self.assertExpansion(actual,expected)

    def test_disk_and_memory_cache_rebind_force_constants_and_frequencies(self):
        with tempfile.TemporaryDirectory() as folder:
            evaluator = self.model(ndim=2,folder=folder)
            first,_ = evaluator.get_action_expansion()
            stats = evaluator.compilation_statistics()
            evaluator.expansions[1][0] *= 1.7
            evaluator.expansions[2][0] *= .3
            evaluator.freqs *= 1.2
            evaluator.expansions[0] = [np.diag(evaluator.freqs)/2]*2
            changed,_ = evaluator.get_action_expansion()
            self.assertEqual(evaluator.compilation_statistics()['compilations'],stats['compilations'])
            self.assertFalse(np.allclose(first[2],changed[2]))
            evaluator.compile_expressions = False
            direct,_ = evaluator.get_action_expansion()
            self.assertExpansion(changed,direct)
            warm = self.model(ndim=2,folder=folder)
            warm.expansions,warm.freqs = evaluator.expansions,evaluator.freqs
            resumed,_ = warm.get_action_expansion()
            self.assertExpansion(resumed,changed)
            self.assertEqual(warm.compilation_statistics()['compilations'],0)

    def test_zero_coefficients_are_not_baked_into_symbolic_kernels(self):
        evaluator = self.model(ndim=2)
        saved = evaluator.expansions[1][0].copy()
        evaluator.expansions[1][0][:] = 0
        evaluator.get_action_expansion()
        compiled = evaluator.compilation_statistics()['compilations']
        evaluator.expansions[1][0][:] = saved
        actual,_ = evaluator.get_action_expansion()
        self.assertEqual(evaluator.compilation_statistics()['compilations'],compiled)
        fresh,_ = self.model(ndim=2).get_action_expansion()
        self.assertExpansion(actual,fresh)

    def test_symbolic_derivative_matches_parameter_rebinding(self):
        evaluator = self.model(ndim=1)
        kernels = evaluator.precompile_action_expansion()
        coefficient = kernels[2,(2,)]
        prepared = evaluator._prep_expansions(evaluator.expansions)
        w = evaluator.freqs[0]
        a = prepared[1][0][0,0,0]
        epsilon = 1e-7
        values = []
        for delta in (-epsilon,epsilon):
            trial = [[x.copy() if hasattr(x,'copy') else x for x in row] for row in prepared]
            trial[1][0][0,0,0] = a+delta
            values.append(coefficient.evaluate_coefficient([trial],np.array([w]),indices=(0,))[0])
        derivative = (values[1]-values[0])/(2*epsilon)
        self.assertAlmostEqual(derivative,-15*a/(2*w),places=10)
        values = [coefficient.evaluate_coefficient([prepared],np.array([w+delta]),indices=(0,))[0]
                  for delta in (-epsilon,epsilon)]
        self.assertAlmostEqual((values[1]-values[0])/(2*epsilon),15*a*a/(4*w*w),places=10)

    def test_interfaces_and_reject_unsupported_or_bound_definitions(self):
        evaluator = self.model(ndim=2)
        runner = AnalyticVPTRunner(evaluator,order=2,logger=False,dipole_expansion=[])
        result,_ = runner.get_action_expansion()
        reference,_ = evaluator.get_action_expansion()
        self.assertExpansion(result,reference)
        with self.assertRaises(ValueError):evaluator.get_action_expansion(order=-1)
        with self.assertRaises(ValueError):evaluator.get_action_expansion(correction='wavefunction')
        with self.assertRaises(ValueError):evaluator.get_action_expansion(degenerate_changes=[(1,-2)])
        compiled = evaluator.precompile_action_expansion()
        prepared = evaluator._prep_expansions(evaluator.expansions)
        for indices in ((.5,),(2,),(-1,),()):
            with self.assertRaises(ValueError):
                compiled[2,(2,)].evaluate_coefficient([prepared],evaluator.freqs,indices=indices)
        # compile=False must also hold inside an outer compiled context.
        from Psience.VPT2.Evaluators import compiled_evaluation_context
        evaluator.compile_expressions = False
        with compiled_evaluation_context(evaluator.expression_store):
            with mock.patch.object(evaluator.expression_store,'compile',side_effect=AssertionError('opt-out ignored')):
                uncompiled,_ = evaluator.get_action_expansion()
        self.assertExpansion(uncompiled,reference)
        evaluator.solver = A.AnalyticPerturbationTheorySolver.from_order(
            2,internals=False,logger=False,polynomial_representation='eager')
        eager,_ = evaluator.get_action_expansion()
        self.assertExpansion(eager,reference)
        evaluator.solver.allowed_energy_changes = [(1,)]
        with self.assertRaises(ValueError):evaluator.get_action_expansion()
        evaluator = self.model(ndim=2)
        evaluator.solver.hamiltonian_expansion[1] = A.HamiltonianExpansionTerm(
            [['x']*2],order=1,identities=[0],logger=False)
        with self.assertRaises(ValueError):evaluator.get_action_expansion()
        projection = _OccupationCoefficientProjection((1,))
        transform = projection.compiler()
        mask = A._ChangeDomainIndex(True,0)
        bound = A.PTEnergyChangeProductSum({((mask,0),):A.PolyPath.from_coeffs([[1]],steps=0)})
        with self.assertRaises(ValueError):transform(((1,0,0,0,0),),bound)
        evaluator = self.model(ndim=2)
        evaluator.solver.hamiltonian_expansion[2] = A.HamiltonianExpansionTerm(
            [['x']*6],order=2,identities=[0],logger=False)
        with self.assertRaises(ValueError):evaluator.get_action_expansion()


class AnyOrderActionExpansionTests(unittest.TestCase):
    @staticmethod
    def model(terms, dimensions=1, compile=True, folder=None):
        # Each entry is (PT order, operator degree, sparse symmetric tensor data).
        w = np.array([.043,.067,.091,.119][:dimensions])
        base = A.AnalyticPerturbationTheorySolver.from_order(2,internals=False,logger=False)
        maximum = max((order for order,_,_ in terms),default=0)
        H = [base.hamiltonian_expansion[0]]+[A.ZeroTerm(logger=False) for _ in range(maximum)]
        coefficients = [[np.diag(w)/2]*2]+[[0]*5 for _ in range(maximum)]
        for order,degree,data in terms:
            H[order] = A.HamiltonianExpansionTerm([['x']*degree],order=order,identities=[0],logger=False)
            tensor = np.zeros((dimensions,)*degree)
            for indices,value in data:
                for permutation in A.UniquePermutations(indices).permutations():
                    tensor[tuple(permutation)] = value
            coefficients[order][0] = tensor
        solver = A.AnalyticPerturbationTheorySolver(H,logger=False)
        return A.PerturbationTheoryEvaluator(solver,coefficients,w,compile=compile,expressions_folder=folder)

    @staticmethod
    def reconstruct(tensors, states):
        actions = np.asarray(states)+.5
        values = np.full(len(actions),tensors[0],dtype=float)
        for degree,tensor in enumerate(tensors[1:],1):
            for index in itertools.combinations_with_replacement(range(actions.shape[1]),degree):
                values += tensor[index]*np.prod(actions[:,index],axis=1)
        return values

    def assertTensors(self, left, right):
        self.assertEqual(len(left),len(right))
        for a,b in zip(left,right):np.testing.assert_allclose(a,b,atol=1e-13,rtol=1e-12)

    def test_primitive_construction_matches_explicit_assignments(self):
        class RepeatedSigns:
            def __init__(self, values):self.values = values
            def permutations(self):return np.array(list(itertools.permutations(self.values)))

        state = np.array([2,3,1])
        frequencies = np.array([.043,.067,.091])
        for degree in range(1,7):
            fast = A.OperatorExpansionTerm([['x']*degree],order=0,identities=[0],logger=False)
            expected = fast.get_changes()
            with mock.patch.object(A,'UniquePermutations',RepeatedSigns):
                explicit = fast.get_changes()
            self.assertEqual(list(expected),list(explicit))
            # An equivalent user symmetrizer retains the general assignment
            # algorithm, giving an independent check on multinomial counting.
            slow = A.OperatorExpansionTerm([['x']*degree],order=0,identities=[0],
                symmetrizers=[lambda indices:tuple(sorted(indices))],logger=False)
            tensor = (1+np.sum(np.indices((3,)*degree)**2,axis=0))*.001
            changes = [c for c in ((),(1,),(2,),(1,-1),(3,),(2,-1))
                       if sum(map(abs,c))<=degree and sum(map(abs,c))%2==degree%2]
            for change in changes:
                values = []
                for term in (fast,slow):
                    wrapper = A.PerturbationTheoryExpressionEvaluator(term,term.get_poly_terms(change),change)
                    values.append(wrapper.evaluate([(state,np.arange(3)[None])],[[tensor]],frequencies,
                        evaluation_mode='dag',evaluation_backend='indexed',zero_cutoff=0.)[0])
                np.testing.assert_allclose(values[0],values[1],atol=1e-12,rtol=1e-12)

    def test_arbitrary_order_patterns_and_validation(self):
        for order in range(13):
            powers = A.AnalyticPerturbationTheorySolver.action_expansion_powers(order)
            self.assertEqual(len(powers),2**(order//2+1))
            self.assertEqual(len(powers),len(set(powers)))
            self.assertIn((),powers)
            for pattern in powers:
                self.assertTrue(all(p>0 for p in pattern))
                self.assertLessEqual(sum(pattern),order//2+1)
        fourth = A.AnalyticPerturbationTheorySolver.action_expansion_powers(4)
        self.assertTrue({(3,),(2,1),(1,2),(1,1,1)} <= set(fourth))
        sixth = A.AnalyticPerturbationTheorySolver.action_expansion_powers(6)
        self.assertTrue({(3,1),(1,3),(2,2),(2,1,1),(1,2,1),(1,1,2),(1,1,1,1)} <= set(sixth))
        evaluator = self.model([])
        for invalid in (-1,2.5,True,None):
            with self.assertRaises(ValueError):evaluator.get_action_expansion(order=invalid)

    def test_harmonic_and_odd_requested_orders(self):
        evaluator = self.model([])
        for source in ('energy','hamiltonian'):
            reference,k = evaluator.get_action_expansion(order=0,correction=source)
            self.assertTensors(reference,[0.,evaluator.freqs])
            self.assertEqual(set(k),{(0,()),(0,(1,))})
            actual,k = evaluator.get_action_expansion(order=1,correction=source)
            self.assertTensors(actual,reference)
            self.assertTrue(all(n%2==0 for n,_ in k))
        evaluator = self.model([(2,2,[((0,0),.001)])])
        even,_ = evaluator.get_action_expansion(order=4)
        odd,k = evaluator.get_action_expansion(order=5)
        self.assertTensors(even,odd)
        self.assertEqual({n for n,_ in k},{0,2,4})

    def test_quadratic_recurrence_through_eighth_order_has_no_spurious_actions(self):
        g = .001
        for compiled in (True,False):
            evaluator = self.model([(2,2,[((0,0),g)])],compile=compiled)
            w = evaluator.freqs[0]
            binomial = 1.
            expected = w
            for k in range(1,5):
                binomial *= (.5-(k-1))/k
                expected += w*binomial*(2*g/w)**k
            for source in ('energy','hamiltonian'):
                tensors,kernels = evaluator.get_action_expansion(order=8,correction=source)
                self.assertEqual(len(tensors),6)
                self.assertEqual({n for n,_ in kernels},{0,2,4,6,8})
                self.assertAlmostEqual(tensors[0],0.,places=13)
                self.assertAlmostEqual(tensors[1][0],expected,places=13)
                for tensor in tensors[2:]:np.testing.assert_allclose(tensor,0.,atol=1e-13,rtol=0)

    def test_cubic_action_tensor_has_distinct_mixed_powers_and_three_anchors(self):
        a,b,c = .000002,.000003,.0000004
        terms = [(4,6,[((0,0,0,0,1,1),a),((0,0,1,1,1,1),b),((0,0,1,1,2,2),c)])]
        for compiled in (True,False):
            evaluator = self.model(terms,dimensions=3,compile=compiled)
            expected = [0.,evaluator.freqs+np.array([45*b/8,45*a/8,0]),
                        np.zeros((3,3)),np.zeros((3,3,3))]
            expected[3][0,0,1] = 45*a/2
            expected[3][0,1,1] = 45*b/2
            expected[3][0,1,2] = 90*c
            for source in ('energy','hamiltonian'):
                tensors,k = evaluator.get_action_expansion(order=4,correction=source)
                self.assertTensors(tensors,expected)
                self.assertTrue({(4,(2,1)),(4,(1,2)),(4,(1,1,1))} <= set(k))

    def test_quartic_action_tensor_has_four_anchors_and_quantum_lower_terms(self):
        a,b = .0000002,.00000003
        terms = [(6,8,[((0,0,0,0,1,1,2,2),a),((0,0,1,1,2,2,3,3),b)])]
        for compiled in (True,False):
            evaluator = self.model(terms,dimensions=4,compile=compiled)
            expected = [0.,evaluator.freqs,np.zeros((4,)*2),np.zeros((4,)*3),np.zeros((4,)*4)]
            expected[2][1,2] = 315*a/2
            expected[4][0,0,1,2] = 630*a
            expected[4][0,1,2,3] = 2520*b
            for source in ('energy','hamiltonian'):
                tensors,k = evaluator.get_action_expansion(order=6,correction=source)
                self.assertTensors(tensors,expected)
                self.assertIn((6,(2,1,1)),k)
                self.assertIn((6,(1,1,1,1)),k)

    def test_twelfth_order_direct_term_generates_seventh_rank(self):
        from Psience.BasisReps import HarmonicOscillatorBasis
        strength = 2e-12
        evaluator = self.model([(12,14,[((0,)*14,strength)])])
        q = HarmonicOscillatorBasis.qmatrix_ho(40).toarray()
        states = np.arange(4)[:,None]
        expected = (states[:,0]+.5)*evaluator.freqs[0]+strength*np.diag(np.linalg.matrix_power(q,14))[:4]
        for source in ('energy','hamiltonian'):
            tensors,k = evaluator.get_action_expansion(order=12,correction=source)
            self.assertEqual(len(tensors),8)
            self.assertAlmostEqual(float(tensors[7].ravel()[0]),strength*math.comb(14,7)/2**7,places=18)
            np.testing.assert_allclose(self.reconstruct(tensors,states),expected,atol=1e-12,rtol=1e-12)
            self.assertEqual({n for n,_ in k},set(range(0,13,2)))

    def test_mixed_cubic_quartic_recurrence_matches_native_matrix_vpt(self):
        from Psience.BasisReps import HarmonicOscillatorProductBasis,BasisStateSpace,ArrayStateIndexer
        from Psience.VPT2.Solver import PerturbationTheorySolver
        from McUtils.Scaffolding import NullLogger,NullCheckpointer
        a,b = .0003,.00002
        w = .043
        basis = HarmonicOscillatorProductBasis([40],indexer=ArrayStateIndexer([40]))
        def representation(word,value):
            rank = len(word)
            return basis.representation(*word,coeffs=np.full((1,)*rank,value),axes=[list(range(rank))]*2)
        h0 = representation(['x']*2,w/2)+representation(['p']*2,w/2)
        states = np.arange(4)[:,None]
        for order,cubic in ((4,a),(6,0.)):
            terms = [(2,4,[((0,0,0,0),b)])]
            if cubic:terms.append((1,3,[((0,0,0),cubic)]))
            evaluator = self.model(terms)
            matrix = PerturbationTheorySolver(
                [h0,representation(['x']*3,cubic),representation(['x']*4,b)]+[0]*(order-2),
                BasisStateSpace(basis,states),order=order,handle_strong_couplings=False,
                logger=NullLogger(),checkpointer=NullCheckpointer(None),results=NullCheckpointer(None),
                zero_element_warning=False)
            expected = matrix.apply_VPT().energy_corrs.sum(axis=1)
            reference = None
            for source in ('energy','hamiltonian'):
                tensors,_ = evaluator.get_action_expansion(order=order,correction=source)
                np.testing.assert_allclose(self.reconstruct(tensors,states),expected,atol=1e-12,rtol=1e-12)
                if reference is None:reference=tensors
                self.assertTensors(tensors,reference)

    def test_higher_order_cache_rebinding_and_warm_disk(self):
        with tempfile.TemporaryDirectory() as folder:
            evaluator = self.model([(2,2,[((0,0),.001)])],folder=folder)
            evaluator.get_action_expansion(order=6)
            count = evaluator.compilation_statistics()['compilations']
            evaluator.freqs *= 1.2
            evaluator.expansions[0] = [np.diag(evaluator.freqs)/2]*2
            evaluator.expansions[2][0] *= 1.7
            changed,_ = evaluator.get_action_expansion(order=6)
            self.assertEqual(evaluator.compilation_statistics()['compilations'],count)
            evaluator.compile_expressions = False
            direct,_ = evaluator.get_action_expansion(order=6)
            self.assertTensors(changed,direct)
            warm = self.model([(2,2,[((0,0),.0017)])],folder=folder)
            warm.freqs,warm.expansions = evaluator.freqs,evaluator.expansions
            result,_ = warm.get_action_expansion(order=6)
            self.assertTensors(result,changed)
            self.assertEqual(warm.compilation_statistics()['compilations'],0)

    def test_all_requested_hamiltonian_orders_are_validated(self):
        evaluator = self.model([(4,8,[((0,)*8,.0001)])])
        # The unused invalid fourth-order term does not affect VPT2.
        evaluator.get_action_expansion(order=2)
        with self.assertRaisesRegex(ValueError,'degree bound'):evaluator.get_action_expansion(order=4)
        evaluator = self.model([(4,5,[((0,)*5,.0001)])])
        with self.assertRaisesRegex(ValueError,'parity'):evaluator.get_action_expansion(order=4)

    def test_higher_order_construction_is_symbolic_and_resonance_independent(self):
        evaluator = self.model([(4,6,[((0,)*6,.000001)])])
        evaluator.solver.resonance_condition = A.ResonanceBlockCondition()
        evaluator.solver.filtered_block_recurrence = True
        with (mock.patch.object(A.PolyPath,'evaluate_polynomial',side_effect=AssertionError('state sampling')),
              mock.patch.object(A.PolyPath,'to_eager',side_effect=AssertionError('polynomial flattening')),
              mock.patch.object(A.PTTensorCoeffProductDAG,'to_eager',side_effect=AssertionError('tensor flattening')),
              mock.patch.object(A.ConditionalTransformationStack,'extend',side_effect=AssertionError('Conditional route'))):
            for source in ('energy','hamiltonian'):
                kernels = evaluator.solver.action_expansion_terms(order=4,correction=source)
                self.assertEqual(len(kernels),24)
                self.assertTrue(all(not kernel.poly_obj.has_change_domains for kernel in kernels.values()))


if __name__ == '__main__':
    unittest.main()
