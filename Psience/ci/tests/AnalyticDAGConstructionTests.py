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


if __name__ == '__main__':
    unittest.main()
