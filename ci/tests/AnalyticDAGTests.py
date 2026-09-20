import inspect
import os
import tempfile
import unittest

import numpy as np

try:
    import Psience.VPT2.Analytic as Analytic
    from Psience.VPT2 import AnalyticVPTRunner
except ModuleNotFoundError:
    import Psience.Psience.VPT2.Analytic as Analytic
    from Psience.Psience.VPT2 import AnalyticVPTRunner


class AnalyticDAGTests(unittest.TestCase):

    @staticmethod
    def _evaluate(poly, states):
        if Analytic.nput.is_numeric(poly):
            return np.full(len(states), poly, dtype=float)
        return Analytic.PolyPath.from_polynomial(poly).evaluate_polynomial(states)

    @classmethod
    def _flatten(cls, poly, key=()):
        if Analytic.nput.is_numeric(poly):
            return {} if poly == 0 else {key: poly}
        if isinstance(poly, Analytic.SqrtChangePoly):
            return cls._flatten(
                poly.poly_obj,
                key + (('sqrt', tuple(poly.poly_change), tuple(poly.shift_start)),)
            )
        if isinstance(poly, Analytic.PTTensorCoeffProductSum):
            flattened = {}
            for subkey, value in poly.terms.items():
                flattened.update(cls._flatten(value, key + (('tensor', subkey),)))
            return flattened
        if isinstance(poly, Analytic.PTEnergyChangeProductSum):
            flattened = {}
            for subkey, value in poly.terms.items():
                flattened.update(cls._flatten(value, key + (('energy', subkey),)))
            return flattened
        return {key: poly}

    def assertPolynomialEqual(self, eager, path, seed=12981):
        eager = self._flatten(eager)
        path = self._flatten(path)
        rng = np.random.default_rng(seed)
        for key in eager.keys() | path.keys():
            eager_poly = eager.get(key, 0)
            path_poly = path.get(key, 0)
            ndim = max(
                0 if Analytic.nput.is_numeric(eager_poly) else eager_poly.ndim,
                0 if Analytic.nput.is_numeric(path_poly) else path_poly.ndim
            )
            states = rng.integers(0, 8, size=(12, ndim))
            np.testing.assert_allclose(
                self._evaluate(eager_poly, states),
                self._evaluate(path_poly, states),
                rtol=2e-12,
                atol=2e-12,
                err_msg=str(key)
            )

    def test_poly_path_algebra(self):
        left = Analytic.ProductPTPolynomial([
            np.array([1.0, 2.0]),
            np.array([3.0, -1.0, 0.5])
        ], prefactor=2, steps=0)
        right = Analytic.ProductPTPolynomial([
            np.array([0.5, 1.0]),
            np.array([2.0, 1.0])
        ], prefactor=-3, steps=0)
        path_left = Analytic.PolyPath.from_polynomial(left)
        path_right = Analytic.PolyPath.from_polynomial(right)
        cases = [
            (left + right, path_left + path_right),
            (left.shift([1, -2]), path_left.shift([1, -2])),
            (left.permute([1, 0]), path_left.permute([1, 0])),
            (left.mul_simple(right), path_left.mul_simple(path_right)),
            (
                left.mul_along(right, [[0], [1]]),
                path_left.mul_along(path_right, [[0], [1]])
            ),
            (
                left.permute([0, 1]) + left.permute([1, 0]),
                path_left.permutation_sum([[0, 1], [1, 0]])
            )
        ]
        for eager, path in cases:
            ndim = max(eager.ndim, path.ndim)
            states = np.arange(3 * ndim, dtype=int).reshape(3, ndim) % 6
            np.testing.assert_allclose(
                self._evaluate(eager, states),
                self._evaluate(path, states),
                rtol=2e-12,
                atol=2e-12
            )

        full = Analytic.PolyPath.from_coeffs(
            [[2.0], [1.0, 2.0], [3.0]], prefactor=4, steps=0
        )
        condensed, projected = full.condense(return_inds=True)
        np.testing.assert_array_equal(condensed, [0, 2])
        projected_states = np.array([[0], [2], [5]])
        full_states = np.column_stack([
            np.zeros(3, dtype=int),
            projected_states[:, 0],
            np.zeros(3, dtype=int)
        ])
        np.testing.assert_allclose(
            full.evaluate_polynomial(full_states),
            projected.evaluate_polynomial(projected_states)
        )
        np.testing.assert_allclose(
            self._evaluate(projected.to_eager(), projected_states),
            projected.evaluate_polynomial(projected_states)
        )

    def test_materialized_path_cache_reuses_shared_nodes_and_is_bounded(self):
        shared = Analytic.PolyPath.from_coeffs([
            [1.0, 2.0, -0.25],
            [0.5, -1.0]
        ])
        left = shared.scale(2.0)
        right = shared.scale(-3.0)
        states = np.array([
            [0, 1],
            [2, 3],
            [4, 2],
            [1, 5]
        ])
        perm_substates = states[np.newaxis, :, :]
        tuple_states = [[tuple(state) for state in states]]
        cache = Analytic._MaterializedEvaluationCache(
            {}, max_items=4, max_bytes=4096
        )

        values = []
        for poly in (left, right):
            values.append(
                Analytic.PerturbationTheoryExpressionEvaluator._eval_poly(
                    cache,
                    tuple_states,
                    perm_substates,
                    None,
                    poly,
                    [],
                    None,
                    False,
                    Analytic.Logger.lookup(None)
                )
            )

        expected = shared.evaluate_polynomial(states)
        np.testing.assert_allclose(values[0][:, 0], 2 * expected)
        np.testing.assert_allclose(values[1][:, 0], -3 * expected)
        stats = cache.stats()
        self.assertGreater(stats['cache_hits'], 0)
        self.assertLessEqual(stats['cache_peak_items'], 4)
        self.assertLessEqual(stats['cache_peak_bytes'], 4096)

    def test_vectorized_prefactors_match_scalar_reference(self):
        rng = np.random.default_rng(9182)
        permutations = np.array([
            rng.permutation(5) for _ in range(12)
        ])
        coefficient_indices = (
            ((1, 0), (0, 2)),
            ((2, 0), (1,)),
            ((0, 0), ())
        )
        tensors = [
            [rng.normal(size=(5, 5)), rng.normal(size=5), 1.25],
            [rng.normal(size=(5, 5)), rng.normal(size=5), -.75]
        ]
        tensors[0][0][0, 0] = 0
        tensors[1][1][2] = 0
        cutoff = 1e-12

        expected_products = np.ones((len(tensors), len(permutations)))
        expected_good = np.full(expected_products.shape, True)
        for coefficient_pos, (_, indices) in enumerate(coefficient_indices):
            for tensor_pos, tensor_list in enumerate(tensors):
                tensor = tensor_list[coefficient_pos]
                if len(indices) == 0:
                    expected_products[tensor_pos] *= tensor
                else:
                    for perm_pos in range(len(permutations)):
                        if expected_good[tensor_pos, perm_pos]:
                            index = tuple(
                                permutations[perm_pos, axis] for axis in indices
                            )
                            expected_products[tensor_pos, perm_pos] *= tensor[index]
                            if abs(expected_products[tensor_pos, perm_pos]) < cutoff:
                                expected_good[tensor_pos, perm_pos] = False

        products, good = (
            Analytic.PerturbationTheoryExpressionEvaluator
            ._compute_coefficient_prefactor_products(
                permutations, coefficient_indices, tensors, cutoff
            )
        )
        np.testing.assert_allclose(products, expected_products)
        np.testing.assert_array_equal(good, np.any(expected_good, axis=0))

    def test_polynomial_block_cache_bypasses_scalar_lookup(self):
        poly = Analytic.PolyPath.from_coeffs([
            [1.0, 2.0, -0.25],
            [0.5, -1.0]
        ])
        states = np.array([
            [0, 1],
            [2, 3],
            [4, 2],
            [1, 5]
        ])
        perm_substates = states[np.newaxis, :, :]
        tuple_states = [[tuple(state) for state in states]]
        cache = Analytic._MaterializedEvaluationCache(
            {}, max_items=128, max_bytes=64 * 1024
        )
        block_key = (
            'full', Analytic._StatePermutationBlockIdentity(tuple_states)
        )

        first = Analytic.PerturbationTheoryExpressionEvaluator._eval_poly(
            cache, tuple_states, perm_substates, None,
            poly, [], None, False, Analytic.Logger.lookup(None),
            block_key=block_key
        )
        exact_hits = cache.exact_hits
        second = Analytic.PerturbationTheoryExpressionEvaluator._eval_poly(
            cache, tuple_states, perm_substates, None,
            poly, [], None, False, Analytic.Logger.lookup(None),
            block_key=(
                'full', Analytic._StatePermutationBlockIdentity(tuple_states)
            )
        )

        self.assertIs(first, second)
        self.assertEqual(cache.exact_hits, exact_hits)
        self.assertEqual(cache.block_hits, 1)
        self.assertEqual(cache.block_misses, 1)

    def test_state_permutation_gathers_are_cached(self):
        state = np.array([[0, 1, 2, 3, 4]])
        frequencies = np.linspace(.5, 1.5, 5)
        permutations = np.array([
            [0, 2, 4],
            [4, 1, 3],
            [2, 3, 0]
        ])
        cache = {}
        args = (
            0, state, frequencies, 0, (0, 2, 4), permutations, cache,
            (0, 1, 2), (0, 2, 4)
        )
        first = Analytic.PerturbationTheoryExpressionEvaluator._get_state_perms(*args)
        second = Analytic.PerturbationTheoryExpressionEvaluator._get_state_perms(*args)

        self.assertIs(first, second)
        np.testing.assert_array_equal(
            first[0], np.moveaxis(Analytic.nput.vector_take(state, permutations), 0, 1)
        )
        np.testing.assert_array_equal(
            first[2], Analytic.nput.vector_take(frequencies, permutations)
        )
        self.assertEqual(len(cache), 1)
        self.assertIsInstance(first[3], Analytic._StatePermutationBlockIdentity)
        self.assertEqual(
            first[3].states,
            tuple(tuple(block) for block in first[1])
        )

    def test_state_permutation_cache_is_byte_bounded(self):
        state = np.array([[0, 1, 2, 3, 4]])
        frequencies = np.linspace(.5, 1.5, 5)
        permutations = np.array([
            [0, 2, 4],
            [4, 1, 3],
            [2, 3, 0]
        ])
        cache = Analytic._StatePermutationCache(
            max_items=100,
            max_bytes=4096
        )
        for state_idx in range(12):
            Analytic.PerturbationTheoryExpressionEvaluator._get_state_perms(
                state_idx, state, frequencies, 0, (0, 2, 4),
                permutations, cache, (0, 1, 2), (0, 2, 4)
            )

        stats = cache.stats()
        self.assertLessEqual(stats['cache_peak_items'], 100)
        self.assertLessEqual(stats['cache_peak_bytes'], 4096)
        self.assertGreater(stats['cache_evictions'], 0)

    def test_correction_backend_parity(self):
        builders = [
            ('energy', 4, lambda solver: solver.energy_correction(2)([])),
            ('wavefunction+', 3, lambda solver: solver.wavefunction_correction(1)([1])),
            ('wavefunction-', 3, lambda solver: solver.wavefunction_correction(1)([-1])),
            ('wavefunction-2/+2', 4, lambda solver: solver.wavefunction_correction(2)([2])),
            ('wavefunction-2/+1+1', 4, lambda solver: solver.wavefunction_correction(2)([1, 1])),
            ('overlap', 4, lambda solver: solver.overlap_correction(2)([]))
        ]
        for name, order, builder in builders:
            expressions = {}
            for representation in ('eager', 'path'):
                Analytic.AnalyticPerturbationTheorySolver.clear_caches()
                solver = Analytic.AnalyticPerturbationTheorySolver.from_order(
                    order,
                    polynomial_representation=representation
                )
                expressions[representation] = builder(solver).expr
                if representation == 'path' and isinstance(
                        expressions[representation], Analytic.SqrtChangePoly
                ):
                    self.assertIsInstance(
                        expressions[representation].poly_obj,
                        Analytic.PTTensorCoeffProductDAG
                    )
                    self.assertEqual(
                        Analytic.PTTensorCoeffProductDAG.cache_info()['tensor_materializations'],
                        0
                    )
            with self.subTest(name=name):
                self.assertPolynomialEqual(expressions['eager'], expressions['path'])

    def test_full_energy_evaluator_backend_parity(self):
        evaluators = {}
        coefficient_keys = set()
        for representation in ('eager', 'path'):
            Analytic.AnalyticPerturbationTheorySolver.clear_caches()
            solver = Analytic.AnalyticPerturbationTheorySolver.from_order(
                4,
                polynomial_representation=representation
            )
            evaluators[representation] = solver.energy_correction(2)([])
            coefficient_keys.update(
                coefficient
                for product in evaluators[representation].expr.poly_obj.terms
                for coefficient in product
            )

        rng = np.random.default_rng(42)
        coefficient_expansion = [
            [] for _ in range(max(key[0] for key in coefficient_keys) + 1)
        ]
        for coefficient_type, expansion in enumerate(coefficient_expansion):
            orders = [key[1] for key in coefficient_keys if key[0] == coefficient_type]
            for order in range(max(orders, default=0) + 1):
                ranks = [
                    len(key) - 2
                    for key in coefficient_keys
                    if key[:2] == (coefficient_type, order)
                ]
                expansion.append(
                    0 if len(ranks) == 0 else rng.normal(size=(2,) * max(ranks))
                )

        state_permutations = [
            np.array([1, 2]),
            np.array([[0, 1], [1, 0]])
        ]
        frequencies = np.array([1.1, 1.7])
        values = {}
        for representation, evaluator in evaluators.items():
            Analytic.PerturbationTheoryExpressionEvaluator._cached_expansion = None
            values[representation] = evaluator.evaluate(
                state_permutations,
                coefficient_expansion,
                frequencies
            )
        np.testing.assert_allclose(
            values['eager'], values['path'], rtol=2e-12, atol=2e-12
        )

        tensor_dag = evaluators['path'].expr.poly_obj
        materialized = tensor_dag.to_eager()
        materialization_count = Analytic.PTTensorCoeffProductDAG.cache_info()[
            'tensor_materializations'
        ]
        self.assertIs(tensor_dag.to_eager(), materialized)
        self.assertEqual(
            Analytic.PTTensorCoeffProductDAG.cache_info()['tensor_materializations'],
            materialization_count
        )

        excluded_operator = next(iter(tensor_dag.operator_keys))
        pruned_eager = materialized.prune_operators([excluded_operator])
        pruned_dag = tensor_dag.prune_operators([excluded_operator]).to_eager()
        self.assertPolynomialEqual(pruned_eager, pruned_dag)

    def test_direct_evaluator_is_bounded_and_matches_materialization(self):
        Analytic.AnalyticPerturbationTheorySolver.clear_caches()
        solver = Analytic.AnalyticPerturbationTheorySolver.from_order(
            4, polynomial_representation='path'
        )
        evaluator = solver.energy_correction(2)([])
        coefficient_keys = {
            coefficient
            for product in evaluator.expr.poly_obj.to_eager().terms
            for coefficient in product
        }
        rng = np.random.default_rng(812)
        coefficient_expansion = [
            [] for _ in range(max(key[0] for key in coefficient_keys) + 1)
        ]
        for coefficient_type, expansion in enumerate(coefficient_expansion):
            orders = [key[1] for key in coefficient_keys if key[0] == coefficient_type]
            for order in range(max(orders, default=0) + 1):
                ranks = [
                    len(key) - 2
                    for key in coefficient_keys
                    if key[:2] == (coefficient_type, order)
                ]
                expansion.append(
                    0 if len(ranks) == 0 else rng.normal(size=(3,) * max(ranks))
                )

        state_permutations = [
            np.array([1, 2, 0]),
            np.array([[0, 1, 2], [1, 0, 2], [2, 1, 0]])
        ]
        frequencies = np.array([0.8, 1.3, 1.9])
        materialized_legacy = evaluator.evaluate(
            state_permutations, coefficient_expansion, frequencies,
            evaluation_mode='materialized_legacy'
        )

        Analytic.PerturbationTheoryExpressionEvaluator._poly_cache = (
            Analytic.PerturbationTheoryExpressionEvaluator.get_cache()
        )
        materialized = evaluator.evaluate(
            state_permutations, coefficient_expansion, frequencies,
            evaluation_mode='materialized',
            path_cache_size=8,
            path_cache_bytes=4096
        )
        np.testing.assert_allclose(
            materialized, materialized_legacy, rtol=2e-12, atol=2e-12
        )
        materialized_stats = (
            Analytic.PerturbationTheoryExpressionEvaluator
            .get_last_materialized_evaluation_stats()
        )
        self.assertLessEqual(materialized_stats['cache_peak_items'], 8)
        self.assertLessEqual(materialized_stats['cache_peak_bytes'], 4096)

        Analytic.PerturbationTheoryExpressionEvaluator._poly_cache = (
            Analytic.PerturbationTheoryExpressionEvaluator.get_cache()
        )
        direct_legacy = evaluator.evaluate(
            state_permutations, coefficient_expansion, frequencies,
            evaluation_mode='dag_legacy',
            dag_cache_size=8,
            dag_cache_bytes=4096,
            dag_chunk_size=3
        )

        # Force both item- and byte-pressure so this exercises eviction rather
        # than merely checking the configured limits.  Reset the exact cache so
        # DAG reuse is measured within this evaluation batch rather than being
        # inherited from the materialized or legacy reference runs above.
        Analytic.PerturbationTheoryExpressionEvaluator._poly_cache = (
            Analytic.PerturbationTheoryExpressionEvaluator.get_cache()
        )
        direct = evaluator.evaluate(
            state_permutations, coefficient_expansion, frequencies,
            evaluation_mode='dag',
            dag_cache_size=8,
            dag_cache_bytes=4096,
            dag_chunk_size=3
        )
        np.testing.assert_allclose(direct_legacy, materialized_legacy, rtol=2e-12, atol=2e-12)
        np.testing.assert_allclose(direct, materialized_legacy, rtol=2e-12, atol=2e-12)
        stats = Analytic.PerturbationTheoryExpressionEvaluator.get_last_dag_evaluation_stats()
        self.assertLessEqual(stats['cache_peak_items'], 8)
        self.assertLessEqual(stats['cache_peak_bytes'], 4096)
        self.assertGreater(stats['cache_evictions'], 0)
        self.assertTrue(stats['dag_path_cache_enabled'])
        self.assertGreater(stats['dag_exact_cache_hits'], 0)
        self.assertEqual(stats['dag_materializations'], 0)

    def test_fourth_order_derivation_stays_lazy(self):
        Analytic.AnalyticPerturbationTheorySolver.clear_caches()
        solver = Analytic.AnalyticPerturbationTheorySolver.from_order(
            6,
            polynomial_representation='path'
        )
        expression = solver.energy_correction(4)([]).expr
        self.assertIsInstance(
            expression.poly_obj,
            Analytic.PTTensorCoeffProductDAG
        )
        cache_info = solver.polynomial_cache_info()
        self.assertEqual(cache_info['tensor_materializations'], 0)
        self.assertEqual(cache_info['axis_materializations'], 0)
        self.assertLess(cache_info['tensor_nodes'], 50000)

    def test_path_is_default_and_eager_remains_selectable(self):
        solver = Analytic.AnalyticPerturbationTheorySolver.from_order(3)
        self.assertEqual(solver.polynomial_representation, 'path')
        eager = Analytic.AnalyticPerturbationTheorySolver.from_order(
            3,
            polynomial_representation='eager'
        )
        self.assertEqual(eager.polynomial_representation, 'eager')
        self.assertEqual(
            inspect.signature(AnalyticVPTRunner.construct)
                .parameters['polynomial_representation'].default,
            'path'
        )

    def test_path_checkpoint_roundtrip(self):
        Analytic.AnalyticPerturbationTheorySolver.clear_caches()
        with tempfile.TemporaryDirectory() as tmpdir:
            checkpoint = os.path.join(tmpdir, 'analytic-path.hdf5')
            solver = Analytic.AnalyticPerturbationTheorySolver.from_order(
                4,
                checkpoint=checkpoint,
                polynomial_representation='path'
            )
            original = solver.energy_correction(2)([])
            self.assertEqual(
                Analytic.PTTensorCoeffProductDAG.cache_info()['tensor_materializations'],
                0
            )

            restored_solver = Analytic.AnalyticPerturbationTheorySolver.from_order(
                4,
                checkpoint=checkpoint,
                polynomial_representation='path'
            )
            restored = restored_solver.energy_correction(2)([])
            self.assertTrue(Analytic.polynomial_uses_path(restored.expr))
            self.assertIsInstance(
                restored.expr.poly_obj,
                Analytic.PTTensorCoeffProductDAG
            )
            self.assertPolynomialEqual(original.expr, restored.expr)


if __name__ == '__main__':
    unittest.main()
