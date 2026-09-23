import inspect
import os
import tempfile
import unittest
import warnings
import weakref
from unittest import mock

import numpy as np

try:
    import Psience.VPT2.Analytic as Analytic
    from Psience.VPT2 import AnalyticVPTRunner
    from McUtils.Parallelizers import MultiprocessingParallelizer
except ModuleNotFoundError:
    import Psience.Psience.VPT2.Analytic as Analytic
    from Psience.Psience.VPT2 import AnalyticVPTRunner
    from McUtils.McUtils.Parallelizers import MultiprocessingParallelizer


class AnalyticDAGTests(unittest.TestCase):

    def test_perturbation_operator_projects_out_reference_state(self):
        Analytic.AnalyticPerturbationTheorySolver.clear_caches()
        solver = Analytic.AnalyticPerturbationTheorySolver.from_order(
            4, polynomial_representation='path'
        )
        first_order = solver.wavefunction_correction(1)
        resolvent = Analytic.PerturbationOperator.lookup(first_order)
        self.assertNotIn((), resolvent.changes)
        self.assertEqual(resolvent.get_poly_terms(()), 0)
        self.assertEqual(resolvent.get_poly_terms((1,), shift=(-1,)), 0)
        self.assertTrue(
            Analytic.nput.is_zero(solver.wavefunction_correction(2)([]).expr)
        )
        self.assertFalse(
            Analytic.nput.is_zero(solver.overlap_correction(2)([]).expr)
        )

    def test_small_energy_denominator_is_pruned_after_polynomial_evaluation(self):
        evaluator = Analytic.PerturbationTheoryExpressionEvaluator
        numerator = np.array([[1e-20, 1e-6, 1e-20]])
        denominator = np.array([0., 1e-20, 2.])
        states = [[(0,)], [(1,)], [(2,)]]
        permutation_states = np.array([[[0]], [[1]], [[2]]])
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            result = evaluator._divide_polynomial_by_energy(
                numerator, denominator, 1e-18,
                energy_changes=((0, 0),),
                evaluation_context={'operator': 'W[2]', 'change': ()},
                tuple_states=states,
                permutation_substates=permutation_states
            )

        self.assertEqual(result[0, 0], 0)
        self.assertEqual(result[0, 2], 0)
        self.assertEqual(result[0, 1], 1e14)
        self.assertFalse(np.any(np.isnan(result)))
        self.assertEqual(len(caught), 2)
        self.assertIs(caught[0].category, Analytic.SmallEnergyDenominatorWarning)
        self.assertIs(
            caught[1].category, Analytic.SurvivingEnergyDenominatorWarning
        )
        self.assertIn("'operator': 'W[2]'", str(caught[0].message))
        self.assertIn("'states':", str(caught[1].message))
        self.assertIn("'permutation_states':", str(caught[1].message))

        with warnings.catch_warnings(record=True) as ordinary_warnings:
            warnings.simplefilter('always')
            ordinary = evaluator._divide_polynomial_by_energy(
                np.array([[2., 1e-20]]), np.array([4., 2.]), 1e-18,
                evaluation_context='ordinary'
            )
        np.testing.assert_array_equal(ordinary, np.array([[.5, 0.]]))
        self.assertEqual(ordinary_warnings, [])

    def test_block_energy_division_matches_row_evaluation(self):
        evaluator = Analytic.PerturbationTheoryExpressionEvaluator
        numerators = np.array([
            [[2., 1e-20, 6.], [4., 8., 1e-20]],
            [[1e-20, 3., 5.], [1e-20, 9., 10.]]
        ])
        denominators = np.array([
            [4., 2., 3.],
            [0., 3., 5.]
        ])
        diagnostic_calls = []

        def diagnostics(row):
            diagnostic_calls.append(row)
            return {
                'evaluation_context': 'row-{}'.format(row),
                'energy_changes': ((row, -row),)
            }

        with warnings.catch_warnings(record=True) as block_warnings:
            warnings.simplefilter('always')
            block = evaluator._divide_polynomial_block_by_energy(
                numerators, denominators, 1e-18,
                row_diagnostics=diagnostics
            )
        with warnings.catch_warnings(record=True) as row_warnings:
            warnings.simplefilter('always')
            rows = np.stack([
                evaluator._divide_polynomial_by_energy(
                    numerator, denominator, 1e-18,
                    evaluation_context='row-{}'.format(row),
                    energy_changes=((row, -row),)
                )
                for row, (numerator, denominator) in enumerate(zip(
                    numerators, denominators
                ))
            ])

        np.testing.assert_array_equal(block, rows)
        self.assertEqual(diagnostic_calls, [1])
        self.assertEqual(
            [warning.category for warning in block_warnings],
            [warning.category for warning in row_warnings]
        )
        self.assertIn('row-1', str(block_warnings[0].message))

    def test_parallel_evaluation_blocks_are_balanced_and_complete(self):
        partition = (
            Analytic.PerturbationTheoryExpressionEvaluator
            ._partition_evaluation_blocks
        )
        for combinations, processes in ((0, 4), (3, 5), (66, 4), (741, 4)):
            blocks = partition(combinations, processes)
            cursor = 0
            sizes = []
            for start, stop in blocks:
                self.assertEqual(start, cursor)
                self.assertGreater(stop, start)
                sizes.append(stop - start)
                cursor = stop
            self.assertEqual(cursor, combinations)
            if sizes:
                self.assertLessEqual(max(sizes) - min(sizes), 1)

        self.assertEqual(
            partition(66, 4),
            [(0, 17), (17, 34), (34, 50), (50, 66)]
        )

        aggregate = (
            Analytic.PerturbationTheoryExpressionEvaluator
            ._aggregate_indexed_evaluation_stats
        )
        worker_stats = [
            {'calls': 3, 'routing_seconds': 2.0, 'max_rows': 7},
            {'calls': 5, 'routing_seconds': 3.0, 'max_rows': 4}
        ]
        self.assertEqual(
            aggregate(worker_stats, elapsed_mode='max'),
            {'calls': 8, 'routing_seconds': 3.0, 'max_rows': 7}
        )
        self.assertEqual(
            aggregate(worker_stats, elapsed_mode='sum'),
            {'calls': 8, 'routing_seconds': 5.0, 'max_rows': 7}
        )

    def test_parallel_dag_evaluation_matches_serial_with_remainder(self):
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
        rng = np.random.default_rng(9417)
        nmodes = 4
        coefficient_expansion = [
            [] for _ in range(max(key[0] for key in coefficient_keys) + 1)
        ]
        for coefficient_type, expansion in enumerate(coefficient_expansion):
            orders = [
                key[1] for key in coefficient_keys
                if key[0] == coefficient_type
            ]
            for order in range(max(orders, default=0) + 1):
                ranks = [
                    len(key) - 2
                    for key in coefficient_keys
                    if key[:2] == (coefficient_type, order)
                ]
                expansion.append(
                    0 if len(ranks) == 0 else
                    rng.normal(scale=.01, size=(nmodes,) * max(ranks))
                )

        state_permutations = [
            np.array([1, 0, 2, 0]),
            np.array([
                [0, 1, 2, 3],
                [1, 0, 3, 2],
                [2, 3, 0, 1]
            ])
        ]
        frequencies = np.array([.8, 1.1, 1.6, 2.0])

        evaluator_cls = Analytic.PerturbationTheoryExpressionEvaluator
        for backend, evaluation_mode in (
                ('legacy', 'dag'),
                ('indexed', 'materialized'),
                ('indexed', 'dag')
        ):
            evaluator_cls._cached_expansion = None
            evaluator_cls._poly_cache = evaluator_cls.get_cache()
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                serial = evaluator.evaluate(
                    state_permutations, coefficient_expansion, frequencies,
                    evaluation_mode=evaluation_mode,
                    evaluation_backend=backend
                )

            evaluator_cls._cached_expansion = None
            evaluator_cls._poly_cache = evaluator_cls.get_cache()
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                with MultiprocessingParallelizer(
                        processes=4, stall_timeout=30
                ) as parallelizer:
                    parallel = evaluator.evaluate(
                        state_permutations, coefficient_expansion, frequencies,
                        evaluation_mode=evaluation_mode,
                        evaluation_backend=backend,
                        parallelizer=parallelizer
                    )
            np.testing.assert_allclose(
                parallel, serial, rtol=2e-12, atol=2e-12
            )
            if backend == 'indexed':
                stats = evaluator_cls.get_last_indexed_evaluation_stats()
                self.assertGreater(stats['calls'], 0)
                self.assertGreater(stats['work_items'], 0)
                self.assertGreater(stats['contraction_rows'], 0)

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
        evaluator = Analytic.PerturbationTheoryExpressionEvaluator
        state = np.array([[0, 1, 2, 3, 4]])
        frequencies = np.linspace(.5, 1.5, 5)
        permutations = np.array([
            [0, 2, 4],
            [4, 1, 3],
            [2, 3, 0]
        ])
        old_state_mode = evaluator.use_hashable_numpy_state_indices
        old_mask_mode = evaluator.use_hashable_numpy_mask_indices
        old_threshold = evaluator.hashable_numpy_index_threshold
        try:
            for use_numpy in (False, True):
                evaluator.use_hashable_numpy_state_indices = use_numpy
                evaluator.use_hashable_numpy_mask_indices = use_numpy
                evaluator.hashable_numpy_index_threshold = 0
                cache = {}
                args = (
                    0, state, frequencies, 0, (0, 2, 4), permutations,
                    cache, np.array([0, 1, 2]), (0, 2, 4)
                )
                first = evaluator._get_state_perms(*args)
                second = evaluator._get_state_perms(*args)

                self.assertIs(first, second)
                np.testing.assert_array_equal(
                    first[0],
                    np.moveaxis(
                        Analytic.nput.vector_take(state, permutations), 0, 1
                    )
                )
                np.testing.assert_array_equal(
                    first[2],
                    Analytic.nput.vector_take(frequencies, permutations)
                )
                self.assertEqual(len(cache), 1)
                self.assertIsInstance(
                    first[3], Analytic._StatePermutationBlockIdentity
                )
                if use_numpy:
                    self.assertTrue(all(
                        isinstance(key, Analytic._HashableNumPyIndex)
                        for block in first[1] for key in block
                    ))
                    self.assertEqual(
                        [[key.tolist() for key in block] for block in first[1]],
                        first[0].tolist()
                    )
                else:
                    self.assertEqual(
                        first[3].states,
                        tuple(tuple(block) for block in first[1])
                    )
        finally:
            evaluator.use_hashable_numpy_state_indices = old_state_mode
            evaluator.use_hashable_numpy_mask_indices = old_mask_mode
            evaluator.hashable_numpy_index_threshold = old_threshold

    def test_hashable_numpy_index_has_immutable_content_semantics(self):
        base = np.arange(12, dtype=np.int16).reshape(3, 4)
        first = Analytic._HashableNumPyIndex(base[:, ::2])
        same = Analytic._HashableNumPyIndex(
            np.ascontiguousarray(base[:, ::2])
        )
        different_dtype = Analytic._HashableNumPyIndex(
            base[:, ::2].astype(np.int32)
        )

        self.assertEqual(first, same)
        self.assertEqual(hash(first), hash(same))
        self.assertNotEqual(first, different_dtype)
        self.assertEqual({first: 'cached'}[same], 'cached')
        expected = first.tolist()
        base[:, ::2] = -1
        self.assertEqual(first.tolist(), expected)
        self.assertEqual(first.nbytes, 6 * np.dtype(np.int16).itemsize)

    def test_hashable_numpy_index_threshold_keeps_short_state_tuples(self):
        evaluator = Analytic.PerturbationTheoryExpressionEvaluator
        old_threshold = evaluator.hashable_numpy_index_threshold
        try:
            evaluator.hashable_numpy_index_threshold = 8
            short = evaluator._index_state_permutations(
                np.arange(8, dtype=np.int64).reshape(1, 1, 8)
            )[0][0]
            long = evaluator._index_state_permutations(
                np.arange(9, dtype=np.int64).reshape(1, 1, 9)
            )[0][0]
        finally:
            evaluator.hashable_numpy_index_threshold = old_threshold

        self.assertIsInstance(short, tuple)
        self.assertIsInstance(long, Analytic._HashableNumPyIndex)
        self.assertEqual(long.tolist(), list(range(9)))

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

    def test_state_permutation_cache_sizing_uses_fixed_layout(self):
        class NonIterableList(list):
            def __iter__(self):
                raise AssertionError('state cache sizing must not recurse')

        states = np.zeros((2, 3, 4), dtype=int)
        frequencies = np.zeros((2, 4), dtype=float)
        tuple_states = NonIterableList([
            [(0, 0, 0, 0)] * 3,
            [(0, 0, 0, 0)] * 3
        ])
        identity = Analytic._StatePermutationBlockIdentity([
            [(0, 0, 0, 0)] * 3,
            [(0, 0, 0, 0)] * 3
        ])
        size = Analytic._StatePermutationCache._value_size((
            states, tuple_states, frequencies, identity
        ))
        self.assertGreater(size, states.nbytes + frequencies.nbytes)

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

    def test_indexed_backend_matches_legacy_and_compare_mode(self):
        from Psience.VPT2.IndexedEvaluator import IndexedBlockEvaluator

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
        rng = np.random.default_rng(9182)
        nmodes = 3
        coefficient_expansion = []
        for coefficient_type in range(max(key[0] for key in coefficient_keys) + 1):
            expansion = []
            orders = [
                key[1] for key in coefficient_keys
                if key[0] == coefficient_type
            ]
            for order in range(max(orders, default=0) + 1):
                ranks = [
                    len(key) - 2
                    for key in coefficient_keys
                    if key[:2] == (coefficient_type, order)
                ]
                expansion.append(
                    0 if len(ranks) == 0 else
                    rng.normal(scale=.02, size=(nmodes,) * max(ranks))
                )
            coefficient_expansion.append(expansion)

        state_permutations = [
            np.array([1, 2, 0]),
            np.array([
                [0, 1, 2],
                [1, 0, 2],
                [2, 1, 0]
            ])
        ]
        frequencies = np.array([.8, 1.3, 1.9])
        old_routing = IndexedBlockEvaluator.use_compiled_permutation_routing
        try:
            for evaluation_mode in ('materialized', 'dag'):
                legacy = evaluator.evaluate(
                    state_permutations,
                    coefficient_expansion,
                    frequencies,
                    evaluation_mode=evaluation_mode,
                    evaluation_backend='legacy'
                )
                IndexedBlockEvaluator.use_compiled_permutation_routing = False
                mask_routed = evaluator.evaluate(
                    state_permutations,
                    coefficient_expansion,
                    frequencies,
                    evaluation_mode=evaluation_mode,
                    evaluation_backend='indexed'
                )
                IndexedBlockEvaluator.use_compiled_permutation_routing = True
                indexed = evaluator.evaluate(
                    state_permutations,
                    coefficient_expansion,
                    frequencies,
                    evaluation_mode=evaluation_mode,
                    evaluation_backend='indexed'
                )
                np.testing.assert_allclose(
                    indexed, legacy, rtol=2e-12, atol=2e-12
                )
                np.testing.assert_array_equal(indexed, mask_routed)
        finally:
            IndexedBlockEvaluator.use_compiled_permutation_routing = old_routing

        compared = evaluator.evaluate(
            state_permutations,
            coefficient_expansion,
            frequencies,
            evaluation_mode='dag',
            evaluation_backend='compare'
        )
        np.testing.assert_allclose(
            compared, legacy, rtol=2e-12, atol=2e-12
        )
        self.assertEqual(
            Analytic.PerturbationTheoryExpressionEvaluator.default_evaluation_backend,
            'indexed'
        )
        defaulted = evaluator.evaluate(
            state_permutations,
            coefficient_expansion,
            frequencies
        )
        np.testing.assert_array_equal(defaulted, indexed)
        stats = (
            Analytic.PerturbationTheoryExpressionEvaluator
            .get_last_indexed_evaluation_stats()
        )
        self.assertGreater(stats['plan_terms'], 0)
        self.assertGreater(stats['work_items'], 0)
        self.assertGreater(stats['contraction_rows'], 0)
        self.assertGreater(stats['contraction_segments'], 0)
        self.assertGreater(stats['degeneracy_selection_hits'], 0)
        self.assertGreater(stats['routing_plans_compiled'], 0)
        self.assertGreater(stats['routing_matches'], 0)
        self.assertLessEqual(stats['max_rows_per_block'], 64)

    def test_indexed_contraction_segments_stream_in_bounded_order(self):
        from Psience.VPT2.IndexedEvaluator import (
            ConcreteBlockTable,
            IndexedContractionSegmentTable
        )

        blocks = ConcreteBlockTable()
        block_index = blocks.append(
            3, np.array([1, 4]), np.array([1, 4]),
            np.ones((1, 2)), ('key',), object(),
            [(0,), (1,)], np.zeros((2, 1)), [(0,), (1,)],
            np.zeros((2, 1)), ('block',), np.zeros((1, 2))
        )
        self.assertEqual(block_index, 0)
        self.assertEqual(len(blocks), 1)
        self.assertEqual(blocks.state_indices[0], 3)

        segments = IndexedContractionSegmentTable()
        first = [object(), object(), object()]
        second = [object(), object()]
        selection = (0, 2)
        segments.append(4, None, first)
        segments.append(7, selection, second)
        streamed = list(segments.iter_blocks(2))
        self.assertTrue(all(len(terms) <= 2 for _, _, terms in streamed))
        self.assertEqual(
            [int(index) for indices, _, _ in streamed for index in indices],
            [4, 4, 4, 7, 7]
        )
        self.assertEqual(
            [term for _, _, terms in streamed for term in terms],
            first + second
        )
        self.assertEqual(
            [value for _, selections, _ in streamed for value in selections],
            [None, None, None, selection, selection]
        )

        identified = list(segments.iter_blocks(2, return_segment_ids=True))
        self.assertEqual(
            [int(segment) for _, _, _, ids in identified for segment in ids],
            [0, 0, 0, 1, 1]
        )
        self.assertEqual(
            [term for _, _, terms, _ in identified for term in terms],
            first + second
        )

    def test_indexed_polynomial_block_uses_one_bounded_materialization(self):
        from Psience.VPT2.IndexedEvaluator import (
            IndexedBlockEvaluator,
            IndexedEvaluationStats
        )

        stats = IndexedEvaluationStats()
        evaluator = object.__new__(IndexedBlockEvaluator)
        evaluator.stats = stats
        full = np.arange(6., dtype=float).reshape(2, 3)
        factors = [2., np.array([1., 2., 3.]), full]
        block = evaluator._materialize_polynomial_block(factors, 2, 3)

        expected = np.empty((3, 2, 3), dtype=float)
        expected[0] = 2.
        expected[1] = np.array([1., 2., 3.])
        expected[2] = full
        np.testing.assert_array_equal(block, expected)
        self.assertTrue(block.flags.owndata)
        self.assertEqual(stats.max_polynomial_scratch_bytes, block.nbytes)
        self.assertEqual(
            block.nbytes,
            len(factors) * 2 * 3 * np.dtype(float).itemsize
        )

        smaller = evaluator._materialize_polynomial_block([1.], 1, 2)
        self.assertEqual(smaller.nbytes, 2 * np.dtype(float).itemsize)
        self.assertEqual(stats.max_polynomial_scratch_bytes, block.nbytes)

    def test_indexed_permutation_routing_matches_mask_split_routing(self):
        from Psience.VPT2.IndexedEvaluator import (
            IndexedPermutationRoutingPlan
        )

        permutation_map = np.array([
            2, 0, 4,
            1, 2,
            4, 3, 2, 0
        ])
        state_lengths = [3, 2, 4]
        eval_permutations = np.array([0, 2, 4, 1, 4, 3])
        eval_coefficients = np.array([0, 0, 0, 2, 2, 5])
        plan = IndexedPermutationRoutingPlan.compile(
            permutation_map,
            state_lengths,
            permutation_count=5
        )
        routed = list(plan.route(eval_permutations, eval_coefficients))

        state_offsets = np.concatenate(([0], np.cumsum(state_lengths)))
        expected = []
        for coefficient in np.unique(eval_coefficients):
            value_group = np.flatnonzero(eval_coefficients == coefficient)
            row_for_permutation = {
                int(eval_permutations[row]): int(row)
                for row in value_group
            }
            for state_index in range(len(state_lengths)):
                state_map = permutation_map[
                    state_offsets[state_index]:state_offsets[state_index + 1]
                ]
                mask_positions = np.array([
                    local_index
                    for local_index, permutation in enumerate(state_map)
                    if int(permutation) in row_for_permutation
                ], dtype=np.intp)
                if len(mask_positions) == 0:
                    continue
                value_positions = np.array([
                    row_for_permutation[int(state_map[local_index])]
                    for local_index in mask_positions
                ], dtype=np.intp)
                expected.append((
                    int(coefficient), state_index,
                    mask_positions, value_positions
                ))

        self.assertEqual(len(routed), len(expected))
        for actual, reference in zip(routed, expected):
            self.assertEqual(actual[:2], reference[:2])
            np.testing.assert_array_equal(actual[2], reference[2])
            np.testing.assert_array_equal(actual[3], reference[3])

    def test_array_degeneracy_plan_matches_legacy_predicate_tree(self):
        tests = [
            Analytic.PerturbationTheoryExpressionEvaluator._deg_test(
                (1, -1, 3, -1)
            ),
            Analytic.PerturbationTheoryExpressionEvaluator._deg_test(
                (-1, 2, 4, -1)
            ),
            Analytic.PerturbationTheoryExpressionEvaluator._deg_test(
                (0, 1, -1, 3)
            )
        ]
        legacy = (
            Analytic.PerturbationTheoryExpressionEvaluator
            ._make_full_deg_test(tests, array_oriented=False)
        )
        compiled = (
            Analytic.PerturbationTheoryExpressionEvaluator
            ._make_full_deg_test(tests, array_oriented=True)
        )
        states = np.array([
            [1, 8, 3, 2],
            [7, 2, 4, 1],
            [0, 1, 9, 3],
            [1, 2, 5, 3],
            [0, 0, 0, 0]
        ])
        expected = np.array([bool(legacy(state)) for state in states])
        np.testing.assert_array_equal(compiled.evaluate(states), expected)
        np.testing.assert_array_equal(
            np.array([compiled(state) for state in states]),
            expected
        )

        exact = (
            Analytic.PerturbationTheoryExpressionEvaluator
            ._make_full_deg_test(
                {(1, 2, 3), (3, 2, 1)},
                array_oriented=True
            )
        )
        np.testing.assert_array_equal(
            exact.evaluate(np.array([
                [1, 2, 3],
                [3, 2, 1],
                [1, 2, 1]
            ])),
            np.array([True, True, False])
        )

    def test_shared_degeneracy_pool_is_exact_and_byte_bounded(self):
        from Psience.VPT2.IndexedEvaluator import (
            IndexedDegeneracyEvaluationContext,
            IndexedEvaluationStats
        )

        left = Analytic.DegeneracyTestPlan(np.array([
            [1, -1, 3],
            [2, -1, 4]
        ]))
        right = Analytic.DegeneracyTestPlan(np.array([
            [-1, 5, 3],
            [-1, 6, 4]
        ]))
        pool = np.array([
            [1, 5, 3],
            [1, 6, 3],
            [2, 6, 4],
            [2, 5, 4],
            [7, 5, 3]
        ])
        stats = IndexedEvaluationStats()
        context = IndexedDegeneracyEvaluationContext(pool, stats)
        for use_left in (False, True):
            for use_right in (False, True):
                for join in (np.logical_and, np.logical_or):
                    actual = context.selection_mask(
                        left, right, use_left, use_right, join
                    )
                    left_mask = np.array([left(state) for state in pool])
                    right_mask = np.array([right(state) for state in pool])
                    if not use_left:
                        left_mask = np.logical_not(left_mask)
                    if not use_right:
                        right_mask = np.logical_not(right_mask)
                    np.testing.assert_array_equal(
                        actual,
                        join(left_mask, right_mask)
                    )
        self.assertEqual(stats.degeneracy_pool_predicate_misses, 2)
        self.assertGreater(stats.degeneracy_pool_predicate_hits, 0)

        bounded_stats = IndexedEvaluationStats()
        bounded = IndexedDegeneracyEvaluationContext(
            np.tile(np.arange(2), (128, 1)),
            bounded_stats,
            max_cache_bytes=1
        )
        bounded.predicate_mask(
            Analytic.DegeneracyTestPlan(np.array([[1, -1]])),
            True
        )
        self.assertEqual(bounded.cache_bytes, 0)
        self.assertEqual(len(bounded.predicate_masks), 0)
        self.assertEqual(bounded_stats.degeneracy_pool_cache_skips, 1)

    def test_compiled_degeneracy_selection_matches_legacy_and_reuses_plans(self):
        import collections

        from Psience.VPT2.IndexedEvaluator import (
            ConcreteBlockTable,
            IndexedBlockEvaluator,
            IndexedDegeneracyEvaluationContext,
            IndexedEnergyTerm,
            IndexedEvaluationStats,
            IndexedSubexpressionPlan
        )

        first = Analytic.DegeneracyTestPlan(np.array([[1, -1, 3]]))
        second = Analytic.DegeneracyTestPlan(np.array([[2, -1, 4]]))
        equivalent_first = Analytic.DegeneracyTestPlan(
            np.array([[1, -1, 3]])
        )
        terms = (
            IndexedEnergyTerm(0, ('first',), object(), (1,)),
            IndexedEnergyTerm(1, ('second',), object(), (1,)),
            IndexedEnergyTerm(2, ('equivalent-first',), object(), (1,))
        )
        subexpression = IndexedSubexpressionPlan(
            object(), False, 1, (((1,), terms),)
        )
        changes = ({
            'expression': {
                ('first',): first,
                ('second',): second,
                ('equivalent-first',): equivalent_first
            }
        }, {})
        pool = np.array([
            [1, 8, 3],
            [2, 7, 4],
            [1, 0, 3],
            [9, 0, 9]
        ])
        global_indices = np.array([0, 1, 2, 3], dtype=np.intp)
        blocks = ConcreteBlockTable()
        blocks.plans.append(subexpression)
        blocks.expression_keys.append('expression')
        blocks.permutation_subsets.append(pool)
        blocks.global_permutation_indices.append(global_indices)
        mode = (True, False, np.logical_and)

        def make_backend():
            backend = object.__new__(IndexedBlockEvaluator)
            backend.stats = IndexedEvaluationStats()
            backend._degeneracy_selection_plans = collections.OrderedDict()
            backend._degeneracy_group_plans = collections.OrderedDict()
            return backend

        legacy_backend = make_backend()
        legacy_context = IndexedDegeneracyEvaluationContext(
            pool, legacy_backend.stats
        )
        legacy = legacy_backend._select_energy_groups_legacy(
            blocks, 0, changes, mode,
            degeneracy_context=legacy_context
        )

        compiled_backend = make_backend()
        compiled_context = IndexedDegeneracyEvaluationContext(
            pool, compiled_backend.stats
        )
        compiled = compiled_backend._select_energy_groups_compiled(
            blocks, 0, changes, mode,
            degeneracy_context=compiled_context
        )

        def normalize(groups):
            return [
                (
                    shape,
                    None if selected is None else tuple(selected),
                    tuple(term.ordinal for term in selected_terms)
                )
                for shape, selected, selected_terms in groups
            ]

        self.assertEqual(normalize(compiled), normalize(legacy))
        compiled_backend._select_energy_groups_compiled(
            blocks, 0, changes, mode,
            degeneracy_context=compiled_context
        )
        self.assertEqual(
            compiled_backend.stats.degeneracy_selection_plan_misses, 1
        )
        self.assertEqual(
            compiled_backend.stats.degeneracy_selection_plan_hits, 1
        )
        self.assertEqual(compiled_backend.stats.degeneracy_group_plan_misses, 1)
        self.assertEqual(compiled_backend.stats.degeneracy_group_plan_hits, 1)
        self.assertEqual(
            compiled_backend.stats.degeneracy_block_selection_misses, 1
        )
        self.assertEqual(
            compiled_backend.stats.degeneracy_block_selection_hits, 1
        )

    def test_linear_and_perfect_degeneracy_plans_match(self):
        import itertools

        evaluator_type = Analytic.PerturbationTheoryExpressionEvaluator
        solver = Analytic.AnalyticPerturbationTheorySolver.from_order(
            4,
            polynomial_representation='path'
        )
        expression = solver.energy_correction(2)([]).expr.poly_obj.to_eager()
        changes = [((0,), (-1,))]
        linear = evaluator_type._identify_possible_degeneracies(
            None, None, expression, changes, 5, method='linear'
        )
        perfect = evaluator_type._identify_possible_degeneracies(
            None, None, expression, changes, 5, method='perfect'
        )
        for linear_side, perfect_side in zip(linear, perfect):
            self.assertEqual(linear_side.keys(), perfect_side.keys())
            for coefficient_key in linear_side:
                self.assertEqual(
                    linear_side[coefficient_key].keys(),
                    perfect_side[coefficient_key].keys()
                )
                for energy_key, linear_test in linear_side[
                        coefficient_key].items():
                    perfect_test = perfect_side[coefficient_key][energy_key]
                    arity = linear_test.patterns.shape[1]
                    states = np.asarray(list(itertools.permutations(
                        range(5),
                        arity
                    )))
                    np.testing.assert_array_equal(
                        linear_test.evaluate(states),
                        perfect_test.evaluate(states)
                    )

    def test_indexed_degeneracy_identification_matches_old_plan(self):
        import itertools

        evaluator = Analytic.PerturbationTheoryExpressionEvaluator
        energy = Analytic.PTEnergyChangeProductSum({
            ((0, -1, 0, 1), (1, 1, -1, 0)): 1,
            ((0, -1, -1, 0), (1, 0, 1, -1)): 1,
            ((0, 2, -2, 0),): 1,
            ((1, 3, -3, 0),): 1
        }, canonicalize=False)
        expression = Analytic.PTTensorCoeffProductSum(
            {(): energy}, canonicalize=False, reduced=False
        )
        changes = [
            ((0, 1), (-1, 1)),
            ((2, 3), (-1, 1)),
            ((0, 1), (-1, 1)),  # duplicate input change
            ((0, 2), (-1, -1)),  # equal-quanta permutation block
            ((1, 3), (2, -2)),
            ((3, 4), (9, -9))  # does not match an energy term
        ]
        for method in ('linear', 'perfect'):
            with self.subTest(method=method):
                old = evaluator._identify_possible_degeneracies_legacy(
                    None, None, expression, changes, 5, method=method
                )
                context = Analytic.DegeneracyIdentificationContext(
                    changes, max_items=128, max_bytes=32 * 1024 ** 2
                )
                new = evaluator._identify_possible_degeneracies(
                    None, None, expression, changes, 5,
                    method=method, context=context
                )
                self.assertEqual(
                    [set(side) for side in new],
                    [set(side) for side in old]
                )
                for old_side, new_side in zip(old, new):
                    for coefficient_key, old_energy in old_side.items():
                        self.assertEqual(
                            set(old_energy), set(new_side[coefficient_key])
                        )
                        for energy_key, old_test in old_energy.items():
                            new_test = new_side[coefficient_key][energy_key]
                            arity = old_test.patterns.shape[1]
                            states = np.asarray(list(itertools.permutations(
                                range(5), arity
                            )), dtype=np.intp)
                            np.testing.assert_array_equal(
                                old_test.evaluate(states),
                                new_test.evaluate(states)
                            )
                again = context.bind(expression, 5, method)
                self.assertEqual(
                    [set(side) for side in again],
                    [set(side) for side in new]
                )
                stats = context.stats()
                self.assertGreater(stats['bound_plan_hits'], 0)
                self.assertLessEqual(stats['cache_items'], 128)
                self.assertLessEqual(stats['cache_bytes'], 32 * 1024 ** 2)
                limited = Analytic.DegeneracyIdentificationContext(
                    changes, max_items=12, max_bytes=4096
                )
                limited.bind(expression, 5, method)
                stats = limited.stats()
                self.assertLessEqual(stats['cache_bytes'], 4096)
                self.assertLessEqual(stats['cache_items'], 12)
                self.assertGreater(
                    stats.get('cache_evictions', 0)
                    + stats.get('cache_skips', 0), 0
                )

        unmatched = [((3, 4), (9, -9))]
        old = evaluator._identify_possible_degeneracies_legacy(
            None, None, expression, unmatched, 5, method='linear'
        )
        new = evaluator._identify_possible_degeneracies(
            None, None, expression, unmatched, 5, method='linear'
        )
        self.assertEqual(tuple(old), tuple(new))

    def test_compiled_identification_interns_equal_predicates(self):
        evaluator = Analytic.PerturbationTheoryExpressionEvaluator
        first = ((0, -1, 1),)
        second = ((0, -1, 1), (0, 2, -2))
        energy = Analytic.PTEnergyChangeProductSum(
            {first: 1, second: 1}, canonicalize=False
        )
        expression = Analytic.PTTensorCoeffProductSum(
            {(): energy}, canonicalize=False, reduced=False
        )
        changes = [((0, 1), (-1, 1))]
        legacy = evaluator._identify_possible_degeneracies_legacy(
            None, None, expression, changes, 4, method='linear'
        )
        context = Analytic.DegeneracyIdentificationContext(changes)
        compiled = evaluator._identify_possible_degeneracies(
            None, None, expression, changes, 4,
            method='linear', context=context
        )
        self.assertIs(legacy[0][()][first], legacy[0][()][second])
        self.assertIs(compiled[0][()][first], compiled[0][()][second])
        self.assertGreater(context.stats()['predicate_intern_content_hits'], 0)
        for key in (first, second):
            states = np.asarray([(0, 1), (1, 0), (2, 3)], dtype=np.intp)
            np.testing.assert_array_equal(
                compiled[0][()][key].evaluate(states),
                legacy[0][()][key].evaluate(states)
            )
        again = context.bind(expression, 4, 'linear')
        self.assertIs(again[0][()][first], again[0][()][second])
        self.assertGreater(context.stats()['bound_plan_hits'], 0)

    def test_degeneracy_predicate_interning_verifies_digest_and_is_bounded(self):
        context = Analytic.DegeneracyIdentificationContext(
            [], max_items=1, max_bytes=4096
        )
        first = Analytic.DegeneracyTestPlan([[0, 1]])
        second = Analytic.DegeneracyTestPlan([[1, 0]])
        by_content = {}
        by_identity = weakref.WeakKeyDictionary()
        with mock.patch.object(Analytic.hashlib, 'blake2b') as digest:
            digest.return_value.digest.return_value = b'collision'
            self.assertIs(
                context._intern_predicate(first, by_content, by_identity),
                first
            )
            self.assertIs(
                context._intern_predicate(second, by_content, by_identity),
                second
            )
        self.assertEqual(len(by_content), 1)
        self.assertEqual(context.stats()['predicate_intern_skips'], 1)

    def test_indexed_coefficient_gather_matches_legacy(self):
        from Psience.VPT2.IndexedEvaluator import (
            IndexedCoefficientGatherPlan,
            IndexedEvaluationStats
        )

        permutations = np.array([
            [0, 1, 2],
            [1, 2, 0],
            [2, 0, 1]
        ])
        coefficient_specs = [
            ((0, 0, 0), (1, 0, 1)),
            ((0, 0, 0, 1),),
            ((2, 0),)
        ]
        rank_one_a = np.array([.5, 0, 1.5])
        rank_one_b = np.array([2., 3., 4.])
        rank_two = np.arange(1, 10, dtype=float).reshape(3, 3) / 10
        coefficient_tensors = [
            [
                [rank_one_a, rank_one_b],
                [rank_one_a * 2, rank_one_b / 2]
            ],
            [
                [rank_two],
                [rank_two * 3]
            ],
            [
                [2.],
                [3.]
            ]
        ]
        subset = (2,)
        num_fixed = 1
        remapped = [
            tuple(
                (
                    coefficient[:2],
                    tuple(
                        subset[index - num_fixed]
                        if index >= num_fixed else index
                        for index in coefficient[2:]
                    )
                )
                for coefficient in product
            )
            for product in coefficient_specs
        ]
        factorials = [1, 1, 2]
        cutoff = .1
        legacy = (
            Analytic.PerturbationTheoryExpressionEvaluator
            ._get_prefacs(
                permutations,
                remapped,
                coefficient_tensors,
                {},
                factorials,
                cutoff
            )
        )
        stats = IndexedEvaluationStats()
        indexed = IndexedCoefficientGatherPlan(
            coefficient_specs,
            coefficient_tensors,
            2,
            max_workspace_bytes=128
        ).gather(
            permutations,
            subset,
            num_fixed,
            {},
            factorials,
            cutoff,
            stats
        )
        np.testing.assert_array_equal(indexed[0], legacy[0])
        np.testing.assert_array_equal(indexed[1], legacy[1])
        np.testing.assert_allclose(indexed[2], legacy[2], rtol=0, atol=0)
        self.assertGreater(stats.coefficient_product_chunks, 1)
        self.assertGreater(stats.coefficient_gather_requests, 0)

    def test_indexed_coefficient_pattern_cache_reuses_and_evicts(self):
        from Psience.VPT2.IndexedEvaluator import (
            IndexedCoefficientGatherPlan,
            IndexedEvaluationStats
        )

        plan = IndexedCoefficientGatherPlan(
            (),
            (),
            0,
            max_pattern_cache_items=2,
            max_pattern_cache_bytes=1024
        )
        patterns = [(0, 2), (0, 2), (2, 0)]
        stats = IndexedEvaluationStats()
        unique, inverse = plan._get_request_pattern_map(patterns, stats)
        repeated_unique, repeated_inverse = plan._get_request_pattern_map(
            patterns,
            stats
        )
        np.testing.assert_array_equal(repeated_unique, unique)
        np.testing.assert_array_equal(repeated_inverse, inverse)
        self.assertEqual(stats.coefficient_pattern_cache_misses, 1)
        self.assertEqual(stats.coefficient_pattern_cache_hits, 1)
        self.assertEqual(stats.coefficient_pattern_compilations, 1)

        np.testing.assert_array_equal(
            plan._remap_request_patterns(unique, (3, 4), 1),
            np.array([[0, 4], [4, 0]])
        )
        np.testing.assert_array_equal(
            plan._remap_request_patterns(unique, (5, 6), 1),
            np.array([[0, 6], [6, 0]])
        )

        plan._get_request_pattern_map([(0,), (1,)], stats)
        before_evictions = plan.request_pattern_cache.evictions
        plan._get_request_pattern_map([(1,), (2,)], stats)
        self.assertGreater(
            plan.request_pattern_cache.evictions,
            before_evictions
        )
        self.assertLessEqual(len(plan.request_pattern_cache), 2)
        self.assertLessEqual(plan.request_pattern_cache.bytes, 1024)
        self.assertLessEqual(
            plan.request_pattern_cache.peak_items,
            plan.max_pattern_cache_items
        )
        self.assertLessEqual(
            plan.request_pattern_cache.peak_bytes,
            plan.max_pattern_cache_bytes
        )

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
