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
        with tempfile.TemporaryDirectory() as tmpdir:
            checkpoint = os.path.join(tmpdir, 'analytic-path.hdf5')
            solver = Analytic.AnalyticPerturbationTheorySolver.from_order(
                4,
                checkpoint=checkpoint,
                polynomial_representation='path'
            )
            original = solver.energy_correction(2)([])

            restored_solver = Analytic.AnalyticPerturbationTheorySolver.from_order(
                4,
                checkpoint=checkpoint,
                polynomial_representation='path'
            )
            restored = restored_solver.energy_correction(2)([])
            self.assertTrue(Analytic.polynomial_uses_path(restored.expr))
            self.assertPolynomialEqual(original.expr, restored.expr)


if __name__ == '__main__':
    unittest.main()
