"""Check partition assignments independently of the FCF combinatorial implementation."""
import itertools
import math
import unittest
import numpy as np
from scipy.special import eval_hermite
from McUtils.Zachary import DensePolynomial
from Psience.Vibronic.FCFs import FranckCondonModel, HermiteProductPolynomial


def monomial_poly(powers):
    factors = {}
    for coordinate, power in enumerate(powers):
        coefficients = np.zeros(power + 1); coefficients[-1] = 1
        factors[coordinate] = DensePolynomial(coefficients)
    return HermiteProductPolynomial(factors, ndim=len(powers))


def monomial_quadrature(powers, Q, alphas):
    degree = 1 + sum(powers) // 2
    nodes, weights = np.polynomial.hermite.hermgauss(degree)
    indices = np.array(list(itertools.product(range(degree), repeat=len(alphas))))
    coordinates = np.sqrt(2 / alphas) * nodes[indices]
    values = np.prod((coordinates @ Q) ** np.asarray(powers), axis=-1)
    return np.sum(np.prod(weights[indices] / np.sqrt(np.pi), axis=-1) * values)


def harmonic_quadrature(wg, we, mg, me, cg, ce, ground_states, excited_states):
    zg, ze = mg @ np.diag(wg) @ mg.T, me @ np.diag(we) @ me.T
    precision = zg + ze
    covariance = np.linalg.inv(precision)
    mean = covariance @ (zg @ cg + ze @ ce)
    constant = cg @ zg @ cg + ce @ ze @ ce - mean @ precision @ mean
    baseline = (2**(len(wg) / 2) * np.prod(wg * we)**.25
                / np.sqrt(np.linalg.det(precision)) * np.exp(-constant / 2))
    degree = 1 + (max(map(sum, ground_states)) + max(map(sum, excited_states))) // 2
    nodes, weights = np.polynomial.hermite.hermgauss(degree)
    indices = np.array(list(itertools.product(range(degree), repeat=len(wg))))
    coords = mean + np.sqrt(2) * nodes[indices] @ np.linalg.cholesky(covariance).T
    qg, qe = (coords - cg) @ mg, (coords - ce) @ me
    weight = np.prod(weights[indices] / np.sqrt(np.pi), axis=-1)
    values = []
    for ng in ground_states:
        for ne in excited_states:
            value = np.ones(len(coords))
            for quanta, frequency, coordinate in zip(ng, wg, qg.T):
                value *= eval_hermite(quanta, np.sqrt(frequency) * coordinate) / np.sqrt(2**quanta * math.factorial(quanta))
            for quanta, frequency, coordinate in zip(ne, we, qe.T):
                value *= eval_hermite(quanta, np.sqrt(frequency) * coordinate) / np.sqrt(2**quanta * math.factorial(quanta))
            values.append(baseline * np.sum(weight * value))
    return np.array(values)


class FCFPartitionRegressionTests(unittest.TestCase):
    def test_sixth_moment_missing_reverse_assignment(self):
        Q = np.ones((3, 1)) / np.sqrt(3)
        observed = FranckCondonModel.evaluate_shifted_poly_overlap(monomial_poly([6]), Q, np.ones(3), 1.)
        self.assertAlmostEqual(observed, 15., places=12)

    def test_linear_form_moments_through_order_twelve(self):
        rng = np.random.default_rng(79214)
        for ndim in (1, 2, 3, 5):
            alphas = rng.uniform(.7, 1.8, size=ndim)
            Q = rng.normal(size=(ndim, 1)) / np.sqrt(ndim)
            variance = np.sum(Q[:, 0]**2 / alphas)
            for order in (2, 4, 6, 8, 10, 12):
                expected = math.prod(range(1, order, 2)) * variance**(order // 2)
                with self.subTest(dimensions=ndim, order=order):
                    observed = FranckCondonModel.evaluate_shifted_poly_overlap(monomial_poly([order]), Q, alphas, 1.)
                    np.testing.assert_allclose(observed, expected, atol=2e-11, rtol=2e-12)

    def test_mixed_powers_against_quadrature(self):
        rng = np.random.default_rng(79215)
        for powers in ([4, 2], [3, 3], [4, 2, 2], [2, 2, 2, 2], [6, 2], [5, 3, 2]):
            Q = rng.normal(size=(3, len(powers))) / np.sqrt(3)
            alphas = np.array([.7, 1.3, 1.9])
            expected = monomial_quadrature(powers, Q, alphas)
            with self.subTest(powers=powers):
                observed = FranckCondonModel.evaluate_shifted_poly_overlap(monomial_poly(powers), Q, alphas, 1.)
                np.testing.assert_allclose(observed, expected, atol=1e-11, rtol=2e-12)

    def test_coordinate_reordering_preserves_amplitudes(self):
        rng = np.random.default_rng(79216)
        Q = rng.normal(size=(4, 3)) / 2
        alphas = np.array([.7, 1.1, 1.5, 1.9])
        poly = monomial_poly([4, 2, 2])
        expected = monomial_quadrature([4, 2, 2], Q, alphas)
        for permutation in itertools.permutations(range(4)):
            indices = np.array(permutation)
            observed = FranckCondonModel.evaluate_shifted_poly_overlap(poly, Q[indices], alphas[indices], 1.)
            np.testing.assert_allclose(observed, expected, atol=1e-11, rtol=2e-12)

    def test_repeated_powers_are_not_counted_twice(self):
        # Independent even powers each have a single coordinate assignment.
        for powers in ([2, 2], [2, 2, 2], [4, 4]):
            Q = np.eye(len(powers))
            expected = math.prod(math.prod(range(1, order, 2)) for order in powers)
            observed = FranckCondonModel.evaluate_shifted_poly_overlap(monomial_poly(powers), Q, np.ones(len(powers)), 1.)
            self.assertAlmostEqual(observed, expected, places=12)

    def test_batched_and_precomputed_plans(self):
        rng = np.random.default_rng(79217)
        alphas = np.array([.8, .9, 1.1, 1.4, 1.9])
        matrices = [rng.normal(size=(5, 1)) / 2 for _ in range(2)]
        exponentiated = [Q[:, :, None] ** np.arange(9) for Q in matrices]
        baseline = .37
        expected = np.array([105 * np.sum(Q[:, 0]**2 / alphas)**4 * baseline for Q in matrices])
        for cached in (False, True):
            plan = FranckCondonModel.get_poly_evaluation_plan(
                (8,), alphas=alphas if cached else None, zpe_prod=baseline if cached else None)
            for batch_size in (1, 2, 3, 4, 6, 10, 10**7):
                with self.subTest(cached=cached, batch_size=batch_size):
                    observed = plan(alphas, baseline, exponentiated, batch_size=batch_size)
                    np.testing.assert_allclose(observed, expected, atol=2e-11, rtol=2e-12)

    def test_rotated_high_quanta_and_hot_bands(self):
        rng = np.random.default_rng(79218)
        wg, we = np.array([.9, 1.3, 1.8]), np.array([1.1, 1.5, 2.])
        mg, me = np.eye(3), np.linalg.qr(rng.normal(size=(3, 3)))[0]
        cg, ce = np.zeros(3), np.array([.3, -.2, .1])
        ground = [[0, 0, 0], [1, 1, 0], [2, 0, 1]]
        excited = [[0, 0, 0], [6, 0, 0], [3, 2, 1], [4, 0, 0], [1, 1, 2], [4, 2, 2]]
        expected = harmonic_quadrature(wg, we, mg, me, cg, ce, ground, excited)
        for order in ('gs', 'es'):
            for center in ('gs', 'es'):
                with self.subTest(rotation_order=order, rotation_center=center):
                    observed = np.asarray(FranckCondonModel.eval_fcf_overlaps(
                        ground, wg, mg, mg.T, cg, excited, we, me, me.T, ce,
                        rotation_order=order, rotation_center=center))
                    np.testing.assert_allclose(observed, expected, atol=2e-11, rtol=2e-10)
                    np.testing.assert_allclose(observed**2, expected**2, atol=2e-12, rtol=4e-10)


if __name__ == '__main__':
    unittest.main()
