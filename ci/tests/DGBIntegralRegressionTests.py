"""Check polynomial DGB integrals against moments and direct quadrature."""
import itertools
import math
import unittest

import numpy as np
from Psience.DGB import DGB
from Psience.DGB.Evaluators import DGBEvaluator


def monomial(powers):
    powers = np.asarray(powers)
    def evaluate(coords, deriv_order=None):
        coords = np.asarray(coords)
        value = np.prod(coords**powers, axis=-1)
        if deriv_order is None:
            return value
        derivatives = [value]
        for order in range(1, deriv_order + 1):
            tensor = np.zeros(coords.shape[:-1] + (len(powers),) * order)
            for indices in itertools.product(range(len(powers)), repeat=order):
                counts = np.bincount(indices, minlength=len(powers))
                if np.any(counts > powers):
                    continue
                coefficient = math.prod(math.factorial(int(p)) // math.factorial(int(p - k))
                                        for p, k in zip(powers, counts))
                tensor[(...,) + indices] = coefficient * np.prod(coords**(powers - counts), axis=-1)
            derivatives.append(tensor)
        return derivatives
    return evaluate


class DGBIntegralRegressionTests(unittest.TestCase):
    def test_momentum_moments_through_order_eight(self):
        # Completing the square gives a complex Gaussian with mean i*p/(2*a).
        # Its moment recurrence supplies an independent polynomial oracle.
        p = np.array([-.8, 0., .3, 1.1])
        a = np.array([.6, 1., 1.4, 2.])
        mean = 1j * p / (2 * a)
        variance = 1 / (2 * a)
        moments = [np.ones_like(mean), mean]
        for order in range(2, 9):
            moments.append(mean * moments[-1] + (order - 1) * variance * moments[-2])
        for order, moment in enumerate(moments):
            with self.subTest(order=order):
                expected = moment.imag if order % 2 else moment.real
                np.testing.assert_allclose(DGBEvaluator.momentum_integral(p, a, order),
                                           expected, atol=1e-12, rtol=1e-12)

    def test_mixed_quartic_normalization(self):
        dgb = DGB.construct(np.zeros((1, 2)), monomial([2, 2]), alphas=1,
                            masses=[1], expansion_degree=4, logger=False)
        # |psi|^2 has independent coordinate variances 1/(4*alpha).
        self.assertAlmostEqual(dgb.V[0, 0], 1 / 16, places=13)

    def check_diagonal(self, powers, phased):
        options = {'momenta': np.array([[.5, .7]])} if phased else {}
        dgb = DGB.construct(np.zeros((1, 2)), monomial(powers), alphas=1,
                            masses=[1], expansion_degree=sum(powers), logger=False, **options)
        expansion = dgb.V.copy()
        dgb.pot.expansion_degree = None
        dgb.pot.quadrature_degree = 15 if phased else 4
        dgb._V = None
        np.testing.assert_allclose(expansion, dgb.V, atol=1e-11, rtol=1e-11)

    def test_unphased_diagonal_moments(self):
        for powers in ([4, 0], [2, 2], [2, 4]):
            with self.subTest(powers=powers):
                self.check_diagonal(powers, False)

    def test_phased_diagonal_moments(self):
        for powers in ([1, 1], [4, 0], [3, 1], [2, 2], [3, 3]):
            with self.subTest(powers=powers):
                self.check_diagonal(powers, True)

    def check_matrix(self, phased):
        centers = np.array([[0, 0], [.2, -.3], [-.3, .1], [.4, .3]])
        alphas = np.array([[1, .7], [.8, 1.3], [1.2, .9], [.6, 1.1]])
        angles = [.2, .4, -.3, .5]
        rotations = np.array([[[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]] for t in angles])
        options = {'momenta': np.array([[.5, .7], [.6, -.2], [-.2, .3], [.2, -.5]])} if phased else {}
        for powers in ([1, 1], [2, 2], [4, 0], [3, 1], [2, 4]):
            with self.subTest(powers=powers):
                dgb = DGB.construct(centers, monomial(powers), alphas=alphas,
                                    transformations=rotations, masses=[1],
                                    expansion_degree=sum(powers), logger=False, **options)
                expansion = dgb.V.copy()
                dgb.pot.expansion_degree = None
                dgb.pot.quadrature_degree = 15 if phased else 4
                dgb._V = None
                np.testing.assert_allclose(expansion, dgb.V, atol=1e-10, rtol=1e-10)

    def test_unphased_displaced_rotated_unequal_width_matrices(self):
        self.check_matrix(False)

    def test_phased_displaced_rotated_unequal_width_matrices(self):
        self.check_matrix(True)


if __name__ == '__main__':
    unittest.main()
