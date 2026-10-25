"""Compare grouped DGB integration with the original path and independent quadrature."""
import itertools
import math
import unittest
import numpy as np
from McUtils.Scaffolding import NullLogger
from McUtils.Zachary import SparsePolynomial
from Psience.DGB import DGB
from Psience.DGB.Evaluators import DGBEvaluator, DGBPotentialEnergyEvaluator


TERMS = {(0, 0): .2, (2, 0): .5, (0, 2): .7, (1, 1): .1,
         (4, 0): .12, (3, 1): -.05, (2, 2): .08, (8, 0): 1e-4,
         (4, 4): 2e-4, (16, 0): 1e-10, (8, 8): 2e-10}


def polynomial(terms):
    def evaluate(coords, deriv_order=None):
        coords = np.asarray(coords)
        derivatives = []
        for order in range((0 if deriv_order is None else deriv_order) + 1):
            # For symmetric bivariate derivatives, only order+1 unique values are needed.
            values = np.zeros(coords.shape[:-1] + (order + 1,), dtype=np.result_type(*terms.values(), float))
            for ny in range(order + 1):
                nx = order - ny
                for (a, b), coefficient in terms.items():
                    if nx <= a and ny <= b:
                        factor = math.factorial(a) / math.factorial(a - nx) * math.factorial(b) / math.factorial(b - ny)
                        values[..., ny] += coefficient * factor * coords[..., 0]**(a - nx) * coords[..., 1]**(b - ny)
            if order:
                ny = np.array([sum(slots) for slots in itertools.product(range(2), repeat=order)])
                derivatives.append(values[..., ny].reshape(coords.shape[:-1] + (2,) * order))
            else:
                derivatives.append(values[..., 0])
        return derivatives[0] if deriv_order is None else derivatives
    return evaluate


def overlap(phased=False):
    centers = np.array([[.2, -.3], [.5, .1], [-.4, .4]])
    angles = [.2, -.6, .9]
    rotations = np.array([[[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]] for t in angles])
    alphas = np.array([[1., 1.2], [.8, 1.4], [1.3, .9]])
    momenta = np.array([[.5, .7], [.6, -.2], [-.2, .3]]) if phased else None
    return DGBEvaluator.get_overlap_gaussians(centers, alphas, (rotations, rotations.transpose(0, 2, 1)),
                                              momenta, logger=NullLogger())


def sparse(terms):
    return SparsePolynomial({(0,) * a + (1,) * b: c for (a, b), c in terms.items()}, ndim=2)


class DGBGroupedIntegralTests(unittest.TestCase):
    def test_orders_two_four_eight_sixteen(self):
        od = overlap()
        derivatives = polynomial(TERMS)(od.centers, deriv_order=16)
        copies = [d.copy() for d in derivatives]
        for degree in (2, 4, 8, 16):
            with self.subTest(degree=degree):
                direct = DGBEvaluator.tensor_expansion_integrate(3, derivatives[:degree + 1], od, method='direct')
                default = DGBEvaluator.tensor_expansion_integrate(3, derivatives[:degree + 1], od)
                grouped = DGBEvaluator.tensor_expansion_integrate(3, derivatives[:degree + 1], od, method='grouped')
                cutoff = DGBEvaluator.tensor_expansion_integrate(3, derivatives, od, method='grouped', expansion_degree=degree)
                global_poly = DGBEvaluator.tensor_expansion_integrate(3, sparse(TERMS), od, method='grouped',
                                                                     expansion_type='taylor', expansion_degree=degree)
                np.testing.assert_array_equal(default, direct)
                np.testing.assert_allclose(grouped, direct, atol=2e-12, rtol=2e-12)
                np.testing.assert_allclose(cutoff, direct, atol=2e-12, rtol=2e-12)
                np.testing.assert_allclose(global_poly, direct, atol=2e-12, rtol=2e-12)
        for original, derivative in zip(copies, derivatives):
            np.testing.assert_array_equal(original, derivative)

    def test_phased_partial_expansions_and_quadrature(self):
        terms = {p: c for p, c in TERMS.items() if sum(p) <= 8}
        od = overlap(True)
        derivatives = polynomial(terms)(od.centers, deriv_order=8)
        for degree in (0, 1, 2, 3, 4, 8):
            with self.subTest(degree=degree):
                direct = DGBEvaluator.tensor_expansion_integrate(3, derivatives[:degree + 1], od)
                grouped = DGBEvaluator.tensor_expansion_integrate(3, derivatives[:degree + 1], od, method='grouped')
                global_poly = DGBEvaluator.tensor_expansion_integrate(3, sparse(terms), od, method='grouped',
                                                                     expansion_type='taylor', expansion_degree=degree)
                np.testing.assert_allclose(grouped, direct, atol=2e-12)
                np.testing.assert_allclose(global_poly, direct, atol=2e-12)
        quadrature = DGBEvaluator.quad_integrate(polynomial(terms), od, degree=22, logger=NullLogger())
        np.testing.assert_allclose(grouped, quadrature, atol=2e-11, rtol=2e-11)

    def test_vector_matrix_and_complex_operators(self):
        od = overlap(True)
        terms = {(0, 0): .2, (3, 1): -.05, (2, 2): .08}
        derivatives = polynomial(terms)(od.centers, deriv_order=4)
        reference = DGBEvaluator.tensor_expansion_integrate(3, derivatives, od)
        for scale in (np.array([1., -.5]), np.array([[1 + 2j, -.5j], [.3, -.2 + .7j]])):
            expanded = [d.reshape((len(d),) + (1,) * scale.ndim + d.shape[1:])
                        * scale.reshape((1,) + scale.shape + (1,) * order)
                        for order, d in enumerate(derivatives)]
            grouped = DGBEvaluator.tensor_expansion_integrate(3, expanded, od, method='grouped')
            np.testing.assert_allclose(grouped, reference.reshape(reference.shape + (1,) * scale.ndim) * scale,
                                       atol=2e-12)
            poly = sparse({p: c * scale[None] for p, c in terms.items()})
            grouped_sparse = DGBEvaluator.tensor_expansion_integrate(3, poly, od, method='grouped', expansion_type='taylor')
            np.testing.assert_allclose(grouped_sparse, grouped, atol=2e-12)

    def test_sparse_local_polynomial_pair_coefficients_and_scaling(self):
        for phased in (False, True):
            od = overlap(phased)
            coefficient = np.arange(len(od.centers)) / 4 + 1
            poly = SparsePolynomial({(): coefficient, (1,): .3 * coefficient, (0, 0): .2 * coefficient}, prefactor=2)
            derivatives = [2 * coefficient, np.stack([np.zeros_like(coefficient), .6 * coefficient], axis=-1),
                           np.zeros((len(coefficient), 2, 2))]
            derivatives[2][:, 0, 0] = .8 * coefficient
            for degree in (1, 2):
                reference = DGBEvaluator.tensor_expansion_integrate(3, derivatives[:degree + 1], od)
                grouped = DGBEvaluator.tensor_expansion_integrate(3, poly, od, method='grouped', expansion_degree=degree)
                np.testing.assert_allclose(grouped, reference, atol=1e-12)

    def test_global_derivative_tensors_at_origin(self):
        od = overlap(True)
        terms = {(0, 0): .2, (4, 0): .5, (3, 1): -.1}
        derivatives = polynomial(terms)(np.zeros_like(od.centers), deriv_order=4)
        for degree in (1, 2, 4):
            grouped = DGBEvaluator.tensor_expansion_integrate(3, derivatives, od, method='grouped',
                                                             expansion_type='taylor', expansion_degree=degree)
            local = polynomial(terms)(od.centers, deriv_order=degree)
            direct = DGBEvaluator.tensor_expansion_integrate(3, local, od)
            np.testing.assert_allclose(grouped, direct, atol=1e-12)

    def test_unweighted_legacy_rotated_coefficients(self):
        for phased in (False, True):
            od = overlap(phased)
            derivatives = polynomial({(2, 2): .3, (3, 1): -.2, (1, 0): .4})(od.centers, deriv_order=4)
            direct = DGBEvaluator.tensor_expansion_integrate(3, derivatives, od, reweight=False)
            grouped = DGBEvaluator.tensor_expansion_integrate(3, derivatives, od, reweight=False, method='grouped')
            np.testing.assert_allclose(grouped, direct, atol=2e-12)

    def test_six_dimensional_derivative_tensors(self):
        rng = np.random.default_rng(981)
        centers = rng.normal(size=(3, 6)) / 3
        rotations = np.linalg.qr(rng.normal(size=(3, 6, 6)))[0]
        od = DGBEvaluator.get_overlap_gaussians(centers, rng.uniform(.5, 2, (3, 6)),
                                                (rotations, rotations.transpose(0, 2, 1)), None, logger=NullLogger())
        # (v . delta)^4 produces a symmetric rank-four tensor coupling all axes.
        vector = rng.normal(size=6)
        fourth = vector[:, None, None, None] * vector[None, :, None, None] * vector[None, None, :, None] * vector[None, None, None, :]
        derivatives = [np.zeros(len(od.centers)), 0, 0, 0, np.broadcast_to(24 * fourth, (len(od.centers),) + fourth.shape)]
        grouped = DGBEvaluator.tensor_expansion_integrate(3, derivatives, od, method='grouped')
        variance = np.sum(vector * (od.covariance_matrices @ vector), axis=-1)
        np.testing.assert_allclose(grouped[od.indices], 3 * variance**2, atol=1e-11)
        dense = [np.zeros((len(od.centers),) + (6,) * n) if np.ndim(d) == 0 else d for n, d in enumerate(derivatives)]
        np.testing.assert_allclose(grouped, DGBEvaluator.tensor_expansion_integrate(3, dense, od), atol=1e-11)

    def test_factory_potential_and_operator_dispatch(self):
        terms = {(2, 0): .5, (0, 2): .7, (2, 2): .08}
        dgb = DGB.construct(np.array([[0., 0.], [.2, -.3]]), polynomial(terms), alphas=1,
                            masses=[1], expansion_degree=4, expansion_method='grouped', logger=False)
        grouped = dgb.V.copy()
        od = dgb.gaussians.overlap_data
        dgb.pot.expansion_method = 'direct'; dgb._V = None
        np.testing.assert_allclose(grouped, dgb.V, atol=1e-13)
        sparse_value = DGBPotentialEnergyEvaluator.expansion_integrate(sparse(terms), od, 'multicenter',
                                                                      expansion_degree=4, method='grouped')
        dense_value = dgb.pot.evaluate_op(polynomial(terms), od, expansion_method='grouped')
        np.testing.assert_allclose(sparse_value, dense_value, atol=1e-13)

    def test_zero_placeholders_constants_and_invalid_method(self):
        od = overlap()
        derivative = polynomial({(0, 0): 2., (2, 0): .5})(od.centers, deriv_order=2)
        derivative[1] = 0
        grouped = DGBEvaluator.tensor_expansion_integrate(3, derivative, od, method='grouped')
        expected = DGBEvaluator.tensor_expansion_integrate(3, polynomial({(0, 0): 2., (2, 0): .5})(od.centers, deriv_order=2), od)
        np.testing.assert_allclose(grouped, expected)
        np.testing.assert_array_equal(DGBEvaluator.tensor_expansion_integrate(3, SparsePolynomial({}, ndim=2), od,
                                                                            method='grouped'), np.zeros((3, 3)))
        with self.assertRaises(ValueError):
            DGBEvaluator.tensor_expansion_integrate(3, derivative, od, method='unknown')
        with self.assertRaises(ValueError):
            DGBEvaluator.tensor_expansion_integrate(3, [], od, method='grouped', expansion_degree=-1)


if __name__ == '__main__':
    unittest.main()
