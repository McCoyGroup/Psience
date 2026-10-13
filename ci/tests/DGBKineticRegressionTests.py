"""DGB kinetic integrals checked against direct Gaussian gradient quadrature."""
import itertools
import unittest

import numpy as np

from McUtils.Scaffolding import NullLogger
from Psience.DGB import DGB, DGBWavefunctions
from Psience.DGB.Evaluators import DGBEvaluator, DGBWatsonEvaluator


class DGBKineticRegressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.centers = np.array([[0., 0., 0.], [.2, -.3, .1],
                                [-.3, .1, -.2], [.4, .3, -.1]])
        cls.alphas = np.array([[1., .7, .8], [.8, 1.3, 1.1],
                               [1.2, .9, .6], [.6, 1.1, .9]])
        cls.rotations = np.array([
            [[np.cos(t), -np.sin(t), 0], [np.sin(t), np.cos(t), 0], [0, 0, 1]]
            for t in (.2, .4, -.3, .5)
        ])
        cls.momenta = np.array([[0., 0., 0.], [.6, -.2, .3],
                                [-.2, .3, .4], [.2, -.5, -.1]])
        cls.zeta = np.array([[[0., .8, -.3], [-.8, 0., .2], [.3, -.2, 0.]],
                             [[0., -.2, .7], [.2, 0., -.4], [-.7, .4, 0.]]])
        cls.rotational_constants = np.array([.1, .07])
        cls.tensor = np.tensordot(cls.rotational_constants[:, None, None]*cls.zeta,
                                  cls.zeta, axes=[0, 0])
        nodes, weights = np.polynomial.hermite.hermgauss(13)
        slots = np.array(list(itertools.product(range(13), repeat=3)))
        cls.nodes = nodes[slots]
        cls.weights = np.prod(weights[slots], axis=-1) / np.pi**1.5

    def parameters(self, phased, rotated=True):
        rotations = (self.rotations if rotated else
                     np.broadcast_to(np.eye(3), self.rotations.shape))
        precision = (rotations*self.alphas[:, None, :]) @ rotations.transpose(0, 2, 1)
        phases = ((rotations @ self.momenta[:, :, None])[:, :, 0] if phased else
                  np.zeros_like(self.centers))
        normalization = np.prod((2*self.alphas/np.pi)**.25, axis=-1)
        if phased:
            phase = np.sum(phases*self.centers, axis=-1)
            decay = np.sum(self.momenta**2/(2*self.alphas), axis=-1)
            normalization /= np.sqrt((1+np.exp(-decay)*np.cos(2*phase))/2)
        return rotations, precision, phases, normalization

    def quadrature(self, i, j, precision):
        # Form product means and quadrature axes independently of overlap_data.
        total = precision[i]+precision[j]
        mean = np.linalg.solve(total, precision[i]@self.centers[i] +
                               precision[j]@self.centers[j])
        values, vectors = np.linalg.eigh(total)
        q = mean + (self.nodes/np.sqrt(values))@vectors.T
        density = np.prod(np.sqrt(values/np.pi))*np.exp(-np.sum(self.nodes**2, axis=-1))
        return q, self.weights/density

    def basis_values_and_gradients(self, q, i, precision, phases, normalization):
        delta = q-self.centers[i]
        envelope = normalization[i]*np.exp(-np.sum((delta@precision[i])*delta, axis=-1))
        phase = q@phases[i]
        values = envelope*np.cos(phase)
        gradients = envelope[:, None]*(-2*(delta@precision[i])*np.cos(phase)[:, None]
                                        - phases[i]*np.sin(phase)[:, None])
        return values, gradients

    def check_kinetic_integrals(self, phased, rotated=True):
        rotations, precision, phases, normalization = self.parameters(phased, rotated)
        transforms = (rotations, rotations.transpose(0, 2, 1))
        overlap = DGBEvaluator.get_overlap_gaussians(
            self.centers, self.alphas, transforms, self.momenta if phased else None,
            logger=NullLogger())
        prefactor, S = DGBEvaluator.evaluate_overlap(
            overlap, logger=NullLogger(), return_prefactor=True)
        expected_S = np.zeros_like(S)
        expected_T = np.zeros_like(S)
        expected_C = np.zeros_like(S)
        for i, j in zip(*overlap.indices):
            q, weights = self.quadrature(i, j, precision)
            phi_i, grad_i = self.basis_values_and_gradients(q, i, precision, phases, normalization)
            phi_j, grad_j = self.basis_values_and_gradients(q, j, precision, phases, normalization)
            directions = np.tensordot(q, self.zeta, axes=[1, 1])
            angular_i = np.sum(directions*grad_i[:, None, :], axis=-1)
            angular_j = np.sum(directions*grad_j[:, None, :], axis=-1)
            expected_S[i, j] = weights@(phi_i*phi_j)
            expected_T[i, j] = weights@np.sum(grad_i*grad_j, axis=-1)/2
            expected_C[i, j] = weights@np.sum(
                angular_i*angular_j*self.rotational_constants, axis=-1)
        rows, cols = overlap.indices
        for matrix in (expected_S, expected_T, expected_C):
            matrix[cols, rows] = matrix[rows, cols]
        C = prefactor*DGBWatsonEvaluator.evaluate_coriolis_contrib(
            np.broadcast_to(self.tensor, (len(rows),)+self.tensor.shape), overlap)
        T = prefactor*DGBWatsonEvaluator.evaluate_diagonal_rotated_momentum_contrib(
            overlap, np.ones(3))
        U = prefactor*DGBWatsonEvaluator.evaluate_watson_term(
            np.full(len(rows), np.sum(self.rotational_constants)), overlap)
        np.testing.assert_allclose(S, expected_S, atol=2e-12, rtol=2e-12)
        np.testing.assert_allclose(T, expected_T, atol=2e-12, rtol=2e-12)
        np.testing.assert_allclose(C, expected_C, atol=2e-12, rtol=2e-12)
        np.testing.assert_allclose(U, -np.sum(self.rotational_constants)*expected_S/4,
                                   atol=2e-12, rtol=2e-12)

    def test_unphased_kinetic_gradient_quadrature(self):
        self.check_kinetic_integrals(False)

    def test_phased_kinetic_gradient_quadrature(self):
        self.check_kinetic_integrals(True)

    def test_phased_unrotated_kinetic_gradient_quadrature(self):
        self.check_kinetic_integrals(True, rotated=False)

    def wavefunctions(self, phased):
        # construct accepts momenta in the coordinate frame; the overlap
        # evaluator accepts momenta already aligned with each Gaussian.
        momenta = self.parameters(True)[2] if phased else None
        dgb = DGB.construct(self.centers, lambda q: np.sum(q*q, axis=-1)/2,
                            alphas=self.alphas, masses=[1], transformations=self.rotations,
                            momenta=momenta, logger=False)
        return DGBWavefunctions(np.arange(4.), np.eye(4), hamiltonian=dgb)

    def test_rotated_phased_basis_evaluation(self):
        states = self.wavefunctions(True)
        _, precision, phases, normalization = self.parameters(True)
        for i in range(4):
            q, weights = self.quadrature(i, i, precision)
            expected, _ = self.basis_values_and_gradients(q, i, precision, phases, normalization)
            actual = states[i].evaluate(q)
            np.testing.assert_allclose(actual, expected, atol=2e-12, rtol=2e-12)
            self.assertAlmostEqual(float(weights@(actual*actual)), 1., places=11)

    def test_multidimensional_evaluation_grid(self):
        q = np.arange(36, dtype=float).reshape(3, 4, 3)/20-.5
        for phased in (False, True):
            with self.subTest(phased=phased):
                states = self.wavefunctions(phased)
                expected = states[1].evaluate(q.reshape(-1, 3)).reshape(3, 4)
                np.testing.assert_allclose(states[1].evaluate(q), expected, atol=1e-13)


if __name__ == '__main__':
    unittest.main()
