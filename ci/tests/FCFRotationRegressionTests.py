"""Independent harmonic-integral checks for Duschinsky FCF amplitudes."""

from itertools import product
from math import factorial, pi, sqrt
from unittest import TestCase

import numpy as np
from scipy.special import eval_hermite

from McUtils.Coordinerds import CartesianCoordinates3D
from Psience.Modes import NormalModes
from Psience.Vibronic import FranckCondonModel


class FCFRotationRegressionTests(TestCase):
    @staticmethod
    def direct_overlap(wg, we, mg, me, cg, ce, ng, ne):
        """Integrate the two normalized Cartesian wavefunctions by Gauss-Hermite quadrature."""
        zg = mg @ np.diag(wg) @ mg.T
        ze = me @ np.diag(we) @ me.T
        z = zg + ze
        covariance = np.linalg.inv(z)
        mean = covariance @ (zg @ cg + ze @ ce)
        constant = cg @ zg @ cg + ce @ ze @ ce - mean @ z @ mean
        baseline = (2 ** (len(wg) / 2) * np.prod(wg * we) ** .25
                    / sqrt(np.linalg.det(z)) * np.exp(-constant / 2))

        # Five points per dimension exactly integrate the degree <= 4 polynomials here.
        nodes, weights = np.polynomial.hermite_e.hermegauss(5)
        root = np.linalg.cholesky(covariance)
        expectation = 0.0
        for ids in product(range(5), repeat=len(wg)):
            x = mean + root @ nodes[list(ids)]
            qg = mg.T @ (x - cg)
            qe = me.T @ (x - ce)
            polynomial = 1.0
            for n, w, q in zip(ng, wg, qg):
                polynomial *= eval_hermite(n, sqrt(w) * q) / sqrt(2 ** n * factorial(n))
            for n, w, q in zip(ne, we, qe):
                polynomial *= eval_hermite(n, sqrt(w) * q) / sqrt(2 ** n * factorial(n))
            expectation += np.prod(weights[list(ids)]) * polynomial
        return baseline * expectation / (2 * pi) ** (len(wg) / 2)

    def test_rotated_signed_amplitudes_and_squares(self):
        wg = np.array([1.0, 2.0, 1.6])
        we = np.array([3.0, 4.0, 2.2])
        c, s = np.cos(.5), np.sin(.5)
        c2, s2 = np.cos(.3), np.sin(.3)
        mg = np.eye(3)
        me = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]]) @ np.array(
            [[c2, 0, -s2], [0, 1, 0], [s2, 0, c2]])
        cg = np.zeros(3)
        ce = np.array([1.0, .2, -.35])
        ground_states = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
        excited_states = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [2, 0, 0]]
        expected = np.array([
            self.direct_overlap(wg, we, mg, me, cg, ce, ng, ne)
            for ng in ground_states for ne in excited_states
        ])
        for order in ('gs', 'es'):
            for center in ('gs', 'es'):
                with self.subTest(order=order, center=center):
                    actual = np.asarray(FranckCondonModel.eval_fcf_overlaps(
                        ground_states, wg, mg, mg.T, cg,
                        excited_states, we, me, me.T, ce,
                        rotation_order=order, rotation_center=center,
                    ))
                    np.testing.assert_allclose(actual, expected, rtol=1e-10, atol=1e-12)
                    np.testing.assert_allclose(actual ** 2, expected ** 2, rtol=1e-10, atol=1e-12)

    def test_two_mode_zero_zero_reference(self):
        wg, we = np.array([1.0, 2.0]), np.array([3.0, 4.0])
        c, s = np.cos(.5), np.sin(.5)
        mg, me = np.eye(2), np.array([[c, -s], [s, c]])
        cg, ce = np.zeros(2), np.array([1.0, .2])
        expected = self.direct_overlap(wg, we, mg, me, cg, ce, [0, 0], [0, 0])
        self.assertAlmostEqual(expected, 0.6028123095072196, places=12)
        actual = FranckCondonModel.eval_fcf_overlaps(
            [[0, 0]], wg, mg, mg.T, cg,
            [[0, 0]], we, me, me.T, ce,
        )[0]
        self.assertAlmostEqual(actual, expected, places=12)
        self.assertAlmostEqual(actual ** 2, 0.363382680493428, places=12)

        gs = NormalModes(CartesianCoordinates3D, mg, inverse=mg.T,
                         freqs=wg, origin=cg)
        es = NormalModes(CartesianCoordinates3D, me, inverse=me.T,
                         freqs=we, origin=ce)
        model = FranckCondonModel(gs, es, embed=False, mass_weight=False)
        states = [[0, 0], [1, 0]]
        amplitudes = np.asarray(model.get_overlaps(
            states, return_states=False, embed=False, mass_weight=False))
        self.assertLess(amplitudes[1], 0)
        np.testing.assert_allclose(amplitudes, [
            self.direct_overlap(wg, we, mg, me, cg, ce, [0, 0], state)
            for state in states
        ], rtol=1e-10, atol=1e-12)
        spectrum = model.get_spectrum(states)
        np.testing.assert_allclose(spectrum.intensities, amplitudes ** 2,
                                   rtol=1e-10, atol=1e-12)
