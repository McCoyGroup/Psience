import unittest
import numpy as np
from scipy.linalg import eigh
from Psience.DGB import DGB


class DGBGeneralizedEigenproblemTests(unittest.TestCase):
    def test_rotated_harmonic_solver_modes(self):
        centers = np.random.default_rng(881).uniform(-2, 2, (32, 2))
        theta = .37
        rotation = np.array([[np.cos(theta), -np.sin(theta)],
                             [np.sin(theta), np.cos(theta)]])
        potential = lambda q, deriv_order=None: np.sum(q**2, axis=-1)/2
        ham = DGB.construct(centers, potential, alphas=np.tile([.7, 1.3], (32, 1)), masses=[1],
                            transformations=np.broadcast_to(rotation, (32, 2, 2)), quadrature_degree=3)
        h, s = ham.T + ham.V, ham.S
        reference = eigh(h, s, eigvals_only=True)[:4]
        for method in ['pivoted-cholesky', 'random-cholesky', 'accelerated-cholesky',
                       'lobpcg', 'fix-heiberger', 'rank-projection', 'rank-completion']:
            with self.subTest(method=method):
                options = dict(states=4)
                if method.startswith('rank-'):
                    options.update(normal_rank=len(s), seed=17)
                elif method == 'lobpcg':
                    options.update(seed=17, preconditioner='shifted', tolerance=1e-9)
                elif method == 'fix-heiberger':
                    options['epsilon'] = 1e-12
                elif method in ['random-cholesky', 'accelerated-cholesky']:
                    options['seed'] = 17
                wavefunctions = ham.get_wavefunctions(mode=method, eigensolver_options=options)
                np.testing.assert_allclose(wavefunctions.energies, reference, atol=1e-8)
                identity = lambda q, deriv_order=None: np.ones(q.shape[:-1])
                operator = lambda q, deriv_order=None: q
                np.testing.assert_allclose(wavefunctions.expectation(identity), np.eye(4), atol=2e-8)
                moments = wavefunctions.expectation(operator)
                np.testing.assert_allclose(moments, moments.swapaxes(0,1).conj(), atol=2e-8)
                phases = np.exp(1j*np.array([.4,-.7,1.,2.]))
                wavefunctions.wavefunctions = wavefunctions.wavefunctions * phases[None,:]
                np.testing.assert_allclose(wavefunctions.expectation(identity), np.eye(4), atol=2e-8)
                np.testing.assert_allclose(abs(wavefunctions.expectation(operator)), abs(moments), atol=2e-8)
                e, c = ham.diagonalize(mode=method, eigensolver_options=options)
                np.testing.assert_allclose(c.conj().T @ s @ c, np.eye(4), atol=2e-8)
                np.testing.assert_allclose(h @ c, (s @ c)*e[None, :], atol=2e-8)


if __name__ == '__main__':
    unittest.main()
