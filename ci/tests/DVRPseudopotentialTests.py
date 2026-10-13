"""Regression tests for the optional internal-coordinate DVR V' operator."""
import itertools
import math
import unittest
from unittest.mock import patch

import numpy as np
from McUtils.Scaffolding import NullLogger
from Psience.DVR import DVR, CartesianDVR, DVRException, SelfConsistentDVR
from Psience.Molecools import Molecule
from Psience.Molecools.Hamiltonian import PseudopotentialExpansion
from Psience.Modes import NormalModes
from Psience.VPT2.Terms import KineticTerms, PotentialLikeTerm


class DVRPseudopotentialTests(unittest.TestCase):
    @staticmethod
    def molecule(zmatrix=False):
        internals = ([[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]]
                     if zmatrix else {'specs': [[0, 1], [0, 2], [1, 0, 2]]})
        return Molecule(['O', 'H', 'H'],
                        [[0., 0., 0.], [1.8, 0., 0.], [-.45, 1.74, 0.]],
                        internals=internals).get_embedded_molecule(load_properties=False)

    @staticmethod
    def cartesian_coordinates(q):
        q = np.asarray(q)
        r1, r2, angle = np.moveaxis(q, -1, 0)
        carts = np.zeros(q.shape[:-1] + (3, 3))
        carts[..., 1, 0] = r1
        carts[..., 2, 0] = r2 * np.cos(angle)
        carts[..., 2, 1] = r2 * np.sin(angle)
        return carts

    @staticmethod
    def physical_vprime(q, masses):
        # Independent flattening of J = r1**2 r2**2 sin(theta):
        # V' = div(G grad(log J))/4 + grad(log J).G.grad(log J)/8.
        r1, r2, theta = np.moveaxis(q, -1, 0)
        oxygen = 1 / masses[0]
        inverse_reduced_mass = oxygen + 1 / masses[1]
        return (-inverse_reduced_mass/4 * (1/r1**2 + 1/r2**2)
                - inverse_reduced_mass/8/np.tan(theta)**2 * (1/r1**2 + 1/r2**2)
                + oxygen*np.cos(theta)/(4*r1*r2*np.tan(theta)**2))

    def test_hamiltonian_matches_vpt_terms(self):
        # Both raw internals and a nontrivial internal-mode rotation/scaling.
        for zmatrix in (False, True):
            mol = self.molecule(zmatrix)
            mol.internal_coordinates.system.name = 'WaterInternals'
            rotation = np.array([[.8, -.6, 0.], [.6, .8, 0.], [0., 0., 1.]])
            matrix = rotation * np.array([.7, 1.2, .9])[:, None]
            modes = NormalModes(mol.internal_coordinates.system, matrix,
                                inverse=np.linalg.inv(matrix), freqs=np.ones(3),
                                mass_weighted=True, frequency_scaled=True)
            ham = mol.get_hamiltonian(modes=modes)
            expected = PotentialLikeTerm(
                mol, modes=modes, undimensionalize=False,
                expansion_handling_mode='analytical', gmatrix_tolerance=None,
                check_input_gmatrix=False, logger=NullLogger()).get_terms(2)
            actual = ham.pseudopotential_expansion(2)
            for a, b in zip(actual, expected):
                np.testing.assert_allclose(a, b, atol=2e-12, rtol=2e-9)
            via_vpt = ham.get_VPT_expansions({'pseudopotential': 2})['pseudopotential']
            for a, b in zip(via_vpt, expected):
                np.testing.assert_allclose(a, b, atol=2e-12, rtol=2e-9)

    def test_raw_tensors_and_derivatives_match_vpt(self):
        rng = np.random.default_rng(37)
        def expansion(base):
            result = [base]
            for order in range(1, 5):
                tensor = rng.normal(size=(3,)*order + (3, 3)) * .02
                tensor = (tensor + tensor.swapaxes(-1, -2))/2
                tensor = sum(tensor.transpose(p+(order, order+1))
                             for p in itertools.permutations(range(order))) / math.factorial(order)
                result.append(tensor)
            return result
        g = expansion(np.eye(3)*.5)
        inertia = expansion(np.eye(3)*2)
        mol = self.molecule()
        modes = NormalModes(mol.coords.system, np.eye(9)[:, :3], freqs=np.ones(3),
                            mass_weighted=True, frequency_scaled=True)
        term = PotentialLikeTerm(mol, modes=modes, undimensionalize=False, logger=NullLogger())
        term.get_cartesians_by_modes = lambda order: [np.eye(3)] + [0]*(order-1)
        term.inertial_frame_derivatives = lambda: inertia[1:]
        with patch.object(type(mol), 'inertia_tensor', property(lambda self: inertia[0])), \
                patch.object(KineticTerms, 'get_terms', lambda self, order, **kw: g[:order+1]):
            expected = term.get_terms(2)
        actual = PseudopotentialExpansion.get_U(g, inertia)
        for a, b in zip(actual, expected):
            np.testing.assert_allclose(a, b, atol=1e-13, rtol=1e-12)

        def shifted(exp, q):
            result = []
            for derivative in range(len(exp)):
                total = 0
                for order in range(derivative, len(exp)):
                    value = exp[order]
                    for _ in range(order-derivative):
                        value = np.tensordot(q, value, axes=[0, 0])
                    total = total + value/math.factorial(order-derivative)
                result.append(total)
            return result
        def evaluate(q):
            return PseudopotentialExpansion.get_U(shifted(g, q)[:3], shifted(inertia, q)[:3])[0]
        step = 1e-4
        gradient = np.empty(3)
        hessian = np.empty((3, 3))
        for i in range(3):
            x = np.eye(3)[i]*step
            gradient[i] = (evaluate(x)-evaluate(-x))/(2*step)
            for j in range(3):
                y = np.eye(3)[j]*step
                hessian[i, j] = (evaluate(x+y)-evaluate(x-y)-evaluate(-x+y)+evaluate(-x-y))/(4*step**2)
        np.testing.assert_allclose(actual[1], gradient, atol=1e-9)
        np.testing.assert_allclose(actual[2], hessian, atol=3e-8)
        batch_g = [np.stack([t, t]) for t in g]
        batch_i = [np.stack([t, t]) for t in inertia]
        for batch, single in zip(PseudopotentialExpansion.get_U(batch_g, batch_i), actual):
            np.testing.assert_allclose(batch, np.stack([single, single]), atol=1e-13)

    def test_molecular_grid_function_matches_volume_element(self):
        q = np.array([[1.8, 1.8, 1.8], [1.9, 1.7, 1.7], [2., 1.6, 1.9]])
        for zmatrix in (False, True):
            mol = self.molecule(zmatrix)
            expected = self.physical_vprime(q, mol.atomic_masses)
            for batch_size in (1, 2, 64):
                func = mol.hamiltonian.get_pseudopotential_function(
                    self.cartesian_coordinates, batch_size=batch_size)
                np.testing.assert_allclose(func(q), expected, atol=1e-15)
                np.testing.assert_allclose(func(q.reshape((1, 3, 3))), expected[None], atol=1e-15)
            np.testing.assert_allclose(mol.hamiltonian.get_pseudopotential_function()(q), expected, atol=1e-14)

    def test_constant_tensors_have_zero_pseudopotential(self):
        g = [np.eye(2), 0, 0, 0, 0]
        inertia = [np.eye(3), 0, 0, 0, 0]
        for term in PseudopotentialExpansion.get_U(g, inertia):
            np.testing.assert_allclose(term, 0)

    def test_molecular_derivatives_match_volume_element(self):
        mol = self.molecule()
        q = np.array([1.8, np.hypot(.45, 1.74), np.arctan2(1.74, -.45)])
        actual = mol.hamiltonian.pseudopotential_expansion(2, modes=False)
        evaluate = lambda r: 8*self.physical_vprime(r, mol.atomic_masses)
        step = 1e-4
        gradient = np.empty(3)
        hessian = np.empty((3, 3))
        for i in range(3):
            x = np.eye(3)[i]*step
            gradient[i] = (evaluate(q+x)-evaluate(q-x))/(2*step)
            for j in range(3):
                y = np.eye(3)[j]*step
                hessian[i, j] = (evaluate(q+x+y)-evaluate(q+x-y)-evaluate(q-x+y)+evaluate(q-x-y))/(4*step**2)
        np.testing.assert_allclose(actual[0], evaluate(q), atol=1e-15)
        np.testing.assert_allclose(actual[1], gradient, atol=1e-11)
        np.testing.assert_allclose(actual[2], hessian, atol=1e-10)

    def test_default_and_constant_shift(self):
        calls = []
        def vprime(q):
            calls.append(len(q))
            return -.2
        dvr = DVR(domain=(-5, 5), divs=60, mass=1, logger=False,
                  potential_function=lambda q: q**2/2, pseudopotential_function=vprime,
                  num_wfns=4, diag_mode='dense')
        base = dvr.run()
        self.assertEqual(calls, [])
        self.assertIsNone(base.pseudopotential_energy)
        shifted = dvr.run(include_pseudopotential=True)
        np.testing.assert_allclose(shifted.wavefunctions.energies, base.wavefunctions.energies-.2, atol=1e-12)
        np.testing.assert_allclose(shifted.wavefunctions.frequencies(), base.wavefunctions.frequencies(), atol=1e-12)
        np.testing.assert_allclose(shifted.potential_energy, base.potential_energy)
        reused = dvr.run(include_pseudopotential=True, grid=base.grid,
                         potential_energy=base.potential_energy, kinetic_energy=base.kinetic_energy,
                         pseudopotential_energy=shifted.pseudopotential_energy)
        np.testing.assert_allclose(reused.wavefunctions.energies, shifted.wavefunctions.energies, atol=1e-12)
        scaled = dvr.run(result='pseudopotential_energy', include_pseudopotential=True, hb=2)
        np.testing.assert_allclose(scaled.pseudopotential_energy.diagonal(), -.8)

    def test_grid_dependent_term_changes_wavefunctions(self):
        dvr = CartesianDVR(domain=(-6, 6), divs=80, mass=1, logger=False,
                           potential_function=lambda q: q**2/2, num_wfns=4)
        corrected = dvr.run(include_pseudopotential=True, pseudopotential_function=lambda q: .5*q**2)
        # The correction doubles the oscillator's force constant.
        np.testing.assert_allclose(corrected.wavefunctions.energies,
                                   np.sqrt(2)*(np.arange(4)+.5), atol=2e-7)

    def test_direct_product_and_potential_optimization(self):
        for optimize in (False, True):
            options = dict(domain=[(-5, 5)]*2, divs=[24]*2, mass=[1., 1.],
                           potential_function=lambda q: np.sum(q**2, axis=-1)/2,
                           potential_optimize=optimize, po_divs=[7, 7],
                           num_wfns=4, diag_mode='dense', logger=False)
            base = DVR(**options).run()
            corrected = DVR(**options, include_pseudopotential=True,
                            pseudopotential_function=lambda q: -.12).run()
            np.testing.assert_allclose(corrected.wavefunctions.energies,
                                       base.wavefunctions.energies-.12, atol=2e-9)
            grid = corrected.grid[0] if isinstance(corrected.grid, tuple) else corrected.grid
            self.assertEqual(corrected.pseudopotential_energy.shape[0], np.prod(grid.shape[:-1]))

    def test_scf_sees_effective_potential_once(self):
        dvr = DVR(domain=[(-4, 4)]*2, divs=[15]*2, mass=[1., 1.],
                  potential_function=lambda q: np.sum(q**2, axis=-1)/2, logger=False,
                  include_pseudopotential=True, pseudopotential_function=lambda q: -.3)
        result = dvr.run(result='potential_energy')
        scf = SelfConsistentDVR(dvr)
        np.testing.assert_allclose(scf.vals.ravel(),
                                   result.potential_energy.diagonal()-.3, atol=1e-13)

    def test_invalid_pseudopotential_is_reported(self):
        dvr = CartesianDVR(domain=(-3, 3), divs=20, mass=1, logger=False,
                           potential_function=lambda q: q**2/2)
        for function in (None, lambda q: np.full(len(q), np.nan), lambda q: q[:, None]):
            with self.assertRaises(DVRException):
                dvr.run(include_pseudopotential=True, pseudopotential_function=function)


if __name__ == '__main__':
    unittest.main()
