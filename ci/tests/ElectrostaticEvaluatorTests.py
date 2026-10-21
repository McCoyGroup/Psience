"""
Tests for the `Psience.Molecools` electrostatic-potential evaluators.

MACE-POLAR itself isn't required: the MACE-POLAR evaluator is exercised with a stand-in
calculator result (`density_coefficients` + smearing width), which checks the full
conversion pipeline (real e3nn-ordered moments in e Angstrom^l, GTO smearing, units).
"""

import unittest

import numpy as np
import scipy.special

from McUtils.Data import UnitsData
from Psience.Molecools import Molecule
from Psience.Molecools.Evaluator import (
    ElectrostaticPotentialEvaluator,
    ElectrostaticMultipoleExpansionEvaluator,
    ElectrostaticPotentialFunctionEvaluator,
    MACEPolarElectrostaticPotentialEvaluator,
    ChargeFunctionChargeEvaluator,
)

A2B = UnitsData.convert("Angstroms", "BohrRadius")


def direct_point_potential(sites_bohr, charges, pts_bohr):
    d = pts_bohr[:, None, :] - sites_bohr[None, :, :]
    return np.sum(charges[None, :] / np.linalg.norm(d, axis=-1), axis=-1)


def toy_charges(coords):
    """smooth geometry-dependent charges (Bohr in, e out) summing to zero"""
    coords = np.asanyarray(coords)
    base = np.array([-.8, .4, .4])
    d = coords[..., :, None, :] - coords[..., None, :, :]
    w = np.exp(-np.sum(d ** 2, axis=-1) / 8)
    q = base * np.sum(w, axis=-1)
    return q - np.mean(q, axis=-1, keepdims=True)


def toy_charge_function(coords, order=0, **opts):
    coords = np.asanyarray(coords, dtype=float)
    if order == 0:
        return toy_charges(coords)
    elif order == 1:
        # central differences: "analytic" from the evaluator framework's point of view
        h = 1e-5
        base_shape = coords.shape[:-2]
        flat = coords.reshape(base_shape + (-1,))
        n = flat.shape[-1]
        jac = []
        for k in range(n):
            e = np.zeros(n); e[k] = h
            jac.append(
                (toy_charges((flat + e).reshape(coords.shape)) - toy_charges((flat - e).reshape(coords.shape))) / (2 * h)
            )
        return np.stack(jac, axis=-2)
    raise ValueError(order)


class ElectrostaticEvaluatorTests(unittest.TestCase):

    def setUp(self):
        self.atoms = ['O', 'H', 'H']
        self.coords_ang = np.array([[0., 0., .117], [0., .757, -.469], [0., -.757, -.469]])
        rng = np.random.default_rng(3)
        pts = rng.normal(size=(12, 3))
        self.points_ang = 3.5 * pts / np.linalg.norm(pts, axis=-1)[:, None]

    def test_Registry(self):
        base, opts = ElectrostaticPotentialEvaluator.resolve_evaluator('multipoles')
        self.assertIs(base, ElectrostaticMultipoleExpansionEvaluator)
        base, opts = ElectrostaticPotentialEvaluator.resolve_evaluator('mace-polar:medium')
        self.assertIs(base, MACEPolarElectrostaticPotentialEvaluator)
        self.assertEqual(opts, {'model': 'medium'})
        base, _ = ElectrostaticPotentialEvaluator.resolve_evaluator(lambda c, order=0, **kw: 0)
        self.assertIs(base, ElectrostaticPotentialFunctionEvaluator)

    def test_FixedMultipoles(self):
        q = np.array([-.8, .4, .4])
        mu = np.array([[0., 0., .1], [0., .05, 0.], [0., -.05, 0.]])  # e Angstrom
        real = np.concatenate([q[:, None], mu[:, [1, 2, 0]]], axis=-1)
        ev = ElectrostaticMultipoleExpansionEvaluator(multipoles=real, points=self.points_ang)
        V = ev.evaluate(self.coords_ang)[0]
        sites = self.coords_ang * A2B
        pts = self.points_ang * A2B
        d = pts[:, None, :] - sites[None, :, :]
        r = np.linalg.norm(d, axis=-1)
        exact = np.sum(q / r + np.einsum('ia,nia->ni', mu * A2B, d) / r ** 3, axis=-1)
        self.assertTrue(np.allclose(V, exact))
        E = ev.evaluate(self.coords_ang, quantity='field')[0]
        self.assertEqual(E.shape, (12, 3))
        # field = -grad V, check against finite differences in the probe positions
        h = 1e-4
        fd = np.stack([
            -(ev.evaluate(self.coords_ang, points=self.points_ang + h * e)[0]
              - ev.evaluate(self.coords_ang, points=self.points_ang - h * e)[0]) / (2 * h * A2B)
            for e in np.eye(3)
        ], axis=-1)
        self.assertTrue(np.allclose(E, fd, atol=1e-7))
        self.assertEqual(ev.evaluate(self.coords_ang, quantity='field_gradient')[0].shape, (12, 3, 3))

    def test_CollapsedExpansion(self):
        q = np.array([-.8, .4, .4])
        far = self.points_ang * 10
        exact = ElectrostaticMultipoleExpansionEvaluator(multipoles=q[:, None], points=far)
        coll = ElectrostaticMultipoleExpansionEvaluator(multipoles=q[:, None], points=far, expansion_order=10)
        self.assertTrue(np.allclose(coll.evaluate(self.coords_ang)[0], exact.evaluate(self.coords_ang)[0]))

    def test_ChargeEvaluatorDerivatives(self):
        charges = ChargeFunctionChargeEvaluator(toy_charge_function, distance_units='BohrRadius')
        ev = ElectrostaticMultipoleExpansionEvaluator(charge_evaluator=charges, points=self.points_ang * A2B)
        self.assertEqual(ev.analytic_derivative_order, 1)
        self.assertEqual(ev.distance_units, 'BohrRadius')
        coords = self.coords_ang * A2B
        V, dV = ev.evaluate(coords, order=1)
        self.assertTrue(np.allclose(V, direct_point_potential(coords, toy_charges(coords), self.points_ang * A2B)))
        self.assertEqual(dV.shape, (9, 12))
        fd_ev = ElectrostaticMultipoleExpansionEvaluator(charge_evaluator=charges, points=self.points_ang * A2B)
        fd_ev.analytic_derivative_order = 0
        _, dV_fd = fd_ev.evaluate(coords, order=1)
        self.assertTrue(np.allclose(dV, dV_fd, atol=1e-6))
        # brute force
        h = 1e-4
        flat = coords.reshape(-1)
        brute = []
        for k in range(9):
            e = np.zeros(9); e[k] = h
            brute.append(
                (ev.evaluate((flat + e).reshape(3, 3))[0] - ev.evaluate((flat - e).reshape(3, 3))[0]) / (2 * h)
            )
        self.assertTrue(np.allclose(dV, np.array(brute), atol=1e-6))

    def test_MoleculeHook(self):
        mol = Molecule(self.atoms, self.coords_ang * A2B)
        charges = ChargeFunctionChargeEvaluator(toy_charge_function, distance_units='BohrRadius')
        pts = self.points_ang * A2B
        V = mol.calculate_electrostatic_potential(pts, evaluator='multipoles', charge_evaluator=charges)
        self.assertTrue(np.allclose(V, direct_point_potential(mol.coords, toy_charges(mol.coords), pts)))
        E = mol.calculate_electrostatic_potential(pts, evaluator='multipoles', charge_evaluator=charges, quantity='field')
        self.assertEqual(E.shape, (12, 3))

    def test_MACEPolarPipeline(self):
        sigma = .9  # Angstrom
        q = np.array([-.7, .35, .35])
        mu = np.array([[0., 0., .08], [0., .03, -.02], [0., -.03, -.02]])  # e Angstrom
        dens = np.concatenate([q[:, None], mu[:, [1, 2, 0]]], axis=-1)  # e3nn (y, z, x) ordering

        class _Model:
            atomic_multipoles_smearing_width = sigma

        class _Calc:
            models = [_Model()]

        class _Base:
            eval = _Calc()

        class StubMACEPolar(MACEPolarElectrostaticPotentialEvaluator):
            def evaluate_polar_properties(self, coords, properties=None):
                coords = np.asanyarray(coords)
                return {'density_coefficients': np.broadcast_to(dens, coords.shape[:-2] + dens.shape)}

        ev = StubMACEPolar(self.atoms, energy_evaluator=_Base(), points=self.points_ang)
        self.assertEqual(ev.model_smearing, sigma)
        V = ev.evaluate(self.coords_ang)[0]
        # independent reference: smeared charges (erf) + smeared dipoles (-p . grad of a smeared unit charge)
        s = sigma * A2B
        sites = self.coords_ang * A2B
        pts = self.points_ang * A2B
        d = pts[:, None, :] - sites[None, :, :]
        r = np.linalg.norm(d, axis=-1)
        u = r / (np.sqrt(2) * s)
        V0 = scipy.special.erf(u) / r
        dV0_dr = (2 / np.sqrt(np.pi) * np.exp(-u ** 2) / (np.sqrt(2) * s) * r - scipy.special.erf(u)) / r ** 2
        grad_V0 = dV0_dr[..., None] * d / r[..., None]
        # a smeared dipole is -p . grad of a smeared unit charge: V_dip = -p . grad_r V0
        exact = np.sum(q * V0 - np.einsum('ia,nia->ni', mu * A2B, grad_V0), axis=-1)
        self.assertTrue(np.allclose(V, exact))
        # finite at the nuclei
        V_nuc = ev.evaluate(self.coords_ang, points=self.coords_ang)[0]
        self.assertTrue(np.all(np.isfinite(V_nuc)))
        # nuclear-coordinate derivatives fall back to finite differences
        V1, dV1 = ev.evaluate(self.coords_ang, order=1)
        self.assertEqual(dV1.shape, (9, 12))


if __name__ == '__main__':
    unittest.main()
