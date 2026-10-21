"""
Tests for `Psience.ElectronicStructure`.

`TestData/mace_polar_states_20261020.json` holds real MACE-POLAR-1-S outputs (Entropy job 157007:
charge and spin-channel monopoles; the model's l = 1 coefficients were exactly zero) for an
eMolecules conformer in three charge states and for several radicals. The MACE-POLAR state model is
exercised by replaying those outputs through `MACEPolarBackend.evaluate_polar_properties`, so
nothing here needs MACE (or anything from `Psience.Molecools` beyond `Molecule`).
"""

import json
import os
import unittest

import numpy as np

from McUtils.Data import AtomData, UnitsData
from Psience.Molecools import Molecule
from Psience.Molecools.Evaluator import ChargeFunctionChargeEvaluator
from Psience.ElectronicStructure import (
    AtomCenteredDensity, MACEPolarStateModel, ChargeEvaluatorStateModel, BackendChargeStateModel,
    ElectronicStateModel, ElectronicStructureBackend, MACEPolarBackend, SpinDensity, FukuiFunctions
)

A2B = UnitsData.convert("Angstroms", "BohrRadius")
DATA = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "TestData", "mace_polar_states_20261020.json")))
SIGMA = DATA["smearing_width_angstrom"]


class _Model:
    atomic_multipoles_smearing_width = SIGMA
class _Calc:
    models = [_Model()]


def replay_model(name):
    system = DATA["systems"][name]
    x = np.array(system["positions_angstrom"])
    states = {(s["charge"], s["multiplicity"]): s for s in system["states"]}

    class Replay(MACEPolarBackend):
        def evaluate_polar_properties(self, coords, properties=None, charge=None, multiplicity=None):
            assert np.allclose(np.asanyarray(coords), x, atol=1e-8)
            s = states[(charge, multiplicity)]
            up = np.concatenate([np.array(s["q_up"])[:, None], np.zeros((len(x), 3))], axis=-1)
            down = np.concatenate([np.array(s["q_down"])[:, None], np.zeros((len(x), 3))], axis=-1)
            return {"energy": s["energy_eV"], "density_coefficients": up + down,
                    "spin_charge_density": np.stack([up, down], axis=-2)}

    backend = Replay(system["symbols"], calculator=_Calc())
    return MACEPolarStateModel(backend), system["symbols"], x, states


def grid(lim=9., n=55):
    g = np.linspace(-lim, lim, n)
    pts = np.stack(np.meshgrid(g, g, g, indexing='ij'), axis=-1).reshape(-1, 3)
    return pts, (g[1] - g[0]) ** 3


class ElectronicStructureTests(unittest.TestCase):

    def test_AtomCenteredDensity(self):
        rng = np.random.default_rng(1)
        centers = rng.uniform(-1, 1, (5, 3))
        uniform = AtomCenteredDensity.from_populations(np.full(5, .2), centers, smearing=1.)
        self.assertTrue(np.isclose(uniform.participation_ratio(), 5))
        single = AtomCenteredDensity.from_populations(np.eye(5)[2], centers, smearing=1.)
        self.assertTrue(np.isclose(single.participation_ratio(), 1))
        real = rng.normal(size=(5, 4)) * .3
        d = AtomCenteredDensity(real, centers, smearing=1.2)
        pts, dV = grid(8., 49)
        rho = d.density(pts)
        self.assertTrue(np.isclose(np.sum(rho) * dV, d.total, atol=1e-6))
        # the density's dipole about the origin = site dipoles + populations * centers
        self.assertTrue(np.allclose((rho * dV) @ pts, d.dipoles.sum(0) + d.populations @ centers, atol=1e-6))
        self.assertTrue(np.allclose((d - d).multipoles, 0))
        with self.assertRaises(ValueError):
            d - AtomCenteredDensity(real, centers + 1, smearing=1.2)
        with self.assertRaises(ValueError):
            AtomCenteredDensity(real, centers).density(pts)

    def test_MACEPolarSpinDensity(self):
        for name in ("benzyl_radical", "allyl_radical", "phenoxyl_radical"):
            model, atoms, x, _ = replay_model(name)
            sd = SpinDensity.compute(model, x, charge=0, multiplicity=2)
            self.assertTrue(np.isclose(sd.density.total, 1, atol=1e-6), name)
            self.assertEqual(sd.density.centers.shape, (len(atoms), 3))
            self.assertTrue(np.allclose(sd.density.centers, x * A2B))
        model, atoms, x, states = replay_model("water_dimer")
        trip = SpinDensity.compute(model, x, charge=0, multiplicity=3)
        self.assertTrue(np.isclose(trip.density.total, 2, atol=1e-6))
        pts, dV = grid(20., 61)  # sigma = 1.5 A = 2.8 Bohr, so the grid must reach well past the atoms
        self.assertTrue(np.isclose(np.sum(trip.density.density(pts)) * dV, 2, atol=1e-4))
        diff = SpinDensity.compute(model, x, charge=0, multiplicity=3, reference=(0, 1))
        s3, s1 = states[(0, 3)], states[(0, 1)]
        expected = (np.array(s3["q_up"]) - s3["q_down"]) - (np.array(s1["q_up"]) - s1["q_down"])
        self.assertTrue(np.allclose(diff.populations, expected))
        self.assertTrue(np.isclose(diff.density.total, 2, atol=1e-6))

    def test_MACEPolarFukui(self):
        model, atoms, x, states = replay_model("emols7")
        f = FukuiFunctions.compute(model, x, charge=0, multiplicity=1)
        self.assertEqual(f.cation.multiplicity, 2)
        self.assertEqual(f.anion.multiplicity, 2)
        q = {k: np.array(s["q_up"]) + np.array(s["q_down"]) for k, s in states.items()}
        self.assertTrue(np.allclose(f.f_minus.populations, q[(1, 2)] - q[(0, 1)]))
        self.assertTrue(np.allclose(f.f_plus.populations, q[(0, 1)] - q[(-1, 2)]))
        for g in (f.f_plus, f.f_minus, f.f_zero):
            self.assertTrue(np.isclose(g.total, 1, atol=1e-6))
        self.assertTrue(np.isclose(f.dual_descriptor.total, 0, atol=1e-6))
        Z = np.array([AtomData[a, "Number"] for a in atoms])
        self.assertTrue(np.allclose(f.reference.electron_density.populations, Z - q[(0, 1)]))
        s = f.summary()
        self.assertGreater(s["f_minus"]["participation_ratio"], 1)
        self.assertTrue(np.isclose(f.reference.energy, states[(0, 1)]["energy_eV"] * UnitsData.convert("ElectronVolts", "Hartrees")))

    def test_FromMolecule(self):
        model, atoms, x, states = replay_model("emols7")
        mol = Molecule(atoms, x * A2B)
        self.assertIs(ElectronicStateModel.from_molecule(mol, model=model), model)
        f = FukuiFunctions.from_molecule(mol, model=model, which=('-',))
        self.assertIsNone(f.anion)
        self.assertTrue(np.isclose(f.f_minus.total, 1, atol=1e-6))
        sd = SpinDensity.from_molecule(mol, model=model, charge=1, multiplicity=2, reference=(0, 1))
        self.assertTrue(np.isclose(sd.density.total, 1, atol=1e-6))

    def test_ChargeEvaluatorModel(self):
        n = 4
        x = np.random.default_rng(3).uniform(-2, 2, (n, 3))
        atoms = ['C', 'O', 'N', 'C']

        def factory(charge, multiplicity):
            base = np.array([.1, -.3, .25, -.05])
            fn = lambda coords, order=0, **kw: base + charge * np.array([.4, .1, .3, .2])
            return ChargeFunctionChargeEvaluator(fn, distance_units='Angstroms')

        model = ChargeEvaluatorStateModel(factory, atoms)
        f = FukuiFunctions.compute(model, x)
        self.assertTrue(np.allclose(f.f_minus.populations, [.4, .1, .3, .2]))
        self.assertTrue(np.allclose(f.f_plus.populations, [.4, .1, .3, .2]))
        with self.assertRaises(ValueError):
            SpinDensity.compute(model, x)
        smeared = ChargeEvaluatorStateModel(factory, atoms, smearing=1.)
        self.assertEqual(FukuiFunctions.compute(smeared, x).f_minus.density(x * A2B).shape, (n,))

    def test_BackendChargeModel(self):
        # any backend providing charges, evaluated per state through its charge/multiplicity overrides
        class LinearCharges(ElectronicStructureBackend):
            property_units = {'charges': 'ElementaryCharge', 'energy': 'Hartrees'}
            analytic_derivative_orders = {'charges': 0, 'energy': 0}
            calls = []
            def evaluate_properties(self, coords, orders, charge, multiplicity, **opts):
                self.calls.append((charge, multiplicity, tuple(orders)))
                base = np.array([.1, -.3, .25, -.05])
                q = np.broadcast_to(base + charge * np.array([.4, .1, .3, .2]), coords.shape[:-1])
                return {'charges': [q], 'energy': [np.full(len(coords), -1. * charge)]}

        atoms = ['C', 'O', 'N', 'C']
        x = np.random.default_rng(3).uniform(-2, 2, (4, 3))
        model = BackendChargeStateModel(LinearCharges(atoms))
        f = FukuiFunctions.compute(model, x)
        self.assertTrue(np.allclose(f.f_minus.populations, [.4, .1, .3, .2]))
        self.assertTrue(np.allclose(f.f_plus.populations, [.4, .1, .3, .2]))
        self.assertTrue(np.isclose(f.cation.energy, -1))
        self.assertEqual(sorted(LinearCharges.calls), [(-1, 2, ('charges', 'energy')), (0, 1, ('charges', 'energy')), (1, 2, ('charges', 'energy'))])
        self.assertTrue(np.allclose(f.reference.charge_density.centers, x * A2B))


if __name__ == '__main__':
    unittest.main()
