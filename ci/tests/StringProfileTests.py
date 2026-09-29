"""Check the Psience string profile adapter with an analytic evaluator."""

import unittest

import numpy as np
from numpy.testing import assert_allclose

from McUtils.Zachary.DifferentiableFunctions import GaussianFunction
from Psience.Reactions.ProfileGenerator import GrowingString, ProfileGenerator


class Image:
    def __init__(self, coords):
        self.coords = np.asarray(coords)
        self.masses = np.ones(len(self.coords))

    def modify(self, coords):
        return Image(coords)


class Evaluator:
    distance_units = 'BohrRadius'

    def __init__(self, potential):
        self.potential = potential

    def evaluate(self, coords, order=0):
        points = coords[..., 0, :2]
        if order == 0:
            return [self.potential(points, order=0)[0]]
        if order == [1]:
            grad = np.zeros_like(coords)
            grad[..., 0, :2] = self.potential(points, order=1)[1]
            return [grad]
        raise ValueError(order)


class StringProfileTests(unittest.TestCase):
    def test_registry_and_profile_generation(self):
        self.assertIs(ProfileGenerator.get_profile_generators()['string'], GrowingString)
        potential = GaussianFunction(-10, [.5, .8], [[-2., 0.], [0., -4.]])
        evaluator = Evaluator(potential)
        base = [Image([[x, x, 0.]]) for x in np.linspace(0., 1., 6)]
        profile = GrowingString.__new__(GrowingString)
        profile.reactants = base[0]
        # A single atom has no independent translation-free gradient; use the
        # analytic Cartesian derivative to exercise the profile plumbing.
        profile._jacobian = lambda energy_evaluator: (
            lambda coords, mask: energy_evaluator.evaluate(
                coords.reshape(coords.shape[:-1] + (-1, 3)), order=[1]
            )[0].reshape(coords.shape)
        )
        before, after = profile.generate(
            base_images=base, energy_evaluator=evaluator, return_preopt=True,
            embedding_options={}, reembed=False, step_size=.002,
            max_iterations=30, tol=.1,
        )
        self.assertIs(before, base)
        self.assertEqual(len(after), len(base))
        assert_allclose(after[0].coords, base[0].coords)
        assert_allclose(after[-1].coords, base[-1].coords)
        self.assertTrue(all(isinstance(image, Image) for image in after))
        self.assertGreater(np.max([
            np.linalg.norm(new.coords - old.coords)
            for old, new in zip(base[1:-1], after[1:-1])
        ]), .01)

    def test_butane_conformers_with_rdkit_evaluator(self):
        try:
            from rdkit import Chem
            from rdkit.Chem import rdMolTransforms
        except ImportError:
            self.skipTest('RDKit is not installed')
        from McUtils.ExternalPrograms import RDMolecule
        from Psience.Molecools import Molecule

        reactant = Molecule.from_string(
            'CCCC', 'smi', optimize=True,
            confgen_opts={'random_seed': 7}, energy_evaluator='rdkit'
        )
        rd_product = Chem.Mol(reactant.rdmol.rdmol)
        conformer = rd_product.GetConformer()
        angle = rdMolTransforms.GetDihedralDeg(conformer, 0, 1, 2, 3)
        rotated_angle = angle + 120
        if rotated_angle > 180:
            rotated_angle -= 360
        rdMolTransforms.SetDihedralDeg(conformer, 0, 1, 2, 3, rotated_angle)
        product = Molecule.from_rdmol(
            RDMolecule.from_rdmol(rd_product), energy_evaluator='rdkit'
        )

        generator = GrowingString(
            reactant, product, energy_evaluator='rdkit', num_images=7
        )
        initial, relaxed = generator.generate(
            return_preopt=True, max_iterations=60, tol=.01,
            step_size=.01, max_displacement_norm=.05
        )
        self.assertEqual(len(relaxed), 7)
        self.assertTrue(all(isinstance(image, Molecule) for image in relaxed))
        assert_allclose(np.asarray(relaxed[0].coords), np.asarray(reactant.coords))
        assert_allclose(np.asarray(relaxed[-1].coords), np.asarray(product.coords))
        initial_energy = np.asarray(generator.evaluate_profile_energies(initial))
        relaxed_energy = np.asarray(generator.evaluate_profile_energies(relaxed))
        self.assertTrue(np.all(np.isfinite(relaxed_energy)))
        self.assertLess(np.max(relaxed_energy), np.max(initial_energy))
        self.assertGreater(np.max([
            np.linalg.norm(new.coords - old.coords)
            for old, new in zip(initial[1:-1], relaxed[1:-1])
        ]), .05)


if __name__ == '__main__':
    unittest.main()
