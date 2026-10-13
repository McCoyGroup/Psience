"""`Molecule.to_string('inchi', include_tag=True)` -> `Molecule.from_string(..., 'inchi')` round trips, with and without hydrogens in the tag"""
import unittest
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from McUtils.ExternalPrograms import RDMolecule
from Psience.Molecools import Molecule


def distances(x):
    return np.linalg.norm(x[:, None] - x[None], axis=-1)


class MoleculeInChITagRoundTripTests(unittest.TestCase):
    SMILES = ["CCO", "CC(=O)O", "CC(=O)N", "CC.O", "CO"]

    @staticmethod
    def molecule(smiles):
        rd = Chem.AddHs(Chem.MolFromSmiles(smiles))
        AllChem.EmbedMolecule(rd, randomSeed=7)
        return Molecule.from_rdmol(RDMolecule(rd.GetConformer()))

    def roundtrip(self, smiles, keep_hydrogens):
        mol = self.molecule(smiles)
        tag_opts = dict(remove_implicit_hydrogens=False) if keep_hydrogens else {}
        tagged, order = mol.to_string('inchi', include_tag=True, return_reordering=True, **tag_opts)
        # whether the tag holds hydrogens is detected
        loaded = Molecule.from_string(tagged, 'inchi', add_implicit_hydrogens=False)
        return mol, loaded, mol.rdmol.coords[order,]

    def test_heavy_atom_tag(self):
        for smiles in self.SMILES:
            with self.subTest(smiles):
                mol, loaded, X = self.roundtrip(smiles, False)
                Y = loaded.rdmol.coords[:len(X)]
                self.assertTrue(np.allclose(distances(X), distances(Y), atol=.05))

    def test_tag_with_hydrogens(self):
        for smiles in self.SMILES:
            with self.subTest(smiles):
                mol, loaded, X = self.roundtrip(smiles, True)
                self.assertEqual(len(loaded.rdmol.coords), len(mol.rdmol.coords))
                self.assertEqual(sorted(loaded.atoms), sorted(mol.atoms))
                self.assertTrue(np.allclose(distances(X), distances(loaded.rdmol.coords), atol=.05))

    def test_inchi_itself_is_unchanged_by_the_tag(self):
        mol = self.molecule("CCO")
        plain = mol.to_string('inchi')
        for opts in [dict(), dict(remove_implicit_hydrogens=False)]:
            self.assertEqual(mol.to_string('inchi', include_tag=True, **opts).partition('_')[0], plain)


if __name__ == "__main__":
    unittest.main()
