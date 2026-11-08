"""Isotope identity and masses survive Psience/RDKit conversion."""

import os as _test_os
import unittest as _test_unittest
_RUN_EXTRAS = _test_os.environ.get("MCUTILS_TEST_EXTRAS") == "1"

import unittest

import numpy as np

from McUtils.Data import AtomData, UnitsData
from McUtils.ExternalPrograms import RDMolecule
from Psience.Molecools import Molecule


class MoleculeIsotopeRoundTripTests(unittest.TestCase):
    coords = np.array([[0., 0., 0.], [.96, 0., 0.], [-.24, .93, 0.]])
    bonds = [[0, 1, 1], [0, 2, 1]]

    def wrapped(self):
        return RDMolecule.from_coords(['O18', 'D', 'T'], self.coords, bonds=self.bonds)

    @_test_unittest.skipUnless(_RUN_EXTRAS, "optional dependencies; set MCUTILS_TEST_EXTRAS=1")
    def test_psience_to_rdkit_and_back(self):
        bohr = self.coords * UnitsData.convert('Angstroms', 'BohrRadius')
        mol = Molecule(['O18', 'D', 'T'], bohr, bonds=self.bonds)
        self.assertEqual(mol.rdmol.atoms, ['O', 'H', 'H'])
        self.assertEqual(mol.rdmol.isotopes, [18, 2, 3])
        restored = Molecule.from_rdmol(mol.rdmol)
        self.assertEqual(restored.atoms, ('O18', 'D', 'T'))
        np.testing.assert_array_equal(restored.masses, mol.masses)
        np.testing.assert_allclose(np.asarray(restored.coords), np.asarray(mol.coords), atol=1e-14)

    @_test_unittest.skipUnless(_RUN_EXTRAS, "optional dependencies; set MCUTILS_TEST_EXTRAS=1")
    def test_wrapped_raw_molecule_and_conformer_imports(self):
        wrapped = self.wrapped()
        for source in [wrapped, wrapped.rdmol, wrapped.rdmol.GetConformer()]:
            with self.subTest(source=type(source).__name__):
                mol = Molecule.from_rdmol(source)
                self.assertEqual(mol.atoms, ('O18', 'D', 'T'))
                self.assertEqual(mol.rdmol.isotopes, [18, 2, 3])
                np.testing.assert_array_equal(mol.masses, [AtomData[a, 'Mass'] for a in mol.atoms])

    @_test_unittest.skipUnless(_RUN_EXTRAS, "optional dependencies; set MCUTILS_TEST_EXTRAS=1")
    def test_rebuilding_from_psience_preserves_heavy_isotopes(self):
        mol = Molecule.from_rdmol(self.wrapped())
        rebuilt = RDMolecule.from_mol(mol, coord_unit='BohrRadius')
        self.assertEqual(rebuilt.isotopes, [18, 2, 3])
        np.testing.assert_allclose(rebuilt.coords, self.coords, atol=1e-14)
        np.testing.assert_array_equal(rebuilt.masses, mol.masses)

    @_test_unittest.skipUnless(_RUN_EXTRAS, "optional dependencies; set MCUTILS_TEST_EXTRAS=1")
    def test_explicit_mass_override_is_respected(self):
        mol = Molecule.from_rdmol(self.wrapped(), masses=[18., 2., 3.])
        np.testing.assert_array_equal(mol.masses, [18., 2., 3.])
        self.assertEqual(mol.atoms, ('O18', 'D', 'T'))

    @_test_unittest.skipUnless(_RUN_EXTRAS, "optional dependencies; set MCUTILS_TEST_EXTRAS=1")
    def test_bond_guessing_accepts_d_and_t(self):
        bohr = self.coords * UnitsData.convert('Angstroms', 'BohrRadius')
        mol = Molecule(['O', 'D', 'T'], bohr)
        self.assertEqual({tuple(sorted(b[:2])) for b in mol.bonds}, {(0, 1), (0, 2)})
        self.assertEqual(mol.rdmol.isotopes, [0, 2, 3])

    @_test_unittest.skipUnless(_RUN_EXTRAS, "optional dependencies; set MCUTILS_TEST_EXTRAS=1")
    def test_carbon13_and_explicit_primary_isotopes(self):
        wrapped = RDMolecule.from_coords(['C13', 'C12', 'H1'], self.coords, bonds=[])
        mol = Molecule.from_rdmol(wrapped)
        self.assertEqual(mol.atoms, ('C13', 'C', 'H'))
        self.assertEqual(mol.rdmol.isotopes, [13, 12, 1])
        np.testing.assert_array_equal(mol.masses, [AtomData[a, 'Mass'] for a in ('C13', 'C12', 'H1')])


if __name__ == '__main__':
    unittest.main()
