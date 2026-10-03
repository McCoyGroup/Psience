"""Regressions for small `Molecule` bugs found while building readouts (2026-10)."""
import os
import unittest

from Psience.Molecools import Molecule

TEST_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "TestData")


class PointGroupCacheTests(unittest.TestCase):
    def test_point_group_is_computed_lazily_and_reset_with_coords(self):
        # `_pg` was never initialized, so `mol.point_group` raised AttributeError
        mol = Molecule.from_file(os.path.join(TEST_DATA, "HOH_freq.fchk"))
        self.assertEqual(str(mol.point_group), "PointGroup<C2v>")
        self.assertIsNotNone(mol._pg)
        mol.coords = mol.coords
        self.assertIsNone(mol._pg)
        mol.point_group = "custom"
        self.assertEqual(mol.point_group, "custom")


if __name__ == "__main__":
    unittest.main()
