"""Display radii stay independent of optional scientific radius fields."""

import unittest

from McUtils.Data import AtomData
from Psience.Molecools import Molecule
from Psience.Molecools.Visualizations import MoleculePlotter


class AtomicRadiusPlotTests(unittest.TestCase):
    def plotter(self, molecule=None):
        return MoleculePlotter(molecule, [], mode='x3d', backend='x3d', atoms=['H'])

    def test_small_icon_never_requires_or_uses_vdw(self):
        standalone = self.plotter()
        attached = self.plotter(Molecule)
        for icon in (.3, .5, .79, 1.4):
            for record in ({"IconRadius": icon}, {"IconRadius": icon, "VanDerWaalsRadius": 9}):
                with self.subTest(icon=icon, record=record):
                    self.assertEqual(Molecule._get_atomic_radius(record), icon)
                    self.assertEqual(standalone._get_atomic_radius(record), icon)
                    self.assertEqual(attached._get_atomic_radius(record), icon)

    def test_explicit_radius_fields_are_honored_and_missing_values_raise(self):
        for resolve in (Molecule._get_atomic_radius, self.plotter()._get_atomic_radius,
                        self.plotter(Molecule)._get_atomic_radius):
            self.assertEqual(resolve({"IconRadius": .5, "CovalentRadius": 1.6}, "CovalentRadius"), 1.6)
            with self.assertRaises(KeyError):
                resolve({"IconRadius": .5}, "VanDerWaalsRadius")
            with self.assertRaises(KeyError):
                resolve({"VanDerWaalsRadius": 2})

    def test_all_elements_resolve_to_icon_without_scientific_fields(self):
        plotter = self.plotter()
        for number in range(1, 119):
            data = AtomData[number]
            record = {"IconRadius": data["IconRadius"]}
            with self.subTest(number=number):
                self.assertEqual(Molecule._get_atomic_radius(record), data["IconRadius"])
                self.assertEqual(plotter._get_atomic_radius(record), data["IconRadius"])


if __name__ == "__main__":
    unittest.main()
