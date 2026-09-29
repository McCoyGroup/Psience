"""Regression checks for interpolation with geometry-local auto internals."""

import unittest

import numpy as np
from numpy.testing import assert_allclose

from McUtils.Coordinerds import (
    CoordinateSet, CoordinateSystem, CoordinateSystemConverters,
    SimpleCoordinateSystemConverter
)
from Psience.Molecools import Molecule
from Psience.Reactions.ProfileGenerator import InterpolatingProfileGenerator


class AutoInternalsInterpolationTests(unittest.TestCase):
    def test_temporary_registered_converter_is_released(self):
        source = CoordinateSystem()
        target = CoordinateSystem()
        converter = SimpleCoordinateSystemConverter(
            (source, target), lambda coords, **kw: coords + 1
        )
        source.registered_converters = (converter,)
        before = len(CoordinateSystemConverters.converters)
        result = CoordinateSet([1.0], source).convert(target)
        assert_allclose(result, [2.0])
        self.assertEqual(len(CoordinateSystemConverters.converters), before)
        self.assertFalse(source._preregistered)

        def fail(coords, **kw):
            raise RuntimeError('conversion failed')

        source.registered_converters = (
            SimpleCoordinateSystemConverter((source, target), fail),
        )
        with self.assertRaisesRegex(RuntimeError, 'conversion failed'):
            CoordinateSet([1.0], source).convert(target)
        self.assertEqual(len(CoordinateSystemConverters.converters), before)
        self.assertFalse(source._preregistered)

    def test_butane_local_internals_do_not_grow_converter_registry(self):
        try:
            import rdkit  # noqa: F401
        except ImportError:
            self.skipTest('RDKit is not installed')

        reactant, product = Molecule.from_string(
            'CCCC', 'smi', num_confs=2, take_min=False, optimize=True,
            confgen_opts={'random_seed': 7}, energy_evaluator='rdkit'
        )
        endpoints = [np.asarray(mol.coords).copy() for mol in (reactant, product)]
        for choice in ('auto', 'zmatrix'):
            with self.subTest(internals=choice):
                generator = InterpolatingProfileGenerator(
                    reactant, product, internals=choice, num_images=5
                )
                # Each frame retains its own reference geometry and internal system.
                internals = generator.interpolator.interpolator.internals
                self.assertIsNot(internals[0].system, internals[1].system)
                carts = generator.interpolator.interpolator.coords
                for cart, internal in zip(carts, internals):
                    self.assertIs(internal.system.coords.system, cart.system)
                    self.assertTrue(any(
                        conv.types[0] is cart.system and conv.types[1] is internal.system
                        for conv in internal.system.registered_converters
                    ))

                registered_before = len(CoordinateSystemConverters.converters)
                images = generator.generate()
                registered_after = len(CoordinateSystemConverters.converters)
                self.assertEqual(registered_after, registered_before)
                self.assertEqual(len(images), 5)
                assert_allclose(np.asarray(images[0].coords), endpoints[0], atol=1e-8)
                assert_allclose(np.asarray(images[-1].coords), endpoints[1], atol=1e-8)
                self.assertTrue(np.all(np.isfinite([
                    np.asarray(image.coords) for image in images
                ])))


if __name__ == '__main__':
    unittest.main()
