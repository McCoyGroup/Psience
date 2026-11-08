"""Public animation API coverage for amplitude adjustment bounds."""

import os as _test_os
import unittest as _test_unittest
_RUN_EXTRAS = _test_os.environ.get("MCUTILS_TEST_EXTRAS") == "1"


from html.parser import HTMLParser
import unittest

import numpy as np
from Psience.Molecools import *


class _ScaleInput(HTMLParser):
    def __init__(self, component):
        super().__init__()
        self.attributes = None
        for element in component.elements:
            self.feed(element.tostring())

    def handle_starttag(self, tag, attributes):
        attributes = dict(attributes)
        if tag == 'input' and attributes.get('data-role') == 'scale':
            self.attributes = attributes


@_test_unittest.skipUnless(_RUN_EXTRAS, "optional dependencies; set MCUTILS_TEST_EXTRAS=1")
class AnimationAdjustmentTests(unittest.TestCase):
    def setUp(self):
        self.molecule = Molecule(['H', 'H'], [[0, 0, 0], [1.4, 0, 0]], bonds=[(0, 1, 1)])
        self.modes = np.eye(6)[:, :2]

    def animate(self, which=(0, 1), **options):
        return self.molecule.animate_mode(
            which, modes=self.modes, mass_scale=False,
            animation_options={'running': False}, include_save_buttons=False, **options)

    def check_bounds(self, scene, expected, initial=1):
        component = scene.components[0]
        self.assertTrue(component.data['allowAdjustments'])
        self.assertEqual(component.data['scaleBounds'], list(expected))
        attributes = _ScaleInput(component).attributes
        self.assertEqual(float(attributes['min']), expected[0])
        self.assertEqual(float(attributes['max']), expected[1])
        self.assertEqual(float(attributes['value']), initial)

    def test_boolean_settings(self):
        self.check_bounds(self.animate(allow_adjustments=True), (0, 3))
        component = self.animate(allow_adjustments=False).components[0]
        self.assertFalse(component.data['allowAdjustments'])
        self.assertIsNone(_ScaleInput(component).attributes)

    def test_custom_bounds_single_and_multiple_modes(self):
        for which in (0, [0, 1]):
            with self.subTest(which=which):
                self.check_bounds(self.animate(which, allow_adjustments=[0, 50]), (0, 50))

    def test_coordinate_bounds_and_array_bounds(self):
        for bounds in ([0, 50], (0, 50), np.array([0., 50.])):
            with self.subTest(bounds=bounds):
                scene = self.molecule.animate_coordinate(
                    np.array([0., 0., 0., 1., 0., 0.]), allow_adjustments=bounds,
                    animation_options={'running': False}, include_save_buttons=False)
                self.check_bounds(scene, (0, 50))

    def test_initial_scale_within_bounds(self):
        self.check_bounds(self.animate(allow_adjustments=[2, 50]), (2, 50), 2)
        self.check_bounds(self.animate(allow_adjustments=[0, .5]), (0, .5), .5)

    def test_invalid_bounds(self):
        for bounds in ([], [0], [0, 1, 2], [2, 1], [1, 1], [-1, 2],
                       [0, np.inf], [np.nan, 1], 'invalid', 5):
            for which in (0, [0, 1]):
                with self.subTest(bounds=bounds, which=which):
                    with self.assertRaisesRegex(ValueError, 'allow_adjustments'):
                        self.animate(which, allow_adjustments=bounds)

    def test_playback_controls_independent_of_amplitude(self):
        scene = self.molecule.animate_mode(
            [0, 1], modes=self.modes, mass_scale=False, allow_adjustments=[0, 50],
            animation_options={'running': False, 'playback_controls': False},
            include_save_buttons=False)
        self.check_bounds(scene, (0, 50))
        self.assertNotIn('mcutils-x3d-playback', ''.join(
            element.tostring() for element in scene.components[0].elements))


if __name__ == '__main__':
    unittest.main()
