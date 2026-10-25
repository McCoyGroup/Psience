import unittest
import numpy as np
from McUtils.Coordinerds import CartesianCoordinates3D
from Psience.Modes import NormalModes
from Psience.Vibronic.FCFs import FranckCondonModel


class FCFPairingRegressionTests(unittest.TestCase):
    @staticmethod
    def model():
        rng = np.random.default_rng(23462)
        a = np.eye(3); b = np.linalg.qr(rng.normal(size=(3, 3)))[0]
        gs = NormalModes(CartesianCoordinates3D, a, inverse=a.T, freqs=[.8, 1.3, 1.7], origin=[0., 0., 0.])
        es = NormalModes(CartesianCoordinates3D, b, inverse=b.T, freqs=[1., 1.5, 1.9], origin=[.2, -.3, .1])
        return FranckCondonModel(gs, es, embed=False, mass_weight=False, logger=False)

    def test_existing_partition_pair_counts(self):
        for degree, count in ((2, 1), (4, 3), (6, 15), (8, 105)):
            self.assertEqual(len(FranckCondonModel.get_pairing_evaluation_plan((degree,))), count)

    def test_high_level_methods_rotations_and_hot_bands(self):
        model = self.model()
        excited = [[0, 0, 0], [6, 0, 0], [3, 2, 1], [2, 2, 2], [4, 2, 2]]
        ground = [[0, 0, 0], [1, 1, 0]]
        for order in ('gs', 'es'):
            for center in ('gs', 'es'):
                options = dict(ground_states=ground, return_states=False, rotation_order=order, rotation_center=center)
                expected = model.get_overlaps(excited, integration_method='partitions', **options)
                observed = model.get_overlaps(excited, integration_method='pairings', **options)
                np.testing.assert_allclose(observed, expected, atol=2e-12, rtol=2e-11)
                default = model.get_overlaps(excited, **options)
                np.testing.assert_array_equal(default, expected)

    def test_moment_cache_belongs_to_each_geometry(self):
        model = self.model()
        states = [[6, 0, 0], [4, 2, 0]]
        for displacement in ([.2, -.3, .1], [.5, .1, -.2], [.2, -.3, .1]):
            model.es_nms.origin = np.array(displacement)
            expected = model.get_overlaps(states, return_states=False, integration_method='partitions')
            actual = model.get_overlaps(states, return_states=False, integration_method='pairings')
            np.testing.assert_allclose(actual, expected, atol=2e-12, rtol=2e-11)

    def test_spectrum_interface_and_cutoff(self):
        model = self.model()
        states = [[0, 0, 0], [6, 0, 0], [2, 2, 2]]
        for cutoff in (1e-20, 1e-5, .01):
            expected = model.get_spectrum(states, integration_method='partitions', duschinsky_cutoff=cutoff)
            actual = model.get_spectrum(states, integration_method='pairings', duschinsky_cutoff=cutoff)
            np.testing.assert_allclose(actual.frequencies, expected.frequencies, atol=1e-12)
            np.testing.assert_allclose(actual.intensities, expected.intensities, atol=2e-12, rtol=2e-11)

    def test_invalid_method(self):
        with self.assertRaises(ValueError):
            self.model().get_overlaps([[0, 0, 0]], integration_method='invalid')

    def test_bounded_moment_cache(self):
        model = self.model()
        states = [[6, 0, 0], [3, 2, 1], [2, 2, 2], [6, 0, 0]]
        expected = model.get_overlaps(states, return_states=False, integration_method='partitions')
        old_size = FranckCondonModel.pairing_cache_size
        try:
            FranckCondonModel.pairing_cache_size = 2
            actual = model.get_overlaps(states, return_states=False, integration_method='pairings')
        finally:
            FranckCondonModel.pairing_cache_size = old_size
        np.testing.assert_allclose(actual, expected, atol=2e-12, rtol=2e-11)
