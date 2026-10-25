import itertools
import unittest
import numpy as np
from Psience.BasisReps import BasisStateSpace


class StateQuantaBoundsTests(unittest.TestCase):
    def test_inclusive_quanta_and_energy_bounds(self):
        frequencies = np.array([.29, .1, .17])
        for minimum, maximum in ((None, 0), (None, 6), (6, 6), (3, 4)):
            actual = BasisStateSpace.states_under_freq_threshold(frequencies, 1., min_quanta=minimum, max_quanta=maximum)
            expected = [state for state in itertools.product(range(maximum + 1), repeat=3)
                        if (minimum is None or sum(state) >= minimum) and sum(state) <= maximum
                        and np.dot(state, frequencies) <= 1. + 1e-14]
            self.assertEqual(set(map(tuple, actual)), set(expected))

    def test_initial_state_respects_both_bounds(self):
        cases = [dict(initial_state=[1, 1], min_quanta=3, max_quanta=3),
                 dict(initial_state=[2, 1], min_quanta=2, max_quanta=2),
                 dict(initial_state=[2, 1], min_quanta=3, max_quanta=3)]
        for options in cases:
            actual = BasisStateSpace.states_under_freq_threshold(np.array([.1, .2]), 1., **options)
            for state in actual:
                self.assertGreaterEqual(sum(state), options['min_quanta'])
                self.assertLessEqual(sum(state), options['max_quanta'])

    def test_initial_state_above_energy_window_is_excluded(self):
        actual = BasisStateSpace.states_in_windows(np.array([1., 2.]), [[0., 1.]], initial_state=[0, 1], max_quantas=[2])[0]
        self.assertEqual(len(actual), 0)
