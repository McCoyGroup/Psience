import unittest
import numpy as np
from Psience.Nonlinear import (nonlinear_response_generators, orientational_averaging_function,
                               four_wave_averaging_function)
from Psience.Nonlinear.NonlinearResponse import prep_nonlinear_transition_data, expand_transition_data


class OrientationalResponseTests(unittest.TestCase):
    @staticmethod
    def transition_data(moment=None):
        return {((0,), (1,)): {'frequency': 1., 'coupling': 0.,
                              'transition_moment': [0., 0., 1.] if moment is None else moment}}

    def test_six_and_eight_factor_end_to_end_path_response(self):
        for order in (4, 6, 8):
            opts = dict(paths=[[0]*(order-1)], num_interactions=order, driving_frequency=False)
            baseline = nonlinear_response_generators(self.transition_data(), **opts)
            delays = [np.array([0., .0001])]*(order-1)
            base = baseline.evaluate_times(*delays)
            self.assertGreater(np.max(np.abs(base)), .1)
            for method in ('wigner', 'bonvicini', 'young') if order <= 6 else ('wigner',):
                response = nonlinear_response_generators(self.transition_data(), polarization='Z'*order,
                            orientational_averaging_method=method, **opts)
                np.testing.assert_allclose(response.evaluate_times(*delays), base/(order+1), atol=1e-14)
                with self.assertRaises(ValueError):
                    response.evaluate_times(*delays[:-1])

    def test_existing_four_factor_formula(self):
        rng = np.random.default_rng(574)
        pairings = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
        for _ in range(20):
            lab, mol = rng.normal(size=(4, 3)), rng.normal(size=(4, 3))
            a = np.array([np.prod([lab[i]@lab[j] for i, j in p]) for p in pairings])
            b = np.array([np.prod([mol[i]@mol[j] for i, j in p]) for p in pairings])
            expected = a@(5*np.eye(3)-np.ones((3, 3)))@b/30
            np.testing.assert_allclose(four_wave_averaging_function(lab)(mol), expected, atol=3e-14)

    def test_complex_transition_dtype_normalization_and_field_conjugations(self):
        moment = np.array([1., 1j, 0.])
        data = self.transition_data(moment)
        td = prep_nonlinear_transition_data(data)
        np.testing.assert_array_equal(td.transition_moments[0, 1], moment)
        np.testing.assert_array_equal(td.transition_moments[1, 0], moment.conj())
        expanded = expand_transition_data(td, td.states)
        np.testing.assert_array_equal(expanded.transition_moments, td.transition_moments)
        opts = dict(paths=[[0, 0, 0]], num_interactions=4, driving_frequency=False)
        response = nonlinear_response_generators(data, polarization=np.tile(moment, (4, 1)),
                   polarization_conjugations=[False, True, False, True], **opts)
        self.assertTrue(np.all(np.isfinite(response.evaluate_times(0., 0., 0.))))
        captured = {}
        def capture(space, blocks, frequencies, magnitudes, directions, *args, **kwargs):
            captured.update(space=space, magnitudes=magnitudes, directions=directions)
            return lambda *times: 0
        nonlinear_response_generators(data, response_function_generator=capture, **opts)
        a, b = captured['space'].find(td.states)
        expected = np.zeros_like(captured['magnitudes'])
        expected[a, b] = expected[b, a] = np.sqrt(2)
        np.testing.assert_allclose(captured['magnitudes'], expected)
        np.testing.assert_array_equal(captured['directions'][a, a], np.zeros(3))
        np.testing.assert_allclose(captured['directions'][a, b], moment/np.sqrt(2))
        for method in ('wigner', 'bonvicini', 'young'):
            average = orientational_averaging_function(np.tile(moment, (2, 1)), method, [False, True])
            self.assertAlmostEqual(average(np.tile(moment, (2, 1)), conjugations=[False, True]), 4/3)

    def test_arity_and_labels(self):
        with self.assertRaises(ValueError):
            four_wave_averaging_function('ZZZZZZ')
        with self.assertRaises(ValueError):
            orientational_averaging_function('XXQZ')
        with self.assertRaises(ValueError):
            nonlinear_response_generators(self.transition_data(), paths=[[0]*5], num_interactions=6,
                                           polarization='ZZZZ', driving_frequency=False)


if __name__ == '__main__':
    unittest.main()
