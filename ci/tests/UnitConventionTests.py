"""Psience caller contracts and numerical checks for explicit energy equivalences.

The checked conversion boundary verifies explicit caller conventions and delegates
to the real McUtils conversion pipeline, including its equivalence policy.
"""
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

from McUtils.Data import UnitsData
from Psience.AIMD.Simulator import AIMDSimulator, LangevinThermostat
from Psience.BasisReps import BasisStateSpace, HarmonicOscillatorProductBasis
from Psience.Nonlinear import propagate_density_matrix_paths
from Psience.Nonlinear.NonlinearResponse import NonlinearResponseFunction


class UnitConventionTests(unittest.TestCase):
    def setUp(self):
        self.original_convert = UnitsData.convert
        self.calls = []
        self.hartree_joules = self.original_convert('Hartrees', 'Joules')
        self.planck = UnitsData.constant('PlanckConstant')
        self.boltzmann = UnitsData.constant('BoltzmannConstant')
        self.speed_of_light = UnitsData.constant('SpeedOfLightInVacuum')
        conversion_patch = patch.object(UnitsData, 'convert', side_effect=self.checked_conversion)
        conversion_patch.start()
        self.addCleanup(conversion_patch.stop)

    def checked_conversion(self, source, target, *, equivalences=None):
        families = {
            'Hartrees': 'energy', 'Joules': 'energy', 'ElectronVolts': 'energy',
            'Kelvins': 'thermal', 'Wavenumbers': 'wavenumber', 'Hertz': 'frequency',
        }
        source_family, target_family = families.get(source), families.get(target)
        required = set()
        if source_family is not None and target_family is not None and source_family != target_family:
            required = {family for family in (source_family, target_family) if family != 'energy'}
        enabled = {'wavenumber'} if equivalences is None else set(equivalences)
        self.assertTrue(required <= enabled,
                        f'{source} -> {target} needs explicit {sorted(required - enabled)}')
        self.calls.append((source, target, equivalences))
        return self.original_convert(source, target, equivalences=equivalences)

    def test_temperature_uses_equipartition_in_kelvin(self):
        simulator = object.__new__(AIMDSimulator)
        simulator._mass = np.array([[2.], [3.]])
        simulator._n_dof = 3
        velocities = np.array([[.1, -.2, .3], [.2, .1, -.1]])
        kinetic_energy = .5 * np.sum(simulator._mass * velocities**2)
        expected = 2 * kinetic_energy * self.hartree_joules / (3 * self.boltzmann)
        np.testing.assert_allclose(simulator._temperature(velocities), expected, rtol=2e-12)
        self.assertIn(('Kelvins', 'Hartrees', ('thermal',)), self.calls)

    def test_langevin_noise_uses_thermal_energy(self):
        thermostat = LangevinThermostat(target_temperature=300., friction=.02, seed=314)
        mass = np.array([[2.], [3.]])
        simulator = SimpleNamespace(_mass=mass, dt=.5)
        positions = np.zeros((2, 3))
        velocities = np.full((2, 3), .01)
        kbt = self.boltzmann * 300 / self.hartree_joules
        normal = np.random.default_rng(314).normal(size=velocities.shape)
        expected = ((1 - .02 * .5) * velocities
                    + np.sqrt(.5) * np.sqrt(2 * .02 * kbt * mass) * normal / mass)
        result_positions, result_velocities = thermostat.apply_thermostat(positions, velocities, simulator)
        np.testing.assert_array_equal(result_positions, positions)
        np.testing.assert_allclose(result_velocities, expected, rtol=2e-12)

    def test_phase_rates_for_energy_wavenumber_kelvin_and_hertz(self):
        cyclic_rates = {
            'Hartrees': self.hartree_joules / self.planck,
            'ElectronVolts': 1.602176634e-19 / self.planck,
            'Wavenumbers': 100 * self.speed_of_light,
            'Kelvins': self.boltzmann / self.planck,
            'Hertz': 1.,
        }
        for frequency_unit, rate in cyclic_rates.items():
            for time_unit, seconds in [('Seconds', 1.), ('PicoSeconds', 1e-12), ('FemtoSeconds', 1e-15)]:
                with self.subTest(frequency_unit=frequency_unit, time_unit=time_unit):
                    actual = NonlinearResponseFunction._freq_conv(frequency_unit, time_unit)
                    np.testing.assert_allclose(actual, 2 * np.pi * rate * seconds, rtol=1e-8)
        frequency_calls = [call for call in self.calls if call[1] == 'Hertz']
        self.assertTrue(all(call[2] == ('wavenumber', 'thermal', 'frequency') for call in frequency_calls))
        self.assertTrue(all(call[2] == () for call in self.calls if call[1] == 'Seconds'))

    def test_response_evaluation_accumulates_one_angular_phase(self):
        response = NonlinearResponseFunction(lambda time: np.exp(-1j * time),
                                             frequency_unit='Hartrees', time_unit='Seconds')
        times = np.array([0., 1e-17, 3e-17])
        expected = np.exp(-1j * times * 2 * np.pi * self.hartree_joules / self.planck)
        np.testing.assert_allclose(response.evaluate_times(times), expected, atol=1e-12)

    def test_sampling_floor_matches_the_phase_convention(self):
        frequency_scales = {
            'Hertz': 100 * self.speed_of_light,
            'Kelvins': 100 * self.planck * self.speed_of_light / self.boltzmann,
        }
        for unit, frequency_scale in frequency_scales.items():
            with self.subTest(unit=unit):
                response = NonlinearResponseFunction(lambda *times: 0, frequency_unit=unit)
                with patch.object(response, '_determine_fft_samples', return_value=(8, .1)) as samples:
                    response.prep_sampling_ranges(np.array([0., 10 * frequency_scale]), 0,
                                                  np.array([0., 20 * frequency_scale]))
                self.assertEqual(samples.call_count, 2)
                for call in samples.call_args_list:
                    np.testing.assert_allclose(call.args[0], 2 * np.pi * 100 * self.speed_of_light * 1e-12,
                                               rtol=1e-8)

    def test_density_matrix_phase_is_consistent_across_unit_representations(self):
        # No pulse kicks: a known initial coherence evolves only by exp(-i E t / hbar).
        gap_hertz = 1000 * 100 * self.speed_of_light
        gaps = {
            'Wavenumbers': 1000.,
            'Hartrees': gap_hertz * self.planck / self.hartree_joules,
            'Kelvins': gap_hertz * self.planck / self.boltzmann,
            'Hertz': gap_hertz,
        }
        filter_space = BasisStateSpace(HarmonicOscillatorProductBasis(1), [(0,), (1,)])
        delays = (.01, .02)
        for unit, gap in gaps.items():
            with self.subTest(unit=unit):
                transitions = {((0,), (1,)): {'frequency': gap, 'transition_moment': [1., 0., 0.]}}
                result = propagate_density_matrix_paths(
                    transitions, frequency_unit=unit, initial_density=np.full((2, 2), .5),
                    pulse_areas=(0., 0., 0.), delay_times=delays,
                    state_filter_opts={'filter_space': filter_space},
                )
                states = [tuple(state) for state in result.total_space.excitations]
                ground, excited = states.index((0,)), states.index((1,))
                expected = .5 * np.exp(1j * 2 * np.pi * gap_hertz * sum(delays) * 1e-12)
                np.testing.assert_allclose(result.density_matrices[-1, ground, excited], expected, atol=1e-8)
                np.testing.assert_allclose(np.trace(result.density_matrices[-1]), 1., atol=1e-12)


if __name__ == '__main__':
    unittest.main()
