"""MACE calculator reuse without importing torch, MACE, or downloading weights."""

import concurrent.futures
import contextlib
import os
import threading
import time
import types
import unittest
from unittest import mock

import numpy as np

from Psience.ElectronicStructure import MACEBackend, MACEPolarBackend
from Psience.Molecools.Evaluator import MACEEnergyEvaluator


class MACEBackendCacheTests(unittest.TestCase):
    def setUp(self):
        self.loads = []

        def load(**settings):
            self.loads.append(settings)
            return types.SimpleNamespace(settings=settings)

        self.load = load
        self.cache_patch = mock.patch.object(MACEBackend, '_mace_cache', {})
        self.cache_patch.start()
        self.addCleanup(self.cache_patch.stop)

    def setup_calculator(self, **settings):
        opts = dict(model='polar-1-l', model_type=self.load, device='cpu')
        return MACEBackend.setup_calculator(**(opts | settings))

    def test_repeated_dtype_settings_reuse_one_calculator(self):
        first = self.setup_calculator(default_dtype='float64')
        for _ in range(50):
            self.assertIs(self.setup_calculator(default_dtype='float64'), first)
        self.assertEqual(len(self.loads), 1)
        self.assertEqual(len(MACEBackend._mace_cache), 1)

    def test_empty_settings_still_reuse(self):
        self.assertIs(self.setup_calculator(), self.setup_calculator())
        self.assertEqual(len(self.loads), 1)

    def test_distinct_configurations_do_not_alias(self):
        configurations = [
            dict(default_dtype='float64'),
            dict(default_dtype='float32'),
            dict(default_dtype='float64', device='cuda:0'),
            dict(default_dtype='float64', model='polar-1-s'),
            dict(default_dtype='float64', model_dir='/synthetic/mace/cache-a'),
            dict(default_dtype='float64', model_dir='/synthetic/mace/cache-b'),
        ]
        calculators = [self.setup_calculator(**settings) for settings in configurations]
        self.assertEqual(len({id(calc) for calc in calculators}), len(configurations))
        for settings, calculator in zip(configurations, calculators):
            self.assertIs(self.setup_calculator(**settings), calculator)

    def test_distinct_model_factories_do_not_alias(self):
        first = self.setup_calculator(default_dtype='float64')
        other_factory = lambda **settings: types.SimpleNamespace(settings=settings)
        second = self.setup_calculator(model_type=other_factory, default_dtype='float64')
        self.assertIsNot(first, second)

    def test_nested_settings_ignore_mapping_order(self):
        first = self.setup_calculator(
            default_dtype='float64',
            compile_options={'layers': [1, 2], 'options': {'mode': 'fast', 'flags': {'a', 'b'}}}
        )
        second = self.setup_calculator(
            compile_options={'options': {'flags': {'b', 'a'}, 'mode': 'fast'}, 'layers': [1, 2]},
            default_dtype='float64'
        )
        self.assertIs(first, second)
        self.assertEqual(len(self.loads), 1)

    def test_setting_types_and_sequence_order_are_preserved(self):
        values = [True, 1, 1., [1, 2], (1, 2), [2, 1]]
        calculators = [self.setup_calculator(option=value) for value in values]
        self.assertEqual(len({id(calc) for calc in calculators}), len(values))
        for value, calculator in zip(values, calculators):
            self.assertIs(self.setup_calculator(option=value), calculator)

    def test_numeric_arrays_reuse_by_value_shape_and_dtype(self):
        array = np.arange(6., dtype=np.float64).reshape(2, 3)
        first = self.setup_calculator(option=array)
        self.assertIs(self.setup_calculator(option=np.asfortranarray(array)), first)
        self.assertIsNot(self.setup_calculator(option=array.reshape(3, 2)), first)
        self.assertIsNot(self.setup_calculator(option=array.astype(np.float32)), first)
        self.assertIsNot(self.setup_calculator(option=array + 1), first)

    def test_mutated_settings_get_a_new_calculator(self):
        settings = {'layers': [1], 'array': np.array([2.])}
        first = self.setup_calculator(option=settings)
        settings['layers'].append(3)
        second = self.setup_calculator(option=settings)
        settings['array'][0] = 4.
        third = self.setup_calculator(option=settings)
        self.assertEqual(len({id(first), id(second), id(third)}), 3)
        self.assertIs(self.setup_calculator(option={'layers': [1], 'array': np.array([2.])}), first)

    def test_unsupported_settings_bypass_without_false_reuse(self):
        class Unhashable:
            __hash__ = None

        cyclic = []
        cyclic.append(cyclic)
        for value in [Unhashable(), np.array([{}], dtype=object), cyclic]:
            with self.subTest(type=type(value)):
                self.assertIsNot(self.setup_calculator(option=value), self.setup_calculator(option=value))
        self.assertEqual(len(MACEBackend._mace_cache), 0)

    def test_failed_load_is_not_cached(self):
        result = object()
        loader = mock.Mock(side_effect=[RuntimeError('failed load'), result, result])
        with self.assertRaisesRegex(RuntimeError, 'failed load'):
            self.setup_calculator(model_type=loader, default_dtype='float64')
        self.assertEqual(len(MACEBackend._mace_cache), 0)
        self.assertIs(self.setup_calculator(model_type=loader, default_dtype='float64'), result)
        self.assertIs(self.setup_calculator(model_type=loader, default_dtype='float64'), result)
        self.assertEqual(loader.call_count, 2)

    def test_concurrent_construction_loads_once(self):
        start = threading.Barrier(8)
        loads = []

        def load(**settings):
            loads.append(settings)
            time.sleep(.03)  # let competing lookups run before this load completes
            return object()

        def construct(_):
            start.wait(timeout=5)
            return self.setup_calculator(model_type=load, default_dtype='float64')

        # Avoid modifying process stdout while testing the old implementation's race.
        with mock.patch('Psience.ElectronicStructure.MLIPBackends.quiet_mode', contextlib.nullcontext):
            with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
                calculators = list(pool.map(construct, range(8)))
        self.assertTrue(all(calc is calculators[0] for calc in calculators))
        self.assertEqual(len(loads), 1)

    def test_polar_backend_and_energy_evaluator_share_default_dtype_model(self):
        mace = types.ModuleType('mace')
        mace.__path__ = []
        calculators = types.ModuleType('mace.calculators')
        calculators.mace_polar = self.load
        mace.calculators = calculators
        with mock.patch.dict('sys.modules', {'mace': mace, 'mace.calculators': calculators}):
            first = MACEPolarBackend(['O', 'H', 'H'], device='cpu')
            second = MACEPolarBackend(['O', 'H', 'H'], device='cpu', charge=1)
            energy = MACEEnergyEvaluator(['O', 'H', 'H'], model='polar', device='cpu', default_dtype='float64')
        self.assertIs(first.model, second.model)
        self.assertIs(first.model, energy.eval)
        self.assertEqual(len(self.loads), 1)

    def test_model_cache_directory_override_is_restored(self):
        with mock.patch.dict(os.environ, {'XDG_CACHE_HOME': '/synthetic/original-cache'}):
            first = self.setup_calculator(model_dir='/synthetic/model-cache', default_dtype='float64')
            self.assertEqual(os.environ['XDG_CACHE_HOME'], '/synthetic/original-cache')
            self.assertIs(self.setup_calculator(model_dir='/synthetic/model-cache', default_dtype='float64'), first)
            self.assertEqual(os.environ['XDG_CACHE_HOME'], '/synthetic/original-cache')
