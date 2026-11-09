"""Exercise complete analytic timing scopes without changing numerical work."""
import io
import tempfile
import unittest
from unittest import mock

import numpy as np

from McUtils.Scaffolding import Logger
from Psience.VPT2 import AnalyticVPTRunner
from Psience.VPT2.Analytic import (
    AnalyticPerturbationTheorySolver, HamiltonianExpansionTerm,
    PerturbationTheoryEvaluator,
)
from Psience.VPT2.Analytic import SameStateCondition, ZeroTerm


class AnalyticStageTimingTests(unittest.TestCase):
    def setUp(self):
        AnalyticPerturbationTheorySolver.clear_caches()

    @staticmethod
    def capture(level=None):
        stream = io.StringIO()
        return Logger(log_file=stream, log_level=level), stream

    @staticmethod
    def evaluator(route, logger):
        w = np.array([.04])
        h0 = AnalyticPerturbationTheorySolver.from_order(
            2, internals=False, logger=logger).hamiltonian_expansion[0]
        h1 = HamiltonianExpansionTerm([['x']], order=1, identities=[0], logger=logger)
        solver = AnalyticPerturbationTheorySolver([h0, h1, ZeroTerm(logger=logger)],
            logger=logger, resonance_condition=SameStateCondition() if route=='conditional' else None)
        evaluator = PerturbationTheoryEvaluator(solver,
            [[np.diag(w)/2, np.diag(w)/2], [np.array([.004])], [0.]], freqs=w)
        return evaluator

    @staticmethod
    def runner(logger):
        evaluator = AnalyticStageTimingTests.evaluator('conditional', logger)
        dipole = [[0., np.array([scale]), np.zeros((1,1)), np.zeros((1,1,1))]
                  for scale in [1., 2., 3.]]
        return AnalyticVPTRunner(evaluator, order=2, dipole_expansion=dipole)

    def test_energy_totals_cover_ordinary_conditional_and_projected_routes(self):
        for route in ['ordinary', 'conditional', 'projected']:
            with self.subTest(route=route):
                logger, stream = self.capture()
                evaluator = self.evaluator(route, logger)
                options = {'degenerate_blocks': [[[0], [1]]]} if route=='projected' else {}
                values = evaluator.get_energy_corrections([[0], [1]], order=2, **options)
                text = stream.getvalue()
                self.assertIn('Energy corrections (total)', text)
                self.assertRegex(text, r'took \d+\.\d{3}s')
                np.testing.assert_allclose(values[0], [.02, .06])
                quiet = evaluator.get_energy_corrections([[0], [1]], order=2, logger=False, **options)
                np.testing.assert_allclose(values, quiet, atol=1e-14)
                if route == 'ordinary':
                    # The existing scalar term builder uses its solver logger.
                    # Per-call False must still suppress the newly added timers.
                    added = stream.getvalue()[len(text):]
                    self.assertNotIn('Energy corrections (total)', added)
                else:
                    self.assertEqual(stream.getvalue(), text)

    def test_grouping_logger_override_and_detail_level(self):
        default, default_output = self.capture()
        override, output = self.capture()
        evaluator = self.evaluator('ordinary', default)
        changes = evaluator.get_diff_map([[[0], [[1], [2]]]], logger=override)
        self.assertEqual(set(changes), {(1,), (2,)})
        self.assertIn('Grouping state changes and permutations', output.getvalue())
        self.assertEqual(default_output.getvalue(), '')
        before = output.getvalue()
        quiet = evaluator.get_diff_map([[[0], [[1], [2]]]], logger=override,
                                       log_level=Logger.LogLevel.Debug)
        self.assertEqual(output.getvalue(), before)
        self.assertEqual(set(changes), set(quiet))

    def test_runner_state_preparation_and_complete_run(self):
        logger, output = self.capture()
        runner = self.runner(logger)
        result = runner.run_VPT([[0], [1]], order=2, calculate_intensities=False,
            handle_degeneracies=False, clear_caches=False, force_return_on_crash=False)
        self.assertTrue(np.isfinite(result.energy_corrections).all())
        text = output.getvalue()
        for label in ['Analytic VPT run (total)', 'Preparing analytic states',
                      'Energy correction call (total)', 'Energy corrections (total)']:
            self.assertIn(label, text)

    def test_dipole_axes_keep_one_evaluation_batch_and_current_logger(self):
        default, default_output = self.capture()
        override, output = self.capture()
        runner = self.runner(default)
        with mock.patch.object(runner.eval, 'get_operator_corrections',
                               wraps=runner.eval.get_operator_corrections) as evaluate:
            values = runner.get_transition_moment_corrections([[[0], [[1], [2]]]],
                order=2, axes=(axis for axis in [2,0]), logger=override)
            self.assertEqual(evaluate.call_count, 1)
            self.assertEqual(len(evaluate.call_args.args[0]), 2)
        text = output.getvalue()
        for label in ['Preparing transition moment axis z', 'Preparing transition moment axis x',
                      'Transition moment corrections (axes: z, x; batched)',
                      'Transition moment correction call (total)', 'Preparing operator expansions',
                      'Preparing operator coefficients', 'Regrouping operator correction results']:
            self.assertIn(label, text)
        self.assertEqual(default_output.getvalue(), '')
        quiet = runner.get_transition_moment_corrections([[[0], [[1], [2]]]], order=2,
                                                       axes=[2,0], logger=False)
        self.assertEqual(output.getvalue(), text)
        self.assertEqual(default_output.getvalue(), '')
        for component in range(2):
            for order in range(3):
                np.testing.assert_allclose(values[component][0][order], quiet[component][0][order])
                np.testing.assert_allclose(values[0][0][order], 3*values[1][0][order], atol=1e-14)

    def test_remaining_correction_totals_and_result_assembly(self):
        logger, output = self.capture()
        evaluator = self.evaluator('conditional', logger)
        states = [[[0], [[0], [1], [2]]]]
        evaluator.get_overlap_corrections([[0], [1]], order=2)
        evaluator.get_matrix_corrections(states, order=2)
        evaluator.get_full_wavefunction_corrections(states, order=2)
        evaluator.get_wavefunction_corrections(states, order=2)
        evaluator.get_reexpressed_hamiltonian([[0], [1]], order=2)
        text = output.getvalue()
        for label in ['Overlap corrections (total)', 'Matrix corrections (total)',
                      'State-by-state corrections (total)', 'Full wavefunction corrections (total)',
                      'Wavefunction corrections (total)', 'Effective Hamiltonian corrections (total)',
                      'Assembling and sorting correction results', 'Grouping state changes and permutations']:
            self.assertIn(label, text)

    def test_positional_logger_false_suppresses_new_timers(self):
        logger, output = self.capture()
        evaluator = self.evaluator('conditional', logger)
        # The logger is positional after states/order/expansions/freqs/cutoff/degs/verbose.
        values = evaluator.get_energy_corrections([[0]], 2, None, None, None, None, False, False)
        self.assertEqual(output.getvalue(), '')
        np.testing.assert_allclose(values, [[.02], [-.0002]], atol=1e-14)

    def test_timers_preserve_debug_verbosity_and_restore_after_errors(self):
        logger, output = self.capture(Logger.LogLevel.MoreDebug)
        evaluator = self.evaluator('ordinary', logger)
        before = (logger.verbosity, logger.block_level, logger.padding)
        with self.assertRaises((TypeError, ValueError, IndexError)):
            evaluator.get_diff_map([[[0], None]])
        self.assertEqual((logger.verbosity, logger.block_level, logger.padding), before)
        evaluator.get_diff_map([[[0], [[1]]]])
        self.assertEqual((logger.verbosity, logger.block_level, logger.padding), before)
        self.assertIn('Grouping state changes and permutations', output.getvalue())

    def test_compiled_correction_is_included_in_total(self):
        logger, output = self.capture()
        evaluator = self.evaluator('conditional', logger)
        with tempfile.TemporaryDirectory() as folder:
            evaluator.expressions_folder = folder
            values = evaluator.get_energy_corrections([[0], [1]], order=2)
        text = output.getvalue()
        self.assertIn('Energy corrections (total)', text)
        self.assertIn('Preparing compiled expression...', text)
        self.assertIn('evaluating compiled expression...', text)
        np.testing.assert_allclose(values, [[.02, .06], [-.0002, -.0002]], atol=1e-14)


if __name__ == '__main__':
    unittest.main()
