"""Logging and numerical regressions for both analytic resonance routes."""
import contextlib
import io
import tempfile
import unittest
import warnings
from unittest import mock

import numpy as np

from McUtils.Scaffolding import Logger
from Psience.VPT2.Analytic import (
    AnalyticPerturbationTheorySolver, HamiltonianExpansionTerm,
    PerturbationTheoryEvaluator, PerturbationTheoryExpressionEvaluator,
)
from Psience.VPT2.Conditional import SameStateCondition, ZeroTerm


class AnalyticResonanceLoggingTests(unittest.TestCase):
    def setUp(self):
        AnalyticPerturbationTheorySolver.clear_caches()

    @staticmethod
    def capture():
        stream = io.StringIO()
        return Logger(log_file=stream), stream

    @staticmethod
    def evaluator(conditional, logger):
        w = np.array([.04])
        h0 = AnalyticPerturbationTheorySolver.from_order(
            2, internals=False, logger=logger).hamiltonian_expansion[0]
        h1 = HamiltonianExpansionTerm([['x']], order=1, identities=[0], logger=logger)
        solver = AnalyticPerturbationTheorySolver(
            [h0, h1, ZeroTerm(logger=logger)], logger=logger,
            resonance_condition=SameStateCondition() if conditional else None)
        return PerturbationTheoryEvaluator(
            solver, [[np.diag(w)/2, np.diag(w)/2], [np.array([.004])], [0.]], freqs=w)

    @staticmethod
    def options(conditional):
        return {} if conditional else {'degenerate_blocks': [[[0], [1]]]}

    def assert_timings(self, text):
        self.assertIn('Getting corrections at order 0', text)
        self.assertIn('Getting corrections at order 2', text)
        self.assertIn('evaluating...', text)
        self.assertRegex(text, r'took \d+\.\d{3}s')

    def test_conditional_energy_restores_standard_timers_and_values(self):
        logger, stream = self.capture()
        evaluator = self.evaluator(True, logger)
        values = evaluator.get_energy_corrections([[0], [1]], order=2)
        self.assert_timings(stream.getvalue())
        self.assertIn('Building evaluator', stream.getvalue())
        self.assertIn('Building conditional transformation at order 1', stream.getvalue())
        np.testing.assert_allclose(values, [[.02, .06], [-.0002, -.0002]], atol=1e-14)

    def test_conditional_overlap_has_build_and_evaluation_timers(self):
        logger, stream = self.capture()
        evaluator = self.evaluator(True, logger)
        values = evaluator.get_overlap_corrections([[0], [1]], order=2)
        self.assertEqual(values.shape, (2, 3))
        self.assertIn('Getting overlap corrections at order 2', stream.getvalue())
        self.assertIn('Building evaluator', stream.getvalue())
        self.assertIn('evaluating...', stream.getvalue())
        quiet = self.evaluator(True, False).get_overlap_corrections([[0], [1]], order=2)
        np.testing.assert_allclose(values, quiet, atol=1e-14)

    def test_compiled_expression_preparation_and_execution_have_timers(self):
        logger, stream = self.capture()
        evaluator = self.evaluator(True, logger)
        with tempfile.TemporaryDirectory() as folder:
            evaluator.expressions_folder = folder
            values = evaluator.get_energy_corrections([[0], [1]], order=2)
            self.assertIn('Preparing compiled expression...', stream.getvalue())
            self.assertIn('evaluating compiled expression...', stream.getvalue())
        np.testing.assert_allclose(values, [[.02, .06], [-.0002, -.0002]], atol=1e-14)

    def test_projected_energy_times_new_contraction_stages(self):
        logger, stream = self.capture()
        evaluator = self.evaluator(False, logger)
        values = evaluator.get_energy_corrections([[0], [1]], order=2, **self.options(False))
        self.assert_timings(stream.getvalue())
        for tag in ['Building projected block (2 states) at order 1',
                    'Contracting block normalization...', 'Contracting block Hamiltonian...',
                    'Building projected resolvents and wavefunctions...']:
            self.assertIn(tag, stream.getvalue())
        np.testing.assert_allclose(values[0], [.02, .06])
        quiet = self.evaluator(False, False).get_energy_corrections(
            [[0], [1]], order=2, **self.options(False))
        np.testing.assert_allclose(values, quiet, atol=1e-14)

    def test_projected_overlap_retains_state_by_order_layout(self):
        logger, stream = self.capture()
        evaluator = self.evaluator(False, logger)
        values = evaluator.get_overlap_corrections([[0], [1]], order=2, **self.options(False))
        self.assertEqual(values.shape, (2, 3))
        self.assert_timings(stream.getvalue())
        quiet = self.evaluator(False, False).get_overlap_corrections(
            [[0], [1]], order=2, **self.options(False))
        np.testing.assert_allclose(values, quiet, atol=1e-14)

    def test_energy_per_call_logger_overrides_solver_and_can_be_disabled(self):
        for conditional in [True, False]:
            with self.subTest(conditional=conditional):
                logger, solver_output = self.capture()
                override, call_output = self.capture()
                evaluator = self.evaluator(conditional, logger)
                options = self.options(conditional)
                values = evaluator.get_energy_corrections([[0], [1]], order=2,
                    logger=override, **options)
                self.assert_timings(call_output.getvalue())
                self.assertEqual(solver_output.getvalue(), '')
                before = call_output.getvalue()
                quiet = evaluator.get_energy_corrections([[0], [1]], order=2,
                    logger=False, **options)
                np.testing.assert_allclose(values, quiet, atol=1e-14)
                self.assertEqual(call_output.getvalue(), before)
                self.assertEqual(solver_output.getvalue(), '')

    def test_operator_override_and_verbose_reach_expression_evaluators(self):
        original = PerturbationTheoryExpressionEvaluator.evaluate
        for conditional in [True, False]:
            with self.subTest(conditional=conditional):
                logger, solver_output = self.capture()
                override, call_output = self.capture()
                evaluator = self.evaluator(conditional, logger)
                seen = []
                def evaluate(expression, *args, **kwargs):
                    seen.append((expression.logger, kwargs.get('verbose')))
                    return original(expression, *args, **kwargs)
                with mock.patch.object(PerturbationTheoryExpressionEvaluator, 'evaluate', evaluate):
                    result = evaluator.get_operator_corrections(
                        [[np.array([1.]), np.zeros((1, 1)), np.zeros((1, 1, 1))]],
                        [[[0], [[1], [2]]]], order=2, logger=override,
                        check_single=False, verbose=True, **self.options(conditional))[0]
                self.assertTrue(seen)
                self.assertTrue(all(log is override and verbose is True for log, verbose in seen))
                self.assert_timings(call_output.getvalue())
                self.assertIn('Gathering indexed coefficients', call_output.getvalue())
                self.assertIn('Routed', call_output.getvalue())
                self.assertIn('Evaluated indexed blocks; took', call_output.getvalue())
                self.assertEqual(solver_output.getvalue(), '')
                quiet = evaluator.get_operator_corrections(
                    [[np.array([1.]), np.zeros((1, 1)), np.zeros((1, 1, 1))]],
                    [[[0], [[1], [2]]]], order=2, logger=False, check_single=False,
                    **self.options(conditional))[0]
                for n in range(3):
                    np.testing.assert_allclose(result.corrections[0][n], quiet.corrections[0][n], atol=1e-14)

    def test_hamiltonian_override_reaches_transformations_and_evaluation(self):
        for conditional in [True, False]:
            with self.subTest(conditional=conditional):
                logger, solver_output = self.capture()
                override, call_output = self.capture()
                evaluator = self.evaluator(conditional, logger)
                result = evaluator.get_reexpressed_hamiltonian([[0], [1]], order=2,
                    logger=override, **self.options(conditional))
                self.assert_timings(call_output.getvalue())
                self.assertEqual(solver_output.getvalue(), '')
                self.assertTrue(np.isfinite(result.corrections[0][2]).all())

    def test_cached_blocks_use_current_logger_for_uncached_elements(self):
        first, first_output = self.capture()
        second, second_output = self.capture()
        evaluator = self.evaluator(False, False)
        evaluator.get_energy_corrections([[0], [1]], order=2,
            logger=first, **self.options(False))
        before = first_output.getvalue()
        result = evaluator.get_operator_corrections(
            [[np.array([1.]), np.zeros((1, 1)), np.zeros((1, 1, 1))]],
            [[[0], [[2]]]], order=2, logger=second, verbose=True,
            check_single=False, **self.options(False))[0]
        self.assertEqual(first_output.getvalue(), before)
        self.assertIn('Building projected evaluator', second_output.getvalue())
        self.assertNotIn('Building projected block (2 states)', second_output.getvalue())
        self.assertTrue(np.isfinite(result.corrections[0][2]).all())

    def test_wavefunctions_and_intermediate_normalization_have_timings(self):
        for conditional in [True, False]:
            with self.subTest(conditional=conditional):
                logger, stream = self.capture()
                evaluator = self.evaluator(conditional, logger)
                evaluator.solver.intermediate_normalization = True
                result = evaluator.get_full_wavefunction_corrections(
                    [[[0], [[0], [1], [2]]]], order=2, **self.options(conditional))
                self.assert_timings(stream.getvalue())
                self.assertIn('intermediate normalization', stream.getvalue())
                self.assertTrue(np.isfinite(result.corrections[0][2]).all())

    def test_target_orders_preserve_skips_and_standard_timers(self):
        for conditional in [True, False]:
            with self.subTest(conditional=conditional):
                logger, stream = self.capture()
                evaluator = self.evaluator(conditional, logger)
                result = evaluator.get_wavefunction_corrections([[[0], [[2]]]],
                    order=2, target_orders=[2], **self.options(conditional))
                if conditional:
                    # The native result formatter fills skipped orders with zeros.
                    np.testing.assert_array_equal(result.corrections[0][0], 0.)
                    np.testing.assert_array_equal(result.corrections[0][1], 0.)
                else:
                    self.assertIsNone(result.corrections[0][0])
                    self.assertIsNone(result.corrections[0][1])
                self.assertIn('Getting corrections at order 2', stream.getvalue())
                self.assertNotIn('Getting corrections at order 0', stream.getvalue())
                self.assertNotIn('Getting corrections at order 1', stream.getvalue())

    def test_logging_disabled_and_legacy_energy_remain_supported(self):
        output = io.StringIO()
        with contextlib.redirect_stdout(output), contextlib.redirect_stderr(output):
            for conditional in [True, False]:
                evaluator = self.evaluator(conditional, False)
                evaluator.get_energy_corrections([[0], [1]], order=2, **self.options(conditional))
        self.assertEqual(output.getvalue(), '')
        logger, stream = self.capture()
        legacy = self.evaluator(False, logger).get_energy_corrections([[0], [1]], order=2)
        self.assert_timings(stream.getvalue())
        np.testing.assert_allclose(legacy, [[.02, .06], [-.0002, -.0002]], atol=1e-14)


    def test_existing_denominator_warnings_reach_normal_log_without_verbose(self):
        E = PerturbationTheoryExpressionEvaluator
        logger, output = self.capture()
        with warnings.catch_warnings(), mock.patch.object(E,'_eval_poly',return_value=np.array([[2.]])):
            warnings.simplefilter('ignore')
            value, _, _ = E._eval_subcontrib(
                ((0,1),), np.array([[1e-30]]), {}, {},
                None, None, None, None, None, None, False, logger,
                zero_cutoff=1e-24)
        np.testing.assert_allclose(value,[[2e30]])
        text = output.getvalue()
        self.assertIn('WARNING: energy denominator below zero_cutoff',text)
        self.assertIn('term survived numerator pruning',text)

    def test_existing_block_warning_and_pruned_zero_remain_visible(self):
        E = PerturbationTheoryExpressionEvaluator
        logger, output = self.capture()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            value = E._divide_polynomial_block_by_energy(
                np.array([[[0.,2.]]]),np.array([[0.,1e-30]]),1e-24,logger=logger)
        np.testing.assert_allclose(value,[[[0.,2e30]]])
        self.assertIn('energy denominator below zero_cutoff',output.getvalue())
        self.assertIn('term survived numerator pruning',output.getvalue())

    def test_projected_warning_recovers_actual_reference_and_intermediate_pair(self):
        from Psience.VPT2.Projected import ProjectionFactor
        E = PerturbationTheoryExpressionEvaluator
        reference = [1,0]
        factor = ProjectionFactor('R',0)
        bindings = ((((1,0),),(.01,.0100000000001),(1,0),False),)
        # Concrete mode permutation reverses the two local delta positions.
        rows = ((factor,1,0),)
        context = dict(operator='fixture',_state_pair_context=lambda:(
            np.array([[0,0]]),np.array([[1,0]]),np.array([[.0100000000001,.01]]),None,bindings))
        details = E._format_small_denominator_context(context,np.array([0]),np.array([1e-30]),energy_changes=rows)
        pair = details['denominator_state_pairs'][0]['state_pairs'][0]
        self.assertEqual(pair['left_state'],(1,0))
        self.assertEqual(pair['right_state'],(0,1))
        self.assertEqual(pair['energy_difference'],'E(left)-E(right)')
        self.assertNotIn('_state_pair_context',details['evaluation'])
        logger, output = self.capture()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            E._divide_polynomial_by_energy(np.array([[2.]]),np.array([1e-30]),1e-24,
                energy_changes=rows,evaluation_context=context,logger=logger)
        self.assertIn("'left_state': (1, 0)",output.getvalue())
        self.assertIn("'right_state': (0, 1)",output.getvalue())

    def test_conditional_warning_recovers_shifted_endpoint_states(self):
        from Psience.VPT2.Conditional import ConditionEndpoint
        E = PerturbationTheoryExpressionEvaluator
        condition = SameStateCondition()
        rows = ((ConditionEndpoint(condition,'resolvent',False,0),0),
                (ConditionEndpoint(condition,'resolvent',False,1),1))
        result = E._denominator_state_pairs(rows,np.array([0]),lambda:(
            np.array([[2]]),np.array([[0]]),np.array([[.01]]),np.array([1])))
        pair = result[0]['state_pairs'][0]
        self.assertEqual(pair['left_state'],(3,))
        self.assertEqual(pair['right_state'],(4,))
        self.assertAlmostEqual(pair['gap_hartree'],-.01)

    def test_relative_energy_rows_do_not_invent_unique_intermediate_states(self):
        result = PerturbationTheoryExpressionEvaluator._denominator_state_pairs(
            ((0,1,-1),),np.array([0]),lambda:(np.array([[1,0]]),np.array([[0,1]]),np.array([[.01,.01]]),None))
        self.assertEqual(result[0]['state_pairs'],[])
        self.assertEqual(result[0]['source_state'],(1,0))
        self.assertEqual(result[0]['unresolved_factors'][0]['energy_change'],(1,-1))

    def test_requested_zero_cutoff_default(self):
        self.assertEqual(PerturbationTheoryExpressionEvaluator.default_zero_cutoff,1e-24)


if __name__ == '__main__':
    unittest.main()
