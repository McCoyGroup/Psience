"""Keep matrix VPT state columns aligned and energy terms on one order schema."""
import contextlib
import io
import re
import tempfile
import unittest
from types import SimpleNamespace
from unittest import mock

import numpy as np

from McUtils.Data import UnitsData
from McUtils.Scaffolding import Logger, NullLogger, NullCheckpointer
from Psience.BasisReps import HarmonicOscillatorProductBasis, BasisStateSpace, StateMaker
from Psience.VPT2.Solver import PerturbationTheorySolver
from Psience.VPT2.Corrections import PerturbationTheoryCorrections
from Psience.VPT2.Wavefunctions import PerturbationTheoryWavefunctions
from Psience.VPT2 import VPTRunner, VPTSystem
from Psience.VPT2.Runner import VPTRuntimeOptions
from Psience.VPT2.DegeneracySpecs import DegenerateMultiStateSpace
from Psience.VPT2.Analyzer import VPTAnalyzerLogParser


class MatrixEnergyTermLoggingTests(unittest.TestCase):
    @staticmethod
    def render(states, terms, real_fmt="{:>12.5f}", state_formatting='quanta'):
        # Exercise the public formatter without constructing unrelated spectra.
        wavefunctions = object.__new__(PerturbationTheoryWavefunctions)
        basis = HarmonicOscillatorProductBasis(len(states[0]))
        wavefunctions.corrs = SimpleNamespace(
            states=BasisStateSpace(basis, states), total_basis=basis,
            all_energy_corrs=terms)
        return wavefunctions.format_energy_corrections_table(real_fmt=real_fmt, state_formatting=state_formatting)

    @staticmethod
    def headings(text):
        return re.findall(r'<0\|dH\(\d+\)\|\d+>', text.splitlines()[0])

    @staticmethod
    def expected_headings(orders):
        return [f'<0|dH({k-i})|{i}>' for k in orders for i in range(k)]

    def test_reported_second_and_third_order_mixed_rows(self):
        h2w = UnitsData.convert('Hartrees', 'Wavenumbers')
        for order in [2, 3]:
            with self.subTest(order=order):
                terms = [
                    [None, None, np.array([3390.96135, -3776.80635])/h2w],
                    [None, np.array([0.]), np.array([3215.83297, -3960.57419])/h2w],
                    [None, np.array([0.]), np.array([3173.39455, -3684.48529])/h2w],
                ]
                if order == 3:
                    terms[0].append(None)
                    terms[1].append(np.zeros(3))
                    terms[2].append(np.zeros(3))
                text = self.render([[0], [1], [2]], terms)
                self.assertEqual(self.headings(text), self.expected_headings(range(1, order+1)))
                rows = [line.split() for line in text.splitlines()[1:]]
                width = 1 + sum(range(1, order+1))
                self.assertTrue(all(len(row) == width for row in rows))
                self.assertEqual(rows[0][1:4], ['--', '3390.96135', '-3776.80635'])
                self.assertEqual(rows[1][1:4], ['0.00000', '3215.83297', '-3960.57419'])
                if order == 3:
                    self.assertEqual(rows[0][4:], ['--'] * 3)
                    self.assertEqual(rows[1][4:], ['0.00000'] * 3)

    def test_higher_orders_not_retained_by_first_state_still_have_headings(self):
        terms = [[None, None, np.array([1., 2.]), None, None],
                 [None, None, np.array([3., 4.]), None, np.array([5., 6., 7., 8.])]]
        text = self.render([[0], [1]], terms)
        self.assertEqual(self.headings(text), self.expected_headings([2, 4]))
        self.assertEqual(text.splitlines()[1].split()[-4:], ['--'] * 4)

    def test_empty_correction_row_does_not_drop_or_reassign_state(self):
        h2w = UnitsData.convert('Hartrees', 'Wavenumbers')
        text = self.render([[0], [2]], [[None, None, None],
            [None, None, np.array([12., 34.])/h2w]])
        rows = [line.split() for line in text.splitlines()[1:]]
        self.assertEqual(rows, [['()', '--', '--'], ['1(2)', '12.00000', '34.00000']])

    def test_zero_order_is_excluded_and_narrow_format_stays_aligned(self):
        h2w = UnitsData.convert('Hartrees', 'Wavenumbers')
        text = self.render([[10], [2]], [
            [np.array([100.]), None, np.array([1., 2.])/h2w],
            [np.array([200.]), None, np.array([3., 4.])/h2w]], real_fmt='{:.1f}')
        self.assertEqual(self.headings(text), self.expected_headings([2]))
        self.assertNotIn('100.', text)
        self.assertNotIn('200.', text)
        rows = [line.split() for line in text.splitlines()[1:]]
        self.assertEqual(rows, [['1(10)', '1.0', '2.0'], ['1(2)', '3.0', '4.0']])

    def test_compact_columns_fit_long_labels_and_keep_numeric_columns_aligned(self):
        states = np.zeros((4, 72), dtype=int)
        states[1, -1] = 1
        states[2, [-3, -1]] = [1, 2]
        states[3, [0, -1]] = [12, 10]
        labels = ['()', '1(1)', '1(2)3(1)', '1(10)72(12)']
        values = np.arange(1, 5, dtype=float)
        wfns = object.__new__(PerturbationTheoryWavefunctions)
        wfns.corrs = SimpleNamespace(total_basis=SimpleNamespace(ndim=72))
        wfns.initial_states = None
        wfns._init_inds = (0,)
        tables = [
            (wfns._format_energies_table(states[0], [1., 2.], states[1:],
                np.column_stack([values[1:], 2*values[1:]])), 2, 12),
            (wfns._format_intensities_table(states, values, 2*values, 3*values, 4*values), 2, 12),
            (wfns._format_operator_table(states, values, 2*values), 1, 12),
            (wfns.format_property_matrices(states, [[values[None, :]]]), 1, 16),
        ]
        width = len(labels[-1]) + 1
        for table, headers, numeric_width in tables:
            with self.subTest(table=table.splitlines()[0]):
                lines = table.splitlines()
                self.assertTrue(any(line[:width].rstrip() == 'State' for line in lines[:headers]))
                self.assertEqual(len(lines[headers:]), len(labels))
                for label, line in zip(labels, lines[headers:]):
                    self.assertEqual(line[:width], label.ljust(width))
                    self.assertEqual(len(line[width:width+numeric_width]), numeric_width)
                    self.assertTrue(line[width:width+numeric_width].strip())
                # Right edges are stable even when the printed numbers have different lengths.
                ends = [[m.end() for m in re.finditer(r'-?\d+\.\d+(?:e[+-]\d+)?|(?<!\S)-(?!\S)',
                    line[width:])] for line in lines[headers:]]
                self.assertTrue(all(row == ends[0] for row in ends))
        correction_table = self.render(states, [[None, np.array([v])] for v in values])
        self.assertEqual([line.split()[0] for line in correction_table.splitlines()[1:]], labels)
        self.assertEqual(len({line.rindex('.') for line in correction_table.splitlines()[1:]}), 1)

    def test_short_label_column_also_fits_state_heading_and_vector_flag(self):
        table = PerturbationTheoryWavefunctions._format_energies_table(
            [0], [1., 2.], [[1]], [[3., 4.]])
        self.assertEqual(table.splitlines()[0][:6], 'State ')
        self.assertTrue(table.splitlines()[2].startswith('()    '))
        vector = self.render([[0], [2]], [[None, np.array([1.])], [None, np.array([2.])]],
            state_formatting='vector')
        self.assertEqual([line.split()[0] for line in vector.splitlines()[1:]], ['0', '2'])

    def test_runner_forwards_state_format_to_all_final_tables(self):
        for state_formatting in ['quanta', 'vector']:
            with self.subTest(state_formatting=state_formatting):
                stream = io.StringIO()
                wfns = mock.Mock()
                wfns.format_dipole_contribs_tables.return_value = ['x', 'y', 'z']
                wfns.format_deperturbed_dipole_contribs_tables.return_value = ['x', 'y', 'z']
                VPTRunner.print_output_tables(wfns, operators={'O': []},
                    logger=Logger(log_file=stream), state_formatting=state_formatting)
                for name in ['format_energy_corrections_table', 'format_energies_table',
                        'format_deperturbed_energies_table', 'format_intensities_table',
                        'format_deperturbed_intensities_table', 'format_dipole_contribs_tables',
                        'format_deperturbed_dipole_contribs_tables']:
                    getattr(wfns, name).assert_called_once_with(state_formatting=state_formatting)
                wfns.format_operator_table.assert_called_once_with({'O': []}, state_formatting=state_formatting)

    def test_multiple_initial_state_headings_and_tables(self):
        states = np.array([[0, 0], [1, 0], [0, 2]])
        wfns = object.__new__(PerturbationTheoryWavefunctions)
        engs = np.arange(3, dtype=float)
        matrix = np.ones((3, 3))
        wfns.corrs = SimpleNamespace(states=SimpleNamespace(excitations=states),
            total_basis=SimpleNamespace(ndim=2), energy_corrs=engs[:, None], degenerate_transf=None)
        wfns.initial_states = states[:2]
        wfns._init_inds = (0, 1)
        ints = {0: np.ones(3), 1: np.ones(3)}
        with mock.patch.object(PerturbationTheoryWavefunctions, 'zero_order_intensities',
                new_callable=mock.PropertyMock, return_value=ints):
            ir = wfns._format_group_intensities_table(engs, ints)
        operator = wfns._format_group_operator_table([(matrix, [[matrix]])])
        contributions = wfns.format_property_matrices(states, [[np.ones((2, 3))]])
        for table in [ir, operator, contributions]:
            self.assertIn('Initial State: ()\n', table)
            self.assertIn('Initial State: 2(1)\n', table)
            self.assertIn('1(2)', table)

    def test_logger_pads_table_headers_and_rows_equally(self):
        table = PerturbationTheoryWavefunctions._format_energies_table(
            [0, 0], [1., 2.], [[0, 1], [1, 2]], [[3., 4.], [5., 6.]])
        wfns = mock.Mock(degenerate_transformation=None)
        wfns.format_energies_table.return_value = table
        stream = io.StringIO()
        logger = Logger(log_file=stream)
        with logger.block(tag='Output'):
            VPTRunner.print_output_tables(wfns, print_energy_corrections=False,
                print_intensities=False, logger=logger)
        printed_rows = [line for line in stream.getvalue().splitlines() if line.startswith('  > ')]
        self.assertEqual([line[4:] for line in printed_rows], table.splitlines())

    @staticmethod
    def matrix_corrections(order, ignore_odd, state_formatting='quanta', logger=None):
        basis = HarmonicOscillatorProductBasis(3)
        states = BasisStateSpace(basis, [[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        def rep(labels, coefficient):
            return basis.representation(*labels, coeffs=coefficient,
                axes=[list(range(len(labels))), list(range(len(labels)))])
        w = np.diag([.04, .071, .103])
        perturbations = [rep(['x', 'x'], w/2) + rep(['p', 'p'], w/2)]
        for k in range(1, order+1):
            coefficient = np.zeros((3,)*(k+2))
            coefficient[(0,)*(k+2)] = 1e-4 / (10**k)
            perturbations.append(rep(['x']*(k+2), coefficient))
        solver = PerturbationTheorySolver(
            perturbations, states, order=order,
            degenerate_states=[np.array([[1, 0, 0], [0, 1, 0]])], handle_strong_couplings=False,
            ignore_odd_order_energies=ignore_odd, check_overlap=False,
            logger=NullLogger() if logger is None else logger, checkpointer=NullCheckpointer(None),
            state_formatting=state_formatting,
            results=NullCheckpointer(None), zero_element_warning=False)
        return solver.apply_VPT()

    def test_matrix_block_and_scalar_routes_agree_on_skipped_orders(self):
        for order in [2, 3, 4]:
            with self.subTest(order=order):
                corrections = self.matrix_corrections(order, True)
                for state_terms, totals in zip(corrections.all_energy_corrs, corrections.energy_corrs):
                    self.assertIsNone(state_terms[0])
                    for k in range(1, order+1):
                        if k % 2:
                            self.assertIsNone(state_terms[k])
                            self.assertEqual(totals[k], 0.)
                        else:
                            self.assertEqual(len(state_terms[k]), k)
                            self.assertAlmostEqual(sum(state_terms[k]), totals[k], places=14)
                text = self.render(corrections.states.excitations, corrections.all_energy_corrs)
                self.assertEqual(self.headings(text), self.expected_headings(range(2, order+1, 2)))

    def test_explicit_odd_orders_keep_all_third_order_terms(self):
        corrections = self.matrix_corrections(3, False)
        for state_terms, totals in zip(corrections.all_energy_corrs, corrections.energy_corrs):
            for k in range(1, 4):
                self.assertEqual(len(state_terms[k]), k)
                self.assertAlmostEqual(sum(state_terms[k]), totals[k], places=14)
        text = self.render(corrections.states.excitations, corrections.all_energy_corrs)
        self.assertEqual(self.headings(text), self.expected_headings([1, 2, 3]))

    def test_deperturbed_diagnostic_fits_labels_and_preserves_frequencies(self):
        basis = HarmonicOscillatorProductBasis(72)
        vectors = np.zeros((3, 72), dtype=int)
        vectors[0, -1] = 1
        vectors[2, [0, -1]] = [12, 2]
        states = BasisStateSpace(basis, vectors)
        energies = np.array([.02, .01, .04])
        before = energies.copy()
        formatter = PerturbationTheorySolver._fmt_depert_engs(states, energies)
        lines = formatter()
        self.assertEqual(lines[0], 'Deperturbed States/Frequencies:')
        self.assertEqual([line.split()[0] for line in lines[1:]], ['1(1)', '()', '1(2)72(12)'])
        self.assertEqual(len({line.rindex('.') for line in lines[1:]}), 1)
        h2w = UnitsData.convert('Hartrees', 'Wavenumbers')
        self.assertEqual([line.split()[-1] for line in lines[1:]],
            [f'{e*h2w:.4f}' for e in [.01, .01, .03]])
        np.testing.assert_array_equal(energies, before)
        excited = states.take_subspace([0, 2])
        self.assertEqual(PerturbationTheorySolver._fmt_depert_engs(excited, energies[[0, 2]])()[0],
            'Deperturbed States/Energies:')
        with mock.patch.object(StateMaker, 'parse_state', side_effect=AssertionError('eager formatting')):
            NullLogger().log_print(None, message_prepper=formatter)

    def test_format_flag_reaches_real_solver_and_resonance_corrections(self):
        results = []
        for mode in ['quanta', 'vector']:
            stream = io.StringIO()
            corrs = self.matrix_corrections(2, True, state_formatting=mode, logger=Logger(log_file=stream))
            results.append(corrs.energies)
            self.assertEqual(corrs.state_formatting, mode)
            subset = corrs.take_subspace(corrs.states.take_subspace([0, 1]))
            self.assertEqual(subset.state_formatting, mode)
            self.assertEqual(VPTRuntimeOptions(state_formatting=mode).solver_opts['state_formatting'], mode)
            output = stream.getvalue()
            if mode == 'quanta':
                self.assertIn('3(1)', output)
                self.assertIn('Deperturbed States/Frequencies:', output)
            else:
                self.assertIn('1 0 0', output)
        np.testing.assert_allclose(results[0], results[1], atol=1e-14)

    def test_both_group_filters_format_warning_states(self):
        basis = HarmonicOscillatorProductBasis(3)
        group = BasisStateSpace(basis, [[1, 0, 0], [0, 1, 0]])
        uncoupled = group.take_subspace([0])
        for filter in [DegenerateMultiStateSpace.default_group_filter, PerturbationTheorySolver.filter_degenerate_group]:
            for mode in ['quanta', 'vector']:
                stream = io.StringIO()
                filter(group, uncoupled_states=uncoupled, state_formatting=mode, logger=Logger(log_file=stream))
                self.assertIn('3(1)' if mode == 'quanta' else '1 0 0', stream.getvalue())

    def test_nielsen_table_defaults_to_occupations_and_supports_vectors(self):
        runner = object.__new__(VPTRunner)
        runner.states = SimpleNamespace(state_list=[[0, 0], [1, 0]])
        runner.get_Nielsen_energies = mock.Mock(return_value=(np.array([.01, .03]), np.array([.011, .029])))
        compact, vector = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(compact):
            runner.print_Nielsen_frequencies()
        with contextlib.redirect_stdout(vector):
            runner.print_Nielsen_frequencies(state_formatting='vector')
        self.assertIn('2(1)', compact.getvalue())
        self.assertNotIn('2(1)', vector.getvalue())

    def test_missing_correction_messages_in_all_storage_routes(self):
        basis = HarmonicOscillatorProductBasis(3)
        states = BasisStateSpace(basis, [[0, 0, 0], [1, 0, 0]])
        for storage in ['dense', 'sparse', 'compact']:
            for mode in ['quanta', 'vector']:
                values = PerturbationTheoryCorrections.create_correction_array(storage, 2, 1, 2, 1024)
                values[:] = np.zeros((2, 2, 2))
                stream = io.StringIO()
                PerturbationTheoryCorrections.create_coupling_matrix(values, states, states, 2, 1,
                    logger=Logger(log_file=stream), state_formatting=mode)
                self.assertIn('state 3(1)' if mode == 'quanta' else 'state [1 0 0]', stream.getvalue())

    def test_run_simple_formats_state_options_and_forwards_flag_before_solve(self):
        vectors = np.array([[0, 0, 0], [1, 0, 0]])
        groups = np.array([vectors])
        for mode in ['quanta', 'vector']:
            stream = io.StringIO()
            logger = Logger(log_file=stream)
            wfns = SimpleNamespace(logger=logger)
            runner = SimpleNamespace(hamiltonian=SimpleNamespace(logger=logger),
                ham_opts=SimpleNamespace(opts={}),
                pt_opts=SimpleNamespace(opts={'degenerate_states': groups,
                    'zero_order_energy_corrections': [(vectors[1], .01)]}),
                runtime_opts=SimpleNamespace(ham_opts={}, solver_opts={}),
                get_wavefunctions=mock.Mock(return_value=wfns), print_tables=mock.Mock())
            space = SimpleNamespace(state_list=vectors, degeneracy_specs=None, degenerate_states=groups)
            with mock.patch.object(VPTRunner, 'construct', return_value=(runner, (None, space, None, None, None))) as construct:
                VPTRunner.run_simple(None, vectors, state_formatting=mode)
            self.assertEqual(construct.call_args.kwargs['state_formatting'], mode)
            self.assertEqual(runner.print_tables.call_args.kwargs['state_formatting'], mode)
            output = stream.getvalue()
            self.assertIn('3(1)' if mode == 'quanta' else '[1 0 0]', output)
            if mode == 'quanta':
                self.assertNotIn('[1 0 0]', output)

    def test_construct_preserves_format_for_state_preparation_and_solver_factories(self):
        system = object.__new__(VPTSystem)
        system.mode_selection = [0, 1]
        system.use_local_modes = False
        system.mol = SimpleNamespace(normal_modes=SimpleNamespace(modes=SimpleNamespace(freqs=np.array([.01, .02]))))
        for mode in ['quanta', 'vector']:
            runner, _ = VPTRunner.construct(system, [[0, 0], [1, 0]], state_formatting=mode,
                degeneracy_specs={'groups': [[[1, 0], [0, 1]]]})
            self.assertEqual(runner.states.state_formatting, mode)
            runner._ham = mock.Mock()
            runner.get_solver()
            runner.get_wavefunctions()
            self.assertEqual(runner._ham.get_solver.call_args.kwargs['state_formatting'], mode)
            self.assertEqual(runner._ham.get_wavefunctions.call_args.kwargs['state_formatting'], mode)

    def test_strong_coupling_report_uses_occupations_and_keeps_vector_option(self):
        basis = HarmonicOscillatorProductBasis(3)
        states = BasisStateSpace(basis, [[0, 0, 0], [1, 0, 0]])
        corrs = object.__new__(PerturbationTheoryCorrections)
        corrs.states, corrs.total_basis = states, states
        couplings = {int(states.indices[1]): [states]}
        compact = corrs.format_strong_couplings_report(couplings)
        self.assertIn('state:  3(1)', compact)
        self.assertIn('()', compact)
        vector = corrs.format_strong_couplings_report(couplings, state_formatting='vector')
        self.assertNotIn('3(1)', vector)
        self.assertIn('1   0   0', vector)

    def test_log_readers_accept_varying_occupation_labels_and_legacy_vectors(self):
        for parser_type in [VPTAnalyzerLogParser.EnergiesBlockParser,
                VPTAnalyzerLogParser.SpectrumBlockParser, VPTAnalyzerLogParser.TransitionMomentBlockParser]:
            parser = parser_type('')
            for state in ['()', '1(1)', '1(2)72(12)']:
                label, values = parser.handle_block_line('()', state.ljust(13) + '  1.00000  -2.00000')
                self.assertEqual(label, state)
                self.assertEqual(values, ['1.00000', '-2.00000'])
            legacy = parser_type('')
            label, values = legacy.handle_block_line('0 0 0', '1 0 0    1.00000  -2.00000')
            self.assertEqual(label, '1 0 0')
            self.assertEqual(values, ['1.00000', '-2.00000'])

    def test_compact_logged_tables_round_trip_through_analyzer(self):
        states = np.array([[0, 0, 0], [1, 0, 0], [1, 0, 2]])
        wavefunctions = object.__new__(PerturbationTheoryWavefunctions)
        wavefunctions.corrs = SimpleNamespace(total_basis=SimpleNamespace(ndim=3))
        wavefunctions.initial_states = states[:1]
        wavefunctions._init_inds = (0,)
        contributions = wavefunctions.format_property_matrices(states, [[np.array([[1., 2., 3.]])]])
        energies = wavefunctions._format_energies_table(states[0], [1., 2.], states[1:], [[3., 4.], [5., 6.]])
        intensities = 'Initial State: ()\n' + wavefunctions._format_intensities_table(
            states[1:], [4., 6.], [10., 20.], [3., 5.], [11., 21.])
        wfns = mock.Mock(degenerate_transformation=None)
        wfns.format_energies_table.return_value = energies
        wfns.format_dipole_contribs_tables.return_value = [contributions] * 3
        wfns.format_intensities_table.return_value = intensities
        stream = io.StringIO()
        logger = Logger(log_file=stream)
        with logger.block(tag='Starting Perturbation Theory Runner'):
            VPTRunner.print_output_tables(wfns, print_energy_corrections=False, logger=logger)
        with tempfile.NamedTemporaryFile(mode='w+') as log_file:
            log_file.write(stream.getvalue())
            log_file.flush()
            parser = VPTAnalyzerLogParser(log_file.name)
            labels, harm = parser.harmonic_energies
            self.assertEqual(labels, ['()', '3(1)', '1(2)3(1)'])
            np.testing.assert_allclose(harm, [1., 3., 5.])
            self.assertIn('1(2)3(1)', str(parser.spectra))
            self.assertIn('1(2)3(1)', str(parser.transition_moment_corrections))


if __name__ == '__main__':
    unittest.main()
