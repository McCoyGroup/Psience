"""Keep matrix VPT energy term headings and rows on one order schema."""
import re
import unittest
from types import SimpleNamespace

import numpy as np

from McUtils.Data import UnitsData
from McUtils.Scaffolding import NullLogger, NullCheckpointer
from Psience.BasisReps import HarmonicOscillatorProductBasis, BasisStateSpace
from Psience.VPT2.Solver import PerturbationTheorySolver
from Psience.VPT2.Wavefunctions import PerturbationTheoryWavefunctions


class MatrixEnergyTermLoggingTests(unittest.TestCase):
    @staticmethod
    def render(states, terms, real_fmt="{:>12.5f}"):
        # Exercise the public formatter without constructing unrelated spectra.
        wavefunctions = object.__new__(PerturbationTheoryWavefunctions)
        basis = HarmonicOscillatorProductBasis(len(states[0]))
        wavefunctions.corrs = SimpleNamespace(
            states=BasisStateSpace(basis, states), total_basis=basis,
            all_energy_corrs=terms)
        return wavefunctions.format_energy_corrections_table(real_fmt=real_fmt)

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
        self.assertEqual(rows, [['0', '--', '--'], ['2', '12.00000', '34.00000']])

    def test_zero_order_is_excluded_and_narrow_format_stays_aligned(self):
        h2w = UnitsData.convert('Hartrees', 'Wavenumbers')
        text = self.render([[10], [2]], [
            [np.array([100.]), None, np.array([1., 2.])/h2w],
            [np.array([200.]), None, np.array([3., 4.])/h2w]], real_fmt='{:.1f}')
        self.assertEqual(self.headings(text), self.expected_headings([2]))
        self.assertNotIn('100.', text)
        self.assertNotIn('200.', text)
        rows = [line.split() for line in text.splitlines()[1:]]
        self.assertEqual(rows, [['10', '1.0', '2.0'], ['2', '3.0', '4.0']])

    @staticmethod
    def matrix_corrections(order, ignore_odd):
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
            logger=NullLogger(), checkpointer=NullCheckpointer(None),
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


if __name__ == '__main__':
    unittest.main()
