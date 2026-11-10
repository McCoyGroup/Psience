"""
Keeps the original (pre-September 2026) analytic VPT flows reachable.

Reference energies were produced by Psience f1669cc6 / McUtils 40956d88
(`AnalyticVPTRunner.run_simple`, default options), the last state before
LLM-assisted development of VPT2 began.
"""
import os, unittest
import numpy as np

from Psience.VPT2 import AnalyticVPTRunner
from Psience.VPT2.Analytic import PerturbationTheoryExpressionEvaluator
from Psience.VPT2.Common import Diagnostics

TEST_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'TestData')
h2w = 219474.6313632

JULY_HOH_POLYAD = np.array([4605.9531, 8350.6953, 12021.3409, 8227.9324, 11806.3164,
                            11737.8357, 6178.6604, 9900.3319, 9780.5667, 7723.344])
JULY_OCHH_BLOCKS = np.array([5785.9524, 8635.4185, 8606.502, 11351.5096, 10296.6796,
                             8753.6573, 11534.0148])


class VPT2LegacyFlowTests(unittest.TestCase):

    def run_hoh(self, **opts):
        return AnalyticVPTRunner.run_simple(
            os.path.join(TEST_DIR, 'HOH_freq.fchk'), 2,
            degeneracy_specs={'polyads': [[[0, 0, 1], [0, 1, 0]]]},
            logger=False, calculate_intensities=False, **opts
        )

    def test_filtered_degeneracy_handling_reproduces_original_energies(self):
        corrs = self.run_hoh(degeneracy_handling='filtered')
        np.testing.assert_allclose(corrs.energies * h2w, JULY_HOH_POLYAD, atol=1e-3)
        # the 001/010 Fermi coupling of the effective Hamiltonian
        ham = [np.asarray(h) for h in corrs.degenerate_hamiltonians if np.shape(h)[0] == 3][0]
        self.assertAlmostEqual(ham[0, 1] * h2w, -78.8757, places=3)

    def test_filtered_and_projected_agree_at_second_order(self):
        filtered = self.run_hoh(degeneracy_handling='filtered')
        projected = self.run_hoh()
        np.testing.assert_allclose(filtered.energies, projected.energies, atol=1e-9)

    def test_filtered_explicit_state_blocks(self):
        corrs = AnalyticVPTRunner.run_simple(
            os.path.join(TEST_DIR, 'OCHH_freq.fchk'),
            [[0, [[0, 0, 0, 0, 0, 1], [0, 1, 0, 1, 0, 0], [0, 0, 0, 1, 1, 0], [0, 0, 0, 0, 1, 0]]],
             [[0, 0, 0, 0, 1, 0], [[0, 0, 0, 0, 1, 1], [0, 1, 0, 1, 1, 0]]]],
            degeneracy_specs=[[[0, 0, 0, 0, 0, 1], [0, 1, 0, 1, 0, 0]]],
            degeneracy_handling='filtered', logger=False, calculate_intensities=False
        )
        np.testing.assert_allclose(corrs.energies * h2w, JULY_OCHH_BLOCKS, atol=1e-3)

    def test_unknown_degeneracy_handling_is_rejected(self):
        with self.assertRaises(ValueError):
            self.run_hoh(degeneracy_handling='nonsense')

    def test_statistics_are_not_recorded_by_default(self):
        self.assertFalse(Diagnostics.enabled)
        PerturbationTheoryExpressionEvaluator._last_indexed_evaluation_stats = None
        self.run_hoh()
        self.assertIsNone(PerturbationTheoryExpressionEvaluator.get_last_indexed_evaluation_stats())
        with Diagnostics.collect():
            self.run_hoh()
        self.assertIsNotNone(PerturbationTheoryExpressionEvaluator.get_last_indexed_evaluation_stats())


if __name__ == '__main__':
    unittest.main()
