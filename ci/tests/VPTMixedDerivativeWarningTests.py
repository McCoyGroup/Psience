import unittest
import numpy as np
from Psience.VPT2.Terms import PotentialTerms


class RecordingLogger:
    def __init__(self):
        self.messages = []

    def log_print(self, template, **values):
        self.messages.append(template.format(**values))


class VPTMixedDerivativeWarningTests(unittest.TestCase):
    def run_handling(self, raw, mode="analytical", **opts):
        logger = RecordingLogger()
        result = PotentialTerms._symmetrize_mixed_derivatives(
            raw, mode, 2, logger=logger, term_id="v4_cart",
            warning_diff=opts.pop("warning_diff", 0.5), **opts
        )
        return result, logger.messages

    def test_resolution_message_precedes_warning_and_names_handling(self):
        raw = np.zeros((2,) * 4)
        raw[0, 0, 1, 1], raw[1, 1, 0, 0] = 2, 8
        for mode, expected in [("analytical", 2), ("numerical", 8), ("averaged", 5)]:
            with self.subTest(mode=mode):
                result, messages = self.run_handling(raw, mode)
                self.assertIn("before symmetry completion", messages[0])
                self.assertIn("will be resolved using " + mode + " handling", messages[0])
                self.assertTrue(messages[1].startswith("WARNING: large difference"))
                self.assertEqual(result[0, 1, 0, 1], expected)
        self.assertEqual(raw[0, 0, 1, 1], 2)
        self.assertEqual(raw[1, 1, 0, 0], 8)

    def test_no_resolution_claim_for_unhandled_data(self):
        raw = np.zeros((2,) * 4)
        raw[0, 0, 1, 1], raw[1, 1, 0, 0] = 2, 8
        result, messages = self.run_handling(raw, "unhandled")
        self.assertTrue(any(m.startswith("WARNING:") for m in messages))
        self.assertFalse(any("will be resolved" in m for m in messages))
        np.testing.assert_array_equal(result, raw)

    def test_no_extra_message_without_a_discrepancy_warning(self):
        raw = np.zeros((2,) * 4)
        raw[0, 0, 1, 1], raw[1, 1, 0, 0] = 2, 2.1
        _, messages = self.run_handling(raw)
        self.assertEqual(messages, [])
        _, messages = self.run_handling(raw, warning_diff=-1)
        self.assertEqual(messages, [])

    def test_resolution_announced_once_across_value_axes(self):
        raw = np.zeros((3, 2, 2, 2, 2))
        raw[:, 0, 0, 1, 1], raw[:, 1, 1, 0, 0] = 2, 8
        _, messages = self.run_handling(raw, val_axes=1)
        self.assertEqual(sum("will be resolved" in m for m in messages), 1)
        self.assertEqual(sum(m.startswith("WARNING:") for m in messages), 3)


if __name__ == "__main__":
    unittest.main()
