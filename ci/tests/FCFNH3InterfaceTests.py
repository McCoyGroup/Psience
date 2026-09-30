"""Regression checks for NH3 FCF embedding and ezFCF mode output."""

from pathlib import Path
from unittest import TestCase
import xml.etree.ElementTree as ET

import numpy as np

from Psience.Vibronic import FranckCondonModel


class FCFNH3InterfaceTests(TestCase):
    @classmethod
    def setUpClass(cls):
        data = Path(__file__).parent / "TestData"
        cls.model = FranckCondonModel.from_files(
            str(data / "nh3_s0.fchk"),
            str(data / "nh3_s1.fchk"),
            logger=False,
        )

    def test_overlaps_do_not_embed_twice(self):
        states = [[0] * 6] + [[int(i == j) for i in range(6)] for j in range(6)]
        default = self.model.get_overlaps(states, return_states=False)
        prepared = self.model.get_overlaps(
            states, return_states=False, embed=False, mass_weight=False)
        np.testing.assert_allclose(default, prepared, atol=1e-12, rtol=1e-12)

    def test_ezfcf_input_has_same_mass_weighted_modes(self):
        job = self.model.get_ezFCF_input(
            2, always_run_parallel=True, print_all=False)
        root = ET.fromstring(job.format().tostring())
        for tag, nms in [
            ("initial_state", self.model.gs_nms),
            ("target_state", self.model.es_nms),
        ]:
            node = root.find(f"{tag}/normal_modes")
            self.assertEqual(node.get("if_mass_weighted"), "false")
            raw = np.fromstring(node.get("text"), sep=" ").reshape(2, 4, 9)
            modes = np.column_stack([
                raw[block, :, 3 * mode:3 * (mode + 1)].reshape(-1)
                for block in range(2) for mode in range(3)
            ])
            np.testing.assert_allclose(modes, nms.modes_by_coords, atol=5e-7)

    def test_parallel_zero_zero_against_ezfcf_v12(self):
        overlap = self.model.get_overlaps(
            [[0] * 6], return_states=False, include_rotation=False)[0]
        self.assertAlmostEqual(overlap, 0.1210281, delta=2e-6)
        self.assertAlmostEqual(overlap ** 2, 0.0146478, delta=2e-6)
