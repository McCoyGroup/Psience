"""The `Molecule` report template: found automatically, fixed outline, no model calls when structure-only."""

import os as _test_os
import unittest as _test_unittest
_RUN_EXTRAS = _test_os.environ.get("MCUTILS_TEST_EXTRAS") == "1"

import os
import unittest
import warnings

from Psience.Molecools import Molecule
from Psience.Molecools.Readouts.MoleculeReports import MoleculeReportTemplate
from McUtils.JobDrivers import ScriptedClient, DecisionRunner
from McUtils.Readouts.Reports import SourceCatalog, ReportHarness, WritingPolicy, find_report_template

TEST_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "TestData")


@_test_unittest.skipUnless(_RUN_EXTRAS, "optional dependencies; set MCUTILS_TEST_EXTRAS=1")
class MoleculeReportTemplateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.filterwarnings("ignore")
        cls.mol = Molecule.from_file(os.path.join(TEST_DATA, "HOH_freq.fchk"))
        cls.catalog = SourceCatalog()
        cls.catalog.add("hoh", cls.mol, title="Water")

    def test_template_is_found_from_the_readout_interface(self):
        t = find_report_template(self.catalog.get("hoh"))
        self.assertIsInstance(t, MoleculeReportTemplate)
        self.assertTrue(t.is_fixed)

    def test_every_section_has_a_clear_description(self):
        rows = self.catalog.describe("hoh")["sections"]
        for r in rows:
            d = r["description"]
            self.assertTrue(d.endswith("."), (r["name"], d))
            self.assertNotIn("`", d, r["name"])
            self.assertGreater(len(d), 25, r["name"])

    def test_structure_only_needs_no_model_calls(self):
        runner = DecisionRunner({"small": ScriptedClient([])})
        res = ReportHarness(self.catalog, runner, policy=WritingPolicy.STRUCTURE_ONLY).compose("Water", ["hoh"])
        self.assertEqual(res.cost["calls"], 0)
        md = res.doc.to_markdown()
        self.assertEqual(res.doc.validate(), [])
        self.assertLess(md.index("Structure"), md.index("Vibrations"))
        includes = [b.include for b in res.doc.blocks if getattr(b, "include", None)]
        self.assertEqual(includes, [["overview"], ["structure"], ["normal_modes.frequencies"]])
        self.assertIn("Harmonic vibrational frequencies", md)
        self.assertNotIn("animations", md)

    def test_the_vibrations_section_is_only_the_frequency_table(self):
        runner = DecisionRunner({"small": ScriptedClient([])})
        res = ReportHarness(self.catalog, runner, policy=WritingPolicy.STRUCTURE_ONLY).compose("Water", ["hoh"])
        txt = res.doc.to_markdown()
        self.assertEqual(txt.count("normal_modes.frequencies"), 1)


if __name__ == "__main__":
    unittest.main()
