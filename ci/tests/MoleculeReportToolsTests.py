"""Tool mode on a real Molecule: the template menu, `apply_template`, and grounded text against its data."""

import os as _test_os
import unittest as _test_unittest
_RUN_EXTRAS = _test_os.environ.get("MCUTILS_TEST_EXTRAS") == "1"

import os
import unittest
import warnings

from Psience.Molecools import Molecule
from McUtils.Readouts.Reports import SourceCatalog, ReportDocument, report_toolset

TEST_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "TestData")


def ok(result):
    assert result.ok, result.error
    return result.value


@_test_unittest.skipUnless(_RUN_EXTRAS, "optional dependencies; set MCUTILS_TEST_EXTRAS=1")
class MoleculeReportToolsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.filterwarnings("ignore")
        cls.mol = Molecule.from_file(os.path.join(TEST_DATA, "HOH_freq.fchk"))

    def setUp(self):
        self.catalog = SourceCatalog()
        self.catalog.add("hoh", self.mol, title="Water")
        self.doc = ReportDocument("Water", self.catalog)
        self.tools = report_toolset(self.doc, self.catalog, policy="grounded")

    def test_the_menu_tells_the_model_what_it_can_add(self):
        t = ok(self.tools.call("describe_source", {"handle": "hoh"}))["template"]
        self.assertEqual(t["name"], "molecule")
        self.assertEqual([h["heading"] for h in t["outline"]], ["Structure", "Vibrations"])
        opt = {r["name"]: r for r in t["optional"]}
        self.assertEqual(set(opt), {"source", "rdkit", "cartesians", "internals", "bonds",
                                    "normal_modes.animations", "normal_modes.displacements",
                                    "spectrum.lines", "spectrum.stick", "spectrum.levels", "spectrum.broadened"})
        self.assertTrue(all(r["description"] for r in opt.values()))
        self.assertEqual(opt["normal_modes.animations"]["part"], "animations")
        self.assertTrue(opt["cartesians"]["available"])

    def test_apply_template_with_an_optional_section(self):
        v = ok(self.tools.call("apply_template", {"handle": "hoh", "optional": ["cartesians"]}))
        inc = [b.include for b in self.doc.blocks if b.kind == "readout"]
        self.assertEqual(inc, [["overview"], ["structure"], ["normal_modes.frequencies"], ["cartesians"]])
        self.assertEqual(self.doc.validate(), [])
        self.assertEqual(v["skipped_unavailable"], [])

    def test_grounded_text_about_the_frequencies(self):
        ok(self.tools.call("apply_template", {"handle": "hoh"}))
        freq = [b.id for b in self.doc.blocks if b.kind == "readout"][2]
        ok(self.tools.call("add_text", {"text": "The bend is at 1622.3 and the stretches at 3803.3 and 3937.5.",
                                        "cites": [freq]}))
        r = self.tools.call("add_text", {"text": "The bend is at 1700.0.", "cites": [freq]})
        self.assertEqual((r.ok, r.error["code"]), (False, "ungrounded"))


if __name__ == "__main__":
    unittest.main()
