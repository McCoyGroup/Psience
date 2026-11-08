"""A report built by hand from a Molecule: Markdown, HTML (with an X3DOM scene) and PowerPoint."""

import os as _test_os
import unittest as _test_unittest
_RUN_EXTRAS = _test_os.environ.get("MCUTILS_TEST_EXTRAS") == "1"

import os
import re
import tempfile
import unittest
import warnings
import zipfile

import numpy as np

from Psience.Molecools import Molecule
from McUtils.Data import UnitsData
from McUtils.Readouts.Reports import SourceCatalog, ReportDocument, ReportRenderer

TEST_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "TestData")


def slide_texts(path):
    with zipfile.ZipFile(path) as z:
        names = sorted((n for n in z.namelist() if re.match(r"ppt/slides/slide\d+\.xml$", n)),
                       key=lambda n: int(re.findall(r"\d+", n.rsplit("/", 1)[1])[0]))
        return [re.findall(r"<a:t>([^<]*)</a:t>", z.read(n).decode("utf-8")) for n in names]


@_test_unittest.skipUnless(_RUN_EXTRAS, "optional dependencies; set MCUTILS_TEST_EXTRAS=1")
class MoleculeReportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.filterwarnings("ignore")
        cls.mol = Molecule.from_file(os.path.join(TEST_DATA, "HOH_freq.fchk"))
        cls.catalog = SourceCatalog()
        cls.catalog.add("hoh", cls.mol, title="Water")
        cls.freqs = np.asarray(cls.mol.normal_modes.modes.freqs) * UnitsData.convert("Hartrees", "Wavenumbers")

    def build(self):
        doc = ReportDocument("Water, harmonic analysis", self.catalog)
        doc.add_heading("Structure")
        doc.add_readout("hoh", include=["overview"], caption="Formula, size and charge of the molecule.")
        doc.add_readout("hoh", include=["structure"], caption="The equilibrium geometry.")
        doc.add_heading("Vibrations")
        doc.add_readout("hoh", include=["normal_modes.frequencies"], caption="Harmonic frequencies in wavenumbers.")
        return doc

    def test_the_catalog_describes_the_molecule(self):
        d = self.catalog.describe("hoh")
        rows = {r["name"]: r for r in d["sections"]}
        self.assertTrue({"overview", "structure", "normal_modes"} <= set(rows))
        self.assertEqual(rows["structure"]["kind"], "figure")
        self.assertEqual(rows["overview"]["kind"], "fields")
        self.assertTrue(all(r["description"] for r in rows.values()))
        self.assertTrue(d["type"].endswith("Molecule"))
        self.assertIn("frequenc", self.catalog.peek("hoh", "normal_modes", max_chars=3000).lower())

    def test_markdown_round_trip_with_a_molecule(self):
        doc = self.build()
        md = doc.to_markdown()
        again = ReportDocument.from_markdown(md, self.catalog)
        self.assertEqual(again.to_markdown(), md)
        self.assertEqual(again.validate(), [])
        self.assertEqual(again.get("r3").include, ["normal_modes.frequencies"])

    def test_html_and_pptx(self):
        doc = self.build()
        with tempfile.TemporaryDirectory() as tmp:
            res = ReportRenderer(self.catalog, timeout=300).render(doc, tmp)
            self.assertTrue(res.ok, res.problems)
            self.assertEqual(res.statuses, {"r1": "ok", "r2": "ok", "r3": "ok"})
            with open(res.paths["html"], encoding="utf-8") as f:
                html = f.read()
            for f_ in self.freqs:                       # the frequencies, in the display units
                self.assertIn(f"{f_:.1f}", html)
            head = html[:html.index("</head>")]
            self.assertIn("x3dom", head.lower())        # the scene block brought the runtime into <head>
            self.assertIn("<X3D", html)
            self.assertIn("Harmonic frequencies in wavenumbers.", html)
            self.assertEqual(html.count('class="mcr-root"'), 3)
            self.assertNotIn("report_source(", html)                  # no code, and no repeated molecule title blocks
            self.assertNotIn('class="mcr-header"', html.split("<body")[1])
            slides = slide_texts(res.paths["pptx"])
            self.assertGreaterEqual(len(slides), 2)
            flat = " ".join(t for s in slides for t in s)
            for heading in ("Structure", "Vibrations"):
                self.assertIn(heading, flat)
            self.assertIn("Harmonic frequencies in wavenumbers.", flat)
            for f_ in self.freqs:
                self.assertIn(f"{f_:.1f}", flat)
            with zipfile.ZipFile(res.paths["pptx"]) as z:
                self.assertTrue(any(n.endswith(".glb") for n in z.namelist()))   # the structure travelled as a model

    def test_a_second_render_runs_nothing(self):
        doc = self.build()
        with tempfile.TemporaryDirectory() as tmp:
            r = ReportRenderer(self.catalog, timeout=300)
            self.assertEqual(len(r.render(doc, tmp, formats=("html",)).executed), 3)
            self.assertEqual(r.render(doc, tmp, formats=("html",)).executed, [])

    def test_full_normal_modes_section_renders_as_markup(self):
        """The mode animations carry playback controls (void elements such as <input>); the page must keep them as markup."""
        doc = ReportDocument("Modes", self.catalog)
        doc.add_readout("hoh", include=["normal_modes"], caption="All the modes.")
        with tempfile.TemporaryDirectory() as d:
            out = ReportRenderer(self.catalog, timeout=300).render(doc, os.path.join(d, "r"), formats=("html",))
            self.assertTrue(out.ok, out.problems)
            html = open(out.paths["html"], encoding="utf-8").read()
        body = html.split("<body", 1)[1]
        self.assertNotIn("&lt;", body)
        self.assertIn('data-role="playback-position"', body)
        self.assertIn("<h1>Modes</h1>", body)



if __name__ == "__main__":
    unittest.main()
