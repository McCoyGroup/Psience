"""The `DiscreteSpectrum` readout and the molecule's (non-default) ``spectrum`` section, alone and in reports.

The tests that build a `Molecule` (and so need the optional packages) run with ``MCUTILS_TEST_EXTRAS=1``."""
import os as _test_os
import unittest as _test_unittest
_RUN_EXTRAS = _test_os.environ.get("MCUTILS_TEST_EXTRAS") == "1"

import os
import tempfile
import unittest
import warnings

import numpy as np

from McUtils.Readouts import to_readout
from McUtils.Readouts.Reports import SourceCatalog, ReportDocument, ReportRenderer, report_toolset
from Psience.Molecools import Molecule
from Psience.Spectra import DiscreteSpectrum, DiscreteSpectrumReadoutInterface

TEST_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "TestData")


def ok(result):
    assert result.ok, result.error
    return result.value


def col(ro, name, section=None):
    """The exported column ``name`` (of ``section`` if several sections have one; a lone section's name is dropped)."""
    d = ro.to_flat_data()[0]
    hits = [k for k in d if k.endswith("/" + name)]
    if section is not None and len(hits) > 1:
        hits = [k for k in hits if f"/{section}/" in k]
    return np.asarray(d[hits[0]])


class DiscreteSpectrumReadoutTests(unittest.TestCase):
    def setUp(self):
        self.f = np.array([1648.2, 3832.7, 3943.9])
        self.i = np.array([60.0, 5.0, 70.0])
        self.spec = DiscreteSpectrum(self.f, self.i, kind="ir")

    def test_adapter_is_found_for_the_class_and_subclasses(self):
        from Psience.Spectra import HarmonicSpectrum
        self.assertIsInstance(to_readout(self.spec).to_text(), str)
        hs = HarmonicSpectrum(self.f, self.i)
        self.assertIn("Lines", to_readout(hs).to_text())

    def test_default_sections_are_lines_and_stick(self):
        ro = to_readout(self.spec)
        keys = set(ro.to_flat_data()[0])
        self.assertTrue(any(k.endswith("/frequency") for k in keys))
        self.assertTrue(any(k.endswith("stick/intensity") for k in keys))
        self.assertFalse(any("broadened" in k for k in keys))

    def test_values_are_the_spectrum(self):
        ro = to_readout(self.spec, include=["lines", "stick"])
        for sec in ("lines", "stick"):
            np.testing.assert_allclose(col(ro, "frequency", sec), self.f)
            np.testing.assert_allclose(col(ro, "intensity", sec), self.i)

    def test_min_intensity_and_which(self):
        ro = to_readout(self.spec, include=["lines"], lines={"min_intensity": 0.5})
        np.testing.assert_allclose(col(ro, "line"), [1, 3])
        ro = to_readout(self.spec, include=["lines"], lines={"which": [2]})
        np.testing.assert_allclose(col(ro, "frequency"), [3943.9])

    def test_broadened_is_a_long_series(self):
        n = len(col(to_readout(self.spec, include=["broadened"]), "intensity"))
        self.assertGreater(n, 1000)
        coarse = to_readout(self.spec, include=["broadened"], broadened={"step_size": 10.0})
        self.assertLess(len(col(coarse, "intensity")), n / 5)

    def test_kind_changes_captions_only(self):
        raman = DiscreteSpectrum(self.f, self.i, kind="raman")
        self.assertIn("Raman", to_readout(raman).to_text())
        self.assertNotIn("Raman", to_readout(self.spec).to_text())

    def test_normalized_spectrum_is_labelled_relative_and_unitless(self):
        text = to_readout(self.spec.normalize()).to_text()
        self.assertIn("Relative intensity", text)
        self.assertNotIn("km/mol", text)

    def test_ladder_of_a_spectrum(self):
        ro = to_readout(self.spec, include=["levels"])
        np.testing.assert_allclose(col(ro, "frequency"), self.f)
        self.assertIn("ladder", ro.to_text().lower())

    def test_empty_spectrum_has_no_broadened_curve(self):
        empty = DiscreteSpectrum(np.zeros(0), np.zeros(0))
        text = to_readout(empty, include=["broadened"]).to_text()
        self.assertIn("no lines", text)


@_test_unittest.skipUnless(_RUN_EXTRAS, "optional dependencies; set MCUTILS_TEST_EXTRAS=1")
class MoleculeSpectrumSectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.filterwarnings("ignore")
        cls.mol = Molecule.from_file(os.path.join(TEST_DATA, "HOH_freq.fchk"))
        cls.spec = cls.mol.get_harmonic_spectrum()

    def test_not_a_default_section(self):
        rows = {r["name"]: r for r in self.mol.get_readout_interface().list_readout_sections()}
        self.assertIn("spectrum", rows)
        self.assertFalse(rows["spectrum"]["default"])
        self.assertNotIn("Harmonic spectrum", self.mol.to_readout().to_text())

    def test_included_on_request_with_the_molecules_spectrum(self):
        ro = self.mol.to_readout(include=["spectrum.lines"])
        np.testing.assert_allclose(col(ro, "frequency"), self.spec.frequencies, rtol=1e-6)
        np.testing.assert_allclose(col(ro, "intensity"), self.spec.intensities, rtol=1e-6)

    def test_whole_section_is_lines_and_stick(self):
        keys = list(self.mol.to_readout(include=["spectrum"]).to_flat_data()[0])
        self.assertTrue(any("spectrum/lines" in k for k in keys))
        self.assertTrue(any("spectrum/stick" in k for k in keys))
        self.assertFalse(any("spectrum/broadened" in k for k in keys))

    def test_options_reach_the_sections(self):
        ro = self.mol.to_readout(include=["spectrum.lines"], spectrum={"lines": {"min_intensity": 0.5}})
        np.testing.assert_allclose(col(ro, "line"), [1, 3])
        coarse = self.mol.to_readout(include=["spectrum.broadened"], spectrum={"broadened": {"step_size": 5.0}})
        fine = self.mol.to_readout(include=["spectrum.broadened"])
        self.assertLess(len(col(coarse, "intensity")), len(col(fine, "intensity")))

    def test_unknown_kind_is_an_error(self):
        with self.assertRaises(Exception):
            self.mol.to_readout(include=["spectrum"], spectrum={"kind": "uv"}, strict=True)

    def test_raman_intensities_use_the_spectrums_own_unit_and_label(self):
        text = self.mol.to_readout(include=["spectrum.lines"], spectrum={"kind": "raman"}).to_text()
        self.assertIn("Raman lines", text)
        self.assertIn("Raman activity", text)
        self.assertNotIn("km/mol", text)
        self.assertIn("km/mol", self.mol.to_readout(include=["spectrum.lines"]).to_text())

    def test_raman_cross_sections_are_written_in_scientific_notation(self):
        ro = self.mol.to_readout(include=["spectrum.lines"],
                                 spectrum={"kind": "raman", "spectrum_options": {"pump_frequency": 0.0656}})
        text = ro.to_text()
        self.assertIn("Raman differential cross section", text)
        self.assertRegex(text, r"\de-\d\d")

    def test_lines_match_the_normal_mode_table(self):
        ro = self.mol.to_readout(include=["spectrum.lines", "normal_modes.frequencies"])
        d = ro.to_flat_data()[0]
        sp = np.asarray(next(v for k, v in d.items() if "spectrum" in k and k.endswith("/frequency")))
        nm = np.asarray(next(v for k, v in d.items() if "normal_modes" in k and k.endswith("/frequency")))
        np.testing.assert_allclose(sp, nm, rtol=1e-5)


@_test_unittest.skipUnless(_RUN_EXTRAS, "optional dependencies; set MCUTILS_TEST_EXTRAS=1")
@_test_unittest.skipUnless(_RUN_EXTRAS, "optional dependencies; set MCUTILS_TEST_EXTRAS=1")
class MoleculeLadderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.filterwarnings("ignore")
        cls.mol = Molecule.from_file(os.path.join(TEST_DATA, "HOH_freq.fchk"))
        cls.spec = cls.mol.get_harmonic_spectrum()

    def test_ladder_is_the_frequencies(self):
        ro = self.mol.to_readout(include=["spectrum.levels"])
        np.testing.assert_allclose(col(ro, "frequency"), self.spec.frequencies, rtol=1e-6)
        np.testing.assert_allclose(col(ro, "level"), [1, 2, 3])
        self.assertIn("<svg", ro.to_html_fragment().html)

    def test_ladder_is_not_in_the_default_section(self):
        keys = list(self.mol.to_readout(include=["spectrum"]).to_flat_data()[0])
        self.assertFalse(any("level" in k for k in keys))

    def test_ladder_options_filter_the_levels(self):
        ro = self.mol.to_readout(include=["spectrum.levels"], spectrum={"levels": {"min_intensity": 0.5}})
        np.testing.assert_allclose(col(ro, "level"), [1, 2])

    def test_a_ladder_block_is_summarised_for_a_report(self):
        cat = SourceCatalog()
        cat.add("hoh", self.mol)
        doc = ReportDocument("t", cat)
        tools = report_toolset(doc, cat)
        r = ok(tools.call("add_readout", {"handle": "hoh", "include": ["spectrum.levels"]}))["block_id"]
        ok(tools.call("add_text", {"text": "The top level is at {%s.frequency[max]}." % r}))
        self.assertIn(f"{self.spec.frequencies.max():.1f}", doc.get(doc.blocks[-1].id).text)


class MoleculeSpectrumReportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.filterwarnings("ignore")
        cls.mol = Molecule.from_file(os.path.join(TEST_DATA, "HOH_freq.fchk"))
        cls.cat = SourceCatalog()
        cls.cat.add("hoh", cls.mol, title="Water")
        cls.spec = cls.mol.get_harmonic_spectrum()

    def tools(self, **kw):
        doc = ReportDocument("Water spectrum", self.cat)
        return doc, report_toolset(doc, self.cat, **kw)

    def test_the_catalog_lists_spectrum_sections_as_non_default_with_kinds(self):
        rows = {r["name"]: r for r in self.cat.describe("hoh")["sections"]}
        self.assertIn("spectrum", rows)
        self.assertFalse(rows["spectrum"]["default"])

    def test_strongest_line_by_selector(self):
        doc, tools = self.tools()
        r = ok(tools.call("add_readout", {"handle": "hoh", "include": ["spectrum.stick"]}))["block_id"]
        ok(tools.call("add_text", {"text": f"The strongest line is at {{{r}.frequency[max:intensity]}} and reaches "
                                           f"{{{r}.intensity[max]}}."}))
        text = doc.get(doc.blocks[-1].id).text
        top = int(np.argmax(self.spec.intensities))
        self.assertIn(f"{self.spec.frequencies[top]:.1f}", text)
        self.assertIn(f"{self.spec.intensities[top]:.4f}", text)

    def test_stick_block_peek_lists_placeholders(self):
        doc, tools = self.tools()
        out = ok(tools.call("peek", {"handle": "hoh", "section": "spectrum"}))
        self.assertIn("frequency", str(out))

    def test_broadened_block_is_summarised_not_listed(self):
        doc, tools = self.tools()
        rec = ok(tools.call("add_readout", {"handle": "hoh", "include": ["spectrum.broadened"]}))
        joined = "\n".join(rec["placeholders"])
        self.assertIn("[max]", joined)
        self.assertLess(len(joined), 1500)

    def test_optional_template_sections_are_accepted(self):
        from Psience.Molecools.Readouts.MoleculeReports import MoleculeReportTemplate
        self.assertIn("spectrum.stick", MoleculeReportTemplate.optional)
        doc, tools = self.tools()
        out = ok(tools.call("apply_template", {"handle": "hoh", "optional": ["spectrum.stick"]}))
        self.assertTrue(out["added"], out)

    def test_rendered_html_and_pptx_contain_the_figure(self):
        doc, tools = self.tools()
        ok(tools.call("add_heading", {"text": "Spectrum", "level": 2}))
        ok(tools.call("add_readout", {"handle": "hoh", "include": ["spectrum.lines", "spectrum.stick"]}))
        with tempfile.TemporaryDirectory() as d:
            try:
                res = ReportRenderer(self.cat, strict=False).render(doc, d, formats=("html", "pptx"))
            except ImportError:
                self.skipTest("an optional package for the HTML report is missing")
            self.assertIn("pptx", res.paths, res.problems)
            if "html" in res.paths:
                html = open(res.paths["html"], encoding="utf-8").read()
                self.assertIn("Intensity", html)
                self.assertIn("<svg", html)


if __name__ == "__main__":
    unittest.main()
