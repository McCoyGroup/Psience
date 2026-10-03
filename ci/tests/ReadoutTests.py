"""Molecule readouts (McUtils.Jupyter.Readouts) with normal modes, exported to HTML and PowerPoint."""
import io
import json
import os
import re
import tempfile
import unittest
import zipfile

import numpy as np

from Psience.Molecools import Molecule
from McUtils.Data import UnitsData

TEST_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "TestData")


class MoleculeReadoutTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mol = Molecule.from_file(os.path.join(TEST_DATA, "HOH_freq.fchk"))
        cls.readout = cls.mol.to_readout(include="all")

    def test_sections(self):
        ids = [c.id for c in self.readout.children]
        self.assertEqual(ids, ["overview", "source", "rdkit", "structure", "cartesians", "internals",
                               "bonds", "normal_modes"])
        self.assertEqual([c.id for c in self.readout["normal_modes"].children],
                         ["frequencies", "animations", "displacements"])
        self.assertEqual(self.readout["overview/overview"].data["charge"].value, 0)
        self.assertEqual(self.readout["source/source"].data["format"].value, "fchk")

    def test_frequencies_match_across_exports(self):
        expected = np.asarray(self.mol.normal_modes.modes.freqs) * UnitsData.convert("Hartrees", "Wavenumbers")
        flat, meta = self.readout.to_flat_data()
        np.testing.assert_allclose(flat["molecule/normal_modes/frequencies/frequency"], expected)
        self.assertEqual(meta["molecule/normal_modes/frequencies/frequency"]["unit"], "Wavenumbers")
        html = self.readout.to_html()
        for f in expected:
            self.assertIn(f"{f:.1f}", html)

    def test_nested_options_and_dotted_includes(self):
        ro = self.mol.to_readout(include=["normal_modes.frequencies"])
        self.assertEqual([c.id for c in ro["normal_modes"].children], ["frequencies"])
        ro = self.mol.to_readout(include=["normal_modes"], normal_modes={"animations": {"which": [0, 2]}})
        self.assertEqual([c.id for c in ro["normal_modes/animations/gallery"].children], ["mode_1", "mode_3"])

    def test_html_and_powerpoint(self):
        ro = self.mol.to_readout(include=["overview", "structure", "normal_modes"])
        with tempfile.TemporaryDirectory() as tmp:
            html_file = ro.to_html(os.path.join(tmp, "water.html"))
            with open(html_file) as f:
                html = f.read()
            self.assertEqual(html.count("<X3D"), 4)  # structure + 3 modes
            self.assertIn("x3dom-full.js", html)
            pptx = ro.to_pptx(os.path.join(tmp, "water.pptx"))
            with zipfile.ZipFile(pptx) as z:
                models = [n for n in z.namelist() if n.endswith(".glb")]
                self.assertEqual(len(models), 4)
                slides = "".join(z.read(n).decode() for n in z.namelist() if re.match(r"ppt/slides/slide\d+\.xml$", n))
                self.assertEqual(slides.count("<am3d:model3d"), 4)
                for name in models:
                    data = z.read(name)
                    length = int.from_bytes(data[12:16], "little")
                    gltf = json.loads(data[20:20 + length])
                    if gltf.get("animations"):
                        self.assertTrue(gltf.get("skins"))   # PowerPoint only plays skinned animation

    def test_data_exports(self):
        with tempfile.TemporaryDirectory() as tmp:
            f = os.path.join(tmp, "water.npz")
            self.readout.to_npz(f)
            data = np.load(f, allow_pickle=False)
            self.assertEqual(data["molecule/normal_modes/displacements"].shape, (3, 3, 3))
            self.assertEqual(str(data["molecule/rdkit/inchi_key"]), "XLYOFNOQVPJJNP-UHFFFAOYSA-N")
        frames = self.readout.to_pandas()
        self.assertEqual(list(frames["molecule/cartesians"].columns), ["index", "atom", "mass", "x", "y", "z"])


if __name__ == "__main__":
    unittest.main()
