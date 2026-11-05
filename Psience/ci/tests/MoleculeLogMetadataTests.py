"""Molecule imports log electronic states through the McUtils parser stack."""
import os
import tempfile
import unittest

import numpy as np

from McUtils.Data import UnitsData
from Psience.Molecools import Molecule


# Small format fixtures matching the coordinate and metadata layouts in the
# repository's water_freq.log, qchem_samp.out and proplybenz.out examples.
GAUSSIAN_COORDS = """ Input orientation:
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
      1          1           0        0.000000    0.000000    0.000000
      2          1           0        0.000000    0.000000    0.740000
 ---------------------------------------------------------------------
"""
QCHEM_COORDS = """ Standard Nuclear Orientation (Angstroms)
    I     Atom           X                Y                Z
 ----------------------------------------------------------------
    1      H       0.0000000000     0.0000000000     0.0000000000
    2      H       0.0000000000     0.0000000000     0.7400000000
 ----------------------------------------------------------------

"""
ORCA_COORDS = """CARTESIAN COORDINATES (A.U.)
----------------------------
  NO LB      ZA    FRAG     MASS         X           Y           Z
   0 H     1.0000    0     1.008    0.000000    0.000000    0.000000
   1 H     1.0000    0     1.008    0.000000    0.000000    1.400000

"""


def fixture(mode, charge=1, multiplicity=2):
    if mode == 'log':
        metadata = '' if charge is None else f' Charge = {charge:+d} Multiplicity = {multiplicity}\n NAtoms= 2\n'
        return metadata + GAUSSIAN_COORDS
    if mode == 'qchem':
        metadata = '' if charge is None else f'$molecule\n{charge:+d} {multiplicity}\nH 0 0 0\nH 0 0 .74\n$end\n'
        return metadata + QCHEM_COORDS + 'Thank you very much for using Q-Chem\n'
    metadata = '' if charge is None else (
        f' Total Charge           Charge          ....    {charge}\n'
        f' Multiplicity           Mult            ....    {multiplicity}\n\n'
    )
    # ORCA prints the computed state after its initial geometry.
    return ORCA_COORDS + metadata + 'ORCA TERMINATED NORMALLY\n'


class MoleculeLogMetadataTests(unittest.TestCase):
    modes = ['log', 'qchem', 'orca']

    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tempdir.cleanup)

    def load(self, mode, fixture_charge=1, fixture_multiplicity=2, text=None, **opts):
        path = os.path.join(self.tempdir.name, 'job.' + ('out' if mode != 'log' else 'log'))
        with open(path, 'w') as stream:
            stream.write(fixture(mode, fixture_charge, fixture_multiplicity) if text is None else text)
        return Molecule.from_file(path, mode=mode, energy=0, **opts)

    def test_cation_doublet_in_molecule_and_overview(self):
        for mode in self.modes:
            with self.subTest(mode=mode):
                mol = self.load(mode)
                self.assertEqual((mol.charge, mol.spin), (1, 2))
                overview = mol.to_readout(include=['overview'])['overview/overview'].data
                self.assertEqual(overview['charge'].value, 1)
                self.assertEqual(overview['spin'].value, 2)

    def test_neutral_anion_and_triplet(self):
        for mode in self.modes:
            for charge, multiplicity in [(0, 1), (-1, 2), (0, 3)]:
                with self.subTest(mode=mode, charge=charge, multiplicity=multiplicity):
                    mol = self.load(mode, charge, multiplicity)
                    self.assertEqual((mol.charge, mol.spin), (charge, multiplicity))

    def test_independent_explicit_overrides(self):
        for mode in self.modes:
            for opts, expected in [({'charge': 0}, (0, 2)), ({'spin': 1}, (1, 1)),
                                   ({'charge': -1, 'spin': 3}, (-1, 3)),
                                   ({'charge': None, 'spin': None}, (1, 2))]:
                with self.subTest(mode=mode, opts=opts):
                    mol = self.load(mode, **{'text': fixture(mode), **opts})
                    self.assertEqual((mol.charge, mol.spin), expected)

    def test_format_options_precedence(self):
        for mode in self.modes:
            with self.subTest(mode=mode):
                mol = self.load(mode, text=fixture(mode), format_options={mode: {'charge': -1, 'spin': 3}}, spin=1)
                self.assertEqual((mol.charge, mol.spin), (-1, 1))

    def test_missing_metadata_preserves_defaults_and_overrides(self):
        for mode in self.modes:
            with self.subTest(mode=mode):
                text = fixture(mode, None)
                mol = self.load(mode, text=text)
                self.assertEqual((mol.charge, mol.spin), (0, None))
                mol = self.load(mode, text=text, **{'charge': -1, 'spin': 2})
                self.assertEqual((mol.charge, mol.spin), (-1, 2))

    def test_geometry_units_masses_and_source_are_preserved(self):
        for mode in self.modes:
            with self.subTest(mode=mode):
                mol = self.load(mode)
                z = 1.4 if mode == 'orca' else .74 * UnitsData.convert('Angstroms', 'BohrRadius')
                self.assertEqual(mol.atoms, ('H', 'H'))
                np.testing.assert_allclose(mol.coords, [[0, 0, 0], [0, 0, z]])
                self.assertTrue(mol.source_file.endswith('job.' + ('log' if mode == 'log' else 'out')))
                if mode == 'orca':
                    np.testing.assert_allclose(mol.masses, [1.008, 1.008])

    def test_qchem_coordinate_filters_keep_metadata(self):
        text = fixture('qchem') + QCHEM_COORDS.replace('0.7400000000', '0.8400000000')
        for coordinate_filter in [None, lambda coords: coords[:1], lambda coords: coords[0]]:
            with self.subTest(coordinate_filter=coordinate_filter):
                result = self.load('qchem', text=text, coordinate_filter=coordinate_filter)
                molecules = result if isinstance(result, list) else [result]
                self.assertGreater(len(molecules), 0)
                self.assertTrue(all((mol.charge, mol.spin) == (1, 2) for mol in molecules))

    def test_gaussian_num_keeps_metadata(self):
        text = fixture('log') + GAUSSIAN_COORDS.replace('0.740000', '0.840000')
        mol = self.load('log', text=text, num=1)
        self.assertEqual((mol.charge, mol.spin), (1, 2))
        self.assertAlmostEqual(mol.coords[-1, -1], .74 * UnitsData.convert('Angstroms', 'BohrRadius'))

    def test_gaussian_linked_job_metadata_follows_selected_geometry(self):
        route = ' --------------------\n #p uhf/sto-3g\n --------------------\n'
        text = route + fixture('log', 0, 1) + ' Normal termination of Gaussian\n'
        text += route + fixture('log', 1, 2).replace('0.740000', '0.840000')
        first = self.load('log', text=text, num=1)
        last = self.load('log', text=text)
        self.assertEqual((first.charge, first.spin), (0, 1))
        self.assertEqual((last.charge, last.spin), (1, 2))
        self.assertAlmostEqual(last.coords[-1, -1], .84 * UnitsData.convert('Angstroms', 'BohrRadius'))
        # An unfinished later job with no geometry must not supply the state.
        text += ' Normal termination of Gaussian\n' + route + ' Charge = -1 Multiplicity = 2\n'
        last = self.load('log', text=text)
        self.assertEqual((last.charge, last.spin), (1, 2))


if __name__ == '__main__':
    unittest.main()
