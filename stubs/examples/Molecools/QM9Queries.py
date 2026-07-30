"""Extracted from MolecoolsTests.test_QM9Queries via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest MolecoolsTests.test_QM9Queries"""

import itertools
import os.path
import pprint
from Peeves.TestUtils import *
from unittest import TestCase
from Peeves import BlockProfiler
from Psience.Molecools import Molecule, MolecularNormalModes
from Psience.Data import DipoleSurface
from McUtils.GaussianInterface import GaussianFChkReader, GaussianLogReader
from McUtils.Plots import *
import McUtils.Plots as plt
from McUtils.Coordinerds import cartesian_to_zmatrix
from McUtils.Data import UnitsData
import numpy as np, scipy
import McUtils.Numputils as nput
import McUtils.Profilers as prof
from McUtils.Formatters import TableFormatter
import McUtils.Formatters as mfmt

class MolecoolsTests(TestCase):

    def setUp(self):
        self.test_log_water = TestManager.test_data('water_OH_scan.log')
        self.test_log_freq = TestManager.test_data('water_freq.log')
        self.test_HOD = TestManager.test_data('HOD_freq.fchk')
        self.test_fchk = TestManager.test_data('water_freq.fchk')
        self.test_log_h2 = TestManager.test_data('outer_H2_scan_new.log')

    def tearDown(self):
        ...

    @classmethod
    def setup_OCHH(cls, optimize=True):
        from McUtils.Extensions import ModuleLoader
        loader = ModuleLoader(os.path.expanduser('~/Documents/Postdoc/Projects/DGB'))
        h2co_mod = loader.load('H2COPot')

        def internal_pot(coords, order=None):
            coords = coords[..., (0, 1, 3, 2, 4, 5)]
            vals = h2co_mod.InternalsPotential.get_pot(coords)
            return vals
        ochh = Molecule.from_file(TestManager.test_data('OCHH_freq.fchk'), energy_evaluator={'potential_function': internal_pot, 'distance_units': 'Angstroms', 'energy_units': 'Wavenumbers', 'strip_embedding': True}, internals=[[0, -1, -1, -1], [1, 0, -1, -1], [2, 1, 0, -1], [3, 1, 2, 0]])
        if optimize:
            base_dip = ochh.dipole_derivatives
            ochh = ochh.optimize(method='conjugate-gradient', max_iterations=50, stencil=3, prevent_oscillations=3, restart_interval=15).modify(dipole_derivatives=base_dip)
        return ochh

    @validationTest
    def test_QM9Queries(self):
        from McUtils.ExternalPrograms import QM9, RDMolecule
        import rdkit.Chem as Chem
        qm9_path = os.path.expanduser('~/Documents/Postdoc/datasets/qm9.npz')
        supplier = QM9(qm9_path)
        hmm = supplier.smiles_query('CCC', upto=500, sanitize=True)
        bond_lengths = {}
        for idx in hmm:
            data = supplier.load_data(idx, ['coords', 'atoms', 'smiles'])
            rdmol = RDMolecule.parse_smiles(data['smiles'])
            huh = Molecule(data['atoms'], coords=data['coords'] * UnitsData.convert('Angstroms', 'BohrRadius'), charge=Chem.GetFormalCharge(rdmol))
            atoms = huh.atoms
            try:
                bonds = huh.bonds
            except ValueError:
                continue
            i, j = np.array([b[:2] for b in bonds]).T
            bls = np.linalg.norm(huh.coords[i, :] - huh.coords[j, :], axis=-1) * UnitsData.bohr_to_angstroms
            for (b1, b2, t), r in zip(bonds, bls):
                a1, a2 = sorted([atoms[b1], atoms[b2]])
                bond_lengths.setdefault((a1, a2, t), []).append(r)
        import pprint
        pprint.pprint({k: np.average(v) for k, v in bond_lengths.items()})
