"""Extracted from MolecoolsTests.test_FragEmbedding via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest MolecoolsTests.test_FragEmbedding"""

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
    def test_FragEmbedding(self):
        bits = Molecule.from_string('FC1CCC(C(=C)C)CC1.C1OCC(C(OC)O)C1C(OC)O', 'smi', confgen_opts=dict(verbose=True, random_seed=12321))
        s3 = bits.modify(internals={'specs': [{'transrot': bits.fragment_indices[1]}]})
        s3.calculate_energy(order=1)
        return
        bits.to_file('/Users/Mark/Desktop/struct.mol')
        ugh = bits.modify(internals=bits.get_canonical_zmatrix()).to_string('zmat')
        new_bits = Molecule.from_string(ugh)
        print(np.linalg.norm(bits.fragments[0].center_of_mass - bits.fragments[1].center_of_mass))
        print(np.linalg.norm(new_bits.fragments[0].center_of_mass - new_bits.fragments[1].center_of_mass))
        enc_smi = bits.to_string('smi', remove_hydrogens=True, include_tag=True)
        print(enc_smi)
        uuuh = Molecule.from_string(enc_smi, 'smi', add_implicit_hydrogens=True, confgen_opts=dict(verbose=True, ignore_smoothing_failures=True))
        print(enc_smi)
        print(np.linalg.norm(bits.fragments[0].center_of_mass - bits.fragments[1].center_of_mass))
        print(np.linalg.norm(uuuh.fragments[0].center_of_mass - uuuh.fragments[1].center_of_mass))
