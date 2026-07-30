"""Extracted from MolecoolsTests.test_FragInternalsScan via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest MolecoolsTests.test_FragInternalsScan"""

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
    def test_FragInternalsScan(self):
        import McUtils.Coordinerds as coordops
        mol = Molecule.from_string('OC=O', 'smi', confgen_opts=dict(verbose=True, random_seed=12321))
        mol = Molecule.from_string('COC=O', 'smi', confgen_opts=dict(verbose=True, random_seed=12321)).get_embedded_molecule(sel=[0, 1, 2, 3])
        coord = [(0, 4), (1, 0, 4), (2, 1, 0, 4)]
        x = 2
        zm = mol.get_bond_zmatrix(root_coordinates=coord)
        which = coordops.zmatrix_indices(zm, coord)[x]
        mol.modify(internals=zm).animate_coordinate(which, 0.2, highlight_atoms=coord[x]).show()
        return
        print(len(mol.atoms) * 3 - 6)
        specs = [(0, 1), (1, 2), (2, 3), (0, 4), (0, 1, 2), (1, 0, 4), (1, 2, 3), (4, 0, 1, 2), (0, 1, 2, 3)]
        mol.internals = {'specs': specs}
        internal_spec = coordops.InternalSpec(specs)
        ints = internal_spec.cartesians_to_internals(mol.coords)
        carts = internal_spec.internals_to_cartesians(ints, reference_cartesians=mol.coords)
        mode = -2
        geoms = mol.get_scan_coordinates([[-0.5, 0.5, 7]], which=[mode], internals='reembed')
        mol.plot(geoms, backend='x3d', highlight_atoms=specs[mode], image_size=800, include_save_buttons=True).show()
