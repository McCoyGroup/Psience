"""Extracted from MolecoolsTests.test_BackboneChains via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest MolecoolsTests.test_BackboneChains"""

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
    def test_BackboneChains(self):
        from Psience.Molecools import Molecule
        import McUtils.Coordinerds as coordops
        from Psience.Reactions import Reaction
        woof = Reaction.from_smiles('C=C.C=CC=C(c1c2ccccc2ccc1)>>C1CCC=CC1(c1c2ccccc2ccc1)', fragment_expansion_method='centroid', optimize=True, min_distance=0.1, add_radius=False, expansion_factor=0.01)
        reactant_complex = woof.reactant_complex
        full_zmat = reactant_complex.get_bond_zmatrix()
        int_comp = reactant_complex.modify(internals=full_zmat)
        return
        woof = Molecule.construct('CCCC')
        zm = coordops.chain_zmatrix(4)
        print(woof.atoms)
        print(coordops.add_missing_zmatrix_bonds(zm, [b[:2] for b in woof.bonds]))
        return
        napthalene = Molecule.construct('CCCCC(c1c2ccccc2ccc1)CCCC')
        chains = napthalene.edge_graph.segment_by_chains()
        zm = coordops.bond_graph_zmatrix([b[:2] for b in napthalene.bonds], chains)
        print(zm)
        return
        backbone, (side_chain,) = napthalene.edge_graph.segment_by_chains()
        atom_styles = {backbone[0]: {'color': 'white', 'glow': 'red'}, backbone[-1]: {'color': 'white', 'glow': 'blue'}}
        for a in side_chain:
            atom_styles[a] = {'color': 'white', 'glow': 'purple'}
        bond_style = {k: {'color': 'white', 'glow': 'red'} for i in range(len(backbone) - 1) for k in [(backbone[i], backbone[i + 1]), (backbone[i + 1], backbone[i])]}
        for i in range(len(side_chain) - 1):
            bond_style[side_chain[i], side_chain[i + 1]] = {'color': 'white', 'glow': 'purple'}
            bond_style[side_chain[i + 1], side_chain[i]] = {'color': 'white', 'glow': 'purple'}
        napthalene.plot(highlight_atoms=backbone[1:-1], atom_style=atom_styles, bond_style=bond_style, include_save_buttons=True).show()
