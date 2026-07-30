"""Extracted from MolecoolsTests.test_PointGroups via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest MolecoolsTests.test_PointGroups"""

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
    def test_PointGroups(self):
        print()
        mol = Molecule(['C', 'C', 'H', 'H', 'H', 'H', 'H', 'H'], [[1.3054152479869834, 0.24647130925656535, -0.5530793729770965], [-1.3054057993563604, -0.24642973528182355, 0.553026460645607], [2.7502129154741146, -0.7828587313748425, 0.5769125988608768], [1.3606991857629105, -0.4278075384495296, -2.5445918158535847], [1.7132634991082158, 2.309094146202721, -0.4999931866841094], [-2.751762490896308, 0.7821255176384876, -0.5756408131790036], [-1.3610034316689752, 0.4288638953531954, 2.544251665151152], [-1.7114191264105811, -2.3094531941664, 0.49911446403615845]])
        mol = Molecule(['C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H'], [[1.5314834, -2.07901109, 0.03362855], [2.66130207, 0.27021048, -0.00918031], [1.09977014, 2.36405638, -0.04301921], [-1.49251245, 2.07781621, -0.03367612], [-2.66353391, -0.24015647, 0.00865698], [-1.08442635, -2.34605922, 0.04267503], [2.65427297, -3.82363066, 0.06211177], [4.65827828, 0.4512349, -0.01567418], [1.93012887, 4.24957282, -0.07725653], [-2.65023509, 3.80315125, -0.06177008], [-4.68709859, -0.49931772, 0.01657436], [-1.95742935, -4.22786689, 0.07692975]]).get_embedded_molecule(load_properties=False)
        (coords, atoms), pg = mol.symmetrize(grouping_tol=0.8, tol=0.3, return_point_group=True)
        self.assertEquals(coords.shape, mol.coords.shape)
        base = mol.plot(backend='x3d', include_save_buttons=True)
        (coords, atoms), pg2 = Molecule(atoms, coords).symmetrize(pg, return_point_group=True)
        Molecule(atoms, coords).plot(figure=base, highlight_atoms=True)
        pg2.plot(figure=base, origin=mol.center_of_mass * UnitsData.bohr_to_angstroms)
        return
        mol = Molecule(['C', 'H', 'H', 'H', 'H'], [[0.000149273393, -8.31939792e-05, -3.01752454e-05], [-0.129672979, 2.0272439, 0.364678828], [-0.990194225, -0.448290077, -1.7545249], [-0.863649504, -1.03676833, 1.5617178], [1.98411318, -0.54251937, -0.171993834]])
        pg = mol.get_point_group(grouping_tol=0.8, tol=0.3, mom_tol=5, verbose=False)
        print(pg)
        print(pg.elements)
        base = mol.plot(backend='x3d', principle_axes=True)
        pg.plot(figure=base, origin=mol.center_of_mass * UnitsData.bohr_to_angstroms)
        base.show()
