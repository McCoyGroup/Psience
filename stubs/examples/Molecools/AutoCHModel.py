"""Extracted from MolecoolsTests.test_AutoCHModel via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest MolecoolsTests.test_AutoCHModel"""

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
    def test_AutoCHModel(self):
        import McUtils.Coordinerds as coordops
        from Psience.BasisReps import LocalHarmonicModel, StateMaker, TaborCHModel
        propylbenzene = Molecule.from_file(TestManager.test_data('proplybenz.hess'))
        lhm = LocalHarmonicModel
        model = LocalHarmonicModel.from_molecule(propylbenzene, oblique=False, coordinate_filter=lambda coords: {c: l for c, l in coords.items() if l.atoms in {'CH', 'HCH'}}, anharmonic_couplings={lhm.state_pair(((2, 1),), ('CH', 'stretch'), (('HCH', 'bend'), ('HCH', 'bend'))): 5.6 / UnitsData.hartrees_to_wavenumbers}, anharmonic_shifts={lhm.state('benzene', None, 'CH', 'stretch'): -30 / UnitsData.hartrees_to_wavenumbers, lhm.state('methyl', 'CH', 'stretch'): -8 / UnitsData.hartrees_to_wavenumbers, lhm.state('ethyl', 'CH', 'stretch'): -5 / UnitsData.hartrees_to_wavenumbers, lhm.state('HCH', 'bend'): -2 * 9.6 / UnitsData.hartrees_to_wavenumbers})
        dim = model.basis.ndim
        state = StateMaker(dim, mode='high-low')
        wfns = model.get_wavefunctions({'max_freq': 3250 / UnitsData.hartrees_to_wavenumbers, 'min_freq': 3050 / UnitsData.hartrees_to_wavenumbers, 'max_quanta': 3})
        print(model.internals)
        print(mfmt.TableFormatter('{:.3f}').format(wfns.hamiltonian * UnitsData.hartrees_to_wavenumbers))
        spec = wfns.get_spectrum()
        print(spec.intensities)
        return
        spec.plot().show()
        ham = model.get_hamiltonian([state([dim - 3, 2]), state(dim - 4, dim - 3)])
        print()
        print(TableFormatter('{:.0f}').format(ham * 219474.63))
        return
        model.get_hamiltonian([state(1)])
        stretch, angles, dihedrals = coordops.get_stretch_coordinate_system([tuple(s[:2]) for s in pb.bonds])
        labels = pb.edge_graph.get_label_types()
        stretch_types = [coordops.get_coordinate_label(c, labels) for c in stretch]
        bend_types = [coordops.get_coordinate_label(c, labels) for c in angles]
        good_coords = {c: l for c, l in zip(stretch, stretch_types) if l.atoms == 'CH'}
        good_coords.update({c: l for c, l in zip(angles, bend_types) if l.atoms == 'HCH'})
        nms = pb.get_normal_modes()
        loc_modes = nms.localize(internals=good_coords).make_oblique()
        base_hess = loc_modes.compute_hessian()
        print('Base Oblique Hessian')
        print(TableFormatter('{:.0f}').format(base_hess * 219474.63))
        print('Scaled Oblique Hessian')
        scaled_hess = modify_internal_hamiltonian(base_hess, good_coords, scaling_types={('CH', 'stretch'): 0.96, (('methyl', 'CH', 'stretch'), ('methyl', 'CH', 'stretch')): 0.9, (('ethyl', 'CH', 'stretch'), ('ethyl', 'CH', 'stretch')): 0.96})
        print(TableFormatter('{:.0f}').format(scaled_hess * 219474.63))
