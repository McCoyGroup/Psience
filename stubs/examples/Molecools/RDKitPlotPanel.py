"""Extracted from MolecoolsTests.test_RDKitPlotPanel via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest MolecoolsTests.test_RDKitPlotPanel"""

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
    def test_RDKitPlotPanel(self):
        systems = ['CON=[C:3]([C:2](=O)[NH:1]C)/[C:4]1=[CH:5]/[C:6]2=[N+:7](C(C)C(C)=NN(c3ccccc3)c3ccccc3)OC1=C(C)C2', 'CO/N=C(\\C(O)=N[C@@H]1C(=O)N2C(C(=O)O)=C(C[N+:7]3=[C:6]=[C:5]4CC=C/[C:3]([S:2](=O)(=O)[OH:1])=[C:4]\\4C=C3)CS[C@H]12)c1csc(=N)[nH]1', 'C=C/C(=C\\N)C/C=C(\\C=N)c1ccc2c(ccc3c4ccc(C5=C[CH:2]([C:3]6=[N+:7]=[CH:6]/[CH:5]=[CH:4]/6)[NH:1]C(/C(C=C)=C/N)=C5)cc4n(-c4ccccc4)c23)c1', 'C=C/C(=C\\N)C/C=C(\\C=N)c1ccc2c(ccc3c4ccc(C5=C[CH:2]([C:3]6=[N+:7]=[CH:6]/[CH:5]=[CH:4]/6)[NH:1]C(/C(C=C)=C/N)=C5)cc4n(-c4ccccc4)c23)c1', 'CCCCCOc1cc(C2C(=O)C(=[C:3]3[C:2]([NH:1]S(C)(=O)=O)=C[C:6](=[N+:7](C)Cc4ccccc4)/[C:5](C)=[CH:4]/3)C2[O-])c(NS(C)(=O)=O)cc1N(C)Cc1ccccc1', 'C=C/C(CNC(=O)c1ccnc(NC[C:6]2=[N+:7](C)[C:3]([C:2]3=CCN=C[NH:1]3)/[N:4]=[N:5]/2)c1)=C(\\N=C)C(F)(F)F', 'C=C/C(CNC(=O)c1ccnc(NC[C:6]2=[N+:7](C)[C:3]([C:2]3=CCN=C[NH:1]3)/[N:4]=[N:5]/2)c1)=C(\\N=C)C(F)(F)F', 'CCCCCCCCCCCCC/C=C/C(O)C(CO)NC(=O)Cn1nnc2c1CC[C@@H]1[C@H](CC2)[C@@H]1COC(=O)CC[NH:1][C:2]1=[C:3]2/[CH:4]=[CH:5]/[CH:6]=[N+:7]2[B-](F)(F)n2cccc21', 'COc1c(Nc2ccc(C)cc2)nc(N2CCOCC2)nc1C(=N)/C=[C:3]([CH:2]=[NH:1])/[C:4]1=[CH:5]/[CH:6]=[N+:7]1C', 'CCCCCCCC[N:3]1[CH2:2][NH:1]C[C:5](/[CH:6]=[N+:7](/C)[O-])=[CH:4]\\1', 'CCCC(CNC)OC1C=CC(C(C)[NH:1][CH2:2][C:3]2=[N+:7]=[C:6]3C=C(C(=C(C)NC4CCCC4)C(C)N)C(=O)CCC(CC)/[C:5]3=[N:4]\\2)CC1', 'CC1(CCCCCC(=O)O)c2cc(S(=O)(=O)O)ccc2[N+:7](CCCCS(=O)(=O)O)=[C:6]1/[CH:5]=[CH:4]/[CH:3]=[CH:2][NH:1]c1ccccc1', 'C=[N+:7]=[C:6]1C[N:3]([C:2](=Nc2cccc(OC(F)F)c2)[NH:1]C#N)/[N:4]=[C:5]/1c1ccc(Cl)c(C)c1', 'C#CC(=NC1C(C)=C(OCCN2CCCOCC2)C=CC1Cl)[NH:1][C:2]1=[N+:7]=[C:6]=[CH:5]/[C:4](C(F)(F)F)=[CH:3]/1', 'O/C=C/C=C/CC/C=C/CCCC[N+:7]1=[CH:6]/[CH:5]=[CH:4]/[C@:3]2(CC/C=C/CCCCCC[NH:1][CH2:2]2)C1']
        mols = [Molecule.from_string(b, 'smi').get_embedded_molecule() for b in systems]
        vd = 100
        offset = 25
        lineskip = 30
        nrows = len(mols) // 2
        baseline = nrows // 2 * -lineskip
        fs = 8
        fig = mols[0].plot(backend='rdkit', plot_range=10.1 * np.array([[-1, 1], [-1, 1], [-1, 1]]), image_size=[1500, 1500], view_settings={'view_distance': vd, 'view_center': [-offset, baseline, 0]}, font_size=fs)
        for i, m in enumerate(mols[1:]):
            i = i + 1
            i, j = (i // 2, i % 2)
            m.plot(figure=fig, backend='rdkit', view_settings={'view_distance': vd, 'view_center': [-offset if j == 0 else offset, baseline + lineskip * i, 0]}, font_size=fs)
        fig.show()
