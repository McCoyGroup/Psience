"""Extracted from MolecoolsTests.test_BondGraphZMatrixIssues via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest MolecoolsTests.test_BondGraphZMatrixIssues"""

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
    def test_BondGraphZMatrixIssues(self):
        horizontal_structure = '\n  C   -0.0531542   -0.2705697   -0.4038580\n  C   -1.2672917   -0.9889447   -0.5655902\n  C   -2.5115504   -0.3220336   -0.4206904\n  C   -3.7326727   -1.0553190   -0.5317414\n  C   -3.7004761   -2.4495575   -0.8014742\n  C   -2.4508476   -3.1180666   -0.9547299\n  C   -2.4145647   -4.5238734   -1.2198701\n  C   -3.6621332   -5.2283272   -1.3328590\n  C   -4.8569665   -4.5892997   -1.1846209\n  H   -5.7886526   -5.1430158   -1.2712961\n  C   -4.9257958   -3.1802597   -0.9094725\n  C   -6.1422696   -2.5061569   -0.7417504\n  H   -7.0716942   -3.0658565   -0.8238723\n  C   -6.1964402   -1.1330449   -0.4675717\n  C   -7.4373714   -0.4323877   -0.2841295\n  H   -8.3634594   -0.9960396   -0.3668735\n  C   -7.4675572    0.9027752   -0.0104014\n  H   -8.4176460    1.4127807    0.1290235\n  C   -6.2596814    1.6721378    0.1065559\n  C   -6.2671014    3.0419631    0.4027137\n  H   -7.2201538    3.5465627    0.5477102\n  C   -5.0827029    3.7776905    0.5303886\n  C   -5.0759603    5.1765715    0.8640348\n  H   -6.0312166    5.6798428    0.9938500\n  C   -3.9096314    5.8651178    1.0239392\n  H   -3.9280356    6.9205177    1.2854065\n  C   -2.6313319    5.2261265    0.8675107\n  C   -2.6065567    3.8427713    0.5018790\n  C   -3.8255173    3.1208176    0.3413015\n  C   -3.7956779    1.7351059    0.0306407\n  C   -5.0084680    1.0043192   -0.0784111\n  C   -4.9769052   -0.3925329   -0.3620967\n  C   -2.5425203    1.0693727   -0.1391634\n  C   -1.3290713    1.7891599    0.0179298\n  C   -0.0839449    1.1156505   -0.0953746\n  C    1.1355172    1.8214190    0.1444458\n  C    1.1010171    3.1946003    0.5049308\n  C   -0.1474428    3.8757259    0.6003097\n  C   -1.3601817    3.1780730    0.3542430\n  C   -0.1858234    5.2499558    0.9952833\n  C   -1.4242191    5.8941573    1.1061277\n  H   -1.4495383    6.9350662    1.4194308\n  C    1.0532032    5.9119919    1.2917556\n  C    2.2465413    5.2622367    1.1941631\n  C    2.3197271    3.8832904    0.7977264\n  C    3.5352228    3.1943918    0.6970674\n  C    3.5930945    1.8423682    0.3340199\n  C    2.3778386    1.1384405    0.0625407\n  C    2.4102323   -0.2468567   -0.2728220\n  C    1.1993952   -0.9519230   -0.5029697\n  C    1.2283278   -2.3432451   -0.7890186\n  C    0.0098448   -3.0637208   -0.9542032\n  C   -1.2361383   -2.3926339   -0.8313463\n  C    0.0349492   -4.4697870   -1.2200387\n  C   -1.1751012   -5.1636932   -1.3514989\n  H   -1.1509809   -6.2327171   -1.5521573\n  C    1.3123634   -5.1187284   -1.3301511\n  H    1.3318089   -6.1850993   -1.5412707\n  C    2.4774035   -4.4302726   -1.1665430\n  C    2.4838369   -3.0225639   -0.8768600\n  C    3.6676819   -2.3042610   -0.6642272\n  C    3.6602742   -0.9382196   -0.3533601\n  C    4.8654746   -0.1983466   -0.0980494\n  H    5.8151245   -0.7232786   -0.1667677\n  C    4.8334523    1.1238831    0.2322146\n  H    5.7571926    1.6620172    0.4300728\n  H    4.6199581   -2.8265000   -0.7288371\n  H    3.4322477   -4.9442495   -1.2457993\n  H    4.4607529    3.7230634    0.9146432\n  H    3.1731065    5.7801519    1.4281981\n  H    1.0166878    6.9498086    1.6110402\n  H   -3.6334126   -6.2956952   -1.5377463\n  N   -0.9768306   -0.8645727    2.6960567\n  C   -1.1153261    0.4844050    3.0234710\n  C   -0.0778126    1.3920491    3.2650545\n  C   -0.3249244    2.7228872    3.5872077\n  C   -1.6325156    3.2294720    3.7094695\n  C   -2.6582688    2.3098377    3.4306428\n  C   -2.4217818    0.9854865    3.0975439\n  F   -3.4554562    0.1841828    2.8127935\n  F   -3.9567056    2.6815385    3.4512927\n  C   -1.8441701    4.6700822    4.1229961\n  O   -0.9060656    5.3915137    4.4314908\n  N   -3.1391561    5.1091075    4.1948416\n  H   -3.8828970    4.6489540    3.6932711\n  H   -3.2290700    6.1045076    4.3462867\n  F    0.7511767    3.4882777    3.7556598\n  F    1.1959563    0.9736924    3.1611571\n  N    0.1556565   -1.3754611    2.6214998\n  N    1.0867607   -2.0159052    2.5061938'
        horp = Molecule.from_string(horizontal_structure, units='Angstroms')
        f1: Molecule = horp.fragments[1]
        zmat = f1.get_bond_zmatrix()
        hmm = [z[0] for z in zmat[:8]]
        segs = f1.find_backbone_segments()
        f1.plot(highlight_atoms=segs[0]).show()
