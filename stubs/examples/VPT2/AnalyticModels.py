"""Extracted from VPT2Tests.test_AnalyticModels via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest VPT2Tests.test_AnalyticModels"""

import itertools
import tempfile
try:
    from Peeves.TestUtils import *
    from Peeves import BlockProfiler
except:
    pass
from unittest import TestCase
from Psience.VPT2 import *
from Psience.Molecools import Molecule
from Psience.BasisReps import HarmonicOscillatorProductBasis, BasisStateSpace
from McUtils.Data import UnitsData
import McUtils.Plots as plt
import McUtils.Numputils as nput
from McUtils.Scaffolding import *
from McUtils.Profilers import Timer
from McUtils.Parallelizers import Parallelizer, SerialNonParallelizer, MultiprocessingParallelizer
from McUtils.Zachary import FiniteDifferenceDerivative
import sys, os, numpy as np, itertools as ip

class VPT2Tests(TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        np.set_printoptions(linewidth=int(100000000.0))

    @validationTest
    def test_AnalyticModels(self):
        from Psience.AnalyticModels import AnalyticModel as Model
        from McUtils.Data import AtomData, UnitsData
        order = 4
        expansion_order = {'potential': 0, 'gmatrix': 4, 'pseudopotential': 4}
        hoh_params = {}
        hoh_params['mH'] = AtomData['H']['Mass'] * UnitsData.convert('AtomicMassUnits', 'AtomicUnitOfMass')
        hoh_params['mO'] = AtomData['O']['Mass'] * UnitsData.convert('AtomicMassUnits', 'AtomicUnitOfMass')
        cm2borh = UnitsData.convert('Angstroms', 'BohrRadius')
        hoh_params['re'] = 0.9575 * cm2borh
        erg2h = UnitsData.convert('Ergs', 'Hartrees')
        invcm2borh = UnitsData.convert('InverseAngstroms', 'InverseBohrRadius')
        hoh_params['De'] = 8.84e-12 * erg2h
        hoh_params['a'] = 2.175 * invcm2borh
        hoh_params['b_e'] = np.deg2rad(104.5)
        hoh_params['k_b'] = 3.2 ** 2 * 1600 * UnitsData.convert('Wavenumbers', 'Hartrees')
        model = Model([Model.r(1, 2), Model.r(2, 3), Model.a(1, 2, 3)], Model.Potential.morse(1, 2, De=hoh_params['De'], a=hoh_params['a'], re=hoh_params['re']) + Model.Potential.morse(2, 3, De=hoh_params['De'], a=hoh_params['a'], re=hoh_params['re']) + Model.Potential.harmonic(1, 2, 3, k=hoh_params['k_b'], qe=hoh_params['b_e']), values={Model.m(1): hoh_params['mH'], Model.m(2): hoh_params['mO'], Model.m(3): hoh_params['mH'], Model.r(1, 2): hoh_params['re'], Model.r(2, 3): hoh_params['re'], Model.a(1, 2, 3): hoh_params['b_e']})
        model.run_VPT(order=order, return_analyzer=False, expansion_order=expansion_order)

        class harmonically_coupled_morse:

            def __init__(self, De_1, a_1, re_1, De_2, a_2, re_2, kb, b_e):
                self.De_1 = De_1
                self.a_1 = a_1
                self.re_1 = re_1
                self.De_2 = De_2
                self.a_2 = a_2
                self.re_2 = re_2
                self.kb = kb
                self.b_e = b_e

            def __call__(self, carts):
                v1 = carts[..., 1, :] - carts[..., 0, :]
                v2 = carts[..., 2, :] - carts[..., 0, :]
                r1 = nput.vec_norms(v1) - self.re_1
                r2 = nput.vec_norms(v2) - self.re_2
                bend, _ = nput.vec_angles(v1, v2)
                bend = bend - self.b_e
                return self.De_1 * (1 - np.exp(-self.a_1 * r1)) ** 2 + self.De_2 * (1 - np.exp(-self.a_2 * r2)) ** 2 + self.kb * bend ** 2
        atoms = ['O', 'H', 'H']
        coords = np.array([[0.0, 0.0, 0.0], [hoh_params['re'], 0.0, 0.0], np.dot(nput.rotation_matrix([0, 0, 1], hoh_params['b_e']), [hoh_params['re'], 0.0, 0.0])])
        masses = np.array([AtomData[x]['Mass'] for x in atoms]) * UnitsData.convert('AtomicMassUnits', 'AtomicUnitOfMass')
        pot_file = os.path.expanduser('~/Desktop/water_pot.hdf5')
        water_chk = Checkpointer.from_file(pot_file)
        if expansion_order['potential'] > -1:
            with water_chk as wat:
                try:
                    potential_derivatives = wat['potential_derivatives']
                except (OSError, KeyError):
                    potential_function = harmonically_coupled_morse(hoh_params['De'], hoh_params['a'], hoh_params['re'], hoh_params['De'], hoh_params['a'], hoh_params['re'], hoh_params['k_b'], hoh_params['b_e'])
                    deriv_gen = FiniteDifferenceDerivative(potential_function, function_shape=((None, None), 0), stencil=5 + expansion_order['potential'], mesh_spacing=0.001).derivatives(coords)
                    potential_derivatives = deriv_gen.derivative_tensor(list(range(1, order + 3)))
                    wat['potential_derivatives'] = potential_derivatives
        else:
            potential_derivatives = []
        analyzer = VPTRunner.run_simple([atoms, coords, dict(masses=masses)], 2, potential_derivatives=potential_derivatives, calculate_intensities=False, internals=[[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]], order=order, internal_fd_mesh_spacing=0.01, cartesian_fd_mesh_spacing=0.01, checkpoint=os.path.expanduser('~/Desktop/water_analyt.hdf5'), expansion_order=expansion_order)
