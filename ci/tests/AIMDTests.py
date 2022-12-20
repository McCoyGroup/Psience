# import itertools
#
# import scipy.linalg
#
# import McUtils.Zachary
from Peeves.TestUtils import *
from unittest import TestCase

# from McUtils.Data import UnitsData, PotentialData
# from McUtils.Zachary import Interpolator, Symbols
# import McUtils.Plots as plt
# from McUtils.GaussianInterface import GaussianLogReader

from McUtils.Extensions import ModuleLoader
from McUtils.Zachary import FiniteDifferenceDerivative

from Psience.AIMD import *
from Psience.Molecools import Molecule
import numpy as np

class AIMDTests(TestCase):

    # @debugTest
    # def test_Morse(self):
    #
    #     sym = Symbols('r')
    #     fn = sym.morse(sym.r, de=10, a=1)
    #
    #     pot = PairwisePotential(fn)
    #
    #     coords = np.random.rand(7, 3)
    #
    #     r = np.arange(7)
    #     ri, ci = np.triu_indices(7, 1)
    #
    #     forces = lambda c:pot.forces(coords)
    #
    #     sim = AIMDSimulator()

    @debugTest
    def test_mbpol(self):
        # sym = Symbols('r')
        # fn = sym.morse(sym.r, de=10, a=1)
        #
        # pot = PairwisePotential(fn)

        loader = ModuleLoader(TestManager.current_manager().test_data_dir)
        mbpol = loader.load("LegacyMBPol").MBPol

        r1 = np.random.normal(1.8, .3, 1000)
        r2 = np.random.normal(1.8, .3, 1000)
        coords = np.array([
            [
                [0, 0, 0],
                [a, 0, 0],
                [0, b, 0]
            ] for a, b in zip(r1, r2)
        ])
        forces = lambda c: mbpol.get_pot_grad(nwaters=1, coords=c.reshape(-1, 3, 3), threading_vars=['energy', 'grad', 'coords'], threading_mode='omp')['grad'].reshape(
            c.shape
        )

        sim = AIMDSimulator(["O", "H", "H"], coords, force_function=forces)

        traj = np.array(sim.propagate(1000))
        engs = mbpol.get_pot(nwaters=1, coords=traj.reshape((-1, 3, 3)), threading_vars=['energy', 'coords']).reshape(traj.shape[:2])
        grads = forces(traj)
        hessians = FiniteDifferenceDerivative(forces, function_shape=((3, 3), (3, 3)) ).derivatives(traj).compute_derivatives(1)

        raise Exception(hessians.shape)
        # raise Exception(traj.shape, engs)



