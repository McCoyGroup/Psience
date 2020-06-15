from Peeves.TestUtils import *
from unittest import TestCase
from PyVPT import *
from McUtils.Plots import *
from McUtils.Data import UnitsData
import sys, os, numpy as np

class VPTTests(TestCase):

    @classmethod
    def setUpClass(cls):
        TestManager.test_data_dir = os.path.join(
            os.path.dirname(TestManager.test_dir),
            "PyVPT",
            "Tests",
            "TestData"
        )


    @debugTest
    def test_WaterVPT(self):

        hammer = PerturbationTheoryHamiltonian.from_fchk(TestManager.test_data("HOD_freq.fchk"), n_quanta=5)
        h2w = UnitsData.convert("Hartrees", "Wavenumbers")
        get_ind = hammer.get_state_indices
        get_qn = hammer.get_state_quantum_numbers
        np.set_printoptions(suppress=True)

        # np.savetxt("/Users/Mark/Desktop/H0.txt", hammer.H0[:, :])
        # np.savetxt("/Users/Mark/Desktop/H1.txt", hammer.H1[:, :])
        # np.savetxt("/Users/Mark/Desktop/H2.txt", hammer.H2[:, :])
        # coeffs, corrs = hammer.get_corrections([0, 1])

        states = ((0, 0, 0), (0, 0, 1), (0, 1, 0), (1, 0, 0))
        "100 300 010 210 120 030 201 001 111 021 102 012 003"
        c_states = (
            (0, 0, 0),
            (0, 0, 1),
            (0, 0, 3),
            (0, 1, 0),
            (0, 1, 2),
            (0, 2, 1),
            (0, 3, 0),
            (1, 0, 2),
            (1, 0, 0),
            (1, 1, 1),
            (1, 2, 0),
            (2, 0, 1),
            (2, 1, 0),
            (3, 0, 0)
            )
        # print(hammer.H1[get_ind([c_states[0]]), get_ind([c_states[4]])])
        # cf, cr = hammer.get_corrections(c_states[:2], coupled_states=c_states)
        # print(h2w*cr[1])

        hammer.H1[get_ind([0, 1, 0]), get_ind([0, 4, 0])]

        thresh = 300 / h2w
        coeffs, corrs = hammer.get_corrections(states, coupled_states=None)
        energies = h2w*sum(corrs)
        corrs = [h2w*x for x in corrs]

        print("State Energies:\n",
              *(
                  "{:<1.0f} {:<1.0f} {:<1.0f} {:>5.0f} {:>5.0f} {:>5.0f} {:>5.0f}\n".format(*x) for x in np.round(
                  np.column_stack((states, corrs[0], energies, corrs[0] - corrs[0][0], (energies - energies[0])))
              ))
              )


        # print(h2w*hammer.martin_test())

""" 
 n m l  E(Ha) E(An) F(Ha) F(An)
------------------------------- 
        Gaussian:
 0 0 0  4053  3995     0     0
 0 0 1     -     -  3874  3685
 0 1 0     -     -  2810  2706
 1 0 0     -     -  1422  1383
        Mark
 0 0 0  4053  4003     0     0
 0 0 1  7927  7682  3874  3679
 0 1 0  6863  6707  2810  2704
 1 0 0  5475  5377  1422  1374
        Anne
 0 0 0  4053  4003     0     0
 0 0 1  7927  7682  3874  3660
 0 1 0  6863  6707  2810  2687
 1 0 0  5475  5377  1422  1374
 """

        # print(corrs[0])
        # print(corrs[1], coeffs)
        # print(corrs[2])



