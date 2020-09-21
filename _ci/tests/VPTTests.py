from Peeves.TestUtils import *
from unittest import TestCase
from Psience.VPT2 import *
from McUtils.Data import UnitsData
import sys, os, numpy as np

class VPTTests(TestCase):

    @inactiveTest
    def test_WaterVPTInternals(self):

        internals = [
            [0, -1, -1, -1],
            [1,  0, -1, -1],
            [2,  0,  1, -1]
        ]

        hammer = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("HOD_freq.fchk"),
            n_quanta=6,
            internals=internals
        )
        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        states = (
            (0, 0, 0),
            (0, 0, 1), (0, 1, 0), (1, 0, 0),
            (0, 0, 2), (0, 2, 0), (2, 0, 0),
            (0, 1, 1), (1, 0, 1), (1, 1, 0)
        )
        coeffs, corrs = hammer.get_corrections(states, coupled_states=None)
        energies = h2w * sum(corrs)
        corrs = [h2w * x for x in corrs]

        gaussian_states = [(0, 0, 1), (0, 1, 0), (1, 0, 0),
                           (0, 0, 2), (0, 2, 0), (2, 0, 0),
                           (0, 1, 1), (1, 0, 1), (1, 1, 0)]
        gaussian_energies = np.array([4052.912, 3994.844])
        gaussian_freqs = np.array([
            [3873.846, 3685.815],
            [2810.031, 2706.132],
            [1421.946, 1383.391],
            [7747.692, 7202.835],
            [5620.062, 5323.917],
            [2843.893, 2749.027],
            [6683.877, 6377.958],
            [5295.792, 5044.721],
            [4231.977, 4072.407]
        ])

        my_energies = np.array([corrs[0][0], energies[0]])
        my_freqs = np.column_stack([
            corrs[0][1:] - corrs[0][0],
            energies[1:] - energies[0]
        ])

        print("Gaussian:\n",
              "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*gaussian_energies, "-", "-"),
              *(
                  "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                  zip(gaussian_states, gaussian_freqs)
              )
              )
        print("State Energies:\n",
              "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*my_energies, "-", "-"),
              *(
                  "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                  zip(states[1:], my_freqs)
              )
              )
        print("Difference Energies:\n",
              "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*(my_energies-gaussian_energies), "-", "-"),
              *(
                  "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                  zip(states[1:], my_freqs-gaussian_freqs)
              )
              )

    @inactiveTest
    def test_WaterVPTInternalsHOT(self):
        internals = [
            [0, -1, -1, -1],
            [1, 0, -1, -1],
            [2, 0, 1, -1]
        ]

        hammer = PerturbationTheoryHamiltonian.from_fchk(
            TestManager.test_data("HOT_freq.fchk"),
            n_quanta=6,
            internals=internals
        )
        h2w = UnitsData.convert("Hartrees", "Wavenumbers")

        states = (
            (0, 0, 0),
            (0, 0, 1), (0, 1, 0), (1, 0, 0),
            (0, 0, 2), (0, 2, 0), (2, 0, 0),
            (0, 1, 1), (1, 0, 1), (1, 1, 0)
        )
        coeffs, corrs = hammer.get_corrections(states, coupled_states=None)
        energies = h2w * sum(corrs)
        corrs = [h2w * x for x in corrs]

        gaussian_states = [(0, 0, 1), (0, 1, 0), (1, 0, 0),
            (0, 0, 2), (0, 2, 0), (2, 0, 0),
            (0, 1, 1), (1, 0, 1), (1, 1, 0)]
        gaussian_energies = np.array([3789.117, 3736.855])
        gaussian_freqs = np.array([
            [3872.325, 3692.047],
            [2357.491, 2285.849],
            [1348.418, 1313.165],
            [7744.651, 7214.810],
            [4714.982, 4509.256],
            [2696.835, 2607.130],
            [6229.816, 5973.755],
            [5220.743, 4987.363],
            [3705.909, 3584.756]
        ])

        my_energies = np.array([corrs[0][0], energies[0]])
        my_freqs = np.column_stack([
            corrs[0][1:] - corrs[0][0],
            energies[1:] - energies[0]
        ])

        print("Gaussian:\n",
              "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*gaussian_energies, "-", "-"),
              *(
                  "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s,e in
                  zip(gaussian_states, gaussian_freqs)
              )
              )
        print("States:\n",
              "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*my_energies, "-", "-"),
              *(
                  "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                  zip(states[1:], my_freqs)
              )
              )
        print("Difference:\n",
              "0 0 0 {:>8.3f} {:>8.3f} {:>8} {:>8}\n".format(*(my_energies - gaussian_energies), "-", "-"),
              *(
                  "{:<1.0f} {:<1.0f} {:<1.0f} {:>8} {:>8} {:>8.3f} {:>8.3f}\n".format(*s, "-", "-", *e) for s, e in
                  zip(states[1:], my_freqs - gaussian_freqs)
              )
              )

    @inactiveTest
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

        # hammer.H1[get_ind([0, 1, 0]), get_ind([0, 4, 0])]

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
        print("Target:\n",
              """
 0 0 0  4053  3995     0     0
 0 0 1     -     -  3874  3685
 0 1 0     -     -  2810  2706
 1 0 0     -     -  1422  1383
 """
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



