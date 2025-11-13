from Peeves.TestUtils import *
from unittest import TestCase
from Psience.Nonlinear import *
import McUtils.Plots as plt
from McUtils.Data import UnitsData
import sys, os, numpy as np
import copy

class NonlinearTests(TestCase):

    @debugTest
    def test_BasicPathways(self):
        from Psience.BasisReps import BasisStateSpace, HarmonicOscillatorProductBasis
        initial_basis = BasisStateSpace.from_quanta(
            HarmonicOscillatorProductBasis(1),
            0
        )
        # pathways = liouville_pathways(3)

        responses = experimental_response_generator(
            {
                ((0,0), (0,1)):{
                    'frequency':1610,
                    'transition_moment':[0, 0, .1]
                },
                ((0,0), (1,0)):{
                    'frequency':1570,
                    'transition_moment':[0, .1, 0]
                },
                ((0,1), (0,2)): {
                    'frequency': 1600,
                    'transition_moment': [0, 0, .12]
                },
                ((1,0), (2, 0)): {
                    'frequency': 1560,
                    'transition_moment': [0, .12, 0]
                },
                ((0, 1), (1, 1)): {
                    'frequency': 1570,
                    'transition_moment': [0, .1, 0]
                },
                ((1, 0), (1, 1)): {
                    'frequency': 1610,
                    'transition_moment': [0,  0, .1]
                }
            },
            band_coherences={
                (0, 1):2,
                (1, 2):2
            },
            frequency_unit="Wavenumbers",
            application_domain="frequency",
            polarization=None,
            # response_function_generator="alt2d",
            response_function_generator="simple2dir",
            # included_signals=['non-rephasing']
        )

        spec = responses.get_spectrum(
            # [1550, 1620], 0, [1550, 1620],
            [1400, 1700], 0, [1400, 1700],
            num_samples=1024,
            time_step=.0625,
            default_frequency_divisions=500
        )
        # subspec = spec.frequency_filter([1400, 1700], [1400, 1700])
        subspec = spec.frequency_filter([1550, 1620], [1550, 1620])
        # print(subspec.clip(1e-6, 1))
        subspec.clip(1e-8, 1).plot(levels=15).show()


        # for p in paths[1]:
        #     print(p)

        # for cur, nxt in enumerate_basis_path(initial_basis, pathways[0]):
        #     print("-"*50)
        #     print(cur.excitations)
        #     print("."*10)
        #     print(nxt.excitations)