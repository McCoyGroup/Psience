
import os as _test_os
import unittest as _test_unittest
_RUN_EXTRAS = _test_os.environ.get("MCUTILS_TEST_EXTRAS") == "1"

import itertools
import os.path
import pprint
import numpy as np

from Peeves.TestUtils import *
from unittest import TestCase
import McUtils.Plots as plt
from Psience.Plots import *

class PlotTests(TestCase):

    @_test_unittest.skipUnless(_RUN_EXTRAS, "optional dependencies; set MCUTILS_TEST_EXTRAS=1")
    @validationTest
    def test_EnergyLevels(self):
        ploot = plot_energy_levels(
            {
                "A": {
                    "x":0,
                    "y":[0, 1, 2],
                    'color':'red'
                },
                "AC": {
                    "x":1,
                    "y":[.2, .8, 2.2],
                    'color':'green'
                },
                "B": {
                    "x":4,
                    "y":[0, 1, 4],
                    'color':'blue'
                }
            }
        )

        # ploot = plot_energy_levels([
        #     [0, 1, 2],
        #     [0, 1, 3]
        # ])

        ploot.show()

    @_test_unittest.skipUnless(_RUN_EXTRAS, "optional dependencies; set MCUTILS_TEST_EXTRAS=1")
    @validationTest
    def test_ConnectEnergyLevels(self):
        plot_energy_levels({
            'r': [0., 1.85532628],
            'ts': [19.7746029, 19.80986744],
            'prod': [-25.01321145379916]
        }, connect=True).show()

    @_test_unittest.skipUnless(_RUN_EXTRAS, "optional dependencies; set MCUTILS_TEST_EXTRAS=1")
    @validationTest
    def test_PlotOnOtherPlot(self):

        fun = lambda x: -(np.cos((1 - np.array(x)) * 2*np.pi) - (1 - np.array(x)))

        pts = [0, .475, 1]
        vals = fun(pts)

        f = plt.Plot(
            np.linspace(0, 1, 20),
            fun(np.linspace(0, 1, 20)),
        )
        plot_energy_levels({
            'r': {'x': pts[0], 'y':[vals[0]]},
            'ts': {'x': pts[1], 'y':[vals[1]]},
            'prod': {'x': pts[2], 'y':[vals[2]]},
        }, connect=True, figure=f, bar_spacing=.4).show()


    @_test_unittest.skipUnless(_RUN_EXTRAS, "optional dependencies; set MCUTILS_TEST_EXTRAS=1")
    def test_DisorederedXRange(self):
        fun = lambda x: -(np.cos((1 - np.array(x)) * 2*np.pi) - (1 - np.array(x)))

        pts = [1, 0, .475, 1]
        vals = fun(pts)

        f = plt.Plot(
            np.linspace(0, 1, 20),
            fun(np.linspace(0, 1, 20)),
        )
        plot_energy_levels({
            'r': {'x': pts[0], 'y':[vals[0]]},
            'ts': {'x': pts[1], 'y':[vals[1]]},
            'prod': {'x': pts[2], 'y':[vals[2]]},
            'prod2': {'x': pts[3], 'y':[vals[3]]},
        }, connect=True, figure=f, bar_spacing=.4).show()