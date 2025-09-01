import itertools
import os.path
import pprint

from Peeves.TestUtils import *
from unittest import TestCase
from Psience.Plots import *

class PlotTests(TestCase):

    @debugTest
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