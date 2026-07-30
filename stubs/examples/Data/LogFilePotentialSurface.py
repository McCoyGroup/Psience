"""Extracted from DataTests.test_LogFilePotentialSurface via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest DataTests.test_LogFilePotentialSurface"""

from Peeves.TestUtils import *
from Psience.Data import *
from McUtils.Coordinerds import cartesian_to_zmatrix
from McUtils.Plots import *
import McUtils.Numputils as nput
from unittest import TestCase
import sys, h5py, math, numpy as np

class DataTests(TestCase):
    maxDiff = None

    @validationTest
    def test_LogFilePotentialSurface(self):
        log = TestManager.test_data('water_OH_scan.log')
        conv = lambda x: np.linalg.norm(x[:, 0] - x[:, 1], axis=1)
        surf = PotentialSurface.from_log_file(log, conv)
        pots = surf(np.arange(0.5, 2, 0.1))
        self.assertEquals(pots.shape, ((2 - 0.5) / 0.1,))
