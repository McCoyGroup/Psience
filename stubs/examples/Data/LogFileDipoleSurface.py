"""Extracted from DataTests.test_LogFileDipoleSurface via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest DataTests.test_LogFileDipoleSurface"""

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
    def test_LogFileDipoleSurface(self):
        log = TestManager.test_data('water_OH_scan.log')
        conv = lambda x: cartesian_to_zmatrix(x, ordering=[[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]]).coords[:, 0, 0]
        surf = DipoleSurface.from_log_file(log, conv)
        dips = surf(np.arange(0.5, 2, 0.1))
        self.assertEquals(dips.shape, ((2 - 0.5) / 0.1, 3))
