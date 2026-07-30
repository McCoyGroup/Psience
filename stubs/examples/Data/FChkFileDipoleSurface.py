"""Extracted from DataTests.test_FChkFileDipoleSurface via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest DataTests.test_FChkFileDipoleSurface"""

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
    def test_FChkFileDipoleSurface(self):
        fchk = TestManager.test_data('HOD_freq.fchk')
        surf = DipoleSurface.from_fchk_file(fchk)
        surf_center = surf.surfs[0].base.data['center']
        self.assertIsInstance(surf_center, np.ndarray)
        self.assertTrue(np.allclose(surf(surf_center) - np.array([s.base.data['ref'] for s in surf.surfs]), 0.0))
        self.assertEquals(surf([[0, 0, 0], [1, 0, 0], [0, 1, 0]]).shape, (1, 3))
        self.assertEquals(surf([[[0, 0, 0], [1, 0, 0], [0, 1, 0]], [[0, 0, 0], [1, 0, 0], [0, 1, 0]]]).shape, (2, 3))
