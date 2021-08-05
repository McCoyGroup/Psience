# <a id="Psience.Data">Psience.Data</a>
    
Provides core Data-related types and data structures.
Intended to be high-level data as opposed to the lower-level stuff in `McUtils.Data`.
That means including stuff like dipole and potential energy surfaces that know how to compute their own properties.
Currently...well that's all we have. But wrappers for commonly-used potentials & bases could well come.
Not sure at this point, though.

### Members:

  - [DipoleSurface](Data/Surfaces/DipoleSurface.md)
  - [PotentialSurface](Data/Surfaces/PotentialSurface.md)

### Examples:



```python

from Peeves.TestUtils import *
from Psience.Data import *
from McUtils.Coordinerds import cartesian_to_zmatrix
from McUtils.Plots import *
from unittest import TestCase
import sys, h5py, math, numpy as np

class PsienceDataTests(TestCase):

    @debugTest
    def test_FChkFileDipoleSurface(self):
        fchk = TestManager.test_data("HOD_freq.fchk")
        surf = DipoleSurface.from_fchk_file(fchk)
        surf_center = surf.surfs[0].base.data['center']
        self.assertIsInstance(surf_center, np.ndarray)
        self.assertTrue(
            np.allclose(surf(surf_center) - np.array([s.base.data['ref'] for s in surf.surfs]), 0.)
        )
        self.assertEquals(surf([[0, 0, 0], [1, 0, 0], [0, 1, 0]]).shape, (1, 3))
        self.assertEquals(surf([
            [[0, 0, 0], [1, 0, 0], [0, 1, 0]],
            [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
        ]).shape, (2, 3))

    @debugTest
    def test_LogFileDipoleSurface(self):
        log = TestManager.test_data("water_OH_scan.log")
        conv = lambda x: cartesian_to_zmatrix(
            x, ordering=[[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]]
        ).coords[:, 0, 0]
        surf = DipoleSurface.from_log_file(log, conv)
        dips = surf(np.arange(.5, 2, .1))
        self.assertEquals(dips.shape, ((2-.5)/.1, 3))

        # surf_center = surf.surfs[0].base.data['center']
        # self.assertIsInstance(surf_center, np.ndarray)
        # self.assertTrue(
        #     np.allclose(surf(surf_center) - np.array([s.base.data['ref'] for s in surf.surfs]), 0.)
        # )
        # self.assertEquals(surf([[0, 0, 0], [1, 0, 0], [0, 1, 0]]).shape, (1, 3))
        # self.assertEquals(surf([
        #     [[0, 0, 0], [1, 0, 0], [0, 1, 0]],
        #     [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
        # ]).shape, (2, 3))

    @debugTest
    def test_LogFilePotentialSurface(self):
        log = TestManager.test_data("water_OH_scan.log")
        conv = lambda x: np.linalg.norm(x[:, 0] - x[:, 1], axis=1)
        surf = PotentialSurface.from_log_file(log, conv)
        pots = surf(np.arange(.5, 2, .1))
        self.assertEquals(pots.shape, ((2-.5)/.1,))
```