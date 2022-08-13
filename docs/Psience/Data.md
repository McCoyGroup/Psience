# <a id="Psience.Data">Psience.Data</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/tree/master/Psience/Data)]
</div>
    
Provides core Data-related types and data structures.
Intended to be high-level data as opposed to the lower-level stuff in `McUtils.Data`.
That means including stuff like dipole and potential energy surfaces that know how to compute their own properties.
Currently...well that's all we have. But wrappers for commonly-used potentials & bases could well come.
Not sure at this point, though.

<div class="container alert alert-secondary bg-light">
  <div class="row">
   <div class="col" markdown="1">
[DipoleSurface](Data/Surfaces/DipoleSurface.md)   
</div>
   <div class="col" markdown="1">
[PotentialSurface](Data/Surfaces/PotentialSurface.md)   
</div>
   <div class="col" markdown="1">
[KEData](Data/KEData/KEData.md)   
</div>
</div>
  <div class="row">
   <div class="col" markdown="1">
[KEDataHandler](Data/KEData/KEDataHandler.md)   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
</div>
</div>






<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
### <a class="collapse-link" data-toggle="collapse" href="#tests">Tests</a> <a class="float-right" data-toggle="collapse" href="#tests"><i class="fa fa-chevron-down"></i></a>
 </div>
<div class="collapsible-section collapsible-section-body collapse show" id="tests" markdown="1">

- [FChkFileDipoleSurface](#FChkFileDipoleSurface)
- [LogFileDipoleSurface](#LogFileDipoleSurface)
- [LogFilePotentialSurface](#LogFilePotentialSurface)

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
#### <a class="collapse-link" data-toggle="collapse" href="#test-setup">Setup</a> <a class="float-right" data-toggle="collapse" href="#test-setup"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="test-setup" markdown="1">

Before we can run our examples we should get a bit of setup out of the way.
Since these examples were harvested from the unit tests not all pieces
will be necessary for all situations.
```python
from Peeves.TestUtils import *
from Psience.Data import *
from McUtils.Coordinerds import cartesian_to_zmatrix
from McUtils.Plots import *
from unittest import TestCase
import sys, h5py, math, numpy as np
```

All tests are wrapped in a test class
```python
class DataTests(TestCase):
    maxDiff = None
```

 </div>
</div>

#### <a name="FChkFileDipoleSurface">FChkFileDipoleSurface</a>
```python
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
```
#### <a name="LogFileDipoleSurface">LogFileDipoleSurface</a>
```python
    def test_LogFileDipoleSurface(self):
        log = TestManager.test_data("water_OH_scan.log")
        conv = lambda x: cartesian_to_zmatrix(
            x, ordering=[[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]]
        ).coords[:, 0, 0]
        surf = DipoleSurface.from_log_file(log, conv)
        dips = surf(np.arange(.5, 2, .1))
        self.assertEquals(dips.shape, ((2-.5)/.1, 3))
```
#### <a name="LogFilePotentialSurface">LogFilePotentialSurface</a>
```python
    def test_LogFilePotentialSurface(self):
        log = TestManager.test_data("water_OH_scan.log")
        conv = lambda x: np.linalg.norm(x[:, 0] - x[:, 1], axis=1)
        surf = PotentialSurface.from_log_file(log, conv)
        pots = surf(np.arange(.5, 2, .1))
        self.assertEquals(pots.shape, ((2-.5)/.1,))
```

 </div>
</div>

___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/master/ci/examples/Psience/Data.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/master/?filename=ci/examples/Psience/Data.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/master/ci/docs/Psience/Data.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/master/?filename=ci/docs/templates/Psience/Data.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/__init__.py?message=Update%20Docs)