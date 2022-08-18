## <a id="Psience.Data.Surfaces.DipoleSurface">DipoleSurface</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Data/Surfaces.py#L15)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/Surfaces.py#L15?message=Update%20Docs)]
</div>

Provides a unified interface to working with dipole surfaces.
Currently basically no fancier than a regular surface (although with convenient loading functions), but dipole-specific
stuff could come

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

<a id="Psience.Data.Surfaces.DipoleSurface.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, mu_x, mu_y, mu_z): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Data/Surfaces.py#L21)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/Surfaces.py#L21?message=Update%20Docs)]
</div>


- `mu_z`: `Surface`
    >Z-component of dipole moment
- `mu_y`: `Surface`
    >Y-component of dipole moment
- `mu_x`: `Surface`
    >X-component of dipole moment

<a id="Psience.Data.Surfaces.DipoleSurface.get_log_values" class="docs-object-method">&nbsp;</a> 
```python
get_log_values(log_file, keys=('StandardCartesianCoordinates', 'DipoleMoments')): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Data/Surfaces.py#L41)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/Surfaces.py#L41?message=Update%20Docs)]
</div>

<a id="Psience.Data.Surfaces.DipoleSurface.from_log_file" class="docs-object-method">&nbsp;</a> 
```python
from_log_file(log_file, coord_transf, keys=('StandardCartesianCoordinates', 'DipoleMoments'), tol=0.001, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Data/Surfaces.py#L51)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/Surfaces.py#L51?message=Update%20Docs)]
</div>

Loads dipoles from a Gaussian log file and builds a dipole surface by interpolating.
Obviously this only really works if we have a subset of "scan" coordinates, so at this stage the user is obligated
to furnish a function that'll take a set of Cartesian coordinates and convert them to "scan" coordinates.
Coordinerds can be helpful with this, as it provides a convenient syntax for Cartesian <-> ZMatrix conversions
- `:returns`: `_`
    >
- `log_file`: `str`
    >a Gaussian log file to pull from

<a id="Psience.Data.Surfaces.DipoleSurface.get_fchk_values" class="docs-object-method">&nbsp;</a> 
```python
get_fchk_values(fchk_file): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Data/Surfaces.py#L106)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/Surfaces.py#L106?message=Update%20Docs)]
</div>

<a id="Psience.Data.Surfaces.DipoleSurface.from_fchk_file" class="docs-object-method">&nbsp;</a> 
```python
from_fchk_file(fchk_file, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Data/Surfaces.py#L118)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/Surfaces.py#L118?message=Update%20Docs)]
</div>

Loads dipoles from a Gaussian formatted checkpoint file and builds a dipole surface via a linear approximation
- `:returns`: `_`
    >
- `log_file`: `str`
    >
- `fchk_file`: `Any`
    >a Gaussian fchk file to pull from

<a id="Psience.Data.Surfaces.DipoleSurface.__call__" class="docs-object-method">&nbsp;</a> 
```python
__call__(self, gridpoints, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Data/Surfaces.py#L146)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/Surfaces.py#L146?message=Update%20Docs)]
</div>

Explicitly overrides the Surface-level evaluation because we know the Taylor surface needs us to flatten our gridpoints
- `:returns`: `_`
    >
- `opts`: `Any`
    >
- `gridpoints`: `Any`
    >

 </div>
</div>





<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
### <a class="collapse-link" data-toggle="collapse" href="#tests">Tests</a> <a class="float-right" data-toggle="collapse" href="#tests"><i class="fa fa-chevron-down"></i></a>
 </div>
<div class="collapsible-section collapsible-section-body collapse show" id="tests" markdown="1">

- [FChkFileDipoleSurface](#FChkFileDipoleSurface)
- [LogFileDipoleSurface](#LogFileDipoleSurface)

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

 </div>
</div>

___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Data/Surfaces/DipoleSurface.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Data/Surfaces/DipoleSurface.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Data/Surfaces/DipoleSurface.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Data/Surfaces/DipoleSurface.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/Surfaces.py#L15?message=Update%20Docs)