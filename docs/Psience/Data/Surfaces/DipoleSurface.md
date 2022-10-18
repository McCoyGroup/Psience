## <a id="Psience.Data.Surfaces.DipoleSurface">DipoleSurface</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Data/Surfaces.py#L15)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Data/Surfaces.py#L15?message=Update%20Docs)]
</div>

Provides a unified interface to working with dipole surfaces.
Currently basically no fancier than a regular surface (although with convenient loading functions), but dipole-specific
stuff could come







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.Data.Surfaces.DipoleSurface.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, mu_x, mu_y, mu_z): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Data/Surfaces/DipoleSurface.py#L21)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Data/Surfaces/DipoleSurface.py#L21?message=Update%20Docs)]
</div>

  - `mu_x`: `Surface`
    > X-component of dipole moment
  - `mu_y`: `Surface`
    > Y-component of dipole moment
  - `mu_z`: `Surface`
    > Z-component of dipole moment


<a id="Psience.Data.Surfaces.DipoleSurface.get_log_values" class="docs-object-method">&nbsp;</a> 
```python
get_log_values(log_file, keys=('StandardCartesianCoordinates', 'DipoleMoments')): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Data/Surfaces/DipoleSurface.py#L41)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Data/Surfaces/DipoleSurface.py#L41?message=Update%20Docs)]
</div>


<a id="Psience.Data.Surfaces.DipoleSurface.from_log_file" class="docs-object-method">&nbsp;</a> 
```python
from_log_file(log_file, coord_transf, keys=('StandardCartesianCoordinates', 'DipoleMoments'), tol=0.001, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Data/Surfaces/DipoleSurface.py#L51)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Data/Surfaces/DipoleSurface.py#L51?message=Update%20Docs)]
</div>
Loads dipoles from a Gaussian log file and builds a dipole surface by interpolating.
Obviously this only really works if we have a subset of "scan" coordinates, so at this stage the user is obligated
to furnish a function that'll take a set of Cartesian coordinates and convert them to "scan" coordinates.
Coordinerds can be helpful with this, as it provides a convenient syntax for Cartesian <-> ZMatrix conversions
  - `log_file`: `str`
    > a Gaussian log file to pull from
  - `:returns`: `_`
    >


<a id="Psience.Data.Surfaces.DipoleSurface.get_fchk_values" class="docs-object-method">&nbsp;</a> 
```python
get_fchk_values(fchk_file): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Data/Surfaces/DipoleSurface.py#L106)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Data/Surfaces/DipoleSurface.py#L106?message=Update%20Docs)]
</div>


<a id="Psience.Data.Surfaces.DipoleSurface.from_fchk_file" class="docs-object-method">&nbsp;</a> 
```python
from_fchk_file(fchk_file, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Data/Surfaces/DipoleSurface.py#L118)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Data/Surfaces/DipoleSurface.py#L118?message=Update%20Docs)]
</div>
Loads dipoles from a Gaussian formatted checkpoint file and builds a dipole surface via a linear approximation
  - `fchk_file`: `Any`
    > a Gaussian fchk file to pull from
  - `log_file`: `str`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Data.Surfaces.DipoleSurface.__call__" class="docs-object-method">&nbsp;</a> 
```python
__call__(self, gridpoints, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Data/Surfaces/DipoleSurface.py#L146)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Data/Surfaces/DipoleSurface.py#L146?message=Update%20Docs)]
</div>
Explicitly overrides the Surface-level evaluation because we know the Taylor surface needs us to flatten our gridpoints
  - `gridpoints`: `Any`
    > 
  - `opts`: `Any`
    > 
  - `:returns`: `_`
    >
 </div>
</div>




## Examples













<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#Tests-ce3856" markdown="1"> Tests</a> <a class="float-right" data-toggle="collapse" href="#Tests-ce3856"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Tests-ce3856" markdown="1">
 - [FChkFileDipoleSurface](#FChkFileDipoleSurface)
- [LogFileDipoleSurface](#LogFileDipoleSurface)

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
### <a class="collapse-link" data-toggle="collapse" href="#Setup-d26cf3" markdown="1"> Setup</a> <a class="float-right" data-toggle="collapse" href="#Setup-d26cf3"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Setup-d26cf3" markdown="1">
 
Before we can run our examples we should get a bit of setup out of the way.
Since these examples were harvested from the unit tests not all pieces
will be necessary for all situations.

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






---


<div markdown="1" class="text-secondary">
<div class="container">
  <div class="row">
   <div class="col" markdown="1">
**Feedback**   
</div>
   <div class="col" markdown="1">
**Examples**   
</div>
   <div class="col" markdown="1">
**Templates**   
</div>
   <div class="col" markdown="1">
**Documentation**   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
</div>
  <div class="row">
   <div class="col" markdown="1">
[Bug](https://github.com/McCoyGroup/Psience/issues/new?title=Documentation%20Improvement%20Needed)/[Request](https://github.com/McCoyGroup/Psience/issues/new?title=Example%20Request)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Data/Surfaces/DipoleSurface.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Data/Surfaces/DipoleSurface.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Data/Surfaces/DipoleSurface.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Data/Surfaces/DipoleSurface.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Data/Surfaces.py#L15?message=Update%20Docs)   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
</div>
</div>
</div>