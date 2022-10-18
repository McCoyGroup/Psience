## <a id="Psience.Data.Surfaces.PotentialSurface">PotentialSurface</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Data/Surfaces.py#L167)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Data/Surfaces.py#L167?message=Update%20Docs)]
</div>

A potential surface structure to go along with the DipoleSurface.
Provides convenient access to dipole data + a unified interface to things like energy minimization







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.Data.Surfaces.PotentialSurface.get_log_values" class="docs-object-method">&nbsp;</a> 
```python
get_log_values(log_file, keys=('StandardCartesianCoordinates', 'ScanEnergies')): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Data/Surfaces/PotentialSurface.py#L173)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Data/Surfaces/PotentialSurface.py#L173?message=Update%20Docs)]
</div>


<a id="Psience.Data.Surfaces.PotentialSurface.from_log_file" class="docs-object-method">&nbsp;</a> 
```python
from_log_file(log_file, coord_transf, keys=('StandardCartesianCoordinates', 'ScanEnergies'), tol=0.001, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Data/Surfaces/PotentialSurface.py#L189)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Data/Surfaces/PotentialSurface.py#L189?message=Update%20Docs)]
</div>
Loads dipoles from a Gaussian log file and builds a potential surface by interpolating.
Obviously this only really works if we have a subset of "scan" coordinates, so at this stage the user is obligated
to furnish a function that'll take a set of Cartesian coordinates and convert them to "scan" coordinates.
Coordinerds can be helpful with this, as it provides a convenient syntax for Cartesian <-> ZMatrix conversions.
  - `log_file`: `str`
    > a Gaussian log file to pull from
  - `:returns`: `_`
    >


<a id="Psience.Data.Surfaces.PotentialSurface.get_fchk_values" class="docs-object-method">&nbsp;</a> 
```python
get_fchk_values(fchk_file): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Data/Surfaces/PotentialSurface.py#L242)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Data/Surfaces/PotentialSurface.py#L242?message=Update%20Docs)]
</div>


<a id="Psience.Data.Surfaces.PotentialSurface.from_fchk_file" class="docs-object-method">&nbsp;</a> 
```python
from_fchk_file(fchk_file, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Data/Surfaces/PotentialSurface.py#L255)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Data/Surfaces/PotentialSurface.py#L255?message=Update%20Docs)]
</div>
Loads potential from a Gaussian formatted checkpoint file and builds a potential surface via a quartic approximation
  - `fchk_file`: `Any`
    > a Gaussian fchk file to pull from
  - `log_file`: `str`
    > 
  - `:returns`: `_`
    >
 </div>
</div>




## Examples













<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#Tests-48be42" markdown="1"> Tests</a> <a class="float-right" data-toggle="collapse" href="#Tests-48be42"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Tests-48be42" markdown="1">
 - [LogFilePotentialSurface](#LogFilePotentialSurface)

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
### <a class="collapse-link" data-toggle="collapse" href="#Setup-10193d" markdown="1"> Setup</a> <a class="float-right" data-toggle="collapse" href="#Setup-10193d"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Setup-10193d" markdown="1">
 
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Data/Surfaces/PotentialSurface.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Data/Surfaces/PotentialSurface.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Data/Surfaces/PotentialSurface.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Data/Surfaces/PotentialSurface.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Data/Surfaces.py#L167?message=Update%20Docs)   
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