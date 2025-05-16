## <a id="Psience.Data.Surfaces.PotentialSurface">PotentialSurface</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Data/Surfaces.py#L209)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/Surfaces.py#L209?message=Update%20Docs)]
</div>

A potential surface structure to go along with the DipoleSurface.
Provides convenient access to potential data + a unified interface to things like energy minimization







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.Data.Surfaces.PotentialSurface.get_log_values" class="docs-object-method">&nbsp;</a> 
```python
@staticmethod
get_log_values(log_file, keys=('StandardCartesianCoordinates', 'ScanEnergies')): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/staticmethod.py#L215)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/staticmethod.py#L215?message=Update%20Docs)]
</div>


<a id="Psience.Data.Surfaces.PotentialSurface.from_log_file" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_log_file(cls, log_file, coord_transf, keys=('StandardCartesianCoordinates', 'ScanEnergies'), tol=0.001, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L231)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L231?message=Update%20Docs)]
</div>
Loads dipoles from a Gaussian log file and builds a potential surface by interpolating.
Obviously this only really works if we have a subset of "scan" coordinates, so at this stage the user is obligated
to furnish a function that'll take a set of Cartesian coordinates and convert them to "scan" coordinates.
Coordinerds can be helpful with this, as it provides a convenient syntax for Cartesian <-> ZMatrix conversions.
  - `log_file`: `str`
    > a Gaussian log file to pull from
  - `:returns`: `_`
    >


<a id="Psience.Data.Surfaces.PotentialSurface.from_fchk_file" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_fchk_file(cls, fchk_file, ref=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L299)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L299?message=Update%20Docs)]
</div>
Loads potential from a Gaussian formatted checkpoint file and builds a potential surface via a quartic approximation
  - `fchk_file`: `Any`
    > a Gaussian fchk file to pull from
  - `log_file`: `str`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Data.Surfaces.PotentialSurface.from_mol" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_mol(cls, mol, expansion=None, center=None, transforms=None, transformed_derivatives=False, use_internals=True, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L322)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L322?message=Update%20Docs)]
</div>


<a id="Psience.Data.Surfaces.PotentialSurface.from_derivatives" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_derivatives(cls, expansion, center=None, ref=None, transforms=None, transformed_derivatives=False, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L344)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L344?message=Update%20Docs)]
</div>


<a id="Psience.Data.Surfaces.PotentialSurface.__call__" class="docs-object-method">&nbsp;</a> 
```python
__call__(self, gridpoints, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Data/Surfaces/PotentialSurface.py#L360)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/Surfaces/PotentialSurface.py#L360?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/Surfaces.py#L209?message=Update%20Docs)   
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