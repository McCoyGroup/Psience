## <a id="Psience.Data.Surfaces.DipoleSurface">DipoleSurface</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Data/Surfaces.py#L16)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Data/Surfaces.py#L16?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Data/Surfaces/DipoleSurface.py#L22)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Data/Surfaces/DipoleSurface.py#L22?message=Update%20Docs)]
</div>

  - `mu_x`: `Surface`
    > X-component of dipole moment
  - `mu_y`: `Surface`
    > Y-component of dipole moment
  - `mu_z`: `Surface`
    > Z-component of dipole moment


<a id="Psience.Data.Surfaces.DipoleSurface.center" class="docs-object-method">&nbsp;</a> 
```python
@property
center(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Data/Surfaces/DipoleSurface.py#L42)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Data/Surfaces/DipoleSurface.py#L42?message=Update%20Docs)]
</div>


<a id="Psience.Data.Surfaces.DipoleSurface.ref" class="docs-object-method">&nbsp;</a> 
```python
@property
ref(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Data/Surfaces/DipoleSurface.py#L45)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Data/Surfaces/DipoleSurface.py#L45?message=Update%20Docs)]
</div>


<a id="Psience.Data.Surfaces.DipoleSurface.expansion_tensors" class="docs-object-method">&nbsp;</a> 
```python
@property
expansion_tensors(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Data/Surfaces/DipoleSurface.py#L48)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Data/Surfaces/DipoleSurface.py#L48?message=Update%20Docs)]
</div>


<a id="Psience.Data.Surfaces.DipoleSurface.get_log_values" class="docs-object-method">&nbsp;</a> 
```python
@staticmethod
get_log_values(log_file, keys=('StandardCartesianCoordinates', 'DipoleMoments')): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L57)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L57?message=Update%20Docs)]
</div>


<a id="Psience.Data.Surfaces.DipoleSurface.from_log_file" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_log_file(cls, log_file, coord_transf, keys=('StandardCartesianCoordinates', 'DipoleMoments'), tol=0.001, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L67)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L67?message=Update%20Docs)]
</div>
Loads dipoles from a Gaussian log file and builds a dipole surface by interpolating.
Obviously this only really works if we have a subset of "scan" coordinates, so at this stage the user is obligated
to furnish a function that'll take a set of Cartesian coordinates and convert them to "scan" coordinates.
Coordinerds can be helpful with this, as it provides a convenient syntax for Cartesian <-> ZMatrix conversions
  - `log_file`: `str`
    > a Gaussian log file to pull from
  - `:returns`: `_`
    >


<a id="Psience.Data.Surfaces.DipoleSurface.from_fchk_file" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_fchk_file(cls, fchk_file, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L134)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L134?message=Update%20Docs)]
</div>
Loads dipoles from a Gaussian formatted checkpoint file and builds a dipole surface via a linear approximation
  - `fchk_file`: `Any`
    > a Gaussian fchk file to pull from
  - `log_file`: `str`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Data.Surfaces.DipoleSurface.from_derivatives" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_derivatives(cls, expansion, center=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L151)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L151?message=Update%20Docs)]
</div>


<a id="Psience.Data.Surfaces.DipoleSurface.from_mol" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_mol(cls, mol, expansion=None, center=None, transforms=None, use_internals=True, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L175)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L175?message=Update%20Docs)]
</div>


<a id="Psience.Data.Surfaces.DipoleSurface.__call__" class="docs-object-method">&nbsp;</a> 
```python
__call__(self, gridpoints, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Data/Surfaces/DipoleSurface.py#L190)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Data/Surfaces/DipoleSurface.py#L190?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Data/Surfaces/DipoleSurface.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Data/Surfaces/DipoleSurface.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Data/Surfaces/DipoleSurface.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Data/Surfaces/DipoleSurface.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Data/Surfaces.py#L16?message=Update%20Docs)   
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