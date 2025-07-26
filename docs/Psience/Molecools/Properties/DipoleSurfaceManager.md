## <a id="Psience.Molecools.Properties.DipoleSurfaceManager">DipoleSurfaceManager</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties.py#L1644)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties.py#L1644?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
name: str
```
<a id="Psience.Molecools.Properties.DipoleSurfaceManager.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, mol, surface=None, derivatives=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties.py#L1646)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties.py#L1646?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.DipoleSurfaceManager.from_data" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_data(cls, mol, data): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1661)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1661?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.DipoleSurfaceManager.surface" class="docs-object-method">&nbsp;</a> 
```python
@property
surface(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L1665)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L1665?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.DipoleSurfaceManager.numerical_derivatives" class="docs-object-method">&nbsp;</a> 
```python
@property
numerical_derivatives(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L1670)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L1670?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.DipoleSurfaceManager.get_derivatives" class="docs-object-method">&nbsp;</a> 
```python
get_derivatives(self, quiet=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L1681)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L1681?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.DipoleSurfaceManager.derivatives" class="docs-object-method">&nbsp;</a> 
```python
@property
derivatives(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L1691)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L1691?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.DipoleSurfaceManager.load" class="docs-object-method">&nbsp;</a> 
```python
load(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L1703)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L1703?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.DipoleSurfaceManager.update" class="docs-object-method">&nbsp;</a> 
```python
update(self, val): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L1708)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L1708?message=Update%20Docs)]
</div>
Updates the held values
  - `val`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.DipoleSurfaceManager.load_dipole_surface" class="docs-object-method">&nbsp;</a> 
```python
load_dipole_surface(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L1719)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L1719?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.DipoleSurfaceManager.load_dipole_derivatives" class="docs-object-method">&nbsp;</a> 
```python
load_dipole_derivatives(self, file=None, quiet=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L1778)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L1778?message=Update%20Docs)]
</div>
Loads dipole derivatives from a file (or from `source_file` if set)
  - `file`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.DipoleSurfaceManager.apply_transformation" class="docs-object-method">&nbsp;</a> 
```python
apply_transformation(self, transf): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L1823)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L1823?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.DipoleSurfaceManager.insert_atoms" class="docs-object-method">&nbsp;</a> 
```python
insert_atoms(self, atoms, coords, where): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L1850)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L1850?message=Update%20Docs)]
</div>
Handles the insertion of new atoms into the structure
  - `atoms`: `tuple[str]`
    > 
  - `coords`: `CoordinateSet`
    > 
  - `where`: `tuple[int]`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.DipoleSurfaceManager.delete_atoms" class="docs-object-method">&nbsp;</a> 
```python
delete_atoms(self, where): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L1873)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L1873?message=Update%20Docs)]
</div>
Handles the deletion from the structure
  - `atoms`: `tuple[str]`
    > 
  - `coords`: `CoordinateSet`
    > 
  - `where`: `tuple[int]`
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Molecools/Properties/DipoleSurfaceManager.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Molecools/Properties/DipoleSurfaceManager.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Molecools/Properties/DipoleSurfaceManager.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Molecools/Properties/DipoleSurfaceManager.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties.py#L1644?message=Update%20Docs)   
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