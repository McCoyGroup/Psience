## <a id="Psience.Molecools.Properties.PotentialSurfaceManager">PotentialSurfaceManager</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties.py#L1878)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties.py#L1878?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
name: str
```
<a id="Psience.Molecools.Properties.PotentialSurfaceManager.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, mol, surface=None, derivatives=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties.py#L1880)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties.py#L1880?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.from_data" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_data(cls, mol, data): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1892)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1892?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.surface" class="docs-object-method">&nbsp;</a> 
```python
@property
surface(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L1896)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L1896?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.surface_coords" class="docs-object-method">&nbsp;</a> 
```python
@property
surface_coords(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L1901)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L1901?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.get_derivs" class="docs-object-method">&nbsp;</a> 
```python
get_derivs(self, quiet=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L1908)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L1908?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.derivatives" class="docs-object-method">&nbsp;</a> 
```python
@property
derivatives(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L1912)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L1912?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.force_constants" class="docs-object-method">&nbsp;</a> 
```python
@property
force_constants(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L1919)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L1919?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.load_potential_derivatives" class="docs-object-method">&nbsp;</a> 
```python
load_potential_derivatives(self, file=None, quiet=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L1923)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L1923?message=Update%20Docs)]
</div>
Loads potential derivatives from a file (or from `source_file` if set)
  - `file`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.load_potential_surface" class="docs-object-method">&nbsp;</a> 
```python
load_potential_surface(self, coordinates): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L1975)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L1975?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.load" class="docs-object-method">&nbsp;</a> 
```python
load(self, coordinates=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2004)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2004?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.update" class="docs-object-method">&nbsp;</a> 
```python
update(self, val): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2011)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2011?message=Update%20Docs)]
</div>
Updates the held values
  - `val`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.apply_transformation" class="docs-object-method">&nbsp;</a> 
```python
apply_transformation(self, transf): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2022)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2022?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.insert_atoms" class="docs-object-method">&nbsp;</a> 
```python
insert_atoms(self, atoms, coords, where): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2040)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2040?message=Update%20Docs)]
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


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.delete_atoms" class="docs-object-method">&nbsp;</a> 
```python
delete_atoms(self, where): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2059)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2059?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Molecools/Properties/PotentialSurfaceManager.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Molecools/Properties/PotentialSurfaceManager.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Molecools/Properties/PotentialSurfaceManager.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Molecools/Properties/PotentialSurfaceManager.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties.py#L1878?message=Update%20Docs)   
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