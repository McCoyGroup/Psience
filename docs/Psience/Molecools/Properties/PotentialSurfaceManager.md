## <a id="Psience.Molecools.Properties.PotentialSurfaceManager">PotentialSurfaceManager</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties.py#L1968)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties.py#L1968?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties.py#L1970)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties.py#L1970?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.from_data" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_data(cls, mol, data): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1982)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1982?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.surface" class="docs-object-method">&nbsp;</a> 
```python
@property
surface(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L1986)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L1986?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.surface_coords" class="docs-object-method">&nbsp;</a> 
```python
@property
surface_coords(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L1991)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L1991?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.get_derivs" class="docs-object-method">&nbsp;</a> 
```python
get_derivs(self, quiet=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L1998)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L1998?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.derivatives" class="docs-object-method">&nbsp;</a> 
```python
@property
derivatives(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2002)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2002?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.force_constants" class="docs-object-method">&nbsp;</a> 
```python
@property
force_constants(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2009)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2009?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.load_potential_derivatives" class="docs-object-method">&nbsp;</a> 
```python
load_potential_derivatives(self, file=None, mode=None, quiet=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2013)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2013?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2111)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2111?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.load" class="docs-object-method">&nbsp;</a> 
```python
load(self, coordinates=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2140)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2140?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.update" class="docs-object-method">&nbsp;</a> 
```python
update(self, val): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2147)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2147?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2158)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2158?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.insert_atoms" class="docs-object-method">&nbsp;</a> 
```python
insert_atoms(self, atoms, coords, where): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2176)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2176?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2195)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2195?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties.py#L1968?message=Update%20Docs)   
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