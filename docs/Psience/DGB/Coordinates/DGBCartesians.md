## <a id="Psience.DGB.Coordinates.DGBCartesians">DGBCartesians</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Coordinates.py#L200)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates.py#L200?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.DGB.Coordinates.DGBCartesians.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, coords, masses, *, natoms=None, atom_sel=None, ndim=None, xyz_sel=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Coordinates.py#L201)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates.py#L201?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBCartesians.centers" class="docs-object-method">&nbsp;</a> 
```python
@property
centers(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Coordinates/DGBCartesians.py#L214)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates/DGBCartesians.py#L214?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBCartesians.cart_shape" class="docs-object-method">&nbsp;</a> 
```python
@property
cart_shape(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Coordinates/DGBCartesians.py#L217)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates/DGBCartesians.py#L217?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBCartesians.kinetic_energy_evaluator" class="docs-object-method">&nbsp;</a> 
```python
@property
kinetic_energy_evaluator(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Coordinates/DGBCartesians.py#L220)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates/DGBCartesians.py#L220?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBCartesians.pairwise_potential_evaluator_type" class="docs-object-method">&nbsp;</a> 
```python
@property
pairwise_potential_evaluator_type(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Coordinates/DGBCartesians.py#L225)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates/DGBCartesians.py#L225?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBCartesians.resolve_masses" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
resolve_masses(cls, coords, masses=None, atoms=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L228)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L228?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBCartesians.from_cartesians" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_cartesians(cls, centers, masses=None, atoms=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L245)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L245?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBCartesians.infer_shape_sel" class="docs-object-method">&nbsp;</a> 
```python
infer_shape_sel(self, selector): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Coordinates/DGBCartesians.py#L251)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates/DGBCartesians.py#L251?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBCartesians.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item) -> 'DGBCartesians': 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Coordinates/DGBCartesians.py#L275)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates/DGBCartesians.py#L275?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBCartesians.take_indices" class="docs-object-method">&nbsp;</a> 
```python
take_indices(self, subinds): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Coordinates/DGBCartesians.py#L302)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates/DGBCartesians.py#L302?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBCartesians.embed_function" class="docs-object-method">&nbsp;</a> 
```python
embed_function(self, function): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Coordinates/DGBCartesians.py#L307)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates/DGBCartesians.py#L307?message=Update%20Docs)]
</div>
Embeds assuming we got a function in Cartesians _before_ any selections happened
  - `function`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.DGB.Coordinates.DGBCartesians.gmatrix" class="docs-object-method">&nbsp;</a> 
```python
gmatrix(self, coords: numpy.ndarray) -> numpy.ndarray: 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Coordinates/DGBCartesians.py#L326)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates/DGBCartesians.py#L326?message=Update%20Docs)]
</div>
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/DGB/Coordinates/DGBCartesians.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/DGB/Coordinates/DGBCartesians.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/DGB/Coordinates/DGBCartesians.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/DGB/Coordinates/DGBCartesians.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates.py#L200?message=Update%20Docs)   
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