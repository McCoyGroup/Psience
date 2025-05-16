## <a id="Psience.DGB.Coordinates.DGBWatsonModes">DGBWatsonModes</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Coordinates.py#L333)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates.py#L333?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.DGB.Coordinates.DGBWatsonModes.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, coords, modes, *, coriolis_inertia_function=None, masses=None, subselection=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Coordinates.py#L334)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates.py#L334?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBWatsonModes.centers" class="docs-object-method">&nbsp;</a> 
```python
@property
centers(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Coordinates/DGBWatsonModes.py#L346)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates/DGBWatsonModes.py#L346?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBWatsonModes.kinetic_energy_evaluator" class="docs-object-method">&nbsp;</a> 
```python
@property
kinetic_energy_evaluator(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Coordinates/DGBWatsonModes.py#L349)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates/DGBWatsonModes.py#L349?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBWatsonModes.pairwise_potential_evaluator_type" class="docs-object-method">&nbsp;</a> 
```python
@property
pairwise_potential_evaluator_type(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Coordinates/DGBWatsonModes.py#L352)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates/DGBWatsonModes.py#L352?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBWatsonModes.zeta_momi" class="docs-object-method">&nbsp;</a> 
```python
@staticmethod
zeta_momi(watson_coords, modes, masses): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/staticmethod.py#L355)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/staticmethod.py#L355?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBWatsonModes.default_coriolis_inertia_function" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
default_coriolis_inertia_function(cls, modes, masses): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L377)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L377?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBWatsonModes.embed_coords" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
embed_coords(cls, carts, modes, shift=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L398)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L398?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBWatsonModes.unembed_coords" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
unembed_coords(cls, mode_coords, modes, masses=None, shift=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L409)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L409?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBWatsonModes.embed_derivs" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
embed_derivs(cls, derivs, modes): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L425)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L425?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBWatsonModes.from_cartesians" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_cartesians(cls, coords, modes, masses=None, coriolis_inertia_function=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L434)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L434?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBWatsonModes.as_cartesians" class="docs-object-method">&nbsp;</a> 
```python
as_cartesians(self, masses=None) -> 'tuple[DGBCartesians, tuple[np.ndarray, np.ndarray]]': 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Coordinates/DGBWatsonModes.py#L449)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates/DGBWatsonModes.py#L449?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBWatsonModes.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Coordinates/DGBWatsonModes.py#L462)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates/DGBWatsonModes.py#L462?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBWatsonModes.gmatrix" class="docs-object-method">&nbsp;</a> 
```python
gmatrix(self, coords: numpy.ndarray) -> numpy.ndarray: 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Coordinates/DGBWatsonModes.py#L492)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates/DGBWatsonModes.py#L492?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBWatsonModes.embed_function" class="docs-object-method">&nbsp;</a> 
```python
embed_function(self, fn): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Coordinates/DGBWatsonModes.py#L496)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates/DGBWatsonModes.py#L496?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/DGB/Coordinates/DGBWatsonModes.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/DGB/Coordinates/DGBWatsonModes.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/DGB/Coordinates/DGBWatsonModes.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/DGB/Coordinates/DGBWatsonModes.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates.py#L333?message=Update%20Docs)   
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