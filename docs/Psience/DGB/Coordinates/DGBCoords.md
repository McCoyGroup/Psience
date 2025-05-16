## <a id="Psience.DGB.Coordinates.DGBCoords">DGBCoords</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Coordinates.py#L19)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates.py#L19?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
DGBEmbeddedFunction: DGBEmbeddedFunction
```
<a id="Psience.DGB.Coordinates.DGBCoords.centers" class="docs-object-method">&nbsp;</a> 
```python
@property
centers(self) -> 'np.ndarray': 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Coordinates/DGBCoords.py#L21)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates/DGBCoords.py#L21?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBCoords.shape" class="docs-object-method">&nbsp;</a> 
```python
@property
shape(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Coordinates/DGBCoords.py#L25)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates/DGBCoords.py#L25?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBCoords.kinetic_energy_evaluator" class="docs-object-method">&nbsp;</a> 
```python
@property
kinetic_energy_evaluator(self) -> 'DGBKineticEnergyEvaluator': 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Coordinates/DGBCoords.py#L29)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates/DGBCoords.py#L29?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBCoords.pairwise_potential_evaluator_type" class="docs-object-method">&nbsp;</a> 
```python
@property
pairwise_potential_evaluator_type(self) -> 'type[DGBPairwisePotentialEvaluator]': 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Coordinates/DGBCoords.py#L34)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates/DGBCoords.py#L34?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBCoords.pairwise_potential_evaluator" class="docs-object-method">&nbsp;</a> 
```python
pairwise_potential_evaluator(self, potential_functions) -> 'DGBPairwisePotentialEvaluator': 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Coordinates/DGBCoords.py#L38)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates/DGBCoords.py#L38?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBCoords.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item) -> 'DGBCoords': 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Coordinates/DGBCoords.py#L48)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates/DGBCoords.py#L48?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBCoords.take_indices" class="docs-object-method">&nbsp;</a> 
```python
take_indices(self, subinds) -> 'DGBCoords': 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Coordinates/DGBCoords.py#L51)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates/DGBCoords.py#L51?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBCoords.drop_indices" class="docs-object-method">&nbsp;</a> 
```python
drop_indices(self, subinds) -> 'DGBCoords': 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Coordinates/DGBCoords.py#L53)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates/DGBCoords.py#L53?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBCoords.gmatrix" class="docs-object-method">&nbsp;</a> 
```python
gmatrix(self, coords: numpy.ndarray) -> numpy.ndarray: 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Coordinates/DGBCoords.py#L56)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates/DGBCoords.py#L56?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBCoords.embedded_mode_function" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
embedded_mode_function(cls, func, modes, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L60)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L60?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBCoords.embedded_subcoordinate_function" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
embedded_subcoordinate_function(cls, func, sel, ndim): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L79)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L79?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBCoords.embedded_cartesian_function" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
embedded_cartesian_function(cls, func, atom_sel, xyz_sel, natoms, ndim): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L104)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L104?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBCoords.embed_function" class="docs-object-method">&nbsp;</a> 
```python
embed_function(self, fn) -> 'DGBEmbeddedFunction': 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Coordinates/DGBCoords.py#L191)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates/DGBCoords.py#L191?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Coordinates.DGBCoords.as_cartesians" class="docs-object-method">&nbsp;</a> 
```python
as_cartesians(self) -> 'tuple[DGBCartesians, tuple[np.ndarray, np.ndarray]]': 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Coordinates/DGBCoords.py#L195)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates/DGBCoords.py#L195?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/DGB/Coordinates/DGBCoords.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/DGB/Coordinates/DGBCoords.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/DGB/Coordinates/DGBCoords.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/DGB/Coordinates/DGBCoords.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Coordinates.py#L19?message=Update%20Docs)   
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