## <a id="Psience.Molecools.CoordinateSystems.ModeEmbedding">ModeEmbedding</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/CoordinateSystems.py#L823)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/CoordinateSystems.py#L823?message=Update%20Docs)]
</div>

Provides a specialization on a `MoleculaEmbedding` to express all properties
in terms of the attendant normal modes







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.Molecools.CoordinateSystems.ModeEmbedding.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, embedding: Psience.Molecools.CoordinateSystems.MolecularEmbedding, modes, mass_weight=None, dimensionless=None, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/CoordinateSystems/ModeEmbedding.py#L828)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/CoordinateSystems/ModeEmbedding.py#L828?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.ModeEmbedding.mw_conversion" class="docs-object-method">&nbsp;</a> 
```python
mw_conversion(self, strip_dummies=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/CoordinateSystems/ModeEmbedding.py#L853)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/CoordinateSystems/ModeEmbedding.py#L853?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.ModeEmbedding.mw_inverse" class="docs-object-method">&nbsp;</a> 
```python
mw_inverse(self, strip_dummies=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/CoordinateSystems/ModeEmbedding.py#L864)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/CoordinateSystems/ModeEmbedding.py#L864?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.ModeEmbedding.get_mw_cartesians_by_internals" class="docs-object-method">&nbsp;</a> 
```python
get_mw_cartesians_by_internals(self, order=None, strip_embedding=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/CoordinateSystems/ModeEmbedding.py#L876)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/CoordinateSystems/ModeEmbedding.py#L876?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.ModeEmbedding.get_internals_by_mw_cartesians" class="docs-object-method">&nbsp;</a> 
```python
get_internals_by_mw_cartesians(self, order=None, strip_embedding=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/CoordinateSystems/ModeEmbedding.py#L888)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/CoordinateSystems/ModeEmbedding.py#L888?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.ModeEmbedding.get_internals_by_cartesians" class="docs-object-method">&nbsp;</a> 
```python
get_internals_by_cartesians(self, order=None, strip_embedding=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/CoordinateSystems/ModeEmbedding.py#L904)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/CoordinateSystems/ModeEmbedding.py#L904?message=Update%20Docs)]
</div>
expresses raw internals or modes (internals or Cartesian) in terms of mass-weighted Cartesians
  - `order`: `Any`
    > 
  - `strip_embedding`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.CoordinateSystems.ModeEmbedding.get_cartesians_by_internals" class="docs-object-method">&nbsp;</a> 
```python
get_cartesians_by_internals(self, order=None, strip_embedding=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/CoordinateSystems/ModeEmbedding.py#L951)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/CoordinateSystems/ModeEmbedding.py#L951?message=Update%20Docs)]
</div>
expresses raw internals or modes (internals or Cartesian) in terms of mass-weighted Cartesians
  - `order`: `Any`
    > 
  - `strip_embedding`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.CoordinateSystems.ModeEmbedding.get_inertia_tensor_expansion" class="docs-object-method">&nbsp;</a> 
```python
get_inertia_tensor_expansion(self, order=None, strip_embedding=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/CoordinateSystems/ModeEmbedding.py#L999)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/CoordinateSystems/ModeEmbedding.py#L999?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.ModeEmbedding.get_inertial_frame" class="docs-object-method">&nbsp;</a> 
```python
get_inertial_frame(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/CoordinateSystems/ModeEmbedding.py#L1004)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/CoordinateSystems/ModeEmbedding.py#L1004?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Molecools/CoordinateSystems/ModeEmbedding.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Molecools/CoordinateSystems/ModeEmbedding.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Molecools/CoordinateSystems/ModeEmbedding.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Molecools/CoordinateSystems/ModeEmbedding.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/CoordinateSystems.py#L823?message=Update%20Docs)   
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