## <a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticKineticEnergyConstructor">AnalyticKineticEnergyConstructor</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors.py#L266)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors.py#L266?message=Update%20Docs)]
</div>

Provides G and V' elements from Frederick and Woywood







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticKineticEnergyConstructor.kinetic_exprs" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
kinetic_exprs(cls, inds1: 'Iterable[int]', inds2: 'Iterable[int]', coord_types=None, target_symbols=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L364)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L364?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticKineticEnergyConstructor.kinetic_exprs_direct" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
kinetic_exprs_direct(cls, inds1: 'Iterable[int]', inds2: 'Iterable[int]', coord_types=None, return_vp=True, target_symbols=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L487)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L487?message=Update%20Docs)]
</div>
Evaluated using the simple expressions in Table 1 from Frederick and Woywood
  - `inds1`: `Any`
    > 
  - `inds2`: `Any`
    > 
  - `coord_types`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticKineticEnergyConstructor.g" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
g(cls, inds1: 'Iterable[int]', inds2: 'Iterable[int]', coord_types=None, target_symbols=None, return_function=False, method='lookup'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L529)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L529?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticKineticEnergyConstructor.vp" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
vp(cls, inds1: 'Iterable[int]', inds2: 'Iterable[int]', coord_types=None, target_symbols=None, return_function=False, method='lookup'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L539)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L539?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.AnalyticKineticEnergyConstructor.g_matrix" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
g_matrix(cls, coord_specs, return_function=False, return_matrix=True, method='lookup'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L550)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L550?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/AnalyticModels/AnalyticModelConstructors/AnalyticKineticEnergyConstructor.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/AnalyticModels/AnalyticModelConstructors/AnalyticKineticEnergyConstructor.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/AnalyticModels/AnalyticModelConstructors/AnalyticKineticEnergyConstructor.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/AnalyticModels/AnalyticModelConstructors/AnalyticKineticEnergyConstructor.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors.py#L266?message=Update%20Docs)   
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