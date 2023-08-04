## <a id="Psience.BasisReps.Representations.ExpansionRepresentation">ExpansionRepresentation</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations.py#L906)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations.py#L906?message=Update%20Docs)]
</div>

Provides support for terms that look like `1/2 pGp + 1/2 dV/dQdQ QQ` by computing each term on its own







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.BasisReps.Representations.ExpansionRepresentation.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, coeffs, computers, basis, name=None, logger=None, memory_constrained=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/ExpansionRepresentation.py#L910)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/ExpansionRepresentation.py#L910?message=Update%20Docs)]
</div>

  - `coeffs`: `Iterable[float]`
    > The expansion coefficients
  - `compute`: `Iterable[callable | Operator]`
    > the functions that turns indices into values
  - `n_quanta`: `tuple[int]`
    > the total quanta used in the representations (necessary for shape reasons)


<a id="Psience.BasisReps.Representations.ExpansionRepresentation.is_diagonal" class="docs-object-method">&nbsp;</a> 
```python
@property
is_diagonal(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/ExpansionRepresentation.py#L923)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/ExpansionRepresentation.py#L923?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Representations.ExpansionRepresentation.is_zero" class="docs-object-method">&nbsp;</a> 
```python
@property
is_zero(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/ExpansionRepresentation.py#L929)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/ExpansionRepresentation.py#L929?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Representations.ExpansionRepresentation.skipped_indices" class="docs-object-method">&nbsp;</a> 
```python
@property
skipped_indices(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/ExpansionRepresentation.py#L936)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/ExpansionRepresentation.py#L936?message=Update%20Docs)]
</div>

  - `:returns`: `_`
    >


<a id="Psience.BasisReps.Representations.ExpansionRepresentation.clear_cache" class="docs-object-method">&nbsp;</a> 
```python
clear_cache(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/ExpansionRepresentation.py#L961)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/ExpansionRepresentation.py#L961?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Representations.ExpansionRepresentation.__rmul__" class="docs-object-method">&nbsp;</a> 
```python
__rmul__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/ExpansionRepresentation.py#L965)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/ExpansionRepresentation.py#L965?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Representations.ExpansionRepresentation.__mul__" class="docs-object-method">&nbsp;</a> 
```python
__mul__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/ExpansionRepresentation.py#L976)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/ExpansionRepresentation.py#L976?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Representations.ExpansionRepresentation.__add__" class="docs-object-method">&nbsp;</a> 
```python
__add__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/ExpansionRepresentation.py#L988)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/ExpansionRepresentation.py#L988?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Representations.ExpansionRepresentation.get_brakets" class="docs-object-method">&nbsp;</a> 
```python
get_brakets(self, states, check_orthogonality=True, memory_constrained=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/ExpansionRepresentation.py#L1066)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/ExpansionRepresentation.py#L1066?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Representations.ExpansionRepresentation.get_element" class="docs-object-method">&nbsp;</a> 
```python
get_element(self, n, m): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/ExpansionRepresentation.py#L1070)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/ExpansionRepresentation.py#L1070?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Representations.ExpansionRepresentation.chunk_size" class="docs-object-method">&nbsp;</a> 
```python
@property
chunk_size(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/ExpansionRepresentation.py#L1072)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/ExpansionRepresentation.py#L1072?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Representations.ExpansionRepresentation.selection_rules" class="docs-object-method">&nbsp;</a> 
```python
@property
selection_rules(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/ExpansionRepresentation.py#L1080)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/ExpansionRepresentation.py#L1080?message=Update%20Docs)]
</div>

  - `:returns`: `_`
    >


<a id="Psience.BasisReps.Representations.ExpansionRepresentation.selection_rule_steps" class="docs-object-method">&nbsp;</a> 
```python
@property
selection_rule_steps(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/ExpansionRepresentation.py#L1097)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/ExpansionRepresentation.py#L1097?message=Update%20Docs)]
</div>

  - `:returns`: `_`
    >


<a id="Psience.BasisReps.Representations.ExpansionRepresentation.get_transformed_space" class="docs-object-method">&nbsp;</a> 
```python
get_transformed_space(self, space, rules=None, parallelizer=None, logger=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/ExpansionRepresentation.py#L1114)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/ExpansionRepresentation.py#L1114?message=Update%20Docs)]
</div>
Returns the state space obtained by using the
held operators to transform `space`
  - `space`: `BasisStateSpace`
    > 
  - `:returns`: `SelectionRuleStateSpace`
    >


<a id="Psience.BasisReps.Representations.ExpansionRepresentation.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Representations/ExpansionRepresentation.py#L1149)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations/ExpansionRepresentation.py#L1149?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/Representations/ExpansionRepresentation.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/Representations/ExpansionRepresentation.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/Representations/ExpansionRepresentation.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/Representations/ExpansionRepresentation.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Representations.py#L906?message=Update%20Docs)   
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