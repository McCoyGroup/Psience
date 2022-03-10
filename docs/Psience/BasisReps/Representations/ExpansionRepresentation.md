## <a id="Psience.BasisReps.Representations.ExpansionRepresentation">ExpansionRepresentation</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L871)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L871?message=Update%20Docs)]
</div>

Provides support for terms that look like `1/2 pGp + 1/2 dV/dQdQ QQ` by computing each term on its own

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

<a id="Psience.BasisReps.Representations.ExpansionRepresentation.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, coeffs, computers, basis, name=None, logger=None, memory_constrained=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L875)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L875?message=Update%20Docs)]
</div>


- `coeffs`: `Iterable[float]`
    >The expansion coefficients
- `compute`: `Iterable[callable | Operator]`
    >the functions that turns indices into values
- `n_quanta`: `tuple[int]`
    >the total quanta used in the representations (necessary for shape reasons)

<a id="Psience.BasisReps.Representations.ExpansionRepresentation.clear_cache" class="docs-object-method">&nbsp;</a> 
```python
clear_cache(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L888)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L888?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.ExpansionRepresentation.__rmul__" class="docs-object-method">&nbsp;</a> 
```python
__rmul__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L892)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L892?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.ExpansionRepresentation.__mul__" class="docs-object-method">&nbsp;</a> 
```python
__mul__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L903)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L903?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.ExpansionRepresentation.__add__" class="docs-object-method">&nbsp;</a> 
```python
__add__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L915)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L915?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.ExpansionRepresentation.get_brakets" class="docs-object-method">&nbsp;</a> 
```python
get_brakets(self, states, check_orthogonality=True, memory_constrained=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L993)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L993?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.ExpansionRepresentation.get_element" class="docs-object-method">&nbsp;</a> 
```python
get_element(self, n, m): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L997)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L997?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.ExpansionRepresentation.chunk_size" class="docs-object-method">&nbsp;</a> 
```python
@property
chunk_size(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.ExpansionRepresentation.selection_rules" class="docs-object-method">&nbsp;</a> 
```python
@property
selection_rules(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L?message=Update%20Docs)]
</div>


- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Representations.ExpansionRepresentation.selection_rule_steps" class="docs-object-method">&nbsp;</a> 
```python
@property
selection_rule_steps(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L?message=Update%20Docs)]
</div>


- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Representations.ExpansionRepresentation.get_transformed_space" class="docs-object-method">&nbsp;</a> 
```python
get_transformed_space(self, space, rules=None, parallelizer=None, logger=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L1041)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L1041?message=Update%20Docs)]
</div>

Returns the state space obtained by using the
        held operators to transform `space`
- `space`: `BasisStateSpace`
    >No description...
- `:returns`: `SelectionRuleStateSpace`
    >No description...

<a id="Psience.BasisReps.Representations.ExpansionRepresentation.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L1068)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L1068?message=Update%20Docs)]
</div>

 </div>
</div>




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/Representations/ExpansionRepresentation.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/Representations/ExpansionRepresentation.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/Representations/ExpansionRepresentation.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/Representations/ExpansionRepresentation.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L871?message=Update%20Docs)