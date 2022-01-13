## <a id="Psience.BasisReps.Representations.ExpansionRepresentation">ExpansionRepresentation</a>
Provides support for terms that look like `1/2 pGp + 1/2 dV/dQdQ QQ` by computing each term on its own

### Properties and Methods
<a id="Psience.BasisReps.Representations.ExpansionRepresentation.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, coeffs, computers, basis, name=None, logger=None, memory_constrained=False): 
```

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

<a id="Psience.BasisReps.Representations.ExpansionRepresentation.__rmul__" class="docs-object-method">&nbsp;</a>
```python
__rmul__(self, other): 
```

<a id="Psience.BasisReps.Representations.ExpansionRepresentation.__mul__" class="docs-object-method">&nbsp;</a>
```python
__mul__(self, other): 
```

<a id="Psience.BasisReps.Representations.ExpansionRepresentation.__add__" class="docs-object-method">&nbsp;</a>
```python
__add__(self, other): 
```

<a id="Psience.BasisReps.Representations.ExpansionRepresentation.get_brakets" class="docs-object-method">&nbsp;</a>
```python
get_brakets(self, states, check_orthogonality=True): 
```

<a id="Psience.BasisReps.Representations.ExpansionRepresentation.get_element" class="docs-object-method">&nbsp;</a>
```python
get_element(self, n, m): 
```

<a id="Psience.BasisReps.Representations.ExpansionRepresentation.selection_rules" class="docs-object-method">&nbsp;</a>
```python
@property
selection_rules(self): 
```

- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Representations.ExpansionRepresentation.selection_rule_steps" class="docs-object-method">&nbsp;</a>
```python
@property
selection_rule_steps(self): 
```

- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Representations.ExpansionRepresentation.get_transformed_space" class="docs-object-method">&nbsp;</a>
```python
get_transformed_space(self, space, rules=None, parallelizer=None, logger=None, **opts): 
```
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

### Examples




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/edit/ci/examples/ci/docs/Psience/BasisReps/Representations/ExpansionRepresentation.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/examples/ci/docs/Psience/BasisReps/Representations/ExpansionRepresentation.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/edit/ci/docs/ci/docs/Psience/BasisReps/Representations/ExpansionRepresentation.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/docs/templates/ci/docs/Psience/BasisReps/Representations/ExpansionRepresentation.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py?message=Update%20Docs)