## <a id="Psience.BasisReps.Terms.ExpansionRepresentation">ExpansionRepresentation</a>
Provides support for terms that look like `1/2 pGp + 1/2 dV/dQdQ QQ` by computing each term on its own

### Properties and Methods
<a id="Psience.BasisReps.Terms.ExpansionRepresentation.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, coeffs, computers, basis, logger=None): 
```

- `coeffs`: `Iterable[float]`
    >The expansion coefficients
- `compute`: `Iterable[callable | Operator]`
    >the functions that turns indices into values
- `n_quanta`: `tuple[int]`
    >the total quanta used in the representations (necessary for shape reasons)

<a id="Psience.BasisReps.Terms.ExpansionRepresentation.clear_cache" class="docs-object-method">&nbsp;</a>
```python
clear_cache(self): 
```

<a id="Psience.BasisReps.Terms.ExpansionRepresentation.__rmul__" class="docs-object-method">&nbsp;</a>
```python
__rmul__(self, other): 
```

<a id="Psience.BasisReps.Terms.ExpansionRepresentation.__mul__" class="docs-object-method">&nbsp;</a>
```python
__mul__(self, other): 
```

<a id="Psience.BasisReps.Terms.ExpansionRepresentation.__add__" class="docs-object-method">&nbsp;</a>
```python
__add__(self, other): 
```

<a id="Psience.BasisReps.Terms.ExpansionRepresentation.get_brakets" class="docs-object-method">&nbsp;</a>
```python
get_brakets(self, states): 
```

<a id="Psience.BasisReps.Terms.ExpansionRepresentation.get_element" class="docs-object-method">&nbsp;</a>
```python
get_element(self, n, m): 
```

<a id="Psience.BasisReps.Terms.ExpansionRepresentation.__repr__" class="docs-object-method">&nbsp;</a>
```python
__repr__(self): 
```

### Examples


