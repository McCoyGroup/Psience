## <a id="Psience.BasisReps.Terms.Representation">Representation</a>
A `Representation` provides a simple interface to compute only some elements of high-dimensional tensors.
It takes a tensor shape and a function to compute tensor elements.
The `compute` function should be able to take a block of indices and return all the matrix elements.

### Properties and Methods
<a id="Psience.BasisReps.Terms.Representation.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, compute, basis, logger=None): 
```

- `compute`: `callable | Operator`
    >the function that turns indices into values
- `n_quanta`: `RepresentationBasis`
    >the basis quanta used in the representations (necessary for shape reasons)
- `logger`: `None | Logger`
    >logger for printing out debug info

<a id="Psience.BasisReps.Terms.Representation.clear_cache" class="docs-object-method">&nbsp;</a>
```python
clear_cache(self): 
```

<a id="Psience.BasisReps.Terms.Representation.diag" class="docs-object-method">&nbsp;</a>
```python
@property
diag(self): 
```

<a id="Psience.BasisReps.Terms.Representation.ndims" class="docs-object-method">&nbsp;</a>
```python
@property
ndims(self): 
```

<a id="Psience.BasisReps.Terms.Representation.dim_inds" class="docs-object-method">&nbsp;</a>
```python
@property
dim_inds(self): 
```

<a id="Psience.BasisReps.Terms.Representation.get_brakets" class="docs-object-method">&nbsp;</a>
```python
get_brakets(self, states): 
```
Computes term elements based on getting a BraKetSpace.
        Can directly pass element specs through, since the shape management shit
        is managed by the BraKetSpace
- `states`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Terms.Representation.get_element" class="docs-object-method">&nbsp;</a>
```python
get_element(self, n, m): 
```
Computes term elements.
        Determines first whether it needs to pull single elements or blocks of them.
- `n`: `Any`
    >No description...
- `m`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Terms.Representation.__getitem__" class="docs-object-method">&nbsp;</a>
```python
__getitem__(self, item): 
```

<a id="Psience.BasisReps.Terms.Representation.__rmul__" class="docs-object-method">&nbsp;</a>
```python
__rmul__(self, other): 
```

<a id="Psience.BasisReps.Terms.Representation.__mul__" class="docs-object-method">&nbsp;</a>
```python
__mul__(self, other): 
```

<a id="Psience.BasisReps.Terms.Representation.__add__" class="docs-object-method">&nbsp;</a>
```python
__add__(self, other): 
```

<a id="Psience.BasisReps.Terms.Representation.__repr__" class="docs-object-method">&nbsp;</a>
```python
__repr__(self): 
```

### Examples


