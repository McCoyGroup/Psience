## <a id="Psience.BasisReps.Terms.Representation">Representation</a>
A `Representation` provides a simple interface to build matrix representations of operators expressed
in high-dimensional spaces.

### Properties and Methods
<a id="Psience.BasisReps.Terms.Representation.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, compute, basis, name=None, logger=None, selection_rules=None): 
```

- `compute`: `callable | Operator`
    >the function that turns indices into values
- `basis`: `RepresentationBasis`
    >the basis quanta used in the representations
- `logger`: `None | Logger`
    >logger for printing out debug info

<a id="Psience.BasisReps.Terms.Representation.parallelizer" class="docs-object-method">&nbsp;</a>
```python
@property
parallelizer(self): 
```

<a id="Psience.BasisReps.Terms.Representation.compute" class="docs-object-method">&nbsp;</a>
```python
compute(self, inds, **kwargs): 
```

<a id="Psience.BasisReps.Terms.Representation.compute_cached" class="docs-object-method">&nbsp;</a>
```python
compute_cached(self, inds): 
```

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
get_brakets(self, states, check_orthogonality=True): 
```
Computes term elements based on getting a BraKetSpace.
        Can directly pass element specs through, since the shape management shit
        is managed by the BraKetSpace
- `states`: `BraKetSpace | Tuple[np.ndarray, np.ndarray]`
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

<a id="Psience.BasisReps.Terms.Representation.selection_rules" class="docs-object-method">&nbsp;</a>
```python
@property
selection_rules(self): 
```

- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Terms.Representation.get_transformed_space" class="docs-object-method">&nbsp;</a>
```python
get_transformed_space(self, space, parallelizer=None, logger=None, **opts): 
```
Returns the state space obtained by using the
        held operator to transform `space`
- `space`: `Any`
    >No description...
- `:returns`: `SelectionRuleStateSpace`
    >connected state spaces

<a id="Psience.BasisReps.Terms.Representation.apply" class="docs-object-method">&nbsp;</a>
```python
apply(self, other): 
```

<a id="Psience.BasisReps.Terms.Representation.get_representation_matrix" class="docs-object-method">&nbsp;</a>
```python
get_representation_matrix(self, coupled_space, total_space, filter_space=None, diagonal=False, logger=None, zero_element_warning=True, clear_sparse_caches=True, clear_operator_caches=True): 
```
Actively constructs a perturbation theory Hamiltonian representation
- `h`: `Any`
    >No description...
- `cs`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Terms.Representation.get_diagonal_representation" class="docs-object-method">&nbsp;</a>
```python
get_diagonal_representation(self, coupled_space, total_space, logger=None, zero_element_warning=True, clear_sparse_caches=True): 
```
Actively constructs a perturbation theory Hamiltonian representation
- `h`: `Any`
    >No description...
- `cs`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

### Examples


