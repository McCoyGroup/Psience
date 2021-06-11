## <a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix">StateSpaceMatrix</a>
A `SparseArray` that holds onto a `BasisStateSpace` that keeps track of the
total set of states involved.
By default is assumed real-symmetric. This can be relaxed in the future.

TODO: The idea is good, but calculating what is "in" the array and what is "out"
        every single time this is applied could be slow...
      We'll need to test to see how slow

### Properties and Methods
```python
identity_from_space: method
```
<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, initial_basis, initial_vals=None, column_space=None, symmetric=True): 
```

- `initial_basis`: `BasisStateSpace | RepresentationBasis`
    >No description...

<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.shape" class="docs-object-method">&nbsp;</a>
```python
@property
shape(self): 
```

<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.basis" class="docs-object-method">&nbsp;</a>
```python
@property
basis(self): 
```
Returns the basis for the matrix rep
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.brakets" class="docs-object-method">&nbsp;</a>
```python
@property
brakets(self): 
```
Returns the BraKetSpace for the held indices
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.extend_basis" class="docs-object-method">&nbsp;</a>
```python
extend_basis(self, states, extend_columns=True): 
```
Extends the held state space and resizes the held array if need be
- `states`: `BasisStateSpace`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.compute_values" class="docs-object-method">&nbsp;</a>
```python
compute_values(self, func, brakets): 
```
Computes new values into the held `SparseArray` based on the function and brakets provided
        and returns the entire array of values
- `func`: `Any`
    >A function that can take a braket spec and compute values
- `brakets`: `Any`
    >A set of brakets to compute values for
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.dot" class="docs-object-method">&nbsp;</a>
```python
dot(self, other): 
```
Performs a dot product between the held SparseArray and another
        StateSpaceMatrix
- `other`: `StateSpaceMatrix`
    >other matrix
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.__getitem__" class="docs-object-method">&nbsp;</a>
```python
__getitem__(self, item): 
```

<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.__setitem__" class="docs-object-method">&nbsp;</a>
```python
__setitem__(self, item, vals): 
```

<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.__repr__" class="docs-object-method">&nbsp;</a>
```python
__repr__(self): 
```

### Examples


