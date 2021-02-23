## <a id="Psience.BasisReps.StateSpaces.BraKetSpace">BraKetSpace</a>
Represents a set of pairs of states that can be fed into a `Representation` or `Operator`
to efficiently tell it what terms it need to calculate.
This basically just implements a bunch of stuff for generating a Graph defining
the connections between states.

### Properties and Methods
```python
from_indices: method
```
<a id="Psience.BasisReps.StateSpaces.BraKetSpace.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, bra_space, ket_space): 
```

- `bra_space`: `AbstractStateSpace`
    >No description...
- `ket_space`: `AbstractStateSpace`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.__len__" class="docs-object-method">&nbsp;</a>
```python
__len__(self): 
```

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.__repr__" class="docs-object-method">&nbsp;</a>
```python
__repr__(self): 
```

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.load_non_orthog" class="docs-object-method">&nbsp;</a>
```python
load_non_orthog(self): 
```

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.get_non_orthog" class="docs-object-method">&nbsp;</a>
```python
get_non_orthog(self, inds, assume_unique=False): 
```
Returns whether the states are non-orthogonal under the set of indices.
- `inds`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.get_sel_rules_from1d" class="docs-object-method">&nbsp;</a>
```python
get_sel_rules_from1d(self, inds, rules): 
```

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.get_sel_rule_filter" class="docs-object-method">&nbsp;</a>
```python
get_sel_rule_filter(self, rules): 
```

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.take_subspace" class="docs-object-method">&nbsp;</a>
```python
take_subspace(self, sel): 
```

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.take_subdimensions" class="docs-object-method">&nbsp;</a>
```python
take_subdimensions(self, inds): 
```

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.apply_non_orthogonality" class="docs-object-method">&nbsp;</a>
```python
apply_non_orthogonality(self, inds, assume_unique=False): 
```
Takes the bra-ket pairs that are non-orthogonal under the
        indices `inds`
- `inds`: `Any`
    >No description...
- `assume_unique`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.apply_sel_rules" class="docs-object-method">&nbsp;</a>
```python
apply_sel_rules(self, rules): 
```
Applies selections rules
- `rules`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.adjacency_matrix" class="docs-object-method">&nbsp;</a>
```python
adjacency_matrix(self, total_space=None): 
```
Generates the (sparse) unweighted adjacency matrix for the bras & kets
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.split" class="docs-object-method">&nbsp;</a>
```python
split(self, chunksize): 
```
splits the brakets into blocks of at max chunksize
- `chunksize`: `int`
    >No description...
- `:returns`: `Iterable[BraKetSpace]`
    >No description...

### Examples


