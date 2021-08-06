## <a id="Psience.BasisReps.StateSpaces.AbstractStateSpace">AbstractStateSpace</a>
Represents a generalized state space which will provide core
methods to index into a basis and generate representations

### Properties and Methods
```python
keep_excitations: bool
keep_indices: bool
StateSpaceSpec: EnumMeta
StateSpaceCache: type
excitations_dtype: dtype[int8]
indices_dtype: dtype[uint64]
from_state: method
get_states_with_quanta: method
```
<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, basis): 
```

- `basis`: `RepresentationBasis`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.to_state" class="docs-object-method">&nbsp;</a>
```python
to_state(self, serializer=None): 
```
Provides just the state that is needed to
        serialize the object
- `serializer`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.ndim" class="docs-object-method">&nbsp;</a>
```python
@property
ndim(self): 
```

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.excitations" class="docs-object-method">&nbsp;</a>
```python
@property
excitations(self): 
```

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.mode" class="docs-object-method">&nbsp;</a>
```python
@property
mode(self): 
```

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.get_mode" class="docs-object-method">&nbsp;</a>
```python
get_mode(self): 
```
Returns the mode (indices or excitations) for the held states
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.has_indices" class="docs-object-method">&nbsp;</a>
```python
@property
has_indices(self): 
```

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.has_excitations" class="docs-object-method">&nbsp;</a>
```python
@property
has_excitations(self): 
```

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.indices" class="docs-object-method">&nbsp;</a>
```python
@property
indices(self): 
```

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.indexer" class="docs-object-method">&nbsp;</a>
```python
@property
indexer(self): 
```

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.exc_indexer" class="docs-object-method">&nbsp;</a>
```python
@property
exc_indexer(self): 
```

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.find" class="docs-object-method">&nbsp;</a>
```python
find(self, to_search, check=True): 
```
Finds the indices of a set of indices inside the space
- `to_search`: `np.ndarray | AbstractStateSpace`
    >array of ints
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.__len__" class="docs-object-method">&nbsp;</a>
```python
__len__(self): 
```

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.unique_len" class="docs-object-method">&nbsp;</a>
```python
@property
unique_len(self): 
```

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.unique_indices" class="docs-object-method">&nbsp;</a>
```python
@property
unique_indices(self): 
```
Returns the unique indices
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.unique_excitations" class="docs-object-method">&nbsp;</a>
```python
@property
unique_excitations(self): 
```
Returns the unique excitations
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.as_indices" class="docs-object-method">&nbsp;</a>
```python
as_indices(self): 
```
Returns the index version of the stored states
- `:returns`: `np.ndarray`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.as_unique_indices" class="docs-object-method">&nbsp;</a>
```python
as_unique_indices(self, sort=False): 
```
Returns unique indices
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.as_excitations" class="docs-object-method">&nbsp;</a>
```python
as_excitations(self): 
```
Returns the excitation version of the stored states
- `:returns`: `np.ndarray`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.as_unique_excitations" class="docs-object-method">&nbsp;</a>
```python
as_unique_excitations(self, sort=False): 
```
Returns unique excitations
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.get_representation_indices" class="docs-object-method">&nbsp;</a>
```python
get_representation_indices(self, other=None, selection_rules=None, freqs=None, freq_threshold=None): 
```
Returns bra and ket indices that can be used as indices to generate representations
- `other`: `Any`
    >No description...
- `selection_rules`: `Any`
    >No description...
- `freqs`: `Any`
    >No description...
- `freq_threshold`: `Any`
    >No description...
- `:returns`: `(np.ndarray, np.ndarray)`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.get_representation_brakets" class="docs-object-method">&nbsp;</a>
```python
get_representation_brakets(self, other=None, selection_rules=None, freqs=None, freq_threshold=None): 
```
Returns a BraKetSpace that can be used as generate representations
- `other`: `Any`
    >No description...
- `selection_rules`: `Any`
    >No description...
- `freqs`: `Any`
    >No description...
- `freq_threshold`: `Any`
    >No description...
- `:returns`: `BraKetSpace`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.take_states" class="docs-object-method">&nbsp;</a>
```python
take_states(self, states): 
```
Takes the intersection of self and the specified states
- `states`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.take_subspace" class="docs-object-method">&nbsp;</a>
```python
take_subspace(self, sel): 
```
Takes a subset of the states
- `sel`: `Any`
    >No description...
- `:returns`: `AbstractStateSpace`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.take_subdimensions" class="docs-object-method">&nbsp;</a>
```python
take_subdimensions(self, inds): 
```
Takes a subset of the state dimensions
- `sel`: `Any`
    >No description...
- `:returns`: `AbstractStateSpace`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.drop_states" class="docs-object-method">&nbsp;</a>
```python
drop_states(self, states): 
```
Takes the difference of self and the specified states
- `states`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.drop_subspace" class="docs-object-method">&nbsp;</a>
```python
drop_subspace(self, sel): 
```
Drops a subset of the states
- `sel`: `Any`
    >No description...
- `:returns`: `AbstractStateSpace`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.drop_subdimensions" class="docs-object-method">&nbsp;</a>
```python
drop_subdimensions(self, inds): 
```
Drops a subset of the state dimensions
- `sel`: `Any`
    >No description...
- `:returns`: `AbstractStateSpace`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.to_single" class="docs-object-method">&nbsp;</a>
```python
to_single(self, track_excitations=True, track_indices=True): 
```
Flattens any complicated state space structure into a
        single space like a `BasisStateSpace`
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.split" class="docs-object-method">&nbsp;</a>
```python
split(self, chunksize): 
```
Subclass overridable function to allow for spaces to be
        split up into chunks
- `chunksize`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

### Examples




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/edit/ci/examples/ci/docs/Psience/BasisReps/StateSpaces/AbstractStateSpace.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/examples/ci/docs/Psience/BasisReps/StateSpaces/AbstractStateSpace.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/edit/ci/docs/ci/docs/Psience/BasisReps/StateSpaces/AbstractStateSpace.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/docs/templates/ci/docs/Psience/BasisReps/StateSpaces/AbstractStateSpace.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py?message=Update%20Docs)