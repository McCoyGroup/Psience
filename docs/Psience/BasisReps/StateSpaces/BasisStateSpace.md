## <a id="Psience.BasisReps.StateSpaces.BasisStateSpace">BasisStateSpace</a>
Represents a subspace of states inside a representation basis.
Useful largely to provide consistent, unambiguous representations of multiple states across
the different representation-generating methods in the code base.

### Properties and Methods
```python
from_state: method
from_quanta: method
```
<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, basis, states, mode=None): 
```

- `basis`: `RepresentationBasis`
    >No description...
- `states`: `Iterable[int]`
    >No description...
- `mode`: `None | str | StateSpaceSpec`
    >whether the states were supplied as indices or as excitations

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.check_indices" class="docs-object-method">&nbsp;</a>
```python
check_indices(self): 
```

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.to_state" class="docs-object-method">&nbsp;</a>
```python
to_state(self, serializer=None): 
```

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.indices" class="docs-object-method">&nbsp;</a>
```python
@property
indices(self): 
```
Returns held indices
- `inds`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.excitations" class="docs-object-method">&nbsp;</a>
```python
@property
excitations(self): 
```
Returns held excitations
- `inds`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.get_mode" class="docs-object-method">&nbsp;</a>
```python
get_mode(self): 
```

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.infer_state_inds_type" class="docs-object-method">&nbsp;</a>
```python
infer_state_inds_type(self): 
```

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.as_excitations" class="docs-object-method">&nbsp;</a>
```python
as_excitations(self): 
```
Returns states as sets of excitations, rather than indices indo the basis functions.
        For 1D, this just means converting a list of states into tuples of length 1.
- `states`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.as_indices" class="docs-object-method">&nbsp;</a>
```python
as_indices(self): 
```
Returns states as sets of excitations, rather than indices indo the basis functions.
        For 1D, this just means converting a list of states into tuples of length 1.
- `states`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.to_single" class="docs-object-method">&nbsp;</a>
```python
to_single(self): 
```
Basically a no-op
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.is_unique" class="docs-object-method">&nbsp;</a>
```python
is_unique(self): 
```
Returns `True` if the number of states is equal to number of unique states
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.is_sorted" class="docs-object-method">&nbsp;</a>
```python
is_sorted(self, allow_indeterminate=True): 
```
Checks and then sets a flag
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.take_unique" class="docs-object-method">&nbsp;</a>
```python
take_unique(self, sort=False, use_indices=False): 
```
Returns only the unique states, but preserves
        ordering and all of that unless explicitly allowed not to
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.as_sorted" class="docs-object-method">&nbsp;</a>
```python
as_sorted(self): 
```
Returns only the unique states, but preserves
        ordering and all of that unless explicitly allowed not to
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.apply_selection_rules" class="docs-object-method">&nbsp;</a>
```python
apply_selection_rules(self, selection_rules, target_dimensions=None, filter_space=None, parallelizer=None, logger=None, iterations=1, new_state_space_class=None): 
```
Generates a new state space from the application of `selection_rules` to the state space.
        Returns a `BasisMultiStateSpace` where each state tracks the effect of the application of the selection rules
        up to the number of iteration specified.
- `basis`: `Any`
    >No description...
- `selection_rules`: `Any`
    >No description...
- `states`: `Any`
    >No description...
- `iterations`: `Any`
    >No description...
- `filter_space`: `Any`
    >No description...
- `:returns`: `SelectionRuleStateSpace`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.permutationally_reduce" class="docs-object-method">&nbsp;</a>
```python
permutationally_reduce(self): 
```

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.get_representation_indices" class="docs-object-method">&nbsp;</a>
```python
get_representation_indices(self, other=None, selection_rules=None, freqs=None, freq_threshold=None, filter=None, return_filter=False, parallelizer=None): 
```
Generates a set of indices that can be fed into a `Representation` to provide a sub-representation
        in this state space.
        Basically just takes all pairs of indices.
        Only returns the upper-triangle indices
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.get_representation_brakets" class="docs-object-method">&nbsp;</a>
```python
get_representation_brakets(self, other=None, selection_rules=None, freqs=None, freq_threshold=None, filter=None, return_filter=False): 
```
Generates a `BraKetSpace` that can be fed into a `Representation`
        Only returns the upper-triangle pairs because we assume symmetry
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.take_subspace" class="docs-object-method">&nbsp;</a>
```python
take_subspace(self, sel, assume_sorted=False): 
```
Returns a subsample of the space.
        Intended to be a cheap operation, so samples
        along either the indices or the excitations, depending
        on which we have
        If we know the subsample is sorted then we can actually reuse more information
        and so we make use of that
- `sel`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.take_subdimensions" class="docs-object-method">&nbsp;</a>
```python
take_subdimensions(self, inds): 
```
Returns a subsample of the space with some dimensions
        dropped
- `inds`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.take_states" class="docs-object-method">&nbsp;</a>
```python
take_states(self, states, sort=False, assume_sorted=False): 
```
Takes the set of specified states from the space.
        A lot like take_subspace, but operates on states, not indices
- `states`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.drop_subspace" class="docs-object-method">&nbsp;</a>
```python
drop_subspace(self, sel): 
```
Returns a subsample of the space.
        Intended to be a cheap operation, so samples
        along either the indices or the excitations, depending
        on which we have
- `sel`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.drop_subdimensions" class="docs-object-method">&nbsp;</a>
```python
drop_subdimensions(self, inds): 
```
Returns a subsample of the space with some dimensions
        dropped
- `inds`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.drop_states" class="docs-object-method">&nbsp;</a>
```python
drop_states(self, states): 
```
Takes the set of specified states from the space.
        A lot like take_subspace, but operates on states, not indices
- `states`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.split" class="docs-object-method">&nbsp;</a>
```python
split(self, chunksize): 
```
Splits the space up into chunks of at max chunksize
- `chunksize`: `int`
    >No description...
- `:returns`: `Iterable[BasisStateSpace]`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.concatenate" class="docs-object-method">&nbsp;</a>
```python
concatenate(self, other): 
```
Just does a direct concatenation with no unions or any
        of that
- `other`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.union" class="docs-object-method">&nbsp;</a>
```python
union(self, other, sort=False, use_indices=False): 
```
Returns a merged version of self and other, making
        use of as much of the information inherent in both as is possible
- `other`: `BasisStateSpace`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.intersection" class="docs-object-method">&nbsp;</a>
```python
intersection(self, other, sort=False, use_indices=False): 
```
Returns an intersected self and other
- `other`: `BasisStateSpace`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.difference" class="docs-object-method">&nbsp;</a>
```python
difference(self, other, sort=False, use_indices=False): 
```
Returns an diff'ed self and other
- `other`: `BasisStateSpace`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.__repr__" class="docs-object-method">&nbsp;</a>
```python
__repr__(self): 
```

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.__eq__" class="docs-object-method">&nbsp;</a>
```python
__eq__(self, other): 
```

- `other`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

### Examples


