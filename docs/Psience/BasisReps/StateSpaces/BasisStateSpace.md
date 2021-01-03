## <a id="Psience.BasisReps.StateSpaces.BasisStateSpace">BasisStateSpace</a>
Represents a subspace of states inside a representation basis.
Useful largely to provide consistent, unambiguous representations of multiple states across
the different representation-generating methods in the code base.

### Properties and Methods
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

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.apply_selection_rules" class="docs-object-method">&nbsp;</a>
```python
apply_selection_rules(self, selection_rules, filter_space=None, iterations=1): 
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

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.get_representation_indices" class="docs-object-method">&nbsp;</a>
```python
get_representation_indices(self, other=None, selection_rules=None, freqs=None, freq_threshold=None): 
```
Generates a set of indices that can be fed into a `Representation` to provide a sub-representation
        in this state space.
        Basically just takes all pairs of indices.
        Only returns the upper-triangle indices
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.get_representation_brakets" class="docs-object-method">&nbsp;</a>
```python
get_representation_brakets(self, other=None, selection_rules=None, freqs=None, freq_threshold=None): 
```
Generates a `BraKetSpace` that can be fed into a `Representation`
        Basically just takes all pairs of indices.
        Only returns the upper-triangle indices
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.take_intersection" class="docs-object-method">&nbsp;</a>
```python
take_intersection(self, states): 
```
Takes the intersection of self and the specified states
- `states`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.take_subspace" class="docs-object-method">&nbsp;</a>
```python
take_subspace(self, sel): 
```

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.take_subdimensions" class="docs-object-method">&nbsp;</a>
```python
take_subdimensions(self, inds): 
```

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.__repr__" class="docs-object-method">&nbsp;</a>
```python
__repr__(self): 
```

### Examples


