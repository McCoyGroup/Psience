## <a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace">SelectionRuleStateSpace</a>
A `BasisMultiStateSpace` subclass that is only built from applying selection rules to an initial space

### Properties and Methods
```python
from_state: method
from_rules: method
```
<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, init_space, excitations, selection_rules=None): 
```

- `init_space`: `Any`
    >No description...
- `excitations`: `Any`
    >No description...
- `selection_rules`: `Any`
    >No description...

<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.to_state" class="docs-object-method">&nbsp;</a>
```python
to_state(self, serializer=None): 
```

<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.representative_space" class="docs-object-method">&nbsp;</a>
```python
@property
representative_space(self): 
```

<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.take_states" class="docs-object-method">&nbsp;</a>
```python
take_states(self, states): 
```
Takes the intersection of each held space and the specified states
- `states`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.take_subdimensions" class="docs-object-method">&nbsp;</a>
```python
take_subdimensions(self, inds): 
```
Takes the subdimensions from each space
- `inds`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.get_representation_indices" class="docs-object-method">&nbsp;</a>
```python
get_representation_indices(self, other=None, freqs=None, freq_threshold=None, selection_rules=None): 
```
This is where this pays dividends, as we know that only the init_space and the held excitations can couple
        which reduces the combinatoric work by a factor of like 2.
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.filter_representation_inds" class="docs-object-method">&nbsp;</a>
```python
filter_representation_inds(self, ind_pairs, q_changes): 
```
Filters representation indices by the allowed #quantum changes.
        Not sure I'll even need this, if `get_representation_indices` is tight enough.
- `ind_pairs`: `Any`
    >No description...
- `q_changes`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.__getitem__" class="docs-object-method">&nbsp;</a>
```python
__getitem__(self, item): 
```

### Examples


