## <a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace">SelectionRuleStateSpace</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L2182)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L2182?message=Update%20Docs)]
</div>

A `BasisMultiStateSpace` subclass that is only built from applying selection rules to an initial space
This really should have been called `TransformedStateSpace` but I am dumb

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

```python
direct_sum_chunk_size: int
```
<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, init_space, excitations, selection_rules=None, ignore_shapes=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L2187)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L2187?message=Update%20Docs)]
</div>


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
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L2212)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L2212?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.from_state" class="docs-object-method">&nbsp;</a> 
```python
from_state(data, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L2218)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L2218?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.as_indices" class="docs-object-method">&nbsp;</a> 
```python
as_indices(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L2237)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L2237?message=Update%20Docs)]
</div>

Pulls the full set indices out of all of the
        held spaces and returns them as a flat vector
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.as_excitations" class="docs-object-method">&nbsp;</a> 
```python
as_excitations(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L2249)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L2249?message=Update%20Docs)]
</div>

Pulls the full set excitations out of all of the
        held spaces and returns them as a flat vector
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.representative_space" class="docs-object-method">&nbsp;</a> 
```python
@property
representative_space(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.nstates" class="docs-object-method">&nbsp;</a> 
```python
@property
nstates(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.check_indices" class="docs-object-method">&nbsp;</a> 
```python
check_indices(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L2271)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L2271?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.take_states" class="docs-object-method">&nbsp;</a> 
```python
take_states(self, states, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L2275)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L2275?message=Update%20Docs)]
</div>

Takes the intersection of each held space and the specified states
- `states`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.take_subspace" class="docs-object-method">&nbsp;</a> 
```python
take_subspace(self, states): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L2303)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L2303?message=Update%20Docs)]
</div>

Takes the intersection of each held space and the specified states
- `states`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.take_subdimensions" class="docs-object-method">&nbsp;</a> 
```python
take_subdimensions(self, inds): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L2330)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L2330?message=Update%20Docs)]
</div>

Takes the subdimensions from each space
- `inds`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.drop_states" class="docs-object-method">&nbsp;</a> 
```python
drop_states(self, states): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L2354)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L2354?message=Update%20Docs)]
</div>

Takes the intersection of each held space and the specified states
- `states`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.drop_subspace" class="docs-object-method">&nbsp;</a> 
```python
drop_subspace(self, inds): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L2375)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L2375?message=Update%20Docs)]
</div>

Takes the intersection of each held space and the specified states
- `states`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.drop_subdimensions" class="docs-object-method">&nbsp;</a> 
```python
drop_subdimensions(self, inds): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L2396)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L2396?message=Update%20Docs)]
</div>

Takes the subdimensions from each space
- `inds`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.get_representation_indices" class="docs-object-method">&nbsp;</a> 
```python
get_representation_indices(self, other=None, freqs=None, freq_threshold=None, selection_rules=None, filter=None, return_filter=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L2410)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L2410?message=Update%20Docs)]
</div>

This is where this pays dividends, as we know that only the init_space and the held excitations can couple
        which reduces the combinatoric work by a factor of like 2.
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.filter_representation_inds" class="docs-object-method">&nbsp;</a> 
```python
filter_representation_inds(self, ind_pairs, q_changes): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L2488)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L2488?message=Update%20Docs)]
</div>

Filters representation indices by the allowed #quantum changes.
        Not sure I'll even need this, if `get_representation_indices` is tight enough.
- `ind_pairs`: `Any`
    >No description...
- `q_changes`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.from_rules" class="docs-object-method">&nbsp;</a> 
```python
from_rules(space, selection_rules, target_dimensions=None, filter_space=None, iterations=1, method='new', parallelizer=None, chunk_size=None, logger=None, track_excitations=True, track_indices=True, full_basis=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L2682)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L2682?message=Update%20Docs)]
</div>


- `space`: `BasisStateSpace | BasisMultiStateSpace`
    >initial space to which to apply the transformations
- `selection_rules`: `Iterable[Iterable[int]]`
    >different possible transformations
- `iterations`: `int`
    >number of times to apply the transformations
- `filter_space`: `BasisStateSpace | None`
    >a space within which all generated `BasisStateSpace` objects must be contained
- `:returns`: `SelectionRuleStateSpace`
    >No description...

<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.union" class="docs-object-method">&nbsp;</a> 
```python
union(self, other, handle_subspaces=True, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L2955)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L2955?message=Update%20Docs)]
</div>

Returns a merged version of self and other, adding
        any states in other to self and merging where they intersect
- `other`: `SelectionRuleStateSpace`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.intersection" class="docs-object-method">&nbsp;</a> 
```python
intersection(self, other, handle_subspaces=True, use_indices=False, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L3060)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L3060?message=Update%20Docs)]
</div>

Returns an intersected self and other
- `other`: `SelectionRuleStateSpace`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.difference" class="docs-object-method">&nbsp;</a> 
```python
difference(self, other, handle_subspaces=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L3129)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L3129?message=Update%20Docs)]
</div>

Returns an diff'ed self and other.
        We get fundamentally different behaviour for `handle_subspaces` than without it.
        If we have it _on_ then differences are computed for each states in the intersection of
          the primary (key) states.
        If we have it off then the difference in the key states is computed and nothing more is
        done.
- `other`: `SelectionRuleStateSpace`
    >No description...
- `:returns`: `SelectionRuleStateSpace`
    >No description...

<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L3191)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L3191?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.__setitem__" class="docs-object-method">&nbsp;</a> 
```python
__setitem__(self, item, vals): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L3199)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L3199?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L3206)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L3206?message=Update%20Docs)]
</div>

 </div>
</div>




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/StateSpaces/SelectionRuleStateSpace.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/StateSpaces/SelectionRuleStateSpace.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/StateSpaces/SelectionRuleStateSpace.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/StateSpaces/SelectionRuleStateSpace.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L2182?message=Update%20Docs)