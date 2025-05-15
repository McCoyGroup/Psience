## <a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace">SelectionRuleStateSpace</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces.py#L2328)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces.py#L2328?message=Update%20Docs)]
</div>

A `BasisMultiStateSpace` subclass that is only built from applying selection rules to an initial space
This really should have been called `TransformedStateSpace` but I am dumb







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
track_change_positions: bool
direct_sum_chunk_size: int
```
<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, init_space, excitations, selection_rules=None, ignore_shapes=False, changes=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L2333)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L2333?message=Update%20Docs)]
</div>

  - `init_space`: `Any`
    > 
  - `excitations`: `Any`
    > 
  - `selection_rules`: `Any`
    >


<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.to_state" class="docs-object-method">&nbsp;</a> 
```python
to_state(self, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L2363)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L2363?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.from_state" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_state(cls, data, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L2369)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L2369?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.as_indices" class="docs-object-method">&nbsp;</a> 
```python
as_indices(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L2388)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L2388?message=Update%20Docs)]
</div>
Pulls the full set indices out of all of the
held spaces and returns them as a flat vector
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.as_excitations" class="docs-object-method">&nbsp;</a> 
```python
as_excitations(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L2400)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L2400?message=Update%20Docs)]
</div>
Pulls the full set excitations out of all of the
held spaces and returns them as a flat vector
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.representative_space" class="docs-object-method">&nbsp;</a> 
```python
@property
representative_space(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L2414)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L2414?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.nstates" class="docs-object-method">&nbsp;</a> 
```python
@property
nstates(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L2418)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L2418?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.check_indices" class="docs-object-method">&nbsp;</a> 
```python
check_indices(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L2422)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L2422?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.take_states" class="docs-object-method">&nbsp;</a> 
```python
take_states(self, states, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L2426)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L2426?message=Update%20Docs)]
</div>
Takes the intersection of each held space and the specified states
  - `states`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.take_subspace" class="docs-object-method">&nbsp;</a> 
```python
take_subspace(self, states): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L2457)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L2457?message=Update%20Docs)]
</div>
Takes the intersection of each held space and the specified states
  - `states`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.take_subdimensions" class="docs-object-method">&nbsp;</a> 
```python
take_subdimensions(self, inds): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L2488)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L2488?message=Update%20Docs)]
</div>
Takes the subdimensions from each space
  - `inds`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.drop_states" class="docs-object-method">&nbsp;</a> 
```python
drop_states(self, states): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L2515)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L2515?message=Update%20Docs)]
</div>
Takes the intersection of each held space and the specified states
  - `states`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.drop_subspace" class="docs-object-method">&nbsp;</a> 
```python
drop_subspace(self, inds): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L2539)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L2539?message=Update%20Docs)]
</div>
Takes the intersection of each held space and the specified states
  - `states`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.drop_subdimensions" class="docs-object-method">&nbsp;</a> 
```python
drop_subdimensions(self, inds): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L2564)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L2564?message=Update%20Docs)]
</div>
Takes the subdimensions from each space
  - `inds`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.get_representation_indices" class="docs-object-method">&nbsp;</a> 
```python
get_representation_indices(self, other=None, freqs=None, freq_threshold=None, selection_rules=None, filter=None, return_filter=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L2581)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L2581?message=Update%20Docs)]
</div>
This is where this pays dividends, as we know that only the init_space and the held excitations can couple
which reduces the combinatoric work by a factor of like 2.
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.filter_representation_inds" class="docs-object-method">&nbsp;</a> 
```python
filter_representation_inds(self, ind_pairs, q_changes): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L2658)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L2658?message=Update%20Docs)]
</div>
Filters representation indices by the allowed #quantum changes.
Not sure I'll even need this, if `get_representation_indices` is tight enough.
  - `ind_pairs`: `Any`
    > 
  - `q_changes`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.get_representation_brakets" class="docs-object-method">&nbsp;</a> 
```python
get_representation_brakets(self, freqs=None, freq_threshold=None, other=None, selection_rules=None, filter=None, return_filter=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L2687)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L2687?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.from_rules" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_rules(cls, space, selection_rules, target_dimensions=None, filter_space=None, iterations=1, method='new', parallelizer=None, chunk_size=None, logger=None, track_excitations=True, track_indices=True, full_basis=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L2907)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L2907?message=Update%20Docs)]
</div>

  - `space`: `BasisStateSpace | BasisMultiStateSpace`
    > initial space to which to apply the transformations
  - `selection_rules`: `Iterable[Iterable[int]]`
    > different possible transformations
  - `iterations`: `int`
    > number of times to apply the transformations
  - `filter_space`: `BasisStateSpace | None`
    > a space within which all generated `BasisStateSpace` objects must be contained
  - `:returns`: `SelectionRuleStateSpace`
    >


<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.filter_transitions" class="docs-object-method">&nbsp;</a> 
```python
filter_transitions(self, excluded_transitions, in_place=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L3135)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L3135?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.generate_change_indices" class="docs-object-method">&nbsp;</a> 
```python
generate_change_indices(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L3229)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L3229?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.union" class="docs-object-method">&nbsp;</a> 
```python
union(self, other, handle_subspaces=True, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L3241)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L3241?message=Update%20Docs)]
</div>
Returns a merged version of self and other, adding
any states in other to self and merging where they intersect
  - `other`: `SelectionRuleStateSpace`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.intersection" class="docs-object-method">&nbsp;</a> 
```python
intersection(self, other, handle_subspaces=True, use_indices=False, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L3406)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L3406?message=Update%20Docs)]
</div>
Returns an intersected self and other
  - `other`: `SelectionRuleStateSpace`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.difference" class="docs-object-method">&nbsp;</a> 
```python
difference(self, other, handle_subspaces=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L3515)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L3515?message=Update%20Docs)]
</div>
Returns an diff'ed self and other.
We get fundamentally different behaviour for `handle_subspaces` than without it.
If we have it _on_ then differences are computed for each states in the intersection of
the primary (key) states.
If we have it off then the difference in the key states is computed and nothing more is
done.
  - `other`: `SelectionRuleStateSpace`
    > 
  - `:returns`: `SelectionRuleStateSpace`
    >


<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L3584)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L3584?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.__setitem__" class="docs-object-method">&nbsp;</a> 
```python
__setitem__(self, item, vals): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L3592)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L3592?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.map" class="docs-object-method">&nbsp;</a> 
```python
map(self, f): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L3599)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L3599?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.SelectionRuleStateSpace.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L3610)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/SelectionRuleStateSpace.py#L3610?message=Update%20Docs)]
</div>
 </div>
</div>












---


<div markdown="1" class="text-secondary">
<div class="container">
  <div class="row">
   <div class="col" markdown="1">
**Feedback**   
</div>
   <div class="col" markdown="1">
**Examples**   
</div>
   <div class="col" markdown="1">
**Templates**   
</div>
   <div class="col" markdown="1">
**Documentation**   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
</div>
  <div class="row">
   <div class="col" markdown="1">
[Bug](https://github.com/McCoyGroup/Psience/issues/new?title=Documentation%20Improvement%20Needed)/[Request](https://github.com/McCoyGroup/Psience/issues/new?title=Example%20Request)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/StateSpaces/SelectionRuleStateSpace.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/StateSpaces/SelectionRuleStateSpace.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/StateSpaces/SelectionRuleStateSpace.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/StateSpaces/SelectionRuleStateSpace.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces.py#L2328?message=Update%20Docs)   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
</div>
</div>
</div>