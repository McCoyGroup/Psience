## <a id="Psience.BasisReps.StateSpaces.BasisStateSpace">BasisStateSpace</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces.py#L494)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces.py#L494?message=Update%20Docs)]
</div>

Represents a subspace of states inside a representation basis.
Useful largely to provide consistent, unambiguous representations of multiple states across
the different representation-generating methods in the code base.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, basis, states, full_basis=None, mode=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BasisStateSpace.py#L500)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BasisStateSpace.py#L500?message=Update%20Docs)]
</div>

  - `basis`: `RepresentationBasis`
    > 
  - `states`: `Iterable[int]`
    > 
  - `mode`: `None | str | StateSpaceSpec`
    > whether the states were supplied as indices or as excitations


<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.check_indices" class="docs-object-method">&nbsp;</a> 
```python
check_indices(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BasisStateSpace.py#L553)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BasisStateSpace.py#L553?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.to_state" class="docs-object-method">&nbsp;</a> 
```python
to_state(self, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BasisStateSpace.py#L563)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BasisStateSpace.py#L563?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.from_state" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_state(cls, data, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L569)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L569?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.indices" class="docs-object-method">&nbsp;</a> 
```python
@property
indices(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BasisStateSpace.py#L581)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BasisStateSpace.py#L581?message=Update%20Docs)]
</div>
Returns held indices
  - `inds`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.excitations" class="docs-object-method">&nbsp;</a> 
```python
@property
excitations(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BasisStateSpace.py#L603)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BasisStateSpace.py#L603?message=Update%20Docs)]
</div>
Returns held excitations
  - `inds`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.from_quanta" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_quanta(cls, basis, quants): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L624)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L624?message=Update%20Docs)]
</div>
Returns states with `quants` quanta of excitation
using the basis `basis`
  - `basis`: `RepresentationBasis`
    > 
  - `quants`: `int | Iterable[int]`
    > set of integers
  - `:returns`: `_`
    > B
a
s
i
s
S
t
a
t
e
S
p
a
c
e


<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.states_in_windows" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
states_in_windows(cls, freqs, windows: 'list[[int,int]]', max_state=None, min_quantas=None, max_quantas=None, initial_state=None, fixed_modes=None, basis=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L650)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L650?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.states_under_freq_threshold" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
states_under_freq_threshold(cls, freqs, thresh, min_freq=None, max_state=None, min_quanta=None, max_quanta=None, basis=None, fixed_modes=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L730)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L730?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.get_mode" class="docs-object-method">&nbsp;</a> 
```python
get_mode(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BasisStateSpace.py#L747)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BasisStateSpace.py#L747?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.infer_state_inds_type" class="docs-object-method">&nbsp;</a> 
```python
infer_state_inds_type(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BasisStateSpace.py#L750)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BasisStateSpace.py#L750?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.as_excitations" class="docs-object-method">&nbsp;</a> 
```python
as_excitations(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BasisStateSpace.py#L761)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BasisStateSpace.py#L761?message=Update%20Docs)]
</div>
Returns states as sets of excitations, rather than indices indo the basis functions.
For 1D, this just means converting a list of states into tuples of length 1.
  - `states`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.as_indices" class="docs-object-method">&nbsp;</a> 
```python
as_indices(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BasisStateSpace.py#L797)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BasisStateSpace.py#L797?message=Update%20Docs)]
</div>
Returns states as sets of excitations, rather than indices indo the basis functions.
For 1D, this just means converting a list of states into tuples of length 1.
  - `states`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.to_single" class="docs-object-method">&nbsp;</a> 
```python
to_single(self, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BasisStateSpace.py#L827)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BasisStateSpace.py#L827?message=Update%20Docs)]
</div>
Basically a no-op
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.is_unique" class="docs-object-method">&nbsp;</a> 
```python
is_unique(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BasisStateSpace.py#L843)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BasisStateSpace.py#L843?message=Update%20Docs)]
</div>
Returns `True` if the number of states is equal to number of unique states
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.is_sorted" class="docs-object-method">&nbsp;</a> 
```python
is_sorted(self, allow_indeterminate=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BasisStateSpace.py#L854)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BasisStateSpace.py#L854?message=Update%20Docs)]
</div>
Checks and then sets a flag
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.take_unique" class="docs-object-method">&nbsp;</a> 
```python
take_unique(self, sort=False, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BasisStateSpace.py#L867)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BasisStateSpace.py#L867?message=Update%20Docs)]
</div>
Returns only the unique states, but preserves
ordering and all of that unless explicitly allowed not to
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.as_sorted" class="docs-object-method">&nbsp;</a> 
```python
as_sorted(self, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BasisStateSpace.py#L911)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BasisStateSpace.py#L911?message=Update%20Docs)]
</div>
Returns a sorted version of the state space
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.apply_selection_rules" class="docs-object-method">&nbsp;</a> 
```python
apply_selection_rules(self, selection_rules, target_dimensions=None, filter_space=None, parallelizer=None, logger=None, iterations=1, new_state_space_class=None, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BasisStateSpace.py#L949)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BasisStateSpace.py#L949?message=Update%20Docs)]
</div>
Generates a new state space from the application of `selection_rules` to the state space.
Returns a `BasisMultiStateSpace` where each state tracks the effect of the application of the selection rules
up to the number of iteration specified.
  - `basis`: `Any`
    > 
  - `selection_rules`: `Any`
    > 
  - `states`: `Any`
    > 
  - `iterations`: `Any`
    > 
  - `filter_space`: `Any`
    > 
  - `:returns`: `SelectionRuleStateSpace`
    >


<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.permutationally_reduce" class="docs-object-method">&nbsp;</a> 
```python
permutationally_reduce(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BasisStateSpace.py#L982)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BasisStateSpace.py#L982?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.get_representation_indices" class="docs-object-method">&nbsp;</a> 
```python
get_representation_indices(self, other=None, selection_rules=None, freqs=None, freq_threshold=None, filter=None, return_filter=False, parallelizer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BasisStateSpace.py#L985)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BasisStateSpace.py#L985?message=Update%20Docs)]
</div>
Generates a set of indices that can be fed into a `Representation` to provide a sub-representation
in this state space.
Basically just takes all pairs of indices.
Only returns the upper-triangle indices
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.get_representation_brakets" class="docs-object-method">&nbsp;</a> 
```python
get_representation_brakets(self, other=None, selection_rules=None, freqs=None, freq_threshold=None, filter=None, return_filter=False, track_excitations=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BasisStateSpace.py#L1068)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BasisStateSpace.py#L1068?message=Update%20Docs)]
</div>
Generates a `BraKetSpace` that can be fed into a `Representation`
Only returns the upper-triangle pairs because we assume symmetry
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.take_subspace" class="docs-object-method">&nbsp;</a> 
```python
take_subspace(self, sel, assume_sorted=False, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BasisStateSpace.py#L1141)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BasisStateSpace.py#L1141?message=Update%20Docs)]
</div>
Returns a subsample of the space.
Intended to be a cheap operation, so samples
along either the indices or the excitations, depending
on which we have
If we know the subsample is sorted then we can actually reuse more information
and so we make use of that
  - `sel`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.take_subdimensions" class="docs-object-method">&nbsp;</a> 
```python
take_subdimensions(self, inds, exc=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BasisStateSpace.py#L1198)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BasisStateSpace.py#L1198?message=Update%20Docs)]
</div>
Returns a subsample of the space with some dimensions
dropped
  - `inds`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.take_states" class="docs-object-method">&nbsp;</a> 
```python
take_states(self, states, sort=False, assume_sorted=False, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BasisStateSpace.py#L1216)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BasisStateSpace.py#L1216?message=Update%20Docs)]
</div>
Takes the set of specified states from the space.
A lot like take_subspace, but operates on states, not indices
  - `states`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.drop_subspace" class="docs-object-method">&nbsp;</a> 
```python
drop_subspace(self, sel, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BasisStateSpace.py#L1238)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BasisStateSpace.py#L1238?message=Update%20Docs)]
</div>
Returns a subsample of the space.
Intended to be a cheap operation, so samples
along either the indices or the excitations, depending
on which we have
  - `sel`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.drop_subdimensions" class="docs-object-method">&nbsp;</a> 
```python
drop_subdimensions(self, inds, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BasisStateSpace.py#L1262)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BasisStateSpace.py#L1262?message=Update%20Docs)]
</div>
Returns a subsample of the space with some dimensions
dropped
  - `inds`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.drop_states" class="docs-object-method">&nbsp;</a> 
```python
drop_states(self, states, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BasisStateSpace.py#L1280)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BasisStateSpace.py#L1280?message=Update%20Docs)]
</div>
Takes the set of specified states from the space.
A lot like take_subspace, but operates on states, not indices
  - `states`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.split" class="docs-object-method">&nbsp;</a> 
```python
split(self, chunksize): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BasisStateSpace.py#L1299)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BasisStateSpace.py#L1299?message=Update%20Docs)]
</div>
Splits the space up into chunks of at max chunksize
  - `chunksize`: `int`
    > 
  - `:returns`: `Iterable[BasisStateSpace]`
    >


<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.concatenate" class="docs-object-method">&nbsp;</a> 
```python
concatenate(self, other, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BasisStateSpace.py#L1333)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BasisStateSpace.py#L1333?message=Update%20Docs)]
</div>
Just does a direct concatenation with no unions or any
of that
  - `other`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.union" class="docs-object-method">&nbsp;</a> 
```python
union(self, other, sort=False, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BasisStateSpace.py#L1397)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BasisStateSpace.py#L1397?message=Update%20Docs)]
</div>
Returns a merged version of self and other, making
use of as much of the information inherent in both as is possible
  - `other`: `BasisStateSpace`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.intersection" class="docs-object-method">&nbsp;</a> 
```python
intersection(self, other, sort=False, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BasisStateSpace.py#L1520)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BasisStateSpace.py#L1520?message=Update%20Docs)]
</div>
Returns an intersected self and other
  - `other`: `BasisStateSpace`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.difference" class="docs-object-method">&nbsp;</a> 
```python
difference(self, other, sort=False, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BasisStateSpace.py#L1648)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BasisStateSpace.py#L1648?message=Update%20Docs)]
</div>
Returns an diff'ed self and other
  - `other`: `BasisStateSpace`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BasisStateSpace.py#L1756)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BasisStateSpace.py#L1756?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.__eq__" class="docs-object-method">&nbsp;</a> 
```python
__eq__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BasisStateSpace.py#L1774)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BasisStateSpace.py#L1774?message=Update%20Docs)]
</div>

  - `other`: `Any`
    > 
  - `:returns`: `_`
    >
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/StateSpaces/BasisStateSpace.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/StateSpaces/BasisStateSpace.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/StateSpaces/BasisStateSpace.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/StateSpaces/BasisStateSpace.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces.py#L494?message=Update%20Docs)   
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