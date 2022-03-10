## <a id="Psience.BasisReps.StateSpaces.BasisStateSpace">BasisStateSpace</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L470)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L470?message=Update%20Docs)]
</div>

Represents a subspace of states inside a representation basis.
Useful largely to provide consistent, unambiguous representations of multiple states across
the different representation-generating methods in the code base.

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, basis, states, full_basis=None, mode=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L476)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L476?message=Update%20Docs)]
</div>


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
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L517)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L517?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.to_state" class="docs-object-method">&nbsp;</a> 
```python
to_state(self, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L527)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L527?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.from_state" class="docs-object-method">&nbsp;</a> 
```python
from_state(data, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L533)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L533?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.indices" class="docs-object-method">&nbsp;</a> 
```python
@property
indices(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L?message=Update%20Docs)]
</div>

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
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L?message=Update%20Docs)]
</div>

Returns held excitations
- `inds`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.from_quanta" class="docs-object-method">&nbsp;</a> 
```python
from_quanta(basis, quants): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L588)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L588?message=Update%20Docs)]
</div>

Returns states with `quants` quanta of excitation
        using the basis `basis`
- `basis`: `RepresentationBasis`
    >No description...
- `quants`: `int | Iterable[int]`
    >set of integers
- `:returns`: `_`
    >BasisStateSpace

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.get_mode" class="docs-object-method">&nbsp;</a> 
```python
get_mode(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L611)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L611?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.infer_state_inds_type" class="docs-object-method">&nbsp;</a> 
```python
infer_state_inds_type(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L614)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L614?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.as_excitations" class="docs-object-method">&nbsp;</a> 
```python
as_excitations(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L625)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L625?message=Update%20Docs)]
</div>

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
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L659)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L659?message=Update%20Docs)]
</div>

Returns states as sets of excitations, rather than indices indo the basis functions.
        For 1D, this just means converting a list of states into tuples of length 1.
- `states`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.to_single" class="docs-object-method">&nbsp;</a> 
```python
to_single(self, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L691)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L691?message=Update%20Docs)]
</div>

Basically a no-op
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.is_unique" class="docs-object-method">&nbsp;</a> 
```python
is_unique(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L707)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L707?message=Update%20Docs)]
</div>

Returns `True` if the number of states is equal to number of unique states
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.is_sorted" class="docs-object-method">&nbsp;</a> 
```python
is_sorted(self, allow_indeterminate=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L718)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L718?message=Update%20Docs)]
</div>

Checks and then sets a flag
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.take_unique" class="docs-object-method">&nbsp;</a> 
```python
take_unique(self, sort=False, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L731)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L731?message=Update%20Docs)]
</div>

Returns only the unique states, but preserves
        ordering and all of that unless explicitly allowed not to
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.as_sorted" class="docs-object-method">&nbsp;</a> 
```python
as_sorted(self, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L774)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L774?message=Update%20Docs)]
</div>

Returns a sorted version of the state space
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.apply_selection_rules" class="docs-object-method">&nbsp;</a> 
```python
apply_selection_rules(self, selection_rules, target_dimensions=None, filter_space=None, parallelizer=None, logger=None, iterations=1, new_state_space_class=None, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L810)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L810?message=Update%20Docs)]
</div>

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
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L843)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L843?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.get_representation_indices" class="docs-object-method">&nbsp;</a> 
```python
get_representation_indices(self, other=None, selection_rules=None, freqs=None, freq_threshold=None, filter=None, return_filter=False, parallelizer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L846)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L846?message=Update%20Docs)]
</div>

Generates a set of indices that can be fed into a `Representation` to provide a sub-representation
        in this state space.
        Basically just takes all pairs of indices.
        Only returns the upper-triangle indices
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.get_representation_brakets" class="docs-object-method">&nbsp;</a> 
```python
get_representation_brakets(self, other=None, selection_rules=None, freqs=None, freq_threshold=None, filter=None, return_filter=False, track_excitations=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L929)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L929?message=Update%20Docs)]
</div>

Generates a `BraKetSpace` that can be fed into a `Representation`
        Only returns the upper-triangle pairs because we assume symmetry
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.take_subspace" class="docs-object-method">&nbsp;</a> 
```python
take_subspace(self, sel, assume_sorted=False, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L1002)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L1002?message=Update%20Docs)]
</div>

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
take_subdimensions(self, inds, exc=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L1059)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L1059?message=Update%20Docs)]
</div>

Returns a subsample of the space with some dimensions
        dropped
- `inds`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.take_states" class="docs-object-method">&nbsp;</a> 
```python
take_states(self, states, sort=False, assume_sorted=False, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L1077)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L1077?message=Update%20Docs)]
</div>

Takes the set of specified states from the space.
        A lot like take_subspace, but operates on states, not indices
- `states`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.drop_subspace" class="docs-object-method">&nbsp;</a> 
```python
drop_subspace(self, sel, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L1099)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L1099?message=Update%20Docs)]
</div>

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
drop_subdimensions(self, inds, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L1123)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L1123?message=Update%20Docs)]
</div>

Returns a subsample of the space with some dimensions
        dropped
- `inds`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.drop_states" class="docs-object-method">&nbsp;</a> 
```python
drop_states(self, states, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L1141)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L1141?message=Update%20Docs)]
</div>

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
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L1160)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L1160?message=Update%20Docs)]
</div>

Splits the space up into chunks of at max chunksize
- `chunksize`: `int`
    >No description...
- `:returns`: `Iterable[BasisStateSpace]`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.concatenate" class="docs-object-method">&nbsp;</a> 
```python
concatenate(self, other, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L1194)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L1194?message=Update%20Docs)]
</div>

Just does a direct concatenation with no unions or any
        of that
- `other`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.union" class="docs-object-method">&nbsp;</a> 
```python
union(self, other, sort=False, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L1253)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L1253?message=Update%20Docs)]
</div>

Returns a merged version of self and other, making
        use of as much of the information inherent in both as is possible
- `other`: `BasisStateSpace`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.intersection" class="docs-object-method">&nbsp;</a> 
```python
intersection(self, other, sort=False, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L1368)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L1368?message=Update%20Docs)]
</div>

Returns an intersected self and other
- `other`: `BasisStateSpace`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.difference" class="docs-object-method">&nbsp;</a> 
```python
difference(self, other, sort=False, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L1496)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L1496?message=Update%20Docs)]
</div>

Returns an diff'ed self and other
- `other`: `BasisStateSpace`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L1604)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L1604?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.__eq__" class="docs-object-method">&nbsp;</a> 
```python
__eq__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L1622)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L1622?message=Update%20Docs)]
</div>


- `other`: `Any`
    >No description...
- `:returns`: `_`
    >No description...



___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/ci/docs/Psience/BasisReps/StateSpaces/BasisStateSpace.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/ci/docs/Psience/BasisReps/StateSpaces/BasisStateSpace.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/ci/docs/Psience/BasisReps/StateSpaces/BasisStateSpace.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/ci/docs/Psience/BasisReps/StateSpaces/BasisStateSpace.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L470?message=Update%20Docs)