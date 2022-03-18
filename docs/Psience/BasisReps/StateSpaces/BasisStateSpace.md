## <a id="Psience.BasisReps.StateSpaces.BasisStateSpace">BasisStateSpace</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L470)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L470?message=Update%20Docs)]
</div>

Represents a subspace of states inside a representation basis.
Useful largely to provide consistent, unambiguous representations of multiple states across
the different representation-generating methods in the code base.

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, basis, states, full_basis=None, mode=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L476)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L476?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L517)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L517?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.to_state" class="docs-object-method">&nbsp;</a> 
```python
to_state(self, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L527)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L527?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.from_state" class="docs-object-method">&nbsp;</a> 
```python
from_state(data, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L533)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L533?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.indices" class="docs-object-method">&nbsp;</a> 
```python
@property
indices(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L588)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L588?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L611)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L611?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.infer_state_inds_type" class="docs-object-method">&nbsp;</a> 
```python
infer_state_inds_type(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L614)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L614?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.as_excitations" class="docs-object-method">&nbsp;</a> 
```python
as_excitations(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L625)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L625?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L659)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L659?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L691)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L691?message=Update%20Docs)]
</div>

Basically a no-op
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.is_unique" class="docs-object-method">&nbsp;</a> 
```python
is_unique(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L707)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L707?message=Update%20Docs)]
</div>

Returns `True` if the number of states is equal to number of unique states
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.is_sorted" class="docs-object-method">&nbsp;</a> 
```python
is_sorted(self, allow_indeterminate=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L718)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L718?message=Update%20Docs)]
</div>

Checks and then sets a flag
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.take_unique" class="docs-object-method">&nbsp;</a> 
```python
take_unique(self, sort=False, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L731)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L731?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L774)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L774?message=Update%20Docs)]
</div>

Returns a sorted version of the state space
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.apply_selection_rules" class="docs-object-method">&nbsp;</a> 
```python
apply_selection_rules(self, selection_rules, target_dimensions=None, filter_space=None, parallelizer=None, logger=None, iterations=1, new_state_space_class=None, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L810)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L810?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L843)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L843?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.get_representation_indices" class="docs-object-method">&nbsp;</a> 
```python
get_representation_indices(self, other=None, selection_rules=None, freqs=None, freq_threshold=None, filter=None, return_filter=False, parallelizer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L846)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L846?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L929)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L929?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L1002)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L1002?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L1059)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L1059?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L1077)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L1077?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L1099)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L1099?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L1123)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L1123?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L1141)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L1141?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L1160)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L1160?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L1194)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L1194?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L1253)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L1253?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L1368)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L1368?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L1496)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L1496?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L1604)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L1604?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.BasisStateSpace.__eq__" class="docs-object-method">&nbsp;</a> 
```python
__eq__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L1622)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L1622?message=Update%20Docs)]
</div>


- `other`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

 </div>
</div>



<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
### <a class="collapse-link" data-toggle="collapse" href="#tests">Tests</a> <a class="float-right" data-toggle="collapse" href="#tests"><i class="fa fa-chevron-down"></i></a>
 </div>
<div class="collapsible-section collapsible-section-body collapse show" id="tests" markdown="1">

- [HOBasis2DPP](#HOBasis2DPP)
- [HarmHam](#HarmHam)
- [HOBasis3DPXP](#HOBasis3DPXP)
- [HOBasis3DXXX2D](#HOBasis3DXXX2D)
- [HOBasis3DXXX2DContracted](#HOBasis3DXXX2DContracted)
- [HOSelRuleTerms](#HOSelRuleTerms)
- [GenerateSelectionRuleSpace](#GenerateSelectionRuleSpace)
- [GenerateFilteredSelectionRuleSpace](#GenerateFilteredSelectionRuleSpace)
- [StateIndexing](#StateIndexing)
- [FindIndices](#FindIndices)
- [PermIndexingChange](#PermIndexingChange)
- [NewOrthogonalityCalcs](#NewOrthogonalityCalcs)
- [StateSpaceIntersections](#StateSpaceIntersections)
- [BasisRepMatrixOps](#BasisRepMatrixOps)
- [OperatorAdjacencyGraph](#OperatorAdjacencyGraph)

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
#### <a class="collapse-link" data-toggle="collapse" href="#test-setup">Setup</a> <a class="float-right" data-toggle="collapse" href="#test-setup"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="test-setup" markdown="1">

Before we can run our examples we should get a bit of setup out of the way.
Since these examples were harvested from the unit tests not all pieces
will be necessary for all situations.
```python
from Peeves import Timer, BlockProfiler
from McUtils.Scaffolding import *
import McUtils.Plots as plt
from McUtils.Combinatorics import CompleteSymmetricGroupSpace
from Peeves.TestUtils import *
from unittest import TestCase
from Psience.BasisReps import *
import sys, os, numpy as np
```

All tests are wrapped in a test class
```python
class BasisSetTests(TestCase):
    def get_states(self, n_quanta, n_modes, max_quanta=None):
        return [np.flip(x) for x in BasisStateSpace.from_quanta(
            HarmonicOscillatorProductBasis(n_modes),
            range(n_quanta)
        ).excitations]
```

 </div>
</div>

#### <a name="HOBasis2DPP">HOBasis2DPP</a>
```python
    def test_HOBasis2DPP(self):
        from Peeves import Timer, BlockProfiler

        n = 10
        m = 2
        oppo = HarmonicOscillatorProductBasis((n,) * m)
        oppo2 = SimpleProductBasis(HarmonicOscillatorBasis, (n,) * m)

        term = ['p', 'p']
        iphase = (-1) ** (term.count("p") // 2)
        n_terms = len(term)

        g1 = np.array(
            [[-1.81146079e-04, 3.97836803e-05],
             [3.97836803e-05, 2.63572358e-05]])
        xxpp1 = 2 * oppo.representation(*term, coeffs=g1, axes=[[0, 1], [1, 0]])
        xxpp1 = xxpp1 + xxpp1
        xxpp2 = 2 * oppo2.representation(*term, coeffs=g1, axes=[[0, 1], [1, 0]])
        xxpp2 = xxpp2 + xxpp2

        usr = os.path.expanduser('~')
        job_is_dumb = [
            os.path.join(usr, "Documents/Python/config/python3.7/lib/python3.7/"),
            os.path.join(usr, "Documents/UW/Research/Development")
        ]

        quant_states = BasisStateSpace(
            oppo,
            self.get_states(9, m, max_quanta=10)
        )
        brakets = quant_states.get_representation_brakets()

        # with Timer("New style"):
            # with BlockProfiler("New Style", strip_dirs=job_is_dumb, num_lines=10, sort_by='tottime', filter="Psience"):
        vals1 = xxpp1[brakets]

        # with Timer("Old style"):
            # with BlockProfiler("Old Style", strip_dirs=job_is_dumb, num_lines=10, sort_by='tottime', filter="Psience"):
        vals2 = xxpp2[brakets]

        v1 = vals1
        v2 = iphase * vals2

        # n = len(quant_states)
        # plt.ArrayPlot(v1.reshape((n, n)))
        # plt.ArrayPlot(v2.reshape((n, n))).show()

        self.assertLess(np.max(np.abs(v1 - v2)), 1.0e-14)
```
#### <a name="HarmHam">HarmHam</a>
```python
    def test_HarmHam(self):

        n = 10
        m = 3
        basis = HarmonicOscillatorProductBasis((n,) * m)
        G, V = [
            np.array([[6.47886479e-03, 5.17641431e-12, -1.12922679e-12],
                      [5.17641431e-12, 1.28034398e-02, -3.15629792e-12],
                      [-1.12922679e-12, -3.15629792e-12, 1.76505371e-02]]),
            np.array([[6.47886478e-03, -8.45595180e-13, -1.01327126e-11],
                      [-8.45595549e-13, 1.28034398e-02, -4.72136245e-12],
                      [-1.01327124e-11, -4.72136255e-12, 1.76505372e-02]])]

        mommy = (1 / 2) * basis.representation('p', 'p', coeffs=G)
        possy = (1 / 2) * basis.representation('x', 'x', coeffs=V)
        H0 = ( mommy + possy )

        states = BasisStateSpace(basis, self.get_states(2, 3, max_quanta=10), mode='excitations')

        diag_inds = BraKetSpace(states, states)

        # raise Exception(diag_inds.state_pairs)

        diags = H0[diag_inds]

        self.assertEquals(np.average(diags), 0.036932841734999985)
```
#### <a name="HOBasis3DPXP">HOBasis3DPXP</a>
```python
    def test_HOBasis3DPXP(self):
        from Peeves import Timer, BlockProfiler

        n = 10
        m = 3
        oppo = HarmonicOscillatorProductBasis((n,) * m)
        oppo2 = SimpleProductBasis(HarmonicOscillatorBasis, (n,) * m)

        term = ['p', 'x', 'p']
        iphase = (-1) ** (term.count("p") // 2)
        n_terms = len(term)

        g1 = np.array([
            [[-1.81146079e-04,  3.97836803e-05,  2.91649691e-05],
             [ 3.97836803e-05,  2.63572358e-05,  2.37597837e-04],
             [ 2.91649691e-05,  2.37597837e-04, -3.38457268e-05]],

            [[-4.36589189e-04,  2.79004059e-05, -1.50059967e-05],
             [ 2.79004059e-05, -1.44188965e-06,  3.49657651e-06],
             [-1.50059967e-05,  3.49657652e-06,  3.11501367e-06]],

            [[-8.10821036e-04,  6.31615150e-06,  5.50255712e-05],
             [ 6.31615151e-06,  4.05569426e-06,  3.51303496e-08],
             [ 5.50255712e-05,  3.51303696e-08, -3.80070492e-06]]])
        xxpp1 = oppo.representation(*term, coeffs=g1, axes=[[0, 1, 2], [1, 0, 2]])
        # xxpp1 = xxpp1 + xxpp1
        xxpp2 = oppo2.representation(*term, coeffs=g1,  axes=[[0, 1, 2], [1, 0, 2]])
        # xxpp2 = xxpp2 + xxpp2

        quant_states = BasisStateSpace(
            oppo,
            self.get_states(3, 3, max_quanta=10)
        )
        inds = quant_states.get_representation_brakets()

        #
        # # raise Exception(inds.bras.indices)
        #
        # quant_states = BasisStateSpace(
        #     oppo,
        #     self.get_states(3, 3, max_quanta=10)
        # )
        # new_stuff = quant_states.apply_selection_rules([[-1, 1]])
        # inds2 = new_stuff.get_representation_brakets()
        #
        # plt.ArrayPlot(inds2.adjacency_matrix().toarray()).show()


        # inds = BasisStateSpace(
        #     oppo,
        #     (
        #         [0, 0, 0],
        #         [1, 0, 0],
        #     )
        # ).get_representation_brakets()


        # with Timer("New style"):
        vals1 = xxpp1[inds]
        # with Timer("Old style"):
        vals2 = xxpp2[inds]

        v1 = vals1
        v2 = iphase * vals2

        # n = len(quant_states)
        # plt.ArrayPlot(v1.reshape((n, n)))
        # plt.ArrayPlot(v2.reshape((n, n)))
        # plt.ArrayPlot(v1.reshape((n, n)) - v1.reshape((n, n)).T,
        #               plot_style=dict(vmin=-1.0e-5, vmax=1.0e-5))
        # plt.ArrayPlot(v2.reshape((n, n)) - v2.reshape((n, n)).T,
        #               plot_style=dict(vmin=-1.0e-5, vmax=1.0e-5))
        # plt.ArrayPlot((v1 - v2).reshape((n, n)),
        #               plot_style=dict(vmin=-1.0e-5, vmax=1.0e-5)).show()

        self.assertTrue(
            np.allclose(
                v1[:15],
                [0.00000000e+00, -2.86578374e-04, 0.00000000e+00, 3.29150701e-06,
                 -1.53766049e-04, 0.00000000e+00, -1.59263719e-06, 0.00000000e+00,
                 -5.52442364e-06, 1.24871307e-06, -6.66923918e-05, 0.00000000e+00,
                 -3.81027078e-05, 0.00000000e+00, -1.61862393e-04],
                atol=1.0e-5
            ))

        self.assertLess(np.max(np.abs(v1 - v2)), 1.0e-14)
```
#### <a name="HOBasis3DXXX2D">HOBasis3DXXX2D</a>
```python
    def test_HOBasis3DXXX2D(self):

        n = 15
        m = 2
        oppo = HarmonicOscillatorProductBasis((n,) * m)
        oppo2 = SimpleProductBasis(HarmonicOscillatorBasis, (n,) * m)

        term = ['x', 'x', 'x']
        iphase = (-1) ** (term.count("p") // 2)

        xxpp1 = oppo.representation(*term)
        xxpp2 = oppo2.representation(*term)

        states = BasisStateSpace.from_quanta(oppo, range(10))
        brakets = states.get_representation_brakets()
        vals1 = xxpp1[brakets]
        vals2 = xxpp2[brakets]

        v1 = vals1.asarray()
        v2 = iphase * vals2.asarray()

        # with JSONCheckpointer(os.path.expanduser("~/Desktop/test_terms.json")) as chk:
        #     chk['XXX_exc'] = states.excitations
        #     chk['XXX_3D_new'] = v1
        #     chk['XXX_3D_old'] = v2

        self.assertLess(np.max(np.abs(v1 - v2)), 2.0e-14)
```
#### <a name="HOBasis3DXXX2DContracted">HOBasis3DXXX2DContracted</a>
```python
    def test_HOBasis3DXXX2DContracted(self):
        n = 15
        m = 2
        oppo = HarmonicOscillatorProductBasis((n,) * m)
        oppo2 = SimpleProductBasis(HarmonicOscillatorBasis, (n,) * m)

        term = ['x', 'x', 'x']
        iphase = (-1) ** (term.count("p") // 2)

        xxpp1 = oppo.representation(*term, coeffs=np.ones((m, m, m)))
        xxpp2 = oppo2.representation(*term, coeffs=np.ones((m, m, m)))


        states = BasisStateSpace.from_quanta(oppo, range(10))
        brakets = states.get_representation_brakets()
        vals1 = xxpp1[brakets]
        vals2 = xxpp2[brakets]

        v1 = vals1
        v2 = iphase * vals2

        # with JSONCheckpointer(os.path.expanduser("~/Desktop/test_terms.json")) as chk:
        #     chk['XXX_exc'] = states.excitations
        #     chk['XXX_3D_new'] = v1
        #     chk['XXX_3D_old'] = v2

        self.assertLess(np.max(np.abs(v1 - v2)), 2.0e-14)
```
#### <a name="HOSelRuleTerms">HOSelRuleTerms</a>
```python
    def test_HOSelRuleTerms(self):
        """
        Profiler to see how quickly terms can be generated

        :return:
        :rtype:
        """

        n = 15
        m = 6
        basis = HarmonicOscillatorProductBasis((n,) * m)

        states = BasisStateSpace(
            basis,
            self.get_states(2, m)
        )

        transitions_h1 = [
            [-1],
            [1],
            [-3],
            [3],
            [-1, -1, -1],
            [-1, -1, 1],
            [-1, 1, 1],
            [1, 1, 1],
            [1, 2],
            [-1, 2],
            [1, -2],
            [-1, -2]
        ]

        with BlockProfiler("Selection Rules"):
            h1_space = states.apply_selection_rules(
                transitions_h1,
                1
            )
```
#### <a name="GenerateSelectionRuleSpace">GenerateSelectionRuleSpace</a>
```python
    def test_GenerateSelectionRuleSpace(self):
        """
        Tests (and profiles) the generation of a state
        space from a set of selection rules and initial states.
        Mostly here to more easily speed up state space generation
        for use in VPT2.

        :return:
        :rtype:
        """

        basis = HarmonicOscillatorProductBasis(8)
        rules = basis.selection_rules("x", "x", "x", "x")

        states = BasisStateSpace.from_quanta(basis, 3)

        # with BlockProfiler(""):
        h2_space = states.apply_selection_rules(rules, iterations=1)

        self.assertEquals(h2_space.nstates, 120)
```
#### <a name="GenerateFilteredSelectionRuleSpace">GenerateFilteredSelectionRuleSpace</a>
```python
    def test_GenerateFilteredSelectionRuleSpace(self):
        """
        Tests (and profiles) the generation of a state
        space from a set of selection rules and initial states.
        Mostly here to more easily speed up state space generation
        for use in VPT2.

        :return:
        :rtype:
        """

        basis = HarmonicOscillatorProductBasis(8)
        rules = basis.selection_rules("x", "x", "x", "x")

        states = BasisStateSpace.from_quanta(basis, 3)

        h2_space = states.apply_selection_rules(rules, iterations=1)

        sub_h2_space = h2_space.spaces[0].take_subspace(np.arange(10))

        h2_space2 = states.apply_selection_rules(rules,
                                                 filter_space=sub_h2_space,
                                                 iterations=1
                                                 )

        uinds, counts = np.unique(h2_space2.indices, return_counts=True)
        sorting = np.argsort(h2_space2.indices)
        ind_tag = (hash(tuple(uinds)), hash(tuple(counts)), hash(tuple(sorting)))
        # raise Exception(ind_tag)
        self.assertEquals(ind_tag, (320425735722628681, 4044592283957769633))
```
#### <a name="StateIndexing">StateIndexing</a>
```python
    def test_StateIndexing(self):
        """
        Tests indexing state specs through a more
        intelligent lexicographic order
        :return:
        :rtype:
        """

        ndim = 6
        indexer = PermutationStateIndexer(ndim)

        states = BasisStateSpace.from_quanta(HarmonicOscillatorProductBasis(ndim), range(10)).excitations
        # print(states)
        inds = indexer.to_indices(states)

        # print(states[44:])

        checks = inds != np.arange(len(states))
        self.assertFalse(
            checks.any()
            , msg="{} no good ({} out of {})".format(states[checks], inds[checks], inds)
        )

        # np.random.seed(0)
        # some_sel = np.arange(len(states))
        some_sel = np.unique(np.random.choice(np.arange(len(states)), 100))
        rev = indexer.from_indices(inds[some_sel,])
        self.assertTrue((states[some_sel,] == rev).all(),
                        msg="{} != {}".format(states[some_sel,], rev))
```
#### <a name="FindIndices">FindIndices</a>
```python
    def test_FindIndices(self):
        ndim = 6
        states = BasisStateSpace.from_quanta(HarmonicOscillatorProductBasis(ndim), range(5))
        test_1 = states.find(states)
        ntest = np.arange(len(test_1))
        self.assertEquals(tuple(test_1), tuple(ntest))

        sel = np.random.choice(ntest, 15)
        _, upos = np.unique(sel, return_index=True)
        sel = sel[np.sort(upos)]
        states2 = states.take_subspace(sel)
        test_2 = states2.find(states2)
        self.assertEquals(tuple(test_2), tuple(np.arange(len(sel))))
```
#### <a name="PermIndexingChange">PermIndexingChange</a>
```python
    def test_PermIndexingChange(self):
        import json
        ndim = 5
        basis = HarmonicOscillatorProductBasis(ndim, indexer=PermutationStateIndexer(ndim))
        rules = basis.selection_rules("x", "x", "x", "x")

        full_states = BasisStateSpace.from_quanta(basis, 1)
        for x in range(4):
            states = full_states.take_subspace([x])
            # if len(states) == 0:
            #     raise ValueError(states)

            # with BlockProfiler(""):
            h2_space = states.apply_selection_rules(rules, iterations=1)

            print(h2_space)

            print(states.indices.tolist(), np.sort(h2_space.indices).tolist())
```
#### <a name="NewOrthogonalityCalcs">NewOrthogonalityCalcs</a>
```python
    def test_NewOrthogonalityCalcs(self):

        n = 15
        m = 4
        oppo = HarmonicOscillatorProductBasis((n,) * m)

        states = BasisStateSpace.from_quanta(oppo, range(10))
        brakets = states.get_representation_brakets()
        orthog_1 = brakets.get_non_orthog([0, 0, 1])

        brakets2 = states.get_representation_brakets()
        brakets2.preindex_trie_enabled=False
        brakets2.aggressive_caching_enabled = False
        orthog_2 = brakets2.get_non_orthog([0, 0, 1])

        self.assertTrue( (orthog_1==orthog_2).all() )

        m = 2
        oppo = HarmonicOscillatorProductBasis((n,) * m)

        states = BasisStateSpace.from_quanta(oppo, range(10))

        brakets = states.get_representation_brakets()
        orthog_1 = brakets.get_non_orthog([0, 0, 1])

        brakets2 = states.get_representation_brakets()
        brakets2.preindex_trie_enabled = False
        brakets2.aggressive_caching_enabled = False
        orthog_2 = brakets2.get_non_orthog([0, 0, 1])

        # raise Exception(orthog_1, orthog_2)

        self.assertTrue((orthog_1 == orthog_2).all())
```
#### <a name="StateSpaceIntersections">StateSpaceIntersections</a>
```python
    def test_StateSpaceIntersections(self):

        basis = HarmonicOscillatorProductBasis(6)

        np.random.seed(0)
        subinds = np.random.random_integers(0, 100, 20)

        subspace = BasisStateSpace(basis, subinds)

        subsubinds = np.random.choice(subinds, 5, replace=False)
        filter_inds = np.unique(
            np.concatenate([
                subsubinds,
                np.random.random_integers(0, 100, 30)
            ])
        )

        filter_space = BasisStateSpace(basis, filter_inds)

        inter_space = subspace.intersection(filter_space)

        self.assertEquals(list(np.sort(inter_space.indices)), list(np.intersect1d(filter_inds, subinds)))

        subinds = np.array([12, 78, 0, 11, 10, 9])
        subspace = BasisStateSpace(basis, subinds, mode=BasisStateSpace.StateSpaceSpec.Indices)
        exc = subspace.excitations
        subspace = BasisStateSpace(basis, exc, mode=BasisStateSpace.StateSpaceSpec.Excitations)

        filter_inds = np.array([0, 6, 5, 4, 3, 2, 1])
        filter_space = BasisStateSpace(basis, filter_inds)

        inter_space = subspace.intersection(filter_space)

        self.assertEquals(
            list(np.sort(inter_space.indices)),
            list(np.intersect1d(filter_inds, subinds))
        )
```
#### <a name="BasisRepMatrixOps">BasisRepMatrixOps</a>
```python
    def test_BasisRepMatrixOps(self):

        n = 15 # totally meaningless these days
        m = 4
        basis = HarmonicOscillatorProductBasis((n,) * m)

        mat = StateSpaceMatrix(basis)

        self.assertEquals(mat.array.shape[0], 0)

        states = BasisStateSpace.from_quanta(basis, range(10))
        brakets = BraKetSpace(states, states)
        vals = np.ones(len(brakets))
        mat_2 = StateSpaceMatrix(brakets, vals)

        def wat(state_space):
            return np.ones(len(state_space))
        sub_brakets = BasisStateSpace.from_quanta(basis, range(4)).get_representation_brakets()
        mat2_vals = mat_2.compute_values(wat, sub_brakets)

        self.assertEquals(mat2_vals.tolist(), mat_2[sub_brakets].tolist())
```
#### <a name="OperatorAdjacencyGraph">OperatorAdjacencyGraph</a>
```python
    def test_OperatorAdjacencyGraph(self):
        """
        Tests building an adjacency graph for an operator
        under an initial set of states

        :return:
        :rtype:
        """

        from McUtils.Numputils import SparseArray

        ndim = 4
        basis = HarmonicOscillatorProductBasis(ndim)
        oppo = basis.representation("x", "x", "x", "x", coeffs=np.ones((ndim,)*4))
        rules = basis.selection_rules("x", "x", "x", "x")

        states = BasisStateSpace.from_quanta(basis, 3)
        h2_space = states.apply_selection_rules(rules, iterations=1)
        bk = h2_space.get_representation_brakets()
        flat_total_space = h2_space.to_single().take_unique()

        # pull from get_vpt2_reps or whatever
        # sub = oppo[bk]
        # flat_total_space = h2_space.to_single().take_unique()
        # N = len(flat_total_space)
        #
        # row_inds = flat_total_space.find(bk.bras)
        # col_inds = flat_total_space.find(bk.kets)
        #
        # up_tri = np.array([row_inds, col_inds]).T
        # low_tri = np.array([col_inds, row_inds]).T
        # # but now we need to remove the duplicates, because many sparse matrix implementations
        # # will sum up any repeated elements
        # full_inds = np.concatenate([up_tri, low_tri])
        # full_dat = np.concatenate([sub, sub])
        #
        # _, idx = np.unique(full_inds, axis=0, return_index=True)
        # sidx = np.sort(idx)
        # full_inds = full_inds[sidx]
        # full_dat = full_dat[sidx]
        # adj_mat = SparseArray((full_dat, full_inds.T), shape=(N, N))

        adj_arr = bk.adjacency_matrix(total_space=flat_total_space).toarray()
```

 </div>
</div>

___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/StateSpaces/BasisStateSpace.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/StateSpaces/BasisStateSpace.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/StateSpaces/BasisStateSpace.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/StateSpaces/BasisStateSpace.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L470?message=Update%20Docs)