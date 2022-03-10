## <a id="Psience.BasisReps.StateIndexers.PermutationStateIndexer">PermutationStateIndexer</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateIndexers.py#L90)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateIndexers.py#L90?message=Update%20Docs)]
</div>

A sophisticated indexer that takes a state dimension and provides
indices based on the `shortlex` ordering, where ordering is defined
first by # of quanta of excitation, then by which partitioning of the #quanta,
 it represents, and finally by which permutation of that paritioning it is.
Is likely about as stable as an indexer can be expected to be over large
numbers of states. Unlikely to exhaust the max integers available for most
systems.

<a id="Psience.BasisReps.StateIndexers.PermutationStateIndexer.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, ndim): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateIndexers.py#L100)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateIndexers.py#L100?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateIndexers.PermutationStateIndexer.to_state" class="docs-object-method">&nbsp;</a> 
```python
to_state(self, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateIndexers.py#L105)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateIndexers.py#L105?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateIndexers.PermutationStateIndexer.from_state" class="docs-object-method">&nbsp;</a> 
```python
from_state(data, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateIndexers.py#L109)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateIndexers.py#L109?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateIndexers.PermutationStateIndexer.to_indices" class="docs-object-method">&nbsp;</a> 
```python
to_indices(self, states): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateIndexers.py#L113)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateIndexers.py#L113?message=Update%20Docs)]
</div>

Finds the appropriate integer partitioning for each state
- `states`: `np.ndarray`
    >2D array of states as excitations
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateIndexers.PermutationStateIndexer.from_indices" class="docs-object-method">&nbsp;</a> 
```python
from_indices(self, indices, check=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateIndexers.py#L124)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateIndexers.py#L124?message=Update%20Docs)]
</div>

Inverts the index calculation.
        First determines what number of quanta the index corresponds to,
        then which integer partition, and finally just loops through the unique
        permutations of the partition to get the right one.
        This is not assured to be a fast process in any way.
- `indices`: `Iterable[int]`
    >No description...
- `:returns`: `_`
    >No description...



___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/ci/docs/Psience/BasisReps/StateIndexers/PermutationStateIndexer.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/ci/docs/Psience/BasisReps/StateIndexers/PermutationStateIndexer.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/ci/docs/Psience/BasisReps/StateIndexers/PermutationStateIndexer.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/ci/docs/Psience/BasisReps/StateIndexers/PermutationStateIndexer.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateIndexers.py#L90?message=Update%20Docs)