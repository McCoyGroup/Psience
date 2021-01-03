## <a id="Psience.BasisReps.StateIndexers.PermutationStateIndexer">PermutationStateIndexer</a>
A sophisticated indexer that takes a state dimension and provides
indices based on the `shortlex` ordering, where ordering is defined
first by # of quanta of excitation, then by which partitioning of the #quanta,
 it represents, and finally by which permutation of that paritioning it is.
Is likely about as stable as an indexer can be expected to be over large
numbers of states. Unlikely to exhaust the max integers available for most
systems.

### Properties and Methods
```python
PartitionPermutationIndexer: type
```
<a id="Psience.BasisReps.StateIndexers.PermutationStateIndexer.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, ndim): 
```

<a id="Psience.BasisReps.StateIndexers.PermutationStateIndexer.integer_partitions" class="docs-object-method">&nbsp;</a>
```python
integer_partitions(self, num): 
```
Gives basically the second sort criterion by calculating integer partitions
        (i.e. number of way to split up num quanta across modes)
        Sorted by default, which is our saving grace
- `num`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateIndexers.PermutationStateIndexer.to_indices" class="docs-object-method">&nbsp;</a>
```python
to_indices(self, states): 
```
Finds the appropriate integer partitioning for each state
- `states`: `np.ndarray`
    >2D array of states as excitations
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateIndexers.PermutationStateIndexer.from_indices" class="docs-object-method">&nbsp;</a>
```python
from_indices(self, indices): 
```
Inverts the index calculation...not implemented yet
- `indices`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

### Examples


