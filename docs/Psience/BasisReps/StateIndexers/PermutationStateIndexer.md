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
max_quants: int
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
- `:returns`: `Tuple[int, int, Iterable[PermutationStateIndexer.PartitionPermutationIndexer]]`
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
from_indices(self, indices, check=False): 
```
Inverts the index calculation.
        First determines what number of quanta the index corresponds to,
        then which integer partition, and finally just loops through the unique
        permutations of the partition to get the right one.
        This is not assured to be a fast process in any way.
- `indices`: `Iterable[int]`
    >No description...
- `:returns`: `_`
    >No description...

### Examples


