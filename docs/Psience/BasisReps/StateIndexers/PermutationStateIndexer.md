## <a id="Psience.BasisReps.StateIndexers.PermutationStateIndexer">PermutationStateIndexer</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateIndexers.py#L90)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateIndexers.py#L90?message=Update%20Docs)]
</div>

A sophisticated indexer that takes a state dimension and provides
indices based on the `shortlex` ordering, where ordering is defined
first by # of quanta of excitation, then by which partitioning of the #quanta,
 it represents, and finally by which permutation of that paritioning it is.
Is likely about as stable as an indexer can be expected to be over large
numbers of states. Unlikely to exhaust the max integers available for most
systems.

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

<a id="Psience.BasisReps.StateIndexers.PermutationStateIndexer.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, ndim): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateIndexers.py#L100)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateIndexers.py#L100?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateIndexers.PermutationStateIndexer.to_state" class="docs-object-method">&nbsp;</a> 
```python
to_state(self, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateIndexers.py#L105)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateIndexers.py#L105?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateIndexers.PermutationStateIndexer.from_state" class="docs-object-method">&nbsp;</a> 
```python
from_state(data, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateIndexers.py#L109)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateIndexers.py#L109?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateIndexers.PermutationStateIndexer.to_indices" class="docs-object-method">&nbsp;</a> 
```python
to_indices(self, states): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateIndexers.py#L113)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateIndexers.py#L113?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateIndexers.py#L124)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateIndexers.py#L124?message=Update%20Docs)]
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

 </div>
</div>



<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
### <a class="collapse-link" data-toggle="collapse" href="#tests">Tests</a> <a class="float-right" data-toggle="collapse" href="#tests"><i class="fa fa-chevron-down"></i></a>
 </div>
<div class="collapsible-section collapsible-section-body collapse show" id="tests" markdown="1">

- [StateIndexing](#StateIndexing)
- [PermIndices](#PermIndices)

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
#### <a name="PermIndices">PermIndices</a>
```python
    def test_PermIndices(self):
        ndim = 3
        indexer = PermutationStateIndexer(ndim)

        states = [[0, 0, 0],
         [0, 0, 1],
         [0, 0, 2],
         [0, 0, 3],
         [0, 0, 4],
         [0, 0, 5],
         [0, 0, 6],
         [0, 1, 0],
         [0, 1, 1],
         [0, 1, 2],
         [0, 2, 3],
         [0, 2, 4],
         [3, 1, 0],
         [3, 1, 1],
         [3, 1, 2],
         [3, 2, 0],
         [3, 2, 1],
         [3, 3, 0],
         [4, 0, 0],
         [4, 0, 1],
         [4, 0, 2],
         [4, 1, 0],
         [4, 1, 1],
         [4, 2, 0],
         [5, 0, 0],
         [5, 0, 1],
         [5, 1, 0],
         [6, 0, 0]]
        raise Exception(indexer.to_indices(states))
```

 </div>
</div>

___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/StateIndexers/PermutationStateIndexer.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/StateIndexers/PermutationStateIndexer.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/StateIndexers/PermutationStateIndexer.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/StateIndexers/PermutationStateIndexer.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateIndexers.py#L90?message=Update%20Docs)