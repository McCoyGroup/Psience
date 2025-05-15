## <a id="Psience.BasisReps.StateIndexers.PermutationStateIndexer">PermutationStateIndexer</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateIndexers.py#L90)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateIndexers.py#L90?message=Update%20Docs)]
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
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.BasisReps.StateIndexers.PermutationStateIndexer.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, ndim): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateIndexers/PermutationStateIndexer.py#L100)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateIndexers/PermutationStateIndexer.py#L100?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateIndexers.PermutationStateIndexer.to_state" class="docs-object-method">&nbsp;</a> 
```python
to_state(self, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateIndexers/PermutationStateIndexer.py#L105)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateIndexers/PermutationStateIndexer.py#L105?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateIndexers.PermutationStateIndexer.from_state" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_state(cls, data, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L109)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L109?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateIndexers.PermutationStateIndexer.to_indices" class="docs-object-method">&nbsp;</a> 
```python
to_indices(self, states): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateIndexers/PermutationStateIndexer.py#L113)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateIndexers/PermutationStateIndexer.py#L113?message=Update%20Docs)]
</div>
Finds the appropriate integer partitioning for each state
  - `states`: `np.ndarray`
    > 2D array of states as excitations
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateIndexers.PermutationStateIndexer.from_indices" class="docs-object-method">&nbsp;</a> 
```python
from_indices(self, indices, check=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateIndexers/PermutationStateIndexer.py#L124)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateIndexers/PermutationStateIndexer.py#L124?message=Update%20Docs)]
</div>
Inverts the index calculation.
First determines what number of quanta the index corresponds to,
then which integer partition, and finally just loops through the unique
permutations of the partition to get the right one.
This is not assured to be a fast process in any way.
  - `indices`: `Iterable[int]`
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/StateIndexers/PermutationStateIndexer.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/StateIndexers/PermutationStateIndexer.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/StateIndexers/PermutationStateIndexer.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/StateIndexers/PermutationStateIndexer.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateIndexers.py#L90?message=Update%20Docs)   
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