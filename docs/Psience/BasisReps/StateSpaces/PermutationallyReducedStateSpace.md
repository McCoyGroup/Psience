## <a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace">PermutationallyReducedStateSpace</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces.py#L1686)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces.py#L1686?message=Update%20Docs)]
</div>

Defines a basis state space where terms are reduced over their
permutationally equivalent operations, making many operations
dramatically faster







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, basis, class_reps, perms): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/PermutationallyReducedStateSpace.py#L1693)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/PermutationallyReducedStateSpace.py#L1693?message=Update%20Docs)]
</div>

  - `original_space`: `BasisStateSpace`
    >


<a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace.to_equivalence_class_space" class="docs-object-method">&nbsp;</a> 
```python
to_equivalence_class_space(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/PermutationallyReducedStateSpace.py#L1701)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/PermutationallyReducedStateSpace.py#L1701?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace.from_space" class="docs-object-method">&nbsp;</a> 
```python
from_space(original_space): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/PermutationallyReducedStateSpace.py#L1704)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/PermutationallyReducedStateSpace.py#L1704?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace.get_equivalent_permutations" class="docs-object-method">&nbsp;</a> 
```python
get_equivalent_permutations(exc): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/PermutationallyReducedStateSpace.py#L1709)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/PermutationallyReducedStateSpace.py#L1709?message=Update%20Docs)]
</div>

  - `exc`: `np.ndarray`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace.permutation_direct_product" class="docs-object-method">&nbsp;</a> 
```python
permutation_direct_product(self, perms): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/PermutationallyReducedStateSpace.py#L1732)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/PermutationallyReducedStateSpace.py#L1732?message=Update%20Docs)]
</div>
Creates a new space by taking permutation products
  - `perms`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace.apply_selection_rules" class="docs-object-method">&nbsp;</a> 
```python
apply_selection_rules(self, selection_rules, target_dimensions=None, filter_space=None, parallelizer=None, logger=None, iterations=1, new_state_space_class=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/PermutationallyReducedStateSpace.py#L1752)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/PermutationallyReducedStateSpace.py#L1752?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace.permutationally_reduce" class="docs-object-method">&nbsp;</a> 
```python
permutationally_reduce(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/PermutationallyReducedStateSpace.py#L1772)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/PermutationallyReducedStateSpace.py#L1772?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace.representative_space" class="docs-object-method">&nbsp;</a> 
```python
representative_space(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/PermutationallyReducedStateSpace.py#L1774)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/PermutationallyReducedStateSpace.py#L1774?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace.permutationally_expand" class="docs-object-method">&nbsp;</a> 
```python
permutationally_expand(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/PermutationallyReducedStateSpace.py#L1776)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/PermutationallyReducedStateSpace.py#L1776?message=Update%20Docs)]
</div>

  - `:returns`: `BasisStateSpace`
    >


<a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace.take_permutations" class="docs-object-method">&nbsp;</a> 
```python
take_permutations(self, *p): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/PermutationallyReducedStateSpace.py#L1787)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/PermutationallyReducedStateSpace.py#L1787?message=Update%20Docs)]
</div>
Takes subsets of the stored permutations.
This function is subject to change as the held structure of the permutations
changes.
Since permutation structure is stored like a direct product to maintain equivalence
class relations we index from the bottom out, i.e. asking for `take_permutations(i, j)`
will give you the states where the original state was `i` and the first product was in `j`
  - `p`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace.take_subspace" class="docs-object-method">&nbsp;</a> 
```python
take_subspace(self, sel, assume_sorted=False, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/PermutationallyReducedStateSpace.py#L1803)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/PermutationallyReducedStateSpace.py#L1803?message=Update%20Docs)]
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


<a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace.take_subdimensions" class="docs-object-method">&nbsp;</a> 
```python
take_subdimensions(self, inds): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/PermutationallyReducedStateSpace.py#L1828)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/PermutationallyReducedStateSpace.py#L1828?message=Update%20Docs)]
</div>
Returns a subsample of the space with some dimensions
dropped
  - `inds`: `Any`
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/StateSpaces/PermutationallyReducedStateSpace.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/StateSpaces/PermutationallyReducedStateSpace.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/StateSpaces/PermutationallyReducedStateSpace.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/StateSpaces/PermutationallyReducedStateSpace.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces.py#L1686?message=Update%20Docs)   
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