## <a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace">PermutationallyReducedStateSpace</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L1664)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L1664?message=Update%20Docs)]
</div>

Defines a basis state space where terms are reduced over their
permutationally equivalent operations, making many operations
dramatically faster

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

<a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, basis, class_reps, perms): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L1671)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L1671?message=Update%20Docs)]
</div>


- `original_space`: `BasisStateSpace`
    >No description...

<a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace.to_equivalence_class_space" class="docs-object-method">&nbsp;</a> 
```python
to_equivalence_class_space(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L1679)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L1679?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace.from_space" class="docs-object-method">&nbsp;</a> 
```python
from_space(original_space): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L1682)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L1682?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace.get_equivalent_permutations" class="docs-object-method">&nbsp;</a> 
```python
get_equivalent_permutations(exc): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L1687)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L1687?message=Update%20Docs)]
</div>


- `exc`: `np.ndarray`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace.permutation_direct_product" class="docs-object-method">&nbsp;</a> 
```python
permutation_direct_product(self, perms): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L1710)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L1710?message=Update%20Docs)]
</div>

Creates a new space by taking permutation products
- `perms`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace.apply_selection_rules" class="docs-object-method">&nbsp;</a> 
```python
apply_selection_rules(self, selection_rules, target_dimensions=None, filter_space=None, parallelizer=None, logger=None, iterations=1, new_state_space_class=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L1730)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L1730?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace.permutationally_reduce" class="docs-object-method">&nbsp;</a> 
```python
permutationally_reduce(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L1750)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L1750?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace.representative_space" class="docs-object-method">&nbsp;</a> 
```python
representative_space(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L1752)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L1752?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace.permutationally_expand" class="docs-object-method">&nbsp;</a> 
```python
permutationally_expand(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L1754)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L1754?message=Update%20Docs)]
</div>


- `:returns`: `BasisStateSpace`
    >No description...

<a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace.take_permutations" class="docs-object-method">&nbsp;</a> 
```python
take_permutations(self, *p): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L1765)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L1765?message=Update%20Docs)]
</div>

Takes subsets of the stored permutations.
        This function is subject to change as the held structure of the permutations
        changes.
        Since permutation structure is stored like a direct product to maintain equivalence
        class relations we index from the bottom out, i.e. asking for `take_permutations(i, j)`
        will give you the states where the original state was `i` and the first product was in `j`
- `p`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace.take_subspace" class="docs-object-method">&nbsp;</a> 
```python
take_subspace(self, sel, assume_sorted=False, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L1781)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L1781?message=Update%20Docs)]
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

<a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace.take_subdimensions" class="docs-object-method">&nbsp;</a> 
```python
take_subdimensions(self, inds): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L1806)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L1806?message=Update%20Docs)]
</div>

Returns a subsample of the space with some dimensions
        dropped
- `inds`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

 </div>
</div>




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/StateSpaces/PermutationallyReducedStateSpace.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/StateSpaces/PermutationallyReducedStateSpace.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/StateSpaces/PermutationallyReducedStateSpace.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/StateSpaces/PermutationallyReducedStateSpace.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L1664?message=Update%20Docs)