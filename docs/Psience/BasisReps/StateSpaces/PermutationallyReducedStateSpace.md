## <a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace">PermutationallyReducedStateSpace</a>
Defines a basis state space where terms are reduced over their
permutationally equivalent operations, making many operations
dramatically faster

### Properties and Methods
```python
from_space: method
get_equivalent_permutations: method
```
<a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, basis, class_reps, perms): 
```

- `original_space`: `BasisStateSpace`
    >No description...

<a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace.to_equivalence_class_space" class="docs-object-method">&nbsp;</a>
```python
to_equivalence_class_space(self): 
```

<a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace.permutation_direct_product" class="docs-object-method">&nbsp;</a>
```python
permutation_direct_product(self, perms): 
```
Creates a new space by taking permutation products
- `perms`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace.apply_selection_rules" class="docs-object-method">&nbsp;</a>
```python
apply_selection_rules(self, selection_rules, target_dimensions=None, filter_space=None, parallelizer=None, logger=None, iterations=1, new_state_space_class=None): 
```

<a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace.permutationally_reduce" class="docs-object-method">&nbsp;</a>
```python
permutationally_reduce(self): 
```

<a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace.representative_space" class="docs-object-method">&nbsp;</a>
```python
representative_space(self): 
```

<a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace.permutationally_expand" class="docs-object-method">&nbsp;</a>
```python
permutationally_expand(self): 
```

- `:returns`: `BasisStateSpace`
    >No description...

<a id="Psience.BasisReps.StateSpaces.PermutationallyReducedStateSpace.take_permutations" class="docs-object-method">&nbsp;</a>
```python
take_permutations(self, *p): 
```
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
Returns a subsample of the space with some dimensions
        dropped
- `inds`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

### Examples


___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/edit/ci/examples/ci/docs/Psience/BasisReps/StateSpaces/PermutationallyReducedStateSpace.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/examples/ci/docs/Psience/BasisReps/StateSpaces/PermutationallyReducedStateSpace.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/edit/ci/docs/ci/docs/Psience/BasisReps/StateSpaces/PermutationallyReducedStateSpace.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/docs/templates/ci/docs/Psience/BasisReps/StateSpaces/PermutationallyReducedStateSpace.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py?message=Update%20Docs)