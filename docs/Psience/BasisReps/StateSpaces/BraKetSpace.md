## <a id="Psience.BasisReps.StateSpaces.BraKetSpace">BraKetSpace</a>
Represents a set of pairs of states that can be fed into a `Representation` or `Operator`
to efficiently tell it what terms it need to calculate.
This basically just implements a bunch of stuff for generating a Graph defining
the connections between states.

### Properties and Methods
```python
OrthogoIndexerTrie: type
CachingOrthogonalIndexCalculator: type
OrthogonalIndexSparseCalculator: type
OrthogonalIndexCalculator: type
aggressive_caching_enabled: bool
preindex_trie_enabled: bool
```
<a id="Psience.BasisReps.StateSpaces.BraKetSpace.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, bra_space, ket_space): 
```

- `bra_space`: `BasisStateSpace`
    >No description...
- `ket_space`: `BasisStateSpace`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.state_pairs" class="docs-object-method">&nbsp;</a>
```python
@property
state_pairs(self): 
```

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.from_indices" class="docs-object-method">&nbsp;</a>
```python
from_indices(inds, basis=None, quanta=None): 
```

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.__len__" class="docs-object-method">&nbsp;</a>
```python
__len__(self): 
```

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.__repr__" class="docs-object-method">&nbsp;</a>
```python
__repr__(self): 
```

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.remove_duplicates" class="docs-object-method">&nbsp;</a>
```python
remove_duplicates(self, assume_symmetric=True): 
```

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.load_space_diffs" class="docs-object-method">&nbsp;</a>
```python
load_space_diffs(self): 
```

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.load_non_orthog" class="docs-object-method">&nbsp;</a>
```python
load_non_orthog(self, use_aggressive_caching=None, use_preindex_trie=None, preindex_trie_depth=None): 
```

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.clear_cache" class="docs-object-method">&nbsp;</a>
```python
clear_cache(self): 
```

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.get_non_orthog" class="docs-object-method">&nbsp;</a>
```python
get_non_orthog(self, inds, assume_unique=False, use_aggressive_caching=None, use_preindex_trie=None, preindex_trie_depth=None): 
```
Returns whether the states are non-orthogonal under the set of indices.
- `inds`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.get_sel_rules_from1d" class="docs-object-method">&nbsp;</a>
```python
get_sel_rules_from1d(self, inds, rules): 
```

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.get_sel_rule_filter" class="docs-object-method">&nbsp;</a>
```python
get_sel_rule_filter(self, rules): 
```

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.take_subspace" class="docs-object-method">&nbsp;</a>
```python
take_subspace(self, sel): 
```

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.take_subdimensions" class="docs-object-method">&nbsp;</a>
```python
take_subdimensions(self, inds): 
```

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.apply_non_orthogonality" class="docs-object-method">&nbsp;</a>
```python
apply_non_orthogonality(self, inds, use_aggressive_caching=None, use_preindex_trie=None, preindex_trie_depth=None, assume_unique=False): 
```
Takes the bra-ket pairs that are non-orthogonal under the indices `inds`
- `inds`: `Any`
    >No description...
- `assume_unique`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.apply_sel_rules_along" class="docs-object-method">&nbsp;</a>
```python
apply_sel_rules_along(self, rules, inds, permute=True, dim=None): 
```

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.apply_sel_sums" class="docs-object-method">&nbsp;</a>
```python
apply_sel_sums(self, rules, inds): 
```
We reckon it's fast enough to just determine if the number
        of quanta in the bra is compatible with the number of
        quanta in the ket...
- `rules`: `Any`
    >No description...
- `inds`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.apply_sel_rules" class="docs-object-method">&nbsp;</a>
```python
apply_sel_rules(self, rules): 
```
Applies selections rules
- `rules`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.adjacency_matrix" class="docs-object-method">&nbsp;</a>
```python
adjacency_matrix(self, total_space=None): 
```
Generates the (sparse) unweighted adjacency matrix for the bras & kets
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.split" class="docs-object-method">&nbsp;</a>
```python
split(self, chunksize): 
```
splits the brakets into blocks of at max chunksize
- `chunksize`: `int`
    >No description...
- `:returns`: `Iterable[BraKetSpace]`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.concatenate" class="docs-object-method">&nbsp;</a>
```python
concatenate(self, other): 
```

### Examples




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/edit/ci/examples/ci/docs/Psience/BasisReps/StateSpaces/BraKetSpace.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/examples/ci/docs/Psience/BasisReps/StateSpaces/BraKetSpace.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/edit/ci/docs/ci/docs/Psience/BasisReps/StateSpaces/BraKetSpace.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/docs/templates/ci/docs/Psience/BasisReps/StateSpaces/BraKetSpace.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py?message=Update%20Docs)