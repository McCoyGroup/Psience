## <a id="Psience.BasisReps.StateSpaces.BraKetSpace">BraKetSpace</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L3212)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L3212?message=Update%20Docs)]
</div>

Represents a set of pairs of states that can be fed into a `Representation` or `Operator`
to efficiently tell it what terms it need to calculate.
This basically just implements a bunch of stuff for generating a Graph defining
the connections between states.

```python
aggressive_caching_enabled: bool
preindex_trie_enabled: bool
OrthogoIndexerTrie: type
CachingOrthogonalIndexCalculator: type
OrthogonalIndexSparseCalculator: type
OrthogonalIndexCalculator: type
```
<a id="Psience.BasisReps.StateSpaces.BraKetSpace.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, bra_space, ket_space): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L3219)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L3219?message=Update%20Docs)]
</div>


- `bra_space`: `BasisStateSpace`
    >No description...
- `ket_space`: `BasisStateSpace`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.state_pairs" class="docs-object-method">&nbsp;</a> 
```python
@property
state_pairs(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.from_indices" class="docs-object-method">&nbsp;</a> 
```python
from_indices(inds, basis=None, quanta=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L3254)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L3254?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.__len__" class="docs-object-method">&nbsp;</a> 
```python
__len__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L3286)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L3286?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L3294)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L3294?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.remove_duplicates" class="docs-object-method">&nbsp;</a> 
```python
remove_duplicates(self, assume_symmetric=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L3298)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L3298?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.load_space_diffs" class="docs-object-method">&nbsp;</a> 
```python
load_space_diffs(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L3329)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L3329?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.load_non_orthog" class="docs-object-method">&nbsp;</a> 
```python
load_non_orthog(self, use_aggressive_caching=None, use_preindex_trie=None, preindex_trie_depth=None, shared_memory_manager=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L3335)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L3335?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.share" class="docs-object-method">&nbsp;</a> 
```python
share(self, shared_memory_manager): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L3670)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L3670?message=Update%20Docs)]
</div>

Creates a shared memory version of the `BraKetSpace`
- `shared_memory_manager`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.unshare" class="docs-object-method">&nbsp;</a> 
```python
unshare(self, shared_memory_manager): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L3691)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L3691?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.clear_cache" class="docs-object-method">&nbsp;</a> 
```python
clear_cache(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L3696)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L3696?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.free" class="docs-object-method">&nbsp;</a> 
```python
free(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L3701)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L3701?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.get_non_orthog" class="docs-object-method">&nbsp;</a> 
```python
get_non_orthog(self, inds, assume_unique=False, use_aggressive_caching=None, use_preindex_trie=None, preindex_trie_depth=None, shared_memory_manager=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L3737)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L3737?message=Update%20Docs)]
</div>

Returns whether the states are non-orthogonal under the set of indices.
- `inds`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.get_sel_rules_from1d" class="docs-object-method">&nbsp;</a> 
```python
get_sel_rules_from1d(self, inds, rules): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L3767)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L3767?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.get_sel_rule_filter" class="docs-object-method">&nbsp;</a> 
```python
get_sel_rule_filter(self, rules): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L3789)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L3789?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.take_subspace" class="docs-object-method">&nbsp;</a> 
```python
take_subspace(self, sel): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L3814)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L3814?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.take_subdimensions" class="docs-object-method">&nbsp;</a> 
```python
take_subdimensions(self, inds): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L3826)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L3826?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.apply_non_orthogonality" class="docs-object-method">&nbsp;</a> 
```python
apply_non_orthogonality(self, inds, use_aggressive_caching=None, use_preindex_trie=None, preindex_trie_depth=None, assume_unique=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L3837)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L3837?message=Update%20Docs)]
</div>

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
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L3900)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L3900?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.apply_sel_sums" class="docs-object-method">&nbsp;</a> 
```python
apply_sel_sums(self, rules, inds): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L3951)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L3951?message=Update%20Docs)]
</div>

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
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L3973)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L3973?message=Update%20Docs)]
</div>

Applies selections rules
- `rules`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.adjacency_matrix" class="docs-object-method">&nbsp;</a> 
```python
adjacency_matrix(self, total_space=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L3985)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L3985?message=Update%20Docs)]
</div>

Generates the (sparse) unweighted adjacency matrix for the bras & kets
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.split" class="docs-object-method">&nbsp;</a> 
```python
split(self, chunksize): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L4014)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L4014?message=Update%20Docs)]
</div>

splits the brakets into blocks of at max chunksize
- `chunksize`: `int`
    >No description...
- `:returns`: `Iterable[BraKetSpace]`
    >No description...

<a id="Psience.BasisReps.StateSpaces.BraKetSpace.concatenate" class="docs-object-method">&nbsp;</a> 
```python
concatenate(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/StateSpaces.py#L4030)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L4030?message=Update%20Docs)]
</div>



___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/ci/docs/Psience/BasisReps/StateSpaces/BraKetSpace.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/ci/docs/Psience/BasisReps/StateSpaces/BraKetSpace.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/ci/docs/Psience/BasisReps/StateSpaces/BraKetSpace.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/ci/docs/Psience/BasisReps/StateSpaces/BraKetSpace.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/StateSpaces.py#L3212?message=Update%20Docs)