## <a id="Psience.BasisReps.StateSpaces.BraKetSpace">BraKetSpace</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces.py#L3640)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces.py#L3640?message=Update%20Docs)]
</div>

Represents a set of pairs of states that can be fed into a `Representation` or `Operator`
to efficiently tell it what terms it need to calculate.
This basically just implements a bunch of stuff for generating a Graph defining
the connections between states.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
aggressive_caching_enabled: bool
preindex_trie_enabled: bool
OrthogoIndexerTrie: OrthogoIndexerTrie
CachingOrthogonalIndexCalculator: CachingOrthogonalIndexCalculator
OrthogonalIndexSparseCalculator: OrthogonalIndexSparseCalculator
OrthogonalIndexCalculator: OrthogonalIndexCalculator
use_change_indices: bool
```
<a id="Psience.BasisReps.StateSpaces.BraKetSpace.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, bra_space, ket_space, changes=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L3647)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L3647?message=Update%20Docs)]
</div>

  - `bra_space`: `BasisStateSpace`
    > 
  - `ket_space`: `BasisStateSpace`
    >


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.state_pairs" class="docs-object-method">&nbsp;</a> 
```python
@property
state_pairs(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L3670)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L3670?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.from_indices" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_indices(cls, inds, basis=None, quanta=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L3685)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L3685?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.__len__" class="docs-object-method">&nbsp;</a> 
```python
__len__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L3717)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L3717?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L3725)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L3725?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.remove_duplicates" class="docs-object-method">&nbsp;</a> 
```python
remove_duplicates(self, assume_symmetric=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L3729)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L3729?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.load_space_diffs" class="docs-object-method">&nbsp;</a> 
```python
load_space_diffs(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L3760)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L3760?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.load_non_orthog" class="docs-object-method">&nbsp;</a> 
```python
load_non_orthog(self, use_aggressive_caching=None, use_preindex_trie=None, preindex_trie_depth=None, shared_memory_manager=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L3766)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L3766?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.share" class="docs-object-method">&nbsp;</a> 
```python
share(self, shared_memory_manager): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L4102)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L4102?message=Update%20Docs)]
</div>
Creates a shared memory version of the `BraKetSpace`
  - `shared_memory_manager`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.unshare" class="docs-object-method">&nbsp;</a> 
```python
unshare(self, shared_memory_manager): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L4123)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L4123?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.clear_cache" class="docs-object-method">&nbsp;</a> 
```python
clear_cache(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L4128)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L4128?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.free" class="docs-object-method">&nbsp;</a> 
```python
free(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L4133)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L4133?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.get_non_orthog" class="docs-object-method">&nbsp;</a> 
```python
get_non_orthog(self, inds, assume_unique=False, use_aggressive_caching=None, use_preindex_trie=None, preindex_trie_depth=None, shared_memory_manager=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L4169)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L4169?message=Update%20Docs)]
</div>
Returns whether the states are non-orthogonal under the set of indices.
  - `inds`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.get_sel_rules_from1d" class="docs-object-method">&nbsp;</a> 
```python
get_sel_rules_from1d(self, inds, rules): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L4199)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L4199?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.get_sel_rule_filter" class="docs-object-method">&nbsp;</a> 
```python
get_sel_rule_filter(self, rules): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L4221)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L4221?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.take_subspace" class="docs-object-method">&nbsp;</a> 
```python
take_subspace(self, sel): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L4246)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L4246?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.take_subdimensions" class="docs-object-method">&nbsp;</a> 
```python
take_subdimensions(self, inds): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L4263)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L4263?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.apply_non_orthogonality" class="docs-object-method">&nbsp;</a> 
```python
apply_non_orthogonality(self, inds, use_aggressive_caching=None, use_preindex_trie=None, preindex_trie_depth=None, assume_unique=False, use_change_indices=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L4275)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L4275?message=Update%20Docs)]
</div>
Takes the bra-ket pairs that are non-orthogonal under the indices `inds`
  - `inds`: `Any`
    > 
  - `assume_unique`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.apply_sel_rules_along" class="docs-object-method">&nbsp;</a> 
```python
apply_sel_rules_along(self, rules, inds, permute=True, dim=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L4391)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L4391?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.apply_sel_sums" class="docs-object-method">&nbsp;</a> 
```python
apply_sel_sums(self, rules, inds): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L4442)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L4442?message=Update%20Docs)]
</div>
We reckon it's fast enough to just determine if the number
of quanta in the bra is compatible with the number of
quanta in the ket...
  - `rules`: `Any`
    > 
  - `inds`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.apply_sel_rules" class="docs-object-method">&nbsp;</a> 
```python
apply_sel_rules(self, rules): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L4464)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L4464?message=Update%20Docs)]
</div>
Applies selections rules
  - `rules`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.adjacency_matrix" class="docs-object-method">&nbsp;</a> 
```python
adjacency_matrix(self, total_space=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L4476)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L4476?message=Update%20Docs)]
</div>
Generates the (sparse) unweighted adjacency matrix for the bras & kets
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.split" class="docs-object-method">&nbsp;</a> 
```python
split(self, chunksize): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L4505)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L4505?message=Update%20Docs)]
</div>
splits the brakets into blocks of at max chunksize
  - `chunksize`: `int`
    > 
  - `:returns`: `Iterable[BraKetSpace]`
    >


<a id="Psience.BasisReps.StateSpaces.BraKetSpace.concatenate" class="docs-object-method">&nbsp;</a> 
```python
concatenate(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/BraKetSpace.py#L4521)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/BraKetSpace.py#L4521?message=Update%20Docs)]
</div>
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/StateSpaces/BraKetSpace.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/StateSpaces/BraKetSpace.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/StateSpaces/BraKetSpace.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/StateSpaces/BraKetSpace.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces.py#L3640?message=Update%20Docs)   
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