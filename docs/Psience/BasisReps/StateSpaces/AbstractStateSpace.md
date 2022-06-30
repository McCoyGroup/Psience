## <a id="Psience.BasisReps.StateSpaces.AbstractStateSpace">AbstractStateSpace</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L25)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L25?message=Update%20Docs)]
</div>

Represents a generalized state space which will provide core
methods to index into a basis and generate representations

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

```python
keep_excitations: bool
keep_indices: bool
StateSpaceSpec: EnumMeta
StateSpaceCache: type
excitations_dtype: dtype[int8]
indices_dtype: dtype[uint64]
```
<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, basis): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L61)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L61?message=Update%20Docs)]
</div>


- `basis`: `RepresentationBasis`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.to_state" class="docs-object-method">&nbsp;</a> 
```python
to_state(self, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L77)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L77?message=Update%20Docs)]
</div>

Provides just the state that is needed to
        serialize the object
- `serializer`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.from_state" class="docs-object-method">&nbsp;</a> 
```python
from_state(state, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L90)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L90?message=Update%20Docs)]
</div>

Loads from the stored state
- `serializer`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.ndim" class="docs-object-method">&nbsp;</a> 
```python
@property
ndim(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.excitations" class="docs-object-method">&nbsp;</a> 
```python
@property
excitations(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.mode" class="docs-object-method">&nbsp;</a> 
```python
@property
mode(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.get_mode" class="docs-object-method">&nbsp;</a> 
```python
get_mode(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L129)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L129?message=Update%20Docs)]
</div>

Returns the mode (indices or excitations) for the held states
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.has_indices" class="docs-object-method">&nbsp;</a> 
```python
@property
has_indices(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.has_excitations" class="docs-object-method">&nbsp;</a> 
```python
@property
has_excitations(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.indices" class="docs-object-method">&nbsp;</a> 
```python
@property
indices(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.indexer" class="docs-object-method">&nbsp;</a> 
```python
@property
indexer(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.exc_indexer" class="docs-object-method">&nbsp;</a> 
```python
@property
exc_indexer(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.find" class="docs-object-method">&nbsp;</a> 
```python
find(self, to_search, check=True, minimal_dtype=False, dtype=None, missing_val='raise'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L187)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L187?message=Update%20Docs)]
</div>

Finds the indices of a set of indices inside the space
- `to_search`: `np.ndarray | AbstractStateSpace`
    >array of ints
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.__len__" class="docs-object-method">&nbsp;</a> 
```python
__len__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L221)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L221?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.unique_len" class="docs-object-method">&nbsp;</a> 
```python
@property
unique_len(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.unique_indices" class="docs-object-method">&nbsp;</a> 
```python
@property
unique_indices(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L?message=Update%20Docs)]
</div>

Returns the unique indices
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.unique_excitations" class="docs-object-method">&nbsp;</a> 
```python
@property
unique_excitations(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L?message=Update%20Docs)]
</div>

Returns the unique excitations
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.as_indices" class="docs-object-method">&nbsp;</a> 
```python
as_indices(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L252)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L252?message=Update%20Docs)]
</div>

Returns the index version of the stored states
- `:returns`: `np.ndarray`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.as_unique_indices" class="docs-object-method">&nbsp;</a> 
```python
as_unique_indices(self, sort=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L261)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L261?message=Update%20Docs)]
</div>

Returns unique indices
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.as_excitations" class="docs-object-method">&nbsp;</a> 
```python
as_excitations(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L285)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L285?message=Update%20Docs)]
</div>

Returns the excitation version of the stored states
- `:returns`: `np.ndarray`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.as_unique_excitations" class="docs-object-method">&nbsp;</a> 
```python
as_unique_excitations(self, sort=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L294)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L294?message=Update%20Docs)]
</div>

Returns unique excitations
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.get_representation_indices" class="docs-object-method">&nbsp;</a> 
```python
get_representation_indices(self, other=None, selection_rules=None, freqs=None, freq_threshold=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L314)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L314?message=Update%20Docs)]
</div>

Returns bra and ket indices that can be used as indices to generate representations
- `other`: `Any`
    >No description...
- `selection_rules`: `Any`
    >No description...
- `freqs`: `Any`
    >No description...
- `freq_threshold`: `Any`
    >No description...
- `:returns`: `(np.ndarray, np.ndarray)`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.get_representation_brakets" class="docs-object-method">&nbsp;</a> 
```python
get_representation_brakets(self, other=None, selection_rules=None, freqs=None, freq_threshold=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L338)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L338?message=Update%20Docs)]
</div>

Returns a BraKetSpace that can be used as generate representations
- `other`: `Any`
    >No description...
- `selection_rules`: `Any`
    >No description...
- `freqs`: `Any`
    >No description...
- `freq_threshold`: `Any`
    >No description...
- `:returns`: `BraKetSpace`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.take_states" class="docs-object-method">&nbsp;</a> 
```python
take_states(self, states): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L361)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L361?message=Update%20Docs)]
</div>

Takes the intersection of self and the specified states
- `states`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.take_subspace" class="docs-object-method">&nbsp;</a> 
```python
take_subspace(self, sel): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L371)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L371?message=Update%20Docs)]
</div>

Takes a subset of the states
- `sel`: `Any`
    >No description...
- `:returns`: `AbstractStateSpace`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.take_subdimensions" class="docs-object-method">&nbsp;</a> 
```python
take_subdimensions(self, inds): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L382)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L382?message=Update%20Docs)]
</div>

Takes a subset of the state dimensions
- `sel`: `Any`
    >No description...
- `:returns`: `AbstractStateSpace`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.drop_states" class="docs-object-method">&nbsp;</a> 
```python
drop_states(self, states): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L394)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L394?message=Update%20Docs)]
</div>

Takes the difference of self and the specified states
- `states`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.drop_subspace" class="docs-object-method">&nbsp;</a> 
```python
drop_subspace(self, sel): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L404)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L404?message=Update%20Docs)]
</div>

Drops a subset of the states
- `sel`: `Any`
    >No description...
- `:returns`: `AbstractStateSpace`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.drop_subdimensions" class="docs-object-method">&nbsp;</a> 
```python
drop_subdimensions(self, inds): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L415)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L415?message=Update%20Docs)]
</div>

Drops a subset of the state dimensions
- `sel`: `Any`
    >No description...
- `:returns`: `AbstractStateSpace`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.get_states_with_quanta" class="docs-object-method">&nbsp;</a> 
```python
get_states_with_quanta(n, ndim): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L427)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L427?message=Update%20Docs)]
</div>

Returns the states with number of quanta equal to n
- `quanta`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.to_single" class="docs-object-method">&nbsp;</a> 
```python
to_single(self, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L449)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L449?message=Update%20Docs)]
</div>

Flattens any complicated state space structure into a
        single space like a `BasisStateSpace`
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.split" class="docs-object-method">&nbsp;</a> 
```python
split(self, chunksize): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L463)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L463?message=Update%20Docs)]
</div>

Subclass overridable function to allow for spaces to be
        split up into chunks
- `chunksize`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.share" class="docs-object-method">&nbsp;</a> 
```python
share(self, shared_memory_manager): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L474)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L474?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.unshare" class="docs-object-method">&nbsp;</a> 
```python
unshare(self, shared_memory_manager): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/StateSpaces.py#L476)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L476?message=Update%20Docs)]
</div>

 </div>
</div>




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/StateSpaces/AbstractStateSpace.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/StateSpaces/AbstractStateSpace.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/StateSpaces/AbstractStateSpace.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/StateSpaces/AbstractStateSpace.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/StateSpaces.py#L25?message=Update%20Docs)