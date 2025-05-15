## <a id="Psience.BasisReps.StateSpaces.AbstractStateSpace">AbstractStateSpace</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces.py#L26)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces.py#L26?message=Update%20Docs)]
</div>

Represents a generalized state space which will provide core
methods to index into a basis and generate representations







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
keep_excitations: bool
keep_indices: bool
StateSpaceSpec: StateSpaceSpec
StateSpaceCache: StateSpaceCache
excitations_dtype: Int8DType
indices_dtype: UInt64DType
```
<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, basis): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/AbstractStateSpace.py#L62)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/AbstractStateSpace.py#L62?message=Update%20Docs)]
</div>

  - `basis`: `RepresentationBasis`
    >


<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.to_state" class="docs-object-method">&nbsp;</a> 
```python
to_state(self, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/AbstractStateSpace.py#L78)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/AbstractStateSpace.py#L78?message=Update%20Docs)]
</div>
Provides just the state that is needed to
serialize the object
  - `serializer`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.from_state" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_state(cls, state, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L91)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L91?message=Update%20Docs)]
</div>
Loads from the stored state
  - `serializer`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.ndim" class="docs-object-method">&nbsp;</a> 
```python
@property
ndim(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/AbstractStateSpace.py#L104)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/AbstractStateSpace.py#L104?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.excitations" class="docs-object-method">&nbsp;</a> 
```python
@property
excitations(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/AbstractStateSpace.py#L120)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/AbstractStateSpace.py#L120?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.mode" class="docs-object-method">&nbsp;</a> 
```python
@property
mode(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/AbstractStateSpace.py#L127)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/AbstractStateSpace.py#L127?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.get_mode" class="docs-object-method">&nbsp;</a> 
```python
get_mode(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/AbstractStateSpace.py#L130)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/AbstractStateSpace.py#L130?message=Update%20Docs)]
</div>
Returns the mode (indices or excitations) for the held states
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.has_indices" class="docs-object-method">&nbsp;</a> 
```python
@property
has_indices(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/AbstractStateSpace.py#L139)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/AbstractStateSpace.py#L139?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.has_excitations" class="docs-object-method">&nbsp;</a> 
```python
@property
has_excitations(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/AbstractStateSpace.py#L142)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/AbstractStateSpace.py#L142?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.indices" class="docs-object-method">&nbsp;</a> 
```python
@property
indices(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/AbstractStateSpace.py#L159)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/AbstractStateSpace.py#L159?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.indexer" class="docs-object-method">&nbsp;</a> 
```python
@property
indexer(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/AbstractStateSpace.py#L166)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/AbstractStateSpace.py#L166?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.exc_indexer" class="docs-object-method">&nbsp;</a> 
```python
@property
exc_indexer(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/AbstractStateSpace.py#L179)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/AbstractStateSpace.py#L179?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.find" class="docs-object-method">&nbsp;</a> 
```python
find(self, to_search, check=True, minimal_dtype=False, dtype=None, missing_val='raise'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/AbstractStateSpace.py#L188)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/AbstractStateSpace.py#L188?message=Update%20Docs)]
</div>
Finds the indices of a set of indices inside the space
  - `to_search`: `np.ndarray | AbstractStateSpace`
    > array of ints
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.__len__" class="docs-object-method">&nbsp;</a> 
```python
__len__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/AbstractStateSpace.py#L222)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/AbstractStateSpace.py#L222?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.unique_len" class="docs-object-method">&nbsp;</a> 
```python
@property
unique_len(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/AbstractStateSpace.py#L228)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/AbstractStateSpace.py#L228?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.unique_indices" class="docs-object-method">&nbsp;</a> 
```python
@property
unique_indices(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/AbstractStateSpace.py#L235)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/AbstractStateSpace.py#L235?message=Update%20Docs)]
</div>
Returns the unique indices
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.unique_excitations" class="docs-object-method">&nbsp;</a> 
```python
@property
unique_excitations(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/AbstractStateSpace.py#L244)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/AbstractStateSpace.py#L244?message=Update%20Docs)]
</div>
Returns the unique excitations
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.as_indices" class="docs-object-method">&nbsp;</a> 
```python
as_indices(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/AbstractStateSpace.py#L253)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/AbstractStateSpace.py#L253?message=Update%20Docs)]
</div>
Returns the index version of the stored states
  - `:returns`: `np.ndarray`
    >


<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.as_unique_indices" class="docs-object-method">&nbsp;</a> 
```python
as_unique_indices(self, sort=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/AbstractStateSpace.py#L262)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/AbstractStateSpace.py#L262?message=Update%20Docs)]
</div>
Returns unique indices
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.as_excitations" class="docs-object-method">&nbsp;</a> 
```python
as_excitations(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/AbstractStateSpace.py#L286)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/AbstractStateSpace.py#L286?message=Update%20Docs)]
</div>
Returns the excitation version of the stored states
  - `:returns`: `np.ndarray`
    >


<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.as_unique_excitations" class="docs-object-method">&nbsp;</a> 
```python
as_unique_excitations(self, sort=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/AbstractStateSpace.py#L295)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/AbstractStateSpace.py#L295?message=Update%20Docs)]
</div>
Returns unique excitations
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.get_representation_indices" class="docs-object-method">&nbsp;</a> 
```python
get_representation_indices(self, other=None, selection_rules=None, freqs=None, freq_threshold=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/AbstractStateSpace.py#L315)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/AbstractStateSpace.py#L315?message=Update%20Docs)]
</div>
Returns bra and ket indices that can be used as indices to generate representations
  - `other`: `Any`
    > 
  - `selection_rules`: `Any`
    > 
  - `freqs`: `Any`
    > 
  - `freq_threshold`: `Any`
    > 
  - `:returns`: `(np.ndarray, np.ndarray)`
    >


<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.get_representation_brakets" class="docs-object-method">&nbsp;</a> 
```python
get_representation_brakets(self, other=None, selection_rules=None, freqs=None, freq_threshold=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/AbstractStateSpace.py#L339)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/AbstractStateSpace.py#L339?message=Update%20Docs)]
</div>
Returns a BraKetSpace that can be used as generate representations
  - `other`: `Any`
    > 
  - `selection_rules`: `Any`
    > 
  - `freqs`: `Any`
    > 
  - `freq_threshold`: `Any`
    > 
  - `:returns`: `BraKetSpace`
    >


<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.take_states" class="docs-object-method">&nbsp;</a> 
```python
take_states(self, states): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/AbstractStateSpace.py#L362)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/AbstractStateSpace.py#L362?message=Update%20Docs)]
</div>
Takes the intersection of self and the specified states
  - `states`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.take_subspace" class="docs-object-method">&nbsp;</a> 
```python
take_subspace(self, sel): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/AbstractStateSpace.py#L372)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/AbstractStateSpace.py#L372?message=Update%20Docs)]
</div>
Takes a subset of the states
  - `sel`: `Any`
    > 
  - `:returns`: `AbstractStateSpace`
    >


<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.take_subdimensions" class="docs-object-method">&nbsp;</a> 
```python
take_subdimensions(self, inds): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/AbstractStateSpace.py#L383)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/AbstractStateSpace.py#L383?message=Update%20Docs)]
</div>
Takes a subset of the state dimensions
  - `sel`: `Any`
    > 
  - `:returns`: `AbstractStateSpace`
    >


<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.drop_states" class="docs-object-method">&nbsp;</a> 
```python
drop_states(self, states): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/AbstractStateSpace.py#L395)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/AbstractStateSpace.py#L395?message=Update%20Docs)]
</div>
Takes the difference of self and the specified states
  - `states`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.drop_subspace" class="docs-object-method">&nbsp;</a> 
```python
drop_subspace(self, sel): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/AbstractStateSpace.py#L405)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/AbstractStateSpace.py#L405?message=Update%20Docs)]
</div>
Drops a subset of the states
  - `sel`: `Any`
    > 
  - `:returns`: `AbstractStateSpace`
    >


<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.drop_subdimensions" class="docs-object-method">&nbsp;</a> 
```python
drop_subdimensions(self, inds): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/AbstractStateSpace.py#L416)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/AbstractStateSpace.py#L416?message=Update%20Docs)]
</div>
Drops a subset of the state dimensions
  - `sel`: `Any`
    > 
  - `:returns`: `AbstractStateSpace`
    >


<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.get_states_with_quanta" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_states_with_quanta(cls, n, ndim): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L428)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L428?message=Update%20Docs)]
</div>
Returns the states with number of quanta equal to n
  - `quanta`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.num_states_with_quanta" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
num_states_with_quanta(cls, n, ndim): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L450)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L450?message=Update%20Docs)]
</div>
Returns the states with number of quanta equal to n
  - `quanta`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.to_single" class="docs-object-method">&nbsp;</a> 
```python
to_single(self, track_excitations=True, track_indices=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/AbstractStateSpace.py#L464)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/AbstractStateSpace.py#L464?message=Update%20Docs)]
</div>
Flattens any complicated state space structure into a
single space like a `BasisStateSpace`
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.split" class="docs-object-method">&nbsp;</a> 
```python
split(self, chunksize): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/AbstractStateSpace.py#L478)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/AbstractStateSpace.py#L478?message=Update%20Docs)]
</div>
Subclass overridable function to allow for spaces to be
split up into chunks
  - `chunksize`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.share" class="docs-object-method">&nbsp;</a> 
```python
share(self, shared_memory_manager): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/AbstractStateSpace.py#L489)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/AbstractStateSpace.py#L489?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.AbstractStateSpace.unshare" class="docs-object-method">&nbsp;</a> 
```python
unshare(self, shared_memory_manager): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/AbstractStateSpace.py#L491)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/AbstractStateSpace.py#L491?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/StateSpaces/AbstractStateSpace.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/StateSpaces/AbstractStateSpace.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/StateSpaces/AbstractStateSpace.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/StateSpaces/AbstractStateSpace.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces.py#L26?message=Update%20Docs)   
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