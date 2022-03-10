## <a id="Psience.BasisReps.Bases.SimpleProductBasis">SimpleProductBasis</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Bases.py#L306)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Bases.py#L306?message=Update%20Docs)]
</div>

Defines a direct product basis from a 1D basis.
Mixed product bases aren't currently supported, but this provides
at least a sample for how that kind of things could be
generated.

```python
array_indexer_cutoff: int
```
<a id="Psience.BasisReps.Bases.SimpleProductBasis.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, basis_type, n_quanta, indexer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Bases.py#L316)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Bases.py#L316?message=Update%20Docs)]
</div>


- `basis_type`: `type`
    >the type of basis to do a product over
- `n_quanta`: `Iterable[int]`
    >the number of quanta for the representations
- `indexer`: `BaseStateIndexer`
    >an object that can turn state specs into indices and indices into state specs

<a id="Psience.BasisReps.Bases.SimpleProductBasis.to_state" class="docs-object-method">&nbsp;</a> 
```python
to_state(self, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Bases.py#L344)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Bases.py#L344?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Bases.SimpleProductBasis.from_state" class="docs-object-method">&nbsp;</a> 
```python
from_state(data, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Bases.py#L350)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Bases.py#L350?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Bases.SimpleProductBasis.ndim" class="docs-object-method">&nbsp;</a> 
```python
@property
ndim(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Bases.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Bases.py#L?message=Update%20Docs)]
</div>

Provides the number of dimensions of the basis
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Bases.SimpleProductBasis.dimensions" class="docs-object-method">&nbsp;</a> 
```python
@property
dimensions(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Bases.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Bases.py#L?message=Update%20Docs)]
</div>

Provides the dimensions of the basis
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Bases.SimpleProductBasis.__eq__" class="docs-object-method">&nbsp;</a> 
```python
__eq__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Bases.py#L376)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Bases.py#L376?message=Update%20Docs)]
</div>


- `other`: `SimpleProductBasis`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Bases.SimpleProductBasis.quanta" class="docs-object-method">&nbsp;</a> 
```python
@property
quanta(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Bases.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Bases.py#L?message=Update%20Docs)]
</div>

Provides the quanta in each dimension of the basis
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Bases.SimpleProductBasis.selection_rules_mapping" class="docs-object-method">&nbsp;</a> 
```python
@property
selection_rules_mapping(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Bases.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Bases.py#L?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Bases.SimpleProductBasis.ravel_state_inds" class="docs-object-method">&nbsp;</a> 
```python
ravel_state_inds(self, idx): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Bases.py#L406)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Bases.py#L406?message=Update%20Docs)]
</div>

Converts state indices from an array of quanta to an array of indices
- `idx`: `Iterable[Iterable[int]]`
    >indices
- `:returns`: `tuple[int]`
    >array of state indices in the basis

<a id="Psience.BasisReps.Bases.SimpleProductBasis.unravel_state_inds" class="docs-object-method">&nbsp;</a> 
```python
unravel_state_inds(self, idx): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Bases.py#L418)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Bases.py#L418?message=Update%20Docs)]
</div>

Converts state indices from an array of ints to an array of quanta
- `idx`: `Iterable[int]`
    >indices
- `:returns`: `tuple[tuple[int]]`
    >array of state tuples in the basis

<a id="Psience.BasisReps.Bases.SimpleProductBasis.get_function" class="docs-object-method">&nbsp;</a> 
```python
get_function(self, idx): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Bases.py#L432)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Bases.py#L432?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Bases.SimpleProductBasis.operator" class="docs-object-method">&nbsp;</a> 
```python
operator(self, *terms, coeffs=None, axes=None, parallelizer=None, logger=None, chunk_size=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Bases.py#L436)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Bases.py#L436?message=Update%20Docs)]
</div>

Builds an operator based on supplied terms, remapping names where possible.
        If `coeffs` or `axes` are supplied, a `ContractedOperator` is built.
- `terms`: `Any`
    >No description...
- `coeffs`: `Any`
    >No description...
- `axes`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Bases.SimpleProductBasis.representation" class="docs-object-method">&nbsp;</a> 
```python
representation(self, *terms, coeffs=None, name=None, axes=None, logger=None, parallelizer=None, chunk_size=None, memory_constrained=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Bases.py#L475)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Bases.py#L475?message=Update%20Docs)]
</div>

Provides a representation of a product operator specified by _terms_.
        If `coeffs` or `axes` are supplied, a `ContractedOperator` is built.
- `terms`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Bases.SimpleProductBasis.x" class="docs-object-method">&nbsp;</a> 
```python
x(self, n): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Bases.py#L492)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Bases.py#L492?message=Update%20Docs)]
</div>

Returns the representation of x in the multi-dimensional basis with every term evaluated up to n quanta
        Whether this is what we want or not is still TBD
- `n`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Bases.SimpleProductBasis.p" class="docs-object-method">&nbsp;</a> 
```python
p(self, n): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Bases.py#L502)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Bases.py#L502?message=Update%20Docs)]
</div>

Returns the representation of p in the multi-dimensional basis with every term evaluated up to n quanta
        Whether this is what we want or not is still TBD
- `n`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Bases.SimpleProductBasis.take_subdimensions" class="docs-object-method">&nbsp;</a> 
```python
take_subdimensions(self, dims): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Bases.py#L513)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Bases.py#L513?message=Update%20Docs)]
</div>

Casts down to lower dimensional space
- `dims`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Bases.SimpleProductBasis.get_state_space" class="docs-object-method">&nbsp;</a> 
```python
get_state_space(self, quanta): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Bases.py#L525)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Bases.py#L525?message=Update%20Docs)]
</div>



___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/ci/docs/Psience/BasisReps/Bases/SimpleProductBasis.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/ci/docs/Psience/BasisReps/Bases/SimpleProductBasis.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/ci/docs/Psience/BasisReps/Bases/SimpleProductBasis.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/ci/docs/Psience/BasisReps/Bases/SimpleProductBasis.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Bases.py#L306?message=Update%20Docs)