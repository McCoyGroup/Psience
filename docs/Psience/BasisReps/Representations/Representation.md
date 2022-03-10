## <a id="Psience.BasisReps.Representations.Representation">Representation</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L24)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L24?message=Update%20Docs)]
</div>

A `Representation` provides a simple interface to build matrix representations of operators expressed
in high-dimensional spaces.

<a id="Psience.BasisReps.Representations.Representation.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, compute, basis, name=None, logger=None, selection_rules=None, selection_rule_steps=None, memory_constrained=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L31)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L31?message=Update%20Docs)]
</div>


- `compute`: `callable | Operator`
    >the function that turns indices into values
- `basis`: `RepresentationBasis`
    >the basis quanta used in the representations
- `logger`: `None | Logger`
    >logger for printing out debug info

<a id="Psience.BasisReps.Representations.Representation.parallelizer" class="docs-object-method">&nbsp;</a> 
```python
@property
parallelizer(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.Representation.compute" class="docs-object-method">&nbsp;</a> 
```python
compute(self, inds, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L70)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L70?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.Representation.compute_cached" class="docs-object-method">&nbsp;</a> 
```python
compute_cached(self, inds): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L76)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L76?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.Representation.chunk_size" class="docs-object-method">&nbsp;</a> 
```python
@property
chunk_size(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.Representation.clear_cache" class="docs-object-method">&nbsp;</a> 
```python
clear_cache(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L89)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L89?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.Representation.diag" class="docs-object-method">&nbsp;</a> 
```python
@property
diag(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.Representation.ndims" class="docs-object-method">&nbsp;</a> 
```python
@property
ndims(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.Representation.dim_inds" class="docs-object-method">&nbsp;</a> 
```python
@property
dim_inds(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.Representation.get_brakets" class="docs-object-method">&nbsp;</a> 
```python
get_brakets(self, states, check_orthogonality=True, memory_constrained=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L138)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L138?message=Update%20Docs)]
</div>

Computes term elements based on getting a BraKetSpace.
        Can directly pass element specs through, since the shape management shit
        is managed by the BraKetSpace
- `states`: `BraKetSpace | Tuple[np.ndarray, np.ndarray]`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Representations.Representation.get_element" class="docs-object-method">&nbsp;</a> 
```python
get_element(self, n, m): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L159)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L159?message=Update%20Docs)]
</div>

Computes term elements.
        Determines first whether it needs to pull single elements or blocks of them.
- `n`: `Any`
    >No description...
- `m`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Representations.Representation.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L249)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L249?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.Representation.__rmul__" class="docs-object-method">&nbsp;</a> 
```python
__rmul__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L266)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L266?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.Representation.__mul__" class="docs-object-method">&nbsp;</a> 
```python
__mul__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L278)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L278?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.Representation.__add__" class="docs-object-method">&nbsp;</a> 
```python
__add__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L291)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L291?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.Representation.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L324)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L324?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.Representation.selection_rules" class="docs-object-method">&nbsp;</a> 
```python
@property
selection_rules(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L?message=Update%20Docs)]
</div>


- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Representations.Representation.selection_rule_steps" class="docs-object-method">&nbsp;</a> 
```python
@property
selection_rule_steps(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L?message=Update%20Docs)]
</div>


- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Representations.Representation.get_transformed_space" class="docs-object-method">&nbsp;</a> 
```python
get_transformed_space(self, space, parallelizer=None, logger=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L372)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L372?message=Update%20Docs)]
</div>

Returns the state space obtained by using the
        held operator to transform `space`
- `space`: `Any`
    >No description...
- `:returns`: `SelectionRuleStateSpace`
    >connected state spaces

<a id="Psience.BasisReps.Representations.Representation.apply" class="docs-object-method">&nbsp;</a> 
```python
apply(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L398)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L398?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Representations.Representation.get_representation_matrix" class="docs-object-method">&nbsp;</a> 
```python
get_representation_matrix(self, coupled_space, total_space, filter_space=None, diagonal=False, logger=None, zero_element_warning=True, clear_sparse_caches=True, clear_operator_caches=True, assume_symmetric=True, remove_duplicates=True, memory_constrained=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L559)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L559?message=Update%20Docs)]
</div>

Actively constructs a perturbation theory Hamiltonian representation
- `h`: `Any`
    >No description...
- `cs`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Representations.Representation.get_diagonal_representation" class="docs-object-method">&nbsp;</a> 
```python
get_diagonal_representation(self, coupled_space, total_space, logger=None, zero_element_warning=True, clear_sparse_caches=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Representations.py#L796)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L796?message=Update%20Docs)]
</div>

Actively constructs a perturbation theory Hamiltonian representation
- `h`: `Any`
    >No description...
- `cs`: `Any`
    >No description...
- `:returns`: `_`
    >No description...



___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/ci/docs/Psience/BasisReps/Representations/Representation.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/ci/docs/Psience/BasisReps/Representations/Representation.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/ci/docs/Psience/BasisReps/Representations/Representation.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/ci/docs/Psience/BasisReps/Representations/Representation.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Representations.py#L24?message=Update%20Docs)