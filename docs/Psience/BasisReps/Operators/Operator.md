## <a id="Psience.BasisReps.Operators.Operator">Operator</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Operators.py#L24)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Operators.py#L24?message=Update%20Docs)]
</div>

Provides a (usually) _lazy_ representation of an operator, which allows things like
QQQ and pQp to be calculated block-by-block.
Crucially, the underlying basis for the operator is assumed to be orthonormal.

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

<a id="Psience.BasisReps.Operators.Operator.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, funcs, quanta, prod_dim=None, symmetries=None, selection_rules=None, selection_rule_steps=None, parallelizer=None, logger=None, zero_threshold=None, skipped_indices=None, chunk_size=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Operators.py#L30)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Operators.py#L30?message=Update%20Docs)]
</div>


- `funcs`: `callable | Iterable[callable]`
    >The functions use to calculate representation
- `quanta`: `int | Iterable[int]`
    >The number of quanta to do the deepest-level calculations up to (also tells us dimension)
- `prod_dim`: `int | None`
    >The number of functions in `funcs`, if `funcs` is a direct term generator
- `symmetries`: `Iterable[int] | None`
    >Labels for the funcs where if two funcs share a label they are symmetry equivalent

<a id="Psience.BasisReps.Operators.Operator.clear_cache" class="docs-object-method">&nbsp;</a> 
```python
clear_cache(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Operators.py#L81)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Operators.py#L81?message=Update%20Docs)]
</div>


- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Operators.Operator.ndim" class="docs-object-method">&nbsp;</a> 
```python
@property
ndim(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Operators.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Operators.py#L?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Operators.Operator.shape" class="docs-object-method">&nbsp;</a> 
```python
@property
shape(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Operators.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Operators.py#L?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Operators.Operator.load_parallelizer" class="docs-object-method">&nbsp;</a> 
```python
load_parallelizer(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Operators.py#L101)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Operators.py#L101?message=Update%20Docs)]
</div>

Loads the held parallelizer if needed
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Operators.Operator.parallelizer" class="docs-object-method">&nbsp;</a> 
```python
@property
parallelizer(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Operators.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Operators.py#L?message=Update%20Docs)]
</div>

Loads a parallelizer that can be used to speed up various bits of the calculation
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Operators.Operator.get_inner_indices" class="docs-object-method">&nbsp;</a> 
```python
get_inner_indices(self, reduced_inds=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Operators.py#L128)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Operators.py#L128?message=Update%20Docs)]
</div>

Gets the n-dimensional array of ijkl (e.g.) indices that functions will map over
        Basically returns the indices of the inner-most tensor
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Operators.Operator.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Operators.py#L146)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Operators.py#L146?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Operators.Operator.__getstate__" class="docs-object-method">&nbsp;</a> 
```python
__getstate__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Operators.py#L166)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Operators.py#L166?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Operators.Operator.filter_symmetric_indices" class="docs-object-method">&nbsp;</a> 
```python
filter_symmetric_indices(self, inds): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Operators.py#L198)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Operators.py#L198?message=Update%20Docs)]
</div>

Determines which inds are symmetry unique.
        For something like `qqq` all permutations are equivalent, but for `pqp` we have `pi qj pj` distinct from `pj qj pi`.
        This means for `qqq` we have `(1, 0, 0) == (0, 1, 0)` but for `pqp` we only have stuff like `(2, 0, 1) == (1, 0, 2)` .
- `inds`: `np.ndarray`
    >indices to filter symmetric bits out of
- `:returns`: `_`
    >symmetric indices & inverse map

<a id="Psience.BasisReps.Operators.Operator.get_elements" class="docs-object-method">&nbsp;</a> 
```python
get_elements(self, idx, parallelizer=None, check_orthogonality=True, memory_constrained=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Operators.py#L756)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Operators.py#L756?message=Update%20Docs)]
</div>

Calculates a subset of elements
- `idx`: `BraKetSpace`
    >bra and ket states as tuples of elements
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Operators.Operator.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Operators.py#L797)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Operators.py#L797?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Operators.Operator.get_transformed_space" class="docs-object-method">&nbsp;</a> 
```python
get_transformed_space(self, base_space, rules=None, parallelizer=None, logger=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Operators.py#L804)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Operators.py#L804?message=Update%20Docs)]
</div>

Returns the space one would get from applying
        the selection rules from this operator
- `base_space`: `BasisStateSpace`
    >No description...
- `:returns`: `SelectionRuleStateSpace`
    >No description...

<a id="Psience.BasisReps.Operators.Operator.apply_reduced" class="docs-object-method">&nbsp;</a> 
```python
apply_reduced(self, base_space, parallelizer=None, logger=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Operators.py#L1011)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Operators.py#L1011?message=Update%20Docs)]
</div>

Takes a base space as input and applies the held selection rules in semi-efficient
        fashion only on the indices that can change and then uses this to compute all matrix
        elements, returning then the final generated space
- `base_space`: `BasisStateSpace | PermutationallyReducedStateSpace`
    >No description...
- `:returns`: `tuple[SparseArray, BraKetSpace]`
    >No description...

 </div>
</div>




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/Operators/Operator.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/Operators/Operator.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/Operators/Operator.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/Operators/Operator.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Operators.py#L24?message=Update%20Docs)