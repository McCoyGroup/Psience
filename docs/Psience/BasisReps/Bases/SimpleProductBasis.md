## <a id="Psience.BasisReps.Bases.SimpleProductBasis">SimpleProductBasis</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Bases.py#L310)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases.py#L310?message=Update%20Docs)]
</div>

Defines a direct product basis from a 1D basis.
Mixed product bases aren't currently supported, but this provides
at least a sample for how that kind of things could be
generated.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
array_indexer_cutoff: int
```
<a id="Psience.BasisReps.Bases.SimpleProductBasis.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, basis_type, n_quanta, indexer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Bases/SimpleProductBasis.py#L320)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases/SimpleProductBasis.py#L320?message=Update%20Docs)]
</div>

  - `basis_type`: `type`
    > the type of basis to do a product over
  - `n_quanta`: `Iterable[int]`
    > the number of quanta for the representations
  - `indexer`: `BaseStateIndexer`
    > an object that can turn state specs into indices and indices into state specs


<a id="Psience.BasisReps.Bases.SimpleProductBasis.to_state" class="docs-object-method">&nbsp;</a> 
```python
to_state(self, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Bases/SimpleProductBasis.py#L348)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases/SimpleProductBasis.py#L348?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Bases.SimpleProductBasis.from_state" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_state(cls, data, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L354)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L354?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Bases.SimpleProductBasis.ndim" class="docs-object-method">&nbsp;</a> 
```python
@property
ndim(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Bases/SimpleProductBasis.py#L360)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases/SimpleProductBasis.py#L360?message=Update%20Docs)]
</div>
Provides the number of dimensions of the basis
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.Bases.SimpleProductBasis.dimensions" class="docs-object-method">&nbsp;</a> 
```python
@property
dimensions(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Bases/SimpleProductBasis.py#L370)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases/SimpleProductBasis.py#L370?message=Update%20Docs)]
</div>
Provides the dimensions of the basis
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.Bases.SimpleProductBasis.__eq__" class="docs-object-method">&nbsp;</a> 
```python
__eq__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Bases/SimpleProductBasis.py#L380)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases/SimpleProductBasis.py#L380?message=Update%20Docs)]
</div>

  - `other`: `SimpleProductBasis`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.Bases.SimpleProductBasis.quanta" class="docs-object-method">&nbsp;</a> 
```python
@property
quanta(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Bases/SimpleProductBasis.py#L389)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases/SimpleProductBasis.py#L389?message=Update%20Docs)]
</div>
Provides the quanta in each dimension of the basis
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.Bases.SimpleProductBasis.selection_rules_mapping" class="docs-object-method">&nbsp;</a> 
```python
@property
selection_rules_mapping(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Bases/SimpleProductBasis.py#L406)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases/SimpleProductBasis.py#L406?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Bases.SimpleProductBasis.ravel_state_inds" class="docs-object-method">&nbsp;</a> 
```python
ravel_state_inds(self, idx): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Bases/SimpleProductBasis.py#L410)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases/SimpleProductBasis.py#L410?message=Update%20Docs)]
</div>
Converts state indices from an array of quanta to an array of indices
  - `idx`: `Iterable[Iterable[int]]`
    > indices
  - `:returns`: `tuple[int]`
    > a
r
r
a
y
 
o
f
 
s
t
a
t
e
 
i
n
d
i
c
e
s
 
i
n
 
t
h
e
 
b
a
s
i
s


<a id="Psience.BasisReps.Bases.SimpleProductBasis.unravel_state_inds" class="docs-object-method">&nbsp;</a> 
```python
unravel_state_inds(self, idx): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Bases/SimpleProductBasis.py#L422)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases/SimpleProductBasis.py#L422?message=Update%20Docs)]
</div>
Converts state indices from an array of ints to an array of quanta
  - `idx`: `Iterable[int]`
    > indices
  - `:returns`: `tuple[tuple[int]]`
    > a
r
r
a
y
 
o
f
 
s
t
a
t
e
 
t
u
p
l
e
s
 
i
n
 
t
h
e
 
b
a
s
i
s


<a id="Psience.BasisReps.Bases.SimpleProductBasis.get_function" class="docs-object-method">&nbsp;</a> 
```python
get_function(self, idx): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Bases/SimpleProductBasis.py#L436)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases/SimpleProductBasis.py#L436?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Bases.SimpleProductBasis.operator" class="docs-object-method">&nbsp;</a> 
```python
operator(self, *terms, coeffs=None, axes=None, parallelizer=None, logger=None, chunk_size=None, **operator_settings): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Bases/SimpleProductBasis.py#L440)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases/SimpleProductBasis.py#L440?message=Update%20Docs)]
</div>
Builds an operator based on supplied terms, remapping names where possible.
If `coeffs` or `axes` are supplied, a `ContractedOperator` is built.
  - `terms`: `Any`
    > 
  - `coeffs`: `Any`
    > 
  - `axes`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.Bases.SimpleProductBasis.representation" class="docs-object-method">&nbsp;</a> 
```python
representation(self, *terms, coeffs=None, name=None, axes=None, logger=None, parallelizer=None, chunk_size=None, memory_constrained=False, **operator_settings): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Bases/SimpleProductBasis.py#L484)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases/SimpleProductBasis.py#L484?message=Update%20Docs)]
</div>
Provides a representation of a product operator specified by _terms_.
If `coeffs` or `axes` are supplied, a `ContractedOperator` is built.
  - `terms`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.Bases.SimpleProductBasis.x" class="docs-object-method">&nbsp;</a> 
```python
x(self, n): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Bases/SimpleProductBasis.py#L506)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases/SimpleProductBasis.py#L506?message=Update%20Docs)]
</div>
Returns the representation of x in the multi-dimensional basis with every term evaluated up to n quanta
Whether this is what we want or not is still TBD
  - `n`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.Bases.SimpleProductBasis.p" class="docs-object-method">&nbsp;</a> 
```python
p(self, n): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Bases/SimpleProductBasis.py#L516)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases/SimpleProductBasis.py#L516?message=Update%20Docs)]
</div>
Returns the representation of p in the multi-dimensional basis with every term evaluated up to n quanta
Whether this is what we want or not is still TBD
  - `n`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.Bases.SimpleProductBasis.take_subdimensions" class="docs-object-method">&nbsp;</a> 
```python
take_subdimensions(self, dims): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Bases/SimpleProductBasis.py#L527)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases/SimpleProductBasis.py#L527?message=Update%20Docs)]
</div>
Casts down to lower dimensional space
  - `dims`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.Bases.SimpleProductBasis.get_state_space" class="docs-object-method">&nbsp;</a> 
```python
get_state_space(self, quanta): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Bases/SimpleProductBasis.py#L539)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases/SimpleProductBasis.py#L539?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/Bases/SimpleProductBasis.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/Bases/SimpleProductBasis.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/Bases/SimpleProductBasis.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/Bases/SimpleProductBasis.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases.py#L310?message=Update%20Docs)   
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