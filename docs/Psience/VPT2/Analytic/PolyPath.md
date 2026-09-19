## <a id="Psience.VPT2.Analytic.PolyPath">PolyPath</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic.py#L1595)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic.py#L1595?message=Update%20Docs)]
</div>

Canonical DAG-backed representation of a sum of separable polynomial products.

Each key is an interned :class:`PolyTerm`; the corresponding value is its
scalar coefficient.  Multiplication joins interned axis paths rather than
convolving coefficient arrays.  ``polys`` and ``to_eager`` provide a
compatibility boundary for the legacy implementation.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.VPT2.Analytic.PolyPath.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, terms, reduced=False, node=None, ndim=None, order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic.py#L1617)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic.py#L1617?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.is_zero" class="docs-object-method">&nbsp;</a> 
```python
@property
is_zero(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L1658)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L1658?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.from_coeffs" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_coeffs(cls, coeffs, prefactor=1, idx=None, steps=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1662)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1662?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.from_polynomial" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_polynomial(cls, poly): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1667)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1667?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.clear_caches" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
clear_caches(cls): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1685)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1685?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.cache_info" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
cache_info(cls): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1696)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1696?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.path_terms" class="docs-object-method">&nbsp;</a> 
```python
@property
path_terms(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L1709)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L1709?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.polys" class="docs-object-method">&nbsp;</a> 
```python
@property
polys(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L1715)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L1715?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.to_eager" class="docs-object-method">&nbsp;</a> 
```python
to_eager(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L1734)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L1734?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.prep_serialization_dict" class="docs-object-method">&nbsp;</a> 
```python
prep_serialization_dict(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L1816)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L1816?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.mutate" class="docs-object-method">&nbsp;</a> 
```python
mutate(self, polynomials=<DefaultValues.DEFAULT: 'default'>, prefactor=<DefaultValues.DEFAULT: 'default'>, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L1824)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L1824?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.ndim" class="docs-object-method">&nbsp;</a> 
```python
@property
ndim(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L1834)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L1834?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.order" class="docs-object-method">&nbsp;</a> 
```python
@property
order(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L1840)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L1840?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.__hash__" class="docs-object-method">&nbsp;</a> 
```python
__hash__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L1849)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L1849?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.__eq__" class="docs-object-method">&nbsp;</a> 
```python
__eq__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L1854)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L1854?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L1861)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L1861?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.format_expr" class="docs-object-method">&nbsp;</a> 
```python
format_expr(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L1867)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L1867?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.audit" class="docs-object-method">&nbsp;</a> 
```python
audit(self, target=None, ignore_constants=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L1871)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L1871?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.ensure_dimension" class="docs-object-method">&nbsp;</a> 
```python
ensure_dimension(self, ndim): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L1883)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L1883?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.pad" class="docs-object-method">&nbsp;</a> 
```python
pad(self, left_right_pads): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L1893)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L1893?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.permute" class="docs-object-method">&nbsp;</a> 
```python
permute(self, new_inds, check_perm=True, allow_padding=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L1905)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L1905?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.permutation_sum" class="docs-object-method">&nbsp;</a> 
```python
permutation_sum(self, permutations, check_perm=True, allow_padding=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L1924)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L1924?message=Update%20Docs)]
</div>
Represent a symmetry sum without constructing each remapped child.


<a id="Psience.VPT2.Analytic.PolyPath.shift" class="docs-object-method">&nbsp;</a> 
```python
shift(self, shift): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L1960)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L1960?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.scale" class="docs-object-method">&nbsp;</a> 
```python
scale(self, scaling): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L1967)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L1967?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.evaluate_polynomial" class="docs-object-method">&nbsp;</a> 
```python
evaluate_polynomial(self, substates, node_cache=None, axis_cache=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L1975)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L1975?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.combine" class="docs-object-method">&nbsp;</a> 
```python
combine(self, *args, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L2136)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L2136?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.condense" class="docs-object-method">&nbsp;</a> 
```python
condense(self, inds=None, return_inds=False, check_inds=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L2177)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L2177?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.constant_rescale" class="docs-object-method">&nbsp;</a> 
```python
constant_rescale(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L2199)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L2199?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.mul_simple" class="docs-object-method">&nbsp;</a> 
```python
mul_simple(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L2205)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L2205?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.rmul_simple" class="docs-object-method">&nbsp;</a> 
```python
rmul_simple(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L2228)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L2228?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.mul_along" class="docs-object-method">&nbsp;</a> 
```python
mul_along(self, other, inds, remainder=None, mapping=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L2233)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L2233?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.rmul_along" class="docs-object-method">&nbsp;</a> 
```python
rmul_along(self, other, inds, remainder=None, mapping=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L2276)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L2276?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.__mul__" class="docs-object-method">&nbsp;</a> 
```python
__mul__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L2283)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L2283?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.__rmul__" class="docs-object-method">&nbsp;</a> 
```python
__rmul__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L2286)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L2286?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.__add__" class="docs-object-method">&nbsp;</a> 
```python
__add__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L2289)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L2289?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.__radd__" class="docs-object-method">&nbsp;</a> 
```python
__radd__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L2307)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L2307?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Analytic/PolyPath.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Analytic/PolyPath.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Analytic/PolyPath.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Analytic/PolyPath.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic.py#L1595?message=Update%20Docs)   
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