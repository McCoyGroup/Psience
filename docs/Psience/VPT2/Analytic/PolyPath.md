## <a id="Psience.VPT2.Analytic.PolyPath">PolyPath</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic.py#L2491)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic.py#L2491?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic.py#L2513)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic.py#L2513?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.is_zero" class="docs-object-method">&nbsp;</a> 
```python
@property
is_zero(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L2554)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L2554?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.from_coeffs" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_coeffs(cls, coeffs, prefactor=1, idx=None, steps=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L2558)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L2558?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.from_polynomial" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_polynomial(cls, poly): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L2563)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L2563?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.clear_caches" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
clear_caches(cls): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L2581)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L2581?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.cache_info" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
cache_info(cls): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L2592)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L2592?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.path_terms" class="docs-object-method">&nbsp;</a> 
```python
@property
path_terms(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L2605)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L2605?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.polys" class="docs-object-method">&nbsp;</a> 
```python
@property
polys(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L2611)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L2611?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.to_eager" class="docs-object-method">&nbsp;</a> 
```python
to_eager(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L2630)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L2630?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.prep_serialization_dict" class="docs-object-method">&nbsp;</a> 
```python
prep_serialization_dict(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L2712)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L2712?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.mutate" class="docs-object-method">&nbsp;</a> 
```python
mutate(self, polynomials=<DefaultValues.DEFAULT: 'default'>, prefactor=<DefaultValues.DEFAULT: 'default'>, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L2720)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L2720?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.ndim" class="docs-object-method">&nbsp;</a> 
```python
@property
ndim(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L2730)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L2730?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.order" class="docs-object-method">&nbsp;</a> 
```python
@property
order(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L2736)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L2736?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.__hash__" class="docs-object-method">&nbsp;</a> 
```python
__hash__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L2745)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L2745?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.__eq__" class="docs-object-method">&nbsp;</a> 
```python
__eq__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L2750)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L2750?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L2757)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L2757?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.format_expr" class="docs-object-method">&nbsp;</a> 
```python
format_expr(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L2763)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L2763?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.audit" class="docs-object-method">&nbsp;</a> 
```python
audit(self, target=None, ignore_constants=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L2767)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L2767?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.ensure_dimension" class="docs-object-method">&nbsp;</a> 
```python
ensure_dimension(self, ndim): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L2779)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L2779?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.pad" class="docs-object-method">&nbsp;</a> 
```python
pad(self, left_right_pads): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L2789)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L2789?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.permute" class="docs-object-method">&nbsp;</a> 
```python
permute(self, new_inds, check_perm=True, allow_padding=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L2801)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L2801?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.permutation_sum" class="docs-object-method">&nbsp;</a> 
```python
permutation_sum(self, permutations, check_perm=True, allow_padding=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L2820)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L2820?message=Update%20Docs)]
</div>
Represent a symmetry sum without constructing each remapped child.


<a id="Psience.VPT2.Analytic.PolyPath.shift" class="docs-object-method">&nbsp;</a> 
```python
shift(self, shift): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L2856)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L2856?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.scale" class="docs-object-method">&nbsp;</a> 
```python
scale(self, scaling): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L2863)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L2863?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.evaluate_polynomial" class="docs-object-method">&nbsp;</a> 
```python
evaluate_polynomial(self, substates, node_cache=None, axis_cache=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L2871)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L2871?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.combine" class="docs-object-method">&nbsp;</a> 
```python
combine(self, *args, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L3037)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L3037?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.condense" class="docs-object-method">&nbsp;</a> 
```python
condense(self, inds=None, return_inds=False, check_inds=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L3078)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L3078?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.constant_rescale" class="docs-object-method">&nbsp;</a> 
```python
constant_rescale(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L3100)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L3100?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.mul_simple" class="docs-object-method">&nbsp;</a> 
```python
mul_simple(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L3106)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L3106?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.rmul_simple" class="docs-object-method">&nbsp;</a> 
```python
rmul_simple(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L3129)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L3129?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.mul_along" class="docs-object-method">&nbsp;</a> 
```python
mul_along(self, other, inds, remainder=None, mapping=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L3134)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L3134?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.rmul_along" class="docs-object-method">&nbsp;</a> 
```python
rmul_along(self, other, inds, remainder=None, mapping=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L3177)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L3177?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.__mul__" class="docs-object-method">&nbsp;</a> 
```python
__mul__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L3184)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L3184?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.__rmul__" class="docs-object-method">&nbsp;</a> 
```python
__rmul__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L3187)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L3187?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.__add__" class="docs-object-method">&nbsp;</a> 
```python
__add__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L3190)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L3190?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PolyPath.__radd__" class="docs-object-method">&nbsp;</a> 
```python
__radd__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PolyPath.py#L3208)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PolyPath.py#L3208?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic.py#L2491?message=Update%20Docs)   
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