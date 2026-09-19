## <a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG">PTTensorCoeffProductDAG</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic.py#L3861)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic.py#L3861?message=Update%20Docs)]
</div>

Lazy, canonical operation DAG for tensor-coefficient expressions.

Leaves retain the existing dictionary-backed ``PTTensorCoeffProductSum``.
Algebra builds immutable nodes and materializes the legacy representation
only at compatibility boundaries such as serialization or the current
tensor evaluator.  Materialization is cached on every node, so repeated
evaluation never replays the derivation tree.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, node, ndim=None, reduced=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic.py#L3876)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic.py#L3876?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.from_sum" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_sum(cls, expression): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L3906)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L3906?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.clear_caches" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
clear_caches(cls): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L3927)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L3927?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.cache_info" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
cache_info(cls): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L3934)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L3934?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.terms" class="docs-object-method">&nbsp;</a> 
```python
@property
terms(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L3943)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L3943?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.ndim" class="docs-object-method">&nbsp;</a> 
```python
@property
ndim(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L3948)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L3948?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.operator_keys" class="docs-object-method">&nbsp;</a> 
```python
@property
operator_keys(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L3954)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L3954?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.__hash__" class="docs-object-method">&nbsp;</a> 
```python
__hash__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L3972)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L3972?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.__eq__" class="docs-object-method">&nbsp;</a> 
```python
__eq__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L3980)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L3980?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L3989)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L3989?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.format_expr" class="docs-object-method">&nbsp;</a> 
```python
format_expr(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L3992)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L3992?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.to_eager" class="docs-object-method">&nbsp;</a> 
```python
to_eager(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4000)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4000?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.prep_serialization_dict" class="docs-object-method">&nbsp;</a> 
```python
prep_serialization_dict(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4134)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4134?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.mutate" class="docs-object-method">&nbsp;</a> 
```python
mutate(self, terms=<DefaultValues.DEFAULT: 'default'>, *, prefactor=<DefaultValues.DEFAULT: 'default'>, ndim=<DefaultValues.DEFAULT: 'default'>, inds_map=<DefaultValues.DEFAULT: 'default'>, canonicalize=<DefaultValues.DEFAULT: 'default'>, reduced=<DefaultValues.DEFAULT: 'default'>): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4137)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4137?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.audit" class="docs-object-method">&nbsp;</a> 
```python
audit(self, target=None, required_dimension=None, ignore_constants=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4154)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4154?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.get_inds" class="docs-object-method">&nbsp;</a> 
```python
get_inds(self, key): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4159)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4159?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.ensure_dimension" class="docs-object-method">&nbsp;</a> 
```python
ensure_dimension(self, ndim): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4164)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4164?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.sort" class="docs-object-method">&nbsp;</a> 
```python
sort(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4171)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4171?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.permute" class="docs-object-method">&nbsp;</a> 
```python
permute(self, new_inds, check_perm=True, allow_padding=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4174)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4174?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.free_up_indices" class="docs-object-method">&nbsp;</a> 
```python
free_up_indices(self, start, stop): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4181)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4181?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.shift" class="docs-object-method">&nbsp;</a> 
```python
shift(self, shift): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4186)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4186?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.shift_energies" class="docs-object-method">&nbsp;</a> 
```python
shift_energies(self, change): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4192)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4192?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.scale" class="docs-object-method">&nbsp;</a> 
```python
scale(self, scaling): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4198)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4198?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.flip_energy_terms" class="docs-object-method">&nbsp;</a> 
```python
flip_energy_terms(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4206)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4206?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.filter_coefficients" class="docs-object-method">&nbsp;</a> 
```python
filter_coefficients(self, terms, mode='match'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4224)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4224?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.filter_energies" class="docs-object-method">&nbsp;</a> 
```python
filter_energies(self, terms, mode='match'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4229)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4229?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.prune_operators" class="docs-object-method">&nbsp;</a> 
```python
prune_operators(self, ops): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4234)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4234?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.combine" class="docs-object-method">&nbsp;</a> 
```python
combine(self, combine_coeffs=False, combine_subterms=True, combine_energies=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4237)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4237?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.mul_along" class="docs-object-method">&nbsp;</a> 
```python
mul_along(self, other, inds, remainder=None, index_classes=None, mapping=None, baseline=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4247)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4247?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.rmul_along" class="docs-object-method">&nbsp;</a> 
```python
rmul_along(self, other, inds, remainder=None, mapping=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4263)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4263?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.mul_simple" class="docs-object-method">&nbsp;</a> 
```python
mul_simple(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4272)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4272?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.rmul_simple" class="docs-object-method">&nbsp;</a> 
```python
rmul_simple(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4280)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4280?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.__add__" class="docs-object-method">&nbsp;</a> 
```python
__add__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4285)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4285?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.__radd__" class="docs-object-method">&nbsp;</a> 
```python
__radd__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4295)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4295?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.__mul__" class="docs-object-method">&nbsp;</a> 
```python
__mul__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4298)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4298?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.__rmul__" class="docs-object-method">&nbsp;</a> 
```python
__rmul__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4301)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4301?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic.py#L3861?message=Update%20Docs)   
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