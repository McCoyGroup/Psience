## <a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG">PTTensorCoeffProductDAG</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic.py#L4762)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic.py#L4762?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic.py#L4777)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic.py#L4777?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.from_sum" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_sum(cls, expression): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L4807)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L4807?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.clear_caches" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
clear_caches(cls): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L4828)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L4828?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.cache_info" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
cache_info(cls): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L4835)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L4835?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.terms" class="docs-object-method">&nbsp;</a> 
```python
@property
terms(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4844)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4844?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.ndim" class="docs-object-method">&nbsp;</a> 
```python
@property
ndim(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4849)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4849?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.operator_keys" class="docs-object-method">&nbsp;</a> 
```python
@property
operator_keys(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4855)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4855?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.__hash__" class="docs-object-method">&nbsp;</a> 
```python
__hash__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4873)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4873?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.__eq__" class="docs-object-method">&nbsp;</a> 
```python
__eq__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4881)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4881?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4890)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4890?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.format_expr" class="docs-object-method">&nbsp;</a> 
```python
format_expr(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4893)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4893?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.to_eager" class="docs-object-method">&nbsp;</a> 
```python
to_eager(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4901)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L4901?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.prep_serialization_dict" class="docs-object-method">&nbsp;</a> 
```python
prep_serialization_dict(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5035)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5035?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.mutate" class="docs-object-method">&nbsp;</a> 
```python
mutate(self, terms=<DefaultValues.DEFAULT: 'default'>, *, prefactor=<DefaultValues.DEFAULT: 'default'>, ndim=<DefaultValues.DEFAULT: 'default'>, inds_map=<DefaultValues.DEFAULT: 'default'>, canonicalize=<DefaultValues.DEFAULT: 'default'>, reduced=<DefaultValues.DEFAULT: 'default'>): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5038)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5038?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.audit" class="docs-object-method">&nbsp;</a> 
```python
audit(self, target=None, required_dimension=None, ignore_constants=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5055)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5055?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.get_inds" class="docs-object-method">&nbsp;</a> 
```python
get_inds(self, key): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5060)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5060?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.ensure_dimension" class="docs-object-method">&nbsp;</a> 
```python
ensure_dimension(self, ndim): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5065)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5065?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.sort" class="docs-object-method">&nbsp;</a> 
```python
sort(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5072)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5072?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.permute" class="docs-object-method">&nbsp;</a> 
```python
permute(self, new_inds, check_perm=True, allow_padding=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5075)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5075?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.free_up_indices" class="docs-object-method">&nbsp;</a> 
```python
free_up_indices(self, start, stop): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5082)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5082?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.shift" class="docs-object-method">&nbsp;</a> 
```python
shift(self, shift): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5087)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5087?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.shift_energies" class="docs-object-method">&nbsp;</a> 
```python
shift_energies(self, change): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5093)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5093?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.scale" class="docs-object-method">&nbsp;</a> 
```python
scale(self, scaling): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5099)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5099?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.flip_energy_terms" class="docs-object-method">&nbsp;</a> 
```python
flip_energy_terms(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5107)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5107?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.filter_coefficients" class="docs-object-method">&nbsp;</a> 
```python
filter_coefficients(self, terms, mode='match'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5125)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5125?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.filter_energies" class="docs-object-method">&nbsp;</a> 
```python
filter_energies(self, terms, mode='match'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5130)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5130?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.prune_operators" class="docs-object-method">&nbsp;</a> 
```python
prune_operators(self, ops): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5135)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5135?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.combine" class="docs-object-method">&nbsp;</a> 
```python
combine(self, combine_coeffs=False, combine_subterms=True, combine_energies=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5138)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5138?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.mul_along" class="docs-object-method">&nbsp;</a> 
```python
mul_along(self, other, inds, remainder=None, index_classes=None, mapping=None, baseline=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5148)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5148?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.rmul_along" class="docs-object-method">&nbsp;</a> 
```python
rmul_along(self, other, inds, remainder=None, mapping=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5164)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5164?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.mul_simple" class="docs-object-method">&nbsp;</a> 
```python
mul_simple(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5173)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5173?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.rmul_simple" class="docs-object-method">&nbsp;</a> 
```python
rmul_simple(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5181)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5181?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.__add__" class="docs-object-method">&nbsp;</a> 
```python
__add__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5186)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5186?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.__radd__" class="docs-object-method">&nbsp;</a> 
```python
__radd__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5196)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5196?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.__mul__" class="docs-object-method">&nbsp;</a> 
```python
__mul__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5199)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5199?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PTTensorCoeffProductDAG.__rmul__" class="docs-object-method">&nbsp;</a> 
```python
__rmul__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5202)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PTTensorCoeffProductDAG.py#L5202?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic.py#L4762?message=Update%20Docs)   
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