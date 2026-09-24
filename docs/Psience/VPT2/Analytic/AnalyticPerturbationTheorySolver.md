## <a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver">AnalyticPerturbationTheorySolver</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic.py#L940)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic.py#L940?message=Update%20Docs)]
</div>

A re-attempt at using the recursive expressions
to provide simpler code for getting APT expressions







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
operator_expansion_index: int
```
<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, hamiltonian_expansion, logger=None, checkpoint=None, allowed_terms=None, allowed_coefficients=None, disallowed_coefficients=None, allowed_energy_changes=None, intermediate_normalization=None, polynomial_representation='path'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic.py#L945)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic.py#L945?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.from_order" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_order(cls, order, internals=True, logger=None, checkpoint=None, allowed_terms=None, allowed_coefficients=None, disallowed_coefficients=None, allowed_energy_changes=None, intermediate_normalization=None, polynomial_representation='path'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L966)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L966?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.modify_hamiltonian" class="docs-object-method">&nbsp;</a> 
```python
modify_hamiltonian(self, hamiltonian_corrections): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L1034)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L1034?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.get_correction" class="docs-object-method">&nbsp;</a> 
```python
get_correction(self, key, cls, order, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L1054)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L1054?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.shifted_hamiltonian_correction" class="docs-object-method">&nbsp;</a> 
```python
shifted_hamiltonian_correction(self, order, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L1074)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L1074?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.energy_correction" class="docs-object-method">&nbsp;</a> 
```python
energy_correction(self, order, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L1077)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L1077?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.wavefunction_correction" class="docs-object-method">&nbsp;</a> 
```python
wavefunction_correction(self, order, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L1080)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L1080?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.overlap_correction" class="docs-object-method">&nbsp;</a> 
```python
overlap_correction(self, order, degenerate_changes=None, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L1083)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L1083?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.full_wavefunction_correction" class="docs-object-method">&nbsp;</a> 
```python
full_wavefunction_correction(self, order, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L1090)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L1090?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.operator_correction" class="docs-object-method">&nbsp;</a> 
```python
operator_correction(self, order, operator_type=None, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L1093)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L1093?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.operator_degenerate_correction" class="docs-object-method">&nbsp;</a> 
```python
operator_degenerate_correction(self, order, /, degenerate_changes, operator_type=None, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L1096)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L1096?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.reexpressed_hamiltonian" class="docs-object-method">&nbsp;</a> 
```python
reexpressed_hamiltonian(self, order, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L1101)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L1101?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.reexpressed_hamiltonian_degenerate_correction" class="docs-object-method">&nbsp;</a> 
```python
reexpressed_hamiltonian_degenerate_correction(self, order, /, degenerate_changes, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L1103)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L1103?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.operator_expansion_terms" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
operator_expansion_terms(cls, order, logger=None, base_index=None, operator_type=None, polynomial_representation='path'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1108)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1108?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.clear_caches" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
clear_caches(cls): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1162)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1162?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.polynomial_cache_info" class="docs-object-method">&nbsp;</a> 
```python
polynomial_cache_info(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L1184)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L1184?message=Update%20Docs)]
</div>
Return backend-specific counters useful for path/eager timing comparisons.
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic.py#L940?message=Update%20Docs)   
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