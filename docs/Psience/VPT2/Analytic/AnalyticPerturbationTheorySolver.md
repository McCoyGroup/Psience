## <a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver">AnalyticPerturbationTheorySolver</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Analytic.py#L41)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Analytic.py#L41?message=Update%20Docs)]
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
__init__(self, hamiltonian_expansion, logger=None, checkpoint=None, allowed_terms=None, allowed_coefficients=None, disallowed_coefficients=None, allowed_energy_changes=None, intermediate_normalization=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L46)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L46?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.from_order" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_order(cls, order, internals=True, logger=None, checkpoint=None, allowed_terms=None, allowed_coefficients=None, disallowed_coefficients=None, allowed_energy_changes=None, intermediate_normalization=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L62)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L62?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.modify_hamiltonian" class="docs-object-method">&nbsp;</a> 
```python
modify_hamiltonian(self, hamiltonian_corrections): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L126)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L126?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.get_correction" class="docs-object-method">&nbsp;</a> 
```python
get_correction(self, key, cls, order, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L145)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L145?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.shifted_hamiltonian_correction" class="docs-object-method">&nbsp;</a> 
```python
shifted_hamiltonian_correction(self, order, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L164)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L164?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.energy_correction" class="docs-object-method">&nbsp;</a> 
```python
energy_correction(self, order, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L167)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L167?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.wavefunction_correction" class="docs-object-method">&nbsp;</a> 
```python
wavefunction_correction(self, order, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L170)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L170?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.overlap_correction" class="docs-object-method">&nbsp;</a> 
```python
overlap_correction(self, order, degenerate_changes=None, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L173)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L173?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.full_wavefunction_correction" class="docs-object-method">&nbsp;</a> 
```python
full_wavefunction_correction(self, order, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L180)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L180?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.operator_correction" class="docs-object-method">&nbsp;</a> 
```python
operator_correction(self, order, operator_type=None, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L183)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L183?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.operator_degenerate_correction" class="docs-object-method">&nbsp;</a> 
```python
operator_degenerate_correction(self, order, /, degenerate_changes, operator_type=None, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L186)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L186?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.reexpressed_hamiltonian" class="docs-object-method">&nbsp;</a> 
```python
reexpressed_hamiltonian(self, order, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L191)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L191?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.reexpressed_hamiltonian_degenerate_correction" class="docs-object-method">&nbsp;</a> 
```python
reexpressed_hamiltonian_degenerate_correction(self, order, /, degenerate_changes, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L193)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L193?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.operator_expansion_terms" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
operator_expansion_terms(cls, order, logger=None, base_index=None, operator_type=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L198)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L198?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.clear_caches" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
clear_caches(cls): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L249)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L249?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Analytic.py#L41?message=Update%20Docs)   
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