## <a id="Psience.VPT2.Analytic.CompactDegeneracyTestPlan">CompactDegeneracyTestPlan</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic.py#L345)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic.py#L345?message=Update%20Docs)]
</div>

Union of linear resonance clauses without expanded pattern rows.

For one energy-change requirement, each compiled inverse mode map projects a
concrete permutation back into the sorted mode tuple held by
``DegeneracyChangeIndex``. Membership in that small source index is exactly
the test that the expanded pattern union performs. In particular, keeping
the actual mode maps preserves the current all-equal-quanta behavior.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
default_scalar_state_cutoff: int
```
<a id="Psience.VPT2.Analytic.CompactDegeneracyTestPlan.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, clauses, width, mode_index=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic.py#L357)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic.py#L357?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.CompactDegeneracyTestPlan.__call__" class="docs-object-method">&nbsp;</a> 
```python
__call__(self, state): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/CompactDegeneracyTestPlan.py#L426)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/CompactDegeneracyTestPlan.py#L426?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.CompactDegeneracyTestPlan.evaluate" class="docs-object-method">&nbsp;</a> 
```python
evaluate(self, states): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/CompactDegeneracyTestPlan.py#L440)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/CompactDegeneracyTestPlan.py#L440?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.CompactDegeneracyTestPlan.evaluate_array" class="docs-object-method">&nbsp;</a> 
```python
evaluate_array(self, states): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/CompactDegeneracyTestPlan.py#L451)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/CompactDegeneracyTestPlan.py#L451?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.CompactDegeneracyTestPlan.evaluate_with_value_cache" class="docs-object-method">&nbsp;</a> 
```python
evaluate_with_value_cache(self, states, cache): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/CompactDegeneracyTestPlan.py#L476)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/CompactDegeneracyTestPlan.py#L476?message=Update%20Docs)]
</div>
Reuse predicate results across pools with equal checked values.

Only the positions this plan inspects enter the key.  In particular,
different unchecked mode assignments share one result.  Large pools
keep the existing vectorized kernel to avoid Python row-key overhead.
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Analytic/CompactDegeneracyTestPlan.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Analytic/CompactDegeneracyTestPlan.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Analytic/CompactDegeneracyTestPlan.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Analytic/CompactDegeneracyTestPlan.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic.py#L345?message=Update%20Docs)   
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