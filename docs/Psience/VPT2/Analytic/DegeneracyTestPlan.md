## <a id="Psience.VPT2.Analytic.DegeneracyTestPlan">DegeneracyTestPlan</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic.py#L117)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic.py#L117?message=Update%20Docs)]
</div>

Array-oriented union of mode-pattern degeneracy predicates.

Negative pattern entries are wildcards.  The object remains callable for
the legacy evaluator, while indexed evaluation can test a complete
permutation block with :meth:`evaluate` and avoid Python predicate/tree
traversal for each state.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
default_compare_workspace_bytes: int
default_scalar_state_cutoff: int
```
<a id="Psience.VPT2.Analytic.DegeneracyTestPlan.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, patterns=None, checks=None, fallback_tests=(), compare_workspace_bytes=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic.py#L129)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic.py#L129?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.DegeneracyTestPlan.cache_key" class="docs-object-method">&nbsp;</a> 
```python
@property
cache_key(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/DegeneracyTestPlan.py#L174)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/DegeneracyTestPlan.py#L174?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.DegeneracyTestPlan.evaluate" class="docs-object-method">&nbsp;</a> 
```python
evaluate(self, states): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/DegeneracyTestPlan.py#L183)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/DegeneracyTestPlan.py#L183?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.DegeneracyTestPlan.evaluate_array" class="docs-object-method">&nbsp;</a> 
```python
evaluate_array(self, states): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/DegeneracyTestPlan.py#L229)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/DegeneracyTestPlan.py#L229?message=Update%20Docs)]
</div>
Evaluate an explicitly batched state/permutation pool.

Unlike :meth:`evaluate`, this never falls back to the scalar trie for a
small batch.  The indexed evaluator uses it after coalescing work-item
subsets into one shared permutation pool.


<a id="Psience.VPT2.Analytic.DegeneracyTestPlan.__call__" class="docs-object-method">&nbsp;</a> 
```python
__call__(self, state): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/DegeneracyTestPlan.py#L335)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/DegeneracyTestPlan.py#L335?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Analytic/DegeneracyTestPlan.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Analytic/DegeneracyTestPlan.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Analytic/DegeneracyTestPlan.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Analytic/DegeneracyTestPlan.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic.py#L117?message=Update%20Docs)   
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