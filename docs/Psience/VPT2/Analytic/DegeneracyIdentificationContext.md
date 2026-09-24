## <a id="Psience.VPT2.Analytic.DegeneracyIdentificationContext">DegeneracyIdentificationContext</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic.py#L631)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic.py#L631?message=Update%20Docs)]
</div>

Compile and bind energy-change predicates within a bounded batch cache.

Only symbolic requirements and predicate plans are retained.  Concrete
state/permutation masks belong to the indexed evaluator's separate batch
cache.  Direct callers can create a short-lived context; high-level
correction workflows share one across their expressions and DAG chunks.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
default_max_items: int
default_max_bytes: int
use_compact_index: bool
```
<a id="Psience.VPT2.Analytic.DegeneracyIdentificationContext.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, changes, max_items=None, max_bytes=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic.py#L645)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic.py#L645?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.DegeneracyIdentificationContext.bind" class="docs-object-method">&nbsp;</a> 
```python
bind(self, expr, nmodes, method): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/DegeneracyIdentificationContext.py#L903)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/DegeneracyIdentificationContext.py#L903?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.DegeneracyIdentificationContext.stats" class="docs-object-method">&nbsp;</a> 
```python
stats(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/DegeneracyIdentificationContext.py#L923)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/DegeneracyIdentificationContext.py#L923?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Analytic/DegeneracyIdentificationContext.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Analytic/DegeneracyIdentificationContext.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Analytic/DegeneracyIdentificationContext.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Analytic/DegeneracyIdentificationContext.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic.py#L631?message=Update%20Docs)   
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