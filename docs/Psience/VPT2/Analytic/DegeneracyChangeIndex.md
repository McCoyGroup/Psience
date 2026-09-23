## <a id="Psience.VPT2.Analytic.DegeneracyChangeIndex">DegeneracyChangeIndex</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic.py#L529)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic.py#L529?message=Update%20Docs)]
</div>

An ordered, immutable index of resonance changes by sorted quanta.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
canonical_mode_cache_max_bytes: int
```
<a id="Psience.VPT2.Analytic.DegeneracyChangeIndex.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, changes): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic.py#L534)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic.py#L534?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.DegeneracyChangeIndex.get" class="docs-object-method">&nbsp;</a> 
```python
get(self, signature): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/DegeneracyChangeIndex.py#L563)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/DegeneracyChangeIndex.py#L563?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.DegeneracyChangeIndex.get_canonical_mode_rows" class="docs-object-method">&nbsp;</a> 
```python
get_canonical_mode_rows(self, signature, ordering): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/DegeneracyChangeIndex.py#L566)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/DegeneracyChangeIndex.py#L566?message=Update%20Docs)]
</div>
Canonical sparse row bytes, shared across checked-position maps.

The bounded cache holds only distinct mode-tuple reorderings.  Thus
interning equivalent predicates does not reconstruct the full
resonance-mode x projection-map Cartesian product for every plan.


<a id="Psience.VPT2.Analytic.DegeneracyChangeIndex.get_mode_membership" class="docs-object-method">&nbsp;</a> 
```python
get_mode_membership(self, signature): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/DegeneracyChangeIndex.py#L594)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/DegeneracyChangeIndex.py#L594?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Analytic/DegeneracyChangeIndex.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Analytic/DegeneracyChangeIndex.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Analytic/DegeneracyChangeIndex.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Analytic/DegeneracyChangeIndex.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic.py#L529?message=Update%20Docs)   
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