## <a id="Psience.VPT2.Conditional.IntegerLatticeCondition">IntegerLatticeCondition</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Conditional.py#L290)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Conditional.py#L290?message=Update%20Docs)]
</div>

A displacement belongs to the integer lattice of resonance relations.

Hermite reduction is done once symbolically. Runtime predicates use exact
integer divisibility and reconstruction, not a floating null-space test.
The family is valid at every occupation with nonnegative endpoints.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
symmetric: bool
diagonal: bool
equivalence: bool
may_retain_change: _lru_cache_wrapper
```
<a id="Psience.VPT2.Conditional.IntegerLatticeCondition.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, relations): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Conditional.py#L300)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Conditional.py#L300?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Conditional.IntegerLatticeCondition.from_pairs" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_pairs(cls, pairs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L335)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L335?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Conditional.IntegerLatticeCondition.evaluate_shifted" class="docs-object-method">&nbsp;</a> 
```python
evaluate_shifted(self, state, permutations, start, end, frequencies): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Conditional/IntegerLatticeCondition.py#L341)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Conditional/IntegerLatticeCondition.py#L341?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Conditional.IntegerLatticeCondition.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Conditional/IntegerLatticeCondition.py#L350)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Conditional/IntegerLatticeCondition.py#L350?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/master/ci/docs/examples/Psience/VPT2/Conditional/IntegerLatticeCondition.md)/[New](https://github.com/McCoyGroup/Psience/new/master/?filename=ci/docs/examples/Psience/VPT2/Conditional/IntegerLatticeCondition.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/ci/docs/templates/Psience/VPT2/Conditional/IntegerLatticeCondition.md)/[New](https://github.com/McCoyGroup/Psience/new/master/?filename=ci/docs/templates/Psience/VPT2/Conditional/IntegerLatticeCondition.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Conditional.py#L290?message=Update%20Docs)   
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