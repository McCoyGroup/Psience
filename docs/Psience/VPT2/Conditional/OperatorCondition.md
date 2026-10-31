## <a id="Psience.VPT2.Conditional.OperatorCondition">OperatorCondition</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Conditional.py#L60)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Conditional.py#L60?message=Update%20Docs)]
</div>

A serializable relation on the two occupations of an operator step.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
symmetric: bool
diagonal: NoneType
equivalence: bool
requires_frequencies: bool
```
<a id="Psience.VPT2.Conditional.OperatorCondition.__setattr__" class="docs-object-method">&nbsp;</a> 
```python
__setattr__(self, name, value): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Conditional/OperatorCondition.py#L67)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Conditional/OperatorCondition.py#L67?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Conditional.OperatorCondition.__hash__" class="docs-object-method">&nbsp;</a> 
```python
__hash__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Conditional/OperatorCondition.py#L72)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Conditional/OperatorCondition.py#L72?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Conditional.OperatorCondition.__eq__" class="docs-object-method">&nbsp;</a> 
```python
__eq__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Conditional/OperatorCondition.py#L75)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Conditional/OperatorCondition.py#L75?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Conditional.OperatorCondition.__and__" class="docs-object-method">&nbsp;</a> 
```python
__and__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Conditional/OperatorCondition.py#L80)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Conditional/OperatorCondition.py#L80?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Conditional.OperatorCondition.__or__" class="docs-object-method">&nbsp;</a> 
```python
__or__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Conditional/OperatorCondition.py#L83)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Conditional/OperatorCondition.py#L83?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Conditional.OperatorCondition.__invert__" class="docs-object-method">&nbsp;</a> 
```python
__invert__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Conditional/OperatorCondition.py#L86)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Conditional/OperatorCondition.py#L86?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Conditional.OperatorCondition.symmetrized" class="docs-object-method">&nbsp;</a> 
```python
symmetrized(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Conditional/OperatorCondition.py#L89)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Conditional/OperatorCondition.py#L89?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Conditional.OperatorCondition.evaluate_shifted" class="docs-object-method">&nbsp;</a> 
```python
evaluate_shifted(self, state, permutations, start, end, frequencies): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Conditional/OperatorCondition.py#L92)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Conditional/OperatorCondition.py#L92?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Conditional.OperatorCondition.matches" class="docs-object-method">&nbsp;</a> 
```python
matches(self, initial, final, frequencies=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Conditional/OperatorCondition.py#L95)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Conditional/OperatorCondition.py#L95?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Conditional.OperatorCondition.matches" class="docs-object-method">&nbsp;</a> 
```python
__call__(self, initial, final, frequencies=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Conditional/OperatorCondition.py#L95)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Conditional/OperatorCondition.py#L95?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Conditional.OperatorCondition.may_retain_change" class="docs-object-method">&nbsp;</a> 
```python
may_retain_change(self, change): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Conditional/OperatorCondition.py#L109)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Conditional/OperatorCondition.py#L109?message=Update%20Docs)]
</div>
False proves a change impossible; None leaves it unresolved.
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
[Edit](https://github.com/McCoyGroup/Psience/edit/master/ci/docs/examples/Psience/VPT2/Conditional/OperatorCondition.md)/[New](https://github.com/McCoyGroup/Psience/new/master/?filename=ci/docs/examples/Psience/VPT2/Conditional/OperatorCondition.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/ci/docs/templates/Psience/VPT2/Conditional/OperatorCondition.md)/[New](https://github.com/McCoyGroup/Psience/new/master/?filename=ci/docs/templates/Psience/VPT2/Conditional/OperatorCondition.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Conditional.py#L60?message=Update%20Docs)   
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