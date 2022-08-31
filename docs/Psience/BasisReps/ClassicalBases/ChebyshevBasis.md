## <a id="Psience.BasisReps.ClassicalBases.ChebyshevBasis">ChebyshevBasis</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/ClassicalBases.py#L45)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/ClassicalBases.py#L45?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
name: str
selection_rules_mapping: dict
```
<a id="Psience.BasisReps.ClassicalBases.ChebyshevBasis.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, n_quanta): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/ClassicalBases/ChebyshevBasis.py#L47)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/ClassicalBases/ChebyshevBasis.py#L47?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.ClassicalBases.ChebyshevBasis.__eq__" class="docs-object-method">&nbsp;</a> 
```python
__eq__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/ClassicalBases/ChebyshevBasis.py#L49)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/ClassicalBases/ChebyshevBasis.py#L49?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.ClassicalBases.ChebyshevBasis.p" class="docs-object-method">&nbsp;</a> 
```python
p(self, n): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/ClassicalBases/ChebyshevBasis.py#L58)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/ClassicalBases/ChebyshevBasis.py#L58?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.ClassicalBases.ChebyshevBasis.p2" class="docs-object-method">&nbsp;</a> 
```python
p2(self, n): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/ClassicalBases/ChebyshevBasis.py#L61)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/ClassicalBases/ChebyshevBasis.py#L61?message=Update%20Docs)]
</div>
0  -> (Which[# \[Equal] 1, 3, # \[Equal] 2, 1/2,
True, (-1 - 8 (# - 1) - 4 (# - 1)^2)]*1/8 &),
2  -> (If[# \[Equal] 1,
15/(8*Sqrt[2]), (15 + 16 (# - 1) + 4 (# - 1)^2)*1/16] &),
-2 -> (If[#2 \[Equal] 1,
15/(8*Sqrt[2]), (15 + 16 (#2 - 1) + 4 (#2 - 1)^2)*1/16] &)
:param n:
:type n:
:return:
:rtype:


<a id="Psience.BasisReps.ClassicalBases.ChebyshevBasis.x" class="docs-object-method">&nbsp;</a> 
```python
x(self, n): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/ClassicalBases/ChebyshevBasis.py#L85)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/ClassicalBases/ChebyshevBasis.py#L85?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.ClassicalBases.ChebyshevBasis.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/ClassicalBases/ChebyshevBasis.py#L93)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/ClassicalBases/ChebyshevBasis.py#L93?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/ClassicalBases/ChebyshevBasis.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/ClassicalBases/ChebyshevBasis.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/ClassicalBases/ChebyshevBasis.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/ClassicalBases/ChebyshevBasis.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/ClassicalBases.py#L45?message=Update%20Docs)   
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