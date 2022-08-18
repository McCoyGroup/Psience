## <a id="Psience.BasisReps.ClassicalBases.ChebyshevBasis">ChebyshevBasis</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/ClassicalBases.py#L45)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/ClassicalBases.py#L45?message=Update%20Docs)]
</div>



<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

```python
name: str
selection_rules_mapping: dict
```
<a id="Psience.BasisReps.ClassicalBases.ChebyshevBasis.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, n_quanta): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/ClassicalBases.py#L47)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/ClassicalBases.py#L47?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.ClassicalBases.ChebyshevBasis.__eq__" class="docs-object-method">&nbsp;</a> 
```python
__eq__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/ClassicalBases.py#L49)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/ClassicalBases.py#L49?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.ClassicalBases.ChebyshevBasis.p" class="docs-object-method">&nbsp;</a> 
```python
p(self, n): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/ClassicalBases.py#L58)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/ClassicalBases.py#L58?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.ClassicalBases.ChebyshevBasis.p2" class="docs-object-method">&nbsp;</a> 
```python
p2(self, n): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/ClassicalBases.py#L61)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/ClassicalBases.py#L61?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/ClassicalBases.py#L85)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/ClassicalBases.py#L85?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.ClassicalBases.ChebyshevBasis.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/ClassicalBases.py#L93)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/ClassicalBases.py#L93?message=Update%20Docs)]
</div>

 </div>
</div>






___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/ClassicalBases/ChebyshevBasis.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/ClassicalBases/ChebyshevBasis.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/ClassicalBases/ChebyshevBasis.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/ClassicalBases/ChebyshevBasis.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/ClassicalBases.py#L45?message=Update%20Docs)