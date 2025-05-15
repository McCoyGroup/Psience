## <a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix">StateSpaceMatrix</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces.py#L4534)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces.py#L4534?message=Update%20Docs)]
</div>

A `SparseArray` that holds onto a `BasisStateSpace` that keeps track of the
total set of states involved.
By default is assumed real-symmetric. This can be relaxed in the future.

TODO: The idea is good, but calculating what is "in" the array and what is "out"
every single time this is applied could be slow...
We'll need to test to see how slow







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, initial_basis, initial_vals=None, column_space=None, symmetric=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4545)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4545?message=Update%20Docs)]
</div>

  - `initial_basis`: `BasisStateSpace | RepresentationBasis`
    >


<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.shape" class="docs-object-method">&nbsp;</a> 
```python
@property
shape(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4604)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4604?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.basis" class="docs-object-method">&nbsp;</a> 
```python
@property
basis(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4610)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4610?message=Update%20Docs)]
</div>
Returns the basis for the matrix rep
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.brakets" class="docs-object-method">&nbsp;</a> 
```python
@property
brakets(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4620)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4620?message=Update%20Docs)]
</div>
Returns the BraKetSpace for the held indices
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.identity_from_space" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
identity_from_space(cls, space, column_space=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L4631)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L4631?message=Update%20Docs)]
</div>
Returns a StateSpaceMatrix where the diagonal is filled with 1s
  - `space`: `Any`
    > 
  - `column_space`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.extend_basis" class="docs-object-method">&nbsp;</a> 
```python
extend_basis(self, states, extend_columns=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4669)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4669?message=Update%20Docs)]
</div>
Extends the held state space and resizes the held array if need be
  - `states`: `BasisStateSpace`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.compute_values" class="docs-object-method">&nbsp;</a> 
```python
compute_values(self, func, brakets): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4752)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4752?message=Update%20Docs)]
</div>
Computes new values into the held `SparseArray` based on the function and brakets provided
and returns the entire array of values
  - `func`: `Any`
    > A function that can take a braket spec and compute values
  - `brakets`: `Any`
    > A set of brakets to compute values for
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.dot" class="docs-object-method">&nbsp;</a> 
```python
dot(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4804)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4804?message=Update%20Docs)]
</div>
Performs a dot product between the held SparseArray and another
StateSpaceMatrix
  - `other`: `StateSpaceMatrix`
    > other matrix
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4840)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4840?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.__setitem__" class="docs-object-method">&nbsp;</a> 
```python
__setitem__(self, item, vals): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4851)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4851?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateSpaces.StateSpaceMatrix.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4862)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces/StateSpaceMatrix.py#L4862?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/StateSpaces/StateSpaceMatrix.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/StateSpaces/StateSpaceMatrix.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/StateSpaces/StateSpaceMatrix.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/StateSpaces/StateSpaceMatrix.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateSpaces.py#L4534?message=Update%20Docs)   
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