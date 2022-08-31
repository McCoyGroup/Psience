## <a id="Psience.BasisReps.Bases.RepresentationBasis">RepresentationBasis</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Bases.py#L18)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases.py#L18?message=Update%20Docs)]
</div>

Metaclass for representations.
Requires concrete implementations of the position and momentum operators.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
name: str
selection_rules_mapping: dict
```
<a id="Psience.BasisReps.Bases.RepresentationBasis.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, function_generator, n_quanta): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Bases/RepresentationBasis.py#L24)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases/RepresentationBasis.py#L24?message=Update%20Docs)]
</div>

  - `function_generator`: `Any`
    > 
  - `n_quanta`: `int`
    > numbers of quanta (hold over from initial implementation)


<a id="Psience.BasisReps.Bases.RepresentationBasis.__eq__" class="docs-object-method">&nbsp;</a> 
```python
__eq__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Bases/RepresentationBasis.py#L35)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases/RepresentationBasis.py#L35?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Bases.RepresentationBasis.dimensions" class="docs-object-method">&nbsp;</a> 
```python
@property
dimensions(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Bases/RepresentationBasis.py#L39)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases/RepresentationBasis.py#L39?message=Update%20Docs)]
</div>
Returns the dimensions of the basis
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.Bases.RepresentationBasis.ndim" class="docs-object-method">&nbsp;</a> 
```python
@property
ndim(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Bases/RepresentationBasis.py#L49)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases/RepresentationBasis.py#L49?message=Update%20Docs)]
</div>
Returns the number of dimensions of the basis
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.Bases.RepresentationBasis.ravel_state_inds" class="docs-object-method">&nbsp;</a> 
```python
ravel_state_inds(self, idx): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Bases/RepresentationBasis.py#L59)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases/RepresentationBasis.py#L59?message=Update%20Docs)]
</div>
Converts state indices from an array of quanta to an array of indices...except in 1D this really isn't doing anything
  - `idx`: `Iterable[Iterable[int]]`
    > indices
  - `:returns`: `tuple[int]`
    > a
r
r
a
y
 
o
f
 
s
t
a
t
e
 
i
n
d
i
c
e
s
 
i
n
 
t
h
e
 
b
a
s
i
s


<a id="Psience.BasisReps.Bases.RepresentationBasis.unravel_state_inds" class="docs-object-method">&nbsp;</a> 
```python
unravel_state_inds(self, idx): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Bases/RepresentationBasis.py#L76)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases/RepresentationBasis.py#L76?message=Update%20Docs)]
</div>
Converts state indices from an array of ints to an array of quanta...except in 1D this really isn't doing anything
  - `idx`: `Iterable[int]`
    > indices
  - `:returns`: `tuple[tuple[int]]`
    > a
r
r
a
y
 
o
f
 
s
t
a
t
e
 
t
u
p
l
e
s
 
i
n
 
t
h
e
 
b
a
s
i
s


<a id="Psience.BasisReps.Bases.RepresentationBasis.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Bases/RepresentationBasis.py#L93)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases/RepresentationBasis.py#L93?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Bases.RepresentationBasis.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Bases/RepresentationBasis.py#L97)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases/RepresentationBasis.py#L97?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Bases.RepresentationBasis.p" class="docs-object-method">&nbsp;</a> 
```python
p(self, n): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Bases/RepresentationBasis.py#L105)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases/RepresentationBasis.py#L105?message=Update%20Docs)]
</div>
Generates the momentum matrix up to n-quanta.
There's one big subtlety to what we're doing here, which is that
for efficiency reasons we return an entirely real matrix
The reason for that is we assumed people would mostly use it in the context
of stuff like pp, pQp, or pQQp, in which case the imaginary part pulls out
and becomes a negative sign
We actually use this assumption across _all_ of our representations.
  - `n`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.Bases.RepresentationBasis.x" class="docs-object-method">&nbsp;</a> 
```python
x(self, n): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Bases/RepresentationBasis.py#L123)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases/RepresentationBasis.py#L123?message=Update%20Docs)]
</div>
Generates the coordinate matrix up to n-quanta
  - `n`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.Bases.RepresentationBasis.I" class="docs-object-method">&nbsp;</a> 
```python
I(self, n): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Bases/RepresentationBasis.py#L135)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases/RepresentationBasis.py#L135?message=Update%20Docs)]
</div>
Generates the identity matrix up to n-quanta
  - `n`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.Bases.RepresentationBasis.operator_mapping" class="docs-object-method">&nbsp;</a> 
```python
@property
operator_mapping(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Bases/RepresentationBasis.py#L146)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases/RepresentationBasis.py#L146?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Bases.RepresentationBasis.operator" class="docs-object-method">&nbsp;</a> 
```python
operator(self, *terms, logger=None, parallelizer=None, chunk_size=None, **operator_settings): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Bases/RepresentationBasis.py#L151)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases/RepresentationBasis.py#L151?message=Update%20Docs)]
</div>
Provides an `Operator` to handle the given terms
  - `terms`: `Any`
    > 
  - `logger`: `Any`
    > 
  - `parallelizer`: `Any`
    > 
  - `chunk_size`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.Bases.RepresentationBasis.representation" class="docs-object-method">&nbsp;</a> 
```python
representation(self, *terms, logger=None, name=None, parallelizer=None, chunk_size=None, memory_constrained=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Bases/RepresentationBasis.py#L175)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases/RepresentationBasis.py#L175?message=Update%20Docs)]
</div>
Provides a representation of a product operator specified by `terms`
  - `terms`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.Bases.RepresentationBasis.selection_rule_steps" class="docs-object-method">&nbsp;</a> 
```python
selection_rule_steps(self, *terms): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Bases/RepresentationBasis.py#L272)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases/RepresentationBasis.py#L272?message=Update%20Docs)]
</div>
Generates the full set of possible selection rules for terms
  - `terms`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.Bases.RepresentationBasis.selection_rules" class="docs-object-method">&nbsp;</a> 
```python
selection_rules(self, *terms): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Bases/RepresentationBasis.py#L294)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases/RepresentationBasis.py#L294?message=Update%20Docs)]
</div>
Generates the full set of possible selection rules for terms
  - `terms`: `Any`
    > 
  - `:returns`: `_`
    >
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/Bases/RepresentationBasis.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/Bases/RepresentationBasis.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/Bases/RepresentationBasis.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/Bases/RepresentationBasis.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Bases.py#L18?message=Update%20Docs)   
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