## <a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorBasis">HarmonicOscillatorBasis</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator.py#L25)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator.py#L25?message=Update%20Docs)]
</div>

Provides a concrete implementation of RepresentationBasis using the H.O.
Need to make it handle units a bit better.
Currently 1D, need to make multidimensional in the future.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
name: str
selection_rules_mapping: dict
p: _lru_cache_wrapper
x: _lru_cache_wrapper
```
<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorBasis.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, n_quanta, m=None, re=None, dimensionless=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorBasis.py#L32)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorBasis.py#L32?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorBasis.__eq__" class="docs-object-method">&nbsp;</a> 
```python
__eq__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorBasis.py#L40)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorBasis.py#L40?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorBasis.pmatrix_ho" class="docs-object-method">&nbsp;</a> 
```python
@staticmethod
pmatrix_ho(n): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L59)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L59?message=Update%20Docs)]
</div>
There's one big subtlety to what we're doing here, which is that
for efficiency reasons we return an entirely real matrix
The reason for that is we assumed people would mostly use it in the context
of stuff like pp, pQp, or pQQp, in which case the imaginary part pulls out
and becomes a negative sign
We actually use this assumption across _all_ of our representations
  - `n`: `Any`
    > 
  - `:returns`: `sp.csr_matrix`
    >


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorBasis.qmatrix_ho" class="docs-object-method">&nbsp;</a> 
```python
@staticmethod
qmatrix_ho(n): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L84)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L84?message=Update%20Docs)]
</div>

  - `n`: `Any`
    > 
  - `:returns`: `sp.csr_matrix`
    >


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorBasis.operator" class="docs-object-method">&nbsp;</a> 
```python
operator(self, *terms, logger=None, parallelizer=None, chunk_size=None, **operator_settings): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorBasis.py#L101)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorBasis.py#L101?message=Update%20Docs)]
</div>
Builds an operator based on supplied terms, remapping names where possible.
If `coeffs` or `axes` are supplied, a `ContractedOperator` is built.
  - `terms`: `Any`
    > 
  - `coeffs`: `Any`
    > 
  - `axes`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorBasis.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorBasis.py#L158)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorBasis.py#L158?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/HarmonicOscillator/HarmonicOscillatorBasis.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/HarmonicOscillator/HarmonicOscillatorBasis.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/HarmonicOscillator/HarmonicOscillatorBasis.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/HarmonicOscillator/HarmonicOscillatorBasis.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator.py#L25?message=Update%20Docs)   
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