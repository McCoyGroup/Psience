## <a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorProductBasis">HarmonicOscillatorProductBasis</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator.py#L163)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator.py#L163?message=Update%20Docs)]
</div>

Tiny, tiny layer on `SimpleProductBasis` that makes use of some analytic work done
to support representations of `p` and `x`.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
nquant_max: int
```
<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorProductBasis.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, n_quanta, indexer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorProductBasis.py#L171)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorProductBasis.py#L171?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorProductBasis.__eq__" class="docs-object-method">&nbsp;</a> 
```python
__eq__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorProductBasis.py#L176)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorProductBasis.py#L176?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorProductBasis.to_state" class="docs-object-method">&nbsp;</a> 
```python
to_state(self, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorProductBasis.py#L182)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorProductBasis.py#L182?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorProductBasis.from_state" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_state(cls, data, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L187)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L187?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorProductBasis.operator" class="docs-object-method">&nbsp;</a> 
```python
operator(self, *terms, coeffs=None, axes=None, parallelizer=None, logger=None, chunk_size=None, **operator_settings): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorProductBasis.py#L191)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorProductBasis.py#L191?message=Update%20Docs)]
</div>
Builds an operator based on supplied terms, remapping names where possible.
If `coeffs` are supplied, a `ContractedOperator` is built.
  - `terms`: `Any`
    > 
  - `coeffs`: `Any`
    > 
  - `axes`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorProductBasis.take_subdimensions" class="docs-object-method">&nbsp;</a> 
```python
take_subdimensions(self, dims): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorProductBasis.py#L237)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorProductBasis.py#L237?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorProductBasis.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/HarmonicOscillator/HarmonicOscillatorProductBasis.py#L241)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator/HarmonicOscillatorProductBasis.py#L241?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/HarmonicOscillator/HarmonicOscillatorProductBasis.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/HarmonicOscillator/HarmonicOscillatorProductBasis.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/HarmonicOscillator/HarmonicOscillatorProductBasis.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/HarmonicOscillator/HarmonicOscillatorProductBasis.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/HarmonicOscillator.py#L163?message=Update%20Docs)   
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