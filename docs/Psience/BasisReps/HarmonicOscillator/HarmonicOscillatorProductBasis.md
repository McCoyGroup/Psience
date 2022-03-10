## <a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorProductBasis">HarmonicOscillatorProductBasis</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/HarmonicOscillator.py#L152)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/HarmonicOscillator.py#L152?message=Update%20Docs)]
</div>

Tiny, tiny layer on `SimpleProductBasis` that makes use of some analytic work done
to support representations of `p` and `x`.

```python
nquant_max: int
```
<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorProductBasis.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, n_quanta, indexer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/HarmonicOscillator.py#L160)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/HarmonicOscillator.py#L160?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorProductBasis.__eq__" class="docs-object-method">&nbsp;</a> 
```python
__eq__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/HarmonicOscillator.py#L165)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/HarmonicOscillator.py#L165?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorProductBasis.to_state" class="docs-object-method">&nbsp;</a> 
```python
to_state(self, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/HarmonicOscillator.py#L171)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/HarmonicOscillator.py#L171?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorProductBasis.from_state" class="docs-object-method">&nbsp;</a> 
```python
from_state(data, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/HarmonicOscillator.py#L176)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/HarmonicOscillator.py#L176?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorProductBasis.operator" class="docs-object-method">&nbsp;</a> 
```python
operator(self, *terms, coeffs=None, axes=None, parallelizer=None, logger=None, chunk_size=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/HarmonicOscillator.py#L180)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/HarmonicOscillator.py#L180?message=Update%20Docs)]
</div>

Builds an operator based on supplied terms, remapping names where possible.
        If `coeffs` are supplied, a `ContractedOperator` is built.
- `terms`: `Any`
    >No description...
- `coeffs`: `Any`
    >No description...
- `axes`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorProductBasis.take_subdimensions" class="docs-object-method">&nbsp;</a> 
```python
take_subdimensions(self, dims): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/HarmonicOscillator.py#L222)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/HarmonicOscillator.py#L222?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorProductBasis.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/HarmonicOscillator.py#L226)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/HarmonicOscillator.py#L226?message=Update%20Docs)]
</div>



___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/ci/docs/Psience/BasisReps/HarmonicOscillator/HarmonicOscillatorProductBasis.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/ci/docs/Psience/BasisReps/HarmonicOscillator/HarmonicOscillatorProductBasis.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/ci/docs/Psience/BasisReps/HarmonicOscillator/HarmonicOscillatorProductBasis.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/ci/docs/Psience/BasisReps/HarmonicOscillator/HarmonicOscillatorProductBasis.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/HarmonicOscillator.py#L152?message=Update%20Docs)