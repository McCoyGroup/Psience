## <a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorProductBasis">HarmonicOscillatorProductBasis</a>
Tiny, tiny layer on `SimpleProductBasis` that makes use of some analytic work done
to support representations of `p` and `x`.

### Properties and Methods
```python
nquant_max: int
from_state: method
```
<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorProductBasis.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, n_quanta, indexer=None): 
```

<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorProductBasis.__eq__" class="docs-object-method">&nbsp;</a>
```python
__eq__(self, other): 
```

<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorProductBasis.to_state" class="docs-object-method">&nbsp;</a>
```python
to_state(self, serializer=None): 
```

<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorProductBasis.operator" class="docs-object-method">&nbsp;</a>
```python
operator(self, *terms, coeffs=None, axes=None, parallelizer=None, logger=None, chunk_size=None): 
```
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

<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorProductBasis.__repr__" class="docs-object-method">&nbsp;</a>
```python
__repr__(self): 
```

### Examples


___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/edit/ci/examples/ci/docs/Psience/BasisReps/HarmonicOscillator/HarmonicOscillatorProductBasis.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/examples/ci/docs/Psience/BasisReps/HarmonicOscillator/HarmonicOscillatorProductBasis.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/edit/ci/docs/ci/docs/Psience/BasisReps/HarmonicOscillator/HarmonicOscillatorProductBasis.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/docs/templates/ci/docs/Psience/BasisReps/HarmonicOscillator/HarmonicOscillatorProductBasis.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/HarmonicOscillator.py?message=Update%20Docs)