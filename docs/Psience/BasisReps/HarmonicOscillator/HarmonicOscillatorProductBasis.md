## <a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorProductBasis">HarmonicOscillatorProductBasis</a>
Tiny, tiny layer on `SimpleProductBasis` that makes use of some analytic work done
to support representations of `p` and `x`.

### Properties and Methods
<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorProductBasis.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, n_quanta): 
```

<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorProductBasis.operator" class="docs-object-method">&nbsp;</a>
```python
operator(self, *terms, coeffs=None, axes=None): 
```
Builds an operator based on supplied terms, remapping names where possible.
        If `coeffs` or `axes` are supplied, a `ContractedOperator` is built.
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

