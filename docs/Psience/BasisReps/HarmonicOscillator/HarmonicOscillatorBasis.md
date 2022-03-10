## <a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorBasis">HarmonicOscillatorBasis</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/HarmonicOscillator.py#L20)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/HarmonicOscillator.py#L20?message=Update%20Docs)]
</div>

Provides a concrete implementation of RepresentationBasis using the H.O.
Need to make it handle units a bit better.
Currently 1D, need to make multidimensional in the future.

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
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/HarmonicOscillator.py#L27)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/HarmonicOscillator.py#L27?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorBasis.__eq__" class="docs-object-method">&nbsp;</a> 
```python
__eq__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/HarmonicOscillator.py#L35)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/HarmonicOscillator.py#L35?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorBasis.pmatrix_ho" class="docs-object-method">&nbsp;</a> 
```python
pmatrix_ho(n): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/HarmonicOscillator.py#L54)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/HarmonicOscillator.py#L54?message=Update%20Docs)]
</div>

There's one big subtlety to what we're doing here, which is that
          for efficiency reasons we return an entirely real matrix
        The reason for that is we assumed people would mostly use it in the context
          of stuff like pp, pQp, or pQQp, in which case the imaginary part pulls out
          and becomes a negative sign
        We actually use this assumption across _all_ of our representations
- `n`: `Any`
    >No description...
- `:returns`: `sp.csr_matrix`
    >No description...

<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorBasis.qmatrix_ho" class="docs-object-method">&nbsp;</a> 
```python
qmatrix_ho(n): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/HarmonicOscillator.py#L79)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/HarmonicOscillator.py#L79?message=Update%20Docs)]
</div>


- `n`: `Any`
    >No description...
- `:returns`: `sp.csr_matrix`
    >No description...

<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorBasis.operator" class="docs-object-method">&nbsp;</a> 
```python
operator(self, *terms, logger=None, parallelizer=None, chunk_size=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/HarmonicOscillator.py#L96)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/HarmonicOscillator.py#L96?message=Update%20Docs)]
</div>

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

<a id="Psience.BasisReps.HarmonicOscillator.HarmonicOscillatorBasis.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/HarmonicOscillator.py#L147)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/HarmonicOscillator.py#L147?message=Update%20Docs)]
</div>



___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/ci/docs/Psience/BasisReps/HarmonicOscillator/HarmonicOscillatorBasis.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/ci/docs/Psience/BasisReps/HarmonicOscillator/HarmonicOscillatorBasis.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/ci/docs/Psience/BasisReps/HarmonicOscillator/HarmonicOscillatorBasis.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/ci/docs/Psience/BasisReps/HarmonicOscillator/HarmonicOscillatorBasis.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/HarmonicOscillator.py#L20?message=Update%20Docs)