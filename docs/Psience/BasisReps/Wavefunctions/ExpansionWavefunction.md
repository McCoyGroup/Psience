## <a id="Psience.BasisReps.Wavefunctions.ExpansionWavefunction">ExpansionWavefunction</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Wavefunctions.py#L123)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Wavefunctions.py#L123?message=Update%20Docs)]
</div>

Simple wave function that takes a set of expansion coefficients alongside its basis.
Technically this should be called a _linear expansion wave function_, but
that was too long for my taste.

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

<a id="Psience.BasisReps.Wavefunctions.ExpansionWavefunction.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, energy, coefficients, basis_wfns): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Wavefunctions.py#L129)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Wavefunctions.py#L129?message=Update%20Docs)]
</div>


- `energy`: `float`
    >energy of the wavefunction
- `coefficients`: `Iterable[float]`
    >expansion coefficients
- `basis_wfns`: `Wavefunctions`
    >basis functions for the expansion

<a id="Psience.BasisReps.Wavefunctions.ExpansionWavefunction.coeffs" class="docs-object-method">&nbsp;</a> 
```python
@property
coeffs(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Wavefunctions.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Wavefunctions.py#L?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Wavefunctions.ExpansionWavefunction.basis" class="docs-object-method">&nbsp;</a> 
```python
@property
basis(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Wavefunctions.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Wavefunctions.py#L?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Wavefunctions.ExpansionWavefunction.evaluate" class="docs-object-method">&nbsp;</a> 
```python
evaluate(self, *args, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Wavefunctions.py#L153)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Wavefunctions.py#L153?message=Update%20Docs)]
</div>

Evaluates the wavecfunction as any other linear expansion.
- `args`: `Any`
    >coordinates + any other args the basis takes
- `kwargs`: `Any`
    >any keyword arguments the basis takes
- `:returns`: `_`
    >values of the wavefunction

<a id="Psience.BasisReps.Wavefunctions.ExpansionWavefunction.expect" class="docs-object-method">&nbsp;</a> 
```python
expect(self, operator): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Wavefunctions.py#L166)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Wavefunctions.py#L166?message=Update%20Docs)]
</div>

Provides the expectation value of the operator `op`.
        Uses the basis to compute the reps and then expands with the expansion coeffs.
- `operator`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Wavefunctions.ExpansionWavefunction.expectation" class="docs-object-method">&nbsp;</a> 
```python
expectation(self, op, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Wavefunctions.py#L185)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Wavefunctions.py#L185?message=Update%20Docs)]
</div>

Computes the expectation value of operator `op` over the wavefunction `other` and `self`.
        **Note**: _the basis of `other`, `self`, and `op` are assumed to be the same_.
- `op`: `Operator`
    >an operator represented in the basis of the expansion
- `other`: `ExpansionWavefunction`
    >the other wavefunction to expand over
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Wavefunctions.ExpansionWavefunction.probability_density" class="docs-object-method">&nbsp;</a> 
```python
probability_density(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Wavefunctions.py#L204)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Wavefunctions.py#L204?message=Update%20Docs)]
</div>

Computes the probability density of the current wavefunction
- `:returns`: `_`
    >No description...

 </div>
</div>




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/Wavefunctions/ExpansionWavefunction.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/Wavefunctions/ExpansionWavefunction.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/Wavefunctions/ExpansionWavefunction.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/Wavefunctions/ExpansionWavefunction.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Wavefunctions.py#L123?message=Update%20Docs)