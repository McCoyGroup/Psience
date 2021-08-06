## <a id="Psience.BasisReps.Wavefunctions.AnalyticWavefunction">AnalyticWavefunction</a>
Little extension to RepresentationBasis so that we can use p and x and stuff
to evaluate out matrix elements and stuff.
This will be more progressively more tightly meshed with RepresentationBasis in the future,
but for now I just want to provide the scaffolding so that I can build off of it.

### Properties and Methods
<a id="Psience.BasisReps.Wavefunctions.AnalyticWavefunction.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, energy, data, **opts): 
```

<a id="Psience.BasisReps.Wavefunctions.AnalyticWavefunction.evaluate" class="docs-object-method">&nbsp;</a>
```python
evaluate(self, *args, **kwargs): 
```

<a id="Psience.BasisReps.Wavefunctions.AnalyticWavefunction.plot" class="docs-object-method">&nbsp;</a>
```python
plot(self, figure=None, plot_class=None, domain=(-5, 5), **opts): 
```
Uses McUtils to plot the wavefunction on the passed figure (makes a new one if none)
- `figure`: `Graphics | Graphics3D`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Wavefunctions.AnalyticWavefunction.expect" class="docs-object-method">&nbsp;</a>
```python
expect(self, operator): 
```
Provides expectation values of operators, but the operators have to be Operator objects...
          basically all the logic is inside the operator, but this is worth it for use in ExpansionWavefunction
        We can also potentially add support for ExpansionOperators or SymbolicOperators in the future that are
          able to very cleanly reuse stuff like the `p` matrix that a RepresentationBasis defines
- `operator`: `Operator`
    >the operator to take the expectation of

<a id="Psience.BasisReps.Wavefunctions.AnalyticWavefunction.expectation" class="docs-object-method">&nbsp;</a>
```python
expectation(self, op, other): 
```
Computes the expectation value of operator op over the wavefunction other and self
- `other`: `AnalyticWavefunction`
    >the other wavefunction
- `op`: `Operator`
    >the operator to take the matrix element of
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Wavefunctions.AnalyticWavefunction.probability_density" class="docs-object-method">&nbsp;</a>
```python
probability_density(self): 
```
Computes the probability density of the current wavefunction
- `:returns`: `_`
    >No description...

### Examples


___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/edit/ci/examples/ci/docs/Psience/BasisReps/Wavefunctions/AnalyticWavefunction.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/examples/ci/docs/Psience/BasisReps/Wavefunctions/AnalyticWavefunction.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/edit/ci/docs/ci/docs/Psience/BasisReps/Wavefunctions/AnalyticWavefunction.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/docs/templates/ci/docs/Psience/BasisReps/Wavefunctions/AnalyticWavefunction.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Wavefunctions.py?message=Update%20Docs)