## <a id="Psience.BasisReps.Wavefunctions.AnalyticWavefunction">AnalyticWavefunction</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Wavefunctions.py#L20)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Wavefunctions.py#L20?message=Update%20Docs)]
</div>

Little extension to RepresentationBasis so that we can use p and x and stuff
to evaluate out matrix elements and stuff.
This will be more progressively more tightly meshed with RepresentationBasis in the future,
but for now I just want to provide the scaffolding so that I can build off of it.

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

<a id="Psience.BasisReps.Wavefunctions.AnalyticWavefunction.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, energy, data, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Wavefunctions.py#L27)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Wavefunctions.py#L27?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Wavefunctions.AnalyticWavefunction.evaluate" class="docs-object-method">&nbsp;</a> 
```python
evaluate(self, *args, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Wavefunctions.py#L30)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Wavefunctions.py#L30?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Wavefunctions.AnalyticWavefunction.plot" class="docs-object-method">&nbsp;</a> 
```python
plot(self, figure=None, plot_class=None, domain=(-5, 5), **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Wavefunctions.py#L33)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Wavefunctions.py#L33?message=Update%20Docs)]
</div>

Uses McUtils to plot the wavefunction on the passed figure (makes a new one if none)
- `figure`: `Graphics | Graphics3D`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Wavefunctions.AnalyticWavefunction.expect" class="docs-object-method">&nbsp;</a> 
```python
expect(self, operator): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Wavefunctions.py#L52)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Wavefunctions.py#L52?message=Update%20Docs)]
</div>

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
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Wavefunctions.py#L64)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Wavefunctions.py#L64?message=Update%20Docs)]
</div>

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
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/BasisReps/Wavefunctions.py#L79)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Wavefunctions.py#L79?message=Update%20Docs)]
</div>

Computes the probability density of the current wavefunction
- `:returns`: `_`
    >No description...

 </div>
</div>




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/Wavefunctions/AnalyticWavefunction.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/Wavefunctions/AnalyticWavefunction.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/Wavefunctions/AnalyticWavefunction.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/Wavefunctions/AnalyticWavefunction.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/BasisReps/Wavefunctions.py#L20?message=Update%20Docs)