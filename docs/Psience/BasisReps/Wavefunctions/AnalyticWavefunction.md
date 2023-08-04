## <a id="Psience.BasisReps.Wavefunctions.AnalyticWavefunction">AnalyticWavefunction</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Wavefunctions.py#L41)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Wavefunctions.py#L41?message=Update%20Docs)]
</div>

Little extension to RepresentationBasis so that we can use p and x and stuff
to evaluate out matrix elements and stuff.
This will be more progressively more tightly meshed with RepresentationBasis in the future,
but for now I just want to provide the scaffolding so that I can build off of it.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.BasisReps.Wavefunctions.AnalyticWavefunction.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, energy, data, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Wavefunctions/AnalyticWavefunction.py#L48)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Wavefunctions/AnalyticWavefunction.py#L48?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Wavefunctions.AnalyticWavefunction.evaluate" class="docs-object-method">&nbsp;</a> 
```python
evaluate(self, *args, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Wavefunctions/AnalyticWavefunction.py#L51)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Wavefunctions/AnalyticWavefunction.py#L51?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.Wavefunctions.AnalyticWavefunction.plot" class="docs-object-method">&nbsp;</a> 
```python
plot(self, figure=None, plot_class=None, domain=(-5, 5), **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Wavefunctions/AnalyticWavefunction.py#L54)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Wavefunctions/AnalyticWavefunction.py#L54?message=Update%20Docs)]
</div>
Uses McUtils to plot the wavefunction on the passed figure (makes a new one if none)
  - `figure`: `Graphics | Graphics3D`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.Wavefunctions.AnalyticWavefunction.expect" class="docs-object-method">&nbsp;</a> 
```python
expect(self, operator): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Wavefunctions/AnalyticWavefunction.py#L73)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Wavefunctions/AnalyticWavefunction.py#L73?message=Update%20Docs)]
</div>
Provides expectation values of operators, but the operators have to be Operator objects...
basically all the logic is inside the operator, but this is worth it for use in ExpansionWavefunction
We can also potentially add support for ExpansionOperators or SymbolicOperators in the future that are
able to very cleanly reuse stuff like the `p` matrix that a RepresentationBasis defines
  - `operator`: `Operator`
    > the operator to take the expectation of


<a id="Psience.BasisReps.Wavefunctions.AnalyticWavefunction.expectation" class="docs-object-method">&nbsp;</a> 
```python
expectation(self, op, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Wavefunctions/AnalyticWavefunction.py#L85)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Wavefunctions/AnalyticWavefunction.py#L85?message=Update%20Docs)]
</div>
Computes the expectation value of operator op over the wavefunction other and self
  - `other`: `AnalyticWavefunction`
    > the other wavefunction
  - `op`: `Operator`
    > the operator to take the matrix element of
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.Wavefunctions.AnalyticWavefunction.probability_density" class="docs-object-method">&nbsp;</a> 
```python
probability_density(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Wavefunctions/AnalyticWavefunction.py#L101)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Wavefunctions/AnalyticWavefunction.py#L101?message=Update%20Docs)]
</div>
Computes the probability density of the current wavefunction
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.Wavefunctions.AnalyticWavefunction.marginalize_out" class="docs-object-method">&nbsp;</a> 
```python
marginalize_out(self, dofs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/Wavefunctions/AnalyticWavefunction.py#L109)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Wavefunctions/AnalyticWavefunction.py#L109?message=Update%20Docs)]
</div>
Computes the projection of the current wavefunction onto a set of degrees
of freedom
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/Wavefunctions/AnalyticWavefunction.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/Wavefunctions/AnalyticWavefunction.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/Wavefunctions/AnalyticWavefunction.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/Wavefunctions/AnalyticWavefunction.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/Wavefunctions.py#L41?message=Update%20Docs)   
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