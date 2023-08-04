## <a id="Psience.Wavefun.Wavefunctions.Wavefunction">Wavefunction</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Wavefun/Wavefunctions.py#L21)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions.py#L21?message=Update%20Docs)]
</div>

Represents a single wavefunction object







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.Wavefun.Wavefunctions.Wavefunction.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, energy, data, parent=None, index=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Wavefun/Wavefunctions/Wavefunction.py#L23)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions/Wavefunction.py#L23?message=Update%20Docs)]
</div>


<a id="Psience.Wavefun.Wavefunctions.Wavefunction.get_dimension" class="docs-object-method">&nbsp;</a> 
```python
get_dimension(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Wavefun/Wavefunctions/Wavefunction.py#L30)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions/Wavefunction.py#L30?message=Update%20Docs)]
</div>


<a id="Psience.Wavefun.Wavefunctions.Wavefunction.ndim" class="docs-object-method">&nbsp;</a> 
```python
@property
ndim(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Wavefun/Wavefunctions/Wavefunction.py#L33)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions/Wavefunction.py#L33?message=Update%20Docs)]
</div>


<a id="Psience.Wavefun.Wavefunctions.Wavefunction.plot" class="docs-object-method">&nbsp;</a> 
```python
plot(self, figure=None, domain=None, grid=None, values=None, plot_points=100, index=0, scaling=1, shift=0, plotter=None, plot_density=False, zero_tol=1e-08, contour_levels=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Wavefun/Wavefunctions/Wavefunction.py#L37)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions/Wavefunction.py#L37?message=Update%20Docs)]
</div>
Plots a single wave function on the grid
  - `figure`: `Any`
    > 
  - `grid`: `Any`
    > 
  - `index`: `Any`
    > 
  - `scaling`: `Any`
    > 
  - `shift`: `Any`
    > 
  - `opts`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Wavefun.Wavefunctions.Wavefunction.projection_plot" class="docs-object-method">&nbsp;</a> 
```python
projection_plot(self, coords, figure=None, **plot_options): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Wavefun/Wavefunctions/Wavefunction.py#L113)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions/Wavefunction.py#L113?message=Update%20Docs)]
</div>
A convenience function to plot multiple projections
on the same set of axes
  - `coords`: `Any`
    > 
  - `figure`: `Any`
    > 
  - `plot_options`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Wavefun.Wavefunctions.Wavefunction.expectation" class="docs-object-method">&nbsp;</a> 
```python
expectation(self, op, other=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Wavefun/Wavefunctions/Wavefunction.py#L146)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions/Wavefunction.py#L146?message=Update%20Docs)]
</div>
Computes the expectation value of operator op over the wavefunction other and self
  - `other`: `Wavefunction`
    > 
  - `op`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Wavefun.Wavefunctions.Wavefunction.overlap" class="docs-object-method">&nbsp;</a> 
```python
overlap(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Wavefun/Wavefunctions/Wavefunction.py#L158)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions/Wavefunction.py#L158?message=Update%20Docs)]
</div>


<a id="Psience.Wavefun.Wavefunctions.Wavefunction.evaluate" class="docs-object-method">&nbsp;</a> 
```python
evaluate(self, points): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Wavefun/Wavefunctions/Wavefunction.py#L160)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions/Wavefunction.py#L160?message=Update%20Docs)]
</div>
Evaluates the current wavefunction
  - `:returns`: `_`
    >


<a id="Psience.Wavefun.Wavefunctions.Wavefunction.probability_density" class="docs-object-method">&nbsp;</a> 
```python
@property
probability_density(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Wavefun/Wavefunctions/Wavefunction.py#L169)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions/Wavefunction.py#L169?message=Update%20Docs)]
</div>
Computes the probability density of the current wavefunction
  - `:returns`: `_`
    >


<a id="Psience.Wavefun.Wavefunctions.Wavefunction.marginalize_out" class="docs-object-method">&nbsp;</a> 
```python
marginalize_out(self, dofs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Wavefun/Wavefunctions/Wavefunction.py#L178)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions/Wavefunction.py#L178?message=Update%20Docs)]
</div>
Integrates out the contributions from the degrees of freedom `dofs`
  - `:returns`: `Wavefunction`
    >


<a id="Psience.Wavefun.Wavefunctions.Wavefunction.project" class="docs-object-method">&nbsp;</a> 
```python
project(self, dofs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Wavefun/Wavefunctions/Wavefunction.py#L187)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions/Wavefunction.py#L187?message=Update%20Docs)]
</div>
Computes the projection of the current wavefunction onto a set of degrees
of freedom, returning a projected wave function object
  - `:returns`: `Wavefunction`
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Wavefun/Wavefunctions/Wavefunction.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Wavefun/Wavefunctions/Wavefunction.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Wavefun/Wavefunctions/Wavefunction.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Wavefun/Wavefunctions/Wavefunction.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions.py#L21?message=Update%20Docs)   
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