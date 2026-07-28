## <a id="Psience.DVR.Wavefunctions.DVRWavefunction">DVRWavefunction</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Wavefunctions.py#L15)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Wavefunctions.py#L15?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.DVR.Wavefunctions.DVRWavefunction.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, energy, data, parent=None, grid=None, index=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Wavefunctions.py#L17)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Wavefunctions.py#L17?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a single DVR wavefunction from its value at the DVR grid points, falling back to the parent collection's grid if none is given directly.
  - `energy`: `float`
    > the wavefunction's energy
  - `data`: `np.ndarray`
    > the wavefunction's values at the grid points
  - `parent`: `DVRWavefunctions | None`
    > the `DVRWavefunctions` collection this wavefunction belongs to
  - `grid`: `np.ndarray | None`
    > the DVR grid this wavefunction is defined on; defaults to `parent.grid`
  - `index`: `object | None`
    > this wavefunction's index within its parent collection
  - `opts`: `dict`
    > extra options forwarded to the base `Wavefunction.__init__`
  - `:returns`: `None`
    > None


<a id="Psience.DVR.Wavefunctions.DVRWavefunction.get_dimension" class="docs-object-method">&nbsp;</a> 
```python
get_dimension(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Wavefunctions/DVRWavefunction.py#L43)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Wavefunctions/DVRWavefunction.py#L43?message=Update%20Docs)]
</div>
**LLM Docstring**

The number of degrees of freedom this wavefunction is defined over, inferred from the trailing dimension of its grid.
  - `:returns`: `int`
    > the dimensionality


<a id="Psience.DVR.Wavefunctions.DVRWavefunction.plot" class="docs-object-method">&nbsp;</a> 
```python
plot(self, figure=None, grid=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Wavefunctions/DVRWavefunction.py#L54)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Wavefunctions/DVRWavefunction.py#L54?message=Update%20Docs)]
</div>
**LLM Docstring**

Plot the wavefunction using its own DVR grid and values, delegating to the base `Wavefunction.plot`.
  - `figure`: `object | None`
    > an existing figure to draw into
  - `grid`: `np.ndarray | None`
    > the grid to plot over; defaults to `self.grid`
  - `opts`: `dict`
    > extra options forwarded to the base `plot`
  - `:returns`: `object`
    > the resulting figure


<a id="Psience.DVR.Wavefunctions.DVRWavefunction.expectation" class="docs-object-method">&nbsp;</a> 
```python
expectation(self, op, other=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Wavefunctions/DVRWavefunction.py#L73)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Wavefunctions/DVRWavefunction.py#L73?message=Update%20Docs)]
</div>
Computes the expectation value of operator op over the wavefunction other and self
  - `other`: `Wavefunction | np.ndarray`
    > 
  - `op`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.DVR.Wavefunctions.DVRWavefunction.interp" class="docs-object-method">&nbsp;</a> 
```python
@property
interp(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Wavefunctions/DVRWavefunction.py#L88)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Wavefunctions/DVRWavefunction.py#L88?message=Update%20Docs)]
</div>
**LLM Docstring**

A (lazily built and cached) continuous interpolant of the wavefunction's grid values, used by `evaluate` to evaluate the wavefunction off-grid.
  - `:returns`: `Interpolator`
    > the cached (or newly built) interpolator


<a id="Psience.DVR.Wavefunctions.DVRWavefunction.evaluate" class="docs-object-method">&nbsp;</a> 
```python
evaluate(self, points): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Wavefunctions/DVRWavefunction.py#L101)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Wavefunctions/DVRWavefunction.py#L101?message=Update%20Docs)]
</div>
Evaluates the functions at the given points
  - `:returns`: `_`
    >


<a id="Psience.DVR.Wavefunctions.DVRWavefunction.marginalize_out" class="docs-object-method">&nbsp;</a> 
```python
marginalize_out(self, dofs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Wavefunctions/DVRWavefunction.py#L110)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Wavefunctions/DVRWavefunction.py#L110?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/DVR/Wavefunctions/DVRWavefunction.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/DVR/Wavefunctions/DVRWavefunction.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/DVR/Wavefunctions/DVRWavefunction.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/DVR/Wavefunctions/DVRWavefunction.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Wavefunctions.py#L15?message=Update%20Docs)   
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