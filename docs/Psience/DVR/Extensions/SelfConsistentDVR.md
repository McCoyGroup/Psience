## <a id="Psience.DVR.Extensions.SelfConsistentDVR">SelfConsistentDVR</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Extensions.py#L54)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Extensions.py#L54?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.DVR.Extensions.SelfConsistentDVR.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, base_dvr: 'DirectProductDVR', **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Extensions.py#L55)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Extensions.py#L55?message=Update%20Docs)]
</div>
**LLM Docstring**

Set up a self-consistent-field (SCF) treatment of a multi-dimensional `DirectProductDVR`, precomputing the grid and potential-energy values and building the per-dimension SCF wavefunction generators before delegating to `GridSCF.__init__`.
  - `base_dvr`: `DirectProductDVR`
    > the multi-dimensional direct-product DVR to run the SCF procedure on
  - `opts`: `dict`
    > extra options, filtered and forwarded to the base `GridSCF.__init__`
  - `:returns`: `None`
    > None


<a id="Psience.DVR.Extensions.SelfConsistentDVR.create_grid_vals" class="docs-object-method">&nbsp;</a> 
```python
create_grid_vals(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Extensions/SelfConsistentDVR.py#L80)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Extensions/SelfConsistentDVR.py#L80?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the full grid and potential-energy values for the underlying multi-dimensional DVR, used to initialize the SCF procedure.
  - `:returns`: `tuple[np.ndarray, np.ndarray]`
    > `(grid, pe)` -- the DVR grid and the potential-energy values reshaped to match the grid's spatial shape


<a id="Psience.DVR.Extensions.SelfConsistentDVR.create_solvers" class="docs-object-method">&nbsp;</a> 
```python
create_solvers(self, grid, pe): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Extensions/SelfConsistentDVR.py#L98)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Extensions/SelfConsistentDVR.py#L98?message=Update%20Docs)]
</div>
**LLM Docstring**

Build the per-dimension SCF wavefunction generators, first rebinding each 1D DVR's `g`/`g_deriv` kinetic-coupling functions (if callable) so they can be evaluated at the fixed "other-dimension" grid point (from the initial guess) while varying only their own coordinate.
  - `grid`: `np.ndarray`
    > the full multi-dimensional grid
  - `pe`: `np.ndarray`
    > the potential-energy values on the grid
  - `:returns`: `list[SCFWavefunctionGenerator]`
    > the list of per-dimension `SCFWavefunctionGenerator` objects


<a id="Psience.DVR.Extensions.SelfConsistentDVR.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Extensions/SelfConsistentDVR.py#L185)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Extensions/SelfConsistentDVR.py#L185?message=Update%20Docs)]
</div>
**LLM Docstring**

Debug string representation showing the class name and the wrapped base DVR.
  - `:returns`: `str`
    > string of the form `ClassName(base_dvr)`
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/DVR/Extensions/SelfConsistentDVR.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/DVR/Extensions/SelfConsistentDVR.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/DVR/Extensions/SelfConsistentDVR.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/DVR/Extensions/SelfConsistentDVR.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Extensions.py#L54?message=Update%20Docs)   
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