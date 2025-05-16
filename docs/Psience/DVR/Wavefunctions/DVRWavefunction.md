## <a id="Psience.DVR.Wavefunctions.DVRWavefunction">DVRWavefunction</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/Wavefunctions.py#L15)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/Wavefunctions.py#L15?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/Wavefunctions/DVRWavefunction.py#L17)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/Wavefunctions/DVRWavefunction.py#L17?message=Update%20Docs)]
</div>


<a id="Psience.DVR.Wavefunctions.DVRWavefunction.get_dimension" class="docs-object-method">&nbsp;</a> 
```python
get_dimension(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/Wavefunctions/DVRWavefunction.py#L23)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/Wavefunctions/DVRWavefunction.py#L23?message=Update%20Docs)]
</div>


<a id="Psience.DVR.Wavefunctions.DVRWavefunction.plot" class="docs-object-method">&nbsp;</a> 
```python
plot(self, figure=None, grid=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/Wavefunctions/DVRWavefunction.py#L26)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/Wavefunctions/DVRWavefunction.py#L26?message=Update%20Docs)]
</div>


<a id="Psience.DVR.Wavefunctions.DVRWavefunction.expectation" class="docs-object-method">&nbsp;</a> 
```python
expectation(self, op, other=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/Wavefunctions/DVRWavefunction.py#L31)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/Wavefunctions/DVRWavefunction.py#L31?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/Wavefunctions/DVRWavefunction.py#L46)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/Wavefunctions/DVRWavefunction.py#L46?message=Update%20Docs)]
</div>


<a id="Psience.DVR.Wavefunctions.DVRWavefunction.evaluate" class="docs-object-method">&nbsp;</a> 
```python
evaluate(self, points): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/Wavefunctions/DVRWavefunction.py#L51)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/Wavefunctions/DVRWavefunction.py#L51?message=Update%20Docs)]
</div>
Evaluates the functions at the given points
  - `:returns`: `_`
    >


<a id="Psience.DVR.Wavefunctions.DVRWavefunction.marginalize_out" class="docs-object-method">&nbsp;</a> 
```python
marginalize_out(self, dofs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/Wavefunctions/DVRWavefunction.py#L60)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/Wavefunctions/DVRWavefunction.py#L60?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/Wavefunctions.py#L15?message=Update%20Docs)   
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