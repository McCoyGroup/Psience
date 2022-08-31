## <a id="Psience.DVR.Wavefunctions.DVRWavefunctions">DVRWavefunctions</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/Wavefunctions.py#L80)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/Wavefunctions.py#L80?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, energies=None, wavefunctions=None, grid=None, wavefunction_class=<class 'Psience.DVR.Wavefunctions.DVRWavefunction'>, results: Psience.DVR.BaseDVR.DVRResults = None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/Wavefunctions/DVRWavefunctions.py#L82)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/Wavefunctions/DVRWavefunctions.py#L82?message=Update%20Docs)]
</div>


<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/Wavefunctions/DVRWavefunctions.py#L86)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/Wavefunctions/DVRWavefunctions.py#L86?message=Update%20Docs)]
</div>


<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.__len__" class="docs-object-method">&nbsp;</a> 
```python
__len__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/Wavefunctions/DVRWavefunctions.py#L93)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/Wavefunctions/DVRWavefunctions.py#L93?message=Update%20Docs)]
</div>


<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.__iter__" class="docs-object-method">&nbsp;</a> 
```python
__iter__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/Wavefunctions/DVRWavefunctions.py#L95)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/Wavefunctions/DVRWavefunctions.py#L95?message=Update%20Docs)]
</div>


<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/Wavefunctions/DVRWavefunctions.py#L98)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/Wavefunctions/DVRWavefunctions.py#L98?message=Update%20Docs)]
</div>
Provides a single `DVRWavefunction` or slice of `DVRWavefunctions`
  - `item`: `Any`
    > 
  - `:returns`: `DVRWavefunction | DVRWavefunctions`
    >


<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.plot" class="docs-object-method">&nbsp;</a> 
```python
plot(self, figure=None, graphics_class=None, plot_style=None, scaling=1, shift=0, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/Wavefunctions/DVRWavefunctions.py#L123)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/Wavefunctions/DVRWavefunctions.py#L123?message=Update%20Docs)]
</div>
Plots the held wavefunctions
  - `figure`: `Any`
    > 
  - `graphics_class`: `Any`
    > 
  - `plot_style`: `Any`
    > 
  - `scaling`: `Any`
    > 
  - `shift`: `Any`
    > 
  - `opts`: `Any`
    > 
  - `:returns`: `Graphics`
    >


<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.expectation" class="docs-object-method">&nbsp;</a> 
```python
expectation(self, op, other=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/Wavefunctions/DVRWavefunctions.py#L159)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/Wavefunctions/DVRWavefunctions.py#L159?message=Update%20Docs)]
</div>
Computes the expectation value of operator op over the wavefunction other and self
  - `other`: `DVRWavefunctions | np.ndarray`
    > 
  - `op`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.transform_operator" class="docs-object-method">&nbsp;</a> 
```python
transform_operator(self, M): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/Wavefunctions/DVRWavefunctions.py#L187)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/Wavefunctions/DVRWavefunctions.py#L187?message=Update%20Docs)]
</div>


<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.probability_density" class="docs-object-method">&nbsp;</a> 
```python
probability_density(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/Wavefunctions/DVRWavefunctions.py#L192)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/Wavefunctions/DVRWavefunctions.py#L192?message=Update%20Docs)]
</div>
Computes the probability density of the set of wavefunctions
  - `:returns`: `_`
    >


<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.coordinate" class="docs-object-method">&nbsp;</a> 
```python
coordinate(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/Wavefunctions/DVRWavefunctions.py#L200)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/Wavefunctions/DVRWavefunctions.py#L200?message=Update%20Docs)]
</div>


<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.momentum" class="docs-object-method">&nbsp;</a> 
```python
momentum(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/Wavefunctions/DVRWavefunctions.py#L202)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/Wavefunctions/DVRWavefunctions.py#L202?message=Update%20Docs)]
</div>


<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.laplacian" class="docs-object-method">&nbsp;</a> 
```python
laplacian(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/Wavefunctions/DVRWavefunctions.py#L206)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/Wavefunctions/DVRWavefunctions.py#L206?message=Update%20Docs)]
</div>


<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.kinetic_energy" class="docs-object-method">&nbsp;</a> 
```python
kinetic_energy(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/Wavefunctions/DVRWavefunctions.py#L211)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/Wavefunctions/DVRWavefunctions.py#L211?message=Update%20Docs)]
</div>


<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.potential_energy" class="docs-object-method">&nbsp;</a> 
```python
potential_energy(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/Wavefunctions/DVRWavefunctions.py#L217)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/Wavefunctions/DVRWavefunctions.py#L217?message=Update%20Docs)]
</div>
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/DVR/Wavefunctions/DVRWavefunctions.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/DVR/Wavefunctions/DVRWavefunctions.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/DVR/Wavefunctions/DVRWavefunctions.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/DVR/Wavefunctions/DVRWavefunctions.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/Wavefunctions.py#L80?message=Update%20Docs)   
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