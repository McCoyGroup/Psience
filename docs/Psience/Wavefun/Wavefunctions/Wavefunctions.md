## <a id="Psience.Wavefun.Wavefunctions.Wavefunctions">Wavefunctions</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Wavefun/Wavefunctions.py#L229)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions.py#L229?message=Update%20Docs)]
</div>

An object representing a set of wavefunctions.
Provides concrete, but potentially inefficient methods for doing all the wavefunction ops.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
wavefunction_class: Wavefunction
```
<a id="Psience.Wavefun.Wavefunctions.Wavefunctions.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, energies=None, wavefunctions=None, indices=None, wavefunction_class=None, dipole_function=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Wavefun/Wavefunctions/Wavefunctions.py#L236)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions/Wavefunctions.py#L236?message=Update%20Docs)]
</div>


<a id="Psience.Wavefun.Wavefunctions.Wavefunctions.get_wavefunctions" class="docs-object-method">&nbsp;</a> 
```python
get_wavefunctions(self, which): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Wavefun/Wavefunctions/Wavefunctions.py#L248)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions/Wavefunctions.py#L248?message=Update%20Docs)]
</div>


<a id="Psience.Wavefun.Wavefunctions.Wavefunctions.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Wavefun/Wavefunctions/Wavefunctions.py#L269)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions/Wavefunctions.py#L269?message=Update%20Docs)]
</div>
Returns a single Wavefunction object


<a id="Psience.Wavefun.Wavefunctions.Wavefunctions.__len__" class="docs-object-method">&nbsp;</a> 
```python
__len__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Wavefun/Wavefunctions/Wavefunctions.py#L273)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions/Wavefunctions.py#L273?message=Update%20Docs)]
</div>


<a id="Psience.Wavefun.Wavefunctions.Wavefunctions.__iter__" class="docs-object-method">&nbsp;</a> 
```python
__iter__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Wavefun/Wavefunctions/Wavefunctions.py#L275)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions/Wavefunctions.py#L275?message=Update%20Docs)]
</div>


<a id="Psience.Wavefun.Wavefunctions.Wavefunctions.frequencies" class="docs-object-method">&nbsp;</a> 
```python
frequencies(self, start_at=0): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Wavefun/Wavefunctions/Wavefunctions.py#L279)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions/Wavefunctions.py#L279?message=Update%20Docs)]
</div>


<a id="Psience.Wavefun.Wavefunctions.Wavefunctions.get_spectrum" class="docs-object-method">&nbsp;</a> 
```python
get_spectrum(self, dipole_function=None, *, start_at=0, **options): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Wavefun/Wavefunctions/Wavefunctions.py#L282)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions/Wavefunctions.py#L282?message=Update%20Docs)]
</div>


<a id="Psience.Wavefun.Wavefunctions.Wavefunctions.plot" class="docs-object-method">&nbsp;</a> 
```python
plot(self, figure=None, graphics_class=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Wavefun/Wavefunctions/Wavefunctions.py#L305)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions/Wavefunctions.py#L305?message=Update%20Docs)]
</div>
Plots all of the wavefunctions on one set of axes
  - `opts`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Wavefun.Wavefunctions.Wavefunctions.expectation" class="docs-object-method">&nbsp;</a> 
```python
expectation(self, op, other=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Wavefun/Wavefunctions/Wavefunctions.py#L340)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions/Wavefunctions.py#L340?message=Update%20Docs)]
</div>
Computes the expectation value of operator op over the wavefunction other and self
  - `other`: `Wavefunctions`
    > 
  - `op`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Wavefun.Wavefunctions.Wavefunctions.overlap" class="docs-object-method">&nbsp;</a> 
```python
overlap(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Wavefun/Wavefunctions/Wavefunctions.py#L361)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions/Wavefunctions.py#L361?message=Update%20Docs)]
</div>


<a id="Psience.Wavefun.Wavefunctions.Wavefunctions.coordinate" class="docs-object-method">&nbsp;</a> 
```python
coordinate(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Wavefun/Wavefunctions/Wavefunctions.py#L364)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions/Wavefunctions.py#L364?message=Update%20Docs)]
</div>
Provides the coordinate operator in the wavefunction basis
  - `:returns`: `_`
    >


<a id="Psience.Wavefun.Wavefunctions.Wavefunctions.momentum" class="docs-object-method">&nbsp;</a> 
```python
momentum(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Wavefun/Wavefunctions/Wavefunctions.py#L372)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions/Wavefunctions.py#L372?message=Update%20Docs)]
</div>
Provides the real part of the representation of the momentum operator in the wavefunction basis
  - `:returns`: `_`
    >


<a id="Psience.Wavefun.Wavefunctions.Wavefunctions.laplacian" class="docs-object-method">&nbsp;</a> 
```python
laplacian(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Wavefun/Wavefunctions/Wavefunctions.py#L380)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions/Wavefunctions.py#L380?message=Update%20Docs)]
</div>
Provides the representation of the laplacian in the wavefunction basis
  - `:returns`: `_`
    >


<a id="Psience.Wavefun.Wavefunctions.Wavefunctions.kinetic_energy" class="docs-object-method">&nbsp;</a> 
```python
kinetic_energy(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Wavefun/Wavefunctions/Wavefunctions.py#L388)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions/Wavefunctions.py#L388?message=Update%20Docs)]
</div>
Provides the representation of the KE in the wavefunction basis
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Wavefun/Wavefunctions/Wavefunctions.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Wavefun/Wavefunctions/Wavefunctions.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Wavefun/Wavefunctions/Wavefunctions.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Wavefun/Wavefunctions/Wavefunctions.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions.py#L229?message=Update%20Docs)   
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