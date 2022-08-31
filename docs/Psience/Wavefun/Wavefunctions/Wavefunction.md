## <a id="Psience.Wavefun.Wavefunctions.Wavefunction">Wavefunction</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Wavefun/Wavefunctions.py#L17)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions.py#L17?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Wavefun/Wavefunctions/Wavefunction.py#L19)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions/Wavefunction.py#L19?message=Update%20Docs)]
</div>


<a id="Psience.Wavefun.Wavefunctions.Wavefunction.plot" class="docs-object-method">&nbsp;</a> 
```python
plot(self, figure=None, index=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Wavefun/Wavefunctions/Wavefunction.py#L25)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions/Wavefunction.py#L25?message=Update%20Docs)]
</div>
Uses McUtils to plot the wavefunction on the passed figure (makes a new one if none)
  - `figure`: `Graphics | Graphics3D`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Wavefun.Wavefunctions.Wavefunction.expectation" class="docs-object-method">&nbsp;</a> 
```python
expectation(self, op, other=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Wavefun/Wavefunctions/Wavefunction.py#L35)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions/Wavefunction.py#L35?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Wavefun/Wavefunctions/Wavefunction.py#L47)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions/Wavefunction.py#L47?message=Update%20Docs)]
</div>


<a id="Psience.Wavefun.Wavefunctions.Wavefunction.probability_density" class="docs-object-method">&nbsp;</a> 
```python
probability_density(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Wavefun/Wavefunctions/Wavefunction.py#L49)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions/Wavefunction.py#L49?message=Update%20Docs)]
</div>
Computes the probability density of the current wavefunction
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Wavefun/Wavefunctions/Wavefunction.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Wavefun/Wavefunctions/Wavefunction.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Wavefun/Wavefunctions/Wavefunction.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Wavefun/Wavefunctions/Wavefunction.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Wavefun/Wavefunctions.py#L17?message=Update%20Docs)   
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