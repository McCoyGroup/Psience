## <a id="Psience.DGB.Wavefunctions.DGBWavefunctions">DGBWavefunctions</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Wavefunctions.py#L289)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Wavefunctions.py#L289?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
wavefunction_class: DGBWavefunction
```
<a id="Psience.DGB.Wavefunctions.DGBWavefunctions.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, energies=None, wavefunctions=None, hamiltonian=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Wavefunctions.py#L291)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Wavefunctions.py#L291?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Wavefunctions.DGBWavefunctions.as_cartesian_wavefunction" class="docs-object-method">&nbsp;</a> 
```python
as_cartesian_wavefunction(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Wavefunctions/DGBWavefunctions.py#L295)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Wavefunctions/DGBWavefunctions.py#L295?message=Update%20Docs)]
</div>
Projects the wavefunction back to Cartesians
  - `:returns`: `_`
    >


<a id="Psience.DGB.Wavefunctions.DGBWavefunctions.hamiltonian" class="docs-object-method">&nbsp;</a> 
```python
@property
hamiltonian(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Wavefunctions/DGBWavefunctions.py#L320)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Wavefunctions/DGBWavefunctions.py#L320?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Wavefunctions.DGBWavefunctions.operator_representation" class="docs-object-method">&nbsp;</a> 
```python
operator_representation(self, op, embed=True, expansion_degree=None, quadrature_degree=None, expansion_type=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Wavefunctions/DGBWavefunctions.py#L334)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Wavefunctions/DGBWavefunctions.py#L334?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Wavefunctions.DGBWavefunctions.expectation" class="docs-object-method">&nbsp;</a> 
```python
expectation(self, op, expansion_degree=None, quadrature_degree=None, expansion_type=None, embed=True, other=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Wavefunctions/DGBWavefunctions.py#L349)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Wavefunctions/DGBWavefunctions.py#L349?message=Update%20Docs)]
</div>
Computes the expectation value of operator op over the wavefunction other and self
  - `other`: `Wavefunction | np.ndarray`
    > 
  - `op`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.DGB.Wavefunctions.DGBWavefunctions.localize" class="docs-object-method">&nbsp;</a> 
```python
localize(self, criterion, which=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Wavefunctions/DGBWavefunctions.py#L395)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Wavefunctions/DGBWavefunctions.py#L395?message=Update%20Docs)]
</div>
Find a transformation that maximally localizes the wavefunctions in the Boys' sense
by minimizing <r^2> - <r>^2 over unitary transformations
  - `criterion`: `Any`
    > 
  - `which`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.DGB.Wavefunctions.DGBWavefunctions.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Wavefunctions/DGBWavefunctions.py#L419)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Wavefunctions/DGBWavefunctions.py#L419?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/DGB/Wavefunctions/DGBWavefunctions.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/DGB/Wavefunctions/DGBWavefunctions.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/DGB/Wavefunctions/DGBWavefunctions.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/DGB/Wavefunctions/DGBWavefunctions.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Wavefunctions.py#L289?message=Update%20Docs)   
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