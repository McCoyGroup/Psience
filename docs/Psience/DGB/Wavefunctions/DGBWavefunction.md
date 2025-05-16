## <a id="Psience.DGB.Wavefunctions.DGBWavefunction">DGBWavefunction</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Wavefunctions.py#L18)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Wavefunctions.py#L18?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.DGB.Wavefunctions.DGBWavefunction.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, energy, data, gaussians: Psience.DGB.Gaussians.DGBGaussians = None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Wavefunctions.py#L19)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Wavefunctions.py#L19?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Wavefunctions.DGBWavefunction.get_dimension" class="docs-object-method">&nbsp;</a> 
```python
get_dimension(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Wavefunctions/DGBWavefunction.py#L25)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Wavefunctions/DGBWavefunction.py#L25?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Wavefunctions.DGBWavefunction.plot" class="docs-object-method">&nbsp;</a> 
```python
plot(self, figure=None, domain=None, plot_centers=False, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Wavefunctions/DGBWavefunction.py#L28)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Wavefunctions/DGBWavefunction.py#L28?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Wavefunctions.DGBWavefunction.to_cartesian_wavefunction" class="docs-object-method">&nbsp;</a> 
```python
to_cartesian_wavefunction(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Wavefunctions/DGBWavefunction.py#L64)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Wavefunctions/DGBWavefunction.py#L64?message=Update%20Docs)]
</div>
Projects the wavefunction back to Cartesians
  - `:returns`: `_`
    >


<a id="Psience.DGB.Wavefunctions.DGBWavefunction.plot_cartesians" class="docs-object-method">&nbsp;</a> 
```python
plot_cartesians(self, xyz_sel=None, *, atom_sel=None, figure=None, plot_centers=False, atom_styles=None, **plot_styles): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Wavefunctions/DGBWavefunction.py#L80)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Wavefunctions/DGBWavefunction.py#L80?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Wavefunctions.DGBWavefunction.evaluate" class="docs-object-method">&nbsp;</a> 
```python
evaluate(self, points): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Wavefunctions/DGBWavefunction.py#L146)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Wavefunctions/DGBWavefunction.py#L146?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Wavefunctions.DGBWavefunction.marginalize_out" class="docs-object-method">&nbsp;</a> 
```python
marginalize_out(self, dofs, rescale=True) -> 'DGBWavefunction': 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DGB/Wavefunctions/DGBWavefunction.py#L203)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Wavefunctions/DGBWavefunction.py#L203?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/DGB/Wavefunctions/DGBWavefunction.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/DGB/Wavefunctions/DGBWavefunction.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/DGB/Wavefunctions/DGBWavefunction.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/DGB/Wavefunctions/DGBWavefunction.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DGB/Wavefunctions.py#L18?message=Update%20Docs)   
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