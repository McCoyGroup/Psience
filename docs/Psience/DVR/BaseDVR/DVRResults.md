## <a id="Psience.DVR.BaseDVR.DVRResults">DVRResults</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/BaseDVR.py#L380)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/BaseDVR.py#L380?message=Update%20Docs)]
</div>

A subclass that can wrap all of the DVR run parameters and results into a clean interface for reuse and extension







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.DVR.BaseDVR.DVRResults.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, grid=None, kinetic_energy=None, potential_energy=None, hamiltonian=None, wavefunctions=None, parent=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/BaseDVR.py#L384)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/BaseDVR.py#L384?message=Update%20Docs)]
</div>


<a id="Psience.DVR.BaseDVR.DVRResults.dimension" class="docs-object-method">&nbsp;</a> 
```python
@property
dimension(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/BaseDVR/DVRResults.py#L403)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/BaseDVR/DVRResults.py#L403?message=Update%20Docs)]
</div>


<a id="Psience.DVR.BaseDVR.DVRResults.plot_potential" class="docs-object-method">&nbsp;</a> 
```python
plot_potential(self, plot_class=None, figure=None, plot_units=None, energy_threshold=None, zero_shift=False, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/BaseDVR/DVRResults.py#L410)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/BaseDVR/DVRResults.py#L410?message=Update%20Docs)]
</div>
Simple plotting function for the potential.
Should be updated to deal with higher dimensional cases
  - `plot_class`: `McUtils.Plots.Graphics`
    > the graphics class to use for the plot
  - `opts`: `Any`
    > plot styling options
  - `:returns`: `McUtils.Plots.Graphics`
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/DVR/BaseDVR/DVRResults.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/DVR/BaseDVR/DVRResults.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/DVR/BaseDVR/DVRResults.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/DVR/BaseDVR/DVRResults.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/BaseDVR.py#L380?message=Update%20Docs)   
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