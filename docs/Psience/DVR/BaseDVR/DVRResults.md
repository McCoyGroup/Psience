## <a id="Psience.DVR.BaseDVR.DVRResults">DVRResults</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/BaseDVR.py#L380)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/BaseDVR.py#L380?message=Update%20Docs)]
</div>

A subclass that can wrap all of the DVR run parameters and results into a clean interface for reuse and extension

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

<a id="Psience.DVR.BaseDVR.DVRResults.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, grid=None, kinetic_energy=None, potential_energy=None, hamiltonian=None, wavefunctions=None, parent=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/BaseDVR.py#L384)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/BaseDVR.py#L384?message=Update%20Docs)]
</div>

<a id="Psience.DVR.BaseDVR.DVRResults.dimension" class="docs-object-method">&nbsp;</a> 
```python
@property
dimension(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/BaseDVR.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/BaseDVR.py#L?message=Update%20Docs)]
</div>

<a id="Psience.DVR.BaseDVR.DVRResults.plot_potential" class="docs-object-method">&nbsp;</a> 
```python
plot_potential(self, plot_class=None, figure=None, plot_units=None, energy_threshold=None, zero_shift=False, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/BaseDVR.py#L410)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/BaseDVR.py#L410?message=Update%20Docs)]
</div>

Simple plotting function for the potential.
        Should be updated to deal with higher dimensional cases
- `plot_class`: `McUtils.Plots.Graphics`
    >the graphics class to use for the plot
- `opts`: `Any`
    >plot styling options
- `:returns`: `McUtils.Plots.Graphics`
    >No description...

 </div>
</div>




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/DVR/BaseDVR/DVRResults.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/DVR/BaseDVR/DVRResults.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/DVR/BaseDVR/DVRResults.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/DVR/BaseDVR/DVRResults.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/BaseDVR.py#L380?message=Update%20Docs)