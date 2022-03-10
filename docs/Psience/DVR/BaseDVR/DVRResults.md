## <a id="Psience.DVR.BaseDVR.DVRResults">DVRResults</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/DVR/BaseDVR.py#L315)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/DVR/BaseDVR.py#L315?message=Update%20Docs)]
</div>

A subclass that can wrap all of the DVR run parameters and results into a clean interface for reuse and extension

<a id="Psience.DVR.BaseDVR.DVRResults.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, grid=None, kinetic_energy=None, potential_energy=None, hamiltonian=None, wavefunctions=None, parent=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/DVR/BaseDVR.py#L319)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/DVR/BaseDVR.py#L319?message=Update%20Docs)]
</div>

<a id="Psience.DVR.BaseDVR.DVRResults.dimension" class="docs-object-method">&nbsp;</a> 
```python
@property
dimension(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/DVR/BaseDVR.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/DVR/BaseDVR.py#L?message=Update%20Docs)]
</div>

<a id="Psience.DVR.BaseDVR.DVRResults.plot_potential" class="docs-object-method">&nbsp;</a> 
```python
plot_potential(self, plot_class=None, figure=None, plot_units=None, energy_threshold=None, zero_shift=False, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/DVR/BaseDVR.py#L345)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/DVR/BaseDVR.py#L345?message=Update%20Docs)]
</div>

Simple plotting function for the potential.
        Should be updated to deal with higher dimensional cases
- `plot_class`: `McUtils.Plots.Graphics`
    >the graphics class to use for the plot
- `opts`: `Any`
    >plot styling options
- `:returns`: `McUtils.Plots.Graphics`
    >No description...



___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/ci/docs/Psience/DVR/BaseDVR/DVRResults.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/ci/docs/Psience/DVR/BaseDVR/DVRResults.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/ci/docs/Psience/DVR/BaseDVR/DVRResults.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/ci/docs/Psience/DVR/BaseDVR/DVRResults.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/DVR/BaseDVR.py#L315?message=Update%20Docs)