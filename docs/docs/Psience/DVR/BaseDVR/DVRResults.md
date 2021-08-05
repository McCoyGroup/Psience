## <a id="Psience.DVR.BaseDVR.DVRResults">DVRResults</a>
A subclass that can wrap all of the DVR run parameters and results into a clean interface for reuse and extension

### Properties and Methods
<a id="Psience.DVR.BaseDVR.DVRResults.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, grid=None, kinetic_energy=None, potential_energy=None, hamiltonian=None, wavefunctions=None, parent=None, **opts): 
```

<a id="Psience.DVR.BaseDVR.DVRResults.dimension" class="docs-object-method">&nbsp;</a>
```python
@property
dimension(self): 
```

<a id="Psience.DVR.BaseDVR.DVRResults.plot_potential" class="docs-object-method">&nbsp;</a>
```python
plot_potential(self, plot_class=None, figure=None, plot_units=None, energy_threshold=None, zero_shift=False, **opts): 
```
Simple plotting function for the potential.
        Should be updated to deal with higher dimensional cases
- `plot_class`: `McUtils.Plots.Graphics`
    >the graphics class to use for the plot
- `opts`: `Any`
    >plot styling options
- `:returns`: `McUtils.Plots.Graphics`
    >No description...

### Examples


