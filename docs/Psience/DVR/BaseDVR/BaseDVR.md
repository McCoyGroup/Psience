## <a id="Psience.DVR.BaseDVR.BaseDVR">BaseDVR</a>
Provides the abstract interface for creating a
convenient runnable DVR that can be cleanly subclassed to provide
extensions

### Properties and Methods
<a id="Psience.DVR.BaseDVR.BaseDVR.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, domain=None, divs=None, potential_function=None, **base_opts): 
```

- `base_opts`: `Any`
    >base opts to use when running

<a id="Psience.DVR.BaseDVR.BaseDVR.__repr__" class="docs-object-method">&nbsp;</a>
```python
__repr__(self): 
```

<a id="Psience.DVR.BaseDVR.BaseDVR.get_grid" class="docs-object-method">&nbsp;</a>
```python
get_grid(self, domain=None, divs=None, **kwargs): 
```

<a id="Psience.DVR.BaseDVR.BaseDVR.grid" class="docs-object-method">&nbsp;</a>
```python
grid(self, domain=None, divs=None, **kwargs): 
```

<a id="Psience.DVR.BaseDVR.BaseDVR.get_kinetic_energy" class="docs-object-method">&nbsp;</a>
```python
get_kinetic_energy(self, grid=None, mass=None, hb=1, **kwargs): 
```

<a id="Psience.DVR.BaseDVR.BaseDVR.kinetic_energy" class="docs-object-method">&nbsp;</a>
```python
kinetic_energy(self, grid=None, mass=None, hb=1, g=None, g_deriv=None, **kwargs): 
```

<a id="Psience.DVR.BaseDVR.BaseDVR.real_momentum" class="docs-object-method">&nbsp;</a>
```python
real_momentum(self, grid=None, mass=None, hb=1, **kwargs): 
```

<a id="Psience.DVR.BaseDVR.BaseDVR.potential_energy" class="docs-object-method">&nbsp;</a>
```python
potential_energy(self, grid=None, potential_function=None, potential_values=None, potential_grid=None, **pars): 
```
Calculates the potential energy at the grid points based
        on dispatching on the input form of the potential
- `grid`: `Any`
    >the grid of points built earlier in the DVR
- `potential_function`: `Any`
    >a function to evaluate the potential energy at the points
- `potential_values`: `Any`
    >the values of the potential at the DVR points
- `potential_grid`: `Any`
    >a grid of points and values to be interpolated
- `pars`: `Any`
    >ignored keyword arguments
- `:returns`: `_`
    >No description...

<a id="Psience.DVR.BaseDVR.BaseDVR.hamiltonian" class="docs-object-method">&nbsp;</a>
```python
hamiltonian(self, kinetic_energy=None, potential_energy=None, potential_threshold=None, **pars): 
```
Calculates the total Hamiltonian from the kinetic and potential matrices
- `kinetic_energy`: `Any`
    >No description...
- `potential_energy`: `np.ndarray | sp.spmatrix`
    >No description...
- `potential_threshold`: `Any`
    >No description...
- `pars`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.DVR.BaseDVR.BaseDVR.wavefunctions" class="docs-object-method">&nbsp;</a>
```python
wavefunctions(self, hamiltonian=None, num_wfns=25, nodeless_ground_state=False, diag_mode=None, **pars): 
```
Calculates the wavefunctions for the given Hamiltonian.
        Doesn't support any kind of pruning based on potential values although that might be a good feature
        to support explicitly in the future
- `hamiltonian`: `Any`
    >No description...
- `num_wfns`: `Any`
    >No description...
- `nodeless_ground_state`: `Any`
    >No description...
- `diag_mode`: `Any`
    >No description...
- `pars`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.DVR.BaseDVR.BaseDVR.run" class="docs-object-method">&nbsp;</a>
```python
run(self, result='wavefunctions', **opts): 
```

- `:returns`: `DVRResults`
    >No description...

### Examples




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/edit/ci/examples/ci/docs/Psience/DVR/BaseDVR/BaseDVR.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/examples/ci/docs/Psience/DVR/BaseDVR/BaseDVR.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/edit/ci/docs/ci/docs/Psience/DVR/BaseDVR/BaseDVR.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/docs/templates/ci/docs/Psience/DVR/BaseDVR/BaseDVR.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/DVR/BaseDVR.py?message=Update%20Docs)