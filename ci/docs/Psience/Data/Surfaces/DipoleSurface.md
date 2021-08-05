## <a id="Psience.Data.Surfaces.DipoleSurface">DipoleSurface</a>
Provides a unified interface to working with dipole surfaces.
Currently basically no fancier than a regular surface (although with convenient loading functions), but dipole-specific
stuff could come

### Properties and Methods
```python
from_log_file: method
from_fchk_file: method
```
<a id="Psience.Data.Surfaces.DipoleSurface.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, mu_x, mu_y, mu_z): 
```

- `mu_x`: `Surface`
    >X-component of dipole moment
- `mu_y`: `Surface`
    >Y-component of dipole moment
- `mu_z`: `Surface`
    >Z-component of dipole moment

<a id="Psience.Data.Surfaces.DipoleSurface.get_log_values" class="docs-object-method">&nbsp;</a>
```python
get_log_values(log_file, keys=('StandardCartesianCoordinates', 'DipoleMoments')): 
```

<a id="Psience.Data.Surfaces.DipoleSurface.get_fchk_values" class="docs-object-method">&nbsp;</a>
```python
get_fchk_values(fchk_file): 
```

<a id="Psience.Data.Surfaces.DipoleSurface.__call__" class="docs-object-method">&nbsp;</a>
```python
__call__(self, gridpoints, **opts): 
```
Explicitly overrides the Surface-level evaluation because we know the Taylor surface needs us to flatten our gridpoints
- `gridpoints`: `Any`
    >No description...
- `opts`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

### Examples


