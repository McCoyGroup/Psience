## <a id="Psience.Data.Surfaces.PotentialSurface">PotentialSurface</a>
A potential surface structure to go along with the DipoleSurface.
Provides convenient access to dipole data + a unified interface to things like energy minimization

### Properties and Methods
<a id="Psience.Data.Surfaces.PotentialSurface.get_log_values" class="docs-object-method">&nbsp;</a>
```python
get_log_values(log_file, keys=('StandardCartesianCoordinates', 'ScanEnergies')): 
```

<a id="Psience.Data.Surfaces.PotentialSurface.from_log_file" class="docs-object-method">&nbsp;</a>
```python
from_log_file(log_file, coord_transf, keys=('StandardCartesianCoordinates', 'ScanEnergies'), tol=0.001, **opts): 
```
Loads dipoles from a Gaussian log file and builds a potential surface by interpolating.
        Obviously this only really works if we have a subset of "scan" coordinates, so at this stage the user is obligated
        to furnish a function that'll take a set of Cartesian coordinates and convert them to "scan" coordinates.
        Coordinerds can be helpful with this, as it provides a convenient syntax for Cartesian <-> ZMatrix conversions.
- `log_file`: `str`
    >a Gaussian log file to pull from
- `:returns`: `_`
    >No description...

<a id="Psience.Data.Surfaces.PotentialSurface.get_fchk_values" class="docs-object-method">&nbsp;</a>
```python
get_fchk_values(fchk_file): 
```

<a id="Psience.Data.Surfaces.PotentialSurface.from_fchk_file" class="docs-object-method">&nbsp;</a>
```python
from_fchk_file(fchk_file, **opts): 
```
Loads potential from a Gaussian formatted checkpoint file and builds a potential surface via a quartic approximation
- `fchk_file`: `Any`
    >a Gaussian fchk file to pull from
- `log_file`: `str`
    >No description...
- `:returns`: `_`
    >No description...

### Examples




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/edit/ci/examples/ci/docs/Psience/Data/Surfaces/PotentialSurface.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/examples/ci/docs/Psience/Data/Surfaces/PotentialSurface.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/edit/ci/docs/ci/docs/Psience/Data/Surfaces/PotentialSurface.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/docs/templates/ci/docs/Psience/Data/Surfaces/PotentialSurface.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Data/Surfaces.py?message=Update%20Docs)