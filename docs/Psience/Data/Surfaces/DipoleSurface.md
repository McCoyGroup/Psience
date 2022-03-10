## <a id="Psience.Data.Surfaces.DipoleSurface">DipoleSurface</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Data/Surfaces.py#L15)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Data/Surfaces.py#L15?message=Update%20Docs)]
</div>

Provides a unified interface to working with dipole surfaces.
Currently basically no fancier than a regular surface (although with convenient loading functions), but dipole-specific
stuff could come

<a id="Psience.Data.Surfaces.DipoleSurface.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, mu_x, mu_y, mu_z): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Data/Surfaces.py#L21)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Data/Surfaces.py#L21?message=Update%20Docs)]
</div>


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
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Data/Surfaces.py#L41)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Data/Surfaces.py#L41?message=Update%20Docs)]
</div>

<a id="Psience.Data.Surfaces.DipoleSurface.from_log_file" class="docs-object-method">&nbsp;</a> 
```python
from_log_file(log_file, coord_transf, keys=('StandardCartesianCoordinates', 'DipoleMoments'), tol=0.001, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Data/Surfaces.py#L51)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Data/Surfaces.py#L51?message=Update%20Docs)]
</div>

Loads dipoles from a Gaussian log file and builds a dipole surface by interpolating.
        Obviously this only really works if we have a subset of "scan" coordinates, so at this stage the user is obligated
        to furnish a function that'll take a set of Cartesian coordinates and convert them to "scan" coordinates.
        Coordinerds can be helpful with this, as it provides a convenient syntax for Cartesian <-> ZMatrix conversions
- `log_file`: `str`
    >a Gaussian log file to pull from
- `:returns`: `_`
    >No description...

<a id="Psience.Data.Surfaces.DipoleSurface.get_fchk_values" class="docs-object-method">&nbsp;</a> 
```python
get_fchk_values(fchk_file): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Data/Surfaces.py#L106)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Data/Surfaces.py#L106?message=Update%20Docs)]
</div>

<a id="Psience.Data.Surfaces.DipoleSurface.from_fchk_file" class="docs-object-method">&nbsp;</a> 
```python
from_fchk_file(fchk_file, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Data/Surfaces.py#L118)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Data/Surfaces.py#L118?message=Update%20Docs)]
</div>

Loads dipoles from a Gaussian formatted checkpoint file and builds a dipole surface via a linear approximation
- `fchk_file`: `Any`
    >a Gaussian fchk file to pull from
- `log_file`: `str`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.Data.Surfaces.DipoleSurface.__call__" class="docs-object-method">&nbsp;</a> 
```python
__call__(self, gridpoints, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Data/Surfaces.py#L146)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Data/Surfaces.py#L146?message=Update%20Docs)]
</div>

Explicitly overrides the Surface-level evaluation because we know the Taylor surface needs us to flatten our gridpoints
- `gridpoints`: `Any`
    >No description...
- `opts`: `Any`
    >No description...
- `:returns`: `_`
    >No description...



___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/ci/docs/Psience/Data/Surfaces/DipoleSurface.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/ci/docs/Psience/Data/Surfaces/DipoleSurface.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/ci/docs/Psience/Data/Surfaces/DipoleSurface.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/ci/docs/Psience/Data/Surfaces/DipoleSurface.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Data/Surfaces.py#L15?message=Update%20Docs)