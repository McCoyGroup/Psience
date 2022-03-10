## <a id="Psience.Molecools.Vibrations.MolecularNormalModes">MolecularNormalModes</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Molecools/Vibrations.py#L241)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Molecools/Vibrations.py#L241?message=Update%20Docs)]
</div>

A Coordinerds CoordinateSystem object that manages all of the data needed to
work with normal mode coordinates + some convenience functions for generating and whatnot

```python
name: str
```
<a id="Psience.Molecools.Vibrations.MolecularNormalModes.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, molecule, coeffs, name=None, freqs=None, internal=False, origin=None, basis=None, inverse=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Molecools/Vibrations.py#L247)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Molecools/Vibrations.py#L247?message=Update%20Docs)]
</div>


- `molecule`: `AbstractMolecule`
    >No description...
- `coeffs`: `Any`
    >No description...
- `name`: `Any`
    >No description...
- `freqs`: `Any`
    >No description...
- `internal`: `Any`
    >No description...
- `origin`: `Any`
    >No description...
- `basis`: `Any`
    >No description...
- `inverse`: `Any`
    >No description...

<a id="Psience.Molecools.Vibrations.MolecularNormalModes.molecule" class="docs-object-method">&nbsp;</a> 
```python
@property
molecule(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Molecools/Vibrations.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Molecools/Vibrations.py#L?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Vibrations.MolecularNormalModes.to_internals" class="docs-object-method">&nbsp;</a> 
```python
to_internals(self, intcrds=None, dYdR=None, dRdY=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Molecools/Vibrations.py#L312)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Molecools/Vibrations.py#L312?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Vibrations.MolecularNormalModes.origin" class="docs-object-method">&nbsp;</a> 
```python
@property
origin(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Molecools/Vibrations.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Molecools/Vibrations.py#L?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Vibrations.MolecularNormalModes.embed" class="docs-object-method">&nbsp;</a> 
```python
embed(self, frame): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Molecools/Vibrations.py#L352)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Molecools/Vibrations.py#L352?message=Update%20Docs)]
</div>


- `frame`: `MolecularTransformation`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.Molecools.Vibrations.MolecularNormalModes.insert" class="docs-object-method">&nbsp;</a> 
```python
insert(self, val, where): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Molecools/Vibrations.py#L416)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Molecools/Vibrations.py#L416?message=Update%20Docs)]
</div>

Inserts values into the appropriate positions in the mode matrix
- `val`: `Any`
    >No description...
- `where`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.Molecools.Vibrations.MolecularNormalModes.rescale" class="docs-object-method">&nbsp;</a> 
```python
rescale(self, scaling_factors): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Molecools/Vibrations.py#L448)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Molecools/Vibrations.py#L448?message=Update%20Docs)]
</div>

Rescales each mode in the expansion matrix
        by the passed `scaling_factors`
- `scaling_factors`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.Molecools.Vibrations.MolecularNormalModes.from_force_constants" class="docs-object-method">&nbsp;</a> 
```python
from_force_constants(molecule, fcs, atoms=None, masses=None, mass_units='AtomicMassUnits', inverse_mass_matrix=False, remove_transrot=True, normalize=False, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Molecools/Vibrations.py#L469)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Molecools/Vibrations.py#L469?message=Update%20Docs)]
</div>

Generates normal modes from the specified force constants
- `molecule`: `AbstractMolecule`
    >No description...
- `fcs`: `np.ndarray`
    >force constants array
- `atoms`: `Iterable[str]`
    >atom list
- `masses`: `Iterable[float]`
    >mass list
- `mass_units`: `str`
    >units for the masses...not clear if this is useful or a distraction
- `inverse_mass_matrix`: `bool`
    >whether or not we have G or G^-1 (default: `False`)
- `remove_transrot`: `bool`
    >whether or not to remove the translations and rotations (default: `True`)
- `normalize`: `bool`
    >whether or not to normalize the modes (default: `True`)
- `opts`: `Any`
    >No description...
- `:returns`: `MolecularNormalModes`
    >No description...

<a id="Psience.Molecools.Vibrations.MolecularNormalModes.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Molecools/Vibrations.py#L547)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Molecools/Vibrations.py#L547?message=Update%20Docs)]
</div>

Takes a slice of the modes
- `item`: `Any`
    >No description...
- `:returns`: `_`
    >No description...



___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/ci/docs/Psience/Molecools/Vibrations/MolecularNormalModes.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/ci/docs/Psience/Molecools/Vibrations/MolecularNormalModes.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/ci/docs/Psience/Molecools/Vibrations/MolecularNormalModes.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/ci/docs/Psience/Molecools/Vibrations/MolecularNormalModes.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Molecools/Vibrations.py#L241?message=Update%20Docs)