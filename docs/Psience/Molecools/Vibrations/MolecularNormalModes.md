## <a id="Psience.Molecools.Vibrations.MolecularNormalModes">MolecularNormalModes</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/Vibrations.py#L311)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/Vibrations.py#L311?message=Update%20Docs)]
</div>

A Coordinerds CoordinateSystem object that manages all of the data needed to
work with normal mode coordinates + some convenience functions for generating and whatnot







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
name: str
```
<a id="Psience.Molecools.Vibrations.MolecularNormalModes.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, molecule, coeffs, name=None, freqs=None, internal=False, origin=None, basis=None, inverse=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/Vibrations/MolecularNormalModes.py#L317)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/Vibrations/MolecularNormalModes.py#L317?message=Update%20Docs)]
</div>

  - `molecule`: `AbstractMolecule`
    > 
  - `coeffs`: `Any`
    > 
  - `name`: `Any`
    > 
  - `freqs`: `Any`
    > 
  - `internal`: `Any`
    > 
  - `origin`: `Any`
    > 
  - `basis`: `Any`
    > 
  - `inverse`: `Any`
    >


<a id="Psience.Molecools.Vibrations.MolecularNormalModes.molecule" class="docs-object-method">&nbsp;</a> 
```python
@property
molecule(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/Vibrations/MolecularNormalModes.py#L372)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/Vibrations/MolecularNormalModes.py#L372?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Vibrations.MolecularNormalModes.change_mol" class="docs-object-method">&nbsp;</a> 
```python
change_mol(self, mol): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/Vibrations/MolecularNormalModes.py#L381)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/Vibrations/MolecularNormalModes.py#L381?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Vibrations.MolecularNormalModes.coords_by_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
coords_by_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/Vibrations/MolecularNormalModes.py#L393)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/Vibrations/MolecularNormalModes.py#L393?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Vibrations.MolecularNormalModes.modes_by_coords" class="docs-object-method">&nbsp;</a> 
```python
@property
modes_by_coords(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/Vibrations/MolecularNormalModes.py#L396)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/Vibrations/MolecularNormalModes.py#L396?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Vibrations.MolecularNormalModes.to_internals" class="docs-object-method">&nbsp;</a> 
```python
to_internals(self, intcrds=None, dYdR=None, dRdY=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/Vibrations/MolecularNormalModes.py#L401)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/Vibrations/MolecularNormalModes.py#L401?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Vibrations.MolecularNormalModes.origin" class="docs-object-method">&nbsp;</a> 
```python
@property
origin(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/Vibrations/MolecularNormalModes.py#L434)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/Vibrations/MolecularNormalModes.py#L434?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Vibrations.MolecularNormalModes.embed" class="docs-object-method">&nbsp;</a> 
```python
embed(self, frame): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/Vibrations/MolecularNormalModes.py#L444)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/Vibrations/MolecularNormalModes.py#L444?message=Update%20Docs)]
</div>

  - `frame`: `MolecularTransformation`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Vibrations.MolecularNormalModes.insert" class="docs-object-method">&nbsp;</a> 
```python
insert(self, val, where): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/Vibrations/MolecularNormalModes.py#L507)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/Vibrations/MolecularNormalModes.py#L507?message=Update%20Docs)]
</div>
Inserts values into the appropriate positions in the mode matrix
  - `val`: `Any`
    > 
  - `where`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Vibrations.MolecularNormalModes.to_new_modes" class="docs-object-method">&nbsp;</a> 
```python
to_new_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/Vibrations/MolecularNormalModes.py#L560)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/Vibrations/MolecularNormalModes.py#L560?message=Update%20Docs)]
</div>
Converts to the new generalized normal modes
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Vibrations.MolecularNormalModes.from_new_modes" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_new_modes(cls, mol, modes): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L580)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L580?message=Update%20Docs)]
</div>
Converts to the new generalized normal modes
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Vibrations.MolecularNormalModes.from_force_constants" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_force_constants(cls, molecule, fcs, *, atoms=None, masses=None, mass_units='AtomicMassUnits', inverse_mass_matrix=False, remove_transrot=True, dimensionless=False, mass_weighted=False, normalize=False, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L596)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L596?message=Update%20Docs)]
</div>
Generates normal modes from the specified force constants
  - `molecule`: `AbstractMolecule`
    > 
  - `fcs`: `np.ndarray`
    > force constants array
  - `atoms`: `Iterable[str]`
    > atom list
  - `masses`: `Iterable[float]`
    > mass list
  - `mass_units`: `str`
    > units for the masses...not clear if this is useful or a distraction
  - `inverse_mass_matrix`: `bool`
    > whether or not we have G or G^-1 (default: `False`)
  - `remove_transrot`: `bool`
    > whether or not to remove the translations and rotations (default: `True`)
  - `normalize`: `bool`
    > whether or not to normalize the modes (default: `True`)
  - `opts`: `Any`
    > 
  - `:returns`: `MolecularNormalModes`
    >


<a id="Psience.Molecools.Vibrations.MolecularNormalModes.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/Vibrations/MolecularNormalModes.py#L653)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/Vibrations/MolecularNormalModes.py#L653?message=Update%20Docs)]
</div>
Takes a slice of the modes
  - `item`: `Any`
    > 
  - `:returns`: `_`
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Molecools/Vibrations/MolecularNormalModes.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Molecools/Vibrations/MolecularNormalModes.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Molecools/Vibrations/MolecularNormalModes.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Molecools/Vibrations/MolecularNormalModes.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/Vibrations.py#L311?message=Update%20Docs)   
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