## <a id="Psience.Molecools.Vibrations.MolecularNormalModes">MolecularNormalModes</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations.py#L413)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations.py#L413?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations.py#L419)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations.py#L419?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularNormalModes.py#L474)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularNormalModes.py#L474?message=Update%20Docs)]
</div>
**LLM Docstring**

Property getter/setter for the molecule these normal modes belong to. Setting it also propagates the new molecule to `self.basis.molecule` before updating `self._molecule`.
  - `mol`: `AbstractMolecule`
    > (setter only) the new molecule to associate
  - `:returns`: `AbstractMolecule`
    > (getter) the stored molecule


<a id="Psience.Molecools.Vibrations.MolecularNormalModes.change_mol" class="docs-object-method">&nbsp;</a> 
```python
change_mol(self, mol): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularNormalModes.py#L503)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularNormalModes.py#L503?message=Update%20Docs)]
</div>
**LLM Docstring**

Rebind these normal modes to a different molecule, keeping the mode matrix, name, frequencies, internal-coordinate flag, origin, basis, and inverse the same.
  - `mol`: `AbstractMolecule`
    > the new molecule to associate with the normal modes
  - `:returns`: `MolecularNormalModes`
    > a new `MolecularNormalModes` for `mol` built from this object's `matrix`, `name`, `freqs`, `in_internals`, `_origin`, `basis`, and `inverse`


<a id="Psience.Molecools.Vibrations.MolecularNormalModes.coords_by_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
coords_by_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularNormalModes.py#L525)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularNormalModes.py#L525?message=Update%20Docs)]
</div>
**LLM Docstring**

Matrix mapping mode displacements back to coordinate displacements (i.e. the inverse transformation).
  - `:returns`: `np.ndarray`
    > the stored inverse matrix, `self.inverse`


<a id="Psience.Molecools.Vibrations.MolecularNormalModes.modes_by_coords" class="docs-object-method">&nbsp;</a> 
```python
@property
modes_by_coords(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularNormalModes.py#L536)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularNormalModes.py#L536?message=Update%20Docs)]
</div>
**LLM Docstring**

Matrix mapping coordinate displacements onto normal-mode displacements.
  - `:returns`: `np.ndarray`
    > the stored mode matrix, `self.matrix`


<a id="Psience.Molecools.Vibrations.MolecularNormalModes.to_internals" class="docs-object-method">&nbsp;</a> 
```python
to_internals(self, intcrds=None, dYdR=None, dRdY=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularNormalModes.py#L549)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularNormalModes.py#L549?message=Update%20Docs)]
</div>
**LLM Docstring**

Intended to convert Cartesian-basis normal modes into an internal-coordinate representation using the supplied Jacobians. The method immediately raises `NotImplementedError` directing callers to use the newer `NormalModes` object instead, so the remaining body (computing `dQdR`/`dRdQ` and constructing a new internal-coordinate `MolecularNormalModes`) is dead code that never executes.
  - `intcrds`: `CoordinateSet | None`
    > internal-coordinate system to convert into; if `None`, taken from `self.molecule.internal_coordinates`
  - `dYdR`: `np.ndarray | None`
    > Jacobian of mass-weighted Cartesians with respect to internal coordinates; computed from `intcrds`/`molecule` if not given (unreachable)
  - `dRdY`: `np.ndarray | None`
    > Jacobian of internal coordinates with respect to mass-weighted Cartesians; computed if not given (unreachable)
  - `:returns`: `None`
    > never returns; always raises


<a id="Psience.Molecools.Vibrations.MolecularNormalModes.origin" class="docs-object-method">&nbsp;</a> 
```python
@property
origin(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularNormalModes.py#L597)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularNormalModes.py#L597?message=Update%20Docs)]
</div>
**LLM Docstring**

Reference geometry the normal modes are expanded about. If no explicit origin was stored, falls back to the molecule's internal coordinates (if `in_internals`) or Cartesian coordinates otherwise.
  - `:returns`: `CoordinateSet`
    > the origin coordinates


<a id="Psience.Molecools.Vibrations.MolecularNormalModes.embed" class="docs-object-method">&nbsp;</a> 
```python
embed(self, frame): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularNormalModes.py#L615)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularNormalModes.py#L615?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularNormalModes.py#L678)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularNormalModes.py#L678?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularNormalModes.py#L731)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularNormalModes.py#L731?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L751)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L751?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L767)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L767?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Vibrations/MolecularNormalModes.py#L824)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations/MolecularNormalModes.py#L824?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Vibrations.py#L413?message=Update%20Docs)   
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