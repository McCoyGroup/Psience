## <a id="Psience.Molecools.Molecule.Molecule">Molecule</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L43)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L43?message=Update%20Docs)]
</div>

General purpose 'Molecule' class where the 'Molecule' need not be a molecule at all







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
bond_guessing_mode: str
default_energy_evalutor: str
highlight_styles: dict
default_display_mode: str
```
<a id="Psience.Molecools.Molecule.Molecule.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, atoms, coords, bonds=None, masses=None, name=None, internals=None, rdmol=None, dipole_surface=None, dipole_derivatives=None, potential_surface=None, potential_derivatives=None, normal_modes=None, source_file=None, guess_bonds=True, charge=None, display_mode=None, energy_evaluator=None, dipole_evaluator=None, charge_evaluator=None, **metadata): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L48)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L48?message=Update%20Docs)]
</div>

  - `atoms`: `Iterable[str]`
    > atoms specified by name, either full name or short
  - `coords`: `np.ndarray`
    > coordinates for the molecule, assumed to be in Bohr by default
  - `bonds`: `Iterable[Iterable[int]] | None`
    > bond specification for the molecule
  - `obmol`: `Any`
    > OpenBabel molecule for doing conversions
  - `charge`: `int | None`
    > Net charge on the molecule
  - `name`: `np.ndarray[int] | None`
    > Name for the molecule
The internal coordinate specification for the molecule
  - `dipole_surface`: `DipoleSurface | None`
    > The dipole surface for the system
  - `dipole_derivatives`: `Iterable[np.ndarray] | None`
    > Derivatives of the dipole surface
  - `potential_surface`: `PotentialSurface | None`
    > The potential surface for the system
  - `potential_derivatives`: `Iterable[np.ndarray] | None`
    > Derivatives of the potential surface
  - `guess_bonds`: `bool`
    > Whether or not to guess the bonding arrangement when that would be used
  - `source_file`: `str`
    > The data file the molecule was loaded from
  - `kw`: `Any`
    > Other bound parameters that might be useful


<a id="Psience.Molecools.Molecule.Molecule.modify" class="docs-object-method">&nbsp;</a> 
```python
modify(self, atoms=<McUtils.Devutils.core.DefaultType instance>, coords=<McUtils.Devutils.core.DefaultType instance>, *, internals=<McUtils.Devutils.core.DefaultType instance>, masses=<McUtils.Devutils.core.DefaultType instance>, bonds=<McUtils.Devutils.core.DefaultType instance>, guess_bonds=<McUtils.Devutils.core.DefaultType instance>, energy_evaluator=<McUtils.Devutils.core.DefaultType instance>, dipole_evaluator=<McUtils.Devutils.core.DefaultType instance>, charge_evaluator=<McUtils.Devutils.core.DefaultType instance>, display_mode=<McUtils.Devutils.core.DefaultType instance>, charge=<McUtils.Devutils.core.DefaultType instance>, normal_modes=<McUtils.Devutils.core.DefaultType instance>, dipole_surface=<McUtils.Devutils.core.DefaultType instance>, potential_surface=<McUtils.Devutils.core.DefaultType instance>, dipole_derivatives=<McUtils.Devutils.core.DefaultType instance>, potential_derivatives=<McUtils.Devutils.core.DefaultType instance>): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L145)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L145?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.canonicalize_internals" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
canonicalize_internals(cls, spec, atoms, coords, bonds, relocalize=True, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L250)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L250?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.prep_internal_spec" class="docs-object-method">&nbsp;</a> 
```python
prep_internal_spec(self, spec, relocalize=True, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L283)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L283?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.coords" class="docs-object-method">&nbsp;</a> 
```python
@property
coords(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L294)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L294?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.masses" class="docs-object-method">&nbsp;</a> 
```python
@property
masses(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L300)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L300?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.internals" class="docs-object-method">&nbsp;</a> 
```python
@property
internals(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L307)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L307?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.charge" class="docs-object-method">&nbsp;</a> 
```python
@property
charge(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L310)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L310?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.charges" class="docs-object-method">&nbsp;</a> 
```python
@property
charges(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L316)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L316?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_charge_evaluator" class="docs-object-method">&nbsp;</a> 
```python
get_charge_evaluator(self, evaluator=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L323)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L323?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.calculate_charges" class="docs-object-method">&nbsp;</a> 
```python
calculate_charges(self, evaluator=None, order=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L333)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L333?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.internal_coordinates" class="docs-object-method">&nbsp;</a> 
```python
@property
internal_coordinates(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L350)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L350?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.redundant_internal_transformation" class="docs-object-method">&nbsp;</a> 
```python
@property
redundant_internal_transformation(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L353)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L353?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.mode_embedding" class="docs-object-method">&nbsp;</a> 
```python
@property
mode_embedding(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L357)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L357?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_internals" class="docs-object-method">&nbsp;</a> 
```python
get_internals(self, strip_embedding=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L362)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L362?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_cartesians_by_internals" class="docs-object-method">&nbsp;</a> 
```python
get_cartesians_by_internals(self, order=None, strip_embedding=False, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L365)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L365?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_internals_by_cartesians" class="docs-object-method">&nbsp;</a> 
```python
get_internals_by_cartesians(self, order=None, strip_embedding=False, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L368)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L368?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_cartesians_by_modes" class="docs-object-method">&nbsp;</a> 
```python
get_cartesians_by_modes(self, order=None, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L371)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L371?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_modes_by_cartesians" class="docs-object-method">&nbsp;</a> 
```python
get_modes_by_cartesians(self, order=None, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L374)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L374?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.dipole_surface" class="docs-object-method">&nbsp;</a> 
```python
@property
dipole_surface(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L380)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L380?message=Update%20Docs)]
</div>

  - `:returns`: `DipoleSurfaceManager`
    >


<a id="Psience.Molecools.Molecule.Molecule.dipole_derivatives" class="docs-object-method">&nbsp;</a> 
```python
@property
dipole_derivatives(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L394)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L394?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_cartesian_dipole_derivatives" class="docs-object-method">&nbsp;</a> 
```python
get_cartesian_dipole_derivatives(self, order=None, evaluator=None, include_constant_term=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L400)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L400?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_internal_dipole_derivatives" class="docs-object-method">&nbsp;</a> 
```python
get_internal_dipole_derivatives(self, order=None, reembed=True, strip_embedding=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L423)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L423?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.potential_surface" class="docs-object-method">&nbsp;</a> 
```python
@property
potential_surface(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L433)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L433?message=Update%20Docs)]
</div>

  - `:returns`: `PotentialSurfaceManager`
    >


<a id="Psience.Molecools.Molecule.Molecule.potential_derivatives" class="docs-object-method">&nbsp;</a> 
```python
@property
potential_derivatives(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L447)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L447?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_cartesian_potential_derivatives" class="docs-object-method">&nbsp;</a> 
```python
get_cartesian_potential_derivatives(self, order=None, evaluator=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L471)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L471?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_internal_potential_derivatives" class="docs-object-method">&nbsp;</a> 
```python
get_internal_potential_derivatives(self, order=None, reembed=True, strip_embedding=True, zero_gradient=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L488)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L488?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.normal_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
normal_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L501)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L501?message=Update%20Docs)]
</div>

  - `:returns`: `NormalModesManager`
    >


<a id="Psience.Molecools.Molecule.Molecule.get_normal_modes" class="docs-object-method">&nbsp;</a> 
```python
get_normal_modes(self, masses=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L516)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L516?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_reaction_path_modes" class="docs-object-method">&nbsp;</a> 
```python
get_reaction_path_modes(self, masses=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L519)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L519?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.metadata" class="docs-object-method">&nbsp;</a> 
```python
@property
metadata(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L522)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L522?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_harmonic_spectrum" class="docs-object-method">&nbsp;</a> 
```python
get_harmonic_spectrum(self, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L534)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L534?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L539)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L539?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.num_atoms" class="docs-object-method">&nbsp;</a> 
```python
@property
num_atoms(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L560)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L560?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.atom_positions" class="docs-object-method">&nbsp;</a> 
```python
@property
atom_positions(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L563)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L563?message=Update%20Docs)]
</div>
A mapping of atom types to positions
  - `spec`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Molecule.Molecule.dummy_positions" class="docs-object-method">&nbsp;</a> 
```python
@property
dummy_positions(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L580)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L580?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.atoms" class="docs-object-method">&nbsp;</a> 
```python
@property
atoms(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L584)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L584?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.atomic_masses" class="docs-object-method">&nbsp;</a> 
```python
@property
atomic_masses(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L598)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L598?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.bonds" class="docs-object-method">&nbsp;</a> 
```python
@property
bonds(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L601)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L601?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.formula" class="docs-object-method">&nbsp;</a> 
```python
@property
formula(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L609)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L609?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.multiconfig" class="docs-object-method">&nbsp;</a> 
```python
@property
multiconfig(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L612)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L612?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.name" class="docs-object-method">&nbsp;</a> 
```python
@property
name(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L615)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L615?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.source_file" class="docs-object-method">&nbsp;</a> 
```python
@property
source_file(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L621)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L621?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.shape" class="docs-object-method">&nbsp;</a> 
```python
@property
shape(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L628)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L628?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.__len__" class="docs-object-method">&nbsp;</a> 
```python
__len__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L631)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L631?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.copy" class="docs-object-method">&nbsp;</a> 
```python
copy(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L637)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L637?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.take_submolecule" class="docs-object-method">&nbsp;</a> 
```python
take_submolecule(self, pos): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L653)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L653?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.prop" class="docs-object-method">&nbsp;</a> 
```python
prop(self, name, *args, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L677)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L677?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_guessed_bonds" class="docs-object-method">&nbsp;</a> 
```python
get_guessed_bonds(self, mode=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L689)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L689?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.edge_graph" class="docs-object-method">&nbsp;</a> 
```python
@property
edge_graph(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L704)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L704?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.fragment_indices" class="docs-object-method">&nbsp;</a> 
```python
@property
fragment_indices(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L708)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L708?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.fragments" class="docs-object-method">&nbsp;</a> 
```python
@property
fragments(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L712)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L712?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.mass_weighted_coords" class="docs-object-method">&nbsp;</a> 
```python
@property
mass_weighted_coords(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L717)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L717?message=Update%20Docs)]
</div>

  - `:returns`: `CoordinateSet`
    >


<a id="Psience.Molecools.Molecule.Molecule.center_of_mass" class="docs-object-method">&nbsp;</a> 
```python
@property
center_of_mass(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L725)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L725?message=Update%20Docs)]
</div>

  - `:returns`: `CoordinateSet`
    >


<a id="Psience.Molecools.Molecule.Molecule.inertia_tensor" class="docs-object-method">&nbsp;</a> 
```python
@property
inertia_tensor(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L732)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L732?message=Update%20Docs)]
</div>

  - `:returns`: `(np.ndarray, np.ndarray)`
    >


<a id="Psience.Molecools.Molecule.Molecule.inertial_eigensystem" class="docs-object-method">&nbsp;</a> 
```python
@property
inertial_eigensystem(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L739)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L739?message=Update%20Docs)]
</div>

  - `:returns`: `(np.ndarray, np.ndarray)`
    >


<a id="Psience.Molecools.Molecule.Molecule.moments_of_inertia" class="docs-object-method">&nbsp;</a> 
```python
@property
moments_of_inertia(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L746)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L746?message=Update%20Docs)]
</div>

  - `:returns`: `np.ndarray`
    >


<a id="Psience.Molecools.Molecule.Molecule.inertial_axes" class="docs-object-method">&nbsp;</a> 
```python
@property
inertial_axes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L753)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L753?message=Update%20Docs)]
</div>

  - `:returns`: `np.ndarray`
    >


<a id="Psience.Molecools.Molecule.Molecule.translation_rotation_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
translation_rotation_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L761)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L761?message=Update%20Docs)]
</div>

  - `:returns`: `np.ndarray`
    >


<a id="Psience.Molecools.Molecule.Molecule.get_translation_rotation_projector" class="docs-object-method">&nbsp;</a> 
```python
get_translation_rotation_projector(self, mass_weighted=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L769)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L769?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_translation_rotation_invariant_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_translation_rotation_invariant_transformation(self, mass_weighted=False, strip_embedding=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L783)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L783?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_energy_evaluator" class="docs-object-method">&nbsp;</a> 
```python
get_energy_evaluator(self, evaluator=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L817)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L817?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_energy_function" class="docs-object-method">&nbsp;</a> 
```python
get_energy_function(self, evaluator=None, *, order=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L830)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L830?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.calculate_energy" class="docs-object-method">&nbsp;</a> 
```python
calculate_energy(self, coords=None, *, evaluator=None, order=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L845)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L845?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.optimize" class="docs-object-method">&nbsp;</a> 
```python
optimize(self, evaluator=None, *, method=None, tol=None, max_iterations=None, logger=None, reembed=True, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L858)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L858?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_dipole_evaluator" class="docs-object-method">&nbsp;</a> 
```python
get_dipole_evaluator(self, evaluator=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L891)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L891?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_dipole_function" class="docs-object-method">&nbsp;</a> 
```python
get_dipole_function(self, evaluator=None, *, order=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L902)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L902?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.calculate_dipole" class="docs-object-method">&nbsp;</a> 
```python
calculate_dipole(self, evaluator=None, order=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L918)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L918?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_reduced_potential_generator" class="docs-object-method">&nbsp;</a> 
```python
get_reduced_potential_generator(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L929)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L929?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_1d_potentials" class="docs-object-method">&nbsp;</a> 
```python
get_1d_potentials(self, spec, evaluator=None, energy_expansion=None, potential_params=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L931)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L931?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.evaluate" class="docs-object-method">&nbsp;</a> 
```python
evaluate(self, func, use_internals=None, order=None, strip_embedding=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L946)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L946?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.evaluate_at" class="docs-object-method">&nbsp;</a> 
```python
evaluate_at(self, func, coords, use_internals=None, order=None, strip_embedding=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L958)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L958?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_displaced_coordinates" class="docs-object-method">&nbsp;</a> 
```python
get_displaced_coordinates(self, displacements, which=None, sel=None, axes=None, use_internals=False, coordinate_expansion=None, strip_embedding=False, shift=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L973)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L973?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_scan_coordinates" class="docs-object-method">&nbsp;</a> 
```python
get_scan_coordinates(self, domains, internals=False, modes=None, order=None, which=None, sel=None, axes=None, shift=True, coordinate_expansion=None, strip_embedding=False, return_displacements=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L989)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L989?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_nearest_displacement_atoms" class="docs-object-method">&nbsp;</a> 
```python
get_nearest_displacement_atoms(self, points, sel=None, axes=None, weighting_function=None, return_distances=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1024)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1024?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_nearest_displacement_coordinates" class="docs-object-method">&nbsp;</a> 
```python
get_nearest_displacement_coordinates(self, points, sel=None, axes=None, weighting_function=None, modes_nearest=False, return_distances=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1034)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1034?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_nearest_scan_coordinates" class="docs-object-method">&nbsp;</a> 
```python
get_nearest_scan_coordinates(self, domains, sel=None, axes=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1046)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1046?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.plot_molecule_function" class="docs-object-method">&nbsp;</a> 
```python
plot_molecule_function(self, function, *, axes, sel=None, embed=False, modes_nearest=False, domain=None, domain_padding=1, plot_points=500, weighting_function=None, mask_function=None, mask_value=0, plot_atoms=False, atom_colors=None, atom_radii=None, plotter=None, epilog=None, **plot_options): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1049)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1049?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_model" class="docs-object-method">&nbsp;</a> 
```python
get_model(self, potential_specs, dipole=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1161)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1161?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.setup_AIMD" class="docs-object-method">&nbsp;</a> 
```python
setup_AIMD(self, potential_function=None, timestep=0.5, seed=None, total_energy=None, total_energy_scaling=None, trajectories=1, sampled_modes=None, initial_energies=None, initial_displacements=None, initial_mode_directions=None, displaced_coords=None, track_kinetic_energy=False, track_velocities=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1286)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1286?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.setup_VPT" class="docs-object-method">&nbsp;</a> 
```python
setup_VPT(self, *, states=2, order=2, use_internals=None, potential_derivatives=None, energy_evaluator=None, dipole_derivatives=None, dipole_evaluator=None, runner='matrix', modes=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1374)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1374?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_gmatrix" class="docs-object-method">&nbsp;</a> 
```python
get_gmatrix(self, masses=None, coords=None, use_internals=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1429)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1429?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.g_matrix" class="docs-object-method">&nbsp;</a> 
```python
@property
g_matrix(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1448)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1448?message=Update%20Docs)]
</div>
Returns the molecular g-matrix for the system
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Molecule.Molecule.coriolis_constants" class="docs-object-method">&nbsp;</a> 
```python
@property
coriolis_constants(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1457)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1457?message=Update%20Docs)]
</div>
Returns the molecular g-matrix for the system
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Molecule.Molecule.bond_length" class="docs-object-method">&nbsp;</a> 
```python
bond_length(self, i, j): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1466)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1466?message=Update%20Docs)]
</div>
Returns the bond length of the coordinates
  - `i`: `Any`
    > 
  - `j`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Molecule.Molecule.bond_angle" class="docs-object-method">&nbsp;</a> 
```python
bond_angle(self, i, j, k): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1478)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1478?message=Update%20Docs)]
</div>
Returns the bond angle of the specified coordinates
  - `i`: `Any`
    > 
  - `j`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Molecule.Molecule.dihedral" class="docs-object-method">&nbsp;</a> 
```python
dihedral(self, i, j, k, l): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1490)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1490?message=Update%20Docs)]
</div>
Returns the dihedral angle of the specified coordinates
  - `i`: `Any`
    > 
  - `j`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Molecule.Molecule.principle_axis_frame" class="docs-object-method">&nbsp;</a> 
```python
principle_axis_frame(self, sel=None, inverse=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1506)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1506?message=Update%20Docs)]
</div>
Gets the principle axis frame(s) for the molecule
  - `mol`: `Any`
    > 
  - `sel`: `Any`
    > selection of atoms to use when getting the Eckart frame
  - `inverse`: `bool`
    > whether to return the inverse of the rotations or not
  - `:returns`: `MolecularTransformation | List[MolecularTransformation]`
    >


<a id="Psience.Molecools.Molecule.Molecule.principle_axis_data" class="docs-object-method">&nbsp;</a> 
```python
@property
principle_axis_data(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1520)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1520?message=Update%20Docs)]
</div>
Gets the principle axis embedded coords and embedding parameters for the molecule
  - `:returns`: `MolecularTransformation | List[MolecularTransformation]`
    >


<a id="Psience.Molecools.Molecule.Molecule.permute_atoms" class="docs-object-method">&nbsp;</a> 
```python
permute_atoms(self, perm): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1544)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1544?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.apply_affine_transformation" class="docs-object-method">&nbsp;</a> 
```python
apply_affine_transformation(self, transformation, load_properties=False, embed_properties=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1556)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1556?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.apply_rotation" class="docs-object-method">&nbsp;</a> 
```python
apply_rotation(self, rotation_matrix, shift_com=None, load_properties=False, embed_properties=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1581)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1581?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.eckart_frame" class="docs-object-method">&nbsp;</a> 
```python
eckart_frame(self, mol, sel=None, inverse=False, planar_ref_tolerance=None, proper_rotation=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1591)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1591?message=Update%20Docs)]
</div>
Gets the Eckart frame(s) for the molecule
  - `mol`: `Any`
    > 
  - `sel`: `Any`
    > selection of atoms to use when getting the Eckart frame
  - `inverse`: `bool`
    > whether to return the inverse of the rotations or not
  - `:returns`: `MolecularTransformation`
    >


<a id="Psience.Molecools.Molecule.Molecule.embed_coords" class="docs-object-method">&nbsp;</a> 
```python
embed_coords(self, crds, sel=None, in_paf=False, planar_ref_tolerance=None, proper_rotation=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1612)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1612?message=Update%20Docs)]
</div>
Embeds coords in the Eckart frame using `self` as a reference
  - `crds`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Molecule.Molecule.get_embedding_data" class="docs-object-method">&nbsp;</a> 
```python
get_embedding_data(self, crds, sel=None, in_paf=False, planar_ref_tolerance=None, proper_rotation=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1626)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1626?message=Update%20Docs)]
</div>
Gets the necessary data to embed crds in the Eckart frame using `self` as a reference
  - `crds`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Molecule.Molecule.get_embedded_molecule" class="docs-object-method">&nbsp;</a> 
```python
get_embedded_molecule(self, ref=None, sel=None, planar_ref_tolerance=None, proper_rotation=False, embed_properties=True, load_properties=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1638)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1638?message=Update%20Docs)]
</div>
Returns a Molecule embedded in an Eckart frame if ref is not None, otherwise returns
a principle-axis embedded Molecule
  - `:returns`: `Molecule`
    >


<a id="Psience.Molecools.Molecule.Molecule.align_molecule" class="docs-object-method">&nbsp;</a> 
```python
align_molecule(self, other: 'typing.Self', reindex_bonds=True, permute_atoms=True, align_structures=True, sel=None, embed_properties=True, load_properties=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1660)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1660?message=Update%20Docs)]
</div>
Aligns `other` with `self` by first finding the reindexing of the bonds of `other` that
lead to the best graph overlap with `self`, then determining which atoms can be permuted based on their graph
structures, then determining which permutation of equivalent atoms leads to the best agreement between the structures,
and then finally finding the Eckart/min-RMSD transformation after this transformation has been applied
  - `other`: `Any`
    > 
  - `reindex_bonds`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Molecule.Molecule.rdmol" class="docs-object-method">&nbsp;</a> 
```python
@property
rdmol(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1752)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1752?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.from_zmat" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_zmat(cls, zmat, internals=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1761)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1761?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.from_openbabel" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_openbabel(cls, mol, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1771)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1771?message=Update%20Docs)]
</div>

  - `mol`: `pybel.mol`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Molecule.Molecule.get_obmol" class="docs-object-method">&nbsp;</a> 
```python
get_obmol(self, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1790)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1790?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.from_rdmol" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_rdmol(cls, rdmol, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1825)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1825?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.from_name" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_name(cls, name, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1892)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1892?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_atom_strings" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_atom_strings(cls): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1897)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1897?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_string_format_dispatchers" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_string_format_dispatchers(cls): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1938)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1938?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.from_string" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_string(cls, string, fmt=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1948)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1948?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_file_format_dispatchers" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_file_format_dispatchers(cls): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1971)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1971?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.from_file" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_file(cls, file, mode=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1981)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1981?message=Update%20Docs)]
</div>
In general we'll delegate to pybel except for like Fchk and Log files
  - `file`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Molecule.Molecule.get_string_export_dispatchers" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_string_export_dispatchers(cls): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L2058)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L2058?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.to_string" class="docs-object-method">&nbsp;</a> 
```python
to_string(self, fmt, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2067)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2067?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_file_export_dispatchers" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_file_export_dispatchers(cls): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L2093)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L2093?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.to_file" class="docs-object-method">&nbsp;</a> 
```python
to_file(self, file, mode=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2097)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2097?message=Update%20Docs)]
</div>

  - `file`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Molecule.Molecule.construct" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
construct(cls, spec, fmt=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L2149)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L2149?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.plot" class="docs-object-method">&nbsp;</a> 
```python
plot(self, *geometries, figure=None, return_objects=False, bond_radius=0.1, atom_radius_scaling=0.25, atom_style=None, bond_style=None, vector_style=None, highlight_atoms=None, highlight_bonds=None, highlight_rings=None, highlight_styles=None, mode_vectors=None, mode_vector_origins=None, mode_vector_origin_mode='set', mode_vector_display_cutoff=0.01, dipole=None, dipole_origin=None, dipole_origin_mode='set', render_multiple_bonds=True, up_vector=None, multiple_bond_spacing=None, mode=None, backend=None, objects=False, graphics_class=None, cylinder_class=None, sphere_class=None, arrow_class=None, animate=None, animation_options=None, jsmol_load_script=None, units='Angstroms', **plot_ops): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2179)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2179?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_animation_geoms" class="docs-object-method">&nbsp;</a> 
```python
get_animation_geoms(self, which, extent=0.35, steps=8, strip_embedding=True, units=None, coordinate_expansion=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2564)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2564?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.animate_coordinate" class="docs-object-method">&nbsp;</a> 
```python
animate_coordinate(self, which, extent=0.5, steps=8, return_objects=False, strip_embedding=True, units='Angstroms', backend=None, mode=None, jsmol_load_script=None, coordinate_expansion=None, **plot_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2586)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2586?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.animate_mode" class="docs-object-method">&nbsp;</a> 
```python
animate_mode(self, which, extent=0.5, steps=8, modes=None, coordinate_expansion=None, order=None, normalize=True, mass_weight=False, mass_scale=True, frequency_scale=True, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2612)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2612?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.format_structs" class="docs-object-method">&nbsp;</a> 
```python
format_structs(self, geoms, format='xyz'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2662)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2662?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.format_animation_file" class="docs-object-method">&nbsp;</a> 
```python
format_animation_file(self, which, format='xyz', extent=0.35, steps=8, strip_embedding=True, units='Angstroms', coordinate_expansion=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2694)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2694?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.jsmol_viz" class="docs-object-method">&nbsp;</a> 
```python
jsmol_viz(self, xyz=None, animate=False, vibrate=False, script=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2711)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2711?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.jupyter_viz" class="docs-object-method">&nbsp;</a> 
```python
jupyter_viz(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2721)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2721?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.to_widget" class="docs-object-method">&nbsp;</a> 
```python
to_widget(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2729)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2729?message=Update%20Docs)]
</div>
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Molecools/Molecule/Molecule.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Molecools/Molecule/Molecule.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Molecools/Molecule/Molecule.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Molecools/Molecule/Molecule.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L43?message=Update%20Docs)   
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