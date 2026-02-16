## <a id="Psience.Molecools.Molecule.Molecule">Molecule</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L50)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L50?message=Update%20Docs)]
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
vector_style: dict
principle_axes_style: list
draw_coords_style: dict
draw_coords_label_style: dict
backend_options_resolution: dict
default_display_mode: str
```
<a id="Psience.Molecools.Molecule.Molecule.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, atoms, coords, bonds=None, masses=None, name=None, internals=None, rdmol=None, dipole_surface=None, dipole_derivatives=None, potential_surface=None, potential_derivatives=None, normal_modes=None, source_file=None, guess_bonds=True, charge=None, formal_charges=None, spin=None, display_mode=None, display_settings=None, energy=None, energy_evaluator=None, dipole_evaluator=None, charge_evaluator=None, polarizability_evaluator=None, polarizability_derivatives=None, checkpoint_file=None, **metadata): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L55)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L55?message=Update%20Docs)]
</div>

  - `atoms`: `Iterable[str]`
    > atoms specified by name, either full name or short
  - `coords`: `np.ndarray | Iterable[Iterable[float]]`
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
modify(self, atoms=<McUtils.Devutils.core.DefaultType instance>, coords=<McUtils.Devutils.core.DefaultType instance>, *, internals=<McUtils.Devutils.core.DefaultType instance>, masses=<McUtils.Devutils.core.DefaultType instance>, bonds=<McUtils.Devutils.core.DefaultType instance>, guess_bonds=<McUtils.Devutils.core.DefaultType instance>, energy=<McUtils.Devutils.core.DefaultType instance>, energy_evaluator=<McUtils.Devutils.core.DefaultType instance>, dipole_evaluator=<McUtils.Devutils.core.DefaultType instance>, charge_evaluator=<McUtils.Devutils.core.DefaultType instance>, polarizability_evaluator=<McUtils.Devutils.core.DefaultType instance>, charge=<McUtils.Devutils.core.DefaultType instance>, spin=<McUtils.Devutils.core.DefaultType instance>, rdmol=<McUtils.Devutils.core.DefaultType instance>, display_mode=<McUtils.Devutils.core.DefaultType instance>, display_settings=<McUtils.Devutils.core.DefaultType instance>, normal_modes=<McUtils.Devutils.core.DefaultType instance>, dipole_surface=<McUtils.Devutils.core.DefaultType instance>, potential_surface=<McUtils.Devutils.core.DefaultType instance>, dipole_derivatives=<McUtils.Devutils.core.DefaultType instance>, potential_derivatives=<McUtils.Devutils.core.DefaultType instance>, polarizability_derivatives=<McUtils.Devutils.core.DefaultType instance>, meta=<McUtils.Devutils.core.DefaultType instance>, source_file=<McUtils.Devutils.core.DefaultType instance>): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L174)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L174?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.to_state" class="docs-object-method">&nbsp;</a> 
```python
to_state(self, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L283)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L283?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.from_state" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_state(cls, data, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L331)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L331?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.cached_eval" class="docs-object-method">&nbsp;</a> 
```python
cached_eval(self, key, generator, *, condition=None, args=(), kwargs=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L335)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L335?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.canonicalize_internals" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
canonicalize_internals(cls, spec, atoms, coords, bonds, relocalize=True, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L416)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L416?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.prep_internal_spec" class="docs-object-method">&nbsp;</a> 
```python
prep_internal_spec(self, spec, relocalize=True, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L466)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L466?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.coords" class="docs-object-method">&nbsp;</a> 
```python
@property
coords(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L477)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L477?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.masses" class="docs-object-method">&nbsp;</a> 
```python
@property
masses(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L483)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L483?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.internals" class="docs-object-method">&nbsp;</a> 
```python
@property
internals(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L490)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L490?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.charge" class="docs-object-method">&nbsp;</a> 
```python
@property
charge(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L493)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L493?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.spin" class="docs-object-method">&nbsp;</a> 
```python
@property
spin(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L499)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L499?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.charges" class="docs-object-method">&nbsp;</a> 
```python
@property
charges(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L505)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L505?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.formal_charges" class="docs-object-method">&nbsp;</a> 
```python
@property
formal_charges(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L511)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L511?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_charge_evaluator" class="docs-object-method">&nbsp;</a> 
```python
get_charge_evaluator(self, evaluator=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L518)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L518?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.calculate_charges" class="docs-object-method">&nbsp;</a> 
```python
calculate_charges(self, evaluator=None, order=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L528)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L528?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.internal_coordinates" class="docs-object-method">&nbsp;</a> 
```python
@property
internal_coordinates(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L545)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L545?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.redundant_internal_transformation" class="docs-object-method">&nbsp;</a> 
```python
@property
redundant_internal_transformation(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L548)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L548?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_coordinate_filer" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_coordinate_filer(cls, allowed_coordinate_types=None, excluded_coordinate_types=None, allowed_ring_types=None, excluded_ring_types=None, allowed_group_types=None, excluded_group_types=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L575)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L575?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_bond_graph_internals" class="docs-object-method">&nbsp;</a> 
```python
get_bond_graph_internals(self, include_stretches=True, include_bends=True, include_dihedrals=True, include_fragments=True, pruning=None, fragment=None, base_internals=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L600)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L600?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.prune_internals" class="docs-object-method">&nbsp;</a> 
```python
prune_internals(self, coords, method='b_matrix'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L650)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L650?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_labeled_internals" class="docs-object-method">&nbsp;</a> 
```python
get_labeled_internals(self, coordinate_filter=None, allowed_coordinate_types=None, excluded_coordinate_types=None, allowed_ring_types=None, excluded_ring_types=None, allowed_group_types=None, excluded_group_types=None, include_stretches=True, include_bends=True, include_dihedrals=True, coordinate_sorting=None, pruning=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L674)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L674?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_mode_labels" class="docs-object-method">&nbsp;</a> 
```python
get_mode_labels(self, internals=None, modes=None, use_redundants=True, expansions=None, return_modes=False, **internals_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L725)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L725?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.mode_embedding" class="docs-object-method">&nbsp;</a> 
```python
@property
mode_embedding(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L798)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L798?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_internals" class="docs-object-method">&nbsp;</a> 
```python
get_internals(self, strip_embedding=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L803)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L803?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_cartesians_by_internals" class="docs-object-method">&nbsp;</a> 
```python
get_cartesians_by_internals(self, order=None, strip_embedding=False, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L806)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L806?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_internals_by_cartesians" class="docs-object-method">&nbsp;</a> 
```python
get_internals_by_cartesians(self, order=None, strip_embedding=False, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L809)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L809?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_cartesians_by_modes" class="docs-object-method">&nbsp;</a> 
```python
get_cartesians_by_modes(self, order=None, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L812)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L812?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_modes_by_cartesians" class="docs-object-method">&nbsp;</a> 
```python
get_modes_by_cartesians(self, order=None, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L815)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L815?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.dipole_surface" class="docs-object-method">&nbsp;</a> 
```python
@property
dipole_surface(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L821)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L821?message=Update%20Docs)]
</div>

  - `:returns`: `DipoleSurfaceManager`
    >


<a id="Psience.Molecools.Molecule.Molecule.dipole_derivatives" class="docs-object-method">&nbsp;</a> 
```python
@property
dipole_derivatives(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L835)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L835?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_cartesian_dipole_derivatives" class="docs-object-method">&nbsp;</a> 
```python
get_cartesian_dipole_derivatives(self, order=None, evaluator=None, include_constant_term=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L841)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L841?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_internal_dipole_derivatives" class="docs-object-method">&nbsp;</a> 
```python
get_internal_dipole_derivatives(self, order=None, reembed=True, strip_embedding=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L867)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L867?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_cartesian_polarizability_derivatives" class="docs-object-method">&nbsp;</a> 
```python
get_cartesian_polarizability_derivatives(self, order=None, evaluator=None, include_constant_term=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L876)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L876?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_hamiltonian" class="docs-object-method">&nbsp;</a> 
```python
get_hamiltonian(self, embedding=None, potential_derivatives=None, modes=None, dipole_derivatives=None, **etc): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L904)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L904?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.potential_surface" class="docs-object-method">&nbsp;</a> 
```python
@property
potential_surface(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L926)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L926?message=Update%20Docs)]
</div>

  - `:returns`: `PotentialSurfaceManager`
    >


<a id="Psience.Molecools.Molecule.Molecule.potential_derivatives" class="docs-object-method">&nbsp;</a> 
```python
@property
potential_derivatives(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L940)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L940?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_cartesian_potential_derivatives" class="docs-object-method">&nbsp;</a> 
```python
get_cartesian_potential_derivatives(self, order=None, evaluator=None, use_cached=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L964)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L964?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_internal_potential_derivatives" class="docs-object-method">&nbsp;</a> 
```python
get_internal_potential_derivatives(self, order=None, reembed=True, strip_embedding=True, zero_gradient=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L981)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L981?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.normal_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
normal_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L994)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L994?message=Update%20Docs)]
</div>

  - `:returns`: `NormalModesManager`
    >


<a id="Psience.Molecools.Molecule.Molecule.get_normal_modes" class="docs-object-method">&nbsp;</a> 
```python
get_normal_modes(self, masses=None, potential_derivatives=None, use_internals=None, project_transrot=True, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1009)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1009?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_reaction_path_modes" class="docs-object-method">&nbsp;</a> 
```python
get_reaction_path_modes(self, masses=None, potential_derivatives=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1022)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1022?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.metadata" class="docs-object-method">&nbsp;</a> 
```python
@property
metadata(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1030)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1030?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_harmonic_spectrum" class="docs-object-method">&nbsp;</a> 
```python
get_harmonic_spectrum(self, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1042)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1042?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_harmonic_raman_spectrum" class="docs-object-method">&nbsp;</a> 
```python
get_harmonic_raman_spectrum(self, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1046)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1046?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1051)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1051?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.num_atoms" class="docs-object-method">&nbsp;</a> 
```python
@property
num_atoms(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1072)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1072?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.atom_positions" class="docs-object-method">&nbsp;</a> 
```python
@property
atom_positions(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1075)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1075?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.dummy_positions" class="docs-object-method">&nbsp;</a> 
```python
@property
dummy_positions(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1084)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1084?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.atoms" class="docs-object-method">&nbsp;</a> 
```python
@property
atoms(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1088)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1088?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.atomic_masses" class="docs-object-method">&nbsp;</a> 
```python
@property
atomic_masses(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1102)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1102?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.bonds" class="docs-object-method">&nbsp;</a> 
```python
@property
bonds(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1105)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1105?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.break_bonds" class="docs-object-method">&nbsp;</a> 
```python
break_bonds(self, bonds): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1113)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1113?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.formula" class="docs-object-method">&nbsp;</a> 
```python
@property
formula(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1119)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1119?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.multiconfig" class="docs-object-method">&nbsp;</a> 
```python
@property
multiconfig(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1122)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1122?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.name" class="docs-object-method">&nbsp;</a> 
```python
@property
name(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1125)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1125?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.source_file" class="docs-object-method">&nbsp;</a> 
```python
@property
source_file(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1131)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1131?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.source_mode" class="docs-object-method">&nbsp;</a> 
```python
@property
source_mode(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1148)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1148?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.shape" class="docs-object-method">&nbsp;</a> 
```python
@property
shape(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1153)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1153?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.__len__" class="docs-object-method">&nbsp;</a> 
```python
__len__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1156)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1156?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.copy" class="docs-object-method">&nbsp;</a> 
```python
copy(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1162)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1162?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.take_submolecule" class="docs-object-method">&nbsp;</a> 
```python
take_submolecule(self, pos): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1180)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1180?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.prop" class="docs-object-method">&nbsp;</a> 
```python
prop(self, name, *args, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1204)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1204?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_guessed_bonds" class="docs-object-method">&nbsp;</a> 
```python
get_guessed_bonds(self, mode=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1216)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1216?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.edge_graph" class="docs-object-method">&nbsp;</a> 
```python
@property
edge_graph(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1232)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1232?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.find_substructure" class="docs-object-method">&nbsp;</a> 
```python
find_substructure(self, pattern): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1236)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1236?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.find_heavy_atom_backbone" class="docs-object-method">&nbsp;</a> 
```python
find_heavy_atom_backbone(self, root=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1239)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1239?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.find_backbone_segments" class="docs-object-method">&nbsp;</a> 
```python
find_backbone_segments(self, root=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1242)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1242?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_backbone_zmatrix" class="docs-object-method">&nbsp;</a> 
```python
get_backbone_zmatrix(self, root=None, segments=None, return_remainder=False, return_segments=False, validate=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1245)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1245?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_canonical_zmatrix" class="docs-object-method">&nbsp;</a> 
```python
get_canonical_zmatrix(self, ordering=None, validate=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1280)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1280?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_bond_zmatrix" class="docs-object-method">&nbsp;</a> 
```python
get_bond_zmatrix(self, fragments=None, segments=None, root=None, attachment_points=None, check_attachment_points=True, validate=True, for_fragment=None, fragment_ordering=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1285)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1285?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.fragment_indices" class="docs-object-method">&nbsp;</a> 
```python
@property
fragment_indices(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1386)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1386?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.fragments" class="docs-object-method">&nbsp;</a> 
```python
@property
fragments(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1390)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1390?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.mass_weighted_coords" class="docs-object-method">&nbsp;</a> 
```python
@property
mass_weighted_coords(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1395)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1395?message=Update%20Docs)]
</div>

  - `:returns`: `CoordinateSet`
    >


<a id="Psience.Molecools.Molecule.Molecule.center_of_mass" class="docs-object-method">&nbsp;</a> 
```python
@property
center_of_mass(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1403)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1403?message=Update%20Docs)]
</div>

  - `:returns`: `CoordinateSet`
    >


<a id="Psience.Molecools.Molecule.Molecule.inertia_tensor" class="docs-object-method">&nbsp;</a> 
```python
@property
inertia_tensor(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1410)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1410?message=Update%20Docs)]
</div>

  - `:returns`: `(np.ndarray, np.ndarray)`
    >


<a id="Psience.Molecools.Molecule.Molecule.inertial_eigensystem" class="docs-object-method">&nbsp;</a> 
```python
@property
inertial_eigensystem(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1417)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1417?message=Update%20Docs)]
</div>

  - `:returns`: `(np.ndarray, np.ndarray)`
    >


<a id="Psience.Molecools.Molecule.Molecule.moments_of_inertia" class="docs-object-method">&nbsp;</a> 
```python
@property
moments_of_inertia(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1424)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1424?message=Update%20Docs)]
</div>

  - `:returns`: `np.ndarray`
    >


<a id="Psience.Molecools.Molecule.Molecule.inertial_axes" class="docs-object-method">&nbsp;</a> 
```python
@property
inertial_axes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1431)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1431?message=Update%20Docs)]
</div>

  - `:returns`: `np.ndarray`
    >


<a id="Psience.Molecools.Molecule.Molecule.translation_rotation_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
translation_rotation_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1439)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1439?message=Update%20Docs)]
</div>

  - `:returns`: `np.ndarray`
    >


<a id="Psience.Molecools.Molecule.Molecule.get_translation_rotation_projector" class="docs-object-method">&nbsp;</a> 
```python
get_translation_rotation_projector(self, mass_weighted=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1447)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1447?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_translation_rotation_invariant_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_translation_rotation_invariant_transformation(self, mass_weighted=False, strip_embedding=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1461)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1461?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.energy" class="docs-object-method">&nbsp;</a> 
```python
@property
energy(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1512)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1512?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_energy_evaluator" class="docs-object-method">&nbsp;</a> 
```python
get_energy_evaluator(self, evaluator=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1521)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1521?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_energy_function" class="docs-object-method">&nbsp;</a> 
```python
get_energy_function(self, evaluator=None, *, order=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1534)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1534?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.calculate_energy" class="docs-object-method">&nbsp;</a> 
```python
calculate_energy(self, coords=None, *, evaluator=None, order=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1549)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1549?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.partial_force_field" class="docs-object-method">&nbsp;</a> 
```python
partial_force_field(self, coords=None, modes=None, *, evaluator=None, order=4, mesh_spacing=1, analytic_derivative_order=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1561)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1561?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.optimize" class="docs-object-method">&nbsp;</a> 
```python
optimize(self, evaluator=None, *, method=None, tol=None, max_iterations=None, logger=None, reembed=True, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1582)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1582?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_dipole_evaluator" class="docs-object-method">&nbsp;</a> 
```python
get_dipole_evaluator(self, evaluator=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1619)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1619?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_dipole_function" class="docs-object-method">&nbsp;</a> 
```python
get_dipole_function(self, evaluator=None, *, order=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1630)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1630?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.calculate_dipole" class="docs-object-method">&nbsp;</a> 
```python
calculate_dipole(self, evaluator=None, order=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1646)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1646?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_polarizability_evaluator" class="docs-object-method">&nbsp;</a> 
```python
get_polarizability_evaluator(self, evaluator=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1657)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1657?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.calculate_dipole_polarizability" class="docs-object-method">&nbsp;</a> 
```python
calculate_dipole_polarizability(self, evaluator=None, order=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1667)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1667?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.polarizability_derivatives" class="docs-object-method">&nbsp;</a> 
```python
@property
polarizability_derivatives(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1677)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1677?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_reduced_potential_generator" class="docs-object-method">&nbsp;</a> 
```python
get_reduced_potential_generator(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1684)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1684?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_1d_potentials" class="docs-object-method">&nbsp;</a> 
```python
get_1d_potentials(self, spec, evaluator=None, energy_expansion=None, potential_params=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1686)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1686?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.evaluate" class="docs-object-method">&nbsp;</a> 
```python
evaluate(self, func, use_internals=None, order=None, strip_embedding=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1701)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1701?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.evaluate_at" class="docs-object-method">&nbsp;</a> 
```python
evaluate_at(self, func, coords, use_internals=None, order=None, strip_embedding=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1713)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1713?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_displaced_coordinates" class="docs-object-method">&nbsp;</a> 
```python
get_displaced_coordinates(self, displacements, which=None, sel=None, axes=None, use_internals=False, coordinate_expansion=None, strip_embedding=False, shift=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1728)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1728?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_scan_coordinates" class="docs-object-method">&nbsp;</a> 
```python
get_scan_coordinates(self, domains, internals=False, modes=None, order=None, which=None, sel=None, axes=None, shift=True, coordinate_expansion=None, strip_embedding=False, return_displacements=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1744)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1744?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_nearest_displacement_atoms" class="docs-object-method">&nbsp;</a> 
```python
get_nearest_displacement_atoms(self, points, sel=None, axes=None, weighting_function=None, return_distances=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1779)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1779?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_nearest_displacement_coordinates" class="docs-object-method">&nbsp;</a> 
```python
get_nearest_displacement_coordinates(self, points, sel=None, axes=None, weighting_function=None, modes_nearest=False, return_distances=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1789)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1789?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_nearest_scan_coordinates" class="docs-object-method">&nbsp;</a> 
```python
get_nearest_scan_coordinates(self, domains, sel=None, axes=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1801)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1801?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.plot_molecule_function" class="docs-object-method">&nbsp;</a> 
```python
plot_molecule_function(self, function, *, axes, sel=None, embed=False, modes_nearest=False, domain=None, domain_padding=1, plot_points=500, weighting_function=None, mask_function=None, mask_value=0, plot_atoms=False, atom_colors=None, atom_radii=None, plotter=None, epilog=None, **plot_options): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1813)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1813?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_model" class="docs-object-method">&nbsp;</a> 
```python
get_model(self, potential_specs, dipole=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1929)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1929?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_point_group" class="docs-object-method">&nbsp;</a> 
```python
get_point_group(self, coords=None, masses=None, *, sel=None, verbose=False, return_identifier=False, **tols): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2053)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2053?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.point_group" class="docs-object-method">&nbsp;</a> 
```python
@property
point_group(self) -> McUtils.Symmetry.PointGroups.PointGroup: 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2072)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2072?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_point_group_embedded_coordinates" class="docs-object-method">&nbsp;</a> 
```python
get_point_group_embedded_coordinates(self, pg=None, sel=None, return_point_group=False, return_identifier=False, **tols): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2171)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2171?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.symmetrize" class="docs-object-method">&nbsp;</a> 
```python
symmetrize(self, pg=None, return_identifier=False, tol=0.1, sel=None, return_coordinates=None, return_point_group=False, **tols): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2206)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2206?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_symmetrized_internals" class="docs-object-method">&nbsp;</a> 
```python
get_symmetrized_internals(self, point_group=None, *, internals=None, extra_internals=None, masses=None, return_expansions=False, atom_selection=None, as_characters=True, normalize=None, drop_empty_modes=None, perms=None, return_base_expansion=False, return_point_group=False, reduce_redundant_coordinates=None, ops=None, permutation_tol=0.01, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2263)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2263?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_surface" class="docs-object-method">&nbsp;</a> 
```python
get_surface(self, radius_type='VanDerWaalsRadius', *, surface_type=None, radius_units='Angstroms', samples=50, radius_scaling=1, **etc): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2329)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2329?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_surface_mesh" class="docs-object-method">&nbsp;</a> 
```python
get_surface_mesh(self, radius_type='VanDerWaalsRadius', *, surface_type=None, radius_units='Angstroms', samples=50, expansion=0.01, mesh_options=None, **etc): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2352)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2352?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.setup_AIMD" class="docs-object-method">&nbsp;</a> 
```python
setup_AIMD(self, potential_function=None, timestep=0.5, seed=None, total_energy=None, total_energy_scaling=None, trajectories=1, sampled_modes=None, initial_energies=None, initial_displacements=None, initial_mode_directions=None, displaced_coords=None, track_kinetic_energy=False, track_velocities=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2377)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2377?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.setup_VPT" class="docs-object-method">&nbsp;</a> 
```python
setup_VPT(self, *, states=2, order=2, use_internals=None, potential_derivatives=None, energy_evaluator=None, dipole_derivatives=None, dipole_evaluator=None, runner='matrix', use_reaction_path=False, modes=None, projected_modes=None, mode_transformation=None, potential_terms=None, dipole_terms=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2474)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2474?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_gmatrix" class="docs-object-method">&nbsp;</a> 
```python
get_gmatrix(self, masses=None, coords=None, use_internals=None, power=None, **internals_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2569)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2569?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.g_matrix" class="docs-object-method">&nbsp;</a> 
```python
@property
g_matrix(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2607)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2607?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2616)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2616?message=Update%20Docs)]
</div>
Returns the molecular g-matrix for the system
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Molecule.Molecule.principle_axis_frame" class="docs-object-method">&nbsp;</a> 
```python
principle_axis_frame(self, sel=None, inverse=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2625)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2625?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2639)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2639?message=Update%20Docs)]
</div>
Gets the principle axis embedded coords and embedding parameters for the molecule
  - `:returns`: `MolecularTransformation | List[MolecularTransformation]`
    >


<a id="Psience.Molecools.Molecule.Molecule.permute_atoms" class="docs-object-method">&nbsp;</a> 
```python
permute_atoms(self, perm): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2649)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2649?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.apply_affine_transformation" class="docs-object-method">&nbsp;</a> 
```python
apply_affine_transformation(self, transformation, load_properties=False, embed_properties=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2662)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2662?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.apply_rotation" class="docs-object-method">&nbsp;</a> 
```python
apply_rotation(self, rotation_matrix, shift_com=None, load_properties=None, embed_properties=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2706)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2706?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.eckart_frame" class="docs-object-method">&nbsp;</a> 
```python
eckart_frame(self, mol, sel=None, inverse=False, planar_ref_tolerance=None, proper_rotation=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2716)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2716?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2737)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2737?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2751)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2751?message=Update%20Docs)]
</div>
Gets the necessary data to embed crds in the Eckart frame using `self` as a reference
  - `crds`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Molecule.Molecule.get_embedded_molecule" class="docs-object-method">&nbsp;</a> 
```python
get_embedded_molecule(self, ref=None, sel=None, planar_ref_tolerance=None, proper_rotation=False, embed_properties=True, load_properties=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2763)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2763?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2786)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2786?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2876)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2876?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.from_zmat" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_zmat(cls, zmat, internals=None, axes=None, origin=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L2887)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L2887?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.from_openbabel" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_openbabel(cls, mol, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L2899)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L2899?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2918)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2918?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.from_rdmol" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_rdmol(cls, rdmol, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L2995)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L2995?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.from_name" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_name(cls, name, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L3255)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L3255?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_atom_strings" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_atom_strings(cls): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L3260)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L3260?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_string_format_dispatchers" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_string_format_dispatchers(cls): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L3317)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L3317?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.from_string" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_string(cls, string, fmt=None, format_options=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L3333)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L3333?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_file_format_dispatchers" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_file_format_dispatchers(cls): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L3358)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L3358?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.from_file" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_file(cls, file, mode=None, format_options=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L3372)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L3372?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L3497)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L3497?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.to_string" class="docs-object-method">&nbsp;</a> 
```python
to_string(self, fmt, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L3509)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L3509?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_file_export_dispatchers" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_file_export_dispatchers(cls): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L3535)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L3535?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.to_file" class="docs-object-method">&nbsp;</a> 
```python
to_file(self, file, mode=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L3539)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L3539?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L3591)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L3591?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.plot" class="docs-object-method">&nbsp;</a> 
```python
plot(self, *geometries, figure=None, return_objects=False, bonds=None, bond_radius=0.1, atom_radius_scaling=0.25, atom_style=None, atom_radii=None, atom_text=None, display_atom_numbers=False, radius_type=None, bond_style=None, capped_bonds=False, reflectiveness=None, vector_style=None, highlight_atoms=None, highlight_bonds=None, highlight_rings=None, highlight_styles=None, mode_vectors=None, mode_vector_origins=None, mode_vector_origin_mode='set', mode_vector_display_cutoff=0.01, principle_axes=None, principle_axes_origin=None, principle_axes_origin_mode='set', principle_axes_style=None, dipole=None, dipole_origin=None, dipole_origin_mode='set', render_multiple_bonds=True, render_fractional_bonds=True, draw_coords=None, draw_coords_style=None, up_vector=None, multiple_bond_spacing=None, mode=None, backend=None, include_save_buttons=None, objects=False, graphics_class=None, cylinder_class=None, sphere_class=None, arrow_class=None, animate=None, recording_options=None, animation_options=None, jsmol_load_script=None, include_jsmol_script_interface=False, dynamic_loading=None, units='Angstroms', **plot_ops): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L4085)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L4085?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.get_animation_geoms" class="docs-object-method">&nbsp;</a> 
```python
get_animation_geoms(self, which, extent=0.35, steps=8, strip_embedding=True, units=None, coordinate_expansion=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L5015)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L5015?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.animate_coordinate" class="docs-object-method">&nbsp;</a> 
```python
animate_coordinate(self, which, extent=0.5, steps=3, return_objects=False, strip_embedding=True, units='Angstroms', backend=None, mode=None, jsmol_load_script=None, coordinate_expansion=None, **plot_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L5037)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L5037?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.animate_mode" class="docs-object-method">&nbsp;</a> 
```python
animate_mode(self, which, extent=0.5, steps=3, modes=None, coordinate_expansion=None, order=None, normalize=True, mass_weight=False, mass_scale=True, frequency_scale=False, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L5063)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L5063?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.format_structs" class="docs-object-method">&nbsp;</a> 
```python
format_structs(self, geoms, format='xyz'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L5113)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L5113?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.format_animation_file" class="docs-object-method">&nbsp;</a> 
```python
format_animation_file(self, which, format='xyz', extent=0.35, steps=8, strip_embedding=True, units='Angstroms', coordinate_expansion=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L5145)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L5145?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.jsmol_viz" class="docs-object-method">&nbsp;</a> 
```python
jsmol_viz(self, xyz=None, animate=False, vibrate=False, script=None, include_script_interface=False, image_size=None, width=None, height=None, **etc): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L5162)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L5162?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.jupyter_viz" class="docs-object-method">&nbsp;</a> 
```python
jupyter_viz(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L5198)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L5198?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Molecule.Molecule.to_widget" class="docs-object-method">&nbsp;</a> 
```python
to_widget(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L5206)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L5206?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L50?message=Update%20Docs)   
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