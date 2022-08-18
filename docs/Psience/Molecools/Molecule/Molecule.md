## <a id="Psience.Molecools.Molecule.Molecule">Molecule</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L30)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L30?message=Update%20Docs)]
</div>

General purpose 'Molecule' class where the 'Molecule' need not be a molecule at all



<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

<a id="Psience.Molecools.Molecule.Molecule.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, atoms, coords, bonds=None, masses=None, name=None, internals=None, obmol=None, dipole_surface=None, dipole_derivatives=None, potential_surface=None, potential_derivatives=None, normal_modes=None, source_file=None, guess_bonds=True, charge=None, **metadata): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L35)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L35?message=Update%20Docs)]
</div>


- `kw`: `Any`
    >Other bound parameters that might be useful
- `source_file`: `str`
    >The data file the molecule was loaded from
- `guess_bonds`: `bool`
    >Whether or not to guess the bonding arrangement when that would be used
- `potential_derivatives`: `Iterable[np.ndarray] | None`
    >Derivatives of the potential surface
- `potential_surface`: `PotentialSurface | None`
    >The potential surface for the system
- `dipole_derivatives`: `Iterable[np.ndarray] | None`
    >Derivatives of the dipole surface
- `dipole_surface`: `DipoleSurface | None`
    >The dipole surface for the system
- `name`: `np.ndarray[int] | None`
    >Name for the molecule
The internal coordinate specification for the molecule
- `charge`: `int | None`
    >Net charge on the molecule
- `obmol`: `Any`
    >OpenBabel molecule for doing conversions
- `bonds`: `Iterable[Iterable[int]] | None`
    >bond specification for the molecule
- `coords`: `np.ndarray`
    >coordinates for the molecule, assumed to be in Bohr by default
- `atoms`: `Iterable[str]`
    >atoms specified by name, either full name or short

<a id="Psience.Molecools.Molecule.Molecule.dipole_surface" class="docs-object-method">&nbsp;</a> 
```python
@property
dipole_surface(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>


- `:returns`: `DipoleSurfaceManager`
    >

<a id="Psience.Molecools.Molecule.Molecule.dipole_derivatives" class="docs-object-method">&nbsp;</a> 
```python
@property
dipole_derivatives(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.potential_surface" class="docs-object-method">&nbsp;</a> 
```python
@property
potential_surface(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>


- `:returns`: `PotentialSurfaceManager`
    >

<a id="Psience.Molecools.Molecule.Molecule.potential_derivatives" class="docs-object-method">&nbsp;</a> 
```python
@property
potential_derivatives(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.normal_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
normal_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>


- `:returns`: `NormalModesManager`
    >

<a id="Psience.Molecools.Molecule.Molecule.metadata" class="docs-object-method">&nbsp;</a> 
```python
@property
metadata(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L190)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L190?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.num_atoms" class="docs-object-method">&nbsp;</a> 
```python
@property
num_atoms(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.atom_positions" class="docs-object-method">&nbsp;</a> 
```python
@property
atom_positions(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>

A mapping of atom types to positions
- `:returns`: `_`
    >
- `spec`: `Any`
    >

<a id="Psience.Molecools.Molecule.Molecule.dummy_positions" class="docs-object-method">&nbsp;</a> 
```python
@property
dummy_positions(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.atoms" class="docs-object-method">&nbsp;</a> 
```python
@property
atoms(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.masses" class="docs-object-method">&nbsp;</a> 
```python
@property
masses(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.bonds" class="docs-object-method">&nbsp;</a> 
```python
@property
bonds(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.coords" class="docs-object-method">&nbsp;</a> 
```python
@property
coords(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.sys" class="docs-object-method">&nbsp;</a> 
```python
@property
sys(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.formula" class="docs-object-method">&nbsp;</a> 
```python
@property
formula(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.multiconfig" class="docs-object-method">&nbsp;</a> 
```python
@property
multiconfig(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.name" class="docs-object-method">&nbsp;</a> 
```python
@property
name(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.source_file" class="docs-object-method">&nbsp;</a> 
```python
@property
source_file(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.insert_atoms" class="docs-object-method">&nbsp;</a> 
```python
insert_atoms(self, atoms, coords, where, handle_properties=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L274)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L274?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.delete_atoms" class="docs-object-method">&nbsp;</a> 
```python
delete_atoms(self, where, handle_properties=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L300)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L300?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.take_submolecule" class="docs-object-method">&nbsp;</a> 
```python
take_submolecule(self, spec): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L316)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L316?message=Update%20Docs)]
</div>

Takes a 'slice' of a molecule if working with Cartesian coords.
If not, need to do some corner case handling for that.
- `:returns`: `_`
    >
- `spec`: `Any`
    >

<a id="Psience.Molecools.Molecule.Molecule.shape" class="docs-object-method">&nbsp;</a> 
```python
@property
shape(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.__len__" class="docs-object-method">&nbsp;</a> 
```python
__len__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L344)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L344?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.__iter__" class="docs-object-method">&nbsp;</a> 
```python
__iter__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L349)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L349?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L355)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L355?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.copy" class="docs-object-method">&nbsp;</a> 
```python
copy(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L358)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L358?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.prop" class="docs-object-method">&nbsp;</a> 
```python
prop(self, name, *args, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L374)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L374?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.mass_weighted_coords" class="docs-object-method">&nbsp;</a> 
```python
@property
mass_weighted_coords(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>


- `:returns`: `CoordinateSet`
    >

<a id="Psience.Molecools.Molecule.Molecule.center_of_mass" class="docs-object-method">&nbsp;</a> 
```python
@property
center_of_mass(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>


- `:returns`: `CoordinateSet`
    >

<a id="Psience.Molecools.Molecule.Molecule.inertia_tensor" class="docs-object-method">&nbsp;</a> 
```python
@property
inertia_tensor(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>


- `:returns`: `(np.ndarray, np.ndarray)`
    >

<a id="Psience.Molecools.Molecule.Molecule.inertial_eigensystem" class="docs-object-method">&nbsp;</a> 
```python
@property
inertial_eigensystem(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>


- `:returns`: `(np.ndarray, np.ndarray)`
    >

<a id="Psience.Molecools.Molecule.Molecule.moments_of_inertia" class="docs-object-method">&nbsp;</a> 
```python
@property
moments_of_inertia(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>


- `:returns`: `np.ndarray`
    >

<a id="Psience.Molecools.Molecule.Molecule.inertial_axes" class="docs-object-method">&nbsp;</a> 
```python
@property
inertial_axes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>


- `:returns`: `np.ndarray`
    >

<a id="Psience.Molecools.Molecule.Molecule.translation_rotation_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
translation_rotation_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>


- `:returns`: `np.ndarray`
    >

<a id="Psience.Molecools.Molecule.Molecule.canonicalize_internal_coordinate_spec" class="docs-object-method">&nbsp;</a> 
```python
canonicalize_internal_coordinate_spec(spec): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L456)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L456?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.internals" class="docs-object-method">&nbsp;</a> 
```python
@property
internals(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.zmatrix" class="docs-object-method">&nbsp;</a> 
```python
@property
zmatrix(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.internal_coordinates" class="docs-object-method">&nbsp;</a> 
```python
@property
internal_coordinates(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.g_matrix" class="docs-object-method">&nbsp;</a> 
```python
@property
g_matrix(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>

Returns the molecular g-matrix for the system
- `:returns`: `_`
    >

<a id="Psience.Molecools.Molecule.Molecule.bond_length" class="docs-object-method">&nbsp;</a> 
```python
bond_length(self, i, j): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L588)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L588?message=Update%20Docs)]
</div>

Returns the bond length of the coordinates
- `:returns`: `_`
    >
- `j`: `Any`
    >
- `i`: `Any`
    >

<a id="Psience.Molecools.Molecule.Molecule.bond_angle" class="docs-object-method">&nbsp;</a> 
```python
bond_angle(self, i, j, k): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L600)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L600?message=Update%20Docs)]
</div>

Returns the bond angle of the specified coordinates
- `:returns`: `_`
    >
- `j`: `Any`
    >
- `i`: `Any`
    >

<a id="Psience.Molecools.Molecule.Molecule.dihedral" class="docs-object-method">&nbsp;</a> 
```python
dihedral(self, i, j, k, l): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L612)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L612?message=Update%20Docs)]
</div>

Returns the dihedral angle of the specified coordinates
- `:returns`: `_`
    >
- `j`: `Any`
    >
- `i`: `Any`
    >

<a id="Psience.Molecools.Molecule.Molecule.principle_axis_frame" class="docs-object-method">&nbsp;</a> 
```python
principle_axis_frame(self, sel=None, inverse=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L628)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L628?message=Update%20Docs)]
</div>

Gets the principle axis frame(s) for the molecule
- `:returns`: `MolecularTransformation | List[MolecularTransformation]`
    >
- `inverse`: `bool`
    >whether to return the inverse of the rotations or not
- `sel`: `Any`
    >selection of atoms to use when getting the Eckart frame
- `mol`: `Any`
    >

<a id="Psience.Molecools.Molecule.Molecule.principle_axis_data" class="docs-object-method">&nbsp;</a> 
```python
@property
principle_axis_data(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>

Gets the principle axis embedded coords and embedding parameters for the molecule
- `:returns`: `MolecularTransformation | List[MolecularTransformation]`
    >

<a id="Psience.Molecools.Molecule.Molecule.eckart_frame" class="docs-object-method">&nbsp;</a> 
```python
eckart_frame(self, mol, sel=None, inverse=False, planar_ref_tolerance=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L652)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L652?message=Update%20Docs)]
</div>

Gets the Eckart frame(s) for the molecule
- `:returns`: `MolecularTransformation`
    >
- `inverse`: `bool`
    >whether to return the inverse of the rotations or not
- `sel`: `Any`
    >selection of atoms to use when getting the Eckart frame
- `mol`: `Any`
    >

<a id="Psience.Molecools.Molecule.Molecule.embed_coords" class="docs-object-method">&nbsp;</a> 
```python
embed_coords(self, crds, sel=None, planar_ref_tolerance=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L666)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L666?message=Update%20Docs)]
</div>

Embeds coords in the Eckart frame using `self` as a reference
- `:returns`: `_`
    >
- `crds`: `Any`
    >

<a id="Psience.Molecools.Molecule.Molecule.get_embedding_data" class="docs-object-method">&nbsp;</a> 
```python
get_embedding_data(self, crds, sel=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L676)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L676?message=Update%20Docs)]
</div>

Gets the necessary data to embed crds in the Eckart frame using `self` as a reference
- `:returns`: `tuple[np.ndarray, tuple[np.ndarray], tuple[np.ndarray]]`
    >
- `crds`: `Any`
    >

<a id="Psience.Molecools.Molecule.Molecule.get_embedded_molecule" class="docs-object-method">&nbsp;</a> 
```python
get_embedded_molecule(self, ref=None, embed_properties=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L685)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L685?message=Update%20Docs)]
</div>

Returns a Molecule embedded in an Eckart frame if ref is not None, otherwise returns
a principle-axis embedded Molecule
- `:returns`: `Molecule`
    >

<a id="Psience.Molecools.Molecule.Molecule.from_zmat" class="docs-object-method">&nbsp;</a> 
```python
from_zmat(zmat, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L718)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L718?message=Update%20Docs)]
</div>

Little z-matrix importer
- `:returns`: `Molecule`
    >
- `zmat`: `str | tuple`
    >

<a id="Psience.Molecools.Molecule.Molecule.from_pybel" class="docs-object-method">&nbsp;</a> 
```python
from_pybel(mol, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L735)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L735?message=Update%20Docs)]
</div>


- `:returns`: `_`
    >
- `mol`: `pybel.mol`
    >

<a id="Psience.Molecools.Molecule.Molecule.from_file" class="docs-object-method">&nbsp;</a> 
```python
from_file(file, mode=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L783)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L783?message=Update%20Docs)]
</div>

In general we'll delegate to pybel except for like Fchk and Log files
- `:returns`: `_`
    >
- `file`: `Any`
    >

<a id="Psience.Molecools.Molecule.Molecule.from_spec" class="docs-object-method">&nbsp;</a> 
```python
from_spec(spec): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L836)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L836?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.plot" class="docs-object-method">&nbsp;</a> 
```python
plot(self, *geometries, figure=None, bond_radius=0.1, atom_radius_scaling=0.25, atom_style=None, bond_style=None, mode='fast', objects=False, **plot_ops): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L870)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L870?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.jupyter_viz" class="docs-object-method">&nbsp;</a> 
```python
jupyter_viz(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L973)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L973?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.to_widget" class="docs-object-method">&nbsp;</a> 
```python
to_widget(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L980)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L980?message=Update%20Docs)]
</div>

 </div>
</div>




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Molecools/Molecule/Molecule.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Molecools/Molecule/Molecule.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Molecools/Molecule/Molecule.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Molecools/Molecule/Molecule.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L30?message=Update%20Docs)