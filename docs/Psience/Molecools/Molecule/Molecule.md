## <a id="Psience.Molecools.Molecule.Molecule">Molecule</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L25)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L25?message=Update%20Docs)]
</div>

General purpose 'Molecule' class where the 'Molecule' need not be a molecule at all

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

<a id="Psience.Molecools.Molecule.Molecule.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, atoms, coords, bonds=None, masses=None, name=None, zmatrix=None, obmol=None, dipole_surface=None, dipole_derivatives=None, potential_surface=None, potential_derivatives=None, normal_modes=None, source_file=None, guess_bonds=True, charge=None, **metadata): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L30)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L30?message=Update%20Docs)]
</div>


- `atoms`: `Iterable[str]`
    >atoms specified by name, either full name or short
- `coords`: `np.ndarray`
    >coordinates for the molecule, assumed to be in Bohr by default
- `bonds`: `Iterable[Iterable[int]] | None`
    >bond specification for the molecule
- `obmol`: `Any`
    >OpenBabel molecule for doing conversions
- `charge`: `int | None`
    >Net charge on the molecule
- `name`: `str | None`
    >Name for the molecule
- `dipole_surface`: `DipoleSurface | None`
    >The dipole surface for the system
- `dipole_derivatives`: `Iterable[np.ndarray] | None`
    >Derivatives of the dipole surface
- `potential_surface`: `PotentialSurface | None`
    >The potential surface for the system
- `potential_derivatives`: `Iterable[np.ndarray] | None`
    >Derivatives of the potential surface
- `guess_bonds`: `bool`
    >Whether or not to guess the bonding arrangement when that would be used
- `source_file`: `str`
    >The data file the molecule was loaded from
- `kw`: `Any`
    >Other bound parameters that might be useful

<a id="Psience.Molecools.Molecule.Molecule.dipole_surface" class="docs-object-method">&nbsp;</a> 
```python
@property
dipole_surface(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>


- `:returns`: `DipoleSurfaceManager`
    >No description...

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
    >No description...

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
    >No description...

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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L185)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L185?message=Update%20Docs)]
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
- `spec`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L264)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L264?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.delete_atoms" class="docs-object-method">&nbsp;</a> 
```python
delete_atoms(self, where, handle_properties=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L290)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L290?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.take_submolecule" class="docs-object-method">&nbsp;</a> 
```python
take_submolecule(self, spec): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L306)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L306?message=Update%20Docs)]
</div>

Takes a 'slice' of a molecule if working with Cartesian coords.
        If not, need to do some corner case handling for that.
- `spec`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L334)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L334?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.__iter__" class="docs-object-method">&nbsp;</a> 
```python
__iter__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L339)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L339?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L345)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L345?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.copy" class="docs-object-method">&nbsp;</a> 
```python
copy(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L348)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L348?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.prop" class="docs-object-method">&nbsp;</a> 
```python
prop(self, name, *args, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L364)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L364?message=Update%20Docs)]
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
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.center_of_mass" class="docs-object-method">&nbsp;</a> 
```python
@property
center_of_mass(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>


- `:returns`: `CoordinateSet`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.inertia_tensor" class="docs-object-method">&nbsp;</a> 
```python
@property
inertia_tensor(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>


- `:returns`: `(np.ndarray, np.ndarray)`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.inertial_eigensystem" class="docs-object-method">&nbsp;</a> 
```python
@property
inertial_eigensystem(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>


- `:returns`: `(np.ndarray, np.ndarray)`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.moments_of_inertia" class="docs-object-method">&nbsp;</a> 
```python
@property
moments_of_inertia(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>


- `:returns`: `np.ndarray`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.inertial_axes" class="docs-object-method">&nbsp;</a> 
```python
@property
inertial_axes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>


- `:returns`: `np.ndarray`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.translation_rotation_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
translation_rotation_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>


- `:returns`: `np.ndarray`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.zmatrix" class="docs-object-method">&nbsp;</a> 
```python
@property
zmatrix(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L?message=Update%20Docs)]
</div>


- `:returns`: `_`
    >No description...

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
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.bond_length" class="docs-object-method">&nbsp;</a> 
```python
bond_length(self, i, j): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L466)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L466?message=Update%20Docs)]
</div>

Returns the bond length of the coordinates
- `i`: `Any`
    >No description...
- `j`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.bond_angle" class="docs-object-method">&nbsp;</a> 
```python
bond_angle(self, i, j, k): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L478)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L478?message=Update%20Docs)]
</div>

Returns the bond angle of the specified coordinates
- `i`: `Any`
    >No description...
- `j`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.dihedral" class="docs-object-method">&nbsp;</a> 
```python
dihedral(self, i, j, k, l): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L490)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L490?message=Update%20Docs)]
</div>

Returns the dihedral angle of the specified coordinates
- `i`: `Any`
    >No description...
- `j`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.principle_axis_frame" class="docs-object-method">&nbsp;</a> 
```python
principle_axis_frame(self, sel=None, inverse=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L506)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L506?message=Update%20Docs)]
</div>

Gets the principle axis frame(s) for the molecule
- `mol`: `Any`
    >No description...
- `sel`: `Any`
    >selection of atoms to use when getting the Eckart frame
- `inverse`: `bool`
    >whether to return the inverse of the rotations or not
- `:returns`: `MolecularTransformation | List[MolecularTransformation]`
    >No description...

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
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.eckart_frame" class="docs-object-method">&nbsp;</a> 
```python
eckart_frame(self, mol, sel=None, inverse=False, planar_ref_tolerance=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L530)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L530?message=Update%20Docs)]
</div>

Gets the Eckart frame(s) for the molecule
- `mol`: `Any`
    >No description...
- `sel`: `Any`
    >selection of atoms to use when getting the Eckart frame
- `inverse`: `bool`
    >whether to return the inverse of the rotations or not
- `:returns`: `MolecularTransformation`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.embed_coords" class="docs-object-method">&nbsp;</a> 
```python
embed_coords(self, crds, sel=None, planar_ref_tolerance=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L544)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L544?message=Update%20Docs)]
</div>

Embeds coords in the Eckart frame using `self` as a reference
- `crds`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.get_embedding_data" class="docs-object-method">&nbsp;</a> 
```python
get_embedding_data(self, crds, sel=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L554)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L554?message=Update%20Docs)]
</div>

Gets the necessary data to embed crds in the Eckart frame using `self` as a reference
- `crds`: `Any`
    >No description...
- `:returns`: `tuple[np.ndarray, tuple[np.ndarray], tuple[np.ndarray]]`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.get_embedded_molecule" class="docs-object-method">&nbsp;</a> 
```python
get_embedded_molecule(self, ref=None, embed_properties=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L563)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L563?message=Update%20Docs)]
</div>

Returns a Molecule embedded in an Eckart frame if ref is not None, otherwise returns
        a principle-axis embedded Molecule
- `:returns`: `Molecule`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.from_zmat" class="docs-object-method">&nbsp;</a> 
```python
from_zmat(zmat, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L596)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L596?message=Update%20Docs)]
</div>

Little z-matrix importer
- `zmat`: `str | tuple`
    >No description...
- `:returns`: `Molecule`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.from_pybel" class="docs-object-method">&nbsp;</a> 
```python
from_pybel(mol, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L613)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L613?message=Update%20Docs)]
</div>


- `mol`: `pybel.mol`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.from_file" class="docs-object-method">&nbsp;</a> 
```python
from_file(file, mode=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L661)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L661?message=Update%20Docs)]
</div>

In general we'll delegate to pybel except for like Fchk and Log files
- `file`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.from_spec" class="docs-object-method">&nbsp;</a> 
```python
from_spec(spec): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L714)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L714?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.plot" class="docs-object-method">&nbsp;</a> 
```python
plot(self, *geometries, figure=None, bond_radius=0.1, atom_radius_scaling=0.25, atom_style=None, bond_style=None, mode='fast', objects=False, **plot_ops): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L751)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L751?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.Molecule.Molecule.jupyter_viz" class="docs-object-method">&nbsp;</a> 
```python
jupyter_viz(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L849)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L849?message=Update%20Docs)]
</div>

 </div>
</div>




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Molecools/Molecule/Molecule.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Molecools/Molecule/Molecule.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Molecools/Molecule/Molecule.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Molecools/Molecule/Molecule.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L25?message=Update%20Docs)