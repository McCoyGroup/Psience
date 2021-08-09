## <a id="Psience.Molecools.Molecule.Molecule">Molecule</a>
General purpose 'Molecule' class where the 'Molecule' need not be a molecule at all

### Properties and Methods
<a id="Psience.Molecools.Molecule.Molecule.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, atoms, coords, bonds=None, masses=None, name=None, zmatrix=None, obmol=None, dipole_surface=None, dipole_derivatives=None, potential_surface=None, potential_derivatives=None, normal_modes=None, source_file=None, guess_bonds=True, charge=None, **metadata): 
```

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

- `:returns`: `DipoleSurfaceManager`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.dipole_derivatives" class="docs-object-method">&nbsp;</a>
```python
@property
dipole_derivatives(self): 
```

<a id="Psience.Molecools.Molecule.Molecule.potential_surface" class="docs-object-method">&nbsp;</a>
```python
@property
potential_surface(self): 
```

- `:returns`: `PotentialSurfaceManager`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.potential_derivatives" class="docs-object-method">&nbsp;</a>
```python
@property
potential_derivatives(self): 
```

<a id="Psience.Molecools.Molecule.Molecule.normal_modes" class="docs-object-method">&nbsp;</a>
```python
@property
normal_modes(self): 
```

- `:returns`: `NormalModesManager`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.metadata" class="docs-object-method">&nbsp;</a>
```python
@property
metadata(self): 
```

<a id="Psience.Molecools.Molecule.Molecule.__repr__" class="docs-object-method">&nbsp;</a>
```python
__repr__(self): 
```

<a id="Psience.Molecools.Molecule.Molecule.num_atoms" class="docs-object-method">&nbsp;</a>
```python
@property
num_atoms(self): 
```

<a id="Psience.Molecools.Molecule.Molecule.atom_positions" class="docs-object-method">&nbsp;</a>
```python
@property
atom_positions(self): 
```
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

<a id="Psience.Molecools.Molecule.Molecule.atoms" class="docs-object-method">&nbsp;</a>
```python
@property
atoms(self): 
```

<a id="Psience.Molecools.Molecule.Molecule.masses" class="docs-object-method">&nbsp;</a>
```python
@property
masses(self): 
```

<a id="Psience.Molecools.Molecule.Molecule.bonds" class="docs-object-method">&nbsp;</a>
```python
@property
bonds(self): 
```

<a id="Psience.Molecools.Molecule.Molecule.coords" class="docs-object-method">&nbsp;</a>
```python
@property
coords(self): 
```

<a id="Psience.Molecools.Molecule.Molecule.sys" class="docs-object-method">&nbsp;</a>
```python
@property
sys(self): 
```

<a id="Psience.Molecools.Molecule.Molecule.formula" class="docs-object-method">&nbsp;</a>
```python
@property
formula(self): 
```

<a id="Psience.Molecools.Molecule.Molecule.multiconfig" class="docs-object-method">&nbsp;</a>
```python
@property
multiconfig(self): 
```

<a id="Psience.Molecools.Molecule.Molecule.name" class="docs-object-method">&nbsp;</a>
```python
@property
name(self): 
```

<a id="Psience.Molecools.Molecule.Molecule.source_file" class="docs-object-method">&nbsp;</a>
```python
@property
source_file(self): 
```

<a id="Psience.Molecools.Molecule.Molecule.insert_atoms" class="docs-object-method">&nbsp;</a>
```python
insert_atoms(self, atoms, coords, where, handle_properties=True): 
```

<a id="Psience.Molecools.Molecule.Molecule.delete_atoms" class="docs-object-method">&nbsp;</a>
```python
delete_atoms(self, where, handle_properties=True): 
```

<a id="Psience.Molecools.Molecule.Molecule.take_submolecule" class="docs-object-method">&nbsp;</a>
```python
take_submolecule(self, spec): 
```
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

<a id="Psience.Molecools.Molecule.Molecule.__len__" class="docs-object-method">&nbsp;</a>
```python
__len__(self): 
```

<a id="Psience.Molecools.Molecule.Molecule.__iter__" class="docs-object-method">&nbsp;</a>
```python
__iter__(self): 
```

<a id="Psience.Molecools.Molecule.Molecule.__getitem__" class="docs-object-method">&nbsp;</a>
```python
__getitem__(self, item): 
```

<a id="Psience.Molecools.Molecule.Molecule.copy" class="docs-object-method">&nbsp;</a>
```python
copy(self): 
```

<a id="Psience.Molecools.Molecule.Molecule.prop" class="docs-object-method">&nbsp;</a>
```python
prop(self, name, *args, **kwargs): 
```

<a id="Psience.Molecools.Molecule.Molecule.mass_weighted_coords" class="docs-object-method">&nbsp;</a>
```python
@property
mass_weighted_coords(self): 
```

- `:returns`: `CoordinateSet`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.center_of_mass" class="docs-object-method">&nbsp;</a>
```python
@property
center_of_mass(self): 
```

- `:returns`: `CoordinateSet`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.inertia_tensor" class="docs-object-method">&nbsp;</a>
```python
@property
inertia_tensor(self): 
```

- `:returns`: `(np.ndarray, np.ndarray)`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.inertial_eigensystem" class="docs-object-method">&nbsp;</a>
```python
@property
inertial_eigensystem(self): 
```

- `:returns`: `(np.ndarray, np.ndarray)`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.moments_of_inertia" class="docs-object-method">&nbsp;</a>
```python
@property
moments_of_inertia(self): 
```

- `:returns`: `np.ndarray`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.inertial_axes" class="docs-object-method">&nbsp;</a>
```python
@property
inertial_axes(self): 
```

- `:returns`: `np.ndarray`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.translation_rotation_modes" class="docs-object-method">&nbsp;</a>
```python
@property
translation_rotation_modes(self): 
```

- `:returns`: `np.ndarray`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.zmatrix" class="docs-object-method">&nbsp;</a>
```python
@property
zmatrix(self): 
```

- `:returns`: `_`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.internal_coordinates" class="docs-object-method">&nbsp;</a>
```python
@property
internal_coordinates(self): 
```

<a id="Psience.Molecools.Molecule.Molecule.g_matrix" class="docs-object-method">&nbsp;</a>
```python
@property
g_matrix(self): 
```
Returns the molecular g-matrix for the system
- `:returns`: `_`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.bond_length" class="docs-object-method">&nbsp;</a>
```python
bond_length(self, i, j): 
```
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
Gets the principle axis embedded coords and embedding parameters for the molecule
- `:returns`: `MolecularTransformation | List[MolecularTransformation]`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.eckart_frame" class="docs-object-method">&nbsp;</a>
```python
eckart_frame(self, mol, sel=None, inverse=False): 
```
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
embed_coords(self, crds, sel=None): 
```
Embeds coords in the Eckart frame using `self` as a reference
- `crds`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.get_embedding_data" class="docs-object-method">&nbsp;</a>
```python
get_embedding_data(self, crds, sel=None): 
```
Gets the necessary data to embed crds in the Eckart frame using `self` as a reference
- `crds`: `Any`
    >No description...
- `:returns`: `tuple[np.ndarray, tuple[np.ndarray], tuple[np.ndarray]]`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.get_embedded_molecule" class="docs-object-method">&nbsp;</a>
```python
get_embedded_molecule(self, ref=None, embed_properties=True): 
```
Returns a Molecule embedded in an Eckart frame if ref is not None, otherwise returns
        a principle-axis embedded Molecule
- `:returns`: `Molecule`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.from_zmat" class="docs-object-method">&nbsp;</a>
```python
from_zmat(zmat, **opts): 
```
Little z-matrix importer
- `zmat`: `str | tuple`
    >No description...
- `:returns`: `Molecule`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.from_pybel" class="docs-object-method">&nbsp;</a>
```python
from_pybel(mol, **opts): 
```

- `mol`: `pybel.mol`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.from_file" class="docs-object-method">&nbsp;</a>
```python
from_file(file, mode=None, **opts): 
```
In general we'll delegate to pybel except for like Fchk and Log files
- `file`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.Molecools.Molecule.Molecule.from_spec" class="docs-object-method">&nbsp;</a>
```python
from_spec(spec): 
```

<a id="Psience.Molecools.Molecule.Molecule.plot" class="docs-object-method">&nbsp;</a>
```python
plot(self, *geometries, figure=None, bond_radius=0.1, atom_radius_scaling=0.25, atom_style=None, bond_style=None, mode='fast', objects=False, **plot_ops): 
```

### Examples




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/edit/ci/examples/ci/docs/Psience/Molecools/Molecule/Molecule.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/examples/ci/docs/Psience/Molecools/Molecule/Molecule.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/edit/ci/docs/ci/docs/Psience/Molecools/Molecule/Molecule.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/docs/templates/ci/docs/Psience/Molecools/Molecule/Molecule.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Molecools/Molecule.py?message=Update%20Docs)