## <a id="Psience.Molecools.Molecule.Molecule">Molecule</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L51)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L51?message=Update%20Docs)]
</div>

General purpose 'Molecule' class where the 'Molecule' need not be a molecule at all







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
default_coordinate_pruning: str
bond_guessing_mode: str
default_energy_evalutor: str
default_display_mode: str
```
<a id="Psience.Molecools.Molecule.Molecule.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, atoms, coords, bonds=None, masses=None, name=None, internals=None, rdmol=None, dipole_surface=None, dipole_derivatives=None, potential_surface=None, potential_derivatives=None, normal_modes=None, source_file=None, guess_bonds=True, charge=None, formal_charges=None, spin=None, display_mode=None, display_settings=None, energy=None, energy_evaluator=None, dipole_evaluator=None, charge_evaluator=None, polarizability_evaluator=None, polarizability_derivatives=None, checkpoint_file=None, **metadata): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule.py#L56)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L56?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L184)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L184?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a new `Molecule` that is a copy of this one with the given fields overridden, treating any argument left at its `dev.default` sentinel as "keep the current value" (with some fields, like masses/energy/surfaces, only carried over automatically when the arguments they depend on -- e.g. `atoms`, `coords`, `energy_evaluator` -- are also left unspecified).
  - `atoms`: `Iterable[str] | object`
    > replacement atoms, or `dev.default` to keep the current ones
  - `coords`: `np.ndarray | object`
    > replacement coordinates, or `dev.default` to keep the current ones
  - `internals`: `object`
    > replacement internal-coordinate specification
  - `masses`: `np.ndarray | object`
    > replacement masses
  - `bonds`: `object`
    > replacement bonds
  - `guess_bonds`: `bool | object`
    > replacement bond-guessing flag
  - `energy`: `float | object`
    > replacement cached energy value
  - `energy_evaluator`: `object`
    > replacement energy evaluator
  - `dipole_evaluator`: `object`
    > replacement dipole evaluator
  - `charge_evaluator`: `object`
    > replacement charge evaluator
  - `polarizability_evaluator`: `object`
    > replacement polarizability evaluator
  - `charge`: `int | object`
    > replacement net charge
  - `spin`: `object`
    > replacement spin
  - `rdmol`: `object`
    > replacement RDKit molecule
  - `display_mode`: `str | object`
    > replacement display mode
  - `display_settings`: `dict | object`
    > replacement display settings
  - `normal_modes`: `object`
    > replacement normal modes
  - `dipole_surface`: `object`
    > replacement dipole surface
  - `potential_surface`: `object`
    > replacement potential surface
  - `dipole_derivatives`: `object`
    > replacement dipole derivatives
  - `potential_derivatives`: `object`
    > replacement potential derivatives
  - `polarizability_derivatives`: `object`
    > replacement polarizability derivatives
  - `meta`: `dict | object`
    > replacement/merged metadata
  - `source_file`: `str | object`
    > replacement source file path
  - `:returns`: `Molecule`
    > the new, modified `Molecule`


<a id="Psience.Molecools.Molecule.Molecule.__del__" class="docs-object-method">&nbsp;</a> 
```python
__del__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L349)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L349?message=Update%20Docs)]
</div>
**LLM Docstring**

Clean up the molecule's coordinate embedding (deregistering any converters it registered) when the object is garbage collected.
  - `:returns`: `None`
    > None


<a id="Psience.Molecools.Molecule.Molecule.to_state" class="docs-object-method">&nbsp;</a> 
```python
to_state(self, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L361)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L361?message=Update%20Docs)]
</div>
**LLM Docstring**

Serialize this molecule's essential data (atoms, masses, coordinates, bonds, internal-coordinate spec, evaluators, potential/dipole derivatives, charge, spin) into a plain dict, stripping out non-serializable embedding-specific converter options first.
  - `serializer`: `object | None`
    > accepted for interface consistency but not used in this method's body
  - `:returns`: `dict`
    > the serialized state dict


<a id="Psience.Molecools.Molecule.Molecule.from_state" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_state(cls, data, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L419)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L419?message=Update%20Docs)]
</div>
**LLM Docstring**

Reconstruct a `Molecule` from a previously serialized state dict, by passing its entries directly as constructor keyword arguments.
  - `data`: `dict`
    > the serialized state, as produced by `to_state`
  - `serializer`: `object | None`
    > accepted for interface consistency but not used in this method's body
  - `:returns`: `Molecule`
    > the reconstructed molecule


<a id="Psience.Molecools.Molecule.Molecule.cached_eval" class="docs-object-method">&nbsp;</a> 
```python
cached_eval(self, key, generator, *, condition=None, args=(), kwargs=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L435)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L435?message=Update%20Docs)]
</div>
**LLM Docstring**

Evaluate (and cache, via this molecule's on-disk/in-memory checkpoint) a value under `key`, delegating to `self.checkpoint.cached_eval`.
  - `key`: `str`
    > the cache key to look up or populate
  - `generator`: `callable`
    > callable used to compute the value when it is not already cached
  - `condition`: `callable | None`
    > optional predicate controlling whether the cached value should be recomputed
  - `args`: `tuple`
    > positional arguments passed to `generator` if it is called
  - `kwargs`: `dict | None`
    > keyword arguments passed to `generator` if it is called
  - `:returns`: `object`
    > the cached or newly computed value


<a id="Psience.Molecools.Molecule.Molecule.canonicalize_internals" class="docs-object-method">&nbsp;</a> 
```python
canonicalize_internals(self, spec, atoms, coords, bonds, relocalize=True, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L629)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L629?message=Update%20Docs)]
</div>
**LLM Docstring**

Normalize the many accepted forms of an internal-coordinate specification (the strings `'auto'`/`'zmatrix'`, a dict with `'primitives'`/`'specs'`/`'zmatrix'` keys where `'specs'` may itself be `'auto'`/`'natural'`, a bare Z-matrix-like array, or a bare list of primitive specs) down into the canonical dict form expected by `MolecularEmbedding`, recursively re-dispatching as needed.
  - `spec`: `str | dict | Iterable | None`
    > the internal-coordinate specification to canonicalize
  - `atoms`: `Iterable[str]`
    > the atom labels
  - `coords`: `np.ndarray`
    > the Cartesian coordinates
  - `bonds`: `Iterable[Iterable[int]] | None`
    > the bonds to use when auto-generating coordinates
  - `relocalize`: `bool`
    > whether redundant coordinates should be relocalized by default
  - `masses`: `np.ndarray | None`
    > atomic masses, forwarded to the auto-generation routines
  - `:returns`: `dict | None`
    > the canonicalized specification


<a id="Psience.Molecools.Molecule.Molecule.prep_internal_spec" class="docs-object-method">&nbsp;</a> 
```python
prep_internal_spec(self, spec, relocalize=True, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L707)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L707?message=Update%20Docs)]
</div>
**LLM Docstring**

Canonicalize an internal-coordinate specification against this molecule's own atoms, coordinates, bonds, and masses, via `canonicalize_internals`.
  - `spec`: `object`
    > the internal-coordinate specification to canonicalize
  - `relocalize`: `bool`
    > whether redundant coordinates should be relocalized by default
  - `masses`: `np.ndarray | None`
    > atomic masses to use instead of this molecule's own
  - `:returns`: `dict | None`
    > the canonicalized specification


<a id="Psience.Molecools.Molecule.Molecule.embedding" class="docs-object-method">&nbsp;</a> 
```python
@property
embedding(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L731)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L731?message=Update%20Docs)]
</div>
**LLM Docstring**

Property getter/setter for the molecule's `MolecularEmbedding`. Setting it also resets the cached evaluator and Hamiltonian, since both depend on the embedding.
  - `e`: `MolecularEmbedding`
    > (setter only) the new embedding
  - `:returns`: `MolecularEmbedding`
    > (getter) the current embedding


<a id="Psience.Molecools.Molecule.Molecule.get_evaluator" class="docs-object-method">&nbsp;</a> 
```python
get_evaluator(self, embedding=None, normal_modes=<McUtils.Devutils.core.DefaultType instance>): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L759)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L759?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a `MolecularEvaluator` for this molecule's embedding (or an alternate one), using either the given `normal_modes` or this molecule's own.
  - `embedding`: `MolecularEmbedding | None`
    > an alternate embedding to build the evaluator for; defaults to `self.embedding`
  - `normal_modes`: `object`
    > normal modes to use instead of `self._normal_modes`
  - `:returns`: `MolecularEvaluator`
    > the constructed evaluator


<a id="Psience.Molecools.Molecule.Molecule.evaluator" class="docs-object-method">&nbsp;</a> 
```python
@property
evaluator(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L775)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L775?message=Update%20Docs)]
</div>
**LLM Docstring**

Property getter/setter for the `MolecularEvaluator` used for energy/derivative calculations. The getter lazily builds one via `get_evaluator` the first time it's needed.
  - `e`: `MolecularEvaluator`
    > (setter only) the new evaluator
  - `:returns`: `MolecularEvaluator`
    > (getter) the cached (or newly built) evaluator


<a id="Psience.Molecools.Molecule.Molecule.coords" class="docs-object-method">&nbsp;</a> 
```python
@property
coords(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L804)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L804?message=Update%20Docs)]
</div>
**LLM Docstring**

Property getter/setter for the Cartesian coordinates, delegating to `self.embedding.coords`.
  - `coords`: `np.ndarray`
    > (setter only) the new Cartesian coordinates
  - `:returns`: `CoordinateSet`
    > (getter) the Cartesian coordinates


<a id="Psience.Molecools.Molecule.Molecule.masses" class="docs-object-method">&nbsp;</a> 
```python
@property
masses(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L830)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L830?message=Update%20Docs)]
</div>
**LLM Docstring**

Property getter/setter for the atomic masses. The setter also updates the embedding's masses (via `atomic_masses`, i.e. in atomic units).
  - `masses`: `np.ndarray`
    > (setter only) the new masses
  - `:returns`: `np.ndarray`
    > (getter) the atomic masses


<a id="Psience.Molecools.Molecule.Molecule.internals" class="docs-object-method">&nbsp;</a> 
```python
@property
internals(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L857)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L857?message=Update%20Docs)]
</div>
**LLM Docstring**

Getter for the raw (canonicalized) internal-coordinate specification, delegating to `self.embedding.internals`. (A setter with the same name separately rebuilds the embedding from a new specification.)
  - `:returns`: `dict | None`
    > the internal-coordinate specification, or `None` if none is set


<a id="Psience.Molecools.Molecule.Molecule.charge" class="docs-object-method">&nbsp;</a> 
```python
@property
charge(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L868)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L868?message=Update%20Docs)]
</div>
**LLM Docstring**

Property getter/setter for the molecule's net charge, stored in its metadata (getter defaults to `0` if unset).
  - `c`: `int`
    > (setter only) the new net charge
  - `:returns`: `int`
    > (getter) the net charge


<a id="Psience.Molecools.Molecule.Molecule.spin" class="docs-object-method">&nbsp;</a> 
```python
@property
spin(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L894)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L894?message=Update%20Docs)]
</div>
**LLM Docstring**

Property getter/setter for the molecule's spin, stored in its metadata.
  - `c`: `object`
    > (setter only) the new spin value
  - `:returns`: `object | None`
    > (getter) the spin, or `None` if unset


<a id="Psience.Molecools.Molecule.Molecule.charges" class="docs-object-method">&nbsp;</a> 
```python
@property
charges(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L920)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L920?message=Update%20Docs)]
</div>
**LLM Docstring**

Property getter/setter for the molecule's per-atom partial charges, stored in its metadata.
  - `c`: `np.ndarray`
    > (setter only) the new per-atom charges
  - `:returns`: `np.ndarray | None`
    > (getter) the per-atom charges, or `None` if unset


<a id="Psience.Molecools.Molecule.Molecule.formal_charges" class="docs-object-method">&nbsp;</a> 
```python
@property
formal_charges(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L946)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L946?message=Update%20Docs)]
</div>
**LLM Docstring**

Property getter/setter for the molecule's per-atom formal charges, stored in its metadata.
  - `c`: `np.ndarray`
    > (setter only) the new per-atom formal charges
  - `:returns`: `np.ndarray | None`
    > (getter) the per-atom formal charges, or `None` if unset


<a id="Psience.Molecools.Molecule.Molecule.get_charge_evaluator" class="docs-object-method">&nbsp;</a> 
```python
get_charge_evaluator(self, evaluator=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L973)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L973?message=Update%20Docs)]
</div>
**LLM Docstring**

Resolve (and, if needed, instantiate from this molecule) a charge-evaluator object, defaulting to `self.charge_evaluator` if none is given explicitly.
  - `evaluator`: `object | None`
    > an explicit evaluator (or evaluator-type specification) to resolve; defaults to `self.charge_evaluator`
  - `opts`: `dict`
    > extra options forwarded to the evaluator's `from_mol` constructor, if applicable
  - `:returns`: `object`
    > the resolved charge-evaluator instance


<a id="Psience.Molecools.Molecule.Molecule.calculate_charges" class="docs-object-method">&nbsp;</a> 
```python
calculate_charges(self, evaluator=None, order=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L996)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L996?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute partial-charge values (and, optionally, their derivatives) using the resolved charge evaluator at this molecule's current geometry.
  - `evaluator`: `object | None`
    > an explicit evaluator (or evaluator-type specification) to use; defaults to `self.charge_evaluator`
  - `order`: `int | None`
    > the highest derivative order to compute; if `None`, only the charges themselves are returned
  - `opts`: `dict`
    > extra options forwarded to `get_charge_evaluator`
  - `:returns`: `np.ndarray | list[np.ndarray]`
    > the charges (if `order` is `None`) or the full charge/derivative expansion


<a id="Psience.Molecools.Molecule.Molecule.internal_coordinates" class="docs-object-method">&nbsp;</a> 
```python
@property
internal_coordinates(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1035)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1035?message=Update%20Docs)]
</div>
**LLM Docstring**

Getter for the internal coordinates at the molecule's current geometry, delegating to `self.embedding.internal_coordinates`.
  - `:returns`: `CoordinateSet | None`
    > the internal coordinates, or `None` if none are defined


<a id="Psience.Molecools.Molecule.Molecule.redundant_internal_transformation" class="docs-object-method">&nbsp;</a> 
```python
@property
redundant_internal_transformation(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1046)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1046?message=Update%20Docs)]
</div>
**LLM Docstring**

Getter for the redundant-to-non-redundant internal-coordinate transformation, delegating to `self.embedding.redundant_internal_transformation`.
  - `:returns`: `np.ndarray | None`
    > the redundant transformation, or `None` if not applicable


<a id="Psience.Molecools.Molecule.Molecule.get_coordinate_filer" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_coordinate_filer(cls, allowed_coordinate_types=None, excluded_coordinate_types=None, allowed_ring_types=None, excluded_ring_types=None, allowed_group_types=None, excluded_group_types=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1103)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1103?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a filter function (closing over the given allow/exclude criteria) that, given a dict of coordinate-to-label mappings, returns only the entries whose label passes `_check_label`.
  - `allowed_coordinate_types`: `Iterable | None`
    > forwarded to `_check_label`
  - `excluded_coordinate_types`: `Iterable | None`
    > forwarded to `_check_label`
  - `allowed_ring_types`: `Iterable | None`
    > forwarded to `_check_label`
  - `excluded_ring_types`: `Iterable | None`
    > forwarded to `_check_label`
  - `allowed_group_types`: `Iterable | None`
    > forwarded to `_check_label`
  - `excluded_group_types`: `Iterable | None`
    > forwarded to `_check_label`
  - `:returns`: `callable`
    > the constructed coordinate-filtering function


<a id="Psience.Molecools.Molecule.Molecule.get_bond_graph_internals" class="docs-object-method">&nbsp;</a> 
```python
get_bond_graph_internals(self, include_stretches=True, include_bends=True, include_dihedrals=True, include_fragments=True, pruning=None, fragment=None, base_internals=None, use_distance_matrix=True, concatenate=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1159)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1159?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a set of internal coordinates (bond stretches, bends, dihedrals, and/or inter-fragment coordinates) directly from the bonding graph, optionally restricted to a single fragment (recursively, with the result permuted back into the full atom indexing) and/or pruned down to a well-conditioned subset.
  - `include_stretches`: `bool`
    > whether to include bond-stretch coordinates
  - `include_bends`: `bool`
    > whether to include bond-angle coordinates
  - `include_dihedrals`: `bool`
    > whether to include dihedral-angle coordinates
  - `include_fragments`: `bool`
    > whether to include coordinates connecting separate molecular fragments
  - `pruning`: `bool | str | dict | None`
    > whether/how to prune the resulting coordinates (`True` for the default method, or an explicit method spec), forwarded to `prune_internals`
  - `fragment`: `int | Iterable[int] | None`
    > restrict to a single fragment, given as a fragment index or an explicit list of atom indices
  - `base_internals`: `object | None`
    > accepted and forwarded when recursing on a fragment, but not otherwise used directly in this method's own body
  - `use_distance_matrix`: `bool`
    > whether to precompute a distance matrix for the fragment-coordinate generation
  - `concatenate`: `bool`
    > whether to concatenate the different coordinate categories (stretches/bends/dihedrals/fragments) into a single list, or return them as separate groups
  - `:returns`: `list`
    > the generated internal coordinates, as a single concatenated list or a list of category groups depending on `concatenate`


<a id="Psience.Molecools.Molecule.Molecule.prune_internals" class="docs-object-method">&nbsp;</a> 
```python
prune_internals(self, coords, method='b_matrix', check_rigidity=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1256)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1256?message=Update%20Docs)]
</div>
**LLM Docstring**

Reduce a set of internal coordinates down to a non-redundant, well-conditioned subset, defaulting to a B-matrix-rank-based method (building the necessary translation/rotation-projected B-matrix generator and a sensible `max_coords` cap) if no custom method is supplied.
  - `coords`: `list`
    > the internal-coordinate specs to prune
  - `method`: `str | dict`
    > the pruning method: a method-name string, or a dict of method options (with a `'method'` key defaulting to `'b_matrix'`)
  - `check_rigidity`: `bool`
    > whether to check that the pruned coordinate set spans a rigid (non-redundant) representation
  - `:returns`: `list`
    > the pruned coordinate specs


<a id="Psience.Molecools.Molecule.Molecule.get_labeled_internals" class="docs-object-method">&nbsp;</a> 
```python
get_labeled_internals(self, coordinate_filter=None, allowed_coordinate_types=None, excluded_coordinate_types=None, allowed_ring_types=None, excluded_ring_types=None, allowed_group_types=None, excluded_group_types=None, include_stretches=True, include_bends=True, include_dihedrals=True, include_fragments=True, coordinate_sorting=None, pruning=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1307)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1307?message=Update%20Docs)]
</div>
**LLM Docstring**

Build the internal coordinates from the bonding graph (via `get_bond_graph_internals`) and label each one by its atom types/ring/functional-group membership (via `edge_graph.get_label_types` and `coordops.get_coordinate_label`), then filter and sort them.
  - `coordinate_filter`: `callable | None`
    > an explicit filter function to apply instead of building one from the allow/exclude arguments
  - `allowed_coordinate_types`: `Iterable | None`
    > forwarded to `get_coordinate_filer` if `coordinate_filter` is not given
  - `excluded_coordinate_types`: `Iterable | None`
    > forwarded to `get_coordinate_filer`
  - `allowed_ring_types`: `Iterable | None`
    > forwarded to `get_coordinate_filer`
  - `excluded_ring_types`: `Iterable | None`
    > forwarded to `get_coordinate_filer`
  - `allowed_group_types`: `Iterable | None`
    > forwarded to `get_coordinate_filer`
  - `excluded_group_types`: `Iterable | None`
    > forwarded to `get_coordinate_filer`
  - `include_stretches`: `bool`
    > whether to include bond-stretch coordinates
  - `include_bends`: `bool`
    > whether to include bond-angle coordinates
  - `include_dihedrals`: `bool`
    > whether to include dihedral-angle coordinates
  - `include_fragments`: `bool`
    > whether to include inter-fragment coordinates
  - `coordinate_sorting`: `callable | bool | None`
    > a custom sorting function to apply to the labeled coordinates instead of the default `coordops.sort_internal_coordinates`; pass a falsy value to skip sorting
  - `pruning`: `bool | str | dict`
    > whether/how to prune the coordinates, forwarded to `get_bond_graph_internals`
  - `:returns`: `dict`
    > a mapping from coordinate spec to its label, filtered and sorted


<a id="Psience.Molecools.Molecule.Molecule.get_mode_labels" class="docs-object-method">&nbsp;</a> 
```python
get_mode_labels(self, internals=None, modes=None, use_redundants=True, expansions=None, return_modes=False, **internals_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1394)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1394?message=Update%20Docs)]
</div>
**LLM Docstring**

Assign human-readable labels (e.g. "C-H stretch") to a set of normal modes by projecting them onto labeled internal coordinates, handling both redundant and non-redundant internal-coordinate expansions and both Cartesian- and internal-coordinate-basis modes.
  - `internals`: `dict | None`
    > the labeled internal coordinates to project onto; computed via `get_labeled_internals` if not given
  - `modes`: `object | None`
    > the normal modes to label; computed via `get_normal_modes` if not given
  - `use_redundants`: `bool`
    > whether to build a redundant-coordinate expansion (with relocalization) for the projection, rather than using the internal coordinates directly
  - `expansions`: `tuple | None`
    > precomputed `(expansions, inverse_expansion)` internal-coordinate Jacobian data to reuse instead of recomputing it
  - `return_modes`: `bool`
    > whether to also return the internal-coordinate-basis mode matrix alongside the labels
  - `internals_opts`: `dict`
    > extra options forwarded to `get_labeled_internals` if `internals` is not given
  - `:returns`: `list | tuple`
    > the mode labels, or `(internal_modes, labels)` if `return_modes` is set


<a id="Psience.Molecools.Molecule.Molecule.mode_embedding" class="docs-object-method">&nbsp;</a> 
```python
@property
mode_embedding(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1487)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1487?message=Update%20Docs)]
</div>
**LLM Docstring**

The (cached) `ModeEmbedding` combining this molecule's coordinate embedding with its normal modes, built lazily the first time it's needed.
  - `:returns`: `ModeEmbedding`
    > the mode embedding


<a id="Psience.Molecools.Molecule.Molecule.get_internals" class="docs-object-method">&nbsp;</a> 
```python
get_internals(self, coords=None, *, strip_embedding=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1500)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1500?message=Update%20Docs)]
</div>
**LLM Docstring**

Fetch internal coordinates, either the molecule's own cached ones or those for an alternate set of Cartesian `coords`, via `self.embedding.get_internals`.
  - `coords`: `np.ndarray | None`
    > alternate Cartesian coordinates to convert instead of using the cached internal coordinates
  - `strip_embedding`: `bool`
    > whether to strip the fixed embedding coordinates from the result
  - `:returns`: `CoordinateSet | None`
    > the internal coordinates, or `None` if none are defined


<a id="Psience.Molecools.Molecule.Molecule.get_cartesians_by_internals" class="docs-object-method">&nbsp;</a> 
```python
get_cartesians_by_internals(self, order=None, coords=None, *, strip_embedding=False, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1515)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1515?message=Update%20Docs)]
</div>
**LLM Docstring**

Fetch the Cartesians-by-internals Jacobian expansion, via `self.embedding.get_cartesians_by_internals`.
  - `order`: `int | None`
    > the highest derivative order to compute
  - `coords`: `np.ndarray | None`
    > alternate coordinates to compute the Jacobian at
  - `strip_embedding`: `bool`
    > whether to strip the fixed embedding coordinates from the result
  - `kw`: `dict`
    > extra options forwarded to the embedding
  - `:returns`: `list[np.ndarray]`
    > the Cartesians-by-internals Jacobian tensors


<a id="Psience.Molecools.Molecule.Molecule.get_internals_by_cartesians" class="docs-object-method">&nbsp;</a> 
```python
get_internals_by_cartesians(self, order=None, *, coords=None, strip_embedding=False, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1534)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1534?message=Update%20Docs)]
</div>
**LLM Docstring**

Fetch the internals-by-Cartesians Jacobian expansion, via `self.embedding.get_internals_by_cartesians`.
  - `order`: `int | None`
    > the highest derivative order to compute
  - `coords`: `np.ndarray | None`
    > alternate coordinates to compute the Jacobian at
  - `strip_embedding`: `bool`
    > whether to strip the fixed embedding coordinates from the result
  - `kw`: `dict`
    > extra options forwarded to the embedding
  - `:returns`: `list[np.ndarray]`
    > the internals-by-Cartesians Jacobian tensors


<a id="Psience.Molecools.Molecule.Molecule.get_cartesians_by_modes" class="docs-object-method">&nbsp;</a> 
```python
get_cartesians_by_modes(self, order=None, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1553)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1553?message=Update%20Docs)]
</div>
**LLM Docstring**

Fetch the Cartesians-by-normal-modes Jacobian expansion, via `self.mode_embedding.get_cartesians_by_internals`.
  - `order`: `int | None`
    > the highest derivative order to compute
  - `kw`: `dict`
    > extra options forwarded to the mode embedding
  - `:returns`: `list[np.ndarray]`
    > the Cartesians-by-modes Jacobian tensors


<a id="Psience.Molecools.Molecule.Molecule.get_modes_by_cartesians" class="docs-object-method">&nbsp;</a> 
```python
get_modes_by_cartesians(self, order=None, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1568)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1568?message=Update%20Docs)]
</div>
**LLM Docstring**

Fetch the normal-modes-by-Cartesians Jacobian expansion, via `self.mode_embedding.get_internals_by_cartesians`.
  - `order`: `int | None`
    > the highest derivative order to compute
  - `kw`: `dict`
    > extra options forwarded to the mode embedding
  - `:returns`: `list[np.ndarray]`
    > the modes-by-Cartesians Jacobian tensors


<a id="Psience.Molecools.Molecule.Molecule.dipole_surface" class="docs-object-method">&nbsp;</a> 
```python
@property
dipole_surface(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1586)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1586?message=Update%20Docs)]
</div>

  - `:returns`: `DipoleSurfaceManager`
    >


<a id="Psience.Molecools.Molecule.Molecule.dipole_derivatives" class="docs-object-method">&nbsp;</a> 
```python
@property
dipole_derivatives(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1611)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1611?message=Update%20Docs)]
</div>
**LLM Docstring**

Property getter/setter for the dipole derivative tensors. The getter delegates to `self.dipole_surface.get_derivatives(quiet=True)`; the setter assigns to `self.dipole_surface.derivatives`.
  - `derivs`: `list[np.ndarray]`
    > (setter only) the new dipole derivative tensors
  - `:returns`: `list[np.ndarray] | None`
    > (getter) the dipole derivative tensors, or `None` if unavailable


<a id="Psience.Molecools.Molecule.Molecule.get_cartesian_dipole_derivatives" class="docs-object-method">&nbsp;</a> 
```python
get_cartesian_dipole_derivatives(self, order=None, evaluator=None, include_constant_term=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1637)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1637?message=Update%20Docs)]
</div>
**LLM Docstring**

Fetch the dipole derivatives in Cartesian coordinates, computing them via `calculate_dipole` (and caching the result on the molecule, if the same evaluator is configured as the default) if not already available to the requested order.
  - `order`: `int | None`
    > the highest derivative order needed; if `None`, whatever is available is returned
  - `evaluator`: `object | None`
    > an explicit dipole evaluator to use instead of `self.dipole_evaluator`
  - `include_constant_term`: `bool`
    > whether to include the zeroth-order (reference dipole) term in the result
  - `:returns`: `list[np.ndarray] | None`
    > the dipole derivative tensors (from first order, or zeroth if `include_constant_term`), or `None` if unavailable


<a id="Psience.Molecools.Molecule.Molecule.get_internal_dipole_derivatives" class="docs-object-method">&nbsp;</a> 
```python
get_internal_dipole_derivatives(self, order=None, reembed=True, strip_embedding=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1677)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1677?message=Update%20Docs)]
</div>
**LLM Docstring**

Fetch the dipole derivatives re-expressed in internal coordinates, by re-expanding the Cartesian dipole derivatives through the Cartesians-by-internals Jacobian.
  - `order`: `int | None`
    > the highest derivative order needed
  - `reembed`: `bool`
    > whether to use the Eckart-reembedded Cartesians-by-internals Jacobian
  - `strip_embedding`: `bool`
    > whether to strip the fixed embedding coordinates from the Jacobian
  - `:returns`: `list[np.ndarray]`
    > the internal-coordinate dipole derivative tensors


<a id="Psience.Molecools.Molecule.Molecule.get_cartesian_polarizability_derivatives" class="docs-object-method">&nbsp;</a> 
```python
get_cartesian_polarizability_derivatives(self, order=None, evaluator=None, include_constant_term=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1700)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1700?message=Update%20Docs)]
</div>
**LLM Docstring**

Fetch the dipole-polarizability derivatives in Cartesian coordinates, computing them via `calculate_dipole_polarizability` (and caching the result, if the same evaluator is configured as the default) if not already available to the requested order, or re-expanding them through the normal modes if they were stored in a mode basis smaller than the full Cartesian space.
  - `order`: `int | None`
    > the highest derivative order needed; if `None`, whatever is available is returned
  - `evaluator`: `object | None`
    > an explicit polarizability evaluator to use instead of `self.polarizability_evaluator`
  - `include_constant_term`: `bool`
    > accepted for interface consistency with `get_cartesian_dipole_derivatives` but not used in this method's body
  - `:returns`: `list[np.ndarray]`
    > the polarizability derivative tensors


<a id="Psience.Molecools.Molecule.Molecule.get_hamiltonian" class="docs-object-method">&nbsp;</a> 
```python
get_hamiltonian(self, embedding=None, potential_derivatives=None, modes=None, dipole_derivatives=None, **etc): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1742)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1742?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a `MolecularHamiltonian` for this molecule, defaulting to its own embedding, potential derivatives, normal modes, and dipole derivatives wherever not explicitly overridden.
  - `embedding`: `MolecularEmbedding | None`
    > an alternate coordinate embedding to use
  - `potential_derivatives`: `list[np.ndarray] | None`
    > alternate potential-energy derivative tensors to use
  - `modes`: `object | None`
    > alternate normal modes to use
  - `dipole_derivatives`: `list[np.ndarray] | None`
    > alternate dipole derivative tensors to use
  - `etc`: `dict`
    > extra options forwarded to the `MolecularHamiltonian` constructor
  - `:returns`: `MolecularHamiltonian`
    > the constructed Hamiltonian


<a id="Psience.Molecools.Molecule.Molecule.hamiltonian" class="docs-object-method">&nbsp;</a> 
```python
@property
hamiltonian(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1781)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1781?message=Update%20Docs)]
</div>
**LLM Docstring**

Property getter/setter for the `MolecularHamiltonian`. The getter lazily builds one via `get_hamiltonian` the first time it's needed.
  - `e`: `MolecularHamiltonian`
    > (setter only) the new Hamiltonian
  - `:returns`: `MolecularHamiltonian`
    > (getter) the cached (or newly built) Hamiltonian


<a id="Psience.Molecools.Molecule.Molecule.potential_surface" class="docs-object-method">&nbsp;</a> 
```python
@property
potential_surface(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1810)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1810?message=Update%20Docs)]
</div>

  - `:returns`: `PotentialSurfaceManager`
    >


<a id="Psience.Molecools.Molecule.Molecule.potential_derivatives" class="docs-object-method">&nbsp;</a> 
```python
@property
potential_derivatives(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1835)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1835?message=Update%20Docs)]
</div>
**LLM Docstring**

Property getter/setter for the potential-energy derivative tensors. The getter fetches them from `self.potential_surface`, normalizing missing/placeholder entries to `0` and trimming any trailing zero-padded (unset) higher-order terms; the setter assigns to `self.potential_surface.derivatives`.
  - `derivs`: `list[np.ndarray]`
    > (setter only) the new potential derivative tensors
  - `:returns`: `list[np.ndarray] | None`
    > (getter) the potential derivative tensors with trailing unset orders trimmed, or `None` if unavailable


<a id="Psience.Molecools.Molecule.Molecule.get_cartesian_potential_derivatives" class="docs-object-method">&nbsp;</a> 
```python
get_cartesian_potential_derivatives(self, order=None, evaluator=None, use_cached=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1879)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1879?message=Update%20Docs)]
</div>
**LLM Docstring**

Fetch the potential-energy derivatives in Cartesian coordinates, computing them via `calculate_energy` (and caching the result on the molecule, if the same evaluator is configured as the default) if not already available to the requested order.
  - `order`: `int | None`
    > the highest derivative order needed; defaults to `2` if a fresh calculation is required
  - `evaluator`: `object | None`
    > an explicit energy evaluator to use instead of `self.energy_evaluator`
  - `use_cached`: `bool`
    > accepted for interface consistency but not used in this method's body
  - `:returns`: `list[np.ndarray] | None`
    > the potential derivative tensors, truncated to `order` if given, or `None` if unavailable


<a id="Psience.Molecools.Molecule.Molecule.get_internal_potential_derivatives" class="docs-object-method">&nbsp;</a> 
```python
get_internal_potential_derivatives(self, order=None, reembed=True, strip_embedding=True, zero_gradient=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1910)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1910?message=Update%20Docs)]
</div>
**LLM Docstring**

Fetch the potential-energy derivatives re-expressed in internal coordinates, by re-expanding the Cartesian potential derivatives through the Cartesians-by-internals Jacobian.
  - `order`: `int | None`
    > the highest derivative order needed
  - `reembed`: `bool`
    > whether to use the Eckart-reembedded Cartesians-by-internals Jacobian
  - `strip_embedding`: `bool`
    > whether to strip the fixed embedding coordinates from the Jacobian
  - `zero_gradient`: `bool`
    > whether to zero out the first-order (gradient) term before re-expanding
  - `:returns`: `list[np.ndarray] | None`
    > the internal-coordinate potential derivative tensors, or `None` if unavailable


<a id="Psience.Molecools.Molecule.Molecule.normal_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
normal_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1939)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1939?message=Update%20Docs)]
</div>

  - `:returns`: `NormalModesManager`
    >


<a id="Psience.Molecools.Molecule.Molecule.get_normal_modes" class="docs-object-method">&nbsp;</a> 
```python
get_normal_modes(self, masses=None, potential_derivatives=None, use_internals=None, project_transrot=True, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1964)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1964?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute this molecule's normal modes (via `NormalModes.from_molecule`), optionally using alternate masses/potential derivatives and controlling whether internal coordinates and translation/rotation projection are used.
  - `masses`: `np.ndarray | None`
    > masses to use instead of this molecule's own
  - `potential_derivatives`: `list[np.ndarray] | None`
    > potential derivatives to use instead of this molecule's own
  - `use_internals`: `bool | None`
    > whether to compute the modes in internal coordinates rather than Cartesians
  - `project_transrot`: `bool`
    > whether to project out translational/rotational degrees of freedom
  - `opts`: `dict`
    > extra options forwarded to `NormalModes.from_molecule`
  - `:returns`: `NormalModes`
    > the computed normal modes


<a id="Psience.Molecools.Molecule.Molecule.get_reaction_path_modes" class="docs-object-method">&nbsp;</a> 
```python
get_reaction_path_modes(self, masses=None, potential_derivatives=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L1995)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L1995?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute reaction-path-following normal modes for this molecule, via `ReactionPathModes.from_molecule`.
  - `masses`: `np.ndarray | None`
    > masses to use instead of this molecule's own
  - `potential_derivatives`: `list[np.ndarray] | None`
    > potential derivatives to use instead of this molecule's own
  - `opts`: `dict`
    > extra options forwarded to `ReactionPathModes.from_molecule`
  - `:returns`: `ReactionPathModes`
    > the computed reaction-path modes


<a id="Psience.Molecools.Molecule.Molecule.metadata" class="docs-object-method">&nbsp;</a> 
```python
@property
metadata(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2017)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2017?message=Update%20Docs)]
</div>
**LLM Docstring**

Property getter/setter for the molecule's metadata dict. The setter requires the new value to already be a `dict`.
  - `val`: `dict`
    > (setter only) the new metadata dict
  - `:returns`: `dict`
    > (getter) the metadata dict


<a id="Psience.Molecools.Molecule.Molecule.get_harmonic_spectrum" class="docs-object-method">&nbsp;</a> 
```python
get_harmonic_spectrum(self, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2051)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2051?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a harmonic IR spectrum for this molecule, via `HarmonicSpectrum.from_mol`.
  - `opts`: `dict`
    > extra options forwarded to `HarmonicSpectrum.from_mol`
  - `:returns`: `HarmonicSpectrum`
    > the constructed harmonic spectrum


<a id="Psience.Molecools.Molecule.Molecule.get_harmonic_raman_spectrum" class="docs-object-method">&nbsp;</a> 
```python
get_harmonic_raman_spectrum(self, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2065)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2065?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a harmonic Raman spectrum for this molecule, via `HarmonicSpectrum.raman_from_mol`.
  - `opts`: `dict`
    > extra options forwarded to `HarmonicSpectrum.raman_from_mol`
  - `:returns`: `HarmonicSpectrum`
    > the constructed harmonic Raman spectrum


<a id="Psience.Molecools.Molecule.Molecule.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2080)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2080?message=Update%20Docs)]
</div>
**LLM Docstring**

Debug string representation using the molecule's name (falling back to the source-file basename, then an SMILES string from any attached RDKit molecule, then the chemical formula) and its coordinate shape.
  - `:returns`: `str`
    > string of the form `ClassName('name', shape=coord_shape)`


<a id="Psience.Molecools.Molecule.Molecule.num_atoms" class="docs-object-method">&nbsp;</a> 
```python
@property
num_atoms(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2109)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2109?message=Update%20Docs)]
</div>
**LLM Docstring**

The number of atoms in the molecule.
  - `:returns`: `int`
    > the atom count


<a id="Psience.Molecools.Molecule.Molecule.atom_positions" class="docs-object-method">&nbsp;</a> 
```python
@property
atom_positions(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2120)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2120?message=Update%20Docs)]
</div>
**LLM Docstring**

A mapping from element symbol to the list of atom indices having that symbol.
  - `:returns`: `dict`
    > the symbol-to-indices mapping


<a id="Psience.Molecools.Molecule.Molecule.dummy_positions" class="docs-object-method">&nbsp;</a> 
```python
@property
dummy_positions(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2137)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2137?message=Update%20Docs)]
</div>
**LLM Docstring**

The indices of any dummy (`'X'`) atoms in the molecule.
  - `:returns`: `list[int]`
    > the dummy-atom indices, or an empty list if there are none


<a id="Psience.Molecools.Molecule.Molecule.atoms" class="docs-object-method">&nbsp;</a> 
```python
@property
atoms(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2149)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2149?message=Update%20Docs)]
</div>
**LLM Docstring**

The element symbols of the molecule's atoms.
  - `:returns`: `tuple[str]`
    > the atom symbols


<a id="Psience.Molecools.Molecule.Molecule.atomic_masses" class="docs-object-method">&nbsp;</a> 
```python
@property
atomic_masses(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2179)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2179?message=Update%20Docs)]
</div>
**LLM Docstring**

The atomic masses in atomic units (electron masses), via `_atomic_masses`.
  - `:returns`: `np.ndarray`
    > the masses in atomic units


<a id="Psience.Molecools.Molecule.Molecule.bonds" class="docs-object-method">&nbsp;</a> 
```python
@property
bonds(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2190)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2190?message=Update%20Docs)]
</div>
**LLM Docstring**

Property getter/setter for the molecule's bonds. The getter lazily guesses the bonds (via `get_guessed_bonds`) if none are set and `self.guess_bonds` is enabled.
  - `b`: `list[tuple] | None`
    > (setter only) the new bonds
  - `:returns`: `list[tuple] | None`
    > (getter) the bonds, or `None` if unset and bond-guessing is disabled


<a id="Psience.Molecools.Molecule.Molecule.break_bonds" class="docs-object-method">&nbsp;</a> 
```python
break_bonds(self, bonds, use_rdkit=False, **rdopts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2218)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2218?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a copy of this molecule with the specified bonds removed, either by delegating to the attached RDKit molecule's `break_bonds` or by filtering the bond list directly.
  - `bonds`: `Iterable[tuple]`
    > the bonds to remove, each an atom-index pair
  - `use_rdkit`: `bool`
    > whether to perform the break via the RDKit molecule instead of filtering `self.bonds` directly
  - `rdopts`: `dict`
    > extra options forwarded to the RDKit `break_bonds` call
  - `:returns`: `Molecule`
    > the new molecule with the given bonds removed


<a id="Psience.Molecools.Molecule.Molecule.formula" class="docs-object-method">&nbsp;</a> 
```python
@property
formula(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2241)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2241?message=Update%20Docs)]
</div>
**LLM Docstring**

The molecule's chemical formula, via `self.prop('chemical_formula')`.
  - `:returns`: `str`
    > the chemical formula


<a id="Psience.Molecools.Molecule.Molecule.multiconfig" class="docs-object-method">&nbsp;</a> 
```python
@property
multiconfig(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2252)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2252?message=Update%20Docs)]
</div>
**LLM Docstring**

Whether this molecule holds multiple geometry configurations at once, delegating to `self.coords.multiconfig`.
  - `:returns`: `bool`
    > whether multiple configurations are stored


<a id="Psience.Molecools.Molecule.Molecule.name" class="docs-object-method">&nbsp;</a> 
```python
@property
name(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2263)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2263?message=Update%20Docs)]
</div>
**LLM Docstring**

The molecule's name, falling back to `"Unnamed"` if none was set.
  - `:returns`: `str`
    > the molecule's name


<a id="Psience.Molecools.Molecule.Molecule.source_file" class="docs-object-method">&nbsp;</a> 
```python
@property
source_file(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2277)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2277?message=Update%20Docs)]
</div>
**LLM Docstring**

Setter for the source file: if given a plain path string, infers the source `mode` from its extension and wraps both into a `{'file':..., 'mode':...}` dict; an already-structured dict is stored as-is.
  - `src`: `str | dict | None`
    > the new source file, as a path string or a `{'file', 'mode'}` dict
  - `:returns`: `None`
    > None


<a id="Psience.Molecools.Molecule.Molecule.source_mode" class="docs-object-method">&nbsp;</a> 
```python
@property
source_mode(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2314)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2314?message=Update%20Docs)]
</div>
**LLM Docstring**

The inferred/stored source-file mode (e.g. `'fchk'`, `'log'`), if a source file is set.
  - `:returns`: `str | None`
    > the source mode, or `None` if no source file is set or no mode was recorded


<a id="Psience.Molecools.Molecule.Molecule.shape" class="docs-object-method">&nbsp;</a> 
```python
@property
shape(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2327)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2327?message=Update%20Docs)]
</div>
**LLM Docstring**

The shape of the molecule's Cartesian coordinates.
  - `:returns`: `tuple[int]`
    > the coordinate array's shape


<a id="Psience.Molecools.Molecule.Molecule.__len__" class="docs-object-method">&nbsp;</a> 
```python
__len__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2338)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2338?message=Update%20Docs)]
</div>
**LLM Docstring**

The number of geometry configurations held by this molecule: the leading dimension of its coordinates if `multiconfig`, otherwise `1`.
  - `:returns`: `int`
    > the number of configurations


<a id="Psience.Molecools.Molecule.Molecule.copy" class="docs-object-method">&nbsp;</a> 
```python
copy(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2352)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2352?message=Update%20Docs)]
</div>
**LLM Docstring**

Make a copy of this molecule by calling `modify()` with no overrides.
  - `:returns`: `Molecule`
    > the copied molecule


<a id="Psience.Molecools.Molecule.Molecule.take_submolecule" class="docs-object-method">&nbsp;</a> 
```python
take_submolecule(self, pos): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2378)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2378?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a new molecule containing only the atoms at the given positions, remapping bonds to the new (sub)indexing and dropping any bonds that reference atoms outside the subset.
  - `pos`: `Iterable[int]`
    > the atom indices to keep, in the desired order for the submolecule
  - `:returns`: `Molecule`
    > the constructed submolecule


<a id="Psience.Molecools.Molecule.Molecule.prop" class="docs-object-method">&nbsp;</a> 
```python
prop(self, name, *args, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2412)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2412?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute a named derived molecular property by dispatching to the corresponding function on `MolecularProperties`.
  - `name`: `str`
    > the property name, matching an attribute of `MolecularProperties`
  - `args`: `tuple`
    > positional arguments forwarded to the property function
  - `kwargs`: `dict`
    > keyword arguments forwarded to the property function
  - `:returns`: `object`
    > the computed property value


<a id="Psience.Molecools.Molecule.Molecule.get_guessed_bonds" class="docs-object-method">&nbsp;</a> 
```python
get_guessed_bonds(self, mode=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2439)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2439?message=Update%20Docs)]
</div>
**LLM Docstring**

Guess the bonding arrangement for this molecule, either via RDKit (from the Cartesian geometry) or via `MolecularProperties.guessed_bonds`, depending on `mode`.
  - `mode`: `str | None`
    > the bond-guessing strategy to use (`'rdkit'` or another mode understood by `MolecularProperties.guessed_bonds`); defaults to `self.bond_guessing_mode`
  - `opts`: `dict`
    > extra options forwarded to the underlying bond-guessing routine
  - `:returns`: `list[tuple]`
    > the guessed bonds


<a id="Psience.Molecools.Molecule.Molecule.edge_graph" class="docs-object-method">&nbsp;</a> 
```python
@property
edge_graph(self) -> 'EdgeGraph': 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2467)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2467?message=Update%20Docs)]
</div>
**LLM Docstring**

The (cached) `EdgeGraph` representation of the molecule's bonding structure, built lazily via `MolecularProperties.edge_graph`.
  - `:returns`: `EdgeGraph`
    > the edge graph


<a id="Psience.Molecools.Molecule.Molecule.find_path" class="docs-object-method">&nbsp;</a> 
```python
find_path(self, atom1, atom2): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2481)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2481?message=Update%20Docs)]
</div>
**LLM Docstring**

Find a path between two atoms through the bonding graph.
  - `atom1`: `int`
    > the starting atom index
  - `atom2`: `int`
    > the ending atom index
  - `:returns`: `list[int]`
    > the path between the two atoms


<a id="Psience.Molecools.Molecule.Molecule.find_substructure" class="docs-object-method">&nbsp;</a> 
```python
find_substructure(self, pattern): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2496)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2496?message=Update%20Docs)]
</div>
**LLM Docstring**

Find matches of a substructure pattern within the molecule, delegating to the attached RDKit molecule.
  - `pattern`: `str`
    > the substructure pattern to search for (e.g. a SMARTS string)
  - `:returns`: `object`
    > the matching substructures


<a id="Psience.Molecools.Molecule.Molecule.apply_smarts" class="docs-object-method">&nbsp;</a> 
```python
apply_smarts(self, pattern): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2509)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2509?message=Update%20Docs)]
</div>
**LLM Docstring**

Apply a SMARTS reaction/transformation pattern to the molecule via RDKit, returning each resulting product as a new `Molecule`.
  - `pattern`: `str`
    > the SMARTS pattern to apply
  - `:returns`: `list[Molecule]`
    > the resulting molecules


<a id="Psience.Molecools.Molecule.Molecule.neighborhood" class="docs-object-method">&nbsp;</a> 
```python
neighborhood(self, loc, size=1): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2526)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2526?message=Update%20Docs)]
</div>
**LLM Docstring**

Find the atoms within a given graph-distance of a location in the bonding graph.
  - `loc`: `int`
    > the atom index to center the neighborhood on
  - `size`: `int`
    > the neighborhood radius (in bond-graph steps)
  - `:returns`: `tuple[int]`
    > the neighboring atom indices


<a id="Psience.Molecools.Molecule.Molecule.remove_hydrogens" class="docs-object-method">&nbsp;</a> 
```python
remove_hydrogens(self, positions=None, max=None, *, hydrogen_types=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2541)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2541?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a copy of this molecule with hydrogen-type atoms (`H`/`D`/`T` by default) removed, either all of them or only those neighboring specified positions (optionally capped at `max` per position).
  - `positions`: `int | Iterable[int] | None`
    > atom(s) whose neighboring hydrogens should be removed; if `None`, every hydrogen-type atom in the molecule is removed
  - `max`: `int | None`
    > maximum number of hydrogens to remove per position in `positions`
  - `hydrogen_types`: `set[str] | None`
    > the set of element symbols treated as hydrogen isotopes; defaults to `{'H', 'D', 'T'}`
  - `:returns`: `Molecule`
    > the molecule with the selected hydrogens removed


<a id="Psience.Molecools.Molecule.Molecule.fragment_embedding" class="docs-object-method">&nbsp;</a> 
```python
fragment_embedding(self, fragment_indices, ref=None, return_axes=False, view_inds=(1, 2), use_moments=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2572)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2572?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute a local coordinate frame (origin, offset vector, and an up-vector or full axis set) anchored at a fragment of the molecule, used as the reference frame for attaching or orienting substituents; falls back to center-of-mass/principal-axis reference points (encoded as indices `-1`/`-2`/`-3`) when the fragment doesn't have enough atoms of its own to define a frame.
  - `fragment_indices`: `int | Iterable[int]`
    > the atom index (or indices) defining the fragment to embed
  - `ref`: `Iterable[int] | None`
    > reference atom(s) (outside the fragment) used to anchor the origin/frame; computed from the local neighborhood if not given
  - `return_axes`: `bool`
    > whether to return a full 3x3 axis frame instead of just an up-vector
  - `view_inds`: `tuple[int, int]`
    > which two fragment-atom positions define the "view" direction used to build the axis frame
  - `use_moments`: `bool`
    > whether to derive the up-vector from the fragment's moments of inertia rather than from its first three atom positions
  - `:returns`: `tuple`
    > `(origin, offset, up_or_axes)` -- the reference origin point, the offset from origin to the fragment's first atom, and either an up-vector or (if `return_axes`) a full axis frame


<a id="Psience.Molecools.Molecule.Molecule.attach_functional_group" class="docs-object-method">&nbsp;</a> 
```python
attach_functional_group(self, target_fragment, atoms, new_coords, bonds='recompute', ref=None, masses=None, distance='auto', angle=0, dihedral=0, embedding='auto', bond_order=None, use_absolue_posititions=False, group_site=None) -> "'typing.Self'": 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2650)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2650?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a copy of this molecule with a new group of atoms (`atoms`/`new_coords`) attached at `target_fragment`, positioning and orienting the new group using the fragment's local reference frame (bond distance/angle/dihedral, or an explicit embedding), and splicing the corresponding bonds into the result; supports designating a `group_site` atom within the new group as its attachment point, in which case the method recurses after re-deriving the embedding/bonds relative to that site.
  - `target_fragment`: `int | Iterable[int]`
    > the atom(s) of this molecule the new group attaches to/replaces
  - `atoms`: `Iterable[str]`
    > the element symbols of the atoms in the new group
  - `new_coords`: `np.ndarray`
    > the (local) coordinates of the new group's atoms
  - `bonds`: `str | list | None`
    > bonds within the new group; `'recompute'` to guess them fresh, `None` to reuse `self.bonds` remapped, or an explicit bond list
  - `ref`: `Iterable[int] | None`
    > reference atom(s) used to anchor the attachment frame; computed automatically if not given
  - `masses`: `np.ndarray | None`
    > masses for the new group's atoms; looked up from `atoms` if not given
  - `distance`: `str | float | None`
    > the bond distance to place the new group at; `'auto'` to look it up from `BondData`, or `None`/a number
  - `angle`: `float`
    > rotation angle (about the up-vector) to apply to the new group
  - `dihedral`: `float`
    > rotation angle (about the offset axis) to apply to the new group
  - `embedding`: `str | tuple | np.ndarray | None`
    > the reference orientation for the new group; `'auto'` to derive it from moments of inertia, or an explicit `(origin, axes)`/axes specification
  - `bond_order`: `float | None`
    > the bond order connecting the new group to the target fragment; defaults to `1` (or inferred when `group_site` is used)
  - `use_absolue_posititions`: `bool`
    > whether `new_coords` should be used as absolute coordinates rather than being repositioned relative to the fragment frame
  - `group_site`: `int | None`
    > index (within `atoms`/`new_coords`) of the atom that should serve as the attachment point; if given, the method recurses with the group re-anchored at this site and that atom excluded from the final group
  - `:returns`: `Molecule`
    > the molecule with the new group attached


<a id="Psience.Molecools.Molecule.Molecule.find_heavy_atom_backbone" class="docs-object-method">&nbsp;</a> 
```python
find_heavy_atom_backbone(self, root=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2830)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2830?message=Update%20Docs)]
</div>
**LLM Docstring**

Find the longest chain of atoms in the bonding graph (the heavy-atom backbone), via `edge_graph.find_longest_chain`.
  - `root`: `int | None`
    > an atom index to force as one end of the chain
  - `:returns`: `list[int]`
    > the backbone atom indices, in chain order


<a id="Psience.Molecools.Molecule.Molecule.find_backbone_segments" class="docs-object-method">&nbsp;</a> 
```python
find_backbone_segments(self, root=None, initial_backbone=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2843)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2843?message=Update%20Docs)]
</div>
**LLM Docstring**

Split the bonding graph into backbone-connected segments, via `edge_graph.segment_by_chains`.
  - `root`: `int | None`
    > an atom index to anchor the segmentation at
  - `initial_backbone`: `Iterable[int] | None`
    > an initial backbone chain to build the segmentation around
  - `:returns`: `list`
    > the resulting chain segments


<a id="Psience.Molecools.Molecule.Molecule.get_backbone_zmatrix" class="docs-object-method">&nbsp;</a> 
```python
get_backbone_zmatrix(self, root=None, segments=None, return_remainder=False, return_segments=False, required_coordinates=None, isolated_coordinates=None, root_coordinates=None, initial_backbone=None, validate=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2858)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2858?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a Z-matrix for a (typically single-fragment) molecule by first segmenting its bonding graph into backbone chains (via `find_backbone_segments`, validating there are no duplicate atoms across segments), constructing the base Z-matrix graph from the bonds and segments, filling in any bonds missing from the initial graph, and (if requested) enforcing required/isolated/root coordinate constraints.
  - `root`: `int | None`
    > an atom index to anchor the backbone segmentation at
  - `segments`: `list | None`
    > precomputed backbone segments to use instead of calling `find_backbone_segments`
  - `return_remainder`: `bool`
    > whether to also return the bonds that had to be added beyond the base backbone graph
  - `return_segments`: `bool`
    > whether to also return the backbone segments used
  - `required_coordinates`: `Iterable | None`
    > internal coordinates that must appear in the resulting Z-matrix
  - `isolated_coordinates`: `Iterable | None`
    > coordinates that must be built in isolation from others
  - `root_coordinates`: `Iterable | None`
    > coordinates to anchor at the root of the Z-matrix
  - `initial_backbone`: `Iterable[int] | None`
    > an initial backbone chain to seed the segmentation with
  - `validate`: `bool`
    > whether to validate the Z-matrix construction at each step (duplicate atoms, valid additions)
  - `:returns`: `np.ndarray | tuple`
    > the Z-matrix, or a tuple additionally including the segments and/or remainder bonds depending on the `return_*` flags


<a id="Psience.Molecools.Molecule.Molecule.get_canonical_zmatrix" class="docs-object-method">&nbsp;</a> 
```python
get_canonical_zmatrix(self, ordering=None, validate=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2941)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2941?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a canonical Z-matrix ordering for the molecule from a canonical fragmentation of the bonding graph.
  - `ordering`: `np.ndarray | None`
    > the atom ordering to use as the basis for canonicalization; defaults to the natural `0..N` ordering
  - `validate`: `bool`
    > whether to validate each Z-matrix addition
  - `:returns`: `np.ndarray`
    > the canonical Z-matrix


<a id="Psience.Molecools.Molecule.Molecule.get_bond_zmatrix" class="docs-object-method">&nbsp;</a> 
```python
get_bond_zmatrix(self, fragments=None, segments=None, root=None, required_coordinates=None, isolated_coordinates=None, root_coordinates=None, attachment_points=None, check_attachment_points=True, validate=True, for_fragment=None, fragment_ordering=None, connect_fragments=True, initial_backbone=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L2996)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L2996?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a full Z-matrix for the molecule from its bonding graph, handling the single-fragment case via `get_backbone_zmatrix` directly and the multi-fragment case by building a per-fragment Z-matrix for each fragment (optionally reordering/rooting/filtering required coordinates per fragment) and then splicing them together into one connected Z-matrix (via `coordops.complex_zmatrix`) unless `connect_fragments` is `False`; can also be restricted to build the Z-matrix for just one fragment (`for_fragment`), in which case it recurses on the corresponding submolecule and reindexes the result back to the full atom numbering.
  - `fragments`: `list[Iterable[int]] | None`
    > explicit fragment atom-index groups to use instead of `self.fragment_indices`
  - `segments`: `list | None`
    > precomputed backbone segments for the single-fragment case
  - `root`: `int | Iterable | None`
    > root atom(s) to anchor the Z-matrix construction at (per fragment, in the multi-fragment case)
  - `required_coordinates`: `Iterable | None`
    > internal coordinates that must appear in the resulting Z-matrix
  - `isolated_coordinates`: `Iterable | None`
    > coordinates that must be built in isolation from others
  - `root_coordinates`: `Iterable | None`
    > coordinates to anchor at the root of the Z-matrix
  - `attachment_points`: `dict | Iterable | None`
    > explicit inter-fragment attachment points to use when connecting fragments
  - `check_attachment_points`: `bool`
    > whether to validate the attachment points used when connecting fragments
  - `validate`: `bool`
    > whether to validate each Z-matrix addition
  - `for_fragment`: `int | Iterable[int] | None`
    > restrict the Z-matrix construction to just this fragment (an index into `self.fragment_indices`, or an explicit list of atom indices), returning the result reindexed to the full molecule
  - `fragment_ordering`: `Iterable[int] | None`
    > explicit ordering to apply to the fragments before connecting them
  - `connect_fragments`: `bool`
    > whether to splice the per-fragment Z-matrices into one connected Z-matrix, or return them as a list of separate Z-matrices
  - `initial_backbone`: `Iterable[int] | None`
    > an initial backbone chain to seed the segmentation of each fragment with
  - `:returns`: `np.ndarray | list[np.ndarray]`
    > the (connected) Z-matrix, or a list of per-fragment Z-matrices if `connect_fragments` is `False`


<a id="Psience.Molecools.Molecule.Molecule.fragment_indices" class="docs-object-method">&nbsp;</a> 
```python
@property
fragment_indices(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L3201)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L3201?message=Update%20Docs)]
</div>
**LLM Docstring**

The (cached) grouping of atom indices into connected molecular fragments, computed lazily via `MolecularProperties.fragment_indices`.
  - `:returns`: `list[np.ndarray]`
    > the list of per-fragment atom-index groups


<a id="Psience.Molecools.Molecule.Molecule.fragments" class="docs-object-method">&nbsp;</a> 
```python
@property
fragments(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L3215)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L3215?message=Update%20Docs)]
</div>
**LLM Docstring**

The molecule split into its connected fragments (as separate `Molecule` objects), via `MolecularProperties.fragments`.
  - `:returns`: `list[Molecule]`
    > the list of fragment molecules


<a id="Psience.Molecools.Molecule.Molecule.mass_weighted_coords" class="docs-object-method">&nbsp;</a> 
```python
@property
mass_weighted_coords(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L3228)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L3228?message=Update%20Docs)]
</div>

  - `:returns`: `CoordinateSet`
    >


<a id="Psience.Molecools.Molecule.Molecule.center_of_mass" class="docs-object-method">&nbsp;</a> 
```python
@property
center_of_mass(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L3236)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L3236?message=Update%20Docs)]
</div>

  - `:returns`: `CoordinateSet`
    >


<a id="Psience.Molecools.Molecule.Molecule.inertia_tensor" class="docs-object-method">&nbsp;</a> 
```python
@property
inertia_tensor(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L3243)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L3243?message=Update%20Docs)]
</div>

  - `:returns`: `(np.ndarray, np.ndarray)`
    >


<a id="Psience.Molecools.Molecule.Molecule.inertial_eigensystem" class="docs-object-method">&nbsp;</a> 
```python
@property
inertial_eigensystem(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L3250)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L3250?message=Update%20Docs)]
</div>

  - `:returns`: `(np.ndarray, np.ndarray)`
    >


<a id="Psience.Molecools.Molecule.Molecule.moments_of_inertia" class="docs-object-method">&nbsp;</a> 
```python
@property
moments_of_inertia(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L3257)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L3257?message=Update%20Docs)]
</div>

  - `:returns`: `np.ndarray`
    >


<a id="Psience.Molecools.Molecule.Molecule.inertial_axes" class="docs-object-method">&nbsp;</a> 
```python
@property
inertial_axes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L3264)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L3264?message=Update%20Docs)]
</div>

  - `:returns`: `np.ndarray`
    >


<a id="Psience.Molecools.Molecule.Molecule.translation_rotation_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
translation_rotation_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L3272)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L3272?message=Update%20Docs)]
</div>

  - `:returns`: `np.ndarray`
    >


<a id="Psience.Molecools.Molecule.Molecule.get_translation_rotation_projector" class="docs-object-method">&nbsp;</a> 
```python
get_translation_rotation_projector(self, mass_weighted=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L3280)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L3280?message=Update%20Docs)]
</div>
**LLM Docstring**

Build the projector matrix that removes translational and rotational displacement components from a Cartesian displacement.
  - `mass_weighted`: `bool`
    > whether the projector should act on mass-weighted coordinates
  - `:returns`: `np.ndarray`
    > the translation/rotation projector


<a id="Psience.Molecools.Molecule.Molecule.get_translation_rotation_invariant_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_translation_rotation_invariant_transformation(self, mass_weighted=False, strip_embedding=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L3304)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L3304?message=Update%20Docs)]
</div>
**LLM Docstring**

Build the transformation (and its inverse) that projects out the translational and rotational degrees of freedom from this molecule's Cartesian coordinates.
  - `mass_weighted`: `bool`
    > whether the transformation should act on mass-weighted coordinates
  - `strip_embedding`: `bool`
    > whether to strip the fixed embedding coordinates from the result
  - `:returns`: `tuple`
    > the translation/rotation-invariant transformation and its inverse


<a id="Psience.Molecools.Molecule.Molecule.energy" class="docs-object-method">&nbsp;</a> 
```python
@property
energy(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L3375)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L3375?message=Update%20Docs)]
</div>
**LLM Docstring**

The molecule's total energy: loaded from its source file if possible, else computed via `calculate_energy` if an energy evaluator is available, and cached thereafter.
  - `:returns`: `float | None`
    > the total energy, or `None` if neither a source-file energy nor an evaluator is available


<a id="Psience.Molecools.Molecule.Molecule.get_energy_evaluator" class="docs-object-method">&nbsp;</a> 
```python
get_energy_evaluator(self, evaluator=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L3392)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L3392?message=Update%20Docs)]
</div>
**LLM Docstring**

Resolve (and, if needed, instantiate from this molecule) an energy-evaluator object, defaulting to `self.energy_evaluator` (then `self.default_energy_evalutor`) if none is given explicitly.
  - `evaluator`: `object | None`
    > an explicit evaluator (or evaluator-type specification) to resolve
  - `opts`: `dict`
    > extra options forwarded to the evaluator's `from_mol` constructor, if applicable
  - `:returns`: `object`
    > the resolved energy-evaluator instance


<a id="Psience.Molecools.Molecule.Molecule.get_energy_function" class="docs-object-method">&nbsp;</a> 
```python
get_energy_function(self, evaluator=None, *, order=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L3418)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L3418?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a standalone callable that evaluates this molecule's energy (and derivatives) at arbitrary coordinates using the resolved energy evaluator.
  - `evaluator`: `object | None`
    > an explicit evaluator (or evaluator-type specification) to use
  - `order`: `int | None`
    > the derivative order the returned function should evaluate to by default
  - `opts`: `dict`
    > extra options forwarded to `get_energy_evaluator`
  - `:returns`: `callable`
    > a function `(coords, order=order) -> energy_or_expansion`


<a id="Psience.Molecools.Molecule.Molecule.calculate_energy" class="docs-object-method">&nbsp;</a> 
```python
calculate_energy(self, coords=None, *, evaluator=None, order=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L3459)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L3459?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the energy (and, optionally, its derivatives) at the given (or this molecule's own) coordinates using the resolved energy evaluator.
  - `coords`: `np.ndarray | None`
    > the coordinates to evaluate at, in Bohr; defaults to `self.coords`
  - `evaluator`: `object | None`
    > an explicit evaluator (or evaluator-type specification) to use
  - `order`: `int | None`
    > the highest derivative order to compute; if `None`, only the energy itself is returned
  - `opts`: `dict`
    > extra options forwarded to `get_energy_evaluator`
  - `:returns`: `float | list`
    > the energy, or the full energy/derivative expansion


<a id="Psience.Molecools.Molecule.Molecule.partial_force_field" class="docs-object-method">&nbsp;</a> 
```python
partial_force_field(self, coords=None, modes=None, *, evaluator=None, order=4, mesh_spacing=1, analytic_derivative_order=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L3487)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L3487?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute a partial (mode-selected) force-field expansion of the potential, evaluated in the given normal modes via the resolved energy evaluator's `partial_force_field` method.
  - `coords`: `np.ndarray | None`
    > the coordinates to evaluate at; defaults to `self.coords`
  - `modes`: `object | None`
    > the normal modes to expand in; computed via `get_normal_modes` if not given
  - `evaluator`: `object | None`
    > an explicit evaluator (or evaluator-type specification) to use
  - `order`: `int`
    > the highest derivative order to compute
  - `mesh_spacing`: `float`
    > finite-difference step size to use for the underlying force-field evaluation
  - `analytic_derivative_order`: `int | None`
    > order up to which derivatives should be computed analytically rather than numerically
  - `opts`: `dict`
    > extra options forwarded to `get_energy_evaluator`
  - `:returns`: `object`
    > the partial force-field expansion


<a id="Psience.Molecools.Molecule.Molecule.optimize" class="docs-object-method">&nbsp;</a> 
```python
optimize(self, evaluator=None, *, method=None, tol=None, max_iterations=None, logger=None, reembed=True, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L3530)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L3530?message=Update%20Docs)]
</div>
**LLM Docstring**

Geometry-optimize the molecule using the resolved energy evaluator, splitting `opts` into optimizer-specific and evaluator-specific options, re-embedding the optimized (and any trajectory) coordinates back onto the original geometry via an Eckart embedding, and returning a modified copy of the molecule at the optimized geometry.
  - `evaluator`: `object | None`
    > an explicit evaluator (or evaluator-type specification) to use
  - `method`: `str | None`
    > the optimization method/algorithm to use
  - `tol`: `float | None`
    > convergence tolerance for the optimizer
  - `max_iterations`: `int | None`
    > maximum number of optimizer iterations
  - `logger`: `Logger | str | None`
    > logger to report optimization progress/warnings to
  - `reembed`: `bool`
    > whether to Eckart-reembed the optimized geometry (and trajectory) onto the original one
  - `opts`: `dict`
    > extra options, split between optimizer options and evaluator-construction options
  - `:returns`: `Molecule | tuple[Molecule, list]`
    > the optimized molecule, or `(optimized_molecule, trajectory)` if the optimizer returned a trajectory


<a id="Psience.Molecools.Molecule.Molecule.relaxed_scan" class="docs-object-method">&nbsp;</a> 
```python
relaxed_scan(self, scan_values, scan_coordinates, evaluator=None, *, method=None, tol=None, max_iterations=None, logger=None, reembed=False, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L3603)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L3603?message=Update%20Docs)]
</div>
**LLM Docstring**

Perform a relaxed (constrained) potential-energy scan over the given coordinates/values using the resolved energy evaluator, optionally Eckart-reembedding the resulting trajectory onto the original geometry.
  - `scan_values`: `object`
    > the coordinate value(s) to scan over
  - `scan_coordinates`: `object`
    > the coordinate(s) to hold fixed/scan along
  - `evaluator`: `object | None`
    > an explicit evaluator (or evaluator-type specification) to use
  - `method`: `str | None`
    > the optimization method/algorithm to use for the constrained relaxations
  - `tol`: `float | None`
    > convergence tolerance for the optimizer
  - `max_iterations`: `int | None`
    > maximum number of optimizer iterations
  - `logger`: `Logger | str | None`
    > logger to report progress/warnings to
  - `reembed`: `bool`
    > whether to Eckart-reembed the resulting trajectory onto the original geometry
  - `opts`: `dict`
    > extra options, split between optimizer/scan options and evaluator-construction options
  - `:returns`: `tuple`
    > `(opts, traj, meta)` -- the resolved options, the resulting geometry trajectory, and scan metadata


<a id="Psience.Molecools.Molecule.Molecule.get_dipole_evaluator" class="docs-object-method">&nbsp;</a> 
```python
get_dipole_evaluator(self, evaluator=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L3676)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L3676?message=Update%20Docs)]
</div>
**LLM Docstring**

Resolve (and, if needed, instantiate from this molecule) a dipole-evaluator object, defaulting to `self.dipole_evaluator` if none is given explicitly.
  - `evaluator`: `object | None`
    > an explicit evaluator (or evaluator-type specification) to resolve
  - `opts`: `dict`
    > extra options forwarded to the evaluator's `from_mol` constructor, if applicable
  - `:returns`: `object`
    > the resolved dipole-evaluator instance


<a id="Psience.Molecools.Molecule.Molecule.get_dipole_function" class="docs-object-method">&nbsp;</a> 
```python
get_dipole_function(self, evaluator=None, *, order=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L3700)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L3700?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a standalone callable that evaluates this molecule's dipole (and derivatives) at arbitrary coordinates using the resolved dipole evaluator.
  - `evaluator`: `object | None`
    > an explicit evaluator (or evaluator-type specification) to use
  - `order`: `int | None`
    > the derivative order the returned function should evaluate to by default
  - `opts`: `dict`
    > extra options forwarded to `get_dipole_evaluator`
  - `:returns`: `callable`
    > a function `(coords, order=order) -> dipole_or_expansion`


<a id="Psience.Molecools.Molecule.Molecule.calculate_dipole" class="docs-object-method">&nbsp;</a> 
```python
calculate_dipole(self, evaluator=None, order=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L3742)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L3742?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the dipole moment (and, optionally, its derivatives) at this molecule's current coordinates using the resolved dipole evaluator.
  - `evaluator`: `object | None`
    > an explicit evaluator (or evaluator-type specification) to use
  - `order`: `int | None`
    > the highest derivative order to compute; if `None`, only the dipole itself is returned
  - `opts`: `dict`
    > extra options forwarded to `get_dipole_evaluator`
  - `:returns`: `np.ndarray | list`
    > the dipole, or the full dipole/derivative expansion


<a id="Psience.Molecools.Molecule.Molecule.get_polarizability_evaluator" class="docs-object-method">&nbsp;</a> 
```python
get_polarizability_evaluator(self, evaluator=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L3767)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L3767?message=Update%20Docs)]
</div>
**LLM Docstring**

Resolve (and, if needed, instantiate from this molecule) a dipole-polarizability-evaluator object, defaulting to `self.polarizability_evaluator` if none is given explicitly.
  - `evaluator`: `object | None`
    > an explicit evaluator (or evaluator-type specification) to resolve
  - `opts`: `dict`
    > extra options forwarded to the evaluator's `from_mol` constructor, if applicable
  - `:returns`: `object`
    > the resolved polarizability-evaluator instance


<a id="Psience.Molecools.Molecule.Molecule.get_dipole_polarizability_function" class="docs-object-method">&nbsp;</a> 
```python
get_dipole_polarizability_function(self, evaluator=None, *, order=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L3790)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L3790?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a standalone callable that evaluates this molecule's dipole polarizability (and derivatives) at arbitrary coordinates using the resolved polarizability evaluator.
  - `evaluator`: `object | None`
    > an explicit evaluator (or evaluator-type specification) to use
  - `order`: `int | None`
    > the derivative order the returned function should evaluate to by default
  - `opts`: `dict`
    > extra options forwarded to `get_polarizability_evaluator`
  - `:returns`: `callable`
    > a function `(coords, order=order) -> polarizability_or_expansion`


<a id="Psience.Molecools.Molecule.Molecule.calculate_dipole_polarizability" class="docs-object-method">&nbsp;</a> 
```python
calculate_dipole_polarizability(self, evaluator=None, order=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L3831)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L3831?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the dipole polarizability (and, optionally, its derivatives) at this molecule's current coordinates using the resolved polarizability evaluator.
  - `evaluator`: `object | None`
    > an explicit evaluator (or evaluator-type specification) to use
  - `order`: `int | None`
    > the highest derivative order to compute; if `None`, only the zeroth-order term of each returned quantity is kept
  - `opts`: `dict`
    > extra options forwarded to `get_polarizability_evaluator`
  - `:returns`: `list`
    > the polarizability expansion


<a id="Psience.Molecools.Molecule.Molecule.polarizability_derivatives" class="docs-object-method">&nbsp;</a> 
```python
@property
polarizability_derivatives(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L3855)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L3855?message=Update%20Docs)]
</div>
**LLM Docstring**

Property getter/setter for the dipole-polarizability derivative tensors, delegating to `self.dipole_surface.polarizability_derivatives`.
  - `derivs`: `list[np.ndarray]`
    > (setter only) the new polarizability derivative tensors
  - `:returns`: `list[np.ndarray]`
    > (getter) the polarizability derivative tensors


<a id="Psience.Molecools.Molecule.Molecule.get_reduced_potential_generator" class="docs-object-method">&nbsp;</a> 
```python
get_reduced_potential_generator(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L3882)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L3882?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a `ReducedDimensionalPotentialHandler` for generating reduced-dimensionality potential slices of this molecule.
  - `:returns`: `ReducedDimensionalPotentialHandler`
    > the constructed handler


<a id="Psience.Molecools.Molecule.Molecule.get_1d_potentials" class="docs-object-method">&nbsp;</a> 
```python
get_1d_potentials(self, spec, evaluator=None, energy_expansion=None, potential_params=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L3892)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L3892?message=Update%20Docs)]
</div>
**LLM Docstring**

Generate 1D potential-energy slices along the specified coordinate(s), via a `ReducedDimensionalPotentialHandler`.
  - `spec`: `object`
    > the coordinate specification(s) to slice along
  - `evaluator`: `object | None`
    > an explicit energy evaluator to use for the potential evaluations
  - `energy_expansion`: `list[np.ndarray] | None`
    > a precomputed energy expansion to use instead of evaluating fresh
  - `potential_params`: `dict | None`
    > extra parameters controlling the potential generation
  - `opts`: `dict`
    > extra options forwarded to `ReducedDimensionalPotentialHandler.get_1d_potentials`
  - `:returns`: `object`
    > the generated 1D potentials


<a id="Psience.Molecools.Molecule.Molecule.evaluate" class="docs-object-method">&nbsp;</a> 
```python
evaluate(self, func, use_internals=None, order=None, strip_embedding=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L3925)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L3925?message=Update%20Docs)]
</div>
**LLM Docstring**

Evaluate an arbitrary function of the molecule's coordinates (and its derivatives) via the molecule's `MolecularEvaluator`.
  - `func`: `callable`
    > the function to evaluate
  - `use_internals`: `bool | None`
    > whether to evaluate/differentiate in internal coordinates rather than Cartesians
  - `order`: `int | None`
    > the highest derivative order to compute
  - `strip_embedding`: `bool`
    > whether to strip the fixed embedding coordinates from the result
  - `:returns`: `object`
    > the evaluated function value(s)/derivatives


<a id="Psience.Molecools.Molecule.Molecule.evaluate_at" class="docs-object-method">&nbsp;</a> 
```python
evaluate_at(self, func, coords, use_internals=None, order=None, strip_embedding=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L3953)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L3953?message=Update%20Docs)]
</div>
**LLM Docstring**

Evaluate an arbitrary function of the molecule's coordinates (and its derivatives) at a specific set of coordinates, via the molecule's `MolecularEvaluator`.
  - `func`: `callable`
    > the function to evaluate
  - `coords`: `np.ndarray`
    > the coordinates to evaluate at
  - `use_internals`: `bool | None`
    > whether to evaluate/differentiate in internal coordinates rather than Cartesians
  - `order`: `int | None`
    > the highest derivative order to compute
  - `strip_embedding`: `bool`
    > whether to strip the fixed embedding coordinates from the result
  - `:returns`: `object`
    > the evaluated function value(s)/derivatives


<a id="Psience.Molecools.Molecule.Molecule.get_displaced_coordinates" class="docs-object-method">&nbsp;</a> 
```python
get_displaced_coordinates(self, displacements, which=None, sel=None, axes=None, use_internals=False, coordinate_expansion=None, strip_embedding=False, shift=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L3986)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L3986?message=Update%20Docs)]
</div>
**LLM Docstring**

Build displaced copies of the molecule's coordinates along specified atoms/axes, via the molecule's `MolecularEvaluator`.
  - `displacements`: `np.ndarray`
    > the displacement values to apply
  - `which`: `object | None`
    > which atoms/coordinates to displace
  - `sel`: `Iterable[int] | None`
    > a selection of atoms to restrict the displacement to
  - `axes`: `object | None`
    > which axes to displace along
  - `use_internals`: `bool`
    > whether the displacements are given in internal coordinates rather than Cartesians
  - `coordinate_expansion`: `list[np.ndarray] | None`
    > a coordinate-transformation expansion to apply to the displacements before applying them
  - `strip_embedding`: `bool`
    > whether the coordinate expansion has had its embedding coordinates stripped
  - `shift`: `bool`
    > whether the displacements are relative shifts (`True`) or absolute target values
  - `:returns`: `np.ndarray`
    > the displaced coordinates


<a id="Psience.Molecools.Molecule.Molecule.get_scan_coordinates" class="docs-object-method">&nbsp;</a> 
```python
get_scan_coordinates(self, domains, internals=False, modes=None, order=None, which=None, sel=None, axes=None, shift=True, coordinate_expansion=None, strip_embedding=False, return_displacements=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L4026)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L4026?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a grid of displaced coordinates spanning the given coordinate domains, optionally re-expressed through a normal-mode (or other) coordinate expansion, via the molecule's `MolecularEvaluator`.
  - `domains`: `object`
    > the coordinate ranges to scan over
  - `internals`: `bool`
    > whether the scan coordinates are internal coordinates rather than Cartesians
  - `modes`: `object | bool | None`
    > normal modes to scan along instead of raw coordinates; `True` to use this molecule's own modes (with translation/rotation projection disabled)
  - `order`: `int | None`
    > the Jacobian order to use when building the mode-based coordinate expansion
  - `which`: `object | None`
    > which atoms/coordinates to scan
  - `sel`: `Iterable[int] | None`
    > a selection of atoms to restrict the scan to
  - `axes`: `object | None`
    > which axes to scan along
  - `shift`: `bool`
    > whether the scan values are relative shifts or absolute target values
  - `coordinate_expansion`: `list[np.ndarray] | None`
    > an explicit coordinate-transformation expansion to combine with the mode-based one (if `modes` is given) or use directly
  - `strip_embedding`: `bool`
    > whether the coordinate expansion has had its embedding coordinates stripped
  - `return_displacements`: `bool`
    > whether to also return the raw displacement values used
  - `:returns`: `np.ndarray | tuple`
    > the scan coordinates, or additionally the displacements if `return_displacements` is set


<a id="Psience.Molecools.Molecule.Molecule.get_nearest_displacement_atoms" class="docs-object-method">&nbsp;</a> 
```python
get_nearest_displacement_atoms(self, points, sel=None, axes=None, weighting_function=None, return_distances=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L4091)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L4091?message=Update%20Docs)]
</div>
**LLM Docstring**

Find the atoms nearest to a set of query points (optionally restricted to a selection/axes, with a custom distance weighting), via the molecule's `MolecularEvaluator`.
  - `points`: `np.ndarray`
    > the query points to find nearest atoms for
  - `sel`: `Iterable[int] | None`
    > a selection of atoms to restrict the search to
  - `axes`: `object | None`
    > which coordinate axes to compute distances over
  - `weighting_function`: `callable | None`
    > a custom function to weight distances by
  - `return_distances`: `bool`
    > whether to also return the computed distances
  - `:returns`: `np.ndarray | tuple`
    > the nearest atom indices, or additionally the distances if `return_distances` is set


<a id="Psience.Molecools.Molecule.Molecule.get_nearest_displacement_coordinates" class="docs-object-method">&nbsp;</a> 
```python
get_nearest_displacement_coordinates(self, points, sel=None, axes=None, weighting_function=None, modes_nearest=False, return_distances=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L4119)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L4119?message=Update%20Docs)]
</div>
**LLM Docstring**

Find the displaced coordinates nearest to a set of query points, optionally in normal-mode space, via the molecule's `MolecularEvaluator`.
  - `points`: `np.ndarray`
    > the query points to find nearest coordinates for
  - `sel`: `Iterable[int] | None`
    > a selection of atoms to restrict the search to
  - `axes`: `object | None`
    > which coordinate axes to compute distances over
  - `weighting_function`: `callable | None`
    > a custom function to weight distances by
  - `modes_nearest`: `bool`
    > whether to find the nearest point in normal-mode space rather than Cartesian/internal space
  - `return_distances`: `bool`
    > whether to also return the computed distances
  - `:returns`: `np.ndarray | tuple`
    > the nearest coordinates, or additionally the distances if `return_distances` is set


<a id="Psience.Molecools.Molecule.Molecule.get_nearest_scan_coordinates" class="docs-object-method">&nbsp;</a> 
```python
get_nearest_scan_coordinates(self, domains, sel=None, axes=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L4151)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L4151?message=Update%20Docs)]
</div>
**LLM Docstring**

Find the scan-grid coordinates nearest to a set of query domains, via the molecule's `MolecularEvaluator`.
  - `domains`: `object`
    > the coordinate domains to build the nearest scan grid for
  - `sel`: `Iterable[int] | None`
    > a selection of atoms to restrict the search to
  - `axes`: `object | None`
    > which coordinate axes to compute distances over
  - `:returns`: `np.ndarray`
    > the nearest scan coordinates


<a id="Psience.Molecools.Molecule.Molecule.plot_molecule_function" class="docs-object-method">&nbsp;</a> 
```python
plot_molecule_function(self, function, *, axes, sel=None, embed=False, modes_nearest=False, domain=None, domain_padding=1, plot_points=500, weighting_function=None, mask_function=None, mask_value=0, plot_atoms=False, atom_colors=None, atom_radii=None, plotter=None, epilog=None, **plot_options): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L4189)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L4189?message=Update%20Docs)]
</div>
**LLM Docstring**

Plot an arbitrary scalar function of the molecule's geometry over a 1D or 2D grid of displaced coordinates, evaluating the function at the grid points nearest to each displacement, optionally masking out-of-range values and overlaying the atoms' positions.
  - `function`: `callable`
    > the function to evaluate, taking a batch of coordinates and returning scalar values
  - `axes`: `int | Iterable[int]`
    > which atom/coordinate axis (or axes, up to 2) to build the plot grid over
  - `sel`: `Iterable[int] | None`
    > a selection of atoms to restrict the displacement search to
  - `embed`: `bool`
    > whether to Eckart-embed the evaluation points before calling `function`
  - `modes_nearest`: `bool`
    > whether to find the nearest displacement in normal-mode space
  - `domain`: `np.ndarray | None`
    > explicit `(min, max)` bounds per axis; computed from the geometry's bounding box (plus `domain_padding`) if not given
  - `domain_padding`: `float | Iterable[float] | None`
    > padding to add around the automatically computed domain
  - `plot_points`: `int | Iterable[int]`
    > number of grid points per axis
  - `weighting_function`: `callable | None`
    > custom distance-weighting function used when finding nearest displacements
  - `mask_function`: `callable | None`
    > a function `(values, eval_points, dists) -> mask` used to blank out certain grid values
  - `mask_value`: `float`
    > the value to substitute where `mask_function` returns `True`
  - `plot_atoms`: `bool`
    > whether to overlay disks at the atoms' projected positions
  - `atom_colors`: `list | None`
    > per-atom overlay colors; defaults to each atom's icon color
  - `atom_radii`: `list | None`
    > per-atom overlay radii; defaults to each atom's resolved radius
  - `plotter`: `type | callable | None`
    > the plotting class/function to use; defaults to `plt.Plot` (1D) or `plt.TriContourPlot` (2D)
  - `epilog`: `list | None`
    > extra plot elements to draw on top of the function plot
  - `plot_options`: `dict`
    > extra options forwarded to the plotting call
  - `:returns`: `object`
    > the resulting plot


<a id="Psience.Molecools.Molecule.Molecule.get_model" class="docs-object-method">&nbsp;</a> 
```python
get_model(self, potential_specs, dipole=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L4367)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L4367?message=Update%20Docs)]
</div>
**LLM Docstring**

Build an analytic `MolecularModel` (a symbolic potential/dipole surface expressed in terms of this molecule's Z-matrix internal coordinates) from a specification of per-coordinate potential (and optional dipole) function contributions.
  - `potential_specs`: `dict | object`
    > either a raw potential expression, or a dict mapping internal-coordinate index (or index tuples, for coupled terms) to a function-type specification (a dict with a `'function_type'` key or a single `{function_type: params}` entry) describing that coordinate's contribution
  - `dipole`: `Iterable | None`
    > an optional 3-element (x/y/z) sequence of dipole-component specifications, each following the same per-coordinate specification format as `potential_specs`
  - `:returns`: `MolecularModel`
    > the constructed analytic molecular model


<a id="Psience.Molecools.Molecule.Molecule.get_point_group" class="docs-object-method">&nbsp;</a> 
```python
get_point_group(self, coords=None, masses=None, *, sel=None, verbose=False, return_identifier=False, **tols): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L4527)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L4527?message=Update%20Docs)]
</div>
**LLM Docstring**

Identify the molecular point group at the given (or this molecule's own) geometry, via `McUtils.Symmetry.PointGroupIdentifier`.
  - `coords`: `np.ndarray | None`
    > the coordinates to identify the point group for; defaults to `self.coords`
  - `masses`: `np.ndarray | None`
    > the masses to use; defaults to `self.atomic_masses`
  - `sel`: `Iterable[int] | None`
    > a selection of atoms to restrict the identification to
  - `verbose`: `bool`
    > whether the identifier should print diagnostic output
  - `return_identifier`: `bool`
    > whether to also return the underlying `PointGroupIdentifier` object
  - `tols`: `dict`
    > extra tolerance options forwarded to `PointGroupIdentifier`
  - `:returns`: `object | tuple`
    > the identified point group, or `(identifier, point_group)` if `return_identifier` is set


<a id="Psience.Molecools.Molecule.Molecule.point_group" class="docs-object-method">&nbsp;</a> 
```python
@property
point_group(self) -> 'symm.PointGroup': 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L4565)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L4565?message=Update%20Docs)]
</div>
**LLM Docstring**

Property getter/setter for the molecule's (cached) point group. The getter computes it lazily via `get_point_group` the first time it's needed.
  - `point_group`: `symm.PointGroup`
    > (setter only) the new point group to cache
  - `:returns`: `symm.PointGroup`
    > (getter) the cached (or newly computed) point group


<a id="Psience.Molecools.Molecule.Molecule.get_point_group_embedded_coordinates" class="docs-object-method">&nbsp;</a> 
```python
get_point_group_embedded_coordinates(self, pg=None, sel=None, return_point_group=False, return_identifier=False, **tols): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L4684)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L4684?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute this molecule's coordinates re-embedded into the standard reference frame of its molecular point group (translating to the center of mass, aligning to the principal axes, then aligning to the point group's own axis convention).
  - `pg`: `symm.PointGroup | None`
    > an explicit point group to embed into; identified automatically via `get_point_group` if not given
  - `sel`: `Iterable[int] | None`
    > a selection of atoms to restrict the point-group identification to
  - `return_point_group`: `bool`
    > whether to also return the (embedded) point group
  - `return_identifier`: `bool`
    > whether to also return the underlying `PointGroupIdentifier`
  - `tols`: `dict`
    > extra tolerance options forwarded to `PointGroupIdentifier`
  - `:returns`: `np.ndarray | tuple`
    > the point-group-embedded coordinates, or additionally the point group and/or identifier depending on the `return_*` flags


<a id="Psience.Molecools.Molecule.Molecule.symmetrize" class="docs-object-method">&nbsp;</a> 
```python
symmetrize(self, pg=None, return_identifier=False, tol=0.1, sel=None, return_coordinates=None, return_point_group=False, **tols): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L4737)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L4737?message=Update%20Docs)]
</div>
**LLM Docstring**

Symmetrize the molecule's geometry to conform exactly to a (identified or given) point group, adding/relabeling atoms as needed if the symmetrization changes the atom count, and matching the symmetrized geometry back onto the original atom ordering when the count is unchanged.
  - `pg`: `symm.PointGroup | None`
    > an explicit point group to symmetrize to; identified automatically via `get_point_group` if not given
  - `return_identifier`: `bool`
    > whether to also return the underlying `PointGroupIdentifier`
  - `tol`: `float`
    > tolerance used for point-group identification and symmetrization
  - `sel`: `Iterable[int] | None`
    > a selection of atoms to restrict the symmetrization to
  - `return_coordinates`: `bool | None`
    > whether to return raw coordinates instead of a modified `Molecule` (when the atom count is unchanged)
  - `return_point_group`: `bool`
    > whether to also return the point group used
  - `tols`: `dict`
    > extra tolerance options forwarded to `PointGroupIdentifier`
  - `:returns`: `tuple`
    > a tuple starting with the symmetrized molecule/coordinates (or `None` if the atom count changed and no direct match was possible), followed by `(new_atoms, new_coords)`, and optionally the identifier and/or point group


<a id="Psience.Molecools.Molecule.Molecule.get_symmetrized_internals" class="docs-object-method">&nbsp;</a> 
```python
get_symmetrized_internals(self, point_group=None, *, internals=None, extra_internals=None, masses=None, return_expansions=False, atom_selection=None, as_characters=True, normalize=None, drop_empty_modes=None, perms=None, return_base_expansion=False, return_point_group=False, reduce_redundant_coordinates=None, ops=None, permutation_tol=0.01, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L4816)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L4816?message=Update%20Docs)]
</div>
**LLM Docstring**

Build symmetry-adapted linear combinations of internal coordinates for the molecule under a (identified or given) point group, via `McUtils.Symmetry.symmetrize_internals`.
  - `point_group`: `symm.PointGroup | None`
    > an explicit point group to symmetrize under; identified automatically via `get_point_group` (or on a submolecule, if `atom_selection` is given) if not given
  - `internals`: `Iterable | None`
    > the internal coordinates to symmetrize; extracted from `self.internals` (Z-matrix or generic spec) if not given
  - `extra_internals`: `Iterable | None`
    > extra internal coordinates to prepend to the extracted set (deduplicated)
  - `masses`: `np.ndarray | None`
    > masses to use; defaults to `self.atomic_masses`
  - `return_expansions`: `bool`
    > whether to also return the coordinate-derivative expansions used
  - `atom_selection`: `Iterable[int] | None`
    > restrict the point-group identification to a submolecule
  - `as_characters`: `bool`
    > whether to express the symmetrized coordinates in terms of irrep characters
  - `normalize`: `bool | None`
    > whether to normalize the symmetrized coordinate combinations
  - `drop_empty_modes`: `bool | None`
    > whether to drop symmetry combinations that come out identically zero
  - `perms`: `Iterable | None`
    > explicit atom permutations to use for the symmetry operations
  - `return_base_expansion`: `bool`
    > whether to also return the base (pre-symmetrization) coordinate expansion
  - `return_point_group`: `bool`
    > whether to also return the point group used
  - `reduce_redundant_coordinates`: `bool | None`
    > whether to reduce a redundant coordinate set down to a non-redundant symmetrized basis
  - `ops`: `Iterable | None`
    > explicit symmetry operations to use instead of the full point-group operation set
  - `permutation_tol`: `float`
    > tolerance used when matching atom permutations to symmetry operations
  - `opts`: `dict`
    > extra options forwarded to point-group identification/`symmetrize_internals`
  - `:returns`: `object | tuple`
    > the symmetrized internal coordinates, plus any additionally requested return values, prefixed with the point group if `return_point_group` is set


<a id="Psience.Molecools.Molecule.Molecule.get_surface" class="docs-object-method">&nbsp;</a> 
```python
get_surface(self, radius_type='VanDerWaalsRadius', *, surface_type=None, radius_units='Angstroms', samples=100, radius_scaling=1, **etc): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L4923)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L4923?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a surface representation of the molecule (by default, a union of per-atom spheres) using per-atom radii of the requested type.
  - `radius_type`: `str`
    > the `AtomData` radius field to use for each atom
  - `surface_type`: `type | None`
    > the surface class to construct; defaults to `zach.SphereUnionSurface`
  - `radius_units`: `str`
    > the units the radii are given in before converting to Bohr
  - `samples`: `int`
    > number of samples used to represent the surface
  - `radius_scaling`: `float`
    > uniform scale factor applied to all radii
  - `etc`: `dict`
    > extra options forwarded to the surface constructor
  - `:returns`: `object`
    > the constructed surface


<a id="Psience.Molecools.Molecule.Molecule.get_surface_mesh" class="docs-object-method">&nbsp;</a> 
```python
get_surface_mesh(self, radius_type='VanDerWaalsRadius', *, surface_type=None, radius_units='Angstroms', samples=50, expansion=0.01, mesh_options=None, **etc): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L4966)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L4966?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a triangulated mesh of the molecule's surface, via `get_surface` followed by `generate_mesh`.
  - `radius_type`: `str`
    > the `AtomData` radius field to use for each atom
  - `surface_type`: `type | None`
    > the surface class to construct; defaults to `zach.SphereUnionSurface`
  - `radius_units`: `str`
    > the units the radii are given in before converting to Bohr
  - `samples`: `int`
    > number of samples used to represent the surface
  - `expansion`: `float`
    > expansion factor applied to the surface before meshing
  - `mesh_options`: `dict | None`
    > extra options forwarded to `generate_mesh`
  - `etc`: `dict`
    > extra options forwarded to `get_surface`
  - `:returns`: `object`
    > the generated surface mesh


<a id="Psience.Molecools.Molecule.Molecule.setup_AIMD" class="docs-object-method">&nbsp;</a> 
```python
setup_AIMD(self, potential_function=None, timestep=0.5, seed=None, total_energy=None, total_energy_scaling=None, trajectories=1, sampled_modes=None, initial_energies=None, initial_displacements=None, initial_mode_directions=None, displaced_coords=None, track_kinetic_energy=False, track_velocities=False, modes=None, **etc): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L5013)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L5013?message=Update%20Docs)]
</div>
**LLM Docstring**

Set up an `AIMDSimulator` (ab initio / classical molecular dynamics) for this molecule, either propagating from an explicit set of displaced initial positions, or sampling initial normal-mode velocities/energies (from explicit mode directions, explicit per-trajectory energies, or randomly-sampled directions distributing a target total energy) when starting from the equilibrium geometry.
  - `potential_function`: `callable | None`
    > the energy/gradient function to use; built via `get_energy_function` if not given
  - `timestep`: `float`
    > the simulation timestep
  - `seed`: `int | None`
    > random seed for sampling initial mode directions
  - `total_energy`: `float | None`
    > the total vibrational energy to distribute across sampled trajectories/modes
  - `total_energy_scaling`: `float | None`
    > scale factor applied to the default total energy (sum of mode frequencies) if `total_energy` isn't given
  - `trajectories`: `int`
    > number of trajectories to sample when randomly generating mode directions
  - `sampled_modes`: `Iterable[int] | None`
    > which normal modes to sample energy into; defaults to all of them
  - `initial_energies`: `np.ndarray | None`
    > explicit per-trajectory, per-mode energies to use instead of sampling them
  - `initial_displacements`: `np.ndarray | None`
    > explicit initial Cartesian displacements to start the trajectories from, bypassing the normal-mode energy-sampling path entirely
  - `initial_mode_directions`: `np.ndarray | None`
    > explicit per-trajectory mode-direction vectors to convert into initial energies
  - `displaced_coords`: `object | None`
    > which coordinates `initial_displacements` applies to
  - `track_kinetic_energy`: `bool`
    > whether the simulator should track kinetic energy over the trajectory
  - `track_velocities`: `bool`
    > whether the simulator should track velocities over the trajectory
  - `modes`: `object | None`
    > normal modes to use for the energy-to-velocity conversion; computed via `get_normal_modes` if not given
  - `etc`: `dict`
    > extra options forwarded to the `AIMDSimulator` constructor
  - `:returns`: `AIMDSimulator`
    > the constructed simulator


<a id="Psience.Molecools.Molecule.Molecule.setup_VPT" class="docs-object-method">&nbsp;</a> 
```python
setup_VPT(self, *, states=2, order=2, use_internals=None, potential_derivatives=None, energy_evaluator=None, dipole_derivatives=None, dipole_evaluator=None, runner='matrix', use_reaction_path=False, modes=None, projected_modes=None, mode_transformation=None, potential_terms=None, dipole_terms=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L5152)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L5152?message=Update%20Docs)]
</div>
**LLM Docstring**

Set up a vibrational perturbation theory (VPT) calculation for this molecule, resolving the potential/dipole derivative expansions (computing them via the configured evaluators if not already available), the normal modes (optionally localized against reaction-path-projected modes), and dispatching to the appropriate VPT runner (`VPTRunner` for the default matrix-based approach, or `AnalyticVPTRunner`).
  - `states`: `int | object`
    > the target vibrational states specification (e.g. max quanta) to compute
  - `order`: `int`
    > the perturbation-theory order to compute the potential/dipole expansions to
  - `use_internals`: `bool | None`
    > whether to run VPT in internal coordinates (constructing from `self.modify()`) rather than bare atoms/coordinates
  - `potential_derivatives`: `list[np.ndarray] | None`
    > explicit potential-energy derivative tensors to use instead of computing them
  - `energy_evaluator`: `object | None`
    > an explicit energy evaluator to use when computing the potential derivatives
  - `dipole_derivatives`: `list[np.ndarray] | None`
    > explicit dipole derivative tensors to use instead of computing them
  - `dipole_evaluator`: `object | None`
    > an explicit dipole evaluator to use when computing the dipole derivatives
  - `runner`: `str | object`
    > which VPT runner to use: `'matrix'` for `VPTRunner`, or an explicit runner class/object
  - `use_reaction_path`: `bool`
    > whether to project the normal modes against reaction-path modes
  - `modes`: `object | None`
    > explicit normal modes to use instead of computing them
  - `projected_modes`: `object | None`
    > explicit modes to localize the normal modes against
  - `mode_transformation`: `object | None`
    > an explicit mode-coordinate transformation to use
  - `potential_terms`: `object | None`
    > precomputed potential expansion terms, bypassing derivative computation entirely
  - `dipole_terms`: `object | None`
    > precomputed dipole expansion terms, bypassing derivative computation entirely
  - `opts`: `dict`
    > extra options forwarded to the VPT runner's `construct` method
  - `:returns`: `object`
    > the constructed VPT runner/calculation


<a id="Psience.Molecools.Molecule.Molecule.setup_job" class="docs-object-method">&nbsp;</a> 
```python
setup_job(self, job_type, *args, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L5285)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L5285?message=Update%20Docs)]
</div>
**LLM Docstring**

Set up an external-program computational job (e.g. a Gaussian/ORCA job) for this molecule, via `ExternalProgramJob.resolve` and `from_mol`.
  - `job_type`: `str | object`
    > the job type to resolve (a name or job-type object)
  - `args`: `tuple`
    > positional arguments forwarded to the resolved job type's `from_mol`
  - `kwargs`: `dict`
    > keyword arguments forwarded to the resolved job type's `from_mol`, merged over its default options
  - `:returns`: `object`
    > the constructed job


<a id="Psience.Molecools.Molecule.Molecule.get_gmatrix" class="docs-object-method">&nbsp;</a> 
```python
get_gmatrix(self, masses=None, coords=None, use_internals=None, power=None, **internals_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L5307)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L5307?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the G-matrix (inverse effective-mass metric tensor) for this molecule, either the trivial diagonal inverse-mass form (in plain Cartesians) or the projected `B G0 B^T`-style form built from the internals-by-Cartesians Jacobian (when using internal coordinates), optionally raised to a fractional power.
  - `masses`: `np.ndarray | None`
    > masses to use instead of this molecule's own
  - `coords`: `np.ndarray | None`
    > alternate coordinates to compute the G-matrix at
  - `use_internals`: `bool | None`
    > whether to compute the G-matrix in internal coordinates; defaults to `True` if the molecule has internal coordinates defined
  - `power`: `float | None`
    > an optional power to raise the resulting G-matrix to
  - `internals_opts`: `dict`
    > extra options forwarded to `get_internals_by_cartesians` when using internal coordinates
  - `:returns`: `np.ndarray`
    > the G-matrix


<a id="Psience.Molecools.Molecule.Molecule.g_matrix" class="docs-object-method">&nbsp;</a> 
```python
@property
g_matrix(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L5363)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L5363?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L5372)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L5372?message=Update%20Docs)]
</div>
Returns the molecular g-matrix for the system
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Molecule.Molecule.principle_axis_frame" class="docs-object-method">&nbsp;</a> 
```python
principle_axis_frame(self, sel=None, inverse=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L5381)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L5381?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L5395)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L5395?message=Update%20Docs)]
</div>
Gets the principle axis embedded coords and embedding parameters for the molecule
  - `:returns`: `MolecularTransformation | List[MolecularTransformation]`
    >


<a id="Psience.Molecools.Molecule.Molecule.permute_atoms" class="docs-object-method">&nbsp;</a> 
```python
permute_atoms(self, perm): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L5405)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L5405?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a copy of this molecule with its atoms reordered according to a permutation, remapping bonds to the new indexing accordingly.
  - `perm`: `Iterable[int]`
    > the new atom ordering, given as the sequence of old indices in their new order
  - `:returns`: `Molecule`
    > the permuted molecule


<a id="Psience.Molecools.Molecule.Molecule.apply_affine_transformation" class="docs-object-method">&nbsp;</a> 
```python
apply_affine_transformation(self, transformation, load_properties=False, embed_properties=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L5428)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L5428?message=Update%20Docs)]
</div>
**LLM Docstring**

Apply an affine (rotation/translation) transformation to the molecule's coordinates, optionally propagating the transformation onto its normal modes, potential surface, and dipole surface as well.
  - `transformation`: `np.ndarray | object`
    > the affine transformation to apply, either a raw transformation matrix or an object exposing an `apply` method (e.g. a `MolecularTransformation`)
  - `load_properties`: `bool | None`
    > whether to force-load the normal modes/potential/dipole derivatives before transforming (so they get carried over even if not already computed); `None` to load only if cheaply available
  - `embed_properties`: `bool`
    > whether to propagate the transformation onto already-loaded (or newly loaded) normal modes/potential/dipole surfaces
  - `:returns`: `Molecule`
    > the transformed molecule (with `source_file` cleared)


<a id="Psience.Molecools.Molecule.Molecule.apply_rotation" class="docs-object-method">&nbsp;</a> 
```python
apply_rotation(self, rotation_matrix, shift_com=None, load_properties=None, embed_properties=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L5485)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L5485?message=Update%20Docs)]
</div>
**LLM Docstring**

Apply a rotation (optionally shifted to rotate about the center of mass) to the molecule, via `apply_affine_transformation`.
  - `rotation_matrix`: `np.ndarray`
    > the rotation matrix (or full affine matrix) to apply
  - `shift_com`: `bool | None`
    > whether to first shift so the rotation is performed about the center of mass (only applied if `rotation_matrix` isn't already a 4x4 affine matrix)
  - `load_properties`: `bool | None`
    > forwarded to `apply_affine_transformation`
  - `embed_properties`: `bool`
    > forwarded to `apply_affine_transformation`
  - `:returns`: `Molecule`
    > the rotated molecule


<a id="Psience.Molecools.Molecule.Molecule.eckart_frame" class="docs-object-method">&nbsp;</a> 
```python
eckart_frame(self, mol, sel=None, inverse=False, planar_ref_tolerance=None, proper_rotation=False, reset_com=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L5511)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L5511?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L5538)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L5538?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L5552)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L5552?message=Update%20Docs)]
</div>
Gets the necessary data to embed crds in the Eckart frame using `self` as a reference
  - `crds`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Molecule.Molecule.get_embedded_molecule" class="docs-object-method">&nbsp;</a> 
```python
get_embedded_molecule(self, ref=None, sel=None, planar_ref_tolerance=None, proper_rotation=False, embed_properties=True, load_properties=None, reset_com=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L5564)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L5564?message=Update%20Docs)]
</div>
Returns a Molecule embedded in an Eckart frame if ref is not None, otherwise returns
a principle-axis embedded Molecule
  - `:returns`: `Molecule`
    >


<a id="Psience.Molecools.Molecule.Molecule.get_rmsd" class="docs-object-method">&nbsp;</a> 
```python
get_rmsd(self, other: "'typing.Self | np.ndarray'", sel=None, embed=True, embedding_sel=None, mass_weighted=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L5588)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L5588?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the (optionally mass-weighted) root-mean-square displacement between this molecule's geometry and another geometry/molecule, optionally Eckart-embedding both onto a common reference frame first and restricting the comparison to a subset of atoms.
  - `other`: `Molecule | np.ndarray`
    > the geometry (or molecule) to compare against
  - `sel`: `Iterable[int] | None`
    > a selection of atoms to restrict the RMSD computation to
  - `embed`: `bool`
    > whether to Eckart-embed both geometries onto a common frame before comparing
  - `embedding_sel`: `Iterable[int] | None`
    > a selection of atoms to use for the embedding fit, if different from `sel`
  - `mass_weighted`: `bool`
    > whether to mass-weight the displacement before computing the norm
  - `:returns`: `float | np.ndarray`
    > the RMSD value(s)


<a id="Psience.Molecools.Molecule.Molecule.align_molecule" class="docs-object-method">&nbsp;</a> 
```python
align_molecule(self, other: "'typing.Self'", reindex_bonds=True, permute_atoms=True, align_structures=True, sel=None, embed_properties=True, load_properties=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L5639)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L5639?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L5727)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L5727?message=Update%20Docs)]
</div>
**LLM Docstring**

The (cached) attached RDKit molecule representation, built lazily via `RDMolecule.from_mol` the first time it's needed (returning `None` if RDKit isn't available), and kept in sync with the current coordinates on subsequent access.
  - `:returns`: `object | None`
    > the RDKit molecule, or `None` if RDKit isn't available


<a id="Psience.Molecools.Molecule.Molecule.to_ase" class="docs-object-method">&nbsp;</a> 
```python
to_ase(self, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L5748)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L5748?message=Update%20Docs)]
</div>
**LLM Docstring**

Convert this molecule to an ASE (Atomic Simulation Environment) molecule object, via `ASEMolecule.from_mol`.
  - `kwargs`: `dict`
    > extra options forwarded to `ASEMolecule.from_mol`
  - `:returns`: `object`
    > the constructed ASE molecule


<a id="Psience.Molecools.Molecule.Molecule.from_ase" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_ase(cls, ase_mol, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L5762)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L5762?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a `Molecule` from an ASE molecule object.
  - `ase_mol`: `object`
    > the ASE molecule to convert
  - `kwargs`: `dict`
    > extra keyword arguments merged over the ASE molecule's metadata and forwarded to the constructor
  - `:returns`: `Molecule`
    > the constructed molecule


<a id="Psience.Molecools.Molecule.Molecule.from_zmat" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_zmat(cls, zmat, internals=None, axes=None, origin=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L5783)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L5783?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a `Molecule` from a Z-matrix specification (either a Z-matrix string or an explicit `(atoms, (ordering, coords))` tuple), converting the internal Z-matrix coordinates to Cartesians.
  - `zmat`: `str | tuple`
    > the Z-matrix specification, as a string or an `(atoms, (ordering, coords))` tuple
  - `internals`: `object | None`
    > the internal-coordinate specification to store on the resulting molecule; defaults to the Z-matrix's own atom ordering
  - `axes`: `np.ndarray | None`
    > reference axes to use when converting to Cartesians
  - `origin`: `np.ndarray | None`
    > reference origin to use when converting to Cartesians
  - `opts`: `dict`
    > extra options forwarded to the constructor
  - `:returns`: `Molecule`
    > the constructed molecule


<a id="Psience.Molecools.Molecule.Molecule.from_openbabel" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_openbabel(cls, mol, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L5813)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L5813?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L5832)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L5832?message=Update%20Docs)]
</div>
**LLM Docstring**

Convert this molecule to an OpenBabel molecule object, via `OBMolecule.from_mol`.
  - `opts`: `dict`
    > extra options forwarded to `OBMolecule.from_mol`
  - `:returns`: `object`
    > the constructed OpenBabel molecule


<a id="Psience.Molecools.Molecule.Molecule.from_rdmol" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_rdmol(cls, rdmol, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L5980)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L5980?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a `Molecule` from an RDKit molecule (or a raw RDKit `Mol`/owning-mol object, which is first wrapped in an `RDMolecule`), carrying over its atoms, coordinates, bonds, charges, and metadata.
  - `rdmol`: `object`
    > the RDKit molecule (or wrappable RDKit object) to convert
  - `opts`: `dict`
    > extra options merged over the RDKit molecule's metadata and forwarded to the constructor
  - `:returns`: `Molecule`
    > the constructed molecule


<a id="Psience.Molecools.Molecule.Molecule.from_name" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_name(cls, name, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L6513)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L6513?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a `Molecule` by looking up a compound name, via `from_string` with format `'name'`.
  - `name`: `str`
    > the compound name to look up
  - `opts`: `dict`
    > extra options forwarded to `from_string`
  - `:returns`: `Molecule | list[Molecule]`
    > the constructed molecule (or list of molecules, if multiple matches are found)


<a id="Psience.Molecools.Molecule.Molecule.get_atom_strings" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_atom_strings(cls): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L6530)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L6530?message=Update%20Docs)]
</div>
**LLM Docstring**

The (cached) set of up-to-2-character atomic element symbols known to `AtomData`, used for heuristically detecting SMILES/Z-matrix strings.
  - `:returns`: `set[str]`
    > the set of element symbol prefixes


<a id="Psience.Molecools.Molecule.Molecule.get_string_format_dispatchers" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_string_format_dispatchers(cls): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L6624)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L6624?message=Update%20Docs)]
</div>
**LLM Docstring**

The mapping from string-format key to the constructor method that parses that format, used by `from_string`.
  - `:returns`: `dict`
    > the format-to-constructor mapping


<a id="Psience.Molecools.Molecule.Molecule.from_string" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_string(cls, string, fmt=None, allow_names=False, format_options=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L6648)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L6648?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a `Molecule` from a raw string in any supported structural format, inferring the format automatically if not given, dispatching to the matching in-memory parser if one exists, or otherwise writing the string to a temporary file and dispatching through `from_file`.
  - `string`: `str`
    > the structural string to parse
  - `fmt`: `str | None`
    > the format to parse as; inferred via `_infer_str_format` if not given
  - `allow_names`: `bool`
    > whether to allow inferring the format as a compound name
  - `format_options`: `dict | None`
    > per-format parsing options, keyed by format, merged under `opts`
  - `opts`: `dict`
    > extra options forwarded to the format-specific parser
  - `:returns`: `Molecule`
    > the constructed molecule


<a id="Psience.Molecools.Molecule.Molecule.get_file_format_dispatchers" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_file_format_dispatchers(cls): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L6705)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L6705?message=Update%20Docs)]
</div>
**LLM Docstring**

The mapping from file-format key (typically a file extension) to the constructor method that parses that format, used by `from_file`.
  - `:returns`: `dict`
    > the format-to-constructor mapping


<a id="Psience.Molecools.Molecule.Molecule.from_file" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_file(cls, file, mode=None, format_options=None, use_ob_fallback=False, **opts) -> 'Molecule': 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L6727)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L6727?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L6969)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L6969?message=Update%20Docs)]
</div>
**LLM Docstring**

The mapping from string-export format key to the exporter method that produces that format, used by `to_string`.
  - `:returns`: `dict`
    > the format-to-exporter mapping


<a id="Psience.Molecools.Molecule.Molecule.to_string" class="docs-object-method">&nbsp;</a> 
```python
to_string(self, fmt, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L6989)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L6989?message=Update%20Docs)]
</div>
**LLM Docstring**

Export this molecule to a string in the given format, dispatching to an in-memory string exporter if one exists, otherwise round-tripping through a temporary file via `to_file`, or falling back to OpenBabel's string-export machinery.
  - `fmt`: `str`
    > the export format key
  - `opts`: `dict`
    > extra options forwarded to the format-specific exporter
  - `:returns`: `str`
    > the exported string


<a id="Psience.Molecools.Molecule.Molecule.get_file_export_dispatchers" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_file_export_dispatchers(cls): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L7027)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L7027?message=Update%20Docs)]
</div>
**LLM Docstring**

The mapping from file-export format key to the exporter method that writes that format to disk, used by `to_file`. Currently empty (no dedicated file exporters beyond the string-based ones and OpenBabel's fallback).
  - `:returns`: `dict`
    > the (currently empty) format-to-exporter mapping


<a id="Psience.Molecools.Molecule.Molecule.to_file" class="docs-object-method">&nbsp;</a> 
```python
to_file(self, file, mode=None, use_ob_fallback=False, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L7039)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L7039?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L7106)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L7106?message=Update%20Docs)]
</div>
**LLM Docstring**

Universal `Molecule` constructor: builds a molecule from essentially any reasonable input (an existing `Molecule` to copy/modify, an RDKit/ASE object, a file path, a raw structural string, a dict of constructor kwargs, or an `(atoms, coords)`/Z-matrix pair), inferring the format automatically if not given.
  - `spec`: `object`
    > the specification to build the molecule from
  - `fmt`: `str | tuple | None`
    > an explicit format to use instead of inferring one; can also be an `(atoms, coords)` pair for direct construction
  - `opts`: `dict`
    > extra options forwarded to the resolved constructor
  - `:returns`: `Molecule`
    > the constructed molecule


<a id="Psience.Molecools.Molecule.Molecule.plot" class="docs-object-method">&nbsp;</a> 
```python
plot(self, *geometries, figure=None, return_objects=False, bonds=None, bond_radius=None, atom_radius_scaling=None, atom_style=None, atom_radii=None, atom_text=None, display_atom_numbers=False, radius_type=None, bond_style=None, reconcile_bonds=True, capped_bonds=None, reflectiveness=None, vector_style=None, highlight_atoms=None, highlight_bonds=None, highlight_rings=None, highlight_styles=None, comparison_styles=None, animation_frame_styles=None, mode_vectors=None, mode_vector_origins=None, mode_vector_origin_mode='set', mode_vector_display_cutoff=0.01, principle_axes=None, principle_axes_origin=None, principle_axes_origin_mode='set', principle_axes_style=None, dipole=None, dipole_origin=None, dipole_origin_mode='set', render_multiple_bonds=None, render_fractional_bonds=None, fractional_bond_offset=None, bond_center_radius_offset=None, draw_coords=None, draw_coords_style=None, up_vector=None, multiple_bond_spacing=None, mode=None, backend=None, include_save_buttons=None, objects=False, graphics_class=None, cylinder_class=None, cylinder_options=None, sphere_class=None, sphere_options=None, arrow_class=None, arrow_options=None, line_class=None, line_options=None, disk_class=None, disk_options=None, animate=None, recording_options=None, animation_options=None, jsmol_load_script=None, include_jsmol_script_interface=None, dynamic_loading=None, units='Angstroms', label_style=None, theme='default', theme_function=None, plot_range_padding='auto', annotation_function=None, **plot_ops): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L7158)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L7158?message=Update%20Docs)]
</div>
Dispatches to the appropriate `MoleculePlotter` for the resolved backend/mode.

The full keyword surface now lives on `MoleculePlotter.plot_molecule`; calls,
defaults, and return conventions are unchanged.


<a id="Psience.Molecools.Molecule.Molecule.get_animation_geoms" class="docs-object-method">&nbsp;</a> 
```python
get_animation_geoms(self, which, extent=0.35, steps=8, strip_embedding=True, units=None, coordinate_expansion=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L7306)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L7306?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a back-and-forth looping sequence of displaced geometries for animating a single coordinate (or an arbitrary coordinate-expansion direction), via `get_scan_coordinates`.
  - `which`: `int | np.ndarray | list[np.ndarray]`
    > the coordinate index to animate (if `coordinate_expansion` is not given directly), or an explicit displacement-direction array/list of arrays (in which case this becomes the coordinate index within that expansion, defaulting to `0`)
  - `extent`: `float`
    > the maximum displacement magnitude in each direction
  - `steps`: `int`
    > number of steps from the equilibrium geometry out to `extent`
  - `strip_embedding`: `bool`
    > whether to strip the fixed embedding coordinates from the default Cartesians-by-internals expansion
  - `units`: `str | None`
    > units to convert the resulting geometries into; left in Bohr if `None`
  - `coordinate_expansion`: `list[np.ndarray] | np.ndarray | None`
    > an explicit coordinate-transformation expansion to displace along, instead of the default internal-coordinate Jacobian
  - `:returns`: `np.ndarray`
    > the looping sequence of displaced geometries (out to `extent` and back)


<a id="Psience.Molecools.Molecule.Molecule.animate_coordinate" class="docs-object-method">&nbsp;</a> 
```python
animate_coordinate(self, which, extent=0.5, steps=3, return_objects=False, strip_embedding=True, units='Angstroms', backend=None, mode=None, jsmol_load_script=None, coordinate_expansion=None, **plot_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L7347)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L7347?message=Update%20Docs)]
</div>
**LLM Docstring**

Build an animation of a displaced coordinate and return it as a plottable/displayable object, either as a JSMol vibration animation or as a sequence of rendered frames via `plot`.
  - `which`: `int | np.ndarray | list[np.ndarray]`
    > the coordinate index (or explicit direction) to animate, forwarded to `get_animation_geoms`/`format_animation_file`
  - `extent`: `float`
    > the maximum displacement magnitude
  - `steps`: `int`
    > number of steps out to `extent`
  - `return_objects`: `bool`
    > whether to return the constructed graphics objects alongside the figure (non-JSMol path only)
  - `strip_embedding`: `bool`
    > whether to strip the fixed embedding coordinates from the default coordinate expansion
  - `units`: `str`
    > units to display the geometries in
  - `backend`: `str | None`
    > the rendering backend to use; defaults to `self.display_mode`
  - `mode`: `str | None`
    > the display mode to use; defaults to `backend`
  - `jsmol_load_script`: `str | list[str] | None`
    > extra JSMol load-script text, used only in the JSMol path
  - `coordinate_expansion`: `list[np.ndarray] | np.ndarray | None`
    > an explicit coordinate-transformation expansion to animate along
  - `plot_opts`: `dict`
    > extra options forwarded to `plot`
  - `:returns`: `object`
    > the resulting animation figure/widget


<a id="Psience.Molecools.Molecule.Molecule.animate_mode" class="docs-object-method">&nbsp;</a> 
```python
animate_mode(self, which, extent=0.5, steps=3, modes=None, coordinate_expansion=None, order=None, normalize=True, mass_weight=False, mass_scale=True, frequency_scale=False, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L7407)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L7407?message=Update%20Docs)]
</div>
**LLM Docstring**

Build an animation of a normal mode's displacement, converting the mode into a coordinate-expansion direction (with optional normalization, mass-weighting/scaling, and frequency scaling of the displacement extent) and delegating to `animate_coordinate`.
  - `which`: `int`
    > the mode index to animate
  - `extent`: `float`
    > the base displacement extent, before any mass/frequency scaling
  - `steps`: `int`
    > number of steps out to `extent`
  - `modes`: `object | None`
    > the normal modes to animate; defaults to this molecule's own (converted to a fresh mode basis)
  - `coordinate_expansion`: `list[np.ndarray] | None`
    > an additional coordinate-transformation expansion to combine with the mode-displacement direction
  - `order`: `int | None`
    > the Jacobian order to use when building the mode-based coordinate expansion (only used if given)
  - `normalize`: `bool`
    > whether to normalize the mode-displacement direction before use
  - `mass_weight`: `bool`
    > whether to keep the modes mass-weighted rather than removing the mass-weighting
  - `mass_scale`: `bool`
    > whether to scale the displacement extent by the mode's effective mass (only applied when not mass-weighted)
  - `frequency_scale`: `bool`
    > whether to scale the displacement extent relative to the mode's frequency
  - `opts`: `dict`
    > extra options forwarded to `animate_coordinate`
  - `:returns`: `object`
    > the resulting animation figure/widget


<a id="Psience.Molecools.Molecule.Molecule.format_structs" class="docs-object-method">&nbsp;</a> 
```python
format_structs(self, geoms, format='xyz'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L7505)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L7505?message=Update%20Docs)]
</div>
**LLM Docstring**

Format a batch of geometries as a concatenated multi-frame string in the given format.
  - `geoms`: `np.ndarray`
    > the geometries to format, reshaped to `(-1, natoms, 3)`
  - `format`: `str`
    > the output format; currently only `'xyz'` is supported
  - `:returns`: `str`
    > the concatenated multi-frame string


<a id="Psience.Molecools.Molecule.Molecule.format_animation_file" class="docs-object-method">&nbsp;</a> 
```python
format_animation_file(self, which, format='xyz', extent=0.35, steps=8, strip_embedding=True, units='Angstroms', coordinate_expansion=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L7564)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L7564?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a formatted animation string/block for a displaced coordinate, either as a single JSMol vibration block (`'jmol'` format) or as a looping multi-frame structure string via `get_animation_geoms`/`format_structs`.
  - `which`: `int | np.ndarray | list[np.ndarray]`
    > the coordinate index (or explicit direction) to animate
  - `format`: `str`
    > the output format (`'jmol'` for a JSMol vibration block, otherwise forwarded to `format_structs`)
  - `extent`: `float`
    > the maximum displacement magnitude
  - `steps`: `int`
    > number of steps out to `extent` (non-`'jmol'` formats only)
  - `strip_embedding`: `bool`
    > whether to strip the fixed embedding coordinates from the default coordinate expansion
  - `units`: `str`
    > units to format the geometries in
  - `coordinate_expansion`: `list[np.ndarray] | None`
    > an explicit coordinate-transformation expansion to animate along
  - `:returns`: `str`
    > the formatted animation string


<a id="Psience.Molecools.Molecule.Molecule.to_widget" class="docs-object-method">&nbsp;</a> 
```python
to_widget(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Molecule/Molecule.py#L7604)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule/Molecule.py#L7604?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a displayable widget for this molecule, either a JSMol applet or an X3D scene, depending on `self.display_mode`.
  - `:returns`: `object`
    > the constructed widget/scene object
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
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Molecule.py#L51?message=Update%20Docs)   
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