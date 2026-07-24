## <a id="Psience.Molecools.Properties.PotentialSurfaceManager">PotentialSurfaceManager</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties.py#L2375)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties.py#L2375?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
name: str
```
<a id="Psience.Molecools.Properties.PotentialSurfaceManager.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, mol, surface=None, derivatives=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties.py#L2377)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties.py#L2377?message=Update%20Docs)]
</div>
**LLM Docstring**

Set up the potential-surface manager, either copying the surface/derivatives/surface-coordinates from an existing surface-like object, or storing the supplied `surface`/`derivatives` directly.
  - `mol`: `AbstractMolecule`
    > the molecule this manager is attached to
  - `surface`: `object | None`
    > an existing potential surface, or another object exposing `_surf`/`_derivs`/`_surface_coords` to copy from
  - `derivatives`: `list[np.ndarray] | None`
    > raw potential derivative tensors
  - `:returns`: `None`
    > None


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.from_data" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_data(cls, mol, data): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L2403)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L2403?message=Update%20Docs)]
</div>
**LLM Docstring**

Abstract constructor for building a potential-surface manager directly from raw data. Not implemented.
  - `mol`: `AbstractMolecule`
    > the molecule the manager will be attached to
  - `data`: `object`
    > the raw data to build from
  - `:returns`: `PotentialSurfaceManager`
    > never returns


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.surface" class="docs-object-method">&nbsp;</a> 
```python
@property
surface(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2420)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2420?message=Update%20Docs)]
</div>
**LLM Docstring**

The potential surface object, lazily loaded via `load_potential_surface(self.surface_coords)` if not already set.
  - `:returns`: `object`
    > the potential surface


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.surface_coords" class="docs-object-method">&nbsp;</a> 
```python
@property
surface_coords(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2433)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2433?message=Update%20Docs)]
</div>
**LLM Docstring**

Property getter/setter for the coordinate specification used when loading the potential surface (e.g. bond/angle/dihedral tuples).
  - `coords`: `object`
    > (setter only) the new coordinate specification
  - `:returns`: `object`
    > (getter) the stored coordinate specification


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.get_derivs" class="docs-object-method">&nbsp;</a> 
```python
get_derivs(self, quiet=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2460)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2460?message=Update%20Docs)]
</div>
**LLM Docstring**

Fetch the potential-energy derivative tensors, lazily loading them via `load_potential_derivatives` if not already cached.
  - `quiet`: `bool`
    > if `True`, suppresses errors from the underlying loader when the derivatives can't be found
  - `:returns`: `tuple`
    > the potential derivative tensors


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.derivatives" class="docs-object-method">&nbsp;</a> 
```python
@property
derivatives(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2474)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2474?message=Update%20Docs)]
</div>
**LLM Docstring**

Property getter/setter for the potential-energy derivative tensors. The getter delegates to `get_derivs()`.
  - `v`: `tuple`
    > (setter only) the new derivative tensors to store
  - `:returns`: `tuple`
    > (getter) the potential derivative tensors


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.force_constants" class="docs-object-method">&nbsp;</a> 
```python
@property
force_constants(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2501)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2501?message=Update%20Docs)]
</div>
**LLM Docstring**

The force-constant (second-derivative/Hessian) tensor from the potential derivative expansion.
  - `:returns`: `np.ndarray`
    > `self.derivatives[1]`


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.load_potential_derivatives" class="docs-object-method">&nbsp;</a> 
```python
load_potential_derivatives(self, file=None, mode=None, quiet=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2513)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2513?message=Update%20Docs)]
</div>
Loads potential derivatives from a file (or from `source_file` if set)
  - `file`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.load_potential_surface" class="docs-object-method">&nbsp;</a> 
```python
load_potential_surface(self, coordinates): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2611)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2611?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a `PotentialSurface` object by loading a Gaussian scan log file (`self.mol.source_file`) and, if `coordinates` gives an iterable of atom-index tuples (bond/angle/dihedral specs), constructing the coordinate-transform function used to interpret the scan coordinates; otherwise `coordinates` itself is treated as that transform function.
  - `coordinates`: `callable | Iterable[tuple]`
    > either a callable that maps Cartesian coordinates to scan coordinates, or an iterable of 2/3/4-atom index tuples specifying bond lengths/angles/dihedrals to use as scan coordinates
  - `:returns`: `PotentialSurface`
    > the loaded potential surface


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.load" class="docs-object-method">&nbsp;</a> 
```python
load(self, coordinates=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2661)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2661?message=Update%20Docs)]
</div>
**LLM Docstring**

Load the potential surface using `coordinates` if one is given and no surface is cached yet, then return whichever representation (surface or derivatives) is available.
  - `coordinates`: `object | None`
    > coordinate specification forwarded to `load_potential_surface` if the surface still needs to be loaded
  - `:returns`: `object`
    > the potential surface or its derivatives


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.update" class="docs-object-method">&nbsp;</a> 
```python
update(self, val): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2678)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2678?message=Update%20Docs)]
</div>
Updates the held values
  - `val`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.apply_transformation" class="docs-object-method">&nbsp;</a> 
```python
apply_transformation(self, transf): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2689)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2689?message=Update%20Docs)]
</div>
**LLM Docstring**

Apply a coordinate transformation to a copy of this potential-surface manager, transforming the surface object (if set) and/or the derivative tensors (if set) via either a direct linear transform or a general `TensorDerivativeConverter` re-expansion.
  - `transf`: `object | np.ndarray | list`
    > the transformation to apply, either an object exposing `transformation_function`, a `3x3` transformation matrix, or a general derivative-tensor sequence
  - `:returns`: `PotentialSurfaceManager`
    > a new `PotentialSurfaceManager` with the transformation applied


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.insert_atoms" class="docs-object-method">&nbsp;</a> 
```python
insert_atoms(self, atoms, coords, where): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2717)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2717?message=Update%20Docs)]
</div>
Handles the insertion of new atoms into the structure
  - `atoms`: `tuple[str]`
    > 
  - `coords`: `CoordinateSet`
    > 
  - `where`: `tuple[int]`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.PotentialSurfaceManager.delete_atoms" class="docs-object-method">&nbsp;</a> 
```python
delete_atoms(self, where): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2736)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/PotentialSurfaceManager.py#L2736?message=Update%20Docs)]
</div>
Handles the deletion from the structure
  - `atoms`: `tuple[str]`
    > 
  - `coords`: `CoordinateSet`
    > 
  - `where`: `tuple[int]`
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Molecools/Properties/PotentialSurfaceManager.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Molecools/Properties/PotentialSurfaceManager.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Molecools/Properties/PotentialSurfaceManager.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Molecools/Properties/PotentialSurfaceManager.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties.py#L2375?message=Update%20Docs)   
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