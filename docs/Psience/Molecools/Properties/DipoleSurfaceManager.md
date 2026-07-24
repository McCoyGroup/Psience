## <a id="Psience.Molecools.Properties.DipoleSurfaceManager">DipoleSurfaceManager</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties.py#L1934)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties.py#L1934?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
name: str
```
<a id="Psience.Molecools.Properties.DipoleSurfaceManager.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, mol, surface=None, derivatives=None, polarizability_derivatives=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties.py#L1936)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties.py#L1936?message=Update%20Docs)]
</div>
**LLM Docstring**

Set up the dipole-surface manager, either copying the surface/derivatives/analytic-flag from an existing surface-like object, or storing the supplied `surface`/`derivatives`/`polarizability_derivatives` directly (unpacking `derivatives` into numerical and analytic parts if given as a dict).
  - `mol`: `AbstractMolecule`
    > the molecule this manager is attached to
  - `surface`: `object | None`
    > an existing dipole surface, or another object exposing `_surf`/`_derivs`/`_analytic_derivatives` to copy from
  - `derivatives`: `list[np.ndarray] | dict | None`
    > raw dipole derivative tensors, or a dict with `'numerical'`/`'analytic'` keys
  - `polarizability_derivatives`: `list[np.ndarray] | None`
    > raw polarizability derivative tensors
  - `:returns`: `None`
    > None


<a id="Psience.Molecools.Properties.DipoleSurfaceManager.from_data" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_data(cls, mol, data): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1969)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1969?message=Update%20Docs)]
</div>
**LLM Docstring**

Abstract constructor for building a dipole-surface manager directly from raw data. Not implemented.
  - `mol`: `AbstractMolecule`
    > the molecule the manager will be attached to
  - `data`: `object`
    > the raw data to build from
  - `:returns`: `DipoleSurfaceManager`
    > never returns


<a id="Psience.Molecools.Properties.DipoleSurfaceManager.surface" class="docs-object-method">&nbsp;</a> 
```python
@property
surface(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L1986)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L1986?message=Update%20Docs)]
</div>
**LLM Docstring**

The dipole surface object, lazily loaded via `load_dipole_surface` if not already set.
  - `:returns`: `object`
    > the dipole surface


<a id="Psience.Molecools.Properties.DipoleSurfaceManager.numerical_derivatives" class="docs-object-method">&nbsp;</a> 
```python
@property
numerical_derivatives(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L1999)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L1999?message=Update%20Docs)]
</div>
**LLM Docstring**

The numerically-computed dipole derivatives, lazily loaded (together with the analytic derivatives, if provided by the same source) via `load_dipole_derivatives` if neither is already cached.
  - `:returns`: `tuple | None`
    > the numerical dipole derivatives, or `None` if unavailable


<a id="Psience.Molecools.Properties.DipoleSurfaceManager.get_derivatives" class="docs-object-method">&nbsp;</a> 
```python
get_derivatives(self, quiet=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L2018)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L2018?message=Update%20Docs)]
</div>
**LLM Docstring**

Fetch the analytic dipole derivatives, lazily loading them (together with the numerical derivatives, if bundled) via `load_dipole_derivatives` if not already cached.
  - `quiet`: `bool`
    > if `True`, suppresses errors from the underlying loader when the derivatives can't be found
  - `:returns`: `tuple`
    > the analytic dipole derivatives


<a id="Psience.Molecools.Properties.DipoleSurfaceManager.derivatives" class="docs-object-method">&nbsp;</a> 
```python
@property
derivatives(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L2038)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L2038?message=Update%20Docs)]
</div>
**LLM Docstring**

Property getter/setter for the analytic dipole derivatives. The getter delegates to `get_derivatives()`; the setter accepts either a raw derivative tuple or a dict with `'numerical'`/`'analytic'` keys.
  - `derivatives`: `tuple | dict`
    > (setter only) the new derivatives to store
  - `:returns`: `tuple`
    > (getter) the analytic dipole derivatives


<a id="Psience.Molecools.Properties.DipoleSurfaceManager.get_pol_derivatives" class="docs-object-method">&nbsp;</a> 
```python
get_pol_derivatives(self, quiet=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L2070)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L2070?message=Update%20Docs)]
</div>
**LLM Docstring**

Fetch the polarizability derivatives, lazily loading them via `load_polarizability_derivatives` if not already cached.
  - `quiet`: `bool`
    > if `True`, suppresses errors from the underlying loader when the derivatives can't be found
  - `:returns`: `list`
    > the polarizability derivatives


<a id="Psience.Molecools.Properties.DipoleSurfaceManager.polarizability_derivatives" class="docs-object-method">&nbsp;</a> 
```python
@property
polarizability_derivatives(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L2085)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L2085?message=Update%20Docs)]
</div>
**LLM Docstring**

Property getter/setter for the polarizability derivatives. The getter delegates to `get_pol_derivatives()`.
  - `derivatives`: `list`
    > (setter only) the new polarizability derivatives to store
  - `:returns`: `list`
    > (getter) the polarizability derivatives


<a id="Psience.Molecools.Properties.DipoleSurfaceManager.load" class="docs-object-method">&nbsp;</a> 
```python
load(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L2112)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L2112?message=Update%20Docs)]
</div>
**LLM Docstring**

Return whichever representation of the dipole is available: the surface object if one is set, otherwise the derivative expansion.
  - `:returns`: `object`
    > the dipole surface or its derivatives


<a id="Psience.Molecools.Properties.DipoleSurfaceManager.update" class="docs-object-method">&nbsp;</a> 
```python
update(self, val): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L2125)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L2125?message=Update%20Docs)]
</div>
Updates the held values
  - `val`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.DipoleSurfaceManager.load_dipole_surface" class="docs-object-method">&nbsp;</a> 
```python
load_dipole_surface(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L2136)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L2136?message=Update%20Docs)]
</div>
**LLM Docstring**

Not currently implemented: general (non-derivative-expansion) dipole surfaces aren't supported yet.
  - `:returns`: `object`
    > never returns


<a id="Psience.Molecools.Properties.DipoleSurfaceManager.load_dipole_derivatives" class="docs-object-method">&nbsp;</a> 
```python
load_dipole_derivatives(self, file=None, quiet=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L2218)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L2218?message=Update%20Docs)]
</div>
Loads dipole derivatives from a file (or from `source_file` if set)
  - `file`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.DipoleSurfaceManager.load_polarizability_derivatives" class="docs-object-method">&nbsp;</a> 
```python
load_polarizability_derivatives(self, file=None, quiet=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L2263)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L2263?message=Update%20Docs)]
</div>
Loads dipole derivatives from a file (or from `source_file` if set)
  - `file`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.DipoleSurfaceManager.apply_transformation" class="docs-object-method">&nbsp;</a> 
```python
apply_transformation(self, transf): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L2299)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L2299?message=Update%20Docs)]
</div>
**LLM Docstring**

Apply a coordinate transformation to a copy of this dipole-surface manager, transforming the surface object (if set) and/or the derivative tensors (if set); non-linear transformations of the derivatives are not yet supported.
  - `transf`: `object | np.ndarray`
    > the transformation to apply, either an object exposing `transformation_function` or a `3x3` transformation matrix
  - `:returns`: `DipoleSurfaceManager`
    > a new `DipoleSurfaceManager` with the transformation applied


<a id="Psience.Molecools.Properties.DipoleSurfaceManager.insert_atoms" class="docs-object-method">&nbsp;</a> 
```python
insert_atoms(self, atoms, coords, where): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L2337)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L2337?message=Update%20Docs)]
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


<a id="Psience.Molecools.Properties.DipoleSurfaceManager.delete_atoms" class="docs-object-method">&nbsp;</a> 
```python
delete_atoms(self, where): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L2360)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/DipoleSurfaceManager.py#L2360?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Molecools/Properties/DipoleSurfaceManager.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Molecools/Properties/DipoleSurfaceManager.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Molecools/Properties/DipoleSurfaceManager.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Molecools/Properties/DipoleSurfaceManager.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties.py#L1934?message=Update%20Docs)   
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