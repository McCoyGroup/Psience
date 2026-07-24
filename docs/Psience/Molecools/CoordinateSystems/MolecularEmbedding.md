## <a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding">MolecularEmbedding</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems.py#L32)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems.py#L32?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
internal_fd_defaults: dict
cart_fd_defaults: dict
cartesian_by_internals_method: str
```
<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, masses, coords, internals, internal_fd_opts=None, cartesian_fd_opts=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems.py#L34)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems.py#L34?message=Update%20Docs)]
</div>
**LLM Docstring**

Set up a molecule's coordinate embedding: wraps the Cartesian coordinates in a `MolecularCartesianCoordinateSystem`, canonicalizes the internal-coordinate specification (or stores it directly if already a `CoordinateSet`), and initializes the Jacobian cache and finite-difference option overrides.
  - `masses`: `np.ndarray`
    > the atomic masses
  - `coords`: `np.ndarray`
    > the Cartesian coordinates
  - `internals`: `object`
    > the internal-coordinate specification (Z-matrix, generic-internal specs, a callable conversion, or an already-built `CoordinateSet`)
  - `internal_fd_opts`: `dict | None`
    > overrides for the internal-coordinate finite-difference defaults
  - `cartesian_fd_opts`: `dict | None`
    > overrides for the Cartesian finite-difference defaults
  - `:returns`: `None`
    > None


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.get_direct_converter" class="docs-object-method">&nbsp;</a> 
```python
get_direct_converter(self, target): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L75)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L75?message=Update%20Docs)]
</div>
**LLM Docstring**

Provide a converter from this embedding's coordinate system directly to plain (non-molecular) 3D Cartesian coordinates, if `target` is Cartesian-compatible.
  - `target`: `object`
    > the coordinate system being converted to
  - `:returns`: `MolecularCartesianToRegularCartesianConverter | None`
    > a `MolecularCartesianToRegularCartesianConverter`, or `None` if `target` isn't Cartesian-compatible


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.get_inverse_converter" class="docs-object-method">&nbsp;</a> 
```python
get_inverse_converter(self, target): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L89)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L89?message=Update%20Docs)]
</div>
**LLM Docstring**

Provide a converter from plain (non-molecular) 3D Cartesian coordinates into this embedding's coordinate system, if `target` is Cartesian-compatible.
  - `target`: `object`
    > the coordinate system being converted from
  - `:returns`: `RegularCartesianToMolecularCartesianConverter | None`
    > a `RegularCartesianToMolecularCartesianConverter`, or `None` if `target` isn't Cartesian-compatible


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.__del__" class="docs-object-method">&nbsp;</a> 
```python
__del__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L104)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L104?message=Update%20Docs)]
</div>
**LLM Docstring**

Deregister any coordinate converters this embedding registered, via `cleanup`, when the object is garbage collected.
  - `:returns`: `None`
    > None


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.cleanup" class="docs-object-method">&nbsp;</a> 
```python
cleanup(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L115)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L115?message=Update%20Docs)]
</div>
**LLM Docstring**

Deregister every converter previously registered via `register`, if any.
  - `:returns`: `None`
    > None


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.register" class="docs-object-method">&nbsp;</a> 
```python
register(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L128)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L128?message=Update%20Docs)]
</div>
**LLM Docstring**

Register the Cartesian coordinate converters for this embedding's coordinate system with the global converter registry, if not already registered.
  - `:returns`: `None`
    > None


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.coords" class="docs-object-method">&nbsp;</a> 
```python
@property
coords(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L144)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L144?message=Update%20Docs)]
</div>
**LLM Docstring**

Property getter/setter for the Cartesian coordinates. The getter registers the coordinate converters (via `register`) before returning them. The setter accepts a raw array or an already-systemed `CoordinateSet`, invalidates the Jacobian cache and inertial-frame cache, and marks the converters as needing re-registration if the coordinate system changed.
  - `coords`: `np.ndarray | CoordinateSet`
    > (setter only) the new Cartesian coordinates
  - `:returns`: `CoordinateSet`
    > (getter) the Cartesian coordinates


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.masses" class="docs-object-method">&nbsp;</a> 
```python
@property
masses(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L180)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L180?message=Update%20Docs)]
</div>
**LLM Docstring**

The atomic masses, taken from the Cartesian coordinate system.
  - `:returns`: `np.ndarray`
    > the atomic masses


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.canonicalize_internal_coordinate_spec" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
canonicalize_internal_coordinate_spec(cls, spec): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L230)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L230?message=Update%20Docs)]
</div>
**LLM Docstring**

Normalize the many accepted forms of an internal-coordinate specification (an options dict with `'zmatrix'`/`'specs'`/`'conversion'` keys, a bare callable conversion function, a raw Z-matrix-like array, or a list of generic-internal-coordinate specs) into the single canonical dict form (`'specs'`, `'zmatrix'`, `'conversion'`, `'inverse'`, `'converter_options'`) used internally, wrapping any conversion callables via `_wrap_conv` and filling in default embedding/jacobian-prep converter options for Z-matrices.
  - `spec`: `dict | callable | Iterable | None`
    > the internal-coordinate specification to canonicalize
  - `:returns`: `dict | None`
    > the canonicalized specification dict, or `None`/the original value if `spec` is `None`


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.internals" class="docs-object-method">&nbsp;</a> 
```python
@property
internals(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L291)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L291?message=Update%20Docs)]
</div>
**LLM Docstring**

Property getter/setter for the raw (canonicalized) internal-coordinate specification. The setter re-canonicalizes the given specification and invalidates any already-computed internal coordinates.
  - `internals`: `object`
    > (setter only) the new internal-coordinate specification, in any form accepted by `canonicalize_internal_coordinate_spec`
  - `:returns`: `dict | None`
    > (getter) the canonicalized specification dict, or `None` if none is set


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.zmatrix" class="docs-object-method">&nbsp;</a> 
```python
@property
zmatrix(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L321)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L321?message=Update%20Docs)]
</div>
**LLM Docstring**

Property getter/setter for just the Z-matrix ordering array out of the internal-coordinate specification. The setter validates the Z-matrix shape, builds a fresh specification if none exists yet, and invalidates any already-computed internal coordinates.
  - `zmat`: `np.ndarray | None`
    > (setter only) the new Z-matrix ordering array
  - `:returns`: `np.ndarray | None`
    > (getter) the stored Z-matrix array, or `None`


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.convert_to_internals" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
convert_to_internals(cls, coords, masses, spec): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L360)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L360?message=Update%20Docs)]
</div>
**LLM Docstring**

Build the internal-coordinate `CoordinateSet` described by `spec`: constructs (and registers converters for) a generic-internal, Z-matrix, or iterative-Z-matrix coordinate system as appropriate, converts `coords` into it, layers on any extra custom `conversion`/`inverse` via a `CompositeCoordinateSystem`, and returns the resulting coordinates together with the (possibly updated, e.g. with redundant-transformation info) spec.
  - `coords`: `CoordinateSet`
    > the Cartesian coordinates to convert
  - `masses`: `np.ndarray`
    > the atomic masses
  - `spec`: `dict`
    > the canonicalized internal-coordinate specification (as produced by `canonicalize_internal_coordinate_spec`)
  - `:returns`: `tuple[CoordinateSet, dict]`
    > `(coords, spec)` -- the internal coordinates and the (possibly updated) specification


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.internal_coordinates_from_spec" class="docs-object-method">&nbsp;</a> 
```python
internal_coordinates_from_spec(self, spec: dict): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L437)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L437?message=Update%20Docs)]
</div>
**LLM Docstring**

Build the internal coordinates for this embedding's current Cartesian coordinates and masses from a given specification, via `convert_to_internals`.
  - `spec`: `dict`
    > the canonicalized internal-coordinate specification
  - `:returns`: `tuple[CoordinateSet, dict]`
    > `(coords, spec)`, as returned by `convert_to_internals`


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.internal_coordinates" class="docs-object-method">&nbsp;</a> 
```python
@property
internal_coordinates(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L454)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L454?message=Update%20Docs)]
</div>
**LLM Docstring**

Property getter/setter for the internal coordinates. The getter lazily computes them from the stored specification (via `internal_coordinates_from_spec`) the first time they're needed. The setter requires an already-built `CoordinateSet`.
  - `ics`: `CoordinateSet`
    > (setter only) the new internal coordinates
  - `:returns`: `CoordinateSet | None`
    > (getter) the internal coordinates, or `None` if no internal-coordinate specification is set


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.strip_embedding_coordinates" class="docs-object-method">&nbsp;</a> 
```python
strip_embedding_coordinates(self, coords): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L493)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L493?message=Update%20Docs)]
</div>
**LLM Docstring**

Drop the fixed embedding coordinates (e.g. the 6 translation/rotation degrees of freedom implied by the Z-matrix embedding) from a coordinate array or list of derivative tensors, if the underlying internal-coordinate system defines any.
  - `coords`: `np.ndarray | list[np.ndarray]`
    > the coordinates (or list of derivative tensors) to strip
  - `:returns`: `np.ndarray | list[np.ndarray]`
    > the coordinates with embedding coordinates removed, or unchanged if there are none to strip or they're already stripped


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.strip_derivative_embedding" class="docs-object-method">&nbsp;</a> 
```python
strip_derivative_embedding(self, derivs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L521)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L521?message=Update%20Docs)]
</div>
**LLM Docstring**

Drop the fixed embedding coordinates from every axis of each tensor in a list of Cartesian-derivative tensors, if the underlying internal-coordinate system defines any.
  - `derivs`: `list[np.ndarray]`
    > the list of Cartesian-derivative tensors (order-`n` tensor at index `n-1`) to strip
  - `:returns`: `list[np.ndarray]`
    > the derivative tensors with embedding coordinates removed from every relevant axis, or unchanged if there are none to strip or they're already stripped


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.restore_embedding_coordinates" class="docs-object-method">&nbsp;</a> 
```python
restore_embedding_coordinates(self, coords): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L544)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L544?message=Update%20Docs)]
</div>
**LLM Docstring**

Reinsert the fixed embedding coordinates (filled in from the reference internal coordinates) back into a stripped coordinate array or list, undoing `strip_embedding_coordinates`.
  - `coords`: `np.ndarray | list[np.ndarray]`
    > the stripped coordinates (or list) to restore
  - `:returns`: `np.ndarray | list[np.ndarray]`
    > the coordinates with embedding coordinates reinserted, or unchanged if there are none to restore or they're already present


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.restore_derivative_embedding" class="docs-object-method">&nbsp;</a> 
```python
restore_derivative_embedding(self, derivs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L585)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L585?message=Update%20Docs)]
</div>
**LLM Docstring**

Reinsert zeroed-out placeholder entries for the fixed embedding coordinates back into every axis of each tensor in a list of stripped Cartesian-derivative tensors, undoing `strip_derivative_embedding`.
  - `derivs`: `list[np.ndarray]`
    > the stripped list of Cartesian-derivative tensors to restore
  - `:returns`: `list[np.ndarray]`
    > the derivative tensors with embedding-coordinate axes reinserted (as zeros), or unchanged if there are none to restore or they're already present


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.get_internals" class="docs-object-method">&nbsp;</a> 
```python
get_internals(self, *, coords=None, strip_embedding=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L614)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L614?message=Update%20Docs)]
</div>
**LLM Docstring**

Fetch the internal coordinates, either the cached ones for this embedding's own geometry or freshly computed ones for an alternate set of Cartesian `coords`, optionally stripping the fixed embedding coordinates.
  - `coords`: `np.ndarray | None`
    > alternate Cartesian coordinates to convert instead of using the cached internal coordinates
  - `strip_embedding`: `bool`
    > whether to strip the fixed embedding coordinates from the result
  - `:returns`: `CoordinateSet | None`
    > the internal coordinates, or `None` if none are defined


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.get_cartesians" class="docs-object-method">&nbsp;</a> 
```python
get_cartesians(self, *, coords=None, strip_embedding=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L637)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L637?message=Update%20Docs)]
</div>
**LLM Docstring**

Fetch Cartesian coordinates, either this embedding's own cached ones or the Cartesian coordinates corresponding to a given set of internal `coords`, optionally restoring any stripped embedding coordinates first.
  - `coords`: `np.ndarray | None`
    > internal-coordinate values to convert to Cartesians instead of returning the cached Cartesian coordinates
  - `strip_embedding`: `bool`
    > whether `coords` has had its embedding coordinates stripped and needs them restored before conversion
  - `:returns`: `CoordinateSet`
    > the Cartesian coordinates


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.redundant_internal_transformation" class="docs-object-method">&nbsp;</a> 
```python
@property
redundant_internal_transformation(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L663)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L663?message=Update%20Docs)]
</div>
**LLM Docstring**

The redundant-to-non-redundant transformation matrix associated with the current internal coordinates, if the internal-coordinate system used a redundant coordinate generator.
  - `:returns`: `np.ndarray | None`
    > the redundant transformation, or `None` if not applicable


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.embedding_coords" class="docs-object-method">&nbsp;</a> 
```python
@property
embedding_coords(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L1012)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L1012?message=Update%20Docs)]
</div>
**LLM Docstring**

The indices of the internal-coordinate system's fixed embedding coordinates (translation/rotation degrees of freedom), if any are defined.
  - `:returns`: `np.ndarray | None`
    > the embedding-coordinate indices, or `None`


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.get_cartesians_by_internals" class="docs-object-method">&nbsp;</a> 
```python
get_cartesians_by_internals(self, order=None, strip_embedding=False, reembed=True, method=None, coords=None, **fd_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L1042)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L1042?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the Cartesians-by-internals Jacobian expansion up to the requested `order`, either via the fast route (inverting the internals-by-Cartesians Jacobian through a translation/rotation-invariant reduction, with caching) or the classic finite-difference/analytic route, depending on `method` (auto-selected based on the internal-coordinate system type) and whether `reembed`/`strip_embedding`/explicit `coords` are requested.
  - `order`: `int | None`
    > the highest derivative order to compute; if `None`, returns whatever is already cached
  - `strip_embedding`: `bool`
    > whether to strip the fixed embedding coordinates from the result
  - `reembed`: `bool`
    > whether to use the Eckart-reembedded (translation/rotation-invariant) formulation, for the `'fast'` method
  - `method`: `str | None`
    > which computation strategy to use (`'fast'` or `'classic'`); auto-selected if `None`
  - `coords`: `np.ndarray | None`
    > alternate Cartesian coordinates to compute the Jacobian at, instead of this embedding's own geometry
  - `fd_opts`: `dict`
    > extra finite-difference options forwarded to the underlying Jacobian computation
  - `:returns`: `list[np.ndarray]`
    > the Cartesians-by-internals Jacobian tensors, one per order


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.get_internals_by_cartesians" class="docs-object-method">&nbsp;</a> 
```python
get_internals_by_cartesians(self, order=None, strip_embedding=False, coords=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L1145)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L1145?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the internals-by-Cartesians Jacobian expansion up to the requested `order`, via finite difference/analytic derivatives (through `_get_cart_jacobs`), reshaping the results to `(..., ncart, ncart, ..., nint)`-style tensors and optionally stripping the fixed embedding coordinates.
  - `order`: `int | None`
    > the highest derivative order to compute; if `None`, returns whatever is already cached
  - `strip_embedding`: `bool`
    > whether to strip the fixed embedding coordinates from the result
  - `coords`: `np.ndarray | None`
    > alternate Cartesian coordinates to compute the Jacobian at
  - `opts`: `dict`
    > extra finite-difference options forwarded to `_get_cart_jacobs`
  - `:returns`: `list[np.ndarray]`
    > the internals-by-Cartesians Jacobian tensors, one per order


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.embed_coords" class="docs-object-method">&nbsp;</a> 
```python
embed_coords(self, coords, sel=None, in_paf=False, planar_ref_tolerance=None, proper_rotation=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L1195)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L1195?message=Update%20Docs)]
</div>
**LLM Docstring**

Eckart-embed a set of Cartesian coordinates onto this embedding's reference geometry.
  - `coords`: `np.ndarray`
    > the coordinates to embed
  - `sel`: `Iterable[int] | None`
    > subset of atoms to use for the embedding fit
  - `in_paf`: `bool`
    > whether to embed into the principal-axis frame
  - `planar_ref_tolerance`: `float | None`
    > tolerance for detecting a (near-)planar reference structure
  - `proper_rotation`: `bool`
    > whether to restrict the embedding to proper (determinant +1) rotations
  - `:returns`: `np.ndarray`
    > the Eckart-embedded coordinates


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.unembed_derivs" class="docs-object-method">&nbsp;</a> 
```python
unembed_derivs(self, coords, derivs, sel=None, in_paf=False, planar_ref_tolerance=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L1224)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L1224?message=Update%20Docs)]
</div>
**LLM Docstring**

Undo an Eckart embedding's rotation on a set of Cartesian derivative tensors, transforming them back by the combination of the embedding's axis frame and rotation.
  - `coords`: `np.ndarray`
    > the (embedded) coordinates the derivatives were computed at
  - `derivs`: `list[np.ndarray]`
    > the Cartesian derivative tensors to un-rotate
  - `sel`: `Iterable[int] | None`
    > subset of atoms used for the embedding fit
  - `in_paf`: `bool`
    > whether the embedding used the principal-axis frame
  - `planar_ref_tolerance`: `float | None`
    > tolerance for detecting a (near-)planar reference structure
  - `:returns`: `list[np.ndarray]`
    > the un-rotated derivative tensors


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.inertia_tensor" class="docs-object-method">&nbsp;</a> 
```python
@property
inertia_tensor(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L1260)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L1260?message=Update%20Docs)]
</div>
**LLM Docstring**

The molecule's inertia tensor at its current Cartesian coordinates.
  - `:returns`: `np.ndarray`
    > the inertia tensor


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.inertial_frame" class="docs-object-method">&nbsp;</a> 
```python
@property
inertial_frame(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L1275)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L1275?message=Update%20Docs)]
</div>
Provides the inertial axis frame
  - `:returns`: `_`
    >


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.inertial_frame_derivatives" class="docs-object-method">&nbsp;</a> 
```python
inertial_frame_derivatives(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L1297)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L1297?message=Update%20Docs)]
</div>
**LLM Docstring**

The first and second derivatives of the inertia tensor with respect to mass-weighted Cartesian displacements.
  - `:returns`: `list[np.ndarray]`
    > `[I0Y, I0YY]`, as returned by `StructuralProperties.get_prop_inertial_frame_derivatives`


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.translation_rotation_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
translation_rotation_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L1312)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L1312?message=Update%20Docs)]
</div>
**LLM Docstring**

The (cached) translation and rotation eigenvectors of the molecule at its current geometry.
  - `:returns`: `tuple`
    > the translation/rotation eigenvectors


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.get_translation_rotation_invariant_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_translation_rotation_invariant_transformation(self, order=0, mass_weighted=True, strip_embedding=True, coords=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L1331)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L1331?message=Update%20Docs)]
</div>
**LLM Docstring**

Build the transformation (and its inverse) that projects out the translational and rotational degrees of freedom from a set of Cartesian coordinates.
  - `order`: `int`
    > the derivative order of the transformation to build
  - `mass_weighted`: `bool`
    > whether the transformation should act on mass-weighted coordinates
  - `strip_embedding`: `bool`
    > whether to strip the fixed embedding coordinates from the result
  - `coords`: `np.ndarray | None`
    > alternate Cartesian coordinates to build the transformation at, instead of this embedding's own geometry
  - `:returns`: `tuple`
    > the translation/rotation-invariant transformation and its inverse
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Molecools/CoordinateSystems/MolecularEmbedding.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Molecools/CoordinateSystems/MolecularEmbedding.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Molecools/CoordinateSystems/MolecularEmbedding.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Molecools/CoordinateSystems/MolecularEmbedding.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems.py#L32?message=Update%20Docs)   
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