## <a id="Psience.Modes.MixtureModes.MixtureModes">MixtureModes</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes.py#L18)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes.py#L18?message=Update%20Docs)]
</div>

A `McUtils.Coordinerds.CoordinateSystem` object that expresses coordinates as
a rotation on some base set of coordinates with some associated frequencies.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
name: str
zero_freq_cutoff: float
ModeData: ModeData
default_zero_freq_cutoff: float
localization_type: str
localization_zero_freq_cutoff: float
LocalizationMethods: LocalizationMethods
localization_options: tuple
```
<a id="Psience.Modes.MixtureModes.MixtureModes.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, basis, coeffs, freqs=None, origin=None, masses=None, inverse=None, mass_weighted=False, frequency_scaled=False, g_matrix=None, name=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes.py#L24)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes.py#L24?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a `MixtureModes` object -- a `CoordinateSystem` expressing coordinates as a (possibly non-orthogonal) linear rotation of some base coordinate system, with associated frequencies, masses, and mass-weighting/frequency-scaling/G-matrix bookkeeping.
  - `basis`: `CoordinateSystem`
    > the base coordinate system the modes are expressed in
  - `coeffs`: `np.ndarray | list[np.ndarray]`
    > the mode coefficient (coordinates-by-modes) matrix, or a list of such matrices/tensors representing higher-order terms of a nonlinear coordinate expansion (the first entry is used as the linear matrix)
  - `freqs`: `np.ndarray | None`
    > the vibrational frequencies associated with each mode
  - `origin`: `np.ndarray | None`
    > the reference geometry the modes are expanded about
  - `masses`: `np.ndarray | None`
    > the atomic masses
  - `inverse`: `np.ndarray | None`
    > the modes-by-coordinates (inverse) matrix
  - `mass_weighted`: `bool`
    > whether the modes are mass-weighted
  - `frequency_scaled`: `bool`
    > whether the modes are frequency-scaled (dimensionless)
  - `g_matrix`: `np.ndarray | None`
    > the G-matrix associated with the modes
  - `name`: `str | None`
    > a name for the mode set; defaults to the class's `name` attribute
  - `:returns`: `None`
    > None


<a id="Psience.Modes.MixtureModes.MixtureModes.to_state" class="docs-object-method">&nbsp;</a> 
```python
to_state(self, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L98)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L98?message=Update%20Docs)]
</div>
**LLM Docstring**

Serialize this mode set's essential data (basis, matrix, inverse, frequencies, masses, mass-weighting/frequency-scaling flags, G-matrix) into a plain dict.
  - `serializer`: `object`
    > the serializer used to recursively serialize the `basis` object
  - `:returns`: `dict`
    > the serialized state dict


<a id="Psience.Modes.MixtureModes.MixtureModes.from_state" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_state(cls, data, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L119)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L119?message=Update%20Docs)]
</div>
**LLM Docstring**

Reconstruct a `MixtureModes` from a previously serialized state dict.
  - `data`: `dict`
    > the serialized state, as produced by `to_state`
  - `serializer`: `object`
    > the serializer used to recursively deserialize the `basis` object
  - `:returns`: `MixtureModes`
    > the reconstructed mode set


<a id="Psience.Modes.MixtureModes.MixtureModes.prep_modes" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
prep_modes(cls, modes): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L144)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L144?message=Update%20Docs)]
</div>
**LLM Docstring**

Coerce a variety of mode representations (an already-built `MixtureModes`-like object, a dict of matrix/inverse/basis/freqs, an object exposing `.matrix`, or a raw coefficient array) into a proper `MixtureModes` instance, inferring a sensible default basis (Cartesian or generic internal) if none is given.
  - `modes`: `object | dict | np.ndarray`
    > the mode data to coerce
  - `:returns`: `MixtureModes`
    > the constructed (or passed-through) mode set


<a id="Psience.Modes.MixtureModes.MixtureModes.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L207)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L207?message=Update%20Docs)]
</div>
Takes a slice of the modes
  - `item`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Modes.MixtureModes.MixtureModes.modify" class="docs-object-method">&nbsp;</a> 
```python
modify(self, matrix=None, *, freqs=None, origin=None, masses=None, inverse=None, name=None, mass_weighted=None, frequency_scaled=None, g_matrix=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L235)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L235?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a new `MixtureModes` (of the same concrete subclass) with the given fields overridden, defaulting each unspecified field to this object's own current value.
  - `matrix`: `np.ndarray | None`
    > replacement modes-by-coordinates matrix; defaults to `self.matrix`
  - `freqs`: `np.ndarray | None`
    > replacement frequencies; defaults to `self.freqs`
  - `origin`: `np.ndarray | None`
    > replacement reference geometry; defaults to `self.origin`
  - `masses`: `np.ndarray | None`
    > replacement masses; defaults to `self.masses`
  - `inverse`: `np.ndarray | None`
    > replacement coords-by-modes matrix; defaults to `self.inverse`
  - `name`: `str | None`
    > replacement name; defaults to `self.name`
  - `mass_weighted`: `bool | None`
    > replacement mass-weighting flag; defaults to `self.mass_weighted`
  - `frequency_scaled`: `bool | None`
    > replacement frequency-scaling flag; defaults to `self.frequency_scaled`
  - `g_matrix`: `np.ndarray | None`
    > replacement G-matrix; defaults to `self.g_matrix`
  - `:returns`: `MixtureModes`
    > the new, modified mode set


<a id="Psience.Modes.MixtureModes.MixtureModes.rotate" class="docs-object-method">&nbsp;</a> 
```python
rotate(self, rot, in_place=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L286)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L286?message=Update%20Docs)]
</div>
**LLM Docstring**

Intended to rotate the mode set by a given rotation. Not implemented -- the intended semantics were judged too ambiguous.
  - `rot`: `object`
    > the rotation to apply
  - `in_place`: `bool`
    > whether to apply the rotation in place
  - `:returns`: `None`
    > never returns


<a id="Psience.Modes.MixtureModes.MixtureModes.transform" class="docs-object-method">&nbsp;</a> 
```python
transform(self, tf, inv=None, origin=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L302)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L302?message=Update%20Docs)]
</div>
**LLM Docstring**

Intended to apply a linear transformation (and its inverse) to the mode set. Not implemented -- immediately raises before reaching the (dead) implementation code below it, which is retained but unreachable.
  - `tf`: `np.ndarray`
    > the transformation matrix
  - `inv`: `np.ndarray | None`
    > the inverse transformation; computed from `tf` if `tf` is square and `inv` isn't given (in the unreachable code)
  - `origin`: `np.ndarray | None`
    > the reference geometry to transform; defaults to `self.origin` (in the unreachable code)
  - `:returns`: `None`
    > never returns


<a id="Psience.Modes.MixtureModes.MixtureModes.cartesian_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
cartesian_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L347)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L347?message=Update%20Docs)]
</div>
**LLM Docstring**

Whether the modes are expressed relative to a Cartesian (as opposed to internal-coordinate) origin, inferred from the origin's dimensionality.
  - `:returns`: `bool`
    > `True` if the origin is 2-dimensional (`(natoms, 3)`)


<a id="Psience.Modes.MixtureModes.MixtureModes.embed_coords" class="docs-object-method">&nbsp;</a> 
```python
embed_coords(self, carts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L359)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L359?message=Update%20Docs)]
</div>
**LLM Docstring**

Convert a batch of Cartesian displacement structures into mode coordinates, by subtracting the reference origin and projecting through the coords-by-modes (inverse) matrix.
  - `carts`: `np.ndarray`
    > the Cartesian structures to embed
  - `:returns`: `np.ndarray`
    > the mode coordinates


<a id="Psience.Modes.MixtureModes.MixtureModes.unembed_coords" class="docs-object-method">&nbsp;</a> 
```python
unembed_coords(self, mode_coords): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L375)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L375?message=Update%20Docs)]
</div>
**LLM Docstring**

Convert a batch of mode coordinates back into Cartesian structures, by projecting through the modes-by-coordinates matrix and adding back the reference origin.
  - `mode_coords`: `np.ndarray`
    > the mode coordinates to unembed
  - `:returns`: `np.ndarray`
    > the Cartesian structures


<a id="Psience.Modes.MixtureModes.MixtureModes.total_transformation" class="docs-object-method">&nbsp;</a> 
```python
@property
total_transformation(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L393)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L393?message=Update%20Docs)]
</div>
**LLM Docstring**

The full (possibly multi-term, nonlinear) forward coordinate-expansion tensors this mode set was constructed with.
  - `:returns`: `list[np.ndarray]`
    > the list of expansion-order tensors, starting with the linear modes-by-coordinates matrix


<a id="Psience.Modes.MixtureModes.MixtureModes.inverse_transformation" class="docs-object-method">&nbsp;</a> 
```python
@property
inverse_transformation(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L404)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L404?message=Update%20Docs)]
</div>
**LLM Docstring**

The (lazily computed and cached) inverse of `total_transformation`, i.e. the Taylor-series expansion mapping mode coordinates back to base coordinates.
  - `:returns`: `list[np.ndarray]`
    > the inverse expansion-order tensors


<a id="Psience.Modes.MixtureModes.MixtureModes.embed_derivs" class="docs-object-method">&nbsp;</a> 
```python
embed_derivs(self, derivs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L419)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L419?message=Update%20Docs)]
</div>
**LLM Docstring**

Re-express a set of derivative tensors (with respect to the base coordinates) in terms of mode coordinates, by re-expanding through `total_transformation`.
  - `derivs`: `list[np.ndarray]`
    > the derivative tensors to re-express
  - `:returns`: `list[np.ndarray]`
    > the re-expressed derivative tensors


<a id="Psience.Modes.MixtureModes.MixtureModes.unembed_derivs" class="docs-object-method">&nbsp;</a> 
```python
unembed_derivs(self, derivs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L431)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L431?message=Update%20Docs)]
</div>
**LLM Docstring**

Re-express a set of derivative tensors (with respect to mode coordinates) back in terms of the base coordinates, by re-expanding through `inverse_transformation`.
  - `derivs`: `list[np.ndarray]`
    > the derivative tensors to re-express
  - `:returns`: `list[np.ndarray]`
    > the re-expressed derivative tensors


<a id="Psience.Modes.MixtureModes.MixtureModes.is_cartesian" class="docs-object-method">&nbsp;</a> 
```python
@property
is_cartesian(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L444)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L444?message=Update%20Docs)]
</div>
**LLM Docstring**

Whether the modes are expressed over a Cartesian coordinate space, either inferred from the mode-matrix row count matching `3 * len(masses)` (if masses are known) or from the base coordinate system's name.
  - `:returns`: `bool`
    > whether the modes are Cartesian-basis


<a id="Psience.Modes.MixtureModes.MixtureModes.coords_by_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
coords_by_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L458)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L458?message=Update%20Docs)]
</div>
**LLM Docstring**

The coordinates-by-modes (inverse) transformation matrix.
  - `:returns`: `np.ndarray`
    > `self.inverse`


<a id="Psience.Modes.MixtureModes.MixtureModes.modes_by_coords" class="docs-object-method">&nbsp;</a> 
```python
@property
modes_by_coords(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L469)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L469?message=Update%20Docs)]
</div>
**LLM Docstring**

The modes-by-coordinates transformation matrix.
  - `:returns`: `np.ndarray`
    > `self.matrix`


<a id="Psience.Modes.MixtureModes.MixtureModes.compute_local_transformations" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
compute_local_transformations(cls, f, g): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L548)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L548?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the pair of diagonal rescaling transformations (for the force-constant matrix and its inverse-G counterpart) that define the "local mode" representation from diagonal `f`/`g` values.
  - `f`: `np.ndarray`
    > the diagonal force-constant values
  - `g`: `np.ndarray`
    > the diagonal G-matrix values
  - `:returns`: `list[np.ndarray]`
    > `[local_tf(f, g), local_tf(g, f)]`


<a id="Psience.Modes.MixtureModes.MixtureModes.compute_local_hessian" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
compute_local_hessian(cls, f, g): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L567)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L567?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the rescaled ("local mode") force-constant matrix from full `f`/`g` matrices, using the diagonal rescaling from `_local_tf`.
  - `f`: `np.ndarray`
    > the force-constant (Hessian) matrix
  - `g`: `np.ndarray`
    > the G-matrix
  - `:returns`: `np.ndarray`
    > the rescaled local Hessian


<a id="Psience.Modes.MixtureModes.MixtureModes.compute_local_gmatrix" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
compute_local_gmatrix(cls, f, g): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L584)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L584?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the rescaled ("local mode") G-matrix from full `f`/`g` matrices, using the diagonal rescaling from `_local_tf`.
  - `f`: `np.ndarray`
    > the force-constant (Hessian) matrix
  - `g`: `np.ndarray`
    > the G-matrix
  - `:returns`: `np.ndarray`
    > the rescaled local G-matrix


<a id="Psience.Modes.MixtureModes.MixtureModes.compute_hessian" class="docs-object-method">&nbsp;</a> 
```python
compute_hessian(self, system='modes'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L601)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L601?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the (signed-frequency-squared) Hessian matrix in either the mode basis or the underlying coordinate basis, by reconstructing it from the stored frequencies and the mode transformation.
  - `system`: `str`
    > which basis to compute the Hessian in, `'modes'` or `'coords'`
  - `:returns`: `np.ndarray`
    > the Hessian matrix


<a id="Psience.Modes.MixtureModes.MixtureModes.compute_gmatrix" class="docs-object-method">&nbsp;</a> 
```python
compute_gmatrix(self, system='modes', return_fractional=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L622)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L622?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the G-matrix in either the mode basis or the underlying coordinate basis, using whichever information is available (mass-weighting, a stored G-matrix, or Cartesian masses).
  - `system`: `str`
    > which basis to compute the G-matrix in, `'modes'` or `'coords'`
  - `return_fractional`: `bool`
    > whether to also return the G-matrix's square root and inverse square root
  - `:returns`: `np.ndarray | tuple | None`
    > the G-matrix, or `(g, g12, gi12)` if `return_fractional` is set; `None` if no G-matrix information is available in the `'modes'` case


<a id="Psience.Modes.MixtureModes.MixtureModes.compute_freqs" class="docs-object-method">&nbsp;</a> 
```python
compute_freqs(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L691)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L691?message=Update%20Docs)]
</div>
**LLM Docstring**

Recompute the vibrational frequencies for this mode set from its coordinate-basis Hessian and G-matrix, via a generalized eigenvalue solve.
  - `:returns`: `np.ndarray`
    > the signed square-root frequencies


<a id="Psience.Modes.MixtureModes.MixtureModes.local_hessian" class="docs-object-method">&nbsp;</a> 
```python
@property
local_hessian(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L708)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L708?message=Update%20Docs)]
</div>
**LLM Docstring**

The "local mode" (diagonally rescaled) force-constant matrix for this mode set.
  - `:returns`: `np.ndarray`
    > the local Hessian


<a id="Psience.Modes.MixtureModes.MixtureModes.local_gmatrix" class="docs-object-method">&nbsp;</a> 
```python
@property
local_gmatrix(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L723)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L723?message=Update%20Docs)]
</div>
**LLM Docstring**

The "local mode" (diagonally rescaled) G-matrix for this mode set.
  - `:returns`: `np.ndarray`
    > the local G-matrix


<a id="Psience.Modes.MixtureModes.MixtureModes.local_freqs" class="docs-object-method">&nbsp;</a> 
```python
@property
local_freqs(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L738)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L738?message=Update%20Docs)]
</div>
**LLM Docstring**

The diagonal entries of the local-mode Hessian, giving an approximate per-mode "local" force constant/frequency.
  - `:returns`: `np.ndarray`
    > the local-mode diagonal values


<a id="Psience.Modes.MixtureModes.MixtureModes.local_mode_transformations" class="docs-object-method">&nbsp;</a> 
```python
@property
local_mode_transformations(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L750)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L750?message=Update%20Docs)]
</div>
**LLM Docstring**

The pair of diagonal rescaling transformations mapping between this mode set's coordinate-basis Hessian/G-matrix and their local-mode counterparts.
  - `:returns`: `list[np.ndarray]`
    > `[hessian_scaling, gmatrix_scaling]`


<a id="Psience.Modes.MixtureModes.MixtureModes.get_nearest_mode_transform" class="docs-object-method">&nbsp;</a> 
```python
get_nearest_mode_transform(self, alternate_modes: numpy.ndarray, mass_weighted=None, atoms=None, maximum_similarity=True, unitarize=True, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L765)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L765?message=Update%20Docs)]
</div>
**LLM Docstring**

Find the transformation (and its inverse) that best aligns this mode set with an externally supplied set of `alternate_modes`, either via a maximum-similarity (orthogonal Procrustes-style) alignment or via direct projection (optionally unitarized), with optional mass-(un)weighting and atom-restriction beforehand.
  - `alternate_modes`: `np.ndarray`
    > the external mode matrix to align/project against
  - `mass_weighted`: `bool | None`
    > if given, coerces this mode set to (or from) mass-weighted form before comparing
  - `atoms`: `Iterable[int] | None`
    > restrict the comparison to a subset of atoms via an atom projector (Cartesian modes only)
  - `maximum_similarity`: `bool`
    > whether to use a maximum-similarity (orthogonal Procrustes) alignment rather than direct projection
  - `unitarize`: `bool`
    > whether to unitarize the resulting transformation (direct-projection path only)
  - `masses`: `np.ndarray | None`
    > masses to use for the mass-(un)weighting conversion
  - `:returns`: `tuple[np.ndarray, np.ndarray]`
    > `(tf, inv)` -- the forward and inverse transformation matrices


<a id="Psience.Modes.MixtureModes.MixtureModes.get_normal_modes" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_normal_modes(cls, f_matrix, mass_spec, remove_transrot=True, dimensionless=False, mass_weighted=None, zero_freq_cutoff=None, return_gmatrix=False, projector=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L834)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L834?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute ordinary harmonic normal modes from a force-constant matrix and mass specification (masses, element symbols, or a full G-matrix), via a generalized eigenvalue solve on the mass-weighted Hessian, optionally removing near-zero-frequency (translation/rotation) modes, applying a custom constraint projector, mass-weighting and/or frequency-scaling (dimensionless-izing) the result, and supporting batched inputs (including "ragged" batches where a different number of modes survives per entry).
  - `f_matrix`: `np.ndarray`
    > the force-constant (Hessian) matrix/matrices
  - `mass_spec`: `np.ndarray | list[str]`
    > the atomic masses (or element symbols), or an already-built G-matrix/mass matrix
  - `remove_transrot`: `bool`
    > whether to discard near-zero-frequency modes
  - `dimensionless`: `bool`
    > whether the resulting modes should be frequency-scaled (dimensionless)
  - `mass_weighted`: `bool | None`
    > whether the resulting modes should be mass-weighted; defaults to `dimensionless` when `mass_spec` is a bare mass vector
  - `zero_freq_cutoff`: `float | None`
    > the frequency cutoff below which modes are discarded; defaults to `cls.default_zero_freq_cutoff`
  - `return_gmatrix`: `bool`
    > whether to also return the (broadcast) mass/G-matrix used
  - `projector`: `np.ndarray | None`
    > an additional constraint projector to apply to the mass-weighted Hessian before diagonalizing
  - `:returns`: `object | tuple`
    > the mode data (`ModeData(freqs, modes, inverse)`), or `(mode_data, mass_spec)` if `return_gmatrix` is set


<a id="Psience.Modes.MixtureModes.MixtureModes.get_projected_localized_mode_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_projected_localized_mode_transformation(self, projectors, masses=None, origin=None, localization_type=None, allow_mode_mixing=False, maximum_similarity=False, unitarize=True, zero_freq_cutoff=None, orthogonal_projection=False, project_zero_gmatrix_modes=True, project_zero_gmatrix_cutoff=1e-08, atoms=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1042)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1042?message=Update%20Docs)]
</div>
**LLM Docstring**

Shared driver behind the various coordinate/atom/fragment-localized-mode constructors: diagonalizes the (optionally translation/rotation-projected) mass-weighted Hessian restricted to (or blocked by) a set of projectors, either mixing all projected subspaces together (`allow_mode_mixing`) or diagonalizing each independently and concatenating the results, then finds the transformation aligning the current modes onto the resulting localized modes via `get_nearest_mode_transform`.
  - `projectors`: `list[np.ndarray]`
    > the projection matrix (or matrices) defining the coordinate subspace(s) to localize within
  - `masses`: `np.ndarray | None`
    > masses to use instead of `self.masses`
  - `origin`: `np.ndarray | None`
    > reference geometry to use for the translation/rotation-invariant transformation (`'direct'` localization type only)
  - `localization_type`: `str | None`
    > `'ned'`/other for the standard projected-Hessian approach, or `'direct'` to first factor out translation/rotation via an explicit invariant transformation (Cartesian modes only)
  - `allow_mode_mixing`: `bool`
    > whether to mix all projected subspaces together into a single diagonalization rather than treating each independently
  - `maximum_similarity`: `bool`
    > forwarded to `get_nearest_mode_transform`
  - `unitarize`: `bool`
    > forwarded to `get_nearest_mode_transform`
  - `zero_freq_cutoff`: `float | None`
    > the frequency cutoff used when diagonalizing each subspace; defaults to `self.localization_zero_freq_cutoff`
  - `orthogonal_projection`: `bool`
    > whether non-square entries in `projectors` should be treated as bases for an orthogonal (rather than oblique) projection
  - `project_zero_gmatrix_modes`: `bool`
    > accepted for interface consistency but not used in this method's body
  - `project_zero_gmatrix_cutoff`: `float`
    > accepted for interface consistency but not used in this method's body
  - `atoms`: `object | None`
    > accepted for interface consistency with sibling localizers but not used in this method's body
  - `:returns`: `tuple[np.ndarray, np.ndarray]`
    > `(tf, inv)`, the localizing transformation and its inverse


<a id="Psience.Modes.MixtureModes.MixtureModes.get_atom_localized_mode_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_atom_localized_mode_transformation(self, atoms, masses=None, origin=None, localization_type='ned', allow_mode_mixing=False, maximum_similarity=False, orthogonal_projection=False, unitarize=True, zero_freq_cutoff=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1180)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1180?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a localized-mode transformation that concentrates vibrational character onto (or, if `orthogonal_projection`, away from) a specified set of atoms, via `get_projected_localized_mode_transformation` with per-atom (or combined) atom-selection projectors.
  - `atoms`: `int | Iterable[int]`
    > the atom index (or indices) to localize modes onto
  - `masses`: `np.ndarray | None`
    > masses to use instead of `self.masses`
  - `origin`: `np.ndarray | None`
    > reference geometry, forwarded to `get_projected_localized_mode_transformation`
  - `localization_type`: `str`
    > the localization scheme to use
  - `allow_mode_mixing`: `bool`
    > whether to build one combined projector spanning all `atoms` rather than one per atom
  - `maximum_similarity`: `bool`
    > forwarded to `get_projected_localized_mode_transformation`
  - `orthogonal_projection`: `bool`
    > whether to localize onto the complement of `atoms` (all other atoms) instead of `atoms` itself
  - `unitarize`: `bool`
    > forwarded to `get_projected_localized_mode_transformation`
  - `zero_freq_cutoff`: `float | None`
    > forwarded to `get_projected_localized_mode_transformation`
  - `:returns`: `tuple[np.ndarray, np.ndarray]`
    > `(tf, inv)`, the localizing transformation and its inverse


<a id="Psience.Modes.MixtureModes.MixtureModes.get_fragment_localized_mode_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_fragment_localized_mode_transformation(self, fragment, masses=None, origin=None, localization_type='ned', allow_mode_mixing=True, maximum_similarity=False, orthogonal_projection=False, unitarize=True, zero_freq_cutoff=None, **etc): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1245)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1245?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a localized-mode transformation for a molecular fragment (a set of atoms), via `get_atom_localized_mode_transformation`.
  - `fragment`: `Iterable[int]`
    > the atom indices making up the fragment
  - `masses`: `np.ndarray | None`
    > masses to use instead of `self.masses`
  - `origin`: `np.ndarray | None`
    > reference geometry
  - `localization_type`: `str`
    > the localization scheme to use
  - `allow_mode_mixing`: `bool`
    > whether to build one combined projector for the whole fragment
  - `maximum_similarity`: `bool`
    > forwarded through to the underlying localizer
  - `orthogonal_projection`: `bool`
    > whether to localize onto the complement of the fragment instead
  - `unitarize`: `bool`
    > forwarded through to the underlying localizer
  - `zero_freq_cutoff`: `float | None`
    > forwarded through to the underlying localizer
  - `etc`: `dict`
    > extra options forwarded to `get_atom_localized_mode_transformation`
  - `:returns`: `tuple[np.ndarray, np.ndarray]`
    > `(tf, inv)`, the localizing transformation and its inverse


<a id="Psience.Modes.MixtureModes.MixtureModes.get_coordinate_projected_localized_mode_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_coordinate_projected_localized_mode_transformation(self, coordinate_constraints, atoms=None, masses=None, origin=None, localization_type='ned', allow_mode_mixing=False, maximum_similarity=False, orthogonal_projection=True, unitarize=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1296)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1296?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a localized-mode transformation that concentrates vibrational character along a set of internal-coordinate directions (bond/angle/dihedral bases for Cartesian modes, or by coordinate index for internal-coordinate modes), optionally restricted to a subset of atoms, via `get_projected_localized_mode_transformation`.
  - `coordinate_constraints`: `object | Iterable[int]`
    > the internal-coordinate specification(s) (for Cartesian modes) or coordinate index/indices (for internal-coordinate modes) to localize along
  - `atoms`: `Iterable[int] | None`
    > restrict the resulting projector(s) to this subset of atoms (Cartesian modes only)
  - `masses`: `np.ndarray | None`
    > masses to use instead of `self.masses`
  - `origin`: `np.ndarray | None`
    > reference geometry, forwarded to `get_projected_localized_mode_transformation`
  - `localization_type`: `str`
    > the localization scheme to use
  - `allow_mode_mixing`: `bool`
    > whether to combine all coordinate projectors into one before localizing
  - `maximum_similarity`: `bool`
    > forwarded to `get_projected_localized_mode_transformation`
  - `orthogonal_projection`: `bool`
    > whether the coordinate bases define orthogonal (rather than oblique) projections
  - `unitarize`: `bool`
    > forwarded to `get_projected_localized_mode_transformation`
  - `:returns`: `tuple[np.ndarray, np.ndarray]`
    > `(tf, inv)`, the localizing transformation and its inverse


<a id="Psience.Modes.MixtureModes.MixtureModes.get_internal_localized_mode_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_internal_localized_mode_transformation(self, expansion_coordinates: 'Iterable[Iterable[int]|dict]', fixed_atoms=None, mass_weighted=False, project_transrot=True, atoms=None, maximum_similarity=False, orthogonal_projection=False, projection=False, allow_mode_mixing=False, unitarize=True, origin=None, masses=None, localization_type='ned'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1414)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1414?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a localized-mode transformation aligned with a set of internal-coordinate displacement directions (computed at the reference geometry), either by projecting the Hessian onto (or away from) those directions and re-diagonalizing (`projection`/`orthogonal_projection`), or by directly finding the nearest-mode alignment to the raw coordinate-derivative directions (optionally mass-weighted and/or translation/rotation-projected).
  - `expansion_coordinates`: `Iterable[Iterable[int]] | dict`
    > the internal coordinate(s) whose Cartesian derivative directions define the localization target
  - `fixed_atoms`: `Iterable[int] | None`
    > atoms to hold fixed when computing the internal-coordinate derivative tensors
  - `mass_weighted`: `bool`
    > whether to mass-weight the coordinate-derivative directions
  - `project_transrot`: `bool`
    > whether to project out translational/rotational components from the derivative directions (direct-alignment path only)
  - `atoms`: `Iterable[int] | None`
    > restrict the projector(s)/alignment to a subset of atoms
  - `maximum_similarity`: `bool`
    > forwarded to the underlying localizer
  - `orthogonal_projection`: `bool`
    > whether to use orthogonal (rather than oblique) projections in the projection-based path
  - `projection`: `bool`
    > whether to use the projection-and-rediagonalize approach rather than direct nearest-mode alignment
  - `allow_mode_mixing`: `bool`
    > whether to combine all coordinate directions into a single projected subspace (projection-based path)
  - `unitarize`: `bool`
    > forwarded to the underlying localizer
  - `origin`: `np.ndarray | None`
    > reference geometry to use
  - `masses`: `np.ndarray | None`
    > masses to use instead of `self.masses`
  - `localization_type`: `str`
    > the localization scheme to use for the projection-based path
  - `:returns`: `tuple[np.ndarray, np.ndarray]`
    > `(tf, inv)`, the localizing transformation and its inverse


<a id="Psience.Modes.MixtureModes.MixtureModes.get_displacement_localized_mode_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_displacement_localized_mode_transformation(self, mode_blocks=None, atoms=None, mass_weighted=True, unitarize=True, **maximizer_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1567)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1567?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a unitary localized-mode transformation by iteratively (Jacobi-style) rotating blocks of modes to maximize a displacement-localization criterion (via `nput.jacobi_maximize`/`nput.displacement_localizing_rotation_generator`), optionally restricted to a subset of atoms and/or applied independently within specified mode blocks.
  - `mode_blocks`: `Iterable[int] | list[Iterable[int]] | None`
    > groups of mode indices to localize independently; a single flat list of indices localizes all modes together, and `None` localizes all modes as one block
  - `atoms`: `Iterable[int] | None`
    > restrict the localization criterion to a subset of atoms
  - `mass_weighted`: `bool`
    > whether to localize the mass-weighted (rather than un-mass-weighted) mode matrix
  - `unitarize`: `bool`
    > must be `True`; the Jacobi-rotation scheme only produces unitary transformations
  - `maximizer_opts`: `dict`
    > extra options forwarded to `nput.jacobi_maximize`
  - `:returns`: `tuple[np.ndarray, np.ndarray]`
    > `(tf_inv, tf_inv.T)`, the (unitary) localizing transformation's inverse and its transpose (used as the forward transformation)


<a id="Psience.Modes.MixtureModes.MixtureModes.get_mass_scaled_mode_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_mass_scaled_mode_transformation(self, mass_scaling, *, atoms, localization_cutoff=0.8, num_modes=None, project_transrot=False, unitarize=True, **diag_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1619)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1619?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a localized-mode transformation by artificially scaling the masses of a set of atoms (making them effectively heavier or lighter) and re-diagonalizing, then selecting the resulting modes that are most concentrated on those atoms (by projected-displacement norm), optionally applying a localization cutoff or a fixed mode count.
  - `mass_scaling`: `float | np.ndarray`
    > the scale factor (or per-atom scale factors) to apply to the selected atoms' masses
  - `atoms`: `int | Iterable[int]`
    > the atom(s) whose masses should be scaled
  - `localization_cutoff`: `float | None`
    > minimum normalized displacement-on-atoms required for a mode to be selected; if `None`, the top `num_modes` (or `3*len(atoms)`) modes by that metric are taken instead
  - `num_modes`: `int | None`
    > the number of modes to select; defaults to `3 * len(atoms)`
  - `project_transrot`: `bool`
    > whether to project out translation/rotation using the mass-scaled masses before diagonalizing
  - `unitarize`: `bool`
    > accepted for interface consistency with sibling localizers but not used in this method's body
  - `diag_opts`: `dict`
    > extra options forwarded to `NormalModes.get_normal_modes`
  - `:returns`: `tuple[np.ndarray, np.ndarray] | tuple[None, None]`
    > `(tf, tf.T)`, the selected mode transformation and its transpose, or `(None, None)` if no modes clear `localization_cutoff`


<a id="Psience.Modes.MixtureModes.MixtureModes.localizer_dispatch" class="docs-object-method">&nbsp;</a> 
```python
@property
localizer_dispatch(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1718)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1718?message=Update%20Docs)]
</div>
**LLM Docstring**

The mapping from `LocalizationMethods` value to the `(constructor_method, primary_argument_name)` pair used by `localize` to dispatch a localization request to the right underlying method.
  - `:returns`: `dict`
    > the method-name-to-`(constructor, arg_name)` mapping


<a id="Psience.Modes.MixtureModes.MixtureModes.localize" class="docs-object-method">&nbsp;</a> 
```python
localize(self, method=None, *, atoms=None, fragment=None, target_modes=None, internals=None, mode_blocks=None, coordinate_constraints=None, projections=None, reorthogonalize=None, mass_scaling=None, unitarize=True, allow_mode_mixing=False, project_zero_gmatrix_modes=None, project_zero_gmatrix_cutoff=1e-08, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1758)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1758?message=Update%20Docs)]
</div>
**LLM Docstring**

Top-level entry point for building a `LocalizedModes` object: infers which localization method to use from whichever keyword argument was supplied (if `method` isn't given explicitly), dispatches to the corresponding `get_*_localized_mode_transformation` method, optionally projects out (and, if mixing is allowed, re-diagonalizes within) any modes with a near-zero effective G-matrix eigenvalue, optionally re-orthogonalizes the result, and wraps the resulting transformation in a `LocalizedModes` object.
  - `method`: `str | LocalizationMethods | None`
    > which localization method to use; inferred from the other arguments if not given
  - `atoms`: `int | Iterable[int] | None`
    > atom(s) to localize onto, for the `'atoms'` method (or as a restriction for others)
  - `fragment`: `Iterable[int] | None`
    > fragment atoms to localize onto, for the `'fragment'` method
  - `target_modes`: `np.ndarray | None`
    > external modes to align to, for the `'target_modes'` method
  - `internals`: `object | None`
    > internal coordinates to localize along, for the `'coordinates'` method
  - `mode_blocks`: `object | None`
    > mode index groupings, for the `'displacements'` method
  - `coordinate_constraints`: `object | None`
    > coordinate constraints, for the `'constraints'` method
  - `projections`: `list[np.ndarray] | None`
    > explicit projectors, for the `'projections'` method
  - `reorthogonalize`: `bool | None`
    > whether to re-orthogonalize the localized modes' mass-weighted representation after localization; defaults to `not unitarize`
  - `mass_scaling`: `float | np.ndarray | None`
    > mass-scaling factor(s), for the `'mass_scaling'` method
  - `unitarize`: `bool`
    > whether the underlying localizer should produce a unitary transformation
  - `allow_mode_mixing`: `bool`
    > whether the underlying localizer is allowed to mix modes across different projected subspaces
  - `project_zero_gmatrix_modes`: `bool | None`
    > whether to project out (and handle) modes with a near-zero effective G-matrix eigenvalue after localizing; defaults to `allow_mode_mixing`
  - `project_zero_gmatrix_cutoff`: `float`
    > the G-matrix eigenvalue magnitude below which a mode is treated as zero
  - `opts`: `dict`
    > extra options forwarded to the dispatched localizer method
  - `:returns`: `LocalizedModes | None`
    > the constructed `LocalizedModes` object, or `None` if the underlying localizer returned no transformation


<a id="Psience.Modes.MixtureModes.MixtureModes.make_mass_weighted" class="docs-object-method">&nbsp;</a> 
```python
make_mass_weighted(self, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1908)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1908?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a mass-weighted version of this mode set by rescaling the mode/inverse matrices and origin through the G-matrix's square root/inverse-square-root; returns `self` unchanged if already mass-weighted.
  - `masses`: `np.ndarray | None`
    > masses to use instead of `self.masses`
  - `:returns`: `MixtureModes`
    > the mass-weighted mode set


<a id="Psience.Modes.MixtureModes.MixtureModes.remove_mass_weighting" class="docs-object-method">&nbsp;</a> 
```python
remove_mass_weighting(self, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1931)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1931?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a non-mass-weighted version of this mode set by undoing the mass-weighting rescaling; returns `self` unchanged if not already mass-weighted.
  - `masses`: `np.ndarray | None`
    > masses to use instead of `self.masses`
  - `:returns`: `MixtureModes`
    > the non-mass-weighted mode set


<a id="Psience.Modes.MixtureModes.MixtureModes.make_frequency_scaled" class="docs-object-method">&nbsp;</a> 
```python
make_frequency_scaled(self, freqs=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1972)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1972?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a frequency-scaled (dimensionless) version of this mode set by rescaling the mode/inverse matrices by the per-mode frequency-scaling factor; returns `self` unchanged if already frequency-scaled.
  - `freqs`: `np.ndarray | None`
    > the frequencies to compute the scaling from; defaults to `self.freqs`
  - `:returns`: `MixtureModes`
    > the frequency-scaled mode set


<a id="Psience.Modes.MixtureModes.MixtureModes.remove_frequency_scaling" class="docs-object-method">&nbsp;</a> 
```python
remove_frequency_scaling(self, freqs=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1994)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1994?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a non-frequency-scaled version of this mode set by undoing the frequency-scaling rescaling; returns `self` unchanged if not already frequency-scaled.
  - `freqs`: `np.ndarray | None`
    > the frequencies to compute the scaling from; defaults to `self.freqs`
  - `:returns`: `MixtureModes`
    > the non-frequency-scaled mode set


<a id="Psience.Modes.MixtureModes.MixtureModes.make_dimensionless" class="docs-object-method">&nbsp;</a> 
```python
make_dimensionless(self, freqs=None, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L2016)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L2016?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a fully dimensionless version of this mode set: mass-weighted and frequency-scaled.
  - `freqs`: `np.ndarray | None`
    > the frequencies to compute the frequency scaling from
  - `masses`: `np.ndarray | None`
    > the masses to compute the mass-weighting from
  - `:returns`: `MixtureModes`
    > the dimensionless mode set


<a id="Psience.Modes.MixtureModes.MixtureModes.make_dimensioned" class="docs-object-method">&nbsp;</a> 
```python
make_dimensioned(self, freqs=None, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L2035)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L2035?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a fully dimensioned version of this mode set: not mass-weighted and not frequency-scaled.
  - `freqs`: `np.ndarray | None`
    > the frequencies to compute the frequency-descaling from
  - `masses`: `np.ndarray | None`
    > the masses to compute the mass-unweighting from
  - `:returns`: `MixtureModes`
    > the dimensioned mode set


<a id="Psience.Modes.MixtureModes.MixtureModes.apply_projection" class="docs-object-method">&nbsp;</a> 
```python
apply_projection(self, proj, project_transrot=True, masses=None, origin=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L2067)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L2067?message=Update%20Docs)]
</div>
**LLM Docstring**

Apply an arbitrary projection matrix to the mode/inverse matrices (optionally first combining it with a translation/rotation projector), returning a new mode set restricted to (or excluding) the projected subspace.
  - `proj`: `np.ndarray`
    > the projection matrix to apply
  - `project_transrot`: `bool`
    > whether to additionally combine `proj` with a translation/rotation projector before applying
  - `masses`: `np.ndarray | None`
    > masses to use for the translation/rotation projector and to store on the result
  - `origin`: `np.ndarray | None`
    > reference geometry to use for the translation/rotation projector and to store on the result
  - `:returns`: `MixtureModes`
    > the projected mode set


<a id="Psience.Modes.MixtureModes.MixtureModes.apply_constraints" class="docs-object-method">&nbsp;</a> 
```python
apply_constraints(self, coordinate_constraints, atoms=None, masses=None, origin=None, orthogonal_projection=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L2138)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L2138?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a new mode set with the given internal coordinates constrained (projected out or onto, depending on `orthogonal_projection`), optionally restricted to a subset of atoms, via `apply_projection`.
  - `coordinate_constraints`: `object`
    > the internal coordinate specification(s) to constrain
  - `atoms`: `Iterable[int] | None`
    > restrict the constraint projector(s) to a subset of atoms
  - `masses`: `np.ndarray | None`
    > masses to use instead of `self.masses`
  - `origin`: `np.ndarray | None`
    > reference geometry to use instead of `self.origin` (un-mass-weighted)
  - `orthogonal_projection`: `bool`
    > whether the constraint bases define orthogonal (rather than oblique) projections, and whether multiple constraints are combined by sequential orthogonal projection (`True`) or simple summation (`False`)
  - `:returns`: `MixtureModes`
    > the constrained mode set


<a id="Psience.Modes.MixtureModes.MixtureModes.apply_transformation" class="docs-object-method">&nbsp;</a> 
```python
apply_transformation(self, tf, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L2221)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L2221?message=Update%20Docs)]
</div>
**LLM Docstring**

Apply an arbitrary linear transformation to this mode set, returning the result as a `LocalizedModes` object.
  - `tf`: `np.ndarray | tuple`
    > the transformation to apply, in any form accepted by `LocalizedModes`
  - `opts`: `dict`
    > extra options forwarded to the `LocalizedModes` constructor
  - `:returns`: `LocalizedModes`
    > the transformed modes


<a id="Psience.Modes.MixtureModes.MixtureModes.make_oblique" class="docs-object-method">&nbsp;</a> 
```python
make_oblique(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L2238)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L2238?message=Update%20Docs)]
</div>
**LLM Docstring**

Build an "oblique" mode representation (see `ObliqueModeGenerator`) that makes the mode transformation as close to orthogonal as possible while still diagonalizing the mode-basis Hessian/G-matrix.
  - `:returns`: `LocalizedModes`
    > the oblique-transformed modes
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Modes/MixtureModes/MixtureModes.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Modes/MixtureModes/MixtureModes.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Modes/MixtureModes/MixtureModes.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Modes/MixtureModes/MixtureModes.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes.py#L18?message=Update%20Docs)   
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