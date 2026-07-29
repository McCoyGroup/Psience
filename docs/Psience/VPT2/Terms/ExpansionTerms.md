## <a id="Psience.VPT2.Terms.ExpansionTerms">ExpansionTerms</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L324)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L324?message=Update%20Docs)]
</div>

Base class for kinetic, potential, and dipole derivative terms







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.VPT2.Terms.ExpansionTerms.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, molecule, modes=None, mode_selection=None, mode_transformation=None, use_internal_modes=None, logger=None, parallelizer=None, checkpointer=None, undimensionalize=None, numerical_jacobians=True, eckart_embed_derivatives=True, eckart_embed_planar_ref_tolerance=None, strip_dummies=False, strip_embedding=True, mixed_derivative_handling_mode=None, mixed_derivative_warning_threshold=0.00025, mixed_derivative_handle_zeros=False, backpropagate_internals=False, direct_propagate_cartesians=False, zero_mass_term=10000000.0, expansion_handling_mode='old', internal_fd_mesh_spacing=0.01, internal_fd_stencil=None, cartesian_fd_mesh_spacing=0.01, cartesian_fd_stencil=None, cartesian_analytic_deriv_order=None, cartesian_by_internal_derivative_method=None, internal_by_cartesian_order=3, cartesian_by_internal_order=4, jacobian_warning_threshold=10000.0, coordinate_transformations=None, coordinate_derivatives=None, imaginary_frequency_handling_mode='abs'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L370)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L370?message=Update%20Docs)]
</div>

  - `molecule`: `Molecule`
    > the molecule we're doing the expansion for
  - `modes`: `MolecularVibrations`
    > normal modes in Cartesian coordinates
  - `mode_selection`: `None | Iterable[int]`
    > the selection of modes to use
  - `undimensionalize`: `bool`
    > whether or not we need to do some units fuckery on the modes


<a id="Psience.VPT2.Terms.ExpansionTerms.num_atoms" class="docs-object-method">&nbsp;</a> 
```python
@property
num_atoms(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L520)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L520?message=Update%20Docs)]
</div>
Gets the number of atoms (excluding dummies if `strip_dummies` is `True`)
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Terms.ExpansionTerms.modes" class="docs-object-method">&nbsp;</a> 
```python
@property
modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L592)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L592?message=Update%20Docs)]
</div>
**LLM Docstring**

The stored mode object (in whatever basis -- Cartesian or internal -- it was constructed with).
  - `:returns`: `MixtureModes`
    > the mode object


<a id="Psience.VPT2.Terms.ExpansionTerms.get_terms" class="docs-object-method">&nbsp;</a> 
```python
get_terms(self, order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L655)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L655?message=Update%20Docs)]
</div>
Gets the terms up to the given order
  - `order`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Terms.ExpansionTerms.get_term" class="docs-object-method">&nbsp;</a> 
```python
get_term(self, t): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L666)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L666?message=Update%20Docs)]
</div>
Provides the term at order `t`
  - `t`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Terms.ExpansionTerms.terms" class="docs-object-method">&nbsp;</a> 
```python
@property
terms(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L679)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L679?message=Update%20Docs)]
</div>
**LLM Docstring**

The (cached) full set of expansion terms, computed lazily via `get_terms()` the first time they're needed.
  - `:returns`: `list[np.ndarray]`
    > the expansion terms


<a id="Psience.VPT2.Terms.ExpansionTerms.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L693)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L693?message=Update%20Docs)]
</div>
**LLM Docstring**

Fetch a single term at the given order, via `get_term`.
  - `item`: `int`
    > the order to fetch
  - `:returns`: `np.ndarray`
    > the term at that order


<a id="Psience.VPT2.Terms.ExpansionTerms.get_int_jacobs" class="docs-object-method">&nbsp;</a> 
```python
get_int_jacobs(self, jacs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L762)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L762?message=Update%20Docs)]
</div>
Gets the specified Internal->Cartesian Jacobians
  - `jacs`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Terms.ExpansionTerms.get_cart_jacobs" class="docs-object-method">&nbsp;</a> 
```python
get_cart_jacobs(self, jacs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L813)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L813?message=Update%20Docs)]
</div>
Gets the specified Cartesian->Internal Jacobians
  - `jacs`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Terms.ExpansionTerms.inertial_frame" class="docs-object-method">&nbsp;</a> 
```python
@property
inertial_frame(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L859)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L859?message=Update%20Docs)]
</div>
Provides the inertial axis frame
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Terms.ExpansionTerms.inertial_frame_derivatives" class="docs-object-method">&nbsp;</a> 
```python
inertial_frame_derivatives(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L881)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L881?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the first and second derivatives of the (mass-weighted) inertia tensor with respect to mass-weighted Cartesian displacements, using closed-form tensor expressions rather than finite differences.
  - `:returns`: `list[np.ndarray]`
    > `[I0Y, I0YY]`, the first derivative tensor (shape `(3*nAt, 3, 3)`) and second derivative tensor (shape `(3*nAt, 3*nAt, 3, 3)`) of the inertia tensor with respect to mass-weighted Cartesian displacements


<a id="Psience.VPT2.Terms.ExpansionTerms.moment_of_inertia_derivs" class="docs-object-method">&nbsp;</a> 
```python
moment_of_inertia_derivs(self, order): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L931)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L931?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the Taylor-series derivatives of the inverse inertia tensor with respect to the normal-mode coordinates, up to the requested order, via the recursive relation built from the first-order Cartesian inertia-tensor derivative re-expressed in mode coordinates.
  - `order`: `int`
    > the highest derivative order to compute
  - `:returns`: `list[np.ndarray]`
    > the reciprocal-inertia-tensor expansion terms, `order + 1` entries starting from the inertia tensor itself


<a id="Psience.VPT2.Terms.ExpansionTerms.get_coordinate_transforms" class="docs-object-method">&nbsp;</a> 
```python
get_coordinate_transforms(self, internal_by_cartesian_order=None, cartesian_by_internal_order=None, current_cache=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L983)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L983?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute (and cache, per-molecule, both in memory and via the checkpointer) the full set of Jacobians relating Cartesian coordinates, internal coordinates, Cartesian normal modes, and internal-coordinate-basis normal modes to each other, up to the requested derivative orders: computes the internals-by-Cartesians Jacobians (mass-weighting them and warning about/zeroing any anomalously large entries), then chains them together with the mode transformation via `TensorDerivativeConverter` to populate every entry of `JacobianKeys`.
  - `internal_by_cartesian_order`: `int | None`
    > derivative order (number of Cartesian derivatives) to compute for internals-by-Cartesians Jacobians; defaults to `self.internal_by_cartesian_order`
  - `cartesian_by_internal_order`: `int | None`
    > derivative order (number of internal derivatives) to compute for Cartesians-by-internals Jacobians; defaults to `self.cartesian_by_internal_order`
  - `current_cache`: `dict | None`
    > an existing partial cache to extend, instead of the per-molecule cache (in-memory or checkpointed) that would otherwise be looked up
  - `:returns`: `dict`
    > the (possibly newly extended) cache mapping each `JacobianKeys` member to a list of Jacobian tensors by order


<a id="Psience.VPT2.Terms.ExpansionTerms.cartesian_L_matrix" class="docs-object-method">&nbsp;</a> 
```python
@property
cartesian_L_matrix(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1410)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1410?message=Update%20Docs)]
</div>
**LLM Docstring**

First-order Cartesians-by-Cartesian-normal-modes transformation matrix.
  - `:returns`: `np.ndarray`
    > the leading term of `get_cartesians_by_cartesian_modes(1)`


<a id="Psience.VPT2.Terms.ExpansionTerms.get_cartesians_by_cartesian_modes" class="docs-object-method">&nbsp;</a> 
```python
get_cartesians_by_cartesian_modes(self, order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1421)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1421?message=Update%20Docs)]
</div>
**LLM Docstring**

Fetch the Cartesians-by-Cartesian-normal-modes Jacobians up to the requested order, computing them (via `get_coordinate_transforms`) if not already cached.
  - `order`: `int | None`
    > number of derivative orders to return; if `None`, all currently cached orders are returned
  - `:returns`: `list[np.ndarray]`
    > list of Jacobian tensors, one per derivative order


<a id="Psience.VPT2.Terms.ExpansionTerms.cartesian_L_inverse" class="docs-object-method">&nbsp;</a> 
```python
@property
cartesian_L_inverse(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1446)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1446?message=Update%20Docs)]
</div>
**LLM Docstring**

First-order Cartesian-normal-modes-by-Cartesians transformation matrix.
  - `:returns`: `np.ndarray`
    > the leading term of `get_cartesian_modes_by_cartesians(1)`


<a id="Psience.VPT2.Terms.ExpansionTerms.get_cartesian_modes_by_cartesians" class="docs-object-method">&nbsp;</a> 
```python
get_cartesian_modes_by_cartesians(self, order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1457)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1457?message=Update%20Docs)]
</div>
**LLM Docstring**

Fetch the Cartesian-normal-modes-by-Cartesians Jacobians up to the requested order, computing them if not already cached.
  - `order`: `int | None`
    > number of derivative orders to return; if `None`, all currently cached orders are returned
  - `:returns`: `list[np.ndarray]`
    > list of Jacobian tensors, one per derivative order


<a id="Psience.VPT2.Terms.ExpansionTerms.internal_L_matrix" class="docs-object-method">&nbsp;</a> 
```python
@property
internal_L_matrix(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1483)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1483?message=Update%20Docs)]
</div>
**LLM Docstring**

First-order internal-normal-modes-by-internals transformation matrix.
  - `:returns`: `np.ndarray`
    > the leading term of `get_internal_modes_by_internals(1)`


<a id="Psience.VPT2.Terms.ExpansionTerms.get_internal_modes_by_internals" class="docs-object-method">&nbsp;</a> 
```python
get_internal_modes_by_internals(self, order=None, strip_embedding=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1494)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1494?message=Update%20Docs)]
</div>
**LLM Docstring**

Fetch the internal-normal-modes-by-internals Jacobians up to the requested order, optionally stripping embedding-coordinate rows from the result.
  - `order`: `int | None`
    > number of derivative orders to return; if `None`, all currently cached orders are returned
  - `strip_embedding`: `bool`
    > whether to strip embedding-coordinate rows from the result (only applied if not already stripped globally via `self.strip_embedding`)
  - `:returns`: `list[np.ndarray]`
    > list of Jacobian tensors, one per derivative order


<a id="Psience.VPT2.Terms.ExpansionTerms.internal_L_inverse" class="docs-object-method">&nbsp;</a> 
```python
@property
internal_L_inverse(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1530)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1530?message=Update%20Docs)]
</div>
**LLM Docstring**

First-order internals-by-internal-normal-modes transformation matrix.
  - `:returns`: `np.ndarray`
    > the leading term of `get_internals_by_internal_modes(1)`


<a id="Psience.VPT2.Terms.ExpansionTerms.get_internals_by_internal_modes" class="docs-object-method">&nbsp;</a> 
```python
get_internals_by_internal_modes(self, order=None, strip_embedding=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1541)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1541?message=Update%20Docs)]
</div>
**LLM Docstring**

Fetch the internals-by-internal-normal-modes Jacobians up to the requested order, optionally stripping embedding-coordinate columns from the result.
  - `order`: `int | None`
    > number of derivative orders to return; if `None`, all currently cached orders are returned
  - `strip_embedding`: `bool`
    > whether to strip embedding-coordinate columns from the result (only applied if not already stripped globally via `self.strip_embedding`)
  - `:returns`: `list[np.ndarray]`
    > list of Jacobian tensors, one per derivative order


<a id="Psience.VPT2.Terms.ExpansionTerms.cartesians_by_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
cartesians_by_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1572)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1572?message=Update%20Docs)]
</div>
**LLM Docstring**

All cached Cartesians-by-internal-modes Jacobians, computing the default set if not already cached.
  - `:returns`: `list[np.ndarray]`
    > the `JacobianKeys.CartesiansByInternalModes` entry from `get_cartesians_by_modes()`


<a id="Psience.VPT2.Terms.ExpansionTerms.get_cartesians_by_modes" class="docs-object-method">&nbsp;</a> 
```python
get_cartesians_by_modes(self, order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1583)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1583?message=Update%20Docs)]
</div>
**LLM Docstring**

Fetch the Cartesians-by-internal-normal-modes Jacobians up to the requested order, computing them if not already cached.
  - `order`: `int | None`
    > number of derivative orders to return; if `None`, all currently cached orders are returned
  - `:returns`: `list[np.ndarray]`
    > list of Jacobian tensors, one per derivative order


<a id="Psience.VPT2.Terms.ExpansionTerms.modes_by_cartesians" class="docs-object-method">&nbsp;</a> 
```python
@property
modes_by_cartesians(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1612)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1612?message=Update%20Docs)]
</div>
**LLM Docstring**

All cached internal-normal-modes-by-Cartesians Jacobians, computing the default set if not already cached.
  - `:returns`: `list[np.ndarray]`
    > the `JacobianKeys.InternalModesByCartesians` entry from `get_coordinate_transforms()`


<a id="Psience.VPT2.Terms.ExpansionTerms.get_modes_by_cartesians" class="docs-object-method">&nbsp;</a> 
```python
get_modes_by_cartesians(self, order=None, strip_embedding=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1623)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1623?message=Update%20Docs)]
</div>
**LLM Docstring**

Fetch the internal-normal-modes-by-Cartesians Jacobians up to the requested order.
  - `order`: `int | None`
    > number of derivative orders to return; if `None`, all currently cached orders are returned
  - `strip_embedding`: `bool`
    > accepted for interface consistency with sibling methods but not used in this method's body
  - `:returns`: `list[np.ndarray]`
    > list of Jacobian tensors, one per derivative order


<a id="Psience.VPT2.Terms.ExpansionTerms.cartesians_by_internals" class="docs-object-method">&nbsp;</a> 
```python
@property
cartesians_by_internals(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1650)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1650?message=Update%20Docs)]
</div>
**LLM Docstring**

All cached Cartesians-by-internals Jacobians, computing the default set if not already cached.
  - `:returns`: `list[np.ndarray]`
    > the `JacobianKeys.CartesiansByInternals` entry from `get_coordinate_transforms()`


<a id="Psience.VPT2.Terms.ExpansionTerms.get_cartesians_by_internals" class="docs-object-method">&nbsp;</a> 
```python
get_cartesians_by_internals(self, order=None, strip_embedding=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1661)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1661?message=Update%20Docs)]
</div>
**LLM Docstring**

Fetch the Cartesians-by-internals Jacobians up to the requested order, optionally stripping embedding coordinates from every axis but the first.
  - `order`: `int | None`
    > number of derivative orders to return; if `None`, all currently cached orders are returned
  - `strip_embedding`: `bool`
    > whether to strip embedding coordinates from the trailing axes of the result
  - `:returns`: `list[np.ndarray]`
    > list of Jacobian tensors, one per derivative order


<a id="Psience.VPT2.Terms.ExpansionTerms.internals_by_cartesians" class="docs-object-method">&nbsp;</a> 
```python
@property
internals_by_cartesians(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1693)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1693?message=Update%20Docs)]
</div>
**LLM Docstring**

All cached internals-by-Cartesians Jacobians, computing the default set if not already cached.
  - `:returns`: `list[np.ndarray]`
    > the `JacobianKeys.InternalsByCartesians` entry from `get_coordinate_transforms()`


<a id="Psience.VPT2.Terms.ExpansionTerms.get_internals_by_cartesians" class="docs-object-method">&nbsp;</a> 
```python
get_internals_by_cartesians(self, order=None, strip_embedding=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1704)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1704?message=Update%20Docs)]
</div>
**LLM Docstring**

Fetch the internals-by-Cartesians Jacobians up to the requested order, optionally stripping embedding coordinates from the trailing axis.
  - `order`: `int | None`
    > number of derivative orders to return; if `None`, all currently cached orders are returned
  - `strip_embedding`: `bool`
    > whether to strip embedding coordinates from the trailing axis of the result
  - `:returns`: `list[np.ndarray]`
    > list of Jacobian tensors, one per derivative order


<a id="Psience.VPT2.Terms.ExpansionTerms.cartesian_modes_by_internal_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
cartesian_modes_by_internal_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1736)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1736?message=Update%20Docs)]
</div>
**LLM Docstring**

All cached Cartesian-normal-modes-by-internal-normal-modes Jacobians, computing the default set if not already cached.
  - `:returns`: `list[np.ndarray]`
    > the `JacobianKeys.CartesianModesByInternalModes` entry from `get_coordinate_transforms()`


<a id="Psience.VPT2.Terms.ExpansionTerms.get_cartesian_modes_by_internal_modes" class="docs-object-method">&nbsp;</a> 
```python
get_cartesian_modes_by_internal_modes(self, order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1747)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1747?message=Update%20Docs)]
</div>
**LLM Docstring**

Fetch the Cartesian-normal-modes-by-internal-normal-modes Jacobians up to the requested order.
  - `order`: `int | None`
    > number of derivative orders to return; if `None`, all currently cached orders are returned
  - `:returns`: `list[np.ndarray]`
    > list of Jacobian tensors, one per derivative order


<a id="Psience.VPT2.Terms.ExpansionTerms.internal_modes_by_cartesian_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
internal_modes_by_cartesian_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1773)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1773?message=Update%20Docs)]
</div>
**LLM Docstring**

All cached internal-normal-modes-by-Cartesian-normal-modes Jacobians, computing the default set if not already cached.
  - `:returns`: `list[np.ndarray]`
    > the `JacobianKeys.InternalModesByCartesianModes` entry from `get_coordinate_transforms()`


<a id="Psience.VPT2.Terms.ExpansionTerms.get_internal_modes_by_cartesian_modes" class="docs-object-method">&nbsp;</a> 
```python
get_internal_modes_by_cartesian_modes(self, order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1785)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1785?message=Update%20Docs)]
</div>
**LLM Docstring**

Fetch the internal-normal-modes-by-Cartesian-normal-modes Jacobians up to the requested order.
  - `order`: `int | None`
    > number of derivative orders to return; if `None`, all currently cached orders are returned
  - `:returns`: `list[np.ndarray]`
    > list of Jacobian tensors, one per derivative order
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Terms/ExpansionTerms.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Terms/ExpansionTerms.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Terms/ExpansionTerms.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Terms/ExpansionTerms.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L324?message=Update%20Docs)   
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