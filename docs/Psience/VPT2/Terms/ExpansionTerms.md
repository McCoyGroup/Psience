## <a id="Psience.VPT2.Terms.ExpansionTerms">ExpansionTerms</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Terms.py#L153)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Terms.py#L153?message=Update%20Docs)]
</div>

Base class for kinetic, potential, and dipole derivative terms

<a id="Psience.VPT2.Terms.ExpansionTerms.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, molecule, modes=None, mode_selection=None, logger=None, parallelizer=None, checkpointer=None, undimensionalize=True, numerical_jacobians=True, eckart_embed_derivatives=True, eckart_embed_planar_ref_tolerance=None, strip_dummies=False, strip_embedding=False, mixed_derivative_handling_mode='numerical', backpropagate_internals=False, direct_propagate_cartesians=False, zero_mass_term=10000000.0, internal_fd_mesh_spacing=0.001, internal_fd_stencil=9, cartesian_fd_mesh_spacing=0.001, cartesian_fd_stencil=9, cartesian_analytic_deriv_order=1, internal_by_cartesian_order=3, cartesian_by_internal_order=4, jacobian_warning_threshold=10000.0, coordinate_transformations=None, coordinate_derivatives=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Terms.py#L194)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Terms.py#L194?message=Update%20Docs)]
</div>


- `molecule`: `Molecule`
    >the molecule we're doing the expansion for
- `modes`: `MolecularVibrations`
    >normal modes in Cartesian coordinates
- `mode_selection`: `None | Iterable[int]`
    >the selection of modes to use
- `undimensionalize`: `bool`
    >whether or not we need to do some units fuckery on the modes

<a id="Psience.VPT2.Terms.ExpansionTerms.num_atoms" class="docs-object-method">&nbsp;</a> 
```python
@property
num_atoms(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Terms.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Terms.py#L?message=Update%20Docs)]
</div>

Gets the number of atoms (excluding dummies if `strip_dummies` is `True`)
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Terms.ExpansionTerms.undimensionalize" class="docs-object-method">&nbsp;</a> 
```python
undimensionalize(self, masses, modes): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Terms.py#L307)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Terms.py#L307?message=Update%20Docs)]
</div>

Removes units from normal modes
- `masses`: `Any`
    >No description...
- `modes`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Terms.ExpansionTerms.get_terms" class="docs-object-method">&nbsp;</a> 
```python
get_terms(self, order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Terms.py#L335)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Terms.py#L335?message=Update%20Docs)]
</div>

Gets the terms up to the given order
- `order`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Terms.ExpansionTerms.get_term" class="docs-object-method">&nbsp;</a> 
```python
get_term(self, t): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Terms.py#L346)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Terms.py#L346?message=Update%20Docs)]
</div>

Provides the term at order `t`
- `t`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Terms.ExpansionTerms.terms" class="docs-object-method">&nbsp;</a> 
```python
@property
terms(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Terms.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Terms.py#L?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Terms.py#L365)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Terms.py#L365?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.get_int_jacobs" class="docs-object-method">&nbsp;</a> 
```python
get_int_jacobs(self, jacs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Terms.py#L388)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Terms.py#L388?message=Update%20Docs)]
</div>

Gets the specified Internal->Cartesian Jacobians
- `jacs`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Terms.ExpansionTerms.get_cart_jacobs" class="docs-object-method">&nbsp;</a> 
```python
get_cart_jacobs(self, jacs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Terms.py#L434)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Terms.py#L434?message=Update%20Docs)]
</div>

Gets the specified Cartesian->Internal Jacobians
- `jacs`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Terms.ExpansionTerms.inertial_frame" class="docs-object-method">&nbsp;</a> 
```python
@property
inertial_frame(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Terms.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Terms.py#L?message=Update%20Docs)]
</div>

Provides the inertial axis frame
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Terms.ExpansionTerms.inertial_frame_derivatives" class="docs-object-method">&nbsp;</a> 
```python
inertial_frame_derivatives(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Terms.py#L498)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Terms.py#L498?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.moment_of_inertia_derivs" class="docs-object-method">&nbsp;</a> 
```python
moment_of_inertia_derivs(self, order): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Terms.py#L540)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Terms.py#L540?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.get_coordinate_transforms" class="docs-object-method">&nbsp;</a> 
```python
get_coordinate_transforms(self, internal_by_cartesian_order=None, cartesian_by_internal_order=None, current_cache=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Terms.py#L741)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Terms.py#L741?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.cartesians_by_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
cartesians_by_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Terms.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Terms.py#L?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.get_cartesians_by_modes" class="docs-object-method">&nbsp;</a> 
```python
get_cartesians_by_modes(self, order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Terms.py#L1029)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Terms.py#L1029?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.modes_by_cartesians" class="docs-object-method">&nbsp;</a> 
```python
@property
modes_by_cartesians(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Terms.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Terms.py#L?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.get_modes_by_cartesians" class="docs-object-method">&nbsp;</a> 
```python
get_modes_by_cartesians(self, order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Terms.py#L1051)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Terms.py#L1051?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.cartesians_by_internals" class="docs-object-method">&nbsp;</a> 
```python
@property
cartesians_by_internals(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Terms.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Terms.py#L?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.get_cartesians_by_internals" class="docs-object-method">&nbsp;</a> 
```python
get_cartesians_by_internals(self, order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Terms.py#L1068)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Terms.py#L1068?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.internals_by_cartesians" class="docs-object-method">&nbsp;</a> 
```python
@property
internals_by_cartesians(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Terms.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Terms.py#L?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.get_internals_by_cartesians" class="docs-object-method">&nbsp;</a> 
```python
get_internals_by_cartesians(self, order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Terms.py#L1085)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Terms.py#L1085?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.cartesian_modes_by_internal_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
cartesian_modes_by_internal_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Terms.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Terms.py#L?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.get_cartesian_modes_by_internal_modes" class="docs-object-method">&nbsp;</a> 
```python
get_cartesian_modes_by_internal_modes(self, order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Terms.py#L1103)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Terms.py#L1103?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.internal_modes_by_cartesian_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
internal_modes_by_cartesian_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Terms.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Terms.py#L?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.get_internal_modes_by_cartesian_modes" class="docs-object-method">&nbsp;</a> 
```python
get_internal_modes_by_cartesian_modes(self, order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Terms.py#L1122)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Terms.py#L1122?message=Update%20Docs)]
</div>



___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/ci/docs/Psience/VPT2/Terms/ExpansionTerms.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/ci/docs/Psience/VPT2/Terms/ExpansionTerms.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/ci/docs/Psience/VPT2/Terms/ExpansionTerms.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/ci/docs/Psience/VPT2/Terms/ExpansionTerms.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Terms.py#L153?message=Update%20Docs)