## <a id="Psience.VPT2.Terms.ExpansionTerms">ExpansionTerms</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L158)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L158?message=Update%20Docs)]
</div>

Base class for kinetic, potential, and dipole derivative terms

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

<a id="Psience.VPT2.Terms.ExpansionTerms.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, molecule, modes=None, mode_selection=None, use_internal_modes=None, logger=None, parallelizer=None, checkpointer=None, undimensionalize=None, numerical_jacobians=True, eckart_embed_derivatives=True, eckart_embed_planar_ref_tolerance=None, strip_dummies=False, strip_embedding=True, mixed_derivative_handling_mode='numerical', backpropagate_internals=False, direct_propagate_cartesians=False, zero_mass_term=10000000.0, internal_fd_mesh_spacing=0.001, internal_fd_stencil=9, cartesian_fd_mesh_spacing=0.001, cartesian_fd_stencil=9, cartesian_analytic_deriv_order=0, internal_by_cartesian_order=3, cartesian_by_internal_order=4, jacobian_warning_threshold=10000.0, coordinate_transformations=None, coordinate_derivatives=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L199)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L199?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L?message=Update%20Docs)]
</div>

Gets the number of atoms (excluding dummies if `strip_dummies` is `True`)
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Terms.ExpansionTerms.modes" class="docs-object-method">&nbsp;</a> 
```python
@property
modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.undimensionalize" class="docs-object-method">&nbsp;</a> 
```python
undimensionalize(self, masses, modes): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L361)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L361?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L395)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L395?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L406)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L406?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L425)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L425?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.get_int_jacobs" class="docs-object-method">&nbsp;</a> 
```python
get_int_jacobs(self, jacs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L448)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L448?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L494)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L494?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L?message=Update%20Docs)]
</div>

Provides the inertial axis frame
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Terms.ExpansionTerms.inertial_frame_derivatives" class="docs-object-method">&nbsp;</a> 
```python
inertial_frame_derivatives(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L558)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L558?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.moment_of_inertia_derivs" class="docs-object-method">&nbsp;</a> 
```python
moment_of_inertia_derivs(self, order): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L600)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L600?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.get_coordinate_transforms" class="docs-object-method">&nbsp;</a> 
```python
get_coordinate_transforms(self, internal_by_cartesian_order=None, cartesian_by_internal_order=None, current_cache=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L801)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L801?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.cartesian_L_matrix" class="docs-object-method">&nbsp;</a> 
```python
@property
cartesian_L_matrix(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.get_cartesians_by_cartesian_modes" class="docs-object-method">&nbsp;</a> 
```python
get_cartesians_by_cartesian_modes(self, order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L1161)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L1161?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.cartesian_L_inverse" class="docs-object-method">&nbsp;</a> 
```python
@property
cartesian_L_inverse(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.get_cartesian_modes_by_cartesians" class="docs-object-method">&nbsp;</a> 
```python
get_cartesian_modes_by_cartesians(self, order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L1182)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L1182?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.internal_L_matrix" class="docs-object-method">&nbsp;</a> 
```python
@property
internal_L_matrix(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.get_internal_modes_by_internals" class="docs-object-method">&nbsp;</a> 
```python
get_internal_modes_by_internals(self, order=None, strip_embedding=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L1200)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L1200?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.internal_L_inverse" class="docs-object-method">&nbsp;</a> 
```python
@property
internal_L_inverse(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.get_internals_by_internal_modes" class="docs-object-method">&nbsp;</a> 
```python
get_internals_by_internal_modes(self, order=None, strip_embedding=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L1225)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L1225?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.cartesians_by_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
cartesians_by_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.get_cartesians_by_modes" class="docs-object-method">&nbsp;</a> 
```python
get_cartesians_by_modes(self, order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L1246)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L1246?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.modes_by_cartesians" class="docs-object-method">&nbsp;</a> 
```python
@property
modes_by_cartesians(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.get_modes_by_cartesians" class="docs-object-method">&nbsp;</a> 
```python
get_modes_by_cartesians(self, order=None, strip_embedding=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L1267)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L1267?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.cartesians_by_internals" class="docs-object-method">&nbsp;</a> 
```python
@property
cartesians_by_internals(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.get_cartesians_by_internals" class="docs-object-method">&nbsp;</a> 
```python
get_cartesians_by_internals(self, order=None, strip_embedding=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L1284)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L1284?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.internals_by_cartesians" class="docs-object-method">&nbsp;</a> 
```python
@property
internals_by_cartesians(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.get_internals_by_cartesians" class="docs-object-method">&nbsp;</a> 
```python
get_internals_by_cartesians(self, order=None, strip_embedding=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L1306)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L1306?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.cartesian_modes_by_internal_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
cartesian_modes_by_internal_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.get_cartesian_modes_by_internal_modes" class="docs-object-method">&nbsp;</a> 
```python
get_cartesian_modes_by_internal_modes(self, order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L1328)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L1328?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.internal_modes_by_cartesian_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
internal_modes_by_cartesian_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.ExpansionTerms.get_internal_modes_by_cartesian_modes" class="docs-object-method">&nbsp;</a> 
```python
get_internal_modes_by_cartesian_modes(self, order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L1347)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L1347?message=Update%20Docs)]
</div>

 </div>
</div>




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Terms/ExpansionTerms.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Terms/ExpansionTerms.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Terms/ExpansionTerms.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Terms/ExpansionTerms.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L158?message=Update%20Docs)