## <a id="Psience.VPT2.Terms.ExpansionTerms">ExpansionTerms</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L161)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L161?message=Update%20Docs)]
</div>

Base class for kinetic, potential, and dipole derivative terms







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.VPT2.Terms.ExpansionTerms.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, molecule, modes=None, mode_selection=None, mode_transformation=None, use_internal_modes=None, logger=None, parallelizer=None, checkpointer=None, undimensionalize=None, numerical_jacobians=True, eckart_embed_derivatives=True, eckart_embed_planar_ref_tolerance=None, strip_dummies=False, strip_embedding=True, mixed_derivative_handling_mode='old', backpropagate_internals=False, direct_propagate_cartesians=False, zero_mass_term=10000000.0, internal_fd_mesh_spacing=0.001, internal_fd_stencil=None, cartesian_fd_mesh_spacing=0.01, cartesian_fd_stencil=None, cartesian_analytic_deriv_order=0, internal_by_cartesian_derivative_method='old', internal_by_cartesian_order=3, cartesian_by_internal_order=4, jacobian_warning_threshold=10000.0, coordinate_transformations=None, coordinate_derivatives=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L203)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L203?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L327)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L327?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L379)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L379?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Terms.ExpansionTerms.get_terms" class="docs-object-method">&nbsp;</a> 
```python
get_terms(self, order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L424)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L424?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L435)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L435?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L448)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L448?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Terms.ExpansionTerms.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L454)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L454?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Terms.ExpansionTerms.get_int_jacobs" class="docs-object-method">&nbsp;</a> 
```python
get_int_jacobs(self, jacs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L477)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L477?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L528)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L528?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L574)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L574?message=Update%20Docs)]
</div>
Provides the inertial axis frame
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Terms.ExpansionTerms.inertial_frame_derivatives" class="docs-object-method">&nbsp;</a> 
```python
inertial_frame_derivatives(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L596)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L596?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Terms.ExpansionTerms.moment_of_inertia_derivs" class="docs-object-method">&nbsp;</a> 
```python
moment_of_inertia_derivs(self, order): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L638)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L638?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Terms.ExpansionTerms.get_coordinate_transforms" class="docs-object-method">&nbsp;</a> 
```python
get_coordinate_transforms(self, internal_by_cartesian_order=None, cartesian_by_internal_order=None, current_cache=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L672)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L672?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Terms.ExpansionTerms.cartesian_L_matrix" class="docs-object-method">&nbsp;</a> 
```python
@property
cartesian_L_matrix(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1071)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1071?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Terms.ExpansionTerms.get_cartesians_by_cartesian_modes" class="docs-object-method">&nbsp;</a> 
```python
get_cartesians_by_cartesian_modes(self, order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1074)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1074?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Terms.ExpansionTerms.cartesian_L_inverse" class="docs-object-method">&nbsp;</a> 
```python
@property
cartesian_L_inverse(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1088)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1088?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Terms.ExpansionTerms.get_cartesian_modes_by_cartesians" class="docs-object-method">&nbsp;</a> 
```python
get_cartesian_modes_by_cartesians(self, order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1091)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1091?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Terms.ExpansionTerms.internal_L_matrix" class="docs-object-method">&nbsp;</a> 
```python
@property
internal_L_matrix(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1106)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1106?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Terms.ExpansionTerms.get_internal_modes_by_internals" class="docs-object-method">&nbsp;</a> 
```python
get_internal_modes_by_internals(self, order=None, strip_embedding=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1109)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1109?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Terms.ExpansionTerms.internal_L_inverse" class="docs-object-method">&nbsp;</a> 
```python
@property
internal_L_inverse(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1132)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1132?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Terms.ExpansionTerms.get_internals_by_internal_modes" class="docs-object-method">&nbsp;</a> 
```python
get_internals_by_internal_modes(self, order=None, strip_embedding=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1135)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1135?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Terms.ExpansionTerms.cartesians_by_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
cartesians_by_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1153)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1153?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Terms.ExpansionTerms.get_cartesians_by_modes" class="docs-object-method">&nbsp;</a> 
```python
get_cartesians_by_modes(self, order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1156)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1156?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Terms.ExpansionTerms.modes_by_cartesians" class="docs-object-method">&nbsp;</a> 
```python
@property
modes_by_cartesians(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1174)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1174?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Terms.ExpansionTerms.get_modes_by_cartesians" class="docs-object-method">&nbsp;</a> 
```python
get_modes_by_cartesians(self, order=None, strip_embedding=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1177)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1177?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Terms.ExpansionTerms.cartesians_by_internals" class="docs-object-method">&nbsp;</a> 
```python
@property
cartesians_by_internals(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1191)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1191?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Terms.ExpansionTerms.get_cartesians_by_internals" class="docs-object-method">&nbsp;</a> 
```python
get_cartesians_by_internals(self, order=None, strip_embedding=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1194)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1194?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Terms.ExpansionTerms.internals_by_cartesians" class="docs-object-method">&nbsp;</a> 
```python
@property
internals_by_cartesians(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1213)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1213?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Terms.ExpansionTerms.get_internals_by_cartesians" class="docs-object-method">&nbsp;</a> 
```python
get_internals_by_cartesians(self, order=None, strip_embedding=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1216)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1216?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Terms.ExpansionTerms.cartesian_modes_by_internal_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
cartesian_modes_by_internal_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1235)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1235?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Terms.ExpansionTerms.get_cartesian_modes_by_internal_modes" class="docs-object-method">&nbsp;</a> 
```python
get_cartesian_modes_by_internal_modes(self, order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1238)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1238?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Terms.ExpansionTerms.internal_modes_by_cartesian_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
internal_modes_by_cartesian_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1253)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1253?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Terms.ExpansionTerms.get_internal_modes_by_cartesian_modes" class="docs-object-method">&nbsp;</a> 
```python
get_internal_modes_by_cartesian_modes(self, order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/ExpansionTerms.py#L1257)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/ExpansionTerms.py#L1257?message=Update%20Docs)]
</div>
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
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L161?message=Update%20Docs)   
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