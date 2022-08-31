## <a id="Psience.VPT2.Runner.VPTHamiltonianOptions">VPTHamiltonianOptions</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Runner.py#L479)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Runner.py#L479?message=Update%20Docs)]
</div>

Provides a helper to keep track of the levers available for
setting up the Hamiltonian







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.VPT2.Runner.VPTHamiltonianOptions.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, mode_selection=None, include_potential=None, include_gmatrix=None, include_coriolis_coupling=None, include_pseudopotential=None, potential_terms=None, kinetic_terms=None, coriolis_terms=None, pseudopotential_terms=None, dipole_terms=None, dipole_derivatives=None, undimensionalize_normal_modes=None, use_numerical_jacobians=None, eckart_embed_derivatives=None, eckart_embed_planar_ref_tolerance=None, strip_dummy_atoms=None, strip_embedding_coordinates=None, mixed_derivative_handling_mode=None, backpropagate_internals=None, direct_propagate_cartesians=None, zero_mass_term=None, internal_fd_mesh_spacing=None, internal_fd_stencil=None, cartesian_fd_mesh_spacing=None, cartesian_fd_stencil=None, cartesian_analytic_deriv_order=None, internal_by_cartesian_order=None, cartesian_by_internal_order=None, jacobian_warning_threshold=None, check_input_force_constants=None, hessian_tolerance=None, grad_tolerance=None, freq_tolerance=None, g_derivative_threshold=None, gmatrix_tolerance=None, use_internal_modes=None, use_cartesian_kinetic_energy=None, operator_coefficient_threshold=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Runner/VPTHamiltonianOptions.py#L526)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Runner/VPTHamiltonianOptions.py#L526?message=Update%20Docs)]
</div>

  - `mode_selection`: `Iterable[int]|None`
    > the set of the supplied normal modes to do perturbation theory on
  - `include_coriolis_coupling`: `bool`
    > whether or not to include Coriolis coupling in Cartesian normal mode calculation
  - `include_pseudopotential`: `bool`
    > whether or not to include the pseudopotential/Watson term
  - `potential_terms`: `Iterable[np.ndarray]`
    > explicit values for the potential terms (e.g. from analytic models), should be a list of tensors starting with the Hessian with each axis of length `nmodes`
  - `kinetic_terms`: `Iterable[np.ndarray]`
    > explicit values for the kinetic terms (e.g. from analytic models), same format as for the potential
  - `coriolis_terms`: `Iterable[np.ndarray]`
    > explicit values for the Coriolis terms
  - `pseudopotential_terms`: `Iterable[np.ndarray]`
    > explicit values for the psuedopotential terms
  - `undimensionalize_normal_modes`: `bool`
    > whether or not to convert normal modes into dimensional coordinates
  - `use_numerical_jacobians`: `bool`
    > whether or not to use numerical differentiation when getting coordinate transformations
  - `eckart_embed_derivatives`: `bool`
    > whether or not to use Eckart embedding when getting Cartesian to internal transformations (needed for proper results)
  - `strip_dummy_atoms`: `bool`
    > whether or not to strip off dummy atoms when doing transformations
  - `strip_embedding_coordinates`: `bool`
    > whether or not to strip off translation/rotation embedding coordinates when doing transformations
  - `mixed_derivative_handling_mode`: `bool`
    > how to handle differences between numerical/analytical mixed derivatives of potential/dipole terms
  - `backpropagate_internals`: `bool`
    > whether or not to do Cartesian coordinate calculations with values backpropagated from internals
  - `internal_fd_mesh_spacing`: `float`
    > mesh spacing for finite difference of Cartesian coordinates with internals
  - `internal_fd_stencil`: `int`
    > stencil for finite difference of Cartesian coordinates with internals
  - `cartesian_fd_mesh_spacing`: `float`
    > mesh spacing for finite difference of internal coordinates with respect to Cartesians
  - `cartesian_fd_stencil`: `int`
    > stencil for finite difference of internal coordinates with respect to Cartesians
  - `cartesian_analytic_deriv_order`: `int`
    > order of analytic derivatives to use for derivatives of internal coordinates with respect to Cartesians (supports `0` or `1`)
  - `jacobian_warning_threshold`: `float`
    > the value at which to warn that the Jacobian is ill-conditions
  - `check_input_force_constants`: `bool`
    > whether or not to check that the input force constants match the input frequencies
  - `hessian_tolerance`: `float`
    > the deviation to allow when transforming from Cartesian to internal Hessian
  - `grad_tolerance`: `float`
    > the size of the norm of the gradient above which to print a warning
  - `freq_tolerance`: `float`
    > the deviation from the input frequencies to allow when transforming from Cartesians to internals
  - `g_derivative_threshold`: `float`
    > the size of the norm of any G-matrix derivative above which to print a warning
  - `operator_coefficient_threshold`: `float|None`
    > the minimum size of a coefficient to keep when evaluating representation terms
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Runner/VPTHamiltonianOptions.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Runner/VPTHamiltonianOptions.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Runner/VPTHamiltonianOptions.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Runner/VPTHamiltonianOptions.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Runner.py#L479?message=Update%20Docs)   
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