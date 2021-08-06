## <a id="Psience.VPT2.Runner.VPTHamiltonianOptions">VPTHamiltonianOptions</a>
Provides a helper to keep track of the levers available for
setting up the Hamiltonian

### Properties and Methods
<a id="Psience.VPT2.Runner.VPTHamiltonianOptions.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, coriolis_coupling=None, include_pseudopotential=None, potential_terms=None, kinetic_terms=None, coriolis_terms=None, pseudopotential_terms=None, undimensionalize_normal_modes=None, use_numerical_jacobians=None, eckart_embed_derivatives=None, strip_dummy_atoms=None, strip_embedding_coordinates=None, mixed_derivative_handling_mode=None, backpropagate_internals=None, zero_mass_term=None, internal_fd_mesh_spacing=None, internal_fd_stencil=None, cartesian_fd_mesh_spacing=None, cartesian_fd_stencil=None, cartesian_analytic_deriv_order=None, internal_by_cartesian_order=None, cartesian_by_internal_order=None, jacobian_warning_threshold=None, check_input_force_constants=None, hessian_tolerance=None, grad_tolerance=None, freq_tolerance=None, g_derivative_threshold=None): 
```

- `coriolis_coupling`: `bool`
    >whether or not to include Coriolis coupling in Cartesian normal mode calculation
- `include_pseudopotential`: `bool`
    >whether or not to include the pseudopotential/Watson term
- `potential_terms`: `Iterable[np.ndarray]`
    >explicit values for the potential terms (e.g. from analytic models)
- `kinetic_terms`: `Iterable[np.ndarray]`
    >explicit values for the kinetic terms (e.g. from analytic models)
- `coriolis_terms`: `Iterable[np.ndarray]`
    >explicit values for the Coriolis terms
- `pseudopotential_terms`: `Iterable[np.ndarray]`
    >explicit values for the psuedopotential terms
- `undimensionalize_normal_modes`: `bool`
    >whether or not to convert normal modes into dimensional coordinates
- `use_numerical_jacobians`: `bool`
    >whether or not to use numerical differentiation when getting coordinate transformations
- `eckart_embed_derivatives`: `bool`
    >whether or not to use Eckart embedding when getting Cartesian to internal transformations (needed for proper results)
- `strip_dummy_atoms`: `bool`
    >whether or not to strip off dummy atoms when doing transformations
- `strip_embedding_coordinates`: `bool`
    >whether or not to strip off translation/rotation embedding coordinates when doing transformations
- `mixed_derivative_handling_mode`: `bool`
    >how to handle differences between numerical/analytical mixed derivatives of potential/dipole terms
- `backpropagate_internals`: `bool`
    >whether or not to do Cartesian coordinate calculations with values backpropagated from internals
- `zero_mass_term`: `float`
    >a placeholder value for dummy atom masses
- `internal_fd_mesh_spacing`: `float`
    >mesh spacing for finite difference of Cartesian coordinates with internals
- `internal_fd_stencil`: `int`
    >stencil for finite difference of Cartesian coordinates with internals
- `cartesian_fd_mesh_spacing`: `float`
    >mesh spacing for finite difference of internal coordinates with respect to Cartesians
- `cartesian_fd_stencil`: `int`
    >stencil for finite difference of internal coordinates with respect to Cartesians
- `cartesian_analytic_deriv_order`: `int`
    >order of analytic derivatives to use for derivatives of internal coordinates with respect to Cartesians (supports `0` or `1`)
- `jacobian_warning_threshold`: `float`
    >the value at which to warn that the Jacobian is ill-conditions
- `check_input_force_constants`: `bool`
    >whether or not to check that the input force constants match the input frequencies
- `hessian_tolerance`: `float`
    >the deviation to allow when transforming from Cartesian to internal Hessian
- `grad_tolerance`: `float`
    >the size of the norm of the gradient above which to print a warning
- `freq_tolerance`: `float`
    >the deviation from the input frequencies to allow when transforming from Cartesians to internals
- `g_derivative_threshold`: `float`
    >the size of the norm of any G-matrix derivative above which to print a warning

### Examples




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/edit/ci/examples/ci/docs/Psience/VPT2/Runner/VPTHamiltonianOptions.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/examples/ci/docs/Psience/VPT2/Runner/VPTHamiltonianOptions.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/edit/ci/docs/ci/docs/Psience/VPT2/Runner/VPTHamiltonianOptions.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/docs/templates/ci/docs/Psience/VPT2/Runner/VPTHamiltonianOptions.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Runner.py?message=Update%20Docs)