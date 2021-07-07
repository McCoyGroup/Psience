## <a id="Psience.VPT2.Terms.ExpansionTerms">ExpansionTerms</a>
Base class for kinetic, potential, and dipole derivative terms

### Properties and Methods
```python
internal_fd_mesh_spacing: float
internal_fd_stencil: int
cartesian_fd_mesh_spacing: float
cartesian_fd_stencil: int
cartesian_analytic_deriv_order: int
internal_by_cartesian_order: int
cartesian_by_internal_order: int
jacobian_warning_threshold: float
```
<a id="Psience.VPT2.Terms.ExpansionTerms.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, molecule, modes=None, mode_selection=None, undimensionalize=True, logger=None, parallelizer=None, checkpointer=None, numerical_jacobians=True, eckart_embed=True): 
```

- `molecule`: `Molecule`
    >the molecule we're doing the expansion for
- `modes`: `MolecularVibrations`
    >normal modes in Cartesian coordinates
- `mode_selection`: `None | Iterable[int]`
    >the selection of modes to use
- `undimensionalize`: `bool`
    >whether or not we need to do some units fuckery on the modes

<a id="Psience.VPT2.Terms.ExpansionTerms.undimensionalize" class="docs-object-method">&nbsp;</a>
```python
undimensionalize(self, masses, modes): 
```

<a id="Psience.VPT2.Terms.ExpansionTerms.get_terms" class="docs-object-method">&nbsp;</a>
```python
get_terms(self, order=None): 
```

<a id="Psience.VPT2.Terms.ExpansionTerms.get_term" class="docs-object-method">&nbsp;</a>
```python
get_term(self, t): 
```

<a id="Psience.VPT2.Terms.ExpansionTerms.terms" class="docs-object-method">&nbsp;</a>
```python
@property
terms(self): 
```

<a id="Psience.VPT2.Terms.ExpansionTerms.__getitem__" class="docs-object-method">&nbsp;</a>
```python
__getitem__(self, item): 
```

<a id="Psience.VPT2.Terms.ExpansionTerms.get_int_jacobs" class="docs-object-method">&nbsp;</a>
```python
get_int_jacobs(self, jacs): 
```

<a id="Psience.VPT2.Terms.ExpansionTerms.get_cart_jacobs" class="docs-object-method">&nbsp;</a>
```python
get_cart_jacobs(self, jacs): 
```

<a id="Psience.VPT2.Terms.ExpansionTerms.inertial_frame" class="docs-object-method">&nbsp;</a>
```python
@property
inertial_frame(self): 
```

<a id="Psience.VPT2.Terms.ExpansionTerms.inertial_frame_derivatives" class="docs-object-method">&nbsp;</a>
```python
inertial_frame_derivatives(self): 
```

<a id="Psience.VPT2.Terms.ExpansionTerms.get_coordinate_transforms" class="docs-object-method">&nbsp;</a>
```python
get_coordinate_transforms(self, internal_by_cartesian_order=None, cartesian_by_internal_order=None, current_cache=None): 
```

<a id="Psience.VPT2.Terms.ExpansionTerms.cartesians_by_modes" class="docs-object-method">&nbsp;</a>
```python
@property
cartesians_by_modes(self): 
```

<a id="Psience.VPT2.Terms.ExpansionTerms.get_cartesians_by_modes" class="docs-object-method">&nbsp;</a>
```python
get_cartesians_by_modes(self, order=None): 
```

<a id="Psience.VPT2.Terms.ExpansionTerms.modes_by_cartesians" class="docs-object-method">&nbsp;</a>
```python
@property
modes_by_cartesians(self): 
```

<a id="Psience.VPT2.Terms.ExpansionTerms.get_modes_by_cartesians" class="docs-object-method">&nbsp;</a>
```python
get_modes_by_cartesians(self, order=None): 
```

<a id="Psience.VPT2.Terms.ExpansionTerms.cartesians_by_internals" class="docs-object-method">&nbsp;</a>
```python
@property
cartesians_by_internals(self): 
```

<a id="Psience.VPT2.Terms.ExpansionTerms.get_cartesians_by_internals" class="docs-object-method">&nbsp;</a>
```python
get_cartesians_by_internals(self, order=None): 
```

<a id="Psience.VPT2.Terms.ExpansionTerms.internals_by_cartesians" class="docs-object-method">&nbsp;</a>
```python
@property
internals_by_cartesians(self): 
```

<a id="Psience.VPT2.Terms.ExpansionTerms.get_internals_by_cartesians" class="docs-object-method">&nbsp;</a>
```python
get_internals_by_cartesians(self, order=None): 
```

<a id="Psience.VPT2.Terms.ExpansionTerms.cartesian_modes_by_internal_modes" class="docs-object-method">&nbsp;</a>
```python
@property
cartesian_modes_by_internal_modes(self): 
```

<a id="Psience.VPT2.Terms.ExpansionTerms.get_cartesian_modes_by_internal_modes" class="docs-object-method">&nbsp;</a>
```python
get_cartesian_modes_by_internal_modes(self, order=None): 
```

### Examples


