## <a id="Psience.VPT2.Terms.ExpansionTerms">ExpansionTerms</a>
Base class for kinetic, potential, and dipole derivative terms

### Properties and Methods
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
get_terms(self): 
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

<a id="Psience.VPT2.Terms.ExpansionTerms.get_coordinate_transforms" class="docs-object-method">&nbsp;</a>
```python
get_coordinate_transforms(self): 
```

<a id="Psience.VPT2.Terms.ExpansionTerms.cartesians_by_modes" class="docs-object-method">&nbsp;</a>
```python
@property
cartesians_by_modes(self): 
```

<a id="Psience.VPT2.Terms.ExpansionTerms.modes_by_cartesians" class="docs-object-method">&nbsp;</a>
```python
@property
modes_by_cartesians(self): 
```

<a id="Psience.VPT2.Terms.ExpansionTerms.cartesians_by_internals" class="docs-object-method">&nbsp;</a>
```python
@property
cartesians_by_internals(self): 
```

<a id="Psience.VPT2.Terms.ExpansionTerms.internals_by_cartesians" class="docs-object-method">&nbsp;</a>
```python
@property
internals_by_cartesians(self): 
```

<a id="Psience.VPT2.Terms.ExpansionTerms.cartesian_modes_by_internal_modes" class="docs-object-method">&nbsp;</a>
```python
@property
cartesian_modes_by_internal_modes(self): 
```

### Examples


