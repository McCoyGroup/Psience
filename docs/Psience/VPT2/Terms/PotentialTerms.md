## <a id="Psience.VPT2.Terms.PotentialTerms">PotentialTerms</a>
A helper class that can transform the derivatives of the potential from Cartesian to normal coordinates

### Properties and Methods
```python
check_input_force_constants: bool
hessian_tolerance: float
grad_tolerance: float
freq_tolerance: float
```
<a id="Psience.VPT2.Terms.PotentialTerms.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, molecule, mixed_derivs=None, modes=None, potential_derivatives=None, mode_selection=None, logger=None, parallelizer=None, checkpointer=None): 
```

- `molecule`: `Molecule`
    >the molecule that will supply the potential derivatives
- `mixed_derivs`: `bool`
    >whether or not the pulled derivatives are partially derivatives along the normal coords
- `modes`: `None | MolecularVibrations`
    >the normal modes to use when doing calculations
- `mode_selection`: `None | Iterable[int]`
    >the subset of normal modes to use

<a id="Psience.VPT2.Terms.PotentialTerms.get_terms" class="docs-object-method">&nbsp;</a>
```python
get_terms(self, order=None, logger=None): 
```

### Examples


