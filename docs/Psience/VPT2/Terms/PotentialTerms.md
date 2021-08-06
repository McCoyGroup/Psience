## <a id="Psience.VPT2.Terms.PotentialTerms">PotentialTerms</a>
A helper class that can transform the derivatives of the potential from Cartesian to normal coordinates

### Properties and Methods
<a id="Psience.VPT2.Terms.PotentialTerms.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, molecule, mixed_derivs=None, modes=None, potential_derivatives=None, mode_selection=None, logger=None, parallelizer=None, checkpointer=None, check_input_force_constants=True, hessian_tolerance=0.0001, grad_tolerance=0.0001, freq_tolerance=0.002, **opts): 
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




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/edit/ci/examples/ci/docs/Psience/VPT2/Terms/PotentialTerms.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/examples/ci/docs/Psience/VPT2/Terms/PotentialTerms.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/edit/ci/docs/ci/docs/Psience/VPT2/Terms/PotentialTerms.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/docs/templates/ci/docs/Psience/VPT2/Terms/PotentialTerms.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Terms.py?message=Update%20Docs)