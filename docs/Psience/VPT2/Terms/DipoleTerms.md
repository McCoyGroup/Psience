## <a id="Psience.VPT2.Terms.DipoleTerms">DipoleTerms</a>


### Properties and Methods
<a id="Psience.VPT2.Terms.DipoleTerms.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, molecule, derivatives=None, mixed_derivs=None, modes=None, mode_selection=None, logger=None, parallelizer=None, checkpointer=None): 
```

- `molecule`: `Molecule`
    >the molecule that will supply the dipole derivatives
- `mixed_derivs`: `bool`
    >whether or not the pulled derivatives are partially derivatives along the normal coords
- `modes`: `None | MolecularNormalModes`
    >the normal modes to use when doing calculations
- `mode_selection`: `None | Iterable[int]`
    >the subset of normal modes to use

<a id="Psience.VPT2.Terms.DipoleTerms.get_terms" class="docs-object-method">&nbsp;</a>
```python
get_terms(self, order=None): 
```

### Examples


