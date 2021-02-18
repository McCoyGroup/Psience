## <a id="Psience.BasisReps.Operators.ContractedOperator">ContractedOperator</a>
Provides support for terms that look like `pGp` or `p(dG/dQ)Qp` by
expanding them out as the pure operator component that depends on the basis states (i.e. `pp` or `pQp`)
and doing the appropriate tensor contractions with the expansion coefficients (i.e. `G` or `dG/dQ`)

### Properties and Methods
<a id="Psience.BasisReps.Operators.ContractedOperator.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, coeffs, funcs, quanta, prod_dim=None, axes=None, symmetries=None, selection_rules=None, parallelizer=None, logger=None, zero_threshold=1e-14): 
```

- `coeffs`: `np.ndarray | int`
    >The tensor of coefficients contract with the operator representation (`0` means no term)
- `funcs`: `callable | Iterable[callable]`
    >The functions use to calculate representation
- `quanta`: `int | Iterable[int]`
    >The number of quanta to do the deepest-level calculations up to
- `axes`: `Iterable[int] | None`
    >The axes to use when doing the contractions
- `symmetries`: `Iterable[int] | None`
    >The symmetries to pass through to `Operator`

<a id="Psience.BasisReps.Operators.ContractedOperator.get_elements" class="docs-object-method">&nbsp;</a>
```python
get_elements(self, idx, parallelizer=None): 
```
Computes the operator values over the specified indices
- `idx`: `Iterable[int]`
    >which elements of H0 to compute
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Operators.ContractedOperator.__repr__" class="docs-object-method">&nbsp;</a>
```python
__repr__(self): 
```

### Examples


