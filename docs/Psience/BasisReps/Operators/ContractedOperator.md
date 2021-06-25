## <a id="Psience.BasisReps.Operators.ContractedOperator">ContractedOperator</a>
Provides support for terms that look like `pGp` or `p(dG/dQ)Qp` by
expanding them out as the pure operator component that depends on the basis states (i.e. `pp` or `pQp`)
and doing the appropriate tensor contractions with the expansion coefficients (i.e. `G` or `dG/dQ`)

### Properties and Methods
<a id="Psience.BasisReps.Operators.ContractedOperator.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, coeffs, funcs, quanta, prod_dim=None, axes=None, symmetries=None, selection_rules=None, zero_threshold=1e-14, chunk_size=None, parallelizer=None, logger=None): 
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
- `prod_dim`: `Any`
    >No description...
- `selection_rules`: `Any`
    >No description...
- `parallelizer`: `Any`
    >No description...
- `logger`: `Any`
    >No description...
- `zero_threshold`: `Any`
    >No description...
- `chunk_size`: `int | None`
    >number of elements that can be evaluated at once (for memory reasons)

<a id="Psience.BasisReps.Operators.ContractedOperator.get_elements" class="docs-object-method">&nbsp;</a>
```python
get_elements(self, idx, parallelizer=None, check_orthogonality=True): 
```
Computes the operator values over the specified indices
- `idx`: `Iterable[int]`
    >which elements of H0 to compute
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Operators.ContractedOperator.apply_reduced" class="docs-object-method">&nbsp;</a>
```python
apply_reduced(self, base_space, parallelizer=None, logger=None): 
```

<a id="Psience.BasisReps.Operators.ContractedOperator.__repr__" class="docs-object-method">&nbsp;</a>
```python
__repr__(self): 
```

### Examples


