## <a id="Psience.BasisReps.Operators.ContractedOperator">ContractedOperator</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Operators.py#L1063)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Operators.py#L1063?message=Update%20Docs)]
</div>

Provides support for terms that look like `pGp` or `p(dG/dQ)Qp` by
expanding them out as the pure operator component that depends on the basis states (i.e. `pp` or `pQp`)
and doing the appropriate tensor contractions with the expansion coefficients (i.e. `G` or `dG/dQ`)

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

<a id="Psience.BasisReps.Operators.ContractedOperator.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, coeffs, funcs, quanta, prod_dim=None, axes=None, symmetries=None, selection_rules=None, selection_rule_steps=None, zero_threshold=1e-14, chunk_size=None, parallelizer=None, logger=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Operators.py#L1070)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Operators.py#L1070?message=Update%20Docs)]
</div>


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
get_elements(self, idx, parallelizer=None, check_orthogonality=True, memory_constrained=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Operators.py#L1157)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Operators.py#L1157?message=Update%20Docs)]
</div>

Computes the operator values over the specified indices
- `idx`: `Iterable[int]`
    >which elements of H0 to compute
- `:returns`: `_`
    >No description...

<a id="Psience.BasisReps.Operators.ContractedOperator.apply_reduced" class="docs-object-method">&nbsp;</a> 
```python
apply_reduced(self, base_space, parallelizer=None, logger=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Operators.py#L1188)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Operators.py#L1188?message=Update%20Docs)]
</div>

<a id="Psience.BasisReps.Operators.ContractedOperator.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Operators.py#L1231)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Operators.py#L1231?message=Update%20Docs)]
</div>

 </div>
</div>




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/Operators/ContractedOperator.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/Operators/ContractedOperator.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/Operators/ContractedOperator.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/Operators/ContractedOperator.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Operators.py#L1063?message=Update%20Docs)