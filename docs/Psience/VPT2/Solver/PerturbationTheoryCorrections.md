## <a id="Psience.VPT2.Solver.PerturbationTheoryCorrections">PerturbationTheoryCorrections</a>
Represents a set of corrections from perturbation theory.
Can be used to correct other operators in the basis of the original calculation.

### Properties and Methods
<a id="Psience.VPT2.Solver.PerturbationTheoryCorrections.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, states, coupled_states, total_basis, energy_corrs, wfn_corrections, all_energy_corrections=None, degenerate_states=None, degenerate_transformation=None, degenerate_energies=None, logger=None): 
```

- `states`: `BasisStateSpace`
    >No description...
- `coupled_states`: `BasisMultiStateSpace`
    >No description...
- `total_basis`: `BasisMultiStateSpace`
    >No description...
- `energy_corrs`: `np.ndarray`
    >No description...
- `wfn_corrections`: `Iterable[SparseArray]`
    >No description...
- `degenerate_states`: `None | np.ndarray`
    >No description...
- `degenerate_transformation`: `None | np.ndarray`
    >No description...
- `degenerate_energies`: `None | np.ndarray`
    >No description...

<a id="Psience.VPT2.Solver.PerturbationTheoryCorrections.from_dicts" class="docs-object-method">&nbsp;</a>
```python
from_dicts(states, corrections, hamiltonians, **opts): 
```

- `states`: `dict`
    >a dict with the states described by the corrections, the set of states coupled, and the size of the overall basis
- `corrections`: `dict`
    >the corrections generated, including the corrections for the energies, wavefunctions, and a transformation from degenerate PT

<a id="Psience.VPT2.Solver.PerturbationTheoryCorrections.degenerate" class="docs-object-method">&nbsp;</a>
```python
@property
degenerate(self): 
```

<a id="Psience.VPT2.Solver.PerturbationTheoryCorrections.energies" class="docs-object-method">&nbsp;</a>
```python
@property
energies(self): 
```

<a id="Psience.VPT2.Solver.PerturbationTheoryCorrections.order" class="docs-object-method">&nbsp;</a>
```python
@property
order(self): 
```

<a id="Psience.VPT2.Solver.PerturbationTheoryCorrections.take_subspace" class="docs-object-method">&nbsp;</a>
```python
take_subspace(self, space): 
```
Takes only those elements that are in space
- `space`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Solver.PerturbationTheoryCorrections.operator_representation" class="docs-object-method">&nbsp;</a>
```python
operator_representation(self, operator_expansion, order=None, subspace=None, contract=True, logger_symbol='A', logger_conversion=None): 
```
Generates the representation of the operator in the basis of stored states
- `operator_expansion`: `Iterable[float] | Iterable[np.ndarray]`
    >the expansion of the operator
- `order`: `Iterable[float] | Iterable[np.ndarray]`
    >the order of correction to go up to
- `subspace`: `None | BasisStateSpace`
    >the subspace of terms in which the operator expansion is defined
- `:returns`: `Iterable[np.ndarray]`
    >the set of representation matrices for this operator

<a id="Psience.VPT2.Solver.PerturbationTheoryCorrections.get_overlap_matrices" class="docs-object-method">&nbsp;</a>
```python
get_overlap_matrices(self): 
```
Returns the overlap matrices for the set of corrections
        at each order of correction
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Solver.PerturbationTheoryCorrections.savez" class="docs-object-method">&nbsp;</a>
```python
savez(self, file): 
```

<a id="Psience.VPT2.Solver.PerturbationTheoryCorrections.loadz" class="docs-object-method">&nbsp;</a>
```python
loadz(file): 
```

<a id="Psience.VPT2.Solver.PerturbationTheoryCorrections.to_state" class="docs-object-method">&nbsp;</a>
```python
to_state(self, serializer=None): 
```

<a id="Psience.VPT2.Solver.PerturbationTheoryCorrections.from_state" class="docs-object-method">&nbsp;</a>
```python
from_state(data, serializer=None): 
```

### Examples




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/edit/ci/examples/ci/docs/Psience/VPT2/Solver/PerturbationTheoryCorrections.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/examples/ci/docs/Psience/VPT2/Solver/PerturbationTheoryCorrections.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/edit/ci/docs/ci/docs/Psience/VPT2/Solver/PerturbationTheoryCorrections.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/docs/templates/ci/docs/Psience/VPT2/Solver/PerturbationTheoryCorrections.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Solver.py?message=Update%20Docs)