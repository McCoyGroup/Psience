## <a id="Psience.VPT2.Hamiltonian.PerturbationTheoryCorrections">PerturbationTheoryCorrections</a>
Represents a set of corrections from perturbation theory.
Can be used to correct other operators in the basis of the original calculation.

### Properties and Methods
```python
from_dicts: method
loadz: method
```
<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryCorrections.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, hamiltonians, states, coupled_states, total_basis, energy_corrs, wfn_corrections, degenerate_states=None, degenerate_transformation=None, degenerate_energies=None): 
```

- `hamiltonians`: `Iterable[SparseArray]`
    >No description...
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

<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryCorrections.degenerate" class="docs-object-method">&nbsp;</a>
```python
@property
degenerate(self): 
```

<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryCorrections.energies" class="docs-object-method">&nbsp;</a>
```python
@property
energies(self): 
```

<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryCorrections.order" class="docs-object-method">&nbsp;</a>
```python
@property
order(self): 
```

<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryCorrections.take_subspace" class="docs-object-method">&nbsp;</a>
```python
take_subspace(self, space): 
```
Takes only those elements that are in space
- `space`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryCorrections.operator_representation" class="docs-object-method">&nbsp;</a>
```python
operator_representation(self, operator_expansion, order=None, subspace=None): 
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

<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryCorrections.savez" class="docs-object-method">&nbsp;</a>
```python
savez(self, file): 
```

### Examples


