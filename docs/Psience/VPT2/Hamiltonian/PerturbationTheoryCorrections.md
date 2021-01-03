## <a id="Psience.VPT2.Hamiltonian.PerturbationTheoryCorrections">PerturbationTheoryCorrections</a>
Represents a set of corrections from perturbation theory.
Can be used to correct other operators in the basis of the original calculation.

### Properties and Methods
```python
loadz: method
```
<a id="Psience.VPT2.Hamiltonian.PerturbationTheoryCorrections.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, states, corrections, hamiltonians): 
```

- `states`: `dict`
    >a dict with the states described by the corrections, the set of states coupled, and the size of the overall basis
- `corrections`: `dict`
    >the corrections generated, including the corrections for the energies, wavefunctions, and a transformation from degenerate PT
- `hamiltonians`: `Iterable[np.ndarray]`
    >the set of Hamiltonian matrices used as an expansion

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


