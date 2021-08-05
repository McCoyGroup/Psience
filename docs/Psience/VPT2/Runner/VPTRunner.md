## <a id="Psience.VPT2.Runner.VPTRunner">VPTRunner</a>
A helper class to make it easier to run jobs by making the inputs/options
clear and making it easier to customize run options

### Properties and Methods
```python
InputSystem: type
HamiltonianOptions: type
RuntimeOptions: type
PerturbationTheoryOptions: type
get_states: method
get_degenerate_polyad_space: method
get_state_space_filter: method
```
<a id="Psience.VPT2.Runner.VPTRunner.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, system, hamiltonian_options=None, perturbation_theory_options=None, runtime_options=None): 
```

- `system`: `InputSystem`
    >the system to run perturbation theory on
- `hamiltonian_options`: `HamiltonianOptions`
    >options to configure the Hamiltonian
- `perturbation_theory_options`: `PerturbationTheoryOptions`
    >options to configure the way the perturbation theory is applied
- `runtime_options`: `RuntimeOptions`
    >options to configure the way the code runs

<a id="Psience.VPT2.Runner.VPTRunner.get_Hamiltonian" class="docs-object-method">&nbsp;</a>
```python
get_Hamiltonian(self): 
```

<a id="Psience.VPT2.Runner.VPTRunner.get_wavefunctions" class="docs-object-method">&nbsp;</a>
```python
get_wavefunctions(self): 
```

<a id="Psience.VPT2.Runner.VPTRunner.print_tables" class="docs-object-method">&nbsp;</a>
```python
print_tables(self, wfns=None): 
```
Prints a bunch of formatted output data from a PT run
- `wfns`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

### Examples


