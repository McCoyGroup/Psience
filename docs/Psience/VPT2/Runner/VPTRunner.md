## <a id="Psience.VPT2.Runner.VPTRunner">VPTRunner</a>
A helper class to make it easier to run jobs by making the inputs/options
clear and making it easier to customize run options

### Properties and Methods
<a id="Psience.VPT2.Runner.VPTRunner.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, system, states, hamiltonian_options=None, solver_options=None, runtime_options=None): 
```

- `system`: `VPTSystem`
    >the system to run perturbation theory on
- `hamiltonian_options`: `VPTHamiltonianOptions`
    >options to configure the Hamiltonian
- `solver_options`: `VPTSolverOptions`
    >options to configure the way the perturbation theory is applied
- `runtime_options`: `VPTRuntimeOptions`
    >options to configure the way the code runs

<a id="Psience.VPT2.Runner.VPTRunner.get_Hamiltonian" class="docs-object-method">&nbsp;</a>
```python
get_Hamiltonian(self): 
```

<a id="Psience.VPT2.Runner.VPTRunner.hamiltonian" class="docs-object-method">&nbsp;</a>
```python
@property
hamiltonian(self): 
```

<a id="Psience.VPT2.Runner.VPTRunner.get_wavefunctions" class="docs-object-method">&nbsp;</a>
```python
get_wavefunctions(self): 
```

<a id="Psience.VPT2.Runner.VPTRunner.print_tables" class="docs-object-method">&nbsp;</a>
```python
print_tables(self, wfns=None, file=None, print_intensities=True, sep_char='=', sep_len=100): 
```
Prints a bunch of formatted output data from a PT run
- `wfns`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Runner.VPTRunner.run_simple" class="docs-object-method">&nbsp;</a>
```python
run_simple(system, states, target_property=None, **opts): 
```

### Examples




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/edit/ci/examples/ci/docs/Psience/VPT2/Runner/VPTRunner.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/examples/ci/docs/Psience/VPT2/Runner/VPTRunner.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/edit/ci/docs/ci/docs/Psience/VPT2/Runner/VPTRunner.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/docs/templates/ci/docs/Psience/VPT2/Runner/VPTRunner.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Runner.py?message=Update%20Docs)