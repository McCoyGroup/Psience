## <a id="Psience.VPT2.Analyzer.VPTAnalyzer">VPTAnalyzer</a>
Provides analysis tools on VPT results

### Properties and Methods
<a id="Psience.VPT2.Analyzer.VPTAnalyzer.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, res): 
```

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.potential_terms" class="docs-object-method">&nbsp;</a>
```python
@property
potential_terms(self): 
```
Returns the expansion of the potential
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.kinetic_terms" class="docs-object-method">&nbsp;</a>
```python
@property
kinetic_terms(self): 
```
Returns the expansion of the kinetic energy
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.dipole_terms" class="docs-object-method">&nbsp;</a>
```python
@property
dipole_terms(self): 
```
Returns the expansion of the dipole moment
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.basis" class="docs-object-method">&nbsp;</a>
```python
@property
basis(self): 
```
Returns the basis for the calculation
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.target_states" class="docs-object-method">&nbsp;</a>
```python
@property
target_states(self): 
```
Returns the target states for the calculation
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.spectrum" class="docs-object-method">&nbsp;</a>
```python
@property
spectrum(self): 
```
Returns the IR spectrum calculated from perturbation theory
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.energy_corrections" class="docs-object-method">&nbsp;</a>
```python
@property
energy_corrections(self): 
```
Returns the corrections to the energies
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.energies" class="docs-object-method">&nbsp;</a>
```python
@property
energies(self): 
```

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.frequencies" class="docs-object-method">&nbsp;</a>
```python
@property
frequencies(self): 
```

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.deperturbed_spectrum" class="docs-object-method">&nbsp;</a>
```python
@property
deperturbed_spectrum(self): 
```

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.deperturbed_frequencies" class="docs-object-method">&nbsp;</a>
```python
@property
deperturbed_frequencies(self): 
```

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.wavefunction_corrections" class="docs-object-method">&nbsp;</a>
```python
@property
wavefunction_corrections(self): 
```
Returns the corrections to the wavefunctions
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.transition_moment_corrections" class="docs-object-method">&nbsp;</a>
```python
@property
transition_moment_corrections(self): 
```
Returns the corrections to the wavefunctions
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.transition_moments" class="docs-object-method">&nbsp;</a>
```python
@property
transition_moments(self): 
```

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.deperturbed_transition_moment_corrections" class="docs-object-method">&nbsp;</a>
```python
@property
deperturbed_transition_moment_corrections(self): 
```
Returns the corrections to the wavefunctions
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.deperturbed_transition_moments" class="docs-object-method">&nbsp;</a>
```python
@property
deperturbed_transition_moments(self): 
```

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.deperturbed_hamiltonians" class="docs-object-method">&nbsp;</a>
```python
@property
deperturbed_hamiltonians(self): 
```
Returns the deperturbed Hamiltonians used to make the degenerate transform
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.deperturbed_energies" class="docs-object-method">&nbsp;</a>
```python
@property
deperturbed_energies(self): 
```

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.degenerate_states" class="docs-object-method">&nbsp;</a>
```python
@property
degenerate_states(self): 
```
Returns the deperturbed states used to make the degenerate transform
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.degenerate_energies" class="docs-object-method">&nbsp;</a>
```python
@property
degenerate_energies(self): 
```
Returns the deperturbed states used to make the degenerate transform
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.shift_and_transform_hamiltonian" class="docs-object-method">&nbsp;</a>
```python
shift_and_transform_hamiltonian(self, hams, shifts): 
```

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.get_shifted_transformed_transition_moments" class="docs-object-method">&nbsp;</a>
```python
get_shifted_transformed_transition_moments(self, deg_states, target_states, hams, shifts, tmoms, handling_mode='transpose'): 
```

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.get_shifted_transformed_spectrum" class="docs-object-method">&nbsp;</a>
```python
get_shifted_transformed_spectrum(self, zpe, deg_states, target_states, hams, shifts, tmoms, handling_mode='transpose'): 
```

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.shifted_transformed_spectrum" class="docs-object-method">&nbsp;</a>
```python
shifted_transformed_spectrum(self, deg_states, hams, shifts, return_transformation=False, handling_mode='transpose'): 
```

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.transition_data" class="docs-object-method">&nbsp;</a>
```python
transition_data(self, states, keys=['frequency', 'transition_moment', 'intensity'], data='deperturbed'): 
```

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.transition_moment_term_sums" class="docs-object-method">&nbsp;</a>
```python
transition_moment_term_sums(self, states, terms=None, data='deperturbed'): 
```

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.transition_moment_term_sums_first_order" class="docs-object-method">&nbsp;</a>
```python
transition_moment_term_sums_first_order(self, states, data='deperturbed'): 
```

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.intensity_breakdown" class="docs-object-method">&nbsp;</a>
```python
intensity_breakdown(self, states, terms=None, data='deperturbed'): 
```

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.degenerate_coupling_element" class="docs-object-method">&nbsp;</a>
```python
degenerate_coupling_element(self, state1, state2): 
```

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.format_deperturbed_hamiltonian" class="docs-object-method">&nbsp;</a>
```python
format_deperturbed_hamiltonian(self, which): 
```

### Examples




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/edit/ci/examples/ci/docs/Psience/VPT2/Analyzer/VPTAnalyzer.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/examples/ci/docs/Psience/VPT2/Analyzer/VPTAnalyzer.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/edit/ci/docs/ci/docs/Psience/VPT2/Analyzer/VPTAnalyzer.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/docs/templates/ci/docs/Psience/VPT2/Analyzer/VPTAnalyzer.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Analyzer.py?message=Update%20Docs)