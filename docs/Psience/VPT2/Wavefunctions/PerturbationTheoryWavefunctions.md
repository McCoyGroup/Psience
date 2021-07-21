## <a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions">PerturbationTheoryWavefunctions</a>
These things are fed the first and second order corrections

### Properties and Methods
```python
TermHolder: type
DipolePartitioningMethod: EnumMeta
write_CSV_breakdown: method
```
<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, mol, basis, corrections, modes=None, mode_selection=None, logger=None, operator_settings=None): 
```

- `mol`: `Molecule`
    >the molecule the wavefunction is for
- `basis`: `SimpleProductBasis`
    >the basis the expansion is being done in
- `corrections`: `PerturbationTheoryCorrections`
    >the corrections to the terms

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.energies_to_order" class="docs-object-method">&nbsp;</a>
```python
energies_to_order(self, order): 
```

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.order" class="docs-object-method">&nbsp;</a>
```python
@property
order(self): 
```

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.expectation" class="docs-object-method">&nbsp;</a>
```python
expectation(self, operator, other): 
```

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.zero_order_energies" class="docs-object-method">&nbsp;</a>
```python
@property
zero_order_energies(self): 
```

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.get_M0" class="docs-object-method">&nbsp;</a>
```python
get_M0(self, mu_0): 
```

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.get_M1" class="docs-object-method">&nbsp;</a>
```python
get_M1(self, mu_1): 
```

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.get_M2" class="docs-object-method">&nbsp;</a>
```python
get_M2(self, mu_2): 
```

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.get_M3" class="docs-object-method">&nbsp;</a>
```python
get_M3(self, mu_3): 
```

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.dipole_terms" class="docs-object-method">&nbsp;</a>
```python
@property
dipole_terms(self): 
```

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.dipole_partitioning" class="docs-object-method">&nbsp;</a>
```python
@property
dipole_partitioning(self): 
```

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.transition_moments" class="docs-object-method">&nbsp;</a>
```python
@property
transition_moments(self): 
```
Computes the transition moments between wavefunctions stored in the object
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.transition_moment_corrections" class="docs-object-method">&nbsp;</a>
```python
@property
transition_moment_corrections(self): 
```
Computes the transition moment corrections between wavefunctions stored in the object
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.oscillator_strengths" class="docs-object-method">&nbsp;</a>
```python
@property
oscillator_strengths(self): 
```
Computes the oscillator strengths for transitions from the ground state to the other states
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.oscillator_strengths_to_order" class="docs-object-method">&nbsp;</a>
```python
oscillator_strengths_to_order(self, order): 
```

- `tms`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.intensities" class="docs-object-method">&nbsp;</a>
```python
@property
intensities(self): 
```
Computes the intensities for transitions from the ground state to the other states
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.intensities_to_order" class="docs-object-method">&nbsp;</a>
```python
intensities_to_order(self, order): 
```
Computes the intensities for transitions from the ground state to the other states
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.zero_order_intensities" class="docs-object-method">&nbsp;</a>
```python
@property
zero_order_intensities(self): 
```
Computes the harmonic intensities for transitions from the ground state to the other states
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.generate_intensity_breakdown" class="docs-object-method">&nbsp;</a>
```python
generate_intensity_breakdown(self, include_wavefunctions=True): 
```
Generates a breakdown of the terms that contribute to the intensity
        Returns in a format that can be directly exported to JSON if desired.
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.format_energies_table" class="docs-object-method">&nbsp;</a>
```python
format_energies_table(self, states=None, zpe=None, freqs=None, real_fmt='{:>12.5f}', dash_fmt='{:>12}'): 
```

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.format_intensities_table" class="docs-object-method">&nbsp;</a>
```python
format_intensities_table(self): 
```

### Examples


