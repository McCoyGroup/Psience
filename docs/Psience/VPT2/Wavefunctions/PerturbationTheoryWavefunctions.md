## <a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions">PerturbationTheoryWavefunctions</a>
These things are fed the first and second order corrections

### Properties and Methods
```python
TermHolder: type
DipolePartitioningMethod: EnumMeta
```
<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, mol, basis, corrections, modes=None, mode_selection=None, logger=None, checkpoint=None, results=None, operator_settings=None, expansion_options=None, degenerate_transformation_layout=None): 
```

- `mol`: `Molecule`
    >the molecule the wavefunction is for
- `basis`: `SimpleProductBasis`
    >the basis the expansion is being done in
- `corrections`: `PerturbationTheoryCorrections`
    >the corrections to the terms

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.energies" class="docs-object-method">&nbsp;</a>
```python
@property
energies(self): 
```

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.to_state" class="docs-object-method">&nbsp;</a>
```python
to_state(self, serializer=None): 
```

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.from_state" class="docs-object-method">&nbsp;</a>
```python
from_state(data, serializer=None): 
```

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.degenerate_transformation" class="docs-object-method">&nbsp;</a>
```python
@property
degenerate_transformation(self): 
```

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.deperturbed_energies" class="docs-object-method">&nbsp;</a>
```python
@property
deperturbed_energies(self): 
```

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

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.get_Mi" class="docs-object-method">&nbsp;</a>
```python
get_Mi(self, i, mu): 
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

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.deperturbed_transition_moments" class="docs-object-method">&nbsp;</a>
```python
@property
deperturbed_transition_moments(self): 
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

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.deperturbed_transition_moment_corrections" class="docs-object-method">&nbsp;</a>
```python
@property
deperturbed_transition_moment_corrections(self): 
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

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.deperturbed_oscillator_strengths" class="docs-object-method">&nbsp;</a>
```python
@property
deperturbed_oscillator_strengths(self): 
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

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.deperturbed_oscillator_strengths_to_order" class="docs-object-method">&nbsp;</a>
```python
deperturbed_oscillator_strengths_to_order(self, order): 
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

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.deperturbed_intensities" class="docs-object-method">&nbsp;</a>
```python
@property
deperturbed_intensities(self): 
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

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.deperturbed_intensities_to_order" class="docs-object-method">&nbsp;</a>
```python
deperturbed_intensities_to_order(self, order): 
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

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.write_CSV_breakdown" class="docs-object-method">&nbsp;</a>
```python
write_CSV_breakdown(file, intensity_breakdown, padding=None): 
```
Writes an intensity breakdown to a CSV by annoyingly flattening all the arrays
- `file`: `Any`
    >No description...
- `intensity_breakdown`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.format_energies_table" class="docs-object-method">&nbsp;</a>
```python
format_energies_table(self, states=None, zpe=None, freqs=None, real_fmt='{:>12.5f}'): 
```

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.format_deperturbed_energies_table" class="docs-object-method">&nbsp;</a>
```python
format_deperturbed_energies_table(self, states=None, zpe=None, freqs=None, real_fmt='{:>12.5f}'): 
```

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.format_property_matrices" class="docs-object-method">&nbsp;</a>
```python
format_property_matrices(self, states, prop_corrs, real_fmt='{:>.8e}', padding_fmt='{:>16}'): 
```

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.format_dipole_contribs_tables" class="docs-object-method">&nbsp;</a>
```python
format_dipole_contribs_tables(self): 
```

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.format_deperturbed_dipole_contribs_tables" class="docs-object-method">&nbsp;</a>
```python
format_deperturbed_dipole_contribs_tables(self): 
```

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.format_energy_corrections_table" class="docs-object-method">&nbsp;</a>
```python
format_energy_corrections_table(self, real_fmt='{:>12.5f}'): 
```

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.format_intensities_table" class="docs-object-method">&nbsp;</a>
```python
format_intensities_table(self, real_fmt='{:>12.5f}'): 
```

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.format_deperturbed_intensities_table" class="docs-object-method">&nbsp;</a>
```python
format_deperturbed_intensities_table(self, real_fmt='{:>12.5f}'): 
```

### Examples




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/edit/ci/examples/ci/docs/Psience/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/examples/ci/docs/Psience/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/edit/ci/docs/ci/docs/Psience/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/docs/templates/ci/docs/Psience/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Wavefunctions.py?message=Update%20Docs)