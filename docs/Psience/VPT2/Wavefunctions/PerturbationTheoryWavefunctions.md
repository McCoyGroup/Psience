## <a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions">PerturbationTheoryWavefunctions</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L55)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L55?message=Update%20Docs)]
</div>

These things are fed the first and second order corrections

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

```python
TermHolder: type
DipolePartitioningMethod: EnumMeta
```
<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, mol, basis, corrections, modes=None, mode_selection=None, logger=None, checkpoint=None, results=None, operator_settings=None, expansion_options=None, degenerate_transformation_layout=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L60)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L60?message=Update%20Docs)]
</div>


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
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.to_state" class="docs-object-method">&nbsp;</a> 
```python
to_state(self, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L110)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L110?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.from_state" class="docs-object-method">&nbsp;</a> 
```python
from_state(data, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L123)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L123?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.degenerate_transformation" class="docs-object-method">&nbsp;</a> 
```python
@property
degenerate_transformation(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.energies_to_order" class="docs-object-method">&nbsp;</a> 
```python
energies_to_order(self, order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L140)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L140?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.deperturbed_energies" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_energies(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.deperturbed_frequencies" class="docs-object-method">&nbsp;</a> 
```python
deperturbed_frequencies(self, order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L147)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L147?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.order" class="docs-object-method">&nbsp;</a> 
```python
@property
order(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.expectation" class="docs-object-method">&nbsp;</a> 
```python
expectation(self, operator, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L158)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L158?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.zero_order_energies" class="docs-object-method">&nbsp;</a> 
```python
@property
zero_order_energies(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.get_M0" class="docs-object-method">&nbsp;</a> 
```python
get_M0(self, mu_0): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L165)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L165?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.get_M1" class="docs-object-method">&nbsp;</a> 
```python
get_M1(self, mu_1): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L170)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L170?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.get_M2" class="docs-object-method">&nbsp;</a> 
```python
get_M2(self, mu_2): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L175)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L175?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.get_M3" class="docs-object-method">&nbsp;</a> 
```python
get_M3(self, mu_3): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L181)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L181?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.get_Mi" class="docs-object-method">&nbsp;</a> 
```python
get_Mi(self, i, mu): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L187)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L187?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.dipole_terms" class="docs-object-method">&nbsp;</a> 
```python
@property
dipole_terms(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.dipole_partitioning" class="docs-object-method">&nbsp;</a> 
```python
@property
dipole_partitioning(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.transition_moments" class="docs-object-method">&nbsp;</a> 
```python
@property
transition_moments(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L?message=Update%20Docs)]
</div>

Computes the transition moments between wavefunctions stored in the object
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.deperturbed_transition_moments" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_transition_moments(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L?message=Update%20Docs)]
</div>

Computes the transition moments between wavefunctions stored in the object
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.transition_moment_corrections" class="docs-object-method">&nbsp;</a> 
```python
@property
transition_moment_corrections(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L?message=Update%20Docs)]
</div>

Computes the transition moment corrections between wavefunctions stored in the object
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.deperturbed_transition_moment_corrections" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_transition_moment_corrections(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L?message=Update%20Docs)]
</div>

Computes the transition moment corrections between wavefunctions stored in the object
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.oscillator_strengths" class="docs-object-method">&nbsp;</a> 
```python
@property
oscillator_strengths(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L?message=Update%20Docs)]
</div>

Computes the oscillator strengths for transitions from the ground state to the other states
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.deperturbed_oscillator_strengths" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_oscillator_strengths(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L?message=Update%20Docs)]
</div>

Computes the oscillator strengths for transitions from the ground state to the other states
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.oscillator_strengths_to_order" class="docs-object-method">&nbsp;</a> 
```python
oscillator_strengths_to_order(self, order): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L820)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L820?message=Update%20Docs)]
</div>


- `tms`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.deperturbed_oscillator_strengths_to_order" class="docs-object-method">&nbsp;</a> 
```python
deperturbed_oscillator_strengths_to_order(self, order): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L834)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L834?message=Update%20Docs)]
</div>


- `tms`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.intensities" class="docs-object-method">&nbsp;</a> 
```python
@property
intensities(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L?message=Update%20Docs)]
</div>

Computes the intensities for transitions from the ground state to the other states
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.deperturbed_intensities" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_intensities(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L?message=Update%20Docs)]
</div>

Computes the intensities for transitions from the ground state to the other states
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.intensities_to_order" class="docs-object-method">&nbsp;</a> 
```python
intensities_to_order(self, order): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L882)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L882?message=Update%20Docs)]
</div>

Computes the intensities for transitions from the ground state to the other states
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.deperturbed_intensities_to_order" class="docs-object-method">&nbsp;</a> 
```python
deperturbed_intensities_to_order(self, order): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L891)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L891?message=Update%20Docs)]
</div>

Computes the intensities for transitions from the ground state to the other states
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.zero_order_intensities" class="docs-object-method">&nbsp;</a> 
```python
@property
zero_order_intensities(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L?message=Update%20Docs)]
</div>

Computes the harmonic intensities for transitions from the ground state to the other states
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.generate_intensity_breakdown" class="docs-object-method">&nbsp;</a> 
```python
generate_intensity_breakdown(self, include_wavefunctions=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L939)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L939?message=Update%20Docs)]
</div>

Generates a breakdown of the terms that contribute to the intensity
        Returns in a format that can be directly exported to JSON if desired.
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.write_CSV_breakdown" class="docs-object-method">&nbsp;</a> 
```python
write_CSV_breakdown(file, intensity_breakdown, padding=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L1024)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L1024?message=Update%20Docs)]
</div>

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
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L1150)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L1150?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.format_deperturbed_energies_table" class="docs-object-method">&nbsp;</a> 
```python
format_deperturbed_energies_table(self, states=None, zpe=None, freqs=None, real_fmt='{:>12.5f}'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L1175)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L1175?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.format_property_matrices" class="docs-object-method">&nbsp;</a> 
```python
format_property_matrices(self, states, prop_corrs, real_fmt='{:>.8e}', padding_fmt='{:>16}'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L1200)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L1200?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.format_dipole_contribs_tables" class="docs-object-method">&nbsp;</a> 
```python
format_dipole_contribs_tables(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L1239)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L1239?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.format_deperturbed_dipole_contribs_tables" class="docs-object-method">&nbsp;</a> 
```python
format_deperturbed_dipole_contribs_tables(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L1250)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L1250?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.format_energy_corrections_table" class="docs-object-method">&nbsp;</a> 
```python
format_energy_corrections_table(self, real_fmt='{:>12.5f}'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L1261)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L1261?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.format_intensities_table" class="docs-object-method">&nbsp;</a> 
```python
format_intensities_table(self, real_fmt='{:>12.5f}'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L1324)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L1324?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.format_deperturbed_intensities_table" class="docs-object-method">&nbsp;</a> 
```python
format_deperturbed_intensities_table(self, real_fmt='{:>12.5f}'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Wavefunctions.py#L1339)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L1339?message=Update%20Docs)]
</div>

 </div>
</div>




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Wavefunctions.py#L55?message=Update%20Docs)