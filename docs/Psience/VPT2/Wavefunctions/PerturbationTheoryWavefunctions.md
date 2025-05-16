## <a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions">PerturbationTheoryWavefunctions</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions.py#L56)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions.py#L56?message=Update%20Docs)]
</div>

These things are fed the first and second order corrections







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
TermHolder: TermHolder
DipolePartitioningMethod: DipolePartitioningMethod
```
<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, mol, basis, corrections, initial_states=None, modes=None, mode_selection=None, logger=None, checkpoint=None, results=None, operator_settings=None, expansion_options=None, degenerate_transformation_layout=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L61)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L61?message=Update%20Docs)]
</div>

  - `mol`: `Molecule`
    > the molecule the wavefunction is for
  - `basis`: `SimpleProductBasis`
    > the basis the expansion is being done in
  - `corrections`: `PerturbationTheoryCorrections`
    > the corrections to the terms


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.get_dimension" class="docs-object-method">&nbsp;</a> 
```python
get_dimension(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L118)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L118?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.energies" class="docs-object-method">&nbsp;</a> 
```python
@property
energies(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L121)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L121?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.to_state" class="docs-object-method">&nbsp;</a> 
```python
to_state(self, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L129)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L129?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.from_state" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_state(cls, data, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L142)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L142?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.degenerate_transformation" class="docs-object-method">&nbsp;</a> 
```python
@property
degenerate_transformation(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L154)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L154?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.initial_state_indices" class="docs-object-method">&nbsp;</a> 
```python
@property
initial_state_indices(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L158)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L158?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.energies_to_order" class="docs-object-method">&nbsp;</a> 
```python
energies_to_order(self, order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L167)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L167?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.deperturbed_energies" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_energies(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L171)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L171?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.deperturbed_frequencies" class="docs-object-method">&nbsp;</a> 
```python
deperturbed_frequencies(self, order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L174)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L174?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.order" class="docs-object-method">&nbsp;</a> 
```python
@property
order(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L180)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L180?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.expectation" class="docs-object-method">&nbsp;</a> 
```python
expectation(self, operator, other=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L187)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L187?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.zero_order_energies" class="docs-object-method">&nbsp;</a> 
```python
@property
zero_order_energies(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L190)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L190?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.get_M0" class="docs-object-method">&nbsp;</a> 
```python
get_M0(self, mu_0): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L194)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L194?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.get_M1" class="docs-object-method">&nbsp;</a> 
```python
get_M1(self, mu_1): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L199)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L199?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.get_M2" class="docs-object-method">&nbsp;</a> 
```python
get_M2(self, mu_2): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L204)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L204?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.get_M3" class="docs-object-method">&nbsp;</a> 
```python
get_M3(self, mu_3): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L210)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L210?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.get_Mi" class="docs-object-method">&nbsp;</a> 
```python
get_Mi(self, i, mu, base_sym='M'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L215)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L215?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.dipole_terms" class="docs-object-method">&nbsp;</a> 
```python
@property
dipole_terms(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1044)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1044?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.dipole_partitioning" class="docs-object-method">&nbsp;</a> 
```python
@property
dipole_partitioning(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1074)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1074?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.transition_moments" class="docs-object-method">&nbsp;</a> 
```python
@property
transition_moments(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1092)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1092?message=Update%20Docs)]
</div>
Computes the transition moments between wavefunctions stored in the object
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.deperturbed_transition_moments" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_transition_moments(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1107)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1107?message=Update%20Docs)]
</div>
Computes the transition moments between wavefunctions stored in the object
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.transition_moment_corrections" class="docs-object-method">&nbsp;</a> 
```python
@property
transition_moment_corrections(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1119)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1119?message=Update%20Docs)]
</div>
Computes the transition moment corrections between wavefunctions stored in the object
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.deperturbed_transition_moment_corrections" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_transition_moment_corrections(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1134)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1134?message=Update%20Docs)]
</div>
Computes the transition moment corrections between wavefunctions stored in the object
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.oscillator_strengths" class="docs-object-method">&nbsp;</a> 
```python
@property
oscillator_strengths(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1146)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1146?message=Update%20Docs)]
</div>
Computes the oscillator strengths for transitions from the ground state to the other states
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.deperturbed_oscillator_strengths" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_oscillator_strengths(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1158)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1158?message=Update%20Docs)]
</div>
Computes the oscillator strengths for transitions from the ground state to the other states
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.oscillator_strengths_to_order" class="docs-object-method">&nbsp;</a> 
```python
oscillator_strengths_to_order(self, order): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1170)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1170?message=Update%20Docs)]
</div>

  - `tms`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.deperturbed_oscillator_strengths_to_order" class="docs-object-method">&nbsp;</a> 
```python
deperturbed_oscillator_strengths_to_order(self, order): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1184)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1184?message=Update%20Docs)]
</div>

  - `tms`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.intensities" class="docs-object-method">&nbsp;</a> 
```python
@property
intensities(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1209)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1209?message=Update%20Docs)]
</div>
Computes the intensities for transitions from the ground state to the other states
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.deperturbed_intensities" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_intensities(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1224)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1224?message=Update%20Docs)]
</div>
Computes the intensities for transitions from the ground state to the other states
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.intensities_to_order" class="docs-object-method">&nbsp;</a> 
```python
intensities_to_order(self, order, return_freqs=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1234)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1234?message=Update%20Docs)]
</div>
Computes the intensities for transitions from the ground state to the other states
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.deperturbed_intensities_to_order" class="docs-object-method">&nbsp;</a> 
```python
deperturbed_intensities_to_order(self, order, return_freqs=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1246)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1246?message=Update%20Docs)]
</div>
Computes the intensities for transitions from the ground state to the other states
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.zero_order_intensities" class="docs-object-method">&nbsp;</a> 
```python
@property
zero_order_intensities(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1296)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1296?message=Update%20Docs)]
</div>
Computes the harmonic intensities for transitions from the ground state to the other states
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.prep_operator_terms" class="docs-object-method">&nbsp;</a> 
```python
prep_operator_terms(self, coeffs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1311)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1311?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.operator_correction_data" class="docs-object-method">&nbsp;</a> 
```python
operator_correction_data(self, operator_coeffs, order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1338)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1338?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.generate_intensity_breakdown" class="docs-object-method">&nbsp;</a> 
```python
generate_intensity_breakdown(self, include_wavefunctions=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1358)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1358?message=Update%20Docs)]
</div>
Generates a breakdown of the terms that contribute to the intensity
Returns in a format that can be directly exported to JSON if desired.
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.write_CSV_breakdown" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
write_CSV_breakdown(cls, file, intensity_breakdown, padding=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L1443)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L1443?message=Update%20Docs)]
</div>
Writes an intensity breakdown to a CSV by annoyingly flattening all the arrays
  - `file`: `Any`
    > 
  - `intensity_breakdown`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.format_energies_table" class="docs-object-method">&nbsp;</a> 
```python
format_energies_table(self, states=None, zpe=None, freqs=None, real_fmt='{:>12.5f}'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1620)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1620?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.format_deperturbed_energies_table" class="docs-object-method">&nbsp;</a> 
```python
format_deperturbed_energies_table(self, states=None, zpe=None, freqs=None, real_fmt='{:>12.5f}'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1629)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1629?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.format_property_matrices" class="docs-object-method">&nbsp;</a> 
```python
format_property_matrices(self, states, prop_corrs, real_fmt='{:>.8e}', padding_fmt='{:>16}'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1638)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1638?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.format_dipole_contribs_tables" class="docs-object-method">&nbsp;</a> 
```python
format_dipole_contribs_tables(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1694)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1694?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.format_deperturbed_dipole_contribs_tables" class="docs-object-method">&nbsp;</a> 
```python
format_deperturbed_dipole_contribs_tables(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1705)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1705?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.format_energy_corrections_table" class="docs-object-method">&nbsp;</a> 
```python
format_energy_corrections_table(self, real_fmt='{:>12.5f}'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1716)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1716?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.format_intensities_table" class="docs-object-method">&nbsp;</a> 
```python
format_intensities_table(self, real_fmt='{:>12.5f}'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1832)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1832?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.format_deperturbed_intensities_table" class="docs-object-method">&nbsp;</a> 
```python
format_deperturbed_intensities_table(self, real_fmt='{:>12.5f}'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1840)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1840?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.get_spectrum" class="docs-object-method">&nbsp;</a> 
```python
get_spectrum(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1880)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1880?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.get_deperturbed_spectrum" class="docs-object-method">&nbsp;</a> 
```python
get_deperturbed_spectrum(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1885)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1885?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunctions.format_operator_table" class="docs-object-method">&nbsp;</a> 
```python
format_operator_table(self, operators, real_fmt='{:>12.5f}'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1961)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.py#L1961?message=Update%20Docs)]
</div>
 </div>
</div>












---


<div markdown="1" class="text-secondary">
<div class="container">
  <div class="row">
   <div class="col" markdown="1">
**Feedback**   
</div>
   <div class="col" markdown="1">
**Examples**   
</div>
   <div class="col" markdown="1">
**Templates**   
</div>
   <div class="col" markdown="1">
**Documentation**   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
</div>
  <div class="row">
   <div class="col" markdown="1">
[Bug](https://github.com/McCoyGroup/Psience/issues/new?title=Documentation%20Improvement%20Needed)/[Request](https://github.com/McCoyGroup/Psience/issues/new?title=Example%20Request)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Wavefunctions/PerturbationTheoryWavefunctions.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Wavefunctions.py#L56?message=Update%20Docs)   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
</div>
</div>
</div>