## <a id="Psience.VPT2.Analyzer.VPTAnalyzer">VPTAnalyzer</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L418)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L418?message=Update%20Docs)]
</div>

Provides analysis tools on VPT results







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.VPT2.Analyzer.VPTAnalyzer.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, res): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L423)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L423?message=Update%20Docs)]
</div>

  - `res`: `Any`
    >


<a id="Psience.VPT2.Analyzer.VPTAnalyzer.run_VPT" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
run_VPT(cls, *args, logger=None, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L433)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L433?message=Update%20Docs)]
</div>
Runs a VPT calculation through `VPTRunner.run_simple` and
stores the output wave functions to use
  - `args`: `Any`
    > 
  - `kwargs`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTAnalyzer.potential_terms" class="docs-object-method">&nbsp;</a> 
```python
@property
potential_terms(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L457)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L457?message=Update%20Docs)]
</div>
Returns the expansion of the potential
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTAnalyzer.kinetic_terms" class="docs-object-method">&nbsp;</a> 
```python
@property
kinetic_terms(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L465)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L465?message=Update%20Docs)]
</div>
Returns the expansion of the kinetic energy
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTAnalyzer.dipole_terms" class="docs-object-method">&nbsp;</a> 
```python
@property
dipole_terms(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L473)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L473?message=Update%20Docs)]
</div>
Returns the expansion of the dipole moment
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTAnalyzer.basis" class="docs-object-method">&nbsp;</a> 
```python
@property
basis(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L482)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L482?message=Update%20Docs)]
</div>
Returns the basis for the calculation
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTAnalyzer.target_states" class="docs-object-method">&nbsp;</a> 
```python
@property
target_states(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L490)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L490?message=Update%20Docs)]
</div>
Returns the target states for the calculation
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTAnalyzer.spectrum" class="docs-object-method">&nbsp;</a> 
```python
@property
spectrum(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L499)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L499?message=Update%20Docs)]
</div>
Returns the IR spectrum calculated from perturbation theory
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTAnalyzer.energy_corrections" class="docs-object-method">&nbsp;</a> 
```python
@property
energy_corrections(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L507)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L507?message=Update%20Docs)]
</div>
Returns the corrections to the energies
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTAnalyzer.energies" class="docs-object-method">&nbsp;</a> 
```python
@property
energies(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L515)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L515?message=Update%20Docs)]
</div>

  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTAnalyzer.frequencies" class="docs-object-method">&nbsp;</a> 
```python
@property
frequencies(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L526)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L526?message=Update%20Docs)]
</div>

  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTAnalyzer.zero_order_spectrum" class="docs-object-method">&nbsp;</a> 
```python
@property
zero_order_spectrum(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L535)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L535?message=Update%20Docs)]
</div>
Returns the zero-order IR spectrum calculated from perturbation theory
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTAnalyzer.deperturbed_spectrum" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_spectrum(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L543)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L543?message=Update%20Docs)]
</div>

  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTAnalyzer.deperturbed_frequencies" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_frequencies(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L558)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L558?message=Update%20Docs)]
</div>

  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTAnalyzer.wavefunction_corrections" class="docs-object-method">&nbsp;</a> 
```python
@property
wavefunction_corrections(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L566)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L566?message=Update%20Docs)]
</div>
Returns the corrections to the wavefunctions
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTAnalyzer.transition_moment_corrections" class="docs-object-method">&nbsp;</a> 
```python
@property
transition_moment_corrections(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L574)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L574?message=Update%20Docs)]
</div>
Returns the corrections to the wavefunctions
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTAnalyzer.transition_moments" class="docs-object-method">&nbsp;</a> 
```python
@property
transition_moments(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L582)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L582?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analyzer.VPTAnalyzer.deperturbed_transition_moment_corrections" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_transition_moment_corrections(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L590)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L590?message=Update%20Docs)]
</div>
Returns the corrections to the wavefunctions
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTAnalyzer.deperturbed_transition_moments" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_transition_moments(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L598)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L598?message=Update%20Docs)]
</div>

  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTAnalyzer.deperturbed_hamiltonians" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_hamiltonians(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L606)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L606?message=Update%20Docs)]
</div>
Returns the deperturbed Hamiltonians used to make the degenerate transform
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTAnalyzer.deperturbed_energies" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_energies(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L614)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L614?message=Update%20Docs)]
</div>

  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTAnalyzer.degenerate_states" class="docs-object-method">&nbsp;</a> 
```python
@property
degenerate_states(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L622)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L622?message=Update%20Docs)]
</div>
Returns the deperturbed states used to make the degenerate transform
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTAnalyzer.degenerate_energies" class="docs-object-method">&nbsp;</a> 
```python
@property
degenerate_energies(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L630)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L630?message=Update%20Docs)]
</div>
Returns the deperturbed states used to make the degenerate transform
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTAnalyzer.shift_and_transform_hamiltonian" class="docs-object-method">&nbsp;</a> 
```python
shift_and_transform_hamiltonian(self, hams, shifts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L639)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L639?message=Update%20Docs)]
</div>

  - `hams`: `Any`
    > 
  - `shifts`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTAnalyzer.get_shifted_transformed_transition_moments" class="docs-object-method">&nbsp;</a> 
```python
get_shifted_transformed_transition_moments(self, deg_states, target_states, hams, shifts, tmoms, handling_mode='transpose'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L692)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L692?message=Update%20Docs)]
</div>

  - `deg_states`: `Any`
    > 
  - `target_states`: `Any`
    > 
  - `hams`: `Any`
    > 
  - `shifts`: `Any`
    > 
  - `tmoms`: `Any`
    > 
  - `handling_mode`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTAnalyzer.get_shifted_transformed_spectrum" class="docs-object-method">&nbsp;</a> 
```python
get_shifted_transformed_spectrum(self, zpe, deg_states, target_states, hams, shifts, tmoms, handling_mode='transpose'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L722)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L722?message=Update%20Docs)]
</div>

  - `zpe`: `Any`
    > 
  - `deg_states`: `Any`
    > 
  - `target_states`: `Any`
    > 
  - `hams`: `Any`
    > 
  - `shifts`: `Any`
    > 
  - `tmoms`: `Any`
    > 
  - `handling_mode`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTAnalyzer.shifted_transformed_spectrum" class="docs-object-method">&nbsp;</a> 
```python
shifted_transformed_spectrum(self, deg_states, hams, shifts, return_transformation=False, handling_mode='transpose'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L755)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L755?message=Update%20Docs)]
</div>

  - `deg_states`: `Any`
    > 
  - `hams`: `Any`
    > 
  - `shifts`: `Any`
    > 
  - `return_transformation`: `Any`
    > 
  - `handling_mode`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTAnalyzer.transition_data" class="docs-object-method">&nbsp;</a> 
```python
transition_data(self, states, keys=['frequency', 'transition_moment', 'intensity'], data='deperturbed'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L781)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L781?message=Update%20Docs)]
</div>

  - `states`: `Any`
    > 
  - `keys`: `Any`
    > 
  - `data`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTAnalyzer.transition_moment_term_sums" class="docs-object-method">&nbsp;</a> 
```python
transition_moment_term_sums(self, states, terms=None, rotation=None, data='deperturbed'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L846)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L846?message=Update%20Docs)]
</div>

  - `states`: `Any`
    > 
  - `terms`: `Any`
    > 
  - `rotation`: `Any`
    > 
  - `data`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTAnalyzer.transition_moment_term_sums_first_order" class="docs-object-method">&nbsp;</a> 
```python
transition_moment_term_sums_first_order(self, states, rotation=None, data='deperturbed'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L893)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L893?message=Update%20Docs)]
</div>

  - `states`: `Any`
    > 
  - `rotation`: `Any`
    > 
  - `data`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTAnalyzer.intensity_breakdown" class="docs-object-method">&nbsp;</a> 
```python
intensity_breakdown(self, states, terms=None, data='deperturbed'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L916)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L916?message=Update%20Docs)]
</div>

  - `states`: `Any`
    > 
  - `terms`: `Any`
    > 
  - `data`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTAnalyzer.degenerate_coupling_element" class="docs-object-method">&nbsp;</a> 
```python
degenerate_coupling_element(self, state1, state2): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L941)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L941?message=Update%20Docs)]
</div>

  - `state1`: `Any`
    > 
  - `state2`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTAnalyzer.format_deperturbed_hamiltonian" class="docs-object-method">&nbsp;</a> 
```python
format_deperturbed_hamiltonian(self, which): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L964)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L964?message=Update%20Docs)]
</div>

  - `which`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTAnalyzer.log_parser" class="docs-object-method">&nbsp;</a> 
```python
@property
log_parser(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L975)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L975?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analyzer.VPTAnalyzer.print_output_tables" class="docs-object-method">&nbsp;</a> 
```python
print_output_tables(self, print_energy_corrections=False, print_energies=False, print_transition_moments=False, print_intensities=True, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L993)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTAnalyzer.py#L993?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Analyzer/VPTAnalyzer.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Analyzer/VPTAnalyzer.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Analyzer/VPTAnalyzer.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Analyzer/VPTAnalyzer.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L418?message=Update%20Docs)   
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