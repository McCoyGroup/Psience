## <a id="Psience.VPT2.Runner.AnalyticVPTRunner">AnalyticVPTRunner</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L2514)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L2514?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
hamiltonian_correction_modification_type: str
```
<a id="Psience.VPT2.Runner.AnalyticVPTRunner.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, expansions, order=None, expansion_order=None, freqs=None, internals=True, logger=None, hamiltonian=None, checkpoint=None, dipole_expansion=None, allowed_terms=None, allowed_coefficients=None, disallowed_coefficients=None, allowed_energy_changes=None, intermediate_normalization=None, local_mode_couplings=None, local_mode_coupling_order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L2515)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L2515?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.from_hamiltonian" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_hamiltonian(cls, ham, order, expansion_order=None, logger=None, checkpoint=None, allowed_terms=None, allowed_coefficients=None, disallowed_coefficients=None, allowed_energy_changes=None, take_diagonal_v4_terms=True, intermediate_normalization=None, corrected_fundamental_frequencies=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L2563)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L2563?message=Update%20Docs)]
</div>
A driver powered by a classic PerturbationTheoryHamiltonian object
  - `ham`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.construct" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
construct(cls, system, states=None, *, order=2, expressions_file=None, allowed_terms=None, allowed_coefficients=None, disallowed_coefficients=None, allowed_energy_changes=None, mixed_derivative_handling_mode='analytical', degeneracy_specs=None, corrected_fundamental_frequencies=None, **settings) -> '(AnalyticVPTRunner, VPTMultiStateSpace)': 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L2645)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L2645?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.from_file" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_file(cls, file_name, order=2, allowed_terms=None, allowed_coefficients=None, disallowed_coefficients=None, allowed_energy_changes=None, expressions_file=None, **settings): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L2706)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L2706?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.construct_classic_runner" class="docs-object-method">&nbsp;</a> 
```python
construct_classic_runner(self, states, system=None, logger=None, corrected_fundamental_frequencies=None, potential_terms=None, kinetic_terms=None, coriolis_terms=None, pseudopotential_terms=None, dipole_terms=None, initial_states=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L2733)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L2733?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.clear_caches" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
clear_caches(cls): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L2826)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L2826?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.prep_multispace" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
prep_multispace(self, states, freqs, system=None, degeneracy_specs=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L2830)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L2830?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.prep_states" class="docs-object-method">&nbsp;</a> 
```python
prep_states(self, states, degeneracy_specs=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L2852)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L2852?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.evaluate_expressions" class="docs-object-method">&nbsp;</a> 
```python
evaluate_expressions(self, states, exprs, zero_cutoff=None, operator_expansions=None, degeneracy_specs=None, verbose=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L2861)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L2861?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.get_matrix_corrections" class="docs-object-method">&nbsp;</a> 
```python
get_matrix_corrections(self, states, order=None, degeneracy_specs=None, zero_cutoff=None, verbose=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L2873)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L2873?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.get_energy_corrections" class="docs-object-method">&nbsp;</a> 
```python
get_energy_corrections(self, states, order=None, degeneracy_specs=None, zero_cutoff=None, verbose=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L2880)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L2880?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.get_overlap_corrections" class="docs-object-method">&nbsp;</a> 
```python
get_overlap_corrections(self, states, order=None, degeneracy_specs=None, zero_cutoff=None, verbose=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L2889)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L2889?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.prep_eval_state_pairs" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
prep_eval_state_pairs(cls, states): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L2901)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L2901?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.get_full_wavefunction_corrections" class="docs-object-method">&nbsp;</a> 
```python
get_full_wavefunction_corrections(self, states, order=None, degeneracy_specs=None, zero_cutoff=None, verbose=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L2908)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L2908?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.get_wavefunction_corrections" class="docs-object-method">&nbsp;</a> 
```python
get_wavefunction_corrections(self, states, order=None, degeneracy_specs=None, zero_cutoff=None, verbose=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L2921)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L2921?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.unflatten_corr" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
unflatten_corr(cls, states, corrs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L2935)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L2935?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.get_operator_corrections" class="docs-object-method">&nbsp;</a> 
```python
get_operator_corrections(self, operator_expansion, states, order=None, terms=None, degeneracy_specs=None, verbose=False, operator_type=None, check_single=True, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L2963)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L2963?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.construct_corrections_vectors" class="docs-object-method">&nbsp;</a> 
```python
construct_corrections_vectors(self, states, corrs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3009)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3009?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.construct_corrections_matrix" class="docs-object-method">&nbsp;</a> 
```python
construct_corrections_matrix(self, group, corrs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3039)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3039?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.get_transition_moment_corrections" class="docs-object-method">&nbsp;</a> 
```python
get_transition_moment_corrections(self, states, dipole_expansion=None, order=None, degeneracy_specs=None, axes=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3066)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3066?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.get_freqs" class="docs-object-method">&nbsp;</a> 
```python
get_freqs(self, states, order=None, degeneracy_specs=None, return_corrections=False, verbose=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3096)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3096?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.get_reexpressed_hamiltonian" class="docs-object-method">&nbsp;</a> 
```python
get_reexpressed_hamiltonian(self, states, order=None, degeneracy_specs=None, only_degenerate_terms=True, verbose=False, hamiltonian_corrections=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3107)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3107?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.get_wfc_test_states" class="docs-object-method">&nbsp;</a> 
```python
get_wfc_test_states(self, input_states: Psience.BasisReps.StateSpaces.BasisStateSpace, energy_window): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3135)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3135?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.get_test_wfn_corrs" class="docs-object-method">&nbsp;</a> 
```python
get_test_wfn_corrs(self, input_states: Psience.BasisReps.StateSpaces.BasisStateSpace, energy_window): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3162)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3162?message=Update%20Docs)]
</div>
We take the expansions and frequencies that we have and at find the possible terms
that could possibly lead to a correction greater than the specified threshold
To do this, we first determine from the expansions what magnitude of energy difference
could possible lead to terms above this threshold


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.format_energies_table" class="docs-object-method">&nbsp;</a> 
```python
format_energies_table(self, states, energies, energy_corrections, zpe_pos, number_format='.3f'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3179)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3179?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.format_degenerate_energies_table" class="docs-object-method">&nbsp;</a> 
```python
format_degenerate_energies_table(self, states, energies, deperturbed_energies, zpe_pos, number_format='.3f'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3212)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3212?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.format_transition_moment_table" class="docs-object-method">&nbsp;</a> 
```python
format_transition_moment_table(self, states, transition_moments, transition_moment_corrections, number_format='.8f'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3240)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3240?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.format_operators_table" class="docs-object-method">&nbsp;</a> 
```python
format_operators_table(self, states, keys, operator_values, operator_corrections, number_format='.8f'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3289)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3289?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.format_spectrum_table" class="docs-object-method">&nbsp;</a> 
```python
format_spectrum_table(self, states, harmonic_spectra, spectra, deperturbed_spectra=None, number_format='.3f'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3349)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3349?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.prep_operators" class="docs-object-method">&nbsp;</a> 
```python
prep_operators(self, operator_expansions, operator_terms, order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3397)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3397?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.format_matrix" class="docs-object-method">&nbsp;</a> 
```python
format_matrix(self, ham): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3438)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3438?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.modify_hamiltonian" class="docs-object-method">&nbsp;</a> 
```python
modify_hamiltonian(self, hamiltonian_corrections): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3443)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3443?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.run_VPT" class="docs-object-method">&nbsp;</a> 
```python
run_VPT(self, states, calculate_intensities=True, operator_expansions=None, operator_terms=None, operator_type=None, order=None, verbose=False, degeneracy_specs=None, handle_degeneracies=True, zero_cutoff=None, transition_moment_terms=None, hamiltonian_corrections=None, clear_caches=True, hamiltonian_correction_type=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3453)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3453?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.run_simple" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
run_simple(cls, system, states, calculate_intensities=True, operator_expansions=None, operator_terms=None, operator_type=None, verbose=False, return_runner=False, degeneracy_specs=None, degeneracy_states=None, handle_degeneracies=True, zero_cutoff=None, clear_caches=True, hamiltonian_correction_type=None, hamiltonian_corrections=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L3684)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L3684?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Runner/AnalyticVPTRunner.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Runner/AnalyticVPTRunner.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Runner/AnalyticVPTRunner.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Runner/AnalyticVPTRunner.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L2514?message=Update%20Docs)   
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