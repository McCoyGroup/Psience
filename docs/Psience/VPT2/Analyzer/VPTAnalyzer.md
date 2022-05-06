## <a id="Psience.VPT2.Analyzer.VPTAnalyzer">VPTAnalyzer</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L411)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L411?message=Update%20Docs)]
</div>

Provides analysis tools on VPT results

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, res): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L416)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L416?message=Update%20Docs)]
</div>


- `res`: `Any`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.run_VPT" class="docs-object-method">&nbsp;</a> 
```python
run_VPT(*args, logger=None, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L426)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L426?message=Update%20Docs)]
</div>

Runs a VPT calculation through `VPTRunner.run_simple` and
        stores the output wave functions to use
- `args`: `Any`
    >No description...
- `kwargs`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.potential_terms" class="docs-object-method">&nbsp;</a> 
```python
@property
potential_terms(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L?message=Update%20Docs)]
</div>

Returns the expansion of the potential
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.kinetic_terms" class="docs-object-method">&nbsp;</a> 
```python
@property
kinetic_terms(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L?message=Update%20Docs)]
</div>

Returns the expansion of the kinetic energy
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.dipole_terms" class="docs-object-method">&nbsp;</a> 
```python
@property
dipole_terms(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L?message=Update%20Docs)]
</div>

Returns the expansion of the dipole moment
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.basis" class="docs-object-method">&nbsp;</a> 
```python
@property
basis(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L?message=Update%20Docs)]
</div>

Returns the basis for the calculation
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.target_states" class="docs-object-method">&nbsp;</a> 
```python
@property
target_states(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L?message=Update%20Docs)]
</div>

Returns the target states for the calculation
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.spectrum" class="docs-object-method">&nbsp;</a> 
```python
@property
spectrum(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L?message=Update%20Docs)]
</div>

Returns the IR spectrum calculated from perturbation theory
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.energy_corrections" class="docs-object-method">&nbsp;</a> 
```python
@property
energy_corrections(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L?message=Update%20Docs)]
</div>

Returns the corrections to the energies
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.energies" class="docs-object-method">&nbsp;</a> 
```python
@property
energies(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L?message=Update%20Docs)]
</div>


- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.frequencies" class="docs-object-method">&nbsp;</a> 
```python
@property
frequencies(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L?message=Update%20Docs)]
</div>


- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.zero_order_spectrum" class="docs-object-method">&nbsp;</a> 
```python
@property
zero_order_spectrum(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L?message=Update%20Docs)]
</div>

Returns the zero-order IR spectrum calculated from perturbation theory
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.deperturbed_spectrum" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_spectrum(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L?message=Update%20Docs)]
</div>


- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.deperturbed_frequencies" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_frequencies(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L?message=Update%20Docs)]
</div>


- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.wavefunction_corrections" class="docs-object-method">&nbsp;</a> 
```python
@property
wavefunction_corrections(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L?message=Update%20Docs)]
</div>

Returns the corrections to the wavefunctions
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.transition_moment_corrections" class="docs-object-method">&nbsp;</a> 
```python
@property
transition_moment_corrections(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L?message=Update%20Docs)]
</div>

Returns the corrections to the wavefunctions
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.transition_moments" class="docs-object-method">&nbsp;</a> 
```python
@property
transition_moments(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.deperturbed_transition_moment_corrections" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_transition_moment_corrections(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L?message=Update%20Docs)]
</div>

Returns the corrections to the wavefunctions
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.deperturbed_transition_moments" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_transition_moments(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L?message=Update%20Docs)]
</div>


- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.deperturbed_hamiltonians" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_hamiltonians(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L?message=Update%20Docs)]
</div>

Returns the deperturbed Hamiltonians used to make the degenerate transform
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.deperturbed_energies" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_energies(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L?message=Update%20Docs)]
</div>


- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.degenerate_states" class="docs-object-method">&nbsp;</a> 
```python
@property
degenerate_states(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L?message=Update%20Docs)]
</div>

Returns the deperturbed states used to make the degenerate transform
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.degenerate_energies" class="docs-object-method">&nbsp;</a> 
```python
@property
degenerate_energies(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L?message=Update%20Docs)]
</div>

Returns the deperturbed states used to make the degenerate transform
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.shift_and_transform_hamiltonian" class="docs-object-method">&nbsp;</a> 
```python
shift_and_transform_hamiltonian(self, hams, shifts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L630)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L630?message=Update%20Docs)]
</div>


- `hams`: `Any`
    >No description...
- `shifts`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.get_shifted_transformed_transition_moments" class="docs-object-method">&nbsp;</a> 
```python
get_shifted_transformed_transition_moments(self, deg_states, target_states, hams, shifts, tmoms, handling_mode='transpose'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L679)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L679?message=Update%20Docs)]
</div>


- `deg_states`: `Any`
    >No description...
- `target_states`: `Any`
    >No description...
- `hams`: `Any`
    >No description...
- `shifts`: `Any`
    >No description...
- `tmoms`: `Any`
    >No description...
- `handling_mode`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.get_shifted_transformed_spectrum" class="docs-object-method">&nbsp;</a> 
```python
get_shifted_transformed_spectrum(self, zpe, deg_states, target_states, hams, shifts, tmoms, handling_mode='transpose'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L709)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L709?message=Update%20Docs)]
</div>


- `zpe`: `Any`
    >No description...
- `deg_states`: `Any`
    >No description...
- `target_states`: `Any`
    >No description...
- `hams`: `Any`
    >No description...
- `shifts`: `Any`
    >No description...
- `tmoms`: `Any`
    >No description...
- `handling_mode`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.shifted_transformed_spectrum" class="docs-object-method">&nbsp;</a> 
```python
shifted_transformed_spectrum(self, deg_states, hams, shifts, return_transformation=False, handling_mode='transpose'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L742)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L742?message=Update%20Docs)]
</div>


- `deg_states`: `Any`
    >No description...
- `hams`: `Any`
    >No description...
- `shifts`: `Any`
    >No description...
- `return_transformation`: `Any`
    >No description...
- `handling_mode`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.transition_data" class="docs-object-method">&nbsp;</a> 
```python
transition_data(self, states, keys=['frequency', 'transition_moment', 'intensity'], data='deperturbed'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L768)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L768?message=Update%20Docs)]
</div>


- `states`: `Any`
    >No description...
- `keys`: `Any`
    >No description...
- `data`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.transition_moment_term_sums" class="docs-object-method">&nbsp;</a> 
```python
transition_moment_term_sums(self, states, terms=None, rotation=None, data='deperturbed'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L833)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L833?message=Update%20Docs)]
</div>


- `states`: `Any`
    >No description...
- `terms`: `Any`
    >No description...
- `rotation`: `Any`
    >No description...
- `data`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.transition_moment_term_sums_first_order" class="docs-object-method">&nbsp;</a> 
```python
transition_moment_term_sums_first_order(self, states, rotation=None, data='deperturbed'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L880)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L880?message=Update%20Docs)]
</div>


- `states`: `Any`
    >No description...
- `rotation`: `Any`
    >No description...
- `data`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.intensity_breakdown" class="docs-object-method">&nbsp;</a> 
```python
intensity_breakdown(self, states, terms=None, data='deperturbed'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L903)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L903?message=Update%20Docs)]
</div>


- `states`: `Any`
    >No description...
- `terms`: `Any`
    >No description...
- `data`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.degenerate_coupling_element" class="docs-object-method">&nbsp;</a> 
```python
degenerate_coupling_element(self, state1, state2): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L928)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L928?message=Update%20Docs)]
</div>


- `state1`: `Any`
    >No description...
- `state2`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.format_deperturbed_hamiltonian" class="docs-object-method">&nbsp;</a> 
```python
format_deperturbed_hamiltonian(self, which): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L951)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L951?message=Update%20Docs)]
</div>


- `which`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.log_parser" class="docs-object-method">&nbsp;</a> 
```python
@property
log_parser(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Analyzer.VPTAnalyzer.print_output_tables" class="docs-object-method">&nbsp;</a> 
```python
print_output_tables(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L980)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L980?message=Update%20Docs)]
</div>

 </div>
</div>




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Analyzer/VPTAnalyzer.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Analyzer/VPTAnalyzer.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Analyzer/VPTAnalyzer.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Analyzer/VPTAnalyzer.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L411?message=Update%20Docs)