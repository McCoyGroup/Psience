## <a id="Psience.VPT2.Analyzer.VPTResultsLoader">VPTResultsLoader</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L50)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L50?message=Update%20Docs)]
</div>

Provides tools for loading results into canonical
sources from a simulation, both from checkpoint files and from
`PerturbationTheoryWavefunctions` objects and potentially more







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.VPT2.Analyzer.VPTResultsLoader.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, res, res_type=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer.py#L57)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L57?message=Update%20Docs)]
</div>

  - `res`: `Any`
    > 
  - `res_type`: `Any`
    >


<a id="Psience.VPT2.Analyzer.VPTResultsLoader.get_res_type" class="docs-object-method">&nbsp;</a> 
```python
get_res_type(self, res): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L69)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L69?message=Update%20Docs)]
</div>

  - `res`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTResultsLoader.potential_terms" class="docs-object-method">&nbsp;</a> 
```python
potential_terms(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L43)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L43?message=Update%20Docs)]
</div>
Returns the expansion of the potential
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTResultsLoader.kinetic_terms" class="docs-object-method">&nbsp;</a> 
```python
kinetic_terms(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L43)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L43?message=Update%20Docs)]
</div>
Returns the expansion of the kinetic energy
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTResultsLoader.dipole_terms" class="docs-object-method">&nbsp;</a> 
```python
dipole_terms(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L43)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L43?message=Update%20Docs)]
</div>
Returns the expansion of the dipole moment
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTResultsLoader.basis" class="docs-object-method">&nbsp;</a> 
```python
basis(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L43)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L43?message=Update%20Docs)]
</div>
Returns the basis for the calculation
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTResultsLoader.target_states" class="docs-object-method">&nbsp;</a> 
```python
target_states(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L43)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L43?message=Update%20Docs)]
</div>
Returns the target states for the calculation
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTResultsLoader.spectrum" class="docs-object-method">&nbsp;</a> 
```python
spectrum(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L43)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L43?message=Update%20Docs)]
</div>
Returns the IR spectrum calculated from perturbation theory
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTResultsLoader.zero_order_spectrum" class="docs-object-method">&nbsp;</a> 
```python
zero_order_spectrum(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L43)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L43?message=Update%20Docs)]
</div>
Returns the zero-order IR spectrum calculated from perturbation theory
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTResultsLoader.energy_corrections" class="docs-object-method">&nbsp;</a> 
```python
energy_corrections(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L43)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L43?message=Update%20Docs)]
</div>
Returns the corrections to the energies
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTResultsLoader.energies" class="docs-object-method">&nbsp;</a> 
```python
energies(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L236)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L236?message=Update%20Docs)]
</div>

  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTResultsLoader.wavefunction_corrections" class="docs-object-method">&nbsp;</a> 
```python
wavefunction_corrections(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L43)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L43?message=Update%20Docs)]
</div>
Returns the corrections to the wavefunctions
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTResultsLoader.transition_moment_corrections" class="docs-object-method">&nbsp;</a> 
```python
transition_moment_corrections(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L43)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L43?message=Update%20Docs)]
</div>
Returns the corrections to the wavefunctions
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTResultsLoader.transition_moments" class="docs-object-method">&nbsp;</a> 
```python
transition_moments(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L275)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L275?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analyzer.VPTResultsLoader.deperturbed_transition_moment_corrections" class="docs-object-method">&nbsp;</a> 
```python
deperturbed_transition_moment_corrections(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L43)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L43?message=Update%20Docs)]
</div>
Returns the corrections to the wavefunctions
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTResultsLoader.deperturbed_transition_moments" class="docs-object-method">&nbsp;</a> 
```python
deperturbed_transition_moments(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L304)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L304?message=Update%20Docs)]
</div>

  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTResultsLoader.degenerate_states" class="docs-object-method">&nbsp;</a> 
```python
degenerate_states(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L43)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L43?message=Update%20Docs)]
</div>
Returns the deperturbed states used to make the degenerate transform
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTResultsLoader.deperturbed_hamiltonians" class="docs-object-method">&nbsp;</a> 
```python
deperturbed_hamiltonians(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L43)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L43?message=Update%20Docs)]
</div>
Returns the deperturbed Hamiltonians used to make the degenerate transform
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTResultsLoader.degenerate_energies" class="docs-object-method">&nbsp;</a> 
```python
degenerate_energies(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L43)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L43?message=Update%20Docs)]
</div>
Returns the deperturbed states used to make the degenerate transform
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTResultsLoader.degenerate_rotations" class="docs-object-method">&nbsp;</a> 
```python
degenerate_rotations(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L43)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L43?message=Update%20Docs)]
</div>
Returns the deperturbed states used to make the degenerate transform
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analyzer.VPTResultsLoader.log_file" class="docs-object-method">&nbsp;</a> 
```python
log_file(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L43)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer/VPTResultsLoader.py#L43?message=Update%20Docs)]
</div>
Returns the log_file for the run
  - `:returns`: `_`
    >
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Analyzer/VPTResultsLoader.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Analyzer/VPTResultsLoader.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Analyzer/VPTResultsLoader.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Analyzer/VPTResultsLoader.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analyzer.py#L50?message=Update%20Docs)   
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