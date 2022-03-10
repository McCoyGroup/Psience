## <a id="Psience.VPT2.Analyzer.VPTResultsLoader">VPTResultsLoader</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Analyzer.py#L46)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Analyzer.py#L46?message=Update%20Docs)]
</div>

Provides tools for loading results into canonical
sources from a simulation, both from checkpoint files and from
`PerturbationTheoryWavefunctions` objects and potentially more

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

<a id="Psience.VPT2.Analyzer.VPTResultsLoader.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, res, res_type=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Analyzer.py#L53)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Analyzer.py#L53?message=Update%20Docs)]
</div>


- `res`: `Any`
    >No description...
- `res_type`: `Any`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTResultsLoader.get_res_type" class="docs-object-method">&nbsp;</a> 
```python
get_res_type(self, res): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Analyzer.py#L65)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Analyzer.py#L65?message=Update%20Docs)]
</div>


- `res`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTResultsLoader.potential_terms" class="docs-object-method">&nbsp;</a> 
```python
potential_terms(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Analyzer.py#L39)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Analyzer.py#L39?message=Update%20Docs)]
</div>

Returns the expansion of the potential
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTResultsLoader.kinetic_terms" class="docs-object-method">&nbsp;</a> 
```python
kinetic_terms(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Analyzer.py#L39)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Analyzer.py#L39?message=Update%20Docs)]
</div>

Returns the expansion of the kinetic energy
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTResultsLoader.dipole_terms" class="docs-object-method">&nbsp;</a> 
```python
dipole_terms(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Analyzer.py#L39)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Analyzer.py#L39?message=Update%20Docs)]
</div>

Returns the expansion of the dipole moment
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTResultsLoader.basis" class="docs-object-method">&nbsp;</a> 
```python
basis(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Analyzer.py#L39)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Analyzer.py#L39?message=Update%20Docs)]
</div>

Returns the basis for the calculation
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTResultsLoader.target_states" class="docs-object-method">&nbsp;</a> 
```python
target_states(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Analyzer.py#L39)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Analyzer.py#L39?message=Update%20Docs)]
</div>

Returns the target states for the calculation
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTResultsLoader.spectrum" class="docs-object-method">&nbsp;</a> 
```python
spectrum(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Analyzer.py#L39)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Analyzer.py#L39?message=Update%20Docs)]
</div>

Returns the IR spectrum calculated from perturbation theory
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTResultsLoader.zero_order_spectrum" class="docs-object-method">&nbsp;</a> 
```python
zero_order_spectrum(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Analyzer.py#L39)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Analyzer.py#L39?message=Update%20Docs)]
</div>

Returns the zero-order IR spectrum calculated from perturbation theory
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTResultsLoader.energy_corrections" class="docs-object-method">&nbsp;</a> 
```python
energy_corrections(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Analyzer.py#L39)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Analyzer.py#L39?message=Update%20Docs)]
</div>

Returns the corrections to the energies
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTResultsLoader.energies" class="docs-object-method">&nbsp;</a> 
```python
energies(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Analyzer.py#L220)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Analyzer.py#L220?message=Update%20Docs)]
</div>


- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTResultsLoader.wavefunction_corrections" class="docs-object-method">&nbsp;</a> 
```python
wavefunction_corrections(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Analyzer.py#L39)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Analyzer.py#L39?message=Update%20Docs)]
</div>

Returns the corrections to the wavefunctions
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTResultsLoader.transition_moment_corrections" class="docs-object-method">&nbsp;</a> 
```python
transition_moment_corrections(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Analyzer.py#L39)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Analyzer.py#L39?message=Update%20Docs)]
</div>

Returns the corrections to the wavefunctions
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTResultsLoader.transition_moments" class="docs-object-method">&nbsp;</a> 
```python
transition_moments(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Analyzer.py#L259)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Analyzer.py#L259?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Analyzer.VPTResultsLoader.deperturbed_transition_moment_corrections" class="docs-object-method">&nbsp;</a> 
```python
deperturbed_transition_moment_corrections(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Analyzer.py#L39)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Analyzer.py#L39?message=Update%20Docs)]
</div>

Returns the corrections to the wavefunctions
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTResultsLoader.deperturbed_transition_moments" class="docs-object-method">&nbsp;</a> 
```python
deperturbed_transition_moments(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Analyzer.py#L288)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Analyzer.py#L288?message=Update%20Docs)]
</div>


- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTResultsLoader.degenerate_states" class="docs-object-method">&nbsp;</a> 
```python
degenerate_states(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Analyzer.py#L39)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Analyzer.py#L39?message=Update%20Docs)]
</div>

Returns the deperturbed states used to make the degenerate transform
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTResultsLoader.deperturbed_hamiltonians" class="docs-object-method">&nbsp;</a> 
```python
deperturbed_hamiltonians(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Analyzer.py#L39)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Analyzer.py#L39?message=Update%20Docs)]
</div>

Returns the deperturbed Hamiltonians used to make the degenerate transform
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTResultsLoader.degenerate_energies" class="docs-object-method">&nbsp;</a> 
```python
degenerate_energies(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Analyzer.py#L39)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Analyzer.py#L39?message=Update%20Docs)]
</div>

Returns the deperturbed states used to make the degenerate transform
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Analyzer.VPTResultsLoader.degenerate_rotations" class="docs-object-method">&nbsp;</a> 
```python
degenerate_rotations(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Analyzer.py#L39)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Analyzer.py#L39?message=Update%20Docs)]
</div>

Returns the deperturbed states used to make the degenerate transform
- `:returns`: `_`
    >No description...

 </div>
</div>




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Analyzer/VPTResultsLoader.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Analyzer/VPTResultsLoader.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Analyzer/VPTResultsLoader.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Analyzer/VPTResultsLoader.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Analyzer.py#L46?message=Update%20Docs)