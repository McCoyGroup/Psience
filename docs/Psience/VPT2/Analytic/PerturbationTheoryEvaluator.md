## <a id="Psience.VPT2.Analytic.PerturbationTheoryEvaluator">PerturbationTheoryEvaluator</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic.py#L6436)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic.py#L6436?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.VPT2.Analytic.PerturbationTheoryEvaluator.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, solver: Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver, expansion, freqs=None, parallelizer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic.py#L6438)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic.py#L6438?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PerturbationTheoryEvaluator.modify_hamiltonian" class="docs-object-method">&nbsp;</a> 
```python
modify_hamiltonian(self, hamiltonian_corrections): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PerturbationTheoryEvaluator.py#L6445)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PerturbationTheoryEvaluator.py#L6445?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PerturbationTheoryEvaluator.get_energy_corrections" class="docs-object-method">&nbsp;</a> 
```python
get_energy_corrections(self, states, order=None, expansions=None, freqs=None, zero_cutoff=None, degenerate_states=None, verbose=False, logger=None, parallelizer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PerturbationTheoryEvaluator.py#L6457)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PerturbationTheoryEvaluator.py#L6457?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PerturbationTheoryEvaluator.is_single_expansion" class="docs-object-method">&nbsp;</a> 
```python
@staticmethod
is_single_expansion(expansion, min_order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/staticmethod.py#L6492)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/staticmethod.py#L6492?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PerturbationTheoryEvaluator.get_overlap_corrections" class="docs-object-method">&nbsp;</a> 
```python
get_overlap_corrections(self, states, order=None, expansions=None, degenerate_states=None, freqs=None, zero_cutoff=None, verbose=False, parallelizer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PerturbationTheoryEvaluator.py#L6523)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PerturbationTheoryEvaluator.py#L6523?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PerturbationTheoryEvaluator.get_diff_map" class="docs-object-method">&nbsp;</a> 
```python
get_diff_map(self, state_map): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PerturbationTheoryEvaluator.py#L6558)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PerturbationTheoryEvaluator.py#L6558?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PerturbationTheoryEvaluator.get_finals" class="docs-object-method">&nbsp;</a> 
```python
@staticmethod
get_finals(initial, change, perms): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/staticmethod.py#L6580)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/staticmethod.py#L6580?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PerturbationTheoryEvaluator.get_degenerate_changes" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_degenerate_changes(cls, degenerate_pairs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L6827)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L6827?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PerturbationTheoryEvaluator.get_state_by_state_corrections" class="docs-object-method">&nbsp;</a> 
```python
get_state_by_state_corrections(self, generator, states, order=None, terms=None, epaths=None, expansions=None, freqs=None, verbose=False, allowed_coefficients=None, disallowed_coefficients=None, degenerate_states=None, only_degenerate_terms=False, degenerate_correction_generator=None, include_degenerate_correction_terms=True, log_scaled=False, zero_cutoff=None, return_sorted=False, logger=None, parallelizer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PerturbationTheoryEvaluator.py#L6860)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PerturbationTheoryEvaluator.py#L6860?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PerturbationTheoryEvaluator.get_matrix_corrections" class="docs-object-method">&nbsp;</a> 
```python
get_matrix_corrections(self, states, order=None, expansions=None, freqs=None, zero_cutoff=None, verbose=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PerturbationTheoryEvaluator.py#L6910)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PerturbationTheoryEvaluator.py#L6910?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PerturbationTheoryEvaluator.get_full_wavefunction_corrections" class="docs-object-method">&nbsp;</a> 
```python
get_full_wavefunction_corrections(self, states, order=None, expansions=None, freqs=None, zero_cutoff=None, degenerate_states=None, verbose=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PerturbationTheoryEvaluator.py#L6915)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PerturbationTheoryEvaluator.py#L6915?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PerturbationTheoryEvaluator.get_wavefunction_corrections" class="docs-object-method">&nbsp;</a> 
```python
get_wavefunction_corrections(self, states, order=None, expansions=None, freqs=None, zero_cutoff=None, degenerate_states=None, verbose=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PerturbationTheoryEvaluator.py#L6924)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PerturbationTheoryEvaluator.py#L6924?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PerturbationTheoryEvaluator.get_reexpressed_hamiltonian" class="docs-object-method">&nbsp;</a> 
```python
get_reexpressed_hamiltonian(self, states, order=None, expansions=None, freqs=None, degenerate_states=None, only_degenerate_terms=False, verbose=False, include_diagonal=False, hamiltonian_corrections=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PerturbationTheoryEvaluator.py#L6946)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PerturbationTheoryEvaluator.py#L6946?message=Update%20Docs)]
</div>

  - `states`: `Any`
    > 
  - `order`: `Any`
    > 
  - `expansions`: `Any`
    > 
  - `freqs`: `Any`
    > 
  - `degenerate_states`: `Any`
    > 
  - `only_degenerate_terms`: `Any`
    > 
  - `verbose`: `Any`
    > 
  - `include_diagonal`: `Any`
    > 
  - `hamiltonian_corrections`: `Any`
    > `[[(order, terms), expansion], ...]`
  - `opts`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Analytic.PerturbationTheoryEvaluator.get_operator_corrections" class="docs-object-method">&nbsp;</a> 
```python
get_operator_corrections(self, operator_expansion, states, order=None, expansions=None, freqs=None, degenerate_states=None, operator_type=None, check_single=True, terms=None, min_order=1, verbose=False, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PerturbationTheoryEvaluator.py#L7060)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PerturbationTheoryEvaluator.py#L7060?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.PerturbationTheoryEvaluator.evaluate_expressions" class="docs-object-method">&nbsp;</a> 
```python
evaluate_expressions(self, states, exprs, expansions=None, operator_expansions=None, degenerate_states=None, zero_cutoff=None, verbose=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Analytic/PerturbationTheoryEvaluator.py#L7098)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic/PerturbationTheoryEvaluator.py#L7098?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Analytic/PerturbationTheoryEvaluator.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Analytic/PerturbationTheoryEvaluator.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Analytic/PerturbationTheoryEvaluator.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Analytic/PerturbationTheoryEvaluator.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Analytic.py#L6436?message=Update%20Docs)   
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