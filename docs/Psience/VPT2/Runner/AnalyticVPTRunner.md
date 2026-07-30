## <a id="Psience.VPT2.Runner.AnalyticVPTRunner">AnalyticVPTRunner</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L3366)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L3366?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
matrix_formatting_options: dict
hamiltonian_correction_modification_type: str
```
<a id="Psience.VPT2.Runner.AnalyticVPTRunner.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, expansions, order=None, expansion_order=None, freqs=None, internals=True, logger=None, hamiltonian=None, checkpoint=None, dipole_expansion=None, allowed_terms=None, allowed_coefficients=None, disallowed_coefficients=None, allowed_energy_changes=None, intermediate_normalization=None, local_mode_couplings=None, local_mode_coupling_order=None, parallelizer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L3367)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L3367?message=Update%20Docs)]
</div>
**LLM Docstring**

Set up an analytic (symbolic) VPT evaluator, either wrapping an already-built `PerturbationTheoryEvaluator` directly or constructing one from a raw set of Taylor-expansion coefficients and an `AnalyticPerturbationTheorySolver` built for the requested order, and resolving the dipole expansion (from the attached classic Hamiltonian, if any and not given directly) needed for transition-moment/intensity calculations.
  - `expansions`: `list | PerturbationTheoryEvaluator`
    > the raw per-order `[V, G, ...]` expansion coefficients, or an already-built `PerturbationTheoryEvaluator`
  - `order`: `int | None`
    > the perturbation-theory order to build the solver for
  - `expansion_order`: `dict | None`
    > the per-term expansion orders used to resolve the dipole expansion's default order
  - `freqs`: `np.ndarray | None`
    > the mode frequencies
  - `internals`: `bool`
    > whether the expansion is expressed in internal coordinates
  - `logger`: `Logger | None`
    > logger for diagnostics
  - `hamiltonian`: `PerturbationTheoryHamiltonian | None`
    > an associated classic `PerturbationTheoryHamiltonian`, used as a fallback source for the dipole expansion and other molecule-specific data
  - `checkpoint`: `str | Checkpointer | None`
    > checkpoint file/object for caching intermediate symbolic expressions
  - `dipole_expansion`: `list[np.ndarray] | None`
    > explicit dipole-derivative expansion terms to use instead of deriving them from `hamiltonian`
  - `allowed_terms`: `object | None`
    > restrict the analytic solver to only these perturbation term types
  - `allowed_coefficients`: `object | None`
    > restrict the analytic solver to only these expansion coefficients
  - `disallowed_coefficients`: `object | None`
    > exclude these expansion coefficients from the analytic solver
  - `allowed_energy_changes`: `object | None`
    > restrict the analytic solver to only these energy-change patterns
  - `intermediate_normalization`: `bool | None`
    > whether to use intermediate (rather than full) normalization in the analytic solver
  - `local_mode_couplings`: `list | None`
    > extra local-mode coupling terms to inject into the Hamiltonian when running VPT
  - `local_mode_coupling_order`: `int | None`
    > the perturbative order local-mode couplings should be injected at
  - `parallelizer`: `Parallelizer | None`
    > parallelization backend for the evaluator
  - `:returns`: `None`
    > None


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.from_hamiltonian" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_hamiltonian(cls, ham, order, expansion_order=None, logger=None, checkpoint=None, parallelizer=None, allowed_terms=None, allowed_coefficients=None, disallowed_coefficients=None, allowed_energy_changes=None, take_diagonal_v4_terms=True, intermediate_normalization=None, corrected_fundamental_frequencies=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L3463)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L3463?message=Update%20Docs)]
</div>
A driver powered by a classic PerturbationTheoryHamiltonian object
  - `ham`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.construct" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
construct(cls, system, states=None, *, order=2, expressions_file=None, allowed_terms=None, allowed_coefficients=None, disallowed_coefficients=None, allowed_energy_changes=None, mixed_derivative_handling_mode='analytical', degeneracy_specs=None, corrected_fundamental_frequencies=None, parallelizer=None, **settings) -> '(AnalyticVPTRunner, VPTMultiStateSpace)': 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L3547)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L3547?message=Update%20Docs)]
</div>
**LLM Docstring**

Build an `AnalyticVPTRunner` (and, if target states are given, its resolved `MultiVPTStateSpace`) directly from a molecule/system specification: builds a throwaway single-state `VPTRunner` to assemble the classic Hamiltonian, then derives the analytic evaluator from it via `from_hamiltonian`.
  - `system`: `str | list | Molecule | VPTSystem`
    > the molecule or system specification to run VPT on
  - `states`: `object | None`
    > the target states (possibly as `(initial, target)` pairs)
  - `order`: `int`
    > the perturbation-theory order
  - `expressions_file`: `str | None`
    > checkpoint file for caching the analytic expressions
  - `allowed_terms`: `object | None`
    > restrict the analytic solver to only these perturbation term types
  - `allowed_coefficients`: `object | None`
    > restrict the analytic solver to only these expansion coefficients
  - `disallowed_coefficients`: `object | None`
    > exclude these expansion coefficients from the analytic solver
  - `allowed_energy_changes`: `object | None`
    > restrict the analytic solver to only these energy-change patterns
  - `mixed_derivative_handling_mode`: `str`
    > the mixed-derivative handling mode to use when building the classic Hamiltonian
  - `degeneracy_specs`: `object | None`
    > the degeneracy specification(s) to build the `MultiVPTStateSpace` with
  - `corrected_fundamental_frequencies`: `np.ndarray | None`
    > replacement fundamental frequencies to use
  - `parallelizer`: `Parallelizer | None`
    > parallelization backend for the evaluator
  - `settings`: `dict`
    > extra options forwarded to `VPTRunner.construct`
  - `:returns`: `AnalyticVPTRunner | tuple[AnalyticVPTRunner, MultiVPTStateSpace]`
    > the constructed runner, or `(runner, states)` if `states` was given


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.from_file" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_file(cls, file_name, order=2, allowed_terms=None, allowed_coefficients=None, disallowed_coefficients=None, allowed_energy_changes=None, expressions_file=None, **settings): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L3644)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L3644?message=Update%20Docs)]
</div>
**LLM Docstring**

Build an `AnalyticVPTRunner` from a molecule file (e.g. an FChk), via a throwaway single-state `VPTRunner.construct` followed by `from_hamiltonian`.
  - `file_name`: `str`
    > the molecule file to load
  - `order`: `int`
    > the perturbation-theory order
  - `allowed_terms`: `object | None`
    > restrict the analytic solver to only these perturbation term types
  - `allowed_coefficients`: `object | None`
    > restrict the analytic solver to only these expansion coefficients
  - `disallowed_coefficients`: `object | None`
    > exclude these expansion coefficients from the analytic solver
  - `allowed_energy_changes`: `object | None`
    > restrict the analytic solver to only these energy-change patterns
  - `expressions_file`: `str | None`
    > checkpoint file for caching the analytic expressions
  - `settings`: `dict`
    > extra options forwarded to `VPTRunner.construct`
  - `:returns`: `AnalyticVPTRunner`
    > the constructed runner


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.construct_classic_runner" class="docs-object-method">&nbsp;</a> 
```python
construct_classic_runner(self, states, system=None, logger=None, corrected_fundamental_frequencies=None, potential_terms=None, kinetic_terms=None, coriolis_terms=None, pseudopotential_terms=None, dipole_terms=None, initial_states=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3695)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3695?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a classic `VPTRunner` reproducing this analytic evaluator's expansion data (rescaling the stored symbolic-coefficient conventions back into ordinary Taylor-expansion derivative tensors for potential/kinetic/Coriolis/pseudopotential/dipole terms), for cross-validation against the classic (matrix) VPT solver, or for use where a classic runner interface is required.
  - `states`: `object | MultiVPTStateSpace`
    > the target states, or a `MultiVPTStateSpace` (in which case its flattened state space and per-pair initial states are used)
  - `system`: `Molecule | None`
    > the molecule/system to attach; derived from the associated Hamiltonian, or built as a dummy placeholder molecule, if not given
  - `logger`: `Logger | None`
    > logger for the constructed runner; defaults to the associated Hamiltonian's logger
  - `corrected_fundamental_frequencies`: `np.ndarray | None`
    > replacement fundamental frequencies; defaults to this evaluator's own frequencies
  - `potential_terms`: `list[np.ndarray] | None`
    > explicit potential-expansion tensors, bypassing the default rescaling
  - `kinetic_terms`: `list[np.ndarray] | None`
    > explicit kinetic-expansion tensors, bypassing the default rescaling
  - `coriolis_terms`: `list[np.ndarray] | None`
    > explicit Coriolis-expansion tensors, bypassing the default rescaling
  - `pseudopotential_terms`: `list[np.ndarray] | None`
    > explicit pseudopotential-expansion tensors, bypassing the default rescaling
  - `dipole_terms`: `list[np.ndarray] | None`
    > explicit dipole-expansion tensors; defaults to this evaluator's own dipole expansion
  - `initial_states`: `object | None`
    > the initial states for the classic runner; derived from `states` if it's a `MultiVPTStateSpace`
  - `opts`: `dict`
    > extra options forwarded to `VPTRunner.construct`
  - `:returns`: `tuple`
    > `(runner, states)`, the constructed classic runner and its resolved state space


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.clear_caches" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
clear_caches(cls): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L3818)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L3818?message=Update%20Docs)]
</div>
**LLM Docstring**

Clear the global caches used by the underlying `AnalyticPerturbationTheorySolver` (e.g. cached symbolic expressions shared across instances).
  - `:returns`: `None`
    > None


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.prep_multispace" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
prep_multispace(self, states, freqs, system=None, degeneracy_specs=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L3830)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L3830?message=Update%20Docs)]
</div>
**LLM Docstring**

Coerce a raw state specification into a `MultiVPTStateSpace`, passing an already-built one through unchanged and normalizing a bare `[[raise, lower], ...]` polyad-pair specification into the `{'polyads': ...}` dict form expected by `DegeneracySpec`.
  - `states`: `object | MultiVPTStateSpace`
    > the raw state specification, or an already-built `MultiVPTStateSpace`
  - `freqs`: `np.ndarray`
    > the mode frequencies used to build the underlying `VPTStateSpace`(s)
  - `system`: `object | None`
    > the system (molecule/mode context) to associate with the state space
  - `degeneracy_specs`: `object | None`
    > the degeneracy specification(s) to apply
  - `:returns`: `MultiVPTStateSpace`
    > the resolved `MultiVPTStateSpace`


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.prep_states" class="docs-object-method">&nbsp;</a> 
```python
prep_states(self, states, degeneracy_specs=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3878)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3878?message=Update%20Docs)]
</div>
**LLM Docstring**

Coerce a raw state specification into a `MultiVPTStateSpace` using this evaluator's own frequencies and a lightweight dummy system object, via `prep_multispace`.
  - `states`: `object`
    > the raw state specification
  - `degeneracy_specs`: `object | None`
    > the degeneracy specification(s) to apply
  - `:returns`: `MultiVPTStateSpace`
    > the resolved `MultiVPTStateSpace`


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.evaluate_expressions" class="docs-object-method">&nbsp;</a> 
```python
evaluate_expressions(self, states, exprs, zero_cutoff=None, operator_expansions=None, degeneracy_specs=None, verbose=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3899)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3899?message=Update%20Docs)]
</div>
**LLM Docstring**

Evaluate a set of arbitrary symbolic perturbation-theory expressions numerically for the given target states, via the underlying `PerturbationTheoryEvaluator`.
  - `states`: `object`
    > the target states to evaluate the expressions for
  - `exprs`: `object`
    > the symbolic expression(s) to evaluate
  - `zero_cutoff`: `float | None`
    > the magnitude below which a term is treated as exactly zero
  - `operator_expansions`: `object | None`
    > extra operator expansions the expressions may reference
  - `degeneracy_specs`: `object | None`
    > the degeneracy specification(s) to apply when resolving the state space
  - `verbose`: `bool`
    > whether to log detailed evaluation progress
  - `:returns`: `object`
    > the evaluated expression results


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.get_matrix_corrections" class="docs-object-method">&nbsp;</a> 
```python
get_matrix_corrections(self, states, order=None, degeneracy_specs=None, zero_cutoff=None, verbose=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3931)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3931?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the perturbative matrix-element corrections for the given target states, via the underlying `PerturbationTheoryEvaluator`.
  - `states`: `object`
    > the target states to compute corrections for
  - `order`: `int | None`
    > the perturbation-theory order to compute to
  - `degeneracy_specs`: `object | None`
    > the degeneracy specification(s) to apply
  - `zero_cutoff`: `float | None`
    > the magnitude below which a term is treated as exactly zero
  - `verbose`: `bool`
    > whether to log detailed evaluation progress
  - `:returns`: `object`
    > the matrix-element corrections


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.get_energy_corrections" class="docs-object-method">&nbsp;</a> 
```python
get_energy_corrections(self, states, order=None, degeneracy_specs=None, zero_cutoff=None, verbose=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3956)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3956?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the perturbative energy corrections for the given target states, via the underlying `PerturbationTheoryEvaluator`.
  - `states`: `object`
    > the target states to compute energy corrections for
  - `order`: `int | None`
    > the perturbation-theory order to compute to
  - `degeneracy_specs`: `object | None`
    > the degeneracy specification(s) to apply
  - `zero_cutoff`: `float | None`
    > the magnitude below which a term is treated as exactly zero
  - `verbose`: `bool`
    > whether to log detailed evaluation progress
  - `:returns`: `np.ndarray`
    > the per-order energy corrections


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.get_overlap_corrections" class="docs-object-method">&nbsp;</a> 
```python
get_overlap_corrections(self, states, order=None, degeneracy_specs=None, zero_cutoff=None, verbose=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3983)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L3983?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the perturbative wavefunction-overlap corrections for the given target states, via the underlying `PerturbationTheoryEvaluator`.
  - `states`: `object`
    > the target states to compute overlap corrections for
  - `order`: `int | None`
    > the perturbation-theory order to compute to
  - `degeneracy_specs`: `object | None`
    > the degeneracy specification(s) to apply
  - `zero_cutoff`: `float | None`
    > the magnitude below which a term is treated as exactly zero
  - `verbose`: `bool`
    > whether to log detailed evaluation progress
  - `:returns`: `object`
    > the overlap corrections


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.prep_eval_state_pairs" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
prep_eval_state_pairs(cls, states): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L4013)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L4013?message=Update%20Docs)]
</div>
**LLM Docstring**

Flatten a `MultiVPTStateSpace`'s `(initial, final)` state-list pairs into a flat list of `[single_initial_state, final_states]` pairs, one per individual initial state (rather than grouped by pair-block), for direct consumption by the underlying evaluator.
  - `states`: `MultiVPTStateSpace`
    > the multi-state space to flatten
  - `:returns`: `list[list]`
    > the flattened `[initial_state, final_states]` pairs


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.get_full_wavefunction_corrections" class="docs-object-method">&nbsp;</a> 
```python
get_full_wavefunction_corrections(self, states, order=None, degeneracy_specs=None, zero_cutoff=None, verbose=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4030)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4030?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the full (all-component) perturbative wavefunction corrections for the given target states, via the underlying `PerturbationTheoryEvaluator`.
  - `states`: `object`
    > the target states to compute wavefunction corrections for
  - `order`: `int | None`
    > the perturbation-theory order to compute to
  - `degeneracy_specs`: `object | None`
    > the degeneracy specification(s) to apply
  - `zero_cutoff`: `float | None`
    > the magnitude below which a term is treated as exactly zero
  - `verbose`: `bool`
    > whether to log detailed evaluation progress
  - `:returns`: `object`
    > the full wavefunction corrections


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.get_wavefunction_corrections" class="docs-object-method">&nbsp;</a> 
```python
get_wavefunction_corrections(self, states, order=None, degeneracy_specs=None, zero_cutoff=None, verbose=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4061)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4061?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the perturbative wavefunction corrections for the given target states, via the underlying `PerturbationTheoryEvaluator`.
  - `states`: `object`
    > the target states to compute wavefunction corrections for
  - `order`: `int | None`
    > the perturbation-theory order to compute to
  - `degeneracy_specs`: `object | None`
    > the degeneracy specification(s) to apply
  - `zero_cutoff`: `float | None`
    > the magnitude below which a term is treated as exactly zero
  - `verbose`: `bool`
    > whether to log detailed evaluation progress
  - `:returns`: `object`
    > the wavefunction corrections


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.unflatten_corr" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
unflatten_corr(cls, states, corrs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L4093)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L4093?message=Update%20Docs)]
</div>
**LLM Docstring**

Regroup a flat correction result (expressed over the combined initial/final state spaces) back into the per-`(initial, final)`-block structure implied by a `MultiVPTStateSpace`'s state-list pairs.
  - `states`: `MultiVPTStateSpace`
    > the multi-state space whose block structure the flat corrections should be regrouped into
  - `corrs`: `object`
    > the flat correction data (with `initial_states`/`final_states`/`corrections` attributes) to regroup
  - `:returns`: `list[list[np.ndarray]]`
    > the per-block regrouped correction arrays


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.get_operator_corrections" class="docs-object-method">&nbsp;</a> 
```python
get_operator_corrections(self, operator_expansion, states, order=None, terms=None, degeneracy_specs=None, verbose=False, operator_type=None, check_single=True, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4133)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4133?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the perturbative corrections to one or more arbitrary operator expansions for the given target states, via the underlying `PerturbationTheoryEvaluator`, then regroup the flat results back into per-block structure via `unflatten_corr`.
  - `operator_expansion`: `list`
    > the raw operator expansion(s) to evaluate, each a list of derivative tensors
  - `states`: `object`
    > the target states to compute corrections for
  - `order`: `int | None`
    > the perturbation-theory order to compute to
  - `terms`: `object | None`
    > restrict the evaluation to specific term contributions
  - `degeneracy_specs`: `object | None`
    > the degeneracy specification(s) to apply
  - `verbose`: `bool`
    > whether to log detailed evaluation progress
  - `operator_type`: `str | int | None`
    > the operator's symmetry/rank type; inferred from the expansion tensor shapes if not given
  - `check_single`: `bool`
    > whether to auto-detect and wrap a single (rather than multiple) operator expansion
  - `opts`: `dict`
    > extra options forwarded to the underlying evaluator
  - `:returns`: `list`
    > the per-block regrouped operator corrections


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.construct_corrections_vectors" class="docs-object-method">&nbsp;</a> 
```python
construct_corrections_vectors(self, states, corrs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4205)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4205?message=Update%20Docs)]
</div>
**LLM Docstring**

Assemble a set of flat per-order correction matrices spanning the full flat state space, from one or more raw per-block correction results, indexing rows by the union of every block's initial states.
  - `states`: `object`
    > the state space whose flat state list defines the column indexing
  - `corrs`: `object | list`
    > one (or a list of) raw correction result(s), each with `initial_states`/`final_states`/`corrections`
  - `:returns`: `tuple`
    > `(init_states, op_corr_mats)` -- the union of initial states used for row indexing, and the assembled per-order correction matrices (single set if `corrs` was a single result)


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.construct_corrections_matrix" class="docs-object-method">&nbsp;</a> 
```python
construct_corrections_matrix(self, group, corrs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4247)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4247?message=Update%20Docs)]
</div>
**LLM Docstring**

Assemble a set of square per-order correction matrices restricted to a single group of states, from one or more raw per-block correction results.
  - `group`: `object`
    > the group of states (typically a degenerate block) to build the matrix over
  - `corrs`: `object | list`
    > one (or a list of) raw correction result(s), each with `initial_states`/`final_states`/`corrections`
  - `:returns`: `list`
    > the assembled per-order correction matrices, restricted to `group` (single set if `corrs` was a single result)


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.get_transition_moment_corrections" class="docs-object-method">&nbsp;</a> 
```python
get_transition_moment_corrections(self, states, dipole_expansion=None, order=None, degeneracy_specs=None, axes=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4286)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4286?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the perturbative transition-dipole-moment corrections for the given target states, using this evaluator's dipole expansion (or an explicitly given one), scaled by the appropriate factorial normalization, via `get_operator_corrections`.
  - `states`: `object`
    > the target states to compute transition moments for
  - `dipole_expansion`: `list[np.ndarray] | None`
    > an explicit dipole-derivative expansion to use instead of `self.dipole_expansion`
  - `order`: `int | None`
    > the perturbation-theory order to compute to; inferred from `self.expansion_order`/the dipole expansion length if not given
  - `degeneracy_specs`: `object | None`
    > the degeneracy specification(s) to apply
  - `axes`: `Iterable[int] | None`
    > which Cartesian axes to compute transition moments for; defaults to all three
  - `opts`: `dict`
    > extra options forwarded to `get_operator_corrections`
  - `:returns`: `list`
    > the per-block regrouped transition-moment corrections


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.get_freqs" class="docs-object-method">&nbsp;</a> 
```python
get_freqs(self, states, order=None, degeneracy_specs=None, return_corrections=False, verbose=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4336)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4336?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the vibrational transition frequencies (in wavenumbers) for the given target states, relative to the ground/reference state.
  - `states`: `object`
    > the target states to compute frequencies for
  - `order`: `int | None`
    > the perturbation-theory order to compute to
  - `degeneracy_specs`: `object | None`
    > the degeneracy specification(s) to apply
  - `return_corrections`: `bool`
    > whether to also return the raw per-order energy corrections
  - `verbose`: `bool`
    > whether to log detailed evaluation progress
  - `:returns`: `tuple`
    > `(zpe, freqs)`, or `((zpe, freqs), corrections)` if `return_corrections` is set


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.get_reexpressed_hamiltonian" class="docs-object-method">&nbsp;</a> 
```python
get_reexpressed_hamiltonian(self, states, order=None, degeneracy_specs=None, only_degenerate_terms=True, verbose=False, hamiltonian_corrections=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4365)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4365?message=Update%20Docs)]
</div>
**LLM Docstring**

Build the deperturbed (degenerate-block) Hamiltonian matrices for each degenerate group of states, re-expressing the analytic Hamiltonian corrections in the basis of each group and summing them into square correction matrices.
  - `states`: `object`
    > the target states (used to resolve degenerate groupings)
  - `order`: `int | None`
    > the perturbation-theory order to compute to
  - `degeneracy_specs`: `object | None`
    > the degeneracy specification(s) to apply
  - `only_degenerate_terms`: `bool`
    > whether to include only the strictly degenerate-coupling terms in the reexpressed Hamiltonian
  - `verbose`: `bool`
    > whether to log detailed evaluation progress
  - `hamiltonian_corrections`: `object | None`
    > extra Hamiltonian corrections to fold in
  - `opts`: `dict`
    > extra options forwarded to the underlying evaluator
  - `:returns`: `tuple | None`
    > `(all_mats, all_corrs)` -- the summed Hamiltonian matrix and the individual per-order correction matrices, one entry per degenerate group; `None` if there are no degenerate states


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.get_wfc_test_states" class="docs-object-method">&nbsp;</a> 
```python
get_wfc_test_states(self, input_states: Psience.BasisReps.StateSpaces.BasisStateSpace, energy_window): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4415)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4415?message=Update%20Docs)]
</div>
**LLM Docstring**

Identify the candidate states that could plausibly be strongly coupled (via wavefunction-correction magnitude) to a given set of input states, based purely on energy proximity (within `energy_window`) and the maximum possible change in quantum numbers implied by the Hamiltonian expansion's term structure.
  - `input_states`: `BasisStateSpace`
    > the states to find coupling candidates for
  - `energy_window`: `float`
    > the energy window (around each input state's own energy) to search within
  - `:returns`: `list[BasisStateSpace]`
    > the candidate coupling states for each input state


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.get_test_wfn_corrs" class="docs-object-method">&nbsp;</a> 
```python
get_test_wfn_corrs(self, input_states: Psience.BasisReps.StateSpaces.BasisStateSpace, energy_window): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4454)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4454?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4471)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4471?message=Update%20Docs)]
</div>
**LLM Docstring**

Format a table of state energies/frequencies alongside their per-order corrections, converting to wavenumbers and displaying every non-zero-point-energy state as a frequency shift relative to the ZPE.
  - `states`: `MultiVPTStateSpace`
    > the states the energies correspond to
  - `energies`: `np.ndarray`
    > the total energies (in Hartrees), one per state
  - `energy_corrections`: `np.ndarray`
    > the per-order energy corrections (in Hartrees)
  - `zpe_pos`: `int`
    > the index of the zero-point-energy (reference) state
  - `number_format`: `str`
    > the format spec used for each numeric column
  - `:returns`: `str`
    > the formatted table


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.format_degenerate_energies_table" class="docs-object-method">&nbsp;</a> 
```python
format_degenerate_energies_table(self, states, energies, deperturbed_energies, zpe_pos, number_format='.3f'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4522)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4522?message=Update%20Docs)]
</div>
**LLM Docstring**

Format a table comparing each state's degenerate-perturbation-theory-corrected energy/frequency against its deperturbed counterpart, both relative to the zero-point energy.
  - `states`: `MultiVPTStateSpace`
    > the states the energies correspond to
  - `energies`: `np.ndarray`
    > the degenerate-corrected total energies (in Hartrees)
  - `deperturbed_energies`: `np.ndarray`
    > the deperturbed total energies (in Hartrees)
  - `zpe_pos`: `int`
    > the index of the zero-point-energy (reference) state
  - `number_format`: `str`
    > the format spec used for each numeric column
  - `:returns`: `str`
    > the formatted table


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.format_transition_moment_table" class="docs-object-method">&nbsp;</a> 
```python
format_transition_moment_table(self, states, transition_moments, transition_moment_corrections, number_format='.8f'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4568)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4568?message=Update%20Docs)]
</div>
**LLM Docstring**

Format a table of transition-dipole moments (and their per-order corrections) for each initial-state block in a `MultiVPTStateSpace`, with a labeled header separating each initial state's sub-table when there's more than one.
  - `states`: `MultiVPTStateSpace`
    > the initial/final state pairs the transition moments correspond to
  - `transition_moments`: `list`
    > the per-axis, per-block final transition moments
  - `transition_moment_corrections`: `list`
    > the per-axis, per-block, per-order transition-moment corrections
  - `number_format`: `str`
    > the format spec used for each numeric column
  - `:returns`: `str`
    > the formatted table(s), concatenated with initial-state header separators


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.format_operators_table" class="docs-object-method">&nbsp;</a> 
```python
format_operators_table(self, states, keys, operator_values, operator_corrections, number_format='.8f'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4633)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4633?message=Update%20Docs)]
</div>
**LLM Docstring**

Format a table of arbitrary operator expectation values (and their per-order corrections) for each initial-state block in a `MultiVPTStateSpace`, with a labeled header separating each initial state's sub-table when there's more than one.
  - `states`: `MultiVPTStateSpace`
    > the initial/final state pairs the operator values correspond to
  - `keys`: `Iterable`
    > the operator names/labels, used as column headers
  - `operator_values`: `list`
    > the per-operator, per-block final operator values
  - `operator_corrections`: `list`
    > the per-operator, per-block, per-order operator corrections
  - `number_format`: `str`
    > the format spec used for each numeric column
  - `:returns`: `str`
    > the formatted table(s), concatenated with initial-state header separators


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.format_spectrum_table" class="docs-object-method">&nbsp;</a> 
```python
format_spectrum_table(self, states, harmonic_spectra, spectra, deperturbed_spectra=None, number_format='.3f'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4711)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4711?message=Update%20Docs)]
</div>
**LLM Docstring**

Format a table of harmonic, anharmonic, and (optionally) deperturbed IR spectra for each initial-state block in a `MultiVPTStateSpace`, with a labeled header separating each initial state's sub-table when there's more than one.
  - `states`: `MultiVPTStateSpace`
    > the initial/final state pairs the spectra correspond to
  - `harmonic_spectra`: `list`
    > the per-block harmonic `DiscreteSpectrum` objects
  - `spectra`: `list`
    > the per-block anharmonic `DiscreteSpectrum` objects
  - `deperturbed_spectra`: `list | None`
    > the per-block deperturbed `DiscreteSpectrum` objects, if degenerate treatment was applied
  - `number_format`: `str`
    > the format spec used for each numeric column
  - `:returns`: `str`
    > the formatted table(s), concatenated with initial-state header separators


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.prep_operators" class="docs-object-method">&nbsp;</a> 
```python
prep_operators(self, operator_expansions, operator_terms, order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4777)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4777?message=Update%20Docs)]
</div>
**LLM Docstring**

Normalize a user-supplied operator specification (raw expansion coefficients, in either list or named-dict form) into the fully-expanded per-mode-basis operator terms needed for correction calculations, using `self.ham.prep_operator_terms` to expand raw finite coefficients if pre-expanded `operator_terms` weren't already given directly.
  - `operator_expansions`: `object | None`
    > raw operator expansion coefficients (a single operator's coefficient list, a list of such lists, or a name-keyed dict of them), used if `operator_terms` isn't given
  - `operator_terms`: `object | None`
    > already-expanded per-mode-basis operator terms (a single list, a list of lists, or a name-keyed dict), used directly if given
  - `order`: `int | None`
    > the expansion order to expand `operator_expansions` to; inferred from the coefficient list length if not given
  - `:returns`: `tuple`
    > `(keys, operator_terms)` -- the resolved operator names/indices and their expanded per-mode-basis terms


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.format_matrix" class="docs-object-method">&nbsp;</a> 
```python
format_matrix(self, ham): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4835)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4835?message=Update%20Docs)]
</div>
**LLM Docstring**

Format a matrix as a plain-text string using this class's standard print options (fixed precision, no truncation, suppressed scientific notation), stripping the outer bracket characters for a cleaner look.
  - `ham`: `np.ndarray | object`
    > the matrix to format
  - `:returns`: `str`
    > the formatted matrix string


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.modify_hamiltonian" class="docs-object-method">&nbsp;</a> 
```python
modify_hamiltonian(self, hamiltonian_corrections): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4850)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4850?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a new `AnalyticVPTRunner` whose underlying evaluator has extra Hamiltonian corrections applied, via `PerturbationTheoryEvaluator.modify_hamiltonian`, preserving this runner's order/expansion-order/dipole-expansion/logger settings.
  - `hamiltonian_corrections`: `object`
    > the Hamiltonian corrections to apply (e.g. local-mode coupling terms)
  - `:returns`: `AnalyticVPTRunner`
    > the new runner with the modified Hamiltonian


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.run_VPT" class="docs-object-method">&nbsp;</a> 
```python
run_VPT(self, states, calculate_intensities=True, operator_expansions=None, operator_terms=None, operator_type=None, order=None, verbose=False, degeneracy_specs=None, handle_degeneracies=True, zero_cutoff=None, transition_moment_terms=None, hamiltonian_corrections=None, clear_caches=True, hamiltonian_correction_type=None, only_degenerate_terms=True, force_return_on_crash=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4909)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/AnalyticVPTRunner.py#L4909?message=Update%20Docs)]
</div>
**LLM Docstring**

Top-level entry point for running a full analytic VPT calculation on a set of target states: optionally injects local-mode coupling terms (and/or applies other Hamiltonian corrections, either up front via `modify_hamiltonian` or as post-hoc corrections during the solve), resolves the degenerate-state structure, computes energies/wavefunction corrections/transition moments/operator values (with strong-coupling-based degenerate-state extension if requested), and assembles the results into formatted output tables and spectra.
  - `states`: `object`
    > the target states (and, optionally, initial states) to compute
  - `calculate_intensities`: `bool`
    > whether to compute transition moments/IR intensities
  - `operator_expansions`: `object | None`
    > extra operator expansions to additionally evaluate corrections for
  - `operator_terms`: `object | None`
    > pre-expanded operator terms to additionally evaluate corrections for
  - `operator_type`: `str | int | None`
    > the operator symmetry/rank type, forwarded to the correction evaluator
  - `order`: `int | None`
    > the perturbation-theory order to compute to
  - `verbose`: `bool`
    > whether to log detailed evaluation progress
  - `degeneracy_specs`: `object | None`
    > the degeneracy specification(s) to apply
  - `handle_degeneracies`: `bool`
    > whether to apply degenerate-perturbation-theory handling at all
  - `zero_cutoff`: `float | None`
    > the magnitude below which a term is treated as exactly zero
  - `transition_moment_terms`: `object | None`
    > restrict the transition-moment evaluation to specific term contributions
  - `hamiltonian_corrections`: `object | None`
    > extra Hamiltonian corrections to apply, either up front (`'primary'`) or as post-hoc adjustments (`'degenerate'`), per `hamiltonian_correction_type`
  - `clear_caches`: `bool`
    > whether to clear the global analytic-solver caches before running
  - `hamiltonian_correction_type`: `str | None`
    > `'primary'` to fold `hamiltonian_corrections` into the Hamiltonian itself before solving, or `'degenerate'`/other to apply them only within degenerate blocks; defaults to `self.hamiltonian_correction_modification_type`
  - `only_degenerate_terms`: `bool`
    > whether the reexpressed Hamiltonian should include only strictly degenerate-coupling terms
  - `force_return_on_crash`: `bool`
    > whether to catch exceptions during the run and still return whatever partial results were computed, rather than propagating the error
  - `:returns`: `object`
    > the computed VPT results (energies, wavefunction data, spectra, and formatted tables)


<a id="Psience.VPT2.Runner.AnalyticVPTRunner.run_simple" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
run_simple(cls, system, states, calculate_intensities=True, operator_expansions=None, operator_terms=None, operator_type=None, verbose=False, return_runner=False, degeneracy_specs=None, degeneracy_states=None, handle_degeneracies=True, zero_cutoff=None, clear_caches=True, hamiltonian_correction_type=None, hamiltonian_corrections=None, only_degenerate_terms=True, force_return_on_crash=True, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L5188)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L5188?message=Update%20Docs)]
</div>
**LLM Docstring**

Convenience one-shot entry point: builds an `AnalyticVPTRunner` for the given system (via `construct`) and immediately runs the full VPT calculation on it (via `run_VPT`).
  - `system`: `str | list | Molecule | VPTSystem`
    > the molecule or system specification to run VPT on
  - `states`: `object`
    > the target states (and, optionally, initial states) to compute
  - `calculate_intensities`: `bool`
    > whether to compute transition moments/IR intensities
  - `operator_expansions`: `object | None`
    > extra operator expansions to additionally evaluate corrections for
  - `operator_terms`: `object | None`
    > pre-expanded operator terms to additionally evaluate corrections for
  - `operator_type`: `str | int | None`
    > the operator symmetry/rank type
  - `verbose`: `bool`
    > whether to log detailed evaluation progress
  - `return_runner`: `bool`
    > whether to also return the constructed runner alongside the results
  - `degeneracy_specs`: `object | None`
    > the degeneracy specification(s) to apply
  - `handle_degeneracies`: `bool`
    > whether to apply degenerate-perturbation-theory handling at all
  - `zero_cutoff`: `float | None`
    > the magnitude below which a term is treated as exactly zero
  - `clear_caches`: `bool`
    > whether to clear the global analytic-solver caches before running
  - `hamiltonian_correction_type`: `str | None`
    > how to apply `hamiltonian_corrections`, forwarded to `run_VPT`
  - `hamiltonian_corrections`: `object | None`
    > extra Hamiltonian corrections to apply
  - `only_degenerate_terms`: `bool`
    > whether the reexpressed Hamiltonian should include only strictly degenerate-coupling terms
  - `force_return_on_crash`: `bool`
    > whether to catch exceptions during the run and still return partial results
  - `opts`: `dict`
    > extra options forwarded to `construct`
  - `:returns`: `object | tuple`
    > the computed VPT results, or `(runner, results)` if `return_runner` is set
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
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L3366?message=Update%20Docs)   
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