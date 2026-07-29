## <a id="Psience.VPT2.Solver.PerturbationTheorySolver">PerturbationTheorySolver</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver.py#L27)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver.py#L27?message=Update%20Docs)]
</div>

A solver that applies perturbation theory
given a series of corrections and population of states.
Supports degenerate and non-degenerate PT.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
PastIndexableTuple: PastIndexableTuple
StateSpaceWrapper: StateSpaceWrapper
ProjectionOperatorWrapper: ProjectionOperatorWrapper
ProjectedOperator: ProjectedOperator
default_strong_coupling_threshold: float
PTResults: PTResults
```
<a id="Psience.VPT2.Solver.PerturbationTheorySolver.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, perturbations, states, coupled_states=None, order=2, total_space=None, flat_total_space=None, state_space_iterations=None, state_space_terms=None, state_space_filters=None, extended_state_space_filter_generator=None, extended_state_space_postprocessor=None, target_property_rules=None, allow_sakurai_degs=False, allow_post_PT_calc=True, modify_degenerate_perturbations=False, gaussian_resonance_handling=False, ignore_odd_order_energies=False, intermediate_normalization=False, reorthogonalize_degenerate_states=None, check_overlap=True, zero_element_warning=True, degenerate_states=None, degeneracy_handlers=None, handle_strong_couplings=False, local_coupling_hamiltonian=None, local_coupling_order=None, low_frequency_mode_cutoff=0.00115, zero_order_energy_corrections=None, nondeg_hamiltonian_precision=3, memory_constrained=False, keep_hamiltonians=None, logger=None, parallelizer=None, checkpointer=None, results=None, checkpoint_keys=None, use_cached_representations=False, use_cached_basis=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver.py#L34)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver.py#L34?message=Update%20Docs)]
</div>

  - `perturbations`: `Iterable[Representation]`
    > 
  - `states`: `BasisStateSpace`
    > 
  - `coupled_states`: `BasisMultiStateSpace`
    > 
  - `order`: `Any`
    > 
  - `degenerate_states`: `Any`
    > 
  - `degeneracy_mode`: `Any`
    > 
  - `logger`: `Any`
    > 
  - `parallelizer`: `Any`
    > 
  - `checkpointer`: `Any`
    > 
  - `results`: `Any`
    >


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.coupled_states" class="docs-object-method">&nbsp;</a> 
```python
@property
coupled_states(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L208)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L208?message=Update%20Docs)]
</div>

  - `:returns`: `_`
    >


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.total_space_dim" class="docs-object-method">&nbsp;</a> 
```python
@property
total_space_dim(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L218)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L218?message=Update%20Docs)]
</div>

  - `:returns`: `_`
    >


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.flat_total_space" class="docs-object-method">&nbsp;</a> 
```python
@property
flat_total_space(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L228)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L228?message=Update%20Docs)]
</div>

  - `:returns`: `_`
    >


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.total_state_space" class="docs-object-method">&nbsp;</a> 
```python
@property
total_state_space(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L238)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L238?message=Update%20Docs)]
</div>

  - `:returns`: `_`
    >


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.representations" class="docs-object-method">&nbsp;</a> 
```python
@property
representations(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L282)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L282?message=Update%20Docs)]
</div>

  - `:returns`: `Iterable[SparseArray]`
    >


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.degenerate_spaces" class="docs-object-method">&nbsp;</a> 
```python
@property
degenerate_spaces(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L303)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L303?message=Update%20Docs)]
</div>

  - `:returns`: `_`
    >


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.merge_deg_spaces" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
merge_deg_spaces(cls, new_states): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L350)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L350?message=Update%20Docs)]
</div>
**LLM Docstring**

Combine several independently-identified sets of degenerate state groups into one consistent set, by flattening every group (from either raw `BasisStateSpace`s or `DegenerateMultiStateSpace`s) down to its excitation vectors and merging any that share states via `DegeneracySpec.merge_state_blocks`.
  - `new_states`: `list`
    > the separate sets of degenerate groups to merge, each either a `DegenerateMultiStateSpace` or an iterable of raw state blocks
  - `:returns`: `list[np.ndarray]`
    > the merged degenerate groups


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.zero_order_energies" class="docs-object-method">&nbsp;</a> 
```python
@property
zero_order_energies(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L382)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L382?message=Update%20Docs)]
</div>

  - `:returns`: `_`
    >


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.apply_VPT" class="docs-object-method">&nbsp;</a> 
```python
apply_VPT(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L426)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L426?message=Update%20Docs)]
</div>
Applies perturbation theory to the held basis of states using the
built representations and degenerate state spaces
  - `:returns`: `PerturbationTheoryCorrections`
    >


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.get_VPT_representations" class="docs-object-method">&nbsp;</a> 
```python
get_VPT_representations(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L457)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L457?message=Update%20Docs)]
</div>
Gets the sparse representations of the passed perturbation inside the basis of coupled states.
  - `:returns`: `Iterable[SparseArray]`
    >


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.extend_VPT_representations" class="docs-object-method">&nbsp;</a> 
```python
extend_VPT_representations(self, new_flat_space, new_states): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L533)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L533?message=Update%20Docs)]
</div>
**LLM Docstring**

Extend the cached Hamiltonian perturbation representation matrices to cover a newly added block of coupled states, by computing just the new matrix elements (between the new states and the full flat space) and splicing them into the existing sparse representation data rather than recomputing everything from scratch.
  - `new_flat_space`: `BasisStateSpace`
    > the newly added states, flattened, that the zeroth-order representation needs new diagonal elements for
  - `new_states`: `list[BasisStateSpace]`
    > the newly added coupled-state blocks, one per higher-order perturbation, that need new off-diagonal elements computed
  - `:returns`: `list[SparseArray]`
    > the extended list of representation matrices


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.load_state_spaces" class="docs-object-method">&nbsp;</a> 
```python
load_state_spaces(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L734)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L734?message=Update%20Docs)]
</div>

  - `:returns`: `_`
    >


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.extend_state_spaces" class="docs-object-method">&nbsp;</a> 
```python
extend_state_spaces(self, new_targets, degenerate_states=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L870)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L870?message=Update%20Docs)]
</div>
**LLM Docstring**

Extend the solver's state spaces to additionally include a set of new target states: (re)computes the coupled-state blocks needed for the new targets (applying any configured extended-state-space filter/postprocessor), merges them into the existing per-order coupled-state spaces and flattened total space, and updates `self.states`/`self._total_space` accordingly.
  - `new_targets`: `BasisStateSpace`
    > the new target states to add to the solve
  - `degenerate_states`: `object | None`
    > the degenerate-state groupings driving this extension, forwarded to `extended_state_space_postprocessor` if one is configured
  - `:returns`: `tuple | None`
    > `(flat_space, new_spaces)` -- the newly added (deduplicated) flat states and the per-order newly added coupled-state blocks, or `None` if nothing new was actually found


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.load_coupled_spaces" class="docs-object-method">&nbsp;</a> 
```python
load_coupled_spaces(self, degenerate_spaces=None, spaces=None, wavefunction_terms=None, property_filter=None, filter_spaces=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L951)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L951?message=Update%20Docs)]
</div>
Determines which states need to be coupled at which levels of correction
to handle the PT
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.get_coupled_space" class="docs-object-method">&nbsp;</a> 
```python
get_coupled_space(self, input_state_space, degenerate_space, use_second_deg, allow_PT_degs=True, wavefunction_terms=None, spaces=None, property_filter=None, filter_spaces=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L1505)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L1505?message=Update%20Docs)]
</div>
Applies the VPT equations semi-symbolically, dispatching based on how many
degeneracies we need to handle
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.get_nondeg_coupled_space" class="docs-object-method">&nbsp;</a> 
```python
get_nondeg_coupled_space(self, input_state_space, degenerate_space=None, spaces=None, wavefunction_terms=None, property_filter=None, filter_spaces=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L1538)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L1538?message=Update%20Docs)]
</div>
Applies the non-degenerate equations in semi-symbolic form to determine
which states needs to be calculated.
This will always be the initial input to a calculation and then
certain additional states may be calculated on the fly if they are needed to handle
truly degenerate stuff.
The only difference there will be to add on
  - `input_state_space`: `BasisStateSpace`
    > 
  - `degenerate_space`: `BasisStateSpace`
    > 
  - `spaces`: `Any`
    > 
  - `wavefunction_terms`: `None | Iterable[Iterable[int]]`
    > which terms to include when calculating corrections
  - `property_filter`: `None | Iterable[Iterable[int], Iterable[Iterable[int]]]`
    > a set of states and selection rules to allow for being targeted in state to calculate
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.get_corrections" class="docs-object-method">&nbsp;</a> 
```python
get_corrections(self, non_zero_cutoff=None, handle_strong_couplings=None, check_overlap=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L1649)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L1649?message=Update%20Docs)]
</div>
Applies the perturbation theory equations to obtain
corrections to the wave functions and energies
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.high_frequency_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
high_frequency_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L1697)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L1697?message=Update%20Docs)]
</div>
**LLM Docstring**

The indices of the vibrational modes whose fundamental transition frequency exceeds `self.low_frequency_mode_cutoff`, used as the default set of modes considered for strong-coupling/degeneracy testing.
  - `:returns`: `list[int]`
    > the high-frequency mode indices


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.identify_strong_couplings" class="docs-object-method">&nbsp;</a> 
```python
identify_strong_couplings(self, corrs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L1716)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L1716?message=Update%20Docs)]
</div>
**LLM Docstring**

Find the strongly-coupled state pairs among the solved corrections, delegating to whichever configured `degeneracy_handlers` spec exposes a `wfc_threshold` (or, if none do, a freshly built default `'auto'` spec).
  - `corrs`: `PerturbationTheoryCorrections`
    > the perturbation-theory corrections to search for strong couplings in
  - `:returns`: `tuple`
    > `(strong_couplings, threshold)` -- the identified strong-coupling data and the threshold used to find it


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.construct_strong_coupling_spaces" class="docs-object-method">&nbsp;</a> 
```python
construct_strong_coupling_spaces(self, spec, sc, corrs, states, threshold): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L1740)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L1740?message=Update%20Docs)]
</div>
**LLM Docstring**

Build the degenerate state groups implied by a set of identified strong couplings, and, if the spec allows extending the state space (`spec.extend_spaces`), extend the solver's states/representations to cover any newly implicated states before returning the resulting (possibly larger) state space and representation data.
  - `spec`: `DegeneracySpec`
    > the degeneracy spec (typically a `StronglyCoupledDegeneracySpec`) whose group filter/extension settings govern this construction
  - `sc`: `dict`
    > the raw strong-coupling data (per-state, per-order coupled-state lists) to build groups from
  - `corrs`: `PerturbationTheoryCorrections`
    > the perturbation-theory corrections the coupling data came from, used to collapse it into per-state coupled-state spaces
  - `states`: `BasisStateSpace`
    > the current target state space
  - `threshold`: `float`
    > the coupling-strength threshold to build the group filter with
  - `:returns`: `tuple`
    > `(degenerate_states, (states, perturbations, flat_total_space, N))` -- the identified degenerate groups, and the (possibly extended) states/representations/flat space/dimension to use going forward


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.apply_VPT_equations" class="docs-object-method">&nbsp;</a> 
```python
apply_VPT_equations(self, state_index, degenerate_space_indices, degenerate_energies, zero_order_state, degenerate_subspace, degenerate_subsubspace, perturbations=None, allow_PT_degs=None, ignore_odd_orders=None, intermediate_normalization=None, non_zero_cutoff=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L2128)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L2128?message=Update%20Docs)]
</div>
Applies VPT equations, dispatching based on how many
degeneracies we need to handle
  - `state_index`: `int`
    > the index of the primary state being treated using the PT
  - `degenerate_space_indices`: `np.ndarray[int]`
    > the indices corresponding to degeneracies with the primary state in the zero-order picture
  - `degenerate_energies`: `Iterable[float | None]`
    > the first and (possibly) second order correction to the energies
  - `zero_order_states`: `np.ndarray[float]`
    > the vector for the proper zero-order state corresponding ot state_index
  - `degenerate_subsubspace`: `tuple[np.ndarray[float], np.ndarray[int]]`
    > the set of vectors for the zero-order states in the secondary degenerate subspace
  - `non_zero_cutoff`: `float`
    > cutoff for when a term can be called zero for performance reasons
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.apply_VPT_nondeg_equations" class="docs-object-method">&nbsp;</a> 
```python
apply_VPT_nondeg_equations(self, state_index, deg_group, perturbations=None, non_zero_cutoff=None, check_overlap=None, intermediate_normalization=False, ignore_odd_orders=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L2175)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L2175?message=Update%20Docs)]
</div>
Does the dirty work of doing the VPT iterative equations.
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.apply_VPT_2k1_rules" class="docs-object-method">&nbsp;</a> 
```python
apply_VPT_2k1_rules(self, existing_corrs, perturbations=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L2333)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L2333?message=Update%20Docs)]
</div>
Apply expressions allowing for obtaining higher-order
corrections to the energies from lower-order corrections to the
wavefunctions
  - `existing_corrs`: `Any`
    > 
  - `perturbations`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.apply_post_PT_variational_calc" class="docs-object-method">&nbsp;</a> 
```python
apply_post_PT_variational_calc(self, degenerate_states, corrs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L2407)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L2407?message=Update%20Docs)]
</div>
Applies degenerate perturbation theory by building a representation
for the degenerate terms in the Hamiltonian.
This is then diagonalized, allowing the degenerate states to be expressed
in the basis of non-degenerate states
  - `H`: `Iterable[SparseArray]`
    > 
  - `corrs`: `PerturbationTheoryCorrections`
    > the standard PerturbationTheory Corrections object that comes out of the application of non-deg PT
  - `degenerate_states`: `Any`
    > population of degenerate states
  - `logger`: `Logger`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.drop_deg_pert_els" class="docs-object-method">&nbsp;</a> 
```python
drop_deg_pert_els(self, perts, deg_groups): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L2509)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L2509?message=Update%20Docs)]
</div>

  - `perts`: `Any`
    > 
  - `deg_groups`: `Any`
    > 
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Solver/PerturbationTheorySolver.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Solver/PerturbationTheorySolver.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Solver/PerturbationTheorySolver.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Solver/PerturbationTheorySolver.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver.py#L27?message=Update%20Docs)   
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