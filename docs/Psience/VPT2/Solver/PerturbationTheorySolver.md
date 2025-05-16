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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L260)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L260?message=Update%20Docs)]
</div>

  - `:returns`: `Iterable[SparseArray]`
    >


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.degenerate_spaces" class="docs-object-method">&nbsp;</a> 
```python
@property
degenerate_spaces(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L273)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L273?message=Update%20Docs)]
</div>

  - `:returns`: `_`
    >


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.merge_deg_spaces" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
merge_deg_spaces(cls, new_states): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L320)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L320?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.zero_order_energies" class="docs-object-method">&nbsp;</a> 
```python
@property
zero_order_energies(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L342)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L342?message=Update%20Docs)]
</div>

  - `:returns`: `_`
    >


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.apply_VPT" class="docs-object-method">&nbsp;</a> 
```python
apply_VPT(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L386)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L386?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L417)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L417?message=Update%20Docs)]
</div>
Gets the sparse representations of the passed perturbation inside the basis of coupled states.
  - `:returns`: `Iterable[SparseArray]`
    >


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.extend_VPT_representations" class="docs-object-method">&nbsp;</a> 
```python
extend_VPT_representations(self, new_flat_space, new_states): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L493)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L493?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.load_state_spaces" class="docs-object-method">&nbsp;</a> 
```python
load_state_spaces(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L665)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L665?message=Update%20Docs)]
</div>

  - `:returns`: `_`
    >


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.extend_state_spaces" class="docs-object-method">&nbsp;</a> 
```python
extend_state_spaces(self, new_targets, degenerate_states=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L793)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L793?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.load_coupled_spaces" class="docs-object-method">&nbsp;</a> 
```python
load_coupled_spaces(self, degenerate_spaces=None, spaces=None, wavefunction_terms=None, property_filter=None, filter_spaces=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L862)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L862?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L1278)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L1278?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L1311)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L1311?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L1422)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L1422?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L1470)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L1470?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.identify_strong_couplings" class="docs-object-method">&nbsp;</a> 
```python
identify_strong_couplings(self, corrs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L1481)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L1481?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.construct_strong_coupling_spaces" class="docs-object-method">&nbsp;</a> 
```python
construct_strong_coupling_spaces(self, spec, sc, corrs, states, threshold): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L1495)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L1495?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.apply_VPT_equations" class="docs-object-method">&nbsp;</a> 
```python
apply_VPT_equations(self, state_index, degenerate_space_indices, degenerate_energies, zero_order_state, degenerate_subspace, degenerate_subsubspace, perturbations=None, allow_PT_degs=None, ignore_odd_orders=None, intermediate_normalization=None, non_zero_cutoff=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L1853)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L1853?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L1900)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L1900?message=Update%20Docs)]
</div>
Does the dirty work of doing the VPT iterative equations.
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Solver.PerturbationTheorySolver.apply_VPT_2k1_rules" class="docs-object-method">&nbsp;</a> 
```python
apply_VPT_2k1_rules(self, existing_corrs, perturbations=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L2058)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L2058?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L2108)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L2108?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L2210)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Solver/PerturbationTheorySolver.py#L2210?message=Update%20Docs)]
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