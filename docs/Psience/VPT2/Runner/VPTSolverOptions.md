## <a id="Psience.VPT2.Runner.VPTSolverOptions">VPTSolverOptions</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L1010)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L1010?message=Update%20Docs)]
</div>

Provides a helper to keep track of the options available
for configuring the way the perturbation theory is applied







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.VPT2.Runner.VPTSolverOptions.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, order=2, expansion_order=None, coupled_states=None, total_space=None, flat_total_space=None, state_space_iterations=None, state_space_terms=None, state_space_filters=None, extended_state_space_filter_generator=None, extended_state_space_postprocessor=None, allow_post_PT_calc=None, modify_degenerate_perturbations=None, gaussian_resonance_handling=None, ignore_odd_order_energies=None, intermediate_normalization=None, zero_element_warning=None, degenerate_states=None, handle_strong_couplings=None, strong_coupling_test_modes=None, strong_couplings_state_filter=None, strongly_coupled_group_filter=None, extend_strong_coupling_spaces=None, strong_coupling_zero_order_energy_cutoff=None, low_frequency_mode_cutoff=None, zero_order_energy_corrections=None, check_overlap=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L1070)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L1070?message=Update%20Docs)]
</div>

  - `order`: `int`
    > the order of perturbation theory to apply
  - `expansion_order`: `int | dict`
    > the order to go to in the expansions of the perturbations, this can be supplied for different properties independently, like
```python
expansion_order = {
    'potential':some_int,
    'kinetic':some_int,
    'dipole':some_int
}
```
  - `degenerate_states`: `Iterable[BasisStateSpace]`
    > the set of degeneracies to handle
  - `coupled_states`: `Iterable[SelectionRuleStateSpace]`
    > explicit bases of states to use at each order in the perturbation theory
  - `total_space`: `Iterable[BasisStateSpace]`
    > the total state spaces at each order in the perturbation theory
  - `flat_total_space`: `BasisStateSpace`
    > the union of all of the total state spaces
  - `state_space_iterations`: `int`
    > the order to go to when getting the `coupled_states`
  - `state_space_terms`: `Iterable[(int, int)]`
    > the explicit set of terms to include, as a tuple `(i, j)` which indicates `(H(i), |n(j)>)`
  - `state_space_filters`: `dict`
    > filters that can be used to cut down on the size of bases (see `VPTRunner.get_state_space_filter`)
  - `allow_post_PT_calc`: `bool`
    > whether to do the post-perturbation theory variational calculation for degeneracy handling
  - `modify_degenerate_perturbations`: `bool`
    > whether to modify the perturbation representation matrices themselves when doing degeneracy handling
  - `gaussian_resonance_handling`: `bool`
    > whether or not to skip the post-PT variational calculation for states with more than two quanta of excitation
  - `ignore_odd_order_energies`: `bool`
    > whether or not to skip actually calculating the energy corrections for odd orders
  - `intermediate_normalization`: `bool`
    > whether or not to use 'intermediate normalization' in the wavefunctions
  - `zero_element_warning`: `bool`
    > whether or not to warn if an element of the representations evaluated to zero (i.e. we wasted effort)
  - `low_frequency_mode_cutoff`: `float (default:500 cm-1)`
    > the energy below which to consider a mode to be "low frequency"
  - `zero_order_energy_corrections`: `dict`
    > energies to use for the zero-order states instead of the diagonal of `H(0)`
  - `check_overlap`: `bool default:True`
    > whether or not to ensure states are normalized in the VPT


<a id="Psience.VPT2.Runner.VPTSolverOptions.get_zero_order_energies" class="docs-object-method">&nbsp;</a> 
```python
@staticmethod
get_zero_order_energies(corrected_fundamental_freqs, states): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/staticmethod.py#L1199)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/staticmethod.py#L1199?message=Update%20Docs)]
</div>

  - `corrected_fundamental_freqs`: `Any`
    > 
  - `states`: `Any`
    > 
  - `:returns`: `_`
    >
 </div>
</div>



<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#Details-a9d7f5" markdown="1"> Details</a> <a class="float-right" data-toggle="collapse" href="#Details-a9d7f5"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Details-a9d7f5" markdown="1">
 The `basis_postfilters` have multiple possible values.
Here are the currently supported cases

```python
{
    'max_quanta': [2, -1, 1, -1, ...] # the max number of quanta allowed in a given mode in the basis (-1 means infinity)
}
```

- for excluding transitions

```python
{
    'excluded_transitions': [[0, 0, 1, 0, ...], [1, 0, 0, 0, ...], ...] # a set of transitions that are forbidden on the input states
}
```

- for excluding based on a test

```python
{
    'test': func # a function that takes the basis and tests if states should be allowed
}
```
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Runner/VPTSolverOptions.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Runner/VPTSolverOptions.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Runner/VPTSolverOptions.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Runner/VPTSolverOptions.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L1010?message=Update%20Docs)   
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