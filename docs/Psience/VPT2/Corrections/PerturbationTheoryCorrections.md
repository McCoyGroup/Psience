## <a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections">PerturbationTheoryCorrections</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Corrections.py#L17)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Corrections.py#L17?message=Update%20Docs)]
</div>

Represents a set of corrections from perturbation theory.
Can be used to correct other operators in the basis of the original calculation.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, states, coupled_states, total_basis, energy_corrs, wfn_corrections, all_energy_corrections=None, degenerate_states=None, degenerate_transformation=None, degenerate_energies=None, degenerate_hamiltonians=None, nondeg_hamiltonian_precision=3, logger=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L23)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L23?message=Update%20Docs)]
</div>

  - `states`: `BasisStateSpace`
    > 
  - `coupled_states`: `BasisMultiStateSpace`
    > 
  - `total_basis`: `BasisMultiStateSpace`
    > 
  - `energy_corrs`: `np.ndarray`
    > 
  - `wfn_corrections`: `Iterable[SparseArray]`
    > 
  - `degenerate_states`: `None | np.ndarray`
    > 
  - `degenerate_transformation`: `None | np.ndarray`
    > 
  - `degenerate_energies`: `None | np.ndarray`
    >


<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.from_dicts" class="docs-object-method">&nbsp;</a> 
```python
from_dicts(states, corrections, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L68)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L68?message=Update%20Docs)]
</div>

  - `states`: `dict`
    > a dict with the states described by the corrections, the set of states coupled, and the size of the overall basis
  - `corrections`: `dict`
    > the corrections generated, including the corrections for the energies, wavefunctions, and a transformation from degenerate PT


<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.degenerate" class="docs-object-method">&nbsp;</a> 
```python
@property
degenerate(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L114)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L114?message=Update%20Docs)]
</div>

  - `:returns`: `_`
    >


<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.energies" class="docs-object-method">&nbsp;</a> 
```python
@property
energies(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L123)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L123?message=Update%20Docs)]
</div>

  - `:returns`: `_`
    >


<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.order" class="docs-object-method">&nbsp;</a> 
```python
@property
order(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L135)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L135?message=Update%20Docs)]
</div>

  - `:returns`: `_`
    >


<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.take_subspace" class="docs-object-method">&nbsp;</a> 
```python
take_subspace(self, space): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L144)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L144?message=Update%20Docs)]
</div>
Takes only those elements that are in space
  - `space`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.create_coupling_matrix" class="docs-object-method">&nbsp;</a> 
```python
create_coupling_matrix(corrs, states: Psience.BasisReps.StateSpaces.BasisStateSpace, flat_total_space: Psience.BasisReps.StateSpaces.BasisStateSpace, nstates, order, filters=None, non_zero_cutoff=1e-14, logger=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L168)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L168?message=Update%20Docs)]
</div>

  - `corrs`: `Any`
    > 
  - `states`: `Any`
    > 
  - `flat_total_space`: `Any`
    > 
  - `nstates`: `Any`
    > 
  - `order`: `Any`
    > 
  - `filters`: `Any`
    > 
  - `non_zero_cutoff`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.prune" class="docs-object-method">&nbsp;</a> 
```python
prune(self, threshold=0.1, in_place=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L309)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L309?message=Update%20Docs)]
</div>
Returns corrections with couplings less than the given cutoff set to zero
  - `threshold`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.get_transformed_Hamiltonians" class="docs-object-method">&nbsp;</a> 
```python
get_transformed_Hamiltonians(self, hams, deg_group=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L348)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L348?message=Update%20Docs)]
</div>

  - `corrs`: `Any`
    > 
  - `deg_group`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.get_degenerate_rotation" class="docs-object-method">&nbsp;</a> 
```python
get_degenerate_rotation(self, deg_group, hams, label=None, zero_point_energy=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L378)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L378?message=Update%20Docs)]
</div>

  - `deg_group`: `Any`
    > 
  - `corrs`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.get_degenerate_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_degenerate_transformation(self, group, hams, gaussian_resonance_handling=False, label=None, zero_point_energy=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L482)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L482?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.default_state_filter" class="docs-object-method">&nbsp;</a> 
```python
default_state_filter(state, couplings, energy_cutoff=None, energies=None, basis=None, target_modes=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L512)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L512?message=Update%20Docs)]
</div>
Excludes modes that differ in only one position, prioritizing states with fewer numbers of quanta
(potentially add restrictions to high frequency modes...?)
  - `input_state`: `Any`
    > 
  - `couplings`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.find_strong_couplings" class="docs-object-method">&nbsp;</a> 
```python
find_strong_couplings(self, threshold=0.1, state_filter=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L556)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L556?message=Update%20Docs)]
</div>
Finds positions in the expansion matrices where the couplings are too large
  - `threshold`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.format_strong_couplings_report" class="docs-object-method">&nbsp;</a> 
```python
format_strong_couplings_report(self, couplings=None, threshold=0.1, int_fmt='{:>3.0f}', padding='{:<8}', join=True, use_excitations=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L589)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L589?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.collapse_strong_couplings" class="docs-object-method">&nbsp;</a> 
```python
collapse_strong_couplings(self, sc: dict): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L609)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L609?message=Update%20Docs)]
</div>

  - `sc`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.operator_representation" class="docs-object-method">&nbsp;</a> 
```python
operator_representation(self, operator_expansion, order=None, subspace=None, contract=True, logger_symbol='A', logger_conversion=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L679)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L679?message=Update%20Docs)]
</div>
Generates the representation of the operator in the basis of stored states
  - `operator_expansion`: `Iterable[float] | Iterable[np.ndarray]`
    > the expansion of the operator
  - `order`: `Iterable[float] | Iterable[np.ndarray]`
    > the order of correction to go up to
  - `subspace`: `None | BasisStateSpace`
    > the subspace of terms in which the operator expansion is defined
  - `:returns`: `Iterable[np.ndarray]`
    > t
h
e
 
s
e
t
 
o
f
 
r
e
p
r
e
s
e
n
t
a
t
i
o
n
 
m
a
t
r
i
c
e
s
 
f
o
r
 
t
h
i
s
 
o
p
e
r
a
t
o
r


<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.get_overlap_matrices" class="docs-object-method">&nbsp;</a> 
```python
get_overlap_matrices(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L761)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L761?message=Update%20Docs)]
</div>
Returns the overlap matrices for the set of corrections
at each order of correction
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.savez" class="docs-object-method">&nbsp;</a> 
```python
savez(self, file): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L836)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L836?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.loadz" class="docs-object-method">&nbsp;</a> 
```python
loadz(file): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L852)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L852?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.to_state" class="docs-object-method">&nbsp;</a> 
```python
to_state(self, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L872)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L872?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.from_state" class="docs-object-method">&nbsp;</a> 
```python
from_state(data, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L884)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Corrections/PerturbationTheoryCorrections.py#L884?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Corrections/PerturbationTheoryCorrections.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Corrections/PerturbationTheoryCorrections.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Corrections/PerturbationTheoryCorrections.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Corrections/PerturbationTheoryCorrections.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Corrections.py#L17?message=Update%20Docs)   
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