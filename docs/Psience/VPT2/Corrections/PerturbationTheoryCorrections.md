## <a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections">PerturbationTheoryCorrections</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections.py#L18)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections.py#L18?message=Update%20Docs)]
</div>

Represents a set of corrections from perturbation theory.
Can be used to correct other operators in the basis of the original calculation.

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, states, coupled_states, total_basis, energy_corrs, wfn_corrections, all_energy_corrections=None, degenerate_states=None, degenerate_transformation=None, degenerate_energies=None, degenerate_hamiltonians=None, nondeg_hamiltonian_precision=3, logger=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections.py#L24)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections.py#L24?message=Update%20Docs)]
</div>


- `degenerate_energies`: `None | np.ndarray`
    >
- `degenerate_transformation`: `None | np.ndarray`
    >
- `degenerate_states`: `None | np.ndarray`
    >
- `wfn_corrections`: `Iterable[SparseArray]`
    >
- `energy_corrs`: `np.ndarray`
    >
- `total_basis`: `BasisMultiStateSpace`
    >
- `coupled_states`: `BasisMultiStateSpace`
    >
- `states`: `BasisStateSpace`
    >

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.from_dicts" class="docs-object-method">&nbsp;</a> 
```python
from_dicts(states, corrections, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections.py#L69)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections.py#L69?message=Update%20Docs)]
</div>


- `corrections`: `dict`
    >the corrections generated, including the corrections for the energies, wavefunctions, and a transformation from degenerate PT
- `states`: `dict`
    >a dict with the states described by the corrections, the set of states coupled, and the size of the overall basis

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.degenerate" class="docs-object-method">&nbsp;</a> 
```python
@property
degenerate(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections.py#L?message=Update%20Docs)]
</div>


- `:returns`: `_`
    >

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.energies" class="docs-object-method">&nbsp;</a> 
```python
@property
energies(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections.py#L?message=Update%20Docs)]
</div>


- `:returns`: `_`
    >

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.order" class="docs-object-method">&nbsp;</a> 
```python
@property
order(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections.py#L?message=Update%20Docs)]
</div>


- `:returns`: `_`
    >

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.take_subspace" class="docs-object-method">&nbsp;</a> 
```python
take_subspace(self, space): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections.py#L145)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections.py#L145?message=Update%20Docs)]
</div>

Takes only those elements that are in space
- `:returns`: `_`
    >
- `space`: `Any`
    >

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.create_coupling_matrix" class="docs-object-method">&nbsp;</a> 
```python
create_coupling_matrix(corrs, states: Psience.BasisReps.StateSpaces.BasisStateSpace, flat_total_space: Psience.BasisReps.StateSpaces.BasisStateSpace, nstates, order, filters=None, non_zero_cutoff=1e-14, logger=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections.py#L169)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections.py#L169?message=Update%20Docs)]
</div>


- `:returns`: `_`
    >
- `non_zero_cutoff`: `Any`
    >
- `filters`: `Any`
    >
- `order`: `Any`
    >
- `nstates`: `Any`
    >
- `flat_total_space`: `Any`
    >
- `states`: `Any`
    >
- `corrs`: `Any`
    >

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.prune" class="docs-object-method">&nbsp;</a> 
```python
prune(self, threshold=0.1, in_place=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections.py#L310)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections.py#L310?message=Update%20Docs)]
</div>

Returns corrections with couplings less than the given cutoff set to zero
- `:returns`: `_`
    >
- `threshold`: `Any`
    >

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.get_transformed_Hamiltonians" class="docs-object-method">&nbsp;</a> 
```python
get_transformed_Hamiltonians(self, hams, deg_group=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections.py#L349)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections.py#L349?message=Update%20Docs)]
</div>


- `:returns`: `_`
    >
- `deg_group`: `Any`
    >
- `corrs`: `Any`
    >

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.get_degenerate_rotation" class="docs-object-method">&nbsp;</a> 
```python
get_degenerate_rotation(self, deg_group, hams, zero_point_energy=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections.py#L379)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections.py#L379?message=Update%20Docs)]
</div>


- `:returns`: `_`
    >
- `corrs`: `Any`
    >
- `deg_group`: `Any`
    >

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.get_degenerate_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_degenerate_transformation(self, group, hams, gaussian_resonance_handling=False, zero_point_energy=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections.py#L483)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections.py#L483?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.default_state_filter" class="docs-object-method">&nbsp;</a> 
```python
default_state_filter(state, couplings, energy_cutoff=None, energies=None, basis=None, target_modes=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections.py#L513)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections.py#L513?message=Update%20Docs)]
</div>

Excludes modes that differ in only one position, prioritizing states with fewer numbers of quanta
(potentially add restrictions to high frequency modes...?)
- `:returns`: `_`
    >
- `couplings`: `Any`
    >
- `input_state`: `Any`
    >

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.find_strong_couplings" class="docs-object-method">&nbsp;</a> 
```python
find_strong_couplings(self, threshold=0.1, state_filter=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections.py#L557)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections.py#L557?message=Update%20Docs)]
</div>

Finds positions in the expansion matrices where the couplings are too large
- `:returns`: `_`
    >
- `threshold`: `Any`
    >

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.format_strong_couplings_report" class="docs-object-method">&nbsp;</a> 
```python
format_strong_couplings_report(self, couplings=None, threshold=0.1, int_fmt='{:>3.0f}', padding='{:<8}', join=True, use_excitations=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections.py#L590)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections.py#L590?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.collapse_strong_couplings" class="docs-object-method">&nbsp;</a> 
```python
collapse_strong_couplings(self, sc: dict): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections.py#L610)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections.py#L610?message=Update%20Docs)]
</div>


- `:returns`: `_`
    >
- `sc`: `Any`
    >

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.operator_representation" class="docs-object-method">&nbsp;</a> 
```python
operator_representation(self, operator_expansion, order=None, subspace=None, contract=True, logger_symbol='A', logger_conversion=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections.py#L680)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections.py#L680?message=Update%20Docs)]
</div>

Generates the representation of the operator in the basis of stored states
- `:returns`: `Iterable[np.ndarray]`
    >t
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
- `subspace`: `None | BasisStateSpace`
    >the subspace of terms in which the operator expansion is defined
- `order`: `Iterable[float] | Iterable[np.ndarray]`
    >the order of correction to go up to
- `operator_expansion`: `Iterable[float] | Iterable[np.ndarray]`
    >the expansion of the operator

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.get_overlap_matrices" class="docs-object-method">&nbsp;</a> 
```python
get_overlap_matrices(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections.py#L761)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections.py#L761?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections.py#L836)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections.py#L836?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.loadz" class="docs-object-method">&nbsp;</a> 
```python
loadz(file): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections.py#L852)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections.py#L852?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.to_state" class="docs-object-method">&nbsp;</a> 
```python
to_state(self, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections.py#L871)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections.py#L871?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.from_state" class="docs-object-method">&nbsp;</a> 
```python
from_state(data, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections.py#L883)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections.py#L883?message=Update%20Docs)]
</div>

 </div>
</div>






___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Corrections/PerturbationTheoryCorrections.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Corrections/PerturbationTheoryCorrections.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Corrections/PerturbationTheoryCorrections.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Corrections/PerturbationTheoryCorrections.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections.py#L18?message=Update%20Docs)