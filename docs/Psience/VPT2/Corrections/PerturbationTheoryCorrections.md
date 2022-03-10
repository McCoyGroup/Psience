## <a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections">PerturbationTheoryCorrections</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Corrections.py#L15)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Corrections.py#L15?message=Update%20Docs)]
</div>

Represents a set of corrections from perturbation theory.
Can be used to correct other operators in the basis of the original calculation.

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, states, coupled_states, total_basis, energy_corrs, wfn_corrections, all_energy_corrections=None, degenerate_states=None, degenerate_transformation=None, degenerate_energies=None, logger=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Corrections.py#L21)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Corrections.py#L21?message=Update%20Docs)]
</div>


- `states`: `BasisStateSpace`
    >No description...
- `coupled_states`: `BasisMultiStateSpace`
    >No description...
- `total_basis`: `BasisMultiStateSpace`
    >No description...
- `energy_corrs`: `np.ndarray`
    >No description...
- `wfn_corrections`: `Iterable[SparseArray]`
    >No description...
- `degenerate_states`: `None | np.ndarray`
    >No description...
- `degenerate_transformation`: `None | np.ndarray`
    >No description...
- `degenerate_energies`: `None | np.ndarray`
    >No description...

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.from_dicts" class="docs-object-method">&nbsp;</a> 
```python
from_dicts(states, corrections, hamiltonians, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Corrections.py#L62)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Corrections.py#L62?message=Update%20Docs)]
</div>


- `states`: `dict`
    >a dict with the states described by the corrections, the set of states coupled, and the size of the overall basis
- `corrections`: `dict`
    >the corrections generated, including the corrections for the energies, wavefunctions, and a transformation from degenerate PT

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.degenerate" class="docs-object-method">&nbsp;</a> 
```python
@property
degenerate(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Corrections.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Corrections.py#L?message=Update%20Docs)]
</div>


- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.energies" class="docs-object-method">&nbsp;</a> 
```python
@property
energies(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Corrections.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Corrections.py#L?message=Update%20Docs)]
</div>


- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.order" class="docs-object-method">&nbsp;</a> 
```python
@property
order(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Corrections.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Corrections.py#L?message=Update%20Docs)]
</div>


- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.take_subspace" class="docs-object-method">&nbsp;</a> 
```python
take_subspace(self, space): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Corrections.py#L139)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Corrections.py#L139?message=Update%20Docs)]
</div>

Takes only those elements that are in space
- `space`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.create_coupling_matrix" class="docs-object-method">&nbsp;</a> 
```python
create_coupling_matrix(corrs, states: Psience.BasisReps.StateSpaces.BasisStateSpace, flat_total_space: Psience.BasisReps.StateSpaces.BasisStateSpace, nstates, order, filters=None, non_zero_cutoff=1e-14): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Corrections.py#L163)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Corrections.py#L163?message=Update%20Docs)]
</div>


- `corrs`: `Any`
    >No description...
- `states`: `Any`
    >No description...
- `flat_total_space`: `Any`
    >No description...
- `nstates`: `Any`
    >No description...
- `order`: `Any`
    >No description...
- `filters`: `Any`
    >No description...
- `non_zero_cutoff`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.prune" class="docs-object-method">&nbsp;</a> 
```python
prune(self, threshold=0.1): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Corrections.py#L296)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Corrections.py#L296?message=Update%20Docs)]
</div>

Returns corrections with couplings less than the given cutoff set to zero
- `threshold`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.find_strong_couplings" class="docs-object-method">&nbsp;</a> 
```python
find_strong_couplings(self, threshold=0.1): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Corrections.py#L308)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Corrections.py#L308?message=Update%20Docs)]
</div>

Finds positions in the expansion matrices where the couplings are too large
- `threshold`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.format_strong_couplings_report" class="docs-object-method">&nbsp;</a> 
```python
format_strong_couplings_report(self, couplings=None, threshold=0.1, int_fmt='{:>3.0f}', padding='{:<8}', join=True, use_excitations=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Corrections.py#L336)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Corrections.py#L336?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.collapse_strong_couplings" class="docs-object-method">&nbsp;</a> 
```python
collapse_strong_couplings(self, sc: dict): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Corrections.py#L356)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Corrections.py#L356?message=Update%20Docs)]
</div>


- `sc`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.operator_representation" class="docs-object-method">&nbsp;</a> 
```python
operator_representation(self, operator_expansion, order=None, subspace=None, contract=True, logger_symbol='A', logger_conversion=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Corrections.py#L423)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Corrections.py#L423?message=Update%20Docs)]
</div>

Generates the representation of the operator in the basis of stored states
- `operator_expansion`: `Iterable[float] | Iterable[np.ndarray]`
    >the expansion of the operator
- `order`: `Iterable[float] | Iterable[np.ndarray]`
    >the order of correction to go up to
- `subspace`: `None | BasisStateSpace`
    >the subspace of terms in which the operator expansion is defined
- `:returns`: `Iterable[np.ndarray]`
    >the set of representation matrices for this operator

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.get_overlap_matrices" class="docs-object-method">&nbsp;</a> 
```python
get_overlap_matrices(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Corrections.py#L504)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Corrections.py#L504?message=Update%20Docs)]
</div>

Returns the overlap matrices for the set of corrections
        at each order of correction
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.savez" class="docs-object-method">&nbsp;</a> 
```python
savez(self, file): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Corrections.py#L579)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Corrections.py#L579?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.loadz" class="docs-object-method">&nbsp;</a> 
```python
loadz(file): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Corrections.py#L595)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Corrections.py#L595?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.to_state" class="docs-object-method">&nbsp;</a> 
```python
to_state(self, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Corrections.py#L614)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Corrections.py#L614?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.from_state" class="docs-object-method">&nbsp;</a> 
```python
from_state(data, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Corrections.py#L626)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Corrections.py#L626?message=Update%20Docs)]
</div>



___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/ci/docs/Psience/VPT2/Corrections/PerturbationTheoryCorrections.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/ci/docs/Psience/VPT2/Corrections/PerturbationTheoryCorrections.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/ci/docs/Psience/VPT2/Corrections/PerturbationTheoryCorrections.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/ci/docs/Psience/VPT2/Corrections/PerturbationTheoryCorrections.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Corrections.py#L15?message=Update%20Docs)