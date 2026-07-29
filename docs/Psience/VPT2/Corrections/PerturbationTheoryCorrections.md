## <a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections">PerturbationTheoryCorrections</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections.py#L23)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections.py#L23?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections.py#L29)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections.py#L29?message=Update%20Docs)]
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
@classmethod
from_dicts(cls, states, corrections, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L74)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L74?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/PerturbationTheoryCorrections.py#L120)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/PerturbationTheoryCorrections.py#L120?message=Update%20Docs)]
</div>

  - `:returns`: `_`
    >


<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.energies" class="docs-object-method">&nbsp;</a> 
```python
@property
energies(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/PerturbationTheoryCorrections.py#L129)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/PerturbationTheoryCorrections.py#L129?message=Update%20Docs)]
</div>

  - `:returns`: `_`
    >


<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.order" class="docs-object-method">&nbsp;</a> 
```python
@property
order(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/PerturbationTheoryCorrections.py#L141)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/PerturbationTheoryCorrections.py#L141?message=Update%20Docs)]
</div>

  - `:returns`: `_`
    >


<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.take_subspace" class="docs-object-method">&nbsp;</a> 
```python
take_subspace(self, space): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/PerturbationTheoryCorrections.py#L150)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/PerturbationTheoryCorrections.py#L150?message=Update%20Docs)]
</div>
Takes only those elements that are in space
  - `space`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.create_coupling_matrix" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
create_coupling_matrix(cls, corrs, states: Psience.BasisReps.StateSpaces.BasisStateSpace, flat_total_space: Psience.BasisReps.StateSpaces.BasisStateSpace, nstates, order, filters=None, non_zero_cutoff=1e-14, logger=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L174)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L174?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/PerturbationTheoryCorrections.py#L315)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/PerturbationTheoryCorrections.py#L315?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/PerturbationTheoryCorrections.py#L354)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/PerturbationTheoryCorrections.py#L354?message=Update%20Docs)]
</div>

  - `corrs`: `Any`
    > 
  - `deg_group`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.get_degenerate_rotation" class="docs-object-method">&nbsp;</a> 
```python
get_degenerate_rotation(self, deg_group, hams, label=None, zero_point_energy=None, local_coupling_hamiltonian=None, local_coupling_order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/PerturbationTheoryCorrections.py#L384)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/PerturbationTheoryCorrections.py#L384?message=Update%20Docs)]
</div>

  - `deg_group`: `Any`
    > 
  - `corrs`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.get_degenerate_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_degenerate_transformation(self, group, hams, gaussian_resonance_handling=False, label=None, zero_point_energy=None, local_coupling_hamiltonian=None, local_coupling_order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/PerturbationTheoryCorrections.py#L508)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/PerturbationTheoryCorrections.py#L508?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the degenerate-perturbation-theory rotation for a single group of resonant/degenerate states: finds the group's positions within this object's overall state space (warning about any states in the group that aren't present), and, if the group has more than one matched state (and isn't skipped due to Gaussian-style high-order resonance handling), diagonalizes the corresponding non-degenerate Hamiltonian block via `get_degenerate_rotation`.
  - `group`: `BasisStateSpace`
    > the group of mutually resonant/degenerate states to build the transformation for
  - `hams`: `list[np.ndarray]`
    > the Hamiltonian correction matrices to build the block from
  - `gaussian_resonance_handling`: `bool`
    > whether to skip building a rotation for groups whose states have more than 2 quanta of excitation (mimicking Gaussian's resonance-handling convention)
  - `label`: `str | None`
    > a label used for logging
  - `zero_point_energy`: `float | None`
    > the zero-point energy, used when building/logging the non-degenerate Hamiltonian block
  - `local_coupling_hamiltonian`: `np.ndarray | None`
    > an explicit local coupling Hamiltonian to use instead of building one from `hams`
  - `local_coupling_order`: `int | None`
    > the perturbative order to build the local coupling Hamiltonian to, if not given explicitly
  - `:returns`: `tuple`
    > `(deg_inds, H_nd, deg_rot, deg_engs)` -- the group's indices within this state space, the non-degenerate Hamiltonian block, the diagonalizing rotation, and the resulting (sorted) degenerate energies; `H_nd`/`deg_rot`/`deg_engs` are all `None` if the group doesn't need (or isn't eligible for) degenerate treatment


<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.default_state_filter" class="docs-object-method">&nbsp;</a> 
```python
@staticmethod
default_state_filter(state, couplings, energy_cutoff=None, energies=None, basis=None, target_modes=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/staticmethod.py#L568)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/staticmethod.py#L568?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/PerturbationTheoryCorrections.py#L612)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/PerturbationTheoryCorrections.py#L612?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/PerturbationTheoryCorrections.py#L645)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/PerturbationTheoryCorrections.py#L645?message=Update%20Docs)]
</div>
**LLM Docstring**

Format a human-readable report of the states found by `find_strong_couplings` (or an explicitly supplied `couplings` dict), listing each state alongside the other states it's strongly coupled to at each perturbative order.
  - `couplings`: `dict | None`
    > the strong-coupling data to format; computed via `find_strong_couplings` if not given
  - `threshold`: `float`
    > the coupling-strength threshold forwarded to `find_strong_couplings` if `couplings` isn't given
  - `int_fmt`: `str`
    > the format string used for each quantum-number column
  - `padding`: `str`
    > the format string used for row labels/indentation
  - `join`: `bool`
    > whether to join the report lines into a single string, or return them as a list
  - `use_excitations`: `bool`
    > whether to display states as their excitation-quantum-number vectors (rather than raw basis indices)
  - `:returns`: `str | list[str]`
    > the formatted report, as a single string or a list of lines


<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.collapse_strong_couplings" class="docs-object-method">&nbsp;</a> 
```python
collapse_strong_couplings(self, sc: dict): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/PerturbationTheoryCorrections.py#L685)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/PerturbationTheoryCorrections.py#L685?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/PerturbationTheoryCorrections.py#L774)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/PerturbationTheoryCorrections.py#L774?message=Update%20Docs)]
</div>
Generates the representation of the operator in the basis of stored states
  - `operator_expansion`: `Iterable[float] | Iterable[np.ndarray]`
    > the expansion of the operator
  - `order`: `Iterable[float] | Iterable[np.ndarray]`
    > the order of correction to go up to
  - `subspace`: `None | BasisStateSpace`
    > the subspace of terms in which the operator expansion is defined
  - `:returns`: `Iterable[np.ndarray]`
    > the set of representation matrices for this operator


<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.get_overlap_matrices" class="docs-object-method">&nbsp;</a> 
```python
get_overlap_matrices(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/PerturbationTheoryCorrections.py#L856)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/PerturbationTheoryCorrections.py#L856?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/PerturbationTheoryCorrections.py#L931)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/PerturbationTheoryCorrections.py#L931?message=Update%20Docs)]
</div>
**LLM Docstring**

Intended to serialize the corrections to an `npz` file, but currently disabled -- immediately raises, noting the implementation is outdated; use `to_state`/`from_state` instead.
  - `file`: `str`
    > the target file path
  - `:returns`: `None`
    > never returns


<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.loadz" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
loadz(cls, file): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L958)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L958?message=Update%20Docs)]
</div>
**LLM Docstring**

Intended to reconstruct corrections from an `npz` file previously written by `savez`, but currently disabled -- immediately raises, noting the implementation is outdated; use `to_state`/`from_state` instead.
  - `file`: `str`
    > the source file path
  - `:returns`: `PerturbationTheoryCorrections`
    > never returns


<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.to_state" class="docs-object-method">&nbsp;</a> 
```python
to_state(self, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/PerturbationTheoryCorrections.py#L989)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/PerturbationTheoryCorrections.py#L989?message=Update%20Docs)]
</div>
**LLM Docstring**

Serialize this object's core data (states, coupled states, total basis, energy/wavefunction corrections, and any degenerate-perturbation-theory data) into a plain dict.
  - `serializer`: `object | None`
    > accepted for interface consistency but not used in this method's body
  - `:returns`: `dict`
    > the serialized state dict


<a id="Psience.VPT2.Corrections.PerturbationTheoryCorrections.from_state" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_state(cls, data, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1011)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1011?message=Update%20Docs)]
</div>
**LLM Docstring**

Reconstruct a `PerturbationTheoryCorrections` object from a previously serialized state dict, deserializing the state-space objects via the given `serializer` and delegating to `from_dicts`.
  - `data`: `dict`
    > the serialized state, as produced by `to_state`
  - `serializer`: `object`
    > the serializer used to deserialize the state-space objects
  - `:returns`: `PerturbationTheoryCorrections`
    > the reconstructed corrections object
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
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections.py#L23?message=Update%20Docs)   
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