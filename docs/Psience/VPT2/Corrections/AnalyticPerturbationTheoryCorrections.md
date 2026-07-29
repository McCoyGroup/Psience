## <a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections">AnalyticPerturbationTheoryCorrections</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections.py#L1046)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections.py#L1046?message=Update%20Docs)]
</div>

AnalyticPerturbationTheoryCorrections(states: Psience.BasisReps.StateSpaces.BasisStateSpace, state_lists: 'list[tuple[np.ndarray, np.ndarray]]', _energies: numpy.ndarray = None, _transition_moments: 'Iterable[np.ndarray]' = None, _spectra: 'Iterable[DiscreteSpectrum]' = None, _deperturbed_energies: numpy.ndarray = None, _deperturbed_transition_moments: 'Iterable[np.ndarray]' = None, _deperturbed_spectra: Psience.Spectra.BaseSpectrum.DiscreteSpectrum = None, degenerate_states: 'Iterable[BasisStateSpace]' = None, only_degenerate_terms: 'bool' = True, _degenerate_hamiltonians: 'Iterable[np.ndarray]' = None, _degenerate_coefficients: 'Iterable[np.ndarray]' = None, _degenerate_state_list_transformations: 'Iterable[list[np.ndarray, np.ndarray]]' = None, energy_corrections: Psience.VPT2.Corrections.PTCorrections = None, transition_moment_corrections: 'Iterable[BasicAPTCorrections]' = None, degenerate_hamiltonian_corrections: 'Iterable[BasicAPTCorrections]' = None, operator_corrections: 'Iterable[BasicAPTCorrections]' = None, _deperturbed_operator_values: 'Iterable[np.ndarray]' = None, _operator_values: 'Iterable[np.ndarray]' = None, operator_keys: 'Iterable[Any]' = None, logger: 'Logger' = None, _zpe_pos: int = None)







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
degenerate_states: NoneType
only_degenerate_terms: bool
energy_corrections: NoneType
transition_moment_corrections: NoneType
degenerate_hamiltonian_corrections: NoneType
operator_corrections: NoneType
operator_keys: NoneType
logger: NoneType
```
<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.get_zpe_pos" class="docs-object-method">&nbsp;</a> 
```python
get_zpe_pos(self) -> int: 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1072)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1072?message=Update%20Docs)]
</div>
**LLM Docstring**

Find (and cache) the index of the zero-point-energy (ground) state within `self.states`, falling back to index `0` if the all-zeros excitation vector can't be found.
  - `:returns`: `int`
    > the ZPE state's index


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.energies" class="docs-object-method">&nbsp;</a> 
```python
@property
energies(self) -> numpy.ndarray: 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1092)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1092?message=Update%20Docs)]
</div>
**LLM Docstring**

The final state energies: the sum of the (per-order) energy corrections if there are no degenerate states, otherwise the degenerate-perturbation-theory-corrected energies (computed via `get_degenerate_transformations`, which also populates the cached degenerate Hamiltonians/coefficients).
  - `:returns`: `np.ndarray`
    > the state energies


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.deperturbed_energies" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_energies(self) -> numpy.ndarray: 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1112)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1112?message=Update%20Docs)]
</div>
**LLM Docstring**

The (cached) deperturbed state energies -- the sum of the per-order energy corrections, without any degenerate-perturbation-theory rotation applied.
  - `:returns`: `np.ndarray`
    > the deperturbed energies


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.handle_degenerate_transformation" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
handle_degenerate_transformation(cls, degenerate_ham): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1126)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1126?message=Update%20Docs)]
</div>
**LLM Docstring**

Diagonalize a degenerate-block Hamiltonian and reorder its eigenvalues/eigenvectors so that each output state is matched (by maximum overlap) to a distinct input state, avoiding two output states mapping to the same input state.
  - `degenerate_ham`: `np.ndarray`
    > the (Hermitian) degenerate-block Hamiltonian matrix to diagonalize
  - `:returns`: `tuple[np.ndarray, np.ndarray]`
    > `(deg_engs, deg_transf)` -- the reordered eigenvalues and eigenvector (transformation) matrix


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.get_degenerate_transformations" class="docs-object-method">&nbsp;</a> 
```python
get_degenerate_transformations(self, basis, energies): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1167)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1167?message=Update%20Docs)]
</div>
**LLM Docstring**

Apply degenerate perturbation theory block by block: for each group of degenerate states, build its Hamiltonian block (summing the per-order corrections, with the diagonal set to the current deperturbed energies), diagonalize it via `handle_degenerate_transformation`, and write the resulting energies back into the running energy array.
  - `basis`: `BasisStateSpace`
    > the state space the energies are indexed against
  - `energies`: `np.ndarray`
    > the deperturbed energies to correct, indexed against `basis`
  - `:returns`: `tuple[np.ndarray, tuple[list[np.ndarray], list[np.ndarray]]]`
    > `(energies, (hams, transf))` -- the corrected energies, and the list of per-block Hamiltonians and diagonalizing transformations


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.degenerate_hamiltonians" class="docs-object-method">&nbsp;</a> 
```python
@property
degenerate_hamiltonians(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1218)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1218?message=Update%20Docs)]
</div>
**LLM Docstring**

The (cached) per-block degenerate Hamiltonians, computed as a side effect of evaluating `energies` if not already cached.
  - `:returns`: `list[np.ndarray]`
    > the list of degenerate-block Hamiltonians


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.degenerate_coefficients" class="docs-object-method">&nbsp;</a> 
```python
@property
degenerate_coefficients(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1231)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1231?message=Update%20Docs)]
</div>
**LLM Docstring**

The (cached) per-block degenerate-perturbation-theory mixing coefficients, computed as a side effect of evaluating `energies` if not already cached.
  - `:returns`: `list[np.ndarray]`
    > the list of degenerate-block mixing-coefficient matrices


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.get_freqs" class="docs-object-method">&nbsp;</a> 
```python
get_freqs(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1245)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1245?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the vibrational transition frequencies (final energies) relative to the zero-point energy.
  - `:returns`: `np.ndarray`
    > the frequencies


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.get_deperturbed_freqs" class="docs-object-method">&nbsp;</a> 
```python
get_deperturbed_freqs(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1256)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1256?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the deperturbed vibrational transition frequencies relative to the deperturbed zero-point energy, if degenerate states are present; otherwise falls back to `get_freqs`.
  - `:returns`: `np.ndarray`
    > the deperturbed frequencies


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.degenerate_transformation_pairs" class="docs-object-method">&nbsp;</a> 
```python
@property
degenerate_transformation_pairs(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1270)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1270?message=Update%20Docs)]
</div>
**LLM Docstring**

The (cached) per-block `(row_tf, col_tf)` transformation pairs mapping each `state_lists` block's raw initial/final states onto their degenerate-perturbation-theory-mixed counterparts, computed lazily via `_get_degenerate_tfs_mats`.
  - `:returns`: `list[list[np.ndarray]]`
    > the list of per-block `[row_tf, col_tf]` transformation pairs


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.transition_moments" class="docs-object-method">&nbsp;</a> 
```python
@property
transition_moments(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1433)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1433?message=Update%20Docs)]
</div>
**LLM Docstring**

The (cached) final transition-dipole moments: the deperturbed transition moments if there are no degenerate states, otherwise those moments rotated by the degenerate-perturbation-theory transformation (via `_apply_degs_to_corrs`).
  - `:returns`: `list`
    > the per-axis, per-block transition moments


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.harmonic_transition_moments" class="docs-object-method">&nbsp;</a> 
```python
@property
harmonic_transition_moments(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1458)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1458?message=Update%20Docs)]
</div>
**LLM Docstring**

The purely harmonic (zeroth-order) transition-dipole moments, extracted from the first entry of each block's transition-moment corrections.
  - `:returns`: `list`
    > the per-axis, per-block harmonic transition moments


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.deperturbed_transition_moments" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_transition_moments(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1475)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1475?message=Update%20Docs)]
</div>
**LLM Docstring**

The (cached) deperturbed transition-dipole moments, summing each block's transition-moment corrections over all perturbative orders.
  - `:returns`: `list`
    > the per-axis, per-block deperturbed transition moments


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.get_spectra" class="docs-object-method">&nbsp;</a> 
```python
get_spectra(self, energies, transition_moments): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1494)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1494?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a per-block list of discrete IR spectra from a set of state energies and transition moments, one spectrum per initial state in each `state_lists` block, with transition frequencies computed relative to that initial state's energy.
  - `energies`: `np.ndarray`
    > the state energies to compute transition frequencies from
  - `transition_moments`: `list`
    > the per-axis, per-block transition-moment data (as returned by `transition_moments`/`deperturbed_transition_moments`/`harmonic_transition_moments`)
  - `:returns`: `list[list[DiscreteSpectrum]]`
    > the per-block lists of per-initial-state discrete spectra


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.harmonic_spectra" class="docs-object-method">&nbsp;</a> 
```python
@property
harmonic_spectra(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1522)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1522?message=Update%20Docs)]
</div>
**LLM Docstring**

The purely harmonic (zeroth-order) IR spectra, built from the zeroth-order energy corrections and harmonic transition moments via `get_spectra`.
  - `:returns`: `list[list[DiscreteSpectrum]]`
    > the per-block lists of per-initial-state harmonic spectra


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.deperturbed_spectra" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_spectra(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1537)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1537?message=Update%20Docs)]
</div>
**LLM Docstring**

The (cached) deperturbed IR spectra, built from the deperturbed energies and transition moments via `get_spectra`.
  - `:returns`: `list[list[DiscreteSpectrum]]`
    > the per-block lists of per-initial-state deperturbed spectra


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.spectra" class="docs-object-method">&nbsp;</a> 
```python
@property
spectra(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1554)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1554?message=Update%20Docs)]
</div>
**LLM Docstring**

The (cached) final IR spectra: the deperturbed spectra if there are no degenerate states, otherwise the spectra built from the degenerate-perturbation-theory-corrected energies and transition moments.
  - `:returns`: `list[list[DiscreteSpectrum]]`
    > the per-block lists of per-initial-state spectra


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.deperturbed_operator_values" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_operator_values(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1574)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1574?message=Update%20Docs)]
</div>
**LLM Docstring**

The (cached) deperturbed values of any extra tracked operators, summing each operator's per-block corrections over all perturbative orders.
  - `:returns`: `list`
    > the per-operator, per-block deperturbed operator values


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.operator_values" class="docs-object-method">&nbsp;</a> 
```python
@property
operator_values(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1591)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1591?message=Update%20Docs)]
</div>
**LLM Docstring**

The (cached) final values of any extra tracked operators: the deperturbed operator values if there are no degenerate states, otherwise those values rotated by the degenerate-perturbation-theory transformation (via `_apply_degs_to_corrs`).
  - `:returns`: `list`
    > the per-operator, per-block operator values


<a id="Psience.VPT2.Corrections.__create_fn__.<locals>.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, states: Psience.BasisReps.StateSpaces.BasisStateSpace, state_lists: 'list[tuple[np.ndarray, np.ndarray]]', _energies: numpy.ndarray = None, _transition_moments: 'Iterable[np.ndarray]' = None, _spectra: 'Iterable[DiscreteSpectrum]' = None, _deperturbed_energies: numpy.ndarray = None, _deperturbed_transition_moments: 'Iterable[np.ndarray]' = None, _deperturbed_spectra: Psience.Spectra.BaseSpectrum.DiscreteSpectrum = None, degenerate_states: 'Iterable[BasisStateSpace]' = None, only_degenerate_terms: 'bool' = True, _degenerate_hamiltonians: 'Iterable[np.ndarray]' = None, _degenerate_coefficients: 'Iterable[np.ndarray]' = None, _degenerate_state_list_transformations: 'Iterable[list[np.ndarray, np.ndarray]]' = None, energy_corrections: Psience.VPT2.Corrections.PTCorrections = None, transition_moment_corrections: 'Iterable[BasicAPTCorrections]' = None, degenerate_hamiltonian_corrections: 'Iterable[BasicAPTCorrections]' = None, operator_corrections: 'Iterable[BasicAPTCorrections]' = None, _deperturbed_operator_values: 'Iterable[np.ndarray]' = None, _operator_values: 'Iterable[np.ndarray]' = None, operator_keys: 'Iterable[Any]' = None, logger: 'Logger' = None, _zpe_pos: int = None) -> None: 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/__create_fn__.py#L)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/__create_fn__.py#L?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Corrections.__create_fn__.<locals>.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/__create_fn__/<locals>.py#L363)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/__create_fn__/<locals>.py#L363?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Corrections.__create_fn__.<locals>.__eq__" class="docs-object-method">&nbsp;</a> 
```python
__eq__(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/__create_fn__/<locals>.py#L)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/__create_fn__/<locals>.py#L?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections.py#L1046?message=Update%20Docs)   
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