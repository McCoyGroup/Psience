## <a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections">AnalyticPerturbationTheoryCorrections</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections.py#L941)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections.py#L941?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L967)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L967?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.energies" class="docs-object-method">&nbsp;</a> 
```python
@property
energies(self) -> numpy.ndarray: 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L979)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L979?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.deperturbed_energies" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_energies(self) -> numpy.ndarray: 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L991)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L991?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.handle_degenerate_transformation" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
handle_degenerate_transformation(cls, degenerate_ham): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L997)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L997?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.get_degenerate_transformations" class="docs-object-method">&nbsp;</a> 
```python
get_degenerate_transformations(self, basis, energies): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1028)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1028?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.degenerate_hamiltonians" class="docs-object-method">&nbsp;</a> 
```python
@property
degenerate_hamiltonians(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1067)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1067?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.degenerate_coefficients" class="docs-object-method">&nbsp;</a> 
```python
@property
degenerate_coefficients(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1072)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1072?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.get_freqs" class="docs-object-method">&nbsp;</a> 
```python
get_freqs(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1078)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1078?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.get_deperturbed_freqs" class="docs-object-method">&nbsp;</a> 
```python
get_deperturbed_freqs(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1081)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1081?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.degenerate_transformation_pairs" class="docs-object-method">&nbsp;</a> 
```python
@property
degenerate_transformation_pairs(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1087)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1087?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.transition_moments" class="docs-object-method">&nbsp;</a> 
```python
@property
transition_moments(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1219)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1219?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.harmonic_transition_moments" class="docs-object-method">&nbsp;</a> 
```python
@property
harmonic_transition_moments(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1236)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1236?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.deperturbed_transition_moments" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_transition_moments(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1245)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1245?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.get_spectra" class="docs-object-method">&nbsp;</a> 
```python
get_spectra(self, energies, transition_moments): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1256)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1256?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.harmonic_spectra" class="docs-object-method">&nbsp;</a> 
```python
@property
harmonic_spectra(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1272)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1272?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.deperturbed_spectra" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_spectra(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1279)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1279?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.spectra" class="docs-object-method">&nbsp;</a> 
```python
@property
spectra(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1288)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1288?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.deperturbed_operator_values" class="docs-object-method">&nbsp;</a> 
```python
@property
deperturbed_operator_values(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1300)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1300?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Corrections.AnalyticPerturbationTheoryCorrections.operator_values" class="docs-object-method">&nbsp;</a> 
```python
@property
operator_values(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1309)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections/AnalyticPerturbationTheoryCorrections.py#L1309?message=Update%20Docs)]
</div>


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
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Corrections.py#L941?message=Update%20Docs)   
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