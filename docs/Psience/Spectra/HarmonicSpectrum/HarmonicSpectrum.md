## <a id="Psience.Spectra.HarmonicSpectrum.HarmonicSpectrum">HarmonicSpectrum</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/HarmonicSpectrum.py#L9)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/HarmonicSpectrum.py#L9?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.Spectra.HarmonicSpectrum.HarmonicSpectrum.from_normal_modes" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_normal_modes(cls, nms, dipole_derivatives, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L11)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L11?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a harmonic (double-harmonic-approximation) IR spectrum from a set of normal modes and their Cartesian dipole derivatives, converting the modes to a non-mass-weighted, frequency-scaled basis and computing the transition moments from the dipole derivatives projected onto the modes.
  - `nms`: `MixtureModes`
    > the normal modes to build the spectrum from
  - `dipole_derivatives`: `np.ndarray`
    > the Cartesian dipole-derivative tensor (first derivatives of the dipole with respect to Cartesian displacements)
  - `opts`: `dict`
    > extra options forwarded to `from_transition_moments`
  - `:returns`: `HarmonicSpectrum`
    > the constructed harmonic spectrum


<a id="Psience.Spectra.HarmonicSpectrum.HarmonicSpectrum.raman_from_modes" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
raman_from_modes(cls, nms, polarizability_derivatives, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L36)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L36?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a harmonic Raman spectrum from a set of normal modes and their Cartesian polarizability derivatives, converting the modes to a non-mass-weighted, frequency-scaled basis and computing the polarizability transition moments projected onto the modes.
  - `nms`: `MixtureModes`
    > the normal modes to build the spectrum from
  - `polarizability_derivatives`: `np.ndarray`
    > the Cartesian polarizability-derivative tensor (first derivatives of the polarizability with respect to Cartesian displacements)
  - `opts`: `dict`
    > extra options forwarded to `from_raman_moments`
  - `:returns`: `HarmonicSpectrum`
    > the constructed harmonic Raman spectrum


<a id="Psience.Spectra.HarmonicSpectrum.HarmonicSpectrum.from_mol" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_mol(cls, mol, modes=None, dipole_derivatives=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L61)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L61?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a harmonic IR spectrum for a molecule, using its normal modes and Cartesian dipole derivatives (computed from the molecule if not given directly), via `from_normal_modes`.
  - `mol`: `Molecule`
    > the molecule to build the spectrum for
  - `modes`: `MixtureModes | None`
    > the normal modes to use; computed via `mol.get_normal_modes()` if not given
  - `dipole_derivatives`: `object | None`
    > the dipole derivatives to use; if given as a full derivative expansion, the first-order term (`[1]`) is used, otherwise computed via `mol.get_cartesian_dipole_derivatives(1)[0]`
  - `opts`: `dict`
    > extra options forwarded to `from_normal_modes`
  - `:returns`: `HarmonicSpectrum`
    > the constructed harmonic spectrum


<a id="Psience.Spectra.HarmonicSpectrum.HarmonicSpectrum.raman_from_mol" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
raman_from_mol(cls, mol, modes=None, polarizability_derivatives=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L85)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L85?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a harmonic Raman spectrum for a molecule, using its normal modes and Cartesian polarizability derivatives (computed from the molecule if not given directly), via `raman_from_modes`.
  - `mol`: `Molecule`
    > the molecule to build the spectrum for
  - `modes`: `MixtureModes | None`
    > the normal modes to use; computed via `mol.get_normal_modes()` if not given
  - `polarizability_derivatives`: `list[np.ndarray] | None`
    > the full polarizability-derivative expansion to use (the first-order term is used); computed via `mol.get_cartesian_polarizability_derivatives(1)` if not given
  - `opts`: `dict`
    > extra options forwarded to `raman_from_modes`
  - `:returns`: `HarmonicSpectrum`
    > the constructed harmonic Raman spectrum
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Spectra/HarmonicSpectrum/HarmonicSpectrum.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Spectra/HarmonicSpectrum/HarmonicSpectrum.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Spectra/HarmonicSpectrum/HarmonicSpectrum.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Spectra/HarmonicSpectrum/HarmonicSpectrum.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/HarmonicSpectrum.py#L9?message=Update%20Docs)   
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