## <a id="Psience.Spectra.BaseSpectrum.DiscreteSpectrum">DiscreteSpectrum</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/BaseSpectrum.py#L152)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum.py#L152?message=Update%20Docs)]
</div>

Concrete implementation of `BaseSpectrum` that exists
solely to allow for plotting and broadening.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.Spectra.BaseSpectrum.DiscreteSpectrum.from_transition_moments" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_transition_moments(cls, frequencies, transition_moments, **meta): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L158)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L158?message=Update%20Docs)]
</div>
Assumes frequencies and transition moments in a.u.
  - `frequencies`: `Any`
    > 
  - `transition_moments`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Spectra.BaseSpectrum.DiscreteSpectrum.from_raman_moments" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_raman_moments(cls, frequencies, transition_polarizabilities, pump_frequency=0, **meta): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L180)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L180?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a discrete Raman spectrum from transition frequencies and polarizability transition moments (in atomic units), converting frequencies to wavenumbers and computing intensities from the transition-polarizability magnitudes scaled by the (pump-shifted) frequency to the fourth power.
  - `frequencies`: `np.ndarray`
    > the transition frequencies, in Hartrees
  - `transition_polarizabilities`: `np.ndarray`
    > the polarizability transition-moment tensors
  - `pump_frequency`: `float`
    > the pump laser frequency to add before the frequency^4 scaling
  - `meta`: `dict`
    > extra metadata stored on the spectrum
  - `:returns`: `DiscreteSpectrum`
    > the constructed Raman spectrum


<a id="Psience.Spectra.BaseSpectrum.DiscreteSpectrum.normalize" class="docs-object-method">&nbsp;</a> 
```python
normalize(self, which=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/BaseSpectrum/DiscreteSpectrum.py#L209)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum/DiscreteSpectrum.py#L209?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a copy of the spectrum with intensities normalized by either the overall maximum intensity or a specific reference intensity.
  - `which`: `int | None`
    > the index of a specific transition to normalize against; if `None`, normalizes by the maximum intensity
  - `:returns`: `DiscreteSpectrum`
    > the normalized spectrum


<a id="Psience.Spectra.BaseSpectrum.DiscreteSpectrum.plot" class="docs-object-method">&nbsp;</a> 
```python
plot(self, figure=None, plot_style=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/BaseSpectrum/DiscreteSpectrum.py#L226)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum/DiscreteSpectrum.py#L226?message=Update%20Docs)]
</div>
Plots a spectrum using `McUtils.Plots.StickSpectrum`
  - `figure`: `None | McUtils.Plots.Graphics`
    > figure to plot the spectrum on
  - `opts`: `Any`
    > any of the many, many options supported by `McUtils.Plots.Graphics`
  - `:returns`: `_`
    >


<a id="Psience.Spectra.BaseSpectrum.DiscreteSpectrum.broaden" class="docs-object-method">&nbsp;</a> 
```python
broaden(self, breadth=10, *, broadening_type='gaussian', noising_function=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/BaseSpectrum/DiscreteSpectrum.py#L249)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum/DiscreteSpectrum.py#L249?message=Update%20Docs)]
</div>
Applies a broadening to the spectrum
  - `broadening_type`: `Any`
    > 
  - `breadth`: `Any`
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Spectra/BaseSpectrum/DiscreteSpectrum.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Spectra/BaseSpectrum/DiscreteSpectrum.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Spectra/BaseSpectrum/DiscreteSpectrum.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Spectra/BaseSpectrum/DiscreteSpectrum.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum.py#L152?message=Update%20Docs)   
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