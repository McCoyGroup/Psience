## <a id="Psience.Spectra.BaseSpectrum.DiscreteSpectrum">DiscreteSpectrum</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Spectra/BaseSpectrum.py#L142)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Spectra/BaseSpectrum.py#L142?message=Update%20Docs)]
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
from_transition_moments(frequencies, transition_moments, **meta): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Spectra/BaseSpectrum/DiscreteSpectrum.py#L148)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Spectra/BaseSpectrum/DiscreteSpectrum.py#L148?message=Update%20Docs)]
</div>
Assumes frequencies and transition moments in a.u.
  - `frequencies`: `Any`
    > 
  - `transition_moments`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Spectra.BaseSpectrum.DiscreteSpectrum.normalize" class="docs-object-method">&nbsp;</a> 
```python
normalize(self, which=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Spectra/BaseSpectrum/DiscreteSpectrum.py#L170)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Spectra/BaseSpectrum/DiscreteSpectrum.py#L170?message=Update%20Docs)]
</div>


<a id="Psience.Spectra.BaseSpectrum.DiscreteSpectrum.plot" class="docs-object-method">&nbsp;</a> 
```python
plot(self, figure=None, plot_style=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Spectra/BaseSpectrum/DiscreteSpectrum.py#L177)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Spectra/BaseSpectrum/DiscreteSpectrum.py#L177?message=Update%20Docs)]
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
broaden(self, broadening_type='gaussian', breadth=10): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Spectra/BaseSpectrum/DiscreteSpectrum.py#L200)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Spectra/BaseSpectrum/DiscreteSpectrum.py#L200?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Spectra/BaseSpectrum.py#L142?message=Update%20Docs)   
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