## <a id="Psience.Spectra.BaseSpectrum.BroadenedSpectrum">BroadenedSpectrum</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Spectra/BaseSpectrum.py#L222)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Spectra/BaseSpectrum.py#L222?message=Update%20Docs)]
</div>

A stick spectrum with associated broadening function







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.Spectra.BaseSpectrum.BroadenedSpectrum.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, frequencies, intensities, broadening_type='gaussian', breadth=10, **meta): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Spectra/BaseSpectrum/BroadenedSpectrum.py#L227)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Spectra/BaseSpectrum/BroadenedSpectrum.py#L227?message=Update%20Docs)]
</div>

  - `frequencies`: `Any`
    > 
  - `intensities`: `Any`
    > 
  - `broadening_type`: `"gaussian" | "lorentzian" | function`
    > the type of broadening to apply (can be any function)
  - `breadth`: `Any`
    > the breadth or list of breads for the peaks in the spectrum
  - `meta`: `Any`
    >


<a id="Psience.Spectra.BaseSpectrum.BroadenedSpectrum.plot" class="docs-object-method">&nbsp;</a> 
```python
plot(self, step_size=0.5, freq_min=None, freq_max=None, figure=None, plot_style=None, filled=False, adjust_width=True, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Spectra/BaseSpectrum/BroadenedSpectrum.py#L338)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Spectra/BaseSpectrum/BroadenedSpectrum.py#L338?message=Update%20Docs)]
</div>
Applies the broadening then plots it using `McUtils.Plots.Plot`
  - `step_size`: `Any`
    > step size to use when getting evaluation points for evaluating the broadening
  - `freq_min`: `Any`
    > min freq for evaluation
  - `freq_max`: `Any`
    > max freq for evaluation
  - `figure`: `Any`
    > 
  - `plot_style`: `Any`
    > 
  - `opts`: `Any`
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Spectra/BaseSpectrum/BroadenedSpectrum.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Spectra/BaseSpectrum/BroadenedSpectrum.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Spectra/BaseSpectrum/BroadenedSpectrum.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Spectra/BaseSpectrum/BroadenedSpectrum.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Spectra/BaseSpectrum.py#L222?message=Update%20Docs)   
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