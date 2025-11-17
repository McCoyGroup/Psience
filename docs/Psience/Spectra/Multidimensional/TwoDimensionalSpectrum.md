## <a id="Psience.Spectra.Multidimensional.TwoDimensionalSpectrum">TwoDimensionalSpectrum</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/Multidimensional.py#L10)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/Multidimensional.py#L10?message=Update%20Docs)]
</div>

Base class to support spectral operation







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
default_styles: dict
default_line_style: dict
```
<a id="Psience.Spectra.Multidimensional.TwoDimensionalSpectrum.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, freq1, freq2, intensities, **meta): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/Multidimensional.py#L14)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/Multidimensional.py#L14?message=Update%20Docs)]
</div>

  - `frequencies`: `np.ndarray`
    > frequency list
  - `intensities`: `np.ndarray`
    > intensity list
  - `meta`: `Any`
    > metadata


<a id="Psience.Spectra.Multidimensional.TwoDimensionalSpectrum.take_subspectrum" class="docs-object-method">&nbsp;</a> 
```python
take_subspectrum(self, sample_x, sample_y): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/Multidimensional/TwoDimensionalSpectrum.py#L28)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/Multidimensional/TwoDimensionalSpectrum.py#L28?message=Update%20Docs)]
</div>
Takes a subset of frequencies/intensities specified by `pos`
  - `pos`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Spectra.Multidimensional.TwoDimensionalSpectrum.frequency_filter" class="docs-object-method">&nbsp;</a> 
```python
frequency_filter(self, freq_span_x, freq_span_y): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/Multidimensional/TwoDimensionalSpectrum.py#L45)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/Multidimensional/TwoDimensionalSpectrum.py#L45?message=Update%20Docs)]
</div>


<a id="Psience.Spectra.Multidimensional.TwoDimensionalSpectrum.intensity_filter" class="docs-object-method">&nbsp;</a> 
```python
intensity_filter(self, int_min, int_max): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/Multidimensional/TwoDimensionalSpectrum.py#L53)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/Multidimensional/TwoDimensionalSpectrum.py#L53?message=Update%20Docs)]
</div>


<a id="Psience.Spectra.Multidimensional.TwoDimensionalSpectrum.clip" class="docs-object-method">&nbsp;</a> 
```python
clip(self, int_min, int_max, clip_abs=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/Multidimensional/TwoDimensionalSpectrum.py#L62)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/Multidimensional/TwoDimensionalSpectrum.py#L62?message=Update%20Docs)]
</div>


<a id="Psience.Spectra.Multidimensional.TwoDimensionalSpectrum.plot" class="docs-object-method">&nbsp;</a> 
```python
plot(self, plot_filled=True, contour_line_style=None, figure=None, symmetric_range=True, remove_baseline=True, vmin=None, vmax=None, levels=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/Multidimensional/TwoDimensionalSpectrum.py#L86)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/Multidimensional/TwoDimensionalSpectrum.py#L86?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Spectra/Multidimensional/TwoDimensionalSpectrum.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Spectra/Multidimensional/TwoDimensionalSpectrum.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Spectra/Multidimensional/TwoDimensionalSpectrum.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Spectra/Multidimensional/TwoDimensionalSpectrum.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/Multidimensional.py#L10?message=Update%20Docs)   
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