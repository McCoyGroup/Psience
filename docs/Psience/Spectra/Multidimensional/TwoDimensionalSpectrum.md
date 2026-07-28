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
**LLM Docstring**

Restrict the spectrum to a rectangular window of frequency values along both axes, via `take_subspectrum`.
  - `freq_span_x`: `tuple[float, float]`
    > the `(min, max)` frequency bounds to keep along the `freq1` axis
  - `freq_span_y`: `tuple[float, float]`
    > the `(min, max)` frequency bounds to keep along the `freq2` axis
  - `:returns`: `TwoDimensionalSpectrum`
    > the frequency-restricted subspectrum


<a id="Psience.Spectra.Multidimensional.TwoDimensionalSpectrum.intensity_filter" class="docs-object-method">&nbsp;</a> 
```python
intensity_filter(self, int_min, int_max): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/Multidimensional/TwoDimensionalSpectrum.py#L65)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/Multidimensional/TwoDimensionalSpectrum.py#L65?message=Update%20Docs)]
</div>
**LLM Docstring**

Restrict the spectrum to the rectangular index range spanning the points whose intensity falls within `[int_min, int_max]`. Note the intensity test uses the bitwise `&` operator rather than `and`, which for numeric `int_min <= self.intensities` (a scalar comparison producing a bare bool) does not behave as intended for filtering an array by both bounds.
  - `int_min`: `float`
    > the minimum intensity to include
  - `int_max`: `float`
    > the maximum intensity to include
  - `:returns`: `TwoDimensionalSpectrum`
    > the intensity-restricted subspectrum


<a id="Psience.Spectra.Multidimensional.TwoDimensionalSpectrum.clip" class="docs-object-method">&nbsp;</a> 
```python
clip(self, int_min, int_max, clip_abs=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/Multidimensional/TwoDimensionalSpectrum.py#L86)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/Multidimensional/TwoDimensionalSpectrum.py#L86?message=Update%20Docs)]
</div>
**LLM Docstring**

Clip the spectrum's intensities to a `[int_min, int_max]` range, either clipping the absolute magnitude while preserving sign (and zeroing out points below `int_min` entirely) or clipping the raw signed values directly.
  - `int_min`: `float`
    > the minimum intensity magnitude (or raw value, if `clip_abs` is `False`) to keep
  - `int_max`: `float`
    > the maximum intensity magnitude (or raw value) to clip to
  - `clip_abs`: `bool`
    > whether to clip based on absolute magnitude (preserving sign) rather than the raw signed values
  - `:returns`: `TwoDimensionalSpectrum`
    > the clipped spectrum


<a id="Psience.Spectra.Multidimensional.TwoDimensionalSpectrum.plot" class="docs-object-method">&nbsp;</a> 
```python
plot(self, plot_filled=True, contour_line_style=None, figure=None, symmetric_range=True, remove_baseline=True, vmin=None, vmax=None, levels=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/Multidimensional/TwoDimensionalSpectrum.py#L124)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/Multidimensional/TwoDimensionalSpectrum.py#L124?message=Update%20Docs)]
</div>
**LLM Docstring**

Render the 2D spectrum as a filled and/or line contour plot, optionally removing the median baseline from the intensities first and using a symmetric (about zero) color/level range by default.
  - `plot_filled`: `bool`
    > whether to draw a filled contour plot (via `ContourPlot`)
  - `contour_line_style`: `dict | bool | None`
    > style options for an overlaid contour-line plot; `True`/`None` for default line styling, `False` to omit the line plot entirely (only meaningful together with `plot_filled`)
  - `figure`: `object | None`
    > an existing figure to draw into
  - `symmetric_range`: `bool`
    > whether to force the color/level range to be symmetric about zero
  - `remove_baseline`: `bool`
    > whether to subtract the median intensity before plotting
  - `vmin`: `float | None`
    > the minimum value for the color scale; derived from the data if not given
  - `vmax`: `float | None`
    > the maximum value for the color scale; derived from the data if not given
  - `levels`: `int | Iterable[float] | None`
    > the number of contour levels, or explicit level values
  - `opts`: `dict`
    > extra plotting options forwarded to `ContourPlot`/`ContourLinePlot`
  - `:returns`: `object`
    > the resulting figure
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