## <a id="Psience.Spectra.BaseSpectrum.BaseSpectrum">BaseSpectrum</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/BaseSpectrum.py#L19)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum.py#L19?message=Update%20Docs)]
</div>

Base class to support spectral operation







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.Spectra.BaseSpectrum.BaseSpectrum.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, frequencies, intensities, **meta): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/BaseSpectrum.py#L23)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum.py#L23?message=Update%20Docs)]
</div>

  - `frequencies`: `np.ndarray`
    > frequency list
  - `intensities`: `np.ndarray`
    > intensity list
  - `meta`: `Any`
    > metadata


<a id="Psience.Spectra.BaseSpectrum.BaseSpectrum.take_subspectrum" class="docs-object-method">&nbsp;</a> 
```python
take_subspectrum(self, pos): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/BaseSpectrum/BaseSpectrum.py#L36)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum/BaseSpectrum.py#L36?message=Update%20Docs)]
</div>
Takes a subset of frequencies/intensities specified by `pos`
  - `pos`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Spectra.BaseSpectrum.BaseSpectrum.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/BaseSpectrum/BaseSpectrum.py#L48)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum/BaseSpectrum.py#L48?message=Update%20Docs)]
</div>
**LLM Docstring**

Index into the spectrum: a scalar index returns the corresponding `(frequency, intensity)` pair, while any other index (slice, array, etc.) returns a subspectrum via `take_subspectrum`.
  - `item`: `int | slice | np.ndarray`
    > the index or slice to apply
  - `:returns`: `tuple[float, float] | BaseSpectrum`
    > a `(frequency, intensity)` pair for a scalar index, or a subspectrum otherwise


<a id="Psience.Spectra.BaseSpectrum.BaseSpectrum.frequency_filter" class="docs-object-method">&nbsp;</a> 
```python
frequency_filter(self, freq_min, freq_max): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/BaseSpectrum/BaseSpectrum.py#L66)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum/BaseSpectrum.py#L66?message=Update%20Docs)]
</div>
Filters by frequencies >= `freq_min` and <= `freq_max`
  - `freq_min`: `float`
    > min frequency
  - `freq_max`: `float`
    > max frequency
  - `:returns`: `BaseSpectrum`
    > subspectrum


<a id="Psience.Spectra.BaseSpectrum.BaseSpectrum.intensity_filter" class="docs-object-method">&nbsp;</a> 
```python
intensity_filter(self, int_min, int_max): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/BaseSpectrum/BaseSpectrum.py#L81)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum/BaseSpectrum.py#L81?message=Update%20Docs)]
</div>
Filters by intensities >= `int_min` and <= `int_max`
  - `int_min`: `float`
    > min intensity
  - `int_max`: `float`
    > max intensity
  - `:returns`: `BaseSpectrum`
    > subspectrum


<a id="Psience.Spectra.BaseSpectrum.BaseSpectrum.save" class="docs-object-method">&nbsp;</a> 
```python
save(self, file): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/BaseSpectrum/BaseSpectrum.py#L97)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum/BaseSpectrum.py#L97?message=Update%20Docs)]
</div>
Saves the spectrum in JSON format
  - `file`: `Any`
    > str | file-like object
  - `:returns`: `_`
    >


<a id="Psience.Spectra.BaseSpectrum.BaseSpectrum.load" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
load(cls, file): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L118)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L118?message=Update%20Docs)]
</div>
Saves a spectrum from a JSON file
  - `file`: `Any`
    > str | file-like object
  - `:returns`: `_`
    >


<a id="Psience.Spectra.BaseSpectrum.BaseSpectrum.plot" class="docs-object-method">&nbsp;</a> 
```python
plot(self, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/BaseSpectrum/BaseSpectrum.py#L137)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum/BaseSpectrum.py#L137?message=Update%20Docs)]
</div>
A stub so that subclasses can implement their own `plot` methods
  - `opts`: `Any`
    > plotting options to be fed through to whatever the plotting function uses
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Spectra/BaseSpectrum/BaseSpectrum.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Spectra/BaseSpectrum/BaseSpectrum.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Spectra/BaseSpectrum/BaseSpectrum.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Spectra/BaseSpectrum/BaseSpectrum.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum.py#L19?message=Update%20Docs)   
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