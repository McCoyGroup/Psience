## <a id="Psience.Spectra.BaseSpectrum.MultiSpectrum">MultiSpectrum</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/BaseSpectrum.py#L411)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum.py#L411?message=Update%20Docs)]
</div>

A wrapper for multiple spectra, really just for the plotting convenience







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.Spectra.BaseSpectrum.MultiSpectrum.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, spectra: 'Iterable[BaseSpectrum]', **meta): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/BaseSpectrum.py#L419)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum.py#L419?message=Update%20Docs)]
</div>

  - `frequencies`: `np.ndarray`
    > frequency list
  - `intensities`: `np.ndarray`
    > intensity list
  - `meta`: `Any`
    > metadata


<a id="Psience.Spectra.BaseSpectrum.MultiSpectrum.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/BaseSpectrum/MultiSpectrum.py#L431)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum/MultiSpectrum.py#L431?message=Update%20Docs)]
</div>


<a id="Psience.Spectra.BaseSpectrum.MultiSpectrum.frequency_filter" class="docs-object-method">&nbsp;</a> 
```python
frequency_filter(self, freq_min, freq_max): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/BaseSpectrum/MultiSpectrum.py#L439)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum/MultiSpectrum.py#L439?message=Update%20Docs)]
</div>
Filters by frequencies >= `freq_min` and <= `freq_max`
  - `freq_min`: `float`
    > min frequency
  - `freq_max`: `float`
    > max frequency
  - `:returns`: `MultiSpectrum`
    > s
u
b
s
p
e
c
t
r
u
m


<a id="Psience.Spectra.BaseSpectrum.MultiSpectrum.intensity_filter" class="docs-object-method">&nbsp;</a> 
```python
intensity_filter(self, int_min, int_max): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/BaseSpectrum/MultiSpectrum.py#L456)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum/MultiSpectrum.py#L456?message=Update%20Docs)]
</div>
Filters by intensities >= `int_min` and <= `int_max`
  - `int_min`: `float`
    > min intensity
  - `int_max`: `float`
    > max intensity
  - `:returns`: `BaseSpectrum`
    > s
u
b
s
p
e
c
t
r
u
m


<a id="Psience.Spectra.BaseSpectrum.MultiSpectrum.plot" class="docs-object-method">&nbsp;</a> 
```python
plot(self, figure=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/BaseSpectrum/MultiSpectrum.py#L473)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum/MultiSpectrum.py#L473?message=Update%20Docs)]
</div>
A just plots all the spectra on the same figure
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Spectra/BaseSpectrum/MultiSpectrum.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Spectra/BaseSpectrum/MultiSpectrum.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Spectra/BaseSpectrum/MultiSpectrum.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Spectra/BaseSpectrum/MultiSpectrum.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum.py#L411?message=Update%20Docs)   
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