## <a id="Psience.Spectra.BaseSpectrum.BaseSpectrum">BaseSpectrum</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/BaseSpectrum.py#L17)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum.py#L17?message=Update%20Docs)]
</div>

Base class to support spectral operation



<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

<a id="Psience.Spectra.BaseSpectrum.BaseSpectrum.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, frequencies, intensities, **meta): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/BaseSpectrum.py#L21)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum.py#L21?message=Update%20Docs)]
</div>


- `meta`: `Any`
    >metadata
- `intensities`: `np.ndarray`
    >intensity list
- `frequencies`: `np.ndarray`
    >frequency list

<a id="Psience.Spectra.BaseSpectrum.BaseSpectrum.take_subspectrum" class="docs-object-method">&nbsp;</a> 
```python
take_subspectrum(self, pos): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/BaseSpectrum.py#L34)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum.py#L34?message=Update%20Docs)]
</div>

Takes a subset of frequencies/intensities specified by `pos`
- `:returns`: `_`
    >
- `pos`: `Any`
    >

<a id="Psience.Spectra.BaseSpectrum.BaseSpectrum.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/BaseSpectrum.py#L46)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum.py#L46?message=Update%20Docs)]
</div>

<a id="Psience.Spectra.BaseSpectrum.BaseSpectrum.frequency_filter" class="docs-object-method">&nbsp;</a> 
```python
frequency_filter(self, freq_min, freq_max): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/BaseSpectrum.py#L54)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum.py#L54?message=Update%20Docs)]
</div>

Filters by frequencies >= `freq_min` and <= `freq_max`
- `:returns`: `BaseSpectrum`
    >s
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
- `freq_max`: `float`
    >max frequency
- `freq_min`: `float`
    >min frequency

<a id="Psience.Spectra.BaseSpectrum.BaseSpectrum.intensity_filter" class="docs-object-method">&nbsp;</a> 
```python
intensity_filter(self, int_min, int_max): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/BaseSpectrum.py#L69)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum.py#L69?message=Update%20Docs)]
</div>

Filters by intensities >= `int_min` and <= `int_max`
- `:returns`: `BaseSpectrum`
    >s
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
- `int_max`: `float`
    >max intensity
- `int_min`: `float`
    >min intensity

<a id="Psience.Spectra.BaseSpectrum.BaseSpectrum.save" class="docs-object-method">&nbsp;</a> 
```python
save(self, file): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/BaseSpectrum.py#L85)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum.py#L85?message=Update%20Docs)]
</div>

Saves the spectrum in JSON format
- `:returns`: `_`
    >
- `file`: `Any`
    >str | file-like object

<a id="Psience.Spectra.BaseSpectrum.BaseSpectrum.load" class="docs-object-method">&nbsp;</a> 
```python
load(file): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/BaseSpectrum.py#L106)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum.py#L106?message=Update%20Docs)]
</div>

Saves a spectrum from a JSON file
- `:returns`: `_`
    >
- `file`: `Any`
    >str | file-like object

<a id="Psience.Spectra.BaseSpectrum.BaseSpectrum.plot" class="docs-object-method">&nbsp;</a> 
```python
plot(self, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/BaseSpectrum.py#L125)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum.py#L125?message=Update%20Docs)]
</div>

A stub so that subclasses can implement their own `plot` methods
- `:returns`: `_`
    >
- `opts`: `Any`
    >plotting options to be fed through to whatever the plotting function uses

 </div>
</div>




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Spectra/BaseSpectrum/BaseSpectrum.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Spectra/BaseSpectrum/BaseSpectrum.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Spectra/BaseSpectrum/BaseSpectrum.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Spectra/BaseSpectrum/BaseSpectrum.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum.py#L17?message=Update%20Docs)