## <a id="Psience.Spectra.BaseSpectrum.BroadenedSpectrum">BroadenedSpectrum</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/BaseSpectrum.py#L214)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum.py#L214?message=Update%20Docs)]
</div>

A stick spectrum with associated broadening function

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

<a id="Psience.Spectra.BaseSpectrum.BroadenedSpectrum.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, frequencies, intensities, broadening_type='gaussian', breadth=10, **meta): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/BaseSpectrum.py#L219)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum.py#L219?message=Update%20Docs)]
</div>


- `frequencies`: `Any`
    >No description...
- `intensities`: `Any`
    >No description...
- `broadening_type`: `"gaussian" | "lorentzian" | function`
    >the type of broadening to apply (can be any function)
- `breadth`: `Any`
    >the breadth or list of breads for the peaks in the spectrum
- `meta`: `Any`
    >No description...

<a id="Psience.Spectra.BaseSpectrum.BroadenedSpectrum.plot" class="docs-object-method">&nbsp;</a> 
```python
plot(self, step_size=0.5, freq_min=None, freq_max=None, figure=None, plot_style=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Spectra/BaseSpectrum.py#L324)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum.py#L324?message=Update%20Docs)]
</div>

Applies the broadening then plots it using `McUtils.Plots.Plot`
- `step_size`: `Any`
    >step size to use when getting evaluation points for evaluating the broadening
- `freq_min`: `Any`
    >min freq for evaluation
- `freq_max`: `Any`
    >max freq for evaluation
- `figure`: `Any`
    >No description...
- `plot_style`: `Any`
    >No description...
- `opts`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

 </div>
</div>




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Spectra/BaseSpectrum/BroadenedSpectrum.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Spectra/BaseSpectrum/BroadenedSpectrum.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Spectra/BaseSpectrum/BroadenedSpectrum.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Spectra/BaseSpectrum/BroadenedSpectrum.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/Spectra/BaseSpectrum.py#L214?message=Update%20Docs)