## <a id="Psience.Spectra.BaseSpectrum.DiscreteSpectrum">DiscreteSpectrum</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Spectra/BaseSpectrum.py#L140)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Spectra/BaseSpectrum.py#L140?message=Update%20Docs)]
</div>

Concrete implementation of `BaseSpectrum` that exists
solely to allow for plotting and broadening.

<a id="Psience.Spectra.BaseSpectrum.DiscreteSpectrum.plot" class="docs-object-method">&nbsp;</a> 
```python
plot(self, figure=None, plot_style=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Spectra/BaseSpectrum.py#L146)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Spectra/BaseSpectrum.py#L146?message=Update%20Docs)]
</div>

Plots a spectrum using `McUtils.Plots.StickSpectrum`
- `figure`: `None | McUtils.Plots.Graphics`
    >figure to plot the spectrum on
- `opts`: `Any`
    >any of the many, many options supported by `McUtils.Plots.Graphics`
- `:returns`: `_`
    >No description...

<a id="Psience.Spectra.BaseSpectrum.DiscreteSpectrum.broaden" class="docs-object-method">&nbsp;</a> 
```python
broaden(self, broadening_type='gaussian', breadth=10): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Spectra/BaseSpectrum.py#L169)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Spectra/BaseSpectrum.py#L169?message=Update%20Docs)]
</div>

Applies a broadening to the spectrum
- `broadening_type`: `Any`
    >No description...
- `breadth`: `Any`
    >No description...
- `:returns`: `_`
    >No description...



___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/ci/docs/Psience/Spectra/BaseSpectrum/DiscreteSpectrum.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/ci/docs/Psience/Spectra/BaseSpectrum/DiscreteSpectrum.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/ci/docs/Psience/Spectra/BaseSpectrum/DiscreteSpectrum.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/ci/docs/Psience/Spectra/BaseSpectrum/DiscreteSpectrum.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Spectra/BaseSpectrum.py#L140?message=Update%20Docs)