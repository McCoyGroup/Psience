## <a id="Psience.DVR.Wavefunctions.DVRWavefunctions">DVRWavefunctions</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Wavefunctions.py#L86)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Wavefunctions.py#L86?message=Update%20Docs)]
</div>



<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, energies=None, wavefunctions=None, wavefunction_class=<class 'Psience.DVR.Wavefunctions.DVRWavefunction'>, results: Psience.DVR.BaseDVR.DVRResults = None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Wavefunctions.py#L88)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Wavefunctions.py#L88?message=Update%20Docs)]
</div>

<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Wavefunctions.py#L91)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Wavefunctions.py#L91?message=Update%20Docs)]
</div>

<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.__len__" class="docs-object-method">&nbsp;</a> 
```python
__len__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Wavefunctions.py#L98)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Wavefunctions.py#L98?message=Update%20Docs)]
</div>

<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.__iter__" class="docs-object-method">&nbsp;</a> 
```python
__iter__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Wavefunctions.py#L100)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Wavefunctions.py#L100?message=Update%20Docs)]
</div>

<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Wavefunctions.py#L103)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Wavefunctions.py#L103?message=Update%20Docs)]
</div>

Provides a single `DVRWavefunction` or slice of `DVRWavefunctions`
- `item`: `Any`
    >No description...
- `:returns`: `DVRWavefunction | DVRWavefunctions`
    >No description...

<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.plot" class="docs-object-method">&nbsp;</a> 
```python
plot(self, figure=None, graphics_class=None, plot_style=None, scaling=1, shift=0, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Wavefunctions.py#L127)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Wavefunctions.py#L127?message=Update%20Docs)]
</div>

Plots the held wavefunctions
- `figure`: `Any`
    >No description...
- `graphics_class`: `Any`
    >No description...
- `plot_style`: `Any`
    >No description...
- `scaling`: `Any`
    >No description...
- `shift`: `Any`
    >No description...
- `opts`: `Any`
    >No description...
- `:returns`: `Graphics`
    >No description...

<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.expectation" class="docs-object-method">&nbsp;</a> 
```python
expectation(self, op, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Wavefunctions.py#L163)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Wavefunctions.py#L163?message=Update%20Docs)]
</div>

Computes the expectation value of operator op over the wavefunction other and self
- `other`: `DVRWavefunctions | np.ndarray`
    >No description...
- `op`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.DVR.Wavefunctions.DVRWavefunctions.probability_density" class="docs-object-method">&nbsp;</a> 
```python
probability_density(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Wavefunctions.py#L179)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Wavefunctions.py#L179?message=Update%20Docs)]
</div>

Computes the probability density of the set of wavefunctions
- `:returns`: `_`
    >No description...

 </div>
</div>




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/DVR/Wavefunctions/DVRWavefunctions.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/DVR/Wavefunctions/DVRWavefunctions.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/DVR/Wavefunctions/DVRWavefunctions.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/DVR/Wavefunctions/DVRWavefunctions.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Wavefunctions.py#L86?message=Update%20Docs)