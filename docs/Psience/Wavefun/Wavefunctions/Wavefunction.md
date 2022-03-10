## <a id="Psience.Wavefun.Wavefunctions.Wavefunction">Wavefunction</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Wavefun/Wavefunctions.py#L17)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Wavefun/Wavefunctions.py#L17?message=Update%20Docs)]
</div>

Represents a single wavefunction object

<a id="Psience.Wavefun.Wavefunctions.Wavefunction.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, energy, data, parent=None, index=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Wavefun/Wavefunctions.py#L19)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Wavefun/Wavefunctions.py#L19?message=Update%20Docs)]
</div>

<a id="Psience.Wavefun.Wavefunctions.Wavefunction.plot" class="docs-object-method">&nbsp;</a> 
```python
plot(self, figure=None, index=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Wavefun/Wavefunctions.py#L25)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Wavefun/Wavefunctions.py#L25?message=Update%20Docs)]
</div>

Uses McUtils to plot the wavefunction on the passed figure (makes a new one if none)
- `figure`: `Graphics | Graphics3D`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.Wavefun.Wavefunctions.Wavefunction.expectation" class="docs-object-method">&nbsp;</a> 
```python
expectation(self, op, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Wavefun/Wavefunctions.py#L35)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Wavefun/Wavefunctions.py#L35?message=Update%20Docs)]
</div>

Computes the expectation value of operator op over the wavefunction other and self
- `other`: `Wavefunction`
    >No description...
- `op`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.Wavefun.Wavefunctions.Wavefunction.probability_density" class="docs-object-method">&nbsp;</a> 
```python
probability_density(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Wavefun/Wavefunctions.py#L47)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Wavefun/Wavefunctions.py#L47?message=Update%20Docs)]
</div>

Computes the probability density of the current wavefunction
- `:returns`: `_`
    >No description...



___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/ci/docs/Psience/Wavefun/Wavefunctions/Wavefunction.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/ci/docs/Psience/Wavefun/Wavefunctions/Wavefunction.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/ci/docs/Psience/Wavefun/Wavefunctions/Wavefunction.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/ci/docs/Psience/Wavefun/Wavefunctions/Wavefunction.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Wavefun/Wavefunctions.py#L17?message=Update%20Docs)