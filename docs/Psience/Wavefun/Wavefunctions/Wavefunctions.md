## <a id="Psience.Wavefun.Wavefunctions.Wavefunctions">Wavefunctions</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Wavefun/Wavefunctions.py#L58)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Wavefun/Wavefunctions.py#L58?message=Update%20Docs)]
</div>

An object representing a set of wavefunctions.
Provides concrete, but potentially inefficient methods for doing all the wavefunction ops.

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

```python
wavefunction_class: type
```
<a id="Psience.Wavefun.Wavefunctions.Wavefunctions.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, energies=None, wavefunctions=None, indices=None, wavefunction_class=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Wavefun/Wavefunctions.py#L64)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Wavefun/Wavefunctions.py#L64?message=Update%20Docs)]
</div>

<a id="Psience.Wavefun.Wavefunctions.Wavefunctions.get_wavefunctions" class="docs-object-method">&nbsp;</a> 
```python
get_wavefunctions(self, which): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Wavefun/Wavefunctions.py#L73)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Wavefun/Wavefunctions.py#L73?message=Update%20Docs)]
</div>

<a id="Psience.Wavefun.Wavefunctions.Wavefunctions.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Wavefun/Wavefunctions.py#L87)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Wavefun/Wavefunctions.py#L87?message=Update%20Docs)]
</div>

Returns a single Wavefunction object

<a id="Psience.Wavefun.Wavefunctions.Wavefunctions.__iter__" class="docs-object-method">&nbsp;</a> 
```python
__iter__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Wavefun/Wavefunctions.py#L91)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Wavefun/Wavefunctions.py#L91?message=Update%20Docs)]
</div>

<a id="Psience.Wavefun.Wavefunctions.Wavefunctions.frequencies" class="docs-object-method">&nbsp;</a> 
```python
frequencies(self, start_at=0): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Wavefun/Wavefunctions.py#L97)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Wavefun/Wavefunctions.py#L97?message=Update%20Docs)]
</div>

<a id="Psience.Wavefun.Wavefunctions.Wavefunctions.plot" class="docs-object-method">&nbsp;</a> 
```python
plot(self, figure=None, graphics_class=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Wavefun/Wavefunctions.py#L100)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Wavefun/Wavefunctions.py#L100?message=Update%20Docs)]
</div>

Plots all of the wavefunctions on one set of axes
- `:returns`: `_`
    >
- `opts`: `Any`
    >

<a id="Psience.Wavefun.Wavefunctions.Wavefunctions.expectation" class="docs-object-method">&nbsp;</a> 
```python
expectation(self, op, other=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Wavefun/Wavefunctions.py#L136)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Wavefun/Wavefunctions.py#L136?message=Update%20Docs)]
</div>

Computes the expectation value of operator op over the wavefunction other and self
- `:returns`: `_`
    >
- `op`: `Any`
    >
- `other`: `Wavefunctions`
    >

<a id="Psience.Wavefun.Wavefunctions.Wavefunctions.overlap" class="docs-object-method">&nbsp;</a> 
```python
overlap(self, other): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Wavefun/Wavefunctions.py#L157)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Wavefun/Wavefunctions.py#L157?message=Update%20Docs)]
</div>

<a id="Psience.Wavefun.Wavefunctions.Wavefunctions.coordinate" class="docs-object-method">&nbsp;</a> 
```python
coordinate(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Wavefun/Wavefunctions.py#L160)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Wavefun/Wavefunctions.py#L160?message=Update%20Docs)]
</div>

Provides the coordinate operator in the wavefunction basis
- `:returns`: `_`
    >

<a id="Psience.Wavefun.Wavefunctions.Wavefunctions.momentum" class="docs-object-method">&nbsp;</a> 
```python
momentum(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Wavefun/Wavefunctions.py#L168)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Wavefun/Wavefunctions.py#L168?message=Update%20Docs)]
</div>

Provides the real part of the representation of the momentum operator in the wavefunction basis
- `:returns`: `_`
    >

<a id="Psience.Wavefun.Wavefunctions.Wavefunctions.laplacian" class="docs-object-method">&nbsp;</a> 
```python
laplacian(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Wavefun/Wavefunctions.py#L176)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Wavefun/Wavefunctions.py#L176?message=Update%20Docs)]
</div>

Provides the representation of the laplacian in the wavefunction basis
- `:returns`: `_`
    >

<a id="Psience.Wavefun.Wavefunctions.Wavefunctions.kinetic_energy" class="docs-object-method">&nbsp;</a> 
```python
kinetic_energy(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Wavefun/Wavefunctions.py#L184)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Wavefun/Wavefunctions.py#L184?message=Update%20Docs)]
</div>

Provides the representation of the KE in the wavefunction basis
- `:returns`: `_`
    >

 </div>
</div>






___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Wavefun/Wavefunctions/Wavefunctions.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Wavefun/Wavefunctions/Wavefunctions.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Wavefun/Wavefunctions/Wavefunctions.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Wavefun/Wavefunctions/Wavefunctions.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/Wavefun/Wavefunctions.py#L58?message=Update%20Docs)