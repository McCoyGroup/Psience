## <a id="Psience.DVR.DirectProduct.DirectProductDVR">DirectProductDVR</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/DirectProduct.py#L17)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/DirectProduct.py#L17?message=Update%20Docs)]
</div>





<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

<a id="Psience.DVR.DirectProduct.DirectProductDVR.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, dvrs_1D, zero_threshold=1e-14, **base_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/DirectProduct.py#L18)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/DirectProduct.py#L18?message=Update%20Docs)]
</div>


- `base_opts`: `Any`
    >
- `dvrs_1D`: `Iterable[AbstractDVR]`
    >a series of 1D DVRs that can provide the inputs we'll product together

<a id="Psience.DVR.DirectProduct.DirectProductDVR.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/DirectProduct.py#L32)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/DirectProduct.py#L32?message=Update%20Docs)]
</div>

<a id="Psience.DVR.DirectProduct.DirectProductDVR.get_grid" class="docs-object-method">&nbsp;</a> 
```python
get_grid(self, domain=None, divs=None, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/DirectProduct.py#L39)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/DirectProduct.py#L39?message=Update%20Docs)]
</div>

<a id="Psience.DVR.DirectProduct.DirectProductDVR.grid" class="docs-object-method">&nbsp;</a> 
```python
grid(self, domain=None, divs=None, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/DirectProduct.py#L55)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/DirectProduct.py#L55?message=Update%20Docs)]
</div>

<a id="Psience.DVR.DirectProduct.DirectProductDVR.get_kinetic_energy" class="docs-object-method">&nbsp;</a> 
```python
get_kinetic_energy(self, grid=None, mass=None, hb=1, g=None, g_deriv=None, logger=None, include_kinetic_coupling=True, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/DirectProduct.py#L68)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/DirectProduct.py#L68?message=Update%20Docs)]
</div>

<a id="Psience.DVR.DirectProduct.DirectProductDVR.kinetic_energy" class="docs-object-method">&nbsp;</a> 
```python
kinetic_energy(self, grid=None, mass=None, hb=1, g=None, g_deriv=None, logger=None, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/DirectProduct.py#L246)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/DirectProduct.py#L246?message=Update%20Docs)]
</div>

Computes the N-dimensional kinetic energy
- `:returns`: `_`
    >
- `g_deriv`: `Any`
    >
- `g`: `Any`
    >
- `hb`: `Any`
    >
- `mass`: `Any`
    >
- `grid`: `Any`
    >

 </div>
</div>




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/DVR/DirectProduct/DirectProductDVR.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/DVR/DirectProduct/DirectProductDVR.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/DVR/DirectProduct/DirectProductDVR.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/DVR/DirectProduct/DirectProductDVR.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/DirectProduct.py#L17?message=Update%20Docs)