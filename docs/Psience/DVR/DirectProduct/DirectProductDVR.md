## <a id="Psience.DVR.DirectProduct.DirectProductDVR">DirectProductDVR</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/DVR/DirectProduct.py#L17)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/DVR/DirectProduct.py#L17?message=Update%20Docs)]
</div>



<a id="Psience.DVR.DirectProduct.DirectProductDVR.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, dvrs_1D, **base_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/DVR/DirectProduct.py#L19)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/DVR/DirectProduct.py#L19?message=Update%20Docs)]
</div>


- `dvrs_1D`: `Iterable[AbstractDVR]`
    >a series of 1D DVRs that can provide the inputs we'll product together
- `base_opts`: `Any`
    >No description...

<a id="Psience.DVR.DirectProduct.DirectProductDVR.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/DVR/DirectProduct.py#L31)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/DVR/DirectProduct.py#L31?message=Update%20Docs)]
</div>

<a id="Psience.DVR.DirectProduct.DirectProductDVR.get_grid" class="docs-object-method">&nbsp;</a> 
```python
get_grid(self, domain=None, divs=None, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/DVR/DirectProduct.py#L38)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/DVR/DirectProduct.py#L38?message=Update%20Docs)]
</div>

<a id="Psience.DVR.DirectProduct.DirectProductDVR.grid" class="docs-object-method">&nbsp;</a> 
```python
grid(self, domain=None, divs=None, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/DVR/DirectProduct.py#L49)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/DVR/DirectProduct.py#L49?message=Update%20Docs)]
</div>

<a id="Psience.DVR.DirectProduct.DirectProductDVR.get_kinetic_energy" class="docs-object-method">&nbsp;</a> 
```python
get_kinetic_energy(self, grid=None, mass=None, hb=1, g=None, g_deriv=None, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/DVR/DirectProduct.py#L61)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/DVR/DirectProduct.py#L61?message=Update%20Docs)]
</div>

<a id="Psience.DVR.DirectProduct.DirectProductDVR.kinetic_energy" class="docs-object-method">&nbsp;</a> 
```python
kinetic_energy(self, grid=None, mass=None, hb=1, g=None, g_deriv=None, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/DVR/DirectProduct.py#L214)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/DVR/DirectProduct.py#L214?message=Update%20Docs)]
</div>

Computes the N-dimensional kinetic energy
- `grid`: `Any`
    >No description...
- `mass`: `Any`
    >No description...
- `hb`: `Any`
    >No description...
- `g`: `Any`
    >No description...
- `g_deriv`: `Any`
    >No description...
- `:returns`: `_`
    >No description...



___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/ci/docs/Psience/DVR/DirectProduct/DirectProductDVR.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/ci/docs/Psience/DVR/DirectProduct/DirectProductDVR.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/ci/docs/Psience/DVR/DirectProduct/DirectProductDVR.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/ci/docs/Psience/DVR/DirectProduct/DirectProductDVR.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/DVR/DirectProduct.py#L17?message=Update%20Docs)