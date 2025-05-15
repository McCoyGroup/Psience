## <a id="Psience.DVR.DirectProduct.DirectProductDVR">DirectProductDVR</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/DirectProduct.py#L17)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/DirectProduct.py#L17?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.DVR.DirectProduct.DirectProductDVR.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, dvrs_1D, zero_threshold=1e-14, **base_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/DirectProduct/DirectProductDVR.py#L18)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/DirectProduct/DirectProductDVR.py#L18?message=Update%20Docs)]
</div>

  - `dvrs_1D`: `Iterable[AbstractDVR]`
    > a series of 1D DVRs that can provide the inputs we'll product together
  - `base_opts`: `Any`
    >


<a id="Psience.DVR.DirectProduct.DirectProductDVR.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/DirectProduct/DirectProductDVR.py#L32)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/DirectProduct/DirectProductDVR.py#L32?message=Update%20Docs)]
</div>


<a id="Psience.DVR.DirectProduct.DirectProductDVR.get_grid" class="docs-object-method">&nbsp;</a> 
```python
get_grid(self, domain=None, divs=None, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/DirectProduct/DirectProductDVR.py#L39)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/DirectProduct/DirectProductDVR.py#L39?message=Update%20Docs)]
</div>


<a id="Psience.DVR.DirectProduct.DirectProductDVR.grid" class="docs-object-method">&nbsp;</a> 
```python
grid(self, domain=None, divs=None, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/DirectProduct/DirectProductDVR.py#L55)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/DirectProduct/DirectProductDVR.py#L55?message=Update%20Docs)]
</div>


<a id="Psience.DVR.DirectProduct.DirectProductDVR.get_kinetic_energy" class="docs-object-method">&nbsp;</a> 
```python
get_kinetic_energy(self, grid=None, mass=None, hb=1, g=None, g_deriv=None, logger=None, include_kinetic_coupling=True, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/DirectProduct/DirectProductDVR.py#L68)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/DirectProduct/DirectProductDVR.py#L68?message=Update%20Docs)]
</div>


<a id="Psience.DVR.DirectProduct.DirectProductDVR.kinetic_energy" class="docs-object-method">&nbsp;</a> 
```python
kinetic_energy(self, grid=None, mass=None, hb=1, g=None, g_deriv=None, logger=None, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DVR/DirectProduct/DirectProductDVR.py#L246)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/DirectProduct/DirectProductDVR.py#L246?message=Update%20Docs)]
</div>
Computes the N-dimensional kinetic energy
  - `grid`: `Any`
    > 
  - `mass`: `Any`
    > 
  - `hb`: `Any`
    > 
  - `g`: `Any`
    > 
  - `g_deriv`: `Any`
    > 
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/DVR/DirectProduct/DirectProductDVR.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/DVR/DirectProduct/DirectProductDVR.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/DVR/DirectProduct/DirectProductDVR.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/DVR/DirectProduct/DirectProductDVR.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/DVR/DirectProduct.py#L17?message=Update%20Docs)   
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