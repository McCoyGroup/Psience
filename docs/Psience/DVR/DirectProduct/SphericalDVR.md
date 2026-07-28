## <a id="Psience.DVR.DirectProduct.SphericalDVR">SphericalDVR</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/DirectProduct.py#L384)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/DirectProduct.py#L384?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.DVR.DirectProduct.SphericalDVR.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, r_max, divs, **base_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/DirectProduct.py#L385)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/DirectProduct.py#L385?message=Update%20Docs)]
</div>
**LLM Docstring**

Set up a 3D spherical-coordinate DVR by combining a `RadialDVR` (over `(0, r_max)`), a `RingDVR`, and a `PolarDVR` via `DirectProductDVR`. Note the `PolarDVR` component is constructed with `domain=(0, 2*np.pi)`, matching the `RingDVR` above it rather than the `(0, pi)` domain `PolarDVR` normally defaults to and is designed for -- this looks like a copy-paste bug.
  - `r_max`: `float`
    > the maximum radial extent
  - `divs`: `tuple[int, int]`
    > the `(radial_divs, angular_divs)` number of grid points for the radial and angular dimensions (the same `angular_divs` is used for both the ring and polar dimensions)
  - `base_opts`: `dict`
    > extra options forwarded to the base `DirectProductDVR.__init__`
  - `:returns`: `None`
    > None
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/DVR/DirectProduct/SphericalDVR.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/DVR/DirectProduct/SphericalDVR.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/DVR/DirectProduct/SphericalDVR.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/DVR/DirectProduct/SphericalDVR.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/DirectProduct.py#L384?message=Update%20Docs)   
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