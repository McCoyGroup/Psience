## <a id="Psience.DVR.DirectProduct.RingNDDVR">RingNDDVR</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/DirectProduct.py#L357)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/DirectProduct.py#L357?message=Update%20Docs)]
</div>

Provides an ND-DVR for products of periodic (0, 2Pi) ranges







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.DVR.DirectProduct.RingNDDVR.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, divs, **base_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/DirectProduct.py#L362)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/DirectProduct.py#L362?message=Update%20Docs)]
</div>
**LLM Docstring**

Set up an N-dimensional DVR over a product of periodic `(0, 2*pi)` ring domains by building an independent `RingDVR` for each requested division count and combining them via `DirectProductDVR`.
  - `divs`: `Iterable[int]`
    > the per-dimension number of grid points
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/DVR/DirectProduct/RingNDDVR.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/DVR/DirectProduct/RingNDDVR.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/DVR/DirectProduct/RingNDDVR.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/DVR/DirectProduct/RingNDDVR.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/DirectProduct.py#L357?message=Update%20Docs)   
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