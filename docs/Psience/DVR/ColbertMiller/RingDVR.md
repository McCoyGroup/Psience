## <a id="Psience.DVR.ColbertMiller.RingDVR">RingDVR</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/ColbertMiller.py#L107)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/ColbertMiller.py#L107?message=Update%20Docs)]
</div>

Provides a DVR for working on the (0, 2Pi) range with periodicity from Colbert and Miller







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.DVR.ColbertMiller.RingDVR.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, domain=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/ColbertMiller.py#L112)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/ColbertMiller.py#L112?message=Update%20Docs)]
</div>
**LLM Docstring**

Set up a Colbert-Miller DVR on the periodic `(0, 2*pi)` ring domain, defaulting the domain to `(0, 2*pi)` if not given.
  - `domain`: `tuple | None`
    > the coordinate domain; defaults to `(0, 2*pi)`
  - `opts`: `dict`
    > extra options forwarded to the base `BaseDVR.__init__`
  - `:returns`: `None`
    > None


<a id="Psience.DVR.ColbertMiller.RingDVR.get_grid" class="docs-object-method">&nbsp;</a> 
```python
get_grid(self, domain=None, divs=None, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/ColbertMiller/RingDVR.py#L129)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/ColbertMiller/RingDVR.py#L129?message=Update%20Docs)]
</div>
Provides the Colbert-Miller 1D grid for the [0, 2Pi] range
  - `domain`: `Any`
    > 
  - `divs`: `Any`
    > 
  - `kw`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.DVR.ColbertMiller.RingDVR.get_kinetic_energy" class="docs-object-method">&nbsp;</a> 
```python
get_kinetic_energy(self, grid=None, mass=1, hb=1, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/ColbertMiller/RingDVR.py#L147)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/ColbertMiller/RingDVR.py#L147?message=Update%20Docs)]
</div>
Colbert-Miller kinetic energy for the [0, 2pi] range
  - `grid`: `Any`
    > 
  - `mass`: `Any`
    > 
  - `hb`: `Any`
    > 
  - `kw`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.DVR.ColbertMiller.RingDVR.real_momentum" class="docs-object-method">&nbsp;</a> 
```python
real_momentum(self, grid=None, hb=1, **kw): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/ColbertMiller/RingDVR.py#L181)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/ColbertMiller/RingDVR.py#L181?message=Update%20Docs)]
</div>
Provides the real part of the momentum for the [0, 2pi] range
  - `grid`: `Any`
    > 
  - `hb`: `Any`
    > 
  - `kw`: `Any`
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/DVR/ColbertMiller/RingDVR.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/DVR/ColbertMiller/RingDVR.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/DVR/ColbertMiller/RingDVR.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/DVR/ColbertMiller/RingDVR.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/ColbertMiller.py#L107?message=Update%20Docs)   
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