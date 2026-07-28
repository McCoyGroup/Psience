## <a id="Psience.DVR.Extensions.PotentialOptimizedDVR">PotentialOptimizedDVR</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Extensions.py#L199)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Extensions.py#L199?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.DVR.Extensions.PotentialOptimizedDVR.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, wfns_1D: 'Iterable[DVRWavefunctions]', **base_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/DVR/Extensions.py#L200)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Extensions.py#L200?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a multi-dimensional DVR whose per-dimension bases are potential-optimized (built from a given set of 1D wavefunctions rather than a fixed grid), via `WavefunctionBasisDVR` and `DirectProductDVR`.
  - `wfns_1D`: `Iterable[DVRWavefunctions]`
    > the per-dimension 1D wavefunctions to use as the optimized basis
  - `base_opts`: `dict`
    > extra options forwarded to the base `DirectProductDVR.__init__`
  - `:returns`: `None`
    > None


<a id="Psience.DVR.Extensions.PotentialOptimizedDVR.from_minimum" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_minimum(cls, base_dvr: 'DirectProductDVR|SelfConsistentDVR', **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L222)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L222?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a `PotentialOptimizedDVR` using the SCF wavefunctions computed from an initial (unconverged, single-iteration) guess at the potential minimum, wrapping `base_dvr` in a `SelfConsistentDVR` first if it isn't already one.
  - `base_dvr`: `DirectProductDVR | SelfConsistentDVR`
    > the base multi-dimensional DVR (or an already-built `SelfConsistentDVR`) to optimize the basis from
  - `opts`: `dict`
    > extra options, overriding the base DVR's own stored options, forwarded to the constructor
  - `:returns`: `PotentialOptimizedDVR`
    > the potential-optimized DVR


<a id="Psience.DVR.Extensions.PotentialOptimizedDVR.from_scf" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_scf(cls, scf_dvr: 'DirectProductDVR|SelfConsistentDVR', wfns=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L247)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L247?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a `PotentialOptimizedDVR` using the wavefunctions from a fully-converged SCF run (or an explicitly supplied set of wavefunctions), wrapping `scf_dvr` in a `SelfConsistentDVR` first if it isn't already one.
  - `scf_dvr`: `DirectProductDVR | SelfConsistentDVR`
    > the base multi-dimensional DVR (or an already-built `SelfConsistentDVR`) to run/use for the optimized basis
  - `wfns`: `Iterable[DVRWavefunctions] | None`
    > explicit wavefunctions to use instead of running the SCF procedure
  - `opts`: `dict`
    > extra options, overriding the base DVR's own stored options, forwarded to the constructor
  - `:returns`: `PotentialOptimizedDVR`
    > the potential-optimized DVR
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/DVR/Extensions/PotentialOptimizedDVR.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/DVR/Extensions/PotentialOptimizedDVR.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/DVR/Extensions/PotentialOptimizedDVR.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/DVR/Extensions/PotentialOptimizedDVR.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/DVR/Extensions.py#L199?message=Update%20Docs)   
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