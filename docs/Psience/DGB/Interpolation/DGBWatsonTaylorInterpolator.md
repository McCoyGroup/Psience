## <a id="Psience.DGB.Interpolation.DGBWatsonTaylorInterpolator">DGBWatsonTaylorInterpolator</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DGB/Interpolation.py#L202)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DGB/Interpolation.py#L202?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.DGB.Interpolation.DGBWatsonTaylorInterpolator.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, centers, potential_derivatives, modes, power=4, include_harmonic_basis=False, harmonic_distance_cutoff=None, pairwise_potential_functions=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DGB/Interpolation/DGBWatsonTaylorInterpolator.py#L203)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DGB/Interpolation/DGBWatsonTaylorInterpolator.py#L203?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Interpolation.DGBWatsonTaylorInterpolator.take_remainder_potential" class="docs-object-method">&nbsp;</a> 
```python
take_remainder_potential(self, centers, potential_derivatives, modes): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DGB/Interpolation/DGBWatsonTaylorInterpolator.py#L258)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DGB/Interpolation/DGBWatsonTaylorInterpolator.py#L258?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Interpolation.DGBWatsonTaylorInterpolator.take_ppf_remainder" class="docs-object-method">&nbsp;</a> 
```python
take_ppf_remainder(self, centers, potential_derivatives): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DGB/Interpolation/DGBWatsonTaylorInterpolator.py#L269)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DGB/Interpolation/DGBWatsonTaylorInterpolator.py#L269?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Interpolation.DGBWatsonTaylorInterpolator.taylor_interp" class="docs-object-method">&nbsp;</a> 
```python
taylor_interp(self, points, dists, neighbors, derivs, power=None, deriv_order=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DGB/Interpolation/DGBWatsonTaylorInterpolator.py#L277)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DGB/Interpolation/DGBWatsonTaylorInterpolator.py#L277?message=Update%20Docs)]
</div>


<a id="Psience.DGB.Interpolation.DGBWatsonTaylorInterpolator.__call__" class="docs-object-method">&nbsp;</a> 
```python
__call__(self, points, deriv_order=None, *, neighborhood_size=None, power=None, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/DGB/Interpolation/DGBWatsonTaylorInterpolator.py#L340)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/DGB/Interpolation/DGBWatsonTaylorInterpolator.py#L340?message=Update%20Docs)]
</div>
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/DGB/Interpolation/DGBWatsonTaylorInterpolator.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/DGB/Interpolation/DGBWatsonTaylorInterpolator.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/DGB/Interpolation/DGBWatsonTaylorInterpolator.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/DGB/Interpolation/DGBWatsonTaylorInterpolator.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/DGB/Interpolation.py#L202?message=Update%20Docs)   
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