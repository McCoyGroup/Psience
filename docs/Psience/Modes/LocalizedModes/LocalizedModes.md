## <a id="Psience.Modes.LocalizedModes.LocalizedModes">LocalizedModes</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes.py#L13)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes.py#L13?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.Modes.LocalizedModes.LocalizedModes.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, normal_modes: Psience.Modes.MixtureModes.MixtureModes, transformation, inverse=None, origin=None, masses=None, freqs=None, mass_weighted=None, frequency_scaled=None, **etc): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes.py#L15)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes.py#L15?message=Update%20Docs)]
</div>


<a id="Psience.Modes.LocalizedModes.LocalizedModes.freqs" class="docs-object-method">&nbsp;</a> 
```python
@property
freqs(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L58)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L58?message=Update%20Docs)]
</div>


<a id="Psience.Modes.LocalizedModes.LocalizedModes.mass_weighted" class="docs-object-method">&nbsp;</a> 
```python
@property
mass_weighted(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L70)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L70?message=Update%20Docs)]
</div>


<a id="Psience.Modes.LocalizedModes.LocalizedModes.g_matrix" class="docs-object-method">&nbsp;</a> 
```python
@property
g_matrix(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L90)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L90?message=Update%20Docs)]
</div>


<a id="Psience.Modes.LocalizedModes.LocalizedModes.modify" class="docs-object-method">&nbsp;</a> 
```python
modify(self, base_modes=None, *, transformation=None, freqs=None, origin=None, masses=None, inverse=None, name=None, mass_weighted=None, frequency_scaled=None, g_matrix=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L106)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L106?message=Update%20Docs)]
</div>


<a id="Psience.Modes.LocalizedModes.LocalizedModes.make_mass_weighted" class="docs-object-method">&nbsp;</a> 
```python
make_mass_weighted(self, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L142)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L142?message=Update%20Docs)]
</div>


<a id="Psience.Modes.LocalizedModes.LocalizedModes.remove_mass_weighting" class="docs-object-method">&nbsp;</a> 
```python
remove_mass_weighting(self, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L144)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L144?message=Update%20Docs)]
</div>


<a id="Psience.Modes.LocalizedModes.LocalizedModes.make_frequency_scaled" class="docs-object-method">&nbsp;</a> 
```python
make_frequency_scaled(self, freqs=None, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L146)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L146?message=Update%20Docs)]
</div>


<a id="Psience.Modes.LocalizedModes.LocalizedModes.remove_frequency_scaling" class="docs-object-method">&nbsp;</a> 
```python
remove_frequency_scaling(self, freqs=None, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L155)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L155?message=Update%20Docs)]
</div>


<a id="Psience.Modes.LocalizedModes.LocalizedModes.compute_hessian" class="docs-object-method">&nbsp;</a> 
```python
compute_hessian(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L173)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L173?message=Update%20Docs)]
</div>


<a id="Psience.Modes.LocalizedModes.LocalizedModes.apply_transformation" class="docs-object-method">&nbsp;</a> 
```python
apply_transformation(self, transformation, inverse=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L183)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L183?message=Update%20Docs)]
</div>


<a id="Psience.Modes.LocalizedModes.LocalizedModes.get_complement" class="docs-object-method">&nbsp;</a> 
```python
get_complement(self, concatenate=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L206)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L206?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Modes/LocalizedModes/LocalizedModes.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Modes/LocalizedModes/LocalizedModes.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Modes/LocalizedModes/LocalizedModes.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Modes/LocalizedModes/LocalizedModes.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes.py#L13?message=Update%20Docs)   
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