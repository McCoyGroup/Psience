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
__init__(self, normal_modes: Psience.Modes.NormalModes.NormalModes, transformation, inverse=None, origin=None, masses=None, freqs=None, mass_weighted=None, frequency_scaled=None, **etc): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes.py#L15)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes.py#L15?message=Update%20Docs)]
</div>


<a id="Psience.Modes.LocalizedModes.LocalizedModes.mass_weighted" class="docs-object-method">&nbsp;</a> 
```python
@property
mass_weighted(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L50)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L50?message=Update%20Docs)]
</div>


<a id="Psience.Modes.LocalizedModes.LocalizedModes.frequency_scaled" class="docs-object-method">&nbsp;</a> 
```python
@property
frequency_scaled(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L60)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L60?message=Update%20Docs)]
</div>


<a id="Psience.Modes.LocalizedModes.LocalizedModes.g_matrix" class="docs-object-method">&nbsp;</a> 
```python
@property
g_matrix(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L70)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L70?message=Update%20Docs)]
</div>


<a id="Psience.Modes.LocalizedModes.LocalizedModes.modify" class="docs-object-method">&nbsp;</a> 
```python
modify(self, base_modes=None, *, transformation=None, freqs=None, origin=None, masses=None, inverse=None, name=None, mass_weighted=None, frequency_scaled=None, g_matrix=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L77)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L77?message=Update%20Docs)]
</div>


<a id="Psience.Modes.LocalizedModes.LocalizedModes.make_mass_weighted" class="docs-object-method">&nbsp;</a> 
```python
make_mass_weighted(self, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L101)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L101?message=Update%20Docs)]
</div>


<a id="Psience.Modes.LocalizedModes.LocalizedModes.remove_mass_weighting" class="docs-object-method">&nbsp;</a> 
```python
remove_mass_weighting(self, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L103)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L103?message=Update%20Docs)]
</div>


<a id="Psience.Modes.LocalizedModes.LocalizedModes.make_frequency_scaled" class="docs-object-method">&nbsp;</a> 
```python
make_frequency_scaled(self, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L105)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L105?message=Update%20Docs)]
</div>


<a id="Psience.Modes.LocalizedModes.LocalizedModes.remove_frequency_scaling" class="docs-object-method">&nbsp;</a> 
```python
remove_frequency_scaling(self, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L107)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L107?message=Update%20Docs)]
</div>


<a id="Psience.Modes.LocalizedModes.LocalizedModes.local_hessian" class="docs-object-method">&nbsp;</a> 
```python
@property
local_hessian(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L110)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L110?message=Update%20Docs)]
</div>


<a id="Psience.Modes.LocalizedModes.LocalizedModes.localize" class="docs-object-method">&nbsp;</a> 
```python
localize(self, method=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L118)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L118?message=Update%20Docs)]
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