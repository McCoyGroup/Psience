## <a id="Psience.Modes.MixtureModes.MixtureModes">MixtureModes</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes.py#L15)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes.py#L15?message=Update%20Docs)]
</div>

A `McUtils.Coordinerds.CoordinateSystem` object that expresses coordinates as
a rotation on some base set of coordinates with some associated frequencies.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
name: str
```
<a id="Psience.Modes.MixtureModes.MixtureModes.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, basis, coeffs, freqs=None, origin=None, masses=None, inverse=None, mass_weighted=False, frequency_scaled=False, g_matrix=None, name=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes.py#L21)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes.py#L21?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.prep_modes" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
prep_modes(cls, modes): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L58)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L58?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L110)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L110?message=Update%20Docs)]
</div>
Takes a slice of the modes
  - `item`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Modes.MixtureModes.MixtureModes.modify" class="docs-object-method">&nbsp;</a> 
```python
modify(self, matrix=None, *, freqs=None, origin=None, masses=None, inverse=None, name=None, mass_weighted=None, frequency_scaled=None, g_matrix=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L135)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L135?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.rotate" class="docs-object-method">&nbsp;</a> 
```python
rotate(self, rot, in_place=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L160)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L160?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.transform" class="docs-object-method">&nbsp;</a> 
```python
transform(self, tf, inv=None, origin=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L163)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L163?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.cartesian_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
cartesian_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L192)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L192?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.embed_coords" class="docs-object-method">&nbsp;</a> 
```python
embed_coords(self, carts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L196)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L196?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.unembed_coords" class="docs-object-method">&nbsp;</a> 
```python
unembed_coords(self, mode_coords): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L202)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L202?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.total_transformation" class="docs-object-method">&nbsp;</a> 
```python
@property
total_transformation(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L210)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L210?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.inverse_transformation" class="docs-object-method">&nbsp;</a> 
```python
@property
inverse_transformation(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L213)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L213?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.embed_derivs" class="docs-object-method">&nbsp;</a> 
```python
embed_derivs(self, derivs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L220)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L220?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.unembed_derivs" class="docs-object-method">&nbsp;</a> 
```python
unembed_derivs(self, derivs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L222)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L222?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.is_cartesian" class="docs-object-method">&nbsp;</a> 
```python
@property
is_cartesian(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L225)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L225?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.coords_by_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
coords_by_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L231)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L231?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.modes_by_coords" class="docs-object-method">&nbsp;</a> 
```python
@property
modes_by_coords(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L234)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L234?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.compute_local_hessian" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
compute_local_hessian(cls, f, g): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L266)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L266?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.compute_hessian" class="docs-object-method">&nbsp;</a> 
```python
compute_hessian(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L271)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L271?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.compute_gmatrix" class="docs-object-method">&nbsp;</a> 
```python
compute_gmatrix(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L274)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L274?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.local_hessian" class="docs-object-method">&nbsp;</a> 
```python
@property
local_hessian(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L282)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L282?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.local_freqs" class="docs-object-method">&nbsp;</a> 
```python
@property
local_freqs(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L289)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L289?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.make_mass_weighted" class="docs-object-method">&nbsp;</a> 
```python
make_mass_weighted(self, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L293)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L293?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.remove_mass_weighting" class="docs-object-method">&nbsp;</a> 
```python
remove_mass_weighting(self, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L306)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L306?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.make_frequency_scaled" class="docs-object-method">&nbsp;</a> 
```python
make_frequency_scaled(self, freqs=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L327)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L327?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.remove_frequency_scaling" class="docs-object-method">&nbsp;</a> 
```python
remove_frequency_scaling(self, freqs=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L339)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L339?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.make_dimensionless" class="docs-object-method">&nbsp;</a> 
```python
make_dimensionless(self, freqs=None, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L351)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L351?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.make_dimensioned" class="docs-object-method">&nbsp;</a> 
```python
make_dimensioned(self, freqs=None, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L358)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L358?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.apply_projection" class="docs-object-method">&nbsp;</a> 
```python
apply_projection(self, proj, project_transrot=True, masses=None, origin=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L378)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L378?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.apply_constraints" class="docs-object-method">&nbsp;</a> 
```python
apply_constraints(self, coordinate_constraints, atoms=None, masses=None, origin=None, orthogonal_projection=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L419)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L419?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Modes/MixtureModes/MixtureModes.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Modes/MixtureModes/MixtureModes.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Modes/MixtureModes/MixtureModes.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Modes/MixtureModes/MixtureModes.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes.py#L15?message=Update%20Docs)   
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