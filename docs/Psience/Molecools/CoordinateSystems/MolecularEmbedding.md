## <a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding">MolecularEmbedding</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems.py#L30)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems.py#L30?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
internal_fd_defaults: dict
cart_fd_defaults: dict
cartesian_by_internals_method: str
```
<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, masses, coords, internals, internal_fd_opts=None, cartesian_fd_opts=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems.py#L32)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems.py#L32?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.coords" class="docs-object-method">&nbsp;</a> 
```python
@property
coords(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L56)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L56?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.masses" class="docs-object-method">&nbsp;</a> 
```python
@property
masses(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L72)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L72?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.canonicalize_internal_coordinate_spec" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
canonicalize_internal_coordinate_spec(cls, spec): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L92)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L92?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.internals" class="docs-object-method">&nbsp;</a> 
```python
@property
internals(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L142)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L142?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.zmatrix" class="docs-object-method">&nbsp;</a> 
```python
@property
zmatrix(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L152)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L152?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.convert_to_internals" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
convert_to_internals(cls, coords, masses, spec): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L169)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L169?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.internal_coordinates_from_spec" class="docs-object-method">&nbsp;</a> 
```python
internal_coordinates_from_spec(self, spec: dict): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L217)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L217?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.internal_coordinates" class="docs-object-method">&nbsp;</a> 
```python
@property
internal_coordinates(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L224)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L224?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.strip_embedding_coordinates" class="docs-object-method">&nbsp;</a> 
```python
strip_embedding_coordinates(self, coords): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L241)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L241?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.strip_derivative_embedding" class="docs-object-method">&nbsp;</a> 
```python
strip_derivative_embedding(self, derivs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L259)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L259?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.restore_embedding_coordinates" class="docs-object-method">&nbsp;</a> 
```python
restore_embedding_coordinates(self, coords): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L272)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L272?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.restore_derivative_embedding" class="docs-object-method">&nbsp;</a> 
```python
restore_derivative_embedding(self, derivs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L301)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L301?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.get_internals" class="docs-object-method">&nbsp;</a> 
```python
get_internals(self, *, coords=None, strip_embedding=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L320)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L320?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.redundant_internal_transformation" class="docs-object-method">&nbsp;</a> 
```python
@property
redundant_internal_transformation(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L332)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L332?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.embedding_coords" class="docs-object-method">&nbsp;</a> 
```python
@property
embedding_coords(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L611)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L611?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.get_cartesians_by_internals" class="docs-object-method">&nbsp;</a> 
```python
get_cartesians_by_internals(self, order=None, strip_embedding=False, reembed=True, method=None, coords=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L625)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L625?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.get_internals_by_cartesians" class="docs-object-method">&nbsp;</a> 
```python
get_internals_by_cartesians(self, order=None, strip_embedding=False, coords=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L692)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L692?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.embed_coords" class="docs-object-method">&nbsp;</a> 
```python
embed_coords(self, coords, sel=None, in_paf=False, planar_ref_tolerance=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L721)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L721?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.unembed_derivs" class="docs-object-method">&nbsp;</a> 
```python
unembed_derivs(self, coords, derivs, sel=None, in_paf=False, planar_ref_tolerance=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L731)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L731?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.inertia_tensor" class="docs-object-method">&nbsp;</a> 
```python
@property
inertia_tensor(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L749)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L749?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.inertial_frame" class="docs-object-method">&nbsp;</a> 
```python
@property
inertial_frame(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L756)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L756?message=Update%20Docs)]
</div>
Provides the inertial axis frame
  - `:returns`: `_`
    >


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.inertial_frame_derivatives" class="docs-object-method">&nbsp;</a> 
```python
inertial_frame_derivatives(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L778)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L778?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.translation_rotation_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
translation_rotation_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L785)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L785?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.MolecularEmbedding.get_translation_rotation_invariant_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_translation_rotation_invariant_transformation(self, order=0, mass_weighted=True, strip_embedding=True, coords=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L796)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularEmbedding.py#L796?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Molecools/CoordinateSystems/MolecularEmbedding.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Molecools/CoordinateSystems/MolecularEmbedding.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Molecools/CoordinateSystems/MolecularEmbedding.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Molecools/CoordinateSystems/MolecularEmbedding.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems.py#L30?message=Update%20Docs)   
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