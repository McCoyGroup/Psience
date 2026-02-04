## <a id="Psience.Molecools.Properties.StructuralProperties">StructuralProperties</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties.py#L38)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties.py#L38?message=Update%20Docs)]
</div>

The set of molecular properties that depend on its coordinates/configuration.
Slowly trying to move code out of this and into numputils/Hamiltonian/Evaluator







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
planar_ref_tolerance: float
EmbeddingData: PrincipleAxisData
EckartData: EckartData
```
<a id="Psience.Molecools.Properties.StructuralProperties.get_prop_mass_weighted_coords" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_prop_mass_weighted_coords(cls, coords, masses): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L44)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L44?message=Update%20Docs)]
</div>
Gets the mass-weighted coordinates for the system
  - `coords`: `CoordinateSet`
    > 
  - `masses`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.StructuralProperties.get_prop_center_of_mass" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_prop_center_of_mass(cls, coords, masses): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L63)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L63?message=Update%20Docs)]
</div>
Gets the center of mass for the coordinates
  - `coords`: `CoordinateSet`
    > 
  - `masses`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.StructuralProperties.get_prop_inertia_tensors" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_prop_inertia_tensors(cls, coords, masses): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L80)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L80?message=Update%20Docs)]
</div>
Computes the moment of intertia tensors for the walkers with coordinates coords (assumes all have the same masses)
  - `coords`: `CoordinateSet`
    > 
  - `masses`: `np.ndarray`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.StructuralProperties.get_prop_inertial_frame_derivatives" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_prop_inertial_frame_derivatives(cls, crds, mass, sel=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L109)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L109?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.StructuralProperties.get_prop_moments_of_inertia" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_prop_moments_of_inertia(cls, coords, masses): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L168)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L168?message=Update%20Docs)]
</div>
Computes the moment of inertia tensor for the walkers with coordinates coords (assumes all have the same masses)
  - `coords`: `CoordinateSet`
    > 
  - `masses`: `np.ndarray`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.StructuralProperties.get_prop_principle_axis_rotation" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_prop_principle_axis_rotation(cls, coords, masses, sel=None, inverse=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L213)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L213?message=Update%20Docs)]
</div>
Generates the principle axis transformation for a set of coordinates and positions
  - `coords`: `CoordinateSet`
    > 
  - `masses`: `np.ndarray`
    > 
  - `:returns`: `MolecularTransformation | List[MolecularTransformation]`
    >


<a id="Psience.Molecools.Properties.StructuralProperties.get_principle_axis_embedded_coords" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_principle_axis_embedded_coords(cls, coords, masses, sel=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L253)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L253?message=Update%20Docs)]
</div>
Returns coordinate embedded in the principle axis frame
  - `coords`: `Any`
    > 
  - `masses`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.StructuralProperties.get_prop_principle_axis_data" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_prop_principle_axis_data(cls, coords, masses): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L286)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L286?message=Update%20Docs)]
</div>
Generates the principle axis transformation for a set of coordinates and positions
  - `coords`: `CoordinateSet`
    > 
  - `masses`: `np.ndarray`
    > 
  - `:returns`: `MolecularTransformation | List[MolecularTransformation]`
    >


<a id="Psience.Molecools.Properties.StructuralProperties.get_prop_translation_rotation_eigenvectors" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_prop_translation_rotation_eigenvectors(cls, coords, masses): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L302)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L302?message=Update%20Docs)]
</div>
Returns the eigenvectors corresponding to translations and rotations
in the system
  - `coords`: `Any`
    > 
  - `masses`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.StructuralProperties.get_eckart_rotations" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_eckart_rotations(cls, masses, ref, coords, sel=None, in_paf=False, planar_ref_tolerance=None, proper_rotation=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L429)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L429?message=Update%20Docs)]
</div>
Generates the Eckart rotation that will align ref and coords, assuming initially that `ref` and `coords` are
in the principle axis frame
  - `masses`: `Any`
    > 
  - `ref`: `Any`
    > 
  - `coords`: `np.ndarray`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.StructuralProperties.get_eckart_embedding_data" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_eckart_embedding_data(cls, masses, ref, coords, sel=None, in_paf=False, planar_ref_tolerance=None, proper_rotation=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L541)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L541?message=Update%20Docs)]
</div>
Embeds a set of coordinates in the reference frame
  - `masses`: `np.ndarray`
    > 
  - `ref`: `CoordinateSet`
    > 
  - `coords`: `CoordinateSet`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.StructuralProperties.get_prop_eckart_transformation" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_prop_eckart_transformation(cls, masses, ref, coords, sel=None, inverse=False, reset_com=False, planar_ref_tolerance=None, proper_rotation=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L564)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L564?message=Update%20Docs)]
</div>
Computes Eckart transformations for a set of coordinates
  - `masses`: `np.ndarray`
    > 
  - `ref`: `CoordinateSet`
    > 
  - `coords`: `CoordinateSet`
    > 
  - `:returns`: `MolecularTransformation | List[MolecularTransformation]`
    >


<a id="Psience.Molecools.Properties.StructuralProperties.get_eckart_embedded_coords" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_eckart_embedded_coords(cls, masses, ref, coords, reset_com=False, in_paf=False, sel=None, planar_ref_tolerance=None, proper_rotation=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L618)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L618?message=Update%20Docs)]
</div>
Embeds a set of coordinates in the reference frame
  - `masses`: `np.ndarray`
    > 
  - `ref`: `CoordinateSet`
    > 
  - `coords`: `CoordinateSet`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.StructuralProperties.get_prop_g_matrix" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_prop_g_matrix(cls, masses, coords, internal_coords): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L679)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L679?message=Update%20Docs)]
</div>
Gets the molecular g-matrix
  - `masses`: `np.ndarray`
    > 
  - `coords`: `CoordinateSet`
    > 
  - `internal_coords`: `CoordinateSet`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.StructuralProperties.get_prop_coriolis_constants" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_prop_coriolis_constants(cls, carts, modes, masses): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L719)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L719?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Molecools/Properties/StructuralProperties.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Molecools/Properties/StructuralProperties.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Molecools/Properties/StructuralProperties.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Molecools/Properties/StructuralProperties.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties.py#L38?message=Update%20Docs)   
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