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
**LLM Docstring**

Compute the first and second derivatives of the (mass-weighted) inertia tensor with respect to mass-weighted Cartesian displacements, using closed-form tensor expressions rather than finite differences. Handles batched inputs (multiple geometries at once) and can restrict to a subset of atoms via `sel`.
  - `crds`: `np.ndarray`
    > Cartesian coordinates of the atoms, with an optional leading batch dimension
  - `mass`: `np.ndarray`
    > atomic masses (only positive entries are treated as real atoms), with an optional leading batch dimension
  - `sel`: `np.ndarray | None`
    > optional subset of atom indices to restrict the calculation to; intersected with the real (mass > 0) atoms
  - `:returns`: `list[np.ndarray]`
    > `[I0Y, I0YY]`, the first derivative tensor (shape `(..., 3*nAt, 3, 3)`) and second derivative tensor (shape `(..., 3*nAt, 3*nAt, 3, 3)`) of the inertia tensor with respect to mass-weighted Cartesian displacements


<a id="Psience.Molecools.Properties.StructuralProperties.get_prop_moments_of_inertia" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_prop_moments_of_inertia(cls, coords, masses): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L182)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L182?message=Update%20Docs)]
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
get_prop_principle_axis_rotation(cls, coords, masses, sel=None, inverse=False, shift_sel=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L227)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L227?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L270)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L270?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L303)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L303?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L319)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L319?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L446)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L446?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L558)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L558?message=Update%20Docs)]
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
get_prop_eckart_transformation(cls, masses, ref, coords, sel=None, inverse=False, reset_com=False, in_paf=False, planar_ref_tolerance=None, proper_rotation=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L590)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L590?message=Update%20Docs)]
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
get_eckart_embedded_coords(cls, masses, ref, coords, reset_com=True, in_paf=False, sel=None, planar_ref_tolerance=None, proper_rotation=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L660)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L660?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L736)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L736?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L776)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L776?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the Coriolis zeta constants for a set of normal modes, by projecting the mode displacement vectors into the molecule's inertial (principal-axis) frame and forming the antisymmetric combination that gives the coupling between each pair of modes about each principal axis.
  - `carts`: `np.ndarray`
    > Cartesian coordinates of the atoms, with an optional leading batch dimension
  - `modes`: `np.ndarray`
    > the (mass-weighted) normal-mode displacement vectors
  - `masses`: `np.ndarray`
    > the atomic masses
  - `:returns`: `tuple[np.ndarray, tuple[np.ndarray, np.ndarray]]`
    > a `(zeta, (mom_i, eigs))` pair, where `zeta` has shape `(..., 3, nmodes, nmodes)` giving the zeta constant about each of the 3 principal axes for every mode pair, `mom_i` is the moments of inertia, and `eigs` is the principal-axis frame
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