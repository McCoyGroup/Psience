## <a id="Psience.Molecools.CoordinateSystems.ModeEmbedding">ModeEmbedding</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems.py#L1359)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems.py#L1359?message=Update%20Docs)]
</div>

Provides a specialization on a `MoleculaEmbedding` to express all properties
in terms of the attendant normal modes







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.Molecools.CoordinateSystems.ModeEmbedding.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, embedding: Psience.Molecools.CoordinateSystems.MolecularEmbedding, modes, mass_weight=None, dimensionless=None, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems.py#L1365)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems.py#L1365?message=Update%20Docs)]
</div>
**LLM Docstring**

Set up a mode-aware specialization of a `MolecularEmbedding`: resolves whatever form `modes` was given in (a manager, a `MolecularVibrations`, or a raw normal-modes object) down to a plain normal-modes object, optionally converting it to a dimensionless or mass-weighted basis, and records whether the resulting modes are mass-weighted.
  - `embedding`: `MolecularEmbedding`
    > the underlying molecular embedding to specialize
  - `modes`: `object`
    > the normal modes (or a manager/vibrations object wrapping them) to express properties in terms of
  - `mass_weight`: `bool | None`
    > whether to convert the modes to a mass-weighted basis
  - `dimensionless`: `bool | None`
    > whether to convert the modes to a dimensionless basis
  - `masses`: `np.ndarray | None`
    > masses to use for the mass-weighting/dimensionless conversions; defaults to the embedding's own masses
  - `:returns`: `None`
    > None


<a id="Psience.Molecools.CoordinateSystems.ModeEmbedding.mw_conversion" class="docs-object-method">&nbsp;</a> 
```python
mw_conversion(self, strip_dummies=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/ModeEmbedding.py#L1411)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/ModeEmbedding.py#L1411?message=Update%20Docs)]
</div>
**LLM Docstring**

Build the diagonal mass-weighting matrix (`sqrt(mass)` per Cartesian coordinate) used to convert plain Cartesian derivatives into mass-weighted ones.
  - `strip_dummies`: `bool | None`
    > whether to exclude dummy (non-positive-mass) atoms from the mass vector
  - `:returns`: `np.ndarray`
    > the diagonal mass-weighting matrix


<a id="Psience.Molecools.CoordinateSystems.ModeEmbedding.mw_inverse" class="docs-object-method">&nbsp;</a> 
```python
mw_inverse(self, strip_dummies=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/ModeEmbedding.py#L1432)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/ModeEmbedding.py#L1432?message=Update%20Docs)]
</div>
**LLM Docstring**

Build the diagonal inverse-mass-weighting matrix (`1/sqrt(mass)` per Cartesian coordinate) used to convert mass-weighted Cartesian derivatives back into plain ones.
  - `strip_dummies`: `bool | None`
    > whether to exclude dummy (non-positive-mass) atoms from the mass vector
  - `:returns`: `np.ndarray`
    > the diagonal inverse-mass-weighting matrix


<a id="Psience.Molecools.CoordinateSystems.ModeEmbedding.get_mw_cartesians_by_internals" class="docs-object-method">&nbsp;</a> 
```python
get_mw_cartesians_by_internals(self, order=None, mass_weighted=None, coords=None, strip_embedding=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/ModeEmbedding.py#L1454)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/ModeEmbedding.py#L1454?message=Update%20Docs)]
</div>
**LLM Docstring**

Fetch the Cartesians-by-internals Jacobian expansion, converted to mass-weighted Cartesians if `mass_weighted` (or `self.mass_weighted` by default) is set.
  - `order`: `int | None`
    > the highest derivative order to compute
  - `mass_weighted`: `bool | None`
    > whether to mass-weight the result; defaults to `self.mass_weighted`
  - `coords`: `np.ndarray | None`
    > alternate coordinates to compute the Jacobian at
  - `strip_embedding`: `bool`
    > whether to strip the fixed embedding coordinates
  - `:returns`: `list[np.ndarray]`
    > the (optionally mass-weighted) Cartesians-by-internals Jacobian tensors


<a id="Psience.Molecools.CoordinateSystems.ModeEmbedding.get_internals_by_mw_cartesians" class="docs-object-method">&nbsp;</a> 
```python
get_internals_by_mw_cartesians(self, order=None, mass_weighted=None, coords=None, strip_embedding=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/ModeEmbedding.py#L1485)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/ModeEmbedding.py#L1485?message=Update%20Docs)]
</div>
**LLM Docstring**

Fetch the internals-by-Cartesians Jacobian expansion, converted to be with respect to mass-weighted Cartesians if `mass_weighted` (or `self.mass_weighted` by default) is set.
  - `order`: `int | None`
    > the highest derivative order to compute
  - `mass_weighted`: `bool | None`
    > whether to mass-weight the result; defaults to `self.mass_weighted`
  - `coords`: `np.ndarray | None`
    > alternate coordinates to compute the Jacobian at
  - `strip_embedding`: `bool`
    > whether to strip the fixed embedding coordinates
  - `:returns`: `list[np.ndarray]`
    > the (optionally mass-weighted) internals-by-Cartesians Jacobian tensors


<a id="Psience.Molecools.CoordinateSystems.ModeEmbedding.get_internals_by_cartesians" class="docs-object-method">&nbsp;</a> 
```python
get_internals_by_cartesians(self, order=None, coords=None, strip_embedding=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/ModeEmbedding.py#L1526)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/ModeEmbedding.py#L1526?message=Update%20Docs)]
</div>
expresses raw internals or modes (internals or Cartesian) in terms of mass-weighted Cartesians
  - `order`: `Any`
    > 
  - `strip_embedding`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.CoordinateSystems.ModeEmbedding.get_cartesians_by_internals" class="docs-object-method">&nbsp;</a> 
```python
get_cartesians_by_internals(self, order=None, coords=None, strip_embedding=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/ModeEmbedding.py#L1578)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/ModeEmbedding.py#L1578?message=Update%20Docs)]
</div>
expresses raw internals or modes (internals or Cartesian) in terms of mass-weighted Cartesians
  - `order`: `Any`
    > 
  - `strip_embedding`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.CoordinateSystems.ModeEmbedding.get_inertia_tensor_expansion" class="docs-object-method">&nbsp;</a> 
```python
get_inertia_tensor_expansion(self, order=None, strip_embedding=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/ModeEmbedding.py#L1633)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/ModeEmbedding.py#L1633?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the Taylor expansion of the inertia tensor in terms of this embedding's coordinates (internal coordinates or normal modes), by re-expanding the inertial-frame derivatives through the Cartesians-by-internals Jacobian.
  - `order`: `int | None`
    > the highest derivative order to compute
  - `strip_embedding`: `bool`
    > whether to strip the fixed embedding coordinates from the underlying Jacobian
  - `:returns`: `list[np.ndarray]`
    > `[I0] + [derivative terms...]`, the inertia tensor and its derivatives


<a id="Psience.Molecools.CoordinateSystems.ModeEmbedding.get_inertial_frame" class="docs-object-method">&nbsp;</a> 
```python
get_inertial_frame(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/ModeEmbedding.py#L1650)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/ModeEmbedding.py#L1650?message=Update%20Docs)]
</div>
**LLM Docstring**

The molecule's inertial (principal-axis) frame, delegated to the underlying embedding.
  - `:returns`: `tuple`
    > the inertial frame, as returned by `MolecularEmbedding.inertial_frame`


<a id="Psience.Molecools.CoordinateSystems.ModeEmbedding.get_modes_by_coords" class="docs-object-method">&nbsp;</a> 
```python
get_modes_by_coords(self, mass_weighted=None, frequency_scaled=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/ModeEmbedding.py#L1661)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/ModeEmbedding.py#L1661?message=Update%20Docs)]
</div>
**LLM Docstring**

The modes-by-coordinates transformation matrix, optionally adjusting the mass-weighting/frequency-scaling convention of the modes first, and (if the modes are Cartesian but expressed relative to an internal-coordinate embedding) re-expressing them in terms of internal coordinates.
  - `mass_weighted`: `bool | None`
    > `True`/`False` to force mass-weighting on/off before extracting the matrix; `None` to leave the modes' current convention
  - `frequency_scaled`: `bool | None`
    > `True`/`False` to adjust frequency scaling before extracting the matrix (both branches currently call `remove_frequency_scaling`); `None` to leave it unchanged
  - `:returns`: `np.ndarray | None`
    > the modes-by-coordinates matrix, or `None` if no modes are set


<a id="Psience.Molecools.CoordinateSystems.ModeEmbedding.get_coords_by_modes" class="docs-object-method">&nbsp;</a> 
```python
get_coords_by_modes(self, mass_weighted=None, frequency_scaled=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/ModeEmbedding.py#L1691)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/ModeEmbedding.py#L1691?message=Update%20Docs)]
</div>
**LLM Docstring**

The coordinates-by-modes transformation matrix, optionally adjusting the mass-weighting/frequency-scaling convention of the modes first, and (if the modes are Cartesian but expressed relative to an internal-coordinate embedding) re-expressing them in terms of internal coordinates.
  - `mass_weighted`: `bool | None`
    > `True`/`False` to force mass-weighting on/off before extracting the matrix; `None` to leave the modes' current convention
  - `frequency_scaled`: `bool | None`
    > `True`/`False` to adjust frequency scaling before extracting the matrix (both branches currently call `remove_frequency_scaling`); `None` to leave it unchanged
  - `:returns`: `np.ndarray | None`
    > the coordinates-by-modes matrix, or `None` if no modes are set
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Molecools/CoordinateSystems/ModeEmbedding.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Molecools/CoordinateSystems/ModeEmbedding.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Molecools/CoordinateSystems/ModeEmbedding.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Molecools/CoordinateSystems/ModeEmbedding.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems.py#L1359?message=Update%20Docs)   
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