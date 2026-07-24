## <a id="Psience.Molecools.CoordinateSystems.MolecularZMatrixCoordinateSystem">MolecularZMatrixCoordinateSystem</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems.py#L1876)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems.py#L1876?message=Update%20Docs)]
</div>

Mirrors the standard ZMatrix coordinate system in _almost_ all regards, but forces an embedding







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
name: str
embedding_coords: list
```
<a id="Psience.Molecools.CoordinateSystems.MolecularZMatrixCoordinateSystem.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, masses, coords, converter_options=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems.py#L1882)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems.py#L1882?message=Update%20Docs)]
</div>

  - `molecule`: `AbstractMolecule`
    > 
  - `converter_options`: `Any`
    > 
  - `opts`: `Any`
    >


<a id="Psience.Molecools.CoordinateSystems.MolecularZMatrixCoordinateSystem.origins" class="docs-object-method">&nbsp;</a> 
```python
@property
origins(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.py#L1903)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.py#L1903?message=Update%20Docs)]
</div>
**LLM Docstring**

The Z-matrix embedding's origin points (typically the reference center of mass), from `converter_options['origins']`.
  - `:returns`: `np.ndarray`
    > the origin points


<a id="Psience.Molecools.CoordinateSystems.MolecularZMatrixCoordinateSystem.axes" class="docs-object-method">&nbsp;</a> 
```python
@property
axes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.py#L1914)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.py#L1914?message=Update%20Docs)]
</div>
**LLM Docstring**

The Z-matrix embedding's reference axes (typically two principal axes), from `converter_options['axes']`.
  - `:returns`: `np.ndarray`
    > the reference axes


<a id="Psience.Molecools.CoordinateSystems.MolecularZMatrixCoordinateSystem.get_direct_converter" class="docs-object-method">&nbsp;</a> 
```python
get_direct_converter(self, target): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.py#L1926)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.py#L1926?message=Update%20Docs)]
</div>
**LLM Docstring**

Provide a converter from this molecular Z-matrix system directly to the plain (non-molecular) `ZMatrix` coordinate system, if `target` is one.
  - `target`: `object`
    > the coordinate system being converted to
  - `:returns`: `MolecularZMatrixToRegularZMatrixConverter | None`
    > a `MolecularZMatrixToRegularZMatrixConverter`, or `None` if `target` isn't a `ZMatrix` system


<a id="Psience.Molecools.CoordinateSystems.MolecularZMatrixCoordinateSystem.get_inverse_converter" class="docs-object-method">&nbsp;</a> 
```python
get_inverse_converter(self, target): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.py#L1940)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.py#L1940?message=Update%20Docs)]
</div>
**LLM Docstring**

Provide a converter from the plain `ZMatrix` coordinate system into this molecular Z-matrix system, if `target` is one.
  - `target`: `object`
    > the coordinate system being converted from
  - `:returns`: `RegularZMatrixToMolecularZMatrixConverter | None`
    > a `RegularZMatrixToMolecularZMatrixConverter`, or `None` if `target` isn't a `ZMatrix` system


<a id="Psience.Molecools.CoordinateSystems.MolecularZMatrixCoordinateSystem.pre_convert_to" class="docs-object-method">&nbsp;</a> 
```python
pre_convert_to(self, system, opts=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.py#L1955)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.py#L1955?message=Update%20Docs)]
</div>
**LLM Docstring**

Re-establish the embedding options (via `set_embedding`) before delegating to the base class's `pre_convert_to`, preserving the existing atom `ordering` if the caller supplied its own options dict.
  - `system`: `object`
    > the coordinate system being converted to
  - `opts`: `dict | None`
    > explicit conversion options to use instead of `self.converter_options`
  - `:returns`: `dict`
    > the resolved conversion options


<a id="Psience.Molecools.CoordinateSystems.MolecularZMatrixCoordinateSystem.set_embedding" class="docs-object-method">&nbsp;</a> 
```python
set_embedding(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.py#L1976)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.py#L1976?message=Update%20Docs)]
</div>
**LLM Docstring**

(Re)compute and store this Z-matrix system's embedding options -- the reference origin, reference axes (chosen via `_get_best_axes` to avoid ill-conditioned choices), axis labels, masses, dummy-atom positions, and reference coordinates -- based on the current center of mass and inertial axes, fixing up the Z-matrix `ordering`'s first three rows to reference the embedding's dummy origin/axis points if an ordering is present.
  - `:returns`: `None`
    > None


<a id="Psience.Molecools.CoordinateSystems.MolecularZMatrixCoordinateSystem.jacobian" class="docs-object-method">&nbsp;</a> 
```python
jacobian(self, coords, *args, reembed=None, strip_dummies=None, converter_options=None, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.py#L2012)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.py#L2012?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the Jacobian of this Z-matrix system with respect to Cartesian coordinates, handling batched/multi-frame inputs, optional dummy-atom stripping, and temporarily overriding the `reembed` converter option for the duration of the call.
  - `coords`: `np.ndarray`
    > the Cartesian coordinates to compute the Jacobian at
  - `args`: `tuple`
    > extra positional arguments forwarded to the base class's `jacobian`
  - `reembed`: `bool | None`
    > whether to re-embed (Eckart-align) during the Jacobian calculation; falls back to the converter options, then defaults to `True`
  - `strip_dummies`: `bool | None`
    > whether to exclude dummy-atom coordinates from the Jacobian; falls back to the converter options, then defaults to `False`
  - `converter_options`: `dict | None`
    > extra converter options merged with `self.converter_options`
  - `kwargs`: `dict`
    > extra keyword arguments forwarded to the base class's `jacobian`
  - `:returns`: `list[np.ndarray]`
    > the computed Jacobian tensor(s)
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems.py#L1876?message=Update%20Docs)   
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