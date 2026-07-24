## <a id="Psience.Molecools.CoordinateSystems.MolecularCartesianCoordinateSystem">MolecularCartesianCoordinateSystem</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems.py#L2106)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems.py#L2106?message=Update%20Docs)]
</div>

Mirrors the standard Cartesian coordinate system in _almost_ all regards, but forces an embedding







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
name: str
```
<a id="Psience.Molecools.CoordinateSystems.MolecularCartesianCoordinateSystem.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, masses, coords, dummy_positions=None, converter_options=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems.py#L2111)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems.py#L2111?message=Update%20Docs)]
</div>

  - `molecule`: `AbstractMolecule`
    > 
  - `converter_options`: `Any`
    > 
  - `opts`: `Any`
    >


<a id="Psience.Molecools.CoordinateSystems.MolecularCartesianCoordinateSystem.to_state" class="docs-object-method">&nbsp;</a> 
```python
to_state(self, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularCartesianCoordinateSystem.py#L2133)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularCartesianCoordinateSystem.py#L2133?message=Update%20Docs)]
</div>
**LLM Docstring**

Serialize this coordinate system's state, adding the masses, coordinates, and dummy-atom positions on top of whatever the base class's `to_state` produces.
  - `serializer`: `object | None`
    > the serializer to use, forwarded to the base class
  - `:returns`: `dict`
    > the serialized state dict


<a id="Psience.Molecools.CoordinateSystems.MolecularCartesianCoordinateSystem.from_state" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_state(cls, data, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L2149)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L2149?message=Update%20Docs)]
</div>
**LLM Docstring**

Reconstruct a `MolecularCartesianCoordinateSystem` from a previously serialized state dict.
  - `data`: `dict`
    > the serialized state, as produced by `to_state`
  - `serializer`: `object | None`
    > the serializer to use, accepted for interface consistency but not used in this method's body
  - `:returns`: `MolecularCartesianCoordinateSystem`
    > the reconstructed coordinate system


<a id="Psience.Molecools.CoordinateSystems.MolecularCartesianCoordinateSystem.pre_convert_to" class="docs-object-method">&nbsp;</a> 
```python
pre_convert_to(self, system, opts=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularCartesianCoordinateSystem.py#L2171)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularCartesianCoordinateSystem.py#L2171?message=Update%20Docs)]
</div>
**LLM Docstring**

Ensure the masses are up to date in `converter_options`, and, if converting to a Z-matrix-family system, re-establish the embedding options (via `set_embedding`) before delegating to the base class's `pre_convert_to`.
  - `system`: `object`
    > the coordinate system being converted to
  - `opts`: `dict | None`
    > explicit conversion options to use instead of `self.converter_options`
  - `:returns`: `dict`
    > the resolved conversion options


<a id="Psience.Molecools.CoordinateSystems.MolecularCartesianCoordinateSystem.set_embedding" class="docs-object-method">&nbsp;</a> 
```python
set_embedding(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularCartesianCoordinateSystem.py#L2194)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularCartesianCoordinateSystem.py#L2194?message=Update%20Docs)]
</div>
Sets up the embedding options...
  - `:returns`: `_`
    >


<a id="Psience.Molecools.CoordinateSystems.MolecularCartesianCoordinateSystem.jacobian" class="docs-object-method">&nbsp;</a> 
```python
jacobian(self, coords, system, order=None, strip_dummies=None, converter_options=None, analytic_deriv_order=None, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularCartesianCoordinateSystem.py#L2225)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularCartesianCoordinateSystem.py#L2225?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the Jacobian of these Cartesian coordinates with respect to `system`, resolving the analytic-derivative order (defaulting to purely numerical for Z-matrix targets) and optionally excluding dummy-atom coordinates.
  - `coords`: `np.ndarray`
    > the coordinates to compute the Jacobian at
  - `system`: `object`
    > the target coordinate system
  - `order`: `int | list[int] | None`
    > the derivative order(s) to compute
  - `strip_dummies`: `bool | None`
    > whether to exclude dummy-atom coordinates; falls back to the converter options, then defaults to `False`
  - `converter_options`: `dict | None`
    > extra converter options merged with `self.converter_options`
  - `analytic_deriv_order`: `int | None`
    > the order up to which to compute the Jacobian analytically rather than numerically; falls back to the converter options, then defaults based on whether `system` is a Z-matrix
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Molecools/CoordinateSystems/MolecularCartesianCoordinateSystem.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Molecools/CoordinateSystems/MolecularCartesianCoordinateSystem.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Molecools/CoordinateSystems/MolecularCartesianCoordinateSystem.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Molecools/CoordinateSystems/MolecularCartesianCoordinateSystem.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems.py#L2106?message=Update%20Docs)   
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