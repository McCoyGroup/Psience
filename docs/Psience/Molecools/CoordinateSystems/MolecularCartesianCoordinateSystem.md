## <a id="Psience.Molecools.CoordinateSystems.MolecularCartesianCoordinateSystem">MolecularCartesianCoordinateSystem</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems.py#L1395)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems.py#L1395?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems.py#L1400)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems.py#L1400?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularCartesianCoordinateSystem.py#L1422)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularCartesianCoordinateSystem.py#L1422?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.MolecularCartesianCoordinateSystem.from_state" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_state(cls, data, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1428)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1428?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.MolecularCartesianCoordinateSystem.pre_convert" class="docs-object-method">&nbsp;</a> 
```python
pre_convert(self, system): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularCartesianCoordinateSystem.py#L1438)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularCartesianCoordinateSystem.py#L1438?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.MolecularCartesianCoordinateSystem.set_embedding" class="docs-object-method">&nbsp;</a> 
```python
set_embedding(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularCartesianCoordinateSystem.py#L1443)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularCartesianCoordinateSystem.py#L1443?message=Update%20Docs)]
</div>
Sets up the embedding options...
  - `:returns`: `_`
    >


<a id="Psience.Molecools.CoordinateSystems.MolecularCartesianCoordinateSystem.jacobian" class="docs-object-method">&nbsp;</a> 
```python
jacobian(self, coords, system, order=None, strip_dummies=None, converter_options=None, analytic_deriv_order=None, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems/MolecularCartesianCoordinateSystem.py#L1474)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems/MolecularCartesianCoordinateSystem.py#L1474?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Molecools/CoordinateSystems/MolecularCartesianCoordinateSystem.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Molecools/CoordinateSystems/MolecularCartesianCoordinateSystem.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Molecools/CoordinateSystems/MolecularCartesianCoordinateSystem.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Molecools/CoordinateSystems/MolecularCartesianCoordinateSystem.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems.py#L1395?message=Update%20Docs)   
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