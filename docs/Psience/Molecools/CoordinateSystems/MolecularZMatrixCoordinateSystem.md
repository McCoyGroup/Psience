## <a id="Psience.Molecools.CoordinateSystems.MolecularZMatrixCoordinateSystem">MolecularZMatrixCoordinateSystem</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/CoordinateSystems.py#L1129)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/CoordinateSystems.py#L1129?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.py#L1135)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.py#L1135?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.py#L1156)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.py#L1156?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.MolecularZMatrixCoordinateSystem.axes" class="docs-object-method">&nbsp;</a> 
```python
@property
axes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.py#L1159)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.py#L1159?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.MolecularZMatrixCoordinateSystem.pre_convert" class="docs-object-method">&nbsp;</a> 
```python
pre_convert(self, system): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.py#L1163)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.py#L1163?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.MolecularZMatrixCoordinateSystem.set_embedding" class="docs-object-method">&nbsp;</a> 
```python
set_embedding(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.py#L1167)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.py#L1167?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.CoordinateSystems.MolecularZMatrixCoordinateSystem.jacobian" class="docs-object-method">&nbsp;</a> 
```python
jacobian(self, coords, *args, reembed=None, strip_dummies=None, converter_options=None, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.py#L1195)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.py#L1195?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/CoordinateSystems.py#L1129?message=Update%20Docs)   
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