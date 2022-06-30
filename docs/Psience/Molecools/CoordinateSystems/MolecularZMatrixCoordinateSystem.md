## <a id="Psience.Molecools.CoordinateSystems.MolecularZMatrixCoordinateSystem">MolecularZMatrixCoordinateSystem</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems.py#L76)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems.py#L76?message=Update%20Docs)]
</div>

Mirrors the standard ZMatrix coordinate system in _almost_ all regards, but forces an embedding

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

```python
name: str
embedding_coords: list
```
<a id="Psience.Molecools.CoordinateSystems.MolecularZMatrixCoordinateSystem.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, molecule, converter_options=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems.py#L82)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems.py#L82?message=Update%20Docs)]
</div>


- `molecule`: `AbstractMolecule`
    >No description...
- `converter_options`: `Any`
    >No description...
- `opts`: `Any`
    >No description...

<a id="Psience.Molecools.CoordinateSystems.MolecularZMatrixCoordinateSystem.origins" class="docs-object-method">&nbsp;</a> 
```python
@property
origins(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems.py#L?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.CoordinateSystems.MolecularZMatrixCoordinateSystem.axes" class="docs-object-method">&nbsp;</a> 
```python
@property
axes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems.py#L?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.CoordinateSystems.MolecularZMatrixCoordinateSystem.pre_convert" class="docs-object-method">&nbsp;</a> 
```python
pre_convert(self, system): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems.py#L106)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems.py#L106?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.CoordinateSystems.MolecularZMatrixCoordinateSystem.set_embedding" class="docs-object-method">&nbsp;</a> 
```python
set_embedding(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems.py#L110)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems.py#L110?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.CoordinateSystems.MolecularZMatrixCoordinateSystem.jacobian" class="docs-object-method">&nbsp;</a> 
```python
jacobian(self, *args, reembed=None, strip_dummies=None, converter_options=None, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems.py#L134)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems.py#L134?message=Update%20Docs)]
</div>

 </div>
</div>




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems.py#L76?message=Update%20Docs)