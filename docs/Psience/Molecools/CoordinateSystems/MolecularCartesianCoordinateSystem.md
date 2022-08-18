## <a id="Psience.Molecools.CoordinateSystems.MolecularCartesianCoordinateSystem">MolecularCartesianCoordinateSystem</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems.py#L186)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems.py#L186?message=Update%20Docs)]
</div>

Mirrors the standard Cartesian coordinate system in _almost_ all regards, but forces an embedding

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

```python
name: str
```
<a id="Psience.Molecools.CoordinateSystems.MolecularCartesianCoordinateSystem.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, molecule, converter_options=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems.py#L191)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems.py#L191?message=Update%20Docs)]
</div>


- `opts`: `Any`
    >
- `converter_options`: `Any`
    >
- `molecule`: `AbstractMolecule`
    >

<a id="Psience.Molecools.CoordinateSystems.MolecularCartesianCoordinateSystem.pre_convert" class="docs-object-method">&nbsp;</a> 
```python
pre_convert(self, system): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems.py#L208)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems.py#L208?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.CoordinateSystems.MolecularCartesianCoordinateSystem.set_embedding" class="docs-object-method">&nbsp;</a> 
```python
set_embedding(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems.py#L212)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems.py#L212?message=Update%20Docs)]
</div>

Sets up the embedding options...
- `:returns`: `_`
    >

<a id="Psience.Molecools.CoordinateSystems.MolecularCartesianCoordinateSystem.jacobian" class="docs-object-method">&nbsp;</a> 
```python
jacobian(self, coords, system, strip_dummies=None, converter_options=None, analytic_deriv_order=None, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/CoordinateSystems.py#L241)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems.py#L241?message=Update%20Docs)]
</div>

 </div>
</div>






___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Molecools/CoordinateSystems/MolecularCartesianCoordinateSystem.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Molecools/CoordinateSystems/MolecularCartesianCoordinateSystem.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Molecools/CoordinateSystems/MolecularCartesianCoordinateSystem.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Molecools/CoordinateSystems/MolecularCartesianCoordinateSystem.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/CoordinateSystems.py#L186?message=Update%20Docs)