## <a id="Psience.Molecools.CoordinateSystems.MolecularCartesianCoordinateSystem">MolecularCartesianCoordinateSystem</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Molecools/CoordinateSystems.py#L190)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Molecools/CoordinateSystems.py#L190?message=Update%20Docs)]
</div>

Mirrors the standard Cartesian coordinate system in _almost_ all regards, but forces an embedding

```python
name: str
```
<a id="Psience.Molecools.CoordinateSystems.MolecularCartesianCoordinateSystem.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, molecule, converter_options=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Molecools/CoordinateSystems.py#L195)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Molecools/CoordinateSystems.py#L195?message=Update%20Docs)]
</div>


- `molecule`: `AbstractMolecule`
    >No description...
- `converter_options`: `Any`
    >No description...
- `opts`: `Any`
    >No description...

<a id="Psience.Molecools.CoordinateSystems.MolecularCartesianCoordinateSystem.pre_convert" class="docs-object-method">&nbsp;</a> 
```python
pre_convert(self, system): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Molecools/CoordinateSystems.py#L212)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Molecools/CoordinateSystems.py#L212?message=Update%20Docs)]
</div>

<a id="Psience.Molecools.CoordinateSystems.MolecularCartesianCoordinateSystem.set_embedding" class="docs-object-method">&nbsp;</a> 
```python
set_embedding(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Molecools/CoordinateSystems.py#L215)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Molecools/CoordinateSystems.py#L215?message=Update%20Docs)]
</div>

Sets up the embedding options...
- `:returns`: `_`
    >No description...

<a id="Psience.Molecools.CoordinateSystems.MolecularCartesianCoordinateSystem.jacobian" class="docs-object-method">&nbsp;</a> 
```python
jacobian(self, coords, system, strip_dummies=None, converter_options=None, analytic_deriv_order=None, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/Molecools/CoordinateSystems.py#L244)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Molecools/CoordinateSystems.py#L244?message=Update%20Docs)]
</div>



___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/ci/docs/Psience/Molecools/CoordinateSystems/MolecularCartesianCoordinateSystem.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/ci/docs/Psience/Molecools/CoordinateSystems/MolecularCartesianCoordinateSystem.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/ci/docs/Psience/Molecools/CoordinateSystems/MolecularCartesianCoordinateSystem.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/ci/docs/Psience/Molecools/CoordinateSystems/MolecularCartesianCoordinateSystem.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Molecools/CoordinateSystems.py#L190?message=Update%20Docs)