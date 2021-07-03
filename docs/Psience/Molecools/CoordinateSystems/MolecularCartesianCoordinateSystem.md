## <a id="Psience.Molecools.CoordinateSystems.MolecularCartesianCoordinateSystem">MolecularCartesianCoordinateSystem</a>
Mirrors the standard Cartesian coordinate system in _almost_ all regards, but forces an embedding

### Properties and Methods
```python
name: str
```
<a id="Psience.Molecools.CoordinateSystems.MolecularCartesianCoordinateSystem.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, molecule, converter_options=None, **opts): 
```

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

<a id="Psience.Molecools.CoordinateSystems.MolecularCartesianCoordinateSystem.set_embedding" class="docs-object-method">&nbsp;</a>
```python
set_embedding(self): 
```
Sets up the embedding options...
- `:returns`: `_`
    >No description...

<a id="Psience.Molecools.CoordinateSystems.MolecularCartesianCoordinateSystem.jacobian" class="docs-object-method">&nbsp;</a>
```python
jacobian(self, coords, system, strip_dummies=None, converter_options=None, analytic_deriv_order=None, **kwargs): 
```

### Examples


