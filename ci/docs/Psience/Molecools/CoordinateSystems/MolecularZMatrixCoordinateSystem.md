## <a id="Psience.Molecools.CoordinateSystems.MolecularZMatrixCoordinateSystem">MolecularZMatrixCoordinateSystem</a>
Mirrors the standard ZMatrix coordinate system in _almost_ all regards, but forces an embedding

### Properties and Methods
```python
name: str
```
<a id="Psience.Molecools.CoordinateSystems.MolecularZMatrixCoordinateSystem.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, molecule, converter_options=None, **opts): 
```

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

<a id="Psience.Molecools.CoordinateSystems.MolecularZMatrixCoordinateSystem.axes" class="docs-object-method">&nbsp;</a>
```python
@property
axes(self): 
```

<a id="Psience.Molecools.CoordinateSystems.MolecularZMatrixCoordinateSystem.pre_convert" class="docs-object-method">&nbsp;</a>
```python
pre_convert(self, system): 
```

<a id="Psience.Molecools.CoordinateSystems.MolecularZMatrixCoordinateSystem.set_embedding" class="docs-object-method">&nbsp;</a>
```python
set_embedding(self): 
```

<a id="Psience.Molecools.CoordinateSystems.MolecularZMatrixCoordinateSystem.jacobian" class="docs-object-method">&nbsp;</a>
```python
jacobian(self, *args, reembed=None, strip_dummies=None, converter_options=None, **kwargs): 
```

### Examples


