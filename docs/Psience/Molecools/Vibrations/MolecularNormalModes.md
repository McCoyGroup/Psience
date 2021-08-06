## <a id="Psience.Molecools.Vibrations.MolecularNormalModes">MolecularNormalModes</a>
A Coordinerds CoordinateSystem object that manages all of the data needed to
work with normal mode coordinates + some convenience functions for generating and whatnot

### Properties and Methods
```python
name: str
from_force_constants: method
```
<a id="Psience.Molecools.Vibrations.MolecularNormalModes.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, molecule, coeffs, name=None, freqs=None, internal=False, origin=None, basis=None, inverse=None): 
```

- `molecule`: `AbstractMolecule`
    >No description...
- `coeffs`: `Any`
    >No description...
- `name`: `Any`
    >No description...
- `freqs`: `Any`
    >No description...
- `internal`: `Any`
    >No description...
- `origin`: `Any`
    >No description...
- `basis`: `Any`
    >No description...
- `inverse`: `Any`
    >No description...

<a id="Psience.Molecools.Vibrations.MolecularNormalModes.molecule" class="docs-object-method">&nbsp;</a>
```python
@property
molecule(self): 
```

<a id="Psience.Molecools.Vibrations.MolecularNormalModes.to_internals" class="docs-object-method">&nbsp;</a>
```python
to_internals(self, intcrds=None, dYdR=None, dRdY=None): 
```

<a id="Psience.Molecools.Vibrations.MolecularNormalModes.origin" class="docs-object-method">&nbsp;</a>
```python
@property
origin(self): 
```

<a id="Psience.Molecools.Vibrations.MolecularNormalModes.embed" class="docs-object-method">&nbsp;</a>
```python
embed(self, frame): 
```

- `frame`: `MolecularTransformation`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.Molecools.Vibrations.MolecularNormalModes.insert" class="docs-object-method">&nbsp;</a>
```python
insert(self, val, where): 
```
Inserts values into the appropriate positions in the mode matrix
- `val`: `Any`
    >No description...
- `where`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.Molecools.Vibrations.MolecularNormalModes.rescale" class="docs-object-method">&nbsp;</a>
```python
rescale(self, scaling_factors): 
```
Rescales each mode in the expansion matrix
        by the passed `scaling_factors`
- `scaling_factors`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.Molecools.Vibrations.MolecularNormalModes.__getitem__" class="docs-object-method">&nbsp;</a>
```python
__getitem__(self, item): 
```
Takes a slice of the modes
- `item`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

### Examples




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/edit/ci/examples/ci/docs/Psience/Molecools/Vibrations/MolecularNormalModes.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/examples/ci/docs/Psience/Molecools/Vibrations/MolecularNormalModes.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/edit/ci/docs/ci/docs/Psience/Molecools/Vibrations/MolecularNormalModes.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/docs/templates/ci/docs/Psience/Molecools/Vibrations/MolecularNormalModes.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/Molecools/Vibrations.py?message=Update%20Docs)