## <a id="Psience.Molecools.Vibrations.MolecularVibrations">MolecularVibrations</a>


### Properties and Methods
<a id="Psience.Molecools.Vibrations.MolecularVibrations.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, molecule, basis, freqs=None, init=None): 
```
Sets up a vibration for a Molecule object over the CoordinateSystem basis
- `molecule`: `Molecule`
    >No description...
- `init`: `None | CoordinateSet`
    >No description...
- `basis`: `MolecularNormalModes`
    >No description...

<a id="Psience.Molecools.Vibrations.MolecularVibrations.basis" class="docs-object-method">&nbsp;</a>
```python
@property
basis(self): 
```

<a id="Psience.Molecools.Vibrations.MolecularVibrations.molecule" class="docs-object-method">&nbsp;</a>
```python
@property
molecule(self): 
```

<a id="Psience.Molecools.Vibrations.MolecularVibrations.coords" class="docs-object-method">&nbsp;</a>
```python
@property
coords(self): 
```

<a id="Psience.Molecools.Vibrations.MolecularVibrations.__len__" class="docs-object-method">&nbsp;</a>
```python
__len__(self): 
```

<a id="Psience.Molecools.Vibrations.MolecularVibrations.displace" class="docs-object-method">&nbsp;</a>
```python
displace(self, displacements=None, amt=0.1, n=1, which=0): 
```

<a id="Psience.Molecools.Vibrations.MolecularVibrations.visualize" class="docs-object-method">&nbsp;</a>
```python
visualize(self, step_size=0.1, steps=(5, 5), which=0, anim_opts=None, mode='fast', **plot_args): 
```

<a id="Psience.Molecools.Vibrations.MolecularVibrations.__getitem__" class="docs-object-method">&nbsp;</a>
```python
__getitem__(self, item): 
```
Takes a slice of the modes
- `item`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.Molecools.Vibrations.MolecularVibrations.embed" class="docs-object-method">&nbsp;</a>
```python
embed(self, frame): 
```

- `frame`: `MolecularTransformation`
    >No description...
- `:returns`: `_`
    >No description...

### Examples

