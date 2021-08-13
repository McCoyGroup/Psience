## <a id="Psience.VPT2.Runner.VPTSystem">VPTSystem</a>
Provides a little helper for setting up the input
system for a VPT job

### Properties and Methods
<a id="Psience.VPT2.Runner.VPTSystem.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, mol, internals=None, dummy_atoms=None, modes=None, mode_selection=None, potential_derivatives=None, potential_function=None, dipole_derivatives=None): 
```

- `mol`: `str | Molecule`
    >the molecule or system specification to use (doesn't really even need to be a molecule)
- `internals`: `Any`
    >the Z-matrix for the internal coordinates (in the future will support a general function for this too)
- `modes`: `Any`
    >the normal modes to use if not already supplied by the Molecule
- `mode_selection`: `Any`
    >the subset of normal modes to do perturbation theory on
- `potential_derivatives`: `Iterable[np.ndarray]`
    >the derivatives of the potential to use for expansions
- `dipole_derivatives`: `Iterable[np.ndarray]`
    >the set of dipole derivatives to use for expansions

<a id="Psience.VPT2.Runner.VPTSystem.nmodes" class="docs-object-method">&nbsp;</a>
```python
@property
nmodes(self): 
```

<a id="Psience.VPT2.Runner.VPTSystem.get_potential_derivatives" class="docs-object-method">&nbsp;</a>
```python
get_potential_derivatives(self, potential_function, order=2, **fd_opts): 
```

<a id="Psience.VPT2.Runner.VPTSystem.from_harmonic_scan" class="docs-object-method">&nbsp;</a>
```python
from_harmonic_scan(scan_array): 
```

### Examples




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/edit/ci/examples/ci/docs/Psience/VPT2/Runner/VPTSystem.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/examples/ci/docs/Psience/VPT2/Runner/VPTSystem.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/edit/ci/docs/ci/docs/Psience/VPT2/Runner/VPTSystem.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/docs/templates/ci/docs/Psience/VPT2/Runner/VPTSystem.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Runner.py?message=Update%20Docs)