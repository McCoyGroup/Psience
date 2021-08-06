## <a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunction">PerturbationTheoryWavefunction</a>
These things are fed the first and second order corrections

### Properties and Methods
<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunction.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, mol, basis, corrections): 
```

- `mol`: `Molecule`
    >the molecule the wavefunction is for
- `basis`: `RepresentationBasis`
    >the basis the expansion is being done in
- `corrections`: `PerturbationTheoryCorrections`
    >the corrections to the terms

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunction.order" class="docs-object-method">&nbsp;</a>
```python
@property
order(self): 
```

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunction.expectation" class="docs-object-method">&nbsp;</a>
```python
expectation(self, operator, other): 
```

<a id="Psience.VPT2.Wavefunctions.PerturbationTheoryWavefunction.zero_order_energy" class="docs-object-method">&nbsp;</a>
```python
@property
zero_order_energy(self): 
```

### Examples


___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/edit/ci/examples/ci/docs/Psience/VPT2/Wavefunctions/PerturbationTheoryWavefunction.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/examples/ci/docs/Psience/VPT2/Wavefunctions/PerturbationTheoryWavefunction.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/edit/ci/docs/ci/docs/Psience/VPT2/Wavefunctions/PerturbationTheoryWavefunction.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/docs/templates/ci/docs/Psience/VPT2/Wavefunctions/PerturbationTheoryWavefunction.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Wavefunctions.py?message=Update%20Docs)