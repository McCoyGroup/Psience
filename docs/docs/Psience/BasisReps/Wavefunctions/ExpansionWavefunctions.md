## <a id="Psience.BasisReps.Wavefunctions.ExpansionWavefunctions">ExpansionWavefunctions</a>
Simple expansion wave function rep that takes multiple sets of coefficients.
As with all objects deriving from `Wavefunctions`, can be iterated through to
provide a manifold of standalone `ExpansionWavefunction` objects.
Currently there are major conceptual issues, as I need this to _both_ support `AnalyticWavefunctions`
and `RepresentationBasis` as the basis...
which means `AnalyticWavefunctions` needs to track a basis...
but `AnalyticWavefunctions` wasn't really designed for that, so I need to go back and figure out how
that binding is going to be managed.

### Properties and Methods
<a id="Psience.BasisReps.Wavefunctions.ExpansionWavefunctions.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, energies, coefficients, basis_wfns, **ops): 
```

- `energies`: `Iterable[float]`
    >energies for the stored wavefunctions
- `coefficients`: `Iterable[Iterable[float]]`
    >expansion coefficients
- `basis_wfns`: `Wavefunctions`
    >wavefunctions to use as the basis for the expansion
- `ops`: `Any`
    >extra options for feeding through to `Wavefunctions`

### Examples


