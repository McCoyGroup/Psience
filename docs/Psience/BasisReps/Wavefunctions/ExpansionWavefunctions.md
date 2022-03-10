## <a id="Psience.BasisReps.Wavefunctions.ExpansionWavefunctions">ExpansionWavefunctions</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Wavefunctions.py#L191)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Wavefunctions.py#L191?message=Update%20Docs)]
</div>

Simple expansion wave function rep that takes multiple sets of coefficients.
As with all objects deriving from `Wavefunctions`, can be iterated through to
provide a manifold of standalone `ExpansionWavefunction` objects.
Currently there are major conceptual issues, as I need this to _both_ support `AnalyticWavefunctions`
and `RepresentationBasis` as the basis...
which means `AnalyticWavefunctions` needs to track a basis...
but `AnalyticWavefunctions` wasn't really designed for that, so I need to go back and figure out how
that binding is going to be managed.

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

<a id="Psience.BasisReps.Wavefunctions.ExpansionWavefunctions.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, energies, coefficients, basis_wfns, **ops): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/BasisReps/Wavefunctions.py#L202)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Wavefunctions.py#L202?message=Update%20Docs)]
</div>


- `energies`: `Iterable[float]`
    >energies for the stored wavefunctions
- `coefficients`: `Iterable[Iterable[float]]`
    >expansion coefficients
- `basis_wfns`: `Wavefunctions`
    >wavefunctions to use as the basis for the expansion
- `ops`: `Any`
    >extra options for feeding through to `Wavefunctions`

 </div>
</div>




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/Wavefunctions/ExpansionWavefunctions.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/Wavefunctions/ExpansionWavefunctions.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/Wavefunctions/ExpansionWavefunctions.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/Wavefunctions/ExpansionWavefunctions.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/BasisReps/Wavefunctions.py#L191?message=Update%20Docs)