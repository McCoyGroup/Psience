## <a id="Psience.VPT2.Terms.PotentialTerms">PotentialTerms</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L1231)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L1231?message=Update%20Docs)]
</div>

A helper class that can transform the derivatives of the potential from Cartesian to normal coordinates

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

<a id="Psience.VPT2.Terms.PotentialTerms.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, molecule, mixed_derivs=None, modes=None, potential_derivatives=None, mode_selection=None, logger=None, parallelizer=None, checkpointer=None, check_input_force_constants=True, allow_higher_potential_terms=False, hessian_tolerance=0.0001, grad_tolerance=0.0001, freq_tolerance=0.002, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L1242)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L1242?message=Update%20Docs)]
</div>


- `mode_selection`: `None | Iterable[int]`
    >the subset of normal modes to use
- `modes`: `None | MolecularVibrations`
    >the normal modes to use when doing calculations
- `mixed_derivs`: `bool`
    >whether or not the pulled derivatives are partially derivatives along the normal coords
- `molecule`: `Molecule`
    >the molecule that will supply the potential derivatives

<a id="Psience.VPT2.Terms.PotentialTerms.v_derivs" class="docs-object-method">&nbsp;</a> 
```python
@property
v_derivs(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Terms.PotentialTerms.get_terms" class="docs-object-method">&nbsp;</a> 
```python
get_terms(self, order=None, logger=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L1548)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L1548?message=Update%20Docs)]
</div>

 </div>
</div>






___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Terms/PotentialTerms.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Terms/PotentialTerms.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Terms/PotentialTerms.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Terms/PotentialTerms.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L1231?message=Update%20Docs)