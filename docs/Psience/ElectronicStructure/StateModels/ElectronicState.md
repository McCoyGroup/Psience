## <a id="Psience.ElectronicStructure.StateModels.ElectronicState">ElectronicState</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/StateModels.py#L24)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/StateModels.py#L24?message=Update%20Docs)]
</div>

The densities a model predicts for one charge/spin state at one geometry.

`charge_density` is the net charge density (nuclei included, so its populations are atomic
partial charges); `spin_density` is majority-minus-minority spin (integrating to
`multiplicity - 1`), or `None` if the model has no spin resolution.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.ElectronicStructure.StateModels.ElectronicState.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, atoms, coords, charge, multiplicity, charge_density, spin_density=None, energy=None, model=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/StateModels.py#L33)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/StateModels.py#L33?message=Update%20Docs)]
</div>

  - `coords`: `np.ndarray`
    > nuclear coordinates in Bohr
  - `charge_density`: `AtomCenteredDensity`
    > 
  - `spin_density`: `AtomCenteredDensity | None`
    > 
  - `energy`: `float | np.ndarray | None`
    > total energy in Hartree, if the model provides it


<a id="Psience.ElectronicStructure.StateModels.ElectronicState.partial_charges" class="docs-object-method">&nbsp;</a> 
```python
@property
partial_charges(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/StateModels/ElectronicState.py#L51)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/StateModels/ElectronicState.py#L51?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.StateModels.ElectronicState.electron_density" class="docs-object-method">&nbsp;</a> 
```python
@property
electron_density(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/StateModels/ElectronicState.py#L55)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/StateModels/ElectronicState.py#L55?message=Update%20Docs)]
</div>
The electron-number density implied by the charge density: populations `Z_i - q_i`,
higher moments `-p_ilm` (the nuclei carry no higher moments)


<a id="Psience.ElectronicStructure.StateModels.ElectronicState.__repr__" class="docs-object-method">&nbsp;</a> 
```python
__repr__(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/StateModels/ElectronicState.py#L66)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/StateModels/ElectronicState.py#L66?message=Update%20Docs)]
</div>
 </div>
</div>












---


<div markdown="1" class="text-secondary">
<div class="container">
  <div class="row">
   <div class="col" markdown="1">
**Feedback**   
</div>
   <div class="col" markdown="1">
**Examples**   
</div>
   <div class="col" markdown="1">
**Templates**   
</div>
   <div class="col" markdown="1">
**Documentation**   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
</div>
  <div class="row">
   <div class="col" markdown="1">
[Bug](https://github.com/McCoyGroup/Psience/issues/new?title=Documentation%20Improvement%20Needed)/[Request](https://github.com/McCoyGroup/Psience/issues/new?title=Example%20Request)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/ElectronicStructure/StateModels/ElectronicState.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/ElectronicStructure/StateModels/ElectronicState.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/ElectronicStructure/StateModels/ElectronicState.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/ElectronicStructure/StateModels/ElectronicState.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/StateModels.py#L24?message=Update%20Docs)   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
</div>
</div>
</div>