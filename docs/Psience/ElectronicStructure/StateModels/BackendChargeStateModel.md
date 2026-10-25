## <a id="Psience.ElectronicStructure.StateModels.BackendChargeStateModel">BackendChargeStateModel</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/StateModels.py#L207)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/StateModels.py#L207?message=Update%20Docs)]
</div>

Partial charges (and energies) from any `ElectronicStructureBackend` that provides
`'charges'` (AIMNet2, xTB, PySCF, ASE calculators that report charges), evaluated per state
with the backend's per-call charge/multiplicity overrides. No spin resolution; sites are
point charges unless a `smearing` (Bohr) is supplied for visualization.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.ElectronicStructure.StateModels.BackendChargeStateModel.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, backend, smearing=None, with_energy=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/StateModels.py#L215)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/StateModels.py#L215?message=Update%20Docs)]
</div>

  - `backend`: `ElectronicStructureBackend`
    > 
  - `with_energy`: `Any`
    > also evaluate the state energy


<a id="Psience.ElectronicStructure.StateModels.BackendChargeStateModel.from_molecule" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_molecule(cls, mol, tag=None, backend=None, smearing=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L226)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L226?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.StateModels.BackendChargeStateModel.evaluate_state" class="docs-object-method">&nbsp;</a> 
```python
evaluate_state(self, coords, charge=0, multiplicity=1): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/StateModels/BackendChargeStateModel.py#L233)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/StateModels/BackendChargeStateModel.py#L233?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/master/ci/docs/examples/Psience/ElectronicStructure/StateModels/BackendChargeStateModel.md)/[New](https://github.com/McCoyGroup/Psience/new/master/?filename=ci/docs/examples/Psience/ElectronicStructure/StateModels/BackendChargeStateModel.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/ci/docs/templates/Psience/ElectronicStructure/StateModels/BackendChargeStateModel.md)/[New](https://github.com/McCoyGroup/Psience/new/master/?filename=ci/docs/templates/Psience/ElectronicStructure/StateModels/BackendChargeStateModel.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/StateModels.py#L207?message=Update%20Docs)   
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