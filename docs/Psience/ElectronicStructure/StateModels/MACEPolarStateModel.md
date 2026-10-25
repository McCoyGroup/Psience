## <a id="Psience.ElectronicStructure.StateModels.MACEPolarStateModel">MACEPolarStateModel</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/StateModels.py#L146)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/StateModels.py#L146?message=Update%20Docs)]
</div>

MACE-POLAR atomic multipoles (`density_coefficients`) and spin-channel multipoles
(`spin_charge_density`, spin density = up channel minus down channel, which integrates to
`multiplicity - 1`), Gaussian-smeared with the model's width.

Benchmark note (MACE-POLAR-1-S, 2026-10-20): the spin channels carry per-atom "spin" of order
+-1 even for closed-shell singlets (it integrates to zero) and miss the textbook spin patterns
of allyl, benzyl and phenoxyl radicals, so they should not be read as a physical spin density
without validation. Total-charge differences (Fukui functions) behave sensibly.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
supports_spin: bool
```
<a id="Psience.ElectronicStructure.StateModels.MACEPolarStateModel.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, backend=None, atoms=None, **backend_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/StateModels.py#L159)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/StateModels.py#L159?message=Update%20Docs)]
</div>

  - `backend`: `MACEPolarBackend | None`
    > a `MACEPolarBackend` (built from `atoms` and `backend_opts` if not given)


<a id="Psience.ElectronicStructure.StateModels.MACEPolarStateModel.from_molecule" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_molecule(cls, mol, tag=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L172)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L172?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.StateModels.MACEPolarStateModel.smearing" class="docs-object-method">&nbsp;</a> 
```python
@property
smearing(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/StateModels/MACEPolarStateModel.py#L178)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/StateModels/MACEPolarStateModel.py#L178?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.StateModels.MACEPolarStateModel.evaluate_state" class="docs-object-method">&nbsp;</a> 
```python
evaluate_state(self, coords, charge=0, multiplicity=1): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/StateModels/MACEPolarStateModel.py#L185)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/StateModels/MACEPolarStateModel.py#L185?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/master/ci/docs/examples/Psience/ElectronicStructure/StateModels/MACEPolarStateModel.md)/[New](https://github.com/McCoyGroup/Psience/new/master/?filename=ci/docs/examples/Psience/ElectronicStructure/StateModels/MACEPolarStateModel.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/ci/docs/templates/Psience/ElectronicStructure/StateModels/MACEPolarStateModel.md)/[New](https://github.com/McCoyGroup/Psience/new/master/?filename=ci/docs/templates/Psience/ElectronicStructure/StateModels/MACEPolarStateModel.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/StateModels.py#L146?message=Update%20Docs)   
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