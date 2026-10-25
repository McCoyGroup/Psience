## <a id="Psience.ElectronicStructure.StateModels.ElectronicStateModel">ElectronicStateModel</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/StateModels.py#L73)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/StateModels.py#L73?message=Update%20Docs)]
</div>

Predicts `ElectronicState`s. Concrete models take coordinates in their `distance_units`;
`evaluate_state_bohr` takes Bohr.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
distance_units: str
supports_spin: bool
```
<a id="Psience.ElectronicStructure.StateModels.ElectronicStateModel.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, atoms): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/StateModels.py#L81)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/StateModels.py#L81?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.StateModels.ElectronicStateModel.evaluate_state" class="docs-object-method">&nbsp;</a> 
```python
evaluate_state(self, coords, charge=0, multiplicity=1): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/StateModels/ElectronicStateModel.py#L84)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/StateModels/ElectronicStateModel.py#L84?message=Update%20Docs)]
</div>

  - `coords`: `np.ndarray`
    > nuclear coordinates, `(natoms, 3)` in `distance_units`
  - `:returns`: `ElectronicState`
    >


<a id="Psience.ElectronicStructure.StateModels.ElectronicStateModel.evaluate_state_bohr" class="docs-object-method">&nbsp;</a> 
```python
evaluate_state_bohr(self, coords, charge=0, multiplicity=1): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/StateModels/ElectronicStateModel.py#L93)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/StateModels/ElectronicStateModel.py#L93?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.StateModels.ElectronicStateModel.evaluate_states" class="docs-object-method">&nbsp;</a> 
```python
evaluate_states(self, coords, states): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/StateModels/ElectronicStateModel.py#L99)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/StateModels/ElectronicStateModel.py#L99?message=Update%20Docs)]
</div>

  - `states`: `Iterable[tuple[int, int]]`
    > `(charge, multiplicity)` pairs
  - `:returns`: `list[ElectronicState]`
    >


<a id="Psience.ElectronicStructure.StateModels.ElectronicStateModel.get_models" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_models(cls): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L107)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L107?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.StateModels.ElectronicStateModel.from_molecule" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_molecule(cls, mol, model='mace-polar', **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L117)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L117?message=Update%20Docs)]
</div>
Resolve a model for a `Molecule`.
  - `model`: `Any`
    > a model instance, `'mace-polar[:size]'`, `'charges:backend[:tag]'` (any
    `ElectronicStructureBackend` that provides charges, e.g. `'charges:aimnet2'`,
    `'charges:xtb'`, `'charges:pyscf:b3lyp/def2-svp'`), or `'charge-evaluator:name'` (a
    charge evaluator from `mol.get_charge_evaluator`)
  - `:returns`: `ElectronicStateModel`
    >
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
[Edit](https://github.com/McCoyGroup/Psience/edit/master/ci/docs/examples/Psience/ElectronicStructure/StateModels/ElectronicStateModel.md)/[New](https://github.com/McCoyGroup/Psience/new/master/?filename=ci/docs/examples/Psience/ElectronicStructure/StateModels/ElectronicStateModel.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/ci/docs/templates/Psience/ElectronicStructure/StateModels/ElectronicStateModel.md)/[New](https://github.com/McCoyGroup/Psience/new/master/?filename=ci/docs/templates/Psience/ElectronicStructure/StateModels/ElectronicStateModel.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/StateModels.py#L73?message=Update%20Docs)   
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