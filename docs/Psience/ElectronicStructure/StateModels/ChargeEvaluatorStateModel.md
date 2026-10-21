## <a id="Psience.ElectronicStructure.StateModels.ChargeEvaluatorStateModel">ChargeEvaluatorStateModel</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/StateModels.py#L247)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/StateModels.py#L247?message=Update%20Docs)]
</div>

Partial charges from a factory of charge evaluators, one per `(charge, multiplicity)`: any
object with `evaluate(coords, order=0) -> [charges, ...]`, `distance_units`, and
`target_property_units` (e.g. a `Psience.Molecools` `ChargeEvaluator`). No spin resolution;
sites are point charges unless a `smearing` (Bohr) is supplied for visualization.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.ElectronicStructure.StateModels.ChargeEvaluatorStateModel.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, evaluator_factory, atoms, smearing=None, distance_units=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/StateModels.py#L255)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/StateModels.py#L255?message=Update%20Docs)]
</div>

  - `evaluator_factory`: `callable`
    > `(charge, multiplicity) -> ChargeEvaluator`


<a id="Psience.ElectronicStructure.StateModels.ChargeEvaluatorStateModel.from_molecule" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_molecule(cls, mol, tag=None, evaluator=None, smearing=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L267)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L267?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.StateModels.ChargeEvaluatorStateModel.get_evaluator" class="docs-object-method">&nbsp;</a> 
```python
get_evaluator(self, charge, multiplicity): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/StateModels/ChargeEvaluatorStateModel.py#L278)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/StateModels/ChargeEvaluatorStateModel.py#L278?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.StateModels.ChargeEvaluatorStateModel.evaluate_state" class="docs-object-method">&nbsp;</a> 
```python
evaluate_state(self, coords, charge=0, multiplicity=1): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/StateModels/ChargeEvaluatorStateModel.py#L284)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/StateModels/ChargeEvaluatorStateModel.py#L284?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/ElectronicStructure/StateModels/ChargeEvaluatorStateModel.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/ElectronicStructure/StateModels/ChargeEvaluatorStateModel.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/ElectronicStructure/StateModels/ChargeEvaluatorStateModel.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/ElectronicStructure/StateModels/ChargeEvaluatorStateModel.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/StateModels.py#L247?message=Update%20Docs)   
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