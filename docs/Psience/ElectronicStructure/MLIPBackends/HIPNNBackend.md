## <a id="Psience.ElectronicStructure.MLIPBackends.HIPNNBackend">HIPNNBackend</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends.py#L565)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends.py#L565?message=Update%20Docs)]
</div>

A `hippynn` model loaded from a training directory, with energies (`energy_key`, gradients
from `gradient_key`, higher derivatives by autograd) and dipoles (`dipole_key`).







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
distance_units: str
property_units: dict
analytic_derivative_orders: dict
property_shapes: dict
energy_key: str
gradient_key: str
dipole_key: str
dipole_gradient_key: str
```
<a id="Psience.ElectronicStructure.MLIPBackends.HIPNNBackend.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, atoms, model_dir, energy_key=<McUtils.Devutils.core.DefaultType instance>, gradient_key=<McUtils.Devutils.core.DefaultType instance>, dipole_key=<McUtils.Devutils.core.DefaultType instance>, charge=None, multiplicity=None, device=None, quiet=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends.py#L587)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends.py#L587?message=Update%20Docs)]
</div>

  - `model_dir`: `Any`
    > the directory `hippynn` saved the experiment to
  - `energy_key`: `Any`
    > the graph node holding the energy (default `'energies'`)
  - `gradient_key`: `Any`
    > the graph node holding its gradient (default `'Grad'`; `None` to use autograd)
  - `dipole_key`: `Any`
    > the graph node holding the dipole (default `'dipole'`)


<a id="Psience.ElectronicStructure.MLIPBackends.HIPNNBackend.handle_specialization" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
handle_specialization(cls, tag): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L608)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L608?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.MLIPBackends.HIPNNBackend.setup_model" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
setup_model(cls, model_dir, *, device=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L612)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L612?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.MLIPBackends.HIPNNBackend.print_model_properties" class="docs-object-method">&nbsp;</a> 
```python
print_model_properties(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends/HIPNNBackend.py#L630)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends/HIPNNBackend.py#L630?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.MLIPBackends.HIPNNBackend.enable_charge_equilibration_outputs" class="docs-object-method">&nbsp;</a> 
```python
enable_charge_equilibration_outputs(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends/HIPNNBackend.py#L633)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends/HIPNNBackend.py#L633?message=Update%20Docs)]
</div>
Rebuild the predictor so it also exposes the charge-equilibration (`ChEQ`) Coulomb energy
and its inputs, for fluctuating-charge responses.


<a id="Psience.ElectronicStructure.MLIPBackends.HIPNNBackend.autodiff" class="docs-object-method">&nbsp;</a> 
```python
@staticmethod
autodiff(expr, coord_list, pad_dim, create_graph=False, retain_graph=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/staticmethod.py#L648)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/staticmethod.py#L648?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.MLIPBackends.HIPNNBackend.process_derivatives" class="docs-object-method">&nbsp;</a> 
```python
process_derivatives(self, expr, coords, order): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends/HIPNNBackend.py#L688)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends/HIPNNBackend.py#L688?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.MLIPBackends.HIPNNBackend.process_results" class="docs-object-method">&nbsp;</a> 
```python
process_results(self, base_shape, coords, outputs, order, property_key=None, gradient_key=None, property_shape=(), output_shape=None, **etc): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends/HIPNNBackend.py#L705)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends/HIPNNBackend.py#L705?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.MLIPBackends.HIPNNBackend.prep_eval" class="docs-object-method">&nbsp;</a> 
```python
prep_eval(self, coords, order, property_key=<McUtils.Devutils.core.DefaultType instance>, gradient_key=<McUtils.Devutils.core.DefaultType instance>, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends/HIPNNBackend.py#L739)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends/HIPNNBackend.py#L739?message=Update%20Docs)]
</div>
Run the predictor on `coords` (`(..., natoms, 3)`, Angstrom), adding `property_key` and
`gradient_key` to its outputs if needed.
  - `:returns`: `_`
    > `(base_shape, torch_coords, outputs)`


<a id="Psience.ElectronicStructure.MLIPBackends.HIPNNBackend.evaluate_energy" class="docs-object-method">&nbsp;</a> 
```python
evaluate_energy(self, coords, order, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends/HIPNNBackend.py#L800)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends/HIPNNBackend.py#L800?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.MLIPBackends.HIPNNBackend.evaluate_dipole" class="docs-object-method">&nbsp;</a> 
```python
evaluate_dipole(self, coords, order, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends/HIPNNBackend.py#L807)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends/HIPNNBackend.py#L807?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.MLIPBackends.HIPNNBackend.evaluate_field_response" class="docs-object-method">&nbsp;</a> 
```python
evaluate_field_response(self, coords, coord_order=0, field_order=2, charge=None, multiplicity=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends/HIPNNBackend.py#L814)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends/HIPNNBackend.py#L814?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.MLIPBackends.HIPNNBackend.evaluate_properties" class="docs-object-method">&nbsp;</a> 
```python
evaluate_properties(self, coords, orders, charge, multiplicity, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends/HIPNNBackend.py#L821)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends/HIPNNBackend.py#L821?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/master/ci/docs/examples/Psience/ElectronicStructure/MLIPBackends/HIPNNBackend.md)/[New](https://github.com/McCoyGroup/Psience/new/master/?filename=ci/docs/examples/Psience/ElectronicStructure/MLIPBackends/HIPNNBackend.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/ci/docs/templates/Psience/ElectronicStructure/MLIPBackends/HIPNNBackend.md)/[New](https://github.com/McCoyGroup/Psience/new/master/?filename=ci/docs/templates/Psience/ElectronicStructure/MLIPBackends/HIPNNBackend.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends.py#L565?message=Update%20Docs)   
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