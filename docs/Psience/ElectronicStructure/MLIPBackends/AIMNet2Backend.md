## <a id="Psience.ElectronicStructure.MLIPBackends.AIMNet2Backend">AIMNet2Backend</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends.py#L38)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends.py#L38?message=Update%20Docs)]
</div>

AIMNet2 (and AIMNet2-NSE) energies, charges, and fluctuating-charge polarizabilities, with
derivatives by autograd. Borrows structure from `AIMNet2ASE` to call the model directly.







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
```
<a id="Psience.ElectronicStructure.MLIPBackends.AIMNet2Backend.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, atoms, model='aimnet2', charge=None, device=None, model_dir=None, multiplicity=None, quiet=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends.py#L60)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends.py#L60?message=Update%20Docs)]
</div>

  - `model`: `Any`
    > an AIMNet2 model name (`'aimnet2'`, `'aimnet2nse'`, ...) or a loaded `AIMNet2Calculator`
  - `device`: `Any`
    > the torch device (default: CUDA if available)
  - `model_dir`: `Any`
    > a cache directory for downloaded weights (default `$MODEL_CACHE_DIR`)


<a id="Psience.ElectronicStructure.MLIPBackends.AIMNet2Backend.handle_specialization" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
handle_specialization(cls, tag): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L70)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L70?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.MLIPBackends.AIMNet2Backend.setup_model" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
setup_model(cls, model, device=None, model_dir=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L154)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L154?message=Update%20Docs)]
</div>
Load an `AIMNet2Calculator` (or return `model` if it already is one).


<a id="Psience.ElectronicStructure.MLIPBackends.AIMNet2Backend.autodiff" class="docs-object-method">&nbsp;</a> 
```python
@staticmethod
autodiff(expr, coord, pad_dim, base_shape=None, create_graph=False, retain_graph=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/staticmethod.py#L190)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/staticmethod.py#L190?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.MLIPBackends.AIMNet2Backend.prep_eval" class="docs-object-method">&nbsp;</a> 
```python
prep_eval(self, coords, order, charge=None, multiplicity=None, keep_graph=None, forces=None, hessian=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.py#L250)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.py#L250?message=Update%20Docs)]
</div>
Run the model on `coords` (`(..., natoms, 3)`, Angstrom), setting up autograd for
`order` energy derivatives.
  - `:returns`: `_`
    > `(base_shape, flat_coords, data)`


<a id="Psience.ElectronicStructure.MLIPBackends.AIMNet2Backend.process_aimnet_derivs" class="docs-object-method">&nbsp;</a> 
```python
process_aimnet_derivs(self, base_shape, coords, data, root_key, order, extra_deriv_coord=None, deriv_key=None, property_shape=(), pad_coord=1, replace_nans=True, output_shape=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.py#L381)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.py#L381?message=Update%20Docs)]
</div>
Collect `data[root_key]` and its coordinate derivatives through `order` (and, with
`extra_deriv_coord=(var, n)`, mixed derivatives through order `n` in `var`) as numpy arrays.


<a id="Psience.ElectronicStructure.MLIPBackends.AIMNet2Backend.evaluate_energy" class="docs-object-method">&nbsp;</a> 
```python
evaluate_energy(self, coords, order, charge=None, multiplicity=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.py#L481)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.py#L481?message=Update%20Docs)]
</div>
`[E, dE/dx, ...]` on a flat batch `(nstruct, natoms, 3)` (Hessians are taken one structure at a time)


<a id="Psience.ElectronicStructure.MLIPBackends.AIMNet2Backend.evaluate_charges" class="docs-object-method">&nbsp;</a> 
```python
evaluate_charges(self, coords, order, charge=None, multiplicity=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.py#L493)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.py#L493?message=Update%20Docs)]
</div>
`[q, dq/dx, ...]` on a flat batch `(nstruct, natoms, 3)` (derivatives are taken one
structure at a time: the padded per-atom autograd doesn't separate batched molecules)


<a id="Psience.ElectronicStructure.MLIPBackends.AIMNet2Backend.evaluate_electrostatic_potential_derivatives" class="docs-object-method">&nbsp;</a> 
```python
evaluate_electrostatic_potential_derivatives(self, coords, coord_order, electrostatic_order, charge=None, multiplicity=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.py#L514)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.py#L514?message=Update%20Docs)]
</div>
The model's Coulomb energy and its derivatives with respect to the atomic charges (through
`electrostatic_order`) and coordinates (through `coord_order`), as
`[[E, dE/dx, ...], [dE/dq, ...], [d2E/dq2, ...]]`


<a id="Psience.ElectronicStructure.MLIPBackends.AIMNet2Backend.evaluate_field_response" class="docs-object-method">&nbsp;</a> 
```python
evaluate_field_response(self, coords, coord_order=0, field_order=2, charge=None, multiplicity=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.py#L535)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.py#L535?message=Update%20Docs)]
</div>
The fluctuating-charge dipole and polarizability expansions,
`[[mu, dmu/dx, ...], [alpha, dalpha/dx, ...]][:field_order]`.
  - `coords`: `Any`
    > `(..., natoms, 3)`, Angstrom


<a id="Psience.ElectronicStructure.MLIPBackends.AIMNet2Backend.evaluate_properties" class="docs-object-method">&nbsp;</a> 
```python
evaluate_properties(self, coords, orders, charge, multiplicity, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.py#L549)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.py#L549?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/master/ci/docs/examples/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.md)/[New](https://github.com/McCoyGroup/Psience/new/master/?filename=ci/docs/examples/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/ci/docs/templates/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.md)/[New](https://github.com/McCoyGroup/Psience/new/master/?filename=ci/docs/templates/Psience/ElectronicStructure/MLIPBackends/AIMNet2Backend.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends.py#L38?message=Update%20Docs)   
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