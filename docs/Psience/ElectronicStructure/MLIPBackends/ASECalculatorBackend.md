## <a id="Psience.ElectronicStructure.MLIPBackends.ASECalculatorBackend">ASECalculatorBackend</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends.py#L835)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends.py#L835?message=Update%20Docs)]
</div>

Any ASE calculator. Energies come with forces and (if the calculator has `get_hessian`)
Hessians; any other result the calculator produces (`'charges'`, `'dipole'`, model-specific
keys) is available without derivatives, and is computed in the same calculator pass as the
energy and forces.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
distance_units: str
default_multiplicity: int
property_units: dict
analytic_derivative_orders: dict
property_shapes: dict
hessian_function_attribute: str
```
<a id="Psience.ElectronicStructure.MLIPBackends.ASECalculatorBackend.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, atoms, charge=None, multiplicity=None, quiet=True, calculator=None, info=None, calc=None, **calculator_settings): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends.py#L859)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends.py#L859?message=Update%20Docs)]
</div>

  - `calculator`: `Any`
    > an ASE calculator (built by `setup_calculator(**calculator_settings)` if not given)
  - `info`: `dict | None`
    > extra entries for the ASE `Atoms.info` dict (charge and spin are set from the state)
  - `calc`: `Any`
    > an alias for `calculator`


<a id="Psience.ElectronicStructure.MLIPBackends.ASECalculatorBackend.setup_calculator" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
setup_calculator(cls, **settings): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L877)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L877?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.MLIPBackends.ASECalculatorBackend.calculator" class="docs-object-method">&nbsp;</a> 
```python
@property
calculator(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends/ASECalculatorBackend.py#L881)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends/ASECalculatorBackend.py#L881?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.MLIPBackends.ASECalculatorBackend.get_analytic_derivative_order" class="docs-object-method">&nbsp;</a> 
```python
get_analytic_derivative_order(self, prop): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends/ASECalculatorBackend.py#L885)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends/ASECalculatorBackend.py#L885?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.MLIPBackends.ASECalculatorBackend.prep_ase" class="docs-object-method">&nbsp;</a> 
```python
prep_ase(self, coords, charge=None, multiplicity=None, info=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends/ASECalculatorBackend.py#L889)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends/ASECalculatorBackend.py#L889?message=Update%20Docs)]
</div>
An `ASEMolecule` at `coords` with this backend's calculator attached.


<a id="Psience.ElectronicStructure.MLIPBackends.ASECalculatorBackend.calculate_results" class="docs-object-method">&nbsp;</a> 
```python
calculate_results(self, coords, keys, charge=None, multiplicity=None, info=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends/ASECalculatorBackend.py#L904)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends/ASECalculatorBackend.py#L904?message=Update%20Docs)]
</div>
The calculator's raw `results` for `keys` at each structure of `coords`.
  - `coords`: `Any`
    > `(..., natoms, 3)`, Angstrom
  - `:returns`: `dict`
    > arrays keyed by result name, with the leading batch shape of `coords`


<a id="Psience.ElectronicStructure.MLIPBackends.ASECalculatorBackend.evaluate_properties" class="docs-object-method">&nbsp;</a> 
```python
evaluate_properties(self, coords, orders, charge, multiplicity, info=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends/ASECalculatorBackend.py#L919)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends/ASECalculatorBackend.py#L919?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/master/ci/docs/examples/Psience/ElectronicStructure/MLIPBackends/ASECalculatorBackend.md)/[New](https://github.com/McCoyGroup/Psience/new/master/?filename=ci/docs/examples/Psience/ElectronicStructure/MLIPBackends/ASECalculatorBackend.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/ci/docs/templates/Psience/ElectronicStructure/MLIPBackends/ASECalculatorBackend.md)/[New](https://github.com/McCoyGroup/Psience/new/master/?filename=ci/docs/templates/Psience/ElectronicStructure/MLIPBackends/ASECalculatorBackend.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends.py#L835?message=Update%20Docs)   
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