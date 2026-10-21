## <a id="Psience.ElectronicStructure.QuantumChemistryBackends.XTBBackend">XTBBackend</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/QuantumChemistryBackends.py#L20)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/QuantumChemistryBackends.py#L20?message=Update%20Docs)]
</div>

GFN-xTB energies, gradients, Mulliken charges, and dipoles through the `xtb` python bindings
(conda-forge: `xtb-python`). Every requested property comes from one single point.
  - `atoms`: `Any`
    > the atom symbols
  - `method`: `Any`
    > ``'GFN2-xTB'`` (default), ``'GFN1-xTB'`` or ``'GFN0-xTB'``
  - `charge`: `Any`
    > the total charge (default 0)
  - `multiplicity`: `Any`
    > the spin multiplicity (default 1); passed to xtb as ``multiplicity - 1`` unpaired electrons
  - `quiet`: `Any`
    > silence xtb's output
  - `accuracy`: `Any`
    > xtb's numerical accuracy parameter (smaller is tighter), if not the default
  - `max_iterations`: `Any`
    > the SCF iteration limit, if not the default
  - `electronic_temperature`: `Any`
    > the electronic temperature in Kelvin, if not the default
  - `solvent`: `Any`
    > an ALPB solvent name (e.g. ``'water'``)
  - `retry_iterations`: `Any`
    > if the SCF fails to converge, try once more with this iteration limit
    (the electronic temperature and accuracy are never changed behind your back); `None` disables the retry






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
```
<a id="Psience.ElectronicStructure.QuantumChemistryBackends.XTBBackend.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, atoms, method='GFN2-xTB', charge=None, multiplicity=None, quiet=True, accuracy=None, max_iterations=None, electronic_temperature=None, solvent=None, retry_iterations=1000): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/QuantumChemistryBackends.py#L55)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/QuantumChemistryBackends.py#L55?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.QuantumChemistryBackends.XTBBackend.handle_specialization" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
handle_specialization(cls, tag): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L75)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L75?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.QuantumChemistryBackends.XTBBackend.discover_dependencies" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
discover_dependencies(cls, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L79)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L79?message=Update%20Docs)]
</div>
The `xtb` python bindings come from conda-forge (``xtb-python``); the PyPI ``xtb`` package needs a
separately installed libxtb.


<a id="Psience.ElectronicStructure.QuantumChemistryBackends.XTBBackend.get_calculator" class="docs-object-method">&nbsp;</a> 
```python
get_calculator(self, coords, charge=None, multiplicity=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/QuantumChemistryBackends/XTBBackend.py#L94)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/QuantumChemistryBackends/XTBBackend.py#L94?message=Update%20Docs)]
</div>
Build an xtb `Calculator` at the given Cartesians (Bohr) with this backend's settings.


<a id="Psience.ElectronicStructure.QuantumChemistryBackends.XTBBackend.get_single_point" class="docs-object-method">&nbsp;</a> 
```python
get_single_point(self, coords, charge=None, multiplicity=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/QuantumChemistryBackends/XTBBackend.py#L124)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/QuantumChemistryBackends/XTBBackend.py#L124?message=Update%20Docs)]
</div>
Run one xTB single point, retrying once with a larger iteration limit if the SCF doesn't converge.


<a id="Psience.ElectronicStructure.QuantumChemistryBackends.XTBBackend.evaluate_properties" class="docs-object-method">&nbsp;</a> 
```python
evaluate_properties(self, coords, orders, charge, multiplicity, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/QuantumChemistryBackends/XTBBackend.py#L145)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/QuantumChemistryBackends/XTBBackend.py#L145?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/ElectronicStructure/QuantumChemistryBackends/XTBBackend.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/ElectronicStructure/QuantumChemistryBackends/XTBBackend.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/ElectronicStructure/QuantumChemistryBackends/XTBBackend.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/ElectronicStructure/QuantumChemistryBackends/XTBBackend.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/QuantumChemistryBackends.py#L20?message=Update%20Docs)   
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