## <a id="Psience.ElectronicStructure.MLIPBackends.MACEPolarBackend">MACEPolarBackend</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends.py#L1007)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends.py#L1007?message=Update%20Docs)]
</div>

MACE-POLAR: energies plus the Gaussian-smeared atomic multipoles it predicts
(`density_coefficients`, real `m = -l..l` moments in e Angstrom^l in the `graph_longrange`
GTO normalization) and their spin channels (`spin_charge_density`), with smearing width
`smearing` (Angstrom).







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
model_types: dict
default_model: str
property_units: dict
property_shapes: dict
polar_properties: tuple
```
<a id="Psience.ElectronicStructure.MLIPBackends.MACEPolarBackend.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, atoms, model=None, charge=None, multiplicity=None, external_field=None, calculator=None, device=None, model_dir=None, default_dtype='float64', quiet=True, smearing=None, info=None, **settings): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends.py#L1030)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends.py#L1030?message=Update%20Docs)]
</div>

  - `model`: `Any`
    > `'polar-1-s'`/`'polar-1-m'`/`'polar-1-l'` (or `'small'`/`'medium'`/`'large'`)
  - `external_field`: `Any`
    > uniform external field (`atoms.info["external_field"]`), in the model's native units
  - `calculator`: `Any`
    > an existing MACE-POLAR ASE calculator to reuse
  - `smearing`: `Any`
    > override for the model's smearing width (Angstrom)


<a id="Psience.ElectronicStructure.MLIPBackends.MACEPolarBackend.handle_specialization" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
handle_specialization(cls, tag): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1051)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1051?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.MLIPBackends.MACEPolarBackend.model_smearing" class="docs-object-method">&nbsp;</a> 
```python
@property
model_smearing(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends/MACEPolarBackend.py#L1055)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends/MACEPolarBackend.py#L1055?message=Update%20Docs)]
</div>
the smearing width the loaded model uses (Angstrom)


<a id="Psience.ElectronicStructure.MLIPBackends.MACEPolarBackend.smearing" class="docs-object-method">&nbsp;</a> 
```python
@property
smearing(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends/MACEPolarBackend.py#L1063)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends/MACEPolarBackend.py#L1063?message=Update%20Docs)]
</div>
the smearing width in use (Angstrom): the override, else the model's


<a id="Psience.ElectronicStructure.MLIPBackends.MACEPolarBackend.evaluate_polar_properties" class="docs-object-method">&nbsp;</a> 
```python
evaluate_polar_properties(self, coords, properties=None, charge=None, multiplicity=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/MLIPBackends/MACEPolarBackend.py#L1070)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends/MACEPolarBackend.py#L1070?message=Update%20Docs)]
</div>
Runs the MACE-POLAR calculator and returns its raw results.
  - `coords`: `Any`
    > nuclear coordinates, `(..., natoms, 3)`, in Angstrom
  - `properties`: `Any`
    > calculator result keys (e.g. `'electrostatic_energy'`, `'spin_charge_density'`, `'dipole'`)
  - `charge`: `Any`
    > total charge override for this call
  - `multiplicity`: `Any`
    > spin multiplicity override for this call
  - `:returns`: `dict`
    > arrays keyed by property, each with the leading batch shape of `coords`
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/ElectronicStructure/MLIPBackends/MACEPolarBackend.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/ElectronicStructure/MLIPBackends/MACEPolarBackend.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/ElectronicStructure/MLIPBackends/MACEPolarBackend.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/ElectronicStructure/MLIPBackends/MACEPolarBackend.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/MLIPBackends.py#L1007?message=Update%20Docs)   
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