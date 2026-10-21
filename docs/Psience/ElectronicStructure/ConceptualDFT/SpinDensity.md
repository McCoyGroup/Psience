## <a id="Psience.ElectronicStructure.ConceptualDFT.SpinDensity">SpinDensity</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/ConceptualDFT.py#L26)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/ConceptualDFT.py#L26?message=Update%20Docs)]
</div>

The spin density (majority minus minority) of one state, optionally relative to a reference
state at the same geometry.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.ElectronicStructure.ConceptualDFT.SpinDensity.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, state, reference=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/ConceptualDFT.py#L32)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/ConceptualDFT.py#L32?message=Update%20Docs)]
</div>

  - `state`: `ElectronicState`
    > 
  - `reference`: `ElectronicState | None`
    > a state whose spin density is subtracted (e.g. the closed-shell neutral),
    a heuristic for models whose spin channels carry a closed-shell partition; the
    difference still integrates to `multiplicity - 1` if the reference is a singlet


<a id="Psience.ElectronicStructure.ConceptualDFT.SpinDensity.compute" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
compute(cls, model, coords, charge=0, multiplicity=2, reference=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L45)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L45?message=Update%20Docs)]
</div>

  - `model`: `ElectronicStateModel`
    > the state model
  - `coords`: `Any`
    > coordinates in the model's distance units
  - `reference`: `tuple[int, int] | None`
    > `(charge, multiplicity)` of a reference state to subtract, or `None`
  - `:returns`: `SpinDensity`
    >


<a id="Psience.ElectronicStructure.ConceptualDFT.SpinDensity.from_molecule" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_molecule(cls, mol, model='mace-polar', charge=None, multiplicity=None, reference=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L59)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L59?message=Update%20Docs)]
</div>

  - `mol`: `Molecule`
    > the molecule (its charge/spin are the defaults)


<a id="Psience.ElectronicStructure.ConceptualDFT.SpinDensity.density" class="docs-object-method">&nbsp;</a> 
```python
@property
density(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/ConceptualDFT/SpinDensity.py#L73)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/ConceptualDFT/SpinDensity.py#L73?message=Update%20Docs)]
</div>
the spin density as an `AtomCenteredDensity`


<a id="Psience.ElectronicStructure.ConceptualDFT.SpinDensity.populations" class="docs-object-method">&nbsp;</a> 
```python
@property
populations(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/ConceptualDFT/SpinDensity.py#L81)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/ConceptualDFT/SpinDensity.py#L81?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.ConceptualDFT.SpinDensity.participation_ratio" class="docs-object-method">&nbsp;</a> 
```python
participation_ratio(self, heavy_atoms_only=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/ConceptualDFT/SpinDensity.py#L85)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/ConceptualDFT/SpinDensity.py#L85?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.ConceptualDFT.SpinDensity.summary" class="docs-object-method">&nbsp;</a> 
```python
summary(self, top=8): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/ConceptualDFT/SpinDensity.py#L88)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/ConceptualDFT/SpinDensity.py#L88?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/ElectronicStructure/ConceptualDFT/SpinDensity.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/ElectronicStructure/ConceptualDFT/SpinDensity.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/ElectronicStructure/ConceptualDFT/SpinDensity.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/ElectronicStructure/ConceptualDFT/SpinDensity.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/ConceptualDFT.py#L26?message=Update%20Docs)   
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