## <a id="Psience.ElectronicStructure.ConceptualDFT.FukuiFunctions">FukuiFunctions</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/ConceptualDFT.py#L102)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/ConceptualDFT.py#L102?message=Update%20Docs)]
</div>

Finite-difference (vertical) Fukui functions in the electron-density convention,

    f+ = n(N+1) - n(N)        (electron addition, LUMO-like)
    f- = n(N) - n(N-1)        (electron removal, HOMO-like)
    f0 = (f+ + f-) / 2,  dual descriptor = f+ - f-

from the predicted charge densities of the reference state and its anion/cation at the same
geometry. Each integrates to 1 (the dual descriptor to 0) and is an `AtomCenteredDensity`
(condensed Fukui indices are its populations; the continuous form uses the model's smearing).







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.ElectronicStructure.ConceptualDFT.FukuiFunctions.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, reference, cation=None, anion=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/ConceptualDFT.py#L115)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/ConceptualDFT.py#L115?message=Update%20Docs)]
</div>

  - `reference`: `ElectronicState`
    > 
  - `cation`: `ElectronicState | None`
    > 
  - `anion`: `ElectronicState | None`
    >


<a id="Psience.ElectronicStructure.ConceptualDFT.FukuiFunctions.compute" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
compute(cls, model, coords, charge=0, multiplicity=1, cation_multiplicity=None, anion_multiplicity=None, which=('+', '-')): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L125)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L125?message=Update%20Docs)]
</div>

  - `model`: `ElectronicStateModel`
    > the state model
  - `coords`: `Any`
    > coordinates in the model's distance units
  - `which`: `Any`
    > which ions to evaluate (`'+'` needs the anion, `'-'` the cation)
  - `:returns`: `FukuiFunctions`
    >


<a id="Psience.ElectronicStructure.ConceptualDFT.FukuiFunctions.from_molecule" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_molecule(cls, mol, model='mace-polar', charge=None, multiplicity=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L145)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L145?message=Update%20Docs)]
</div>

  - `mol`: `Molecule`
    > the molecule (its charge/spin are the defaults)


<a id="Psience.ElectronicStructure.ConceptualDFT.FukuiFunctions.f_plus" class="docs-object-method">&nbsp;</a> 
```python
@property
f_plus(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/ConceptualDFT/FukuiFunctions.py#L161)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/ConceptualDFT/FukuiFunctions.py#L161?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.ConceptualDFT.FukuiFunctions.f_minus" class="docs-object-method">&nbsp;</a> 
```python
@property
f_minus(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/ConceptualDFT/FukuiFunctions.py#L167)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/ConceptualDFT/FukuiFunctions.py#L167?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.ConceptualDFT.FukuiFunctions.f_zero" class="docs-object-method">&nbsp;</a> 
```python
@property
f_zero(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/ConceptualDFT/FukuiFunctions.py#L173)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/ConceptualDFT/FukuiFunctions.py#L173?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.ConceptualDFT.FukuiFunctions.dual_descriptor" class="docs-object-method">&nbsp;</a> 
```python
@property
dual_descriptor(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/ConceptualDFT/FukuiFunctions.py#L177)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/ConceptualDFT/FukuiFunctions.py#L177?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.ConceptualDFT.FukuiFunctions.summary" class="docs-object-method">&nbsp;</a> 
```python
summary(self, top=8): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/ConceptualDFT/FukuiFunctions.py#L181)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/ConceptualDFT/FukuiFunctions.py#L181?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/ElectronicStructure/ConceptualDFT/FukuiFunctions.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/ElectronicStructure/ConceptualDFT/FukuiFunctions.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/ElectronicStructure/ConceptualDFT/FukuiFunctions.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/ElectronicStructure/ConceptualDFT/FukuiFunctions.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/ConceptualDFT.py#L102?message=Update%20Docs)   
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