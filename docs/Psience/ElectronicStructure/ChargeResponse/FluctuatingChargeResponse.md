## <a id="Psience.ElectronicStructure.ChargeResponse.FluctuatingChargeResponse">FluctuatingChargeResponse</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/ChargeResponse.py#L15)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/ChargeResponse.py#L15?message=Update%20Docs)]
</div>

For a model whose charges minimize `E(q) = chi . q + q . J q / 2` subject to a fixed total
charge, the induced dipole and the dipole polarizability follow from the electronegativities
`chi = dE/dq` and hardness matrix `J = d^2E/dq^2`. Given those (and their nuclear-coordinate
derivatives), `expand` returns the coordinate expansions of both.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.ElectronicStructure.ChargeResponse.FluctuatingChargeResponse.mu_factor" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
mu_factor(cls, total_charge, X, J_inv): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L23)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L23?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.ChargeResponse.FluctuatingChargeResponse.dipole_factor_expansion" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
dipole_factor_expansion(cls, coords, total_charge, X, J_inv): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L27)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L27?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.ChargeResponse.FluctuatingChargeResponse.polarizability_factor" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
polarizability_factor(cls, coords, J_inv): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L32)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L32?message=Update%20Docs)]
</div>


<a id="Psience.ElectronicStructure.ChargeResponse.FluctuatingChargeResponse.expand" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
expand(cls, total_charge, coords, potential_by_charges_by_coords): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L41)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L41?message=Update%20Docs)]
</div>

  - `total_charge`: `Any`
    > the total charge
  - `coords`: `Any`
    > nuclear coordinates, `(..., natoms, 3)`
  - `potential_by_charges_by_coords`: `Any`
    > `[[E, dE/dx, ...], [chi, dchi/dx, ...], [J, dJ/dx, ...]]`
    (charge-derivative order outermost, coordinate-derivative order innermost)
  - `:returns`: `list[list[np.ndarray]]`
    > `[dipole_expansion, polarizability_expansion]` (each a list over coordinate orders),
    truncated to the available charge-derivative orders


<a id="Psience.ElectronicStructure.ChargeResponse.FluctuatingChargeResponse.distorted_harmonic_polarizabilities" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
distorted_harmonic_polarizabilities(cls, charges, dipole_expansion, coords, hessians): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L110)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L110?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/ElectronicStructure/ChargeResponse/FluctuatingChargeResponse.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/ElectronicStructure/ChargeResponse/FluctuatingChargeResponse.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/ElectronicStructure/ChargeResponse/FluctuatingChargeResponse.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/ElectronicStructure/ChargeResponse/FluctuatingChargeResponse.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/ChargeResponse.py#L15?message=Update%20Docs)   
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