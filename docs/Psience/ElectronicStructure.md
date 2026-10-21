# <a id="Psience.ElectronicStructure">Psience.ElectronicStructure</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/ElectronicStructure/__init__.py#L1)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/__init__.py#L1?message=Update%20Docs)]
</div>
    
Electronic-structure models and descriptors, kept separate from `Molecule`.

`ElectronicStructureBackend`s are the raw evaluation layer for MLIPs (AIMNet2, HIPNN, MACE,
MACE-POLAR, UMA, UPET, Orb), xTB, and PySCF: energies, charges, dipoles, and polarizabilities
with their analytic Cartesian derivatives, in each model's native units. The
`Psience.Molecools` evaluators are thin wrappers over them that add unit handling, finite
differences, embeddings, and optimization.

On top of the backends, `ElectronicStateModel`s predict atomic charge (and spin) densities for a
given charge and multiplicity, represented as atom-centered Gaussian-smeared multipole
expansions (`AtomCenteredDensity`, on `McUtils.Zachary.MultipoleExpansion`), from which spin
densities, Fukui functions, and related conceptual-DFT quantities are built.

### Members
<div class="container alert alert-secondary bg-light">
  <div class="row">
   <div class="col" markdown="1">
[ElectronicStructureBackend](ElectronicStructure/Backends/ElectronicStructureBackend.md)   
</div>
   <div class="col" markdown="1">
[FluctuatingChargeResponse](ElectronicStructure/ChargeResponse/FluctuatingChargeResponse.md)   
</div>
   <div class="col" markdown="1">
[AIMNet2Backend](ElectronicStructure/MLIPBackends/AIMNet2Backend.md)   
</div>
</div>
  <div class="row">
   <div class="col" markdown="1">
[HIPNNBackend](ElectronicStructure/MLIPBackends/HIPNNBackend.md)   
</div>
   <div class="col" markdown="1">
[ASECalculatorBackend](ElectronicStructure/MLIPBackends/ASECalculatorBackend.md)   
</div>
   <div class="col" markdown="1">
[MACEBackend](ElectronicStructure/MLIPBackends/MACEBackend.md)   
</div>
</div>
  <div class="row">
   <div class="col" markdown="1">
[MACEPolarBackend](ElectronicStructure/MLIPBackends/MACEPolarBackend.md)   
</div>
   <div class="col" markdown="1">
[UMABackend](ElectronicStructure/MLIPBackends/UMABackend.md)   
</div>
   <div class="col" markdown="1">
[UPETBackend](ElectronicStructure/MLIPBackends/UPETBackend.md)   
</div>
</div>
  <div class="row">
   <div class="col" markdown="1">
[OrbBackend](ElectronicStructure/MLIPBackends/OrbBackend.md)   
</div>
   <div class="col" markdown="1">
[XTBBackend](ElectronicStructure/QuantumChemistryBackends/XTBBackend.md)   
</div>
   <div class="col" markdown="1">
[PySCFBackend](ElectronicStructure/QuantumChemistryBackends/PySCFBackend.md)   
</div>
</div>
  <div class="row">
   <div class="col" markdown="1">
[AtomCenteredDensity](ElectronicStructure/Densities/AtomCenteredDensity.md)   
</div>
   <div class="col" markdown="1">
[ElectronicState](ElectronicStructure/StateModels/ElectronicState.md)   
</div>
   <div class="col" markdown="1">
[ElectronicStateModel](ElectronicStructure/StateModels/ElectronicStateModel.md)   
</div>
</div>
  <div class="row">
   <div class="col" markdown="1">
[MACEPolarStateModel](ElectronicStructure/StateModels/MACEPolarStateModel.md)   
</div>
   <div class="col" markdown="1">
[BackendChargeStateModel](ElectronicStructure/StateModels/BackendChargeStateModel.md)   
</div>
   <div class="col" markdown="1">
[ChargeEvaluatorStateModel](ElectronicStructure/StateModels/ChargeEvaluatorStateModel.md)   
</div>
</div>
  <div class="row">
   <div class="col" markdown="1">
[SpinDensity](ElectronicStructure/ConceptualDFT/SpinDensity.md)   
</div>
   <div class="col" markdown="1">
[FukuiFunctions](ElectronicStructure/ConceptualDFT/FukuiFunctions.md)   
</div>
   <div class="col" markdown="1">
   
</div>
</div>
</div>





## Examples


















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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/ElectronicStructure.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/ElectronicStructure.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/ElectronicStructure.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/ElectronicStructure.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/ElectronicStructure/__init__.py#L1?message=Update%20Docs)   
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