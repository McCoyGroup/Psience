## <a id="Psience.AnalyticModels.AnalyticModelConstructors.MolecularModel">MolecularModel</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors.py#L1337)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors.py#L1337?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.AnalyticModels.AnalyticModelConstructors.MolecularModel.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, mol, coords, potential, dipole=None, values=None, rotation=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/MolecularModel.py#L1338)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/MolecularModel.py#L1338?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.MolecularModel.potential" class="docs-object-method">&nbsp;</a> 
```python
@property
potential(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/MolecularModel.py#L1341)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/MolecularModel.py#L1341?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.MolecularModel.gmatrix" class="docs-object-method">&nbsp;</a> 
```python
@property
gmatrix(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/MolecularModel.py#L1344)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/MolecularModel.py#L1344?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.MolecularModel.vprime" class="docs-object-method">&nbsp;</a> 
```python
@property
vprime(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/MolecularModel.py#L1347)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/MolecularModel.py#L1347?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.MolecularModel.dipole" class="docs-object-method">&nbsp;</a> 
```python
@property
dipole(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/MolecularModel.py#L1350)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/MolecularModel.py#L1350?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.MolecularModel.setup_AIMD" class="docs-object-method">&nbsp;</a> 
```python
setup_AIMD(self, timestep=0.5, initial_energies=None, initial_displacements=None, displaced_coords=None, track_kinetic_energy=False, track_velocities=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/MolecularModel.py#L1354)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/MolecularModel.py#L1354?message=Update%20Docs)]
</div>


<a id="Psience.AnalyticModels.AnalyticModelConstructors.MolecularModel.setup_DGB" class="docs-object-method">&nbsp;</a> 
```python
setup_DGB(self, centers, *, masses=None, modes='normal', transformations=None, alphas='auto', cartesians=None, potential_function=None, dipole_function=None, optimize_centers=None, quadrature_degree=None, expansion_degree=None, pairwise_potential_functions=None, internals=False, logger=True, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/AnalyticModels/AnalyticModelConstructors/MolecularModel.py#L1377)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors/MolecularModel.py#L1377?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/AnalyticModels/AnalyticModelConstructors/MolecularModel.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/AnalyticModels/AnalyticModelConstructors/MolecularModel.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/AnalyticModels/AnalyticModelConstructors/MolecularModel.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/AnalyticModels/AnalyticModelConstructors/MolecularModel.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/AnalyticModels/AnalyticModelConstructors.py#L1337?message=Update%20Docs)   
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