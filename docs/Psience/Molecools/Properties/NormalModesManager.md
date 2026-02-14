## <a id="Psience.Molecools.Properties.NormalModesManager">NormalModesManager</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties.py#L2204)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties.py#L2204?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
recalc_normal_mode_tolerance: float
```
<a id="Psience.Molecools.Properties.NormalModesManager.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, mol, normal_modes=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties.py#L2205)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties.py#L2205?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.NormalModesManager.from_data" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_data(cls, mol, data): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L2228)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L2228?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.NormalModesManager.set_molecule" class="docs-object-method">&nbsp;</a> 
```python
set_molecule(self, mol): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/NormalModesManager.py#L2252)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/NormalModesManager.py#L2252?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.NormalModesManager.get_modes" class="docs-object-method">&nbsp;</a> 
```python
get_modes(self, quiet=False, allow_compute=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/NormalModesManager.py#L2256)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/NormalModesManager.py#L2256?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.NormalModesManager.modes" class="docs-object-method">&nbsp;</a> 
```python
@property
modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/NormalModesManager.py#L2260)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/NormalModesManager.py#L2260?message=Update%20Docs)]
</div>

  - `:returns`: `MolecularVibrations`
    >


<a id="Psience.Molecools.Properties.NormalModesManager.construct_normal_modes" class="docs-object-method">&nbsp;</a> 
```python
construct_normal_modes(self, modes): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/NormalModesManager.py#L2278)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/NormalModesManager.py#L2278?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.NormalModesManager.load" class="docs-object-method">&nbsp;</a> 
```python
load(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/NormalModesManager.py#L2299)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/NormalModesManager.py#L2299?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.NormalModesManager.update" class="docs-object-method">&nbsp;</a> 
```python
update(self, modes): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/NormalModesManager.py#L2301)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/NormalModesManager.py#L2301?message=Update%20Docs)]
</div>

  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.NormalModesManager.load_normal_modes" class="docs-object-method">&nbsp;</a> 
```python
load_normal_modes(self, file=None, mode=None, rephase=True, recalculate=False, quiet=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/NormalModesManager.py#L2431)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/NormalModesManager.py#L2431?message=Update%20Docs)]
</div>
Loads potential derivatives from a file (or from `source_file` if set)
  - `file`: `Any`
    > 
  - `rephase`: `bool`
    > whether to rephase FChk normal modes or not
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.NormalModesManager.get_normal_modes" class="docs-object-method">&nbsp;</a> 
```python
get_normal_modes(self, quiet=False, compute_force_constants=True, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/NormalModesManager.py#L2522)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/NormalModesManager.py#L2522?message=Update%20Docs)]
</div>
Loads normal modes from file or calculates
from force constants
  - `kwargs`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.NormalModesManager.get_force_constants" class="docs-object-method">&nbsp;</a> 
```python
get_force_constants(self, compute_force_constants=True, quiet=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/NormalModesManager.py#L2550)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/NormalModesManager.py#L2550?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.NormalModesManager.get_dipole_derivative_based_rephasing" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_dipole_derivative_based_rephasing(cls, modes, analytic_dipoles, numerical_dipoles, strict=True, allow_swaps=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L2564)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L2564?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.NormalModesManager.get_fchk_normal_mode_rephasing" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_fchk_normal_mode_rephasing(cls, mol, modes, use_dipoles=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L2646)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L2646?message=Update%20Docs)]
</div>
Returns the necessary rephasing to make the numerical dipole derivatives
agree with the analytic dipole derivatives as pulled from a Gaussian FChk file
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.NormalModesManager.get_partial_cubic_based_rephasing" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_partial_cubic_based_rephasing(cls, modes, partial_cubics: numpy.ndarray, degenerate_freq_threshold=1e-06, equivalent_threshold=0.01, nonzero_threshold=4.5e-06, ignore_single_zeros=False, frequency_scale=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L2711)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L2711?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.NormalModesManager.apply_transformation" class="docs-object-method">&nbsp;</a> 
```python
apply_transformation(self, transf): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/NormalModesManager.py#L3194)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/NormalModesManager.py#L3194?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.Properties.NormalModesManager.insert_atoms" class="docs-object-method">&nbsp;</a> 
```python
insert_atoms(self, atoms, coords, where): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/NormalModesManager.py#L3204)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/NormalModesManager.py#L3204?message=Update%20Docs)]
</div>
Handles the insertion of new atoms into the structure
  - `atoms`: `tuple[str]`
    > 
  - `coords`: `CoordinateSet`
    > 
  - `where`: `tuple[int]`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Molecools.Properties.NormalModesManager.delete_atoms" class="docs-object-method">&nbsp;</a> 
```python
delete_atoms(self, where): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/Properties/NormalModesManager.py#L3227)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties/NormalModesManager.py#L3227?message=Update%20Docs)]
</div>
Handles the deletion from the structure
  - `atoms`: `tuple[str]`
    > 
  - `coords`: `CoordinateSet`
    > 
  - `where`: `tuple[int]`
    > 
  - `:returns`: `_`
    >
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Molecools/Properties/NormalModesManager.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Molecools/Properties/NormalModesManager.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Molecools/Properties/NormalModesManager.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Molecools/Properties/NormalModesManager.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/Properties.py#L2204?message=Update%20Docs)   
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