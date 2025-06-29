## <a id="Psience.Modes.NormalModes.NormalModes">NormalModes</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/NormalModes.py#L19)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/NormalModes.py#L19?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
name: str
ModeData: ModeData
default_zero_freq_cutoff: float
default_projected_zero_freq_cutoff: NoneType
localization_type: str
zero_freq_cutoff: float
LocalizationMethods: LocalizationMethods
```
<a id="Psience.Modes.NormalModes.NormalModes.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, basis, coeffs, freqs=None, origin=None, masses=None, inverse=None, name=None, mass_weighted=False, frequency_scaled=False, g_matrix=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/NormalModes.py#L22)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/NormalModes.py#L22?message=Update%20Docs)]
</div>


<a id="Psience.Modes.NormalModes.NormalModes.get_normal_modes" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_normal_modes(cls, f_matrix, mass_spec, remove_transrot=True, dimensionless=False, mass_weighted=None, zero_freq_cutoff=None, return_gmatrix=False, projector=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L49)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L49?message=Update%20Docs)]
</div>


<a id="Psience.Modes.NormalModes.NormalModes.from_fg" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_fg(cls, basis, f_matrix, mass_spec, remove_transrot=True, dimensionless=False, zero_freq_cutoff=None, mass_weighted=None, origin=None, projector=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L232)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L232?message=Update%20Docs)]
</div>
Generates normal modes from the specified F and G matrices
  - `basis`: `Any`
    > 
  - `f_matrix`: `Any`
    > second derivatives of the potential
  - `mass_spec`: `Any`
    > 
  - `mass_units`: `Any`
    > 
  - `remove_transrot`: `Any`
    > 
  - `opts`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Modes.NormalModes.NormalModes.from_molecule" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_molecule(cls, mol, dimensionless=False, use_internals=None, potential_derivatives=None, project_transrot=True, zero_freq_cutoff=None, masses=None, energy_evaluator=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L281)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L281?message=Update%20Docs)]
</div>


<a id="Psience.Modes.NormalModes.NormalModes.get_nearest_mode_transform" class="docs-object-method">&nbsp;</a> 
```python
get_nearest_mode_transform(self, alternate_modes: numpy.ndarray, mass_weighted=None, atoms=None, maximum_similarity=True, unitarize=True, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/NormalModes/NormalModes.py#L352)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/NormalModes/NormalModes.py#L352?message=Update%20Docs)]
</div>


<a id="Psience.Modes.NormalModes.NormalModes.get_projected_localized_mode_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_projected_localized_mode_transformation(self, projectors, masses=None, origin=None, localization_type=None, allow_mode_mixing=False, maximum_similarity=False, unitarize=True, zero_freq_cutoff=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/NormalModes/NormalModes.py#L400)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/NormalModes/NormalModes.py#L400?message=Update%20Docs)]
</div>


<a id="Psience.Modes.NormalModes.NormalModes.get_atom_localized_mode_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_atom_localized_mode_transformation(self, atoms, masses=None, origin=None, localization_type='ned', allow_mode_mixing=False, maximum_similarity=False, unitarize=True, zero_freq_cutoff=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/NormalModes/NormalModes.py#L484)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/NormalModes/NormalModes.py#L484?message=Update%20Docs)]
</div>


<a id="Psience.Modes.NormalModes.NormalModes.get_coordinate_projected_localized_mode_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_coordinate_projected_localized_mode_transformation(self, coordinate_constraints, atoms=None, masses=None, origin=None, localization_type='ned', allow_mode_mixing=False, maximum_similarity=False, orthogonal_projection=True, unitarize=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/NormalModes/NormalModes.py#L520)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/NormalModes/NormalModes.py#L520?message=Update%20Docs)]
</div>


<a id="Psience.Modes.NormalModes.NormalModes.get_internal_localized_mode_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_internal_localized_mode_transformation(self, expansion_coordinates: 'Iterable[Iterable[int]|dict]', fixed_atoms=None, mass_weighted=False, project_transrot=True, atoms=None, maximum_similarity=False, orthogonal_projection=False, projection=False, allow_mode_mixing=False, unitarize=True, origin=None, masses=None, localization_type='ned'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/NormalModes/NormalModes.py#L611)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/NormalModes/NormalModes.py#L611?message=Update%20Docs)]
</div>


<a id="Psience.Modes.NormalModes.NormalModes.get_displacement_localized_mode_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_displacement_localized_mode_transformation(self, mode_blocks=None, atoms=None, mass_weighted=True, unitarize=True, **maximizer_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/NormalModes/NormalModes.py#L728)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/NormalModes/NormalModes.py#L728?message=Update%20Docs)]
</div>


<a id="Psience.Modes.NormalModes.NormalModes.localizer_dispatch" class="docs-object-method">&nbsp;</a> 
```python
@property
localizer_dispatch(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/NormalModes/NormalModes.py#L769)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/NormalModes/NormalModes.py#L769?message=Update%20Docs)]
</div>


<a id="Psience.Modes.NormalModes.NormalModes.localize" class="docs-object-method">&nbsp;</a> 
```python
localize(self, method=None, *, atoms=None, target_modes=None, internals=None, mode_blocks=None, coordinate_constraints=None, projections=None, reorthogonalize=None, unitarize=True, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/NormalModes/NormalModes.py#L780)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/NormalModes/NormalModes.py#L780?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Modes/NormalModes/NormalModes.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Modes/NormalModes/NormalModes.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Modes/NormalModes/NormalModes.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Modes/NormalModes/NormalModes.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/NormalModes.py#L19?message=Update%20Docs)   
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