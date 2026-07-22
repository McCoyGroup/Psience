## <a id="Psience.Modes.MixtureModes.MixtureModes">MixtureModes</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes.py#L18)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes.py#L18?message=Update%20Docs)]
</div>

A `McUtils.Coordinerds.CoordinateSystem` object that expresses coordinates as
a rotation on some base set of coordinates with some associated frequencies.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
name: str
zero_freq_cutoff: float
ModeData: ModeData
default_zero_freq_cutoff: float
localization_type: str
localization_zero_freq_cutoff: float
LocalizationMethods: LocalizationMethods
localization_options: tuple
```
<a id="Psience.Modes.MixtureModes.MixtureModes.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, basis, coeffs, freqs=None, origin=None, masses=None, inverse=None, mass_weighted=False, frequency_scaled=False, g_matrix=None, name=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes.py#L24)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes.py#L24?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.to_state" class="docs-object-method">&nbsp;</a> 
```python
to_state(self, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L70)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L70?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.from_state" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_state(cls, data, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L81)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L81?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.prep_modes" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
prep_modes(cls, modes): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L94)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L94?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L146)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L146?message=Update%20Docs)]
</div>
Takes a slice of the modes
  - `item`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Modes.MixtureModes.MixtureModes.modify" class="docs-object-method">&nbsp;</a> 
```python
modify(self, matrix=None, *, freqs=None, origin=None, masses=None, inverse=None, name=None, mass_weighted=None, frequency_scaled=None, g_matrix=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L174)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L174?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.rotate" class="docs-object-method">&nbsp;</a> 
```python
rotate(self, rot, in_place=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L199)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L199?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.transform" class="docs-object-method">&nbsp;</a> 
```python
transform(self, tf, inv=None, origin=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L202)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L202?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.cartesian_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
cartesian_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L231)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L231?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.embed_coords" class="docs-object-method">&nbsp;</a> 
```python
embed_coords(self, carts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L235)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L235?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.unembed_coords" class="docs-object-method">&nbsp;</a> 
```python
unembed_coords(self, mode_coords): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L241)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L241?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.total_transformation" class="docs-object-method">&nbsp;</a> 
```python
@property
total_transformation(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L249)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L249?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.inverse_transformation" class="docs-object-method">&nbsp;</a> 
```python
@property
inverse_transformation(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L252)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L252?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.embed_derivs" class="docs-object-method">&nbsp;</a> 
```python
embed_derivs(self, derivs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L259)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L259?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.unembed_derivs" class="docs-object-method">&nbsp;</a> 
```python
unembed_derivs(self, derivs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L261)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L261?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.is_cartesian" class="docs-object-method">&nbsp;</a> 
```python
@property
is_cartesian(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L264)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L264?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.coords_by_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
coords_by_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L270)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L270?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.modes_by_coords" class="docs-object-method">&nbsp;</a> 
```python
@property
modes_by_coords(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L273)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L273?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.compute_local_transformations" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
compute_local_transformations(cls, f, g): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L311)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L311?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.compute_local_hessian" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
compute_local_hessian(cls, f, g): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L318)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L318?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.compute_local_gmatrix" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
compute_local_gmatrix(cls, f, g): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L323)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L323?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.compute_hessian" class="docs-object-method">&nbsp;</a> 
```python
compute_hessian(self, system='modes'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L328)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L328?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.compute_gmatrix" class="docs-object-method">&nbsp;</a> 
```python
compute_gmatrix(self, system='modes', return_fractional=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L338)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L338?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.compute_freqs" class="docs-object-method">&nbsp;</a> 
```python
compute_freqs(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L382)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L382?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.local_hessian" class="docs-object-method">&nbsp;</a> 
```python
@property
local_hessian(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L391)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L391?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.local_gmatrix" class="docs-object-method">&nbsp;</a> 
```python
@property
local_gmatrix(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L398)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L398?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.local_freqs" class="docs-object-method">&nbsp;</a> 
```python
@property
local_freqs(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L405)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L405?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.local_mode_transformations" class="docs-object-method">&nbsp;</a> 
```python
@property
local_mode_transformations(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L409)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L409?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.get_nearest_mode_transform" class="docs-object-method">&nbsp;</a> 
```python
get_nearest_mode_transform(self, alternate_modes: numpy.ndarray, mass_weighted=None, atoms=None, maximum_similarity=True, unitarize=True, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L416)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L416?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.get_normal_modes" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_normal_modes(cls, f_matrix, mass_spec, remove_transrot=True, dimensionless=False, mass_weighted=None, zero_freq_cutoff=None, return_gmatrix=False, projector=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L465)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L465?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.get_projected_localized_mode_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_projected_localized_mode_transformation(self, projectors, masses=None, origin=None, localization_type=None, allow_mode_mixing=False, maximum_similarity=False, unitarize=True, zero_freq_cutoff=None, orthogonal_projection=False, project_zero_gmatrix_modes=True, project_zero_gmatrix_cutoff=1e-08, atoms=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L649)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L649?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.get_atom_localized_mode_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_atom_localized_mode_transformation(self, atoms, masses=None, origin=None, localization_type='ned', allow_mode_mixing=False, maximum_similarity=False, orthogonal_projection=False, unitarize=True, zero_freq_cutoff=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L755)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L755?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.get_fragment_localized_mode_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_fragment_localized_mode_transformation(self, fragment, masses=None, origin=None, localization_type='ned', allow_mode_mixing=True, maximum_similarity=False, orthogonal_projection=False, unitarize=True, zero_freq_cutoff=None, **etc): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L794)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L794?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.get_coordinate_projected_localized_mode_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_coordinate_projected_localized_mode_transformation(self, coordinate_constraints, atoms=None, masses=None, origin=None, localization_type='ned', allow_mode_mixing=False, maximum_similarity=False, orthogonal_projection=True, unitarize=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L817)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L817?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.get_internal_localized_mode_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_internal_localized_mode_transformation(self, expansion_coordinates: 'Iterable[Iterable[int]|dict]', fixed_atoms=None, mass_weighted=False, project_transrot=True, atoms=None, maximum_similarity=False, orthogonal_projection=False, projection=False, allow_mode_mixing=False, unitarize=True, origin=None, masses=None, localization_type='ned'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L908)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L908?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.get_displacement_localized_mode_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_displacement_localized_mode_transformation(self, mode_blocks=None, atoms=None, mass_weighted=True, unitarize=True, **maximizer_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1027)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1027?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.get_mass_scaled_mode_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_mass_scaled_mode_transformation(self, mass_scaling, *, atoms, localization_cutoff=0.8, num_modes=None, project_transrot=False, unitarize=True, **diag_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1060)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1060?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.localizer_dispatch" class="docs-object-method">&nbsp;</a> 
```python
@property
localizer_dispatch(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1137)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1137?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.localize" class="docs-object-method">&nbsp;</a> 
```python
localize(self, method=None, *, atoms=None, fragment=None, target_modes=None, internals=None, mode_blocks=None, coordinate_constraints=None, projections=None, reorthogonalize=None, mass_scaling=None, unitarize=True, allow_mode_mixing=False, project_zero_gmatrix_modes=None, project_zero_gmatrix_cutoff=1e-08, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1169)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1169?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.make_mass_weighted" class="docs-object-method">&nbsp;</a> 
```python
make_mass_weighted(self, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1281)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1281?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.remove_mass_weighting" class="docs-object-method">&nbsp;</a> 
```python
remove_mass_weighting(self, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1294)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1294?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.make_frequency_scaled" class="docs-object-method">&nbsp;</a> 
```python
make_frequency_scaled(self, freqs=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1315)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1315?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.remove_frequency_scaling" class="docs-object-method">&nbsp;</a> 
```python
remove_frequency_scaling(self, freqs=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1327)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1327?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.make_dimensionless" class="docs-object-method">&nbsp;</a> 
```python
make_dimensionless(self, freqs=None, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1339)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1339?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.make_dimensioned" class="docs-object-method">&nbsp;</a> 
```python
make_dimensioned(self, freqs=None, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1346)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1346?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.apply_projection" class="docs-object-method">&nbsp;</a> 
```python
apply_projection(self, proj, project_transrot=True, masses=None, origin=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1366)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1366?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.apply_constraints" class="docs-object-method">&nbsp;</a> 
```python
apply_constraints(self, coordinate_constraints, atoms=None, masses=None, origin=None, orthogonal_projection=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1407)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1407?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.apply_transformation" class="docs-object-method">&nbsp;</a> 
```python
apply_transformation(self, tf, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1472)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1472?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.make_oblique" class="docs-object-method">&nbsp;</a> 
```python
make_oblique(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1477)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1477?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Modes/MixtureModes/MixtureModes.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Modes/MixtureModes/MixtureModes.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Modes/MixtureModes/MixtureModes.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Modes/MixtureModes/MixtureModes.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes.py#L18?message=Update%20Docs)   
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