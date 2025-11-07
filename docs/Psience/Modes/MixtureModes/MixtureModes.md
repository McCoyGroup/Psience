## <a id="Psience.Modes.MixtureModes.MixtureModes">MixtureModes</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes.py#L17)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes.py#L17?message=Update%20Docs)]
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
localization_type: str
zero_freq_cutoff: float
LocalizationMethods: LocalizationMethods
```
<a id="Psience.Modes.MixtureModes.MixtureModes.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, basis, coeffs, freqs=None, origin=None, masses=None, inverse=None, mass_weighted=False, frequency_scaled=False, g_matrix=None, name=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes.py#L23)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes.py#L23?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.to_state" class="docs-object-method">&nbsp;</a> 
```python
to_state(self, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L69)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L69?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.from_state" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_state(cls, data, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L80)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L80?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.prep_modes" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
prep_modes(cls, modes): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L93)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L93?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L145)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L145?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L173)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L173?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.rotate" class="docs-object-method">&nbsp;</a> 
```python
rotate(self, rot, in_place=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L198)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L198?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.transform" class="docs-object-method">&nbsp;</a> 
```python
transform(self, tf, inv=None, origin=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L201)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L201?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.cartesian_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
cartesian_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L230)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L230?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.embed_coords" class="docs-object-method">&nbsp;</a> 
```python
embed_coords(self, carts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L234)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L234?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.unembed_coords" class="docs-object-method">&nbsp;</a> 
```python
unembed_coords(self, mode_coords): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L240)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L240?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.total_transformation" class="docs-object-method">&nbsp;</a> 
```python
@property
total_transformation(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L248)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L248?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.inverse_transformation" class="docs-object-method">&nbsp;</a> 
```python
@property
inverse_transformation(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L251)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L251?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.embed_derivs" class="docs-object-method">&nbsp;</a> 
```python
embed_derivs(self, derivs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L258)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L258?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.unembed_derivs" class="docs-object-method">&nbsp;</a> 
```python
unembed_derivs(self, derivs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L260)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L260?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.is_cartesian" class="docs-object-method">&nbsp;</a> 
```python
@property
is_cartesian(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L263)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L263?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.coords_by_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
coords_by_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L269)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L269?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.modes_by_coords" class="docs-object-method">&nbsp;</a> 
```python
@property
modes_by_coords(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L272)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L272?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.compute_local_transformations" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
compute_local_transformations(cls, f, g): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L310)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L310?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.compute_local_hessian" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
compute_local_hessian(cls, f, g): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L317)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L317?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.compute_local_gmatrix" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
compute_local_gmatrix(cls, f, g): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L322)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L322?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.compute_hessian" class="docs-object-method">&nbsp;</a> 
```python
compute_hessian(self, system='modes'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L327)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L327?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.compute_gmatrix" class="docs-object-method">&nbsp;</a> 
```python
compute_gmatrix(self, system='modes', return_fractional=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L337)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L337?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.compute_freqs" class="docs-object-method">&nbsp;</a> 
```python
compute_freqs(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L362)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L362?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.local_hessian" class="docs-object-method">&nbsp;</a> 
```python
@property
local_hessian(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L368)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L368?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.local_gmatrix" class="docs-object-method">&nbsp;</a> 
```python
@property
local_gmatrix(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L375)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L375?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.local_freqs" class="docs-object-method">&nbsp;</a> 
```python
@property
local_freqs(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L382)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L382?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.local_mode_transformations" class="docs-object-method">&nbsp;</a> 
```python
@property
local_mode_transformations(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L386)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L386?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.get_nearest_mode_transform" class="docs-object-method">&nbsp;</a> 
```python
get_nearest_mode_transform(self, alternate_modes: numpy.ndarray, mass_weighted=None, atoms=None, maximum_similarity=True, unitarize=True, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L393)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L393?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.get_projected_localized_mode_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_projected_localized_mode_transformation(self, projectors, masses=None, origin=None, localization_type=None, allow_mode_mixing=False, maximum_similarity=False, unitarize=True, zero_freq_cutoff=None, orthogonal_projection=False, atoms=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L441)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L441?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.get_atom_localized_mode_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_atom_localized_mode_transformation(self, atoms, masses=None, origin=None, localization_type='ned', allow_mode_mixing=False, maximum_similarity=False, orthogonal_projection=False, unitarize=True, zero_freq_cutoff=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L545)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L545?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.get_coordinate_projected_localized_mode_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_coordinate_projected_localized_mode_transformation(self, coordinate_constraints, atoms=None, masses=None, origin=None, localization_type='ned', allow_mode_mixing=False, maximum_similarity=False, orthogonal_projection=True, unitarize=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L584)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L584?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.get_internal_localized_mode_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_internal_localized_mode_transformation(self, expansion_coordinates: 'Iterable[Iterable[int]|dict]', fixed_atoms=None, mass_weighted=False, project_transrot=True, atoms=None, maximum_similarity=False, orthogonal_projection=False, projection=False, allow_mode_mixing=False, unitarize=True, origin=None, masses=None, localization_type='ned'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L675)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L675?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.get_displacement_localized_mode_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_displacement_localized_mode_transformation(self, mode_blocks=None, atoms=None, mass_weighted=True, unitarize=True, **maximizer_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L792)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L792?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.get_mass_scaled_mode_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_mass_scaled_mode_transformation(self, mass_scaling, *, atoms, localization_cutoff=0.8, num_modes=None, project_transrot=False, unitarize=True, **diag_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L825)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L825?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.localizer_dispatch" class="docs-object-method">&nbsp;</a> 
```python
@property
localizer_dispatch(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L901)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L901?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.localize" class="docs-object-method">&nbsp;</a> 
```python
localize(self, method=None, *, atoms=None, target_modes=None, internals=None, mode_blocks=None, coordinate_constraints=None, projections=None, reorthogonalize=None, mass_scaling=None, unitarize=True, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L913)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L913?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.make_mass_weighted" class="docs-object-method">&nbsp;</a> 
```python
make_mass_weighted(self, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L991)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L991?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.remove_mass_weighting" class="docs-object-method">&nbsp;</a> 
```python
remove_mass_weighting(self, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1004)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1004?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.make_frequency_scaled" class="docs-object-method">&nbsp;</a> 
```python
make_frequency_scaled(self, freqs=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1025)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1025?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.remove_frequency_scaling" class="docs-object-method">&nbsp;</a> 
```python
remove_frequency_scaling(self, freqs=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1037)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1037?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.make_dimensionless" class="docs-object-method">&nbsp;</a> 
```python
make_dimensionless(self, freqs=None, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1049)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1049?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.make_dimensioned" class="docs-object-method">&nbsp;</a> 
```python
make_dimensioned(self, freqs=None, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1056)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1056?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.apply_projection" class="docs-object-method">&nbsp;</a> 
```python
apply_projection(self, proj, project_transrot=True, masses=None, origin=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1076)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1076?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.apply_constraints" class="docs-object-method">&nbsp;</a> 
```python
apply_constraints(self, coordinate_constraints, atoms=None, masses=None, origin=None, orthogonal_projection=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1117)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1117?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.apply_transformation" class="docs-object-method">&nbsp;</a> 
```python
apply_transformation(self, tf, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1182)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1182?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.make_oblique" class="docs-object-method">&nbsp;</a> 
```python
make_oblique(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1187)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1187?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes.py#L17?message=Update%20Docs)   
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