## <a id="Psience.Modes.MixtureModes.MixtureModes">MixtureModes</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes.py#L16)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes.py#L16?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes.py#L22)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes.py#L22?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.to_state" class="docs-object-method">&nbsp;</a> 
```python
to_state(self, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L68)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L68?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.from_state" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_state(cls, data, serializer=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L79)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L79?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.prep_modes" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
prep_modes(cls, modes): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L92)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L92?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L144)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L144?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L169)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L169?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.rotate" class="docs-object-method">&nbsp;</a> 
```python
rotate(self, rot, in_place=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L194)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L194?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.transform" class="docs-object-method">&nbsp;</a> 
```python
transform(self, tf, inv=None, origin=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L197)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L197?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.cartesian_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
cartesian_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L226)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L226?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.embed_coords" class="docs-object-method">&nbsp;</a> 
```python
embed_coords(self, carts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L230)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L230?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.unembed_coords" class="docs-object-method">&nbsp;</a> 
```python
unembed_coords(self, mode_coords): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L236)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L236?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.total_transformation" class="docs-object-method">&nbsp;</a> 
```python
@property
total_transformation(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L244)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L244?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.inverse_transformation" class="docs-object-method">&nbsp;</a> 
```python
@property
inverse_transformation(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L247)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L247?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.embed_derivs" class="docs-object-method">&nbsp;</a> 
```python
embed_derivs(self, derivs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L254)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L254?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.unembed_derivs" class="docs-object-method">&nbsp;</a> 
```python
unembed_derivs(self, derivs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L256)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L256?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.is_cartesian" class="docs-object-method">&nbsp;</a> 
```python
@property
is_cartesian(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L259)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L259?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.coords_by_modes" class="docs-object-method">&nbsp;</a> 
```python
@property
coords_by_modes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L265)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L265?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.modes_by_coords" class="docs-object-method">&nbsp;</a> 
```python
@property
modes_by_coords(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L268)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L268?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.compute_local_hessian" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
compute_local_hessian(cls, f, g): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L300)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L300?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.compute_local_gmatrix" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
compute_local_gmatrix(cls, f, g): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L305)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L305?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.compute_hessian" class="docs-object-method">&nbsp;</a> 
```python
compute_hessian(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L310)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L310?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.compute_gmatrix" class="docs-object-method">&nbsp;</a> 
```python
compute_gmatrix(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L313)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L313?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.compute_freqs" class="docs-object-method">&nbsp;</a> 
```python
compute_freqs(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L323)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L323?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.local_hessian" class="docs-object-method">&nbsp;</a> 
```python
@property
local_hessian(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L329)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L329?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.local_gmatrix" class="docs-object-method">&nbsp;</a> 
```python
@property
local_gmatrix(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L336)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L336?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.local_freqs" class="docs-object-method">&nbsp;</a> 
```python
@property
local_freqs(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L343)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L343?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.get_nearest_mode_transform" class="docs-object-method">&nbsp;</a> 
```python
get_nearest_mode_transform(self, alternate_modes: numpy.ndarray, mass_weighted=None, atoms=None, maximum_similarity=True, unitarize=True, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L349)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L349?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.get_projected_localized_mode_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_projected_localized_mode_transformation(self, projectors, masses=None, origin=None, localization_type=None, allow_mode_mixing=False, maximum_similarity=False, unitarize=True, zero_freq_cutoff=None, orthogonal_projection=False, atoms=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L397)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L397?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.get_atom_localized_mode_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_atom_localized_mode_transformation(self, atoms, masses=None, origin=None, localization_type='ned', allow_mode_mixing=False, maximum_similarity=False, unitarize=True, zero_freq_cutoff=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L493)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L493?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.get_coordinate_projected_localized_mode_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_coordinate_projected_localized_mode_transformation(self, coordinate_constraints, atoms=None, masses=None, origin=None, localization_type='ned', allow_mode_mixing=False, maximum_similarity=False, orthogonal_projection=True, unitarize=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L529)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L529?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.get_internal_localized_mode_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_internal_localized_mode_transformation(self, expansion_coordinates: 'Iterable[Iterable[int]|dict]', fixed_atoms=None, mass_weighted=False, project_transrot=True, atoms=None, maximum_similarity=False, orthogonal_projection=False, projection=False, allow_mode_mixing=False, unitarize=True, origin=None, masses=None, localization_type='ned'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L620)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L620?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.get_displacement_localized_mode_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_displacement_localized_mode_transformation(self, mode_blocks=None, atoms=None, mass_weighted=True, unitarize=True, **maximizer_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L737)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L737?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.get_mass_scaled_mode_transformation" class="docs-object-method">&nbsp;</a> 
```python
get_mass_scaled_mode_transformation(self, mass_scaling, *, atoms, localization_cutoff=0.8, num_modes=None, project_transrot=False, unitarize=True, **diag_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L770)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L770?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.localizer_dispatch" class="docs-object-method">&nbsp;</a> 
```python
@property
localizer_dispatch(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L846)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L846?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.localize" class="docs-object-method">&nbsp;</a> 
```python
localize(self, method=None, *, atoms=None, target_modes=None, internals=None, mode_blocks=None, coordinate_constraints=None, projections=None, reorthogonalize=None, mass_scaling=None, unitarize=True, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L858)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L858?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.make_mass_weighted" class="docs-object-method">&nbsp;</a> 
```python
make_mass_weighted(self, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L935)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L935?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.remove_mass_weighting" class="docs-object-method">&nbsp;</a> 
```python
remove_mass_weighting(self, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L948)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L948?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.make_frequency_scaled" class="docs-object-method">&nbsp;</a> 
```python
make_frequency_scaled(self, freqs=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L969)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L969?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.remove_frequency_scaling" class="docs-object-method">&nbsp;</a> 
```python
remove_frequency_scaling(self, freqs=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L981)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L981?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.make_dimensionless" class="docs-object-method">&nbsp;</a> 
```python
make_dimensionless(self, freqs=None, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L993)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L993?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.make_dimensioned" class="docs-object-method">&nbsp;</a> 
```python
make_dimensioned(self, freqs=None, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1000)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1000?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.apply_projection" class="docs-object-method">&nbsp;</a> 
```python
apply_projection(self, proj, project_transrot=True, masses=None, origin=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1020)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1020?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.apply_constraints" class="docs-object-method">&nbsp;</a> 
```python
apply_constraints(self, coordinate_constraints, atoms=None, masses=None, origin=None, orthogonal_projection=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1061)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1061?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.apply_transformation" class="docs-object-method">&nbsp;</a> 
```python
apply_transformation(self, tf, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1126)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1126?message=Update%20Docs)]
</div>


<a id="Psience.Modes.MixtureModes.MixtureModes.make_oblique" class="docs-object-method">&nbsp;</a> 
```python
make_oblique(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/MixtureModes/MixtureModes.py#L1131)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes/MixtureModes.py#L1131?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/MixtureModes.py#L16?message=Update%20Docs)   
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