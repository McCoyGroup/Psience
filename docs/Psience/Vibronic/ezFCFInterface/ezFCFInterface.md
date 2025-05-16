## <a id="Psience.Vibronic.ezFCFInterface.ezFCFInterface">ezFCFInterface</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Vibronic/ezFCFInterface.py#L147)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Vibronic/ezFCFInterface.py#L147?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
ezFCFRunner: ezFCFRunner
```
<a id="Psience.Vibronic.ezFCFInterface.ezFCFInterface.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, atoms, gs_nms, es_nms, excitations, masses=None, ground_states=None, include_rotation=True, rotation_order='gs', rotation_blocks=None, rotation_center=None, logger=None, mode_reordering=None, rotation_method='duschinsky', mass_weight=False, dimensionless=False, always_run_parallel=True, print_all=True, embed=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Vibronic/ezFCFInterface.py#L148)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Vibronic/ezFCFInterface.py#L148?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.ezFCFInterface.ezFCFInterface.format" class="docs-object-method">&nbsp;</a> 
```python
format(self, job_type='harmonic_pes', temperature=0, spectrum_intensity_threshold=1e-08): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Vibronic/ezFCFInterface/ezFCFInterface.py#L196)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Vibronic/ezFCFInterface/ezFCFInterface.py#L196?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.ezFCFInterface.ezFCFInterface.prep_masses" class="docs-object-method">&nbsp;</a> 
```python
prep_masses(self, atoms, masses, units=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Vibronic/ezFCFInterface/ezFCFInterface.py#L235)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Vibronic/ezFCFInterface/ezFCFInterface.py#L235?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.ezFCFInterface.ezFCFInterface.format_masses_file" class="docs-object-method">&nbsp;</a> 
```python
format_masses_file(self, atom_map): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Vibronic/ezFCFInterface/ezFCFInterface.py#L273)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Vibronic/ezFCFInterface/ezFCFInterface.py#L273?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.ezFCFInterface.ezFCFInterface.run" class="docs-object-method">&nbsp;</a> 
```python
run(self, ezFCF_binary, dir=None, dir_prefix=None, dir_suffix=None, mode='w', prefix='ezFCF-', suffix='.xml', delete=True, raise_errors=True, **job_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Vibronic/ezFCFInterface/ezFCFInterface.py#L304)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Vibronic/ezFCFInterface/ezFCFInterface.py#L304?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.ezFCFInterface.ezFCFInterface.canonicalize_excitation_options" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
canonicalize_excitation_options(cls, nms, threshold=None, fixed_modes=None, ground_states=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L357)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L357?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.ezFCFInterface.ezFCFInterface.prep_excitations" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
prep_excitations(cls, nms, ground_states, excitations): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L410)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L410?message=Update%20Docs)]
</div>
Dispatcher to get appropriate state spaces
  - `excitations`: `Any`
    > 
  - `check`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Vibronic.ezFCFInterface.ezFCFInterface.prep_state_str" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
prep_state_str(cls, nms, state): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L508)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L508?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.ezFCFInterface.ezFCFInterface.prep_parallel" class="docs-object-method">&nbsp;</a> 
```python
prep_parallel(self, *elems, rotation_order='gs', max_vibr_excitations_in_initial_el_state=0, max_vibr_excitations_in_target_el_state=0, **etc): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Vibronic/ezFCFInterface/ezFCFInterface.py#L520)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Vibronic/ezFCFInterface/ezFCFInterface.py#L520?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.ezFCFInterface.ezFCFInterface.prep_duschinsky" class="docs-object-method">&nbsp;</a> 
```python
prep_duschinsky(self, *elems, rotation_order='gs', max_vibr_excitations_in_initial_el_state=0, max_vibr_excitations_in_target_el_state=0, **etc): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Vibronic/ezFCFInterface/ezFCFInterface.py#L539)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Vibronic/ezFCFInterface/ezFCFInterface.py#L539?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.ezFCFInterface.ezFCFInterface.prep_initial_state" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
prep_initial_state(cls, atoms, nms, autoembed=False, excitation_energy_units='cm-1'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L557)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L557?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.ezFCFInterface.ezFCFInterface.prep_target_state" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
prep_target_state(cls, atoms, nms, autoembed=False, mode_reordering=None, excitation_energy_units='cm-1'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L574)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L574?message=Update%20Docs)]
</div>

  - `atoms`: `Any`
    > 
  - `nms`: `Any`
    > 
  - `excitation_energy_units`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Vibronic.ezFCFInterface.ezFCFInterface.prep_state" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
prep_state(cls, tag, atoms, nms, *elems): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L605)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L605?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.ezFCFInterface.ezFCFInterface.prep_geometry" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
prep_geometry(cls, atoms, nms, linear=False, units='BohrRadius'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L614)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L614?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.ezFCFInterface.ezFCFInterface.format_modes_block" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
format_modes_block(cls, modes): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L629)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L629?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.ezFCFInterface.ezFCFInterface.parse_modes_block" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
parse_modes_block(cls, modes_str, num_atoms=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L652)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L652?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.ezFCFInterface.ezFCFInterface.prep_normal_modes" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
prep_normal_modes(cls, atoms, nms): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L671)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L671?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.ezFCFInterface.ezFCFInterface.format_freqs_block" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
format_freqs_block(cls, freqs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L682)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L682?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.ezFCFInterface.ezFCFInterface.parse_freqs_block" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
parse_freqs_block(cls, freqs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L691)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L691?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.ezFCFInterface.ezFCFInterface.prep_frequencies" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
prep_frequencies(cls, nms): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L695)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L695?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.ezFCFInterface.ezFCFInterface.parse_state" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
parse_state(self, state_xml): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L701)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L701?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.ezFCFInterface.ezFCFInterface.parse_fc_model" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
parse_fc_model(cls, input_xml, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L760)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L760?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Vibronic/ezFCFInterface/ezFCFInterface.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Vibronic/ezFCFInterface/ezFCFInterface.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Vibronic/ezFCFInterface/ezFCFInterface.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Vibronic/ezFCFInterface/ezFCFInterface.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Vibronic/ezFCFInterface.py#L147?message=Update%20Docs)   
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