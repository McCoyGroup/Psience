## <a id="Psience.Vibronic.FCFs.FranckCondonModel">FranckCondonModel</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Vibronic/FCFs.py#L56)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Vibronic/FCFs.py#L56?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
default_rotation_method: str
default_rotation_order: State
default_rotation_center: State
default_embedding_ref: State
default_include_rotation: bool
OverlapData: OverlapData
Embedding: Embedding
duschinsky_cutoff: float
evaluator_plans: dict
integral_block_size: int
state_space_prep_registry: dict
```
<a id="Psience.Vibronic.FCFs.FranckCondonModel.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, gs_nms: Psience.Modes.NormalModes.NormalModes, es_nms, atoms=None, *, logger=None, embed=True, embedding_ref=None, masses=None, mass_weight=True, dimensionless=False, mode_selection=None, mode_reordering=None, rotation_method=None, rotation_order=None, rotation_center=None, include_rotation=None, rotation_blocks=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Vibronic/FCFs.py#L74)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Vibronic/FCFs.py#L74?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.FCFs.FranckCondonModel.prep_modes" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
prep_modes(cls, gs_nms, es_nms, embed=True, embedding_ref=None, masses=None, mass_weight=False, dimensionless=False, mode_selection=None, mode_reordering=None, **rotation_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L108)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L108?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.FCFs.FranckCondonModel.from_files" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_files(cls, gs_file, es_file, logger=None, mode_selection=None, mode_reordering=None, internals=None, internals_ref='gs', **rotation_embedding_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L141)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L141?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.FCFs.FranckCondonModel.convert_internal_modes" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
convert_internal_modes(cls, mol, nms): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L166)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L166?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.FCFs.FranckCondonModel.from_mols" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_mols(cls, gs, es, logger=None, remove_transrot=True, use_internals=True, embed=True, mass_weight=True, **rotation_embedding_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L187)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L187?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.FCFs.FranckCondonModel.get_overlaps" class="docs-object-method">&nbsp;</a> 
```python
get_overlaps(self, excitations, *, duschinsky_cutoff=None, ground_states=None, return_states=True, **rotation_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Vibronic/FCFs/FranckCondonModel.py#L221)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Vibronic/FCFs/FranckCondonModel.py#L221?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.FCFs.FranckCondonModel.get_overlap_data" class="docs-object-method">&nbsp;</a> 
```python
get_overlap_data(self, **rotation_embedding_opts) -> Psience.Vibronic.FCFs.OverlapData: 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Vibronic/FCFs/FranckCondonModel.py#L251)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Vibronic/FCFs/FranckCondonModel.py#L251?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.FCFs.FranckCondonModel.prep_overlap_args" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
prep_overlap_args(self, gs_nms, es_nms): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L271)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L271?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.FCFs.FranckCondonModel.get_poly_evaluation_plan" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_poly_evaluation_plan(self, exponents, alphas=None, zpe_prod=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L294)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L294?message=Update%20Docs)]
</div>
Provides a function that can take a set of indices and rotation matrices
from the gs and es bases to the shared basis of the central Gaussian and compute
the corresponding term contributions by considering every even
permutation of indices that could lead to a non-zero contribution
  - `tg`: `Any`
    > 
  - `te`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Vibronic.FCFs.FranckCondonModel.zero_point_alpha_contrib" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
zero_point_alpha_contrib(cls, alphas): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L363)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L363?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.FCFs.FranckCondonModel.term_evaluator" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
term_evaluator(self, exponents_list, splits_list, splits_inds_list, weights_list, gammas, alphas=None, zpe_prod=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L366)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L366?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.FCFs.FranckCondonModel.evaluate_poly_chunks" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
evaluate_poly_chunks(cls, poly_coeffs, exps, splits, split_inds, weights, alphas, include_baseline=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L452)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L452?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.FCFs.FranckCondonModel.evaluate_poly_contrib_chunk" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
evaluate_poly_contrib_chunk(self, inds, exponents_list, splits_list, splits_inds_list, weights_list, alphas, coeffs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L471)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L471?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.FCFs.FranckCondonModel.evaluate_shifted_poly_overlap" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
evaluate_shifted_poly_overlap(self, poly: 'HermiteProductPolynomial', Q, alphas, zpe_prod, duschinsky_cutoff=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L514)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L514?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.FCFs.FranckCondonModel.df_weights" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
df_weights(cls, n): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L572)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L572?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.FCFs.FranckCondonModel.get_overlap_gaussian_data" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_overlap_gaussian_data(cls, freqs_gs, modes_gs, inv_gs, center_gs, freqs_es, modes_es, inv_es, center_es, rotation_method=None, rotation_order=None, rotation_center=None, include_rotation=None, rotation_blocks=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L581)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L581?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.FCFs.FranckCondonModel.eval_fcf_overlaps" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
eval_fcf_overlaps(self, excitations_gs, freqs_gs, modes_gs, inv_gs, center_gs, excitations_es, freqs_es, modes_es, inv_es, center_es, duschinsky_cutoff=None, logger=None, **rotation_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L703)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L703?message=Update%20Docs)]
</div>
Evaluates the Gaussian overlaps between two H.O. wave functions defined by
a set of polynomial coefficients, broadening factors, and centers, assuming
the modes and centers are in an Eckart fream


<a id="Psience.Vibronic.FCFs.FranckCondonModel.embed_modes" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
embed_modes(cls, gs_nms: 'NormalModes', es_nms, ref=None, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L838)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L838?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.FCFs.FranckCondonModel.mass_weight_nms" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
mass_weight_nms(cls, nms, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L928)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L928?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.FCFs.FranckCondonModel.make_dimensionless" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
make_dimensionless(cls, nms, freqs=None, masses=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L934)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L934?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.FCFs.FranckCondonModel.prep_states_from_threshold_and_quanta" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
prep_states_from_threshold_and_quanta(cls, nms, *, threshold=None, min_freq=None, max_state=None, min_quanta=None, max_quanta=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L940)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L940?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.FCFs.FranckCondonModel.prep_states_from_excitations" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
prep_states_from_excitations(cls, nms, *, states, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L963)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L963?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.FCFs.FranckCondonModel.state_space_prep_dispatchers" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
state_space_prep_dispatchers(cls): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L973)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L973?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.FCFs.FranckCondonModel.dispatch_state_space_prep" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
dispatch_state_space_prep(cls, spec, nms): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L982)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L982?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.FCFs.FranckCondonModel.prep_state_space" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
prep_state_space(cls, excitations, nms, check=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1012)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1012?message=Update%20Docs)]
</div>
Dispatcher to get appropriate state spaces
  - `excitations`: `Any`
    > 
  - `check`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Vibronic.FCFs.FranckCondonModel.get_fcfs" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_fcfs(cls, gs_nms: 'NormalModes', es_nms: 'NormalModes', excitations, ground_states=None, duschinsky_cutoff=None, logger=None, **rotation_embedding_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1056)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1056?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.FCFs.FranckCondonModel.format_overlap_tables" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
format_overlap_tables(cls, es, overlaps, include_headers=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1082)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1082?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.FCFs.FranckCondonModel.get_fcf_spectrum" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_fcf_spectrum(self, gs_nms: 'NormalModes', es_nms: 'NormalModes', excitations, ground_states=None, logger=None, duschinsky_cutoff=None, return_states=False, **rotation_embedding_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L1093)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L1093?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.FCFs.FranckCondonModel.prep_opts" class="docs-object-method">&nbsp;</a> 
```python
prep_opts(self, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Vibronic/FCFs/FranckCondonModel.py#L1141)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Vibronic/FCFs/FranckCondonModel.py#L1141?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.FCFs.FranckCondonModel.get_spectrum" class="docs-object-method">&nbsp;</a> 
```python
get_spectrum(self, excitations, *, ground_states=None, return_states=False, duschinsky_cutoff=None, **rotation_embedding_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Vibronic/FCFs/FranckCondonModel.py#L1160)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Vibronic/FCFs/FranckCondonModel.py#L1160?message=Update%20Docs)]
</div>


<a id="Psience.Vibronic.FCFs.FranckCondonModel.get_ezFCF_input" class="docs-object-method">&nbsp;</a> 
```python
get_ezFCF_input(self, excitations, atoms=None, ground_states=None, **rotation_embedding_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Vibronic/FCFs/FranckCondonModel.py#L1178)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Vibronic/FCFs/FranckCondonModel.py#L1178?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Vibronic/FCFs/FranckCondonModel.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Vibronic/FCFs/FranckCondonModel.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Vibronic/FCFs/FranckCondonModel.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Vibronic/FCFs/FranckCondonModel.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Vibronic/FCFs.py#L56?message=Update%20Docs)   
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