## <a id="Psience.Modes.NormalModes.ReactionPathModes">ReactionPathModes</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/NormalModes.py#L390)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/NormalModes.py#L390?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
zero_gradient_cutoff: float
```
<a id="Psience.Modes.NormalModes.ReactionPathModes.get_rp_modes" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_rp_modes(cls, gradient, f_matrix, mass_spec, remove_transrot=True, dimensionless=False, mass_weighted=None, zero_freq_cutoff=None, return_gmatrix=False, projector=None, zero_gradient_cutoff=None, use_max_gradient_cutoff=True, gradient_check_transformation=None, return_indices=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L393)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L393?message=Update%20Docs)]
</div>


<a id="Psience.Modes.NormalModes.ReactionPathModes.from_grad_fg" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_grad_fg(cls, basis, gradient, f_matrix, mass_spec, remove_transrot=True, dimensionless=False, zero_freq_cutoff=None, mass_weighted=None, origin=None, projector=None, zero_gradient_cutoff=None, gradient_check_transformation=None, return_status=False, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L669)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L669?message=Update%20Docs)]
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


<a id="Psience.Modes.NormalModes.ReactionPathModes.from_molecule" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_molecule(cls, mol, dimensionless=False, use_internals=None, potential_derivatives=None, project_transrot=True, zero_freq_cutoff=None, masses=None, zero_gradient_cutoff=None, return_status=False, gradient_check_internals=None, gradient_check_transformation=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L735)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L735?message=Update%20Docs)]
</div>


<a id="Psience.Modes.NormalModes.ReactionPathModes.from_modes_and_grad" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_modes_and_grad(cls, modes: Psience.Modes.MixtureModes.MixtureModes, grad: numpy.ndarray, zero_gradient_cutoff=None, use_max_gradient_cutoff=True, return_status=False, mass_weighted=None, **projection_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L818)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L818?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Modes/NormalModes/ReactionPathModes.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Modes/NormalModes/ReactionPathModes.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Modes/NormalModes/ReactionPathModes.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Modes/NormalModes/ReactionPathModes.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/NormalModes.py#L390?message=Update%20Docs)   
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