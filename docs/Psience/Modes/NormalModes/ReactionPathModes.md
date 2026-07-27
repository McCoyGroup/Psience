## <a id="Psience.Modes.NormalModes.ReactionPathModes">ReactionPathModes</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/NormalModes.py#L278)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/NormalModes.py#L278?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L281)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L281?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute reaction-path-following vibrational modes from a gradient and force-constant/mass specification: mass-weights and normalizes the gradient to form the leading "reaction coordinate" mode, projects the remaining Hessian (and any existing constraint projector) orthogonal to that direction to get the transverse vibrational modes via `get_normal_modes`, and (unless `remove_transrot`) drops the transverse mode most similar to the (near-zero-frequency) projected-out translation/rotation modes so the gradient mode takes its place. Points whose gradient magnitude falls below `zero_gradient_cutoff` are instead treated as ordinary stationary-point normal modes with no reaction coordinate.
  - `gradient`: `np.ndarray`
    > the gradient vector(s) defining the reaction-path direction(s)
  - `f_matrix`: `np.ndarray`
    > the force-constant (Hessian) matrix/matrices
  - `mass_spec`: `np.ndarray | list[str]`
    > the atomic masses (or element symbols), or an already-built G-matrix/mass matrix
  - `remove_transrot`: `bool`
    > whether to remove translation/rotation modes as usual (dropping the one most similar to the reaction-coordinate mode)
  - `dimensionless`: `bool`
    > whether the resulting modes should be frequency-scaled (dimensionless)
  - `mass_weighted`: `bool | None`
    > whether the resulting modes should be mass-weighted
  - `zero_freq_cutoff`: `float | None`
    > frequency cutoff for discarding near-zero-frequency modes
  - `return_gmatrix`: `bool`
    > whether to also return the (broadcast) mass/G-matrix used
  - `projector`: `np.ndarray | None`
    > an existing constraint projector to combine with the gradient-orthogonal projector
  - `zero_gradient_cutoff`: `float | None`
    > the gradient-magnitude cutoff below which a point is treated as having no reaction-path direction; defaults to `cls.zero_gradient_cutoff`
  - `use_max_gradient_cutoff`: `bool`
    > whether to test the cutoff against the maximum gradient component rather than its norm
  - `gradient_check_transformation`: `list[np.ndarray] | None`
    > an optional coordinate-transformation expansion used to re-express the gradient before testing it against the cutoff (without changing the gradient actually used to build the mode)
  - `return_indices`: `bool`
    > whether to also return which batch positions used the reaction-path treatment versus the plain stationary-point treatment
  - `:returns`: `object | tuple`
    > the mode data (`ModeData(freqs, modes, inv)`), or a tuple additionally including the mass/G-matrix and/or the `(regular_mode_pos, rem_pos)` index split depending on the `return_*` flags


<a id="Psience.Modes.NormalModes.ReactionPathModes.from_grad_fg" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_grad_fg(cls, basis, gradient, f_matrix, mass_spec, remove_transrot=True, dimensionless=False, zero_freq_cutoff=None, mass_weighted=None, origin=None, projector=None, zero_gradient_cutoff=None, gradient_check_transformation=None, return_status=False, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L591)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L591?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L657)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L657?message=Update%20Docs)]
</div>
**LLM Docstring**

Build reaction-path-following normal modes for a molecule at a (generally non-stationary) point, using its gradient together with either its internal-coordinate Hessian/G-matrix (if `use_internals`) or its (optionally translation/rotation-projected) Cartesian Hessian/masses, dispatching to `from_grad_fg`.
  - `mol`: `Molecule`
    > the molecule to build reaction-path modes for
  - `dimensionless`: `bool`
    > whether the resulting modes should be dimensionless
  - `use_internals`: `bool | None`
    > whether to build the modes in internal coordinates
  - `potential_derivatives`: `list[np.ndarray] | None`
    > explicit potential derivative tensors (including the gradient) to use instead of the molecule's own/computed ones
  - `project_transrot`: `bool`
    > whether to project out translational/rotational degrees of freedom (Cartesian case only)
  - `zero_freq_cutoff`: `float | None`
    > the frequency cutoff below which modes are discarded as translation/rotation
  - `masses`: `np.ndarray | None`
    > masses to use instead of the molecule's own
  - `zero_gradient_cutoff`: `float | None`
    > the gradient-magnitude cutoff below which the point is treated as stationary
  - `return_status`: `bool`
    > whether to also return whether the reaction-path treatment was actually used
  - `gradient_check_internals`: `object | None`
    > an alternate internal-coordinate specification used to build `gradient_check_transformation`, if that isn't given directly
  - `gradient_check_transformation`: `list[np.ndarray] | None`
    > an optional coordinate-transformation expansion used only to test the gradient against the zero-gradient cutoff
  - `opts`: `dict`
    > extra options forwarded to `from_grad_fg`
  - `:returns`: `ReactionPathModes | tuple`
    > the constructed reaction-path modes, or `(modes, status)` if `return_status` is set


<a id="Psience.Modes.NormalModes.ReactionPathModes.from_modes_and_grad" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_modes_and_grad(cls, modes: Psience.Modes.MixtureModes.MixtureModes, grad: numpy.ndarray, zero_gradient_cutoff=None, use_max_gradient_cutoff=True, return_status=False, mass_weighted=None, **projection_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L772)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L772?message=Update%20Docs)]
</div>
**LLM Docstring**

Given an existing mode set and a gradient vector, build reaction-path-following modes by localizing the modes against the (mass-weighted, normalized) gradient direction as an orthogonal projection constraint, unless the gradient falls below `zero_gradient_cutoff`, in which case the original modes are returned unchanged.
  - `modes`: `MixtureModes`
    > the mode set to re-localize around the gradient direction
  - `grad`: `np.ndarray`
    > the gradient vector
  - `zero_gradient_cutoff`: `float | None`
    > the gradient-magnitude cutoff below which the modes are left unchanged
  - `use_max_gradient_cutoff`: `bool`
    > whether to test the cutoff against the maximum gradient component rather than its norm
  - `return_status`: `bool`
    > whether to also return whether the reaction-path treatment was actually applied
  - `mass_weighted`: `bool | None`
    > whether the resulting modes should be mass-weighted; defaults to `modes.mass_weighted`
  - `projection_opts`: `dict`
    > extra options forwarded to `MixtureModes.localize`
  - `:returns`: `MixtureModes | tuple`
    > the reaction-path modes, or `(modes, status)` if `return_status` is set
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
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/NormalModes.py#L278?message=Update%20Docs)   
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