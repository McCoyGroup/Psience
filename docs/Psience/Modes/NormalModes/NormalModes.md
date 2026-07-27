## <a id="Psience.Modes.NormalModes.NormalModes">NormalModes</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/NormalModes.py#L16)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/NormalModes.py#L16?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
name: str
default_projected_zero_freq_cutoff: NoneType
```
<a id="Psience.Modes.NormalModes.NormalModes.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, basis, coeffs, freqs=None, origin=None, masses=None, inverse=None, name=None, mass_weighted=False, frequency_scaled=False, g_matrix=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/NormalModes.py#L19)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/NormalModes.py#L19?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a `NormalModes` object (a `MixtureModes` specialized for harmonic normal modes) from a mode basis, coefficient matrix, and associated frequency/mass/origin data.
  - `basis`: `CoordinateSystem`
    > the coordinate system the modes are expressed in
  - `coeffs`: `np.ndarray`
    > the mode coefficient (coordinates-by-modes) matrix
  - `freqs`: `np.ndarray | None`
    > the vibrational frequencies associated with each mode
  - `origin`: `np.ndarray | None`
    > the reference geometry the modes are expanded about
  - `masses`: `np.ndarray | None`
    > the atomic masses
  - `inverse`: `np.ndarray | None`
    > the modes-by-coordinates (inverse) matrix
  - `name`: `str | None`
    > a name for the mode set
  - `mass_weighted`: `bool`
    > whether the modes are mass-weighted
  - `frequency_scaled`: `bool`
    > whether the modes are frequency-scaled (dimensionless)
  - `g_matrix`: `np.ndarray | None`
    > the G-matrix associated with the modes
  - `:returns`: `None`
    > None


<a id="Psience.Modes.NormalModes.NormalModes.from_fg" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_fg(cls, basis, f_matrix, mass_spec, remove_transrot=True, dimensionless=False, zero_freq_cutoff=None, mass_weighted=None, origin=None, projector=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L72)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L72?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L121)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L121?message=Update%20Docs)]
</div>
**LLM Docstring**

Build normal modes for a molecule, dispatching to `from_fg` with either the internal-coordinate Hessian/G-matrix (if `use_internals`) or the (optionally translation/rotation-projected) Cartesian Hessian/masses, computing missing potential derivatives via the molecule's energy evaluator or falling back to its existing normal modes where possible.
  - `mol`: `Molecule`
    > the molecule to build normal modes for
  - `dimensionless`: `bool`
    > whether the resulting modes should be dimensionless (frequency-scaled)
  - `use_internals`: `bool | None`
    > whether to build the modes in internal coordinates; defaults to whether the molecule has internal coordinates defined
  - `potential_derivatives`: `list[np.ndarray] | None`
    > explicit potential derivative tensors to use instead of the molecule's own/computed ones
  - `project_transrot`: `bool`
    > whether to project out translational/rotational degrees of freedom before diagonalizing (Cartesian case only)
  - `zero_freq_cutoff`: `float | None`
    > the frequency cutoff below which modes are discarded as translation/rotation; defaults to `default_projected_zero_freq_cutoff` or 1 wavenumber when projecting
  - `masses`: `np.ndarray | None`
    > masses to use instead of the molecule's own
  - `energy_evaluator`: `object | None`
    > an explicit energy evaluator to use for computing missing potential derivatives
  - `opts`: `dict`
    > extra options forwarded to `from_fg`
  - `:returns`: `NormalModes`
    > the constructed normal modes


<a id="Psience.Modes.NormalModes.NormalModes.get_reaction_path_modes" class="docs-object-method">&nbsp;</a> 
```python
get_reaction_path_modes(self, grad: numpy.ndarray, zero_gradient_cutoff=None, return_status=False, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/NormalModes/NormalModes.py#L250)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/NormalModes/NormalModes.py#L250?message=Update%20Docs)]
</div>
**LLM Docstring**

Build reaction-path-following normal modes from this mode set and a gradient vector, via `ReactionPathModes.from_modes_and_grad`.
  - `grad`: `np.ndarray`
    > the gradient vector defining the reaction path direction
  - `zero_gradient_cutoff`: `float | None`
    > the gradient-magnitude cutoff below which the gradient direction is treated as zero (falling back to ordinary normal modes)
  - `return_status`: `bool`
    > whether to also return whether the gradient direction was actually used
  - `kwargs`: `dict`
    > extra options forwarded to `ReactionPathModes.from_modes_and_grad`
  - `:returns`: `ReactionPathModes | tuple`
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Modes/NormalModes/NormalModes.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Modes/NormalModes/NormalModes.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Modes/NormalModes/NormalModes.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Modes/NormalModes/NormalModes.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/NormalModes.py#L16?message=Update%20Docs)   
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