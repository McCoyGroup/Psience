## <a id="Psience.VPT2.Terms.PotentialTerms">PotentialTerms</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L1811)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L1811?message=Update%20Docs)]
</div>

A helper class that can transform the derivatives of the potential from Cartesian to normal coordinates







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.VPT2.Terms.PotentialTerms.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, molecule, mixed_derivs=None, modes=None, potential_derivatives=None, mode_selection=None, mode_transformation=None, full_surface_mode_selection=None, logger=None, parallelizer=None, checkpointer=None, check_input_force_constants=True, allow_higher_potential_terms=False, hessian_tolerance=0.0001, grad_tolerance=0.0001, freq_tolerance=0.002, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L1822)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L1822?message=Update%20Docs)]
</div>

  - `molecule`: `Molecule`
    > the molecule that will supply the potential derivatives
  - `mixed_derivs`: `bool`
    > whether or not the pulled derivatives are partially derivatives along the normal coords
  - `modes`: `None | MolecularVibrations`
    > the normal modes to use when doing calculations
  - `mode_selection`: `None | Iterable[int]`
    > the subset of normal modes to use


<a id="Psience.VPT2.Terms.PotentialTerms.v_derivs" class="docs-object-method">&nbsp;</a> 
```python
@property
v_derivs(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/PotentialTerms.py#L1870)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/PotentialTerms.py#L1870?message=Update%20Docs)]
</div>
**LLM Docstring**

Property getter/setter for the canonicalized potential-energy derivative tensors. The getter lazily pulls the raw derivatives from `self.molecule.potential_surface.derivatives` (if none were supplied at construction) and canonicalizes them via `_canonicalize_derivs`, caching the result.
  - `v`: `list[np.ndarray]`
    > (setter only) the new canonicalized derivative tensors to store directly
  - `:returns`: `list[np.ndarray]`
    > (getter) the canonicalized potential derivative tensors


<a id="Psience.VPT2.Terms.PotentialTerms.get_terms" class="docs-object-method">&nbsp;</a> 
```python
get_terms(self, order=None, logger=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/PotentialTerms.py#L2431)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/PotentialTerms.py#L2431?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the potential-energy expansion terms in the molecule's normal-mode coordinates, zeroing out the (should-be-vanishing) gradient term, re-expanding the Cartesian (or internal-coordinate) potential derivatives through the appropriate coordinate Jacobians (handling mixed-derivative-basis terms via `_symmetrize_mixed_derivatives` where relevant), and caching the result via the checkpointer.
  - `order`: `int | None`
    > the highest derivative order to compute; if `None`, uses however many terms are available in `v_derivs`
  - `logger`: `Logger | None`
    > logger to report progress/timing/warnings to; defaults to `self.logger`
  - `:returns`: `list[np.ndarray]`
    > the potential-energy expansion terms, in mode-basis coordinates, from the (zeroed) gradient upward


<a id="Psience.VPT2.Terms.PotentialTerms.get_potential_optimized_coordinates" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_potential_optimized_coordinates(cls, V_expansion, order=2): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L2765)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L2765?message=Update%20Docs)]
</div>
**LLM Docstring**

Find, order by order, the coordinate transformation that eliminates as much of the potential-energy expansion as possible (a "potential-optimized" coordinate system), by solving at each order for the new-coordinate term that cancels the corresponding remainder of the transformed potential.
  - `V_expansion`: `list[np.ndarray]`
    > the potential-energy expansion terms (from the quadratic/Hessian term upward, i.e. without the zeroth-order energy)
  - `order`: `int`
    > the highest order to optimize the transformation to
  - `:returns`: `tuple[list[np.ndarray], list[np.ndarray]]`
    > `(forward_derivs, reverse_derivs)` -- the forward and reverse coordinate-transformation derivative tensors


<a id="Psience.VPT2.Terms.PotentialTerms.optimize_coordinates" class="docs-object-method">&nbsp;</a> 
```python
optimize_coordinates(self, order=2): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/PotentialTerms.py#L2797)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/PotentialTerms.py#L2797?message=Update%20Docs)]
</div>
**LLM Docstring**

Build the potential-optimized coordinate transformation for this potential expansion and re-express it in terms of both the internal coordinates and the Cartesians, via `get_potential_optimized_coordinates` combined with the Cartesians-by-internals/internals-by-Cartesians Jacobians.
  - `order`: `int`
    > the highest order to optimize the transformation to
  - `:returns`: `tuple[tuple[list[np.ndarray], list[np.ndarray]], tuple[list[np.ndarray], list[np.ndarray]]]`
    > `((QR, RQ), (QX, XQ))` -- the forward/reverse transformations between the optimized coordinates and both the internal coordinates (`R`) and the Cartesians (`X`)
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Terms/PotentialTerms.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Terms/PotentialTerms.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Terms/PotentialTerms.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Terms/PotentialTerms.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L1811?message=Update%20Docs)   
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