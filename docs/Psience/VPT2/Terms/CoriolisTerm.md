## <a id="Psience.VPT2.Terms.CoriolisTerm">CoriolisTerm</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L3264)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L3264?message=Update%20Docs)]
</div>

Calculates the Coriolis coupling term







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.VPT2.Terms.CoriolisTerm.get_zetas_and_momi" class="docs-object-method">&nbsp;</a> 
```python
get_zetas_and_momi(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/CoriolisTerm.py#L3268)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/CoriolisTerm.py#L3268?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the Coriolis zeta constants and the inertial-frame moments of inertia: mass-weights and frequency-descales the mode matrix, reshapes it into an atom/mode/Cartesian-axis tensor, rotates it into the inertial (principal-axis) frame, and forms the antisymmetric (Levi-Civita) combination across atoms that gives the zeta constants.
  - `:returns`: `tuple[np.ndarray, np.ndarray]`
    > `(zeta, B_e)` -- the zeta-constant tensor (mode x mode x 3 x 3) and the inertial-frame moments of inertia


<a id="Psience.VPT2.Terms.CoriolisTerm.get_zetas" class="docs-object-method">&nbsp;</a> 
```python
get_zetas(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/CoriolisTerm.py#L3308)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/CoriolisTerm.py#L3308?message=Update%20Docs)]
</div>
**LLM Docstring**

The Coriolis zeta-constant tensor alone, via `get_zetas_and_momi`.
  - `:returns`: `np.ndarray`
    > the zeta-constant tensor


<a id="Psience.VPT2.Terms.CoriolisTerm.get_terms" class="docs-object-method">&nbsp;</a> 
```python
get_terms(self, order=None, J=0): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/CoriolisTerm.py#L3322)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/CoriolisTerm.py#L3322?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the Coriolis rotational-coupling operator's Taylor-expansion terms by combining the frequency-dimensioned zeta constants with the expansion of the reciprocal moment-of-inertia tensor, adding one extra coordinate axis per order.
  - `order`: `int | None`
    > the highest derivative order to compute
  - `J`: `int`
    > the total rotational angular momentum quantum number; only `J=0` (pure vibration-rotation coupling with no external rotation) is currently supported
  - `:returns`: `list[np.ndarray]`
    > the Coriolis expansion terms, one per order
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Terms/CoriolisTerm.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Terms/CoriolisTerm.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Terms/CoriolisTerm.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Terms/CoriolisTerm.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L3264?message=Update%20Docs)   
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