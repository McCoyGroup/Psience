## <a id="Psience.Modes.LocalizedModes.LocalizedModes">LocalizedModes</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes.py#L13)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes.py#L13?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.Modes.LocalizedModes.LocalizedModes.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, normal_modes: Psience.Modes.MixtureModes.MixtureModes, transformation, inverse=None, origin=None, masses=None, freqs=None, mass_weighted=None, frequency_scaled=None, **etc): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes.py#L15)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes.py#L15?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a `LocalizedModes` object by applying a (possibly non-square) linear transformation to an existing set of normal modes, storing the transformation and the base modes for later use (e.g. by `modify`/`apply_transformation`/`get_complement`).
  - `normal_modes`: `MixtureModes`
    > the base mode set to localize
  - `transformation`: `np.ndarray | tuple`
    > the localizing transformation to apply to the mode-by-coordinates matrix, either a single 2D array (with its transpose used as the inverse unless `inverse` is given) or a `(forward, inverse)` pair
  - `inverse`: `np.ndarray | None`
    > an explicit inverse for `transformation`, overriding the default
  - `origin`: `np.ndarray | None`
    > the reference geometry; defaults to `normal_modes.origin`
  - `masses`: `np.ndarray | None`
    > the atomic masses; defaults to `normal_modes.masses`
  - `freqs`: `np.ndarray | None`
    > explicit frequencies to associate with the localized modes; if `None`, computed lazily on first access
  - `mass_weighted`: `bool | None`
    > accepted for interface consistency but not used directly (mass-weighting is taken from `normal_modes.mass_weighted`)
  - `frequency_scaled`: `bool | None`
    > whether the localized modes are frequency-scaled
  - `etc`: `dict`
    > extra options forwarded to the base `MixtureModes.__init__`
  - `:returns`: `None`
    > None


<a id="Psience.Modes.LocalizedModes.LocalizedModes.__getitem__" class="docs-object-method">&nbsp;</a> 
```python
__getitem__(self, item): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L84)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L84?message=Update%20Docs)]
</div>
Takes a slice of the modes
  - `item`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.Modes.LocalizedModes.LocalizedModes.freqs" class="docs-object-method">&nbsp;</a> 
```python
@property
freqs(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L112)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L112?message=Update%20Docs)]
</div>
**LLM Docstring**

Property getter/setter for the localized modes' frequencies. The getter lazily computes them via `compute_freqs` if not already cached; the setter is currently a no-op (frequencies can't be set directly).
  - `freqs`: `np.ndarray`
    > (setter only) ignored; setting `freqs` directly is not supported
  - `:returns`: `np.ndarray`
    > (getter) the (cached or newly computed) frequencies


<a id="Psience.Modes.LocalizedModes.LocalizedModes.mass_weighted" class="docs-object-method">&nbsp;</a> 
```python
@property
mass_weighted(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L144)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L144?message=Update%20Docs)]
</div>
**LLM Docstring**

Property getter/setter for whether the modes are mass-weighted, delegated to (and required to match) the base modes' own `mass_weighted` state, since a `LocalizedModes` can't independently change mass-weighting.
  - `new`: `bool`
    > (setter only) the requested mass-weighting state; must match `self.base_modes.mass_weighted`
  - `:returns`: `bool`
    > (getter) `self.base_modes.mass_weighted`


<a id="Psience.Modes.LocalizedModes.LocalizedModes.g_matrix" class="docs-object-method">&nbsp;</a> 
```python
@property
g_matrix(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L186)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L186?message=Update%20Docs)]
</div>
**LLM Docstring**

Property getter/setter for the G-matrix, delegated to `self.base_modes.g_matrix`. The setter is currently a no-op.
  - `g`: `np.ndarray`
    > (setter only) ignored
  - `:returns`: `np.ndarray | None`
    > (getter) `self.base_modes.g_matrix`


<a id="Psience.Modes.LocalizedModes.LocalizedModes.modify" class="docs-object-method">&nbsp;</a> 
```python
modify(self, base_modes=None, *, transformation=None, freqs=None, origin=None, masses=None, inverse=None, name=None, mass_weighted=None, frequency_scaled=None, g_matrix=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L222)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L222?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a new `LocalizedModes` with the given fields overridden, defaulting each unspecified field to this object's own value (base modes, localizing transformation/inverse, name), and re-normalizing a `(forward, inverse)`-tuple `transformation` argument the same way the constructor does.
  - `base_modes`: `MixtureModes | None`
    > replacement base mode set; defaults to `self.base_modes`
  - `transformation`: `np.ndarray | tuple | None`
    > replacement localizing transformation; defaults to `self.localizing_transformation[0]`
  - `freqs`: `np.ndarray | None`
    > replacement frequencies
  - `origin`: `np.ndarray | None`
    > replacement reference geometry
  - `masses`: `np.ndarray | None`
    > replacement masses
  - `inverse`: `np.ndarray | None`
    > replacement inverse transformation; defaults to `self.localizing_transformation[1]`
  - `name`: `str | None`
    > replacement name; defaults to `self.name`
  - `mass_weighted`: `bool | None`
    > accepted but not forwarded to the new instance's constructor
  - `frequency_scaled`: `bool | None`
    > replacement frequency-scaling flag
  - `g_matrix`: `np.ndarray | None`
    > replacement G-matrix
  - `:returns`: `LocalizedModes`
    > the new, modified `LocalizedModes`


<a id="Psience.Modes.LocalizedModes.LocalizedModes.make_mass_weighted" class="docs-object-method">&nbsp;</a> 
```python
make_mass_weighted(self, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L303)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L303?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a mass-weighted version of these localized modes by mass-weighting the underlying base modes and re-localizing on top of them.
  - `kwargs`: `dict`
    > extra options forwarded to `self.base_modes.make_mass_weighted`
  - `:returns`: `LocalizedModes`
    > the mass-weighted localized modes


<a id="Psience.Modes.LocalizedModes.LocalizedModes.remove_mass_weighting" class="docs-object-method">&nbsp;</a> 
```python
remove_mass_weighting(self, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L315)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L315?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a non-mass-weighted version of these localized modes by removing the mass-weighting from the underlying base modes and re-localizing on top of them.
  - `kwargs`: `dict`
    > extra options forwarded to `self.base_modes.remove_mass_weighting`
  - `:returns`: `LocalizedModes`
    > the non-mass-weighted localized modes


<a id="Psience.Modes.LocalizedModes.LocalizedModes.make_frequency_scaled" class="docs-object-method">&nbsp;</a> 
```python
make_frequency_scaled(self, freqs=None, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L327)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L327?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a frequency-scaled (dimensionless) version of these localized modes by rescaling the localizing transformation by the per-mode frequency-scaling factors; returns `self` unchanged if already frequency-scaled.
  - `freqs`: `np.ndarray | None`
    > the frequencies to compute the scaling from; defaults to `self.local_freqs`
  - `kwargs`: `dict`
    > accepted but not used in this method's body
  - `:returns`: `LocalizedModes`
    > the frequency-scaled localized modes


<a id="Psience.Modes.LocalizedModes.LocalizedModes.remove_frequency_scaling" class="docs-object-method">&nbsp;</a> 
```python
remove_frequency_scaling(self, freqs=None, **kwargs): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L348)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L348?message=Update%20Docs)]
</div>
**LLM Docstring**

Build a non-frequency-scaled version of these localized modes by undoing the frequency-scaling factor in the localizing transformation; returns `self` unchanged if not already frequency-scaled.
  - `freqs`: `np.ndarray | None`
    > the frequencies to compute the scaling from; defaults to `self.local_freqs`
  - `kwargs`: `dict`
    > accepted but not used in this method's body
  - `:returns`: `LocalizedModes`
    > the non-frequency-scaled localized modes


<a id="Psience.Modes.LocalizedModes.LocalizedModes.compute_hessian" class="docs-object-method">&nbsp;</a> 
```python
compute_hessian(self, system='modes'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L378)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L378?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the Hessian in either the localized-mode basis (via the localizing transformation and the base modes' frequencies) or the Cartesian/coordinate basis (delegating to the base class).
  - `system`: `str`
    > which basis to compute the Hessian in, `'modes'` or `'coords'`
  - `:returns`: `np.ndarray`
    > the Hessian matrix


<a id="Psience.Modes.LocalizedModes.LocalizedModes.apply_transformation" class="docs-object-method">&nbsp;</a> 
```python
apply_transformation(self, transformation, inverse=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L404)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L404?message=Update%20Docs)]
</div>
**LLM Docstring**

Apply an additional linear transformation on top of the existing localizing transformation, returning a new `LocalizedModes` built from the same base modes with the combined transformation.
  - `transformation`: `np.ndarray | tuple`
    > the additional transformation to apply, either a single 2D array (transpose used as inverse unless `inverse` is given) or a `(forward, inverse)` pair
  - `inverse`: `np.ndarray | None`
    > an explicit inverse for `transformation`
  - `opts`: `dict`
    > extra options forwarded to the constructor
  - `:returns`: `LocalizedModes`
    > the new `LocalizedModes` with the combined transformation


<a id="Psience.Modes.LocalizedModes.LocalizedModes.get_complement" class="docs-object-method">&nbsp;</a> 
```python
get_complement(self, concatenate=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L441)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes/LocalizedModes.py#L441?message=Update%20Docs)]
</div>
**LLM Docstring**

Build the mode set complementary to this localized set (i.e. the remaining degrees of freedom not spanned by the localizing transformation), either as a freshly re-diagonalized (locally harmonic) mode set, or -- if `concatenate` -- as a single combined mode set spanning both this localized subspace and its complement.
  - `concatenate`: `bool`
    > whether to return the complementary modes standalone, or concatenated together with this localized set into one combined mode set spanning the full space
  - `:returns`: `MixtureModes`
    > the complementary modes, or (if `concatenate`) the combined full-space mode set
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Modes/LocalizedModes/LocalizedModes.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Modes/LocalizedModes/LocalizedModes.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Modes/LocalizedModes/LocalizedModes.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Modes/LocalizedModes/LocalizedModes.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Modes/LocalizedModes.py#L13?message=Update%20Docs)   
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