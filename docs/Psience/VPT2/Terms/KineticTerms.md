## <a id="Psience.VPT2.Terms.KineticTerms">KineticTerms</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L2818)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L2818?message=Update%20Docs)]
</div>

Represents the KE coefficients







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.VPT2.Terms.KineticTerms.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, molecule, g_derivative_threshold=0.001, gmatrix_tolerance=1e-06, use_cartesian_kinetic_energy=False, check_input_gmatrix=True, freq_tolerance=0.002, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L2828)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L2828?message=Update%20Docs)]
</div>
**LLM Docstring**

Set up a G-matrix (kinetic-energy coefficient) expansion generator for a molecule, storing the tolerances/thresholds used to validate and warn about the computed G-matrix and its derivatives.
  - `molecule`: `Molecule`
    > the molecule to compute the G-matrix expansion for
  - `g_derivative_threshold`: `float`
    > the magnitude above which a G-matrix derivative term triggers a logged warning
  - `gmatrix_tolerance`: `float`
    > the tolerance used when checking that the zeroth-order G-matrix is diagonal
  - `use_cartesian_kinetic_energy`: `bool`
    > whether to force computing the kinetic energy directly in Cartesian coordinates rather than internal coordinates
  - `check_input_gmatrix`: `bool`
    > whether to validate the reconstructed G-matrix frequencies against the nominal mode frequencies
  - `freq_tolerance`: `float`
    > the tolerance used for that frequency validation
  - `opts`: `dict`
    > extra options forwarded to the base `ExpansionTerms.__init__`
  - `:returns`: `None`
    > None


<a id="Psience.VPT2.Terms.KineticTerms.get_terms" class="docs-object-method">&nbsp;</a> 
```python
get_terms(self, order=None, logger=None, return_expressions=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/KineticTerms.py#L2866)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/KineticTerms.py#L2866?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the G-matrix (kinetic-energy coefficient) Taylor-expansion terms in the molecule's normal-mode coordinates, either via a simplified direct Cartesian-mode-matrix contraction (when working purely in Cartesians) or by chaining together the internals/Cartesians-by-modes Jacobians via `TensorDerivativeConverter` and differentiating iteratively, validating that the zeroth-order term is diagonal and matches the nominal mode frequencies, and warning about anomalously large higher-order derivatives.
  - `order`: `int | None`
    > the highest derivative order to compute
  - `logger`: `Logger | None`
    > logger to report progress/timing/warnings to; defaults to `self.logger`
  - `return_expressions`: `bool`
    > whether to also return the underlying symbolic `TensorExpansionTerms` objects alongside the numeric arrays
  - `:returns`: `list[np.ndarray] | tuple`
    > the G-matrix expansion terms (or `(terms, expressions)` if `return_expressions` is set)


<a id="Psience.VPT2.Terms.KineticTerms.reexpress_G" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
reexpress_G(self, G_expansion, forward_derivs, reverse_derivs=None, order=2): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L3135)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L3135?message=Update%20Docs)]
</div>
Apply a coordinate transformation to the G-matrix
  - `forward_derivs`: `Any`
    > 
  - `reverse_derivs`: `Any`
    > 
  - `order`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Terms.KineticTerms.reexpress" class="docs-object-method">&nbsp;</a> 
```python
reexpress(self, forward_derivs, reverse_derivs=None, order=2): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/KineticTerms.py#L3176)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/KineticTerms.py#L3176?message=Update%20Docs)]
</div>
Finds a coordinate transformation the give 0 contribution to the G-matrix
  - `forward_derivs`: `Any`
    > 
  - `reverse_derivs`: `Any`
    > 
  - `order`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Terms.KineticTerms.get_kinetic_optimized_coordinates" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_kinetic_optimized_coordinates(cls, G_expansion, order=2): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L3189)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L3189?message=Update%20Docs)]
</div>
**LLM Docstring**

Intended to iteratively find the coordinate transformation that eliminates as much of the G-matrix expansion as possible order by order. As written, after computing the first-order correction `R2` and the resulting re-expressed G-matrix, the method unconditionally executes `raise Exception(new_G[1])` -- it never returns normally, so it does not currently function, and the remaining (otherwise complete-looking) loop below that line is unreachable dead code.
  - `G_expansion`: `list[np.ndarray]`
    > the G-matrix expansion terms to optimize
  - `order`: `int`
    > the target expansion order to optimize up to
  - `:returns`: `None`
    > never returns normally; always raises


<a id="Psience.VPT2.Terms.KineticTerms.optimize_coordinates" class="docs-object-method">&nbsp;</a> 
```python
optimize_coordinates(self, order=2): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/KineticTerms.py#L3241)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/KineticTerms.py#L3241?message=Update%20Docs)]
</div>
**LLM Docstring**

Build the kinetic-energy-optimized coordinate transformation for this G-matrix expansion and re-express it in terms of both the internal coordinates and the Cartesians. Note this method calls `get_kinetic_optimized_coordinates`, which as currently implemented always raises rather than returning a transformation (see that method's docstring), so calling this method will likewise fail.
  - `order`: `int`
    > the highest order to optimize the transformation to
  - `:returns`: `tuple`
    > intended to be `((QR, RQ), (QX, XQ))`, the forward/reverse transformations between the optimized coordinates and both the internal coordinates and the Cartesians; never actually returned due to the exception raised by `get_kinetic_optimized_coordinates`
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Terms/KineticTerms.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Terms/KineticTerms.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Terms/KineticTerms.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Terms/KineticTerms.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L2818?message=Update%20Docs)   
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