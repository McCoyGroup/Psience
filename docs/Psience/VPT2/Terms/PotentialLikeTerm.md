## <a id="Psience.VPT2.Terms.PotentialLikeTerm">PotentialLikeTerm</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms.py#L3370)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L3370?message=Update%20Docs)]
</div>

This accounts for the potential-like term.
In Cartesian diplacement modes this is the Watson U.
In proper internals, this is the V' term.







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.VPT2.Terms.PotentialLikeTerm.get_terms" class="docs-object-method">&nbsp;</a> 
```python
get_terms(self, order=None, logger=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Terms/PotentialLikeTerm.py#L3377)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms/PotentialLikeTerm.py#L3377?message=Update%20Docs)]
</div>
**LLM Docstring**

Compute the potential-like (Watson `U` / internal-coordinate `V'`) correction term's Taylor-expansion terms: either directly from the trace of the reciprocal-inertia-tensor derivatives (when working purely in Cartesian modes), or, for internal coordinates, via the standard Watson pseudopotential derivation combining the G-matrix and inertia-tensor log-determinant derivatives (`d/dQ[ln(detI) - ln(detG)]`, using a matrix-cookbook identity for the log-determinant derivative).
  - `order`: `int | None`
    > the highest derivative order to compute
  - `logger`: `Logger | None`
    > logger to report progress to; only used when computing the G-matrix terms via the base class
  - `:returns`: `list[np.ndarray]`
    > the potential-like-term expansion, one per order
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Terms/PotentialLikeTerm.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Terms/PotentialLikeTerm.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Terms/PotentialLikeTerm.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Terms/PotentialLikeTerm.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Terms.py#L3370?message=Update%20Docs)   
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