## <a id="Psience.VPT2.Terms.KineticTerms">KineticTerms</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Terms.py#L1975)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Terms.py#L1975?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Terms/KineticTerms.py#L1985)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Terms/KineticTerms.py#L1985?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Terms.KineticTerms.get_terms" class="docs-object-method">&nbsp;</a> 
```python
get_terms(self, order=None, logger=None, return_expressions=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Terms/KineticTerms.py#L2001)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Terms/KineticTerms.py#L2001?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Terms.KineticTerms.reexpress_G" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
reexpress_G(self, G_expansion, forward_derivs, reverse_derivs=None, order=2): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L2214)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L2214?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Terms/KineticTerms.py#L2255)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Terms/KineticTerms.py#L2255?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L2268)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L2268?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Terms.KineticTerms.optimize_coordinates" class="docs-object-method">&nbsp;</a> 
```python
optimize_coordinates(self, order=2): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Terms/KineticTerms.py#L2307)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Terms/KineticTerms.py#L2307?message=Update%20Docs)]
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Terms/KineticTerms.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Terms/KineticTerms.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Terms/KineticTerms.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Terms/KineticTerms.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Terms.py#L1975?message=Update%20Docs)   
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