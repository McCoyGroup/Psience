## <a id="Psience.BasisReps.StateFilters.BasisStateSpaceFilter">BasisStateSpaceFilter</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateFilters.py#L82)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateFilters.py#L82?message=Update%20Docs)]
</div>

Provides an easier constructor for the VPT state space filters







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.BasisReps.StateFilters.BasisStateSpaceFilter.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, input_space, prefilters, postfilters): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateFilters/BasisStateSpaceFilter.py#L87)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateFilters/BasisStateSpaceFilter.py#L87?message=Update%20Docs)]
</div>

  - `input_space`: `BasisStateSpace`
    > 
  - `prefilters`: `Any`
    > 
  - `postfilters`: `Any`
    >


<a id="Psience.BasisReps.StateFilters.BasisStateSpaceFilter.from_data" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_data(cls, input_space, data): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L102)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L102?message=Update%20Docs)]
</div>
Works to canonicalize inputs and initialize appropriately from there
  - `data`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateFilters.BasisStateSpaceFilter.from_rules" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_rules(cls, input_space, *rules): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L143)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L143?message=Update%20Docs)]
</div>
Builds a set of filter spaces from dicts of rules
  - `rules`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateFilters.BasisStateSpaceFilter.prefilters" class="docs-object-method">&nbsp;</a> 
```python
@property
prefilters(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateFilters/BasisStateSpaceFilter.py#L177)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateFilters/BasisStateSpaceFilter.py#L177?message=Update%20Docs)]
</div>

  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateFilters.BasisStateSpaceFilter.postfilters" class="docs-object-method">&nbsp;</a> 
```python
@property
postfilters(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateFilters/BasisStateSpaceFilter.py#L189)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateFilters/BasisStateSpaceFilter.py#L189?message=Update%20Docs)]
</div>

  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateFilters.BasisStateSpaceFilter.canonicalize_postfilters" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
canonicalize_postfilters(self, input_space, filters): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L200)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L200?message=Update%20Docs)]
</div>


<a id="Psience.BasisReps.StateFilters.BasisStateSpaceFilter.canonicalize_prefilters" class="docs-object-method">&nbsp;</a> 
```python
canonicalize_prefilters(self, basis, prefilters): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/BasisReps/StateFilters/BasisStateSpaceFilter.py#L271)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateFilters/BasisStateSpaceFilter.py#L271?message=Update%20Docs)]
</div>
Puts the prefilters in canonical form...
  - `basis`: `Any`
    > 
  - `prefilters`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateFilters.BasisStateSpaceFilter.from_property_rules" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_property_rules(cls, initial_space, target_space, perturbation_rules, property_rules, order=2, postfilters=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L318)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L318?message=Update%20Docs)]
</div>

  - `initial_space`: `Any`
    > 
  - `target_space`: `Any`
    > 
  - `perturbation_rules`: `Any`
    > 
  - `property_rules`: `Any`
    > 
  - `order`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.BasisReps.StateFilters.BasisStateSpaceFilter.generate_nquanta_filter" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
generate_nquanta_filter(cls, initials, rules, finals): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/__init__.py#L665)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/__init__.py#L665?message=Update%20Docs)]
</div>
Takes the initial number of quanta, a set of possible rules, and
a set of final numbers of quanta and determines which rules apply
  - `initial`: `Any`
    > 
  - `rules`: `Any`
    > 
  - `:returns`: `_`
    >
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/BasisReps/StateFilters/BasisStateSpaceFilter.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/BasisReps/StateFilters/BasisStateSpaceFilter.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/BasisReps/StateFilters/BasisStateSpaceFilter.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/BasisReps/StateFilters/BasisStateSpaceFilter.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/BasisReps/StateFilters.py#L82?message=Update%20Docs)   
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