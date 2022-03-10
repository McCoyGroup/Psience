## <a id="Psience.VPT2.StateFilters.PerturbationTheoryStateSpaceFilter">PerturbationTheoryStateSpaceFilter</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/StateFilters.py#L9)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/StateFilters.py#L9?message=Update%20Docs)]
</div>

Provides an easier constructor for the VPT state space filters

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

<a id="Psience.VPT2.StateFilters.PerturbationTheoryStateSpaceFilter.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, input_space, prefilters, postfilters): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/StateFilters.py#L13)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/StateFilters.py#L13?message=Update%20Docs)]
</div>


- `input_space`: `BasisStateSpace`
    >No description...
- `prefilters`: `Any`
    >No description...
- `postfilters`: `Any`
    >No description...

<a id="Psience.VPT2.StateFilters.PerturbationTheoryStateSpaceFilter.from_data" class="docs-object-method">&nbsp;</a> 
```python
from_data(input_space, data): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/StateFilters.py#L28)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/StateFilters.py#L28?message=Update%20Docs)]
</div>

Works to canonicalize inputs and initialize appropriately from there
- `data`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.StateFilters.PerturbationTheoryStateSpaceFilter.from_rules" class="docs-object-method">&nbsp;</a> 
```python
from_rules(input_space, *rules): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/StateFilters.py#L69)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/StateFilters.py#L69?message=Update%20Docs)]
</div>

Builds a set of filter spaces from dicts of rules
- `rules`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.StateFilters.PerturbationTheoryStateSpaceFilter.prefilters" class="docs-object-method">&nbsp;</a> 
```python
@property
prefilters(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/StateFilters.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/StateFilters.py#L?message=Update%20Docs)]
</div>


- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.StateFilters.PerturbationTheoryStateSpaceFilter.postfilter" class="docs-object-method">&nbsp;</a> 
```python
@property
postfilter(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/StateFilters.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/StateFilters.py#L?message=Update%20Docs)]
</div>


- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.StateFilters.PerturbationTheoryStateSpaceFilter.canonicalize_prefilters" class="docs-object-method">&nbsp;</a> 
```python
canonicalize_prefilters(self, basis, prefilters): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/StateFilters.py#L180)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/StateFilters.py#L180?message=Update%20Docs)]
</div>

Puts the prefilters in canonical form...
- `basis`: `Any`
    >No description...
- `prefilters`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.StateFilters.PerturbationTheoryStateSpaceFilter.from_property_rules" class="docs-object-method">&nbsp;</a> 
```python
from_property_rules(initial_space, target_space, perturbation_rules, property_rules, order=2): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/StateFilters.py#L227)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/StateFilters.py#L227?message=Update%20Docs)]
</div>


- `initial_space`: `Any`
    >No description...
- `target_space`: `Any`
    >No description...
- `perturbation_rules`: `Any`
    >No description...
- `property_rules`: `Any`
    >No description...
- `order`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.StateFilters.PerturbationTheoryStateSpaceFilter.generate_nquanta_filter" class="docs-object-method">&nbsp;</a> 
```python
generate_nquanta_filter(initials, rules, finals): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/StateFilters.py#L435)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/StateFilters.py#L435?message=Update%20Docs)]
</div>

Takes the initial number of quanta, a set of possible rules, and
        a set of final numbers of quanta and determines which rules apply
- `initial`: `Any`
    >No description...
- `rules`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

 </div>
</div>




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/StateFilters/PerturbationTheoryStateSpaceFilter.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/StateFilters/PerturbationTheoryStateSpaceFilter.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/StateFilters/PerturbationTheoryStateSpaceFilter.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/StateFilters/PerturbationTheoryStateSpaceFilter.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/StateFilters.py#L9?message=Update%20Docs)