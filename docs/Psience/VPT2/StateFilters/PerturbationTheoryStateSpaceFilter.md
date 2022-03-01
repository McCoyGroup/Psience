## <a id="Psience.VPT2.StateFilters.PerturbationTheoryStateSpaceFilter">PerturbationTheoryStateSpaceFilter</a>
Provides an easier constructor for the VPT state space filters

### Properties and Methods
<a id="Psience.VPT2.StateFilters.PerturbationTheoryStateSpaceFilter.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, input_space, prefilters, postfilters): 
```

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
Works to canonicalize inputs and initialize appropriately from there
- `data`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.StateFilters.PerturbationTheoryStateSpaceFilter.from_rules" class="docs-object-method">&nbsp;</a>
```python
from_rules(input_space, *rules): 
```
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

<a id="Psience.VPT2.StateFilters.PerturbationTheoryStateSpaceFilter.postfilter" class="docs-object-method">&nbsp;</a>
```python
@property
postfilter(self): 
```

<a id="Psience.VPT2.StateFilters.PerturbationTheoryStateSpaceFilter.canonicalize_prefilters" class="docs-object-method">&nbsp;</a>
```python
canonicalize_prefilters(self, basis, prefilters): 
```
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

<a id="Psience.VPT2.StateFilters.PerturbationTheoryStateSpaceFilter.generate_nquanta_filter" class="docs-object-method">&nbsp;</a>
```python
generate_nquanta_filter(initials, rules, finals): 
```
Takes the initial number of quanta, a set of possible rules, and
        a set of final numbers of quanta and determines which rules apply
- `initial`: `Any`
    >No description...
- `rules`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

### Examples




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/edit/ci/examples/ci/docs/Psience/VPT2/StateFilters/PerturbationTheoryStateSpaceFilter.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/examples/ci/docs/Psience/VPT2/StateFilters/PerturbationTheoryStateSpaceFilter.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/edit/ci/docs/ci/docs/Psience/VPT2/StateFilters/PerturbationTheoryStateSpaceFilter.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/docs/templates/ci/docs/Psience/VPT2/StateFilters/PerturbationTheoryStateSpaceFilter.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/StateFilters.py?message=Update%20Docs)