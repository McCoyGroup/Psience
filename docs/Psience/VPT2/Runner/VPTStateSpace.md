## <a id="Psience.VPT2.Runner.VPTStateSpace">VPTStateSpace</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L278)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L278?message=Update%20Docs)]
</div>

Provides a helper to make it easier to set up the input
state spaces/degenerate spaces to run the perturbation theory







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.VPT2.Runner.VPTStateSpace.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, states, degeneracy_specs=None, system=None, frequencies=None, evaluator=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L365)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L365?message=Update%20Docs)]
</div>

  - `states`: `list | int`
    > A list of states or a number of quanta to target
  - `degeneracy_specs`: `'auto' | list | dict`
    > A specification of degeneracies, either as polyads, explicit groups of states, or parameters to a method. (see Details for more info)


<a id="Psience.VPT2.Runner.VPTStateSpace.from_system_and_spec" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_system_and_spec(cls, system, spec, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L446)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L446?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.VPTStateSpace.from_system_and_quanta" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_system_and_quanta(cls, system, quanta, target_modes=None, only_target_modes=False, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L455)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L455?message=Update%20Docs)]
</div>
Takes a system and a number of quanta and constructs a state space
based on that
  - `system`: `Any`
    > 
  - `quanta`: `Any`
    > 
  - `opts`: `Any`
    > any of the options supported by
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Runner.VPTStateSpace.get_state_list_from_quanta" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_state_list_from_quanta(cls, n_quanta, n_modes, target_modes=None, only_target_modes=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L484)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L484?message=Update%20Docs)]
</div>
Gets states up to `n_quanta` over `n_modes`
  - `n_quanta`: `int | Iterable[int]`
    > the number of quanta to provide excitations for
  - `n_modes`: `int`
    > the number of modes in the system
  - `target_modes`: `Iterable[int]`
    > modes that must be excited
  - `only_target_modes`: `bool`
    > whether or not to _only_ support excitations in the `target_modes`
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Runner.VPTStateSpace.build_degenerate_state_spaces" class="docs-object-method">&nbsp;</a> 
```python
build_degenerate_state_spaces(self, degeneracy_specs, states, system=None, evaluator=None, freqs=None) -> '(None|DegeneracySpec, None|list[np.ndarray])': 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/VPTStateSpace.py#L514)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/VPTStateSpace.py#L514?message=Update%20Docs)]
</div>

  - `degeneracy_specs`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Runner.VPTStateSpace.filter_generator" class="docs-object-method">&nbsp;</a> 
```python
filter_generator(self, target_property, order=2, initial_states=None, postfilters=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/VPTStateSpace.py#L561)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/VPTStateSpace.py#L561?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.VPTStateSpace.get_filter" class="docs-object-method">&nbsp;</a> 
```python
get_filter(self, target_property, order=2, initial_states=None, postfilters=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner/VPTStateSpace.py#L570)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner/VPTStateSpace.py#L570?message=Update%20Docs)]
</div>
Obtains a state space filter for the given target property
using the states we want to get corrections for
  - `target_property`: `Any`
    > 
  - `order`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Runner.VPTStateSpace.get_state_space_filter" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
get_state_space_filter(cls, states, initial_states=None, n_modes=None, order=2, target='wavefunctions', postfilters=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L588)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L588?message=Update%20Docs)]
</div>
Gets `state_space_filters` for the input `states` targeting some property
  - `states`: `Any`
    > the input states
  - `n_modes`: `int`
    > 
  - `target`: `str`
    > the property to target, one of `('frequencies', 'intensities', 'wavefunctions')`
  - `:returns`: `_`
    >
 </div>
</div>



<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#Details-29e580" markdown="1"> Details</a> <a class="float-right" data-toggle="collapse" href="#Details-29e580"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Details-29e580" markdown="1">
 
There are multiple possible values for the `degeneracy_specs`.
The simplest is to use the automatic approach, in which we supply a numeric type (`int`, `float`, etc.) to use as the `WFC` threshold.
The next simplest is to explicitly supply the groups we want, like

```python
[
    [ # the first resonant space
        state_1,
        state_2,
        state_3
    ],
    [ # the second
        state_5, state_11, ...
    ],
    ...
]
```

We can also supply pairs of relations for determining resonances, like

```python
[
    [state_1,  state_2], # A first relation
    [state_3,  state_4],  # another relation
    ...
]
```

To allow for extra options, you can also supply a `dict`. If you wanted to have a different `wfc_threshold` and you wanted to do the secondary resonant space splitting step with a very large threshold, you could do that by supplying

```python
{
    'wfc_threshold':.1,
    'energy_cutoff':1.0 # in Hartree
}
```

or you can explicitly add extra groups to the pairs of polyad rules by saying

```python
{
    'polyads':[
            [state_1,  state_2], # A first relation
            [state_3,  state_4],  # another relation
            ...
        ],
    'extra_groups': [
        [ # the first resonant space
            state_a,
            state_b,
            state_c
        ],
        [ # the second
            state_d, state_e, ...
        ],
        ...
    ]
}
```

This also allows us to define more resonance handling strategies.

The Martin Test is supported,
```python
{
    'martin_threshold':.1/219465, #in Hartree
}
```

As are total quanta vectors/polyads
```python
{
    'nT': [1, 1, 1, 0, 2, 2, 0] # e.g.
}
```
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Runner/VPTStateSpace.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Runner/VPTStateSpace.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Runner/VPTStateSpace.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Runner/VPTStateSpace.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L278?message=Update%20Docs)   
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