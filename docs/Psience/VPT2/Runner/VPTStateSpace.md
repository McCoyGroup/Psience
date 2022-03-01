## <a id="Psience.VPT2.Runner.VPTStateSpace">VPTStateSpace</a>
Provides a helper to make it easier to set up the input
state spaces/degenerate spaces to run the perturbation theory

### Properties and Methods
<a id="Psience.VPT2.Runner.VPTStateSpace.__init__" class="docs-object-method">&nbsp;</a>
```python
__init__(self, states, degeneracy_specs=None): 
```

<a id="Psience.VPT2.Runner.VPTStateSpace.from_system_and_quanta" class="docs-object-method">&nbsp;</a>
```python
from_system_and_quanta(system, quanta, target_modes=None, only_target_modes=False, **opts): 
```
Takes a system and a number of quanta and constructs a state space
        based on that
- `system`: `Any`
    >No description...
- `quanta`: `Any`
    >No description...
- `opts`: `Any`
    >any of the options supported by
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Runner.VPTStateSpace.get_state_list_from_quanta" class="docs-object-method">&nbsp;</a>
```python
get_state_list_from_quanta(n_quanta, n_modes, target_modes=None, only_target_modes=False): 
```
Gets states up to `n_quanta` over `n_modes`
- `n_quanta`: `int | Iterable[int]`
    >the number of quanta to provide excitations for
- `n_modes`: `int`
    >the number of modes in the system
- `target_modes`: `Iterable[int]`
    >modes that must be excited
- `only_target_modes`: `bool`
    >whether or not to _only_ support excitations in the `target_modes`
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Runner.VPTStateSpace.build_degenerate_state_spaces" class="docs-object-method">&nbsp;</a>
```python
build_degenerate_state_spaces(self, degeneracy_specs): 
```

- `degeneracy_specs`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Runner.VPTStateSpace.get_filter" class="docs-object-method">&nbsp;</a>
```python
get_filter(self, target_property, order=2): 
```

<a id="Psience.VPT2.Runner.VPTStateSpace.get_state_space_filter" class="docs-object-method">&nbsp;</a>
```python
get_state_space_filter(states, n_modes=None, order=2, target='wavefunctions'): 
```
Gets `state_space_filters` for the input `states` targeting some property
- `states`: `Any`
    >the input states
- `n_modes`: `int`
    >No description...
- `target`: `str`
    >the property to target, one of `('frequencies', 'intensities', 'wavefunctions')`
- `:returns`: `_`
    >No description...

### Examples




___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/edit/ci/examples/ci/docs/Psience/VPT2/Runner/VPTStateSpace.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/examples/ci/docs/Psience/VPT2/Runner/VPTStateSpace.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/edit/ci/docs/ci/docs/Psience/VPT2/Runner/VPTStateSpace.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/docs/templates/ci/docs/Psience/VPT2/Runner/VPTStateSpace.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Runner.py?message=Update%20Docs)