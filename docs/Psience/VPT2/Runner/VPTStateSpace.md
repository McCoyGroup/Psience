## <a id="Psience.VPT2.Runner.VPTStateSpace">VPTStateSpace</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Runner.py#L138)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Runner.py#L138?message=Update%20Docs)]
</div>

Provides a helper to make it easier to set up the input
state spaces/degenerate spaces to run the perturbation theory

<a id="Psience.VPT2.Runner.VPTStateSpace.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, states, degeneracy_specs=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Runner.py#L147)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Runner.py#L147?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Runner.VPTStateSpace.from_system_and_quanta" class="docs-object-method">&nbsp;</a> 
```python
from_system_and_quanta(system, quanta, target_modes=None, only_target_modes=False, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Runner.py#L162)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Runner.py#L162?message=Update%20Docs)]
</div>

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
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Runner.py#L190)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Runner.py#L190?message=Update%20Docs)]
</div>

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
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Runner.py#L220)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Runner.py#L220?message=Update%20Docs)]
</div>


- `degeneracy_specs`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Runner.VPTStateSpace.get_filter" class="docs-object-method">&nbsp;</a> 
```python
get_filter(self, target_property, order=2): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Runner.py#L257)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Runner.py#L257?message=Update%20Docs)]
</div>

Obtains a state space filter for the given target property
        using the states we want to get corrections for
- `target_property`: `Any`
    >No description...
- `order`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Runner.VPTStateSpace.get_state_space_filter" class="docs-object-method">&nbsp;</a> 
```python
get_state_space_filter(states, n_modes=None, order=2, target='wavefunctions'): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Runner.py#L273)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Runner.py#L273?message=Update%20Docs)]
</div>

Gets `state_space_filters` for the input `states` targeting some property
- `states`: `Any`
    >the input states
- `n_modes`: `int`
    >No description...
- `target`: `str`
    >the property to target, one of `('frequencies', 'intensities', 'wavefunctions')`
- `:returns`: `_`
    >No description...



___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/ci/docs/Psience/VPT2/Runner/VPTStateSpace.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/ci/docs/Psience/VPT2/Runner/VPTStateSpace.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/ci/docs/Psience/VPT2/Runner/VPTStateSpace.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/ci/docs/Psience/VPT2/Runner/VPTStateSpace.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Runner.py#L138?message=Update%20Docs)