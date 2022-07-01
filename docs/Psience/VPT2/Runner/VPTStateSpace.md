## <a id="Psience.VPT2.Runner.VPTStateSpace">VPTStateSpace</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L142)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L142?message=Update%20Docs)]
</div>

Provides a helper to make it easier to set up the input
state spaces/degenerate spaces to run the perturbation theory

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

<a id="Psience.VPT2.Runner.VPTStateSpace.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, states, degeneracy_specs=None, system=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L151)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L151?message=Update%20Docs)]
</div>


- `states`: `list | int`
    >A list of states or a number of quanta to target
- `degeneracy_specs`: `list | dict`
    >A specification of degeneracies, either as polyads or explicit groups of states

<a id="Psience.VPT2.Runner.VPTStateSpace.from_system_and_quanta" class="docs-object-method">&nbsp;</a> 
```python
from_system_and_quanta(system, quanta, target_modes=None, only_target_modes=False, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L190)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L190?message=Update%20Docs)]
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
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L219)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L219?message=Update%20Docs)]
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
build_degenerate_state_spaces(self, degeneracy_specs, states, system=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L249)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L249?message=Update%20Docs)]
</div>


- `degeneracy_specs`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Runner.VPTStateSpace.filter_generator" class="docs-object-method">&nbsp;</a> 
```python
filter_generator(self, target_property, order=2, postfilters=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L291)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L291?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Runner.VPTStateSpace.get_filter" class="docs-object-method">&nbsp;</a> 
```python
get_filter(self, target_property, order=2, postfilters=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L295)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L295?message=Update%20Docs)]
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
get_state_space_filter(states, n_modes=None, order=2, target='wavefunctions', postfilters=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L312)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L312?message=Update%20Docs)]
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

 </div>
</div>



<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
### <a class="collapse-link" data-toggle="collapse" href="#tests">Tests</a> <a class="float-right" data-toggle="collapse" href="#tests"><i class="fa fa-chevron-down"></i></a>
 </div>
<div class="collapsible-section collapsible-section-body collapse show" id="tests" markdown="1">

- [HOHVPTRunnerFlow](#HOHVPTRunnerFlow)
- [GetDegenerateSpaces](#GetDegenerateSpaces)

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
#### <a class="collapse-link" data-toggle="collapse" href="#test-setup">Setup</a> <a class="float-right" data-toggle="collapse" href="#test-setup"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="test-setup" markdown="1">

Before we can run our examples we should get a bit of setup out of the way.
Since these examples were harvested from the unit tests not all pieces
will be necessary for all situations.
```python
try:
    from Peeves.TestUtils import *
    from Peeves import BlockProfiler
except:
    pass
from unittest import TestCase
from Psience.VPT2 import *
from Psience.Molecools import Molecule
from Psience.BasisReps import HarmonicOscillatorProductBasis, BasisStateSpace
from McUtils.Data import UnitsData
import McUtils.Plots as plt
import McUtils.Numputils as nput
from McUtils.Scaffolding import *
from McUtils.Parallelizers import SerialNonParallelizer, MultiprocessingParallelizer
from McUtils.Zachary import FiniteDifferenceDerivative
import sys, os, numpy as np, itertools as ip
```

All tests are wrapped in a test class
```python
class VPT2Tests(TestCase):
```

 </div>
</div>

#### <a name="HOHVPTRunnerFlow">HOHVPTRunnerFlow</a>
```python
    def test_HOHVPTRunnerFlow(self):

        file_name = "HOH_freq.fchk"
        VPTRunner.run_simple(
            TestManager.test_data(file_name),
            3,
            memory_constrained=True,
            logger=True
        )

        system = VPTSystem(TestManager.test_data(file_name))
        states = VPTStateSpace.from_system_and_quanta(system, 3)
        pt_opts = VPTSolverOptions(state_space_filters=states.get_filter("intensities"))
        run_opts = VPTRuntimeOptions(logger=True)
        runner = VPTRunner(system, states, runtime_options=run_opts, solver_options=pt_opts)
        runner.print_tables()
```
#### <a name="GetDegenerateSpaces">GetDegenerateSpaces</a>
```python
    def test_GetDegenerateSpaces(self):

        base_states = [
            [0, 0, 1],
            [0, 1, 0],
            [0, 2, 1],
            [0, 4, 0]
        ]

        degenerate_states = VPTStateSpace.get_degenerate_polyad_space(
            base_states,
            [
                [
                    [0, 2, 0],
                    [0, 0, 1]
                ]
            ],
        )
```

 </div>
</div>

___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Runner/VPTStateSpace.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Runner/VPTStateSpace.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Runner/VPTStateSpace.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Runner/VPTStateSpace.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L142?message=Update%20Docs)