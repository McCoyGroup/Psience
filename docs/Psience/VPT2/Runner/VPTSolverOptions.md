## <a id="Psience.VPT2.Runner.VPTSolverOptions">VPTSolverOptions</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L643)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L643?message=Update%20Docs)]
</div>

Provides a helper to keep track of the options available
for configuring the way the perturbation theory is applied

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

<a id="Psience.VPT2.Runner.VPTSolverOptions.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, order=2, expansion_order=None, coupled_states=None, total_space=None, flat_total_space=None, state_space_iterations=None, state_space_terms=None, state_space_filters=None, allow_post_PT_calc=None, modify_degenerate_perturbations=None, gaussian_resonance_handling=None, ignore_odd_order_energies=None, intermediate_normalization=None, zero_element_warning=None, degenerate_states=None, handle_strong_couplings=None, strong_coupling_test_modes=None, strong_couplings_state_filter=None, strongly_coupled_group_filter=None, extend_strong_coupling_spaces=None, strong_coupling_zero_order_energy_cutoff=None, low_frequency_mode_cutoff=None, zero_order_energy_corrections=None, check_overlap=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L675)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L675?message=Update%20Docs)]
</div>


- `order`: `int`
    >the order of perturbation theory to apply
- `expansion_order`: `int`
    >the order to go to in the expansions of the perturbations
- `degenerate_states`: `Iterable[BasisStateSpace]`
    >the set of degeneracies to handle
- `coupled_states`: `Iterable[SelectionRuleStateSpace]`
    >explicit bases of states to use at each order in the perturbation theory
- `total_space`: `Iterable[BasisStateSpace]`
    >the total state spaces at each order in the perturbation theory
- `flat_total_space`: `BasisStateSpace`
    >the union of all of the total state spaces
- `state_space_iterations`: `int`
    >the order to go to when getting the `coupled_states`
- `state_space_terms`: `Iterable[(int, int)]`
    >the explicit set of terms to include, as a tuple `(i, j)` which indicates `(H(i), |n(j)>)`
- `state_space_filters`: `dict`
    >filters that can be used to cut down on the size of bases (see `VPTRunner.get_state_space_filter`)
- `allow_post_PT_calc`: `bool`
    >whether to do the post-perturbation theory variational calculation for degeneracy handling
- `modify_degenerate_perturbations`: `bool`
    >whether to modify the perturbation representation matrices themselves when doing degeneracy handling
- `gaussian_resonance_handling`: `bool`
    >whether or not to skip the post-PT variational calculation for states with more than two quanta of excitation
- `ignore_odd_order_energies`: `bool`
    >whether or not to skip actually calculating the energy corrections for odd orders
- `intermediate_normalization`: `bool`
    >whether or not to use 'intermediate normalization' in the wavefunctions
- `zero_element_warning`: `bool`
    >whether or not to warn if an element of the representations evaluated to zero (i.e. we wasted effort)
- `zero_order_energy_corrections`: `dict`
    >energies to use for the zero-order states instead of the diagonal of `H(0)`

<a id="Psience.VPT2.Runner.VPTSolverOptions.get_zero_order_energies" class="docs-object-method">&nbsp;</a> 
```python
get_zero_order_energies(corrected_fundamental_freqs, states): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L769)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L769?message=Update%20Docs)]
</div>


- `corrected_fundamental_freqs`: `Any`
    >No description...
- `states`: `Any`
    >No description...
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

 </div>
</div>

___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Runner/VPTSolverOptions.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Runner/VPTSolverOptions.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Runner/VPTSolverOptions.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Runner/VPTSolverOptions.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L643?message=Update%20Docs)