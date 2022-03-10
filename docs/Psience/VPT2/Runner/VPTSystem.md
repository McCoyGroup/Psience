## <a id="Psience.VPT2.Runner.VPTSystem">VPTSystem</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Runner.py#L27)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Runner.py#L27?message=Update%20Docs)]
</div>

Provides a little helper for setting up the input
system for a VPT job

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

<a id="Psience.VPT2.Runner.VPTSystem.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, mol, internals=None, dummy_atoms=None, modes=None, mode_selection=None, potential_derivatives=None, potential_function=None, order=2, dipole_derivatives=None, eckart_embed=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Runner.py#L44)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Runner.py#L44?message=Update%20Docs)]
</div>


- `mol`: `str | Molecule`
    >the molecule or system specification to use (doesn't really even need to be a molecule)
- `internals`: `Any`
    >the Z-matrix for the internal coordinates (in the future will support a general function for this too)
- `modes`: `Any`
    >the normal modes to use if not already supplied by the Molecule
- `mode_selection`: `Any`
    >the subset of normal modes to do perturbation theory on
- `potential_derivatives`: `Iterable[np.ndarray]`
    >the derivatives of the potential to use for expansions
- `dipole_derivatives`: `Iterable[np.ndarray]`
    >the set of dipole derivatives to use for expansions

<a id="Psience.VPT2.Runner.VPTSystem.nmodes" class="docs-object-method">&nbsp;</a> 
```python
@property
nmodes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Runner.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Runner.py#L?message=Update%20Docs)]
</div>

Provides the number of modes in the system
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Runner.VPTSystem.get_potential_derivatives" class="docs-object-method">&nbsp;</a> 
```python
get_potential_derivatives(self, potential_function, order=2, **fd_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Runner.py#L113)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Runner.py#L113?message=Update%20Docs)]
</div>

Computes potential derivatives for the given function through finite difference
- `potential_function`: `Any`
    >No description...
- `order`: `Any`
    >No description...
- `fd_opts`: `Any`
    >No description...
- `:returns`: `_`
    >No description...

<a id="Psience.VPT2.Runner.VPTSystem.from_harmonic_scan" class="docs-object-method">&nbsp;</a> 
```python
from_harmonic_scan(scan_array): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2/Runner.py#L134)/[edit](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Runner.py#L134?message=Update%20Docs)]
</div>

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

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Runner/VPTSystem.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Runner/VPTSystem.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Runner/VPTSystem.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Runner/VPTSystem.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/Runner.py#L27?message=Update%20Docs)