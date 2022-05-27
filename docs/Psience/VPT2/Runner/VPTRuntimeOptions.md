## <a id="Psience.VPT2.Runner.VPTRuntimeOptions">VPTRuntimeOptions</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L550)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L550?message=Update%20Docs)]
</div>

Provides a helper to keep track of the options available
for configuring the way the code runs

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

<a id="Psience.VPT2.Runner.VPTRuntimeOptions.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, operator_chunk_size=None, logger=None, verbose=None, checkpoint=None, results=None, parallelizer=None, memory_constrained=None, checkpoint_keys=None, use_cached_representations=None, use_cached_basis=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L567)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L567?message=Update%20Docs)]
</div>


- `operator_chunk_size`: `int`
    >the number of representation matrix elements to calculate at once
- `logger`: `Logger`
    >the `Logger` object to use when logging the status of the calculation
- `verbose`: `bool`
    >whether or not to be verbose in log output
- `checkpoint`: `str`
    >the checkpoint file or `Checkpointer` object to use
- `parallelizer`: `Parallelizer`
    >the `Parallelizer` object to use when parallelizing pieces of the calculation
- `memory_constrained`: `bool`
    >whether or not to attempt memory optimizations
- `checkpoint_keys`: `Iterable[str]`
    >the keys to write to the checkpoint file
- `use_cached_representations`: `bool`
    >whether or not to try to load representation matrices from the checkpoint
- `use_cached_basis`: `bool`
    >whether or not to try to load the bases to use from the checkpoint

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

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Runner/VPTRuntimeOptions.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Runner/VPTRuntimeOptions.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Runner/VPTRuntimeOptions.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Runner/VPTRuntimeOptions.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L550?message=Update%20Docs)