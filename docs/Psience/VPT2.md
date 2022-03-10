# <a id="Psience.VPT2">Psience.VPT2</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/edit/Psience/VPT2)]
</div>
    
An implementation of vibrational perturbation theory (VPT) that uses sparse matrix methods to obtain
corrections.
Makes heavy use of the `BasisReps` package as well as `McUtils.Coordinerds` to obtain representations
of the corrections to the vibrational Hamiltonian.
Is technically not restricted to VPT in a harmonic basis, but no other forms of PT are likely to
be implemented in the near future.
For the purposes of papers, we've been calling this implementation `PyVibPTn`

The code flow is detailed below

![pt design](/Psience/img/PyVibPTnDesign.png){:width=600px}

<div class="container alert alert-secondary bg-light">
  <div class="row">
   <div class="col" markdown="1">
[VPTRunner](VPT2/Runner/VPTRunner.md)   
</div>
   <div class="col" markdown="1">
[VPTSystem](VPT2/Runner/VPTSystem.md)   
</div>
   <div class="col" markdown="1">
[VPTStateSpace](VPT2/Runner/VPTStateSpace.md)   
</div>
</div>
  <div class="row">
   <div class="col" markdown="1">
[VPTStateMaker](VPT2/Runner/VPTStateMaker.md)   
</div>
   <div class="col" markdown="1">
[VPTHamiltonianOptions](VPT2/Runner/VPTHamiltonianOptions.md)   
</div>
   <div class="col" markdown="1">
[VPTRuntimeOptions](VPT2/Runner/VPTRuntimeOptions.md)   
</div>
</div>
  <div class="row">
   <div class="col" markdown="1">
[VPTSolverOptions](VPT2/Runner/VPTSolverOptions.md)   
</div>
   <div class="col" markdown="1">
[VPTResultsLoader](VPT2/Analyzer/VPTResultsLoader.md)   
</div>
   <div class="col" markdown="1">
[VPTResultsSource](VPT2/Analyzer/VPTResultsSource.md)   
</div>
</div>
  <div class="row">
   <div class="col" markdown="1">
[VPTAnalyzer](VPT2/Analyzer/VPTAnalyzer.md)   
</div>
   <div class="col" markdown="1">
[PerturbationTheoryHamiltonian](VPT2/Hamiltonian/PerturbationTheoryHamiltonian.md)   
</div>
   <div class="col" markdown="1">
[PerturbationTheoryCorrections](VPT2/Corrections/PerturbationTheoryCorrections.md)   
</div>
</div>
  <div class="row">
   <div class="col" markdown="1">
[PerturbationTheorySolver](VPT2/Solver/PerturbationTheorySolver.md)   
</div>
   <div class="col" markdown="1">
[PerturbationTheoryCorrections](VPT2/Corrections/PerturbationTheoryCorrections.md)   
</div>
   <div class="col" markdown="1">
[PerturbationTheoryWavefunctions](VPT2/Wavefunctions/PerturbationTheoryWavefunctions.md)   
</div>
</div>
  <div class="row">
   <div class="col" markdown="1">
[ExpansionTerms](VPT2/Terms/ExpansionTerms.md)   
</div>
   <div class="col" markdown="1">
[KineticTerms](VPT2/Terms/KineticTerms.md)   
</div>
   <div class="col" markdown="1">
[PotentialTerms](VPT2/Terms/PotentialTerms.md)   
</div>
</div>
  <div class="row">
   <div class="col" markdown="1">
[DipoleTerms](VPT2/Terms/DipoleTerms.md)   
</div>
   <div class="col" markdown="1">
[CoriolisTerm](VPT2/Terms/CoriolisTerm.md)   
</div>
   <div class="col" markdown="1">
[PotentialLikeTerm](VPT2/Terms/PotentialLikeTerm.md)   
</div>
</div>
  <div class="row">
   <div class="col" markdown="1">
[PerturbationTheoryStateSpaceFilter](VPT2/StateFilters/PerturbationTheoryStateSpaceFilter.md)   
</div>
</div>
</div>

## Examples

### A note on input formats

A very common use case for this package is to extract more information from a VPT2 calculation
performed by some electronic structure package. To support this use case the [`Molecule`](../Molecools/Molecule.md)
object can read in formatted checkpoint (`.fchk`) files and most of these examples will run off of those.

### Running Jobs

The way to run jobs with the least boilerplate is to use the `run_simple` method of the `VPTRunner` class.

Run second-order vibrational perturbation theory for all states with up to 
three quanta of excitation and print their infrared spectroscopic data

<div class="card in-out-block" markdown="1" id="Markdown_code">

```python
VPTRunner.run_simple(
        "HOH_freq.fchk",
        3,
        logger=True
    )
```

<div class="card-body out-block" markdown="1">

```lang-none
...
                   Harmonic                  Anharmonic
State       Frequency    Intensity       Frequency    Intensity
  0 0 1    3937.52466     67.02051      3744.74223     64.17167
  0 1 0    3803.29960      4.14283      3621.97931      3.11401
  1 0 0    1622.30302     67.45626      1572.70734     68.32367
  0 0 2    7875.04932      0.00000      7391.41648      0.01483
  0 2 0    7606.59919      0.00000      7155.85397      0.31496
  2 0 0    3244.60604      0.00000      3117.39090      0.55473
  0 1 1    7740.82426      0.00000      7200.36337      2.20979
  1 0 1    5559.82768      0.00000      5294.37886      3.76254
  1 1 0    5425.60262      0.00000      5174.61359      0.06232
  0 0 3   11812.57398      0.00000     10940.02275      0.04985
  0 3 0   11409.89879      0.00000     10601.62396      0.00898
  3 0 0    4866.90906      0.00000      4634.05068      0.00350
  0 1 2   11678.34892      0.00000     10680.67944      0.00001
  1 0 2    9497.35234      0.00000      8917.98240      0.00333
  0 2 1   11544.12385      0.00000     10567.87984      0.08362
  2 0 1    7182.13070      0.00000      6815.99171      0.16303
  1 2 0    9228.90221      0.00000      8688.41518      0.00427
  2 1 0    7047.90564      0.00000      6699.22408      0.00661
  1 1 1    9363.12728      0.00000      8729.92693      0.09713
```

</div>
</div>

You can do the same by using the provided objects directly like

<div class="card in-out-block" markdown="1" id="Markdown_code">

```python
system = VPTSystem("HOH_freq.fchk")
states = VPTStateSpace.from_system_and_quanta(system, 3)
pt_opts = VPTSolverOptions(state_space_filters=states.get_filter("intensities"))
run_opts = VPTRuntimeOptions(logger=True)
runner = VPTRunner(system, states, runtime_options=run_opts, solver_options=pt_opts)
runner.print_tables()
``` 

<div class="card-body out-block" markdown="1">

```lang-none
...
                   Harmonic                  Anharmonic
State       Frequency    Intensity       Frequency    Intensity
  0 0 1    3937.52466     67.02051      3744.74223     64.17167
  0 1 0    3803.29960      4.14283      3621.97931      3.11401
  1 0 0    1622.30302     67.45626      1572.70734     68.32367
  0 0 2    7875.04932      0.00000      7391.41648      0.01483
  0 2 0    7606.59919      0.00000      7155.85397      0.31496
  2 0 0    3244.60604      0.00000      3117.39090      0.55473
  0 1 1    7740.82426      0.00000      7200.36337      2.20979
  1 0 1    5559.82768      0.00000      5294.37886      3.76254
  1 1 0    5425.60262      0.00000      5174.61359      0.06232
  0 0 3   11812.57398      0.00000     10940.02275      0.04985
  0 3 0   11409.89879      0.00000     10601.62396      0.00898
  3 0 0    4866.90906      0.00000      4634.05068      0.00350
  0 1 2   11678.34892      0.00000     10680.67944      0.00001
  1 0 2    9497.35234      0.00000      8917.98240      0.00333
  0 2 1   11544.12385      0.00000     10567.87984      0.08362
  2 0 1    7182.13070      0.00000      6815.99171      0.16303
  1 2 0    9228.90221      0.00000      8688.41518      0.00427
  2 1 0    7047.90564      0.00000      6699.22408      0.00661
  1 1 1    9363.12728      0.00000      8729.92693      0.09713
```

</div>
</div>




<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
### <a class="collapse-link" data-toggle="collapse" href="#tests">Tests</a> <a class="float-right" data-toggle="collapse" href="#tests"><i class="fa fa-chevron-down"></i></a>
 </div>
<div class="collapsible-section collapsible-section-body collapse show" id="tests" markdown="1">

- [HOHVPTRunner](#HOHVPTRunner)
- [HOHVPTRunnerFlow](#HOHVPTRunnerFlow)
- [HOHVPTRunnerShifted](#HOHVPTRunnerShifted)
- [HOHVPTRunner3rd](#HOHVPTRunner3rd)
- [GetDegenerateSpaces](#GetDegenerateSpaces)
- [ClHOClRunner](#ClHOClRunner)

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

#### <a name="HOHVPTRunner">HOHVPTRunner</a>
```python
    def test_HOHVPTRunner(self):

        file_name = "HOH_freq.fchk"
        VPTRunner.run_simple(
            TestManager.test_data(file_name),
            3,
            memory_constrained=True,
            logger=True
        )
```
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
#### <a name="HOHVPTRunnerShifted">HOHVPTRunnerShifted</a>
```python
    def test_HOHVPTRunnerShifted(self):

        file_name = "HOH_freq.fchk"
        VPTRunner.run_simple(
            TestManager.test_data(file_name),
            3,
            logger=True,
            corrected_fundamental_frequencies=np.array([1600, 3775, 3880])/UnitsData.convert("Hartrees", "Wavenumbers")
        )
```
#### <a name="HOHVPTRunner3rd">HOHVPTRunner3rd</a>
```python
    def test_HOHVPTRunner3rd(self):
        """
        test that runner works for 3rd order PT, too

        :return:
        :rtype:
        """

        file_name = "HOH_freq.fchk"

        handling_mode="unhandled"

        logger=Logger()
        with logger.block(tag="Internals 2nd Order triad"):
            VPTRunner.run_simple(
                TestManager.test_data(file_name),
                3,
                logger=logger,
                order=2,
                internals=[
                    [0, -1, -1, -1],
                    [1,  0, -1, -1],
                    [2,  0,  1, -1]
                ],
                expansion_order=2,
                degeneracy_specs=[
                    [[0, 0, 1], [2, 0, 0]],
                    [[0, 1, 0], [2, 0, 0]],
                ]
            )

        with logger.block(tag="Internals 3rd Order triad"):
            VPTRunner.run_simple(
                TestManager.test_data(file_name),
                3,
                logger=logger,
                order=3,
                internals=[
                    [0, -1, -1, -1],
                    [1,  0, -1, -1],
                    [2,  0,  1, -1]
                ],
                expansion_order=2,
                degeneracy_specs=[
                    [[0, 0, 1], [2, 0, 0]],
                    [[0, 1, 0], [2, 0, 0]],
                ]
            )

        with logger.block(tag="Internals 2nd Order dyad"):
            VPTRunner.run_simple(
                TestManager.test_data(file_name),
                3,
                logger=logger,
                order=2,
                internals=[
                    [0, -1, -1, -1],
                    [1, 0, -1, -1],
                    [2, 0, 1, -1]
                ],
                expansion_order=2,
                degeneracy_specs=[
                    [[0, 1, 0], [2, 0, 0]],
                ]
            )

        with logger.block(tag="Internals 3rd Order dyad"):
            VPTRunner.run_simple(
                TestManager.test_data(file_name),
                3,
                logger=logger,
                order=3,
                internals=[
                    [0, -1, -1, -1],
                    [1, 0, -1, -1],
                    [2, 0, 1, -1]
                ],
                expansion_order=2,
                degeneracy_specs=[
                    [[0, 1, 0], [2, 0, 0]],
                ]
            )

        with logger.block(tag="Cartesians 2nd Order triad"):
            VPTRunner.run_simple(
                TestManager.test_data(file_name),
                3,
                logger=logger,
                order=2,
                expansion_order=2,
                degeneracy_specs=[
                    [[0, 0, 1], [2, 0, 0]],
                    [[0, 1, 0], [2, 0, 0]],
                ]
            )

        with logger.block(tag="Cartesians 3rd Order triad"):
            VPTRunner.run_simple(
                TestManager.test_data(file_name),
                3,
                logger=logger,
                order=3,
                expansion_order=2,
                mixed_derivative_handling_mode=handling_mode,
                degeneracy_specs=[
                    [[0, 0, 1], [2, 0, 0]],
                    [[0, 1, 0], [2, 0, 0]],
                ]
            )

        with logger.block(tag="Cartesians 2nd Order dyad"):
            VPTRunner.run_simple(
                TestManager.test_data(file_name),
                3,
                logger=logger,
                order=2,
                expansion_order=2,
                degeneracy_specs=[
                    [[0, 1, 0], [2, 0, 0]],
                ]
            )

        with logger.block(tag="Cartesians 3rd Order dyad"):
            VPTRunner.run_simple(
                TestManager.test_data(file_name),
                3,
                logger=logger,
                order=3,
                expansion_order=2,
                degeneracy_specs=[
                    [[0, 1, 0], [2, 0, 0]],
                ]
            )
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
#### <a name="ClHOClRunner">ClHOClRunner</a>
```python
    def test_ClHOClRunner(self):
        file_name = "cl_hocl.fchk"
        state = VPTStateMaker(6)
        COM = -3
        A  = -2
        C  = -1
        _  = 1000
        O  = 0
        H  = 1
        Cl = 2
        X  = 3
        VPTRunner.run_simple(
            TestManager.test_data(file_name),
            [
                state(),
                state([1, 1]),
                state([1, 2]),
                state([1, 3]),
                state([1, 2], [5, 1]),
                state([1, 1], [2, 2]),
            ],
            degenerate_states=[
                [
                    [0, 0, 0, 0, 0, 2],
                    [0, 0, 0, 0, 0, 3],
                    [0, 1, 0, 0, 0, 2],
                    [0, 0, 0, 0, 2, 1]
                ]
            ],
            logger=True,
            handle_strong_couplings=False
            , internals=[
                    [Cl,    _,    _,     _],
                    [ O,   Cl,    _,     _],
                    [ X,    O,   Cl,     _],
                    [ H,    O,   Cl,    X],
                ]
        )
```

 </div>
</div>

___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/edit/ci/examples/Psience/VPT2.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/examples/Psience/VPT2.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/edit/ci/docs/Psience/VPT2.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/edit/?filename=ci/docs/templates/Psience/VPT2.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/edit/Psience/VPT2/__init__.py?message=Update%20Docs)