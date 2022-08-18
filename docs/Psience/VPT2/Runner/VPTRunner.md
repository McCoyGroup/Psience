## <a id="Psience.VPT2.Runner.VPTRunner">VPTRunner</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L942)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L942?message=Update%20Docs)]
</div>

A helper class to make it easier to run jobs by making the inputs/options
clear and making it easier to customize run options



<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
 
### <a class="collapse-link" data-toggle="collapse" href="#methods">Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>

 </div>
 <div class="collapsible-section collapsible-section-body collapse" id="methods" markdown="1">

```python
helpers: type
```
<a id="Psience.VPT2.Runner.VPTRunner.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, system, states, hamiltonian_options=None, solver_options=None, runtime_options=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L948)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L948?message=Update%20Docs)]
</div>


- `runtime_options`: `VPTRuntimeOptions`
    >options to configure the way the code runs
- `solver_options`: `VPTSolverOptions`
    >options to configure the way the perturbation theory is applied
- `hamiltonian_options`: `VPTHamiltonianOptions`
    >options to configure the Hamiltonian
- `system`: `VPTSystem`
    >the system to run perturbation theory on

<a id="Psience.VPT2.Runner.VPTRunner.get_Hamiltonian" class="docs-object-method">&nbsp;</a> 
```python
get_Hamiltonian(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L989)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L989?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Runner.VPTRunner.hamiltonian" class="docs-object-method">&nbsp;</a> 
```python
@property
hamiltonian(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Runner.VPTRunner.get_wavefunctions" class="docs-object-method">&nbsp;</a> 
```python
get_wavefunctions(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L1002)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L1002?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Runner.VPTRunner.print_output_tables" class="docs-object-method">&nbsp;</a> 
```python
print_output_tables(wfns=None, file=None, print_intensities=True, print_energies=True, print_energy_corrections=True, print_transition_moments=True, logger=None, sep_char='=', sep_len=100): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L1012)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L1012?message=Update%20Docs)]
</div>

Prints a bunch of formatted output data from a PT run
- `:returns`: `_`
    >
- `wfns`: `Any`
    >

<a id="Psience.VPT2.Runner.VPTRunner.print_tables" class="docs-object-method">&nbsp;</a> 
```python
print_tables(self, wfns=None, file=None, print_intensities=True, print_energy_corrections=True, print_transition_moments=True, sep_char='=', sep_len=100): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L1083)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L1083?message=Update%20Docs)]
</div>

Prints a bunch of formatted output data from a PT run
- `:returns`: `_`
    >
- `wfns`: `Any`
    >

<a id="Psience.VPT2.Runner.VPTRunner.construct" class="docs-object-method">&nbsp;</a> 
```python
construct(system, states, target_property=None, basis_filters=None, corrected_fundamental_frequencies=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L1109)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L1109?message=Update%20Docs)]
</div>

<a id="Psience.VPT2.Runner.VPTRunner.run_simple" class="docs-object-method">&nbsp;</a> 
```python
run_simple(system, states, target_property=None, corrected_fundamental_frequencies=None, calculate_intensities=True, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/VPT2/Runner.py#L1183)/[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L1183?message=Update%20Docs)]
</div>

The standard runner for VPT.
Makes a runner using the `construct` method and then calls that
runner's `print_tables` method after printing out run info.
- `opts`: `Any`
    >options that work for a `VPTSystem`, `VPTStateSpace`, `VPTRuntimeOptions`, `VPTSolverOptions`, or `VPTHamiltonianOptions` object which will be filtered automatically
- `calculate_intensities`: `bool default:True`
    >whether or not to calculate energies
- `corrected_fundamental_frequencies`: `Iterable[float]|None`
    >a set of fundamental frequencies to use to get new zero-order energies
- `target_property`: `str`
    >the target property to get corrections for (one of 'frequencies', 'intensities', 'wavefunctions')
- `states`: `int|list`
    >the states to get corrections for either an `int` (up to that many quanta) or an explicit state list
- `system`: `list|str|Molecule`
    >the system spec, either as a `Molecule`, molecule spec (atoms, coords, opts) or a file to construct a `Molecule`

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
- [ClHOClRunner](#ClHOClRunner)
- [AnalyticModels](#AnalyticModels)
- [HOHCorrectedDegeneracies](#HOHCorrectedDegeneracies)
- [WaterSkippedCouplings](#WaterSkippedCouplings)
- [H2COPolyads](#H2COPolyads)
- [H2COModeSel](#H2COModeSel)
- [HODRephase](#HODRephase)
- [HOHRephase](#HOHRephase)
- [NH3](#NH3)
- [HOONO](#HOONO)
- [H2COSkippedCouplings](#H2COSkippedCouplings)
- [WaterDimerSkippedCouplings](#WaterDimerSkippedCouplings)

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
    """ Threshold = 0
    
    ::> building ExpansionRepresentation<H(0)>
      ::> in Representation<T(0)>
        > evaluating in BraKet space BraKetSpace(nstates=210)
        > evaluating 210 elements over 21 unique indices sequentially
        > took 0.127s
      <::
      ::> in Representation<V(0)>
        > evaluating in BraKet space BraKetSpace(nstates=210)
        > evaluating 210 elements over 21 unique indices sequentially
        > took 0.184s
      <::
      > took 0.520s
    <::
    ::> building ExpansionRepresentation<H(1)>
      ::> in Representation<T(1)>
        > evaluating in BraKet space BraKetSpace(nstates=1190)
        > took 0.000s
      <::
      ::> in Representation<V(1)>
        > evaluating in BraKet space BraKetSpace(nstates=1190)
        > evaluating 1190 elements over 56 unique indices sequentially
        > took 0.204s
      <::
      > took 0.499s
    <::
    ::> building ExpansionRepresentation<H(2)>
      ::> in Representation<T(2)>
        > evaluating in BraKet space BraKetSpace(nstates=43)
        > took 0.000s
      <::
      ::> in Representation<V(2)>
        > evaluating in BraKet space BraKetSpace(nstates=43)
        > evaluating 43 elements over 126 unique indices sequentially
        > took 0.657s
      <::
      ::> in Representation<Coriolis(0)>
        > evaluating in BraKet space BraKetSpace(nstates=43)
        > evaluating 43 elements over 906 unique indices sequentially
        > took 1.005s
      <::
      ::> in Representation<V'(0)>
        > evaluating in BraKet space BraKetSpace(nstates=43)
        > evaluating identity tensor over 43 elements
        > took 0.036s
      <::
      > took 2.244s
    <::
  <::
  
  ::> Energy Corrections
  > State          <0|dH(2)|0>  <0|dH(1)|1> 
  0 0 0 0 0 0     -5.49879    -75.41416
  0 0 0 0 0 1     48.22641   -294.08843
  0 0 0 0 1 0     46.72377   -284.71337
  0 0 0 1 0 0      2.02819   -114.86847
  0 0 1 0 0 0    -32.88017    -81.93283
  0 1 0 0 0 0    -34.90506    -66.80892
  1 0 0 0 0 0    -46.81545    -55.50575
  0 1 0 1 0 0    -25.19818   -114.53080

  0 0 0 0 0 1    3061.70147     95.24194      2849.45769     63.44723
  0 0 0 0 1 0    2977.64050     69.29750      2820.56385     64.99539
  0 0 0 1 0 0    1727.08265     63.79277      1695.15532     65.08450
  0 0 1 0 0 0    1527.04079     11.12160      1493.14075      9.28519
  0 1 0 0 0 0    1252.16397      9.69252      1231.36294     10.18280
  1 0 0 0 0 0    1188.11375      7.01998      1166.70551      7.08986
  0 1 0 1 0 0    2979.24662      0.00000      2967.72530     43.44534
  """
    """  Threshold = 0.05 cm^-1
      0 0 0 0 0 1    3061.70147     95.24194      2849.45684     63.44730
      0 0 0 0 1 0    2977.64050     69.29750      2820.56243     64.99418
      0 0 0 1 0 0    1727.08265     63.79277      1695.15532     65.08644
      0 0 1 0 0 0    1527.04080     11.12160      1493.14075      9.28423
      0 1 0 0 0 0    1252.16397      9.69252      1231.36294     10.18282
      1 0 0 0 0 0    1188.11376      7.01998      1166.70551      7.08986
      0 1 0 1 0 0    2979.24662      0.00000      2967.72474     43.44434

      ::> building ExpansionRepresentation<H(0)>
          ::> in Representation<T(0)>
            > evaluating in BraKet space BraKetSpace(nstates=144)
            > evaluating 144 elements over 21 unique indices sequentially
            > took 0.089s
          <::
          ::> in Representation<V(0)>
            > evaluating in BraKet space BraKetSpace(nstates=144)
            > evaluating 144 elements over 21 unique indices sequentially
            > took 0.176s
          <::
          > took 0.449s
        <::
        ::> building ExpansionRepresentation<H(1)>
          ::> in Representation<T(1)>
            > evaluating in BraKet space BraKetSpace(nstates=287)
            > took 0.000s
          <::
          ::> in Representation<V(1)>
            > evaluating in BraKet space BraKetSpace(nstates=287)
            > evaluating 287 elements over 56 unique indices sequentially
            > took 0.238s
          <::
          > took 0.559s
        <::
        ::> building ExpansionRepresentation<H(2)>
          ::> in Representation<T(2)>
            > evaluating in BraKet space BraKetSpace(nstates=21)
            > took 0.000s
          <::
          ::> in Representation<V(2)>
            > evaluating in BraKet space BraKetSpace(nstates=21)
            > evaluating 21 elements over 126 unique indices sequentially
            > took 0.415s
          <::
          ::> in Representation<Coriolis(0)>
            > evaluating in BraKet space BraKetSpace(nstates=21)
            > evaluating 21 elements over 906 unique indices sequentially
            > took 0.506s
          <::
          ::> in Representation<V'(0)>
            > evaluating in BraKet space BraKetSpace(nstates=21)
            > evaluating identity tensor over 21 elements
            > took 0.118s
          <::
          > took 1.760s
        <::
      """
    """ Threshold = 1.0 cm^-1
    
    ::> building ExpansionRepresentation<H(0)>
      ::> in Representation<T(0)>
        > evaluating in BraKet space BraKetSpace(nstates=144)
        > evaluating 144 elements over 21 unique indices sequentially
        > took 0.063s
      <::
      ::> in Representation<V(0)>
        > evaluating in BraKet space BraKetSpace(nstates=144)
        > evaluating 144 elements over 21 unique indices sequentially
        > took 0.142s
      <::
      > took 0.582s
    <::
    ::> building ExpansionRepresentation<H(1)>
      ::> in Representation<T(1)>
        > evaluating in BraKet space BraKetSpace(nstates=287)
        > took 0.000s
      <::
      ::> in Representation<V(1)>
        > evaluating in BraKet space BraKetSpace(nstates=287)
        > evaluating 287 elements over 56 unique indices sequentially
        > took 0.262s
      <::
      > took 0.901s
    <::
    ::> building ExpansionRepresentation<H(2)>
      ::> in Representation<T(2)>
        > evaluating in BraKet space BraKetSpace(nstates=19)
        > took 0.000s
      <::
      ::> in Representation<V(2)>
        > evaluating in BraKet space BraKetSpace(nstates=19)
        > evaluating 19 elements over 126 unique indices sequentially
        > took 0.336s
      <::
      ::> in Representation<Coriolis(0)>
        > evaluating in BraKet space BraKetSpace(nstates=19)
        > evaluating 19 elements over 906 unique indices sequentially
        > took 0.601s
      <::
      ::> in Representation<V'(0)>
        > evaluating in BraKet space BraKetSpace(nstates=19)
        > evaluating identity tensor over 19 elements
        > took 0.064s
      <::
      > took 1.756s
    <::
  <::
  
  0 0 0 0 0 0     -4.96621    -75.41416
  0 0 0 0 0 1     48.17888   -294.08843
  0 0 0 0 1 0     46.58555   -284.71337
  0 0 0 1 0 0      1.52477   -114.86847
  0 0 1 0 0 0    -33.06100    -81.93283
  0 1 0 0 0 0    -34.75406    -66.80892
  1 0 0 0 0 0    -47.74137    -55.50575
  0 1 0 1 0 0    -26.31829   -114.53080

  0 0 0 0 0 1    3061.70147     95.24194      2848.44632     62.90510
  0 0 0 0 1 0    2977.64050     69.29750      2819.89305     64.85348
  0 0 0 1 0 0    1727.08265     63.79277      1694.11932     65.38942
  0 0 1 0 0 0    1527.04080     11.12160      1492.42734      9.04394
  0 1 0 0 0 0    1252.16397      9.69252      1230.98136     10.06742
  1 0 0 0 0 0    1188.11376      7.01998      1165.24700      7.08479
  0 1 0 1 0 0    2979.24662      0.00000      2966.50387     43.86153
    """
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
            logger=True,
            handle_strong_couplings=True,
            strong_coupling_test_modes=list(range(3, 6))
        )
        """
        > [[0 0 0 0 0 2]
        >  [1 0 0 0 0 2]
        >  [0 1 0 0 0 2]
        >  [0 0 0 0 2 1]
        >  [0 0 0 0 4 0]
        >  [0 2 0 0 0 2]
        >  [0 1 0 0 2 1]]
                                 Harmonic                  Anharmonic
        State             Frequency    Intensity       Frequency    Intensity
          0 0 0 0 0 1    2709.16096   2782.25433      2285.38768   2611.96281
          0 0 0 0 0 2    5418.32192      0.00000      4140.45935     19.13726
          0 0 0 0 0 3    8127.48288      0.00000      5353.97008      0.16004
          0 1 0 0 0 2    5699.88024      0.00000      4605.01290      3.93389
          0 0 0 0 2 1    5592.48466      0.00000      5023.09956      7.99053
        """
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
            logger=True,
            handle_strong_couplings=True
        )
        """        
        > [[0 0 0 0 0 2]
        >  [0 0 0 0 2 1]
        >  [0 0 0 0 4 0]]
                                 Harmonic                  Anharmonic
        State             Frequency    Intensity       Frequency    Intensity
          0 0 0 0 0 1    2709.16096   2782.25433      2285.38768   2611.96281
          0 0 0 0 0 2    5418.32192      0.00000      4096.10015     21.49976
          0 0 0 0 0 3    8127.48288      0.00000      5353.97008      0.16004
          0 1 0 0 0 2    5699.88024      0.00000      4374.63719      2.66632
          0 0 0 0 2 1    5592.48466      0.00000      4964.16010      6.86648
        """
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
            logger=True,
            handle_strong_couplings=True
            , strong_coupling_test_modes=list(range(3, 6))
            , internals=[
                [Cl, _, _, _],
                [O, Cl, _, _],
                [X, O, Cl, _],
                [H, O, Cl, X],
            ]
        )
        """        
        > [[0 0 0 0 0 2]
        >  [1 0 0 0 0 2]
        >  [0 1 0 0 0 2]
        >  [0 0 0 0 2 1]
        >  [0 0 0 0 4 0]
        >  [0 2 0 0 0 2]
        >  [0 1 0 0 2 1]]
                                 Harmonic                  Anharmonic
        State             Frequency    Intensity       Frequency    Intensity
          0 0 0 0 0 1    2709.16096   2782.25433      2305.91708   2613.72238
          0 0 0 0 0 2    5418.32192      0.00000      4174.53466     17.82679
          0 0 0 0 0 3    8127.48288      0.00000      5353.97246      0.16004
          0 1 0 0 0 2    5699.88024      0.00000      4645.21985      6.56429
          0 0 0 0 2 1    5592.48466      0.00000      5114.27886      9.19226
        """
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
            logger=True,
            handle_strong_couplings=True
            , internals=[
                [Cl, _, _, _],
                [O, Cl, _, _],
                [X, O, Cl, _],
                [H, O, Cl, X],
            ]
        )

        """
        > [[0 0 0 0 0 2]
        >  [0 0 0 0 2 1]
        >  [0 0 0 0 4 0]]
                                 Harmonic                  Anharmonic
        State             Frequency    Intensity       Frequency    Intensity
          0 0 0 0 0 1    2709.16096   2782.25433      2305.91708   2613.72238
          0 0 0 0 0 2    5418.32192      0.00000      4130.15930     22.20869
          0 0 0 0 0 3    8127.48288      0.00000      5353.97246      0.16004
          0 1 0 0 0 2    5699.88024      0.00000      4374.63638      2.66634
          0 0 0 0 2 1    5592.48466      0.00000      5053.08138      7.79496
        """
```
#### <a name="AnalyticModels">AnalyticModels</a>
```python
    def test_AnalyticModels(self):
        from Psience.AnalyticModels import AnalyticModel as Model
        from McUtils.Data import AtomData, UnitsData

        order = 4
        # include_potential=False
        # include_gmatrix=True
        # include_pseudopotential=True
        expansion_order={
            'potential':0,
            'gmatrix':4,
            'pseudopotential':4
        }

        hoh_params = {}

        hoh_params["mH"] = AtomData["H"]["Mass"] * UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
        hoh_params["mO"] = AtomData["O"]["Mass"] * UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")

        # morse stretch parameters
        cm2borh = UnitsData.convert("Angstroms", "BohrRadius")
        hoh_params['re'] = 0.9575 * cm2borh
        erg2h = UnitsData.convert("Ergs", "Hartrees")
        invcm2borh = UnitsData.convert("InverseAngstroms", "InverseBohrRadius")
        hoh_params['De'] = 8.84e-12 * erg2h
        hoh_params['a'] = 2.175 * invcm2borh

        # harmonic bend parameters
        hoh_params['b_e'] = np.deg2rad(104.5)
        hoh_params['k_b'] = 3.2 ** 2 * 1600 * UnitsData.convert("Wavenumbers",
                                                                "Hartrees")  # the 3.2 is some deep magic I don't understand


        model = Model(
            [
                Model.r(1, 2),
                Model.r(2, 3),
                Model.a(1, 2, 3)
            ],
            Model.Potential.morse(1, 2,
                                  De=hoh_params["De"],
                                  a=hoh_params["a"],
                                  re=hoh_params["re"]
                                  )
            + Model.Potential.morse(2, 3,
                                    De=hoh_params["De"],
                                    a=hoh_params["a"],
                                    re=hoh_params["re"]
                                    )
            + Model.Potential.harmonic(1, 2, 3,
                                       k=hoh_params["k_b"],
                                       qe=hoh_params["b_e"]
                                       ),
            # dipole=Model.r(1, 2) + Model.r(2, 3),
            values={
                Model.m(1): hoh_params["mH"],
                Model.m(2): hoh_params["mO"],
                Model.m(3): hoh_params["mH"],
                Model.r(1, 2): hoh_params["re"],
                Model.r(2, 3): hoh_params["re"],
                Model.a(1, 2, 3): hoh_params["b_e"]
            }
        )

        model.run_VPT(order=order, return_analyzer=False, expansion_order=expansion_order)

        class harmonically_coupled_morse:
            # mass_weights = masses[:2] / np.sum(masses[:2])
            def __init__(self,
                         De_1, a_1, re_1,
                         De_2, a_2, re_2,
                         kb, b_e
                         ):
                self.De_1 = De_1
                self.a_1 = a_1
                self.re_1 = re_1
                self.De_2 = De_2
                self.a_2 = a_2
                self.re_2 = re_2
                self.kb = kb
                self.b_e = b_e

            def __call__(self, carts):
                v1 = carts[..., 1, :] - carts[..., 0, :]
                v2 = carts[..., 2, :] - carts[..., 0, :]
                r1 = nput.vec_norms(v1) - self.re_1
                r2 = nput.vec_norms(v2) - self.re_2
                bend, _ = nput.vec_angles(v1, v2)
                bend = bend - self.b_e

                return (
                        self.De_1 * (1 - np.exp(-self.a_1 * r1)) ** 2
                        + self.De_2 * (1 - np.exp(-self.a_2 * r2)) ** 2
                        + self.kb * bend ** 2
                )

        atoms = ["O", "H", "H"]
        coords = np.array([
            [0.000000, 0.000000, 0.000000],
            [hoh_params["re"], 0.000000, 0.000000],
            np.dot(
                nput.rotation_matrix([0, 0, 1], hoh_params["b_e"]),
                [hoh_params["re"], 0.000000, 0.000000]
            )
        ])
        masses = np.array([AtomData[x]["Mass"] for x in atoms]) * UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
        pot_file = os.path.expanduser('~/Desktop/water_pot.hdf5')
        water_chk = Checkpointer.from_file(pot_file)
        if expansion_order['potential'] > -1:
            with water_chk as wat:
                try:
                    potential_derivatives = wat['potential_derivatives']
                except (OSError, KeyError):
                    potential_function = harmonically_coupled_morse(
                        hoh_params["De"], hoh_params["a"], hoh_params["re"],
                        hoh_params["De"], hoh_params["a"], hoh_params["re"],
                        hoh_params["k_b"], hoh_params["b_e"]
                    )
                    deriv_gen = FiniteDifferenceDerivative(potential_function,
                                                           function_shape=((None, None), 0),
                                                           stencil=5 + expansion_order['potential'],
                                                           mesh_spacing=1e-3,
                                                           ).derivatives(coords)
                    potential_derivatives = deriv_gen.derivative_tensor(list(range(1, order + 3)))
                    wat['potential_derivatives'] = potential_derivatives
        else:
            potential_derivatives = []


        # analyzer = VPTRunner.run_simple(
        #     [atoms, coords, dict(masses=masses)],
        #     2,
        #     potential_derivatives=potential_derivatives,
        #     calculate_intensities=False,
        #     order=order,
        #     include_potential=include_potential,
        #     include_gmatrix=include_gmatrix,
        #     include_pseudopotential=include_pseudopotential,
        #     include_coriolis_coupling=include_gmatrix
        # )
        # analyzer.print_output_tables(print_intensities=False, print_energies=True)

        # try:
        #     os.remove(os.path.expanduser('~/Desktop/water_analyt.hdf5'))
        # except:
        #     pass
        analyzer = VPTRunner.run_simple(
            [atoms, coords, dict(masses=masses)],
            2,
            potential_derivatives=potential_derivatives,
            calculate_intensities=False,
            internals=[
                [0, -1, -1, -1],
                [1,  0, -1, -1],
                [2,  0,  1, -1]
            ],
            order=order,
            internal_fd_mesh_spacing=1e-2,
            cartesian_fd_mesh_spacing=1e-2,
            checkpoint=os.path.expanduser('~/Desktop/water_analyt.hdf5'),
            expansion_order=expansion_order
        )
```
#### <a name="HOHCorrectedDegeneracies">HOHCorrectedDegeneracies</a>
```python
    def test_HOHCorrectedDegeneracies(self):
        VPTRunner.run_simple(
            TestManager.test_data('HOH_freq.fchk'),
            2,
            zero_order_energy_corrections=[
                [(0, 1, 0), (4681.56364+3800) * UnitsData.convert("Wavenumbers", "Hartrees")],
                [(0, 2, 0), (4681.56364+7800) * UnitsData.convert("Wavenumbers", "Hartrees")],
                [(0, 0, 2), (4681.56364+7801) * UnitsData.convert("Wavenumbers", "Hartrees")],
                [(0, 3, 1), (4681.56364+3801) * UnitsData.convert("Wavenumbers", "Hartrees")],
            ],
            degeneracy_specs={
                'wfc_threshold':.3,
                'extra_groups':[
                    [[0, 2, 0], [0, 1, 1]],
                    [[0, 0, 2], [0, 3, 1]],
                ]
            }
            # operator_chunk_size=int(12)
        )
```
#### <a name="WaterSkippedCouplings">WaterSkippedCouplings</a>
```python
    def test_WaterSkippedCouplings(self):

        VPTRunner.run_simple(
            TestManager.test_data('water_freq.fchk'),
            1,
            degeneracy_specs='auto',
            operator_coefficient_threshold=(1.0e-8)
        )
```
#### <a name="H2COPolyads">H2COPolyads</a>
```python
    def test_H2COPolyads(self):

        internals = [
            [0, -1, -1, -1],
            [1, 0, -1, -1],
            [2, 1, 0, -1],
            [3, 1, 0, 2]
        ]

        VPTRunner.run_simple(
            TestManager.test_data('OCHH_freq.fchk'),
            2,
            degeneracy_specs=True
            # degeneracy_specs={
            #     'nT':[1, 1, 1, 1, 2, 2]
            # }
        )
```
#### <a name="H2COModeSel">H2COModeSel</a>
```python
    def test_H2COModeSel(self):

        internals = [
            [0, -1, -1, -1],
            [1,  0, -1, -1],
            [2,  1,  0, -1],
            [3,  1,  0,  2]
        ]

        VPTRunner.run_simple(
            TestManager.test_data('OCHH_freq.fchk'),
            1,
            degeneracy_specs=None,
            mode_selection=[1, 2, 3, 4, 5]
        )
```
#### <a name="HODRephase">HODRephase</a>
```python
    def test_HODRephase(self):
        VPTRunner.run_simple(
            TestManager.test_data('HOD_freq_16.fchk'),
            1,
            degeneracy_specs=None,
            # order=4,
            # expansion_order=2
        )
```
#### <a name="HOHRephase">HOHRephase</a>
```python
    def test_HOHRephase(self):
        VPTRunner.run_simple(
            TestManager.test_data('HOH_freq.fchk'),
            1,
            degeneracy_specs=None,
            # order=4,
            # expansion_order=2
        )
```
#### <a name="NH3">NH3</a>
```python
    def test_NH3(self):

        VPTRunner.run_simple(
            TestManager.test_data('NH3_freq.fchk'),
            2,
            # degeneracy_specs=False,
            order=4,
            expansion_order=2,

            # basis_filters={
            #     'max_quanta':[2, -1, -1, -1, -1, -1]
            # }
        )
```
#### <a name="HOONO">HOONO</a>
```python
    def test_HOONO(self):

        VPTRunner.run_simple(
            TestManager.test_data('HOONO_freq.fchk'),
            1,
            degeneracy_specs=None,
            # order=4,
            # expansion_order=2
        )
```
#### <a name="H2COSkippedCouplings">H2COSkippedCouplings</a>
```python
    def test_H2COSkippedCouplings(self):

        VPTRunner.run_simple(
            TestManager.test_data('OCHH_freq.fchk'),
            1,
            degeneracy_specs='auto',
            operator_coefficient_threshold=1.00 / 219475
        )
```
#### <a name="WaterDimerSkippedCouplings">WaterDimerSkippedCouplings</a>
```python
    def test_WaterDimerSkippedCouplings(self):

        COM = -3
        A = -2
        C = -1
        X = 1000
        LHF = 0
        LO = 1
        SH = 2
        RO = 3
        RH1 = 4
        RH2 = 5

        internals = [
            [LHF, X, X, X],
            [LO, LHF, X, X],
            [SH, LO, LHF, X],
            [RH2, SH, LO, LHF],  # get out of plane
            [RO, LO, RH2, LHF],
            [RH1, RO, RH2, LHF]
        ]

        VPTRunner.run_simple(
            TestManager.test_data('water_dimer_freq.fchk'),
            1,
            degeneracy_specs='auto',
            # operator_coefficient_threshold=0.005/219475,
            # basis_filters=[
            #     {'max_quanta': [3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]}
            # ],
            # internals=internals,
            mode_selection=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
        )
```

 </div>
</div>

___

[Edit Examples](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Runner/VPTRunner.md) or 
[Create New Examples](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Runner/VPTRunner.md) <br/>
[Edit Template](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Runner/VPTRunner.md) or 
[Create New Template](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Runner/VPTRunner.md) <br/>
[Edit Docstrings](https://github.com/McCoyGroup/Psience/edit/master/Psience/VPT2/Runner.py#L942?message=Update%20Docs)