## <a id="Psience.VPT2.Runner.VPTStateSpace">VPTStateSpace</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Runner.py#L169)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Runner.py#L169?message=Update%20Docs)]
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
__init__(self, states, degeneracy_specs=None, system=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Runner/VPTStateSpace.py#L256)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Runner/VPTStateSpace.py#L256?message=Update%20Docs)]
</div>

  - `states`: `list | int`
    > A list of states or a number of quanta to target
  - `degeneracy_specs`: `'auto' | list | dict`
    > A specification of degeneracies, either as polyads, explicit groups of states, or parameters to a method. (see Details for more info)


<a id="Psience.VPT2.Runner.VPTStateSpace.from_system_and_quanta" class="docs-object-method">&nbsp;</a> 
```python
from_system_and_quanta(system, quanta, target_modes=None, only_target_modes=False, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Runner/VPTStateSpace.py#L295)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Runner/VPTStateSpace.py#L295?message=Update%20Docs)]
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
get_state_list_from_quanta(n_quanta, n_modes, target_modes=None, only_target_modes=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Runner/VPTStateSpace.py#L324)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Runner/VPTStateSpace.py#L324?message=Update%20Docs)]
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
build_degenerate_state_spaces(self, degeneracy_specs, states, system=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Runner/VPTStateSpace.py#L354)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Runner/VPTStateSpace.py#L354?message=Update%20Docs)]
</div>

  - `degeneracy_specs`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Runner.VPTStateSpace.filter_generator" class="docs-object-method">&nbsp;</a> 
```python
filter_generator(self, target_property, order=2, postfilters=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Runner/VPTStateSpace.py#L396)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Runner/VPTStateSpace.py#L396?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Runner.VPTStateSpace.get_filter" class="docs-object-method">&nbsp;</a> 
```python
get_filter(self, target_property, order=2, postfilters=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Runner/VPTStateSpace.py#L400)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Runner/VPTStateSpace.py#L400?message=Update%20Docs)]
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
get_state_space_filter(states, n_modes=None, order=2, target='wavefunctions', postfilters=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Runner/VPTStateSpace.py#L417)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Runner/VPTStateSpace.py#L417?message=Update%20Docs)]
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
## <a class="collapse-link" data-toggle="collapse" href="#Details-54f472" markdown="1"> Details</a> <a class="float-right" data-toggle="collapse" href="#Details-54f472"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Details-54f472" markdown="1">
 
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


## Examples













<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#Tests-f13b8a" markdown="1"> Tests</a> <a class="float-right" data-toggle="collapse" href="#Tests-f13b8a"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Tests-f13b8a" markdown="1">
 - [HOHVPTRunnerFlow](#HOHVPTRunnerFlow)
- [GetDegenerateSpaces](#GetDegenerateSpaces)

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
### <a class="collapse-link" data-toggle="collapse" href="#Setup-a66fcc" markdown="1"> Setup</a> <a class="float-right" data-toggle="collapse" href="#Setup-a66fcc"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Setup-a66fcc" markdown="1">
 
Before we can run our examples we should get a bit of setup out of the way.
Since these examples were harvested from the unit tests not all pieces
will be necessary for all situations.

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
[Edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Runner.py#L169?message=Update%20Docs)   
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