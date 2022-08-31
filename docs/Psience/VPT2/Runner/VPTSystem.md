## <a id="Psience.VPT2.Runner.VPTSystem">VPTSystem</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Runner.py#L30)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Runner.py#L30?message=Update%20Docs)]
</div>

Provides a little helper for setting up the input
system for a VPT job







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.VPT2.Runner.VPTSystem.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, mol, internals=None, dummy_atoms=None, modes=None, mode_selection=None, potential_derivatives=None, potential_function=None, order=2, dipole_derivatives=None, eckart_embed=False): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Runner/VPTSystem.py#L60)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Runner/VPTSystem.py#L60?message=Update%20Docs)]
</div>

  - `mol`: `str | list | Molecule`
    > the molecule or system specification to use (doesn't really even need to be a molecule)
  - `internals`: `list | dict`
    > the Z-matrix for the internal coordinates optionally with a specification of a `conversion` and `inverse`
To supply a conversion function, provide a `dict` like so
```python
{
    'zmatrix': [[atom1, bond1, angle1, dihed1], [atom2, bond2, angle2, dihed2], ...] or None,
    'conversion': 'a function to convert from Z-matrix coordinates to desired coordinates',
    'inverse': 'the inverse conversion'
}
```
  - `modes`: `MolecularVibrations|dict`
    > the normal modes to use if not already supplied by the Molecule
  - `potential_derivatives`: `Iterable[np.ndarray]`
    > the derivatives of the potential to use for expansions
  - `dipole_derivatives`: `Iterable[np.ndarray]`
    > the set of dipole derivatives to use for expansions


<a id="Psience.VPT2.Runner.VPTSystem.nmodes" class="docs-object-method">&nbsp;</a> 
```python
@property
nmodes(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Runner/VPTSystem.py#L131)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Runner/VPTSystem.py#L131?message=Update%20Docs)]
</div>
Provides the number of modes in the system
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Runner.VPTSystem.get_potential_derivatives" class="docs-object-method">&nbsp;</a> 
```python
get_potential_derivatives(self, potential_function, order=2, **fd_opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Runner/VPTSystem.py#L144)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Runner/VPTSystem.py#L144?message=Update%20Docs)]
</div>
Computes potential derivatives for the given function through finite difference
  - `potential_function`: `Any`
    > 
  - `order`: `Any`
    > 
  - `fd_opts`: `Any`
    > 
  - `:returns`: `_`
    >


<a id="Psience.VPT2.Runner.VPTSystem.from_harmonic_scan" class="docs-object-method">&nbsp;</a> 
```python
from_harmonic_scan(scan_array): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Runner/VPTSystem.py#L165)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Runner/VPTSystem.py#L165?message=Update%20Docs)]
</div>
 </div>
</div>



<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#Details-261327" markdown="1"> Details</a> <a class="float-right" data-toggle="collapse" href="#Details-261327"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Details-261327" markdown="1">
 When using functions of internal (Z-matrix/polyspherical) coordinates, a sample form of the conversion function is
```python
def conv(r, t, f, **kwargs):
    '''
    Takes the bond lengths (`r`), angles `(t)` and dihedrals `(f)`,
    and returns new coordinates that are functions of these coordinates
    '''
    ... # convert the coordinates
    return np.array([r, t, f])
```
and then the inverse function will take the output of `conv` and return the original Z-matrix/polyspherical coordinates.
 </div>
</div>


## Examples













<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#Tests-9eceb7" markdown="1"> Tests</a> <a class="float-right" data-toggle="collapse" href="#Tests-9eceb7"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Tests-9eceb7" markdown="1">
 - [HOHVPTRunnerFlow](#HOHVPTRunnerFlow)

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
### <a class="collapse-link" data-toggle="collapse" href="#Setup-bca53c" markdown="1"> Setup</a> <a class="float-right" data-toggle="collapse" href="#Setup-bca53c"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Setup-bca53c" markdown="1">
 
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Runner/VPTSystem.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Runner/VPTSystem.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Runner/VPTSystem.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Runner/VPTSystem.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Runner.py#L30?message=Update%20Docs)   
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