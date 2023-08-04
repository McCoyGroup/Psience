## <a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver">AnalyticPerturbationTheorySolver</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Analytic.py#L17)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Analytic.py#L17?message=Update%20Docs)]
</div>

A re-attempt at using the recursive expressions
to provide simpler code for getting APT expressions







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, hamiltonian_expansion): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L22)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L22?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.from_order" class="docs-object-method">&nbsp;</a> 
```python
from_order(order, internals=True): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L25)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L25?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.energy_correction" class="docs-object-method">&nbsp;</a> 
```python
energy_correction(self, order): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L72)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L72?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.wavefunction_correction" class="docs-object-method">&nbsp;</a> 
```python
wavefunction_correction(self, order): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L75)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L75?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.overlap_correction" class="docs-object-method">&nbsp;</a> 
```python
overlap_correction(self, order): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L78)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L78?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.operator_correction" class="docs-object-method">&nbsp;</a> 
```python
operator_correction(self, order, operator_type=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L81)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L81?message=Update%20Docs)]
</div>


<a id="Psience.VPT2.Analytic.AnalyticPerturbationTheorySolver.operator_expansion_terms" class="docs-object-method">&nbsp;</a> 
```python
operator_expansion_terms(order, operator_type=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L84)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Analytic/AnalyticPerturbationTheorySolver.py#L84?message=Update%20Docs)]
</div>
 </div>
</div>




## Examples













<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#Tests-991a6c" markdown="1"> Tests</a> <a class="float-right" data-toggle="collapse" href="#Tests-991a6c"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Tests-991a6c" markdown="1">
 - [AnalyticPTOperators](#AnalyticPTOperators)

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
### <a class="collapse-link" data-toggle="collapse" href="#Setup-20a4ee" markdown="1"> Setup</a> <a class="float-right" data-toggle="collapse" href="#Setup-20a4ee"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Setup-20a4ee" markdown="1">
 
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

#### <a name="AnalyticPTOperators">AnalyticPTOperators</a>
```python
    def test_AnalyticPTOperators(self):

        internals = False
        vpt2 = AnalyticPerturbationTheorySolver.from_order(2, internals=internals)
        # print(
        #     vpt2.hamiltonian_expansion[2].get_poly_terms([])
        # )
        # raise Exception(...)
        #
        # H1PH1 = vpt2.energy_correction(2).expressions[1]
        #
        # H1 = vpt2.hamiltonian_expansion[1]
        # H1H1 = H1*H1

        """
        ==================== V[1](0, 0, 0)V[1](0, 0, 0) 1 ====================
  :: 1
   > 1 [array([0.   , 0.   , 0.   , 1.125])]
   > 1.1249999999999993 [array([1., 3., 3., 1.])]
   > 0.7499999999999996 [array([1.        , 1.83333333, 1.        , 0.16666667])]
   > 1 [array([ 0.   ,  0.25 , -0.375,  0.125])]
   """

        # H1H1.get_poly_terms([],
        #                     # allowed_paths=[
        #                     #     ((1,), (-1,)),
        #                     #     ((3,), (-3,))
        #                     # ]
        #                     ).prune_operators([(1, 1)]).print_tree()
        # raise Exception(...)

        PH1 = vpt2.wavefunction_correction(1).expressions[0]
        E2 = vpt2.energy_correction(2)

        # # pt_shit = H1PH1.get_poly_terms((4,)).combine()
        # pt_shit = PH1((1,), simplify=False).expr.combine()

        # raise Exception(
        #     PH1.get_poly_terms((1,)).terms[((1,0, 0, 0, 0),)].terms[((1,),)].coeffs
        # )

        # pt_shit = E2.expressions[1]([], simplify=False).expr
        # raise Exception(pt_shit)

        # subpoly = pt_shit.terms[((1, 1, 0, 2, 1), (1, 1, 1, 0, 2))].combine(combine_energies=False).terms[((1, 1, 1),)]
        # for p in subpoly.polys:
        #     print(p.prefactor, p.coeffs)
        # raise Exception(...)

        # h3_poly = vpt2.hamiltonian_expansion[1]([3]).expr.terms[((1, 0, 0, 0, 0),)]
        # raise Exception(
        #     h3_poly.prefactor,
        #     h3_poly.coeffs
        # )
        # with np.printoptions(linewidth=1e8):
        #     pt_shit = H1H1([]).expr
        #     pt_shit.print_tree()
        #
        # raise Exception(...)

        # # raise Exception(
        # #     Y1.get_poly_terms((1,))
        # #     # [
        # #     #     [p.prefactor for p in psum.polys]
        # #     #     for psum in h1y1.get_poly_terms(()).terms.values()
        # #     # ]
        # # )
        #
        # h1y1 = E2.expressions[1]
        # # raise Exception(h1y1, h1y1.gen2, h1y1.gen2.expressions)
        #
        # raise Exception(
        #     h1y1.get_poly_terms(()).combine()
        #     # [
        #     #     [p.prefactor for p in psum.polys]
        #     #     for psum in h1y1.get_poly_terms(()).terms.values()
        #     # ]
        # )

        # for _ in vpt2.energy_correction(2).expressions[1].changes[()]:
        #     print(_)

        # test_sum = list(jesus_fuck.terms.values())[1]
        # # for poly in test_sum.polys:
        # #     print(poly.coeffs)
        # raise Exception(test_sum, test_sum.combine())
        # raise Exception(test_sum.polys[0].combine(test_sum.polys[1]).coeffs)


        # load H20 parameters...
        file_name = "HOH_freq.fchk"
        from Psience.BasisReps import HarmonicOscillatorMatrixGenerator
        HarmonicOscillatorMatrixGenerator.default_evaluator_mode = 'rho'
        runner, _ = VPTRunner.construct(
            TestManager.test_data(file_name),
            1,
            internals=(
                [[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]]
                    if internals else
                None
            ),
            logger=False,
            zero_element_warning=False
            # mode_selection=[0, 1]
        )
        ham = runner.hamiltonian
        V = ham.V_terms
        G = ham.G_terms
        U = ham.pseudopotential_term
        g2 = G[2]

        if internals:
            water_expansion = [
                [V[0]/2,  G[0]/2],
                [
                    # np.zeros(V[1].shape),
                    np.sum([V[1].transpose(p) for p in itertools.permutations([0, 1, 2])], axis=0)/np.math.factorial(3)/6,
                    # np.zeros(V[1].shape),
                    -np.moveaxis(G[1], -1, 0)/2
                        if isinstance(G[1], np.ndarray) else
                    np.zeros(V[1].shape)
                ],
                [
                    # np.zeros(V[2].shape),
                    V[2]/24,
                    # np.zeros(V[2].shape),
                    -np.moveaxis(G[2], -1, 0)/4
                        if isinstance(G[2], np.ndarray) else
                    np.zeros(V[2].shape),
                    # 0
                    U[0]/8
                ]
            ]
        else:
            Z = ham.coriolis_terms
            water_expansion = [
                [V[0] / 2, G[0] / 2],
                [
                    np.zeros(V[1].shape),
                    # np.sum(
                    #     [V[1].transpose(p) for p in itertools.permutations([0, 1, 2])],
                    #     axis=0
                    # ) / np.math.factorial(3) / 6,
                    0 # G
                ],
                [
                    # np.zeros(V[2].shape),
                    V[2] / 24,
                    0, # G
                    0, # V'
                    # np.zeros(V[2].shape),
                    -Z[0],
                    # 0
                    U[0] / 8 # Watson
                ]
            ]
        # raise Exception(
        #     -.25 * np.array([
        #         Z[0][0, 0, 2, 2],
        #         Z[0][0, 2, 0, 2],
        #         Z[0][0, 2, 2, 0],
        #         Z[0][2, 0, 0, 2],
        #         Z[0][2, 0, 2, 0],
        #         Z[0][2, 2, 0, 0]
        #     ]) * 219475
        # )
        # raise Exception(
        #     -.25*np.array([
        #             Z[0][0, 0, 2, 2],
        #             Z[0][0, 2, 0, 2],
        #             Z[0][0, 2, 2, 0],
        #             Z[0][2, 0, 0, 2],
        #             Z[0][2, 0, 2, 0],
        #             Z[0][2, 2, 0, 0]
        #     ]) * 219475
        # )

        water_freqs = ham.modes.freqs

        # raise Exception(list(sorted([
        #     [p, np.linalg.norm((np.transpose(water_expansion[2][3], p) - water_expansion[2][3]).flatten())]
        #     for p in itertools.permutations([0, 1, 2, 3])
        #     ],
        #     key=lambda x:x[1]
        # )))

        # solver = runner.hamiltonian.get_solver(runner.states.state_list)
        # raise Exception(
        #     solver.representations[1].todense()[0][:4],
        #     solver.flat_total_space.excitations[:4]
        # )

        # wfns = runner.get_wavefunctions()

        # ham.G_terms = [
        #     G[0],
        #     # G[1],
        #     np.zeros(G[1].shape),
        #     # G[2]
        #     np.zeros(G[2].shape),
        # ]
        # ham.V_terms = [
        #     V[0],
        #     # np.sum([V[1].transpose(p) for p in itertools.permutations([0, 1, 2])], axis=0) / np.math.factorial(3),
        #     np.zeros(V[1].shape),
        #     V[2],
        #     # np.zeros(V[2].shape)
        # ]
        # ham.coriolis_terms = [
        #     np.zeros(V[2].shape)
        #     # (Z[0] + np.transpose(Z[0], [0, 3, 2, 1])) / 2
        # ]
        # ham.pseudopotential_term = [0]
        runner.print_tables(print_intensities=False)
        """
        :: State    <0|dH(2)|0>  <0|dH(1)|1> 
          0 0 0      0.00000   -141.74965
          0 0 1      0.00000   -514.12544
          0 1 0      0.00000   -509.59810
          1 0 0      0.00000   -165.17359
  
        :: State    <0|dH(2)|0>  <0|dH(1)|1> 
          0 0 0     81.09533   -156.70145
          0 0 1    276.68965   -545.08549
          0 1 0    281.61581   -538.53511
          1 0 0     88.49014   -213.67502
          
        :: State    <0|dH(2)|0>  <0|dH(1)|1>  # Cartesians
          0 0 0     42.02414   -117.62974
          0 0 1    199.53514   -467.93034
          0 1 0    195.87334   -452.78667
          1 0 0    -32.13943    -93.04841
        >>--------------------------------------------------<<
        >>------------------------- States Energies -------------------------
        :: State     Harmonic   Anharmonic     Harmonic   Anharmonic
                       ZPE          ZPE    Frequency    Frequency
        0 0 0   4681.56362   4605.95750            -            - 
        0 0 1            -            -   3937.52464   3744.73491 
        0 1 0            -            -   3803.29958   3621.98639 
        1 0 0            -            -   1622.30304   1572.72428 
        >>--------------------------------------------------<<
        """

        # V1 = water_expansion[1][0]
        # G1 = water_expansion[1][1]
        # print(G[1]/2)
        # print(G1)
        # raise Exception(...)

        # raise Exception(
        #     wfns.corrs.wfn_corrections[1].todense()[0][:4],
        #     wfns.corrs.total_basis.excitations[:4],
        #     -7.14635718e-03 * water_freqs[0],
        #     -5.83925920e-03 * water_freqs[0],
        #     0.02061049 * water_freqs[0],
        #     (
        #             V1[0, 0, 0]*1.060660171779821
        #             +(V1[0, 1, 1] + V1[0, 2, 2])*1.060660171779821
        #             + G1[0, 0, 0]*-0.35355339
        #             + (G1[1, 0, 1] + G1[2, 0, 2])*-0.35355339
        #     ) / water_freqs[0]
        #     # # wfns.corrs.wfn_corrections[1].todense()[0, 1],
        #     # # should be .005839259198189573
        #     # V[1][0, 0, 0]/water_freqs[0]*1.06066017
        #     # # + G[1][0, 0, 0]/water_freqs[0]*-0.35355339
        #     # + (V[1][0, 1, 1] + V[1][0, 2, 2])/water_freqs[0]*0.70710678*.5
        #     # # + (G[1][0, 1, 1] + G[1][0, 2, 2])/water_freqs[0]*0.70710678*-.5
        # )
        # 0.02819074222400865, -0.035337103975354986, [-0.00758025182870669, -0.027756852139690837]

        # wfns = runner.get_wavefunctions()
        # Y1 = vpt2.wavefunction_correction(1)
        # raise Exception(
        #     wfns.corrs.wfn_corrections[1].todense()[0][:10],
        #     wfns.corrs.total_basis.excitations[:10],
        #     Y1([3,]).evaluate([0, 0, 0], water_expansion, water_freqs),
        #     Y1([2, 1]).evaluate([0, 0, 0], water_expansion, water_freqs),
        # )

        with np.printoptions(linewidth=1e8):
            jesus_fuck = E2([])
            # jesus_fuck.expr.print_tree()
            corr = jesus_fuck.evaluate([0, 0, 0], water_expansion, water_freqs, verbose=True) * UnitsData.convert("Hartrees", "Wavenumbers")
        print(corr)
        raise Exception(corr)
            # vpt2.energy_correction(2).expressions[1].changes[()]

        raise Exception(
            AnalyticPerturbationTheoryDriver.from_order(2).energy_correction_driver(
                2
            ).get_poly_evaluator(
                [ # tensor coeffs
                    [
                        np.eye(3),
                        np.eye(3),
                    ],
                    [
                        np.ones((3, 3, 3)),
                        np.ones((3, 3, 3))
                    ],
                    [
                        np.ones((3, 3, 3, 3)),
                        np.ones((3, 3, 3, 3)),
                        1
                    ]
                ],
                [1, 1, 1] # freqs
            )
        )

        corrections = AnalyticPTCorrectionGenerator([
            [
                ['x', 'x', 'x'],
                ['p', 'x', 'p']
            ],
            # ['x'],
            [
                ['x', 'x', 'x'],
                ['p', 'x', 'p']
            ]
        ]).get_correction([2])

        v1 = np.ones((3, 3, 3))
        g1 = np.ones((3, 3, 3))

        return corrections.evaluate([0, 0, 0], [[v1, g1], [v1, g1]], [1, 1, 1], 1)

        raise Exception(corrections)

        coeffs = np.array([
            TensorCoeffPoly({((1, 0, 0),):2, ((0, 1, 0),):1}),
            TensorCoeffPoly({((0, 0, 1),):1}),
        ], dtype=object)

        # new_b_poly = np.dot(
        #             [[1, 3], [2, 1]],
        #             coeffs
        #         )
        # raise Exception(
        #     np.dot(
        #         [[-1, 3], [2, -1]],
        #         np.dot(
        #             [[1, 3], [2, 1]],
        #             coeffs
        #         )
        #     )/5
        # )

        # from Psience.AnalyticModels import AnalyticModel
        #
        # AnalyticModel(
        #     [
        #         AnalyticModel.r(0, 1),
        #
        #         ]
        # ).g()




        from McUtils.Zachary import DensePolynomial

        shifted = DensePolynomial([
            [1, 2, 3, 4, 5],
            [1, -2, 0, 6, 8]
        ]).shift([2, 3])

        # raise Exception(shifted.deriv(1).coeffs)

        self.assertTrue(
            np.allclose(
                shifted.coeffs,
                [
                    [2157., 2716., 1281., 268., 21.],
                    [ 805., 1024.,  486., 102.,  8.]
                ]
            )
        )

        #
        new = DensePolynomial(coeffs)*DensePolynomial(coeffs)
        raise Exception(new)

 #        """
 #        DensePolynomial([TensorCoeffPoly({((1, 0, 0), (1, 0, 0)): 4, ((0, 1, 0), (1, 0, 0)): 4, ((0, 1, 0), (0, 1, 0)): 1},1)
 # TensorCoeffPoly({((0, 0, 1), (0, 1, 0)): 2, ((0, 0, 1), (1, 0, 0)): 4},1)
 # TensorCoeffPoly({((0, 0, 1), (0, 0, 1)): 1},1)], 1)"""

        # raise Exception(
        #     PTPoly(coeffs)*
        #     PTPoly(coeffs)
        # )

        base_classes = RaisingLoweringClasses(14, [2, 4, 2])
        print(list(base_classes))
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/VPT2/Analytic/AnalyticPerturbationTheorySolver.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/VPT2/Analytic.py#L17?message=Update%20Docs)   
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