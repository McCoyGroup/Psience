
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

class VPT2Tests(TestCase):

    # region Water Analogs

    @validationTest
    def test_HOHVPTRunner(self):

        file_name = "HOH_freq.fchk"
        VPTRunner.run_simple(
            TestManager.test_data(file_name),
            3,
            memory_constrained=True,
            logger=True
        )

    @validationTest
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

    @validationTest
    def test_HOHVPTRunnerShifted(self):

        file_name = "HOH_freq.fchk"
        VPTRunner.run_simple(
            TestManager.test_data(file_name),
            3,
            logger=True,
            corrected_fundamental_frequencies=np.array([1600, 3775, 3880])/UnitsData.convert("Hartrees", "Wavenumbers")
        )

    @validationTest
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

    @validationTest
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

    #endregion

    @validationTest
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

    @debugTest
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

        # model.run_VPT(order=order, return_analyzer=False, expansion_order=expansion_order)
        """
        ::> Energy Corrections (4)
  > State    <0|dH(2)|0>  <0|dH(1)|1>  <0|dH(4)|0>  <0|dH(3)|1>  <0|dH(2)|2>  <0|dH(1)|3> 
  0 0 0     58.77245   -129.09822      0.74196     -4.07776      2.09670      0.72134
  0 0 1    212.70254   -460.34848      4.15867    -20.89442      1.50185     14.23673
  0 1 0    211.20631   -449.30152      4.26017    -20.40002     -5.93737     22.03242
  1 0 0     70.58602   -158.64843      0.61308     -5.00042      2.19406      1.11387
  0 0 2    441.32022   -953.29021     11.32184    -56.56544    -93.85290    200.31423
  0 2 0    436.14405   -925.36339     11.41748    -54.75476    -49.40033     31.53122
  2 0 0     82.17782   -197.39848      0.48565     -6.18033      1.37629      1.82257
  0 1 1    512.27353  -1094.10494     15.05825    -73.95472   -115.61485    174.62554
  1 0 1    235.86395   -515.46315      4.64152    -23.88543     -1.33125     18.74167
  1 1 0    235.74270   -493.98804      4.71744    -23.09038    -16.06563     34.66543
<::
::> States Energies (4)
  > State     Harmonic   Anharmonic     Harmonic   Anharmonic
               ZPE          ZPE    Frequency    Frequency
0 0 0   4671.98118   4601.13767            -            - 
0 0 1            -            -   3898.39551   3720.59591 
0 1 0            -            -   3840.97965   3673.68316 
1 0 0            -            -   1604.58720   1586.28889 
0 0 2            -            -   7796.79102   7416.88228 
0 2 0            -            -   7681.95930   7202.37707 
2 0 0            -            -   3209.17441   3162.30144 
0 1 1            -            -   7739.37516   7228.50149 
1 0 1            -            -   5502.98272   5292.39354 
1 1 0            -            -   5445.56685   5258.39189 
<::
"""

        """"
::> Energy Corrections (4 - int)
  > State    <0|dH(2)|0>  <0|dH(1)|1>  <0|dH(4)|0>  <0|dH(3)|1>  <0|dH(2)|2>  <0|dH(1)|3> 
  0 0 0     59.35897   -129.04702     -3.54396     -3.82024      2.09509      0.97740
  0 0 1    213.30534   -460.29036     -0.87381    -20.39105      1.52004     14.71216
  0 1 0    211.73506   -448.90679    -10.81118    -20.02599     -6.03063     22.48989
  1 0 0     71.25746   -158.57990     -7.36846     -4.71729      2.19070      1.39454
  0 0 2    441.94349   -953.21967     10.03720    -55.85125    -95.05966    203.83185
  0 2 0    436.56498   -924.67248    -15.47662    -54.42478    -48.42560     29.36111
  2 0 0     83.00399   -197.50757    -13.63095     -5.88012      1.33945      2.13295
  0 1 1    512.80780  -1093.52329    -13.02580    -73.42807   -115.70497    175.14432
  1 0 1    236.58663   -515.56485     -4.32149    -23.38510     -1.24587     19.13576
  1 1 0    236.35167   -492.97424    -22.53337    -22.74588    -16.44501     35.35035
<::
::> States Energies
  > State     Harmonic   Anharmonic     Harmonic   Anharmonic
               ZPE          ZPE    Frequency    Frequency
0 0 0   4671.50834   4597.52858            -            - 
0 0 1            -            -   3898.07519   3720.03727 
0 1 0            -            -   3841.98586   3664.41598 
1 0 0            -            -   1602.95562   1581.11244 
0 0 2            -            -   7796.15039   7421.81210 
0 2 0            -            -   7683.97173   7180.87810 
2 0 0            -            -   3205.91124   3149.34874 
0 1 1            -            -   7740.06106   7206.31081 
1 0 1            -            -   5501.03081   5286.21566 
1 1 0            -            -   5444.94148   5235.92476 
<::
        """


        """
=== No U Anal ===
::> Energy Corrections
  > State    <0|dH(2)|0>  <0|dH(1)|1>  <0|dH(4)|0>  <0|dH(3)|1>  <0|dH(2)|2>  <0|dH(1)|3> 
  0 0 0     79.40046   -129.09822      1.33074     -4.13575      2.09670      0.66336
  0 0 1    233.33055   -460.34848      5.20152    -20.98900      1.50185     14.14215
  0 1 0    231.83432   -449.30152      5.13964    -20.53773     -5.93737     21.89471
  1 0 0     91.21403   -158.64843      1.63462     -5.05808      2.19406      1.05622
  0 0 2    461.94823   -953.29021     12.81878    -56.69661    -93.85290    200.18307
  0 2 0    456.77206   -925.36339     12.58764    -54.97220    -49.40033     31.31378
  2 0 0    102.80583   -197.39848      1.93994     -6.23764      1.37629      1.76525
  0 1 1    532.90154  -1094.10494     16.39181    -74.12902   -115.61485    174.45124
  1 0 1    256.49196   -515.46315      6.11714    -23.97967     -1.33125     18.64743
  1 1 0    256.37071   -493.98804      6.02967    -23.22776    -16.06563     34.52806
<::
::> States Energies
  > State     Harmonic   Anharmonic     Harmonic   Anharmonic
               ZPE          ZPE    Frequency    Frequency
0 0 0   4671.98118   4622.23847            -            - 
0 0 1            -            -   3898.39551   3720.97682 
0 1 0            -            -   3840.97965   3673.81440 
1 0 0            -            -   1604.58720   1586.72233 
0 0 2            -            -   7796.79102   7417.64409 
0 2 0            -            -   7681.95930   7202.63956 
2 0 0            -            -   3209.17441   3163.16832 
0 1 1            -            -   7739.37516   7229.01364 
1 0 1            -            -   5502.98272   5293.20789 
1 1 0            -            -   5445.56685   5258.95657 
<::
=== No U Int ===
::> Energy Corrections
  > State    <0|dH(2)|0>  <0|dH(1)|1>  <0|dH(4)|0>  <0|dH(3)|1>  <0|dH(2)|2>  <0|dH(1)|3> 
  0 0 0     79.41128   -129.04702      6.65695     -4.16068      2.09509      0.63696
  0 0 1    233.35765   -460.29036     14.67717    -21.06033      1.52004     14.04288
  0 1 0    231.78737   -448.90679      4.94286    -20.67249     -6.03063     21.84339
  1 0 0     91.30978   -158.57990     12.33106     -5.10371      2.19070      1.00812
  0 0 2    461.99580   -953.21967     30.93825    -56.84937    -95.05966    202.83373
  0 2 0    456.61729   -924.67249      5.83055    -55.37734    -48.42560     28.40856
  2 0 0    103.05630   -197.50757     15.56717     -6.31253      1.33945      1.70055
  0 1 1    532.86011  -1093.52329      8.07831    -74.40341   -115.70497    174.16899
  1 0 1    256.63894   -515.56485     20.72810    -24.10036     -1.24587     18.42050
  1 1 0    256.40398   -492.97424      2.71928    -23.43836    -16.44501     34.65787
<::
::> States Energies
  > State     Harmonic   Anharmonic     Harmonic   Anharmonic
               ZPE          ZPE    Frequency    Frequency
0 0 0   4671.50834   4627.10092            -            - 
0 0 1            -            -   3898.07519   3724.72966 
0 1 0            -            -   3841.98586   3669.35699 
1 0 0            -            -   1602.95562   1590.51908 
0 0 2            -            -   7796.15039   7431.19689 
0 2 0            -            -   7683.97173   7190.76012 
2 0 0            -            -   3205.91124   3168.16203 
0 1 1            -            -   7740.06106   7215.94422 
1 0 1            -            -   5501.03081   5300.31470 
1 1 0            -            -   5444.94148   5250.27242 

"""

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
        """
        ::> States Energies
  > State     Harmonic   Anharmonic     Harmonic   Anharmonic
               ZPE          ZPE    Frequency    Frequency
0 0 0   4671.50834   4601.82020            -            - 
0 0 1            -            -   3898.07520   3720.77819 
0 1 0            -            -   3841.98587   3674.50193 
1 0 0            -            -   1602.95562   1585.32112 
0 0 2            -            -   7796.15039   7354.56223 
0 2 0            -            -   7683.97173   7265.55168 
2 0 0            -            -   3205.91124   3161.09541 
0 1 1            -            -   7740.06106   7229.03319 
1 0 1            -            -   5501.03081   5291.74053 
1 1 0            -            -   5444.94148   5258.00648

::> Energy Corrections
  > State    <0|dH(2)|0>  <0|dH(1)|1>  <0|dH(4)|0>  <0|dH(3)|1>  <0|dH(2)|2>  <0|dH(1)|3> 
  0 0 0     24.77485    -94.46299     -0.18449     -1.25972      1.33962      0.33167
  0 0 1    141.19064   -388.17578      0.24609     -9.72724      1.89569      8.29192
  0 1 0    128.27090   -365.44296     -0.06078     -7.76138     -7.82200     16.86981
  1 0 0    -25.47576    -61.84687     -0.83659      1.32999     -0.78748     -0.16486
  0 0 2    331.92505   -843.20134      1.64458    -32.53866    -85.19641    181.56445
  0 2 0    303.68749   -791.79568      0.84710    -24.50108    -40.63612      3.14733
  2 0 0    -52.15963    -62.34432     -0.22926      0.91216     -4.71500      0.87207
  0 1 1    392.17671   -972.89271      2.52002    -40.26134   -102.89989    142.86106
  1 0 1     27.64436   -306.62277     -2.34347      4.01703     -2.45830      1.21836
  1 1 0     -6.31875   -250.30439     -3.12355      5.71557    -16.28739     15.41413
<::
::> States Energies (4)
  > State     Harmonic   Anharmonic     Harmonic   Anharmonic
               ZPE          ZPE    Frequency    Frequency
0 0 0   4671.50834   4602.04729            -            - 
0 0 1            -            -   3898.07520   3721.25757 
0 1 0            -            -   3841.98586   3675.50050 
1 0 0            -            -   1602.95562   1584.63509 
0 0 2            -            -   7796.15039   7419.80911 
0 2 0            -            -   7683.97173   7204.18183 
2 0 0            -            -   3205.91124   3157.70830 
0 1 1            -            -   7740.06106   7231.02596 
1 0 1            -            -   5501.03081   5291.94707 
1 1 0            -            -   5444.94148   5259.49816 
"""

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
        # analyzer.print_output_tables(print_intensities=False, print_energies=True)

"""
::> Energy Corrections
  > State    <0|dH(2)|0>  <0|dH(1)|1>  <0|dH(4)|0>  <0|dH(3)|1>  <0|dH(2)|2>  <0|dH(1)|3> 
  0 0 0    -14.81786     -0.89037     -0.43869     -0.02506     -0.02596     -0.01256
  0 0 1     -9.14842     -1.16766     -0.58993     -0.06620     -0.10732     -0.04020
  0 1 0     -8.45971      3.86545     -0.43585      0.11212     -0.16681      0.00906
  1 0 0     -3.02919     -5.27783     -0.56817     -0.14851     -0.08173     -0.06229
  0 0 2     -3.47507     -1.78907     -0.58928     -0.10337     -0.24315     -0.05697
  0 2 0     -2.09690      8.27186     -0.28334      0.39796     -0.33809      0.65282
  2 0 0      8.53771    -18.87146     -0.69620     -0.53363     -0.21704      0.00091
  0 1 1     -2.78171      3.94151     -0.28711      0.20989     -0.34169      0.27442
  1 0 1     13.96273     -5.77479     -0.10890     -0.30804     -0.34613     -0.14126
  1 1 0     16.02735      9.33511      0.01964      0.14542     -0.51555     -1.06421
<::
::> States Energies
  > State     Harmonic   Anharmonic     Harmonic   Anharmonic
               ZPE          ZPE    Frequency    Frequency
0 0 0   4671.98118   4655.77067            -            - 
0 0 1            -            -   3898.39551   3903.48630 
0 1 0            -            -   3840.97965   3852.11442 
1 0 0            -            -   1604.58720   1611.62999 
0 0 2            -            -   7796.79102   7806.74462 
0 2 0            -            -   7681.95930   7704.77413 
2 0 0            -            -   3209.17441   3213.60520 
0 1 1            -            -   7739.37516   7756.60098 
1 0 1            -            -   5502.98272   5526.47684 
1 1 0            -            -   5445.56685   5485.72512 
<::
::> Energy Corrections
  > State    <0|dH(2)|0>  <0|dH(1)|1>  <0|dH(4)|0>  <0|dH(3)|1>  <0|dH(2)|2>  <0|dH(1)|3> 
  0 0 0    -14.22611     -0.85451     -3.04382     -0.00213     -0.02637      0.00985
  0 0 1     -8.55243     -1.12596     -3.05109     -0.01464     -0.10762      0.00764
  0 1 0     -7.86609      4.19469     -8.10474      0.11033     -0.17409      0.00344
  1 0 0     -2.37918     -5.20713     -3.46249     -0.07967     -0.08368      0.00432
  0 0 2     -2.87876     -1.72562     -2.90610     -0.02325     -0.24382      0.01572
  0 2 0     -1.48771      8.91729    -13.01472      0.37676     -0.34492      0.65374
  2 0 0      9.31569    -18.96712     -3.88051     -0.42166     -0.22527      0.11450
  0 1 1     -2.17606      4.30105     -7.81164      0.24299     -0.35306      0.29907
  1 0 1     14.62549     -5.74288     -2.85819     -0.20987     -0.34703     -0.06733
  1 1 0     16.64778     10.21587     -7.93680      0.19185     -0.53549     -1.06923
<::
::> States Energies
  > State     Harmonic   Anharmonic     Harmonic   Anharmonic
               ZPE          ZPE    Frequency    Frequency
0 0 0   4671.50834   4653.36525            -            - 
0 0 1            -            -   3898.07519   3903.37418 
0 1 0            -            -   3841.98586   3848.29249 
1 0 0            -            -   1602.95562   1609.89087 
0 0 2            -            -   7796.15039   7806.53165 
0 2 0            -            -   7683.97173   7697.21526 
2 0 0            -            -   3205.91124   3209.98994 
0 1 1            -            -   7740.06106   7752.70650 
1 0 1            -            -   5501.03081   5524.57409 
1 1 0            -            -   5444.94148   5480.59855 
<::
"""