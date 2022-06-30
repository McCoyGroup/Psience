
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

    @validationTest
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
        # analyzer.print_output_tables(print_intensities=False, print_energies=True)

    @debugTest
    def test_HOHCorrected(self):
        VPTRunner.run_simple(
            TestManager.test_data('HOH_freq.fchk'),
            2,
            zero_order_energy_corrections=[
                [(0, 1, 0), (4681.56364+3800) * UnitsData.convert("Wavenumbers", "Hartrees")],
                [(0, 2, 0), (4681.56364+7800) * UnitsData.convert("Wavenumbers", "Hartrees")],
                [(0, 0, 2), (4681.56364+7801) * UnitsData.convert("Wavenumbers", "Hartrees")],
            ],
            degeneracy_specs='auto'
        )

    @validationTest
    def test_OCHHInternals(self):

        tag = 'OCHH Internals'
        file_name = "OCHH_freq.fchk"

        zmatrix = [
            [0, -1, -1, -1],
            [1,  0, -1, -1],
            [2,  1,  0, -1],
            [3,  1,  0,  2]
        ]

        def conv(r, t, f, **kwargs):
            return np.array([r**2+1, t, f])
        def inv(r2, t, f, **kwargs):
            return np.array([np.sqrt(r2-1), t, f])

        # def conv(crds, **kwargs):
        #     return crds
        # def inv(crds, **kwargs):
        #     return crds

        internals = {
            'zmatrix':zmatrix,
            'conversion':conv,
            'inverse':inv,
            # 'converter_options':{
            #     'pointwise':False,
            #     # 'jacobian_prep':ZMatrixCoordinateSystem.jacobian_prep_coordinates
            # }
        }

        print("Cart:")
        # VPTAnalyzer.run_VPT(TestManager.test_data(file_name), 2).print_output_tables()
        print("ZMat:")
        # VPTAnalyzer.run_VPT(TestManager.test_data(file_name), 2, internals=zmatrix).print_output_tables()
        print("Custom:")
        # VPTRunner.run_simple(TestManager.test_data(file_name), 2, internals=internals)
        VPTAnalyzer.run_VPT(TestManager.test_data(file_name), 2, internals=internals).print_output_tables()

    @debugTest
    def test_NH3Units(self):

        tag = 'NH3 Internals'
        file_name = "nh3.fchk"

        _ = -1
        zmatrix = [
            [0, _, _, _],
            [1, 0, _, _],
            [2, 0, 1, _],
            [3, 0, 1, 2]
        ]

        def conv(r, t, f, **kwargs):
            cp1 = np.cos(f[..., 3])  # skip three for embedding
            ct1 = np.cos(t[..., 2])  # skip two for embedding
            ct2 = np.cos(t[..., 3])
            st1 = np.sin(t[..., 2])
            st2 = np.sin(t[..., 3])
            f[..., 3] = np.arccos(st1*st2*cp1+ct1*ct2)
            return np.array([r, t, f])
        def inv(r, t, f, **kwargs):
            cp1 = np.cos(f[..., 3])
            ct1 = np.cos(t[..., 2])
            ct2 = np.cos(t[..., 3])
            st1 = np.sin(t[..., 2])
            st2 = np.sin(t[..., 3])
            f[..., 3] = np.arccos((cp1-ct1*ct2)/(st1*st2))
            return np.array([r, t, f])

        # def conv(crds, **kwargs):
        #     return crds
        # def inv(crds, **kwargs):
        #     return crds

        # internals = {
        #     'zmatrix':zmatrix,
        #     'conversion':conv,
        #     'inverse':inv,
        #     # 'converter_options':{
        #     #     'pointwise':False,
        #     #     # 'jacobian_prep':ZMatrixCoordinateSystem.jacobian_prep_coordinates
        #     # }
        # }

        file = TestManager.test_data(file_name)
        print("Cart:")
        VPTAnalyzer.run_VPT(file, 2).print_output_tables()
        print("ZMat:")
        # VPTAnalyzer.run_VPT(TestManager.test_data(file_name), 2, internals=zmatrix).print_output_tables()
        print("Custom:")
        # VPTRunner.run_simple(TestManager.test_data(file_name), 2, internals=internals)
        # VPTAnalyzer.run_VPT(TestManager.test_data(file_name), 2, internals=internals, handle_strong_couplings=True).print_output_tables()

        """ With Cartesian coordinates
                         Harmonic                  Anharmonic
State             Frequency    Intensity       Frequency    Intensity
  0 0 0 0 0 1    3649.84753      8.40063      3673.97101      8.16596
  0 0 0 0 1 0    3649.84749      8.40063      3673.99153      8.16597
  0 0 0 1 0 0    3502.88652      3.39049      3504.89231      4.29459
  0 0 1 0 0 0    1668.90366     14.31874      1611.87387     15.15334
  0 1 0 0 0 0    1668.90366     14.31874      1611.87470     15.15358
  1 0 0 0 0 0    1037.51781    139.18086       907.20372    146.77249
  0 0 0 0 0 2    7299.69506      0.00000      7358.46413      0.23938
  0 0 0 0 2 0    7299.69499      0.00000      7322.01149      0.00313
  0 0 0 2 0 0    7005.77304      0.00000      6948.64326      0.01480
  0 0 2 0 0 0    3337.80732      0.00000      3216.64730      0.29740
  0 2 0 0 0 0    3337.80733      0.00000      3191.38576      0.04236
  2 0 0 0 0 0    2075.03561      0.00000      1716.91900      0.05990
  0 0 0 0 1 1    7299.69502      0.00000      7541.05092      0.25778
  0 0 0 1 0 1    7152.73405      0.00000      7166.06224      0.21966
  0 0 1 0 0 1    5318.75119      0.00000      5240.10874      0.00432
  0 1 0 0 0 1    5318.75119      0.00000      5279.00970      1.21013
  1 0 0 0 0 1    4687.36534      0.00000      4547.87059      0.68588
  0 0 0 1 1 0    7152.73402      0.00000      7202.90789      0.20904
  0 0 1 0 1 0    5318.75115      0.00000      5274.40898      0.00033
  0 1 0 0 1 0    5318.75116      0.00000      5235.65311      1.20013
  1 0 0 0 1 0    4687.36530      0.00000      4547.88276      0.68588
  0 0 1 1 0 0    5171.79018      0.00000      5099.83410      0.04409
  0 1 0 1 0 0    5171.79019      0.00000      5099.83542      0.04422
  1 0 0 1 0 0    4540.40433      0.00000      4425.74209      0.84746
  0 1 1 0 0 0    3337.80732      0.00000      3216.43526      0.29740
  1 0 1 0 0 0    2706.42147      0.00000      2512.73682      0.00177
  1 1 0 0 0 0    2706.42147      0.00000      2512.73850      0.00177
  0 2 0 1 0 0    6840.69385      0.00000      6713.79845      0.00011
  0 0 2 1 0 0    6840.69384      0.00000      6698.61681      0.00040
  0 3 0 0 0 0    5006.71099      0.00000      4789.39177      0.00714
  0 0 3 0 0 0    5006.71098      0.00000      4789.38260      0.00626
                """
        """ With regular Z-matrix coordinates
============================================= IR Data ==============================================
                         Harmonic                  Anharmonic
State             Frequency    Intensity       Frequency    Intensity
  0 0 0 0 0 1    3649.84753      8.40063      3726.14107      8.01000
  0 0 0 0 1 0    3649.84749      8.40063      3720.43203      8.14962
  0 0 0 1 0 0    3502.88652      3.39049      3504.98059      4.27590
  0 0 1 0 0 0    1668.90366     14.31874      1635.12777     15.11316
  0 1 0 0 0 0    1668.90367     14.31874      1634.58245     15.13005
  1 0 0 0 0 0    1037.51781    139.18086       951.20470    148.04683
  0 0 0 0 0 2    7299.69506      0.00000      7443.94606      0.24253
  0 0 0 0 2 0    7299.69499      0.00000      7420.05073      0.00705
  0 0 0 2 0 0    7005.77304      0.00000      6958.79881      0.01981
  0 0 2 0 0 0    3337.80732      0.00000      3263.93164      0.30176
  0 2 0 0 0 0    3337.80733      0.00000      3234.19765      0.04307
  2 0 0 0 0 0    2075.03562      0.00000      1804.93029      0.06297
  0 0 0 0 1 1    7299.69502      0.00000      7636.09541      0.26285
  0 0 0 1 0 1    7152.73405      0.00000      7227.52952      0.22936
  0 0 1 0 0 1    5318.75119      0.00000      5338.82560      0.28104
  0 1 0 0 0 1    5318.75120      0.00000      5378.12318      1.09984
  1 0 0 0 0 1    4687.36534      0.00000      4697.63552      0.70847
  0 0 0 1 1 0    7152.73402      0.00000      7257.29782      0.19690
  0 0 1 0 1 0    5318.75115      0.00000      5374.55183      0.13318
  0 1 0 0 1 0    5318.75116      0.00000      5333.95854      0.94634
  1 0 0 0 1 0    4687.36530      0.00000      4675.38562      0.70511
  0 0 1 1 0 0    5171.79018      0.00000      5128.93044      0.04608
  0 1 0 1 0 0    5171.79019      0.00000      5129.48827      0.04659
  1 0 0 1 0 0    4540.40433      0.00000      4469.82653      0.85590
  0 1 1 0 0 0    3337.80732      0.00000      3257.94962      0.30112
  1 0 1 0 0 0    2706.42147      0.00000      2577.72518      0.00182
  1 1 0 0 0 0    2706.42147      0.00000      2579.04652      0.00182
  0 3 0 0 0 0    5006.71100      0.00000      4848.28115      0.00859
  0 0 3 0 0 0    5006.71098      0.00000      4848.52584      0.00458
        """
        """ With symmetrized coordinates
State             Frequency    Intensity       Frequency    Intensity
  0 0 0 0 0 1    3649.84753      8.40063      3723.70074      8.02349
  0 0 0 0 1 0    3649.84749      8.40063      3723.72163      8.02350
  0 0 0 1 0 0    3502.88652      3.39049      3504.97874      4.27470
  0 0 1 0 0 0    1668.90366     14.31874      1634.65206     15.13159
  0 1 0 0 0 0    1668.90367     14.31874      1634.64634     15.13156
  1 0 0 0 0 0    1037.51781    139.18086       952.40922    148.07587
  0 0 0 0 0 2    7299.69506      0.00000      7444.16390      0.24658
  0 0 0 0 2 0    7299.69499      0.00000      7421.28673      0.00317
  0 0 0 2 0 0    7005.77305      0.00000      6958.79513      0.01981
  0 0 2 0 0 0    3337.80732      0.00000      3263.57753      0.30173
  0 2 0 0 0 0    3337.80733      0.00000      3233.82922      0.04292
  2 0 0 0 0 0    2075.03561      0.00000      1807.33816      0.06305
  0 0 0 0 1 1    7299.69503      0.00000      7637.00847      0.26180
  0 0 0 1 0 1    7152.73405      0.00000      7229.66109      0.21733
  0 0 1 0 0 1    5318.75119      0.00000      5339.24402      0.00449
  0 1 0 0 0 1    5318.75120      0.00000      5376.53040      1.23248
  1 0 0 0 0 1    4687.36534      0.00000      4689.39755      0.70723
  0 0 0 1 1 0    7152.73402      0.00000      7256.19539      0.20990
  0 0 1 0 1 0    5318.75115      0.00000      5373.43275      0.00026
  0 1 0 0 1 0    5318.75116      0.00000      5336.29143      1.22319
  1 0 0 0 1 0    4687.36530      0.00000      4689.40994      0.70723
  0 0 1 1 0 0    5171.79018      0.00000      5129.06995      0.04504
  0 1 0 1 0 0    5171.79019      0.00000      5129.06553      0.04502
  1 0 0 1 0 0    4540.40433      0.00000      4471.02511      0.85613
  0 1 1 0 0 0    3337.80733      0.00000      3257.49832      0.30121
  1 0 1 0 0 0    2706.42147      0.00000      2579.33175      0.00182
  1 1 0 0 0 0    2706.42147      0.00000      2579.32464      0.00182
  0 3 0 0 0 0    5006.71100      0.00000      4847.85818      0.00744
  0 0 3 0 0 0    5006.71098      0.00000      4847.85797      0.00578
  """