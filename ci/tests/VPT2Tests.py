
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

    @debugTest
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