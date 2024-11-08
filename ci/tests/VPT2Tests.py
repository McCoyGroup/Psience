import itertools

import McUtils.Numputils

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
    @classmethod
    def setUpClass(cls) -> None:
        np.set_printoptions(linewidth=int(1e8))

    # region Water Analogs

    @validationTest
    def test_MultdiDegHOH(self):

        # space = VPTStateSpace(
        #     [
        #         [0, 1, 0],
        #         [0, 0, 1]
        #     ],
        #     degeneracy_specs=[
        #         {
        #             'polyads': [
        #                 [[0, 0, 1], [1, 1, 0]],
        #             ]
        #         },
        #         {
        #             'polyads': [
        #                 [[0, 1, 0], [1, 0, 0]]
        #             ]
        #         }
        #     ]
        # )
        #
        # space = VPTRunner.run_simple(
        #     TestManager.test_data('HOH_freq.fchk'),
        #     1,
        #     degeneracy_specs=[
        #          {
        #             'polyads':[
        #                 [[0, 0, 1], [1, 1, 0]],
        #             ]
        #         },
        #         {
        #             'polyads': [
        #                 [[0, 1, 0], [1, 0, 0]]
        #             ]
        #         }
        #     ]
        # )

        space = VPTRunner.run_simple(
            TestManager.test_data('HOH_freq.fchk'),
            2,
            degeneracy_specs=[
                {
                    'polyads': [
                        [[0, 0, 1], [1, 1, 0]],
                    ]
                },
                {
                    'wfc_threshold': .3
                }
            ]
        )

        raise Exception(...)

    @validationTest
    def test_HOHAnalytic(self):

        file_name = "HOH_freq.fchk"
        runner, states = AnalyticVPTRunner.construct(
            TestManager.test_data(file_name),
            [[
                0,
                [[0, 0, 1]]
            ]],
            # expressions_file=os.path.expanduser("~/Desktop/exprs.hdf5")
        )
        classic, _ = runner.construct_classic_runner(states)
        classic.print_tables()
        # corrs_a = runner.get_transition_moment_corrections(
        #     states,
        #     # axes=[2],
        #     # dipole_expansion=dips,
        #     # terms=[(0, 0, 0)]
        # )  # .corrections[0]
        # print(corrs_a)
        runner.run_VPT(states)

    @validationTest
    def test_AnalyticWFC(self):

        file_name = "OCHH_freq.fchk"
        AnalyticVPTRunner.run_simple(
            TestManager.test_data(file_name),
            [
                [
                    0,
                    [
                        [0, 0, 0, 0, 0, 1],
                        [0, 1, 0, 1, 0, 0],
                        [0, 0, 0, 1, 1, 0],
                        [0, 0, 0, 0, 1, 0],
                    ],
                ],
                [
                    [0, 0, 0, 0, 1, 0],
                    [
                        [0, 0, 0, 0, 1, 1],
                        [0, 1, 0, 1, 1, 0]
                    ]
                ]
            ],
            expressions_file=os.path.expanduser("~/Desktop/exprs.hdf5"),
            degeneracy_specs='auto',
            handle_degeneracies=True
        )

    @inactiveTest
    def test_PartialRebuild(self):
        state = VPTStateMaker(7)
        corrs = AnalyticVPTRunner.run_simple(
            mol,
            # [
            #     state(),
            #     state(21),
            #     state(20),
            #     state(19),
            #     state([21, 2]),
            #     state([20, 2]),
            #     state([19, 2]),
            #     state(21, 19),
            #     state(21, 20),
            #     state(20, 19)
            # ],
            [
                [
                    [state()],
                    [
                        state(1),
                        state(2),
                        state(3),
                        state([1, 2]),
                        state([2, 2]),
                        state([3, 2]),
                        state(1, 2),
                        state(1, 3),
                        state(2, 3),
                    ]
                ]
            ],
            full_surface_mode_selection=[108 - 108, 108 - 107, 108 - 106, 108 - 105, 108 - 21, 108 - 20, 108 - 19,
                                         108 - 1],
            # degeneracy_specs = [
            #   {'polyads':[[state(19), state(20, 107)]]},
            #   {'polyads':[[state(106, 106), state(108, 108)]]},
            #   {'polyads':[[state(106), state(107)], [state(105), state(107)], [state(105), state(108)]]}
            # ],
            mode_selection=[108 - 108, 108 - 107, 108 - 106, 108 - 105, 108 - 21, 108 - 20, 108 - 19],
            # degeneracy_specs = [
            #   {'polyads':[[state(1), state(2, 6)]]},
            #   {'polyads':[[state([5, 2]), state([7, 2])]]},
            #   {
            #       'polyads':
            #         [
            #           [state(4), state(6)],
            #           [state(4), state(7)],
            #           [state(5), state(6)]
            #         ]
            #     }
            # ],
            # logger=output_file,
            expressions_file=os.path.expanduser("exprs.hdf5")
        )

    @debugTest
    def test_AnalyticOCHHMultiple(self):

        """

::       1(1) | 2908.616 | 90.480 | 2689.933 | 68.349 | 2732.141 | 105.713
::       2(1) | 2828.758 | 65.833 | 2660.651 | 60.476 | 2660.651 |  60.476
::       3(2) | 3281.457 |  0.000 | 3193.397 |  5.802 | 3193.397 |   5.802
::   3(1)5(1) | 2830.284 |  0.000 | 2811.622 | 37.349 | 2769.415 |   0.001
        :return:
        """

        file_name = "OCHH_freq.fchk"
        runner, states = AnalyticVPTRunner.construct(
            TestManager.test_data(file_name),
            [
                [
                    0,
                    [
                        [0, 0, 0, 0, 0, 1],
                        [0, 0, 0, 0, 1, 0],
                        [0, 0, 0, 2, 0, 0]
                    ],
                ],
                [
                    [
                        [0, 0, 0, 0, 0, 1],
                    ],
                    [
                        [0, 0, 0, 0, 1, 0],
                        [0, 0, 0, 2, 0, 0]
                    ],
                ]
            ],
            degeneracy_specs='auto'
            # degeneracy_specs={
            #     'polyads': [
            #         [
            #             [0, 0, 0, 0, 0, 1],
            #             [0, 0, 0, 0, 1, 0]
            #         ],
            #         [
            #             [0, 0, 0, 0, 0, 1],
            #             [0, 1, 0, 1, 0, 0]
            #         ]
            #     ]
            # }
        )
        classic, _ = runner.construct_classic_runner(states,
                                                     zero_element_warning=False
                                                     )
        classic.print_tables()
        runner.run_VPT(states)
        raise Exception(...)

    @validationTest
    def test_AnalyticOCHH(self):

        file_name = "OCHH_freq.fchk"
        AnalyticVPTRunner.run_simple(
            TestManager.test_data(file_name),
            [
                [
                    0,
                    [
                        [0, 0, 0, 0, 0, 1],
                        [0, 1, 0, 1, 0, 0],
                        [0, 0, 0, 1, 1, 0],
                        [0, 0, 0, 0, 1, 0],
                    ],
                ],
                [
                    [0, 0, 0, 0, 1, 0],
                    [
                        [0, 0, 0, 0, 1, 1],
                        [0, 1, 0, 1, 1, 0]
                    ]
                ]
            ],
            # expressions_file=os.path.expanduser("~/Desktop/exprs.hdf5"),
            degeneracy_specs=[
                [[0, 0, 0, 0, 0, 1], [0, 1, 0, 1, 0, 0]]
            ]
        )

    @validationTest
    def test_AnalyticHOONO(self):

        file_name = "HOONO_freq.fchk"
        state = VPTStateMaker(9)
        with BlockProfiler():
            AnalyticVPTRunner.run_simple(
                TestManager.test_data(file_name),
                2,
                degeneracy_specs=[
                    [state(1), state(8, 7)]
                ],
                expressions_file=os.path.expanduser("~/Documents/Postdoc/exprs.hdf5")
            )
        raise Exception(...)

        # VPTRunner.run_simple(
        #     TestManager.test_data(file_name),
        #     2,
        #     # mode_selection=mode_selection,
        #     # internals=internals,
        #     logger=True,
        #     mixed_derivative_handling_mode='averaged'
        # )

        runner, states = AnalyticVPTRunner.construct(
            TestManager.test_data(file_name),
            2,
            # mode_selection=[8, 7, 6],
            # internals=internals,
            logger=True,
            # expansion_order = {
            #     'coriolis':0,
            #     'pseudopotential':0
            # },
            # return_runner=True,
            # return_states=True,
            # calculate_intensities=False
        )

        # raise Exception(
        #     runner.get_freqs(states)[1] * UnitsData.convert("Hartrees", "Wavenumbers")
        # )

        # runner.eval.expansions[2][0][:] = 0
        # v3 = runner.eval.expansions[1][0]
        # for i in itertools.combinations(range(3), 1):
        #     for p in itertools.permutations([i, i, i]):
        #         v3[p] = 0#1e-5
        # for i,j in itertools.combinations(range(3), 2):
        #     for p in itertools.permutations([i, i, j]):
        #         v3[p] = 0#1e-5
        #     for p in itertools.permutations([i, j, j]):
        #         v3[p] = 0
        # for i,j,k in itertools.combinations(range(3), 3):
        #     for p in itertools.permutations([i, j, k]):
        #         v3[p] = 1e-4

        og, _ = runner.construct_classic_runner(
            TestManager.test_data(file_name),
            states,
            mode_selection=np.arange(len(states[0])),
            zero_element_warning=False,
            logger=False
            # internals=internals,
            # degeneracy_specs='auto'
            # dipole_terms=dips
        )


        # from Peeves import Timer
        # with Timer():
        # og.print_tables(print_intensities=True)
        # print(freqs[0], freqs[1:]-freqs[0])
        # with BlockProfiler():
        # with Timer():
        # with BlockProfiler():
        #     spec = runner.get_spectrum(states, verbose=False)
        # print(np.array(spec).T)

    @validationTest
    def test_TrimerAnalytic(self):

        file_name = "water_trimer_freq.fchk"
        # VPTRunner.run_simple(
        #     TestManager.test_data(file_name),
        #     2,
        #     # mode_selection=mode_selection,
        #     # internals=internals,
        #     logger=True,
        #     mixed_derivative_handling_mode='averaged'
        # )

        runner, states = AnalyticVPTRunner.construct(
            TestManager.test_data(file_name),
            2,
            logger=True
        )

        og, _ = runner.construct_classic_runner(
            TestManager.test_data(file_name),
            states,
            mode_selection=np.arange(len(states[0])),
            logger=False
        )

        # og.print_tables(print_intensities=True)
        with BlockProfiler(print_options={'show_all':True}):
            spec = runner.get_spectrum(states, verbose=False)
        print(np.array(spec).T)

    @validationTest
    def test_SelAnharmAnalytic(self):

        file_name ="freq_anion.fchk"

        state = VPTStateMaker(108)
        # with BlockProfiler():
        corrs = AnalyticVPTRunner.run_simple(
            file_name,
            [
                state(),
                state(21),
                state(20),
                state(19),
            ],
            full_surface_mode_selection=[108-21, 108-20, 108-19, 108-1],
            logger=True,
            expressions_file=os.path.expanduser("exprs.hdf5")
        )

        print(corrs)

        raise Exception(...)

    @validationTest
    def test_PyreneAnalytic(self):
        file_name = "pyrene_try3.fchk"
        state = VPTStateMaker(72)
        degs = [
            [state(8), state(9)],
            [state(15), state(69)],
            [state(16), state(39)],
            [state(19), state(20)],
            [state(19), state(15)],
            [state(20), state(16)],
            [state(21), state(70)],
            [state(22), state(18)],
            [state(23), state(24)],
            [state(24), state(55)],
            [state(27), state(28)],
            [state(35), state(66)],
            [state(36), state(37)],
            [state(37), state(15)],
            [state(38), state(69)],
            [state(39), state(40)],
            [state(41), state(71)],
            [state(42), state(43)],
            [state(67), state(68)],
            [state(67), state(37)],
            [state(68), state(69)],
            [state(71), state(18)],
            [state(71), state(22)]
        ]
        # with BlockProfiler():
        corrs = AnalyticVPTRunner.run_simple(
            file_name,
            # [
            #     state(8), state(9),
            #     state(10), state(11)
            # ],
            1,
            expressions_file="exprs.hdf5",
            degeneracy_specs=degs,
            logger=True,
            zero_cutoff=1e-12
        )

    @validationTest
    def test_NewEmbedding(self):
        file_name = "HOH_freq.fchk"
        test_internals = [[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]]
        # file_name = "water_dimer_freq.fchk"
        # test_internals = dimer_internals
        # file_name = "HOONO_freq.fchk"
        # test_internals = hoono_internals
        # file_name = "OCHH_freq.fchk"

        mol = Molecule.from_file(
            TestManager.test_data(file_name),
            internals=test_internals
        )#.get_embedded_molecule()
        h1 = mol.hamiltonian

        runner1, _ = VPTRunner.construct(
            mol,
            2,
            internals=test_internals,
            logger=False
        )
        h2 = runner1.hamiltonian

        t1 = h1.gmatrix_expansion()
        t2 = h2.G_terms

        # print(t1.embedding.embedding.coords)
        # print(t2.base_terms.molecule.coords)
        # raise Exception(...)

        # print(...)
        # # print(
        # #     t1.embedding.get_internals_by_mw_cartesians(order=2, strip_embedding=True)[0]
        # # )
        # # print(
        # #     t2.base_terms.get_internals_by_cartesians(order=2)[0]
        # # )
        # print(
        #     np.round(
        #         t1.embedding.get_mw_cartesians_by_internals(order=2, strip_embedding=True)[0],
        #         8
        #     )
        # )
        # print("="*50)
        # print(
        #     np.round(
        #         t2.base_terms.get_cartesians_by_internals(order=1)[0],
        #         8
        #     )
        # )

        # print(
        #     np.round(
        #         t1.embedding.get_internals_by_mw_cartesians(order=2, strip_embedding=True)[0],
        #         8
        #     )
        # )
        # print(
        #     np.round(
        #         t2.base_terms.get_internals_by_cartesians(order=1)[0],
        #         8
        #     )
        # )

        # print(
        #     np.round(
        #         t1.embedding.get_internals_by_cartesians(order=2, strip_embedding=True)[1][0],
        #         8
        #     )
        # )
        # print(
        #     np.round(
        #         t2.base_terms.get_modes_by_cartesians(order=2)[1][0],
        #         8
        #     )
        # )
        # raise Exception(...)
        # QY_derivs = self.get_modes_by_cartesians(order=order + 1)
        # YQ_derivs = self.get_cartesians_by_modes(order=order + 1)
        # g1 = h1.get_VPT_expansions({'kinetic':2})['kinetic']
        # g2 = h2.G_terms.base_terms
        # t2.base_terms.gmatrix_tolerance = None
        # t2.base_terms.freq_tolerance = None

        # print(np.round(g1[1][0] - g2[1][0], 8))
        # print("_"*20)
        # print(np.round(g1[2][0, 0], 8))
        # print("_"*20)
        # print(np.round(g2[2][0, 0], 8))
        # # print(np.round(g2[1][0], 8))
        # # g2.get_terms(2)
        # # print(g1[2] - g2[2])
        #
        # raise Exception(...)

        # G = g1
        # VPTRunner.construct(
        #     mol,
        #     2,
        #     internals=test_internals,
        #     kinetic_terms=G,
        #     # potential_terms=[G[0], 0, 0],
        #     include_pseudopotential=False,
        #     include_coriolis_coupling=False,
        #     logger=False
        # )[0].print_tables(print_intensities=False)
        # VPTRunner.construct(
        #     mol,
        #     2,
        #     internals=test_internals,
        #     kinetic_terms=[h2.G_terms[0], h2.G_terms[1], h2.G_terms[2]],
        #     # potential_terms=[G[0], 0, 0],
        #     include_pseudopotential=False,
        #     include_coriolis_coupling=False,
        #     logger=False
        # )[0].print_tables(print_intensities=False)

        import itertools
        from McUtils.Zachary import TensorDerivativeConverter
        # def symm(a, k):
        #     n = 0
        #     t = 0
        #     for p in itertools.permutations(range(k)):
        #         n += 1
        #         t += a.transpose(p + tuple(range(k, a.ndim)))
        #     return t / n
        #
        # u1 = h1.get_VPT_expansions({'potential': 3})['potential']
        # print(u1[-1])
        # u2 = h2.V_terms[1]
        # print(u2)

        gg = h1.gmatrix_expansion(dimensionless=True)
        ke, ki = h1.get_kinetic_optmizing_transformation(2, dimensionless=True)
        woof = gg.get_terms(2, transformation=(ke, ki))
        print(woof[0])
        print(woof[1])

        # ki = nput.inverse_transformation(ke, 2)
        # for t in nput.tensor_reexpand(ki, ke, 2):
        #     print(np.round(t, 8))
        # for t in nput.tensor_reexpand(ke, ki, 2):
        #     print(np.round(t, 8))

        # woof = gg.reexpress(2, ke, ki)
        # print(woof[0])
        # print(woof[1])

    @inactiveTest
    def test_HOHNoKE(self):
        print()

        # VPTRunner.run_simple(
        #     TestManager.test_data(file_name),
        #     2,
        #     # include_pseudopotential=False,
        #     # include_coriolis_coupling=False,
        #     # kinetic_terms=[np.zeros((3,)*(i+2)) for i in range(3)],
        #     potential_terms=[np.zeros((3,)*(i+2)) for i in range(3)],
        #     zero_element_warning=False
        # )



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

        dimer_internals = [
            [LHF, X, X, X],
            [LO, LHF, X, X],
            [SH, LO, LHF, X],
            [RH2, SH, LO, LHF],  # get out of plane
            [RO, LO, RH2, LHF],
            [RH1, RO, RH2, LHF]
        ]

        hoono_internals = [
            [1, -1, -1, -1],
            [2,  1, -1, -1],
            [3,  2,  1, -1],
            [0,  1,  2,  3],
            [4,  3,  2,  1]
        ]


        test_internals = None
        file_name = "HOH_freq.fchk"
        # test_internals = [[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]]
        # file_name = "water_dimer_freq.fchk"
        # test_internals = dimer_internals
        # file_name = "HOONO_freq.fchk"
        # test_internals = hoono_internals
        # file_name = "OCHH_freq.fchk"
        # file_name = "nh3.fchk"
        # test_internals = [[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1], [3, 0, 1, 2]]

        runner1, _ = VPTRunner.construct(
            TestManager.test_data(file_name),
            2,
            internals=test_internals,
            # expansion_order={'potential':0},
            # include_pseudopotential=False,
            logger=False
        )
        # runner1.print_tables()
        # R2 = np.zeros((3, 3, 3))
        # for i,j,k in itertools.permutations()
        # R = [np.eye(3), ]
        # Q = [np.random.rand(3, 3), np.random.rand(3, 3, 3), np.random.rand(3, 3, 3, 3)]

        G_new = runner1.hamiltonian.G_terms.base_terms.optimize_coordinates()
        raise Exception([G.shape for G in G_new])


        runner2, _ = VPTRunner.construct(
            TestManager.test_data(file_name),
            2,
            internals=[[1, -1, -1, -1], [2, 1, -1, -1], [0, 1, 2, -1]],
            logger=False
        )

        # print(runner1.hamiltonian.G_terms[1]+2/3*runner1.hamiltonian.V_terms[1])
        # print(runner2.hamiltonian.G_terms[1]+2/3*runner2.hamiltonian.V_terms[1])
        #
        # raise Exception(...)

        from McUtils.Zachary import TensorDerivativeConverter

        def symm_terms(v_terms):
            test_V = [0, 0, 0, 0]
            for i in [3, 2, 1]:
                test_V[i] = v_terms[i-1]
            test_V[0] = np.zeros(len(test_V[1]))
            n = len(test_V[0])
            for i in range(n):
                for j in range(n):
                    for k in range(n):
                        perms = list(itertools.permutations([i, j, k]))
                        v3 = sum(test_V[2][p] for p in perms) / len(perms)
                        for p in perms:
                            test_V[2][p] = v3
                        for l in range(n):
                            perms = list(itertools.permutations([i, j, k, l]))
                            v4 = sum(test_V[3][p] for p in perms) / len(perms)
                            for p in perms:
                                test_V[3][p] = v4

            return test_V

        test_V1 = symm_terms(runner1.hamiltonian.V_terms)
        n = len(test_V1[0])
        test_G = [
            np.zeros((n,) * (i + 2)) if isinstance(g, int) and g == 0 else g
            for i, g in enumerate(reversed([runner1.hamiltonian.G_terms[n] for n in [2, 1, 0]]))
        ]
        # test_V2 = symm_terms(runner2.hamiltonian.V_terms)

        V3 = test_V1[2]
        freqs = np.diag(test_V1[1])

        V3 = np.tensordot(V3, np.tensordot(V3, np.diag(1/freqs), axes=[2, -1]), axes=[2, -1])
        V3 = V3 + np.moveaxis(V3, 2, 1) + np.moveaxis(V3, 2, 0)

        ders = [
            np.eye(len(test_V1[1])),
            -test_V1[2] / (3*freqs[np.newaxis, np.newaxis, :]),
            -(test_V1[3] - 5/9*V3) / (4*freqs[np.newaxis, np.newaxis, np.newaxis, :]),
            0
        ]
        # test_V[2] = 0
        # test_V[3] = 0
        # corr = wtf + np.moveaxis(wtf, 2, 0) + np.moveaxis(wtf, 1, 0)
        # print(corr[0])# + test_V[2])
        # raise Exception(...)
        # convs = TensorDerivativeConverter(ders, test_V1).convert()#print_transformations=True)
        convs = TensorDerivativeConverter.convert_fast(ders, test_V1, order=len(ders))
        # print(convs[2])
        # print(convs[3])

        R2 = test_V1[2] / (3*freqs[np.newaxis, np.newaxis, :])
        # R3 = (test_V1[3] - 6/9*V3) / (4*freqs[np.newaxis, np.newaxis, np.newaxis, :])
        # R2Q2 = np.tensordot(R2, ders[1], axes=[2, 0])
        # R2Q2 = R2Q2 + R2Q2.transpose(0,2,1,3) + R2Q2.transpose(2,1,0,3)
        # R3 = -(ders[2] + R2Q2)
        R3 = (test_V1[3] - 1/9*V3) / (4*freqs[np.newaxis, np.newaxis, np.newaxis, :])

        dersR = [
            np.eye(len(test_V1[1])),
            R2,
            R3,
            0
        ]
        # convs_R = TensorDerivativeConverter(dersR, convs).convert()  # print_transformations=True)

        # print(np.round(convs[2], 7))
        # print(np.round(convs[3], 7))
        # print(np.round(convs_R[1] - test_V1[1], 7))
        # print(np.round(convs_R[2] - test_V1[2], 7))
        # print(np.round(convs_R[3] - test_V1[3], 7))
        # raise Exception(...)

        G0 = np.diag(freqs)

        Q1 = ders[0]
        Q2 = ders[1]
        Q3 = ders[2]

        R1 = dersR[0]

        G1_Q = np.tensordot(
            R2,
            G0,
            axes=[1, 0]
        )
        # G1 = test_G[1] + (G1_Q + G1_Q.transpose(0, 2, 1))

        W = freqs[:, np.newaxis] / freqs[np.newaxis, :]
        V1 = test_V1[2]
        VV = test_V1[2][:, :, :, np.newaxis, np.newaxis] * test_V1[2][np.newaxis, np.newaxis, :, :, :]
        G1 = test_G[1] + test_V1[2] / 3 * (W + W.T)[np.newaxis]

        G2_QR = 2 * np.tensordot(
            ders[1],
            test_V1[2]/3,
            axes=[2,0]
        )

        # G2_QR_exp = -sum(
        #     (
        #         VV[:, :, a, :, :] / 3 * W[np.newaxis, np.newaxis, :, :]
        #         #+ V1[:, :, a, np.newaxis, np.newaxis] * test_G[1][np.newaxis, np.newaxis, a, :, :]
        #     ) / (3 * freqs[a])
        #     for a in range(3)
        # )
        #
        # print(G2_QR - G2_QR_exp)
        #
        # raise Exception(...)

        G2_RR = np.tensordot(
            R2,
            R2 * freqs[np.newaxis, :, np.newaxis],
            axes=[1,1]
        )
        G2_RR = sum(
            G2_RR.transpose(*p)
            for p in [
                [0, 2, 1, 3],
                [0, 2, 3, 1]
            ]
        )

        G2_R3 = R3 * freqs[np.newaxis, np.newaxis, :, np.newaxis]
        G2_R3 = G2_R3 + G2_R3.transpose(0, 1, 3, 2)

        G2_RG = np.tensordot(
            R2,
            test_G[1],
            axes=[1,1]
        )

        G2_RG = sum(
            G2_RG.transpose(*p)
            for p in [
                [0,2,1,3],
                [2,0,1,3],
                [0,2,3,1],
                [2,0,3,1]
            ]
        )

        G2_QG = np.tensordot(
            Q2,
            test_G[1],
            axes=[2,0]
        )

        G2 = (
            G2_QR
            + G2_RR
            + G2_RG
            + G2_QG
            + G2_R3
            + test_G[2]
        )


        # print((G1_Q[0] + G1_Q.transpose(0, 2, 1))[0])
        # print(G1_Q[0])
        # print(test_G[1][0])
        # print(R2[0])
        # print(G1[0] - test_G[1][0])

        # print(
        #     (
        #             + G2_RR
        #             # + G2_RG
        #             + G2_R3
        #             + test_G[2]
        #     )[0, 0]
        # )

        ggg = [G0, G1, G2]
        ggg_test = runner1.hamiltonian.reexpress_G([Q1, Q2, Q3], [R1, R2, R3])

        # print(ggg[1] - ggg_test[1])
        # print(ggg[1][0])
        # print(ggg_test[1][0])
        # print(ggg[2][0, 0])
        # print(ggg_test[2][0, 0])
        # runner1.print_tables()
        # R = [np.random.rand(3, 3), np.random.rand(3, 3, 3), np.random.rand(3, 3, 3, 3)]
        # Q = [np.random.rand(3, 3), np.random.rand(3, 3, 3), np.random.rand(3, 3, 3, 3)]
        #
        # G_new = runner1.hamiltonian.reexpress_G(R, Q)
        # raise Exception([G.shape for G in G_new])

        (Q, R), _ = runner1.hamiltonian.get_potential_optimized_coordinates()
        # print(Q[1][0])
        # print(Q2[0])
        # print(R[1] - R2)
        # print(Q[2][:, :, 0, 0])
        # print(Q3[:, :, 0, 0])

        # raise Exception(...)

        # print(np.sum(np.abs(G2 - np.moveaxis(G2, 1, 0))))
        # print(np.sum(np.abs(G2 - np.moveaxis(G2, 3, 2))))
        # # print(np.max(G2), np.min(G2), np.max(test_G[2]), np.min(test_G[2]))
        # raise Exception(...)


        # runner1.print_tables(print_intensities=False)
        rrr, _ = VPTRunner.construct(
            TestManager.test_data(file_name),
            2,
            internals=test_internals,
            kinetic_terms=[G0, test_G[1], test_G[2]],
            potential_terms=[G0, test_V1[2], test_V1[3]],
            include_pseudopotential=False,
            include_coriolis_coupling=True,
            logger=False
        )
        rrr.print_tables(print_intensities=False)
        """
        :: State    <0|dH(2)|0>  <0|dH(1)|1> 
              0 0 0    100.74075   -156.70145
              0 0 1    296.33507   -545.08549
              0 1 0    301.26122   -538.53511
              1 0 0    108.13556   -213.67502
              0 0 2    587.47448  -1127.09543
              0 2 0    597.90621  -1104.59990
              2 0 0    109.32398   -292.47319
              0 1 1    687.69956  -1284.13087
              1 0 1    312.98434   -634.36364
              1 1 0    326.60393   -633.50201
        """
        """
              0 0 0    105.49844   -157.99275
              0 0 1    335.26400   -580.55272
              0 1 0    279.87458   -513.67920
              1 0 0    105.42497   -207.50600
              0 0 2    703.29509  -1239.45907
              0 2 0    554.24618  -1057.46833
              2 0 0    108.34454   -288.04967
              0 1 1    678.58740  -1271.55370
              1 0 1    349.24337   -667.17858
              1 1 0    259.61523   -563.04545
              """

        """
            0 0 0     10.85263   -144.38370
            0 0 1    103.73295   -428.81417
            0 1 0    102.30966   -416.98116
            1 0 0    -35.81457   -148.46241
            0 0 2    249.91709   -862.39340
            0 2 0    252.94485   -836.86328
            2 0 0    -83.05123   -180.01127
            0 1 1    313.14314   -985.39758
            1 0 1     18.26562   -419.47037
            1 1 0      2.24708   -387.69952
            """

        # for i in range(n):
        #     for j in range(n):
        #         for k in range(n):
        #             if not (i==j and j==k):
        #                 G1[i,j,k] = 0

        rrr, _ = VPTRunner.construct(
            TestManager.test_data(file_name),
            2,
            internals=test_internals,
            kinetic_terms=ggg_test,#[G0, G1, G2],
            potential_terms=[G0, 0, 0],
            include_pseudopotential=False,
            include_coriolis_coupling=True,
            logger=False,
            zero_element_warning=False
        )
        rrr.print_tables(print_intensities=False)

        # print(
        #     # H2       H1
        #     # 113.097, -305.887
        #     219475 * runner1.hamiltonian.get_2nd_order_freqs(
        #         runner1.states.state_list,
        #         V_terms=[G0, 0, 0],
        #         G_terms=[G0, G1, G2]
        #     )
        # )

        """
          0 0 0    176.43652   -163.08293
          0 0 1    289.53377   -468.96990
          0 1 0    298.51621   -466.47580
          1 0 0    170.71509   -206.94025
          """
        raise Exception(...)


        """
              0 0 0    176.43652   -163.08293
              0 0 1    289.53377   -468.96990
              0 1 0    298.51621   -466.47580
              1 0 0    170.71509   -206.94025
              0 0 2    455.91832   -926.22498
              0 2 0    479.53208   -916.91148
              2 0 0    153.91367   -267.74859
              0 1 1    529.20837  -1056.32539
              1 0 1    285.83732   -537.90233
              1 1 0    301.48690   -539.07068
        """
        """
          0 0 0     88.73554   -172.09145
          0 0 1    239.96002   -516.11032
          0 1 0    191.24982   -455.91603
          1 0 0    108.69731   -241.63993
          0 0 2    461.67666  -1028.70222
          0 2 0    353.91552   -887.99927
          2 0 0    138.27476   -348.84148
          0 1 1    455.22336  -1079.05125
          1 0 1    308.63734   -657.43414
          1 1 0    183.18822   -517.48003
          """

        raise Exception(...)

        # print(runner1.hamiltonian.V_terms[0])
        # print(runner2.hamiltonian.V_terms[0])
        # print(runner1.hamiltonian.G_terms[1])
        # print(runner1.hamiltonian.G_terms[1])
        # print(runner1.hamiltonian.V_terms[1])
        # print(runner2.hamiltonian.G_terms[1])
        # print(runner2.hamiltonian.V_terms[1])
        raise Exception(...)

        VPTRunner.run_simple(
            TestManager.test_data(file_name),
            2,
            # include_pseudopotential=False,
            # kinetic_terms=[np.zeros((3,)*(i+2)) for i in range(3)],
            # potential_terms=[np.zeros((3,) * (i + 2)) for i in range(3)],
            internals=[[1, -1, -1, -1], [2, 1, -1, -1], [0, 1, 2, -1]],
            # potential_derivatives=[mol.potential_derivatives[0], mol.potential_derivatives[1]] + [
            #     np.zeros_like(mol.potential_derivatives[2]),
            #     np.zeros_like(mol.potential_derivatives[3])
            # ],
            # expansion_order={'kinetic':0},
            include_pseudopotential=False,
            # potential_terms=[0, 0, 0],
            potential_derivatives=[np.zeros_like(m) if i != 1 else m for i,m in enumerate(mol.potential_derivatives)],
            zero_element_warning=False,
            calculate_intensities=False
        )
        VPTRunner.run_simple(
            TestManager.test_data(file_name),
            2,
            # include_pseudopotential=False,
            # kinetic_terms=[np.zeros((3,)*(i+2)) for i in range(3)],
            # potential_terms=[np.zeros((3,)*(i+2)) for i in range(3)],
            # potential_derivatives=[mol.potential_derivatives[0], mol.potential_derivatives[1]] + [
            #     np.zeros_like(mol.potential_derivatives[2]),
            #     np.zeros_like(mol.potential_derivatives[3])
            # ],
            # expansion_order={'kinetic':0},
            include_pseudopotential=False,
            # potential_terms=[0, 0, 0],
            potential_derivatives=[np.zeros_like(m) if i != 1 else m for i,m in enumerate(mol.potential_derivatives)],
            internals=[[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]],
            zero_element_warning=False,
            calculate_intensities=False
        )
        # VPTRunner.run_simple(
        #     TestManager.test_data(file_name),
        #     3,
        #     memory_constrained=True,
        #     include_pseudopotential=False,
        #     include_coriolis_coupling=False,
        #     logger=True
        # )

    @validationTest
    def test_HOHVPTRunner(self):

        file_name = "HOH_freq.fchk"

        VPTRunner.run_simple(
            TestManager.test_data(file_name),
            3,
            degeneracy_specs='auto'
        )

    @validationTest
    def test_HODVPTRunner(self):

        file_name = "HOD_freq.fchk"

        VPTRunner.run_simple(
            TestManager.test_data(file_name),
            3
        )

    @validationTest
    def test_HOHVPTSubstates(self):

        file_name = "HOH_freq.fchk"
        mol = Molecule.from_file(TestManager.test_data(file_name),
                                 internals=[[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]])
        # raise Exception(mol.internal_coordinates)
        ics = mol.internal_coordinates
        derivs = mol.coords.jacobian(
            ics.system,
            [1, 2],
            all_numerical=True,
            converter_options={'reembed': True}
        )
        derivs = [x.reshape((9,) * (i + 1) + (3, 3)) for i, x in enumerate(derivs)]

        VPTRunner.run_simple(
            # TestManager.test_data(file_name),
            mol,
            # 1,
            [[0, 0, 0], [1, 0, 0]],
            # [[0, 0, 0], [0, 2, 1]],
            # 2,
            order=2,
            expansion_order=2,
            # 3,
            # expansion_order={'default':1, 'dipole':2},
            # target_property='wavefunctions',
            # internals=mol.zmatrix,
            # initial_states=1,
            # operators={
            #     'OH1': [ics[1, 0], derivs[0][:, 1, 0], derivs[1][:, :, 1, 0]],
            #     'OH2': [ics[2, 0], derivs[0][:, 2, 0], derivs[1][:, :, 2, 0]],
            #     'HOH': [ics[2, 1], derivs[0][:, 2, 1], derivs[1][:, :, 2, 1]]
            # },
            logger=True
        )

    @validationTest
    def test_BlockLabels(self):
        VPTRunner.run_simple(
            TestManager.test_data("i_hoh_opt.fchk"),
            VPTStateSpace.get_state_list_from_quanta(4, 6) + [
                [0, 1, 2, 2, 0, 0]
            ],
            initial_states=[
                [0, 0, 0, 0, 0, 0],
                [0, 0, 0, 2, 0, 0],
                [0, 1, 0, 2, 0, 0],
                [0, 0, 0, 0, 1, 0]
            ],
            # degeneracy_specs='auto',
            degeneracy_specs={
                'wfc_threshold':.3
                # "polyads": [
                #     [
                #         [0, 0, 0, 0, 1, 0],
                #         [0, 0, 0, 2, 0, 0]
                #     ],
                #     [
                #         [0, 0, 0, 1, 0, 0],
                #         [0, 0, 2, 0, 0, 0]
                #     ]
                # ]
            },
            # target_property='wavefunctions',
            logger=True,
            # logger=os.path.expanduser("~/Desktop/specks/run.txt"),
            plot_spectrum=False
        )

    @validationTest
    def test_ResultsFileAnalysis(self):

        temp_file = os.path.expanduser('~/Desktop/test_results.hdf5')
        log_file = os.path.expanduser('~/Desktop/test_results.txt')
        # os.remove(temp_file)

        # wfns = VPTRunner.run_simple(
        #     TestManager.test_data("i_hoh_opt.fchk"),
        #     2,
        #     plot_spectrum=False
        #     # initial_states=[
        #     #     [0, 0, 0, 0, 0, 0],
        #     #     [0, 0, 0, 2, 0, 0],
        #     #     [0, 1, 0, 2, 0, 0],
        #     #     [0, 0, 0, 0, 1, 0]
        #     # ]
        # )

        if not os.path.exists(temp_file):
            VPTRunner.run_simple(
                TestManager.test_data("i_hoh_opt.fchk"),
                2,
                # initial_states=[
                #     [0, 0, 0, 0, 0, 0],
                #     [0, 0, 0, 2, 0, 0],
                #     [0, 1, 0, 2, 0, 0],
                #     [0, 0, 0, 0, 1, 0]
                # ],
                # degeneracy_specs='auto',
                degeneracy_specs={
                    "polyads": [
                        [
                            [0, 0, 0, 0, 1, 0],
                            [0, 0, 0, 2, 0, 0]
                        ],
                        [
                            [0, 0, 0, 1, 0, 0],
                            [0, 0, 2, 0, 0, 0]
                        ]
                    ],
                    "extra_groups": [
                        [
                            [0, 0, 0, 0, 1, 0],
                            [0, 1, 0, 0, 1, 0],
                            [1, 0, 0, 0, 1, 0],
                            [0, 0, 0, 2, 0, 0],
                            [0, 1, 0, 2, 0, 0],
                            [1, 0, 0, 2, 0, 0],
                            [0, 0, 2, 1, 0, 0],
                            [0, 1, 2, 1, 0, 0],
                            [1, 0, 2, 1, 0, 0],
                            [0, 0, 4, 0, 0, 0],
                            [0, 1, 4, 0, 0, 0],
                            [1, 0, 4, 0, 0, 0]
                        ]
                    ]
                },
                # target_property='wavefunctions',
                # logger=os.path.expanduser("~/Desktop/specks/run_wfns.txt"),
                results=temp_file,
                logger=log_file,
                plot_spectrum=False
            )

        analyzer = VPTAnalyzer(temp_file)
        # target_state = [0, 0, 2, 0, 0, 0]
        # for i,block in analyzer.degenerate_states:
        #     # intersect([target_state], block) -> check non-empty
        #     ...
        # McUtils.Numputils.intersection()

        shifted_spec = analyzer.shifted_transformed_spectrum(
            analyzer.degenerate_states[4], # 2 states, bend and OOP overtone
            analyzer.deperturbed_hamiltonians[4], # 2x2 matrix
            [0, -50 / UnitsData.hartrees_to_wavenumbers] # the shifts I want to add onto the diagonal
        )
        shifted_spec.plot()#.show()
        print(shifted_spec.frequencies, shifted_spec.intensities)

        with analyzer.log_parser as parser:
            for i, block in enumerate(parser.get_blocks()):
                for subblock in block.lines:
                    print(subblock.tag)

        from McUtils.Scaffolding import LogParser
        with LogParser(log_file) as parser:
            for i, block in enumerate(parser.get_blocks()):
                for subblock in block.lines:
                    print(subblock.tag)

    @validationTest
    def test_IHOHExcited(self):
        wfns = VPTRunner.run_simple(
            TestManager.test_data("i_hoh_opt.fchk"),
            VPTStateSpace.get_state_list_from_quanta(4, 6) + [
                [0, 1, 2, 2, 0, 0]
            ],
            initial_states=[
                [0, 0, 0, 0, 0, 0],
                [0, 0, 0, 2, 0, 0],
                [0, 1, 0, 2, 0, 0],
                [0, 0, 0, 0, 1, 0]
            ],
            # degeneracy_specs='auto',
            degeneracy_specs = {
                "polyads":[
                    [
                        [0, 0, 0, 0, 1, 0],
                        [0, 0, 0, 2, 0, 0]
                    ],
                    [
                        [0, 0, 0, 1, 0, 0],
                        [0, 0, 2, 0, 0, 0]
                    ]
                ],
                "extra_groups": [
                    [
                        [0, 0, 0, 0, 1, 0],
                        [0, 1, 0, 0, 1, 0],
                        [1, 0, 0, 0, 1, 0],
                        [0, 0, 0, 2, 0, 0],
                        [0, 1, 0, 2, 0, 0],
                        [1, 0, 0, 2, 0, 0],
                        [0, 0, 2, 1, 0, 0],
                        [0, 1, 2, 1, 0, 0],
                        [1, 0, 2, 1, 0, 0],
                        [0, 0, 4, 0, 0, 0],
                        [0, 1, 4, 0, 0, 0],
                        [1, 0, 4, 0, 0, 0]
                    ]
                    # [
                    #     [0, 0, 0, 1, 1, 0],
                    #     [0, 1, 0, 1, 1, 0],
                    #     [1, 0, 0, 1, 1, 0],
                    #     [0, 0, 0, 3, 0, 0],
                    #     [0, 1, 0, 3, 0, 0],
                    #     [1, 0, 0, 3, 0, 0],
                    #     [0, 0, 2, 2, 0, 0],
                    #     [0, 1, 2, 2, 0, 0],
                    #     [1, 0, 2, 2, 0, 0],
                    #     [0, 0, 4, 1, 0, 0],
                    #     [0, 1, 4, 1, 0, 0],
                    #     [1, 0, 4, 1, 0, 0]
                    # ]
                ]
            },
            target_property='wavefunctions',
            logger=os.path.expanduser("~/Desktop/specks/run_wfns.txt"),
            # logger=os.path.expanduser("~/Desktop/specks/run.txt"),
            plot_spectrum=False
        )
        # raise Exception(wfns.initial_states, wfns.initial_state_indices)

        multispec = wfns.get_spectrum().frequency_filter(600, 4400)
        multispec.plot().savefig(
                os.path.expanduser(f"~/Desktop/specks/full.pdf"),
                transparent=True
            )
        for state,spec in zip(wfns.initial_states, multispec):
            s = "".join(str(s) for s in state)
            spec.plot(plot_range=[[600, 4400], [0, 500]], padding=[[0, 0], [0, 0]],
                      image_size=[.75 * 20, 5.48 * 20]
                      ).savefig(
                os.path.expanduser(f"~/Desktop/specks/state_{s}.pdf"),
                transparent=True
            )
        multispec = wfns.get_deperturbed_spectrum().frequency_filter(600, 4400)
        for state,spec in zip(wfns.initial_states, multispec):
            s = "".join(str(s) for s in state)
            spec.plot(plot_range=[[600, 4400], [0, 500]], padding=[[0, 0], [0, 0]],
                      image_size=[.75 * 20, 5.48 * 20]
                      ).savefig(
                os.path.expanduser(f"~/Desktop/specks/state_depert_{s}.pdf"),
                transparent=True
            )

    @validationTest
    def test_HOHVPTAnneManip(self):

        runner, _ = VPTRunner.helpers.run_anne_job(
            TestManager.test_data("vpt2_helpers_api/hod/r"),
            return_runner=True,
            order=2,
            expansion_order=2
        )
        # runner, opts = VPTRunner.construct('HOH', 3)

        # Collect expansion data from runner
        H = runner.hamiltonian
        V = H.V_terms
        freqs = np.diag(V[0]) # how they actually get fed into the code...
        G = H.G_terms
        U = H.pseudopotential_term
        D = H.expansion_options['dipole_terms'] # these usually get fed forward to the wave functions

        # raise Exception([[x.shape if isinstance(x, np.ndarray) else x for x in D[a]] for a in range(3)])

        # Define new shifted frequencies
        frequency_shift = np.array([-1, 0, 0]) * UnitsData.convert("Wavenumbers", "Hartrees")
        new_freqs = freqs + frequency_shift

        # Rescale parameters
        scaling_factor = np.sqrt(new_freqs) / np.sqrt(freqs)
        # we use NumPy broadcasting tricks to rescale everything

        G[2] # this requires the highest-order derivatives, so by doing it first and letting everything
             # cache less junk gets printed to screen

        v_expansion = [
            # d V /dq_i dq_j
            V[0] * (scaling_factor[:, np.newaxis] * scaling_factor[np.newaxis, :]),
            # d V /dq_i dq_j dq_k
            V[1] * (
                    scaling_factor[:, np.newaxis, np.newaxis] *
                    scaling_factor[np.newaxis, :, np.newaxis] *
                    scaling_factor[np.newaxis, np.newaxis, :]
            ),
            # d V /dq_i dq_j dq_k dq_l
            V[2] * (
                    scaling_factor[:, np.newaxis, np.newaxis, np.newaxis] *
                    scaling_factor[np.newaxis, :, np.newaxis, np.newaxis] *
                    scaling_factor[np.newaxis, np.newaxis, :, np.newaxis] *
                    scaling_factor[np.newaxis, np.newaxis, np.newaxis, :]
            ),
        ]

        # For the momentum axes we _divide_ by the scaling factor
        g_expansion = [
            # Formally we should be dividing by this scaling factor, but to make sure
            # V[0] == G[0] we multiply
            # G_i,j
            G[0] * (scaling_factor[:, np.newaxis] * scaling_factor[np.newaxis, :]),
            # For the derivatives, we do the scaling correct
            # d G_j,k / dq_i (i.e. q-index corresponds to axis 0)
            G[1] * (
                    scaling_factor[:, np.newaxis, np.newaxis] /
                    scaling_factor[np.newaxis, :, np.newaxis] /
                    scaling_factor[np.newaxis, np.newaxis, :]
            ),
            # d G_k,l / dq_i dq_j (i.e. q-indices correspond to axes 0,1)
            G[2] * (
                    scaling_factor[:, np.newaxis, np.newaxis, np.newaxis] * # Note that we multiply for the first two axes
                    scaling_factor[np.newaxis, :, np.newaxis, np.newaxis] /
                    scaling_factor[np.newaxis, np.newaxis, :, np.newaxis] /
                    scaling_factor[np.newaxis, np.newaxis, np.newaxis, :]
            )
        ]

        # I haven't done the math to figure out how exactly u should transform
        u_expansion = [U[0]]

        d_expansion = [
            [
                D[a][0],  # mu_a
                # d mu_a / dq_i
                D[a][1] * scaling_factor[:],
                # d mu_a / dq_i dq_j
                D[a][2] * (
                        scaling_factor[:, np.newaxis] *
                        scaling_factor[np.newaxis, :]
                ),
                # d mu_a / dq_i dq_j dq_k
                D[a][3] * (
                        scaling_factor[:, np.newaxis, np.newaxis] *
                        scaling_factor[np.newaxis, :, np.newaxis] *
                        scaling_factor[np.newaxis, np.newaxis, :]
                )
            ]
            for a in range(3)  # loop over the x, y, and z axes
        ]

        new_runner, _ = VPTRunner.construct(
            runner.system.mol, # not actually used
            runner.states.state_list,
            potential_terms=v_expansion,
            kinetic_terms=g_expansion,
            pseudopotential_terms=u_expansion,
            dipole_terms=d_expansion
        )

        runner.print_tables() # runs the code and prints the IR tables
        new_runner.print_tables()

    @validationTest
    def test_HOHPartialQuartic(self):

        VPTRunner.helpers.run_anne_job(
            # os.path.expanduser("~/Desktop/r_as"),
            TestManager.test_data("vpt2_helpers_api/hod/x"),
            # states=2, # max quanta to be focusing on
            states=[[0, 0, 0], [0, 0, 1], [0, 1, 0]],
            # max quanta to be focusing on
            order=2,  # None, # orderr of VPT
            expansion_order=2,  # None, # order of expansion of H can use {
            # 'potential':int,
            # 'kinetic':int,
            # 'dipole':int}
            # logger=filename
            # mode_selection=[1, 2],
            calculate_intensities=False,
            zero_element_warning=False,
            include_only_mode_couplings=[1, 2],
            include_coriolis_coupling=False
            # target_property='wavefunctions'
            # return_runner=True
        )  # output file name

        VPTRunner.helpers.run_anne_job(
            # os.path.expanduser("~/Desktop/r_as"),
            TestManager.test_data("vpt2_helpers_api/hod/x_no_bend"),
            # states=2, # max quanta to be focusing on
            states=[[0, 0, 0], [0, 0, 1], [0, 1, 0]],
            # max quanta to be focusing on
            order=2,  # None, # orderr of VPT
            expansion_order=2,  # None, # order of expansion of H can use {
            # 'potential':int,
            # 'kinetic':int,
            # 'dipole':int}
            # logger=filename
            # mode_selection=[1, 2],
            calculate_intensities=False,
            zero_element_warning=False,
            # include_only_mode_couplings=[1, 2],
            include_coriolis_coupling=True
            # target_property='wavefunctions'
            # return_runner=True
        )  # output file name

        raise Exception(...)

        VPTRunner.helpers.run_anne_job(
            # os.path.expanduser("~/Desktop/r_as"),
            TestManager.test_data("vpt2_helpers_api/hod/x_decoupled"),
            # states=2, # max quanta to be focusing on
            states=1,
            # max quanta to be focusing on
            order=2,  # None, # orderr of VPT
            expansion_order=2,  # None, # order of expansion of H can use {
            # 'potential':int,
            # 'kinetic':int,
            # 'dipole':int}
            # logger=filename
            # mode_selection=[1, 2],
            calculate_intensities=False,
            zero_element_warning=False,
            include_coriolis_coupling=False
            # target_property='wavefunctions'
            # return_runner=True
        )  # output file name

        """
::> States Energies
  > State     Harmonic   Anharmonic     Harmonic   Anharmonic
               ZPE          ZPE    Frequency    Frequency
0 0 0   4052.91097   4001.04707            -            - 
0 0 1            -            -   3873.84521   3688.94993 
0 1 0            -            -   2810.03028   2723.42011 
1 0 0            -            -   1421.94645   1383.13454 
"""

        VPTRunner.helpers.run_anne_job(
            # os.path.expanduser("~/Desktop/r_as"),
            TestManager.test_data("vpt2_helpers_api/hod/x"),
            # states=2, # max quanta to be focusing on
            states=1,
            # max quanta to be focusing on
            order=2,  # None, # orderr of VPT
            expansion_order=2,  # None, # order of expansion of H can use {
            # 'potential':int,
            # 'kinetic':int,
            # 'dipole':int}
            # logger=filename
            mode_selection=[0, 2],
            calculate_intensities=False,
            zero_element_warning=False,
            include_coriolis_coupling=False
            # target_property='wavefunctions'
            # return_runner=True
        )  # output file name

        # runner3 = VPTRunner.helpers.run_anne_job(
        #     # os.path.expanduser("~/Desktop/r_as"),
        #     TestManager.test_data("vpt2_helpers_api/hod/x"),
        #     # states=2, # max quanta to be focusing on
        #     states=1,
        #     # max quanta to be focusing on
        #     order=2,  # None, # orderr of VPT
        #     expansion_order=2,  # None, # order of expansion of H can use {
        #     # 'potential':int,
        #     # 'kinetic':int,
        #     # 'dipole':int}
        #     # logger=filename
        #     # mode_selection=[1, 2],
        #     calculate_intensities=False,
        #     operator_coefficient_threshold=1e-12,
        #     zero_element_warning=False,
        #     # target_property='wavefunctions'
        #     return_runner=True
        # )  # output file name
        #
        # runner4 = VPTRunner.helpers.run_anne_job(
        #     # os.path.expanduser("~/Desktop/r_as"),
        #     TestManager.test_data("vpt2_helpers_api/hod/x_sub"),
        #     # states=2, # max quanta to be focusing on
        #     states=1,
        #     # max quanta to be focusing on
        #     order=2,  # None, # orderr of VPT
        #     expansion_order=2,  # None, # order of expansion of H can use {
        #     # 'potential':int,
        #     # 'kinetic':int,
        #     # 'dipole':int}
        #     # logger=filename
        #     # mode_selection=[1, 2],
        #     calculate_intensities=False,
        #     operator_coefficient_threshold=-1,
        #     zero_element_warning=False
        #     # target_property='wavefunctions'
        #     , return_runner=True
        # )  # output file name

        """
0 0 0   4052.91097   3994.84632            -            - 
0 0 1            -            -   3873.84521   3685.79215 
0 1 0            -            -   2810.03028   2706.14630 
1 0 0            -            -   1421.94645   1383.40761 

0 0 0   4052.91097   4033.89584            -            - 
0 0 1            -            -   3873.84521   3697.97351 
0 1 0            -            -   2810.03028   2902.77195 
1 0 0            -            -   1421.94645   1380.70462
"""
        # split1 = runner3[0].hamiltonian.get_Nielsen_energies([[0, 0, 0], [0, 0, 1]], return_split=True)
        # split2 = runner4[0].hamiltonian.get_Nielsen_energies([[0, 0, 0], [0, 0, 1]], return_split=True)
        # raise Exception(
        #     split1[2][0] - split2[2][0] # cubic contributions to X matrix
        # )
        # nie_full = np.sum(runner3[0].hamiltonian.get_Nielsen_energies([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]), axis=0)
        # nie_sub = np.sum(runner4[0].hamiltonian.get_Nielsen_energies([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]), axis=0)
        """
        array([1383.40761359, 2706.14629421, 3685.79215281])
        array([1380.70462067, 2902.77194309, 3697.97351711]))
        """
        raise Exception(...)
        #     (nie_full[1:] - nie_full[0]) * UnitsData.hartrees_to_wavenumbers,
        #     (nie_sub[1:] - nie_sub[0]) * UnitsData.hartrees_to_wavenumbers
        # )


        raise Exception(...)
        """
        0   1973.85245   1992.46835            -            - 
        1            -            -   3947.70491   4042.24072
        """

        runner1 = VPTRunner.helpers.run_anne_job(
            os.path.expanduser("~/Desktop/r_as"),
            # states=2, # max quanta to be focusing on
            states=1,
            # max quanta to be focusing on
            order=2,  # None, # orderr of VPT
            expansion_order=2,  # None, # order of expansion of H can use {
            # 'potential':int,
            # 'kinetic':int,
            # 'dipole':int}
            # logger=filename
            calculate_intensities=False,
            operator_coefficient_threshold=-1
            # target_property='wavefunctions'
            # return_runner=True
        )  # output file name

        """
  State       Frequency    Intensity       Frequency    Intensity
  0 0 1    3947.69802     75.46200      3761.26977     70.72838
  0 1 0    3821.87392      5.56098      3652.31667      4.82837
  1 0 0    1628.37574      0.00000      1590.97610      0.00483
  State       Frequency    Intensity       Frequency    Intensity
  0 0 1    3947.69802     75.46200      3874.95435     71.78714
  0 1 0    3821.87392      0.00000      3864.33291      0.00111
  1 0 0    1628.37574      0.00000      1620.53058      0.00003
"""

        """
        ::> States Energies
  > State     Harmonic   Anharmonic     Harmonic   Anharmonic
               ZPE          ZPE    Frequency    Frequency
0 0 0   4698.97384   4628.17891            -            - 
0 0 1            -            -   3947.69802   3761.26977 
0 1 0            -            -   3821.87392   3652.31667 
1 0 0            -            -   1628.37574   1590.97610 
"""

        runner2 = VPTRunner.helpers.run_anne_job(
            os.path.expanduser("~/Desktop/r_a"),
            # states=2, # max quanta to be focusing on
            states=1,
            # max quanta to be focusing on
            order=2,  # None, # orderr of VPT
            expansion_order=2,  # None, # order of expansion of H can use {
            # 'potential':int,
            # 'kinetic':int,
            # 'dipole':int}
            # logger=filename
            calculate_intensities=False,
            operator_coefficient_threshold=-1
            # target_property='wavefunctions'
            # return_runner=True
        )  # output file name

        raise Exception(...)

        raise Exception(
            runner1[0].ham_opts.opts['potential_terms'][1],
            runner2[0].ham_opts.opts['potential_terms'][1]
            # runner2.ham_opts['potential_derivatives'],
        )

    @validationTest
    def test_HOHVPTNonGSRunner(self):

        file_name = "HOH_freq.fchk"
        mol = Molecule.from_file(TestManager.test_data(file_name),
                                 internals=[[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]])
        # raise Exception(mol.internal_coordinates)
        ics = mol.internal_coordinates
        derivs = mol.coords.jacobian(
            ics.system,
            [1, 2],
            all_numerical=True,
            converter_options={'reembed':True}
        )
        derivs = [x.reshape((9,)*(i+1) + (3, 3)) for i,x in enumerate(derivs)]

        VPTRunner.run_simple(
            TestManager.test_data(file_name),
            2,
            # target_property='wavefunctions',
            # internals=mol.zmatrix,
            initial_states=1,
            operators={
                'OH1':[ics[1, 0], derivs[0][:, 1, 0], derivs[1][:, :, 1, 0]],
                'OH2':[ics[2, 0], derivs[0][:, 2, 0], derivs[1][:, :, 2, 0]],
                'HOH':[ics[2, 1], derivs[0][:, 2, 1], derivs[1][:, :, 2, 1]]
            },
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
            degeneracy_specs='auto',
            corrected_fundamental_frequencies=np.array([1600, 3775, 3880])/UnitsData.convert("Hartrees", "Wavenumbers")
        )

    @validationTest
    def test_OCHHSubspaceTargetProps(self):

        file_name = "OCHH_freq.fchk"
        VPTRunner.run_simple(
            TestManager.test_data(file_name),
            [
                [0, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 1],
                [0, 1, 2, 0, 0, 0]
            ],
            # expansion_order=1,
            logger=True,
            target_property='wavefunctions'
            # , calculate_intensities=False
            , degeneracy_specs='auto',
            # extended_space_target_property='frequencies'
        )
        """
State           Harmonic   Anharmonic     Harmonic   Anharmonic
                     ZPE          ZPE    Frequency    Frequency
0 0 0 0 0 0   5866.87156   5785.95861            -            - 
0 0 1 0 0 1            -            -   4588.74227   4415.63305 
0 1 2 0 0 0            -            -   4306.24556   4191.96329 
"""
        """
::> IR Data
  > Initial State: 0 0 0 0 0 0 
                         Harmonic                  Anharmonic
State             Frequency    Intensity       Frequency    Intensity
  0 0 1 0 0 1    4588.74227      0.00000      4446.79734      0.00164
  0 1 2 0 0 0    4306.24556      0.00000      4167.25991      0.30489
  0 1 1 1 0 0    4506.28742      0.00000      4345.62126      0.48928
  """

        VPTRunner.run_simple(
            TestManager.test_data(file_name),
            [
                [0, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 1],
                [0, 1, 2, 0, 0, 0]
            ],
            # expansion_order=1,
            logger=True,
            target_property='intensities'
            # , calculate_intensities=False
            , degeneracy_specs='auto'
            # extended_space_target_property='frequencies'
        )
        """
::> IR Data
  > Initial State: 0 0 0 0 0 0 
                         Harmonic                  Anharmonic
State             Frequency    Intensity       Frequency    Intensity
  0 0 1 0 0 1    4588.74227      0.00000      4446.79734      0.00164
  0 1 2 0 0 0    4306.24556      0.00000      4167.25991      0.30489
  0 1 1 1 0 0    4506.28742      0.00000      4345.62126      0.48928
  """

        VPTRunner.run_simple(
            TestManager.test_data(file_name),
            [
                [0, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 1],
                [0, 1, 2, 0, 0, 0]
            ],
            # expansion_order=1,
            logger=True,
            target_property='frequencies'
            # , calculate_intensities=False
            , degeneracy_specs='auto'
        )
        """
::> IR Data
  > Initial State: 0 0 0 0 0 0 
                         Harmonic                  Anharmonic
State             Frequency    Intensity       Frequency    Intensity
  0 0 1 0 0 1    4588.74227      0.00000      4458.35353      0.00402
  0 1 2 0 0 0    4306.24556      0.00000      4170.69772      0.23905
  0 1 1 1 0 0    4506.28742      0.00000      4330.62726      0.55570
  """

    @validationTest
    def test_OCHHFasterDegen(self):

        file_name = "OCHH_freq.fchk"
        # VPTRunner.run_simple(
        #     TestManager.test_data(file_name),
        #     2,
        #     logger=True,
        #     degeneracy_specs='auto',
        #     extended_space_target_property='frequencies'
        # )

        VPTRunner.run_simple(
            TestManager.test_data(file_name),
            2,
            logger=True,
            degeneracy_specs='auto'
        )

    @validationTest
    def test_OCHHFasterDegenSubspace(self):

        file_name = "OCHH_freq.fchk"
        VPTRunner.run_simple(
            TestManager.test_data(file_name),
            [
                [0, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 1]
            ],
            logger=True,
            degeneracy_specs='auto',
            extended_space_target_property='frequencies'
        )
        """
::> IR Data
  > Initial State: 0 0 0 0 0 0 
                         Harmonic                  Anharmonic
State             Frequency    Intensity       Frequency    Intensity
  0 0 1 0 0 1    4588.74227      0.00000      4446.79734      0.00134
  0 1 2 0 0 0    4306.24556      0.00000      4167.25991      0.29964
  0 1 1 1 0 0    4506.28742      0.00000      4345.62126      0.48923
  """

        VPTRunner.run_simple(
            TestManager.test_data(file_name),
            [
                [0, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 1]
            ],
            logger=True,
            degeneracy_specs='auto'
        )
        """
::> IR Data
  > Initial State: 0 0 0 0 0 0 
                         Harmonic                  Anharmonic
State             Frequency    Intensity       Frequency    Intensity
  0 0 1 0 0 1    4588.74227      0.00000      4446.79734      0.00000
  0 1 2 0 0 0    4306.24556      0.00000      4167.25991      0.26680
  0 1 1 1 0 0    4506.28742      0.00000      4345.62126      0.47084
  """

    @validationTest
    def test_HOONOFasterDegen(self):

        file_name = "HOONO_freq.fchk"
        VPTRunner.run_simple(
            TestManager.test_data(file_name),
            2,
            logger=True,
            degeneracy_specs='auto',
            extended_space_target_property='frequencies'
        )

        VPTRunner.run_simple(
            TestManager.test_data(file_name),
            2,
            logger=True,
            degeneracy_specs='auto'
        )

    @validationTest
    def test_HOONOFasterDegenSubspace(self):

        file_name = "HOONO_freq.fchk"
        VPTRunner.run_simple(
            TestManager.test_data(file_name),
            [[0, 0, 0, 0, 0, 0, 0, 0, 0],
             [1, 0, 0, 0, 0, 0, 1, 0, 0],
             [1, 0, 0, 0, 0, 0, 0, 1, 0],
             [0, 0, 1, 0, 0, 0, 1, 0, 0],
             [0, 0, 1, 0, 0, 0, 0, 1, 0]],
            logger=True,
            degeneracy_specs='auto',
            extended_space_target_property='frequencies'
        )
        """
State                   Frequency    Intensity       Frequency    Intensity
  1 0 0 0 0 0 1 0 0    1789.34288      0.00000      1669.11650      0.00850
  1 0 0 0 0 0 0 1 0    1924.23563      0.00000      1825.06905      0.30108
  0 0 1 0 0 0 1 0 0    1957.70875      0.00000      1886.51516      0.05466
  0 0 1 0 0 0 0 1 0    2092.60150      0.00000      2045.13170      0.02941
  1 0 0 2 0 0 0 0 0    1787.50010      0.00000      2413.94904     80.96353
  0 0 1 2 0 0 0 0 0    1955.86597      0.00000      1971.94764      2.37340
  1 0 0 1 1 0 0 0 0    1908.01157      0.00000      2131.62528      0.72881
  0 0 1 1 1 0 0 0 0    2076.37744      0.00000      2080.76165      2.72838
  """

        VPTRunner.run_simple(
            TestManager.test_data(file_name),
            [[0, 0, 0, 0, 0, 0, 0, 0, 0],
             [1, 0, 0, 0, 0, 0, 1, 0, 0],
             [1, 0, 0, 0, 0, 0, 0, 1, 0],
             [0, 0, 1, 0, 0, 0, 1, 0, 0],
             [0, 0, 1, 0, 0, 0, 0, 1, 0]],
            logger=True,
            degeneracy_specs='auto'
        )
        """
State                   Frequency    Intensity       Frequency    Intensity
  1 0 0 0 0 0 1 0 0    1789.34288      0.00000      1669.11650      0.04924
  1 0 0 0 0 0 0 1 0    1924.23563      0.00000      1825.06905      0.07528
  0 0 1 0 0 0 1 0 0    1957.70875      0.00000      1886.51516      0.00766
  0 0 1 0 0 0 0 1 0    2092.60150      0.00000      2045.13170      0.01768
  1 0 0 2 0 0 0 0 0    1787.50010      0.00000      2413.94904      0.08504
  0 0 1 2 0 0 0 0 0    1955.86597      0.00000      1971.94764      0.00159
  1 0 0 1 1 0 0 0 0    1908.01157      0.00000      2131.62528      0.00035
  0 0 1 1 1 0 0 0 0    2076.37744      0.00000      2080.76165      0.00332
        """

    @validationTest
    def test_OCHHVPTRunnerShifted(self):
        file_name = "OCHH_freq.fchk"
        VPTRunner.run_simple(
            TestManager.test_data(file_name),
            2,
            logger=True,
            degeneracy_specs='auto',
            corrected_fundamental_frequencies=np.array([1188, 1252, 1527, 1727, 2977, 3070]) / UnitsData.convert("Hartrees", "Wavenumbers")
        )

    @validationTest
    def test_HOONOVPTRunnerShifted(self):
        file_name = "HOONO_freq.fchk"
        VPTRunner.run_simple(
            TestManager.test_data(file_name),
            2,
            logger=True,
            degeneracy_specs='auto',
            corrected_fundamental_frequencies=np.array([
                355.73348, 397.16760, 524.09935,
                715.88331, 836.39478, 970.87676,
                1433.60940, 1568.50215, 3486.85528
            ]) / UnitsData.convert("Hartrees", "Wavenumbers")
        )

    @inactiveTest
    def test_CrieegeeVPTRunnerShifted(self):
        # with BlockProfiler('Crieegee', print_res=True):
        freqs = VPTSystem('criegee_eq_anh.fchk').mol.normal_modes.modes.freqs
        freqs = freqs.copy()
        freqs[1] += 10/UnitsData.convert("Hartrees", "Wavenumbers")
        VPTRunner.run_simple(
            # 'criegee_eq_anh.fchk',
            2,
            logger=True,
            degeneracy_specs='auto',
            corrected_fundamental_frequencies=freqs
            # corrected_fundamental_frequencies=np.array([
            #     200.246, 301.985 + 10, 462.536, 684.792, 736.234, 961.474, 984.773, 1038.825, 1120.260, 1327.450, 1402.397,
            #     1449.820, 1472.576, 1519.875, 3037.286, 3078.370, 3174.043, 3222.828
            # ])/UnitsData.convert("Hartrees", "Wavenumbers")
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

        # raise Exception(model.g())

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

    @validationTest
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

    @validationTest
    def test_HOHCorrectedPostfilters(self):

        # VPTAnalyzer.run_VPT(
        #     TestManager.test_data('HOH_freq.fchk'),
        #     2,
        #     zero_order_energy_corrections=[
        #         [(0, 1, 0), (4681.56364 + 3800) * UnitsData.convert("Wavenumbers", "Hartrees")],
        #         [(0, 2, 0), (4681.56364 + 7800) * UnitsData.convert("Wavenumbers", "Hartrees")],
        #         [(0, 0, 2), (4681.56364 + 7801) * UnitsData.convert("Wavenumbers", "Hartrees")],
        #     ],
        #     degeneracy_specs='auto',
        #     basis_filters={'max_quanta':[0, -1, -1]}
        # ).print_output_tables()

        VPTAnalyzer.run_VPT(
            TestManager.test_data('HOH_freq.fchk'),
            2,
            zero_order_energy_corrections=[
                [(0, 1, 0), (4681.56364 + 3800) * UnitsData.convert("Wavenumbers", "Hartrees")],
                [(0, 2, 0), (4681.56364 + 7800) * UnitsData.convert("Wavenumbers", "Hartrees")],
                [(0, 0, 2), (4681.56364 + 7801) * UnitsData.convert("Wavenumbers", "Hartrees")],
            ],
            degeneracy_specs='auto',
            basis_filters=[
                {'max_quanta': [2, -1, -1]},
                {'excluded_transitions': [[1, 0, 0], [0, 1, 0], [0, 0, 1]]}
            ]
        ).print_output_tables()

        # VPTAnalyzer.run_VPT(
        #     TestManager.test_data('HOH_freq.fchk'),
        #     2,
        #     zero_order_energy_corrections=[
        #         [(0, 1, 0), (4681.56364 + 3800) * UnitsData.convert("Wavenumbers", "Hartrees")],
        #         [(0, 2, 0), (4681.56364 + 7800) * UnitsData.convert("Wavenumbers", "Hartrees")],
        #         [(0, 0, 2), (4681.56364 + 7801) * UnitsData.convert("Wavenumbers", "Hartrees")],
        #     ],
        #     degeneracy_specs='auto',
        #     basis_filters={'max_quanta': [4, -1, -1]}
        # ).print_output_tables()

        # VPTAnalyzer.run_VPT(
        #     TestManager.test_data('HOH_freq.fchk'),
        #     2,
        #     zero_order_energy_corrections=[
        #         [(0, 1, 0), (4681.56364 + 3800) * UnitsData.convert("Wavenumbers", "Hartrees")],
        #         [(0, 2, 0), (4681.56364 + 7800) * UnitsData.convert("Wavenumbers", "Hartrees")],
        #         [(0, 0, 2), (4681.56364 + 7801) * UnitsData.convert("Wavenumbers", "Hartrees")],
        #     ],
        #     degeneracy_specs='auto',
        #     basis_filters={'max_quanta': [6, -1, -1]}
        # ).print_output_tables()
        #
        # VPTAnalyzer.run_VPT(
        #     TestManager.test_data('HOH_freq.fchk'),
        #     2,
        #     zero_order_energy_corrections=[
        #         [(0, 1, 0), (4681.56364 + 3800) * UnitsData.convert("Wavenumbers", "Hartrees")],
        #         [(0, 2, 0), (4681.56364 + 7800) * UnitsData.convert("Wavenumbers", "Hartrees")],
        #         [(0, 0, 2), (4681.56364 + 7801) * UnitsData.convert("Wavenumbers", "Hartrees")],
        #     ],
        #     degeneracy_specs='auto',
        #     basis_filters={'max_quanta': [-1, -1, -1]}
        # ).print_output_tables()

    @validationTest
    def test_WaterSkippedCouplings(self):

        VPTRunner.run_simple(
            TestManager.test_data('water_freq.fchk'),
            1,
            degeneracy_specs='auto',
            operator_coefficient_threshold=(1.0e-8)
        )  # .print_output_tables()

    @validationTest
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

    @validationTest
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

        # VPTRunner.run_simple(
        #     TestManager.test_data('OCHH_freq.fchk'),
        #     1,
        #     degeneracy_specs='auto',
        #     internals=internals,
        #     mode_selection=[1, 2, 3, 4, 5]
        # )

    @validationTest
    def test_HODRephase(self):
        VPTRunner.run_simple(
            TestManager.test_data('HOD_freq_16.fchk'),
            1,
            degeneracy_specs=None,
            # order=4,
            # expansion_order=2
        )

    @validationTest
    def test_HOHRephase(self):
        VPTRunner.run_simple(
            TestManager.test_data('HOH_freq.fchk'),
            1,
            degeneracy_specs=None,
            # order=4,
            # expansion_order=2
        )

    @validationTest
    def test_NH3(self):

        VPTRunner.run_simple(
            TestManager.test_data('nh3.fchk'),
            2,
            # degeneracy_specs=False,
            order=4,
            expansion_order=2,

            # basis_filters={
            #     'max_quanta':[2, -1, -1, -1, -1, -1]
            # }
        )

    @validationTest
    def test_HOONO(self):

        VPTRunner.run_simple(
            TestManager.test_data('HOONO_freq.fchk'),
            1,
            degeneracy_specs=None,
            # order=4,
            # expansion_order=2
        )

        # VPTRunner.run_simple(
        #     TestManager.test_data('OCHH_freq.fchk'),
        #     1,
        #     degeneracy_specs='auto',
        #     internals=internals,
        #     mode_selection=[1, 2, 3, 4, 5]
        # )

    @validationTest
    def test_H2COSkippedCouplings(self):

        VPTRunner.run_simple(
            TestManager.test_data('OCHH_freq.fchk'),
            1,
            degeneracy_specs='auto',
            operator_coefficient_threshold=1.00 / 219475
        )

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

    @validationTest
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
        )#.print_output_tables()

        # VPTAnalyzer.run_VPT(
        #     TestManager.test_data('HOH_freq.fchk'),
        #     2,
        #     zero_order_energy_corrections=[
        #         [(0, 1, 0), (4681.56364 + 3800) * UnitsData.convert("Wavenumbers", "Hartrees")],
        #         [(0, 2, 0), (4681.56364 + 7800) * UnitsData.convert("Wavenumbers", "Hartrees")],
        #         [(0, 0, 2), (4681.56364 + 7801) * UnitsData.convert("Wavenumbers", "Hartrees")],
        #     ],
        #     degeneracy_specs='auto',
        #     basis_filters={'max_quanta': [4, -1, -1]}
        # ).print_output_tables()

        # VPTAnalyzer.run_VPT(
        #     TestManager.test_data('HOH_freq.fchk'),
        #     2,
        #     zero_order_energy_corrections=[
        #         [(0, 1, 0), (4681.56364 + 3800) * UnitsData.convert("Wavenumbers", "Hartrees")],
        #         [(0, 2, 0), (4681.56364 + 7800) * UnitsData.convert("Wavenumbers", "Hartrees")],
        #         [(0, 0, 2), (4681.56364 + 7801) * UnitsData.convert("Wavenumbers", "Hartrees")],
        #     ],
        #     degeneracy_specs='auto',
        #     basis_filters={'max_quanta': [6, -1, -1]}
        # ).print_output_tables()
        #
        # VPTAnalyzer.run_VPT(
        #     TestManager.test_data('HOH_freq.fchk'),
        #     2,
        #     zero_order_energy_corrections=[
        #         [(0, 1, 0), (4681.56364 + 3800) * UnitsData.convert("Wavenumbers", "Hartrees")],
        #         [(0, 2, 0), (4681.56364 + 7800) * UnitsData.convert("Wavenumbers", "Hartrees")],
        #         [(0, 0, 2), (4681.56364 + 7801) * UnitsData.convert("Wavenumbers", "Hartrees")],
        #     ],
        #     degeneracy_specs='auto',
        #     basis_filters={'max_quanta': [-1, -1, -1]}
        # ).print_output_tables()

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

    @validationTest
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

    @validationTest
    def test_AnneAPI(self):
        VPTRunner.run_simple(
            TestManager.test_data('HOD_freq_16.fchk'), 2,
            # calculate_intensities=False
        )
        # raise Exception("...")
        # test_folder = TestManager.test_data('vpt2_helpers_api/hod')
        # runner, dat = VPTRunner.construct(TestManager.test_data('HOD_freq_16.fchk'), 2)

        # fuck = dat[0]
        # raise Exception(fuck.mol.normal_modes.modes.basis.matrix.T)

        # fuck = runner.hamiltonian
        # shit_strings = []
        # cub = fuck.V_terms[1]
        # h2w = UnitsData.convert("Hartrees", "Wavenumbers")
        # for i in range(cub.shape[0]):
        #     for j in range(i, cub.shape[1]):
        #         for k in range(j, cub.shape[2]):
        #             shit_strings.append(f"{i+1} {j+1} {k+1} {cub[i, j ,k]*h2w}")

        # fuck = runner.hamiltonian
        # shit_strings = []
        # quart = fuck.V_terms[2]
        # h2w = UnitsData.convert("Hartrees", "Wavenumbers")
        # for i in range(quart.shape[0]):
        #     for j in range(i, quart.shape[2]):
        #         for k in range(j, quart.shape[2]):
        #             for l in range(k, quart.shape[3]):
        #                 shit_strings.append(f"{i+1} {j+1} {k+1} {l+1} {quart[i, j, k, l]*h2w}")
        # #
        # raise Exception("\n".join(shit_strings))


        # fuck = runner.get_wavefunctions()
        # shit_strings = []
        # dts = fuck.dipole_terms
        #
        # lints = np.array([d[1] for d in dts])
        # for a in range(lints.shape[0]):
        #     for i in range(lints.shape[1]):
        #         shit_strings.append(f"{a+1} {i+1} {lints[a, i]}")
        # shit_strings.append("")
        #
        # dints = np.array([d[2] for d in dts])
        # for a in range(dints.shape[0]):
        #     for i in range(dints.shape[1]):
        #         for j in range(i, dints.shape[2]):
        #             shit_strings.append(f"{a + 1} {i + 1} {j + 1} {dints[a, i, j]}")
        # shit_strings.append("")
        #
        # cints = np.array([d[3] for d in dts])
        # for a in range(cints.shape[0]):
        #     for i in range(cints.shape[1]):
        #         for j in range(i, cints.shape[2]):
        #             for k in range(j, cints.shape[3]):
        #                 shit_strings.append(f"{a + 1} {i + 1} {j + 1} {k + 1} {cints[a, i, j, k]}")
        # #
        # raise Exception("\n".join(shit_strings))

        VPTRunner.helpers.run_anne_job(
            TestManager.test_data('vpt2_helpers_api/hod/x'),
            # calculate_intensities=False,
            # expansion_order=2
        )
