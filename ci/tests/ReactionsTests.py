import os.path

from Peeves.TestUtils import *
from unittest import TestCase
from Peeves import BlockProfiler

import numpy as np
from Psience.Reactions import Reaction

class ReactionsTests(TestCase):

    @debugTest
    def test_SetupReactants(self):

        # import McUtils.Plots as plt
        # from McUtils.Data import UnitsData
        #
        # plt.Plot(
        #     np.linspace(-1, 1, 10), np.linspace(-1, 1, 10) ** 2 * UnitsData.convert("Hartrees", "Kilojoules/Mole"),
        #     marker='.', linestyle='dashed', color="red",
        #     axes_labels=['Reaction Coordinate', "Energy (kJ mol$^{-1}$)"],
        #     padding=[[60, 20], [60, 0]]
        # ).show()
        # return


        # def get_specs(zmat):
        #     specs = []
        #     for row in zmat:
        #         if row[1] < 0: continue
        #         specs.append(tuple(row[:2]))
        #         if row[2] < 0: continue
        #         specs.append(tuple(row[:3]))
        #         if row[3] < 0: continue
        #         specs.append(tuple(row[:4]))
        #     return specs
        #
        # from Psience.Molecools import Molecule
        # hoono3 = Molecule.from_file(
        #     TestManager.test_data('HOONO_freq.fchk'),
        #     internals={'specs': get_specs([
        #         [0, -1, -1, -1],  # H
        #         [1,  0, -1, -1],  # O
        #         [2,  1, 0, -1],  # O
        #         [3,  2, 1, 0],  # N
        #         [4,  3, 2, 0]  # O
        #     ])},
        # ).get_embedded_molecule()
        # hoono = Molecule.from_file(
        #     TestManager.test_data('HOONO_freq.fchk'),
        #     internals={'zmatrix': [
        #         [0, -1, -1, -1],  # H
        #         [1,  0, -1, -1],  # O
        #         [2,  1,  0, -1],  # O
        #         [3,  2,  1,  0],  # N
        #         [4,  3,  2,  0]  # O
        #     ], 'iterative': True},
        # ).get_embedded_molecule()
        # hoono2 = Molecule.from_file(
        #     TestManager.test_data('HOONO_freq.fchk'),
        #     internals=[
        #         [0, -1, -1, -1],  # H
        #         [1,  0, -1, -1],  # O
        #         [2,  1,  0, -1],  # O
        #         [3,  2,  1,  0],  # N
        #         [4,  3,  2,  0]  # O
        #     ]
        # ).get_embedded_molecule()
        # # h1 = hoono.get_internals_by_cartesians(1)[0]
        # # h2 = hoono2.get_internals_by_cartesians(1)[0]
        # ho2_ics = hoono2.internal_coordinates.copy()
        # ho2_ics[3, 2] += 1
        # ho2_new_carts = ho2_ics.convert(hoono2.coords.system)
        #
        # print("="*100)
        #
        # ho3_ics = hoono3.internal_coordinates.copy()
        # print(hoono3.internal_coordinates.converter_options['specs'])
        # ho3_ics[5] += 1
        # ho3_new_carts = ho3_ics.convert(hoono.coords.system)
        #
        # print("="*100)
        #
        # ho1_ics = hoono.internal_coordinates.copy()
        # ho1_ics[3, 2] += 1
        # ho1_new_carts = ho1_ics.convert(hoono.coords.system)
        #
        # from McUtils.Numputils import eckart_embedding
        # hoono_og_carts = eckart_embedding(
        #     ho2_new_carts,
        #     hoono.coords,
        #     masses=hoono.masses
        # ).coordinates
        #
        # ho1_new_carts = eckart_embedding(
        #     ho2_new_carts,
        #     ho1_new_carts,
        #     masses=hoono.masses
        # ).coordinates
        #
        # ho3_new_carts = eckart_embedding(
        #     ho2_new_carts,
        #     ho3_new_carts,
        #     masses=hoono.masses
        # ).coordinates
        # print("="*20)
        # print(ho2_new_carts - ho1_new_carts)
        # print(ho2_new_carts - ho3_new_carts)
        # print("_"*20)
        #
        # return

        # ploot = hoono2.plot(hoono_og_carts, backend='x3d')
        # hoono.plot(ho1_new_carts, backend='x3d', figure=ploot,
        #            atom_style={'color': 'white', 'glow': '#00ff00'},
        #            bond_style={'color': 'white', 'glow': '#00ff00'}
        #            )
        # hoono.plot(ho2_new_carts, backend='x3d', figure=ploot,
        #            atom_style={'color': 'white', 'glow': '#f0f'},
        #            bond_style={'color': 'white', 'glow': '#f0f'}
        #            )
        # ploot.show()
        #
        # return
        #
        # ho1_ics = hoono.internal_coordinates.copy()
        # ho1_ics[2, 2] += .1
        # ho1_new_carts = ho1_ics.convert(hoono.coords.system)
        #
        # ploot = hoono2.plot(ho2_new_carts, backend='x3d')
        # hoono.plot(ho1_new_carts, backend='x3d', figure=ploot,
        #            atom_style={'color': 'white', 'glow': '#00ff00'},
        #            bond_style={'color': 'white', 'glow': '#00ff00'}
        #            )
        # hoono.plot(backend='x3d', figure=ploot,
        #            atom_style={'color': 'white', 'glow': '#f0f'},
        #            bond_style={'color': 'white', 'glow': '#f0f'}
        #            )
        # ploot.show()
        #
        # return

        r = Reaction.from_smiles(
            "C=C.C=CC=C>>C1CCC=CC1",
            fragment_expansion_method='centroid',
            min_distance=.4,
            add_radius=False,
            expansion_factor=.01,
            # energy_evaluator='rdkit'
            energy_evaluator='rdkit',
            optimize=True
        )

        # for mol in r.reactant_complex.fragments:
        #     print(mol.coords)
        # rxn_plot = r.reactant_complex.plot(backend='x3d')
        # r.product_complex.plot(backend='x3d', figure=rxn_plot,
        #                        highlight_atoms=list(range(len(r.reactant_complex.atoms))))
        # rxn_plot.show()
        zmat = [
                [0, -1, -1, -1],
                [1,  0, -1, -1],
                [2,  1,  0, -1],
                [3,  2,  1,  0],
                [4,  3,  2,  1],
                [5,  4,  3,  2],
                [6,  0,  5,  1],
                [7,  0,  5,  1],
                [8,  1,  0,  2],
                [9,  1,  0,  2],
                [10, 2,  1,  3],
                [11, 2,  1,  3],
                [12, 3,  2,  4],
                [13, 4,  3,  5],
                [14, 5,  4,  0],
                [15, 5,  4,  0]
            ]
        def get_specs(zmat):
            specs = []
            for row in zmat:
                if row[1] < 0: continue
                specs.append(tuple(row[:2]))
                if row[2] < 0: continue
                specs.append(tuple(row[:3]))
                if row[3] < 0: continue
                specs.append(tuple(row[:4]))
            return specs

        # gen_cart = r.get_profile_generator(
        #     'neb',
        #     energy_evaluator='aimnet2',
        #     # initial_image_positions=[0, .2, .6, .7, .8, .9, .95, 1],
        #     # initial_image_positions=[0, .1, 1],
        #     num_images=25,
        #     # interpolation_gradient_scaling=.3,
        #     # internals=zmat,
        #     # internals={'zmatrix':zmat, 'iterative':True},
        #     # internals={'specs':get_specs(zmat)},
        #     # internals={'primitives':get_specs(zmat)}
        #     # internals={'specs':'auto', 'prune_coordinates':False},
        #     # internals={'primitives':'auto', 'prune_coordinates':False},
        #     spring_constant=1
        # )
        #
        # gen = r.get_profile_generator(
        #     'neb',
        #     energy_evaluator='aimnet2',
        #     initial_image_positions=[0, .2, .6, .7, .8, .9, .95, 1],
        #     # initial_image_positions=[0, .1, 1],
        #     num_images=25,
        #     # interpolation_gradient_scaling=.3,
        #     # internals=zmat,
        #     internals={'zmatrix':zmat, 'iterative':True},
        #     # internals={'specs':get_specs(zmat)},
        #     # internals={'primitives':get_specs(zmat)}
        #     # internals={'specs':'auto', 'prune_coordinates':False},
        #     # internals={'primitives':'auto', 'prune_coordinates':False},
        #     spring_constant=1
        # )

        # gen_gen = r.get_profile_generator(
        #     'neb',
        #     energy_evaluator='aimnet2',
        #     # initial_image_positions=[0, .2, .6, .7, .8, .9, .95, 1],
        #     # initial_image_positions=[0, .1, 1],
        #     num_images=8,
        #     # interpolation_gradient_scaling=.3,
        #     # internals=zmat,
        #     # internals={'zmatrix': zmat, 'iterative': True},
        #     # internals={'specs':get_specs(zmat)},
        #     # internals={'primitives':get_specs(zmat)}
        #     internals={'specs':'auto', 'prune_coordinates':False},
        #     # internals={'primitives':'auto', 'prune_coordinates':False},
        #     spring_constant=1
        # )

        # int_vals = [0, .2,  .4, .6, .8, 1]
        # crds = gen.interpolator.interpolate(int_vals)
        # img = [
        #     r.product_complex.modify(coords=struct)
        #         if p > .7 else
        #     r.reactant_complex.modify(coords=struct)
        #         for p,struct in zip(int_vals, crds)
        #     ]

        from McUtils.Data import UnitsData

        generators = [
            r.get_profile_generator(
                'neb',
                energy_evaluator='aimnet2',
                # initial_image_positions=[0, .2, .4, .6, .625, .65, .675, .7, .725, .75, .775, .8, .9, 1],
                # initial_image_positions=[0, .1, 1],
                # initial_image_positions = np.concatenate([
                #     np.arange(0, .6, .1),
                #     np.arange(.6, .8, .05),
                #     np.arange(.8, 1, .1),
                #     [.1]
                # ]),
                num_images=50,
                # interpolation_gradient_scaling=.3,
                internals=coords,
                # internals=zmat,
                # internals={'zmatrix': zmat, 'iterative': True},
                # internals={'specs':get_specs(zmat)},
                # internals={'primitives':get_specs(zmat)}
                # internals={'specs':'auto', 'prune_coordinates':False},
                # internals={'primitives':'auto', 'prune_coordinates':False},
                spring_constant=.5
            )
            for coords in [
                # {'zmatrix': zmat, 'iterative': True},
                zmat,
                None,
                # zmat
                # {'specs':'auto', 'prune_coordinates':False}
            ]
        ]

        profiles = []
        data = []
        for g in generators:
            pre, post = g.generate(return_preopt=True,
                                   max_iterations=25,
                                   max_displacement=.1
                                   )
            profiles.append([pre, post])
            pre_eng = g.evaluate_profile_energies(pre)
            post_eng = g.evaluate_profile_energies(post)
            pre_dist = g.evaluate_profile_distances(pre)
            post_dist = g.evaluate_profile_distances(post)
            data.append([pre_dist,
                         (pre_eng - pre_eng[0]) * UnitsData.convert("Hartrees", "Kilojoules/Mole"),
                         post_dist,
                         (post_eng - pre_eng[0]) * UnitsData.convert("Hartrees", "Kilojoules/Mole")])

        import McUtils.Plots as plt
        eng_plot = None
        colors = ["#5e81b5", "#e19c24", "#8fb032"]
        for i, (pre_dist, pre_eng, post_dist, post_eng) in enumerate(data):
            eng_plot = plt.Plot(pre_dist, pre_eng, figure=eng_plot, marker='.', linestyle='dashed', color=colors[i % len(colors)],
                                axes_labels=['Reaction Coordinate', "Energy (kJ mol$^{-1}$)"],
                                padding=[[80, 20], [60, 0]],
                                image_size=[350, 350]
                                )
            eng_plot = plt.Plot(post_dist, post_eng, figure=eng_plot, marker='.', color=colors[i % len(colors)])
        eng_plot.savefig(os.path.expanduser("~/Desktop/comp_plot_expanded.pdf"))
        eng_plot.show()

        return

        img = profiles[0][1]

        from McUtils.Jupyter import ExamplesManager
        vm = ExamplesManager.parse_x3d_view_matrix("{\"_00\":-0.2344741117584002,\"_01\":-0.9704258662186166,\"_02\":-0.05740669899034502,\"_03\":-0.5905212871236172,\"_10\":0.927679505284065,\"_11\":-0.2057128258010835,\"_12\":-0.3115974466789714,\"_13\":-0.3163106911448954,\"_20\":0.2905729278357282,\"_21\":-0.12631655265559674,\"_22\":0.948478519595548,\"_23\":-14.848572311368708,\"_30\":0,\"_31\":0,\"_32\":0,\"_33\":1.0000000000000002}")

        # start = r.reactant_complex.modify(coords=crds[0])
        rxn_plot = img[0].plot(backend='x3d', image_size=(800, 800),
                              view_settings=vm
                              )
        rxn_plot.show()
        return

        colors = ([
            '#00ff00',
            '#ff00ff',
            '#0000ff',
            '#dddd00',
            '#00cc00',
            '#cc00cc',
            '#0000cc',
            '#f0f000',
        ]*(len(img)//8 + 1))[:len(img) - 2]
        for mid, color in zip(
                img[1:],
                colors + [None]
        ):
            mid.plot(backend='x3d',
                     figure=rxn_plot,
                     atom_style={'color': 'white', 'glow': color} if color is not None else None,
                     bond_style={'color': 'white', 'glow': color} if color is not None else None
                     )
        rxn_plot.show()
