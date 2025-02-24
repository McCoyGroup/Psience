from Peeves.TestUtils import *
from unittest import TestCase
from Peeves import BlockProfiler

import numpy as np
from Psience.Reactions import Reaction

class ReactionsTests(TestCase):

    @debugTest
    def test_SetupReactants(self):

        r = Reaction.from_smiles(
            "C=C.C=CC=C>>C1CCC=CC1",
            fragment_expansion_method='centroid',
            # energy_evaluator='rdkit'
            energy_evaluator='aimnet2'
        )
        # for mol in r.reactant_complex.fragments:
        #     print(mol.coords)
        # rxn_plot = r.reactant_complex.plot(backend='x3d')
        # r.product_complex.plot(backend='x3d', figure=rxn_plot,
        #                        highlight_atoms=list(range(len(r.reactant_complex.atoms))))
        # rxn_plot.show()
        zmat = [
                [0, -1, -1, -1],
                [1, 0, -1, -1],
                [2, 1, 0, -1],
                [3, 2, 1, 0],
                [4, 3, 2, 1],
                [5, 4, 3, 2],
                [6, 0, 5, 1],
                [7, 0, 5, 1],
                [8, 1, 0, 2],
                [9, 1, 0, 2],
                [10, 2, 1, 3],
                [11, 2, 1, 3],
                [12, 3, 2, 4],
                [13, 4, 3, 5],
                [14, 5, 4, 0],
                [15, 5, 4, 0]
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

        gen = r.get_profile_generator(
            'neb',
            energy_evaluator='aimnet2',
            initial_image_positions=[0,  .4, .6, .7, .8, .9, .95, 1],
            # internals=zmat,
            # internals={'specs':get_specs(zmat)},
            # internals={'primitives':get_specs(zmat)}
            # internals={'specs':'auto', 'prune_coordinates':False}
            internals={'primitives':'auto', 'prune_coordinates':False},
            spring_constant=.5
        )
        # int_vals = [0, .2,  .4, .6, .8, 1]
        # crds = gen.interpolator.interpolate(int_vals)
        # img = [
        #     r.product_complex.modify(coords=struct)
        #         if p > .7 else
        #     r.reactant_complex.modify(coords=struct)
        #         for p,struct in zip(int_vals, crds)
        #     ]

        pre, post = gen.generate(return_preopt=True, max_iterations=15)
        # return

        pre_eng = gen.evaluate_profile_energies(pre)
        post_eng = gen.evaluate_profile_energies(post)
        pre_dist = gen.evaluate_profile_distances(post)
        post_dist = gen.evaluate_profile_distances(post)

        import McUtils.Plots as plt
        eng_plot = plt.Plot(pre_dist, pre_eng, marker='.')
        eng_plot = plt.Plot(post_dist, post_eng, figure=eng_plot, marker='.')
        eng_plot.show()
        # return

        img = post

        from McUtils.Jupyter import ExamplesManager
        vm = ExamplesManager.parse_x3d_view_matrix("{\"_00\":-0.2344741117584002,\"_01\":-0.9704258662186166,\"_02\":-0.05740669899034502,\"_03\":-0.5905212871236172,\"_10\":0.927679505284065,\"_11\":-0.2057128258010835,\"_12\":-0.3115974466789714,\"_13\":-0.3163106911448954,\"_20\":0.2905729278357282,\"_21\":-0.12631655265559674,\"_22\":0.948478519595548,\"_23\":-14.848572311368708,\"_30\":0,\"_31\":0,\"_32\":0,\"_33\":1.0000000000000002}")

        # start = r.reactant_complex.modify(coords=crds[0])
        rxn_plot = img[0].plot(backend='x3d', image_size=(800, 800),
                              view_settings=vm
                              )
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
