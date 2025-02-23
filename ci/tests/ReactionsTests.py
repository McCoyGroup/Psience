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
            fragment_expansion_method='centroid'
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

        """
        ==> [[[ 0.77190252  0.         -0.0719183  ...  0.          0.
    0.04552509]
  [-0.16557176  0.         -0.38929849 ...  0.          0.
   -0.01881451]
  [ 0.61380168  0.         -0.01456972 ...  0.          0.
    0.08105092]
  ...
  [ 0.          0.          0.         ... -0.55669355 -0.40498602
    0.00877491]
  [ 0.          0.          0.         ... -0.7321073   0.24340891
    0.26671989]
  [ 0.          0.          0.         ... -0.39257     0.12036494
   -0.50985179]]]
        """
        import McUtils.Numputils as nput
        import McUtils.McUtils.Numputils.CoordOps as coops
        import itertools
        coords = r.reactant_complex.coords
        # for c in itertools.combinations(range(coords.shape[0]), 4):
        #     for p in itertools.permutations(c):
        #         coops.fast_proj = True
        #         new = nput.dihed_vec(coords, *p, order=1)
        #         coops.fast_proj = False
        #         old = nput.dihed_vec(coords, *p, order=1)#, method='classic')
        #         if not np.allclose(new[1], old[1]):
        #             raise ValueError(
        #                 p, new[1], old[1]
        #             )
        # coops.fast_proj = True
        # new = nput.angle_vec(coords, 15, 5, 4, angle_ordering='ijk', order=1)
        # coops.fast_proj = False
        # old = nput.angle_vec(coords, 15, 5, 4, angle_ordering='ijk', order=1)
        # print(new[1][-3:])
        # print(old[1][-3:])
        # return

        gen = r.get_profile_generator(
            'interpolate',
            # internals='auto'
            internals=zmat
            # internals={'specs':get_specs(zmat)}
        )

        # from McUtils.Profilers import BlockProfiler
        # with BlockProfiler.profiler(mode='deterministic'):
        #     for c in itertools.combinations(range(coords.shape[0]), 3):
        #         for p in itertools.permutations(c):
        #             new = nput.angle_vec(coords, *p, order=2)
        # # structs = gen.generate()
        # return
        int_vals = [0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1]
        # int_vals = [0, .1, 1]
        crds = gen.interpolator.interpolator(int_vals)
        # return

        start = r.reactant_complex.modify(coords=crds[0])
        rxn_plot = start.plot(backend='x3d')
        colors = [
            '#00ff00',
            '#ff00ff',
            '#0000ff',
            '#dddd00',
            '#00cc00',
            '#cc00cc',
            '#0000cc',
            '#f0f000',
        ][:len(crds) - 2]
        for struct, color, p in zip(
                crds[1:],
                colors + [None],
                int_vals
        ):
            if p > .7:
                mid = r.product_complex.modify(coords=struct)
            else:
                mid = r.reactant_complex.modify(coords=struct)

            mid.plot(backend='x3d',
                     figure=rxn_plot,
                     atom_style={'color': 'white', 'glow': color} if color is not None else None,
                     bond_style={'color': 'white', 'glow': color} if color is not None else None
                     )
        rxn_plot.show()
