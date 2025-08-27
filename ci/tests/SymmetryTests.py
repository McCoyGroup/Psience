# import itertools
#
# import scipy.linalg
#
# import McUtils.Zachary
from Peeves.TestUtils import *
from unittest import TestCase

# from McUtils.Data import UnitsData, PotentialData
# from McUtils.Zachary import Interpolator, Symbols
# import McUtils.Plots as plt
# from McUtils.GaussianInterface import GaussianLogReader

from Psience.Symmetry import *
from Psience.Molecools import Molecule
import numpy as np
import McUtils.Numputils as nput

class SymmetryTests(TestCase):

    @validationTest
    def test_SymmetryElements(self):

        for s in [
            InversionElement(),
            RotationElement(5, [0, 1, 1], 2),
            ImproperRotationElement(11, [1, 0, 0]),
            ReflectionElement([1, 0, 1])
        ]:
            self.assertEquals(
                SymmetryElement.from_transformation_matrix(
                    s.get_transformation()
                ),
                s
            )

    @validationTest
    def test_Transformations(self):

        np.random.seed(12223)
        tf = nput.view_matrix(np.random.rand(3))
        ax = np.random.rand(3)
        for a in [
            InversionElement(),
            RotationElement(2, ax),
            RotationElement(5, ax),
            RotationElement(7, ax, 2),
            RotationElement(22, ax, 11),
            ReflectionElement(ax),
            ImproperRotationElement(3, ax),
            ImproperRotationElement(6, ax, 5),
            ImproperRotationElement(7, ax, 5),
        ]:
            s = a.transform(tf)

            x1 = s.get_transformation()
            x2 = tf @ a.get_transformation() @ tf.T
            # print(x1)
            # print(x2)

            self.assertTrue(
                np.allclose(x1, x2)
            )

    @validationTest
    def test_Composition(self):
        np.random.seed(123)
        ax = np.random.rand(3)
        ax2 = np.cross(ax, [0, 1, 0])
        ax3 = RotationElement(6, ax).get_transformation() @ ax
        for a,b in [
            [InversionElement(), RotationElement(7, ax2, root=5)],
            [RotationElement(2, ax), RotationElement(7, ax2, root=5)],
            [RotationElement(2, ax), ImproperRotationElement(7, ax2, root=5)],
            [RotationElement(3, ax), RotationElement(7, ax, root=5)],
            [RotationElement(6, ax), RotationElement(11, ax3)],
            [RotationElement(6, ax), ImproperRotationElement(12, ax3)],
        ]:
            s1 = a @ b
            s2 = SymmetryElement.compose(a, b)
            print(a, "@", b, "=>", s1)

            if hasattr(s1, 'bits'):
                print(":", SymmetryElement.from_transformation_matrix(s1.get_transformation()))


            self.assertTrue(
                np.allclose(
                    s1.get_transformation(),
                    s2.get_transformation()
                )
            )

    @validationTest
    def test_PointGroupElements(self):
        for name in [
            # ("Cv", 2),
            # ("S", 4),
            # ("Ch", 4),
            # ("Ch", 6),
            # ("Dh", 4),
            # ("Dh", 6),
            # ("Dd", 7),
            # ("D", 7),
            ("T",),
            ("Td",),
            ("Th",),
            ("O",),
            ("Oh",),
            ("I",),
            ("Ih",),
        ]:
            pg = PointGroup.from_name(*name)
            ct = pg.character_table
            elements = ct.permutations
            classes = ct.classes
            print(pg, pg.elements)
            self.assertIsNot(classes, None)
            self.assertEquals(
                np.sort(np.concatenate(classes)).tolist(),
                np.arange(sum(len(l) for l in classes)).tolist()
            )
            self.assertEquals(len(elements), sum(len(l) for l in classes))


    @validationTest
    def test_PointGroupAlignments(self):
        pg = PointGroup.from_name("Ch", 6)
        print(pg.axes @ pg.axes.T)
        new_pg = pg.align(np.eye(3))
        print(new_pg.get_axes())

    @validationTest
    def test_Visualization(self):
        pg = PointGroup.from_name("Dd", 3)
        # base = pg.plot() #.show()
        print(pg.axes)
        new_pg = pg.align(np.eye(3))
        print(new_pg.axes)
        new_pg.plot().show()
        for e in pg.elements:
            if hasattr(e, 'axis'):
                print(e, e.axis)

    # @debugTest
    # def test_SymmetryReduction(self):
    #     perms = nput.permutation_indices(5, 5)[:5]
    #     print(perms)
    #     cycles = nput.permutation_cycles(perms, return_groups=True)
    #     print(cycles)
    #     return
    #
    #     tf = nput.rotation_matrix([0, 0, 1], 2*np.pi/3)
    #     np.random.seed(5)
    #     coords = np.round(np.random.rand(2, 3), 3)
    #     coords = np.concatenate([coords, coords@tf, coords@tf@tf], axis=0)
    #     print(nput.symmetry_permutation(coords, tf))