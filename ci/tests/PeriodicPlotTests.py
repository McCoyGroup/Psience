"""Directional periodicity and native plotting overlays for PeriodicSystem."""
import copy
import io
import unittest

import numpy as np
from numpy.testing import assert_allclose

from McUtils.Data import UnitsData
from McUtils.Plots import Arrow, Cylinder, Line, Text
import McUtils.Numputils as nput
from Psience.Molecools import Molecule, PeriodicSystem
from PeriodicSystemTests import CIF


class PeriodicPlotTests(unittest.TestCase):
    def silicon(self):
        fcc = np.array([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
        frac = np.concatenate([fcc, fcc + .25])
        lattice = np.eye(3) * 5.431 * UnitsData.convert("Angstroms", "BohrRadius")
        return PeriodicSystem(Molecule(["Si"] * 8, frac @ lattice), lattice)

    def system(self, periodicity=(True, True, True)):
        return PeriodicSystem(Molecule(["C"], [[.2, .3, .4]], guess_bonds=False),
                              np.diag([4., 5., 6.]), periodicity=periodicity)

    def plot(self, system, **options):
        return system.plot(bonds=False, rings=False, objects=True, return_objects=True, units=None, **options)

    def test_cell_and_lattice_arrows_follow_periodic_dimensions(self):
        for mask, edge_count, labels in (((True, False, False), 1, "a"),
                                          ((True, True, False), 4, "ab"),
                                          ((True, True, True), 12, "abc"),
                                          ((True, False, True), 4, "ac")):
            with self.subTest(mask=mask):
                p = self.system(mask)
                _, atoms, bonds, annotations = self.plot(p, unit_axes=True)
                primitives = annotations[0]
                self.assertEqual(sum(isinstance(o, Line) for o in primitives), edge_count)
                self.assertEqual("".join(o.txt for o in primitives if isinstance(o, Text)), labels)
                arrows = [o for o in primitives if isinstance(o, Arrow)]
                self.assertEqual(len(arrows), sum(mask))
                for arrow, i in zip(arrows, np.flatnonzero(mask)):
                    assert_allclose(arrow.pos1, p.origin)
                    assert_allclose(arrow.pos2, p.origin + p.lattice[i])

    def test_cell_direction_styles_cylinders_and_axes_display_units(self):
        p = self.system()
        style = {"radius": .025, "transparency": .25, "a": {"color": "red"},
                 "b": {"color": "green"}, "c": {"color": "blue"}}
        original = copy.deepcopy(style)
        _, _, _, annotations = p.plot(bonds=False, rings=False, objects=True, return_objects=True,
                                     cell_style=style, unit_axes=True, unit_axes_origin=[1., 2., 3.],
                                     unit_axes_scale=.75, unit_axes_style={"radius": .04, "a": {"color": "purple"}},
                                     unit_axes_label_style={"a": {"font_size": 18}})
        edges = [o for o in annotations[0] if isinstance(o, Cylinder)]
        self.assertEqual(len(edges), 12)
        self.assertEqual(sorted(o.opts["color"] for o in edges), ["blue"] * 4 + ["green"] * 4 + ["red"] * 4)
        self.assertTrue(all(o.opts["transparency"] == .25 for o in edges))
        arrows = [o for o in annotations[0] if isinstance(o, Arrow)]
        factor = UnitsData.convert("BohrRadius", "Angstroms")
        assert_allclose(arrows[0].pos1, np.array([1., 2., 3.]) * factor)
        assert_allclose(arrows[0].pos2, (np.array([1., 2., 3.]) + .75 * p.lattice[0]) * factor)
        self.assertEqual(arrows[0].opts["color"], "purple")
        self.assertEqual([o for o in annotations[0] if isinstance(o, Text)][0].opts["font_size"], 18)
        self.assertEqual(style, original)

    def test_wire_tiling_and_wrapping_preserve_finite_coordinates(self):
        p = PeriodicSystem(Molecule(["C"], [[4.2, 8., -3.]], guess_bonds=False),
                           np.diag([4., 5., 6.]), periodicity=(True, False, False))
        tiled = p.tile(2)
        assert_allclose(tiled.molecule.coords, [[.2, 8., -3.], [4.2, 8., -3.]])
        assert_allclose(tiled.lattice, np.diag([8., 5., 6.]))
        assert_allclose(tiled.periodicity, [True, False, False])
        assert_allclose(tiled.copy().periodicity, p.periodicity)
        self.assertEqual(len(tiled.reduce_symmetry().molecule.atoms), 1)
        with self.assertRaisesRegex(ValueError, "nonperiodic"):
            p.tile((2, 2, 1))

    def test_slab_expansion_does_not_identify_sites_separated_in_finite_direction(self):
        molecule = Molecule(["C", "C"], [[.2, .3, .4], [.2, .3, 6.4]], guess_bonds=False)
        p = PeriodicSystem(molecule, np.diag([4., 5., 6.]), periodicity=(True, True, False), symmetry_reduced=True)
        full = p.expand_symmetry()
        self.assertEqual(len(full.molecule.atoms), 2)
        assert_allclose(np.asarray(full.molecule.coords), np.asarray(molecule.coords))
        self.assertEqual(len(full.reduce_symmetry().molecule.atoms), 2)

    def test_finite_symmetry_translation_is_not_reduced_modulo_one(self):
        mirror = np.eye(4)
        mirror[2, 2] = -1
        mirror[2, 3] = 4
        p = PeriodicSystem(Molecule(["C"], [[.2, .3, 2.3]], guess_bonds=False), np.eye(3),
                           periodicity=(True, True, False), symmetry_reduced=True, symmetry_operations=[mirror])
        full = p.expand_symmetry()
        assert_allclose(full.molecule.coords, [[.2, .3, 2.3], [.2, .3, 1.7]])
        assert_allclose(np.asarray(full.reduce_symmetry().expand_symmetry().molecule.coords), np.asarray(full.molecule.coords))

    def test_cif_periodicity_and_incompatible_symmetry(self):
        p = PeriodicSystem.from_cif(io.StringIO(CIF), guess_bonds=False, periodicity=(True, True, False))
        assert_allclose(p.fractional_coordinates, [[.1, .2, .3], [.9, .8, -.3], [0, 0, 0]], atol=1e-14)
        self.assertEqual(len(p.reduce_symmetry().molecule.atoms), 2)
        swap = np.eye(4)
        swap[:3, :3] = [[0, 1, 0], [1, 0, 0], [0, 0, 1]]
        with self.assertRaisesRegex(ValueError, "periodic directions"):
            PeriodicSystem(Molecule(["C"], [[0., 0., 0.]], guess_bonds=False), np.eye(3),
                           periodicity=(True, False, True), symmetry_operations=[swap])

    def test_default_images_are_outside_central_axes_with_half_transparency(self):
        p = self.system()
        extra = Line([[0, 0, 0], [1, 0, 0]])
        _, atoms, _, annotations = self.plot(p, show_periodic_images=True, unit_axes=True,
                                             annotation_function=lambda *args: [extra])
        self.assertEqual(len(atoms), 27)
        self.assertNotIn("transparency", atoms[0][0].opts)
        self.assertTrue(all(a[0].opts["transparency"] == .5 for a in atoms[1:]))
        positions = np.array([a[0].pos for a in atoms])
        assert_allclose(positions.min(axis=0), np.asarray(p.molecule.coords)[0] - p.lattice.sum(axis=0))
        assert_allclose(positions.max(axis=0), np.asarray(p.molecule.coords)[0] + p.lattice.sum(axis=0))
        self.assertEqual(len(annotations[0]), 19)  # twelve edges, three arrows/labels, user annotation
        self.assertTrue(all(not group for group in annotations[1:]))
        self.assertIs(annotations[0][-1], extra)
        self.assertEqual(len(p.molecule.atoms), 1)

    def test_boolean_and_integer_image_counts_and_explicit_disable(self):
        p = self.system()
        for counts, expected in (([True, False, True], 9), ([2, 0, 1], 15), ([0, 0, 0], 1)):
            with self.subTest(counts=counts):
                _, atoms, _, _ = self.plot(p, periodic_images=counts)
                self.assertEqual(len(atoms), expected)
        _, atoms, _, _ = self.plot(p, periodic_images=[2, 0, 1], show_periodic_images=False)
        self.assertEqual(len(atoms), 1)

    def test_asymmetric_image_pairs_set_negative_and_positive_layers(self):
        p = self.system()
        counts = [2, 0, [1, 0]]
        original = copy.deepcopy(counts)
        _, atoms, _, annotations = self.plot(p, periodic_images=counts, unit_axes=True)
        self.assertEqual(len(atoms), 10)
        centers = np.array([a[0].pos for a in atoms])
        fractions = np.linalg.solve(p.lattice.T, (centers - np.asarray(p.molecule.coords)[0]).T).T
        actual = {tuple(np.round(f).astype(int)) for f in fractions}
        expected = {(a, 0, c) for a in range(-2, 3) for c in (-1, 0)}
        self.assertEqual(actual, expected)
        self.assertEqual(counts, original)
        self.assertTrue(all(image[0].opts["transparency"] == .5 for image in atoms[1:]))
        self.assertTrue(all(not group for group in annotations[1:]))
        for counts, expected in (([[0, 2], 0, 0], 3),
                                 ([[1, 2], [0, 1], [2, 0]], 24),
                                 ([True, False, [True, False]], 6),
                                 (np.array([[0, 1], [1, 0], [0, 0]]), 4)):
            with self.subTest(counts=counts):
                _, atoms, _, _ = self.plot(p, periodic_images=counts)
                self.assertEqual(len(atoms), expected)

    def test_asymmetric_counts_respect_finite_directions_and_validate_pairs(self):
        p = self.system((True, False, True))
        _, atoms, _, _ = self.plot(p, periodic_images=[[0, 1], [3, 4], [1, 0]])
        self.assertEqual(len(atoms), 4)
        for counts in ([2, 0, [-1, 0]], [2, 0, [1.5, 0]], [2, 0, [1]],
                       [2, 0, [1, 0, 0]], [2, 0, [[1, 0]]]):
            with self.subTest(counts=counts), self.assertRaisesRegex(ValueError, "count pairs"):
                self.plot(p, periodic_images=counts)

    def test_images_never_translate_along_finite_directions(self):
        for mask, expected in (((True, False, False), 3), ((True, True, False), 9), (1, 3), (2, 9), (3, 27)):
            p = self.system(mask)
            _, atoms, _, _ = self.plot(p, show_periodic_images=True)
            self.assertEqual(len(atoms), expected)
            positions = np.array([a[0].pos for a in atoms])
            for i in np.flatnonzero(~p.periodicity):
                assert_allclose(positions[:, i], np.asarray(p.molecule.coords)[0, i])

    def test_image_material_styles_preserve_multiple_bonds_and_user_inputs(self):
        molecule = Molecule(["C", "O"], [[.2, .3, .4], [1.2, .3, .4]], bonds=[[0, 1, 2]], guess_bonds=False)
        p = PeriodicSystem(molecule, np.diag([4., 5., 6.]))
        style = {"transparency": .7, "color": "cyan", "atom_style": {"C": {"transparency": .25}}}
        original = copy.deepcopy(style)
        _, atoms, bonds, _ = p.plot(rings=False, objects=True, return_objects=True, units=None,
                                   periodic_images=[1, 0, 0], periodic_image_styles=style,
                                   atom_style={"transparency": 0., "O": {"color": "red"}},
                                   bond_style={"transparency": 0.})
        self.assertEqual(len(atoms), 3)
        self.assertEqual(atoms[0][1].opts["color"], "red")
        for image in atoms[1:]:
            self.assertEqual(image[0].opts["transparency"], .25)
            self.assertEqual(image[1].opts["transparency"], .7)
            self.assertTrue(all(o.opts["color"] == "cyan" for o in image))
        self.assertEqual(len(bonds[1][0]), len(bonds[0][0]))
        self.assertGreater(len(bonds[1][0]), 2)
        self.assertTrue(all(o.opts["transparency"] == .7 for image in bonds[1:] for bond in image for o in bond))
        self.assertEqual(style, original)

    def test_actual_native_x3d_axes_and_images_export(self):
        p = self.system((True, False, True))
        fig = p.plot(bonds=False, rings=False, unit_axes=True, periodic_images=[1, 0, 1],
                     cell_style={"radius": .025, "color": "navy", "transparency": .2})
        html = fig.figure.to_x3d().to_html().tostring()
        self.assertEqual(html.count("<Sphere "), 9)
        self.assertEqual(html.lower().count("<cone "), 2)
        self.assertIn('transparency="0.500"', html)
        self.assertIn('string="a"', html)
        self.assertIn('string="c"', html)

    def test_lattice_and_principal_axes_coexist_without_copying_them_to_images(self):
        p = PeriodicSystem.from_cif(io.StringIO(CIF), guess_bonds=False)
        _, atoms, _, annotations = self.plot(p, unit_axes=True, principle_axes=True, periodic_images=[1, 0, 0])
        self.assertEqual(sum(isinstance(o, Arrow) for o in annotations[0]), 6)
        self.assertEqual("".join(o.txt for o in annotations[0] if isinstance(o, Text)), "abc")
        self.assertTrue(all(not group for group in annotations[1:]))
        _, atoms, _, _ = self.plot(p, atom_style=False, periodic_images=[1, 0, 0])
        self.assertTrue(all(not group for group in atoms))

    def test_silicon_image_network_matches_full_supercell_and_nearest_neighbors(self):
        p = self.silicon()
        original = np.asarray(p.molecule.coords).copy()
        network, bridges, cells = p._image_bond_network(p._image_translations([2, 0, [1, 0]]), None)
        distances = nput.distance_matrix(np.asarray(network.coords))
        nearest = np.sqrt(3) * p.lattice[0, 0] / 4
        expected = set(map(tuple, np.argwhere(np.triu(np.isclose(distances, nearest), 1))))
        actual = {tuple(sorted(b[:2])) for b in network.bonds}
        self.assertEqual(actual, expected)
        self.assertEqual(len(actual), 100)
        self.assertEqual(len(bridges), 30)
        self.assertEqual(sum(cells[i] == 0 or cells[j] == 0 for i, j, *_ in bridges), 7)
        _, atoms, bonds, _ = p.plot(periodic_images=[2, 0, [1, 0]], rings=False,
                                     objects=True, return_objects=True)
        self.assertEqual(sum(len(group or []) for group in bonds), 100)
        self.assertEqual(sum(len(group or []) for group in atoms), 80)
        tiled = p.tile((5, 1, 2))
        offset = -2 * p.lattice[0] - p.lattice[2]

        def physical_edges(molecule, shift):
            coords = np.asarray(molecule.coords) + shift
            return {tuple(sorted((tuple(np.round(coords[i], 8)), tuple(np.round(coords[j], 8)))))
                    for i, j, *_ in molecule.bonds}

        self.assertEqual(physical_edges(tiled.molecule, offset), physical_edges(network, 0))
        self.assertEqual(len(p.molecule.bonds), 7)
        assert_allclose(np.asarray(p.molecule.coords), original)

    def test_silicon_boundary_bond_halves_follow_cell_materials(self):
        p = self.silicon()
        _, _, bonds, _ = p.plot(periodic_images=[2, 0, [1, 0]], rings=False,
                                 objects=True, return_objects=True,
                                 bond_style={"transparency": .1, "color": "orange"},
                                 periodic_image_styles={"transparency": .6, "color": "cyan"})
        styles = [tuple(part.opts.get("transparency", 0) for part in bond)
                  for group in bonds for bond in (group or [])]
        self.assertEqual(styles.count((.1, .1)), 7)
        self.assertEqual(styles.count((.1, .6)), 7)
        self.assertEqual(styles.count((.6, .6)), 86)
        mixed = [bond for group in bonds for bond in (group or [])
                 if tuple(part.opts.get("transparency", 0) for part in bond) == (.1, .6)]
        self.assertTrue(all(tuple(part.opts["color"] for part in bond) == ("orange", "cyan") for bond in mixed))

    def test_bond_inference_and_rendering_can_be_disabled(self):
        p = self.silicon()
        p.molecule.bonds = p.molecule.get_guessed_bonds()
        p.molecule.guess_bonds = False
        _, _, bonds, _ = p.plot(periodic_images=[2, 0, [1, 0]], rings=False, objects=True, return_objects=True)
        self.assertEqual(sum(len(group or []) for group in bonds), 70)
        self.assertEqual(len(p.tile((5, 1, 2)).molecule.bonds), 70)
        p.molecule.guess_bonds = True
        _, _, bonds, _ = p.plot(periodic_images=[2, 0, [1, 0]], bonds=False, rings=False,
                                 objects=True, return_objects=True)
        self.assertTrue(all(not group for group in bonds))
        _, _, bonds, _ = p.plot(periodic_images=[2, 0, [1, 0]], rings=False,
                                 periodic_image_styles={"bond_style": False}, objects=True, return_objects=True)
        self.assertEqual(sum(len(group or []) for group in bonds), 7)

    def test_wire_connects_adjacent_cells_and_preserves_supplied_orders(self):
        factor = UnitsData.convert("Angstroms", "BohrRadius")
        molecule = Molecule(["C", "C"], np.array([[.5, 0, 0], [1.92, 0, 0]]) * factor,
                            bonds=[[0, 1, 2]])
        p = PeriodicSystem(molecule, np.diag([2.84, 8., 8.]) * factor, periodicity=1)
        tiled = p.tile(3)
        self.assertEqual(len(tiled.molecule.bonds), 5)
        for offset in (0, 2, 4):
            self.assertIn([offset, offset + 1, 2], tiled.molecule.bonds)
        _, _, bonds, _ = p.plot(show_periodic_images=True, rings=False, objects=True, return_objects=True)
        self.assertEqual(sum(len(group or []) for group in bonds), 5)

    def test_new_plot_validation(self):
        p = self.system()
        for kwargs in ({"periodic_images": [-1, 0, 1]}, {"periodic_images": [1.2, 0, 0]},
                       {"periodic_images": [1, 0]}, {"unit_axes": [True, False]},
                       {"unit_axes": True, "unit_axes_scale": 0}, {"cell_style": {"radius": -1}},
                       {"cell_style": [{}, {}]}, {"show_periodic_images": "yes"}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                self.plot(p, **kwargs)
        with self.assertRaisesRegex(ValueError, "periodicity"):
            PeriodicSystem(p.molecule, p.lattice, periodicity=[True, False])
        with self.assertRaisesRegex(ValueError, "3D"):
            p.plot(unit_axes=True, backend="jsmol")
        with self.assertRaisesRegex(ValueError, "static"):
            p.plot(show_periodic_images=True, animate=True)


if __name__ == "__main__":
    unittest.main()
