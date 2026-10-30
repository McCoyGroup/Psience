"""Regressions for periodic geometry, CIF import, and native molecular plotting."""
import io
import os
import tempfile
import unittest

import numpy as np
from numpy.testing import assert_allclose

from McUtils.Data import UnitsData
from McUtils.Jupyter.Readouts import cell_matrix
from McUtils.Plots import Line
from Psience.Molecools import Molecule, PeriodicSystem


CIF = """data_inversion
_cell_length_a 4.0(2)
_cell_length_b 5.0
_cell_length_c 6.0
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_space_group_name_H-M_alt 'P -1'
_space_group_IT_number 2
loop_
_space_group_symop_id
_space_group_symop_operation_xyz
1 'x, y, z'
2 '-x, -y, -z'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
C1 C 0.1(1) 0.2 0.3 1
O1 O 0 0 0 .5
"""


def inversion():
    operation = np.eye(4)
    operation[:3, :3] *= -1
    return [np.eye(4), operation]


def records(system):
    return sorted((a, *np.round(f % 1, 8), round(o, 8)) for a, f, o in zip(
        system.molecule.atoms, system.fractional_coordinates, system.occupancies))


class PeriodicSystemTests(unittest.TestCase):
    def read(self, text=CIF, **options):
        return PeriodicSystem.from_cif(io.StringIO(text), guess_bonds=False, **options)

    def test_wrapper_units_origin_and_independent_copy(self):
        scale = UnitsData.convert("Angstroms", "BohrRadius")
        lattice = np.diag([4., 5., 6.])
        origin = np.array([1., 2., 3.])
        frac = np.array([[.1, .2, .3]])
        molecule = Molecule(["C"], (origin + frac @ lattice) * scale, guess_bonds=False)
        p = PeriodicSystem(molecule, lattice, origin=origin, lattice_units="Angstroms")
        self.assertIs(p.molecule, molecule)
        self.assertNotIsInstance(p, Molecule)
        self.assertFalse(p.symmetry_reduced)
        assert_allclose(p.fractional_coordinates, frac)
        assert_allclose(p.cell_vertices[[0, -1]], [origin * scale, (origin + lattice.sum(axis=0)) * scale])
        copied = p.copy()
        copied.lattice[0, 0] *= 2
        copied.occupancies[0] = .5
        copied.molecule.coords[0, 0] += 1
        assert_allclose(p.lattice, lattice * scale)
        assert_allclose(p.occupancies, [1])
        assert_allclose(p.fractional_coordinates, frac)

    def test_fractional_boundary_tolerance_and_translated_operation(self):
        shift = np.eye(4)
        shift[0, 3] = .5
        p = PeriodicSystem(Molecule(["C", "O"], [[1 - 1e-9, 0, 0], [.25, .2, .3]], guess_bonds=False),
                           np.eye(3), symmetry_reduced=True, symmetry_operations=[shift])
        full = p.expand_symmetry()
        assert_allclose(full.fractional_coordinates, [[0, 0, 0], [.5, 0, 0], [.25, .2, .3], [.75, .2, .3]], atol=2e-9)
        self.assertEqual(records(full.reduce_symmetry().expand_symmetry()), records(full))

    def test_trigonal_centering_operations_and_disordered_sites(self):
        # A self-contained trigonal R-centered cell with special positions.
        rotations = []
        threefold = np.array([[0, -1, 0], [1, -1, 0], [0, 0, 1]])
        for center in ([0, 0, 0], [1/3, 2/3, 1/3], [2/3, 1/3, 2/3]):
            for sign in (1, -1):
                for power in range(3):
                    op = np.eye(4)
                    op[:3, :3] = sign * np.linalg.matrix_power(threefold, power)
                    op[:3, 3] = center
                    rotations.append(op)
        lattice = cell_matrix(3.97, 3.97, 28.59, 90, 90, 120)
        frac = np.array([[0, 0, .1], [0, 0, .2], [0, 0, .2]])
        p = PeriodicSystem(Molecule(["O", "Bi", "Sr"], frac @ lattice, guess_bonds=False), lattice,
                           symmetry_reduced=True, symmetry_operations=rotations, occupancies=[1, .74, .26])
        realized = p.expand_symmetry()
        self.assertEqual(len(realized.molecule.atoms), 18)
        self.assertEqual(len(realized.reduce_symmetry().molecule.atoms), 3)
        self.assertEqual(records(realized.reduce_symmetry().expand_symmetry()), records(realized))

    def test_cif_modern_symmetry_uncertainties_occupancies(self):
        p = self.read()
        self.assertEqual(p.molecule.atoms, ("C", "C", "O"))
        assert_allclose(p.fractional_coordinates, [[.1, .2, .3], [.9, .8, .7], [0, 0, 0]], atol=1e-14)
        assert_allclose(p.occupancies, [1, 1, .5])
        self.assertEqual(p.site_labels, ["C1", "C1", "O1"])
        self.assertEqual(p.space_group["space_group_it_number"], "2")
        self.assertIn("cif", p.metadata)
        assert_allclose(p.lattice, np.diag([4, 5, 6]) * UnitsData.convert("Angstroms", "BohrRadius"), atol=1e-14)

    def test_cif_reduced_legacy_and_dotted_tags(self):
        legacy = CIF.replace("space_group_symop_operation_xyz", "symmetry_equiv_pos_as_xyz")
        dotted = CIF.replace("space_group_symop_operation_xyz", "space_group_symop.operation_xyz")
        for text in (CIF, legacy, dotted):
            with self.subTest(text=text):
                p = self.read(text, symmetry_reduced=True)
                self.assertTrue(p.symmetry_reduced)
                self.assertEqual(p.molecule.atoms, ("C", "O"))
                self.assertEqual(records(p.expand_symmetry()), records(self.read()))

    def test_cif_readonly_path_and_multiple_blocks(self):
        text = CIF + CIF.replace("data_inversion", "data_second").replace("4.0(2)", "8.0")
        p = self.read(text, block=1)
        assert_allclose(p.lattice[0, 0], 8 * UnitsData.convert("Angstroms", "BohrRadius"))
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "sample.cif")
            with open(path, "w") as stream:
                stream.write(CIF)
            os.chmod(path, 0o444)
            try:
                p = PeriodicSystem.from_cif(path, guess_bonds=False)
                self.assertEqual(p.molecule.source_file, path)
                self.assertEqual(records(p), records(self.read()))
            finally:
                os.chmod(path, 0o644)

    def test_cif_loop_row_can_span_lines_and_symbols_can_use_labels(self):
        text = CIF.replace("C1 C 0.1(1) 0.2 0.3 1", "C1 C\n0.1(1) 0.2\n0.3 1")
        self.assertEqual(records(self.read(text)), records(self.read()))
        text = CIF.replace("_atom_site_type_symbol\n", "").replace("C1 C ", "C1 ").replace("O1 O ", "O1 ")
        self.assertEqual(records(self.read(text)), records(self.read()))

    def test_cif_triclinic_geometry(self):
        p = self.read(CIF.replace("_cell_angle_alpha 90", "_cell_angle_alpha 75")
                      .replace("_cell_angle_beta 90", "_cell_angle_beta 80")
                      .replace("_cell_angle_gamma 90", "_cell_angle_gamma 65"))
        expected = cell_matrix(4, 5, 6, 75, 80, 65) * UnitsData.convert("Angstroms", "BohrRadius")
        assert_allclose(p.lattice, expected)
        assert_allclose(np.asarray(p.molecule.coords), p.fractional_coordinates @ expected)

    def test_expansion_deduplicates_sites_and_special_positions(self):
        p = PeriodicSystem(Molecule(["C", "C", "O"], [[.1, .2, .3], [.9, .8, .7], [0, 0, 0]], guess_bonds=False),
                           np.eye(3), symmetry_reduced=True, symmetry_operations=inversion())
        p = p.expand_symmetry()
        self.assertEqual(p.molecule.atoms, ("C", "C", "O"))
        self.assertFalse(p.symmetry_reduced)

    def test_reduction_roundtrip_and_inconsistent_symmetry(self):
        p = self.read()
        reduced = p.reduce_symmetry()
        self.assertTrue(reduced.symmetry_reduced)
        self.assertEqual(reduced.molecule.atoms, ("C", "O"))
        self.assertEqual(records(reduced.expand_symmetry()), records(p))
        self.assertEqual(len(p.molecule.atoms), 3)
        bad = PeriodicSystem(Molecule(["C"], [[.1, .2, .3]], guess_bonds=False), np.eye(3),
                             symmetry_operations=inversion())
        with self.assertRaisesRegex(ValueError, "do not preserve"):
            bad.reduce_symmetry()

    def test_reduction_distinguishes_elements_isotopes_and_occupancies(self):
        molecule = Molecule(["C", "C", "C", "C", "O"],
                            [[.1, 0, 0], [.9, 0, 0], [.1, 0, 0], [.9, 0, 0], [0, 0, 0]],
                            masses=[12, 12, 13, 13, 16], guess_bonds=False)
        p = PeriodicSystem(molecule, np.eye(3), symmetry_operations=inversion(), occupancies=[1, 1, .5, .5, 1])
        r = p.reduce_symmetry()
        self.assertEqual(r.molecule.atoms, ("C", "C", "O"))
        assert_allclose(r.molecule.masses, [12, 13, 16])
        assert_allclose(r.occupancies, [1, .5, 1])
        self.assertEqual(len(r.expand_symmetry().molecule.atoms), 5)

    def test_tile_skew_cell_translation_order_and_origin(self):
        lattice = np.array([[3., 0, 0], [1., 4., 0], [.5, 1., 5.]])
        origin = np.array([2., 3., 4.])
        frac = np.array([[.1, .2, .3]])
        p = PeriodicSystem(Molecule(["C"], origin + frac @ lattice, guess_bonds=False), lattice, origin=origin)
        s = p.tile((2, 3, 1))
        assert_allclose(s.lattice, lattice * np.array([2, 3, 1])[:, None])
        assert_allclose(s.molecule.coords, [origin + (frac[0] + [i, j, 0]) @ lattice for i in range(2) for j in range(3)])
        self.assertEqual(len(s.molecule.atoms), 6)
        self.assertEqual(len(s.reduce_symmetry().molecule.atoms), 1)
        self.assertEqual(len(p.molecule.atoms), 1)
        assert_allclose(p.lattice, lattice)

    def test_tile_reduced_cell_and_transformed_symmetry_roundtrip(self):
        p = self.read(symmetry_reduced=True)
        s = p.tile((2, 1, 1))
        self.assertEqual(len(s.molecule.atoms), 6)
        self.assertFalse(s.symmetry_reduced)
        self.assertEqual(len(s.symmetry_operations), 4)
        self.assertEqual(records(s.reduce_symmetry().expand_symmetry()), records(s))
        self.assertEqual(len(p.molecule.atoms), 2)
        self.assertEqual(len(p.tile(2).molecule.atoms), 24)

    def test_anisotropic_tiling_drops_incompatible_rotations(self):
        swap = np.eye(4)
        swap[:3, :3] = [[0, 1, 0], [1, 0, 0], [0, 0, 1]]
        p = PeriodicSystem(Molecule(["C"], [[.1, .2, .3]], guess_bonds=False), np.eye(3),
                           symmetry_reduced=True, symmetry_operations=[swap])
        cubic = p.tile((2, 2, 1))
        stretched = p.tile((2, 1, 1))
        self.assertEqual(len(cubic.symmetry_operations), 8)
        self.assertEqual(len(stretched.symmetry_operations), 2)
        self.assertEqual(records(cubic.reduce_symmetry().expand_symmetry()), records(cubic))
        self.assertEqual(records(stretched.reduce_symmetry().expand_symmetry()), records(stretched))

    def test_explicit_bonds_and_masses_survive_expansion_and_tiling(self):
        molecule = Molecule(["C", "O"], [[.1, .2, .3], [.2, .3, .4]], bonds=[[0, 1, 2]], masses=[13, 18],
                            guess_bonds=False)
        p = PeriodicSystem(molecule, np.eye(3), symmetry_reduced=True, symmetry_operations=inversion())
        full = p.expand_symmetry()
        self.assertEqual(full.molecule.bonds, [[0, 2, 2], [1, 3, 2]])
        assert_allclose(full.molecule.masses, [13, 13, 18, 18])
        tiled = full.tile((2, 1, 1))
        self.assertEqual(tiled.molecule.bonds, [[0, 2, 2], [1, 3, 2], [4, 6, 2], [5, 7, 2]])

    def test_native_plot_cell_edges_units_and_user_annotations(self):
        p = self.read(symmetry_reduced=True)
        p.origin = np.array([2., 3., 4.])
        extra = Line([[0, 0, 0], [1, 0, 0]])
        for units in (None, "Angstroms"):
            fig, atoms, bonds, annotations = p.plot(bonds=False, rings=False, objects=True, return_objects=True,
                                                   units=units, cell_style={"color": "red"},
                                                   annotation_function=lambda *args: [extra])
            self.assertEqual(len(atoms[0]), 2)  # reduced molecule is plotted as stored
            self.assertEqual(len(annotations[0]), 13)
            self.assertIs(annotations[0][-1], extra)
            vertices = np.concatenate([line.points for line in annotations[0][:-1]])
            factor = 1 if units is None else UnitsData.convert("BohrRadius", units)
            assert_allclose(vertices.min(axis=0), p.cell_vertices.min(axis=0) * factor)
            assert_allclose(vertices.max(axis=0), p.cell_vertices.max(axis=0) * factor)
            self.assertTrue(all(line.opts["color"] == "red" for line in annotations[0][:-1]))

    def test_native_x3d_plot_and_disable_cell(self):
        p = self.read()
        fig, atoms, bonds, lines = p.plot(bonds=False, rings=False, return_objects=True)
        self.assertEqual(len(atoms[0]), 3)
        self.assertEqual(len(lines[0]), 12)
        self.assertIn("X3D", type(fig.figure).__name__)
        fig, atoms, bonds, lines = p.plot(bonds=False, rings=False, draw_cell=False, return_objects=True)
        self.assertFalse(lines[0])

    def test_cif_missing_symmetry_requires_non_p1_operations(self):
        text = CIF[:CIF.index("loop_")] + CIF[CIF.index("loop_\n_atom_site_label"):]
        with self.assertRaisesRegex(ValueError, "no explicit symmetry"):
            self.read(text)
        self.assertEqual(len(self.read(text, symmetry_operations=inversion()).molecule.atoms), 3)
        p1 = text.replace("'P -1'", "'P 1'").replace("_space_group_IT_number 2", "_space_group_IT_number 1")
        self.assertEqual(len(self.read(p1).molecule.atoms), 2)

    def test_validation(self):
        molecule = Molecule(["C"], [[0., 0., 0.]], guess_bonds=False)
        for lattice in (np.zeros((3, 3)), np.ones((2, 3)), np.eye(3) * np.nan):
            with self.subTest(lattice=lattice), self.assertRaises(ValueError):
                PeriodicSystem(molecule, lattice)
        for options in ({"occupancies": [.5, 1]}, {"occupancies": [np.nan]}, {"occupancies": [2]},
                        {"site_labels": []}, {"origin": [1, 2]}, {"symmetry_operations": np.eye(3)}):
            with self.subTest(options=options), self.assertRaises(ValueError):
                PeriodicSystem(molecule, np.eye(3), **options)
        p = PeriodicSystem(molecule, np.eye(3))
        for repeats in ((0, 1, 1), (1, 2), (-1, 1, 1), (1.5, 1, 1), True):
            with self.subTest(repeats=repeats), self.assertRaises(ValueError):
                p.tile(repeats)
        for tol in (0, -1, np.nan, .5):
            with self.subTest(tol=tol), self.assertRaises(ValueError):
                p.reduce_symmetry(tol=tol)
        with self.assertRaisesRegex(ValueError, "block"):
            self.read(block=1)
        with self.assertRaisesRegex(ValueError, "missing required"):
            self.read(CIF.replace("_cell_length_a 4.0(2)\n", ""))
        with self.assertRaisesRegex(ValueError, "cell angles"):
            self.read(CIF.replace("_cell_angle_alpha 90", "_cell_angle_alpha 170")
                      .replace("_cell_angle_beta 90", "_cell_angle_beta 10")
                      .replace("_cell_angle_gamma 90", "_cell_angle_gamma 10"))


if __name__ == "__main__":
    unittest.main()
