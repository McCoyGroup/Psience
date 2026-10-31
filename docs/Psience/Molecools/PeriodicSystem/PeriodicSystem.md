## <a id="Psience.Molecools.PeriodicSystem.PeriodicSystem">PeriodicSystem</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/PeriodicSystem.py#L81)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/PeriodicSystem.py#L81?message=Update%20Docs)]
</div>

A molecule plus a periodic lattice, without inheriting from Molecule.

``lattice`` contains the three Cartesian lattice vectors as **rows**. Its
lengths and ``origin`` use ``lattice_units`` (Bohr by default), while the
stored molecule follows Molecule's usual Bohr convention. Internally all
lengths are Bohr. Fractional coordinates obey ``xyz = origin + frac @ lattice``.

Set ``symmetry_reduced=True`` when the molecule contains asymmetric-unit
representatives instead of a realized cell. ``symmetry_operations`` is the
complete list of fractional-coordinate 4x4 affine matrices, acting on column
vectors. Identity is included automatically. Tolerances are fractional.
``periodicity`` enables translation/wrapping along a, b, and c independently;
retain independent lattice rows for finite directions (e.g. slab vacuum).







<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 
<a id="Psience.Molecools.PeriodicSystem.PeriodicSystem.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, molecule, lattice, *, symmetry_reduced=False, symmetry_operations=None, origin=None, lattice_units='BohrRadius', space_group=None, site_labels=None, occupancies=None, metadata=None, periodicity=(True, True, True)): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/PeriodicSystem.py#L97)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/PeriodicSystem.py#L97?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.PeriodicSystem.PeriodicSystem.fractional_coordinates" class="docs-object-method">&nbsp;</a> 
```python
@property
fractional_coordinates(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/PeriodicSystem/PeriodicSystem.py#L155)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/PeriodicSystem/PeriodicSystem.py#L155?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.PeriodicSystem.PeriodicSystem.cell_vertices" class="docs-object-method">&nbsp;</a> 
```python
@property
cell_vertices(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/PeriodicSystem/PeriodicSystem.py#L159)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/PeriodicSystem/PeriodicSystem.py#L159?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.PeriodicSystem.PeriodicSystem.copy" class="docs-object-method">&nbsp;</a> 
```python
copy(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/PeriodicSystem/PeriodicSystem.py#L218)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/PeriodicSystem/PeriodicSystem.py#L218?message=Update%20Docs)]
</div>


<a id="Psience.Molecools.PeriodicSystem.PeriodicSystem.expand_symmetry" class="docs-object-method">&nbsp;</a> 
```python
expand_symmetry(self, tol=1e-06): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/PeriodicSystem/PeriodicSystem.py#L227)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/PeriodicSystem/PeriodicSystem.py#L227?message=Update%20Docs)]
</div>
Realize a reduced cell, wrapping sites and removing equivalent duplicates.


<a id="Psience.Molecools.PeriodicSystem.PeriodicSystem.reduce_symmetry" class="docs-object-method">&nbsp;</a> 
```python
reduce_symmetry(self, tol=1e-06): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/PeriodicSystem/PeriodicSystem.py#L252)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/PeriodicSystem/PeriodicSystem.py#L252?message=Update%20Docs)]
</div>
Keep one representative per supplied symmetry orbit, preserving species/occupancy.

The molecule must contain a complete realized cell for these operations.
This reduces the atom list, not the lattice to a primitive cell.


<a id="Psience.Molecools.PeriodicSystem.PeriodicSystem.tile" class="docs-object-method">&nbsp;</a> 
```python
tile(self, repetitions, tol=1e-06): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/PeriodicSystem/PeriodicSystem.py#L278)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/PeriodicSystem/PeriodicSystem.py#L278?message=Update%20Docs)]
</div>
Return a realized n_a x n_b x n_c supercell (or n x n x n for an integer).

Retain operations compatible with the supercell lattice and add the
internal cell translations. Supplied bonds are replicated, and Molecule's
normal bond inference connects copies when guess_bonds is enabled.


<a id="Psience.Molecools.PeriodicSystem.PeriodicSystem.from_cif" class="docs-object-method">&nbsp;</a> 
```python
@classmethod
from_cif(cls, file, *, block=0, symmetry_reduced=False, symmetry_operations=None, periodicity=(True, True, True), **molecule_options): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/classmethod.py#L317)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/classmethod.py#L317?message=Update%20Docs)]
</div>
Read one CIF block using McUtils (cell lengths in Angstroms, angles in degrees).

Preserve labels, occupancies, space-group fields, and the source CIF block.
Non-P1 files without explicit operations require ``symmetry_operations``;
space-group names alone are not used to guess a setting or origin.


<a id="Psience.Molecools.PeriodicSystem.PeriodicSystem.plot" class="docs-object-method">&nbsp;</a> 
```python
plot(self, *, draw_cell=True, cell_style=None, unit_axes=False, unit_axes_style=None, unit_axes_label_style=None, unit_axes_origin=None, unit_axes_scale=1, show_periodic_images=None, periodic_images=None, periodic_image_styles=None, backend='x3d', units='Angstroms', **options): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Molecools/PeriodicSystem/PeriodicSystem.py#L563)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/PeriodicSystem/PeriodicSystem.py#L563?message=Update%20Docs)]
</div>
Plot the stored molecule with its periodic boundary and optional neighbors.

``cell_style`` and ``unit_axes_style`` accept common styles, a/b/c overrides,
or three style dictionaries. A cell ``radius`` draws cylindrical edges;
otherwise edges are lines. Radii use the display units. ``unit_axes=True``
adds labeled lattice-vector arrows along periodic directions only.

Images form a full neighborhood outside the central cell: counts [2,0,1]
mean +/-2 cells along a, none along b, +/-1 along c, including corners.
An entry [negative,positive] sets the two sides independently: [2,0,[1,0]]
adds one layer in -c and none in +c. Booleans mean counts 1/0.
Finite directions never generate images. Supplying
counts enables images unless ``show_periodic_images=False``. Shared image
material styles default to transparency=.5; nested atom_style/bond_style
entries override those defaults. Only the central cell gets axes/annotations.
With bond inference enabled, copies are connected using the same Molecule
rules as tile(); native split bonds inherit each endpoint's cell styling.
With return_objects=True, image objects follow the central geometry's objects.
 </div>
</div>




## Examples













<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#Tests-68a1b3" markdown="1"> Tests</a> <a class="float-right" data-toggle="collapse" href="#Tests-68a1b3"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Tests-68a1b3" markdown="1">
 - [wrapper_units_origin_and_independent_copy](#wrapper_units_origin_and_independent_copy)
- [fractional_boundary_tolerance_and_translated_operation](#fractional_boundary_tolerance_and_translated_operation)
- [trigonal_centering_operations_and_disordered_sites](#trigonal_centering_operations_and_disordered_sites)
- [cif_modern_symmetry_uncertainties_occupancies](#cif_modern_symmetry_uncertainties_occupancies)
- [cif_reduced_legacy_and_dotted_tags](#cif_reduced_legacy_and_dotted_tags)
- [cif_readonly_path_and_multiple_blocks](#cif_readonly_path_and_multiple_blocks)
- [cif_loop_row_can_span_lines_and_symbols_can_use_labels](#cif_loop_row_can_span_lines_and_symbols_can_use_labels)
- [cif_triclinic_geometry](#cif_triclinic_geometry)
- [expansion_deduplicates_sites_and_special_positions](#expansion_deduplicates_sites_and_special_positions)
- [reduction_roundtrip_and_inconsistent_symmetry](#reduction_roundtrip_and_inconsistent_symmetry)
- [reduction_distinguishes_elements_isotopes_and_occupancies](#reduction_distinguishes_elements_isotopes_and_occupancies)
- [tile_skew_cell_translation_order_and_origin](#tile_skew_cell_translation_order_and_origin)
- [tile_reduced_cell_and_transformed_symmetry_roundtrip](#tile_reduced_cell_and_transformed_symmetry_roundtrip)
- [anisotropic_tiling_drops_incompatible_rotations](#anisotropic_tiling_drops_incompatible_rotations)
- [explicit_bonds_and_masses_survive_expansion_and_tiling](#explicit_bonds_and_masses_survive_expansion_and_tiling)
- [native_plot_cell_edges_units_and_user_annotations](#native_plot_cell_edges_units_and_user_annotations)
- [native_x3d_plot_and_disable_cell](#native_x3d_plot_and_disable_cell)
- [cif_missing_symmetry_requires_non_p1_operations](#cif_missing_symmetry_requires_non_p1_operations)
- [validation](#validation)

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
### <a class="collapse-link" data-toggle="collapse" href="#Setup-8643b1" markdown="1"> Setup</a> <a class="float-right" data-toggle="collapse" href="#Setup-8643b1"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Setup-8643b1" markdown="1">
 
Before we can run our examples we should get a bit of setup out of the way.
Since these examples were harvested from the unit tests not all pieces
will be necessary for all situations.

All tests are wrapped in a test class
```python
class PeriodicSystemTests(unittest.TestCase):
    def read(self, text=CIF, **options):
        return PeriodicSystem.from_cif(io.StringIO(text), guess_bonds=False, **options)
```

 </div>
</div>

#### <a name="wrapper_units_origin_and_independent_copy">wrapper_units_origin_and_independent_copy</a>
```python
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
```

#### <a name="fractional_boundary_tolerance_and_translated_operation">fractional_boundary_tolerance_and_translated_operation</a>
```python
    def test_fractional_boundary_tolerance_and_translated_operation(self):
        shift = np.eye(4)
        shift[0, 3] = .5
        p = PeriodicSystem(Molecule(["C", "O"], [[1 - 1e-9, 0, 0], [.25, .2, .3]], guess_bonds=False),
                           np.eye(3), symmetry_reduced=True, symmetry_operations=[shift])
        full = p.expand_symmetry()
        assert_allclose(full.fractional_coordinates, [[0, 0, 0], [.5, 0, 0], [.25, .2, .3], [.75, .2, .3]], atol=2e-9)
        self.assertEqual(records(full.reduce_symmetry().expand_symmetry()), records(full))
```

#### <a name="trigonal_centering_operations_and_disordered_sites">trigonal_centering_operations_and_disordered_sites</a>
```python
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
```

#### <a name="cif_modern_symmetry_uncertainties_occupancies">cif_modern_symmetry_uncertainties_occupancies</a>
```python
    def test_cif_modern_symmetry_uncertainties_occupancies(self):
        p = self.read()
        self.assertEqual(p.molecule.atoms, ("C", "C", "O"))
        assert_allclose(p.fractional_coordinates, [[.1, .2, .3], [.9, .8, .7], [0, 0, 0]], atol=1e-14)
        assert_allclose(p.occupancies, [1, 1, .5])
        self.assertEqual(p.site_labels, ["C1", "C1", "O1"])
        self.assertEqual(p.space_group["space_group_it_number"], "2")
        self.assertIn("cif", p.metadata)
        assert_allclose(p.lattice, np.diag([4, 5, 6]) * UnitsData.convert("Angstroms", "BohrRadius"), atol=1e-14)
```

#### <a name="cif_reduced_legacy_and_dotted_tags">cif_reduced_legacy_and_dotted_tags</a>
```python
    def test_cif_reduced_legacy_and_dotted_tags(self):
        legacy = CIF.replace("space_group_symop_operation_xyz", "symmetry_equiv_pos_as_xyz")
        dotted = CIF.replace("space_group_symop_operation_xyz", "space_group_symop.operation_xyz")
        for text in (CIF, legacy, dotted):
            with self.subTest(text=text):
                p = self.read(text, symmetry_reduced=True)
                self.assertTrue(p.symmetry_reduced)
                self.assertEqual(p.molecule.atoms, ("C", "O"))
                self.assertEqual(records(p.expand_symmetry()), records(self.read()))
```

#### <a name="cif_readonly_path_and_multiple_blocks">cif_readonly_path_and_multiple_blocks</a>
```python
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
```

#### <a name="cif_loop_row_can_span_lines_and_symbols_can_use_labels">cif_loop_row_can_span_lines_and_symbols_can_use_labels</a>
```python
    def test_cif_loop_row_can_span_lines_and_symbols_can_use_labels(self):
        text = CIF.replace("C1 C 0.1(1) 0.2 0.3 1", "C1 C\n0.1(1) 0.2\n0.3 1")
        self.assertEqual(records(self.read(text)), records(self.read()))
        text = CIF.replace("_atom_site_type_symbol\n", "").replace("C1 C ", "C1 ").replace("O1 O ", "O1 ")
        self.assertEqual(records(self.read(text)), records(self.read()))
```

#### <a name="cif_triclinic_geometry">cif_triclinic_geometry</a>
```python
    def test_cif_triclinic_geometry(self):
        p = self.read(CIF.replace("_cell_angle_alpha 90", "_cell_angle_alpha 75")
                      .replace("_cell_angle_beta 90", "_cell_angle_beta 80")
                      .replace("_cell_angle_gamma 90", "_cell_angle_gamma 65"))
        expected = cell_matrix(4, 5, 6, 75, 80, 65) * UnitsData.convert("Angstroms", "BohrRadius")
        assert_allclose(p.lattice, expected)
        assert_allclose(np.asarray(p.molecule.coords), p.fractional_coordinates @ expected)
```

#### <a name="expansion_deduplicates_sites_and_special_positions">expansion_deduplicates_sites_and_special_positions</a>
```python
    def test_expansion_deduplicates_sites_and_special_positions(self):
        p = PeriodicSystem(Molecule(["C", "C", "O"], [[.1, .2, .3], [.9, .8, .7], [0, 0, 0]], guess_bonds=False),
                           np.eye(3), symmetry_reduced=True, symmetry_operations=inversion())
        p = p.expand_symmetry()
        self.assertEqual(p.molecule.atoms, ("C", "C", "O"))
        self.assertFalse(p.symmetry_reduced)
```

#### <a name="reduction_roundtrip_and_inconsistent_symmetry">reduction_roundtrip_and_inconsistent_symmetry</a>
```python
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
```

#### <a name="reduction_distinguishes_elements_isotopes_and_occupancies">reduction_distinguishes_elements_isotopes_and_occupancies</a>
```python
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
```

#### <a name="tile_skew_cell_translation_order_and_origin">tile_skew_cell_translation_order_and_origin</a>
```python
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
```

#### <a name="tile_reduced_cell_and_transformed_symmetry_roundtrip">tile_reduced_cell_and_transformed_symmetry_roundtrip</a>
```python
    def test_tile_reduced_cell_and_transformed_symmetry_roundtrip(self):
        p = self.read(symmetry_reduced=True)
        s = p.tile((2, 1, 1))
        self.assertEqual(len(s.molecule.atoms), 6)
        self.assertFalse(s.symmetry_reduced)
        self.assertEqual(len(s.symmetry_operations), 4)
        self.assertEqual(records(s.reduce_symmetry().expand_symmetry()), records(s))
        self.assertEqual(len(p.molecule.atoms), 2)
        self.assertEqual(len(p.tile(2).molecule.atoms), 24)
```

#### <a name="anisotropic_tiling_drops_incompatible_rotations">anisotropic_tiling_drops_incompatible_rotations</a>
```python
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
```

#### <a name="explicit_bonds_and_masses_survive_expansion_and_tiling">explicit_bonds_and_masses_survive_expansion_and_tiling</a>
```python
    def test_explicit_bonds_and_masses_survive_expansion_and_tiling(self):
        molecule = Molecule(["C", "O"], [[.1, .2, .3], [.2, .3, .4]], bonds=[[0, 1, 2]], masses=[13, 18],
                            guess_bonds=False)
        p = PeriodicSystem(molecule, np.eye(3), symmetry_reduced=True, symmetry_operations=inversion())
        full = p.expand_symmetry()
        self.assertEqual(full.molecule.bonds, [[0, 2, 2], [1, 3, 2]])
        assert_allclose(full.molecule.masses, [13, 13, 18, 18])
        tiled = full.tile((2, 1, 1))
        self.assertEqual(tiled.molecule.bonds, [[0, 2, 2], [1, 3, 2], [4, 6, 2], [5, 7, 2]])
```

#### <a name="native_plot_cell_edges_units_and_user_annotations">native_plot_cell_edges_units_and_user_annotations</a>
```python
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
```

#### <a name="native_x3d_plot_and_disable_cell">native_x3d_plot_and_disable_cell</a>
```python
    def test_native_x3d_plot_and_disable_cell(self):
        p = self.read()
        fig, atoms, bonds, lines = p.plot(bonds=False, rings=False, return_objects=True)
        self.assertEqual(len(atoms[0]), 3)
        self.assertEqual(len(lines[0]), 12)
        self.assertIn("X3D", type(fig.figure).__name__)
        fig, atoms, bonds, lines = p.plot(bonds=False, rings=False, draw_cell=False, return_objects=True)
        self.assertFalse(lines[0])
```

#### <a name="cif_missing_symmetry_requires_non_p1_operations">cif_missing_symmetry_requires_non_p1_operations</a>
```python
    def test_cif_missing_symmetry_requires_non_p1_operations(self):
        text = CIF[:CIF.index("loop_")] + CIF[CIF.index("loop_\n_atom_site_label"):]
        with self.assertRaisesRegex(ValueError, "no explicit symmetry"):
            self.read(text)
        self.assertEqual(len(self.read(text, symmetry_operations=inversion()).molecule.atoms), 3)
        p1 = text.replace("'P -1'", "'P 1'").replace("_space_group_IT_number 2", "_space_group_IT_number 1")
        self.assertEqual(len(self.read(p1).molecule.atoms), 2)
```

#### <a name="validation">validation</a>
```python
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
```

 </div>
</div>






---


<div markdown="1" class="text-secondary">
<div class="container">
  <div class="row">
   <div class="col" markdown="1">
**Feedback**   
</div>
   <div class="col" markdown="1">
**Examples**   
</div>
   <div class="col" markdown="1">
**Templates**   
</div>
   <div class="col" markdown="1">
**Documentation**   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
</div>
  <div class="row">
   <div class="col" markdown="1">
[Bug](https://github.com/McCoyGroup/Psience/issues/new?title=Documentation%20Improvement%20Needed)/[Request](https://github.com/McCoyGroup/Psience/issues/new?title=Example%20Request)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/ci/docs/examples/Psience/Molecools/PeriodicSystem/PeriodicSystem.md)/[New](https://github.com/McCoyGroup/Psience/new/master/?filename=ci/docs/examples/Psience/Molecools/PeriodicSystem/PeriodicSystem.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/ci/docs/templates/Psience/Molecools/PeriodicSystem/PeriodicSystem.md)/[New](https://github.com/McCoyGroup/Psience/new/master/?filename=ci/docs/templates/Psience/Molecools/PeriodicSystem/PeriodicSystem.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Molecools/PeriodicSystem.py#L81?message=Update%20Docs)   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
</div>
</div>
</div>