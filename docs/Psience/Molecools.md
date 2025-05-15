# <a id="Psience.Molecools">Psience.Molecools</a> 
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Molecools/__init__.py#L1)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/__init__.py#L1?message=Update%20Docs)]
</div>
    
Molecules provides wrapper utilities for working with and visualizing molecular systems

### Members
<div class="container alert alert-secondary bg-light">
  <div class="row">
   <div class="col" markdown="1">
[MolecularVibrations](Molecools/Vibrations/MolecularVibrations.md)   
</div>
   <div class="col" markdown="1">
[MolecularNormalModes](Molecools/Vibrations/MolecularNormalModes.md)   
</div>
   <div class="col" markdown="1">
[Molecule](Molecools/Molecule/Molecule.md)   
</div>
</div>
  <div class="row">
   <div class="col" markdown="1">
[MolecoolException](Molecools/Molecule/MolecoolException.md)   
</div>
   <div class="col" markdown="1">
[StructuralProperties](Molecools/Properties/StructuralProperties.md)   
</div>
   <div class="col" markdown="1">
[BondingProperties](Molecools/Properties/BondingProperties.md)   
</div>
</div>
  <div class="row">
   <div class="col" markdown="1">
[MolecularProperties](Molecools/Properties/MolecularProperties.md)   
</div>
   <div class="col" markdown="1">
[MolecularPropertyError](Molecools/Properties/MolecularPropertyError.md)   
</div>
   <div class="col" markdown="1">
[OpenBabelMolManager](Molecools/Properties/OpenBabelMolManager.md)   
</div>
</div>
  <div class="row">
   <div class="col" markdown="1">
[DipoleSurfaceManager](Molecools/Properties/DipoleSurfaceManager.md)   
</div>
   <div class="col" markdown="1">
[PotentialSurfaceManager](Molecools/Properties/PotentialSurfaceManager.md)   
</div>
   <div class="col" markdown="1">
[NormalModesManager](Molecools/Properties/NormalModesManager.md)   
</div>
</div>
  <div class="row">
   <div class="col" markdown="1">
[MolecularEmbedding](Molecools/CoordinateSystems/MolecularEmbedding.md)   
</div>
   <div class="col" markdown="1">
[ModeEmbedding](Molecools/CoordinateSystems/ModeEmbedding.md)   
</div>
   <div class="col" markdown="1">
[MolecularZMatrixCoordinateSystem](Molecools/CoordinateSystems/MolecularZMatrixCoordinateSystem.md)   
</div>
</div>
  <div class="row">
   <div class="col" markdown="1">
[MolecularCartesianCoordinateSystem](Molecools/CoordinateSystems/MolecularCartesianCoordinateSystem.md)   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
</div>
</div>





## Examples













<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#Tests-e56515" markdown="1"> Tests</a> <a class="float-right" data-toggle="collapse" href="#Tests-e56515"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Tests-e56515" markdown="1">
 - [NormalModeRephasing](#NormalModeRephasing)
- [MolecularGMatrix](#MolecularGMatrix)
- [ImportMolecule](#ImportMolecule)
- [PrincipleAxisEmbedding](#PrincipleAxisEmbedding)
- [EckartEmbed](#EckartEmbed)
- [Eckart](#Eckart)
- [HOONODihedral](#HOONODihedral)
- [EckartEmbedDipoles](#EckartEmbedDipoles)
- [EckartEmbedMolecule](#EckartEmbedMolecule)
- [EmbeddedMolecule](#EmbeddedMolecule)
- [AddDummyAtoms](#AddDummyAtoms)
- [AddDummyAtomProperties](#AddDummyAtomProperties)
- [AddDummyAtomJacobians](#AddDummyAtomJacobians)
- [InternalCoordOrder](#InternalCoordOrder)
- [Plotting](#Plotting)
- [BondGuessing](#BondGuessing)
- [Frags](#Frags)
- [AutoZMat](#AutoZMat)
- [HODModes](#HODModes)
- [H2OModes](#H2OModes)
- [RenormalizeGaussianModes](#RenormalizeGaussianModes)
- [VisualizeNormalModes](#VisualizeNormalModes)
- [InternalCartesianJacobians](#InternalCartesianJacobians)
- [CompositeCoordinates](#CompositeCoordinates)
- [RDKitSpectrum](#RDKitSpectrum)
- [ExpansionPotential](#ExpansionPotential)
- [OpenBabel](#OpenBabel)
- [1DPotentialReps](#1DPotentialReps)
- [Constructors](#Constructors)
- [InternalConv](#InternalConv)
- [AutomaticConversion](#AutomaticConversion)
- [FastInternals](#FastInternals)

<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
### <a class="collapse-link" data-toggle="collapse" href="#Setup-bad481" markdown="1"> Setup</a> <a class="float-right" data-toggle="collapse" href="#Setup-bad481"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="Setup-bad481" markdown="1">
 
Before we can run our examples we should get a bit of setup out of the way.
Since these examples were harvested from the unit tests not all pieces
will be necessary for all situations.

All tests are wrapped in a test class
```python
class MolecoolsTests(TestCase):
    def setUp(self):
        self.test_log_water = TestManager.test_data("water_OH_scan.log")
        self.test_log_freq = TestManager.test_data("water_freq.log")
        self.test_HOD = TestManager.test_data("HOD_freq.fchk")
        self.test_fchk = TestManager.test_data("water_freq.fchk")
        self.test_log_h2 = TestManager.test_data("outer_H2_scan_new.log")
    def setup_OCHH(cls, optimize=True):
        from McUtils.Extensions import ModuleLoader

        loader = ModuleLoader(os.path.expanduser("~/Documents/Postdoc/Projects/DGB"))
        h2co_mod = loader.load("H2COPot")

        def internal_pot(coords, order=None):
            coords = coords[..., (0, 1, 3, 2, 4, 5)]
            vals = h2co_mod.InternalsPotential.get_pot(coords)
            return vals

        ochh = Molecule.from_file(
            TestManager.test_data('OCHH_freq.fchk'),
            energy_evaluator={
                'potential_function': internal_pot,
                "distance_units": "Angstroms",
                "energy_units": "Wavenumbers",
                "strip_embedding": True,
            },
            internals=[
                [0, -1, -1, -1],
                [1, 0, -1, -1],
                [2, 1, 0, -1],
                [3, 1, 2, 0]
            ]
        )
        if optimize:
            base_dip = ochh.dipole_derivatives
            ochh = ochh.optimize(
                # method='quasi-newton'
                method='conjugate-gradient'
                # method='gradient-descent'
                , max_iterations=50
                , stencil=3
                # , logger=True
                # , max_displacement=.01
                , prevent_oscillations=3
                , restart_interval=15
            ).modify(dipole_derivatives=base_dip)
            # ochh = Molecule(
            #     ochh.atoms,
            #     ochh_opt.coords,
            #     energy_evaluator={
            #         'potential_function': internal_pot,
            #         "distance_units": "Angstroms",
            #         "energy_units": "Wavenumbers",
            #         "strip_embedding": True,
            #     },
            #     internals=[
            #         [0, -1, -1, -1],
            #         [1,  0, -1, -1],
            #         [2,  1,  0, -1],
            #         [3,  1,  0,  2],
            #     ]
            # )

        return ochh
```

 </div>
</div>

#### <a name="NormalModeRephasing">NormalModeRephasing</a>
```python
    def test_NormalModeRephasing(self):
        m_16 = Molecule.from_file(TestManager.test_data('CH2DT_freq_16.fchk'))
        m_09 = Molecule.from_file(TestManager.test_data('CH2DT_freq.fchk'))
        modes_09 = m_09.normal_modes.modes
        # modes_16 = m_16.normal_modes

        modes_09 = np.array([x / np.linalg.norm(x) for x in modes_09.basis.matrix.T])
        # modes_16 = np.array([x / np.linalg.norm(x) for x in modes_16.basis.matrix.T])

        phases = m_16.normal_modes.get_fchk_normal_mode_rephasing()
        rescaled = m_16.normal_modes.modes.rescale(phases)

        rescaled_16 = np.array([x / np.linalg.norm(x) for x in rescaled.basis.matrix.T])

        phase_test = np.sign(np.diag(np.dot(modes_09, rescaled_16.T)))

        self.assertEquals(np.sum(np.diff(phase_test)), 0)
```

#### <a name="MolecularGMatrix">MolecularGMatrix</a>
```python
    def test_MolecularGMatrix(self):
        mol = Molecule.from_file(self.test_fchk)
        mol.zmatrix = [
            [0, -1, -1, -1],
            [1,  0, -1, -1],
            [2,  0,  1, -1]
        ]
        g = mol.g_matrix

        self.assertEquals(g.shape, (3, 3))
```

#### <a name="ImportMolecule">ImportMolecule</a>
```python
    def test_ImportMolecule(self):

        n = 3 # water
        m = Molecule.from_file(self.test_fchk)
        self.assertEquals(m.atoms, ("O", "H", "H"))
```

#### <a name="PrincipleAxisEmbedding">PrincipleAxisEmbedding</a>
```python
    def test_PrincipleAxisEmbedding(self):
        ref_file = TestManager.test_data("tbhp_180.fchk")

        ref = Molecule.from_file(ref_file)
        self.assertEquals(ref.center_of_mass.tolist(),
                          [-0.10886336323443993, -7.292720327263524e-05, -0.04764041570644441]
                          )

        ref_inerts = [
                [-0.9998646051394727,    -1.6944914059526497e-5, -1.6455123887595957e-2],
                [-1.6455007408638932e-2,  4.930578772682442e-3,   0.9998524501765987],
                [-6.419087070136426e-5,  -0.9999878444790397,     4.930190026343585e-3]
            ]
        inerts = ref.inertial_axes
        test_inerts = (inerts * np.array([-1, 1, 1])).T
        self.assertTrue(np.allclose(
            test_inerts,
            ref_inerts
        ),
            msg="principle axes {} and {} don't align".format(
                test_inerts,
                ref_inerts
            )
        )

        pax_rot = ref.principle_axis_frame()  # type: MolecularTransformation
        self.assertTrue(np.allclose(
            pax_rot.transformation_function.transform,
            inerts.T
        ))
        rot_ref = pax_rot.apply(ref)

        # g, _, _ = ref.plot(atom_style=dict(color='black'))
        # rot_ref.plot(figure=g)
        # g.show()

        self.assertTrue(np.allclose(
            rot_ref.center_of_mass,
            [0., 0., 0.]
        ),
            msg="COM: {} was {}".format(rot_ref.center_of_mass, ref.center_of_mass))

        test_coords = np.matmul(
                    (ref.coords - ref.center_of_mass[np.newaxis])[:, np.newaxis, :],
                    inerts
                ).squeeze()
        # raise Exception(rot_ref.coords, test_coords)
        self.assertTrue(
            np.allclose(
                rot_ref.coords,
                test_coords
            )
        )

        mathematica_coords = np.array([
                [ 2.094928525160645e-4,   8.85212868882308e-2,  -0.8400509406910139],
                [ 2.389396575506497,      1.697491740062459,    -0.8428256390972853],
                [ 2.435043833038253,      2.934952064361808,     0.7950074811481486],
                [ 4.0560845074996985,     0.4921123166233054,   -0.8003781737352631],
                [-4.983484850171475e-3,  -1.5885626031388058,    1.2992229461755922],
                [-1.7490151872158886e-4, -1.8815600632167903e-3, 3.5774728125123842],
                [-4.314406779968471e-3,  -1.3424852433777361,    4.810480604689872],
                [-4.312429484356625e-3,  -1.7659250558813848,   -3.0429810385290326],
                [-1.6805757842711242,    -2.9559004963767235,   -2.984461679814903],
                [ 1.663962078887355,      -2.9669237481136603,   -2.9820756778710344],
                [ 4.171884239172418e-4,   -0.7242576512048614,   -4.816727043081511],
                [-2.3797319162701913,     1.7110998385574014,   -0.8442221100234485],
                [-4.053502667206945,      0.5153958278660512,   -0.8051208327551433],
                [-2.439171179603177,      2.871593767591361,    -2.543401568931165],
                [-2.419963556488472,      2.947396453869957,     0.7945604672548087],
                [2.4576648430627377,      2.8566629998551765,   -2.5425989365331256]
            ])
        self.assertTrue(np.allclose(
            rot_ref.coords,
            -mathematica_coords[:, (2, 1, 0)]
        ),
        msg="{} but mathematica {}".format(
            rot_ref.coords,
            -mathematica_coords[:, (2, 1, 0)]
        ))
```

#### <a name="EckartEmbed">EckartEmbed</a>
```python
    def test_EckartEmbed(self):
        m = Molecule.from_file(TestManager.test_data('HOH_freq.fchk'))
        crd = m.embed_coords(m.coords)
        self.assertTrue(np.allclose(m.coords, crd))
```

#### <a name="Eckart">Eckart</a>
```python
    def test_Eckart(self):
        scan_file = TestManager.test_data("tbhp_030.log")
        ref_file = TestManager.test_data("tbhp_180.fchk")

        scan = Molecule.from_file(scan_file)
        ref = Molecule.from_file(ref_file)
        sel = np.where(ref.masses > 3)[0]

        pax_rot = ref.principle_axis_frame(sel=sel) #type: MolecularTransformation
        rot_ref = pax_rot.apply(ref)

        self.assertTrue(np.allclose(
            rot_ref.center_of_mass,
            [0., 0., 0.]
        ))

        #transf = scan.principle_axis_frame(sel=sel)
        transf = scan.eckart_frame(rot_ref, sel=sel)
        tf_test = transf[0].transformation_function

        tf_mat = tf_test.transform
        self.assertTrue(np.allclose(tf_mat@tf_mat.T - np.eye(3), 0.))
        self.assertEquals(tf_test.transf.shape, (4, 4))

        for t, m in zip(transf, scan):
            # t = m.principle_axis_frame(sel=sel)  # type: MolecularTransformation

            new_mol = t(m)
            # rot_ref.guess_bonds = False
            # ref.guess_bonds = False
            # m.guess_bonds = False
            # new_mol.guess_bonds = False
            # m = m #type: Molecule
            # # g1, a, b = ref.plot()
            # # ref.plot(figure=g1)
            # # rot_ref.plot(figure=g1)
            # g, a, b = new_mol.plot()
            # rot_ref.plot(figure=g, atom_style=dict(color='black'))
            # g.show()

            fuckup = np.linalg.norm(new_mol.coords[sel] - rot_ref.coords[sel])
            self.assertLess(fuckup / len(sel), .1,
                            msg="new: {}\nref: {}".format(
                                new_mol.coords,
                                rot_ref.coords
                            )
                            )

            # transf = scan.principle_axis_frame(sel=sel)
        transf = scan.eckart_frame(ref, sel=sel)
        for t, m in zip(transf, scan):
            # t = m.principle_axis_frame(sel=sel)  # type: MolecularTransformation

            new_mol = t(m)
            # rot_ref.guess_bonds = False
            # ref.guess_bonds = False
            # m.guess_bonds = False
            # new_mol.guess_bonds = False
            # m = m #type: Molecule
            # # g1, a, b = ref.plot()
            # # ref.plot(figure=g1)
            # # rot_ref.plot(figure=g1)
            # g, a, b = new_mol.plot()
            # rot_ref.plot(figure=g, atom_style=dict(color='black'))
            # g.show()

            fuckup = np.linalg.norm(new_mol.coords[sel] - ref.coords[sel])
            self.assertLess(fuckup / len(sel), .1,
                            msg="new: {}\nref: {}".format(
                                new_mol.coords,
                                ref.coords
                            )
                            )
```

#### <a name="HOONODihedral">HOONODihedral</a>
```python
    def test_HOONODihedral(self):
        # should be broken

        mol = Molecule.from_file(TestManager.test_data('HOONO_freq.fchk'))
        mol.zmatrix = [
            [1, -1, -1, -1],
            [2,  1, -1, -1],
            [3,  2,  1, -1],
            [0,  1,  2,  3],
            [4,  3,  2,  1]
        ]

        intcds = mol.internal_coordinates
        ccoords = mol.coords
        carts = ccoords.system
        internals = intcds.system


        raise Exception(nput.dihed_deriv(
            ccoords,
            0, 1, 2, 3
        ))

        new_jacs_anal, = ccoords.jacobian(internals, [1],
                                          mesh_spacing=1.0e-3,
                                          stencil=5,
                                          # all_numerical=True,
                                          analytic_deriv_order=1,
                                          converter_options=dict(strip_dummies=True)
                                          )

        raise Exception(new_jacs_anal.shape)
        raise Exception(new_jacs_anal[1][2][3], np.deg2rad(45))


        new_jacs_num, = ccoords.jacobian(internals, [1],
                                    mesh_spacing=1.0e-3,
                                    stencil=5,
                                    # all_numerical=True,
                                    analytic_deriv_order=0,
                                    converter_options=dict(strip_dummies=True)
                                   )


        raise Exception(new_jacs_num[1][2][3], np.deg2rad(45))



        raise Exception(new_jacs_num[1][2], new_jacs_anal[1][2])
```

#### <a name="EckartEmbedDipoles">EckartEmbedDipoles</a>
```python
    def test_EckartEmbedDipoles(self):
        scan_file = TestManager.test_data("tbhp_030.log")
        ref_file = TestManager.test_data("tbhp_180.fchk")

        scan = Molecule.from_file(scan_file)
        ref = Molecule.from_file(ref_file)
        sel = np.where(ref.masses>3)[0]
        pax_rot = ref.principle_axis_frame(sel=sel, inverse=True)  # type: MolecularTransformation
        rot_ref = pax_rot.apply(ref)

        transf = scan.eckart_frame(rot_ref, sel=sel)

        carts, dips = DipoleSurface.get_log_values(scan_file, keys=("StandardCartesianCoordinates", "OptimizedDipoleMoments"))
        rot_dips = np.array([ np.dot(t.transformation_function.transform, d) for t,d in zip(transf, dips) ])
        self.assertTrue(np.allclose(np.linalg.norm(dips, axis=1)-np.linalg.norm(rot_dips, axis=1), 0.))
```

#### <a name="EckartEmbedMolecule">EckartEmbedMolecule</a>
```python
    def test_EckartEmbedMolecule(self):

        ref_file = TestManager.test_data("tbhp_180.fchk")
        ref = Molecule.from_file(ref_file)
        new = ref.get_embedded_molecule()
```

#### <a name="EmbeddedMolecule">EmbeddedMolecule</a>
```python
    def test_EmbeddedMolecule(self):

        file_name = TestManager.test_data("HOH_freq.fchk")

        mol1 = Molecule.from_file(file_name)
        # init_mat1 = mol1.normal_modes.modes
        mol = mol1.get_embedded_molecule()
        init_mat = mol1.normal_modes.modes.basis.matrix
        self.assertTrue(np.allclose(mol.moments_of_inertia, mol1.moments_of_inertia),
                        msg="(HOH) Moments of inertia changed post-rotation: {} to {}".format(mol1.moments_of_inertia, mol.moments_of_inertia)
                        )
        self.assertTrue(np.allclose(mol.inertial_axes, np.eye(3)),
                        msg="(HOH) Principle axes are not identity matrix post-rotation: {}".format(mol.inertial_axes)
                        )

        norms_1 = np.linalg.norm(mol.normal_modes.modes.basis.matrix, axis=0)
        norms_2 = np.linalg.norm(init_mat, axis=0)
        self.assertTrue(np.allclose(norms_1, norms_2),
                        msg="(HOH) Normal modes renomalized:{} different from {}".format(norms_1, norms_2)
                        )

        # try on TBHP
        file_name = TestManager.test_data("tbhp_180.fchk")
        mol1 = Molecule.from_file(file_name)
        # init_mat1 = mol1.normal_modes.modes
        mol = mol1.get_embedded_molecule()
        init_mat = mol1.normal_modes.modes.basis.matrix
        self.assertTrue(np.allclose(mol.moments_of_inertia, mol1.moments_of_inertia),
                        msg="(TBHP) Moments of inertia changed post-rotation: {} to {}".format(mol1.moments_of_inertia,
                                                                                        mol.moments_of_inertia)
                        )
        self.assertTrue(np.allclose(mol.inertial_axes, np.eye(3)),
                        msg="(TBHP) Principle axes are not identity matrix post-rotation: {}".format(mol.inertial_axes)
                        )

        norms_1 = np.linalg.norm(mol.normal_modes.modes.basis.matrix, axis=0)
        norms_2 = np.linalg.norm(init_mat, axis=0)
        self.assertTrue(np.allclose(norms_1, norms_2),
                        msg="(TBHP) Normal modes renomalized: {} different from {}".format(norms_1, norms_2)
                        )


        # try on HOONO
        file_name = TestManager.test_data("HOONO_freq.fchk")
        mol1 = Molecule.from_file(file_name)
        # init_mat1 = mol1.normal_modes.modes
        mol = mol1.get_embedded_molecule()
        init_mat = mol1.normal_modes.modes.basis.matrix
        self.assertTrue(np.allclose(mol.moments_of_inertia, mol1.moments_of_inertia),
                        msg="(HOONO) Moments of inertia changed post-rotation: {} to {}".format(mol1.moments_of_inertia,
                                                                                        mol.moments_of_inertia)
                        )
        self.assertTrue(np.allclose(mol.inertial_axes, np.eye(3)),
                        msg="(HOONO) Principle axes are not identity matrix post-rotation: {}".format(mol.inertial_axes)
                        )

        norms_1 = np.linalg.norm(mol.normal_modes.modes.basis.matrix, axis=0)
        norms_2 = np.linalg.norm(init_mat, axis=0)
        self.assertTrue(np.allclose(norms_1, norms_2),
                        msg="(HOONO) Normal modes renomalized: {} different from {}".format(norms_1, norms_2)
        )
```

#### <a name="AddDummyAtoms">AddDummyAtoms</a>
```python
    def test_AddDummyAtoms(self):

        file_name = TestManager.test_data("HOONO_freq.fchk")

        mol = Molecule.from_file(file_name)
        n_pos = mol.atom_positions["N"]
        o_pos = mol.atom_positions["O"]

        normal = nput.vec_crosses(
            mol.coords[o_pos[0]] - mol.coords[o_pos[1]],
            mol.coords[n_pos[0]] - mol.coords[o_pos[1]],
            normalize=True
        )

        mol2 = mol.insert_atoms("X", mol.coords[o_pos[1]] + 5 * normal, 5, handle_properties=False)
        del mol # to elim hard to debug errors

        self.assertEquals(mol2.atoms,
                          ("H", "O", "O", "N", "O", "X")
                          )
        self.assertEquals(np.linalg.norm(mol2.coords[o_pos[1]] - mol2.coords[-1]), 5.0)

        mol2.zmatrix = [
            [1, -1, -1, -1], #O
            [2,  1, -1, -1], #O
            [3,  2,  1, -1], #N
            [5,  2,  1,  3], #X
            [0,  1,  2,  5], #H
            [4,  3,  2,  5], #O
        ]

        self.assertEquals(
            mol2.internal_coordinates[3, 0], 5.0
        )
        self.assertEquals(
            mol2.internal_coordinates[3, 1], np.pi/2
        )
        self.assertEquals(
            mol2.internal_coordinates[3, 2], np.pi/2
        )
        self.assertEquals(
            mol2.internal_coordinates[4, 2], np.pi/2
        )
        self.assertEquals(
            mol2.internal_coordinates[5, 2], -np.pi/2
        )
```

#### <a name="AddDummyAtomProperties">AddDummyAtomProperties</a>
```python
    def test_AddDummyAtomProperties(self):

        file_name = TestManager.test_data("HOONO_freq.fchk")

        mol = Molecule.from_file(file_name)
        n_pos = mol.atom_positions["N"]
        o_pos = mol.atom_positions["O"]

        normal = nput.vec_crosses(
            mol.coords[o_pos[0]] - mol.coords[o_pos[1]],
            mol.coords[n_pos[0]] - mol.coords[o_pos[1]],
            normalize=True
        )

        mol2 = mol.insert_atoms("X", mol.coords[o_pos[1]] + 5 * normal, 5, handle_properties=False)

        self.assertEquals(
            mol2.moments_of_inertia.tolist(),
            mol.moments_of_inertia.tolist()
        )

        self.assertEquals(
            mol2.inertial_axes.tolist(),
            mol.inertial_axes.tolist()
        )
```

#### <a name="AddDummyAtomJacobians">AddDummyAtomJacobians</a>
```python
    def test_AddDummyAtomJacobians(self):

        file_name = TestManager.test_data("HOONO_freq.fchk")

        mol = Molecule.from_file(file_name)
        n_pos = mol.atom_positions["N"]
        o_pos = mol.atom_positions["O"]

        normal = nput.vec_crosses(
            mol.coords[o_pos[0]] - mol.coords[o_pos[1]],
            mol.coords[n_pos[0]] - mol.coords[o_pos[1]],
            normalize=True
        )
        mol2 = mol.insert_atoms("X", mol.coords[o_pos[1]] + 5 * normal, 3, handle_properties=False)
        mol2.zmatrix = [
            [1, -1, -1, -1],  # O
            [2,  1, -1, -1],  # O
            [3,  2,  1, -1],  # N
            [5,  2,  1,  3],  # X
            [0,  1,  2,  5],  # H
            [4,  3,  2,  5],  # O
        ]

        jacobians_no_dummy = mol2.coords.jacobian(mol2.internal_coordinates.system,
                                                  [1, 2],
                                                  stencil=3,
                                                  all_numerical=True,
                                                  converter_options=dict(strip_dummies=True),
                                                  )
        self.assertEquals(jacobians_no_dummy[0].shape, (5, 3, 5, 3))
        self.assertEquals(jacobians_no_dummy[1].shape, (5, 3, 5, 3, 5, 3))
        jacobians = mol2.coords.jacobian(mol2.internal_coordinates.system,
                                                       [1, 2, 3],
                                                       stencil=5,
                                                       all_numerical=True,
                                                       converter_options=dict(strip_dummies=False),
                                                       )
        self.assertEquals(jacobians[0].shape, (6, 3, 6, 3))
        self.assertEquals(jacobians[1].shape, (6, 3, 6, 3, 6, 3))
        self.assertEquals(jacobians[2].shape, (6, 3, 6, 3, 6, 3, 6, 3))
        jacobians_analytic = mol2.coords.jacobian(mol2.internal_coordinates.system,
                                                       [1, 2],
                                                       stencil=5,
                                                       analytic_deriv_order=1,
                                                       converter_options=dict(strip_dummies=False),
                                                       )
        self.assertEquals(jacobians_analytic[0].shape, (6, 3, 6, 3))
        self.assertEquals(jacobians_analytic[1].shape, (6, 3, 6, 3, 6, 3))
        jacobians_no_dummy_analytic = mol2.coords.jacobian(mol2.internal_coordinates.system,
                                                  [1, 2],
                                                  stencil=3,
                                                  analytic_deriv_order=1,
                                                  converter_options=dict(strip_dummies=True),
                                                  )
        self.assertEquals(jacobians_no_dummy_analytic[0].shape, (5, 3, 5, 3))
        self.assertEquals(jacobians_no_dummy_analytic[1].shape, (5, 3, 5, 3, 5, 3))

        self.assertTrue(np.allclose(
            jacobians[0][0, 0][:2],
            jacobians_no_dummy[0][0, 0][:2]
        ))

        self.assertTrue(np.allclose(
            jacobians[1][0, 0, 0, 0][:2], jacobians_no_dummy[1][0, 0, 0, 0][:2]
        ))

        # with BlockProfiler():
        jacobians_no_dummy = mol2.internal_coordinates.jacobian(mol2.coords.system,
                                                       [1, 2],
                                                       stencil=3,
                                                       all_numerical=True,
                                                       converter_options=dict(strip_dummies=True),
                                                       )
        self.assertEquals(jacobians_no_dummy[0].shape, (5, 3, 5, 3))
        self.assertEquals(jacobians_no_dummy[1].shape, (5, 3, 5, 3, 5, 3))
        jacobians = mol2.internal_coordinates.jacobian(mol2.coords.system,
                                                       [1, 2],
                                                       stencil=3,
                                                       all_numerical=True,
                                                       converter_options=dict(strip_dummies=False),
                                                       )
        self.assertEquals(jacobians[0].shape, (6, 3, 6, 3))
        self.assertEquals(jacobians[1].shape, (6, 3, 6, 3, 6, 3))
```

#### <a name="InternalCoordOrder">InternalCoordOrder</a>
```python
    def test_InternalCoordOrder(self):
        file_name = TestManager.test_data("HOONO_freq.fchk")

        mol = Molecule.from_file(file_name)
        mol.zmatrix = [
            [1, -1, -1, -1],
            [2,  1, -1, -1],
            [3,  2,  1, -1],
            [0,  1,  2,  3],
            [4,  3,  2,  1]
        ]
        mol_ics = mol.internal_coordinates

        mol2 = Molecule.from_file(file_name)
        mol2.zmatrix = [
            [0, -1, -1, -1],  # H
            [1,  0, -1, -1],  # O
            [2,  1,  0, -1],  # O
            [3,  2,  1,  0],  # N
            [4,  3,  2,  0]   # O
        ]
        mol2_ics = mol2.internal_coordinates

        self.assertEquals(mol_ics[1, 0], mol2_ics[2, 0])
        self.assertEquals(mol_ics[3, 0], mol2_ics[1, 0])
        self.assertEquals(mol_ics[3, 2], mol2_ics[3, 2])

        jacs = mol.coords.jacobian(mol_ics.system, [1])[0]
        jacs2 = mol2.coords.jacobian(mol2_ics.system, [1])[0]

        self.assertEquals(jacs[0, 0][3, 0], jacs2[0, 0][1, 0])
        self.assertEquals(jacs[0, 0][1, 0], jacs2[0, 0][2, 0])
        self.assertEquals(jacs[0, 0][3, 2], jacs2[0, 0][3, 2])

        remade_carts = np.round(mol_ics.convert(mol.coords.system), 4)
        remade_carts2 = np.round(mol2_ics.convert(mol2.coords.system), 4)
        # raise Exception(remade_carts, remade_carts2)

        jacs = mol_ics.jacobian(mol.coords.system, [1], all_numerical=True)[0]
        jacs2 = mol2_ics.jacobian(mol2.coords.system, [1], all_numerical=True)[0]

        self.assertTrue(np.allclose(jacs[3, 0], jacs2[1, 0]))
```

#### <a name="Plotting">Plotting</a>
```python
    def test_Plotting(self):

        # g = Graphics3D(
        #     image_size=[1500, 1500],
        #     plot_range=[[-10, 10]]*3,
        #     backend="VTK"
        #     )
        # h5 = Molecule.from_file(
        #     self.test_log_h2,
        #     # self.test_fchk,
        #     # bonds = [
        #     #     [0, 1, 1],
        #     #     [0, 2, 1]
        #     # ]
        # )
        # h5.plot(
        #     figure=g
        #     # mode='3D',
        #     # bond_style= { "circle_points": 24 },
        #     # atom_style= { "sphere_points": 24 }
        # )
        m = Molecule.from_file(
            self.test_fchk,
            bonds = [
                [0, 1, 1],
                [0, 2, 1]
            ]
        )
        m.plot(
            # figure=g
            # mode='3D',
            # bond_style= { "circle_points": 24 },
            # atom_style= { "sphere_points": 24 }
            )
```

#### <a name="BondGuessing">BondGuessing</a>
```python
    def test_BondGuessing(self):
        m = Molecule.from_file(self.test_fchk)
        self.assertEquals(m.bonds, [[0, 1, 1], [0, 2, 1]])
```

#### <a name="Frags">Frags</a>
```python
    def test_Frags(self):
        m = Molecule.from_file(self.test_fchk)
        self.assertEquals(len(m.prop("fragments")), 1)
```

#### <a name="AutoZMat">AutoZMat</a>
```python
    def test_AutoZMat(self):
        raise NotImplementedError("saddy")
        m = Molecule.from_file(self.test_fchk)
```

#### <a name="HODModes">HODModes</a>
```python
    def test_HODModes(self):
        # oops fucked up getting D out
        m = Molecule.from_file(self.test_HOD, bonds=[[0, 1, 1], [0, 2, 1]])
        modes = m.normal_modes
        self.assertEquals(m.atoms, ("O", "H", "D"))
        self.assertEquals(
            tuple(np.round(modes.freqs*UnitsData.convert("Hartrees", "Wavenumbers"))),
            (1422.0, 2810.0, 3874.0)
        )
```

#### <a name="H2OModes">H2OModes</a>
```python
    def test_H2OModes(self):
        m = Molecule.from_file(self.test_fchk, bonds=[[0, 1, 1], [0, 2, 1]])
        modes = m.normal_modes
        self.assertEquals(m.atoms, ("O", "H", "H"))
        self.assertEquals(
            tuple(np.round(modes.freqs*UnitsData.convert("Hartrees", "Wavenumbers"))),
            (1622.0, 3803.0, 3938.0)
        )
```

#### <a name="RenormalizeGaussianModes">RenormalizeGaussianModes</a>
```python
    def test_RenormalizeGaussianModes(self):

        with GaussianFChkReader(self.test_HOD) as gr:
            parse = gr.parse(["Coordinates", "Gradient", "AtomicMasses",
                              "ForceConstants", "ForceDerivatives", "VibrationalModes", "VibrationalData"])

        coords = UnitsData.convert("Angstroms", "AtomicUnitOfLength") * parse["Coordinates"]
        masses = UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass") * parse["AtomicMasses"]
        modes = parse["VibrationalModes"].T
        freqs = parse["VibrationalData"]["Frequencies"]
        fcs = parse["ForceConstants"].array
        sad = UnitsData.convert("Hartrees", "Wavenumbers") * np.sqrt(np.diag(np.dot(np.dot(modes.T, fcs), modes)))
        modes = modes * freqs/sad
        print( UnitsData.convert("Hartrees", "Wavenumbers") * np.sqrt(np.diag(np.dot(np.dot(modes.T, fcs), modes))))

        masses = np.broadcast_to(masses, (len(masses), 3)).T.flatten()
        # print(modes-np.linalg.pinv(modes).T)
        print(np.dot(np.dot(modes.T, np.diag(masses)), modes))

        modes_2 = Molecule.from_file(self.test_HOD).get_normal_modes(normalize=False)
        mm = modes_2._basis.matrix

        print(np.dot(np.dot(mm.T, np.diag(masses)), mm))
        print(UnitsData.convert("Hartrees", "Wavenumbers") * np.sqrt(np.diag(np.dot(np.dot(mm.T, fcs), mm))))
```

#### <a name="VisualizeNormalModes">VisualizeNormalModes</a>
```python
    def test_VisualizeNormalModes(self):

        from Psience.Molecools.Vibrations import MolecularVibrations, MolecularNormalModes
        from McUtils.Plots import GraphicsGrid, Graphics3D

        m = Molecule.from_file(self.test_fchk, bonds = [[0, 1, 1], [0, 2, 1]])

        with GaussianFChkReader(self.test_fchk) as reader:
            parse = reader.parse(("VibrationalModes", "VibrationalData"))
        modes = parse["VibrationalModes"].T

        test_freqs = parse["VibrationalData"]["Frequencies"]

        nms = m.normal_modes
        realvibs = MolecularVibrations(m, basis=MolecularNormalModes(m, modes, freqs=test_freqs))

        realvibs.visualize(mode='jupyter') # get no bugs

        plot_vibrations = False
        if plot_vibrations:
            nmodes = 1
            mode_start = 0
            g = GraphicsGrid(nrows=2, ncols=nmodes,
                             graphics_class=Graphics3D,
                             plot_range = [[-2, 2], [-2, 2], [-2, 2]],
                             fig_kw = dict(figsize = (17, 5)),
                             tighten = True
                             )

            for i in range(nmodes):
                nms.visualize(step_size=.1, figure = g[0, i], which=mode_start + i,
                              anim_opts= dict(interval = 10)
                              )

            for i in range(nmodes):
                realvibs.visualize(step_size=.1, figure = g[1, i], which= mode_start+i,
                                   anim_opts= dict(interval = 10)
                                   )

            g.show()

        self.assertEquals(
            tuple(np.round(UnitsData.convert("Hartrees", "Wavenumbers")*nms.modes.freqs, 4)),
            tuple(np.round(test_freqs, 4))
        )
```

#### <a name="InternalCartesianJacobians">InternalCartesianJacobians</a>
```python
    def test_InternalCartesianJacobians(self):
        import McUtils.Plots as plt
        m = Molecule.from_file(TestManager.test_data('HOH_freq.fchk'),
                               zmatrix=[
                                   [0, -1, -1, -1],
                                   [1,  0, -1, -1],
                                   [2,  0,  1, -1]
                               ]
                               )
        # m = m.get_embedded_molecule()
        intcds = m.internal_coordinates
        carts = m.coords
        # ijacsnum, ijacs2num = intcds.jacobian(carts.system, [1, 2], analytic_deriv_order=0, mesh_spacing=1.0e-2)
        ijacsnum, ijacs2num = intcds.jacobian(carts.system, [1, 2], all_numerical=True, mesh_spacing=1.0e-2)
        ijacs, ijacs2 = intcds.jacobian(carts.system, [1, 2], analytic_deriv_order=1,
                                        converter_options=dict(reembed=False)
                                        )#, mesh_spacing=1.0e-2)
        jacs, jacs2 = carts.jacobian(intcds.system, [1, 2], mesh_spacing=1.0e-5)

        meh1 = ijacs.squeeze().reshape(9, 9)
        meh0 = ijacsnum.squeeze().reshape(9, 9)
        meh2 = jacs.squeeze().reshape(9, 9)

        itest = np.dot(meh1, meh2)
        itest2 = np.dot(meh2, meh1)

        # plt.ArrayPlot(meh1)
        # plt.ArrayPlot(meh1)
        # plt.ArrayPlot(meh0).show()
        # plt.ArrayPlot(np.dot(meh0, meh2)).show()
        self.assertTrue(np.allclose(np.eye(9), itest))


        good_sel = (...,) + np.ix_((3, 5, 6), (3, 5, 6))
        meh12 = ijacs2.squeeze().reshape(9, 9, 9)
        meh12 = meh12.transpose(2, 0, 1).reshape(3, 3, 9, 9)
        meh22 = ijacs2num.squeeze().reshape(9, 9, 9)
        meh22 = meh22.transpose(2, 0, 1).reshape(3, 3, 9, 9)
        meh12 = meh12[good_sel]
        meh22 = meh22[good_sel]
        ps = dict(vmin=-.05, vmax=.05)
        plt.TensorPlot(meh12, plot_style=ps)
        plt.TensorPlot(meh22, plot_style=ps).show()
        # plt.TensorPlot(meh22-meh12, plot_style=ps).show()

        self.assertAlmostEquals(meh22[1, 1, 0, 0], .009235, places=6)
        self.assertTrue(np.allclose(meh12, meh22))
```

#### <a name="CompositeCoordinates">CompositeCoordinates</a>
```python
    def test_CompositeCoordinates(self):
        def conv(r, t, f, **kwargs):
            return [r**2, np.cos(t), np.sin(f)]
        def inv(r2, t, f, **kwargs):
            return [np.sqrt(r2), np.arccos(t), np.arcsin(f)]

        mol = Molecule.from_file(
            TestManager.test_data('HOONO_freq.fchk'),
            internals = {
                'zmatrix':[
                    [1, -1, -1, -1],
                    [2,  1, -1, -1],
                    [3,  2,  1, -1],
                    [0,  1,  2,  3],
                    [4,  3,  2,  1]
                ],
                'conversion':conv,
                'inverse':inv,
                'converter_options':{'pointwise':True}
            }
        )

        mol2 = Molecule.from_file(
            TestManager.test_data('HOONO_freq.fchk'),
            internals = {
                'zmatrix':[
                    [1, -1, -1, -1],
                    [2,  1, -1, -1],
                    [3,  2,  1, -1],
                    [0,  1,  2,  3],
                    [4,  3,  2,  1]
                ]
            }
        )

        ic1 = mol.internal_coordinates
        ic2 = mol2.internal_coordinates

        self.assertAlmostEquals(np.sum(ic1.convert(ic2.system)-ic2)[()], 0.)
```

#### <a name="RDKitSpectrum">RDKitSpectrum</a>
```python
    def test_RDKitSpectrum(self):
        # propane = Molecule.from_string('CCC', 'smi', energy_evaluator='rdkit').optimize()
        # propane.get_harmonic_spectrum().plot().show()

        propanol = Molecule.from_string('CCCO', 'smi',
                                        energy_evaluator='rdkit',
                                        charge_evaluator='aimnet2'
                                        ).optimize()
        propanol.get_harmonic_spectrum().plot(
            plot_range=[
                None,
                [0, 170]
            ],
            axes_labels=['Freq. (cm$^{-1}$)', 'Int. (km mol$^{-1}$)'],
            plot_label='AIMNet2 Charges w/ MMFF Modes',
            padding=[[55, 0], [40, 20]]
        ).savefig(os.path.expanduser('~/Desktop/aimnet2_propanol_rdkit_opt.png'))

        propanol = propanol.modify(charge_evaluator='rdkit')
        propanol.get_harmonic_spectrum().plot(
            plot_range=[
                None,
                [0, 170]
            ],
            axes_labels=['Freq. (cm$^{-1}$)', 'Int. (km mol$^{-1}$)'],
            plot_label='Gasteiger Charges w/ MMFF Modes',
            padding=[[55, 0], [40, 20]]
        ).savefig(os.path.expanduser('~/Desktop/gasteiger_propanol_rdkit_opt.png'))

        return
        # water = Molecule.from_file(TestManager.test_data("HOH_freq.fchk"))

        water = Molecule.from_string('O', 'smi',
                                     energy_evaluator='rdkit',
                                     charge_evaluator='aimnet2'
                                     ).optimize()
        water.get_harmonic_spectrum().plot()

        water = water.modify(charge_evaluator='rdkit')
        spec = water.get_harmonic_spectrum()
        spec.plot().show()
```

#### <a name="ExpansionPotential">ExpansionPotential</a>
```python
    def test_ExpansionPotential(self):
        h2co = Molecule.from_file(TestManager.test_data('OCHH_freq.fchk'))
        disps, scan_coords = h2co.get_scan_coordinates(
            [[-.5, .5, 1000]],
            which=[[2, 2]],
            return_displacements=True
        )
        [grad, hess, cubes, quarts] = h2co.potential_derivatives
        vals_2 = h2co.calculate_energy(scan_coords,
                                       evaluator=(
                                           'expansion',
                                           {
                                               'expansion': [grad, hess]
                                           }
                                       )
                                       ) * 219475.6
        vals_3 = h2co.calculate_energy(scan_coords,
                                     evaluator=(
                                         'expansion',
                                         {
                                             'expansion':[grad, hess, cubes]
                                         }
                                     )
                                     ) * 219475.6
        vals_4 = h2co.calculate_energy(scan_coords,
                                     evaluator=(
                                         'expansion',
                                         {
                                             'expansion':[grad, hess, cubes, quarts]
                                         }
                                     )
                                     ) * 219475.6
        # h2co.plot(scan_coords).show()
        base_fig = plt.Plot(disps, vals_2, plot_range=[None, [0, 20000]])
        plt.Plot(disps, vals_3, figure=base_fig)
        plt.Plot(disps, vals_4, figure=base_fig)
        base_fig.show()
```

#### <a name="OpenBabel">OpenBabel</a>
```python
    def test_OpenBabel(self):
        mol = Molecule.from_file(TestManager.test_data("nh3.fchk"))
        print(mol.to_string("pdb"))
        return
```

#### <a name="1DPotentialReps">1DPotentialReps</a>
```python
    def test_1DPotentialReps(self):
        ochh = self.setup_OCHH(optimize=True)
        int_ochh = ochh.modify(internals=[
            [0, -1, -1, -1],
            [1, 0, -1, -1],
            [2, 1, 0, -1],
            [3, 1, 2, 0],
        ])

        scan_disps = [-1.0, 1.0, 51]
        scan_angles = int_ochh.get_scan_coordinates(
            [scan_disps],
            which=[[3, 1]],
            internals='reembed'
        )
        # scan_disp_dist = [-.4, .7, 25]
        # # int_ochh.plot(scan_angles).show()
        # scan_dists = int_ochh.get_scan_coordinates(
        #     [scan_disp_dist],
        #     which=[[1, 0]],
        #     internals='reembed'
        # )

        pot_vals = int_ochh.calculate_energy(scan_angles)

        # pot_vals_dists = int_ochh.calculate_energy(scan_dists)

        # woof = ochh.get_anharmonic_parameters(
        #     [(0, 1), (1, 2), (1, 3), (2, 1, 3), (0, 1, 2, 3)]
        # )
        woof_pots = ochh.get_1d_potentials(
            [(0, 1), (1, 2), (1, 3), (2, 1, 3), (0, 1, 2, 3)]
        )
        woof_pots_4 = ochh.get_1d_potentials(
            [(0, 1), (1, 2), (1, 3), (2, 1, 3), (0, 1, 2, 3)],
            poly_expansion_order=4
        )
        woof_pots_morse = ochh.get_1d_potentials(
            [(0, 1), (1, 2), (1, 3), (2, 1, 3), (0, 1, 2, 3)],
            quartic_potential_cutoff=0
        )
        # for method, (params, re) in woof:
        #     print(
        #         [
        #             p * 219475.6
        #             for n, p in enumerate(params)
        #         ],
        #         re
        #     )

        # _, pot_vals_dists_appx = woof_pots[0](
        #     scan_dists
        # )
        angs, pot_vals_angs_appx = woof_pots[3](
            scan_angles
        )
        angs, pot_vals_angs_appx_4 = woof_pots_4[3](
            scan_angles
        )
        angs, pot_vals_angs_appx_morse = woof_pots_morse[3](
            scan_angles
        )

        # uh = nput.internal_coordinate_tensors(
        #     scan_angles,
        #     [
        #         [0, 1],
        #         [1, 2],
        #         [1, 3],
        #         [2, 1, 3]
        #     ],
        #     order=0
        # )
        # print(uh[0])

        disp_vals = np.linspace(*scan_disps)
        ploots = plt.Plot(disp_vals, pot_vals * 219475.6, color='black', plot_range=[None, [None, 35000]])
        # plt.Plot(np.linspace(*scan_disp_dist), pot_vals_dists * 219475.6, figure=ploots)
        # plt.Plot(np.linspace(*scan_disp_dist), pot_vals_dists_appx[0] * 219475.6, figure=ploots,
        #          linestyle='dashed')
        # k = len(disp_vals) // 2
        # f2_alt = (pot_vals[k+1] + pot_vals[k-1] - 2*pot_vals[k]) / (2*(disp_vals[1] - disp_vals[0])**2)
        # print(f2_alt)
        # print(angs[0] - angs[0][k])
        # print(disp_vals)
        shang_vals = (ochh.calculate_energy() + pot_vals_angs_appx[0])
        shang_vals_4 = (ochh.calculate_energy() + pot_vals_angs_appx_4[0])
        shang_vals_morse = (ochh.calculate_energy() + pot_vals_angs_appx_morse[0])
        plt.Plot(np.linspace(*scan_disps), shang_vals * 219475.6, figure=ploots,
                 linestyle='dashed')
        plt.Plot(np.linspace(*scan_disps), shang_vals_4 * 219475.6, figure=ploots,
                 linestyle='dashed')
        plt.Plot(np.linspace(*scan_disps), shang_vals_morse * 219475.6, figure=ploots,
                 linestyle='dashed')
        ploots.show()
        return
```

#### <a name="Constructors">Constructors</a>
```python
    def test_Constructors(self):

        a = Molecule.construct(TestManager.test_data('OCHH_freq.fchk'))
        b = Molecule.construct([
            a.atoms, a.coords,
            dict(internals=[
                [0, -1, -1, -1],
                [1,  0, -1, -1],
                [2,  1,  0, -1],
                [3,  1,  0,  2]
            ])
        ])
        c = Molecule.construct([b.atoms, (b.internals['zmatrix'], b.internal_coordinates[1:])])
        b = Molecule.construct('formaldehyde')
        c = Molecule.construct('OC')
```

#### <a name="InternalConv">InternalConv</a>
```python
    def test_InternalConv(self):

        gggg = Molecule(
            ['C', 'C', 'H', 'H', 'H', 'C', 'C', 'H', 'H', 'H', 'C', 'C', 'H', 'H', 'H', 'C', 'H', 'H', 'H'],
            [[0., 0., 0.],
             [2.5221177, 0., 0.],
             [-1.074814, 1.74867225, 0.],
             [-1.06173943, -1.75324789, -0.02644452],
             [3.50660496, -1.80459232, -0.09993741],
             [4.12940149, 2.25082363, 0.07118096],
             [3.53603628, 4.43884306, 1.1781188],
             [5.96989258, 2.06662578, -0.83065536],
             [4.82076346, 6.03519987, 1.12058693],
             [1.77236274, 4.68977978, 2.19737026],
             [0.05431618, 1.81919915, 6.66571583],
             [0.10361848, 4.11457325, 7.68790077],
             [1.66160335, 1.19828653, 5.53987418],
             [-1.46076767, 4.82280651, 8.81630128],
             [1.71219257, 5.36239542, 7.43802933],
             [-2.06783807, -0.02189412, 6.91272006],
             [-2.80613635, -0.54136474, 5.04934629],
             [-1.42168595, -1.77751975, 7.80094268],
             [-3.61859877, 0.74452369, 8.04209746]],
            internals={
                'primitives': 'auto',
                'nonredundant_coordinates': [(15, 10, 11, 13)]
            }
        )

        raise Exception(gggg.internal_coordinates.shape)


        nh3 = Molecule.from_file(
            TestManager.test_data("nh3.fchk"),
            internals=[
                [0, -1, -1, -1],
                [1, 0, -1, -1],
                [2, 0, 1, -1],
                [3, 0, 1, 2]
            ]
        )
        emb_nh3 = nh3.get_embedded_molecule()
        emb_test = emb_nh3.internal_coordinates.convert(emb_nh3.coords.system)
        # conv_1 = nh3.internal_coordinates.convert(nh3.coords.system)
        # print(emb_nh3.internal_coordinates)
        # print(emb_nh3.internal_coordinates.converter_options)
        # print(np.round(emb_nh3.coords, 8))
        # print(np.round(emb_test, 8))
        # raise Exception(...)
        # conv_2 = ...

        nh3 = Molecule.from_file(
            TestManager.test_data("nh3.fchk"),
            internals=[
                [0, -1, -1, -1],
                [1, 0, -1, -1],
                [2, 0, 1, -1],
                [3, 0, 1, 2]
            ]
        )
        ders1 = nh3.get_cartesians_by_internals(1, strip_embedding=True, reembed=True)[0]
        ders_inv1 = nh3.get_internals_by_cartesians(1, strip_embedding=True)[0]

        # nh3 = Molecule.from_file(
        #     TestManager.test_data("nh3.fchk"),
        #     internals=[
        #         [0, -1, -1, -1],
        #         [1,  0, -1, -1],
        #         [2,  0,  1, -1],
        #         [3,  0,  1,  2]
        #     ]
        # )
        ders2 = nh3.get_cartesians_by_internals(1, method='classic', strip_embedding=True, reembed=True)[0]
        # print(np.round(ders1, 7)[0])
        # print(np.round(ders2, 7)[0])

        # print(np.round(ders2 @ ders_inv1, 6))
        # print(np.round(ders1 @ ders_inv1, 8))

        nh3_derivs_internal = nput.tensor_reexpand(
            nh3.get_cartesians_by_internals(2, strip_embedding=True, reembed=True),
            [0, nh3.potential_derivatives[1]]
        )

        nh3_gmatrix = nh3.g_matrix

        freqs_int, _ = scipy.linalg.eigh(nh3_derivs_internal[1], nh3_gmatrix, type=2)
        freqs_cart, _ = scipy.linalg.eigh(nh3.potential_derivatives[1], nh3.get_gmatrix(use_internals=False), type=2)

        print(np.sqrt(freqs_int) * UnitsData.convert("Hartrees", "Wavenumbers"))
        print(np.sqrt(freqs_cart[6:]) * UnitsData.convert("Hartrees", "Wavenumbers"))
```

#### <a name="AutomaticConversion">AutomaticConversion</a>
```python
    def test_AutomaticConversion(self):
        sys = 'benzene.sdf'
        mol = Molecule.from_file(
            TestManager.test_data(sys),
            internals='auto'
        )

        with BlockProfiler():
            disps = mol.get_displaced_coordinates(
                np.linspace(-1, 1, 25)[:, np.newaxis],
                which=[0],
                use_internals=True
            )
```

#### <a name="FastInternals">FastInternals</a>
```python
    def test_FastInternals(self):

        sys = 'nh3.fchk'
        mol = Molecule.from_file(
            TestManager.test_data(sys),
            internals=[
                (1, 0),
                (2, 0),
                (0, 1, 2),
                (0, 3),
                (0, 1, 3),
                (3, 0, 1, 2)
            ]
        )
        # plt, _, _ = mol.plot(backend='vpython')
        # plt.show()
        #
        # raise Exception(...)

        mol2 = Molecule.from_file(
            TestManager.test_data(sys),
            internals=[
                [0, -1, -1, -1],
                [1,  0, -1, -1],
                [2,  0,  1, -1],
                [3,  0,  1,  2]
            ]
        )

        # disp_carts = mol.get_displaced_coordinates([.2], [3], use_internals=True).convert(
        #     mol.coords.system
        # )
        # raise Exception([
        #     s.shape for s in mol.get_cartesians_by_internals(2)
        # ])

        # raise Exception([
        #     s.shape for s in mol2.get_cartesians_by_internals(2, strip_embedding=True)
        # ])

        # with ...:
        #     mol.get_cartesians_by_internals(1)
        # int = mol.get_internals_by_cartesians(2)
        # int2 = mol2.get_internals_by_cartesians(2, strip_embedding=True)
        #
        # print(int[0][0])
        # print("-"*10)
        # print(int2[0][0])
        #
        # print("="*10)
        #
        # print(int[1][0, 0])
        # print("-"*10)
        # print(int2[1][0, 0])
        # raise Exception(...)

        cart = mol.get_cartesians_by_internals(1)[0][0]
        # cart2 = mol.get_cartesians_by_internals(1, method='og')[0][0]
        # print(cart2 / np.linalg.norm(cart2))
        cart3 = mol2.get_cartesians_by_internals(1, strip_embedding=True)[0][0]

        # print()
        # print(cart / np.linalg.norm(cart))
        # print(cart3 / np.linalg.norm(cart3))
        # print(
        #     (np.abs(cart) > 1e-14) * (cart) / (cart3)
        # )
        #
        # raise Exception(...)

        raise Exception(
            [
                np.round(s1 - s2, 6) for s1, s2 in zip(
                    mol.get_cartesians_by_internals(2),
                    mol2.get_cartesians_by_internals(2, method='classic', strip_embedding=True),
                    # mol2.get_cartesians_by_internals(2, strip_embedding=True)
                )
            ]
        )

        raise Exception(mol.internal_coordinates, mol.coords - disp_carts)
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
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Molecools.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Molecools.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Molecools.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Molecools.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Molecools/__init__.py#L1?message=Update%20Docs)   
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