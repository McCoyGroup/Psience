import itertools
import os.path
import pprint

from Peeves.TestUtils import *
from unittest import TestCase
from Peeves import BlockProfiler

from Psience.Molecools import Molecule, MolecularNormalModes
# from Psience.Molecools.Transformations import MolecularTransformation
from Psience.Data import DipoleSurface # this will be leaving Zachary very soon I think...
from McUtils.GaussianInterface import GaussianFChkReader, GaussianLogReader
from McUtils.Plots import *
import McUtils.Plots as plt
from McUtils.Coordinerds import cartesian_to_zmatrix
from McUtils.Data import UnitsData
import numpy as np, scipy
import McUtils.Numputils as nput
import McUtils.Profilers as prof
from McUtils.Formatters import TableFormatter
import McUtils.Formatters as mfmt

class MolecoolsTests(TestCase):
    def setUp(self):
        self.test_log_water = TestManager.test_data("water_OH_scan.log")
        self.test_log_freq = TestManager.test_data("water_freq.log")
        self.test_HOD = TestManager.test_data("HOD_freq.fchk")
        self.test_fchk = TestManager.test_data("water_freq.fchk")
        self.test_log_h2 = TestManager.test_data("outer_H2_scan_new.log")
        # np.seterr(all='ignore')

    def tearDown(self):
        ...

    @validationTest
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

    @validationTest
    def test_MolecularGMatrix(self):
        mol = Molecule.from_file(self.test_fchk)
        mol.zmatrix = [
            [0, -1, -1, -1],
            [1,  0, -1, -1],
            [2,  0,  1, -1]
        ]
        g = mol.g_matrix

        self.assertEquals(g.shape, (3, 3))

    @validationTest
    def test_ImportMolecule(self):

        n = 3 # water
        m = Molecule.from_file(self.test_fchk)
        self.assertEquals(m.atoms, ("O", "H", "H"))

    @inactiveTest
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

    @validationTest
    def test_EckartEmbed(self):
        m = Molecule.from_file(TestManager.test_data('HOH_freq.fchk'))
        crd = m.embed_coords(m.coords)
        self.assertTrue(np.allclose(m.coords, crd))

    @validationTest
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

    @inactiveTest
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

    @inactiveTest
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

        # ### Visualize dipole surface
        # dists = np.linalg.norm(carts[1:, 5] - carts[1:, 6], axis=1)
        # Graphics.default_style['image_size'] = 575
        # g = GraphicsGrid(nrows=1, ncols=2, padding=((.075, 0), (0, .45)))
        # p = Plot(dists, rot_dips[:, 0], figure=g[0, 0])
        # Plot(dists, rot_dips[:, 1], figure=p)
        # Plot(dists, rot_dips[:, 2], figure=p)
        # p2= Plot(dists, dips[:, 0], figure=g[0, 1])
        # Plot(dists, dips[:, 1], figure=p2)
        # Plot(dists, dips[:, 2], figure=p2)
        # g.show()

    @validationTest
    def test_EckartEmbedMolecule(self):

        ref_file = TestManager.test_data("tbhp_180.fchk")
        ref = Molecule.from_file(ref_file)
        new = ref.get_embedded_molecule()

    @validationTest
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

        # raise Exception(mol1.coords, mol1.normal_modes.modes.basis.matrix.T)

        # raise Exception(mol.coords, mol.normal_modes.modes.basis.matrix.T)

    @validationTest
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
        del mol # to elim hard to find errors

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

    @validationTest
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

    @validationTest
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

        # # now test against regular stuff
        # mol = Molecule.from_file(file_name)
        #

    @validationTest
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

    @validationTest
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
        # g.show()

    @inactiveTest
    def test_BondGuessing(self):
        m = Molecule.from_file(self.test_fchk)
        self.assertEquals(m.bonds, [[0, 1, 1], [0, 2, 1]])

    @inactiveTest
    def test_Frags(self):
        m = Molecule.from_file(self.test_fchk)
        self.assertEquals(len(m.prop("fragments")), 1)

    @inactiveTest
    def test_AutoZMat(self):
        raise NotImplementedError("saddy")
        m = Molecule.from_file(self.test_fchk)

    @validationTest
    def test_HODModes(self):
        # oops fucked up getting D out
        m = Molecule.from_file(self.test_HOD, bonds=[[0, 1, 1], [0, 2, 1]])
        modes = m.normal_modes
        self.assertEquals(m.atoms, ("O", "H", "D"))
        self.assertEquals(
            tuple(np.round(modes.freqs*UnitsData.convert("Hartrees", "Wavenumbers"))),
            (1422.0, 2810.0, 3874.0)
        )

    @validationTest
    def test_H2OModes(self):
        m = Molecule.from_file(self.test_fchk, bonds=[[0, 1, 1], [0, 2, 1]])
        modes = m.normal_modes
        self.assertEquals(m.atoms, ("O", "H", "H"))
        self.assertEquals(
            tuple(np.round(modes.freqs*UnitsData.convert("Hartrees", "Wavenumbers"))),
            (1622.0, 3803.0, 3938.0)
        )

    @inactiveTest
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
        # print(modes._basis.matrix.T.dot(m.force_constants).shape)

    @validationTest
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

    @inactiveTest
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

    @validationTest
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

    @validationTest
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

    @validationTest
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
        # print(vals)

    @classmethod
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

    @validationTest
    def test_OpenBabel(self):
        mol = Molecule.from_file(TestManager.test_data("nh3.fchk"))
        print(mol.to_string("pdb"))
        return

    @validationTest
    def test_MethanolTorsionScan(self):
        import McUtils.Coordinerds as coordops

        methanol_zmatrix = coordops.reindex_zmatrix(
            coordops.functionalized_zmatrix(
                3,
                {
                    (2, 1, 0): [
                        [0, -1, -2, -3],
                        [1, -1, 0, -2],
                        [2, -1, 0, 1],
                    ]
                }
            ),
            [5, 1, 0, 2, 3, 4]
        )
        methanol = Molecule(
            ['C', 'O', 'H', 'H', 'H', 'H'],
            [[-0.6988896,  0.00487717, 0.00528378],
                    [ 1.69694605, -1.08628154, -0.22505594],
                    [-1.27384807, -0.22559494, 1.98568702],
                    [-0.59371792,  2.01534286, -0.59617633],
                    [-2.04665278, -0.99128091, -1.21504054],
                    [ 2.91616232,  0.28293736, 0.04530201]],
            energy_evaluator='aimnet2'
        )


        # # rpnms = methanol.get_reaction_path_modes(zero_gradient_cutoff=zero_gradient_cutoff=.1)
        # rpnms = methanol.get_reaction_path_modes()
        # print(rpnms.freqs * 219474.6)
        # # print(np.round(rpnms.coords_by_modes @ rpnms.modes_by_coords, 4))
        #
        # return


        methanol = methanol.optimize()

        # me_opt_scipy = methanol_1.optimize(mode='scipy')
        # methanol = methanol_1.optimize()

        # print(me_opt_scipy.coords)

        # print(methanol.calculate_energy() - me_opt_scipy.calculate_energy())
        # return

        #.optimize(mode='scipy')
        me_int = methanol.modify(
            internals=methanol_zmatrix,
        )
        scan_spec = [-np.deg2rad(180), np.deg2rad(180), 72]
        me_coords = me_int.get_scan_coordinates(
            [scan_spec],
            which=[11],
            internals='reembed'
        )
        me_pot = methanol.get_energy_function()(me_coords)

        # me_aimnet2 = methanol.get_energy_function()
        # me_pot = np.array([me_aimnet2(x) for x in me_coords])

        # base_pot = plt.Plot(
        #     np.rad2deg(np.linspace(*scan_spec)),
        #     (me_pot - methanol.calculate_energy()) * UnitsData.convert("Hartrees", "Wavenumbers")
        # )
        # subpot = plt.Plot(
        #     [-60, -60],
        #     [0, 500],
        #     linestyle='dashed',
        #     figure=base_pot,
        #     color='green'
        # )
        # subpot = plt.Plot(
        #     [60, 60],
        #     [0, 500],
        #     linestyle='dashed',
        #     figure=base_pot,
        #     color='green'
        # )
        #
        # base_pot.show()

    @validationTest
    def test_MethanolConstrainedOpt(self):
        import McUtils.Coordinerds as coordops

        # methanol = Molecule.from_string(
        #     'methanol',
        #     # energy_evaluator='aimnet2'
        # )


        methanol = Molecule(
            ['C', 'O', 'H', 'H', 'H', 'H'],
            [[-0.6988896 ,  0.00487717,  0.00528378],
             [ 1.69694605, -1.08628154, -0.22505594],
             [-1.27384807, -0.22559494,  1.98568702],
             [-0.59371792,  2.01534286, -0.59617633],
             [-2.04665278, -0.99128091, -1.21504054],
             [ 2.91616232,  0.28293736,  0.04530201]],
            # energy_evaluator='aimnet2'
        )

        eng0 = methanol.calculate_energy()

        methanol_zmatrix = coordops.reindex_zmatrix(
            coordops.functionalized_zmatrix(
                3,
                {
                    (2, 1, 0): [
                        [0, -1, -2, -3],
                        [1, -1, 0, -2],
                        [2, -1, 0, 1],
                    ]
                }
            ),
            [5, 1, 0, 2, 3, 4]
        )

        meth_int = methanol.modify(
            # evaluator=lambda c:
            # energy_evaluator='aimnet2',
            internals=methanol_zmatrix,
        )
        ints0 = meth_int.internal_coordinates
        meth_int = meth_int.optimize(max_displacement=.5, max_iterations=500,
                   coordinate_constraints=[(0,1)],
                   # logger=True
                   )
        eng1 = meth_int.calculate_energy()
        ints1 = meth_int.internal_coordinates

        # meth_opt = methanol.optimize(max_displacement=.5, max_iterations=500,
        #                              coordinate_constraints=[(0, 1)],
        #                              # logger=True,
        #                              # method='scipy'
        #                              )
        # eng2 = meth_opt.calculate_energy()
        # ints2 = meth_opt.modify(
        #     internals=methanol_zmatrix,
        # ).internal_coordinates
        # ints3 = huh2.modify(
        #     internals=methanol_zmatrix,
        # ).internal_coordinates
        #
        # # print(methanol.coords)
        # print(eng1, eng2, eng0, eng2 - eng1)
        # print(ints0)
        # print(ints1)
        # print(ints2)
        # print(ints3)

    @validationTest
    def test_ProjectedLocalModes(self):
        import McUtils.Coordinerds as coordops

        # methanol = Molecule.from_string(
        #     'methanol',
        #     # energy_evaluator='aimnet2'
        # )

        methanol = Molecule(
            ['C', 'O', 'H', 'H', 'H', 'H'],
            [[-0.6988896, 0.00487717, 0.00528378],
             [1.69694605, -1.08628154, -0.22505594],
             [-1.27384807, -0.22559494, 1.98568702],
             [-0.59371792, 2.01534286, -0.59617633],
             [-2.04665278, -0.99128091, -1.21504054],
             [2.91616232, 0.28293736, 0.04530201]],
            energy_evaluator='rdkit'
        ).optimize(max_iterations=50)
        # methanol.optimize(mode='scipy')
        nms = methanol.get_normal_modes(zero_freq_cutoff=1e-4)

        def do_the_thing(lms, freqs):
            print(len(freqs))
            disps = []
            print(freqs * 219474.56)
            for i in range(len(freqs)):
                scans = methanol.get_scan_coordinates(
                    [[-75, 75, 3]],
                    which=[i],
                    modes=lms
                )
                disps.append(
                    np.diff(nput.internal_coordinate_tensors(scans, [(0, 1)], order=0)[0][:, 0], axis=0)
                )
            print(np.round(np.array(disps), 4))

        print()
        # lms = nms.localize(coordinate_constraints=[(0,1)], orthogonal_projection=False)
        lfs = nms.freqs
        do_the_thing(nms, lfs)


        print()
        lms = nms.localize(coordinate_constraints=[(0,1)], orthogonal_projection=False)
        lfs = lms.local_freqs
        do_the_thing(lms, lfs)

        print()
        lms = nms.localize(coordinate_constraints=[(0,1)], orthogonal_projection=True)
        lfs = lms.local_freqs
        do_the_thing(lms, lfs)
        # print(nms.localize(atoms=[0,1]).local_freqs * 219474.56)

        print()
        lms = nms.localize(internals=[(0, 1)])
        lfs = lms.local_freqs
        do_the_thing(lms, lfs)

        print()
        lms = nms.localize(internals=[(0, 1)], projection=True, orthogonal_projection=False)
        lfs = lms.local_freqs
        do_the_thing(lms, lfs)

        print()
        lms = nms.localize(internals=[(0, 1)], projection=True, orthogonal_projection=True)
        lfs = lms.local_freqs
        do_the_thing(lms, lfs)

    @validationTest
    def test_ProjectedConstrainedModes(self):
        methanol = Molecule(
            ['C', 'O', 'H', 'H', 'H', 'H'],
            [[-0.6988896, 0.00487717, 0.00528378],
             [1.69694605, -1.08628154, -0.22505594],
             [-1.27384807, -0.22559494, 1.98568702],
             [-0.59371792, 2.01534286, -0.59617633],
             [-2.04665278, -0.99128091, -1.21504054],
             [2.91616232, 0.28293736, 0.04530201]],
            energy_evaluator='rdkit'
        ).optimize(max_iterations=50)

        # fixed_coord = [(0, 1, 5)]
        fixed_coord = [(2, 0, 1, 5)]

        # b, i, _ = nput.internal_basis(methanol.coords, fixed_coord,
        #                               masses=methanol.atomic_masses)
        # import McUtils.McUtils.Numputils.CoordOps as coop
        # ii = coop.fixed_angle_basis(methanol.coords, *fixed_coord[0])
        # _, i = nput.translation_rotation_eigenvectors(methanol.coords,
        #                                               masses=methanol.atomic_masses,
        #                                               mass_weighted=False)
        # i = [nput.find_basis(i, method='svd')]
        # i = [i]
        # def do_the_other_thing(expansion):
        #     print()
        #     # print(len(freqs))
        #     disps = []
        #     # print(freqs * 219474.56)
        #     for i in range(expansion[0].shape[0]):
        #         scans = methanol.get_scan_coordinates(
        #             [[-.5, .5, 3]],
        #             which=[i],
        #             coordinate_expansion=expansion
        #         )
        #         ints = nput.internal_coordinate_tensors(scans, fixed_coord, order=0)[0][:, 0]
        #         diffs = np.diff(ints, axis=0)
        #         disps.append(ints)
        #     print(np.round(np.array(disps), 4))

        # do_the_other_thing([i[0].T])
        # return

        # methanol.animate_coordinate(1, coordinate_expansion=[i[0].T])

        # b, i = nput.internal_coordinate_tensors(methanol.coords, [(0, 1, 2)], return_inverse=True)
        # return
        # methanol.animate_coordinate(0, coordinate_expansion=[i.T])

        # methanol.optimize(mode='scipy')
        nms = methanol.get_normal_modes(zero_freq_cutoff=1e-4)#.make_mass_weighted()
        cnms = nms.apply_constraints(fixed_coord, orthogonal_projection=True)
        # print(nms.local_freqs * UnitsData.hartrees_to_wavenumbers)
        # print(cnms.local_freqs * UnitsData.hartrees_to_wavenumbers)

        def do_the_thing(lms, freqs):
            print()
            print(len(freqs))
            disps = []
            print(freqs * 219474.56)
            for i in range(len(freqs)):
                scans = methanol.get_scan_coordinates(
                    [[-75, 75, 3]],
                    which=[i],
                    modes=lms
                )
                disps.append(
                    np.diff(nput.internal_coordinate_tensors(scans, fixed_coord, order=0)[0][:, 0], axis=0)
                )
            print(np.round(np.array(disps), 4))

        do_the_thing(nms, nms.local_freqs)
        do_the_thing(cnms, cnms.local_freqs)

        # scans = methanol.get_scan_coordinates(
        #     [[-100, 100, 3]],
        #     which=[4],
        #     modes=cnms
        # )
        # print(nput.internal_coordinate_tensors(scans, [(0, 1)], order=0))

    @validationTest
    def test_DihedConstrainedOpt(self):
        import McUtils.Coordinerds as coordops

        # methanol = Molecule.from_string(
        #     'methanol',
        #     # energy_evaluator='aimnet2'
        # )


        methanol = Molecule(
            ['C', 'O', 'H', 'H', 'H', 'H'],
            [[-0.6988896, 0.00487717, 0.00528378],
             [1.69694605, -1.08628154, -0.22505594],
             [-1.27384807, -0.22559494, 1.98568702],
             [-0.59371792, 2.01534286, -0.59617633],
             [-2.04665278, -0.99128091, -1.21504054],
             [2.91616232, 0.28293736, 0.04530201]],
            # energy_evaluator='aimnet2'
        )

        constraints = [(2, 0, 1, 5)]
        meth_opt = methanol.optimize(
            max_displacement=.5,
            max_iterations=10,
            coordinate_constraints=constraints,
            # logger=True
        )

        eng0 = methanol.calculate_energy() * UnitsData.hartrees_to_wavenumbers
        eng1 = meth_opt.calculate_energy() * UnitsData.hartrees_to_wavenumbers

        print(eng1 - eng0)

        ints0 = nput.internal_coordinate_tensors(methanol.coords, constraints, order=0)[0]
        ints1 = nput.internal_coordinate_tensors(meth_opt.coords, constraints, order=0)[0]
        print(ints0)
        print(ints1)

        methanol_zmatrix = coordops.reindex_zmatrix(
            coordops.functionalized_zmatrix(
                3,
                {
                    (2, 1, 0): [
                        [0, -1, -2, -3],
                        [1, -1, 0, -2],
                        [2, -1, 0, 1],
                    ]
                }
            ),
            [5, 1, 0, 2, 3, 4]
        )

        print(
            coordops.zmatrix_unit_convert(
                coordops.cartesian_to_zmatrix(methanol.coords, methanol_zmatrix).coords,
                UnitsData.convert("BohrRadius", "Angstroms"),
                rad2deg=True
            )
        )
        print(
            coordops.zmatrix_unit_convert(
                coordops.cartesian_to_zmatrix(meth_opt.coords, methanol_zmatrix).coords,
                UnitsData.convert("BohrRadius", "Angstroms"),
                rad2deg=True
            )
        )

    @validationTest
    def test_DihedOptRPNMScan(self):
        import McUtils.Coordinerds as coordops
        from Psience.Modes import ReactionPathModes

        methanol = Molecule.from_string(
            'methanol',
            energy_evaluator='aimnet2'
        ).optimize(max_iterations=250)

        methanol_zmatrix = coordops.reindex_zmatrix(
            coordops.functionalized_zmatrix(
                3,
                {
                    (2, 1, 0): [
                        [0, -1, -2, -3],
                        [1, -1, 0, -2],
                        [2, -1, 0, 1],
                    ]
                }
            ),
            [5, 1, 0, 2, 3, 4]
        )
        me_int = methanol.modify(
            internals=methanol_zmatrix,
        )

        scan_coord = [(2, 0, 1, 5)]

        scan_spec = [-np.deg2rad(180), np.deg2rad(0), 6]
        me_coords = me_int.get_scan_coordinates(
            [scan_spec],
            which=coordops.zmatrix_indices(methanol_zmatrix, scan_coord),
            internals='reembed',
            strip_embedding=True
        )

        opt_mols = [
            methanol.modify(coords=c).optimize(coordinate_constraints=scan_coord)
            for c in me_coords
        ]
        opt_coords = np.array([o.coords for o in opt_mols])

        me_aimnet2 = methanol.get_energy_function('aimnet2')
        # me_base_expansions = [me_aimnet2(c, order=2) for c in me_coords]
        me_expansions = [me_aimnet2(c, order=2) for c in opt_coords]

        me_grads = [exp[1] for exp in me_expansions]
        me_hesss = [exp[2] for exp in me_expansions]
        me_proj = nput.translation_rotation_projector(
            opt_coords,
            methanol.atomic_masses,
            mass_weighted=True
        )
        # me_proj_grads = me_proj @ np.array(me_grads)[:, :, np.newaxis]

        rpnms_opt = ReactionPathModes.get_rp_modes(
            me_grads,
            me_hesss,
            methanol.atomic_masses,
            projector=me_proj
        )

    @validationTest
    def test_AIMNetExpansions(self):
        methanol = Molecule(
            ['C', 'O', 'H', 'H', 'H', 'H'],
            [[-0.71174571, 0.0161939, 0.02050266],
             [1.71884591, -1.07310118, -0.2778059],
             [-1.30426891, 0.02589585, 1.99632677],
             [-0.77962613, 1.94036941, -0.7197672],
             [-2.02413643, -1.14525287, -1.05166036],
             [2.91548382, -0.08353621, 0.65084457]],
            energy_evaluator='aimnet2'
        )

        expansion = methanol.calculate_energy(order=3)

    @validationTest
    def test_Raman(self):
        water = Molecule.from_file(TestManager.test_data("water_freq.fchk"))
        print(water.get_harmonic_raman_spectrum())

    @validationTest
    def test_AIMNetDipoles(self):
        water = Molecule(
            ["O", "H", "H"],
            [[-3.09880964e-09, 1.23091261e-01, 0.00000000e+00],
             [-1.43810501e+00, -9.76773835e-01, 0.00000000e+00],
             [1.43810505e+00, -9.76773797e-01, 0.00000000e+00]],
            energy_evaluator='aimnet2',
            charge_evaluator='aimnet2',
            dipole_evaluator='aimnet2',
            polarizability_evaluator='aimnet2'
        )

        print()
        print(water.calculate_dipole())

        pol = water.calculate_dipole_polarizability(order=1)
        for p in pol:
            print([pp.shape for pp in p])

        print(pol[1][0])
        print(pol[1][1].shape)

        return

        water2 = water.modify(dipole_evaluator='expansion')
        water2.get_normal_modes()

        print(water.get_cartesian_dipole_derivatives(include_constant_term=True))
        print(water2.get_cartesian_dipole_derivatives(include_constant_term=True))

        return

        methanol = Molecule(
            ['C', 'O', 'H', 'H', 'H', 'H'],
            [[-0.71174571, 0.0161939, 0.02050266],
             [1.71884591, -1.07310118, -0.2778059],
             [-1.30426891, 0.02589585, 1.99632677],
             [-0.77962613, 1.94036941, -0.7197672],
             [-2.02413643, -1.14525287, -1.05166036],
             [2.91548382, -0.08353621, 0.65084457]],
            energy_evaluator='aimnet2',
            dipole_evaluator='aimnet2'
        )
        # expansion = methanol.calculate_energy(order=2)
        expansion = methanol.calculate_dipole(order=1)
        print([e.shape for e in expansion])

    @validationTest
    def test_RPNMVPT(self):
        methanol = Molecule.from_file(TestManager.test_data('methanol_vpt_3.fchk'))
        runner, _ = methanol.setup_VPT(use_reaction_path=True, degeneracy_specs='auto')
        runner.print_tables()

    @validationTest
    def test_LocalModeCHModel(self):
        from Psience.BasisReps import modify_internal_hamiltonian

        methanol = Molecule.from_file(TestManager.test_data('methanol_vpt_3.fchk'))
        base_nms = methanol.get_normal_modes()

        internals = {
                (0, 1):"OH",  # OH
                (1, 2):"CO",  # CO
                (2, 3):"CH3_stretch",  # CH1
                (2, 4):"CH3_stretch",  # CH2
                (2, 5):"CH3_stretch",  # CH3
                (1, 2, 3):"CH3_bend",  # OCH1
                (1, 2, 4):"CH3_bend",  # OCH2
                (1, 2, 5):"CH3_bend",  # OCH3
        }
        loc_modes = base_nms.localize(
            internals=internals
        )

        ob_modes = loc_modes.make_oblique()
        # f_ob = ob_modes.compute_hessian()
        # g_ob = ob_modes.compute_gmatrix()
        # h_ob = ob_modes.local_hessian
        # print(TableFormatter("{:.0f}").format(f_ob * UnitsData.hartrees_to_wavenumbers))
        # print(TableFormatter("{:.0f}").format(g_ob * UnitsData.hartrees_to_wavenumbers))
        # print(TableFormatter("{:.0f}").format(h_ob * UnitsData.hartrees_to_wavenumbers))
        # return
        base_hess = ob_modes.local_hessian

        new_hess = modify_internal_hamiltonian(
            base_hess,
            {
                (0, 1): "OH",
                (1, 2): "CO",
                (2, 3): "CH3_stretch",
                (2, 4): "CH3_stretch",
                (2, 5): "CH3_stretch",
                (1, 2, 3): "CH3_bend",
                (1, 2, 4): "CH3_bend",
                (1, 2, 5): "CH3_bend",
            },
            scaling_types={
                "OH":.93,
                "CH3_stretch":.96,
            },
            coupling_types={
                (1, "CH3_stretch", "CH3_stretch"):-22 * UnitsData.convert("Wavenumbers", "Hartrees")
            }
        )


        g_mat = loc_modes.local_gmatrix
        def print_arr(header, array=None):
            if array is None:
                array = header
                header = []
            elif isinstance(header, str):
                header = [header]
            if array.ndim == 1:
                array = array[np.newaxis]
            print(*header,  TableFormatter('{:.0f}').format(array * UnitsData.hartrees_to_wavenumbers) )

        print()
        print("="*25, "Unscaled", "="*25)
        freqs_old, _ = scipy.linalg.eigh(base_hess, g_mat, type=3)
        print_arr("Freqs:", np.sqrt(freqs_old))
        print("Hessian:")
        print_arr(base_hess)

        print("="*25, "Scaled", "="*25)
        freqs_new, _ = scipy.linalg.eigh(new_hess, g_mat, type=3)
        print_arr("Freqs:", np.sqrt(freqs_new))
        print("Hessian:")
        print_arr(new_hess)

    @validationTest
    def test_Caching(self):
        from McUtils.Scaffolding import Checkpointer
        import McUtils.Coordinerds as coordops
        import tempfile as tf

        methanol = Molecule(
            ['C', 'O', 'H', 'H', 'H', 'H'],
            [[-0.71174571,  0.0161939, 0.02050266],
             [ 1.71884591, -1.07310118, -0.2778059],
             [-1.30426891,  0.02589585, 1.99632677],
             [-0.77962613,  1.94036941, -0.7197672],
             [-2.02413643, -1.14525287, -1.05166036],
             [ 2.91548382, -0.08353621, 0.65084457]],
            energy_evaluator='aimnet2',
            # internals=methanol_zmatrix
        )

        methanol_zmatrix = coordops.reindex_zmatrix(
            coordops.functionalized_zmatrix(
                3,
                {
                    (2, 1, 0): [
                        [0, -1, -2, -3],
                        [1, -1, 0, -2],
                        [2, -1, 0, 1],
                    ]
                }
            ),
            [5, 1, 0, 2, 3, 4]
        )

        with tf.NamedTemporaryFile(suffix='.json') as temp:
            temp.close()

            cache = Checkpointer.from_file(temp.name)
            meth_int = methanol.modify(internals=methanol_zmatrix)
            with cache:
                cache['mol'] = meth_int

            cache = Checkpointer.from_file(temp.name)
            meth_int2 = cache['mol']

        # return

        with tf.NamedTemporaryFile(suffix='.json') as temp:
            temp.close()

            cache = Checkpointer.from_file(temp.name)
            modes = methanol.get_normal_modes(zero_freq_cutoff=50 / UnitsData.hartrees_to_wavenumbers)
            with cache:
                cache['modes'] = modes
                with prof.Timer("base"):
                    potential_derivatives = cache.cached_eval(
                        'potential_derivatives',
                        lambda: methanol.partial_force_field(
                            order=4,
                            modes=modes
                        )
                    )



            cache = Checkpointer.from_file(temp.name)
            modes2 = cache['modes']
            with prof.Timer("cached"):
                potential_derivatives = cache.cached_eval(
                    'potential_derivatives',
                    lambda: methanol.partial_force_field(
                        order=4,
                        modes=modes
                    )
                )
            print(
                np.asanyarray(potential_derivatives[-1]).shape
            )
        # methanol.animate_mode(0, modes=modes2).show()
        # print(modes2 is modes)

    @validationTest
    def test_OrcaImport(self):
        # propyl = Molecule.from_file(
        #     TestManager.test_data('proplybenz.out'),
        #     'orca'
        # )
        #
        # print(
        #     propyl.get_normal_modes(
        #         # project_transrot=False
        #     ).freqs * UnitsData.hartrees_to_wavenumbers
        # )

        propyl = Molecule.from_file(
            TestManager.test_data('proplybenz.hess')
        )

        print(
            propyl.get_normal_modes(
                # project_transrot=False
            ).freqs * UnitsData.hartrees_to_wavenumbers
        )

    @validationTest
    def test_PartialQuartic(self):
        import McUtils.Coordinerds as coordops

        methanol_zmatrix = coordops.reindex_zmatrix(
            coordops.functionalized_zmatrix(
                3,
                {
                    (2, 1, 0): [
                        [0, -1, -2, -3],
                        [1, -1, 0, -2],
                        [2, -1, 0, 1],
                    ]
                }
            ),
            [5, 1, 0, 2, 3, 4]
        )

        crds_aimnet_opt = np.array(
            [[-0.71174571,  0.0161939, 0.02050266],
             [ 1.71884591, -1.07310118, -0.2778059],
             [-1.30426891,  0.02589585, 1.99632677],
             [-0.77962613,  1.94036941, -0.7197672],
             [-2.02413643, -1.14525287, -1.05166036],
             [ 2.91548382, -0.08353621, 0.65084457]]
        )
        crds_bad = np.array(
            [[[[0.99603661, -0.03076131, 0.317282],
               [2.28225031, -0.60719144, 0.1594239],
               [0.68719378, -0.0280038, 1.36425206],
               [0.95774141, 0.98826296, -0.07215449],
               [0.29937143, -0.64460867, -0.24823604],
               [2.91548382, -0.08353621, 0.65084457]]],
             [[[0.99603661, -0.03076131, 0.317282],
               [2.28225031, -0.60719144, 0.1594239],
               [0.68248684, -0.02562726, 1.36284309],
               [0.96011584, 0.98746852, -0.07445194],
               [0.30154935, -0.64537247, -0.25008224],
               [2.91548382, -0.08353621, 0.65084457]]],
             [[[0.99603661, -0.03076131, 0.317282],
               [2.28225031, -0.60719144, 0.1594239],
               [0.67778774, -0.02325085, 1.36140798],
               [0.9625001, 0.98666785, -0.07673702],
               [0.3037351, -0.64614144, -0.25191702],
               [2.91548382, -0.08353621, 0.65084457]]]]
        ) / UnitsData.bohr_to_angstroms

        methanol = Molecule(
            ['C', 'O', 'H', 'H', 'H', 'H'],
            crds_aimnet_opt,
            energy_evaluator='rdkit',
            # internals=methanol_zmatrix
        ).optimize()

        woof = methanol.partial_force_field(
            order=3,
            modes=methanol.get_normal_modes(project_transrot=False)
        )

        print(woof[-1][0, 0])
        print(woof[-1][0, 1])

        return

        # der = methanol.calculate_energy(order=1, coords=crds_bad)
        # print(der[1])
        # return

        # print(methanol.calculate_energy() - methanol.optimize().calculate_energy())
        # return

        # der = methanol.calculate_energy(order=1)
        # print(der[1])
        # return

        proj_dir = os.path.expanduser("~/Documents/Postdoc/Projects/CoordinatePaper/ml_fd_tests/")
        os.makedirs(proj_dir, exist_ok=True)
        from McUtils.Formatters import TableFormatter

        modes = methanol.get_normal_modes()#zero_freq_cutoff=50 / UnitsData.hartrees_to_wavenumbers)

        # from Psience.Psience.Molecools.Evaluator import AIMNet2EnergyEvaluator
        # AIMNet2EnergyEvaluator.analytic_derivative_order = 3
        # expansion = methanol.calculate_energy(order=3, mesh_spacing=0.01, stencil=3)
        # terms = np.tensordot(modes.coords_by_modes, expansion[-1], axes=[-1, 0])
        #
        # with open(os.path.join(proj_dir, f"aimnet_analytic.txt"), 'w+') as outs:
        #     for term in terms:
        #         print("=" * 100, file=outs)
        #         print(TableFormatter("{:.3f}").format(
        #             term * UnitsData.hartrees_to_wavenumbers
        #         ), file=outs)
        #
        #
        # AIMNet2EnergyEvaluator.analytic_derivative_order = 2
        #
        # return



        # for step_size in [0.1, 0.25, 0.5, 1]:
        #     for stencil in [3, 5, 7]:
        for step_size in [1]:
            for stencil in [5]:
                # with open(os.path.join(proj_dir, f"aimnet_fd_{step_size*100:.0f}_{stencil}.txt"), 'w+') as outs:
                print(f"Mesh spacing: {step_size}")
                print(f"Stencil: {stencil}")
                print(f"Freqs: {modes.freqs * UnitsData.hartrees_to_wavenumbers}")
                expansion = methanol.partial_force_field(order=3, modes=modes,
                                                         mesh_spacing=step_size,
                                                         stencil=stencil)
                terms = expansion[-1]
                # term = (expansion[-1])[0]
                # term = expansion[-1][0].reshape(-1, 18)
                # print([np.asanyarray(e).shape for e in expansion])
                for term in terms:
                    print("="*100)
                    print(TableFormatter("{:.3f}").format(
                        term * UnitsData.hartrees_to_wavenumbers
                    ))
                # print(expansion[-1][0, 0])

    @validationTest
    def test_InternalProjectedModes(self):
        import McUtils.Coordinerds as crops

        methanol_zmatrix = crops.functionalized_zmatrix(
            3,
            {
                (2, 1, 0): [
                    [0, -1, -2, -3],
                    [1, -1, 0, -2],
                    [2, -1, 0, 1],
                ]
            }
        )
        methanol_zmatrix = crops.set_zmatrix_embedding(methanol_zmatrix)

        me_ints = Molecule.from_file(
            TestManager.test_data('methanol_vpt_1.fchk'),
            internals=methanol_zmatrix
        )
        nms = me_ints.get_normal_modes(project_transrot=False)
        locs = nms.localize(coordinate_constraints=crops.zmatrix_indices(
            methanol_zmatrix,
            [(3, 2, 1, 0)]
        ))
        print()

        from McUtils.Formatters import TableFormatter
        print(TableFormatter('{:.0f}').format(nms.freqs[np.newaxis] * UnitsData.hartrees_to_wavenumbers))
        print(
            TableFormatter('{:.0f}').format(locs.local_hessian * UnitsData.hartrees_to_wavenumbers)
        )

        me_carts = Molecule.from_file(
            TestManager.test_data('methanol_vpt_1.fchk')
        )
        nms_carts = me_carts.get_normal_modes(project_transrot=False, use_internals=False)

        # locs = nms.localize(internals=[(3, 2, 1, 0)], orthogonal_projection=True)
        # print(TableFormatter('{:.0f}').format(nms.freqs[np.newaxis] * UnitsData.hartrees_to_wavenumbers))
        # print(
        #     TableFormatter('{:.0f}').format(locs.local_hessian * UnitsData.hartrees_to_wavenumbers)
        # )
        # loc_2 = nms_carts.apply_transformation(locs.localizing_transformation[0])
        # print(TableFormatter('{:.0f}').format(nms_carts.freqs[np.newaxis] * UnitsData.hartrees_to_wavenumbers))
        # print(
        #     TableFormatter('{:.0f}').format(loc_2.local_hessian * UnitsData.hartrees_to_wavenumbers)
        # )
        # print(
        #     TableFormatter('{:.0f}').format(loc_2.compute_gmatrix())
        # )
        # return

        loc_2 = nms_carts.apply_transformation(locs.localizing_transformation).make_dimensionless()
        # cart_udim = nms_carts.make_dimensionless()
        f_nmw = me_carts.potential_derivatives[1]
        g12 = nput.fractional_power(me_carts.g_matrix, 1/2)
        f_mw = g12 @ f_nmw @ g12
        f_cart = nput.tensor_reexpand(
            [loc_2.coords_by_modes],
            [0, f_mw]
        )[-1]
        # g_cart = nms_carts.compute_gmatrix()
        # print(
        #     TableFormatter('{:.3f}').format(locs.localizing_transformation[1] @ locs.localizing_transformation[0])
        # )
        # f_loc = locs.localizing_transformation[1] @ f_cart @ locs.localizing_transformation[1].T
        # print(f_loc.shape)
        print(
            TableFormatter('{:.3f}').format(
                f_cart * (UnitsData.hartrees_to_wavenumbers)
            )
        )
        # print(locs.localizing_transformation[1] @ locs.localizing_transformation[0])
        return

        runner, _ = me_ints.setup_VPT(states=2,
                                       degeneracy_specs='auto',
                                       cartesian_analytic_deriv_order=-1,
                                       internal_by_cartesian_derivative_method='fast',
                                       modes=nms_carts,
                                       mode_transformation=locs.localizing_transformation
                                       )
        runner.print_tables()






    @validationTest
    def test_MultiGMatrix(self):


        mol = Molecule.from_file(
            TestManager.test_data("nh3.fchk"),
            internals=[
                [0, -1, -1, -1],
                [1,  0, -1, -1],
                [2,  0,  1, -1],
                [3,  0,  1,  2]
            ]
        )
        other_structs = mol.get_displaced_coordinates(
            np.random.rand(5, 4, 3) / 5
        )
        gms = mol.get_gmatrix(coords=other_structs)

    @validationTest
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

    @validationTest
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

    @validationTest
    def test_FormatImports(self):
        woof = Molecule.from_file(
            TestManager.test_data('water_dimer_freq_unopt.log'),
            'gspec'
        )
        print(len(woof.potential_derivatives))
        # woof.plot().show()

    @validationTest
    def test_ModeSelectedNMs(self):
        propylbenzene = Molecule.from_file(
            TestManager.test_data('proplybenz.hess')
        )
        modes = propylbenzene.get_normal_modes()
        loc_1 = modes.localize(
            internals=[(19, 8)]
        )
        stretches = modes[tuple(t-7 for t in [57, 58, 59, 60, 61, 62, 63])]
        loc_2 = stretches.localize(
            internals=[(19, 8)]
        )
        print(loc_1.local_freqs * UnitsData.hartrees_to_wavenumbers)
        print(loc_2.local_freqs * UnitsData.hartrees_to_wavenumbers)

    @inactiveTest
    def test_NMFiniteDifference(self):

        jobs_dir = ...
        def check_directory():
            ...
        def create_jobs(displacements):
            ...
        def prep_jobs_directory(spec, displacements, jobs):
            ...
        propylbenzene_setup = Molecule.from_file(
            TestManager.test_data('proplybenz.hess'),
            energy_evaluator=create_jobs
        )
        setup = propylbenzene_setup.partial_force_field()

        # run jobs
        def parse_jobs():
            ...
        propylbenzene_parse = Molecule.from_file(
            TestManager.test_data('proplybenz.hess'),
            energy_evaluator=create_jobs
        )
        setup = propylbenzene_setup.partial_force_field()

    @validationTest
    def test_CoordinateSystems(self):
        import McUtils.Coordinerds as coordops

        propylbenzene = Molecule.from_file(
            TestManager.test_data('proplybenz.hess')
        )
        coords = propylbenzene.get_bond_graph_internals(pruning=True)
        # g12 = propylbenzene.get_gmatrix(power=1/2)
        # def b_gen(pos, crds):
        #     return g12 @ nput.internal_coordinate_tensors(propylbenzene.coords, crds, order=1)[1]
        #
        # pruned = coordops.prune_internal_coordinates(
        #     coords,
        #     b_gen,
        #     method='b_matrix'
        # )
        # # pruned = coordops.prune_internal_coordinates(coords)
        # print(coords)

    @validationTest
    def test_PySCFEnergy(self):
        # ethanol = Molecule.from_string('CCO', energy_evaluator='xtb')
        # print(ethanol.calculate_energy())

        ethanol = Molecule.from_string('CCO',
            energy_evaluator=(
                'pyscf',
                {'level_of_theory': 'b3lyp', 'basis_set': 'ccpvdz'}
            )
        )
        print(ethanol.calculate_energy())

    @validationTest
    def test_MACEEnergy(self):
        ethanol = Molecule.from_string('CCO', energy_evaluator='mace')
        print(ethanol.calculate_energy())

        # ethanol = ethanol.modify(
        #     energy_evaluator=(
        #         'pyscf',
        #         {'level_of_theory': 'b3lyp', 'basis_set': 'ccpvdz'}
        #     )
        # )
        # print(ethanol.calculate_energy())

    @validationTest
    def test_UMAEnergy(self):
        ethanol = Molecule.from_string('CCO', energy_evaluator='uma')
        print(ethanol.calculate_energy())

    @validationTest
    def test_BackboneChains(self):
        from Psience.Molecools import Molecule
        import McUtils.Coordinerds as coordops
        from Psience.Reactions import Reaction

        # woof = Reaction.from_smiles("C=C.C=CC=C>>C1CCC=CC1",
        #                             fragment_expansion_method='centroid',
        #                             optimize=True,
        #                             min_distance=.1,
        #                             add_radius=False,
        #                             expansion_factor=.01,
        #                             )

        woof = Reaction.from_smiles("C=C.C=CC=C(c1c2ccccc2ccc1)>>C1CCC=CC1(c1c2ccccc2ccc1)",
                                    fragment_expansion_method='centroid',
                                    optimize=True,
                                    min_distance=.1,
                                    add_radius=False,
                                    expansion_factor=.01,
                                    )

        reactant_complex = woof.reactant_complex
        full_zmat = reactant_complex.get_bond_zmatrix()
        int_comp = reactant_complex.modify(internals=full_zmat)

        # int_comp.animate_coordinate(30-6).show()

        return

        woof = Molecule.construct('CCCC')
        zm = coordops.chain_zmatrix(4)

        print(woof.atoms)
        print(
            coordops.add_missing_zmatrix_bonds(
                zm,
                [b[:2] for b in woof.bonds]
            )
        )




        return


        napthalene = Molecule.construct('CCCCC(c1c2ccccc2ccc1)CCCC')
        # backbone = napthalene.find_heavy_atom_backbone()



        chains = napthalene.edge_graph.segment_by_chains()
        zm = coordops.bond_graph_zmatrix(
            [b[:2] for b in napthalene.bonds],
            chains
        )

        print(zm)
        return

        backbone, (side_chain,) = napthalene.edge_graph.segment_by_chains()
        atom_styles = {
            backbone[0]: {'color': 'white', 'glow': 'red'},
            backbone[-1]: {'color': 'white', 'glow': 'blue'}
        }
        for a in side_chain:
            atom_styles[a] = {'color': 'white', 'glow': 'purple'}
        bond_style = {
            k: {'color': 'white', 'glow': 'red'}
            for i in range(len(backbone) - 1)
            for k in [
                (backbone[i], backbone[i + 1]),
                (backbone[i + 1], backbone[i])
            ]
        }

        for i in range(len(side_chain) - 1):
            bond_style[(side_chain[i], side_chain[i + 1])] = {'color': 'white', 'glow': 'purple'}
            bond_style[(side_chain[i + 1], side_chain[i])] = {'color': 'white', 'glow': 'purple'}
        napthalene.plot(
            highlight_atoms=backbone[1:-1],
            atom_style=atom_styles,
            bond_style=bond_style,
            include_save_buttons=True
        ).show()

    @validationTest
    def test_EasyZMatrices(self):
        # cpmo = Molecule.from_file(
        #     TestManager.test_data('cpmo3m_opt.xyz'),
        #     units='Angstroms'
        # )
        #
        # cpmo_split = cpmo.modify(
        #     bonds=[
        #         b for b in cpmo.bonds
        #         if tuple(sorted(b[:2])) not in {
        #             (0, 4),
        #             (0, 5),
        #             (0, 6),
        #             (0, 7),
        #             (0, 8)
        #         }
        #     ]
        # )
        #
        # pprint.pprint(
        #     cpmo_split.get_bond_zmatrix(
        #         attachment_points={0:(4, 6, 8)}
        #     )
        # )

        import McUtils.Coordinerds as coordops

        alt_mol = Molecule.from_string("O=C(COc1cc(Cl)ccc1Cl)N/N=C/c1c[nH]c2c1cccc2", "smi")
        zmat = alt_mol.get_bond_zmatrix()
        coordops.validate_zmatrix(zmat, raise_exception=True)
        pprint.pprint(zmat)

    @validationTest
    def test_PointGroups(self):
        print()
        # from Psience.Symmetry import (
        #     PointGroupIdentifier, PointGroup,
        #     SymmetryElement,
        #     InversionElement, RotationElement, ReflectionElement, ImproperRotationElement
        # )

        # water = Molecule.from_file(TestManager.test_data('water_freq.fchk'))
        # mol = Molecule.construct("CC", energy_evaluator='rdkit').optimize()
        # print(mol.atoms)
        # print(mol.coords.tolist())
        mol = Molecule(
            ['C', 'C', 'H', 'H', 'H', 'H', 'H', 'H'],
            [
                [1.3054152479869834, 0.24647130925656535, -0.5530793729770965],
                [-1.3054057993563604, -0.24642973528182355, 0.553026460645607],
                [2.7502129154741146, -0.7828587313748425, 0.5769125988608768],
                [1.3606991857629105, -0.4278075384495296, -2.5445918158535847],
                [1.7132634991082158, 2.309094146202721, -0.4999931866841094],
                [-2.751762490896308, 0.7821255176384876, -0.5756408131790036],
                [-1.3610034316689752, 0.4288638953531954, 2.544251665151152],
                [-1.7114191264105811, -2.3094531941664, 0.49911446403615845]
            ]
        )
        mol = Molecule(
            ['C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H'],
            [[ 1.5314834 , -2.07901109,  0.03362855],
             [ 2.66130207,  0.27021048, -0.00918031],
             [ 1.09977014,  2.36405638, -0.04301921],
             [-1.49251245,  2.07781621, -0.03367612],
             [-2.66353391, -0.24015647,  0.00865698],
             [-1.08442635, -2.34605922,  0.04267503],
             [ 2.65427297, -3.82363066,  0.06211177],
             [ 4.65827828,  0.4512349 , -0.01567418],
             [ 1.93012887,  4.24957282, -0.07725653],
             [-2.65023509,  3.80315125, -0.06177008],
             [-4.68709859, -0.49931772,  0.01657436],
             [-1.95742935, -4.22786689,  0.07692975]]
        ).get_embedded_molecule(load_properties=False)
        (coords, atoms), pg = mol.symmetrize(grouping_tol=.8, tol=.3, return_point_group=True)
        self.assertEquals(coords.shape, mol.coords.shape)
        base = mol.plot(backend='x3d', include_save_buttons=True)#, principle_axes=True)
        (coords, atoms), pg2 = Molecule(atoms, coords).symmetrize(pg, return_point_group=True)
        # print(pg2)
        # print(coords, atoms)
        Molecule(atoms, coords).plot(figure=base, highlight_atoms=True)
        # print(nput.distance_matrix(mol.coords[:6]))
        # print(atoms)
        # print(nput.distance_matrix(coords[-6:]))
        pg2.plot(figure=base, origin=mol.center_of_mass * UnitsData.bohr_to_angstroms)
        # base.show()
        # print(coords)
        return

        mol = Molecule(
            ["C", "H", "H", "H", "H"],
            [[ 1.49273393e-04, -8.31939792e-05, -3.01752454e-05],
             [-1.29672979e-01,  2.02724390e+00,  3.64678828e-01],
             [-9.90194225e-01, -4.48290077e-01, -1.75452490e+00],
             [-8.63649504e-01, -1.03676833e+00,  1.56171780e+00],
             [ 1.98411318e+00, -5.42519370e-01, -1.71993834e-01]]
        )
        # print(mol.coords)
        pg = mol.get_point_group(grouping_tol=.8, tol=.3, mom_tol=5, verbose=False)
        print(pg)
        print(pg.elements)
        base = mol.plot(backend='x3d', principle_axes=True)
        pg.plot(figure=base, origin=mol.center_of_mass * UnitsData.bohr_to_angstroms)
        base.show()

        # print(id.coord_data.rotor_type, id.coord_data.planar)
        # symm, pg = id.identify_point_group()
        # print(PointGroup.from_symmetry_elements(symm))
        # print(pg.get_character_table().format())

    @validationTest
    def test_Symmetrization(self):
        cpmo3m = Molecule.from_file(
            TestManager.test_data('cpmo3m_opt.xyz'),
            units='Angstroms'
        )
        new_struct, new_coords, pg = cpmo3m.symmetrize(sel=[0, 1, 2, 3, 4, 10], tol=.8, return_point_group=True)
        self.assertIsNot(new_struct, None)
        print(
            nput.distance_matrix(cpmo3m.coords[(0, 1, 2, 3, 4, 10,),]) -
            nput.distance_matrix(new_struct.coords[(0, 1, 2, 3, 4, 10,),])
        )
        # base_pol = new_struct.plot(backend='x3d')
        # cpmo3m.plot(highlight_atoms='all', figure=base_pol)
        # base_pol.show()

    @validationTest
    def test_PartialSymmetries(self):
        import McUtils.Coordinerds as coordops

        nh3 = Molecule.from_file(
            TestManager.test_data('nh3.fchk'),
            internals=coordops.spoke_zmatrix(3)
        ).get_embedded_molecule()

        # base = nh3.plot(backend='x3d', include_save_buttons=True)#, principle_axes=True)
        # pg = nh3.get_point_group()
        # pg.plot(figure=base,
        #         origin=nh3.center_of_mass * UnitsData.bohr_to_angstroms,
        #         elements='all'
        #         )
        # base.show()
        #
        # return

        coeffs, basis, expansions = nh3.get_symmetrized_internals(return_expansions=True)
        # print(coeffs[0])
        # print(basis)

        # return

        cpmo3m = Molecule.from_file(
            TestManager.test_data('cpmo3m_opt.xyz'),
            units='Angstroms'
        )
        zmat = coordops.reindex_zmatrix(
            coordops.functionalized_zmatrix(
                1,
                {
                    (0, -1, -2): coordops.chain_zmatrix(2),
                    (0, 1, 2): coordops.chain_zmatrix(2),
                    (0, 3, 4): coordops.chain_zmatrix(2),
                    (0, 5, 6): coordops.chain_zmatrix(2),
                    (0, 7, 8): coordops.chain_zmatrix(2),
                    (0, 9, 10): coordops.chain_zmatrix(2),
                    (0, 11, 12): coordops.chain_zmatrix(2),
                    (0, 13, 14): coordops.chain_zmatrix(2)
                }
            ),
            [10, 11, 12, 13, 14, 15, 16, 0, 9, 1, 8, 2, 7, 3, 6, 4, 5]
        )
        cpmo3m_int = cpmo3m.modify(
            # cpmo3m.atoms,  # [10:],
            # cpmo3m.coords,  # [10:],
            internals=zmat
        )

        coeffs, intenrals, expansions = cpmo3m_int.get_symmetrized_internals(
            atom_selection=list(range(11)),
            tol=.8,
            permutation_tol=.9,
            return_expansions=True,
            return_point_group=False,
            extra_internals=[
                (0, 1) # add in one CC bond length to get symmetrized
            ]
        )

        a1_modes, a1_inv = expansions[0]
        cpmo3m_int.animate_coordinate(0, coordinate_expansion=[
            a1_inv[0][(15,),]
        ]).show()




    @validationTest
    def test_SufaceArea(self):
        # propylbenzene = Molecule.from_file(
        #     TestManager.test_data('proplybenz.hess')
        # )
        # propylbenzene = Molecule.from_file(
        #     TestManager.test_data('methanol_vpt_1.fchk')
        # )
        propylbenzene = Molecule.from_file(
            TestManager.test_data('tbhp_180.fchk')
        )
        # propylbenzene = Molecule(
        #     np.asanyarray(propylbenzene.atoms)[(0, 6, 14),],
        #     propylbenzene.coords[(0, 6, 14),]
        # )

        scale = 1.0
        surf = propylbenzene.get_surface(samples=100, tolerance=1e-6, expansion=1e-12, radius_scaling=scale)
        print()
        # print(surf.surface_area(include_doubles=False, include_triples=False, include_quadruples=False))
        # print(surf.surface_area(include_triples=False, include_quadruples=False))
        # print(surf.surface_area(include_triples=True, include_quadruples=False))
        # print("Trips:", surf.surface_area(include_triples=True, include_quadruples=False))
        s1 = surf.surface_area(include_triples=True, include_quadruples=True) #* UnitsData.convert("BohrRadius", "Angstroms")**2
        s2 = surf.surface_area(method='sampling') #* UnitsData.convert("BohrRadius", "Angstroms")**2
        print("Quads:", s1)
        print("Sampled:", s2)
        print("Diff:", s1 - s2)
        # print(surf.surface_area(method='mesh'))
        # print(surf.surface_area(include_triples=True, include_quadruples=True))
        mol_plot = propylbenzene.plot(backend='x3d',
                                      background=['white', 'blue'],
                                      # highlight_atoms=[0, 2],
                                      atom_style={'transparency':.6},
                                      atom_radius_scaling=scale,
                                      capped_bonds=True,
                                      include_save_buttons=True,
                                      image_size=800
                                      )
        surf.plot(figure=mol_plot,
                  # color=None,
                  sphere_color=None,
                  plot_intersections=True
                  )
        # mesh = propylbenzene.get_surface_mesh(mesh_options={'depth':15})
        # mesh.plot(figure=mol_plot)

        mol_plot.show()

    @validationTest
    def test_SufaceTriangulation(self):
        from McUtils.Plots import ColorPalette

        # g = np.linspace(0, 6.28, 100)
        # fig = None
        # palette = ColorPalette('viridis').lighten(-.2)
        # for i in range(0, 6, 2):
        #     fig = plt.Plot(
        #         g, np.sin(i*g),
        #         color=palette(i/3, return_color_code=True),
        #         figure=fig
        #     )
        # fig.show()
        # return

        propylbenzene = Molecule.from_file(
            TestManager.test_data('proplybenz.hess')
        )
        mesh = propylbenzene.get_surface_mesh()
        mol_plot = propylbenzene.plot(backend='x3d', include_save_buttons=True)
        mesh.plot(
            # line_color=None,
            function=lambda pts:pts[:, 2],
            transparency=.2,
            figure=mol_plot
        )

        mol_plot.show()
        return

        # surf = SphereUnionSurface.from_xyz(
        #     propylbenzene.atoms,
        #     propylbenzene.coords,
        #     expansion=.01,
        #     samples=100
        # )

        # print([s.shape for s in surf.generate_points(preserve_origins=True)])
        #
        # return

        # pts = surf.sampling_points
        # verts, tris = surf.generate_mesh(pts)
        # mol_plot = propylbenzene.plot(backend='x3d', include_save_buttons=True)
        # plt.Line(pts[:15], glow='red', line_style='2', line_thickness=2).plot(mol_plot)
        # print(
        #     mol_plot.figure.to_x3d().to_x3d().tostring()
        # )
        # mol_plot.show()
        # return
        # delaunay = scipy.spatial.Delaunay(pts)
        # subtris = delaunay.points[delaunay.simplices][:, :, :3]
        # tri_points = verts[tris] * UnitsData.convert("BohrRadius", "Angstroms")
        # plt.Point(verts * UnitsData.convert("BohrRadius", "Angstroms"), color='red').plot(mol_plot)
        # plt.Triangle(tri_points, transparency=.8, color='black').plot(mol_plot)
        # plt.Line(tri_points, color='red').plot(mol_plot)
        # plt.Sphere(pts * UnitsData.convert("BohrRadius", "Angstroms"), .1, color='purple').plot(mol_plot)
        # print(
        #     mol_plot.figure.to_x3d().to_x3d().tostring()
        # )

        # mol_plot.show()
        return

        pts = surf.sampling_points
        dm = nput.distance_matrix(pts)
        np.fill_diagonal(dm, 100)
        print(np.min(dm))
        print(np.max(np.min(dm, axis=1)))

        pts2 = SphereUnionSurface.adjust_point_cloud_density(pts,
                                                             centers=surf.centers,
                                                             radii=surf.radii,
                                                             min_component=.6,
                                                             max_iterations=250
                                                             )
        dm = nput.distance_matrix(pts2)
        np.fill_diagonal(dm, 100)
        print(np.min(dm))
        print(np.max(np.min(dm, axis=1)))

        pts3 = SphereUnionSurface.point_cloud_repulsion(pts,
                                                        surf.centers,
                                                        surf.radii,
                                                        max_iterations=25
                                                        )
        dm = nput.distance_matrix(pts3)
        np.fill_diagonal(dm, 100)
        print(np.min(dm))
        print(np.max(np.min(dm, axis=1)))

        # pts4 = SphereUnionSurface.sphere_boundary_pruning(pts,
        #                                                   surf.centers
        #                                                   )
        # dm = nput.distance_matrix(pts4)
        # np.fill_diagonal(dm, 100)
        # print(np.min(dm))
        # print(np.max(np.min(dm, axis=1)))


        mol_plot = propylbenzene.plot(backend='x3d', image_size=[950, 700], include_save_buttons=True)
        plt.Sphere(pts * UnitsData.convert("BohrRadius", "Angstroms"), .1, color='teal').plot(mol_plot)
        mol_plot.show()

        mol_plot = propylbenzene.plot(backend='x3d', image_size=[950, 700], include_save_buttons=True)
        plt.Sphere(pts2 * UnitsData.convert("BohrRadius", "Angstroms"), .1, color='purple').plot(mol_plot)
        mol_plot.show()

        mol_plot = propylbenzene.plot(backend='x3d', image_size=[950, 700], include_save_buttons=True)
        plt.Sphere(pts3 * UnitsData.convert("BohrRadius", "Angstroms"), .1, color='red').plot(mol_plot)
        mol_plot.show()

        # pts = surf.sampling_points
        # print(np.min(pts))

    @validationTest
    def test_ModeLabels(self):
        import McUtils.Coordinerds as coordops

        propylbenzene = Molecule.from_file(
            TestManager.test_data('proplybenz.hess')
        )

        # for c,l in propylbenzene.get_labeled_internals().items():
        #     print(c, l)
        # return

        modes = propylbenzene.get_normal_modes()
        mode_labs = propylbenzene.get_mode_labels(pruning=False, use_redundants=True)

        for i,(freq,lab) in enumerate(zip(reversed(modes.freqs), reversed(mode_labs))):
            print(
                "Mode {}: {:.0f} {}".format(i+1, freq * UnitsData.hartrees_to_wavenumbers,
                                            "mixed"
                                                if lab.type is None else
                                            lab.type
                                            )
            )

        return

    @validationTest
    def test_HamiltonianExpansions(self):
        from Psience.BasisReps import TaborCHModel

        print()
        tbhp = Molecule.from_file(
            TestManager.test_data('tbhp_180.fchk')
        )

        print(
            mfmt.format_mode_labels(
                tbhp.get_mode_labels(),
                tbhp.get_normal_modes().freqs * UnitsData.hartrees_to_wavenumbers
            )
        )
        return

        print(
            mfmt.format_symmetric_tensor_elements(
                tbhp.potential_derivatives[1] * UnitsData.hartrees_to_wavenumbers,
                cutoff=1000
            )
        )

        return

        model = TaborCHModel.from_molecule(
            tbhp,
            oblique=True
        )

        ham = tbhp.get_hamiltonian(modes=model.modes)

        print(
            mfmt.TableFormatter("{:.0f}").format(model.f * UnitsData.hartrees_to_wavenumbers)
        )

        v_exp = ham.potential_expansion(2)
        print(len(v_exp))

    @validationTest
    def test_AutoCHModel(self):
        # import McUtils.Devutils as dev
        #
        # print(dev.merge_dicts(
        #     {},
        #     {'a':{"1":{1}}}
        # ))
        # return

        import McUtils.Coordinerds as coordops
        from Psience.BasisReps import LocalHarmonicModel, StateMaker, TaborCHModel

        propylbenzene = Molecule.from_file(
            TestManager.test_data('proplybenz.hess')
        )

        # print(
        #     TaborCHModel.from_molecule(propylbenzene).internals
        # )
        # return


        lhm = LocalHarmonicModel
        model = LocalHarmonicModel.from_molecule(
            propylbenzene,
            oblique=False,
            # coordinate_filter=lambda coords: {
            #     c: l
            #     for c, l in coords.items()
            #     if l.atoms in {'CH', 'HCH'}
            # },
            # localization_mode_spaces = {
            #     ("CH", "stretch"):[-12, -11, -10, -9, -8, -7, -6],
            #     ("HCH", "bend"):[-21, -19, -18, -17, -16]
            # },
            # localization_mode_spaces={
            #     ("CH", "stretch"): (("methyl", "ethyl"), "CH", "stretch"),
            #     ("HCH", "bend"): [
            #         [1390 / UnitsData.hartrees_to_wavenumbers, 1505 / UnitsData.hartrees_to_wavenumbers],
            #         ("bend",)
            #     ]
            # },
            # mode_labels=True,
            coordinate_filter=lambda coords: {
                c:l
                # c: (l,
                #     [-12, -11, -10, -9, -8, -7, -6]
                #         if l.atoms == "CH" else
                #     [-21, -19, -18, -17, -16]
                #     )
                for c, l in coords.items()
                if l.atoms in {'CH', 'HCH'} #and l.ring != 'benzene'
            },
            # anharmonic_scalings={
            #     # lhm.state("CH", "stretch"): .96,  # diagonal scaling
            #     # lhm.state_pair(
            #     #     ("methyl", "CH", "stretch"),
            #     #     ("methyl", "CH", "stretch")
            #     # ): .9,  # methyl-methyl stretch scaling
            #     # lhm.state("HCH", "bend"): .96,
            #     # lhm.state_pair(
            #     #     ("ethyl", "CH", "stretch"),
            #     #     ("ethyl", "CH", "stretch")
            #     # ): .96,  # ethyl-ethyl stretch scaling
            #     lhm.state(
            #         ("HCH", "bend"),
            #         ("HCH", "bend")
            #     ): 0.975,
            # },
            anharmonic_couplings={
                # lhm.state_pair(
                #     2, # shared atoms
                #     ("CH", "stretch"), # stretch fundamental
                #     (2, "HCH", "bend") # bend overtone
                #  ): 22 / UnitsData.hartrees_to_wavenumbers,
                lhm.state_pair(
                    ((2, 1),),  # shared atoms
                    ("CH", "stretch"),  # stretch fundamental
                    (("HCH", "bend"), ("HCH", "bend"))  # bend overtone
                ): 5.6 / UnitsData.hartrees_to_wavenumbers
            },
            anharmonic_shifts = {
                lhm.state("benzene", None, "CH", "stretch"): -30 / UnitsData.hartrees_to_wavenumbers,
                lhm.state("methyl", "CH", "stretch"): -8 / UnitsData.hartrees_to_wavenumbers,
                lhm.state("ethyl", "CH", "stretch"): -5 / UnitsData.hartrees_to_wavenumbers,
                lhm.state("HCH", "bend"): -2*9.6 / UnitsData.hartrees_to_wavenumbers
            }
        )

        # import pprint
        # pprint.pprint(model.scalings)
        # return


        dim = model.basis.ndim
        state = StateMaker(dim, mode='high-low')
        wfns = model.get_wavefunctions({
            'max_freq': 3250 / UnitsData.hartrees_to_wavenumbers,
            'min_freq': 3050 / UnitsData.hartrees_to_wavenumbers,
            'max_quanta': 3
        })
        # print(wfns.basis.excitations)
        #     state(1),
        #     state(2),
        #     state(3),
        #     state(4),
        #     state(5),
        #     state(6), # methyl CH stretch
        #     state(7), # methyl CH stretch
        #     state([dim-4, 2]), # methyl HCH bend
        #     state([dim - 3, 2]),  # methyl HCH bend
        #     state(dim - 4, dim - 3),  # methyl combination
        # ])

        print(model.internals)
        print(
            mfmt.TableFormatter("{:.3f}").format(wfns.hamiltonian * UnitsData.hartrees_to_wavenumbers)
        )

        spec = wfns.get_spectrum()
        # spec.plot().show()
        print(spec.intensities)
        return

        spec.plot().show()

        ham = model.get_hamiltonian(
            [
                # state(1), # benzene CH stretch
                # state(6), # methyl CH stretch
                # state(7), # methyl CH stretch
                # state([dim-4, 2]), # methyl HCH bend
                state([dim-3, 2]), # methyl HCH bend
                state(dim-4, dim-3), # methyl combination
            ]
        )
        print()
        print(
            TableFormatter("{:.0f}").format(ham * 219474.63)
        )
        return

        model.get_hamiltonian([
            state(1)

        ])

        stretch, angles, dihedrals = coordops.get_stretch_coordinate_system([tuple(s[:2]) for s in pb.bonds])
        labels = pb.edge_graph.get_label_types()
        stretch_types = [
            coordops.get_coordinate_label(
                c,
                labels
            )
            for c in stretch
        ]
        bend_types = [
            coordops.get_coordinate_label(
                c,
                labels
            )
            for c in angles
        ]

        good_coords = {
            c:l
            for c, l in zip(stretch, stretch_types)
            if l.atoms == 'CH'
        }

        good_coords.update({
            c:l
            for c,l in zip(angles, bend_types)
            if l.atoms == 'HCH'
        })

        nms = pb.get_normal_modes()
        loc_modes = nms.localize(internals=good_coords).make_oblique()
        base_hess = loc_modes.compute_hessian()

        # print(stretch_types)
        # print(bend_types)


        print("Base Oblique Hessian")
        print(TableFormatter("{:.0f}").format(base_hess * 219474.63))

        print("Scaled Oblique Hessian")
        scaled_hess = modify_internal_hamiltonian(
            base_hess,
            good_coords,
            scaling_types={
                ("CH", "stretch"):.96,
                (("methyl", "CH", "stretch"), ("methyl", "CH", "stretch")):.9,
                (("ethyl", "CH", "stretch"), ("ethyl", "CH", "stretch")):.96,
            }
        )
        print(TableFormatter("{:.0f}").format(scaled_hess * 219474.63))


    @validationTest
    def test_AtomTypeMap(self):
        import McUtils.Coordinerds as coordops
        # pb = Molecule.from_file(
        #     TestManager.test_data('proplybenz.hess')
        # )
        pb = Molecule.construct('c1ccccc1CCOCC(=O)O')
        # pb = Molecule.construct('CC(C)(C)O')
        labels = pb.edge_graph.get_label_types()

        for b in pb.bonds:
            print(
                b[:2],
                coordops.get_coordinate_label(
                    b[:2],
                    labels
                )
            )
        return

        g = pb.edge_graph
        print(g.find_functional_groups())
        print([
            g.categorize_ring(r)
            for r in g.rings
        ])
        # print(
        #     pb.edge_graph.get_label_types(neighbor_depth=2)
        # )
        return

        r = g.get_rings()
        # return
        # woof = (pb.edge_graph.get_rings())
        # print(pb.atoms)
        # print([len(w) for w in woof])
        # return


        print(
            pb.edge_graph.get_label_types(neighbor_depth=1)
        )
        pb = Molecule.from_file(
            TestManager.test_data('methanol_vpt_3.fchk')
        )
        print(
            pb.edge_graph.get_label_types(neighbor_depth=1)
        )
        print(
            pb.edge_graph.get_label_types(neighbor_depth=2)
        )

        pb = Molecule.construct('acetic acid')
        print(
            pb.edge_graph.get_label_types(neighbor_depth=2)
        )

    @validationTest
    def test_RedundantConversion(self):

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
                # 'primitives': 'auto',
                'primitives': [(0,1), (0, 2), (0, 3), (15, 10, 11, 13)],
                'untransformed_coordinates': [(15, 10, 11, 13)]
            }
        )

        print(
            gggg.internal_coordinates.converter_options['redundant_transformation']
        )

    @validationTest
    def test_InternalConv(self):
        test_inverse = False
        test_vpt = False
        test_direct_expansions = False
        test_numerical_expansions = ['all']

        # nh3 = Molecule.from_file(
        #     TestManager.test_data("HOH_freq.fchk"),
        #     internals=[
        #         [0, -1, -1, -1],
        #         [1,  0, -1, -1],
        #         [2,  0,  1, -1],
        #         # [3,  0,  1,  2]
        #     ]
        # )
        # nh3 = Molecule(
        #     atoms=nh3.atoms[:2],
        #     coords=nh3.coords[:2],
        #     internals=[
        #         [0, -1, -2, -3],
        #         [1,  0, -1, -2]
        #     ]
        # )
        #
        nh3 = Molecule.from_file(
            TestManager.test_data("nh3.fchk"),
            internals=[
                [0, -1, -1, -1],
                [1,  0, -1, -1],
                [2,  0,  1, -1],
                [3,  0,  1,  2]
            ]
        )
        nh3_2 = nh3.copy()

        # nh3 = Molecule.from_file(
        #     TestManager.test_data("OCHH_freq.fchk"),
        #     internals=[
        #         [0, -1, -1, -1],
        #         [1,  0, -1, -1],
        #         [2,  1,  0, -1],
        #         [3,  1,  2,  0]
        #     ]
        # )

        #
        test_coords = (
                # np.random.rand(1, 3) +
                nh3.coords
                # * np.array([2, 1, 2, .5])[:, np.newaxis]
        )

        if test_direct_expansions:
            with np.printoptions(linewidth=1e8, suppress=True):
                e1 = nput.angle_vec(test_coords, 0, 1, 2, method='old', angle_ordering='jik', order=2)
                e2 = nput.angle_vec(test_coords, 0, 1, 2, method='expansion', angle_ordering='jik', order=2)
                for i1, i2 in zip(e1, e2):
                    print(":::" * 50)
                    print(":::" * 50)
                    print(i1.shape)
                    w = np.where(np.abs(i1 - i2) > 1e-6)
                    print(w)
                    if len(w[0]) > 0:
                        print("=" * 50)
                        print(np.round(i1, 6))
                        print("-" * 50)
                        print(np.round(i2, 6))
                        print("=" * 50)
                    # i1[np.abs(i1) < 1e-8] = 1
                    # print(np.round(i2 / i1, 8))

            # print("=" * 20)
            # print(nput.pts_dihedrals(test_coords[3], test_coords[0], test_coords[1], test_coords[2]))
            e1 = nput.dihed_vec(test_coords, 3, 0, 1, 2, method='old', order=2)
            e2 = nput.dihed_vec(test_coords, 3, 0, 1, 2, method='expansion', order=2)
            with np.printoptions(linewidth=1e8, suppress=True):
                for i1, i2 in zip(e1, e2):
                    print(":::" * 50)
                    print(":::" * 50)
                    print(i1.shape)
                    w = np.where(np.abs(i1 - i2) > 1e-6)
                    print(w)
                    if len(w[0]) > 0:
                        print("=" * 50)
                        print(np.round(i1, 6))
                        print("-" * 50)
                        print(np.round(i2, 6))
                        print("=" * 50)
                    # i1[np.abs(i1) < 1e-8] = 1
                    # print(np.round(i2 / i1, 8))

        # #
        # print("=" * 20)
        # e1 = np.round(nput.dihed_vec(test_coords, 1, 3, 2, 0, method='old'), 8)
        # e2 = np.round(nput.dihed_vec(test_coords, 1, 3, 2, 0, method='expansion'), 8)
        # print(e1.reshape(-1, 3))
        # print(e2.reshape(-1, 3))
        # e1[e1 == 0] = 1
        # print(e2 / e1)
        # #
        # print("=" * 20)
        # e1 = np.round(nput.dihed_vec(test_coords, 1, 3, 0, 2, method='old'), 8)
        # e2 = np.round(nput.dihed_vec(test_coords, 1, 3, 0, 2, method='expansion'), 8)
        # print(e1.reshape(-1, 3))
        # print(e2.reshape(-1, 3))
        # e1[e1 == 0] = 1
        # print(e2 / e1)
        #

        if test_numerical_expansions:
            if isinstance(test_numerical_expansions, bool):
                test_numerical_expansions = ['dihedrals', 'all']
            exp0 = nh3.get_internals_by_cartesians(4,
                                                   strip_embedding=True,
                                                   analytic_derivative_order=0,
                                                   # stencil=9,
                                                   # mesh_spacing=5e-4,
                                                   all_numerical=True)

            if 'dihedrals' in test_numerical_expansions:
                with np.printoptions(linewidth=1e8, suppress=True):
                    # e2 = nput.angle_vec(test_coords, 0, 2, 1, angle_ordering='jik', method='old', order=2)
                    # e1 = [e[..., 5] for e in exp0]
                    # e2 = nput.dihed_vec(test_coords, 3, 0, 1, 2, method='expansion', order=3)[1:]
                    # e1 = [e[..., 2] for e in exp0]
                    # e2 = nput.angle_vec(test_coords, 2, 0, 1, method='expansion', order=3)[1:]

                    # e1 = [e[..., 1] for e in exp0]
                    # e2 = nput.dist_vec(test_coords, 0, 2, method='expansion', order=3)[1:]

                    e1 = [e[..., -1] for e in exp0]
                    # e2 = nput.dihed_vec(test_coords, 3, 0, 1, 2, method='expansion', order=3)[1:]

                    e2 = nput.dihed_vec(test_coords, *nh3.internals['zmatrix'][-1], method='old', order=2)[1:]
                    for i1, i2 in zip(e1, e2):
                        print(":::" * 50)
                        print(":::" * 50)
                        print(i1.shape)
                        i1[np.abs(i1) < 1e-6] = 0
                        i2[np.abs(i2) < 1e-6] = 0

                        d1 = -i1.copy() # add a sign flip because the dihed. derivs I have are flipped
                        d1[np.abs(d1) < 1e-6] = 1
                        r = i2 / d1
                        r[np.logical_and(
                            np.abs(i1) < 1e-6,
                            np.abs(i2) < 1e-6
                        )] = 1
                        w = np.where(np.abs(r - 1) > 1e-2)
                        r[np.logical_and(
                            np.abs(i1) < 1e-6,
                            np.abs(i2) < 1e-6
                        )] = 0
                        print(w)
                        if len(w[0]) > 0:
                            print("=" * 50)
                            print(np.round(i1, 6)[2])
                            print("." * 10)
                            print(np.round(i2, 6)[2])
                            # print("-" * 50)
                            # print(np.round(i1, 6)[5])
                            # print("." * 10)
                            # print(np.round(i2, 6)[5])
                            print("=" * 50)
                            print(np.round(r, 6))
                # return
                #
                # # print(exp1[2][0, 0, 1], exp1[2][0, 1, 0], exp1[2][1, 0, 0])
                # # return
                #
            if 'all' in test_numerical_expansions:
                exp1 = nh3_2.get_internals_by_cartesians(4,
                                                         strip_embedding=True,
                                                         analytic_derivative_order=-1,
                                                         all_numerical=False)
                for x in range(6):
                    print("##" * 50)
                    print("TESTING:", x)
                    print("##" * 50)
                    e1 = [e[..., x] for e in exp0]
                    # e2 = nput.dihed_vec(test_coords, 3, 0, 1, 2, method='expansion', order=3)[1:]

                    e2 = [e[..., x] for e in exp1]
                    for i1, i2 in zip(e1, e2):
                        print(":::" * 50)
                        print(":::" * 50)
                        print(i1.shape)
                        i1[np.abs(i1) < 1e-6] = 0
                        i2[np.abs(i2) < 1e-6] = 0

                        d1 = (-1 if x == 5 else 1 ) * i1.copy()
                        d1[np.abs(d1) < 1e-6] = 1
                        r = i2 / d1
                        r[np.logical_and(
                            np.abs(i1) < 1e-6,
                            np.abs(i2) < 1e-6
                        )] = 1
                        w = np.where(np.abs(r - 1) > 1e-2)
                        r[np.logical_and(
                            np.abs(i1) < 1e-6,
                            np.abs(i2) < 1e-6
                        )] = 0
                        print(w)
                        if len(w[0]) > 0:
                            print("=" * 50)
                            print(np.round(i1, 6)[2])
                            print("." * 10)
                            print(np.round(i2, 6)[2])
                            # print("-" * 50)
                            # print(np.round(i1, 6)[5])
                            # print("." * 10)
                            # print(np.round(i2, 6)[5])
                            print("=" * 50)
                            print(np.round(r, 6))

        if test_inverse:
            ## TEST INVERSE QUALITY
            with np.printoptions(linewidth=1e8, suppress=True):
                # exp0 = nh3.get_internals_by_cartesians(3,
                #                                        strip_embedding=True,
                #                                        analytic_derivative_order=0,
                #                                        all_numerical=True
                #                                        )
                # inv0 = nh3.get_cartesians_by_internals(3,
                #                                        strip_embedding=True,
                #                                        reembed=True,
                #                                        analytic_derivative_order=0,
                #                                        all_numerical=True,
                #                                        method='old')

                exp0 = nh3.get_internals_by_cartesians(3,
                                                       strip_embedding=True
                                                       )
                inv0 = nh3.get_cartesians_by_internals(3,
                                                       strip_embedding=True,
                                                       reembed=True,
                                                       method='fast')

                # inv1 = nh3.get_cartesians_by_internals(3,
                #                                        strip_embedding=True,
                #                                        reembed=True,
                #                                        analytic_derivative_order=-1,
                #                                        all_numerical=False,
                #                                        method='fast')

                #     for i1,i2 in zip(inv0, inv1):
                #         print(":::"*50)
                #         print(":::"*50)
                #         print(np.where(np.abs(i1 - i2) > 1e-4))
                #         print("="*50)
                #         print(np.round(i1[0], 6))
                #         print("-"*50)
                #         print(np.round(i2[0], 6))
                #         print("="*50)
                #         print(np.round(i1 - i2, 6)[0])
                for t in nput.tensor_reexpand(inv0, exp0):
                    print(
                        np.round(t, 5)
                        # np.round(exp0[2] - exp1[2], 8)
                    )
                # for t in nput.tensor_reexpand(inv1, exp1):
                #     print(
                #         np.round(t, 8)
                #         # np.round(exp0[2] - exp1[2], 8)
                #     )

            # return

        if test_vpt:

            nh3.setup_VPT(states=1,
                          logger=False,
                          cartesian_analytic_deriv_order=-1,
                          cartesian_by_internal_derivative_method='fast',
                          )[0].print_tables(print_intensities=False, print_energy_corrections=False)
            nh3.setup_VPT(states=1,
                          logger=False,
                          cartesian_analytic_deriv_order=0,
                          cartesian_by_internal_derivative_method='old',
                          )[0].print_tables(print_intensities=False, print_energy_corrections=False)
            nh3 = Molecule.from_file(nh3.source_file)
            nh3.setup_VPT(states=1, logger=False)[0].print_tables(print_intensities=False, print_energy_corrections=False)

        #
        #
        # emb_nh3 = nh3.get_embedded_molecule()
        # emb_test = emb_nh3.internal_coordinates.convert(emb_nh3.coords.system)
        # # conv_1 = nh3.internal_coordinates.convert(nh3.coords.system)
        # # print(emb_nh3.internal_coordinates)
        # # print(emb_nh3.internal_coordinates.converter_options)
        # # print(np.round(emb_nh3.coords, 8))
        # # print(np.round(emb_test, 8))
        # # raise Exception(...)
        # # conv_2 = ...
        #
        # wtf = nh3.get_cartesians_by_internals(1,
        #                                     strip_embedding=True,
        #                                     reembed=True,
        #                                     analytic_derivative_order=-1,
        #                                     all_numerical=False,
        #                                     method='classic')
        #
        # nh3_derivs_internal = [
        #     0,
        #     wtf[0] @ nh3.potential_derivatives[1] @ wtf[0].T
        # ]
        #
        # nh3_gmatrix = nh3.get_gmatrix(
        #     analytic_derivative_order=-1,
        #     all_numerical=False
        # )
        #
        # freqs_int, _ = scipy.linalg.eigh(nh3_derivs_internal[1], nh3_gmatrix, type=2)
        # freqs_cart, _ = scipy.linalg.eigh(nh3.potential_derivatives[1], nh3.get_gmatrix(use_internals=False), type=2)
        #
        # print(np.sign(freqs_int) * np.sqrt(np.abs(freqs_int)) * UnitsData.convert("Hartrees", "Wavenumbers"))
        # print(np.sqrt(freqs_cart[6:]) * UnitsData.convert("Hartrees", "Wavenumbers"))
        #
        # return

        # methanol = Molecule.from_file(
        #     TestManager.test_data("methanol_vpt_1.fchk"),
        #     internals=[
        #         [0, -1, -2, -3],
        #         [1, 0, -1, -2],
        #         [2, 1, 0, -1],
        #         [3, 2, 1, 0],
        #         [4, 2, 3, 1],
        #         [5, 2, 3, 4]
        #     ]
        # )
        # ders1 = methanol.get_cartesians_by_internals(1, strip_embedding=True, reembed=True)[0]
        # ders_inv1 = methanol.get_internals_by_cartesians(1, strip_embedding=True)[0]
        #
        # # nh3 = Molecule.from_file(
        # #     TestManager.test_data("nh3.fchk"),
        # #     internals=[
        # #         [0, -1, -1, -1],
        # #         [1,  0, -1, -1],
        # #         [2,  0,  1, -1],
        # #         [3,  0,  1,  2]
        # #     ]
        # # )
        # ders2 = methanol.get_cartesians_by_internals(1, method='classic', strip_embedding=True, reembed=True)[0]
        # # print(np.round(ders1, 7)[0])
        # # print(np.round(ders2, 7)[0])
        #
        # print(np.round(ders2 @ ders_inv1, 6))
        # print(np.round(ders1 @ ders_inv1, 8))
        #
        # print(np.round(ders2 - ders1, 5))
        #
        #
        #
        # meth_derivs_internal = nput.tensor_reexpand(
        #     methanol.get_cartesians_by_internals(2, strip_embedding=True, reembed=True),
        #     [0, methanol.potential_derivatives[1]]
        # )
        #
        # meth_gmatrix = methanol.g_matrix
        #
        # freqs_int, _ = scipy.linalg.eigh(meth_derivs_internal[1], meth_gmatrix, type=2)
        # freqs_cart, _ = scipy.linalg.eigh(methanol.potential_derivatives[1], methanol.get_gmatrix(use_internals=False), type=2)
        #
        # print(np.sqrt(freqs_int) * UnitsData.convert("Hartrees", "Wavenumbers"))
        # print(np.sqrt(freqs_cart[6:]) * UnitsData.convert("Hartrees", "Wavenumbers"))
        #
        # return

    @validationTest
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



    @validationTest
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


        # # sys = 'HOONO_freq.fchk'
        # # internals = [
        # #             [1, -1, -1, -1],
        # #             [2, 1, -1, -1],
        # #             [3, 2, 1, -1],
        # #             [0, 1, 2, 3],
        # #             [4, 3, 2, 1]
        # #         ]
        # sys = 'nh3.fchk'
        # internals = [
        #     [0, -1, -1, -1],
        #     [1,  0, -1, -1],
        #     [2,  0,  1, -1],
        #     [3,  0,  1,  2]
        # ]
        # mol2 = Molecule.from_file(
        #     TestManager.test_data(sys),
        #     internals={
        #         'zmatrix': internals
        #     }
        # )
        # mol2.embedding.cartesian_by_internals_method = 'fast'
        # mol3 = Molecule.from_file(
        #     TestManager.test_data(sys),
        #     internals={
        #         'zmatrix': internals
        #     }
        # )
        # mol3.embedding.cartesian_by_internals_method = 'numerical'
        # # print(mol2.hamiltonian.potential_expansion(3)[2][0])
        # # print(mol3.hamiltonian.potential_expansion(3)[2][0])
        # # raise Exception(...)
        # L_base = mol2.embedding.get_translation_rotation_invariant_transformation(mass_weighted=False)
        # jacs_1 = mol2.get_internals_by_cartesians(2, strip_embedding=True)
        # new_tf = nput.tensor_reexpand([L_base.T], jacs_1, 2)
        # inverse_tf = nput.inverse_transformation(new_tf, 2)
        # jacs_2_test = [
        #     np.tensordot(j, L_base, axes=[-1, -1])
        #     for j in inverse_tf
        # ]
        # jacs_2 = mol2.embedding.get_cartesians_by_internals(2, reembed=True, strip_embedding=True, method='og')
        # i,j = [0, 2]
        # print(np.round(jacs_2[1][i, j], 8))
        # print(np.round(jacs_2_test[1][i, j], 8))
        # print(np.round(jacs_2[1][i, j] - jacs_2_test[1][i, j], 8))
        # # print(jacs_2[1].shape)
        # woof = np.array(np.where(np.abs(jacs_2[1] - jacs_2_test[1]) > 1e-7)).T
        # print(np.unique(woof[:, :2], axis=0))
        # # print(len(woof[0]))
        # raise Exception(...)


    @validationTest
    def test_MoreBondZMatrix(self):
        import McUtils.Coordinerds as coordops

        for smi in [
            'CC1(C)C(/C=C/C2=C(O)C=CC3=C2C=CC=C3)=[N+](CCCS(=O)([O-])=O)C4=CC=CC=C41',
            'OC(C=CC=C1)=C1/C=C/C2=[N+](CCCS(=O)([O-])=O)C3=CC=CC=C3C2(C)C',
            'COC1=CC=C(C2=C1)[N+](CCCOS([O-])=O)=C(C2(C)C)/C=C/C3=C(C=C(C=C3)OC)O'
        ]:
            test = Molecule.from_string(smi)
            z = test.get_bond_zmatrix()
            coords = coordops.extract_zmatrix_internals(z)
            if None in coords:
                for c in z:
                    for i,j in itertools.combinations(c, 2):
                        if i == j:
                            raise ValueError(c)
            pprint.pprint(z)

    @validationTest
    def test_EvenMoreZMatrix(self):
        import McUtils.Coordinerds as coordops
        # ts = Molecule.from_file(TestManager.test_data('ts_samp.xyz'))
        # ts.get_bond_zmatrix(for_fragment=0)

        ts = Molecule.from_file(TestManager.test_data('ts_samp2.xyz'))
        # ts.get_bond_zmatrix()

        zm_sub = ts.modify(
            bonds=[b for b in ts.bonds if b[0] not in {19, 18} or b[1] not in {19, 18}]
        ).get_bond_zmatrix(
            for_fragment=ts.fragment_indices[1],
            fragment_ordering=[0, 1],
            attachment_points={19:18}
        )
        # pprint.pprint(zm_sub)

        f1 = ts.fragments[1]
        zm_f1 = f1.modify(
            bonds=[b for b in f1.bonds if b[0] not in {0, 1} or b[1] not in {0, 1}]
        ).get_bond_zmatrix(
            fragment_ordering=[0, 1],
            attachment_points={1: 0}
        )

        f1_int = f1.modify(internals={'specs':coordops.extract_zmatrix_internals(zm_f1)})
        woof = f1_int.get_cartesians_by_internals(1)
        # print(woof[0].shape)

    @validationTest
    def test_RDKitInputFormats(self):
        Molecule.from_string('MDSKGSGS', 'fasta')

    @validationTest
    def test_BreakBondZMat(self):
        import McUtils.Coordinerds as coordops
        # ts = Molecule.from_file(TestManager.test_data('ts_samp.xyz'))
        # ts.get_bond_zmatrix(for_fragment=0)

        ts = Molecule.from_file(TestManager.test_data('ts_samp.xyz'))
        ts_mod = ts.break_bonds([(18, 22)])
        zm = ts_mod.get_bond_zmatrix()
        ts_int = ts.modify(internals=zm)
        # ts_int.animate_coordinate(0)
        # pprint.pprint(ts_mod.get_bond_zmatrix())

    @validationTest
    def test_NeverEndingZMatrix(self):
        # import McUtils.Coordinerds as coordops
        # ts = Molecule.from_file(TestManager.test_data('ts_samp.xyz'))
        # ts.get_bond_zmatrix(for_fragment=0)

        react = Molecule.from_file(TestManager.test_data('test_react.xyz'), units='Angstroms')
        react.get_bond_zmatrix()

    @validationTest
    def test_TRIC(self):
        # import McUtils.Coordinerds as coordops
        # ts = Molecule.from_file(TestManager.test_data('ts_samp.xyz'))
        # ts.get_bond_zmatrix(for_fragment=0)

        react = Molecule.from_file(TestManager.test_data('test_react.xyz'), units='Angstroms')
        tr_react = react.modify(internals={
            'specs':[
                {'transrot':(0, 1, 2, 3)},
                {'orientation': ((0, 1, 4), (2, 3, 5))}
            ]
        })
        uuuh = tr_react.get_internals_by_cartesians(1)
        tr_react.get_cartesians_by_internals(1)
        tr_react.animate_coordinate(0, backend='x3d', highlight_atoms=[0, 1, 2, 3])

    @validationTest
    def test_NewAnim(self):
        ohh = Molecule.from_file(TestManager.test_data('water_freq.fchk'))
        ohh.animate_mode(0, backend='x3d')

    @validationTest
    def test_RDKitIssues(self):
        from Psience.Molecools import Molecule

        Molecule.from_string('COc1cc([OH]C2([N+](c3c(C2(C)C)cc(OC)cc3)CCCOS([O-])=O)C=C4)c4cc1', 'smi')

    @validationTest
    def test_RDKitConfGen(self):
        from Psience.Molecools import Molecule

        mols = Molecule.from_string(
            'Cc1ccc(OC[C:6]([NH:5]/[N:4]=[CH:3]/[c:2]2cc(=O)[nH]c(=O)[nH:1]2)=[O:7])c([N+](=O)[O-])c1',
            'smi',
            confgen_opts={
                'distance_constraints': {(4, 5): [1.2, 1.4]},
                'random_seed': 100
            }
        )

        mols = Molecule.from_string(
            'Cc1ccc(OC[C:6]([NH:5]/[N:4]=[CH:3]/[c:2]2cc(=O)[nH]c(=O)[nH:1]2)=[O:7])c([N+](=O)[O-])c1',
            'smi',
            num_confs=10,
            confgen_opts={
                'distance_constraints': {(4, 5): [1.2, 1.4]},
                'random_seed': 100
            }
        )

        raise Exception(mols)

        print(
            Molecule(
                ["O", "H", "H"],
                [[0, 0, 0], [0, 1, 0], [1, 0, 0]]
            ).bonds
        )

        mol = Molecule.from_string(
            'Cc1ccc(OC[C:6]([NH:5]/[N:4]=[CH:3]/[c:2]2cc(=O)[nH]c(=O)[nH:1]2)=[O:7])c([N+](=O)[O-])c1',
            'smi',
            confgen_opts={
                'distance_constraints': {(4, 5): [1.2, 1.4]},
                'random_seed': 100
            }
        )

        print(
            Molecule(
                mol.atoms,
                mol.coords
            ).bonds[:3]
        )

        print(
            np.linalg.norm(mol.coords[4] - mol.coords[5]) * UnitsData.bohr_to_angstroms
        )

        Molecule.from_string(
            "NC(N)=NC(=O)c1[nH]c(C(=O)O)c2cc3c(cc12)c1c2cc4c(C(=O)O)[nH]c(C(=O)N=C(N)N)c4cc2c2c4cc5c(C(=O)O)[nH:1][c:2]([C:3](=O)[N:4]=[C:5](N)[NH2:6])c5cc4c4c5cc6c(C(=O)O)[nH]c(C(=O)N=C(N)N)c6cc5c5c6cc7c(C(=O)O)[nH]c(C(=O)N=C(N)N)c7cc6c3c3c1c2c4c53",
            "smi"
        )

    @validationTest
    def test_CanonicalZMatrix(self):
        from Psience.Molecools import Molecule

        # mol = Molecule.from_string(
        #     'Cc1ccc(OCC(NN=CC2CC(=O)NC(=O)N2)=O)c(N(=O)O)c1',
        #     add_implicit_hydrogens=False
        # )
        #
        # smi = mol.to_string('smi', include_tag=True, remove_hydrogens=True)
        # mol2 = Molecule.from_string(smi, 'smi', add_implicit_hydrogens=False)

        # from Psience.Molecools import Molecule
        #
        # mol = Molecule.from_string(
        #     'Cc1ccc(OCC(NN=CC2CC(=O)NC(=O)N2)=O)c(N(=O)O)c1',
        #     add_implicit_hydrogens=True
        # )
        #
        # smi = mol.to_string('smi', preserve_atom_order=True, include_tag=True, remove_hydrogens=True)
        # print(smi.partition("_")[0])
        # mol2 = Molecule.from_string(smi, 'smi', add_implicit_hydrogens=True)
        #
        # mol = mol.get_embedded_molecule()
        # mol2 = mol2.get_embedded_molecule(ref=mol)
        import McUtils.Coordinerds as coordops


        test_smi = 'Cc1ccc(OC[C]([NH]/[N]=[CH]/[c]2cc(=O)[nH]c(=O)[nH]2)=[O])c([N+](=O)[O-])c1'
        mol = Molecule.from_string(
            test_smi
        )
        smi = mol.to_string('smi', include_tag=True, remove_hydrogens=True)

        mol2 = mol.from_string(smi, 'smi', add_implicit_hydrogens=True)

        # yeesh2 = mol2.modify(internals=mol2.get_canonical_zmatrix())
        #
        # zmat = mol.get_canonical_zmatrix()
        # yeesh = mol.modify(internals=zmat)
        # ca = yeesh.internal_coordinates
        # cb = yeesh2.internal_coordinates
        #
        # print(coordops.set_zmatrix_embedding(ca[:25] - cb[:25], embedding=np.zeros(6)))


    @validationTest
    def test_FlexiblePlotting(self):
        from Psience.Molecools import Molecule

        mol = Molecule.from_string('O=[C:2]([NH:1]c1cccc(C(F)(F)F)c1)[O:3]/[N:4]=[C:5](\C1COc2ccccc2O1)[NH2:6]', 'smi')
        # smi = mol.to_string('smi', preserve_atom_order=True, include_tag=True, remove_hydrogens=True)
        # mol2 = Molecule.from_string(smi, 'smi')
        fig1 = mol.plot(
            highlight_bonds=[(0, 1), (2, 3), (4, 5)],
            bond_style={(0,1):{'color':'blue'}},
            bond_radius=5,
            use_default_radii=False,
            # include_script_interface=True,
            atom_style={i: {"color": "#FF00FF"} for i in range(5)},
            # background='blue',
            image_size=[800, 500],
            background='blue',
            atom_labels={
                i: {}
                for i in range(27)
            },
            draw_coords={
                (0, 2):{
                    'label':'r<msub>1</msub>'
                },
                (0, 1, 2):{
                    'scaling':.5,
                    # 'label':{'text':'r<msub>1</msub>', 'color':'green'},
                    'styles':{'color':'red'}
                },
                (3, 4, 5): {
                    'scaling': 1,
                    'label': {'text': 'a', 'offset': [3.1, 0]}
                },
                (5, 6, 7):{
                    'color':None,
                    'label':{'text':'6'}
                }
            },
            # atom_note_color='gray',
            label_style={
                'color': 'gray',
                'font_size': 7
            },
            # drawing_extents_include=['all'],
            backend='2d'
        )

        # mol.plot(
        #     figure=fig1,
        #     background='#FFFFFFFF',
        #     backend='2d'
        # )

    @validationTest
    def test_StableInternals(self):
        mol = Molecule.from_file('water_freq.fchk',
                                   internals=[
                                       [0, -1, -2, -3],
                                       [1,  0, -1, -2],
                                       [2,  0,  1, -1],
                                   ])

        print(mol.get_internals_by_cartesians())

    @validationTest
    def test_Opts(self):
        import McUtils.Coordinerds as coordops
        import McUtils.Formatters as mfmt
        import numpy as np
        np.seterr(all='raise')

        # coordops.InternalCoordinateType.resolve(
        #     {"orientation":((0, 1, 2), (3, 4, 5))}
        # )
        #
        # return

        mol = Molecule.from_file(
                TestManager.test_data('OCHH_freq.fchk'),
                energy_evaluator='aimnet2'
            )
        fig = mol.plot(use_default_bonds=False)
        mol.plot(use_default_bonds=False, figure=fig)
        return

        # mol = Molecule.from_file(
        #     TestManager.test_data('OCHH_freq.fchk'),
        #     energy_evaluator='aimnet2'
        # )
        # mol = Molecule.from_file(
        #     TestManager.test_data('tbhp_180.fchk')
        # )
        mol = Molecule.from_file(
            TestManager.test_data('react_samp.xyz'),
            energy_evaluator='rdkit'
        )
        coordops.validate_zmatrix(mol.get_bond_zmatrix())
        zmcs = mol.get_bond_zmatrix()

        n = 4
        constraints = [
                z
                    for z in coordops.extract_zmatrix_internals(zmcs)
                if len(z) == n
            ][:3]

        # bases, _, _ = nput.internal_basis(mol.coords, constraints)
        # # base_tensors = tf_fun(coords)
        # # if mask is not None:
        # #     checks = mask(base_tensors[0])
        # #     #TODO: handle the mixed-shape case
        # #     base_tensors = base_tensors[1][..., :, checks]
        # base_tensors = np.concatenate(bases, axis=-1)
        # a, b = nput.orthogonalize_transformations([
        #     [[base_tensors.T], [base_tensors]]
        # ])
        #
        # # print(a[0] @ b[0])
        # # return
        #
        # # print(bases[0])
        # proj = nput.orthogonal_projection_matrix(base_tensors, orthonormal=False)
        #
        # rdm = np.random.rand(*mol.coords.shape).flatten()[np.newaxis] @ proj
        #
        # disps = mol.get_scan_coordinates(
        #     [[-1, 1, 3]],
        #     which=[0],
        #     coordinate_expansion=[rdm]
        # )
        #
        # di = mol.modify(internals=zmcs).get_internals(coords=disps, strip_embedding=False)
        # print(di.shape)
        # idx = coordops.zmatrix_indices(zmcs, constraints)
        #
        # pop_vals = lambda ints: coordops.extract_zmatrix_values(ints, idx)
        #
        # print(pop_vals(di[1]))
        # print(pop_vals(di[0]))
        #
        # return


        opt = mol.optimize(
            coordinate_constraints=constraints
        )
        opt2 = mol.optimize(
            # coordinate_contraints=[
            #     z
            #     for z in coordops.extract_zmatrix_internals(zmcs)
            #     if len(z) == 4
            # ][:3]
        )
        print(constraints)
        idx = coordops.zmatrix_indices(zmcs, constraints)
        pop_vals = lambda m:coordops.extract_zmatrix_values(m.modify(internals=zmcs).internal_coordinates, idx)

        print(pop_vals(mol))
        print(pop_vals(opt))
        print(pop_vals(opt2))

        # hmm = mol.modify(internals=mol.get_bond_zmatrix()).get_scan_coordinates(
        #     [[-.05, .05, 3]],
        #     which=[5],
        #     internals='reembed',
        #     shift=True,
        #     strip_embedding=True
        # )
        # print(hmm[1] - mol.coords)
        return

        # zmat = mol.get_bond_zmatrix()
        # print()
        # print(mfmt.format_zmatrix(zmat))
        # print(mol.fragment_indices)
        # return

        # int_tbhp = mol.modify(internals=mol.get_bond_zmatrix())
        # dx = int_tbhp.get_cartesians_by_internals(order=1)[0]

        # u_traj = []
        # ugh = mol.optimize(
        #     # func=lambda s:(-np.dot(s, dx[6])),
        #     gradient_modification_function=lambda c,g:g+.2*dx[6].reshape(g.shape),
        #     # initialization_function=lambda c:c+5*dx[6].reshape(-1, 3),
        #     return_trajectory=True,
        #     line_search=False,
        #     mode='scipy'
        #     # logger=True
        # )
        # print(
        #     len(ugh[1])
        # )

        zmat = mol.get_bond_zmatrix(validate=True)

        # spec = coordops.InternalSpec(coordops.extract_zmatrix_internals(zmat))
        # _, inv = spec.get_expansion(mol.coords, order=1, return_inverse=True, orthogonalize=False)
        # print(inv[0])
        # return

        int_mol = mol.modify(internals=zmat)
        # print(int_mol.internal_coordinates)
        print(mol.coords[:5])
        exp = int_mol.get_cartesians_by_internals(
            method='classic',
            use_direct_expansions=True,#[2],
            orthogonalize_derivatives=False,
            allow_fd=False,
            order=1,
            strip_embedding=False
        )
        # print(exp[0][0])
        exp = mol.modify(internals=zmat).get_cartesians_by_internals(
            # method='classic',
            # use_direct_expansions=True,
            # orthogonalize_derivatives=False,
            allow_fd=False,
            order=1,
            strip_embedding=True
        )
        # print(exp[0][0])
        return
        # raise Exception(exp[0].shape)

        mol = Molecule.from_file(
            TestManager.test_data('tbhp_180.fchk'),
            energy_evaluator='aimnet2'
        )

        ugh =  mol.modify(internals=zmat).optimize(
            coordinate_constraints=[
                c
                for c in coordops.extract_zmatrix_internals(zmat)
                if len(c) == 4
            ],
            track_best=True,
            max_iterations=50
        )

        print(
            np.concatenate([
                ugh.modify(internals=zmat).internal_coordinates,
                mol.modify(internals=zmat).internal_coordinates
                ], axis=-1),
        )
        print(
            (ugh.calculate_energy() - mol.calculate_energy())*UnitsData.convert("Hartrees", "Kilocalories/Mole")
        )

    @validationTest
    def test_FragBaseDraw(self):
        bits = Molecule.from_string(
            'FC1CCC(C(=C)C)CC1.C1OCC(C(OC)O)C1C(OC)O',
            # 'c1ccccn1',
            'smi',
            confgen_opts=dict(verbose=True, random_seed=12321)
        ).fragments[0].get_embedded_molecule(
            # sel=[5, 6, 7]
        )
        bits.plot(
            mode=(
                # 'matplotlib3D'
                'x3d'
            ),
            principle_axes=True,
            image_size=[500, 500],
            theme='flat',
            draw_coords={
                (0, 1, 2):{'label':r'$\theta$',
                           'label_style':{
                               'offset_magnitude':1.1,
                               'billboard':True,
                               'color':'red', 'fontsize':24
                           }}
            },
            highlight_atoms=[0, 1, 2],
            view_settings={
                # 'view_distance':200,
                'view_vector': [0, 0, 1],
            }).show()#.savefig("/Users/Mark/Desktop/view_xy_simp_bonds.svg")

    @validationTest
    def test_FragInternalsSpec(self):
        import McUtils.Coordinerds as coordops

        ints3 = [[2.59092749, 0., 0.],
             [2.89495928, 1.93458804, 0.],
             [2.89963503, 1.8239319, 3.37828743],
             [2.83053209, 1.83815421, 5.06413309],
             [2.79207623, 2.00065824, 3.29813193],
             [2.52079324, 2.15594405, 4.10504516],
             [2.82036909, 2.05789597, 0.9634783],
             [2.82837755, 1.91974668, 1.11715569],
             [2.8622616, 1.93798376, 5.26578887],
             [2.09534423, 1.97234113, 4.22162045],
             [2.0556097, 1.98211883, 5.54531358],
             [2.12037191, 1.97468473, 4.16514767],
             [2.10238981, 1.86667929, 0.80257191],
             [2.07437367, 1.97295648, 4.17603043],
             [2.12921129, 1.8901408, 5.36779409],
             [2.04513453, 2.1175231, 3.14158509],
             [2.05011789, 2.05580992, 3.1416388],
             [2.06074897, 2.453391, 0.31696486],
             [2.09040482, 1.89999553, 2.33908385],
             [2.11914653, 1.84889244, 1.97498365],
             [2.13523637, 1.74785729, 5.16708695],
             [2.08252609, 1.9355175, 4.49096479],
             [2.07995171, 1.95503155, 5.23240844],
             [2.09065558, 1.93768086, 2.16018763]]
        zm3 = [[0, -1, -2, -3],
               [1, 0, -1, -2],
               [2, 1, 0, -1],
               [3, 2, 1, 0],
               [4, 3, 2, 1],
               [5, 4, 3, 2],
               [6, 5, 4, 3],
               [7, 5, 4, 3],
               [8, 4, 3, 2],
               [9, 8, 4, 3],
               [10, 1, 0, 2],
               [11, 2, 1, 0],
               [12, 2, 11, 1],
               [13, 3, 2, 1],
               [14, 3, 13, 2],
               [15, 4, 3, 2],
               [16, 6, 5, 4],
               [17, 6, 16, 5],
               [18, 7, 6, 5],
               [19, 7, 18, 6],
               [20, 7, 18, 19],
               [21, 8, 7, 6],
               [22, 8, 21, 7],
               [23, 9, 8, 7],
               [24, 9, 23, 8]]

        bits = Molecule.from_string(
            'FC1CCC(C(=C)C)CC1.C1OCC(C(OC)O)C1C(OC)O',
            # 'c1ccccn1',
            'smi',
            confgen_opts=dict(verbose=True, random_seed=12321)
        )#.fragments[0]

        # bits.coords[bits.fragment_indices[1], :] += 10
        # pprint.pprint(bits.get_canonical_zmatrix())
        # print(np.array(bits.get_canonical_zmatrix()[:6]))
        # print()
        # for f in bits.get_bond_zmatrix(validate=True, connect_fragments=False):
        #     print(np.array(f))
        zm = bits.get_bond_zmatrix(validate=True, connect_fragments=False)
        spec = coordops.InternalSpec.from_zmatrix(
            *(np.array(z).tolist() for z in zm)
        )

        # carts3 = coordops.zmatrix_to_cartesian(
        #     ints3,
        #     zm3
        # )
        # print(coordops.cartesian_to_zmatrix(
        #     carts3,
        #     zm3
        # ).coords)
        # raise Exception(...)

        ints = spec.cartesians_to_internals(bits.coords)
        # zints = coordops.cartesian_to_zmatrix(bits.coords, zm[0])
        # print(coordops.zmatrix_from_values(ints, partial_embedding=True))
        # print(zints.coords)
        # return

        # print(np.array(zm[0])[:4])
        # conv = coordops.find_internal_conversion(spec.rad_set, [(0, 1, 2, 3)],
        #                                          triangles_and_dihedrons=spec.get_triangulation(),
        #                                          allow_completion=False,
        #                                          prep_conversions=False)
        # print(conv(ints))
        # print(ints[:6])
        # return
        # intzm = coordops.zmatrix_from_values(ints, partial_embedding=True)
        # zmstuff3 = coordops.cartesian_to_zmatrix(bits.coords, zm[0])
        # intzm3 = zmstuff3.coords
        # c3 = coordops.zmatrix_to_cartesian(intzm3, zm[0])
        # ints2 = spec.cartesians_to_internals(c3)
        # print(ints[:9])
        # print(ints2[:9])
        # return
        zm2, _ = spec.get_zmat_conv()
        carts, exp = spec.internals_to_cartesians(ints, reference_cartesians=bits.coords, order=1)
        # print(carts - bits.coords)
        # print(exp[0][1].shape)
        uuuh = spec.get_direct_inverses(carts, order=1)
        # print(uuuh)
        print(uuuh[1][0].shape, exp[0][1].shape)
        print(np.round(uuuh[1] @ exp[0][1], 8)[:6, :6])

        # ints2 = spec.cartesians_to_internals(carts)
        # print(np.array(zm[0])[:4])
        # print(np.array(zm2)[:4])
        # print(ints[:9])
        # print(ints2[:9])
        #
        # print(np.round(ints - ints2, 8))


        return

    @validationTest
    def test_FragInternalsScan(self):
        import McUtils.Coordinerds as coordops

        mol = Molecule.from_string(
            # 'FC1CCC(C(=C)C)CC1.C1OCC(C(OC)O)C1C(OC)O',
            # 'FC1CCC(C(=C)C)CC1',
            'OC=O',
            # "N",
            'smi',
            confgen_opts=dict(verbose=True, random_seed=12321)
        )

        # mol.plot(highlight_atoms=[0, 1, 3, 4]).show()

        spec = coordops.InternalSpec.from_zmatrix(mol.get_bond_zmatrix())
        pprint.pprint(spec.rad_set)
        pprint.pprint(spec.get_triangulation())

        # graph = coordops.InternalCoordinateGraph(spec.rad_set)
        # conv = graph.find_conversions([(3, 4)])
        # print(...)
        # print(conv)
        # return


        mol = Molecule.from_string(
            # 'FC1CCC(C(=C)C)CC1.C1OCC(C(OC)O)C1C(OC)O',
            'FC1CCC(C(=C)C)CC1',
            # 'COC=O',
            # "N",
            'smi',
            confgen_opts=dict(verbose=True, random_seed=12321)
        )
        # mol.plot(backend='2d',
        #          include_save_buttons=True,
        #          draw_coords=[
        #              (0, 2),
        #              (1, 2, 3)
        #          ]).show()
        # return
        # dists, bends, diheds =  mol.get_bond_graph_internals(include_fragments=False, concatenate=False)
        # import pprint
        # pprint.pprint(diheds[:50])
        print(len(mol.atoms) * 3 - 6)


        # from McUtils.Graphs import analyze_internal_coords
        # pprint.pprint(
        #     analyze_internal_coords(dists, bends, diheds)
        # )
        # return
        # mol.internals = {
        #     'primitives':dists + bends + diheds#[:50],
        #     # 'method':'iterative'
        # }
        # mol.internals = coordops.extract_zmatrix_internals(
        #     mol.get_bond_zmatrix()
        # )

        with BlockProfiler():
            specs = mol.get_bond_graph_internals(include_fragments=False, pruning='b_matrix')
        print(len(specs))
        # specs = mol.get_bond_graph_internals(include_fragments=False, pruning='b_matrix')
        # specs = mol.get_bond_graph_internals(include_fragments=False)
        # import pprint
        # pprint.pprint([
        #     s for s in specs
        #     if all(k in [1, 2, 3, 4, 8, 9] for k in s)
        # ])

        # pprint.pprint(specs2)
        # pprint.pprint([
        #     x for x in specs if x not in specs2
        # ])
        # pprint.pprint([
        #     x for x in specs2 if x not in specs
        # ])
        pprint.pprint(specs)
        # return
        # raise Exception(len(specs))
        mol.internals = {
            'specs':specs,
            # 'method':'iterative'
        }

        # print(mol.internal_coordinates.tolist())
        # return

        # mol.animate_coordinate(-5).show()
        # return
        geoms = mol.get_scan_coordinates(
            [[-.5, .5, 5]],
            which=[-5],
            internals='reembed'
        )

        mol.plot(geoms,
                 backend='x3d',
                 highlight_atoms=[0, 2, 3],
                 image_size=800,
                 include_save_buttons=True).show()

    @validationTest
    def test_MultiXYZParsing(self):
        mol = Molecule.from_file(
            TestManager.test_data('traj.xyz'),
            units='Angstroms'
        )
        self.assertEquals(mol.coords.shape, (27, 3))


        traj = Molecule.from_file(
            TestManager.test_data('traj.xyz'),
            units='Angstroms',
            max_blocks=-1
        )
        traj[0].plot([t.coords for t in traj], image_size=800, include_save_buttons=True).show()

    @validationTest
    def test_ZMatOpt(self):
        import McUtils.Coordinerds as coordops

#         z1 = Molecule.from_string("""C
# C  1   1.4008
# C  2   1.3827  1 120.2081
# C  3   1.3804  2 121.0304  1    1.0784
# C  4   1.3971  3 119.2584  2   -2.3141
# C  5   1.3898  4 118.4788  3   -3.0201
# S  6   1.7465  5 118.8810  4  176.9933
# C  7   1.7497  6 103.1779  5   74.6801
# C  8   1.4215  7 128.2225  6   48.8667
# C  9   1.5373  8 119.6962  7   -2.2336
# O 10   1.1998  9 122.7413  8   81.1035
# O  7   1.4245  6 106.4017  5  -38.8122
# O 10   3.4864  9  72.5418  8   25.4569
# N 10   1.3508  9 112.3805  8 -107.5567
# H  8   1.0854  7 103.6846  6  -71.8987
# H  9   1.0965  8 108.2016  7 -120.1102
# H  5   1.0809  4 122.4201  3  176.9501
# H  4   1.0751  3 121.5289  2  177.5960
# H  3   1.0834  2 119.4871  1 -179.3770
# H  2   1.0660  1 119.4576  3 -179.8877
# H  1   1.0766  2 122.3884  3  177.4101
# H 14   1.0057 13 109.6920 12  131.3817
# H 14   1.0048 22 116.7792 13  113.3685
# C  8   3.1451  7 118.3616  9  131.5411
# C 24   1.5207  8  86.1233  7  -92.2703
# C 25   1.5467 24  94.5056  8  -63.7068
# C 26   1.4443 25 100.1371 24  -48.5602
# C 27   1.3944 26 105.4206 25   43.0025
# H 24   1.0929 25 119.7568 26  174.9601
# H 26   1.0919 25 117.8925 24  179.0337
# H 28   1.0990 27 127.8184 26  175.0305
# H 27   1.0754 26 124.8731 25 -150.8246
# H 25   1.0816 24 111.1788 26 -119.2753
# H 25   1.1009 33 109.9263 24  128.2950""", 'zmat', energy_evaluator='rdkit')

        import warnings
        warnings.filterwarnings("ignore", category=RuntimeWarning) # new mac annoyance

        z1 = Molecule.from_string(
            'C(=O)O', 'smi',
            energy_evaluator='aimnet2',
            conf_get_options={'random_seed': 12321}
        ).optimize(max_iterations=80)
        z1.internals = z1.get_bond_zmatrix()

        disp_coords = z1.get_displaced_coordinates(
            [[np.pi/6]],
            which=[-1],
            use_internals='reembed',
            strip_embedding=True
        )[0]
        z1 = z1.modify(coords=disp_coords)


        # int_scan = z1.get_internals(
        #     z1.get_scan_coordinates(
        #         [[0, np.pi, 10]],
        #         which=[-1],
        #         internals='reembed'
        #     )
        # )[:, -1]
        #
        # raise Exception(int_scan)

        # woof2 = scan_ahh.modify(internals={'specs': coordops.extract_zmatrix_internals(zm0)})
        # scan_ahh.animate_coordinate(1, .05, coordinate_expansion=bbb)

        opt = z1.optimize(max_iterations=25)#, mode='scipy', method='nelder-mead')

        opt2 = z1.modify(internals=None).optimize(max_iterations=15)
        print(opt.calculate_energy() - z1.calculate_energy())
        print(opt.calculate_energy() - opt2.calculate_energy())

        print(
            coordops.set_zmatrix_embedding(
                np.round(z1.internal_coordinates - z1.modify(coords=opt.coords).internal_coordinates, 2)
            )
        )
        print(
            coordops.set_zmatrix_embedding(
                np.round(z1.internal_coordinates - z1.modify(coords=opt2.coords).internal_coordinates, 2)
            )
        )

    @validationTest
    def test_RDKitNumberIssues(self):
        Molecule.from_string("""CC(N1)=NC2=C1C(/C=C/C3=[N+](CCC[S-](=O)(=O)=O)C4=CC=CC=C4S3)=CC=C2""", 'smi')

    @validationTest
    def test_FragEmbedding(self):
        bits = Molecule.from_string(
            'FC1CCC(C(=C)C)CC1.C1OCC(C(OC)O)C1C(OC)O',
            # 'c1ccccn1',
            'smi',
            confgen_opts=dict(verbose=True, random_seed=12321)
        )

        s3 = bits.modify(internals={'specs': [
            {"transrot": bits.fragment_indices[1]}
        ]})

        s3.calculate_energy(order=1)

        # bits.modify(energy_evaluator='aimnet2').get_embedded_molecule()

        return

        bits.to_file("/Users/Mark/Desktop/struct.mol")

        ugh = bits.modify(
            internals=bits.get_canonical_zmatrix()
        ).to_string("zmat")

        # print(ugh)

        new_bits = Molecule.from_string(ugh)
        print(np.linalg.norm(bits.fragments[0].center_of_mass - bits.fragments[1].center_of_mass))
        print(np.linalg.norm(new_bits.fragments[0].center_of_mass - new_bits.fragments[1].center_of_mass))
        # return

        enc_smi = bits.to_string('smi', remove_hydrogens=True, include_tag=True)
        print(enc_smi)
        uuuh = Molecule.from_string(
            enc_smi,
            'smi',
            add_implicit_hydrogens=True,
            confgen_opts=dict(
                verbose=True,
                ignore_smoothing_failures=True
            )
        )
        print(enc_smi)

        print(np.linalg.norm(bits.fragments[0].center_of_mass - bits.fragments[1].center_of_mass))
        print(np.linalg.norm(uuuh.fragments[0].center_of_mass - uuuh.fragments[1].center_of_mass))

        # print(nput.distance_matrix(uuuh.coords, return_triu=True))
        # print(nput.distance_matrix(bits.coords, return_triu=True))
        # print(bits.center_of_mass, uuuh.center_of_mass)

        # print(np.linalg.norm(bits.fragments[0].center_of_mass - bits.fragments[1].center_of_mass))
        # print(np.linalg.norm(uuuh.fragments[0].center_of_mass - uuuh.fragments[1].center_of_mass))

    @validationTest
    def test_PlotlyBackend(self):

        from Psience.Molecools import Molecule
        tests = [
            [[2.38332071e+00, -4.96430360e-01, 5.55111512e-17],
             [6.58769136e-03, 9.88766830e-01, 3.33066907e-16],
             [-2.38990841e+00, -4.92336470e-01, -4.99600361e-16],
             [-2.08988725e+00, -3.00745588e+00, 7.41033629e-01],
             [-5.23349285e-01, -4.08233895e+00, -1.09291015e+00],
             [2.14118284e+00, -3.26708995e+00, -4.17890874e-01],
             [3.66973327e+00, 2.55947571e-01, -1.50495176e+00],
             [3.43121803e+00, -1.92532829e-01, 1.82150928e+00],
             [-6.08841279e-02, 2.25686365e+00, 1.71104297e+00],
             [-1.14242408e-01, 2.33347006e+00, -1.62987350e+00],
             [-3.20094923e+00, -6.04268579e-01, -1.95106407e+00],
             [-3.84381689e+00, 3.36703601e-01, 1.30762207e+00],
             [-8.78968145e-01, -3.32719027e+00, -3.00837486e+00],
             [-6.11932834e-01, -6.15747003e+00, -1.09757998e+00],
             [3.41377283e+00, -4.02279921e+00, -1.90838061e+00],
             [2.60551826e+00, -4.35184035e+00, 1.33902027e+00]],
            [[2.26805368e+00, 5.92936157e-01, 0.00000000e+00],
             [2.57086576e-03, -1.18386057e+00, -2.22044605e-16],
             [-2.27062454e+00, 5.90924417e-01, 1.11022302e-16],
             [-2.30791921e+00, 1.88012754e+00, -2.29124159e+00],
             [-2.54755984e-01, 3.46313899e+00, -2.50768436e+00],
             [2.15182481e+00, 1.90602459e+00, -2.56714591e+00],
             [4.06416242e+00, -3.86408559e-01, 3.13668950e-01],
             [1.84673260e+00, 2.01434259e+00, 1.49524714e+00],
             [-2.23771964e-02, -2.45799578e+00, 1.63765593e+00],
             [-8.37319822e-02, -2.19336565e+00, -1.85026755e+00],
             [-3.97754998e+00, -5.75479502e-01, 2.12293297e-01],
             [-2.05680660e+00, 1.91328775e+00, 1.59105223e+00],
             [-3.49214054e-01, 4.52908630e+00, -4.28173613e+00],
             [-1.30106083e-01, 4.71540063e+00, -8.38457033e-01],
             [1.91592871e+00, 3.49950245e-01, -3.95763078e+00],
             [3.83541515e+00, 3.03601575e+00, -2.92206855e+00]],
            [[-2.37496538e+00, -5.16913078e-01, -2.22044605e-16],
             [-1.51564096e-03, 1.03283778e+00, 4.44089210e-16],
             [2.37648102e+00, -5.15924698e-01, 0.00000000e+00],
             [2.36112435e+00, -2.49907484e+00, 1.67326377e+00],
             [2.67785361e-01, -4.03453305e+00, 1.54551068e+00],
             [-2.20303996e+00, -2.64979149e+00, 1.86935105e+00],
             [-3.92449341e+00, 7.52133206e-01, 6.81237472e-01],
             [-2.88710909e+00, -1.17905511e+00, -1.90284281e+00],
             [-2.53225012e-03, 2.11878872e+00, -1.81438423e+00],
             [3.79870511e-02, 2.45152654e+00, 1.54550417e+00],
             [3.91483037e+00, 7.97233384e-01, 6.25195636e-01],
             [2.89487215e+00, -1.07093185e+00, -1.97610252e+00],
             [1.33465406e-01, -5.19598512e+00, -1.68889241e-01],
             [3.92752782e-01, -5.37624976e+00, 3.17017578e+00],
             [-2.57208138e+00, -2.04623188e+00, 3.82601529e+00],
             [-3.70060088e+00, -4.03037183e+00, 1.36342150e+00]],

            ## 2.196092835668157 0.9490992792321991 [ 9.34821571e-17 -1.05211987e-15 -9.49099279e-01]
            [[-2.25547156e+00, 6.00809059e-01, -8.88178420e-16],
             [2.22534122e-04, -1.20179599e+00, 1.33226763e-15],
             [2.25524902e+00, 6.00986929e-01, -4.44089210e-16],
             [1.75066984e+00, 2.75425742e+00, -1.34530337e+00],
             [1.83205907e-01, 2.48701259e+00, -3.38252635e+00],
             [-2.40251269e+00, 1.42533820e+00, -2.75825525e+00],
             [-3.98808226e+00, -3.09181461e-01, 6.73014560e-01],
             [-1.82228951e+00, 2.31487316e+00, 1.16334253e+00],
             [-4.07583757e-02, -2.46686804e+00, 1.63773950e+00],
             [4.28858789e-02, -2.32777829e+00, -1.78495144e+00],
             [4.01707996e+00, -3.90110707e-01, -5.51681364e-01],
             [2.41274739e+00, 1.21229589e+00, 2.03414022e+00],
             [-1.64379751e-01, 4.40744415e+00, -4.13914932e+00],
             [1.14647365e+00, 1.44677024e+00, -4.91357472e+00],
             [-3.96825875e+00, 2.74648401e+00, -3.11420515e+00],
             [-2.73132791e+00, -3.93039213e-01, -3.79281507e+00]]
        ]
        uuh2 = Molecule(
            ['C', 'C', 'C', 'O', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'],
            tests[3]
        )
        # uuh2 = Molecule(
        #     ('C', 'C', 'C', 'O', 'C', 'C', 'H', 'H',
        #      'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'),
        #     [[-2.37496538e+00, -5.16913078e-01, -2.22044605e-16],
        #      [-1.51564096e-03, 1.03283778e+00, 4.44089210e-16],
        #      [2.37648102e+00, -5.15924698e-01, 0.00000000e+00],
        #      [2.36112435e+00, -2.49907484e+00, 1.67326377e+00],
        #      [2.67785361e-01, -4.03453305e+00, 1.54551068e+00],
        #      [-2.20303996e+00, -2.64979149e+00, 1.86935105e+00],
        #      [-3.92449341e+00, 7.52133206e-01, 6.81237472e-01],
        #      [-2.88710909e+00, -1.17905511e+00, -1.90284281e+00],
        #      [-2.53225012e-03, 2.11878872e+00, -1.81438423e+00],
        #      [3.79870511e-02, 2.45152654e+00, 1.54550417e+00],
        #      [3.91483037e+00, 7.97233384e-01, 6.25195636e-01],
        #      [2.89487215e+00, -1.07093185e+00, -1.97610252e+00],
        #      [1.33465406e-01, -5.19598512e+00, -1.68889241e-01],
        #      [3.92752782e-01, -5.37624976e+00, 3.17017578e+00],
        #      [-2.57208138e+00, -2.04623188e+00, 3.82601529e+00],
        #      [-3.70060088e+00, -4.03037183e+00, 1.36342150e+00]]
        # ).get_embedded_molecule()
        # uuh2 = Molecule.from_string("C1CCOCC1=O", "smi").get_embedded_molecule(sel=[0, 1, 2])
        # print(uuh2.bonds)
        # uuh2 = Molecule(
        #     ('C', 'C', 'C', 'O', 'C', 'C', 'O', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'),
        #     [[2.40248768e+00, 4.87574030e-01, 1.11022302e-16],
        #      [-3.73056809e-02, -9.98589102e-01, 2.22044605e-16],
        #      [-2.36518200e+00, 5.11015071e-01, 3.33066907e-16],
        #      [-2.33820339e+00, 2.97459390e+00, 6.19461561e-01],
        #      [-2.97642862e-01, 4.08330967e+00, 1.67448987e+00],
        #      [1.93759442e+00, 2.44053814e+00, 2.03326492e+00],
        #      [3.36587018e+00, 2.59335666e+00, 3.83564101e+00],
        #      [2.51685464e+00, 1.55960350e+00, -1.78585388e+00],
        #      [3.99395479e+00, -7.47614214e-01, 4.05268307e-01],
        #      [1.86910428e-02, -2.15188852e+00, 1.76521607e+00],
        #      [9.52341488e-02, -2.39849452e+00, -1.57135667e+00],
        #      [-3.20701964e+00, 3.66560107e-01, -1.95258316e+00],
        #      [-3.84198611e+00, -4.27302551e-01, 1.19337845e+00],
        #      [-7.10626249e-01, 4.97676320e+00, 3.52639815e+00],
        #      [3.37684581e-01, 5.64952636e+00, 4.10428772e-01]]
        #     ##
        #     # [[-0.98978292, 2.13944996, -0.74785071],
        #     #  [-2.37421443, -0.17229396, 0.13897264],
        #     #  [-0.84146122, -2.57197311, 0.02119341],
        #     #  [1.54231381, -2.41614628, 1.01377972],
        #     #  [2.99186019, -0.39568311, 0.32402585],
        #     #  [1.58760012, 2.05609658, 0.39871331],
        #     #  [2.47689799, 3.92228997, 1.34505311],
        #     #  [-1.99443943, 3.78284084, 0.00833893],
        #     #  [-0.78881228, 2.09402206, -2.83188136],
        #     #  [-4.00185885, -0.44419277, -1.201296],
        #     #  [-3.20938307, 0.10260383, 2.02073282],
        #     #  [-1.91966954, -3.9299204, 1.24531937],
        #     #  [-0.92135447, -3.41409647, -1.89883835],
        #     #  [4.56567818, -0.2361968, 1.7029855],
        #     #  [3.87662593, -0.51680037, -1.53924822]]
        #     ## -2.9405216419521696 [-0.3256606  -0.16229313 -0.93145376]
        #     # [[-1.8048029, -1.6340269, 0.20826272],
        #     # [-2.11990191, 1.18155187, -0.06155707],
        #     # [0.25521018, 2.62375034, 0.5352225],
        #     # [2.40426119, 1.73354464, -0.5846305],
        #     # [2.89836166, -0.81862443, -0.28512281],
        #     # [0.6907083, -2.43287548, -0.77753825],
        #     # [0.95061245, -4.38814714, -1.98079781],
        #     # [-3.40719786, -2.51631386, -0.76390484],
        #     # [-1.87647756, -2.18293024, 2.22453257],
        #     # [-3.61347563, 1.74187993, 1.2523594],
        #     # [-2.83159418, 1.61058188, -2.0158157],
        #     # [-0.04698966, 4.56279462, -0.30967277],
        #     # [0.33894847, 2.93750686, 2.59499621],
        #     # [3.74411763, -1.14405503, 1.62293356],
        #     # [4.41821982, -1.27463704, -1.65926721]]
        #     ## -1.171263910156821 [-0.27779964  0.06078523  0.95871399]
        #     ,
        #     bonds=[[0, 1, 1.0], [1, 2, 1.0], [2, 3, 1.0], [3, 4, 1.0], [4, 5, 1.0], [5, 6, 2.0], [5, 0, 1.0], [0, 7, 1.0], [0, 8, 1.0], [1, 9, 1.0], [1, 10, 1.0], [2, 11, 1.0], [2, 12, 1.0], [4, 13, 1.0], [4, 14, 1.0]]
        # )#.get_embedded_molecule()


        view_vector, right_vector, up_vector = nput.view_matrix(
            [1, 0, 0],
            [0, 1, 0],
            output_order=['x', 'y', 'z']
        ).T
        fig = uuh2.plot(backend='x3d',
                        image_size=[500, 500],
                        highlight_atoms=[0, 1, 2],
                        draw_coords={
                            (0, 4): {
                                'label': "r",
                                'line_color':'pink',
                                'label_style': {'font_size': 22, 'color': 'red'}
                            },
                            (0, 1, 2): {
                                'label': "O",
                                'line_color':'pink',
                                'label_style': {'font_size': 40, 'color': 'blue'}
                            }
                        },
                        view_settings={
                            'up_vector': up_vector,
                            'view_vector': view_vector,
                            # 'right_vector': uuh2.coords[0] - uuh2.coords[1],
                            'view_distance': 9,
                            # 'view_center': 2.5 * right_vector
                        },
                        # plot_range=[[-1, 1], [-1, 1], [-1, 1]],
                        # font_family='Helvetica',
                        # postdraw=[
                        #     {
                        #         'pattern': 'labels_1',
                        #         'classes': True,
                        #         'replacement':{'text':"θ",  'mode':'svg',
                        #                        # 'font_options':{
                        #                        #     'family':"Times"
                        #                        # }
                        #                        }
                        #     }
                        # ]
                        )
        # uuh2 = Molecule.from_string("C1CCOCC1=O", "smi").get_embedded_molecule(sel=[0, 1, 2])
        # uuh2.plot(backend='rdkit',
        #           figure=fig,
        #           image_size=[500, 500],
        #           highlight_atoms=[0, 1, 2],
        #           draw_coords={
        #               (0, 4): {
        #                   'label': "r",
        #                   'label_style': {'font_size': 22, 'color': 'red'}
        #               },
        #               (0, 1, 2): {
        #                   'label': "\u03B8",
        #                   'label_style': {'font_size': 40, 'font_family': 'Symbol', 'color': 'blue'}
        #               }
        #           },
        #           view_settings={
        #               'up_vector': up_vector,
        #               'view_vector': view_vector,
        #               # 'right_vector': uuh2.coords[0] - uuh2.coords[1],
        #               'view_distance': 14,
        #               'view_center': 3*up_vector
        #           },
        #           font_family='Helvetica',
        #           no_free_type=False,
        #           postdraw='annotate'
        #           )
        # uuh2.plot(backend='rdkit',
        #           figure=fig,
        #           image_size=[500, 500],
        #           highlight_atoms=[0, 1, 2],
        #           draw_coords={
        #               (0, 4): {
        #                   'label': "r",
        #                   'label_style': {'font_size': 22, 'color': 'red'}
        #               },
        #               (0, 1, 2): {
        #                   'label': "\u03B8",
        #                   'label_style': {'font_size': 40, 'font_family': 'Symbol', 'color': 'blue'}
        #               }
        #           },
        #           view_settings={
        #               'up_vector': nput.rotation_matrix(view_vector, np.deg2rad(45)) @ up_vector,
        #               'view_vector': view_vector,
        #               # 'right_vector': uuh2.coords[0] - uuh2.coords[1],
        #               'view_distance': 14,
        #               'view_center': -3 * up_vector
        #           },
        #           font_family='Helvetica',
        #           no_free_type=False,
        #           postdraw='annotate'
        #           )
        fig.show()
        return

        # uuh2 = Molecule.from_string("C1CCOCC1", "smi")
        uuh3 = Molecule(
            ['C', 'C', 'C', 'O', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'],
            # [[2.38332071e+00, -4.96430360e-01, 5.55111512e-17],
            #  [6.58769136e-03, 9.88766830e-01, 3.33066907e-16],
            #  [-2.38990841e+00, -4.92336470e-01, -4.99600361e-16],
            #  [-2.08988725e+00, -3.00745588e+00, 7.41033629e-01],
            #  [-5.23349285e-01, -4.08233895e+00, -1.09291015e+00],
            #  [2.14118284e+00, -3.26708995e+00, -4.17890874e-01],
            #  [3.66973327e+00, 2.55947571e-01, -1.50495176e+00],
            #  [3.43121803e+00, -1.92532829e-01, 1.82150928e+00],
            #  [-6.08841279e-02, 2.25686365e+00, 1.71104297e+00],
            #  [-1.14242408e-01, 2.33347006e+00, -1.62987350e+00],
            #  [-3.20094923e+00, -6.04268579e-01, -1.95106407e+00],
            #  [-3.84381689e+00, 3.36703601e-01, 1.30762207e+00],
            #  [-8.78968145e-01, -3.32719027e+00, -3.00837486e+00],
            #  [-6.11932834e-01, -6.15747003e+00, -1.09757998e+00],
            #  [3.41377283e+00, -4.02279921e+00, -1.90838061e+00],
            #  [2.60551826e+00, -4.35184035e+00, 1.33902027e+00]]
            ## u/noflip (array([ 0.,  0., -1.]), array([-0.848,  0.53 ,  0.   ]), array([0.851, 0.526, 0.   ]), -0.897)
            # [[2.26805368e+00, 5.92936157e-01, 0.00000000e+00],
            #   [2.57086576e-03, -1.18386057e+00, -2.22044605e-16],
            #   [-2.27062454e+00, 5.90924417e-01, 1.11022302e-16],
            #   [-2.30791921e+00, 1.88012754e+00, -2.29124159e+00],
            #   [-2.54755984e-01, 3.46313899e+00, -2.50768436e+00],
            #   [2.15182481e+00, 1.90602459e+00, -2.56714591e+00],
            #   [4.06416242e+00, -3.86408559e-01, 3.13668950e-01],
            #   [1.84673260e+00, 2.01434259e+00, 1.49524714e+00],
            #   [-2.23771964e-02, -2.45799578e+00, 1.63765593e+00],
            #   [-8.37319822e-02, -2.19336565e+00, -1.85026755e+00],
            #   [-3.97754998e+00, -5.75479502e-01, 2.12293297e-01],
            #   [-2.05680660e+00, 1.91328775e+00, 1.59105223e+00],
            #   [-3.49214054e-01, 4.52908630e+00, -4.28173613e+00],
            #   [-1.30106083e-01, 4.71540063e+00, -8.38457033e-01],
            #   [1.91592871e+00, 3.49950245e-01, -3.95763078e+00],
            #   [3.83541515e+00, 3.03601575e+00, -2.92206855e+00]]
            ## v/no flip: (array([ 0., -0.,  1.]), array([0.787, 0.617, 0.   ]), array([-0.788,  0.615,  0.   ]), 0.971)
            # [[-2.37496538e+00, -5.16913078e-01, -2.22044605e-16],
            # [-1.51564096e-03, 1.03283778e+00, 4.44089210e-16],
            # [2.37648102e+00, -5.15924698e-01, 0.00000000e+00],
            # [2.36112435e+00, -2.49907484e+00, 1.67326377e+00],
            # [2.67785361e-01, -4.03453305e+00, 1.54551068e+00],
            # [-2.20303996e+00, -2.64979149e+00, 1.86935105e+00],
            # [-3.92449341e+00, 7.52133206e-01, 6.81237472e-01],
            # [-2.88710909e+00, -1.17905511e+00, -1.90284281e+00],
            # [-2.53225012e-03, 2.11878872e+00, -1.81438423e+00],
            # [3.79870511e-02, 2.45152654e+00, 1.54550417e+00],
            # [3.91483037e+00, 7.97233384e-01, 6.25195636e-01],
            # [2.89487215e+00, -1.07093185e+00, -1.97610252e+00],
            # [1.33465406e-01, -5.19598512e+00, -1.68889241e-01],
            # [3.92752782e-01, -5.37624976e+00, 3.17017578e+00],
            # [-2.57208138e+00, -2.04623188e+00, 3.82601529e+00],
            # [-3.70060088e+00, -4.03037183e+00, 1.36342150e+00]]
            ## v/flip: (array([-0., -0.,  1.]), array([-0.837, -0.547, -0.   ]), array([ 0.838, -0.546, -0.   ]), 0.915)
            [[-2.25547156e+00, 6.00809059e-01, -8.88178420e-16],
             [2.22534122e-04, -1.20179599e+00, 1.33226763e-15],
             [2.25524902e+00, 6.00986929e-01, -4.44089210e-16],
             [1.75066984e+00, 2.75425742e+00, -1.34530337e+00],
             [1.83205907e-01, 2.48701259e+00, -3.38252635e+00],
             [-2.40251269e+00, 1.42533820e+00, -2.75825525e+00],
             [-3.98808226e+00, -3.09181461e-01, 6.73014560e-01],
             [-1.82228951e+00, 2.31487316e+00, 1.16334253e+00],
             [-4.07583757e-02, -2.46686804e+00, 1.63773950e+00],
             [4.28858789e-02, -2.32777829e+00, -1.78495144e+00],
             [4.01707996e+00, -3.90110707e-01, -5.51681364e-01],
             [2.41274739e+00, 1.21229589e+00, 2.03414022e+00],
             [-1.64379751e-01, 4.40744415e+00, -4.13914932e+00],
             [1.14647365e+00, 1.44677024e+00, -4.91357472e+00],
             [-3.96825875e+00, 2.74648401e+00, -3.11420515e+00],
             [-2.73132791e+00, -3.93039213e-01, -3.79281507e+00]]
            ## u/flip: (array([ 0., -0., -1.]), array([ 0.781, -0.624,  0.   ]), array([-0.781, -0.624,  0.   ]), -0.975)
        )
        ploot = uuh3.plot(backend='matplotlib3D',
                          # include_save_buttons=True,
                          # include_script_interface=True,
                          highlight_atoms=[0, 1, 2],
                          draw_coords={
                              (0, 4): {
                                  'label': 'r',
                                  'label_style': {'color': 'blue', 'font_size': 20}
                              },
                              (0, 1, 2): {
                                  'label': "q",
                                  'label_style': {'color': 'red', 'font_size':32, 'font_family':'serif'}
                              }
                          },
                          view_settings={
                              'view_vector': [0, 0, 1],
                              # 'view_distance': 5
                          })
        # print(ploot.tostring())
        ploot.write("/Users/Mark/Desktop/why.html")
        # ploot.show()

        # bits = Molecule.from_file(TestManager.test_data("water_freq.fchk"))
        # bits.plot(
        #     backend=(
        #         # 'matplotlib3D'
        #         # 'x3d'
        #         'plotly'
        #     ),
        #     # principle_axes=True,
        #     image_size=[500, 500],
        #     # theme='flat',
        #     # draw_coords={
        #     #     (0, 1, 2): {'label': r'$\theta$',
        #     #                 'label_style': {
        #     #                     'offset_magnitude': 1.1,
        #     #                     'billboard': True,
        #     #                     'color': 'red', 'fontsize': 24
        #     #                 }}
        #     # },
        #     # highlight_atoms=[0, 1, 2],
        #     # view_settings={
        #     #     # 'view_distance':200,
        #     #     'view_vector': [0, 0, 1],
        #     # },
        #     include_save_buttons=True
        # ).show()  # .savefig("/Users/Mark/Desktop/view_xy_simp_bonds.svg")
        # # return

    @validationTest
    def test_SVGBackend(self):
        from Psience.Molecools import Molecule
        tests = [
            [[2.38332071e+00, -4.96430360e-01, 5.55111512e-17],
             [6.58769136e-03, 9.88766830e-01, 3.33066907e-16],
             [-2.38990841e+00, -4.92336470e-01, -4.99600361e-16],
             [-2.08988725e+00, -3.00745588e+00, 7.41033629e-01],
             [-5.23349285e-01, -4.08233895e+00, -1.09291015e+00],
             [2.14118284e+00, -3.26708995e+00, -4.17890874e-01],
             [3.66973327e+00, 2.55947571e-01, -1.50495176e+00],
             [3.43121803e+00, -1.92532829e-01, 1.82150928e+00],
             [-6.08841279e-02, 2.25686365e+00, 1.71104297e+00],
             [-1.14242408e-01, 2.33347006e+00, -1.62987350e+00],
             [-3.20094923e+00, -6.04268579e-01, -1.95106407e+00],
             [-3.84381689e+00, 3.36703601e-01, 1.30762207e+00],
             [-8.78968145e-01, -3.32719027e+00, -3.00837486e+00],
             [-6.11932834e-01, -6.15747003e+00, -1.09757998e+00],
             [3.41377283e+00, -4.02279921e+00, -1.90838061e+00],
             [2.60551826e+00, -4.35184035e+00, 1.33902027e+00]],
            [[2.26805368e+00, 5.92936157e-01, 0.00000000e+00],
             [2.57086576e-03, -1.18386057e+00, -2.22044605e-16],
             [-2.27062454e+00, 5.90924417e-01, 1.11022302e-16],
             [-2.30791921e+00, 1.88012754e+00, -2.29124159e+00],
             [-2.54755984e-01, 3.46313899e+00, -2.50768436e+00],
             [2.15182481e+00, 1.90602459e+00, -2.56714591e+00],
             [4.06416242e+00, -3.86408559e-01, 3.13668950e-01],
             [1.84673260e+00, 2.01434259e+00, 1.49524714e+00],
             [-2.23771964e-02, -2.45799578e+00, 1.63765593e+00],
             [-8.37319822e-02, -2.19336565e+00, -1.85026755e+00],
             [-3.97754998e+00, -5.75479502e-01, 2.12293297e-01],
             [-2.05680660e+00, 1.91328775e+00, 1.59105223e+00],
             [-3.49214054e-01, 4.52908630e+00, -4.28173613e+00],
             [-1.30106083e-01, 4.71540063e+00, -8.38457033e-01],
             [1.91592871e+00, 3.49950245e-01, -3.95763078e+00],
             [3.83541515e+00, 3.03601575e+00, -2.92206855e+00]],
            [[-2.37496538e+00, -5.16913078e-01, -2.22044605e-16],
             [-1.51564096e-03, 1.03283778e+00, 4.44089210e-16],
             [2.37648102e+00, -5.15924698e-01, 0.00000000e+00],
             [2.36112435e+00, -2.49907484e+00, 1.67326377e+00],
             [2.67785361e-01, -4.03453305e+00, 1.54551068e+00],
             [-2.20303996e+00, -2.64979149e+00, 1.86935105e+00],
             [-3.92449341e+00, 7.52133206e-01, 6.81237472e-01],
             [-2.88710909e+00, -1.17905511e+00, -1.90284281e+00],
             [-2.53225012e-03, 2.11878872e+00, -1.81438423e+00],
             [3.79870511e-02, 2.45152654e+00, 1.54550417e+00],
             [3.91483037e+00, 7.97233384e-01, 6.25195636e-01],
             [2.89487215e+00, -1.07093185e+00, -1.97610252e+00],
             [1.33465406e-01, -5.19598512e+00, -1.68889241e-01],
             [3.92752782e-01, -5.37624976e+00, 3.17017578e+00],
             [-2.57208138e+00, -2.04623188e+00, 3.82601529e+00],
             [-3.70060088e+00, -4.03037183e+00, 1.36342150e+00]],

            ## 2.196092835668157 0.9490992792321991 [ 9.34821571e-17 -1.05211987e-15 -9.49099279e-01]
            [[-2.25547156e+00, 6.00809059e-01, -8.88178420e-16],
             [2.22534122e-04, -1.20179599e+00, 1.33226763e-15],
             [2.25524902e+00, 6.00986929e-01, -4.44089210e-16],
             [1.75066984e+00, 2.75425742e+00, -1.34530337e+00],
             [1.83205907e-01, 2.48701259e+00, -3.38252635e+00],
             [-2.40251269e+00, 1.42533820e+00, -2.75825525e+00],
             [-3.98808226e+00, -3.09181461e-01, 6.73014560e-01],
             [-1.82228951e+00, 2.31487316e+00, 1.16334253e+00],
             [-4.07583757e-02, -2.46686804e+00, 1.63773950e+00],
             [4.28858789e-02, -2.32777829e+00, -1.78495144e+00],
             [4.01707996e+00, -3.90110707e-01, -5.51681364e-01],
             [2.41274739e+00, 1.21229589e+00, 2.03414022e+00],
             [-1.64379751e-01, 4.40744415e+00, -4.13914932e+00],
             [1.14647365e+00, 1.44677024e+00, -4.91357472e+00],
             [-3.96825875e+00, 2.74648401e+00, -3.11420515e+00],
             [-2.73132791e+00, -3.93039213e-01, -3.79281507e+00]]
        ]
        uuh2 = Molecule(
            ['C', 'C', 'C', 'O', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'],
            tests[3]
        )
        # uuh2 = Molecule(
        #     ('C', 'C', 'C', 'O', 'C', 'C', 'H', 'H',
        #      'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'),
        #     [[-2.37496538e+00, -5.16913078e-01, -2.22044605e-16],
        #      [-1.51564096e-03, 1.03283778e+00, 4.44089210e-16],
        #      [2.37648102e+00, -5.15924698e-01, 0.00000000e+00],
        #      [2.36112435e+00, -2.49907484e+00, 1.67326377e+00],
        #      [2.67785361e-01, -4.03453305e+00, 1.54551068e+00],
        #      [-2.20303996e+00, -2.64979149e+00, 1.86935105e+00],
        #      [-3.92449341e+00, 7.52133206e-01, 6.81237472e-01],
        #      [-2.88710909e+00, -1.17905511e+00, -1.90284281e+00],
        #      [-2.53225012e-03, 2.11878872e+00, -1.81438423e+00],
        #      [3.79870511e-02, 2.45152654e+00, 1.54550417e+00],
        #      [3.91483037e+00, 7.97233384e-01, 6.25195636e-01],
        #      [2.89487215e+00, -1.07093185e+00, -1.97610252e+00],
        #      [1.33465406e-01, -5.19598512e+00, -1.68889241e-01],
        #      [3.92752782e-01, -5.37624976e+00, 3.17017578e+00],
        #      [-2.57208138e+00, -2.04623188e+00, 3.82601529e+00],
        #      [-3.70060088e+00, -4.03037183e+00, 1.36342150e+00]]
        # ).get_embedded_molecule()
        # uuh2 = Molecule.from_string("C1CCOCC1=O", "smi").get_embedded_molecule(sel=[0, 1, 2])
        # print(uuh2.bonds)
        # uuh2 = Molecule(
        #     ('C', 'C', 'C', 'O', 'C', 'C', 'O', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'),
        #     [[2.40248768e+00, 4.87574030e-01, 1.11022302e-16],
        #      [-3.73056809e-02, -9.98589102e-01, 2.22044605e-16],
        #      [-2.36518200e+00, 5.11015071e-01, 3.33066907e-16],
        #      [-2.33820339e+00, 2.97459390e+00, 6.19461561e-01],
        #      [-2.97642862e-01, 4.08330967e+00, 1.67448987e+00],
        #      [1.93759442e+00, 2.44053814e+00, 2.03326492e+00],
        #      [3.36587018e+00, 2.59335666e+00, 3.83564101e+00],
        #      [2.51685464e+00, 1.55960350e+00, -1.78585388e+00],
        #      [3.99395479e+00, -7.47614214e-01, 4.05268307e-01],
        #      [1.86910428e-02, -2.15188852e+00, 1.76521607e+00],
        #      [9.52341488e-02, -2.39849452e+00, -1.57135667e+00],
        #      [-3.20701964e+00, 3.66560107e-01, -1.95258316e+00],
        #      [-3.84198611e+00, -4.27302551e-01, 1.19337845e+00],
        #      [-7.10626249e-01, 4.97676320e+00, 3.52639815e+00],
        #      [3.37684581e-01, 5.64952636e+00, 4.10428772e-01]]
        #     ##
        #     # [[-0.98978292, 2.13944996, -0.74785071],
        #     #  [-2.37421443, -0.17229396, 0.13897264],
        #     #  [-0.84146122, -2.57197311, 0.02119341],
        #     #  [1.54231381, -2.41614628, 1.01377972],
        #     #  [2.99186019, -0.39568311, 0.32402585],
        #     #  [1.58760012, 2.05609658, 0.39871331],
        #     #  [2.47689799, 3.92228997, 1.34505311],
        #     #  [-1.99443943, 3.78284084, 0.00833893],
        #     #  [-0.78881228, 2.09402206, -2.83188136],
        #     #  [-4.00185885, -0.44419277, -1.201296],
        #     #  [-3.20938307, 0.10260383, 2.02073282],
        #     #  [-1.91966954, -3.9299204, 1.24531937],
        #     #  [-0.92135447, -3.41409647, -1.89883835],
        #     #  [4.56567818, -0.2361968, 1.7029855],
        #     #  [3.87662593, -0.51680037, -1.53924822]]
        #     ## -2.9405216419521696 [-0.3256606  -0.16229313 -0.93145376]
        #     # [[-1.8048029, -1.6340269, 0.20826272],
        #     # [-2.11990191, 1.18155187, -0.06155707],
        #     # [0.25521018, 2.62375034, 0.5352225],
        #     # [2.40426119, 1.73354464, -0.5846305],
        #     # [2.89836166, -0.81862443, -0.28512281],
        #     # [0.6907083, -2.43287548, -0.77753825],
        #     # [0.95061245, -4.38814714, -1.98079781],
        #     # [-3.40719786, -2.51631386, -0.76390484],
        #     # [-1.87647756, -2.18293024, 2.22453257],
        #     # [-3.61347563, 1.74187993, 1.2523594],
        #     # [-2.83159418, 1.61058188, -2.0158157],
        #     # [-0.04698966, 4.56279462, -0.30967277],
        #     # [0.33894847, 2.93750686, 2.59499621],
        #     # [3.74411763, -1.14405503, 1.62293356],
        #     # [4.41821982, -1.27463704, -1.65926721]]
        #     ## -1.171263910156821 [-0.27779964  0.06078523  0.95871399]
        #     ,
        #     bonds=[[0, 1, 1.0], [1, 2, 1.0], [2, 3, 1.0], [3, 4, 1.0], [4, 5, 1.0], [5, 6, 2.0], [5, 0, 1.0], [0, 7, 1.0], [0, 8, 1.0], [1, 9, 1.0], [1, 10, 1.0], [2, 11, 1.0], [2, 12, 1.0], [4, 13, 1.0], [4, 14, 1.0]]
        # )#.get_embedded_molecule()

        # uuh2 = Molecule.from_file(
        #     TestManager.test_data("water_freq.fchk")
        # )

        right_vector, up_vector, view_vector = nput.view_matrix(
            [0, 1, 0],
            [0, 0, 1],
            output_order=['x', 'y', 'z']
        ).T #@ nput.rotation_matrix([1, 1, 0], np.pi/1.5)
        uuh2.coords = uuh2.coords - uuh2.center_of_mass[np.newaxis]
        fig = uuh2.plot(backend='x3d',
                        image_size=[500, 500],
                        highlight_atoms=[0, 1, 2],
                        draw_coords={
                            (0, 4): {
                                'label': "r",
                                'line_color': 'pink',
                                'label_style': {'color': 'red'}
                            },
                            (0, 1, 2): {
                                'label': "O",
                                'line_color': 'pink',
                                'label_style': {'color': 'blue'}
                            }
                        },
                        view_settings={
                            'up_vector': up_vector,
                            'view_vector': view_vector,
                            # 'viewAll':True,
                            # 'right_vector': uuh2.coords[0] - uuh2.coords[1],
                            'view_distance': 15,
                            # 'view_center': 2.5 * right_vector
                        },
                        # plot_range=np.array([[-1, 1], [-1, 1], [-1, 1]]),
                        # font_family='Helvetica',
                        # postdraw=[
                        #     {
                        #         'pattern': 'labels_1',
                        #         'classes': True,
                        #         'replacement':{'text':"θ",  'mode':'svg',
                        #                        # 'font_options':{
                        #                        #     'family':"Times"
                        #                        # }
                        #                        }
                        #     }
                        # ]
                        )
        # uuh2 = Molecule.from_string("C1CCOCC1=O", "smi").get_embedded_molecule(sel=[0, 1, 2])
        # uuh2.plot(backend='rdkit',
        #           figure=fig,
        #           image_size=[500, 500],
        #           highlight_atoms=[0, 1, 2],
        #           draw_coords={
        #               (0, 4): {
        #                   'label': "r",
        #                   'label_style': {'font_size': 22, 'color': 'red'}
        #               },
        #               (0, 1, 2): {
        #                   'label': "\u03B8",
        #                   'label_style': {'font_size': 40, 'font_family': 'Symbol', 'color': 'blue'}
        #               }
        #           },
        #           view_settings={
        #               'up_vector': up_vector,
        #               'view_vector': view_vector,
        #               # 'right_vector': uuh2.coords[0] - uuh2.coords[1],
        #               'view_distance': 14,
        #               'view_center': 3*up_vector
        #           },
        #           font_family='Helvetica',
        #           no_free_type=False,
        #           postdraw='annotate'
        #           )
        # uuh2.plot(backend='rdkit',
        #           figure=fig,
        #           image_size=[500, 500],
        #           highlight_atoms=[0, 1, 2],
        #           draw_coords={
        #               (0, 4): {
        #                   'label': "r",
        #                   'label_style': {'font_size': 22, 'color': 'red'}
        #               },
        #               (0, 1, 2): {
        #                   'label': "\u03B8",
        #                   'label_style': {'font_size': 40, 'font_family': 'Symbol', 'color': 'blue'}
        #               }
        #           },
        #           view_settings={
        #               'up_vector': nput.rotation_matrix(view_vector, np.deg2rad(45)) @ up_vector,
        #               'view_vector': view_vector,
        #               # 'right_vector': uuh2.coords[0] - uuh2.coords[1],
        #               'view_distance': 14,
        #               'view_center': -3 * up_vector
        #           },
        #           font_family='Helvetica',
        #           no_free_type=False,
        #           postdraw='annotate'
        #           )

        # print(
        #     fig.to_widget().tostring(prettify=True)
        # )
        # fig.to_widget().write(os.path.expanduser("~/Desktop/mol_as_svg.html"))
        fig.show()
        return

        # uuh2 = Molecule.from_string("C1CCOCC1", "smi")
        uuh3 = Molecule(
            ['C', 'C', 'C', 'O', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'],
            # [[2.38332071e+00, -4.96430360e-01, 5.55111512e-17],
            #  [6.58769136e-03, 9.88766830e-01, 3.33066907e-16],
            #  [-2.38990841e+00, -4.92336470e-01, -4.99600361e-16],
            #  [-2.08988725e+00, -3.00745588e+00, 7.41033629e-01],
            #  [-5.23349285e-01, -4.08233895e+00, -1.09291015e+00],
            #  [2.14118284e+00, -3.26708995e+00, -4.17890874e-01],
            #  [3.66973327e+00, 2.55947571e-01, -1.50495176e+00],
            #  [3.43121803e+00, -1.92532829e-01, 1.82150928e+00],
            #  [-6.08841279e-02, 2.25686365e+00, 1.71104297e+00],
            #  [-1.14242408e-01, 2.33347006e+00, -1.62987350e+00],
            #  [-3.20094923e+00, -6.04268579e-01, -1.95106407e+00],
            #  [-3.84381689e+00, 3.36703601e-01, 1.30762207e+00],
            #  [-8.78968145e-01, -3.32719027e+00, -3.00837486e+00],
            #  [-6.11932834e-01, -6.15747003e+00, -1.09757998e+00],
            #  [3.41377283e+00, -4.02279921e+00, -1.90838061e+00],
            #  [2.60551826e+00, -4.35184035e+00, 1.33902027e+00]]
            ## u/noflip (array([ 0.,  0., -1.]), array([-0.848,  0.53 ,  0.   ]), array([0.851, 0.526, 0.   ]), -0.897)
            # [[2.26805368e+00, 5.92936157e-01, 0.00000000e+00],
            #   [2.57086576e-03, -1.18386057e+00, -2.22044605e-16],
            #   [-2.27062454e+00, 5.90924417e-01, 1.11022302e-16],
            #   [-2.30791921e+00, 1.88012754e+00, -2.29124159e+00],
            #   [-2.54755984e-01, 3.46313899e+00, -2.50768436e+00],
            #   [2.15182481e+00, 1.90602459e+00, -2.56714591e+00],
            #   [4.06416242e+00, -3.86408559e-01, 3.13668950e-01],
            #   [1.84673260e+00, 2.01434259e+00, 1.49524714e+00],
            #   [-2.23771964e-02, -2.45799578e+00, 1.63765593e+00],
            #   [-8.37319822e-02, -2.19336565e+00, -1.85026755e+00],
            #   [-3.97754998e+00, -5.75479502e-01, 2.12293297e-01],
            #   [-2.05680660e+00, 1.91328775e+00, 1.59105223e+00],
            #   [-3.49214054e-01, 4.52908630e+00, -4.28173613e+00],
            #   [-1.30106083e-01, 4.71540063e+00, -8.38457033e-01],
            #   [1.91592871e+00, 3.49950245e-01, -3.95763078e+00],
            #   [3.83541515e+00, 3.03601575e+00, -2.92206855e+00]]
            ## v/no flip: (array([ 0., -0.,  1.]), array([0.787, 0.617, 0.   ]), array([-0.788,  0.615,  0.   ]), 0.971)
            # [[-2.37496538e+00, -5.16913078e-01, -2.22044605e-16],
            # [-1.51564096e-03, 1.03283778e+00, 4.44089210e-16],
            # [2.37648102e+00, -5.15924698e-01, 0.00000000e+00],
            # [2.36112435e+00, -2.49907484e+00, 1.67326377e+00],
            # [2.67785361e-01, -4.03453305e+00, 1.54551068e+00],
            # [-2.20303996e+00, -2.64979149e+00, 1.86935105e+00],
            # [-3.92449341e+00, 7.52133206e-01, 6.81237472e-01],
            # [-2.88710909e+00, -1.17905511e+00, -1.90284281e+00],
            # [-2.53225012e-03, 2.11878872e+00, -1.81438423e+00],
            # [3.79870511e-02, 2.45152654e+00, 1.54550417e+00],
            # [3.91483037e+00, 7.97233384e-01, 6.25195636e-01],
            # [2.89487215e+00, -1.07093185e+00, -1.97610252e+00],
            # [1.33465406e-01, -5.19598512e+00, -1.68889241e-01],
            # [3.92752782e-01, -5.37624976e+00, 3.17017578e+00],
            # [-2.57208138e+00, -2.04623188e+00, 3.82601529e+00],
            # [-3.70060088e+00, -4.03037183e+00, 1.36342150e+00]]
            ## v/flip: (array([-0., -0.,  1.]), array([-0.837, -0.547, -0.   ]), array([ 0.838, -0.546, -0.   ]), 0.915)
            [[-2.25547156e+00, 6.00809059e-01, -8.88178420e-16],
             [2.22534122e-04, -1.20179599e+00, 1.33226763e-15],
             [2.25524902e+00, 6.00986929e-01, -4.44089210e-16],
             [1.75066984e+00, 2.75425742e+00, -1.34530337e+00],
             [1.83205907e-01, 2.48701259e+00, -3.38252635e+00],
             [-2.40251269e+00, 1.42533820e+00, -2.75825525e+00],
             [-3.98808226e+00, -3.09181461e-01, 6.73014560e-01],
             [-1.82228951e+00, 2.31487316e+00, 1.16334253e+00],
             [-4.07583757e-02, -2.46686804e+00, 1.63773950e+00],
             [4.28858789e-02, -2.32777829e+00, -1.78495144e+00],
             [4.01707996e+00, -3.90110707e-01, -5.51681364e-01],
             [2.41274739e+00, 1.21229589e+00, 2.03414022e+00],
             [-1.64379751e-01, 4.40744415e+00, -4.13914932e+00],
             [1.14647365e+00, 1.44677024e+00, -4.91357472e+00],
             [-3.96825875e+00, 2.74648401e+00, -3.11420515e+00],
             [-2.73132791e+00, -3.93039213e-01, -3.79281507e+00]]
            ## u/flip: (array([ 0., -0., -1.]), array([ 0.781, -0.624,  0.   ]), array([-0.781, -0.624,  0.   ]), -0.975)
        )
        ploot = uuh3.plot(backend='matplotlib3D',
                          # include_save_buttons=True,
                          # include_script_interface=True,
                          highlight_atoms=[0, 1, 2],
                          draw_coords={
                              (0, 4): {
                                  'label': 'r',
                                  'label_style': {'color': 'blue', 'font_size': 20}
                              },
                              (0, 1, 2): {
                                  'label': "q",
                                  'label_style': {'color': 'red', 'font_size': 32, 'font_family': 'serif'}
                              }
                          },
                          view_settings={
                              'view_vector': [0, 0, 1],
                              # 'view_distance': 5
                          })
        # print(ploot.tostring())
        ploot.write("/Users/Mark/Desktop/why.html")
        # ploot.show()

        # bits = Molecule.from_file(TestManager.test_data("water_freq.fchk"))
        # bits.plot(
        #     backend=(
        #         # 'matplotlib3D'
        #         # 'x3d'
        #         'plotly'
        #     ),
        #     # principle_axes=True,
        #     image_size=[500, 500],
        #     # theme='flat',
        #     # draw_coords={
        #     #     (0, 1, 2): {'label': r'$\theta$',
        #     #                 'label_style': {
        #     #                     'offset_magnitude': 1.1,
        #     #                     'billboard': True,
        #     #                     'color': 'red', 'fontsize': 24
        #     #                 }}
        #     # },
        #     # highlight_atoms=[0, 1, 2],
        #     # view_settings={
        #     #     # 'view_distance':200,
        #     #     'view_vector': [0, 0, 1],
        #     # },
        #     include_save_buttons=True
        # ).show()  # .savefig("/Users/Mark/Desktop/view_xy_simp_bonds.svg")
        # # return

    @validationTest
    def test_BondGraphZMatrixIssues(self):
        horizontal_structure = '''
  C   -0.0531542   -0.2705697   -0.4038580
  C   -1.2672917   -0.9889447   -0.5655902
  C   -2.5115504   -0.3220336   -0.4206904
  C   -3.7326727   -1.0553190   -0.5317414
  C   -3.7004761   -2.4495575   -0.8014742
  C   -2.4508476   -3.1180666   -0.9547299
  C   -2.4145647   -4.5238734   -1.2198701
  C   -3.6621332   -5.2283272   -1.3328590
  C   -4.8569665   -4.5892997   -1.1846209
  H   -5.7886526   -5.1430158   -1.2712961
  C   -4.9257958   -3.1802597   -0.9094725
  C   -6.1422696   -2.5061569   -0.7417504
  H   -7.0716942   -3.0658565   -0.8238723
  C   -6.1964402   -1.1330449   -0.4675717
  C   -7.4373714   -0.4323877   -0.2841295
  H   -8.3634594   -0.9960396   -0.3668735
  C   -7.4675572    0.9027752   -0.0104014
  H   -8.4176460    1.4127807    0.1290235
  C   -6.2596814    1.6721378    0.1065559
  C   -6.2671014    3.0419631    0.4027137
  H   -7.2201538    3.5465627    0.5477102
  C   -5.0827029    3.7776905    0.5303886
  C   -5.0759603    5.1765715    0.8640348
  H   -6.0312166    5.6798428    0.9938500
  C   -3.9096314    5.8651178    1.0239392
  H   -3.9280356    6.9205177    1.2854065
  C   -2.6313319    5.2261265    0.8675107
  C   -2.6065567    3.8427713    0.5018790
  C   -3.8255173    3.1208176    0.3413015
  C   -3.7956779    1.7351059    0.0306407
  C   -5.0084680    1.0043192   -0.0784111
  C   -4.9769052   -0.3925329   -0.3620967
  C   -2.5425203    1.0693727   -0.1391634
  C   -1.3290713    1.7891599    0.0179298
  C   -0.0839449    1.1156505   -0.0953746
  C    1.1355172    1.8214190    0.1444458
  C    1.1010171    3.1946003    0.5049308
  C   -0.1474428    3.8757259    0.6003097
  C   -1.3601817    3.1780730    0.3542430
  C   -0.1858234    5.2499558    0.9952833
  C   -1.4242191    5.8941573    1.1061277
  H   -1.4495383    6.9350662    1.4194308
  C    1.0532032    5.9119919    1.2917556
  C    2.2465413    5.2622367    1.1941631
  C    2.3197271    3.8832904    0.7977264
  C    3.5352228    3.1943918    0.6970674
  C    3.5930945    1.8423682    0.3340199
  C    2.3778386    1.1384405    0.0625407
  C    2.4102323   -0.2468567   -0.2728220
  C    1.1993952   -0.9519230   -0.5029697
  C    1.2283278   -2.3432451   -0.7890186
  C    0.0098448   -3.0637208   -0.9542032
  C   -1.2361383   -2.3926339   -0.8313463
  C    0.0349492   -4.4697870   -1.2200387
  C   -1.1751012   -5.1636932   -1.3514989
  H   -1.1509809   -6.2327171   -1.5521573
  C    1.3123634   -5.1187284   -1.3301511
  H    1.3318089   -6.1850993   -1.5412707
  C    2.4774035   -4.4302726   -1.1665430
  C    2.4838369   -3.0225639   -0.8768600
  C    3.6676819   -2.3042610   -0.6642272
  C    3.6602742   -0.9382196   -0.3533601
  C    4.8654746   -0.1983466   -0.0980494
  H    5.8151245   -0.7232786   -0.1667677
  C    4.8334523    1.1238831    0.2322146
  H    5.7571926    1.6620172    0.4300728
  H    4.6199581   -2.8265000   -0.7288371
  H    3.4322477   -4.9442495   -1.2457993
  H    4.4607529    3.7230634    0.9146432
  H    3.1731065    5.7801519    1.4281981
  H    1.0166878    6.9498086    1.6110402
  H   -3.6334126   -6.2956952   -1.5377463
  N   -0.9768306   -0.8645727    2.6960567
  C   -1.1153261    0.4844050    3.0234710
  C   -0.0778126    1.3920491    3.2650545
  C   -0.3249244    2.7228872    3.5872077
  C   -1.6325156    3.2294720    3.7094695
  C   -2.6582688    2.3098377    3.4306428
  C   -2.4217818    0.9854865    3.0975439
  F   -3.4554562    0.1841828    2.8127935
  F   -3.9567056    2.6815385    3.4512927
  C   -1.8441701    4.6700822    4.1229961
  O   -0.9060656    5.3915137    4.4314908
  N   -3.1391561    5.1091075    4.1948416
  H   -3.8828970    4.6489540    3.6932711
  H   -3.2290700    6.1045076    4.3462867
  F    0.7511767    3.4882777    3.7556598
  F    1.1959563    0.9736924    3.1611571
  N    0.1556565   -1.3754611    2.6214998
  N    1.0867607   -2.0159052    2.5061938'''
        horp = Molecule.from_string(horizontal_structure, units='Angstroms')
        f1:Molecule = horp.fragments[1]
        zmat = f1.get_bond_zmatrix()
        hmm = [z[0] for z in zmat[:8]]
        segs = f1.find_backbone_segments()
        f1.plot(highlight_atoms=segs[0]).show()

    @validationTest
    def test_SimpleRelaxedScan(self):
        # import McUtils.Iterators as itut
        #
        # for x in itut.zigzag_product(['a', 'b', 'c'], [0, 1, 2], [-2, 1, 0]):
        #     print(x)
        # return

        import warnings
        warnings.filterwarnings("ignore", category=RuntimeWarning) # new mac annoyance

        mol = Molecule.from_string('[CH3:1][O:2][H:3]',
                                   energy_evaluator='aimnet2',
                                   internals='zmatrix'
                                   )
        # mol.plot(highlight_atoms=[1, 2]).show()


        _, scan_geoms, _ = mol.relaxed_scan(
            [0, np.pi/6, 10],
            (2, 1, 0, 3),
            max_iterations=50
        )

        _, rigid_geoms, _ = mol.relaxed_scan(
            [0, np.pi / 6, 10],
            (2, 1, 0, 3),
            max_iterations=0
        )

        import McUtils.Plots as plt

        fig = plt.Plot(
            np.linspace(0, np.pi/6, 10),
            mol.calculate_energy(coords=scan_geoms, use_internals=False)
        )
        plt.Plot(
            np.linspace(0, np.pi / 6, 10),
            mol.calculate_energy(coords=rigid_geoms, use_internals=False),
            figure=fig
        )
        fig.show()

        # mol.plot(scan_geoms, include_save_buttons=True).show()

        return


        # import rdkit.Chem
        # import numpy as np

        # mol1 = rdkit.Chem.MolFromSmiles("[CH2:1]1[CH2:3]C=C[CH2:4][CH2:2]1")
        # mol1 = rdkit.Chem.AddHs(mol1)
        # mol2 = rdkit.Chem.MolFromSmiles("[CH2:1]1[CH2:3]C=C[CH2:4][CH2:2]1")
        # mol2 = rdkit.Chem.AddHs(mol2)
        # mol1 = rdkit.Chem.RemoveHs(mol1)
        # mol2 = rdkit.Chem.RemoveHs(mol2)
        # mol = rdkit.Chem.CombineMols(mol1, mol2)
        # mol = rdkit.Chem.AddHs(mol)
        #
        # coords = np.random.rand(mol.GetNumAtoms(), 3)
        # conf = rdkit.Chem.Conformer(coords.shape[0])
        # conf.SetPositions(coords)
        # conf.SetId(0)
        # conf_id = mol.AddConformer(conf)
        # conf = mol.GetConformer(conf_id)
        # print(conf_id, conf.GetPositions())
        # return

    @validationTest
    def test_RequiredCoordinateZMatrix(self):

        import warnings
        warnings.filterwarnings("ignore", category=RuntimeWarning)  # new mac annoyance

        mol: Molecule = Molecule.from_string('[CH2:1]1[CH2:3]C=C[CH2:4][CH2:2]1',
                                             energy_evaluator='aimnet2',
                                             internals='zmatrix'
                                             )

        zm = mol.get_bond_zmatrix(required_coordinates=[
            # (0, 2),
            (1, 3),
            (0, 1),
            # (3, 2, 4)
        ])

    @validationTest
    def test_ProblemTransformedZMatrix(self):
        import McUtils.Coordinerds as coordops

        prod: Molecule = Molecule.from_string("[CH2:1]1[CH2:3]C=C[CH2:4][CH2:2]1")
        react = prod.apply_smarts("[C:1]1[C:3][C:5]=[C:6][C:4][C:2]1 >> [C:1]=[C:2].[C:3]=[C:5]-[C:6]=[C:4]")

        # print(react[0].fragment_indices)
        zzz = react[0].get_bond_zmatrix()
        print(np.array(zzz))
        int_react = react[0].modify(internals=zzz)
        pos = coordops.zmatrix_indices(zzz, [(2, 0)], strip_embedding=True)
        # react[0].plot(highlight_atoms=[3, 5, 4, 2]).show()
        int_react.animate_coordinate(pos,
                                     strip_embedding=True,
                                     draw_coords=[(2, 1)]).show()


    @validationTest
    def test_ComplexRelaxedScan(self):
        import McUtils.Coordinerds as coordops

        import warnings
        warnings.filterwarnings("ignore", category=RuntimeWarning)  # new mac annoyance

        mol:Molecule = Molecule.from_string('[CH2:1]1[CH2:3]C=C[CH2:4][CH2:2]1',
                                   energy_evaluator='rdkit',
                                   internals='zmatrix'
                                   )

        zm = mol.get_bond_zmatrix(required_coordinates=[
            (0, 2),
            (1, 3),
            (2, 4, 5, 3)
        ])

        int_mol = mol.modify(internals=zm)
        _, geoms, _ = int_mol.relaxed_scan(
            [0, 2, 15],
            {
                (0,2):1,
                (1,3):1
            },
            max_iterations=0,
            coordinate_constraints=[(2, 4, 5, 3)]
            # region_constraints={
            #     (0,1,23)
            # }
        )

        wtf = int_mol.get_internals(geoms, strip_embedding=True)
        idx = coordops.zmatrix_indices(zm,
                                       [
                                           (0, 2),
                                           (1, 3),
                                           (2, 4, 5, 3)
                                       ], strip_embedding=True)
        print(idx)
        # print(
        #     wtf
        # )
        print(
            wtf[..., idx]
        )


        # int_mol.plot(highlight_atoms=[2, 4, 5, 3]).show()
        int_mol.plot(geoms,
                     highlight_atoms=[2, 4, 5, 3]
                     , bonds="recompute"
                     , include_save_buttons=True
                     ).show()

    @validationTest
    def test_ASEProfileGenerator(self):
        import warnings
        warnings.filterwarnings("ignore", category=RuntimeWarning)  # new mac annoyance

        import McUtils.Devutils as dev

        # prod:Molecule = Molecule.from_string("[CH2:1]1[CH2:3]C=C[CH2:4][CH2:2]1")
        # react = prod.apply_smarts("[C:1]1[C:2][C:3]=[C:4][C:5][C:6]1 >> [C:1]=[C:6].[C:2]=[C:3]-[C:4]=[C:5]")
        traj:list[Molecule] = Molecule.from_file(
            TestManager.test_data('traj.xyz'),
            max_blocks=-1,
            units='Angstroms'
        )

        init_js = dev.read_json(TestManager.test_data('product.json'))
        new_js = dev.read_json(TestManager.test_data('trajectory.json'))
        traj = [
            Molecule(init_js['atoms'], c)
            for c in new_js['final_trajectory']
        ]

        from Psience.Reactions import Reaction
        rxn = Reaction([traj[0]], [traj[-1]])
        prof = rxn.get_profile_generator('pys-neb',
                                         energy_evaluator='aimnet2')
        new_images = prof.generate(base_images=traj)
        traj[0].plot([i.coords for i in new_images]).show()

        f1 = plt.Plot(
            prof.evaluate_profile_distances(new_images, normalize=False),
            prof.evaluate_profile_energies(new_images),
        )
        plt.Plot(
            new_js['final_rmsds'],
            new_js['final_energies'],
            figure=f1
        )
        f1.show()

        # traj[0].plot(
        #     [t.coords for t in traj]
        # ).show()
        # prod: Molecule = Molecule.from_string("[CH2:1]1[CH2:3]C=C[CH2:4][CH2:2]1")
        # new = prod.break_bonds([(0, 2), (1, 3)], use_rdkit=True)
        # new.plot().show()
        # return
        # react = prod.apply_smarts("[C:1]1[C:3][C:5]=[C:6][C:4][C:2]1 >> [C:1]=[C:2].[C:3]=[C:5]-[C:6]=[C:4]")
        #
        # # raise Exception(react[0].coords)
        # # react[0].plot().show()

    @validationTest
    def test_PysisGSMProfileGenerator(self):
        import warnings
        warnings.filterwarnings("ignore", category=RuntimeWarning)  # new mac annoyance

        import McUtils.Devutils as dev

        # prod:Molecule = Molecule.from_string("[CH2:1]1[CH2:3]C=C[CH2:4][CH2:2]1")
        # react = prod.apply_smarts("[C:1]1[C:2][C:3]=[C:4][C:5][C:6]1 >> [C:1]=[C:6].[C:2]=[C:3]-[C:4]=[C:5]")
        traj: list[Molecule] = Molecule.from_file(
            TestManager.test_data('traj.xyz'),
            max_blocks=-1,
            units='Angstroms'
        )

        init_js = dev.read_json(TestManager.test_data('product.json'))
        new_js = dev.read_json(TestManager.test_data('trajectory.json'))
        traj = [
            Molecule(init_js['atoms'], c)
            for c in new_js['final_trajectory']
        ]

        from Psience.Reactions import Reaction
        rxn = Reaction([traj[0]], [traj[-1]])
        prof = rxn.get_profile_generator('ase-neb',
                                         energy_evaluator='aimnet2',
                                         spring_constant=1,
                                         climb=True)
        new_images = prof.generate(base_images=traj)
        traj[0].plot([i.coords for i in new_images]).show()

        new_structs = np.array([
            i.coords for i in new_images
        ])

        old_structs = np.array(new_js['final_trajectory'])


        f1 = plt.Plot(
            np.linalg.norm(new_structs[:, 2] - new_structs[:, 3], axis=-1),
            prof.evaluate_profile_energies(new_images),
        )
        plt.Plot(
            np.linalg.norm(old_structs[:, 2] - old_structs[:, 3], axis=-1),
            new_js['final_energies'],
            figure=f1
        )
        f1.show()

        # traj[0].plot(
        #     [t.coords for t in traj]
        # ).show()
        # prod: Molecule = Molecule.from_string("[CH2:1]1[CH2:3]C=C[CH2:4][CH2:2]1")
        # new = prod.break_bonds([(0, 2), (1, 3)], use_rdkit=True)
        # new.plot().show()
        # return
        # react = prod.apply_smarts("[C:1]1[C:3][C:5]=[C:6][C:4][C:2]1 >> [C:1]=[C:2].[C:3]=[C:5]-[C:6]=[C:4]")
        #
        # # raise Exception(react[0].coords)
        # # react[0].plot().show()

    @validationTest
    def test_RMSDChecks(self):
        import warnings
        warnings.filterwarnings("ignore", category=RuntimeWarning)  # new mac annoyance

        import McUtils.Numputils as nput

        # prod:Molecule = Molecule.from_string("[CH2:1]1[CH2:3]C=C[CH2:4][CH2:2]1")
        # react = prod.apply_smarts("[C:1]1[C:2][C:3]=[C:4][C:5][C:6]1 >> [C:1]=[C:6].[C:2]=[C:3]-[C:4]=[C:5]")
        traj: list[Molecule] = Molecule.from_file(
            TestManager.test_data('traj.xyz'),
            max_blocks=-1,
            units='Angstroms'
        )

        print(traj[0].get_rmsd(traj[1]))
        print(
            nput.eckart_rmsd(traj[0].coords, traj[1].coords,
                             masses=traj[0].masses,
                               mass_weighted=True)
        )

    @validationTest
    def test_ASEDimerProfile(self):
        import warnings
        warnings.filterwarnings("ignore", category=RuntimeWarning)  # new mac annoyance

        import McUtils.Devutils as dev

        # prod:Molecule = Molecule.from_string("[CH2:1]1[CH2:3]C=C[CH2:4][CH2:2]1")
        # react = prod.apply_smarts("[C:1]1[C:2][C:3]=[C:4][C:5][C:6]1 >> [C:1]=[C:6].[C:2]=[C:3]-[C:4]=[C:5]")
        init_js = dev.read_json(TestManager.test_data('product.json'))
        new_js = dev.read_json(TestManager.test_data('trajectory.json'))
        traj = [
            Molecule(init_js['atoms'], c)
            for c in new_js['final_trajectory']
        ]

        from Psience.Reactions import Reaction
        rxn = Reaction([traj[0]], [traj[-1]])
        prof = rxn.get_profile_generator('ase-dimer',
                                         energy_evaluator='aimnet2',
                                         climb=True)
        new_images = prof.generate(base_images=traj)
        traj[0].plot([i.coords for i in new_images]).show()

        new_structs = np.array([
            i.coords for i in new_images
        ])

        old_structs = np.array(new_js['final_trajectory'])

        f1 = plt.Plot(
            np.linalg.norm(new_structs[:, 2] - new_structs[:, 3], axis=-1),
            prof.evaluate_profile_energies(new_images),
        )
        plt.Plot(
            np.linalg.norm(old_structs[:, 2] - old_structs[:, 3], axis=-1),
            new_js['final_energies'],
            figure=f1
        )
        f1.show()

        # traj[0].plot(
        #     [t.coords for t in traj]
        # ).show()
        # prod: Molecule = Molecule.from_string("[CH2:1]1[CH2:3]C=C[CH2:4][CH2:2]1")
        # new = prod.break_bonds([(0, 2), (1, 3)], use_rdkit=True)
        # new.plot().show()
        # return
        # react = prod.apply_smarts("[C:1]1[C:3][C:5]=[C:6][C:4][C:2]1 >> [C:1]=[C:2].[C:3]=[C:5]-[C:6]=[C:4]")
        #
        # # raise Exception(react[0].coords)
        # # react[0].plot().show()

    @validationTest
    def test_PysisyphusZTSProfile(self):
        import warnings
        warnings.filterwarnings("ignore", category=RuntimeWarning)  # new mac annoyance

        import McUtils.Devutils as dev

        # prod:Molecule = Molecule.from_string("[CH2:1]1[CH2:3]C=C[CH2:4][CH2:2]1")
        # react = prod.apply_smarts("[C:1]1[C:2][C:3]=[C:4][C:5][C:6]1 >> [C:1]=[C:6].[C:2]=[C:3]-[C:4]=[C:5]")
        init_js = dev.read_json(TestManager.test_data('product.json'))
        new_js = dev.read_json(TestManager.test_data('trajectory.json'))
        traj = [
            Molecule(init_js['atoms'], c)
            for c in new_js['final_trajectory']
        ]

        from Psience.Reactions import Reaction
        rxn = Reaction([traj[0]], [traj[-1]])
        prof = rxn.get_profile_generator('pys-string',
                                         energy_evaluator='aimnet2',
                                         climb=True)
        new_images = prof.generate(base_images=traj)
        traj[0].plot([i.coords for i in new_images]).show()

        new_structs = np.array([
            i.coords for i in new_images
        ])

        old_structs = np.array(new_js['final_trajectory'])

        f1 = plt.Plot(
            np.linalg.norm(new_structs[:, 2] - new_structs[:, 3], axis=-1),
            prof.evaluate_profile_energies(new_images),
        )
        plt.Plot(
            np.linalg.norm(old_structs[:, 2] - old_structs[:, 3], axis=-1),
            new_js['final_energies'],
            figure=f1
        )
        f1.show()

    @validationTest
    def test_PysisyphusDimerProfile(self):
        import warnings
        warnings.filterwarnings("ignore", category=RuntimeWarning)  # new mac annoyance

        import McUtils.Devutils as dev

        # prod:Molecule = Molecule.from_string("[CH2:1]1[CH2:3]C=C[CH2:4][CH2:2]1")
        # react = prod.apply_smarts("[C:1]1[C:2][C:3]=[C:4][C:5][C:6]1 >> [C:1]=[C:6].[C:2]=[C:3]-[C:4]=[C:5]")
        init_js = dev.read_json(TestManager.test_data('product.json'))
        new_js = dev.read_json(TestManager.test_data('trajectory.json'))
        traj = [
            Molecule(init_js['atoms'], c)
            for c in new_js['final_trajectory']
        ]

        from Psience.Reactions import Reaction
        rxn = Reaction([traj[0]], [traj[-1]])
        prof = rxn.get_profile_generator('pys-dimer',
                                         energy_evaluator='aimnet2',
                                         climb=True)
        new_images = prof.generate(base_images=traj)
        traj[0].plot([i.coords for i in new_images]).show()

        new_structs = np.array([
            i.coords for i in new_images
        ])

        old_structs = np.array(new_js['final_trajectory'])

        f1 = plt.Plot(
            np.linalg.norm(new_structs[:, 2] - new_structs[:, 3], axis=-1),
            prof.evaluate_profile_energies(new_images),
        )
        plt.Plot(
            np.linalg.norm(old_structs[:, 2] - old_structs[:, 3], axis=-1),
            new_js['final_energies'],
            figure=f1
        )
        f1.show()

    @validationTest
    def test_ProblemCanonicalZMatrix(self):

        mol = Molecule.from_file(TestManager.test_data('frame_0000.xyz'), units='Angstroms')
        mol.get_canonical_zmatrix()

    @validationTest
    def test_QM9Loading(self):
        from McUtils.ExternalPrograms import QM9

        qm9_path = os.path.expanduser("~/Documents/Postdoc/datasets/qm9.npz")
        supplier = QM9(qm9_path)

        index = np.random.randint(130_000)
        # print(index)

        data = supplier.load_data(index, ['coords', 'atoms'])
        # huh = Molecule.from_string(data['smiles'], coords=data['coords']) # an alternate loader, for strange SMILES strings this can break
        huh = Molecule(data['atoms'], coords=data['coords'] * UnitsData.convert("Angstroms", "BohrRadius"))

        huh.plot(include_save_buttons=True).show()

    @validationTest
    def test_QM9Queries(self):
        from McUtils.ExternalPrograms import QM9, RDMolecule
        import rdkit.Chem as Chem

        qm9_path = os.path.expanduser("~/Documents/Postdoc/datasets/qm9.npz")
        supplier = QM9(qm9_path)

        hmm = supplier.smiles_query('CCC',
                                    upto=500,
                                    # track_failures=True,
                                    # quiet=False,
                                    sanitize=True
                                    )

        bond_lengths = {}
        for idx in hmm:
            data = supplier.load_data(idx, ['coords', 'atoms', 'smiles'])
            rdmol = RDMolecule.parse_smiles(data['smiles'])
            huh = Molecule(
                data['atoms'],
                coords=data['coords'] * UnitsData.convert("Angstroms", "BohrRadius"),
                charge=Chem.GetFormalCharge(rdmol)
            )
            atoms = huh.atoms
            try:
                bonds = huh.bonds
            except ValueError:
                continue

            i,j = np.array([b[:2] for b in bonds]).T
            bls = np.linalg.norm(huh.coords[i, :] - huh.coords[j, :], axis=-1) * UnitsData.bohr_to_angstroms
            for (b1,b2,t),r in zip(bonds, bls):
                a1, a2 = sorted([atoms[b1], atoms[b2]])
                bond_lengths.setdefault((a1, a2, t), []).append(r)

        import pprint
        pprint.pprint({
            k:np.average(v)
            for k,v in bond_lengths.items()
        })

    @validationTest
    def test_QM9Iter(self):
        from McUtils.ExternalPrograms import QM9, RDMolecule
        import rdkit.Chem as Chem

        qm9_path = os.path.expanduser("~/Documents/Postdoc/datasets/qm9.npz")
        supplier = QM9(qm9_path)

        bond_lengths = {}
        update_interval = 500

        for n,data in enumerate(supplier.data_iter(['coords', 'atoms', 'smiles'], start_at=65500)):
            rdmol = RDMolecule.parse_smiles(data['smiles'], quiet=True, sanitize=True)
            if rdmol is None: continue
            if n % update_interval == 0:
                import pprint
                pprint.pprint({
                    k: np.average(v).tolist()
                    for k, v in bond_lengths.items()
                })

            huh = Molecule(
                data['atoms'],
                coords=data['coords'] * UnitsData.convert("Angstroms", "BohrRadius"),
                charge=Chem.GetFormalCharge(rdmol)
            )
            atoms = huh.atoms
            try:
                bonds = huh.bonds
            except ValueError:
                continue

            i, j = np.array([b[:2] for b in bonds]).T
            bls = np.linalg.norm(huh.coords[i, :] - huh.coords[j, :], axis=-1) * UnitsData.bohr_to_angstroms
            for (b1, b2, t), r in zip(bonds, bls):
                a1, a2 = sorted([atoms[b1], atoms[b2]])
                bond_lengths.setdefault((a1, a2, t), []).append(r)

        import pprint
        pprint.pprint({
            k: np.average(v)
            for k, v in bond_lengths.items()
        })

    @validationTest
    def test_RDKitPlotErrors(self):
        Molecule.from_string(
            r'''CO/N=C(\\C(O)=N[C@@H]1C(=O)N2C(C(=O)O)=C(C[N+:7]3=[C:6]=[C:5]4CC=C/[C:3]([S:2](=O)(=O)[OH:1])=[C:4]\\4C=C3)CS[C@H]12)c1csc(=N)[nH]1'''
        ).plot(backend='rdkit').show()

    @validationTest
    def test_RDKitPlotPanel(self):
        systems = [
            'CON=[C:3]([C:2](=O)[NH:1]C)/[C:4]1=[CH:5]/[C:6]2=[N+:7](C(C)C(C)=NN(c3ccccc3)c3ccccc3)OC1=C(C)C2',
            'CO/N=C(\\C(O)=N[C@@H]1C(=O)N2C(C(=O)O)=C(C[N+:7]3=[C:6]=[C:5]4CC=C/[C:3]([S:2](=O)(=O)[OH:1])=[C:4]\\4C=C3)CS[C@H]12)c1csc(=N)[nH]1',
            'C=C/C(=C\\N)C/C=C(\\C=N)c1ccc2c(ccc3c4ccc(C5=C[CH:2]([C:3]6=[N+:7]=[CH:6]/[CH:5]=[CH:4]/6)[NH:1]C(/C(C=C)=C/N)=C5)cc4n(-c4ccccc4)c23)c1',
            'C=C/C(=C\\N)C/C=C(\\C=N)c1ccc2c(ccc3c4ccc(C5=C[CH:2]([C:3]6=[N+:7]=[CH:6]/[CH:5]=[CH:4]/6)[NH:1]C(/C(C=C)=C/N)=C5)cc4n(-c4ccccc4)c23)c1',
            'CCCCCOc1cc(C2C(=O)C(=[C:3]3[C:2]([NH:1]S(C)(=O)=O)=C[C:6](=[N+:7](C)Cc4ccccc4)/[C:5](C)=[CH:4]/3)C2[O-])c(NS(C)(=O)=O)cc1N(C)Cc1ccccc1',
            'C=C/C(CNC(=O)c1ccnc(NC[C:6]2=[N+:7](C)[C:3]([C:2]3=CCN=C[NH:1]3)/[N:4]=[N:5]/2)c1)=C(\\N=C)C(F)(F)F',
            'C=C/C(CNC(=O)c1ccnc(NC[C:6]2=[N+:7](C)[C:3]([C:2]3=CCN=C[NH:1]3)/[N:4]=[N:5]/2)c1)=C(\\N=C)C(F)(F)F',
            'CCCCCCCCCCCCC/C=C/C(O)C(CO)NC(=O)Cn1nnc2c1CC[C@@H]1[C@H](CC2)[C@@H]1COC(=O)CC[NH:1][C:2]1=[C:3]2/[CH:4]=[CH:5]/[CH:6]=[N+:7]2[B-](F)(F)n2cccc21',
            'COc1c(Nc2ccc(C)cc2)nc(N2CCOCC2)nc1C(=N)/C=[C:3]([CH:2]=[NH:1])/[C:4]1=[CH:5]/[CH:6]=[N+:7]1C',
            'CCCCCCCC[N:3]1[CH2:2][NH:1]C[C:5](/[CH:6]=[N+:7](/C)[O-])=[CH:4]\\1',
            'CCCC(CNC)OC1C=CC(C(C)[NH:1][CH2:2][C:3]2=[N+:7]=[C:6]3C=C(C(=C(C)NC4CCCC4)C(C)N)C(=O)CCC(CC)/[C:5]3=[N:4]\\2)CC1',
            'CC1(CCCCCC(=O)O)c2cc(S(=O)(=O)O)ccc2[N+:7](CCCCS(=O)(=O)O)=[C:6]1/[CH:5]=[CH:4]/[CH:3]=[CH:2][NH:1]c1ccccc1',
            'C=[N+:7]=[C:6]1C[N:3]([C:2](=Nc2cccc(OC(F)F)c2)[NH:1]C#N)/[N:4]=[C:5]/1c1ccc(Cl)c(C)c1',
            'C#CC(=NC1C(C)=C(OCCN2CCCOCC2)C=CC1Cl)[NH:1][C:2]1=[N+:7]=[C:6]=[CH:5]/[C:4](C(F)(F)F)=[CH:3]/1',
            'O/C=C/C=C/CC/C=C/CCCC[N+:7]1=[CH:6]/[CH:5]=[CH:4]/[C@:3]2(CC/C=C/CCCCCC[NH:1][CH2:2]2)C1'
        ]
        mols = [Molecule.from_string(b, 'smi').get_embedded_molecule() for b in systems]
        vd = 100
        offset = 25
        lineskip = 30
        nrows = len(mols) // 2
        baseline = nrows // 2 * (-lineskip)
        fs = 8
        fig = mols[0].plot(
            backend='rdkit',
            plot_range=10.1 * np.array([[-1, 1], [-1, 1], [-1, 1]]),
            image_size=[1500, 1500],
            # background='pink',
            view_settings={
                'view_distance': vd,
                'view_center': [-offset, baseline, 0]
            },
            font_size=fs
        )
        for i, m in enumerate(mols[1:]):
            i = i + 1
            i, j = i // 2, i % 2
            m.plot(figure=fig,
                   backend='rdkit',
                   view_settings={
                       'view_distance': vd,
                       'view_center': [-offset if j == 0 else offset, baseline + lineskip * i, 0]
                   },
                   font_size=fs
                   )
        fig.show()
    @inactiveTest
    def test_Manipulator(self):
        from Psience.Molecools import Molecule
        import McUtils.Numputils as nput
        import numpy as np
        bits = Molecule.from_string("C1CCOCC1=O", "smi", confgen_opts=dict(random_seed=21232)).get_embedded_molecule()

        vv = [0, 0, 1]
        uv = [1, 0, 0]

        def prep_up_vector(view_angle):
            if view_angle is None or isinstance(view_angle, str):
                view_angle = 0
            return nput.rotation_matrix(vv, view_angle) @ np.array(uv)

        def plot_mol(view_angle, view_dist, view_x=0, view_z=0, **etc):
            uuuh = bits.plot(backend='matplotlib',
                             image_size=[500, 500],
                             highlight_atoms=[0, 1, 2],
                             draw_coords={
                                 (0, 4): {
                                     'label': "r",
                                     'line_color': 'pink',
                                     'label_style': {'font_size': 20, 'color': 'red'}
                                 },
                                 (0, 1, 2): {
                                     'label': "a",
                                     'line_color': 'pink',
                                     'label_style': {'font_size': 20, 'color': 'blue'}
                                 }
                             },
                             plot_range=[[-3, 3], [-3, 3], [-3, 3]],
                             view_settings={
                                 'view_distance': (
                                     5
                                     if view_dist is None or isinstance(view_dist, str) else
                                     view_dist
                                 ),
                                 'up_vector': prep_up_vector(view_angle),
                                 'view_vector': vv
                             },
                             principle_axes=True,
                             include_save_buttons=True)
            return uuuh.to_widget()

        import McUtils.Jupyter as interactive
        woof = interactive.Manipulator(
            plot_mol,
            ['view_angle', {'range': [-np.pi, np.pi], 'value': 0, 'continuous_update': True}],
            ['view_dist', {'range': [1, 20], 'value': 10}],
        )
        woof

    @inactiveTest
    def test_ManipulatorRDKit(self):
        from Psience.Molecools import Molecule
        import McUtils.Numputils as nput
        import numpy as np
        bits = Molecule.from_string("C1CCOCC1=O", "smi").get_embedded_molecule()

        def prep_up_vector(view_angle):
            if view_angle is None or isinstance(view_angle, str):
                view_angle = 0
            return nput.rotation_matrix([0, 0, 1], view_angle) @ np.array([0, 1, 0])

        def plot_mol(view_angle, view_dist, view_x=0, view_z=0, **etc):
            uuuh = bits.plot(backend='rdkit',
                             image_size=[500, 500],
                             highlight_atoms=[0, 1, 2],
                             draw_coords={
                                 (0, 4): "r",
                                 (0, 1, 2): {
                                     'label': "t",
                                     'label_style': {'font_size': 22, 'font_family': 'Arial', 'color': 'red'}
                                 }
                             },
                             postdraw=[
                                 {
                                     'pattern': 'labels_1',
                                     'classes': True,
                                     'replacement': {'text': "θ", 'mode': 'svg'}
                                 }
                             ],
                             view_settings={
                                 'view_distance': (
                                     5
                                     if view_dist is None or isinstance(view_dist, str) else
                                     view_dist
                                 ),
                                 'up_vector': prep_up_vector(view_angle),
                                 'view_vector': [0, 0, 1]
                             })
            return uuuh.to_widget()