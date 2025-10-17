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

    @debugTest
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

