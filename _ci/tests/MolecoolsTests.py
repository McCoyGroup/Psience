
from Peeves.TestUtils import *
from unittest import TestCase
from Psience.Molecools import Molecule, MolecularNormalModes
from Psience.Molecools.Transformations import MolecularTransformation
from Psience.Data import DipoleSurface # this will be leaving Zachary very soon I think...
from McUtils.GaussianInterface import GaussianFChkReader, GaussianLogReader
from McUtils.Plots import *
from McUtils.Coordinerds import cartesian_to_zmatrix
from McUtils.Data import UnitsData
import numpy as np

class MolecoolsTests(TestCase):

    test_log_water = TestManager.test_data("water_OH_scan.log")
    test_log_freq = TestManager.test_data("water_freq.log")
    test_HOD = TestManager.test_data("HOD_freq.fchk")
    test_fchk = TestManager.test_data("water_freq.fchk")
    test_log_h2 = TestManager.test_data("outer_H2_scan_new.log")

    @validationTest
    def test_ImportMolecule(self):

        n = 3 # water
        m = Molecule.from_file(self.test_fchk)
        self.assertEquals(m.atoms, ("O", "H", "H"))

    @debugTest
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

        #transf = scan.principle_axis_frame(sel=sel)
        transf = scan.eckart_frame(ref, sel=sel)
        tf_test = transf[0].transformation_function

        tf_mat = tf_test.transform
        self.assertTrue(np.allclose(tf_mat@tf_mat.T - np.eye(3), 0.))
        self.assertEquals(tf_test.transf.shape, (4, 4))

        for t, m in zip(transf, scan):
            new_mol = t(m)
            # rot_ref.guess_bonds = False
            # ref.guess_bonds = False
            # m.guess_bonds = False
            # new_mol.guess_bonds = False
            # m = m #type: Molecule
            # g1, a, b = ref.plot()
            # # ref.plot(figure=g1)
            # rot_ref.plot(figure=g1)
            # g, a, b = new_mol.plot()
            # rot_ref.plot(figure=g, atom_style=dict(color='black'))
            # g.show()
            fuckup = np.linalg.norm(new_mol.coords[sel] - rot_ref.coords[sel])
            self.assertLess(fuckup/len(sel), .1)

    @validationTest
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
            tuple(np.round(UnitsData.convert("Hartrees", "Wavenumbers")*nms.freqs, 4)),
            tuple(np.round(test_freqs, 4))
        )

    @debugTest
    def test_InternalCartesianJacobians(self):
        import McUtils.Plots as plt
        m = Molecule.from_file(TestManager.test_data('HOH_freq.fchk'),
                               zmatrix=[
                                   [0, -1, -1, -1],
                                   [1,  0, -1, -1],
                                   [2,  0,  1, -1]
                               ]
                               )
        m = m.get_embedded_molecule()
        intcds = m.internal_coordinates
        carts = m.coords
        ijacs = intcds.jacobian(carts.system, 1, mesh_spacing=1.0e-5)
        ijacs2num = intcds.jacobian(carts.system, 2, mesh_spacing=1.0e-2, all_numerical=True)
        jacs, jacs2 = carts.jacobian(intcds.system, [1, 2], mesh_spacing=1.0e-5)

        meh1 = ijacs.squeeze().reshape(9, 9)
        meh2 = jacs.squeeze().reshape(9, 9)

        itest = np.dot(meh1, meh2)
        itest2 = np.dot(meh2, meh1)

        # plt.ArrayPlot(meh1)
        # plt.ArrayPlot(meh2)
        # plt.ArrayPlot(itest2)
        # plt.ArrayPlot(itest).show()
        self.assertTrue(np.allclose(np.eye(9), itest))

        good_sel = (...,) + np.ix_((3, 5, 6), (3, 5, 6))
        # meh12 = ijacs2.squeeze().reshape(9, 9, 9)
        # meh12 = meh12.transpose(2, 0, 1).reshape(3, 3, 9, 9)
        # meh12 = meh12[good_sel]
        meh22 = ijacs2num.squeeze().reshape(9, 9, 9)
        meh22 = meh22.transpose(2, 0, 1).reshape(3, 3, 9, 9)
        meh22 = meh22[good_sel]
        ps = dict(vmin=-.05, vmax=.05)
        # plt.TensorPlot(meh12, plot_style=ps)
        plt.TensorPlot(meh22, plot_style=ps).show()
        # plt.TensorPlot(meh22-meh12, plot_style=ps).show()
