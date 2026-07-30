"""Extracted from DGBTests.test_ModelPotentialAIMD2D via McUtils.Docs.ExamplesParser -- not the original file, and may reference test-only setup/state. Run with: python -m unittest DGBTests.test_ModelPotentialAIMD2D"""

import collections
import gc
import itertools
import os
import shutil
import scipy.linalg
import McUtils.Zachary
from Peeves.TestUtils import *
from Peeves import BlockProfiler
from unittest import TestCase
from McUtils.Data import UnitsData, PotentialData, AtomData
from McUtils.Zachary import Interpolator, FiniteDifferenceDerivative
import McUtils.Plots as plt
import McUtils.Numputils as nput
import McUtils.Devutils as dev
from McUtils.GaussianInterface import GaussianFChkReader, GaussianLogReader
from McUtils.Extensions import ModuleLoader
from McUtils.Scaffolding import Checkpointer
from Psience.DGB import *
from Psience.Molecools import Molecule
from Psience.AIMD import AIMDSimulator
import numpy as np

class DGBTests(TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        np.set_printoptions(linewidth=int(100000000.0))
    w2h = UnitsData.convert('Wavenumbers', 'Hartrees')

    @staticmethod
    def multiply_model_functions(func1: 'dict[int|tuple, dict|tuple]', func2: 'dict[int|tuple, dict|tuple]'):
        new_func = {}
        for k, f1 in func1.items():
            for k2, f2 in func2.items():
                if isinstance(k, int):
                    k = (k,)
                    f1 = (f1,)
                if isinstance(k2, int):
                    k2 = (k2,)
                    f2 = (f2,)
                new_func[k + k2] = f1 + f2
        return new_func

    @classmethod
    def buildWaterModel(cls, *, oh_model=False, atoms=None, w=3869.47 * w2h, wx=84 * w2h, w2=3869.47 * w2h, wx2=84 * w2h, ka=1600 ** 2 / 150 * w2h, dudr1=None, dudr2=None, duda=None, dipole=None, dipole_magnitude=None, dipole_direction='auto'):
        base_water = Molecule.from_file(TestManager.test_data('water_freq.fchk'), internals=[[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]])
        if oh_model:
            w2 = None
            wx2 = None
            ka = None
            mol = Molecule(base_water.atoms[:2], base_water.coords[:2], internals=[[0, -1, -1, -1], [1, 0, -1, -1]]).get_embedded_molecule(load_properties=False)
        elif atoms is not None:
            mol = Molecule(atoms, base_water.coords, internals=[[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]]).get_embedded_molecule(load_properties=False)
        else:
            mol = base_water.get_embedded_molecule(load_properties=False)
        r1 = 0
        r2 = 1
        a12 = 2
        potential_params = {}
        if wx is not None:
            potential_params[r1] = {'morse': {'w': w, 'wx': wx}}
        elif w is not None:
            potential_params[r1] = {'harmonic': {'k': w}}
        if wx2 is not None:
            potential_params[r2] = {'morse': {'w': w2, 'wx': wx2}}
        elif w2 is not None:
            potential_params[r2] = {'harmonic': {'k': w2}}
        elif w2 is not None:
            potential_params[r2] = {'harmonic': {'k': w2}}
        if ka is not None:
            potential_params[a12] = {'harmonic': {'k': ka}}
        if dipole is None and dudr1 is not None or dudr2 is not None or duda is not None:
            if oh_model:
                dipole = [{r1: {'linear': {'eq': 0, 'scaling': dudr1}}}, 0, 0]
            else:
                dipole = [{(r1, a12): ({'linear': {'eq': 0, 'scaling': dudr1}}, {'sin': {'eq': 0}}), (r2, a12): ({'linear': {'eq': 0, 'scaling': dudr2}}, {'sin': {'eq': 0, 'scaling': -1}})}, {(r1, a12): ({'linear': {'eq': 0, 'scaling': dudr1}}, {'cos': {'eq': 0, 'scaling': 1 / 2}}), (r2, a12): ({'linear': {'eq': 0, 'scaling': dudr2}}, {'cos': {'eq': 0, 'scaling': 1 / 2}})}, 0]
        return (mol, mol.get_model(potential_params, dipole=dipole))

    @classmethod
    def buildTrajectory(cls, mol, cart_pot_func, steps, timestep=0.5, initial_energies=None, initial_displacements=None, displaced_coords=None, seed=0):
        if initial_displacements is not None:
            init_pos = mol.get_displaced_coordinates(initial_displacements, which=displaced_coords, use_internals='reembed')
            sim = AIMDSimulator(mol.masses, init_pos, lambda c: -cart_pot_func(c, deriv_order=1)[1].reshape(c.shape), timestep=timestep, track_kinetic_energy=True)
        else:
            mol.potential_derivatives = cart_pot_func(mol.coords, deriv_order=2)[1:]
            nms = mol.normal_modes.modes.basis
            sim = AIMDSimulator(mol.atomic_masses, [mol.coords] * len(initial_energies), lambda c: -cart_pot_func(c, deriv_order=1)[1].reshape(c.shape), velocities=AIMDSimulator.mode_energies_to_velocities(nms.inverse.T, mol.atomic_masses, initial_energies, inverse=nms.matrix.T), timestep=timestep, track_kinetic_energy=True)
        sim.propagate(steps)
        raise Exception(np.array(sim.trajectory).shape)
        coords = np.array(sim.trajectory).reshape((-1,) + mol.coords.shape)
        coords = mol.embed_coords(coords)
        return (coords, sim)

    @staticmethod
    def buildRunDGB(coords, pot, dipole, *, logger=True, plot_wavefunctions=True, plot_spectrum=True, **opts):
        dgb = DGB.construct(np.round(coords, 8), pot, logger=logger, **opts)
        logger = dgb.logger
        with logger.block(tag='Running DGB'):
            logger.log_print('num coords: {c}', c=len(dgb.gaussians.coords.centers))
            with logger.block(tag='S'):
                logger.log_print(logger.prep_array(dgb.S[:5, :5]))
            with logger.block(tag='T'):
                logger.log_print(logger.prep_array(dgb.T[:5, :5]))
            with logger.block(tag='V'):
                logger.log_print(logger.prep_array(dgb.V[:5, :5]))
            wfns_cart = dgb.get_wavefunctions()
            with logger.block(tag='Energies'):
                logger.log_print(logger.prep_array(wfns_cart.energies[:5] * UnitsData.convert('Hartrees', 'Wavenumbers')))
            with logger.block(tag='Frequencies'):
                logger.log_print(logger.prep_array(wfns_cart.frequencies()[:5] * UnitsData.convert('Hartrees', 'Wavenumbers')))
            if plot_wavefunctions:
                for i in range(4):
                    wfns_cart[i].plot().show()
            spec = wfns_cart[:4].get_spectrum(dipole)
            with logger.block(tag='Intensities'):
                logger.log_print(logger.prep_array(spec.intensities))
            if plot_spectrum:
                spec.plot().show()

    @classmethod
    def setupCartesianModelDGB(cls):
        ...

    @classmethod
    def plot_dgb_potential(cls, dgb, mol, potential, coordinate_sel=None, domain=None, domain_padding=1, potential_cutoff=17000, potential_min=0, plot_cartesians=None, plot_atoms=True, cmap=None, modes_nearest=False, plot_points=100, levels=24, **plot_styles):

        def cutoff_pot(points, cutoff=potential_cutoff / UnitsData.hartrees_to_wavenumbers, cutmin=potential_min / UnitsData.hartrees_to_wavenumbers):
            values = potential(points)
            values[np.logical_or(values < cutmin, values > cutoff)] = cutoff
            return values
        if isinstance(dgb, DGBCoords):
            coords = dgb
        else:
            coords = dgb.gaussians.coords
        if plot_cartesians is None:
            plot_cartesians = isinstance(coords, DGBCartesians)
            if coordinate_sel is None:
                coordinate_sel = [0, 1]
        if plot_cartesians:
            figure = mol.plot_molecule_function(cutoff_pot, axes=coordinate_sel, domain=domain, modes_nearest=modes_nearest, domain_padding=domain_padding, cmap=cmap, plot_points=plot_points, levels=levels, plot_atoms=plot_atoms, **plot_styles)
        else:
            if domain is None:
                from McUtils.Zachary import Mesh
                domain = Mesh(coords.centers).bounding_box
            points = DGBWavefunction.prep_plot_grid(domain, domain_padding=domain_padding)
            if points.shape[-1] == 1:
                figure = plt.Plot(*np.moveaxis(points, -1, 0), potential(points), **plot_styles)
            else:
                figure = plt.TriContourPlot(*np.moveaxis(points, -1, 0), cutoff_pot(points), cmap=cmap, levels=levels, **plot_styles)
        return figure

    @classmethod
    def plot_gaussians(cls, dgb, mol, *, domain=None, domain_padding=1, cmap='RdBu', plot_dir=None, plot_name='gaussian_{i}.pdf', **plot_options):
        n = len(dgb.gaussians.coords.centers)
        wfns = DGBWavefunctions(np.zeros(n), np.eye(n), dgb)
        cls.plot_wavefunctions(wfns, dgb, mol, cmap=cmap, plot_name=plot_name, plot_dir=plot_dir, domain=domain, domain_padding=domain_padding, potential_styles=dict(domain=domain, domain_padding=domain_padding), **plot_options)
    default_num_plot_wfns = 5

    @classmethod
    def plot_wavefunctions(cls, wfns, dgb, mol, which=True, coordinate_sel=None, cartesians=None, plot_dir=None, plot_name='wfn_{i}.pdf', plot_label='{e} cm-1', plot_potential=True, plot_atoms=None, plot_centers=True, potential_styles=None, **plot_options):
        figure = None
        if cartesians:
            wfns = wfns.as_cartesian_wavefunction()
            dgb = wfns.hamiltonian
        if coordinate_sel is None:
            coordinate_sel = list(range(dgb.gaussians.alphas.shape[-1]))
        figs = []
        if which is True:
            which = cls.default_num_plot_wfns
        if isinstance(which, int):
            which = list(range(which))
        if isinstance(which, list) and len(which) > 0:
            pot = dgb.pot.potential_function
            if plot_potential:
                if plot_atoms is None:
                    plot_atoms = bool(plot_centers)
                if potential_styles is None:
                    potential_styles = {}
                pot_figure = cls.plot_dgb_potential(dgb, mol, pot, coordinate_sel=coordinate_sel, plot_atoms=plot_atoms, **potential_styles)
            else:
                pot_figure = None
            if plot_dir is not None:
                os.makedirs(plot_dir, exist_ok=True)
            if coordinate_sel is None and wfns.gaussians.alphas.shape[-1] == 1 or len(coordinate_sel) == 1:
                for k in ['cmap', 'levels', 'plotter', 'contour_levels']:
                    if k in plot_options:
                        del plot_options[k]
                    if k in potential_styles:
                        del potential_styles[k]
            for i in which:
                if i < len(wfns):
                    if pot_figure is not None:
                        figure = pot_figure.copy()
                    w = wfns[i]
                    if i == 0:
                        e = w.energy
                    else:
                        e = w.energy - wfns.energies[0]
                    e = e * UnitsData.hartrees_to_wavenumbers
                    lab = plot_label.format(i=i, e=e) if plot_label is not None else None
                    if isinstance(dgb.gaussians.coords, DGBCartesians):
                        figs.append(w.plot_cartesians(coordinate_sel, plot_centers=plot_centers, figure=figure, plot_label=lab, **plot_options))
                    elif wfns.gaussians.alphas.shape[-1] > 1:
                        if coordinate_sel is not None:
                            w = w.project(coordinate_sel)
                        figs.append(w.plot(figure=figure, plot_centers=plot_centers, plot_label=lab, **plot_options))
                    else:
                        figs.append(w.plot(figure=figure, plot_centers=plot_centers, plot_label=lab, **plot_options))
                    if plot_dir is not None:
                        fig = figs.pop()
                        fig.savefig(os.path.join(plot_dir, plot_name.format(i=i)))
                        fig.close()
            if plot_dir is None:
                figs[0].show()

    @classmethod
    def runDGB(cls, dgb: DGB, mol, plot_centers=True, plot_wavefunctions=True, plot_spectrum=False, pot_cmap='viridis', wfn_cmap='RdBu', wfn_points=100, wfn_contours=12, plot_dir=None, plot_potential=True, pot_points=100, domain=None, domain_padding=1, potential_cutoff=15000, mode=None, nodeless_ground_state=None, min_singular_value=None, subspace_size=None, plot_similarity=False, similarity_cutoff=None, similarity_chunk_size=None, similar_det_cutoff=None, **plot_options):
        print('--->', len(dgb.gaussians.coords.centers))
        print(dgb.S[:5, :5])
        print(dgb.T[:5, :5])
        print(dgb.V[:5, :5])
        try:
            wfns, spec = dgb.run(calculate_spectrum=plot_spectrum, mode=mode, nodeless_ground_state=nodeless_ground_state, subspace_size=subspace_size, min_singular_value=min_singular_value, similarity_cutoff=similarity_cutoff, similarity_chunk_size=similarity_chunk_size, similar_det_cutoff=similar_det_cutoff)
        except Exception as e:
            if plot_wavefunctions is not False:
                print(e)
                if isinstance(plot_wavefunctions, str) and plot_wavefunctions == 'cartesians':
                    plot_wavefunctions = {'cartesians': None}
                cartesian_plot_axes = None
                if isinstance(plot_wavefunctions, dict):
                    if 'cartesians' in plot_wavefunctions:
                        cartesian_plot_axes = plot_wavefunctions['cartesians']
                        dgb = dgb.as_cartesian_dgb()
                    else:
                        raise ValueError(plot_wavefunctions)
                pot = dgb.pot.potential_function
                figure = cls.plot_dgb_potential(dgb, mol, pot, domain=domain, domain_padding=domain_padding, potential_cutoff=potential_cutoff)
                dgb.gaussians.plot_centers(figure, xyz_sel=cartesian_plot_axes)
                figure.show()
            raise e
        else:
            if plot_similarity:
                plt.ArrayPlot(dgb.get_similarity_matrix()).show()
            use_cartesians = False
            if isinstance(plot_wavefunctions, str) and plot_wavefunctions == 'cartesians':
                plot_wavefunctions = {'cartesians': None}
            coordinate_sel = None
            if isinstance(plot_wavefunctions, dict):
                if 'cartesians' in plot_wavefunctions:
                    use_cartesians = True
                    coordinate_sel = plot_wavefunctions['cartesians']
                    plot_wavefunctions = plot_wavefunctions.get('num', True)
                elif 'modes' in plot_wavefunctions:
                    coordinate_sel = plot_wavefunctions['modes']
                    plot_wavefunctions = plot_wavefunctions.get('num', True)
                    plot_potential = False
                else:
                    raise ValueError(plot_wavefunctions)
            if plot_spectrum:
                which = plot_wavefunctions
                if which is True:
                    which = cls.default_num_plot_wfns
                if isinstance(which, int):
                    which = list(range(which))
                spec_plot = spec[which].plot()
                if plot_dir is not None:
                    os.makedirs(plot_dir, exist_ok=True)
                    spec_plot.savefig(os.path.join(plot_dir, 'spec.png'))
                    spec_plot.close()
            if plot_wavefunctions:
                cls.plot_wavefunctions(wfns, dgb, mol, which=plot_wavefunctions, cartesians=use_cartesians, coordinate_sel=coordinate_sel, plot_dir=plot_dir, contour_levels=wfn_contours, cmap=wfn_cmap, plot_points=wfn_points, plot_centers={'color': 'red'} if plot_centers else False, domain=domain, domain_padding=domain_padding, plot_potential=plot_potential, scaling=0.2, plotter=plt.TriContourLinesPlot, potential_styles=dict(domain=domain, domain_padding=domain_padding, cmap=pot_cmap, plot_points=pot_points), **plot_options)
            if plot_spectrum:
                return (wfns, spec)
            else:
                return wfns

    @classmethod
    def getMorseParameters(cls, w=None, wx=None, m1=None, m2=None, re=None):
        if w is None:
            w = 3869.47 * cls.w2h
        freq = w
        if wx is None:
            wx = 84 * cls.w2h
        anh = wx
        De = freq ** 2 / (4 * anh)
        if m1 is None:
            m1 = AtomData['O', 'Mass'] * UnitsData.convert('AtomicMassUnits', 'AtomicUnitOfMass')
        if m2 is None:
            m2 = AtomData['H', 'Mass'] * UnitsData.convert('AtomicMassUnits', 'AtomicUnitOfMass')
        muv = 1 / m1 + 1 / m2
        a = np.sqrt(2 * anh / muv)
        if re is None:
            re = 1.82534
        return (De, a, re)

    @classmethod
    def setupMorseFunction(cls, model, i, j, w=None, wx=None):
        from McUtils.Data import PotentialData
        if isinstance(model, float):
            m1 = model
            m2 = i
            re = j
        else:
            m1 = model.vals[model.m(i)]
            m2 = model.vals[model.m(j)]
            re = model.vals[model.r(i, j)]
        De, a, re = cls.getMorseParameters(w=w, wx=wx, m1=m1, m2=m2, re=re)

        def morse_basic(r, re=re, alpha=a, De=De, deriv_order=None, _morse=PotentialData['MorsePotential']):
            return _morse(r, re=re, alpha=alpha, De=De, deriv_order=deriv_order)
        return morse_basic

    @classmethod
    def plot_interpolation_error(cls, dgb, pot):
        sel = slice(None)
        centers = dgb.gaussians.overlap_data.centers[sel]
        embpot = dgb.gaussians.coords.embed_function(pot)
        realpots, realgrad, real_hess = [d * 219475 for d in embpot(centers, deriv_order=2)]
        interpots, intergrad, inter_hess = [d * 219475 for d in dgb.pot.potential_function(centers, deriv_order=2)]
        ords = np.argsort(realpots)
        devs = interpots[ords] - realpots[ords]
        rows, cols = np.triu_indices_from(dgb.S)
        utris = dgb.S[rows, cols]
        unscaled_devs = devs
        devs = devs * utris
        max_dev_pos = np.flip(np.argsort(np.abs(devs)))[:5]
        inter_gnorm = np.sum(np.abs(intergrad[ords]), axis=-1)
        real_gnorm = np.sum(np.abs(realgrad[ords]), axis=-1)
        dev_gnorm = inter_gnorm - real_gnorm
        grad_plot = plt.ScatterPlot(realpots[ords], dev_gnorm)
        inter_trace = np.sum(np.sum(np.abs(inter_hess[ords]), axis=-1), axis=-1)
        real_trace = np.sum(np.sum(np.abs(real_hess[ords]), axis=-1), axis=-1)
        dev_trace = inter_trace - real_trace
        dev_hess = inter_hess - real_hess
        hess_plot = plt.ScatterPlot(realpots[ords], dev_trace)
        print('Mean Absolute Error:', np.mean(np.abs(unscaled_devs)), 'Std:', np.std(unscaled_devs))
        print('Mean Scaled Error:', np.mean(np.abs(devs)), 'Std:', np.std(devs))
        print('Mean Hessian Error:', np.mean(np.abs(dev_hess.flatten())), 'Std:', np.std(dev_hess.flatten()), 'Max:', np.max(np.abs(dev_hess.flatten())))
        print('Mean Summed Hessian Error:', np.mean(np.abs(dev_trace.flatten())), 'Std:', np.std(dev_trace.flatten()), 'Max:', np.max(np.abs(dev_trace.flatten())))
        print('Maximum (Scaled) Error:', devs[max_dev_pos])
        print('Maximum (Scaled) Interpolation Error:')
        for l, r, c, tt, ii, ov in zip(dgb.gaussians.coords.centers[rows[sel][ords[max_dev_pos]]], dgb.gaussians.coords.centers[cols[sel][ords[max_dev_pos]]], dgb.gaussians.overlap_data.centers[sel][ords[max_dev_pos]], realpots[ords[max_dev_pos]], interpots[ords[max_dev_pos]], utris[sel][ords[max_dev_pos]]):
            print(f'Centers: {c} ({ov}) <- {l} {r}')
            print(f'  Error: {ii - tt} <- {tt} {ii}')
        bad_bad = np.abs(devs) > 50
        dev_plot = plt.ScatterPlot(realpots[ords], unscaled_devs)
        dev_plot.show()

    def symmetrizePoints(self, coords, equivalent_atoms, unique_tol=7):
        coords = np.asanyarray(coords)
        all_perms = [np.unique(list(itertools.permutations(a)), axis=0) for a in equivalent_atoms]
        stacks = [coords]
        for swaps in itertools.product(*all_perms):
            perm = np.arange(coords.shape[1])
            for og, sw in zip(equivalent_atoms, swaps):
                perm[og] = perm[sw]
            stacks.append(coords[:, perm])
        coords = np.concatenate(stacks, axis=0)
        if unique_tol is not None:
            _, upos = np.unique(np.round(coords, unique_tol), return_index=True, axis=0)
            coords = coords[np.sort(upos)]
        return coords

    @classmethod
    def declusterPoints(cls, points, radius):
        pivots = np.arange(len(points))
        dec_pts = points.reshape(points.shape[0], -1)
        for i in range(len(points)):
            cur_pos = pivots[i]
            dists = np.linalg.norm(dec_pts[cur_pos][np.newaxis, :] - dec_pts[pivots[i + 1:], :], axis=1)
            good_pos = np.where(dists > radius)
            if len(good_pos) == 0 or len(good_pos[0]) == 0:
                break
            pivots = np.concatenate([pivots[:i + 1], pivots[i + 1:][good_pos]])
        return points[pivots]

    @classmethod
    def plot_quadratic_opt_potentials(cls, interp_data, mol, local_coordinates=None):
        modes = mol.normal_modes.modes.basis.to_new_modes()
        print(modes.matrix.shape)
        print([v.shape for v in interp_data['values']])
        embedded_coords = modes.embed_coords(interp_data['centers'])
        embedded_vals = modes.embed_derivs(interp_data['values'])
        fitting_vals = embedded_vals[0] * 219475.6 < 3600
        cents = interp_data['centers'][fitting_vals]
        embedded_coords = embedded_coords[fitting_vals]
        embedded_vals = [v[fitting_vals] for v in embedded_vals]

        def bond_vec(coords, i, j):
            bond_tf = nput.dist_deriv(coords, i, j)[1]
            bond_vectors = np.zeros(coords.shape)
            bond_vectors[:, i] = bond_tf[0]
            bond_vectors[:, j] = bond_tf[1]
            bond_vectors = bond_vectors.reshape(bond_vectors.shape[0], -1)
            return bond_vectors

        def angle_vec(coords, i, j, k):
            ang_tf = nput.angle_deriv(coords, i, j, k)[1]
            ang_vectors = np.zeros(coords.shape)
            ang_vectors[:, i] = ang_tf[0]
            ang_vectors[:, j] = ang_tf[1]
            ang_vectors[:, k] = ang_tf[2]
            ang_vectors = ang_vectors.reshape(ang_vectors.shape[0], -1)
            return ang_vectors

        def poly_kernel(centers, points, m=3):
            vals = 1 + np.sum(centers[:, np.newaxis, :] * points[np.newaxis, :, :], axis=-1)
            return vals ** m

        def poly_krr(coords, vals, l=1, m=3):
            kernel_mat = poly_kernel(coords, coords, m=m)
            regularizer = np.eye(len(coords))
            np.fill_diagonal(regularizer, l)
            weights = np.dot(np.linalg.inv(kernel_mat + regularizer), vals)

            def regression(pts, *, coords=coords, kernel=poly_kernel, m=m, weights=weights):
                return np.dot(weights, kernel(coords, pts, m=m))
            return regression
        if local_coordinates is not None:
            localizers = np.array([bond_vec(cents, *inds) if len(inds) == 2 else angle_vec(cents, *inds) for inds in local_coordinates]).transpose((1, 2, 0))
            pseudos = np.linalg.inv(modes.matrix.T @ modes.matrix) @ modes.matrix.T
            vec_ls = pseudos[np.newaxis] @ localizers[:, :]
            svd_bits = np.linalg.svd(vec_ls)
            polar_rots = svd_bits[0] @ svd_bits[2]
            embedded_coords = nput.vec_tensordot(polar_rots, embedded_coords, shared=1, axes=[2, 1])
            _ = []
            for n, d in enumerate(embedded_vals):
                for __ in range(n):
                    d = nput.vec_tensordot(d, polar_rots, shared=1, axes=[1, 2])
                _.append(d)
            embedded_vals = _
        figs = []
        for crd in range(embedded_coords.shape[-1]):
            rem = np.setdiff1d(np.arange(embedded_coords.shape[-1]), crd)
            inverse_sub_hess = np.linalg.inv(embedded_vals[2][:, rem, :][:, :, rem])
            sub_grad = embedded_vals[1][:, rem]
            pot_contrib = sub_grad[:, np.newaxis, :] @ inverse_sub_hess @ sub_grad[:, :, np.newaxis]
            crd_sort = np.argsort(embedded_coords[:, crd])
            emb_pot = embedded_vals[0] - 1 / 2 * pot_contrib.flatten()
            crds = embedded_coords[crd_sort, crd]
            regression = poly_krr(crds[:, np.newaxis], emb_pot[crd_sort])
            regressed = plt.ScatterPlot(crds, emb_pot[crd_sort] * 219475.6)
            plt.ScatterPlot(crds, regression(crds[:, np.newaxis]) * 219475.6, figure=regressed)
            figs.append(regressed)
            figs.append(plt.ScatterPlot(embedded_coords[crd_sort, crd], embedded_vals[0][crd_sort] * 219475.6))
        figs[0].show()
        raise Exception(...)

    @classmethod
    def getMBPolModel(cls, atoms=None, ref=None, embed=True):
        loader = ModuleLoader(TestManager.current_manager().test_data_dir)
        mbpol = loader.load('LegacyMBPol').MBPol
        b2a = UnitsData.convert('BohrRadius', 'Angstroms')

        def potential(coords, deriv_order=0, chunk_size=int(500000.0)):
            coords = coords.reshape(-1, 9)
            just_vals = deriv_order is None
            if just_vals:
                deriv_order = 0
            chunks = [[] for _ in range(deriv_order + 1)]
            num_chunks = int(len(coords) / chunk_size) + 1
            for coords in np.array_split(coords, num_chunks):
                energies = mbpol.get_pot(coords=coords.reshape(-1, 3, 3) * b2a, nwaters=1, threading_vars=['energy', 'coords'], threading_mode='omp')
                if deriv_order > 0:
                    derivs = []
                    grads = lambda c: mbpol.get_pot_grad(nwaters=1, coords=c.reshape(-1, 3, 3) * b2a, threading_vars=['energy', 'grad', 'coords'], threading_mode='omp')['grad'].reshape(c.shape) * b2a
                    derivs.append(grads(coords))
                    if deriv_order > 1:
                        hess_fun = FiniteDifferenceDerivative(grads, function_shape=(9, 9))
                        new_derivs = hess_fun.derivatives(coords).derivative_tensor(list(range(1, deriv_order)))
                        derivs.extend((np.moveaxis(d, -2, 0) for i, d in enumerate(new_derivs)))
                    for i, d in enumerate([energies] + derivs):
                        chunks[i].append(d)
                else:
                    chunks[0].append(energies)
            for i, c in enumerate(chunks):
                chunks[i] = np.concatenate(c, axis=0)
            if just_vals:
                chunks = chunks[0]
            return chunks
        if ref is None:
            base_mol = Molecule.from_file(TestManager.test_data('water_freq.fchk'))
            ref = base_mol.embed_coords(np.array([[0.0, 0.0656215885, 0.0], [0.757391014, -0.520731105, 0.0], [-0.757391014, -0.520731105, 0.0]]) * UnitsData.convert('Angstroms', 'BohrRadius'))
        if atoms is None:
            atoms = ['O', 'H', 'H']
        ref_mol = Molecule(atoms, ref, internals=[[0, -1, -1, -1], [1, 0, -1, -1], [2, 0, 1, -1]])
        if embed:
            ref_mol = ref_mol.get_embedded_molecule(load_properties=False)
        return (potential, ref_mol)

    @classmethod
    def setupOCHHModel(cls):
        loader = ModuleLoader(os.path.expanduser('~/Documents/Postdoc'))
        h2co_mod = loader.load('H2COPot')
        h2co = h2co_mod.Potential
        h2co_i = h2co_mod.InternalsPotential
        struct = Molecule.from_file(TestManager.test_data('OCHH_freq.fchk'))
        pds = struct.potential_derivatives
        ints = h2co.get_internals(struct.coords * UnitsData.convert('BohrRadius', 'Angstroms'), order=[3, 2, 1, 0])

        def pot_i(c):
            return h2co_i.get_pot(c)
        import scipy.optimize as opt
        opt_ints = opt.minimize(pot_i, ints, method='nelder-mead', tol=1e-08)
        ref_crds = opt_ints.x
        ref_val = opt_ints.fun
        ref_struct = h2co.from_internals(ref_crds)
        ref_struct = struct.embed_coords(ref_struct[(3, 2, 0, 1),] * UnitsData.convert('Angstroms', 'BohrRadius'))
        ref_val = h2co.get_pot(ref_struct * UnitsData.convert('BohrRadius', 'Angstroms'), order=(3, 2, 1, 0))

        def pot(crds, deriv_order=None):
            crds = crds * UnitsData.convert('BohrRadius', 'Angstroms')
            pot_vals = h2co.get_pot(crds, order=(3, 2, 1, 0)) - ref_val
            pot_vals *= UnitsData.convert('Wavenumbers', 'Hartrees')
            if deriv_order is not None and deriv_order > 0:
                from McUtils.Zachary import FiniteDifferenceDerivative
                fd = FiniteDifferenceDerivative(lambda c: h2co.get_pot(c, order=(3, 2, 1, 0)), function_shape=((4, 3), None)).derivatives(crds)
                derivs = fd.derivative_tensor(list(range(1, deriv_order + 1)))
                a2b = UnitsData.convert('BohrRadius', 'Angstroms')
                w2h = UnitsData.convert('Wavenumbers', 'Hartrees')
                for n, d in enumerate(derivs):
                    derivs[n] = np.moveaxis(d * (w2h * a2b ** (n + 1)), -1, 0)
                pot_vals = [pot_vals] + derivs
            return pot_vals
        ref_dip0 = np.array([1.94509054e-10, -1.61826704e-10, -0.928736197])
        ref_dip = np.array([[-0.318, -0.0, -0.0], [-0.0, -0.434, 0.0], [-0.0, 0.0, -0.745], [0.126, 0.0, 0.0], [-0.0, 0.755, 0.0], [-0.0, -0.0, 0.874], [0.096, -0.0, 0.0], [-0.0, -0.161, 0.144], [-0.0, 0.083, -0.065], [0.096, -0.0, -0.0], [0.0, -0.161, -0.144], [0.0, -0.083, -0.065]])

        def dip(crds, deriv_order=None, ref_dip=ref_dip):
            crds = mol.embed_coords(crds)
            dx = (crds - mol.coords[np.newaxis]).reshape(-1, 12)
            mu = ref_dip0[np.newaxis] + dx @ ref_dip
            if deriv_order is not None:
                vals = [mu]
                if deriv_order > 0:
                    der = ref_dip
                    for _ in crds.shape[:-2]:
                        der = np.expand_dims(der, 0)
                    der = np.broadcast_to(der, crds.shape[:-2] + der.shape)
                    vals.append(der)
                if deriv_order > 1:
                    for n in range(2, deriv_order + 1):
                        vals.append(np.zeros(crds.shape[:-2] + (12,) * n + (3,)))
            else:
                vals = mu
            return vals
            vOH = crds[..., 0, :] - crds[..., 1, :]
            vCH1 = crds[..., 2, :] - crds[..., 1, :]
            vCH2 = crds[..., 3, :] - crds[..., 1, :]
            mCH1 = muCH * vCH1
            mCH2 = muCH * vCH2
            muCH1 = muCH * np.diag(nput.vec_normalize(vCH1))
            muCH2 = muCH * np.diag(nput.vec_normalize(vCH2))
            mOH = muOH * vOH
            mu = mCH1 + mCH2 + mOH
            muuOH = mOH * nput.vec_normalize(vOH)
            if deriv_order is not None:
                vals = [mu]
                if deriv_order > 0:
                    for _ in crds.shape[:-2]:
                        der = np.expand_dims(der, 0)
                    der = np.broadcast_to(der, crds.shape[:-2] + der.shape)
                    vals.append(der)
                if deriv_order > 1:
                    for n in range(2, deriv_order + 1):
                        vals.append(np.zeros(crds.shape[:-2] + (12,) * n + (3,)))
            else:
                vals = mu
            return vals
        mol = Molecule(['O', 'C', 'H', 'H'], ref_struct)
        return (pot, dip, mol)
    default_OH_freq = 3869.47
    default_OH_anh = 84
    default_HOH_freq = 1624

    @classmethod
    def setup_Water(cls, use_mbpol=False, atoms=None, oh_model=False, stretch_model=False, bend_model=False, full_model=True, w=None, wx=None, w2=None, wx2=None, w_bend=None, freq_units='Wavenumbers', dipole_model='linear', potential_params=None, reoptimize=True):
        mol = Molecule.from_file(TestManager.test_data('water_freq.fchk'))
        if atoms is not None:
            mol = mol.modify(atoms=atoms).get_embedded_molecule()
        h2w = UnitsData.convert(freq_units, 'Hartrees')
        base_params = potential_params
        potential_params = {} if potential_params is None else potential_params
        if oh_model:
            mol = mol.modify(atoms=mol.atoms[:2], coords=mol.coords[:2], dipole_evaluator='expansion', dipole_derivatives=[mol.dipole_derivatives[0], mol.dipole_derivatives[1][:6, :]]).get_embedded_molecule(load_properties=False)
            if w is None:
                w = cls.default_OH_freq
            if wx is None:
                wx = cls.default_OH_anh
        elif stretch_model:
            if w is None:
                w = cls.default_OH_freq
            if wx is None:
                wx = cls.default_OH_anh
            if w2 is None:
                w2 = cls.default_OH_freq
            if wx2 is None:
                wx2 = cls.default_OH_anh
        elif full_model:
            if w is None:
                w = cls.default_OH_freq
            if wx is None:
                wx = cls.default_OH_anh
            if w2 is None:
                w2 = cls.default_OH_freq
            if wx2 is None:
                wx2 = cls.default_OH_anh
            if w_bend is None:
                w_bend = cls.default_HOH_freq
        if not use_mbpol and base_params is None:
            if w is not None:
                if wx is None:
                    potential_params[0, 1] = {'w_coeffs': [w * h2w]}
                else:
                    potential_params[0, 1] = {'w': w * h2w, 'wx': wx * h2w}
            if w2 is not None:
                if wx2 is None:
                    potential_params[0, 2] = {'w_coeffs': [w2 * h2w]}
                else:
                    potential_params[0, 2] = {'w': w2 * h2w, 'wx': wx2 * h2w}
            if w_bend is not None:
                potential_params[1, 0, 2] = {'w_coeffs': [w_bend * h2w]}
        if len(potential_params) > 0:
            pot = mol.get_1d_potentials(potential_params)
            base_pot = sum(pot)

            def potential(coords, order=None, base_pot=base_pot):
                none_ord = order is None
                if none_ord:
                    order = 0
                pot_exp = base_pot(coords, order=order)[1]
                if none_ord:
                    pot_exp = pot_exp[0]
                return pot_exp
            mol = mol.modify(energy_evaluator={'potential_function': potential, 'analytic_derivative_order': 6, 'batched_orders': True}, dipole_derivatives=mol.dipole_derivatives[:2] if dev.str_is(dipole_model, 'linear') else mol.dipole_derivatives, dipole_evaluator='expansion')
        else:
            loader = ModuleLoader(TestManager.current_manager().test_data_dir)
            mbpol = loader.load('LegacyMBPol').MBPol

            def potential(coords, order=None, chunk_size=int(500000.0)):
                coords = np.asanyarray(coords)
                base_shape = coords.shape[:-1] if coords.shape[-1] == 9 else coords.shape[:-2]
                coords = coords.reshape(-1, 9)
                just_vals = order is None
                if just_vals:
                    order = 0
                chunks = [[] for _ in range(order + 1)]
                num_chunks = int(len(coords) / chunk_size) + 1
                for coords in np.array_split(coords, num_chunks):
                    if order == 0:
                        energies = mbpol.get_pot(coords=coords.reshape(-1, 3, 3), nwaters=1, threading_vars=['energy', 'coords'], threading_mode='omp')
                        chunks[0].append(energies)
                    else:
                        grad_vals = mbpol.get_pot_grad(nwaters=1, coords=coords.reshape(-1, 3, 3), threading_vars=['energy', 'grad', 'coords'], threading_mode='omp')
                        chunks[0].append(grad_vals['energy'])
                        chunks[1].append(grad_vals['grad'].reshape((-1, 9)))
                chunks = [np.concatenate(c, axis=0) for c in chunks]
                chunks = [c.reshape(base_shape + c.shape[1:]) for c in chunks]
                if just_vals:
                    chunks = chunks[0]
                return chunks
            mol = mol.modify(coords=np.array([[0.0, 0.12400683, 0.0], [1.43126159, -0.98403917, 0.0], [-1.43126159, -0.98403917, 0.0]]), energy_evaluator={'potential_function': potential, 'analytic_derivative_order': 1, 'distance_units': 'Angstroms', 'energy_units': 'Hartrees', 'batched_orders': True, 'mesh_spacing': 0.0001}, dipole_derivatives=mol.dipole_derivatives, dipole_evaluator='expansion').get_embedded_molecule(load_properties=False)
            if reoptimize:
                opt_mol = mol.optimize()
                mol = opt_mol
        return mol

    @classmethod
    def setup_OCHH(cls, optimize=False, use_internals=True, load_potential=True, embed=True, **fd_opts):
        if load_potential:
            loader = ModuleLoader(os.path.expanduser('~/Documents/Postdoc/Projects/DGB'))
            h2co_mod = loader.load('H2COPot')

        def pot(coords, order=None):
            if coords.shape[-1] == 12:
                base_shape = coords.shape[:-1]
            else:
                base_shape = coords.shape[:-2]
            vals = h2co_mod.Potential.get_pot(coords.reshape(-1, 4, 3), order=(3, 2, 1, 0))
            return vals.reshape(base_shape)
        if not use_internals:
            ochh = Molecule.from_file(TestManager.test_data('OCHH_freq.fchk'), energy_evaluator=dict({'potential_function': pot, 'distance_units': 'Angstroms', 'energy_units': 'Wavenumbers'}, **fd_opts))
        else:

            def internal_pot(coords, order=None):
                coords = coords[..., (0, 1, 3, 2, 4, 5)]
                vals = h2co_mod.InternalsPotential.get_pot(coords)
                return vals
            ochh = Molecule.from_file(TestManager.test_data('OCHH_freq.fchk'), energy_evaluator=dict({'potential_function': internal_pot, 'distance_units': 'Angstroms', 'energy_units': 'Wavenumbers', 'strip_embedding': True, 'supports_internals': True}, **fd_opts), internals=[[0, -1, -1, -1], [1, 0, -1, -1], [2, 1, 0, -1], [3, 1, 0, 2]])
        base_dip = ochh.dipole_derivatives[:2]
        ochh = ochh.modify(coords=[[0, 0, 1.27603352], [0, 0, -0.998625259], [0, 1.77174246, -2.09630576], [0, -1.77174246, -2.09630576]], dipole_derivatives=base_dip, dipole_evaluator='expansion')
        if embed:
            ochh = ochh.get_embedded_molecule(embed_properties=True)
        if optimize:
            ochh = ochh.optimize(method='conjugate-gradient', max_iterations=50, stencil=3, logger=True, prevent_oscillations=3, restart_interval=15)
            print(ochh.coords)
        return ochh

    @validationTest
    def test_ModelPotentialAIMD2D(self):
        mol, model = self.buildWaterModel(ka=None, dudr1=1 / 5.5, dudr2=1 / 5.5)
        check_freqs = False
        if check_freqs:
            freqs = model.normal_modes()[0]
            raise Exception(freqs * UnitsData.convert('Hartrees', 'Wavenumbers'))
        check_anh = False
        if check_anh:
            from Psience.VPT2 import VPTRunner
            VPTRunner.run_simple([mol.atoms, mol.coords], potential_derivatives=model.potential(mol.coords, deriv_order=4)[1:], dipole_derivatives=model.dipole(mol.coords, deriv_order=3), order=2, states=3, logger=True, degeneracy_specs='auto', calculate_intensities=True)
            '\n            ZPE: 3869.37229   3814.75070 \n            State     Frequency    Intensity       Frequency    Intensity\n              0 1    3896.87028     64.98650      3726.72077     63.58462\n              1 0    3841.87433      0.26781      3675.96177      0.25708\n              0 2    7793.74057      0.00000      7415.60644      0.00317\n              2 0    7683.74866      0.00000      7220.11984      0.00861\n              1 1    7738.74461      0.00000      7236.25926      1.32009\n              0 3   11690.61085      0.00000     10986.10686      0.00048\n              3 0   11525.62299      0.00000     10895.67849      0.00008\n              1 2   11635.61490      0.00000     10593.72292      0.00056\n              2 1   11580.61894      0.00000     10596.33867      0.05435\n            '
            raise Exception(...)
        '\n        >>------------------------- Running distributed Gaussian basis calculation -------------------------\n        :: diagonalizing in the space of 69 S functions\n        :: ZPE: 3815.1602161487854\n        :: Frequencies: [ 3676.15271753  3727.03078901  7220.80144039  7236.73919476  7416.68728081 10595.20948421 10597.32265836 10897.3034496  10990.64643556 13795.47060442 13797.44564618 14285.28218497 14321.1505273  14573.69425921 16830.3614471  16835.00173818 17518.48942695 17543.56990564 17847.28555492 18095.84541425]\n        >>--------------------------------------------------<<\n        '
        check_dvr = False
        if check_dvr:
            raise Exception('getting model DVRs working in normal mode subspaces unfinished')
            nm_model, freqs = model.to_normal_modes()
            raise Exception(nm_model.coords, nm_model.rotation)
            dvr = nm_model.setup_DVR(domain=[[-10, 10], [-10, 10]], divs=[65, 65], logger=True)
            po_data = dvr.run()
            raise Exception(po_data.wavefunctions.energies[0] * UnitsData.hartrees_to_wavenumbers, po_data.wavefunctions.frequencies() * UnitsData.hartrees_to_wavenumbers)
        symm, asymm = model.normal_modes()[0]
        sim = model.setup_AIMD(initial_energies=np.array([[symm, asymm], [-symm, asymm], [0, -asymm], [-symm, 0]]) * 2, timestep=15, track_velocities=True)
        sim.propagate(25)
        coords, velocities = sim.extract_trajectory(flatten=True, embed=mol.coords)
        momenta = 1872 * velocities * mol.masses[np.newaxis, :, np.newaxis]
        cartesians = False
        with BlockProfiler(inactive=True):
            dgb = model.setup_DGB(coords, optimize_centers=1e-14, modes=None if cartesians else 'normal', cartesians=[0, 1] if cartesians else None, quadrature_degree=3, kinetic_options=dict(include_diagonal_contribution=True, include_coriolis_coupling=True, include_watson_term=True), expansion_degree=2, pairwise_potential_functions={'functions': {(0, 1): self.setupMorseFunction(model, 0, 1), (0, 2): self.setupMorseFunction(model, 0, 2)}, 'quadrature_degree': 7})
            '\n            >>------------------------- Running distributed Gaussian basis calculation -------------------------\n            :: diagonalizing in the space of 80 S functions\n            :: ZPE: 3830.3764717334275\n            :: Frequencies: [ 3675.84806189  3726.63119184  7220.21030118  7236.04910246  7416.59927689 10594.35262638 10596.41056233 10898.18829694 10992.57836132 13795.01411084 13795.28408024 14287.32637531 14325.71150395 14607.27796273 16827.61503786 16828.53168599 17521.47359305 17552.51517707 17879.83078742 18202.01573779]\n            >>--------------------------------------------------<<\n            '
            plot_gaussians = False
            if plot_gaussians:
                self.plot_gaussians(dgb, mol)
                raise Exception(...)
            '\n            >>------------------------- Running distributed Gaussian basis calculation -------------------------\n            :: diagonalizing in the space of 59 S functions\n            :: ZPE: 3815.0254181251166\n            :: Frequencies: [ 3676.14207931  3726.99361301  7221.08377964  7236.78275062  7418.52358929 10595.91399038 10597.43855691 10901.20911454 11001.80376921 13797.34317942 13799.82687628 14293.18090478 14329.25062911 14640.80903924 16841.23176496 16849.01710935 17531.49955265 17572.95249918 17891.12454688 18319.17520319]\n            >>--------------------------------------------------<<\n            '
            '\n            >>------------------------- Running distributed Gaussian basis calculation -------------------------\n            :: diagonalizing in the space of 41 S functions\n            :: ZPE: 3815.1110786450845\n            :: Frequencies: [ 3676.77263767  3727.70483878  7227.82029744  7237.90664201  7440.76152723 10601.58120583 10618.10732722 10920.95275887 11104.66728657 13808.26706889 13827.90646153 14416.22589799 14471.29425293 15204.74679781 16858.16199974 16902.38702413 17668.82463041 17864.53996371 18809.9243302  18939.68118366]\n            >>--------------------------------------------------<<\n            '
            self.runDGB(dgb, mol, similarity_cutoff=0.9, plot_spectrum=False, plot_wavefunctions=True)
