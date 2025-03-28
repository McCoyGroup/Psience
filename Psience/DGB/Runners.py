
import numpy as np, os
from McUtils.Data import AtomData, UnitsData
import McUtils.Plots as plt
import McUtils.Numputils as nput

from ..Molecools import Molecule

from .DGB import DGB
from .Coordinates import DGBCoords, DGBCartesians
from .Wavefunctions import DGBWavefunctions, DGBWavefunction

__all__ = [
    "DGBRunner"
]

__reload_hook__ = ['.DGB', '.Coordinates', '.Wavefunctions']

class DGBRunner:

    # @classmethod
    # def prep_pairwise_functions(cls, coord, ...):

    @classmethod
    def construct_from_mol_simulation(cls,
                                      sim, mol,
                                      *,
                                      potential_function=None,
                                      dipole_function=None,
                                      use_cartesians=False,
                                      use_momenta=False,
                                      quadrature_degree=3,
                                      expansion_degree=2,
                                      use_interpolation=True,
                                      use_quadrature=False,
                                      symmetrizations=None,
                                      momentum_scaling=None,
                                      skip_initial_configurations=True,
                                      modes='normal',
                                      transformations='diag',
                                      logger=True,
                                      **opts
                                      ):

        coords, velocities = sim.extract_trajectory(flatten=True, embed=mol.coords)
        if skip_initial_configurations:
            nconf = len(sim.coords)
            coords = coords[nconf - 1:]
            velocities = velocities[nconf - 1:]
        if use_momenta:
            momenta = velocities * mol.masses[np.newaxis, :, np.newaxis]
            if momentum_scaling is not None:
                momenta = momentum_scaling * momenta
        else:
            momenta = None

        if potential_function is not None:
            mol = mol.modify(energy_evaluator=potential_function)
            mol.potential_derivatives = mol.calculate_energy(mol.coords, order=2)[1:]
        else:
            mol = mol.modify(
                potential_derivatives=mol.calculate_energy(mol.coords, order=2)[1:]
            )
        if dipole_function is not None:
            mol = mol.modify(dipole_evaluator=dipole_function)
        potential_function = mol.get_energy_function()
        nms = mol.get_normal_modes()

        if use_interpolation:
            if symmetrizations is not None:
                emb_coords = nms.embed_coords(coords)
                embs = [coords]
                for symm in symmetrizations:
                    if nput.is_numeric(symm) or all(s >= 0 for s in symm):
                        s = [1] * len(nms.freqs)
                        for i in symm: s[i] = -1
                        symm = s
                    new_emb = emb_coords * np.array(symm)[np.newaxis, :]
                    new_coords = nms.unembed_coords(new_emb)
                    embs.append(new_coords)
                interp_coords = np.concatenate(embs, axis=0)
            else:
                interp_coords = coords
            potential_data = {
                'centers': interp_coords,
                'values': mol.calculate_energy(interp_coords, order=2)
            }
        else:
            potential_data = potential_function

        if use_quadrature and use_interpolation:
            raise ValueError("don't use interpolation with quadrature...")

        if use_cartesians:
            diffs = [np.ptp(coords[:, i]) for i in range(3)]
            axes = [i for i, d in enumerate(diffs) if d > 1e-8]
        else:
            axes = None

        return DGB.construct(
            coords,
            **dict(
                dict(
                    potential_function=potential_data,
                    dipole_function=mol.get_dipole_function(),
                    modes=None if use_cartesians else modes,
                    cartesians=axes
                    , quadrature_degree=quadrature_degree
                    , expansion_degree=expansion_degree if not use_quadrature else None
                    , transformations=transformations
                    , momenta=momenta
                    , logger=logger
                ),
                **opts
            )
        )

    @classmethod
    def construct_from_model_simulation(cls,
                                        sim, model, mol=None,
                                        *,
                                        use_cartesians=False,
                                        use_momenta=False,
                                        quadrature_degree=3,
                                        expansion_degree=2,
                                        use_interpolation=True,
                                        use_quadrature=False,
                                        symmetrizations=None,
                                        momentum_scaling=None,
                                        skip_initial_configurations=True,
                                        modes='normal',
                                        transformations='diag',
                                        **opts
                                        ):
        if mol is None:
            mol = model.mol.get_embedded_molecule()

        coords, velocities = sim.extract_trajectory(flatten=True, embed=mol.coords)
        if skip_initial_configurations:
            nconf = len(sim.coords)
            coords = coords[nconf - 1:]
            velocities = velocities[nconf - 1:]
        if use_momenta:
            momenta = velocities * mol.masses[np.newaxis, :, np.newaxis]
            if momentum_scaling is not None:
                momenta = momentum_scaling * momenta
        else:
            momenta = None

        mol.potential_derivatives = model.potential(mol.coords, deriv_order=2)[1:]
        nms = mol.normal_modes.modes.basis.to_new_modes()

        if use_interpolation:
            if symmetrizations is not None:

                emb_coords = nms.embed_coords(coords)
                embs = [coords]
                for symm in symmetrizations:
                    if nput.is_numeric(symm) or all(s >= 0 for s in symm):
                        s = [1]*len(nms.freqs)
                        for i in symm: s[i] = -1
                        symm = s
                    new_emb = emb_coords * np.array(symm)[np.newaxis, :]
                    new_coords = nms.unembed_coords(new_emb)
                    embs.append(new_coords)
                interp_coords = np.concatenate(embs, axis=0)
            else:
                interp_coords = coords
            potential_data = {
                'centers': interp_coords,
                'values': model.potential(interp_coords, deriv_order=2)
            }
        else:
            potential_data = None

        if use_quadrature and use_interpolation:
            raise ValueError("don't use interpolation with quadrature...")

        if use_cartesians:
            diffs = [np.ptp(coords[:, i]) for i in range(3)]
            axes = [i for i,d in enumerate(diffs) if d > 1e-8]
        else:
            axes = None

        return model.setup_DGB(
            coords,
            **dict(
                dict(
                    potential_function=potential_data,
                    modes=None if use_cartesians else modes,
                    cartesians=axes
                    , quadrature_degree=quadrature_degree
                    , expansion_degree=expansion_degree if not use_quadrature else None
                    , transformations=transformations
                    , momenta=momenta
                ),
                **opts
            )
        )

    @classmethod
    def construct_from_model(cls,
                             model,
                             trajectories=10,
                             *,
                             sim=None,
                             propagation_time=10,
                             timestep=50,
                             use_cartesians=False,
                             use_momenta=False,
                             use_pairwise=True,
                             use_interpolation=True,
                             use_quadrature=False,
                             symmetrizations=None,
                             momentum_scaling=None,
                             track_velocities=True,
                             **aimd_options
                             ):

        if sim is None:
            sim = model.setup_AIMD(
                trajectories=trajectories,
                timestep=timestep,
                track_velocities=True,
                **aimd_options
            )
            sim.propagate(propagation_time)

        return cls.construct_from_model_simulation(
            sim, model,
            use_cartesians=use_cartesians,
            use_momenta=use_momenta,
            use_pairwise=use_pairwise,
            use_interpolation=use_interpolation,
            use_quadrature=use_quadrature,
            momentum_scaling=momentum_scaling,
            symmetrizations=symmetrizations,
        )

    @classmethod
    def from_mol(cls,
                 mol, sim=None,
                 *,
                 potential_function=None,
                 dipole_function=None,
                 trajectories=10,
                 propagation_time=10,
                 timestep=50,
                 use_cartesians=False,
                 use_momenta=False,
                 use_pairwise=True,
                 use_interpolation=True,
                 use_quadrature=False,
                 symmetrizations=None,
                 momentum_scaling=None,
                 track_velocities=True,
                 **aimd_options
                 ):

        if potential_function is not None:
            mol = mol.modify(energy_evaluator=potential_function)
        if sim is None:
            sim = mol.setup_AIMD(
                trajectories=trajectories,
                timestep=timestep,
                track_velocities=True,
                **aimd_options
            )
            sim.propagate(propagation_time)

        return cls.construct_from_mol_simulation(
            sim, mol,
            potential_function=None,
            use_cartesians=use_cartesians,
            use_momenta=use_momenta,
            use_pairwise=use_pairwise,
            use_interpolation=use_interpolation,
            use_quadrature=use_quadrature,
            symmetrizations=symmetrizations,
            momentum_scaling=momentum_scaling,
            dipole_function=dipole_function,
        )

    # @classmethod
    # def from_coord(cls,
    #             coords,
    #             potential=None,
    #             dipole=None,
    #             *,
    #             logger=True,
    #             plot_wavefunctions=True,
    #             plot_spectrum=True,
    #             **opts
    #     ):
    #
    #         dgb = DGB.construct(
    #             np.round(coords, 8),  # this ends up really mattering to keep optimize_centers stable
    #             pot,
    #             logger=logger,
    #             **opts
    #         )

    @classmethod
    def run_dgb(cls,
                dgb,
                *,
                plot_wavefunctions=True,
                plot_spectrum=True
                ):

            logger = dgb.logger
            with logger.block(tag="Running DGB"):
                logger.log_print("num coords: {c}", c=len(dgb.gaussians.coords.centers))
                with logger.block(tag="S"):
                    logger.log_print(logger.prep_array(dgb.S[:5, :5]))
                with logger.block(tag="T"):
                    logger.log_print(logger.prep_array(dgb.T[:5, :5]))
                with logger.block(tag="V"):
                    logger.log_print(logger.prep_array(dgb.V[:5, :5]))

                wfns_cart = dgb.get_wavefunctions()
                with logger.block(tag="Energies"):
                    logger.log_print(
                        logger.prep_array(wfns_cart.energies[:5] * UnitsData.convert("Hartrees", "Wavenumbers"))
                    )
                with logger.block(tag="Frequencies"):
                    logger.log_print(
                        logger.prep_array(wfns_cart.frequencies()[:5] * UnitsData.convert("Hartrees", "Wavenumbers"))
                    )

                if plot_wavefunctions:
                    for i in range(4):
                        wfns_cart[i].plot().show()

                if plot_spectrum:
                    spec = wfns_cart[:4].get_spectrum()
                    with logger.block(tag="Intensities"):
                        logger.log_print(
                            logger.prep_array(spec.intensities)
                        )
                    spec.plot().show()
                else:
                    spec = None
            return wfns_cart, spec
    @classmethod
    def run_simple(cls,
                   system_spec,
                   sim=None,
                   plot_wavefunctions=True,
                   plot_spectrum=True,
                   **opts
                   ):
        if isinstance(system_spec, (str, dict, tuple)):
            system_spec = Molecule.construct(system_spec)
            system_spec.get_dipole_evaluator()
        if hasattr(system_spec, 'potential'): # analyt
            dgb = cls.construct_from_model(system_spec, sim=sim, **opts)
        else:
            dgb = cls.from_mol(system_spec, sim=sim, **opts)

        return dgb, cls.run_dgb(
            dgb,
            plot_wavefunctions=plot_wavefunctions,
            plot_spectrum=plot_spectrum
        )

    @classmethod
    def plot_dgb_potential(cls,
                           dgb, mol, potential,
                           coordinate_sel=None,
                           domain=None, domain_padding=1,
                           potential_cutoff=17000,
                           potential_min=0,
                           plot_cartesians=None,
                           plot_atoms=True,
                           cmap=None,
                           modes_nearest=False,
                           plot_points=100,
                           levels=24,
                           **plot_styles
                           ):
        def cutoff_pot(points, cutoff=potential_cutoff / UnitsData.hartrees_to_wavenumbers,
                       cutmin=potential_min / UnitsData.hartrees_to_wavenumbers
                       ):
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
            figure = mol.plot_molecule_function(
                cutoff_pot,
                axes=coordinate_sel,
                domain=domain,
                modes_nearest=modes_nearest,
                domain_padding=domain_padding,
                cmap=cmap,
                plot_points=plot_points,
                levels=levels,
                plot_atoms=plot_atoms,
                **plot_styles
            )
        else:
            if domain is None:
                from McUtils.Zachary import Mesh
                domain = Mesh(coords.centers).bounding_box
            points = DGBWavefunction.prep_plot_grid(
                domain,
                domain_padding=domain_padding,
            )
            if points.shape[-1] == 1:
                figure = plt.Plot(
                    *np.moveaxis(points, -1, 0),
                    potential(points),
                    **plot_styles
                )
            else:
                figure = plt.TriContourPlot(
                    *np.moveaxis(points, -1, 0),
                    cutoff_pot(points),
                    cmap=cmap,
                    levels=levels,
                    **plot_styles
                )

        return figure

    @classmethod
    def plot_gaussians(cls,
                       dgb, mol,
                       *,
                       domain=None,
                       domain_padding=1,
                       cmap='RdBu',
                       plot_dir=None,
                       plot_name='gaussian_{i}.pdf',
                       **plot_options
                       ):
        n = len(dgb.gaussians.coords.centers)
        wfns = DGBWavefunctions(
            np.zeros(n),
            np.eye(n),
            dgb
        )
        cls.plot_wavefunctions(
            wfns, dgb, mol,
            cmap=cmap,
            plot_name=plot_name,
            plot_dir=plot_dir,
            domain=domain,
            domain_padding=domain_padding,
            potential_styles=dict(
                domain=domain,
                domain_padding=domain_padding
            ),
            **plot_options
        )

    default_num_plot_wfns = 5
    @classmethod
    def plot_wavefunctions(cls,
                           wfns, dgb, mol,
                           which=True,
                           coordinate_sel=None,
                           cartesians=None,
                           plot_dir=None,
                           plot_name='wfn_{i}.pdf',
                           plot_label='{e} cm-1',
                           plot_potential=True,
                           plot_atoms=None,
                           plot_centers=True,
                           potential_styles=None,
                           **plot_options
                           ):

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
                pot_figure = cls.plot_dgb_potential(
                    dgb, mol, pot,
                    coordinate_sel=coordinate_sel,
                    plot_atoms=plot_atoms,
                    **potential_styles
                )
            else:
                pot_figure = None

            if plot_dir is not None:
                os.makedirs(plot_dir, exist_ok=True)
            if (
                    coordinate_sel is None and wfns.gaussians.alphas.shape[-1] == 1
                    or len(coordinate_sel) == 1
            ):
                for k in ['cmap', 'levels', 'plotter', 'contour_levels']:
                    if k in plot_options: del plot_options[k]
                    if k in potential_styles: del potential_styles[k]
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

                        figs.append(
                            w.plot_cartesians(
                                coordinate_sel,
                                plot_centers=plot_centers,
                                figure=figure,
                                plot_label=lab,
                                **plot_options
                            )
                        )
                    else:
                        if wfns.gaussians.alphas.shape[-1] > 1:
                            if coordinate_sel is not None:
                                w = w.project(coordinate_sel)
                            figs.append(
                                w.plot(
                                    figure=figure,
                                    plot_centers=plot_centers,
                                    plot_label=lab,
                                    **plot_options
                                )
                            )
                        else:
                            figs.append(
                                w.plot(
                                    figure=figure,
                                    plot_centers=plot_centers,
                                    plot_label=lab,
                                    **plot_options
                                )
                            )

                    if plot_dir is not None:
                        fig = figs.pop()
                        fig.savefig(os.path.join(plot_dir, plot_name.format(i=i)))
                        fig.close()

            if plot_dir is None:
                figs[0].show()

    @classmethod
    def run_DGB(cls,
               dgb: DGB,
               mol,
               plot_centers=True,
               plot_wavefunctions=True,
               plot_spectrum=False,
               pot_cmap='viridis',
               wfn_cmap='RdBu',
               wfn_points=100,
               wfn_contours=12,
               plot_dir=None,
               plot_potential=True,
               pot_points=100,
               domain=None,
               domain_padding=1,
               potential_cutoff=15000,
               mode=None,
               nodeless_ground_state=None,
               min_singular_value=None,
               subspace_size=None,
               plot_similarity=False,
               similarity_cutoff=None,
               similarity_chunk_size=None,
               similar_det_cutoff=None,
               **plot_options
               ):

        print("--->", len(dgb.gaussians.coords.centers))
        print(dgb.S[:5, :5])
        print(dgb.T[:5, :5])
        print(dgb.V[:5, :5])

        try:
            wfns, spec = dgb.run(
                calculate_spectrum=plot_spectrum,
                mode=mode,
                nodeless_ground_state=nodeless_ground_state,
                subspace_size=subspace_size,
                min_singular_value=min_singular_value,
                similarity_cutoff=similarity_cutoff,
                similarity_chunk_size=similarity_chunk_size,
                similar_det_cutoff=similar_det_cutoff
            )
        except Exception as e:

            if plot_wavefunctions is not False:
                print(e)

                if isinstance(plot_wavefunctions, str) and plot_wavefunctions == 'cartesians':
                    plot_wavefunctions = {'cartesians': None}
                elif plot_wavefunctions is True and dgb.gaussians.coords.centers.shape[-1] > 2:
                    plot_wavefunctions = {'cartesians': None}

                if plot_wavefunctions.get('cartesians', False) is None:
                    embed_coords = dgb.gaussians.coords.as_cartesians()[0].centers
                    diffs = [np.ptp(embed_coords[:, i]) for i in range(3)]
                    axes = [i for i, d in enumerate(diffs) if d > 1e-8]
                    plot_wavefunctions['cartesians'] = axes

                cartesian_plot_axes = None
                if isinstance(plot_wavefunctions, dict):
                    if 'cartesians' in plot_wavefunctions:
                        cartesian_plot_axes = plot_wavefunctions['cartesians']
                        dgb = dgb.as_cartesian_dgb()
                    else:
                        raise ValueError(plot_wavefunctions)

                pot = dgb.pot.potential_function
                figure = cls.plot_dgb_potential(
                    dgb, mol, pot,
                    domain=domain, domain_padding=domain_padding,
                    potential_cutoff=potential_cutoff,
                )

                dgb.gaussians.plot_centers(
                    figure,
                    xyz_sel=cartesian_plot_axes
                )

                figure.show()

            raise e
        else:
            if plot_similarity:
                plt.ArrayPlot(dgb.get_similarity_matrix()).show()


            if isinstance(plot_wavefunctions, str) and plot_wavefunctions == 'cartesians':
                plot_wavefunctions = {'cartesians': None}
            elif plot_wavefunctions is True and dgb.gaussians.coords.centers.shape[-1] > 2:
                plot_wavefunctions = {'cartesians': None}

            if plot_wavefunctions.get('cartesians', False) is None:
                embed_coords = dgb.gaussians.coords.as_cartesians()[0].centers
                diffs = [np.ptp(embed_coords[:, i]) for i in range(3)]
                axes = [i for i, d in enumerate(diffs) if d > 1e-8]
                plot_wavefunctions['cartesians'] = axes

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
                # which = plot_wavefunctions
                # if which is True:
                #     which = cls.default_num_plot_wfns
                # if isinstance(which, int):
                #     which = list(range(which))
                # print(which, spec.intensities, spec.frequencies, spec[which])
                if plot_dir is not None:
                    spec_plot = spec.plot()
                    os.makedirs(plot_dir, exist_ok=True)
                    spec_plot.savefig(os.path.join(plot_dir, 'spec.png'))
                    spec_plot.close()

            if plot_wavefunctions:
                cls.plot_wavefunctions(
                    wfns,
                    dgb,
                    mol,
                    which=plot_wavefunctions,
                    cartesians=use_cartesians,
                    coordinate_sel=coordinate_sel,
                    plot_dir=plot_dir,
                    contour_levels=wfn_contours,
                    cmap=wfn_cmap,
                    plot_points=wfn_points,
                    plot_centers={'color': 'red'} if plot_centers else False,
                    domain=domain,
                    domain_padding=domain_padding,
                    plot_potential=plot_potential,
                    scaling=.2,
                    plotter=plt.TriContourLinesPlot,
                    potential_styles=dict(
                        domain=domain,
                        domain_padding=domain_padding,
                        cmap=pot_cmap,
                        plot_points=pot_points
                    ),
                    **plot_options
                )

            if plot_spectrum:
                return wfns, spec
            else:
                return wfns


    @classmethod
    def getMorseParameters(cls, w=None, wx=None, m1=None, m2=None, re=None):
        if w is None:
            w = 3869.47 / UnitsData.hartrees_to_wavenumbers
        freq = w
        if wx is None:
            wx = 84 / UnitsData.hartrees_to_wavenumbers
        anh = wx
        De = (freq ** 2) / (4 * anh)
        if m1 is None:
            m1 = AtomData["O", "Mass"] * UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
        if m2 is None:
            m2 = AtomData["H", "Mass"] * UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
        muv = (1 / m1 + 1 / m2)
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

        def morse_basic(r,
                        re=re,
                        alpha=a,
                        De=De,
                        deriv_order=None,
                        _morse=PotentialData["MorsePotential"]
                        ):
            return _morse(r, re=re, alpha=alpha, De=De, deriv_order=deriv_order)

        return morse_basic

    @classmethod
    def plot_interpolation_error(cls, dgb, pot):
        sel = slice(None)  # slice(15,30)

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
        # dev_hess = inter_hess - real_hess
        grad_plot = plt.ScatterPlot(realpots[ords], dev_gnorm)

        inter_trace = np.sum(np.sum(np.abs(inter_hess[ords]), axis=-1), axis=-1)
        real_trace = np.sum(np.sum(np.abs(real_hess[ords]), axis=-1), axis=-1)
        dev_trace = inter_trace - real_trace
        dev_hess = inter_hess - real_hess
        hess_plot = plt.ScatterPlot(realpots[ords], dev_trace)

        print("Mean Absolute Error:", np.mean(np.abs(unscaled_devs)), "Std:", np.std(unscaled_devs))
        print("Mean Scaled Error:", np.mean(np.abs(devs)), "Std:", np.std(devs))
        print("Mean Hessian Error:", np.mean(np.abs(dev_hess.flatten())),
              "Std:", np.std(dev_hess.flatten()),
              "Max:", np.max(np.abs(dev_hess.flatten()))
              )
        print("Mean Summed Hessian Error:", np.mean(np.abs(dev_trace.flatten())),
              "Std:", np.std(dev_trace.flatten()),
              "Max:", np.max(np.abs(dev_trace.flatten())),
              )
        print("Maximum (Scaled) Error:", devs[max_dev_pos])
        print("Maximum (Scaled) Interpolation Error:")
        for l, r, c, tt, ii, ov in zip(
                dgb.gaussians.coords.centers[rows[sel][ords[max_dev_pos]]],
                dgb.gaussians.coords.centers[cols[sel][ords[max_dev_pos]]],
                dgb.gaussians.overlap_data.centers[sel][ords[max_dev_pos]],
                realpots[ords[max_dev_pos]],
                interpots[ords[max_dev_pos]],
                utris[sel][ords[max_dev_pos]]
        ):
            print(f"Centers: {c} ({ov}) <- {l} {r}")
            print(f"  Error: {ii - tt} <- {tt} {ii}")

        bad_bad = np.abs(devs) > 50

        # center_plot=plt.ScatterPlot(
        #     centers[:, 0],
        #     centers[:, 1],
        #     c=unscaled_devs
        # )
        # center_plot = plt.ScatterPlot(
        #     centers[:, 0][bad_bad],
        #     centers[:, 1][bad_bad],
        #     c='blue'
        # )
        # plt.ScatterPlot(
        #     dgb.gaussians.coords.centers[:, 0],
        #     dgb.gaussians.coords.centers[:, 1],
        #     c='red',
        #     figure=center_plot
        # )
        # woof = plt.ScatterPlot(realpots[ords], unscaled_devs / realpots[ords])
        # dev_plot = plt.ScatterPlot(realpots[ords], unscaled_devs)
        dev_plot = plt.ScatterPlot(realpots[ords], unscaled_devs)
        # dev_plot = plt.ScatterPlot(realpots[ords], devs
        #                            # plot_range=[None, [-100, 100]]
        #                            )
        dev_plot.show()