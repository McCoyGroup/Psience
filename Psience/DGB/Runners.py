
import numpy as np, os
from McUtils.Data import AtomData, UnitsData
import McUtils.Devutils as dev
import McUtils.Plots as plt
import McUtils.Numputils as nput
from McUtils.Scaffolding import Logger

from ..Molecools import Molecule
from ..AIMD import AIMDSimulator

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
    def prep_interpolation(cls, nms, coords, potential_function, symmetrizations=None):
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
        return {
            'centers': interp_coords,
            'values': potential_function(interp_coords, order=2)
        }

    @classmethod
    def construct_from_mol_simulation(cls,
                                      sim, mol,
                                      *,
                                      potential_function=None,
                                      dipole_function=None,
                                      use_dipole_embedding=True,
                                      use_cartesians=False,
                                      use_momenta=False,
                                      quadrature_degree=3,
                                      expansion_degree=2,
                                      use_interpolation=True,
                                      use_quadrature=False,
                                      symmetrizations=None,
                                      momentum_scaling=None,
                                      skip_initial_configurations=True,
                                      alphas='virial',
                                      modes='normal',
                                      transformations='diag',
                                      pairwise_potential_functions=None,
                                      logger=True,
                                      **opts
                                      ):

        coords, velocities = sim.extract_trajectory(flatten=True, embed=mol.coords)
        if skip_initial_configurations:
            nconf = len(sim.coords)
            coords = coords[nconf - 1:]
            velocities = velocities[nconf - 1:]
        if use_momenta:
            momenta = velocities * sim.masses[np.newaxis, :, np.newaxis]
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
        if dipole_function is not None and use_dipole_embedding:
            mol = mol.modify(dipole_evaluator=dipole_function)
        potential_function = mol.get_energy_function()
        nms = mol.get_normal_modes(project_transrot=False, use_internals=False).remove_mass_weighting()

        if use_interpolation:
            potential_data = cls.prep_interpolation(
                nms,
                coords,
                mol.calculate_energy,
                symmetrizations=symmetrizations
            )
        else:
            def pot(coords, deriv_order=None):
                exp = potential_function(coords, order=deriv_order)
                return exp
            potential_data = pot

        if use_quadrature and use_interpolation:
            raise ValueError("don't use interpolation with quadrature...")

        if use_cartesians:
            diffs = [np.ptp(coords[:, i]) for i in range(3)]
            axes = [i for i, d in enumerate(diffs) if d > 1e-8]
        else:
            axes = None

        if pairwise_potential_functions is not None:
            auto_keys = [
                k for k,v in pairwise_potential_functions.items()
                if dev.str_is(v, 'auto')
            ]
            if len(auto_keys) > 0:
                pots = mol.get_1d_potentials(auto_keys)
                pairwise_potential_functions = pairwise_potential_functions.copy()
                for k,p in zip(auto_keys, pots):
                    def ppf(c, deriv_order=None, _p=p):
                        ord = 0 if deriv_order is None else deriv_order
                        _, conv = _p(c.reshape(-1, 1), preconverted=True, order=ord)
                        conv = [d.reshape(c.shape + d.shape[1:]) for d in conv]
                        # print(deriv_order, c.shape, [d.shape for d in conv])
                        if deriv_order is None:
                            conv = conv[0]
                        return conv
                    pairwise_potential_functions[k] = ppf

        if use_dipole_embedding or dipole_function is None:
            dipole = mol.get_dipole_function()
        else:
            dipole = dipole_function
        def dipole_data(coords, deriv_order=None):
            exp = dipole(coords, order=deriv_order)
            return exp

        if dev.str_is(modes, 'normal'):
            modes = nms
            # modes = mol.get_normal_modes(use_internals=False, project_transrot=True).remove_mass_weighting()
        return DGB.construct(
            coords,
            **dict(
                dict(
                    potential_function=potential_data,
                    dipole_function=dipole_data,
                    modes=None if use_cartesians else modes,
                    cartesians=axes
                    , alphas=alphas
                    , masses=mol.atomic_masses
                    , pairwise_potential_functions=pairwise_potential_functions
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

        coords = sim.extract_trajectory(flatten=True, embed=mol.coords)
        if sim.track_velocities:
            coords, velocities = coords
        else:
            velocities = None
        if skip_initial_configurations:
            nconf = len(sim.coords)
            coords = coords[nconf - 1:]
            if velocities is not None:
                velocities = velocities[nconf - 1:]
        if use_momenta and velocities is not None:
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
                             pairwise_potential_functions=None,
                             use_interpolation=True,
                             use_quadrature=False,
                             symmetrizations=None,
                             momentum_scaling=None,
                             total_energy=None,
                             total_energy_scaling=None,
                             sampled_modes=None,
                             initial_energies=None,
                             initial_displacements=None,
                             initial_mode_directions=None,
                             displaced_coords=None,
                             track_velocities=True,
                             logger=None,
                             **aimd_options
                             ):
        from ..AnalyticModels import MolecularModel
        model: MolecularModel

        logger = Logger.lookup(logger, construct=True)
        if sim is None:
            with logger.block(tag="Running AIMD"):
                logger.log_print(
                    "Normal mode frequencies: {freqs}",
                    freqs=model.normal_modes()[0] * UnitsData.convert("Hartrees", "Wavenumbers")
                )
                aimd_opts = dict(
                    trajectories=trajectories,
                    timestep=timestep,
                    track_velocities=True,
                    total_energy=total_energy,
                    total_energy_scaling=total_energy_scaling,
                    sampled_modes=sampled_modes,
                    initial_energies=initial_energies,
                    initial_displacements=initial_displacements,
                    initial_mode_directions=initial_mode_directions,
                    displaced_coords=displaced_coords,
                    **aimd_options
                )
                for k,v in aimd_opts.items():
                    if v is not None:
                        logger.log_print("{k} = {v}", k=k, v=v)
                sim = model.setup_AIMD(
                    trajectories=trajectories,
                    timestep=timestep,
                    track_velocities=True,
                    total_energy=total_energy,
                    total_energy_scaling=total_energy_scaling,
                    sampled_modes=sampled_modes,
                    initial_energies=initial_energies,
                    initial_displacements=initial_displacements,
                    initial_mode_directions=initial_mode_directions,
                    displaced_coords=displaced_coords,
                    **aimd_options
                )
                sim.propagate(propagation_time)

        with logger.block(tag="Constructing DGB"):
            print_opts = dict(
                use_cartesians=use_cartesians,
                use_momenta=use_momenta,
                pairwise_potential_functions=pairwise_potential_functions,
                use_interpolation=use_interpolation,
                use_quadrature=use_quadrature,
                momentum_scaling=momentum_scaling,
                symmetrizations=symmetrizations,
            )
            for k, v in print_opts.items():
                if v is not None:
                    logger.log_print("{k} = {v}", k=k, v=v)
            dgb = cls.construct_from_model_simulation(
                sim, model,
                use_cartesians=use_cartesians,
                use_momenta=use_momenta,
                pairwise_potential_functions=pairwise_potential_functions,
                use_interpolation=use_interpolation,
                use_quadrature=use_quadrature,
                momentum_scaling=momentum_scaling,
                symmetrizations=symmetrizations,
                logger=logger
            )
        return dgb

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
                 pairwise_potential_functions=None,
                 use_interpolation=True,
                 use_quadrature=False,
                 symmetrizations=None,
                 momentum_scaling=None,
                 trajectory_seed=None,
                 total_energy=None,
                 total_energy_scaling=None,
                 sampled_modes=None,
                 initial_energies=None,
                 initial_displacements=None,
                 initial_mode_directions=None,
                 displaced_coords=None,
                 track_velocities=True,
                 logger=None,
                 **aimd_options
                 ):

        aimd_options, dgb_opts = dev.OptionsSet(aimd_options).split(AIMDSimulator)

        if potential_function is not None:
            mol = mol.modify(energy_evaluator=potential_function)

        logger = Logger.lookup(logger, construct=True)
        if sim is None:
            with logger.block(tag="Running AIMD"):
                logger.log_print(
                    "Normal mode frequencies: {freqs}",
                    freqs=mol.get_normal_modes().freqs * UnitsData.convert("Hartrees", "Wavenumbers")
                )
                aimd_opts = dict(
                    trajectories=trajectories,
                    timestep=timestep,
                    seed=trajectory_seed,
                    total_energy=total_energy,
                    total_energy_scaling=total_energy_scaling,
                    sampled_modes=sampled_modes,
                    initial_energies=initial_energies,
                    initial_displacements=initial_displacements,
                    initial_mode_directions=initial_mode_directions,
                    displaced_coords=displaced_coords,
                    **aimd_options
                )
                for k,v in aimd_opts.items():
                    if v is not None:
                        logger.log_print("{k} = {v}", k=k, v=v)
            sim = mol.setup_AIMD(
                trajectories=trajectories,
                timestep=timestep,
                track_velocities=True,
                seed=trajectory_seed,
                total_energy=total_energy,
                total_energy_scaling=total_energy_scaling,
                sampled_modes=sampled_modes,
                initial_energies=initial_energies,
                initial_displacements=initial_displacements,
                initial_mode_directions=initial_mode_directions,
                displaced_coords=displaced_coords,
                **aimd_options
            )
            sim.propagate(propagation_time)

        with logger.block(tag="Constructing DGB"):
            print_opts = dict(
                use_cartesians=use_cartesians,
                use_momenta=use_momenta,
                pairwise_potential_functions=pairwise_potential_functions,
                use_interpolation=use_interpolation,
                use_quadrature=use_quadrature,
                symmetrizations=symmetrizations,
                momentum_scaling=momentum_scaling,
                dipole_function=dipole_function,
                **dgb_opts
            )
            for k, v in print_opts.items():
                if v is not None:
                    logger.log_print("{k} = {v}", k=k, v=v)

            dgb = cls.construct_from_mol_simulation(
                sim, mol,
                potential_function=None,
                use_cartesians=use_cartesians,
                use_momenta=use_momenta,
                pairwise_potential_functions=pairwise_potential_functions,
                use_interpolation=use_interpolation,
                use_quadrature=use_quadrature,
                symmetrizations=symmetrizations,
                momentum_scaling=momentum_scaling,
                dipole_function=dipole_function,
                logger=logger,
                **dgb_opts
            )
        return dgb

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

    # @classmethod
    # def run_dgb(cls,
    #             dgb,
    #             *,
    #             plot_wavefunctions=True,
    #             plot_spectrum=True
    #             ):
    #
    #         logger = dgb.logger
    #         with logger.block(tag="Running DGB"):
    #             logger.log_print("num coords: {c}", c=len(dgb.gaussians.coords.centers))
    #             with logger.block(tag="S"):
    #                 logger.log_print(logger.prep_array(dgb.S[:5, :5]))
    #             with logger.block(tag="T"):
    #                 logger.log_print(logger.prep_array(dgb.T[:5, :5]))
    #             with logger.block(tag="V"):
    #                 logger.log_print(logger.prep_array(dgb.V[:5, :5]))
    #
    #             wfns_cart = dgb.get_wavefunctions()
    #             with logger.block(tag="Energies"):
    #                 logger.log_print(
    #                     logger.prep_array(wfns_cart.energies[:5] * UnitsData.convert("Hartrees", "Wavenumbers"))
    #                 )
    #             with logger.block(tag="Frequencies"):
    #                 logger.log_print(
    #                     logger.prep_array(wfns_cart.frequencies()[:5] * UnitsData.convert("Hartrees", "Wavenumbers"))
    #                 )
    #
    #             if plot_wavefunctions:
    #                 for i in range(4):
    #                     wfns_cart[i].plot().show()
    #
    #             if plot_spectrum:
    #                 spec = wfns_cart[:4].get_spectrum()
    #                 with logger.block(tag="Intensities"):
    #                     logger.log_print(
    #                         logger.prep_array(spec.intensities)
    #                     )
    #                 spec.plot().show()
    #             else:
    #                 spec = None
    #         return wfns_cart, spec
    @classmethod
    def run_simple(cls,
                   system_spec,
                   sim=None,
                   plot_wavefunctions=True,
                   plot_spectrum=True,
                   trajectories=10,
                   propagation_time=10,
                   timestep=50,
                   use_cartesians=False,
                   use_momenta=False,
                   pairwise_potential_functions=None,
                   use_interpolation=True,
                   use_quadrature=False,
                   symmetrizations=None,
                   momentum_scaling=None,
                   trajectory_seed=None,
                   total_energy=None,
                   total_energy_scaling=None,
                   sampled_modes=None,
                   initial_energies=None,
                   initial_mode_directions=None,
                   initial_displacements=None,
                   displaced_coords=None,
                   **opts
                   ):

        base_opts = dev.OptionsSet(opts)
        opts, run_opts = base_opts.split(None, props=AIMDSimulator.__props__ + DGB.__props__ + (
            'use_dipole_embedding',
        ))
        opts.update(
            trajectories=trajectories,
            propagation_time=propagation_time,
            timestep=timestep,
            use_cartesians=use_cartesians,
            use_momenta=use_momenta,
            pairwise_potential_functions=pairwise_potential_functions,
            use_interpolation=use_interpolation,
            use_quadrature=use_quadrature,
            symmetrizations=symmetrizations,
            momentum_scaling=momentum_scaling,
            trajectory_seed=trajectory_seed,
            total_energy=total_energy,
            total_energy_scaling=total_energy_scaling,
            sampled_modes=sampled_modes,
            initial_energies=initial_energies,
            initial_mode_directions=initial_mode_directions,
            initial_displacements=initial_displacements,
            displaced_coords=displaced_coords
        )

        if isinstance(system_spec, (str, dict, tuple)):
            system_spec = Molecule.construct(system_spec)
            # system_spec.get_dipole_evaluator()
        if hasattr(system_spec, 'potential'): # analyt
            dgb = cls.construct_from_model(system_spec, sim=sim, **opts)
            mol = system_spec.mol
        else:
            dgb = cls.from_mol(system_spec, sim=sim, **opts)
            mol = system_spec

        return dgb, cls.run_dgb(
            dgb,
            mol,
            plot_wavefunctions=plot_wavefunctions,
            plot_spectrum=plot_spectrum,
            **run_opts
        )

    plot_potential_cutoff = 17000
    plot_potential_units = 'Wavenumbers'
    @classmethod
    def plot_dgb_potential(cls,
                           dgb, mol, potential,
                           coordinate_sel=None,
                           domain=None, domain_padding=1,
                           potential_cutoff=None,
                           potential_units=None,
                           potential_min=0,
                           plot_cartesians=None,
                           plot_atoms=True,
                           cmap=None,
                           modes_nearest=False,
                           plot_points=100,
                           levels=24,
                           **plot_styles
                           ):

        if potential_cutoff is None:
            potential_cutoff = cls.plot_potential_cutoff
        if potential_units is None:
            potential_units = cls.plot_potential_units
        def cutoff_pot(points,
                       cutoff=potential_cutoff,
                       cutmin=potential_min
                       ):
            values = potential(points) * UnitsData.convert('Hartrees', potential_units)
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
                    cutoff_pot(points),
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

    gaussian_plot_name = 'gaussian_{i}.pdf'
    @classmethod
    def plot_gaussians(cls,
                       dgb, mol,
                       *,
                       domain=None,
                       domain_padding=1,
                       cmap='RdBu',
                       plot_dir=None,
                       plot_name=None,
                       **plot_options
                       ):
        n = len(dgb.gaussians.coords.centers)
        wfns = DGBWavefunctions(
            np.zeros(n),
            np.eye(n),
            dgb
        )
        if plot_name is None:
            plot_name = cls.gaussian_plot_name
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
    wavefunction_plot_name = 'wfn_{i}.pdf'
    potential_plot_name = 'pot.png'
    @classmethod
    def plot_wavefunctions(cls,
                           wfns, dgb, mol,
                           which=True,
                           coordinate_sel=None,
                           cartesians=None,
                           plot_dir=None,
                           plot_name=None,
                           plot_label='{e:.2f} cm-1',
                           plot_potential=True,
                           separate_potential=False,
                           potential_plot_name=None,
                           potential_units='Wavenumbers',
                           plot_atoms=None,
                           plot_centers=True,
                           potential_styles=None,
                           scaling=None,
                           ticks=None,
                           padding=None,
                           aspect_ratio=None,
                           plot_range=None,
                           image_size=None,
                           **plot_options
                           ):

        figure = None
        if cartesians:
            wfns = wfns.as_cartesian_wavefunction()
            dgb = wfns.hamiltonian

        if dgb.gaussians.coords.shape[-1] == 1:
            if scaling is None: scaling = .05
            scaling *= UnitsData.convert("Hartrees", potential_units)
        elif scaling is None:
            scaling = .2

        if coordinate_sel is None:
            coordinate_sel = list(range(dgb.gaussians.alphas.shape[-1]))

        if plot_name is None:
            plot_name = cls.wavefunction_plot_name

        layout_opts = dict(
            ticks=ticks,
            padding=padding,
            aspect_ratio=aspect_ratio,
            plot_range=plot_range,
            image_size=image_size
        )
        plot_options.update(layout_opts)
        if potential_styles is None:
            potential_styles = {}
        for k,v in layout_opts.items():
            if k not in potential_styles:
                potential_styles[k] = v
        figs = []
        pot_figure = None
        if which is True:
            which = cls.default_num_plot_wfns
        if isinstance(which, int):
            which = list(range(which))
        if isinstance(which, list) and len(which) > 0:
            pot = dgb.pot.potential_function
            if plot_potential:
                if plot_atoms is None:
                    plot_atoms = bool(plot_centers)
                pot_figure = cls.plot_dgb_potential(
                    dgb, mol, pot,
                    coordinate_sel=coordinate_sel,
                    plot_atoms=plot_atoms,
                    potential_units=potential_units,
                    **potential_styles
                )
                if separate_potential and plot_dir is not None:
                    if potential_plot_name is None:
                        potential_plot_name = cls.potential_plot_name
                    pot_figure.savefig(os.path.join(plot_dir, potential_plot_name))
                    pot_figure.close()

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
                    if separate_potential:
                        figure = None
                    elif pot_figure is not None:
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
                                    scaling=scaling,
                                    **plot_options
                                )
                            )
                        else:
                            figs.append(
                                w.plot(
                                    figure=figure,
                                    plot_centers=plot_centers,
                                    plot_label=lab,
                                    scaling=scaling,
                                    **plot_options
                                )
                            )

                    if plot_dir is not None:
                        fig = figs.pop()
                        fig.background = "#ffffff00"
                        fig.savefig(
                            os.path.join(plot_dir, plot_name.format(i=i)),
                            transparent=True,
                            # facecolor="#ffffff00"
                        )
                        fig.close()

            # if plot_dir is None:
            #     figs[0].show()

        if separate_potential:
            return pot_figure, figs
        else:
            return figs

    @classmethod
    def plot_potential_from_spec(cls, dgb, mol, spec,
                                 plot_centers=True,
                                 **opts
                                 ):
        spec = cls.prep_plot_wavefunctions_spec(dgb, spec)
        cartesian_plot_axes = None
        if isinstance(spec, dict):
            if 'cartesians' in spec:
                cartesian_plot_axes = spec['cartesians']
                dgb = dgb.as_cartesian_dgb()

        pot = dgb.pot.potential_function
        figure = cls.plot_dgb_potential(
            dgb, mol, pot,
            **opts
        )

        if plot_centers:
            dgb.gaussians.plot_centers(
                figure,
                xyz_sel=cartesian_plot_axes
            )

        return figure

    @classmethod
    def prep_plot_wavefunctions_spec(cls, dgb, spec):
        if dev.str_is(spec, 'cartesians'):
            spec = {'cartesians': None}
        elif spec is True:
            if dgb.gaussians.coords.centers.shape[-1] > 2:
                spec = {'cartesians': None}
            else:
                spec = {}
        elif dev.is_dict_like(spec) and spec.get('cartesians', False) is None:
            embed_coords = dgb.gaussians.coords.as_cartesians()[0].centers
            diffs = [np.ptp(embed_coords[:, i]) for i in range(3)]
            axes = [i for i, d in enumerate(diffs) if d > 1e-8]
            spec['cartesians'] = axes

        return spec

    similarity_plot_name = 'similarity.png'
    spectrum_plot_name = 'spec.png'
    @classmethod
    def run_dgb(cls,
                dgb: DGB,
                mol,
                plot_centers=True,
                plot_wavefunctions=True,
                plot_spectrum=False,
                spectrum_plot_name=None,
                pot_cmap='viridis',
                wfn_cmap='RdBu',
                wfn_points=100,
                wfn_contours=12,
                plot_dir=None,
                plot_potential=True,
                pot_points=100,
                domain=None,
                domain_padding=1,
                wavefunction_scaling=None,
                potential_cutoff=None,
                potential_units='Wavenumbers',
                mode=None,
                nodeless_ground_state=None,
                min_singular_value=None,
                subspace_size=None,
                plot_similarity=False,
                similarity_plot_name=None,
                similarity_cutoff=None,
                similarity_chunk_size=None,
                similar_det_cutoff=None,
                num_print=None,
                **plot_options
                ):

        logger = dgb.logger
        with logger.block(tag="Running DGB"):
            # print("--->", len(dgb.gaussians.coords.centers))
            # print(dgb.S[:5, :5])
            # print(dgb.T[:5, :5])
            # print(dgb.V[:5, :5])
            dgb_opts = dict(
                mode=mode,
                nodeless_ground_state=nodeless_ground_state,
                subspace_size=subspace_size,
                min_singular_value=min_singular_value,
                similarity_cutoff=similarity_cutoff,
                similarity_chunk_size=similarity_chunk_size,
                similar_det_cutoff=similar_det_cutoff
            )
            for k,v in dgb_opts.items():
                if v is not None:
                    logger.log_print("{k} = {v}", k=k, v=v)

            try:
                wfns, spec = dgb.run(
                    calculate_spectrum=plot_spectrum,
                    mode=mode,
                    nodeless_ground_state=nodeless_ground_state,
                    subspace_size=subspace_size,
                    min_singular_value=min_singular_value,
                    similarity_cutoff=similarity_cutoff,
                    similarity_chunk_size=similarity_chunk_size,
                    similar_det_cutoff=similar_det_cutoff,
                    num_print=num_print
                )
            except Exception as e:
                if plot_wavefunctions is not False:
                    figure = cls.plot_potential_from_spec(
                        dgb, mol, plot_wavefunctions,
                        domain=domain,
                        domain_padding=domain_padding,
                        potential_cutoff=potential_cutoff
                    )
                    figure.show()
                raise e
            else:
                if plot_similarity:
                    sim_plot = plt.ArrayPlot(dgb.get_similarity_matrix())
                    if plot_dir is not None:
                        os.makedirs(plot_dir, exist_ok=True)
                        if similarity_plot_name is None:
                            similarity_plot_name = cls.similarity_plot_name
                        sim_plot.savefig(os.path.join(plot_dir, similarity_plot_name))
                        sim_plot.close()
                    else:
                        sim_plot.show()

                if plot_spectrum:
                    # which = plot_wavefunctions
                    # if which is True:
                    #     which = cls.default_num_plot_wfns
                    # if isinstance(which, int):
                    #     which = list(range(which))
                    # print(which, spec.intensities, spec.frequencies, spec[which])
                    spec_plot = spec.plot()
                    if plot_dir is not None:
                        os.makedirs(plot_dir, exist_ok=True)
                        if spectrum_plot_name is None:
                            spectrum_plot_name = cls.spectrum_plot_name
                        spec_plot.savefig(os.path.join(plot_dir, spectrum_plot_name))
                        spec_plot.close()
                    else:
                        spec_plot.show()


                plot_wavefunctions = cls.prep_plot_wavefunctions_spec(dgb, plot_wavefunctions)
                use_cartesians = False
                if plot_wavefunctions or dev.is_dict_like(plot_wavefunctions):
                    coordinate_sel = None
                    if dev.is_dict_like(plot_wavefunctions):
                        if 'cartesians' in plot_wavefunctions:
                            use_cartesians = True
                            coordinate_sel = plot_wavefunctions['cartesians']
                            plot_wavefunctions = plot_wavefunctions.get('num', True)
                        elif 'modes' in plot_wavefunctions:
                            coordinate_sel = plot_wavefunctions['modes']
                            plot_wavefunctions = plot_wavefunctions.get('num', True)
                            plot_potential = len(coordinate_sel) < 2
                        else:
                            coordinate_sel = plot_wavefunctions.get('coords')
                            plot_wavefunctions = plot_wavefunctions.get('num', True)
                            # plot_potential = False
                        # else:
                        #     raise ValueError(plot_wavefunctions)

                    plot_wfns_opts = dict(
                        which=plot_wavefunctions,
                        cartesians=use_cartesians,
                        coordinate_sel=coordinate_sel,
                        potential_units=potential_units,
                        plot_dir=plot_dir,
                        contour_levels=wfn_contours,
                        cmap=wfn_cmap,
                        plot_points=wfn_points,
                        plot_centers={'color': 'red'} if plot_centers is True else plot_centers,
                        domain=domain,
                        domain_padding=domain_padding,
                        plot_potential=plot_potential,
                        scaling=wavefunction_scaling,
                        plotter=plt.TriContourLinesPlot,
                        potential_styles=dict(
                            domain=domain,
                            domain_padding=domain_padding,
                            cmap=pot_cmap,
                            plot_points=pot_points
                        ),
                        **plot_options
                    )
                    with dgb.logger.block(tag="Plotting wavefunctions"):
                        for k,v in plot_wfns_opts.items():
                            dgb.logger.log_print("{k} = {v}", k=k,v=v)

                        wfns_plots = cls.plot_wavefunctions(
                            wfns,
                            dgb,
                            mol,
                            **plot_wfns_opts
                        )
                else:
                    wfns_plots = None

                if plot_spectrum:
                    return (wfns, wfns_plots), spec
                else:
                    return (wfns, wfns_plots)


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