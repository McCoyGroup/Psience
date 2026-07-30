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
__all__ = ['DGBRunner']
__reload_hook__ = ['.DGB', '.Coordinates', '.Wavefunctions']

class DGBRunner:

    @classmethod
    def prep_interpolation(cls, nms, coords, potential_function, symmetrizations=None):
        ...

    @classmethod
    def construct_from_mol_simulation(cls, sim, mol, *, potential_function=None, dipole_function=None, use_dipole_embedding=True, use_cartesians=False, use_momenta=False, quadrature_degree=3, expansion_degree=2, use_interpolation=True, use_quadrature=False, symmetrizations=None, momentum_scaling=None, skip_initial_configurations=True, alphas='virial', modes='normal', transformations='diag', pairwise_potential_functions=None, logger=True, **opts):
        ...

    @classmethod
    def construct_from_model_simulation(cls, sim, model, mol=None, *, use_cartesians=False, use_momenta=False, quadrature_degree=3, expansion_degree=2, use_interpolation=True, use_quadrature=False, symmetrizations=None, momentum_scaling=None, skip_initial_configurations=True, modes='normal', transformations='diag', **opts):
        ...

    @classmethod
    def construct_from_model(cls, model, trajectories=10, *, sim=None, propagation_time=10, timestep=50, use_cartesians=False, use_momenta=False, pairwise_potential_functions=None, use_interpolation=True, use_quadrature=False, symmetrizations=None, momentum_scaling=None, total_energy=None, total_energy_scaling=None, sampled_modes=None, initial_energies=None, initial_displacements=None, initial_mode_directions=None, displaced_coords=None, track_velocities=True, logger=None, **aimd_options):
        ...

    @classmethod
    def from_mol(cls, mol, sim=None, *, potential_function=None, dipole_function=None, trajectories=10, propagation_time=10, timestep=50, use_cartesians=False, use_momenta=False, pairwise_potential_functions=None, use_interpolation=True, use_quadrature=False, symmetrizations=None, momentum_scaling=None, trajectory_seed=None, total_energy=None, total_energy_scaling=None, sampled_modes=None, initial_energies=None, initial_displacements=None, initial_mode_directions=None, displaced_coords=None, track_velocities=True, logger=None, **aimd_options):
        ...

    @classmethod
    def run_simple(cls, system_spec, sim=None, plot_wavefunctions=True, plot_spectrum=True, trajectories=10, propagation_time=10, timestep=50, use_cartesians=False, use_momenta=False, pairwise_potential_functions=None, use_interpolation=True, use_quadrature=False, symmetrizations=None, momentum_scaling=None, trajectory_seed=None, total_energy=None, total_energy_scaling=None, sampled_modes=None, initial_energies=None, initial_mode_directions=None, initial_displacements=None, displaced_coords=None, **opts):
        ...
    plot_potential_cutoff = 17000
    plot_potential_units = 'Wavenumbers'

    @classmethod
    def plot_dgb_potential(cls, dgb, mol, potential, coordinate_sel=None, domain=None, domain_padding=1, potential_cutoff=None, potential_units=None, potential_min=0, plot_cartesians=None, plot_atoms=True, cmap=None, modes_nearest=False, plot_points=100, levels=24, **plot_styles):
        ...
    gaussian_plot_name = 'gaussian_{i}.pdf'

    @classmethod
    def plot_gaussians(cls, dgb, mol, *, domain=None, domain_padding=1, cmap='RdBu', plot_dir=None, plot_name=None, **plot_options):
        ...
    default_num_plot_wfns = 5
    wavefunction_plot_name = 'wfn_{i}.pdf'
    potential_plot_name = 'pot.png'

    @classmethod
    def plot_wavefunctions(cls, wfns, dgb, mol, which=True, coordinate_sel=None, cartesians=None, plot_dir=None, plot_name=None, plot_label='{e:.2f} cm-1', plot_potential=True, separate_potential=False, potential_plot_name=None, potential_units='Wavenumbers', plot_atoms=None, plot_centers=True, potential_styles=None, scaling=None, ticks=None, padding=None, aspect_ratio=None, plot_range=None, image_size=None, **plot_options):
        ...

    @classmethod
    def plot_potential_from_spec(cls, dgb, mol, spec, plot_centers=True, **opts):
        ...

    @classmethod
    def prep_plot_wavefunctions_spec(cls, dgb, spec):
        ...
    similarity_plot_name = 'similarity.png'
    spectrum_plot_name = 'spec.png'

    @classmethod
    def run_dgb(cls, dgb: DGB, mol, plot_centers=True, plot_wavefunctions=True, plot_spectrum=False, spectrum_plot_name=None, pot_cmap='viridis', wfn_cmap='RdBu', wfn_points=100, wfn_contours=12, plot_dir=None, plot_potential=True, pot_points=100, domain=None, domain_padding=1, wavefunction_scaling=None, potential_cutoff=None, potential_units='Wavenumbers', mode=None, nodeless_ground_state=None, min_singular_value=None, subspace_size=None, plot_similarity=False, similarity_plot_name=None, similarity_cutoff=None, similarity_chunk_size=None, similar_det_cutoff=None, num_print=None, spectrum_plotting_options=None, **plot_options):
        ...

    @classmethod
    def getMorseParameters(cls, w=None, wx=None, m1=None, m2=None, re=None):
        ...

    @classmethod
    def setupMorseFunction(cls, model, i, j, w=None, wx=None):
        ...

    @classmethod
    def plot_interpolation_error(cls, dgb, pot):
        ...