import numpy as np, math
from McUtils.Coordinerds import CartesianCoordinates1D, CartesianCoordinates2D, CartesianCoordinates3D
from McUtils.Scaffolding import Logger
from McUtils.Parallelizers import Parallelizer
from McUtils.Data import AtomData, UnitsData
import McUtils.Numputils as nput
from ..Modes import NormalModes
from .Coordinates import *
from .Evaluators import *
from .Solvers import *
from .Interpolation import *
__all__ = ['DGBGaussians']

class DGBCovarianceTransformations:

    def __init__(self, tfs, invs=None, inverse_sum=None):
        ...

    def transform(self, derivative_operators):
        ...

    def inverse(self):
        ...

    def combine(self, other: 'DGBCovarianceTransformations'):
        ...

class DGBGaussians:
    """
    A class to set up the actual N-dimensional Gaussians used in a DGB
    """

    def __init__(self, coords, alphas, transformations=None, *, momenta=None, poly_coeffs=None, kinetic_options=None, logger=None, parallelizer=None):
        ...

    @property
    def overlap_data(self):
        ...

    def get_S(self, return_prefactor=False):
        ...

    def get_T(self):
        ...

    def optimize(self, optimizer_options, potential_function=None, logger=None, **opts):
        ...

    def take_gaussian_selection(self, full_good_pos):
        ...

    @staticmethod
    def _calc_gs_norm(new, S, old):
        ...

    @staticmethod
    def _calc_gs_next(vec, ortho):
        ...
    _fast_outer = None

    @classmethod
    def _gs_off_diag_norms(cls, S, rows, cols, vec):
        ...

    @staticmethod
    def _gs_diag_norms(vec):
        ...

    @staticmethod
    def _gs_extra_norm(old, vec):
        ...

    @classmethod
    def _gs_direct_norm(cls, S, vec):
        ...

    @classmethod
    def _calc_gs_norm_block(cls, new, S, old):
        ...

    @staticmethod
    def _calc_gs_next_block(vecs, ortho):
        ...

    @staticmethod
    def _calc_gs_next_block_single(cur, overlaps, prev):
        ...

    @classmethod
    def _optimize_gs_iter(cls, S, overlap_cutoff, norm_truncation_cutoff, logger=None):
        ...

    @staticmethod
    def _pivot_vector(vec, i, pivot):
        ...

    @classmethod
    def _pivot_matrix(cls, mat, i, pivot):
        ...

    @classmethod
    def _optimize_gs_block(cls, S, overlap_cutoff, max_overlap_component, logger=None):
        ...
    gs_optimization_overlap_cutoff = 0.001

    def _optimize_gs(self, *, S=None, norm_cutoff=None, norm_truncation_cutoff=0, max_overlap_cutoff=1, allow_pivoting=True, chunk_size=None, logger=None):
        ...

    def _optimize_svd(self, *, S=None, min_value=1e-12, num_vectors=None, contrib_cutoff=0.001, logger=None):
        ...
    default_energy_cutoff = 1600 / UnitsData.hartrees_to_wavenumbers

    def _prune_energy(self, pot, *, cutoff=None, probabilities=None, logger=None, potential_values=None):
        ...

    def _prune_dists(self, *, cluster_radius, logger=None):
        ...

    @classmethod
    def construct(cls, coords, alphas, *, potential_expansion=None, potential_function=None, transformations=None, masses=None, atoms=None, modes=None, kinetic_options=None, internals=None, coordinate_selection=None, cartesians=None, gmat_function=None, momenta=None, poly_coeffs=None, logger=None, pairwise_potential_functions=None, parallelizer=None):
        ...

    @classmethod
    def get_normal_modes(cls, coords, potential_function, masses=None, atoms=None, internals=None, gmat_function=None, reference_structure=None, stationary_point_norm=0.01, project_transrot=True):
        ...

    @classmethod
    def get_reaction_path_transformations(cls, coords, potential_function, gmat_function, stationary_point_norm=0.0001, sort_alphas=True):
        ...

    @classmethod
    def get_hessian_diagonalizing_transformations(cls, coords, potential_function, gmat_function, *, masses=None, project_transrot=True):
        ...

    @staticmethod
    def _filter_alpha_method_keys(method, opts, necessary_keys, optional_keys):
        ...

    @classmethod
    def dispatch_get_alphas(self, alphas, centers, **extra_opts):
        ...

    @classmethod
    def get_mass_alphas(cls, centers, *, masses, scaling=10, use_mean=False):
        ...

    @classmethod
    def get_min_distance_alphas(cls, masses, centers, scaling=1 / 4, use_mean=False):
        ...

    @classmethod
    def get_virial_alphas(cls, coords, *, potential_function, gmat_function, transformations, scaling=1 / 2):
        ...

    @staticmethod
    def _get_hermite_poly(coeff_dict, alphas):
        ...

    @classmethod
    def canonicalize_poly_coeffs(cls, coeffs, alphas):
        ...

    @property
    def transformations(self):
        ...

    @transformations.setter
    def transformations(self, tf):
        ...

    @classmethod
    def canonicalize_transforms(self, coords, tfs):
        ...

    @property
    def prefactor(self):
        ...

    @property
    def S(self):
        ...

    @S.setter
    def S(self, smat):
        ...

    @property
    def T(self):
        ...

    @T.setter
    def T(self, tmat):
        ...
    bad_alpha_limit = 1e-15
    bad_scaling_limit = 0.001

    def marginalize_out(self, indices, *, bad_alpha_limit=None, bad_scaling_limit=None):
        ...

    def as_cartesians(self):
        ...

    def plot_centers(self, figure=None, xyz_sel=None, **plot_styles):
        ...