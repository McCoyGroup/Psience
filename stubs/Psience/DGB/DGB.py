import numpy as np, scipy as sp
from McUtils.Data import UnitsData
from McUtils.Zachary import DensePolynomial
from McUtils.Scaffolding import Logger, NullLogger
from McUtils.Parallelizers import Parallelizer
from .Gaussians import *
from .Evaluators import *
from .Coordinates import *
from .Interpolation import *
from .Wavefunctions import DGBWavefunctions
from .Solvers import DGBEigensolver
__all__ = ['DGB']
__reload_hook__ = ['.Wavefunctions', '..Molecools', '..MixtureModes']

class DGB:
    """

    """

    @classmethod
    def run_simple(cls, centers, potential_function, *, masses=None, atoms=None, modes=None, transformations=None, alphas=None, logger=True, clustering_radius=None, min_singular_value=None, num_svd_vectors=None, svd_contrib_cutoff=None, optimize_centers=None, quadrature_degree=None, expansion_degree=None, expansion_type=None, pairwise_potential_functions=None, dipole_function=None):
        ...
    printed_frequencies_number = 20

    def run(self, quiet=False, calculate_spectrum=True, dipole_degree=0, num_print=None, **wavefunction_options):
        """
        The default case...

        :return:
        """
        ...
    __props__ = ('potential_function', 'gmat_function', 'masses', 'atoms', 'alphas', 'transformations', 'internals', 'modes', 'coordinate_selection', 'cartesians', 'logger', 'parallelizer', 'optimize_centers', 'quadrature_degree', 'expansion_degree', 'expansion_type', 'momenta', 'poly_coeffs', 'pairwise_potential_functions', 'dipole_function', 'kinetic_options')

    @classmethod
    def construct(cls, centers, potential_function, gmat_function=None, masses=None, atoms=None, alphas=None, transformations=None, internals=None, modes=None, coordinate_selection=None, cartesians=None, logger=False, parallelizer=None, optimize_centers=False, quadrature_degree=3, expansion_degree=None, expansion_type='multicenter', momenta=None, poly_coeffs=None, pairwise_potential_functions=None, dipole_function=None, kinetic_options=None):
        ...

    @classmethod
    def construct_gaussians(cls, centers, alphas, potential_spec, gmat_function=None, masses=None, atoms=None, internals=None, coordinate_selection=None, cartesians=None, modes=None, kinetic_options=None, transformations=None, momenta=None, pairwise_potential_functions=None, poly_coeffs=None, logger=None, parallelizer=None):
        ...

    @classmethod
    def construct_potential(cls, potential_function, coords, quadrature_degree=None, expansion_degree=None, expansion_type=None, pairwise_potential_functions=None, logger=None, parallelizer=None):
        ...

    def __init__(self, gaussians: DGBGaussians, potential: DGBPotentialEnergyEvaluator, logger=None, parallelizer=None, wavefunction_options=None):
        ...

    def as_cartesian_dgb(self):
        ...

    def get_S(self):
        ...

    @property
    def S(self):
        ...

    def get_T(self):
        ...

    @property
    def T(self):
        ...

    @T.setter
    def T(self, T):
        ...

    def get_V(self):
        ...

    @property
    def V(self):
        ...

    @V.setter
    def V(self, V):
        ...

    def evaluate_multiplicative_operator(self, function, embed=True, expansion_degree=None, expansion_type=None, quadrature_degree=None, pairwise_functions=None):
        ...
    default_solver_mode = 'similarity'

    def diagonalize(self, *, mode=None, similarity_cutoff=None, similarity_chunk_size=None, similar_det_cutoff=None, similarity_shift=None, subspace_size=None, min_singular_value=None, nodeless_ground_state=True, low_rank_energy_cutoff=None, low_rank_overlap_cutoff=None, low_rank_shift=None, stable_eigenvalue_epsilon=None):
        ...

    def get_similarity_matrix(self):
        ...

    def get_wavefunctions(self, mode=None, similarity_cutoff=None, similarity_chunk_size=None, similar_det_cutoff=None, subspace_size=None, min_singular_value=None, nodeless_ground_state=None, stable_eigenvalue_epsilon=None, **wfn_opts):
        ...