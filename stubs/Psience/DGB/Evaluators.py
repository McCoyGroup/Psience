import abc, math, gc, time
import numpy as np, itertools, functools
from McUtils.Zachary import DensePolynomial, TensorDerivativeConverter
from McUtils.Scaffolding import Logger
from McUtils.Parallelizers import Parallelizer
import McUtils.Numputils as nput
__all__ = ['DGBEvaluator', 'DGBKineticEnergyEvaluator', 'DGBCartesianEvaluator', 'DGBWatsonEvaluator', 'DGBPotentialEnergyEvaluator', 'DGBPairwisePotentialEvaluator', 'DGBCartesianPairwiseEvaluator', 'DGBWatsonPairwiseEvaluator']

class OverlapGaussianData:

    def __init__(self, input_data, product_data):
        ...

    @property
    def npts(self):
        ...

    @property
    def ndim(self):
        ...

    @property
    def centers(self):
        ...

    @property
    def alphas(self):
        ...

    @property
    def rotations(self):
        ...

    @property
    def inverse_rotations(self):
        ...

    @property
    def precision_matrices(self):
        ...

    @property
    def covariance_matrices(self):
        ...

    @property
    def shift_matrices(self):
        ...

    @property
    def phases_diff(self):
        ...

    @property
    def phases_sum(self):
        ...

    @property
    def momenta_sum(self):
        ...

    @property
    def momenta_diff(self):
        ...

    @property
    def rho_sum(self):
        ...

    @property
    def rho_diff(self):
        ...

    @property
    def correlation_factors_sum(self):
        ...

    @property
    def correlation_factors_diff(self):
        ...

    @property
    def cos_correlation_diff(self):
        ...

    @property
    def sin_correlation_diff(self):
        ...

    @property
    def cos_correlation_sum(self):
        ...

    @property
    def sin_correlation_sum(self):
        ...

    @property
    def decay_factor_diff(self):
        ...

    @property
    def decay_factor_sum(self):
        ...

    @property
    def center_difference(self):
        ...

    @property
    def delta_position(self):
        ...

    @property
    def delta_phase_sum(self):
        ...

    @property
    def delta_phase_diff(self):
        ...

    @property
    def indices(self):
        ...

    @property
    def row_indices(self):
        ...

    @property
    def col_indices(self):
        ...

    @property
    def initial_centers(self):
        ...

    @property
    def initial_alphas(self):
        ...

    @property
    def initial_precision_matrices(self):
        ...

    @property
    def initial_momenta(self):
        ...

    @property
    def initial_phases(self):
        ...

    @classmethod
    def from_gaussian_parameters(cls, centers, alphas, transformations, momenta, *, chunk_size=None, rows_cols=None, logger=None, parallelizer=None):
        ...

    @classmethod
    def get_overlap_data(clget_s, centers, alphas, transformations, aligned_momenta, *, chunk_size=None, rows_cols=None, logger=None, parallelizer=None):
        ...

    def take_subselection(self, positions):
        ...

class DGBEvaluator:
    """
    An object that supports evaluating matrix elements in a distributed Gaussian basis.
    Provides support for integrating a function via quadrature or as an expansion in a polynomial tensors
    """

    @classmethod
    def get_inverse_covariances(cls, alphas, transformations):
        """
        Transforms the alphas into proper inverse covariance matrices.
        Chosen so that in the case that the transformations, Q, diagonalize S we can write
            QT S Q = A

        :return:
        :rtype:
        """
        ...

    @classmethod
    def get_covariances(cls, alphas, transformations):
        """
        Transforms the alphas into proper inverse covariance matrices.
        Chosen so that in the case that the transformations, Q, diagonalize S we can write
            QT S Q = A

        :return:
        :rtype:
        """
        ...

    @classmethod
    def get_momentum_vectors(cls, phases, transformations):
        """
        Transforms the momenta so that they're aligned along the Gaussian axes

        :return:
        :rtype:
        """
        ...

    @classmethod
    def get_phase_vectors(cls, momenta, transformations):
        """
        Transforms the momenta so that they're aligned along the Gaussian axes

        :return:
        :rtype:
        """
        ...

    @classmethod
    def get_overlap_gaussians(cls, centers, alphas, transformations, momenta, *, chunk_size=None, rows_cols=None, logger=None, parallelizer=None) -> 'OverlapGaussianData':
        ...

    @classmethod
    def poch(cls, n, m):
        ...

    @classmethod
    def polyint_1D(cls, centers, alphas, n):
        ...

    @classmethod
    def momentum_coeffient(cls, k, n):
        ...

    @classmethod
    def momentum_integral(cls, p, a, k):
        ...

    @classmethod
    def simple_poly_int(cls, n):
        ...

    @classmethod
    def tensor_expansion_integrate(cls, npts, derivs, overlap_data: 'OverlapGaussianData', expansion_type='multicenter', logger=None, reweight=True):
        """
        provides an integral from a polynomial expansion with derivs as an expansion in tensors

        :param npts:
        :param derivs:
        :param centers:
        :param alphas:
        :param inds:
        :param rot_data:
        :param expansion_type:
        :param logger:
        :return:
        """
        ...

    @classmethod
    def quad_weight_eval(cls, function, d_chunk, w_chunk, ndim, centers, squa):
        ...

    @classmethod
    def quad_nd(cls, centers, alphas, function, flatten=False, degree=3, chunk_size=int(1000000.0), normalize=True):
        """
        N-dimensional quadrature

        :param centers:
        :param alphas:
        :param function:
        :param degree:
        :return:
        """
        ...

    @classmethod
    def _wrap_rotated_function(cls, func, rotations, inverse, momenta):
        ...

    @classmethod
    def rotated_gaussian_quadrature(cls, function, alphas, centers, rotations, inverse, momenta, normalize=True, degree=2):
        ...

    @classmethod
    def quad_integrate(cls, function, overlap_data: 'OverlapGaussianData', degree=2, logger=None):
        """
        Integrate potential over all pairs of Gaussians at once

        :param degree:
        :type degree:
        :return:
        :rtype:
        """
        ...

    @classmethod
    def _rot_base_S_components(cls, overlap_data: 'OverlapGaussianData'):
        ...

    @classmethod
    def _rot_base_S(cls, overlap_data: 'OverlapGaussianData', logger=None, return_prefactor=False):
        ...

    @classmethod
    def evaluate_overlap(cls, overlap_data: 'OverlapGaussianData', logger=None, return_prefactor=False):
        ...

class DGBKineticEnergyEvaluator(DGBEvaluator):
    """
    """

    @classmethod
    def _evaluate_polynomial_ke(cls, overlap_data: 'OverlapGaussianData', terms, prefactors):
        ...

    @classmethod
    def _induced_polys(cls, terms, sigs, initial_centers, final_centers):
        ...

    @abc.abstractmethod
    def evaluate_ke(self, overlap_data: 'OverlapGaussianData', logger=None, **opts):
        ...

    @classmethod
    def evaluate_diagonal_rotated_momentum_contrib(self, overlap_data: 'OverlapGaussianData', masses):
        ...

class DGBCartesianEvaluator(DGBKineticEnergyEvaluator):

    def __init__(self, masses):
        ...

    def evaluate_ke(self, overlap_data: 'OverlapGaussianData', logger=None):
        ...

class DGBWatsonEvaluator(DGBKineticEnergyEvaluator):
    """
    """

    def __init__(self, modes, coriolis_inertia_function):
        ...

    @classmethod
    def _embed_watson_modes(cls, watson_data, centers):
        ...

    @staticmethod
    def annoying_coriolis_term(n, u, m, v, Xc, Dx, Sc, Sp, Gi, Gj, DG):
        ...

    @staticmethod
    def annoying_coriolis_momentum_term(n, u, m, v, Xc, r, Jp, Dx, Sc, Sp, DG):
        ...

    @staticmethod
    def annoying_imaginary_momentum_term(n, u, m, v, Xc, r, Jp, Dx, Sc, Sp, DG):
        ...

    @classmethod
    def evaluate_coriolis_contrib(cls, coriolis_tensors, overlap_data: 'OverlapGaussianData'):
        ...

    @classmethod
    def evaluate_watson_term(cls, B_e, overlap_data: 'OverlapGaussianData'):
        ...

    def evaluate_ke(self, overlap_data: 'OverlapGaussianData', logger=None, include_diagonal_contribution=True, include_coriolis_coupling=True, include_watson_term=True):
        ...

class DGBPotentialEnergyEvaluator(DGBEvaluator):
    """
    An evaluator designed
    """

    def __init__(self, potential_function, integral_handler=None, expansion_degree=None, expansion_type=None, quadrature_degree=None, pairwise_functions=None, logger=None):
        ...

    def analytic_integrate(self):
        ...

    @classmethod
    def expansion_integrate(cls, function, overlap_data: 'OverlapGaussianData', expansion_type, expansion_degree=2, pairwise_functions=None, logger=None):
        ...

    @classmethod
    def evaluate_multiplicative(cls, function, overlap_data: 'OverlapGaussianData', integral_handler=None, expansion_degree=None, expansion_type=None, quadrature_degree=None, pairwise_functions=None, logger=None):
        ...

    def evaluate_pe(self, overlap_data: 'OverlapGaussianData', logger=None):
        ...

    def evaluate_op(self, operator, overlap_data: 'OverlapGaussianData', integral_handler=None, expansion_degree=None, expansion_type=None, quadrature_degree=None, pairwise_functions=None, logger=None):
        ...

class DGBPairwisePotentialEvaluator(DGBEvaluator, metaclass=abc.ABCMeta):

    def __init__(self, coords, pairwise_potential_functions, quadrature_degree=3, use_with_interpolation='ignored'):
        ...

    @classmethod
    def get_bond_length_deltas(cls, natoms, ndim, i, j, full=False):
        ...

    @abc.abstractmethod
    def get_coordinate_bond_length_projection(self, i, j) -> 'tuple[np.ndarray, np.ndarray]':
        ...

    def get_coordinate_change_transformation(self, coordinate_projection_data) -> np.ndarray:
        ...

    def get_bond_length_change_transformation(self, overlap_data: 'OverlapGaussianData', i, j) -> np.ndarray:
        ...

    def wrap_distance_function(self, i, j, overlap_data: 'OverlapGaussianData', transformations, pairwise_function):
        ...

    def evaluate_pairwise_contrib(self, overlap_data: 'OverlapGaussianData', quadrature_degree=None, expansion_degree=2):
        ...

class DGBCartesianPairwiseEvaluator(DGBPairwisePotentialEvaluator):

    def __init__(self, coords: 'DGBCartesians', pairwise_functions, **opts):
        ...

    def get_coordinate_bond_length_projection(self, i, j):
        ...

class DGBWatsonPairwiseEvaluator(DGBPairwisePotentialEvaluator):

    def __init__(self, coords: 'DGBWatsonModes', pairwise_functions, **opts):
        ...

    def get_coordinate_bond_length_projection(self, i, j):
        ...