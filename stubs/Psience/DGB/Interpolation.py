import abc
import numpy as np, scipy as sp, math
from McUtils.Zachary import InverseDistanceWeightedInterpolator, TensorDerivativeConverter
import McUtils.Numputils as nput
__all__ = ['DGBInterpolator', 'DGBGenericInterpolator', 'DGBWatsonInterpolator', 'DGBCartesianWatsonInterpolator']

class DGBInterpolator(metaclass=abc.ABCMeta):
    default_neighborhood_size = 10
    default_declustering_alpha = 0
    default_declustering_overlap = 0.1

    def __init__(self, centers, potential_derivatives, declustering_alpha=None, declustering_overlap=None, neighborhood_size=None, logger=None, pairwise_potential_functions=None, **opts):
        ...

    @abc.abstractmethod
    def __call__(self, points, deriv_order=None, **kwargs) -> 'np.ndarray | list[np.ndarray]':
        ...

class DGBGenericInterpolator(DGBInterpolator):

    def __init__(self, centers, potential_derivatives, **opts):
        ...

    def __call__(self, points, deriv_order=None, **kwargs):
        ...

class DGBWatsonInterpolatorSingleRef(DGBInterpolator):

    def __init__(self, centers, potential_derivatives, modes, **opts):
        ...

    def take_remainder_potential(self, centers, potential_derivatives, modes):
        ...

    def __call__(self, points, deriv_order=None, **kwargs):
        ...

class WatsonPairwisePotential:
    """
    Partial mimic of the PairwisePotentialEvaluator, but just taking
    the bits needed to actually convert f(q) -> f(delta) -> f(r) and get the
    derivatives cleanly
    """

    def __init__(self, ppfs: 'dict[(int, int), callable]', modes):
        ...

    @classmethod
    def get_bond_length_deltas(cls, natoms, ndim, i, j, full=False):
        ...

    def get_coordinate_bond_length_projection(self, i, j, ndim=3):
        ...

    def get_bond_lengths(self, coords, i, j, deriv_order=0):
        ...

    def eval_ppf(self, f, coords, i, j, deriv_order=None):
        ...

    def __call__(self, normal_mode_coords, deriv_order=None):
        ...

class DGBWatsonTaylorInterpolator(DGBInterpolator):

    def __init__(self, centers, potential_derivatives, modes, power=4, include_harmonic_basis=False, harmonic_distance_cutoff=None, pairwise_potential_functions=None, **opts):
        ...

    def take_remainder_potential(self, centers, potential_derivatives, modes):
        ...

    def take_ppf_remainder(self, centers, potential_derivatives):
        ...

    def taylor_interp(self, points, dists, neighbors, derivs, power=None, deriv_order=None):
        ...

    def __call__(self, points, deriv_order=None, *, neighborhood_size=None, power=None, **kwargs):
        ...

class DGBWatsonLeastSquaresInterpolator(DGBInterpolator):
    ...
DGBWatsonInterpolator = DGBWatsonTaylorInterpolator

class DGBCartesianWatsonInterpolator(DGBWatsonInterpolator):

    def __init__(self, centers, potential_derivatives, modes, **opts):
        ...

    def __call__(self, points, deriv_order=None, **kwargs):
        ...