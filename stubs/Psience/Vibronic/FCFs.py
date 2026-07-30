import numpy as np, scipy as sp, itertools as it, collections, math, time, enum
import scipy.linalg
import McUtils.Numputils as nput
from McUtils.Scaffolding import Logger
from McUtils.Formatters import TableFormatter
from McUtils.Zachary import DensePolynomial
from McUtils.Combinatorics import IntegerPartitioner, IntegerPartitioner2D, UniquePermutations
from ..BasisReps import StateMaker
from ..Modes import NormalModes
__all__ = ['FranckCondonModel']

class State(enum.Enum):
    GroundState = 'gs'
    ExcitedState = 'es'

    @classmethod
    def get_aliases(cls):
        ...

    @classmethod
    def resolve(cls, method_name):
        ...

class RotationMethods(enum.Enum):
    LeastSquares = 'least-squares'
    Duschinsky = 'duschinsky'

    @classmethod
    def get_aliases(cls):
        ...

    @classmethod
    def resolve(cls, method_name):
        ...

class FranckCondonModel:
    default_rotation_method = 'default'
    default_rotation_order = State.GroundState
    default_rotation_center = State.ExcitedState
    default_embedding_ref = State.ExcitedState
    default_include_rotation = True

    def __init__(self, gs_nms: NormalModes, es_nms, atoms=None, *, logger=None, embed=True, embedding_ref=None, masses=None, mass_weight=True, dimensionless=False, mode_selection=None, mode_reordering=None, rotation_method=None, rotation_order=None, rotation_center=None, include_rotation=None, rotation_blocks=None):
        ...

    @classmethod
    def prep_modes(cls, gs_nms, es_nms, embed=True, embedding_ref=None, masses=None, mass_weight=False, dimensionless=False, mode_selection=None, mode_reordering=None, **rotation_opts):
        ...

    @classmethod
    def from_files(cls, gs_file, es_file, logger=None, mode_selection=None, mode_reordering=None, internals=None, internals_ref='gs', **rotation_embedding_opts):
        ...

    @classmethod
    def convert_internal_modes(cls, mol, nms):
        ...

    @classmethod
    def from_mols(cls, gs, es, logger=None, remove_transrot=True, use_internals=True, embed=True, mass_weight=True, **rotation_embedding_opts):
        ...

    def get_overlaps(self, excitations, *, duschinsky_cutoff=None, ground_states=None, return_states=True, **rotation_opts):
        ...
    OverlapData = collections.namedtuple('OverlapData', ['alphas', 'scaling', 'center', 'gs', 'es'])
    Embedding = collections.namedtuple('Embedding', ['modes', 'center'])

    def get_overlap_data(self, **rotation_embedding_opts) -> OverlapData:
        ...
    _rot_opts = {'rotation_center', 'rotation_order', 'include_rotation', 'rotation_method'}

    @classmethod
    def prep_overlap_args(self, gs_nms, es_nms):
        ...

    @classmethod
    def get_poly_evaluation_plan(self, exponents, alphas=None, zpe_prod=None):
        """
        Provides a function that can take a set of indices and rotation matrices
        from the gs and es bases to the shared basis of the central Gaussian and compute
        the corresponding term contributions by considering every even
        permutation of indices that could lead to a non-zero contribution

        :param tg:
        :param te:
        :return:
        """
        ...

    @classmethod
    def zero_point_alpha_contrib(cls, alphas):
        ...

    @classmethod
    def term_evaluator(self, exponents_list, splits_list, splits_inds_list, weights_list, gammas, alphas=None, zpe_prod=None):
        ...

    @staticmethod
    def _expand(pc, splits):
        ...

    @staticmethod
    def _contract(poly_exps, weights):
        ...

    @classmethod
    def _index_contract(cls, pc, weights, split_inds):
        ...

    @classmethod
    def evaluate_poly_chunks(cls, poly_coeffs, exps, splits, split_inds, weights, alphas, include_baseline=False):
        ...

    @classmethod
    def evaluate_poly_contrib_chunk(self, inds, exponents_list, splits_list, splits_inds_list, weights_list, alphas, coeffs):
        ...
    duschinsky_cutoff = 1e-20
    evaluator_plans = {}

    @classmethod
    def evaluate_shifted_poly_overlap(self, poly: 'HermiteProductPolynomial', Q, alphas, zpe_prod, duschinsky_cutoff=None):
        ...

    @classmethod
    def df_weights(cls, n):
        ...

    @classmethod
    def get_overlap_gaussian_data(cls, freqs_gs, modes_gs, inv_gs, center_gs, freqs_es, modes_es, inv_es, center_es, rotation_method=None, rotation_order=None, rotation_center=None, include_rotation=None, rotation_blocks=None):
        ...
    integral_block_size = 1000

    @classmethod
    def eval_fcf_overlaps(self, excitations_gs, freqs_gs, modes_gs, inv_gs, center_gs, excitations_es, freqs_es, modes_es, inv_es, center_es, duschinsky_cutoff=None, logger=None, **rotation_opts):
        """
        Evaluates the Gaussian overlaps between two H.O. wave functions defined by
        a set of polynomial coefficients, broadening factors, and centers, assuming
        the modes and centers are in an Eckart fream
        """
        ...

    @classmethod
    def embed_modes(cls, gs_nms: 'NormalModes', es_nms, ref=None, masses=None):
        ...

    @classmethod
    def mass_weight_nms(cls, nms, masses=None):
        ...

    @classmethod
    def make_dimensionless(cls, nms, freqs=None, masses=None):
        ...

    @classmethod
    def prep_states_from_threshold_and_quanta(cls, nms, *, threshold=None, min_freq=None, max_state=None, min_quanta=None, max_quanta=None):
        ...

    @classmethod
    def prep_states_from_excitations(cls, nms, *, states, **opts):
        ...
    state_space_prep_registry = {}

    @classmethod
    def state_space_prep_dispatchers(cls):
        ...

    @classmethod
    def dispatch_state_space_prep(cls, spec, nms):
        ...

    @classmethod
    def _check_listable(cls, excitations):
        ...

    @classmethod
    def prep_state_space(cls, excitations, nms, check=True):
        """
        Dispatcher to get appropriate state spaces
        :param excitations:
        :param check:
        :return:
        """
        ...

    @classmethod
    def get_fcfs(cls, gs_nms: 'NormalModes', es_nms: 'NormalModes', excitations, ground_states=None, duschinsky_cutoff=None, logger=None, **rotation_embedding_opts):
        ...

    @classmethod
    def format_overlap_tables(cls, es, overlaps, include_headers=True):
        ...

    @classmethod
    def get_fcf_spectrum(self, gs_nms: 'NormalModes', es_nms: 'NormalModes', excitations, ground_states=None, logger=None, duschinsky_cutoff=None, return_states=False, **rotation_embedding_opts):
        ...

    def prep_opts(self, **opts):
        ...

    def get_spectrum(self, excitations, *, ground_states=None, return_states=False, duschinsky_cutoff=None, **rotation_embedding_opts):
        ...

    def get_ezFCF_input(self, excitations, atoms=None, ground_states=None, **rotation_embedding_opts):
        ...

class HermiteProductPolynomial:

    def __init__(self, poly_dict: dict, ndim):
        ...

    def shift(self, s):
        ...

    def concat(self, other: 'HermiteProductPolynomial'):
        ...

    def term_iter(self, filter=None):
        ...

    @classmethod
    def from_quanta(cls, quanta, alphas):
        ...
    _hermite_cache = {}

    @classmethod
    def _hermite_coeffs(cls, n):
        ...

    @classmethod
    def get_1D_hermite_poly(cls, n, a):
        ...