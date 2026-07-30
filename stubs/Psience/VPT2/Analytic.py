"""
Provides a symbolic approach to vibrational perturbation theory based on a Harmonic description
"""
import abc, itertools, collections, enum, math
import contextlib
import functools
import numpy as np, scipy.signal, time
from McUtils.Zachary import DensePolynomial, TensorCoefficientPoly
import McUtils.Numputils as nput
import McUtils.Devutils as dev
import McUtils.Combinatorics as mcomb
from McUtils.Combinatorics import SymmetricGroupGenerator, IntegerPartitioner, UniquePartitions, UniquePermutations
from McUtils.Scaffolding import Logger, Checkpointer, MaxSizeCache
from McUtils.Parallelizers import Parallelizer
from ..BasisReps import BasisStateSpace, HarmonicOscillatorProductBasis, HarmonicOscillatorMatrixGenerator, HarmonicOscillatorRaisingLoweringPolyTerms
from .Corrections import BasicAPTCorrections
__all__ = ['PerturbationTheoryEvaluator', 'AnalyticPerturbationTheorySolver']
_DEBUG_PRINT = False
_DEBUG_AUDIT = False
_PERMUTE_CHANGES = False
_PERMUTE_FINALS = False
_TAKE_UNIQUE_CHANGES = False

class DefaultValues(enum.Enum):
    DEFAULT = 'default'
default = DefaultValues.DEFAULT

class AnalyticPerturbationTheorySolver:
    """
    A re-attempt at using the recursive expressions
    to provide simpler code for getting APT expressions
    """

    def __init__(self, hamiltonian_expansion, logger=None, checkpoint=None, allowed_terms=None, allowed_coefficients=None, disallowed_coefficients=None, allowed_energy_changes=None, intermediate_normalization=None):
        ...

    @classmethod
    def from_order(cls, order, internals=True, logger=None, checkpoint=None, allowed_terms=None, allowed_coefficients=None, disallowed_coefficients=None, allowed_energy_changes=None, intermediate_normalization=None):
        ...

    def modify_hamiltonian(self, hamiltonian_corrections):
        ...
    _op_maps = {}

    def get_correction(self, key, cls, order, **kw):
        ...

    def shifted_hamiltonian_correction(self, order, **kw):
        ...

    def energy_correction(self, order, **kw):
        ...

    def wavefunction_correction(self, order, **kw):
        ...

    def overlap_correction(self, order, degenerate_changes=None, **kw):
        ...

    def full_wavefunction_correction(self, order, **kw):
        ...

    def operator_correction(self, order, operator_type=None, **kw):
        ...

    def operator_degenerate_correction(self, order, /, degenerate_changes, operator_type=None, **kw):
        ...

    def reexpressed_hamiltonian(self, order, **kw):
        ...

    def reexpressed_hamiltonian_degenerate_correction(self, order, /, degenerate_changes, **kw):
        ...
    operator_expansion_index = 5

    @classmethod
    def operator_expansion_terms(cls, order, logger=None, base_index=None, operator_type=None):
        ...

    @classmethod
    def clear_caches(cls):
        ...

class PolynomialInterface(metaclass=abc.ABCMeta):
    """
    Provides a basic interface to allow for the uniform manipulation
    of objects that dispatch down to some form of scalar multiplied by
    a sum of polynomials
    """

    def format_expr(self) -> str:
        ...

    @property
    @abc.abstractmethod
    def ndim(self) -> int:
        ...

    @abc.abstractmethod
    def audit(self, target, ignore_constants=True):
        ...

    @abc.abstractmethod
    def ensure_dimension(self, ndim) -> 'Self':
        ...

    def align_dimensions(self, other) -> 'tuple[Self, Self]':
        ...

    @abc.abstractmethod
    def shift(self, shift) -> 'Self':
        ...

    @abc.abstractmethod
    def scale(self, scaling) -> 'Self':
        ...

    @abc.abstractmethod
    def permute(self, perm, check_perm=True, allow_padding=False) -> 'Self':
        ...

    @abc.abstractmethod
    def mul_simple(self, other: 'PolynomialInterface') -> 'PolynomialInterface':
        ...

    @abc.abstractmethod
    def rmul_simple(self, other: 'PolynomialInterface') -> 'PolynomialInterface':
        ...

    @abc.abstractmethod
    def mul_along(self, other: 'PolynomialInterface', inds, remainder=None, mapping=None) -> 'PolynomialInterface':
        ...

    @abc.abstractmethod
    def rmul_along(self, other: 'PolynomialInterface', inds, remainder=None, mapping=None) -> 'PolynomialInterface':
        ...

    @abc.abstractmethod
    def combine(self, **kwargs) -> 'Self':
        ...

    @abc.abstractmethod
    def mutate(self, *args, **kwargs) -> 'Self':
        ...

class PolyPath:
    """
    A simple holder class that contains the set of modifications along each dimension
    used to build a given polynomial
    """

    def __init__(self, paths, scaling):
        ...

    def as_tuple(self):
        ...

    def __add__(self, other):
        ...

    def mul_along(self, other, inds, remainder):
        ...

class TreeSerializer:

    @classmethod
    def serialize_iterable(cls, iterable, primitive_test, concat, track_shapes=True):
        ...

    @classmethod
    def build_dict_trees(cls, dict_obj):
        ...

    @classmethod
    def serialize_tree_dict(cls, dict_obj, key_primitive_test=None, vals_primitive_test=None, concat=None):
        ...

    @staticmethod
    def default_concat(arrays):
        ...

    @staticmethod
    def default_prim(obj):
        ...

    @classmethod
    def deserialize_subiterable(cls, shape, flat, i, pad, depth):
        ...

    @classmethod
    def deserialize_iterable(cls, shape, flat, depth):
        ...

    @classmethod
    def _tuplate(cls, key):
        ...

    @classmethod
    def stitch_dict(cls, iterables):
        ...

    @classmethod
    def rebuild_tree_dict(cls, key_trees, tree_shape, val_buffers):
        ...

class ProductPTPolynomial(PolynomialInterface):
    """
    TODO: include prefactor term so we can divide out energy changes
    """

    def __init__(self, coeffs, prefactor=1, idx=None, steps=None):
        ...
    _cache = {}

    @classmethod
    def lookup(cls, idx):
        ...

    def mutate(self, coeffs: 'list[np.ndarray]'=default, prefactor: 'Any'=default, idx=None, steps=default):
        ...

    def to_state(self, serializer=None):
        ...

    @property
    def ndim(self):
        ...

    @property
    def order(self):
        ...

    def audit(self, target, ignore_constants=True):
        ...

    def ensure_dimension(self, ndim):
        ...

    def pad(self, left_right_pads):
        ...

    def __repr__(self):
        ...

    @classmethod
    def prep_float(cls, c):
        ...

    @classmethod
    def format_simple_poly(cls, coeffs, idx):
        ...

    def format_expr(self):
        ...

    def permute(self, new_inds, check_perm=True, allow_padding=False):
        ...

    def constant_rescale(self):
        """
        rescales so constant term is 1

        :return:
        """
        ...

    @property
    def monic_coeffs(self):
        ...

    @monic_coeffs.setter
    def monic_coeffs(self, coeffs):
        ...

    @staticmethod
    def _monify(coeffs, zero_thresh=1e-08):
        ...

    @staticmethod
    def _find_off_pos(self_coeffs, other_coeffs):
        ...

    def combine(self, other: 'ProductPTPolynomial'):
        ...

    def condense(self, inds=None, return_inds=False, check_inds=True):
        ...

    def shift(self, shift):
        ...

    def __mul__(self, other):
        ...

    def __rmul__(self, other):
        ...

    def scale(self, scalar):
        ...

    def mul_simple(self, other: 'PolynomialInterface'):
        ...

    def rmul_simple(self, other: 'PolynomialInterface') -> 'PolynomialInterface':
        ...
    _prod_poly_cache = {}

    @classmethod
    def _poly_mul(cls, self, other):
        ...

    @classmethod
    def fast_ind_remainder(cls, n, diff):
        ...

    @classmethod
    def get_index_mapping(cls, dim_1, dim_2, inds, return_remainder=True):
        """
        Computes the corresponding and indices remainders for multdimensional polynomial multiplications

        :param dim_1:
        :param dim_2:
        :param inds:
        :return:
        """
        ...

    def mul_along(self, other: 'PolynomialInterface', inds, remainder=None, mapping=None):
        ...

    def rmul_along(self, other, inds, remainder=None, mapping=None):
        """
        We multiply every subpoly along the given indices, transposing the appropriate tensor indices

        :param other:
        :param inds:
        :param remainder:
        :return:
        """
        ...

    def __add__(self, other):
        ...

    def __radd__(self, other):
        ...

class ProductPTPolynomialSum(PolynomialInterface):

    def __init__(self, polynomials, prefactor=1, ndim=None, order=None, reduced=False):
        ...

    def prep_serialization_dict(self):
        ...

    @classmethod
    def from_serialization_dict(cls, big_dict):
        ...

    def to_state(self, serializer=None):
        ...

    def mutate(self, polynomials: 'list[PolynomialInterface]'=default, prefactor: 'Any'=default, ndim: 'Any'=default, order: 'Any'=default, reduced: 'Any'=default):
        ...

    @property
    def ndim(self):
        ...

    @property
    def order(self):
        ...

    def audit(self, target=None, ignore_constants=True):
        ...

    def ensure_dimension(self, ndim):
        ...

    def __repr__(self):
        ...

    def format_expr(self):
        ...

    def permute(self, new_inds, check_perm=True, allow_padding=False):
        ...

    @classmethod
    def combine_polys(cls, poly_set, cache):
        ...

    def combine(self):
        ...

    def condense(self, inds=None, return_inds=False, check_inds=True):
        ...

    def shift(self, shift):
        ...

    def __mul__(self, other):
        ...

    def __rmul__(self, other):
        ...

    def scale(self, scalar):
        ...

    def mul_simple(self, other: 'ProductPTPolynomial'):
        ...

    def rmul_simple(self, other: 'PolynomialInterface') -> 'PolynomialInterface':
        ...

    def mul_along(self, other: 'PolynomialInterface', inds, remainder=None, mapping=None):
        ...

    def rmul_along(self, other: 'PolynomialInterface', inds, remainder=None, mapping=None):
        """
        We multiply every subpoly along the given indices, transposing the appropriate tensor indices

        :param other:
        :param inds:
        :param remainder:
        :return:
        """
        ...

    def __add__(self, other):
        ...

    def __radd__(self, other):
        ...

class PTEnergyChangeProductSum(TensorCoefficientPoly, PolynomialInterface):
    """
    A representation of a sum of 1/energy * poly sums
    which is here so we can transpose energy change indices intelligently
    """

    def __init__(self, terms: dict, prefactor=1, canonicalize=True, reduced=False):
        ...

    def mutate(self, terms: dict=default, prefactor: 'Any'=default, canonicalize: 'Any'=default, reduced: 'Any'=default):
        ...

    def flip_energy_terms(self):
        ...

    def prep_serialization_dict(self):
        ...

    @classmethod
    def from_serialization_dict(cls, big_dict):
        ...

    def to_state(self, serializer=None):
        ...

    @classmethod
    def from_state(cls, state, serializer=None):
        ...

    def filter(self, terms, mode='match'):
        ...

    @classmethod
    def canonical_key(cls, monomial_tuple):
        ...

    @classmethod
    def side_change_iter(cls, key):
        ...

    @classmethod
    def format_key(self, key):
        ...

    def __repr__(self):
        ...

    @classmethod
    def format_energy_prod_key(cls, key):
        ...

    def format_expr(self):
        """
        Formats in a Mathematica-ingestible format
        :return:
        """
        ...

    @classmethod
    def shift_key(cls, key, shift):
        ...

    def shift_energies(self, shift):
        ...

    def shift(self, shift, shift_energies=False):
        ...

    def scale(self, scaling):
        ...

    @classmethod
    def get_key_ndim(cls, terms: dict):
        ...

    @property
    def ndim(self):
        ...

    def audit(self, target=None, ignore_constants=True):
        ...

    def ensure_dimension(self, ndim):
        ...

    def sort(self):
        ...

    def _permute_changes(self, changes, new_inds):
        ...

    def permute(self, new_inds, check_perm=True, allow_padding=False):
        ...

    @staticmethod
    def _check_neg(t1, t2):
        ...

    @staticmethod
    def find_term_scaling(key):
        ...

    def combine_energies(self):
        ...

    def combine(self, combine_subterms=True, combine_energies=False):
        ...

    @staticmethod
    def _permute_idx(idx, inds, rem, ndim):
        """
        Builds a new change index by aligning everything that matches in inds and
        incrementing everything that doesn't in remainder
        """
        ...

    @classmethod
    def _build_new_echange_key(cls, k1, k2, inds, remainder):
        ...

    @classmethod
    def _pad_echange_key_right(cls, k, ndim, inds, remainder):
        ...

    @classmethod
    def _pad_echange_key_left(cls, k, ndim, inds, remainder):
        ...

    def mul_along(self, other: 'PolynomialInterface', inds, remainder=None, mapping=None):
        """
        We multiply every subpoly along the given indices, transposing the appropriate tensor indices

        :param other:
        :param inds:
        :param remainder:
        :return:
        """
        ...

    def rmul_along(self, other, inds, remainder=None, mapping=None):
        """
        We multiply every subpoly along the given indices, transposing the appropriate tensor indices

        :param other:
        :param inds:
        :param remainder:
        :return:
        """
        ...

    def mul_simple(self, other: 'ProductPTPolynomial'):
        """
        We multiply every subpoly along the given indices, transposing the appropriate tensor indices

        :param other:
        :param inds:
        :param remainder:
        :return:
        """
        ...

    def rmul_simple(self, other: 'PolynomialInterface') -> 'PolynomialInterface':
        ...

class PTTensorCoeffProductSum(TensorCoefficientPoly, PolynomialInterface):
    """
    A representation for a sum of tensor coefficients * polynomial sums
    which is primarily here so we can transpose tensor coefficients intelligently
    """
    tensor_ids = ['V', 'G', 'U', 'Z', 'u', 'M', 'LV', 'LG']

    def __init__(self, terms, prefactor=1, canonicalize=True, ndim=None, inds_map=None, reduced=False):
        ...

    def mutate(self, terms: 'Any'=default, *, prefactor: 'Any'=default, ndim: 'Any'=default, inds_map: 'Any'=default, canonicalize: 'Any'=default, reduced: 'Any'=default):
        ...

    def prep_serialization_dict(self):
        ...

    @classmethod
    def from_serialization_dict(cls, big_dict):
        ...

    def to_state(self, serializer=None):
        ...

    @classmethod
    def from_state(cls, state, serializer=None):
        ...

    def flip_energy_terms(self):
        ...

    def filter(self, terms, mode='match'):
        ...

    def filter_coefficients(self, terms, mode='match'):
        ...

    def filter_energies(self, terms, mode='match'):
        ...

    @classmethod
    def format_key(cls, key):
        ...

    def __repr__(self):
        ...

    @classmethod
    def format_tensor_key(cls, key):
        ...

    def format_expr(self):
        """
        Formats in a Mathematica-ingestible format
        :return:
        """
        ...

    @classmethod
    def coeff_product_inds(cls, key, return_counts=False):
        ...

    def get_inds(self, key):
        ...

    def prune_operators(self, ops):
        ...

    @property
    def ndim(self):
        ...

    def audit(self, target=None, required_dimension=None, ignore_constants=True):
        """
        Checks to ensure that the number of dimensions aligns with
        the number of indices in the tensor coefficients

        :return:
        """
        ...

    def ensure_dimension(self, ndim):
        ...

    @classmethod
    def ensure_subpoly_dim(cls, terms):
        ...

    def sort(self):
        ...

    def permute(self, new_inds, check_perm=True, allow_padding=False):
        ...

    def _check_equiv(self, k1, k2):
        ...

    def combine_terms(self):
        ...

    @classmethod
    def reindex_terms(cls, terms):
        ...

    def combine(self, combine_coeffs=False, combine_subterms=True, combine_energies=False):
        ...

    @classmethod
    def _get_uperms(cls, p_list):
        ...

    def free_up_indices(self, start, stop):
        ...

    def shift_energies(self, change):
        ...

    def shift(self, shift):
        ...

    def scale(self, scaling):
        ...

    @staticmethod
    def _potential_symmetrizer(idx):
        ...

    @staticmethod
    def _kinetic_symmetrizer(idx):
        ...

    @staticmethod
    def _coriolis_symmetrizer(idx):
        ...
    default_symmetrizers = None

    @classmethod
    def symmetrizers(cls):
        ...

    @classmethod
    def _symmetrize(cls, idx, symmetrizers=None):
        ...

    @classmethod
    def canonical_key(cls, monomial_tuple, symmetrizers=None):
        ...

    @classmethod
    def _handle_mul(cls, new_key, total_order, num_fixed, mul_size, poly_1, poly_2, mul_inds, mul_remainder, og_remainder, left_baseline, baseline, left_choice_x, right_perm_x, left_rem_x, right_rem_x, num_left, num_right, logger, log_level):
        ...

    def _generate_direct_product_values(self, inds, remainder, index_classes, baseline, k1, k2, poly_1, poly_2):
        ...

    def _adjust_key_right(self, k, other, inds, remainder):
        ...

    def _adjust_key_left(self, k, other, inds, remainder):
        ...

    def mul_along(self, other: 'PolynomialInterface', inds, remainder=None, index_classes=None, mapping=None, baseline=None):
        """
        We multiply every subpoly along the given indices, transposing the appropriate tensor indices

        :param other:
        :param inds:
        :param remainder:
        :return:
        """
        ...

    def rmul_along(self, other: 'ProductPTPolynomial', inds, remainder=None, mapping=None):
        """
        We multiply every subpoly along the given indices, transposing the appropriate tensor indices

        :param other:
        :param inds:
        :param remainder:
        :return:
        """
        ...

    def mul_simple(self, other: 'ProductPTPolynomial'):
        """
        We multiply every subpoly along the given indices, transposing the appropriate tensor indices

        :param other:
        :param inds:
        :param remainder:
        :return:
        """
        ...

    def rmul_simple(self, other: 'PolynomialInterface') -> 'PolynomialInterface':
        ...

class SqrtChangePoly(PolynomialInterface):

    def __init__(self, poly_obj: 'PolynomialInterface', change, shift, canonicalize=False):
        ...

    def mutate(self, poly_obj: 'Any'=default, poly_change: 'Any'=default, shift_start: 'Any'=default, canonicalize=False):
        ...

    def to_state(self, serializer=None):
        ...

    @classmethod
    def from_state(cls, state, serializer=None):
        ...

    def strip(self):
        ...

    @property
    def ndim(self):
        ...

    def audit(self, target=None, ignore_constants=True):
        ...

    def ensure_dimension(self, ndim):
        ...

    def __repr__(self):
        ...

    def short_repr(self):
        ...

    def format_sqrt_expr(self):
        ...

    def format_expr(self):
        """
        Formats in a Mathematica-ingestible format
        :return:
        """
        ...

    def __add__(self, other):
        ...

    def combine(self, combine_coeffs=False, combine_subterms=True, combine_energies=False):
        ...

    def sort(self):
        ...

    @classmethod
    def canonicalize(cls, poly_obj, poly_change, shift_start):
        ...

    def shift_energies(self, changes):
        ...

    def shift(self, shift, shift_energies=False):
        ...

    def scale(self, scaling):
        ...

    def permute(self, perm, check_perm=True, allow_padding=False):
        ...

    def filter_coefficients(self, terms, mode='match'):
        ...

    def filter_energies(self, terms, mode='match'):
        ...

    def __radd__(self, other):
        ...

    def __mul__(self, other):
        ...

    def __rmul__(self, other):
        ...

    def get_change_poly(self, initial_changes, extra_changes, initial_shift, extra_shift):
        ...

    def mul_simple(self, other: 'PolynomialInterface'):
        ...

    def rmul_simple(self, other: 'PolynomialInterface') -> 'PolynomialInterface':
        ...

    def mul_along(self, other, inds, remainder=None, mapping=None, baseline=None):
        ...

    def rmul_along(self, other, inds, remainder=None, mapping=None):
        ...

class PerturbationTheoryTerm(metaclass=abc.ABCMeta):
    """
    A generic version of one of the three terms in
    PT that will generate a correction polynomial
    """
    use_intermediate_normalization = False

    def __init__(self, logger=None, checkpoint=None, allowed_terms=None, allowed_energy_changes=None, intermediate_normalization=None, allowed_coefficients=None, disallowed_coefficients=None):
        ...

    def get_subexpressions(self) -> 'Iterable[PerturbationTheoryTerm]':
        ...

    def __mul__(self, other):
        ...

    def __rmul__(self, other):
        ...

    def __neg__(self):
        ...

    @property
    def expressions(self):
        ...

    @staticmethod
    def change_sort_key(changes):
        ...

    @classmethod
    def change_sort(cls, changes):
        ...

    @classmethod
    def sorted_changes(cls, changes):
        ...

    @abc.abstractmethod
    def get_changes(self) -> 'dict[tuple[int], Any]':
        ...

    @property
    def changes(self):
        ...

    def get_serializer_key(self):
        ...

    @property
    def serializer_key(self):
        ...
    debug_logger = None
    _logger_stack = []

    @contextlib.contextmanager
    def debug_logging(self):
        ...

    @classmethod
    def default_logger(cls) -> Logger:
        ...
    _default_filters = None
    _filter_stack = []

    @contextlib.contextmanager
    def default_filtering(self):
        ...

    @classmethod
    def default_filters(cls) -> Logger:
        ...

    def get_core_poly(self, changes, shift=None):
        ...
    simplify_by_default = True
    caching_enabled = True

    def get_poly_terms(self, changes, simplify=None, shift=None) -> 'SqrtChangePoly':
        ...

    def __call__(self, changes, shift=None, coeffs=None, freqs=None, check_sorting=True, simplify=None, return_evaluator=True):
        ...

class OperatorExpansionTerm(PerturbationTheoryTerm):

    def __init__(self, terms, order=None, identities=None, symmetrizers=None, index=None, allowed_terms=None, change_rules=None, **opts):
        ...

    def __repr__(self):
        ...
    _generator_cache = {}

    @classmethod
    def _get_generator(cls, term_list, inds):
        ...

    @classmethod
    def _evaluate_poly_coeffs(cls, terms, inds, delta, shift, sqrt_scale=False):
        ...

    @property
    def change_rules(self):
        ...

    @change_rules.setter
    def change_rules(self, rules):
        ...

    def get_changes(self) -> 'dict[tuple[int], Any]':
        ...

    def get_subexpressions(self) -> 'Iterable[PerturbationTheoryTerm]':
        ...

    def get_core_poly(self, changes, shift=None) -> 'SqrtChangePoly':
        ...
    keep_ints = True
    _poly_cache = {}

    @classmethod
    def _resolve_poly(cls, term_list, partition_sizes, p_index, c_key, s_key, p_vec, changes, shift):
        ...

class HamiltonianExpansionTerm(OperatorExpansionTerm):

    def __init__(self, terms, order=None, identities=None, symmetrizers=None, change_rules=None, **opts):
        ...

    def __repr__(self):
        ...

class PerturbationOperator(PerturbationTheoryTerm):
    _energy_baseline = None

    def __init__(self, subterm):
        ...
    _cache = {}

    @classmethod
    def lookup(cls, subterm):
        ...

    def __repr__(self):
        ...

    def get_changes(self) -> 'dict[tuple[int], Any]':
        ...

    def get_core_poly(self, changes, shift=None):
        ...

class _ShiftedEnergyBaseline(PerturbationTheoryTerm):

    def __init__(self, base_term: 'PerturbationTheoryTerm'):
        ...

    def __repr__(self):
        ...

    def get_changes(self):
        ...

    def get_core_poly(self, changes, shift=None) -> 'SqrtChangePoly':
        ...

class ShiftedEnergyBaseline(PerturbationTheoryTerm):
    """
    Represents a term that will be multipled by on the left rather than the right
    for evaluating things like Y[1]M[0]Y[1], essentially changing raising operations to lowering
    """

    def __init__(self, base_term):
        ...

    def __repr__(self):
        ...
    _cache = {}

    @classmethod
    def lookup(cls, base_term):
        ...

    def get_changes(self):
        ...

    def get_poly_terms(self, changes, shift=None, **opts) -> 'SqrtChangePoly':
        ...

class ShiftedHamiltonianCorrection(PerturbationTheoryTerm):
    """
       Adds the wave function correction and the overlap term
       """

    def __init__(self, parent, order, allowed_terms=None, **opts):
        ...

    def get_serializer_key(self):
        ...

    def __repr__(self):
        ...

    def get_changes(self):
        ...

    def get_subexpressions(self) -> 'Iterable[PerturbationTheoryTerm]':
        ...

class WavefunctionCorrection(PerturbationTheoryTerm):

    def __init__(self, parent, order, allowed_terms=None, **opts):
        ...

    def get_serializer_key(self):
        ...

    def __repr__(self):
        ...

    def get_changes(self):
        ...

    def get_subexpressions(self):
        ...

class EnergyCorrection(PerturbationTheoryTerm):

    def __init__(self, parent, order, allowed_terms=None, **opts):
        ...

    def get_serializer_key(self):
        ...

    def __repr__(self):
        ...

    def get_changes(self) -> 'dict[tuple[int], Any]':
        ...

    def get_subexpressions(self) -> 'Iterable[PerturbationTheoryTerm]':
        ...

    def get_poly_terms(self, changes, shift=None, **opts):
        ...

    def __call__(self, changes, shift=None, coeffs=None, freqs=None, check_sorting=None, simplify=True):
        ...

class WavefunctionOverlapCorrection(PerturbationTheoryTerm):

    def __init__(self, parent, order, allowed_terms=None, degenerate_changes=None, **opts):
        ...

    def get_serializer_key(self):
        ...

    def __repr__(self):
        ...

    def get_changes(self) -> 'dict[tuple[int], Any]':
        ...

    def get_subexpressions(self) -> 'Iterable[PerturbationTheoryTerm]':
        ...

class FullWavefunctionCorrection(PerturbationTheoryTerm):
    """
    Adds the wave function correction and the overlap term
    """

    def __init__(self, parent, order, allowed_terms=None, **opts):
        ...

    def get_serializer_key(self):
        ...

    def __repr__(self):
        ...

    def get_changes(self):
        ...

    def get_subexpressions(self) -> 'Iterable[PerturbationTheoryTerm]':
        ...

class OperatorCorrection(PerturbationTheoryTerm):

    def __init__(self, parent, order, operator_type=None, allowed_terms=None, wavefunction_generator=None, base_index=None, **opts):
        ...

    @classmethod
    def get_type_key(cls, operator_type):
        ...

    def get_serializer_key(self):
        ...
    repr_key = 'M'

    def get_repr_key(self):
        ...

    def __repr__(self):
        ...

    def get_changes(self):
        ...

    def default_wavefunction_generator(self, o: int):
        ...

    def wavefunction_generator(self, o: int):
        ...

    def get_subexpressions(self, bra_wavefunction_generator=None, ket_wavefunction_generator=None, bounds=None, const_zeros=None) -> 'Iterable[PerturbationTheoryTerm]':
        ...

class OperatorDegenerateCorrection(OperatorCorrection):

    def __init__(self, parent, order, degenerate_changes=None, operator_type=None, allowed_terms=None, **opts):
        ...

    def __repr__(self):
        ...

    def default_wavefunction_generator(self, o: int):
        ...

    class Left(OperatorCorrection):

        def __init__(self, parent: OperatorCorrection, order, real_parent):
            ...

        def __repr__(self):
            ...

        def get_subexpressions(self):
            ...

    class Right(OperatorCorrection):

        def __init__(self, parent, order, real_parent):
            ...

        def __repr__(self):
            ...

        def get_subexpressions(self):
            ...

    class Both(OperatorCorrection):

        def __init__(self, parent, order, real_parent):
            ...

        def __repr__(self):
            ...

        def get_subexpressions(self):
            ...

    def get_subexpressions(self, bra_wavefunction_generator=None, ket_wavefunction_generator=None) -> 'Iterable[PerturbationTheoryTerm]':
        ...

    def get_right_degenerate_expressions(self, bra_wavefunction_generator=None, ket_wavefunction_generator=None) -> 'Iterable[PerturbationTheoryTerm]':
        ...

    def get_left_degenerate_expressions(self, bra_wavefunction_generator=None, ket_wavefunction_generator=None) -> 'Iterable[PerturbationTheoryTerm]':
        ...

    def get_both_degenerate_expressions(self, bra_wavefunction_generator=None, ket_wavefunction_generator=None) -> 'Iterable[PerturbationTheoryTerm]':
        ...

class DiagonalHamiltonian(OperatorExpansionTerm):

    def __init__(self):
        ...

    def __repr__(self):
        ...

    def get_changes(self):
        ...

class ReexpressedHamiltonian(OperatorCorrection):
    repr_key = 'H'

    def __init__(self, parent, order, allowed_terms=None, degenerate_changes=None, hamiltonian_corrections=None, **opts):
        ...

    @classmethod
    def prep_expansion(self, base_expansion, corrs, *, logger, allowed_coefficients, disallowed_coefficients, base_index=5):
        ...

    def get_repr_key(self):
        ...

class ReexpressedHamiltonianDegenerateCorrection(OperatorDegenerateCorrection):
    repr_key = 'H'

    def __init__(self, parent, order, allowed_terms=None, hamiltonian_corrections=None, **opts):
        ...

    def get_repr_key(self):
        ...

class ScaledPerturbationTheoryTerm(PerturbationTheoryTerm):

    def __init__(self, base_term: 'PerturbationTheoryTerm', scaling):
        ...
    _cache = {}

    @classmethod
    def lookup(cls, base_term, scaling):
        ...

    def __repr__(self):
        ...

    def get_changes(self) -> 'dict[tuple[int], Any]':
        ...

    def get_core_poly(self, changes, shift=None) -> 'SqrtChangePoly':
        ...

class PerturbationTheoryTermSum(PerturbationTheoryTerm):

    def __init__(self, *terms):
        ...

    def __repr__(self):
        ...

    def get_changes(self):
        ...

    def get_subexpressions(self):
        ...

class PerturbationTheoryTermProduct(PerturbationTheoryTerm):
    _cache = {}

    def __init__(self, post_op, pre_op):
        ...

    @classmethod
    def lookup(cls, post_op, pre_op):
        ...

    def __repr__(self):
        ...
    _change_map = {}
    _perms_map = {}

    @classmethod
    def get_change_classes(cls, changes):
        ...
    _combination_inds = {}
    _combination_comps = {}

    @classmethod
    def get_combination_inds(cls, n, r):
        ...

    @classmethod
    def get_combination_comp(cls, n, r):
        ...

    @staticmethod
    def _fill_change_data(full_changes, n, new_changes, init_changes, subsorts, contacts, rems):
        ...
    _permutation_cache = MaxSizeCache(128, cache_type='fifo')

    @classmethod
    def _get_uperms(cls, p_list):
        ...

    @classmethod
    def _fill_class_counts(cls, data, which, c):
        ...

    @classmethod
    def _fill_inds_data(cls, data, which, counts, mapping, m, permute):
        ...

    @classmethod
    def _prep_product_change_data(cls, i1, ch1, inds_1, rind_1, n1, i2, ch2, inds_2, rind_2, n2, os, pad_change, initial_sort, nz, x1, x2, upairs):
        ...

    def get_changes(self):
        ...

    def get_expressions(self):
        ...

    @classmethod
    def get_poly_product_terms(cls, gen1, gen2, change_1, change_2, target_inds, remainder_inds, reorgs, simplify=True):
        ...

    def get_core_poly(self, changes, shift=None):
        ...

class PerturbationTheoryExpressionEvaluator:

    def __init__(self, op, expr: 'SqrtChangePoly', change=None, logger=None, parallelizer=None):
        ...

    def __repr__(self):
        ...

    @staticmethod
    def _extract_coeffs(coeffs, coeff_indices, perm):
        ...

    @classmethod
    def _eval_raw_poly(cls, perm_substates, poly, change, baseline_shift, pows, verbose, logger):
        ...

    @classmethod
    def _eval_poly(cls, cache, tuple_states, perm_substates, pows, poly, change, baseline_shift, verbose, logger):
        ...

    @classmethod
    def _compute_energy_weights(cls, energy_changes, perm_freqs):
        ...

    @classmethod
    def _get_prefacs(cls, perms, cinds_remapped, ctensors, counts_cache, facs, zero_cutoff):
        ...

    @classmethod
    def _eval_subcontrib(cls, echanges, perm_freqs, poly_cache, energy_cache, tuple_states, perm_substates, pows, polys, change, baseline_shift, verbose, logger):
        ...

    @classmethod
    def _test_degs(cls, subchanges, echanges, perm_subsets, only_deg):
        ...

    @classmethod
    def _eval_perm_core(cls, expr, state, tuple_states, perm_substates, which_perms, change, baseline_shift, prefacs, perm_freqs, pows, key, perm_subsets, degenerate_changes, only_degenerate_terms, poly_cache, energy_cache, verbose, logger, log_level, log_scaling):
        ...

    @classmethod
    def _get_state_perms(cls, state_idx, state, freqs, fixed, subset, sub_perms, take_cache, mask_pos, full_set):
        ...

    @classmethod
    def _eval_perm(cls, expr, change, baseline_shift, subset, state_perms, all_perms, perm_map, freqs, cind_specs, ctensors, num_fixed, degenerate_changes, only_degenerate_terms, zero_cutoff, counts_cache, poly_cache, take_cache, energy_cache, facs, pows, split_spec, verbose, logger, log_scaled):
        ...

    @classmethod
    def _prep_coeff_data(cls, cind_sets, coeff_lists):
        ...

    @classmethod
    def _get_max_order(cls, expr):
        ...

    @classmethod
    def _deg_test_direct(cls, state, modes, check_pos):
        ...
    _direct_degs = collections.namedtuple('_direct_degs', ['check_pos', 'modes'])

    @classmethod
    def _deg_test(cls, modes):
        ...

    @classmethod
    def _deg_test_tree(cls, state, tree):
        ...

    @staticmethod
    def _mode_inclusion_direct(modes, tests):
        ...

    @staticmethod
    def _mode_inclusion_multi(modes, tests):
        ...

    @classmethod
    def _compile_deg_tests(cls, tests):
        ...

    @classmethod
    def _make_full_deg_test(cls, tests):
        ...
    default_deg_id_method = 'linear'
    deg_id_method_nmodes_switch = 50

    @classmethod
    def _identify_possible_degeneracies(cls, all_perms, sides, expr, changes, nmodes, method=None):
        ...

    @classmethod
    def set_cache_size(cls, max_size):
        ...

    @classmethod
    def get_cache(self):
        ...
    _max_cache_size = 10000000.0
    _poly_cache = MaxSizeCache(_max_cache_size, cache_type='fifo')
    _ecoeff_cache = MaxSizeCache(_max_cache_size, cache_type='fifo')
    default_zero_cutoff = 1e-18
    _parallel_eval_main_args = None
    _cached_expansion = None

    @classmethod
    def _initialize_main_eval(cls, coeffs, *args, parallelizer=None):
        ...
    _cached_main_args = None

    @classmethod
    def _extract_main_args(cls):
        ...

    @classmethod
    def _run_eval_block(cls, block_inds, contrib_shapes, block_size, free_inds, cind_sets):
        ...

    @classmethod
    def _run_eval_blocks(cls, blocks, contrib_shapes, *args, parallelizer: Parallelizer):
        ...

    @classmethod
    def evaluate_polynomial_expression(cls, state_perms, coeffs, freqs, expr, change, baseline_shift, num_fixed, op=None, logger=None, parallelizer=None, degenerate_changes=None, only_degenerate_terms=False, zero_cutoff=None, verbose=False, log_scaled=True):
        ...

    def evaluate(self, state_perms, coeffs, freqs, degenerate_changes=None, only_degenerate_terms=False, zero_cutoff=None, parallelizer=None, verbose=False, log_scaled=True):
        ...

class PerturbationTheoryEvaluator:

    def __init__(self, solver: AnalyticPerturbationTheorySolver, expansion, freqs=None, parallelizer=None):
        ...

    def modify_hamiltonian(self, hamiltonian_corrections):
        ...

    def get_energy_corrections(self, states, order=None, expansions=None, freqs=None, zero_cutoff=None, degenerate_states=None, verbose=False, logger=None, parallelizer=None):
        ...

    @staticmethod
    def is_single_expansion(expansion, min_order=None):
        ...

    def _prep_expansions(self, expansions):
        ...

    def get_overlap_corrections(self, states, order=None, expansions=None, degenerate_states=None, freqs=None, zero_cutoff=None, verbose=False, parallelizer=None):
        ...

    def get_diff_map(self, state_map):
        ...

    @staticmethod
    def get_finals(initial, change, perms):
        ...

    def _reformat_corrections(self, order, corrs, change_map, num_expansions):
        ...

    @classmethod
    def _compute_corr(cls, gen_corrs, only_degs, key, gen, expr, *args, only_degenerate_terms=None, **kwargs):
        ...

    @classmethod
    def _build_corrections(cls, corr_gen, degenerate_corr_gen, expansions, order, terms, allowed_coefficients, disallowed_coefficients, epaths, change_map, degenerate_changes, only_degenerate_terms, include_degenerate_correction_terms, freqs, verbose, logger, zero_cutoff, log_scaled, parallelizer):
        ...

    @classmethod
    def get_degenerate_changes(cls, degenerate_pairs):
        ...

    def _sort_corrections(self, corrs: BasicAPTCorrections, state_pairs):
        ...

    def get_state_by_state_corrections(self, generator, states, order=None, terms=None, epaths=None, expansions=None, freqs=None, verbose=False, allowed_coefficients=None, disallowed_coefficients=None, degenerate_states=None, only_degenerate_terms=False, degenerate_correction_generator=None, include_degenerate_correction_terms=True, log_scaled=False, zero_cutoff=None, return_sorted=False, logger=None, parallelizer=None):
        ...

    def get_matrix_corrections(self, states, order=None, expansions=None, freqs=None, zero_cutoff=None, verbose=False):
        ...

    def get_full_wavefunction_corrections(self, states, order=None, expansions=None, freqs=None, zero_cutoff=None, degenerate_states=None, verbose=False):
        ...

    def get_wavefunction_corrections(self, states, order=None, expansions=None, freqs=None, zero_cutoff=None, degenerate_states=None, verbose=False):
        ...

    @classmethod
    def _prep_ham_corrs(cls, base_corrs):
        ...

    def get_reexpressed_hamiltonian(self, states, order=None, expansions=None, freqs=None, degenerate_states=None, only_degenerate_terms=False, verbose=False, include_diagonal=False, hamiltonian_corrections=None, **opts):
        """

        :param states:
        :param order:
        :param expansions:
        :param freqs:
        :param degenerate_states:
        :param only_degenerate_terms:
        :param verbose:
        :param include_diagonal:
        :param hamiltonian_corrections:  `[[(order, terms), expansion], ...]`
        :param opts:
        :return:
        """
        ...

    def _prep_operator_expansion(self, expansions, operator_expansion):
        ...

    def get_operator_corrections(self, operator_expansion, states, order=None, expansions=None, freqs=None, degenerate_states=None, operator_type=None, check_single=True, terms=None, min_order=1, verbose=False, **opts):
        ...

    def evaluate_expressions(self, states, exprs, expansions=None, operator_expansions=None, degenerate_states=None, zero_cutoff=None, verbose=False):
        ...