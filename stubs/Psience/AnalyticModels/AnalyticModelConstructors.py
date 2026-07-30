import itertools
import numpy as np, scipy, math
from McUtils.Data import AtomData, UnitsData
import McUtils.Numputils as nput
from .Helpers import AnalyticModelBase
from .AnalyticJacobianDotCalculator import InternalJacobianDisplacements
from ..Data import KEData
__all__ = ['AnalyticPotentialConstructor', 'AnalyticKineticEnergyConstructor', 'AnalyticModel', 'MolecularModel', 'GeometricFunction']
__reload_hook__ = ['..Data']

class AnalyticPotentialConstructor(AnalyticModelBase):
    """
    Provides a set of symbolic potentials for use in models

    :related:AnalyticModel, AnalyticKineticEnergyConstructor
    """

    @classmethod
    def morse(cls, i, j, De=None, a=None, re=None, eq=None, w=None, wx=None):
        """
        Returns a fully symbolic form of a Morse potential
        :return:
        :rtype:
        """
        ...

    @staticmethod
    def calc_morse(De, a, r, re):
        ...

    @staticmethod
    def harm(k, x, x_e):
        ...

    @classmethod
    def harmonic(cls, *args, k=None, eq=None, qe=None):
        """
        Returns a fully symbolic form of a Morse potential
        :return:
        :rtype:
        """
        ...

    @staticmethod
    def lin(k, x, x_e):
        ...

    @classmethod
    def linear(cls, *args, k=1, eq=None, xe=None):
        """
        Returns a fully symbolic form of a linear function
        :return:
        :rtype:
        """
        ...

    @staticmethod
    def pow(k, x, x_e, n):
        ...

    @classmethod
    def power(cls, *args, k=1, eq=None, n=None, xe=None):
        """
        Returns a fully symbolic form of a linear function
        :return:
        :rtype:
        """
        ...

    @classmethod
    def cos(cls, *args, eq=None, qe=None):
        """
        Returns a fully symbolic form of cos
        :return:
        :rtype:
        """
        ...

    @classmethod
    def sin(cls, *args, eq=None, qe=None):
        """
        Returns a fully symbolic form of sin
        :return:
        :rtype:
        """
        ...

    @classmethod
    def multiwell(cls, *args, turning_points=None, origin=None, eq=None, minimum=0, depth=None):
        """

        :param args:
        :type args:
        :param turning_points:
        :type turning_points:
        :param depth:
        :type depth:
        :return:
        :rtype:
        """
        ...

class GeometricFunction:

    def __init__(self, coords, coord_funcs, val_func):
        ...

    def __call__(self, masses, coords):
        ...

    @classmethod
    def position_function(cls, i):
        ...

    @classmethod
    def normal_function(cls, i, j, k):
        ...

    @classmethod
    def mass_function(cls, i):
        ...

    @classmethod
    def distance_function(cls, i, j):
        ...

    @classmethod
    def angle_function(cls, i, j, k):
        ...

    @classmethod
    def dihedral_function(cls, i, j, k, l):
        ...

    @classmethod
    def book_function(cls, i, j, k, l):
        ...

    @classmethod
    def create_symbol_function(cls, sym):
        ...

    @classmethod
    def sorted_vars(cls, vars):
        ...

    @classmethod
    def from_expr(cls, expr):
        ...

class AnalyticKineticEnergyConstructor(AnalyticModelBase):
    """
    Provides G and V' elements from Frederick and Woywood

    :related:AnalyticModel, AnalyticPotentialConstructor
    """

    @classmethod
    def _remap_inter_inds(cls, inds, intersection, remapping, n_inter, n_diff):
        ...

    @classmethod
    def _canonicalize_inds(cls, t, inds):
        ...

    @classmethod
    def _simple_presorts(cls, coord_types, inds1, inds2):
        ...

    @classmethod
    def _get_coord_key(cls, inds1: 'Iterable[int]', inds2: 'Iterable[int]', coord_types=None):
        ...

    @classmethod
    def kinetic_exprs(cls, inds1: 'Iterable[int]', inds2: 'Iterable[int]', coord_types=None, target_symbols=None):
        ...

    @classmethod
    def _apply_remapping(cls, expr, mapping, target_symbols):
        ...

    @classmethod
    def _g_expr_direct(cls, coord_types, inds1: 'Iterable[int]', inds2: 'Iterable[int]'):
        ...

    @classmethod
    def _vp_expr_direct(cls, coord_types: '[str,str]', inds1: 'Iterable[int]', inds2: 'Iterable[int]', G):
        ...

    @classmethod
    def kinetic_exprs_direct(cls, inds1: 'Iterable[int]', inds2: 'Iterable[int]', coord_types=None, return_vp=True, target_symbols=None):
        """
        Evaluated using the simple expressions in Table 1 from Frederick and Woywood
        :param inds1:
        :param inds2:
        :param coord_types:
        :return:
        """
        ...

    @classmethod
    def g(cls, inds1: 'Iterable[int]', inds2: 'Iterable[int]', coord_types=None, target_symbols=None, return_function=False, method='lookup'):
        ...

    @classmethod
    def vp(cls, inds1: 'Iterable[int]', inds2: 'Iterable[int]', coord_types=None, target_symbols=None, return_function=False, method='lookup'):
        ...

    @classmethod
    def g_matrix(cls, coord_specs, return_function=False, return_matrix=True, method='lookup'):
        ...

    @classmethod
    def _g(cls, inds1: 'Iterable[int]', inds2: 'Iterable[int]'):
        ...

class AnalyticModel:
    """
    Provides a symbolic representation of an analytically evaluatable Hamiltonian
    which can be used to get derived expressions to evaluate.
    Supplies methods to automatically run DVR and VPT calculations from the model
    specifications as well.

    :related:AnalyticPotentialConstructor, AnalyticPotentialConstructor
    """

    def __init__(self, coordinates, potential, dipole=None, values=None, rotation=None):
        ...

    @classmethod
    def from_potential(cls, potential, dipole=None, values=None, rotation=None):
        ...

    @property
    def internal_coordinates(self):
        ...

    @property
    def constants(self):
        ...

    def normal_modes(self, dimensionless=True):
        ...

    def to_normal_modes(self, dimensionless=True):
        ...

    def get_VPT_expansions(self, order=2, expansion_order=None, include_potential=None, include_gmatrix=None, include_pseudopotential=None, evaluate=True):
        ...

    def run_VPT(self, order=2, states=2, return_analyzer=True, expansion_order=None, include_potential=None, include_gmatrix=None, include_pseudopotential=None, atoms=None, coords=None, **kwargs):
        ...

    def _get_rotated_coordinates(self):
        ...

    def _get_inverse_coordinates(self):
        ...

    class SympyExpr:

        def __init__(self, expr, core, ndim):
            ...

        def _broadcast_tree(self, shape, expr_list):
            ...

        def _array_call(self, grid, core):
            ...

        def _rec_call(self, grid, core, depth=0, rec_counter=None):
            ...

        def __call__(self, grid, vector=None, **kwargs):
            ...

        def __repr__(self):
            ...

        @classmethod
        def _lambdiy(cls, vars, expr, lambdify_arrays=False):
            ...

        @classmethod
        def lambdify(cls, vars, expr, ndim, lambdify_arrays=False):
            ...

    @classmethod
    def prep_lambda_expr(cls, base_coords, expr, dummify=False, rewrite_trig=True):
        ...

    @classmethod
    def lambdify(cls, coord_vec, expr, coordinate_transform=None, mode=None, dummify=False, rewrite_trig=True, lambdify_arrays=False):
        ...

    def wrap_function(self, expr, transform_coordinates=True, mode=None):
        ...

    def expand_potential(self, order, lambdify=True, evaluate=True, contract=True):
        ...

    def get_DVR_parameters(self, expansion_order=None, lambdify=True, evaluate='constants'):
        ...

    def setup_DVR(self, domain=None, divs=None, use_normal_modes=False, expansion_order=None, potential_function=None, **params):
        ...

    def evaluate(self, expr, mode='all', numericize=False):
        ...

    def _load_symbols(self):
        ...

    @classmethod
    def parse_symbol(self, sym):
        ...

    def jacobian(self, order=0, evaluate=False, lambdify=False):
        ...

    def jacobian_inverse(self, order=0, evaluate=False):
        ...

    def _base_gmat(self):
        ...

    def g(self, order=0, evaluate=False, lambdify=False):
        ...

    def v(self, order=2, evaluate=False, lambdify=False):
        ...

    def _base_u(self):
        ...

    def vp(self, order=0, evaluate=False, lambdify=False):
        ...

    def mu(self, order=1, evaluate=False, lambdify=False):
        ...
    Potential = AnalyticPotentialConstructor
    morse = AnalyticPotentialConstructor.morse
    harmonic = AnalyticPotentialConstructor.harmonic
    linear = AnalyticPotentialConstructor.linear
    power = AnalyticPotentialConstructor.power
    sin = AnalyticPotentialConstructor.sin
    cos = AnalyticPotentialConstructor.cos
    KE = AnalyticKineticEnergyConstructor

    class NamespaceContext:

        def __init__(self, context=None):
            ...

        def _get_frame_vars(self):
            ...

        def insert_vars(self):
            ...

        def prune_vars(self):
            ...

        def __enter__(self):
            ...

        def __exit__(self, exc_type, exc_val, exc_tb):
            ...

    @classmethod
    def sym(self, base, *args):
        ...

    @classmethod
    def m(self, i):
        ...

    @classmethod
    def r(self, i, j):
        ...

    @classmethod
    def a(self, i, j, k):
        ...

    @classmethod
    def t(self, i, j, k, l):
        ...

    @classmethod
    def y(self, i, j, k, l):
        ...
    convert = UnitsData.convert

    @staticmethod
    def mass(atom):
        ...

    def molecular_potential(self, mol):
        ...

    def molecular_dipole(self, mol):
        ...

    def molecular_gmatrix(self, mol):
        ...

class MolecularModel(AnalyticModel):

    def __init__(self, mol, coords, potential, dipole=None, values=None, rotation=None):
        ...

    @property
    def potential(self):
        ...

    @property
    def gmatrix(self):
        ...

    @property
    def vprime(self):
        ...

    @property
    def dipole(self):
        ...

    def setup_AIMD(self, timestep=0.5, initial_energies=None, initial_displacements=None, displaced_coords=None, track_kinetic_energy=False, track_velocities=False, **aimd_opts):
        ...

    def setup_DGB(self, centers, *, masses=None, modes='normal', transformations=None, alphas='auto', cartesians=None, potential_function=None, dipole_function=None, optimize_centers=None, quadrature_degree=None, expansion_degree=None, pairwise_potential_functions=None, internals=False, logger=True, **opts):
        ...

class MolecularModelFunction:

    def __init__(self, deriv_evaluator, mol):
        ...

    def evaluate(self, carts, deriv_order=None, internals=False, which=None, sel=None, axes=None, derivs=None):
        ...

    def __call__(self, carts, deriv_order=None, internals=False, which=None, sel=None, axes=None):
        ...

class MolecularModelPotentialFunction(MolecularModelFunction):

    def __init__(self, model, mol):
        ...

class MolecularModelDipoleFunction(MolecularModelFunction):

    def __init__(self, model, mol):
        ...

    def evaluate(self, carts, deriv_order=None, internals=False, which=None, sel=None, axes=None):
        """
        This has the added complication of needing to dispatch over the axes...

        :param carts:
        :type carts:
        :param deriv_order:
        :type deriv_order:
        :param internals:
        :type internals:
        :param which:
        :type which:
        :param sel:
        :type sel:
        :param axes:
        :type axes:
        :return:
        :rtype:
        """
        ...

class MolecularModelGMatrixFunction(MolecularModelFunction):

    def __init__(self, model, mol):
        ...

class MolecularModelVPrimeFunction(MolecularModelFunction):

    def __init__(self, model, mol):
        ...