import abc, numpy as np, functools
from McUtils.Data import AtomData, UnitsData
import McUtils.Numputils as nput
from ..Molecools import StructuralProperties
from .Evaluators import *
from .Interpolation import *
__all__ = ['DGBCoords', 'DGBCartesians', 'DGBInternals', 'DGBWatsonModes']

class DGBCoords(metaclass=abc.ABCMeta):

    @property
    @abc.abstractmethod
    def centers(self) -> 'np.ndarray':
        ...

    @property
    def shape(self):
        ...

    @property
    @abc.abstractmethod
    def kinetic_energy_evaluator(self) -> 'DGBKineticEnergyEvaluator':
        ...

    @property
    @abc.abstractmethod
    def pairwise_potential_evaluator_type(self) -> 'type[DGBPairwisePotentialEvaluator]':
        ...

    def pairwise_potential_evaluator(self, potential_functions) -> 'DGBPairwisePotentialEvaluator':
        ...

    @abc.abstractmethod
    def __getitem__(self, item) -> 'DGBCoords':
        ...

    def take_indices(self, subinds) -> 'DGBCoords':
        ...

    def drop_indices(self, subinds) -> 'DGBCoords':
        ...

    @abc.abstractmethod
    def gmatrix(self, coords: np.ndarray) -> np.ndarray:
        ...

    @classmethod
    def embedded_mode_function(cls, func, modes, masses=None):
        ...

    @classmethod
    def embedded_subcoordinate_function(cls, func, sel, ndim):
        ...

    @classmethod
    def embedded_cartesian_function(cls, func, atom_sel, xyz_sel, natoms, ndim):
        ...

    class DGBEmbeddedFunction:

        def __init__(self, embedded_function, original_function, coords):
            ...

        def __call__(self, coords, deriv_order=None):
            ...

    @abc.abstractmethod
    def embed_function(self, fn) -> 'DGBEmbeddedFunction':
        ...

    def as_cartesians(self) -> 'tuple[DGBCartesians, tuple[np.ndarray, np.ndarray]]':
        ...

class DGBCartesians(DGBCoords):

    def __init__(self, coords, masses, *, natoms=None, atom_sel=None, ndim=None, xyz_sel=None):
        ...

    @property
    def centers(self):
        ...

    @property
    def cart_shape(self):
        ...

    @property
    def kinetic_energy_evaluator(self):
        ...

    @property
    def pairwise_potential_evaluator_type(self):
        ...

    @classmethod
    def resolve_masses(cls, coords, masses=None, atoms=None):
        ...

    @classmethod
    def from_cartesians(cls, centers, masses=None, atoms=None):
        ...

    def infer_shape_sel(self, selector):
        ...

    @staticmethod
    def _merge_sel(new_sel, old_sel):
        ...

    def __getitem__(self, item) -> 'DGBCartesians':
        ...

    def take_indices(self, subinds):
        ...

    def embed_function(self, function):
        """
        Embeds assuming we got a function in Cartesians _before_ any selections happened

        :param function:
        :return:
        """
        ...

    def gmatrix(self, coords: np.ndarray) -> np.ndarray:
        ...

class DGBInternals(DGBCoords):

    def __init__(self, coords, gmat_function=None, vprime_function=None):
        ...

class DGBWatsonModes(DGBCoords):

    def __init__(self, coords, modes, *, coriolis_inertia_function=None, masses=None, subselection=None):
        ...

    @property
    def centers(self):
        ...

    @property
    def kinetic_energy_evaluator(self):
        ...

    @property
    def pairwise_potential_evaluator_type(self):
        ...

    @staticmethod
    def zeta_momi(watson_coords, modes, masses):
        ...

    @classmethod
    def default_coriolis_inertia_function(cls, modes, masses):
        ...

    @classmethod
    def embed_coords(cls, carts, modes, shift=True):
        ...

    @classmethod
    def unembed_coords(cls, mode_coords, modes, masses=None, shift=True):
        ...

    @classmethod
    def embed_derivs(cls, derivs, modes):
        ...

    @classmethod
    def from_cartesians(cls, coords, modes, masses=None, coriolis_inertia_function=None):
        ...

    def as_cartesians(self, masses=None) -> 'tuple[DGBCartesians, tuple[np.ndarray, np.ndarray]]':
        ...

    def __getitem__(self, item):
        ...

    def gmatrix(self, coords: np.ndarray) -> np.ndarray:
        ...

    def embed_function(self, fn):
        ...