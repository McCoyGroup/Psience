"""
Temporary storage for an attempt to get G matrix elements purely symbolically
but which lacked sufficient information about the coordinates involved
"""
import numpy as np
from .Helpers import sym, AnalyticModelBase

class OrientationVector:

    def __init__(self, x, y, z, basis):
        ...

    def __neg__(self):
        ...

    def dot(self, other: 'OrientationVector'):
        ...

    @classmethod
    def rotation_matrix(cls, angle, axis: 'OrientationVector'):
        ...

    def __repr__(self):
        ...

class RotationMatrix:

    def __init__(self, v1: OrientationVector, v2: OrientationVector, v3: OrientationVector):
        ...

    def dot(self, v):
        ...

    def __neg__(self):
        ...

class BondVector:
    logger = NullLogger()

    def __init__(self, i, j, norm=1, embedding=None):
        ...

    def one_common_atom(self, bond: 'BondVector', mode='shared'):
        ...

    def shared_atom(self, bond: 'BondVector'):
        ...

    def unshared_atom(self, bond: 'BondVector'):
        ...

    @property
    def direction(self):
        ...

    def __repr__(self):
        ...

    def __neg__(self):
        ...

    def __mul__(self, other):
        ...

    def __rmul__(self, other):
        ...

    def __add__(self, other):
        ...

    def __sub__(self, other):
        ...

    def __eq__(self, other):
        ...

    def as_expr(self):
        ...

    @classmethod
    def get_embedding_angle_cos(self, i, j, k, l):
        ...

    @classmethod
    def get_positioning_cos(cls, i, j, k, l):
        ...

    def angle_cos(self, other, require_int=False):
        ...

    def dot(self, other: 'BondVector'):
        ...

    def cross(self, other):
        ...

    @classmethod
    def polar_components(self, embedding_atoms, vector):
        ...

class BondNormal:

    def __init__(self, a, b, norm=1, embedding=None):
        ...

    @property
    def direction(self):
        ...

    def __repr__(self):
        ...

    def __neg__(self):
        ...

    def __rmul__(self, other):
        ...

    def __mul__(self, other):
        ...

    def __add__(self, other):
        ...

    def __sub__(self, other):
        ...

    def __eq__(self, other):
        ...

    @classmethod
    def embedding_basis(cls, embedding):
        ...

    @classmethod
    def polar_components(cls, embedding, self):
        ...

    def angle_cos(self, other, require_int=True):
        ...

    def as_expr(self):
        ...

    def cross(self, other):
        ...

    def dot(self, other):
        ...

class BondVectorSum:

    def __init__(self, *terms):
        ...

    def __add__(self, other):
        ...

    def __radd__(self, other):
        ...

    def __rmul__(self, other):
        ...

    def __mul__(self, other):
        ...

    def __repr__(self):
        ...

    def __neg__(self):
        ...

    def simplify(self):
        ...

    def as_expr(self):
        ...

    def distribute(self, f):
        ...

    def dot(self, other):
        ...

class InternalJacobianDisplacements:

    @staticmethod
    def symbolic_m(i):
        ...

    @staticmethod
    def symbolic_r(i, j):
        ...

    @staticmethod
    def symbolic_a(i, j, k):
        ...

    @staticmethod
    def symbolic_t(i, j, k, l):
        ...

    @staticmethod
    def symbolic_y(i, j, k, l):
        ...

    @classmethod
    def lam(cls, i, j, k):
        ...

    @classmethod
    def dr(cls, i, j):
        ...

    @classmethod
    def da(cls, i, j, k):
        ...

    @classmethod
    def dt(cls, i, j, k, l):
        ...

    @classmethod
    def dy(cls, i, j, k, l):
        ...

    @classmethod
    def displacement_vectors(cls, inds, coord_type=None):
        ...