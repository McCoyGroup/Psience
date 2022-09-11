"""
Provides core Data-related types and data structures.
Intended to be high-level data as opposed to the lower-level stuff in `McUtils.Data`.
That means including stuff like dipole and potential energy surfaces that know how to compute their own properties.
We also have expressions for G-matrix elements from Frederick and Woywood to use with `sympy`.
"""


__all__ = []
from .Surfaces import *; from .Surfaces import __all__ as exposed
__all__ += exposed
from .KEData import *; from .KEData import __all__ as exposed
__all__ += exposed