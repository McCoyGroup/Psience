"""
Provides core Data-related types and data structures.
Intended to be high-level data as opposed to the lower-level stuff in `McUtils.Data`.
That means including stuff like dipole and potential energy surfaces that know how to compute their own properties.
Currently...well that's all we have. But wrappers for commonly-used potentials & bases could well come.
Not sure at this point, though.
"""


from .Surfaces import *

__all__ = []
__all__ += Surfaces.__all__