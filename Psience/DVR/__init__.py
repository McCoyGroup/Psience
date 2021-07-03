'''
A package for doing generalized DVR in python.
Provides and extensible DVR framework with an easy-to-write structure.
'''

#TODO: migrate Mathematica structure -> more DVR classes, allow for a general direct-product DVR
#       provide hook-in for non-direct-product couplings, allow for coordinate-dependent mass

__all__= [ ]
from .BaseDVR import *; from .BaseDVR import __all__ as exposed
__all__ += exposed
from .ColbertMiller import *; from .ColbertMiller import __all__ as exposed
__all__ += exposed
from .DirectProduct import *; from .DirectProduct import __all__ as exposed
__all__ += exposed
from .Wavefunctions import *; from .Wavefunctions import __all__ as exposed
__all__ += exposed
