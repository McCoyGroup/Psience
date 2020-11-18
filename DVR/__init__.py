'''
A package for doing generalized DVR in python.
Provides and extensible DVR framework with an easy-to-write structure.
'''

#TODO: migrate Mathematica structure -> more DVR classes, allow for a general direct-product DVR
#       provide hook-in for non-direct-product couplings, allow for coordinate-dependent mass

__all__= [ "DVR" ]
from .DVR import *
