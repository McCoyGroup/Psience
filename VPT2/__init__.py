"""
provides a class for doing (2nd order) vibrational perturbation theory in python
builds off of resource packages to handle most of the dirty work and just does the actual potential expansions
and pertubation theory computations
"""

from .PerturbationTheory import *

from .PerturbationTheory import __all__ as PT__all__
__all__ = PT__all__