"""
An implementation of vibrational perturbation theory (VPT) that uses sparse matrix methods to obtain
corrections.
Makes heavy use of the `BasisReps` package as well as `McUtils.Coordinerds` to obtain representations
of the corrections to the vibrational Hamiltonian.
Is technically not restricted to VPT in a harmonic basis, but no other forms of PT are likely to
be implemented in the near future.
For the purposes of papers, we've been calling this implementation `PyVibPTn`

The code flow is detailed below

![pt design](/Psience/img/PyVibPTnDesign.png){:width="100%"}
"""

__all__ = []
from .Runner import *; from .Runner import __all__ as exposed
__all__ += exposed
from .Analyzer import *; from .Analyzer import __all__ as exposed
__all__ += exposed
from .Hamiltonian import *; from .Hamiltonian import __all__ as exposed
__all__ += exposed
from .Solver import *; from .Solver import __all__ as exposed
__all__ += exposed
from .Corrections import *; from .Corrections import __all__ as exposed
__all__ += exposed
from .Wavefunctions import *; from .Wavefunctions import __all__ as exposed
__all__ += exposed
from .Terms import *; from .Terms import __all__ as exposed
__all__ += exposed
from .StateFilters import *; from .StateFilters import __all__ as exposed
__all__ += exposed