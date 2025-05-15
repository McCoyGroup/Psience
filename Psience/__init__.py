"""
The main package for working with core coordinate/molecule/wavefunction stuff.
Provides base classes and common/shared interfaces as well as a number of domain-specific applications.
The packages exposed include: `Molecools`, a handler for molecules (or really systems of atoms) with a
rich set of query-able properties; `DVR`, a package for arbitrary dimensional discrete variable representations
with the ability to generalize easily with full kinetic couplings; and `VPT2`, a misleadingly named package that
implements arbitrary order perturbation theory with flexible coordinate choice.
There are also a number of helper packages implementing things like a wavefunction interface, basis set representations,
and spectrum handling.
"""

__all__ = []

from . import Molecools
__all__ += ["Molecools"]
from . import DVR
__all__ += ["DVR"]
from . import VPT2
__all__ += ["VPT2"]
from . import DGB
__all__ += ["DGB"]
from . import Modes
__all__ += ["Modes"]
from . import Vibronic
__all__ += ["Vibronic"]
from . import Spectra
__all__ += ["Spectra"]
from . import Data
__all__ += ["Data"]
from . import AnalyticModels
__all__ += ["AnalyticModels"]