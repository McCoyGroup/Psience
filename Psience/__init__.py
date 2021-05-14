"""
The main package for working with core coordinate/molecule/wavefunction stuff
Attempts to provide base classes and common/shared interfaces
"""

import Psience.Molecools as Molecools
import Psience.Wavefun as Wavefun
import Psience.DVR as DVR
import Psience.VPT2 as VPT2
import Psience.Data as Data
import Psience.BasisReps as BasisReps

# getting the full list of symbols explicitly in an __all__ variable
__all__ = [
    "Molecools",
    "Wavefun",
    "DVR",
    "VPT2",
    "Data",
    "BasisReps"
]