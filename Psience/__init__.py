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

import Psience.Molecools as Molecools
import Psience.Wavefun as Wavefun
import Psience.DVR as DVR
import Psience.VPT2 as VPT2
import Psience.Data as Data
import Psience.BasisReps as BasisReps
import Psience.AnalyticModels as AnalyticModels

# getting the full list of symbols explicitly in an __all__ variable
__all__ = [
    "Molecools",
    "Wavefun",
    "DVR",
    "VPT2",
    "Data",
    "BasisReps",
    "Spectra",
    "AnalyticModels"
]