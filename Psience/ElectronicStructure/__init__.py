"""
Electronic-structure models and descriptors, kept separate from `Molecule`.

`ElectronicStructureBackend`s are the raw evaluation layer for MLIPs (AIMNet2, HIPNN, MACE,
MACE-POLAR, UMA, UPET, Orb), xTB, and PySCF: energies, charges, dipoles, and polarizabilities
with their analytic Cartesian derivatives, in each model's native units. The
`Psience.Molecools` evaluators are thin wrappers over them that add unit handling, finite
differences, embeddings, and optimization.

On top of the backends, `ElectronicStateModel`s predict atomic charge (and spin) densities for a
given charge and multiplicity, represented as atom-centered Gaussian-smeared multipole
expansions (`AtomCenteredDensity`, on `McUtils.Zachary.MultipoleExpansion`), from which spin
densities, Fukui functions, and related conceptual-DFT quantities are built.
"""

__all__ = []
from .Backends import *; from .Backends import __all__ as exposed
__all__ += exposed
from .ChargeResponse import *; from .ChargeResponse import __all__ as exposed
__all__ += exposed
from .MLIPBackends import *; from .MLIPBackends import __all__ as exposed
__all__ += exposed
from .QuantumChemistryBackends import *; from .QuantumChemistryBackends import __all__ as exposed
__all__ += exposed
from .Densities import *; from .Densities import __all__ as exposed
__all__ += exposed
from .StateModels import *; from .StateModels import __all__ as exposed
__all__ += exposed
from .ConceptualDFT import *; from .ConceptualDFT import __all__ as exposed
__all__ += exposed
