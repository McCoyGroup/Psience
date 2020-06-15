"""
A DMC framework that aims to be flexible, general, and as efficient as possible within those constraints
Provides an abstract DMC class and some concrete implementations of the class
"""

#TODO: Large scale refactor to use the Simulation architecture of RynLib _without_ support for a PotentialCaller
# or MPI or any of that --> Force people to parallelize on their own, if they want it
# But _do_ support the importance sampling/normal DMC stuff

from .AbstractDMC import *
from .DiscreteWeightingDMC import *
from .ContinuousWeightingDMC import *
from .ImportanceSamplingDMC import *
from .PotentialEvaluator import *
from .WalkerSet import *
from .Wavefunction import *

from .AbstractDMC import __all__ as AbstractDMC__all__
from .DiscreteWeightingDMC import __all__ as DiscreteWeightingDMC__all__
from .ContinuousWeightingDMC import __all__ as ContinuousWeightingDMC__all__
from .ImportanceSamplingDMC import __all__ as ImportanceSamplingDMC__all__
from .PotentialEvaluator import __all__ as PotentialEvaluator__all__
from .WalkerSet import __all__ as WalkerSet__all__
from .Wavefunction import __all__ as Wavefunction__all__

__all__ = (
        AbstractDMC__all__ +
        DiscreteWeightingDMC__all__ +
        ContinuousWeightingDMC__all__ +
        ImportanceSamplingDMC__all__ +
        PotentialEvaluator__all__ +
        WalkerSet__all__ +
        Wavefunction__all__
)