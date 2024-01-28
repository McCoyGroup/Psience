"""
Developed to decoupled molecules and vibrations, intended to become the
replacement for the already developed `MolecularNormalModes` / `MolecularVibrations`
"""

__all__ = []
from .NormalModes import *; from .NormalModes import __all__ as exposed
__all__ += exposed