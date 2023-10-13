'''
A package for managing local-mode models.
In the future it should have Sibert-style CH stretch models & coupled Morse oscillator models
'''

__all__= [ ]
from .BlockLocalFG import *; from .BlockLocalFG import __all__ as exposed
__all__ += exposed
