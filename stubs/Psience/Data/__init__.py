"""
Provides core Data-related types and data structures.
Intended to be high-level data as opposed to the lower-level stuff in `McUtils.Data`.
That means including stuff like dipole and potential energy surfaces that know how to compute their own properties.
We also have expressions for G-matrix elements from Frederick and Woywood to use with `sympy`.
"""
__all__ = ['DipoleSurface', 'PotentialSurface', 'KEData', 'KEDataHandler', 'PotentialRegistryAPI', 'ScanManager', 'shape_scan_iterator', 'scan_iterator', 'molecule_scan_geometries_iterator', 'molecule_displaced_geometries_iterator', 'molecule_atom_position_scan_iterator']
from .Surfaces import *
from .KEData import *
from .PotentialRegistry import *
from .ScanManager import *