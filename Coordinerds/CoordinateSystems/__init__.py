"""
Provides a framework for using coordinates with explicit reference to an underlying coordinate system
"""

from .CoordinateSet import CoordinateSet
from .CommonCoordinateSystems import InternalCoordinateSystem, CartesianCoordinateSystem, ZMatrixCoordinateSystem, \
    SphericalCoordinateSystem
from .CommonCoordinateSystems import CartesianCoordinates3D, ZMatrixCoordinates, SphericalCoordinates
from .CoordinateSystem import *
from .CoordinateSystemConverter import *

CoordinateSystemConverters._preload_converters()