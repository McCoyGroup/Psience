from .CoordinateSystem import CoordinateSystem

######################################################################################################
##
##                                   BaseCoordinateSystem Class
##
######################################################################################################

class BaseCoordinateSystem(CoordinateSystem):
    """A CoordinateSystem object that can't be reduced further.
    A common choice might be Cartesian coordinates or internal coordinates

    """

    def __init__(self, name, dimension = None, coeffs = None):
        super().__init__(name=name, dimension=dimension, basis=self, )

    @property
    def basis(self):
        return self._basis
    @property
    def matrix(self):
        return self._matrix
    @property
    def dimension(self):
        return self._dim
