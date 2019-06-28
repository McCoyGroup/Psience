import numpy as np
from .CoordinateSystemConverter import CoordinateSystemConverters as converters

######################################################################################################
##
##                                   CoordinateSystem Class
##
######################################################################################################

class CoordinateSystem:
    """A representation of a coordinate system. It doesn't do much on its own but it *does* provide a way
    to unify internal, cartesian, derived type coordinates

    """
    def __init__(self, name = None, basis = None, matrix = None, dimension = None):
        self.name = name
        self._basis = basis
        self._matrix = matrix
        self._dimension = dimension
        self._validate()

    def _validate(self):
        if self._matrix is None:
            return True
        elif len(self._matrix.shape) != 2:
            raise CoordinateSystemException("{}: expansion matrix must be a matrix".format(type(self).__name__))
        elif self._matrix.shape[0] != self._matrix.shape[1]:
            raise CoordinateSystemException("{}: expansion matrix must square".format(type(self).__name__))
        else:
            return True

    @property
    def basis(self):
        """The basis for the representation of CoordinateSystem.matrix

        :return:
        :rtype: CoordinateSystem
        """
        return self._basis

    @property
    def matrix(self):
        """The matrix representation in the CoordinateSystem.basis
        None is shorthand for the identity matrix

        :return:
        :rtype:
        """
        return self._matrix

    @property
    def dimension(self):
        """The dimension of the coordinate system
        None means unspecified dimension

        :return:
        :rtype: int or None
        """
        return self._dimension

    def converter(self, system):
        """Gets the converter from the current system to a new system

        :param system: the target CoordinateSystem
        :type system: CoordinateSystem
        :return:
        :rtype: CoordinateSystem
        """

        return converters.get_converter(self, system)

    def displacement(self, amts):
        """Generates a displacement or matrix of displacements based on the vector or matrix amts

        :param amts:
        :type amts: np.ndarray
        :return:
        :rtype: np.ndarray
        """
        if self.matrix is None:
            return amts
        else:
            if isinstance(amts, (float, int, np.integer, np.float)):
                amts = np.full(self.matrix.shape[-1], amts)
            return np.matmul(self.matrix, amts)

    jacobian_prep_coordinates = None # used unsurprisingly in prepping data fed into Jacobian calculation

######################################################################################################
##
##                                   CoordinateSystemException Class
##
######################################################################################################
class CoordinateSystemException(Exception):
    pass


######################################################################################################
##
##                                   BaseCoordinateSystem Class
##
######################################################################################################

class BaseCoordinateSystem(CoordinateSystem):
    """A CoordinateSystem object that can't be reduced further.
    A common choice might be Cartesian coordinates or internal coordinates

    """

    def __init__(self, name, dimension = None, matrix = None):
        super().__init__(name=name, dimension=dimension, basis=self, matrix=matrix)