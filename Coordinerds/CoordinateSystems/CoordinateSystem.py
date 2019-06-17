import abc

######################################################################################################
##
##                                   CoordinateSystem Class
##
######################################################################################################

class CoordinateSystem(metaclass=abc.ABCMeta):
    """A representation of a coordinate system. It doesn't do much on its own but it *does* provide a way
    to unify internal, cartesian, derived type coordinates

    """
    @abc.abstractmethod
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
    @abc.abstractmethod
    def basis(self):
        """The basis for the representation of CoordinateSystem.matrix

        :return:
        :rtype: CoordinateSystem
        """
        return self._basis

    @property
    @abc.abstractmethod
    def matrix(self):
        """The matrix representation in the CoordinateSystem.basis
        None is shorthand for the identity matrix

        :return:
        :rtype:
        """
        return self._matrix

    @property
    @abc.abstractmethod
    def dimension(self):
        """The dimension of the coordinate system
        None means unspecified dimension

        :return:
        :rtype: int or None
        """
        return self._dimension

######################################################################################################
##
##                                   CoordinateSystemException Class
##
######################################################################################################
class CoordinateSystemException(Exception):
    pass
