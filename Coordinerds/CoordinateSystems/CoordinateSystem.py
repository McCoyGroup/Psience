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
    def __init__(self):
        pass

    @property
    @abc.abstractmethod
    def basis(self):
        """The basis for the representation of CoordinateSystem.matrix

        :return:
        :rtype: CoordinateSystem
        """
        pass

    @property
    @abc.abstractmethod
    def matrix(self):
        """The matrix representation in the CoordinateSystem.basis
        None is shorthand for the identity matrix

        :return:
        :rtype:
        """
        pass

    @property
    @abc.abstractmethod
    def dimension(self):
        """The dimension of the coordinate system
        None means unspecified dimension

        :return:
        :rtype: int or None
        """
        pass


