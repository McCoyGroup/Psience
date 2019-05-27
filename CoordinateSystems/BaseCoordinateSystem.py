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

    def __init__(self, name, dimension = None):
        self.name = name
        self.__basis = self
        self.__matrix = None
        self.__dim = dimension
        super().__init__()

    @property
    def basis(self):
        return self.__basis
    @property
    def matrix(self):
        return self.__matrix
    @property
    def dimension(self):
        return self.__dim
