from .BaseCoordinateSystem import BaseCoordinateSystem

######################################################################################################
##
##                                   CartesianCoordinateSystem Class
##
######################################################################################################
class CartesianCoordinateSystem(BaseCoordinateSystem):
    """Represents Cartesian coordinates generally

    """
    name = "Cartesian"
    def __init__(self, dimension = None):
        super().__init__(self.name, dimension=dimension)

######################################################################################################
##
##                                   InternalCoordinateSystem Class
##
######################################################################################################
class InternalCoordinateSystem(BaseCoordinateSystem):
    """Represents Internal coordinates generally

    """

    name = "Internal"
    def __init__(self, dimension = None):
        super().__init__(self.name, dimension=dimension)

######################################################################################################
##
##                                   CartesianCoordinates3D Class
##
######################################################################################################
class CartesianCoordinateSystem3D(CartesianCoordinateSystem):
    """Represents Cartesian coordinates in 3D

    """
    name = "Cartesian3D"
    def __init__(self):
        super().__init__(dimension=3)
CartesianCoordinates3D = CartesianCoordinateSystem3D()

######################################################################################################
##
##                                   ZMatrixCoordinateSystem Class
##
######################################################################################################
class ZMatrixCoordinateSystem(InternalCoordinateSystem):
    """Represents Cartesian coordinates generally

    """
    name = "ZMatrix"
    def __init__(self):
        super().__init__(dimension=6)
ZMatrixCoordinates = ZMatrixCoordinateSystem()

######################################################################################################
##
##                                   SphericalCoordinateSystem Class
##
######################################################################################################
class SphericalCoordinateSystem(BaseCoordinateSystem):
    """Represents Cartesian coordinates generally

    """
    name = "SphericalCoordinates"
    def __init__(self):
        super().__init__(self.name, dimension=3)
SphericalCoordinates = SphericalCoordinateSystem()