"""
Provides a CoordinateSet class that acts as a symbolic extension of np.ndarray to provide an explicit basis
"""
import numpy as np
from .CoordinateSystem import CoordinateSystem
from .CommonCoordinateSystems import CartesianCoordinates3D
from .CoordinateUtils import is_multiconfig

######################################################################################################
##
##                                   CoordinateSet Class
##
######################################################################################################


class CoordinateSet(np.ndarray):
    """
    A subclass of np.ndarray that lives in an explicit coordinate system and can convert between them
    """
    # note that the subclassing procedure here is very non-standard because of how ndarray works

    def __new__(cls, coords, system = CartesianCoordinates3D):
        self = np.asarray(coords).view(cls)
        self.system = system
        return self

    def __init__(self, *args, **kwargs):
        # all the heavy lifting is done in _validate
        # and we only want to ever call this once
        self._validate()

    def __array_finalize__(self, coords):
        # basically just a validator...
        if coords is None:
            return None

        self.system = getattr(coords, "system", CartesianCoordinates3D)

    def _validate(self):
        base_dim = self.system.dimension
        if isinstance(base_dim, int):
            core_dim = self.shape[-1]
        else:
            base_dim = tuple(base_dim)
            core_dim = tuple(self.shape[-len(base_dim):])
        if base_dim != core_dim:
            raise TypeError(
                "Dimension of basis {} ({}) and dimension of coordinate set ({}) misaligned".format(
                    self.system,
                    self.system.dimension,
                    core_dim
                )
            )

    def __str__(self):
        return "{}({}, {})".format(type(self).__name__, self.system.name, super().__str__())

    @property
    def multiconfig(self):
        """Determines whether self.coords represents multiple configurations of the coordinates

        :return:
        :rtype:
        """
        return is_multiconfig(self)

    def transform(self, tf):
        """Applies a transformation to the stored coordinates

        :param tf: the transformation function to apply to the system
        :type tf:
        :return:
        :rtype:
        """
        coords = self
        new_coords = tf(coords)
        return type(self)(new_coords, system = self.system)

    def convert(self, system, **kw):
        """Converts across coordinate systems

        :param system: the target coordinate system
        :type system: CoordinateSystem
        :return: new_coords
        :rtype: CoordinateSet
        """
        new_coords = self.system.convert_coords(self, system, **kw)
        return type(self)(new_coords, system = system)

    def jacobian(self, system, order = 1, mesh_spacing = .01, prep = None, coordinates = None, **fd_options):
        """Delegates to the jacobian function of the current coordinate system

        :param system:
        :type system:
        :param order:
        :type order:
        :param mesh_spacing:
        :type mesh_spacing:
        :param prep:
        :type prep:
        :param coordinates:
        :type coordinates:
        :param fd_options:
        :type fd_options:
        :return:
        :rtype:
        """

        return self.system.jacobian(self, system, order=order, mesh_spacing=mesh_spacing, prep=prep, coordinates=coordinates,
                                    **fd_options
                                    )
