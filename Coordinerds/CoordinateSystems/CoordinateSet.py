import numpy as np
from .CoordinateSystemConverter import CoordinateSystemConverters as converters
from .CommonCoordinateSystems import CartesianCoordinates3D

######################################################################################################
##
##                                   CoordinateSet Class
##
######################################################################################################


class CoordinateSet(np.ndarray):
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
        return "{}({})".format(type(self.system).__name__, super().__str__())

    @property
    def multiconfig(self):
        """Determines whether self.coords represents multiple configurations of the coordinates

        :return:
        :rtype:
        """
        return len(self.shape) > 2

    def _mc_safe_apply(self, fun):
        """Applies fun to the coords in such a way that it will apply to an array of valid
        coordinates (as determined by dimension of the basis). This will either be a single configuration
        or multiple configurations

        :param fun:
        :type fun:
        :return:
        :rtype:
        """
        coords = self
        if self.multiconfig:
            base_shape = coords.shape
            new_shape = (np.product(base_shape[:-2]),) + base_shape[-2:]
            coords = np.reshape(coords, new_shape)
            new_coords = fun(coords)
            revert_shape = tuple(base_shape[:-2]) + new_coords.shape[1:]
            new_coords = np.reshape(new_coords, revert_shape)
        else:
            new_coords = fun(coords)
        return new_coords

    def transform(self, tf):
        """Applies a transformation to the stored coordinates

        :param system:
        :type system:
        :return:
        :rtype:
        """
        coords = self
        new_coords = tf(coords)
        return type(self)(new_coords, self.system)

    def convert(self, system, **kw):
        """Converts across coordinate systems

        :param system:
        :type system:
        :return:
        :rtype:
        """
        cosys = self.system
        converter = converters.get_converter(cosys, system)
        if self.multiconfig:
            fun = lambda coords: converter.convert_many(coords, **kw)
        else:
            fun = lambda coords: converter.convert(coords, **kw)
        new_coords = self._mc_safe_apply(fun)
        return type(self)(new_coords, system)
