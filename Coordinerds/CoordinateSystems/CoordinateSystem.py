import numpy as np
from .CoordinateSystemConverter import CoordinateSystemConverters as converters, CoordinateSystemConverter
from .CoordinateUtils import is_multiconfig, mc_safe_apply

######################################################################################################
##
##                                   CoordinateSystem Class
##
######################################################################################################

class CoordinateSystem:
    """A representation of a coordinate system. It doesn't do much on its own but it *does* provide a way
    to unify internal, cartesian, derived type coordinates

    """
    def __init__(self,
                 name = None, basis = None, matrix = None, dimension = None,
                 jacobian_prep = None, coordinate_shape = None
                 ):
        """Sets up the CoordinateSystem object

        :param name: a name to give to the coordinate system
        :type name: str
        :param basis: a basis for the coordinate system
        :type basis:
        :param matrix: an expansion coefficient matrix for the set of coordinates in its basis
        :type matrix: np.ndarray | None
        :param dimension: the dimension of a single coordinate in the coordinate system (for validation)
        :type dimension:
        :param jacobian_prep: a function for preparing coordinates to be used in computing the Jacobian
        :type jacobian_prep: function | None
        :param coordinate_shape: the actual shape of a single coordinate in the coordinate system
        :type coordinate_shape: iterable[int]
        """
        self.name = name
        self._basis = basis
        self._matrix = matrix
        self._dimension = matrix.shape[-1] if matrix is not None else dimension
        self.jacobian_prep = jacobian_prep
        self.coordinate_shape = coordinate_shape
        self._validate()

    def _validate(self):
        if self._matrix is None:
            pass
        elif len(self._matrix.shape) != 2:
            raise CoordinateSystemException("{}: expansion matrix must be a matrix".format(type(self).__name__))
        # elif self._matrix.shape[0] != self._matrix.shape[1]:
        #     raise CoordinateSystemException("{}: expansion matrix must square".format(type(self).__name__))

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
        :rtype: CoordinateSystemConverter
        """

        return converters.get_converter(self, system)

    def convert_coords(self, coords, system, **kw):
        if self.matrix is not None:
            coords = np.tensordot(self.matrix, coords, axes = ((1,), (1,))) # set up the basic expansion coordinates in the more primitive system
            coords = coords.T
            ndim = self.basis.dimension
            coords = np.reshape(coords,
                                coords.shape[:1] + (coords.shape[1] // ndim, ndim) # this currently only handles stuff like 9 -> (3, 3)...
                                )# reshape so as to actually fit the dimension of the basis

            if system is not self.basis:
                coords = self.basis.convert_coords(coords, system, **kw)
            return coords

        else:
            converter = self.converter(system)
            if is_multiconfig(coords):
                fun = lambda coords: converter.convert_many(coords, **kw)
            else:
                fun = lambda coords: converter.convert(coords, **kw)
            new_coords = mc_safe_apply(fun, coords=coords)
            return new_coords

    def displacement(self, amts):
        """Generates a displacement or matrix of displacements based on the vector or matrix amts

        :param amts:
        :type amts: np.ndarray
        :return:
        :rtype: np.ndarray
        """
        return amts # I used to think this would be relevant... but I don't really know now if it is
        # if self.matrix is None:
        #     return amts
        # else:
        #     if isinstance(amts, (float, int, np.integer, np.float)):
        #         amts = np.full((1, self.matrix.shape[-1]), amts)
        #     else:
        #         amts = np.asarray(amts)
        #         if amts.ndim == 1:
        #             amts = np.broadcast_to(np.reshape(amts, amts.shape + (1,)), amts.shape + self.matrix.shape[-1:])
        #     return np.matmul(amts, self.matrix)

    def jacobian(self, coords, system,
                 order = 1,
                 mesh_spacing = .01,
                 prep = None,
                 coordinates = None,
                 return_array = True,
                 **fd_options
                 ):
        """Computes the Jacobian between the current coordinate system and a target coordinate system

        :param system: the target CoordinateSystem
        :type system: CoordinateSystem
        :param order: the order of the Jacobian to compute, 1 for a standard, 2 for the Hessian, etc.
        :type order: int
        :param mesh_spacing: the spacing to use when displacing
        :type mesh_spacing: float | np.ndarray
        :param prep: a function for pre-validating the generated coordinate values and grids
        :type prep: None | function
        :param coordinates: a spec of which coordinates to generate derivatives for (None means all)
        :type coordinates: None | iterable[iterable[int] | None]
        :param fd_options: options to be passed straight through to FiniteDifferenceFunction
        :type fd_options:
        :return:
        :rtype:
        """
        # TODO: make it possible to iterate over elements of the jacobian -- important for high dimensional cases

        from McUtils.Zachary import FiniteDifferenceDerivative

        # as a basic example, let's think about the Jacobian between Cartesians and ZMatrix coordinates
        # we need to take the derivative of every internal coordinate with respect to Cartesian (or vice versa)
        # this means we need to define some set of pertubations for our Cartesians (step them all forward / backward some amount basically)
        # and then we need to compute the new ZMatrix coordinates for each perturbation
        # after that we need to use the mesh spacing for the stepping along each coordinate to compute the FD for that coordinate
        # in general, then, the procedure looks like:
        #
        # To get dq/dx1:
        #     [
        #         [
        #             [x1 - xh, y1, z1]
        #             [x2,      y2, z2]
        #             [x3,      y3, z3]
        #         ]
        #         [
        #             [x1,     y1, z1]
        #             [x2,     y2, z2]
        #             [x3,     y3, z3]
        #         ]
        #         [
        #             [x1 + xh, y1, z1]
        #             [x2,      y2, z2]
        #             [x3,      y3, z3]
        #         ]
        #     ]
        #
        #     ->
        #
        #     [
        #         [
        #             [r11,      0,     0]
        #             [r21,      a21,   0]
        #             [r31,      a31, d31]
        #         ] # for the -xh
        #         [
        #             [r12,      0,     0]
        #             [r22,      a22,   0]
        #             [r32,      a32, d32]
        #         ]
        #         [
        #             [r13,      0,     0]
        #             [r23,      a23,   0]
        #             [r33,      a33, d33]
        #         ]  # for the +xh
        #     ]
        #
        # Then I need to take [ FD([x1-xh, x1, x1+xh], [r11, r12, r13]),  FD([x1-xh, x1, x1+xh], [r21, r22, r23]), ...]
        #
        # Then repeat for y1, z1, x2, y2, z2, ...
        #
        #   In the multidimensional case I need to generate the whole matrix of displaced coordinates. Like I need to take:
        #
        #       Map[
        #           convert@*displace,
        #           {
        #               {x-h, y-h}, {x, y-h}, {x+h, y-h},
        #               {x-h, y},   {x, y},   {x+h, y},
        #               {x-h, y+h}, {x, y+h}, {x+h, y+h},
        #           },
        #           {2}
        #           ]

        convert = self.converter(system)
        self_shape = self.coordinate_shape
        if self_shape is None:
            raise CoordinateSystemException(
                "{}.{}: 'coordinate_shape' {} must be tuple of ints".format(
                    type(self).__name__,
                    'jacobian',
                    self_shape
            ))

        other_shape = system.coordinate_shape
        if other_shape is None:
            raise CoordinateSystemException(
                "{}.{}: 'coordinate_shape' {} must be tuple of ints".format(
                    type(self).__name__,
                    'jacobian',
                    other_shape
                ))

        if prep is None:
            prep = system.jacobian_prep
        deriv = FiniteDifferenceDerivative(
            convert,
            order,
            mesh_spacing = mesh_spacing,
            prep = prep,
            coordinates = coordinates,
            function_shape = (self_shape, other_shape),
            **fd_options
        ).derivatives(coords)

        if return_array:
            return deriv.array
        else:
            return deriv

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