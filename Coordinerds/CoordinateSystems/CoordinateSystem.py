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
        :rtype: CoordinateSystemConverter
        """

        return converters.get_converter(self, system)

    def convert_coords(self, coords, system, **kw):
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
        if self.matrix is None:
            return amts
        else:
            if isinstance(amts, (float, int, np.integer, np.float)):
                amts = np.full(self.matrix.shape[-1], amts)
            return np.matmul(self.matrix, amts)

    jacobian_prep_coordinates = None # used unsurprisingly in prepping data fed into Jacobian calculation

    def jacobian(self, coords, system, order = 1, mesh_spacing = .01, prep = None, coordinates = None, **fd_options):
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
        :type coordinates: None | list[list[int] | None]
        :param fd_options: options to be passed straight through to FiniteDifferenceFunction
        :type fd_options:
        :return:
        :rtype:
        """
        # TODO: make it possible to iterate over elements of the jacobian -- important for high dimensional cases

        from McUtils.Zachary import FiniteDifferenceFunction

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

        # from Peeves import Timer

        shape = coords.shape
        coord_shape = shape[-2:]

        # with Timer("setup code"):
        displacement = self.displacement(mesh_spacing) # provides a default displacement to work with
        if isinstance(displacement, float):
            displacement = np.full(coord_shape, displacement)
        elif displacement.shape == coord_shape[-1:]:
            displacement = np.broadcast_to(displacement, coord_shape)

        if coordinates is not None and coordinates[0] is not None:
            to_gen = coordinates[0]
        else:
            to_gen = np.arange(np.product(np.array(coord_shape)))
        if isinstance(to_gen[0], (int, np.integer)):
            to_gen = np.array(np.meshgrid(*([to_gen] * order)))
            to_gen = np.array(np.unravel_index(to_gen, coord_shape)).transpose()
            to_gen = to_gen.reshape((np.product(to_gen.shape[:-2]), ) + to_gen.shape[-2:] )

        if coordinates is not None and coordinates[1] is not None:
            gen_partners = tuple(coordinates[1])
        else:
            gen_partners = None

        if prep is None:
            try:
                prep = system.jacobian_prep_coordinates
            except AttributeError:
                pass
            if prep is None:
                prep = lambda c, a, b: (a, b)

        fdf = FiniteDifferenceFunction.RegularGridFunction([ 1 ]*order,
                                                           end_point_precision = 0,
                                                           only_core = True,
                                                           **fd_options
                                                           )

        stencil_widths = tuple( len(cf[1]) for cf in fdf.coefficients )
        stencil_shapes = tuple( w[1] for w in fdf.widths )
        finite_difference = fdf.get_FDF(shape = stencil_widths)

        multiconfig = cu.is_multiconfig(coords)
        if not multiconfig:
            coords = np.array([coords])

        full_jacob = None

        for coord_num, coord in enumerate(to_gen):
            # not too worried about looping over coordinates since the number of loops will be like max in the small hundreds

            # with Timer("create displacements ({})".format(order)):
            num_displacements = np.product(stencil_widths)
            displacement_shape = (num_displacements, ) + coords.shape[1:]
            displacements = np.zeros(displacement_shape)
            base_roll = tuple(np.arange(len(stencil_widths)))
            for stencil_shape, c in zip(stencil_shapes, coord): # coord can be like ((1, 1), (1, 2)) in which case we have a 2D derivative
                # creates single displacement matrix
                ctup = tuple(c)
                disp = displacement[ctup]
                steps = np.arange(-stencil_shape[0], stencil_shape[1])
                full_disp = disp * steps
                to_set = np.broadcast_to(full_disp, stencil_widths).transpose(np.roll(base_roll, -2)).flatten()
                # in dimension 1 we want to have this repeat over the slowest moving indices, then the ones after those
                # then after those, etc.
                idx = (...,) + ctup # I'm imposing some assumptions on the form of this for the moment, but we can relax those later...
                displacements[idx] = to_set

            # then we broadcast *this* up to the total number of walkers we have
            # with Timer("create displaced coords ({})".format(order)):
            full_target_shape = (len(coords), ) + displacement_shape
            coords_expanded = np.expand_dims(coords, 1)
            displaced_coords = np.broadcast_to(coords_expanded, full_target_shape) # creates the displaced coordinate sets
            full_full_displacement = np.broadcast_to(displacements, full_target_shape)
            displaced_coords = displaced_coords + full_full_displacement

            # with Timer("convert displaced coords ({})".format(order)):
            # finally we can convert all of our walkers in batch
            conv = self.convert_coords(displaced_coords, system)

            # now we loop over each of these walkers in turn... this could be slow so I probably need a way to
            # optimize this to do just like a single call...
            # with Timer("finite difference stuff ({})".format(order)):
            # for config_num, converted_coords in enumerate(conv):

            disp, converted_coords = prep(coord, displacements, conv)

            # converted_coords is for _all_ of the walkers so we need to reshape it
            num_coords = np.product(converted_coords.shape[-2:])
            new_coord_shape = converted_coords.shape[:-2] + (num_coords, )
            # raise Exception(new_coord_shape)
            converted_coords = converted_coords.reshape(new_coord_shape)

            # now we potentially filter out some coords...
            if gen_partners is not None:
                converted_coords = converted_coords[:, gen_partners]

            # finally we restructure this into a tensor of the appropriate dimension for feeding into the FD code
            converted_coords = converted_coords.reshape((len(conv), ) + stencil_widths +  converted_coords.shape[-1:])

            def _get_diff(c):
                diffs = np.abs(np.diff(disp[(..., ) + tuple(c) ]))
                return np.sort(diffs)[-1] # get the only diff > 0 (assumes we have a regular grid)

            h = [ _get_diff(c) for c in coord ]
            derivs = finite_difference(converted_coords, h = h, axis = 1)

            deriv_center = tuple( int(np.floor(s/2)) for s in derivs.shape[:order] )
            derivs = derivs[deriv_center]

            if full_jacob is None: # we finally know what size it'll be...
                jacob_shape = tuple( len( np.unique( to_gen[:, i], axis=0 ) ) for i in range(order) )
                full_jacob = np.zeros((len(coords), ) + jacob_shape + derivs.shape[1:])

            for config_num in range(len(coords)): # a python loop but hopefully a relatively cheap one...
                set_index = (config_num, ) + tuple( np.ravel_multi_index(c, coord_shape) for c in coord )
                full_jacob[set_index] = derivs[config_num]

        if not multiconfig:
            full_jacob = full_jacob[0]

        return full_jacob

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