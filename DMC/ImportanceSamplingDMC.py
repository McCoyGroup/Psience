from .AbstractDMC import AbstractDMC
import numpy as np

__all__ = ["ImportanceSamplingDMC"]

class ImportanceSamplingDMC(AbstractDMC):

    def __init__(
            self,
            dmc,
            guide_function,
            derivative_function = None,
            mesh_spacing = .001,
            difference_stencil = 3
    ):
        """
        A class that wraps an existing DMC object and layers importance sampling onto it

        :param dmc: Existing DMC object to extend
        :type dmc: AbstractDMC
        :param guide_function: The guiding function for the importance sampling
        :type guide_function: function
        :param derivative_function: The function to take first, second, and mixed derivatives of the `guide_function`
        :type derivative_function: function
        :param mesh_spacing: the mesh spacing to feed in if `derivative_function` is `None`
        :type mesh_spacing: float
        :param difference_stencil: The number of points used in the stencil if `derivative_function` is `None`
        :type difference_stencil: int
        """

        if not isinstance(dmc, AbstractDMC):
            raise TypeError("{cls} extends an existing DMC object, {dmc} is not of a subtype of {pcls}".format(
                cls=type(self).__name__,
                pcls=AbstractDMC.__name__,
                dmc=dmc
            ))

        self.guide_function = guide_function
        self._custom_derivative_function = derivative_function is not None
        if not self._custom_derivative_function:
            from McUtils.Zachary import FiniteDifferenceFunction
            derivative_function = [
                # first derivatives
                FiniteDifferenceFunction.RegularGridFunction(1, mesh_spacing = mesh_spacing, stencil = 3, only_center = True).function,
                # second derivatives
                FiniteDifferenceFunction.RegularGridFunction(2, mesh_spacing = mesh_spacing, stencil = 3, only_center=True).function,
                # mixed derivatives
                FiniteDifferenceFunction.RegularGridFunction([1, 1], mesh_spacing = mesh_spacing, stencil=[3, 3], only_center=True).function
                ]
        self._mesh_spacing = mesh_spacing
        self.derivative_function = derivative_function

        # General optimizations
        self._current_psi_values = None # we'll cache these so we can make use of them inside `drift` and stuff
        self._cached_psi_arrays = None # used to not have to reallocate a block of memory for getting the psi values for the FD

    def displace_walkers(self):
        ...


    def _get_derivs(self, coords):
        """Default method for getting the derivatives using the FiniteDifferenceFunction we have bound

        :param coords:
        :type coords:
        :return:
        :rtype:
        """
        n_walkers = coords.shape[0]
        single_coord_shape = coords.shape[1:]
        n_coords = np.product(single_coord_shape)
        deriv_arrays = [
            np.empty(n_coords, n_walkers),
            np.empty(n_coords, n_walkers),
            np.empty(int(n_coords * (n_coords - 1) / 2), n_walkers)  # we only fill the upper triangle
        ]
        if self._cached_psi_arrays is None:
            # we cache some arrays so we never have to allocate them again -- shouldn't be huge memory hogs but in a memory constrained
            # environment this could be relaxed
            fd_psi_arrays = np.broadcast_to(self._current_psi_values, (3, 3) + self._current_psi_values.shape).copy()
            single_disp_psi_arrays = np.broadcast_to(self._current_psi_values, (n_coords, 2) + self._current_psi_values.shape).copy()
            self._cached_psi_arrays = (fd_psi_arrays, single_disp_psi_arrays)
        else:
            fd_psi_arrays = self._cached_psi_arrays[0]
            fd_psi_arrays[1, 1] = self._current_psi_values
            single_disp_psi_arrays = self._cached_psi_arrays[1]

        # we now have to do all of our displacements for each coordinate pair
        # this means a lot of += h and -= h but it's better to do that than to do a ton of copies
        # on the other hand we'll do a copy for the left-ward displacements and one for the right-ward displacements, just for convenience
        h = self._mesh_spacing
        coords_minus = coords.copy() # do one copy here so that we can cut down on two subtractions per inner loop
        coords_plus = coords.copy() # the first of these populates the 0 column, the second populates the 2 column

        # we'll loop through our pairs of indices in a snaking fashion so that we only ever change one of (i, j) at once
        for i in range(n_coords):
            idx = np.unravel_index(i, single_coord_shape)

            # handle the filling of the psi array if necessary
            # we'll also deal with the requisite displacements if needed
            odd_i = i % 2 == 1
            if i == 0:
                # on the first round through we gotta compute `psi` but every other time through `psi` will already have been calculated
                coords_minus[..., *idx] -= h # displace backwards in the i_th coordinate
                coords_plus[..., *idx] += h # displace forwards in the i_th coordinate
                fd_psi_arrays[1, 0] = single_disp_psi_arrays[i, 0] = self.guide_function(coords_minus)
                fd_psi_arrays[1, 2] = single_disp_psi_arrays[i, 1] = self.guide_function(coords_plus)
            elif odd_i: # only on odd values of i do we need to rebind this and re-displace
                coords_minus[..., *idx] -= h
                coords_plus[..., *idx] += h
                fd_psi_arrays[1, 0] = single_disp_psi_arrays[i, 0]
                fd_psi_arrays[1, 2] = single_disp_psi_arrays[i, 1]
            else:
                # these values were computed in the last step of the snaking
                # but to make the shape of the thing consistent we do this transpose
                fd_psi_arrays = fd_psi_arrays.T # just to make everything look consistent
                # we don't need to do a displacement because we can reuse the displacement from the previous step

            # first derivs and second derivs
            deriv_arrays[0][i] = self.derivative_function[0](fd_psi_arrays[1])  # this a 3 X N so ought to be clean...
            deriv_arrays[1][i] = self.derivative_function[1](fd_psi_arrays[1])

            # we want to snake, so if i is odd, we go in reverse order
            start = n_coords - 1 if odd_i else i + 1
            end = i if odd_i else n_coords
            step = -1 if odd_i else 1
            for j in range(start, end, step):
                jdx = np.unravel_index(j, single_coord_shape)
                if i == 0:
                    # special case where we compute the core psi values -- and hence have to displace the main set of coords, too
                    coords[..., *jdx] -= h
                    fd_psi_arrays[0, 1] = single_disp_psi_arrays[j, 0] = self.guide_function(coords)
                    coords[..., *jdx] += 2*h
                    fd_psi_arrays[2, 1] = single_disp_psi_arrays[j, 1] = self.guide_function(coords)
                    coords[..., *jdx] -= h # gotta reset what we did

                coords_minus[..., *jdx] += h  # this gives us a mixed evaluation corresponding to (-h, +h)
                fd_psi_arrays[2, 0] = self.guide_function(coords_minus)
                coords_minus[..., *jdx] -= 2*h # (-h, -h) -- we do forward then backward because we can sometimes reuse the backward
                fd_psi_arrays[0, 0] = self.guide_function(coords_minus)

                coords_plus[..., *jdx] -= h  # this gives us a mixed evaluation corresponding to (+h, -h)
                fd_psi_arrays[0, 2] = self.guide_function(coords_minus)
                coords_plus[..., *jdx] += 2 * h  # (+h, +h) -- we do backwards then forward because we can sometimes reuse the forward
                fd_psi_arrays[2, 2] = self.guide_function(coords_minus)

                # now the mixed derivative array is full
            coords_minus[..., *idx] += h # reset them
            coords_plus[..., *idx] -= h # reset them

        return deriv_arrays


    def get_derivs(self, coords):
        """Gets the derivatives for the importance sampled calculation

        By default uses a FiniteDifferenceFunction to evaluate the derivatives.
        Currently it's configured to have a single mesh_spacing, but it will be the case some day that we'll want to have multiple mesh_spacings
        This will be particularly relevant if propagating in internal coordinates (not sure if that's a thing or not, but I think it is...)

        :param coords:
        :type coords: np.ndarray
        :return:
        :rtype:
        """
        self._current_psi_values = self.guide_function(coords)
        if self._custom_derivative_function and not isinstance(self.derivative_function, (tuple, list)):
            # this is the kind of thing that'd be used for analytic derivatives
            derivs = self.derivative_function(self.guide_function, coords)
        else:
            # this is the default case where we need to take numerical derivatives
            derivs = self._get_derivs(coords)
        return derivs


    def drift(self, coords):
        """Computes the importance sampling drift term at the current set of coordinates

        :param coords:
        :type coords:
        :return:
        :rtype:
        """

        psi = psi_t(rch, interp, imp_samp_type, coords, multicore=multicore, interp_exp=interp_exp)
        deriv =
        blah = (psi[:, 2] - psi[:, 0]) / dx / psi[:, 1]
        return blah, psi

    # The metropolis step based on those crazy Green's functions
    def metropolis(Fqx, Fqy, x, y, sigmaCH, psi1, psi2):
        psi_1 = psi1[:, 1, 0, 0]
        psi_2 = psi2[:, 1, 0, 0]
        psi_ratio = (psi_2 / psi_1) ** 2
        a = np.exp(1. / 2. * (Fqx + Fqy) * (sigmaCH ** 2 / 4. * (Fqx - Fqy) - (y - x)))
        a = np.prod(np.prod(a, axis=1), axis=1) * psi_ratio
        return a

    def get_energy(self, coords=None, atoms=None):
        # this needs an update here
        raise NotImplemented