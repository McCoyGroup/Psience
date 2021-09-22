"""
Provides implementations of the Colbert-Miller DVR types defined in
http://xbeams.chem.yale.edu/~batista/v572/ColbertMiller.pdf
"""


import numpy as np
from .BaseDVR import BaseDVR

__all__ = [
    "CartesianDVR",
    "RingDVR",
    "RadialDVR",
    "PolarDVR",
]

class CartesianDVR(BaseDVR):
    """
    Provides the Colbert Miller DVR on the Cartesian [-inf, inf] range
    """

    def get_grid(self, domain=None, divs=None, **kw):
        """
        Provides the Colbert-Miller DVR grid for the [-inf, inf] range
        :param domain:
        :type domain:
        :param divs:
        :type divs:
        :param flavor:
        :type flavor:
        :param kw:
        :type kw:
        :return:
        :rtype:
        """

        return domain[0] + (domain[1] - domain[0]) * np.arange(1, divs) / (divs + 1)

    def get_kinetic_energy(self, grid=None, mass=None, hb=1, **kwargs):

        dx = grid[1] - grid[0]  # recomputed here simply to decouple the calling from dvr_grid
        divs = len(grid)
        ke = np.empty((divs, divs))

        coeff = (hb ** 2) / (2 * mass * (dx ** 2))
        # compute the band values for the first row
        b_val_0 = coeff * (np.pi ** 2) / 3
        col_rng = np.arange(1,
                            divs + 1)  # the column indices -- also what will be used for computing the off diagonal bands
        row_rng = np.arange(0, divs)  # the row indices -- computed once and sliced
        b_vals = coeff * ((-1) ** col_rng) * 2 / (col_rng ** 2)

        for i in range(divs):
            if i == 0:
                np.fill_diagonal(ke, b_val_0)
            else:
                col_inds = col_rng[i - 1:-1]  # +(i-1)
                row_inds = row_rng[:-i]
                ke[row_inds, col_inds] = b_vals[i - 1]
                ke[col_inds, row_inds] = b_vals[i - 1]

        return ke

    def real_momentum(self, grid=None, mass=None, hb=1, **kwargs):
        """
        Provides the real part of the momentum for the [0, 2pi] range
        :param grid:
        :type grid:
        :param hb:
        :type hb:
        :param kw:
        :type kw:
        :return:
        :rtype:
        """

        divs = len(grid)
        dx = grid[1] - grid[0]
        p = np.zeros((divs, divs))

        col_rng = np.arange(1, divs + 1)
        row_rng = np.arange(0, divs)
        for i in range(1, divs):
            col_inds = col_rng[i - 1:-1]  # +(i-1)
            row_inds = row_rng[:-i]
            val = hb / dx * (-1) ** (i) / i
            p[row_inds, col_inds] = val
            p[col_inds, row_inds] = -val

        return p

class RingDVR(BaseDVR):
    """
    Provides a DVR for working on the (0, 2Pi) range with periodicity from Colbert and Miller
    """

    def __init__(self, domain=None, **opts):
        if domain is None:
            domain = (0, 2*np.pi)
        super().__init__(domain=domain, **opts)

    def get_grid(self, domain=None, divs=None, **kw):
        """
        Provides the Colbert-Miller 1D grid for the [0, 2Pi] range
        :param domain:
        :type domain:
        :param divs:
        :type divs:
        :param kw:
        :type kw:
        :return:
        :rtype:
        """

        if divs % 2 != 1:
            raise ValueError('number of DVR points must go as (2N + 1), i.e. it must be odd')

        return domain[0] + (domain[1] - domain[0]) * np.arange(1, divs + 1) / divs

    def get_kinetic_energy(self, grid=None, mass=1, hb=1, **kw):
        """
        Colbert-Miller kinetic energy for the [0, 2pi] range
        :param grid:
        :type grid:
        :param mass:
        :type mass:
        :param hb:
        :type hb:
        :param kw:
        :type kw:
        :return:
        :rtype:
        """

        coeff = hb ** 2 / (2 * mass)
        divs = len(grid)
        n = (divs - 1) // 2
        ke = np.empty((divs, divs))

        diag_val = coeff * n * (n + 1) / 3
        np.fill_diagonal(ke, diag_val)

        col_rng = np.arange(1, divs + 1)  # the column indices -- also what will be used for computing the off diagonal bands
        row_rng = np.arange(0, divs)  # the row indices -- computed once and sliced
        for i in range(1, divs):
            col_inds = col_rng[i - 1:-1]  # +(i-1)
            row_inds = row_rng[:-i]
            val = coeff * (-1) ** (i) * np.cos(i * np.pi / divs) / (2 * np.sin(i * np.pi / divs) ** 2)
            ke[row_inds, col_inds] = val
            ke[col_inds, row_inds] = val

        return ke

    def real_momentum(self, grid=None, hb=1, **kw):
        """
        Provides the real part of the momentum for the [0, 2pi] range
        :param grid:
        :type grid:
        :param hb:
        :type hb:
        :param kw:
        :type kw:
        :return:
        :rtype:
        """

        # divs = len(grid)
        # p = np.zeros((divs, divs))
        # for j in range(1, divs+1):
        #     for jp in range(1, divs+1):
        #         i = jp - j
        #         if i == 0:
        #             p[j-1, jp-1] = 0
        #         else:
        #             val = hb / 2 * ((-1) ** i) / np.sin(i * np.pi / divs)
        #             p[j-1, jp-1] = val
                    # p[jp-1, j-1] = val
        # print("FIRST", p)

        coeff = hb / 2
        divs = len(grid)
        p = np.zeros((divs, divs))
        col_rng = np.arange(1, divs + 1)  # the column indices -- also what will be used for computing the off diagonal bands
        row_rng = np.arange(0, divs)  # the row indices -- computed once and sliced
        for i in range(1, divs): # loop over possible values of j - j'
            col_inds = col_rng[i-1:-1]  # +(i-1)
            row_inds = row_rng[:-i]
            val = coeff * ( (-1)**i ) / np.sin(i * np.pi / divs)
            p[row_inds, col_inds] = -val
            p[col_inds, row_inds] =  val

        return p

class PolarDVR(BaseDVR):
    """
    Provides a DVR for working on the (0, pi) range from Colbert and Miller
    """

    def __init__(self, domain=None, **opts):
        if domain is None:
            domain = (0, np.pi)
        super().__init__(domain=domain, **opts)

    def get_grid(self, domain=(0, np.pi), divs=None, **kwargs):
        """
        Provides the grid appropriate for the Colbert-Miller (0, Pi) range

        :param domain:
        :type domain:
        :param divs:
        :type divs:
        :param kwargs:
        :type kwargs:
        :return:
        :rtype:
        """

        if divs % 2 != 1:
            raise ValueError('number of DVR points must go as (2N + 1), i.e. it must be odd')

        return domain[0] + (domain[1] - domain[0]) * np.arange(1, divs) / divs

    def get_kinetic_energy(self, grid=None, mass=None, hb=1, **kwargs):
        """
        Colbert-Miller kinetic energy for the [0, pi] range
        :param grid:
        :type grid:
        :param mass:
        :type mass:
        :param hb:
        :type hb:
        :param kw:
        :type kw:
        :return:
        :rtype:
        """

        coeff = hb ** 2 / (2 * mass)
        divs = len(grid)
        n = divs+1
        ke = np.empty((divs, divs))

        # we'll loop to fill to start
        ke[np.diag_indices_from(ke)] = (2*(n**2) + 1)/3 - 1/(np.sin(grid)^2)

        # then we start with constant diagonal bands like
        # in the previous cases but we also add on terms that
        # are constant along the _antidiagonal_ bands
        col_rng = np.arange(1, divs + 1)
        row_rng = np.arange(0, divs)
        for i in range(1, divs):
            col_inds = col_rng[i - 1:-1]
            row_inds = row_rng[:-i]
            diag_val = 1 /(np.sin(np.pi*i/(2*n))**2)
            ke[row_inds, col_inds] = diag_val
            ke[col_inds, row_inds] = diag_val

        # we fill in the antidiagonal bands
        for i in range(2, 2*divs):
            #antidiag terms
            row_inds = row_rng[:i-1] if i <= divs else row_rng[i-divs-1:]
            antdiag_val = -1 /(np.sin(np.pi*i/(2*n))**2)
            ke[row_inds, np.flip(row_inds)] += antdiag_val
            # ke[flip_cols, row_inds] += val

        # finally we do the phases (hard to get in with the antidiagonals)
        for i in range(1, divs):
            col_inds = col_rng[i - 1:-1]
            row_inds = row_rng[:-i]
            phase = (-1) ** (i)
            ke[row_inds, col_inds] *= phase
            ke[col_inds, row_inds] *= phase

        return coeff * ke

class RadialDVR(BaseDVR):
    """
    Provides a DVR for working on the (0, inf) range from Colbert and Miller
    """

    def get_grid(self, domain=(0, np.pi), divs=None, **kwargs):
        """
        Provides the grid appropriate for the Colbert-Miller (0, Pi) range

        :param domain:
        :type domain:
        :param divs:
        :type divs:
        :param kwargs:
        :type kwargs:
        :return:
        :rtype:
        """

        return domain[0] + (domain[1] - domain[0]) * np.arange(1, divs) / (divs + 1)

    def get_kinetic_energy(self, grid=None, mass=None, hb=1, **kwargs):
        """
        Colbert-Miller kinetic energy for the (0, inf) range

        :param grid:
        :type grid:
        :param mass:
        :type mass:
        :param hb:
        :type hb:
        :param kw:
        :type kw:
        :return:
        :rtype:
        """

        dx = grid[1] - grid[0]
        coeff = (hb ** 2) / (2 * mass * (dx ** 2))
        divs = len(grid)
        ke = np.empty((divs, divs))

        # we'll loop to fill to start
        ke[np.diag_indices_from(ke)] = np.pi**2/3 - 1/np.arange(1, divs+1)

        # then we start with constant diagonal bands like
        # in the previous cases but we also add on terms that
        # are constant along the _antidiagonal_ bands
        col_rng = np.arange(1, divs + 1)
        row_rng = np.arange(0, divs)
        for i in range(1, divs):
            col_inds = col_rng[i - 1:-1]
            row_inds = row_rng[:-i]
            diag_val = 2/(i**2)
            ke[row_inds, col_inds] = diag_val
            ke[col_inds, row_inds] = diag_val

        # we fill in the antidiagonal bands
        for i in range(2, 2*divs):
            #antidiag terms
            row_inds = row_rng[:i-1] if i <= divs else row_rng[i-divs-1:]
            antdiag_val = -2/(i**2)
            ke[row_inds, np.flip(row_inds)] += antdiag_val
            # ke[flip_cols, row_inds] += val

        # finally we do the phases (hard to get in with the antidiagonals)
        for i in range(1, divs):
            col_inds = col_rng[i - 1:-1]
            row_inds = row_rng[:-i]
            phase = (-1) ** (i)
            ke[row_inds, col_inds] *= phase
            ke[col_inds, row_inds] *= phase

        return coeff * ke