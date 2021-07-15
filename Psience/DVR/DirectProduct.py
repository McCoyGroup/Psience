"""
Provides a direct-product extension to a system of 1D DVRs
"""

import numpy as np, scipy.sparse as sp
from functools import reduce
from .BaseDVR import BaseDVR
from .ColbertMiller import *

__all__ = [
    "DirectProductDVR",
    "CartesianNDDVR",
    "RingNDDVR",
    "SphericalDVR"
]

class DirectProductDVR(BaseDVR):

    def __init__(self,
                 dvrs_1D,
                 **base_opts
                 ):
        """
        :param dvrs_1D: a series of 1D DVRs that can provide the inputs we'll product together
        :type dvrs_1D: Iterable[AbstractDVR]
        :param base_opts:
        :type base_opts:
        """
        super().__init__(**base_opts)
        self.dvrs = dvrs_1D
    def __repr__(self):
        return "{}({}, pot={})".format(
            type(self).__name__,
            self.dvrs,
            self.potential_function
        )

    def get_grid(self, domain=None, divs=None, **kwargs):

        if domain is None:
            subgrids = [dvr.grid() for dvr in self.dvrs]
        else:
            subgrids = [dvr.grid(domain=dom, divs=div) for dvr,dom,div in zip(self.dvrs, domain, divs)]
        mesh = np.array(np.meshgrid(*subgrids, indexing='ij'))
        MEHSH = np.moveaxis(mesh, 0, len(subgrids))

        return MEHSH

    def grid(self, domain=None, divs=None, **kwargs):
        if domain is None:
            domain = self.domain
        if divs is None:
            divs = self.divs

        # if divs is None:
        #     raise ValueError("need a value for `domain`")
        # if divs is None:
        #     raise ValueError("need a value for `divs`")
        return self.get_grid(domain=domain, divs=divs, **kwargs)

    def get_kinetic_energy(self, grid=None, mass=None, hb=1, g=None, g_deriv=None, **kwargs):

        ndims = grid.shape[-1]
        try:
            iter(mass); ms = mass
        except TypeError:
            ms = [mass] * ndims

        try:
            iter(hb); hbs = hb
        except TypeError:
            hbs = [hb] * ndims

        ndim = grid.shape[-1]
        grids = [  # build subgrids
            grid[(0,) * i + (...,) + (0,) * (ndim - i - 1) + (i,)]
            for i in range(ndim)
        ]
        if g is not None:
            if g_deriv is None:
                raise ValueError(
                    "if functions for `g` are supplied, also need functions, `g_deriv` for the second derivative of `g`")

            include_coupling = any(
                i != j
                and not (
                        isinstance(g[i][j], (int, float, np.integer, np.floating))
                        and g[i][j] == 0
                )
                for i in range(ndim) for j in range(ndim)
            )

            ms = [1] * len(ms)

        else:
            include_coupling = False

        kes = [
            dvr.kinetic_energy(subg, mass=m, hb=hb)
            for dvr, subg, m, hb in zip(self.dvrs, grids, ms, hbs)
        ]
        kes = [sp.csr_matrix(mat) for mat in kes]
        if g is None:  # we passed constant masses
            def _kron_sum(a, b):
                '''Computes a Kronecker sum to build our Kronecker-Delta tensor product expression'''
                n_1 = a.shape[0]
                n_2 = b.shape[0]
                ident_1 = sp.identity(n_1)
                ident_2 = sp.identity(n_2)

                return sp.kron(a, ident_2) + sp.kron(ident_1, b)

            ke = reduce(_kron_sum, kes)
        else:
            flat_grid = np.reshape(grid, (-1, ndim))
            tot_shape = [len(gr) for gr in grids]  # we'll need this to multiply terms into the right indices
            ke = sp.csr_matrix((len(flat_grid), len(flat_grid)), dtype=kes[0].dtype)
            for i in range(ndim):  # build out all of the coupling term products
                # evaluate g over the terms and average
                if not (
                        isinstance(g[i][i], (int, float, np.integer, np.floating))
                        and g[i][i] == 0
                ):
                    g_vals = np.reshape(g[i][i](flat_grid), grid.shape[:-1])

                    # construct the basic kinetic energy kronecker product
                    sub_kes = [  # set up all the subtensors we'll need for this
                        sp.eye(tot_shape[k]) if k != i else kes[k] for k in range(ndim)
                    ]
                    ke_mat = reduce(sp.kron, sub_kes)

                    # now we need to figure out where to multiply in the g_vals
                    flat_rows, flat_cols, ke_vals = sp.find(ke_mat)
                    # we convert each row and column into its corresponding direct
                    # product index since each is basically a flat index for a multdimensional
                    # array
                    row_inds = np.unravel_index(flat_rows, tot_shape)
                    col_inds = np.unravel_index(flat_cols, tot_shape)

                    # and we pull the G matrix values for the corresponding i and j indices
                    row_vals = g_vals[row_inds]
                    col_vals = g_vals[col_inds]

                    # finally we take the average of the two and put them into a sparse matrix
                    # that can be multiplied by the base kinetic matrix values
                    avg_g_vals = 1 / 2 * (row_vals + col_vals)
                    ke_mat = sp.csr_matrix(
                        (
                            avg_g_vals * ke_vals,
                            (flat_rows, flat_cols)
                        ),
                        shape=ke.shape,
                        dtype=ke_vals.dtype
                    )

                    ke += ke_mat

                    gd_vals = 1 / 4 * g_deriv[i](flat_grid)
                    ke_contrib = sp.diags([gd_vals], [0])
                    ke += ke_contrib

            if include_coupling:
                momenta = [dvr.real_momentum(subg, hb=hb) for dvr,subg in zip(self.dvrs, grids)]
                kinetic_coupling = sp.csr_matrix(ke.shape, dtype=ke.dtype)  # initialize empty tensor
                for i in range(len(momenta)):  # build out all of the coupling term products
                    for j in range(i + 1, len(momenta)):
                        if i!= j and not (isinstance(g[i][j], (int, float, np.integer, np.floating)) and g[i][j] == 0):
                            g_vals = np.reshape(g[i][j](flat_grid), grid.shape[:-1])

                            # construct the basic momenta Kroenecker product
                            sub_momenta = [  # set up all the subtensors we'll need for this
                                sp.eye(tot_shape[k]) if k != i and k != j else sp.csr_matrix(momenta[k])
                                for k in range(len(momenta))
                            ]
                            momentum_mat = reduce(sp.kron, sub_momenta)

                            # now we need to figure out where to multiply in the ij_vals
                            flat_rows, flat_cols, mom_prod_vals = sp.find(momentum_mat)

                            # we convert each row and column into its corresponding direct
                            # product index since each is basically a flat index for a multdimensional
                            # array
                            row_inds = np.unravel_index(flat_rows, tot_shape)
                            col_inds = np.unravel_index(flat_cols, tot_shape)

                            # we do the swap of row/column indices required by the mixed term
                            swap_rows = row_inds[:i] + (col_inds[i],) + row_inds[i+1:]
                            swap_cols = col_inds[:i] + (row_inds[i],) + col_inds[i+1:]

                            # and we pull the G matrix values for the corresponding indices
                            row_vals = g_vals[swap_rows]
                            col_vals = g_vals[swap_cols]

                            # finally we take the sum of the two and put them into a sparse matrix
                            # that can be multiplied by the base momentum matrix values
                            sum_g_vals = (row_vals + col_vals)
                            coupling_term = sp.csr_matrix(
                                (
                                    1/2 * sum_g_vals * mom_prod_vals,
                                    (flat_rows, flat_cols)
                                ),
                                shape=momentum_mat.shape,
                                dtype=momentum_mat.dtype
                            )

                            kinetic_coupling -= coupling_term  # negative sign from the two factors of i
                ke += kinetic_coupling
                # print(ke.getnnz(), np.prod(ke.shape))
                if ke.getnnz() >= 1 / 2 * np.prod(ke.shape):
                    ke = ke.toarray()

        return ke

    def kinetic_energy(self, grid=None, mass=None, hb=1, g=None, g_deriv=None, **kwargs):
        """
        Computes the N-dimensional kinetic energy
        :param grid:
        :type grid:
        :param mass:
        :type mass:
        :param hb:
        :type hb:
        :param g:
        :type g:
        :param g_deriv:
        :type g_deriv:
        :return:
        :rtype:
        """

        if grid is None:
            grid = self.grid()

        return self.get_kinetic_energy(grid=grid, mass=mass, hb=hb, g=g, g_deriv=g_deriv, **kwargs)

class CartesianNDDVR(DirectProductDVR):
    """
    Provides an ND-DVR over different domains
    """
    def __init__(self,
                 domains,
                 **base_opts
                 ):

        sub_dvrs = [CartesianDVR(domain=dom[:2], divs=dom[2]) for dom in domains]
        super().__init__(
            sub_dvrs,
            **base_opts
        )

class RingNDDVR(DirectProductDVR):
    """
    Provides an ND-DVR for products of periodic (0, 2Pi) ranges
    """

    def __init__(self,
                 divs,
                 **base_opts
                 ):
        sub_dvrs = [RingDVR(domain=(0, 2*np.pi), divs=div) for div in divs]
        super().__init__(
            sub_dvrs,
            **base_opts
        )

class SphericalDVR(DirectProductDVR):
    def __init__(self,
                 r_max,
                 divs,
                 **base_opts
                 ):
        sub_dvrs = [
            RadialDVR(domain=(0, r_max), divs=divs[0]),
            RingDVR(domain=(0, 2*np.pi), divs=divs[1]),
            PolarDVR(domain=(0, 2*np.pi), divs=divs[1]),


        ]
        super().__init__(
            sub_dvrs,
            **base_opts
        )

