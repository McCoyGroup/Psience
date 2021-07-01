"""
ND Colbert and Miller DVR on [-inf, inf] range
but the basic template can be directly adapted to the
[0, 2pi] one or really any DVR that is a direct product
of 1D DVRs
"""

import numpy as np
import scipy.sparse as sp
import ColbertMiller1D as cm1D

def grid(domain=None, divs=None, flavor='[-inf,inf]', **kw):
    """

    :param domain:
    :type domain:
    :param divs:
    :type divs:
    :param kw:
    :type kw:
    :return:
    :rtype:
    """

    subgrids = [cm1D.grid(domain=dom, divs=div, flavor=flavor) for dom,div in zip(domain, divs)]
    mesh = np.array(np.meshgrid(*subgrids, indexing='ij'))

    rolly_polly_OLLY = np.roll(np.arange(len(mesh.shape)), -1)
    MEHSH = mesh.transpose(rolly_polly_OLLY)
    # for i in range(mesh.shape[0]):
    #     mesh = mesh.swapaxes(i, i+1)
    return MEHSH

def kinetic_energy(grid=None, mass=1, hb=1, g=None, g_deriv=None, flavor='[-inf,inf]', **kw):
    '''Computes n-dimensional kinetic energy for the grid'''
    from functools import reduce

    ndims = grid.shape[-1]
    try:
        iter(mass); ms = mass
    except TypeError:
        ms = [mass] * ndims

    try:
        iter(hb); hbs = hb
    except TypeError:
        hbs = [hb]*ndims

    ndim = grid.shape[-1]
    grids = [ # build subgrids
        grid[(0, )*i + (...,) + (0, ) * (ndim-i-1) +(i,)]
        for i in range(ndim)
    ]
    if g is not None:
        if g_deriv is None:
            raise ValueError("if functions for `g` are supplied, also need functions, `g_deriv` for the second derivative of `g`")

        include_coupling = any(
            i != j
            and not (
                isinstance(g[i][j], (int, float, np.integer, np.floating))
                and g[i][j] == 0
            )
            for i in range(ndim) for j in range(ndim)
        )

        ms = [1/2]*len(ms)

    else:
        include_coupling = False

    kes = [
        cm1D.kinetic_energy(subg, mass=m, hb=hb, flavor=flavor)
        for subg, m, hb in zip(grids, ms, hbs)
    ]
    kes = [sp.csr_matrix(mat) for mat in kes]
    if g is None: # we passed constant masses
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
                    sp.eye(tot_shape[k]) if k != i else sp.csr_matrix(kes[k]) for k in range(ndim)
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

        if include_coupling:
            momenta = [cm1D.real_momentum(subg, hb=hb, flavor=flavor) for subg in grids]
            kinetic_coupling = sp.csr_matrix(ke.shape, dtype=ke.dtype) # initialize empty tensor
            for i in range(len(momenta)): # build out all of the coupling term products
                for j in range(i+1, len(momenta)):
                    if not (isinstance(g[i][j], (int, float, np.integer, np.floating)) and g[i][j] == 0):
                        # evaluate g over the terms and average
                        g_vals = np.reshape(g[i][i](flat_grid), grid.shape[:-1])

                        # construct the basic momenta kronecker product
                        sub_momenta = [ # set up all the subtensors we'll need for this
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

                        # and we pull the G matrix values for the corresponding i and j indices
                        row_vals = g_vals[row_inds]
                        col_vals = g_vals[col_inds]

                        # finally we take the average of the two and put them into a sparse matrix
                        # that can be multiplied by the base momentum matrix values
                        avg_g_vals = 1/2*(row_vals + col_vals)
                        coupling_term = sp.csr_matrix(
                            (
                                avg_g_vals * mom_prod_vals,
                                (flat_rows, flat_cols)
                            ),
                            shape=momentum_mat.shape,
                            dtype=momentum_mat.dtype
                        )

                        kinetic_coupling -= coupling_term # negative sign from the two factors of i
            ke += kinetic_coupling
            # print(ke.getnnz(), np.prod(ke.shape))
            if ke.getnnz() >= 1/2*np.prod(ke.shape):
                ke = ke.toarray()

    return ke

def wavefunctions(hamiltonian=None, num_wfns=10, diag_mode='sparse', **kw):
    """Computes the wavefunctions using sparse methods"""
    if isinstance(hamiltonian, sp.spmatrix) and diag_mode=='dense':
        hamiltonian = hamiltonian.toarray()
    if isinstance(hamiltonian, sp.spmatrix):
        import scipy.sparse.linalg as la
        return la.eigsh(hamiltonian, num_wfns, which='SM')
    else:
        engs, wfns = np.linalg.eigh(hamiltonian)
        # print(engs[:num_wfns])
        return (engs[:num_wfns], wfns[:, :num_wfns])


