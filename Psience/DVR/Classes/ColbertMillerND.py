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

def kinetic_energy(grid=None, m=1, hb=1, g=None, g_deriv=None, flavor='[-inf,inf]', **kw):
    '''Computes n-dimensional kinetic energy for the grid'''
    from functools import reduce

    ndims = grid.shape[-1]
    try:
        iter(m); ms = m
    except TypeError:
        ms = [m]*ndims

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

        g_diag = [g[i][i] for i in range(len(grids))]

        include_coupling = any(
            i != j
            and not isinstance(g[i][j], (int, float, np.integer, np.floating))
            and g[i][j] != 0
            for i in range(ndim) for j in range(ndim)
        )

    else:
        include_coupling = False
        g_diag = [None] * ndims
        g_deriv = [None] * ndims

    kes = [cm1D.kinetic_energy(subg, m=m, hb=hb, g=gv, g_deriv=gdv, flavor=flavor) for subg, m, hb, gv, gdv in zip(grids, ms, hbs, g_diag, g_deriv)]

    kes = [sp.csr_matrix(mat) for mat in kes]
    def _kron_sum(a, b):
        '''Computes a Kronecker sum to build our Kronecker-Delta tensor product expression'''
        n_1 = a.shape[0]
        n_2 = b.shape[0]
        ident_1 = sp.identity(n_1)
        ident_2 = sp.identity(n_2)

        return sp.kron(a, ident_2) + sp.kron(ident_1, b)

    ke = reduce(_kron_sum, kes)

    # if g_diag[0] is not None:
    #     raise Exception(ke.nnz)

    if include_coupling:
        raise Exception('wat')
        momenta = [sp.csr_matrix(cm1D.real_momentum(subg, hb=hb, flavor=flavor)) for subg in grids]
        kinetic_coupling = sp.csr_matrix(ke.shape, ke.dtype) # initialize empty tensor
        for i in range(len(momenta)): # build out all of the coupling term products
            for j in range(i+1, len(momenta)):
                if not isinstance(g[i][j], (int, float, np.integer, np.floating)) and g == 0:
                    # evaluate g over the terms
                    ij_grid = np.array(np.meshgrid(grids[i], grids[j], indexing='ij')).transpose(2, 0, 1)
                    raise Exception(ij_grid.shape)
                    ij_vals = g[i][j](ij_grid.reshape(-1, 2))

                    sub_momenta = [ # set up all the subtensors we'll need for this
                        sp.eye(len(grids[k])) if k != i and k != j else momenta[k]
                        for k in range(len(momenta))
                    ]
                    kinetic_coupling += -g_mats[i, j] * reduce(sp.kron, sub_momenta) # negative sign from the two factors of i
        ke += kinetic_coupling

    return ke


def potential_energy(grid=None, potential_function=None, potential_values=None, **kw):

    if potential_values is not None:
        return sp.diags([potential_values], [0])
    else:
        from functools import reduce
        from operator import mul

        npts = reduce(mul, grid.shape[:-1], 1)
        gps = np.reshape(grid, (npts, grid.shape[-1]))
        pots = potential_function(gps)
        return sp.diags([pots], [0])

def wavefunctions(hamiltonian=None, num_wfns=10, **kw):
    """Computes the wavefunctions using sparse methods"""
    # if isinstance(hamiltonian, sp.spmatrix):
    #     hamiltonian = hamiltonian.toarray()
    if isinstance(hamiltonian, sp.spmatrix):
        import scipy.sparse.linalg as la
        return la.eigsh(hamiltonian, num_wfns, which='SM')
    else:
        engs, wfns = np.linalg.eigh(hamiltonian)
        # print(engs[:num_wfns])
        return (engs[:num_wfns], wfns[:, :num_wfns])


