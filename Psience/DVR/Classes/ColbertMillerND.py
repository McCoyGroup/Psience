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

def kinetic_energy(grid=None, m=1, hb=1, flavor='[-inf,inf]', **kw):
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
    grids = [
        grid[(0, )*i + (...,) + (0, ) * (ndim-i-1) +(i,)]
        for i in range(ndim)
    ]
    kes = [cm1D.kinetic_energy(subg, m=m, hb=hb, flavor=flavor) for subg, m, hb in zip(grids, ms, hbs)]

    kes = [sp.csr_matrix(mat) for mat in kes]

    def _kron_sum(a, b):
        '''Computes a Kronecker sum to build our Kronecker-Delta tensor product expression'''
        n_1 = a.shape[0]
        n_2 = b.shape[0]
        ident_1 = sp.identity(n_1)
        ident_2 = sp.identity(n_2)

        return sp.kron(a, ident_2) + sp.kron(ident_1, b)

    ke = reduce(_kron_sum, kes)
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
    import scipy.sparse.linalg as la
    return la.eigsh(hamiltonian, num_wfns, which = 'SM')


