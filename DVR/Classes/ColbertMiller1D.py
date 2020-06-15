'''1D Colbert and Miller DVR

Note that we only supply a kinetic energy and a grid as the default in the DVR
class will handle everything else

http://xbeams.chem.yale.edu/~batista/v572/ColbertMiller.pdf
'''

import numpy as np, math

def grid(domain=(-5, 5), divs=10, **kw):
    '''Calculates a 1D grid'''

    return np.linspace(*domain, divs)

def kinetic_energy(grid=None, m=1, hb=1, **kw):
    '''Computes the kinetic energy for the grid'''

    dx=grid[1]-grid[0] # recomputed here simply to decouple the calling from dvr_grid
    divs=len(grid)
    ke=np.empty((divs, divs))

    coeff=(hb**2)/(2*m*(dx**2))
    # compute the band values for the first row
    b_val_0 = coeff*(math.pi**2)/3
    col_rng = np.arange(1, divs+1) # the column indices -- also what will be used for computing the off diagonal bands
    row_rng = np.arange(0, divs) # the row indices -- computed once and sliced
    b_vals = coeff * ((-1)**col_rng) * 2 / (col_rng**2)

    for i in range(divs):
        if i == 0:
            np.fill_diagonal(ke, b_val_0)
        else:
            col_inds = col_rng[i-1:-1]#+(i-1)
            row_inds = row_rng[:-i]
            ke[row_inds, col_inds] = b_vals[i-1]
            ke[col_inds, row_inds] = b_vals[i-1]

    return ke
