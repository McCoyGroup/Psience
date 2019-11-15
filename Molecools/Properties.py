"""
A collection of methods used in computing molecule properties
"""
import numpy as np

def mommy_tensors(coords, masses):
    """Computes the moment of intertia tensors for the walkers with coordinates coords (assumes all have the same masses)

    :param coords:
    :type coords: np.ndarray
    :param masses:
    :type masses: np.ndarray
    :return:
    :rtype:
    """

    x = coords[:, :, 0]
    y = coords[:, :, 1]
    z = coords[:, :, 2]
    o = np.zeros(x.shape)
    # build the skew matrices that we matrix product to get the individual r^2 tensors
    doop = np.array(
        [
            [o, -z, y],
            [z, o, -x],
            [-y, x, o]
        ]
    )

    # current shape (3, 3, N, M) where N is number of walkers and M is number of atoms
    # make it into (N, M, 3, 3)
    doop = doop.transpose(np.roll(np.arange(4), -2)) # might need to be a 2 :)

    # build the inertia products
    # from the np.dot docs:
    #   if a is an N-D array and b is an M-D array (where M>=2), it is a sum product over the last axis of a and the second-to-last axis of b:
    # this means for our (N, M, 3, 3)
    doopdoop = np.matmul(doop, doop)

    # dot in the masses to reduce the dimension of this to (N, 3, 3)
    massy_doop = np.tensordot(masses, doopdoop, axes=(0, 1))

    return massy_doop

def mommies(coords, masses):
    """Computes the moment of intertia tensor for the walkers with coordinates coords (assumes all have the same masses)

    :param coords:
    :type coords: np.ndarray
    :param masses:
    :type masses: np.ndarray
    :return:
    :rtype:
    """

    massy_doop = mommy_tensors(coords, masses)
    return np.linalg.eigh(massy_doop)



