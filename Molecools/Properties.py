"""
A collection of methods used in computing molecule properties
"""
import numpy as np
from McUtils.Coordinerds import CoordinateSet, CartesianCoordinateSystem, CartesianCoordinates3D
from .Molecule import Molecule

__all__ = [
    "MolecularProperties"
]

class MolecularProperties:
    """
    An object whose sole purpose in life is to get molecular properties
    A property should be implemented in two parts:
        1) a classmethod called get_prop_<prop name> that takes raw inputs and uses them to compute a property
        2) a classmethod called get_<prop name> that extracts the property from a passed molecule
        3) an @property called <prop name> that passes the bound molecule to get_<prop name>
    """
    def __init__(self, molecule):
        self.molecule = molecule

    @classmethod
    def get_prop_center_of_mass(cls, coords, masses):
        """Gets the center of mass for the coordinates

        :param coords:
        :type coords: CoordinateSet
        :param masses:
        :type masses:
        :return:
        :rtype:
        """
        return np.tensordot(masses, coords, axes=[0, -2]) / np.sum(masses)

    @classmethod
    def get_prop_inertia_tensors(cls, coords, masses):
        """Computes the moment of intertia tensors for the walkers with coordinates coords (assumes all have the same masses)

            :param coords:
            :type coords: CoordinateSet
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
        doop = doop.transpose(np.roll(np.arange(4), -2))  # might need to be a 2 :)

        # build the inertia products
        # from the np.dot docs:
        #   if a is an N-D array and b is an M-D array (where M>=2), it is a sum product over the last axis of a and the second-to-last axis of b:
        # this means for our (N, M, 3, 3)
        doopdoop = np.matmul(doop, doop)

        # dot in the masses to reduce the dimension of this to (N, 3, 3)
        massy_doop = np.tensordot(masses, doopdoop, axes=(0, 1))

        return massy_doop

    @classmethod
    def get_prop_moments_of_inertia(cls, coords, masses):
        """Computes the moment of inertia tensor for the walkers with coordinates coords (assumes all have the same masses)

        :param coords:
        :type coords: np.ndarray
        :param masses:
        :type masses: np.ndarray
        :return:
        :rtype:
        """

        massy_doop = cls.get_prop_inertia_tensors(coords, masses)
        return np.linalg.eigh(massy_doop)
    @classmethod
    def get_moments_of_inertia(cls, mol):
        """Computes the moments of inertia

        :param mol:
        :type mol:
        :return:
        :rtype:
        """



