"""
A collection of methods used in computing molecule properties
"""
import numpy as np, scipy.sparse as sp
from McUtils.Coordinerds import CoordinateSet, CartesianCoordinateSystem, CartesianCoordinates3D
from .Molecule import Molecule

__all__ = [
    "MolecularProperties",
    "MolecularPropertyError"
]

class MolecularPropertyError(Exception):
    """
    General error class for MolecularProperties
    """

class MolecularProperties:
    """
    An object whose sole purpose in life is to get molecular properties
    A property should be implemented in two parts:
        1) a classmethod called get_prop_<prop name> that takes raw inputs and uses them to compute a property
        2) a classmethod called <prop name> that extracts the property from a passed molecule

    All properties should be appropriately vectorized to work on a single configuration or a set of configurations
    """

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
    def center_of_mass(cls, mol):
        """Computes the moments of inertia

        :param mol:
        :type mol: Molecule
        :return:
        :rtype:
        """

        return cls.get_prop_moments_of_inertia(mol.coords, mol.masses)

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
        :type coords: CoordinateSet
        :param masses:
        :type masses: np.ndarray
        :return:
        :rtype:
        """

        massy_doop = cls.get_prop_inertia_tensors(coords, masses)
        return np.linalg.eigh(massy_doop)

    @classmethod
    def moments_of_inertia(cls, mol):
        """Computes the moments of inertia

        :param mol:
        :type mol: Molecule
        :return:
        :rtype:
        """

        return cls.get_prop_moments_of_inertia(mol.coords, mol.masses)

    @classmethod
    def get_prop_adjacency_matrix(cls, atoms, bonds):
        """
        Returns the adjacency matrix for the molecule

        :param bonds:
        :type bonds: Iterable[int]
        :return:
        :rtype:
        """

        cons = np.array([x[:2] for x in bonds])
        adj = sp.csr_matrix(
            (np.ones((len(bonds), )), (cons[:, 0], cons[:, 1])),
            (len(atoms), len(atoms))
        )
        return adj + adj.T

    @classmethod
    def get_prop_connectivity(cls, atoms, bonds):
        """
        Returns the adjacency matrix for the molecule

        :param bonds:
        :type bonds: Iterable[int]
        :return:
        :rtype:
        """

        if isinstance(bonds, sp.spmatrix):
            adj_mat = bonds
        else:
            adj_mat = cls.get_prop_adjacency_matrix(atoms, bonds)

        return np.array(adj_mat.sum(axis=1))

    @classmethod
    def get_prop_fragments(cls, atoms, bonds):
        """
        Returns the fragments for the molecule

        :param bonds:
        :type bonds: Iterable[int]
        :return:
        :rtype:
        """

        if isinstance(bonds, sp.spmatrix):
            adj_mat = bonds
        else:
            adj_mat = cls.get_prop_adjacency_matrix(atoms, bonds)

        ngroups, labels = sp.csgraph.connected_components(adj_mat)
        def pull(l, i):
            w = np.where(labels==i)
            if len(w) > 0:
                w = w[0]
            return w
        return [pull(labels, i) for i in range(ngroups)]

    @classmethod
    def fragments(cls, mol):
        """

        :param mol:
        :type mol: Molecule
        :return:
        :rtype:
        """

        cds = mol.coords
        bonds = mol.bonds
        ats = mol.atoms

        comps = cls.get_prop_fragments(ats, bonds)
        bond_map = {}
        for k in bonds:
            if k[0] not in bond_map:
                bond_map[k[0]] = [k]
            else:
                bond_map[k[0]].append(k)
        for i in range(len(ats)):
            if i not in bond_map:
                bond_map[i] = []

        frags = [None]*len(comps)
        for i,g in enumerate(comps):
            frag_ats = [ats[x] for x in g]
            frag_cds = cds[g] if not cds.multiconfig else cds[:, g]
            frag_bonds = [bond_map[x] for x in g]
            frags[i] = Molecule(
                frag_ats,
                frag_cds,
                bonds=sum(frag_bonds, [])
            )
        return frags

    @classmethod
    def get_prop_zmat_ordering(cls, atoms, bonds):
        """
        Gets a guessed Z-matrix ordering for the molecule with connectivity defined by bonds based on the following:
            1. Fragments are separated out
            2. The atom with the highest degree of connectivity in each fragment is chosen as the fragment "label"
            3. Fragments are ordered by connectivity of the label from high to low
            4. Fragment labels reference each other with:
                a) the second label on the x-axis
                b) the 3rd in the xy-plane
                c) all others relative to the first three
            5. All other atoms are sorted first by fragment label, then by connection to the fragment label, and then by connectivity
            6. Atoms reference each other based on the following:
                a) if the atom has one bond:
                    i)   the atom it is bound to
                    ii)  the lowest-connectivity atom that one is bound to
                    iii) the second-lowest-connectivity atom OR the next fragment label
                b) if the atom has two bonds:
                    i)   the highest-connectivity atom it is bound to
                    ii)  the lowest-connectivity atom it is bound to
                    iii) the lowest-connectivity atom (i) is bound to
                c) if the atom has three bonds:
                    i)   the highest-connectivity atom it is bound to
                    ii)  the lowest-connectivity atom it is bound to
                    iii) the second-highest connectivity atom it is bound to
              if any of these atoms do not exist, the fragment labels will be used in their place

        :param bonds:
        :type bonds:
        :return:
        :rtype:
        """

        adj = cls.get_prop_adjacency_matrix(atoms, bonds)
        frags = cls.get_prop_fragments(atoms, bonds)
        conn = cls.get_prop_connectivity(atoms, adj)


