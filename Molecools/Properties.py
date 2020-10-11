"""
A collection of methods used in computing molecule properties
"""
import numpy as np, scipy.sparse as sp, itertools as ip
import McUtils.Numputils as nput
from McUtils.Coordinerds import CoordinateSet
from McUtils.Data import AtomData, UnitsData, BondData

from .Molecule import Molecule
from .Transformations import MolecularTransformation

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

        return np.tensordot(masses / np.sum(masses), coords, axes=[0, -2])

    @classmethod
    def center_of_mass(cls, mol):
        """Computes the moments of inertia

        :param mol:
        :type mol: Molecule
        :return:
        :rtype:
        """

        return cls.get_prop_center_of_mass(mol.coords, mol.masses)

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

        d = np.zeros(coords.shape[:-1] + (3, 3), dtype=float)
        diag = nput.vec_dots(coords, coords)
        d[..., (0, 1, 2), (0, 1, 2)] = diag[..., np.newaxis]
        o = nput.vec_outer(coords, coords, axes=[-1, -1])
        tens = np.tensordot(masses, d - o, axes=[0, -3])

        # print(masses)

        return tens

    @classmethod
    def inertia_tensor(cls, mol):
        """
        Computes the inertia tensors for the stored geometries

        :param mol:
        :type mol: Molecule
        :return:
        :rtype:
        """

        return cls.get_prop_inertia_tensors(mol.coords, mol.masses)

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
        from McUtils.Numputils import vec_crosses

        if coords.ndim == 1:
            raise ValueError("can't get moment of inertia for single point (?)")
        elif coords.ndim == 2:
            multiconfig = False
            coords = coords[np.newaxis]
            extra_shape = None
        else:
            multiconfig = True
            extra_shape = coords.shape[:-2]
            coords = coords.reshape((np.product(extra_shape),) + coords.shape[-2:])

        massy_doop = cls.get_prop_inertia_tensors(coords, masses)
        moms, axes = np.linalg.eigh(massy_doop)
        a = axes[:, :, 0]
        c = axes[:, :, 2]
        b = vec_crosses(a, c) # force right-handedness because we can
        axes[:, :, 1] = b # ensure we have true rotation matrices
        dets = np.linalg.det(axes)
        axes[:, :, 1] *= dets # ensure we have true rotation matrices
        if multiconfig:
            moms = moms.reshape(extra_shape + (3,))
            axes = axes.reshape(extra_shape + (3, 3))
        else:
            moms = moms[0]
            axes = axes[0]
        return moms, axes

    @classmethod
    def moments_of_inertia(cls, mol):
        """
        Computes the moments of inertia

        :param mol:
        :type mol: Molecule
        :return:
        :rtype:
        """

        return cls.get_prop_moments_of_inertia(mol.coords, mol.masses)

    @classmethod
    def get_prop_principle_axis_rotation(cls, coords, masses, sel=None, inverse=False):
        """
        Generates the principle axis transformation for a set of coordinates and positions

        :param coords:
        :type coords: CoordinateSet
        :param masses:
        :type masses: np.ndarray
        :return:
        :rtype: MolecularTransformation | List[MolecularTransformation]
        """

        multiconf = coords.multiconfig
        transforms = [None]*(1 if not multiconf else len(coords))
        if sel is not None:
            coords = coords[..., sel, :]
            masses = masses[sel]
        if multiconf:
            coords = list(coords)
        else:
            coords = [coords]
        mass = masses
        for i,c in enumerate(coords):
            com = cls.get_prop_center_of_mass(c, mass)
            transf = MolecularTransformation(-com)
            c = transf(c)
            moms, axes = cls.get_prop_moments_of_inertia(c, mass)
            if inverse:
                axes = axes.T
            transf = MolecularTransformation(axes)(transf)
            transforms[i] = transf

        if not multiconf:
            transforms = transforms[0]

        return transforms

    @classmethod
    def principle_axis_transformation(cls, mol, sel=None, inverse=False):
        """
        Generates the principle axis transformation for a Molecule

        :param mol:
        :type mol: Molecule
        :return:
        :rtype:
        """
        return cls.get_prop_principle_axis_rotation(mol.coords, mol.masses, sel=sel, inverse=inverse)

    @classmethod
    def get_prop_eckart_transformation(cls, masses, ref, coords, sel=None, inverse=False):
        """
        Computes Eckart transformations for a set of coordinates

        :param masses:
        :type masses: np.ndarray
        :param ref:
        :type ref: CoordinateSet
        :param coords:
        :type coords: CoordinateSet
        :return:
        :rtype: MolecularTransformation | List[MolecularTransformation]
        """

        # first we'll put everything in a principle axis frame
        multiconf = coords.multiconfig
        if sel is not None:
            coords = coords[..., sel, :]
            masses = masses[sel]
            ref = ref[..., sel, :]
        transforms = cls.get_prop_principle_axis_rotation(coords, masses)
        if multiconf:
            coords = list(coords)
        else:
            coords = [coords]
        if not multiconf:
            transforms = [transforms]

        ref_transf = cls.get_prop_principle_axis_rotation(ref, masses)
        ref = ref_transf.apply(ref)
        planar_ref = np.allclose(ref[:, 2], 0.)

        for i, ct in enumerate(zip(coords, transforms)):
            c = ct[0]
            t = ct[1] # type: MolecularTransformation
            c = t.apply(c)

            planar_struct = True if planar_ref else np.allclose(c[:, 2], 0.)
            is_planar = planar_ref or planar_struct
            if not is_planar:
                # generate pair-wise product matrix
                A = np.tensordot(masses/np.sum(masses), ref[:, :, np.newaxis] * c[:, np.newaxis, :], axes=[0, 0])
                # take SVD of this
                U, S, V = np.linalg.svd(A)
                rot = U @ V
            else:
                # generate pair-wise product matrix but only in 2D
                A = ref[:, :2, np.newaxis] * c[:, np.newaxis, :2]
                U, S, V = np.linalg.svd(A)
                rot = np.identity(3, np.float)
                rot[:2, :2] = U @ V.T

            if inverse:
                rot = rot.T
            transf = MolecularTransformation(t, rot)
            transforms[i] = transf

        if not multiconf:
            transforms = transforms[0]

        return transforms

    @classmethod
    def eckart_transformation(cls, mol, ref_mol, sel=None, inverse=False):
        """

        :param ref_mol: reference geometry
        :type ref_mol: Molecule
        :param mol: molecules to get Eckart embeddings for
        :type mol: Molecule
        :param sel: coordinate selection to use when doing the Eckart stuff
        :type mol:
        :return:
        :rtype:
        """
        m1 = ref_mol.masses
        m2 = mol.masses
        if not np.all(m1 == m2):
            raise ValueError("Eckart reference has different masses from scan ({}) vs. ({})".format(
                m1,
                m2
            ))
        return cls.get_prop_eckart_transformation(m1, ref_mol.coords, mol.coords, sel=sel, inverse=inverse)

    @classmethod
    def get_prop_translation_rotation_eigenvectors(cls, coords, masses):
        """
        Returns the eigenvectors corresponding to translations and rotations
        in the system

        :param coords:
        :type coords:
        :param masses:
        :type masses:
        :return:
        :rtype:
        """

        n = len(masses)
        mT = np.sum(masses)
        mvec = np.sqrt(masses)
        com = cls.get_prop_center_of_mass(coords, masses)
        shift_crds = coords - com[np.newaxis, :]
        M = np.kron(mvec/mT, np.eye(3)).T # translation eigenvectors
        mom_rot, ax_rot = cls.get_prop_moments_of_inertia(shift_crds, masses)
        inv_rot_2 = np.dot(ax_rot.T, np.diag(1/np.sqrt(mom_rot)), ax_rot)
        e = nput.levi_cevita3
        R = np.dot(shift_crds, np.dot(inv_rot_2, e).transpose((1, 2, 0))).reshape((3*n, 3)) # rotations
        freqs = np.concatenate([
            [1e-14, 1e-14, 1e-14],
            UnitsData.convert("Wavenumbers", "Hartrees")*(1/mom_rot)
            # this isn't right, I'm totally aware, but I think the frequency is supposed to be zero anyway and this
            # will be tiny
        ])
        eigs = np.concatenate([M, R], axis=1)

        np.savetxt("/Users/Mark/Desktop/tr_crd.data", coords)
        np.savetxt("/Users/Mark/Desktop/tr_mass.data", masses)
        np.savetxt("/Users/Mark/Desktop/tr_eigv.data", eigs)

        # import McUtils.Plots as plt
        # plt.ArrayPlot(eigs).show()
        # raise Exception([M, R])
        return freqs, eigs
    @classmethod
    def translation_rotation_eigenvectors(cls, mol, sel=None):
        """

        :param mol: molecules to get eigenvectors for
        :type mol: Molecule
        :param sel: coordinate selection to use when doing the rotation/translation calculations
        :type mol:
        :return:
        :rtype:
        """
        return cls.get_prop_translation_rotation_eigenvectors(mol.coords, mol.masses)

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
        # TODO: finish this off for real

    @classmethod
    def guessed_bonds(cls, mol, tol=1.05, guess_type=True):
        """
        Guesses the bonds for the molecule by finding the ones that are less than some percentage of a single bond for that
        pair of elements

        :return:
        :rtype:
        """

        if mol.multiconfig:
            coords = list(np.asarray(mol.coords))
        else:
            coords = [np.asarray(mol.coords)]
        guessed_bonds = [None] * len(coords)
        for i, coord in enumerate(coords):
            # TODO: generalize this to work for multiple configurations at once
            cds = coord[:, np.newaxis]
            dist_mat = np.linalg.norm(cds - cds.transpose(1, 0, 2), axis=2)
            atoms = np.array([a["ElementSymbol"] for a in mol._ats])
            pair_dists = np.array(
                [BondData.get_distance((a1, a2, 1), default=-1) for a1, a2 in ip.product(atoms, atoms)]).reshape(
                len(atoms), len(atoms)
            )
            pair_dists[np.tril_indices_from(pair_dists)] = -1

            pair_dists *= tol
            pos = np.array(np.where(dist_mat < pair_dists)).T

            if len(pos) == 0:
                return []

            if guess_type:
                bonds = [None] * len(pos)
                for i, ats in enumerate(pos):
                    a1, a2 = ats
                    test_data = BondData[atoms[a1], atoms[a2], None]
                    key = None
                    ref = 100
                    dist = dist_mat[a1, a2]
                    for t, d in test_data.items():
                        if dist < tol * d and d < ref:
                            ref = d
                            key = t
                    if key == "Single":
                        key = 1
                    elif key == "Double":
                        key = 2
                    elif key == "Triple":
                        key = 3
                    bonds[i] = [a1, a2, key]
            else:
                bonds = np.column_stack(pos, np.full((len(pos), 1), 1))

            guessed_bonds[i] = bonds

        # TODO: should maybe put some valence checker in here?
        if not mol.multiconfig:
            guessed_bonds = guessed_bonds[0]

        return guessed_bonds

    @classmethod
    def get_prop_chemical_formula(cls, atoms):
        """

        :param atoms:
        :type atoms: list[str]
        :return:
        :rtype:
        """
        return "".join(a+str(atoms.count(a)) for a in set(atoms))
    @classmethod
    def chemical_formula(cls, mol):
        """

        :param mol:
        :type mol: Molecule
        :return:
        :rtype:
        """
        return cls.get_prop_chemical_formula(mol.atoms)

