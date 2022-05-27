"""
A collection of methods used in computing molecule properties
"""
import numpy as np, scipy.sparse as sp, itertools as ip, os, abc

import McUtils.Numputils as nput
from McUtils.Coordinerds import CoordinateSet
from McUtils.ExternalPrograms import OpenBabelInterface
from McUtils.GaussianInterface import GaussianFChkReader, GaussianFChkReaderException
from McUtils.Data import AtomData, UnitsData, BondData
from McUtils.Zachary import FiniteDifferenceDerivative

from .MoleculeInterface import AbstractMolecule
from .Transformations import MolecularTransformation
from .Vibrations import MolecularVibrations, MolecularNormalModes

__all__ = [
    "MolecularProperties",
    "MolecularPropertyError",
    "OpenBabelMolManager",
    "DipoleSurfaceManager",
    "PotentialSurfaceManager",
    "NormalModesManager"
]

__reload_hook__ = [".MoleculeInterface", '.Vibrations', '.Transformations']

class MolecularPropertyError(Exception):
    """
    General error class for MolecularProperties
    """

class StructuralProperties:
    """
    The set of molecular properties
    that depend on its coordinates/configuration
    """

    @classmethod
    def get_prop_mass_weighted_coords(cls, coords, masses):
        """Gets the mass-weighted coordinates for the system

        :param coords:
        :type coords: CoordinateSet
        :param masses:
        :type masses:
        :return:
        :rtype:
        """

        sel = masses > 0
        comp = masses <= 0
        new_stuff = masses.copy()
        new_stuff[sel] = np.sqrt(masses[sel])
        new_stuff[comp] = 1
        return new_stuff[:, np.newaxis] * coords

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

        masses = masses.copy()
        masses[masses < 0] = 0

        return np.tensordot(masses / np.sum(masses), coords, axes=[0, -2])

    @classmethod
    def get_prop_inertia_tensors(cls, coords, masses):
        """
        Computes the moment of intertia tensors for the walkers with coordinates coords (assumes all have the same masses)

        :param coords:
        :type coords: CoordinateSet
        :param masses:
        :type masses: np.ndarray
        :return:
        :rtype:
        """

        com = cls.get_prop_center_of_mass(coords, masses)
        coords = coords - com[..., np.newaxis, :]

        real_spots = masses > 0 # allow for dropping out dummy atoms
        coords = coords[..., real_spots, :]
        masses = masses[real_spots]

        d = np.zeros(coords.shape[:-1] + (3, 3), dtype=float)
        diag = nput.vec_dots(coords, coords)
        d[..., (0, 1, 2), (0, 1, 2)] = diag[..., np.newaxis]
        # o = np.array([[np.outer(a, a) for a in r] for r in coords])
        o = nput.vec_outer(coords, coords, axes=[-1, -1])
        tens = np.tensordot(masses, d - o, axes=[0, -3])

        return tens

    @classmethod
    def get_prop_moments_of_inertia(cls, coords, masses):
        """
        Computes the moment of inertia tensor for the walkers with coordinates coords (assumes all have the same masses)

        :param coords:
        :type coords: CoordinateSet
        :param masses:
        :type masses: np.ndarray
        :return:
        :rtype:
        """

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
        # a = axes[..., :, 0]
        # c = axes[..., :, 2]
        # b = nput.vec_crosses(a, c)  # force right-handedness to avoid inversions
        # axes[..., :, 1] = b
        a = axes[..., :, 0]
        b = axes[..., :, 1]
        c = nput.vec_crosses(b, a)  # force right-handedness to avoid inversions
        axes[..., :, 2] = c
        dets = np.linalg.det(axes) # ensure we have true rotation matrices to avoid inversions
        axes[..., :, 2] /= dets[..., np.newaxis]

        if multiconfig:
            moms = moms.reshape(extra_shape + (3,))
            axes = axes.reshape(extra_shape + (3, 3))
        else:
            moms = moms[0]
            axes = axes[0]
        return moms, axes

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
        transforms = [None] * (1 if not multiconf else len(coords))
        com_coords = coords
        com_mass = masses
        if sel is not None:
            coords = coords[..., sel, :]
            masses = masses[sel]
        if multiconf:
            coords = list(coords)
        else:
            coords = [coords]
        mass = masses
        for i, (c, c2) in enumerate(zip(coords, com_coords)):
            com = cls.get_prop_center_of_mass(com_coords, com_mass)
            transf = MolecularTransformation(-com)
            c = transf(c)
            moms, axes = cls.get_prop_moments_of_inertia(c, mass)
            if not inverse:
                axes = axes.T
            transf = MolecularTransformation(axes)(transf)
            transforms[i] = transf

        if not multiconf:
            transforms = transforms[0]

        return transforms

    @classmethod
    def get_principle_axis_embedded_coords(cls, coords, masses):
        """
        Returns coordinate embedded in the principle axis frame

        :param coords:
        :type coords:
        :param masses:
        :type masses:
        :return:
        :rtype:
        """
        com = cls.get_prop_center_of_mass(coords, masses)
        # crds_ = coords
        coords = coords - com[..., np.newaxis, :]
        moms, pax_axes = cls.get_prop_moments_of_inertia(coords, masses)
        # pax_axes = np.swapaxes(pax_axes, -2, -1)
        coords = np.matmul(coords, pax_axes)

        return coords, com, pax_axes

    @classmethod
    def get_prop_principle_axis_data(cls, coords, masses):
        """
        Generates the principle axis transformation for a set of coordinates and positions

        :param coords:
        :type coords: CoordinateSet
        :param masses:
        :type masses: np.ndarray
        :return:
        :rtype: MolecularTransformation | List[MolecularTransformation]
        """

        coords, com, axes = cls.get_principle_axis_embedded_coords(coords, masses)
        return coords, com, axes

    planar_ref_tolerance=1e-6
    @classmethod
    def get_eckart_rotations(cls, masses, ref, coords, sel=None, in_paf=False, planar_ref_tolerance=None):
        """
        Generates the Eckart rotation that will align ref and coords, assuming initially that `ref` and `coords` are
        in the principle axis frame

        :param masses:
        :type masses:
        :param ref:
        :type ref:
        :param coords:
        :type coords: np.ndarray
        :return:
        :rtype:
        """

        if coords.ndim == 2:
            coords = np.broadcast_to(coords, (1,) + coords.shape)

        if planar_ref_tolerance is None:
            planar_ref_tolerance = cls.planar_ref_tolerance

        if not in_paf:
            coords, com, pax_axes = cls.get_principle_axis_embedded_coords(coords, masses)
            ref, ref_com, ref_axes = cls.get_principle_axis_embedded_coords(ref, masses)
            # raise ValueError(ref)
        else:
            com = pax_axes = None
            ref_com = ref_axes = None

        if sel is not None:
            coords = coords[..., sel, :]
            masses = masses[sel]
            ref = ref[..., sel, :]

        real_pos = masses > 0
        # print(real_pos)
        og_coords = coords
        coords = coords[..., real_pos, :]
        masses = masses[real_pos,]
        og_ref = ref
        ref = ref[..., real_pos, :]

        planar_ref = np.allclose(ref[:, 2], 0., atol=planar_ref_tolerance)
        if not planar_ref:
            # generate pair-wise product matrix
            A = np.tensordot(masses / np.sum(masses),
                             ref[np.newaxis, :, :, np.newaxis] * coords[:, :, np.newaxis, :],
                             axes=[0, 1])
            # take SVD of this
            U, S, V = np.linalg.svd(A)
            rot = np.matmul(U, V)
        else:
            # generate pair-wise product matrix but only in 2D
            F = ref[np.newaxis, :, :2, np.newaxis] * coords[:, :, np.newaxis, :2]
            A = np.tensordot(masses / np.sum(masses), F, axes=[0, 1])
            U, S, V = np.linalg.svd(A)
            rot = np.broadcast_to(np.eye(3, dtype=float), (len(coords), 3, 3)).copy()
            rot[..., :2, :2] = np.matmul(U, V)

        a = rot[..., :, 0]
        b = rot[..., :, 1]
        c = nput.vec_crosses(a, b, normalize=True)  # force right-handedness because we can
        rot[..., :, 2] = c  # ensure we have true rotation matrices
        dets = np.linalg.det(rot)
        rot[..., :, 2] /= dets[..., np.newaxis]  # ensure we have true rotation matrices

        # dets = np.linalg.det(rot)
        # raise ValueError(dets)

        return rot, (og_ref, ref_com, ref_axes), (og_coords, com, pax_axes)

    @classmethod
    def get_eckart_embedding_data(cls, masses, ref, coords, sel=None):
        """
        Embeds a set of coordinates in the reference frame

        :param masses:
        :type masses: np.ndarray
        :param ref:
        :type ref: CoordinateSet
        :param coords:
        :type coords: CoordinateSet
        :return:
        :rtype:
        """

        return cls.get_eckart_rotations(masses, ref, coords, sel=sel, in_paf=False)

    @classmethod
    def get_prop_eckart_transformation(cls, masses, ref, coords,
                                       sel=None,
                                       inverse=False,
                                       reset_com=False,
                                       planar_ref_tolerance=None
                                       ):
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

        multiconf = coords.multiconfig

        # if multiconf:
        #     coords = list(coords)
        # else:
        #     coords = [coords]

        ek_rot, ref_stuff, coord_stuff = cls.get_eckart_rotations(masses, ref, coords, sel=sel, in_paf=False, planar_ref_tolerance=planar_ref_tolerance)
        ref, ref_com, ref_rot = ref_stuff
        crd, crd_com, crd_rot = coord_stuff

        transforms = [None] * len(coords)
        for i, (rot, com, crd_rot) in enumerate(zip(ek_rot, crd_com, crd_rot)):
            if inverse:
                if reset_com:
                    transf = MolecularTransformation(com, crd_rot, rot.T, ref_rot.T, -ref_com)
                else:
                    transf = MolecularTransformation(crd_rot, rot.T, ref_rot.T)
            else:
                if reset_com:
                    transf = MolecularTransformation(ref_com, ref_rot, rot, crd_rot.T, -com)
                else:
                    transf = MolecularTransformation(ref_rot, rot, crd_rot.T, -com)
            transforms[i] = transf

        if not multiconf:
            transforms = transforms[0]

        return transforms

    @classmethod
    def get_eckart_embedded_coords(cls, masses,
                                   ref, coords,
                                   reset_com=False,
                                   sel=None,
                                   planar_ref_tolerance=None
                                   ):
        """
        Embeds a set of coordinates in the reference frame

        :param masses:
        :type masses: np.ndarray
        :param ref:
        :type ref: CoordinateSet
        :param coords:
        :type coords: CoordinateSet
        :return:
        :rtype:
        """

        multiconf = coords.multiconfig
        if sel is not None:
            coords = coords[..., sel, :]
            masses = masses[sel]
            ref = ref[..., sel, :]

        ek_rot, ref_stuff, coord_stuff = cls.get_eckart_rotations(masses, ref, coords, in_paf=False, planar_ref_tolerance=planar_ref_tolerance)
        ref, ref_com, ref_rot = ref_stuff
        crd, crd_com, crd_rot = coord_stuff

        # crd is in _its_ principle axis frame, so now we transform it using ek_rot
        ek_rot = np.swapaxes(ek_rot, -2, -1)
        crd = crd @ ek_rot
        # now we rotate this back to the reference frame
        crd = crd @ (ref_rot.T)[np.newaxis, :, :]

        if reset_com:
            # and then shift so the COM doesn't change
            crd = crd + ref_com[np.newaxis, np.newaxis, :]

        if not multiconf:
            crd = crd[0]

        return crd

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
        # explicitly put masses in m_e from AMU
        masses = UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass") * masses
        mT = np.sqrt(np.sum(masses))
        mvec = np.sqrt(masses)

        M = np.kron(mvec / mT, np.eye(3)).T  # translation eigenvectors
        mom_rot, ax_rot = cls.get_prop_moments_of_inertia(coords, masses)
        inv_rot_2 = np.dot(np.dot(ax_rot, np.diag(1 / np.sqrt(mom_rot))), ax_rot.T)
        com = cls.get_prop_center_of_mass(coords, masses)
        shift_crds = mvec[:, np.newaxis] * (coords - com[np.newaxis, :])
        e = nput.levi_cevita3
        R = np.tensordot(
            shift_crds,
            np.tensordot(e, inv_rot_2, axes=[0, 1]),
            axes=[1, 0]
        ).reshape((3 * n, 3))  # rotations
        freqs = np.concatenate([
            [1e-14, 1e-14, 1e-14],
            (1 / (2 * mom_rot))
            # this isn't right, I'm totally aware, but I think the frequency is supposed to be zero anyway and this
            # will be tiny
        ])
        eigs = np.concatenate([M, R], axis=1)

        # import McUtils.Plots as plt
        # print(np.round(eigs, 10))
        # plt.ArrayPlot(eigs).show()
        # raise Exception([M, R])
        return freqs, eigs

    @classmethod
    def get_prop_g_matrix(cls, masses, coords, internal_coords):
        """
        Gets the molecular g-matrix
        :param masses:
        :type masses: np.ndarray
        :param coords:
        :type coords: CoordinateSet
        :param internal_coords:
        :type internal_coords: CoordinateSet
        :return:
        :rtype:
        """
        jacobian = coords.jacobian(internal_coords.system, [1])
        if not isinstance(jacobian, np.ndarray):
            jacobian = jacobian[0]
        if coords.multiconfig:
            # strip embedding
            embedding_coords = [0, 1, 2, 4, 5, 8]
            good_coords = np.setdiff1d(np.arange(3*len(masses)), embedding_coords)
            mass_weighting = np.sqrt(masses)[:, np.newaxis, np.newaxis, np.newaxis]
            for i in range(jacobian.ndim - 4):
                jacobian = np.expand_dims(jacobian, 0)
            jacobian *= 1/mass_weighting
            jac_shapes = jacobian.shape[:-4] + (jacobian.shape[-4] * jacobian.shape[-3], jacobian.shape[-2] * jacobian.shape[-1])
            jacobian = jacobian.reshape(jac_shapes)
            jacobian = jacobian[..., :, good_coords]
        else:
            # now mass-weight
            jacobian = jacobian.reshape((len(masses), 3, len(masses), 3))
            # strip embedding
            embedding_coords = [0, 1, 2, 4, 5, 8]
            good_coords = np.setdiff1d(np.arange(3*len(masses)), embedding_coords)
            mass_weighting = np.sqrt(masses)
            jacobian *= 1/mass_weighting[:, np.newaxis, np.newaxis, np.newaxis]
            jacobian = jacobian.reshape(jacobian.shape[0]*jacobian.shape[1], jacobian.shape[2]*jacobian.shape[3])
            jacobian = jacobian[:, good_coords]
        # dot together
        return nput.vec_tensordot(jacobian, jacobian, axes=[-2, -2]).squeeze()

class BondingProperties:
    """
    The set of properties that depend only on bonding
    and that kind of format
    """

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
            (np.ones((len(bonds),)), (cons[:, 0], cons[:, 1])),
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
            w = np.where(labels == i)
            if len(w) > 0:
                w = w[0]
            return w

        return [pull(labels, i) for i in range(ngroups)]

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
        raise NotImplementedError("ran out of time to do this in May ;_;")

    @classmethod
    def get_prop_guessed_bonds(cls, mol, tol=1.05, guess_type=True):
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

class MolecularProperties:
    """
    An object whose sole purpose in life is to get molecular properties
    A property should be implemented in two parts:
        1) a classmethod called get_prop_<prop name> that takes raw inputs and uses them to compute a property
        2) a classmethod called <prop name> that extracts the property from a passed molecule

    All properties should be appropriately vectorized to work on a single configuration or a set of configurations
    """

    @classmethod
    def mass_weighted_coords(cls, mol):
        """
        Computes the moments of inertia

        :param mol:
        :type mol: AbstractMolecule
        :return:
        :rtype:
        """

        return StructuralProperties.get_prop_mass_weighted_coords(mol.coords, mol.masses)


    @classmethod
    def g_matrix(cls, mol):
        """
        :param mol:
        :type mol: AbstractMolecule
        :return:
        :rtype:
        """
        masses = mol.masses * UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
        return StructuralProperties.get_prop_g_matrix(masses, mol.coords, mol.internal_coordinates)

    @classmethod
    def center_of_mass(cls, mol):
        """
        Computes the moments of inertia

        :param mol:
        :type mol: AbstractMolecule
        :return:
        :rtype:
        """

        return StructuralProperties.get_prop_center_of_mass(mol.coords, mol.masses)

    @classmethod
    def inertia_tensor(cls, mol):
        """
        Computes the inertia tensors for the stored geometries

        :param mol:
        :type mol: AbstractMolecule
        :return:
        :rtype:
        """

        return StructuralProperties.get_prop_inertia_tensors(mol.coords, mol.masses)

    @classmethod
    def moments_of_inertia(cls, mol):
        """
        Computes the moments of inertia

        :param mol:
        :type mol: AbstractMolecule
        :return:
        :rtype:
        """

        return StructuralProperties.get_prop_moments_of_inertia(mol.coords, mol.masses)

    @classmethod
    def principle_axis_data(cls, mol, sel=None):
        """
        Generates the center of masses and inertial axes

        :param mol:
        :type mol: AbstractMolecule
        :return:
        :rtype:
        """
        coords = mol.coords
        masses = mol.masses
        if sel is not None:
            coords = coords[..., sel, :]
            masses = masses[sel]
        return StructuralProperties.get_prop_principle_axis_data(coords, masses)

    @classmethod
    def principle_axis_transformation(cls, mol, sel=None, inverse=False):
        """
        Generates the principle axis transformation for a Molecule

        :param mol:
        :type mol: AbstractMolecule
        :return:
        :rtype:
        """
        return StructuralProperties.get_prop_principle_axis_rotation(mol.coords, mol.masses, sel=sel, inverse=inverse)

    @classmethod
    def eckart_embedding_data(cls, mol, coords, sel=None, planar_ref_tolerance=None):
        """

        :param mol:
        :type mol: AbstractMolecule
        :param coords:
        :type coords:
        :param sel:
        :type sel:
        :return:
        :rtype:
        """
        masses = mol.masses
        ref = mol.coords
        coords = CoordinateSet(coords)
        return StructuralProperties.get_eckart_embedding_data(masses, ref, coords, sel=sel, planar_ref_tolerance=planar_ref_tolerance)

    @classmethod
    def eckart_transformation(cls, mol, ref_mol, sel=None, inverse=False, planar_ref_tolerance=None):
        """

        :param ref_mol: reference geometry
        :type ref_mol: AbstractMolecule
        :param mol: molecules to get Eckart embeddings for
        :type mol: AbstractMolecule
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
        return StructuralProperties.get_prop_eckart_transformation(m1, ref_mol.coords, mol.coords, sel=sel, inverse=inverse, planar_ref_tolerance=planar_ref_tolerance)

    @classmethod
    def eckart_embedded_coords(cls, mol, coords, sel=None, planar_ref_tolerance=None):
        """

        :param mol:
        :type mol: AbstractMolecule
        :param coords:
        :type coords:
        :param sel:
        :type sel:
        :return:
        :rtype:
        """
        masses = mol.masses
        ref = mol.coords
        coords = CoordinateSet(coords)
        return StructuralProperties.get_eckart_embedded_coords(masses, ref, coords, sel=sel, planar_ref_tolerance=planar_ref_tolerance)

    @classmethod
    def translation_rotation_eigenvectors(cls, mol, sel=None):
        """

        :param mol: molecules to get eigenvectors for
        :type mol: AbstractMolecule
        :param sel: coordinate selection to use when doing the rotation/translation calculations
        :type mol:
        :return:
        :rtype:
        """
        if sel is not None:
            raise NotImplementedError("Still need to add coordinate subselection")
        return StructuralProperties.get_prop_translation_rotation_eigenvectors(mol.coords, mol.masses)

    @classmethod
    def fragments(cls, mol):
        """

        :param mol:
        :type mol: AbstractMolecule
        :return:
        :rtype:
        """

        from .Molecule import Molecule

        cds = mol.coords
        bonds = mol.bonds
        ats = mol.atoms

        comps = BondingProperties.get_prop_fragments(ats, bonds)
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
    def guessed_bonds(cls, mol, tol=1.05, guess_type=True):
        """
        Guesses the bonds for the molecule by finding the ones that are less than some percentage of a single bond for that
        pair of elements

        :return:
        :rtype:
        """

        return BondingProperties.get_prop_guessed_bonds(mol, tol=tol, guess_type=guess_type)

    @classmethod
    def get_prop_chemical_formula(cls, atoms):
        """

        :param atoms:
        :type atoms: Tuple[str]
        :return:
        :rtype:
        """
        return "".join(a+str(atoms.count(a)) for a in set(atoms))
    @classmethod
    def chemical_formula(cls, mol):
        """

        :param mol:
        :type mol: AbstractMolecule
        :return:
        :rtype:
        """
        return cls.get_prop_chemical_formula(mol.atoms)

class StringFormatHandler(metaclass=abc.ABCMeta):
    """
    A base class to handle converting to/from string formats
    mostly just here to implement SDF so OpenBabel can read off that
    """
    def __init__(self, mol):
        self.mol = mol
    @abc.abstractmethod
    def get_string(self, atoms, coords, bonds, meta):
        """
        Converts the molecular info to string format

        :param atoms:
        :type atoms: tuple[str]
        :param coords:
        :type coords: CoordinateSet
        :param bonds:
        :type bonds: tuple[tuple]
        :param meta:
        :type meta: dict
        :return:
        :rtype: str
        """
        raise NotImplementedError("abstract base class")

    def convert(self):
        return self.get_string(
            self.mol.atoms,
            self.mol.coords,
            self.mol.bonds,
            self.mol.metadata
        )

    @classmethod
    @abc.abstractmethod
    def parse_string(cls, str):
        """

        :param str:
        :type str:
        :return: (atoms, coords, bonds, meta)
        :rtype: (list, np.ndarray, list, dict)
        """
        raise NotImplementedError("abstract base class")

    @classmethod
    def parse(cls, string):
        from .Molecule import Molecule
        atoms, crds, bonds, meta = cls.parse_string(string)
        return Molecule(
            atoms, crds,
            bonds=bonds,
            **meta
        )

class SDFFormatHandler(StringFormatHandler):
    misc_useless_structural_data_header = " 0     0  0  0  0  0  0999 V2000"
    program = 'Psience.Molecools'
    def __init__(self, mol):
        """
        :param mol:
        :type mol: AbstractMolecule
        :param program:
        :type program:
        :param comment:
        :type comment:
        """
        super().__init__(mol)
        self.name = mol.name
        self.comment = '' if 'comment' not in mol.metadata else mol.metadata['comment']

    def convert_header(self, comment=None):
        return "\n".join([
            self.name,
            "  " + self.program,
            " " + self.comment + ("" if comment is None else comment)
        ])

    def convert_counts_line(self, atoms, bonds):
        return "{:>3.0f}{:>3.0f} {}".format(len(atoms), len(bonds),
                                            self.misc_useless_structural_data_header)

    def convert_coordinate_block(self, atoms, coords):
        return "\n".join(
            " {0[0]:>9.5f} {0[1]:>9.5f} {0[2]:>9.5f} {1:<3} 0  0  0  0  0  0  0  0  0  0  0  0".format(
                crd,
                at
            ) for crd, at in zip(coords, atoms)
        )

    def convert_bond_block(self, bonds):
        return "\n".join(
            "{:>3.0f}{:>3.0f}{:>3.0f}  0  0  0  0".format(
                b[0] + 1,
                b[1] + 1,
                b[2] if len(b) > 2 else 1
            ) for b in bonds
        )

    def convert_metadata(self, meta):
        return "\n\n".join(
            ">  <{}>\n{}".format(k, v) for k,v in meta.items()
        )

    def get_single_structure_string(self, atoms, coords, bonds, meta):
        return "{header}\n{counts}\n{atoms}\n{bonds}\nM  END\n{meta}\n$$$$".format(
            header=self.convert_header(),
            counts=self.convert_counts_line(atoms, bonds),
            atoms=self.convert_coordinate_block(atoms, coords),
            bonds=self.convert_bond_block(bonds),
            meta=self.convert_metadata(meta)
        ).strip()

    def get_string(self, atoms, coords, bonds, meta):
        """
        Converts the molecular info to string format

        :param atoms:
        :type atoms: tuple[str]
        :param coords:
        :type coords: CoordinateSet
        :param bonds:
        :type bonds: tuple[tuple]
        :param meta:
        :type meta: dict
        :return:
        :rtype: str
        """
        if coords.multiconfig:
            return "\n".join(
                self.get_single_structure_string(
                    atoms,
                    c,
                    bonds,
                    meta
                ) for c in coords
            )
        else:
            return self.get_single_structure_string(
                atoms,
                coords,
                bonds,
                meta
            )

    @classmethod
    def parse_string(cls, str):
        """

        :param str:
        :type str:
        :return: (atoms, coords, bonds, meta)
        :rtype: (list, np.ndarray, list, dict)
        """
        raise NotImplementedError("SDF parsing not supported/use OpenBabel as interface")

class PropertyManager(metaclass=abc.ABCMeta):
    """
    A utility base class so to make it easier to have a unified way to
    handled derived properties
    """
    name = None
    def __init__(self, mol):
        """
        :param mol:
        :type mol: AbstractMolecule
        """
        self.mol=mol

    def set_molecule(self, mol):
        self.mol = mol

    @abc.abstractmethod
    def load(self):
        """
        Loads in the values

        :return:
        :rtype:
        """
        raise NotImplementedError("abstract base class")
    @abc.abstractmethod
    def update(self, val):
        """
        Updates the held values

        :param val:
        :type val:
        :return:
        :rtype:
        """
        raise NotImplementedError("abstract base class")
    @abc.abstractmethod
    def apply_transformation(self, transf):
        """
        Applies a transformation to the held values

        :param transf:
        :type transf: MolecularTransformation
        :return:
        :rtype: PropertyManager
        """
        raise NotImplementedError("abstract base class")
    @abc.abstractmethod
    def insert_atoms(self, atoms, coords, where):
        """
        Handles the insertion of new atoms into the structure

        :param atoms:
        :type atoms: tuple[str]
        :param coords:
        :type coords: CoordinateSet
        :param where:
        :type where: tuple[int]
        :return:
        :rtype:
        """
        raise NotImplementedError("abstract base class")
    @staticmethod
    def _insert_derivative_zeros(derivs, where, extra_shape=0):
        """
        Inserts zeros into derivative terms to account for
        the insertion of dummy atoms into a system

        :param derivs:
        :type derivs:
        :param where:
        :type where:
        :param extra_shape:
        :type extra_shape:
        :return:
        :rtype:
        """

        n_coords = derivs[0].shape[0]
        if isinstance(where, (int, np.integer)):
            where = [where]
        for d in derivs:
            # we reshape the derivs to the right shape
            # and insert zeros
            new_derivs = []
            for i, d in enumerate(derivs):
                if d is None:
                    new_derivs.append(d)
                elif isinstance(d, (int, float, np.integer, np.floating)):
                    new_derivs.append(d)
                else:
                    shp = d.shape
                    # print(">>", shp)
                    for j in range(d.ndim - extra_shape):
                        # print(j, d.shape[j])
                        if d.shape[j] == n_coords:
                            target_shape = (
                                    d.shape[:j] +
                                    (
                                        d.shape[j] // 3,
                                        3
                                    ) +
                                    d.shape[j + 1:]
                            )
                            new_shape = (
                                    d.shape[:j] +
                                    (
                                        (d.shape[j] // 3 + len(where)) * 3,
                                    ) +
                                    d.shape[j + 1:]
                            )
                            reshape_d = np.reshape(d, target_shape)
                            reshape_d = np.insert(reshape_d, where, 0, axis=j)
                            d = reshape_d.reshape(new_shape)
                    new_derivs.append(d)

        return new_derivs

    @abc.abstractmethod
    def delete_atoms(self, where):
        """
        Handles the deletion from the structure

        :param atoms:
        :type atoms: tuple[str]
        :param coords:
        :type coords: CoordinateSet
        :param where:
        :type where: tuple[int]
        :return:
        :rtype:
        """
        raise NotImplementedError("abstract base class")

    @staticmethod
    def _transform_derivatives(derivs, transf):
        """
        Handles the transformation of derivative tensors
        across coordinate systems

        :param derivs: derivative tensors starting at order 0
        :type derivs: tuple[np.ndarray]
        :param transf:
        :type transf: MolecularTransformation
        :return:
        :rtype: tuple[np.ndarray]
        """
        # if isinstance(derivs[0], (int, float, np.integer, np.floating)):
        #     base_shape = 0
        # else:
        #     base_shape = derivs[0].shape

        n_coords = derivs[0].shape[0]
        if transf.is_affine: # most relevant subcase
            # take inverse?
            tf = transf.transformation_function.transform #type: np.ndarray
            new_derivs = []
            for i, d in enumerate(derivs):
                if d is None:
                    new_derivs.append(d)
                elif isinstance(d, (int, float, np.integer, np.floating)):
                    new_derivs.append(d)
                else:
                    shp = d.shape
                    # print(">>", shp)
                    for j in range(d.ndim):
                        # print(j, d.shape[j])
                        if d.shape[j] == n_coords:
                            target_shape = (
                                    d.shape[:j] +
                                    (
                                        d.shape[j] // 3,
                                        3
                                    ) +
                                    d.shape[j + 1:]
                            )
                            reshape_d = np.reshape(d, target_shape)
                            reshape_d = np.moveaxis(
                                np.tensordot(tf, reshape_d, axes=[[1], [j+1]]),
                                0, j+1
                            )
                            d = reshape_d.reshape(shp)
                    new_derivs.append(d)
            return tuple(new_derivs)
        else:
            raise NotImplementedError('not sure how to transform in non-affine manner')
    def __repr__(self):
        return "{}({})".format(
            self.name if self.name is not None else type(self).__name__,
            self.mol
        )

    def copy(self):
        import copy
        return copy.copy(self)

class OpenBabelMolManager(PropertyManager):
    name="OBMol"
    def __init__(self, mol, obmol=None):
        super().__init__(mol)
        self._obmol = obmol

    def load(self):
        if self._obmol is None:
            pybel = OpenBabelInterface().pybel
            mol_string = SDFFormatHandler(self.mol).convert()
            pymol = pybel.readstring(mol_string, 'sdf')
            self._obmol = pymol.OBMol
        return self._obmol
    def update(self, val):
        self._obmol = val

    def apply_transformation(self, transf):
        """
        Applies a transformation to the held values

        :param transf:
        :type transf: MolecularTransformation
        :return:
        :rtype: PropertyManager
        """
        raise NotImplementedError("Incomplete interface")
    def insert_atoms(self, atoms, coords, where):
        """
        Handles the insertion of new atoms into the structure

        :param atoms:
        :type atoms: tuple[str]
        :param coords:
        :type coords: CoordinateSet
        :param where:
        :type where: tuple[int]
        :return:
        :rtype:
        """
        raise NotImplementedError("Incomplete interface")
    def delete_atoms(self, where):
        """
        Handles the deletion from the structure

        :param atoms:
        :type atoms: tuple[str]
        :param coords:
        :type coords: CoordinateSet
        :param where:
        :type where: tuple[int]
        :return:
        :rtype:
        """
        raise NotImplementedError("Incomplete interface")

class DipoleSurfaceManager(PropertyManager):
    name='DipoleSurface'
    def __init__(self, mol, surface=None, derivatives=None):
        super().__init__(mol)
        self._surf = surface
        if isinstance(derivatives, dict):
            self._derivs = derivatives['numerical']
            self._analytic_derivatives = derivatives['analytic']
        else:
            self._derivs = derivatives
            self._analytic_derivatives = None

    @property
    def surface(self):
        if self._surf is None:
            self._surf = self.load_dipole_surface()
        return self._surf
    @property
    def numerical_derivatives(self):
        if self._numerical_derivs is None and self._derivs is None:
            derivatives = self.load_dipole_derivatives()
            if isinstance(derivatives, dict):
                self._numerical_derivs = derivatives['numerical']
                self._derivs = derivatives['analytic']
            else:
                self._numerical_derivs = None
                self._derivs = derivatives
        return self._numerical_derivs
    @property
    def derivatives(self):
        if self._derivs is None:
            derivatives = self.load_dipole_derivatives()
            if isinstance(derivatives, dict):
                self._numerical_derivs = derivatives['numerical']
                self._derivs = derivatives['analytic']
            else:
                self._numerical_derivs = None
                self._derivs = derivatives
        return self._derivs
    @derivatives.setter
    def derivatives(self, derivatives):
        if isinstance(derivatives, dict):
            self._numerical_derivs = derivatives['numerical']
            self._derivs = derivatives['analytic']
        else:
            self._numerical_derivs = None
            self._derivs = derivatives

    def load(self):
        if self._surf is not None:
            return self.surface
        else:
            return self.derivatives
    def update(self, val):
        """
        Updates the held values

        :param val:
        :type val:
        :return:
        :rtype:
        """
        raise NotImplementedError("Incomplete interface")

    def load_dipole_surface(self):
        raise NotImplementedError("haven't needed general dipole surfaces yet")

    def _load_gaussian_fchk_dipoles(self, file):
        try:
            keys = ['DipoleMoment', 'DipoleDerivatives', 'DipoleHigherDerivatives', 'DipoleNumDerivatives']
            with GaussianFChkReader(file) as gr:
                parse = gr.parse(keys)

            mom, grad, high = tuple(parse[k] for k in keys[:3])
            grad = grad.array
            seconds = high.second_deriv_array
            thirds = high.third_deriv_array
            num_derivs = parse[keys[3]]
            num_grad = num_derivs.first_derivatives
            num_secs = num_derivs.second_derivatives
        except GaussianFChkReaderException:
            keys = ['DipoleMoment', 'DipoleDerivatives']
            with GaussianFChkReader(file) as gr:
                parse = gr.parse(keys)

            mom, grad = tuple(parse[k] for k in keys[:2])
            grad = grad.array
            seconds = thirds = None
            num_grad = num_secs = None

        if seconds is not None:
            amu_conv = UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
            seconds = seconds / np.sqrt(amu_conv)
        if thirds is not None:
            amu_conv = UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
            thirds = thirds / amu_conv

        return {
            "analytic": (mom, grad, seconds, thirds),
            "numerical": (mom, num_grad, num_secs)
        }

    def load_dipole_derivatives(self, file=None):
        """
        Loads dipole derivatives from a file (or from `source_file` if set)

        :param file:
        :type file:
        :return:
        :rtype:
        """

        if file is None:
            file = self.mol.source_file
        if file is None:
            return None

        path, ext = os.path.splitext(file)
        ext = ext.lower()

        if ext == ".fchk":
            return self._load_gaussian_fchk_dipoles(file)
        elif ext == ".log":
            raise NotImplementedError("{}: support for loading dipole derivatives from {} files not there yet".format(
                type(self).__name__,
                ext
            ))
        else:
            raise NotImplementedError("{}: support for loading dipole derivatives from {} files not there yet".format(
                type(self).__name__,
                ext
            ))

    def apply_transformation(self, transf):
        new = self.copy()
        if new._surf is not None:
            new._surf = new._surf.transform(transf)
        if new._derivs is not None:
            tf = transf.transformation_function.transform
            base_derivs = (
                    (new._derivs[0],) +
                    tuple(new._transform_derivatives(new._derivs[1:], transf))
            )
            new._derivs = tuple(np.tensordot(d, tf, axes=[-1, -1]) for d in base_derivs)
        return new

    def insert_atoms(self, atoms, coords, where):
        """
        Handles the insertion of new atoms into the structure

        :param atoms:
        :type atoms: tuple[str]
        :param coords:
        :type coords: CoordinateSet
        :param where:
        :type where: tuple[int]
        :return:
        :rtype:
        """
        if all(a == "X" for a in atoms):
            new = self.copy()
            dip_mom = self.derivatives[0]
            new.derivatives = (
                    (dip_mom,) +
                    tuple(self._insert_derivative_zeros(self.derivatives[1:], where, extra_shape=1))
            )
            return new
        else:
            raise NotImplementedError("don't know how to insert non-dummy atoms")
    def delete_atoms(self, where):
        """
        Handles the deletion from the structure

        :param atoms:
        :type atoms: tuple[str]
        :param coords:
        :type coords: CoordinateSet
        :param where:
        :type where: tuple[int]
        :return:
        :rtype:
        """
        raise NotImplementedError("Incomplete interface")

class PotentialSurfaceManager(PropertyManager):
    name="PotentialSurface"
    def __init__(self, mol, surface=None, derivatives=None):
        super().__init__(mol)
        self._surf = surface
        self._derivs = derivatives

    @property
    def surface(self):
        if self._surf is None:
            self._surf = self.load_potential_surface()
        return self._surf

    @property
    def derivatives(self):
        if self._derivs is None:
            self._derivs = self.load_potential_derivatives()
        return self._derivs
    @derivatives.setter
    def derivatives(self, v):
        self._derivs = v

    @property
    def force_constants(self):
        return self.derivatives[1]

    def load_potential_derivatives(self, file=None):
        """
        Loads potential derivatives from a file (or from `source_file` if set)

        :param file:
        :type file:
        :return:
        :rtype:
        """

        if file is None:
            file = self.mol.source_file

        path, ext = os.path.splitext(file)
        ext = ext.lower()

        if ext == ".fchk":
            keys= ['Gradient', 'ForceConstants', 'ForceDerivatives']
            with GaussianFChkReader(file) as gr:
                parse = gr.parse(keys, default=None)

            seconds = parse["ForceConstants"].array
            if parse["ForceDerivatives"] is not None:
                thirds = parse["ForceDerivatives"].third_deriv_array
                fourths = parse["ForceDerivatives"].fourth_deriv_array
                if isinstance(fourths, nput.SparseArray):
                    fourths = fourths.asarray()

                amu_conv = UnitsData.convert("AtomicMassUnits", "AtomicUnitOfMass")
                thirds = thirds / np.sqrt(amu_conv)
                fourths = fourths / amu_conv
            else:
                thirds = fourths = None

            return (parse["Gradient"], seconds, thirds, fourths)
        elif ext == ".log":
            raise NotImplementedError("{}: support for loading force constants from {} files not there yet".format(
                type(self).__name__,
                ext
            ))
        else:
            raise NotImplementedError("{}: support for loading force constants from {} files not there yet".format(
                type(self).__name__,
                ext
            ))
    def load_potential_surface(self, coordinates):
        from ..Data import PotentialSurface

        ang2bohr = UnitsData.convert("Angstroms", "BohrRadius")

        try:
            iter(coordinates)
        except TypeError:
            coord_transf = coordinates
        else:
            fns = []
            for ctup in coordinates:
                if len(ctup) == 2:
                    fns.append(lambda c, i=ctup[0], j=ctup[1]: ang2bohr*nput.pts_norms(c[:, i], c[:, j]))
                elif len(ctup) == 3:
                    fns.append(lambda c, i=ctup[0], j=ctup[1], k=ctup[2]: nput.pts_angles(c[:, i], c[:, j], c[:, k])[0])
                elif len(ctup) == 4:
                    fns.append(lambda c, i=ctup[0], j=ctup[1], k=ctup[2], l=ctup[3]: nput.pts_dihedrals(c[:, i], c[:, j], c[:, k],  c[:, l])[0])
                else:
                    raise ValueError("don't know how to interpret coordinate spec '{}'".format(ctup))
            def coord_transf(crds, fns=fns):
                return np.array([f(crds) for f in fns]).T

        surf = PotentialSurface.from_log_file(
            self.mol.source_file,
            coord_transf
        )
        return surf

    def load(self, coordinates=None):
        if coordinates is not None and self._surf is None:
            self._surf = self.load_potential_surface(coordinates)
        if self._surf is not None:
            return self.surface
        else:
            return self.derivatives
    def update(self, val):
        """
        Updates the held values

        :param val:
        :type val:
        :return:
        :rtype:
        """
        raise NotImplementedError("Incomplete interface")

    def apply_transformation(self, transf):
        new = self.copy()
        if new._surf is not None:
            new._surf = new._surf.transform(transf)
        if new._derivs is not None:
            # raise Exception(new._derivs)
            # print([x.shape for x in new._derivs])
            new._derivs = new._transform_derivatives(new._derivs, transf)
        return new

    def insert_atoms(self, atoms, coords, where):
        """
        Handles the insertion of new atoms into the structure

        :param atoms:
        :type atoms: tuple[str]
        :param coords:
        :type coords: CoordinateSet
        :param where:
        :type where: tuple[int]
        :return:
        :rtype:
        """
        if all(a == "X" for a in atoms):
            new = self.copy()
            new.derivatives = self._insert_derivative_zeros(self.derivatives, where, extra_shape=0)
            return new
        else:
            raise NotImplementedError("don't know how to insert non-dummy atoms")
    def delete_atoms(self, where):
        """
        Handles the deletion from the structure

        :param atoms:
        :type atoms: tuple[str]
        :param coords:
        :type coords: CoordinateSet
        :param where:
        :type where: tuple[int]
        :return:
        :rtype:
        """
        raise NotImplementedError("Incomplete interface")

class NormalModesManager(PropertyManager):
    def __init__(self, mol, normal_modes=None):
        super().__init__(mol)
        self._modes = normal_modes

    def set_molecule(self, mol):
        super().set_molecule(mol)
        if self._modes is not None:
            self.modes.molecule = mol

    @property
    def modes(self):
        """

        :return:
        :rtype: MolecularVibrations
        """
        if self._modes is None:
            self._modes = self.get_normal_modes()
        return self._modes
    @modes.setter
    def modes(self, modes):
        if not isinstance(modes, MolecularVibrations):
            modes = self.construct_normal_modes(modes)
        if not isinstance(modes, MolecularVibrations):
            raise TypeError("`modes` must be {}".format(
                MolecularVibrations.__name__
            ))
        self._modes = modes

    def construct_normal_modes(self, modes):
        if isinstance(modes, dict):
            modes = modes.copy()
            coeffs = modes['matrix']
            del modes['matrix']
            modes = MolecularNormalModes(self.mol, coeffs, **modes)
        elif isinstance(modes, np.ndarray):
            modes = MolecularNormalModes(self.mol, modes)
        elif not isinstance(modes, MolecularNormalModes):
            raise ValueError("don't know how to construct `{}` from {}".format(
                MolecularNormalModes.__name__,
                modes
            ))
        if isinstance(modes, MolecularNormalModes):
            modes = MolecularVibrations(self.mol, modes)

        return modes

    def load(self):
        return self._modes
    def update(self, modes):
        """

        :return:
        :rtype:
        """
        if not isinstance(modes, MolecularVibrations):
            raise TypeError("{}.{}: '{}' is expected to be a MolecularVibrations object".format(
                type(self).__name__,
                'normal_modes',
                modes
            ))
        self._modes = modes
    #TODO: need to be careful about transition states...
    recalc_normal_mode_tolerance = 1.0e-8
    def load_normal_modes(self, file=None, rephase=True, recalculate=False):
        """
        Loads potential derivatives from a file (or from `source_file` if set)

        :param file:
        :type file:
        :param rephase: whether to rephase FChk normal modes or not
        :type rephase: bool
        :return:
        :rtype:
        """
        from .Vibrations import MolecularNormalModes

        if file is None:
            file = self.mol.source_file
        path, ext = os.path.splitext(file)
        ext = ext.lower()

        if ext == ".fchk":
            with GaussianFChkReader(file) as gr:
                parse = gr.parse(
                    ['Real atomic weights', 'VibrationalModes', 'ForceConstants', 'VibrationalData']
                )

            modes = parse["VibrationalModes"]
            freqs = parse["VibrationalData"]["Frequencies"] * UnitsData.convert("Wavenumbers", "Hartrees")
            amu_conv = UnitsData.convert("AtomicMassUnits", "ElectronMass")
            masses = parse['Real atomic weights'] # * amu_conv

            fcs = self.mol.potential_surface.force_constants

            sqrt_freqs = np.sign(freqs) * np.sqrt(np.abs(freqs))

            # I do a bunch of messing with stuff here for like 0 pay off
            # over the old simple method...not sure what un-dimensioning
            # is really happening here...
            mass_vec = np.broadcast_to(masses[:, np.newaxis], (len(masses), 3)).flatten() * amu_conv
            modes_conv = modes * np.sqrt(mass_vec[np.newaxis, :])
            modes_conv = modes_conv / (
                    np.linalg.norm(modes_conv, axis=1)[:, np.newaxis] *
                    sqrt_freqs[:, np.newaxis] * np.sqrt(mass_vec[np.newaxis, :])
            )

            #
            internal_F2 = np.dot(np.dot(modes_conv, fcs), modes_conv.T)

            final_scale = sqrt_freqs**(3) / np.diag(internal_F2)
            modes_conv *= final_scale[:, np.newaxis]
            modes = modes_conv

            modes = modes.T

            modes = MolecularNormalModes(self.mol, modes, inverse=modes.T, freqs=freqs)

            if rephase:
                try:
                    phases = self.get_fchk_normal_mode_rephasing(modes)
                except NotImplementedError:
                    pass
                else:
                    if phases is not None:
                        modes = modes.rescale(phases)

            if recalculate:
                fcs = self.mol.potential_surface.force_constants
                new_modes = MolecularNormalModes.from_force_constants(self.mol, fcs, self.mol.atoms)
                new_old_phases = np.dot(modes.matrix.T, new_modes.matrix)
                old_old_phases = np.dot(modes.matrix.T, modes.matrix)
                new_new_phases = np.dot(new_modes.matrix.T, new_modes.matrix)

                if not np.allclose(
                        np.diag(new_new_phases),
                        np.diag(old_old_phases),
                    atol=self.recalc_normal_mode_tolerance
                ):
                    raise ValueError("normal modes from Gaussian are normalized differently than normal modes from diagonalizing G...?")

                phases = np.sign(np.diag(new_old_phases) / np.diag(old_old_phases))
                modes = new_modes.rescale(phases)

            return modes
        elif ext == ".log":
            raise NotImplementedError("{}: support for loading normal modes from {} files not there yet".format(
                type(self).__name__,
                ext
            ))
        else:
            raise NotImplementedError("{}: support for loading normal modes from {} files not there yet".format(
                type(self).__name__,
                ext
            ))
    def get_normal_modes(self, **kwargs):
        """
        Loads normal modes from file or calculates
        from force constants

        :param kwargs:
        :type kwargs:
        :return:
        :rtype:
        """

        if self.mol.source_file is not None:
            vibs = MolecularVibrations(self.mol, self.load_normal_modes())
        else:
            fcs = self.mol.potential_surface.force_constants
            vibs = MolecularVibrations(self.mol,
                                       MolecularNormalModes.from_force_constants(self.mol, fcs, self.mol.atoms, **kwargs)
                                       )

        return vibs
    def get_fchk_normal_mode_rephasing(self, modes=None):
        """
        Returns the necessary rephasing to make the numerical dipole derivatives
        agree with the analytic dipole derivatives as pulled from a Gaussian FChk file
        :return:
        :rtype:
        """

        d1_analytic = self.mol.dipole_surface.derivatives
        if d1_analytic is None:
            return None
        else:
            d1_analytic = d1_analytic[1]
            if d1_analytic is None:
                return None
        d1_numerical = self.mol.dipole_surface.numerical_derivatives
        if d1_numerical is None:
            return None
        else:
            d1_numerical = d1_numerical[1]
            if d1_numerical is None:
                return None

        # if not isinstance(d1_analytic, np.ndarray):
        #     d1_analytic = d1_analytic.array
        # if not isinstance(d1_numerical, np.ndarray):
        #     d1_numerical = d1_numerical.first_derivatives

        # raise Exception(d1_analytic.shape, d1_numerical.shape

        if modes is None:
            modes = self.modes.basis

        mode_basis = modes.matrix
        rot_analytic = np.dot(d1_analytic.T, mode_basis)

        rot_analytic = np.nan_to_num(np.array([x/np.linalg.norm(x) for x in rot_analytic]), nan=0.0)
        d1_numerical = np.nan_to_num(np.array([x/np.linalg.norm(x) for x in d1_numerical.T]), nan=0.0)

        # could force some orientation relative to PAX frame?
        # dot each into each other
        phases = np.sign(np.diag(np.dot(d1_numerical.T, rot_analytic)))

        return phases

    def apply_transformation(self, transf):
        # self.modes # load in for some reason?
        new = self.copy()

        modes = new.modes
        modes = modes.embed(transf)
        new.modes = modes

        return new

    def insert_atoms(self, atoms, coords, where):
        """
        Handles the insertion of new atoms into the structure

        :param atoms:
        :type atoms: tuple[str]
        :param coords:
        :type coords: CoordinateSet
        :param where:
        :type where: tuple[int]
        :return:
        :rtype:
        """
        import copy
        new = self.copy()

        modes = new.modes
        new_basis = modes.basis.insert(0, where)
        modes = copy.copy(modes)
        modes.basis = new_basis
        new.modes = modes

        return new
    def delete_atoms(self, where):
        """
        Handles the deletion from the structure

        :param atoms:
        :type atoms: tuple[str]
        :param coords:
        :type coords: CoordinateSet
        :param where:
        :type where: tuple[int]
        :return:
        :rtype:
        """
        raise NotImplementedError("Incomplete interface")