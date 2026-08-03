"""
Provides a `MoleculeBuilder` class that constructs molecules from functional groups or
fragments.

This gathers the construction-oriented routines that used to live directly on `Molecule`:
computing a local attachment frame for a fragment (`fragment_embedding`), attaching a group of
atoms onto a scaffold at that frame (`attach_functional_group`), and assembling a whole
molecule from a scaffold plus a set of replacement fragments driven by a templated-SMILES join
(`from_fragments`), along with the stereochemistry/geometry helpers those depend on.

Every entry point is a `classmethod`, so `MoleculeBuilder` needs no constructor and is never
instantiated -- it's a namespace of build operations. The routines that used to operate on
`self` now take the scaffold `Molecule` as an explicit first argument, and the routines that
used to construct `Molecule` objects via `cls(...)`/`cls.from_string(...)` take the molecule
class as `molecule_type` (defaulting to `Psience.Molecools.Molecule.Molecule`).
"""

from __future__ import annotations

import functools
import numpy as np

from McUtils.Data import AtomData, UnitsData, BondData
import McUtils.Numputils as nput
import McUtils.Devutils as dev
import McUtils.Iterators as itut
from McUtils.ExternalPrograms import build_templated_smiles, parse_smiles_and_atom_map

__all__ = [
    "MoleculeBuilder"
]

__reload_hook__ = ['.Molecule']


class MoleculeBuilder:
    """
    Namespace of classmethods for constructing molecules from functional groups or fragments.

    Not instantiated: all functionality is exposed through classmethods, so no constructor is
    needed. The scaffold molecule (formerly `self`) is passed in explicitly, and the molecule
    class used to build new objects is passed as `molecule_type` (defaulting to the standard
    `Molecule`).
    """

    @classmethod
    def _default_molecule_type(cls):
        """
        Resolve the default `Molecule` class to build with, imported lazily to avoid a circular
        import at module load.

        :return: the standard `Molecule` class
        :rtype: type
        """
        from .Molecule import Molecule
        return Molecule

    #region Fragment frame
    @classmethod
    def fragment_embedding(cls, mol, fragment_indices,
                           ref=None,
                           return_axes=False,
                           order=1,
                           view_inds=(1, 2),
                           excluded=None,
                           use_moments=False):
        """
        Compute a local coordinate frame (origin, offset vector, and an up-vector or full axis
        set) anchored at a fragment of `mol`, used as the reference frame for attaching or
        orienting substituents; falls back to center-of-mass/principal-axis reference points
        (encoded as indices `-1`/`-2`/`-3`) when the fragment doesn't have enough atoms of its
        own to define a frame.

        :param mol: the scaffold molecule the fragment lives in
        :type mol: Molecule
        :param fragment_indices: the atom index (or indices) defining the fragment to embed
        :type fragment_indices: int | Iterable[int]
        :param ref: reference atom(s) (outside the fragment) used to anchor the origin/frame;
            computed from the local neighborhood if not given
        :type ref: Iterable[int] | None
        :param return_axes: whether to return a full 3x3 axis frame instead of just an up-vector
        :type return_axes: bool
        :param view_inds: which two fragment-atom positions define the "view" direction used to
            build the axis frame
        :type view_inds: tuple[int, int]
        :param use_moments: whether to derive the up-vector from the fragment's moments of
            inertia rather than from its first three atom positions
        :type use_moments: bool
        :return: `(origin, offset, up_or_axes)` -- the reference origin point, the offset from
            origin to the fragment's first atom, and either an up-vector or (if `return_axes`) a
            full axis frame
        :rtype: tuple
        """
        if nput.is_int(fragment_indices):
            fragment_indices = [fragment_indices]
        if len(fragment_indices) < order:
            raise ValueError(f"need at least {len(fragment_indices)} indices or embedding of order {order}")
        if excluded is None:
            excluded = ()
        if ref is None:
            ref = list(itut.flatten(zip(*(
                (m for m in mol.neighborhood(fragment_indices[i], size=1, heavy_only=True) if m not in excluded)
                for i in range(order)
            ))))
        if len(fragment_indices) < 3:
            if len(fragment_indices) + len(ref) < 3:
                extra = [m for m in mol.neighborhood(fragment_indices[0], size=2, heavy_only=True) if m not in excluded]
                ref = np.asanyarray(ref, dtype=int)
                extra = np.asanyarray([r for r in extra if r not in ref], dtype=int)
                ref = np.concatenate([ref, [r for r in extra if r not in ref]])
            ref = [r for r in ref if r not in fragment_indices]

            fragment_indices = np.concatenate([fragment_indices, ref])
            if len(fragment_indices) < 3: # these refer to COM and principle axis positions
                fragment_indices = np.concatenate([fragment_indices, [-1, -2, -3]])

        ref = np.asanyarray(ref, dtype=int)[:order]
        r_coords = mol.coords[ref,]
        for i,r0 in enumerate(ref):
            if r0 == -1: r_coords[i] = mol.center_of_mass
            if r0 == -2: r_coords[i] = mol.center_of_mass + mol.inertial_axes[:, 2]
            if r0 == -3: r_coords[i] = mol.center_of_mass + mol.inertial_axes[:, 0]

        f_coords = mol.coords[fragment_indices,]
        com_f = fragment_indices == -1
        if np.any(com_f):
            f_coords[com_f] = mol.center_of_mass
        paxc_f = fragment_indices == -2
        if np.any(paxc_f):
            f_coords[paxc_f] = mol.center_of_mass + mol.inertial_axes[:, 2]
        paxa_f = fragment_indices == -3
        if np.any(paxa_f):
            f_coords[paxa_f] = mol.center_of_mass + mol.inertial_axes[:, 0]

        origin = np.average(r_coords, axis=0)
        offset = np.average(f_coords[:order], axis=0) - origin

        if use_moments:
            _, axes = nput.moments_of_inertia(f_coords, np.asanyarray(mol.masses)[fragment_indices,])
            up = axes[:, 2]
        else:
            # need origin, displacement vector, and up-vector
            fragment_indices = np.asanyarray(fragment_indices, dtype=int)[:3]
            f_coords = mol.coords[fragment_indices,]

            up = nput.pts_normals(*f_coords, normalize=True)
            if return_axes:
                axes = nput.view_matrix(
                    up_vector=up,
                    view_vector=f_coords[view_inds[0]] - f_coords[view_inds[1]]
                )
            else:
                axes = np.eye(3)

        if return_axes:
            return origin, offset, axes
        else:
            return origin, offset, up
    #endregion

    #region Geometry / stereo helpers
    @staticmethod
    def _vdw_clash_metric(frag_coords, other_coords,
                          frag_atoms=None, other_atoms=None,
                          frag_radii=None, other_radii=None):
        """
        Score a candidate placement of a fragment against the rest of a scaffold by the (negated)
        sum of squared van-der-Waals overlaps; higher (less negative) means fewer/softer clashes.

        :param frag_coords: coordinates of the fragment atoms
        :type frag_coords: np.ndarray
        :param other_coords: coordinates of the surviving scaffold atoms
        :type other_coords: np.ndarray
        :param frag_atoms: element symbols of the fragment atoms, used to look up vdW radii
        :type frag_atoms: Iterable[str] | None
        :param other_atoms: element symbols of the scaffold atoms, used to look up vdW radii
        :type other_atoms: Iterable[str] | None
        :param frag_radii: explicit fragment vdW radii (overrides `frag_atoms`)
        :type frag_radii: np.ndarray | None
        :param other_radii: explicit scaffold vdW radii (overrides `other_atoms`)
        :type other_radii: np.ndarray | None
        :return: the clash score (higher is better)
        :rtype: float
        """
        if frag_radii is None:
            if frag_atoms is not None:
                frag_radii = np.array([AtomData[a, "VanDerWaalsRadius"] for a in frag_atoms])
            else:
                frag_radii = np.array([1.5])
        if other_radii is None:
            if other_atoms is not None:
                other_radii = np.array([AtomData[a, "VanDerWaalsRadius"] for a in other_atoms])
            else:
                other_radii = np.array([1.5])
        radii_sum = (frag_radii[:, np.newaxis] + other_radii[np.newaxis, :]) * UnitsData.convert("Angstroms", "BohrRadius")

        diffs = frag_coords[:, np.newaxis, :] - other_coords[np.newaxis, :, :]
        dists = np.linalg.norm(diffs, axis=-1)
        overlap = np.maximum(radii_sum - dists, 0.0)  # only penalize actual clashes
        return -np.sum(overlap ** 2)  # higher (less negative) is better

    @staticmethod
    def _angle_diff(a, b):
        """
        Compute the smallest absolute difference between two angles, accounting for wraparound at
        `2*pi`.

        :param a: the first angle (radians)
        :type a: float
        :param b: the second angle (radians)
        :type b: float
        :return: the wrapped angular difference in `[0, pi]`
        :rtype: float
        """
        d = (a - b) % (2*np.pi)
        return min(d, (2*np.pi) - d)

    @classmethod
    def resolve_stereo_hydrogen(cls, atoms, coords, bonds, stereo_pos, stereos, ref_exclude=()):
        """
        atoms, coords, bonds : the geometry the stereocenter currently lives in
                                (either `mol` or an unattached fragment -- both
                                are static at this point, so their existing
                                torsions are ground truth)
        stereo_pos            : local index of the atom being functionalized
                                 (i.e. one of the two atoms of the double bond)
        stereos                : {(i, j): 'cis'/'trans'}
        ref_exclude             : local indices to exclude when hunting for the
                                  retained reference substituent on the far atom
                                  (pass [target_pos] etc. as needed)

        Returns the local index of the hydrogen to use as the attachment site
        (i.e. the one to pass as `target_fragment` or `group_site`), or None if
        there's no ambiguity to resolve (caller should fall back to its default
        H-picking logic).
        """
        h_candidates = [
            (b[1] if b[0] == stereo_pos else b[0])
            for b in bonds
            if stereo_pos in b[:2]
               and atoms[b[1] if b[0] == stereo_pos else b[0]] in ('H', 'D', 'T')
        ]
        if len(h_candidates) < 2:
            return None  # nothing to disambiguate

        for bb in bonds:
            partner_pos = None
            if bb[0] == stereo_pos:
                partner_pos = bb[1]
            elif bb[1] == stereo_pos:
                partner_pos = bb[0]
            if partner_pos is not None:
                relation = stereos.get((bb[0], bb[1]))
                if relation is None:
                    relation = stereos.get((bb[1], bb[0]))
                if relation is not None:
                    break
        else:
            return None

        # the retained substituent on the far side anchors what cis/trans means
        ref_pos = next(
            (
                (b[1] if b[0] == partner_pos else b[0])
                for b in bonds
                if partner_pos in b[:2]
                   and (b[1] if b[0] == partner_pos else b[0]) not in (stereo_pos, *ref_exclude)
                   and atoms[b[1] if b[0] == partner_pos else b[0]] not in ('H', 'D', 'T')
            ),
            None
        )
        if ref_pos is None:
            return None  # no usable reference substituent found

        target_torsion = 0 if relation == 'cis' else np.pi
        return min(
            h_candidates,
            key=lambda h: cls._angle_diff(
                nput.pts_dihedrals(coords[h], coords[stereo_pos], coords[partner_pos], coords[ref_pos]),
                target_torsion
            )
        )
    #endregion

    #region Functional-group attachment
    @staticmethod
    def _is_multisite_spec(target_fragment, group_site):
        """
        Decide whether the call describes multiple bond sites. Multi-site mode is selected when
        `group_site` is given as a (non-string) sequence of site atoms; a scalar `group_site`
        (or `None`) keeps the legacy single-site behavior.

        :param target_fragment: the target-fragment specification passed to `attach_functional_group`
        :type target_fragment: Any
        :param group_site: the `group_site` specification passed to `attach_functional_group`
        :type group_site: Any
        :return: whether to use the multi-site path
        :rtype: bool
        """
        if group_site is None:
            return False
        if nput.is_int(group_site):
            return False
        # a bare list/array/tuple of ints => one attachment atom per bond site
        try:
            n = len(group_site)
        except TypeError:
            return False
        return n > 0

    @classmethod
    def _group_bonds(cls, atoms, coords, bonds, molecule_type):
        """
        Resolve the group's own bond list, guessing it from `atoms`/`coords` (via `molecule_type`)
        when it is `None` or `'recompute'`.

        :param atoms: the group's element symbols
        :type atoms: Iterable[str]
        :param coords: the group's coordinates
        :type coords: np.ndarray
        :param bonds: the group's bonds; `'recompute'`/`None` to guess, or an explicit list
        :type bonds: str | list | None
        :param molecule_type: the molecule class used to guess bonds
        :type molecule_type: type
        :return: the resolved bond list
        :rtype: list
        """
        if bonds is None or dev.str_is(bonds, 'recompute'):
            return list(molecule_type(list(atoms), np.asanyarray(coords, dtype=float)).bonds)
        return list(bonds)

    @classmethod
    def _attachment_atom_for_placeholder(cls, group_bonds, placeholder):
        """
        Find the heavy attachment atom for a placeholder `group_site`: the atom the placeholder
        is bonded to (its unique neighbor in the group). This is the atom that will bond to the
        scaffold once the placeholder is removed.

        :param group_bonds: the group's bond list
        :type group_bonds: Iterable
        :param placeholder: the local index of the placeholder atom being replaced
        :type placeholder: int
        :return: the local index of the attachment atom
        :rtype: int
        :raises ValueError: if the placeholder has no bonded neighbor to attach through
        """
        for b in group_bonds:
            if b[0] == placeholder:
                return b[1]
            if b[1] == placeholder:
                return b[0]
        raise ValueError(
            f"group_site placeholder {placeholder} has no bonded neighbor to attach through"
        )

    @classmethod
    def attach_functional_group(cls,
                                mol,
                                target_fragment,
                                atoms,
                                new_coords,
                                bonds='recompute',
                                ref=None,
                                masses=None,
                                distance='auto',
                                angle=0,
                                dihedral='auto',
                                bond_sites=None,
                                dihedral_search_steps=36,
                                dihedral_distance_metric=None,
                                embedding='auto',
                                bond_order=None,
                                use_absolue_posititions=False,
                                group_site=None,
                                molecule_type=None
                                ):
        """
        Build a copy of `mol` with a new group of atoms (`atoms`/`new_coords`) attached at
        `target_fragment`, positioning and orienting the new group using the fragment's local
        reference frame (bond distance/angle/dihedral, or an explicit embedding), and splicing
        the corresponding bonds into the result; supports designating a `group_site` atom within
        the new group as its attachment point, in which case the method recurses after
        re-deriving the embedding/bonds relative to that site.

        **Multiple bond sites.** More than one bond can be formed at once by passing a list of
        binding sites: `target_fragment` becomes a list (one target fragment per site), `ref`
        and `bond_order` may be given as matching per-site lists, and `group_site` becomes a
        list of the group atoms to replace at each corresponding site. As in the single-site
        `group_site` path, each `group_site` atom is a *placeholder that is removed* -- typically
        a hydrogen -- and the heavy atom it was bonded to becomes the real attachment atom that
        bonds to the scaffold. This lets a bidentate (or polydentate) ligand be attached to an
        existing scaffold at several points at once, closing a ring or forming a linker. The
        group is oriented so its attachment atoms line up with the scaffold's per-site bond
        directions, with the overall frame anchored on the reference atoms of the *first* binding
        site by default (`embedding='auto'`); an explicit `embedding` overrides this.

        :param mol: the scaffold molecule the group is attached to
        :type mol: Molecule
        :param target_fragment: the atom(s) of `mol` the new group attaches to/replaces; a list
            of target fragments (one per site) selects multi-site mode
        :type target_fragment: int | Iterable[int] | Iterable[Iterable[int]]
        :param atoms: the element symbols of the atoms in the new group
        :type atoms: Iterable[str]
        :param new_coords: the (local) coordinates of the new group's atoms
        :type new_coords: np.ndarray
        :param bonds: bonds within the new group; `'recompute'` to guess them fresh, `None` to
            reuse `mol.bonds` remapped, or an explicit bond list
        :type bonds: str | list | None
        :param ref: reference atom(s) used to anchor the attachment frame; computed automatically
            if not given. In multi-site mode this may be a per-site list of reference-atom lists
        :type ref: Iterable[int] | Iterable[Iterable[int]] | None
        :param masses: masses for the new group's atoms; looked up from `atoms` if not given
        :type masses: np.ndarray | None
        :param distance: the bond distance to place the new group at; `'auto'` to look it up from
            `BondData`, or `None`/a number
        :type distance: str | float | None
        :param angle: rotation angle (about the up-vector) to apply to the new group
        :type angle: float
        :param dihedral: rotation angle (about the offset axis) to apply to the new group;
            `'auto'` to instead scan `dihedral_search_steps` evenly-spaced angles and keep
            whichever maximizes `dihedral_distance_metric` between the new group and the rest of
            the scaffold
        :type dihedral: float | str
        :param dihedral_search_steps: number of evenly-spaced angles (over 360°) to try when
            `dihedral='auto'`
        :type dihedral_search_steps: int
        :param dihedral_distance_metric: `(frag_coords, other_coords) -> float` scoring function
            used when `dihedral='auto'`; higher is better. Defaults to the average pairwise
            distance between the new group's atoms and the surviving scaffold atoms
        :type dihedral_distance_metric: Callable | None
        :param embedding: the reference orientation for the new group; `'auto'` to derive it from
            moments of inertia (single site) or from the per-site bond geometry anchored on the
            first site (multi-site), or an explicit `(origin, axes)`/axes specification
        :type embedding: str | tuple | np.ndarray | None
        :param bond_order: the bond order connecting the new group to the target fragment;
            defaults to `1` (or inferred when `group_site` is used). In multi-site mode may be a
            per-site list
        :type bond_order: float | Iterable[float] | None
        :param use_absolue_posititions: whether `new_coords` should be used as absolute
            coordinates rather than being repositioned relative to the fragment frame
        :type use_absolue_posititions: bool
        :param group_site: index (within `atoms`/`new_coords`) of the placeholder atom that
            marks the attachment point; the placeholder is removed and its heavy neighbor becomes
            the atom that bonds to the scaffold. A single int uses the single-site path; a list
            of ints selects multi-site mode (each placeholder removed and its neighbor bonded at
            the corresponding site), e.g. two hydrogens replaced when a bidentate ligand closes a
            ring or forms a linker
        :type group_site: int | Iterable[int] | None
        :param molecule_type: the molecule class used to build any intermediate molecules;
            defaults to the standard `Molecule`
        :type molecule_type: type | None
        :return: the molecule with the new group attached
        :rtype: Molecule
        """
        if molecule_type is None:
            molecule_type = cls._default_molecule_type()

        if group_site is not None:
            if nput.is_int(group_site): group_site = [(0, group_site)]
            group_site = [(i, g) if nput.is_int(g) else g for i, g in enumerate(group_site)]
            group_roots = [i for i,g in group_site]
            group_inds = [g for i,g in group_site]
            if dev.str_is(embedding, 'auto'):
                if dev.str_is(bonds, 'recompute') or bonds is not None:
                    if dev.str_is(bonds, 'recompute'): bonds = None
                    group_mol = molecule_type(atoms, new_coords, bonds=bonds)
                    if bond_order is None:
                        bond_order = [
                            next(
                                (
                                    b[2] for b in group_mol.bonds
                                    if b[1] in (i, gs) and b[0] in (i, gs)
                                ),
                                None
                            )
                            for i,gs in group_site
                        ]
                    if len(group_site) == 1:
                        embedding = cls.fragment_embedding(
                            group_mol,
                            group_roots,
                            order=len(group_site),
                            ref=group_inds
                        )
                        origin = np.average(new_coords[group_inds,], axis=0)
                        embedding = origin, nput.view_matrix(
                            view_vector=embedding[1],
                            up_vector=embedding[2]
                        )
                    else:
                        embedding = cls.fragment_embedding(
                            group_mol,
                            group_roots,
                            order=len(group_site),
                            # ref=group_inds,
                            excluded=group_inds
                        )
                        og_1 = np.average(new_coords[group_inds,], axis=0)
                        og_2 = np.average(new_coords[group_roots,], axis=0)
                        v0 = og_2 - og_1
                        if np.dot(embedding[2], v0) < 0:
                            embedding = embedding[:2] + (embedding[2] * -1,)

                        origin = embedding[1] + og_2
                        embedding = origin, nput.view_matrix(
                            view_vector=embedding[1],
                            up_vector=embedding[2]
                        )
                elif bonds == None:
                    _, embedding = nput.moments_of_inertia(new_coords, masses=masses)
                    embedding = np.average(new_coords[group_inds,], axis=0), embedding
            rem = np.setdiff1d(np.arange(len(atoms)), group_inds)
            if bonds is not None and not dev.str_is(bonds, 'recompute'):
                remapping = {i: n for n, i in enumerate(rem)}
                bonds = [
                    [remapping[b[0]], remapping[b[1]], 1 if len(b) == 2 else b[2]]
                    for b in bonds
                    if b[0] in remapping and b[1] in remapping
                ]
            if masses is not None:
                masses = np.asanyarray(masses)[rem,]
            if bond_sites is None:
                bond_sites = [i for i,g in group_site]
            return cls.attach_functional_group(
                mol,
                target_fragment,
                [atoms[i] for i in rem],
                new_coords[rem,],
                bonds=bonds,
                ref=ref,
                masses=masses,
                distance=distance,
                angle=angle,
                dihedral=dihedral,
                dihedral_search_steps=dihedral_search_steps,
                dihedral_distance_metric=dihedral_distance_metric,
                embedding=embedding,
                bond_order=bond_order,
                use_absolue_posititions=use_absolue_posititions,
                group_site=None,
                bond_sites=bond_sites,
                molecule_type=molecule_type
            )

        if bond_order is None:
            bond_order = 1
        if bond_sites is None:
            if nput.is_int(bond_order):
                bond_sites = [0]
            else:
                bond_sites = list(range(bond_order))
        if dev.str_is(distance, 'auto'):
            if ref is None:
                ref = mol.neighborhood(target_fragment[0], size=2)
            if nput.is_int(bond_order):
                if ref[0] > -1:
                    ref_type = mol.atoms[ref[0]]
                    distance = BondData[(ref_type, atoms[bond_sites[0]], bond_order)] * UnitsData.convert("Angstroms", "BohrRadius")
                else:
                    distance = None
            else:
                distance = [
                    BondData[(mol.atoms[r], atoms[bs], bo)] * UnitsData.convert("Angstroms", "BohrRadius")
                        if r > -1 else
                    None
                    for r,bo,bs in zip(ref, bond_order, bond_sites)
                ]
        if masses is None:
            masses = np.array([AtomData[a, "Mass"] for a in atoms])

        # needed both for bond remapping below and, now, as the comparison
        # set when `dihedral='auto'` scores candidate orientations
        rem = np.setdiff1d(np.arange(len(mol.atoms)), target_fragment)

        if not use_absolue_posititions:
            if nput.is_int(bond_order):
                attachement_order = 1
            else:
                attachement_order = len(bond_order)
            origin, offset, up = cls.fragment_embedding(mol, target_fragment, ref=ref, order=attachement_order)
            # single bond attachement
            if distance is not None:
                if not nput.is_numeric(distance):
                    distance = np.average(distance)
                offset = nput.vec_normalize(offset) * distance
            if len(new_coords) == 1:
                new_coords = (origin + offset)[np.newaxis]
            else:
                if dev.str_is(embedding, 'auto'):
                    _, embedding = nput.moments_of_inertia(new_coords, masses=masses)
                centroid = np.average(new_coords[bond_sites,], axis=0)
                shift = centroid
                new_coords = new_coords - centroid[np.newaxis]
                if embedding is not None:
                    if len(embedding) == 2:
                        cent, embedding = embedding
                        cent = cent - shift
                    else:
                        cent = nput.center_of_mass(new_coords, masses=masses)
                    u = centroid - cent
                    inv = nput.view_matrix(
                        up_vector=embedding[:, 2],
                        view_vector=u
                    )
                    rot = nput.view_matrix(
                        up_vector=up,
                        view_vector=offset
                    )
                    new_coords = new_coords @ (inv @ rot.T)
                    if angle != 0:
                        new_coords = new_coords @ nput.rotation_matrix(up, angle)

                    if dev.str_is(dihedral, 'auto'):
                        dihedral = 0
                        if attachement_order == 1:
                            # numerically solve for good placement
                            metric = dihedral_distance_metric
                            if metric is None:
                                metric = functools.partial(
                                    cls._vdw_clash_metric,
                                    frag_atoms=atoms,
                                    other_atoms=[mol.atoms[i] for i in rem]
                                )
                            other_coords = mol.coords[rem]
                            base_coords = new_coords
                            best_score = -np.inf
                            best_coords = base_coords
                            for cand in np.linspace(0, np.pi/2, dihedral_search_steps, endpoint=False):
                                cand_coords = base_coords @ nput.rotation_matrix(offset, cand)
                                placed = (origin + offset)[np.newaxis] + cand_coords
                                score = metric(placed, other_coords)
                                if score > best_score:
                                    best_score = score
                                    best_coords = cand_coords
                            new_coords = best_coords
                    if dihedral != 0:
                        # rotate about offset axis, assumed to be perp to up
                        new_coords = new_coords @ nput.rotation_matrix(offset, dihedral)
                new_coords = (origin + offset)[np.newaxis] + new_coords

        remapping = {i:n for n,i in enumerate(rem)}
        if dev.str_is(bonds, 'recompute'):
            total_bonds = None
        elif bonds is None:
            total_bonds = [
                [remapping[b[0]], remapping[b[1]], 1 if len(b) == 2 else b[2]]
                for b in mol.bonds
                if b[0] in remapping and b[1] in remapping
            ]
        else:
            n = len(rem)
            total_bonds = [
                [remapping[b[0]], remapping[b[1]], 1 if len(b) == 2 else b[2]]
                for b in mol.bonds
                if b[0] in remapping and b[1] in remapping
            ] + [
                [remapping[ref[i]], n+bs, bo]
                for i,(bo,bs) in enumerate(zip(bond_order, bond_sites))
            ] + [
                [b[0]+n, b[1]+n, 1 if len(b) == 2 else b[2]]
                for b in bonds
            ]

        return mol.modify(
            atoms=tuple(mol.atoms[i] for i in rem) + tuple(atoms),
            coords=np.concatenate([
                mol.coords[rem],
                new_coords
            ]),
            bonds=total_bonds,
            masses=np.concatenate([[mol.masses[i] for i in rem], masses])
        )
    #endregion

    #region Fragment assembly
    @classmethod
    def from_fragments(
            cls,
            scaffold,
            *replacements,
            active_sites=None,
            chiralities=None,
            stereos=None,
            bond_orders=None,
            atom_replacements=None,
            cache=None,
            add_implicit_hydrogens='full',
            remove_sites=False,
            recompute_properties=True,
            reorder_from_atom_map=False,
            molecule_type=None,
            **opts
    ):
        """
        Construct a molecule from a `scaffold` plus a set of replacement fragments, driving the
        connectivity via a templated-SMILES join (`build_templated_smiles`) and the 3D geometry
        via repeated `attach_functional_group` calls, then reordering the atoms to match the
        final SMILES numbering.

        :param scaffold: the base fragment: a SMILES string, a `Molecule`, or a dict with any of
            `molecule`/`smiles`/`atoms`/`coords`/`bonds`
        :type scaffold: str | Molecule | dict
        :param replacements: the replacement fragments, each in the same accepted forms as
            `scaffold`
        :type replacements: str | Molecule | dict
        :param active_sites: forwarded to `build_templated_smiles`
        :param chiralities: forwarded to `build_templated_smiles`
        :param stereos: `{(i, j): 'cis'/'trans'}` stereochemistry constraints
        :param bond_orders: forwarded to `build_templated_smiles`
        :param atom_replacements: forwarded to `build_templated_smiles`
        :param cache: SMILES-parsing cache forwarded through
        :param add_implicit_hydrogens: implicit-hydrogen handling forwarded through
        :param remove_sites: forwarded to `build_templated_smiles`
        :param recompute_properties: whether to rebuild the final molecule cleanly from the joined
            SMILES (carrying over the assembled coordinates), or just reorder the assembled
            molecule in place
        :type recompute_properties: bool
        :param molecule_type: the molecule class to build with; defaults to the standard
            `Molecule`
        :type molecule_type: type | None
        :param opts: extra options forwarded to the molecule constructor / `from_string`
        :type opts: dict
        :return: the assembled molecule
        :rtype: Molecule
        :raises ValueError: if not every atom could be placed into the final SMILES ordering
        """
        if molecule_type is None:
            molecule_type = cls._default_molecule_type()
        Molecule = molecule_type

        def canonicalize_fragment(scaffold):
            if isinstance(scaffold, str):
                scaffold = {'smiles': scaffold}
            elif hasattr(scaffold, 'to_string'):
                scaffold = {'molecule': scaffold}
            scaffold = scaffold.copy()
            mol = scaffold.pop('molecule', None)
            if mol is None and 'smiles' in scaffold and 'coords' not in scaffold:
                mol = Molecule.from_string(scaffold.pop('smiles'), add_implicit_hydrogens='full')
            if mol is not None:
                ordering = None
                if 'smiles' not in scaffold:
                    scaffold['smiles'], ordering = mol.to_string('smi', remove_hydrogens=True,
                                                                 return_reordering=True)
                    rem = np.setdiff1d(np.arange(len(mol.atoms)), ordering)
                    ordering = np.concatenate([ordering, rem], axis=0)
                if 'atoms' not in scaffold:
                    scaffold['atoms'] = mol.atoms
                    if ordering is not None:
                        a = scaffold['atoms']
                        scaffold['atoms'] = [a[i] for i in ordering]
                if 'coords' not in scaffold:
                    scaffold['coords'] = mol.coords
                    if ordering is not None:
                        scaffold['coords'] =scaffold['coords'][ordering,]
                if 'bonds' not in scaffold:
                    scaffold['bonds'] = mol.bonds
                    if ordering is not None:
                        map = np.argsort(ordering)
                        scaffold['bonds'] = [
                            (map[b[0]], map[b[1]]) + tuple(b[2:])
                            for b in scaffold['bonds']
                        ]
            return scaffold

        frags = [canonicalize_fragment(scaffold)] + [canonicalize_fragment(f) for f in replacements]

        # `build_templated_smiles`/`join_smiles_fragments` only care about the
        # SMILES-level bookkeeping, so strip the 3D-only keys before forwarding
        geom_keys = ('coords', 'bonds', 'atoms')
        smiles_frags = [
            {k: v for k, v in f.items() if k not in geom_keys}
            for f in frags
        ]
        for s in smiles_frags:
            s['functional_group'] = s.pop('smiles')


        smiles_options = dict(
            active_sites=active_sites,
            chiralities=chiralities,
            stereos=stereos,
            bond_orders=bond_orders,
            atom_replacements=atom_replacements,
            cache=cache,
            add_implicit_hydrogens=add_implicit_hydrogens,
            remove_sites=remove_sites,
            reorder_from_atom_map=reorder_from_atom_map,
            return_fragment_indices=True,
            return_new_bonds=True
        )

        smiles, atom_maps, bond_maps = build_templated_smiles(
            smiles_frags[0]['functional_group'],
            *smiles_frags[1:],
            **smiles_options
        )

        def heavy_atom_positions(atoms):
            # SMILES fragments are built with `remove_hydrogens=True`, so
            # `atom_maps`/`bond_maps` only ever reference heavy atoms; this
            # recovers the correspondence with the (H-inclusive) 3D atom list
            return [i for i, a in enumerate(atoms) if a not in ('H', 'D', 'T')]

        def find_bonded_hydrogen(atoms, coords, bonds, heavy_local, exclude=(), return_all=False):
            if bonds is None:
                bonds = Molecule(atoms, coords).bonds
            exclude = set(exclude)
            candidates = []
            for b in bonds:
                a, c = b[0], b[1]
                if a == heavy_local:
                    if atoms[c] in ('H', 'D', 'T') and c not in exclude:
                        if return_all:
                            candidates.append(c)
                        else:
                            return c
                elif c == heavy_local:
                    if atoms[a] in ('H', 'D', 'T') and a not in exclude:
                        if return_all:
                            candidates.append(a)
                        else:
                            return a
            if return_all:
                return candidates
            else:
                return None

        base = frags[0]
        mol = Molecule(base['atoms'], base['coords'], bonds=base.get('bonds'), **opts)

        heavy0 = heavy_atom_positions(base['atoms'])
        # map: final (fully-joined) SMILES heavy-atom index -> current position in `mol`
        current_index = {
            atom_maps[0][i]: heavy0[i]
            for i in range(len(atom_maps[0]))
        }

        for f in frags:
            f['sites'] = parse_smiles_and_atom_map(f['smiles'],
                                                   cache=smiles_options['cache'],
                                                   add_implicit_hydrogens=smiles_options['add_implicit_hydrogens'],
                                                   reorder_from_atom_map=smiles_options['reorder_from_atom_map'])['map']

        cur_sites = frags[0]['sites']
        if stereos is None:
            stereos = {}

        for frag, am, bonds_i in zip(frags[1:], atom_maps[1:], bond_maps):
            if not bonds_i:
                continue

            heavy_frag = heavy_atom_positions(frag['atoms'])
            heavy_to_smiles_local = {loc: k for k, loc in enumerate(heavy_frag)}
            sub_index = {final: i for i, final in enumerate(am)}

            frag_sites = frag['sites']
            frag_bonds_raw = frag.get('bonds') or Molecule(frag['atoms'], frag['coords']).bonds

            # `bonds_i` may contain more than one bond (e.g. a bidentate ligand that
            # closes a ring / forms a linker). Discover a scaffold target atom and a
            # fragment placeholder hydrogen for *each* bond, removing an explicit H at
            # each site where one exists. Hydrogens already claimed by an earlier bond
            # in this same join are not reused.
            base_final_of = {v: k for k, v in current_index.items()}  # local pos -> final idx, base side
            base_substereos = {
                (cur_sites[i + 1], cur_sites[j + 1]): v
                for (i, j), v in stereos.items()
                if (i + 1 in cur_sites) and (j + 1 in cur_sites)
            }
            frag_substereos = {
                (frag_sites[i + 1], frag_sites[j + 1]): v
                for (i, j), v in stereos.items()
                if (i + 1 in frag_sites) and (j + 1 in frag_sites)
            }
            final_of_frag = {loc: am[k] for k, loc in enumerate(heavy_frag)}  # local pos -> final idx, frag side

            used_base_h = set()
            used_frag_h = set()
            sites = []  # one entry per bond: (target, base_heavy, frag_h, frag_local, bond_order)
            if len(bonds_i) == 1:
                bond = bonds_i[0]
                base_final, frag_final, *rest = bond
                bo = rest[0] if rest else 1

                base_heavy = current_index[base_final]
                frag_local = heavy_to_smiles_local[heavy_frag[sub_index[frag_final]]]

                # the atom actually replaced is the explicit hydrogen at this scaffold
                # site (if any) -- NOT the heavy atom -- since the target is deleted
                stereo_h = cls.resolve_stereo_hydrogen(
                    mol.atoms, mol.coords, mol.bonds,
                    base_heavy, base_substereos, base_final_of
                )
                if stereo_h is not None and stereo_h in used_base_h:
                    stereo_h = None  # already claimed by an earlier bond in this join
                base_h = stereo_h if stereo_h is not None else find_bonded_hydrogen(
                    mol.atoms, mol.coords, mol.bonds, base_heavy, exclude=used_base_h
                )
                if base_h is not None:
                    used_base_h.add(base_h)
                target = base_h if base_h is not None else base_heavy

                stereo_h = cls.resolve_stereo_hydrogen(
                    frag['atoms'], frag['coords'], frag_bonds_raw,
                    frag_local, frag_substereos, final_of_frag
                )
                if stereo_h is not None and stereo_h in used_frag_h:
                    stereo_h = None
                frag_h = stereo_h if stereo_h is not None else find_bonded_hydrogen(
                    frag['atoms'], frag['coords'], frag.get('bonds'), frag_local,
                    exclude=used_frag_h
                )
                if frag_h is not None:
                    used_frag_h.add(frag_h)

                sites.append((target, base_heavy, frag_h, frag_local, bo))
            else:
                # find all hydrogen candidates from `bonds_i`
                targets, frag_targs = [], []
                base_candidates, frag_candidates = [], []
                bond_orders = []
                for bond in bonds_i:
                    base_final, frag_final, *rest = bond
                    bo = rest[0] if rest else 1
                    bond_orders.append(bo)

                    base_heavy = current_index[base_final]
                    frag_local = heavy_to_smiles_local[heavy_frag[sub_index[frag_final]]]

                    targets.append(base_heavy)
                    bh = find_bonded_hydrogen(
                            mol.atoms, mol.coords, mol.bonds, base_heavy, exclude=used_base_h,
                            return_all=True
                        )
                    if len(bh) == 0: bh = [None]
                    base_candidates.append(bh)

                    frag_targs.append(frag_local)
                    fh = find_bonded_hydrogen(
                        frag['atoms'], frag['coords'], frag.get('bonds'), frag_local,
                        exclude=used_frag_h,
                        return_all=True
                    )
                    if len(fh) == 0: fh = [None]
                    frag_candidates.append(fh)

                # choose candidates such that their sum of pairwise distances is minimized

                best_base = None
                min_dist = np.inf
                for cands in itut.unique_product(*base_candidates):
                    cands = [c for c in cands if c is not None]
                    if len(cands) == 0:
                        best_base = cands
                        break
                    else:
                        pos = mol.coords[cands,]
                        dists = np.sum(nput.distance_matrix(pos, return_triu=True))
                        if dists < min_dist:
                            min_dist = dists
                            best_base = cands
                if best_base is None:
                    raise ValueError("could not find candidate base hydrogens")

                best_frag = None
                min_dist = np.inf
                for cands in itut.unique_product(*frag_candidates):
                    cands = [c for c in cands if c is not None]
                    if len(cands) == 0:
                        best_frag = cands
                        break
                    else:
                        pos = frag['coords'][cands,]
                        dists = np.sum(nput.distance_matrix(pos, return_triu=True))
                        if dists < min_dist:
                            min_dist = dists
                            best_frag = cands
                if best_frag is None:
                    raise ValueError("could not find candidate fragment hydrogens")

                print(best_base, best_frag)
                print(targets, frag_targs)
                # best_base = [0, 3]
                # best_frag = [5, 8]

                for target, base_heavy, frag_h, frag_local, bo in zip(best_base, targets, best_frag, frag_targs, bond_orders):
                    sites.append((target, base_heavy, frag_h, frag_local, bo))

            n_before = len(mol.atoms)

            if len(sites) == 1:
                # ---- single-bond join: unchanged single-site path ----
                target, base_heavy, frag_h, frag_local, bond_order = sites[0]

                # `group_site` mode assumes the true attachment atom sits at local index 0
                order = [frag_local] + [j for j in range(len(frag['atoms'])) if j != frag_local]
                remap = {old: new for new, old in enumerate(order)}
                frag_atoms = [frag['atoms'][j] for j in order]
                frag_coords = frag['coords'][order,]
                frag_bonds = frag.get('bonds')
                if frag_bonds is not None:
                    frag_bonds = [
                        [remap[b[0]], remap[b[1]]] + list(b[2:])
                        for b in frag_bonds
                    ]
                group_site = remap[frag_h] if frag_h is not None else None

                mol = cls.attach_functional_group(
                    mol,
                    [target],
                    frag_atoms,
                    frag_coords,
                    bonds=frag_bonds if frag_bonds is not None else 'recompute',
                    bond_order=bond_order,
                    group_site=group_site,
                    molecule_type=molecule_type
                )

                removed_targets = [target]
                # the fragment's atoms were appended in `order`, minus whichever
                # position `group_site` removed
                appended = [j for j in order if frag_h is None or j != frag_h]
            else:
                # ---- multi-bond join: multi-site placeholder attachment ----
                # every fragment site must expose an explicit H to remove so it can be
                # cast as a `group_site` placeholder; likewise every scaffold site needs
                # an H to remove (otherwise there's no unambiguous atom to delete)
                if any(frag_h is None for (_, _, frag_h, _, _) in sites):
                    raise ValueError(
                        "multi-bond fragment join needs an explicit hydrogen at every "
                        "fragment bond site to remove; none found for some site"
                    )
                if any(t == bh for (t, bh, _, _, _) in sites):
                    raise ValueError(
                        "multi-bond fragment join needs an explicit hydrogen at every "
                        "scaffold bond site to remove; none found for some site"
                    )

                frag_atoms = list(frag['atoms'])
                frag_coords = np.asanyarray(frag['coords'], dtype=float)
                frag_bonds = frag.get('bonds')

                targets = [t for (t, _, _, _, _) in sites]
                # refs = [bh for (_, bh, _, _, _) in sites]
                group_sites = [(fl, fh) for (_, _, fh, fl, _) in sites]
                bond_orders = [bo for (_, _, _, _, bo) in sites]

                print(group_sites)
                print(frag_atoms)
                mol = cls.attach_functional_group(
                    mol,
                    targets,
                    frag_atoms,
                    frag_coords,
                    bonds=frag_bonds if frag_bonds is not None else 'recompute',
                    # ref=refs,
                    bond_order=bond_orders,
                    group_site=group_sites,
                    molecule_type=molecule_type
                )

                removed_targets = [t for (t, _, _, _, _) in sites]
                # in the multi-site path the surviving fragment atoms are appended in
                # ascending local order, minus the removed placeholders
                order = list(range(len(frag_atoms)))
                appended = [j for j in order if j not in set(group_sites)]

            # --- update current_index after removing the scaffold target(s) ---
            removed_sorted = sorted(removed_targets)
            def _shift(idx):
                return idx - sum(1 for t in removed_sorted if t < idx)
            current_index = {
                final: _shift(idx)
                for final, idx in current_index.items()
                if idx not in removed_targets
            }
            base_offset = n_before - len(removed_targets)  # atom count after removals

            # the fragment heavy atoms were appended after the surviving scaffold atoms,
            # in the order given by `appended`
            for pos, local_idx in enumerate(appended):
                if local_idx in heavy_to_smiles_local:
                    final_idx = am[heavy_to_smiles_local[local_idx]]
                    current_index[final_idx] = base_offset + pos
                # explicit hydrogens carried on the fragment have no counterpart in the
                # (hydrogen-free) SMILES numbering, so they simply aren't tracked here

            # join cur_sites and frag_sites
            site_offset = len(cur_sites)
            cur_sites = cur_sites | {i + site_offset: f + base_offset for i, f in frag_sites.items()}

        def _extend_with_hydrogen_indices(mol, current_index, n_heavy_total):
            """
            current_index: {final heavy-atom index -> current position in `mol`}
            Extends it with hydrogens: {n_heavy_total + rank -> current position},
            where `rank` orders H's by the final index of the heavy atom each is
            bonded to (ties broken by their existing relative order in `mol`).
            """
            pos_to_final_heavy = {pos: final for final, pos in current_index.items()}

            def heavy_neighbor(i):
                return next(
                    (b[1] if b[0] == i else b[0])
                    for b in mol.bonds
                    if i in b[:2]
                    and mol.atoms[b[1] if b[0] == i else b[0]] not in ('H', 'D', 'T')
                )

            h_entries = [
                (pos_to_final_heavy[heavy_neighbor(i)], i)
                for i, a in enumerate(mol.atoms)
                if a in ('H', 'D', 'T')
            ]
            # stable sort: ties (same attachment atom) keep current relative order
            h_entries.sort(key=lambda e: e[0])

            full_index = dict(current_index)
            for rank, (_, pos) in enumerate(h_entries):
                full_index[n_heavy_total + rank] = pos
            return full_index

        n_heavy_total = max(max(am) for am in atom_maps) + 1
        full_index = _extend_with_hydrogen_indices(mol, current_index, n_heavy_total)

        n_total = len(mol.atoms)
        if len(full_index) != n_total:
            missing = set(range(n_total)) - set(full_index.values())
            raise ValueError(
                f"Could not place every atom into the final SMILES ordering; "
                f"unaccounted current positions: {missing}"
            )

        new_ord = [full_index[i] for i in range(n_total)]

        if recompute_properties:
            new_mol = Molecule.from_string(smiles, 'smi', reorder_from_atom_map=False, **opts)
            new_mol.coords = mol.coords[new_ord,]
        else:
            new_mol = mol.take_submolecule(new_ord)

        return new_mol
    #endregion
